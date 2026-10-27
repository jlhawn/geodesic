#!/usr/bin/env python3
"""Writes data/aerosol_mac_2deg.bin, the aerosol climatology that radiation
option aerosol 'mac' reads (decodeAerosol in js/aerosol.module.js).

Source: the Max Planck Institute for Meteorology's MACv2 aerosol climatology,
reference year 2005, anonymous FTP, ftp://ftp-projects.mpimet.mpg.de/aerocom/climatology/
under MACv2_2018/ (the server states no license; the attribution the file is
distributed with is in docs/aerosol-data-notice.md):

  for_ICON/spec_2005/by_wave/gt_f_00550nm.nc        fine mode, 550 nm   (10.1 MB)
  for_ICON/spec_2005/by_wave/gt_c_00550nm.nc        coarse mode, 550 nm (10.1 MB)
  for_ICON/altitude/aeropt_fin_alt_km20.nc          fine mode's AOD fraction in 20 layers (62.5 MB)
  for_ICON/altitude/aeropt_coa_alt_km20.nc          the same for the coarse mode (62.5 MB)

Kinne (2019), The MACv2 aerosol climatology, Tellus B 71, 1623639,
https://doi.org/10.1080/16000889.2019.1623639.

Usage (standard library only):

  python3 -I scripts/packAerosol.py RAW_DIR [out.bin] [FACTOR]

RAW_DIR holds the four files under their names. FACTOR (2 by default) is the
number of source cells per side that make one cell of the output; the source
grid is 1 degree, 360 x 180, and its size is read from the files.

Each altitude file may instead be a record file, aeropt_fin_alt_km20.record or
aeropt_coa_alt_km20.record: one month's record cut from the file by byte range,
a big-endian float32 `time` and then the 20 layers' fractions, nLat x nLon each.
It gives every month the same heights, and, there being no surface altitude in
it, the ground of a cell is taken as the lower edge of its lowest layer that
holds aerosol. A complete file gives each month its own heights over its `asl`.

Per cell and month the file holds eight fields of one byte, each as
low + (high - low) * s, s = q / 255 (linear) or (q / 255)^2 (square):
  0, 1  AOD at 550 nm of the fine and the coarse mode, 0 to 1.2, square scale
  2, 3  1 - single scattering albedo at 550 nm, 0 to 0.3
  4, 5  asymmetry parameter at 550 nm, 0.45 to 0.85
  6, 7  mean height of the mode's aerosol above the ground, km, 0 to 8
Coarsening takes the AOD by area, 1 - SSA weighted by the AOD, the asymmetry
by AOD times SSA and the height by AOD.

Layout, little-endian, after the header of packWoaMonthlySst.py's 'WOAM':
  0   char[4]   'AERM'
  4   uint32    byte offset of the data
  8   uint16    nLon     10  uint16  nLat     12  uint16  nMonth     14  uint16  nField
  16  float32   lon0, dLon, lat0, dLat   (degrees; center of column 0 / row 0)
  32  float32   center[nMonth]  (day of a 365-day year, 0 = January 1 00:00)
      per field: float32 low, float32 high, uint32 scale (0 linear, 1 square)
      uint16    n, then n bytes of ASCII provenance
      zero padding to a multiple of 8
  data: uint8 [field][month][lat][lon], rows south to north, columns west to
        east from lon0.
"""
import array
import math
import struct
import sys
from pathlib import Path

N_FIELD = 8
RANGES = [(0.0, 1.2, 1), (0.0, 1.2, 1), (0.0, 0.3, 0), (0.0, 0.3, 0), (0.45, 0.85, 0), (0.45, 0.85, 0), (0.0, 8.0, 0), (0.0, 8.0, 0)]
LAYER_BOTTOM = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 15000]
LAYER_TOP = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 15000, 20000]
TYPES = {1: ('b', 1), 2: ('c', 1), 3: ('h', 2), 4: ('i', 4), 5: ('f', 4), 6: ('d', 8)}


class Classic:
    """The header and variables of a netCDF-3 classic file (CDF-1 or CDF-2)."""

    def __init__(self, path):
        self.path = path
        with open(path, 'rb') as f:
            self.head = f.read(1 << 20)
        if self.head[:3] != b'CDF' or self.head[3] not in (1, 2):
            raise ValueError(f'{path}: not a netCDF-3 classic file')
        self.version, self.at = self.head[3], 4
        self.numrecs = self.u32()
        self.dims = self.read_dims()
        self.read_atts()
        self.vars = self.read_vars()
        records = [v for v in self.vars.values() if v['record']]
        self.recsize = sum(v['vsize'] for v in records)
        if len(records) == 1:
            self.recsize = records[0]['width']

    def u32(self):
        value = struct.unpack_from('>I', self.head, self.at)[0]
        self.at += 4
        return value

    def u64(self):
        value = struct.unpack_from('>Q', self.head, self.at)[0]
        self.at += 8
        return value

    def name(self):
        n = self.u32()
        text = self.head[self.at:self.at + n].decode('latin1')
        self.at += (n + 3) & ~3
        return text

    def read_dims(self):
        self.u32()
        return [(self.name(), self.u32()) for _ in range(self.u32())]

    def read_atts(self):
        self.u32()
        for _ in range(self.u32()):
            self.name()
            size = TYPES[self.u32()][1]
            count = self.u32()
            self.at += (count * size + 3) & ~3

    def read_vars(self):
        self.u32()
        out = {}
        for _ in range(self.u32()):
            name = self.name()
            ids = [self.u32() for _ in range(self.u32())]
            self.read_atts()
            kind, vsize = self.u32(), self.u32()
            begin = self.u32() if self.version == 1 else self.u64()
            record = bool(ids) and self.dims[ids[0]][1] == 0
            shape = [self.dims[i][1] or self.numrecs for i in ids]
            width = TYPES[kind][1]
            for n in shape[1:] if record else shape:
                width *= n
            out[name] = {'kind': kind, 'shape': shape, 'begin': begin, 'record': record, 'width': width, 'vsize': vsize}
        return out

    def read(self, name, record=0):
        """A variable (one record of a record variable) as an array of its values."""
        v = self.vars[name]
        code, size = TYPES[v['kind']]
        count = v['width'] // size
        start = v['begin'] + (record * self.recsize if v['record'] else 0)
        with open(self.path, 'rb') as f:
            f.seek(start)
            raw = f.read(count * size)
        if len(raw) != count * size:
            raise ValueError(f'{self.path}: {name} is cut short')
        values = array.array(code)
        values.frombytes(raw)
        if sys.byteorder == 'little':
            values.byteswap()
        return values


def read_record(path, layers, cells):
    """A one-month altitude record: float32 time, then layers x cells fractions."""
    with open(path, 'rb') as f:
        raw = f.read()
    if len(raw) != 4 + 4 * layers * cells:
        raise ValueError(f'{path}: {len(raw)} bytes, not {4 + 4 * layers * cells}')
    values = array.array('f')
    values.frombytes(raw[4:])
    if sys.byteorder == 'little':
        values.byteswap()
    return [values[l * cells:(l + 1) * cells] for l in range(layers)]


def finite(x, ceiling=1e20):
    return x == x and abs(x) < ceiling


class Source:
    """The 550 nm optics and the layers of one mode on the source grid."""

    def __init__(self, raw, mode):
        optics = Classic(raw / f'gt_{mode[0]}_00550nm.nc')
        self.optics = optics
        dims = dict(optics.dims)
        self.nLon, self.nLat = dims['lon'], dims['lat']
        self.lat = list(optics.read('lat'))
        self.lon = list(optics.read('lon'))
        self.months = optics.numrecs
        self.centers = [float(optics.read('time', m)[0]) for m in range(self.months)]
        name = 'aeropt_fin_alt_km20' if mode == 'fine' else 'aeropt_coa_alt_km20'
        self.layers_name = 'zkm_aer_fine_mo' if mode == 'fine' else 'zkm_aer_coarse_mo'
        full, record = raw / f'{name}.nc', raw / f'{name}.record'
        cells = self.nLat * self.nLon
        if full.exists():
            self.altitude, self.record = Classic(full), None
            self.top = [float(z) for z in self.altitude.read('ztop')]
            self.bottom = [float(z) for z in self.altitude.read('zbot')]
            self.asl = self.altitude.read('asl')
        elif record.exists():
            self.altitude, self.record = None, read_record(record, len(LAYER_TOP), cells)
            self.top, self.bottom, self.asl = [float(z) for z in LAYER_TOP], [float(z) for z in LAYER_BOTTOM], None
        else:
            raise FileNotFoundError(f'{raw}: neither {name}.nc nor {name}.record')

    def month(self, m):
        """tau, ssa, asy and the mean height above ground (km) of month m, flat over the source grid."""
        cells = self.nLat * self.nLon
        tau, ssa, asy = (self.optics.read(v, m) for v in ('aod', 'ssa', 'asy'))
        if self.record:
            fractions = self.record
        else:
            layers = self.altitude.read(self.layers_name, m)
            fractions = [layers[l * cells:(l + 1) * cells] for l in range(len(self.top))]
        middle = [0.5 * (a + b) / 1000 for a, b in zip(self.top, self.bottom)]
        edges = [b / 1000 for b in self.bottom]
        height = array.array('d', [0.0]) * cells
        for c in range(cells):
            total = weighted = 0.0
            lowest = -1
            for l, layer in enumerate(fractions):
                f = layer[c]
                if f > 0:
                    if lowest < 0:
                        lowest = l
                    total += f
                    weighted += f * middle[l]
            if total > 0:
                ground = self.asl[c] / 1000 if self.asl is not None else edges[lowest]
                height[c] = max(0.0, weighted / total - ground)
        return tau, ssa, asy, height


def coarsen(source, month, factor):
    """One month on the output grid, rows south to north: AOD, 1 - SSA, ASY and height."""
    tau, ssa, asy, height = month
    nLon, nLat = source.nLon, source.nLat
    rows, cols = nLat // factor, nLon // factor
    flip = source.lat[0] > source.lat[-1]
    out = [array.array('d', [0.0]) * (rows * cols) for _ in range(4)]
    for j in range(rows):
        for i in range(cols):
            area = depth = scatter = tilt = lift = 0.0
            plain = [0.0, 0.0, 0.0]
            for dj in range(factor):
                r = (nLat - 1 - (j * factor + dj)) if flip else (j * factor + dj)
                w = max(0.0, math.cos(math.radians(source.lat[r])))
                for di in range(factor):
                    c = r * nLon + i * factor + di
                    a, s, g, h = tau[c], ssa[c], asy[c], height[c]
                    if not (finite(a) and finite(s) and finite(g)):
                        continue
                    a = max(0.0, a)
                    area += w
                    depth += w * a
                    scatter += w * a * s
                    tilt += w * a * s * g
                    lift += w * a * h
                    plain[0] += w * s
                    plain[1] += w * g
                    plain[2] += w * h
            n = j * cols + i
            out[0][n] = depth / area if area > 0 else 0.0
            ssa_mean = scatter / depth if depth > 0 else (plain[0] / area if area > 0 else 1.0)
            out[1][n] = 1.0 - ssa_mean
            out[2][n] = tilt / scatter if scatter > 0 else (plain[1] / area if area > 0 else 0.7)
            out[3][n] = lift / depth if depth > 0 else (plain[2] / area if area > 0 else 1.0)
    return out


def quantize(values, low, high, square, clipped):
    out = bytearray(len(values))
    for n, v in enumerate(values):
        s = (v - low) / (high - low)
        if s < 0 or s > 1:
            clipped[0] += 1
            s = min(1.0, max(0.0, s))
        out[n] = int(round(255 * (math.sqrt(s) if square else s)))
    return out


def main(raw_dir, out_path, factor):
    raw = Path(raw_dir)
    fine, coarse = modes = [Source(raw, 'fine'), Source(raw, 'coarse')]
    if (fine.nLon, fine.nLat, fine.centers) != (coarse.nLon, coarse.nLat, coarse.centers):
        raise ValueError('the fine and coarse files differ in grid or months')
    nLon, nLat, months = fine.nLon, fine.nLat, fine.months
    if nLon % factor or nLat % factor:
        raise ValueError(f'a factor of {factor} does not divide the {nLon} x {nLat} source grid')
    rows, cols = nLat // factor, nLon // factor
    step_lon, step_lat = abs(fine.lon[1] - fine.lon[0]), abs(fine.lat[1] - fine.lat[0])
    lon0, lat0 = min(fine.lon) + 0.5 * (factor - 1) * step_lon, min(fine.lat) + 0.5 * (factor - 1) * step_lat
    planes = [[None] * months for _ in range(N_FIELD)]
    for m in range(months):
        for k, source in enumerate(modes):
            tau, absorption, asymmetry, height = coarsen(source, source.month(m), factor)
            for f, values in ((k, tau), (2 + k, absorption), (4 + k, asymmetry), (6 + k, height)):
                planes[f][m] = values
        print(f'month {m + 1:2d} of {months}: mean fine AOD {sum(planes[0][m]) / len(planes[0][m]):.4f}, coarse {sum(planes[1][m]) / len(planes[1][m]):.4f}')
    clipped = [0] * N_FIELD
    data = bytearray()
    for f in range(N_FIELD):
        low, high, square = RANGES[f]
        count = [0]
        for m in range(months):
            data += quantize(planes[f][m], low, high, square, count)
        clipped[f] = count[0]
    print('values outside the stored range, clipped, per field:', clipped)
    source = ('MACv2 aerosol climatology, Kinne 2019, Tellus B 71, 1623639; reference year 2005; MPI-M '
              'ftp-projects.mpimet.mpg.de/aerocom/climatology/MACv2_2018 (no license stated); '
              f'{factor}-degree cells, AOD by area, 1-SSA by AOD, ASY by AOD*SSA, height by AOD').encode('ascii')
    head = struct.pack('<4sIHHHH', b'AERM', 0, cols, rows, months, N_FIELD)
    head += struct.pack('<4f', lon0, step_lon * factor, lat0, step_lat * factor)
    head += struct.pack(f'<{months}f', *fine.centers)
    for low, high, square in RANGES:
        head += struct.pack('<ffI', low, high, square)
    head += struct.pack('<H', len(source)) + source
    head += b'\0' * (-len(head) % 8)
    head = head[:4] + struct.pack('<I', len(head)) + head[8:]
    with open(out_path, 'wb') as out:
        out.write(head)
        out.write(data)
    print(f'{out_path}: {Path(out_path).stat().st_size} bytes, {cols} x {rows} x {months} months x {N_FIELD} fields, data at byte {len(head)}')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    default = Path(__file__).resolve().parent.parent / 'data' / 'aerosol_mac_2deg.bin'
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else str(default), int(sys.argv[3]) if len(sys.argv) > 3 else 2)
