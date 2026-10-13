#!/usr/bin/env python3
"""Writes data/woa_sst_monthly_1deg.bin, the monthly sea surface temperature
climatology the SST check reads (decodeMonthlySst in
js/ocean/climatology.module.js; scripts/woaSstBias.mjs).

Source: the World Ocean Atlas 2023 (NOAA National Centers for Environmental
Information; public domain), months "01".."12", all decades ("decav"),
objectively analyzed temperature (t_an) on the 1.00 degree grid, in the
atlas's comma-separated form, which holds t_an alone, one file per month
MM = 01..12, woa23_decav_tMMan01.csv.gz (about 5.1 MB each, 61.4 MB for
the twelve):

  https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/temperature/csv/decav/1.00/woa23_decav_tMMan01.csv.gz

Locarnini et al. (2023), WOA23 Volume 1: Temperature, NOAA Atlas NESDIS 89,
https://doi.org/10.25923/54bh-1613.

Usage (standard library only):

  python3 -I scripts/packWoaMonthlySst.py RAW_DIR [out.bin]

RAW_DIR holds the twelve woa23_decav_tMMan01.csv.gz files. Only the 0 m
value of each row is kept; at the surface the in-situ temperature is the
potential temperature, so no conversion is applied.

Layout, little-endian, the header of packWoa.py's 'WOA1' with months in
place of depths and no salinity block:
  0   char[4]   'WOAM'
  4   uint32    byte offset of the data
  8   uint16    nLon     10  uint16  nLat     12  uint16  nMonth (12)
  14  int16     missing sentinel (land, ice shelf, or no analysis)
  16  float32   lon0, dLon, lat0, dLat   (degrees; center of column 0 / row 0)
  32  float32   tOffset, tScale, 0, 0
  48  float32   center[nMonth]  (day of a 365-day year, 0 = January 1 00:00)
      uint16    n, then n bytes of ASCII provenance
      zero padding to a multiple of 8
  data: int16 T[nMonth][nLat][nLon], value = offset + scale * stored (degC),
        rows south to north, columns west to east from lon0.
"""
import array
import gzip
import struct
import sys
from pathlib import Path

N_LON, N_LAT = 360, 180
LON0, LAT0, STEP = -179.5, -89.5, 1.0
T_OFFSET, T_SCALE = 15.0, 0.001
MISSING = -32768
LENGTHS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def centers():
    out, start = [], 0
    for n in LENGTHS:
        out.append(start + 0.5 * n)
        start += n
    return out


def read_month(path):
    grid = array.array('h', [MISSING]) * (N_LON * N_LAT)
    count, low, high = 0, 99.0, -99.0
    with gzip.open(path, 'rt', encoding='ascii') as rows:
        for line in rows:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split(',')
            if len(fields) < 3 or fields[2].strip() == '':
                continue
            lat, lon, t = float(fields[0]), float(fields[1]), float(fields[2])
            row, col = round((lat - LAT0) / STEP), round((lon - LON0) / STEP)
            if not (0 <= row < N_LAT and 0 <= col < N_LON) or abs(LAT0 + row * STEP - lat) > 1e-3 or abs(LON0 + col * STEP - lon) > 1e-3:
                raise ValueError(f'{path}: point {lat}, {lon} is not on the 1 degree grid')
            stored = round((t - T_OFFSET) / T_SCALE)
            if not (MISSING < stored <= 32767):
                raise ValueError(f'{path}: {t} degC outside the int16 range')
            grid[row * N_LON + col] = stored
            count += 1
            low, high = min(low, t), max(high, t)
    return grid, count, low, high


def main(raw_dir, out_path):
    raw = Path(raw_dir)
    months = []
    for m in range(1, 13):
        path = raw / f'woa23_decav_t{m:02d}an01.csv.gz'
        grid, count, low, high = read_month(path)
        print(f'{path.name}: {count} surface points, {low:.3f}..{high:.3f} degC')
        months.append(grid)
    source = ('WOA23 decav monthly (t01-t12) 1.00 degree t_an at 0 m, csv form '
              '(NOAA NCEI, public domain); degC').encode('ascii')
    head = struct.pack('<4sIHHHh', b'WOAM', 0, N_LON, N_LAT, 12, MISSING)
    head += struct.pack('<8f', LON0, STEP, LAT0, STEP, T_OFFSET, T_SCALE, 0.0, 0.0)
    head += struct.pack('<12f', *centers())
    head += struct.pack('<H', len(source)) + source
    head += b'\0' * (-len(head) % 8)
    head = head[:4] + struct.pack('<I', len(head)) + head[8:]
    with open(out_path, 'wb') as out:
        out.write(head)
        for grid in months:
            if sys.byteorder != 'little':
                grid.byteswap()
            out.write(grid.tobytes())
    print(f'{out_path}: {Path(out_path).stat().st_size} bytes, {N_LON} x {N_LAT} x 12, data at byte {len(head)}')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else str(Path(__file__).resolve().parent.parent / 'data' / 'woa_sst_monthly_1deg.bin'))
