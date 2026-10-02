#!/usr/bin/env python3
"""
The fine-scale terrain statistics the subgrid orography schemes need,
from GMTED2010's 30-arc-second mean elevation (Danielson and Gesch 2011,
USGS Open-File Report 2011-1073), filtered as the IFS documents its own
fields (Cy47r3 Part IV §11.3). numpy and Pillow only.

  assemble   the 108 tiles of the cache (20° × 30°, 2400 × 3600 int16,
             RasterPixelIsArea, nodata −32768) into one global
             21600 × 43200 int16 grid, rows from the north, columns from
             −180°, each tile's size, pixel scale, tie point and embedded
             statistics checked
  filter     the 2′30″ grids (4320 × 8640, float32, rows from the north)
             that scripts/subgridTerrain.mjs aggregates per mesh cell:
               h5     the 30″ orography (sea and nodata at 0) smoothed
                      with the IFS's operator (eq. 11.4) at Δ 5 km, δ 1 km,
                      sampled at the block centres (§11.3.2 steps 1–2)
               flt2   the block mean of (h₂ − h₂₀)², h_Δ the 30″
                      orography smoothed at Δ 2 km and Δ 20 km, δ 1 km
                      (§11.3.3 step i)
               mean   the block mean of the 30″ orography
               land   the share of the block's 30″ points above 0 m
  spectrum   the land orography's one-dimensional power spectra along
             rows and columns and their slopes

The smoothing is a convolution with the radial kernel h(r) on the sphere's
local plane, applied by FFT in bands of rows: within a band (at most 1°,
and short enough that cos φ changes by at most 2 % across it) the kernel
is laid out at the band's mean latitude (east–west spacing R cos φ Δλ) and
normalised to unit sum, and the band is periodic in longitude. Within
0.03° of a pole the kernel spans the whole row.

  python3 scripts/subgridTerrain.py assemble CACHE
  python3 scripts/subgridTerrain.py filter CACHE
  python3 scripts/subgridTerrain.py spectrum CACHE
"""
import os
import re
import sys
import time
import numpy as np
from PIL import Image

RADIUS = 6371220.0
ROWS, COLS = 21600, 43200
BLOCK = 5
TILE_ROWS, TILE_COLS = 2400, 3600
NODATA = -32768
LATS = ['90S', '70S', '50S', '30S', '10S', '10N', '30N', '50N', '70N']
LONS = ['180W', '150W', '120W', '090W', '060W', '030W', '000E', '030E', '060E', '090E', '120E', '150E']


def tile_name(lat, lon):
    return f'{lat}{lon}_20101117_gmted_mea300.tif'


def degrees(code):
    value = int(code[:-1])
    return -value if code[-1] in 'SW' else value


def assemble(cache):
    out = np.lib.format.open_memmap(os.path.join(cache, 'gmted_mea300.npy'), mode='w+', dtype=np.int16, shape=(ROWS, COLS))
    report = []
    for lat in LATS:
        for lon in LONS:
            path = os.path.join(cache, tile_name(lat, lon))
            image = Image.open(path)
            tags = image.tag_v2
            assert image.size == (TILE_COLS, TILE_ROWS), (path, image.size)
            assert tags[258] == (16,) and tags[339] == (2,) and tags[259] == 1, path
            scale = tags[33550]
            assert abs(scale[0] - 1 / 120) < 1e-12 and abs(scale[1] - 1 / 120) < 1e-12, (path, scale)
            tie = tags[33922]
            south, west = degrees(lat), degrees(lon)
            offset = 0.5 / 3600
            assert abs(tie[3] - (west - offset)) < 1e-9 and abs(tie[4] - (south + 20 - offset)) < 1e-9, (path, tie)
            keys = tags[34735]
            assert 4326 in keys, path
            data = np.array(image).astype(np.int16)
            meta = tags.get(42112, '')
            stats = {m.group(1): float(m.group(2)) for m in re.finditer(r'STATISTICS_(\w+)" sample="0">([-0-9.eE]+)<', meta)}
            valid = data[data != NODATA].astype(np.float64)
            mine = {'MINIMUM': valid.min(), 'MAXIMUM': valid.max(), 'MEAN': valid.mean(), 'STDDEV': valid.std()} if valid.size else {}
            for key in ('MINIMUM', 'MAXIMUM', 'MEAN', 'STDDEV'):
                if key in stats:
                    assert abs(stats[key] - mine[key]) <= 1e-6 * max(1, abs(stats[key])), (path, key, stats[key], mine[key])
            r0 = (90 - (south + 20)) * 120
            c0 = (west + 180) * 120
            out[r0:r0 + TILE_ROWS, c0:c0 + TILE_COLS] = data
            report.append((lat + lon, os.path.getsize(path), int((data == NODATA).sum()), stats.get('MINIMUM'), stats.get('MAXIMUM'), stats.get('MEAN')))
    out.flush()
    for row in report:
        print(*row)
    print('tiles', len(report), 'bytes', sum(r[1] for r in report), 'nodata', sum(r[2] for r in report))


def load(cache):
    return np.load(os.path.join(cache, 'gmted_mea300.npy'), mmap_mode='r')


def latitude(row, rows=ROWS):
    return 90.0 - (row + 0.5) * 180.0 / rows


def kernel_profile(r, width, edge):
    """IFS Cy47r3 Part IV eq. 11.4, unnormalised (1/Δ dropped)."""
    r = np.abs(r)
    inner, outer = width / 2 - edge, width / 2 + edge
    h = np.where(r < inner, 1.0, np.where(r < outer, 0.5 + 0.5 * np.cos(np.pi * (r - inner) / (2 * edge)), 0.0))
    return h


def band_kernel(phi, width, edge, cols=COLS, rows=ROWS):
    """The kernel on the 30″ grid at latitude φ (degrees): an array of shape (2 m + 1, cols) centred on column 0, unit sum."""
    dy = RADIUS * np.pi / rows
    dx = RADIUS * np.cos(np.radians(phi)) * 2 * np.pi / cols
    reach = width / 2 + edge
    m = int(np.ceil(reach / dy))
    n = min(cols // 2, int(np.ceil(reach / max(dx, 1e-6))))
    yy = np.arange(-m, m + 1)[:, None] * dy
    xx = np.arange(-n, n + 1)[None, :] * dx
    h = kernel_profile(np.hypot(xx, yy), width, edge)
    full = np.zeros((2 * m + 1, cols))
    np.add.at(full, (slice(None), np.arange(-n, n + 1) % cols), h)
    return full / full.sum(), m


def smooth_band(source, r0, r1, width, edge, phi, total_rows=ROWS):
    """Rows r0..r1 of `source` (rows of a global grid of total_rows from the north) smoothed with the kernel at latitude φ."""
    rows, cols = source.shape
    kernel, m = band_kernel(phi, width, edge, cols, total_rows)
    lo, hi = r0 - m, r1 + m
    block = np.zeros((hi - lo, cols))
    a, b = max(0, lo), min(rows, hi)
    block[a - lo:b - lo] = source[a:b]
    if lo < 0:
        block[:-lo] = source[0:1]
    if hi > rows:
        block[b - lo:] = source[rows - 1:rows]
    pad = np.zeros((hi - lo, cols))
    pad[:2 * m + 1] = kernel
    pad = np.roll(pad, -m, axis=0)
    out = np.fft.irfft2(np.fft.rfft2(block) * np.fft.rfft2(pad), s=block.shape)
    return out[m:m + (r1 - r0)]


def clamp(raw):
    h = raw.astype(np.float64)
    h[raw == NODATA] = 0.0
    return np.maximum(h, 0.0)


def filter_fields(source, band_rows=120, block=BLOCK, log=print, row0=0, total_rows=ROWS):
    """The four 2′30″ grids from a 30″ int16 grid (rows from the north, columns from −180°): the global grid, or its rows from row0 on."""
    rows, cols = source.shape
    out_rows, out_cols = rows // block, cols // block
    h5 = np.zeros((out_rows, out_cols), np.float32)
    flt2 = np.zeros((out_rows, out_cols), np.float32)
    mean = np.zeros((out_rows, out_cols), np.float32)
    land = np.zeros((out_rows, out_cols), np.float32)
    start = time.time()
    r0, bands = 0, 0
    while r0 < rows:
        phi = abs(latitude(row0 + r0, total_rows))
        span = 0.02 / max(np.tan(np.radians(min(phi, 89.99))), 1e-9) / (np.pi / total_rows)
        r1 = min(rows, r0 + max(block, min(band_rows, int(span) // block * block)))
        lo, hi = max(0, r0 - 40), min(rows, r1 + 40)
        raw = np.asarray(source[lo:hi])
        local = clamp(raw)
        middle = latitude(row0 + 0.5 * (r0 + r1 - 1), total_rows)
        s5 = smooth_band(local, r0 - lo, r1 - lo, 5000.0, 1000.0, middle, total_rows)
        s2 = smooth_band(local, r0 - lo, r1 - lo, 2000.0, 1000.0, middle, total_rows)
        s20 = smooth_band(local, r0 - lo, r1 - lo, 20000.0, 1000.0, middle, total_rows)
        here = local[r0 - lo:r1 - lo]
        o0, o1 = r0 // block, r1 // block
        centre = block // 2
        h5[o0:o1] = s5[centre::block, centre::block]
        diff = (s2 - s20) ** 2
        flt2[o0:o1] = diff.reshape(o1 - o0, block, out_cols, block).mean(axis=(1, 3))
        mean[o0:o1] = here.reshape(o1 - o0, block, out_cols, block).mean(axis=(1, 3))
        land[o0:o1] = (raw[r0 - lo:r1 - lo] > 0).reshape(o1 - o0, block, out_cols, block).mean(axis=(1, 3))
        if log and bands % 20 == 0:
            log(f'rows {r0}-{r1} of {rows}, {time.time() - start:.0f} s')
        r0, bands = r1, bands + 1
    return h5, flt2, mean, land


def write_fields(cache, fields):
    for name, grid in zip(('h5', 'flt2', 'mean', 'land'), fields):
        path = os.path.join(cache, f'fine_{name}.f32')
        grid.astype('<f4').tofile(path)
        print(name, path, grid.shape, os.path.getsize(path))


def spectrum(cache, segment=1024):
    """Power spectra of land orography along rows (east–west) and columns (north–south), segments wholly on land."""
    source = load(cache)
    results = {}
    window = np.hanning(segment)
    for axis in ('ew', 'ns'):
        total = np.zeros(segment // 2 + 1)
        count = 0
        spacing_sum = 0.0
        for r in range(2400, ROWS - 2400, 120):
            phi = latitude(r)
            if abs(phi) > 60:
                continue
            if axis == 'ew':
                line = np.asarray(source[r]).astype(np.float64)
                dx = RADIUS * np.cos(np.radians(phi)) * 2 * np.pi / COLS
                for c in range(0, COLS - segment, segment):
                    seg = line[c:c + segment]
                    if (seg <= 0).any() or (seg == NODATA).any():
                        continue
                    seg = seg - np.polyval(np.polyfit(np.arange(segment), seg, 1), np.arange(segment))
                    total += np.abs(np.fft.rfft(seg * window)) ** 2 * dx / segment
                    spacing_sum += dx
                    count += 1
            else:
                if r + segment > ROWS:
                    continue
                block = np.asarray(source[r:r + segment:1, ::120]).astype(np.float64)
                dy = RADIUS * np.pi / ROWS
                for j in range(block.shape[1]):
                    seg = block[:, j]
                    if (seg <= 0).any() or (seg == NODATA).any():
                        continue
                    seg = seg - np.polyval(np.polyfit(np.arange(segment), seg, 1), np.arange(segment))
                    total += np.abs(np.fft.rfft(seg * window)) ** 2 * dy / segment
                    spacing_sum += dy
                    count += 1
        spacing = spacing_sum / count
        k = 2 * np.pi * np.fft.rfftfreq(segment, d=spacing)
        power = total / count
        results[axis] = (k, power, count, spacing)
        for lo, hi, label in ((0.000628, 0.003, 'k0-k1 (10-2.1 km)'), (0.00014, 0.00112, 'filter band (45-5.6 km)'), (0.0001, 0.0006, '63-10 km')):
            sel = (k >= lo) & (k <= hi)
            slope = np.polyfit(np.log(k[sel]), np.log(power[sel]), 1)[0]
            print(f'{axis} {label}: slope {slope:.3f} over {sel.sum()} wavenumbers, {count} segments, spacing {spacing:.0f} m')
    return results


def kernel_response(k, width, edge):
    """The continuous kernel's response to a plane wave of wavenumber k (rad/m): ∫ h J₀(kr) r dr / ∫ h r dr."""
    r = np.linspace(0, width / 2 + edge, 20001)[1:]
    t = np.linspace(0, np.pi, 801)
    j0 = np.trapezoid(np.cos(np.outer(k * r, np.sin(t))), t, axis=1) / np.pi
    h = kernel_profile(r, width, edge)
    return np.trapezoid(h * j0 * r, r) / np.trapezoid(h * r, r)


def selftest(amplitude=300.0, base=1000.0):
    """Analytic orography through `filter_fields` on a 5° band at the equator: a ridge east–west at 40 and 8 km and an isotropic field at 40 km."""
    rows, row0 = 600, ROWS // 2 - 300
    phi = latitude(row0 + np.arange(rows), ROWS)[:, None]
    lon = np.radians(-180 + (np.arange(COLS) + 0.5) / 120)[None, :]
    x = RADIUS * np.cos(np.radians(phi)) * lon
    y = RADIUS * np.radians(phi)
    out = {}
    for name, wavelength, iso in (('ridge40', 40000.0, False), ('ridge8', 8000.0, False), ('iso40', 40000.0, True)):
        waves = round(2 * np.pi * RADIUS / wavelength)
        k = waves / RADIUS
        kx = waves * lon
        field = np.sin(kx) * (np.sin(k * y) if iso else np.ones_like(y))
        source = np.round(base + amplitude * field).astype(np.int16)
        h5, flt2, mean, land = filter_fields(source, log=None, row0=row0, total_rows=ROWS)
        inner = slice(30, 90)
        r5 = kernel_response(k * (np.sqrt(2) if iso else 1), 5000.0, 1000.0)
        band = kernel_response(k * (np.sqrt(2) if iso else 1), 2000.0, 1000.0) - kernel_response(k * (np.sqrt(2) if iso else 1), 20000.0, 1000.0)
        variance = 0.25 if iso else 0.5
        out[name] = {
            'k': k, 'h5Mean': float(h5[inner].mean()), 'h5Std': float(h5[inner].std()), 'h5StdAnalytic': float(amplitude * abs(r5) * np.sqrt(variance)),
            'sigmaFlt': float(np.sqrt(flt2[inner].mean())), 'sigmaFltAnalytic': float(amplitude * abs(band) * np.sqrt(variance)),
            'mean': float(mean[inner].mean()), 'land': float(land[inner].mean()),
        }
    return out


def spectrum_selftest(slope=-1.9, size=1 << 15, spacing=926.6, seed=1):
    """A synthetic line with a power-law spectrum through the slope fit of `spectrum`."""
    rng = np.random.default_rng(seed)
    k = 2 * np.pi * np.fft.rfftfreq(size, d=spacing)
    amplitude = np.zeros_like(k)
    amplitude[1:] = k[1:] ** (slope / 2)
    line = np.fft.irfft(amplitude * np.exp(2j * np.pi * rng.random(k.size)) * rng.standard_normal(k.size), n=size)
    segment, window, total, count = 1024, np.hanning(1024), np.zeros(513), 0
    for c in range(0, size - segment, segment):
        seg = line[c:c + segment]
        seg = seg - np.polyval(np.polyfit(np.arange(segment), seg, 1), np.arange(segment))
        total += np.abs(np.fft.rfft(seg * window)) ** 2
        count += 1
    kk = 2 * np.pi * np.fft.rfftfreq(segment, d=spacing)
    sel = (kk >= 0.000628) & (kk <= 0.003)
    return float(np.polyfit(np.log(kk[sel]), np.log(total[sel] / count), 1)[0])


if __name__ == '__main__':
    if sys.argv[1] == 'selftest':
        import json
        print(json.dumps({'filters': selftest(), 'spectrumSlope': spectrum_selftest()}))
        raise SystemExit(0)
    command, cache = sys.argv[1], sys.argv[2]
    if command == 'assemble':
        assemble(cache)
    elif command == 'filter':
        write_fields(cache, filter_fields(load(cache)))
    elif command == 'spectrum':
        spectrum(cache)
    else:
        raise SystemExit(f'unknown command {command}')
