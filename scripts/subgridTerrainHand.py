#!/usr/bin/env python3
"""
An independent check of data/subgrid_N<N>.bin for a few cells, from the
30″ grid of `subgridTerrain.py assemble` directly: every smoothing as a
direct sum over the 30″ points within the kernel's reach, with the
kernel laid out at each target row's own latitude; the cell's members by
brute force against its two rings of neighbours; the resolved orography
on the triangle of cell centres that contains each point; σ_flt over the
cell's 30″ points (the IFS's definition) rather than the 2′30″ blocks.

  node scripts/subgridTerrainCells.mjs 64 39,-99 28,84 -32.6,-70 > cells.json
  python3 scripts/subgridTerrainHand.py CACHE cells.json
"""
import json
import sys
import numpy as np

RADIUS = 6371220.0
ROWS, COLS = 21600, 43200


def profile(r, width, edge):
    inner, outer = width / 2 - edge, width / 2 + edge
    return np.where(r < inner, 1.0, np.where(r < outer, 0.5 + 0.5 * np.cos(np.pi * (r - inner) / (2 * edge)), 0.0))


def smooth_rows(window, rows, width, edge):
    """Each listed row of `window` (30″ rows, `rows` their global indices, whole latitude circles) smoothed by direct sums."""
    dy = RADIUS * np.pi / ROWS
    reach = width / 2 + edge
    m = int(np.ceil(reach / dy))
    out = {}
    for row in rows:
        phi = 90.0 - (row + 0.5) / 120
        dx = RADIUS * np.cos(np.radians(phi)) * 2 * np.pi / COLS
        n = int(np.ceil(reach / dx))
        total = np.zeros(window.shape[1])
        weight = 0.0
        for j in range(-m, m + 1):
            w = profile(np.hypot(np.arange(-n, n + 1) * dx, j * dy), width, edge)
            line = window.row(row + j)
            padded = np.concatenate([line[-n:], line, line[:n]])
            total += np.convolve(padded, w[::-1], mode='valid')
            weight += w.sum()
        out[row] = total / weight
    return out


class Rows:
    def __init__(self, grid):
        self.grid = grid
        self.cache = {}

    def row(self, r):
        r = min(ROWS - 1, max(0, r))
        if r not in self.cache:
            raw = np.asarray(self.grid[r]).astype(np.float64)
            raw[raw == -32768] = 0
            self.cache[r] = np.maximum(raw, 0)
        return self.cache[r]

    @property
    def shape(self):
        return (ROWS, COLS)


def xyz(lat, lon):
    lat, lon = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)


def resolved(points, triangles):
    out = np.zeros(len(points))
    for n, p in enumerate(points):
        best, value = -np.inf, 0.0
        for t in triangles:
            a, b, c = (np.array(v['x']) for v in t)
            weights = np.linalg.solve(np.stack([a, b, c], axis=1), p)
            low = weights.min() / weights.sum()
            if low > best:
                best, value = low, weights @ np.array([v['h'] for v in t]) / weights.sum()
        out[n] = value
    return out


def check(grid, cell):
    rows = Rows(grid)
    centre = np.array(cell['centre'])
    ring = np.array([r['x'] for r in cell['ring']])
    own = [r['i'] for r in cell['ring']].index(cell['cell'])
    lat0, lon0 = cell['lat'], cell['lon']
    half = 1.6
    fine_rows = [j for j in range(int((90 - lat0 - half) * 24), int((90 - lat0 + half) * 24) + 1)]
    span = half / max(np.cos(np.radians(lat0 + np.sign(lat0) * half)), 0.2)
    fine_cols = [(i % 8640) for i in range(int((lon0 - span + 180) * 24), int((lon0 + span + 180) * 24) + 1)]
    flat = lambda j: 90 - (j + 0.5) / 24
    flon = lambda i: -180 + (i + 0.5) / 24
    members = []
    for j in fine_rows:
        for i in fine_cols:
            p = xyz(flat(j), flon(i))
            if np.argmax(ring @ p) == own:
                members.append((j, i))
    need_rows = sorted({5 * j + 2 + d for j, _ in members for d in (-5, 0, 5)})
    s5 = smooth_rows(rows, need_rows, 5000.0, 1000.0)
    h5 = lambda j, i: s5[5 * j + 2][(5 * i + 2) % COLS]
    stencil = []
    for j, i in members:
        stencil += [(j, i), (j, (i + 1) % 8640), (j, (i - 1) % 8640), (j - 1, i), (j + 1, i)]
    stencil = sorted(set(stencil))
    index = {key: n for n, key in enumerate(stencil)}
    points = np.array([xyz(flat(j), flon(i)) for j, i in stencil])
    res = np.array([h5(j, i) for j, i in stencil]) - resolved(points, cell['triangles'])
    w = s = s2 = K = L = M = 0.0
    dlat = np.radians(1 / 24)
    for j, i in members:
        phi = np.radians(flat(j))
        h = res[index[(j, i)]]
        hx = (res[index[(j, (i + 1) % 8640)]] - res[index[(j, (i - 1) % 8640)]]) / (2 * RADIUS * np.cos(phi) * dlat)
        hy = (res[index[(j - 1, i)]] - res[index[(j + 1, i)]]) / (2 * RADIUS * dlat)
        c = np.cos(phi)
        w += c; s += c * h; s2 += c * h * h
        K += c * 0.5 * (hx * hx + hy * hy); L += c * 0.5 * (hx * hx - hy * hy); M += c * hx * hy
    K, L, M = K / w, L / w, M / w
    spread = np.hypot(L, M)
    out = {
        'members': len(members),
        'deviation': float(np.sqrt(s2 / w - (s / w) ** 2)),
        'anisotropy': float(np.sqrt((K - spread) / (K + spread))),
        'orientation': float(0.5 * np.arctan2(M, L)),
        'slope': float(np.sqrt(K + spread)),
    }
    coarse_rows = sorted({r for j, _ in members for r in range(5 * j, 5 * j + 5)})
    a = smooth_rows(rows, coarse_rows, 2000.0, 1000.0)
    b = smooth_rows(rows, coarse_rows, 20000.0, 1000.0)
    total = weight = 0.0
    count = 0
    cols = np.array(sorted({c for _, i in members for c in range(5 * i, 5 * i + 5)}))
    for r in coarse_rows:
        phi = 90.0 - (r + 0.5) / 120
        lons = -180 + (cols + 0.5) / 120
        p = xyz(np.full(cols.size, phi), lons)
        inside = np.argmax(p @ ring.T, axis=1) == own
        d = (a[r][cols % COLS] - b[r][cols % COLS])[inside]
        total += np.cos(np.radians(phi)) * (d * d).sum()
        weight += np.cos(np.radians(phi)) * inside.sum()
        count += inside.sum()
    out['filtered'] = float(np.sqrt(total / weight))
    out['points30'] = int(count)
    return out


if __name__ == '__main__':
    grid = np.load(f'{sys.argv[1]}/gmted_mea300.npy', mmap_mode='r')
    cells = json.load(open(sys.argv[2]))
    for cell in cells:
        hand = check(grid, cell)
        print(json.dumps({'cell': cell['cell'], 'lat': round(cell['lat'], 2), 'lon': round(cell['lon'], 2), 'hand': hand, 'file': cell['file']}))
