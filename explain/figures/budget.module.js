import { Grid } from '../../js/grid.module.js';
import { buildMesh } from '../../js/mesh.module.js';
import { createSigmaCore } from '../../js/dynamics/sigmaCore.module.js';
import { curl, kineticEnergy, gradient, laplacianVelocity } from '../../js/dynamics/operators.module.js';
import { cellVelocity } from '../swCases.module.js';
import { Figure, choice, legend, caption, readout, text, arrow, rampRGB, termColor, paletteVersion, clamp, DARK_NEUTRAL, INK, MUTED, GRID, LINE } from '../runtime.module.js';

const DEG = Math.PI / 180, HOUR = 3600, WIDE = 520, SIDE = 190, PANEL = 170, LON_SPAN = 100, LAT_SPAN = 60, BOX_LAT = 27.5, MAGNIFY = 5, STEP = 3, REACH = 15, KEY = 10, DIM = 0.62, BG = [24, 24, 26];
const LAYER_SIGMAS = [0.85, 0.5, 0.25], TERMS = ['a', 'b', 'c', 'd'], HALO = 'rgba(20, 20, 22, 0.85)', DARK = `rgb(${DARK_NEUTRAL.join(', ')})`;
const SHARED = ['flux', 'dissipation', 'divFlux', 'piSigmaDot', 'exnerLower', 'exnerLayer', 'dExnerDpi', 'thetaLower', 'qLower', 'qcLower', 'thetaV', 'geopotential', 'piVertex'];
const NAMES = { a: 'turning', b: 'height slope', c: 'pressure change', d: 'up and down' };
const VIEWS = {
  a: 'a, the turning: −(ζ + f) ẑ × u', bc: 'b + c, the pressure gradient: −∇(Φ + K) − R T ∇ln π', sum: 'a + b + c + d, the acceleration ∂u/∂t',
  b: 'b, the height slope: −∇(Φ + K)', c: 'c, the pressure change: −R T ∇ln π', d: 'd, up and down: −σ̇ ∂u/∂σ',
};

const states = new Map();
export function loadState(url) {
  if (!states.has(url)) states.set(url, fetch(url).then((r) => r.arrayBuffer()).then((bytes) => {
    if (String.fromCharCode(...new Uint8Array(bytes, 0, 4)) !== 'EXS1') throw new Error(`${url} is not a state file`);
    const length = new DataView(bytes).getUint32(4, true), header = JSON.parse(new TextDecoder().decode(new Uint8Array(bytes, 8, length))), data = new Float32Array(bytes, 8 + length);
    return { header, arrays: Object.fromEntries(header.arrays.map((a) => [a.name, data.subarray(a.offset, a.offset + a.length)])) };
  }));
  return states.get(url);
}

export function momentumBudget({ header, arrays }, sigmas = LAYER_SIGMAS) {
  const { N, a, omega, g, cp, R, p0, nu4, nu4Theta } = header, levels = Float64Array.from(header.levels), K = levels.length - 1;
  const mesh = buildMesh(new Grid(N), { radius: a, omega }), { nCells: C, nEdges: E, nVertices: V, maxEdgesOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge, cellsOnEdge, verticesOnEdge, dcEdge, dvEdge, fVertex } = mesh;
  const sizes = { flux: K * E, dissipation: K * E, piSigmaDot: (K + 1) * C, piVertex: V }, buffers = Object.fromEntries(SHARED.map((name) => [name, new ArrayBuffer(8 * (sizes[name] ?? K * C))]));
  const core = createSigmaCore(mesh, { levels, g, cp, R, p0, nu4, nu4Theta, surfaceGeopotential: Float64Array.from(arrays.surfaceGeopotential), buffers });
  const state = [Float64Array.from(arrays.pi), Float64Array.from(arrays.theta), Float64Array.from(arrays.u)], [pi, , u] = state;
  const tendency = [new Float64Array(C), new Float64Array(K * C), new Float64Array(K * E)];
  core.tendency(state, tendency);
  const { geopotential, exnerLayer, thetaV, piSigmaDot } = core.arrays, dSigma = core.diagnostics.dSigma, sigmaMid = core.sigmaMid;
  const flux = new Float64Array(buffers.flux), piVertex = new Float64Array(buffers.piVertex);
  const piEdge = Float64Array.from({ length: E }, (_, e) => 0.5 * (pi[cellsOnEdge[2 * e]] + pi[cellsOnEdge[2 * e + 1]])), gradLnPi = gradient(mesh, Float64Array.from(pi, Math.log));
  const layers = sigmas.map((target) => {
    let k = 0;
    for (let m = 1; m < K; m++) if (Math.abs(sigmaMid[m] - target) < Math.abs(sigmaMid[k] - target)) k = m;
    const off = k * C, uk = u.subarray(k * E, (k + 1) * E), fk = flux.subarray(k * E, (k + 1) * E), zeta = curl(mesh, uk);
    const qVertex = Float64Array.from(zeta, (z, v) => (z + fVertex[v]) / piVertex[v]), qEdge = Float64Array.from({ length: E }, (_, e) => 0.5 * (qVertex[verticesOnEdge[2 * e]] + qVertex[verticesOnEdge[2 * e + 1]]));
    const kinetic = kineticEnergy(mesh, uk), gradPhi = gradient(mesh, Float64Array.from(kinetic, (x, i) => geopotential[off + i] + x));
    const terms = { a: new Float64Array(E), b: new Float64Array(E), c: new Float64Array(E), d: new Float64Array(E) };
    for (let e = 0; e < E; e++) {
      const i = cellsOnEdge[2 * e], j = cellsOnEdge[2 * e + 1];
      let pv = 0;
      for (let s = 0; s < nEdgesOnEdge[e]; s++) { const slot = maxEdgesOnEdge * e + s, other = edgesOnEdge[slot]; pv += weightsOnEdge[slot] * dvEdge[other] * fk[other] * (0.5 * qEdge[e] + 0.5 * qEdge[other]); }
      const lowerFlow = 0.5 * (piSigmaDot[(k + 1) * C + i] + piSigmaDot[(k + 1) * C + j]), upperFlow = 0.5 * (piSigmaDot[k * C + i] + piSigmaDot[k * C + j]);
      const lowerU = k === K - 1 ? 0 : 0.5 * (uk[e] + u[(k + 1) * E + e]), upperU = k === 0 ? 0 : 0.5 * (uk[e] + u[(k - 1) * E + e]);
      terms.a[e] = pv / dcEdge[e];
      terms.b[e] = -gradPhi[e];
      terms.c[e] = -R * 0.5 * (thetaV[off + i] * exnerLayer[off + i] + thetaV[off + j] * exnerLayer[off + j]) * gradLnPi[e];
      terms.d[e] = -(lowerFlow * lowerU - upperFlow * upperU - uk[e] * (lowerFlow - upperFlow)) / (piEdge[e] * dSigma[k]);
    }
    terms.closure = Float64Array.from(laplacianVelocity(mesh, laplacianVelocity(mesh, uk)), (x) => -nu4 * x);
    return { k, sigma: sigmaMid[k], terms, du: tendency[2].subarray(k * E, (k + 1) * E), geopotential: geopotential.subarray(off, off + C) };
  });
  return { header, mesh, core, K, pi, u, layers };
}

const rgbOf = (color) => (/^#[0-9a-f]{6}$/i.test(color) ? [1, 3, 5].map((n) => parseInt(color.slice(n, n + 2), 16)) : [200, 200, 200]);
const css = (c, alpha = 1) => `rgba(${Math.round(c[0])}, ${Math.round(c[1])}, ${Math.round(c[2])}, ${alpha})`;
const size = (v) => (v < 0.01 ? v.toPrecision(1) : v < 1 ? v.toFixed(2) : v.toFixed(1));
const latName = (lat) => `${Math.abs(Math.round(lat))}°${lat >= 0 ? 'N' : 'S'}`;
const lonName = (lon) => { const l = ((Math.round(lon) % 360) + 540) % 360 - 180; return l === -180 || l === 180 ? '180°' : `${Math.abs(l)}°${l > 0 ? 'E' : 'W'}`; };
const wrap = (x) => Math.atan2(Math.sin(x), Math.cos(x));
const nice = (v) => [1, 2, 5, 10, 20, 50].find((n) => n >= v) ?? 100;

function locator(mesh) {
  const { maxEdges, nEdgesOnCell, cellsOnCell, verticesOnCell, cellsOnVertex, xCell } = mesh;
  const dot = (x, y, z, i) => x * xCell[3 * i] + y * xCell[3 * i + 1] + z * xCell[3 * i + 2];
  const triple = (x, y, z, i, j) => { const ax = xCell[3 * i], ay = xCell[3 * i + 1], az = xCell[3 * i + 2], bx = xCell[3 * j], by = xCell[3 * j + 1], bz = xCell[3 * j + 2]; return x * (ay * bz - az * by) + y * (az * bx - ax * bz) + z * (ax * by - ay * bx); };
  let cell = 0;
  return (x, y, z, cells, weights, slot) => {
    for (let best = dot(x, y, z, cell), moved = true; moved;) {
      moved = false;
      for (let m = 0, from = cell; m < nEdgesOnCell[from]; m++) { const j = cellsOnCell[maxEdges * from + m], d = dot(x, y, z, j); if (d > best) { best = d; cell = j; moved = true; } }
    }
    let bestMin = -Infinity;
    for (let m = 0; m < nEdgesOnCell[cell]; m++) {
      const v = verticesOnCell[maxEdges * cell + m], i = cellsOnVertex[3 * v], j = cellsOnVertex[3 * v + 1], k = cellsOnVertex[3 * v + 2];
      const wi = triple(x, y, z, j, k), wj = triple(x, y, z, k, i), wk = triple(x, y, z, i, j), total = wi + wj + wk, least = Math.min(wi, wj, wk) / total;
      if (least > bestMin) { bestMin = least; cells[3 * slot] = i; cells[3 * slot + 1] = j; cells[3 * slot + 2] = k; weights[3 * slot] = wi / total; weights[3 * slot + 1] = wj / total; weights[3 * slot + 2] = wk / total; }
    }
  };
}

function contours(values, gw, gh, interval) {
  const segments = [], p = [0, 0, 0, 0], cx = [0, 1, 1, 0], cy = [0, 0, 1, 1], cross = [];
  for (let gy = 0; gy < gh - 1; gy++) for (let gx = 0; gx < gw - 1; gx++) {
    p[0] = values[gy * gw + gx]; p[1] = values[gy * gw + gx + 1]; p[2] = values[(gy + 1) * gw + gx + 1]; p[3] = values[(gy + 1) * gw + gx];
    for (let level = Math.ceil(Math.min(...p) / interval) * interval; level < Math.max(...p); level += interval) {
      cross.length = 0;
      for (let m = 0; m < 4; m++) { const n = (m + 1) % 4, a = p[m], b = p[n]; if ((a < level) !== (b < level)) { const t = (level - a) / (b - a); cross.push((gx + cx[m] + t * (cx[n] - cx[m])) * STEP, (gy + cy[m] + t * (cy[n] - cy[m])) * STEP); } }
      segments.push(...cross);
    }
  }
  return segments;
}

export function mountBudget(root) {
  const controls = root.querySelector('.controls');
  const color = { a: termColor('a'), b: termColor('b'), c: termColor('c'), d: termColor('d'), sum: INK }, pressure = rgbOf(color.b).map((x, n) => (x + rgbOf(color.c)[n]) / 2);
  color.bc = css(pressure);
  const faint = { a: css(rgbOf(color.a), 0.55), bc: css(pressure, 0.55) };
  let model = null, layer = 0, view = 'sum', probe = -1, geo = null, painted = null, failed = false;

  const fig = new Figure(root, { height: 420, minHeight: 380, draw });
  const fit = fig.resize.bind(fig);
  fig.resize = () => { const w = fig.stage.clientWidth, narrow = w < WIDE, tall = Math.min(320, Math.round(w * 0.85)) + PANEL; fig.height = narrow ? tall : 420; fig.minHeight = narrow ? tall : 380; fit(); };
  caption(root, 'The baroclinic wave of Part 3 on day 9, as the model holds it: 27 layers on the N = 16 grid. Each arrow is a term of the momentum equation at one cell, computed by the model’s own code, and the terms share one scale so their lengths compare. Tap or hover over a cell, or drag the ring, to take it apart.');
  const key = legend(root, []);
  choice(controls, { label: 'Layer', options: [['near 850 hPa', '0'], ['near 500 hPa', '1'], ['near 250 hPa', '2']], value: '0', onChange: (v) => { layer = Number(v); setKey(); fig.render(); }, span: true });
  choice(controls, { label: 'Show', options: [['the turning, a', 'a'], ['the pressure terms, b + c', 'bc'], ['all four: the acceleration', 'sum'], ['b alone', 'b'], ['c alone', 'c'], ['d alone', 'd']], value: view, onChange: (v) => { view = v; setKey(); fig.render(); }, span: true });
  const out = readout(controls);

  function setKey() {
    const L = model?.layers[layer];
    key.set([
      ['ramp', 'surface pressure, 970 to 1030 hPa', 'cool', 'warm', DARK], ['faint', L ? `lines of equal pressure ${(L.height / 1000).toFixed(1)} km up, about the layer’s height, every ${L.interval} hPa` : 'lines of equal pressure at about the layer’s height', LINE],
      ['force', view === 'sum' ? `${VIEWS.sum}, drawn ${MAGNIFY} times longer than the terms` : VIEWS[view], color[view]], ...(view === 'sum' ? [['force', 'the turning a, faint', faint.a], ['force', 'the pressure terms b + c, faint', faint.bc]] : []),
    ]);
  }
  setKey();

  function prepare(m) {
    const { mesh, pi, u, layers, header, K, core } = m, C = mesh.nCells, v = [0, 0], column = core.arrays.geopotential, sigma = core.sigmaMid;
    let low = 0;
    for (let i = 0; i < C; i++) if (pi[i] < pi[low]) low = i;
    const lat0 = mesh.latCell[low], lon0 = mesh.lonCell[low];
    const box = Array.from({ length: C }, (_, i) => i).filter((i) => Math.abs(wrap(mesh.lonCell[i] - lon0)) <= LON_SPAN / 2 * DEG && Math.abs(mesh.latCell[i] - lat0) <= BOX_LAT * DEG);
    const median = (values) => values.sort((p, q) => p - q)[Math.floor((values.length - 1) / 2)];
    const vectors = (field, offset = 0) => { const out = new Float64Array(2 * C); for (let i = 0; i < C; i++) { cellVelocity(mesh, field, i, v, 0, offset); out[2 * i] = v[0]; out[2 * i + 1] = v[1]; } return out; };
    for (const L of layers) {
      L.cell = Object.fromEntries([...TERMS, 'closure'].map((t) => [t, vectors(L.terms[t])]));
      L.cell.bc = L.cell.b.map((x, n) => x + L.cell.c[n]);
      L.cell.sum = L.cell.a.map((x, n) => x + L.cell.b[n] + L.cell.c[n] + L.cell.d[n]);
      L.wind = vectors(u, L.k * mesh.nEdges);
      L.median = median(box.map((i) => Math.hypot(L.cell.sum[2 * i], L.cell.sum[2 * i + 1]) / Math.hypot(L.cell.bc[2 * i], L.cell.bc[2 * i + 1])));
      L.height = Math.round(median(box.map((i) => L.geopotential[i] / header.g)) / 100) * 100;
      const target = header.g * L.height;
      L.pressure = Float64Array.from({ length: C }, (_, i) => {
        let k = K - 2;
        while (k > 0 && column[k * C + i] < target) k--;
        const w = (target - column[(k + 1) * C + i]) / (column[k * C + i] - column[(k + 1) * C + i]);
        return pi[i] / 100 * Math.exp(Math.log(sigma[k + 1]) + w * (Math.log(sigma[k]) - Math.log(sigma[k + 1])));
      });
      const span = Math.max(...box.map((i) => L.pressure[i])) - Math.min(...box.map((i) => L.pressure[i]));
      L.interval = span / 5 <= 18 ? 5 : 10;
    }
    let start = low, most = 0;
    for (const i of box) {
      const d = Math.acos(clamp(mesh.xCell[3 * i] * mesh.xCell[3 * low] + mesh.xCell[3 * i + 1] * mesh.xCell[3 * low + 1] + mesh.xCell[3 * i + 2] * mesh.xCell[3 * low + 2], -1, 1));
      const s = Math.hypot(layers[0].cell.a[2 * i], layers[0].cell.a[2 * i + 1]);
      if (d < 12 * DEG && s > most) { most = s; start = i; }
    }
    let spacing = 0;
    for (let e = 0; e < mesh.nEdges; e++) spacing += mesh.dcEdge[e];
    const extrema = [];
    for (const i of box) {
      const around = Array.from({ length: mesh.nEdgesOnCell[i] }, (_, n) => pi[mesh.cellsOnCell[mesh.maxEdges * i + n]]);
      if (around.every((p) => p > pi[i])) extrema.push({ i, kind: 'L' });
      else if (around.every((p) => p < pi[i])) extrema.push({ i, kind: 'H' });
    }
    return { ...m, low, lat0, lon0, probe: start, spacing: spacing / mesh.nEdges / header.a, extrema, locate: locator(mesh) };
  }

  function layout(w, h) {
    const narrow = w < WIDE, mapW = narrow ? w : w - SIDE, mapH = narrow ? h - PANEL : h, { lat0, lon0, mesh } = model;
    const s = Math.max(mapW / (LON_SPAN * DEG * Math.cos(lat0)), mapH / (LAT_SPAN * DEG)), cx = mapW / 2, cy = mapH / 2;
    const toX = (lon) => cx + wrap(lon - lon0) * Math.cos(lat0) * s, toY = (lat) => cy - (lat - lat0) * s;
    const cells = [];
    for (let i = 0; i < mesh.nCells; i++) { const x = toX(mesh.lonCell[i]), y = toY(mesh.latCell[i]); if (Math.abs(wrap(mesh.lonCell[i] - lon0)) < Math.PI / 2 && x > 3 && x < mapW - 3 && y > 3 && y < mapH - 3) cells.push([i, x, y]); }
    const gw = Math.ceil(mapW / STEP) + 1, gh = Math.ceil(mapH / STEP) + 1, n = gw * gh, tri = new Int32Array(3 * n), weights = new Float64Array(3 * n);
    for (let gy = 0; gy < gh; gy++) for (let gx = 0; gx < gw; gx++) {
      const lon = lon0 + (gx * STEP - cx) / (s * Math.cos(lat0)), lat = lat0 - (gy * STEP - cy) / s, c = Math.cos(lat);
      model.locate(c * Math.cos(lon), c * Math.sin(lon), Math.sin(lat), tri, weights, gy * gw + gx);
    }
    const sample = (field, scale = 1) => Float64Array.from({ length: n }, (_, m) => (weights[3 * m] * field[tri[3 * m]] + weights[3 * m + 1] * field[tri[3 * m + 1]] + weights[3 * m + 2] * field[tri[3 * m + 2]]) * scale);
    const pxPer = model.spacing * s / REACH * HOUR, square = narrow ? Math.min(PANEL, Math.round(w / 2)) : SIDE;
    const panel = narrow ? { diagram: [0, mapH, square, h], notes: [square, mapH, w, h] } : { diagram: [mapW, 0, w, square + 16], notes: [mapW, square + 16, w, h] };
    return { w, h, narrow, mapW, mapH, toX, toY, cells, gw, gh, surface: sample(model.pi, 0.01), lines: model.layers.map((L) => contours(sample(L.pressure), gw, gh, L.interval)), image: null, pxPer, panel, key: `${w}x${h}` };
  }

  function paint() {
    const image = document.createElement('canvas'), c = image.getContext('2d'), data = c.createImageData(geo.gw, geo.gh), rgb = [0, 0, 0];
    image.width = geo.gw; image.height = geo.gh;
    for (let n = 0; n < geo.surface.length; n++) {
      rampRGB(0.5 + (geo.surface[n] - 1000) / 60, rgb, DARK_NEUTRAL);
      for (let m = 0; m < 3; m++) data.data[4 * n + m] = BG[m] + (rgb[m] * 255 - BG[m]) * DIM;
      data.data[4 * n + 3] = 255;
    }
    c.putImageData(data, 0, 0);
    geo.image = image;
    painted = paletteVersion;
  }

  const vec = (name, i) => { const f = model.layers[layer].cell[name]; return [f[2 * i], f[2 * i + 1]]; };

  function map(ctx) {
    const { mapW, mapH, toX, toY, cells, pxPer } = geo, { lat0, lon0 } = model, lines = geo.lines[layer];
    ctx.save(); ctx.beginPath(); ctx.rect(0, 0, mapW, mapH); ctx.clip();
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(geo.image, -STEP / 2, -STEP / 2, geo.gw * STEP, geo.gh * STEP);
    ctx.strokeStyle = GRID; ctx.lineWidth = 1; ctx.beginPath();
    const latLines = [], lonLines = [];
    for (let lat = Math.ceil((lat0 / DEG - 40) / 10) * 10; lat <= lat0 / DEG + 40; lat += 10) { const y = toY(lat * DEG); if (y > 18 && y < mapH - 20) { ctx.moveTo(0, y); ctx.lineTo(mapW, y); latLines.push([lat, y]); } }
    for (let lon = Math.ceil((lon0 / DEG - 70) / 20) * 20; lon <= lon0 / DEG + 70; lon += 20) { const x = toX(lon * DEG); if (x > 16 && x < mapW - 16) { ctx.moveTo(x, 0); ctx.lineTo(x, mapH); lonLines.push([lon, x]); } }
    ctx.stroke();
    ctx.strokeStyle = LINE; ctx.beginPath();
    for (let n = 0; n < lines.length; n += 4) { ctx.moveTo(lines[n], lines[n + 1]); ctx.lineTo(lines[n + 2], lines[n + 3]); }
    ctx.stroke();
    for (const [lat, y] of latLines) text(ctx, latName(lat), 5, y - 7, { color: MUTED, size: 10, halo: HALO });
    for (const [lon, x] of lonLines) text(ctx, lonName(lon), x + 4, mapH - 9, { color: MUTED, size: 10, halo: HALO });
    const L = model.layers[layer], sum = view === 'sum', shown = sum ? [['a', true, 1], ['bc', true, 1], ['sum', false, MAGNIFY]] : [[view, false, 1]];
    ctx.fillStyle = 'rgba(255,255,255,0.45)';
    for (const [, x, y] of cells) { ctx.beginPath(); ctx.arc(x, y, 1.2, 0, 2 * Math.PI); ctx.fill(); }
    let longest = 0;
    for (const [name, dim, times] of shown) {
      const f = L.cell[name], px = pxPer * times;
      for (const [i, x, y] of cells) {
        arrow(ctx, x, y, x + f[2 * i] * px, y - f[2 * i + 1] * px, { color: dim ? faint[name] : color[name], width: dim ? 1.2 : 1.5, head: 5, dash: [4, 3], open: true });
        if (!dim) longest = Math.max(longest, Math.hypot(f[2 * i], f[2 * i + 1]));
      }
    }
    if (longest * pxPer * (sum ? MAGNIFY : 1) < 3) text(ctx, `${sum ? 'white arrows' : 'arrows'} too short to see: at most ${size(longest * HOUR)} m/s per hour`, mapW / 2, mapH - 30, { color: INK, size: 11, align: 'center', halo: HALO });
    for (const { i, kind } of model.extrema) {
      const x = toX(model.mesh.lonCell[i]), y = toY(model.mesh.latCell[i]);
      if (x < 14 || x > mapW - 14 || y < 14 || y > mapH - 24) continue;
      text(ctx, kind, x, y - 11, { color: INK, size: 14, weight: 700, align: 'center', halo: HALO });
      text(ctx, `${Math.round(model.pi[i] / 100)}`, x, y + 9, { color: INK, size: 10, align: 'center', halo: HALO });
    }
    if (probe >= 0) { const x = toX(model.mesh.lonCell[probe]), y = toY(model.mesh.latCell[probe]); ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(x, y, 7, 0, 2 * Math.PI); ctx.stroke(); }
    const length = KEY * pxPer / HOUR, kx = mapW - 14 - length, rows = sum ? [[`a, b + c: ${KEY} m/s per hour`, MUTED], [`all four: ${KEY / MAGNIFY} m/s per hour`, INK]] : [[`${KEY} m/s per hour`, INK]];
    ctx.font = '400 10px system-ui, sans-serif';
    const label = Math.max(...rows.map(([s]) => ctx.measureText(s).width));
    ctx.fillStyle = HALO; ctx.beginPath(); ctx.roundRect(kx - label - 16, 5, length + label + 24, 4 + 16 * rows.length, 4); ctx.fill();
    rows.forEach(([s, tone], n) => { arrow(ctx, kx, 15 + 16 * n, kx + length, 15 + 16 * n, { color: tone, width: 1.5, head: 5, dash: [4, 3], open: true }); text(ctx, s, kx - 6, 15 + 16 * n, { color: MUTED, size: 10, align: 'right' }); });
    ctx.restore();
  }

  function diagram(ctx, [x0, y0, x1, y1]) {
    const L = model.layers[layer], points = [[0, 0]];
    for (const t of TERMS) { const [e, n] = vec(t, probe), [px, py] = points[points.length - 1]; points.push([px + e * HOUR, py + n * HOUR]); }
    text(ctx, 'at the ring, tip to tail', x0 + 10, y0 + 14, { color: MUTED, size: 11 });
    const top = y0 + 30, bottom = y1 - 32, left = x0 + 16, right = x1 - 16, reach = 0.24 * Math.min(right - left, bottom - top);
    const we = L.wind[2 * probe], wn = L.wind[2 * probe + 1], speed = Math.hypot(we, wn), [ex, ey] = points[4], magnify = (gap) => (gap >= 14 ? 1 : MAGNIFY);
    let scale = 1, ox = 0, oy = 0, times = 1;
    for (let pass = 0; pass < 4; pass++) {
      const all = [...points, [ex * times, ey * times], ...(speed > 0.1 ? [[we / speed * reach / scale, wn / speed * reach / scale]] : [])], xs = all.map((p) => p[0]), ys = all.map((p) => p[1]);
      const xMin = Math.min(...xs), xMax = Math.max(...xs), yMin = Math.min(...ys), yMax = Math.max(...ys);
      scale = Math.min((right - left) / Math.max(xMax - xMin, 1e-3), (bottom - top) / Math.max(yMax - yMin, 1e-3), 40);
      ox = (left + right) / 2 - (xMax + xMin) / 2 * scale; oy = (top + bottom) / 2 + (yMax + yMin) / 2 * scale;
      times = magnify(Math.hypot(ex, ey) * scale);
    }
    const X = (p) => ox + p[0] * scale, Y = (p) => oy - p[1] * scale;
    if (speed > 0.1) arrow(ctx, ox, oy, ox + we / speed * reach, oy - wn / speed * reach, { color: MUTED, width: 1.5, head: 6 });
    const mx = points.reduce((sum, p) => sum + X(p), 0) / points.length, my = points.reduce((sum, p) => sum + Y(p), 0) / points.length;
    TERMS.forEach((t, n) => {
      const ax = X(points[n]), ay = Y(points[n]), bx = X(points[n + 1]), by = Y(points[n + 1]), len = Math.hypot(bx - ax, by - ay), cx = (ax + bx) / 2, cy = (ay + by) / 2;
      arrow(ctx, ax, ay, bx, by, { color: color[t], width: 2, head: 7, dash: [5, 4], open: true });
      const side = (cx - mx) * (by - ay) + (cy - my) * (ax - bx) < 0 ? -10 : 10;
      if (len > 14) text(ctx, t, cx + (by - ay) / len * side, cy + (ax - bx) / len * side, { color: color[t], size: 12, weight: 700, align: 'center', halo: HALO });
    });
    const inside = (m) => { const x = ox + ex * scale * m, y = oy - ey * scale * m; return x > x0 + 6 && x < x1 - 6 && y > y0 + 22 && y < y1 - 26; };
    while (times > 1 && !inside(times)) times = [MAGNIFY, 2, 1].find((m) => m < times);
    arrow(ctx, ox, oy, ox + ex * scale * times, oy - ey * scale * times, { color: INK, width: 2, head: 7, dash: [5, 4], open: true });
    ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(ox, oy, 3, 0, 2 * Math.PI); ctx.fill();
    const unit = nice(30 / scale), length = unit * scale;
    arrow(ctx, x0 + 12, y1 - 14, x0 + 12 + length, y1 - 14, { color: INK, width: 1.5, head: 5, dash: [4, 3], open: true });
    text(ctx, `${unit} m/s per hour`, x0 + 18 + length, y1 - 14, { color: MUTED, size: 10 });
    return { times, speed };
  }

  function notes(ctx, [x0, y0, edge], { times, speed }) {
    const x1 = Math.min(edge, x0 + 230);
    let y = y0 + (geo.narrow ? 16 : 6);
    text(ctx, geo.w < 320 ? 'in m/s per hour' : 'each term, in m/s per hour', x0 + 10, y, { color: MUTED, size: 10 });
    for (const t of [...TERMS, 'sum']) {
      y += 17;
      const [e, n] = vec(t, probe), sum = t === 'sum';
      if (!sum) text(ctx, t, x0 + 12, y, { color: color[t], size: 12, weight: 700 });
      text(ctx, sum ? 'all four' : NAMES[t], x0 + 26, y, { color: sum ? INK : MUTED, size: 11 });
      text(ctx, size(Math.hypot(e, n) * HOUR), x1 - 12, y, { color: color[t], size: 11, align: 'right', weight: sum ? 600 : 400 });
    }
    y += 22;
    text(ctx, `the wind, solid: ${speed.toFixed(0)} m/s`, x0 + 10, y, { color: MUTED, size: 10 });
    if (times > 1) text(ctx, `all four: drawn ${times}× longer`, x0 + 10, y + 15, { color: MUTED, size: 10 });
  }

  function report() {
    const L = model.layers[layer], { mesh, pi } = model, i = probe, mag = (t) => Math.hypot(...vec(t, i)) * HOUR;
    const [ae, an] = vec('a', i), [pe, pn] = vec('bc', i), angle = Math.acos(clamp((ae * pe + an * pn) / (Math.hypot(ae, an) * Math.hypot(pe, pn) || 1), -1, 1)) / DEG;
    out.set([
      ['at the ring', `${latName(mesh.latCell[i] / DEG)} ${lonName(mesh.lonCell[i] / DEG)}, ${Math.round(L.sigma * pi[i] / 100)} hPa, ${(L.geopotential[i] / model.header.g / 1000).toFixed(1)} km up`],
      ['wind', `${Math.hypot(L.wind[2 * i], L.wind[2 * i + 1]).toFixed(1)} m/s`],
      ['in m/s per hour: turning a', size(mag('a'))], ['height slope b', size(mag('b'))], ['pressure change c', size(mag('c'))], ['up and down d', size(mag('d'))], ['all four', size(mag('sum'))],
      ['the smoothing, not drawn', size(mag('closure'))],
      ['a and b + c', `${angle.toFixed(0)}° apart`],
      ['over the map, at the median cell', `the pressure gradient and the turning term cancel to within ${(L.median * 100).toPrecision(2)} percent`],
    ]);
  }

  function draw(ctx, w, h) {
    if (!model) { text(ctx, failed ? 'the model’s state did not load' : 'loading the model’s state…', w / 2, h / 2, { color: MUTED, size: 12, align: 'center' }); return; }
    if (geo?.key !== `${w}x${h}`) geo = layout(w, h);
    if (!geo.image || painted !== paletteVersion) paint();
    map(ctx);
    ctx.strokeStyle = GRID; ctx.lineWidth = 1; ctx.beginPath();
    if (geo.narrow) { ctx.moveTo(0, geo.mapH + 0.5); ctx.lineTo(w, geo.mapH + 0.5); } else { ctx.moveTo(geo.mapW + 0.5, 0); ctx.lineTo(geo.mapW + 0.5, h); }
    ctx.stroke();
    notes(ctx, geo.panel.notes, diagram(ctx, geo.panel.diagram));
    report();
  }

  function pick({ x, y }) {
    if (!geo || !model || x > geo.mapW || y > geo.mapH) return false;
    let best = -1, bestD = Infinity;
    for (const [i, cx, cy] of geo.cells) { const d = (cx - x) ** 2 + (cy - y) ** 2; if (d < bestD) { bestD = d; best = i; } }
    if (best >= 0 && best !== probe) { probe = best; fig.render(); }
    return best >= 0;
  }
  const at = (e) => { const r = fig.canvas.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
  const ring = ({ x, y }) => !!geo && probe >= 0 && Math.hypot(x - geo.toX(model.mesh.lonCell[probe]), y - geo.toY(model.mesh.latCell[probe])) < 24;
  fig.pointer({ hit: ring, down: pick, move: pick });
  fig.canvas.addEventListener('click', (e) => pick(at(e)));
  fig.canvas.addEventListener('pointermove', (e) => { if (e.pointerType === 'mouse' && !e.buttons) pick(at(e)); });

  loadState(new URL('../data/jw06state_N16.bin', import.meta.url).href).then((state) => { model = prepare(momentumBudget(state)); probe = model.probe; setKey(); fig.render(); }).catch(() => { failed = true; fig.render(); });
}
