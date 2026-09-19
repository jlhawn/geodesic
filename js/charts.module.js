import { Grid } from './grid.module.js';
import { createModel } from './model.module.js';
import { cellVector } from './dynamics/operators.module.js';
import { LEVELS, levelFields } from './levels.module.js';

const A1 = 1.340264, A2 = -0.081106, A3 = 0.000893, A4 = 0.003796;
const WIDTH = 1600, MARGIN = 20, TOP = 70, BOTTOM = 90;

function project(lat, lon) {
  const theta = Math.asin(Math.sqrt(3) / 2 * Math.sin(lat));
  const t2 = theta * theta, t6 = t2 * t2 * t2;
  const denom = 3 * (9 * A4 * t6 * t2 + 7 * A3 * t6 + 3 * A2 * t2 + A1);
  return [2 * Math.sqrt(3) * lon * Math.cos(theta) / denom, A4 * t6 * t2 * theta + A3 * t6 * theta + A2 * t2 * theta + A1 * theta];
}
const [xMax] = project(0, Math.PI);
const [, yMax] = project(Math.PI / 2, 0);
const scale = (WIDTH - 2 * MARGIN) / (2 * xMax);
const HEIGHT = Math.round(TOP + 2 * yMax * scale + BOTTOM);
const px = ([x, y]) => [((x + xMax) * scale + MARGIN).toFixed(1), (TOP + (yMax - y) * scale).toFixed(1)];
const lonOf = (x, y) => Math.atan2(y, x);
const latOf = (x, y, z) => Math.atan2(z, Math.hypot(x, y));
const unwrap = (lon, reference) => { while (lon - reference > Math.PI) lon -= 2 * Math.PI; while (lon - reference < -Math.PI) lon += 2 * Math.PI; return lon; };

function outlinePath() {
  const pts = [];
  for (let s = -90; s <= 90; s += 2) pts.push(project(s * Math.PI / 180, Math.PI));
  for (let s = 90; s >= -90; s -= 2) pts.push(project(s * Math.PI / 180, -Math.PI));
  return 'M' + pts.map((p) => px(p).join(' ')).join('L') + 'Z';
}

function colorFor(t) {
  const x = Math.max(0, Math.min(1, t));
  const lo = [232, 240, 251], hi = [8, 48, 107];
  return `rgb(${lo.map((c, i) => Math.round(c + (hi[i] - c) * Math.pow(x, 0.85))).join(',')})`;
}

function niceStep(range, target = 14) {
  const raw = range / target;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  for (const m of [1, 2, 2.5, 4, 5, 10]) if (m * mag >= raw) return m * mag;
  return 10 * mag;
}

const models = new Map();
function modelFor(N) {
  if (!models.has(N)) models.set(N, createModel(new Grid(N)));
  return models.get(N);
}

function contours(mesh, field, step) {
  const { nVertices, cellsOnVertex, xCell } = mesh;
  let fmin = Infinity, fmax = -Infinity;
  for (const v of field) { fmin = Math.min(fmin, v); fmax = Math.max(fmax, v); }
  const polylines = [];
  for (let level = Math.ceil(fmin / step) * step; level < fmax; level += step) {
    const segments = [];
    const point = (a, b) => {
      const t = (level - field[a]) / (field[b] - field[a]);
      const x = xCell[3 * a] + t * (xCell[3 * b] - xCell[3 * a]), y = xCell[3 * a + 1] + t * (xCell[3 * b + 1] - xCell[3 * a + 1]), z = xCell[3 * a + 2] + t * (xCell[3 * b + 2] - xCell[3 * a + 2]);
      const n = Math.hypot(x, y, z);
      return { key: a < b ? `${a},${b}` : `${b},${a}`, lat: latOf(x / n, y / n, z / n), lon: lonOf(x / n, y / n) };
    };
    for (let v = 0; v < nVertices; v++) {
      const cells = [cellsOnVertex[3 * v], cellsOnVertex[3 * v + 1], cellsOnVertex[3 * v + 2]];
      const crossings = [];
      for (let m = 0; m < 3; m++) {
        const a = cells[m], b = cells[(m + 1) % 3];
        if ((field[a] - level) * (field[b] - level) < 0) crossings.push(point(a, b));
      }
      if (crossings.length === 2) segments.push(crossings);
    }
    const byKey = new Map();
    segments.forEach((s, idx) => { for (const p of s) { if (!byKey.has(p.key)) byKey.set(p.key, []); byKey.get(p.key).push(idx); } });
    const used = new Uint8Array(segments.length);
    for (let start = 0; start < segments.length; start++) {
      if (used[start]) continue;
      const chain = [segments[start][0], segments[start][1]];
      used[start] = 1;
      for (const direction of [1, -1]) {
        for (;;) {
          const end = direction === 1 ? chain[chain.length - 1] : chain[0];
          const next = (byKey.get(end.key) || []).find((idx) => !used[idx]);
          if (next === undefined) break;
          used[next] = 1;
          const seg = segments[next];
          const other = seg[0].key === end.key ? seg[1] : seg[0];
          if (direction === 1) chain.push(other); else chain.unshift(other);
        }
      }
      polylines.push({ level, points: chain });
    }
  }
  return { polylines, fmin, fmax };
}

function svgFor(model, grid, title, subtitle, fill, contourField, contourStep, unit) {
  const mesh = model.mesh;
  const C = mesh.nCells;
  let smin = Infinity, smax = -Infinity;
  for (const v of fill) { smin = Math.min(smin, v); smax = Math.max(smax, v); }
  const { polylines, fmin, fmax } = contours(mesh, contourField, contourStep);
  const cellLon = Float64Array.from({ length: C }, (_, i) => lonOf(mesh.xCell[3 * i], mesh.xCell[3 * i + 1]));
  const parts = [`<svg xmlns="http://www.w3.org/2000/svg" width="${WIDTH}" height="${HEIGHT}" viewBox="0 0 ${WIDTH} ${HEIGHT}" font-family="Helvetica, Arial, sans-serif">`,
    `<rect width="${WIDTH}" height="${HEIGHT}" fill="#ffffff"/>`, `<defs><clipPath id="map"><path d="${outlinePath()}"/></clipPath></defs>`,
    `<text x="${MARGIN}" y="30" font-size="22" fill="#111">${title}</text>`, `<text x="${MARGIN}" y="52" font-size="14" fill="#555">${subtitle}</text>`, `<g clip-path="url(#map)">`];
  for (const cell of grid) {
    const i = cell.index;
    const center = cellLon[i];
    const verts = cell.vertices.map((v) => [latOf(v.x, v.y, v.z), unwrap(lonOf(v.x, v.y), center)]);
    const straddles = verts.some(([, lon]) => Math.abs(lon) > Math.PI);
    const color = colorFor((fill[i] - smin) / (smax - smin || 1));
    for (const shift of straddles ? [0, center > 0 ? -2 * Math.PI : 2 * Math.PI] : [0]) {
      parts.push(`<polygon points="${verts.map(([lat, lon]) => px(project(lat, lon + shift)).join(',')).join(' ')}" fill="${color}" stroke="${color}" stroke-width="0.5"/>`);
    }
  }
  for (let lat = -60; lat <= 60; lat += 30) {
    const pts = []; for (let lon = -180; lon <= 180; lon += 2) pts.push(px(project(lat * Math.PI / 180, lon * Math.PI / 180)).join(' '));
    parts.push(`<path d="M${pts.join('L')}" fill="none" stroke="#000" stroke-opacity="0.15" stroke-width="0.8"/>`);
  }
  for (let lon = -150; lon <= 180; lon += 30) {
    const pts = []; for (let lat = -90; lat <= 90; lat += 2) pts.push(px(project(lat * Math.PI / 180, lon * Math.PI / 180)).join(' '));
    parts.push(`<path d="M${pts.join('L')}" fill="none" stroke="#000" stroke-opacity="0.15" stroke-width="0.8"/>`);
  }
  const labels = [];
  for (const { level, points } of polylines) {
    if (points.length < 2) continue;
    const unwrapped = [];
    let previous = points[0].lon;
    for (const p of points) { const lon = unwrap(p.lon, previous); unwrapped.push([p.lat, lon]); previous = lon; }
    for (const shift of [0, 2 * Math.PI, -2 * Math.PI]) {
      const pts = unwrapped.map(([lat, lon]) => px(project(lat, lon + shift)));
      parts.push(`<path d="M${pts.map((p) => p.join(' ')).join('L')}" fill="none" stroke="#1a1a1a" stroke-width="1.2" stroke-linejoin="round"/>`);
      let travelled = 0, nextLabel = 220 + (level * 7919) % 200;
      for (let n = 1; n < pts.length; n++) {
        const dx = pts[n][0] - pts[n - 1][0], dy = pts[n][1] - pts[n - 1][1];
        const length = Math.hypot(dx, dy);
        travelled += length;
        if (travelled >= nextLabel && length > 4) {
          const mx = (Number(pts[n][0]) + Number(pts[n - 1][0])) / 2, my = (Number(pts[n][1]) + Number(pts[n - 1][1])) / 2;
          const lonMid = 0.5 * (unwrapped[n][1] + unwrapped[n - 1][1]) + shift;
          if (Math.abs(lonMid) < 0.97 * Math.PI && my > TOP && my < TOP + 2 * yMax * scale) {
            let angle = Math.atan2(dy, dx) * 180 / Math.PI;
            if (angle > 90) angle -= 180;
            if (angle < -90) angle += 180;
            labels.push(`<text x="${mx.toFixed(1)}" y="${my.toFixed(1)}" transform="rotate(${angle.toFixed(1)} ${mx.toFixed(1)} ${my.toFixed(1)})" font-size="12" fill="#111" text-anchor="middle" dominant-baseline="middle" stroke="#fff" stroke-width="3" paint-order="stroke">${Number.isInteger(level) ? level : level.toFixed(1)}</text>`);
          }
          nextLabel += 420;
        }
      }
    }
  }
  parts.push(`</g>`, `<path d="${outlinePath()}" fill="none" stroke="#333" stroke-width="1"/>`, ...labels);
  const barX = MARGIN, barY = HEIGHT - 52, barW = 420, barH = 14;
  parts.push(`<defs><linearGradient id="bar" x1="0" x2="1" y1="0" y2="0">${[0, 0.25, 0.5, 0.75, 1].map((t) => `<stop offset="${t}" stop-color="${colorFor(t)}"/>`).join('')}</linearGradient></defs>`,
    `<rect x="${barX}" y="${barY}" width="${barW}" height="${barH}" fill="url(#bar)" stroke="#555" stroke-width="0.5"/>`);
  for (const t of [0, 0.25, 0.5, 0.75, 1]) parts.push(`<text x="${(barX + t * barW).toFixed(1)}" y="${barY + barH + 16}" font-size="12" fill="#333" text-anchor="middle">${(smin + t * (smax - smin)).toFixed(t === 0 || t === 1 ? 1 : 0)}</text>`);
  parts.push(`<text x="${barX}" y="${barY - 6}" font-size="13" fill="#333">Wind speed (m/s)</text>`,
    `<text x="${barX + barW + 40}" y="${barY + 11}" font-size="13" fill="#333">Contours every ${contourStep} ${unit}, from ${fmin.toFixed(0)} to ${fmax.toFixed(0)} ${unit}</text>`, `</svg>`);
  return parts.join('\n');
}

/*
 * Builds the chart for one level of a saved state: isobars over the
 * lowest-layer wind at the surface, or geopotential-height contours over
 * the wind interpolated to that pressure.
 */
export function chartFor(state, level) {
  const model = modelFor(state.N);
  const grid = [...new Grid(state.N)];
  const { K, C, E, sigmaMid } = model.core.diagnostics;
  const pi = Float64Array.from(state.pi), theta = Float64Array.from(state.theta), u = Float64Array.from(state.u);
  model.core.diagnose(pi, theta);
  const layerWind = [];
  for (let k = 0; k < K; k++) layerWind.push(cellVector(model.mesh, u.subarray(k * E, (k + 1) * E), new Float64Array(3 * C)));
  if (level === 'surface') {
    const psHpa = Float64Array.from(pi, (p) => p / 100);
    const speed = Float64Array.from({ length: C }, (_, i) => Math.hypot(layerWind[K - 1][3 * i], layerWind[K - 1][3 * i + 1], layerWind[K - 1][3 * i + 2]));
    return svgFor(model, grid, `Day ${state.day} — surface pressure and surface wind`, `Isobars of surface pressure (hPa); fill: wind speed in the lowest layer (${(sigmaMid[K - 1] * 1000).toFixed(0)} hPa at 1000 hPa surface pressure). N=${state.N}, ${(7720 / state.N).toFixed(0)} km cells.`, speed, psHpa, 4, 'hPa');
  }
  const p = Number(level);
  const { speed, height } = levelFields(model.core, pi, theta, (k) => layerWind[k], p);
  const step = niceStep(Math.max(...height) - Math.min(...height));
  return svgFor(model, grid, `Day ${state.day} — ${p} hPa geopotential height and wind speed`, `Contours: height of the ${p} hPa surface (m); fill: wind speed (m/s) interpolated to ${p} hPa in ln p${p >= 990 ? '; extrapolated hydrostatically where the level is below ground' : ''}. N=${state.N}.`, speed, height, step, 'm');
}

async function listStates(directory) {
  const response = await fetch(directory);
  const html = await response.text();
  const names = [...html.matchAll(/href="([^"]+_state_day\d+\.json)"/g)].map((m) => decodeURIComponent(m[1]));
  return names.sort();
}

export default async function runCharts(directory = 'runs/') {
  const stateSelect = document.getElementById('state');
  const levelSelect = document.getElementById('level');
  const status = document.getElementById('status');
  const holder = document.getElementById('chart');
  const svgButton = document.getElementById('download-svg');
  const pngButton = document.getElementById('download-png');
  for (const level of LEVELS) {
    const option = document.createElement('option');
    option.value = String(level);
    option.textContent = level === 'surface' ? 'Surface' : `${level} hPa`;
    levelSelect.appendChild(option);
  }
  let names = [];
  try { names = await listStates(directory); } catch (error) { status.textContent = `could not list ${directory}: ${error.message}`; return; }
  if (!names.length) { status.textContent = `no *_state_day*.json files in ${directory}`; return; }
  for (const name of names) {
    const option = document.createElement('option');
    option.value = name;
    option.textContent = name.replace('_state_day', ' · day ').replace('.json', '');
    stateSelect.appendChild(option);
  }
  const params = new URLSearchParams(location.search);
  if (params.get('state') && names.includes(params.get('state'))) stateSelect.value = params.get('state');
  if (params.get('level')) levelSelect.value = params.get('level');
  const states = new Map();
  let current = '';
  async function render() {
    const name = stateSelect.value, level = levelSelect.value;
    status.textContent = `loading ${name}…`;
    if (!states.has(name)) states.set(name, fetch(directory + name).then((r) => r.json()));
    const state = await states.get(name);
    status.textContent = `rendering ${level === 'surface' ? 'surface' : level + ' hPa'}…`;
    await new Promise((resolve) => setTimeout(resolve, 0));
    const t0 = performance.now();
    current = chartFor(state, level);
    holder.innerHTML = current;
    status.textContent = `${name} · ${level === 'surface' ? 'surface' : level + ' hPa'} · ${(performance.now() - t0).toFixed(0)} ms`;
    history.replaceState(null, '', `?state=${encodeURIComponent(name)}&level=${level}`);
  }
  const download = (blob, filename) => {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  };
  const filename = (ext) => `${stateSelect.value.replace('.json', '')}_${levelSelect.value === 'surface' ? 'surface' : levelSelect.value + 'hPa'}.${ext}`;
  svgButton.addEventListener('click', () => { if (current) download(new Blob([current], { type: 'image/svg+xml' }), filename('svg')); });
  pngButton.addEventListener('click', () => {
    if (!current) return;
    const image = new Image();
    image.onload = () => {
      const canvas = document.createElement('canvas');
      canvas.width = WIDTH * 2;
      canvas.height = HEIGHT * 2;
      const context = canvas.getContext('2d');
      context.drawImage(image, 0, 0, canvas.width, canvas.height);
      canvas.toBlob((blob) => download(blob, filename('png')), 'image/png');
    };
    image.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(current);
  });
  stateSelect.addEventListener('change', render);
  levelSelect.addEventListener('change', render);
  render();
}
