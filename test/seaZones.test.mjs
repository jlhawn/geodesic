import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { readFrames } from '../explain/frames.module.js';
import { layerOverturning, layerSection } from '../explain/oceanRun.module.js';
import { contours, overturningGrid, overturningCells, northwardByLayer, THICK, CELL_DEPTH } from '../explain/seaZones.module.js';

const bytes = readFileSync(new URL('../explain/data/ocean_N16.bin', import.meta.url));
const { header } = readFrames('ocean_N16.bin', bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength));
const L = header.layers, M = L + 1, F = header.fixedLayers;
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const grid = (ni, nj, f) => Float64Array.from({ length: ni * nj }, (_, n) => f(Math.floor(n / nj), n % nj));
const onEdge = (values, nj, [x, y]) => {
  const at = (a, b) => values[a * nj + b];
  if (Number.isInteger(x)) { const j = Math.floor(y), t = y - j; return (1 - t) * at(x, j) + t * at(x, Math.min(nj - 1, j + 1)); }
  const i = Math.floor(x), s = x - i;
  return (1 - s) * at(i, y) + s * at(i + 1, y);
};
const highOnRight = (values, nj, line, level) => {
  const above = ([i, j]) => values[i * nj + j] > level;
  for (let n = 0; n + 1 < line.points.length; n++) {
    const [a, b] = [line.points[n], line.points[n + 1]], d = [b[0] - a[0], b[1] - a[1]];
    const side = ([i, j]) => (i - a[0]) * -d[1] + (j - a[1]) * d[0];
    for (const [x, y] of [a, b]) {
      const ends = Number.isInteger(x) ? [[x, Math.floor(y)], [x, Math.floor(y) + 1]] : [[Math.floor(x), y], [Math.floor(x) + 1, y]];
      for (const corner of ends) if (side(corner) !== 0 && (side(corner) > 0) !== above(corner)) return false;
    }
  }
  return true;
};

test('a contour around a peak is closed, lies on the level and runs with the larger values on its right', () => {
  const ni = 21, nj = 17, peak = grid(ni, nj, (i, j) => Math.exp(-((i - 10) ** 2 + (j - 8) ** 2) / 20));
  const lines = contours(peak, ni, nj, 0.5);
  assert.equal(lines.length, 1);
  const [line] = lines;
  assert.equal(line.closed, true);
  for (const p of line.points) close(onEdge(peak, nj, p), 0.5, 1e-12, `value at ${p}`);
  assert.ok(highOnRight(peak, nj, line, 0.5));
  const top = line.points.reduce((best, p) => (p[1] < best[1] ? p : best)), next = line.points[(line.points.indexOf(top) + 1) % line.points.length];
  assert.ok(next[0] > top[0], 'over a peak drawn with j down the page, the line runs clockwise');
  const trough = grid(ni, nj, (i, j) => -peak[i * nj + j]), [around] = contours(trough, ni, nj, -0.5);
  assert.ok(highOnRight(trough, nj, around, -0.5));
  const low = around.points.reduce((best, p) => (p[1] < best[1] ? p : best)), after = around.points[(around.points.indexOf(low) + 1) % around.points.length];
  assert.ok(after[0] < low[0], 'around a trough it runs counterclockwise');
});

test('an open contour runs edge to edge, and a value growing down the page makes it run to the right', () => {
  const ni = 6, nj = 5, ramp = grid(ni, nj, (i, j) => j);
  const lines = contours(ramp, ni, nj, 2.5);
  assert.equal(lines.length, 1);
  assert.equal(lines[0].closed, false);
  assert.equal(lines[0].points.length, ni);
  lines[0].points.forEach((p, n) => { close(p[0], n, 1e-12, 'i'); close(p[1], 2.5, 1e-12, 'j'); });
  const across = contours(grid(ni, nj, (i) => i), ni, nj, 2.5);
  assert.ok(across[0].points.at(-1)[1] < across[0].points[0][1], 'a value growing to the right makes the line run up the page');
});

test('a saddle splits into two lines, each with the larger values on its right', () => {
  const ni = 4, nj = 4, saddle = grid(ni, nj, (i, j) => (i - 1.5) * (j - 1.5) + 0.01 * i);
  for (const level of [0.3, -0.3]) {
    const lines = contours(saddle, ni, nj, level);
    assert.equal(lines.length, 2, `level ${level}`);
    for (const line of lines) assert.ok(highOnRight(saddle, nj, line, level), `level ${level}`);
  }
});

test('the overturning grid is zero at the poles and the surface, carries only the free surface’s own flow at the floor, and keeps the layer streamfunction between them', () => {
  for (const month of [null, 'March', 'June', 'September', 'December']) {
    const g = overturningGrid(header, month), flow = layerOverturning(header, month);
    assert.equal(g.ni, 37);
    assert.equal(g.nj, M);
    assert.deepEqual([g.lats[0], g.lats[1], g.lats[18], g.lats[36]], [-90, -85, 0, 90]);
    for (let j = 0; j < M; j++) { assert.equal(g.psi[j], 0); assert.equal(g.psi[36 * M + j], 0); }
    for (let i = 0; i < g.ni; i++) { assert.equal(g.psi[i * M], 0); close(g.psi[i * M + L], i && i < 36 ? flow.etaNet[i - 1] : 0, 0.05, `floor at ${g.lats[i]}`); }
    if (!month) for (let i = 0; i < g.ni; i++) assert.ok(Math.abs(g.psi[i * M + L]) < 0.4);
    for (let c = 0; c < 35; c++) for (let k = 0; k <= L; k++) { assert.equal(g.psi[(c + 1) * M + k], flow.psi[c * M + k]); assert.equal(g.depth[(c + 1) * M + k], flow.depth[c * M + k]); }
    const base = layerSection(header, month).base;
    for (let k = 1; k <= L; k++) { assert.equal(g.depth[k], base[k - 1]); assert.equal(g.depth[36 * M + k], base[35 * L + k - 1]); }
    for (let i = 0; i < g.ni; i++) for (let k = 0; k < L; k++) assert.ok(g.depth[i * M + k + 1] >= g.depth[i * M + k] - 1e-9);
  }
});

test('the overturning cells the figure and the page read: the tropical cells, the cells under the westerlies, the deep cells, and where they close', () => {
  const year = overturningCells(header);
  close(year.base, 56.9, 0.05, 'the base of the twelve surface layers');
  close(year.south, -107.3, 0.05, 'south across 5°S'); close(year.north, 102.5, 0.05, 'north across 5°N');
  close(year.westerlySouth, 31.4, 0.05, 'north across 45°S'); close(year.westerlyNorth, -26.9, 0.05, 'south across 45°N');
  close(year.surface(-10), -37.3, 0.05, 'south across 10°S'); close(year.surface(10), 32.1, 0.05, 'north across 10°N'); close(year.surface(0), 13.7, 0.05, 'north across the equator');
  const tropicS = year.peak(-1, -30, 0), tropicN = year.peak(1, 0, 30), westS = year.peak(1, -60, -30), westN = year.peak(-1, 30, 60);
  assert.deepEqual([tropicS.latitude, Math.round(tropicS.depth)], [-5, 57]); assert.deepEqual([tropicN.latitude, Math.round(tropicN.depth)], [5, 57]);
  assert.deepEqual([westS.latitude, Math.round(westS.depth)], [-45, 351]); assert.deepEqual([westN.latitude, Math.round(westN.depth)], [40, 289]);
  close(westS.psi, 40.0, 0.05, 'the southern westerly cell'); close(westN.psi, -36.6, 0.05, 'the northern'); assert.ok(westS.depth < CELL_DEPTH);
  const deepS = year.deep(-1, -90, 0), deepN = year.deep(1, 0, 90);
  assert.deepEqual([deepS.latitude, Math.round(deepS.depth)], [-10, 1107]); assert.deepEqual([deepN.latitude, Math.round(deepN.depth)], [5, 1064]);
  close(deepS.psi, -33.8, 0.05, 'the southern deep cell'); close(deepN.psi, 37.6, 0.05, 'the northern deep cell');
  for (const [z, south, north] of [[57, -107.2, 102.5], [100, -96.2, 92.9], [300, -96.2, 91.5], [500, -6.9, -7.1], [1000, -26.0, 39.0], [3000, -6.7, 5.8]]) { close(year.at(-5, z), south, 0.1, `5°S at ${z} m`); close(year.at(5, z), north, 0.1, `5°N at ${z} m`); }
  close(year.at(-45, 100), 34.8, 0.05, '45°S at 100 m'); close(year.at(-45, 500), 14.7, 0.05, 'at 500 m'); close(year.at(-45, 1000), 1.6, 0.05, 'at 1,000 m');
  const months = Object.fromEntries(['March', 'June', 'September', 'December'].map((m) => [m, overturningCells(header, m)]));
  close(months.September.south, -237.3, 0.05, 'September, 5°S'); close(months.September.north, 46.0, 0.05, 'September, 5°N');
  close(months.March.south, -31.5, 0.05, 'March, 5°S'); close(months.March.north, 230.7, 0.05, 'March, 5°N');
  close(months.June.south, -130.4, 0.05, 'June, 5°S'); close(months.December.north, 139.5, 0.05, 'December, 5°N');
  close(months.September.westerlySouth, 30.8, 0.05, 'September, 45°S'); close(months.September.westerlyNorth, -33.8, 0.05, 'September, 45°N');
  close(months.September.at(0, 300), -90.1, 0.05, 'September across the equator above 300 m'); close(months.March.at(0, 100), 110.5, 0.05, 'March above 100 m');
  const g = overturningGrid(header), inner = contours(g.psi, g.ni, g.nj, -100);
  assert.equal(inner.length, 1);
  assert.ok(inner[0].closed);
  for (const [i] of inner[0].points) assert.ok(g.lats[0] + 5 * i > -7.5 && g.lats[0] + 5 * i < -2.5);
  const outer = contours(g.psi, g.ni, g.nj, 100);
  assert.equal(outer.length, 1); assert.ok(outer[0].closed);
  for (const [i] of outer[0].points) assert.ok(g.lats[0] + 5 * i > 2.5 && g.lats[0] + 5 * i < 7.5);
});

test('the northward current by layer is each layer’s transport over its cross-section, and matches the section where both exist', () => {
  const R = header.constants.radius, v = northwardByLayer(header), base = header.section.base;
  assert.equal(v.values.length, 35 * L);
  const c = header.overturning.caps.indexOf(-10), h = 0.5 * (base[c * L] + base[(c + 1) * L]);
  close(v.values[c * L], 100 * header.overturning.mass[c * L] * 1e6 / (2 * Math.PI * R * Math.cos(10 * Math.PI / 180) * h), 1e-12, 'top layer at 10°S');
  for (const lat of [-35, -20, 20, 35]) for (const k of [0, 5, 11]) {
    const cap = header.overturning.caps.indexOf(lat), sectionMean = 0.5 * (header.section.v[cap * L + k] + header.section.v[(cap + 1) * L + k]);
    close(v.values[cap * L + k], sectionMean, 0.3, `layer ${k} at ${lat}°`);
  }
  close(v.values[c * L], -2.40, 0.005, 'top layer at 10°S, cm/s'); close(v.values[header.overturning.caps.indexOf(10) * L], 2.14, 0.005, 'at 10°N');
  close(v.values[header.overturning.caps.indexOf(-5) * L], -8.35, 0.005, 'at 5°S'); close(v.values[header.overturning.caps.indexOf(5) * L], 8.10, 0.005, 'at 5°N');
  for (let n = 0; n < v.values.length; n++) {
    const cap = Math.floor(n / L), k = n % L, thick = (b) => base[b * L + k] - (k ? base[b * L + k - 1] : 0);
    assert.equal(v.values[n] === null, 0.5 * (thick(cap) + thick(cap + 1)) < THICK);
  }
  const september = northwardByLayer(header, 'September'), flow = header.overturning.months.find((m) => m.name === 'September'), sb = header.section.months.find((m) => m.name === 'September').base;
  close(september.values[c * L], 100 * flow.mass[c * L] * 1e6 / (2 * Math.PI * R * Math.cos(10 * Math.PI / 180) * 0.5 * (sb[c * L] + sb[(c + 1) * L])), 1e-12, 'September, top layer at 10°S');
});

class Node {
  constructor(tag) { Object.assign(this, { tag, children: [], style: {}, dataset: {}, attributes: {}, listeners: {}, className: '', textContent: '', clientWidth: 0, width: 0, height: 0 }); }
  append(...c) { this.children.push(...c); }
  before() {} after() {}
  replaceChildren(...c) { this.children = c; }
  setAttribute(k, v) { this.attributes[k] = v; }
  addEventListener(type, fn) { (this.listeners[type] ??= []).push(fn); }
  querySelector(sel) { return this.children.find((c) => c.className === sel.slice(1)); }
  getBoundingClientRect() { return { left: 0, top: 0, width: this.width, height: this.height }; }
  getContext() { return (this.ctx ??= context(this)); }
}
const text = (node) => (typeof node === 'string' ? node : node.children.length ? node.children.map(text).join('') : node.textContent ?? '').replaceAll('\u00a0', ' ');
function context(canvas) {
  const target = { measureText: (s) => ({ width: 6 * String(s).length }), createImageData: (w, h) => ({ width: w, height: h, data: new Uint8ClampedArray(w * h * 4) }), createLinearGradient: () => ({ addColorStop() {} }) };
  return new Proxy(target, {
    get(t, prop) {
      if (prop in t) return t[prop];
      return (...args) => { if (args.some((a) => typeof a === 'number' && !Number.isFinite(a))) canvas.bad.push(prop); if (prop === 'fillText') canvas.texts.push(String(args[0])); };
    },
    set(t, prop, value) { t[prop] = value; return true; },
  });
}
async function mount(name, mountName, width) {
  globalThis.IntersectionObserver = class { observe() {} }; globalThis.ResizeObserver = class { observe() {} };
  globalThis.matchMedia = () => ({ addEventListener() {} }); globalThis.window = { devicePixelRatio: 1 };
  globalThis.getComputedStyle = () => ({ getPropertyValue: () => '' });
  const canvases = [];
  globalThis.document = { documentElement: {}, createTextNode: (t) => t, createElement: (tag) => { const n = new Node(tag); if (tag === 'canvas') { n.bad = []; n.texts = []; canvases.push(n); } return n; } };
  globalThis.fetch = async (url) => ({ arrayBuffer: async () => bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength) });
  const root = new Node('figure'), stage = new Node('div'), controls = new Node('div');
  stage.className = 'stage'; controls.className = 'controls'; stage.clientWidth = width; root.append(stage, controls);
  const module = await import(`../explain/figures/${name}.module.js`);
  module[mountName](root);
  await new Promise((resolve) => setTimeout(resolve, 600));
  const canvas = canvases[0];
  const buttons = controls.children.flatMap((c) => c.children), readouts = () => controls.children.filter((c) => c.className === 'readout').map(text);
  return {
    canvas, readouts,
    click: (label) => { buttons.find((b) => b.textContent === label).listeners.click[0](); },
    move: (fx, fy) => canvas.listeners.pointermove.at(-1)({ pointerType: 'mouse', buttons: 0, clientX: canvas.width * fx, clientY: canvas.height * fy }),
  };
}

test('the section figure reads what the page quotes: the surface layers’ cells, the probe’s layer and mode, and draws nothing non-finite at either width', async () => {
  for (const width of [640, 360]) {
    const fig = await mount('seazones', 'mountSeazones', width);
    assert.deepEqual(fig.canvas.bad, []);
    const [local, global] = fig.readouts();
    assert.match(local, /latitude 2\.5°N · depth 30 m · in layer 8 of 58, 27 m to 32 m · its label 1018\.25 kg\/m³, fixed-depth · temperature 22\.1 °C · salinity 34\.91 psu · current 77\.0 cm\/s west, 2\.52 cm\/s north · overturning above 47\.3 Sv north/);
    assert.equal(global, 'deepest mixed layer 157 m at 87.5°N · the top 57 m carry 107.3 Sv south across 5°S, 102.5 Sv north across 5°N · and under the westerlies 31.4 Sv north across 45°S, 26.9 Sv south across 45°N');
    fig.click('September'); assert.match(fig.readouts()[1], /the top 57 m carry 237\.3 Sv south across 5°S, 46\.0 Sv north across 5°N · and under the westerlies 30\.8 Sv north across 45°S, 33\.8 Sv south across 45°N/);
    fig.click('March'); assert.match(fig.readouts()[1], /deepest mixed layer 380 m at 87\.5°N · the top 57 m carry 31\.5 Sv south across 5°S, 230\.7 Sv north across 5°N/);
    fig.click('the year');
    for (const view of ['salinity', 'eastward current', 'temperature']) fig.click(view);
    fig.move(0.55, 0.3);
    assert.match(fig.readouts()[0], width === 640 ? /latitude 5°N · depth 281 m · in layer 34 of 58, 264 m to 288 m · its label 1025\.0 kg\/m³, isopycnal · temperature 13\.7 °C · salinity 34\.99 psu · current 45\.8 cm\/s west/ : /in layer \d+ of 58, .* its label 10\d\d\.\d+ kg\/m³, (isopycnal|fixed-depth|between)/);
    fig.move(0.5, 0.05);
    assert.match(fig.readouts()[0], /in layer 1 of 58, 0 m to 3 m · its label 1016\.5 kg\/m³, fixed-depth/);
    assert.deepEqual(fig.canvas.bad, []);
    assert.ok(fig.canvas.texts.includes('1024.0'), 'the 1024.0 isopycnal is labeled at 17.5°N');
  }
});

test('the heat figure reads what the page quotes: the peaks, the sea’s store, the cells’ heat and the surface layers’ flow', async () => {
  const fig = await mount('seaheat', 'mountSeaheat', 640);
  assert.deepEqual(fig.canvas.bad, []);
  const [upper, lower, global] = fig.readouts();
  assert.match(upper, /latitude 10°S/);
  assert.match(lower, /heat carried 3\.07 PW south · by the air, as the heat of dry air 0\.97 PW south · as latent heat 0\.11 PW north · by the sea 2\.21 PW south · by its cells 2\.83 PW south/);
  assert.match(global, /most heat carried, in all 4\.9 PW south at 30°S, 4\.8 PW north at 30°N · by the air 3\.9 PW south at 35°S, 3\.9 PW north at 35°N · by the sea 2\.3 PW south at 15°S, 2\.4 PW north at 15°N · across the equator none: by the sea none, by the air none · heat the sea takes in, and stores 8\.5 W\/m², 8\.5 W\/m² · heat the run’s books cannot place 0\.003 W\/m², spread evenly/);
  fig.click('the sea’s overturning');
  assert.match(fig.readouts()[1], /water carried by the top 57 m 37\.3 Sv south · lighter than 1026, below them 14\.5 Sv north · 1026 to 1026\.8 22\.2 Sv north · from 1026\.8 none/);
  fig.click('the surface’s budget');
  assert.match(fig.readouts()[0], /sunlight absorbed 176\.4 W\/m²/);
  fig.move(0.5, 0.5);
  assert.deepEqual(fig.canvas.bad, []);
});
