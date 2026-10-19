import { Figure, choice, readout, legend, caption, text, rampRGB, clamp, paletteVersion, DARK_NEUTRAL, ACCENT, MUTED, GRID, LINE } from '../runtime.module.js';
import { loadFrames } from '../frames.module.js';
import { DATA_URL, layerSection, zonalMeans, modeWeight, modeName } from '../oceanRun.module.js';
import { overturningGrid, overturningCells, northwardByLayer, contours, SHALLOW, FLOOR } from '../seaZones.module.js';

const DARK = 'rgb(52, 55, 62)', HALO = 'rgba(21, 21, 21, 0.7)', PROBE = 'rgba(255, 255, 255, 0.45)';
const MONTHS = [['the year', ''], ['March', 'March'], ['June', 'June'], ['September', 'September'], ['December', 'December']];
const LEVELS = Array.from({ length: 25 }, (_, n) => 10 + 10 * n).flatMap((v) => [v, -v]);
const signedRoot = (v, top) => 0.5 + 0.5 * Math.sign(v) * Math.sqrt(Math.min(1, Math.abs(v) / top));
const VIEWS = [
  { key: 'T', label: 'temperature', color: (v, rgb) => rampRGB((v + 2) / 26, rgb, DARK_NEUTRAL), legend: ['ramp', 'temperature, from −2 to 24\u00a0°C', 'cool', 'warm', DARK] },
  { key: 'S', label: 'salinity', color: (v, rgb) => rampRGB((v - 33.6) / 2, rgb, DARK_NEUTRAL), legend: ['ramp', 'salinity, from 33.6 to 35.6\u00a0psu', 'cool', 'warm', DARK] },
  { key: 'u', label: 'eastward current', color: (v, rgb) => rampRGB(signedRoot(v, 100), rgb, DARK_NEUTRAL), legend: ['ramp', 'current toward the west in cool colors and the east in warm, up to 100\u00a0cm/s on a square-root scale', 'cool', 'warm', DARK] },
];
const POLAR = 40;
const interval = (step) => (step === 20 ? `20\u00a0Sv between lines, 10\u00a0Sv poleward of ${POLAR}°` : `${step}\u00a0Sv between lines, 20\u00a0Sv poleward of ${POLAR}°`);
const lines = (step) => [['line', 'the base of the mixed layer', ACCENT], ['faint', 'interfaces between layers, brighter where they follow isopycnals', 'rgba(255, 255, 255, 0.5)'], ['dash', 'the 1024.0 isopycnal', '#fff'], ['arrow', `the overturning, ${interval(step)}`, 'rgba(255, 255, 255, 0.9)']];

const place = (lat) => (Math.abs(lat) < 0.05 ? '0°' : `${Math.abs(lat).toFixed(Math.abs(lat) % 1 ? 1 : 0)}°${lat > 0 ? 'N' : 'S'}`);
const meters = (z) => `${Math.round(z).toLocaleString('en-US')}\u00a0m`;
const mix = (a, b, t) => (a === null ? b : b === null ? a : a + t * (b - a));
const toward = (v, plus, minus, digits = 1) => `${Math.abs(v).toFixed(digits)}\u00a0cm/s ${v >= 0 ? plus : minus}`;

export function mountSeazones(root) {
  const controls = root.querySelector('.controls');
  let header = null, failed = false, view = 0, month = '', probe = { lat: 2.5, z: 30 }, geometry = null, image = null, shown = null, step = 20;
  const periods = new Map();

  const fig = new Figure(root, { height: (w) => (w < 520 ? 390 : 440), draw });
  caption(root, 'The sea of the run above, averaged around each circle of latitude over the year or over one month, from the surface to the floor 4,000\u00a0m down; the top 500\u00a0m is drawn seven times as tall for each meter as the water below. Each band between the faint lines is one of the model’s 58 layers. An interface is faint where its layer is held at its minimum thickness, as the layers are down to 234 to 284\u00a0m in the tropics, and brighter where the layer follows a surface of equal density. The mixed layer is not one of the layers; its base and the 1024.0 isopycnal are the model’s own measures, averaged. The white lines are the overturning, the flow north or south of all the water above, the eddy transport included: water moves along them the way their arrows point, 20\u00a0Sv between neighboring lines, or on a narrow screen 40\u00a0Sv, and 10 or 20\u00a0Sv poleward of 40°. “The top 57\u00a0m carry” is that flow above the base of the twelve surface layers; the wind’s flow is not closed there, for most of it goes on to 300\u00a0m and turns back below. The six-hourly samples miss up to 0.5\u00a0Sv over the year and up to 3\u00a0Sv in a month; the figure spreads what they miss through the depth so that every line closes. The probe’s northward current is each layer’s flow across the circles 5° apart over its cross-section. Tap or hover over the section to read one place.');
  const key = legend(root, [VIEWS[view].legend, ...lines(step)]);
  const relabel = () => key.set([VIEWS[view].legend, ...lines(step)]);
  choice(controls, { label: 'Color shows', options: VIEWS.map((v, n) => [v.label, String(n)]), value: '0', onChange: (v) => { view = Number(v); relabel(); fig.render(); }, span: true });
  choice(controls, { label: 'Average over', options: MONTHS, value: month, onChange: (v) => { month = v; summarize(); fig.render(); }, span: true });
  const local = readout(controls), global = readout(controls);

  function period() {
    if (!periods.has(month)) {
      const name = month || null, grid = overturningGrid(header, name);
      periods.set(month, {
        section: layerSection(header, name), zonal: zonalMeans(header, name), grid, north: northwardByLayer(header, name), cells: overturningCells(header, name),
        lines: LEVELS.flatMap((level) => contours(grid.psi, grid.ni, grid.nj, level).map((line) => ({ ...line, level, polar: line.points.every(([i]) => Math.abs(grid.lats[0] + i * (grid.lats[1] - grid.lats[0])) >= POLAR) }))),
      });
    }
    return periods.get(month);
  }

  const L = () => header.layers;
  const binOf = (lat) => { const x = clamp((lat + 87.5) / 5, 0, 35), b = Math.min(34, Math.floor(x)); return [b, x - b]; };
  const capOf = (lat) => { const x = clamp((lat + 85) / 5, 0, 34), c = Math.min(33, Math.floor(x)); return [c, x - c]; };
  const basesAt = (p, lat, out) => { const [b, t] = binOf(lat), n = L(), base = p.section.base; for (let k = 0; k < n; k++) out[k] = (1 - t) * base[b * n + k] + t * base[(b + 1) * n + k]; return out; };
  const valueAt = (p, which, k, lat) => {
    const n = L();
    if (which === 'v') { const [c, t] = capOf(lat); return mix(p.north.values[c * n + k], p.north.values[(c + 1) * n + k], t); }
    const [b, t] = binOf(lat), values = p.section[which];
    return mix(values[b * n + k], values[(b + 1) * n + k], t);
  };
  const layerAt = (bases, z) => { let k = 0; while (k < bases.length - 1 && z > bases[k]) k++; return k; };

  function paint(box) {
    const width = Math.round(box.right - box.left), height = Math.round(box.bottom - box.top), id = `${view}/${month}/${width}x${height}/${paletteVersion}`;
    if (shown === id) return image;
    const p = period(), which = VIEWS[view].key, color = VIEWS[view].color, n = L(), bases = new Float64Array(n), rgb = [0, 0, 0];
    image ??= document.createElement('canvas');
    image.width = width; image.height = height;
    const ctx = image.getContext('2d'), pixels = ctx.createImageData(width, height);
    for (let x = 0; x < width; x++) {
      const lat = box.latitude(box.left + x + 0.5);
      basesAt(p, lat, bases);
      let k = 0, shownLayer = -1, value = null, last = null;
      for (let y = 0; y < height; y++) {
        const z = box.depth(box.top + y + 0.5);
        while (k < n - 1 && z > bases[k]) k++;
        if (k !== shownLayer) { shownLayer = k; value = valueAt(p, which, k, lat) ?? last; if (value !== null) { last = value; color(value, rgb); } }
        const i = 4 * (y * width + x);
        pixels.data[i] = rgb[0] * 255; pixels.data[i + 1] = rgb[1] * 255; pixels.data[i + 2] = rgb[2] * 255; pixels.data[i + 3] = 255;
      }
    }
    ctx.putImageData(pixels, 0, 0);
    shown = id;
    return image;
  }

  function draw(ctx, w, h) {
    if (!header) { text(ctx, failed ? 'The run could not be loaded.' : 'Loading the run…', w / 2, h / 2, { align: 'center', color: MUTED }); return; }
    const narrow = w < 520, size = narrow ? 10 : 11, left = narrow ? 40 : 48, right = w - (narrow ? 10 : 14), top = 22, bottom = h - 26;
    if (step !== (narrow ? 40 : 20)) { step = narrow ? 40 : 20; relabel(); }
    const half = (bottom - top) / 2;
    const box = {
      left, right, top, bottom,
      x: (lat) => left + (lat + 90) / 180 * (right - left),
      y: (z) => (z <= SHALLOW ? top + z / SHALLOW * half : top + half + (Math.min(z, FLOOR) - SHALLOW) / (FLOOR - SHALLOW) * half),
      latitude: (px) => (px - left) / (right - left) * 180 - 90,
      depth: (py) => (py - top <= half ? (py - top) / half * SHALLOW : SHALLOW + (py - top - half) / half * (FLOOR - SHALLOW)),
    };
    geometry = box;
    const p = period(), n = L(), section = p.section, centers = section.centers;
    ctx.drawImage(paint(box), left, top, right - left, bottom - top);

    const sides = [-90, ...centers, 90], bins = [0, ...centers.map((_, b) => b), centers.length - 1];
    const trace = (depthAt, keep, style) => {
      const z = bins.map(depthAt);
      for (let m = 0; m + 1 < sides.length; m++) {
        if (!keep(bins[m], bins[m + 1])) continue;
        if (style) ctx.strokeStyle = style(bins[m], bins[m + 1]);
        ctx.beginPath(); ctx.moveTo(box.x(sides[m]), box.y(z[m])); ctx.lineTo(box.x(sides[m + 1]), box.y(z[m + 1])); ctx.stroke();
      }
    };
    const thickness = (b, k) => section.base[b * n + k] - (k ? section.base[b * n + k - 1] : 0);
    const follows = (b, k) => (k < header.fixedLayers ? 0 : modeWeight(thickness(b, k), header.minimumThickness[k]));
    ctx.lineWidth = 1;
    for (let k = 1; k < n; k++) trace((b) => section.base[b * n + k - 1], (a, c) => thickness(a, k - 1) >= 2 || thickness(c, k - 1) >= 2, (a, c) => `rgba(255, 255, 255, ${(0.06 + 0.26 * 0.5 * (follows(a, k - 1) + follows(c, k - 1))).toFixed(3)})`);
    const { mixedLayerDepth: mixed, thermoclineTop: iso } = p.zonal;
    ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]);
    trace((b) => iso[b], (a, c) => iso[a] >= 1 && iso[c] >= 1);
    ctx.setLineDash([]);
    ctx.strokeStyle = ACCENT; ctx.lineWidth = 2;
    trace((b) => mixed[b], () => true);

    ctx.strokeStyle = GRID; ctx.lineWidth = 1; ctx.setLineDash([2, 3]);
    ctx.beginPath(); ctx.moveTo(left, box.y(SHALLOW)); ctx.lineTo(right, box.y(SHALLOW)); ctx.stroke(); ctx.setLineDash([]);

    const { lats, ni, nj, depth } = p.grid;
    const screen = ([i, j]) => {
      const i0 = Math.min(ni - 2, Math.floor(i)), s = i - i0, j0 = Math.min(nj - 2, Math.floor(j)), t = j - j0;
      const at = (c) => (1 - t) * depth[c * nj + j0] + t * depth[c * nj + j0 + 1];
      return [box.x(lats[i0] + s * (lats[i0 + 1] - lats[i0])), box.y((1 - s) * at(i0) + s * at(i0 + 1))];
    };
    const spacing = narrow ? 120 : 160;
    for (const line of p.lines) {
      if (line.level % step && !line.polar) continue;
      const pts = smooth(line.points.map(screen), line.closed);
      const lengths = [0];
      for (let m = 1; m < pts.length; m++) lengths.push(lengths[m - 1] + Math.hypot(pts[m][0] - pts[m - 1][0], pts[m][1] - pts[m - 1][1]));
      const total = lengths.at(-1);
      if (line.closed && total < 16) continue;
      for (const [color, width] of [[HALO, 3], ['rgba(255, 255, 255, 0.85)', 1.3]]) {
        ctx.strokeStyle = color; ctx.lineWidth = width; ctx.lineJoin = 'round'; ctx.beginPath();
        pts.forEach(([px, py], m) => (m ? ctx.lineTo(px, py) : ctx.moveTo(px, py)));
        if (line.closed) ctx.closePath();
        ctx.stroke();
      }
      if (total < 50) continue;
      for (let at = total < spacing ? total / 2 : spacing / 2; at < total; at += spacing) {
        let m = 1;
        while (m < pts.length - 1 && lengths[m] < at) m++;
        const [ax, ay] = pts[m - 1], [bx, by] = pts[m], len = Math.hypot(bx - ax, by - ay) || 1, u = (at - lengths[m - 1]) / len;
        head(ctx, ax + u * (bx - ax), ay + u * (by - ay), (bx - ax) / len, (by - ay) / len);
      }
    }

    ctx.lineWidth = 1;
    for (const z of [0, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000]) {
      ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(left - 4, box.y(z)); ctx.lineTo(left, box.y(z)); ctx.stroke();
      text(ctx, z.toLocaleString('en-US'), left - 6, box.y(z), { align: 'right', color: MUTED, size });
    }
    text(ctx, 'm', left - 6, top - 12, { align: 'right', color: MUTED, size });
    for (const lat of [-60, -30, 0, 30, 60]) {
      ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(box.x(lat), bottom); ctx.lineTo(box.x(lat), bottom + 4); ctx.stroke();
      text(ctx, place(lat), box.x(lat), bottom + 13, { align: 'center', color: MUTED, size });
    }
    const labelAt = centers.indexOf(17.5);
    if (p.zonal.thermoclineTop[labelAt] >= 1) text(ctx, '1024.0', box.x(17.5) + 4, box.y(p.zonal.thermoclineTop[labelAt]) - 8, { color: '#fff', size, halo: HALO });

    const px = box.x(probe.lat), py = box.y(probe.z);
    ctx.strokeStyle = PROBE; ctx.lineWidth = 1; ctx.setLineDash([2, 3]);
    ctx.beginPath(); ctx.moveTo(px, top); ctx.lineTo(px, bottom); ctx.moveTo(left, py); ctx.lineTo(right, py); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = HALO; ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(px, py, 5, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
    report(p);
  }

  function report(p) {
    const n = L(), bases = basesAt(p, probe.lat, new Float64Array(n)), k = layerAt(bases, probe.z), label = header.labels[k], thick = bases[k] - (k ? bases[k - 1] : 0), weight = k < header.fixedLayers ? 0 : modeWeight(thick, header.minimumThickness[k]);
    const at = (which) => valueAt(p, which, k, probe.lat), T = at('T'), S = at('S'), u = at('u'), v = at('v');
    const { lats, nj, psi, depth } = p.grid, i0 = Math.min(lats.length - 2, Math.floor((probe.lat + 90) / 5)), s = (probe.lat - lats[i0]) / 5;
    const column = (c) => {
      let j = 0;
      while (j < nj - 2 && depth[c * nj + j + 1] < probe.z) j++;
      const d0 = depth[c * nj + j], d1 = depth[c * nj + j + 1], t = d1 > d0 ? clamp((probe.z - d0) / (d1 - d0), 0, 1) : 1;
      return psi[c * nj + j] + t * (psi[c * nj + j + 1] - psi[c * nj + j]);
    };
    const above = (1 - s) * column(i0) + s * column(i0 + 1);
    local.set([
      ['latitude', place(probe.lat)], ['depth', meters(probe.z)],
      ['in', `layer ${k + 1} of ${L()}, ${meters(k ? bases[k - 1] : 0)} to ${meters(bases[k])}`],
      ['its label', `${Number.isInteger(label) ? label.toFixed(1) : label}\u00a0kg/m³, ${modeName(weight)}`],
      ['temperature', T === null ? '–' : `${T.toFixed(1)}\u00a0°C`], ['salinity', S === null ? '–' : `${S.toFixed(2)}\u00a0psu`],
      ['current', `${u === null ? '–' : toward(u, 'east', 'west')}, ${v === null ? '–' : toward(v, 'north', 'south', 2)}`],
      ['overturning above', `${Math.abs(above).toFixed(1)}\u00a0Sv ${above >= 0 ? 'north' : 'south'}`],
    ]);
  }

  function summarize() {
    const p = period(), { zonal, cells } = p, deepest = zonal.mixedLayerDepth.reduce((best, z, b) => (z > best.z ? { z, lat: zonal.centers[b] } : best), { z: 0 });
    const across = (psi, lat) => `${Math.abs(psi).toFixed(1)}\u00a0Sv ${psi >= 0 ? 'north' : 'south'} across ${place(lat)}`;
    global.set([
      ['deepest mixed layer', `${meters(deepest.z)} at ${place(deepest.lat)}`],
      [`the top ${Math.round(cells.base)}\u00a0m carry`, `${across(cells.south, -5)}, ${across(cells.north, 5)}`],
      ['and under the westerlies', `${across(cells.westerlySouth, -45)}, ${across(cells.westerlyNorth, 45)}`],
    ]);
  }

  const pick = (event) => {
    if (!geometry) return;
    const r = fig.canvas.getBoundingClientRect(), x = clamp(event.clientX - r.left, geometry.left, geometry.right), y = clamp(event.clientY - r.top, geometry.top, geometry.bottom);
    const z = geometry.depth(y);
    probe = { lat: clamp(Math.round(geometry.latitude(x) / 2.5) * 2.5, -87.5, 87.5), z: z < SHALLOW ? Math.round(z) : 10 * Math.round(z / 10) };
    fig.render();
  };
  fig.canvas.style.touchAction = 'pan-y';
  fig.canvas.addEventListener('pointerdown', pick);
  fig.canvas.addEventListener('pointermove', (event) => { if (event.pointerType === 'mouse' || event.buttons) pick(event); });

  loadFrames(DATA_URL).then((frames) => { header = frames.header; summarize(); fig.render(); }).catch((error) => { failed = true; console.info('seazones figure', error.message); fig.render(); });
}

function head(ctx, x, y, ux, uy) {
  ctx.save();
  ctx.translate(x, y);
  ctx.beginPath(); ctx.moveTo(ux * 4, uy * 4); ctx.lineTo(-ux * 3.5 - uy * 3.5, -uy * 3.5 + ux * 3.5); ctx.lineTo(-ux * 3.5 + uy * 3.5, -uy * 3.5 - ux * 3.5); ctx.closePath();
  ctx.fillStyle = '#fff'; ctx.strokeStyle = HALO; ctx.lineWidth = 1; ctx.fill(); ctx.stroke();
  ctx.restore();
}

function smooth(points, closed) {
  let pts = points;
  for (let pass = 0; pass < 2; pass++) {
    const out = closed ? [] : [pts[0]], n = pts.length, last = closed ? n : n - 1;
    for (let m = 0; m < last; m++) {
      const [ax, ay] = pts[m], [bx, by] = pts[(m + 1) % n];
      out.push([0.75 * ax + 0.25 * bx, 0.75 * ay + 0.25 * by], [0.25 * ax + 0.75 * bx, 0.25 * ay + 0.75 * by]);
    }
    if (!closed) out.push(pts[n - 1]);
    pts = out;
  }
  return pts;
}
