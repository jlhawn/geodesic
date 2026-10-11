// Shadows the drawing of infrared.module.js and beam.module.js; a fix to either belongs here too.
import { Figure, choice, legend, caption, readout, text, termColor, clamp, INK, MUTED, GRID, LINE, ACCENT } from '../runtime.module.js';
import { loadColumn, K, GLOBAL_BEAM, GLOBAL_MU } from '../radiativeColumn.module.js';
import { slices, ABSORBERS } from '../infraredSlices.module.js';
import { createSplit } from '../sunlightSplit.module.js';
import { SOLAR_CONSTANT } from '../../js/physics/radiation.module.js';

const WIDE = 520, LABELED = 0.02, P_TOP = 10, T_LOW = 180, T_HIGH = 300, RATE = 7, HALO = 'rgba(20, 20, 22, 0.85)', GROUND = '#5a4a36';
const COLORS = { carbonDioxide: '#d55e00', ozone: '#0072b2', methane: '#009e73', nitrousOxide: '#cc79a7', vaporLines: '#5fc4f0', continuum: '#bfe6fa', none: '#5c5c5c' };
const NAMES = { carbonDioxide: 'carbon dioxide', ozone: 'ozone', methane: 'methane', nitrousOxide: 'nitrous oxide', vaporLines: 'vapor’s lines', continuum: 'vapor’s continuum' };
const BUTTONS = { carbonDioxide: 'carbon dioxide', ozone: 'ozone', methane: 'methane', nitrousOxide: 'nitrous oxide', vaporLines: 'vapor lines', continuum: 'continuum' };
const FORMULAS = { carbonDioxide: 'CO₂', ozone: 'O₃', methane: 'CH₄', nitrousOxide: 'N₂O', vaporLines: 'H₂O', continuum: 'continuum' };
const SHORT = { carbonDioxide: 'CO₂', ozone: 'O₃', methane: 'CH₄', nitrousOxide: 'N₂O', vaporLines: 'lines', continuum: 'continuum' };
const VAPORS = [['none', '0'], ['×¼', '0.25'], ['×½', '0.5'], ['×1', '1'], ['×2', '2']];
const AMOUNTS = [['none', '0'], ['×½', '0.5'], ['×1', '1'], ['×2', '2'], ['×4', '4']];
const SUNS = { overhead: { beam: SOLAR_CONSTANT, elevation: 90, ground: 'sea', rate: 30, tick: 10 }, mean: { beam: GLOBAL_BEAM, elevation: Math.asin(GLOBAL_MU) * 180 / Math.PI, ground: 0.07, rate: 10, tick: 2 } };
const FILL = { continuum: 0.8 };
const RGB = Object.fromEntries(Object.entries(COLORS).map(([k, hex]) => [k, [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16))]));
const ZOOM = 3, TOP = 30, PAD = 12, BARS = 184, GAP = 2, SKY = '#dfe6ee', SHINE = '#9aa3ad', SOIL = '#80664a', EARTH = '#5a4a36';
const PRESSURES = [1000, 300, 100, 30, 10, 3, 1, 0.3], HEIGHTS = [10, 20, 30, 40, 50];
const LIGHT = ['vapor', 'ozone', 'gases', 'aerosol', 'byAir', 'byGround', 'ground'];

const rgb = (c, a = 1) => `rgba(${Math.round(c[0])}, ${Math.round(c[1])}, ${Math.round(c[2])}, ${a})`;
const mix = (hex, to, t) => `rgb(${[1, 3, 5].map((n) => Math.round(parseInt(hex.slice(n, n + 2), 16) * (1 - t) + to * t)).join(', ')})`;
const watts = (v) => `${(Math.abs(v) < 0.05 ? 0 : v).toFixed(1)} W/m²`;
const percent = (f) => `${(100 * f).toFixed(1)}%`;
const number = (v, digits = 1) => v.toLocaleString('en-US', { minimumFractionDigits: digits, maximumFractionDigits: digits });
const SUPER = { 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
const pretty = (s) => s.replace(/\^(\d+)/, (_, n) => [...n].map((c) => SUPER[c]).join(''));
const depthText = (d) => (d >= 100 ? d.toExponential(1).replace('e+', ' × 10^') : d >= 10 ? d.toFixed(0) : d >= 1 ? d.toFixed(1) : d.toFixed(d >= 0.01 ? 2 : 4));
const part = (d) => (d >= 100 ? pretty(depthText(d)) : String(Number(d.toPrecision(2))));
const pressureText = (p) => (p >= 10000 ? `${(p / 100).toFixed(0)} hPa` : p >= 1000 ? `${(p / 100).toFixed(1)} hPa` : `${(p / 100).toFixed(2)} hPa`);
const hPa = (p) => { const v = p / 100; return v >= 10 ? v.toFixed(0) : v >= 1 ? v.toFixed(1) : v.toFixed(2); };

function segments(ctx, x, y, height, per, parts) {
  for (const [value, color] of parts) {
    const width = Math.max(0, value) * per;
    if (value > 0) {
      const inset = width > 2 * GAP + 2 ? GAP / 2 : 0, w = Math.max(1, width - 2 * inset);
      ctx.fillStyle = color;
      ctx.beginPath(); ctx.roundRect(x + inset, y, w, height, Math.min(3, w / 2)); ctx.fill();
    }
    x += width;
  }
}

function labels(ctx, items, lo, hi, size, from, to, first) {
  ctx.font = `400 ${size}px system-ui, sans-serif`;
  for (const item of items) item.width = Math.max(ctx.measureText(item.name).width, item.value ? ctx.measureText(item.value).width : 0);
  let edge = lo;
  for (const item of items) { item.left = Math.max(edge, item.anchor - item.width / 2); edge = item.left + item.width + 10; }
  edge = hi;
  for (const item of [...items].reverse()) { item.left = Math.min(item.left, edge - item.width); edge = item.left - 10; }
  const line = size + 2.5;
  for (const item of items) {
    const center = Math.max(lo, item.left) + item.width / 2;
    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(item.anchor, from); ctx.lineTo(center, to); ctx.stroke();
    const rows = item.value ? [[item.name, MUTED], [item.value, INK]] : [[item.name, MUTED]], start = to < from ? first + (2 - rows.length) * line : first;
    rows.forEach(([words, color], n) => text(ctx, words, center, start + n * line, { align: 'center', color, size }));
  }
}

export function mountVaporslices(root) {
  const controls = root.querySelector('.controls'), A = termColor('a'), VAPOR = A, OZONE = mix(A, 255, 0.35), GAS = mix(A, 255, 0.6), HAZE = mix(A, 0, 0.35);
  const on = Object.fromEntries(ABSORBERS.map((g) => [g, true]));
  let profile = null, split = null, view = 'slices', humidity = 1, amount = 1, lastAmount = 1, sun = 'overhead';
  let state = null, shown = null, probe = null, touched = false, geometry = null, first = true;
  let light = null, lit = null, dry = {}, layer = null, plot = null;

  const fig = new Figure(root, { height: (w) => (w < WIDE ? 480 : 460), step, draw });
  const words = caption(root, '');
  const key = legend(root, []);

  choice(controls, { label: 'Show', options: [['infrared slices', 'slices'], ['sunlight', 'sunlight']], value: view, onChange: (v) => { view = v; arrange(); update(); }, span: true });
  choice(controls, { label: 'Water vapor', options: VAPORS, value: '1', onChange: (v) => { humidity = Number(v); arrange(); update(); }, span: true });
  const group = document.createElement('div'), label = document.createElement('span');
  group.className = 'choice'; group.setAttribute('role', 'group'); group.style.gridColumn = '1 / -1';
  label.className = 'caption'; label.textContent = 'Absorbers'; label.id = 'vaporslices-absorbers'; group.setAttribute('aria-labelledby', label.id);
  group.append(label);
  const toggles = Object.fromEntries(ABSORBERS.map((g) => {
    const b = document.createElement('button');
    b.type = 'button'; b.textContent = BUTTONS[g];
    b.addEventListener('click', () => {
      if (g === 'carbonDioxide') { amount = on[g] ? 0 : lastAmount; on[g] = amount > 0; dose.value = String(amount); } else on[g] = !on[g];
      update();
    });
    group.append(b);
    return [g, b];
  }));
  controls.append(group);
  const dose = choice(controls, { label: 'Carbon dioxide', options: AMOUNTS, value: '1', onChange: (v) => { amount = Number(v); if (amount > 0) lastAmount = amount; on.carbonDioxide = amount > 0; update(); }, span: true });
  const doseBox = controls.lastElementChild;
  choice(controls, { label: 'Sun', options: [['overhead, 1,362 W/m²', 'overhead'], ['the global mean, 340.5 W/m²', 'mean']], value: sun, onChange: (v) => { sun = v; update(); }, span: true });
  const sunBox = controls.lastElementChild;
  const totals = readout(controls), detail = readout(controls), detailBox = controls.lastElementChild;

  function arrange() {
    const slicing = view === 'slices';
    group.style.display = slicing ? '' : 'none';
    doseBox.style.display = slicing ? '' : 'none';
    sunBox.style.display = slicing ? 'none' : '';
    detailBox.style.display = slicing ? '' : 'none';
    if (slicing) {
      words.textContent = 'The midlatitude summer over ground at 294.2 K, with its water vapor scaled by the chosen factor. Along the bottom, the ground’s infrared is split into the model’s 34 slices, each as wide as its share, the clearest on the left. Each bar rises to where its slice leaves for space, the level at which the absorbers above it reach an optical depth of 1; where there is no bar, the ground shines through. The row along the top is what each slice sends to space, and the row under the slices is what the air sends back down to the ground, each as a share of the ground’s own emission in that slice. On the right is the air’s temperature. ×1 carbon dioxide is 388.75 ppm. Drag across the slices to read one.';
      key.set([...ABSORBERS.map((g) => ['bar', g === 'vaporLines' ? 'water vapor’s lines, H₂O' : g === 'continuum' ? 'water vapor’s continuum' : `${NAMES[g]}, ${FORMULAS[g]}`, COLORS[g]]), ['bar', 'nearly clear: under 5% absorbed', COLORS.none], ['line', 'temperature of the air', INK]]);
    } else {
      words.textContent = 'The same column in sunlight, through the model’s own radiation code: either a sun straight overhead, over the sea, or the global mean of 340.5 W/m² with the sun 41.8° up, over ground that reflects 7 percent. The top bar splits the arriving sunlight by where it ends up, and the air’s share is drawn again below it, three times larger; the sliver for O₂ and CO₂ also holds what the gases take from light on its way back up. Underneath is how fast the sunlight warms each layer, in K per day, with dry air’s heating dashed; the topmost layer, above 0.3 hPa, is left off. Point at or tap the profile to read another layer.';
      key.set([['bar', 'absorbed by water vapor', VAPOR], ['bar', 'absorbed by ozone', OZONE], ['bar', 'by O₂ and CO₂', GAS], ['bar', 'by aerosol', HAZE], ['bar', 'reflected by the air', SKY], ['bar', 'reflected by the ground', SHINE], ['bar', 'absorbed by the ground', SOIL], ['line', 'heating by sunlight', A], ...(humidity > 0 ? [['dash', 'heating in dry air', MUTED]] : [])]);
    }
  }

  function update() {
    for (const g of ABSORBERS) toggles[g].setAttribute('aria-pressed', String(on[g]));
    if (!profile) return;
    if (view === 'slices') {
      state = slices(profile, { humidity, carbonDioxide: amount, gases: { ...on, carbonDioxide: amount > 0 } });
      if (!shown) shown = state.points.map(() => ({ left: 0, level: 0, kept: 1, back: 0, alpha: 0.55, color: [0, 0, 0] }));
      if (probe === null) probe = state.points.reduce((best, p) => (p.share > state.points[best].share ? p.index : best), 0);
      report();
    } else {
      const s = SUNS[sun];
      dry[sun] ??= Float64Array.from(split({ ...s }).heating);
      const g = split({ ...s, humidity });
      light = { ...g, parts: Float64Array.from(LIGHT, (name) => g[name]), heating: Float64Array.from(g.heating), dry: dry[sun], rate: s.rate, tick: s.tick };
      if (!lit) lit = { parts: Float64Array.from(light.parts), heating: Float64Array.from(light.heating) };
      shine();
    }
    fig.play(true);
    fig.render();
  }

  function report() {
    totals.set([['vapor in the column', humidity > 0 ? `${state.precipitableWater.toFixed(1)} kg/m²` : 'none'], ['the ground emits', watts(state.groundEmission)], ['leaving the top', watts(state.olr)], ['held back, G', watts(state.greenhouse)], ['back down to the ground', watts(state.dlr)], ['straight to space', `${percent(state.escaping)} of the ground’s emission`]]);
    const p = state.points[probe], strongest = ABSORBERS.filter((g) => p.gas[g] >= 0.02 * p.depth && p.gas[g] > 0).sort((a, b) => p.gas[b] - p.gas[a]).slice(0, 4);
    const depth = `${pretty(depthText(p.depth))}${strongest.length > 1 ? `: ${strongest.map((g) => `${SHORT[g]} ${part(p.gas[g])}`).join(', ')}` : ''}`;
    const where = p.level === null ? [['straight to space', `${percent(Math.exp(-p.depth))} of it`]] : [['leaves for space from', `${pressureText(p.level)}, ${(p.height / 1000).toFixed(1)} km, where the air is at ${p.temperature.toFixed(1)} K`]];
    detail.set([[`slice ${p.index + 1} of ${state.points.length}`, `${p.dominant ? `mostly ${NAMES[p.dominant]}` : 'nearly clear'}, ${percent(p.share)} of the ground’s emission`], ['optical depth', depth], ...where, ['to space', `${watts(p.up)} of the ${watts(p.ground)} the ground emits in it`], ['back down', watts(p.down)]]);
  }

  function shine() {
    const g = light;
    totals.set([
      ['arriving', `${number(g.arriving, g.arriving > 1000 ? 0 : 1)} W/m²`],
      ['absorbed in the air', `${number(g.air)} W/m², ${percent(g.air / g.arriving)}`],
      ...(humidity > 0 ? [['by the vapor', `${number(g.vapor)} W/m²`]] : []),
      ['absorbed by the ground', `${number(g.ground)} W/m², ${percent(g.ground / g.arriving)}`],
      ['reflected to space', `${number(g.reflected)} W/m², ${percent(g.reflected / g.arriving)}`],
      ['in the lowest 10 km', `${g.lower.low.toFixed(2)} to ${g.lower.high.toFixed(2)} K/day`],
      ['vapor in the column', humidity > 0 ? `${g.precipitableWater.toFixed(1)} kg/m²` : 'none'],
    ]);
  }

  function targets() {
    const ps = state.surfacePressure, out = [];
    let left = 0;
    for (const j of state.order) {
      const p = state.points[j];
      out[j] = { left, level: Math.log((p.level ?? ps) / P_TOP), kept: p.ground > 0 ? p.up / p.ground : 1, back: p.ground > 0 ? p.down / p.ground : 0, alpha: FILL[p.dominant] ?? 0.55, color: RGB[p.dominant ?? 'none'] };
      left += p.share;
    }
    return out;
  }

  function step(dt) {
    if (view === 'slices') {
      if (!state) return false;
      const goal = targets(), f = first ? 1 : 1 - Math.exp(-RATE * dt);
      let moving = false;
      shown.forEach((s, j) => {
        const g = goal[j];
        s.left += (g.left - s.left) * f; s.level += (g.level - s.level) * f; s.kept += (g.kept - s.kept) * f; s.back += (g.back - s.back) * f; s.alpha += (g.alpha - s.alpha) * f;
        for (let c = 0; c < 3; c++) s.color[c] += (g.color[c] - s.color[c]) * f;
        if (Math.abs(g.left - s.left) > 1e-4 || Math.abs(g.level - s.level) > 1e-3 || Math.abs(g.kept - s.kept) > 1e-4 || Math.abs(g.back - s.back) > 1e-4 || Math.abs(g.color[0] - s.color[0]) + Math.abs(g.color[1] - s.color[1]) + Math.abs(g.color[2] - s.color[2]) > 0.5) moving = true;
        else { s.left = g.left; s.level = g.level; s.kept = g.kept; s.back = g.back; s.alpha = g.alpha; s.color = [...g.color]; }
      });
      first = false;
      if (!moving) fig.play(false);
      return true;
    }
    if (!light) return false;
    const f = 1 - Math.exp(-dt * 12);
    let moved = 0;
    for (const name of ['parts', 'heating']) for (let n = 0; n < light[name].length; n++) {
      moved = Math.max(moved, Math.abs(light[name][n] - lit[name][n]));
      lit[name][n] += (light[name][n] - lit[name][n]) * f;
    }
    if (moved < 1e-3) { lit.parts.set(light.parts); lit.heating.set(light.heating); fig.play(false); }
    return true;
  }

  function drawSlices(ctx, w, h) {
    if (first) step(0);
    const narrow = w < WIDE, ps = state.surfacePressure, left = narrow ? 46 : 52, right = w - (narrow ? 8 : 14), panel = narrow ? 66 : 120, gap = narrow ? 12 : 24;
    const x0 = left, x1 = right - panel - gap, t0 = x1 + gap, t1 = right, top = 50, bottom = h - 80, row = 10, rowTop = 26, rowBack = bottom + row + 6;
    const span = Math.log(ps / P_TOP), y = (lp) => top + lp / span * (bottom - top), yp = (p) => y(Math.log(p / P_TOP));
    const x = (f) => x0 + f * (x1 - x0), tx = (T) => t0 + (T - T_LOW) / (T_HIGH - T_LOW) * (t1 - t0);
    geometry = { x0, x1 };

    ctx.lineWidth = 1;
    for (const level of [0.1, 1, 10, 100, 1000]) {
      const yy = yp(level * 100);
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(x0 - 4, yy); ctx.lineTo(x1, yy); ctx.moveTo(t0, yy); ctx.lineTo(t1, yy); ctx.stroke();
      text(ctx, `${level}${level === 0.1 ? ' hPa' : ''}`, x0 - 7, yy, { align: 'right', color: MUTED, size: 10 });
    }
    for (const T of [200, 250, 300]) {
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(tx(T), top); ctx.lineTo(tx(T), bottom); ctx.stroke();
      if (!narrow || T !== 250) text(ctx, `${T}`, tx(T), bottom + 22, { align: T === 300 ? 'right' : 'center', color: MUTED, size: 10 });
    }
    text(ctx, narrow ? 'air, K' : 'air temperature, K', (t0 + t1) / 2, bottom + 38, { align: 'center', color: MUTED, size: 10 });
    text(ctx, `to space, ${watts(state.olr)}`, x0, rowTop - 10, { color: INK, size: 11 });
    text(ctx, `held back, ${watts(state.greenhouse)}`, t1, rowTop - 10, { align: 'right', color: MUTED, size: 11 });
    text(ctx, `back down to the ground, ${watts(state.dlr)}`, x0, rowBack + row + 12, { color: INK, size: 11 });
    for (const f of [0, 0.25, 0.5, 0.75, 1]) text(ctx, f === 1 ? '100%' : `${100 * f}`, x(f), bottom + 52, { align: f === 1 ? 'right' : f === 0 ? 'left' : 'center', color: MUTED, size: 10 });
    text(ctx, `the ground’s ${watts(state.groundEmission)}, slice by slice`, (x0 + x1) / 2, bottom + 67, { align: 'center', color: MUTED, size: 10 });

    const order = state.order.slice().sort((a, b) => shown[a].left - shown[b].left);
    for (const j of order) {
      const s = shown[j], p = state.points[j], a = x(s.left), b = x(s.left + p.share), width = b - a;
      if (width < 0.05) continue;
      const inner = width >= 3 ? width - 1 : width, color = s.color, yy = Math.min(bottom, y(s.level));
      ctx.fillStyle = rgb(color, s.alpha); ctx.fillRect(a, yy, inner, bottom - yy);
      if (bottom - yy > 0.5) { ctx.fillStyle = rgb(color); ctx.fillRect(a, yy, inner, Math.min(2, bottom - yy)); }
      ctx.fillStyle = rgb(color); ctx.fillRect(a, bottom + 2, inner, row);
      const kept = clamp(s.kept, 0, 1) * inner, back = clamp(s.back, 0, 1) * inner;
      ctx.fillRect(a + (inner - kept) / 2, rowTop, kept, row);
      ctx.fillRect(a + (inner - back) / 2, rowBack, back, row);
    }
    ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(x0, bottom + 0.5); ctx.lineTo(x1, bottom + 0.5); ctx.stroke();

    let clear = 0;
    for (const p of state.points) if (p.level === null) clear += p.share;
    if (clear > 0.3) {
      const words = narrow ? 'the ground shines through' : 'mostly clear: the ground shines through';
      ctx.font = '400 11px system-ui, sans-serif';
      const half = ctx.measureText(words).width / 2;
      text(ctx, words, clamp(x(clear / 2), x0 + half, x1 - half), yp(20000), { align: 'center', color: MUTED, size: 11 });
    }
    const groups = [];
    order.forEach((j, n) => {
      const p = state.points[j], last = groups.at(-1);
      if (!p.dominant) return;
      if (last && last.gas === p.dominant && last.next === n) { last.share += p.share; last.next = n + 1; if (p.share > last.widest.share) last.widest = p; } else groups.push({ gas: p.dominant, share: p.share, widest: p, next: n + 1 });
    });
    const placed = [];
    ctx.font = `500 11px system-ui, sans-serif`;
    if (probe !== null) {
      const s = shown[probe], a = x(s.left), b = Math.max(a + 2, x(s.left + state.points[probe].share)), yy = Math.min(bottom, y(s.level));
      for (let at = yy; at < bottom; at += 13) placed.push([a, at, 0], [b, at, 0]);
      if (!touched) placed.push([clamp((a + b) / 2, x0 + 14, x1 - 14), Math.min(yy, bottom - 4) - 24, ctx.measureText('drag').width / 2]);
    }
    for (const { gas, widest } of groups.filter((group) => group.share >= (narrow ? 2 : 1) * LABELED).sort((a, b) => b.share - a.share)) {
      const s = shown[widest.index], half = ctx.measureText(FORMULAS[gas]).width / 2, cx = clamp(x(s.left + widest.share / 2), x0 + half + 2, x1 - half - 2);
      let yy = Math.min(bottom, y(s.level)) - 9, tries = 0;
      const hits = () => placed.some(([lx, ly, lw]) => Math.abs(lx - cx) < half + lw + 6 && Math.abs(ly - yy) < 13);
      while (hits() && tries < 2) { yy -= 13; tries++; }
      if (hits()) continue;
      placed.push([cx, yy, half]);
      text(ctx, FORMULAS[gas], cx, yy, { align: 'center', color: INK, size: 11, weight: 500, halo: HALO });
    }

    ctx.strokeStyle = INK; ctx.lineWidth = 2; ctx.lineJoin = 'round'; ctx.beginPath();
    profile.temperature.forEach((T, k) => { const px = tx(T), py = yp(profile.pressure[k]); if (k) ctx.lineTo(px, py); else ctx.moveTo(px, py); });
    ctx.stroke();
    ctx.fillStyle = GROUND; ctx.fillRect(t0, bottom + 2, t1 - t0, 5);
    ctx.fillStyle = INK; ctx.beginPath(); ctx.arc(tx(state.groundTemperature), bottom, 3.5, 0, Math.PI * 2); ctx.fill();

    if (probe !== null) {
      const p = state.points[probe], s = shown[probe], a = x(s.left), b = Math.max(a + 2, x(s.left + p.share)), yy = Math.min(bottom, y(s.level));
      ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5; ctx.strokeRect(a - 0.5, Math.min(yy, bottom - 4) - 0.5, b - a, bottom + 2 + row - Math.min(yy, bottom - 4) + 1);
      ctx.strokeRect(a - 0.5, rowTop - 0.5, b - a, row + 1);
      ctx.strokeRect(a - 0.5, rowBack - 0.5, b - a, row + 1);
      const T = p.level === null ? state.groundTemperature : p.temperature, px = tx(T);
      ctx.strokeStyle = 'rgba(255,255,255,0.5)'; ctx.lineWidth = 1; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(b + 2, yy); ctx.lineTo(px, yy); ctx.stroke(); ctx.setLineDash([]);
      ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(px, yy, 3.5, 0, Math.PI * 2); ctx.fill();
      const tag = `${T.toFixed(0)} K`, low = yy > bottom - 14;
      text(ctx, tag, px + (T > 270 ? -7 : 7), low ? yy - 11 : yy - 10, { align: T > 270 ? 'right' : 'left', color: '#fff', size: 11, weight: 500, halo: HALO });
      if (!touched) text(ctx, 'drag', clamp((a + b) / 2, x0 + 14, x1 - 14), Math.min(yy, bottom - 4) - 24, { align: 'center', color: ACCENT, size: 11, halo: HALO });
    }
  }

  function drawBars(ctx, w, size) {
    const x0 = PAD, span = w - 2 * PAD, per = span / light.arriving, zoom = per * ZOOM, yBar = 66, yZoom = 118, bh = 18;
    const [vapor, ozone, gases, aerosol, byAir, byGround, ground] = lit.parts, air = vapor + ozone + gases + aerosol, mine = air + byAir;
    text(ctx, `sunlight arriving at the top: ${number(light.arriving, light.arriving > 1000 ? 0 : 1)} W/m²`, x0, 12, { color: INK, size: size + 0.5, weight: 500 });
    ctx.fillStyle = 'rgba(255, 255, 255, 0.05)'; ctx.strokeStyle = GRID; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.roundRect(x0 + 0.5, yBar + 0.5, span - 1, bh - 1, 3); ctx.fill(); ctx.stroke();
    segments(ctx, x0, yBar, bh, per, [[air, A], [byAir, SKY], [byGround, SHINE], [ground, SOIL]]);
    const xAir = x0 + mine * per, xShine = xAir + byGround * per, xEnd = xShine + ground * per;
    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(x0, yBar - 3); ctx.lineTo(x0, yBar - 6); ctx.lineTo(xAir, yBar - 6); ctx.lineTo(xAir, yBar - 3); ctx.stroke();
    labels(ctx, [
      { anchor: (x0 + xAir) / 2, name: 'the air' },
      { anchor: (xAir + xShine) / 2, name: 'reflected by the ground', value: number(light.byGround) },
      { anchor: (xShine + xEnd) / 2, name: 'absorbed by the ground', value: number(light.ground) },
    ], x0, x0 + span, size, yBar - 7, yBar - 13, yBar - 34);
    const reach = Math.min(span, mine * zoom), squeeze = mine > 0 ? reach / (mine * zoom) : 1, yGap = yBar + bh, mid = (yGap + yZoom) / 2;
    ctx.fillStyle = 'rgba(255, 255, 255, 0.04)'; ctx.beginPath(); ctx.moveTo(x0, yGap); ctx.lineTo(xAir, yGap); ctx.lineTo(x0 + reach, yZoom); ctx.lineTo(x0, yZoom); ctx.closePath(); ctx.fill();
    ctx.strokeStyle = LINE; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(xAir, yGap + 2); ctx.lineTo(x0 + reach, yZoom - 2); ctx.stroke(); ctx.setLineDash([]);
    text(ctx, 'three times larger', x0 + span, mid, { align: 'right', color: MUTED, size });
    segments(ctx, x0, yZoom, bh, zoom * squeeze, [[vapor, VAPOR], [ozone, OZONE], [gases, GAS], [aerosol, HAZE], [byAir, SKY]]);
    let at = x0;
    const items = [['vapor', vapor, light.vapor], ['ozone', ozone, light.ozone], ['O₂ + CO₂', gases, light.gases], ['aerosol', aerosol, light.aerosol], ['reflected by the air', byAir, light.byAir]].filter(([, , exact]) => exact > 0.005).map(([name, value, exact]) => {
      const width = Math.max(0, value) * zoom * squeeze, item = { anchor: at + width / 2, name, value: number(exact) };
      at += width;
      return item;
    });
    labels(ctx, items, x0, x0 + span, size, yZoom + bh + 1, yZoom + bh + 8, yZoom + bh + 15);
  }

  function drawProfile(ctx, w, h, size) {
    const narrow = w < WIDE, left = narrow ? 56 : 62, right = w - (narrow ? 40 : 46), top = BARS + 30, bottom = h - 34, rate = light.rate;
    const ip = light.interfacePressure, ps = ip[K], span = Math.log(ps) - Math.log(TOP);
    const y = (p) => top + (Math.log(p) - Math.log(TOP)) / span * (bottom - top), x = (r) => left + r / rate * (right - left);
    const at = (z) => { const { height: zs, pressure: p } = light; if (z <= zs[K - 1]) return ps * Math.exp(-z / zs[K - 1] * Math.log(ps / p[K - 1])); for (let k = K - 1; k > 0; k--) if (zs[k - 1] >= z) return Math.exp(Math.log(p[k]) + (z - zs[k]) / (zs[k - 1] - zs[k]) * (Math.log(p[k - 1]) - Math.log(p[k]))); return p[0]; };
    plot = { y, top, bottom, left, right };
    text(ctx, 'heating by sunlight, K per day', left, BARS + 10, { color: INK, size: size + 0.5, weight: 500 });
    ctx.lineWidth = 1;
    for (let r = 0; r <= rate; r += light.tick) { ctx.strokeStyle = r ? GRID : LINE; ctx.beginPath(); ctx.moveTo(x(r), top); ctx.lineTo(x(r), bottom); ctx.stroke(); text(ctx, `${r}`, x(r), bottom + 16, { align: 'center', color: MUTED, size }); }
    for (const p of PRESSURES) { ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(left, y(p * 100)); ctx.lineTo(right, y(p * 100)); ctx.stroke(); text(ctx, `${p} hPa`, left - 6, y(p * 100), { align: 'right', color: MUTED, size }); }
    for (const z of HEIGHTS) { const yy = y(at(z * 1000)); ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(right, yy); ctx.lineTo(right + 4, yy); ctx.stroke(); text(ctx, `${z} km`, right + 7, yy, { color: MUTED, size }); }
    ctx.fillStyle = EARTH; ctx.fillRect(left, bottom, right - left, 5);
    const r = lit.heating, d = light.dry, stairs = (values) => { ctx.moveTo(x(values[1]), y(ip[1])); for (let k = 1; k < K; k++) { ctx.lineTo(x(values[k]), y(ip[k])); ctx.lineTo(x(values[k]), y(ip[k + 1])); } };
    ctx.save(); ctx.beginPath(); ctx.rect(left, top, right - left, bottom - top); ctx.clip();
    ctx.beginPath(); ctx.moveTo(x(0), y(ip[1]));
    for (let k = 1; k < K; k++) { ctx.lineTo(x(r[k]), y(ip[k])); ctx.lineTo(x(r[k]), y(ip[k + 1])); }
    ctx.lineTo(x(0), y(ps)); ctx.closePath();
    ctx.globalAlpha = 0.18; ctx.fillStyle = A; ctx.fill(); ctx.globalAlpha = 1;
    if (humidity > 0) { ctx.beginPath(); stairs(d); ctx.strokeStyle = MUTED; ctx.lineWidth = 1.5; ctx.setLineDash([3, 4]); ctx.stroke(); ctx.setLineDash([]); }
    ctx.beginPath(); stairs(r); ctx.strokeStyle = A; ctx.lineWidth = 2; ctx.lineJoin = 'round'; ctx.stroke();
    ctx.restore();
    const low = light.lower;
    text(ctx, `${low.low.toFixed(2)} to ${low.high.toFixed(2)} K per day in the lowest 10 km`, x(low.high) + 8, y(at(4000)), { color: MUTED, size, halo: HALO });
    const mark = (k, strong) => {
      const px = x(r[k]), py = y(light.pressure[k]), value = light.heating[k], rateText = `${value.toFixed(value < 1 ? 2 : 1)} K/day at ${hPa(light.pressure[k])} hPa`;
      const was = humidity > 0 && !strong ? `, ${d[k].toFixed(d[k] < 1 ? 2 : 1)} in dry air` : '';
      ctx.font = `400 ${size}px system-ui, sans-serif`;
      const options = [`${rateText}${was}`, `${rateText}, ${(light.height[k] / 1000).toFixed(0)} km up`, rateText].flatMap((words) => [[words, px + 10, 'left'], [words, px - 10, 'right']]);
      const fits = ([words, from, align]) => { const width = ctx.measureText(words).width; return align === 'left' ? from + width <= right - 2 : from - width >= left + 2; };
      const [words, from, align] = options.find(fits) ?? [rateText, clamp(px - ctx.measureText(rateText).width / 2, left + 2, right - 2 - ctx.measureText(rateText).width), 'left'];
      if (!strong) { ctx.strokeStyle = 'rgba(255, 255, 255, 0.3)'; ctx.setLineDash([2, 4]); ctx.beginPath(); ctx.moveTo(left, py); ctx.lineTo(right, py); ctx.stroke(); ctx.setLineDash([]); }
      ctx.fillStyle = strong ? A : INK; ctx.strokeStyle = HALO; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.arc(px, py, 4, 0, Math.PI * 2); ctx.stroke(); ctx.fill();
      text(ctx, words, from, options.some(fits) ? py : py + 14, { align, color: INK, size, halo: HALO });
    };
    if (layer !== null && layer !== light.peak.k) mark(layer, false);
    mark(light.peak.k, true);
  }

  function draw(ctx, w, h) {
    const size = w < WIDE ? 10.5 : 11;
    if (view === 'slices' ? !state || !shown : !light || !lit) { text(ctx, 'Loading the midlatitude summer…', w / 2, h / 2, { align: 'center', color: MUTED, size: 12 }); return; }
    if (view === 'slices') { drawSlices(ctx, w, h); return; }
    drawBars(ctx, w, size);
    drawProfile(ctx, w, h, size);
  }

  const pickSlice = (ev) => {
    if (!geometry || !state) return;
    const r = fig.canvas.getBoundingClientRect(), px = ev.clientX - r.left, { x0, x1 } = geometry;
    if (px < x0 - 10 || px > x1 + 10) return;
    const f = clamp((px - x0) / (x1 - x0), 0, 1 - 1e-9);
    let at = 0, best = state.order[0];
    for (const j of state.order) { const share = state.points[j].share; if (share > 1e-3) best = j; if (at + share > f && share > 1e-3) { best = j; break; } at += share; }
    touched = true;
    if (best !== probe) { probe = best; report(); }
    fig.render();
  };
  const pickLayer = (ev) => {
    if (!plot || !light) return;
    const r = fig.canvas.getBoundingClientRect(), py = ev.clientY - r.top, px = ev.clientX - r.left;
    if (py < plot.top - 8 || py > plot.bottom + 8 || px < plot.left - 40) { if (layer !== null) { layer = null; fig.render(); } return; }
    let best = 1;
    for (let k = 1; k < K; k++) if (Math.abs(plot.y(light.pressure[k]) - py) < Math.abs(plot.y(light.pressure[best]) - py)) best = k;
    if (best !== layer) { layer = best; fig.render(); }
  };
  const pick = (ev) => (view === 'slices' ? pickSlice(ev) : pickLayer(ev));
  fig.canvas.style.touchAction = 'pan-y';
  fig.canvas.addEventListener('pointerdown', pick);
  fig.canvas.addEventListener('pointermove', (ev) => { if (ev.pointerType === 'mouse' || (view === 'slices' && ev.buttons)) pick(ev); });
  fig.canvas.addEventListener('pointerleave', (ev) => { if (view !== 'slices' && ev.pointerType === 'mouse' && layer !== null) { layer = null; fig.render(); } });

  arrange();
  update();
  loadColumn().then((column) => { profile = column; split = createSplit(column); update(); });
}
