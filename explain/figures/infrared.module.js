import { Figure, choice, legend, caption, readout, text, clamp, INK, MUTED, GRID, LINE, ACCENT } from '../runtime.module.js';
import { loadColumn } from '../radiativeColumn.module.js';
import { slices, GASES } from '../infraredSlices.module.js';

const WIDE = 520, LABELED = 0.02, P_TOP = 10, T_LOW = 180, T_HIGH = 300, RATE = 7, HALO = 'rgba(20, 20, 22, 0.85)', GROUND = '#5a4a36';
const COLORS = { carbonDioxide: '#d55e00', ozone: '#0072b2', methane: '#009e73', nitrousOxide: '#cc79a7', none: '#5c5c5c' };
const NAMES = { carbonDioxide: 'carbon dioxide', ozone: 'ozone', methane: 'methane', nitrousOxide: 'nitrous oxide' };
const FORMULAS = { carbonDioxide: 'CO₂', ozone: 'O₃', methane: 'CH₄', nitrousOxide: 'N₂O' };
const AMOUNTS = [['none', '0'], ['×½', '0.5'], ['×1', '1'], ['×2', '2'], ['×4', '4']];
const RGB = Object.fromEntries(Object.entries(COLORS).map(([k, hex]) => [k, [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16))]));
const rgb = (c, a = 1) => `rgba(${Math.round(c[0])}, ${Math.round(c[1])}, ${Math.round(c[2])}, ${a})`;
const watts = (v) => `${(Math.abs(v) < 0.05 ? 0 : v).toFixed(1)} W/m²`;
const percent = (f) => `${(100 * f).toFixed(1)}%`;
const depthText = (d) => (d >= 100 ? d.toExponential(1).replace('e+', ' × 10^') : d >= 10 ? d.toFixed(0) : d >= 1 ? d.toFixed(1) : d.toFixed(d >= 0.01 ? 2 : 4));
const SUPER = { 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
const pretty = (s) => s.replace(/\^(\d+)/, (_, n) => [...n].map((c) => SUPER[c]).join(''));
const pressureText = (p) => (p >= 10000 ? `${(p / 100).toFixed(0)} hPa` : p >= 1000 ? `${(p / 100).toFixed(1)} hPa` : `${(p / 100).toFixed(2)} hPa`);

export function mountInfrared(root) {
  const controls = root.querySelector('.controls');
  let profile = null, state = null, shown = null, probe = null, touched = false, geometry = null, amount = 1, lastAmount = 1;
  const on = Object.fromEntries(GASES.map((g) => [g, true]));
  const fig = new Figure(root, { height: (w) => (w < WIDE ? 420 : 440), step, draw });
  caption(root, 'The midlatitude summer atmosphere with its water vapor taken out, over ground at 294.2 K. Along the bottom, the ground’s infrared is split into the model’s 34 slices, each as wide as its share, the clearest on the left. Each bar rises to where its slice leaves for space, the level at which the gases above it reach an optical depth of 1; where there is no bar, the ground shines through. The row along the top is what each slice sends to space. On the right is the air’s temperature. ×1 carbon dioxide is 388.75 ppm. Drag across the slices to read one.');
  legend(root, [...GASES.map((g) => ['bar', `${NAMES[g]}, ${FORMULAS[g]}`, COLORS[g]]), ['bar', 'nearly clear: under 5% absorbed', COLORS.none], ['line', 'temperature of the air', INK]]);

  const group = document.createElement('div'), label = document.createElement('span');
  group.className = 'choice'; group.setAttribute('role', 'group'); group.style.gridColumn = '1 / -1';
  label.className = 'caption'; label.textContent = 'Gases'; label.id = 'infrared-gases'; group.setAttribute('aria-labelledby', label.id);
  group.append(label);
  const toggles = Object.fromEntries(GASES.map((g) => {
    const b = document.createElement('button');
    b.type = 'button'; b.textContent = NAMES[g];
    b.addEventListener('click', () => {
      if (g === 'carbonDioxide') { amount = on[g] ? 0 : lastAmount; on[g] = amount > 0; dose.value = String(amount); } else on[g] = !on[g];
      update();
    });
    group.append(b);
    return [g, b];
  }));
  controls.append(group);
  const dose = choice(controls, { label: 'Carbon dioxide', options: AMOUNTS, value: '1', onChange: (v) => { amount = Number(v); if (amount > 0) lastAmount = amount; on.carbonDioxide = amount > 0; update(); }, span: true });
  const totals = readout(controls), detail = readout(controls);

  function update() {
    for (const g of GASES) toggles[g].setAttribute('aria-pressed', String(on[g]));
    if (!profile) return;
    state = slices(profile, { carbonDioxide: amount, gases: { ...on, carbonDioxide: amount > 0 } });
    if (!shown) shown = state.points.map(() => ({ left: 0, level: 0, kept: 1, color: [0, 0, 0] }));
    if (probe === null) probe = state.points.reduce((best, p) => (p.dominant === 'carbonDioxide' && p.level !== null && (best === null || p.share > state.points[best].share) ? p.index : best), null);
    report();
    fig.play(true);
    fig.render();
  }

  function targets() {
    const ps = state.surfacePressure, out = [];
    let left = 0;
    for (const j of state.order) {
      const p = state.points[j];
      out[j] = { left, level: Math.log((p.level ?? ps) / P_TOP), kept: p.ground > 0 ? p.up / p.ground : 1, color: RGB[p.dominant ?? 'none'] };
      left += p.share;
    }
    return out;
  }

  let first = true;
  function step(dt) {
    if (!state) return false;
    const goal = targets(), f = first ? 1 : 1 - Math.exp(-RATE * dt);
    let moving = false;
    shown.forEach((s, j) => {
      const g = goal[j];
      s.left += (g.left - s.left) * f; s.level += (g.level - s.level) * f; s.kept += (g.kept - s.kept) * f;
      for (let c = 0; c < 3; c++) s.color[c] += (g.color[c] - s.color[c]) * f;
      if (Math.abs(g.left - s.left) > 1e-4 || Math.abs(g.level - s.level) > 1e-3 || Math.abs(g.kept - s.kept) > 1e-4 || Math.abs(g.color[0] - s.color[0]) + Math.abs(g.color[1] - s.color[1]) + Math.abs(g.color[2] - s.color[2]) > 0.5) moving = true;
      else { s.left = g.left; s.level = g.level; s.kept = g.kept; s.color = [...g.color]; }
    });
    first = false;
    if (!moving) fig.play(false);
    return true;
  }

  function report() {
    totals.set([['the ground emits', watts(state.groundEmission)], ['leaving the top', watts(state.olr)], ['held back, G', watts(state.greenhouse)], ['back down to the ground', watts(state.dlr)], ['straight to space', `${percent(state.escaping)} of the ground’s emission`]]);
    const p = state.points[probe], name = NAMES[p.dominant];
    const where = p.level === null ? [['optical depth', depthText(p.depth)], ['straight to space', `${percent(Math.exp(-p.depth))} of it`]] : [['optical depth', pretty(depthText(p.depth))], ['leaves for space from', `${pressureText(p.level)}, ${(p.height / 1000).toFixed(1)} km, where the air is at ${p.temperature.toFixed(1)} K`]];
    detail.set([[`slice ${p.index + 1} of ${state.points.length}`, `${p.dominant ? `mostly ${name}` : 'nearly clear'}, ${percent(p.share)} of the ground’s emission`], ...where, ['to space', `${watts(p.up)} of the ${watts(p.ground)} the ground emits in it`]]);
  }

  function draw(ctx, w, h) {
    if (!state || !shown) { text(ctx, 'Loading the midlatitude summer…', w / 2, h / 2, { align: 'center', color: MUTED, size: 12 }); return; }
    if (first) step(0);
    const narrow = w < WIDE, ps = state.surfacePressure, left = narrow ? 46 : 52, right = w - (narrow ? 8 : 14), panel = narrow ? 66 : 120, gap = narrow ? 12 : 24;
    const x0 = left, x1 = right - panel - gap, t0 = x1 + gap, t1 = right, top = 50, bottom = h - 54, row = 10, rowTop = 26;
    const span = Math.log(ps / P_TOP), y = (lp) => top + lp / span * (bottom - top), yp = (p) => y(Math.log(p / P_TOP));
    const x = (f) => x0 + f * (x1 - x0), tx = (T) => t0 + (T - T_LOW) / (T_HIGH - T_LOW) * (t1 - t0);
    geometry = { x0, x1 };

    ctx.lineWidth = 1;
    for (const hPa of [0.1, 1, 10, 100, 1000]) {
      const yy = yp(hPa * 100);
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(x0 - 4, yy); ctx.lineTo(x1, yy); ctx.moveTo(t0, yy); ctx.lineTo(t1, yy); ctx.stroke();
      text(ctx, `${hPa}${hPa === 0.1 ? ' hPa' : ''}`, x0 - 7, yy, { align: 'right', color: MUTED, size: 10 });
    }
    for (const T of [200, 250, 300]) {
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(tx(T), top); ctx.lineTo(tx(T), bottom); ctx.stroke();
      if (!narrow || T !== 250) text(ctx, `${T}`, tx(T), bottom + 22, { align: T === 300 ? 'right' : 'center', color: MUTED, size: 10 });
    }
    text(ctx, narrow ? 'air, K' : 'air temperature, K', (t0 + t1) / 2, bottom + 38, { align: 'center', color: MUTED, size: 10 });
    for (const f of [0, 0.25, 0.5, 0.75, 1]) text(ctx, f === 1 ? '100%' : `${100 * f}`, x(f), bottom + 22, { align: f === 1 ? 'right' : f === 0 ? 'left' : 'center', color: MUTED, size: 10 });
    text(ctx, `the ground’s ${watts(state.groundEmission)}, slice by slice`, (x0 + x1) / 2, bottom + 38, { align: 'center', color: MUTED, size: 10 });
    text(ctx, `to space, ${watts(state.olr)}`, x0, rowTop - 10, { color: INK, size: 11 });
    text(ctx, `held back, ${watts(state.greenhouse)}`, t1, rowTop - 10, { align: 'right', color: MUTED, size: 11 });

    const order = state.order.slice().sort((a, b) => shown[a].left - shown[b].left);
    for (const j of order) {
      const s = shown[j], p = state.points[j], a = x(s.left), b = x(s.left + p.share), width = b - a;
      if (width < 0.05) continue;
      const inner = width >= 3 ? width - 1 : width, color = s.color, yy = Math.min(bottom, y(s.level));
      ctx.fillStyle = rgb(color, 0.55); ctx.fillRect(a, yy, inner, bottom - yy);
      if (bottom - yy > 0.5) { ctx.fillStyle = rgb(color); ctx.fillRect(a, yy, inner, Math.min(2, bottom - yy)); }
      ctx.fillStyle = rgb(color); ctx.fillRect(a, bottom + 2, inner, row);
      const kept = clamp(s.kept, 0, 1) * inner;
      ctx.fillRect(a + (inner - kept) / 2, rowTop, kept, row);
    }
    ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(x0, bottom + 0.5); ctx.lineTo(x1, bottom + 0.5); ctx.stroke();

    let clear = 0;
    for (const p of state.points) if (p.level === null) clear += p.share;
    if (clear > 0.3) text(ctx, narrow ? 'the ground shines through' : 'mostly clear: the ground shines through', x(clear / 2), yp(20000), { align: 'center', color: MUTED, size: 11 });
    const groups = [];
    order.forEach((j, n) => {
      const p = state.points[j], last = groups.at(-1);
      if (!p.dominant) return;
      if (last && last.gas === p.dominant && last.next === n) { last.share += p.share; last.next = n + 1; if (p.share > last.widest.share) last.widest = p; } else groups.push({ gas: p.dominant, share: p.share, widest: p, next: n + 1 });
    });
    const labels = [], wide = narrow ? 24 : 28;
    if (probe !== null) { const s = shown[probe], a = x(s.left), b = Math.max(a + 2, x(s.left + state.points[probe].share)); for (let yy = Math.min(bottom, y(s.level)); yy < bottom; yy += 13) labels.push([a, yy], [b, yy]); }
    for (const { gas, widest } of groups.filter((group) => group.share >= LABELED).sort((a, b) => b.share - a.share)) {
      const s = shown[widest.index], cx = clamp(x(s.left + widest.share / 2), x0 + 12, x1 - 12);
      let yy = Math.min(bottom, y(s.level)) - 9;
      while (labels.some(([lx, ly]) => Math.abs(lx - cx) < wide && Math.abs(ly - yy) < 13)) yy -= 13;
      labels.push([cx, yy]);
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
      const T = p.level === null ? state.groundTemperature : p.temperature, px = tx(T);
      ctx.strokeStyle = 'rgba(255,255,255,0.5)'; ctx.lineWidth = 1; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(b + 2, yy); ctx.lineTo(px, yy); ctx.stroke(); ctx.setLineDash([]);
      ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(px, yy, 3.5, 0, Math.PI * 2); ctx.fill();
      const tag = `${T.toFixed(0)} K`, low = yy > bottom - 14;
      text(ctx, tag, px + (T > 270 ? -7 : 7), low ? yy - 11 : yy - 10, { align: T > 270 ? 'right' : 'left', color: '#fff', size: 11, weight: 500, halo: HALO });
      if (!touched) text(ctx, 'drag', clamp((a + b) / 2, x0 + 14, x1 - 14), Math.min(yy, bottom - 4) - 24, { align: 'center', color: ACCENT, size: 11, halo: HALO });
    }
  }

  const pick = (ev) => {
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
  fig.canvas.style.touchAction = 'pan-y';
  fig.canvas.addEventListener('pointerdown', pick);
  fig.canvas.addEventListener('pointermove', (ev) => { if (ev.pointerType === 'mouse' || ev.buttons) pick(ev); });

  update();
  loadColumn().then((column) => { profile = column; update(); });
}
