import { Figure, slider, choice, buttons, readout, legend, caption, text, arrow, termColor, clamp, INK, MUTED, GRID, LINE, ACCENT } from '../runtime.module.js';
import { K } from '../radiativeColumn.module.js';
import { createPair, tropopause, logScale, pickLayers, daysPerSecond, STEPS_PER_DAY, PPM } from '../radiativeConvective.module.js';
import { GRAVITY, CP_DRY, R_DRY, P0 } from '../../js/dynamics/sigmaCore.module.js';

const KAPPA = R_DRY / CP_DRY, LAPSE = GRAVITY / CP_DRY, TOP_PRESSURE = 10, LOW_TOP = 3000, WIDE = 520;
const SMALL = 0.03, LARGE = 40, BUDGET = 8, MOST = 75, GHOST = 'rgba(255,255,255,0.4)', HALO = 'rgba(20, 20, 22, 0.85)';
const PRESSURES = [1000, 500, 200, 100, 50, 20, 10, 5, 2, 1, 0.5, 0.2, 0.1];
const TERMS = [['sunlight', 'a'], ['infrared', 'b'], ['ground', 'c'], ['convection', 'd']];
const signed = (v, d = 2) => (Math.abs(v) < 0.5 * 10 ** -d ? (0).toFixed(d) : `${v < 0 ? '−' : '+'}${Math.abs(v).toFixed(d)}`);
const rate = (v) => signed(v, Math.abs(v) >= 9.995 ? 1 : 2);
const plain = (v, d = 1) => (Math.abs(v) < 0.5 * 10 ** -d ? (0).toFixed(d) : `${v < 0 ? '−' : ''}${Math.abs(v).toFixed(d)}`);
const hPa = (p) => { const v = p / 100; return v >= 9.95 ? v.toFixed(0) : v >= 0.995 ? v.toFixed(1) : v.toFixed(2); };
const meters = (z) => (z < 995 ? `${Math.round(z / 10) * 10} m` : `${(z / 1000).toFixed(1)} km`);

export function mountEquilibrium(root) {
  const controls = root.querySelector('.controls'), color = Object.fromEntries([...'abcde'].map((c) => [c, termColor(c)]));
  const pair = createPair();
  let view = 'whole', running = true, probe = K - 1, touched = false, owed = 0, range = null, goal = null, geometry = null;

  const fig = new Figure(root, { height: (w) => (w < WIDE ? 480 : 470), step, draw });
  caption(root, 'The model’s own radiation and dry convection in one column of 36 layers over dark ground, under the global average of sunlight, stepped every 3 simulated hours from 288 K everywhere. On the left is the temperature of each layer and, on the brown strip, of the ground. On the right is how fast each process warms or cools a dozen of the layers: each tenfold step in the heating adds the same length to an arrow, so that 0.1 and 10 K a day both show. Tap or drag to read a layer. The run pauses itself once the column stops changing.');
  legend(root, [['line', 'temperature', INK], ['line', 'with radiation alone', GHOST], ['dash', 'dry adiabat', LINE], ['dots', 'tropopause', 'rgba(255,255,255,0.7)'], ['arrow', 'sunlight', color.a], ['arrow', 'infrared', color.b], ['arrow', 'heat from the ground', color.c], ['arrow', 'convection', color.d], ['arrow', 'all together', color.e]]);
  controls.classList.add('two');
  choice(controls, { label: 'Show', options: [['whole column', 'whole'], ['lowest 3 km', 'low']], value: view, onChange: (v) => { view = v; if (view === 'low' && pair.column.height[probe] > LOW_TOP) probe = K - 1; settle(); }, span: true });
  choice(controls, { label: 'Convection', options: [['off', 'off'], ['on', 'on']], value: 'on', onChange: (v) => { pair.set({ convection: v === 'on' }); resume(); } });
  const [, toggle] = buttons(controls, [['Restart', () => { pair.restart(); resume(); }], ['Pause', () => { if (running) { halt(); settle(); } else { if (pair.settled()) pair.since = 0; resume(); } }]]);
  slider(controls, { label: 'Carbon dioxide', min: 0, max: 4, step: 0.25, value: 1, format: (v) => `×${v}, ${Math.round(v * PPM)} ppm`, onInput: (v) => { pair.set({ carbonDioxide: v }); resume(); }, span: true });
  slider(controls, { label: 'Ground reflects', min: 0.07, max: 0.4, step: 0.01, value: 0.07, format: (v) => `${Math.round(v * 100)}% of sunlight`, onInput: (v) => { pair.set({ albedo: v }); resume(); }, span: true });
  const out = readout(controls), layerOut = readout(controls);
  for (const box of controls.querySelectorAll('.readout')) box.style.gridColumn = '1 / -1';

  function resume() { running = true; toggle.textContent = 'Pause'; settle(); }
  function halt() { running = false; toggle.textContent = 'Run'; }
  function settle() { goal = target(); fig.play(running || !range || !same(range, goal)); fig.render(); }
  const same = (a, b) => Math.abs(a[0] - b[0]) < 0.05 && Math.abs(a[1] - b[1]) < 0.05;

  function target() {
    let lo = Infinity, hi = -Infinity;
    for (const c of [pair.column, pair.ghost]) {
      for (let k = 0; k < K; k++) if (view === 'whole' || c.height[k] <= LOW_TOP * 1.05) { lo = Math.min(lo, c.T[k]); hi = Math.max(hi, c.T[k]); }
      lo = Math.min(lo, c.Ts); hi = Math.max(hi, c.Ts);
    }
    const unit = view === 'whole' ? 10 : 5, pad = view === 'whole' ? 6 : 2, least = view === 'whole' ? 60 : 20;
    lo = Math.floor((lo - pad) / unit) * unit; hi = Math.ceil((hi + pad) / unit) * unit;
    if (hi - lo < least) { const mid = (hi + lo) / 2; lo = Math.floor((mid - least / 2) / unit) * unit; hi = lo + least; }
    return [lo, hi];
  }

  function step(dt) {
    if (running) {
      owed = Math.min(owed + dt * daysPerSecond(pair.since) * STEPS_PER_DAY, MOST);
      const start = performance.now();
      for (; owed >= 1 && performance.now() - start < BUDGET; owed--) pair.step();
      owed = Math.min(owed, 1);
      if (pair.settled()) halt();
    }
    goal = target();
    range = !range ? goal : same(range, goal) ? goal : range.map((v, n) => v + (goal[n] - v) * Math.min(1, dt * 6));
    if (!running && same(range, goal)) { range = goal; fig.play(false); }
  }

  function layout(w, h) {
    const narrow = w < WIDE, left = narrow ? 34 : 44, right = w - (narrow ? 6 : 10), top = 32, bottom = h - 46, gap = narrow ? 14 : 26;
    const mid = left + (right - left - gap) * 0.56, r0 = mid + gap, cx = (r0 + right) / 2, half = (right - r0) / 2 - 2;
    const span = Math.log(pair.column.interfacePressure[K] / TOP_PRESSURE);
    const yOf = view === 'whole'
      ? (c, k) => bottom - Math.log(c.interfacePressure[K] / c.pressure[k]) / span * (bottom - top)
      : (c, k) => bottom - c.height[k] / LOW_TOP * (bottom - top);
    const yP = (p) => bottom - Math.log(pair.column.interfacePressure[K] / p) / span * (bottom - top);
    const yZ = (z) => bottom - z / LOW_TOP * (bottom - top);
    return { narrow, left, right, top, bottom, mid, r0, cx, half, yOf, yP, yZ };
  }

  function draw(ctx, w, h) {
    range ??= target();
    const g = geometry = layout(w, h), { narrow, left, top, bottom, mid, r0, right, cx, half, yOf, yP, yZ } = g, { column, ghost } = pair;
    const [lo, hi] = range, x = (T) => left + (T - lo) / (hi - lo) * (mid - left);
    ctx.lineWidth = 1;
    if (view === 'whole') {
      text(ctx, 'hPa', left - 6, top - 14, { align: 'right', color: MUTED, size: 10 });
      for (const p of PRESSURES) {
        const y = yP(p * 100);
        ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(left - 4, y); ctx.lineTo(mid, y); ctx.moveTo(r0, y); ctx.lineTo(right, y); ctx.stroke();
        text(ctx, `${p}`, left - 7, y, { align: 'right', color: MUTED, size: 10 });
      }
    } else {
      text(ctx, 'km', left - 6, top - 14, { align: 'right', color: MUTED, size: 10 });
      for (let z = 0; z <= LOW_TOP; z += 500) {
        const y = yZ(z);
        ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(left - 4, y); ctx.lineTo(mid, y); ctx.moveTo(r0, y); ctx.lineTo(right, y); ctx.stroke();
        text(ctx, `${z / 1000}`, left - 7, y, { align: 'right', color: MUTED, size: 10 });
      }
    }
    const every = (hi - lo) / (narrow ? 4 : 6), tick = [5, 10, 20, 40, 50, 100].find((s) => s >= every) ?? 100;
    for (let T = Math.ceil(lo / tick) * tick; T <= hi + 1e-9; T += tick) {
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(x(T), top); ctx.lineTo(x(T), bottom); ctx.stroke();
      text(ctx, `${T}`, x(T), bottom + 16, { align: 'center', color: MUTED, size: 10 });
    }
    text(ctx, 'temperature, K', (left + mid) / 2, bottom + 31, { align: 'center', color: MUTED, size: 10 });

    ctx.save();
    ctx.beginPath(); ctx.rect(left, top - 6, mid - left, bottom - top + 6); ctx.clip();
    const tp = tropopause(column);
    if (tp) {
      const y = yOf(column, tp.k);
      ctx.strokeStyle = 'rgba(255,255,255,0.7)'; ctx.lineWidth = 2; ctx.lineCap = 'round'; ctx.setLineDash([0.5, 4]);
      ctx.beginPath(); ctx.moveTo(left, y); ctx.lineTo(mid, y); ctx.stroke(); ctx.setLineDash([]); ctx.lineCap = 'butt';
    }
    ctx.strokeStyle = LINE; ctx.lineWidth = 1.2; ctx.setLineDash([4, 4]); ctx.beginPath();
    if (view === 'whole') {
      const theta = column.theta[K - 1], ps = column.interfacePressure[K];
      for (let n = 0; n <= 80; n++) { const p = ps * (TOP_PRESSURE / ps) ** (n / 80), px = x(theta * (p / P0) ** KAPPA), py = yP(p); if (n) ctx.lineTo(px, py); else ctx.moveTo(px, py); }
    } else {
      const z0 = column.height[K - 1], T0 = column.T[K - 1];
      ctx.moveTo(x(T0 + LAPSE * z0), yZ(0)); ctx.lineTo(x(T0 - LAPSE * (LOW_TOP * 1.1 - z0)), yZ(LOW_TOP * 1.1));
    }
    ctx.stroke(); ctx.setLineDash([]);
    profile(ctx, ghost, x, yOf, GHOST, 1.5, false);
    profile(ctx, column, x, yOf, INK, 2, true);
    ctx.restore();

    ctx.fillStyle = '#5a4a36'; ctx.fillRect(left, bottom, mid - left, 5);
    ctx.lineWidth = 1.5;
    ctx.strokeStyle = GHOST; ctx.fillStyle = HALO; ctx.beginPath(); ctx.arc(x(ghost.Ts), bottom + 2.5, 4, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
    ctx.fillStyle = INK; ctx.beginPath(); ctx.arc(x(column.Ts), bottom + 2.5, 3.5, 0, Math.PI * 2); ctx.fill();

    const ticks = narrow ? [-10, -1, 0, 1, 10] : [-10, -1, -0.1, 0, 0.1, 1, 10];
    for (const v of ticks) {
      const px = cx + logScale(v, SMALL, LARGE) * half;
      ctx.strokeStyle = v ? GRID : LINE; ctx.beginPath(); ctx.moveTo(px, top - 4); ctx.lineTo(px, bottom); ctx.stroke();
      text(ctx, v ? `${v < 0 ? '−' : ''}${Math.abs(v)}` : '0', px, bottom + 16, { align: 'center', color: MUTED, size: 10 });
    }
    text(ctx, narrow ? 'K per day' : 'heating, K per day', cx, bottom + 31, { align: 'center', color: MUTED, size: 10 });
    text(ctx, 'cools', cx - 6, top - 14, { align: 'right', color: MUTED, size: 10 });
    text(ctx, 'warms', cx + 6, top - 14, { color: MUTED, size: 10 });

    const spacing = narrow ? 4 : 4.5, rows = (k) => TERMS.filter(([name]) => Math.abs(column.heating[name][k]) >= 0.005).length + 1;
    const place = (k) => { const n = rows(k), y = yOf(column, k); return clamp(y, top + (n - 1) * spacing / 2 + 2, bottom - (n - 1) * spacing / 2 - 3); };
    const centers = Array.from({ length: K }, (_, k) => (view === 'low' && column.height[k] > LOW_TOP ? NaN : place(k)));
    const shown = pickLayers(centers, narrow ? 20 : 22, [probe, K - 1]);
    const py = centers[probe], layerY = yOf(column, probe);
    ctx.strokeStyle = 'rgba(255,255,255,0.3)'; ctx.lineWidth = 1; ctx.setLineDash([2, 4]); ctx.beginPath(); ctx.moveTo(left, layerY); ctx.lineTo(mid, layerY); ctx.lineTo(r0, py); ctx.lineTo(right, py); ctx.stroke(); ctx.setLineDash([]);
    for (const k of shown) {
      const items = [...TERMS.filter(([name]) => Math.abs(column.heating[name][k]) >= 0.005).map(([name, c]) => [column.heating[name][k], color[c], 1.6]), [column.heating.total[k], color.e, 2]];
      const y0 = centers[k] - (items.length - 1) * spacing / 2;
      ctx.globalAlpha = k === probe ? 1 : 0.55;
      items.forEach(([v, c, width], n) => {
        const y = y0 + n * spacing, len = logScale(v, SMALL, LARGE) * half;
        if (Math.abs(len) >= 2.5) arrow(ctx, cx, y, cx + len, y, { color: c, width, head: 5 });
        else { ctx.fillStyle = c; ctx.beginPath(); ctx.arc(cx, y, 1.6, 0, Math.PI * 2); ctx.fill(); }
      });
      ctx.globalAlpha = 1;
    }
    const ry = yOf(column, probe), rx = x(column.T[probe]);
    if (ry >= top - 6 && ry <= bottom + 1) {
      ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(rx, ry, 6, 0, Math.PI * 2); ctx.stroke();
      if (!touched) { const before = rx > (left + mid) / 2; text(ctx, 'drag', before ? rx - 10 : rx + 10, ry - 11, { align: before ? 'right' : 'left', color: ACCENT, size: 11, halo: HALO }); }
    }
    text(ctx, `day ${Math.floor(column.time / 86400).toLocaleString('en-US')}`, mid, top - 14, { align: 'right', color: MUTED, size: 10 });
    report(fig.running ? 4 : 1, tp);
  }

  function profile(ctx, c, x, yOf, stroke, width, dots) {
    ctx.strokeStyle = stroke; ctx.lineWidth = width; ctx.lineJoin = 'round'; ctx.beginPath();
    for (let k = 0; k < K; k++) { const px = x(c.T[k]), py = yOf(c, k); if (k) ctx.lineTo(px, py); else ctx.moveTo(px, py); }
    ctx.stroke();
    if (!dots) return;
    ctx.fillStyle = stroke;
    for (let k = 0; k < K; k++) { ctx.beginPath(); ctx.arc(x(c.T[k]), yOf(c, k), 1.8, 0, Math.PI * 2); ctx.fill(); }
  }

  function report(every, tp) {
    const { column, ghost } = pair, b = column.budget, day = Math.floor(column.time / 86400).toLocaleString('en-US');
    out.set([
      ['day', `${day}${!running && pair.settled() ? ', settled' : ''}`],
      ['ground', `${column.Ts.toFixed(1)} K, or ${ghost.Ts.toFixed(1)} K with radiation alone`],
      ['sunlight absorbed', `${b.absorbed.toFixed(1)} W/m²`],
      ['infrared leaving the top', `${b.olr.toFixed(1)} W/m²`],
      ['imbalance', `${signed(b.imbalance)} W/m²`],
      ['held back, G', `${plain(b.greenhouse)} W/m²`],
      ['tropopause', tp ? `${hPa(tp.pressure)} hPa, ${meters(tp.height)}, ${tp.temperature.toFixed(1)} K` : 'none yet'],
    ], every);
    const k = probe, names = { sunlight: 'sunlight', infrared: 'infrared', ground: 'heat from the ground', convection: 'convection' };
    layerOut.set([
      ['layer', `${K - k} of ${K}, at ${hPa(column.pressure[k])} hPa and ${meters(column.height[k])}, ${column.T[k].toFixed(1)} K`],
      ...TERMS.filter(([name]) => column.heating[name][k] !== 0).map(([name]) => [names[name], rate(column.heating[name][k])]),
      ['all together', `${rate(column.heating.total[k])} K per day`],
    ], every);
  }

  function pick({ x: px, y: py }) {
    if (!geometry || px > geometry.right || py < geometry.top - 10 || py > geometry.bottom + 8) return false;
    let best = -1, bestD = Infinity;
    for (let k = 0; k < K; k++) {
      if (view === 'low' && pair.column.height[k] > LOW_TOP) continue;
      const d = Math.abs(geometry.yOf(pair.column, k) - py);
      if (d < bestD) { bestD = d; best = k; }
    }
    if (best >= 0) { touched = true; if (best !== probe) { probe = best; fig.render(); } }
    return best >= 0;
  }
  const at = (e) => { const r = fig.canvas.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
  const ring = ({ y }) => !!geometry && Math.abs(y - geometry.yOf(pair.column, probe)) < 16;
  fig.pointer({ hit: ring, down: pick, move: pick });
  fig.canvas.addEventListener('click', (e) => pick(at(e)));
  fig.canvas.addEventListener('pointermove', (e) => { if (e.pointerType === 'mouse' && !e.buttons) pick(at(e)); });

  fig.play(true);
}
