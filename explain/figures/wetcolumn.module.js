// Shadows the drawing of equilibrium.module.js; a fix there belongs here too.
import { Figure, slider, choice, buttons, readout, legend, caption, text, arrow, termColor, clamp, INK, MUTED, GRID, LINE, ACCENT } from '../runtime.module.js';
import { K } from '../radiativeColumn.module.js';
import { tropopause, logScale, pickLayers, daysPerSecond, STEPS_PER_DAY, PPM } from '../radiativeConvective.module.js';
import { createWetColumn, convectionTop, troposphericLapse, parcelHumidity, moistAdiabat, nearestLayer } from '../moistConvective.module.js';
import { CP_DRY, R_DRY, P0 } from '../../js/dynamics/sigmaCore.module.js';

const KAPPA = R_DRY / CP_DRY, TOP_PRESSURE = 10, LOW_TOP = 20000, WIDE = 520;
const SMALL = 0.03, LARGE = 40, BUDGET = 10, MOST = 75, GHOST = 'rgba(255,255,255,0.4)', HALO = 'rgba(20, 20, 22, 0.85)', COLD = 'rgba(255,255,255,0.7)';
const PRESSURES = [1000, 500, 200, 100, 50, 20, 10, 5, 2, 1, 0.5, 0.2, 0.1];
const TERMS = [['sunlight', 'a'], ['infrared', 'b'], ['ground', 'c'], ['convection', 'd']];
const NAMES = { sunlight: 'sunlight', infrared: 'infrared', ground: 'sensible heat from the ground', convection: 'convection and condensation' };
const signed = (v, d = 2) => (Math.abs(v) < 0.5 * 10 ** -d ? (0).toFixed(d) : `${v < 0 ? '−' : '+'}${Math.abs(v).toFixed(d)}`);
const rate = (v) => signed(v, Math.abs(v) >= 9.995 ? 1 : 2);
const plain = (v, d = 1) => (Math.abs(v) < 0.5 * 10 ** -d ? (0).toFixed(d) : `${v < 0 ? '−' : ''}${Math.abs(v).toFixed(d)}`);
const hPa = (p) => { const v = p / 100; return v >= 9.95 ? v.toFixed(0) : v >= 0.995 ? v.toFixed(1) : v.toFixed(2); };
const meters = (z) => (z < 995 ? `${Math.round(z / 10) * 10} m` : `${(z / 1000).toFixed(1)} km`);

export function mountWetcolumn(root) {
  const controls = root.querySelector('.controls'), color = Object.fromEntries([...'abcde'].map((c) => [c, termColor(c)]));
  const run = createWetColumn({ albedo: 0.3 }), column = run.column;
  const up = new Float64Array(K), guide = new Float64Array(K);
  let view = 'low', running = true, probe = nearestLayer(column, 500e2), touched = false, owed = 0, range = null, goal = null, geometry = null, compared = false;

  const fig = new Figure(root, { height: (w) => (w < WIDE ? 500 : 480), step, draw });
  caption(root, 'Part 6’s column of 36 layers under the global average of sunlight, now holding water vapor, over ground that reflects 30 percent of the sunlight, stepped every 3 simulated hours from 288 K everywhere. Moist convection, written for this figure and not taken from the climate model, lifts the lowest layer’s air along its moist adiabat and releases the latent heat of the water the ground evaporates; its rate is each step’s adjustment over the 3 hours. With convection off the column steps every hour instead: with this much vapor its lowest layer, only 40 m thick, cools by its own infrared within about an hour, and 3-hour steps would overshoot back and forth without bound. “Frozen” holds each layer’s vapor at its amount when pressed. On the left is the temperature of each layer and, on the brown strip, of the ground; on the right, on Part 6’s tenfold scale, is how fast each process warms or cools a dozen of the layers. Tap or drag to read a layer. The run pauses itself once the column stops changing.');
  const key = legend(root, legendItems());
  controls.classList.add('two');
  choice(controls, { label: 'Show', options: [['whole column', 'whole'], ['lowest 20 km', 'low']], value: view, onChange: (v) => { view = v; if (view === 'low' && column.height[probe] > LOW_TOP) probe = nearestLayer(column, 500e2); settle(); } });
  const [, toggle] = buttons(controls, [['Restart', () => { run.restart(); resume(); }], ['Pause', () => { if (running) { halt(); settle(); } else { if (run.settled()) run.calm = 0; resume(); } }]]);
  choice(controls, { label: 'Water vapor', options: [['none', 'none'], ['frozen', 'frozen'], ['fixed relative humidity', 'relative']], value: run.vapor, onChange: (v) => change({ vapor: v }), span: true });
  choice(controls, { label: 'Convection', options: [['off', 'off'], ['dry', 'dry'], ['moist', 'moist']], value: run.convection, onChange: (v) => change({ convection: v }), span: true });
  slider(controls, { label: 'Carbon dioxide', min: 0, max: 4, step: 0.25, value: 1, format: (v) => `×${v}, ${Math.round(v * PPM)} ppm`, onInput: (v) => change({ carbonDioxide: v }), span: true });
  slider(controls, { label: 'Ground reflects', min: 0.07, max: 0.4, step: 0.01, value: 0.3, format: (v) => `${Math.round(v * 100)}% of sunlight`, onInput: (v) => change({ albedo: v }), span: true });
  const out = readout(controls), layerOut = readout(controls);
  for (const box of controls.querySelectorAll('.readout')) box.style.gridColumn = '1 / -1';

  function legendItems() {
    return [
      ['line', 'temperature', INK], ...(compared ? [['line', 'last settled state', GHOST]] : []), ['dash', 'dry adiabat', LINE], ['dash', 'moist adiabat', color.d],
      ['faint', 'top of convection', color.d], ['dots', 'cold point', COLD],
      ['arrow', 'sunlight', color.a], ['arrow', 'infrared', color.b], ['arrow', 'sensible heat from the ground', color.c], ['arrow', 'convection and condensation', color.d], ['arrow', 'all together', color.e],
    ];
  }

  function change(changes) {
    run.set(changes);
    if (run.before && !compared) { compared = true; key.set(legendItems()); }
    resume();
  }
  function resume() { running = true; toggle.textContent = 'Pause'; settle(); }
  function halt() { running = false; toggle.textContent = 'Run'; }
  function settle() { goal = target(); fig.play(running || !range || !same(range, goal)); fig.render(); }
  const same = (a, b) => Math.abs(a[0] - b[0]) < 0.05 && Math.abs(a[1] - b[1]) < 0.05;
  const shown = (c, k) => view === 'whole' || c.height[k] <= LOW_TOP * 1.05;

  function target() {
    let lo = Infinity, hi = -Infinity;
    for (const c of run.before ? [column, run.before] : [column]) {
      for (let k = 0; k < K; k++) if (shown(c, k)) { lo = Math.min(lo, c.T[k]); hi = Math.max(hi, c.T[k]); }
      lo = Math.min(lo, c.Ts); hi = Math.max(hi, c.Ts);
    }
    lo = Math.floor((lo - 6) / 10) * 10; hi = Math.ceil((hi + 6) / 10) * 10;
    if (hi - lo < 60) { const mid = (hi + lo) / 2; lo = Math.floor((mid - 30) / 10) * 10; hi = lo + 60; }
    return [lo, hi];
  }

  function step(dt) {
    if (running) {
      owed = Math.min(owed + dt * daysPerSecond(run.since) * STEPS_PER_DAY, MOST);
      const start = performance.now();
      for (; owed >= 1 && !run.settled() && performance.now() - start < BUDGET; owed--) run.step();
      owed = Math.min(owed, 1);
      if (run.settled()) halt();
    }
    goal = target();
    range = !range ? goal : same(range, goal) ? goal : range.map((v, n) => v + (goal[n] - v) * Math.min(1, dt * 6));
    if (!running && same(range, goal)) { range = goal; fig.play(false); }
  }

  function layout(w, h) {
    const narrow = w < WIDE, left = narrow ? 34 : 44, right = w - (narrow ? 6 : 10), top = 32, bottom = h - 46, gap = narrow ? 14 : 26;
    const mid = left + (right - left - gap) * 0.56, r0 = mid + gap, cx = (r0 + right) / 2, half = (right - r0) / 2 - 2;
    const span = Math.log(column.interfacePressure[K] / TOP_PRESSURE);
    const yP = (p) => bottom - Math.log(column.interfacePressure[K] / p) / span * (bottom - top);
    const yZ = (z) => bottom - z / LOW_TOP * (bottom - top);
    const yOf = view === 'whole' ? (c, k) => yP(c.pressure[k]) : (c, k) => yZ(c.height[k]);
    return { narrow, left, right, top, bottom, mid, r0, cx, half, yOf, yP, yZ };
  }

  function curve(ctx, points) {
    ctx.beginPath();
    points.forEach(([px, py], n) => (n ? ctx.lineTo(px, py) : ctx.moveTo(px, py)));
    ctx.stroke();
  }

  function draw(ctx, w, h) {
    range ??= target();
    const g = geometry = layout(w, h), { narrow, left, top, bottom, mid, r0, right, cx, half, yOf, yP, yZ } = g, before = run.before;
    const [lo, hi] = range, x = (T) => left + (T - lo) / (hi - lo) * (mid - left);
    const ground = view === 'whole' ? yP(column.interfacePressure[K]) : yZ(0);
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
      for (let z = 0; z <= LOW_TOP; z += 2000) {
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
    const cold = tropopause(column), crest = convectionTop(column), marks = [];
    if (cold) {
      const y = yOf(column, cold.k);
      ctx.strokeStyle = COLD; ctx.lineWidth = 2; ctx.lineCap = 'round'; ctx.setLineDash([0.5, 4]);
      ctx.beginPath(); ctx.moveTo(left, y); ctx.lineTo(mid, y); ctx.stroke(); ctx.setLineDash([]); ctx.lineCap = 'butt';
      marks.push([y, 'cold point', COLD]);
    }
    if (crest && crest.k < K - 1) {
      const y = yOf(column, crest.k);
      ctx.strokeStyle = color.d; ctx.lineWidth = 1.2;
      ctx.beginPath(); ctx.moveTo(left, y); ctx.lineTo(mid, y); ctx.stroke();
      marks.push([y, 'top of convection', color.d]);
    }

    const ps = column.interfacePressure[K], pb = column.pressure[K - 1], Tb = column.T[K - 1], theta = Tb * (P0 / pb) ** KAPPA;
    for (let n = 0; n < K; n++) up[n] = column.pressure[K - 1 - n];
    const along = (T, n) => [x(T), yOf(column, K - 1 - n)];
    ctx.strokeStyle = LINE; ctx.lineWidth = 1.2; ctx.setLineDash([4, 4]);
    curve(ctx, [[x(theta * (ps / P0) ** KAPPA), ground], ...Array.from({ length: K }, (_, n) => along(theta * (up[n] / P0) ** KAPPA, n))]);
    if (parcelHumidity(column) > 0) {
      moistAdiabat(column, up, guide);
      ctx.strokeStyle = color.d; ctx.lineWidth = 1.4; ctx.setLineDash([6, 4]);
      curve(ctx, [[x(Tb * (ps / pb) ** KAPPA), ground], ...Array.from({ length: K }, (_, n) => along(guide[n], n))]);
    }
    ctx.setLineDash([]);
    if (before) profile(ctx, before, x, yOf, GHOST, 1.5, false);
    profile(ctx, column, x, yOf, INK, 2, true);
    marks.sort((a, b) => a[0] - b[0]);
    marks.forEach(([y, label, c], n) => {
      const below = marks.length === 2 && n === 1 && marks[1][0] - marks[0][0] < 15;
      if (y > top - 2 && y < bottom - 4) text(ctx, label, mid - 4, below ? y + 8 : y - 7, { align: 'right', color: c, size: 10, halo: HALO });
    });
    ctx.restore();

    ctx.fillStyle = '#5a4a36'; ctx.fillRect(left, bottom, mid - left, 5);
    ctx.lineWidth = 1.5;
    if (before) { ctx.strokeStyle = GHOST; ctx.fillStyle = HALO; ctx.beginPath(); ctx.arc(x(before.Ts), bottom + 2.5, 4, 0, Math.PI * 2); ctx.fill(); ctx.stroke(); }
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
    const centers = Array.from({ length: K }, (_, k) => (shown(column, k) ? place(k) : NaN));
    const picked = pickLayers(centers, narrow ? 20 : 22, [probe, K - 1]);
    const py = centers[probe], layerY = yOf(column, probe);
    ctx.strokeStyle = 'rgba(255,255,255,0.3)'; ctx.lineWidth = 1; ctx.setLineDash([2, 4]); ctx.beginPath(); ctx.moveTo(left, layerY); ctx.lineTo(mid, layerY); ctx.lineTo(r0, py); ctx.lineTo(right, py); ctx.stroke(); ctx.setLineDash([]);
    for (const k of picked) {
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
      if (!touched) { const after = rx < (left + mid) / 2; text(ctx, 'drag', after ? rx + 10 : rx - 10, ry - 11, { align: after ? 'left' : 'right', color: ACCENT, size: 11, halo: HALO }); }
    }
    text(ctx, `day ${Math.floor(column.time / 86400).toLocaleString('en-US')}`, mid, top - 14, { align: 'right', color: MUTED, size: 10 });
    report(fig.running ? 4 : 1, cold, crest);
  }

  function profile(ctx, c, x, yOf, stroke, width, dots) {
    ctx.strokeStyle = stroke; ctx.lineWidth = width; ctx.lineJoin = 'round'; ctx.beginPath();
    for (let k = 0; k < K; k++) { const px = x(c.T[k]), py = yOf(c, k); if (k) ctx.lineTo(px, py); else ctx.moveTo(px, py); }
    ctx.stroke();
    if (!dots) return;
    ctx.fillStyle = stroke;
    for (let k = 0; k < K; k++) { ctx.beginPath(); ctx.arc(x(c.T[k]), yOf(c, k), 1.8, 0, Math.PI * 2); ctx.fill(); }
  }

  function report(every, cold, crest) {
    const b = column.budget, day = Math.floor(column.time / 86400).toLocaleString('en-US'), wet = run.convection === 'moist' && run.vapor !== 'none';
    out.set([
      ['day', `${day}${!running && run.settled() ? ', settled' : ''}`],
      ['ground', `${column.Ts.toFixed(1)} K`],
      ...(run.before ? [['since the last settled state', `${signed(column.Ts - run.before.Ts)} K`]] : []),
      ['sunlight absorbed', `${b.absorbed.toFixed(1)} W/m²`],
      ['reflected', `${(b.albedo * 100).toFixed(1)}%`],
      ['infrared leaving the top', `${b.olr.toFixed(1)} W/m²`],
      ['held back, G', `${plain(b.greenhouse)} W/m²`],
      ['precipitable water', `${b.precipitableWater.toFixed(1)} kg/m²`],
      ...(wet ? [['evaporation and rain', `${(b.evaporation * 86400).toFixed(2)} mm per day`], ['latent heat, L·E', `${b.latent.toFixed(1)} W/m²`]] : []),
      ...(run.convection !== 'off' ? [['sensible heat', `${plain(b.sensible)} W/m²`]] : []),
      ['top of convection', crest ? `${hPa(crest.pressure)} hPa, ${meters(crest.height)}` : 'none'],
      ['cold point', cold ? `${hPa(cold.pressure)} hPa, ${meters(cold.height)}, ${cold.temperature.toFixed(1)} K` : 'none yet'],
      ['lapse rate from 850 to 300 hPa', `${troposphericLapse(column).toFixed(1)} K per km`],
    ], every);
    const k = probe, vapor = column.q[k] > 0 ? `, ${(column.q[k] * 1000).toFixed(column.q[k] >= 1e-3 ? 2 : 4)} g/kg of vapor` : '';
    layerOut.set([
      ['layer', `${K - k} of ${K}, at ${hPa(column.pressure[k])} hPa and ${meters(column.height[k])}, ${column.T[k].toFixed(1)} K${vapor}`],
      ...TERMS.filter(([name]) => column.heating[name][k] !== 0).map(([name]) => [NAMES[name], rate(column.heating[name][k])]),
      ['all together', `${rate(column.heating.total[k])} K per day`],
    ], every);
  }

  function pick({ x: px, y: py }) {
    if (!geometry || px > geometry.right || py < geometry.top - 10 || py > geometry.bottom + 8) return false;
    let best = -1, bestD = Infinity;
    for (let k = 0; k < K; k++) {
      if (!shown(column, k)) continue;
      const d = Math.abs(geometry.yOf(column, k) - py);
      if (d < bestD) { bestD = d; best = k; }
    }
    if (best >= 0) { touched = true; if (best !== probe) { probe = best; fig.render(); } }
    return best >= 0;
  }
  const at = (e) => { const r = fig.canvas.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
  const ring = ({ y }) => !!geometry && Math.abs(y - geometry.yOf(column, probe)) < 16;
  fig.pointer({ hit: ring, down: pick, move: pick });
  fig.canvas.addEventListener('click', (e) => pick(at(e)));
  fig.canvas.addEventListener('pointermove', (e) => { if (e.pointerType === 'mouse' && !e.buttons) pick(at(e)); });

  fig.play(true);
}
