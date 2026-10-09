import { Figure, slider, choice, buttons, readout, legend, text, arrow, anomalyColor, clamp, ACCENT, INK, MUTED, LINE, GRID, P0 } from '../runtime.module.js';
import { createSlice } from '../sliceCore.module.js';

export function mountCells(root) {
  const controls = root.querySelector('.controls');
  const slice = createSlice({ dt: 10 });
  const { M, K, heated: c, dSigma } = slice;
  const z0 = Float64Array.from({ length: K + 1 }, (_, k) => slice.interfaceHeight(k, 0));
  const T0k = Float64Array.from({ length: K }, (_, k) => slice.temperature(k, 0));
  const EXAG = 20, SPEED = 600, zMax = 18000;
  let heating = 1, mode = 'column', friction = false;

  const fig = new Figure(root, { height: 470, minHeight: 360, step, draw });
  legend(root, [
    ['gradient', 'cooler or warmer than at the start', 'cool', 'warm'],
    ['arrow', 'wind between columns', 'rgba(255,255,255,0.85)'],
    ['varrow', 'rising or sinking air', ACCENT],
    ['bar', 'change in surface pressure', 'cool'],
  ]);
  controls.classList.add('two');
  slider(controls, { label: 'Warm the middle column by', min: 0, max: 3, step: 0.1, value: 1, format: (v) => `${v.toFixed(1)} °C per hour`, onInput: (v) => { heating = v; }, span: true });
  choice(controls, { label: 'Heat goes into', options: [['the whole column', 'column'], ['the ground', 'ground']], value: 'column', onChange: (v) => { mode = v; }, span: true });
  choice(controls, { label: 'Friction', options: [['off', 'off'], ['on', 'on']], value: 'off', onChange: (v) => { friction = v === 'on'; } });
  buttons(controls, [['Restart', () => { slice.reset(); show(); fig.render(); }], ['Pause', (b) => { fig.play(!fig.running); b.textContent = fig.running ? 'Pause' : 'Run'; }]]);
  const out = readout(controls);

  function step(dt) {
    let remaining = dt * SPEED;
    for (let n = 0; remaining > 0 && n < 200; n++) { const h = Math.min(slice.dt, remaining); slice.step({ heating: heating / 3600, mode, friction }, h); remaining -= h; }
    show(30);
  }

  function show(every = 1) {
    let maxU = 0, maxW = 0;
    for (let k = 0; k < K; k++) for (let i = 0; i < M; i++) { maxU = Math.max(maxU, Math.abs(slice.u[k * M + i])); maxW = Math.max(maxW, Math.abs(slice.verticalVelocity(k, i))); }
    const minutes = Math.round(slice.state.time / 60);
    out.set([['time', `${Math.floor(minutes / 60)} h ${String(minutes % 60).padStart(2, '0')} min`], ['strongest wind', `${maxU.toFixed(1)} m/s`], ['strongest updraught', `${(maxW * 100).toFixed(0)} cm/s`]], every);
  }

  function draw(ctx, w, h) {
    const left = 44, right = 56, top = 18, barsH = 90, bottom = h - barsH - 24;
    const plotW = w - left - right, colW = plotW / 5;
    const xc = (i) => left + (i - (c - 2) + 0.5) * colW;
    const yz = (z) => bottom - (z / zMax) * (bottom - top);
    const zi = (k, i) => z0[k] + EXAG * (slice.interfaceHeight(k, i) - z0[k]);
    const yEdge = (k, i, side) => yz(0.5 * (zi(k, i) + zi(k, i + side)));
    const yTop = (k, i, side) => k === 0 ? top - 5 : side === 0 ? yz(zi(k, i)) : yEdge(k, i, side);
    ctx.save();
    ctx.beginPath(); ctx.rect(left, top, plotW, bottom - top); ctx.clip();
    for (let i = c - 3; i <= c + 3; i++) {
      const xl = xc(i) - colW / 2, x = xc(i), xr = xc(i) + colW / 2;
      for (let k = 0; k < K; k++) {
        ctx.fillStyle = anomalyColor((slice.temperature(k, i) - T0k[k]) / 2.5, 0.6);
        ctx.beginPath();
        ctx.moveTo(xl, yTop(k, i, -1)); ctx.lineTo(x, yTop(k, i, 0)); ctx.lineTo(xr, yTop(k, i, 1));
        ctx.lineTo(xr, yEdge(k + 1, i, 1)); ctx.lineTo(x, yz(zi(k + 1, i))); ctx.lineTo(xl, yEdge(k + 1, i, -1));
        ctx.closePath(); ctx.fill();
      }
    }
    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    for (let k = 1; k <= K; k++) { ctx.beginPath(); for (let i = c - 3; i <= c + 3; i++) { const x = xc(i), y = yz(zi(k, i)); if (i === c - 3) ctx.moveTo(x, y); else ctx.lineTo(x, y); } ctx.stroke(); }
    ctx.strokeStyle = GRID;
    for (let i = c - 2; i <= c + 3; i++) { const x = xc(i) - colW / 2; ctx.beginPath(); ctx.moveTo(x, top); ctx.lineTo(x, bottom); ctx.stroke(); }
    const uScale = colW * 0.5 / 3;
    for (let e = c - 3; e <= c + 2; e++) {
      const x = xc(e) + colW / 2;
      for (let k = 1; k < K; k++) {
        const u = slice.u[k * M + e];
        if (Math.abs(u) < 0.05) continue;
        const y = 0.5 * (yEdge(k, e, 1) + yEdge(k + 1, e, 1)), len = clamp(u * uScale, -colW * 0.7, colW * 0.7);
        arrow(ctx, x - len / 2, y, x + len / 2, y, { color: 'rgba(255,255,255,0.85)', width: 1.5, head: 5 });
      }
    }
    for (let i = c - 3; i <= c + 3; i++) for (let k = 2; k < K; k++) {
      const wv = slice.verticalVelocity(k, i);
      if (Math.abs(wv) < 0.004) continue;
      const x = xc(i), y = yz(zi(k, i)), len = clamp(wv * 300, -36, 36);
      arrow(ctx, x, y + len / 2, x, y - len / 2, { color: ACCENT, width: 1.5, head: 5 });
    }
    ctx.restore();
    for (let km = 0; km <= 15; km += 5) text(ctx, `${km} km`, left - 6, yz(km * 1000), { align: 'right', color: MUTED, size: 11 });
    for (let k = 1; k < K; k++) if (z0[k] < zMax) text(ctx, `${Math.round(k * dSigma * 1000)} hPa`, left + plotW + 6, yz(z0[k]), { color: MUTED, size: 10 });
    ctx.fillStyle = '#5a4a36'; ctx.fillRect(left, bottom, plotW, 5);
    const note = 'heights exaggerated ×20', label = mode === 'ground' ? 'warmed ground' : 'warmed column';
    ctx.font = '11px system-ui, sans-serif';
    const collide = heating > 0 && left + 6 + ctx.measureText(note).width > xc(c) - ctx.measureText(label).width / 2;
    if (heating > 0) { ctx.fillStyle = ACCENT; ctx.fillRect(xc(c) - colW * 0.3, bottom, colW * 0.6, 5); text(ctx, label, xc(c), top + 10, { align: 'center', color: ACCENT, size: 11 }); }
    text(ctx, note, left + 6, collide ? top + 26 : top + 10, { color: MUTED, size: 11 });
    const base = bottom + 5 + (h - bottom - 5) / 2, perHPa = 16;
    const narrow = colW < 60;
    ctx.save(); ctx.translate(left - 30, base); ctx.rotate(-Math.PI / 2); text(ctx, narrow ? 'pressure, hPa' : 'surface pressure', 0, 0, { align: 'center', color: MUTED, size: 10 }); ctx.restore();
    ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(left, base); ctx.lineTo(left + plotW, base); ctx.stroke();
    for (let i = c - 2; i <= c + 2; i++) {
      const d = (slice.pi[i] - P0) / 100, len = clamp(d * perHPa, -36, 36);
      ctx.fillStyle = anomalyColor(Math.sign(d) * 0.8, 0.8);
      ctx.fillRect(xc(i) - colW * 0.25, Math.min(base, base - len), colW * 0.5, Math.abs(len));
      text(ctx, `${d >= 0 ? '+' : '−'}${Math.abs(d).toFixed(1)}${narrow ? '' : ' hPa'}`, xc(i), d >= 0 ? base - len - 9 : base - len + 9, { align: 'center', color: INK, size: 10 });
    }
  }

  show();
  fig.play(true);
}
