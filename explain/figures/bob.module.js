import { Figure, slider, buttons, readout, legend, text, rampColor, clamp, ACCENT, MUTED, LINE, GRID, G, CP } from '../runtime.module.js';
import { thetaAt, pressureAt, T0, P0, KAPPA } from '../physics.module.js';

export function mountBob(root) {
  const controls = root.querySelector('.controls');
  const zTop = 10000, zRest = 3000, SPEED = 300, STEP = 2, WINDOW = 40 * 60;
  const lapseFor = (slope) => { let lapse = G / CP; for (let i = 0; i < 6; i++) lapse = G / CP - slope / 1000 / (P0 / pressureAt(zRest, T0, lapse)) ** KAPPA; return lapse; };
  const LID_JUMP = 8, LID_DEPTH = 300;
  let lapse = lapseFor(3.6), delta = 0, lid = 0, z = zRest, wv = 0, holding = false, t = 0, geometry = null;
  const history = [];
  const smooth = (x) => x <= 0 ? 0 : x >= 1 ? 1 : x * x * (3 - 2 * x);
  const thetaEnv = (zz) => thetaAt(zz, T0, lapse) + (lid > 0 ? LID_JUMP * smooth((zz - lid * 1000 + LID_DEPTH / 2) / LID_DEPTH) : 0);
  let thetaParcel = thetaEnv(zRest);
  function neutralHeight() {
    let previous = thetaEnv(0) - thetaParcel;
    for (let zz = 50; zz <= zTop; zz += 50) { const now = thetaEnv(zz) - thetaParcel; if (previous === 0 || previous * now < 0) return zz - 50 * now / (now - previous); previous = now; }
    return null;
  }
  let span = { mid: 300, width: 12 };
  const shade = (theta) => rampColor(0.5 + (theta - span.mid) / span.width);
  function rescale() {
    let lo = Infinity, hi = -Infinity;
    for (let zz = 0; zz <= zTop; zz += 250) { const th = thetaEnv(zz); lo = Math.min(lo, th); hi = Math.max(hi, th); }
    span = { mid: (lo + hi) / 2, width: Math.max(12, hi - lo) };
  }

  const fig = new Figure(root, { height: 400, step, draw });
  legend(root, [['ramp', 'lower to higher θ', 'cool', 'warm', 'neutral'], ['dash', 'where the parcel\u2019s θ matches the air', 'rgba(255,255,255,0.55)'], ['line', 'where the parcel has been', ACCENT]]);
  slider(controls, { label: 'θ of the air changes with height by', min: -3, max: 8, step: 0.1, value: 3.6, format: (v) => `${v > 0 ? '+' : v < 0 ? '−' : ''}${Math.abs(v).toFixed(1)} °C per km`, onInput: (v) => { lapse = lapseFor(v); thetaParcel = thetaEnv(zRest) + delta; rescale(); stability(); } });
  slider(controls, { label: 'Parcel θ compared with the air at 3 km', min: -8, max: 8, step: 0.5, value: 0, format: (v) => `${v > 0 ? '+' : v < 0 ? '−' : ''}${Math.abs(v).toFixed(1)} °C`, onInput: (v) => { delta = v; thetaParcel = thetaEnv(zRest) + delta; } });
  slider(controls, { label: 'Inversion lid at', min: 0, max: 8, step: 0.1, value: 0, format: (v) => v === 0 ? 'none' : `${v.toFixed(1)} km`, onInput: (v) => { lid = v; thetaParcel = thetaEnv(zRest) + delta; rescale(); stability(); } });
  buttons(controls, [['Nudge it up', () => { z = clamp(z + 1500, 0, zTop); wv = 0; }], ['Nudge it down', () => { z = clamp(z - 1500, 0, zTop); wv = 0; }], ['Reset', () => { z = zRest; wv = 0; history.length = 0; t = 0; }]]);
  const out = readout(controls);

  function stability() {
    const dz = 100, n2 = (G / thetaEnv(zRest)) * (thetaEnv(zRest + dz) - thetaEnv(zRest - dz)) / (2 * dz);
    const cools = ['the air cools by', `${(lapse * 1000).toFixed(1)} °C per km`];
    if (n2 > 2e-7) out.set([cools, ['stable', `it bobs with a period of ${(2 * Math.PI / Math.sqrt(n2) / 60).toFixed(1)} min`]]);
    else if (n2 > -2e-7) out.set([cools, ['neutral', 'it stays wherever it is put']]);
    else out.set([cools, ['unstable', 'it runs away: convection']]);
  }

  function step(dt) {
    for (let remaining = dt * SPEED; remaining > 0;) {
      const h = Math.min(STEP, remaining);
      remaining -= h;
      if (holding) continue;
      const te = thetaEnv(z);
      wv += (G * (thetaParcel - te) / te - wv / 1500) * h;
      z += wv * h;
      if (z <= 0) { z = 0; wv = 0; } else if (z >= zTop) { z = zTop; wv = 0; }
    }
    t += dt * SPEED;
    history.push(t, z);
    while (history.length && history[0] < t - WINDOW) history.splice(0, 2);
  }

  function draw(ctx, w, h) {
    const top = 24, bottom = h - 30, left = 46, colW = Math.min(110, w * 0.2), right = left + colW;
    const y = (zz) => bottom - (zz / zTop) * (bottom - top);
    geometry = { y, left, right };
    const bands = 40;
    for (let s = 0; s < bands; s++) { const z1 = (s / bands) * zTop, z2 = ((s + 1) / bands) * zTop; ctx.fillStyle = shade(thetaEnv((z1 + z2) / 2)); ctx.fillRect(left, y(z2), colW, y(z1) - y(z2) + 0.6); }
    ctx.fillStyle = '#5a4a36'; ctx.fillRect(left - 8, bottom, colW + 16, 6);
    for (let km = 0; km <= 10; km += 2) text(ctx, `${km} km`, left - 6, y(km * 1000), { align: 'right', color: MUTED, size: 11 });
    const level = neutralHeight();
    if (level !== null) { ctx.strokeStyle = 'rgba(255,255,255,0.55)'; ctx.setLineDash([3, 4]); ctx.beginPath(); ctx.moveTo(left, y(level)); ctx.lineTo(right, y(level)); ctx.stroke(); ctx.setLineDash([]); }
    text(ctx, 'θ of the surrounding air', left + colW / 2, 11, { align: 'center', color: MUTED, size: 11 });
    if (lid > 0) text(ctx, 'inversion', right + 6, y(lid * 1000), { color: MUTED, size: 10 });
    const cx = left + colW / 2, py = y(z);
    const cl = right + 44, cr = w - 16;
    ctx.strokeStyle = 'rgba(255,232,160,0.3)'; ctx.setLineDash([2, 4]); ctx.beginPath(); ctx.moveTo(cx + 13, py); ctx.lineTo(cr, py); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = shade(thetaParcel); ctx.strokeStyle = '#fff'; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.arc(cx, py, 13, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
    if (!holding && t === 0) text(ctx, 'drag me', cx, py - 22, { align: 'center', color: ACCENT, size: 11 });
    ctx.strokeStyle = GRID; for (let km = 0; km <= 10; km += 2) { ctx.beginPath(); ctx.moveTo(cl, y(km * 1000)); ctx.lineTo(cr, y(km * 1000)); ctx.stroke(); }
    ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(cl, bottom); ctx.lineTo(cr, bottom); ctx.stroke();
    text(ctx, 'height over the last 40 minutes', (cl + cr) / 2, bottom + 14, { align: 'center', color: MUTED, size: 11 });
    const x = (tt) => cr - ((t - tt) / WINDOW) * (cr - cl);
    ctx.strokeStyle = ACCENT; ctx.lineWidth = 2; ctx.beginPath();
    for (let i = 0; i < history.length; i += 2) { const xx = x(history[i]), yy = y(history[i + 1]); if (i === 0) ctx.moveTo(xx, yy); else ctx.lineTo(xx, yy); }
    ctx.stroke();
    ctx.fillStyle = ACCENT; ctx.beginPath(); ctx.arc(cr, py, 4, 0, Math.PI * 2); ctx.fill();
  }

  fig.pointer({
    hit: ({ x, y }) => geometry != null && Math.hypot(x - (geometry.left + geometry.right) / 2, y - geometry.y(z)) <= 24,
    down: () => { holding = true; },
    move: ({ y }) => { const { y: yOf } = geometry; z = clamp(zTop * (yOf(0) - y) / (yOf(0) - yOf(zTop)), 0, zTop); wv = 0; },
    up: () => { holding = false; },
  });
  rescale();
  stability();
  fig.play(true);
}
