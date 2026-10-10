import { Figure, choice, legend, readout, text, rampRGB, sequentialRGB, clamp, MUTED, GRID, INK } from '../runtime.module.js';

const P_TOP = 50;

function contour(ctx, grid, levels, x, y, style) {
  const { lats, pressures, values } = grid, nx = lats.length, ny = pressures.length;
  for (const level of levels) {
    ctx.save(); ctx.strokeStyle = style(level).color; ctx.lineWidth = style(level).width; ctx.setLineDash(style(level).dash ?? []);
    ctx.beginPath();
    for (let i = 0; i + 1 < nx; i++) for (let j = 0; j + 1 < ny; j++) {
      const corners = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]].map(([a, b]) => ({ a, b, v: values[a][b] }));
      const points = [];
      for (let k = 0; k < 4; k++) {
        const p = corners[k], q = corners[(k + 1) % 4];
        if ((p.v - level) * (q.v - level) < 0 || (p.v === level && q.v !== level)) {
          const t = (level - p.v) / (q.v - p.v), la = lats[p.a] + t * (lats[q.a] - lats[p.a]), pr = pressures[p.b] + t * (pressures[q.b] - pressures[p.b]);
          points.push([x(la), y(pr)]);
        }
      }
      if (points.length >= 2) { ctx.moveTo(...points[0]); ctx.lineTo(...points[1]); }
      if (points.length === 4) { ctx.moveTo(...points[2]); ctx.lineTo(...points[3]); }
    }
    ctx.stroke(); ctx.restore();
  }
}

export function drawCrossSection(ctx, w, h, { lats, pressures, temperature, theta, wind = null, tRange = [200, 310] }) {
  const left = 52, right = 14, top = 14, bottom = h - 34;
  const x = (lat) => left + (lat + 90) / 180 * (w - left - right), y = (p) => top + (p - P_TOP) / (1000 - P_TOP) * (bottom - top);
  const rgb = [0, 0, 0];
  for (let i = 0; i < lats.length; i++) for (let j = 0; j < pressures.length; j++) {
    const x0 = x(i ? (lats[i - 1] + lats[i]) / 2 : -90), x1 = x(i + 1 < lats.length ? (lats[i] + lats[i + 1]) / 2 : 90);
    const p0 = j ? (pressures[j - 1] + pressures[j]) / 2 : P_TOP, p1 = j + 1 < pressures.length ? (pressures[j] + pressures[j + 1]) / 2 : 1000;
    if (p1 < P_TOP) continue;
    sequentialRGB((temperature[i][j] - tRange[0]) / (tRange[1] - tRange[0]), rgb);
    ctx.fillStyle = `rgb(${rgb.map((c) => Math.round(c * 255)).join(',')})`;
    ctx.fillRect(x0, y(Math.max(p0, P_TOP)), x1 - x0 + 0.5, y(p1) - y(Math.max(p0, P_TOP)) + 0.5);
  }
  ctx.save(); ctx.beginPath(); ctx.rect(left, top, w - left - right, bottom - top); ctx.clip();
  const thetaLevels = []; for (let v = 260; v <= 600; v += 10) thetaLevels.push(v);
  contour(ctx, { lats, pressures, values: theta }, thetaLevels, x, y, (v) => ({ color: v % 50 === 0 ? 'rgba(255,255,255,0.7)' : 'rgba(255,255,255,0.3)', width: 1 }));
  if (wind) {
    const windLevels = []; for (let v = -30; v <= 50; v += 5) if (v) windLevels.push(v);
    contour(ctx, { lats, pressures, values: wind }, windLevels, x, y, (v) => ({ color: v > 0 ? 'rgba(120, 220, 255, 0.95)' : 'rgba(120, 220, 255, 0.95)', width: v % 10 === 0 ? 2 : 1.2, dash: v < 0 ? [4, 3] : [] }));
  }
  ctx.restore();
  ctx.strokeStyle = GRID;
  for (const p of [200, 400, 600, 800, 1000]) { ctx.beginPath(); ctx.moveTo(left, y(p)); ctx.lineTo(w - right, y(p)); ctx.stroke(); text(ctx, `${p} hPa`, left - 6, y(p), { align: 'right', color: MUTED, size: 10 }); }
  for (const lat of [-90, -60, -30, 0, 30, 60, 90]) text(ctx, lat === 0 ? '0°' : `${Math.abs(lat)}°${lat > 0 ? 'N' : 'S'}`, x(lat), bottom + 14, { align: 'center', color: MUTED, size: 10 });
  void INK;
}

export function mountCrossSection(root, { modes, initial = 0, legendFor }) {
  const controls = root.querySelector('.controls');
  let mode = initial;
  const fig = new Figure(root, { height: 360, minHeight: 280, draw: (ctx, w, h) => { const d = modes[mode].data(); if (d) drawCrossSection(ctx, w, h, d); } });
  const key = legend(root, legendFor(mode));
  if (modes.length > 1) choice(controls, { label: 'Show', options: modes.map((m, k) => [m.label, String(k)]), value: String(mode), onChange: (v) => { mode = Number(v); key.set(legendFor(mode)); show(); fig.render(); }, span: true });
  const out = readout(controls);
  function show() { const r = modes[mode].readout?.(); if (r) out.set(r); }
  show();
  return { refresh: () => { show(); fig.render(); } };
}
