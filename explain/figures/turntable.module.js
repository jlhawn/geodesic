import { Figure, slider, buttons, legend, text, ACCENT, INK, MUTED, LINE } from '../runtime.module.js';

export function mountTurntable(root) {
  const controls = root.querySelector('.controls');
  const CROSSING = 4;
  let spin = 0.12, t = 0, launched = true, angle = 0;
  const fixed = [], riding = [];

  const fig = new Figure(root, { height: 360, step, draw });
  legend(root, [['line', 'the puck’s path', ACCENT]]);
  slider(controls, { label: 'Turntable spin', min: 0, max: 0.3, step: 0.01, value: spin, format: (v) => `${v.toFixed(2)} turns per second`, onInput: (v) => { spin = v; launch(); } });
  buttons(controls, [['Launch again', launch]]);

  function launch() { t = 0; angle = 0; fixed.length = 0; riding.length = 0; launched = true; }

  function step(dt) {
    if (!launched) return;
    t += dt;
    angle += 2 * Math.PI * spin * dt;
    const x = -1 + 2 * t / CROSSING, y = 0.35 * Math.sin(0.7);
    if (x > 1.05) { launched = false; return; }
    const r = Math.hypot(x, y);
    if (r <= 1) {
      fixed.push([x, y]);
      const c = Math.cos(-angle), s = Math.sin(-angle);
      riding.push([x * c - y * s, x * s + y * c]);
    }
  }

  function draw(ctx, w, h) {
    const R = Math.min(w / 4.6, h / 2.6), yc = h / 2 + 8;
    const panels = [[w * 0.27, angle, 'seen from above, standing still', fixed], [w * 0.73, 0, 'seen riding the turntable', riding]];
    for (const [xc, spokeAngle, label, trace] of panels) {
      ctx.fillStyle = 'rgba(255,255,255,0.06)'; ctx.beginPath(); ctx.arc(xc, yc, R, 0, Math.PI * 2); ctx.fill();
      ctx.strokeStyle = LINE; ctx.lineWidth = 1.5; ctx.stroke();
      for (let k = 0; k < 6; k++) { const a = spokeAngle + k * Math.PI / 3; ctx.strokeStyle = 'rgba(255,255,255,0.15)'; ctx.beginPath(); ctx.moveTo(xc, yc); ctx.lineTo(xc + R * Math.cos(a), yc - R * Math.sin(a)); ctx.stroke(); }
      ctx.fillStyle = 'rgba(255,255,255,0.5)'; ctx.beginPath(); ctx.arc(xc + R * 0.92 * Math.cos(spokeAngle), yc - R * 0.92 * Math.sin(spokeAngle), 4, 0, Math.PI * 2); ctx.fill();
      ctx.strokeStyle = ACCENT; ctx.lineWidth = 2; ctx.beginPath();
      trace.forEach(([x, y], i) => { const sx = xc + x * R, sy = yc - y * R; if (i) ctx.lineTo(sx, sy); else ctx.moveTo(sx, sy); });
      ctx.stroke();
      if (trace.length) { const [x, y] = trace[trace.length - 1]; ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(xc + x * R, yc - y * R, 6, 0, Math.PI * 2); ctx.fill(); }
      text(ctx, label, xc, yc + R + 18, { align: 'center', color: MUTED, size: 11 });
    }
    text(ctx, 'the same puck, sliding in a straight line at a steady speed', w / 2, 14, { align: 'center', color: MUTED, size: 11 });
    void INK;
  }

  launch();
  fig.play(true);
}
