import { Figure, slider, choice, readout, text, anomalyColor, thermometer, ACCENT, MUTED, LINE, KAPPA, R } from '../runtime.module.js';

export function mountPiston(root) {
  const controls = root.querySelector('.controls');
  const T0 = 293.15;
  let T = T0, p = 1000, mode = 'insulated';
  const dots = Array.from({ length: 90 }, () => { const a = Math.random() * Math.PI * 2, s = 0.25 + Math.random() * 0.4; return { x: Math.random(), y: Math.random(), vx: Math.cos(a) * s, vy: Math.sin(a) * s }; });
  const volume = () => (T / T0) / (p / 1000);

  const fig = new Figure(root, { height: 380, step, draw });
  slider(controls, { label: 'Weight on the piston', min: 500, max: 1500, step: 10, value: 1000, format: (v) => `${v} hPa`, onInput: (v) => { T *= (v / p) ** KAPPA; p = v; show(); } });
  choice(controls, { label: 'The cylinder', options: [['is insulated', 'insulated'], ['leaks heat', 'leaky']], value: 'insulated', onChange: (v) => { mode = v; } });
  const out = readout(controls);

  function step(dt) {
    const s = Math.sqrt(T / T0);
    for (const d of dots) {
      d.x += d.vx * s * dt; d.y += d.vy * s * dt;
      if (d.x < 0) { d.x = -d.x; d.vx = -d.vx; } else if (d.x > 1) { d.x = 2 - d.x; d.vx = -d.vx; }
      if (d.y < 0) { d.y = -d.y; d.vy = -d.vy; } else if (d.y > 1) { d.y = 2 - d.y; d.vy = -d.vy; }
    }
    if (mode === 'leaky') T += (T0 - T) * (1 - Math.exp(-dt / 1.2));
    show(12);
  }

  function show(every = 1) {
    out.set([['pressure', `${p} hPa`], ['temperature', `${(T - 273.15).toFixed(1)} °C`], ['volume', `×${volume().toFixed(2)}`], ['density', `${(p * 100 / (R * T)).toFixed(2)} kg/m³`]], every);
  }

  function draw(ctx, w, h) {
    const cw = Math.min(220, w * 0.42), cx = w * 0.42, bottom = h - 28, gh0 = (h - 80) / 2.15, gh = gh0 * volume();
    const left = cx - cw / 2, right = cx + cw / 2, pistonY = bottom - gh;
    ctx.fillStyle = anomalyColor((T - T0) / 18, 0.45);
    ctx.fillRect(left, pistonY, cw, gh);
    ctx.fillStyle = '#fff';
    for (const d of dots) { ctx.beginPath(); ctx.arc(left + 4 + d.x * (cw - 8), pistonY + 4 + d.y * (gh - 8), 2.2, 0, Math.PI * 2); ctx.fill(); }
    ctx.strokeStyle = LINE; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(left, 24); ctx.lineTo(left, bottom); ctx.lineTo(right, bottom); ctx.lineTo(right, 24); ctx.stroke();
    ctx.fillStyle = '#b8b8b8';
    ctx.fillRect(left - 1, pistonY - 10, cw + 2, 10);
    const block = 16 + (p - 500) * 0.07, bw = cw * 0.45;
    ctx.fillStyle = '#8c8c8c';
    ctx.fillRect(cx - bw / 2, pistonY - 10 - block, bw, block);
    text(ctx, `${p} hPa`, cx, pistonY - 10 - block / 2, { align: 'center', color: '#151515', size: 12, weight: 600 });
    if (mode === 'insulated') { ctx.strokeStyle = ACCENT; ctx.setLineDash([3, 4]); ctx.lineWidth = 1; ctx.strokeRect(left - 6, 18, cw + 12, bottom - 18 + 6); ctx.setLineDash([]); }
    text(ctx, mode === 'insulated' ? 'insulated' : 'heat leaks to the room', cx, bottom + 15, { align: 'center', color: MUTED, size: 12 });
    thermometer(ctx, w * 0.8, 40, h - 120, T - 273.15, { min: -20, max: 80, label: 'temperature' });
  }

  fig.play(true);
}
