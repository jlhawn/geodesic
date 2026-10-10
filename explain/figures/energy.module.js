import { Figure, legend, readout, text, MUTED, GRID, LINE, INK, anomalyColor } from '../runtime.module.js';
import { loadTracks } from '../frames.module.js';

export function mountEnergy(root) {
  const controls = root.querySelector('.controls');
  let energy = null, hover = null, geometry = null;
  const fig = new Figure(root, { height: 320, minHeight: 260, draw });
  legend(root, [['line', 'change in kinetic energy', 'warm'], ['line', 'change in potential and internal energy', 'cool'], ['line', 'change in the total', INK]]);
  const out = readout(controls);
  const series = () => {
    const e0 = energy[0];
    return [
      { color: anomalyColor(1), values: energy.map((e) => (e.kinetic - e0.kinetic) / 1000) },
      { color: anomalyColor(-1), values: energy.map((e) => (e.internal - e0.internal) / 1000) },
      { color: INK, values: energy.map((e) => (e.kinetic + e.internal - e0.kinetic - e0.internal) / 1000) },
    ];
  };

  function draw(ctx, w, h) {
    if (!energy) return;
    const left = 56, right = 16, top = 16, bottom = h - 30, days = energy.map((e) => e.day), s = series();
    const lo = Math.min(...s.flatMap((x) => x.values)), hi = Math.max(...s.flatMap((x) => x.values)), span = hi - lo || 1;
    const x = (d) => left + d / days[days.length - 1] * (w - left - right), y = (v) => bottom - (v - lo) / span * (bottom - top);
    geometry = { x, left, right, w, days };
    ctx.strokeStyle = GRID; ctx.lineWidth = 1;
    for (let v = Math.ceil(lo / 50) * 50; v <= hi; v += 50) { ctx.beginPath(); ctx.moveTo(left, y(v)); ctx.lineTo(w - right, y(v)); ctx.stroke(); text(ctx, `${v}`, left - 6, y(v), { align: 'right', color: MUTED, size: 10 }); }
    ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(left, y(0)); ctx.lineTo(w - right, y(0)); ctx.stroke();
    for (let d = 0; d <= days[days.length - 1]; d += 2) text(ctx, `day ${d}`, x(d), bottom + 14, { align: 'center', color: MUTED, size: 10 });
    text(ctx, 'kJ per m²', left - 6, top - 6, { align: 'right', color: MUTED, size: 10 });
    for (const line of s) { ctx.strokeStyle = line.color; ctx.lineWidth = 2; ctx.beginPath(); line.values.forEach((v, n) => { if (n) ctx.lineTo(x(days[n]), y(v)); else ctx.moveTo(x(days[n]), y(v)); }); ctx.stroke(); }
    const n = hover ?? days.length - 1;
    ctx.strokeStyle = 'rgba(255,255,255,0.3)'; ctx.beginPath(); ctx.moveTo(x(days[n]), top); ctx.lineTo(x(days[n]), bottom); ctx.stroke();
    for (const line of s) { ctx.fillStyle = line.color; ctx.beginPath(); ctx.arc(x(days[n]), y(line.values[n]), 3.5, 0, Math.PI * 2); ctx.fill(); }
    const e = energy[n], e0 = energy[0];
    out.set([['day', days[n].toFixed(1)], ['kinetic', `${s[0].values[n] >= 0 ? '+' : ''}${s[0].values[n].toFixed(1)} kJ/m²`], ['potential and internal', `${s[1].values[n].toFixed(1)} kJ/m²`], ['total', `${s[2].values[n].toFixed(1)} kJ/m²`], ['mass', `changed by ${Math.abs((e.mass - e0.mass) / e0.mass).toExponential(0)}`]]);
  }

  const pick = (ev) => {
    if (!geometry) return;
    const r = fig.canvas.getBoundingClientRect(), px = ev.clientX - r.left, { days } = geometry;
    let best = 0; days.forEach((d, n) => { if (Math.abs(geometry.x(d) - px) < Math.abs(geometry.x(days[best]) - px)) best = n; });
    hover = best; fig.render();
  };
  fig.canvas.style.touchAction = 'pan-y';
  fig.canvas.addEventListener('pointerdown', pick);
  fig.canvas.addEventListener('pointermove', pick);
  fig.canvas.addEventListener('pointerleave', (ev) => { if (ev.pointerType === 'mouse') { hover = null; fig.render(); } });

  loadTracks(new URL('../data/tracks_N16.bin', import.meta.url).href).then(({ header }) => { energy = header.energy; fig.render(); });
}
