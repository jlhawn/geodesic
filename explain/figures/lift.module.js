import { Figure, slider, readout, text, thermometer, rampColor, ACCENT, MUTED, GRID, KAPPA, P0 } from '../runtime.module.js';
import { pressureAt, heightOf, T0 } from '../physics.module.js';

export function mountLift(root) {
  const controls = root.querySelector('.controls');
  const zMax = 12000;
  let z = 0;
  const parcel = () => { const p = pressureAt(z), T = T0 * (p / P0) ** KAPPA; return { p, T, volume: (T / T0) / (p / P0) }; };

  const fig = new Figure(root, { height: 420, draw });
  slider(controls, { label: 'Lift the parcel to', min: 0, max: zMax, step: 50, value: 0, format: (v) => `${(v / 1000).toFixed(2)} km`, onInput: (v) => { z = v; fig.render(); show(); } });
  const out = readout(controls);
  const show = () => { const { p, T, volume } = parcel(); out.set([['pressure there', `${(p / 100).toFixed(0)} hPa`], ['parcel temperature', `${(T - 273.15).toFixed(1)} °C`], ['volume', `×${volume.toFixed(2)}`]]); };

  function draw(ctx, w, h) {
    const narrow = w < 480, top = 30, bottom = h - 34, left = 44, panelW = Math.min(240, w * (narrow ? 0.34 : 0.45)), right = left + panelW;
    const y = (zz) => bottom - (zz / zMax) * (bottom - top);
    const sky = ctx.createLinearGradient(0, top, 0, bottom);
    sky.addColorStop(0, 'rgba(20, 30, 70, 0.9)'); sky.addColorStop(1, 'rgba(70, 120, 190, 0.6)');
    ctx.fillStyle = sky; ctx.fillRect(left, top, panelW, bottom - top);
    ctx.fillStyle = '#5a4a36'; ctx.fillRect(left, bottom, panelW, 6);
    for (let km = 0; km <= 12; km += 2) { ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(left, y(km * 1000)); ctx.lineTo(right, y(km * 1000)); ctx.stroke(); text(ctx, `${km} km`, left - 6, y(km * 1000), { align: 'right', color: MUTED, size: 11 }); }
    for (const hPa of narrow ? [1000, 700, 500, 300, 200] : [1000, 850, 700, 500, 400, 300, 200]) { const zz = heightOf(hPa * 100); if (zz > zMax) continue; text(ctx, `${hPa} hPa`, right + 5, y(zz), { color: MUTED, size: 10 }); }
    const { T, volume } = parcel(), r = 13 * Math.sqrt(volume), cx = left + panelW / 2, cy = y(z) - r * (1 - z / zMax);
    ctx.strokeStyle = 'rgba(255,232,160,0.5)'; ctx.setLineDash([2, 5]); ctx.beginPath(); ctx.moveTo(cx, bottom); ctx.lineTo(cx, cy + r); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = rampColor((T - 213) / 90); ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.arc(cx, cy, r, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
    const free = w - right - 50, tx1 = right + 50 + free * (narrow ? 0.18 : 0.25), tx2 = right + 50 + free * (narrow ? 0.66 : 0.75);
    thermometer(ctx, tx1, 46, h - 130, T - 273.15, { min: -80, max: 40, label: narrow ? ['temperature', 'now'] : 'temperature now' });
    thermometer(ctx, tx2, 46, h - 130, T0 - 273.15, { min: -80, max: 40, label: narrow ? ['brought', 'back down'] : 'brought back down' });
    text(ctx, narrow ? '≈ −9.8 °C/km' : 'about −9.8 °C per km', tx1, h - 14, { align: 'center', color: MUTED, size: 11 });
    text(ctx, narrow ? 'θ' : 'potential temperature θ', tx2, h - 14, { align: 'center', color: ACCENT, size: 11 });
  }

  show();
}
