import { Figure, slider, readout, text, arrow, ACCENT, INK, MUTED, LINE, GRID } from '../runtime.module.js';
import { heightOf } from '../physics.module.js';

export function mountColumn(root) {
  const controls = root.querySelector('.controls');
  const levels = [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100], zMax = 20000;
  let surfaceT = 288.15, selected = 2, geometry = null;
  const zOf = (hPa) => heightOf(hPa * 100, surfaceT);

  const fig = new Figure(root, { height: 440, draw });
  const slab = slider(controls, { label: 'Slab', min: 0, max: levels.length - 1, step: 1, value: selected, format: (v) => `${levels[v]} to ${levels[v + 1] ?? 0} hPa`, onInput: (v) => { selected = v; fig.render(); } });
  slider(controls, { label: 'Surface temperature', min: -30, max: 40, step: 1, value: 15, format: (v) => `${v} °C`, onInput: (v) => { surfaceT = v + 273.15; fig.render(); show(); } });
  const out = readout(controls);
  const show = () => out.set([['half of the air is below', `${(zOf(500) / 1000).toFixed(1)} km`], ['nine tenths below', `${(zOf(100) / 1000).toFixed(1)} km`]]);

  function draw(ctx, w, h) {
    const narrow = w < 480, size = narrow ? 10 : 11, top = 22, bottom = h - 30, left = 70, colW = Math.min(190, w * (narrow ? 0.26 : 0.34)), right = left + colW;
    const y = (z) => bottom - (z / zMax) * (bottom - top);
    const bounds = (k) => ({ yb: y(zOf(levels[k])), yt: k + 1 < levels.length ? y(zOf(levels[k + 1])) : top });
    geometry = { left, right, bounds };
    for (let k = 0; k < levels.length; k++) {
      const { yb, yt } = bounds(k);
      let fill = k === selected ? 'rgba(255,232,160,0.25)' : k % 2 ? 'rgba(255,255,255,0.05)' : 'rgba(255,255,255,0.1)';
      if (k + 1 === levels.length) { const g = ctx.createLinearGradient(0, yb, 0, yt); g.addColorStop(0, fill); g.addColorStop(1, 'rgba(255,255,255,0)'); fill = g; }
      ctx.fillStyle = fill;
      ctx.fillRect(left, yt, colW, yb - yt);
      ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(left, yb); ctx.lineTo(right, yb); ctx.stroke();
      text(ctx, `${levels[k]} hPa`, left - 8, yb, { align: 'right', color: k === selected ? INK : MUTED, size });
    }
    ctx.strokeStyle = LINE; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(left, bottom); ctx.lineTo(left, top); ctx.moveTo(right, bottom); ctx.lineTo(right, top); ctx.stroke();
    ctx.fillStyle = '#5a4a36'; ctx.fillRect(left - 10, bottom, colW + 20, 6);
    for (let km = 0; km <= 20; km += 5) { ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(right, y(km * 1000)); ctx.lineTo(right + 6, y(km * 1000)); ctx.stroke(); text(ctx, `${km} km`, right + 10, y(km * 1000), { color: MUTED, size }); }
    const { yb, yt } = bounds(selected), below = levels[selected], above = levels[selected + 1] ?? 0, scale = 0.1;
    const px = narrow ? Math.max(right + 50, w - 130) : Math.min(w - 110, right + 110);
    const tail = Math.min(yb + below * scale, h - 4), ly = Math.min((yb + tail) / 2, h - 22);
    arrow(ctx, px, tail, px, yb, { color: ACCENT, width: 2.5, head: 8 });
    text(ctx, 'air below pushes up', px + 10, ly - 7, { color: INK, size });
    text(ctx, `${below} hPa`, px + 10, ly + 7, { color: ACCENT, size });
    if (above) {
      arrow(ctx, px, yt - above * scale, px, yt, { color: ACCENT, width: 2.5, head: 8 });
      text(ctx, 'air above pushes down', px + 10, yt - above * scale * 0.5 - 7, { color: INK, size });
      text(ctx, `${above} hPa`, px + 10, yt - above * scale * 0.5 + 7, { color: ACCENT, size });
    } else text(ctx, 'nothing above', px + 10, yt + 10, { color: MUTED, size });
    const mid = (yb + yt) / 2, weight = (below - above) * scale, wy = Math.min(mid, ly - 28);
    arrow(ctx, left + colW / 2, mid - weight / 2, left + colW / 2, mid + weight / 2, { color: '#fff', width: 2.5, head: 8 });
    text(ctx, 'weight of the slab', px + 10, wy - 7, { color: '#fff', size, weight: 500 });
    text(ctx, `${below - above} hPa`, px + 10, wy + 7, { color: '#fff', size, weight: 500 });
    text(ctx, narrow ? 'a tenth of the air per slab' : 'each slab holds one tenth of the air · click a slab', left, top - 10, { color: MUTED, size });
  }

  fig.canvas.addEventListener('click', (e) => {
    if (!geometry) return;
    const r = fig.canvas.getBoundingClientRect(), x = e.clientX - r.left, py = e.clientY - r.top;
    if (x < geometry.left - 60 || x > geometry.right) return;
    for (let k = 0; k < levels.length; k++) { const { yb, yt } = geometry.bounds(k); if (py <= yb && py >= yt) { selected = k; slab.value = k; fig.render(); break; } }
  });
  show();
}
