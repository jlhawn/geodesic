import { Figure, slider, choice, buttons, legend, caption, readout, text, termColor, clamp, INK, MUTED, LINE, G, P0 } from '../runtime.module.js';
import { heightOf } from '../physics.module.js';

const LEVELS = [0, 0.2, 0.4, 0.6, 0.75, 0.9, 1], K = LEVELS.length - 1, LOW = 960, HIGH = 1040, BASE = 1000, PUSH = 10, RATE = 5, BIG = 5, REACH = 18;
const tonnes = (pi, share) => pi * 100 * share / G / 1000;
const depth = (s) => heightOf(s * P0);

export function mountWeight(root) {
  const controls = root.querySelector('.controls'), A = termColor('a');
  let pi = BASE, target = null, mag = BIG, magShown = BIG, grab = null, geometry = null;
  const drawn = (p) => BASE + magShown * (p - BASE);

  const fig = new Figure(root, { height: 420, minHeight: 400, step, draw });
  caption(root, `Each layer is drawn as thick as the air it holds, not as tall as it is: the top layer, everything above about ${Math.round(depth(LEVELS[1]) / 1000)} km, holds a fifth of the air, and the lowest, only about ${Math.round(depth(LEVELS[K - 1]) / 100) * 100} m deep, holds a tenth. A layer between two levels of σ holds πΔσ/g of air over each square meter.`);
  const key = legend(root, []);
  const control = slider(controls, { label: 'Surface pressure π', min: LOW, max: HIGH, step: 1, value: pi, format: (v) => `${Math.round(v)} hPa`, onInput: (v) => { target = null; set(v); } });
  controls.lastElementChild.querySelector('output').style.color = A;
  buttons(controls, [['Pile air in', () => go(PUSH)], ['Let air out', () => go(-PUSH)]]);
  choice(controls, { label: 'Draw the change', options: [[`${BIG} times larger`, 'big'], ['at true size', 'true']], value: 'big', onChange: (v) => { mag = v === 'big' ? BIG : 1; label(); fig.play(true); } });
  const out = readout(controls), box = controls.querySelector('.readout');
  const label = () => key.set([['bar', mag === BIG ? `all the air in the column, which π measures, with changes from 1000 hPa drawn ${BIG} times larger` : 'all the air in the column, which π measures', A], ['dash', 'where the top of the column sits at 1000 hPa', MUTED]]);

  function show() {
    const p = Math.round(pi), d = p - BASE;
    out.set([['surface pressure', `${p} hPa`], ['air in the whole column', `${tonnes(p, 1).toFixed(2)} metric tons over each square meter`], ['the top layer’s share', `${Math.round(100 * (LEVELS[1] - LEVELS[0]))} percent, at any pressure`], ['compared with 1000 hPa', d ? `${d > 0 ? '+' : '−'}${Math.abs(d)} hPa, ${Math.round(Math.abs(d) * 100 / G)} kg ${d > 0 ? 'more' : 'less'} air over each square meter` : 'the same']]);
    for (const b of [...box.querySelectorAll('b')].slice(0, 2)) b.style.color = A;
  }

  function set(v) { pi = v; control.value = v; show(); fig.render(); }
  function go(d) { target = clamp((target ?? Math.round(pi)) + d, LOW, HIGH); fig.play(true); }

  function step(dt) {
    let changed = false;
    if (target !== null) {
      const gap = target - pi, move = RATE * dt;
      pi = Math.abs(gap) <= move ? target : pi + Math.sign(gap) * move;
      if (pi === target) target = null;
      control.value = pi; show(); changed = true;
    }
    if (magShown !== mag) { magShown = Math.abs(mag - magShown) < 0.01 ? mag : magShown + (mag - magShown) * (1 - Math.exp(-dt * 8)); changed = true; }
    if (target === null && magShown === mag) fig.play(false);
    return changed;
  }

  function draw(ctx, w, h) {
    const narrow = w < 480, size = narrow ? 10 : 11, labels = narrow ? 108 : 120, colW = narrow ? 104 : 170, after = narrow ? 118 : 150, p = Math.round(pi);
    const left = Math.max(0, (w - labels - colW - after) / 2) + labels, right = left + colW, mid = (left + right) / 2;
    const bottom = h - 34, scale = (bottom - 44) / (BASE + BIG * (HIGH - BASE)), y = (s) => bottom - (1 - s) * drawn(pi) * scale, top = y(0);
    geometry = { left, right, bottom, scale, top };
    text(ctx, 'drag the top of the column to add or remove air', w / 2, 14, { align: 'center', color: MUTED, size: 11 });
    for (let k = 0; k < K; k++) {
      const yt = y(LEVELS[k]), yb = y(LEVELS[k + 1]), share = LEVELS[k + 1] - LEVELS[k];
      ctx.fillStyle = k % 2 ? 'rgba(255,255,255,0.05)' : 'rgba(255,255,255,0.1)';
      ctx.fillRect(left, yt, colW, yb - yt);
      text(ctx, `${tonnes(p, share).toFixed(2)} t/m²`, mid, (yt + yb) / 2, { align: 'center', color: INK, size });
      text(ctx, `${Math.round(share * 100)}%`, right + 10, (yt + yb) / 2, { color: MUTED, size });
    }
    if (Math.abs(pi - BASE) >= 0.5) { const yr = bottom - BASE * scale; ctx.strokeStyle = MUTED; ctx.lineWidth = 1; ctx.setLineDash([3, 4]); ctx.beginPath(); ctx.moveTo(left - 8, yr); ctx.lineTo(right + 8, yr); ctx.stroke(); ctx.setLineDash([]); }
    ctx.strokeStyle = LINE; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(left, bottom); ctx.lineTo(left, top); ctx.moveTo(right, bottom); ctx.lineTo(right, top); ctx.stroke();
    LEVELS.forEach((s, k) => {
      const yy = y(s), ground = k === K;
      if (k) { ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(left - 4, yy); ctx.lineTo(right + 4, yy); ctx.stroke(); }
      text(ctx, `${Math.round(s * p)} hPa`, left - 10, yy, { align: 'right', color: ground ? A : INK, size, weight: ground ? 600 : 400 });
      text(ctx, `σ = ${s}`, left - (narrow ? 60 : 66), yy, { align: 'right', color: MUTED, size });
    });
    ctx.strokeStyle = A; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(left, top); ctx.lineTo(right, top); ctx.stroke();
    ctx.fillStyle = A; ctx.beginPath(); ctx.roundRect(mid - (grab === null ? 14 : 18), top - 4, grab === null ? 28 : 36, 8, 4); ctx.fill();
    ctx.fillStyle = '#5a4a36'; ctx.fillRect(left - 6, bottom, colW + 12, 6);
    text(ctx, 'the ground', mid, bottom + 20, { align: 'center', color: MUTED, size });
    const bx = right + (narrow ? 36 : 52), by = (top + bottom) / 2;
    ctx.strokeStyle = A; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(bx - 5, top); ctx.lineTo(bx, top); ctx.lineTo(bx, bottom); ctx.lineTo(bx - 5, bottom); ctx.stroke();
    text(ctx, 'whole column', bx + 8, by - 17, { color: MUTED, size });
    text(ctx, `${tonnes(p, 1).toFixed(2)} t/m²`, bx + 8, by, { color: A, size: size + 2, weight: 600 });
    text(ctx, '= π / g', bx + 8, by + 17, { color: MUTED, size });
  }

  const near = (q) => geometry !== null && Math.abs(q.y - geometry.top) <= REACH && q.x >= geometry.left - 12 && q.x <= geometry.right + 12;
  fig.pointer({
    hit: near,
    down: (q) => { target = null; grab = geometry.top - q.y; fig.render(); },
    move: (q) => { const { bottom, scale } = geometry; set(clamp(Math.round(BASE + ((bottom - q.y - grab) / scale - BASE) / magShown), LOW, HIGH)); },
    up: () => { grab = null; fig.render(); },
  });
  fig.canvas.addEventListener('pointermove', (e) => { if (grab !== null) return; const r = fig.canvas.getBoundingClientRect(); fig.canvas.style.cursor = near({ x: e.clientX - r.left, y: e.clientY - r.top }) ? 'ns-resize' : ''; });
  label();
  show();
}
