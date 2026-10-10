import { Figure, buttons, legend, caption, readout, text, arrow, termColor, clamp, INK, MUTED, GRID, LINE } from '../runtime.module.js';
import { R, G, P0, heightOf, temperatureAt } from '../physics.module.js';

const K = 6, MAX = 20, STEP = 5, SPEED = 48, LIFE = 4, FADE = 0.6, DOTS = 300, GOLDEN = 0.6180339887, WIDE = 480, HALO = 'rgba(20, 20, 22, 0.85)';
const SIGMA = ['0', '1/6', '1/3', '1/2', '2/3', '5/6', '1'];
const NAMES = ['the top layer', 'the second from the top', 'the third from the top', 'the third from the ground', 'the second from the ground', 'the lowest layer'];
const COUNT = ['', 'one', 'two', 'three', 'four', 'five', 'six'];
const PRESETS = { rising: [0, -5, -10, -10, -10, -5, 0], sinking: [0, 5, 10, 10, 10, 5, 0], none: [0, 0, 0, 0, 0, 0, 0] };
const LESSON = ['The whole column’s total is always exactly zero.', 'The column’s total is always exactly zero.', 'Each arrow moves air from one layer into its neighbor, and no air can cross the top or the ground, so this term can never change the surface pressure.'];

export const changes = (flows) => Array.from({ length: K }, (_, k) => flows[k] - flows[k + 1]);
export function climb(i, tenths) {
  const p = i / K * P0, rho = p / (R * temperatureAt(heightOf(p)));
  return Math.abs(tenths) * 10 / (rho * G);
}
const fmt = (n) => (Math.abs(n) / 10).toFixed(1);
const signed = (n) => (n > 0 ? `+${fmt(n)}` : n < 0 ? `−${fmt(n)}` : '0');
const meters = (m) => (m < 9.95 ? m.toFixed(1) : `${Math.round(m)}`);

export function mountBetween(root) {
  const controls = root.querySelector('.controls'), C = termColor('c');
  const flows = [...PRESETS.rising];
  let selected = 3, grab = null, geometry = null, seed = 0;
  const even = () => (seed = (seed + GOLDEN) % 1);
  const dots = Array.from({ length: DOTS }, (_, n) => ({ x: Math.random(), s: even(), age: (n * 0.7548776662) % 1 * LIFE }));
  const moving = () => flows.some((f) => f !== 0);
  const speedAt = (s) => { const u = clamp(s * K, 0, K), i = Math.min(K - 1, Math.floor(u)), t = u - i; return (flows[i] * (1 - t) + flows[i + 1] * t) / 10; };

  const fig = new Figure(root, { height: 460, minHeight: 460, step, draw });
  const fit = fig.resize.bind(fig);
  fig.resize = () => { const narrow = fig.stage.clientWidth < WIDE; fig.height = fig.minHeight = narrow ? 500 : 460; fit(); };
  caption(root, `Drag an arrow up or down to set how much air crosses that level, in steps of half an hPa an hour, up to 2 either way. Each layer holds a sixth of the column’s ${P0 / 100} hPa of air, about ${Math.round(P0 / 100 / K)} hPa, and 1 hPa an hour carries about ${Math.round(100 / G)} kg of air across each square meter every hour. The dots are air carried along by the flow, with two days passing each second, and they crowd into the layers that gain.`);
  legend(root, [['varrow', 'air crossing between layers, πσ̇, rising or sinking', C], ['bar', 'what a layer gains through its top and bottom; hollow, what it loses', C]]);
  buttons(controls, [['Rising through the middle', () => preset('rising')], ['Sinking through the middle', () => preset('sinking')], ['No flow', () => preset('none')]]);
  const out = readout(controls);

  function preset(name) { PRESETS[name].forEach((v, i) => { flows[i] = v; }); changed(); }

  function changed() { report(); fig.play(moving()); fig.render(); }

  function who(ks) {
    const run = ks.every((k, n) => k === ks[0] + n);
    if (ks.length === 1) return NAMES[ks[0]];
    if (run && ks[0] === 0) return `the top ${COUNT[ks.length]} layers`;
    if (run && ks.at(-1) === K - 1) return `the lowest ${COUNT[ks.length]} layers`;
    if (ks.length === 2) return run && ks[0] === K / 2 - 1 ? 'the middle two layers' : `${NAMES[ks[0]]} and ${NAMES[ks[1]]}`;
    return `${COUNT[ks.length]} layers`;
  }

  function extreme(list, sign) {
    const best = Math.max(...list.map((c) => sign * c));
    return best <= 0 ? 'none' : `${fmt(best)} hPa an hour, ${who(list.flatMap((c, k) => (sign * c === best ? [k] : [])))}`;
  }

  function report() {
    const list = changes(flows), total = list.reduce((a, b) => a + b, 0), f = flows[selected];
    out.set([
      ['largest gain', extreme(list, 1)],
      ['largest loss', extreme(list, -1)],
      ['whole column', `${fmt(total)} hPa an hour, always: the surface pressure cannot change`],
      [`at σ = ${SIGMA[selected]} (${Math.round(selected / K * P0 / 100)} hPa)`, f === 0 ? 'no air crossing' : `${fmt(f)} hPa an hour, the air ${f < 0 ? 'rising' : 'sinking'} about ${meters(climb(selected, f))} m an hour`],
    ]);
  }

  function step(dt) {
    for (const d of dots) {
      d.s = clamp(d.s + speedAt(d.s) / (P0 / 100) * SPEED * dt, 0, 1);
      d.age += dt;
      if (d.age > LIFE) Object.assign(d, { x: Math.random(), s: even(), age: 0 });
    }
  }

  function wrap(ctx, str, width, font) {
    ctx.font = font;
    const lines = [''];
    for (const word of str.split(' ')) { const next = lines.at(-1) ? `${lines.at(-1)} ${word}` : word; if (lines.at(-1) && ctx.measureText(next).width > width) lines.push(word); else lines[lines.length - 1] = next; }
    return lines;
  }

  function bar(ctx, x, y, len, thick, gain) {
    if (len <= 0) return;
    ctx.fillStyle = C; ctx.strokeStyle = C; ctx.lineWidth = 1.2;
    if (gain) ctx.fillRect(x, y - thick / 2, len, thick);
    else { ctx.globalAlpha = 0.2; ctx.fillRect(x, y - thick / 2, len, thick); ctx.globalAlpha = 1; ctx.strokeRect(x + 0.6, y - thick / 2 + 0.6, Math.max(0, len - 1.2), thick - 1.2); }
  }

  function draw(ctx, w, h) {
    const narrow = w < WIDE, size = narrow ? 10 : 11, top = narrow ? 56 : 52, ground = h - (narrow ? 162 : 128), lh = (ground - top) / K;
    const x0 = narrow ? 46 : 71, colW = Math.round(Math.min(230, w * (narrow ? 0.32 : 0.34))), x1 = x0 + colW, ax = Math.round(x0 + (narrow ? 26 : colW * 0.28));
    const numW = narrow ? 28 : 76, p0 = x1 + (narrow ? 14 : 28), p1 = w - (narrow ? 8 : 12) - numW, mid = (p0 + p1) / 2, half = (p1 - p0) / 2, unit = half / (2 * MAX), reach = 0.225 * lh / 10;
    const yI = (i) => top + i * lh, thick = Math.min(16, lh * 0.34), list = changes(flows);
    const gains = list.reduce((a, c) => a + Math.max(0, c), 0), sum = Math.min(unit, (half - 3) / Math.max(1, gains));
    geometry = { x0, x1, top, lh, reach };

    for (let k = 0; k < K; k++) { ctx.fillStyle = k % 2 ? 'rgba(255,255,255,0.04)' : 'rgba(255,255,255,0.075)'; ctx.fillRect(x0, yI(k), colW, lh); }
    for (const d of dots) {
      const alpha = 0.5 * clamp(Math.min(d.age / FADE, (LIFE - d.age) / FADE), 0, 1);
      ctx.fillStyle = `rgba(255,255,255,${alpha.toFixed(3)})`; ctx.beginPath(); ctx.arc(x0 + 3 + d.x * (colW - 6), top + d.s * (ground - top), 1.6, 0, Math.PI * 2); ctx.fill();
    }
    ctx.fillStyle = '#5a4a36'; ctx.fillRect(x0 - 3, ground, colW + 6, 6);
    for (let i = 0; i <= K; i++) {
      const edge = i === 0 || i === K;
      ctx.strokeStyle = edge ? 'rgba(255,255,255,0.6)' : i === selected ? 'rgba(255,255,255,0.55)' : LINE; ctx.lineWidth = edge ? 1.5 : 1;
      ctx.beginPath(); ctx.moveTo(x0, yI(i)); ctx.lineTo(x1, yI(i)); ctx.stroke();
      text(ctx, `${edge || !narrow ? 'σ = ' : ''}${SIGMA[i]}`, x0 - 8, yI(i), { align: 'right', color: i === selected ? INK : MUTED, size });
    }
    ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(x0, top); ctx.lineTo(x0, ground); ctx.moveTo(x1, top); ctx.lineTo(x1, ground); ctx.stroke();
    text(ctx, narrow ? 'top' : 'the top', x0 - 8, top - 13, { align: 'right', color: MUTED, size });
    text(ctx, narrow ? 'ground' : 'the ground', x0 - 8, ground + 13, { align: 'right', color: MUTED, size });

    for (let i = 0; i <= K; i++) {
      const y = yI(i), f = flows[i], edge = i === 0 || i === K, label = edge ? 'always zero' : f === 0 ? 'no flow' : `${f < 0 ? 'rising' : 'sinking'} ${fmt(f)}`;
      if (f === 0) { ctx.strokeStyle = C; ctx.globalAlpha = edge ? 0.55 : 1; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(ax, y, 4, 0, Math.PI * 2); ctx.stroke(); ctx.globalAlpha = 1; }
      else arrow(ctx, ax, y - f * reach, ax, y + f * reach, { color: C, width: 2.5, head: 8 });
      text(ctx, label, ax + 11, y + (edge && i === K ? -9 : edge ? 9 : 0), { color: edge ? MUTED : INK, size, weight: i === selected ? 600 : 400, halo: HALO });
    }

    for (const [str, x, y] of [[narrow ? 'flow between layers' : 'air crossing between layers', x0, 12], [narrow ? 'hPa an hour' : 'hPa an hour · drag the arrows', x0, 25], [narrow ? 'each layer’s change' : 'what each layer gains or loses', p0, 12], ['hPa an hour', p0, 25]]) text(ctx, str, x, y, { color: MUTED, size });

    const yT = ground + (narrow ? 34 : 30);
    ctx.strokeStyle = GRID; ctx.lineWidth = 1; ctx.beginPath();
    for (let v = STEP * 2 - 2 * MAX; v <= 2 * MAX; v += STEP * 2) if (v) { ctx.moveTo(mid + v * unit, top); ctx.lineTo(mid + v * unit, ground); }
    ctx.stroke();
    for (const v of [-2 * MAX, -MAX, MAX, 2 * MAX]) text(ctx, `${v > 0 ? '+' : '−'}${Math.abs(v) / 10}`, mid + v * unit, ground + 10, { align: v === 2 * MAX ? 'right' : v === -2 * MAX ? 'left' : 'center', color: MUTED, size: size - 1 });
    ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(mid, top); ctx.lineTo(mid, ground + 4); ctx.moveTo(mid, yT - thick); ctx.lineTo(mid, yT + thick); ctx.stroke();
    text(ctx, '0', mid, ground + 10, { align: 'center', color: MUTED, size: size - 1, halo: HALO });
    ctx.save(); ctx.beginPath(); ctx.rect(p0, 0, p1 - p0, h); ctx.clip();
    for (let k = 0; k < K; k++) {
      const yc = yI(k) + lh / 2, c = list[k];
      ctx.fillStyle = 'rgba(255,255,255,0.04)'; ctx.fillRect(p0, yc - thick / 2, p1 - p0, thick);
      if (c > 0) bar(ctx, mid, yc, c * unit, thick, true); else bar(ctx, mid + c * unit, yc, -c * unit, thick, false);
    }
    ctx.fillStyle = 'rgba(255,255,255,0.04)'; ctx.fillRect(p0, yT - thick / 2, p1 - p0, thick);
    let right = mid, left = mid;
    for (const c of list) { if (c > 0) { bar(ctx, right, yT, c * sum - 1, thick, true); right += c * sum; } else if (c < 0) { left += c * sum; bar(ctx, left + 1, yT, -c * sum - 1, thick, false); } }
    ctx.restore();
    if (gains) for (const [str, x, align] of [[`losses ${fmt(gains)}`, mid - 6, 'right'], [`gains ${fmt(gains)}`, mid + 6, 'left']]) text(ctx, str, x, yT + thick / 2 + 10, { align, color: MUTED, size });
    for (let k = 0; k < K; k++) {
      const c = list[k];
      text(ctx, narrow ? signed(c) : c > 0 ? `gains ${fmt(c)}` : c < 0 ? `loses ${fmt(c)}` : 'no change', p1 + 8, yI(k) + lh / 2, { color: c ? INK : MUTED, size });
    }
    const total = list.reduce((a, b) => a + b, 0);
    text(ctx, 'whole column', p0 - 8, yT, { align: 'right', color: INK, size: narrow ? 11 : 12, weight: 600 });
    text(ctx, narrow ? fmt(total) : `total ${fmt(total)}`, p1 + 8, yT, { color: C, size: narrow ? 13 : 14, weight: 700 });

    const lines = [...wrap(ctx, LESSON[narrow ? 1 : 0], w - 28, `600 ${narrow ? 12 : 13}px system-ui, sans-serif`).map((s) => [s, INK, narrow ? 12 : 13, 600]), ...wrap(ctx, LESSON[2], w - 28, `400 ${size + 1}px system-ui, sans-serif`).map((s) => [s, MUTED, size + 1, 400])];
    lines.forEach(([s, color, sz, weight], n) => text(ctx, s, w / 2, yT + 40 + n * (narrow ? 16 : 15), { align: 'center', color, size: sz, weight }));
  }

  const nearest = ({ x, y }) => {
    if (!geometry) return null;
    const { x0, x1, top, lh } = geometry, i = Math.round((y - top) / lh);
    return x >= x0 && x <= x1 && y >= top + lh / 2 && y <= top + (K - 0.5) * lh ? clamp(i, 1, K - 1) : null;
  };
  const move = ({ y }) => {
    if (!grab) return;
    const v = clamp(Math.round((grab.v + (y - grab.y) / geometry.reach) / STEP) * STEP, -MAX, MAX);
    if (v !== flows[grab.i]) { flows[grab.i] = v; changed(); }
  };
  fig.pointer({ hit: (p) => nearest(p) !== null, down: (p) => { const i = nearest(p); grab = { i, y: p.y, v: flows[i] }; selected = i; report(); fig.render(); }, move, up: () => { grab = null; } });
  fig.canvas.addEventListener('pointermove', (e) => { const r = fig.canvas.getBoundingClientRect(); fig.canvas.style.cursor = grab || nearest({ x: e.clientX - r.left, y: e.clientY - r.top }) !== null ? 'ns-resize' : ''; });

  report();
  fig.play(moving());
}
