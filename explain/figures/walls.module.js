import { Figure, choice, buttons, legend, caption, readout, text, arrow, termColor, clamp, G, INK, MUTED, GRID, LINE } from '../runtime.module.js';

const SPACING = 60e3, WALL = SPACING / Math.sqrt(3), APOTHEM = SPACING / 2, AREA = Math.sqrt(3) / 2 * SPACING ** 2, PER_HOUR = WALL / AREA * 360000;
const MAX = 15, SPEED = 1800, SUB = 120, DOTS = 130, LOOP = 8 * 3600, HOLD = 1.2, ENTER = 0.35, LEAVE = 0.5, WIDE = 560, PANEL = 230, BUDGET = 172, HALO = 'rgba(20, 20, 22, 0.85)';
const ANGLES = Array.from({ length: 6 }, (_, e) => e * Math.PI / 3), NORMALS = ANGLES.map((t) => [Math.cos(t), Math.sin(t)]), TANGENTS = ANGLES.map((t) => [-Math.sin(t), Math.cos(t)]);
const exact = (v) => Math.round(v * 1e6) / 1e6;
const PRESETS = {
  uniform: { winds: NORMALS.map(([nx]) => exact(10 * nx)), spin: 0 },
  converging: { winds: NORMALS.map(() => -1), spin: 0 },
  diverging: { winds: NORMALS.map(() => 1), spin: 0 },
  spinning: { winds: NORMALS.map(() => 0), spin: 10 },
  front: { winds: NORMALS.map(([nx]) => exact((nx > 0 ? 9 : 10) * nx)), spin: 0 },
};
const SUPERSCRIPT = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
const sci = (v) => { const [m, p] = Math.abs(v).toExponential(1).split('e'); return `${v < 0 ? '−' : ''}${m} × 10${[...String(Number(p))].map((c) => SUPERSCRIPT[c]).join('')}`; };
const tons = (pct) => Number((pct / 360000 * AREA * 100 / G / 1000).toPrecision(3)).toLocaleString('en-US');
const speed = (v) => `${Number.isInteger(v) ? v : v.toFixed(1)} m/s`;
const signed = (v) => (v > 0 ? `+${Math.round(v)}%` : v < 0 ? `−${Math.round(-v)}%` : '0');

function hex(x, y) {
  let r = -Infinity, e = 0;
  NORMALS.forEach(([nx, ny], i) => { const d = (x * nx + y * ny) / APOTHEM; if (d > r) { r = d; e = i; } });
  return [r, e];
}

export function mountWalls(root) {
  const controls = root.querySelector('.controls'), B = termColor('b');
  let winds = [...PRESETS.uniform.winds], spin = 0, dots = [], clock = 0, amount = 1, hold = 0, drag = null, geometry = null, seed = 7;
  let ux = 0, uy = 0, c = 0, s1 = 0, s2 = 0, q = 0, rate = 0, gained = 0, inflow = [], peak = [], budget = { into: 0, out: 0 };
  const acc = new Float64Array(6);
  const random = () => { seed = (seed + 0x6d2b79f5) | 0; let t = Math.imul(seed ^ (seed >>> 15), 1 | seed); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
  const velocity = (x, y) => [ux + (c + s1) * x + s2 * y + q * (x * x - y * y), uy + (c - s1) * y + s2 * x - 2 * q * x * y];
  const onWall = (e, tau) => [APOTHEM * NORMALS[e][0] + tau * WALL / 2 * TANGENTS[e][0], APOTHEM * NORMALS[e][1] + tau * WALL / 2 * TANGENTS[e][1]];
  const across = (e, tau) => { const [x, y] = onWall(e, tau), [u, v] = velocity(x, y); return u * NORMALS[e][0] + v * NORMALS[e][1]; };
  const moving = () => spin !== 0 || winds.some((u) => u !== 0);

  const fig = new Figure(root, { height: 380, minHeight: 300, step, draw });
  const fit = fig.resize.bind(fig);
  fig.resize = () => { const tall = fig.stage.clientWidth < WIDE; fig.height = tall ? 460 : 380; fig.minHeight = tall ? 460 : 300; fit(); };
  caption(root, `One of the model’s cells at N = 128, seen from above, with its six neighbors: their centers are ${SPACING / 1000} km apart, each wall is ${(WALL / 1000).toFixed(1)} km long, and the cell covers ${(Math.round(AREA / 1e7) * 10).toLocaleString('en-US')} km². π is the same in all seven cells, so the air crossing a wall is set by the wind across it alone. Drag the ring on any wall, in or out, to set the wind across it. Time runs an hour every ${3600 / SPEED} seconds.`);
  const key = legend(root, []);
  const preset = choice(controls, { label: 'Wind', options: [['a uniform wind', 'uniform'], ['converging', 'converging'], ['diverging', 'diverging'], ['spinning', 'spinning'], ['a front', 'front'], ['your own', 'own']], value: 'uniform', onChange: (v) => { if (v in PRESETS) { winds = [...PRESETS[v].winds]; spin = PRESETS[v].spin; configure(); reset(); fig.render(); } }, span: true });
  buttons(controls, [['Restart', () => { reset(); fig.play(true); fig.render(); }]]);
  const out = readout(controls);

  function configure() {
    const sum = (f) => winds.reduce((s, u, e) => s + u * f(e), 0);
    c = sum(() => 1) / 6 / APOTHEM; ux = sum((e) => NORMALS[e][0]) / 3; uy = sum((e) => NORMALS[e][1]) / 3;
    s1 = sum((e) => Math.cos(2 * ANGLES[e])) / 3 / APOTHEM; s2 = sum((e) => Math.sin(2 * ANGLES[e])) / 3 / APOTHEM; q = 9 / 8 * sum((e) => (e % 2 ? -1 : 1)) / 6 / APOTHEM ** 2;
    inflow = []; peak = [];
    for (let e = 0; e < 6; e++) {
      let total = 0, most = 0;
      for (let j = 0; j < 24; j++) { const v = Math.max(0, -across(e, -1 + (j + 0.5) / 12)); total += v / 24; most = Math.max(most, v); }
      inflow.push(total * WALL); peak.push(most);
    }
    rate = winds.reduce((s, u) => s + u, 0) * WALL / AREA;
    budget = { into: winds.reduce((s, u) => s + Math.max(0, -u), 0) * PER_HOUR, out: winds.reduce((s, u) => s + Math.max(0, u), 0) * PER_HOUR };
    gained = budget.into - budget.out;
    out.set([
      ['air in', `${Math.round(budget.into)}% of the cell’s air an hour`], ['out', `${Math.round(budget.out)}%`],
      ['in a layer 1 hPa thick', `${tons(budget.into)} metric tons a second in and ${tons(budget.out)} out, ${gained === 0 ? 'no difference' : `${tons(Math.abs(gained))} more ${gained > 0 ? 'in than out' : 'out than in'}`}`],
      ['divergence', gained === 0 ? '0: the layer keeps its air' : `${sci(-gained / 360000)} per second: the layer ${gained > 0 ? 'gains' : 'loses'} ${Math.round(Math.abs(gained))}% of its air an hour`],
      ...(spin ? [['along the walls', `${speed(spin)}, crossing none of them`]] : []),
    ]);
    key.set([['arrow', 'air crossing a wall, π u l: π times the wind across the wall times its length', B], ['dots', 'air that was in the cell at the start', INK], ['dots', 'air that came in through the walls', B], ...(spin ? [['arrow', 'wind along the walls', 'rgba(221, 221, 221, 0.5)']] : [])]);
    if (moving()) fig.play(true);
  }

  function reset() {
    dots = []; clock = 0; amount = 1; hold = 0;
    for (let e = 0; e < 6; e++) acc[e] = random();
    while (dots.length < DOTS) { const x = (random() * 2 - 1) * WALL, y = (random() * 2 - 1) * WALL; if (hex(x, y)[0] < 1) dots.push({ x, y, fresh: false, state: 0, age: ENTER }); }
    launch();
  }

  function launch() {
    dots = dots.filter((d) => d.state !== 1);
    for (let e = 0; e < 6; e++) for (let n = Math.floor(DOTS / AREA * amount * inflow[e] * ENTER * SPEED + random()); n > 0; n--) { const back = random() * ENTER; spawn(e, 0, back * SPEED).age = ENTER - back; }
  }

  function rotate(p, walls) {
    const [r, e] = hex(p.x, p.y);
    if (r < 1e-9) return;
    const tau = clamp((p.x * TANGENTS[e][0] + p.y * TANGENTS[e][1]) / (r * WALL / 2), -1, 1), s = e + (tau + 1) / 2 + walls, k = Math.floor(s), [x, y] = onWall(((k % 6) + 6) % 6, 2 * (s - k) - 1);
    p.x = r * x; p.y = r * y;
  }

  function move(p, s) {
    const [u1, v1] = velocity(p.x, p.y), [u2, v2] = velocity(p.x + u1 * s / 2, p.y + v1 * s / 2);
    p.x += u2 * s; p.y += v2 * s;
    if (spin) rotate(p, spin * s / WALL);
  }

  function travel(p, s) { for (let left = Math.abs(s); left > 1e-6; left -= SUB) move(p, Math.sign(s) * Math.min(SUB, left)); }

  function spawn(e, s, back = ENTER * SPEED + random() * s) {
    let tau = 0;
    for (let tries = 0; tries < 30; tries++) { tau = random() * 2 - 1; if (random() * peak[e] < -across(e, tau)) break; }
    const [x, y] = onWall(e, tau), dot = { x, y, fresh: true, state: 1, age: 0 };
    travel(dot, -back);
    dots.push(dot);
    return dot;
  }

  function advance(s) {
    clock += s; amount *= Math.exp(-rate * s);
    for (let e = 0; e < 6; e++) { acc[e] += DOTS / AREA * amount * inflow[e] * s; while (acc[e] >= 1) { acc[e] -= 1; spawn(e, s); } }
    for (const d of dots) {
      move(d, s);
      const inside = hex(d.x, d.y)[0] <= 1;
      if (d.state === 0 && !inside) { d.state = 2; d.age = 0; }
      else if (d.state === 1 && inside) d.state = 0;
    }
  }

  function step(dt) {
    const fading = dots.some((d) => d.state);
    if (!moving() && !fading) { fig.play(false); return false; }
    if (hold > 0) { hold -= dt; if (hold > 0) return false; reset(); return true; }
    if (moving()) for (let left = dt * SPEED; left > 1e-6; left -= SUB) advance(Math.min(SUB, left));
    for (const d of dots) d.age += dt;
    dots = dots.filter((d) => !(d.state === 2 && d.age > LEAVE) && !(d.state === 1 && (d.age > 4 * ENTER || !moving())));
    if (amount > 2 || amount < 0.5 || clock >= LOOP) hold = HOLD;
    return true;
  }

  function outline(ctx, x0, y0, r) {
    ctx.beginPath();
    for (let j = 0; j < 6; j++) { const t = Math.PI / 6 + j * Math.PI / 3, x = x0 + r * Math.cos(t), y = y0 - r * Math.sin(t); if (j) ctx.lineTo(x, y); else ctx.moveTo(x, y); }
    ctx.closePath();
  }

  function cell(ctx, vw, vh) {
    const ap = Math.max(30, Math.min((vw / 2 - 34) / 1.675, (vh / 2 - 22) / 1.45)), scale = ap / APOTHEM, k = 0.09 * ap, cx = vw / 2, cy = vh / 2 + 4, wallPx = WALL * scale;
    const sx = (x) => cx + x * scale, sy = (y) => cy - y * scale;
    geometry = { cx, cy, ap, k, wallPx };
    ctx.save(); ctx.beginPath(); ctx.rect(0, 0, vw, vh); ctx.clip();
    ctx.strokeStyle = GRID; ctx.lineWidth = 1;
    for (const [nx, ny] of NORMALS) { outline(ctx, sx(SPACING * nx), sy(SPACING * ny), wallPx); ctx.stroke(); }
    outline(ctx, cx, cy, wallPx); ctx.fillStyle = 'rgba(255, 255, 255, 0.025)'; ctx.fill(); ctx.strokeStyle = LINE; ctx.lineWidth = 1.5; ctx.stroke();
    const radius = clamp(ap * 0.024, 1.6, 2.8);
    for (const fresh of [false, true]) {
      ctx.fillStyle = fresh ? B : INK;
      for (const d of dots) {
        if (d.fresh !== fresh) continue;
        ctx.globalAlpha = 0.85 * (d.state === 2 ? clamp(1 - d.age / LEAVE, 0, 1) : clamp(d.age / ENTER, 0, 1));
        ctx.beginPath(); ctx.arc(sx(d.x), sy(d.y), radius, 0, Math.PI * 2); ctx.fill();
      }
    }
    ctx.globalAlpha = 1;
    ctx.restore();
    if (spin) for (let e = 0; e < 6; e++) {
      const [nx, ny] = NORMALS[e], [tx, ty] = TANGENTS[e], len = spin * k, mx = cx + nx * (ap - 12), my = cy - ny * (ap - 12);
      arrow(ctx, mx - tx * len / 2, my + ty * len / 2, mx + tx * len / 2, my - ty * len / 2, { color: 'rgba(221, 221, 221, 0.45)', width: 1.5, head: 6 });
    }
    ctx.font = '400 10px system-ui, sans-serif';
    winds.forEach((u, e) => {
      const [nx, ny] = NORMALS[e], mx = cx + nx * ap, my = cy - ny * ap, half = Math.abs(u) * k / 2, sign = Math.sign(u), held = drag?.e === e;
      ctx.beginPath(); ctx.arc(mx, my, held ? 6 : 4, 0, Math.PI * 2); ctx.fillStyle = HALO; ctx.fill();
      if (held) { ctx.globalAlpha = 0.3; ctx.fillStyle = INK; ctx.fill(); ctx.globalAlpha = 1; }
      ctx.strokeStyle = held ? INK : MUTED; ctx.lineWidth = 1.2; ctx.stroke();
      arrow(ctx, mx - sign * nx * half, my + sign * ny * half, mx + sign * nx * half, my - sign * ny * half, { color: B, width: 1.5 + 0.22 * Math.abs(u), head: 6 + 0.45 * Math.abs(u) });
      const label = u ? `${speed(Math.abs(u))} ${u > 0 ? 'out' : 'in'}` : '0 m/s', room = ctx.measureText(label).width / 2 + 3, side = Math.abs(ny) < 0.1;
      const off = side ? room + 8 : Math.max(half, 5) + 19, lx = clamp(mx + nx * off, room + 2, vw - room - 2), ly = side ? my - 15 : clamp(my - ny * off, 8, vh - 8);
      text(ctx, label, lx, ly, { align: 'center', color: held ? INK : MUTED, size: 10, halo: HALO });
    });
    text(ctx, 'north ↑', 10, 14, { color: MUTED, size: 11 });
  }

  function bar(ctx, x, y, width, f, color, alpha) {
    ctx.fillStyle = GRID; ctx.fillRect(x, y - 4.5, width, 9);
    ctx.globalAlpha = alpha; ctx.fillStyle = color; ctx.fillRect(x, y - 4.5, width * clamp(f, 0, 1), 9); ctx.globalAlpha = 1;
  }

  function panel(ctx, [x0, y0, x1, y1]) {
    const pad = 14, rows = 26, height = 18 + 3 * rows + 14 + 42, top = y0 + Math.max(8, (y1 - y0 - height) / 2), left = x0 + pad, right = x1 - pad, bx = left + 34, bw = right - 40 - bx;
    const full = Math.max(100, budget.into, budget.out);
    text(ctx, 'share of the cell’s air, per hour', left, top + 6, { color: MUTED, size: 11 });
    [['in', budget.into, `${Math.round(budget.into)}%`, 0.55], ['out', budget.out, `${Math.round(budget.out)}%`, 0.55], ['net', Math.abs(gained), signed(gained), 1]].forEach(([name, v, label, alpha], n) => {
      const y = top + 18 + rows * (n + 0.5);
      text(ctx, name, left, y, { color: n === 2 ? INK : MUTED, size: 11 });
      bar(ctx, bx, y, bw, v / full, B, alpha);
      text(ctx, label, right, y, { align: 'right', color: n === 2 ? INK : MUTED, size: 11 });
    });
    const gy = top + 18 + 3 * rows + 14;
    text(ctx, 'the cell’s air', left, gy + 6, { color: MUTED, size: 11 });
    text(ctx, `${Math.round(amount * 100)}% after ${(clock / 3600).toFixed(1)} h`, right, gy + 6, { align: 'right', color: INK, size: 11 });
    const gx = left, gw = right - left, by = gy + 24;
    bar(ctx, gx, by, gw, amount / 2, INK, 0.7);
    ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(gx + gw / 2, by - 8); ctx.lineTo(gx + gw / 2, by + 8); ctx.stroke();
    text(ctx, 'start', gx + gw / 2, by + 15, { align: 'center', color: MUTED, size: 10 });
    text(ctx, '0', gx, by + 15, { color: MUTED, size: 10 });
    text(ctx, 'twice', gx + gw, by + 15, { align: 'right', color: MUTED, size: 10 });
  }

  function draw(ctx, w, h) {
    const tall = w < WIDE, vw = tall ? w : w - PANEL, vh = tall ? h - BUDGET : h;
    cell(ctx, vw, vh);
    panel(ctx, tall ? [0, vh, w, h] : [vw, 0, w, h]);
  }

  const along = (e, p) => (p.x - geometry.cx - NORMALS[e][0] * geometry.ap) * NORMALS[e][0] - (p.y - geometry.cy + NORMALS[e][1] * geometry.ap) * NORMALS[e][1];
  function pick(p) {
    if (!geometry) return -1;
    let best = -1, near = 20;
    winds.forEach((u, e) => {
      const [nx, ny] = NORMALS[e], dx = p.x - geometry.cx - nx * geometry.ap, dy = p.y - geometry.cy + ny * geometry.ap, side = dx * ny + dy * nx, d = Math.hypot(side, Math.max(0, Math.abs(along(e, p)) - Math.abs(u) * geometry.k / 2));
      if (d < near) { near = d; best = e; }
    });
    return best;
  }
  fig.pointer({
    hit: (p) => pick(p) >= 0,
    down: (p) => { const e = pick(p); if (e < 0) return false; drag = { e, s0: along(e, p), u0: winds[e] }; fig.canvas.style.cursor = 'grabbing'; fig.render(); },
    move: (p) => { if (!drag) return; const u = clamp(Math.round((drag.u0 + 2 * (along(drag.e, p) - drag.s0) / geometry.k) * 2) / 2, -MAX, MAX); if (u !== winds[drag.e]) { winds[drag.e] = u; preset.value = 'own'; configure(); launch(); } fig.render(); },
    up: (p) => { drag = null; fig.canvas.style.cursor = pick(p) >= 0 ? 'grab' : ''; fig.render(); },
  });
  fig.canvas.addEventListener('pointermove', (e) => { if (drag) return; const r = fig.canvas.getBoundingClientRect(); fig.canvas.style.cursor = pick({ x: e.clientX - r.left, y: e.clientY - r.top }) >= 0 ? 'grab' : ''; });

  configure();
  reset();
  fig.play(true);
}
