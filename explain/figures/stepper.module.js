import { Figure, slider, choice, buttons, legend, caption, readout, text, arrow, anomalyColor, clamp, MUTED, GRID, LINE, ACCENT } from '../runtime.module.js';

const TURN = 2 * Math.PI, LIMIT = 2 * Math.SQRT2, XMAX = 3.5, YMAX = 2, WIDE = 520, TURNS = 3, ESCAPE = 3, HOLD = 1.2, ROOM = 1.4, METHODS = ['euler', 'rk4'], HALO = 'rgba(20, 20, 22, 0.85)';
const NAMES = { euler: 'Euler', rk4: 'Runge–Kutta 4' };
const tendency = (u, v) => [v, -u];

export function stepOnce(method, [u, v], x) {
  const k1 = tendency(u, v);
  if (method === 'euler') return [u + x * k1[0], v + x * k1[1]];
  const k2 = tendency(u + x / 2 * k1[0], v + x / 2 * k1[1]), k3 = tendency(u + x / 2 * k2[0], v + x / 2 * k2[1]), k4 = tendency(u + x * k3[0], v + x * k3[1]);
  return [u + x / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]), v + x / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])];
}

export const amplification = (method, x) => Math.hypot(...stepOnce(method, [1, 0], x));
const SAMPLES = Array.from({ length: 351 }, (_, i) => i * XMAX / 350);
const CURVES = Object.fromEntries(METHODS.map((m) => [m, SAMPLES.map((x) => amplification(m, x))]));
const factor = (a) => `×${a.toFixed(clamp(Math.ceil(-Math.log10(Math.abs(a - 1))) + 1, 3, 7))}`;

export function mountStepper(root, mark = null) {
  const controls = root.querySelector('.controls');
  let x = 0.25, show = 'both', n = 0, clock = 0, hold = 0, reach = 1, zoom = ROOM, paths = {};
  const shown = (m) => show === 'both' || show === m;
  const color = (m, alpha = 1) => anomalyColor(m === 'euler' ? 1 : -1, alpha);
  const tone = (m) => color(m, shown(m) ? 1 : 0.3);
  const last = (m) => paths[m][paths[m].length - 1];
  const speed = (m) => Math.hypot(...last(m));
  const running = (m) => shown(m) && speed(m) <= ESCAPE;
  const perStep = () => Math.min(0.5, 4 * x / TURN);

  const fig = new Figure(root, { height: 340, minHeight: 300, step, draw });
  const fit = fig.resize.bind(fig);
  fig.resize = () => { const tall = fig.stage.clientWidth < WIDE; fig.height = tall ? 560 : 340; fig.minHeight = tall ? 560 : 300; fit(); };
  caption(root, 'Each dot is one step, the only place the computer knows the wind; the straight lines just join the dots, and earlier turns fade. The open ring marks the exact wind after the same time.');
  legend(root, [['dash', 'the exact answer: the wind turns at a steady speed', MUTED], ['line', 'forward Euler', 'warm'], ['line', 'Runge–Kutta 4', 'cool']]);
  slider(controls, { label: 'Time step', min: 0.05, max: 3.2, step: 0.05, value: x, format: (v) => `f Δt = ${v.toFixed(2)}, ${TURN / v >= 6 ? `1/${Math.round(TURN / v)}` : (v / TURN).toFixed(2)} turn`, onInput: (v) => { x = v; restart(); } });
  choice(controls, { label: 'Show', options: [['both', 'both'], ['forward Euler', 'euler'], ['Runge–Kutta 4', 'rk4']], value: show, onChange: (v) => { show = v; restart(); } });
  buttons(controls, [['Restart', restart]]);
  const out = readout(controls);

  function restart() {
    n = 0; clock = 0; hold = 0; reach = 1;
    paths = { euler: [[0, 1]], rk4: [[0, 1]] };
    report();
  }

  function report() {
    out.set([['steps taken', `${n}`], ...METHODS.filter(shown).flatMap((m) => [[`${NAMES[m]}: speed`, `×${speed(m).toFixed(2)}`], ['each step', factor(amplification(m, x))]])]);
  }

  function advance() {
    n++;
    for (const m of METHODS) if (speed(m) <= ESCAPE) paths[m].push(stepOnce(m, last(m), x));
    const live = METHODS.filter(running);
    if (live.length) reach = Math.max(1, ...live.map(speed));
    if (n * x >= TURNS * TURN - 1e-9 || !live.length) hold = HOLD;
    report();
  }

  function step(dt) {
    zoom += (Math.max(ROOM, 1.12 * reach) - zoom) * Math.min(1, dt * 5);
    if (hold > 0) { hold -= dt; if (hold <= 0) restart(); }
    else for (clock += dt; clock >= perStep() && !hold; clock -= perStep()) advance();
  }

  function plane(ctx, [x0, y0, x1, y1]) {
    const cx = (x0 + x1) / 2, cy = (y0 + y1) / 2 + 6, radius = Math.min(x1 - x0, y1 - y0) / 2 - 26, scale = radius / zoom;
    const sx = (u) => cx + u * scale, sy = (v) => cy - v * scale, dot = clamp(x * scale * 0.15, 1.2, 3);
    ctx.strokeStyle = GRID; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(x0 + 8, cy); ctx.lineTo(x1 - 8, cy); ctx.moveTo(cx, cy + radius + 12); ctx.lineTo(cx, cy - radius - 12); ctx.stroke();
    text(ctx, `step ${n}`, x0 + 10, y0 + 14, { color: MUTED, size: 11 });
    ctx.strokeStyle = MUTED; ctx.setLineDash([4, 5]); ctx.beginPath(); ctx.arc(cx, cy, scale, 0, TURN); ctx.stroke(); ctx.setLineDash([]);
    ctx.save(); ctx.beginPath(); ctx.rect(x0, y0, x1 - x0, y1 - y0); ctx.clip();
    for (const m of METHODS.filter(shown)) {
      const cut = Math.max(0, paths[m].length - 1 - Math.ceil(TURN / x));
      for (const [points, faint] of [[paths[m].slice(0, cut + 1), 0.3], [paths[m].slice(cut), 1]]) {
        ctx.strokeStyle = color(m, 0.55 * faint); ctx.lineWidth = 1.5; ctx.beginPath();
        points.forEach(([u, v], i) => { if (i) ctx.lineTo(sx(u), sy(v)); else ctx.moveTo(sx(u), sy(v)); });
        ctx.stroke();
        ctx.fillStyle = color(m, faint);
        for (const [u, v] of points) { ctx.beginPath(); ctx.arc(sx(u), sy(v), dot, 0, TURN); ctx.fill(); }
      }
    }
    for (const m of METHODS.filter(shown)) { const [u, v] = last(m); arrow(ctx, cx, cy, sx(u), sy(v), { color: color(m), width: 2, head: 8 }); }
    ctx.strokeStyle = MUTED; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(sx(Math.sin(n * x)), sy(Math.cos(n * x)), 5, 0, TURN); ctx.stroke();
    ctx.restore();
    text(ctx, 'east wind u', x1 - 8, cy + 11, { align: 'right', color: MUTED, size: 10, halo: HALO });
    text(ctx, 'north wind v', cx + 6, cy - radius - 8, { color: MUTED, size: 10, halo: HALO });
  }

  function chart(ctx, [x0, y0, x1, y1]) {
    const left = x0 + 40, right = x1 - 16, top = y0 + 34, bottom = y1 - 34;
    const px = (s) => left + s / XMAX * (right - left), py = (a) => bottom - a / YMAX * (bottom - top);
    text(ctx, '|G|: how much one step multiplies the speed', x0 + 10, y0 + 14, { color: MUTED, size: 11 });
    ctx.lineWidth = 1;
    for (let a = 0; a <= YMAX; a += 0.5) { ctx.strokeStyle = a === 0 || a === 1 ? LINE : GRID; ctx.beginPath(); ctx.moveTo(left, py(a)); ctx.lineTo(right, py(a)); ctx.stroke(); text(ctx, `${a}`, left - 6, py(a), { align: 'right', color: MUTED, size: 10 }); }
    for (let s = 0; s <= 3; s++) text(ctx, `${s}`, px(s), bottom + 12, { align: 'center', color: MUTED, size: 10 });
    text(ctx, 'time step ω Δt', right, bottom + 26, { align: 'right', color: MUTED, size: 10 });
    ctx.strokeStyle = 'rgba(255,255,255,0.3)'; ctx.lineWidth = 1; ctx.setLineDash([2, 4]); ctx.beginPath(); ctx.moveTo(px(x), top); ctx.lineTo(px(x), bottom); ctx.stroke(); ctx.setLineDash([]);
    ctx.strokeStyle = color('rk4', 0.6); ctx.setLineDash([3, 4]); ctx.beginPath(); ctx.moveTo(px(LIMIT), top); ctx.lineTo(px(LIMIT), bottom); ctx.stroke(); ctx.setLineDash([]);
    ctx.save(); ctx.beginPath(); ctx.rect(left, top - 3, right - left, bottom - top + 3); ctx.clip();
    for (const m of METHODS) { ctx.strokeStyle = tone(m); ctx.lineWidth = 2; ctx.beginPath(); CURVES[m].forEach((a, i) => { if (i) ctx.lineTo(px(SAMPLES[i]), py(a)); else ctx.moveTo(px(SAMPLES[i]), py(a)); }); ctx.stroke(); }
    ctx.restore();
    text(ctx, 'keeps its speed', px(2.72), py(1) - 9, { align: 'right', color: MUTED, size: 10, halo: HALO });
    text(ctx, 'RK4’s limit', px(LIMIT) - 5, top + 16, { align: 'right', color: MUTED, size: 10, halo: HALO });
    text(ctx, '2√2 ≈ 2.83', px(LIMIT) - 5, top + 29, { align: 'right', color: MUTED, size: 10, halo: HALO });
    if (mark) {
      const mx = px(clamp(mark.x, 0, XMAX));
      ctx.strokeStyle = ACCENT; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(mx, bottom + 4); ctx.lineTo(mx, bottom - 8); ctx.stroke();
      text(ctx, mark.label, mx, bottom - 16, { align: mx < left + 50 ? 'left' : mx > right - 50 ? 'right' : 'center', color: ACCENT, size: 10, halo: HALO });
    }
    text(ctx, 'Euler', px(1) - 8, py(1.5), { align: 'right', color: tone('euler'), size: 11, halo: HALO });
    text(ctx, 'RK4', px(2.42), py(0.5) + 14, { align: 'center', color: tone('rk4'), size: 11, halo: HALO });
    for (const m of METHODS) {
      const a = amplification(m, x);
      ctx.fillStyle = ctx.strokeStyle = tone(m); ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(px(x), py(Math.min(a, YMAX)), 4.5, 0, TURN);
      if (a > YMAX) ctx.stroke(); else ctx.fill();
    }
  }

  function draw(ctx, w, h) {
    const tall = w < WIDE, split = Math.round(tall ? h * 0.56 : w * 0.48);
    plane(ctx, tall ? [0, 0, w, split] : [0, 0, split, h]);
    chart(ctx, tall ? [0, split, w, h] : [split, 0, w, h]);
  }

  restart();
  fig.play(true);
}
