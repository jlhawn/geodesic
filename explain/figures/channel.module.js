import { Figure, choice, buttons, legend, caption, readout, text, arrow, termColor, clamp, INK, MUTED, GRID, LINE } from '../runtime.module.js';

const KM = 1000, LENGTH = 1000 * KM, TOP = 40, KMAX = 800, HOUR = 3600, RUN = 6 * HOUR, RATE = HOUR, FRONT = 0.25, SAMPLES = 240, GAP = 40 * KM, SOFT = 2 / HOUR, WIDE = 520, REACH = 18;
const HALO = 'rgba(20, 20, 22, 0.85)', AIR = 'rgba(255, 255, 255, 0.7)', TURN = 2 * Math.PI;
const KNOTS = [0, 1, 2, 3].map((i) => i * LENGTH / 3);
const PRESETS = { up: [10, 50 / 3, 70 / 3, 30], down: [30, 70 / 3, 50 / 3, 10], steady: [20, 20, 20, 20], jet: [10, 35, 35, 10] };
const signed = (v, n) => `${Number(v.toFixed(n)) === 0 ? '' : v < 0 ? '−' : '+'}${Math.abs(v).toFixed(n)}`;

export function speedCurve(values) {
  const h = LENGTH / 3, d = [0, 1, 2].map((i) => (values[i + 1] - values[i]) / h);
  const end = (a, b) => { const p = 1.5 * a - 0.5 * b; return p * a <= 0 ? 0 : Math.abs(p) > 2 * Math.abs(a) ? 2 * a : p; };
  const m = [end(d[0], d[1]), ...[1, 2].map((i) => (Math.sign(d[i - 1]) + Math.sign(d[i])) * Math.min(Math.abs(d[i - 1]), Math.abs(d[i]), Math.abs(d[i - 1] + d[i]) / 4)), end(d[2], d[1])];
  return (x) => {
    const i = clamp(Math.floor(x / h), 0, 2), s = clamp(x / h - i, 0, 1), s2 = s * s, s3 = s2 * s, a = values[i], b = values[i + 1], ma = m[i] * h, mb = m[i + 1] * h;
    return [(2 * s3 - 3 * s2 + 1) * a + (s3 - 2 * s2 + s) * ma + (3 * s2 - 2 * s3) * b + (s3 - s2) * mb, ((6 * s2 - 6 * s) * (a - b) + (3 * s2 - 4 * s + 1) * ma + (3 * s2 - 2 * s) * mb) / h];
  };
}

export function carried(curve, t) {
  const inflow = curve(0)[0];
  return (x) => {
    if (t > 0 && x <= inflow * t) return [inflow, 0];
    let lo = 0, hi = x;
    for (let n = 0; n < 32; n++) { const mid = (lo + hi) / 2; if (mid + curve(mid)[0] * t < x) lo = mid; else hi = mid; }
    const [u, du] = curve((lo + hi) / 2);
    return [u, du / (1 + du * t)];
  };
}

export function frontTime(curve) {
  let best = Infinity;
  for (let i = 0; i <= 1000; i++) { const x0 = i * LENGTH / 1000, [u, du] = curve(x0); if (du < 0) { const t = (FRONT - 1) / du; if (x0 + u * t <= LENGTH) best = Math.min(best, t); } }
  return best;
}

export function mountChannel(root) {
  const controls = root.querySelector('.controls'), B = termColor('b');
  let values = [...PRESETS.up], probe = LENGTH / 2, start = probe, t = 0, end = RUN, curve = null, dots = [], drag = null, offset = 0, geometry = null;
  const fig = new Figure(root, { height: 440, minHeight: 440, step, draw });
  caption(root, `The wind blows from left to right along a straight channel 1,000 km long. Drag the four points to shape it, and drag along the channel to move the probe. “Let the wind carry itself” runs this term on its own, an hour every second, for up to ${RUN / HOUR} hours, and a ring follows the air that was at the probe when the run began.`);
  const key = legend(root, []);
  const preset = choice(controls, { label: 'Wind', options: [['speeding up', 'up'], ['slowing down', 'down'], ['steady', 'steady'], ['a jet streak', 'jet']], value: 'up', onChange: (v) => { values = [...PRESETS[v]]; reset(); }, span: true });
  buttons(controls, [['Let the wind carry itself', run], ['Reset', reset]]);
  const out = readout(controls), box = controls.querySelector('.readout');

  function shape() {
    curve = speedCurve(values);
    end = Math.min(RUN, frontTime(curve));
    const inflow = curve(0)[0];
    dots = [];
    for (const lane of [0, 1]) for (let x0 = -RUN * TOP + lane * GAP / 2; x0 <= LENGTH; x0 += GAP) dots.push({ x0, u0: x0 < 0 ? inflow : curve(x0)[0], lane });
  }

  function labels(started) {
    key.set([['arrow', 'wind', INK], ['dots', 'bits of air, carried by the wind', AIR], ['line', 'kinetic energy per kilogram, K = ½u²', B], ['force', 'the term −∂K/∂x: how it changes the wind, pointing down the slope of K', B], ...(started ? [['dash', 'the wind at the start', MUTED]] : [])]);
  }

  function reset() { t = 0; fig.play(false); shape(); labels(false); fig.render(); }
  function run() { t = 0; start = probe; labels(true); fig.play(true); fig.render(); }
  function step(dt) { t = Math.min(end, t + dt * RATE); if (t >= end) fig.play(false); }

  function words(u, a) {
    if (u < 0.5) return 'barely changes, because the air is almost still and little new air arrives';
    if (Math.abs(a) < 0.005) return 'holds, because the air arriving from upstream is just as fast';
    return a < 0 ? 'slows, because slower air is arriving from upstream' : 'speeds up, because faster air is arriving from upstream';
  }

  function draw(ctx, w, h) {
    const narrow = w < WIDE, size = narrow ? 10 : 11, left = narrow ? 40 : 48, right = 14, span = w - left - right, n = narrow ? 5 : 8, pitch = span / n, pushes = narrow ? 6 : 10, gap = span / pushes;
    const wallTop = 26, wallBottom = 88, lane = (wallTop + wallBottom) / 2, sTop = 128, sBottom = 236, kTop = 274, kBottom = h - 32;
    const x = (m) => left + m / LENGTH * span, uy = (v) => sBottom - v / TOP * (sBottom - sTop), ky = (k) => kBottom - k / KMAX * (kBottom - kTop);
    const line = (x0, y0, x1, y1) => { ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke(); };
    const path = (points) => { ctx.beginPath(); points.forEach(([px, py], i) => { if (i) ctx.lineTo(px, py); else ctx.moveTo(px, py); }); ctx.stroke(); };
    const push = (px, py, a, width) => { const len = clamp(Math.sign(a) * 0.85 * gap * Math.tanh(Math.abs(a) / SOFT), left - 2 - px, w - 4 - px); if (Math.abs(len) >= 4) arrow(ctx, px, py, px + len, py, { color: B, width, head: 7, dash: [5, 4], open: true }); };
    const above = (k) => (ky(k) - 14 > kTop + 2 ? -14 : 14);
    const field = carried(curve, t), along = Array.from({ length: SAMPLES + 1 }, (_, i) => i * LENGTH / SAMPLES), samples = along.map(field);
    const [pu, pdu] = field(probe), pk = pu * pu / 2, pa = -pu * pdu, px = x(probe), u0 = curve(start)[0], xp = start + u0 * t, followed = t > 0 && xp <= LENGTH, stopped = t >= end && end < RUN;
    geometry = { x, left, span, wallTop, wallBottom, sTop, sBottom, uy, marks: [uy(pu), ky(pk)] };

    ctx.lineWidth = 1;
    for (let v = 0; v <= TOP; v += 10) { ctx.strokeStyle = v ? GRID : LINE; line(left, uy(v), w - right, uy(v)); text(ctx, `${v}`, left - 10, uy(v), { align: 'right', color: MUTED, size: 10 }); }
    for (let k = 0; k <= KMAX; k += 200) { ctx.strokeStyle = k ? GRID : LINE; line(left, ky(k), w - right, ky(k)); text(ctx, `${k}`, left - 10, ky(k), { align: 'right', color: MUTED, size: 10 }); }
    for (let m = 0; m <= LENGTH; m += 250 * KM) { ctx.strokeStyle = GRID; line(x(m), sTop, x(m), sBottom); line(x(m), kTop, x(m), kBottom); text(ctx, m < LENGTH ? `${m / KM}` : '1,000 km', x(m), kBottom + 14, { align: m < LENGTH ? 'center' : 'right', color: MUTED, size: 10 }); }
    text(ctx, 'wind speed u, m/s', left, sTop - 16, { color: MUTED, size });
    text(ctx, 'kinetic energy K = ½u², J per kg', left, kTop - 16, { color: MUTED, size });
    if (t > 0) text(ctx, `after ${(t / HOUR).toFixed(1)} hours${stopped ? ', stopped' : ''}`, w - right, sTop - 16, { align: 'right', color: INK, size });

    ctx.fillStyle = 'rgba(255, 255, 255, 0.04)'; ctx.fillRect(left, wallTop, span, wallBottom - wallTop);
    ctx.strokeStyle = LINE; ctx.lineWidth = 1.5; line(left, wallTop, w - right, wallTop); line(left, wallBottom, w - right, wallBottom);
    for (let i = 0; i < n; i++) { const m = (i + 0.5) * LENGTH / n, len = field(m)[0] / TOP * 0.85 * pitch; arrow(ctx, x(m) - len / 2, lane, x(m) + len / 2, lane, { color: INK, width: 1.8, head: 7 }); }
    ctx.fillStyle = AIR;
    for (const d of dots) { const m = d.x0 + d.u0 * t; if (m < 0 || m > LENGTH) continue; ctx.beginPath(); ctx.arc(x(m), d.lane ? wallBottom - 11 : wallTop + 11, 2.2, 0, TURN); ctx.fill(); }

    if (t > 0) { ctx.strokeStyle = MUTED; ctx.lineWidth = 1.5; ctx.setLineDash([4, 4]); path(along.map((m) => [x(m), uy(curve(m)[0])])); ctx.setLineDash([]); }
    ctx.strokeStyle = INK; ctx.lineWidth = 2; path(samples.map(([u], i) => [x(along[i]), uy(u)]));
    ctx.strokeStyle = B; ctx.lineWidth = 2; path(samples.map(([u], i) => [x(along[i]), ky(u * u / 2)]));
    for (let i = 0; i < pushes; i++) { const at = (i + 0.5) * LENGTH / pushes, [u, du] = field(at); if (Math.abs(x(at) - px) > 0.3 * gap) push(x(at), ky(u * u / 2) + above(u * u / 2), -u * du, 1.6); }

    ctx.strokeStyle = 'rgba(255, 255, 255, 0.3)'; ctx.lineWidth = 1; ctx.setLineDash([2, 4]); line(px, sTop, px, sBottom); line(px, kTop, px, kBottom); ctx.setLineDash([]);
    ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; line(px, wallTop, px, wallBottom);
    ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(px, wallTop, 5, 0, TURN); ctx.fill();
    text(ctx, 'probe', clamp(px, left + 14, w - right - 14), wallTop - 15, { align: 'center', color: MUTED, size: 10 });
    for (const py of geometry.marks) { ctx.beginPath(); ctx.arc(px, py, 4, 0, TURN); ctx.fill(); }
    push(px, ky(pk) + Math.sign(above(pk)) * 5, pa, 2.6);

    KNOTS.forEach((m, i) => {
      const kx = x(m), kyy = uy(values[i]);
      ctx.beginPath(); ctx.arc(kx, kyy, t > 0 ? 4.5 : 6, 0, TURN);
      if (t > 0) { ctx.strokeStyle = MUTED; ctx.lineWidth = 1.5; ctx.stroke(); }
      else { ctx.fillStyle = '#fff'; ctx.fill(); ctx.strokeStyle = HALO; ctx.lineWidth = 2; ctx.stroke(); }
      if (drag === i) text(ctx, `${Math.round(values[i])} m/s`, kx + (i === 3 ? -12 : 12), kyy, { align: i === 3 ? 'right' : 'left', color: '#fff', size, halo: HALO });
    });

    if (followed) for (const fy of [lane, uy(u0), ky(u0 * u0 / 2)]) { ctx.beginPath(); ctx.arc(x(xp), fy, 6, 0, TURN); ctx.fillStyle = HALO; ctx.fill(); ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.stroke(); }

    out.set([
      ['at the probe, the wind', `${pu.toFixed(1)} m/s`], ['kinetic energy', `${pk.toFixed(0)} J per kg`], ['speed change along the channel', `${signed(pdu * 1e5, 1)} m/s per 100 km`],
      ['the term −∂K/∂x, which for straight flow is the advection −u ∂u/∂x,', `${signed(pa * HOUR, 2)} m/s per hour`], ['so the wind here', words(pu, pa * HOUR)],
      ...(t > 0 ? [['the ring, the air that was at the probe at the start,', `still ${u0.toFixed(1)} m/s, ${followed ? `now ${((xp - start) / KM).toFixed(0)} km downstream` : 'now past the end of the channel'}`]] : []),
      ...(stopped ? [['the run stopped because', 'faster air is catching up with slower air ahead of it, and the wind is about to jump from fast to slow']] : []),
    ]);
    const bold = box.querySelectorAll('b');
    bold[1].style.color = bold[3].style.color = B;
  }

  function grab({ x: px, y: py }) {
    if (!geometry) return null;
    const { x, uy, wallTop, wallBottom, marks } = geometry;
    let best = null, reach = 22;
    KNOTS.forEach((m, i) => { const d = Math.hypot(x(m) - px, uy(values[i]) - py); if (d < reach) { reach = d; best = i; } });
    if (best !== null) return best;
    return (py > wallTop - 24 && py < wallBottom + 10) || (Math.abs(px - x(probe)) < REACH && marks.some((my) => Math.abs(py - my) < REACH)) ? 'probe' : null;
  }

  function move({ x: px, y: py }) {
    const { left, span, sTop, sBottom } = geometry;
    if (drag === 'probe') { probe = clamp((px - left) / span, 0, 1) * LENGTH; fig.render(); }
    else if (drag !== null) { values[drag] = clamp(Math.round((sBottom - py - offset) / (sBottom - sTop) * TOP), 0, TOP); preset.value = ''; reset(); }
  }

  fig.pointer({
    hit: (p) => grab(p) !== null,
    down: (p) => { drag = grab(p); if (drag === 'probe') move(p); else { offset = geometry.uy(values[drag]) - p.y; fig.render(); } },
    move,
    up: () => { drag = null; fig.render(); },
  });
  fig.canvas.addEventListener('pointermove', (e) => {
    if (drag !== null) return;
    const r = fig.canvas.getBoundingClientRect(), g = grab({ x: e.clientX - r.left, y: e.clientY - r.top });
    fig.canvas.style.cursor = g === null ? '' : g === 'probe' ? 'ew-resize' : 'ns-resize';
  });

  reset();
}
