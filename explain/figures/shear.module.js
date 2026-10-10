import { Figure, slider, buttons, legend, caption, readout, text, arrow, termColor, clamp, INK, MUTED, GRID, LINE, ACCENT } from '../runtime.module.js';

const KM = 1000, HOUR = 3600, TOP = 12 * KM, KINK = 10 * KM, CHECK = 5 * KM, RUN = 6 * HOUR, SPEED = 1.5 * HOUR, FILL = 0.82, ALOFT_MAX = 60, GROUND_MAX = 20, HALO = 'rgba(20, 20, 22, 0.85)';
const LEVELS = Array.from({ length: 12 }, (_, k) => (k + 1) * KM), PARCELS = Array.from({ length: 19 }, (_, k) => (k - 3) * KM);
const num = (v, d) => `${v < -0.5 * 10 ** -d ? '−' : ''}${Math.abs(v).toFixed(d)}`;
const signed = (v, d = 2) => (Math.abs(v) < 0.5 * 10 ** -d ? (0).toFixed(d) : `${v < 0 ? '−' : '+'}${Math.abs(v).toFixed(d)}`);

export function mountShear(root) {
  const controls = root.querySelector('.controls'), D = termColor('d');
  let ground = 5, aloft = 35, rise = 0.05, t = 0, clock = false, held = null, dragged = false, geometry = null;
  const base = (s) => (s <= 0 ? ground : s >= KINK ? aloft : ground + (aloft - ground) * s / KINK);
  const slope = (s) => (s > 0 && s < KINK ? (aloft - ground) / KINK : 0);
  const wind = (z, at = t) => base(z - rise * at);
  const push = (z) => -rise * slope(z - rise * t - Math.sign(rise));
  const fit = () => clamp(Math.ceil(Math.max(ground, aloft) / 10) * 10, 20, ALOFT_MAX);
  let axis = fit(), goal = axis;

  const fig = new Figure(root, { height: 440, minHeight: 400, step, draw });
  caption(root, `A column of air 12 km tall, seen from the side. The solid arrows are the wind blowing east at each kilometer: it changes steadily from the ground to 10 km and holds steady above. All the air above the ground rises or sinks at one speed, w, set below; measured in height rather than σ, the term is −w ∂u/∂z. Drag the rings to reshape the wind. Running it shows 6 hours, ${SPEED / HOUR} hours each second, and the dots ride with the air.`);
  legend(root, [['arrow', 'the wind blowing east at each height', INK], ['force', 'this term’s push, −w ∂u/∂z, drawn as the change it would make in 6 hours', D], ['faint', 'the wind at the start', LINE]]);
  const lower = slider(controls, { label: 'Wind at the ground', min: 0, max: GROUND_MAX, step: 1, value: ground, format: (v) => `${v} m/s`, onInput: (v) => { ground = v; settle(); } });
  const upper = slider(controls, { label: 'Wind at 10 km', min: 0, max: ALOFT_MAX, step: 1, value: aloft, format: (v) => `${v} m/s`, onInput: (v) => { aloft = v; settle(); } });
  slider(controls, { label: 'Vertical motion', min: -10, max: 10, step: 0.5, value: rise * 100, format: (v) => (v > 0 ? `rising ${v.toFixed(1)} cm/s` : v < 0 ? `sinking ${(-v).toFixed(1)} cm/s` : 'none'), onInput: (v) => { rise = v / 100; fig.render(); } });
  buttons(controls, [['Run 6 hours', () => { t = 0; clock = true; fig.play(true); }], ['Reset', () => { t = 0; clock = false; settle(); }]]);
  const out = readout(controls);

  function report(every) {
    const start = base(CHECK), end = wind(CHECK, RUN);
    out.set([
      ['shear below 10 km', `${num((aloft - ground) / 10, 1)} m/s per km`],
      ['push at 5 km', `${signed(push(CHECK) * HOUR)} m/s per hour`],
      ['wind at 5 km', `${num(wind(CHECK), 1)} m/s`],
      ['change at 5 km after 6 hours', `${signed(end - start)} m/s, to ${num(end, 1)} m/s`],
      ['in 6 hours the air', rise ? `${rise > 0 ? 'rises' : 'sinks'} ${num(Math.abs(rise) * RUN / KM, 2)} km` : 'stays at its height'],
    ], every);
  }

  function settle() {
    if (!held) goal = fit();
    fig.play(clock || axis !== goal);
    fig.render();
  }

  function step(dt) {
    const was = [t, axis];
    if (clock) { t = Math.min(RUN, t + dt * SPEED); clock = t < RUN; }
    axis = Math.abs(goal - axis) < 0.05 ? goal : axis + (goal - axis) * Math.min(1, dt * 8);
    if (!clock && axis === goal) fig.play(false);
    return t !== was[0] || axis !== was[1];
  }

  function curve(ctx, x, y, at, color, width) {
    const zs = Array.from({ length: 241 }, (_, i) => i * TOP / 240);
    for (const k of [rise * at, KINK + rise * at]) if (k > 0 && k < TOP) zs.push(k);
    zs.sort((a, b) => a - b);
    ctx.strokeStyle = color; ctx.lineWidth = width; ctx.lineJoin = 'round'; ctx.beginPath(); ctx.moveTo(x(ground), y(0));
    for (const z of zs) ctx.lineTo(x(wind(z, at)), y(z));
    ctx.stroke();
  }

  function draw(ctx, w, h) {
    const left = 40, right = w - 44, top = 32, bottom = h - 40, scale = (right - left) * FILL / axis, tick = axis > 30 ? 10 : 5;
    const x = (u) => left + u * scale, y = (z) => bottom - z / TOP * (bottom - top);
    geometry = { x, y, left, scale };
    ctx.lineWidth = 1;
    for (let u = 0; x(u) <= right && u <= ALOFT_MAX; u += tick) {
      ctx.strokeStyle = u ? GRID : LINE; ctx.beginPath(); ctx.moveTo(x(u), top - 6); ctx.lineTo(x(u), bottom); ctx.stroke();
      text(ctx, `${u}`, x(u), bottom + 16, { align: 'center', color: MUTED, size: 10 });
    }
    text(ctx, 'wind blowing east, m/s', right, bottom + 30, { align: 'right', color: MUTED, size: 10 });
    ctx.fillStyle = '#5a4a36'; ctx.fillRect(left - 8, bottom, right - left + 16, 6);
    if (t > 0) curve(ctx, x, y, 0, LINE, 1.5);
    curve(ctx, x, y, t, INK, 2);
    for (const z of LEVELS) {
      const py = y(z), rate = push(z) * HOUR, flat = Math.abs(rate) < 0.005;
      arrow(ctx, x(0), py, x(wind(z)), py, { color: INK, width: 1.5, head: 6 });
      text(ctx, flat ? '0' : signed(rate), w - 8, py, { align: 'right', color: flat ? MUTED : D, size: 10 });
      text(ctx, `${z / KM} km`, left - 6, py, { align: 'right', color: z === CHECK ? INK : MUTED, size: 10 });
    }
    if (t > 0) for (const s of PARCELS) {
      const z = s + rise * t;
      if (z < 0 || z > TOP) continue;
      ctx.fillStyle = HALO; ctx.strokeStyle = INK; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.arc(x(base(s)), y(z), 3.5, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
    }
    for (const z of LEVELS) {
      const py = y(z) - 6, tip = x(wind(z)), end = clamp(tip + push(z) * RUN * scale, left, right + 6);
      if (Math.abs(push(z) * HOUR) >= 0.005) arrow(ctx, tip, py, end, py, { color: D, width: 2, head: clamp(Math.abs(end - tip) / 3, 5, 7), dash: [5, 4], open: true });
    }
    for (const [u, z] of [[ground, 0], [aloft, KINK]]) { ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(x(u), y(z), 6.5, 0, Math.PI * 2); ctx.stroke(); }
    if (!dragged) { const roomy = x(aloft) + 40 < right; text(ctx, 'drag', x(aloft) + (roomy ? 11 : -10), y(KINK) + (roomy ? 13 : -16), { align: roomy ? 'left' : 'right', color: ACCENT, size: 11, halo: HALO }); }
    text(ctx, `${(t / HOUR).toFixed(1)} h`, 8, 12, { color: MUTED, size: 11 });
    if (rise) arrow(ctx, 52, rise > 0 ? 18 : 6, 52, rise > 0 ? 6 : 18, { color: INK, width: 1.5, head: 5 });
    text(ctx, rise ? `air ${rise > 0 ? 'rising' : 'sinking'} ${num(Math.abs(rise) * 100, 1)} cm/s` : 'no vertical motion', rise ? 60 : 44, 12, { color: MUTED, size: 11 });
    text(ctx, 'push, m/s per hour', w - 8, 12, { align: 'right', color: D, size: 10 });
    report(fig.running ? 3 : 1);
  }

  const pick = ({ x: px, y: py }) => (geometry ? [['ground', ground, 0], ['aloft', aloft, KINK]].find(([, u, z]) => Math.hypot(px - geometry.x(u), py - geometry.y(z)) <= 22)?.[0] : undefined);
  fig.pointer({
    hit: (p) => pick(p) !== undefined,
    down: (p) => { held = pick(p); dragged = true; fig.render(); },
    move: ({ x: px }) => {
      const v = Math.round((px - geometry.left) / geometry.scale);
      if (held === 'ground') { ground = clamp(v, 0, GROUND_MAX); lower.value = ground; } else { aloft = clamp(v, 0, ALOFT_MAX); upper.value = aloft; }
      settle();
    },
    up: () => { held = null; settle(); },
  });
}
