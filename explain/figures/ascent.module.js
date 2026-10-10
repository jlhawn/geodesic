import { Figure, slider, buttons, legend, caption, readout, text, arrow, termColor, clamp, INK, MUTED, GRID, LINE, ACCENT } from '../runtime.module.js';
import { T0 } from '../physics.module.js';

const KM = 1000, HOUR = 3600, DAY = 86400, TOP = 12 * KM, PEAK = 6 * KM, RUN = 6 * HOUR, SPEED = 1.5 * HOUR, WIDE = 520, WAVE = Math.PI / TOP;
const LOW = 280, HIGH = 365, EVERY = 5, NEAR = 0.05, REACH = 1.6 / 5, SURFACE = 'rgba(255,255,255,0.22)', HALO = 'rgba(20, 20, 22, 0.85)';
const LEVELS = [2, 4, 6, 8, 10].map((z) => z * KM), BEADS = Array.from({ length: 11 }, (_, k) => (k + 1) * KM), AXES = [5, 10, 15, 20, 30, 40];
const shape = (z) => Math.sin(WAVE * z);
const minus = (v, d) => `${v < 0 ? '−' : ''}${Math.abs(v).toFixed(d)}`;
const signed = (v, d = 1) => (Math.abs(v) < 0.5 * 10 ** -d ? (0).toFixed(d) : `${v < 0 ? '−' : '+'}${Math.abs(v).toFixed(d)}`);
const motion = (v) => (Math.abs(v) < 0.005 ? 'still' : `${v > 0 ? 'rising' : 'sinking'} ${Math.abs(v).toFixed(Math.abs(v) < 0.995 ? 2 : 1)} cm/s`);
const warming = (v) => (v > 0 ? `${v.toFixed(1)} K per day` : v < 0 ? `cooling, ${(-v).toFixed(1)} K per day` : 'none');

export function column(stability, rise, heating) {
  const start = (z) => T0 + stability * z;
  const parcel = (z0, t) => {
    if (!rise) return [z0, start(z0) + heating * shape(z0) * t];
    const z = 2 / WAVE * Math.atan(Math.tan(WAVE * z0 / 2) * Math.exp(WAVE * rise * t));
    return [z, start(z0) + heating / rise * (z - z0)];
  };
  const theta = (z, t) => {
    if (!rise) return start(z) + heating * shape(z) * t;
    const z0 = 2 / WAVE * Math.atan(Math.tan(WAVE * z / 2) * Math.exp(-WAVE * rise * t));
    return start(z0) + heating / rise * (z - z0);
  };
  const slope = (z, t) => {
    if (!rise) return stability + heating * WAVE * Math.cos(WAVE * z) * t;
    const u = Math.tan(WAVE * z / 2), e = Math.exp(-WAVE * rise * t), d = e * (1 + u * u) / (1 + u * u * e * e);
    return stability * d + heating / rise * (1 - d);
  };
  const height = (th, t) => {
    if (th <= theta(0, t) || th >= theta(TOP, t)) return null;
    let a = 0, b = TOP;
    for (let i = 0; i < 40; i++) { const m = (a + b) / 2; if (theta(m, t) < th) a = m; else b = m; }
    return (a + b) / 2;
  };
  const rates = (z, t) => { const c = -rise * shape(z) * slope(z, t), h = heating * shape(z); return [c, h, c + h]; };
  return { start, parcel, theta, slope, height, rates };
}

export function mountAscent(root) {
  const controls = root.querySelector('.controls'), C = termColor('c'), E = termColor('e');
  let stability = 3.3, rise = 2, heating = 0, probe = PEAK, t = 0, clock = false, dragged = false, geometry = null, model = null;
  const build = () => { model = column(stability / KM, rise / 100, heating / DAY); };
  const perDay = (z) => model.rates(z, t).map((v) => v * DAY);
  const fit = () => { const m = Math.max(...[...LEVELS, probe].flatMap((z) => perDay(z).map(Math.abs))); return AXES.find((a) => a >= m * 1.05) ?? AXES[AXES.length - 1]; };
  const balanced = () => rise !== 0 && Math.abs(heating - rise / 100 * stability / KM * DAY) < NEAR;
  build();
  let axis = fit(), goal = axis;

  const fig = new Figure(root, { height: 440, minHeight: 440, step, draw });
  caption(root, `A column of air 12 km tall, seen from the side. The white line is its θ at each height, and the thin level lines are the surfaces where θ is 290 K, 295 K and so on. The air rises or sinks fastest at 6 km and not at all at the ground or the top. Where its speed changes with height, air flows in from the sides or out to them with the same θ as the column, so only the vertical motion and the heating change θ. The heating has the same shape, and below zero it is cooling, the way air loses heat by radiating to space. On the right, at five heights, is how fast each term changes θ, in K per day, and the two together. Drag the ring to read another height. Running shows 6 hours, ${SPEED / HOUR} hours each second, and the dots ride with the air.`);
  legend(root, [['line', 'θ at each height', INK], ['dash', 'θ at the start', MUTED], ['faint', 'surfaces of equal θ, every 5 K', LINE], ['varrow', 'the air rising or sinking', INK], ['force', 'the vertical term, −w ∂θ/∂z', C], ['force', 'the heating, H', E], ['force', 'the two together', INK]]);
  slider(controls, { label: 'Stability', min: 1, max: 6, step: 0.1, value: stability, format: (v) => `θ climbs ${v.toFixed(1)} K per km`, onInput: (v) => { stability = v; settle(); } });
  slider(controls, { label: 'Vertical motion at 6 km', min: -5, max: 5, step: 0.1, value: rise, format: (v) => (v ? motion(v) : 'none'), onInput: (v) => { rise = v; settle(); } });
  slider(controls, { label: 'Heating at 6 km', min: -2, max: 6, step: 0.1, value: heating, format: warming, onInput: (v) => { heating = v; settle(); } });
  buttons(controls, [['Run 6 hours', () => { t = 0; clock = true; settle(); }], ['Reset', () => { t = 0; clock = false; settle(); }]]);
  const out = readout(controls);

  function report(every) {
    const [c, e, n] = perDay(probe), w = rise * shape(probe), s = model.slope(probe, t) * KM, lift = w / 100 * DAY / KM, th = model.theta(probe, t), d = th - model.start(probe), even = balanced(), at = `${probe / KM} km`;
    const verdict = even ? '' : !rise && !heating ? ', so nothing changes' : Math.abs(n) < NEAR ? ', almost no change' : `, so the air at ${at} ${n < 0 ? 'cools' : 'warms'}`;
    const since = Math.abs(d) < 0.005 ? 'the same as' : `${Math.abs(d).toFixed(2)} K ${d < 0 ? 'lower' : 'higher'} than`;
    out.set([
      [`at ${at}`, `air ${motion(w)}, θ climbing ${s.toFixed(2)} K per km`],
      ['the vertical term', Math.abs(lift) < 0.005 ? `${signed(c)} K per day` : `−(${minus(lift, 2)} km a day × ${s.toFixed(2)} K per km) = ${signed(c)} K per day`],
      ['the heating', `${signed(e)} K per day`],
      ['the two together', `${signed(n)} K per day${verdict}`],
      [`θ at ${at}`, t > 0 ? `${th.toFixed(2)} K after ${(t / HOUR).toFixed(1)} hours, ${since} at the start` : `${th.toFixed(2)} K`],
      ...(even ? [['in balance', rise > 0 ? 'rising air cools each level exactly as fast as the heating warms it: the balance of the tropics, where condensing water heats the rising air' : 'sinking air warms each level exactly as fast as the air cools by radiating to space: the balance of the subtropics, where dry air sinks over the deserts']] : []),
    ], every);
  }

  function settle() {
    build();
    goal = fit();
    fig.play(clock || axis !== goal);
    fig.render();
  }

  function step(dt) {
    const was = [t, axis];
    if (clock) { t = Math.min(RUN, t + dt * SPEED); clock = t < RUN; goal = fit(); }
    axis = Math.abs(goal - axis) < 0.05 ? goal : axis + (goal - axis) * Math.min(1, dt * 8);
    if (!clock && axis === goal) fig.play(false);
    return t !== was[0] || axis !== was[1];
  }

  function curve(ctx, x, y, at, color, width, dash = []) {
    ctx.strokeStyle = color; ctx.lineWidth = width; ctx.lineJoin = 'round'; ctx.setLineDash(dash); ctx.beginPath();
    for (let i = 0; i <= 120; i++) { const z = i * TOP / 120, px = x(model.theta(z, at)), py = y(z); if (i) ctx.lineTo(px, py); else ctx.moveTo(px, py); }
    ctx.stroke(); ctx.setLineDash([]);
  }

  function draw(ctx, w, h) {
    const narrow = w < WIDE, left = narrow ? 36 : 46, right = w - (narrow ? 30 : 38), top = 34, bottom = h - 44, gap = narrow ? 16 : 30;
    const mid = left + (right - left - gap) / 2, r0 = mid + gap, cx = (r0 + right) / 2, half = (right - r0) / 2, per = half / axis, km = (bottom - top) / 12;
    const x = (th) => left + (th - LOW) / (HIGH - LOW) * (mid - left), y = (z) => bottom - z / TOP * (bottom - top);
    geometry = { y, left, right, top, bottom };
    ctx.lineWidth = 1;
    for (let z = 0; z <= TOP; z += 2 * KM) {
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(left - 4, y(z)); ctx.lineTo(left, y(z)); ctx.moveTo(r0, y(z)); ctx.lineTo(right, y(z)); ctx.stroke();
      text(ctx, `${z / KM} km`, left - 7, y(z), { align: 'right', color: z === probe ? INK : MUTED, size: 10 });
    }
    for (let th = LOW; th <= HIGH; th += 20) {
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(x(th), top); ctx.lineTo(x(th), bottom); ctx.stroke();
      text(ctx, `${th}`, x(th), bottom + 16, { align: 'center', color: MUTED, size: 10 });
    }
    text(ctx, narrow ? 'θ, K' : 'potential temperature θ, K', (left + mid) / 2, bottom + 31, { align: 'center', color: MUTED, size: 10 });
    for (let th = Math.ceil(T0 / EVERY) * EVERY; th < HIGH; th += EVERY) {
      const z = model.height(th, t);
      if (z === null) continue;
      ctx.strokeStyle = SURFACE; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(left, y(z)); ctx.lineTo(mid, y(z)); ctx.stroke();
      if (th % 10 === 0 && y(z) > top + 12) { const low = x(th) < (left + mid) / 2; text(ctx, `${th} K`, low ? mid - 3 : left + 3, y(z) - 7, { align: low ? 'right' : 'left', color: MUTED, size: 10, halo: HALO }); }
    }
    for (let n = -2; n <= 2; n++) {
      const v = n * goal / 2, px = cx + v * per;
      if (Math.abs(px - cx) > half + 1) continue;
      ctx.strokeStyle = n ? GRID : LINE; ctx.beginPath(); ctx.moveTo(px, top - 4); ctx.lineTo(px, bottom); ctx.stroke();
      if (!narrow || n % 2 === 0) text(ctx, n ? signed(v, Number.isInteger(v) ? 0 : 1) : '0', px, bottom + 16, { align: 'center', color: MUTED, size: 10 });
    }
    text(ctx, narrow ? 'K per day' : 'change in θ, K per day', cx, bottom + 31, { align: 'center', color: MUTED, size: 10 });
    text(ctx, 'cools', cx - 6, top - 12, { align: 'right', color: MUTED, size: 10 });
    text(ctx, 'warms', cx + 6, top - 12, { color: MUTED, size: 10 });
    ctx.fillStyle = '#5a4a36'; ctx.fillRect(left, bottom, mid - left, 5);
    const py = y(probe), ring = x(model.theta(probe, t));
    ctx.strokeStyle = 'rgba(255,255,255,0.3)'; ctx.lineWidth = 1; ctx.setLineDash([2, 4]); ctx.beginPath(); ctx.moveTo(left, py); ctx.lineTo(right, py); ctx.stroke(); ctx.setLineDash([]);
    if (t > 0) curve(ctx, x, y, 0, MUTED, 1.5, [3, 4]);
    curve(ctx, x, y, t, INK, 2);
    for (const z of LEVELS) {
      const reach = rise * REACH * shape(z), low = z - Math.abs(reach) * KM / 2, len = reach * km, px = Math.max(left + 5, x(Math.min(model.theta(low, t), model.start(low))) - 10);
      if (Math.abs(len) >= 4) arrow(ctx, px, y(z) + len / 2, px, y(z) - len / 2, { color: INK, width: 1.5, head: 5 });
    }
    for (const z0 of BEADS) {
      const [z, th] = model.parcel(z0, t);
      ctx.fillStyle = HALO; ctx.strokeStyle = INK; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.arc(x(th), y(z), 3.2, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
    }
    ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(ring, py, 6.5, 0, Math.PI * 2); ctx.stroke();
    if (!dragged) text(ctx, 'drag', ring + 11, py - 12, { color: ACCENT, size: 11, halo: HALO });
    for (const z of LEVELS.includes(probe) ? LEVELS : [...LEVELS, probe]) {
      const yy = y(z), focus = z === probe, values = perDay(z);
      ctx.globalAlpha = focus ? 1 : 0.6;
      [[values[0], C, -7], [values[1], E, 0], [values[2], INK, 7]].forEach(([v, color, dy]) => {
        const len = clamp(v * per, -half - 4, half + 4);
        if (Math.abs(len) >= 1.5) arrow(ctx, cx, yy + dy, cx + len, yy + dy, { color, width: 2, head: 7, dash: [5, 4], open: true });
        else { ctx.fillStyle = color; ctx.beginPath(); ctx.arc(cx, yy + dy, 1.8, 0, Math.PI * 2); ctx.fill(); }
      });
      ctx.globalAlpha = 1;
      text(ctx, signed(values[2]), w - 4, yy + 7, { align: 'right', color: focus ? INK : MUTED, size: 10, weight: focus ? 600 : 400 });
    }
    text(ctx, `${(t / HOUR).toFixed(1)} hours`, 8, 12, { color: MUTED, size: 11 });
    if (balanced()) text(ctx, narrow ? 'balanced' : 'balanced: θ holds still at every height', cx, bottom - 10, { align: 'center', color: INK, size: 11, halo: HALO });
    report(fig.running ? 3 : 1);
  }

  fig.pointer({
    hit: ({ x: px, y: py }) => !!geometry && px >= geometry.left - 30 && Math.abs(py - geometry.y(probe)) <= 18,
    down: () => { dragged = true; fig.render(); },
    move: ({ y: py }) => { probe = clamp(Math.round((geometry.bottom - py) / (geometry.bottom - geometry.top) * TOP / KM), 1, 11) * KM; settle(); },
  });
}
