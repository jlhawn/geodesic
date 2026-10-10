import { Figure, slider, buttons, legend, caption, readout, text, arrow, termColor, clamp, INK, MUTED, GRID, LINE } from '../runtime.module.js';

const OMEGA = 7.2921e-5, HOUR = 3600, DAY = 86400, KM = 1000, DEG = Math.PI / 180, TURN = 2 * Math.PI, BASE = 4 * HOUR, LEAD = 1.5 * HOUR, SHORTEST = 3, LONGEST = 8, FIT = 0.92, WIDE = 520, STRIP = 76, EDGE = 8, FOOT = 24, CORNER = 190, DIAL = 12, GRAB = 30, MAX = 30, LABEL = 48, REACH = 0.28;
const LEVELS = Array.from({ length: 15 }, (_, k) => 75 * KM * Math.SQRT2 ** (k - 6)), SPACINGS = [2, 5, 10, 25, 50, 100, 200, 500, 1000].map((d) => d * KM);
const HALO = 'rgba(20, 20, 22, 0.85)', TRACE = 'rgba(255,255,255,0.8)';
const per = (x) => (Math.abs(x * 1e5) < 0.05 ? 'zero' : `${x < 0 ? '−' : ''}${Math.abs(x * 1e5).toFixed(1)} × 10⁻⁵ per second`);
const big = (v, below = 10) => (v < below ? v.toFixed(1) : Math.round(v).toLocaleString('en-US'));

function at(s, e, [u0, v0]) {
  const a = e * s, S = Math.abs(a) < 1e-6 ? s : Math.sin(a) / e, C = Math.abs(a) < 1e-6 ? a * s / 2 : (1 - Math.cos(a)) / e;
  return [u0 * S + v0 * C, v0 * S - u0 * C];
}

function wind(s, e, [u0, v0]) {
  const c = Math.cos(e * s), n = Math.sin(e * s);
  return [u0 * c + v0 * n, v0 * c - u0 * n];
}

export function mountTurning(root) {
  const controls = root.querySelector('.controls'), PUSH = termColor('a');
  let latitude = 45, spin = 0, speed = 10, heading = 0, t = 0, moving = false, touched = false, view = null, geo = null, run = null, cached = null, grab = null;
  const f = () => 2 * OMEGA * Math.sin(latitude * DEG);
  const cancels = () => spin !== 0 && Math.abs(spin + f() * 1e5) <= 0.05 + 1e-9;
  const zeta = () => (cancels() ? -f() : spin * 1e-5);
  const eta = () => zeta() + f();
  const start = () => [speed * Math.cos(heading), speed * Math.sin(heading)];
  const radius = (e, s = speed) => (e ? s / Math.abs(e) : Infinity);
  const period = (e) => (e ? TURN / Math.abs(e) : Infinity);
  const level = () => LEVELS.find((H) => 2 * radius(eta(), Math.max(speed, 1)) <= FIT * H) ?? LEVELS[LEVELS.length - 1];

  const fig = new Figure(root, { height: 420, minHeight: 340, step, draw });
  const fit = fig.resize.bind(fig);
  fig.resize = () => { const narrow = fig.stage.clientWidth < WIDE; fig.height = narrow ? 430 : 420; fig.minHeight = narrow ? 430 : 340; fit(); };
  caption(root, 'Seen from above. Suppose the parcel feels this term and nothing else, with ζ + f held at its starting value, as if the parcel never left its latitude. Drag the tip of the wind arrow, then let it go.');
  legend(root, [['arrow', 'the wind u', INK], ['force', 'the push −(ζ + f) ẑ × u', PUSH], ['dash', 'the path it will take', LINE], ['line', 'the path it has taken', TRACE]]);
  slider(controls, { label: 'Latitude', min: -90, max: 90, step: 1, value: latitude, format: (v) => (v > 0 ? `${v}° north` : v < 0 ? `${-v}° south` : 'the equator'), onInput: (v) => { latitude = v; reset(); } });
  slider(controls, { label: 'Spin of the flow around it, ζ', min: -15, max: 15, step: 0.1, value: spin, format: (v) => `${v > 0 ? '+' : v < 0 ? '−' : ''}${Math.abs(v).toFixed(1)} × 10⁻⁵ per second`, onInput: (v) => { spin = v; reset(); } });
  buttons(controls, [['Let it go', launch], ['Reset', reset]]);
  const out = readout(controls);

  function layout(w, h) {
    const top = w < WIDE ? STRIP : 0, mapH = h - top, half = Math.min(w, mapH) / 2 - 12;
    return { w, h, top, mapH, half, cx: w / 2, cy: top + mapH / 2, px: clamp(half / 40, 2.5, 4.5) };
  }

  function plan() {
    if (!geo) return null;
    const e = eta(), u = start(), H = level(), key = [e, u[0], u[1], H, geo.w, geo.h].join();
    if (cached?.key === key) return cached;
    const T = period(e), scale = geo.half / H, { w, h, top, cx, cy } = geo;
    const gone = (s, k) => { const [x, y] = at(s, e, u), [a, b] = wind(s, e, u), X = cx + (x + a * k) * scale, Y = cy - (y + b * k) * scale; return X < EDGE || X > w - EDGE || Y < top + EDGE || Y > h - FOOT || (Y < STRIP && (X < CORNER || X > w - CORNER)); };
    const first = (k) => {
      const ds = H / 150 / speed;
      for (let i = 1; i <= 4000; i++) {
        if (i * ds >= T) return T;
        if (!gone(i * ds, k)) continue;
        let lo = (i - 1) * ds, hi = i * ds;
        for (let n = 0; n < 30; n++) { const mid = (lo + hi) / 2; if (gone(mid, k)) hi = mid; else lo = mid; }
        return lo;
      }
      return 4000 * ds;
    };
    let D = Math.min(T, DAY);
    if (speed) {
      D = first(0);
      if (D < T) D = Math.min(D, first([...[LABEL, LABEL / 2, 0].map((m) => (geo.px + m / speed) / scale), 0].find((k) => !gone(0, k))));
    }
    return (cached = { key, e, u, T, D, closed: D === T && speed > 0, rate: Math.min(Math.max(BASE, D / LONGEST), D / SHORTEST) });
  }

  function report(every = 1) {
    const e = eta(), [u, v] = run ? wind(t, run.e, run.u) : start(), push = Math.abs(e) * speed * HOUR, T = period(e);
    out.set([
      ['the planet’s spin, f', per(f())],
      ['together, ζ + f', cancels() ? 'zero, ζ cancels f' : per(e)],
      ['the push', e && speed ? `${push.toFixed(push < 10 ? 2 : 1)} m/s per hour, to the ${e > 0 ? 'right' : 'left'} of the wind` : 'none'],
      ['the circle’s radius', !speed ? 'none, the air is still' : e ? `${big(radius(e) / KM)} km` : 'none, it goes straight'],
      ['once around', e ? `${big(T / HOUR, 100)} hours${T > 2 * DAY ? `, ${(T / DAY).toFixed(1)} days` : ''}` : 'never'],
      ['speed', `${Math.hypot(u, v).toFixed(1)} m/s${t > 0 ? ', the same as at the start' : ''}`],
    ], every);
  }

  function launch() {
    const p = plan();
    if (!p) return;
    run = p; t = 0; moving = true;
    report();
    fig.play(true);
  }

  function reset() {
    moving = false; run = null; t = 0;
    report();
    fig.render();
    fig.play(true);
  }

  function step(dt) {
    let changed = false;
    const target = level();
    if (view === null) view = target;
    if (view !== target) { const k = Math.log(target / view); view = Math.abs(k) < 0.01 ? target : view * Math.exp(k * Math.min(1, dt * 6)); changed = true; }
    if (moving) {
      t = Math.min(run.D, t + dt * run.rate);
      if (t >= run.D) moving = false;
      report(moving ? 10 : 1);
      changed = true;
    }
    if (!changed) fig.play(false);
    return changed;
  }

  function trace(ctx, p, s, sx, sy) {
    const n = clamp(Math.ceil(240 * s / p.D), 2, 240);
    ctx.beginPath();
    for (let i = 0; i <= n; i++) { const [x, y] = at(s * i / n, p.e, p.u); if (i) ctx.lineTo(sx(x), sy(y)); else ctx.moveTo(sx(x), sy(y)); }
    ctx.stroke();
  }

  function dial(ctx, x, y, value, name, note) {
    ctx.strokeStyle = GRID; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(x, y, DIAL + 5, 0, TURN); ctx.stroke();
    if (Math.abs(value) < 0.05) { ctx.fillStyle = MUTED; ctx.beginPath(); ctx.arc(x, y, 2, 0, TURN); ctx.fill(); }
    else {
      const sweep = (0.12 + 0.73 * Math.min(1, Math.abs(value) / MAX)) * Math.PI, ccw = value > 0;
      for (const a0 of [-Math.PI / 2, Math.PI / 2]) {
        const a1 = ccw ? a0 - sweep : a0 + sweep, ex = x + DIAL * Math.cos(a1), ey = y + DIAL * Math.sin(a1), tx = ccw ? Math.sin(a1) : -Math.sin(a1), ty = ccw ? -Math.cos(a1) : Math.cos(a1);
        ctx.strokeStyle = INK; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(x, y, DIAL, a0, a1, ccw); ctx.stroke();
        arrow(ctx, ex - tx * 6, ey - ty * 6, ex, ey, { color: INK, width: 1.5, head: 5 });
      }
    }
    text(ctx, name, x, y + DIAL + 15, { align: 'center', color: INK, size: 11 });
    text(ctx, note, x, y + DIAL + 28, { align: 'center', color: MUTED, size: 10 });
  }

  function tip() {
    const p = run ?? plan(), s = run ? t : 0, scale = geo.half / view, [x, y] = at(s, p.e, p.u), [u, v] = wind(s, p.e, p.u);
    return [geo.cx + x * scale + u * geo.px, geo.cy - y * scale - v * geo.px];
  }

  function draw(ctx, w, h) {
    geo = layout(w, h);
    if (view === null) view = level();
    const p = plan(), shown = run ?? p, s = run ? t : 0, { cx, cy, top, mapH, half, px } = geo, scale = half / view;
    const sx = (x) => cx + x * scale, sy = (y) => cy - y * scale, spacing = SPACINGS.find((d) => d * scale >= 48) ?? SPACINGS[SPACINGS.length - 1], gap = spacing * scale;
    ctx.save(); ctx.beginPath(); ctx.rect(0, top, w, mapH); ctx.clip();
    ctx.strokeStyle = GRID; ctx.lineWidth = 1; ctx.beginPath();
    for (let k = -Math.floor(cx / gap); k <= Math.floor(cx / gap); k++) { ctx.moveTo(sx(k * spacing), top); ctx.lineTo(sx(k * spacing), h); }
    for (let k = -Math.floor(mapH / 2 / gap); k <= Math.floor(mapH / 2 / gap); k++) { ctx.moveTo(0, sy(k * spacing)); ctx.lineTo(w, sy(k * spacing)); }
    ctx.stroke();
    ctx.strokeStyle = LINE; ctx.lineWidth = 1.5; ctx.setLineDash([4, 5]); trace(ctx, shown, shown.D, sx, sy); ctx.setLineDash([]);
    if (shown.closed) { const [u0, v0] = shown.u, ox = sx(v0 / shown.e), oy = sy(-u0 / shown.e); ctx.strokeStyle = MUTED; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(ox - 4, oy); ctx.lineTo(ox + 4, oy); ctx.moveTo(ox, oy - 4); ctx.lineTo(ox, oy + 4); ctx.stroke(); }
    if (s > 0) { ctx.strokeStyle = TRACE; ctx.lineWidth = 2; trace(ctx, shown, s, sx, sy); }
    const [x, y] = at(s, shown.e, shown.u), [u, v] = wind(s, shown.e, shown.u), X = sx(x), Y = sy(y), pushLength = Math.abs(shown.e) * speed * LEAD * px;
    const squeeze = pushLength > REACH * half ? REACH * half / pushLength : 1, ax = shown.e * v * LEAD * px * squeeze, ay = -shown.e * u * LEAD * px * squeeze;
    ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(X, Y, 5, 0, TURN); ctx.fill();
    if (Math.hypot(ax, ay) > 2) {
      arrow(ctx, X, Y, X + ax, Y - ay, { color: PUSH, width: 2, head: 8, dash: [5, 4], open: true });
      const push = Math.abs(shown.e) * speed * HOUR, bx = -u / speed, by = v / speed;
      text(ctx, `${push.toFixed(push < 10 ? 1 : 0)} m/s per hour`, X + ax + bx * 8, Y - ay + by * 8, { align: bx > 0.3 ? 'left' : bx < -0.3 ? 'right' : 'center', baseline: by > 0.3 ? 'top' : by < -0.3 ? 'bottom' : 'middle', color: PUSH, size: 10, halo: HALO });
    }
    const tx = X + u * px, ty = Y - v * px;
    if (speed) {
      arrow(ctx, X, Y, tx, ty, { color: INK, width: 2, head: 8 });
      text(ctx, `${speed} m/s`, tx + u / speed * 15, ty - v / speed * 15, { align: u > 2 ? 'left' : u < -2 ? 'right' : 'center', color: INK, size: 10, halo: HALO });
    }
    if (!moving) { ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(tx, ty, 11, 0, TURN); ctx.stroke(); }
    if (!touched && !run) text(ctx, 'drag the arrow’s tip', tx, ty + (ay > 1 ? 24 : -22), { align: 'center', color: MUTED, size: 11, halo: HALO });
    ctx.restore();
    const left = 14 + DIAL + 5, y0 = 10 + DIAL + 5;
    dial(ctx, left, y0, f() * 1e5, 'f', 'the planet');
    text(ctx, '+', left + 33, y0, { align: 'center', color: MUTED, size: 13 });
    dial(ctx, left + 66, y0, zeta() * 1e5, 'ζ', 'the flow');
    text(ctx, '=', left + 99, y0, { align: 'center', color: MUTED, size: 13 });
    dial(ctx, left + 132, y0, eta() * 1e5, 'ζ + f', 'together');
    const rate = shown.rate / HOUR, shownRate = rate < 10 ? +rate.toFixed(1) : Math.round(rate);
    text(ctx, `${(s / HOUR).toFixed(1)} h${run && !moving && run.closed ? ', once around' : ''}`, w - 10, 14, { align: 'right', color: INK, size: 11 });
    text(ctx, `time runs ${shownRate} hour${shownRate === 1 ? '' : 's'} per second`, w - 10, 29, { align: 'right', color: MUTED, size: 10 });
    text(ctx, 'north ↑', cx, top + 12, { align: 'center', color: MUTED, size: 11 });
    text(ctx, 'east →', w - 10, h - 12, { align: 'right', color: MUTED, size: 11 });
    text(ctx, `grid lines every ${spacing / KM} km`, 10, h - 12, { color: MUTED, size: 11 });
  }

  function aim({ x, y }) {
    const a = grab.u + (x - grab.x) / geo.px, b = grab.v - (y - grab.y) / geo.px;
    speed = clamp(Math.round(Math.hypot(a, b)), 0, MAX);
    if (speed) heading = Math.round(Math.atan2(b, a) / (5 * DEG)) * 5 * DEG;
    reset();
  }

  fig.pointer({
    hit: (p) => { if (!geo || view === null) return false; const [x, y] = tip(); return Math.hypot(p.x - x, p.y - y) <= GRAB; },
    down: (p) => { touched = true; const [u, v] = run ? wind(t, run.e, run.u) : start(); grab = { x: p.x, y: p.y, u, v }; aim(p); },
    move: aim,
  });
  report();
}
