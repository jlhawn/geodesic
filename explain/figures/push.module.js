import { Figure, slider, choice, buttons, readout, legend, caption, text, arrow, anomalyColor, clamp, ACCENT, INK, MUTED, LINE, GRID } from '../runtime.module.js';
import { R, G, P0, pressureAt, heightOf, temperatureAt } from '../physics.module.js';

const KM = 1000, WIDTH = 200 * KM, START = 15 * KM, END = 185 * KM, SURFACES = [1000, 850, 500, 300], SPEED = 3600, SUBSTEP = 60, LONGEST = 24 * 3600;
const WINDOW = 100, BELOW = 10, RADIUS = 6, LIFT = 14, FORCE = 0.22, PUSH_SOFT = 1.8e-3, VELOCITY = 0.3, SPEED_SOFT = 20;
const HALO = 'rgba(20, 20, 22, 0.85)', LEVEL = 'rgba(255,255,255,0.55)', CHOSEN = 'rgba(255,255,255,0.8)', OTHER = 'rgba(255,255,255,0.22)';
const comma = (n) => Math.round(n).toLocaleString('en-US');
const meters = (n) => `${n.toLocaleString('en-US', { minimumFractionDigits: 1, maximumFractionDigits: 1 })} m`;
const round2 = (n) => { const k = 10 ** Math.max(0, Math.floor(Math.log10(n)) - 1); return (Math.round(n / k) * k).toLocaleString('en-US'); };
const clock = (t) => { const minutes = Math.round(t / 60); return `${Math.floor(minutes / 60)} h ${String(minutes % 60).padStart(2, '0')} min`; };

export function mountPush(root) {
  const controls = root.querySelector('.controls');
  let view = 'height', level = 0, difference = 2, x = START, u = 0, t = 0, flying = false, launched = false, done = false, geometry = null;
  const hours = [];
  const surface = () => SURFACES[level] * 100;
  const levelHeight = () => heightOf(surface());
  const highSide = () => P0 + difference * 100;
  const surfaceHeight = (m) => heightOf(P0 * surface() / (highSide() - difference * 100 * m / WIDTH));

  function physics() {
    const z = levelHeight(), dp = difference * 100 * pressureAt(z) / P0, drop = surfaceHeight(0) - z, push = G * drop / WIDTH;
    return { z, dp, drop, push, density: push > 0 ? dp / (push * WIDTH) : pressureAt(z) / (R * temperatureAt(z)) };
  }

  const fig = new Figure(root, { height: 380, minHeight: 300, step, draw });
  const BASE = 'Two columns of air 200 km apart whose temperature falls with height as the average atmosphere’s does, the same on both sides at every height. At the ground the low side has 1000 hPa and the high side more by the difference you set. The parcel feels only the push of the pressure: no friction, no turning planet. Time runs an hour per second.';
  const captionFor = () => view === 'height' ? BASE : `${BASE} The height axis spans only ${WINDOW} m; drawn to the slab’s true proportions, the tilt would be flat to the pixel.`;
  const note = caption(root, captionFor());
  const columns = ['ramp', 'the two columns of air', 'rgba(255,255,255,0.3)', 'rgba(255,255,255,0.05)'];
  const legendFor = () => view === 'height'
    ? [['force', 'the push of the pressure', 'warm'], ['arrow', 'the parcel’s velocity', INK], ['line', 'its path, with a dot every hour', ACCENT], ['dash', 'heights where the pressures are compared, the chosen one solid', LEVEL], columns]
    : [['force', 'the push: g times the downhill slope', 'warm'], ['arrow', 'the parcel’s velocity', INK], ['line', 'its path, with a dot every hour', ACCENT], ['faint', 'the pressure surface', INK], ['dash', 'a level line, for comparison', LEVEL], columns];
  const key = legend(root, legendFor());
  choice(controls, { label: 'View', options: [['at one height', 'height'], ['on a pressure surface', 'surface']], value: view, onChange: (v) => { view = v; note.textContent = captionFor(); key.set(legendFor()); update(); } });
  slider(controls, { label: 'Pressure difference at the ground', min: 0, max: 10, step: 0.5, value: difference, format: (v) => `${v.toFixed(1)} hPa`, onInput: (v) => { difference = v; reset(); update(); } });
  const heights = choice(controls, { label: 'Height', options: SURFACES.map((p, k) => [k === 0 ? 'the ground' : `${(heightOf(p * 100) / KM).toFixed(1)} km`, String(k)]), value: '0', onChange: (v) => { if (Number(v) === level) return; level = Number(v); reset(); update(); } });
  const [launch] = buttons(controls, [['Launch', () => { reset(); flying = launched = true; launch.textContent = 'Launch again'; fig.play(true); update(); }]]);
  const info = readout(controls), live = readout(controls);

  function reset() { x = START; u = 0; t = 0; hours.length = 0; flying = launched = done = false; launch.textContent = 'Launch'; fig.play(false); }

  function update() { report(); fig.render(); }

  function report() {
    const { z, dp, density, drop, push } = physics(), where = z < 1 ? 'at the ground' : `at ${(z / KM).toFixed(1)} km`;
    const pushText = push > 0 ? `${(push * 1e3).toFixed(2)} × 10⁻³ m/s²` : 'none', hourText = `${(push * 3600).toFixed(1)} m/s`;
    if (view === 'height') info.set([[`pressure difference ${where}`, `${(dp / 100).toFixed(2)} hPa over 200 km`], ['air density there', `${density.toFixed(2)} kg/m³`], ['push per kilogram', pushText], ['speed gained per hour', hourText]]);
    else info.set([[`the ${SURFACES[level]} hPa surface`, drop > 0.05 ? `${drop.toFixed(1)} m lower on the low side` : 'level'], ['slope', drop > 0.05 ? `1 in ${round2(WIDTH / drop)}` : 'none'], ['push per kilogram, g times the slope', pushText], ['speed gained per hour', hourText]]);
    live.set([['time', clock(t)], ['speed now', `${u.toFixed(1)} m/s`], ['distance', `${((x - START) / KM).toFixed(0)} km`]]);
  }

  function step(dt) {
    if (!flying) { fig.play(false); return false; }
    const a = physics().push;
    for (let remaining = dt * SPEED; remaining > 0 && flying;) {
      const h = Math.min(SUBSTEP, remaining), hour = Math.floor(t / 3600);
      remaining -= h;
      const next = x + u * h + 0.5 * a * h * h;
      if (next >= END) {
        const tau = a > 0 ? (Math.sqrt(u * u + 2 * a * (END - x)) - u) / a : (END - x) / u;
        x = END; t += tau; u += a * tau; flying = false; done = true;
      } else { x = next; u += a * h; t += h; }
      if (Math.floor(t / 3600) > hour) hours.push(x);
      if (t >= LONGEST) flying = false;
    }
    report();
  }

  function draw(ctx, w, h) {
    const narrow = w < 480, size = narrow ? 10 : 11, left = 46, right = 12, top = 30, bottom = h - 30, plotW = w - left - right, edge = left + plotW - 3;
    const z = levelHeight(), zMin = view === 'height' ? 0 : Math.max(0, z - BELOW), zMax = view === 'height' ? 10000 : zMin + WINDOW, base = zMin === 0 ? bottom - RADIUS - 1 : bottom;
    const X = (m) => left + m / WIDTH * plotW, Y = (zz) => base - (zz - zMin) / (zMax - zMin) * (base - top);
    const colW = Math.max(14, plotW * 0.05), inner = [X(0) + colW, X(WIDTH) - colW];
    const at = (m) => Y(view === 'height' ? z : surfaceHeight(m));
    const a = physics().push;
    const arrows = (m, speed) => {
      const px = X(m);
      return { px, py: at(m), force: Math.min(FORCE * plotW * Math.tanh(a / PUSH_SOFT), edge - px), velocity: Math.min(VELOCITY * plotW * Math.tanh(speed / SPEED_SOFT), edge - px - RADIUS) };
    };
    const band = (x1, x2) => {
      let y = Infinity;
      const n = Math.ceil(X(END) - X(START)), over = (from, to) => from <= x2 && to >= x1;
      for (let s = 0; s <= n; s++) {
        const m = START + s / n * (END - START), { px, py, force, velocity } = arrows(m, Math.sqrt(2 * a * (m - START)));
        if (over(px - RADIUS, px + RADIUS + Math.max(0, velocity))) y = Math.min(y, py - RADIUS);
        if (force >= 4 && over(px, px + force)) y = Math.min(y, py - LIFT - 5);
      }
      return y;
    };

    ctx.strokeStyle = GRID; ctx.lineWidth = 1; ctx.setLineDash([2, 4]);
    const tick = view === 'height' ? 2000 : 20, ticks = [];
    for (let zz = Math.ceil(zMin / tick) * tick; zz <= zMax; zz += tick) ticks.push(zz);
    for (const zz of ticks) { if (zz > zMin) { ctx.beginPath(); ctx.moveTo(left, Y(zz)); ctx.lineTo(left + plotW, Y(zz)); ctx.stroke(); } }
    ctx.setLineDash([]);
    ticks.forEach((zz, i) => {
      const last = i === ticks.length - 1, label = view === 'height' ? (last ? `${zz / KM} km` : `${zz / KM}`) : (last ? `${comma(zz)} m` : comma(zz));
      text(ctx, label, left - 6, Y(zz), { align: 'right', color: MUTED, size: 10 });
    });
    for (const d of [0, 50, 100, 150, 200]) text(ctx, d === 200 ? '200 km' : `${d}`, X(d * KM), bottom + 18, { align: d === 0 ? 'left' : d === 200 ? 'right' : 'center', color: MUTED, size: 10 });
    if (zMin === 0) { ctx.fillStyle = '#5a4a36'; ctx.fillRect(left, bottom, plotW, 6); }

    const fade = ctx.createLinearGradient(0, bottom, 0, top), thin = pressureAt(zMax) / P0, thick = pressureAt(zMin) / P0;
    fade.addColorStop(0, `rgba(255,255,255,${(0.18 * thick).toFixed(3)})`); fade.addColorStop(1, `rgba(255,255,255,${(0.18 * thin).toFixed(3)})`);
    ctx.fillStyle = fade; ctx.fillRect(X(0), top, colW, bottom - top); ctx.fillRect(inner[1], top, colW, bottom - top);
    ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(inner[0], top); ctx.lineTo(inner[0], bottom); ctx.moveTo(inner[1], top); ctx.lineTo(inner[1], bottom); ctx.stroke();
    if (difference > 0) {
      text(ctx, narrow ? 'higher' : 'higher pressure', X(0), top - 12, { color: INK, size });
      text(ctx, narrow ? 'lower' : 'lower pressure', X(WIDTH), top - 12, { align: 'right', color: INK, size });
    }
    if (!launched && view === 'height') text(ctx, 'tap a dashed line to choose a height', (inner[0] + inner[1]) / 2, top - 12, { align: 'center', color: MUTED, size: 10 });
    if (launched) text(ctx, `${done ? 'arrived after ' : ''}${clock(t)}${done ? ' at ' : ' · '}${u.toFixed(1)} m/s`, (inner[0] + inner[1]) / 2, top - 12, { align: 'center', color: done ? INK : MUTED, size: 10 });

    for (const k of view === 'height' ? SURFACES.keys() : [level]) {
      const chosen = view === 'height' && k === level, y = Y(heightOf(SURFACES[k] * 100));
      ctx.strokeStyle = chosen ? CHOSEN : view === 'height' ? OTHER : LEVEL; ctx.lineWidth = chosen ? 1.5 : 1; ctx.setLineDash(chosen ? [] : [3, 4]);
      ctx.beginPath(); ctx.moveTo(inner[0], y); ctx.lineTo(inner[1], y); ctx.stroke();
    }
    ctx.setLineDash([]);
    if (view === 'surface') {
      ctx.strokeStyle = INK; ctx.lineWidth = 1.5; ctx.beginPath();
      for (let s = 0; s <= 40; s++) { const m = s / 40 * WIDTH, y = Y(surfaceHeight(m)); if (s) ctx.lineTo(X(m), y); else ctx.moveTo(X(m), y); }
      ctx.stroke();
    }
    const line = (px) => at(clamp((px - left) / plotW, 0, 1) * WIDTH);
    const obstacles = view === 'height' ? SURFACES.map((p, k) => (k === level ? null : Y(heightOf(p * 100)))).filter((v) => v !== null) : [Y(z)];
    const tag = (str, x0, align, under = true, avoid = null) => {
      ctx.font = '400 10px system-ui, sans-serif';
      const span = ctx.measureText(str).width, x1 = align === 'right' ? x0 - span : align === 'center' ? x0 - span / 2 : x0, x2 = x1 + span;
      const below = Math.max(line(x1), line(x2)) + 14, down = under && below + 6 < bottom - 2;
      let y = down ? below : Math.min(band(x1, x2), line(x1), line(x2)) - 10;
      for (let i = 0; i < 6 && obstacles.some((oy) => Math.abs(oy - y) < 9); i++) y += down ? 12 : -12;
      if (avoid !== null && Math.abs(y - avoid) < 12) y = avoid - 12;
      text(ctx, str, x0, Math.max(top + 6, y), { align, color: INK, size: 10, halo: HALO });
      return y;
    };
    const high = tag(view === 'height' ? `${(highSide() * pressureAt(z) / P0 / 100).toFixed(1)} hPa` : meters(surfaceHeight(0)), inner[0] + 5, 'left');
    if (view === 'height') tag(`${(pressureAt(z) / 100).toFixed(1)} hPa`, inner[1] - 5, 'right');
    else {
      if (z >= 1) tag(meters(z), inner[1] - 5, 'right');
      tag(`the ${SURFACES[level]} hPa surface`, inner[0] + 5, 'left', false, high);
    }
    geometry = view === 'height' ? { inner, ys: SURFACES.map((p) => Y(heightOf(p * 100))) } : null;

    ctx.save();
    ctx.beginPath(); ctx.rect(left, top - 4, plotW, bottom - top + 4); ctx.clip();
    if (launched) {
      ctx.strokeStyle = ACCENT; ctx.lineWidth = 2; ctx.beginPath();
      for (let s = 0; s <= 40; s++) { const m = START + s / 40 * (x - START); if (s) ctx.lineTo(X(m), at(m)); else ctx.moveTo(X(m), at(m)); }
      ctx.stroke();
      ctx.fillStyle = ACCENT;
      for (const m of hours) { ctx.beginPath(); ctx.arc(X(m), at(m), 3, 0, Math.PI * 2); ctx.fill(); }
    }
    const { px, py, force, velocity } = arrows(x, u);
    if (!done) {
      if (force >= 4) arrow(ctx, px, py - LIFT, px + force, py - LIFT, { color: anomalyColor(1), width: 2, head: clamp(force * 0.6, 4, 8), dash: [5, 4], open: true });
      if (velocity >= 4) arrow(ctx, px + RADIUS, py, px + RADIUS + velocity, py, { color: INK, width: 2, head: clamp(velocity * 0.6, 4, 8) });
    }
    ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(px, py, RADIUS, 0, Math.PI * 2); ctx.fill();
    ctx.restore();
  }

  function pick({ x: px, y: py }) {
    if (!geometry || px < geometry.inner[0] || px > geometry.inner[1]) return null;
    let best = null, gap = 14;
    geometry.ys.forEach((y, k) => { if (Math.abs(y - py) < gap) { gap = Math.abs(y - py); best = k; } });
    return best;
  }
  fig.pointer({ hit: (at) => pick(at) !== null, down: (at) => { const k = pick(at); if (k === null || k === level) return; level = k; heights.value = String(k); reset(); update(); } });

  report();
}
