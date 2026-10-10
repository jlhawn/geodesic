import { Figure, slider, choice, buttons, legend, caption, readout, text, arrow, anomalyColor, ACCENT, INK, MUTED, GRID } from '../runtime.module.js';

const OMEGA = 7.2921e-5, RHO = 1.2, SPEED = 7200, PUSH = 10, NUDGE = 3, EARTH_RADIUS = 6371e3, KM = 1000, CROWD = 12;
const SIZE = 300 * KM, START = 300 * KM, PER_ACCELERATION = 36 / (100 / (RHO * 1e5));

export function mountInertial(root) {
  const controls = root.querySelector('.controls');
  let mode = 'slope', crowd = false, latitude = 45, gradient = 0, center = -8, t = 0;
  let parcels = [];
  const average = [];
  const latitudeAt = (y) => latitude + (y / EARTH_RADIUS) * 180 / Math.PI;
  const f = (y) => 2 * OMEGA * Math.sin(latitudeAt(y) * Math.PI / 180);
  const pressureAt = (py) => 1000 - gradient * py / (100 * KM);
  function force(px, py) {
    if (mode === 'slope') return [0, gradient * 100 / (RHO * 1e5)];
    const k = (center * 100 / (RHO * SIZE * SIZE)) * Math.exp(-(px * px + py * py) / (2 * SIZE * SIZE));
    return [k * px, k * py];
  }
  const ensemble = () => crowd;
  function balance(r) {
    const fc = f(-r), F = Math.hypot(...force(0, -r)), disc = fc * fc * r * r / 4 - r * F;
    if (center < 0) return -fc * r / 2 + Math.sqrt(fc * fc * r * r / 4 + r * F);
    return disc >= 0 ? fc * r / 2 - Math.sqrt(disc) : null;
  }
  function mean() {
    const n = parcels.length, m = { x: 0, y: 0, u: 0, v: 0 };
    for (const p of parcels) { m.x += p.x / n; m.y += p.y / n; m.u += p.u / n; m.v += p.v / n; }
    return m;
  }

  const fig = new Figure(root, { height: 420, step, draw });
  const note = caption(root, '');
  const key = legend(root, []);
  choice(controls, { label: 'Pressure', options: [['slopes evenly', 'slope'], ['circles a low or a high', 'radial']], value: mode, onChange: (val) => { mode = val; reveal(); restart(); }, span: true });
  slider(controls, { label: 'Latitude', min: 5, max: 90, step: 1, value: latitude, format: (val) => `${val}° north`, onInput: (val) => { latitude = val; restart(); } });
  choice(controls, { label: 'Parcels', options: [['one', 'one'], [`${CROWD}, shoved in different directions`, 'crowd']], value: 'one', onChange: (val) => { crowd = val === 'crowd'; restart(); }, span: true });
  const slopeBox = document.createElement('div'), radialBox = document.createElement('div');
  for (const box of [slopeBox, radialBox]) { box.style.gridColumn = '1 / -1'; box.style.gap = '10px'; controls.append(box); }
  const reveal = () => { slopeBox.style.display = mode === 'slope' ? 'grid' : 'none'; radialBox.style.display = mode === 'radial' ? 'grid' : 'none'; };
  reveal();
  slider(slopeBox, { label: 'Pressure rising southward by', min: 0, max: 3, step: 0.1, value: gradient, format: (val) => `${val.toFixed(1)} hPa per 100 km`, onInput: (val) => { gradient = val; restart(); } });
  slider(radialBox, { label: 'Pressure at the center', min: -10, max: 10, step: 1, value: center, format: (val) => val < 0 ? `${-val} hPa lower: a low` : val > 0 ? `${val} hPa higher: a high` : 'the same as around it', onInput: (val) => { center = val; restart(); } });
  buttons(controls, [['Restart', restart]]);
  const out = readout(controls);

  function restart() {
    t = 0; average.length = 0;
    if (mode === 'radial' && !crowd) parcels = [{ x: 0, y: -START, u: 0, v: 0 }];
    else if (mode === 'radial') {
      const speed = balance(START) ?? f(-START) * START / 2, around = center < 0 ? speed : -speed;
      parcels = Array.from({ length: CROWD }, (_, k) => { const a = 2 * Math.PI * k / CROWD; return { x: 0, y: -START, u: around + NUDGE * Math.cos(a), v: NUDGE * Math.sin(a) }; });
    }
    else if (!crowd) parcels = [{ x: 0, y: -PUSH / f(0), u: 0, v: PUSH }];
    else {
      const ug = force(0, 0)[1] / f(0);
      parcels = Array.from({ length: CROWD }, (_, k) => { const a = 2 * Math.PI * k / CROWD; return { x: 0, y: 0, u: ug + PUSH * Math.cos(a), v: PUSH * Math.sin(a) }; });
    }
    for (const p of parcels) p.trace = [];
    note.textContent = mode === 'radial'
      ? crowd
        ? `${CROWD} parcels start together ${START / KM} km south of the center, each moving with the balanced wind there plus ${NUDGE} m/s in a different direction, and the Coriolis force follows each one\u2019s latitude. Time runs ${SPEED / 3600} hours per second.`
        : `The parcel starts at rest ${START / KM} km south of the center, and the Coriolis force follows its latitude as it moves. Time runs ${SPEED / 3600} hours per second.`
      : crowd
        ? `${CROWD} parcels start together on the dashed line, each moving at ${PUSH} m/s relative to the drift but in a different direction, and the Coriolis force follows each one’s latitude. Time runs ${SPEED / 3600} hours per second, and the view follows their average.`
        : `The parcel starts out moving north at ${PUSH} m/s, one loop’s radius south of the dashed line, and the Coriolis force follows its latitude as it moves. Time runs ${SPEED / 3600} hours per second, and the view follows the parcel.`;
    const lines = [
      ['dash', mode === 'slope' ? 'the line of equal pressure the loops circle around' : crowd ? 'the circle of equal pressure they started on' : 'the circle of equal pressure the parcel started on', 'rgba(255,255,255,0.55)'],
      ['faint', mode === 'slope' ? 'other lines of equal pressure, and 100 km marks' : 'circles of equal pressure every hPa, and 100 km marks', 'rgba(255,255,255,0.3)'],
    ];
    key.set(ensemble()
      ? [['faint', 'each parcel’s path', 'rgba(255,232,160,0.5)'], ['line', 'the path of their average', ACCENT], ['arrow', 'their average wind', INK], ...lines]
      : [['line', 'the parcel’s path over the last two days', ACCENT], ['arrow', 'wind', INK], ['force', 'pressure gradient force', 'warm'], ['force', 'Coriolis force', 'cool'], ...lines]);
    show(1);
  }

  function show(every = 6) {
    const m = mean(), fc = f(m.y), [Fx, Fy] = force(m.x, m.y), geostrophic = Math.hypot(Fx, Fy) / fc, p = parcels[0];
    if (mode === 'radial' && crowd) {
      const r = Math.hypot(m.x, m.y), around = r ? (m.v * m.x - m.u * m.y) / r : 0, balanced = center ? balance(START) : 0;
      out.set([['average distance from the center', `${(r / KM).toFixed(0)} km`], ['average wind around the center', `${Math.abs(around).toFixed(1)} m/s ${around >= 0 ? 'counterclockwise' : 'clockwise'}`], ['balanced wind on the starting circle', balanced === null ? 'none possible this close to so strong a high' : `${balanced.toFixed(1)} m/s`]], every);
    } else if (mode === 'radial') out.set([['latitude now', `${latitudeAt(p.y).toFixed(1)}°`], ['distance from the center', `${(Math.hypot(p.x, p.y) / KM).toFixed(0)} km`], ['wind', `${Math.hypot(p.u, p.v).toFixed(1)} m/s`], ['geostrophic wind here', center ? `${geostrophic.toFixed(1)} m/s` : 'none']], every);
    else if (crowd) out.set([['average wind', `${m.u.toFixed(1)} m/s east, ${m.v.toFixed(1)} m/s north`], ['geostrophic drift', gradient > 0 ? `${geostrophic.toFixed(1)} m/s east` : 'none'], ['average drifted so far', `${(m.x / KM).toFixed(0)} km east`]], every);
    else out.set([['latitude now', `${latitudeAt(p.y).toFixed(1)}°`], ['inertial period', `${(2 * Math.PI / fc / 3600).toFixed(1)} h`], ['loop radius', `${(Math.hypot(p.u - Fy / fc, p.v) / fc / KM).toFixed(0)} km`], ['geostrophic drift', gradient > 0 ? `${geostrophic.toFixed(1)} m/s eastward, along the lines` : 'none'], ['drifted so far', `${(p.x / KM).toFixed(0)} km east`]], every);
  }

  function advance(p, s) {
    const fc = f(p.y), [Fx, Fy] = force(p.x, p.y), ug = Fy / fc, vg = -Fx / fc, angle = fc * s, c = Math.cos(angle), sn = Math.sin(angle);
    const ua = p.u - ug, va = p.v - vg;
    p.x += (ug + (ua * sn - va * (c - 1)) / angle) * s; p.y += (vg + (va * sn + ua * (c - 1)) / angle) * s;
    p.u = ug + ua * c + va * sn; p.v = vg + va * c - ua * sn;
  }

  function step(dt) {
    const h = 60;
    for (let remaining = dt * SPEED; remaining > 0; remaining -= h) {
      const s = Math.min(h, remaining);
      for (const p of parcels) { advance(p, s); p.trace.push(t + s, p.x, p.y); }
      t += s;
      const m = mean();
      average.push(t, m.x, m.y);
    }
    for (const path of [average, ...parcels.map((p) => p.trace)]) while (path.length && path[0] < t - 2 * 86400) path.splice(0, 3);
    show();
  }

  function polyline(ctx, path, sx, sy) {
    ctx.beginPath();
    for (let i = 0; i < path.length; i += 3) { if (i) ctx.lineTo(sx(path[i + 1]), sy(path[i + 2])); else ctx.moveTo(sx(path[i + 1]), sy(path[i + 2])); }
    ctx.stroke();
  }

  function draw(ctx, w, h) {
    if (!parcels.length) return;
    const m = mean();
    let reach = 450 * KM;
    if (mode === 'radial') for (const p of parcels) for (let i = 0; i < p.trace.length; i += 3) reach = Math.max(reach, 1.15 * Math.hypot(p.trace[i + 1], p.trace[i + 2]));
    const half = mode === 'slope' ? 260 * KM : reach, scale = h / (2 * half);
    const ox = mode === 'slope' ? m.x : 0, oy = mode === 'slope' ? m.y : 0, cx = w / 2 - ox * scale, cy = h / 2 + oy * scale;
    const sx = (px) => cx + px * scale, sy = (py) => cy - py * scale;
    ctx.strokeStyle = GRID; ctx.lineWidth = 1;
    for (let k = Math.ceil((ox - w / scale / 2) / (100 * KM)); k * 100 * KM < ox + w / scale / 2; k++) { ctx.beginPath(); ctx.moveTo(sx(k * 100 * KM), 0); ctx.lineTo(sx(k * 100 * KM), h); ctx.stroke(); }
    for (let k = Math.floor((oy - half) / (100 * KM)); k * 100 * KM <= oy + half; k++) {
      const start = mode === 'slope' && k === 0;
      ctx.strokeStyle = start ? 'rgba(255,255,255,0.55)' : GRID; ctx.setLineDash(start ? [4, 5] : []);
      ctx.beginPath(); ctx.moveTo(0, sy(k * 100 * KM)); ctx.lineTo(w, sy(k * 100 * KM)); ctx.stroke(); ctx.setLineDash([]);
      if (mode === 'slope' && gradient > 0) text(ctx, `${pressureAt(k * 100 * KM).toFixed(1)} hPa`, w - 8, sy(k * 100 * KM) - 8, { align: 'right', color: MUTED, size: 10 });
    }
    if (mode === 'radial') {
      for (let level = 1; level < Math.abs(center); level++) {
        const r = SIZE * Math.sqrt(2 * Math.log(Math.abs(center) / (Math.abs(center) - level)));
        ctx.strokeStyle = 'rgba(255,255,255,0.3)'; ctx.beginPath(); ctx.arc(sx(0), sy(0), r * scale, 0, Math.PI * 2); ctx.stroke();
        text(ctx, `${(1000 + Math.sign(center) * (Math.abs(center) - level)).toFixed(0)} hPa`, sx(0) + r * scale * 0.71 + 4, sy(0) - r * scale * 0.71, { color: MUTED, size: 10 });
      }
      ctx.strokeStyle = 'rgba(255,255,255,0.55)'; ctx.setLineDash([4, 5]); ctx.beginPath(); ctx.arc(sx(0), sy(0), START * scale, 0, Math.PI * 2); ctx.stroke(); ctx.setLineDash([]);
      if (center) { ctx.fillStyle = anomalyColor(Math.sign(center)); ctx.beginPath(); ctx.arc(sx(0), sy(0), 4, 0, Math.PI * 2); ctx.fill(); text(ctx, center < 0 ? 'L' : 'H', sx(0) + 8, sy(0) - 8, { color: INK, size: 13, weight: 700 }); }
    }
    if (ensemble()) {
      ctx.strokeStyle = 'rgba(255,232,160,0.4)'; ctx.lineWidth = 1.2;
      for (const p of parcels) polyline(ctx, p.trace, sx, sy);
      ctx.strokeStyle = ACCENT; ctx.lineWidth = 3; polyline(ctx, average, sx, sy);
      for (const p of parcels) { ctx.fillStyle = 'rgba(255,255,255,0.85)'; ctx.beginPath(); ctx.arc(sx(p.x), sy(p.y), 3.5, 0, Math.PI * 2); ctx.fill(); }
      const ax = sx(m.x), ay = sy(m.y);
      ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(ax, ay, 7, 0, Math.PI * 2); ctx.stroke();
      arrow(ctx, ax, ay, ax + m.u * 3, ay - m.v * 3, { color: INK, width: 2, head: 7 });
      text(ctx, 'average', ax + 10, ay + 14, { color: MUTED, size: 10 });
    } else {
      const p = parcels[0];
      ctx.strokeStyle = ACCENT; ctx.lineWidth = 2; polyline(ctx, p.trace, sx, sy);
      const px = sx(p.x), py = sy(p.y), fc = f(p.y), [Fx, Fy] = force(p.x, p.y), coriolis = [fc * p.v * PER_ACCELERATION, -fc * p.u * PER_ACCELERATION];
      ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(px, py, 5, 0, Math.PI * 2); ctx.fill();
      arrow(ctx, px, py, px + p.u * 3, py - p.v * 3, { color: INK, width: 1.5, head: 6 });
      if (Math.hypot(...coriolis) > 2) { arrow(ctx, px, py, px + coriolis[0], py - coriolis[1], { color: anomalyColor(-1), width: 2, head: 8, dash: [5, 4], open: true }); text(ctx, 'Coriolis force', px + coriolis[0] * 1.15 + 6, py - coriolis[1] * 1.15, { color: MUTED, size: 10 }); }
      if (Math.hypot(Fx, Fy) > 0) { arrow(ctx, px, py, px + Fx * PER_ACCELERATION, py - Fy * PER_ACCELERATION, { color: anomalyColor(1), width: 2, head: 8, dash: [5, 4], open: true }); text(ctx, 'pressure gradient force', px + Fx * PER_ACCELERATION * 1.15 + 8, py - Fy * PER_ACCELERATION * 1.15 - 8, { color: MUTED, size: 10 }); }
    }
    text(ctx, 'north ↑', w / 2, 12, { align: 'center', color: MUTED, size: 11 });
    text(ctx, 'east →', w - 10, h - 12, { align: 'right', color: MUTED, size: 11 });
    text(ctx, `${(t / 3600).toFixed(0)} h`, 10, 12, { color: MUTED, size: 11 });
  }

  restart();
  fig.play(true);
}
