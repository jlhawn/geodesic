import { Figure, slider, choice, buttons, legend, caption, readout, text, arrow, sequentialRGB, termColor, paletteVersion, clamp, INK, MUTED, GRID, LINE, ACCENT } from '../runtime.module.js';

const KM = 1000, HOUR = 3600, WIDE = 1500 * KM, TALL = 900 * KM, FRONT = 300 * KM, SIGMA = 150 * KM, RUN = 12 * HOUR, SPEED = 1.5 * HOUR, MAX = 30, SPAN = 10, BRIGHT = 0.8, HALO = 'rgba(20, 20, 22, 0.85)';
const NAMES = ['north', 'north-northeast', 'northeast', 'east-northeast', 'east', 'east-southeast', 'southeast', 'south-southeast', 'south', 'south-southwest', 'southwest', 'west-southwest', 'west', 'west-northwest', 'northwest', 'north-northwest'];
const bump = (x, y) => 6 * Math.exp(-((x - WIDE / 2) ** 2 + (y - TALL / 2) ** 2) / (2 * SIGMA * SIGMA));
const FIELDS = {
  front: { theta: (x, y) => 285 - SPAN * clamp((y - TALL / 2) / FRONT, -0.5, 0.5), slope: (x, y) => [0, Math.abs(y - TALL / 2) < FRONT / 2 ? -SPAN / FRONT : 0], levels: [282, 284, 286, 288], apart: '2 K', probe: [500 * KM, 420 * KM], range: [280, 290], ticks: [280, 285, 290] },
  blob: { theta: (x, y) => 284 + bump(x, y), slope: (x, y) => { const b = bump(x, y) / (SIGMA * SIGMA); return [-b * (x - WIDE / 2), -b * (y - TALL / 2)]; }, levels: [285, 286, 287, 288, 289], apart: '1 K', probe: [1000 * KM, 700 * KM], range: [284, 290], ticks: [284, 287, 290] },
};
const direction = (d) => { const k = Math.round(d / 22.5); return `${Math.abs(d - k * 22.5) < 0.25 ? '' : 'about '}the ${NAMES[k % 16]}`; };
const layout = (w) => { const pad = 12, mw = w - 2 * pad, mh = mw * 0.6; return { pad, mw, mh, top: 22, h: Math.round(22 + mh + (w < 500 ? 110 : 128)) }; };

export function mountCarry(root) {
  const controls = root.querySelector('.controls'), blue = termColor('b'), field = document.createElement('canvas'), paint = field.getContext('2d');
  let preset = 'front', speed = 15, from = 225, probe = [...FIELDS.front.probe], t = 0, drag = null, grab = [0, 0], geometry = null, painted = '', shades = null, shaded = -1, keyed = paletteVersion;
  const wind = () => { const a = from * Math.PI / 180; return [-speed * Math.sin(a), -speed * Math.cos(a)]; };
  const thetaAt = (x, y, s = t) => { const [u, v] = wind(); return FIELDS[preset].theta(x - u * s, y - v * s); };
  const slopeAt = (x, y, s = t) => { const [u, v] = wind(); return FIELDS[preset].slope(x - u * s, y - v * s); };
  const rateAt = (x, y, s = t) => { const [u, v] = wind(), [gx, gy] = slopeAt(x, y, s); return -(u * gx + v * gy) * HOUR; };
  const parcel = () => { const [u, v] = wind(); return [probe[0] + u * t, probe[1] + v * t]; };
  const rgb = (f) => `rgb(${sequentialRGB(f).map((c) => Math.round(c * 255)).join(', ')})`;
  const items = () => [['ramp', `θ, from ${FIELDS[preset].range[0]} K to ${FIELDS[preset].range[1]} K`, rgb(0), rgb(BRIGHT), rgb(BRIGHT / 2)], ['faint', `lines of equal θ, ${FIELDS[preset].apart} apart`, LINE], ['arrow', 'the wind, the same everywhere', INK], ['dots', '−u·∇θ: + marks where the wind warms a fixed spot, − marks where it cools it, bigger where faster', blue], ['line', 'θ at the probe, which stays put', INK], ['line', 'θ in the parcel, which moves with the wind', ACCENT], ['faint', 'the slope of the probe’s θ right now, which is −u·∇θ there', blue]];

  const fig = new Figure(root, { height: 520, minHeight: 320, step, draw });
  const fit = fig.resize.bind(fig);
  fig.resize = () => { fig.height = fig.minHeight = layout(fig.stage.clientWidth).h; fit(); };
  caption(root, 'A map 1,500 km wide and 900 km from south to north. The wind is the same everywhere, and so is π, so the term ∇·(πθu) comes down to π u·∇θ. Divided by π and moved to the other side of the equation, it changes θ at each fixed spot at the rate −u·∇θ: the wind bringing in air of a different θ. The blue marks show that rate, + where it warms the spot and − where it cools it. The white dot is the air that sat over the probe when the clock started. Time runs 1.5 hours a second.');
  const key = legend(root, items());
  choice(controls, { label: 'Start with', options: [['a front', 'front'], ['a warm blob', 'blob']], value: preset, onChange: (v) => { preset = v; probe = [...FIELDS[v].probe]; key.set(items()); restart(); }, span: true });
  const speedControl = slider(controls, { label: 'Wind speed', min: 0, max: MAX, step: 1, value: speed, format: (v) => `${v} m/s`, onInput: (v) => { speed = v; restart(); } });
  const fromControl = slider(controls, { label: 'Wind from', min: 0, max: 360, step: 0.5, value: from, format: direction, onInput: (v) => { from = v % 360; restart(); } });
  const [play] = buttons(controls, [['Play', toggle], ['Reset', () => { fig.play(false); restart(); }]]);
  const out = readout(controls), box = controls.querySelector('.readout');

  function label() { play.textContent = fig.running ? 'Pause' : t >= RUN ? 'Play again' : 'Play'; }
  function restart() { t = 0; label(); report(); fig.render(); }
  function toggle() { if (fig.running) fig.play(false); else { if (t >= RUN) t = 0; fig.play(true); } label(); report(); fig.render(); }

  function report(every = 1) {
    const [u, v] = wind(), [x, y] = probe, [gx, gy] = slopeAt(x, y), rate = rateAt(x, y), up = speed ? -(u * gx + v * gy) / speed * 100 * KM : 0;
    const [qx, qy] = parcel(), gone = qx < 0 || qx > WIDE || qy < 0 || qy > TALL;
    out.set([
      ['wind', speed ? `${speed} m/s from ${direction(from)}, ${+(speed * 3.6).toFixed(1)} km an hour` : 'calm'],
      ['θ at the probe', `${thetaAt(x, y).toFixed(1)} K`],
      speed ? ['upwind, the air is', Math.abs(up) < 0.005 ? 'just as warm' : `${Math.abs(up).toFixed(2)} K ${up > 0 ? 'warmer' : 'colder'} every 100 km`] : ['with no wind,', 'no other air arrives'],
      ['so θ at the probe is', Math.abs(rate) < 0.005 ? 'holding steady' : `${rate > 0 ? 'rising' : 'falling'} ${Math.abs(rate).toFixed(2)} K an hour`],
      ['the parcel’s θ', `${thetaAt(x, y, 0).toFixed(1)} K, ${t === 0 ? 'the same as the probe’s for now' : gone ? 'unchanged, though it has left the map' : 'unchanged'}`],
    ], every);
    const b = box.querySelectorAll('b')[3];
    if (b) b.style.color = blue;
  }

  function step(dt) {
    if (t >= RUN) { fig.play(false); label(); return false; }
    t = Math.min(RUN, t + dt * SPEED);
    if (t >= RUN) { fig.play(false); label(); }
    report(t >= RUN ? 1 : 3);
  }

  function shade(cw, ch, ox, oy) {
    const id = `${preset}/${ox}/${oy}/${cw}x${ch}/${paletteVersion}`;
    if (id === painted) return;
    painted = id;
    if (shaded !== paletteVersion) { shades = Array.from({ length: 256 }, (_, i) => sequentialRGB(i / 255).map((c) => Math.round(c * 255))); shaded = paletteVersion; }
    if (field.width !== cw || field.height !== ch) { field.width = cw; field.height = ch; }
    const img = paint.createImageData(cw, ch), data = img.data, { theta, range: [lo, hi] } = FIELDS[preset], scale = BRIGHT * 255 / (hi - lo);
    for (let j = 0; j < ch; j++) {
      const y = (1 - (j + 0.5) / ch) * TALL - oy;
      for (let i = 0; i < cw; i++) { const c = shades[clamp(Math.round((theta((i + 0.5) / cw * WIDE - ox, y) - lo) * scale), 0, 255)], n = 4 * (j * cw + i); data[n] = c[0]; data[n + 1] = c[1]; data[n + 2] = c[2]; data[n + 3] = 255; }
    }
    paint.putImageData(img, 0, 0);
  }

  function draw(ctx, w, h) {
    if (keyed !== paletteVersion) { keyed = paletteVersion; key.set(items()); }
    const { pad, mw, mh, top } = layout(w), x0 = pad, y0 = top, k = mw / WIDE, bottom = y0 + mh;
    const sx = (x) => x0 + x * k, sy = (y) => bottom - y * k, [u, v] = wind(), ox = u * t, oy = v * t;
    const R = clamp(w * 0.085, 28, 50), dial = { x: x0 + mw - R - 14, y: bottom - R - 14, r: R };
    geometry = { x0, y0, mw, mh, k, sx, sy, dial };
    shade(Math.ceil(mw / 2), Math.ceil(mh / 2), ox, oy);
    ctx.save();
    ctx.beginPath(); ctx.rect(x0, y0, mw, mh); ctx.clip();
    ctx.imageSmoothingEnabled = true; ctx.drawImage(field, x0, y0, mw, mh);
    ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.beginPath();
    if (preset === 'front') for (const level of FIELDS.front.levels) { const y = sy(TALL / 2 + FRONT * (285 - level) / SPAN + oy); ctx.moveTo(x0, y); ctx.lineTo(x0 + mw, y); }
    else for (const level of FIELDS.blob.levels) { const r = SIGMA * Math.sqrt(2 * Math.log(6 / (level - 284))) * k, cx = sx(WIDE / 2 + ox), cy = sy(TALL / 2 + oy); ctx.moveTo(cx + r, cy); ctx.arc(cx, cy, r, 0, 2 * Math.PI); }
    ctx.stroke();
    const bx = x0 + 14, by = bottom - 14, bl = FRONT * k, g = w < 500 ? 19 : 24, nx = Math.floor(mw / g), ny = Math.floor(mh / g), gx0 = x0 + (mw - (nx - 1) * g) / 2, gy0 = y0 + (mh - (ny - 1) * g) / 2;
    ctx.beginPath();
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const gx = gx0 + i * g, gy = gy0 + j * g;
      if (gx > dial.x - R - 12 && gy > dial.y - R - 28 || gx < bx + bl + 8 && gy > by - 26) continue;
      const r = rateAt((gx - x0) / k, (bottom - gy) / k);
      if (Math.abs(r) < 0.05) continue;
      const s = 1.5 + 4.5 * Math.sqrt(Math.min(1, Math.abs(r) / 3.6));
      ctx.moveTo(gx - s, gy); ctx.lineTo(gx + s, gy);
      if (r > 0) { ctx.moveTo(gx, gy - s); ctx.lineTo(gx, gy + s); }
    }
    ctx.lineCap = 'round'; ctx.strokeStyle = 'rgba(20, 20, 22, 0.6)'; ctx.lineWidth = 4.5; ctx.stroke(); ctx.strokeStyle = blue; ctx.lineWidth = 2; ctx.stroke();
    for (const [color, width] of [[HALO, 4], [INK, 1.5]]) { ctx.strokeStyle = color; ctx.lineWidth = width; ctx.beginPath(); ctx.moveTo(bx, by - 4); ctx.lineTo(bx, by); ctx.lineTo(bx + bl, by); ctx.lineTo(bx + bl, by - 4); ctx.stroke(); }
    text(ctx, '300 km', bx + bl / 2, by - 10, { align: 'center', color: INK, size: 10, halo: HALO });
    ctx.fillStyle = 'rgba(20, 20, 22, 0.72)'; ctx.beginPath(); ctx.arc(dial.x, dial.y, R + 6, 0, 2 * Math.PI); ctx.fill();
    ctx.lineWidth = 1;
    for (const s of [10, 20, 30]) { ctx.strokeStyle = s === MAX ? LINE : GRID; ctx.beginPath(); ctx.arc(dial.x, dial.y, s / MAX * R, 0, 2 * Math.PI); ctx.stroke(); }
    arrow(ctx, dial.x, dial.y, dial.x + u / MAX * R, dial.y - v / MAX * R, { color: INK, width: 2, head: 8 });
    ctx.fillStyle = INK; ctx.beginPath(); ctx.arc(dial.x, dial.y, 2, 0, 2 * Math.PI); ctx.fill();
    text(ctx, speed ? `wind ${speed} m/s` : 'calm', dial.x, dial.y - R - 16, { align: 'center', color: INK, size: 11, halo: HALO });
    const px = sx(probe[0]), py = sy(probe[1]), [qx, qy] = parcel(), cx = sx(qx), cy = sy(qy), unit = speed ? [u / speed, v / speed] : [0, -1];
    if (t > 0) { ctx.strokeStyle = ACCENT; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(px, py); ctx.lineTo(cx, cy); ctx.stroke(); }
    ctx.fillStyle = '#fff'; ctx.strokeStyle = HALO; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(cx, cy, 4.5, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
    const tag = (str, x, y, color) => text(ctx, str, clamp(x, x0 + 20, x0 + mw - 20), clamp(y, y0 + 8, bottom - 8), { align: 'center', color, size: 11, halo: HALO });
    if (Math.hypot(cx - px, cy - py) > 30 && qx >= 0 && qx <= WIDE && qy >= 0 && qy <= TALL) tag('parcel', cx - unit[1] * 18, cy - unit[0] * 18, ACCENT);
    for (const [color, width] of [[HALO, 4.5], [INK, 2]]) {
      ctx.strokeStyle = color; ctx.lineWidth = width; ctx.beginPath(); ctx.arc(px, py, 8, 0, 2 * Math.PI);
      for (const [dx, dy] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) { ctx.moveTo(px + dx * 11, py + dy * 11); ctx.lineTo(px + dx * 15, py + dy * 15); }
      ctx.stroke();
    }
    tag('probe', px - unit[0] * 30, py + unit[1] * 30, INK);
    ctx.restore();
    ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.strokeRect(x0 + 0.5, y0 + 0.5, mw - 1, mh - 1);
    text(ctx, t > 0 ? `after ${(t / HOUR).toFixed(1)} hours` : 'drag the probe or the wind arrow', t > 0 ? x0 : w / 2, 11, { align: t > 0 ? 'left' : 'center', color: MUTED, size: 11 });
    text(ctx, 'north ↑', w - pad, 11, { align: 'right', color: MUTED, size: 11 });
    chart(ctx, w, h, pad, bottom);
  }

  function chart(ctx, w, h, pad, under) {
    const left = pad + 34, right = w - pad - 44, top = under + 32, low = h - 20;
    const { range: [lo, hi], ticks } = FIELDS[preset], X = (s) => left + s / RUN * (right - left), Y = (th) => low - (th - lo + 1) / (hi - lo + 2) * (low - top);
    const at = (s) => thetaAt(probe[0], probe[1], s), kept = at(0);
    text(ctx, 'θ over 12 hours', pad, under + 18, { color: MUTED, size: 11 });
    ctx.lineWidth = 1;
    for (const th of ticks) { ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(left, Y(th)); ctx.lineTo(right, Y(th)); ctx.stroke(); text(ctx, `${th} K`, left - 6, Y(th), { align: 'right', color: MUTED, size: 10 }); }
    for (let s = 0; s <= 12; s += 3) text(ctx, s === 12 ? '12 hours' : `${s}`, X(s * HOUR), low + 12, { align: 'center', color: MUTED, size: 10 });
    const line = (fn, a, b, color, width) => {
      const n = Math.max(1, Math.ceil(96 * (b - a) / RUN));
      ctx.strokeStyle = color; ctx.lineWidth = width; ctx.beginPath();
      for (let i = 0; i <= n; i++) { const s = a + (b - a) * i / n; if (i) ctx.lineTo(X(s), Y(fn(s))); else ctx.moveTo(X(s), Y(fn(s))); }
      ctx.stroke();
    };
    line(() => kept, 0, RUN, 'rgba(255, 232, 160, 0.3)', 1.5);
    line(at, 0, RUN, 'rgba(221, 221, 221, 0.3)', 1.5);
    if (t > 0) {
      ctx.strokeStyle = 'rgba(255, 255, 255, 0.2)'; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(X(t), top); ctx.lineTo(X(t), low); ctx.stroke();
      line(() => kept, 0, t, ACCENT, 2);
      line(at, 0, t, INK, 2);
    }
    if (speed) {
      const r = rateAt(probe[0], probe[1]) / HOUR, d = 1.2 * HOUR, th = at(t);
      ctx.save(); ctx.beginPath(); ctx.rect(left - 4, top - 4, right - left + 8, low - top + 8); ctx.clip();
      ctx.strokeStyle = blue; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(X(t - d), Y(th - r * d)); ctx.lineTo(X(t + d), Y(th + r * d)); ctx.stroke();
      ctx.restore();
    }
    for (const [th, color] of [[kept, ACCENT], [at(t), INK]]) { ctx.fillStyle = color; ctx.beginPath(); ctx.arc(X(t), Y(th), 3.5, 0, 2 * Math.PI); ctx.fill(); }
    let a = Y(at(RUN)), b = Y(kept);
    if (Math.abs(a - b) < 12) { const m = (a + b) / 2, s = a <= b ? -1 : 1; a = m + s * 6; b = m - s * 6; }
    text(ctx, 'probe', right + 6, a, { color: INK, size: 10 });
    text(ctx, 'parcel', right + 6, b, { color: ACCENT, size: 10 });
  }

  const pick = (p) => {
    if (!geometry) return null;
    const { sx, sy, dial } = geometry;
    if (Math.hypot(p.x - sx(probe[0]), p.y - sy(probe[1])) < 22) return 'probe';
    if (Math.hypot(p.x - dial.x, p.y - dial.y) < dial.r + 10) return 'dial';
    return null;
  };
  function place(p) {
    const { x0, y0, mh, k } = geometry, [u, v] = wind(), peak = [clamp(WIDE / 2 + u * t, 0, WIDE), clamp(TALL / 2 + v * t, 0, TALL)];
    probe = [clamp((p.x + grab[0] - x0) / k, 0, WIDE), clamp((y0 + mh - p.y - grab[1]) / k, 0, TALL)];
    if (preset === 'blob' && Math.hypot(probe[0] - peak[0], probe[1] - peak[1]) * k < 7) probe = peak;
    report(); fig.render();
  }
  function steer(p) {
    const { dial } = geometry, dx = p.x - dial.x, dy = dial.y - p.y;
    speed = clamp(Math.round(Math.hypot(dx, dy) / dial.r * MAX), 0, MAX);
    if (speed) { const a = (Math.atan2(-dx, -dy) * 180 / Math.PI + 360) % 360, near = Math.round(a / 22.5) * 22.5; from = (Math.abs(a - near) < 4 ? near : Math.round(a * 2) / 2) % 360; }
    speedControl.value = speed; fromControl.value = from;
    restart();
  }
  fig.pointer({
    hit: (p) => pick(p) !== null,
    down: (p) => { drag = pick(p); if (drag === 'probe') grab = [geometry.sx(probe[0]) - p.x, geometry.sy(probe[1]) - p.y]; else steer(p); },
    move: (p) => { if (drag === 'probe') place(p); else if (drag === 'dial') steer(p); },
    up: () => { drag = null; },
  });
  fig.canvas.addEventListener('pointermove', (e) => { if (drag) return; const r = fig.canvas.getBoundingClientRect(); fig.canvas.style.cursor = pick({ x: e.clientX - r.left, y: e.clientY - r.top }) ? 'grab' : ''; });

  report();
}
