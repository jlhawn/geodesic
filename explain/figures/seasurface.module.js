import { Figure, slider, choice, buttons, readout, legend, caption, text, arrow, rampRGB, anomalyColor, clamp, INK, MUTED, GRID, LINE, ACCENT, DARK_NEUTRAL } from '../runtime.module.js';
import { Grid } from '../../js/grid.module.js';

const N = 16, RADIUS = 6371e3, DEG = Math.PI / 180, HOUR = 3600, DAY = 86400, WIDE = 540;
const ETA_RANGE = 0.25, SPEED_RANGE = 0.02, ARROW_SPEED = 0.01, BALANCED = 0.95;
const DEPTHS = [0, 200, 1000, 4000], SHARES = [0, 0.45, 0.75, 1], DEPTH_TICKS = [100, 200, 500, 1000, 2000, 4000];
const SURFACE = '#fff', HALO = 'rgba(20, 20, 22, 0.85)', CHIP = 'rgba(30, 32, 36, 0.88)';
const WATER_ALPHA = [0.3, 0, 0.2, 0.13, 0.07];
const BUMP_HEIGHT = 1;

const km = (m) => Math.round(m / 1000).toLocaleString('en-US');
const minus = (s) => s.replace('-', '−');
const signedCm = (m) => { const v = 100 * m; return Math.abs(v) < 0.05 ? '0.0 cm' : `${v < 0 ? '−' : '+'}${Math.abs(v).toFixed(1)} cm`; };
const cms = (v) => (100 * v).toFixed(2);
const clock = (t) => (t < 48 * HOUR ? `${(t / HOUR).toFixed(t < 10 * HOUR && t % HOUR ? 1 : 0)} h` : `day ${(t / DAY).toFixed(1)}`);
const rate = (kind, t) => (kind === 'lens' ? DAY : t < 8 * HOUR ? HOUR : t < DAY ? 4 * HOUR : 12 * HOUR);

function depthShare(z) {
  for (let i = 1; i < DEPTHS.length; i++) if (z <= DEPTHS[i] || i === DEPTHS.length - 1) return SHARES[i - 1] + (SHARES[i] - SHARES[i - 1]) * clamp((z - DEPTHS[i - 1]) / (DEPTHS[i] - DEPTHS[i - 1]), 0, 1);
  return 1;
}

function projector(latitude) {
  const c = Math.cos(latitude * DEG), s = Math.sin(latitude * DEG);
  return (x, y, z) => {
    const r = Math.hypot(x, y, z), along = (x * c + z * s) / r, east = y / r, north = (-x * s + z * c) / r;
    const t = Math.hypot(east, north);
    if (t < 1e-12) return [0, 0];
    const k = Math.acos(clamp(along, -1, 1)) * RADIUS / t;
    return [east * k, north * k];
  };
}

function niceStep(value) {
  for (const s of [0.05, 0.1, 0.2, 0.25, 0.5, 1]) if (value <= s * 1.0001) return s;
  return 1;
}

export function mountSeasurface(root) {
  const controls = root.querySelector('.controls');
  const state = { kind: 'lens', latitude: 45, depth: 60, width: 1.5e6 };
  const polygons = [...new Grid(N)].map((c) => ({ center: c.centerVertex, vertices: c.vertices }));
  let worker = null, run = 0, pending = false, owed = 0, playing = true, timer = 0, failure = null;
  let layout = null, latest = null, waiting = null, patch = null, arrows = null, lastFront = null, requested = state.latitude;

  const fig = new Figure(root, { height: (w) => (w < WIDE ? 660 : 380), step, draw });
  const visible = fig.setVisible.bind(fig);
  fig.setVisible = (on) => { if (on) boot(); visible(on); };
  const words = caption(root, '');
  legend(root, [
    ['line', 'sea surface, drawn far taller than the depths', SURFACE],
    ['line', 'thermocline: the base of the light water', ACCENT],
    ['faint', 'base of each class', LINE],
    ['ramp', `sea surface ${ETA_RANGE * 100} cm below to ${ETA_RANGE * 100} cm above its mean`, 'cool', 'warm', `rgb(${DARK_NEUTRAL.join(', ')})`],
    ['arrow', 'current in the mixed layer', INK],
    ['dots', 'the sea surface at rest', 'rgba(255, 255, 255, 0.5)'],
    ['dash', 'where the slice cuts', 'rgba(255, 255, 255, 0.55)'],
  ]);
  choice(controls, { label: 'Make', options: [['a lens of light water', 'lens'], ['a bump in the sea surface', 'bump']], value: state.kind, onChange: (v) => { state.kind = v; arrange(); restart(); }, span: true });
  slider(controls, { label: 'Latitude', min: 20, max: 70, step: 5, value: state.latitude, format: (v) => `${v}°N`, onInput: (v) => { state.latitude = v; later(); } });
  slider(controls, { label: 'Thermocline pushed down', min: 20, max: 100, step: 10, value: state.depth, format: (v) => `${v} m at the center`, onInput: (v) => { state.depth = v; later(); } });
  const depthRow = controls.lastElementChild;
  slider(controls, { label: 'Width', min: 1000, max: 2000, step: 250, value: state.width / 1000, format: (v) => `${(2 * v).toLocaleString('en-US')} km across`, onInput: (v) => { state.width = v * 1000; later(); } });
  const [toggle] = buttons(controls, [['Pause', () => {
    if (!latest || latest.done) { restart(); return; }
    playing = !playing; toggle.textContent = playing ? 'Pause' : 'Play'; fig.play(playing);
  }], ['Start again', () => restart()]]);
  const timeOut = readout(controls), levelOut = readout(controls), leftOut = readout(controls), flowOut = readout(controls);
  const leftBox = controls.querySelectorAll('.readout')[2];
  arrange();
  fig.play(true);

  function arrange() {
    const lens = state.kind === 'lens';
    depthRow.style.display = lens ? '' : 'none';
    leftBox.style.display = lens ? 'none' : '';
    const across = `about ${km(2 * state.width)} km across`;
    words.textContent = lens
      ? `The model’s own ocean, run live in your browser on the N = ${N} grid with four classes of water, 1023, 1025, 1026.5 and 1026.95 kg/m³, under a mixed layer at 27 °C, over a flat floor 4,000 m deep, with no wind and no heating. At the start, water of the light 1023 class has taken the place of water of the 1026.5 class over a patch ${across}, so the thermocline lies deeper at its center, but the sea surface is flat. The slice runs from west to east through the center, with the depths stretched toward the surface; the sea surface on top of it is drawn far taller still. The circles are the current through the slice in each layer, crossed where it flows north, away from you, and dotted where it flows south, toward you, their area in proportion to its speed up to 2 cm/s, where they stop growing, and smaller in a layer drawn too thin to hold them. The map is the sea surface seen from above. Time runs at up to a day each second.`
      : `The model’s own ocean, run live in your browser on the N = ${N} grid with four classes of water, 1023, 1025, 1026.5 and 1026.95 kg/m³, under a mixed layer at 27 °C, over a flat floor 4,000 m deep, with no wind and no heating. At the start, the deepest class is ${BUMP_HEIGHT} m thicker over a patch ${across}, which raises the sea surface there by up to ${BUMP_HEIGHT} m without making any water lighter or denser, as Part 5 dropped its stone. The slice runs from west to east through the center, with the depths stretched toward the surface; the sea surface on top of it is drawn far taller still. The circles are the current through the slice in each layer, crossed where it flows north, away from you, and dotted where it flows south, toward you, their area in proportion to its speed up to 2 cm/s, where they stop growing, and smaller in a layer drawn too thin to hold them. The map is the sea surface seen from above. Time runs at up to an hour each second for the first 8 hours, then faster.`;
  }

  function later() { clearTimeout(timer); timer = setTimeout(() => { arrange(); restart(); }, 250); }

  function boot() {
    if (worker) return;
    worker = new Worker(new URL('../ocean.worker.js', import.meta.url), { type: 'module' });
    worker.onmessage = receive;
    worker.onerror = (event) => { event.preventDefault(); failure = 'The ocean could not be started in this browser.'; fig.render(); };
    restart();
  }

  function restart() {
    clearTimeout(timer);
    playing = true; toggle.textContent = 'Pause'; lastFront = null;
    fig.play(true);
    if (!worker) return;
    pending = true; owed = 0; failure = null; requested = state.latitude;
    const options = state.kind === 'lens'
      ? { latitude: state.latitude, depth: state.depth, width: state.width }
      : { latitude: state.latitude, width: state.width, height: BUMP_HEIGHT };
    worker.postMessage({ type: 'start', run: ++run, kind: state.kind, options });
  }

  function receive({ data }) {
    if (data.run !== run) return;
    if (data.type === 'layout') { waiting = data; return; }
    pending = false;
    if (data.type === 'error') { failure = `The ocean stopped: ${data.message}`; fig.render(); return; }
    if (waiting) { adopt(waiting); waiting = null; }
    latest = data;
    if (data.kind === 'bump' && data.frontSpeed != null) lastFront = data.frontSpeed;
    report(data);
    if (data.done) { playing = false; toggle.textContent = 'Play'; fig.play(false); }
    fig.render();
  }

  function adopt(next) {
    layout = next;
    const half = Math.max(...Array.from(next.slice.distance, Math.abs)), reach = projector(requested), cells = [];
    polygons.forEach((cell, i) => {
      const [cx, cy] = reach(cell.center.x, cell.center.y, cell.center.z);
      if (Math.abs(cx) > 1.25 * half || Math.abs(cy) > 1.25 * half) return;
      cells.push({ i, x: cx, y: cy, vertices: cell.vertices.map((v) => reach(v.x, v.y, v.z)), basis: localBasis(reach, next.latitude[i], next.longitude[i], cx, cy) });
    });
    patch = { half: half + 0.6 * spacing(next), cells, sliceHalf: half };
    arrows = null;
  }

  function localBasis(reach, lat, lon, cx, cy) {
    const p = lat * DEG, l = lon * DEG, e = 1e-4;
    const at = (dp, dl) => { const q = p + dp, m = l + dl; return reach(Math.cos(q) * Math.cos(m), Math.cos(q) * Math.sin(m), Math.sin(q)); };
    const [ex, ey] = at(0, e / Math.max(1e-6, Math.cos(p))), [nx, ny] = at(e, 0);
    const el = Math.hypot(ex - cx, ey - cy) || 1, nl = Math.hypot(nx - cx, ny - cy) || 1;
    return [(ex - cx) / el, (ey - cy) / el, (nx - cx) / nl, (ny - cy) / nl];
  }

  function spacing(next) {
    const d = next.slice.distance;
    return d.length > 1 ? (d[d.length - 1] - d[0]) / (d.length - 1) : 4.8e5;
  }

  function report(f) {
    const L = f.layerSpeeds.length, g = f.geostrophic;
    const balance = g.correlation >= BALANCED ? `${g.ratio.toFixed(2)} times the geostrophic` : 'not yet in balance';
    if (f.kind === 'lens') {
      timeOut.set([['day', f.day.toFixed(1)]]);
      levelOut.set([['sea surface at the center', signedCm(f.rise)], ['thermocline', `${f.depression.toFixed(1)} m deeper`], ['mirror of the thermocline', `${(100 * f.expected).toFixed(1)} cm`]]);
    } else {
      timeOut.set(f.time < 48 * HOUR ? [['time', clock(f.time)]] : [['day', (f.time / DAY).toFixed(1)]]);
      levelOut.set([['wave’s leading edge', lastFront == null ? 'timing it' : `${lastFront.toFixed(0)} m/s`], ['√(gD)', `${f.gravitySpeed.toFixed(0)} m/s`]]);
      leftOut.set([['sea surface at the center', signedCm(f.centreHeight)], [`still within ${km(2 * f.width)} km`, f.time < 12 * HOUR ? 'the wave is still passing' : f.kept > 0.0005 ? `${(100 * f.kept).toFixed(0)}% of the raised water` : 'none of the raised water']]);
    }
    flowOut.set([['surface current around it', `${cms(f.layerSpeeds[0])} cm/s, ${balance}`], ['deepest class', `${cms(f.layerSpeeds[L - 1])} cm/s`]]);
  }

  function step(dt) {
    if (!worker || !latest || !playing || latest.done) return false;
    owed += dt;
    if (!pending) { pending = true; worker.postMessage({ type: 'advance', seconds: owed * rate(latest.kind, latest.time) }); owed = 0; }
    return false;
  }

  function geometry(w, h) {
    const narrow = w < WIDE, size = narrow ? 10 : 11;
    if (narrow) {
      const left = 48, right = w - 10, top = 30, strip = 46, depth = 200, map = Math.min(w - 24, 250);
      const bottom = top + strip + depth;
      return { narrow, size, section: { left, right, top, strip, surface: top + strip, bottom }, key: { x: left, y: bottom + 44, stacked: false }, map: { x: (w - map) / 2, y: bottom + 96, size: map } };
    }
    const map = Math.min(236, h - 120), mx = w - 14 - map;
    const left = 52, right = mx - 30, top = 30, strip = 56, depth = h - top - strip - 52;
    const bottom = top + strip + depth;
    return { narrow, size, section: { left, right, top, strip, surface: top + strip, bottom }, key: { x: mx, y: top + map + 42, stacked: true }, map: { x: mx, y: top, size: map } };
  }

  function draw(ctx, w, h) {
    const g = geometry(w, h);
    if (failure || !layout || !latest) {
      text(ctx, failure ?? 'Building the sea…', w / 2, h / 2, { align: 'center', color: MUTED, size: 12 });
      return;
    }
    section(ctx, g, latest);
    map(ctx, g, latest);
  }

  function section(ctx, g, f) {
    const { left, right, top, strip, surface, bottom } = g.section, size = g.size;
    const L = layout.layers, d = f.slice.distance, n = d.length, half = patch.sliceHalf + 0.5 * spacing(layout);
    const X = (m) => (left + right) / 2 + (m / half) * (right - left) / 2;
    const Y = (z) => surface + depthShare(z) * (bottom - surface);
    const ref = f.kind === 'bump' ? f.level : 0;
    const scale = f.kind === 'lens' ? niceStep(1.1 * (layout.dense - layout.light) / layout.densities[0] * f.depth) : BUMP_HEIGHT;
    const E = (strip - 8) / scale;
    const etaY = (s) => surface - (f.slice.eta[s] - ref) * E;
    const interfaceY = (s, k) => Y(-f.slice.interfaces[s * L + k]);
    const topOf = (s, k) => (k === 0 ? etaY(s) : interfaceY(s, k - 1));
    const xs = Array.from(d, X);

    text(ctx, 'Slice from west to east', left, 12, { color: INK, size: size + 1, weight: 500 });

    const trace = (yAt, reverse = false) => {
      const order = reverse ? [...xs.keys()].reverse() : [...xs.keys()];
      const points = [];
      if (!reverse) points.push([left, yAt(0)]); else points.push([right, yAt(n - 1)]);
      for (const s of order) points.push([xs[s], yAt(s)]);
      if (!reverse) points.push([right, yAt(n - 1)]); else points.push([left, yAt(0)]);
      return points;
    };
    ctx.save();
    ctx.beginPath(); ctx.rect(left, top, right - left, bottom - top + 1); ctx.clip();
    for (let k = 0; k < L; k++) {
      const upper = trace((s) => topOf(s, k)), lower = trace((s) => interfaceY(s, k), true);
      ctx.beginPath();
      [...upper, ...lower].forEach(([x, y], j) => (j ? ctx.lineTo(x, y) : ctx.moveTo(x, y)));
      ctx.closePath();
      ctx.fillStyle = k === 1 ? anomalyColor(1, 0.32) : `rgba(105, 155, 215, ${WATER_ALPHA[Math.min(k, WATER_ALPHA.length - 1)]})`;
      ctx.fill();
    }
    ctx.setLineDash([1.5, 3.5]); ctx.lineCap = 'round'; ctx.strokeStyle = 'rgba(255, 255, 255, 0.5)'; ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(left, surface); ctx.lineTo(right, surface); ctx.stroke();
    ctx.setLineDash([]); ctx.lineCap = 'butt';
    const line = (points, color, width) => { ctx.strokeStyle = color; ctx.lineWidth = width; ctx.beginPath(); points.forEach(([x, y], j) => (j ? ctx.lineTo(x, y) : ctx.moveTo(x, y))); ctx.stroke(); };
    for (let k = 0; k < L - 1; k++) if (k !== 1) line(trace((s) => interfaceY(s, k)), LINE, 1);
    line(trace((s) => interfaceY(s, 1)), ACCENT, 1.75);
    line(trace(etaY), SURFACE, 2);
    ctx.restore();
    ctx.strokeStyle = GRID; ctx.lineWidth = 1;
    ctx.strokeRect(left + 0.5, top + 0.5, right - left - 1, bottom - top);

    const gap = n > 1 ? (xs[n - 1] - xs[0]) / (n - 1) : 20;
    for (let s = 0; s < n; s++) {
      for (let k = 0; k < L; k++) {
        const y0 = topOf(s, k), y1 = interfaceY(s, k), room = Math.min(0.42 * gap, 0.42 * (y1 - y0), 7);
        const v = f.slice.across[s * L + k], r = room * Math.sqrt(Math.min(1, Math.abs(v) / SPEED_RANGE));
        if (r < 0.9) continue;
        mark(ctx, xs[s], (y0 + y1) / 2, r, v > 0);
      }
    }

    const tickX = left - 5;
    const tick = (y, label) => { ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(left - 3, y); ctx.lineTo(left, y); ctx.stroke(); text(ctx, label, tickX, y, { align: 'right', color: MUTED, size }); };
    tick(surface - scale * E, scale >= 1 ? `+${scale} m` : `+${Math.round(scale * 100)} cm`);
    tick(surface, '0');
    for (const z of DEPTH_TICKS) tick(Y(z), z === 4000 ? '4000 m' : `${z}`);

    const labelX = right - 4;
    for (let k = 0; k < L; k++) {
      const y0 = topOf(n - 1, k), y1 = interfaceY(n - 1, k);
      if (y1 - y0 < size + 2) continue;
      const name = k === 0 ? 'mixed layer' : k === 1 ? `${layout.densities[k]} kg/m³` : `${layout.densities[k]}`;
      ctx.font = `400 ${size - 1}px system-ui, sans-serif`;
      const width = ctx.measureText(name).width + 8, y = (y0 + y1) / 2;
      ctx.fillStyle = CHIP;
      ctx.beginPath(); ctx.roundRect(labelX - width + 4, y - size / 2 - 2, width, size + 4, 3); ctx.fill();
      text(ctx, name, labelX, y, { align: 'right', color: k === 1 ? ACCENT : MUTED, size: size - 1 });
    }

    if (f.kind === 'lens') {
      const c = xs[Math.floor(n / 2)], tip = etaY(Math.floor(n / 2)), deep = interfaceY(Math.floor(n / 2), 1);
      text(ctx, signedCm(f.rise), c + 10, Math.max(top + 7, tip - 9), { color: SURFACE, size, halo: HALO });
      text(ctx, `${minus((-f.depression).toFixed(0))} m`, c + 10, deep + 10, { color: ACCENT, size, halo: HALO });
    }

    const step = half > 5.5e6 ? 3e6 : 2e6;
    for (let m = -Math.floor(half / step) * step; m <= half; m += step) {
      const x = X(m);
      if (x < left + 6 || x > right - 6) continue;
      ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(x, bottom); ctx.lineTo(x, bottom + 3); ctx.stroke();
      text(ctx, minus(`${Math.round(m / 1000)}`), x, bottom + 11, { align: 'center', color: MUTED, size });
    }
    text(ctx, 'km east of the center', (left + right) / 2, bottom + 25, { align: 'center', color: MUTED, size });

    const { x: kx, y: ky, stacked } = g.key;
    ctx.font = `400 ${size}px system-ui, sans-serif`;
    const second = stacked ? [kx, ky + 20] : [kx + 30 + ctx.measureText('north, away from you').width, ky];
    mark(ctx, kx + 5, ky, 5, true);
    text(ctx, 'north, away from you', kx + 14, ky, { color: MUTED, size });
    mark(ctx, second[0] + 5, second[1], 5, false);
    text(ctx, 'south, toward you', second[0] + 14, second[1], { color: MUTED, size });
  }

  function mark(ctx, x, y, r, north) {
    ctx.save();
    ctx.strokeStyle = INK; ctx.fillStyle = INK; ctx.lineWidth = r > 3 ? 1.2 : 1;
    ctx.beginPath(); ctx.arc(x, y, r, 0, 2 * Math.PI); ctx.stroke();
    if (north) {
      const q = r * 0.62;
      ctx.beginPath(); ctx.moveTo(x - q, y - q); ctx.lineTo(x + q, y + q); ctx.moveTo(x + q, y - q); ctx.lineTo(x - q, y + q); ctx.stroke();
    } else {
      ctx.beginPath(); ctx.arc(x, y, Math.max(0.9, r * 0.28), 0, 2 * Math.PI); ctx.fill();
    }
    ctx.restore();
  }

  function map(ctx, g, f) {
    const { x: mx, y: my, size: S } = g.map, size = g.size;
    const k = S / (2 * patch.half), cx = mx + S / 2, cy = my + S / 2;
    const P = ([x, y]) => [cx + x * k, cy - y * k];
    const ref = f.kind === 'bump' ? f.level : 0, rgb = [0, 0, 0];
    text(ctx, 'Seen from above, north up', mx, my - 18, { color: INK, size: size + 1, weight: 500 });
    ctx.save();
    ctx.beginPath(); ctx.rect(mx, my, S, S); ctx.clip();
    ctx.fillStyle = 'rgba(255,255,255,0.03)'; ctx.fillRect(mx, my, S, S);
    for (const cell of patch.cells) {
      rampRGB(0.5 + clamp((f.eta[cell.i] - ref) / ETA_RANGE, -1, 1) / 2, rgb, DARK_NEUTRAL);
      ctx.fillStyle = `rgb(${Math.round(255 * rgb[0])}, ${Math.round(255 * rgb[1])}, ${Math.round(255 * rgb[2])})`;
      ctx.strokeStyle = 'rgba(20, 20, 22, 0.55)'; ctx.lineWidth = 0.6;
      ctx.beginPath();
      cell.vertices.forEach((v, j) => { const [x, y] = P(v); if (j) ctx.lineTo(x, y); else ctx.moveTo(x, y); });
      ctx.closePath(); ctx.fill(); ctx.stroke();
    }
    ctx.setLineDash([5, 4]); ctx.strokeStyle = 'rgba(255, 255, 255, 0.55)'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(cx - patch.sliceHalf * k, cy); ctx.lineTo(cx + patch.sliceHalf * k, cy); ctx.stroke();
    ctx.setLineDash([]);
    if (!arrows || arrows.k !== k) arrows = { k, cells: thin(patch.cells, k, 17) };
    const unit = 1.5 * spacing(layout) * k / ARROW_SPEED;
    for (const cell of arrows.cells) {
      const ue = f.currentEast[cell.i], un = f.currentNorth[cell.i], [ex, ey, nx, ny] = cell.basis;
      const vx = ue * ex + un * nx, vy = ue * ey + un * ny, len = Math.min(Math.hypot(ue, un) * unit, 3 * spacing(layout) * k);
      const norm = Math.hypot(vx, vy);
      if (len < 2 || !norm) continue;
      const [px, py] = P([cell.x, cell.y]), dx = vx / norm * len, dy = -vy / norm * len;
      arrow(ctx, px - dx / 2, py - dy / 2, px + dx / 2, py + dy / 2, { color: INK, width: 1.3, head: Math.min(5.5, 0.45 * len) });
    }
    ctx.restore();
    ctx.strokeStyle = GRID; ctx.lineWidth = 1; ctx.strokeRect(mx + 0.5, my + 0.5, S - 1, S - 1);

    const barMeters = patch.half > 5e6 ? 2e6 : 1e6, bar = barMeters * k, by = my + S + 14;
    ctx.strokeStyle = MUTED; ctx.lineWidth = 1.2;
    ctx.beginPath(); ctx.moveTo(mx, by - 3); ctx.lineTo(mx, by); ctx.lineTo(mx + bar, by); ctx.lineTo(mx + bar, by - 3); ctx.stroke();
    text(ctx, `${km(barMeters)} km`, mx + bar + 6, by - 1, { color: MUTED, size });
    const ref1 = ARROW_SPEED * unit, ax = mx + S - ref1;
    arrow(ctx, ax, by - 1, ax + ref1, by - 1, { color: INK, width: 1.1, head: 5 });
    text(ctx, `${ARROW_SPEED * 100} cm/s`, ax - 6, by - 1, { align: 'right', color: MUTED, size });
  }

  function thin(cells, k, gapPx) {
    const sorted = [...cells].sort((a, b) => Math.hypot(a.x, a.y) - Math.hypot(b.x, b.y)), kept = [];
    for (const cell of sorted) {
      if (Math.abs(cell.x) > patch.half || Math.abs(cell.y) > patch.half) continue;
      if (kept.every((o) => Math.hypot((o.x - cell.x) * k, (o.y - cell.y) * k) >= gapPx)) kept.push(cell);
    }
    return kept;
  }
}
