import { Figure, slider, choice, buttons, readout, legend, caption, text, arrow, rampRGB, termColor, clamp, INK, MUTED, GRID, LINE, DARK_NEUTRAL } from '../runtime.module.js';
import { Grid } from '../../js/grid.module.js';
import { basinWind, DEG, DAY } from '../oceanCases.module.js';
import { gyreReadings, DEPTH, ROW } from '../gyreReadings.module.js';

const N = 16, DAYS = 100, SPEED = 3 * DAY, WIDE = 520, LON = [-44, 44], LAT = [8, 52], COS = Math.cos(ROW * DEG);
const TAU_MAX = 0.2, SPACING = 4.33, ARROW_CAP = 1.4, REFERENCE = { surface: 0.05, depth: 0.015 };
const LADDER = [0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 7.5, 10, 15, 20, 30, 50];
const SURFACE = [2, 3, 5, 7.5, 10, 15, 20, 30, 50, 75, 100, 150, 200];
const GROUND = '#5a4a36', HALO = 'rgba(20, 20, 22, 0.85)', BAR = 'rgba(255, 232, 160, 0.85)', ARROW = 'rgba(255, 255, 255, 0.92)';

const cm = (v) => { const c = Math.abs(v) * 100; return c >= 10 ? c.toFixed(1) : c >= 1 ? c.toFixed(2) : c.toFixed(3); };
const km = (m) => `${(Math.round(m / 1e4) * 10).toLocaleString('en-US')} km`;
const sv = (v) => `${Math.abs(v).toFixed(1)} Sv`;
const way = (v, positive, negative) => (Math.abs(v) < 5e-6 ? 'none' : `${cm(v)} cm/s ${v > 0 ? positive : negative}`);
const windFactor = (tau0) => Math.max(tau0, 0.01) / 0.1;
const pick = (ladder, v) => ladder.find((s) => s >= v * 1.02) ?? ladder.at(-1);
const inside = (cell) => cell.lat > LAT[0] && cell.lat < LAT[1] && cell.lon > LON[0] && cell.lon < LON[1];
const degrees = (lon) => (lon === 0 ? '0°' : `${Math.abs(lon)}°${lon < 0 ? 'W' : 'E'}`);

export function mountGyre(root) {
  const controls = root.querySelector('.controls'), push = termColor('e');
  const state = { tau0: 0.1, beta: true, walls: true, arrows: 'surface' };
  let worker = null, run = 0, pending = false, owed = 0, layout = null, latest = null, reading = null, failed = null;
  let cells = null, coast = null, wallCoast = null, scale = null, surface = null, timer = 0;

  const fig = new Figure(root, { height: (w) => frame(w).total, step, draw });
  const words = caption(root, '');
  const key = legend(root, []);
  slider(controls, { label: 'Wind stress', min: 0, max: TAU_MAX, step: 0.01, value: state.tau0, format: (v) => `${v.toFixed(2)} N/m² at its strongest`, onInput: (v) => { state.tau0 = v; describe(); clearTimeout(timer); timer = setTimeout(restart, 150); } });
  choice(controls, { label: 'Coriolis parameter', options: [['grows toward the pole', 'beta'], ['the same everywhere', 'constant']], value: 'beta', onChange: (v) => { state.beta = v === 'beta'; describe(); restart(); } });
  choice(controls, { label: 'Walls', options: [['a basin', 'on'], ['none', 'off']], value: 'on', onChange: (v) => { state.walls = v === 'on'; describe(); restart(); } });
  choice(controls, { label: 'Arrows', options: [['the top layer’s current', 'surface'], ['the current averaged over the depth', 'depth']], value: state.arrows, onChange: (v) => { state.arrows = v; describe(); fig.render(); } });
  const [toggle] = buttons(controls, [['Pause', () => { if (latest?.done) restart(); else fig.play(!fig.running); label(); }], ['Restart', () => restart()]]);
  const outs = [readout(controls), readout(controls), readout(controls), readout(controls)];
  describe();

  function label() { toggle.textContent = fig.running ? 'Pause' : latest?.done ? 'Play again' : 'Play'; }

  function describe() {
    words.textContent = `A sea ${DEPTH.toLocaleString('en-US')} m deep${state.walls ? ', walled in from 40°W to 40°E and from 12° to 48°N,' : ' that covers the whole planet'} starts at rest under a wind that blows west over 15–30°N and east over 30–45°N, drawn at the right. The model’s own ocean runs in your browser on the N = ${N} grid, in steps of an hour and a half for ${DAYS} days, with all of its 58 layers, the top one 3 m thick. Arrows grow with the square root of the speed, so that slow currents show too; the arrow under the map is the length for the speed beside it. The bars are the northward speed across 30°N, averaged from the surface to the floor; the readouts average the western current over 25–35°N west of 32°W, and the interior over 25–35°N between 20°W and 20°E.`;
    const range = surface?.range ?? pick(SURFACE, 4.5 * windFactor(state.tau0));
    key.set([
      ['ramp', `sea surface ${range} cm below to ${range} cm above its average`, 'cool', 'warm', 'rgb(52, 55, 62)'],
      ['arrow', state.arrows === 'surface' ? 'current in the top layer, 3 m thick' : 'current averaged over the depth', INK],
      ['line', 'wind stress', push],
      ['bar', 'northward speed across 30°N', BAR],
      ...(state.beta && state.walls ? [['dash', 'Sverdrup’s balance', INK]] : []),
      ...(state.walls ? [['ramp', 'land', GROUND, GROUND, GROUND]] : [['dash', 'where the walls were', LINE]]),
    ]);
  }

  function restart() {
    clearTimeout(timer);
    if (!worker) return;
    pending = true; owed = 0; failed = null;
    worker.postMessage({ type: 'start', run: ++run, kind: 'basin', options: { N, days: DAYS, tau0: state.tau0, beta: state.beta, walls: state.walls } });
    fig.play(true);
    label();
  }

  function begin() {
    const grid = [...new Grid(N)];
    cells = [];
    for (const cell of grid) {
      const c = cell.centerVertex, lat = Math.atan2(c.z, Math.hypot(c.x, c.y)) / DEG, lon = Math.atan2(c.y, c.x) / DEG;
      if (lat < LAT[0] - 6 || lat > LAT[1] + 6 || lon < LON[0] - 8 || lon > LON[1] + 8) continue;
      const corner = (v) => [Math.atan2(v.y, v.x) / DEG, Math.atan2(v.z, Math.hypot(v.x, v.y)) / DEG];
      cells.push({ i: cell.index, lat, lon, corners: cell.vertices.map(corner), neighbors: cell.neighbors.map((n) => n.index) });
    }
    worker = new Worker(new URL('../ocean.worker.js', import.meta.url), { type: 'module' });
    worker.onmessage = receive;
    worker.onerror = (event) => { failed = event.message || 'the worker stopped'; fig.render(); };
    restart();
  }

  function coastOf(sea) {
    const segments = [];
    for (const cell of cells) {
      if (!sea[cell.i]) continue;
      cell.neighbors.forEach((n, k) => { if (!sea[n]) segments.push([cell.corners[(k + cell.corners.length - 1) % cell.corners.length], cell.corners[k]]); });
    }
    return segments;
  }

  function receive({ data }) {
    if (data.run !== run) return;
    if (data.type === 'layout') {
      layout = data; latest = null; reading = null; scale = null; surface = null;
      coast = data.walls ? coastOf(data.sea) : null;
      if (data.walls && !wallCoast) wallCoast = coast;
      return;
    }
    pending = false;
    if (data.type === 'error') { failed = data.message; fig.play(false); label(); fig.render(); return; }
    latest = data;
    reading = gyreReadings(data);
    rescale();
    report();
    if (data.done) { fig.play(false); label(); }
    fig.render();
  }

  function rescale() {
    const floor = 0.45 * windFactor(latest.tau0);
    let high = floor, low = floor;
    for (const v of latest.parallel.north) { const c = v / DEPTH * 100; high = Math.max(high, c); low = Math.max(low, -c); }
    if (latest.sverdrup !== null) low = Math.max(low, -reading.sverdrup.speed * 100);
    scale = { high: Math.max(scale?.high ?? 0, pick(LADDER, high)), low: Math.max(scale?.low ?? 0, pick(LADDER, low)) };
    const sea = layout.sea, eta = latest.eta;
    let mean = 0, count = 0, extent = 0;
    for (const cell of cells) if (sea[cell.i] && inside(cell)) { mean += eta[cell.i]; count++; }
    mean = count ? mean / count : 0;
    for (const cell of cells) if (sea[cell.i] && inside(cell)) extent = Math.max(extent, Math.abs(eta[cell.i] - mean) * 100);
    const range = Math.max(surface?.range ?? 0, pick(SURFACE, Math.max(4.5 * windFactor(latest.tau0), extent)));
    const changed = range !== surface?.range;
    surface = { mean, range };
    if (changed) describe();
  }

  function report() {
    const r = reading, every = latest.done ? 1 : 2;
    outs[0].set([['day', r.day.toFixed(1)]], every);
    if (!latest.walls) {
      const { westward, eastward } = r.zonal;
      outs[1].set([['across 30°N', `${sv(r.across.north)} north and ${sv(r.across.south)} south`]], every);
      outs[2].set([['strongest westward flow', westward && westward.speed < 0 ? `${cm(westward.speed)} cm/s at ${westward.latitude}°N` : 'none'], ['strongest eastward', eastward && eastward.speed > 0 ? `${cm(eastward.speed)} cm/s at ${eastward.latitude}°N` : 'none']], every);
      outs[3].set([['Sverdrup’s balance', 'cannot hold: no wall to send the water back north']], every);
      return;
    }
    const current = r.current;
    if (latest.beta) outs[1].set([['western current', way(r.west, 'north', 'south')], ['width at 30°N', current ? km(current.width) : 'none'], ['carrying', current ? sv(current.transport) : 'nothing north']], every);
    else outs[1].set([['along the western wall', way(r.west, 'north', 'south')], ['carrying', current ? sv(current.transport) : 'nothing north']], every);
    if (latest.beta) {
      outs[2].set([['the interior’s return', way(r.interior, 'north', 'south')], ['the western current is', r.ratio ? `${r.ratio.toFixed(1)} times as fast` : 'no faster']], every);
      outs[3].set([['Sverdrup’s balance predicts', Math.abs(r.sverdrup.speed) < 5e-6 ? 'no flow without a wind' : `${way(r.sverdrup.speed, 'north', 'south')}, ${sv(r.sverdrup.transport)} across the basin`]], every);
    } else {
      outs[2].set([['the interior', way(r.interior, 'north', 'south')], ['eastern current', way(r.east, 'north', 'south')]], every);
      outs[3].set([['Sverdrup’s balance', 'none: f does not change with latitude'], ['across 30°N', `${sv(r.across.north)} each way`]], every);
    }
  }

  function step(dt) {
    if (!worker) { begin(); return false; }
    owed += dt;
    if (!pending && latest && !latest.done) { pending = true; worker.postMessage({ type: 'advance', seconds: owed * SPEED }); owed = 0; }
    return false;
  }

  function frame(w) {
    const narrow = w < WIDE, size = narrow ? 10 : 11;
    const left = narrow ? 30 : 40, gap = narrow ? 8 : 14, windW = narrow ? 50 : 96, right = narrow ? 6 : 12, top = narrow ? 24 : 28;
    const mapW = w - left - gap - windW - right, k = mapW / ((LON[1] - LON[0]) * COS), mapH = k * (LAT[1] - LAT[0]);
    const chartTop = top + mapH + (narrow ? 52 : 58), chartH = narrow ? 92 : 120, total = Math.round(chartTop + chartH + (narrow ? 12 : 16));
    const x = (lon) => left + (lon - LON[0]) * COS * k, y = (lat) => top + (LAT[1] - lat) * k;
    return { narrow, size, left, top, mapW, mapH, k, x, y, windLeft: left + mapW + gap, windW, chartTop, chartBottom: chartTop + chartH, total };
  }

  function draw(ctx, w, h) {
    const g = frame(w);
    if (failed) { text(ctx, `The model could not run: ${failed}`, w / 2, h / 2, { align: 'center', color: MUTED, size: 12 }); return; }
    if (!latest || !layout || !cells) { text(ctx, 'Starting the model…', w / 2, h / 2, { align: 'center', color: MUTED, size: 12 }); return; }
    map(ctx, g);
    wind(ctx, g);
    chart(ctx, g);
  }

  function map(ctx, g) {
    const { left, top, mapW, mapH, x, y, size, narrow } = g, sea = layout.sea, eta = latest.eta, { mean } = surface, range = surface.range / 100, rgb = [0, 0, 0];
    ctx.save();
    ctx.beginPath(); ctx.rect(left, top, mapW, mapH); ctx.clip();
    for (const cell of cells) {
      if (sea[cell.i]) { rampRGB(0.5 + clamp((eta[cell.i] - mean) / range, -1, 1) / 2, rgb, DARK_NEUTRAL); ctx.fillStyle = `rgb(${Math.round(rgb[0] * 255)}, ${Math.round(rgb[1] * 255)}, ${Math.round(rgb[2] * 255)})`; }
      else ctx.fillStyle = GROUND;
      ctx.strokeStyle = 'rgba(255, 255, 255, 0.07)'; ctx.lineWidth = 1;
      ctx.beginPath();
      cell.corners.forEach(([lon, lat], n) => (n ? ctx.lineTo(x(lon), y(lat)) : ctx.moveTo(x(lon), y(lat))));
      ctx.closePath(); ctx.fill(); ctx.stroke();
    }
    const outline = latest.walls ? coast : wallCoast;
    if (outline) {
      ctx.strokeStyle = latest.walls ? 'rgba(255, 255, 255, 0.55)' : LINE; ctx.lineWidth = latest.walls ? 1.5 : 1.2; ctx.setLineDash(latest.walls ? [] : [4, 4]);
      ctx.beginPath();
      for (const [[lon0, lat0], [lon1, lat1]] of outline) { ctx.moveTo(x(lon0), y(lat0)); ctx.lineTo(x(lon1), y(lat1)); }
      ctx.stroke(); ctx.setLineDash([]);
    }
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.3)'; ctx.setLineDash([3, 4]); ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(left, y(ROW)); ctx.lineTo(left + mapW, y(ROW)); ctx.stroke(); ctx.setLineDash([]);
    const depth = state.arrows === 'depth', ref = REFERENCE[state.arrows], cellPx = SPACING * g.k, head = narrow ? 3.5 : 4.5;
    const east = depth ? latest.transportEast : latest.currentEast, north = depth ? latest.transportNorth : latest.currentNorth;
    for (const cell of cells) {
      if (!sea[cell.i] || cell.lat < LAT[0] + 1 || cell.lat > LAT[1] - 1 || cell.lon < LON[0] + 1 || cell.lon > LON[1] - 1) continue;
      const u = east[cell.i] / (depth ? DEPTH : 1), v = north[cell.i] / (depth ? DEPTH : 1), speed = Math.hypot(u, v);
      const length = cellPx * Math.min(ARROW_CAP, Math.sqrt(speed / ref));
      if (length < 2) continue;
      const cx = x(cell.lon), cy = y(cell.lat), dx = u / speed * length / 2, dy = -v / speed * length / 2;
      arrow(ctx, cx - dx, cy - dy, cx + dx, cy + dy, { color: ARROW, width: narrow ? 1.1 : 1.3, head });
    }
    ctx.restore();
    ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.strokeRect(left + 0.5, top + 0.5, mapW - 1, mapH - 1);
    for (const lat of [15, 30, 45]) text(ctx, `${lat}°N`, left - 5, y(lat), { align: 'right', color: MUTED, size });
    for (const lon of narrow ? [-40, 0, 40] : [-40, -20, 0, 20, 40]) text(ctx, degrees(lon), x(lon), top + mapH + 11, { align: 'center', color: MUTED, size });
    const ry = top + mapH + 27, reference = `${(ref * 100).toLocaleString('en-US')} cm/s`;
    arrow(ctx, left, ry, left + cellPx, ry, { color: ARROW, width: narrow ? 1.1 : 1.3, head });
    text(ctx, reference, left + cellPx + 6, ry, { color: MUTED, size });
    text(ctx, `day ${Math.floor(latest.day)}`, left, top - 12, { color: INK, size: size + 1, weight: 500 });
    const mode = `${latest.beta ? 'f grows toward the pole' : 'f the same everywhere'}${latest.walls ? '' : ', no walls'}`;
    text(ctx, mode, left + mapW, top - 12, { align: 'right', color: MUTED, size });
  }

  function wind(ctx, g) {
    const { windLeft, windW, top, mapH, y, size, narrow } = g, zero = windLeft + windW / 2, per = (windW / 2) / TAU_MAX, tau0 = latest.tau0;
    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(zero, top); ctx.lineTo(zero, top + mapH); ctx.stroke();
    ctx.strokeStyle = GRID;
    for (const t of [-TAU_MAX, TAU_MAX]) { ctx.beginPath(); ctx.moveTo(zero + t * per, top); ctx.lineTo(zero + t * per, top + mapH); ctx.stroke(); }
    ctx.strokeStyle = push; ctx.lineWidth = 1.5;
    ctx.beginPath();
    for (let lat = LAT[0]; lat <= LAT[1] + 1e-9; lat += 0.25) { const px = zero + tau0 * basinWind(lat * DEG) * per; if (lat === LAT[0]) ctx.moveTo(px, y(lat)); else ctx.lineTo(px, y(lat)); }
    ctx.stroke();
    for (let lat = 16; lat <= 44; lat += 4) {
      const t = tau0 * basinWind(lat * DEG) * per;
      if (Math.abs(t) >= 3) arrow(ctx, zero, y(lat), zero + t, y(lat), { color: push, width: 1.5, head: narrow ? 4 : 5 });
    }
    text(ctx, 'wind', zero, top - 12, { align: 'center', color: push, size: size + 1, weight: 500 });
    const small = narrow ? 9 : size;
    text(ctx, 'west', windLeft, top + mapH + 24, { color: MUTED, size: small });
    text(ctx, 'east', windLeft + windW, top + mapH + 24, { align: 'right', color: MUTED, size: small });
    if (!narrow) {
      text(ctx, 'trades', zero + 5, y(20), { color: MUTED, size });
      text(ctx, 'westerlies', zero - 5, y(40), { align: 'right', color: MUTED, size });
    }
  }

  function chart(ctx, g) {
    const { left, mapW, x, chartTop, chartBottom, size, narrow } = g, { high, low } = scale, row = latest.parallel;
    const yv = (c) => chartTop + (high - c) / (high + low) * (chartBottom - chartTop);
    ctx.strokeStyle = GRID; ctx.lineWidth = 1;
    for (const c of [high, -low]) { ctx.beginPath(); ctx.moveTo(left, yv(c)); ctx.lineTo(left + mapW, yv(c)); ctx.stroke(); }
    ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(left, yv(0)); ctx.lineTo(left + mapW, yv(0)); ctx.stroke();
    const value = (c) => `${c}`.replace('-', '−');
    text(ctx, value(high), left - 5, yv(high), { align: 'right', color: MUTED, size });
    text(ctx, '0', left - 5, yv(0), { align: 'right', color: MUTED, size });
    text(ctx, value(-low), left - 5, yv(-low), { align: 'right', color: MUTED, size });
    text(ctx, 'cm/s', left - 5, chartTop - 14, { align: 'right', color: MUTED, size });
    text(ctx, narrow ? 'northward across 30°N' : 'northward speed across 30°N, surface to floor', left + mapW, chartTop - 14, { align: 'right', color: INK, size, weight: 500 });
    const n = row.longitude.length;
    for (let s = 0; s < n; s++) {
      const lon = row.longitude[s], spacing = s + 1 < n ? row.longitude[s + 1] - lon : lon - row.longitude[s - 1];
      const half = 0.36 * spacing * COS * g.k, c = row.north[s] / DEPTH * 100, top = yv(Math.max(0, Math.min(high, c))), bottom = yv(Math.min(0, Math.max(-low, c)));
      ctx.fillStyle = BAR; ctx.fillRect(x(lon) - half, top, 2 * half, Math.max(1, bottom - top));
    }
    if (reading.sverdrup && latest.walls) {
      const from = reading.current ? reading.current.end : row.longitude[0], to = row.longitude[n - 1] + 0.5 * (row.longitude[n - 1] - row.longitude[n - 2]), yy = yv(reading.sverdrup.speed * 100);
      ctx.strokeStyle = INK; ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]);
      ctx.beginPath(); ctx.moveTo(x(from), yy); ctx.lineTo(x(to), yy); ctx.stroke(); ctx.setLineDash([]);
      text(ctx, 'Sverdrup', x(to) - 2, yy + 10, { align: 'right', color: INK, size, halo: HALO });
    }
    const current = reading.current;
    if (current && latest.beta) {
      const x0 = x(current.coast), x1 = x(current.end), yy = chartTop - 4;
      ctx.strokeStyle = INK; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(x0, yy + 4); ctx.lineTo(x0, yy); ctx.lineTo(x1, yy); ctx.lineTo(x1, yy + 4); ctx.stroke();
      text(ctx, km(current.width), x1 + 5, yy, { color: INK, size });
    }
  }

  fig.play(true);
  label();
}
