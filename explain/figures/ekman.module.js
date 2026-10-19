import { Figure, slider, choice, buttons, readout, legend, caption, text, arrow, termColor, rampRGB, clamp, DARK_NEUTRAL, ACCENT, INK, MUTED, GRID, LINE } from '../runtime.module.js';
import { Grid } from '../../js/grid.module.js';
import { ekmanSpiral, fitViscosity, windForStress, tracks, coriolis, inertialPeriod, columnDays, RHO, DAY } from '../ekmanTheory.module.js';

const N = 8, PUMP_DEPTH = 100, WIDE = 560, LOOP_SECONDS = 4, PLANET_DAYS = 4, PLANET_SPEED = 8 * 3600, CENTER_LONGITUDE = -20 * Math.PI / 180;
const LATITUDES = [...Array.from({ length: 16 }, (_, k) => -80 + 5 * k), ...Array.from({ length: 16 }, (_, k) => 5 + 5 * k)];
const FAINT = 'rgba(255,255,255,0.55)', HALO = 'rgba(20, 20, 22, 0.85)', DARK = `rgb(${DARK_NEUTRAL.join(', ')})`;
const NICE = [1, 1.5, 2, 3, 5, 7];

const place = (latitude) => (latitude ? `${Math.abs(latitude)}°${latitude > 0 ? 'N' : 'S'}` : '0°');
const tidy = (v) => Number(v.toPrecision(3));
const niceStep = (span) => { const p = 10 ** Math.floor(Math.log10(span)), m = span / p; return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p; };
const niceCeil = (v) => { const p = 10 ** Math.floor(Math.log10(v)); return (NICE.find((m) => m * p >= v * (1 - 1e-9)) ?? 10) * p; };
const css = (rgb) => `rgb(${Math.round(rgb[0] * 255)}, ${Math.round(rgb[1] * 255)}, ${Math.round(rgb[2] * 255)})`;

export function mountEkman(root) {
  const controls = root.querySelector('.controls'), color = { a: termColor('a'), e: termColor('e') };
  let view = 'column', tau = 0.1, latitude = 45, wind = 'trades', match = null, matchFrame = null;
  let layout = null, latest = null, pending = false, run = 0, owed = 0, speed = 1, plan = null, scale = 0, path = null, pathFrame = null, cells = null, restartTimer = 0;

  let failed = null;
  const worker = new Worker(new URL('../ocean.worker.js', import.meta.url), { type: 'module' });
  worker.onerror = (event) => { event.preventDefault(); failed = event.message || 'the worker could not start'; fig.play(false); fig.render(); };
  worker.onmessage = ({ data }) => {
    if (data.run !== run) return;
    if (data.type === 'layout') { layout = data; return; }
    pending = false;
    if (data.type === 'error') { failed = data.message; console.error(`ekman: ${data.message}`); fig.render(); return; }
    if (data.steps === 0 && latest && latest.time === data.time) return;
    latest = data;
    if (data.view === 'column' && !plan) plan = planBox(data);
    if (data.view === 'planet') {
      let most = 0;
      for (const w of data.upwelling) most = Math.max(most, Math.abs(w) * DAY);
      if (most > 0 && niceCeil(most) > scale) { scale = niceCeil(most); key.set(legendItems()); }
    }
    if (data.done) fig.play(false);
    report();
    fig.render();
  };

  const fig = new Figure(root, { height: (w) => (w < WIDE ? 490 : 440), step, draw });
  const note = caption(root, '');
  const key = legend(root, []);
  choice(controls, { label: 'Show', options: [['one latitude', 'column'], ['the whole planet', 'planet']], value: view, onChange: (v) => { view = v; reveal(); restart(); }, span: true });
  slider(controls, { label: 'Wind stress', min: 0.02, max: 0.3, step: 0.01, value: tau, format: (v) => `${v.toFixed(2)} N/m² (wind ${windForStress(v).toFixed(1)} m/s)`, onInput: (v) => { tau = v; later(); }, span: true });
  slider(controls, { label: 'Latitude', min: 0, max: LATITUDES.length - 1, step: 1, value: LATITUDES.indexOf(latitude), format: (k) => place(LATITUDES[k]), onInput: (k) => { latitude = LATITUDES[k]; if (view === 'column') later(); else { report(); fig.render(); } }, span: true });
  choice(controls, { label: 'Wind over the planet', options: [['trades and westerlies', 'trades'], ['westerlies only', 'westerlies']], value: wind, onChange: (v) => { wind = v; restart(); }, span: true });
  const planetRow = controls.lastElementChild;
  buttons(controls, [['Restart', restart]]);
  const out = readout(controls), theoryOut = readout(controls), theoryRow = controls.lastElementChild;

  function reveal() {
    planetRow.style.display = view === 'planet' ? '' : 'none';
    theoryRow.style.display = view === 'column' ? '' : 'none';
    note.textContent = view === 'column'
      ? 'The model’s own ocean on the N = 8 grid, with its 58 layers, at rest until a band of eastward wind stress, about 20° wide, starts to blow along the chosen latitude. The ocean steps every 22.5 minutes, and everything shown is averaged over the ring of cells nearest that latitude. The wind’s push goes into the top layer, 3 m thick, and the model’s turbulent mixing carries it down through the layers beneath. The panel of paths shows where the top layer, the mixed layer and the water well below them would carry a float, seen from above. The velocity panel shows the top layer’s and the mixed layer’s velocity measured against the water well below, their tips tracing hodographs as the velocities change, with the two forces that act on the mixed layer. The third panel shows the model’s current in each layer, averaged over the last two inertial periods, beside Ekman’s spiral for the one eddy viscosity that fits it best, which is theory, not the model. Time runs as fast as your computer allows, up to one inertial period every four seconds.'
      : `The same ocean with the wind over the whole planet, ${wind === 'trades' ? 'the trades blowing west within 30° of the equator and the westerlies east from 30° to 65°' : 'westerlies only, centered on 45° in each hemisphere'}, for ${PLANET_DAYS} days. The strip on the left is the wind’s stress at each latitude. On the globe the colors are the water moving up or down ${PUMP_DEPTH} m below the surface and the arrows the Ekman transport, the water above twice the boundary layer’s depth times its current against the water below that, both averaged over the last two days; arrows longer than twice the one at the bottom are drawn at that length. The dashed line is the latitude the readout describes.`;
    key.set(legendItems());
  }

  function legendItems() {
    if (view === 'column') return [
      ['line', 'the top layer, 3 m thick: its path, and its velocity against the water below', INK],
      ['line', 'the mixed layer: the same', ACCENT],
      ['dash', 'the path of the water well below', FAINT],
      ['arrow', 'wind', color.e],
      ['force', 'the wind’s push on the mixed layer', color.e],
      ['force', 'Coriolis force on its drift', color.a],
      ['arrow', 'the model’s current in a layer, in the third panel', ACCENT],
      ['dash', 'Ekman’s spiral, in theory', INK],
    ];
    const s = scale ? `${tidy(scale)}` : '…';
    return [
      ['ramp', `water ${PUMP_DEPTH} m down rising ${s} m a day to sinking ${s} m a day`, 'cool', 'warm', DARK],
      ['arrow', 'Ekman transport: the water above twice the boundary layer’s depth, against the water below', INK],
      ['line', 'wind stress', color.e],
    ];
  }

  function later() { clearTimeout(restartTimer); restartTimer = setTimeout(restart, 120); }

  function restart() {
    clearTimeout(restartTimer);
    latest = null; layout = null; plan = null; scale = 0; path = null; pathFrame = null; owed = 0; pending = true; failed = null;
    const options = view === 'column' ? { view, tau, latitude, days: columnDays(latitude) } : { view, tau, wind, days: PLANET_DAYS };
    speed = view === 'column' ? inertialPeriod(latitude) / LOOP_SECONDS : PLANET_SPEED;
    worker.postMessage({ type: 'start', run: ++run, kind: 'ekman', options });
    reveal();
    report();
    fig.play(true);
    fig.render();
  }

  function step(dt) {
    owed += dt;
    if (!pending && latest && !latest.done) { pending = true; worker.postMessage({ type: 'advance', seconds: owed * speed }); owed = 0; }
    return false;
  }

  function planBox(frame) {
    const reach = frame.theory / 20 * columnDays(frame.latitude) * DAY, south = frame.f > 0;
    return { west: -0.2 * reach, east: 0.45 * reach, south: south ? -1.05 * reach : -0.1 * reach, north: south ? 0.1 * reach : 1.05 * reach };
  }

  function matched() {
    if (matchFrame !== latest) {
      const p = latest.profile;
      const fit = fitViscosity({ tau: latest.stressMean, f: latest.fMean, center: p.center, thickness: p.thickness, east: p.east, north: p.north, reach: latest.reach });
      match = { ...fit, spiral: ekmanSpiral({ tau: latest.stressMean, f: latest.fMean, viscosity: fit.viscosity }) };
      matchFrame = latest;
    }
    return match;
  }

  function report() {
    if (!latest || latest.view !== view) { out.set([['starting', '…']]); theoryOut.set([]); return; }
    if (view === 'column') {
      const f = latest, side = f.f > 0 ? 'right' : 'left', tr = f.transport, cms = (v) => (v * 100).toFixed(1);
      const drift = f.settled
        ? [['Ekman transport', `${tr.magnitude.toFixed(3)} m²/s, ${tr.ratio.toFixed(2)} times τ/(ρ₀f), ${Math.abs(tr.angle).toFixed(0)}° to the ${side} of the wind`], ['top layer', `${cms(f.top.magnitude)} cm/s, ${Math.abs(f.top.angle).toFixed(0)}° to the ${side}`], ['mixed layer', `${f.mixedDepth.toFixed(0)} m deep, drifting ${cms(f.mixed.magnitude)} cm/s, ${Math.abs(f.mixed.angle).toFixed(0)}° to the ${side}`], ['boundary layer', `${f.boundaryDepth.toFixed(0)} m deep`]]
        : [['drift', `measured once two inertial periods have passed, at ${(2 * f.inertialPeriod / 3600).toFixed(0)} h`], ['mixed layer', `${f.mixedDepth.toFixed(0)} m deep`]];
      out.set([['averaged over', ringSpan()], ['inertial period', `${(f.inertialPeriod / 3600).toFixed(1)} h`], ...drift]);
      const m = matched(), s = m.spiral;
      theoryOut.set([['Ekman’s spiral, in theory, at the viscosity that fits the model', `${m.viscosity.toFixed(3)} m²/s: ${(Math.abs(s.surface) * 100).toFixed(1)} cm/s at the surface, 45° to the ${side} of the wind, reversing at ${s.reversal.toFixed(0)} m`], ['the model’s own viscosity, averaged over its boundary layer', `${f.meanViscosity.toFixed(3)} m²/s`]]);
      return;
    }
    if (!layout) return;
    const band = bandMean();
    const rows = [[`at ${place(latitude)}, wind`, Math.abs(band.stress) < 5e-4 ? 'none' : `${Math.abs(band.stress).toFixed(3)} N/m² ${band.stress > 0 ? 'east' : 'west'}`]];
    const flow = (m) => (Math.abs(m) < 5e-3 ? 'none' : `${Math.abs(m).toFixed(2)} m²/s ${m < 0 ? 'south' : 'north'}`);
    rows.push(['Ekman transport', Math.abs(latitude) > 5 ? `${flow(band.ekman)} against τ/(ρ₀f) = ${flow(band.theory)}` : `${flow(band.ekman)}; τ/(ρ₀f) is no guide this close to the equator`]);
    rows.push([`${PUMP_DEPTH} m down the water`, vertical(band.upwelling)]);
    rows.push(['on the equator it', vertical(latest.equatorUpwelling)]);
    out.set(rows);
  }

  function ringSpan() {
    if (!layout?.ring) return '…';
    let low = Infinity, high = -Infinity;
    for (const i of layout.ring) { low = Math.min(low, layout.latitude[i]); high = Math.max(high, layout.latitude[i]); }
    const at = (lat) => `${Math.abs(lat).toFixed(1)}°${lat >= 0 ? 'N' : 'S'}`;
    return `the cells from ${at(latitude > 0 ? low : high)} to ${at(latitude > 0 ? high : low)}`;
  }

  function vertical(w) {
    const perDay = w * DAY;
    return Math.abs(perDay) < 0.005 ? 'barely moves' : `${perDay > 0 ? 'rises' : 'sinks'} ${Math.abs(perDay).toFixed(2)} m a day`;
  }

  function bandMean() {
    const lat = layout.latitude, sum = { stress: 0, ekman: 0, theory: 0, upwelling: 0 };
    let count = 0;
    for (let i = 0; i < lat.length; i++) {
      if (Math.abs(lat[i] - latitude) >= 2.5) continue;
      sum.stress += latest.stress[i]; sum.ekman += latest.ekmanNorth[i]; sum.upwelling += latest.upwelling[i];
      sum.theory += -latest.stress[i] / (RHO * coriolis(lat[i]));
      count++;
    }
    for (const name in sum) sum[name] /= Math.max(1, count);
    return sum;
  }

  function draw(ctx, w, h) {
    if (failed) { text(ctx, `The model’s ocean could not run: ${failed}`, w / 2, h / 2, { align: 'center', color: MUTED }); return; }
    if (!latest || latest.view !== view || !layout) { text(ctx, 'starting the model’s ocean…', w / 2, h / 2, { align: 'center', color: MUTED }); return; }
    if (view === 'planet') { drawPlanet(ctx, w, h); return; }
    if (pathFrame !== latest) { path = tracks(latest.series); pathFrame = latest; }
    const wide = w >= WIDE;
    const left = wide ? Math.round(w * 0.52) : w, lower = wide ? h : h - 205;
    const panels = wide
      ? [{ x: 0, y: 0, w: left, h }, { x: left, y: 0, w: w - left, h: Math.round(h / 2) }, { x: left, y: Math.round(h / 2), w: w - left, h: h - Math.round(h / 2) }]
      : [{ x: 0, y: 0, w, h: lower }, { x: 0, y: lower, w: Math.round(w / 2), h: h - lower }, { x: Math.round(w / 2), y: lower, w: w - Math.round(w / 2), h: h - lower }];
    ctx.strokeStyle = GRID; ctx.lineWidth = 1;
    ctx.beginPath();
    if (wide) { ctx.moveTo(left + 0.5, 0); ctx.lineTo(left + 0.5, h); ctx.moveTo(left, panels[2].y + 0.5); ctx.lineTo(w, panels[2].y + 0.5); }
    else { ctx.moveTo(0, lower + 0.5); ctx.lineTo(w, lower + 0.5); ctx.moveTo(panels[2].x + 0.5, lower); ctx.lineTo(panels[2].x + 0.5, h); }
    ctx.stroke();
    drawPath(ctx, panels[0]);
    drawHodograph(ctx, panels[1], !wide);
    drawSpiral(ctx, panels[2], !wide);
  }

  function clip(ctx, P, body) { ctx.save(); ctx.beginPath(); ctx.rect(P.x, P.y, P.w, P.h); ctx.clip(); body(); ctx.restore(); }

  function polyline(ctx, points, sx, sy, count = points.length / 2) {
    ctx.beginPath();
    for (let s = 0; s < count; s++) { if (s) ctx.lineTo(sx(points[2 * s]), sy(points[2 * s + 1])); else ctx.moveTo(sx(points[0]), sy(points[1])); }
    ctx.stroke();
  }

  function drawPath(ctx, P) {
    const f = latest, n = f.series.time.length, box = { ...plan }, tracks3 = [path.top, path.mixed, path.deep];
    for (const track of tracks3) for (let s = 0; s < n; s++) {
      box.west = Math.min(box.west, track[2 * s]); box.east = Math.max(box.east, track[2 * s]);
      box.south = Math.min(box.south, track[2 * s + 1]); box.north = Math.max(box.north, track[2 * s + 1]);
    }
    const top = P.y + 50, bottom = P.y + P.h - 26, room = { w: P.w - 28, h: bottom - top };
    const k = Math.min(room.w / (box.east - box.west), room.h / (box.north - box.south));
    const cx = P.x + P.w / 2 - k * (box.west + box.east) / 2, cy = (top + bottom) / 2 + k * (box.south + box.north) / 2;
    const sx = (x) => cx + k * x, sy = (y) => cy - k * y;
    const spacing = niceStep(Math.max(box.east - box.west, box.north - box.south) / 6);
    clip(ctx, P, () => {
      ctx.strokeStyle = GRID; ctx.lineWidth = 1;
      ctx.beginPath();
      for (let x = Math.floor((P.x - cx) / k / spacing) * spacing; sx(x) < P.x + P.w; x += spacing) { ctx.moveTo(Math.round(sx(x)) + 0.5, P.y + 40); ctx.lineTo(Math.round(sx(x)) + 0.5, P.y + P.h); }
      for (let y = Math.floor((cy - P.y - P.h) / k / spacing) * spacing; sy(y) > P.y + 40; y += spacing) { ctx.moveTo(P.x, Math.round(sy(y)) + 0.5); ctx.lineTo(P.x + P.w, Math.round(sy(y)) + 0.5); }
      ctx.stroke();
      ctx.lineJoin = 'round';
      ctx.strokeStyle = FAINT; ctx.lineWidth = 1.5; ctx.setLineDash([4, 4]); polyline(ctx, path.deep, sx, sy, n); ctx.setLineDash([]);
      ctx.strokeStyle = ACCENT; ctx.lineWidth = 2; polyline(ctx, path.mixed, sx, sy, n);
      ctx.strokeStyle = INK; ctx.lineWidth = 1.5; polyline(ctx, path.top, sx, sy, n);
      ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(sx(0), sy(0), 4, 0, 2 * Math.PI); ctx.stroke();
      ctx.fillStyle = MUTED; ctx.beginPath(); ctx.arc(sx(path.deep[2 * n - 2]), sy(path.deep[2 * n - 1]), 3.5, 0, 2 * Math.PI); ctx.fill();
      ctx.fillStyle = ACCENT; ctx.beginPath(); ctx.arc(sx(path.mixed[2 * n - 2]), sy(path.mixed[2 * n - 1]), 4, 0, 2 * Math.PI); ctx.fill();
      ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(sx(path.top[2 * n - 2]), sy(path.top[2 * n - 1]), 4.5, 0, 2 * Math.PI); ctx.fill();
      text(ctx, 'start', sx(0) - 8, sy(0), { align: 'right', color: MUTED, size: 10, halo: HALO });
    });
    text(ctx, 'paths, seen from above', P.x + 12, P.y + 14, { color: INK, size: 11.5 });
    text(ctx, f.time < 100 * 3600 ? `${(f.time / 3600).toFixed(0)} h` : `day ${(f.time / DAY).toFixed(1)}`, P.x + P.w - 12, P.y + 14, { align: 'right', color: MUTED, size: 11 });
    arrow(ctx, P.x + 12, P.y + 34, P.x + 62, P.y + 34, { color: color.e, width: 3, head: 9 });
    text(ctx, `wind, ${tau.toFixed(2)} N/m²`, P.x + 70, P.y + 34, { color: MUTED, size: 10.5 });
    text(ctx, 'north ↑', P.x + P.w - 12, P.y + 34, { align: 'right', color: MUTED, size: 10.5 });
    text(ctx, 'east →', P.x + P.w - 12, P.y + P.h - 12, { align: 'right', color: MUTED, size: 10.5 });
    text(ctx, `grid every ${spacing >= 1000 ? `${tidy(spacing / 1000)} km` : `${tidy(spacing)} m`}`, P.x + 12, P.y + P.h - 12, { color: MUTED, size: 10.5, halo: HALO });
  }

  function label(ctx, P, str, x, y, { align = 'left', color: c = MUTED, size = 10 } = {}) {
    ctx.font = `400 ${size}px system-ui, sans-serif`;
    const width = ctx.measureText(str).width, start = align === 'left' ? x : align === 'right' ? x - width : x - width / 2;
    const fitted = clamp(start, P.x + 4, P.x + P.w - 4 - width);
    text(ctx, str, fitted, clamp(y, P.y + 8, P.y + P.h - 8), { color: c, size, halo: HALO });
  }

  function title(ctx, P, lines) {
    lines.forEach((line, k) => text(ctx, line, P.x + 12, P.y + 14 + 13 * k, { color: k ? MUTED : INK, size: k ? 10.5 : 11.5 }));
    return 14 + 13 * lines.length;
  }

  function drawHodograph(ctx, P, narrow) {
    const f = latest, s = f.series, n = s.time.length, sign = Math.sign(f.f), h1 = s.mixedDepth[n - 1], D0 = f.theory / Math.max(h1, 1);
    const head = title(ctx, P, narrow ? ['velocity against', 'the water well below'] : ['velocity against the water well below']);
    const rel = (name, deep) => Float64Array.from(s[name], (v, j) => v - s[deep][j]);
    const topE = rel('topEast', 'deepEast'), topN = rel('topNorth', 'deepNorth'), mixE = rel('mixedEast', 'deepEast'), mixN = rel('mixedNorth', 'deepNorth');
    let D = D0;
    for (let j = 0; j < n; j++) D = Math.max(D, Math.hypot(topE[j], topN[j]));
    const room = { x: P.x + 14, y: P.y + head + 6, w: P.w - 28, h: P.h - head - 18 };
    const k = Math.min(room.w, room.h) / (2.7 * D), cx = room.x + room.w / 2, cy = room.y + room.h / 2 - sign * k * D;
    const sx = (u) => cx + k * u, sy = (v) => cy - k * v;
    const ring = niceStep(1.2 * D);
    clip(ctx, P, () => {
      ctx.save(); ctx.beginPath(); ctx.rect(P.x, room.y - 4, P.w, P.y + P.h - room.y + 4); ctx.clip();
      ctx.strokeStyle = GRID; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(P.x, sy(0)); ctx.lineTo(P.x + P.w, sy(0)); ctx.moveTo(sx(0), room.y - 4); ctx.lineTo(sx(0), P.y + P.h); ctx.stroke();
      for (let r = ring; r < 2.3 * D; r += ring) { ctx.beginPath(); ctx.arc(sx(0), sy(0), r * k, 0, 2 * Math.PI); ctx.stroke(); }
      ctx.restore();
      ctx.font = '400 9.5px system-ui, sans-serif';
      for (let r = ring; r < 2.3 * D; r += ring) {
        const name = `${(r * 100).toFixed(ring * 100 < 1 ? 1 : 0)} cm/s`;
        if (sx(r) + 3 + ctx.measureText(name).width < P.x + P.w - 4) label(ctx, P, name, sx(r) + 3, sy(0) - sign * 8, { size: 9.5 });
      }
      const tx = sx(0), ty = sy(-sign * D0);
      ctx.strokeStyle = INK; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.moveTo(tx - 5, ty - 5); ctx.lineTo(tx + 5, ty + 5); ctx.moveTo(tx - 5, ty + 5); ctx.lineTo(tx + 5, ty - 5); ctx.stroke();
      label(ctx, P, 'τ/(ρ₀fh)', tx, ty + sign * 16, { align: 'center', color: INK, size: 10.5 });
      ctx.lineJoin = 'round'; ctx.globalAlpha = 0.75; ctx.lineWidth = 1.3;
      for (const [east, north, c] of [[topE, topN, INK], [mixE, mixN, ACCENT]]) {
        ctx.strokeStyle = c;
        ctx.beginPath();
        for (let j = 0; j < n; j++) { const x = sx(east[j]), y = sy(north[j]); if (j) ctx.lineTo(x, y); else ctx.moveTo(x, y); }
        ctx.stroke();
      }
      ctx.globalAlpha = 1;
      arrow(ctx, sx(0), sy(0), sx(topE[n - 1]), sy(topN[n - 1]), { color: INK, width: 1.6, head: 6 });
      const u = mixE[n - 1], v = mixN[n - 1], px = sx(u), py = sy(v);
      arrow(ctx, sx(0), sy(0), px, py, { color: ACCENT, width: 2, head: 7 });
      const reach = 0.6 / Math.abs(f.f), push = f.theory * Math.abs(f.f) / Math.max(h1, 1);
      const forces = [[push, 0, color.e, 'wind'], [f.f * v, -f.f * u, color.a, 'Coriolis']];
      for (const [ax, ay, c, name] of forces) {
        const ex = px + ax * reach * k, ey = py - ay * reach * k;
        if (Math.hypot(ex - px, ey - py) < 4) continue;
        arrow(ctx, px, py, ex, ey, { color: c, width: 2, head: 7, dash: [5, 4], open: true });
        const upright = Math.abs(ey - py) > Math.abs(ex - px);
        if (!narrow) label(ctx, P, name, upright ? ex : ex + (ex >= px ? 6 : -6), upright ? ey + (ey > py ? 11 : -11) : ey, { align: upright ? 'center' : ex >= px ? 'left' : 'right', color: c });
      }
    });
  }

  function drawSpiral(ctx, P, narrow) {
    const f = latest, sign = Math.sign(coriolis(latitude)), m = matched(), spiral = m.spiral, V = Math.abs(spiral.surface), p = f.profile;
    const head = title(ctx, P, narrow ? ['the model’s layers (dots)', 'and Ekman’s theory (dashed)'] : ['the model’s layers (dots) and Ekman’s theory (dashed)']);
    const cut = Math.max(1.5 * f.boundaryDepth, 30), tips = [];
    for (let k = 0; k < p.center.length && p.center[k] < cut; k++) tips.push([p.east[k], p.north[k], p.center[k]]);
    const curve = [];
    for (let j = 0; j <= 120; j++) curve.push(spiral.at(-2 * spiral.reversal * j / 120));
    let west = 0, east = 0.35 * V, south = 0, north = 0;
    for (const [u, v] of curve.concat(tips)) { west = Math.min(west, u); east = Math.max(east, u); south = Math.min(south, v); north = Math.max(north, v); }
    const room = { x: P.x + 16, y: P.y + head + 10, w: P.w - 32, h: P.h - head - 22 };
    const k = Math.min(room.w / (east - west), room.h / (north - south)) * 0.92;
    const cx = room.x + room.w / 2 - k * (west + east) / 2, cy = room.y + room.h / 2 + k * (south + north) / 2;
    const sx = (u) => cx + k * u, sy = (v) => cy - k * v;
    clip(ctx, P, () => {
      ctx.strokeStyle = GRID; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(P.x, sy(0)); ctx.lineTo(P.x + P.w, sy(0)); ctx.moveTo(sx(0), room.y - 6); ctx.lineTo(sx(0), P.y + P.h); ctx.stroke();
      arrow(ctx, sx(0), sy(0), sx(0.35 * V), sy(0), { color: color.e, width: 2.5, head: 8 });
      label(ctx, P, 'wind', sx(0.35 * V) + 4, sy(0) + sign * 9, { color: color.e });
      ctx.strokeStyle = 'rgba(255,255,255,0.45)'; ctx.lineWidth = 1; ctx.setLineDash([3, 3]);
      ctx.beginPath(); curve.forEach(([u, v], j) => { if (j) ctx.lineTo(sx(u), sy(v)); else ctx.moveTo(sx(u), sy(v)); }); ctx.stroke(); ctx.setLineDash([]);
      for (let j = 9; j >= 0; j--) {
        const [u, v] = spiral.at(-j * spiral.reversal / 6);
        arrow(ctx, sx(0), sy(0), sx(u), sy(v), { color: `rgba(255,255,255,${(0.7 - 0.05 * j).toFixed(2)})`, width: 1, head: 4, dash: [3, 3] });
      }
      if (tips.length) {
        ctx.strokeStyle = ACCENT; ctx.lineWidth = 1.5; ctx.lineJoin = 'round';
        ctx.beginPath(); tips.forEach(([u, v], j) => { if (j) ctx.lineTo(sx(u), sy(v)); else ctx.moveTo(sx(u), sy(v)); }); ctx.stroke();
        arrow(ctx, sx(0), sy(0), sx(tips[0][0]), sy(tips[0][1]), { color: ACCENT, width: 2.2, head: 7 });
        ctx.fillStyle = ACCENT;
        for (const [u, v] of tips) { ctx.beginPath(); ctx.arc(sx(u), sy(v), 2.6, 0, 2 * Math.PI); ctx.fill(); }
        for (const target of [0, 10, 20, 30]) {
          let best = 0;
          tips.forEach(([, , z], j) => { if (Math.abs(z - target) < Math.abs(tips[best][2] - target)) best = j; });
          const [u, v, z] = tips[best], r = Math.hypot(u, v) || 1;
          label(ctx, P, target === 0 ? `${z.toFixed(0)} m, the top layer` : `${z.toFixed(0)} m`, sx(u) + 7 * u / r, sy(v) - 7 * v / r + (target === 0 ? sign * 4 : 0), { align: u >= 0 ? 'left' : 'right', color: ACCENT });
        }
      }
    });
  }

  function polygons() {
    if (cells) return cells;
    const c = Math.cos(CENTER_LONGITUDE), s = Math.sin(CENTER_LONGITUDE);
    const project = (p) => [-s * p.x + c * p.y, p.z, c * p.x + s * p.y];
    cells = [...new Grid(N)].map((cell, i) => ({ i, center: project(cell.centerVertex), vertices: cell.vertices.map(project), position: [cell.centerVertex.x, cell.centerVertex.y, cell.centerVertex.z] }));
    cells.sort((a, b) => a.center[2] - b.center[2]);
    return cells;
  }

  function drawPlanet(ctx, w, h) {
    const f = latest, strip = clamp(Math.round(w * 0.17), 58, 110), top = 36, bottom = 34;
    const R = Math.min((w - strip - 14) / 2, (h - top - bottom) / 2), cx = strip + (w - strip) / 2, cy = top + (h - top - bottom) / 2;
    const yOf = (lat) => cy - R * Math.sin(lat * Math.PI / 180);
    const x0 = Math.round(strip / 2) + 0.5, half = strip / 2 - 12, peak = Math.max(1e-9, ...Array.from(f.zonal.stress, Math.abs));
    ctx.strokeStyle = GRID; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(x0, yOf(90)); ctx.lineTo(x0, yOf(-90)); ctx.stroke();
    for (const lat of [-60, -30, 0, 30, 60]) { ctx.beginPath(); ctx.moveTo(x0 - half, yOf(lat)); ctx.lineTo(x0 + half, yOf(lat)); ctx.stroke(); text(ctx, lat ? `${Math.abs(lat)}°${lat > 0 ? 'N' : 'S'}` : '0°', 6, yOf(lat), { color: MUTED, size: 9.5, halo: HALO }); }
    ctx.strokeStyle = color.e; ctx.lineWidth = 1.5;
    ctx.beginPath();
    f.zonal.latitude.forEach((lat, b) => { const x = x0 + half * f.zonal.stress[b] / peak, y = yOf(lat); if (b) ctx.lineTo(x, y); else ctx.moveTo(x, y); });
    ctx.stroke();
    for (let b = 1; b < f.zonal.latitude.length; b += 2) {
      const x = x0 + half * f.zonal.stress[b] / peak, y = yOf(f.zonal.latitude[b]);
      if (Math.abs(x - x0) > 4) arrow(ctx, x0, y, x, y, { color: color.e, width: 1.5, head: 5 });
    }
    text(ctx, 'wind', x0, 14, { align: 'center', color: color.e, size: 11 });
    text(ctx, 'east →', x0, 27, { align: 'center', color: MUTED, size: 9.5 });

    const list = polygons(), rgb = [0, 0, 0], S = scale || 1;
    ctx.save();
    ctx.beginPath(); ctx.arc(cx, cy, R, 0, 2 * Math.PI); ctx.fillStyle = DARK; ctx.fill(); ctx.clip();
    ctx.lineWidth = 0.6; ctx.strokeStyle = 'rgba(0,0,0,0.35)'; ctx.lineJoin = 'round';
    for (const cell of list) {
      if (cell.center[2] < -0.05) continue;
      ctx.fillStyle = css(rampRGB(0.5 - 0.5 * clamp(f.upwelling[cell.i] * DAY / S, -1, 1), rgb, DARK_NEUTRAL));
      ctx.beginPath();
      cell.vertices.forEach(([x, y], j) => { if (j) ctx.lineTo(cx + R * x, cy - R * y); else ctx.moveTo(cx + R * x, cy - R * y); });
      ctx.closePath(); ctx.fill(); ctx.stroke();
    }
    const reference = tidy(niceStep(tau / (RHO * coriolis(30)))), k = 0.12 * R / reference, longest = 2 * reference;
    const c = Math.cos(CENTER_LONGITUDE), s = Math.sin(CENTER_LONGITUDE);
    for (const cellInfo of list) {
      if (cellInfo.center[2] < 0.2) continue;
      const i = cellInfo.i, [x, y, z] = cellInfo.position, r = Math.hypot(x, y);
      if (r < 1e-9) continue;
      let east = f.ekmanEast[i], north = f.ekmanNorth[i];
      const m = Math.hypot(east, north);
      if (m > longest) { east *= longest / m; north *= longest / m; }
      const e = [-y / r, x / r, 0], nn = [-z * x / r, -z * y / r, r];
      const v = [east * e[0] + north * nn[0], east * e[1] + north * nn[1], east * e[2] + north * nn[2]];
      const dx = (-s * v[0] + c * v[1]) * k, dy = v[2] * k, px = cx + R * cellInfo.center[0], py = cy - R * cellInfo.center[1];
      if (Math.hypot(dx, dy) > 3) arrow(ctx, px, py, px + dx, py - dy, { color: INK, width: 1.2, head: 4.5 });
    }
    const y = yOf(latitude), chord = R * Math.cos(latitude * Math.PI / 180);
    ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]);
    ctx.beginPath(); ctx.moveTo(cx - chord, y); ctx.lineTo(cx + chord, y); ctx.stroke(); ctx.setLineDash([]);
    ctx.restore();
    ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(cx, cy, R, 0, 2 * Math.PI); ctx.stroke();
    ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]);
    ctx.beginPath(); ctx.moveTo(x0 - half, y); ctx.lineTo(x0 + half, y); ctx.stroke(); ctx.setLineDash([]);
    text(ctx, `day ${f.day.toFixed(1)}`, w - 10, 14, { align: 'right', color: MUTED, size: 11 });
    const ax = w - 14 - k * reference;
    arrow(ctx, ax, h - 14, ax + k * reference, h - 14, { color: INK, width: 1.2, head: 4.5 });
    text(ctx, `${reference} m²/s`, ax - 6, h - 14, { align: 'right', color: MUTED, size: 10 });
  }

  restart();
}
