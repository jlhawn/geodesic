import { Figure, slider, legend, caption, readout, text, sequentialRGB, termColor, paletteVersion, clamp, INK, MUTED, GRID, LINE } from '../runtime.module.js';
import { insolationMap, daylightHours, noonElevation, subsolarLatitude, seasonPhrase, DAYS, MARKERS } from '../insolation.module.js';

const TOP = 550, LEVELS = [100, 200, 300, 400, 500], TAGGED = 183, AXIS = [0, 200, 400], SIDE = 600, SNAP = 7, HALO = 'rgba(20, 20, 22, 0.85)';
const SHORT = ['Mar', 'Jun', 'Sep', 'Dec', 'Mar'], LONG = ['March equinox', 'June solstice', 'Sept. equinox', 'Dec. solstice', 'March'];
const place = (lat) => (lat === 90 ? 'the North Pole' : lat === -90 ? 'the South Pole' : lat === 0 ? 'the equator' : `${Math.abs(lat)}°${lat > 0 ? 'N' : 'S'}`);
const layout = (w) => {
  const x0 = 44, y0 = 30, gap = 16, sw = clamp(Math.round(w * 0.25), 86, 160), mw = w - x0 - gap - sw - 12, mh = clamp(Math.round(w * 0.46), 230, 290);
  return { x0, y0, mw, mh, sx: x0 + mw + gap, sw, h: y0 + mh + 26 };
};

function contours(map) {
  const { rows, latitude, daily } = map, cols = DAYS + 1, value = (i, j) => daily[i * DAYS + (j % DAYS)];
  return LEVELS.map((level) => {
    const segments = [];
    for (let i = 0; i + 1 < rows; i++) for (let j = 0; j + 1 < cols; j++) {
      const corners = [[i, j], [i, j + 1], [i + 1, j + 1], [i + 1, j]], points = [];
      for (let k = 0; k < 4; k++) {
        const [a, b] = corners[k], [c, d] = corners[(k + 1) % 4], p = value(a, b), q = value(c, d);
        if ((p - level) * (q - level) < 0) { const t = (level - p) / (q - p); points.push(b + t * (d - b), latitude[a] + t * (latitude[c] - latitude[a])); }
      }
      if (points.length >= 4) segments.push(...points.slice(0, 4));
      if (points.length === 8) segments.push(...points.slice(4, 8));
    }
    return Float32Array.from(segments);
  });
}

export function mountSeasons(root) {
  const controls = root.querySelector('.controls'), sun = termColor('a'), field = document.createElement('canvas'), paint = field.getContext('2d');
  const map = insolationMap(), lines = contours(map), overhead = Array.from({ length: DAYS + 1 }, (_, d) => subsolarLatitude(d)), row = (lat) => lat + 90;
  let lat = 45, day = 91, geometry = null, painted = -1, dragging = false, grab = [0, 0];
  const rgb = (f) => `rgb(${sequentialRGB(f).map((c) => Math.round(c * 255)).join(', ')})`;
  const items = () => [['ramp', `daily-mean sunlight, 0 to ${TOP} W/m²`, rgb(0), rgb(1), rgb(0.5)], ['faint', 'every 100 W/m²', LINE], ['dash', 'where the sun stands overhead at noon', INK], ['line', 'the probe’s day, by latitude', sun], ['line', 'the year’s average, by latitude', INK]];
  let keyed = paletteVersion;

  const fig = new Figure(root, { height: (w) => layout(w).h, draw });
  caption(root, 'The sunlight reaching the top of the air, averaged over each day, at every latitude through the model’s year. The year starts at the March equinox. Each day’s mean follows the model’s own sun, sampled every five minutes. Drag the probe or tap the map to read one latitude on one day. The panel on the right shares the map’s latitudes. It shows the probe’s day and the average over the whole year.');
  const key = legend(root, items());
  const latitudeControl = slider(controls, { label: 'Latitude', min: -90, max: 90, step: 1, value: lat, format: (v) => (v ? `${Math.abs(v)}°${v > 0 ? 'N' : 'S'}` : '0°'), onInput: (v) => { lat = v; update(); } });
  const dayControl = slider(controls, { label: 'Day', min: 0, max: DAYS - 1, step: 1, value: day, format: (v) => `day ${v}`, onInput: (v) => { day = v; update(); } });
  const out = readout(controls);

  function report() {
    const mean = map.daily[row(lat) * DAYS + day], hours = daylightHours(lat, day), noon = noonElevation(lat, day), phrase = seasonPhrase(day);
    out.set([
      ['at', `${place(lat)}${phrase.startsWith('the ') ? ' on ' : ', '}${phrase}`],
      ['daily mean', `${mean.toFixed(1)} W/m²`],
      ['daylight', hours === 0 ? 'none, polar night' : hours === 24 ? 'all 24 hours' : `${hours.toFixed(1)} hours`],
      ['sun at noon', noon >= 89.95 ? 'straight overhead' : noon > 0 ? `${noon.toFixed(1)}° up` : 'below the horizon'],
      ['annual mean', `${map.annual[row(lat)].toFixed(1)} W/m²`],
    ]);
  }

  function update() { latitudeControl.value = lat; dayControl.value = day; report(); fig.render(); }

  function shade() {
    if (painted === paletteVersion) return;
    painted = paletteVersion;
    const cols = DAYS + 1, shades = Array.from({ length: 256 }, (_, i) => sequentialRGB(i / 255).map((c) => Math.round(c * 255)));
    field.width = cols; field.height = map.rows;
    const img = paint.createImageData(cols, map.rows), data = img.data;
    for (let i = 0; i < map.rows; i++) for (let j = 0; j < cols; j++) {
      const c = shades[clamp(Math.round(map.daily[i * DAYS + (j % DAYS)] / TOP * 255), 0, 255)], n = 4 * ((map.rows - 1 - i) * cols + j);
      data[n] = c[0]; data[n + 1] = c[1]; data[n + 2] = c[2]; data[n + 3] = 255;
    }
    paint.putImageData(img, 0, 0);
  }

  function draw(ctx, w) {
    if (keyed !== paletteVersion) { keyed = paletteVersion; key.set(items()); }
    shade();
    const { x0, y0, mw, mh, sx, sw } = layout(w), bottom = y0 + mh;
    const X = (d) => x0 + d / DAYS * mw, Y = (l) => y0 + (90 - l) / 180 * mh, S = (v) => sx + v / SIDE * sw;
    geometry = { x0, y0, mw, mh, sx, sw, X, Y };

    ctx.save();
    ctx.beginPath(); ctx.rect(x0, y0, mw, mh); ctx.clip();
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(field, X(-0.5), Y(90.5), X(DAYS + 0.5) - X(-0.5), Y(-90.5) - Y(90.5));
    ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.beginPath();
    for (const s of lines) for (let k = 0; k < s.length; k += 4) { ctx.moveTo(X(s[k]), Y(s[k + 1])); ctx.lineTo(X(s[k + 2]), Y(s[k + 3])); }
    ctx.stroke();
    for (const [color, width] of [[HALO, 3.5], [INK, 1.5]]) {
      ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(color === INK ? [5, 4] : []); ctx.beginPath();
      overhead.forEach((l, d) => (d ? ctx.lineTo(X(d), Y(l)) : ctx.moveTo(X(d), Y(l))));
      ctx.stroke();
    }
    ctx.setLineDash([]);
    for (const level of LEVELS.slice(0, 4)) {
      let at = null;
      for (let l = 0; l < 90 && at === null; l++) { const a = map.daily[row(l) * DAYS + TAGGED], b = map.daily[row(l + 1) * DAYS + TAGGED]; if ((a - level) * (b - level) <= 0 && a !== b) at = l + (a - level) / (a - b); }
      if (at !== null) text(ctx, `${level}`, X(TAGGED), Y(at), { align: 'center', color: INK, size: 10, halo: HALO });
    }
    const px = X(day), py = Y(lat);
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.28)'; ctx.lineWidth = 1; ctx.setLineDash([3, 3]);
    ctx.beginPath(); ctx.moveTo(px, y0); ctx.lineTo(px, bottom); ctx.stroke();
    ctx.setLineDash([]);
    ctx.restore();
    ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.strokeRect(x0 + 0.5, y0 + 0.5, mw - 1, mh - 1);
    for (const [color, width] of [[HALO, 4.5], [INK, 2]]) {
      ctx.strokeStyle = color; ctx.lineWidth = width; ctx.beginPath(); ctx.arc(px, py, 8, 0, 2 * Math.PI);
      for (const [dx, dy] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) if (px + dx * 15 >= x0 && px + dx * 15 <= x0 + mw && py + dy * 15 >= y0 && py + dy * 15 <= bottom) { ctx.moveTo(px + dx * 11, py + dy * 11); ctx.lineTo(px + dx * 15, py + dy * 15); }
      ctx.stroke();
    }

    for (const l of [90, 60, 30, 0, -30, -60, -90]) {
      text(ctx, l ? `${Math.abs(l)}°${l > 0 ? 'N' : 'S'}` : '0°', x0 - 6, clamp(Y(l), y0 + 5, bottom - 5), { align: 'right', color: MUTED, size: 10 });
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(sx, Y(l)); ctx.lineTo(sx + sw, Y(l)); ctx.stroke();
    }
    ctx.font = '400 10px system-ui, sans-serif';
    const names = LONG.every((s, n) => n === 0 || n === 4 || ctx.measureText(s).width < mw / 4 - 10) ? LONG : SHORT;
    MARKERS.forEach(([d], n) => {
      ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(X(d), bottom); ctx.lineTo(X(d), bottom + 4); ctx.stroke();
      text(ctx, names[n], X(d), bottom + 14, { align: n === 0 ? 'left' : n === 4 ? 'right' : 'center', color: MUTED, size: 10 });
    });
    text(ctx, dragging ? `${place(lat)}, day ${day}` : 'drag the probe or tap the map', x0 + mw / 2, 11, { align: 'center', color: MUTED, size: 11 });

    for (const v of AXIS) {
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(S(v), y0); ctx.lineTo(S(v), bottom); ctx.stroke();
      text(ctx, `${v}`, S(v), bottom + 14, { align: v ? 'center' : 'left', color: MUTED, size: 10 });
    }
    text(ctx, 'W/m²', sx + sw, 11, { align: 'right', color: MUTED, size: 11 });
    ctx.strokeStyle = LINE; ctx.strokeRect(sx + 0.5, y0 + 0.5, sw - 1, mh - 1);
    ctx.save();
    ctx.beginPath(); ctx.rect(sx, y0 - 4, sw, mh + 8); ctx.clip();
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.28)'; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(sx, py); ctx.lineTo(sx + sw, py); ctx.stroke(); ctx.setLineDash([]);
    const curve = (value, color) => {
      ctx.strokeStyle = color; ctx.lineWidth = 2; ctx.beginPath();
      for (let i = 0; i < map.rows; i++) { const qx = S(value(i)), qy = Y(map.latitude[i]); if (i) ctx.lineTo(qx, qy); else ctx.moveTo(qx, qy); }
      ctx.stroke();
    };
    curve((i) => map.annual[i], INK);
    curve((i) => map.daily[i * DAYS + day], sun);
    for (const [v, color] of [[map.annual[row(lat)], INK], [map.daily[row(lat) * DAYS + day], sun]]) { ctx.fillStyle = color; ctx.strokeStyle = HALO; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(S(v), py, 3.5, 0, 2 * Math.PI); ctx.fill(); ctx.stroke(); }
    ctx.restore();
  }

  function aim(p, latitudeOnly = false) {
    const { x0, y0, mw, mh } = geometry, snap = (v, targets, scale) => targets.find((t) => Math.abs(t - v) * scale < SNAP) ?? Math.round(v);
    lat = clamp(snap(90 - (p.y - y0) / mh * 180, [-90, -45, 0, 45, 90], mh / 180), -90, 90);
    if (!latitudeOnly) day = clamp(snap((p.x - x0) / mw * DAYS, MARKERS.slice(0, 4).map(([m]) => m), mw / DAYS), 0, DAYS - 1);
    update();
  }
  const at = (e) => { const r = fig.canvas.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
  const near = (p) => !!geometry && Math.hypot(p.x - geometry.X(day), p.y - geometry.Y(lat)) < 22;
  const within = (p, left, width) => p.x >= left && p.x <= left + width && p.y >= geometry.y0 - 6 && p.y <= geometry.y0 + geometry.mh + 6;
  fig.pointer({
    hit: near,
    down: (p) => { dragging = true; grab = [geometry.X(day) - p.x, geometry.Y(lat) - p.y]; fig.render(); },
    move: (p) => { if (dragging) aim({ x: p.x + grab[0], y: p.y + grab[1] }); },
    up: () => { dragging = false; fig.render(); },
  });
  fig.canvas.addEventListener('click', (e) => {
    const p = at(e);
    if (!geometry || near(p)) return;
    if (within(p, geometry.x0, geometry.mw)) aim(p); else if (within(p, geometry.sx, geometry.sw)) aim(p, true);
  });
  fig.canvas.addEventListener('pointermove', (e) => {
    if (dragging || !geometry) return;
    const p = at(e);
    fig.canvas.style.cursor = near(p) ? 'grab' : within(p, geometry.x0, geometry.mw) || within(p, geometry.sx, geometry.sw) ? 'crosshair' : '';
  });

  report();
}
