import { Figure, slider, choice, buttons, readout, legend, caption, text, clamp, rampRGB, paletteVersion, INK, MUTED, GRID, LINE, ACCENT } from '../runtime.module.js';
import { seasonPhrase, MARKERS } from '../insolation.module.js';
import { slabCycle } from '../oceanColumn.module.js';
import { sunSwing, startingSea, airWarming, RHO_CP, SLAB, STORE, DAMPING, POLEWARD, YEAR_DAYS } from '../seaStore.module.js';

const WIDE = 560, DEPTH = 450, ROWS = 225, SPAN = 6, SPEED = 30, AIR = [-50, 150], HALO = 'rgba(20, 20, 22, 0.85)';
const SLAB_DEPTH = SLAB / RHO_CP, SNAP = 4, MODE_DEPTH = 5, TOP_LAYER = 3;
const SHORT = ['Mar', 'Jun', 'Sep', 'Dec', 'Mar'], LONG = ['March equinox', 'June solstice', 'Sept. equinox', 'Dec. solstice', 'March'];
const signed = (v, d) => (Math.abs(v) < 0.5 * 10 ** -d ? (0).toFixed(d) : `${v < 0 ? '−' : '+'}${Math.abs(v).toFixed(d)}`);
const celsius = (k, d = 1) => `${(k - 273.15).toFixed(d).replace('-', '−')} °C`;
const rgb = (f) => `rgb(${rampRGB(f).map((c) => Math.round(c * 255)).join(', ')})`;

export function mountSeastore(root) {
  const controls = root.querySelector('.controls'), field = document.createElement('canvas'), paint = field.getContext('2d');
  const settings = { latitude: POLEWARD, heating: Math.round(sunSwing(POLEWARD)), wind: 0.1 };
  let compare = true, worker = null, run = 0, data = null, day = 0, reached = 0, running = true, geometry = null, painted = 0, paintedPalette = -1, keyed = -1, dragging = false, failed = null;

  const fig = new Figure(root, { height: (w) => layout(w).h, step, draw });
  caption(root, 'One column of the model’s own ocean, its 58 layers, through a year from the March equinox, stepped every 3 hours from the model’s starting sea at the chosen latitude. The heating is prescribed for this figure, not computed by the climate model: it swings with the sunlight a sea absorbs at that latitude under a sky that lets 70 percent of it through, strongest at the June solstice, and takes back 10 W/m² for every kelvin the surface is warmer than at the start. Choosing a latitude sets the heating to that latitude’s swing; the slider changes it. The wind’s stress stirs the sea, and its currents are held at rest. Part 7’s slab, a sea 10 m deep, takes the same heating; in the column it is the bar 10 m deep. At the top is the surface temperature through the year, with a dot on each sea’s warmest day so far, and under it the temperature at every depth, the interfaces between the layers and the base of the mixed layer, the depth to which the water is as dense as the surface. The column on the day marked is drawn on its own, the bar at its left showing how far each layer is held at its minimum thickness or follows its isopycnals. Tap the year or drag the marker to pick a day.');
  const key = legend(root, legendItems());
  slider(controls, { label: 'Latitude', min: 0, max: POLEWARD, step: 5, value: settings.latitude, format: (v) => (v ? `${v}°N` : '0°'), onInput: (v) => { settings.latitude = v; settings.heating = Math.round(sunSwing(v)); heatingControl.value = settings.heating; restart(); } });
  const heatingControl = slider(controls, { label: 'Seasonal heating', min: 0, max: 180, step: 1, value: settings.heating, format: (v) => `±${v} W/m²`, onInput: (v) => { settings.heating = v; restart(); } });
  slider(controls, { label: 'Wind stress', min: 0, max: 0.2, step: 0.01, value: settings.wind, format: (v) => `${v.toFixed(2)} N/m²`, onInput: (v) => { settings.wind = v; restart(); } });
  choice(controls, { label: 'Compare with', options: [['nothing', 'none'], ['Part 7’s 10 m slab', 'slab']], value: 'slab', onChange: (v) => { compare = v === 'slab'; key.set(legendItems()); report(); fig.render(); } });
  const [, toggle] = buttons(controls, [['Restart', () => { day = reached = 0; resume(); }], ['Pause', () => { if (running) halt(); else { if (data && day >= YEAR_DAYS) day = reached = 0; resume(); } }]]);
  const out = readout(controls), heatOut = readout(controls), swingOut = readout(controls);

  function legendItems() {
    const mean = startingSea(settings.latitude).temperature;
    return [
      ['ramp', `water temperature, ${celsius(mean - SPAN, 0)} to ${celsius(mean + SPAN, 0)}`, rgb(0), rgb(1), rgb(0.5)],
      ['line', 'the surface temperature and the mixed layer’s base', ACCENT],
      ['line', 'the column’s temperature on the day marked', INK],
      ['dash', 'the column at the start', 'rgba(255, 255, 255, 0.45)'],
      ['gradient', `heat the top ${STORE} m has gained or lost`, rgb(0.08), rgb(0.92)],
      ...(compare ? [['dash', 'Part 7’s slab, 10 m deep', INK]] : []),
      ['gradient', 'a layer: at its minimum thickness (fixed-depth) to following its isopycnals', 'rgba(255, 255, 255, 0.1)', 'rgba(255, 255, 255, 0.8)'],
    ];
  }

  function start() {
    worker ??= new Worker(new URL('../seastore.worker.js', import.meta.url), { type: 'module' });
    worker.onmessage = ({ data: message }) => {
      if (message.run !== run) return;
      if (message.type === 'layout') receiveLayout(message);
      else if (message.type === 'samples') receive(message);
      else if (message.type === 'error') { failed = message.message; console.error(`seastore: ${failed}`); }
      if (!fig.frame) fig.render();
    };
    worker.onerror = (event) => { event.preventDefault(); failed = event.message || 'the worker could not start'; report(); fig.render(); };
    data = null; painted = 0; failed = null;
    worker.postMessage({ type: 'start', run: ++run, settings: { ...settings } });
  }

  function receiveLayout(layout) {
    const n = layout.samples, L = layout.layers;
    data = {
      layout, received: 0, done: false,
      time: new Float64Array(n), sst: new Float32Array(n), slab: new Float32Array(n), depth: new Float32Array(n), boundary: new Float32Array(n), salinity: new Float32Array(n), stored: new Float64Array(n),
      bases: new Float32Array(n * L), temperatures: new Float32Array(n * L), range: null,
    };
    field.width = n; field.height = ROWS;
    painted = 0;
  }

  function receive(batch) {
    if (!data) return;
    const { layout: { layers: L } } = data, from = batch.from;
    for (const name of ['time', 'sst', 'slab', 'depth', 'boundary', 'salinity', 'stored']) data[name].set(batch[name], from);
    data.bases.set(batch.bases, from * L); data.temperatures.set(batch.temperatures, from * L);
    data.received = from + batch.time.length;
    data.done = batch.done;
    data.range ??= ranges();
  }

  function ranges() {
    const { layout, bases, temperatures } = data, L = layout.layers, mean = layout.mean;
    const reach = (capacity) => slabCycle({ capacity, heating: layout.heating, damping: DAMPING }).amplitude + 1.5;
    let coldest = mean;
    for (let k = 0, top = 0; k < L && top < DEPTH; k++) { if (bases[k] - top > 0.1) coldest = Math.min(coldest, temperatures[k]); top = bases[k]; }
    const around = (r) => [Math.floor(mean - 273.15 - Math.max(r, 3)), Math.ceil(mean - 273.15 + Math.max(r, 3))];
    return { slab: around(reach(SLAB)), model: around(reach(20 * RHO_CP)), coldest: Math.floor(coldest - 273.15 - 0.5) };
  }

  function restart() { day = reached = 0; key.set(legendItems()); start(); resume(); }
  function resume() { running = true; toggle.textContent = 'Pause'; fig.play(true); fig.render(); }
  function halt() { running = false; toggle.textContent = 'Play'; fig.play(false); fig.render(); }

  function step(dt) {
    if (!worker) start();
    if (!running || !data || !data.received) return;
    const last = data.received - 1;
    day = Math.min(day + dt * SPEED, last);
    reached = Math.max(reached, day);
    if (data.done && day >= YEAR_DAYS) { day = YEAR_DAYS; running = false; toggle.textContent = 'Play'; fig.play(false); }
  }

  function layout(w) {
    if (w >= WIDE) {
      const colL = 46, colR = colL + Math.round((w - 110) * 0.36), yearL = colR + 24, yearR = w - 44;
      const stripTop = 32, stripBottom = 128, top = 162, bottom = 392;
      return { w, h: bottom + 34, colL, colR, yearL, yearR, stripTop, stripBottom, top, bottom, columnTop: top, columnBottom: bottom, gauge: { x0: colL - 30, x1: colR, y: 74 } };
    }
    const yearL = 12, yearR = w - 40, stripTop = 28, stripBottom = 108, top = 132, bottom = 282;
    const gaugeY = bottom + 70, columnTop = gaugeY + 80, columnBottom = columnTop + 180;
    return { w, h: columnBottom + 14, colL: 44, colR: w - 12, yearL, yearR, stripTop, stripBottom, top, bottom, columnTop, columnBottom, gauge: { x0: 14, x1: w - 14, y: gaugeY } };
  }

  function shade() {
    const L = data.layout.layers, mean = data.layout.mean, upTo = data.received;
    if (paintedPalette !== paletteVersion) { painted = 0; paintedPalette = paletteVersion; }
    if (painted >= upTo) return;
    const image = paint.createImageData(upTo - painted, ROWS), pixels = image.data, color = [0, 0, 0];
    for (let j = painted; j < upTo; j++) {
      let k = 0;
      for (let r = 0; r < ROWS; r++) {
        const z = (r + 0.5) * DEPTH / ROWS;
        while (k < L - 1 && data.bases[j * L + k] <= z) k++;
        rampRGB(0.5 + (data.temperatures[j * L + k] - mean) / (2 * SPAN), color);
        const p = 4 * (r * (upTo - painted) + j - painted);
        pixels[p] = color[0] * 255; pixels[p + 1] = color[1] * 255; pixels[p + 2] = color[2] * 255; pixels[p + 3] = 255;
      }
    }
    paint.putImageData(image, painted, 0);
    painted = upTo;
  }

  function draw(ctx, w) {
    if (keyed !== paletteVersion) { keyed = paletteVersion; key.set(legendItems()); }
    const g = geometry = layout(w), { yearL, yearR, stripTop, stripBottom, top, bottom } = g;
    const X = (d) => yearL + d / YEAR_DAYS * (yearR - yearL), Y = (z) => top + Math.min(z, DEPTH) / DEPTH * (bottom - top);
    g.X = X;
    text(ctx, 'surface temperature, °C', yearL, stripTop - 14, { color: MUTED, size: 10 });
    text(ctx, 'temperature with depth, and the mixed layer’s base', yearL, top - 12, { color: MUTED, size: 10 });
    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    ctx.strokeRect(yearL + 0.5, stripTop + 0.5, yearR - yearL - 1, stripBottom - stripTop - 1);
    ctx.strokeRect(yearL + 0.5, top + 0.5, yearR - yearL - 1, bottom - top - 1);
    for (let z = 0; z <= DEPTH; z += 100) {
      if (z > 0 && z < DEPTH) { ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(yearL, Y(z)); ctx.lineTo(yearR, Y(z)); ctx.stroke(); }
      text(ctx, `${z} m`, yearR + 5, clamp(Y(z), top + 5, bottom - 5), { color: MUTED, size: 10 });
    }
    ctx.font = '400 10px system-ui, sans-serif';
    const quarter = (yearR - yearL) / 4, names = ctx.measureText(LONG[0]).width + ctx.measureText(LONG[1]).width / 2 < quarter - 10 && LONG.slice(1, 4).every((s) => ctx.measureText(s).width < quarter - 10) ? LONG : SHORT;
    MARKERS.forEach(([d], n) => {
      ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(X(d), bottom); ctx.lineTo(X(d), bottom + 4); ctx.stroke();
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(X(d), stripTop); ctx.lineTo(X(d), stripBottom); ctx.stroke();
      text(ctx, names[n], X(d), bottom + 14, { align: n === 0 ? 'left' : n === 4 ? 'right' : 'center', color: MUTED, size: 10 });
    });

    if (!data || !data.received) {
      text(ctx, failed ? 'the column could not run' : 'the model is starting', (yearL + yearR) / 2, (top + bottom) / 2, { align: 'center', color: MUTED, size: 12 });
      report();
      return;
    }
    shade();
    const n = Math.min(data.received - 1, Math.round(day)), L = data.layout.layers, upto = Math.min(data.received - 1, Math.max(n, Math.round(reached)));
    const [lo, hi] = compare ? data.range.slab : data.range.model, S = (k) => stripBottom - (k - 273.15 - lo) / (hi - lo) * (stripBottom - stripTop);
    const tick = hi - lo > 30 ? 10 : hi - lo > 8 ? 5 : 2;
    for (let t = Math.ceil(lo / tick) * tick; t <= hi; t += tick) {
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(yearL, S(t + 273.15)); ctx.lineTo(yearR, S(t + 273.15)); ctx.stroke();
      text(ctx, `${t}`.replace('-', '−'), yearR + 5, S(t + 273.15), { color: MUTED, size: 10 });
    }

    ctx.save();
    ctx.beginPath(); ctx.rect(yearL, top, X(Math.min(upto + 0.5, YEAR_DAYS)) - yearL, bottom - top); ctx.clip();
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(field, X(-0.5), top, X(data.layout.samples - 0.5) - X(-0.5), bottom - top);
    ctx.restore();
    ctx.save();
    ctx.beginPath(); ctx.rect(yearL, top, yearR - yearL, bottom - top); ctx.clip();
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.14)'; ctx.lineWidth = 1; ctx.beginPath();
    for (let k = 0; k < L - 1 && data.bases[k] < DEPTH + 30; k++) {
      let low = Infinity, high = -Infinity;
      for (let j = 0; j <= upto; j++) { const z = data.bases[j * L + k]; low = Math.min(low, z); high = Math.max(high, z); }
      if (high - low < 0.3) { ctx.moveTo(X(0), Y(low)); ctx.lineTo(X(upto), Y(low)); continue; }
      for (let j = 0; j <= upto; j++) (j ? ctx.lineTo(X(j), Y(data.bases[j * L + k])) : ctx.moveTo(X(j), Y(data.bases[j * L + k])));
    }
    ctx.stroke();
    for (const [color, width] of [[HALO, 4], [ACCENT, 2]]) {
      ctx.strokeStyle = color; ctx.lineWidth = width; ctx.lineJoin = 'round'; ctx.beginPath();
      for (let j = 0; j <= upto; j++) (j ? ctx.lineTo(X(j), Y(data.depth[j])) : ctx.moveTo(X(j), Y(data.depth[j])));
      ctx.stroke();
    }
    ctx.restore();

    ctx.save();
    ctx.beginPath(); ctx.rect(yearL, stripTop, yearR - yearL, stripBottom - stripTop); ctx.clip();
    const curve = (values, color, dash) => {
      ctx.strokeStyle = color; ctx.lineWidth = 2; ctx.setLineDash(dash ?? []); ctx.beginPath();
      for (let j = 0; j <= upto; j++) (j ? ctx.lineTo(X(j), S(values[j])) : ctx.moveTo(X(j), S(values[j])));
      ctx.stroke(); ctx.setLineDash([]);
    };
    if (compare) curve(data.slab, INK, [5, 3]);
    curve(data.sst, ACCENT);
    for (const [values, color, on] of [[data.slab, INK, compare], [data.sst, ACCENT, true]]) {
      if (!on) continue;
      const w0 = warmest(values, upto);
      if (w0 !== null) { ctx.fillStyle = color; ctx.strokeStyle = HALO; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(X(w0), S(values[w0]), 3, 0, 2 * Math.PI); ctx.fill(); ctx.stroke(); }
    }
    ctx.restore();

    const px = X(n);
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.55)'; ctx.lineWidth = 1; ctx.setLineDash([3, 3]);
    ctx.beginPath(); ctx.moveTo(px, stripTop); ctx.lineTo(px, stripBottom); ctx.moveTo(px, top); ctx.lineTo(px, bottom); ctx.stroke(); ctx.setLineDash([]);
    for (const [v, color, on] of [[data.slab[n], INK, compare], [data.sst[n], ACCENT, true]]) {
      if (!on) continue;
      ctx.fillStyle = color; ctx.strokeStyle = HALO; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(px, clamp(S(v), stripTop, stripBottom), 3.5, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
    }
    for (const [color, width] of [[HALO, 4.5], [INK, 2]]) { ctx.strokeStyle = color; ctx.lineWidth = width; ctx.beginPath(); ctx.arc(px, Y(data.depth[n]), 6, 0, 2 * Math.PI); ctx.stroke(); }

    drawGauge(ctx, g, n);
    drawColumn(ctx, g, n, L);
    report(n);
  }

  function warmest(values, n) {
    let best = 0;
    for (let j = 1; j <= n; j++) if (values[j] > values[best]) best = j;
    return best > 0 && values[n] < values[best] - 0.3 ? best : null;
  }

  function drawGauge(ctx, g, n) {
    const { x0, x1, y } = g.gauge, A = (k) => x0 + (k - AIR[0]) / (AIR[1] - AIR[0]) * (x1 - x0), stored = data.stored[n], air = clamp(airWarming(stored), ...AIR);
    text(ctx, `heat stored in the top ${STORE} m`, x0, y - 30, { color: MUTED, size: 10 });
    ctx.fillStyle = 'rgba(255, 255, 255, 0.06)'; ctx.fillRect(x0, y - 6, x1 - x0, 12);
    ctx.fillStyle = rgb(air >= 0 ? 0.92 : 0.08); ctx.fillRect(Math.min(A(0), A(air)), y - 6, Math.abs(A(air) - A(0)), 12);
    ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.strokeRect(x0 + 0.5, y - 5.5, x1 - x0 - 1, 11);
    for (let k = AIR[0]; k <= AIR[1]; k += 50) {
      ctx.strokeStyle = k ? LINE : INK; ctx.beginPath(); ctx.moveTo(A(k), y + (k ? 6 : -9)); ctx.lineTo(A(k), y + 10); ctx.stroke();
      text(ctx, `${k}`.replace('-', '−'), A(k), y + 18, { align: k === AIR[0] ? 'left' : k === AIR[1] ? 'right' : 'center', color: MUTED, size: 10 });
    }
    text(ctx, 'as kelvins of the whole air column', x0, y + 31, { color: MUTED, size: 10 });
    text(ctx, `${signed(stored / 1e9, 2)} GJ/m²`, clamp(A(air), x0 + 34, x1 - 34), y - 15, { align: 'center', color: '#fff', size: 11, weight: 500 });
  }

  function drawColumn(ctx, g, n, L) {
    const { colL: plotL, colR: plotR, columnTop: top, columnBottom: bottom } = g, { layout } = data, range = compare ? data.range.slab : data.range.model;
    const Y = (z) => top + Math.min(z, DEPTH) / DEPTH * (bottom - top), lo = Math.min(data.range.coldest, range[0]), hi = range[1];
    const P = (k) => plotL + (k - 273.15 - lo) / (hi - lo) * (plotR - plotL);
    text(ctx, `the column on day ${n}`, plotL, top - 26, { color: MUTED, size: 10 });
    text(ctx, '°C', plotR, top - 26, { align: 'right', color: MUTED, size: 10 });
    for (let z = 0; z <= DEPTH; z += 100) {
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(plotL, Y(z)); ctx.lineTo(plotR, Y(z)); ctx.stroke();
      text(ctx, `${z} m`, plotL - 5, clamp(Y(z), top + 5, bottom - 5), { align: 'right', color: MUTED, size: 10 });
    }
    const tick = hi - lo > 16 ? 5 : 2;
    for (let t = Math.ceil(lo / tick) * tick; t <= hi; t += tick) {
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(P(t + 273.15), top); ctx.lineTo(P(t + 273.15), bottom); ctx.stroke();
      text(ctx, `${t}`.replace('-', '−'), P(t + 273.15), top - 9, { align: 'center', color: MUTED, size: 10 });
    }

    const profile = (j) => {
      const out = [];
      for (let k = 0, z = 0; k < L && z < DEPTH; k++) {
        const base = data.bases[j * L + k];
        if (base - z > 0.05) out.push({ k, top: z, bottom: base, t: data.temperatures[j * L + k], mode: clamp((base - z - layout.minimum[k]) / MODE_DEPTH, 0, 1) });
        z = base;
      }
      return out;
    };
    const now = profile(n), then = profile(0), at = (layers, z) => (layers.find((l) => z < l.bottom) ?? layers[layers.length - 1]).t;
    ctx.save();
    ctx.beginPath(); ctx.rect(plotL, top, plotR - plotL, bottom - top); ctx.clip();
    ctx.fillStyle = 'rgba(255, 232, 160, 0.08)'; ctx.fillRect(plotL, top, plotR - plotL, Y(data.depth[n]) - top);
    for (let y = top; y < Y(STORE); y++) {
      const z = (y + 0.5 - top) / (bottom - top) * DEPTH, a = P(at(then, z)), b = P(at(now, z));
      if (Math.abs(b - a) < 0.5) continue;
      ctx.fillStyle = rgb(b > a ? 0.92 : 0.08); ctx.globalAlpha = 0.45; ctx.fillRect(Math.min(a, b), y, Math.abs(b - a), 1); ctx.globalAlpha = 1;
    }
    ctx.strokeStyle = GRID; ctx.setLineDash([1, 3]); ctx.beginPath(); ctx.moveTo(plotL, Y(STORE)); ctx.lineTo(plotR, Y(STORE)); ctx.stroke(); ctx.setLineDash([]);
    let labelBottom = Y(data.depth[n]) + 4;
    for (const layer of now.slice(1)) {
      const y0 = Y(layer.top), y1 = Math.min(Y(layer.bottom), bottom);
      ctx.strokeStyle = 'rgba(255, 255, 255, 0.2)'; ctx.beginPath(); ctx.moveTo(plotL, y0); ctx.lineTo(plotR, y0); ctx.stroke();
      const mid = (y0 + y1) / 2;
      if (y1 - y0 >= 14 && mid - 6 > labelBottom) {
        text(ctx, `layer ${layer.k}, ${(layer.bottom - layer.top).toFixed(0)} m`, plotR - 5, mid, { align: 'right', color: MUTED, size: 10, halo: HALO });
        labelBottom = mid + 6;
      }
    }
    for (const layer of now) {
      ctx.fillStyle = `rgba(255, 255, 255, ${0.1 + 0.7 * layer.mode})`;
      ctx.fillRect(plotL + 1, Y(layer.top) + 0.5, 5, Math.max(0.5, Y(layer.bottom) - Y(layer.top) - 1));
    }
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.35)'; ctx.lineWidth = 1; ctx.setLineDash([2, 2]); ctx.beginPath();
    then.forEach((layer, i) => { const x = P(layer.t); if (i) ctx.lineTo(x, Y(layer.top)); else ctx.moveTo(x, Y(layer.top)); ctx.lineTo(x, Y(layer.bottom)); });
    ctx.stroke(); ctx.setLineDash([]);
    for (const [stroke, width] of [[HALO, 4], [INK, 2]]) {
      ctx.strokeStyle = stroke; ctx.lineWidth = width; ctx.lineJoin = 'round'; ctx.beginPath();
      now.forEach((layer, i) => { const x = P(layer.t); if (i) ctx.lineTo(x, Y(layer.top)); else ctx.moveTo(x, Y(layer.top)); ctx.lineTo(x, Y(layer.bottom)); });
      ctx.stroke();
    }
    const base = Y(data.depth[n]);
    for (const [stroke, width] of [[HALO, 4], [ACCENT, 2]]) { ctx.strokeStyle = stroke; ctx.lineWidth = width; ctx.beginPath(); ctx.moveTo(plotL, base); ctx.lineTo(plotR, base); ctx.stroke(); }
    text(ctx, `mixed layer, ${data.depth[n].toFixed(0)} m`, plotL + 5, clamp(base - 9, top + 8, bottom - 8), { color: ACCENT, size: 10, halo: HALO });
    if (compare) {
      const sx = P(data.slab[n]), sy = Y(SLAB_DEPTH);
      ctx.fillStyle = HALO; ctx.fillRect(sx - 3.5, top, 7, sy - top + 1.5);
      ctx.fillStyle = INK; ctx.fillRect(sx - 2, top, 4, sy - top);
      const left = sx - plotL > plotR - sx;
      text(ctx, 'slab, 10 m', left ? sx - 7 : sx + 7, top + 7, { align: left ? 'right' : 'left', color: INK, size: 10, halo: HALO });
    }
    ctx.restore();
    ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.strokeRect(plotL + 0.5, top + 0.5, plotR - plotL - 1, bottom - top - 1);
  }

  function report(n = null) {
    if (n === null || !data) { out.set([['day', '0, the March equinox'], ['the model', failed ? `could not start: ${failed}` : 'starting']]); heatOut.set([]); swingOut.set([]); return; }
    const air = airWarming(data.stored[n]), phrase = seasonPhrase(n), every = running ? 3 : 1, { layers: L, minimum } = data.layout;
    let held = 0, following = 0;
    for (let k = 0, z = 0; k < L && z < DEPTH; k++) {
      const base = data.bases[n * L + k];
      if (clamp((base - z - minimum[k]) / MODE_DEPTH, 0, 1) < 0.5) held++; else following++;
      z = base;
    }
    out.set([
      ['day', `${n}, ${phrase}`],
      ['mixed layer', `${data.depth[n].toFixed(0)} m deep, the wind’s boundary layer ${data.boundary[n].toFixed(0)} m`],
      [`top layer, ${TOP_LAYER} m`, `${celsius(data.sst[n])}, ${data.salinity[n].toFixed(2)} psu`],
      [`layers in the top ${DEPTH} m`, `${held} at their minimum thickness, ${following} following their isopycnals`],
      ...(compare ? [['slab', celsius(data.slab[n])]] : []),
    ], every);
    heatOut.set([
      [`heat stored in the top ${STORE} m since the year began`, `${signed(data.stored[n] / 1e9, 2)} GJ/m²`],
      ['as warming of the whole air column above', `${signed(air, 0)} K`],
    ], every);
    const swing = (values) => { let low = Infinity, high = -Infinity; for (let j = 0; j <= n; j++) { low = Math.min(low, values[j]); high = Math.max(high, values[j]); } return high - low; };
    const peak = (values) => { const w0 = warmest(values, n); return w0 === null ? '' : `, warmest on day ${w0}`; };
    swingOut.set([
      [n >= YEAR_DAYS ? 'over the year the surface swung by' : 'so far the surface has swung by', `${swing(data.sst).toFixed(1)} K${peak(data.sst)}`],
      ...(compare ? [['the slab by', `${swing(data.slab).toFixed(1)} K${peak(data.slab)}`]] : []),
    ], every);
  }

  function aim(p) {
    if (!data || !data.received) return;
    const { yearL, yearR } = geometry;
    const raw = (p.x - yearL) / (yearR - yearL) * YEAR_DAYS, snapped = MARKERS.slice(0, 4).map(([m]) => m).find((m) => Math.abs(m - raw) * (yearR - yearL) / YEAR_DAYS < SNAP);
    day = clamp(Math.round(snapped ?? raw), 0, data.received - 1);
    if (running) halt(); else fig.render();
  }
  const at = (e) => { const r = fig.canvas.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
  const inYear = (p) => !!geometry && p.x >= geometry.yearL - 6 && p.x <= geometry.yearR + 6 && p.y >= geometry.stripTop - 6 && p.y <= geometry.bottom + 6;
  const near = (p) => inYear(p) && !!data && Math.abs(p.x - geometry.X(Math.round(day))) < 16;
  fig.pointer({ hit: near, down: (p) => { dragging = true; aim(p); }, move: (p) => { if (dragging) aim(p); }, up: () => { dragging = false; } });
  fig.canvas.addEventListener('click', (e) => { const p = at(e); if (inYear(p) && !near(p)) aim(p); });
  fig.canvas.addEventListener('pointermove', (e) => { if (e.pointerType !== 'mouse' || dragging) return; const p = at(e); fig.canvas.style.cursor = near(p) ? 'grab' : inYear(p) ? 'crosshair' : ''; });

  report();
  fig.play(true);
}
