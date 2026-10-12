import { Figure, slider, choice, buttons, readout, legend, caption, text, arrow, termColor, clamp, INK, MUTED, GRID, LINE } from '../runtime.module.js';
import { seasonPhrase, DAYS } from '../insolation.module.js';
import { seaIceYears, heatAboveFreezing, snowDepth, DEPTHS, DEEP_CLASS, FREEZING_POINT, SINKING_SALINITY, SNOWFALL, YEARS } from '../seaIceYears.module.js';

const WIDE = 520, PACE = 60, HEAD = 16, HALO = 'rgba(20, 20, 22, 0.85)', CHIP = 'rgba(30, 32, 36, 0.88)', PROBE = 'rgba(255, 255, 255, 0.3)', FAINT = 'rgba(255, 255, 255, 0.22)';
const WATER = 'rgb(95, 165, 225)', SALT = 'rgb(214, 190, 140)', ICE = { bright: 'rgb(205, 218, 228)', dark: 'rgb(78, 84, 92)' }, SNOW = 'rgb(246, 248, 251)';
const SEA = 'rgba(55, 115, 175, 0.42)', DEEP = 'rgba(25, 55, 105, 0.6)', SUN_SCALE = 400, ZERO = 273.15;
const MIXED = { slab: DEPTHS.slab, floor: DEPTHS.floor, deep: DEPTHS.deep }, SHARE = { slab: 0.3, floor: 0.5, deep: 0.75 };
const NAMES = ['Sept. equinox', 'Dec. solstice', 'March equinox', 'June solstice', 'Sept. equinox'], SHORT = ['Sep', 'Dec', 'Mar', 'Jun', 'Sep'];
const MARKS = [0, 91, 182, 273, 365];

const zero = (v, d) => Math.abs(v) < 0.5 * 10 ** -d;
const num = (v, d = 1) => `${v < 0 && !zero(v, d) ? '−' : ''}${Math.abs(zero(v, d) ? 0 : v).toFixed(d)}`;
const celsius = (k, d = 1) => `${num(k - ZERO, d)} °C`;
const percent = (f) => `${Math.round(100 * f)}%`;
const niceStep = (span, most) => [0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5, 10, 20, 25, 50].find((s) => span / s <= most) ?? 100;
const decimals = (step) => (step >= 1 ? 0 : step >= 0.1 ? 1 : 2);
const fade = (hex, a) => `rgba(${[1, 3, 5].map((n) => parseInt(hex.slice(n, n + 2), 16)).join(', ')}, ${a})`;
const heat = (joules) => (joules >= 1e9 ? `${(joules / 1e9).toFixed(2)} GJ/m²` : `${Math.round(joules / 1e6)} MJ/m²`);

export function mountSeaice(root) {
  const controls = root.querySelector('.controls'), sun = termColor('a');
  const state = { ice: 'bright', depth: 'floor', convergence: 85, latitude: 75, snow: false };
  let run = null, scales = null, position = 0, playing = true, geometry = null, summaries = new Map(), pending = 0;

  const fig = new Figure(root, { height: (w) => layout(w).h, step, draw });
  caption(root, 'One column of the model’s own sea ice over a mixed layer of seawater, stepped hourly with the model’s ice code for 20 years from open water at 0\u00a0°C on the September equinox. The air above it is a toy written for this figure, not the model’s, and the heat the winds bring is a fixed amount. Snow, when it is on, falls at 100\u00a0kg/m² a year whenever the surface is below 0\u00a0°C. The mixed layer never mixes with the water beneath it and no current brings it heat, so the brine of its ice stays in it until the ice melts. The cartoon draws the ice and snow to a scale of their own, with the leads, the open water between the floes, on the right, and under the mixed layer the 1026.95 class, the densest water the model’s ocean starts with on a planet without land. The gauge places the water’s density against that class: water below the line is the denser, and in the model it would sink. The charts run from one September equinox to the next, with the year being played drawn bright and earlier years faint. Tap the charts to pick a day, or drag the marker.');
  const legendItems = () => [
    ['arrow', 'sunlight reaching the surface', sun], ['arrow', 'sunlight reflected', fade(sun, 0.55)], ['line', 'ice, its volume per square meter of sea', ICE.bright], ...(state.snow ? [['dots', 'snow depth', SNOW]] : []),
    ['line', 'mixed layer’s temperature', WATER], ['line', 'temperature at the surface', INK], ['line', 'mixed layer’s salinity', SALT], ['faint', 'earlier years', 'rgba(255, 255, 255, 0.35)'],
  ];
  const key = legend(root, legendItems());
  choice(controls, { label: 'Ice', options: [['the model’s, bright', 'bright'], ['Part 7’s, dark', 'dark']], value: state.ice, onChange: (v) => change({ ice: v }) });
  choice(controls, { label: 'Mixed layer', options: [['10 m, Part 7’s slab', 'slab'], ['50 m', 'floor'], ['200 m', 'deep']], value: state.depth, onChange: (v) => change({ depth: v }) });
  slider(controls, { label: 'Heat the winds bring', min: 60, max: 110, step: 1, value: state.convergence, format: (v) => `${v} W/m²`, onInput: (v) => change({ convergence: v }) });
  slider(controls, { label: 'Latitude', min: 60, max: 90, step: 1, value: state.latitude, format: (v) => `${v}°N`, onInput: (v) => change({ latitude: v }) });
  choice(controls, { label: 'Snow', options: [['none', 'off'], [`${SNOWFALL} kg/m² a year`, 'on']], value: 'off', onChange: (v) => { change({ snow: v === 'on' }); key.set(legendItems()); } });
  const yearControl = slider(controls, { label: 'Year', min: 1, max: YEARS, step: 1, value: 1, format: (v) => `${v} of ${YEARS}`, onInput: (v) => { position = (v - 1) * DAYS + dayOf(position); halt(); fig.render(); } });
  const [toggle] = buttons(controls, [['Pause', () => (playing ? halt() : resume())], ['Restart', () => { position = 0; resume(); }]]);
  const dayOut = readout(controls), yearOut = readout(controls);

  const yearOf = (p) => Math.min(YEARS - 1, Math.floor(p / DAYS));
  const dayOf = (p) => Math.floor(p) - yearOf(p) * DAYS;

  function compute() {
    pending = 0;
    run = seaIceYears({ latitude: state.latitude, depth: MIXED[state.depth], ice: state.ice, convergence: state.convergence, snow: state.snow });
    summaries = new Map();
    const r = run.record;
    let iceTop = 1, low = FREEZING_POINT, high = FREEZING_POINT + 2, sLow = Math.min(34, SINKING_SALINITY), sHigh = SINKING_SALINITY + 0.1, dLow = DEEP_CLASS - 0.2, dHigh = DEEP_CLASS + 0.1;
    for (let d = 0; d < run.days; d++) {
      iceTop = Math.max(iceTop, r.thickness[d] + snowDepth(r.snow[d]), r.volume[d]);
      low = Math.min(low, r.surface[d]); high = Math.max(high, r.water[d]);
      sLow = Math.min(sLow, r.salinity[d]); sHigh = Math.max(sHigh, r.salinity[d]);
      dLow = Math.min(dLow, r.density[d]); dHigh = Math.max(dHigh, r.density[d]);
    }
    scales = { iceTop, low, high, sLow, sHigh, dLow, dHigh };
    fig.render();
  }
  function change(changes) { Object.assign(state, changes); if (!pending) pending = requestAnimationFrame(compute); }
  function summary(y) { if (!summaries.has(y)) summaries.set(y, run.summary(y)); return summaries.get(y); }

  function resume() { if (position >= run.days - 1) position = 0; playing = true; toggle.textContent = 'Pause'; fig.play(true); }
  function halt() { playing = false; toggle.textContent = position >= run.days - 1 ? 'Play again' : 'Play'; fig.play(false); }

  function step(dt) {
    if (!playing || !run) return false;
    position = Math.min(run.days - 1, position + dt * PACE);
    if (position >= run.days - 1) halt();
    return true;
  }

  function layout(w) {
    const narrow = w < WIDE;
    if (narrow) {
      const cartoon = { x: 0, y: 0, w, h: 270 }, top = cartoon.h + 12, panel = 96, gap = 12;
      return { narrow, cartoon, charts: { x: 40, w: w - 48, top: top + 18, panel, gap }, h: top + 18 + 3 * panel + 2 * gap + 30 };
    }
    const cw = clamp(Math.round(w * 0.36), 200, 250), panel = 108, gap = 14, h = 22 + 3 * panel + 2 * gap + 30;
    return { narrow, cartoon: { x: 0, y: 0, w: cw, h }, charts: { x: cw + 46, w: w - cw - 56, top: 22, panel, gap }, h };
  }

  function draw(ctx, w) {
    if (!run) return;
    const g = geometry = layout(w), d = Math.floor(position), y = yearOf(position);
    drawCartoon(ctx, g, d);
    drawCharts(ctx, g, d, y);
    report(d, y);
  }

  function drawCartoon(ctx, g, d) {
    const { x: x0, y: y0, w: cw, h: ch } = g.cartoon, r = run.record, narrow = g.narrow;
    const gaugeX = x0 + cw - (narrow ? 50 : 44), colX = x0 + (narrow ? 10 : 14), colW = gaugeX - 18 - colX;
    const skyBottom = y0 + ch * 0.22, yW = y0 + ch * 0.52, bottom = y0 + ch - 8, yM = yW + (bottom - yW) * SHARE[state.depth];
    const scale = Math.min(40, (yW - skyBottom) * 0.94 / scales.iceTop);
    const h = r.thickness[d], A = h > 0 ? r.cover[d] : 0, snow = snowDepth(r.snow[d]), iceW = colW * A, light = r.sunlight[d], albedo = r.albedo[d];

    ctx.fillStyle = SEA; ctx.fillRect(colX, yW, colW, yM - yW);
    ctx.fillStyle = DEEP; ctx.fillRect(colX, yM, colW, bottom - yM);
    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(colX, yW + 0.5); ctx.lineTo(colX + colW, yW + 0.5); ctx.stroke();
    ctx.setLineDash([4, 3]); ctx.beginPath(); ctx.moveTo(colX, yM + 0.5); ctx.lineTo(gaugeX, yM + 0.5); ctx.stroke(); ctx.setLineDash([]);
    if (iceW > 0.5) {
      const top = yW - h * scale;
      ctx.fillStyle = ICE[state.ice]; ctx.fillRect(colX, top, iceW, yW - top);
      if (snow > 0) { ctx.fillStyle = SNOW; ctx.fillRect(colX, top - snow * scale, iceW, snow * scale); }
      if (yW - top > 22 && iceW > 24) {
        const c = run.conducted[d], ax = colX + Math.min(iceW / 2, 22), ink = state.ice === 'bright' ? 'rgb(40, 44, 50)' : INK;
        if (Math.abs(c) >= 0.5) arrow(ctx, ax, c > 0 ? yW - 3 : top + 3, ax, c > 0 ? top + 4 : yW - 4, { color: ink, width: 1.5, head: 6 });
        if (iceW > 70 && Math.abs(c) >= 0.5) text(ctx, `${num(Math.abs(c), 0)} W/m² ${c >= 0 ? 'up' : 'down'}`, ax + 8, (top + yW) / 2, { color: ink, size: 10.5 });
      }
    }
    const topY = h > 0 ? yW - (h + snow) * scale : yW, ly = Math.min(topY, yW - 2) - 8;
    const top = h > 0 ? `top ${celsius(run.skin[d])}` : `surface ${celsius(r.surface[d])}`, thick = h > 0 ? `ice ${h.toFixed(2)} m${snow > 0 ? `, snow ${snow.toFixed(2)} m` : ''}` : '';
    ctx.font = '400 10.5px system-ui, sans-serif';
    const crowded = h > 0 && ctx.measureText(top).width + ctx.measureText(thick).width + 12 > colW, topRow = crowded ? ly - 13 : ly;
    const downX = colX + colW * 0.42, upX = colX + colW * 0.6, skyTop = y0 + 30, lowest = Math.max(skyTop + 12, topRow - 9);
    if (light > 0) {
      const lit = 1.2 + 6 * light / SUN_SCALE, back = 1.2 + 6 * light * (albedo || 0) / SUN_SCALE;
      arrow(ctx, downX, skyTop, downX, lowest, { color: sun, width: lit, head: 5 + lit });
      ctx.save(); ctx.globalAlpha = 0.55; arrow(ctx, upX, lowest, upX, skyTop, { color: sun, width: back, head: 5 + back }); ctx.restore();
      text(ctx, `${Math.round(light)} W/m²`, downX - 7, skyTop + 4, { color: sun, size: 10.5, align: 'right', halo: HALO });
      text(ctx, `${percent(albedo)} back`, upX + 7, skyTop + 4, { color: sun, size: 10.5, halo: HALO });
    }
    const sx = colX + 9, sy = y0 + 14, glow = clamp(light / SUN_SCALE, 0, 1);
    ctx.save();
    ctx.strokeStyle = sun; ctx.fillStyle = sun; ctx.globalAlpha = light > 0 ? 0.35 + 0.65 * glow : 0.4;
    ctx.beginPath(); ctx.arc(sx, sy, 5.5, 0, 2 * Math.PI); if (light > 0) ctx.fill(); else ctx.stroke();
    if (light > 0) for (let k = 0; k < 8; k++) { const a = k * Math.PI / 4; ctx.beginPath(); ctx.moveTo(sx + 8 * Math.cos(a), sy + 8 * Math.sin(a)); ctx.lineTo(sx + 11 * Math.cos(a), sy + 11 * Math.sin(a)); ctx.stroke(); }
    ctx.restore();
    if (!(light > 0)) text(ctx, 'polar night', sx + 14, sy, { color: MUTED, size: 10.5 });

    text(ctx, top, colX + colW - 2, topRow, { color: INK, size: 10.5, align: 'right', halo: HALO });
    if (h > 0) text(ctx, thick, colX + 2, ly, { color: INK, size: 10.5, halo: HALO });
    const deep = MIXED[state.depth];
    text(ctx, `mixed layer, ${deep < 20 ? deep.toFixed(1) : Math.round(deep)} m`, colX + 6, yW + 11, { color: INK, size: 10.5 });
    if (yM - yW > 34) text(ctx, `${celsius(r.water[d], 2)} · ${r.salinity[d].toFixed(2)} psu`, colX + 6, yW + 25, { color: INK, size: 10.5 });
    text(ctx, narrow ? 'the deepest class' : 'the deepest class, 1026.95', colX + 6, Math.max(yM + 10, bottom - 9), { color: MUTED, size: 10.5 });

    const rho = r.density[d], above = Math.max(0.05, DEEP_CLASS - scales.dLow), below = Math.max(0.05, scales.dHigh - DEEP_CLASS);
    const k = Math.min((yM - yW - 10) / above, (bottom - yM - 8) / below), Y = (v) => yM + (v - DEEP_CLASS) * k;
    const tick = niceStep(Math.max(above, below) * 2, Math.max(2, (bottom - yW) / 14));
    ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(gaugeX, yW); ctx.lineTo(gaugeX, bottom); ctx.stroke();
    for (let v = Math.ceil((DEEP_CLASS - above) / tick) * tick; v <= DEEP_CLASS + below + 1e-9; v += tick) {
      const ty = Y(v);
      if (ty < yW + 2 || ty > bottom - 1) continue;
      ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(gaugeX, ty); ctx.lineTo(gaugeX + 3, ty); ctx.stroke();
    }
    text(ctx, 'density,', gaugeX - 2, yW - 19, { color: MUTED, size: 9.5 });
    text(ctx, 'kg/m³', gaugeX - 2, yW - 8, { color: MUTED, size: 9.5 });
    const my = clamp(Y(rho), yW + 2, bottom - 2), denser = rho > DEEP_CLASS, apart = Math.abs(my - yM) >= 11;
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.6)'; ctx.beginPath(); ctx.moveTo(gaugeX, yM + 0.5); ctx.lineTo(gaugeX + 5, yM + 0.5); ctx.stroke();
    text(ctx, DEEP_CLASS.toFixed(2), gaugeX + 6, apart ? yM : yM + (my < yM ? 6 : -6), { color: MUTED, size: 9.5 });
    ctx.fillStyle = denser ? SALT : INK;
    ctx.beginPath(); ctx.moveTo(gaugeX - 1, my); ctx.lineTo(gaugeX - 9, my - 4.5); ctx.lineTo(gaugeX - 9, my + 4.5); ctx.closePath(); ctx.fill();
    ctx.strokeStyle = denser ? SALT : 'rgba(255, 255, 255, 0.5)'; ctx.setLineDash([2, 3]); ctx.beginPath(); ctx.moveTo(colX + colW, my); ctx.lineTo(gaugeX - 10, my); ctx.stroke(); ctx.setLineDash([]);
    text(ctx, rho.toFixed(2), gaugeX + 6, apart ? my : my + (my < yM ? -6 : 6), { color: denser ? SALT : INK, size: 9.5, weight: 500 });
    if (denser) {
      const ax = colX + colW * 0.82;
      arrow(ctx, ax, yM - 12, ax, Math.min(bottom - 4, yM + 14), { color: SALT, width: 1.6, head: 6 });
      text(ctx, 'would sink', ax - 6, yM - 8, { color: SALT, size: 10.5, align: 'right', halo: HALO });
    }
  }

  function drawCharts(ctx, g, d, y) {
    const { x, w, top, panel, gap } = g.charts, r = run.record, from = y * DAYS, now = d - from;
    const X = (day) => x + day / (DAYS - 1) * w;
    const panels = [
      { name: 'ice, m', lo: 0, hi: scales.iceTop * 1.05, lines: [['volume', ICE.bright, 2]], extra: state.snow ? [['snow', SNOW]] : [], faint: 'volume' },
      { name: 'temperature, °C', lo: scales.low - ZERO, hi: scales.high - ZERO, lines: [['surface', INK, 1.4], ['water', WATER, 2]], guide: [FREEZING_POINT - ZERO, 'freezing'], faint: 'water', shift: -ZERO },
      { name: 'salinity, psu', lo: scales.sLow, hi: scales.sHigh, lines: [['salinity', SALT, 2]], guide: [SINKING_SALINITY, 'as dense as the deepest class'], faint: 'salinity' },
    ];
    text(ctx, `year ${y + 1}, ${seasonPhrase(r.day[d])}`, x + w / 2, top - 13, { color: MUTED, size: 11, align: 'center' });
    panels.forEach((p, n) => {
      const y0 = top + n * (panel + gap), y1 = y0 + panel, span = p.hi - p.lo || 1, pad = span * 0.06, lo = p.lo - (n ? pad : 0), hi = p.hi + pad;
      const Y = (v) => y1 - 3 - (v - lo) / (hi - lo) * (panel - HEAD - 3), shift = p.shift ?? 0;
      const stepY = niceStep(hi - lo, Math.max(2, Math.floor((panel - HEAD) / 24))), places = decimals(stepY);
      for (let v = Math.ceil(lo / stepY) * stepY; v <= hi + 1e-9; v += stepY) {
        ctx.strokeStyle = GRID; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(x, Y(v)); ctx.lineTo(x + w, Y(v)); ctx.stroke();
        text(ctx, num(v, places), x - 6, Y(v), { color: MUTED, size: 10, align: 'right' });
      }
      ctx.strokeStyle = LINE; ctx.strokeRect(x + 0.5, y0 + 0.5, w - 1, panel - 1);
      const gy = p.guide ? Y(p.guide[0]) : null;
      if (p.guide) { ctx.strokeStyle = LINE; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(x, gy); ctx.lineTo(x + w, gy); ctx.stroke(); ctx.setLineDash([]); }
      ctx.save(); ctx.beginPath(); ctx.rect(x, y0, w, panel); ctx.clip();
      ctx.strokeStyle = FAINT; ctx.lineWidth = 1; ctx.beginPath();
      for (let past = 0; past < y; past++) {
        const series = r[p.faint], start = past * DAYS;
        for (let k = 0; k < DAYS; k++) { const px = X(k), py = Y(series[start + k] + shift); if (k) ctx.lineTo(px, py); else ctx.moveTo(px, py); }
      }
      ctx.stroke();
      for (const [field, color] of p.extra ?? []) {
        ctx.strokeStyle = color; ctx.lineWidth = 2; ctx.setLineDash([1.5, 3]); ctx.lineCap = 'round'; ctx.beginPath();
        for (let k = 0; k <= now; k++) { const px = X(k), py = Y(snowDepth(r[field][from + k])); if (k) ctx.lineTo(px, py); else ctx.moveTo(px, py); }
        ctx.stroke(); ctx.setLineDash([]); ctx.lineCap = 'butt';
      }
      for (const [field, color, width] of p.lines) {
        ctx.strokeStyle = color; ctx.lineWidth = width; ctx.lineJoin = 'round'; ctx.beginPath();
        for (let k = 0; k <= now; k++) { const px = X(k), py = Y(r[field][from + k] + shift); if (k) ctx.lineTo(px, py); else ctx.moveTo(px, py); }
        ctx.stroke();
      }
      ctx.restore();
      if (p.guide) {
        const words = `${p.guide[1]}, ${num(p.guide[0], 2)}`;
        ctx.font = '400 9.5px system-ui, sans-serif';
        const width = ctx.measureText(words).width + 8, start = Math.max(0, Math.floor((w - 4 - width) / w * (DAYS - 1))), marks = [];
        for (const [field] of p.lines) for (let k = start; k <= now; k++) marks.push(Y(r[field][from + k] + shift));
        guideLabel(ctx, words, width, x + w - 4, gy, y0 + HEAD, y1 - 2, marks);
      }
      ctx.strokeStyle = PROBE; ctx.lineWidth = 1; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(X(now), y0); ctx.lineTo(X(now), y1); ctx.stroke(); ctx.setLineDash([]);
      for (const [field, color] of p.lines) { ctx.fillStyle = color; ctx.strokeStyle = HALO; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(X(now), Y(r[field][d] + shift), 3.2, 0, 2 * Math.PI); ctx.fill(); ctx.stroke(); }
      text(ctx, p.name, x + 5, y0 + 9, { color: MUTED, size: 10, halo: HALO });
    });
    const base = top + 3 * panel + 2 * gap;
    ctx.font = '400 10px system-ui, sans-serif';
    const fits = (labels) => labels.every((label, n) => {
      if (!n) return true;
      const before = ctx.measureText(labels[n - 1]).width * (n === 1 ? 1 : 0.5), mine = ctx.measureText(label).width * (n === 4 ? 1 : 0.5);
      return before + mine + 10 < X(Math.min(MARKS[n], DAYS - 1)) - X(MARKS[n - 1]);
    });
    const names = fits(NAMES) ? NAMES : SHORT;
    MARKS.forEach((m, n) => {
      const mx = X(Math.min(m, DAYS - 1));
      ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(mx, base); ctx.lineTo(mx, base + 4); ctx.stroke();
      text(ctx, names[n], mx, base + 14, { align: n === 0 ? 'left' : n === 4 ? 'right' : 'center', color: MUTED, size: 10 });
    });
  }

  function guideLabel(ctx, words, width, right, gy, top, bottom, marks) {
    const above = gy - 8, below = gy + 8;
    const clear = (cy) => cy - 6 >= top && cy + 6 <= bottom && marks.every((my) => Math.abs(my - cy) > 11);
    const cy = clear(above) ? above : clear(below) ? below : gy - 8 >= top ? above : below;
    ctx.fillStyle = CHIP; ctx.beginPath(); ctx.roundRect(right - width, cy - 6.5, width + 2, 13, 3); ctx.fill();
    text(ctx, words, right - 2, cy, { color: MUTED, size: 9.5, align: 'right' });
  }

  function report(d, y) {
    const r = run.record, h = r.thickness[d], every = playing ? 4 : 1, s = summary(y), latent = run.sea.latent;
    const c = run.conducted[d], growth = 100 * 86400 * c / latent, rho = r.density[d] - DEEP_CLASS, light = r.sunlight[d];
    const against = (v) => `${num(Math.abs(v), 2)} kg/m³ ${v > 0 ? 'denser' : 'lighter'} than the deepest class`;
    dayOut.set([
      [`year ${y + 1},`, seasonPhrase(r.day[d])],
      ['ice', h > 0 ? `${h.toFixed(2)} m thick over ${percent(r.cover[d])} of the sea` : 'none'],
      ...(state.snow && r.snow[d] > 0 ? [['snow', `${snowDepth(r.snow[d]).toFixed(2)} m`]] : []),
      ...(h > 0 ? [['conducted', `${num(Math.abs(c), 1)} W/m² ${c >= 0 ? 'up' : 'down'} through it, ${c >= 0 ? 'freezing' : 'melting'} ${Math.abs(growth).toFixed(2)} cm a day at its base`]] : []),
      ['sunlight', light > 0 ? `${Math.round(light)} W/m², ${percent(r.albedo[d])} reflected` : 'none, polar night'],
      ['water', `${celsius(r.water[d], 2)}, ${r.salinity[d].toFixed(2)} psu, ${against(rho)}`],
      ...(h > 0 ? [] : [['heat to lose before it freezes', heat(heatAboveFreezing(r.water[d], MIXED[state.depth]))]]),
    ], every);
    const when = (i) => seasonPhrase(r.day[i]);
    const season = s.none ? 'none' : s.allYear ? `all year, ${s.thinnest.toFixed(2)} to ${s.thickest.toFixed(2)} m for each square meter of sea` : `${s.forms === null ? 'from the year before' : `forms ${when(s.forms)}`}, ${s.melts === null ? 'lasts into the next year' : `gone ${when(s.melts)}`}, up to ${s.thickest.toFixed(2)} m for each square meter of sea`;
    yearOut.set([
      [`year ${y + 1}: ice`, season],
      ['brine of the winter', `${s.brine.toFixed(2)} psu`],
      ['densest water', against(s.againstClass)],
      ...(s.summerHeat >= 1e6 ? [['most heat above freezing', `${heat(s.summerHeat)}, enough to melt ${(s.summerHeat / latent).toFixed(1)} m of ice`]] : []),
    ], every);
    if (yearControl.value !== y + 1) yearControl.value = y + 1;
  }

  function pick({ x }) {
    if (!geometry) return false;
    const { x: cx, w } = geometry.charts, k = clamp(Math.round((x - cx) / w * (DAYS - 1)), 0, DAYS - 1);
    position = yearOf(position) * DAYS + k;
    halt();
    fig.render();
    return true;
  }
  const inCharts = (p) => {
    if (!geometry) return false;
    const { x, w, top, panel, gap } = geometry.charts;
    return p.x >= x - 6 && p.x <= x + w + 6 && p.y >= top - 4 && p.y <= top + 3 * panel + 2 * gap + 4;
  };
  const near = (p) => {
    if (!inCharts(p)) return false;
    const { x, w } = geometry.charts;
    return Math.abs(p.x - (x + dayOf(position) / (DAYS - 1) * w)) < 16;
  };
  const at = (e) => { const b = fig.canvas.getBoundingClientRect(); return { x: e.clientX - b.left, y: e.clientY - b.top }; };
  let dragging = false;
  fig.pointer({ hit: near, down: (p) => { dragging = true; return pick(p); }, move: (p) => { if (dragging) pick(p); }, up: () => { dragging = false; } });
  fig.canvas.addEventListener('click', (e) => { const p = at(e); if (inCharts(p) && !near(p)) pick(p); });
  fig.canvas.addEventListener('pointermove', (e) => { if (e.pointerType !== 'mouse' || dragging) return; const p = at(e); fig.canvas.style.cursor = near(p) ? 'grab' : inCharts(p) ? 'crosshair' : ''; });

  compute();
  fig.play(true);
}
