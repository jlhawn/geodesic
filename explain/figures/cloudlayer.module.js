import { Figure, slider, choice, readout, legend, caption, text, arrow, termColor, anomalyColor, clamp, INK, MUTED, GRID, LINE, ACCENT, VAPOR, CLOUD } from '../runtime.module.js';
import { loadColumn, makeColumn, K } from '../radiativeColumn.module.js';
import { createSweep, zeroCrossings, cloudFromHumidity, LAYERS, PATHS, SUNS } from '../cloudEffects.module.js';
import { liquidShare } from '../../js/physics/radiation.module.js';

const WIDE = 520, Z_TOP = 15000, BUDGET = 6, RATE = 10, HALO = 'rgba(20, 20, 22, 0.85)', GROUND = '#5a4a36';
const ZERO = 'rgba(255, 255, 255, 0.85)';
const LOWEST = LAYERS.indexOf(26), HUNDRED = PATHS.findIndex((w) => Math.abs(w - 0.1) < 1e-12);
const SCALES = [50, 100, 150, 200, 300, 400, 600, 800, 1000];
const R_LO = 50, R_HI = 150;
const SUN_OPTIONS = [['average', 'average'], ['overhead', 'overhead'], ['15° up', 'low'], ['night', 'night']];

const signed = (v, d = 1) => (Math.abs(v) < 0.5 * 10 ** -d ? (0).toFixed(d) : `${v < 0 ? '−' : '+'}${Math.abs(v).toFixed(d)}`);
const hPa = (p) => `${(p / 100).toFixed(0)} hPa`;
const distance = (z) => (z < 995 ? `${Math.round(z / 10) * 10} m` : `${(z / 1000).toFixed(1)} km`);
const grams = (path) => { const g = path * 1000; return `${g < 9.95 ? g.toFixed(1) : g.toFixed(0)} g/m²`; };
const perKilogram = (q) => { const g = q * 1000; return g >= 0.01 ? g.toFixed(3) : g.toPrecision(2); };
const phase = (T) => { const f = liquidShare(T); return f >= 1 ? 'liquid' : f <= 0 ? 'ice' : `${Math.round(100 * f)}% liquid, the rest ice`; };

export function mountCloudlayer(root) {
  const controls = root.querySelector('.controls'), A = termColor('a'), B = termColor('b'), E = termColor('e');
  const state = { view: 'effect', row: LOWEST, col: HUNDRED, cover: 1, humidity: 1, sun: 'average', air: 'moist' };
  let profile = null, column = null, sweep = null, made = null, goal = null, shown = null, scale = 200, touched = false, geometry = null, failed = false, lastSun = null;

  const fig = new Figure(root, { height: (w) => (w < WIDE ? 450 : 440), step, draw });
  const words = caption(root, '');
  const key = legend(root, []);
  choice(controls, { label: 'Show', options: [['what a cloud does', 'effect'], ['how much cloud', 'amount']], value: state.view, onChange: (v) => { state.view = v; arrange(); update(); } });
  const height = slider(controls, { label: 'Cloud’s height', min: 0, max: LAYERS.length - 1, step: 1, value: state.row, format: (r) => where(r), onInput: (r) => { state.row = r; touched = true; update(); } });
  const water = slider(controls, { label: 'Water in the cloud', min: 0, max: PATHS.length - 1, step: 1, value: state.col, format: (c) => grams(PATHS[c]), onInput: (c) => { state.col = c; touched = true; update(); } });
  const waterRow = controls.lastElementChild;
  slider(controls, { label: 'Cover', min: 0, max: 100, step: 1, value: 100, format: (v) => `${v}% of the sky`, onInput: (v) => { state.cover = v / 100; update(); } });
  const coverRow = controls.lastElementChild;
  const humid = slider(controls, { label: 'Humidity of the layer', min: 60, max: 130, step: 1, value: 100, format: (v) => `${v}% of saturation`, onInput: (v) => { state.humidity = v / 100; update(); } });
  const humidRow = controls.lastElementChild;
  choice(controls, { label: 'Sun', options: SUN_OPTIONS, value: state.sun, onChange: (v) => { state.sun = v; restart(); } });
  choice(controls, { label: 'Air', options: [['moist', 'moist'], ['dry', 'dry']], value: state.air, onChange: (v) => { state.air = v; restart(); } });
  const layerOut = readout(controls), fluxOut = readout(controls), effectOut = readout(controls), heatOut = readout(controls);
  arrange();

  function where(r) {
    const k = LAYERS[r];
    if (!column) return `layer ${K - k} of ${K}`;
    return `${distance(column.height[k])}, ${hPa(column.pressure[k])}`;
  }

  function arrange() {
    const effect = state.view === 'effect';
    waterRow.style.display = coverRow.style.display = effect ? '' : 'none';
    humidRow.style.display = effect ? 'none' : '';
    const arrows = [['arrow', 'sunlight', A], ['arrow', 'infrared', B], ['arrow', 'both together', E]];
    if (effect) {
      words.textContent = 'One layer of cloud placed in Part 6’s midlatitude summer, with its water vapor, over ground that reflects 7 percent of the sunlight, under the planet’s average sunlight unless the Sun control changes it. The model’s own radiation runs on the column twice, with the cloud and with the cloud taken out, and the difference at the top of the air is the cloud’s effect on the planet. The map colors every cloud the figure can place, in each of the model’s layers from the ground to 141 hPa and with 5 to 400 g of water on each square meter, by that effect, sunlight and infrared together. The ring is the chosen cloud, and the arrows below split its effect in two. Drag the ring or tap the map to move it.';
      key.set([['gradient', 'net effect, cooling to warming', 'cool', 'warm'], ['dash', 'no net effect', ZERO], ...arrows]);
    } else {
      words.textContent = 'How the model decides how much cloud a layer holds. It takes the water in each layer of a grid cell, vapor and cloud together, to vary across the cell evenly between a driest and a moistest part, and the part above saturation condenses. The ramp is the cell’s air sorted from its driest part to its moistest, and the white wedge is cloud. Its width is the cloud’s cover. The arrows below are that cloud’s effect on the column, as in the other view. Drag the dot or tap the chart to change the layer’s humidity.';
      key.set([['bar', 'water vapor', VAPOR], ['bar', 'cloud: the water above saturation', CLOUD], ['dash', 'saturation', ACCENT], ...arrows]);
    }
  }

  function restart() {
    if (!sweep) return;
    sweep.start({ humidity: state.air === 'moist' ? profile.humidity : null, sun: SUNS[state.sun] });
    scale = null;
    update();
  }

  function update() {
    if (!column) return;
    const k = LAYERS[state.row], effect = state.view === 'effect';
    let humidity = state.air === 'moist' ? profile.humidity : null, temperature = profile.temperature, cloud = null;
    made = null;
    if (effect) cloud = state.cover > 0 ? { pressure: column.pressure[k], path: PATHS[state.col], cover: state.cover } : null;
    else ({ made, humidity, temperature, cloud } = cloudFromHumidity(profile, k, state.humidity, humidity));
    const sunChanged = lastSun !== state.sun;
    lastSun = state.sun;
    column.set({ humidity, temperature, cloud: null, ...(sunChanged ? SUNS[state.sun] : {}) });
    const clear = { dlr: column.budget.dlr, olr: column.budget.olr, reflected: column.budget.reflected, infrared: column.heating.infrared[k], sunlight: column.heating.sunlight[k] };
    column.set({ cloud });
    const b = column.budget;
    goal = { sw: b.shortwaveCloudEffect, lw: b.longwaveCloudEffect, net: b.shortwaveCloudEffect + b.longwaveCloudEffect, height: column.height[k], path: effect ? PATHS[state.col] : made.path };
    if (!shown) shown = { ...goal };
    report(k, b, clear);
    fig.play(true);
    fig.render();
  }

  function report(k, b, clear) {
    const T = column.T[k], effect = state.view === 'effect';
    const layer = [['the cloud’s layer', `${hPa(column.pressure[k])}, ${distance(column.height[k])} up, at ${T.toFixed(1)} K: ${phase(T)}`]];
    if (!effect) layer.push(['critical humidity', `${(100 * made.critical).toFixed(1)}%`], ['cover', `${(100 * made.cover).toFixed(made.cover > 0 && made.cover < 0.1 ? 1 : 0)}%`], ['cloud water', made.condensate > 0 ? `${perKilogram(made.condensate)} g/kg, ${grams(made.path)} where it covers` : 'none']);
    layerOut.set(layer);
    fluxOut.set([
      ['sunlight reflected to space', `${b.reflected.toFixed(1)} W/m², ${clear.reflected.toFixed(1)} without the cloud`],
      ['infrared leaving the top', `${b.olr.toFixed(1)} W/m², ${clear.olr.toFixed(1)} without`],
      ['infrared reaching the ground', `${b.dlr.toFixed(1)} W/m², ${clear.dlr.toFixed(1)} without`],
    ]);
    effectOut.set([['the cloud’s effect: sunlight', `${signed(goal.sw)} W/m²`], ['infrared', `${signed(goal.lw)} W/m²`], ['together', `${signed(goal.net)} W/m²`]]);
    heatOut.set([['the cloud’s layer, by infrared', `${signed(column.heating.infrared[k])} K/day, ${signed(clear.infrared)} without the cloud`], ['by sunlight', `${signed(column.heating.sunlight[k])} K/day, ${signed(clear.sunlight)} without`]]);
  }

  function step(dt) {
    if (!column) return false;
    let busy = false;
    if (!sweep.done) { sweep.advance(BUDGET); busy = true; }
    if (scale === null && sweep.count >= 2) scale = SCALES.find((s) => s >= sweep.largest() * 1.02) ?? SCALES.at(-1);
    const f = 1 - Math.exp(-RATE * dt);
    for (const name of ['sw', 'lw', 'net', 'height']) {
      const d = goal[name] - shown[name];
      if (Math.abs(d) > (name === 'height' ? 1 : 0.02)) { shown[name] += d * f; busy = true; } else shown[name] = goal[name];
    }
    const d = Math.log(goal.path || PATHS[0]) - Math.log(shown.path || PATHS[0]);
    if (Math.abs(d) > 1e-3) { shown.path = Math.exp(Math.log(shown.path || PATHS[0]) + d * f); busy = true; } else shown.path = goal.path;
    if (!busy) fig.play(false);
    return true;
  }

  function layout(w, h) {
    const narrow = w < WIDE, size = narrow ? 10 : 11, left = narrow ? 32 : 40, right = w - (narrow ? 40 : 48), top = 30;
    const effects = 112, bottom = h - effects - 48;
    return { narrow, size, left, right, top, bottom, effectsTop: bottom + 50 };
  }

  function rowEdges() {
    const heights = LAYERS.map((k) => column.height[k]);
    return heights.map((z, r) => [r === 0 ? 0 : (z + heights[r - 1]) / 2, r === heights.length - 1 ? Math.min(Z_TOP, z + (z - heights[r - 1]) / 2) : (z + heights[r + 1]) / 2]);
  }

  function pathEdges() {
    const logs = PATHS.map(Math.log);
    return logs.map((v, c) => [c === 0 ? v - (logs[1] - v) / 2 : (v + logs[c - 1]) / 2, c === logs.length - 1 ? v + (v - logs[c - 1]) / 2 : (v + logs[c + 1]) / 2]);
  }

  function draw(ctx, w, h) {
    const g = layout(w, h), { size } = g;
    if (!column) { text(ctx, failed ? 'The midlatitude summer could not be loaded.' : 'Loading the midlatitude summer…', w / 2, h / 2, { align: 'center', color: MUTED, size: 12 }); return; }
    if (state.view === 'effect') map(ctx, w, g); else ramp(ctx, w, g);
    effects(ctx, w, g);
  }

  function map(ctx, w, g) {
    const { narrow, size, left, right, top, bottom } = g, rows = rowEdges(), cols = pathEdges();
    const lo = cols[0][0], hi = cols.at(-1)[1], x = (v) => left + (v - lo) / (hi - lo) * (right - left), y = (z) => bottom - z / Z_TOP * (bottom - top);
    geometry = { kind: 'map', left, right, top, bottom, x, y, rows, cols };
    const net = Float64Array.from(sweep.shortwave, (s, n) => s + sweep.longwave[n]), s = scale ?? 200;
    text(ctx, narrow ? 'every cloud’s net effect' : 'the net effect of a cloud at every height and amount of water', left, top - 16, { color: INK, size: size + 0.5, weight: 500 });
    if (!sweep.done) text(ctx, `computing ${Math.round(100 * sweep.count / net.length)}%`, right, top - 16, { align: 'right', color: ACCENT, size });
    ctx.fillStyle = 'rgba(255, 255, 255, 0.03)'; ctx.fillRect(left, top, right - left, bottom - top);
    for (let r = 0; r < sweep.rows; r++) for (let c = 0; c < sweep.cols; c++) {
      const v = net[r * sweep.cols + c] * state.cover;
      if (!Number.isFinite(v)) continue;
      const t = Math.sign(v) * Math.min(1, Math.abs(v) / s) ** 0.7;
      const x0 = Math.round(x(cols[c][0])), x1 = Math.round(x(cols[c][1])), y0 = Math.round(y(rows[r][1])), y1 = Math.round(y(rows[r][0]));
      ctx.fillStyle = anomalyColor(t, 0.9); ctx.fillRect(x0, y0, x1 - x0, y1 - y0);
    }
    ctx.lineWidth = 1;
    const zStep = narrow ? 5000 : 3000;
    for (let z = 0; z <= Z_TOP; z += zStep) {
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(left, y(z)); ctx.lineTo(right, y(z)); ctx.stroke();
      text(ctx, `${z / 1000}`, left - 6, y(z), { align: 'right', color: MUTED, size });
      text(ctx, `${airAt(z).toFixed(0)} K`, right + 5, y(z), { color: MUTED, size });
    }
    text(ctx, 'km', left - 6, top - 13, { align: 'right', color: MUTED, size, baseline: 'bottom' });
    text(ctx, 'air', right + 5, top - 13, { color: MUTED, size, baseline: 'bottom' });
    for (const gram of narrow ? [5, 20, 100, 400] : [5, 10, 20, 50, 100, 200, 400]) {
      const px = x(Math.log(gram / 1000));
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(px, top); ctx.lineTo(px, bottom + 4); ctx.stroke();
      text(ctx, `${gram}`, px, bottom + 14, { align: 'center', color: MUTED, size });
    }
    text(ctx, 'water in the cloud, g/m²', (left + right) / 2, bottom + 29, { align: 'center', color: MUTED, size });
    ctx.fillStyle = GROUND; ctx.fillRect(left, bottom, right - left, 3);

    if (sweep.done) {
      const lines = zeroCrossings(net, sweep.rows, sweep.cols), mid = (r) => (rows[r][0] + rows[r][1]) / 2;
      ctx.save(); ctx.strokeStyle = ZERO; ctx.lineWidth = 1.5; ctx.setLineDash([4, 3]);
      let label = null;
      for (let r = 0; r < sweep.rows; r++) {
        const here = lines[r], next = lines[r + 1] ?? [];
        for (let i = 0; i < here.length; i++) {
          const px = x(Math.log(here[i])), py = y(mid(r));
          if (i < next.length) { ctx.beginPath(); ctx.moveTo(px, py); ctx.lineTo(x(Math.log(next[i])), y(mid(r + 1))); ctx.stroke(); }
          if (i === 0) label = [px, py];
        }
      }
      ctx.restore();
      if (label) {
        ctx.font = `400 ${size}px system-ui, sans-serif`;
        const width = ctx.measureText('no net effect').width, [lx, ly] = label, dot = [x(Math.log(shown.path)), y(shown.height)];
        const spots = [[lx + 8, ly + 9], [lx - 8 - width, ly - 10], [lx - 8 - width, ly + 12], [lx + 8, ly - 10]].map(([a, b]) => [a, clamp(b, top + 8, bottom - 8)]);
        const clear = ([a, b]) => a >= left + 2 && a + width <= right - 2 && (dot[0] < a - 10 || dot[0] > a + width + 10 || Math.abs(dot[1] - b) > 16);
        const [a, b] = spots.find(clear) ?? spots[0];
        text(ctx, 'no net effect', a, b, { color: '#fff', size, halo: HALO });
      }
      let least = 0, most = 0;
      for (let n = 0; n < net.length; n++) { if (net[n] < net[least]) least = n; if (net[n] > net[most]) most = n; }
      if (net[least] * state.cover < -0.05 * s) text(ctx, 'cools', x(Math.log(0.15)), y(2600), { align: 'center', color: '#fff', size: size + 1, weight: 500, halo: HALO });
      if (net[most] * state.cover > 0.05 * s) {
        const r = Math.floor(most / sweep.cols), c = most % sweep.cols;
        text(ctx, 'warms', clamp(x(Math.log(PATHS[c])), left + 24, right - 24), clamp(y(mid(r)) + 18, top + 10, bottom - 10), { align: 'center', color: '#fff', size: size + 1, weight: 500, halo: HALO });
      }
    }

    const px = x(Math.log(shown.path)), py = y(shown.height);
    ctx.fillStyle = HALO; ctx.strokeStyle = '#fff'; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.arc(px, py, 6.5, 0, Math.PI * 2); ctx.stroke();
    ctx.fillStyle = CLOUD; ctx.beginPath(); ctx.arc(px, py, 2.5, 0, Math.PI * 2); ctx.fill();
    geometry.dot = [px, py];
    if (!touched) { const before = px > (left + right) / 2; text(ctx, 'drag', before ? px - 12 : px + 12, py - 12, { align: before ? 'right' : 'left', color: ACCENT, size: 11, halo: HALO }); }
  }

  function airAt(z) {
    if (z <= 0) return profile.groundTemperature;
    for (let k = K - 1; k > 0; k--) if (column.height[k - 1] >= z) {
      const z0 = column.height[k], z1 = column.height[k - 1];
      if (z <= z0) return profile.groundTemperature + (column.T[k] - profile.groundTemperature) * z / z0;
      return column.T[k] + (column.T[k - 1] - column.T[k]) * (z - z0) / (z1 - z0);
    }
    return column.T[0];
  }

  function ramp(ctx, w, g) {
    const { narrow, size, left, right, top, bottom } = g, k = LAYERS[state.row], c = made, mean = 100 * state.humidity, half = 100 * c.spread;
    const x = (f) => left + f * (right - left), y = (v) => bottom - (v - R_LO) / (R_HI - R_LO) * (bottom - top);
    geometry = { kind: 'ramp', left, right, top, bottom, y, x };
    text(ctx, narrow ? 'water across one grid cell' : `the ${hPa(column.pressure[k])} layer: its water across one grid cell`, left, top - 16, { color: INK, size: size + 0.5, weight: 500 });
    ctx.fillStyle = 'rgba(255, 255, 255, 0.03)'; ctx.fillRect(left, top, right - left, bottom - top);
    ctx.lineWidth = 1;
    for (let v = R_LO; v <= R_HI; v += 10) {
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(left, y(v)); ctx.lineTo(right, y(v)); ctx.stroke();
      if (!narrow || v % 20 === 0) text(ctx, `${v}`, left - 6, y(v), { align: 'right', color: MUTED, size });
    }
    text(ctx, '%', left - 6, top - 13, { align: 'right', color: MUTED, size, baseline: 'bottom' });
    for (const f of [0, 0.25, 0.5, 0.75, 1]) text(ctx, `${100 * f}%`, x(f), bottom + 14, { align: f === 0 ? 'left' : f === 1 ? 'right' : 'center', color: MUTED, size });
    text(ctx, narrow ? 'the cell, driest part to moistest' : 'the grid cell’s air, from its driest part to its moistest', (left + right) / 2, bottom + 29, { align: 'center', color: MUTED, size });

    const at = (f) => mean - half + 2 * half * f, sat = (100 - (mean - half)) / (2 * half);
    const cut = half > 0 ? clamp(sat, 0, 1) : (mean >= 100 ? 0 : 1);
    ctx.save(); ctx.beginPath(); ctx.rect(left, top, right - left, bottom - top); ctx.clip();
    ctx.fillStyle = VAPOR; ctx.globalAlpha = 0.5; ctx.beginPath();
    ctx.moveTo(x(0), y(R_LO)); ctx.lineTo(x(0), y(Math.min(100, at(0)))); ctx.lineTo(x(cut), y(Math.min(100, at(cut)))); ctx.lineTo(x(1), y(Math.min(100, at(1)))); ctx.lineTo(x(1), y(R_LO)); ctx.closePath(); ctx.fill();
    ctx.globalAlpha = 1;
    if (cut < 1) {
      ctx.fillStyle = CLOUD; ctx.beginPath();
      ctx.moveTo(x(cut), y(Math.max(100, at(cut)))); ctx.lineTo(x(1), y(at(1))); ctx.lineTo(x(1), y(100)); ctx.lineTo(x(cut), y(100)); ctx.closePath(); ctx.fill();
    }
    ctx.strokeStyle = INK; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(x(0), y(at(0))); ctx.lineTo(x(1), y(at(1))); ctx.stroke();
    ctx.restore();
    ctx.strokeStyle = ACCENT; ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]); ctx.beginPath(); ctx.moveTo(left, y(100)); ctx.lineTo(right, y(100)); ctx.stroke(); ctx.setLineDash([]);
    text(ctx, 'saturation', left + 6, y(100) - 9, { color: ACCENT, size, halo: HALO });
    const critical = 100 * c.critical;
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.55)'; ctx.lineWidth = 1; ctx.setLineDash([1.5, 3]); ctx.beginPath(); ctx.moveTo(left, y(critical)); ctx.lineTo(x(0.5), y(critical)); ctx.stroke(); ctx.setLineDash([]);
    text(ctx, `${narrow ? 'cloud begins at' : 'cloud begins at a mean of'} ${critical.toFixed(1)}%`, left + 6, y(critical) + 9, { color: '#fff', size, halo: HALO });
    if (cut < 1) {
      const yb = clamp(y(at(1)) - 8, top + 24, y(100) - 6);
      ctx.strokeStyle = '#fff'; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(x(cut), yb + 4); ctx.lineTo(x(cut), yb); ctx.lineTo(x(1) - 1, yb); ctx.lineTo(x(1) - 1, yb + 4); ctx.stroke();
      const label = `cover ${(100 * c.cover).toFixed(c.cover < 0.1 ? 1 : 0)}%`;
      ctx.font = `400 ${size}px system-ui, sans-serif`;
      const width = ctx.measureText(label).width;
      text(ctx, label, clamp((x(cut) + x(1)) / 2, left + width / 2 + 4, right - width / 2 - 4), yb - 10, { align: 'center', color: '#fff', size, halo: HALO });
    }
    const my = y(mean);
    ctx.fillStyle = '#fff'; ctx.strokeStyle = HALO; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(x(0.5), my, 5, 0, Math.PI * 2); ctx.stroke(); ctx.fill();
    text(ctx, `mean ${mean.toFixed(0)}%`, x(0.5) + 10, my + 13, { color: '#fff', size, halo: HALO });
    geometry.dot = [x(0.5), my];
    if (!touched) text(ctx, 'drag', x(0.5) - 10, my - 12, { align: 'right', color: ACCENT, size: 11, halo: HALO });
  }

  function effects(ctx, w, g) {
    const { narrow, size, effectsTop } = g, labelEnd = narrow ? 62 : 74, x0 = labelEnd + 8, x1 = w - 14, cx = (x0 + x1) / 2, s = scale ?? 200;
    const x = (v) => cx + clamp(v / s, -1, 1) * (x1 - cx), rowY = (n) => effectsTop + 30 + n * 22;
    text(ctx, narrow ? 'the cloud’s effect, W/m²' : 'the cloud’s effect on the planet, W/m² at the top of the air', 8, effectsTop - 2, { color: INK, size: size + 0.5, weight: 500 });
    text(ctx, 'cools', cx - 8, effectsTop + 14, { align: 'right', color: MUTED, size });
    text(ctx, 'warms', cx + 8, effectsTop + 14, { color: MUTED, size });
    ctx.lineWidth = 1;
    for (const f of [-1, -0.5, 0, 0.5, 1]) {
      const px = cx + f * (x1 - cx);
      ctx.strokeStyle = f ? GRID : LINE; ctx.beginPath(); ctx.moveTo(px, effectsTop + 21); ctx.lineTo(px, rowY(2) + 9); ctx.stroke();
      text(ctx, f ? `${f < 0 ? '−' : ''}${Math.abs(f * s)}` : '0', px, rowY(2) + 20, { align: f === -1 ? 'left' : f === 1 ? 'right' : 'center', color: MUTED, size });
    }
    [['sunlight', shown.sw, A], ['infrared', shown.lw, B], ['together', shown.net, E]].forEach(([name, v, color], n) => {
      const py = rowY(n), tip = x(v);
      text(ctx, name, labelEnd, py, { align: 'right', color, size: size + 0.5 });
      if (Math.abs(tip - cx) >= 3) arrow(ctx, cx, py, tip, py, { color, width: n === 2 ? 3 : 2.5, head: 7 });
      else { ctx.fillStyle = color; ctx.beginPath(); ctx.arc(cx, py, 2.5, 0, Math.PI * 2); ctx.fill(); }
      const label = signed(v), outside = v < 0 ? tip - 8 : tip + 8;
      ctx.font = `500 ${size}px system-ui, sans-serif`;
      const width = ctx.measureText(label).width, fits = v < 0 ? outside - width >= x0 - 4 : outside + width <= x1 + 4;
      if (Math.abs(tip - cx) < 3) text(ctx, label, cx + 8, py, { color: '#fff', size, weight: 500, halo: HALO });
      else if (fits) text(ctx, label, outside, py, { align: v < 0 ? 'right' : 'left', color: '#fff', size, weight: 500, halo: HALO });
      else text(ctx, label, v < 0 ? tip + 4 : tip - 4, py - 9, { align: v < 0 ? 'left' : 'right', color: '#fff', size, weight: 500, halo: HALO });
    });
  }

  function pickMap({ x: px, y: py }) {
    const gm = geometry;
    if (!gm || gm.kind !== 'map' || px < gm.left - 8 || px > gm.right + 8 || py < gm.top - 8 || py > gm.bottom + 8) return false;
    const lo = gm.cols[0][0], hi = gm.cols.at(-1)[1], v = lo + (clamp(px, gm.left, gm.right) - gm.left) / (gm.right - gm.left) * (hi - lo);
    const z = (gm.bottom - clamp(py, gm.top, gm.bottom)) / (gm.bottom - gm.top) * Z_TOP;
    let c = 0, r = 0;
    PATHS.forEach((p, n) => { if (Math.abs(Math.log(p) - v) < Math.abs(Math.log(PATHS[c]) - v)) c = n; });
    gm.rows.forEach(([a, b], n) => { if (z >= a && z <= b) r = n; });
    if (z > gm.rows.at(-1)[1]) r = gm.rows.length - 1;
    touched = true;
    if (c !== state.col || r !== state.row) { state.col = c; state.row = r; water.value = c; height.value = r; update(); } else fig.render();
    return true;
  }

  function pickRamp({ y: py }) {
    const gm = geometry;
    if (!gm || gm.kind !== 'ramp' || py < gm.top - 8 || py > gm.bottom + 8) return false;
    const v = Math.round(clamp(R_HI - (py - gm.top) / (gm.bottom - gm.top) * (R_HI - R_LO), 60, 130));
    touched = true;
    if (v !== Math.round(100 * state.humidity)) { state.humidity = v / 100; humid.value = v; update(); } else fig.render();
    return true;
  }

  const pick = (at) => (state.view === 'effect' ? pickMap(at) : pickRamp(at));
  const near = ({ x: px, y: py }) => !!geometry?.dot && Math.hypot(px - geometry.dot[0], py - geometry.dot[1]) < 22;
  fig.pointer({ hit: near, down: pick, move: pick });
  fig.canvas.addEventListener('click', (e) => { const r = fig.canvas.getBoundingClientRect(); pick({ x: e.clientX - r.left, y: e.clientY - r.top }); });

  loadColumn().then((loaded) => {
    profile = loaded;
    column = makeColumn({ profile, humidity: profile.humidity });
    sweep = createSweep(profile);
    height.value = state.row;
    restart();
  }).catch((error) => { failed = true; console.info('cloud figure', error.message); fig.render(); });
}
