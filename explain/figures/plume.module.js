import { Figure, slider, choice, readout, legend, caption, text, arrow, termColor, clamp, INK, MUTED, GRID, LINE, ACCENT, VAPOR, CLOUD } from '../runtime.module.js';
import { plumeColumn, SPACINGS } from '../plumeColumn.module.js';
import { readPlume } from '../plumePath.module.js';
import { dewPoint } from '../moistParcel.module.js';

const WIDE = 560, TOP = 100e2, T_RANGE = [180, 305], RAIN = [0.1, 1000], RATE = 8;
const HUMIDITY = [0.4, 1.1], WARMING = [0, 1];
const HALO = 'rgba(20, 20, 22, 0.85)', COLDER = 'rgba(255, 255, 255, 0.1)';
const PRESSURES = [1000, 850, 700, 500, 400, 300, 200, 150, 100];
const STEPS = [1, 2, 2.5, 5, 10], ROOMS = [1, 1.5, 2, 3, 4, 5, 6, 8, 10];
const decade = (v) => 10 ** Math.floor(Math.log10(v));
const nice = (v) => (v > 0 ? decade(v) * STEPS.find((s) => s * decade(v) >= v * 0.999) : 1);
const room = (v) => (v > 0 ? decade(v) * ROOMS.find((s) => s * decade(v) >= v * 0.999) : 0);
const hPa = (p) => `${Math.round(p / 100)} hPa`;
const km = (z) => `${(z / 1000).toFixed(1)} km`;
const signed = (v, d = 1) => (Math.abs(v) < 0.5 * 10 ** -d ? (0).toFixed(d) : `${v < 0 ? '−' : '+'}${Math.abs(v).toFixed(d)}`);
const minus = (v) => `${+v.toFixed(2)}`.replace('-', '−');
const rainText = (v) => (v >= 99.95 ? v.toFixed(0) : v >= 9.995 ? v.toFixed(1) : v.toFixed(2));

function withAlpha(color, alpha) {
  const hex = /^#([0-9a-f]{6})$/i.exec(color);
  return hex ? `rgba(${[0, 2, 4].map((n) => parseInt(hex[1].slice(n, n + 2), 16)).join(', ')}, ${alpha})` : color;
}

function ticks(lo, hi, count) {
  const step = nice((hi - lo) / count), out = [];
  for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + 1e-9 * step; v += step) out.push(Math.abs(v) < 1e-9 * step ? 0 : v);
  return out;
}

export function mountPlume(root) {
  const controls = root.querySelector('.controls'), D = termColor('d'), FILL = withAlpha(D, 0.32), EDGE = withAlpha(CLOUD, 0.55);
  const settings = { humidity: 1, warming: 0, spacing: SPACINGS[0], ice: true };
  let result = null, read = null, probe = null, touched = false, geometry = null, range = null, goal = null, failed = false;

  const fig = new Figure(root, { height: (w) => (w < WIDE ? 610 : 480), step, draw });
  caption(root, 'One step of the model’s own moist physics on one column, starting from C. L. Jordan’s mean sounding for the West Indies in the hurricane season. All four panels share the pressure axis. The second shows how much warmer the plume is than the air around it, counting its vapor and cloud water; its shaded part makes up the CAPE. The last counts the drying as the latent heat of the vapor removed. Tap or drag to read a layer.');
  legend(root, [['line', 'air temperature', INK], ['dash', 'dew point', VAPOR], ['line', 'the plume, its mass flux and its heating', D], ['bar', 'plume warmer than the air', FILL], ['bar', 'plume colder than the air', COLDER], ['dash', 'downdraft', D], ['line', 'drying', VAPOR], ['dots', 'cloud base and top', EDGE]]);
  controls.classList.add('two');
  slider(controls, { label: 'Humidity above 850 hPa', min: HUMIDITY[0], max: HUMIDITY[1], step: 0.01, value: settings.humidity, format: (v) => `×${v.toFixed(2)}`, onInput: (v) => change({ humidity: v }), span: true });
  slider(controls, { label: 'Warm the lowest kilometer', min: WARMING[0], max: WARMING[1], step: 0.05, value: settings.warming, format: (v) => `+${v.toFixed(2)} K`, onInput: (v) => change({ warming: v }), span: true });
  choice(controls, { label: 'Grid cells', options: [['as at N = 128', '0'], ['as at N = 64', '1']], value: '0', onChange: (v) => change({ spacing: SPACINGS[Number(v)] }) });
  choice(controls, { label: 'Ice in the plume', options: [['on', 'on'], ['off', 'off']], value: 'on', onChange: (v) => change({ ice: v === 'on' }) });
  const summary = readout(controls), layer = readout(controls);
  for (const box of controls.querySelectorAll('.readout')) box.style.gridColumn = '1 / -1';

  function change(changes) {
    Object.assign(settings, changes);
    compute();
    goal = target();
    fig.play(true);
    fig.render();
  }

  function compute() {
    try {
      result = plumeColumn(settings);
      read = readPlume(result);
      if (probe === null) probe = read.peak.k;
    } catch (error) {
      failed = true;
      console.info('plume figure', error.message);
    }
  }

  function target() {
    const finite = read.excess.filter(Number.isFinite), most = Math.max(0, ...finite), least = Math.min(0, ...finite);
    const flux = Math.max(...result.massFlux), sink = -Math.min(...result.downdraft);
    let heat = 0;
    for (let k = 0; k < result.pressure.length; k++) if (result.pressure[k] >= TOP) heat = Math.max(heat, Math.abs(read.heating[k]), Math.abs(read.drying[k]));
    return {
      excess: [Math.max(-10, Math.min(-2, -room(-least))), Math.max(2, room(most))],
      flux: [-Math.max(0.02, room(sink)), Math.max(0.06, room(flux))],
      heat: Math.min(100, Math.max(30, room(heat))),
    };
  }

  const flat = (g) => [...g.excess, ...g.flux, g.heat];
  function step(dt) {
    if (!goal) return false;
    range ??= goal;
    const f = 1 - Math.exp(-RATE * dt), from = flat(range), to = flat(goal);
    const next = from.map((v, n) => (Math.abs(to[n] - v) < 1e-4 * Math.max(1, Math.abs(to[n])) ? to[n] : v + (to[n] - v) * f));
    range = { excess: next.slice(0, 2), flux: next.slice(2, 4), heat: next[4] };
    if (next.every((v, n) => v === to[n])) fig.play(false);
  }

  function layout(w, h) {
    const narrow = w < WIDE, size = narrow ? 10 : 10.5, axis = narrow ? 34 : 40, right = narrow ? 8 : 10, gap = narrow ? 14 : 14;
    const gaugeH = 64, rows = narrow ? 2 : 1, head = 30, foot = 30;
    const usable = h - gaugeH - rows * (head + foot), rowH = usable / rows;
    const rowsAt = Array.from({ length: rows }, (_, n) => { const top = head + n * (rowH + head + foot); return { top, bottom: top + rowH }; });
    const span = (row, x0, x1) => ({ x0, x1, ...rowsAt[row] });
    let panels;
    if (narrow) {
      const inner = w - axis - right - gap, left = Math.round(inner * 0.55);
      panels = { air: span(0, axis, axis + left), excess: span(0, axis + left + gap, w - right), flux: span(1, axis, axis + left - 20), heat: span(1, axis + left - 20 + gap, w - right) };
    } else {
      const inner = w - axis - right - 3 * gap, unit = inner / 3.4;
      const xs = [axis, axis + unit * 1.1, axis + unit * 1.1 + gap, axis + unit * 1.85 + gap, axis + unit * 1.85 + 2 * gap, axis + unit * 2.5 + 2 * gap, axis + unit * 2.5 + 3 * gap, w - right];
      panels = { air: span(0, xs[0], xs[1]), excess: span(0, xs[2], xs[3]), flux: span(0, xs[4], xs[5]), heat: span(0, xs[6], xs[7]) };
    }
    const ground = result.surfacePressure, scale = Math.log(ground / TOP);
    for (const panel of Object.values(panels)) panel.y = (p) => panel.bottom - Math.log(ground / p) / scale * (panel.bottom - panel.top);
    const gauge = { x0: axis, x1: w - right, top: h - gaugeH + 24, bottom: h - gaugeH + 36 };
    return { narrow, size, axis, panels, gauge, rows: rowsAt };
  }

  function frame(ctx, panel, title, lo, hi, values, format, size, label) {
    const x = (v) => panel.x0 + (v - lo) / (hi - lo) * (panel.x1 - panel.x0);
    ctx.lineWidth = 1;
    for (const p of PRESSURES) { const y = panel.y(p * 100); ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(panel.x0, y); ctx.lineTo(panel.x1, y); ctx.stroke(); }
    for (const v of values) {
      ctx.strokeStyle = v === 0 ? LINE : GRID; ctx.beginPath(); ctx.moveTo(x(v), panel.top); ctx.lineTo(x(v), panel.bottom); ctx.stroke();
      text(ctx, format(v), x(v), panel.bottom + 13, { align: 'center', color: MUTED, size });
    }
    text(ctx, title, panel.x0, panel.top - 12, { color: INK, size: size + 0.5, weight: 500 });
    if (label) text(ctx, label, (panel.x0 + panel.x1) / 2, panel.bottom + 26, { align: 'center', color: MUTED, size });
    return x;
  }

  function pressureAxis(ctx, panel, size) {
    for (const p of PRESSURES) text(ctx, `${p}`, panel.x0 - 6, panel.y(p * 100), { align: 'right', color: MUTED, size });
    text(ctx, 'hPa', panel.x0 - 6, panel.top - 12, { align: 'right', color: MUTED, size });
  }

  function stairs(ctx, panel, x, values, lo, hi) {
    const { interfaces } = result, K = values.length, y = panel.y;
    ctx.beginPath();
    let open = false;
    for (let k = K - 1; k >= 0; k--) {
      const v = values[k];
      if (!Number.isFinite(v) || interfaces[k + 1] < TOP) { open = false; continue; }
      const px = x(clamp(v, lo, hi)), y0 = y(interfaces[k + 1]), y1 = y(Math.max(TOP, interfaces[k]));
      if (open) ctx.lineTo(px, y0); else ctx.moveTo(px, y0);
      ctx.lineTo(px, y1);
      open = true;
    }
  }

  function tag(ctx, taken, words, x, y, align, color, size) {
    ctx.font = `400 ${size}px system-ui, sans-serif`;
    const width = ctx.measureText(words).width, left = align === 'right' ? x - width : x;
    let at = y;
    for (let n = 0; n < 6 && taken.some((b) => left < b.right && left + width > b.left && Math.abs(at - b.y) < size + 2); n++) at += size + 2;
    taken.push({ left, right: left + width, y: at });
    text(ctx, words, x, at, { align, color, size, halo: HALO });
  }

  function overflow(ctx, taken, panel, x, values, lo, hi, color, size) {
    const worst = [-1, -1];
    for (let k = 0; k < values.length; k++) {
      const v = values[k];
      if (!Number.isFinite(v) || result.pressure[k] < TOP || (v >= lo && v <= hi)) continue;
      const y = panel.y(result.pressure[k]), side = v > hi ? 1 : 0, edge = x(side ? hi : lo), out = side ? 1 : -1;
      arrow(ctx, edge - out * 10, y, edge, y, { color, width: 1.5, head: 5 });
      if (worst[side] < 0 || Math.abs(v) > Math.abs(values[worst[side]])) worst[side] = k;
    }
    for (const [side, k] of worst.entries()) {
      if (k < 0) continue;
      const edge = x(side ? hi : lo), out = side ? 1 : -1;
      tag(ctx, taken, minus(Math.round(values[k])), edge - out * 13, panel.y(result.pressure[k]) - 8, side ? 'right' : 'left', color, size);
    }
  }

  function draw(ctx, w, h) {
    if (failed || !result) { text(ctx, failed ? 'This figure needs the model, which could not start here.' : 'Starting the model…', w / 2, h / 2, { align: 'center', color: MUTED, size: 12 }); return; }
    range ??= goal ??= target();
    const g = geometry = layout(w, h), { panels, size, narrow } = g, r = result, K = r.pressure.length, deep = r.type === 'deep';

    const air = panels.air, tx = frame(ctx, air, narrow ? 'air and plume, K' : 'temperature, K', T_RANGE[0], T_RANGE[1], ticks(T_RANGE[0], T_RANGE[1], narrow ? 3 : 4).filter((v) => v >= 190 && v <= 300), (v) => `${v}`, size);
    pressureAxis(ctx, air, size);
    const ex = panels.excess, [elo, ehi] = range.excess, xx = frame(ctx, ex, narrow ? 'plume − air, K' : 'plume minus air, K', elo, ehi, ticks(elo, ehi, narrow ? 3 : 4), (v) => minus(v), size);
    const mf = panels.flux, [flo, fhi] = range.flux, fx = frame(ctx, mf, narrow ? 'mass flux' : 'mass flux, kg/m²/s', flo, fhi, ticks(flo, fhi, 2), (v) => minus(+v.toFixed(2)), size, narrow ? 'kg/m²/s' : null);
    if (narrow) pressureAxis(ctx, mf, size);
    const hp = panels.heat, R = range.heat, hx = frame(ctx, hp, narrow ? 'heating and drying' : 'heating and drying, K per day', -R, R, ticks(-R, R, narrow ? 3 : 4), (v) => minus(v), size, narrow ? 'K per day' : null);

    const clip = (panel, fn) => { ctx.save(); ctx.beginPath(); ctx.rect(panel.x0, panel.top - 2, panel.x1 - panel.x0, panel.bottom - panel.top + 4); ctx.clip(); fn(); ctx.restore(); };

    clip(air, () => {
      ctx.lineJoin = 'round';
      ctx.strokeStyle = VAPOR; ctx.lineWidth = 1.5; ctx.setLineDash([4, 3]); ctx.beginPath();
      let started = false;
      for (let k = K - 1; k >= 0; k--) { if (r.pressure[k] < TOP * 0.9) break; const px = tx(dewPoint(r.vapor[k], r.pressure[k])), py = air.y(r.pressure[k]); if (started) ctx.lineTo(px, py); else { ctx.moveTo(px, py); started = true; } }
      ctx.stroke(); ctx.setLineDash([]);
      ctx.strokeStyle = INK; ctx.lineWidth = 2; ctx.beginPath();
      for (let k = K - 1; k >= 0; k--) { if (r.pressure[k] < TOP * 0.9) break; const px = tx(r.temperature[k]), py = air.y(r.pressure[k]); if (k === K - 1) ctx.moveTo(px, py); else ctx.lineTo(px, py); }
      ctx.stroke();
      if (read.path.length) {
        ctx.strokeStyle = D; ctx.lineWidth = 2; ctx.beginPath();
        read.path.forEach(({ pressure, temperature }, n) => { const px = tx(temperature), py = air.y(pressure); if (n) ctx.lineTo(px, py); else ctx.moveTo(px, py); });
        ctx.stroke();
        const lcl = read.source.lcl;
        ctx.fillStyle = D; ctx.strokeStyle = HALO; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(tx(lcl.temperature), air.y(lcl.pressure), 3.5, 0, Math.PI * 2); ctx.stroke(); ctx.fill();
      }
    });

    clip(ex, () => {
      if (!deep) return;
      for (const [sign, fill] of [[1, FILL], [-1, COLDER]]) {
        ctx.fillStyle = fill; ctx.beginPath();
        for (let k = 0; k < K; k++) {
          const v = read.excess[k];
          if (!(sign * v > 0) || r.interfaces[k + 1] < TOP) continue;
          const y0 = ex.y(r.interfaces[k + 1]), y1 = ex.y(Math.max(TOP, r.interfaces[k]));
          ctx.rect(xx(0), y1, xx(clamp(v, elo, ehi)) - xx(0), y0 - y1);
        }
        ctx.fill();
      }
      ctx.strokeStyle = D; ctx.lineWidth = 2; ctx.lineJoin = 'round';
      stairs(ctx, ex, xx, read.excess, elo, ehi); ctx.stroke();
    });
    overflow(ctx, [], ex, xx, read.excess, elo, ehi, D, size);

    clip(mf, () => {
      ctx.lineJoin = 'round';
      ctx.fillStyle = withAlpha(D, 0.12); ctx.beginPath(); ctx.moveTo(fx(0), mf.y(r.interfaces[K]));
      for (let k = K; k >= 0; k--) { if (r.interfaces[k] < TOP * 0.9) break; ctx.lineTo(fx(r.massFlux[k]), mf.y(r.interfaces[k])); }
      ctx.lineTo(fx(0), mf.y(TOP * 0.9)); ctx.closePath(); ctx.fill();
      for (const [values, dash] of [[r.massFlux, []], [r.downdraft, [4, 3]]]) {
        ctx.strokeStyle = D; ctx.lineWidth = 2; ctx.setLineDash(dash); ctx.beginPath();
        for (let k = K; k >= 0; k--) { if (r.interfaces[k] < TOP * 0.9) break; const px = fx(values[k]), py = mf.y(r.interfaces[k]); if (k === K) ctx.moveTo(px, py); else ctx.lineTo(px, py); }
        ctx.stroke(); ctx.setLineDash([]);
      }
    });
    const peakFlux = r.massFlux.reduce((best, v, k) => (v > r.massFlux[best] ? k : best), 0);
    if (deep) text(ctx, 'up', Math.min(fx(r.massFlux[peakFlux]) + 6, mf.x1 - 14), mf.y(r.interfaces[peakFlux]) - 10, { color: D, size, halo: HALO });
    const deepest = r.downdraft.reduce((best, v, k) => (v < r.downdraft[best] ? k : best), 0);
    if (deep && r.downdraft[deepest] < 0) text(ctx, 'down', Math.max(fx(r.downdraft[deepest]) - 4, mf.x0 + 2), mf.y(r.interfaces[deepest]) - 18, { color: D, size, halo: HALO });

    clip(hp, () => {
      ctx.lineJoin = 'round'; ctx.lineWidth = 2;
      ctx.strokeStyle = VAPOR; stairs(ctx, hp, hx, read.drying, -R, R); ctx.stroke();
      ctx.strokeStyle = D; stairs(ctx, hp, hx, read.heating, -R, R); ctx.stroke();
    });
    const taken = [];
    overflow(ctx, taken, hp, hx, read.drying, -R, R, VAPOR, size);
    overflow(ctx, taken, hp, hx, read.heating, -R, R, D, size);
    const lines = [];
    for (const values of [read.drying, read.heating]) {
      let last = null;
      for (let k = K - 1; k >= 0; k--) {
        if (!Number.isFinite(values[k]) || r.interfaces[k + 1] < TOP) { last = null; continue; }
        const px = hx(clamp(values[k], -R, R)), y0 = hp.y(r.interfaces[k + 1]), y1 = hp.y(Math.max(TOP, r.interfaces[k]));
        if (last !== null) lines.push([Math.min(last, px), Math.max(last, px), y0, y0]);
        lines.push([px, px, y1, y0]);
        last = px;
      }
    }
    const crosses = (left, right, y) => lines.some(([x0, x1, top, bottom]) => x1 >= left - 4 && x0 <= right + 4 && bottom >= y - 0.75 * size && top <= y + 0.75 * size)
      || taken.some((b) => left < b.right && right > b.left && Math.abs(y - b.y) < size + 2);
    const label = (values, color, words, below) => {
      let best = -1;
      for (let k = 0; k < K; k++) if (r.pressure[k] >= TOP && values[k] > 0 && (best < 0 || values[k] > values[best])) best = k;
      if (best < 0) return;
      ctx.font = `400 ${size}px system-ui, sans-serif`;
      const width = ctx.measureText(words).width, px = hx(Math.min(values[best], R)), py = hp.y(r.pressure[best]);
      const spots = [];
      for (const dy of [12, 24, 36, 48, 60, 72]) for (const sign of below ? [1, -1] : [-1, 1]) for (const left of [px + 8, px - 8 - width, hp.x0 + 4]) spots.push([left, py + sign * dy]);
      const fits = ([left, y]) => left >= hp.x0 + 2 && left + width <= hp.x1 - 2 && y >= hp.top + size && y <= hp.bottom - size && !crosses(left, left + width, y);
      const [left, y] = spots.find(fits) ?? [clamp(px + 8, hp.x0 + 2, hp.x1 - 2 - width), py + (below ? 12 : -12)];
      taken.push({ left, right: left + width, y });
      text(ctx, words, left, y, { color, size, halo: HALO });
    };
    label(read.drying, VAPOR, 'drying', true);
    label(read.heating, D, 'heating', false);

    for (const panel of Object.values(panels)) {
      for (const p of [r.base, r.top]) {
        if (p === null || p < TOP) continue;
        const y = panel.y(p);
        ctx.strokeStyle = EDGE; ctx.lineWidth = 1.5; ctx.lineCap = 'round'; ctx.setLineDash([0.5, 4]); ctx.beginPath(); ctx.moveTo(panel.x0, y); ctx.lineTo(panel.x1, y); ctx.stroke(); ctx.setLineDash([]); ctx.lineCap = 'butt';
      }
    }
    for (const [p, words] of [[r.top, 'top'], [r.base, 'cloud base']]) {
      if (p === null || p < TOP) continue;
      const k = r.pressure.reduce((best, q, j) => (Math.abs(Math.log(q / p)) < Math.abs(Math.log(r.pressure[best] / p)) ? j : best), 0);
      const right = tx(r.temperature[k]) < (air.x0 + air.x1) / 2;
      text(ctx, words, right ? air.x1 - 3 : air.x0 + 4, air.y(p) - 8, { align: right ? 'right' : 'left', color: INK, size, halo: HALO });
    }

    const yProbe = r.pressure[probe];
    for (const panel of Object.values(panels)) {
      const y = panel.y(yProbe);
      ctx.strokeStyle = 'rgba(255, 255, 255, 0.35)'; ctx.lineWidth = 1; ctx.setLineDash([2, 4]); ctx.beginPath(); ctx.moveTo(panel.x0, y); ctx.lineTo(panel.x1, y); ctx.stroke(); ctx.setLineDash([]);
    }
    const dot = (panel, px, color) => { ctx.fillStyle = color; ctx.strokeStyle = HALO; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(clamp(px, panel.x0, panel.x1), panel.y(yProbe), 3.5, 0, Math.PI * 2); ctx.stroke(); ctx.fill(); };
    dot(air, tx(r.temperature[probe]), INK);
    if (Number.isFinite(read.excess[probe])) dot(ex, xx(clamp(read.excess[probe], elo, ehi)), D);
    dot(hp, hx(clamp(read.heating[probe], -R, R)), D);
    dot(hp, hx(clamp(read.drying[probe], -R, R)), VAPOR);
    if (!touched) text(ctx, 'drag', air.x0 + 6, air.y(yProbe) - 9, { color: ACCENT, size: size + 0.5, halo: HALO });

    gauge(ctx, g.gauge, size);
    report();
  }

  function gauge(ctx, box, size) {
    const span = Math.log10(RAIN[1] / RAIN[0]), x = (v) => box.x0 + Math.log10(clamp(v, RAIN[0], RAIN[1]) / RAIN[0]) / span * (box.x1 - box.x0);
    const rain = result.rain, mid = (box.top + box.bottom) / 2;
    text(ctx, 'rain at the ground, mm per day', box.x0, box.top - 10, { color: INK, size: size + 0.5, weight: 500 });
    ctx.fillStyle = 'rgba(255, 255, 255, 0.05)'; ctx.strokeStyle = GRID; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.roundRect(box.x0 + 0.5, box.top + 0.5, box.x1 - box.x0 - 1, box.bottom - box.top - 1, 3); ctx.fill(); ctx.stroke();
    if (rain > RAIN[0]) { ctx.fillStyle = VAPOR; ctx.beginPath(); ctx.roundRect(box.x0, box.top, x(rain) - box.x0, box.bottom - box.top, 3); ctx.fill(); }
    for (const v of [0.1, 1, 10, 100, 1000]) {
      const px = x(v);
      ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(px, box.bottom); ctx.lineTo(px, box.bottom + 4); ctx.stroke();
      text(ctx, `${v}`, px, box.bottom + 12, { align: v === 0.1 ? 'left' : v === 1000 ? 'right' : 'center', color: MUTED, size });
    }
    const words = `${rainText(rain)} mm/day`, px = x(rain), inside = px - box.x0 > 90;
    text(ctx, words, inside ? px - 6 : px + 6, mid, { align: inside ? 'right' : 'left', color: inside ? '#151515' : INK, size: size + 0.5, weight: 600 });
  }

  function report() {
    const r = result, deep = r.type === 'deep', cell = (r.cellSize / 1000).toFixed(1);
    const flux = `${r.baseFlux.toFixed(4)} kg/m²/s${r.baseFlux >= read.cap * (1 - 1e-9) ? ', the most the model allows' : ''}`;
    summary.set([
      ['plume', deep ? 'deep: its test parcel’s cloud is deeper than 200 hPa' : r.type === 'shallow' ? 'shallow: it never rains' : 'none'],
      ['cloud base', r.base === null ? '—' : hPa(r.base)],
      ['top', r.top === null ? '—' : `${hPa(r.top)}, ${km(read.topHeight)}`],
      ['CAPE', r.cape === null ? '—' : `${Math.round(r.cape)} J/kg`],
      ['base mass flux', flux],
      ['τ', r.tau === null ? '—' : `${Math.round(r.tau).toLocaleString('en-US')} s, ${(r.tau / 3600).toFixed(1)} h, on a cell ${cell} km across`],
      ['rain', `${rainText(r.rain)} mm/day`],
      ['peak heating', `${read.peak.rate.toFixed(1)} K/day at ${hPa(read.peak.pressure)}`],
    ]);
    const k = probe, rh = Math.round(100 * r.relativeHumidity[k]), drying = -(r.moistening.convection[k] + r.moistening.largeScale[k]);
    const rising = (r.massFlux[k] + r.massFlux[k + 1]) / 2, sinking = -(r.downdraft[k] + r.downdraft[k + 1]) / 2, around = rising - sinking;
    layer.set([
      ['at', `${hPa(r.pressure[k])}, ${km(r.height[k])}`],
      ['air', `${r.temperature[k].toFixed(1)} K, ${rh}% relative humidity`],
      ['plume', Number.isFinite(read.excess[k]) ? `${signed(read.excess[k], 2)} K against the air` : 'not here'],
      ['rising', rising > 0 ? `${rising.toFixed(4)} kg/m²/s` : 'none here'],
      ['downdraft', sinking > 0 ? `${sinking.toFixed(4)} kg/m²/s` : 'none here'],
      ['surrounding air', Math.abs(around) < 5e-5 ? 'still' : `${around > 0 ? 'sinking' : 'rising'} ${Math.abs(around).toFixed(4)} kg/m²/s`],
      ['heating', `${signed(read.heating[k])} K/day`],
      ['drying', `${signed(drying, 2)} g/kg/day, ${signed(read.drying[k])} K/day of latent heat`],
    ]);
  }

  function panelAt({ x, y }) {
    return geometry && Object.values(geometry.panels).find((p) => y >= p.top - 8 && y <= p.bottom + 8 && x >= p.x0 - geometry.axis && x <= p.x1 + 8);
  }

  function pick(point) {
    const panel = result && panelAt(point);
    if (!panel) return false;
    let best = -1;
    for (let k = 0; k < result.pressure.length; k++) {
      if (result.pressure[k] < TOP) continue;
      if (best < 0 || Math.abs(panel.y(result.pressure[k]) - point.y) < Math.abs(panel.y(result.pressure[best]) - point.y)) best = k;
    }
    touched = true;
    probe = best;
    fig.render();
    return true;
  }
  const at = (event) => { const box = fig.canvas.getBoundingClientRect(); return { x: event.clientX - box.left, y: event.clientY - box.top }; };
  const onProbe = (point) => { const panel = result && panelAt(point); return Boolean(panel) && Math.abs(point.y - panel.y(result.pressure[probe])) < 16; };
  fig.pointer({ hit: onProbe, down: pick, move: pick });
  fig.canvas.addEventListener('click', (event) => pick(at(event)));
  fig.canvas.addEventListener('pointermove', (event) => { if (event.pointerType === 'mouse' && !event.buttons) pick(at(event)); });

  compute();
  if (result) goal = target();
  fig.render();
}
