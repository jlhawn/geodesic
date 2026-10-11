import { Figure, slider, choice, readout, legend, caption, text, termColor, clamp, MUTED, GRID, LINE, INK } from '../runtime.module.js';
import { bulkFluxes, WATTS_PER_MM_DAY } from '../surfaceFluxes.module.js';

const WIDE = 520, ZERO = 273.15, LOW = 0, HIGH = 32, SAMPLES = 129, HALO = 'rgba(20, 20, 22, 0.85)', PROBE = 'rgba(255, 255, 255, 0.28)', STORED = '#77736b';
const TOPS = [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60];
const DEFAULT = { sea: 27, humidity: 80, wind: 7, colder: 1 };
const GROUNDS = { all: 'The whole planet’s ground over one model day, W/m²', sea: 'The model’s sea and sea ice over one day, W/m²', land: 'The model’s land over one day, W/m²' };

const number = (v, digits = 1) => (Math.abs(v) < 0.5 * 10 ** -digits ? 0 : v).toLocaleString('en-US', { minimumFractionDigits: digits, maximumFractionDigits: digits });
const signed = (v, digits = 1) => `${v < -0.5 * 10 ** -digits ? '−' : ''}${number(Math.abs(v), digits)}`;
const tickStep = (span, most, steps) => steps.find((s) => span / s <= most) ?? steps[steps.length - 1];

function place(ctx, items, lo, hi, size) {
  ctx.font = `400 ${size}px system-ui, sans-serif`;
  for (const item of items) item.width = Math.max(ctx.measureText(item.name).width, ctx.measureText(item.value).width);
  let edge = lo;
  for (const item of items) { item.left = Math.max(edge, item.anchor - item.width / 2); edge = item.left + item.width + 10; }
  edge = hi;
  for (const item of [...items].reverse()) { item.left = Math.min(item.left, edge - item.width); edge = item.left - 10; }
  for (const item of items) item.center = Math.max(lo, item.left) + item.width / 2;
}

function labels(ctx, items, lo, hi, size, from, to, first) {
  place(ctx, items, lo, hi, size);
  const line = size + 2.5, up = to < from;
  for (const item of items) {
    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(item.anchor, from); ctx.lineTo(item.center, to); ctx.stroke();
    const rows = up ? [[item.value, INK], [item.name, MUTED]] : [[item.name, MUTED], [item.value, INK]];
    rows.forEach(([words, color], n) => text(ctx, words, item.center, first + (up ? -n : n) * line, { align: 'center', color, size }));
  }
}

function bar(ctx, x, y, height, per, parts) {
  for (const [value, color] of parts) {
    const width = Math.max(0, value) * per;
    if (width > 0) {
      const inset = width > 4 ? 1 : 0, w = Math.max(1, width - 2 * inset);
      ctx.fillStyle = color;
      ctx.beginPath(); ctx.roundRect(x + inset, y, w, height, Math.min(3, w / 2)); ctx.fill();
    }
    x += width;
  }
}

export function mountEvaporation(root) {
  const controls = root.querySelector('.controls'), A = termColor('a'), B = termColor('b'), C = termColor('c'), D = termColor('d');
  const state = { ...DEFAULT };
  let here = null, curve = null, goal = 8, top = null, budget = null, failed = false, chart = null, part = 'all';

  const fig = new Figure(root, { height: (w) => (w < WIDE ? 430 : 410), step, draw });
  caption(root, 'Evaporation and sensible heat from an open sea into the air 20 m above it, the height of the model’s lowest layer, computed with the model’s own exchange code. Along the curves only the sea’s temperature changes; the air keeps its humidity, its wind and how much colder it is than the sea. Drag across the chart to set the sea’s temperature. Below is the climate model’s own ground over one day of a run from the state the climate page opens with, with the sea held as a slab, averaged over the whole planet, its sea or its land: the sunlight it absorbs, and where that energy goes.');
  legend(root, [['line', 'evaporation, in mm of water a day and as L·E in W/m²', D], ['line', 'sensible heat, W/m²', C], ['dots', 'the model’s mean evaporation over the ground in the bars', MUTED]]);
  controls.classList.add('two');
  const sea = slider(controls, { label: 'Sea temperature', min: LOW, max: HIGH, step: 0.5, value: state.sea, format: (v) => `${v} °C`, onInput: (v) => change({ sea: v }) });
  slider(controls, { label: 'Air humidity', min: 40, max: 100, step: 1, value: state.humidity, format: (v) => `${v}%`, onInput: (v) => change({ humidity: v }) });
  slider(controls, { label: 'Wind', min: 1, max: 15, step: 0.5, value: state.wind, format: (v) => `${v} m/s`, onInput: (v) => change({ wind: v }) });
  slider(controls, { label: 'Air colder than the sea by', min: 0, max: 5, step: 0.1, value: state.colder, format: (v) => `${v.toFixed(1)} K`, onInput: (v) => change({ colder: v }) });
  const out = readout(controls);
  controls.querySelector('.readout').style.gridColumn = '1 / -1';
  choice(controls, { label: 'Bars', options: [['whole planet', 'all'], ['sea', 'sea'], ['land', 'land']], value: part, onChange: (v) => { part = v; fig.render(); }, span: true });

  const fluxes = (celsius) => bulkFluxes({ seaT: celsius + ZERO, airT: celsius + ZERO - state.colder, rh: state.humidity / 100, wind: state.wind });

  function change(changes) {
    const shape = Object.keys(changes).some((key) => key !== 'sea');
    Object.assign(state, changes);
    if (shape || !curve) {
      curve = { sea: new Float64Array(SAMPLES), E: new Float64Array(SAMPLES), SH: new Float64Array(SAMPLES) };
      let most = 0;
      for (let n = 0; n < SAMPLES; n++) {
        const celsius = LOW + (HIGH - LOW) * n / (SAMPLES - 1), f = fluxes(celsius);
        curve.sea[n] = celsius; curve.E[n] = f.mmPerDay; curve.SH[n] = f.SH;
        most = Math.max(most, f.mmPerDay, f.SH / WATTS_PER_MM_DAY);
      }
      goal = TOPS.find((t) => t >= most * 1.08) ?? TOPS[TOPS.length - 1];
      top ??= goal;
    }
    here = fluxes(state.sea);
    const below = fluxes(state.sea - 0.5).E, above = fluxes(state.sea + 0.5).E;
    const rate = below > 0 ? `${number(here.mmPerDay, 2)} mm/day, rising ${number(100 * (above / below - 1), 1)}% per °C of sea` : here.E > 0 ? `${number(here.mmPerDay, 2)} mm/day` : 'none';
    out.set([
      ['evaporation', rate],
      ['L·E', `${number(here.LE)} W/m²`],
      ['sensible heat', `${signed(here.SH)} W/m²`],
      ['Bowen ratio', here.bowen === null ? 'none' : signed(here.bowen, 3)],
      ['exchange coefficient', `${number(here.CH * 1e3, 2)} × 10⁻³, ${number(here.CHneutral * 1e3, 2)} in neutral air`],
      ['vapor', `${number(here.qSea * 1e3)} g/kg saturated at the sea’s surface, ${number(here.qAir * 1e3)} in the air`],
    ]);
    fig.play(top !== goal);
    fig.render();
  }

  function step(dt) {
    if (top === null) return false;
    const f = 1 - Math.exp(-dt * 10);
    top += (goal - top) * f;
    if (Math.abs(goal - top) < goal * 2e-3) { top = goal; fig.play(false); }
  }

  function plot(ctx, w, size, narrow) {
    const left = narrow ? 30 : 40, right = w - (narrow ? 38 : 46), y0 = 28, y1 = y0 + (narrow ? 206 : 196);
    const x = (c) => left + (c - LOW) / (HIGH - LOW) * (right - left), y = (mm) => y1 - mm / top * (y1 - y0);
    chart = { left, right, y0, y1, x };
    const every = tickStep(top, 5, [0.5, 1, 2, 5, 10, 20]), wattsTop = top * WATTS_PER_MM_DAY, wattsEvery = tickStep(wattsTop, 5, [10, 20, 25, 50, 100, 200, 250, 500]);
    ctx.lineWidth = 1;
    for (let v = 0; v <= top + 1e-9; v += every) {
      ctx.strokeStyle = v ? GRID : LINE; ctx.beginPath(); ctx.moveTo(left, y(v)); ctx.lineTo(right, y(v)); ctx.stroke();
      text(ctx, number(v, every < 1 ? 1 : 0), left - 6, y(v), { align: 'right', color: MUTED, size });
    }
    for (let v = 0; v <= wattsTop + 1e-9; v += wattsEvery) {
      const yy = y(v / WATTS_PER_MM_DAY);
      ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(right, yy); ctx.lineTo(right + 4, yy); ctx.stroke();
      text(ctx, number(v, 0), right + 7, yy, { color: MUTED, size });
    }
    text(ctx, 'mm/day', 2, y0 - 16, { color: MUTED, size });
    text(ctx, 'W/m²', w - 2, y0 - 16, { align: 'right', color: MUTED, size });
    for (let c = 0; c <= HIGH; c += narrow ? 10 : 5) text(ctx, `${c}`, x(c), y1 + 14, { align: 'center', color: MUTED, size });
    text(ctx, 'sea temperature, °C', (left + right) / 2, y1 + 30, { align: 'center', color: MUTED, size });

    const px = x(state.sea), at = (v) => clamp(y(v), y0, y1 + 4), yE = at(here.mmPerDay), ySH = at(here.SH / WATTS_PER_MM_DAY);
    ctx.font = `400 ${size}px system-ui, sans-serif`;
    const [upper, lower] = [{ label: `${number(here.mmPerDay, 2)} mm/day, ${number(here.LE, 0)} W/m²`, at: yE }, { label: `sensible heat ${signed(here.SH)} W/m²`, at: ySH }].sort((a, b) => a.at - b.at);
    upper.y = upper.at - 11;
    lower.y = lower.at + 13 <= y1 - 4 ? lower.at + 13 : lower.at - 11;
    if (lower.y - upper.y < 13) upper.y = lower.y - 13;
    const tags = [upper, lower].map((tag) => {
      const width = ctx.measureText(tag.label).width, room = px - 9 - width >= left + 2, from = room ? px - 9 - width : Math.min(px + 9, right - 2 - width);
      return { ...tag, y: Math.max(tag.y, y0 - 8), from, to: from + width };
    });

    ctx.save();
    ctx.beginPath(); ctx.rect(left, y0 - 4, right - left, y1 - y0 + 10); ctx.clip();
    ctx.strokeStyle = PROBE; ctx.beginPath(); ctx.moveTo(px, y0); ctx.lineTo(px, y1); ctx.stroke();
    let mean = null;
    if (budget) {
      const yy = y(budget[part].evaporation), label = `${{ all: 'planet’s', sea: 'sea’s', land: 'land’s' }[part]} mean ${number(budget[part].evaporation, 2)} mm/day`, width = ctx.measureText(label).width;
      ctx.save(); ctx.strokeStyle = MUTED; ctx.lineWidth = 1.5; ctx.setLineDash([1.5, 3.5]); ctx.lineCap = 'round';
      ctx.beginPath(); ctx.moveTo(left, yy); ctx.lineTo(right, yy); ctx.stroke(); ctx.restore();
      const lines = [curve.E, curve.SH.map((v) => v / WATTS_PER_MM_DAY)];
      const clear = ([from, ly]) => !lines.some((values) => values.some((v, n) => x(curve.sea[n]) >= from - 2 && x(curve.sea[n]) <= from + width + 2 && Math.abs(y(v) - ly) < 8))
        && !tags.some((t) => t.to >= from - 4 && t.from <= from + width + 4 && Math.abs(t.y - ly) < 12) && ly > y0 + 2 && ly < y1 - 4;
      const spots = [-9, 9, -21, -33].flatMap((dy) => [[left + 6, yy + dy], [right - 6 - width, yy + dy]]);
      mean = { label, spot: spots.find(clear) ?? spots[0] };
    }
    const trace = (values, scale, color, width) => {
      ctx.strokeStyle = color; ctx.lineWidth = width; ctx.lineJoin = 'round';
      ctx.beginPath(); values.forEach((v, n) => { const yy = y(v * scale); if (n) ctx.lineTo(x(curve.sea[n]), yy); else ctx.moveTo(x(curve.sea[n]), yy); }); ctx.stroke();
    };
    trace(curve.SH, 1 / WATTS_PER_MM_DAY, C, 2);
    trace(curve.E, 1, D, 2.5);
    if (mean) text(ctx, mean.label, mean.spot[0], mean.spot[1], { color: MUTED, size, halo: HALO });
    ctx.restore();

    for (const [yy, color] of [[ySH, C], [yE, D]]) { ctx.fillStyle = color; ctx.strokeStyle = HALO; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(px, yy, 4.5, 0, Math.PI * 2); ctx.stroke(); ctx.fill(); }
    for (const tag of tags) text(ctx, tag.label, tag.from, tag.y, { color: INK, size, halo: HALO });
  }

  function ground(ctx, w, size, yTop) {
    const pad = w < WIDE ? 10 : 14, x0 = pad + 26, span = w - x0 - pad;
    text(ctx, GROUNDS[part], pad, yTop, { color: INK, size: size + 0.5, weight: 500 });
    if (!budget) { text(ctx, failed ? 'The model’s budget could not be loaded.' : 'Loading the model’s budget…', pad, yTop + 40, { color: MUTED, size }); return; }
    const g = budget[part], drawn = Math.max(0, -g.stored), kept = Math.max(0, g.stored);
    const total = Math.max(g.sunlight + drawn, g.infrared + g.latent + g.sensible + kept), per = span / total, height = 16;
    const yIn = yTop + 46, yOut = yIn + height + 6;
    ctx.fillStyle = 'rgba(255, 255, 255, 0.05)'; ctx.strokeStyle = GRID; ctx.lineWidth = 1;
    for (const yy of [yIn, yOut]) { ctx.beginPath(); ctx.roundRect(x0 + 0.5, yy + 0.5, span - 1, height - 1, 3); ctx.fill(); ctx.stroke(); }
    text(ctx, 'in', x0 - 7, yIn + height / 2, { align: 'right', color: MUTED, size });
    text(ctx, 'out', x0 - 7, yOut + height / 2, { align: 'right', color: MUTED, size });
    bar(ctx, x0, yIn, height, per, [[g.sunlight, A], [drawn, STORED]]);
    bar(ctx, x0, yOut, height, per, [[g.infrared, B], [g.latent, D], [g.sensible, C], [kept, STORED]]);
    const items = (parts) => { let at = x0; return parts.filter(([value]) => value > 0).map(([value, name]) => { const width = value * per, item = { anchor: at + width / 2, name, value: number(value) }; at += width; return item; }); };
    labels(ctx, items([[g.sunlight, 'sunlight absorbed'], [drawn, 'drawn from storage']]), x0, x0 + span, size, yIn - 1, yIn - 7, yIn - 12);
    labels(ctx, items([[g.infrared, 'net infrared'], [g.latent, 'evaporation, L·E'], [g.sensible, 'sensible heat'], [kept, 'stored']]), x0, x0 + span, size, yOut + height + 1, yOut + height + 7, yOut + height + 14);
  }

  function draw(ctx, w, h) {
    const narrow = w < WIDE, size = narrow ? 10.5 : 11;
    if (!curve) return;
    plot(ctx, w, size, narrow);
    ground(ctx, w, size, chart.y1 + 58);
  }

  const setSea = (event) => {
    if (!chart) return;
    const r = fig.canvas.getBoundingClientRect(), c = clamp(LOW + (event.clientX - r.left - chart.left) / (chart.right - chart.left) * (HIGH - LOW), LOW, HIGH);
    const v = Math.round(c * 2) / 2;
    if (v !== state.sea) { sea.value = v; change({ sea: v }); }
  };
  fig.canvas.style.touchAction = 'pan-y';
  fig.canvas.addEventListener('pointerdown', (event) => {
    if (!chart) return;
    const r = fig.canvas.getBoundingClientRect(), py = event.clientY - r.top;
    if (py < chart.y0 - 10 || py > chart.y1 + 20) return;
    fig.canvas.setPointerCapture(event.pointerId);
    setSea(event);
  });
  fig.canvas.addEventListener('pointermove', (event) => { if (fig.canvas.hasPointerCapture(event.pointerId)) setSea(event); });
  for (const end of ['pointerup', 'pointercancel']) fig.canvas.addEventListener(end, (event) => { if (fig.canvas.hasPointerCapture(event.pointerId)) fig.canvas.releasePointerCapture(event.pointerId); });

  change({});
  fetch(new URL('../data/surface_budget.json', import.meta.url)).then((response) => {
    if (!response.ok) throw new Error(`${response.status} for the surface budget`);
    return response.json();
  }).then((data) => { budget = data; fig.render(); }).catch((error) => { failed = true; console.info('evaporation figure', error.message); fig.render(); });
}
