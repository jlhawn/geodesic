import { Figure, choice, slider, readout, legend, caption, text, termColor, clamp, ACCENT, MUTED, GRID, LINE } from '../runtime.module.js';
import { loadFrames } from '../frames.module.js';
import { DATA_URL, transportSplit } from '../oceanRun.module.js';

const D = Math.PI / 180, PW_SPAN = 7, SV_SPAN = 70, HALO = 'rgba(21, 21, 21, 0.85)', PROBE = 'rgba(255, 255, 255, 0.28)';
const SEA = 'rgb(120, 140, 255)', STORED = 'rgba(255, 255, 255, 0.75)', PART7 = 'rgba(255, 255, 255, 0.6)', PART6 = 'rgba(255, 255, 255, 0.3)';
const TOP = { sea: { range: [-60, 60], ticks: [-60, -30, 0, 30, 60] }, surface: { range: [-25, 200], ticks: [0, 50, 100, 150, 200] } };
const VIEWS = { heat: { top: 'sea', bottom: 'heat' }, surface: { top: 'surface', bottom: 'heat' }, flow: { top: 'sea', bottom: 'flow' } };
const minus = (v) => `${v}`.replace('-', '−');
const place = (lat) => (Math.abs(lat) < 0.05 ? '0°' : `${Math.abs(lat).toFixed(Math.abs(lat) % 1 ? 1 : 0)}°${lat > 0 ? 'N' : 'S'}`);
const carried = (pw, digits = 2) => (Math.abs(pw) < 0.5 * 10 ** -digits ? 'none' : `${Math.abs(pw).toFixed(digits)} PW ${pw > 0 ? 'north' : 'south'}`);
const flowing = (sv) => (Math.abs(sv) < 0.05 ? 'none' : `${Math.abs(sv).toFixed(1)} Sv ${sv > 0 ? 'north' : 'south'}`);
const flux = (v) => `${v < 0 ? '−' : v > 0 ? '+' : ''}${Math.abs(v).toFixed(1)} W/m²`;
const plain = (v) => `${v.toFixed(1)} W/m²`;

function withAlpha(color, alpha) {
  const hex = /^#([0-9a-f]{6})$/i.exec(color), rgb = /^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/.exec(color);
  if (hex) return `rgba(${[0, 2, 4].map((n) => parseInt(hex[1].slice(n, n + 2), 16)).join(', ')}, ${alpha})`;
  return rgb ? `rgba(${rgb[1]}, ${rgb[2]}, ${rgb[3]}, ${alpha})` : color;
}

function banded(values, lat) {
  const x = clamp((lat + 90) / 5 - 0.5, 0, values.length - 1), b = Math.min(values.length - 2, Math.floor(x)), t = x - b;
  return values[b] + t * (values[b + 1] - values[b]);
}

export function mountSeaheat(root) {
  const controls = root.querySelector('.controls');
  const sun = termColor('a'), infrared = termColor('b'), sensible = termColor('c'), latentColor = termColor('d'), totalColor = termColor('e');
  const surplus = withAlpha(SEA, 0.42), deficit = 'rgba(255, 255, 255, 0.16)', seaFill = withAlpha(SEA, 0.14);
  let run = null, view = 'heat', probe = -10, geometry = null, failed = false;

  const fig = new Figure(root, { height: (w) => (w < 520 ? 440 : 480), draw });
  caption(root, 'Above, the year’s average of the heat going into the sea through its surface, in bands 5° of latitude wide, from the run of the last section, and of the heat the sea stored or gave up; the shading between them is what its currents must carry away or bring. Below, the heat carried north across each latitude, added up from the South Pole, in all and by each carrier, and the sea’s own flow added up class by class, relative to water at 0\u00a0°C. The latitude scale shrinks toward the poles so that equal widths hold equal areas of the planet. Each part is corrected for the heat it stored over the year, so the curves include the heat the polar sea is storing. The heat the run’s books cannot place is spread evenly over the globe. The sea’s overturning is sampled every six hours, and the samples miss up to 8\u00a0Sv, so its layers do not quite cancel. Poleward of 45° the flow of its densest classes is the densest water of the run’s start draining toward the equator. Tap or drag across the chart to read one latitude.');
  const below = legend(root, []), above = legend(root, []);
  const keys = {
    sea: [['line', 'above: heat into the sea through its surface', SEA], ['dash', 'heat the sea stored, below zero where it gave heat up', STORED], ['bar', 'heat left for the currents to carry away', surplus], ['bar', 'heat the currents bring', deficit]],
    surface: [['line', 'above: sunlight absorbed by the sea', sun], ['line', 'infrared it sends up, less what comes back', infrared], ['line', 'sensible heat', sensible], ['line', 'latent heat of the water it evaporates', latentColor], ['line', 'heat into the sea, what is left', totalColor]],
    heat: [['line', 'below: heat carried north, in all', totalColor], ['line', 'by the air, as the heat of dry air', sensible], ['line', 'by the air, as latent heat in its vapor', latentColor], ['line', 'by the sea', SEA], ['dots', 'by the sea, its flow added up class by class', SEA], ['dash', 'in all by Part 7’s air, over a sea that did not move', PART7], ['faint', 'in all on Part 6’s dry sphere', PART6]],
    flow: [['line', 'below: water carried north by the sea’s mixed layer', ACCENT], ['line', 'by the classes lighter than 1026 kg/m³', SEA], ['dash', 'by the classes from 1026 to 1026.8 kg/m³', SEA], ['dots', 'by the classes from 1026.8 kg/m³ down', SEA]],
  };
  const relabel = () => { above.set(keys[VIEWS[view].top]); below.set(keys[VIEWS[view].bottom]); };
  relabel();
  choice(controls, { label: 'Show', options: [['heat carried', 'heat'], ['the surface’s budget', 'surface'], ['the sea’s overturning', 'flow']], value: view, onChange: (v) => { view = v; relabel(); fig.render(); }, span: true });
  const latitude = slider(controls, { label: 'Latitude', min: -90, max: 90, step: 1, value: probe, format: place, onInput: (v) => { probe = v; fig.render(); }, span: true });
  const upper = readout(controls), lower = readout(controls), global = readout(controls);

  function draw(ctx, w, h) {
    if (!run) { text(ctx, failed ? 'The run could not be loaded.' : 'Loading the run…', w / 2, h / 2, { align: 'center', color: MUTED }); return; }
    const { top: topView, bottom: bottomView } = VIEWS[view], scale = TOP[topView], span = bottomView === 'flow' ? SV_SPAN : PW_SPAN;
    const narrow = w < 520, left = narrow ? 34 : 44, right = narrow ? 10 : 16, size = narrow ? 10 : 11;
    const top = 22, split = Math.round(h * 0.4), low = split + 34, bottom = h - 24;
    const x = (lat) => left + (Math.sin(lat * D) + 1) / 2 * (w - left - right);
    const y = (v) => split - (v - scale.range[0]) / (scale.range[1] - scale.range[0]) * (split - top);
    const z = (v) => (low + bottom) / 2 - v / span * (bottom - low) / 2;
    geometry = { left, right, w };

    ctx.lineWidth = 1;
    for (const v of scale.ticks) {
      ctx.strokeStyle = v ? GRID : LINE; ctx.beginPath(); ctx.moveTo(left, y(v)); ctx.lineTo(w - right, y(v)); ctx.stroke();
      text(ctx, minus(v), left - 6, y(v), { align: 'right', color: MUTED, size });
    }
    text(ctx, 'W/m²', left - 6, top - 12, { align: 'right', color: MUTED, size });
    for (const v of bottomView === 'flow' ? [-60, -30, 0, 30, 60] : [-6, -4, -2, 0, 2, 4, 6]) {
      ctx.strokeStyle = v ? GRID : LINE; ctx.beginPath(); ctx.moveTo(left, z(v)); ctx.lineTo(w - right, z(v)); ctx.stroke();
      text(ctx, minus(v), left - 6, z(v), { align: 'right', color: MUTED, size });
    }
    text(ctx, bottomView === 'flow' ? 'Sv' : 'PW', left - 6, low - 12, { align: 'right', color: MUTED, size });
    text(ctx, 'northward', left + 4, z(span) + 2, { color: MUTED, size, baseline: 'top' });
    text(ctx, 'southward', w - right - 4, z(-span) - 2, { align: 'right', color: MUTED, size, baseline: 'bottom' });
    for (const [lat, label] of [[-60, '60°S'], [-30, '30°S'], [0, '0°'], [30, '30°N'], [60, '60°N']]) {
      text(ctx, label, x(lat), split + 13, { align: 'center', color: MUTED, size });
      text(ctx, label, x(lat), bottom + 13, { align: 'center', color: MUTED, size });
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(x(lat), top); ctx.lineTo(x(lat), split); ctx.moveTo(x(lat), low); ctx.lineTo(x(lat), bottom); ctx.stroke();
    }

    const samples = [];
    for (let lat = -90; lat <= 90; lat += 0.5) samples.push(run.at(lat));
    const curve = (pick, scaleOf, color, width, dash = null) => {
      ctx.save(); ctx.strokeStyle = color; ctx.lineWidth = width; ctx.lineJoin = 'round'; ctx.lineCap = dash && dash[0] < 2 ? 'round' : 'butt';
      if (dash) ctx.setLineDash(dash);
      ctx.beginPath();
      samples.forEach((s, n) => (n ? ctx.lineTo(x(s.lat), scaleOf(pick(s))) : ctx.moveTo(x(s.lat), scaleOf(pick(s)))));
      ctx.stroke(); ctx.restore();
    };

    if (topView === 'sea') {
      for (let n = 0; n + 1 < samples.length; n++) {
        const a = samples[n], b = samples[n + 1];
        ctx.fillStyle = a.input - a.stored + b.input - b.stored >= 0 ? surplus : deficit;
        ctx.beginPath(); ctx.moveTo(x(a.lat), y(a.input)); ctx.lineTo(x(b.lat), y(b.input)); ctx.lineTo(x(b.lat), y(b.stored)); ctx.lineTo(x(a.lat), y(a.stored)); ctx.closePath(); ctx.fill();
      }
      curve((s) => s.stored, y, STORED, 1.5, [4, 3]);
      curve((s) => s.input, y, SEA, 2);
    } else {
      curve((s) => s.sensibleHeat, y, sensible, 2);
      curve((s) => s.infrared, y, infrared, 2);
      curve((s) => s.latentHeat, y, latentColor, 2);
      curve((s) => s.sunlight, y, sun, 2);
      curve((s) => s.input, y, totalColor, 2.5);
    }

    if (bottomView === 'heat') {
      ctx.save(); ctx.fillStyle = seaFill; ctx.beginPath(); ctx.moveTo(x(-90), z(0));
      samples.forEach((s) => ctx.lineTo(x(s.lat), z(s.sea)));
      ctx.lineTo(x(90), z(0)); ctx.closePath(); ctx.fill(); ctx.restore();
      curve((s) => s.part6, z, PART6, 1);
      curve((s) => s.part7, z, PART7, 1.5, [4, 4]);
      curve((s) => s.dry, z, sensible, 2);
      curve((s) => s.latent, z, latentColor, 2);
      curve((s) => s.sea, z, SEA, 2);
      curve((s) => s.direct, z, SEA, 2.5, [0.5, 5]);
      curve((s) => s.total, z, totalColor, 2.5);
      for (const sign of [1, -1]) {
        const peak = run.split.peak('total', sign);
        text(ctx, `${Math.abs(peak.transport).toFixed(1)} PW`, x(peak.latitude), z(peak.transport) - sign * 11, { align: 'center', color: totalColor, size, halo: HALO });
      }
    } else {
      curve((s) => s.groups[3], z, SEA, 2.5, [0.5, 5]);
      curve((s) => s.groups[2], z, SEA, 1.6, [5, 4]);
      curve((s) => s.groups[1], z, SEA, 2);
      curve((s) => s.groups[0], z, ACCENT, 2);
    }

    const at = run.at(probe), px = x(probe);
    ctx.strokeStyle = PROBE; ctx.beginPath(); ctx.moveTo(px, top); ctx.lineTo(px, split); ctx.moveTo(px, low); ctx.lineTo(px, bottom); ctx.stroke();
    const dots = topView === 'sea' ? [[at.input, SEA, y], [at.stored, STORED, y]] : [[at.sunlight, sun, y], [at.infrared, infrared, y], [at.sensibleHeat, sensible, y], [at.latentHeat, latentColor, y], [at.input, totalColor, y]];
    if (bottomView === 'heat') dots.push([at.dry, sensible, z], [at.latent, latentColor, z], [at.sea, SEA, z], [at.total, totalColor, z]);
    else dots.push(...at.groups.map((v, g) => [v, g ? SEA : ACCENT, z]));
    for (const [v, color, scaleOf] of dots) { ctx.fillStyle = color; ctx.beginPath(); ctx.arc(px, scaleOf(v), 4, 0, Math.PI * 2); ctx.fill(); }

    upper.set([
      ['latitude', place(probe)],
      ...(topView === 'sea'
        ? [['heat into the sea', flux(at.input)], ['stored', flux(at.stored)], ['left for the currents', flux(at.input - at.stored)]]
        : [['sunlight absorbed', plain(at.sunlight)], ['infrared up, less what comes back', plain(at.infrared)], ['sensible heat', plain(at.sensibleHeat)], ['latent heat', plain(at.latentHeat)], ['into the sea', flux(at.input)]]),
    ]);
    lower.set(bottomView === 'heat'
      ? [['heat carried', carried(at.total)], ['by the air, as the heat of dry air', carried(at.dry)], ['as latent heat', carried(at.latent)], ['by the sea', carried(at.sea)], ['by its classes, added up', carried(at.direct)]]
      : [['water carried by the mixed layer', flowing(at.groups[0])], ['lighter than 1026', flowing(at.groups[1])], ['1026 to 1026.8', flowing(at.groups[2])], ['from 1026.8', flowing(at.groups[3])], ['heat carried by the sea', carried(at.sea)], ['by its classes, added up', carried(at.direct)]]);
  }

  function summarize() {
    const { split } = run, most = (which) => [-1, 1].map((sign) => { const p = split.peak(which, sign); return `${carried(p.transport, 1)} at ${place(p.latitude)}`; }).join(', '), equator = split.at(0);
    global.set([
      ['most heat carried, in all', most('total')], ['by the air', most('air')], ['by the sea', most('sea')],
      ['across the equator', `${carried(equator.total, 1)}: by the sea ${carried(equator.sea, 1)}, by the air ${carried(equator.air, 1)}`],
      ['heat the run’s books cannot place', `${split.books.global.residual.toFixed(1)} W/m², spread evenly`],
    ]);
  }

  const pick = (event) => {
    if (!geometry) return;
    const r = fig.canvas.getBoundingClientRect(), { left, right, w } = geometry;
    const s = clamp(2 * (event.clientX - r.left - left) / (w - left - right) - 1, -1, 1), lat = Math.round(Math.asin(s) / D);
    if (lat !== probe) { probe = lat; latitude.value = lat; fig.render(); }
  };
  fig.canvas.style.touchAction = 'pan-y';
  fig.canvas.addEventListener('pointerdown', pick);
  fig.canvas.addEventListener('pointermove', (event) => { if (event.pointerType === 'mouse' || event.buttons) pick(event); });

  loadFrames(DATA_URL).then(({ header }) => {
    const split = transportSplit(header), { books } = split, budget = header.zonalBudget;
    const bands = {
      input: books.seaInput, stored: books.sea, sunlight: budget.surfaceShortwave, sensibleHeat: budget.sensibleHeat, latentHeat: budget.latentHeat,
      infrared: budget.surfaceShortwave.map((v, b) => v - budget.sensibleHeat[b] - budget.latentHeat[b] - budget.netSurfaceHeat[b]),
    };
    run = {
      split,
      at(lat) {
        const edge = split.at(lat), out = { lat, groups: edge.groups.map((g) => g.mass) };
        for (const name of ['total', 'sea', 'air', 'latent', 'dry', 'direct', 'part7', 'part6']) out[name] = edge[name];
        for (const [name, values] of Object.entries(bands)) out[name] = banded(values, lat);
        return out;
      },
    };
    summarize(); fig.render();
  }).catch((error) => { failed = true; console.info('seaheat figure', error.message); fig.render(); });
}
