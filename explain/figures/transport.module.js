import { Figure, slider, readout, legend, caption, text, termColor, clamp, ACCENT, MUTED, GRID, LINE } from '../runtime.module.js';
import { loadFrames } from '../frames.module.js';
import { Grid } from '../../js/grid.module.js';
import { DATA_URL, heatTransport, storedHeat } from '../sunlitRun.module.js';

const D = Math.PI / 180, FLUX = [100, 400], PW_SPAN = 1.5, HALO = 'rgba(21, 21, 21, 0.85)', PROBE = 'rgba(255, 255, 255, 0.28)';
const place = (lat) => (Math.abs(lat) < 0.05 ? '0°' : `${Math.abs(lat).toFixed(Math.abs(lat) % 1 ? 1 : 0)}°\u00a0${lat > 0 ? 'N' : 'S'}`);
const signed = (v, digits = 1) => `${v >= 0 ? '+' : '−'}${Math.abs(v).toFixed(digits)}`;
const carried = (pw) => (Math.abs(pw) < 0.005 ? 'none' : `${Math.abs(pw).toFixed(2)}\u00a0PW ${pw > 0 ? 'northward' : 'southward'}`);
const flux = (v, digits = 1) => `${v.toFixed(digits)}\u00a0W/m²`;

function withAlpha(color, alpha) {
  const hex = /^#([0-9a-f]{6})$/i.exec(color);
  return hex ? `rgba(${[0, 2, 4].map((n) => parseInt(hex[1].slice(n, n + 2), 16)).join(', ')}, ${alpha})` : color;
}

export function mountTransport(root) {
  const controls = root.querySelector('.controls');
  const sun = termColor('a'), infrared = termColor('b'), kept = withAlpha(sun, 0.45), lost = withAlpha(infrared, 0.45);
  let run = null, probe = 40, geometry = null, failed = false;

  const fig = new Figure(root, { height: (w) => (w < 520 ? 400 : 440), draw });
  caption(root, 'The year’s average of the sunlight absorbed and the infrared leaving the top of the air, in bands 5° of latitude wide, from the run above. The latitude scale shrinks toward the poles so that equal widths hold equal areas of the planet. Below them is the heat that has to cross each latitude to keep every band in balance, added up from the South Pole. It is corrected for the heat the ground and air stored or shed over the year, and the 0.7\u00a0W/m² that the stored heat and the top-of-air budget leave unaccounted for is spread evenly over the globe. Tap or drag across the chart to read one latitude.');
  legend(root, [['line', 'sunlight absorbed', sun], ['line', 'infrared leaving the top', infrared], ['bar', 'more sunlight absorbed than infrared sent out', kept], ['bar', 'more infrared sent out than sunlight absorbed', lost], ['line', 'heat carried north across each latitude, south below zero', ACCENT]]);
  const latitude = slider(controls, { label: 'Latitude', min: -90, max: 90, step: 1, value: probe, format: place, onInput: (v) => { probe = v; fig.render(); }, span: true });
  const local = readout(controls), global = readout(controls);

  function draw(ctx, w, h) {
    if (!run) { text(ctx, failed ? 'The run could not be loaded.' : 'Loading the run…', w / 2, h / 2, { align: 'center', color: MUTED }); return; }
    const narrow = w < 520, left = narrow ? 40 : 50, right = narrow ? 10 : 16, size = narrow ? 10 : 11;
    const top = 22, split = Math.round(h * 0.56), lower = split + 34, bottom = h - 24;
    const x = (lat) => left + (Math.sin(lat * D) + 1) / 2 * (w - left - right);
    const y = (v) => split - (v - FLUX[0]) / (FLUX[1] - FLUX[0]) * (split - top);
    const z = (pw) => (lower + bottom) / 2 - pw / PW_SPAN * (bottom - lower) / 2;
    geometry = { left, right, w };

    ctx.lineWidth = 1;
    for (let v = FLUX[0]; v <= FLUX[1]; v += 100) {
      ctx.strokeStyle = v > FLUX[0] ? GRID : LINE; ctx.beginPath(); ctx.moveTo(left, y(v)); ctx.lineTo(w - right, y(v)); ctx.stroke();
      text(ctx, `${v}`, left - 6, y(v), { align: 'right', color: MUTED, size });
    }
    text(ctx, 'W/m²', left - 6, top - 12, { align: 'right', color: MUTED, size });
    for (const v of [-1, 0, 1]) {
      ctx.strokeStyle = v ? GRID : LINE; ctx.beginPath(); ctx.moveTo(left, z(v)); ctx.lineTo(w - right, z(v)); ctx.stroke();
      text(ctx, `${v}`, left - 6, z(v), { align: 'right', color: MUTED, size });
    }
    text(ctx, 'PW', left - 6, lower - 12, { align: 'right', color: MUTED, size });
    text(ctx, 'northward', left + 4, z(PW_SPAN) + 2, { color: MUTED, size, baseline: 'top' });
    text(ctx, 'southward', w - right - 4, z(-PW_SPAN) - 2, { align: 'right', color: MUTED, size, baseline: 'bottom' });
    for (const [lat, label] of [[-60, '60°S'], [-30, '30°S'], [0, '0°'], [30, '30°N'], [60, '60°N']]) {
      text(ctx, label, x(lat), split + 13, { align: 'center', color: MUTED, size });
      text(ctx, label, x(lat), bottom + 13, { align: 'center', color: MUTED, size });
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(x(lat), top); ctx.lineTo(x(lat), split); ctx.moveTo(x(lat), lower); ctx.lineTo(x(lat), bottom); ctx.stroke();
    }

    const samples = [];
    for (let lat = -90; lat <= 90; lat += 0.5) samples.push({ lat, ...run.at(lat) });
    for (let n = 0; n + 1 < samples.length; n++) {
      const a = samples[n], b = samples[n + 1];
      ctx.fillStyle = a.absorbed - a.outgoing + b.absorbed - b.outgoing >= 0 ? kept : lost;
      ctx.beginPath(); ctx.moveTo(x(a.lat), y(a.absorbed)); ctx.lineTo(x(b.lat), y(b.absorbed)); ctx.lineTo(x(b.lat), y(b.outgoing)); ctx.lineTo(x(a.lat), y(a.outgoing)); ctx.closePath(); ctx.fill();
    }
    const curve = (key, scale, color, width) => {
      ctx.save(); ctx.strokeStyle = color; ctx.lineWidth = width; ctx.lineJoin = 'round'; ctx.beginPath();
      samples.forEach((s, n) => (n ? ctx.lineTo(x(s.lat), scale(s[key])) : ctx.moveTo(x(s.lat), scale(s[key]))));
      ctx.stroke(); ctx.restore();
    };
    curve('outgoing', y, infrared, 2);
    curve('absorbed', y, sun, 2);
    ctx.save(); ctx.fillStyle = withAlpha(ACCENT, 0.12); ctx.beginPath(); ctx.moveTo(x(-90), z(0));
    samples.forEach((s) => ctx.lineTo(x(s.lat), z(s.transport)));
    ctx.lineTo(x(90), z(0)); ctx.closePath(); ctx.fill(); ctx.restore();
    curve('transport', z, ACCENT, 2);

    for (const sign of [1, -1]) {
      const peak = run.peak(sign);
      text(ctx, `${Math.abs(peak.transport).toFixed(2)} PW`, x(peak.latitude), z(peak.transport) - sign * 10, { align: 'center', color: ACCENT, size, halo: HALO });
    }

    const at = run.at(probe), px = x(probe);
    ctx.strokeStyle = PROBE; ctx.beginPath(); ctx.moveTo(px, top); ctx.lineTo(px, split); ctx.moveTo(px, lower); ctx.lineTo(px, bottom); ctx.stroke();
    for (const [v, color, scale] of [[at.absorbed, sun, y], [at.outgoing, infrared, y], [at.transport, ACCENT, z]]) { ctx.fillStyle = color; ctx.beginPath(); ctx.arc(px, scale(v), 4, 0, Math.PI * 2); ctx.fill(); }

    local.set([['latitude', place(probe)], ['sunlight absorbed', flux(at.absorbed)], ['infrared out', flux(at.outgoing)], ['heat stored over the year', `${signed(at.stored)}\u00a0W/m²`], ['net', `${signed(at.net)}\u00a0W/m²`], ['heat carried', carried(at.transport)]]);
  }

  function summarize() {
    const north = run.peak(1), south = run.peak(-1);
    global.set([
      ['whole planet: sunlight absorbed', flux(run.absorbedMean)], ['infrared out', flux(run.outgoingMean)], ['heat stored', `${signed(run.storedMean, 2)}\u00a0W/m²`],
      ['net, after the heat stored, turns from surplus to loss at', run.crossings.map(place).join(' and ')],
      ['most heat carried', `${carried(south.transport)} at ${place(south.latitude)}, ${carried(north.transport)} at ${place(north.latitude)}`],
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

  loadFrames(DATA_URL).then((frames) => {
    const cells = [...new Grid(frames.header.N)].map((c) => ({ lat: Math.asin(c.centerVertex.z) / D, area: c.area }));
    run = heatTransport(frames.header.zonalBudget, storedHeat(frames, cells).total);
    summarize(); fig.render();
  }).catch((error) => { failed = true; console.info('transport figure', error.message); fig.render(); });
}
