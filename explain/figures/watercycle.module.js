import { Figure, slider, readout, legend, caption, text, termColor, clamp, MUTED, GRID, LINE, VAPOR } from '../runtime.module.js';
import { loadFrames } from '../frames.module.js';
import { DATA_URL, waterTransport, recordedYear } from '../wetRun.module.js';

const D = Math.PI / 180, WATER = [0, 5], PW_SPAN = 6, HALO = 'rgba(21, 21, 21, 0.85)', PROBE = 'rgba(255, 255, 255, 0.28)', PART6 = 'rgba(255, 255, 255, 0.45)';
const place = (lat) => (Math.abs(lat) < 0.05 ? '0°' : `${Math.abs(lat).toFixed(Math.abs(lat) % 1 ? 1 : 0)}°\u00a0${lat > 0 ? 'N' : 'S'}`);
const carried = (pw) => (Math.abs(pw) < 0.05 ? 'none' : `${Math.abs(pw).toFixed(1)}\u00a0PW ${pw > 0 ? 'north' : 'south'}`);
const mm = (v) => `${v.toFixed(2)}\u00a0mm a day`;

function withAlpha(color, alpha) {
  const hex = /^#([0-9a-f]{6})$/i.exec(color), rgb = /^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/.exec(color);
  if (hex) return `rgba(${[0, 2, 4].map((n) => parseInt(hex[1].slice(n, n + 2), 16)).join(', ')}, ${alpha})`;
  return rgb ? `rgba(${rgb[1]}, ${rgb[2]}, ${rgb[3]}, ${alpha})` : color;
}

export function mountWatercycle(root) {
  const controls = root.querySelector('.controls');
  const sensible = termColor('c'), latentColor = termColor('d'), totalColor = termColor('e');
  const source = withAlpha(latentColor, 0.4), sink = withAlpha(VAPOR, 0.45);
  let run = null, year = null, probe = -15, geometry = null, failed = false;

  const fig = new Figure(root, { height: (w) => (w < 520 ? 420 : 460), draw });
  caption(root, 'Above, the year’s average evaporation from the sea and rain back into it, in bands 5° of latitude wide, from the run of the last section. Below, the heat the air carries across each latitude, added up from the South Pole: in all, as the latent heat of its vapor, and as the rest, the heat of the dry air. The latitude scale shrinks toward the poles so that equal widths hold equal areas of the planet. The total is corrected for the heat the sea, the ice and the air stored or shed over the year, and the latent part for the vapor the air gained. The 1.0\u00a0W/m² that the stored heat and the top-of-air budget leave unaccounted for is spread evenly over the globe. Tap or drag across the chart to read one latitude.');
  legend(root, [['line', 'below: heat carried north, in all', totalColor], ['line', 'as latent heat in the vapor', latentColor], ['line', 'as the heat of dry air', sensible], ['faint', 'in all on Part 6’s dry sphere', PART6]]);
  legend(root, [['line', 'above: evaporation', latentColor], ['line', 'rain and snow', VAPOR], ['bar', 'more evaporation than rain: the sea loses water', source], ['bar', 'more rain than evaporation: the sea gains water', sink]]);
  const latitude = slider(controls, { label: 'Latitude', min: -90, max: 90, step: 1, value: probe, format: place, onInput: (v) => { probe = v; fig.render(); }, span: true });
  const local = readout(controls), global = readout(controls);

  function draw(ctx, w, h) {
    if (!run) { text(ctx, failed ? 'The run could not be loaded.' : 'Loading the run…', w / 2, h / 2, { align: 'center', color: MUTED }); return; }
    const narrow = w < 520, left = narrow ? 34 : 44, right = narrow ? 10 : 16, size = narrow ? 10 : 11;
    const top = 22, split = Math.round(h * 0.4), lower = split + 34, bottom = h - 24;
    const x = (lat) => left + (Math.sin(lat * D) + 1) / 2 * (w - left - right);
    const y = (v) => split - (v - WATER[0]) / (WATER[1] - WATER[0]) * (split - top);
    const z = (pw) => (lower + bottom) / 2 - pw / PW_SPAN * (bottom - lower) / 2;
    geometry = { left, right, w };

    ctx.lineWidth = 1;
    for (let v = WATER[0]; v <= WATER[1]; v += 1) {
      ctx.strokeStyle = v > WATER[0] ? GRID : LINE; ctx.beginPath(); ctx.moveTo(left, y(v)); ctx.lineTo(w - right, y(v)); ctx.stroke();
      text(ctx, `${v}`, left - 6, y(v), { align: 'right', color: MUTED, size });
    }
    text(ctx, 'mm a day', left - 6, top - 12, { align: 'left', color: MUTED, size });
    for (const v of [-4, -2, 0, 2, 4]) {
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
      ctx.fillStyle = a.evaporation - a.precipitation + b.evaporation - b.precipitation >= 0 ? source : sink;
      ctx.beginPath(); ctx.moveTo(x(a.lat), y(a.evaporation)); ctx.lineTo(x(b.lat), y(b.evaporation)); ctx.lineTo(x(b.lat), y(b.precipitation)); ctx.lineTo(x(a.lat), y(a.precipitation)); ctx.closePath(); ctx.fill();
    }
    const curve = (key, scale, color, width) => {
      ctx.save(); ctx.strokeStyle = color; ctx.lineWidth = width; ctx.lineJoin = 'round'; ctx.beginPath();
      samples.forEach((s, n) => (n ? ctx.lineTo(x(s.lat), scale(s[key])) : ctx.moveTo(x(s.lat), scale(s[key]))));
      ctx.stroke(); ctx.restore();
    };
    curve('precipitation', y, VAPOR, 2);
    curve('evaporation', y, latentColor, 2);
    curve('part6', z, PART6, 1);
    curve('dry', z, sensible, 2);
    curve('latent', z, latentColor, 2);
    curve('total', z, totalColor, 2.5);

    for (const sign of [1, -1]) {
      const peak = run.peak('total', sign);
      text(ctx, `${Math.abs(peak.transport).toFixed(1)} PW`, x(peak.latitude), z(peak.transport) - sign * 11, { align: 'center', color: totalColor, size, halo: HALO });
    }

    const at = run.at(probe), px = x(probe);
    ctx.strokeStyle = PROBE; ctx.beginPath(); ctx.moveTo(px, top); ctx.lineTo(px, split); ctx.moveTo(px, lower); ctx.lineTo(px, bottom); ctx.stroke();
    for (const [v, color, scale] of [[at.evaporation, latentColor, y], [at.precipitation, VAPOR, y], [at.dry, sensible, z], [at.latent, latentColor, z], [at.total, totalColor, z]]) { ctx.fillStyle = color; ctx.beginPath(); ctx.arc(px, scale(v), 4, 0, Math.PI * 2); ctx.fill(); }

    const gap = at.evaporation - at.precipitation;
    local.set([
      ['latitude', place(probe)], ['evaporation', mm(at.evaporation)], ['rain', mm(at.precipitation)],
      ['E − P', `${gap >= 0 ? '+' : '−'}${Math.abs(gap).toFixed(2)}\u00a0mm a day, the sea ${gap >= 0 ? 'loses' : 'gains'} water`],
      ['heat carried', carried(at.total)], ['as latent heat', carried(at.latent)], ['as the heat of dry air', carried(at.dry)],
    ]);
  }

  function summarize() {
    const north = run.peak('total', 1), south = run.peak('total', -1), equator = run.at(0);
    global.set([
      ['whole planet: evaporation', mm(year.evaporation)], ['rain', mm(year.precipitation)], ['vapor in the air', `${year.precipitableWater.toFixed(1)}\u00a0kg/m², ${year.residenceDays.toFixed(1)} days of rain`],
      ['most heat carried', `${carried(south.transport)} at ${place(south.latitude)}, ${carried(north.transport)} at ${place(north.latitude)}`],
      ['across the equator', `${carried(equator.total)}: ${carried(equator.latent)} as latent heat, ${carried(equator.dry)} as the heat of dry air`],
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
    run = waterTransport(frames.header);
    year = recordedYear(frames.header);
    summarize(); fig.render();
  }).catch((error) => { failed = true; console.info('watercycle figure', error.message); fig.render(); });
}
