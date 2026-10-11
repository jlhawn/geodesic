import { Figure, slider, choice, buttons, readout, legend, caption, text, termColor, ACCENT, MUTED, GRID, LINE } from '../runtime.module.js';
import { loadLatitudes, createLatitudes, heldSuarezSurface } from '../latitudeColumns.module.js';
import { GREENHOUSE_GASES } from '../../js/physics/radiation.module.js';

const BUDGET = 10, DEFAULT = { carbonDioxide: 1, albedo: 0.07 }, FLUX_TOP = 450, HALO = 'rgba(21, 21, 21, 0.85)';
const TARGET = 'rgba(255, 255, 255, 0.7)', GHOST = 'rgba(255, 255, 255, 0.3)', PROBE = 'rgba(255, 255, 255, 0.28)';
const degrees = (lat) => `${lat.toFixed(lat % 1 ? 1 : 0)}°`;
const coldest = ({ pressure, height, temperature }) => `${(pressure / 100).toFixed(0)} hPa, ${(height / 1000).toFixed(1)} km, ${temperature.toFixed(1)} K`;

function withAlpha(color, alpha) {
  const hex = /^#([0-9a-f]{6})$/i.exec(color);
  return hex ? `rgba(${[0, 2, 4].map((n) => parseInt(hex[1].slice(n, n + 2), 16)).join(', ')}, ${alpha})` : color;
}

export function mountLatitudes(root) {
  const controls = root.querySelector('.controls');
  const sun = termColor('a'), infrared = termColor('b'), faintSun = withAlpha(sun, 0.85), reflectedSun = withAlpha(sun, 0.2);
  const settings = { ...DEFAULT };
  let data = null, band = null, reference = null, summary = null, view = 'temperature', probe = 0, perFrame = 4, failed = false, geometry = null, ghost = false;

  const fig = new Figure(root, { height: (w) => (w < 520 ? 300 : 340), step, draw });
  caption(root, 'Ten columns of the model’s own radiation and dry convection, at the N = 16 grid’s cells nearest 0°, 10°, 20° and so on to the pole, each under its own year-round average of sunlight and ozone. No heat passes between them. Change a setting and they run again, live, until every column is in balance. Tap or hover over the chart to read one latitude.');
  const key = legend(root, legendFor());
  choice(controls, { label: 'Show', options: [['ground temperature', 'temperature'], ['sunlight and infrared', 'flux']], value: view, onChange: (v) => { view = v; key.set(legendFor()); fig.render(); }, span: true });
  const carbon = slider(controls, { label: 'Carbon dioxide', min: 0.5, max: 4, step: 0.25, value: DEFAULT.carbonDioxide, format: (v) => `×${v}, ${Math.round(v * GREENHOUSE_GASES.carbonDioxide * 1e6)} ppm`, onInput: (v) => change({ carbonDioxide: v }), span: true });
  const ground = slider(controls, { label: 'Ground reflects', min: 0.07, max: 0.4, step: 0.01, value: DEFAULT.albedo, format: (v) => `${Math.round(v * 100)}% of sunlight`, onInput: (v) => change({ albedo: v }), span: true });
  buttons(controls, [['Start every column at 288 K', () => { if (band) { band.restart(288); run(); } }], ['Reset', reset]]);
  const local = readout(controls), global = readout(controls);

  const isDefault = () => settings.carbonDioxide === DEFAULT.carbonDioxide && Math.abs(settings.albedo - DEFAULT.albedo) < 1e-9;
  const showGhost = () => Boolean(reference && summary) && summary.rows.some((r, n) => Math.abs(r.Ts - reference.rows[n].Ts) > 0.05);

  function legendFor() {
    if (view === 'flux') return [['faint', 'sunlight arriving at the top', faintSun], ['bar', 'sunlight reflected', reflectedSun], ['line', 'sunlight absorbed', sun], ['line', 'infrared leaving the top', infrared], ['dots', 'area mean absorbed', MUTED]];
    return [
      ['line', 'ground temperature, each latitude for itself', ACCENT], ['dash', 'Held–Suarez target at the ground, 315 − 60 sin²φ K', TARGET], ['dots', 'area mean', MUTED],
      ...(ghost ? [['faint', 'in balance with carbon dioxide ×1 and the ground reflecting 7%', GHOST]] : []),
    ];
  }

  function change(changes) {
    Object.assign(settings, changes);
    if (!band) return;
    band.set(changes);
    run();
  }

  function reset() {
    Object.assign(settings, DEFAULT);
    carbon.value = DEFAULT.carbonDioxide;
    ground.value = DEFAULT.albedo;
    if (!band) return;
    band.load(data, DEFAULT);
    run();
  }

  function run() {
    summary = band.summary();
    syncGhost();
    fig.play(true);
    fig.render();
  }

  function step() {
    if (!band || band.settled) { fig.play(false); return false; }
    const start = performance.now();
    band.advance(perFrame);
    const ms = performance.now() - start;
    if (ms < BUDGET * 0.7) perFrame = Math.ceil(perFrame * 1.25);
    else if (ms > BUDGET) perFrame = Math.max(1, Math.floor(perFrame * 0.8));
    summary = band.summary();
    syncGhost();
    if (band.settled) fig.play(false);
  }

  function syncGhost() {
    if (showGhost() === ghost) return;
    ghost = !ghost;
    key.set(legendFor());
  }

  function draw(ctx, w, h) {
    if (!summary) { text(ctx, failed ? 'The columns could not be loaded.' : 'Loading the columns…', w / 2, h / 2, { align: 'center', color: MUTED }); return; }
    const narrow = w < 520, left = narrow ? 42 : 50, right = narrow ? 10 : 16, top = 26, bottom = h - 26, size = narrow ? 10 : 11;
    const x = (lat) => left + lat / 90 * (w - left - right), { rows, mean } = summary, lats = rows.map((r) => r.latitude);
    geometry = { x, lats };
    const flux = view === 'flux';
    const lowest = Math.min(...rows.map((r) => r.Ts), ...(reference ? reference.rows.map((r) => r.Ts) : []));
    const lo = flux ? 0 : Math.min(220, Math.floor((lowest - 5) / 20) * 20), hi = flux ? FLUX_TOP : 320, tick = flux ? 100 : 20;
    const y = (v) => bottom - (v - lo) / (hi - lo) * (bottom - top);
    const line = (values, color, width, dash = []) => {
      ctx.save(); ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash); ctx.lineJoin = 'round';
      ctx.beginPath(); values.forEach((v, n) => (n ? ctx.lineTo(x(lats[n]), y(v)) : ctx.moveTo(x(lats[n]), y(v)))); ctx.stroke(); ctx.restore();
    };

    ctx.lineWidth = 1;
    for (let v = lo; v <= hi; v += tick) {
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(left, y(v)); ctx.lineTo(w - right, y(v)); ctx.stroke();
      text(ctx, `${v}`, left - 6, y(v), { align: 'right', color: MUTED, size });
    }
    text(ctx, flux ? 'W/m²' : 'K', left - 6, top - 14, { align: 'right', color: MUTED, size });
    ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(left, bottom); ctx.lineTo(w - right, bottom); ctx.stroke();
    for (const [lat, label, align] of [[0, 'equator', 'left'], [30, '30°', 'center'], [60, '60°', 'center'], [90, 'pole', 'right']]) text(ctx, label, x(lat), bottom + 13, { align, color: MUTED, size });

    const at = rows[probe], px = x(at.latitude);
    ctx.strokeStyle = PROBE; ctx.beginPath(); ctx.moveTo(px, top); ctx.lineTo(px, bottom); ctx.stroke();

    if (flux) {
      ctx.fillStyle = reflectedSun; ctx.beginPath();
      rows.forEach((r, n) => (n ? ctx.lineTo(x(lats[n]), y(r.insolation)) : ctx.moveTo(x(lats[n]), y(r.insolation))));
      for (let n = rows.length - 1; n >= 0; n--) ctx.lineTo(x(lats[n]), y(rows[n].absorbed));
      ctx.closePath(); ctx.fill();
      line(rows.map((r) => r.insolation), faintSun, 1.5);
      line(rows.map((r) => r.absorbed), sun, 3.5);
      line(rows.map((r) => r.olr), infrared, 1.5);
      areaLine(ctx, mean.absorbed, `area mean ${mean.absorbed.toFixed(0)} W/m²`);
      for (const [v, color] of [[at.absorbed, sun], [at.olr, infrared]]) { ctx.fillStyle = color; ctx.beginPath(); ctx.arc(px, y(v), 4, 0, Math.PI * 2); ctx.fill(); }
    } else {
      ctx.save(); ctx.strokeStyle = TARGET; ctx.lineWidth = 1.5; ctx.setLineDash([4, 4]); ctx.beginPath();
      for (let lat = 0; lat <= 90; lat++) (lat ? ctx.lineTo(x(lat), y(heldSuarezSurface(lat))) : ctx.moveTo(x(lat), y(heldSuarezSurface(lat))));
      ctx.stroke(); ctx.restore();
      if (ghost) line(reference.rows.map((r) => r.Ts), GHOST, 1.5);
      areaLine(ctx, mean.Ts, `area mean ${mean.Ts.toFixed(1)} K`);
      line(rows.map((r) => r.Ts), ACCENT, 2);
      ctx.fillStyle = ACCENT;
      rows.forEach((r, n) => { ctx.beginPath(); ctx.arc(x(lats[n]), y(r.Ts), n === probe ? 5 : 3, 0, Math.PI * 2); ctx.fill(); });
      ctx.strokeStyle = TARGET; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(px, y(at.target), 4, 0, Math.PI * 2); ctx.stroke();
      text(ctx, 'Held–Suarez target', x(narrow ? 6 : 14), y(heldSuarezSurface(narrow ? 6 : 14)) - 11, { color: TARGET, size, halo: HALO });
    }
    text(ctx, band.settled ? 'every column in balance' : `running: day ${Math.floor(band.days)}`, w - right, top - 14, { align: 'right', color: band.settled ? MUTED : ACCENT, size });

    function areaLine(c, v, label) {
      c.save(); c.strokeStyle = MUTED; c.lineWidth = 1.5; c.setLineDash([1.5, 3.5]); c.lineCap = 'round';
      c.beginPath(); c.moveTo(left, y(v)); c.lineTo(w - right, y(v)); c.stroke(); c.restore();
      text(c, label, w - right, y(v) - 9, { align: 'right', color: MUTED, size, halo: HALO });
    }

    const settled = band.settled, every = settled ? 1 : 4, watts = (v) => `${v.toFixed(1)} W/m²`;
    const fluxes = (absorbed, olr) => (settled ? [['sunlight absorbed = infrared out', watts(absorbed)]] : [['sunlight absorbed', watts(absorbed)], ['infrared out', watts(olr)]]);
    local.set([['latitude', degrees(at.latitude)], ['ground', `${at.Ts.toFixed(1)} K`], ['Held–Suarez target', `${at.target.toFixed(1)} K`], ['500 hPa', `${at.T500.toFixed(1)} K`], ['tropopause', coldest(at.coldPoint)], ...fluxes(at.absorbed, at.olr)], every);
    global.set([['area mean: ground', `${mean.Ts.toFixed(1)} K`], ['effective temperature', `${mean.Te.toFixed(1)} K, ${(mean.Ts - mean.Te).toFixed(1)} K colder`], ...fluxes(mean.absorbed, mean.olr), ['equator to pole', `${summary.contrast.toFixed(1)} K, Held–Suarez ${summary.targetContrast.toFixed(0)} K`]], every);
  }

  const pick = (event) => {
    if (!geometry) return;
    const r = fig.canvas.getBoundingClientRect(), px = event.clientX - r.left;
    let best = 0;
    geometry.lats.forEach((lat, n) => { if (Math.abs(geometry.x(lat) - px) < Math.abs(geometry.x(geometry.lats[best]) - px)) best = n; });
    if (best !== probe) { probe = best; fig.render(); }
  };
  fig.canvas.style.touchAction = 'pan-y';
  fig.canvas.addEventListener('pointerdown', pick);
  fig.canvas.addEventListener('pointermove', pick);

  loadLatitudes().then((loaded) => {
    data = loaded;
    band = createLatitudes(data, DEFAULT);
    reference = band.summary();
    if (!isDefault()) band.set(settings);
    run();
  }).catch((error) => { failed = true; console.info('latitudes figure', error.message); fig.render(); });
}
