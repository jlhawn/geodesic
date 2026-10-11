import { Figure, slider, readout, legend, caption, text, thermometer, rampColor, clamp, ACCENT, MUTED, GRID, LINE, VAPOR, CLOUD, P0 } from '../runtime.module.js';
import { heightOf } from '../physics.module.js';
import { liftParcel, parcelAt, moistLapse, DRY_LAPSE } from '../moistParcel.module.js';
import { saturationHumidity } from '../../js/physics/moist.module.js';

const Z_MAX = 12000, P_TOP = 190e2, C0 = 273.15, WATER_MAX = 40;
const CLOUD_SHADE = '#aeb8c6', BASE = 'rgba(255, 255, 255, 0.6)', HALO = 'rgba(20, 30, 60, 0.75)';
const minus = (v, d) => (Math.abs(v) < 0.5 * 10 ** -d ? (0).toFixed(d) : `${v < 0 ? '−' : ''}${Math.abs(v).toFixed(d)}`);
const meters = (z) => `${Math.round(z).toLocaleString('en-US')} m`;
const wobble = (n) => { const s = Math.sin(n * 12.9898 + 4.1414) * 43758.5453; return s - Math.floor(s); };

export function reading(path, height) {
  const at = parcelAt(path, height), limit = saturationHumidity(at.T, at.p);
  return { ...at, limit, rh: Math.min(1, at.q / limit), volume: (at.T / path.T0) / (at.p / P0), lapse: at.cloudy ? moistLapse(at.T, at.p) : DRY_LAPSE };
}

function gauge(ctx, x, top, height, { vapor, cloud, limit, label }) {
  const width = 20, bottom = top + height, scale = (v) => bottom - clamp(v / WATER_MAX, 0, 1) * height;
  ctx.save();
  ctx.lineWidth = 1.5; ctx.strokeStyle = LINE; ctx.fillStyle = 'rgba(255,255,255,0.06)';
  ctx.beginPath(); ctx.roundRect(x - width / 2, top, width, height, 3); ctx.fill(); ctx.stroke();
  const inner = width - 6, yv = scale(vapor), yc = scale(vapor + cloud);
  if (vapor > 0) { ctx.fillStyle = VAPOR; ctx.fillRect(x - inner / 2, yv, inner, bottom - 2 - yv); }
  if (cloud > 0) { ctx.fillStyle = CLOUD; ctx.fillRect(x - inner / 2, yc, inner, yv - yc); }
  ctx.strokeStyle = ACCENT; ctx.lineWidth = 2; ctx.setLineDash([4, 3]);
  ctx.beginPath(); ctx.moveTo(x - width / 2 - 5, scale(limit)); ctx.lineTo(x + width / 2 + 5, scale(limit)); ctx.stroke();
  ctx.setLineDash([]);
  for (let v = 10; v <= WATER_MAX; v += 10) { const y = scale(v); ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(x + width / 2 + 2, y); ctx.lineTo(x + width / 2 + 6, y); ctx.stroke(); text(ctx, `${v}`, x + width / 2 + 9, y, { color: MUTED, size: 10 }); }
  ctx.restore();
  text(ctx, `${vapor.toFixed(1)} g/kg`, x, top - 12, { align: 'center', color: '#fff', weight: 500 });
  label.forEach((line, n) => text(ctx, line, x, bottom + 7.2 + n * 12, { align: 'center', color: MUTED, size: 11 }));
}

function tower(ctx, cx, yBase, yTop, radius, left, width, top) {
  if (yTop >= yBase) return;
  ctx.save();
  ctx.beginPath(); ctx.rect(left, top - radius * 2, width, yBase - top + radius * 2); ctx.clip();
  const fill = ctx.createLinearGradient(0, yTop - radius, 0, yBase);
  fill.addColorStop(0, CLOUD); fill.addColorStop(1, CLOUD_SHADE);
  ctx.fillStyle = fill; ctx.globalAlpha = 0.88;
  ctx.beginPath();
  const puff = (x, y, r) => { ctx.moveTo(x + r, y); ctx.arc(x, y, r, 0, Math.PI * 2); };
  for (let n = 0, y = yBase; y > yTop + radius * 0.4; n++, y -= radius * 0.62) {
    const side = n % 2 ? 1 : -1, r = radius * (0.8 + 0.25 * wobble(n));
    puff(cx + radius * 0.25 * (wobble(n + 7) - 0.5), y, r);
    puff(cx + side * radius * (0.55 + 0.2 * wobble(n + 3)), y + radius * 0.15, radius * (0.55 + 0.15 * wobble(n + 11)));
  }
  puff(cx, yTop, radius);
  puff(cx - radius * 0.75, yTop + radius * 0.4, radius * 0.72);
  puff(cx + radius * 0.8, yTop + radius * 0.35, radius * 0.68);
  ctx.fill();
  ctx.restore();
}

export function mountMoistlift(root) {
  const controls = root.querySelector('.controls');
  const state = { z: 0, rh: 0.7, ground: 15 };
  let path = null;
  const lift = () => { path = liftParcel({ T0: C0 + state.ground, rh: state.rh }, P_TOP); };
  const parcel = () => reading(path, state.z);

  const layout = (w) => {
    const narrow = w < 480, left = narrow ? 36 : 44, labels = narrow ? 44 : 56;
    let panel = Math.min(240, Math.round(w * 0.38)), columns = 4;
    if ((w - left - panel - labels - 8) / 4 < 70) { columns = 2; panel = clamp(w - left - labels - 8 - 2 * 76, 90, 240); }
    const start = left + panel + labels, column = (w - start - 4) / columns;
    return { narrow, left, panel, right: left + panel, start, column, columns };
  };

  const fig = new Figure(root, { height: (w) => (layout(w).columns === 4 ? 440 : 520), draw });
  caption(root, 'One parcel of air lifted from the ground at 1000 hPa, with the heights of Part 1’s standard atmosphere. Below its cloud base it cools as a dry parcel does. Above it the parcel holds only as much vapor as its temperature allows, and the rest condenses into cloud, here dropped behind the parcel at once. The gauge on the right holds the parcel’s water: the vapor it still carries and the cloud water it has made.');
  legend(root, [['ramp', 'parcel colder to warmer', 'cool', 'warm', 'neutral'], ['bar', 'water vapor', VAPOR], ['bar', 'cloud water made', CLOUD], ['dash', 'the most vapor it can hold', ACCENT], ['dash', 'cloud base', BASE]]);
  slider(controls, { label: 'Lift the parcel to', min: 0, max: Z_MAX, step: 50, value: state.z, format: (v) => `${(v / 1000).toFixed(2)} km`, onInput: (v) => { state.z = v; update(); } });
  slider(controls, { label: 'Humidity at the ground', min: 0, max: 100, step: 1, value: state.rh * 100, format: (v) => `${v}%`, onInput: (v) => { state.rh = v / 100; lift(); update(); } });
  slider(controls, { label: 'Temperature at the ground', min: 0, max: 35, step: 0.5, value: state.ground, format: (v) => `${v.toFixed(1)} °C`, onInput: (v) => { state.ground = v; lift(); update(); } });
  const out = readout(controls);

  function update() {
    const s = parcel();
    out.set([
      ['pressure there', `${(s.p / 100).toFixed(0)} hPa`],
      ['parcel temperature', `${minus(s.T - C0, 1)} °C`],
      ['dew point', Number.isFinite(s.dew) ? `${minus(s.dew - C0, 1)} °C` : 'none, the air is dry'],
      ['vapor', `${(1000 * s.q).toFixed(2)} g/kg of the ${(1000 * s.limit).toFixed(2)} it can hold`],
      ['relative humidity', `${Math.round(100 * s.rh)}%`],
      ['cloud water made', `${(1000 * s.condensate).toFixed(2)} g/kg`],
      ['warmer than a dry parcel by', `${minus(s.T - s.dry, 1)} K`],
      ['θ', `${s.theta.toFixed(1)} K`],
    ]);
    fig.render();
  }

  function draw(ctx, w, h) {
    const { narrow, left, panel, right, start, column, columns } = layout(w), top = 30, bottom = h - 34;
    const y = (z) => bottom - (z / Z_MAX) * (bottom - top), s = parcel(), base = path.lcl && path.lcl.height <= Z_MAX ? path.lcl.height : null;
    const sky = ctx.createLinearGradient(0, top, 0, bottom);
    sky.addColorStop(0, 'rgba(20, 30, 70, 0.9)'); sky.addColorStop(1, 'rgba(70, 120, 190, 0.6)');
    ctx.fillStyle = sky; ctx.fillRect(left, top, panel, bottom - top);
    ctx.fillStyle = '#5a4a36'; ctx.fillRect(left, bottom, panel, 6);
    for (let km = 0; km <= 12; km += 2) { ctx.strokeStyle = GRID; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(left, y(km * 1000)); ctx.lineTo(right, y(km * 1000)); ctx.stroke(); text(ctx, `${km} km`, left - 6, y(km * 1000), { align: 'right', color: MUTED, size: narrow ? 10 : 11 }); }
    for (const hPa of narrow ? [700, 500, 300, 200] : [1000, 850, 700, 500, 400, 300, 200]) text(ctx, `${hPa} hPa`, right + 5, y(heightOf(hPa * 100)), { color: MUTED, size: 10 });

    const r = 13 * Math.sqrt(s.volume), cx = left + panel / 2, cy = y(s.z) - r * (1 - s.z / Z_MAX);
    if (base !== null && s.z > base) tower(ctx, cx, y(base), cy, Math.min(panel * 0.2, 24), left, panel, top);
    const trail = base !== null ? Math.max(y(base), cy + r) : cy + r;
    ctx.strokeStyle = 'rgba(255,232,160,0.5)'; ctx.lineWidth = 1; ctx.setLineDash([2, 5]); ctx.beginPath(); ctx.moveTo(cx, bottom); ctx.lineTo(cx, Math.min(bottom, trail)); ctx.stroke(); ctx.setLineDash([]);
    if (base !== null) {
      const yb = y(base);
      ctx.strokeStyle = BASE; ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]); ctx.beginPath(); ctx.moveTo(left, yb); ctx.lineTo(right, yb); ctx.stroke(); ctx.setLineDash([]);
      text(ctx, narrow ? meters(base) : `cloud base, ${meters(base)}`, right - 4, Math.max(top + 8, yb - 8), { align: 'right', color: '#fff', size: 10, halo: HALO });
    }
    ctx.fillStyle = rampColor((s.T - 213) / 90); ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.arc(cx, cy, r, 0, Math.PI * 2); ctx.fill(); ctx.stroke();

    const cell = (i) => { const row = columns === 4 ? 0 : Math.floor(i / 2), col = columns === 4 ? i : i % 2; return { x: start + column * (col + 0.5) - 8, row }; };
    const rowTop = (row) => (columns === 4 ? 46 : 42 + row * (h - 20) / 2), tube = columns === 4 ? h - 130 : (h - 20) / 2 - 88;
    const lapse = `${minus(-1000 * s.lapse, 1)} °C/km`, dryLapse = `${minus(-1000 * DRY_LAPSE, 1)} °C/km`;
    const scale = { min: -120, max: 40 };
    let at = cell(0); thermometer(ctx, at.x, rowTop(at.row), tube, s.T - C0, { ...scale, label: ['parcel', 'now', lapse] });
    at = cell(1); thermometer(ctx, at.x, rowTop(at.row), tube, s.dry - C0, { ...scale, label: ['a dry', 'parcel', dryLapse] });
    at = cell(2); thermometer(ctx, at.x, rowTop(at.row), tube, s.theta, { min: 260, max: 420, unit: 'K', label: ['potential', 'temperature θ'] });
    at = cell(3); gauge(ctx, at.x, rowTop(at.row), tube + 14, { vapor: 1000 * s.q, cloud: 1000 * s.condensate, limit: 1000 * s.limit, label: ['vapor and', 'cloud water'] });
  }

  lift();
  update();
}
