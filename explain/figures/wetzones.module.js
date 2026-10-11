import { Figure, choice, readout, legend, caption, text, rampRGB, clamp, paletteVersion, DARK_NEUTRAL, MUTED, GRID, LINE, VAPOR } from '../runtime.module.js';
import { loadFrames } from '../frames.module.js';
import { DATA_URL, overturning, peakLatitude } from '../wetRun.module.js';

const P_BOTTOM = 1000, P_TOP = 30, D = Math.PI / 180, RADIUS = 6.371e6;
const DARK = 'rgb(52, 55, 62)', PALE = [212, 239, 255], PROBE = 'rgba(255, 255, 255, 0.45)';
const MONTHS = [['the year', ''], ['March', 'March'], ['June', 'June'], ['September', 'September'], ['December', 'December']];
const place = (lat) => (Math.abs(lat) < 0.25 ? '0°' : `${Math.abs(lat).toFixed(Math.abs(lat) % 1 ? 1 : 0)}°\u00a0${lat > 0 ? 'N' : 'S'}`);
const parse = (color) => (Array.isArray(color) ? color : color.startsWith('#') ? color.slice(1).match(/../g).map((h) => parseInt(h, 16)) : color.match(/\d+/g).slice(0, 3).map(Number));

function blend(colors) {
  const list = colors.map(parse);
  return (t, rgb) => {
    const x = clamp(t, 0, 1) * (list.length - 1), n = Math.min(list.length - 2, Math.floor(x)), u = x - n;
    for (let c = 0; c < 3; c++) rgb[c] = (list[n][c] + u * (list[n + 1][c] - list[n][c])) / 255;
    return rgb;
  };
}

const watery = blend([DARK_NEUTRAL, VAPOR, PALE]);
const VIEWS = [
  { key: 'RH', label: 'relative humidity', color: (v, rgb) => watery(v, rgb), legend: ['ramp', 'relative humidity over liquid water, from 0 to 100 percent', DARK, `rgb(${PALE})`, VAPOR] },
  { key: 'q', label: 'water vapor', color: (v, rgb) => watery((Math.log10(Math.max(v, 1e-3)) + 2) / 3.3, rgb), legend: ['ramp', 'water vapor, from 0.01 to 20\u00a0g/kg on a logarithmic scale', DARK, `rgb(${PALE})`, VAPOR] },
  { key: 'T', label: 'temperature', color: (v, rgb) => rampRGB((v - 190) / 110, rgb, DARK_NEUTRAL), legend: ['ramp', 'temperature, from 190 to 300\u00a0K', 'cool', 'warm', DARK] },
  { key: 'u', label: 'eastward wind', color: (v, rgb) => rampRGB(0.5 + v / 120, rgb, DARK_NEUTRAL), legend: ['ramp', 'wind from the east in cool colors, from the west in warm ones, up to 60\u00a0m/s', 'cool', 'warm', DARK] },
];
const MOTION = ['arrow', 'the air’s motion north or south and up or down, averaged around each circle of latitude, fainter where it is slower', 'rgba(255, 255, 255, 0.9)'];

export function mountWetzones(root) {
  const controls = root.querySelector('.controls');
  let header = null, failed = false, view = 0, month = '', probe = { lat: -22.5, p: 500 }, geometry = null, image = null, shown = null, flow = null;

  const fig = new Figure(root, { height: (w) => (w < 520 ? 340 : 390), draw });
  caption(root, 'The air of the run above, averaged around each circle of latitude over the year or over one month, from the sea up to 30\u00a0hPa. Pressure falls up the page on a logarithmic scale, so that each step up thins the air by the same factor. Tap or hover over the section to read one place.');
  const key = legend(root, [VIEWS[view].legend, MOTION]);
  choice(controls, { label: 'Show', options: VIEWS.map((v, n) => [v.label, String(n)]), value: '0', onChange: (v) => { view = Number(v); key.set([VIEWS[view].legend, MOTION]); fig.render(); }, span: true });
  choice(controls, { label: 'Average over', options: MONTHS, value: month, onChange: (v) => { month = v; flow = null; summarize(); fig.render(); }, span: true });
  const local = readout(controls), global = readout(controls);

  const zonal = () => (month ? header.zonalMean.months.find((m) => m.name === month) : header.zonalMean);
  const budget = () => (month ? header.zonalBudget.months.find((m) => m.name === month) : header.zonalBudget);

  function sampler(sigma) {
    const K = sigma.length, logs = sigma.map(Math.log);
    const across = (lat) => { const x = clamp((lat + 90) / 5 - 0.5, 0, 35), b = Math.min(34, Math.floor(x)); return [b, x - b]; };
    const down = (p) => { const s = Math.log(p / 1000); if (s <= logs[0]) return [0, 0]; if (s >= logs[K - 1]) return [K - 2, 1]; let k = 0; while (logs[k + 1] < s) k++; return [k, (s - logs[k]) / (logs[k + 1] - logs[k])]; };
    const at = (values, lat, p, stride = K, levels = null) => {
      const [b, tb] = across(lat);
      let k, tk;
      if (levels) { const s = p / 1000; k = 0; while (k + 2 < levels.length && levels[k + 1] < s) k++; tk = clamp((s - levels[k]) / (levels[k + 1] - levels[k]), 0, 1); } else [k, tk] = down(p);
      const v = (bb, kk) => values[bb * stride + kk];
      return (1 - tb) * ((1 - tk) * v(b, k) + tk * v(b, k + 1)) + tb * ((1 - tk) * v(b + 1, k) + tk * v(b + 1, k + 1));
    };
    return { across, down, at };
  }

  function paint(w, h, box) {
    const id = `${view}/${month}/${w}x${h}/${paletteVersion}`;
    if (shown === id) return image;
    const zm = zonal(), values = zm[VIEWS[view].key], { at } = sampler(header.zonalMean.sigma), color = VIEWS[view].color;
    const width = Math.round(box.right - box.left), height = Math.round(box.bottom - box.top);
    image ??= document.createElement('canvas');
    image.width = width; image.height = height;
    const ctx = image.getContext('2d'), pixels = ctx.createImageData(width, height), rgb = [0, 0, 0];
    for (let y = 0; y < height; y++) {
      const p = box.pressure(box.top + y + 0.5);
      for (let x = 0; x < width; x++) {
        color(at(values, box.latitude(box.left + x + 0.5), p), rgb);
        const n = 4 * (y * width + x);
        pixels.data[n] = rgb[0] * 255; pixels.data[n + 1] = rgb[1] * 255; pixels.data[n + 2] = rgb[2] * 255; pixels.data[n + 3] = 255;
      }
    }
    ctx.putImageData(pixels, 0, 0);
    shown = id;
    return image;
  }

  function draw(ctx, w, h) {
    if (!header) { text(ctx, failed ? 'The run could not be loaded.' : 'Loading the run…', w / 2, h / 2, { align: 'center', color: MUTED }); return; }
    const narrow = w < 520, size = narrow ? 10 : 11, left = narrow ? 38 : 46, right = narrow ? 8 : 14, top = 22, bottom = h - 24;
    const span = Math.log(P_BOTTOM / P_TOP);
    const box = {
      left, right: w - right, top, bottom,
      x: (lat) => left + (lat + 90) / 180 * (w - right - left),
      y: (p) => bottom - Math.log(P_BOTTOM / p) / span * (bottom - top),
      latitude: (x) => (x - left) / (w - right - left) * 180 - 90,
      pressure: (y) => P_BOTTOM * Math.exp(-(bottom - y) / (bottom - top) * span),
    };
    geometry = box;
    ctx.drawImage(paint(w, h, box), left, top, w - right - left, bottom - top);

    ctx.lineWidth = 1;
    const pressures = narrow ? [1000, 500, 200, 100, 50] : [1000, 700, 500, 300, 200, 100, 50];
    for (const p of pressures) {
      ctx.strokeStyle = p === 1000 ? LINE : GRID; ctx.beginPath(); ctx.moveTo(left, box.y(p)); ctx.lineTo(w - right, box.y(p)); ctx.stroke();
      text(ctx, `${p}`, left - 5, box.y(p), { align: 'right', color: MUTED, size });
    }
    text(ctx, 'hPa', left - 5, top - 12, { align: 'right', color: MUTED, size });
    for (const lat of [-60, -30, 0, 30, 60]) {
      ctx.strokeStyle = GRID; ctx.beginPath(); ctx.moveTo(box.x(lat), top); ctx.lineTo(box.x(lat), bottom); ctx.stroke();
      text(ctx, lat === 0 ? '0°' : `${Math.abs(lat)}°${lat > 0 ? 'N' : 'S'}`, box.x(lat), bottom + 13, { align: 'center', color: MUTED, size });
    }

    flow ??= overturning(header, month || null);
    const { at } = sampler(header.zonalMean.sigma), K = flow.K, length = narrow ? 12 : 10;
    const levels = narrow ? Array.from({ length: 7 }, (_, n) => 920 * 0.7 ** n) : Array.from({ length: 9 }, (_, n) => 920 * 0.77 ** n), lats = [];
    for (let lat = -80; lat <= 80; lat += narrow ? 20 : 10) lats.push(lat);
    const pxPerSecond = (lat, p, v, omega) => [v / (RADIUS * D) / 180 * (w - right - left), omega / (p * 100) / span * (bottom - top)];
    for (const lat of lats) for (const p of levels) {
      const v = at(flow.v, lat, p), omega = at(flow.omega, lat, p, K + 1, flow.levels);
      const [dx, dy] = pxPerSecond(lat, p, v, omega), norm = Math.hypot(dx, dy);
      if (!norm) continue;
      const strength = clamp(Math.hypot(v / 1, omega / 0.02), 0, 1);
      if (strength < 0.08) continue;
      const ux = dx / norm, uy = dy / norm, cx = box.x(lat), cy = box.y(p), half = length * (0.55 + 0.45 * strength) / 2;
      const x0 = cx - ux * half, y0 = cy - uy * half, x1 = cx + ux * half, y1 = cy + uy * half;
      const head = () => { ctx.beginPath(); ctx.moveTo(x1 + ux * 2, y1 + uy * 2); ctx.lineTo(x1 - ux * 4.5 - uy * 3.5, y1 - uy * 4.5 + ux * 3.5); ctx.lineTo(x1 - ux * 4.5 + uy * 3.5, y1 - uy * 4.5 - ux * 3.5); ctx.closePath(); };
      for (const [color, width, outline] of [[`rgba(21, 21, 21, ${(0.35 + 0.5 * strength).toFixed(2)})`, 3.6, 2.2], [`rgba(255, 255, 255, ${(0.3 + 0.7 * strength).toFixed(2)})`, 1.6, 0]]) {
        ctx.save(); ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = width; ctx.lineCap = 'round'; ctx.lineJoin = 'round';
        ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
        head(); ctx.fill();
        if (outline) { ctx.lineWidth = outline; ctx.stroke(); }
        ctx.restore();
      }
    }

    const px = box.x(probe.lat), py = box.y(probe.p);
    ctx.strokeStyle = PROBE; ctx.beginPath(); ctx.moveTo(px, top); ctx.lineTo(px, bottom); ctx.moveTo(left, py); ctx.lineTo(w - right, py); ctx.stroke();
    ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(px, py, 3.5, 0, Math.PI * 2); ctx.fill();

    const zm = zonal(), get = (name) => at(zm[name], probe.lat, probe.p), [b, tb] = sampler(header.zonalMean.sigma).across(probe.lat);
    const ps = (1 - tb) * zm.surfacePressure[b] + tb * zm.surfacePressure[b + 1];
    const v = at(flow.v, probe.lat, probe.p), omega = at(flow.omega, probe.lat, probe.p, K + 1, flow.levels), perDay = omega * 864;
    const wind = get('u');
    local.set([
      ['latitude', place(probe.lat)], ['pressure', `${(probe.p / 1000 * ps).toFixed(0)}\u00a0hPa`], ['temperature', `${get('T').toFixed(1)}\u00a0K`],
      ['relative humidity', `${Math.round(100 * get('RH'))}%`], ['water vapor', `${get('q').toPrecision(2)}\u00a0g/kg`], ['cloud water', `${Math.max(0, get('qc')).toPrecision(2)}\u00a0mg/kg`],
      ['wind', `${Math.abs(wind).toFixed(1)}\u00a0m/s from the ${wind >= 0 ? 'west' : 'east'}`],
      ['motion', `${Math.abs(v).toFixed(2)}\u00a0m/s ${v >= 0 ? 'north' : 'south'}, ${perDay < 0 ? 'rising' : 'sinking'} ${Math.abs(perDay).toFixed(0)}\u00a0hPa a day`],
    ]);
  }

  function summarize() {
    const zm = zonal(), K = header.K, sigma = header.zonalMean.sigma, rain = peakLatitude(budget().precipitation);
    const k500 = sigma.findIndex((s) => s > 0.5);
    let driest = { rh: Infinity };
    for (let b = 10; b < 26; b++) if (zm.RH[b * K + k500] < driest.rh) driest = { rh: zm.RH[b * K + k500], lat: -90 + (b + 0.5) * 5 };
    let jet = { u: -Infinity };
    for (let b = 0; b < 36; b++) for (let k = 0; k < K; k++) if (sigma[k] > 0.1 && zm.u[b * K + k] > jet.u) jet = { u: zm.u[b * K + k], lat: -90 + (b + 0.5) * 5, p: sigma[k] * zm.surfacePressure[b] };
    global.set([
      ['rain band', `${place(rain.latitude)}, ${rain.value.toFixed(1)}\u00a0mm a day`],
      ['driest air near 500\u00a0hPa', `${Math.round(100 * driest.rh)}% at ${place(driest.lat)}`],
      ['strongest wind below 100\u00a0hPa', `${jet.u.toFixed(0)}\u00a0m/s at ${place(jet.lat)}, ${jet.p.toFixed(0)}\u00a0hPa`],
    ]);
  }

  const pick = (event) => {
    if (!geometry) return;
    const r = fig.canvas.getBoundingClientRect(), x = clamp(event.clientX - r.left, geometry.left, geometry.right), y = clamp(event.clientY - r.top, geometry.top, geometry.bottom);
    probe = { lat: clamp(Math.round(2 * geometry.latitude(x)) / 2, -87.5, 87.5), p: geometry.pressure(y) };
    fig.render();
  };
  fig.canvas.style.touchAction = 'pan-y';
  fig.canvas.addEventListener('pointerdown', pick);
  fig.canvas.addEventListener('pointermove', (event) => { if (event.pointerType === 'mouse' || event.buttons) pick(event); });

  loadFrames(DATA_URL).then((frames) => { header = frames.header; summarize(); fig.render(); }).catch((error) => { failed = true; console.info('wetzones figure', error.message); fig.render(); });
}
