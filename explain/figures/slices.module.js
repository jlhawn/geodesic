import { Figure, slider, choice, legend, caption, readout, text, anomalyColor, INK, MUTED, LINE } from '../runtime.module.js';
import { P0, pressureAt, heightOf } from '../physics.module.js';
import { sigmaInterfaces } from '../../js/dynamics/sigmaCore.module.js';

const KM = 1000, HPA = 100, SPAN = 600 * KM, HALF_WIDTH = 120 * KM, Z_MAX = 22 * KM, GROUND = '#5a4a36', SHADE = 'rgba(255, 255, 255, 0.06)', HALO = 'rgba(20, 20, 22, 0.85)', PANEL = 'rgba(20, 20, 22, 0.94)';
const MODEL = sigmaInterfaces('bl36'), TENTHS = Array.from({ length: 11 }, (_, k) => k / 10);
const HEIGHTS = Array.from({ length: 21 }, (_, k) => (k + 1) * KM), PRESSURES = Array.from({ length: 9 }, (_, k) => (900 - 100 * k) * HPA);
const LABELED = [1000, 800, 600, 400, 200, 100];
const INSETS = {
  height: { span: KM, title: 'the lowest kilometer over the plain', marks: [[0, 'ground'], [500, '500 m'], [KM, '1 km']] },
  pressure: { span: 100 * HPA, title: 'the lowest 100 hPa over the plain', marks: [[0, 'ground'], [50 * HPA, '950 hPa'], [100 * HPA, '900 hPa']] },
};

export function mountSlices(root) {
  const controls = root.querySelector('.controls');
  let peak = 3200, mode = 'height', axis = 'height';
  const ground = (x) => (Math.abs(x) < HALF_WIDTH ? peak * Math.cos(Math.PI * x / (2 * HALF_WIDTH)) ** 2 : 0);
  const buried = (z) => (z < peak ? (2 * HALF_WIDTH / Math.PI) * Math.acos(Math.sqrt(z / peak)) : 0);
  const draped = (sigma, base) => (sigma > 0 ? heightOf(sigma * pressureAt(base)) : Infinity);
  const lowestThickness = (levels, base) => draped(levels[levels.length - 2], base) - base;
  const pz = (z) => (Number.isFinite(z) ? pressureAt(z) : 0);

  const fig = new Figure(root, { height: 400, minHeight: 300, draw });
  const note = caption(root, '');
  const key = legend(root, []);
  slider(controls, { label: 'Mountain height', min: 0, max: 6000, step: 100, value: peak, format: (v) => `${(v / KM).toFixed(1)} km`, onInput: (v) => { peak = v; update(); } });
  choice(controls, { label: 'Levels', options: [['height levels', 'height'], ['pressure levels', 'pressure'], ['σ levels', 'sigma'], ['the model’s 36 layers', 'model']], value: mode, onChange: (v) => { mode = v; update(); }, span: true });
  choice(controls, { label: 'Up the page', options: [['height', 'height'], ['pressure', 'pressure']], value: axis, onChange: (v) => { axis = v; update(); }, span: true });
  const out = readout(controls);

  function update() {
    const surface = pressureAt(peak), summit = `${Math.round(surface / HPA)} hPa`;
    const earth = ['ramp', 'the ground', GROUND, GROUND, GROUND];
    if (mode === 'height' || mode === 'pressure') {
      const height = mode === 'height', cut = height ? HEIGHTS.filter((z) => z < peak).length : PRESSURES.filter((p) => p > surface).length;
      const range = height ? (cut === 1 ? '1 km' : `1 to ${cut} km`) : (cut === 1 ? '900 hPa' : `900 to ${1000 - 100 * cut} hPa`);
      key.set([['faint', height ? 'levels of equal height, every kilometer' : 'levels of equal pressure, every 100 hPa', LINE], ['dash', 'where a level runs inside the mountain', 'warm'], earth]);
      out.set([['levels the mountain cuts', cut ? `${cut} of ${height ? HEIGHTS.length : PRESSURES.length}: ${range}` : 'none'], ...(peak ? [['pressure on the summit', summit]] : [])]);
    } else if (mode === 'sigma') {
      key.set([['faint', 'levels of equal σ, every 0.1, with alternate layers shaded', LINE], earth]);
      out.set([['each layer', 'a tenth of its column'], ['the lowest layer over the plain', `100 hPa of air, ${Math.round(lowestThickness(TENTHS, 0))} m thick`], ...(peak ? [['over the summit', `${Math.round(surface / 1000)} hPa of air, ${Math.round(lowestThickness(TENTHS, peak))} m thick`]] : [])]);
    } else {
      const lowest = 1 - MODEL[MODEL.length - 2], inside = MODEL.slice(1, -1), lowKm = (base) => inside.filter((s) => draped(s, base) - base < KM).length, hidden = axis === 'height' ? inside.filter((s) => draped(s, 0) > Z_MAX).length : 0;
      key.set([['faint', 'the boundaries of the model’s 36 layers, with alternate layers shaded', LINE], earth]);
      out.set([['the lowest layer', `${(lowest * 100).toFixed(2)}% of its column, ${Math.round(lowestThickness(MODEL, 0))} m thick over the plain${peak ? `, ${Math.round(lowestThickness(MODEL, peak))} m on the summit` : ''}`], ['layers in the lowest kilometer', `${lowKm(0)} over the plain${peak ? `, ${lowKm(peak)} above the summit` : ''}`], ...(hidden ? [['layers above the top of the frame', `${hidden}`]] : [])]);
    }
    fig.render();
  }

  function draw(ctx, w, h) {
    const left = 38, right = 54, top = 12, bottom = h - 24, plotW = w - left - right, H = bottom - top;
    const x = (m) => left + (m / SPAN + 0.5) * plotW, y = (z) => Math.max(top - 2, axis === 'height' ? bottom - z / Z_MAX * H : bottom - (1 - pz(z) / P0) * H);
    note.textContent = `The standard atmosphere, with 1000 hPa at sea level, over a mountain ${2 * HALF_WIDTH / KM} km across. ${axis === 'height' ? `Heights are drawn ${Math.round(H / Z_MAX / (plotW / SPAN))} times steeper than distances.` : 'Pressure runs up the page, so the top of the frame is zero pressure, the top of the atmosphere.'}`;
    const samples = Math.ceil(plotW / 2), xs = Array.from({ length: samples + 1 }, (_, s) => (s / samples - 0.5) * SPAN), bases = xs.map(ground);
    const curve = (zs, reverse = false) => { for (let s = 0; s <= samples; s++) { const i = reverse ? samples - s : s; ctx.lineTo(x(xs[i]), y(zs[i])); } };

    ctx.save();
    ctx.beginPath(); ctx.rect(left, top, plotW, h - top); ctx.clip();
    const levels = mode === 'sigma' ? TENTHS : mode === 'model' ? MODEL : null;
    if (levels) {
      const heights = Array.from(levels, (sigma) => bases.map((base) => draped(sigma, base)));
      ctx.fillStyle = SHADE;
      for (let k = levels.length - 2; k >= 0; k -= 2) { ctx.beginPath(); curve(heights[k]); curve(heights[k + 1], true); ctx.fill(); }
      ctx.strokeStyle = mode === 'model' ? 'rgba(255, 255, 255, 0.3)' : LINE; ctx.lineWidth = 1;
      for (let k = 1; k < levels.length - 1; k++) { ctx.beginPath(); curve(heights[k]); ctx.stroke(); }
      if (axis === 'pressure') { ctx.strokeStyle = LINE; ctx.beginPath(); ctx.moveTo(left, top + 0.5); ctx.lineTo(w - right, top + 0.5); ctx.stroke(); if (mode === 'sigma') text(ctx, 'σ = 0, the top of the atmosphere', w - right - 5, top + 10, { align: 'right', color: MUTED, size: 10, halo: HALO }); }
    }
    ctx.fillStyle = GROUND;
    ctx.beginPath(); ctx.moveTo(left, h); curve(bases); ctx.lineTo(w - right, h); ctx.closePath(); ctx.fill();
    if (!levels) {
      const flat = mode === 'height' ? HEIGHTS : PRESSURES.map((p) => heightOf(p));
      for (const z of flat) {
        const reach = buried(z), yy = y(z);
        ctx.strokeStyle = LINE; ctx.lineWidth = 1;
        ctx.beginPath(); ctx.moveTo(left, yy); ctx.lineTo(x(-reach), yy); ctx.moveTo(x(reach), yy); ctx.lineTo(w - right, yy); ctx.stroke();
        if (reach) { ctx.strokeStyle = anomalyColor(1, 0.85); ctx.lineWidth = 2; ctx.setLineDash([4, 3]); ctx.beginPath(); ctx.moveTo(x(-reach), yy); ctx.lineTo(x(reach), yy); ctx.stroke(); ctx.setLineDash([]); }
      }
    }
    if (mode === 'sigma') for (const sigma of TENTHS.slice(1, -1)) text(ctx, `σ = ${sigma.toFixed(1)}`, left + 5, y(heightOf(sigma * P0)) - 7, { color: MUTED, size: 10, halo: HALO });
    if (mode === 'model') inset(ctx, left + 12, top + 24, Math.round(Math.min(150, Math.max(96, plotW * 0.26))), Math.round(Math.min(170, Math.max(116, H * 0.45))));
    const bar = 100 * KM / SPAN * plotW;
    ctx.strokeStyle = INK; ctx.lineWidth = 1.2;
    ctx.beginPath(); ctx.moveTo(left + 8, h - 8); ctx.lineTo(left + 8 + bar, h - 8); ctx.moveTo(left + 8, h - 11); ctx.lineTo(left + 8, h - 5); ctx.moveTo(left + 8 + bar, h - 11); ctx.lineTo(left + 8 + bar, h - 5); ctx.stroke();
    text(ctx, '100 km', left + 14 + bar, h - 8, { color: INK, size: 10 });
    ctx.restore();

    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    let last = Infinity;
    for (let z = 0; z <= 30 * KM; z += KM) {
      const yy = y(z), major = z % (2 * KM) === 0;
      if (yy < top) break;
      ctx.beginPath(); ctx.moveTo(left - (major ? 5 : 3), yy); ctx.lineTo(left, yy); ctx.stroke();
      if (major && last - yy >= 13 && (axis === 'height' || z <= 16 * KM)) { text(ctx, `${z / KM} km`, left - 7, yy, { align: 'right', color: MUTED, size: 10 }); last = yy; }
    }
    for (let p = 1000; p >= 100; p -= 100) { const yy = y(heightOf(p * HPA)), major = LABELED.includes(p); ctx.beginPath(); ctx.moveTo(w - right, yy); ctx.lineTo(w - right + (major ? 5 : 3), yy); ctx.stroke(); if (major) text(ctx, `${p} hPa`, w - right + 7, yy, { color: MUTED, size: 10 }); }
    if (axis === 'pressure') { ctx.beginPath(); ctx.moveTo(w - right, top); ctx.lineTo(w - right + 5, top); ctx.stroke(); text(ctx, '0 hPa', w - right + 7, top, { color: MUTED, size: 10 }); }
  }

  function inset(ctx, x0, y0, iw, ih) {
    const { span, title, marks } = INSETS[axis], pad = 8, iy = (d) => y0 + ih - pad - d / span * (ih - 2 * pad);
    const tops = Array.from(MODEL, (sigma) => (axis === 'height' ? draped(sigma, 0) : (1 - sigma) * P0));
    ctx.fillStyle = PANEL; ctx.fillRect(x0, y0, iw, ih);
    ctx.fillStyle = SHADE;
    for (let k = MODEL.length - 2; k >= 0 && tops[k + 1] < span; k -= 2) { const zt = Math.min(tops[k], span); ctx.fillRect(x0 + 1, iy(zt), iw - 2, iy(tops[k + 1]) - iy(zt)); }
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.3)'; ctx.lineWidth = 1;
    for (let k = 1; k < MODEL.length - 1; k++) if (tops[k] < span) { ctx.beginPath(); ctx.moveTo(x0 + 1, iy(tops[k])); ctx.lineTo(x0 + iw - 1, iy(tops[k])); ctx.stroke(); }
    ctx.fillStyle = GROUND; ctx.fillRect(x0 + 1, iy(0), iw - 2, pad - 1);
    ctx.strokeStyle = LINE; ctx.strokeRect(x0 + 0.5, y0 + 0.5, iw - 1, ih - 1);
    for (const [d, label] of marks) { ctx.beginPath(); ctx.moveTo(x0 + iw, iy(d)); ctx.lineTo(x0 + iw + 4, iy(d)); ctx.stroke(); text(ctx, label, x0 + iw + 7, iy(d), { color: MUTED, size: 10, halo: HALO }); }
    text(ctx, title, x0, y0 - 9, { color: INK, size: 10, halo: HALO });
  }

  update();
}
