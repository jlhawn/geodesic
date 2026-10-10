import { Figure, slider, legend, readout, text, arrow, anomalyColor, clamp, MUTED, LINE, GRID, ACCENT } from '../runtime.module.js';
import { R, CP, KAPPA, G, P0, pressureAt, heightOf, thetaAt } from '../physics.module.js';
import { sigmaInterfaces } from '../../js/dynamics/sigmaCore.module.js';

const LEVELS = sigmaInterfaces('bl36'), K = LEVELS.length - 1, DX = 1e5, COLUMNS = 24, WIDTH = 2.5e5;
const MID = Array.from({ length: K }, (_, k) => 0.5 * (LEVELS[k] + LEVELS[k + 1])), FIRST = MID.findIndex((s) => s * 1013.25 > 130);

export function mountMountain(root) {
  const controls = root.querySelector('.controls');
  let peak = 3000, layer = K - 4, edge = COLUMNS / 2 - 3, geometry = null;
  const sigmaK = Float64Array.from(LEVELS, (s) => s ** KAPPA), sigma1K = Float64Array.from(LEVELS, (s) => s ** (1 + KAPPA));
  const xOf = (i) => (i - COLUMNS / 2 + 0.5) * DX;
  const ground = (x) => peak * Math.exp(-((x / WIDTH) ** 2));

  function columns() {
    return Array.from({ length: COLUMNS }, (_, i) => {
      const h = ground(xOf(i)), pi = pressureAt(h), s = (pi / P0) ** KAPPA;
      const theta = new Float64Array(K), exner = new Float64Array(K), phi = new Float64Array(K), phiI = new Float64Array(K + 1);
      phiI[K] = G * h;
      for (let k = K - 1; k >= 0; k--) {
        const p = MID[k] * pi, dSigma = LEVELS[k + 1] - LEVELS[k];
        theta[k] = thetaAt(heightOf(p));
        const lower = s * sigmaK[k + 1], upper = s * sigmaK[k];
        exner[k] = s * (sigma1K[k + 1] - sigma1K[k]) / ((1 + KAPPA) * dSigma);
        phi[k] = phiI[k + 1] + CP * theta[k] * (lower - exner[k]);
        phiI[k] = phi[k] + CP * theta[k] * (exner[k] - upper);
      }
      return { h, pi, theta, exner, phi, phiI };
    });
  }

  const fig = new Figure(root, { height: 400, minHeight: 300, draw });
  legend(root, [['force', 'the push from the slope of the layer’s height', 'warm'], ['force', 'the push from the change in pressure along the layer, weighted by temperature', 'cool'], ['force', 'what is left over, drawn 100 times larger', ACCENT]]);
  slider(controls, { label: 'Mountain height', min: 0, max: 5000, step: 100, value: peak, format: (v) => `${(v / 1000).toFixed(1)} km`, onInput: (v) => { peak = v; fig.render(); } });
  slider(controls, { label: 'Layer', min: 1, max: K - FIRST, step: 1, value: K - layer, format: (v) => (v === 1 ? 'the lowest' : `${v} up from the ground`), onInput: (v) => { layer = K - v; fig.render(); } });
  const out = readout(controls);

  function draw(ctx, w, h) {
    const cols = columns(), left = 52, right = 14, top = 14, bottom = h - 30, zMax = 16000;
    const x = (m) => left + (m / (COLUMNS * DX) + 0.5) * (w - left - right), y = (z) => bottom - z / zMax * (bottom - top);
    geometry = { x, left, right, w };
    for (const p of [900, 800, 700, 600, 500, 400, 300, 200, 150]) { const z = heightOf(p * 100); ctx.strokeStyle = GRID; ctx.setLineDash([2, 4]); ctx.beginPath(); ctx.moveTo(left, y(z)); ctx.lineTo(w - right, y(z)); ctx.stroke(); ctx.setLineDash([]); text(ctx, `${p} hPa`, left - 6, y(z), { align: 'right', color: MUTED, size: 10 }); }
    for (let k = FIRST; k <= K; k++) {
      ctx.strokeStyle = k === layer || k === layer + 1 ? 'rgba(255,255,255,0.75)' : LINE; ctx.lineWidth = k === layer || k === layer + 1 ? 1.6 : 1;
      ctx.beginPath(); cols.forEach((c, i) => { const yy = y(c.phiI[k] / G); if (i) ctx.lineTo(x(xOf(i)), yy); else ctx.moveTo(x(xOf(i)), yy); }); ctx.stroke();
    }
    ctx.fillStyle = '#5a4a36'; ctx.beginPath(); ctx.moveTo(left, bottom);
    for (let s = 0; s <= 200; s++) { const m = (s / 200 - 0.5) * COLUMNS * DX; ctx.lineTo(x(m), y(ground(m))); }
    ctx.lineTo(w - right, bottom); ctx.closePath(); ctx.fill();
    const a = cols[edge], b = cols[edge + 1], k = layer;
    const termA = -(b.phi[k] - a.phi[k]) / DX, termB = -R * 0.5 * (a.theta[k] * a.exner[k] + b.theta[k] * b.exner[k]) * Math.log(b.pi / a.pi) / DX, sum = termA + termB;
    const ex = x((xOf(edge) + xOf(edge + 1)) / 2), ey = y((a.phi[k] + b.phi[k]) / (2 * G)), scale = 120 / Math.max(0.05, Math.abs(termA), Math.abs(termB));
    arrow(ctx, ex, ey - 6, ex + termA * scale, ey - 6, { color: anomalyColor(1), width: 2, head: 8, dash: [5, 4], open: true });
    arrow(ctx, ex, ey + 6, ex + termB * scale, ey + 6, { color: anomalyColor(-1), width: 2, head: 8, dash: [5, 4], open: true });
    arrow(ctx, ex, ey + 18, ex + clamp(sum * scale * 100, -200, 200), ey + 18, { color: ACCENT, width: 2, head: 8, dash: [5, 4], open: true });
    ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(ex, ey, 3.5, 0, Math.PI * 2); ctx.fill();
    text(ctx, 'drag along the mountain to move the point', w / 2, top + 4, { align: 'center', color: MUTED, size: 11 });
    const hour = (v) => `${(v * 3600).toFixed(Math.abs(v * 3600) < 10 ? 2 : 0)} m/s`;
    out.set([['the layer’s height', `${hour(termA)} an hour`], ['the pressure change', `${hour(termB)} an hour`], ['left over', `${hour(sum)} an hour, ${(Math.abs(sum) * 86400).toFixed(1)} m/s after a day`], ['the point', `${Math.round(MID[k] * (a.pi + b.pi) / 200)} hPa, ${((a.phi[k] + b.phi[k]) / (2 * G * 1000)).toFixed(1)} km up`], ['the right answer', 'zero: the air is at rest and the pressure surfaces are flat']]);
  }

  fig.pointer({ down: move, move });
  function move({ x: px }) {
    if (!geometry) return;
    const m = ((px - geometry.left) / (geometry.w - geometry.left - geometry.right) - 0.5) * COLUMNS * DX;
    edge = clamp(Math.round(m / DX + COLUMNS / 2 - 1), 0, COLUMNS - 2);
    fig.render();
  }
}
