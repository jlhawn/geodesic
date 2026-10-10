import { Figure, choice, buttons, legend, readout, text, arrow, rampRGB, MUTED, INK, ACCENT, anomalyColor } from '../runtime.module.js';

const N = 24;

export function mountInvariant(root) {
  const controls = root.querySelector('.controls');
  const u = new Float64Array(N * N), v = new Float64Array(N * N);
  let show = 'parts', last = null, geometry = null;
  const at = (i, j) => ((j + N) % N) * N + ((i + N) % N);
  const PRESETS = {
    vortex: (x, y) => { const r2 = x * x + y * y, f = Math.exp(-r2 / 18); return [-y * f * 0.5, x * f * 0.5]; },
    jet: (x, y) => [2.2 * Math.exp(-(y * y) / 8), 0],
    pair: (x, y) => { const a = PRESETS.vortex(x - 4, y), b = PRESETS.vortex(x + 4, y); return [a[0] - b[0], a[1] - b[1]]; },
    calm: () => [0, 0],
  };
  function preset(name) { for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) { const [a, b] = PRESETS[name](i - N / 2 + 0.5, j - N / 2 + 0.5); u[at(i, j)] = a; v[at(i, j)] = b; } fig.render(); }

  function terms() {
    const zeta = new Float64Array(N * N), K = new Float64Array(N * N), adv = new Float64Array(2 * N * N), vort = new Float64Array(2 * N * N), grad = new Float64Array(2 * N * N);
    for (let n = 0; n < N * N; n++) K[n] = 0.5 * (u[n] * u[n] + v[n] * v[n]);
    for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
      const n = at(i, j), dx = (f, a, b) => (f[at(a + 1, b)] - f[at(a - 1, b)]) / 2, dy = (f, a, b) => (f[at(a, b + 1)] - f[at(a, b - 1)]) / 2;
      zeta[n] = dx(v, i, j) - dy(u, i, j);
    }
    for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
      const n = at(i, j), dx = (f) => (f[at(i + 1, j)] - f[at(i - 1, j)]) / 2, dy = (f) => (f[at(i, j + 1)] - f[at(i, j - 1)]) / 2;
      adv[2 * n] = -(u[n] * dx(u) + v[n] * dy(u)); adv[2 * n + 1] = -(u[n] * dx(v) + v[n] * dy(v));
      vort[2 * n] = zeta[n] * v[n]; vort[2 * n + 1] = -zeta[n] * u[n];
      grad[2 * n] = -dx(K); grad[2 * n + 1] = -dy(K);
    }
    return { zeta, adv, vort, grad };
  }

  const fig = new Figure(root, { height: 420, minHeight: 320, draw });
  root.classList.add('drag');
  legend(root, [['ramp', 'the flow’s spin, clockwise to counterclockwise', 'cool', 'warm', 'neutral'], ['arrow', 'the wind', INK], ['force', 'the push the flow gives itself, −(u·∇)u', 'rgb(255, 232, 160)'], ['force', 'its vortex part, −ζ k×u', 'warm'], ['force', 'its kinetic-energy part, −∇K', 'cool']]);
  const start = choice(controls, { label: 'Start with', options: [['a vortex', 'vortex'], ['a straight jet', 'jet'], ['two vortices', 'pair'], ['calm air', 'calm']], value: 'vortex', onChange: preset, span: true });
  choice(controls, { label: 'Show', options: [['the two parts', 'parts'], ['the whole push', 'whole'], ['just the wind', 'wind']], value: show, onChange: (s) => { show = s; fig.render(); }, span: true });
  const out = readout(controls);

  function draw(ctx, w, h) {
    const side = Math.min(w - 24, h - 24), x0 = (w - side) / 2, y0 = (h - side) / 2, cell = side / N;
    geometry = { x0, y0, cell };
    const { zeta, adv, vort, grad } = terms(), rgb = [0, 0, 0];
    let zmax = 1e-9; for (const z of zeta) zmax = Math.max(zmax, Math.abs(z));
    for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) { rampRGB(0.5 + 0.5 * zeta[at(i, j)] / Math.max(zmax, 0.3), rgb); ctx.fillStyle = `rgb(${rgb.map((c) => Math.round(c * 255)).join(',')})`; ctx.fillRect(x0 + i * cell, y0 + (N - 1 - j) * cell, cell + 0.5, cell + 0.5); }
    let amax = 1e-9, wmax = 1e-9, mismatch = 0, total = 0;
    for (let n = 0; n < N * N; n++) { amax = Math.max(amax, Math.hypot(adv[2 * n], adv[2 * n + 1]), Math.hypot(vort[2 * n], vort[2 * n + 1]), Math.hypot(grad[2 * n], grad[2 * n + 1])); wmax = Math.max(wmax, Math.hypot(u[n], v[n])); mismatch += Math.hypot(adv[2 * n] - vort[2 * n] - grad[2 * n], adv[2 * n + 1] - vort[2 * n + 1] - grad[2 * n + 1]); total += Math.hypot(vort[2 * n], vort[2 * n + 1]) + Math.hypot(grad[2 * n], grad[2 * n + 1]); }
    const ws = 1.6 * cell / wmax, fs = 1.6 * cell / amax;
    for (let j = 1; j < N; j += 2) for (let i = 1; i < N; i += 2) {
      const n = at(i, j), cx = x0 + (i + 0.5) * cell, cy = y0 + (N - 0.5 - j) * cell;
      if (show === 'wind' || show === 'whole' || show === 'parts') arrow(ctx, cx, cy, cx + u[n] * ws, cy - v[n] * ws, { color: show === 'wind' ? INK : 'rgba(255,255,255,0.45)', width: 1.4, head: 5 });
      if (show === 'whole') arrow(ctx, cx, cy, cx + adv[2 * n] * fs, cy - adv[2 * n + 1] * fs, { color: ACCENT, width: 1.8, head: 6, dash: [4, 3], open: true });
      if (show === 'parts') {
        arrow(ctx, cx, cy, cx + vort[2 * n] * fs, cy - vort[2 * n + 1] * fs, { color: anomalyColor(1), width: 1.8, head: 6, dash: [4, 3], open: true });
        arrow(ctx, cx, cy, cx + grad[2 * n] * fs, cy - grad[2 * n + 1] * fs, { color: anomalyColor(-1), width: 1.8, head: 6, dash: [4, 3], open: true });
      }
    }
    text(ctx, 'drag across the flow to stir it', w / 2, y0 - 6, { align: 'center', color: MUTED, size: 11 });
    let advMax = 0, partMax = 0;
    for (let n = 0; n < N * N; n++) { advMax = Math.max(advMax, Math.hypot(adv[2 * n], adv[2 * n + 1])); partMax = Math.max(partMax, Math.hypot(vort[2 * n], vort[2 * n + 1])); }
    out.set([['largest whole push', advMax.toFixed(3)], ['largest vortex part', partMax.toFixed(3)], ['the parts add up to the whole', total > 1e-9 ? `to within ${(100 * mismatch / total).toFixed(0)}% of their size on this coarse grid, exactly on paper` : 'there is no flow']]);
  }

  fig.pointer({
    down: (p) => { last = p; },
    move: (p) => {
      if (!geometry || !last) return;
      const { x0, y0, cell } = geometry, gx = (p.x - x0) / cell - 0.5, gy = N - 0.5 - (p.y - y0) / cell, dx = (p.x - last.x) / cell, dy = -(p.y - last.y) / cell;
      last = p;
      for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
        let ddx = i - gx, ddy = j - gy; ddx -= N * Math.round(ddx / N); ddy -= N * Math.round(ddy / N);
        const wgt = Math.exp(-(ddx * ddx + ddy * ddy) / 6) * 0.25;
        u[at(i, j)] += wgt * dx; v[at(i, j)] += wgt * dy;
      }
      fig.render();
    },
  });
  buttons(controls, [['Reset', () => preset(start.value)]]);
  preset('vortex');
}
