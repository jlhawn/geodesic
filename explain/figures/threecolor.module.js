import { gradient, tangential } from '../../js/dynamics/operators.module.js';
import { hexTorus, threeColors, broadHigh, linearWater } from '../hexTorus.module.js';
import { Figure, choice, buttons, legend, caption, readout, text, arrow, anomalyColor, ACCENT, INK, MUTED, LINE } from '../runtime.module.js';

const SQRT3 = Math.sqrt(3), SPACING = 60e3, LATITUDE = 45, F = 2 * 7.2921e-5 * Math.sin(LATITUDE * Math.PI / 180), WIDE = 560, TOP = 24, TRACE = 84, GAP = 34, BOTTOM = 22;
const FLUIDS = {
  sea: { c: 2, rho: 1025, dt: 600, every: 1, window: 96 * 3600, speed: 12 * 3600, tick: 24 },
  air: { c: 50, rho: 1.2, dt: 30, every: 4, window: 8 * 3600, speed: 3600, tick: 2 },
};
const fixed = (v, n = 2) => (Math.abs(v) < 0.5 * 10 ** -n ? (0).toFixed(n) : v.toFixed(n).replace('-', '−'));

export function mountThreeColor(root) {
  const controls = root.querySelector('.controls');
  let patternName = 'three', fluidName = 'sea', torus = null, start = null, runs = [], t = 0, carry = 0, samples = [], fastest = 0, frame = null, steps = 0, mode = null;
  const pushes = { across: null, along: null };

  const fig = new Figure(root, { height: (w) => layout(w).h, draw, step });
  caption(root, `A patch of cells ${SPACING / 1000} km apart, the model’s spacing at N = 128, seen from above with north up and wrapped around, so that what leaves one side comes back in at the other. The planet turns as it does at ${LATITUDE}° north. The pushes and the sum along each edge are the model’s own code, run here on a single layer of fluid.`);
  const key = legend(root, []);
  choice(controls, { label: 'Start from', options: [['three colors', 'three'], ['one high', 'one'], ['a broad high', 'high'], ['level', 'level']], value: patternName, onChange: (v) => { patternName = v; fill(); }, span: true });
  choice(controls, { label: 'The patch holds', options: [['the sea', 'sea'], ['the air', 'air']], value: fluidName, onChange: (v) => { fluidName = v; reset(); }, span: true });
  const [go] = buttons(controls, [['Let go', () => { fig.play(!fig.running); label(); }], ['Reset', () => reset()]]);
  const out = readout(controls);
  const label = () => { go.textContent = fig.running ? 'Pause' : t > 0 ? 'Go on' : 'Let go'; };

  function layout(w) {
    const wide = w >= WIDE, cols = wide ? 9 : 6, rows = wide ? 6 : 8, R = Math.min(34, (w - 24) / (SQRT3 * (cols + 0.5)));
    const patch = (1.5 * (rows - 1) + 2) * R;
    return { cols, rows, R, patch, h: Math.round(TOP + patch + GAP + TRACE + BOTTOM) };
  }

  function rebuild(cols, rows) {
    torus = hexTorus(cols, rows, SPACING);
    fill();
  }

  function fill() {
    if (!torus) return;
    const { mesh } = torus, center = (torus.rows >> 1) * torus.cols + (torus.cols >> 1);
    start = patternName === 'three' ? threeColors(torus) : patternName === 'high' ? broadHigh(torus, 1.2 * SPACING) : Float64Array.from({ length: mesh.nCells }, (_, i) => (patternName === 'one' && i === center ? 1 : 0));
    reset();
  }

  function reset() {
    if (!torus) return;
    const { mesh } = torus, fluid = FLUIDS[fluidName];
    fig.play(false);
    runs = [F, 0].map((f) => ({ step: linearWater(mesh, { c: fluid.c, f, rho: fluid.rho }), state: { p: Float64Array.from(start), u: new Float64Array(mesh.nEdges) } }));
    t = 0; carry = 0; steps = 0; samples = [[0, remaining(start), remaining(start)]]; fastest = 0;
    pushes.across = gradient(mesh, start).map((g) => -g / fluid.rho);
    pushes.along = tangential(mesh, pushes.across);
    label(); show(); fig.render();
  }

  const mean = (p) => p.reduce((a, b) => a + b, 0) / p.length;
  function remaining(p) {
    const m0 = mean(start), m = mean(p);
    let num = 0, den = 0;
    for (let i = 0; i < p.length; i++) { num += (p[i] - m) * (start[i] - m0); den += (start[i] - m0) ** 2; }
    return den ? num / den : 0;
  }

  function step(dt) {
    const fluid = FLUIDS[fluidName];
    carry += dt * fluid.speed;
    while (carry >= fluid.dt) {
      carry -= fluid.dt;
      for (const { step: advance, state } of runs) advance(state, fluid.dt);
      t += fluid.dt;
      if (++steps % fluid.every === 0) samples.push([t, remaining(runs[0].state.p), remaining(runs[1].state.p)]);
    }
    while (samples.length > 2 && samples[1][0] < t - fluid.window) samples.shift();
    fastest = Math.max(Math.max(...Array.from(runs[0].state.u, Math.abs)), fastest * Math.exp(-dt * fluid.speed / (fluid.window / 8)));
    show();
  }

  const level = () => { const m = mean(start); return start.every((v) => Math.abs(v - m) < 1e-9); };

  function show() {
    if (!torus) return;
    const fluid = FLUIDS[fluidName], radius = fluid.c / F, ratio = radius / SPACING;
    const deformation = ['deformation radius', `${Math.round(radius / 1000)} km, ${ratio < 1 ? `${fixed(ratio)} of` : `${ratio.toFixed(1)} times`} the spacing`];
    const now = t === 0 ? 'rest' : 'run';
    if (now !== mode && now === 'rest') key.set([['gradient', 'pressure, lower to higher', 'cool', 'warm'], ['force', 'the push across each edge, from the pressures on either side', INK], ['force', 'the push along each edge, from the sum of the ten around it', ACCENT]]);
    else if (now !== mode) key.set([['gradient', 'pressure, lower to higher', 'cool', 'warm'], ['arrow', 'the wind across each edge', INK], ['arrow', 'the wind along each edge, from the sum', ACCENT], ['line', 'how much of the pattern is left, the planet turning', ACCENT], ['dash', 'and the planet stopped', INK]]);
    mode = now;
    if (t === 0) {
      const across = Math.max(...Array.from(pushes.across, Math.abs)), along = Math.max(...Array.from(pushes.along, Math.abs));
      out.set([['the sum finds a push along the edges', level() ? 'none: the pressure is level and nothing pushes' : along < 1e-9 * across ? 'none, on any edge' : `up to ${Math.round(100 * along / across)}% of the strongest push across one`], deformation]);
      return;
    }
    const [, on, off] = samples[samples.length - 1];
    out.set([['after', `${fixed(t / 3600, t < 36000 ? 1 : 0)} hours`], ['pattern left', `${fixed(on)} with the planet turning, ${fixed(off)} with it stopped`], deformation], 3);
  }

  function draw(ctx, w, h) {
    const L = layout(w);
    if (!torus || torus.cols !== L.cols) rebuild(L.cols, L.rows);
    const { mesh, xy, cols, rows } = torus, R = L.R, pitch = SQRT3 * R;
    const x0 = (w - pitch * (cols + 0.5)) / 2 + pitch / 2, y0 = TOP + R;
    const sx = (x) => x0 + x / SPACING * pitch, sy = (y) => y0 - y / SPACING * pitch;
    frame = { sx, sy, R };
    const p = t === 0 ? start : runs[0].state.p, m = mean(start), mp = mean(p), swing = Math.max(...Array.from(start, (v) => Math.abs(v - m))) || 1;
    const hex = (i) => { ctx.beginPath(); for (let k = 0; k < 6; k++) { const a = Math.PI / 6 + k * Math.PI / 3; ctx.lineTo(sx(xy[i][0]) + R * Math.cos(a), sy(xy[i][1]) - R * Math.sin(a)); } ctx.closePath(); };
    for (let i = 0; i < mesh.nCells; i++) { hex(i); ctx.fillStyle = anomalyColor((p[i] - mp) / swing, 0.75); ctx.fill(); ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.stroke(); }
    ctx.fillStyle = '#fff';
    for (const [x, y] of xy) { ctx.beginPath(); ctx.arc(sx(x), sy(y), 2, 0, Math.PI * 2); ctx.fill(); }

    const resting = t === 0, across = resting ? pushes.across : runs[0].state.u, along = resting ? pushes.along : tangential(mesh, across);
    const biggest = resting ? Math.max(...Array.from(across, Math.abs)) : fastest, per = biggest > 0 ? 0.9 * R / biggest : 0;
    const acrossStyle = resting ? { color: INK, width: 1.5, head: 6, dash: [3, 2.5], open: true } : { color: INK, width: 1.5, head: 6 };
    const alongStyle = resting ? { color: ACCENT, width: 2, head: 6, dash: [3, 2.5], open: true } : { color: ACCENT, width: 2, head: 6 };
    for (let i = 0; i < mesh.nCells; i++) for (let k = 0; k < 6; k++) {
      const j = torus.neighbor(i, k), a = k * Math.PI / 3, dx = Math.cos(a), dy = Math.sin(a);
      const wrapped = Math.hypot(xy[j][0] - xy[i][0] - SPACING * dx, xy[j][1] - xy[i][1] - SPACING * dy) > 1;
      if (k >= 3 && !wrapped) continue;
      const e = mesh.edgesOnCell[6 * i + k], nx = mesh.nEdge[3 * e], ny = mesh.nEdge[3 * e + 1];
      const mx = sx(xy[i][0] + SPACING / 2 * dx), my = sy(xy[i][1] + SPACING / 2 * dy);
      const lu = across[e] * per, lt = along[e] * per;
      if (Math.abs(lu) > 1.5) arrow(ctx, mx - nx * lu / 2, my + ny * lu / 2, mx + nx * lu / 2, my - ny * lu / 2, acrossStyle);
      if (Math.abs(lt) > 1.5) arrow(ctx, mx + ny * lt / 2, my + nx * lt / 2, mx - ny * lt / 2, my - nx * lt / 2, alongStyle);
    }
    if (patternName === 'three' && resting && start.some((v) => v)) {
      const c = (rows >> 1) * cols + (cols >> 1), vx = sx(xy[c][0]) + R * Math.cos(Math.PI / 6), vy = sy(xy[c][1]) - R * Math.sin(Math.PI / 6);
      ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(vx, vy, 0.42 * R, 0, Math.PI * 2); ctx.stroke();
      text(ctx, 'one of each color', vx + 0.5 * R, vy - 0.55 * R, { color: '#fff', size: 11, weight: 700, halo: 'rgba(20,20,22,0.8)' });
    }
    if (resting) text(ctx, 'tap a cell to change its pressure', w / 2, 12, { align: 'center', color: MUTED, size: 11 });
    drawTrace(ctx, w, TOP + L.patch + GAP);
  }

  function drawTrace(ctx, w, top) {
    const fluid = FLUIDS[fluidName], left = 34, right = w - 12, from = Math.max(0, t - fluid.window), to = from + fluid.window;
    const px = (s) => left + (s - from) / (to - from) * (right - left), py = (a) => top + (1 - a) / 2 * TRACE;
    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    for (const a of [1, 0, -1]) { ctx.beginPath(); ctx.moveTo(left, py(a)); ctx.lineTo(right, py(a)); ctx.globalAlpha = a ? 0.4 : 1; ctx.stroke(); }
    ctx.globalAlpha = 1;
    text(ctx, '+1', left - 6, py(1), { align: 'right', color: MUTED, size: 10 });
    text(ctx, '0', left - 6, py(0), { align: 'right', color: MUTED, size: 10 });
    text(ctx, '−1', left - 6, py(-1), { align: 'right', color: MUTED, size: 10 });
    text(ctx, 'the pattern as it started, at +1, and turned inside out, at −1', left, top - 14, { color: MUTED, size: 11 });
    for (let hour = Math.ceil(from / 3600 / fluid.tick) * fluid.tick; hour * 3600 <= to; hour += fluid.tick) text(ctx, `${hour} h`, px(hour * 3600), top + TRACE + 11, { align: 'center', color: MUTED, size: 10 });
    for (const [column, style] of [[1, { color: ACCENT, width: 2, dash: null }], [2, { color: INK, width: 1.5, dash: [4, 3] }]]) {
      ctx.save(); ctx.strokeStyle = style.color; ctx.lineWidth = style.width; if (style.dash) ctx.setLineDash(style.dash);
      ctx.beginPath();
      samples.forEach((s, n) => { const x = px(s[0]), y = py(s[column]); if (n) ctx.lineTo(x, y); else ctx.moveTo(x, y); });
      ctx.stroke(); ctx.restore();
    }
  }

  function cellAt({ x, y }) {
    if (!frame || !torus) return null;
    let best = null, near = 0.9 * frame.R;
    torus.xy.forEach(([cx, cy], i) => { const d = Math.hypot(x - frame.sx(cx), y - frame.sy(cy)); if (d < near) { near = d; best = i; } });
    return best;
  }
  fig.pointer({
    hit: (at) => cellAt(at) !== null,
    down: (at) => {
      const i = cellAt(at);
      if (i === null) return false;
      start[i] = start[i] > 0.5 ? -1 : start[i] < -0.5 ? 0 : 1;
      reset();
    },
  });
}
