import * as THREE from '../../js/three.module.js';
import { slider, buttons, legend, readout, text, termColor, clamp, MUTED, rampRGB } from '../runtime.module.js';
import { Scene3D, hexPatch, Arrows, hexTiles, linear, cssColor } from '../scene3d.module.js';

const K = 6, SPACING = 1.3, GAUGE = { x: 2.7, z: 1.4 }, PERIOD = 16, LAG = Math.PI / 2, SWING = [2, 1.5, 0.75, -0.75, -1.5, -2], CYCLE = 'Rising and sinking, out of step', PRESETS = { rising: [2, 1.5, 0, 0, -1.5, -2], piling: [-1, -1, -1, -1, -1, -1], draining: [1, 1, 1, 1, 1, 1], sinking: [-2, -1.5, 0, 0, 1.5, 2] };

export function mountContinuity(root) {
  const controls = root.querySelector('.controls'), A = termColor('a'), B = termColor('b'), C = termColor('c'), weight = cssColor(A), sides = cssColor(B), between = cssColor(C);
  const patch = hexPatch(3, 3), center = patch.cells[4], cells = patch.cells.filter((cell) => Math.hypot(cell.x - center.x, cell.z - center.z) < 1.9), gx = center.x + GAUGE.x, gz = center.z + GAUGE.z;
  const divergence = [...PRESETS.rising];
  const layerY = (k) => (K - 1 - k + 0.5) * SPACING, interfaceY = (k) => (K - k) * SPACING;
  let result = null, cycling = false, phase = 0, toggle = null;

  const scene = new Scene3D(root, { height: 440, distance: 19, target: [center.x, 4.2, center.z], yaw: -0.5, pitch: 0.28, update, draw: labels, step });
  const ground = scene.add(new THREE.Mesh(new THREE.CircleGeometry(4.2, 48).rotateX(-Math.PI / 2), new THREE.MeshLambertMaterial({ color: new THREE.Color(0.05, 0.045, 0.035) })));
  ground.position.set(center.x, -0.02, center.z);
  const tiles = hexTiles(scene.group, { thickness: 0.04, opacity: 0.45, depthWrite: false }).mesh(cells.length * K);
  const flows = new Arrows(scene.group, 6 * K, { radius: 0.04, head: 0.2, headRadius: 0.1 });
  const lifts = new Arrows(scene.group, K + 1, { radius: 0.05, head: 0.22, headRadius: 0.12 });
  const gauge = new Arrows(scene.group, 1, { radius: 0.09, head: 0.36, headRadius: 0.2 });
  legend(root, [['arrow', 'the surface pressure rising or falling, ∂π/∂t: the ground glows as it changes', A], ['arrow', 'air flowing across the column’s six sides in each layer', B], ['arrow', 'air crossing between layers, πσ̇', C], ['ramp', 'the column’s layers, tinted from air flowing in to air flowing out', 'cool', 'warm', 'neutral']]);
  const sliders = divergence.map((value, k) => slider(controls, { label: k === 0 ? 'Top layer' : k === K - 1 ? 'Lowest layer' : `Layer ${k + 1}`, min: -3, max: 3, step: 0.5, value, format: (v) => v === 0 ? 'balanced' : v > 0 ? `${v} out` : `${-v} in`, onInput: (v) => { stop(); divergence[k] = v; compute(); } }));
  [, , , , toggle] = buttons(controls, [['Rising column', () => preset('rising')], ['Sinking column', () => preset('sinking')], ['Air piling in', () => preset('piling')], ['Air draining out', () => preset('draining')], [CYCLE, () => (cycling ? stop() : start())]]);
  const out = readout(controls);

  function preset(name) { stop(); PRESETS[name].forEach((v, k) => { divergence[k] = v; sliders[k].value = v; }); compute(); }
  function start() { cycling = true; phase = 0; toggle.textContent = 'Stop'; scene.fig.play(true); }
  function stop() { if (!cycling) return; cycling = false; toggle.textContent = CYCLE; scene.fig.play(false); }
  function step(dt) {
    if (!cycling) return false;
    phase += 2 * Math.PI * dt / PERIOD;
    SWING.forEach((a, k) => { const v = a * Math.cos(phase - LAG * k / (K - 1)); divergence[k] = v; sliders[k].value = Math.round(v * 2) / 2; });
    compute(false);
  }

  function compute(render = true) {
    const dSigma = 1 / K, tendency = -divergence.reduce((s, d) => s + d * dSigma, 0), flux = [0];
    let sum = 0;
    for (let k = 0; k < K; k++) { sum += divergence[k] * dSigma; flux.push(-sum - (k + 1) * dSigma * tendency); }
    flux[K] = 0;
    result = { tendency, flux };
    let strongest = 0;
    for (const f of flux) if (Math.abs(f) > Math.abs(strongest)) strongest = f;
    out.set([['surface pressure', Math.abs(tendency) < 1e-9 ? 'steady' : `${tendency > 0 ? 'rising' : 'falling'} by ${Math.abs(tendency).toFixed(2)} hPa per hour`], ['strongest flow between layers', Math.abs(strongest) < 1e-9 ? 'none' : `${Math.abs(strongest).toFixed(2)} hPa per hour ${strongest < 0 ? 'upward' : 'downward'}`], ['at the top and the ground', 'zero, always']]);
    if (render) scene.fig.render();
  }

  const m = new THREE.Matrix4(), color = new THREE.Color(), rgb = [0, 0, 0];
  function update() {
    if (!result) return;
    let n = 0;
    for (let k = 0; k < K; k++) for (const cell of cells) {
      m.makeTranslation(cell.x, layerY(k) - SPACING / 2 + 0.02, cell.z);
      tiles.setMatrixAt(n, m);
      const c = cell === center ? rampRGB(0.5 + divergence[k] / 6, rgb) : [0.32, 0.33, 0.36];
      tiles.setColorAt(n, color.setRGB(linear(c[0]), linear(c[1]), linear(c[2])));
      n++;
    }
    tiles.count = n; tiles.instanceMatrix.needsUpdate = true; if (tiles.instanceColor) tiles.instanceColor.needsUpdate = true;
    flows.begin();
    for (let k = 0; k < K; k++) for (const { edge, sign } of center.edges) {
      const nx = edge.nx * sign, nz = edge.nz * sign, len = 0.32 * divergence[k], mx = (center.x + edge.x) / 2 + nx * 0.45, mz = (center.z + edge.z) / 2 + nz * 0.45;
      flows.push(mx - nx * len / 2, layerY(k), mz - nz * len / 2, nx * len, 0, nz * len, sides);
    }
    flows.end();
    lifts.begin();
    for (let k = 1; k < K; k++) { const f = result.flux[k], len = -0.55 * f; lifts.push(center.x, interfaceY(k) - len / 2, center.z, 0, len, 0, between); }
    lifts.end();
    const t = result.tendency, glow = Math.min(1, Math.abs(t) / 1.5) * 0.6, base = [0.05, 0.045, 0.035];
    ground.material.color.setRGB(...base.map((b, i) => b + (linear(weight[i]) - b) * glow));
    gauge.begin();
    if (Math.abs(t) > 1e-9) { const len = clamp(1.1 * t, -3.3, 3.3); gauge.push(gx, len < 0 ? -len : 0, gz, 0, len, 0, weight); }
    gauge.end();
  }

  function labels(ctx) {
    const top = scene.project(center.x, interfaceY(0) + 0.3, center.z), bottom = scene.project(center.x, -0.05, center.z);
    if (top.visible) text(ctx, 'top of the atmosphere, σ = 0', top.x, top.y - 10, { align: 'center', color: MUTED, size: 11 });
    if (bottom.visible) text(ctx, 'the ground, σ = 1', bottom.x, bottom.y + 16, { align: 'center', color: MUTED, size: 11 });
    const t = result?.tendency ?? 0, steady = Math.abs(t) < 1e-9, at = scene.project(gx, steady ? 0.1 : Math.abs(clamp(1.1 * t, -3.3, 3.3)) + 0.45, gz);
    const style = { align: 'center', color: steady ? MUTED : A, size: 11, halo: 'rgba(20, 20, 22, 0.85)' };
    if (at.visible) { text(ctx, 'surface pressure', at.x, at.y - 20, style); text(ctx, steady ? 'steady' : `${t > 0 ? 'rising' : 'falling'} ${Math.abs(t).toFixed(2)} hPa an hour`, at.x, at.y - 6, style); }
  }

  compute();
}
