import * as THREE from '../../js/three.module.js';
import { slider, buttons, legend, readout, text, MUTED, INK, ACCENT, rampRGB } from '../runtime.module.js';
import { Scene3D, hexPatch, Arrows, hexTiles, linear, cssColor } from '../scene3d.module.js';

const K = 6, SPACING = 1.3, PRESETS = { rising: [2, 1.5, 0, 0, -1.5, -2], piling: [-1, -1, -1, -1, -1, -1], sinking: [-2, -1.5, 0, 0, 1.5, 2] };

export function mountContinuity(root) {
  const controls = root.querySelector('.controls');
  const patch = hexPatch(3, 3), center = patch.cells[4];
  const divergence = [...PRESETS.rising];
  const layerY = (k) => (K - 1 - k + 0.5) * SPACING, interfaceY = (k) => (K - k) * SPACING;
  let result = null;

  const scene = new Scene3D(root, { height: 440, distance: 19, target: [0, 4.2, 0], yaw: -0.5, pitch: 0.28, update, draw: labels });
  const ground = scene.add(new THREE.Mesh(new THREE.CircleGeometry(4.2, 48).rotateX(-Math.PI / 2), new THREE.MeshLambertMaterial({ color: new THREE.Color(0.05, 0.045, 0.035) })));
  ground.position.y = -0.02;
  const tiles = hexTiles(scene.group, { thickness: 0.04, opacity: 0.45 }).mesh(patch.cells.length * K);
  const flows = new Arrows(scene.group, 6 * K, { radius: 0.04, head: 0.2, headRadius: 0.1 });
  const lifts = new Arrows(scene.group, K + 1, { radius: 0.05, head: 0.22, headRadius: 0.12 });
  legend(root, [['arrow', 'air flowing across the column’s six sides in each layer', INK], ['arrow', 'air crossing between layers, πσ̇', ACCENT]]);
  const sliders = divergence.map((value, k) => slider(controls, { label: k === 0 ? 'Top layer' : k === K - 1 ? 'Lowest layer' : `Layer ${k + 1}`, min: -3, max: 3, step: 0.5, value, format: (v) => v === 0 ? 'balanced' : v > 0 ? `${v} out` : `${-v} in`, onInput: (v) => { divergence[k] = v; compute(); } }));
  buttons(controls, [['Rising column', () => preset('rising')], ['Air piling in', () => preset('piling')], ['Sinking column', () => preset('sinking')]]);
  const out = readout(controls);

  function preset(name) { PRESETS[name].forEach((v, k) => { divergence[k] = v; sliders[k].value = v; }); compute(); }

  function compute() {
    const dSigma = 1 / K, tendency = -divergence.reduce((s, d) => s + d * dSigma, 0), flux = [0];
    let sum = 0;
    for (let k = 0; k < K; k++) { sum += divergence[k] * dSigma; flux.push(-sum - (k + 1) * dSigma * tendency); }
    flux[K] = 0;
    result = { tendency, flux };
    let strongest = 0;
    for (const f of flux) if (Math.abs(f) > Math.abs(strongest)) strongest = f;
    out.set([['surface pressure', Math.abs(tendency) < 1e-9 ? 'steady' : `${tendency > 0 ? 'rising' : 'falling'} by ${Math.abs(tendency).toFixed(2)} hPa per hour`], ['strongest flow between layers', Math.abs(strongest) < 1e-9 ? 'none' : `${Math.abs(strongest).toFixed(2)} hPa per hour ${strongest < 0 ? 'upward' : 'downward'}`], ['at the top and the ground', 'zero, always']]);
    scene.fig.render();
  }

  const m = new THREE.Matrix4(), color = new THREE.Color(), rgb = [0, 0, 0];
  function update() {
    if (!result) return;
    let n = 0;
    for (let k = 0; k < K; k++) for (const cell of patch.cells) {
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
      flows.push(mx - nx * len / 2, layerY(k), mz - nz * len / 2, nx * len, 0, nz * len, [0.93, 0.93, 0.93]);
    }
    flows.end();
    lifts.begin();
    for (let k = 1; k < K; k++) { const f = result.flux[k], len = -0.55 * f; lifts.push(center.x, interfaceY(k) - len / 2, center.z, 0, len, 0, cssColor('rgb(255, 232, 160)')); }
    lifts.end();
  }

  function labels(ctx) {
    const top = scene.project(center.x, interfaceY(0) + 0.3, center.z), bottom = scene.project(center.x, -0.05, center.z);
    if (top.visible) text(ctx, 'top of the atmosphere, σ = 0', top.x, top.y - 10, { align: 'center', color: MUTED, size: 11 });
    if (bottom.visible) text(ctx, 'the ground, σ = 1', bottom.x, bottom.y + 16, { align: 'center', color: MUTED, size: 11 });
  }

  compute();
}
