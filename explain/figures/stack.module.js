import * as THREE from '../../js/three.module.js';
import { slider, choice, legend, readout, text, rampRGB, MUTED, INK, paletteVersion } from '../runtime.module.js';
import { Scene3D, hexPatch, Arrows, hexTiles, linear } from '../scene3d.module.js';

const COLS = 7, ROWS = 6, K = 6, SPACING = 1.7, ARROW = 0.11;

export function mountStack(root) {
  const controls = root.querySelector('.controls');
  const patch = hexPatch(COLS, ROWS);
  let shown = K, mode = 'both', painted = -1;
  const layerY = (k) => (K - 1 - k + 0.5) * SPACING;
  const theta = (cell, k) => 288 + 9 * (K - 1 - k) * (1 + 0.04 * (K - 1 - k)) + 1.6 * cell.z + 0.9 * Math.sin(cell.x * 0.7);
  const wind = (x, z, k) => { const up = K - 1 - k; return [5 + 3.2 * up + 1.5 * Math.cos(z * 0.6), 0, 2.2 * Math.sin(x * 0.55 + up * 0.3)]; };

  const scene = new Scene3D(root, { height: 460, distance: 24, target: [0, 4.6, 0], yaw: -0.55, pitch: 0.3, update, draw: labels });
  const ground = scene.add(new THREE.Mesh(new THREE.CircleGeometry(9, 48).rotateX(-Math.PI / 2), new THREE.MeshLambertMaterial({ color: new THREE.Color(0.05, 0.045, 0.035) })));
  ground.position.y = -0.02;
  const tiles = hexTiles(scene.group, { thickness: 0.05, opacity: 0.92 }).mesh(COLS * ROWS * K);
  const arrows = new Arrows(scene.group, patch.edges.length * K, { radius: 0.03, head: 0.17, headRadius: 0.08 });
  legend(root, [['ramp', 'potential temperature θ at each cell’s center, 285 to 345 K', 'cool', 'warm', 'neutral'], ['arrow', 'the wind across each edge, the one number the model keeps there', INK]]);
  slider(controls, { label: 'Layers shown', min: 1, max: K, step: 1, value: shown, format: (v) => v === K ? `all ${K}` : `the lowest ${v}`, onInput: (v) => { shown = v; painted = -1; scene.fig.render(); } });
  choice(controls, { label: 'Show', options: [['both', 'both'], ['θ at the centers', 'theta'], ['wind across the edges', 'wind']], value: mode, onChange: (v) => { mode = v; painted = -1; scene.fig.render(); }, span: true });
  const out = readout(controls);
  out.set([['this patch', `${patch.cells.length} cells, ${patch.edges.length} edges, ${K} layers`], ['the model at N = 128', '163,842 cells and 491,520 edges, on 36 layers']]);

  const m = new THREE.Matrix4(), color = new THREE.Color(), rgb = [0, 0, 0];
  function update() {
    if (painted === paletteVersion) return;
    painted = paletteVersion;
    let n = 0;
    for (let k = K - shown; k < K; k++) for (const cell of patch.cells) {
      if (mode === 'wind') break;
      m.makeTranslation(cell.x, layerY(k), cell.z);
      tiles.setMatrixAt(n, m);
      rampRGB((theta(cell, k) - 285) / 60, rgb);
      tiles.setColorAt(n, color.setRGB(linear(rgb[0]), linear(rgb[1]), linear(rgb[2])));
      n++;
    }
    tiles.count = n; tiles.instanceMatrix.needsUpdate = true; if (tiles.instanceColor) tiles.instanceColor.needsUpdate = true;
    arrows.begin();
    if (mode !== 'theta') for (let k = K - shown; k < K; k++) for (const e of patch.edges) {
      const [u, , v] = wind(e.x, e.z, k), normal = u * e.nx + v * e.nz, len = normal * ARROW;
      arrows.push(e.x - e.nx * len / 2, layerY(k) + 0.12, e.z - e.nz * len / 2, e.nx * len, 0, e.nz * len, [0.93, 0.93, 0.93]);
    }
    arrows.end();
  }

  function labels(ctx) {
    const top = scene.project(0, layerY(K - shown) + 0.6, -5.2), bottom = scene.project(0, 0, -5.2);
    if (top.visible) text(ctx, shown === K ? 'top layer' : `layer ${K - shown + 1} of ${K}`, top.x, top.y, { align: 'center', color: MUTED, size: 11 });
    if (bottom.visible) text(ctx, 'the ground', bottom.x, bottom.y + 12, { align: 'center', color: MUTED, size: 11 });
    const north = scene.project(0, 0.1, -6.4);
    if (north.visible) text(ctx, 'north', north.x, north.y, { align: 'center', color: MUTED, size: 11 });
  }
}
