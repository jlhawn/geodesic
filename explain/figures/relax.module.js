import * as THREE from '../../js/three.module.js';
import { Grid } from '../../js/grid.module.js';
import { Globe, rgb } from '../globe.module.js';
import { slider, buttons, legend, readout, rampColor } from '../runtime.module.js';

export function mountRelax(root) {
  const controls = root.querySelector('.controls');
  const globe = new Globe(root, { height: 440, zoom: 0.42, lookAt: new THREE.Vector3(0, 0, 1) });
  legend(root, [['ramp', 'cell edges, from cut at the midpoint to cut 8% or more off it', 'neutral', 'warm'], ['faint', 'lines between neighboring cell centers', 'rgba(255,255,255,0.35)']]);
  let N = 10, grid = new Grid(N, { relax: 0 }), passes = 0;
  const reset = () => { grid = new Grid(N, { relax: 0 }); passes = 0; globe.zoom = 0.42 * Math.sqrt(6 / N); build(); };
  let pending = 0;
  slider(controls, { label: 'Subdivisions N', min: 10, max: 32, step: 1, value: N, format: (v) => `${v}`, onInput: (v) => { N = v; clearTimeout(pending); pending = setTimeout(reset, 150); } });
  buttons(controls, [
    ['Relax once more', () => { pass(1); }],
    ['Ten more passes', () => { pass(10); }],
    ['Reset', reset],
  ]);
  const out = readout(controls);

  function pass(n) {
    for (let i = 0; i < n; i++) grid.relax();
    for (const cell of grid) cell.calculateArea();
    passes += n;
    build();
  }

  function build() {
    const cells = [...grid];
    globe.clear();
    globe.cells(cells.map((c) => ({ center: c.centerVertex, vertices: c.vertices })), () => rgb('rgb(70, 74, 84)'));
    const edges = [], duals = [];
    let worst = 0, sum = 0, count = 0;
    for (const cell of cells) for (const nb of cell.neighbors) {
      if (nb.index < cell.index) continue;
      const [a, b] = cell.vertices.filter((v) => nb.vertices.includes(v));
      const offset = cell.centerVertex.clone().add(nb.centerVertex).normalize().distanceTo(a.clone().add(b).normalize()) / a.distanceTo(b);
      worst = Math.max(worst, offset); sum += offset; count++;
      edges.push({ a, b, color: rgb(rampColor(0.5 + Math.min(1, offset / 0.08) / 2)) });
      duals.push([cell.centerVertex, nb.centerVertex]);
    }
    globe.lines(duals, { color: 0xffffff, opacity: 0.25, lift: 1.002 });
    globe.coloredLines(edges, { lift: 1.004 });
    out.set([['cells', `${cells.length}`], ['Lloyd passes', `${passes}`], ['worst edge', `cut ${(worst * 100).toFixed(1)}% from its midpoint`], ['typical edge', `${(sum / count * 100).toFixed(2)}%`]]);
    globe.fig.render();
  }

  reset();
}
