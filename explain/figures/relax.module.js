import * as THREE from '../../js/three.module.js';
import { Grid } from '../../js/grid.module.js';
import { Globe, rgb } from '../globe.module.js';
import { choice, legend, readout, rampColor } from '../runtime.module.js';

export function mountRelax(root) {
  const controls = root.querySelector('.controls');
  const globe = new Globe(root, { height: 440, zoom: 0.42, lookAt: new THREE.Vector3(0, 0, 1) });
  legend(root, [['ramp', 'cell edges, from cut at the midpoint to cut 8% or more off it', 'neutral', 'warm'], ['faint', 'lines between neighboring cell centers', 'rgba(255,255,255,0.35)']]);
  const grids = { off: new Grid(6, { relax: 0 }), on: new Grid(6, { relax: 8 }) };
  let relaxed = 'off';
  choice(controls, { label: 'Lloyd relaxation', options: [['off', 'off'], ['on', 'on']], value: relaxed, onChange: (v) => { relaxed = v; build(); } });
  const out = readout(controls);

  function build() {
    const grid = grids[relaxed], cells = [...grid];
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
    out.set([['worst edge', `cut ${(worst * 100).toFixed(0)}% from its midpoint`], ['typical edge', `${(sum / count * 100).toFixed(1)}%`]]);
    globe.fig.render();
  }

  build();
}
