import { Grid } from '../../js/grid.module.js';
import { Icosahedron } from '../../js/isea.module.js';
import { Globe, arc, rgb } from '../globe.module.js';
import { slider, choice, legend, readout, rampColor, ACCENT } from '../runtime.module.js';

const RADIUS = 6371;

export function mountIcosa(root) {
  const controls = root.querySelector('.controls');
  const globe = new Globe(root, { height: 440 });
  legend(root, [['ramp', 'cells smaller to larger than the average, on the same scale as above', 'cool', 'warm', 'neutral'], ['faint', 'the twelve pentagons', ACCENT], ['faint', 'edges of the icosahedron', 'rgba(255,255,255,0.9)']]);
  let N = 6, overlay = true;
  slider(controls, { label: 'Subdivisions N', min: 2, max: 20, step: 1, value: N, format: (v) => `${v}`, onInput: (v) => { N = v; build(); } });
  choice(controls, { label: 'Icosahedron', options: [['shown', 'on'], ['hidden', 'off']], value: 'on', onChange: (v) => { overlay = v === 'on'; build(); } });
  const out = readout(controls);
  const ico = new Icosahedron();
  const icoEdges = [];
  for (let i = 0; i < ico.vertices.length; i++) for (let j = i + 1; j < ico.vertices.length; j++) {
    const a = ico.vertices[i].clone().normalize(), b = ico.vertices[j].clone().normalize();
    if (a.distanceTo(b) < 1.2) icoEdges.push(arc(a, b, 16));
  }

  function build() {
    const grid = new Grid(N), cells = [...grid], areas = cells.map((c) => c.area), mean = areas.reduce((a, b) => a + b) / areas.length;
    globe.clear();
    globe.cells(cells.map((c) => ({ center: c.centerVertex, vertices: c.vertices })), (cell, k) => rgb(rampColor(0.5 + Math.log2(areas[k] / mean) / 4)));
    const edges = [], pentagons = [];
    for (const cell of cells) {
      const n = cell.vertices.length;
      for (let k = 0; k < n; k++) {
        const a = cell.vertices[k], b = cell.vertices[(k + 1) % n];
        if (cell.isPentagon) pentagons.push([a, b]);
        else if (!cell.neighbors.some((nb) => nb.index < cell.index && nb.vertices.includes(a) && nb.vertices.includes(b))) edges.push([a, b]);
      }
    }
    globe.lines(edges, { color: 0x000000, opacity: 0.4 });
    globe.lines(pentagons, { color: ACCENT, opacity: 1, lift: 1.004 });
    if (overlay) globe.lines(icoEdges, { color: 0xffffff, opacity: 0.9, lift: 1.005 });
    const width = Math.sqrt(4 * Math.PI * RADIUS * RADIUS / cells.length);
    out.set([['cells', `10 × ${N}² + 2 = ${cells.length}`], ['largest cell', `${(Math.max(...areas) / Math.min(...areas)).toFixed(2)}× the smallest`], ['cell width', `about ${width.toFixed(0)} km`]]);
    globe.fig.render();
  }

  build();
}
