import { Grid } from '../../js/grid.module.js';
import { buildMesh } from '../../js/mesh.module.js';
import { tangential } from '../../js/dynamics/operators.module.js';
import { Figure, slider, legend, readout, text, arrow, ACCENT, INK, MUTED, LINE } from '../runtime.module.js';

export function mountTrisk(root) {
  const controls = root.querySelector('.controls');
  const mesh = buildMesh(new Grid(8), { radius: 1, omega: 0 });
  const { nEdges, maxEdges, maxEdgesOnEdge, cellsOnEdge, nEdgesOnCell, edgesOnCell, verticesOnCell, verticesOnEdge, xEdge, nEdge, tEdge, xVertex, xCell, edgesOnEdge, weightsOnEdge, dvEdge, dcEdge, latEdge } = mesh;
  const at = (array, i) => [array[3 * i], array[3 * i + 1], array[3 * i + 2]];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const unit = (a) => { const n = Math.hypot(...a); return [a[0] / n, a[1] / n, a[2] / n]; };
  let e0 = 0, best = Infinity;
  for (let e = 0; e < nEdges; e++) {
    const [i, j] = [cellsOnEdge[2 * e], cellsOnEdge[2 * e + 1]];
    if (nEdgesOnCell[i] !== 6 || nEdgesOnCell[j] !== 6) continue;
    const p = at(xEdge, e), score = Math.abs(latEdge[e]) + Math.abs(Math.atan2(p[1], p[0]));
    if (score < best) { best = score; e0 = e; }
  }
  const center = at(xEdge, e0), east = unit(cross([0, 0, 1], center)), north = cross(center, east);
  const project = (p) => [dot(p, east), dot(p, north)];
  const cells = [cellsOnEdge[2 * e0], cellsOnEdge[2 * e0 + 1]];
  const ring = [...new Set(cells.flatMap((c) => Array.from({ length: nEdgesOnCell[c] }, (_, k) => edgesOnCell[maxEdges * c + k])))];
  const coefficient = new Map();
  for (let s = 0; s < mesh.nEdgesOnEdge[e0]; s++) coefficient.set(edgesOnEdge[maxEdgesOnEdge * e0 + s], weightsOnEdge[maxEdgesOnEdge * e0 + s] * dvEdge[edgesOnEdge[maxEdgesOnEdge * e0 + s]] / dcEdge[e0]);
  const SPEED = 10;
  let direction = 35, un = new Float64Array(nEdges), ut = 0, exact = 0;

  const fig = new Figure(root, { height: 420, draw });
  legend(root, [['arrow', 'wind across each edge, the number the model stores', INK], ['arrow', 'wind along the shared edge, reconstructed from the ten', ACCENT], ['faint', 'the true wind', 'rgba(255,255,255,0.9)']]);
  slider(controls, { label: 'Wind direction', min: 0, max: 360, step: 1, value: direction, format: (v) => `${v}°`, onInput: (v) => { direction = v; compute(); fig.render(); } });
  const out = readout(controls);

  function velocity(p) {
    const theta = direction * Math.PI / 180, dir = [0, 1, 2].map((k) => Math.cos(theta) * east[k] + Math.sin(theta) * north[k]);
    const axis = cross(center, dir).map((v) => v * SPEED);
    return cross(axis, p);
  }

  function compute() {
    for (let e = 0; e < nEdges; e++) un[e] = dot(velocity(at(xEdge, e)), at(nEdge, e));
    ut = tangential(mesh, un)[e0];
    exact = dot(velocity(center), at(tEdge, e0));
    out.set([['stored across the shared edge', `${un[e0].toFixed(2)} m/s`], ['reconstructed along it', `${ut.toFixed(2)} m/s`], ['true along-edge wind', `${exact.toFixed(2)} m/s`], ['error', `${Math.abs(ut - exact).toFixed(2)} m/s`]]);
  }

  function draw(ctx, w, h) {
    const points = cells.flatMap((c) => Array.from({ length: nEdgesOnCell[c] }, (_, k) => project(at(xVertex, verticesOnCell[maxEdges * c + k]))));
    const xs = points.map((p) => p[0]), ys = points.map((p) => p[1]);
    const span = Math.max(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys)), scale = Math.min(w, h) * 0.78 / span;
    const cx = (Math.max(...xs) + Math.min(...xs)) / 2, cy = (Math.max(...ys) + Math.min(...ys)) / 2;
    const screen = (p) => [w / 2 + (p[0] - cx) * scale, h / 2 - (p[1] - cy) * scale];
    for (const c of cells) {
      ctx.beginPath();
      for (let k = 0; k < nEdgesOnCell[c]; k++) { const [x, y] = screen(project(at(xVertex, verticesOnCell[maxEdges * c + k]))); if (k === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y); }
      ctx.closePath(); ctx.fillStyle = 'rgba(255,255,255,0.05)'; ctx.fill(); ctx.strokeStyle = LINE; ctx.lineWidth = 1; ctx.stroke();
      const [x, y] = screen(project(at(xCell, c))); ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(x, y, 2.5, 0, Math.PI * 2); ctx.fill();
    }
    const perMs = Math.min(w, h) * 0.022;
    for (const e of ring) {
      const [x, y] = screen(project(at(xEdge, e))), n = project(at(nEdge, e)), len = un[e] * perMs;
      if (e === e0) continue;
      arrow(ctx, x - n[0] * len / 2, y + n[1] * len / 2, x + n[0] * len / 2, y - n[1] * len / 2, { color: INK, width: 1.5, head: 6 });
      const k = coefficient.get(e) ?? 0, t = project(at(tEdge, e));
      text(ctx, `×${k.toFixed(2)}`, x + t[0] * 14 - n[0] * 12, y - t[1] * 14 + n[1] * 12, { align: 'center', color: MUTED, size: 10 });
    }
    const [x, y] = screen(project(center)), n = project(at(nEdge, e0)), t = project(at(tEdge, e0)), v = project(velocity(center));
    const ends = [verticesOnEdge[2 * e0], verticesOnEdge[2 * e0 + 1]].map((vertex) => screen(project(at(xVertex, vertex)))).sort((a, b) => a[0] - b[0]);
    ctx.strokeStyle = LINE; ctx.lineWidth = 3.5; ctx.lineCap = 'round'; ctx.beginPath(); ctx.moveTo(...ends[0]); ctx.lineTo(...ends[1]); ctx.stroke(); ctx.lineCap = 'butt';
    const away = [ends[1][0] - ends[0][0], ends[1][1] - ends[0][1]], span2 = Math.hypot(...away);
    text(ctx, 'the shared edge', ends[1][0] + away[0] / span2 * 14, ends[1][1] + away[1] / span2 * 14, { align: 'left', color: INK, size: 12, weight: 700 });
    arrow(ctx, x - n[0] * un[e0] * perMs / 2, y + n[1] * un[e0] * perMs / 2, x + n[0] * un[e0] * perMs / 2, y - n[1] * un[e0] * perMs / 2, { color: INK, width: 2, head: 7 });
    arrow(ctx, x, y, x + t[0] * ut * perMs, y - t[1] * ut * perMs, { color: ACCENT, width: 2.5, head: 8 });
    arrow(ctx, x, y, x + v[0] * perMs, y - v[1] * perMs, { color: 'rgba(255,255,255,0.9)', width: 1, head: 6 });
    text(ctx, 'each number is that edge’s weight in the sum', w / 2, h - 12, { align: 'center', color: MUTED, size: 11 });
  }

  compute();
}
