import { divergence, gradient, curl, cellVector } from '../../js/dynamics/operators.module.js';
import { Figure, choice, legend, caption, readout, text, arrow, anomalyColor, ACCENT, INK, MUTED, LINE } from '../runtime.module.js';

const SQRT3 = Math.sqrt(3), SPACING = 60e3, WALL = SPACING / SQRT3, AREA = SQRT3 / 2 * SPACING ** 2, TRIANGLE = AREA / 2;
const OMEGA = 7.2921e-5, RHO = 1.2, CORE = 100e3, CONVERGE = 30e3, NEGLIGIBLE = 2e-3 * 2 / CONVERGE, WIDTH = 150e3, WIDE = 560, TOP = 26, BOTTOM = 30, HALO = 'rgba(20, 20, 22, 0.6)';
const WINDS = {
  uniform: () => [10, 0],
  converging: (x, y) => [-x / CONVERGE, -y / CONVERGE],
  disk: (x, y) => [-OMEGA * y, OMEGA * x],
  shear: (x, y) => [12 + 5e-5 * y, 0],
  vortex: (x, y) => { const s = 1e-4 * Math.min(1, CORE ** 2 / (x * x + y * y)); return [-s * y, s * x]; },
};
const PRESSURES = {
  slope: (x) => 1e5 - 1e-3 * x,
  low: (x, y) => 1e5 - 400 * Math.exp(-(x * x + y * y) / WIDTH ** 2),
  high: (x, y) => 1e5 + 400 * Math.exp(-(x * x + y * y) / WIDTH ** 2),
};
const COMPASS = ['east', 'east-northeast', 'northeast', 'north-northeast', 'north', 'north-northwest', 'northwest', 'west-northwest', 'west', 'west-southwest', 'southwest', 'south-southwest', 'south', 'south-southeast', 'southeast', 'east-southeast'];
const SUPERSCRIPT = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
const sci = (v) => { const [m, p] = Math.abs(v).toExponential(2).split('e'); return `${v < 0 ? '−' : ''}${m} × 10${[...String(Number(p))].map((c) => SUPERSCRIPT[c]).join('')}`; };
const fixed = (v) => (Math.abs(v) < 0.005 ? '0.00' : v.toFixed(2).replace('-', '−'));
const listed = (values) => `${values.slice(0, -1).map(fixed).join(', ')} and ${fixed(values[values.length - 1])} m/s`;
const compass = (x, y) => COMPASS[((Math.round(Math.atan2(y, x) / (Math.PI / 8)) % 16) + 16) % 16];
const period = (z) => { const hours = 4 * Math.PI / Math.abs(z) / 3600; return hours < 100 ? `${hours.toFixed(1)} hours` : `${Math.round(hours / 24)} days`; };
const turning = (angle) => (angle + 2 * Math.PI + 1e-9) % (2 * Math.PI);

function patch(cols, rows) {
  const c0 = cols >> 1, r0 = rows >> 1, nCells = cols * rows, xy = [];
  for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) xy.push([SPACING * (c - c0 + ((r & 1) - (r0 & 1)) / 2), -1.5 * WALL * (r - r0)]);
  const rings = xy.map(([x, y], i) => xy.flatMap(([u, v], j) => (j !== i && Math.hypot(u - x, v - y) < 1.01 * SPACING ? [j] : [])).sort((a, b) => turning(Math.atan2(xy[a][1] - y, xy[a][0] - x)) - turning(Math.atan2(xy[b][1] - y, xy[b][0] - x))));
  const nEdgesOnCell = new Int32Array(nCells), edgesOnCell = new Int32Array(6 * nCells).fill(-1), edgeSignOnCell = new Int32Array(6 * nCells);
  const pairs = [], mid = [], normal = [], lookup = new Map();
  rings.forEach((ring, i) => {
    nEdgesOnCell[i] = ring.length;
    ring.forEach((j, k) => {
      if (j < i) { edgesOnCell[6 * i + k] = lookup.get(`${j},${i}`); edgeSignOnCell[6 * i + k] = -1; return; }
      const e = pairs.length, n = [(xy[j][0] - xy[i][0]) / SPACING, (xy[j][1] - xy[i][1]) / SPACING];
      pairs.push([i, j]); mid.push([(xy[i][0] + xy[j][0]) / 2, (xy[i][1] + xy[j][1]) / 2]); normal.push(n); lookup.set(`${i},${j}`, e);
      edgesOnCell[6 * i + k] = e; edgeSignOnCell[6 * i + k] = 1;
    });
  });
  const nEdges = pairs.length, corners = new Map();
  mid.forEach(([mx, my], e) => {
    for (const side of [-1, 1]) {
      const p = [mx - side * normal[e][1] * WALL / 2, my + side * normal[e][0] * WALL / 2], key = `${Math.round(p[0])},${Math.round(p[1])}`;
      if (!corners.has(key)) corners.set(key, { p, edges: [] });
      corners.get(key).edges.push([e, side]);
    }
  });
  const vertices = [...corners.values()].filter((v) => v.edges.length === 3).map((v) => ({ ...v, edges: v.edges.sort(([a], [b]) => turning(Math.atan2(mid[a][1] - v.p[1], mid[a][0] - v.p[0])) - turning(Math.atan2(mid[b][1] - v.p[1], mid[b][0] - v.p[0]))) }));
  for (const v of vertices) v.cells = [...new Set(v.edges.flatMap(([e]) => pairs[e]))].sort((a, b) => turning(Math.atan2(xy[a][1] - v.p[1], xy[a][0] - v.p[0])) - turning(Math.atan2(xy[b][1] - v.p[1], xy[b][0] - v.p[0])));
  const mesh = {
    nCells, nEdges, nVertices: vertices.length, maxEdges: 6, nEdgesOnCell, edgesOnCell, edgeSignOnCell,
    cellsOnEdge: Int32Array.from(pairs.flat()), dcEdge: new Float64Array(nEdges).fill(SPACING), dvEdge: new Float64Array(nEdges).fill(WALL),
    nEdge: Float64Array.from(normal.flatMap(([x, y]) => [x, y, 0])), areaCell: new Float64Array(nCells).fill(AREA),
    edgesOnVertex: Int32Array.from(vertices.flatMap((v) => v.edges.map(([e]) => e))), edgeSignOnVertex: Int32Array.from(vertices.flatMap((v) => v.edges.map(([, s]) => s))),
    areaTriangle: new Float64Array(vertices.length).fill(TRIANGLE),
  };
  return { cols, rows, mesh, xy, mid, normal, vertices, whole: Array.from(nEdgesOnCell, (n) => n === 6) };
}

export function mountOperators(root) {
  const controls = root.querySelector('.controls');
  let view = 'gradient', windName = 'converging', pressureName = 'slope', grid = null, cell = 0, vertex = 0, frame = null;
  let u, div, vort, p, g, slope;
  const converging = () => view === 'divergence' && windName === 'converging';
  const cellAt = (x0) => grid.xy.findIndex(([x, y], i) => grid.whole[i] && Math.hypot(x - x0, y) < 1);

  const fig = new Figure(root, { height: (w) => layout(w).h, draw });
  caption(root, `A flat patch of cells the size of the model’s at N = 128, seen from above with north up: centers ${SPACING / 1000} km apart, walls ${(WALL / 1000).toFixed(1)} km long, cells of ${Math.round(AREA / 1e6).toLocaleString('en-US')} km² and triangles of ${Math.round(TRIANGLE / 1e6).toLocaleString('en-US')} km² between three centers. Every number comes from the model’s own code.`);
  const key = legend(root, []);
  choice(controls, { label: 'Show', options: [['the slope of the pressure, ∇p', 'gradient'], ['the spread of the wind, ∇·u', 'divergence'], ['the spin of the wind, ∇×u', 'curl']], value: view, onChange: (v) => { view = v; cell = cellAt(converging() ? 2 * SPACING : 0); if (converging()) cell = cellAt(2 * SPACING); configure(); }, span: true });
  choice(controls, { label: 'Pressure', options: [['a steady slope', 'slope'], ['a low', 'low'], ['a high', 'high']], value: pressureName, onChange: (v) => { pressureName = v; compute(); fig.render(); }, span: true });
  const pressureRow = controls.lastElementChild;
  choice(controls, { label: 'Wind', options: [['a uniform wind', 'uniform'], ['converging', 'converging'], ['a spinning disk', 'disk'], ['a shear', 'shear'], ['a lone vortex', 'vortex']], value: windName, onChange: (v) => { windName = v; if (converging()) cell = cellAt(2 * SPACING); compute(); fig.render(); }, span: true });
  const windRow = controls.lastElementChild;
  const out = readout(controls);

  function layout(w) {
    const wide = w >= WIDE, cols = wide ? 9 : 7, rows = wide ? 7 : 9, R = Math.min(36, (w - 24) / (SQRT3 * (cols + 0.5)));
    return { cols, rows, R, h: Math.round((1.5 * (rows - 1) + 2) * R + TOP + BOTTOM) };
  }

  function rebuild(cols, rows) {
    grid = patch(cols, rows);
    const { mesh, vertices } = grid;
    u = new Float64Array(mesh.nEdges); div = new Float64Array(mesh.nCells); vort = new Float64Array(mesh.nVertices);
    p = new Float64Array(mesh.nCells); g = new Float64Array(mesh.nEdges); slope = new Float64Array(3 * mesh.nCells);
    cell = cellAt(converging() ? 2 * SPACING : 0);
    vertex = vertices.reduce((best, v, k) => (Math.hypot(v.p[0], v.p[1] - WALL) < Math.hypot(vertices[best].p[0], vertices[best].p[1] - WALL) ? k : best), 0);
    compute();
  }

  function compute() {
    const { mesh, xy, mid, normal } = grid, wind = WINDS[windName], pressure = PRESSURES[pressureName];
    for (let e = 0; e < mesh.nEdges; e++) { const [a, b] = wind(...mid[e]); u[e] = a * normal[e][0] + b * normal[e][1]; }
    for (let i = 0; i < mesh.nCells; i++) p[i] = pressure(...xy[i]);
    divergence(mesh, u, div);
    curl(mesh, u, vort);
    gradient(mesh, p, g);
    cellVector(mesh, g, slope);
    show();
  }

  function configure() {
    pressureRow.style.display = view === 'gradient' ? '' : 'none';
    windRow.style.display = view === 'gradient' ? 'none' : '';
    const wind = [['arrow', 'the wind across each edge, the number the model stores', INK], ['arrow', `the winds in the chosen ${view === 'curl' ? 'triangle’s' : 'cell’s'} sum`, ACCENT]];
    key.set(view === 'gradient' ? [['gradient', 'pressure at each cell’s center, lower to higher', 'cool', 'warm'], ['force', 'the push across each edge, toward lower pressure', ACCENT], ['force', 'the push on the chosen cell, from its six edges', '#fff']]
      : view === 'divergence' ? [['gradient', 'each cell filling (converging) to draining (diverging)', 'cool', 'warm'], ...wind]
        : [['gradient', 'each triangle spinning clockwise to counterclockwise', 'cool', 'warm'], ['dash', 'the triangles joining three cell centers', 'rgba(255,255,255,0.5)'], ...wind]);
    show();
    fig.render();
  }

  function show() {
    if (!grid) return;
    const { mesh, vertices } = grid;
    if (view === 'gradient') {
      const gx = slope[3 * cell], gy = slope[3 * cell + 1], size = Math.hypot(gx, gy);
      out.set([['pressure here', `${(p[cell] / 100).toFixed(1)} hPa`], ...(size < 1e-8 ? [['slope', '0: the pressure is level here, and nothing pushes']] : [['slope', `${sci(size)} Pa/m, ${(size * 1000).toFixed(1)} hPa every 100 km`], ['push', `${(size / RHO * 3600).toFixed(1)} m/s an hour toward the ${compass(-gx, -gy)}`]])]);
    } else if (view === 'divergence') {
      const walls = Array.from({ length: 6 }, (_, k) => mesh.edgeSignOnCell[6 * cell + k] * u[mesh.edgesOnCell[6 * cell + k]]), d = div[cell], share = -d * 3600 * 100;
      out.set([['wind out through the six walls', listed(walls)], ['divergence', `${Math.abs(d) < NEGLIGIBLE ? '0' : sci(d)} per second`], ['the cell', Math.abs(share) < 0.05 ? 'keeps its air' : `${share > 0 ? 'gains' : 'loses'} ${Math.abs(share).toFixed(1)}% of its air an hour`]]);
    } else {
      const sides = vertices[vertex].edges.map(([e, s]) => s * u[e]), z = vort[vertex];
      out.set([['wind along the three sides, counterclockwise', listed(sides)], ['curl', `${Math.abs(z) < 1e-9 ? '0' : sci(z)} per second`], ['a paddle wheel here', Math.abs(z) < 1e-9 ? 'would not turn' : `would turn ${z > 0 ? 'counterclockwise' : 'clockwise'} once every ${period(z)}`]]);
    }
  }

  function draw(ctx, w, h) {
    const L = layout(w);
    if (!grid || grid.cols !== L.cols) rebuild(L.cols, L.rows);
    const { mesh, xy, mid, normal, vertices, whole } = grid, scale = L.R / WALL;
    const xs = xy.map(([x]) => x), ys = xy.map(([, y]) => y);
    const cx = (w - (Math.max(...xs) + Math.min(...xs)) * scale) / 2, cy = TOP + L.R + Math.max(...ys) * scale;
    const sx = (x) => cx + x * scale, sy = (y) => cy - y * scale;
    frame = { sx, sy, scale };
    const hex = (i) => { ctx.beginPath(); for (let k = 0; k < 6; k++) { const a = Math.PI / 6 + k * Math.PI / 3; ctx.lineTo(sx(xy[i][0]) + L.R * Math.cos(a), sy(xy[i][1]) + L.R * Math.sin(a)); } ctx.closePath(); };
    const triangle = (v) => { ctx.beginPath(); for (const c of v.cells) ctx.lineTo(sx(xy[c][0]), sy(xy[c][1])); ctx.closePath(); };
    const swing = Math.max(...Array.from(p, (v) => Math.abs(v - 1e5))) || 1;
    for (let i = 0; i < mesh.nCells; i++) {
      hex(i);
      const tone = view === 'gradient' ? anomalyColor((p[i] - 1e5) / swing, 0.75) : view === 'divergence' && whole[i] ? anomalyColor(div[i] / 1e-4, 0.85) : 'rgba(255,255,255,0.03)';
      ctx.fillStyle = tone; ctx.fill();
    }
    if (view === 'curl') vertices.forEach((v, k) => { triangle(v); ctx.fillStyle = anomalyColor(vort[k] / 1.5e-4, 0.85); ctx.fill(); });
    ctx.strokeStyle = LINE; ctx.lineWidth = 1;
    for (let i = 0; i < mesh.nCells; i++) { hex(i); ctx.globalAlpha = whole[i] || view === 'gradient' ? 1 : 0.5; ctx.stroke(); }
    ctx.globalAlpha = 1;
    if (view === 'curl') { ctx.save(); ctx.setLineDash([3, 3]); ctx.strokeStyle = 'rgba(255,255,255,0.5)'; for (const v of vertices) { triangle(v); ctx.stroke(); } ctx.restore(); }
    ctx.fillStyle = '#fff';
    for (const [x, y] of xy) { ctx.beginPath(); ctx.arc(sx(x), sy(y), 2.2, 0, Math.PI * 2); ctx.fill(); }
    ctx.strokeStyle = ACCENT; ctx.lineWidth = 2.5; ctx.lineJoin = 'round';
    if (view === 'curl') triangle(vertices[vertex]); else hex(cell);
    ctx.stroke();
    const reach = 0.85 * SQRT3 * L.R;
    const haloed = (x0, y0, x1, y1, style) => { arrow(ctx, x0, y0, x1, y1, { ...style, color: HALO, width: style.width + 2, head: style.head + 2, dash: null, open: false }); arrow(ctx, x0, y0, x1, y1, style); };
    const across = (e, len, style, halo = true) => { const x = sx(mid[e][0]), y = sy(mid[e][1]), [nx, ny] = normal[e]; (halo ? haloed : arrow.bind(null, ctx))(x - nx * len / 2, y + ny * len / 2, x + nx * len / 2, y - ny * len / 2, style); };
    if (view === 'gradient') {
      const per = reach / Math.max(1e-12, ...Array.from(g, Math.abs));
      for (let e = 0; e < mesh.nEdges; e++) across(e, -g[e] * per, { color: ACCENT, width: 1.5, head: 6, dash: [3, 2.5], open: true }, false);
      const gx = slope[3 * cell] * per, gy = slope[3 * cell + 1] * per, x = sx(xy[cell][0]), y = sy(xy[cell][1]);
      haloed(x, y, x - gx, y + gy, { color: '#fff', width: 2.2, head: 8, dash: [4, 3], open: true });
    } else {
      const per = reach / Math.max(10, ...Array.from(u, Math.abs));
      const chosen = new Set(view === 'curl' ? vertices[vertex].edges.map(([e]) => e) : Array.from({ length: 6 }, (_, k) => mesh.edgesOnCell[6 * cell + k]));
      for (let e = 0; e < mesh.nEdges; e++) if (!chosen.has(e)) across(e, u[e] * per, { color: INK, width: 1.5, head: 6 });
      for (const e of chosen) across(e, u[e] * per, { color: ACCENT, width: 2, head: 6 });
    }
    text(ctx, view === 'curl' ? 'tap a triangle to read its spin' : 'tap a cell to read it', w / 2, 13, { align: 'center', color: MUTED, size: 11 });
    const formula = { gradient: `across each edge: the pressure difference ÷ ${SPACING / 1000} km`, divergence: `each cell: Σ wind out × ${(WALL / 1000).toFixed(1)} km ÷ ${Math.round(AREA / 1e6).toLocaleString('en-US')} km²`, curl: `each triangle: Σ wind along a side × ${SPACING / 1000} km ÷ ${Math.round(TRIANGLE / 1e6).toLocaleString('en-US')} km²` };
    text(ctx, formula[view], w / 2, h - 13, { align: 'center', color: MUTED, size: 11 });
  }

  function pick({ x, y }) {
    if (!frame || !grid) return null;
    const px = (x - frame.sx(0)) / frame.scale, py = (frame.sy(0) - y) / frame.scale;
    if (view === 'curl') {
      let best = null, near = WALL;
      grid.vertices.forEach((v, k) => { const d = Math.hypot(px - v.p[0], py - v.p[1]); if (d < near) { near = d; best = k; } });
      return best;
    }
    let best = null, near = 0.95 * WALL;
    grid.xy.forEach(([cx, cy], i) => { const d = Math.hypot(px - cx, py - cy); if (grid.whole[i] && d < near) { near = d; best = i; } });
    return best;
  }
  function probe(at) {
    const k = pick(at);
    if (k === null) return false;
    if (view === 'curl') vertex = k; else cell = k;
    show();
    fig.render();
  }
  fig.pointer({ hit: (at) => pick(at) !== null, down: probe, move: probe });

  configure();
}
