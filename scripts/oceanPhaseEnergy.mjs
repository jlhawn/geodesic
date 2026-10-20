// The hybrid ocean's kinetic energy changed by each phase of its step, unforced
// (no stress, no surface flux), from a saved coupled state, on the CPU:
//   FROM=<state.bin> STEPS=32 STEP=1350 OCEAN='{"edgeTransport":"remapped"}' node scripts/oceanPhaseEnergy.mjs
// BOX (la0,la1,lo0,lo1 in degrees, 0,12,125,150 by default) is the region
// reported beside the globe, with its barotropic part; EVERY (8) steps
// between lines; WORKERS (12) worker threads. Before each rescale it also
// reports how far the box's columns sum from the barotropic solver's free
// surface, the water the rescale spreads through each column.
import { readFileSync } from 'node:fs';
const R = new URL('../', import.meta.url).pathname;
const [{ decodeState }, { stateGrid }, { topographyFromInt16 }, { createParallelModel }] = await Promise.all([
  import(R + 'js/stateFile.module.js'), import(R + 'js/stateGrid.module.js'), import(R + 'js/geography.module.js'), import(R + 'js/parallel.module.js')]);
const FROM = process.env.FROM, STEPS = Number(process.env.STEPS ?? 32), STEP = Number(process.env.STEP ?? 1350), WORKERS = Number(process.env.WORKERS ?? 12), EVERY = Number(process.env.EVERY ?? 8);
const OCEAN = JSON.parse(process.env.OCEAN ?? '{}');
const [la0, la1, lo0, lo1] = (process.env.BOX ?? '0,12,125,150').split(',').map(Number);
const topography = topographyFromInt16(readFileSync(R + 'data/topography_0p25.bin').buffer);
const saved = await decodeState(new Uint8Array(readFileSync(FROM)));
const mesh = stateGrid(saved, topography).mesh, C = mesh.nCells, E = mesh.nEdges, deg = 180 / Math.PI;
const parallel = await createParallelModel(mesh, { topography, ocean: { vertical: 'hybrid', everySteps: 1, ...OCEAN } }, WORKERS);
const ocean = parallel.ocean, L = ocean.layers, { h, u } = ocean, avgEta = new Float64Array(ocean.shared.baroAvgEta);
const ice = Float64Array.from(saved.ice), concentration = saved.concentration ? Float64Array.from(saved.concentration) : null;
ocean.load(saved.ocean, Float64Array.from(saved.surfaceT), ice, concentration);
ocean.setStress(new Float64Array(E), ice, concentration);
const inside = (lat, lon) => lat >= la0 && lat <= la1 && ((lon % 360) + 360) % 360 >= lo0 && ((lon % 360) + 360) % 360 <= lo1;
const boxEdge = Uint8Array.from({ length: E }, (_, e) => { const [x, y, z] = mesh.xEdge.subarray(3 * e, 3 * e + 3); return inside(Math.asin(z / Math.hypot(x, y, z)) * deg, Math.atan2(y, x) * deg); });
const boxCell = Uint8Array.from({ length: C }, (_, i) => inside(mesh.latCell[i] * deg, mesh.lonCell[i] * deg));

// PJ, with the cells' mean thickness at each edge: box, globe, the box's barotropic part.
function energy() {
  let box = 0, all = 0, barotropic = 0;
  for (let e = 0; e < E; e++) {
    if (!ocean.edgeOcean[e]) continue;
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1], w = 0.5 * 1025 * mesh.dcEdge[e] * mesh.dvEdge[e];
    let sum = 0, transport = 0, depth = 0;
    for (let k = 0; k < L; k++) { const he = 0.5 * (h[k * C + a] + h[k * C + b]), v = u[k * E + e]; sum += he * v * v; transport += he * v; depth += he; }
    all += w * sum;
    if (boxEdge[e]) { box += w * sum; if (depth > 0) barotropic += w * transport * transport / depth; }
  }
  return [box / 1e15, all / 1e15, barotropic / 1e15];
}
const mismatch = { squares: 0, n: 0, largest: 0 };
function columnMismatch() {
  for (let i = 0; i < C; i++) {
    if (!ocean.cellOcean[i] || !boxCell[i]) continue;
    let sum = 0;
    for (let k = 0; k < L; k++) sum += h[k * C + i];
    const m = ocean.D[i] + avgEta[i] - sum;
    mismatch.squares += m * m; mismatch.n++; mismatch.largest = Math.max(mismatch.largest, Math.abs(m));
  }
}
const P = ocean.phases;
const phases = [['dynamics', () => P.dynamics(STEP)], ['rescale', () => { columnMismatch(); P.rescale(); }], ['eddy', () => P.eddy(STEP)], ['shift', () => P.shift()], ['regrid', () => P.regrid(STEP)], ['column', () => P.column(STEP)]];
const sums = Object.fromEntries(phases.map(([name]) => [name, [0, 0, 0]]));
console.log(`kinetic energy by phase from ${FROM} (day ${saved.day}), unforced, step ${STEP} s, ocean ${JSON.stringify(OCEAN)}, box ${la0}..${la1}N ${lo0}..${lo1}E; PJ as box/globe/box barotropic`);
let now = energy();
console.log(`start: box ${now[0].toFixed(3)} (barotropic ${now[2].toFixed(3)}), globe ${now[1].toFixed(3)}`);
for (let n = 1; n <= STEPS; n++) {
  for (const [name, run] of phases) {
    run();
    const after = energy();
    for (let j = 0; j < 3; j++) sums[name][j] += after[j] - now[j];
    now = after;
  }
  if (n % EVERY && n !== STEPS) continue;
  console.log(`after ${n} steps (${(n * STEP / 3600).toFixed(1)} h): box ${now[0].toFixed(3)} (barotropic ${now[2].toFixed(3)}), globe ${now[1].toFixed(3)}; summed by phase ${phases.map(([name]) => `${name} ${sums[name].map((v) => v.toFixed(3)).join('/')}`).join(', ')}`);
  console.log(`  columns against the free surface before the rescale, box: rms ${Math.sqrt(mismatch.squares / Math.max(1, mismatch.n)).toFixed(3)} m, largest ${mismatch.largest.toFixed(3)} m`);
}
await parallel.close();
process.exit(0);
