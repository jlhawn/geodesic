import { Grid, DEFAULT_RELAX } from './grid.module.js';
import { buildMesh } from './mesh.module.js';
import { createGeography, readTopography, savedGeography, TOPOGRAPHIES } from './geography.module.js';
import { SIDEREAL_DAY } from './model.module.js';

/*
 * The grid a saved state was made on. A state records the Lloyd passes of
 * its mesh in its header as `grid: { relax }`. One saved without them is
 * tried on each of STATE_RELAX in turn, the first whose land mask under
 * `topography` matches the saved ocean's dry columns exactly winning, else
 * the one with the fewest mismatches. A state with no ocean, or no dry
 * column, or nothing to match against, is on the default grid. Under node
 * an undefined `topography` is the repository's raster of the state's
 * recorded geography; null matches nothing; a named raster other than the
 * state's is refused. The candidates' land masks are kept per topography,
 * so a second state of the same size is matched without building a grid.
 */
export const STATE_RELAX = [DEFAULT_RELAX, 8];

export function gridHeader(meshOrGrid) {
  const relax = Number.isInteger(meshOrGrid?.relaxPasses) ? meshOrGrid.relaxPasses : meshOrGrid?.relax;
  return Number.isInteger(relax) ? { grid: { relax } } : {};
}

export function savedRelax(saved) {
  const relax = saved?.grid?.relax;
  return Number.isInteger(relax) && relax >= 0 ? relax : null;
}

export function savedDryColumns(saved) {
  const h = saved?.ocean?.h, N = saved?.N;
  if (!h || !N) return null;
  const C = 10 * N * N + 2, L = h.length / C;
  if (!(Number.isInteger(L) && L >= 1)) return null;
  const dry = new Uint8Array(C);
  for (let i = 0; i < C; i++) {
    let sum = 0;
    for (let k = 0; k < L; k++) sum += h[k * C + i];
    dry[i] = sum > 0 ? 0 : 1;
  }
  return dry;
}

export function maskMismatches(dry, land) {
  let count = 0;
  for (let i = 0; i < dry.length; i++) if (!dry[i] !== !land[i]) count++;
  return count;
}

const repository = new Map();
function repositoryTopography(name) {
  if (repository.has(name)) return repository.get(name);
  const fs = globalThis.process?.getBuiltinModule?.('node:fs');
  const url = new URL(`../data/${TOPOGRAPHIES[name]}`, import.meta.url);
  repository.set(name, fs && fs.existsSync(url) ? readTopography(name) : null);
  return repository.get(name);
}

export const gridMesh = (grid) => buildMesh(grid, { omega: 2 * Math.PI / SIDEREAL_DAY });
function build(N, relax) {
  const grid = new Grid(N, { relax });
  return { grid, mesh: gridMesh(grid) };
}

const masks = new WeakMap();
function landMask(N, relax, topography, geography, built) {
  let kept = masks.get(topography);
  if (!kept) masks.set(topography, (kept = new Map()));
  const key = `${N} ${relax} ${JSON.stringify(geography)}`;
  if (!kept.has(key)) {
    const made = build(N, relax);
    built.set(relax, made);
    kept.set(key, Uint8Array.from(createGeography(made.mesh, topography, geography).land, (l) => (l ? 1 : 0)));
  }
  return kept.get(key);
}

function find(saved, topography, { candidates = STATE_RELAX, geography = {} }, built) {
  const named = savedGeography(saved).topography;
  if (topography?.name && topography.name !== named) throw new Error(`the state is on the ${named} topography, not ${topography.name}: read its geography with stateGeography`);
  const recorded = savedRelax(saved);
  if (recorded !== null) return { relax: recorded, source: 'header', mismatches: null };
  const dry = savedDryColumns(saved), land = topography === undefined ? repositoryTopography(named) : topography;
  if (!dry || !land || !dry.some(Boolean)) return { relax: candidates[0], source: 'default', mismatches: null };
  let best = null;
  for (const relax of candidates) {
    const mismatches = maskMismatches(dry, landMask(saved.N, relax, land, geography, built));
    if (!best || mismatches < best.mismatches) best = { relax, source: 'mask', mismatches };
    if (mismatches === 0) break;
  }
  return best;
}

export function stateRelax(saved, topography = undefined, options = {}) {
  return find(saved, topography, options, new Map());
}

export function stateGrid(saved, topography = undefined, options = {}) {
  const built = new Map(), found = find(saved, topography, options, built);
  return { ...(built.get(found.relax) ?? build(saved.N, found.relax)), ...found };
}
