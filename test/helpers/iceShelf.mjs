import assert from 'node:assert/strict';
import { Grid } from '../../js/grid.module.js';
import { createModel } from '../../js/model.module.js';
import { gridMesh } from '../../js/stateGrid.module.js';
import { createOcean, displacedWeight } from '../../js/ocean/layered.module.js';
import { syntheticTopography, encodeLandIce, decodeLandIce, ICE_DENSITY, SEAWATER_DENSITY } from '../../js/geography.module.js';

export const N = 16;
const deg = 180 / Math.PI;
export const sparse = ({ rows, cols, points, count, ...fields }) => {
  const index = Uint32Array.from(Array.from(count.keys()).filter((p) => count[p] > 0)), take = (Kind, field) => Kind.from(index, (p) => (field ? field[p] : 0));
  return decodeLandIce(encodeLandIce({
    rows, cols, points, index, count: take(Uint8Array, count), floating: take(Uint8Array, fields.floating), thickness: take(Uint16Array, fields.thickness),
    floatThickness: take(Uint16Array, fields.floatThickness), floatSurface: take(Uint16Array, fields.floatSurface), floatBed: take(Int16Array, fields.floatBed),
  }));
};

export function shelfBands(slope = 0, { water = 0, shelfSea = 0, thickness: H = 400 } = {}) {
  const mesh = gridMesh(new Grid(N)), rows = 90, cols = 180, points = 225, freeboard = 50, draft = H * ICE_DENSITY / SEAWATER_DENSITY, deep = -800, shallow = -Math.round(draft + 30);
  const zone = (lat) => (lat < -70 ? 'grounded' : lat < -54 && !water ? 'shallow' : lat < -38 ? 'deep' : lat < -22 && shelfSea ? 'shelfSea' : 'sea');
  const raster = syntheticTopography(rows, cols, (lat) => ({ grounded: 500, shallow: freeboard, deep: freeboard, shelfSea: -shelfSea, sea: -4000 })[zone(lat * deg)]);
  const pixels = (f) => Array.from({ length: rows * cols }, (_, n) => { const lat = 90 - (Math.floor(n / cols) + 0.5) * 180 / rows; return f(zone(lat), lat); }), afloat = (z) => z === 'shallow' || z === 'deep';
  const shelf = (z, lat) => (z === 'deep' ? Math.round(H + slope * Math.max(-1, Math.min(1, (-46 - lat) / 8))) : H);
  const ice = sparse({
    rows, cols, points, count: pixels((z) => (z === 'sea' || z === 'shelfSea' ? 0 : points)), thickness: pixels((z, lat) => (z === 'grounded' ? 1500 : afloat(z) ? shelf(z, lat) : 0)),
    floating: pixels((z) => (afloat(z) ? points : 0)), floatThickness: pixels((z, lat) => (afloat(z) ? shelf(z, lat) : 0)), floatSurface: pixels((z) => (afloat(z) ? freeboard : 0)),
    floatBed: pixels((z, lat) => (z === 'deep' ? (water ? -Math.round(shelf(z, lat) * ICE_DENSITY / SEAWATER_DENSITY + water) : deep) : z === 'shallow' ? shallow : 0)),
  });
  return { mesh, raster, zone, draft, deep, options: { landBridges: {}, seaStraits: {}, landIce: ice } };
}

export function warmOcean(mesh, geography, options = {}) {
  const C = mesh.nCells, ocean = createOcean(mesh, { geography, everySteps: 1, ...options }), L = ocean.layers, surfaceT = new Float64Array(C).fill(275), ice = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const saved = { ...ocean.serialize(), ...ocean.restartArrays(), flux: new Float64Array(C) }, isomip = (z) => [271.25 + 2.9 * z / 720, 33.8 + 0.9 * z / 720];
  for (let i = 0; i < C; i++) {
    if (!geography.cavity[i]) continue;
    let z = geography.iceDraft[i];
    for (let k = 0; k < L; k++) { const n = k * C + i, hk = ocean.h[n], [t, sal] = isomip(z + 0.5 * hk); saved.Q[n] = hk * t; saved.W[n] = hk * sal; z += hk; }
    saved.previousT0[i] = saved.Q[i] / ocean.h[i];
    saved.shelfLoad[i] = displacedWeight(geography.iceDraft[i], isomip);
  }
  ocean.load(saved, surfaceT, ice);
  return { ocean, surfaceT, ice };
}

async function restMonth({ raster, options }, iceShelves, profile, engines) {
  const setup = { topography: raster, geography: { ...options, iceShelves }, ocean: { iceDrift: false } }, { mesh, geography } = createModel(new Grid(N), setup), C = mesh.nCells, E = mesh.nEdges;
  const cavity = geography.cavity ?? new Uint8Array(C), cpu = createOcean(mesh, { geography, everySteps: 1, shelfMelt: false }), L = cpu.layers, ice = new Float64Array(C), calm = new Float64Array(E), flux = new Float64Array(C);
  const [surface] = profile(0);
  cpu.initialize(new Float64Array(C).fill(surface), ice);
  const saved = { ...cpu.serialize(), ...cpu.restartArrays(), flux: new Float64Array(C) };
  for (let i = 0; i < C; i++) {
    let sum = 0, z = cavity[i] ? geography.iceDraft[i] : 0;
    for (let k = 0; k < L; k++) sum += saved.h[k * C + i];
    for (let k = 0; k < L; k++) {
      const n = k * C + i;
      if (sum > 0) saved.h[n] *= cpu.D[i] / sum;
      const [t, s] = profile(z + 0.5 * saved.h[n]);
      z += saved.h[n];
      saved.Q[n] = saved.h[n] * t; saved.W[n] = saved.h[n] * s; saved.T[n] = t; saved.S[n] = s;
    }
    saved.eta[i] = 0; saved.previousT0[i] = saved.T[i];
    if (saved.shelfLoad) saved.shelfLoad[i] = cavity[i] ? displacedWeight(geography.iceDraft[i], profile) : 0;
  }
  saved.u.fill(0);
  const onCpu = engines.includes('cpu'), onGpu = engines.includes('gpu'), createGpuModel = onGpu ? (await import('../../js/gpu/model.gpu.js')).createGpuModel : null;
  const gpu = onGpu ? (await createGpuModel(new Grid(N), { ...setup, ocean: { everySteps: 1, shelfMelt: false, iceDrift: false } })).oceanEngine : null, cpuSkin = new Float64Array(C).fill(surface), gpuSkin = new Float64Array(C).fill(surface);
  cpu.load(saved, cpuSkin, ice);
  if (onGpu) gpu.upload(saved, gpuSkin, ice);
  for (let n = 0; n < 720; n++) { if (onCpu) cpu.advance(cpuSkin, ice, flux, calm, 3600); if (onGpu) await gpu.advance(gpuSkin, ice, calm, 3600); }
  const b = onGpu ? await gpu.download() : null, most = (a) => a.reduce((m, v) => Math.max(m, Math.abs(v)), 0), under = (u) => {
    let m = 0;
    for (let e = 0; e < E; e++) if (cavity[mesh.cellsOnEdge[2 * e]] || cavity[mesh.cellsOnEdge[2 * e + 1]]) for (let k = 0; k < L; k++) m = Math.max(m, Math.abs(u[k * E + e]));
    return m;
  };
  return { ...(onCpu ? { cpu: most(cpu.u), cpuEta: most(cpu.eta), cpuCavity: under(cpu.u) } : {}), ...(onGpu ? { gpu: most(b.u), gpuEta: most(b.eta), gpuCavity: under(b.u) } : {}), ocean: cpu, mesh, geography };
}
export const show = (r) => ['cpu', 'gpu', 'cpuEta', 'gpuEta'].filter((k) => k in r).map((k) => `${k} ${r[k].toExponential(1)}`).join(', ');

export async function restingShelf(engines) {
  const bands = shelfBands(250), still = () => [271.65, 34.4], shelves = await restMonth(bands, 'cavity', still, engines), bare = await restMonth(bands, 'land', still, engines);
  console.log(`a month at rest, homogeneous water: with the shelves ${show(shelves)}; with them as land ${show(bare)}`);
  // Rounding walks both oceans off rest; 10× leaves room for the cavity's added cells and catches a load, draft or pin that pushes the water.
  for (const [u, eta] of engines.map((e) => [e, `${e}Eta`])) assert.ok(shelves[u] < 10 * bare[u] && shelves[eta] < 10 * bare[eta], `${u}: ${show(shelves)} against ${show(bare)}`);
}

export async function steppedDraft(engines) {
  const bands = shelfBands(700, { water: 200, thickness: 760 }), still = () => [271.65, 34.4], isomip = (z) => [271.25 + 2.9 * Math.min(z, 720) / 720, 33.8 + 0.9 * Math.min(z, 720) / 720];
  const shelves = await restMonth(bands, 'cavity', still, engines), bare = await restMonth(bands, 'land', still, engines), { ocean, mesh, geography: { cavity, iceDraft } } = shelves, floor = 50, top = (i) => (cavity[i] ? iceDraft[i] : 0);
  let step = 0, empty = 0, thin = 0;
  for (let e = 0; e < mesh.nEdges; e++) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
    if (!ocean.cellOcean[a] || !ocean.cellOcean[b]) continue;
    const shared = Math.min(top(a) + ocean.D[a], top(b) + ocean.D[b]) - Math.max(top(a), top(b));
    if (cavity[a] && cavity[b]) step = Math.max(step, Math.abs(top(a) - top(b)));
    if (shared <= 0) empty++;
    else if (shared < floor) thin++;
    assert.ok(ocean.sillEdge[e] >= floor, `edge ${e} passes ${ocean.sillEdge[e]} m`);
  }
  console.log(`drafts step by up to ${step.toFixed(0)} m between cavity cells: ${empty} edges share no depth and ${thin} share less than ${floor} m; each passes ${floor} m`);
  assert.ok(step >= 300 && empty > 0, `drafts step by ${step} m, ${empty} edges share no depth`);
  console.log(`a month at rest, homogeneous water: with the shelves ${show(shelves)}; with them as land ${show(bare)}`);
  for (const [u, eta] of engines.map((e) => [e, `${e}Eta`])) assert.ok(shelves[u] < 10 * bare[u] && shelves[eta] < 10 * bare[eta], `${u}: ${show(shelves)} against ${show(bare)}`);
  const layered = await restMonth(bands, 'cavity', isomip, engines), flat = await restMonth(bands, 'land', isomip, engines);
  console.log(`a month at rest on ISOMIP+'s warm profile: with the shelves ${show(layered)}, on the cavities' edges ${engines.map((e) => `${e} ${layered[`${e}Cavity`].toExponential(1)}`).join(', ')}; with them as land ${show(flat)}`);
  // On a stratified profile the sloping base's pressure gradient moves the cavity edges at about twice the open sea's speed; 3× still catches a runaway.
  for (const [u, eta, under] of engines.map((e) => [e, `${e}Eta`, `${e}Cavity`])) assert.ok(layered[u] < 10 * flat[u] && layered[eta] < 10 * flat[eta] && layered[under] < 3 * flat[u], `${u}: ${show(layered)}, cavities ${layered[under]}, against ${show(flat)}`);
}
