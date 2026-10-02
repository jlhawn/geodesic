// The day-mean shortwave and longwave heating of the tropical boxes under two
// radiation configurations on the same held state, on the single-thread CPU
// engine with the ocean off:
//   node scripts/radiationProfiles.mjs <state.bin>
// One model step is taken under RADIATION (JSON, the defaults if unset), and
// the stepped state, cumulus, boundary layer and deck are copied into a second
// model under OTHER (JSON, '{"longwaveScheme":"gray","solarGases":"lacisHansen"}'
// if unset); every box column's radiation then runs in both at TIMES (24)
// instants of the day (dt 0, so the deck's carried state does not move). The
// layer heating is the radiation's own: its longwave record, and the rest of
// its layer flux less the lowest layer's sensible heat as the shortwave.
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS } from '../js/physics/regrid.module.js';
import { SEA_DRAG, LAND_DRAG } from '../js/physics/surface.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';
import { BOXES, inLongitudes } from '../js/audit.module.js';

const FILE = process.argv[2];
if (!FILE) { console.error('usage: node scripts/radiationProfiles.mjs <state.bin>'); process.exit(1); }
const RADIATION = JSON.parse(process.env.RADIATION ?? '{}'), OTHER = JSON.parse(process.env.OTHER ?? '{"longwaveScheme":"gray","solarGases":"lacisHansen"}'), TIMES = Number(process.env.TIMES ?? 24);
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
function make(radiation) {
  const model = createModel(new Grid(saved.N), { topography, levels: savedLevels(saved), ocean: false, radiation });
  STATE_NAMES.forEach((name, a) => model.state[a].set(saved[name]));
  model.seaIce.load(model.state[6], saved.concentration ?? null);
  for (const field of Object.keys(DECK_FIELDS)) model.radiation[field].set(savedDeckField(saved, field, model));
  model.land.load(saved.land, model.state[6]);
  model.time = saved.time;
  const { g } = model.core.diagnostics, zs = Float64Array.from({ length: model.mesh.nCells }, (_, i) => (model.surfaceGeopotential ? model.surfaceGeopotential[i] / g : 0));
  if (saved.boundaryDepth) model.boundaryLayer.depth.set(Float64Array.from(saved.boundaryDepth, (z, i) => z + zs[i]));
  if (saved.mixingTop) model.boundaryLayer.mixingTop.set(Float64Array.from(saved.mixingTop, (z, i) => z + zs[i]));
  if (saved.boundaryRegime) model.boundaryLayer.regime.set(saved.boundaryRegime);
  if (saved.boundaryBuoyancy) model.boundaryLayer.buoyancyFlux.set(saved.boundaryBuoyancy);
  return model;
}
const A = make(RADIATION), B = make({ ...RADIATION, ...OTHER });
A.step(1350 * 16 / saved.N);
STATE_NAMES.forEach((_, a) => B.state[a].set(A.state[a]));
B.moist.cumulusCover.set(A.moist.cumulusCover); B.moist.cumulusWater.set(A.moist.cumulusWater);
for (const name of ['depth', 'mixingTop', 'regime', 'buoyancyFlux']) B.boundaryLayer[name].set(A.boundaryLayer[name]);
for (const field of Object.keys(DECK_FIELDS)) B.radiation[field].set(A.radiation[field]);
const { K, C, sigmaMid, dSigma, g, cp, geopotential } = A.core.diagnostics;
const deg = 180 / Math.PI, lat = Float64Array.from(A.mesh.latCell, (x) => x * deg), lon = Float64Array.from(A.mesh.lonCell, (x) => x * deg), land = A.geography.land;
const boxes = [
  ['Pacific ITCZ 5-12N 160E-100W', (i) => lat[i] >= BOXES.itcz[0] && lat[i] <= BOXES.itcz[1] && inLongitudes(lon[i], BOXES.itcz[2], BOXES.itcz[3])],
  ['warm pool 10S-10N 120-170E sea', (i) => !land[i] && Math.abs(lat[i]) <= 10 && inLongitudes(lon[i], 120, 170)],
  ['15S-15N', (i) => Math.abs(lat[i]) <= 15],
];
const acc = boxes.map(() => ({ area: 0, pressure: new Float64Array(K), shortwave: [new Float64Array(K), new Float64Array(K)], longwave: [new Float64Array(K), new Float64Array(K)] }));
for (const M of [A, B]) { for (let i = 0; i < C; i++) M.core.diagnoseColumn(i, M.state[0], M.state[1], M.state[4], M.state[5]); M.surface.lowestWindSpeed(M.state[2]); }
for (let n = 0; n < TIMES; n++) {
  [A, B].forEach((M, m) => {
    M.radiation.setTime(A.time + n * 86400 / TIMES);
    const [pi, theta, , surfaceT, q, qc, ice] = M.state;
    for (let i = 0; i < C; i++) {
      const hit = boxes.map(([, keep]) => keep(i));
      if (!hit.some((x) => x)) continue;
      let skin = surfaceT[i], openSea = 0;
      if (!land[i]) { const cover = M.seaIce.cover(i, ice[i]); openSea = 1 - cover; if (ice[i] > 0 && cover < 1) skin = cover * surfaceT[i] + (1 - cover) * FREEZING_POINT; }
      const albedo = M.surfaceAlbedo[i], bottom = (K - 1) * C + i;
      M.radiation.column(i, pi[i], theta, skin, M.surface.windSpeed[i], undefined, M.radiation.insolation(i), q[bottom], q, qc, albedo, albedo, 1, land[i] ? LAND_DRAG : SEA_DRAG, openSea, M.boundaryLayer.depth[i] - geopotential[bottom] / g, 0, M.boundaryLayer.mixingTop[i] - geopotential[bottom] / g);
      const flux = M.radiation.layerFlux, sensible = M.radiation.budget.sensibleHeat;
      hit.forEach((inside, b) => {
        if (!inside) return;
        const S = acc[b], a = A.mesh.areaCell[i];
        if (m === 0 && n === 0) S.area += a;
        for (let k = 0; k < K; k++) {
          const perDay = 86400 / (cp * pi[i] * dSigma[k] / g), longwave = M.radiation.longwave[k * C + i];
          S.longwave[m][k] += a * longwave * perDay / TIMES;
          S.shortwave[m][k] += a * (flux[k] - longwave - (k === K - 1 ? sensible : 0)) * perDay / TIMES;
          if (m === 0 && n === 0) S.pressure[k] += a * pi[i] * sigmaMid[k];
        }
      });
    }
  });
}
console.log(`radiation profiles of ${FILE.split('/').pop()}: day ${saved.day}, N=${saved.N}; first ${JSON.stringify(RADIATION)}, second ${JSON.stringify({ ...RADIATION, ...OTHER })}; ${TIMES} instants`);
boxes.forEach(([name], b) => {
  const S = acc[b];
  console.log(`${name}, K/day: hPa, shortwave first and second, longwave first and second, the first's net less the second's`);
  for (let k = 0; k < K; k++) {
    const p = S.pressure[k] / S.area / 100;
    if (p < 100) continue;
    const v = [S.shortwave[0][k], S.shortwave[1][k], S.longwave[0][k], S.longwave[1][k]].map((x) => x / S.area);
    console.log(`  ${p.toFixed(0).padStart(5)} ${v.map((x) => x.toFixed(2).padStart(7)).join(' ')} ${(v[0] + v[2] - v[1] - v[3]).toFixed(2).padStart(7)}`);
  }
});
