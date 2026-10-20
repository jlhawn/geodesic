// The model's own CPU radiation on one saved state at HOURS (24) times of
// the state's day, the cloud held as it is after one step (ocean off):
//   RADIATION='{...}' node scripts/dayRadiation.mjs <state.bin>
// Prints the day means of the absorbed sunlight, the OLR, the cloud
// effects (40S-40N apart for the shortwave's) and the clear-sky OLR, area
// weighted over the globe, so that a radiation option's effect on one
// climate state is read without the cloud's response to it.
import { readFileSync } from 'node:fs';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS } from '../js/physics/regrid.module.js';
import { stateGrid } from '../js/stateGrid.module.js';

const FILE = process.argv[2];
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(stateGrid(saved, topography).mesh, { topography, levels: savedLevels(saved), ocean: false, radiation: { clearSkyPass: true, ...JSON.parse(process.env.RADIATION ?? '{}') } });
const { mesh, core, state, radiation, boundaryLayer: bl, seaIce, land, geography } = model;
const { K, C, g } = core.diagnostics;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
for (const field of Object.keys(DECK_FIELDS)) radiation[field].set(savedDeckField(saved, field, model));
land.load({ soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow), ...(saved.land.vegetation ? { vegetation: Float64Array.from(saved.land.vegetation) } : {}), ...(saved.land.surface ? { surface: Float64Array.from(saved.land.surface) } : {}) });
model.time = saved.time;
const zs = Float64Array.from({ length: C }, (_, i) => (model.surfaceGeopotential ? model.surfaceGeopotential[i] / g : 0));
if (saved.boundaryDepth) bl.depth.set(Float64Array.from(saved.boundaryDepth, (z, i) => z + zs[i]));
if (saved.mixingTop) bl.mixingTop.set(Float64Array.from(saved.mixingTop, (z, i) => z + zs[i]));
if (saved.boundaryRegime) bl.regime.set(saved.boundaryRegime);
if (saved.boundaryBuoyancy) bl.buoyancyFlux.set(saved.boundaryBuoyancy);
model.step(1350 * 16 / saved.N);
const [pi, theta, , surfaceT, q, qc, ice] = state;
const landMask = geography.land, area = mesh.areaCell, deg = 180 / Math.PI;
const fluxT = new Float64Array(C), adir = new Float64Array(C), adif = new Float64Array(C), openSea = new Float64Array(C);
const fluxState = [pi, theta, state[2], fluxT, q, qc];
const out = [null, new Float64Array(K * C), null, null, new Float64Array(K * C)];
for (let i = 0; i < C; i++) {
  fluxT[i] = surfaceT[i];
  if (landMask[i]) { adir[i] = adif[i] = land.albedo(i); openSea[i] = 0; continue; }
  const h = ice[i], cover = seaIce.cover(i, h);
  adir[i] = adif[i] = seaIce.albedo(h, null, seaIce.snow[i], cover);
  openSea[i] = 1 - cover;
  if (h > 0 && cover < 1) fluxT[i] = cover * surfaceT[i] + (1 - cover) * 271.35;
}
const t0 = model.time, HOURS = Number(process.env.HOURS ?? 24), start = performance.now();
radiation.restartSums();
for (let h = 0; h < HOURS; h++) {
  radiation.setTime(t0 + h * 86400 / HOURS);
  radiation.apply(fluxState, out, model.surface.windSpeed, null, 0, C, adir, adif, null, openSea, bl.depth, 0);
}
radiation.readMeans(HOURS);
const g2 = { asr: 0, olr: 0, sw: 0, lw: 0, a: 0, swT: 0, aT: 0 };
for (let i = 0; i < C; i++) {
  const a = area[i]; g2.a += a;
  g2.asr += a * radiation.meanAbsorbedSolar[i]; g2.olr += a * radiation.meanOutgoingLongwave[i]; g2.sw += a * radiation.meanShortwaveCloudEffect[i]; g2.lw += a * radiation.meanLongwaveCloudEffect[i];
  if (Math.abs(mesh.latCell[i] * deg) < 40) { g2.swT += a * radiation.meanShortwaveCloudEffect[i]; g2.aT += a; }
}
const f = (x) => (x / g2.a).toFixed(2);
console.log(`${process.env.RADIATION ?? '{}'}: ASR ${f(g2.asr)} OLR ${f(g2.olr)} SWCRE ${f(g2.sw)} (40S-40N ${(g2.swT / g2.aT).toFixed(2)}) LWCRE ${f(g2.lw)} clear-sky OLR ${((g2.olr + g2.lw) / g2.a).toFixed(2)}; ${((performance.now() - start) / 1000).toFixed(0)} s`);
