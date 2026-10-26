// Runs the climate page's N = 64 default state (runs/eleven64_day1825.bin)
// for one day on the CPU and sums the ground's energy budget for
// explain/water.html:
//   [STEPS=n] node scripts/explainSurfaceBudget.mjs [STATE] [OUT]
// writing explain/data/surface_budget.json. About 25 minutes on one core.
import { writeFileSync } from 'node:fs';
import { createModel } from '../js/model.module.js';
import { savedLevels } from '../js/stateFile.module.js';
import { withCadence } from '../js/cadence.module.js';
import { LATENT_HEAT, placeCumulus } from '../js/physics/moist.module.js';
import { readTopography, readState, tagOf } from './figures/figureState.mjs';
import { stateGrid } from '../js/stateGrid.module.js';

const file = process.argv[2] ?? new URL('../runs/eleven64_day1825.bin', import.meta.url).pathname;
const out = process.argv[3] ?? new URL('../explain/data/surface_budget.json', import.meta.url).pathname;
const saved = await readState(file);
const dt = 1350 * 16 / saved.N, steps = Number(process.env.STEPS ?? Math.round(86400 / dt));
const topography = readTopography(saved);
const model = createModel(stateGrid(saved, topography).mesh, { ...withCadence({ ocean: false }, dt, saved.N), ocean: false, levels: savedLevels(saved), topography });
const { moist, boundaryLayer, radiation, mesh, geography } = model, C = mesh.nCells;
['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => model.state[a].set(saved[name]));
if (saved.boundaryDepth) boundaryLayer.depth.set(saved.boundaryDepth);
model.time = saved.time;
model.seaIce.load(model.state[6], saved.concentration ?? null);
if (model.land && saved.land) model.land.load(saved.land);
if (saved.windSpeed) model.surface.windSpeed.set(saved.windSpeed);
if (saved.exchangeWind && model.exchange) model.exchange.wind.set(saved.exchangeWind);
if (saved.mlmGate) radiation.mlmGate.set(saved.mlmGate);
if (saved.mlmHeight) radiation.mlmHeight.set(saved.mlmHeight);
if (saved.subcloudVirtual && saved.subcloudVirtual.length === moist.subcloudVirtual.length) moist.subcloudVirtual.set(saved.subcloudVirtual);
placeCumulus(model, saved);

const onLand = model.land ? geography.land : new Uint8Array(C);
const net = new Float64Array(C), latent = new Float64Array(C), sensible = new Float64Array(C);
const start = performance.now();
model.restartPrecipitation();
for (let n = 0; n < steps; n++) {
  model.step(dt);
  for (let i = 0; i < C; i++) {
    const LE = LATENT_HEAT * radiation.evaporation[i], SH = radiation.sensibleHeat[i];
    net[i] += dt * radiation.surfaceFlux[i]; latent[i] += dt * LE; sensible[i] += dt * SH;
  }
  if (n % 32 === 31) console.log(`step ${n + 1} of ${steps}, ${((performance.now() - start) / 1000).toFixed(0)} s`);
}
const seconds = steps * dt, rain = moist.precipitation;
const round = (x, d = 2) => +x.toFixed(d);
const budget = (keep) => {
  let area = 0, sunlight = 0, gain = 0, LE = 0, SH = 0, P = 0;
  for (let i = 0; i < C; i++) {
    if (!keep(i)) continue;
    const a = mesh.areaCell[i];
    area += a;
    sunlight += a * (radiation.summed.absorbedSolar[i] - radiation.summed.atmosphereSolar[i]) * dt;
    gain += a * net[i]; LE += a * latent[i]; SH += a * sensible[i]; P += a * rain[i];
  }
  const mean = (x) => x / area / seconds, stored = mean(gain);
  return {
    share: area, sunlight: round(mean(sunlight)), infrared: round(mean(sunlight) - stored - mean(LE) - mean(SH)), latent: round(mean(LE)), sensible: round(mean(SH)), stored: round(stored),
    evaporation: round(86400 * mean(LE) / LATENT_HEAT, 3), precipitation: round(86400 * mean(P), 3),
  };
};
const all = budget(() => true), sea = budget((i) => !onLand[i]), land = budget((i) => onLand[i]);
for (const part of [sea, land]) part.share = round(part.share / all.share, 4);
all.share = 1;
const data = {
  name: 'the ground\'s energy budget',
  source: `scripts/explainSurfaceBudget.mjs: ${tagOf(file)} day ${saved.day} on createModel at N = ${saved.N}, ${steps} steps of ${dt} s; day means in W/m², upward positive, net infrared the remainder; ocean: false, a slab sea`,
  day: saved.day,
  N: saved.N,
  steps,
  dt,
  all, sea, land,
};
console.log(JSON.stringify(data, null, 1));
writeFileSync(out, `${JSON.stringify(data)}\n`);
