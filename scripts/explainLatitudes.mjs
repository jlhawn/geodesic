// Runs the ten columns of explain/latitudeColumns.module.js from 288 K to
// equilibrium for explain/radiation.html:
//   node scripts/explainLatitudes.mjs [OUT]
// writing explain/data/latitudes.json.
import { writeFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { meanSunlight, HEAT_CAPACITY } from '../explain/radiativeColumn.module.js';
import { createLatitudes, STEP, STEPS_PER_DAY } from '../explain/latitudeColumns.module.js';

const out = process.argv[2] ?? new URL('../explain/data/latitudes.json', import.meta.url).pathname;
const DAYS = 3000, D = Math.PI / 180;
const mesh = buildMesh(new Grid(16));
const latitudes = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90].map((target) => {
  let best = 0;
  for (let i = 0; i < mesh.nCells; i++) if (Math.abs(mesh.latCell[i] - target * D) < Math.abs(mesh.latCell[best] - target * D)) best = i;
  return Math.abs(mesh.latCell[best]) / D;
});
const sun = latitudes.map((latitude) => meanSunlight(latitude));
const band = createLatitudes({ latitudes, insolation: sun.map((s) => s.insolation), mu: sun.map((s) => s.mu) });
for (const column of band.columns) for (let s = 0; s < DAYS * STEPS_PER_DAY; s++) column.step(STEP);
for (const column of band.columns) column.evaluate();

const round = (values, digits) => Array.from(values, (x) => +x.toFixed(digits));
const data = {
  name: 'every latitude for itself',
  source: 'scripts/explainLatitudes.mjs: explain/radiativeColumn.module.js at Grid(16) cell latitudes, annual-mean sunlight and AFGL ozone, q = 0, CO₂ ×1, ground albedo 0.07, from 288 K',
  grid: 'bl36',
  step: STEP,
  days: DAYS,
  heatCapacity: HEAT_CAPACITY,
  settings: { carbonDioxide: 1, albedo: 0.07 },
  latitudes: round(latitudes, 6),
  insolation: round(sun.map((s) => s.insolation), 6),
  mu: round(sun.map((s) => s.mu), 8),
  daylight: round(sun.map((s) => s.daylight), 6),
  groundTemperature: round(band.columns.map((c) => c.Ts), 4),
  temperature: band.columns.map((c) => round(c.T, 4)),
};
writeFileSync(out, `${JSON.stringify(data)}\n`);
const { rows, mean, contrast } = band.summary();
for (const r of rows) console.log(`${r.latitude.toFixed(2).padStart(5)}°  Ts ${r.Ts.toFixed(2)} K  500 hPa ${r.T500.toFixed(2)} K  ASR ${r.absorbed.toFixed(2)}  OLR ${r.olr.toFixed(2)}  imbalance ${r.imbalance.toExponential(1)}`);
console.log(`${out}: area mean Ts ${mean.Ts.toFixed(2)} K, ASR ${mean.absorbed.toFixed(2)} W/m², Te ${mean.Te.toFixed(2)} K, contrast ${contrast.toFixed(2)} K`);
