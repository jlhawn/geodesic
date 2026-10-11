// Lays the midlatitude-summer standard atmosphere onto the bl36 layers for
// explain/radiation.html:
//   node scripts/explainColumn.mjs [OUT]
// writing explain/data/mls_bl36.json.
import { writeFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore, sigmaInterfaces, GRAVITY } from '../js/dynamics/sigmaCore.module.js';
import { OZONE_CM_ATM } from '../js/physics/shortwaveGases.module.js';
import { BENCHMARK, modelColumn } from './standardAtmospheres.mjs';

const out = process.argv[2] ?? new URL('../explain/data/mls_bl36.json', import.meta.url).pathname;
const atmosphere = BENCHMARK.atmospheres.MLS, levels = sigmaInterfaces('bl36'), column = modelColumn(atmosphere, levels);
const mesh = buildMesh(new Grid(2)), core = createSigmaCore(mesh, { levels });
const { K, C, sigmaMid, dSigma, exnerLayer, geopotential } = core.diagnostics, ps = column.ps;
const pi = new Float64Array(C).fill(ps), theta = new Float64Array(K * C);
core.diagnoseColumn(0, pi, theta);
for (let k = 0; k < K; k++) theta[k * C] = column.T[k] / exnerLayer[k * C];
core.diagnoseColumn(0, pi, theta);

const round = (values, digits) => Array.from(values, (x) => +x.toFixed(digits));
const significant = (values, digits) => Array.from(values, (x) => +x.toPrecision(digits));
const ozone = Float64Array.from(column.o3, (x, k) => x * dSigma[k] * ps / GRAVITY / OZONE_CM_ATM);
const data = {
  name: atmosphere.name,
  source: `${BENCHMARK.sources.afgl}; laid onto bl36 by scripts/standardAtmospheres.mjs modelColumn`,
  grid: 'bl36',
  levels: Array.from(levels),
  surfacePressure: ps,
  groundTemperature: column.Ts,
  pressure: round(Array.from(sigmaMid, (s) => s * ps), 2),
  height: round(Array.from({ length: K }, (_, k) => geopotential[k * C] / GRAVITY), 1),
  temperature: round(column.T, 3),
  ozone: significant(ozone, 6),
  ozoneColumn: +(ozone.reduce((a, b) => a + b, 0) * 1000).toFixed(1),
  humidity: significant(column.q, 6),
};
writeFileSync(out, `${JSON.stringify(data)}\n`);
console.log(`${out}: ${K} layers, Ts ${data.groundTemperature} K, ozone ${data.ozoneColumn} DU, top layer ${data.height[0]} m`);
