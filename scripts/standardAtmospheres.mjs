// The standard atmospheres and their reference fluxes (data/radiationBenchmark.json,
// sources listed in the file), laid onto a model column, for
// scripts/radiationBenchmark.mjs and scripts/longwaveFit.mjs.
import { readFileSync } from 'node:fs';
import { sigmaInterfaces, GRAVITY, CP_DRY } from '../js/dynamics/sigmaCore.module.js';

export const BENCHMARK = JSON.parse(readFileSync(new URL('../data/radiationBenchmark.json', import.meta.url)));
export const MOLAR = { air: 28.964, h2o: 18.015, co2: 44.01, o3: 47.997, ch4: 16.04, n2o: 44.013 };

// Mass-weighted layer means of temperature (linear in ln p within each
// reference layer) and of the mass mixing ratios (kg per kg of moist air)
// over each model layer of the sigma interfaces `levels`.
export function modelColumn(atmosphere, levels = sigmaInterfaces('bl34'), samples = 400) {
  const K = levels.length - 1, ps = atmosphere.surfacePressure, L = atmosphere.layers;
  const field = () => new Float64Array(K);
  const out = { ps, Ts: atmosphere.surfaceT, levels, T: field(), q: field(), o3: field(), co2: field(), ch4: field(), n2o: field() };
  for (let k = 0; k < K; k++) {
    const a = levels[k] * ps, b = levels[k + 1] * ps;
    const sum = { T: 0, q: 0, o3: 0, co2: 0, ch4: 0, n2o: 0 };
    for (let s = 0; s < samples; s++) {
      const p = a + (b - a) * (s + 0.5) / samples;
      let j = L.findIndex((l) => p <= l.pBottom && p >= l.pTop);
      if (j < 0) j = p > L[0].pBottom ? 0 : L.length - 1;
      const l = L[j];
      const x = Math.min(1, Math.max(0, Math.log(l.pBottom / Math.max(p, 1e-6)) / Math.log(l.pBottom / l.pTop)));
      const w = l.h2o * MOLAR.h2o / MOLAR.air, moist = 1 + w;
      sum.T += l.tBottom + (l.tTop - l.tBottom) * x;
      sum.q += w / moist;
      for (const gas of ['o3', 'co2', 'ch4', 'n2o']) sum[gas] += l[gas] * MOLAR[gas] / MOLAR.air / moist;
    }
    for (const key of Object.keys(sum)) out[key][k] = sum[key] / samples;
  }
  return out;
}

// The column's mean volume mixing ratio of a well-mixed gas (dry air).
export function meanMixingRatio(column, gas) {
  let m = 0, w = 0;
  for (let k = 0; k < column.T.length; k++) {
    const dp = column.levels[k + 1] - column.levels[k];
    m += dp * column[gas][k] / (1 - column.q[k]);
    w += dp;
  }
  return m / w * MOLAR.air / MOLAR[gas];
}

// A level quantity of a reference table (levels bottom first, pressure in Pa)
// interpolated linearly in ln p.
export function referenceAt(table, p, field = 'net') {
  for (let j = 0; j < table.length - 1; j++) {
    const a = table[j], b = table[j + 1];
    if (p <= a.p && p >= b.p) return a[field] + (b[field] - a[field]) * Math.log(a.p / p) / Math.log(a.p / b.p);
  }
  return p > table[0].p ? table[0][field] : table[table.length - 1][field];
}

// A profile given at the column's interfaces (top first) interpolated in ln p.
export function interfaceAt(column, values, p) {
  const { levels, ps } = column;
  for (let k = 0; k < levels.length - 1; k++) {
    const a = Math.max(levels[k] * ps, 1e-3), b = levels[k + 1] * ps;
    if (p >= a && p <= b) return values[k] + (values[k + 1] - values[k]) * Math.log(p / a) / Math.log(b / a);
  }
  return p < levels[1] * ps ? values[0] : values[levels.length - 1];
}

// Layer heating (K/day) from a net upward flux at the interfaces, top first.
export function layerHeating(column, netUp) {
  const { levels, ps } = column;
  return Array.from({ length: levels.length - 1 }, (_, k) => (netUp[k + 1] - netUp[k]) * GRAVITY / (CP_DRY * (levels[k + 1] - levels[k]) * ps) * 86400);
}

// The same for a reference table's net upward flux over the column's layers.
export function referenceHeating(column, table, field = 'net') {
  const { levels, ps } = column, top = table[table.length - 1].p;
  const net = Array.from(levels, (s) => referenceAt(table, Math.max(s * ps, top), field));
  return layerHeating(column, net);
}
