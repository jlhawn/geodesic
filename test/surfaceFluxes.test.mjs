import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { saturationHumidity, LATENT_HEAT } from '../js/physics/moist.module.js';
import { bulkFluxes, LOWEST_HEIGHT, WATTS_PER_MM_DAY } from '../explain/surfaceFluxes.module.js';
import { PUBLISHED_SURFACE } from '../explain/published.module.js';

const levels = sigmaInterfaces('bl36');
const near = (a, b, rel, what) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(a), Math.abs(b), 1e-30), `${what}: ${a} vs ${b}`);
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);

function modelFluxes(model, { seaT, airT, rh, wind, ps = 101325, depth = 0 }) {
  const { K, C, sigmaMid, exnerLayer } = model.core.diagnostics;
  const [pi, theta, , surfaceT, q, qc] = model.state, bottom = (K - 1) * C;
  pi[0] = ps;
  surfaceT[0] = seaT;
  model.core.diagnoseColumn(0, pi, theta, q, qc);
  for (let k = 0; k < K; k++) {
    const idx = k * C, p = ps * sigmaMid[k], T = Math.max(210, airT - 6.5e-3 * 287.06 * airT / 9.806 * Math.log(ps * sigmaMid[K - 1] / p));
    theta[idx] = T / exnerLayer[idx];
    q[idx] = k === K - 1 ? rh * saturationHumidity(airT, p) : 0;
    qc[idx] = 0;
  }
  model.core.diagnoseColumn(0, pi, theta, q, qc);
  const cell = model.exchange.cell(0, pi, theta, q, qc, seaT, wind, 0, 0, 0, 0, 0, depth);
  const CH = model.exchange.heat[0], speed = model.exchange.wind[0];
  model.radiation.column(0, ps, theta, seaT, speed, undefined, 0, q[bottom], q, qc, undefined, undefined, 1, CH);
  return { z: cell.z, Ri: cell.ri, CH, CD: model.exchange.drag[0], speed, E: model.radiation.budget.evaporation, SH: model.radiation.budget.sensibleHeat };
}

const CASES = [];
for (const seaT of [271.5, 288, 300, 304]) for (const airColder of [-2, 0, 1, 3]) for (const rh of [0.6, 0.8, 1]) for (const wind of [1, 3, 7, 12]) CASES.push({ seaT, airT: seaT - airColder, rh, wind });

test('bulkFluxes gives the model\'s own C_H, E and sensible heat over the open sea, without the gust and with it', () => {
  const still = createModel(new Grid(16), { ocean: false, levels, surface: { ...PUBLISHED_SURFACE, convectiveGust: false, gustiness: 0 }, radiation: { gustiness: 0 } });
  const gusty = createModel(new Grid(16), { ocean: false, levels, surface: PUBLISHED_SURFACE });
  assert.ok(!still.exchange.gusty && gusty.exchange.gusty);
  let worst = 0;
  for (const c of CASES) for (const [model, depth] of [[still, 0], [gusty, 600]]) {
    const m = modelFluxes(model, { ...c, depth }), f = bulkFluxes({ ...c, z: m.z, gust: model === gusty ? depth : null });
    const what = `${JSON.stringify(c)}${model === gusty ? ', gust' : ''}`;
    near(f.CH, m.CH, 1e-10, `${what} C_H`);
    near(f.CD, m.CD, 1e-10, `${what} C_D`);
    near(f.speed, m.speed, 1e-12, `${what} wind`);
    near(f.Ri, m.Ri, 1e-10, `${what} Ri`);
    near(f.E, m.E, 1e-10, `${what} E`);
    near(f.SH, m.SH, 1e-10, `${what} SH`);
    worst = Math.max(worst, Math.abs(f.CH - m.CH) / m.CH, m.E > 0 ? Math.abs(f.E - m.E) / m.E : 0);
    if (model === still && c.airT === 288) close(m.z, LOWEST_HEIGHT, 0.2, `${what} the lowest bl36 layer's height`);
  }
  console.log(`${CASES.length} cases twice: C_H and E within ${worst.toExponential(1)} of the model's`);
});

test('a 300 K sea under air 1 K cooler at 80 percent and 7 m/s evaporates 4.39 mm/day, 127.0 W/m², with 7.6 W/m² of sensible heat', () => {
  const base = { seaT: 300, airT: 299, rh: 0.8, wind: 7 };
  const f = bulkFluxes(base);
  console.log(`C_H ${f.CH.toExponential(3)} (neutral ${f.CHneutral.toExponential(3)}), Ri ${f.Ri.toFixed(4)}, E ${f.mmPerDay.toFixed(3)} mm/d, L·E ${f.LE.toFixed(2)}, SH ${f.SH.toFixed(2)}, Bowen ${f.bowen.toFixed(4)}`);
  close(f.CH, 1.145e-3, 0.0005e-3, 'C_H');
  close(f.CHneutral, 0.989e-3, 0.0005e-3, 'neutral C_H');
  close(f.mmPerDay, 4.39, 0.005, 'E');
  close(f.LE, 127.0, 0.05, 'L·E');
  close(f.SH, 7.62, 0.005, 'SH');
  close(f.bowen, 0.060, 0.0005, 'Bowen ratio');
  close(f.E * 86400, f.mmPerDay, 1e-12, 'mm/day');
  close(f.LE, LATENT_HEAT * f.E, 1e-9, 'L·E');
  close(WATTS_PER_MM_DAY, 28.94, 0.005, '1 mm/day');
  const rate = (change) => bulkFluxes({ ...base, ...change }).mmPerDay;
  const table = { 288: rate({ seaT: 288, airT: 287 }), 273: rate({ seaT: 273, airT: 272 }), '3 m/s': rate({ wind: 3 }), '12 m/s': rate({ wind: 12 }), '60 %': rate({ rh: 0.6 }), 'air as warm': rate({ airT: 300 }), '3 K cooler': rate({ airT: 297 }) };
  console.log(Object.entries(table).map(([k, v]) => `${k}: ${v.toFixed(3)}`).join(', '));
  close(table[288], 2.14, 0.005, '288 K sea');
  close(table[273], 0.81, 0.005, '273 K sea');
  close(table['3 m/s'], 2.31, 0.005, '3 m/s');
  close(table['12 m/s'], 6.94, 0.005, '12 m/s');
  close(table['60 %'], 8.02, 0.005, 'air at 60 percent');
  close(bulkFluxes({ ...base, seaT: 288, airT: 287 }).LE, 61.9, 0.05, '288 K L·E');
  close(bulkFluxes({ ...base, seaT: 273, airT: 272 }).LE, 23.35, 0.005, '273 K L·E');
});

test('evaporation grows with the sea\'s warmth by about Clausius–Clapeyron, never goes negative, and the sensible heat follows the dry static energy', () => {
  let previous = 0;
  for (let seaT = 271; seaT <= 305; seaT++) {
    const f = bulkFluxes({ seaT, airT: seaT - 1, rh: 0.8, wind: 7 });
    assert.ok(f.E > previous, `${seaT} K: ${f.E}`);
    if (previous > 0) assert.ok(f.E / previous > 1.05 && f.E / previous < 1.09, `${seaT} K grows ${f.E / previous}`);
    previous = f.E;
  }
  const saturatedWarm = bulkFluxes({ seaT: 290, airT: 292, rh: 1, wind: 5 });
  assert.equal(saturatedWarm.E, 0);
  assert.equal(saturatedWarm.bowen, null);
  assert.ok(saturatedWarm.Ri > 0 && saturatedWarm.SH < 0 && saturatedWarm.CH < saturatedWarm.CHneutral, 'stable air over a cooler sea');
  const even = bulkFluxes({ seaT: 300, airT: 300 - 9.806 * LOWEST_HEIGHT / 1003.5, rh: 0.8, wind: 7 });
  close(even.SH, 0, 1e-9, 'no sensible heat when the air\'s dry static energy matches the sea\'s');
  const unstable = bulkFluxes({ seaT: 300, airT: 297, rh: 0.8, wind: 3 });
  assert.ok(unstable.Ri < 0 && unstable.CH > unstable.CHneutral, 'unstable air exchanges faster than neutral');
});
