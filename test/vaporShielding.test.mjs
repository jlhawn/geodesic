import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore, sigmaInterfaces, GRAVITY } from '../js/dynamics/sigmaCore.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { createRadiation, DAY, SOLAR_CONSTANT } from '../js/physics/radiation.module.js';
import { nearInfraredBands, nearInfraredShielded, NEAR_INFRARED_MOLAR, OXYGEN, OZONE_CM_ATM } from '../js/physics/shortwaveGases.module.js';
import { BENCHMARK, modelColumn, meanMixingRatio } from '../scripts/standardAtmospheres.mjs';

const levels = sigmaInterfaces('bl48'), mesh = buildMesh(new Grid(16)), core = createSigmaCore(mesh, { levels });
const { K, C, dSigma, sigmaMid } = core.diagnostics;
const S0 = 1360.85;
const TROPICAL = modelColumn(BENCHMARK.atmospheres.TROP, levels), MIDLATITUDE = modelColumn(BENCHMARK.atmospheres.MLS, levels);

function gasesOf(column) {
  const ozone = Array.from(column.o3, (x, k) => x * (levels[k + 1] - levels[k]) * column.ps / GRAVITY / OZONE_CM_ATM);
  return { carbonDioxide: meanMixingRatio(column, 'co2'), methane: meanMixingRatio(column, 'ch4'), nitrousOxide: meanMixingRatio(column, 'n2o'), ozoneProfile: Float64Array.from(ozone) };
}
const layerAt = (column, p) => { let k = 0; for (let j = 1; j < K; j++) if (Math.abs(sigmaMid[j] * column.ps - p) < Math.abs(sigmaMid[k] * column.ps - p)) k = j; return k; };

const EQUATORIAL = Array.from({ length: C }, (_, i) => i).reduce((best, i) => (Math.abs(mesh.latCell[i]) < Math.abs(mesh.latCell[best]) ? i : best), 0);

function sunAt(radiation, zenith) {
  const mu = Math.cos(zenith * Math.PI / 180);
  let noon = 0, highest = -1;
  for (let s = 0; s < 1440; s++) { radiation.setTime(60 * s); const c = radiation.cosZenith(EQUATORIAL); if (c > highest) { highest = c; noon = 60 * s; } }
  assert.ok(highest > mu, `the sun reaches cosine ${highest} at noon, short of ${mu}`);
  let low = noon - 6 * 3600, high = noon;
  for (let n = 0; n < 50; n++) { const middle = 0.5 * (low + high); radiation.setTime(middle); if (radiation.cosZenith(EQUATORIAL) < mu) low = middle; else high = middle; }
  return 0.5 * (low + high);
}

function columnBudgets(column, { zenith, clouds, albedo = 0.06, options = {} }) {
  const out = {};
  let time = null;
  for (const vaporShielding of ['none', 'cloudTop']) {
    const radiation = createRadiation(mesh, core, { landAerosol: 0, seaAerosol: 0, clearSkyPass: true, cloudCover: 'overcast', cloudColumns: false, solarConstant: S0, ...gasesOf(column), ...options, vaporShielding });
    time ??= sunAt(radiation, zenith);
    radiation.setTime(time);
    const pi = new Float64Array(C).fill(column.ps), theta = new Float64Array(K * C).fill(280), q = new Float64Array(K * C), qc = new Float64Array(K * C);
    core.diagnose(pi, theta);
    const { exnerLayer } = core.diagnostics;
    for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) { theta[k * C + i] = column.T[k] / exnerLayer[k * C + i]; q[k * C + i] = column.q[k]; }
    for (const [p, path] of clouds) { const k = layerAt(column, p); for (let i = 0; i < C; i++) qc[k * C + i] = path * GRAVITY / (column.ps * dSigma[k]); }
    core.diagnose(pi, theta, q, qc);
    radiation.column(EQUATORIAL, column.ps, theta, column.Ts, 5, undefined, undefined, q[(K - 1) * C + EQUATORIAL], q, qc, albedo, albedo, 0, 0, 0, 0, 0, 0);
    let heating = 0;
    for (let k = 0; k < K; k++) heating += radiation.layerFlux[k] - radiation.longwave[k * C + EQUATORIAL];
    out[vaporShielding] = { ...radiation.budget, beam: radiation.budget.insolation, heating };
  }
  return out;
}

const CASES = [
  ['a stratocumulus', [[92500, 0.1]]], ['a mixed-phase layer', [[45000, 0.3]]], ['a thick anvil', [[25000, 0.3]]], ['a thin anvil', [[25000, 0.1]]],
  ['two layers', [[25000, 0.1], [85000, 0.1]]], ['thin cirrus', [[20000, 0.005]]],
];
const gas = (b) => b.vaporSolar + b.oxygenSolar + b.carbonDioxideSolar + b.methaneSolar + b.shieldedGasSolar;

test('vaporShielding is none or cloudTop, and cloudTop needs the band gases, a scattering atmosphere and a cover it can follow', () => {
  const make = (extra) => createRadiation(mesh, core, extra);
  assert.doesNotThrow(() => make({}));
  assert.doesNotThrow(() => make({ vaporShielding: 'none' }));
  assert.doesNotThrow(() => make({ vaporShielding: 'cloudTop' }));
  assert.doesNotThrow(() => make({ vaporShielding: 'cloudTop', cloudCover: 'overcast', cloudColumns: false }));
  for (const bad of [{ vaporShielding: 'top' }, { vaporShielding: 'cloudTop', nearInfrared: 'clirad' }, { vaporShielding: 'cloudTop', solarGases: 'lacisHansen' }, { vaporShielding: 'cloudTop', upwardAbsorption: false },
    { vaporShielding: 'cloudTop', rayleighDepth: 0, nearInfraredRayleigh: 0, landAerosol: 0, seaAerosol: 0 }, { vaporShielding: 'cloudTop', cloudOverlap: 'maximum', cloudColumns: false }]) {
    assert.throws(() => make(bad), /vaporShielding/, JSON.stringify(bad));
  }
});

test('nearInfraredShielded with the diffuse path set to the beam\'s own gives the beam\'s absorption below the top, and its climb from the surface when the top is the surface', () => {
  for (const column of [TROPICAL, MIDLATITUDE]) {
    const air = { K, vapour: new Float64Array(K), dry: new Float64Array(K), p: new Float64Array(K), T: new Float64Array(K), e: new Float64Array(K) };
    for (let k = 0; k < K; k++) {
      const mass = column.ps * dSigma[k] / GRAVITY, q = column.q[k];
      air.vapour[k] = q * mass; air.dry[k] = mass * (1 - q); air.p[k] = column.ps * sigmaMid[k]; air.T[k] = column.T[k]; air.e[k] = q * air.p[k] / (0.622 + 0.378 * q);
    }
    const mixing = { co2: meanMixingRatio(column, 'co2') * NEAR_INFRARED_MOLAR.co2, ch4: meanMixingRatio(column, 'ch4') * NEAR_INFRARED_MOLAR.ch4, o2: OXYGEN.mixingRatio * NEAR_INFRARED_MOLAR.o2 };
    const m = 1.3, taken = new Float64Array(K), upward = new Float64Array(K), shielded = { down: new Float64Array(K), up: new Float64Array(K), climb: new Float64Array(K) };
    nearInfraredBands(air, m, m, mixing, 1.48, taken, upward);
    for (const top of [0, 10, 30, K]) {
      nearInfraredShielded(air, m, m, mixing, 1.48, top, shielded);
      for (let k = top; k < K; k++) {
        assert.ok(Math.abs(shielded.down[k] - (taken[k] - (k > 0 ? taken[k - 1] : 0))) < 1e-12, `top ${top}: down[${k}] ${shielded.down[k]} against ${taken[k] - (k > 0 ? taken[k - 1] : 0)}`);
        assert.ok(Math.abs(shielded.up[k] - upward[k]) < 1e-12, `top ${top}: up[${k}] ${shielded.up[k]} against ${upward[k]}`);
      }
      for (let k = 0; k < top; k++) assert.ok(shielded.climb[k] >= 0, `climb[${k}] ${shielded.climb[k]}`);
    }
    nearInfraredShielded(air, m, 5 / 3, mixing, 1.48, K, shielded);
    nearInfraredBands(air, m, 5 / 3, mixing, 1.48, taken, upward);
    for (let k = 0; k < K; k++) assert.ok(Math.abs(shielded.climb[k] - upward[k]) < 1e-12, `a cloud at the surface: climb[${k}] ${shielded.climb[k]} against ${upward[k]}`);
  }
});

test('under cloudTop the absorbed and the reflected sunlight make the beam, the layers\' heating makes the atmosphere\'s absorption, and the clear-sky pass is the sequential scheme\'s', () => {
  let worstClosure = 0, worstHeating = 0, columns = 0;
  for (const [name, column] of [['tropical', TROPICAL], ['midlatitude', MIDLATITUDE]]) {
    for (const [label, clouds] of CASES) {
      for (const zenith of [10, 30, 60, 80]) {
        const { none, cloudTop } = columnBudgets(column, { zenith, clouds, options: { seaDropletRadius: 10 } });
        for (const b of [none, cloudTop]) {
          worstClosure = Math.max(worstClosure, Math.abs(b.absorbedSolar + b.reflectedSolar - b.beam) / b.beam);
          worstHeating = Math.max(worstHeating, Math.abs(b.heating - b.atmosphereSolar) / b.beam);
          assert.ok(b.reflectedSolar >= 0 && b.reflectedSolar <= b.beam, `${name} ${label} at ${zenith} degrees: reflected ${b.reflectedSolar} of ${b.beam}`);
          assert.ok(b.surfaceSolar >= 0 && b.atmosphereSolar >= 0, `${name} ${label} at ${zenith} degrees: surface ${b.surfaceSolar}, atmosphere ${b.atmosphereSolar}`);
        }
        assert.equal(cloudTop.clearAbsorbedSolar, none.clearAbsorbedSolar, `${name} ${label} at ${zenith} degrees: clear-sky absorption`);
        assert.equal(cloudTop.beam, none.beam);
        columns++;
      }
    }
  }
  console.log(`${columns} overcast columns: absorbed plus reflected differs from the beam by at most ${worstClosure.toExponential(1)} of it and the layers' heating from the atmosphere's absorption by ${worstHeating.toExponential(1)}`);
  assert.ok(worstClosure < 1e-12 && worstHeating < 1e-12);
});

test('a column without cloud is the sequential scheme\'s bit for bit, and a vanishing cloud approaches it', () => {
  for (const column of [TROPICAL, MIDLATITUDE]) for (const zenith of [20, 70]) {
    const { none, cloudTop } = columnBudgets(column, { zenith, clouds: [] });
    for (const key of Object.keys(none)) assert.equal(cloudTop[key], none[key], `${key} without cloud`);
    const faint = columnBudgets(column, { zenith, clouds: [[30000, 1e-9]], options: { nearInfraredRayleigh: 0 } });
    for (const key of ['absorbedSolar', 'atmosphereSolar', 'surfaceSolar', 'reflectedSolar']) assert.ok(Math.abs(faint.cloudTop[key] - faint.none[key]) < 1e-6 * faint.none.beam, `${key} under a cloud of 1e-9 kg/m2: ${faint.cloudTop[key]} against ${faint.none[key]}`);
  }
});

test('cloudTop takes the gas absorption below a thick high cloud away, adds to it under a stratocumulus, and moves the planet\'s reflection by what the gases gave up', () => {
  const rows = [];
  for (const zenith of [30, 60]) {
    const anvil = columnBudgets(TROPICAL, { zenith, clouds: [[25000, 0.3]] }), deck = columnBudgets(TROPICAL, { zenith, clouds: [[92500, 0.1]], options: { seaDropletRadius: 10 } });
    const cut = gas(anvil.none) - gas(anvil.cloudTop), gain = gas(deck.cloudTop) - gas(deck.none);
    rows.push(`${zenith} degrees: anvil gas ${(100 * gas(anvil.none) / anvil.none.beam).toFixed(1)} to ${(100 * gas(anvil.cloudTop) / anvil.none.beam).toFixed(1)} %, reflected ${(100 * anvil.none.reflectedSolar / anvil.none.beam).toFixed(1)} to ${(100 * anvil.cloudTop.reflectedSolar / anvil.none.beam).toFixed(1)} %; stratocumulus gas ${(100 * gas(deck.none) / deck.none.beam).toFixed(1)} to ${(100 * gas(deck.cloudTop) / deck.none.beam).toFixed(1)} %`);
    assert.ok(cut > 0.4 * gas(anvil.none) && cut < 0.9 * gas(anvil.none), `anvil: the gases give up ${cut} of ${gas(anvil.none)}`);
    assert.ok(anvil.cloudTop.reflectedSolar - anvil.none.reflectedSolar > 0.5 * cut, 'most of what the gases give up the planet reflects');
    assert.ok(anvil.cloudTop.surfaceSolar > 0.9 * anvil.none.surfaceSolar && anvil.cloudTop.surfaceSolar < 1.05 * anvil.none.surfaceSolar, 'the surface keeps its sunlight to a few per cent');
    assert.ok(gain > 0 && gain < 0.15 * gas(deck.none), `stratocumulus: the gases gain ${gain} on ${gas(deck.none)}`);
  }
  console.log(rows.join('\n'));
});

test('under cloudTop the cloudy N = 16 model\'s absorbed plus reflected sunlight makes the insolation of every step, held or not, and the shielding acts in its columns', () => {
  for (const every of [1, 4]) {
    const model = createModel(new Grid(16), { ocean: false, radiation: { stratus: true, clearSkyPass: true, radiationEvery: every, vaporShielding: 'cloudTop' }, surface: { exchange: 'roughness' } });
    initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
    const radiation = model.radiation, cells = model.mesh.nCells, dt = 900;
    model.time = 150 * DAY;
    for (let n = 0; n < 40; n++) model.step(dt);
    const absorbed = new Float64Array(cells), reflected = new Float64Array(cells), insolation = new Float64Array(cells);
    let worst = 0, lit = 0, shielded = 0;
    for (let n = 0; n < 9; n++) {
      const time = model.time;
      model.step(dt);
      radiation.setTime(time);
      for (let i = 0; i < cells; i++) {
        const beam = radiation.summed.insolation[i] - insolation[i], top = radiation.summed.absorbedSolar[i] - absorbed[i], lost = radiation.summed.reflectedSolar[i] - reflected[i];
        if (every === 1) { worst = Math.max(worst, Math.abs(top + lost - beam)); if (beam > 0) lit++; }
        else if (radiation.held.mu[i] > 0) { worst = Math.max(worst, Math.abs(radiation.held.top[i] + radiation.held.reflected[i] - SOLAR_CONSTANT * radiation.held.mu[i])); lit++; }
      }
      absorbed.set(radiation.summed.absorbedSolar); reflected.set(radiation.summed.reflectedSolar); insolation.set(radiation.summed.insolation);
      radiation.setTime(model.time);
    }
    const [pi, theta, , surfaceT, q, qc] = model.state, bottom = (K - 1) * cells;
    radiation.setTime(model.time);
    for (let i = 0; i < cells; i++) {
      radiation.column(i, pi[i], theta, surfaceT[i], 5, undefined, undefined, q[bottom + i], q, qc, 0.1, 0.1, 1, 0, 0, 0, 0, 0);
      if (radiation.budget.shieldedGasSolar !== 0) shielded++;
    }
    console.log(`radiationEvery ${every}: ${every === 1 ? 'absorbed plus reflected differs from the insolation of the step' : 'the held call\'s absorbed plus reflected differs from its beam'} by at most ${worst.toExponential(1)} W/m² over ${lit} lit cell-steps; ${shielded} of ${cells} columns shield the gases at the end`);
    assert.ok(lit > 2 * cells);
    assert.ok(worst < 1e-8, `absorbed plus reflected misses the sunlight by ${worst} W/m²`);
    assert.ok(shielded > 0, 'no column was shielded');
  }
});
