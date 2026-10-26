import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore, sigmaInterfaces, GRAVITY } from '../js/dynamics/sigmaCore.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { createRadiation, DAY, SOLAR_CONSTANT } from '../js/physics/radiation.module.js';
import { twoStreamLayer, addLayers, NEAR_INFRARED_POINTS, NEAR_INFRARED_SHARE } from '../js/physics/layeredShortwave.module.js';
import { LAYERED_TABLE } from '../js/physics/layeredTable.module.js';
import { NEAR_INFRARED_TABLE } from '../js/physics/nearInfraredTable.module.js';
import { nearInfraredBands, NEAR_INFRARED_MOLAR, OZONE_CM_ATM } from '../js/physics/shortwaveGases.module.js';
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
  for (const [name, scheme] of [['sequential', { shortwave: 'sequential' }], ['layered', { shortwave: 'layered' }]]) {
    const radiation = createRadiation(mesh, core, { landAerosol: 0, seaAerosol: 0, clearSkyPass: true, cloudCover: 'overcast', cloudColumns: false, solarConstant: S0, ...gasesOf(column), ...options, ...scheme });
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
    out[name] = { ...radiation.budget, beam: radiation.budget.insolation, heating };
  }
  return out;
}

const CASES = [
  ['a stratocumulus', [[92500, 0.1]]], ['a mixed-phase layer', [[45000, 0.3]]], ['a thick anvil', [[25000, 0.3]]], ['a thin anvil', [[25000, 0.1]]],
  ['two layers', [[25000, 0.1], [85000, 0.1]]], ['thin cirrus', [[20000, 0.005]]],
];

const emptyLayers = (n) => ({ reflect: new Float64Array(n), scatter: new Float64Array(n), direct: new Float64Array(n), reflectDiffuse: new Float64Array(n), transmitDiffuse: new Float64Array(n) });

function slab(n, tau, c, g, mu, surfaceAlbedo, diffuseAlbedo) {
  const layers = emptyLayers(n), one = {};
  for (let j = 0; j < n; j++) {
    twoStreamLayer(tau / n, c, g, mu, one);
    for (const key of Object.keys(layers)) layers[key][j] = one[key];
  }
  const out = { heat: new Float64Array(n), surface: 0, reflected: 0, down: 0, direct: 0 };
  addLayers(layers, n, surfaceAlbedo, diffuseAlbedo, 1, out, new Float64Array(n + 1), new Float64Array(n + 1));
  return { ...out, absorbed: out.heat.reduce((a, b) => a + b, 0) };
}

test('shortwave is sequential or layered, and layered needs the band gases, a scattering atmosphere, no shielding and no gray cloud', () => {
  const make = (extra) => createRadiation(mesh, core, extra);
  assert.doesNotThrow(() => make({}));
  assert.doesNotThrow(() => make({ shortwave: 'sequential' }));
  assert.doesNotThrow(() => make({ shortwave: 'layered' }));
  for (const bad of [{ shortwave: 'twoStream' }, { shortwave: 'layered', nearInfrared: 'clirad' }, { shortwave: 'layered', solarGases: 'lacisHansen' }, { shortwave: 'layered', vaporShielding: 'cloudTop' },
    { shortwave: 'layered', cloudScattering: 0.1 }, { shortwave: 'layered', vaporAbsorption: 0.9 }, { shortwave: 'layered', rayleighDepth: 0, nearInfraredRayleigh: 0, landAerosol: 0, seaAerosol: 0 }]) {
    assert.throws(() => make(bad), /shortwave/, JSON.stringify(bad));
  }
});

test('a layer reflects, transmits and absorbs without negative parts or non-finite numbers over the range of the atmosphere\'s optics', () => {
  const one = {};
  let worst = Infinity, count = 0;
  for (const g of [0, 0.3, 0.7, 0.85, 0.95]) for (const c of [0, 1e-9, 1e-6, 1e-3, 0.01, 0.1, 0.3, 0.5, 0.9, 1]) for (const tau of [1e-6, 1e-3, 0.1, 1, 5, 30, 300]) for (const mu of [0.02, 0.1, 0.3, 0.5774, 0.7, 1]) {
    twoStreamLayer(tau, c, g, mu, one);
    for (const key of ['reflect', 'scatter', 'direct', 'reflectDiffuse', 'transmitDiffuse']) assert.ok(Number.isFinite(one[key]) && one[key] >= 0, `${key} ${one[key]} at tau ${tau} c ${c} g ${g} mu ${mu}`);
    worst = Math.min(worst, 1 - one.reflect - one.scatter - one.direct, 1 - one.reflectDiffuse - one.transmitDiffuse);
    count++;
  }
  console.log(`${count} layers: the least absorption of the beam or of diffuse light is ${worst.toExponential(1)}`);
  assert.ok(worst > -1e-9);
});

test('a layer without absorption conserves the beam and diffuse light, and the layer is smooth through the cosine where k mu is 1', () => {
  const one = {}, next = {};
  for (const g of [0, 0.5, 0.85]) for (const tau of [1e-6, 1e-3, 0.1, 5, 300]) for (const mu of [0.03, 0.3, 1]) {
    twoStreamLayer(tau, 0, g, mu, one);
    assert.ok(Math.abs(1 - one.reflect - one.scatter - one.direct) < 1e-12, `beam at tau ${tau} g ${g} mu ${mu}`);
    assert.ok(Math.abs(1 - one.reflectDiffuse - one.transmitDiffuse) < 1e-12, `diffuse at tau ${tau} g ${g}`);
  }
  for (const [tau, c, g] of [[1, 0.5, 0], [3, 0.4, 0.3], [0.5, 0.7, 0.6]]) {
    let step = 0;
    twoStreamLayer(tau, c, g, 0.1, one);
    for (let mu = 0.1 + 5e-5; mu < 1; mu += 5e-5) {
      twoStreamLayer(tau, c, g, mu, next);
      step = Math.max(step, Math.abs(next.reflect - one.reflect), Math.abs(next.scatter - one.scatter));
      Object.assign(one, next);
    }
    assert.ok(step < 5e-5,`the largest step of the layer between cosines 5e-5 apart is ${step} at tau ${tau} c ${c} g ${g}`);
  }
});

test('a layer of cloud reflects as the Monte Carlo of its scattering says, to the accuracy of the delta-Eddington two-stream', () => {
  let seed = 987654321;
  const rand = () => { seed ^= seed << 13; seed >>>= 0; seed ^= seed >>> 17; seed ^= seed << 5; seed >>>= 0; return (seed + 0.5) / 4294967296; };
  const monteCarlo = (tau, c, g, mu0, photons) => {
    let reflected = 0, transmitted = 0;
    for (let n = 0; n < photons; n++) {
      let mu = mu0, z = 0, w = 1;
      for (;;) {
        z += -Math.log(rand()) * mu;
        if (z < 0) { reflected += w; break; }
        if (z > tau) { transmitted += w; break; }
        w *= 1 - c;
        if (w < 1e-3) { if (rand() < 0.1) w *= 10; else break; }
        const q = (1 - g * g) / (1 - g + 2 * g * rand()), cosT = (1 + g * g - q * q) / (2 * g), phi = 2 * Math.PI * rand(), sinT = Math.sqrt(Math.max(0, 1 - cosT * cosT));
        mu = Math.max(-1, Math.min(1, mu * cosT + Math.sqrt(Math.max(0, 1 - mu * mu)) * sinT * Math.cos(phi)));
      }
    }
    return { reflect: reflected / photons, transmit: transmitted / photons };
  };
  const one = {};
  for (const [tau, c, mu] of [[5, 1e-6, 0.7], [20, 1e-6, 0.5], [20, 0.01, 0.9], [50, 0.1, 0.6]]) {
    const mc = monteCarlo(tau, c, 0.85, mu, 60000);
    twoStreamLayer(tau, c, 0.85, mu, one);
    console.log(`tau ${tau} co-albedo ${c} cosine ${mu}: reflected ${one.reflect.toFixed(3)} against ${mc.reflect.toFixed(3)}, transmitted ${(one.scatter + one.direct).toFixed(3)} against ${mc.transmit.toFixed(3)}`);
    assert.ok(Math.abs(one.reflect - mc.reflect) < 0.05 && Math.abs(one.scatter + one.direct - mc.transmit) < 0.05, `tau ${tau} c ${c} mu ${mu}`);
  }
});

test('layers added in any number give one thick layer\'s sunlight, and the top, the surface and the layers account for all of it', () => {
  for (const [tau, c, g, mu, a] of [[10, 0, 0.85, 0.5, 0.1], [30, 0.01, 0.85, 0.8, 0.2], [0.1, 0, 0, 0.6, 0.2], [1.2, 1, 0, 0.7, 0.3], [2, 0.3, 0.5, 0.9, 0.05], [3, 0.9, 0, 0.5774, 0.2]]) {
    const [one, many, most] = [1, 7, 40].map((n) => slab(n, tau, c, g, mu, a, a));
    for (const r of [one, many, most]) {
      assert.ok(Math.abs(1 - r.reflected - r.surface - r.absorbed) < 1e-12, `closure at tau ${tau} c ${c}`);
      assert.ok(r.reflected >= 0 && r.surface >= 0 && r.absorbed >= -1e-12);
    }
    for (const key of ['reflected', 'surface', 'absorbed', 'down', 'direct']) {
      assert.ok(Math.abs(one[key] - many[key]) < 1e-10 && Math.abs(one[key] - most[key]) < 1e-10, `${key} of one layer ${one[key]}, seven ${many[key]}, forty ${most[key]}`);
    }
  }
});

test('the folded near-infrared table has its 32 g-points, weights that make each band and the beam whole, and the direct beam\'s absorption of the full table to a per cent', () => {
  assert.equal(NEAR_INFRARED_POINTS, 32);
  assert.ok(Math.abs(NEAR_INFRARED_SHARE - NEAR_INFRARED_TABLE.bands.reduce((s, b) => s + b.share, 0)) < 1e-12);
  for (const band of LAYERED_TABLE.bands) assert.ok(Math.abs(band.points.reduce((s, p) => s + p.weight, 0) - 1) < 1e-4, `band ${band.band}`);
  let worst = 0, sum = 0, count = 0;
  for (const column of [TROPICAL, MIDLATITUDE]) for (const m of [1, 2]) {
    const air = { K, vapour: new Float64Array(K), dry: new Float64Array(K), p: new Float64Array(K), T: new Float64Array(K), e: new Float64Array(K) };
    for (let k = 0; k < K; k++) {
      const mass = column.ps * dSigma[k] / GRAVITY, q = column.q[k];
      air.vapour[k] = q * mass; air.dry[k] = mass * (1 - q); air.p[k] = column.ps * sigmaMid[k]; air.T[k] = column.T[k]; air.e[k] = q * air.p[k] / (0.622 + 0.378 * q);
    }
    const mixing = { co2: meanMixingRatio(column, 'co2') * NEAR_INFRARED_MOLAR.co2, ch4: meanMixingRatio(column, 'ch4') * NEAR_INFRARED_MOLAR.ch4, o2: 0.2095 * NEAR_INFRARED_MOLAR.o2 };
    const taken = new Float64Array(K), upward = new Float64Array(K);
    nearInfraredBands(air, m, 5 / 3, mixing, 1, taken, upward);
    let folded = 0;
    for (const band of LAYERED_TABLE.bands) {
      const source = NEAR_INFRARED_TABLE.bands.find((b) => b.band === band.band);
      const base = band.gases.map((name) => {
        const par = source.gases[name];
        return par.k.map((_, j) => Float64Array.from({ length: K }, (__, k) => par.k[j] * (name === 'h2o' ? air.vapour[k] : mixing[name] * air.dry[k]) * Math.exp(par.n[j] * Math.log(Math.sqrt(air.p[k] ** 2 + par.pD ** 2) / LAYERED_TABLE.pRef) + par.t[j] * (air.T[k] - 250))));
      });
      const par = source.gases.h2o, self = Float64Array.from({ length: K }, (_, k) => (par.self > 0 ? par.self * air.vapour[k] * air.e[k] / LAYERED_TABLE.pRef * Math.exp(par.selfT * Math.log(296 / air.T[k])) : 0));
      for (const point of band.points) {
        let total = 0;
        for (let k = 0; k < K; k++) { let d = self[k]; band.gases.forEach((_, i) => point.mix[i].forEach((x, j) => { d += x * base[i][j][k]; })); total += point.scale * d; }
        folded += band.share * point.weight * (1 - Math.exp(-m * total));
      }
    }
    const miss = Math.abs(folded - taken[K - 1]) * S0;
    worst = Math.max(worst, miss); sum += miss; count++;
  }
  console.log(`direct-beam near-infrared absorption of the folded table against the full table: mean miss ${(sum / count).toFixed(2)} W/m2, worst ${worst.toFixed(2)} W/m2 of ${S0}`);
  assert.ok(sum / count < 2 && worst < 5);
});

test('under layered the absorbed and the reflected sunlight make the beam, the layers\' heating makes the atmosphere\'s absorption, and the clear sky is the sequential scheme\'s', () => {
  let worstClosure = 0, worstHeating = 0, columns = 0;
  for (const [name, column] of [['tropical', TROPICAL], ['midlatitude', MIDLATITUDE]]) {
    for (const [label, clouds] of CASES) {
      for (const zenith of [10, 30, 60, 80]) {
        const { sequential, layered } = columnBudgets(column, { zenith, clouds, options: { seaDropletRadius: 10 } });
        for (const b of [sequential, layered]) {
          worstClosure = Math.max(worstClosure, Math.abs(b.absorbedSolar + b.reflectedSolar - b.beam) / b.beam);
          worstHeating = Math.max(worstHeating, Math.abs(b.heating - b.atmosphereSolar) / b.beam);
          assert.ok(b.reflectedSolar >= 0 && b.reflectedSolar <= b.beam, `${name} ${label} at ${zenith} degrees: reflected ${b.reflectedSolar} of ${b.beam}`);
          assert.ok(b.surfaceSolar >= 0 && b.atmosphereSolar >= 0, `${name} ${label} at ${zenith} degrees: surface ${b.surfaceSolar}, atmosphere ${b.atmosphereSolar}`);
        }
        assert.equal(layered.clearAbsorbedSolar, sequential.clearAbsorbedSolar, `${name} ${label} at ${zenith} degrees: clear-sky absorption`);
        columns++;
      }
    }
  }
  console.log(`${columns} overcast columns: absorbed plus reflected differs from the beam by at most ${worstClosure.toExponential(1)} of it and the layers' heating from the atmosphere's absorption by ${worstHeating.toExponential(1)}`);
  assert.ok(worstClosure < 1e-12 && worstHeating < 1e-12);
});

test('under thick cloud the surface takes sunlight at every sun, from overhead to the horizon, and the budget stays closed', () => {
  let least = Infinity, count = 0;
  for (const column of [TROPICAL, MIDLATITUDE]) for (const zenith of [1, 30, 60, 80, 85, 88.5]) for (const clouds of [[[75000, 1.5]], [[60000, 1.4], [85000, 2]], [[20000, 0.3], [92500, 0.4]], [[98000, 0.5]]]) {
    const { layered } = columnBudgets(column, { zenith, clouds, options: { seaDropletRadius: 10 } });
    least = Math.min(least, layered.surfaceSolar / layered.beam, layered.surfaceShortwave / layered.beam, layered.surfaceDirect / layered.beam);
    assert.ok(Math.abs(layered.absorbedSolar + layered.reflectedSolar - layered.beam) < 1e-12 * layered.beam, `closure at ${zenith} degrees under ${JSON.stringify(clouds)}`);
    assert.ok(layered.reflectedSolar >= 0 && layered.reflectedSolar <= layered.beam, `reflected ${layered.reflectedSolar} of ${layered.beam}`);
    count++;
  }
  console.log(`${count} thick-cloud columns: the least surface absorption, downward flux or direct flux is ${(100 * least).toFixed(3)} % of the beam`);
  assert.ok(least >= 0);
});

test('a column without cloud is the sequential scheme\'s bit for bit, and a vanishing cloud approaches it', () => {
  for (const column of [TROPICAL, MIDLATITUDE]) for (const zenith of [20, 70]) {
    const { sequential, layered } = columnBudgets(column, { zenith, clouds: [] });
    for (const key of Object.keys(sequential)) assert.equal(layered[key], sequential[key], `${key} without cloud`);
    const faint = columnBudgets(column, { zenith, clouds: [[30000, 1e-12]], options: { nearInfraredRayleigh: 0 } });
    for (const key of ['absorbedSolar', 'atmosphereSolar', 'surfaceSolar', 'reflectedSolar']) assert.ok(Math.abs(faint.layered[key] - faint.sequential[key]) < 1e-6 * faint.sequential.beam, `${key} under a cloud of 1e-12 kg/m2: ${faint.layered[key]} against ${faint.sequential[key]}`);
  }
});

test('layered shields the gases under a thick high cloud, lets a stratocumulus take more than the gray cloud does, and reflects what the gray cloud absorbed in the visible', () => {
  const rows = [];
  for (const zenith of [30, 60]) {
    const anvil = columnBudgets(TROPICAL, { zenith, clouds: [[25000, 0.3]] }), deck = columnBudgets(TROPICAL, { zenith, clouds: [[92500, 0.1]], options: { seaDropletRadius: 10 } });
    const clear = columnBudgets(TROPICAL, { zenith, clouds: [] }).layered;
    const pct = (x, b) => (100 * x / b.beam).toFixed(1);
    rows.push(`${zenith} degrees: anvil atmosphere ${pct(clear.atmosphereSolar, clear)} clear, ${pct(anvil.sequential.atmosphereSolar, clear)} sequential, ${pct(anvil.layered.atmosphereSolar, clear)} layered %; reflected ${pct(anvil.sequential.reflectedSolar, clear)} to ${pct(anvil.layered.reflectedSolar, clear)} %; stratocumulus atmosphere ${pct(deck.sequential.atmosphereSolar, clear)} to ${pct(deck.layered.atmosphereSolar, clear)} %, surface ${pct(deck.sequential.surfaceSolar, clear)} to ${pct(deck.layered.surfaceSolar, clear)} %`);
    assert.ok(anvil.layered.atmosphereSolar < anvil.sequential.atmosphereSolar - 0.03 * clear.beam, 'the anvil shields the gases below it');
    assert.ok(anvil.layered.atmosphereSolar > clear.atmosphereSolar - 0.12 * clear.beam, 'the anvil\'s own absorption makes up part of what it shields');
    assert.ok(anvil.layered.reflectedSolar > anvil.sequential.reflectedSolar, 'the anvil reflects more of the sunlight');
    assert.ok(deck.layered.atmosphereSolar > deck.sequential.atmosphereSolar && deck.layered.atmosphereSolar > clear.atmosphereSolar, 'the stratocumulus adds to the atmosphere\'s absorption');
    assert.ok(deck.layered.surfaceSolar < deck.sequential.surfaceSolar, 'the surface takes less under the stratocumulus');
  }
  console.log(rows.join('\n'));
});

test('under layered the cloudy N = 16 model\'s absorbed plus reflected sunlight makes the insolation of every step, held or not, and the layers act in its columns', () => {
  for (const every of [1, 4]) {
    const model = createModel(new Grid(16), { ocean: false, radiation: { stratus: true, clearSkyPass: true, radiationEvery: every, shortwave: 'layered' }, surface: { exchange: 'roughness' } });
    initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
    const radiation = model.radiation, cells = model.mesh.nCells, dt = 900;
    model.time = 150 * DAY;
    for (let n = 0; n < 40; n++) model.step(dt);
    const absorbed = new Float64Array(cells), reflected = new Float64Array(cells), insolation = new Float64Array(cells);
    let worst = 0, lit = 0;
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
    let acting = 0;
    for (let i = 0; i < cells; i++) {
      radiation.column(i, pi[i], theta, surfaceT[i], 5, undefined, undefined, q[bottom + i], q, qc, 0.1, 0.1, 1, 0, 0, 0, 0, 0);
      if (radiation.budget.cloudSolar > 0) acting++;
    }
    console.log(`radiationEvery ${every}: ${every === 1 ? 'absorbed plus reflected differs from the insolation of the step' : 'the held call\'s absorbed plus reflected differs from its beam'} by at most ${worst.toExponential(1)} W/m² over ${lit} lit cell-steps; ${acting} of ${cells} columns hold absorbing cloud at the end`);
    assert.ok(lit > 2 * cells);
    assert.ok(worst < 1e-8, `absorbed plus reflected misses the sunlight by ${worst} W/m²`);
    assert.ok(acting > 0, 'no column holds cloud');
  }
});
