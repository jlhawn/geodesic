import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore, sigmaInterfaces, GRAVITY } from '../js/dynamics/sigmaCore.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { createRadiation, DAY, SOLAR_CONSTANT } from '../js/physics/radiation.module.js';
import { OZONE_CM_ATM } from '../js/physics/shortwaveGases.module.js';
import { readTopography, checkGeography, FRESH_GEOGRAPHY, geographyOptions } from '../js/geography.module.js';
import { BENCHMARK, modelColumn, meanMixingRatio } from '../scripts/standardAtmospheres.mjs';
import { syntheticAerosol, fillAerosol } from './helpers/aerosol.mjs';

/*
 * The shortwave under radiation option aerosol 'mac' keeps the sunlight: in the columns of the standard
 * atmospheres, clear and cloudy, sequential and layered, over a dark sea and a bright land, at a high and
 * a low sun; in the N = 16 model's steps, held or not; and the aerosol's own absorption is the layers'
 * heating from it. The layered clear region's aerosol against the sequential clear column's is recorded.
 */
const levels = sigmaInterfaces('bl48'), mesh = buildMesh(new Grid(16)), core = createSigmaCore(mesh, { levels });
const { K, C, dSigma, sigmaMid } = core.diagnostics;
const S0 = 1360.85;
const COLUMNS = { tropical: modelColumn(BENCHMARK.atmospheres.TROP, levels), midlatitude: modelColumn(BENCHMARK.atmospheres.MLS, levels) };
const EQUATORIAL = Array.from({ length: C }, (_, i) => i).reduce((best, i) => (Math.abs(mesh.latCell[i]) < Math.abs(mesh.latCell[best]) ? i : best), 0);
const GLOBAL_MEANS = [0.083, 0.059, 0.0703, 0.0268, 0.644, 0.76, 2.2, 1.4], HEAVY = [0.3, 0.5, 0.09, 0.03, 0.62, 0.76, 3.5, 2.6];

function gasesOf(column) {
  const ozone = Array.from(column.o3, (x, k) => x * (levels[k + 1] - levels[k]) * column.ps / GRAVITY / OZONE_CM_ATM);
  return { carbonDioxide: meanMixingRatio(column, 'co2'), methane: meanMixingRatio(column, 'ch4'), nitrousOxide: meanMixingRatio(column, 'n2o'), ozoneProfile: Float64Array.from(ozone) };
}
const layerAt = (column, p) => { let k = 0; for (let j = 1; j < K; j++) if (Math.abs(sigmaMid[j] * column.ps - p) < Math.abs(sigmaMid[k] * column.ps - p)) k = j; return k; };
function sunAt(radiation, zenith) {
  const mu = Math.cos(zenith * Math.PI / 180);
  let noon = 0, highest = -1;
  for (let s = 0; s < 1440; s++) { radiation.setTime(60 * s); const c = radiation.cosZenith(EQUATORIAL); if (c > highest) { highest = c; noon = 60 * s; } }
  assert.ok(highest > mu, `the sun reaches cosine ${highest} at noon, short of ${mu}`);
  let low = noon - 6 * 3600, high = noon;
  for (let n = 0; n < 50; n++) { const middle = 0.5 * (low + high); radiation.setTime(middle); if (radiation.cosZenith(EQUATORIAL) < mu) low = middle; else high = middle; }
  return 0.5 * (low + high);
}

/*
 * One column of a standard atmosphere under `shortwave`, with `aerosol` ('none', 'mac' with the modes' eight
 * fields, or 'flat' with today's 0.12) over a surface of albedo `albedo`, land or sea, at `zenith`: the
 * radiation's budget, the layers' shortwave heating, and the aerosol's profiles.
 */
function column(atmosphere, { shortwave, zenith, clouds = [], albedo = 0.06, land = false, aerosol = 'none', fields = GLOBAL_MEANS, options = {} }) {
  const radiation = createRadiation(mesh, core, {
    aerosol: 'flat', landAerosol: 0, seaAerosol: 0, clearSkyPass: true, cloudCover: 'overcast', cloudColumns: false, solarConstant: S0, shortwave, ...gasesOf(atmosphere), ...(land ? { land: new Uint8Array(C).fill(1) } : {}),
    ...(aerosol === 'mac' ? { aerosol: 'mac', aerosolFields: fillAerosol(new Float32Array(8 * C), fields) } : aerosol === 'flat' ? { landAerosol: 0.12, seaAerosol: 0.12 } : {}), ...options,
  });
  radiation.setTime(sunAt(radiation, zenith));
  const pi = new Float64Array(C).fill(atmosphere.ps), theta = new Float64Array(K * C).fill(280), q = new Float64Array(K * C), qc = new Float64Array(K * C);
  core.diagnose(pi, theta);
  const { exnerLayer } = core.diagnostics;
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) { theta[k * C + i] = atmosphere.T[k] / exnerLayer[k * C + i]; q[k * C + i] = atmosphere.q[k]; }
  for (const [p, path] of clouds) { const k = layerAt(atmosphere, p); for (let i = 0; i < C; i++) qc[k * C + i] = path * GRAVITY / (atmosphere.ps * dSigma[k]); }
  core.diagnose(pi, theta, q, qc);
  radiation.column(EQUATORIAL, atmosphere.ps, theta, atmosphere.Ts, 5, undefined, undefined, q[(K - 1) * C + EQUATORIAL], q, qc, albedo, albedo, 0, 0, 0, 0, 0, 0);
  const heating = Float64Array.from({ length: K }, (_, k) => radiation.layerFlux[k] - radiation.longwave[k * C + EQUATORIAL]);
  return { ...radiation.budget, beam: radiation.budget.insolation, heating, total: heating.reduce((a, b) => a + b, 0), profiles: radiation.aerosolProfiles && { visible: Float64Array.from(radiation.aerosolProfiles.visible), nearInfrared: Float64Array.from(radiation.aerosolProfiles.nearInfrared) }, mu: radiation.cosZenith(EQUATORIAL) };
}

const CLOUDS = { clear: [], 'a thin cirrus': [[20000, 0.005]], 'a stratocumulus': [[92500, 0.1]], 'two layers': [[25000, 0.1], [85000, 0.1]], 'a thick anvil': [[25000, 0.3]] };

test('with aerosol \'mac\' the absorbed and the reflected sunlight make the beam and the layers\' heating the atmosphere\'s absorption, clear and cloudy, in either scheme, over sea and land, at a high and a low sun', () => {
  let columns = 0, worstClosure = 0, worstHeating = 0, least = Infinity;
  for (const [name, atmosphere] of Object.entries(COLUMNS)) for (const shortwave of ['sequential', 'layered']) for (const [label, clouds] of Object.entries(CLOUDS)) for (const [surface, land, albedo] of [['sea', false, 0.06], ['land', true, 0.3]]) for (const zenith of [10, 45, 80]) for (const fields of [GLOBAL_MEANS, HEAVY]) {
    const b = column(atmosphere, { shortwave, zenith, clouds, albedo, land, aerosol: 'mac', fields, options: { seaDropletRadius: 10 } });
    worstClosure = Math.max(worstClosure, Math.abs(b.absorbedSolar + b.reflectedSolar - b.beam) / b.beam);
    worstHeating = Math.max(worstHeating, Math.abs(b.total - b.atmosphereSolar) / b.beam);
    least = Math.min(least, b.reflectedSolar / b.beam, b.surfaceSolar / b.beam, b.atmosphereSolar / b.beam);
    assert.ok(b.aerosolSolar > 0 && b.aerosolSolar < b.atmosphereSolar, `${name} ${shortwave} ${label} over ${surface} at ${zenith} degrees: the aerosol absorbs ${b.aerosolSolar} of the atmosphere's ${b.atmosphereSolar}`);
    assert.ok(b.reflectedSolar >= 0 && b.reflectedSolar <= b.beam && b.surfaceSolar >= 0, `${name} ${shortwave} ${label} over ${surface} at ${zenith} degrees: reflected ${b.reflectedSolar}, surface ${b.surfaceSolar}`);
    columns++;
  }
  console.log(`${columns} columns with the MACv2 modes: absorbed plus reflected differs from the beam by at most ${worstClosure.toExponential(1)} of it, the layers' heating from the atmosphere's absorption by ${worstHeating.toExponential(1)}`);
  assert.ok(worstClosure < 1e-12 && worstHeating < 1e-12, `${worstClosure}, ${worstHeating}`);
});

test('the aerosol\'s absorption is positive and is the layers\' heating from it, which lies where its mean height says', () => {
  const rows = [], low = [0.1, 0.1, 0.07, 0.07, 0.65, 0.75, 0.6, 0.6];
  for (const shortwave of ['sequential', 'layered']) for (const zenith of [10, 45, 75]) for (const [name, atmosphere] of Object.entries(COLUMNS)) {
    const clouds = shortwave === 'layered' ? [[30000, 1e-6]] : [], options = { upwardAbsorption: false };
    const without = column(atmosphere, { shortwave, zenith, clouds, albedo: 0, options }), mac = column(atmosphere, { shortwave, zenith, clouds, albedo: 0, aerosol: 'mac', fields: low, options });
    let change = 0, above = 0;
    for (let k = 0; k < K; k++) { const d = mac.heating[k] - without.heating[k]; change += d; if (sigmaMid[k] * atmosphere.ps < atmosphere.ps - 300e2) above += d; }
    rows.push(`${shortwave} ${name} at ${zenith} degrees: the aerosol absorbs ${mac.aerosolSolar.toFixed(2)} W/m2 and the layers gain ${change.toFixed(2)}, ${(above / change * 100).toFixed(1)} % of it above 3 km`);
    assert.ok(mac.aerosolSolar > 0.5 && without.aerosolSolar === 0, 'the aerosol absorbs');
    assert.ok(Math.abs(change - mac.aerosolSolar) < (shortwave === 'sequential' ? 1e-9 : 1e-4) * mac.aerosolSolar, `the layers gain ${change} W/m2 where the aerosol absorbs ${mac.aerosolSolar}`);
    assert.ok(Math.abs(above) < 0.1 * change, `${above} of ${change} W/m2 above 3 km, from aerosol with a mean height of 0.6 km`);
    assert.ok(mac.absorbedSolar < without.absorbedSolar && mac.surfaceSolar < without.surfaceSolar, 'the aerosol reflects sunlight over a black surface and shades it');
  }
  console.log(rows.join('\n'));
});

test('no aerosol depth, and aerosolScale 0, are the model without aerosol bit for bit, and the depth scales the absorption', () => {
  const same = ['absorbedSolar', 'atmosphereSolar', 'reflectedSolar', 'surfaceSolar', 'surfaceShortwave', 'surfaceDirect', 'clearAbsorbedSolar', 'clearSurfaceSolar', 'aerosolSolar', 'cloudSolar', 'vaporSolar'];
  for (const shortwave of ['sequential', 'layered']) for (const clouds of [[], [[92500, 0.1]]]) for (const zenith of [20, 70]) {
    const none = column(COLUMNS.tropical, { shortwave, zenith, clouds }), zero = column(COLUMNS.tropical, { shortwave, zenith, clouds, aerosol: 'mac', fields: GLOBAL_MEANS.map((v, f) => (f < 2 ? 0 : v)) });
    const scaled = column(COLUMNS.tropical, { shortwave, zenith, clouds, aerosol: 'mac', options: { aerosolScale: 0 } });
    for (const key of same) { assert.equal(zero[key], none[key], `${key} without depth, ${shortwave} at ${zenith} degrees`); assert.equal(scaled[key], none[key], `${key} at scale 0`); }
    assert.deepEqual(Array.from(zero.heating), Array.from(none.heating));
    const one = column(COLUMNS.tropical, { shortwave, zenith, clouds, aerosol: 'mac' }), twice = column(COLUMNS.tropical, { shortwave, zenith, clouds, aerosol: 'mac', options: { aerosolScale: 2 } }), half = column(COLUMNS.tropical, { shortwave, zenith, clouds, aerosol: 'mac', options: { aerosolScale: 0.5 } });
    assert.ok(twice.aerosolSolar / one.aerosolSolar > 1.8 && twice.aerosolSolar / one.aerosolSolar < 2.02 && one.aerosolSolar / half.aerosolSolar > 1.8 && one.aerosolSolar / half.aerosolSolar < 2.02, `absorption ${half.aerosolSolar}, ${one.aerosolSolar}, ${twice.aerosolSolar}`);
    if (clouds.length === 0) assert.ok(twice.surfaceSolar < one.surfaceSolar && one.surfaceSolar < none.surfaceSolar && twice.absorbedSolar < one.absorbedSolar && one.absorbedSolar < none.absorbedSolar, `the surface and the top take less as the depth rises: ${none.surfaceSolar}, ${one.surfaceSolar}, ${twice.surfaceSolar} at ${zenith} degrees`);
  }
});

test('the layered clear region and the sequential clear column take the MACv2 aerosol alike in the atmosphere and, with the sun 45 degrees high, at the top; the rest is recorded', () => {
  const rows = [], ratios = new Map();
  for (const [name, atmosphere] of Object.entries(COLUMNS)) for (const albedo of [0.06, 0.2]) for (const zenith of [2, 30, 45, 60, 70]) {
    const thin = [[30000, 1e-6]];
    const seq = [column(atmosphere, { shortwave: 'sequential', zenith, albedo }), column(atmosphere, { shortwave: 'sequential', zenith, albedo, aerosol: 'mac' })];
    const lay = [column(atmosphere, { shortwave: 'layered', zenith, albedo, clouds: thin }), column(atmosphere, { shortwave: 'layered', zenith, albedo, clouds: thin, aerosol: 'mac' })];
    const top = [seq[1].clearAbsorbedSolar - seq[0].clearAbsorbedSolar, lay[1].layeredClearAbsorbed - lay[0].layeredClearAbsorbed], surface = [seq[1].clearSurfaceSolar - seq[0].clearSurfaceSolar, lay[1].layeredClearSurface - lay[0].layeredClearSurface];
    const atm = [top[0] - surface[0], top[1] - surface[1]];
    ratios.set(`${name} ${albedo} ${zenith}`, { top: top[1] / top[0], atmosphere: atm[1] / atm[0], surface: surface[1] / surface[0] });
    rows.push(`${name.padEnd(12)} albedo ${albedo} zenith ${String(zenith).padStart(2)}: top ${top[0].toFixed(2).padStart(7)} ${top[1].toFixed(2).padStart(7)} (${(top[1] / top[0]).toFixed(2)}) atmosphere ${atm[0].toFixed(2).padStart(6)} ${atm[1].toFixed(2).padStart(6)} (${(atm[1] / atm[0]).toFixed(2)}) surface ${surface[0].toFixed(2).padStart(7)} ${surface[1].toFixed(2).padStart(7)} (${(surface[1] / surface[0]).toFixed(2)})`);
    assert.ok(atm[0] > 5 && atm[1] > 3 && top[0] < 0 + 5, `the aerosol absorbs in both: ${atm}`);
  }
  console.log(`increments of the absorbed sunlight with the MACv2 global mean modes, W/m2, sequential and layered (ratio):\n${rows.join('\n')}`);
  for (const [key, { top, atmosphere, surface }] of ratios) {
    const zenith = Number(key.split(' ')[2]);
    if (zenith <= 60) assert.ok(atmosphere > 0.75 && atmosphere < 1.35, `${key}: the atmosphere's increments differ by ${atmosphere}`);
    if (zenith === 45) assert.ok(Math.abs(top - 1) < 0.1 && Math.abs(surface - 1) < 0.1, `${key}: top ${top}, surface ${surface}`);
  }
});

const held = (every) => (every === 1 ? 'every step' : `every ${every}th step`);
const geography = checkGeography({ ...FRESH_GEOGRAPHY }), sheets = geographyOptions(geography);
for (const [place, extra] of [['an aquaplanet', {}], ['the land and sea of the fresh geography', { topography: readTopography(geography.topography), geography: sheets.geography, land: sheets.land }]]) for (const [shortwave, every] of [['layered', 1], ['layered', 4], ['sequential', 1]]) {
  test(`over ${place} the N = 16 model's absorbed plus reflected sunlight makes the insolation of every step under aerosol 'mac' and shortwave '${shortwave}', radiation ${held(every)}, across a change of the model day`, () => {
    const { bytes } = syntheticAerosol({ nLon: 36, nLat: 18 });
    const model = createModel(new Grid(16), { ocean: false, ...extra, radiation: { stratus: true, clearSkyPass: true, radiationEvery: every, shortwave, aerosol: 'mac', aerosolData: bytes }, surface: { exchange: 'roughness' } });
    initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
    const radiation = model.radiation, cells = model.mesh.nCells, dt = 1350;
    model.time = 150 * DAY - 9 * dt;
    for (let n = 0; n < 6; n++) model.step(dt);
    const days = new Set(), absorbed = new Float64Array(cells), reflected = new Float64Array(cells), insolation = new Float64Array(cells);
    let worst = 0, lit = 0;
    for (let n = 0; n < 6; n++) {
      const time = model.time;
      model.step(dt);
      days.add(model.aerosol.day);
      radiation.setTime(time);
      for (let i = 0; i < cells; i++) {
        const beam = radiation.summed.insolation[i] - insolation[i], top = radiation.summed.absorbedSolar[i] - absorbed[i], lost = radiation.summed.reflectedSolar[i] - reflected[i];
        if (every === 1) { worst = Math.max(worst, Math.abs(top + lost - beam)); if (beam > 0) lit++; }
        else if (radiation.held.mu[i] > 0) { worst = Math.max(worst, Math.abs(radiation.held.top[i] + radiation.held.reflected[i] - SOLAR_CONSTANT * radiation.held.mu[i])); lit++; }
      }
      absorbed.set(radiation.summed.absorbedSolar); reflected.set(radiation.summed.reflectedSolar); insolation.set(radiation.summed.insolation);
      radiation.setTime(model.time);
    }
    const [pi, theta, , surfaceT, q, qc] = model.state, bottom = (model.core.K - 1) * cells;
    radiation.setTime(model.time);
    let aerosolSum = 0, litNow = 0;
    for (let i = 0; i < cells; i++) {
      if (!(radiation.cosZenith(i) > 0)) continue;
      radiation.column(i, pi[i], theta, surfaceT[i], 5, undefined, undefined, q[bottom + i], q, qc, 0.1, 0.1, 1, 0, 0, 0, 0, 0);
      aerosolSum += radiation.budget.aerosolSolar; litNow++;
    }
    console.log(`${place}, ${shortwave}, radiation ${held(every)}: absorbed plus reflected differs from the sunlight by at most ${worst.toExponential(1)} W/m² over ${lit} lit cell-steps; days ${[...days].join(', ')}; the aerosol absorbs a mean ${(aerosolSum / litNow).toFixed(2)} W/m² of ${litNow} lit columns`);
    assert.ok(days.size >= 2, `the model day changed ${[...days]}`);
    assert.ok(lit > cells && worst < 1e-8, `absorbed plus reflected misses the sunlight by ${worst} W/m² in ${lit} cell-steps`);
    assert.ok(litNow > 0.3 * cells && aerosolSum / litNow > 1, 'the aerosol absorbs');
  });
}
