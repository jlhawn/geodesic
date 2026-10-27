import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';
import { createLandSurface, soilColumn, airCycle, SOIL_LAG } from '../js/physics/land.module.js';
import { MELTING_POINT } from '../js/physics/ice.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { regridLand } from '../js/physics/regrid.module.js';
import { recordedOptions } from '../scripts/stateOptions.mjs';

const N = 16, YEAR = 365 * 86400, DT = 1350, LF = 3.34e5;
const mesh = createModel(new Grid(N), { physics: false }).mesh;
const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 || lat < -1.15 ? 400 : -4000));
const geography = createGeography(mesh, topography, { landBridges: {}, seaStraits: {} });
const C = mesh.nCells;
const firnCell = geography.land.findIndex((l, i) => l && geography.iceSheet[i] > 0);
const northCell = geography.land.findIndex((l, i) => l && !geography.iceSheet[i] && mesh.latCell[i] > 0.5 && mesh.latCell[i] < 0.9);
const southCell = geography.land.findIndex((l, i) => l && !geography.iceSheet[i] && mesh.latCell[i] < -0.5 && mesh.latCell[i] > -0.9);
const landCell = geography.land.findIndex((l, i) => l && !geography.iceSheet[i]);

const soilClass = { heat: 2.19e6, conductivity: 0.7 };
const firnClass = { heat: 2.05e6 * 400 / 917, conductivity: 2.2 * (400 / 917) ** 1.88 };
const DEPTHS = [0.07, 0.21, 0.72, 1.89];
const middles = DEPTHS.map((d, k) => DEPTHS.slice(0, k).reduce((a, b) => a + b, 0) + d / 2);

const build = (options = {}, geo = geography) => createLandSurface(mesh, geo, { soil: 'layers', ...options });
const fresh = (options = {}, geo = geography) => { const land = build(options, geo); land.initialize(); land.snow.fill(0); return land; };
const columnOf = (land, i) => (geography.iceSheet[i] > 0 ? land.column.firn : land.column.soil);
const heatOf = (land, surfaceT, i, reference = 270) => {
  const { capacity } = columnOf(land, i);
  return capacity.reduce((sum, c, k) => sum + c * ((k === 0 ? surfaceT[i] : land.deep[k - 1][i]) - reference), 0);
};
const put = (land, surfaceT, i, temperatures) => { surfaceT[i] = temperatures[0]; temperatures.slice(1).forEach((t, k) => { land.deep[k][i] = t; }); };
const spread = (land, surfaceT, i) => [surfaceT[i], ...land.deep.map((layer) => layer[i])];

test('the layers take the IFS thicknesses and heat capacity, the conductances join their mid-depths, and a firn column takes the density\'s properties', () => {
  const column = soilColumn({ soil: 'layers' });
  assert.deepEqual(column.depths, DEPTHS);
  DEPTHS.forEach((d, k) => assert.ok(Math.abs(column.soil.capacity[k] - 2.19e6 * d) < 1e-6 * column.soil.capacity[k]));
  [0.14, 0.465, 1.305].forEach((dz, k) => assert.ok(Math.abs(column.soil.conductance[k] - 0.7 / dz) < 1e-12, `soil g${k + 1} ${column.soil.conductance[k]}`));
  assert.ok(Math.abs(column.soil.capacity[0] - 1.533e5) < 1, 'the skin is 1.533e5 J/m²/K');
  assert.ok(Math.abs(column.firn.heatCapacity - 2.05e6 * 400 / 917) < 1e-6 && Math.abs(column.firn.conductivity - 2.2 * (400 / 917) ** 1.88) < 1e-12);
  assert.ok(column.firn.heatCapacity > 8.9e5 && column.firn.heatCapacity < 9e5 && column.firn.conductivity > 0.46 && column.firn.conductivity < 0.47);
  assert.ok(Math.abs(column.firn.conductance[0] - column.firn.conductivity / 0.14) < 1e-12);
  const tuned = soilColumn({ soil: 'layers', soilDepths: [0.1, 0.3], soilHeatCapacity: 2e6, soilConductivity: 1, firnDensity: 300 });
  assert.deepEqual(tuned.depths, [0.1, 0.3]);
  assert.ok(Math.abs(tuned.soil.conductance[0] - 1 / 0.2) < 1e-12 && Math.abs(tuned.firn.heatCapacity - 2.05e6 * 300 / 917) < 1e-6);
  assert.equal(soilColumn({ soil: 'slab' }), null);
  assert.equal(soilColumn({}), null);
  assert.throws(() => soilColumn({ soil: 'layers', heatCapacity: 1e6 }), /heatCapacity is the slab's/);
  assert.throws(() => soilColumn({ soil: 'deep' }), /'slab' or 'layers'/);
  assert.throws(() => soilColumn({ soil: 'layers', soilDepths: [] }), /layer thicknesses/);
  assert.throws(() => soilColumn({ soil: 'layers', soilDepths: [0.1, 0] }), /layer thicknesses/);
  assert.throws(() => soilColumn({ soil: 'layers', soilConductivity: 0 }), /above 0/);
  assert.throws(() => soilColumn({ soil: 'layers', firnDensity: 1000 }), /firnDensity/);
  assert.throws(() => createLandSurface(mesh, geography, { soil: 'layers', heatCapacity: 1.5e5 }), /heatCapacity is the slab's/);
});

test('the land energy row: the layers\' heat content changes by the surface flux times the step, less the fusion heat of the snow and glacier ice that went, in soil and firn columns', () => {
  const worst = { constant: 0, alternating: 0, melting: 0, glacier: 0 };
  const advance = (land, surfaceT, flux, i, q, dt, evaporation = 0) => {
    const before = heatOf(land, surfaceT, i, 0), snow = land.snow[i], lost = land.iceLoss[i];
    flux[i] = q;
    land.update(i, surfaceT, flux, evaporation, dt);
    const after = heatOf(land, surfaceT, i, 0), expected = q * dt + LF * (land.snow[i] - snow) - LF * (land.iceLoss[i] - lost);
    return Math.abs(after - before - expected) / (Math.abs(before) + Math.abs(after));
  };
  for (const i of [landCell, firnCell]) {
    const land = fresh(), surfaceT = new Float64Array(C), flux = new Float64Array(C);
    put(land, surfaceT, i, [258, 262, 266, 270]);
    for (let n = 0; n < 400; n++) worst.constant = Math.max(worst.constant, advance(land, surfaceT, flux, i, 120, DT));
    for (let n = 0; n < 400; n++) worst.alternating = Math.max(worst.alternating, advance(land, surfaceT, flux, i, 200 * Math.sin(n / 7) * (n % 3 ? 1 : -2), DT * (1 + (n % 5))));
    const warm = fresh(), skin = new Float64Array(C);
    put(warm, skin, i, [272.5, 271, 270, 270]);
    warm.snow[i] = 8;
    for (let n = 0; n < 60; n++) worst.melting = Math.max(worst.melting, advance(warm, skin, flux, i, 400, DT, n % 4 === 0 ? 2e-5 : 0));
    assert.ok(warm.snow[i] < 8, `cell ${i}: the snow went`);
    assert.ok(skin[i] <= MELTING_POINT + 1e-9 || warm.snow[i] <= 0, `cell ${i}: the melting skin is held at ${skin[i]}`);
  }
  const ice = fresh({ iceSheetBalance: true }), skin = new Float64Array(C), flux = new Float64Array(C);
  ice.snow[firnCell] = 0;
  put(ice, skin, firnCell, [273, 272.5, 272, 271.5]);
  for (let n = 0; n < 40; n++) worst.glacier = Math.max(worst.glacier, advance(ice, skin, flux, firnCell, 500, DT));
  assert.ok(ice.iceLoss[firnCell] > 0 && skin[firnCell] <= MELTING_POINT + 1e-9, `the glacier lost ${ice.iceLoss[firnCell]} kg/m²`);
  console.log(`land energy row, worst residual over the column's heat content in J/m²: constant ${worst.constant.toExponential(1)}, alternating ${worst.alternating.toExponential(1)}, melting snow ${worst.melting.toExponential(1)}, glacier melt ${worst.glacier.toExponential(1)}`);
  for (const [name, value] of Object.entries(worst)) assert.ok(value < 1e-14, `${name}: ${value}`);
});

test('a uniform column stays uniform without a flux, whatever the step, and a flux warms the layers it reaches, the heat spreading through the column without a loss', () => {
  const land = fresh(), surfaceT = new Float64Array(C).fill(281), flux = new Float64Array(C);
  put(land, surfaceT, landCell, [281, 281, 281, 281]);
  for (const dt of [DT, 86400, 30 * 86400]) land.update(landCell, surfaceT, flux, 0, dt);
  assert.deepEqual(spread(land, surfaceT, landCell), [281, 281, 281, 281]);
  const before = spread(land, surfaceT, landCell), total = heatOf(land, surfaceT, landCell);
  flux[landCell] = 50;
  land.update(landCell, surfaceT, flux, 0, DT);
  const after = spread(land, surfaceT, landCell);
  assert.ok(after[0] - before[0] > after[1] - before[1] && after[1] - before[1] > after[2] - before[2] && after[2] - before[2] > after[3] - before[3] && after[3] > before[3], 'the heat moves down and decays with depth');
  flux[landCell] = 0;
  for (let n = 0; n < 4000; n++) land.update(landCell, surfaceT, flux, 0, 86400);
  const { capacity } = land.column.soil, whole = capacity.reduce((a, b) => a + b, 0), equal = 270 + (total + 50 * DT) / whole;
  for (const t of spread(land, surfaceT, landCell)) assert.ok(Math.abs(t - equal) < 1e-6, `${t} against the mean ${equal}`);
});

// phasors of the stepped column under a sinusoidal flux into the skin, from the backward step's own difference operator
const complex = (re, im = 0) => ({ re, im });
const add = (a, b) => complex(a.re + b.re, a.im + b.im), mul = (a, b) => complex(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
const div = (a, b) => { const d = b.re * b.re + b.im * b.im; return complex((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d); };
function periodicColumn({ capacity, conductance }, period, dt, discrete = true) {
  const omega = 2 * Math.PI / period, s = discrete ? complex((1 - Math.cos(omega * dt)) / dt, Math.sin(omega * dt) / dt) : complex(0, omega);
  const n = capacity.length, load = [];
  let below = complex(0);
  for (let k = n - 1; k >= 1; k--) {
    const stored = add(mul(s, complex(capacity[k])), below);
    below = div(mul(complex(conductance[k - 1]), stored), add(stored, complex(conductance[k - 1])));
    load[k - 1] = below;
  }
  const skin = n > 1 ? add(mul(s, complex(capacity[0])), load[0]) : mul(s, complex(capacity[0]));
  const phasor = [div(complex(1), skin)];
  for (let k = 1; k < n; k++) {
    const stored = add(mul(s, complex(capacity[k])), k < n - 1 ? load[k] : complex(0));
    phasor.push(mul(phasor[k - 1], div(complex(conductance[k - 1]), add(stored, complex(conductance[k - 1])))));
  }
  return { phasor, omega, skin, admittance: n > 1 ? load[0] : complex(0) };
}

test('under a sinusoidal flux the stepped column follows the backward-step solution of its own equations step by step, at periods from 12 hours to a year, in soil and in firn', () => {
  const rows = [];
  for (const [label, i] of [['soil', landCell], ['firn', firnCell]]) {
    for (const period of [43200, 86400, 10 * 86400, 30 * 86400, YEAR]) {
      const land = fresh(), surfaceT = new Float64Array(C), flux = new Float64Array(C), cls = columnOf(land, i);
      const { phasor, omega } = periodicColumn(cls, period, DT), amplitude = 20 / Math.sqrt(period / 86400), mean = 276;
      const turn = (n) => complex(Math.cos(omega * DT * n), Math.sin(omega * DT * n));
      put(land, surfaceT, i, phasor.map((h) => mean + amplitude * mul(h, turn(-1)).re));
      let worst = 0, peak = 0;
      const steps = Math.round(period / DT);
      for (let n = 0; n < steps; n++) {
        flux[i] = amplitude * Math.cos(omega * DT * n);
        land.update(i, surfaceT, flux, 0, DT);
        const phase = turn(n), now = spread(land, surfaceT, i);
        now.forEach((t, k) => { worst = Math.max(worst, Math.abs(t - mean - amplitude * mul(phasor[k], phase).re)); });
        peak = Math.max(peak, Math.abs(now[0] - mean));
      }
      const continuous = periodicColumn(cls, period, DT, false), stepped = periodicColumn(cls, period, DT, true);
      const modulus = Math.hypot(stepped.phasor[0].re, stepped.phasor[0].im) / Math.hypot(continuous.phasor[0].re, continuous.phasor[0].im);
      rows.push(`${label} ${(period / 86400).toFixed(1)} d: skin swing ${peak.toFixed(2)} K, worst departure from the stepped solution ${worst.toExponential(1)} K, stepped over continuous skin response ${modulus.toFixed(4)}`);
      assert.ok(worst < 1e-9, `${label} ${period}: ${worst}`);
      assert.ok(Math.abs(modulus - 1) < (period < 86400 ? 0.1 : period < 5 * 86400 ? 0.05 : 0.01), `${label} ${period / 86400} d: stepped skin response ${modulus.toFixed(4)} of the continuous column's`);
    }
  }
  console.log(rows.join('\n'));
});

test('a column of one layer is the slab with that layer\'s capacity, in soil and in firn', () => {
  const D = 0.2;
  for (const [i, heat] of [[landCell, 2.19e6], [firnCell, firnClass.heat]]) {
    const layered = build({ soilDepths: [D] }), slab = createLandSurface(mesh, geography, { heatCapacity: heat * D });
    layered.initialize(); slab.initialize();
    assert.equal(layered.deepNames.length, 0);
    const a = new Float64Array(C).fill(275), b = new Float64Array(C).fill(275), flux = new Float64Array(C);
    for (const l of [layered, slab]) l.snow[i] = 3;
    for (let n = 0; n < 200; n++) {
      flux[i] = 300 * Math.sin(n / 9) + 40 - 130 * (i === firnCell);
      for (const [l, T] of [[layered, a], [slab, b]]) l.update(i, T, flux, n % 7 === 0 ? 1e-5 : 0, DT, 270, 1e-5);
    }
    assert.ok(Math.abs(a[i] - b[i]) < 1e-9 && Math.abs(layered.snow[i] - slab.snow[i]) < 1e-9 && Math.abs(layered.soil[i] - slab.soil[i]) < 1e-9, `cell ${i}: ${a[i]} against ${b[i]}`);
  }
});

// the annual wave of a homogeneous half-space, written out from its constants
function wave(cell, depth, time, firn) {
  const { mean, amplitude } = airCycle(mesh.latCell[cell], geography.elevation[cell]);
  const kappa = (firn ? firnClass : soilClass).conductivity / (firn ? firnClass : soilClass).heat, omega = 2 * Math.PI / YEAR, d = Math.sqrt(2 * kappa / omega);
  const peak = (mesh.latCell[cell] >= 0 ? YEAR / 4 : 3 * YEAR / 4) + 30 * 86400;
  return MELTING_POINT + mean + amplitude * Math.exp(-depth / d) * Math.cos(omega * (time - peak) - depth / d);
}

test('a land state without the deep layers starts them on the annual wave of a half-space, whose mean and amplitude by depth match the formula in soil and in firn, in either hemisphere', () => {
  assert.equal(SOIL_LAG, 30 * 86400);
  for (const cell of [northCell, southCell, firnCell]) {
    const firn = geography.iceSheet[cell] > 0, land = build(), samples = 48, series = [0, 1, 2].map(() => []);
    for (let s = 0; s < samples; s++) {
      const time = (s + 0.37) * YEAR / samples;
      land.load({ soil: new Float64Array(C).fill(100), snow: new Float64Array(C) }, null, time);
      for (let k = 1; k < 4; k++) {
        const expected = wave(cell, middles[k], time, firn);
        assert.ok(Math.abs(land.deep[k - 1][cell] - expected) < 1e-9, `cell ${cell} layer ${k + 1} at ${time}: ${land.deep[k - 1][cell]} against ${expected}`);
        series[k - 1].push(land.deep[k - 1][cell]);
      }
    }
    const { mean, amplitude } = airCycle(mesh.latCell[cell], geography.elevation[cell]);
    const kappa = (firn ? firnClass : soilClass).conductivity / (firn ? firnClass : soilClass).heat, d = Math.sqrt(2 * kappa / (2 * Math.PI / YEAR));
    series.forEach((values, n) => {
      const k = n + 1, average = values.reduce((a, b) => a + b, 0) / samples;
      let cosine = 0, sine = 0;
      values.forEach((v, s) => { const w = 2 * Math.PI * (s + 0.37) / samples; cosine += 2 * (v - average) * Math.cos(w) / samples; sine += 2 * (v - average) * Math.sin(w) / samples; });
      assert.ok(Math.abs(average - MELTING_POINT - mean) < 1e-9, `layer ${k + 1} mean ${average}`);
      assert.ok(Math.abs(Math.hypot(cosine, sine) - amplitude * Math.exp(-middles[k] / d)) < 1e-9, `layer ${k + 1} amplitude ${Math.hypot(cosine, sine)} against ${amplitude * Math.exp(-middles[k] / d)}`);
      const peakDay = ((Math.atan2(sine, cosine) / (2 * Math.PI) + 1) % 1) * 365, expectedDay = (((mesh.latCell[cell] >= 0 ? 0.25 : 0.75) * 365 + 30 + 365 * middles[k] / d / (2 * Math.PI)) % 365);
      assert.ok(Math.abs(peakDay - expectedDay) < 0.01 || Math.abs(Math.abs(peakDay - expectedDay) - 365) < 0.01, `layer ${k + 1} peaks on day ${peakDay.toFixed(2)}, expected ${expectedDay.toFixed(2)}`);
    });
  }
});

test('initialize starts the layers on the wave at the given time, the layers the state carries are used as saved, and a state with some layers missing or empty starts them all', () => {
  const land = build(), time = 100 * 86400;
  land.initialize(time);
  for (const cell of [northCell, southCell, firnCell]) for (let k = 1; k < 4; k++) assert.ok(Math.abs(land.deep[k - 1][cell] - wave(cell, middles[k], time, geography.iceSheet[cell] > 0)) < 1e-9);
  for (let i = 0; i < C; i++) if (!geography.land[i]) for (const layer of land.deep) assert.equal(layer[i], 0);
  const saved = { soil: new Float64Array(C).fill(100), snow: new Float64Array(C) };
  const arrays = Object.fromEntries(land.deepNames.map((name, k) => [name, Float64Array.from({ length: C }, (_, i) => (geography.land[i] ? 250 + k + (i % 11) : 7))]));
  const used = build();
  used.load({ ...saved, ...arrays }, null, 12345);
  land.deepNames.forEach((name, k) => { for (let i = 0; i < C; i++) assert.equal(used.deep[k][i], geography.land[i] ? 250 + k + (i % 11) : 0, `${name} at ${i}`); });
  const partial = build();
  partial.load({ ...saved, soilT2: arrays.soilT2 }, null, time);
  const started = build();
  started.initialize(time);
  land.deepNames.forEach((name, k) => assert.deepEqual(Array.from(partial.deep[k]), Array.from(started.deep[k]), `${name} of a partial save`));
  const holed = build();
  holed.load({ ...saved, ...arrays, soilT3: Float64Array.from(arrays.soilT3, (v, i) => (i === northCell ? 0 : v)) }, null, time);
  assert.ok(Math.abs(holed.soilT3[northCell] - wave(northCell, middles[2], time, false)) < 1e-9, 'a cell with an empty layer starts on the wave');
  assert.ok(Math.abs(holed.soilT2[northCell] - wave(northCell, middles[1], time, false)) < 1e-9, 'and takes all its layers from it');
  const slab = createLandSurface(mesh, geography, {});
  slab.load({ ...saved, ...arrays }, null, time);
  assert.equal(slab.soilT2, undefined);
  assert.deepEqual(Object.keys(slab.serialize()).filter((key) => key.startsWith('soilT')), []);
});

test('started on the half-space wave and forced by the wave at its top, the column drifts over a simulated year by a third (soil) or half (firn) of what a start at the skin temperature does', () => {
  const report = [], checks = [];
  for (const cell of [northCell, southCell, firnCell]) {
    const firn = geography.iceSheet[cell] > 0, drift = {};
    for (const start of ['wave', 'skin']) {
      const land = fresh(), surfaceT = new Float64Array(C), flux = new Float64Array(C);
      const begin = 200 * 86400, skin0 = wave(cell, middles[0], begin, firn);
      land.load({ soil: new Float64Array(C).fill(100), snow: new Float64Array(C) }, null, begin);
      if (start === 'skin') land.deep.forEach((layer) => { layer[cell] = skin0; });
      surfaceT[cell] = skin0;
      const initial = spread(land, surfaceT, cell), sums = initial.map(() => 0), steps = Math.round(YEAR / DT);
      for (let n = 0; n < steps; n++) {
        const time = begin + n * DT, target = wave(cell, middles[0], time + 0.5 * DT, firn);
        flux[cell] = 50 * (target - surfaceT[cell]);
        land.update(cell, surfaceT, flux, 0, DT);
        spread(land, surfaceT, cell).forEach((t, k) => { sums[k] += t / steps; });
      }
      const final = spread(land, surfaceT, cell), { mean } = airCycle(mesh.latCell[cell], geography.elevation[cell]);
      drift[start] = {
        mean: sums.map((s) => s - MELTING_POINT - mean),
        end: final.map((t, k) => t - initial[k]),
      };
    }
    const worstMean = (d) => Math.max(...d.mean.map(Math.abs)), worstEnd = (d) => Math.max(...d.end.map(Math.abs));
    report.push(`${firn ? 'firn' : mesh.latCell[cell] > 0 ? 'north' : 'south'} cell: layer means over the year depart from the wave's by ${worstMean(drift.wave).toFixed(3)} K (${worstMean(drift.skin).toFixed(3)} from the skin start), year end minus start ${worstEnd(drift.wave).toFixed(3)} K (${worstEnd(drift.skin).toFixed(3)}), by layer ${drift.wave.end.map((d) => d.toFixed(2)).join(' ')}`);
    checks.push([worstMean(drift.wave), worstMean(drift.skin), worstEnd(drift.wave), worstEnd(drift.skin), report.at(-1), firn]);
  }
  console.log(report.join('\n'));
  for (const [mean, skinMean, end, skinEnd, line, firn] of checks) {
    assert.ok(mean < (firn ? 0.7 : 0.3) && end < (firn ? 2.6 : 0.8), line);
    assert.ok((firn ? 1.8 : 3) * mean < skinMean && (firn ? 1.8 : 3) * end < skinEnd, line);
  }
});

const hashOf = (arrays) => { const hash = createHash('sha256'); for (const a of arrays) hash.update(Buffer.from(Float64Array.from(a).buffer)); return hash.digest('hex'); };
function slabRun() {
  const land = createLandSurface(mesh, geography, { iceSheetBalance: true });
  land.initialize();
  const surfaceT = Float64Array.from({ length: C }, (_, i) => 255 + 40 * Math.cos(mesh.latCell[i]));
  const cells = [landCell, northCell, southCell, firnCell], flux = new Float64Array(C);
  for (const i of cells) land.snow[i] = i % 2 ? 12 : 0;
  let seed = 12345;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  for (let n = 0; n < 200; n++) {
    for (const i of cells) {
      flux[i] = 500 * (random() - 0.45);
      land.deposit(i, random() < 0.2 ? 4 * random() : 0, surfaceT[i] - 1, DT);
      land.update(i, surfaceT, flux, random() < 0.3 ? 2e-5 * random() : 0, DT, surfaceT[i] - 3, 1e-5 * random());
    }
  }
  const saved = land.serialize();
  return { surfaceT, saved, hash: hashOf([surfaceT, ...Object.keys(saved).sort().map((key) => saved[key])]) };
}

test('the default slab leaves land.update bit-identical and the serialized land without a new key', () => {
  const { saved, hash } = slabRun();
  assert.ok(!Object.keys(saved).some((key) => key.startsWith('soilT')));
  assert.deepEqual(Object.keys(saved).sort(), ['canopy', 'decayMean', 'demandMean', 'litterMean', 'rainMean', 'record', 'seasonLength', 'seasonWarmth', 'snow', 'snowAlbedo', 'snowFreeCover', 'soil', 'soilCarbon', 'surface', 'vegetation']);
  assert.equal(hash, '35dd36e1e87b759003e37cf0e3a0e3cacd906ced069f8e643d29eb5c64a56e57');
  assert.equal(recordedOptions().land.soil, 'slab');
  assert.equal(recordedOptions({ land: { soil: 'layers' } }).land.soil, 'layers');
});

test('a layered land state round-trips through the state file, and regridding carries the layers between resolutions', async () => {
  const land = build();
  land.initialize(50 * 86400);
  const saved = land.serialize();
  const header = { N, day: 50, pi: Float64Array.from({ length: 4 }, (_, i) => 1e5 + i), land: saved };
  const decoded = await decodeState(encodeState(header, { f64: [] }));
  for (const name of land.deepNames) {
    assert.ok(decoded.land[name] && decoded.land[name].length === C, `${name} saved`);
    for (let i = 0; i < C; i++) assert.ok(Math.abs(decoded.land[name][i] - saved[name][i]) < 3e-5 * Math.max(1, saved[name][i] / 100), `${name} at ${i}`);
  }
  const loaded = build();
  loaded.load(decoded.land, null, 0);
  for (let i = 0; i < C; i++) if (geography.land[i]) assert.ok(Math.abs(loaded.soilT3[i] - land.soilT3[i]) < 1e-3);
  const fine = createModel(new Grid(32), { physics: false, topography }), coarse = createModel(new Grid(16), { physics: false, topography });
  const fineLand = createLandSurface(fine.mesh, fine.geography, { soil: 'layers' });
  fineLand.initialize(120 * 86400);
  const sameCells = regridLand(coarse, coarse, land.serialize());
  assert.deepEqual(Array.from(sameCells.soilT4), Array.from(land.soilT4));
  const down = regridLand(fine, coarse, fineLand.serialize(), null, { surfaceT: new Float64Array(fine.mesh.nCells).fill(280) });
  assert.ok(down.soilT2 && down.soilT3 && down.soilT4 && down.soilT2.length === C);
  const target = createLandSurface(coarse.mesh, coarse.geography, { soil: 'layers' });
  target.load(down, null, 0);
  let compared = 0;
  for (let i = 0; i < C; i++) {
    if (!coarse.geography.land[i]) { assert.equal(down.soilT3[i], 0); continue; }
    assert.ok(down.soilT3[i] > 150 && down.soilT3[i] < 330, `cell ${i}: ${down.soilT3[i]}`);
    compared++;
  }
  assert.ok(compared > C / 4);
});
