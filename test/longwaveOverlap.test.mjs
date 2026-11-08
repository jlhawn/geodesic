import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { STEFAN_BOLTZMANN, TRACE_CHAIN_FLOOR } from '../js/physics/radiation.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { GRAVITY } from '../js/dynamics/sigmaCore.module.js';
import { encodeState } from '../js/stateFile.module.js';
import { PREVIOUS_RADIATION, PREVIOUS_MOIST } from './helpers/previousDefaults.mjs';

function partlyCloudy() {
  const model = createModel(new Grid(16), { ocean: false, divergenceDamping: 0, radiation: { mixedLayerDeck: false } });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  model.step(900);
  const [pi, theta, , , q, qc] = model.state, { K } = model.core, C = model.mesh.nCells;
  model.core.diagnose(pi, theta, q, qc);
  const { exnerLayer, sigmaMid } = model.core.diagnostics;
  for (let i = 0; i < C; i++) {
    for (let k = 0; k < K; k++) {
      if (!(sigmaMid[k] > 0.3 && sigmaMid[k] < 0.9) || (k + i) % 5 === 0) continue;
      const x = k * C + i, qs = saturationHumidity(theta[x] * exnerLayer[x], pi[i] * sigmaMid[k]);
      q[x] = (0.85 + 0.01 * ((3 * i + k) % 11)) * qs;
      qc[x] = (0.002 + 0.05 * (((i + 2 * k) % 9) / 8)) * qs;
    }
  }
  model.core.diagnose(pi, theta, q, qc);
  return model;
}

function columns(model, radiation) {
  const run = createModel(model.mesh, { ocean: false, radiation: { mixedLayerDeck: false, clearSkyPass: true, ...radiation } });
  run.state.forEach((a, n) => a.set(model.state[n]));
  const [pi, theta, , surfaceT, q, qc] = run.state, { K } = run.core, C = run.mesh.nCells;
  run.core.diagnose(pi, theta, q, qc);
  run.radiation.setTime(model.time);
  const rad = run.radiation, out = { outgoing: new Float64Array(C), back: new Float64Array(C), clear: new Float64Array(C), longwave: new Float64Array(K * C), cover: new Float64Array(K * C), closure: 0, scale: 0 };
  for (let i = 0; i < C; i++) {
    rad.column(i, pi[i], theta, surfaceT[i], 5, undefined, undefined, q[(K - 1) * C + i], q, qc, 0.07, 0.07);
    out.outgoing[i] = rad.budget.outgoingLongwave; out.back[i] = rad.budget.downwardLongwave; out.clear[i] = rad.budget.clearOutgoingLongwave;
    let sum = 0;
    for (let k = 0; k < K; k++) { out.longwave[k * C + i] = rad.longwave[k * C + i]; out.cover[k * C + i] = rad.layerCover[k]; sum += rad.longwave[k * C + i]; }
    const difference = STEFAN_BOLTZMANN * surfaceT[i] ** 4 - out.back[i] - out.outgoing[i];
    out.closure = Math.max(out.closure, Math.abs(sum - difference));
    out.scale = Math.max(out.scale, out.outgoing[i]);
  }
  return out;
}

test('under exponential-random overlap in the longwave the layers\' longwave heating of every column sums to the surface\'s emission less the back radiation and the OLR, as under random overlap; the overlap raises the OLR and lowers the back radiation of partly cloudy columns, and leaves the clear-sky OLR alone', () => {
  const model = partlyCloudy();
  const exponential = columns(model, {}), random = columns(model, { longwaveOverlap: 'random' });
  const C = model.mesh.nCells;
  let raised = 0, lowered = 0, moved = 0;
  for (let i = 0; i < C; i++) {
    raised += exponential.outgoing[i] - random.outgoing[i]; lowered += random.back[i] - exponential.back[i];
    moved = Math.max(moved, Math.abs(exponential.outgoing[i] - random.outgoing[i]));
    assert.equal(exponential.clear[i], random.clear[i], `the clear-sky OLR of cell ${i}`);
  }
  console.log(`${C} partly cloudy columns at N=16: the layers' longwave closes on the flux difference to ${exponential.closure.toExponential(1)} W/m2 (random ${random.closure.toExponential(1)}); exponential-random less random OLR ${(raised / C).toFixed(3)}, back radiation ${(-lowered / C).toFixed(3)} W/m2 in the mean, at most ${moved.toFixed(2)} in a column`);
  assert.ok(exponential.closure < 1e-11 * exponential.scale && random.closure < 1e-11 * random.scale, `closure ${exponential.closure} and ${random.closure} W/m2`);
  assert.ok(raised > 0 && lowered > 0 && moved > 0.5, `the overlap raises the OLR by ${raised / C} and lowers the back radiation by ${lowered / C} W/m2, at most ${moved}`);
});

test('with a decorrelation length so short that adjacent cloudy layers overlap at random (α = 0), the two-region longwave gives the random overlap\'s OLR, back radiation and layer heating to rounding', () => {
  const model = partlyCloudy(), short = { decorrelationLength: 1e-6, decorrelationSlope: 0 };
  const chain = columns(model, short), random = columns(model, { ...short, longwaveOverlap: 'random' });
  const C = model.mesh.nCells, { K } = model.core;
  let flux = 0, heating = 0, scale = 0;
  for (let i = 0; i < C; i++) {
    flux = Math.max(flux, Math.abs(chain.outgoing[i] - random.outgoing[i]) / random.outgoing[i], Math.abs(chain.back[i] - random.back[i]) / random.back[i]);
    for (let k = 0; k < K; k++) { heating = Math.max(heating, Math.abs(chain.longwave[k * C + i] - random.longwave[k * C + i])); scale = Math.max(scale, Math.abs(random.longwave[k * C + i])); }
  }
  assert.ok(flux < 1e-13 && heating < 1e-11 * scale, `OLR and back radiation apart by ${flux} relative, layer longwave by ${heating} W/m2 of ${scale}`);
});

test('longwave overlap rejects an unknown overlap', () => {
  const model = createModel(new Grid(16), { ocean: false });
  assert.throws(() => createModel(model.mesh, { ocean: false, radiation: { longwaveOverlap: 'maximum' } }), /longwaveOverlap must be 'exponentialRandom' or 'random'/);
});

test('scripts/longwaveOverlap.mjs repeats the radiation\'s OLR and surface downward longwave under its own overlap, exponential-random by default and random on request, and finds the exponential-random OLR above the random', () => {
  const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
  const model = createModel(new Grid(16), { topography, ocean: false, radiation: PREVIOUS_RADIATION, moist: PREVIOUS_MOIST });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  model.land.initialize();
  for (let n = 0; n < 6; n++) model.step(1350);
  const [pi, theta, u, surfaceT, q, qc, ice] = model.state, { mlmSubsidence, mlmHeight, mlmGate } = model.radiation;
  const dir = mkdtempSync(join(tmpdir(), 'longwaveOverlap-')), file = join(dir, 'overlap16_day0000.bin');
  writeFileSync(file, encodeState({ N: 16, K: model.core.K, day: 0, time: model.time, terrain: true, levels: model.core.levels, pi, theta, u, surfaceT, q, qc, ice, concentration: model.seaIce.concentration, mlmSubsidence, mlmHeight, mlmGate, land: model.land.serialize() }));
  const script = (radiation) => spawnSync(process.execPath, [new URL('../scripts/longwaveOverlap.mjs', import.meta.url).pathname, file], { env: { ...process.env, RADIATION: JSON.stringify({ ...PREVIOUS_RADIATION, ...radiation }), MOIST: JSON.stringify(PREVIOUS_MOIST) }, encoding: 'utf8' });
  const runs = [script({}), script({ longwaveOverlap: 'random' })], ramped = script({ traceCloudPath: 1e-3 });
  rmSync(dir, { recursive: true, force: true });
  for (const [n, run] of runs.entries()) {
    assert.equal(run.status, 0, run.stderr);
    const [, olr, down] = run.stdout.match(/repeats the radiation's OLR to ([0-9.e+-]+) and its surface downward longwave to ([0-9.e+-]+) relative/);
    assert.ok(Number(olr) < 1e-6 && Number(down) < 1e-6, `${n ? 'random' : 'exponential-random'}: OLR to ${olr}, downward to ${down}`);
    assert.match(run.stdout, n ? /the random recomputation \(the radiation's longwaveOverlap\)/ : /the exponential-random recomputation \(the radiation's longwaveOverlap\)/);
    const [, raised] = run.stdout.match(/^global: .*exponential-random less random: OLR (-?[0-9.]+),/m);
    assert.ok(Number(raised) > 0, `exponential-random less random OLR ${raised} W/m2`);
  }
  assert.equal(ramped.status, 0, ramped.stderr);
  const [, rampedOlr, rampedDown] = ramped.stdout.match(/repeats the radiation's OLR to ([0-9.e+-]+) and its surface downward longwave to ([0-9.e+-]+) relative/);
  console.log(`traceCloudPath 1e-3 kg/m²: the script repeats the radiation's OLR to ${rampedOlr} and its surface downward longwave to ${rampedDown} relative`);
  assert.ok(Number(rampedOlr) < 1e-6 && Number(rampedDown) < 1e-6, `traceCloudPath 1e-3: OLR to ${rampedOlr}, downward to ${rampedDown}`);
});

test('with traceCloudPath a layer whose cloud water is a trace enters the longwave overlap chain with cover W/t: the OLR no longer jumps when the trace appears, and columns with no layer under the trace path are untouched', () => {
  const model = partlyCloudy(), [pi, theta, , , q, qc] = model.state, { K } = model.core, C = model.mesh.nCells;
  const { exnerLayer, sigmaMid } = model.core.diagnostics;
  const gaps = [];
  for (let i = 0; i < C; i++) for (let k = 0; k < K; k++) {
    const x = k * C + i, inside = sigmaMid[k] > 0.3 && sigmaMid[k] < 0.9;
    if (inside && (k + i) % 5 === 0) { q[x] = 0.9 * saturationHumidity(theta[x] * exnerLayer[x], pi[i] * sigmaMid[k]); qc[x] = 0; gaps.push(x); }
    else if (!inside) qc[x] = 0;
  }
  const trace = { traceCloudPath: 1e-6 }, zero = { base: columns(model, {}), trace: columns(model, trace) };
  for (const x of gaps) qc[x] = 1e-15;
  const traced = { base: columns(model, {}), trace: columns(model, trace) };
  let cliff = 0, kept = 0, untouched = 0;
  for (let i = 0; i < C; i++) {
    cliff = Math.max(cliff, Math.abs(traced.base.outgoing[i] - zero.base.outgoing[i]));
    kept = Math.max(kept, Math.abs(traced.trace.outgoing[i] - zero.trace.outgoing[i]));
    if (zero.trace.outgoing[i] === zero.base.outgoing[i] && zero.trace.back[i] === zero.base.back[i]) untouched++;
  }
  console.log(`${gaps.length} gap layers given 1e-15 kg/kg of cloud water at N=16: the OLR moves by up to ${cliff.toFixed(3)} W/m² without traceCloudPath and ${kept.toExponential(1)} with 1e-6 kg/m²; ${untouched} of ${C} columns without a trace untouched`);
  assert.ok(untouched === C, `${C - untouched} columns with no trace changed`);
  assert.ok(cliff > 0.01 && kept < 1e-3 * cliff, `the trace moves the OLR by ${cliff} W/m² without the option, ${kept} with it`);
  assert.ok(traced.trace.closure < 1e-11 * traced.trace.scale, `closure ${traced.trace.closure} W/m²`);
  for (const x of gaps) qc[x] = 1e-19;
  const deeper = columns(model, trace);
  let floored = 0, flooredDeeper = 0;
  for (let i = 0; i < C; i++) { floored += Math.abs(traced.trace.outgoing[i] - zero.trace.outgoing[i]); flooredDeeper += Math.abs(deeper.outgoing[i] - zero.trace.outgoing[i]); }
  assert.ok(floored > 0 && Math.abs(flooredDeeper / floored - 1) < 1e-2, `below TRACE_CHAIN_FLOOR the chain cover stops at 1e-6, so the trace's OLR effect no longer shrinks with W: summed over the columns ${floored} W/m² at 1e-15 kg/kg and ${flooredDeeper} at 1e-19`);
  assert.throws(() => columns(model, { traceCloudPath: 2e-3 }), /traceCloudPath must be null or an in-cloud water path/);
});

test('with traceCloudPath a layer whose in-cloud path is half the trace path enters the chain at half its cover with the same mean emissivity: alone in its column it leaves the OLR and back radiation as they were to rounding, and between two cloudy layers it moves the OLR through the overlap it no longer shares, alike at two trace paths', () => {
  const model = partlyCloudy(), [pi, theta, , , q, qc] = model.state, { K, levels } = model.core, C = model.mesh.nCells;
  const { exnerLayer, sigmaMid } = model.core.diagnostics;
  let k0 = 1;
  for (let k = 2; k < K - 1; k++) if (Math.abs(sigmaMid[k] - 0.6) < Math.abs(sigmaMid[k0] - 0.6)) k0 = k;
  const cells = Array.from({ length: C }, (_, i) => i).filter((i) => (k0 - 1 + i) % 5 && (k0 + 1 + i) % 5);
  for (const i of cells) { const x = k0 * C + i; q[x] = 0.9 * saturationHumidity(theta[x] * exnerLayer[x], pi[i] * sigmaMid[k0]); }
  const clouds = Float64Array.from(qc), path = (i) => pi[i] * (levels[k0 + 1] - levels[k0]) / GRAVITY;
  const share = (cover, i, t) => qc[k0 * C + i] * path(i) / (cover[k0 * C + i] * t);
  const halve = (t) => {
    const prior = new Map();
    qc.set(clouds);
    for (const i of cells) qc[k0 * C + i] = 1e-6;
    for (let n = 0; n < 8; n++) {
      const { cover } = columns(model, {});
      for (const i of cells) {
        const x = k0 * C + i, r = Math.log(share(cover, i, t) / 0.5), lnq = Math.log(qc[x]), last = prior.get(i);
        const slope = last && lnq !== last[0] ? (r - last[1]) / (lnq - last[0]) : 1;
        prior.set(i, [lnq, r]);
        qc[x] = Math.exp(lnq - r / (slope > 0.1 ? slope : 1));
      }
    }
    return cells.map((i) => qc[k0 * C + i]);
  };
  const moves = (t, traced, alone) => {
    qc.set(clouds);
    cells.forEach((i, n) => { if (alone) for (let k = 0; k < K; k++) qc[k * C + i] = 0; qc[k0 * C + i] = traced[n]; });
    const off = columns(model, {}), on = columns(model, { traceCloudPath: t });
    return { half: Math.max(...cells.map((i) => Math.abs(share(off.cover, i, t) - 0.5))), closure: on.closure / on.scale, olr: Math.max(...cells.map((i) => Math.abs(on.outgoing[i] - off.outgoing[i]))), back: Math.max(...cells.map((i) => Math.abs(on.back[i] - off.back[i]))) };
  };
  const results = [1e-3, 1e-4].map((t) => { const traced = halve(t); return { t, lone: moves(t, traced, true), sandwiched: moves(t, traced, false) }; });
  for (const { t, lone, sandwiched } of results) console.log(`t = ${t} kg/m2, W/f = t/2 to ${Math.max(lone.half, sandwiched.half).toExponential(1)} in ${cells.length} columns: the OLR moves by up to ${lone.olr.toExponential(2)} W/m2 and its back radiation by ${lone.back.toExponential(2)} for a lone layer, the OLR by ${sandwiched.olr.toExponential(2)} for one between two cloudy layers; closure ${lone.closure.toExponential(1)}, ${sandwiched.closure.toExponential(1)}`);
  for (const { t, lone, sandwiched } of results) assert.ok(lone.half < 1e-9 && sandwiched.half < 1e-9 && lone.closure < 1e-11 && sandwiched.closure < 1e-11, `t = ${t}: W/f is t/2 to ${lone.half}, ${sandwiched.half}; closure ${lone.closure}, ${sandwiched.closure}`);
  for (const { t, lone } of results) assert.ok(lone.olr < 1e-10 && lone.back < 1e-10, `t = ${t}: a lone layer at half its cover and twice its in-cloud emissivity moves the OLR by ${lone.olr} and the back radiation by ${lone.back} W/m2`);
  const [wide, narrow] = results.map((r) => r.sandwiched.olr);
  assert.ok(wide > 5e-3 && Math.abs(narrow / wide - 1) < 0.05, `between two cloudy layers the halved chain cover moves the OLR by ${wide} W/m2 at t = 1e-3 and ${narrow} at 1e-4: it is the cover the trace no longer shares, not its emissivity, that moves it`);
});

test('the GPU kernels carry the trace-chain floor as one literal that rounds to TRACE_CHAIN_FLOOR in f32', async () => {
  const { physicsConstants } = await import('../js/gpu/physics.gpu.js');
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  const literals = [...physicsConstants({ ...PHYSICS_DEFAULTS, R: 287 }).matchAll(/const TRACE_FLOOR: f32 = ([^;]+);/g)].map((m) => m[1]);
  assert.ok(literals.length === 1 && Math.fround(Number(literals[0])) === Math.fround(TRACE_CHAIN_FLOOR), `TRACE_FLOOR literals ${JSON.stringify(literals)} against TRACE_CHAIN_FLOOR ${TRACE_CHAIN_FLOOR}, ${Math.fround(TRACE_CHAIN_FLOOR)} in f32`);
});
test('the GPU kernels apply the trace-chain ramp and its floor with the CPU expression, max(W/(f t), min(1, TRACE_CHAIN_FLOOR / chain cover)), under the same guard', () => {
  const gpu = readFileSync(new URL('../js/gpu/physics.gpu.js', import.meta.url), 'utf8'), cpu = readFileSync(new URL('../js/physics/radiation.module.js', import.meta.url), 'utf8');
  const wgsl = [...gpu.matchAll(/if \((LW_CHAIN && TRACE_CHAIN && [^{]+)\) \{ chainF\[k\] \*= ([^;]+); \}/g)], js = [...cpu.matchAll(/if \((traceCloudPath !== null && [^{]+)\) \{ const kept = ([^;]+); chainCloud\[k\] \/= kept; chainCover\[k\] \*= kept; \}/g)];
  const names = [[/Math\.(max|min)/g, '$1'], [/layerCover\[k\]/g, 'f'], [/traceCloudPath !== null/g, 'TRACE_CHAIN'], [/traceCloudPath/g, 'TRACE_PATH'], [/TRACE_CHAIN_FLOOR/g, 'TRACE_FLOOR'], [/chainCover\[k\]/g, 'chainF[k]'], [/\b([01])\.0\b/g, '$1'], [/^LW_CHAIN && /, '']];
  const normal = (text) => names.reduce((t, [from, to]) => t.replace(from, to), text).replace(/\s+/g, '');
  assert.ok(wgsl.length === 1 && js.length === 1, `${wgsl.length} WGSL and ${js.length} CPU trace-chain rules`);
  assert.equal(normal(wgsl[0][2]), normal(js[0][2]), 'the ramp and its floor');
  assert.equal(normal(wgsl[0][1]), normal(js[0][1]), 'the guard');
});
