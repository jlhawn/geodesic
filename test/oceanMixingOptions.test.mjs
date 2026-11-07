import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean, backgroundDiffusivities, HENYEY_RATIO, EQUATORIAL_BACKGROUND, restratificationScale, restratificationOptions, deformationRadii, limitedFaceValue, limitedEdgeValue, eddyDiffusivities, shearResolution, shearResolutionOptions, isoneutralOptions, ISONEUTRAL } from '../js/ocean/layered.module.js';
import { interiorMixing, lmd94Shape, LMD94_VISCOSITY } from '../js/ocean/column.module.js';

/*
 * The ocean's optional interior closures on the CPU at N = 16, each against
 * its published form and against the default it leaves alone: the
 * latitude dependence of the background diffusivity (backgroundLatitude)
 * and Large, McWilliams and Doney's interior shear mixing
 * (interiorShear), the restratification's equatorial taper
 * (mixedLayerRestratification.taperLatitude), the temperature and
 * salinity the eddy phase's layer fluxes carry (eddyTracers) and the
 * resolution function on the shear-scaled diffusion
 * (shearDiffusivityResolution), the factor on the ∇⁴ closure
 * (closureScale), and the isoneutral diffusion of the fixed-depth layers
 * (isoneutralDiffusion).
 */
const mesh = buildMesh(new Grid(16)), { nCells: C, nEdges: E, latCell, latEdge, fEdge, dcEdge, areaCell } = mesh;
const DT = 1350, DEGREE = Math.PI / 180;
const latitudes = (degrees) => ({ nCells: degrees.length, latCell: Float64Array.from(degrees, (d) => d * DEGREE) });

test('backgroundLatitude henyey is MOM6\'s Henyey scaling: κ_b at 30°, the floor at the equator, 1.68 κ_b at the poles, symmetric and rising with |latitude|', () => {
  const degrees = [-90, -60, -30, -15, -5, -1, 0, 1, 5, 15, 30, 60, 90], n = degrees.length;
  const kappa = backgroundDiffusivities(latitudes(degrees), 'henyey', 1e-5, 1e-6);
  const G = (d) => { const s = Math.max(1e-10, Math.abs(Math.sin(d * DEGREE))), x = HENYEY_RATIO / s; return s * Math.log(x + Math.sqrt(x * x - 1)); };
  for (let m = 0; m < n; m++) {
    const expected = Math.max(1e-6, 1e-5 * G(degrees[m]) / G(30));
    assert.ok(Math.abs(kappa[m] - expected) <= 1e-12 * expected, `${degrees[m]}°: ${kappa[m]} against ${expected}`);
    assert.equal(kappa[m], kappa[n - 1 - m]);
  }
  assert.ok(Math.abs(kappa[10] - 1e-5) < 1e-17, `30°: ${kappa[10]}`);
  assert.equal(kappa[6], 1e-6);
  assert.ok(Math.abs(kappa[12] / 1e-5 - 1.684) < 1e-3, `pole: ${kappa[12]}`);
  for (let m = 7; m < n; m++) assert.ok(kappa[m] >= kappa[m - 1], `${degrees[m]}°`);
  const om4 = backgroundDiffusivities(latitudes([0, 2, 5]), 'henyey', 1.5e-5, 2e-6);
  assert.deepEqual(Array.from(om4.subarray(0, 2)), [2e-6, 2e-6]);
  assert.ok(om4[2] > 3.6e-6 && om4[2] < 3.7e-6, `OM4's 5°: ${om4[2]}`);
});

test('backgroundLatitude equatorial is MOM025\'s: κ_eq within 5°, κ_b poleward of 15°, and a raised cosine between with the mean at 10°', () => {
  const degrees = [0, 3, 5, 7.5, 10, 12.5, 15.5, 20, 45, 90, -10, -90];
  const kappa = backgroundDiffusivities(latitudes(degrees), 'equatorial', 1e-5, 1e-6);
  assert.equal(kappa[0], 1e-6); assert.equal(kappa[1], 1e-6); assert.equal(kappa[2], 1e-6);
  for (const m of [6, 7, 8, 9, 11]) assert.equal(kappa[m], 1e-5, `${degrees[m]}°`);
  assert.ok(Math.abs(kappa[4] - 5.5e-6) < 1e-18, `10°: ${kappa[4]}`);
  assert.equal(kappa[10], kappa[4]);
  const quarter = 1e-6 + 9e-6 * 0.5 * (1 - Math.cos(Math.PI / 4));
  assert.ok(Math.abs(kappa[3] - quarter) < 1e-18, `7.5°: ${kappa[3]}`);
  for (let m = 1; m < 10; m++) assert.ok(kappa[m] >= kappa[m - 1], `${degrees[m]}°`);
  assert.deepEqual(EQUATORIAL_BACKGROUND, { inner: 5, outer: 15 });
});

test('backgroundLatitude uniform holds κ_b in every cell, and an unknown shape or a negative equatorial diffusivity is refused', () => {
  const kappa = backgroundDiffusivities(mesh, 'uniform', 1e-5, 1e-6);
  assert.equal(kappa.length, C);
  for (const v of kappa) assert.equal(v, 1e-5);
  assert.throws(() => backgroundDiffusivities(mesh, 'tanh', 1e-5, 1e-6), /backgroundLatitude is one of/);
  assert.throws(() => createOcean(mesh, { backgroundLatitude: 'step' }), /backgroundLatitude is one of/);
  assert.throws(() => createOcean(mesh, { backgroundLatitude: 'equatorial', equatorialDiffusivity: -1e-6 }), /equatorialDiffusivity/);
});

function stratified(options) {
  const ocean = createOcean(mesh, { everySteps: 1, iceDrift: false, ...options });
  ocean.initialize(Float64Array.from(latCell, (lat) => 271.5 + 30 * Math.cos(lat) ** 2), new Float64Array(C));
  return ocean;
}

test('one column phase under backgroundLatitude equatorial leaves every column poleward of 15° as the uniform background does, bit for bit, and moves a tenth of the heat within 5°', () => {
  const uniform = stratified({}), equatorial = stratified({ backgroundLatitude: 'equatorial' }), L = uniform.layers;
  assert.deepEqual(Array.from(equatorial.Q), Array.from(uniform.Q));
  const before = Float64Array.from(uniform.Q);
  uniform.phases.column(DT); equatorial.phases.column(DT);
  let near = 0, moved = 0, movedLess = 0;
  for (let i = 0; i < C; i++) {
    if (!uniform.cellOcean[i]) continue;
    const degrees = Math.abs(latCell[i]) / DEGREE;
    let a = 0, b = 0;
    for (let k = 0; k < L; k++) {
      const n = k * C + i;
      if (degrees > 15.01) { assert.equal(equatorial.Q[n], uniform.Q[n]); assert.equal(equatorial.h[n], uniform.h[n]); assert.equal(equatorial.W[n], uniform.W[n]); }
      a += Math.abs(uniform.Q[n] - before[n]); b += Math.abs(equatorial.Q[n] - before[n]);
    }
    if (degrees < 5 && a > 0) { near++; moved += a; movedLess += b; }
  }
  assert.ok(near > 0, 'no stratified column within 5°');
  const ratio = movedLess / moved;
  assert.ok(ratio > 0.08 && ratio < 0.12, `within 5° the column moves ${ratio} of the uniform background's heat`);
});

test('lmd94Shape is Large, McWilliams and Doney\'s interior shape: 1 at Ri 0 and below, (3/4)³ at Ri 0.35, 0 at and above 0.7', () => {
  assert.equal(lmd94Shape(0), 1);
  assert.equal(lmd94Shape(-2), 1);
  assert.ok(Math.abs(lmd94Shape(0.35) - 0.421875) < 1e-15);
  assert.equal(lmd94Shape(0.7), 0);
  assert.equal(lmd94Shape(5), 0);
  for (let r = 0.01; r < 0.7; r += 0.01) assert.ok(lmd94Shape(r) < lmd94Shape(r - 0.01));
});

test('under interiorShear lmd94 the interior shear term is ν₀[1 − (Ri/0.7)²]³ in viscosity and diffusivity alike, over the backgrounds, with ν₀ lmd94Viscosity, 5e-3 unless set, whatever shearViscosity is; pp81, the default, is Pacanowski and Philander\'s', () => {
  const K = 8, h = new Float64Array(K).fill(10), T = Float64Array.from({ length: K }, (_, k) => 293.15 - 0.5 * k), S = new Float64Array(K).fill(35);
  const u = Float64Array.of(1, 0.75, 0.55, 0.4, 0.3, 0.25, 0.24, 0.235), v = new Float64Array(K);
  const scratch = { rho: new Float64Array(K) }, nu0 = 5e-3, nub = 1e-4, kb = 1e-5;
  const run = (extra) => { const kappa = new Float64Array(K), nu = new Float64Array(K); interiorMixing({ h, T, S, u, v, K, shearViscosity: nu0, backgroundViscosity: nub, backgroundDiffusivity: kb, ...extra }, kappa, nu, scratch); return { kappa, nu }; };
  const off = run({}), pp81 = run({ interiorShear: 'pp81' }), lmd = run({ interiorShear: 'lmd94' }), lmdSheared = run({ interiorShear: 'lmd94', shearViscosity: 3e-2 }), own = 2e-3, lmdOwn = run({ interiorShear: 'lmd94', lmd94Viscosity: own });
  assert.equal(LMD94_VISCOSITY, 5e-3);
  assert.deepEqual(Array.from(pp81.nu), Array.from(off.nu));
  assert.deepEqual(Array.from(pp81.kappa), Array.from(off.kappa));
  assert.deepEqual(Array.from(lmdSheared.nu), Array.from(lmd.nu));
  assert.deepEqual(Array.from(lmdSheared.kappa), Array.from(lmd.kappa));
  let below = 0, above = 0;
  for (let k = 0; k < K - 1; k++) {
    const ri = (Math.sqrt(nu0 / (off.nu[k] - nub)) - 1) / 5, shape = Math.max(0, 1 - (ri / 0.7) ** 2) ** 3, shear = LMD94_VISCOSITY * shape;
    assert.ok(Math.abs(off.kappa[k] - (off.nu[k] / (1 + 5 * ri) + kb)) < 1e-9 * off.kappa[k], `PP81 at interface ${k}`);
    assert.ok(Math.abs(lmd.nu[k] - (shear + nub)) < 1e-9 * lmd.nu[k], `interface ${k} at Ri ${ri}: ${lmd.nu[k]} against ${shear + nub}`);
    assert.ok(Math.abs(lmd.kappa[k] - (shear + kb)) < 1e-9 * lmd.nu[k], `interface ${k} at Ri ${ri}: ${lmd.kappa[k]} against ${shear + kb}`);
    assert.ok(Math.abs(lmdOwn.kappa[k] - (own * shape + kb)) < 1e-9 * lmdOwn.nu[k], `interface ${k} at Ri ${ri} under lmd94Viscosity ${own}: ${lmdOwn.kappa[k]} against ${own * shape + kb}`);
    if (ri < 0.7) below++; else above++;
  }
  assert.ok(below >= 2 && above >= 2, `${below} interfaces below Ri 0.7 and ${above} above`);
  assert.throws(() => createOcean(mesh, { interiorShear: 'kpp' }), /interiorShear/);
  assert.throws(() => createOcean(mesh, { interiorShear: 'lmd94', lmd94Viscosity: -1e-3 }), /lmd94Viscosity/);
});

test('interiorShear lmd94 reaches the ocean\'s column phase: under a sheared surface current it raises some interfaces\' viscosity and lowers others\'', () => {
  const sheared = (options) => {
    const ocean = stratified(options), L = ocean.layers;
    for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) ocean.u[k * E + e] = 0.4 * Math.exp(-k / 3) * Math.cos(latEdge[e]) * Math.sin(3 * e);
    ocean.phases.column(DT);
    return ocean;
  };
  const base = sheared({}), lmd = sheared({ interiorShear: 'lmd94' });
  let more = 0, less = 0;
  for (let n = 0; n < base.viscosity.length; n++) {
    if (lmd.viscosity[n] > base.viscosity[n] * (1 + 1e-6)) more++;
    if (lmd.viscosity[n] < base.viscosity[n] * (1 - 1e-6)) less++;
  }
  assert.ok(more > 0 && less > 0, `${more} viscosities raised and ${less} lowered`);
});

test('taperLatitude scales the restratification by |f|/√(f² + f_t²): none without it, 1/√2 at the taper latitude, above 0.98 at 30° and below 0.2 at 1°', () => {
  const plain = restratificationScale(mesh, restratificationOptions(true)), radii = deformationRadii(mesh);
  for (let e = 0; e < E; e++) assert.equal(plain[e], 1 * Math.min(dcEdge[e], radii[e]) / 500);
  const tapered = restratificationScale(mesh, restratificationOptions({ taperLatitude: 5 })), ft = 2 * mesh.omega * Math.sin(5 * DEGREE);
  for (let e = 0; e < E; e++) {
    const expected = Math.abs(fEdge[e]) / Math.sqrt(fEdge[e] ** 2 + ft ** 2);
    assert.ok(Math.abs(tapered[e] / plain[e] - expected) < 1e-12, `edge ${e}: ${tapered[e] / plain[e]} against ${expected}`);
  }
  const at = (degrees) => { const f = 2 * mesh.omega * Math.sin(degrees * DEGREE); return Math.abs(f) / Math.hypot(f, ft); };
  assert.ok(Math.abs(at(5) - Math.SQRT1_2) < 1e-12 && at(30) > 0.98 && at(1) < 0.2);
  assert.throws(() => restratificationOptions({ taperLatitude: -1 }), /taperLatitude/);
  assert.throws(() => restratificationOptions({ taperLatitude: 90 }), /taperLatitude/);
});

test('limitedFaceValue is limitedEdgeValue on the cells\' values: the midpoint on a line, the upwind value at an extremum', () => {
  const field = Float64Array.of(1, 2, 3, 2.5, 0.5);
  for (const [up, down, uu] of [[1, 2, 0], [2, 1, 3], [3, 4, 2], [1, 0, 2], [3, 2, 4]]) {
    const ratio = 1.3 / 1.1;
    assert.ok(Math.abs(limitedFaceValue(field[up], field[down], field[uu], ratio) - limitedEdgeValue(field, up, down, uu, 1.1, 1.3)) < 1e-15);
  }
  assert.equal(limitedFaceValue(2, 3, 1, 1), 2.5);
  assert.equal(limitedFaceValue(3, 2, 2.5, 1), 3);
});

function uniformFront(options) {
  const ocean = createOcean(mesh, { everySteps: 1, iceDrift: false, eddyDiffusivity: 0, ...options }), L = ocean.layers;
  ocean.initialize(new Float64Array(C).fill(288), new Float64Array(C));
  for (let i = 0; i < C; i++) {
    let top = 0;
    for (let k = 0; k < L && top < 150; k++) { const n = k * C + i; ocean.Q[n] = ocean.h[n] * (288 + 4 * Math.tanh((latCell[i] - 0.8) / 0.12)); ocean.W[n] = ocean.h[n] * 35; top += ocean.h[n]; }
  }
  ocean.boundaryDepths.fill(100);
  return ocean;
}
const columnHeat = (ocean) => Float64Array.from({ length: C }, (_, i) => { let s = 0; for (let k = 0; k < ocean.layers; k++) s += ocean.Q[k * C + i]; return s; });

test('across a front uniform with depth the donor-cell eddy fluxes carry heat from the warm columns to the cold, the centered ones none, and the limited ones less than the donor\'s', () => {
  const moved = {};
  for (const eddyTracers of ['donor', 'third', 'centered']) {
    const ocean = uniformFront({ eddyTracers, mixedLayerRestratification: { coefficient: 20 } }), before = columnHeat(ocean);
    ocean.phases.eddy(DT);
    const after = columnHeat(ocean);
    let total = 0, warmed = 0, scale = 0;
    for (let i = 0; i < C; i++) { const d = (after[i] - before[i]) * areaCell[i]; total += Math.abs(d); scale += Math.abs(before[i]) * areaCell[i]; if (latCell[i] < 0.8 && d > 0) warmed += d; }
    moved[eddyTracers] = { total, warmed, scale };
  }
  console.log(`column heat moved across the front (m·K·m²): ${Object.entries(moved).map(([k, v]) => `${k} ${v.total.toExponential(3)}`).join(', ')}`);
  assert.ok(moved.donor.warmed > 0, 'the donor-cell fluxes did not warm the cold side');
  assert.ok(moved.centered.total < 1e-9 * moved.centered.scale, `the centered fluxes moved ${moved.centered.total}`);
  assert.ok(moved.third.total < 0.5 * moved.donor.total, `the limited fluxes moved ${moved.third.total} against the donor's ${moved.donor.total}`);
});

test('under eddyTracers third and centered the eddy phase keeps water, heat and salt to rounding, the limited values every layer\'s temperature within its start\'s range and the centered ones within a kelvin of it, with the coefficient far too large', () => {
  for (const eddyTracers of ['third', 'centered']) {
    const ocean = uniformFront({ eddyTracers, eddyDiffusivity: 1000, mixedLayerRestratification: { coefficient: 4000 } }), L = ocean.layers;
    const sums = () => { let v = 0, q = 0, w = 0; for (let n = 0; n < L * C; n++) { const a = areaCell[n % C]; v += a * ocean.h[n]; q += a * ocean.Q[n]; w += a * ocean.W[n]; } return { v, q, w }; };
    const before = sums();
    let low = Infinity, high = -Infinity;
    for (let n = 0; n < L * C; n++) if (ocean.h[n] > 1e-6) { const t = ocean.Q[n] / ocean.h[n]; low = Math.min(low, t); high = Math.max(high, t); }
    for (let n = 0; n < 200; n++) ocean.phases.eddy(DT);
    const after = sums();
    for (const name of ['v', 'q', 'w']) assert.ok(Math.abs(after[name] / before[name] - 1) < 1e-12, `${eddyTracers} ${name}: ${after[name] / before[name] - 1}`);
    let lowAfter = Infinity, highAfter = -Infinity;
    for (let n = 0; n < L * C; n++) if (ocean.h[n] > 0.1) { const t = ocean.Q[n] / ocean.h[n]; assert.ok(Number.isFinite(t)); lowAfter = Math.min(lowAfter, t); highAfter = Math.max(highAfter, t); }
    console.log(`${eddyTracers}: layer temperatures ${low.toFixed(3)}-${high.toFixed(3)} K before, ${lowAfter.toFixed(3)}-${highAfter.toFixed(3)} K after 200 eddy phases`);
    const slack = eddyTracers === 'third' ? 1e-9 : 1;
    assert.ok(lowAfter > low - slack && highAfter < high + slack, `${eddyTracers}: ${lowAfter}-${highAfter} K from ${low}-${high} K`);
  }
  assert.throws(() => createOcean(mesh, { eddyTracers: 'upwind' }), /eddyTracers is one of/);
  assert.throws(() => createOcean(mesh, { eddyTracers: 'third', tracerAdvection: 'upwind' }), /beyond cells/);
});

test('shearDiffusivityResolution is MOM6\'s resolution function 1/(1 + (s R_d/Δ)^p): the eddy diffusivity\'s form under true, least near the equator, whole at high latitudes, refused without a positive scale and power', () => {
  assert.equal(shearResolutionOptions(false), null);
  assert.deepEqual(shearResolutionOptions(true), { scale: 1, power: 2 });
  assert.deepEqual(shearResolutionOptions({ power: 4 }), { scale: 1, power: 4 });
  const radius = deformationRadii(mesh), F = shearResolution(mesh, shearResolutionOptions(true)), G = shearResolution(mesh, { scale: 0.5, power: 4 }), kappa = eddyDiffusivities(mesh, 1000);
  let least = 0;
  for (let e = 0; e < E; e++) {
    assert.ok(Math.abs(F[e] - kappa[e] / 1000) < 1e-15, `edge ${e}: ${F[e]} against the eddy diffusivity's ${kappa[e] / 1000}`);
    assert.ok(Math.abs(G[e] - 1 / (1 + (0.5 * radius[e] / dcEdge[e]) ** 4)) < 1e-15, `edge ${e}: ${G[e]}`);
    assert.ok(F[e] > 0 && F[e] <= 1);
    if (Math.abs(latEdge[e]) > 45 * DEGREE) assert.ok(F[e] > 0.99, `edge ${e} at ${latEdge[e] / DEGREE}°: ${F[e]}`);
    if (F[e] < F[least]) least = e;
  }
  assert.ok(Math.abs(latEdge[least]) < 5 * DEGREE && F[least] < 0.9, `least ${F[least]} at ${latEdge[least] / DEGREE}°`);
  for (const bad of [{ scale: 0 }, { power: -1 }, { scale: NaN }]) assert.throws(() => createOcean(mesh, { shearDiffusivityResolution: bad }), /shearDiffusivityResolution takes a positive scale and power/);
});

test('shearDiffusivityResolution reaches the layers\' tracer tendency: a vanishing scale keeps the shear-scaled diffusion whole, a huge one takes it all away, true some of it, and the dynamics keep the sea\'s heat to rounding', () => {
  const sum = (ocean) => { let q = 0; for (let n = 0; n < ocean.Q.length; n++) q += areaCell[n % C] * ocean.Q[n]; return q; };
  const run = (options) => {
    const ocean = stratified({ tracerAdvection: 'upwind', ...options }), L = ocean.layers;
    for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) ocean.u[k * E + e] = 0.4 * Math.exp(-k / 3) * Math.cos(latEdge[e]) * Math.sin(3 * e);
    const before = sum(ocean);
    ocean.phases.dynamics(DT);
    assert.ok(Math.abs(sum(ocean) - before) < 1e-12 * Math.abs(before), `${JSON.stringify(options)}: heat ${before} → ${sum(ocean)}`);
    return Float64Array.from(ocean.Q);
  };
  const apart = (x, y) => { let d = 0; for (let n = 0; n < x.length; n++) d += Math.abs(x[n] - y[n]); return d; };
  const none = run({ shearDiffusivity: 0 }), whole = run({}), kept = run({ shearDiffusivityResolution: { scale: 1e-6 } }), removed = run({ shearDiffusivityResolution: { scale: 1e6 } }), tapered = run({ shearDiffusivityResolution: true });
  const span = apart(whole, none);
  assert.ok(span > 0, 'the shear-scaled diffusion moved no heat');
  assert.ok(apart(kept, whole) < 1e-7 * span, `a vanishing scale: ${apart(kept, whole) / span}`);
  assert.ok(apart(removed, none) < 1e-7 * span, `a huge scale: ${apart(removed, none) / span}`);
  const fromWhole = apart(tapered, whole) / span, fromNone = apart(tapered, none) / span;
  console.log(`true moves the heat ${fromWhole.toFixed(4)} of the shear-scaled diffusion's change from the whole and ${fromNone.toFixed(4)} from none`);
  assert.ok(fromWhole > 1e-3 && fromNone > 0.5, `${fromWhole}, ${fromNone}`);
});

test('closureScale scales the ∇⁴ closure: 1 is the default and 0 no closure, bit for bit, and a half takes about half the closure\'s change to a sheared current', () => {
  const run = (options) => {
    const ocean = stratified(options), L = ocean.layers;
    for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) ocean.u[k * E + e] = 0.4 * Math.exp(-k / 3) * Math.cos(latEdge[e]) * Math.sin(3 * e);
    ocean.phases.dynamics(DT);
    return Float64Array.from(ocean.u);
  };
  const apart = (x, y) => { let d = 0; for (let n = 0; n < x.length; n++) d += Math.abs(x[n] - y[n]); return d; };
  const whole = run({}), none = run({ closureHours: 0 });
  assert.deepEqual(run({ closureScale: 1 }), whole);
  assert.deepEqual(run({ closureScale: 0 }), none);
  const half = apart(run({ closureScale: 0.5 }), none) / apart(whole, none);
  console.log(`closureScale 0.5 makes ${half.toFixed(4)} of the closure's change`);
  assert.ok(half > 0.45 && half < 0.55, `${half}`);
  for (const bad of [-1, NaN, Infinity]) assert.throws(() => createOcean(mesh, { closureScale: bad }), /closureScale is a non-negative factor/);
});

const FULL_STEP = ['dynamics', 'rescale', 'eddy', 'shift', 'regrid', 'column'];
const profiled = (options, profile) => {
  const ocean = stratified(options), L = ocean.layers, { h, Q, W } = ocean;
  ocean.u.fill(0);
  for (let i = 0; i < C; i++) {
    let z = 0;
    for (let k = 0; k < L; k++) { const n = k * C + i, hv = h[n], [t, s] = profile(z + 0.5 * hv, latCell[i], i); z += hv; Q[n] = hv * t; W[n] = hv * s; }
  }
  return ocean;
};
const tracerTendency = (ocean) => {
  ocean.tendency(ocean.state, ocean.stages[0]);
  return { dQ: Float64Array.from(ocean.stages[0][2]), dW: Float64Array.from(ocean.stages[0][3]), k33: ocean.isoKappa ? Float64Array.from(ocean.isoKappa) : null };
};
const tracerValues = (ocean, store) => Float64Array.from(store, (q, n) => {
  if (ocean.cellOcean[n % C]) assert.ok(ocean.h[n] > 1e-6, `layer ${Math.floor(n / C)} of cell ${n % C} has vanished`);
  return ocean.h[n] > 1e-6 ? q / ocean.h[n] : 0;
});
const explicitKappa33 = (ocean, k33, values) => {
  const { h } = ocean, F = ocean.table.fixed, out = new Float64Array(ocean.layers * C);
  for (let i = 0; i < C; i++) for (let m = 1; m <= F; m++) {
    const up = (m - 1) * C + i, down = m * C + i, flux = -k33[up] * (values[up] - values[down]) / Math.max(0.5 * (h[up] + h[down]), 1e-3);
    out[up] += flux; out[down] -= flux;
  }
  return out;
};

test('isoneutralDiffusion takes true for a slope cap of 0.01 under the boundary taper, refuses a negative or infinite cap and a non-boolean taper, and under slopeCap 0 is the horizontal diffusion, a step of the ocean bit for bit', () => {
  assert.equal(isoneutralOptions(false), null);
  assert.deepEqual(isoneutralOptions(true), { slopeCap: 0.01, boundaryTaper: true });
  assert.deepEqual(ISONEUTRAL, { slopeCap: 0.01, boundaryTaper: true });
  assert.deepEqual(isoneutralOptions({ slopeCap: 0.02 }), { slopeCap: 0.02, boundaryTaper: true });
  for (const bad of [{ slopeCap: -1 }, { slopeCap: Infinity }, { slopeCap: NaN }, { boundaryTaper: 1 }]) assert.throws(() => createOcean(mesh, { isoneutralDiffusion: bad }), /isoneutralDiffusion takes a non-negative slopeCap/);
  assert.equal(createOcean(mesh, {}).isoKappa, null);
  const run = (options) => {
    const ocean = stratified(options), L = ocean.layers;
    for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) ocean.u[k * E + e] = 0.4 * Math.exp(-k / 3) * Math.cos(latEdge[e]) * Math.sin(3 * e);
    for (let n = 0; n < 2; n++) for (const phase of FULL_STEP) ocean.phases[phase](DT);
    return ocean;
  };
  const off = run({}), flat = run({ isoneutralDiffusion: { slopeCap: 0 } }), on = run({ isoneutralDiffusion: true });
  for (const field of ['h', 'u', 'Q', 'W', 'eta']) {
    let differ = 0;
    for (let n = 0; n < off[field].length; n++) if (off[field][n] !== flat[field][n]) differ++;
    assert.equal(differ, 0, `${field}: ${differ} values differ`);
  }
  let moved = 0;
  for (let n = 0; n < off.Q.length; n++) moved += Math.abs(on.Q[n] - off.Q[n]);
  assert.ok(moved > 0, 'the rotation moved no heat');
});

test('under isoneutralDiffusion the triads move no heat or salt through the interior of a stratification whose temperature and salinity are linear in one tilted coordinate, which the horizontal diffusion crosses', () => {
  const tilt = 1e-5, radius = 6.371e6, profile = (depth, lat) => { const x = depth + tilt * radius * lat; return [288 - 0.003 * x, 34.5 + 2e-4 * x]; };
  const off = profiled({}, profile), on = profiled({ isoneutralDiffusion: { boundaryTaper: false } }, profile), F = on.table.fixed;
  const flat = tracerTendency(off), turned = tracerTendency(on), { cellsOnCell, nEdgesOnCell, maxEdges } = mesh;
  const thick = (k, i) => k >= 0 && on.h[k * C + i] > 1;
  const interior = (k, i) => {
    for (let m = -1; m < nEdgesOnCell[i]; m++) { const j = m < 0 ? i : cellsOnCell[maxEdges * i + m]; if (!(thick(k - 1, j) && thick(k, j) && thick(k + 1, j))) return false; }
    return true;
  };
  let inside = 0;
  for (let k = 1; k < F; k++) for (let i = 0; i < C; i++) if (on.cellOcean[i] && interior(k, i)) inside++;
  assert.ok(inside > 0.5 * (F - 1) * C * 0.6, `${inside} interior cells`);
  for (const [store, change, label] of [['Q', 'dQ', 'T'], ['W', 'dW', 'S']]) {
    const values = tracerValues(on, on[store]);
    const kappa33 = explicitKappa33(on, turned.k33, values);
    let scale = 0, residual = 0;
    for (let k = 1; k <= F; k++) for (let i = 0; i < C; i++) {
      const n = k * C + i, rest = turned[change][n] + kappa33[n] - (k === F ? flat[change][n] : 0);
      if (k < F) scale = Math.max(scale, Math.abs(flat[change][n]));
      if (k === F || interior(k, i)) residual = Math.max(residual, Math.abs(rest));
    }
    for (let n = (F + 1) * C; n < on.layers * C; n++) assert.equal(turned[change][n], flat[change][n]);
    console.log(`${label}: the horizontal diffusion's largest tendency ${scale.toExponential(3)}, the triads' with κS² ${residual.toExponential(3)}`);
    assert.ok(scale > 0 && residual < 1e-9 * scale, `${label}: ${residual} against ${scale}`);
  }
});

test('the triads with their κS² diffusivity applied explicitly never raise the variance of temperature or salinity, and under isoneutralDiffusion the dynamics, eddy, regrid and column phases keep the sea\'s heat and salt to rounding', () => {
  const { lonCell } = mesh, profile = (depth, lat, i) => [288 - 0.01 * depth + 3 * Math.cos(lat) ** 2 + 0.5 * Math.sin(7 * lonCell[i]) * Math.exp(-depth / 150), 34.6 + 2e-4 * depth + 0.3 * Math.sin(lat)];
  const off = profiled({}, profile), on = profiled({ isoneutralDiffusion: true }, profile);
  on.boundaryDepths.fill(60);
  const flat = tracerTendency(off), turned = tracerTendency(on);
  for (const [store, change, mean] of [['Q', 'dQ', 288], ['W', 'dW', 35]]) {
    const values = tracerValues(on, on[store]), kappa33 = explicitKappa33(on, turned.k33, values);
    let rotated = 0, horizontal = 0, scale = 0;
    for (let n = 0; n < values.length; n++) {
      const a = areaCell[n % C] * (values[n] - mean);
      rotated += a * (turned[change][n] + kappa33[n]); horizontal += a * flat[change][n]; scale += Math.abs(a * flat[change][n]);
    }
    console.log(`${store}: variance tendency ${(rotated / scale).toExponential(3)} rotated and ${(horizontal / scale).toExponential(3)} horizontal, of the horizontal's magnitude`);
    assert.ok(rotated <= 1e-12 * scale && horizontal <= 1e-12 * scale, `${store}: ${rotated}, ${horizontal} against ${scale}`);
  }
  const sheared = stratified({ isoneutralDiffusion: true }), L = sheared.layers;
  for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) sheared.u[k * E + e] = 0.4 * Math.exp(-k / 3) * Math.cos(latEdge[e]) * Math.sin(3 * e);
  const total = (field) => { let q = 0; for (let n = 0; n < field.length; n++) q += areaCell[n % C] * field[n]; return q; };
  for (let n = 0; n < 2; n++) for (const phase of FULL_STEP) {
    const heat = total(sheared.Q), salt = total(sheared.W);
    sheared.phases[phase](DT);
    if (phase === 'rescale' || phase === 'shift') continue;
    assert.ok(Math.abs(total(sheared.Q) - heat) < 1e-12 * Math.abs(heat), `${phase}: heat ${heat} → ${total(sheared.Q)}`);
    assert.ok(Math.abs(total(sheared.W) - salt) < 1e-12 * Math.abs(salt), `${phase}: salt ${salt} → ${total(sheared.W)}`);
  }
  let largest = 0;
  for (const v of sheared.isoKappa) largest = Math.max(largest, v);
  assert.ok(largest > 0, 'no interface took a κS² diffusivity');
});

test('under isoneutralDiffusion an unstable stratification takes every slope at the cap, so the κS² diffusivity goes as the cap squared and the triads\' change to the tendency as the cap; the boundary taper scales an interface\'s diffusivity by the square of its depth over the boundary layer\'s', () => {
  const unstable = (depth, lat) => [280 + 0.01 * depth + 5 * Math.cos(lat) ** 2, 35];
  const flat = tracerTendency(profiled({}, unstable));
  const capped = [0.01, 0.02].map((slopeCap) => tracerTendency(profiled({ isoneutralDiffusion: { slopeCap, boundaryTaper: false } }, unstable)));
  let largest = 0, missQ = 0, miss33 = 0;
  for (let n = 0; n < flat.dQ.length; n++) { const change = capped[0].dQ[n] - flat.dQ[n]; largest = Math.max(largest, Math.abs(change)); missQ = Math.max(missQ, Math.abs(capped[1].dQ[n] - flat.dQ[n] - 2 * change)); }
  for (let n = 0; n < capped[0].k33.length; n++) miss33 = Math.max(miss33, Math.abs(capped[1].k33[n] - 4 * capped[0].k33[n]) / Math.max(capped[1].k33[n], 1e-300));
  console.log(`doubling the cap: the triads' change ${largest.toExponential(3)}, its miss of twice ${missQ.toExponential(3)}; the κS² diffusivity's relative miss of four times ${miss33.toExponential(3)}`);
  assert.ok(largest > 0 && missQ < 1e-9 * largest && miss33 < 1e-12, `${largest}, ${missQ}, ${miss33}`);
  const stable = (depth, lat, i) => [288 - 0.01 * depth + 3 * Math.cos(lat) ** 2 + 0.5 * Math.sin(7 * mesh.lonCell[i]) * Math.exp(-depth / 150), 34.6 + 2e-4 * depth];
  const depthOf = (ocean, m, i) => { let z = 0; for (let j = 0; j < m; j++) z += ocean.h[j * C + i]; return z; };
  const untapered = profiled({ isoneutralDiffusion: { boundaryTaper: false } }, stable), F = untapered.table.fixed, base = tracerTendency(untapered).k33;
  for (const layer of [0, 5000]) {
    const tapered = profiled({ isoneutralDiffusion: true }, stable);
    tapered.boundaryDepths.fill(layer);
    const k33 = tracerTendency(tapered).k33;
    let miss = 0, positive = 0;
    for (let m = 1; m <= F; m++) for (let i = 0; i < C; i++) {
      const n = (m - 1) * C + i, share = layer > 0 ? (depthOf(tapered, m, i) / layer) ** 2 : 1;
      if (base[n] > 0) positive++;
      miss = Math.max(miss, Math.abs(k33[n] - base[n] * share) / Math.max(base[n] * share, 1e-300));
    }
    assert.ok(positive > 0.5 * F * C && miss < 1e-12, `boundary layer ${layer} m: ${positive} interfaces, miss ${miss}`);
  }
});
