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

test('isoneutralDiffusion takes true for a taper centered on a slope of 0.01, 0.001 wide, under the boundary taper, refuses a negative or infinite cap, a width that is not positive and finite and a non-boolean taper, and under slopeCap 0 is the horizontal diffusion, a step of the ocean bit for bit', () => {
  assert.equal(isoneutralOptions(false), null);
  assert.deepEqual(isoneutralOptions(true), { slopeCap: 0.01, slopeWidth: 0.001, boundaryTaper: true });
  assert.deepEqual(ISONEUTRAL, { slopeCap: 0.01, slopeWidth: 0.001, boundaryTaper: true });
  assert.deepEqual(isoneutralOptions({ slopeCap: 0.02 }), { slopeCap: 0.02, slopeWidth: 0.001, boundaryTaper: true });
  for (const bad of [{ slopeCap: -1 }, { slopeCap: Infinity }, { slopeCap: NaN }, { slopeWidth: 0 }, { slopeWidth: -1e-3 }, { slopeWidth: Infinity }, { boundaryTaper: 1 }]) assert.throws(() => createOcean(mesh, { isoneutralDiffusion: bad }), /isoneutralDiffusion takes a non-negative slopeCap/);
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

test('under isoneutralDiffusion the triads move through the interior of a stratification whose temperature and salinity are linear in one tilted coordinate no more heat or salt than the 2e-9 of the horizontal diffusion that the slope taper leaves at so gentle a slope, where the horizontal diffusion crosses it whole', () => {
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
    assert.ok(scale > 0 && residual < 1e-8 * scale, `${label}: ${residual} against ${scale}`);
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

test('under isoneutralDiffusion a neutral or unstable interface takes no κS² diffusivity and turns no flux, however small its horizontal gradient, and each triad of a uniformly tilted stratification takes Danabasoglu and McWilliams\'s tanh taper of κ on its slope terms, a κS² diffusivity of κ|S|² where the tilt is gentle', () => {
  const iso = { boundaryTaper: false };
  for (const [label, profile] of [['unstable', (depth, lat) => [280 + 0.01 * depth + 5 * Math.cos(lat) ** 2, 35]], ['neutral', (depth, lat) => [285 + 1e-9 * Math.cos(lat) ** 2, 35]]]) {
    const flat = tracerTendency(profiled({}, profile)), turned = tracerTendency(profiled({ isoneutralDiffusion: iso }, profile));
    let largest = 0, differ = 0;
    for (const v of turned.k33) largest = Math.max(largest, Math.abs(v));
    for (const f of ['dQ', 'dW']) for (let n = 0; n < flat[f].length; n++) if (flat[f][n] !== turned[f][n]) differ++;
    assert.ok(largest === 0 && differ === 0, `${label}: κS² up to ${largest} m²/s, and ${differ} tendencies part from the horizontal diffusion's`);
  }
  const { dvEdge, cellsOnEdge } = mesh, radius = 6.371e6, speed = 0.005, width = 0.001;
  const perEdge = (values) => { const out = new Float64Array(C); for (let e = 0; e < E; e++) { out[cellsOnEdge[2 * e]] += values[e]; out[cellsOnEdge[2 * e + 1]] += values[e]; } return out; };
  const across = perEdge(Float64Array.from(dvEdge, (dv, e) => dv * dcEdge[e])), kappaMean = perEdge(Float64Array.from(dvEdge, (dv, e) => speed * dcEdge[e] * dv * dcEdge[e]));
  for (const [tilt, cap] of [[1e-3, 0.01], [0.01, 0.01], [0.01, 0.005]]) {
    const profile = (depth, lat) => { const x = depth + tilt * radius * lat; return [285 - 2e-5 * x, 35 + 2e-6 * x]; };
    const ocean = profiled({ isoneutralDiffusion: { ...iso, slopeCap: cap }, lateralDiffusionSpeed: speed }, profile), F = ocean.table.fixed, k33 = tracerTendency(ocean).k33;
    const want = perEdge(Float64Array.from(dcEdge, (dc, e) => {
      const slope = tilt * radius * (latCell[cellsOnEdge[2 * e + 1]] - latCell[cellsOnEdge[2 * e]]) / dc;
      return 0.5 * speed * dc * dc * dvEdge[e] * 0.5 * (1 + Math.tanh((cap - Math.abs(slope)) / width)) * slope * slope;
    }));
    let miss = 0;
    const ratios = [];
    for (let i = 0; i < C; i++) {
      for (let m = 1; m < F; m++) miss = Math.max(miss, Math.abs(k33[(m - 1) * C + i] * areaCell[i] - want[i]) / want[i]);
      if (Math.abs(latCell[i]) < 60 * DEGREE) ratios.push(k33[C + i] / (kappaMean[i] / across[i] * tilt * tilt));
    }
    ratios.sort((x, y) => x - y);
    const median = ratios[ratios.length >> 1];
    console.log(`tilt ${tilt}, slopeCap ${cap}: κS² against Σ ½κ dv dc f(S) S² over the cell's edges, largest relative miss ${miss.toExponential(2)}; median κS²/(κ̄ tilt²) within 60° ${median.toFixed(4)}`);
    assert.ok(miss < 1e-7, `tilt ${tilt}, slopeCap ${cap}: κS² parts from the tapered triads' by ${miss}`);
    if (tilt < 0.1 * cap) assert.ok(median > 0.9 && median < 1.1, `tilt ${tilt}: κS² is ${median} of κ̄|S|²`);
  }
});

test('under the boundary taper a triad above the base of its cell\'s boundary layer takes the slope at that base times its depth over the layer\'s, so a mixed boundary layer\'s κS² diffusivity follows the stratification below it and not slopeCap, and with no boundary layer the taper changes nothing', () => {
  const D = 30, profile = (depth, lat) => [288 + 5 * Math.cos(lat) ** 2 - (depth > D ? 1e-6 * D + 0.02 * (depth - D) : 1e-6 * depth), 35];
  const run = (iso, layer) => { const ocean = profiled({ isoneutralDiffusion: iso }, profile); ocean.boundaryDepths.fill(layer); return { ocean, k33: tracerTendency(ocean).k33 }; };
  const { ocean, k33 } = run(true, D), wide = run({ slopeCap: 0.02 }, D).k33, F = ocean.table.fixed;
  const depthOf = (m, i) => { let z = 0; for (let j = 0; j < m; j++) z += ocean.h[j * C + i]; return z; };
  let inside = 0, capMiss = 0, profileMiss = 0, jumpLow = Infinity, jumpHigh = 0, deepest = 0;
  for (let i = 0; i < C; i++) {
    let scale = 0, m = 1;
    for (; m < F && depthOf(m, i) < D; m++) {
      const n = (m - 1) * C + i, z = depthOf(m, i);
      if (!(k33[n] > 0)) continue;
      inside++;
      deepest = Math.max(deepest, k33[n]);
      capMiss = Math.max(capMiss, Math.abs(wide[n] - k33[n]) / k33[n]);
      const flat = k33[n] * (D / z) ** 2;
      if (scale === 0) scale = flat; else profileMiss = Math.max(profileMiss, Math.abs(flat - scale) / scale);
    }
    if (scale > 0) { const jump = scale / k33[(m - 1) * C + i]; jumpLow = Math.min(jumpLow, jump); jumpHigh = Math.max(jumpHigh, jump); }
  }
  console.log(`${inside} interfaces above a ${D} m boundary layer's base, κS² up to ${deepest.toExponential(3)} m²/s: doubling slopeCap moves them by ${capMiss.toExponential(2)}, κS² (D/z)² parts within a column by ${profileMiss.toExponential(2)}, and over the interface at the base by a factor ${jumpLow.toFixed(4)} to ${jumpHigh.toFixed(4)}`);
  assert.ok(inside > 0.5 * C, `${inside} interfaces`);
  assert.ok(capMiss < 1e-6, `doubling slopeCap moves the boundary layer's κS² by ${capMiss}`);
  assert.ok(profileMiss < 1e-9, `the boundary layer's κS² parts from the base slope's (z/D)² by ${profileMiss}`);
  assert.ok(jumpLow > 0.8 && jumpHigh < 1.25, `κS² over the boundary layer's base jumps by ${jumpLow} to ${jumpHigh}`);
  assert.deepEqual(Array.from(run(true, 0).k33), Array.from(run({ boundaryTaper: false }, 0).k33));
});

test('under isoneutralDiffusion the column takes the κS² diffusivity of the dynamics\' four stages weighted 1:2:2:1, as their cross terms are', () => {
  const profile = (depth, lat, i) => [288 - 0.01 * depth + 3 * Math.cos(lat) ** 2 + 0.5 * Math.sin(7 * mesh.lonCell[i]) * Math.exp(-depth / 150), 34.6 + 2e-4 * depth];
  const ocean = profiled({ isoneutralDiffusion: true, tracerAdvection: 'upwind', shearDiffusivity: 0 }, profile), L = ocean.layers;
  ocean.boundaryDepths.fill(30);
  for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) ocean.u[k * E + e] = 0.4 * Math.exp(-k / 3) * Math.cos(latEdge[e]) * Math.sin(3 * e);
  const start = ocean.state.map((a) => Float64Array.from(a));
  ocean.phases.dynamics(DT);
  const step = Float64Array.from(ocean.isoKappa), rates = ocean.stages.map((stage) => stage.map((a) => Float64Array.from(a)));
  const input = (s) => start.map((a, f) => (s === 0 ? Float64Array.from(a) : Float64Array.from(a, (v, n) => v + [0, DT / 2, DT / 2, DT][s] * rates[s - 1][f][n])));
  const inputs = [0, 1, 2, 3].map(input);
  let trialMiss = 0;
  for (const f of [0, 2, 3]) for (let n = 0; n < inputs[3][f].length; n++) trialMiss = Math.max(trialMiss, Math.abs(inputs[3][f][n] - ocean.trial[f][n]));
  const scratch = ocean.stages[0].map((a) => new Float64Array(a.length));
  const stageKappa = inputs.map((x) => { ocean.tendency(x, scratch); return Float64Array.from(ocean.isoKappa); });
  let miss = 0, apart = 0, largest = 0;
  for (let n = 0; n < step.length; n++) {
    const mean = (stageKappa[0][n] + 2 * stageKappa[1][n] + 2 * stageKappa[2][n] + stageKappa[3][n]) / 6;
    largest = Math.max(largest, mean);
    miss = Math.max(miss, Math.abs(step[n] - mean)); apart = Math.max(apart, Math.abs(step[n] - stageKappa[3][n]));
  }
  console.log(`the fourth stage's state rebuilt to ${trialMiss.toExponential(2)}; the step's κS² against the 1:2:2:1 mean of its stages' ${(miss / largest).toExponential(2)} and against the fourth stage's alone ${(apart / largest).toExponential(2)}, of the largest ${largest.toExponential(3)} m²/s`);
  assert.ok(trialMiss === 0, `the fourth stage's state rebuilt to ${trialMiss}`);
  assert.ok(largest > 0 && miss < 1e-12 * largest, `the step's κS² parts from its stages' 1:2:2:1 mean by ${miss} of ${largest}`);
  assert.ok(apart > 1e-6 * largest, `the stages' κS² differ by only ${apart} of ${largest}`);
});

const SHELF_MINIMUM = 0.5;
const shelfSea = (boundary) => {
  const { lonCell, cellsOnEdge } = mesh, distance = Float64Array.from(latCell, (lat, i) => Math.acos(Math.min(1, Math.cos(lat) * Math.cos(lonCell[i]))) / DEGREE);
  const land = Uint8Array.from(distance, (d) => (d < 7 ? 1 : 0)), elevation = Float64Array.from(distance, (d) => (d < 7 ? 100 : -Math.min(4000, 40 * 1.6 ** ((d - 7) / 4))));
  const edgeOcean = Uint8Array.from({ length: E }, (_, e) => (land[cellsOnEdge[2 * e]] || land[cellsOnEdge[2 * e + 1]] ? 0 : 1));
  const ocean = createOcean(mesh, { geography: { land, elevation, edgeOcean }, everySteps: 1, iceDrift: false, minimumLayer: SHELF_MINIMUM, isoneutralDiffusion: true });
  ocean.initialize(Float64Array.from(latCell, (lat) => 271.5 + 30 * Math.cos(lat) ** 2), new Float64Array(C));
  const { h, Q, W } = ocean;
  ocean.u.fill(0);
  for (let i = 0; i < C; i++) {
    let z = 0;
    for (let k = 0; k < ocean.layers; k++) {
      const n = k * C + i, hv = h[n], mid = z + 0.5 * hv;
      z += hv;
      Q[n] = hv * (288 - 0.02 * mid + 2 * Math.cos(latCell[i]) ** 2 + 0.3 * Math.sin(5 * lonCell[i]));
      W[n] = hv * (34.6 + 4e-3 * mid + 0.2 * Math.sin(latCell[i]));
    }
  }
  ocean.boundaryDepths.fill(boundary);
  const sea = [];
  for (let i = 0; i < C; i++) if (!land[i]) sea.push(i);
  return { ocean, land, sea, thin: 2 * SHELF_MINIMUM };
};
const movedRows = (before, after, rows) => {
  let moved = 0;
  for (const k of rows) for (let i = 0; i < C; i++) if (before[k * C + i] !== after[k * C + i]) moved++;
  return moved;
};
const thinWater = (ocean, cells, thin) => {
  let changed = 0;
  for (const i of cells) for (let k = 0; k < ocean.layers; k++) { const n = k * C + i; if (ocean.h[n] <= thin) { ocean.Q[n] += 0.3 * ocean.h[n]; ocean.W[n] += 0.05 * ocean.h[n]; changed++; } }
  return changed;
};

test('under isoneutralDiffusion no triad reads a land cell: water in the land cells\' layers leaves every sea cell\'s tracer tendencies and κS² diffusivity bit for bit, and the land cells take none', () => {
  const { ocean, land } = shelfSea(10), L = ocean.layers, every = Array.from({ length: L }, (_, k) => k), { cellsOnEdge } = mesh;
  const before = tracerTendency(ocean);
  let coast = 0;
  for (let e = 0; e < E; e++) if (land[cellsOnEdge[2 * e]] !== land[cellsOnEdge[2 * e + 1]]) coast++;
  for (let i = 0; i < C; i++) if (land[i]) for (let k = 0; k < L; k++) { const n = k * C + i; ocean.h[n] = 3; ocean.Q[n] = 3 * (289 + 0.01 * k); ocean.W[n] = 3 * (34.7 + 0.01 * k); }
  const after = tracerTendency(ocean);
  let landTendency = 0;
  for (const run of [before, after]) for (let i = 0; i < C; i++) if (land[i]) for (let k = 0; k < L; k++) landTendency = Math.max(landTendency, Math.abs(run.dQ[k * C + i]), Math.abs(run.dW[k * C + i]));
  const moved = movedRows(before.dQ, after.dQ, every) + movedRows(before.dW, after.dW, every) + movedRows(before.k33, after.k33, every.slice(0, before.k33.length / C));
  console.log(`${coast} coastal edges: water in the land cells moves ${moved} sea tendencies and κS² values, and the land cells' largest tendency is ${landTendency}`);
  assert.ok(coast > 0, 'no coast');
  assert.equal(moved, 0);
  assert.equal(landTendency, 0);
});

test('under the boundary taper a column whose thick layers end above its boundary layer\'s base takes no κS² diffusivity, and no triad takes its slope from a base beside a thin layer, in its own cell or the other cell of its edge: the water of every thin layer moves no column\'s top six layers or their κS²', () => {
  const { ocean, sea, thin } = shelfSea(30), F = ocean.table.fixed, { h } = ocean, top = [0, 1, 2, 3, 4, 5];
  const before = tracerTendency(ocean);
  let shelves = 0, shelfKappa = 0;
  for (const i of sea) {
    let z = 0, k = 0;
    while (k < F && h[k * C + i] > thin) z += h[k++ * C + i];
    if (z >= 30) continue;
    shelves++;
    for (let m = 0; m < F; m++) shelfKappa = Math.max(shelfKappa, before.k33[m * C + i]);
  }
  const changed = thinWater(ocean, sea, thin), after = tracerTendency(ocean);
  const moved = movedRows(before.dQ, after.dQ, top) + movedRows(before.dW, after.dW, top) + movedRows(before.k33, after.k33, top);
  console.log(`${shelves} columns thick to less than their 30 m boundary layer take κS² up to ${shelfKappa}; the water of ${changed} thin layers moves ${moved} tendencies and κS² values in the top six layers`);
  assert.ok(shelves > 0 && changed > 0);
  assert.equal(shelfKappa, 0);
  assert.equal(moved, 0);
});

test('under isoneutralDiffusion a layer no thicker than twice minimumLayer sets no slope and takes no triad flux: no interface beside one takes κS² diffusivity, and the water of a shelf column\'s thin layers moves none of its layers thick in every column, nor their κS²', () => {
  const { ocean, sea, thin } = shelfSea(10), F = ocean.table.fixed, { h } = ocean;
  const before = tracerTendency(ocean);
  let beside = 0, besideKappa = 0, thick = 0;
  for (const i of sea) for (let m = 1; m <= F; m++) if (!(h[(m - 1) * C + i] > thin && h[m * C + i] > thin)) { beside++; besideKappa = Math.max(besideKappa, before.k33[(m - 1) * C + i]); }
  while (thick < F && sea.every((i) => h[thick * C + i] > thin)) thick++;
  const shelf = sea.filter((i) => h[thick * C + i] <= thin), changed = thinWater(ocean, shelf, thin), after = tracerTendency(ocean);
  let moved = 0;
  for (let k = 0; k < thick; k++) for (const i of shelf) { const n = k * C + i; if (before.dQ[n] !== after.dQ[n] || before.dW[n] !== after.dW[n] || before.k33[n] !== after.k33[n]) moved++; }
  console.log(`${beside} interfaces beside a thin layer take κS² up to ${besideKappa}; the water of ${changed} thin layers in ${shelf.length} shelf columns moves ${moved} of the ${thick * shelf.length} values in their ${thick} layers thick everywhere`);
  assert.ok(beside > 0 && shelf.length > 0 && thick > 0);
  assert.equal(besideKappa, 0);
  assert.equal(moved, 0);
});
