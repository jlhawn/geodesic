import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean, backgroundDiffusivities, HENYEY_RATIO, EQUATORIAL_BACKGROUND, restratificationScale, restratificationOptions, deformationRadii, limitedFaceValue, limitedEdgeValue, eddyDiffusivities, shearResolution, shearResolutionOptions } from '../js/ocean/layered.module.js';
import { interiorMixing, shearTaper } from '../js/ocean/column.module.js';

/*
 * The ocean's optional interior closures on the CPU at N = 16, each against
 * its published form and against the default it leaves alone: the
 * latitude dependence of the background diffusivity (backgroundLatitude)
 * and the Richardson-number cutoff of the interior shear mixing
 * (shearCutoff), the restratification's equatorial taper
 * (mixedLayerRestratification.taperLatitude), the temperature and
 * salinity the eddy phase's layer fluxes carry (eddyTracers) and the
 * resolution function on the shear-scaled diffusion
 * (shearDiffusivityResolution).
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

test('shearTaper is Large, McWilliams and Doney\'s interior shape: 1 at Ri 0, (3/4)³ at half the cutoff, 0 at and above it', () => {
  assert.equal(shearTaper(0, 0.7), 1);
  assert.ok(Math.abs(shearTaper(0.35, 0.7) - 0.421875) < 1e-15);
  assert.equal(shearTaper(0.7, 0.7), 0);
  assert.equal(shearTaper(5, 0.7), 0);
  for (let r = 0.01; r < 0.7; r += 0.01) assert.ok(shearTaper(r, 0.7) < shearTaper(r - 0.01, 0.7));
});

test('under shearCutoff the interior keeps PP81\'s background and scales only its shear term by the taper, zero above the cutoff; without it the shear term is PP81\'s', () => {
  const K = 8, h = new Float64Array(K).fill(10), T = Float64Array.from({ length: K }, (_, k) => 293.15 - 0.5 * k), S = new Float64Array(K).fill(35);
  const u = Float64Array.of(1, 0.75, 0.55, 0.4, 0.3, 0.25, 0.24, 0.235), v = new Float64Array(K);
  const scratch = { rho: new Float64Array(K) }, nu0 = 5e-3, nub = 1e-4, kb = 1e-5;
  const run = (extra) => { const kappa = new Float64Array(K), nu = new Float64Array(K); interiorMixing({ h, T, S, u, v, K, shearViscosity: nu0, backgroundViscosity: nub, backgroundDiffusivity: kb, ...extra }, kappa, nu, scratch); return { kappa, nu }; };
  const off = run({}), infinite = run({ shearCutoff: Infinity }), cut = run({ shearCutoff: 0.7 });
  assert.deepEqual(Array.from(infinite.nu), Array.from(off.nu));
  assert.deepEqual(Array.from(infinite.kappa), Array.from(off.kappa));
  let below = 0, above = 0;
  for (let k = 0; k < K - 1; k++) {
    const ri = (Math.sqrt(nu0 / (off.nu[k] - nub)) - 1) / 5;
    if (ri >= 0.7 * (1 + 1e-9)) {
      above++;
      assert.equal(cut.nu[k], nub, `interface ${k} at Ri ${ri}`);
      assert.ok(Math.abs(cut.kappa[k] - (nub / (1 + 5 * ri) + kb)) < 1e-12 * cut.kappa[k], `interface ${k} at Ri ${ri}`);
    } else if (ri < 0.7 * (1 - 1e-6)) {
      below++;
      const ratio = (cut.nu[k] - nub) / (off.nu[k] - nub);
      assert.ok(Math.abs(ratio - (1 - (ri / 0.7) ** 2) ** 3) < 1e-6, `interface ${k} at Ri ${ri}: ${ratio}`);
      assert.ok(Math.abs(cut.kappa[k] - (cut.nu[k] / (1 + 5 * ri) + kb)) < 1e-9 * cut.kappa[k], `interface ${k}`);
    }
  }
  assert.ok(below >= 2 && above >= 2, `${below} interfaces below the cutoff and ${above} above`);
  assert.throws(() => createOcean(mesh, { shearCutoff: 0 }), /shearCutoff/);
});

test('shearCutoff reaches the ocean\'s column phase: under a sheared surface current every interface\'s viscosity is at most the default\'s and some are less', () => {
  const sheared = (options) => {
    const ocean = stratified(options), L = ocean.layers;
    for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) ocean.u[k * E + e] = 0.4 * Math.exp(-k / 3) * Math.cos(latEdge[e]) * Math.sin(3 * e);
    ocean.phases.column(DT);
    return ocean;
  };
  const base = sheared({}), cut = sheared({ shearCutoff: 0.7 });
  let less = 0;
  for (let n = 0; n < base.viscosity.length; n++) {
    assert.ok(cut.viscosity[n] <= base.viscosity[n] * (1 + 1e-12), `entry ${n}: ${cut.viscosity[n]} over ${base.viscosity[n]}`);
    if (cut.viscosity[n] < base.viscosity[n] * (1 - 1e-6)) less++;
  }
  assert.ok(less > 0, 'the cutoff changed no viscosity');
});

test('taperLatitude scales the restratification by |f|/√(f² + f_t²): none without it, half at the equator\'s first edges of a wide taper, 1/√2 at the taper latitude', () => {
  const plain = restratificationScale(mesh, restratificationOptions(true)), radii = deformationRadii(mesh);
  for (let e = 0; e < E; e++) assert.equal(plain[e], 1 * Math.min(dcEdge[e], radii[e]) / 500);
  const tapered = restratificationScale(mesh, restratificationOptions({ taperLatitude: 5 })), ft = 2 * mesh.omega * Math.sin(5 * DEGREE);
  let equator = Infinity;
  for (let e = 0; e < E; e++) {
    const expected = Math.abs(fEdge[e]) / Math.sqrt(fEdge[e] ** 2 + ft ** 2);
    assert.ok(Math.abs(tapered[e] / plain[e] - expected) < 1e-12, `edge ${e}: ${tapered[e] / plain[e]} against ${expected}`);
    equator = Math.min(equator, Math.abs(latEdge[e]));
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
