import { test } from 'node:test';
import assert from 'node:assert/strict';
import { cellVector } from '../js/dynamics/operators.module.js';
import { createOcean, EPS, EDDY_SLACK } from '../js/ocean/layered.module.js';
import { mesh, C, E, DEG, RHO, PLAIN_MIXING, uniformOcean, zonalWindOnEdges, deepen } from './helpers/layered.mjs';

const DT = 2700, DEPTH = 4000, DRAG = 2e-4, WIND_DEPTH = 50, OMEGA = 7.292e-5;

test('a uniform zonal stress drives an Ekman transport of τ/(ρ₀f) to the right of the wind in the north and to the left in the south, and downwind on the equator, and the formula holds at every edge', () => {
  const { ocean, ice } = uniformOcean({ ekmanBlock: true });
  const tau = 0.1, stress = zonalWindOnEdges(mesh, () => tau);
  ocean.setStress(stress, ice);
  ocean.ekmanTransport(DT);
  const vector = cellVector(mesh, ocean.ekmanFlux), rows = [];
  for (const [low, high] of [[-60, -15], [-15, -3], [-3, 3], [3, 15], [15, 60]]) {
    let area = 0, north = 0, east = 0, expected = 0;
    for (let i = 0; i < C; i++) {
      const lat = mesh.latCell[i] / DEG;
      if (lat < low || lat >= high) continue;
      const x = mesh.xCell[3 * i], y = mesh.xCell[3 * i + 1], z = mesh.xCell[3 * i + 2], r = Math.hypot(x, y), radius = Math.hypot(x, y, z);
      const a = mesh.areaCell[i], f = 2 * OMEGA * Math.sin(mesh.latCell[i]);
      area += a;
      north += a * (vector[3 * i] * (-z * x / r / radius) + vector[3 * i + 1] * (-z * y / r / radius) + vector[3 * i + 2] * (r / radius));
      east += a * (vector[3 * i] * (-y / r) + vector[3 * i + 1] * (x / r));
      expected += a * (-tau / (RHO * f));
    }
    rows.push({ low, high, north: north / area, east: east / area, expected: expected / area });
  }
  console.log(`Ekman transport for 0.1 N/m² westerly: ${rows.map((r) => `${r.low}..${r.high}°: ${r.north.toFixed(3)} north, ${r.east.toFixed(3)} east (τ/ρ₀f ${r.expected.toFixed(3)})`).join('; ')} m²/s`);
  const [south, southTropic, equator, northTropic, north] = rows;
  assert.ok(Math.abs(north.north / north.expected - 1) < 0.02 && north.north < 0, `${north.north} against ${north.expected}`);
  assert.ok(Math.abs(south.north / south.expected - 1) < 0.02 && south.north > 0, `${south.north} against ${south.expected}`);
  assert.ok(northTropic.north < 0 && southTropic.north > 0 && Math.abs(northTropic.north + southTropic.north) < 1e-9);
  assert.ok(equator.east > 5 && Math.abs(equator.north) < 1e-9 * equator.east, `${equator.east} east, ${equator.north} north on the equator`);
  let nearest = 0;
  for (let e = 0; e < E; e++) {
    const expected = transportOf({ stress }, e, WIND_DEPTH);
    assert.ok(Math.abs(ocean.ekmanFlux[e] - expected) <= 1e-12 * Math.abs(expected), `edge ${e}: ${ocean.ekmanFlux[e]} against ${expected}`);
    if (Math.abs(mesh.fEdge[e]) < Math.abs(mesh.fEdge[nearest])) nearest = e;
  }
  const r = DRAG / WIND_DEPTH, friction = stress[nearest] * WIND_DEPTH / (RHO * DRAG);
  assert.ok(Math.abs(mesh.fEdge[nearest]) < 0.25 * r && Math.abs(ocean.ekmanFlux[nearest] - friction) < 0.1 * Math.abs(friction) + Math.abs(stress[nearest]) * 0, `edge ${nearest}: ${ocean.ekmanFlux[nearest]} against ${friction} (f ${mesh.fEdge[nearest]}, r ${r})`);
});

let pairEdge = 0;
for (let e = 1; e < E; e++) if (Math.abs(mesh.latEdge[e] - 35 * DEG) < Math.abs(mesh.latEdge[pairEdge] - 35 * DEG)) pairEdge = e;

/*
 * A sea of two cells, a and b, joined by one edge in mid-latitudes, with
 * the stress on that edge alone, so that its transport is the normal
 * stress's, and a wind depth of 50 m: every class a token but the
 * mixed layer, the lid (class 1) and the deep class (2) that fills each
 * column to 4000 m.
 */
function pair(options = {}) {
  const e = pairEdge, a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
  const land = new Uint8Array(C).fill(1), edgeOcean = new Uint8Array(E);
  land[a] = 0; land[b] = 0; edgeOcean[e] = 1;
  const geography = { land, elevation: new Float64Array(C).fill(-DEPTH), edgeOcean };
  const ocean = createOcean(mesh, { everySteps: 1, geography, bathymetry: new Float64Array(C).fill(DEPTH), ekmanBlock: true, ekmanDepth: 0, thermoclineTilt: 0, salinityProfile: () => 35, ...options });
  ocean.initialize(new Float64Array(C).fill(288), new Float64Array(C));
  return { ocean, e, a, b };
}
function setColumn(ocean, i, { layer, lid, temperatures: [layerT, lidT, deepT] }) {
  const { h, Q, W } = ocean, L = ocean.layers;
  const put = (k, thickness, t, s = 35) => { const n = k * C + i; h[n] = thickness; Q[n] = thickness * t; W[n] = thickness * s; };
  for (let k = 3; k < L; k++) put(k, EPS, 275, 35);
  put(0, layer, layerT, 34.5);
  put(1, lid, lidT, 35.2);
  put(2, DEPTH - layer - lid - (L - 3) * EPS, deepT, 34.9);
  ocean.eta[i] = 0;
}

/*
 * The edge's Ekman transport written out from the mesh: the stress's normal
 * and tangential parts at the edge, f and the drag over the reach.
 */
function transportOf(ocean, e, reach) {
  let tangential = 0;
  for (let s = 0; s < mesh.nEdgesOnEdge[e]; s++) {
    const o = mesh.edgesOnEdge[mesh.maxEdgesOnEdge * e + s];
    tangential += mesh.weightsOnEdge[mesh.maxEdgesOnEdge * e + s] * mesh.dvEdge[o] * ocean.stress[o];
  }
  tangential /= mesh.dcEdge[e];
  const f = mesh.fEdge[e], r = DRAG / reach;
  return (r * ocean.stress[e] + f * tangential) / (RHO * (r * r + f * f));
}

/*
 * The change of the three stores the pass made at a and b, as the volume
 * flux from a to b (m³/s) that would have made it, per class, and the
 * flux's heat and salt from the donors the change names.
 */
function changes(ocean, before, a, b, dt) {
  const L = ocean.layers, flux = new Float64Array(L), heat = new Float64Array(L), salt = new Float64Array(L);
  for (let k = 0; k < L; k++) {
    flux[k] = (ocean.h[k * C + b] - before.h[k * C + b]) * mesh.areaCell[b] / dt;
    heat[k] = (ocean.Q[k * C + b] - before.Q[k * C + b]) * mesh.areaCell[b] / dt;
    salt[k] = (ocean.W[k * C + b] - before.W[k * C + b]) * mesh.areaCell[b] / dt;
  }
  return { flux, heat, salt };
}
function differing(left, right) {
  let count = 0;
  left.forEach((a, n) => { for (let i = 0; i < a.length; i++) if (!Object.is(a[i], right[n][i])) count++; });
  return count;
}
const snapshot = (ocean) => ({ h: Float64Array.from(ocean.h), Q: Float64Array.from(ocean.Q), W: Float64Array.from(ocean.W) });

/*
 * What the design says one edge does, from its donor's column: the stress
 * shares σ of the edge's own layers over its reach, the block shares ω of
 * the donor's top max(wind depth, h₀), and U_E dv (ω − σ) per class.
 */
function expectedFlux(ocean, e, reach, wind, stressThickness = 20) {
  const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1], L = ocean.layers;
  const column = Array.from({ length: L }, (_, k) => (k === 0 ? 0.5 * (ocean.h[a] + ocean.h[b]) : Math.min(ocean.h[k * C + a], ocean.h[k * C + b])));
  const shares = new Array(L).fill(0);
  let top = 0, depth = 0;
  for (let k = 0; k < L; k++) {
    const part = Math.min(column[k], Math.max(0, reach - top));
    if (k === 0 || column[k] >= stressThickness) { shares[k] = part; depth += part; }
    top += column[k];
    if (top >= reach) break;
  }
  const transport = transportOf(ocean, e, reach), donor = transport > 0 ? a : b;
  const block = new Array(L).fill(0), blockDepth = Math.max(wind, ocean.h[donor]);
  let z = 0, total = 0;
  for (let k = 0; k < L; k++) {
    const hk = ocean.h[k * C + donor];
    block[k] = Math.min(k === 0 ? hk : Math.max(0, hk - EPS), Math.max(0, blockDepth - z));
    total += block[k]; z += hk;
  }
  const volume = transport * mesh.dvEdge[e];
  return { transport, volume, omega: block.map((o) => o / total), sigma: shares.map((s) => s / depth), flux: block.map((o, k) => volume * (o / total - shares[k] / depth)), donor };
}

for (const sign of [1, -1]) {
  test(`${sign > 0 ? 'a → b' : 'b → a'}: a coastal lid under a thin mixed layer crosses with the Ekman transport in the donor's proportions, the pass keeping each class's volume, heat and salt, the column's volume and the upwind cell's temperatures`, () => {
    const { ocean, e, a, b } = pair();
    const coast = sign > 0 ? a : b, offshore = sign > 0 ? b : a;
    ocean.stress[e] = 0.2 * sign;
    const cases = [
      { name: 'mixed layer thick offshore, stress over the layer alone', coast: { layer: 20, lid: 15 + EPS }, offshore: { layer: 80, lid: EPS }, headline: true },
      { name: 'both mixed layers 20 m, stress spread to the deep class', coast: { layer: 20, lid: 15 }, offshore: { layer: 20, lid: EPS }, headline: false },
    ];
    for (const c of cases) {
      setColumn(ocean, coast, { ...c.coast, temperatures: [293.15, 293.15, 283.15] });
      setColumn(ocean, offshore, { ...c.offshore, temperatures: [298.15, 275, 281.15] });
      const before = snapshot(ocean), plan = expectedFlux(ocean, e, WIND_DEPTH, WIND_DEPTH);
      ocean.ekmanTransport(DT);
      assert.ok(Math.abs(ocean.ekmanFlux[e] - plan.transport) <= 1e-12 * Math.abs(plan.transport), `${ocean.ekmanFlux[e]} against ${plan.transport}`);
      assert.equal(plan.donor, coast);
      const { flux, heat, salt } = changes(ocean, before, a, b, DT);
      const moved = Math.abs(plan.volume);
      const rounding = (k, store, scale) => 8 * Number.EPSILON * scale * (Math.abs(store[k * C + a]) * mesh.areaCell[a] + Math.abs(store[k * C + b]) * mesh.areaCell[b]);
      let lost = 0, gained = 0, slack = 0;
      for (let k = 0; k < ocean.layers; k++) {
        const tolerance = 1e-9 * moved + rounding(k, before.h, 1) / DT;
        slack += rounding(k, before.h, 1) / DT;
        assert.ok(Math.abs(flux[k] - plan.flux[k]) <= tolerance, `${c.name}: class ${k} crossed ${flux[k]} m³/s, expected ${plan.flux[k]}`);
        const source = plan.flux[k] > 0 ? a : b, t = before.Q[k * C + source] / before.h[k * C + source], s = before.W[k * C + source] / before.h[k * C + source];
        assert.ok(Math.abs(heat[k] - plan.flux[k] * t) <= 1e-9 * moved * 300 + rounding(k, before.Q, 1) / DT && Math.abs(salt[k] - plan.flux[k] * s) <= 1e-9 * moved * 35 + rounding(k, before.W, 1) / DT, `${c.name}: class ${k} heat ${heat[k]} salt ${salt[k]}`);
        const dha = (ocean.h[k * C + a] - before.h[k * C + a]) * mesh.areaCell[a], dhb = (ocean.h[k * C + b] - before.h[k * C + b]) * mesh.areaCell[b];
        const dqa = (ocean.Q[k * C + a] - before.Q[k * C + a]) * mesh.areaCell[a], dqb = (ocean.Q[k * C + b] - before.Q[k * C + b]) * mesh.areaCell[b];
        const dwa = (ocean.W[k * C + a] - before.W[k * C + a]) * mesh.areaCell[a], dwb = (ocean.W[k * C + b] - before.W[k * C + b]) * mesh.areaCell[b];
        assert.ok(Math.abs(dha + dhb) <= rounding(k, before.h, 1) && Math.abs(dqa + dqb) <= rounding(k, before.Q, 1) && Math.abs(dwa + dwb) <= rounding(k, before.W, 1), `${c.name}: class ${k} is not kept: ${dha + dhb} m³, ${dqa + dqb}, ${dwa + dwb}`);
        if (k >= 3) assert.ok(ocean.h[k * C + a] === before.h[k * C + a] && ocean.h[k * C + b] === before.h[k * C + b], `${c.name}: token class ${k} changed`);
        lost += flux[k];
        gained += Math.abs(flux[k]);
      }
      assert.ok(Math.abs(lost) <= slack, `${c.name}: the classes' fluxes sum to ${lost} m³/s`);
      for (const cell of [a, b]) {
        let column = 0;
        for (let k = 0; k < ocean.layers; k++) column += ocean.h[k * C + cell] - before.h[k * C + cell];
        assert.ok(Math.abs(column) <= 1e-9, `${c.name}: the column of cell ${cell} changed by ${column} m`);
      }
      assert.ok(Math.abs(ocean.ekmanMoved.reduce((s, x) => s + x, 0) - gained * 86400 / 1e9) <= 1e-9 * gained * 86400 / 1e9);
      const share = (k) => plan.flux[k] / plan.volume;
      console.log(`${c.name} (${sign > 0 ? 'a → b' : 'b → a'}): the layer, lid and deep class change by ${[0, 1, 2].map((k) => share(k).toFixed(4)).join(', ')} of U_E dv; ω ${[0, 1, 2].map((k) => plan.omega[k].toFixed(4)).join(', ')}, σ ${[0, 1, 2].map((k) => plan.sigma[k].toFixed(4)).join(', ')}`);
      if (c.headline) {
        const lidShare = Math.abs(flux[1]) / moved, layerShare = Math.abs(flux[0]) / moved;
        assert.ok(Math.abs(lidShare - 15 / 50) < 1e-3 && Math.abs(flux[1] / plan.volume - 15 / 50) < 1e-3, `the lid crosses at ${lidShare} of U_E dv, not 15/B`);
        assert.ok(Math.abs(plan.omega[0] - 20 / 50) < 1e-3 && Math.abs(layerShare - (1 - 20 / 50)) < 1e-3, `the layer's share falls to ${plan.omega[0]}, its flux changing by ${layerShare}`);
        assert.ok(Math.sign(flux[1]) === Math.sign(plan.volume) && Math.sign(flux[0]) === -Math.sign(plan.volume), 'the lid crosses with the transport and the layer\'s water returns');
        const lidT = before.Q[C + coast] / before.h[C + coast], layerT = before.Q[offshore] / before.h[offshore];
        assert.ok(Math.abs(heat[1] / flux[1] - lidT) < 1e-9 * lidT && Math.abs(heat[0] / flux[0] - layerT) < 1e-9 * layerT, 'the lid carries the coast\'s temperature and the returning layer the offshore\'s');
      }
    }
  });
}

test('a donor class at its token thickness never carries water and one a little above it is held to 1/nEdges of its excess, the whole edge scaled together', () => {
  const columns = (layer, lid) => ({ coast: { layer, lid, temperatures: [293.15, 293.15, 283.15] }, offshore: { layer, lid: EPS, temperatures: [298.15, 275, 281.15] } });
  const run = (lid, dt, tau = 0.2) => {
    const { ocean, e, a, b } = pair();
    ocean.stress[e] = tau;
    const { coast, offshore } = columns(20, lid);
    setColumn(ocean, a, coast); setColumn(ocean, b, offshore);
    const before = snapshot(ocean), plan = expectedFlux(ocean, e, WIND_DEPTH, WIND_DEPTH);
    ocean.ekmanTransport(dt);
    return { ocean, a, b, e, before, plan, ...changes(ocean, before, a, b, dt), limited: ocean.diagnostics().oceanEkmanLimited };
  };

  const token = run(EPS, DT);
  assert.equal(token.ocean.h[C + token.a], token.before.h[C + token.a]);
  assert.equal(token.flux[1], 0);
  assert.ok(token.ocean.h[C + token.a] >= EPS && token.limited === 0);

  const gentle = run(0.2, DT), n = mesh.nEdgesOnCell[gentle.a];
  const excess = 0.2 - EPS;
  assert.equal(gentle.limited, 0);
  for (let k = 0; k < gentle.ocean.layers; k++) assert.ok(Math.abs(gentle.flux[k] - gentle.plan.flux[k]) <= 1e-8 * Math.abs(gentle.plan.volume));
  const dt = 3 * excess * mesh.areaCell[gentle.a] / (n * Math.abs(gentle.plan.flux[1]));
  const hard = run(0.2, dt);
  assert.equal(hard.limited, 1);
  const allowed = (k, cell) => (Math.max(0, hard.before.h[k * C + cell] - EPS) + EDDY_SLACK) * mesh.areaCell[cell] / (mesh.nEdgesOnCell[cell] * dt);
  let lambda = 1;
  for (let k = 0; k < hard.ocean.layers; k++) {
    const f = hard.plan.flux[k];
    if (f !== 0) lambda = Math.min(lambda, allowed(k, f > 0 ? hard.a : hard.b) / Math.abs(f));
  }
  console.log(`a lid 0.2 m thick, the pass held to λ = ${lambda.toFixed(4)} of its flux (the lid's allowance ${allowed(1, hard.a).toExponential(3)} m³/s against ${Math.abs(hard.plan.flux[1]).toExponential(3)}); the lid ends ${hard.ocean.h[C + hard.a].toFixed(4)} m thick`);
  assert.ok(lambda < 0.5 && lambda > 0);
  let sum = 0;
  for (let k = 0; k < hard.ocean.layers; k++) {
    assert.ok(Math.abs(hard.flux[k] - lambda * hard.plan.flux[k]) <= 1e-8 * Math.abs(hard.plan.flux[1]) + 1e-9 * Math.abs(hard.plan.volume) * lambda, `class ${k}: ${hard.flux[k]} against λ·${hard.plan.flux[k]}`);
    sum += hard.flux[k];
  }
  assert.ok(Math.abs(sum) <= 1e-8 * Math.abs(hard.plan.volume), `the held fluxes sum to ${sum}`);
  assert.ok(hard.ocean.h[C + hard.a] >= 0.2 - (excess + EDDY_SLACK) / n - 1e-12 && hard.ocean.h[C + hard.a] >= EPS);
  for (let n2 = 0; n2 < hard.ocean.h.length; n2++) assert.ok(hard.ocean.h[n2] >= Math.min(hard.before.h[n2], EPS) - 1e-12, `layer ${Math.floor(n2 / C)} of cell ${n2 % C} thinned to ${hard.ocean.h[n2]}`);
});

/*
 * Mixed layers of 50 to 160 m, a band of each along the meridians, their
 * surface a little lighter than the class under them.
 */
function deepSea(options = {}) {
  const sea = uniformOcean(options), band = [[], [], []];
  for (let i = 0; i < C; i++) band[Math.min(2, Math.floor(3 * (0.5 + 0.5 * Math.sin(3 * mesh.lonCell[i]))))].push(i);
  band.forEach((cells, n) => deepen(sea, cells, [50, 100, 160][n], -0.002));
  return sea;
}

test('with every mixed layer at least the wind depth the pass moves nothing, bit for bit, in a sea whose mixed layers differ from cell to cell; a step with the option off is the default step, and with it on differs from it only where a column dips a fraction of a meter below the wind depth', () => {
  const { ocean, ice } = deepSea({ ekmanBlock: true });
  const stress = zonalWindOnEdges(mesh, (lat) => -0.12 * Math.cos(3 * lat));
  ocean.setStress(stress, ice);
  let thinnest = Infinity, thickest = 0;
  for (let i = 0; i < C; i++) { thinnest = Math.min(thinnest, ocean.h[i]); thickest = Math.max(thickest, ocean.h[i]); }
  assert.ok(thinnest >= 50 - 1e-9 && thickest > 120, `mixed layers ${thinnest} to ${thickest} m`);
  const before = snapshot(ocean);
  ocean.ekmanTransport(DT);
  let strongest = 0;
  for (let e = 0; e < E; e++) strongest = Math.max(strongest, Math.abs(ocean.ekmanFlux[e]));
  console.log(`mixed layers ${thinnest.toFixed(1)} to ${thickest.toFixed(1)} m, the strongest Ekman transport ${strongest.toFixed(2)} m²/s, the pass moved ${ocean.ekmanMoved.reduce((s, x) => s + x, 0)} km³/day`);
  assert.ok(strongest > 1);
  assert.equal(differing([ocean.h, ocean.Q, ocean.W], [before.h, before.Q, before.W]), 0);
  assert.ok(ocean.ekmanMoved.every((x) => x === 0));

  const mixed = (options) => {
    const { ocean: o, surfaceT, ice: ic, flux } = deepSea(options);
    const wind = zonalWindOnEdges(mesh, (lat) => -0.12 * Math.cos(3 * lat));
    for (let n = 0; n < 6; n++) o.advance(surfaceT, ic, flux, wind, DT, undefined, (n + 1) * DT);
    return [o.h, o.u, o.Q, o.W, o.eta].map((x) => Float64Array.from(x));
  };
  const on = mixed({ ekmanBlock: true }), off = mixed({ ekmanBlock: false });
  assert.equal(differing(off, mixed({})), 0);
  const gap = on.map((x, n) => x.reduce((m, v, i) => Math.max(m, Math.abs(v - off[n][i])), 0));
  console.log(`six steps with the option on against off, the largest differences in h, u, Q, W and η: ${gap.map((g) => g.toExponential(1)).join(', ')}`);
  assert.ok(gap[0] < 1e-6 && gap[1] < 1e-9 && gap[4] < 1e-6, `${gap}`);
});

test('with ekmanBlock off the pass does nothing and the ocean has no Ekman arrays', () => {
  const { ocean, ice } = uniformOcean();
  ocean.setStress(zonalWindOnEdges(mesh, () => 0.1), ice);
  const before = snapshot(ocean);
  ocean.phases.ekman(DT);
  assert.equal(ocean.ekmanFlux, null);
  assert.equal(ocean.ekmanMoved, null);
  assert.equal(ocean.diagnostics().oceanEkmanLimited, undefined);
  assert.equal(differing([ocean.h, ocean.Q, ocean.W], [before.h, before.Q, before.W]), 0);
});

/*
 * A continent from 55°S to 55°N over a sea of a 50 m layer on classes
 * 30, 70, 150, 300 and 3400 m thick, the layers cut to 30 m with the
 * water put into the first class, a 100 m wind depth so that the base of
 * the Ekman block is in the third layer, under a wind with both components.
 */
const SEA = { ...PLAIN_MIXING, everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34, mixedDepth: 50, restoreTime: 1e30, densities: [1026.0, 1026.1, 1026.2, 1026.3, 1026.95], salinities: [34.9, 34.9, 34.9, 34.9, 34.8], bottoms: [80, 150, 300, 600, 12000], momentumDepth: 100, ekmanBlock: true };
function continent() {
  const land = new Uint8Array(C), edgeOcean = new Uint8Array(E);
  for (let i = 0; i < C; i++) land[i] = Math.cos(mesh.lonCell[i]) > 0.4 && Math.abs(mesh.latCell[i]) < 55 * DEG ? 1 : 0;
  for (let e = 0; e < E; e++) edgeOcean[e] = land[mesh.cellsOnEdge[2 * e]] || land[mesh.cellsOnEdge[2 * e + 1]] ? 0 : 1;
  return { land, elevation: Float64Array.from(land, (l) => (l ? 300 : -DEPTH)), edgeOcean };
}
function windOnEdges(eastAt, northAt) {
  const u = new Float64Array(E);
  for (let e = 0; e < E; e++) {
    const x = mesh.xEdge[3 * e], y = mesh.xEdge[3 * e + 1], z = mesh.xEdge[3 * e + 2], r = Math.hypot(x, y), R = Math.hypot(x, y, z);
    const east = [-y / r, x / r, 0], north = [-z * x / (r * R), -z * y / (r * R), r / R], n = [mesh.nEdge[3 * e], mesh.nEdge[3 * e + 1], mesh.nEdge[3 * e + 2]];
    const lat = mesh.latEdge[e], lon = Math.atan2(y, x);
    u[e] = eastAt(lat, lon) * (east[0] * n[0] + east[1] * n[1] + east[2] * n[2]) + northAt(lat, lon) * (north[0] * n[0] + north[1] * n[1] + north[2] * n[2]);
  }
  return u;
}
const coastWind = () => windOnEdges((lat) => -0.1 * Math.cos(3 * lat), (lat, lon) => 0.08 * Math.sin(2 * lat) * Math.cos(lon) + 0.05 * Math.sin(lon));
function coastSea(options, { thin = false, shear = 0 } = {}) {
  const ocean = createOcean(mesh, { ...SEA, geography: continent(), bathymetry: new Float64Array(C).fill(DEPTH), ...options });
  const surfaceT = new Float64Array(C).fill(303), ice = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const { h, Q, W } = ocean;
  const give = (i, from, to, amount) => { const f = amount / h[from * C + i]; for (const store of [Q, W]) store[to * C + i] += store[from * C + i] * f; h[to * C + i] += amount; for (const store of [Q, W]) store[from * C + i] -= store[from * C + i] * f; h[from * C + i] -= amount; };
  for (let i = 0; i < C; i++) {
    if (!ocean.cellOcean[i]) continue;
    give(i, 0, 1, h[i] - 30);
    if (thin) give(i, 2, 3, h[2 * C + i] - 4);
  }
  if (shear) { const wind = zonalWindOnEdges(mesh, () => shear); ocean.u.set(wind); }
  ocean.setStress(coastWind(), ice);
  ocean.readSurface(surfaceT, ice);
  ocean.ekmanTransport(DT);
  return { ocean, surfaceT, ice };
}
/*
 * The upwelling budget b_i = upwellingSource max(0, w_i) dt of every sea
 * cell with a land neighbor from the Ekman transport the sea holds, and the
 * class beneath the base of the block from the cell's column.
 */
function budgets(ocean, source, dt) {
  const b = new Float64Array(C), from = new Int32Array(C).fill(-1);
  for (let i = 0; i < C; i++) {
    if (!ocean.cellOcean[i]) continue;
    let divergence = 0, coast = false;
    for (let m = 0; m < mesh.nEdgesOnCell[i]; m++) {
      const e = mesh.edgesOnCell[mesh.maxEdges * i + m];
      divergence += mesh.edgeSignOnCell[mesh.maxEdges * i + m] * ocean.ekmanFlux[e] * mesh.dvEdge[e];
      if (!ocean.cellOcean[mesh.cellsOnCell[mesh.maxEdges * i + m]]) coast = true;
    }
    if (!coast) continue;
    b[i] = source * Math.max(0, divergence / mesh.areaCell[i]) * dt;
    const depth = Math.max(SEA.momentumDepth, ocean.h[i]);
    let top = ocean.h[i];
    for (let k = 1; k < ocean.layers; k++) {
      const hk = ocean.h[k * C + i];
      if (hk > 5 && top + hk > depth) { from[i] = k; break; }
      top += hk;
    }
  }
  return { b, from };
}
const columnSums = (ocean, store, i) => { let s = 0; for (let k = 0; k < ocean.layers; k++) s += store[k * C + i]; return s; };

/*
 * The same sea twice, once with upwellingSource and once without, each
 * taking one mixed-layer step. A column the budget is not spent in comes
 * out bit for bit the reference's; in the others the mixed layer takes in
 * as much as it does without, and the budget or the demand, the lesser,
 * comes from the class beneath the block's base rather than the first,
 * with heat, salt and water kept.
 */
function against(name, options, { thin = false, shear = 0, source = 20 } = {}) {
  const on = coastSea({ ...options, upwellingSource: source }, { thin, shear }), off = coastSea({ ...options, upwellingSource: 0 }, { thin, shear });
  const { ocean } = on, reference = off.ocean, start = snapshot(ocean), { b, from } = budgets(ocean, source, DT);
  assert.equal(differing([ocean.h, ocean.Q, ocean.W], [reference.h, reference.Q, reference.W]), 0);
  ocean.mixedLayer(DT); reference.mixedLayer(DT);
  let spent = 0, bound = 0, demanded = 0, outside = 0, other = 0;
  for (let i = 0; i < C; i++) {
    if (!ocean.cellOcean[i]) continue;
    const demand = reference.h[i] - start.h[i], expected = Math.min(b[i], demand);
    if (!(b[i] > 0)) {
      outside++;
      for (let k = 0; k < ocean.layers; k++) for (const store of ['h', 'Q', 'W']) assert.ok(Object.is(ocean[store][k * C + i], reference[store][k * C + i]), `${name}: cell ${i} has no budget and layer ${k}'s ${store} differs`);
      assert.equal(ocean.upwelled[i], 0);
      continue;
    }
    assert.ok(Math.abs(ocean.upwelled[i] - expected) <= 1e-12 * Math.max(b[i], demand) + 1e-15, `${name}: cell ${i} drew ${ocean.upwelled[i]} m from the source, expected min(${b[i]}, ${demand})`);
    assert.ok(Math.abs(ocean.h[i] - reference.h[i]) <= 1e-12, `${name}: cell ${i} layer ${ocean.h[i]} against ${reference.h[i]}`);
    if (expected > 0) {
      const s = from[i], first = [1, 2, 3, 4, 5].find((k) => start.h[k * C + i] > 5);
      assert.ok(s > 0);
      if (s !== first) other++;
      assert.ok(ocean.upwelled[i] <= b[i] * (1 + 1e-12) + 1e-15);
      for (let k = 1; k < ocean.layers; k++) {
        const change = ocean.h[k * C + i] - reference.h[k * C + i], expectedChange = (k === first ? 1 : 0) * ocean.upwelled[i] - (k === s ? 1 : 0) * ocean.upwelled[i];
        if (s === first) assert.ok(Math.abs(change) <= 1e-12, `${name}: cell ${i} class ${k} moved ${change} m with the source the first class`);
        else assert.ok(Math.abs(change - expectedChange) <= 1e-12, `${name}: cell ${i} class ${k} differs from the reference by ${change} m, expected ${expectedChange}`);
      }
      spent++; if (b[i] < demand) bound++; else demanded++;
    }
    for (const store of ['h', 'Q', 'W']) {
      const was = columnSums(ocean, start[store], i), is = columnSums(ocean, ocean[store], i);
      assert.ok(Math.abs(is - was) <= 1e-12 * Math.abs(was) * 10, `${name}: cell ${i} ${store} ${was} → ${is}`);
    }
  }
  console.log(`${name}: ${spent} coastal columns drew ${Array.from(ocean.upwelled).filter((x) => x > 0).length} budgets, ${bound} held to b_i and ${demanded} to the mixed layer's demand, ${other} from a class other than the first thick one; ${outside} sea columns without a budget are the reference's bit for bit`);
  return { spent, bound, demanded, other, outside, ocean, reference, start, b, from };
}

test('under upwellingSource the mixed layer\'s refill takes the Ekman divergence of a coastal cell from the class that holds the block\'s base and the rest as it does without, other cells untouched', () => {
  const result = against('refill', { refillRate: 3 / 86400, minimumThickness: 50, stirring: 0 });
  assert.ok(result.bound >= 2 && result.demanded >= 2 && result.other >= 4, JSON.stringify([result.bound, result.demanded, result.other]));
  const thin = against('refill, a thin class at the base', { refillRate: 3 / 86400, minimumThickness: 50, stirring: 0 }, { thin: true });
  assert.ok(thin.bound >= 2 && thin.demanded >= 2 && thin.other >= 4);
  let sources = new Set();
  for (let i = 0; i < C; i++) if (thin.b[i] > 0) sources.add(thin.from[i]);
  assert.ok(sources.has(3) && !sources.has(2), `${[...sources]}`);
});

test('under upwellingSource the mixed layer\'s stirring entrainment takes the budget from the same class and the rest from the first thick class', () => {
  const result = against('stirring', { refillRate: 0, minimumThickness: 1, stirring: 100 }, { source: 60 });
  assert.ok(result.bound >= 2 && result.demanded >= 2 && result.other >= 4, JSON.stringify([result.bound, result.demanded, result.other]));
});

test('under upwellingSource the mixed layer\'s shear entrainment takes the budget from the same class and the rest from the first thick class', () => {
  const result = against('shear', { refillRate: 0, minimumThickness: 1, stirring: 0, shearEntrainment: 0.65, carryMomentum: true }, { shear: 3.5, source: 20 });
  assert.ok(result.bound >= 2 && result.other >= 4, JSON.stringify([result.bound, result.demanded, result.other]));
});

test('the shear entrainment, the stirring and the refill share one budget, and a cell without a land neighbor or with a convergence draws nothing from the source', () => {
  const options = { refillRate: 3 / 86400, minimumThickness: 50, stirring: 100, shearEntrainment: 0.65, carryMomentum: true };
  const { ocean } = coastSea({ ...options, upwellingSource: 20 }, { shear: 1.5 }), { b } = budgets(ocean, 20, DT);
  ocean.mixedLayer(DT);
  let spent = 0, whole = 0;
  for (let i = 0; i < C; i++) {
    assert.ok(ocean.upwelled[i] <= b[i] * (1 + 1e-12) + 1e-15, `cell ${i} drew ${ocean.upwelled[i]} of ${b[i]} m`);
    if (b[i] > 0 && ocean.upwelled[i] > 0) spent++;
    if (b[i] > 0 && Math.abs(ocean.upwelled[i] - b[i]) <= 1e-12 * b[i]) whole++;
    if (!(b[i] > 0)) assert.equal(ocean.upwelled[i], 0);
  }
  console.log(`all three processes under one budget: ${spent} cells drew from the source, ${whole} the whole of it`);
  assert.ok(spent >= 8 && whole >= 4);
});

test('the intake ledger splits the mixed layer\'s intake by process and class without changing the ocean, its upwelled share is the source\'s draw, and a mask records only its cells', () => {
  const all = { refillRate: 3 / 86400, minimumThickness: 50, stirring: 100, shearEntrainment: 0.65, carryMomentum: true, upwellingSource: 20 };
  const scenarios = [
    { options: all, fixture: { shear: 1.5 } },
    { options: { ...all, refillRate: 0, minimumThickness: 1, stirring: 0 }, fixture: { shear: 3.5 } },
    { options: all, fixture: { shear: 1.5 }, cool: 29 },
  ];
  const total = { convective: 0, shear: 0, stirring: 0, refill: 0, upwelled: 0 };
  let cold = 0, outside = 0, first = null;
  for (const { options, fixture, cool } of scenarios) {
    const make = (extra) => {
      const { ocean, surfaceT, ice } = coastSea({ ...options, ...extra }, fixture);
      if (cool) { surfaceT.fill(303 - cool); ocean.readSurface(surfaceT, ice); }
      return ocean;
    };
    const plain = make({}), ocean = make({ intakeLedger: true });
    const start = snapshot(ocean), L = ocean.layers, { intake } = ocean;
    assert.equal(plain.intake, null);
    assert.deepEqual(Array.from(intake.processes), ['convective', 'shear', 'stirring', 'refill', 'upwelled']);
    assert.equal(intake.cells.length, Array.from(ocean.cellOcean).filter((x) => x).length);
    plain.mixedLayer(DT); ocean.mixedLayer(DT);
    assert.equal(differing([ocean.h, ocean.Q, ocean.W], [plain.h, plain.Q, plain.W]), 0);
    for (let slot = 0; slot < intake.cells.length; slot++) {
      const i = intake.cells[slot], by = {};
      let sum = 0;
      for (let p = 0; p < intake.processes.length; p++) {
        let water = 0;
        for (let k = 0; k < L; k++) {
          const n = (p * L + k) * intake.cells.length + slot;
          assert.ok(intake.water[n] >= 0);
          if (intake.water[n] > 0) assert.ok(intake.degrees[n] / intake.water[n] > 250 && intake.degrees[n] / intake.water[n] < 320, `class ${k} at ${intake.degrees[n] / intake.water[n]}`);
          if (k === 0) assert.equal(intake.water[n], 0, 'the mixed layer is not its own source');
          water += intake.water[n];
        }
        by[intake.processes[p]] = water; total[intake.processes[p]] += water; sum += water;
      }
      assert.ok(Math.abs(by.upwelled - (ocean.upwelled[i] ?? 0)) <= 1e-12 * Math.max(1, ocean.upwelled[i]), `cell ${i}: ${by.upwelled} against ${ocean.upwelled[i]}`);
      assert.ok(by.shear <= ocean.sheared[i] + 1e-12 && ocean.sheared[i] - by.shear <= by.upwelled + 1e-12);
      assert.ok(ocean.refilled.water[i] - by.refill <= by.upwelled + 1e-12);
      assert.ok(sum >= ocean.h[i] - start.h[i] - 1e-12, `cell ${i} gained ${ocean.h[i] - start.h[i]} m from ${sum} m`);
      if (by.upwelled > 0) cold++; else if (ocean.upwelled[i] === 0) outside++;
    }
    first ??= { make, intake, L };
  }
  console.log(`intake in m over ${scenarios.length} scenarios: ${JSON.stringify(Object.fromEntries(Object.entries(total).map(([name, x]) => [name, +x.toFixed(2)])))}`);
  assert.ok(cold >= 8 && outside > cold);
  for (const [name, x] of Object.entries(total)) assert.ok(x > 0, `no ${name} intake`);
  const { make, intake, L } = first, mask = new Uint8Array(C);
  for (let slot = 0; slot < intake.cells.length; slot += 5) mask[intake.cells[slot]] = 1;
  const masked = make({ intakeLedger: mask }), full = make({ intakeLedger: true });
  masked.mixedLayer(DT); full.mixedLayer(DT);
  assert.equal(masked.intake.cells.length, mask.reduce((a, b) => a + b, 0));
  for (let slot = 0; slot < masked.intake.cells.length; slot++) {
    const i = masked.intake.cells[slot], at = full.intake.slot[i];
    assert.ok(mask[i] && at >= 0);
    for (let p = 0; p < intake.processes.length; p++) for (let k = 0; k < L; k++) {
      assert.equal(masked.intake.water[(p * L + k) * masked.intake.cells.length + slot], full.intake.water[(p * L + k) * full.intake.cells.length + at]);
    }
  }
  for (let i = 0; i < C; i++) if (!mask[i]) assert.equal(masked.intake.slot[i], -1);
});

test('upwellingSource needs ekmanBlock and no less than zero, and at zero the ocean is the same bit for bit with no source array', () => {
  assert.throws(() => createOcean(mesh, { ...SEA, ekmanBlock: false, upwellingSource: 1 }), /upwellingSource/);
  assert.throws(() => createOcean(mesh, { ...SEA, upwellingSource: -1 }), /upwellingSource/);
  const run = (options) => {
    const { ocean, surfaceT, ice } = coastSea(options);
    const flux = new Float64Array(C), wind = coastWind();
    for (let n = 0; n < 6; n++) ocean.advance(surfaceT, ice, flux, wind, DT, undefined, (n + 1) * DT);
    return { ocean, state: [ocean.h, ocean.u, ocean.Q, ocean.W, ocean.eta].map((x) => Float64Array.from(x)) };
  };
  const base = { refillRate: 3 / 86400, minimumThickness: 50, stirring: 100, shearEntrainment: 0.65, carryMomentum: true };
  const absent = run(base), zero = run({ ...base, upwellingSource: 0 }), on = run({ ...base, upwellingSource: 20 });
  assert.equal(differing(absent.state, zero.state), 0);
  assert.equal(zero.ocean.upwelled, null);
  const gap = on.state.map((x, n) => x.reduce((m, v, i) => Math.max(m, Math.abs(v - absent.state[n][i])), 0));
  console.log(`six steps with upwellingSource 20 against none, the largest differences in h, u, Q, W and η: ${gap.map((g) => g.toExponential(1)).join(', ')}`);
  assert.ok(gap[0] > 1e-3 && Array.from(on.ocean.upwelled).some((x) => x > 0));
});
