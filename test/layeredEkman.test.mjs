import { test } from 'node:test';
import assert from 'node:assert/strict';
import { cellVector } from '../js/dynamics/operators.module.js';
import { createOcean, EPS, EDDY_SLACK } from '../js/ocean/layered.module.js';
import { mesh, C, E, DEG, RHO, uniformOcean, zonalWindOnEdges, deepen } from './helpers/layered.mjs';

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
