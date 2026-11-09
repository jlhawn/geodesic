import test from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { F32_CELL, F32_RMS, S9_TIERS, buildSpread, columnLargest, distributionCheck, hemispheres, roundDown, roundUp, spreadStatistics } from './helpers/distribution.mjs';
import { bandLayers, groupChecks, groupLayout, spreadGroups } from './helpers/spreadGroups.mjs';

const mesh = buildMesh(new Grid(16)), n = mesh.nCells, sides = hemispheres(mesh), REFERENCE = { max: 20, rms: 5 };

function generator(seed) {
  const uniform = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; };
  return { uniform, normal: () => Math.sqrt(-2 * Math.log(uniform() + 1e-12)) * Math.cos(2 * Math.PI * uniform()) };
}

function twins(seed, count, size) {
  const { uniform, normal } = generator(seed), knifeEdge = Array.from({ length: 40 }, () => Math.floor(uniform() * size));
  return Array.from({ length: count }, () => {
    const gaps = Float64Array.from({ length: size }, () => 2e-4 * Math.exp(1.2 * normal()) * (uniform() < 0.5 ? 1 : -1)), bias = 1e-5 * normal();
    for (let x = 0; x < size; x++) gaps[x] += bias;
    for (let flip = 0; flip < 3; flip++) gaps[knifeEdge[Math.floor(uniform() * 40)]] += (uniform() < 0.5 ? 0.5 : -0.5) * (0.05 + 0.95 * uniform());
    return gaps;
  });
}

const ensemble = twins(12345, 30, n), build = ensemble.slice(0, 20), fresh = ensemble.slice(20), spread = buildSpread(build, sides, REFERENCE), M = spread.twins.largest;
const failures = (gaps, bounds = spread) => distributionCheck(gaps, sides, bounds, REFERENCE, { require: S9_TIERS }).failed;
const mutant = (change) => { const gaps = Float64Array.from(fresh[0]); change(gaps); return failures(gaps); };
const caughtBy = (failed, pattern) => failed.some((text) => pattern.test(text));

test('columnLargest keeps a column\'s signed largest gap over the layers it is given, and a NaN or an infinity in any of them', () => {
  assert.deepEqual([...columnLargest(Float64Array.from([1, -3, 2]), 1, [0, 1, 2])], [-3]);
  assert.deepEqual([...columnLargest(Float64Array.from([1, 5, -3, 2, 9, -40]), 2, [0, 1])], [-3, 5]);
  for (const bad of [NaN, Infinity, -Infinity]) {
    for (let layer = 0; layer < 3; layer++) {
      const gaps = Float64Array.from([1, -3, 2]);
      gaps[layer] = bad;
      const [worst] = columnLargest(gaps, 1, [0, 1, 2]);
      assert.ok(Number.isNaN(bad) ? Number.isNaN(worst) : worst === bad, `${bad} in layer ${layer} of 3 gives ${worst}`);
    }
  }
});

test('roundUp and roundDown give the nearest three-figure value on their side, at exact decades too', () => {
  for (const [x, up] of [[0.001, 0.001], [1000, 1000], [0.00101, 0.00101], [0.0010001, 0.00101], [9.991, 10], [1.0000001e-5, 1.01e-5]]) assert.equal(roundUp(x), up, `roundUp(${x})`);
  for (const [x, down] of [[0.001, 0.001], [9.99999e-6, 9.99e-6], [0.0019999, 0.00199], [1000.5, 1000]]) assert.equal(roundDown(x), down, `roundDown(${x})`);
  const { uniform } = generator(2718);
  for (let trial = 0; trial < 20000; trial++) {
    const x = Math.exp((uniform() - 0.5) * 120), up = roundUp(x), down = roundDown(x);
    assert.ok(up >= x && +up.toPrecision(3) === up && roundDown(up * (1 - 1e-9)) < x, `roundUp(${x}) = ${up}`);
    assert.ok(down <= x && +down.toPrecision(3) === down && roundUp(down * (1 + 1e-9)) > x, `roundDown(${x}) = ${down}`);
  }
});

test('the ten equatorial cells of the N=16 mesh count in both hemispheres, and a plain array of sides or a latitude array agrees', () => {
  const equator = sides.filter((side) => side === 0).length, ones = new Float64Array(n).fill(1), { north, south } = spreadStatistics(ones, sides).hemisphere;
  assert.equal(equator, 10);
  assert.equal(north.cells + south.cells, n + equator);
  const plain = spreadStatistics(ones, Array.from(sides)).hemisphere;
  assert.deepEqual([plain.north.cells, plain.south.cells], [north.cells, south.cells]);
  const latitudes = spreadStatistics(ones, mesh.latCell).hemisphere;
  assert.equal(latitudes.north.cells + latitudes.south.cells, n);
  const marked = spreadStatistics([1, 1, 1], [0.5, 0, -0.5]).hemisphere;
  assert.deepEqual([marked.north.cells, marked.south.cells], [2, 2]);
});

test('buildSpread bounds the fresh twins of its own ensemble and catches what none of them does', () => {
  console.log(`twenty twins of ${n} gaps: M ${M.toExponential(2)}, ladder ${JSON.stringify(spread.counts)}, hemisphere ${JSON.stringify(spread.hemisphere)}, bulk ${JSON.stringify(spread.bulk)}`);
  assert.ok(spread.counts[0][0] <= M && spread.twins.over[0] >= 1, `the top rung ${spread.counts[0][0]} against M ${M}`);
  [...build, ...fresh].forEach((gaps, t) => assert.deepEqual(failures(gaps), [], `twin ${t}`));
  const bias = mutant((gaps) => gaps.forEach((_, x) => { gaps[x] += 3e-4; }));
  assert.ok(caughtBy(bias, /trimmed signed mean/), `a uniform bias of 3e-4: ${bias.join('; ')}`);
  const columns = mutant((gaps) => { for (let x = 0; x < 17 * 150; x += 150) gaps[x] = 0.97 * M; });
  assert.ok(caughtBy(columns, /at or over/), `17 cells at 0.97 M: ${columns.join('; ')}`);
  assert.ok(caughtBy(mutant((gaps) => { gaps[5] = 2.1 * M; }), /^largest/), 'one cell at 2.1 M');
  assert.deepEqual(mutant((gaps) => { gaps[5] = 1.9 * M; }), [], 'one cell at 1.9 M');
  const opposite = mutant((gaps) => gaps.forEach((_, x) => { gaps[x] += Math.sign(sides[x]) * 3e-4; }));
  assert.ok(caughtBy(opposite, /northern trimmed signed/) && caughtBy(opposite, /southern trimmed signed/), `hemispheres apart by ±3e-4: ${opposite.join('; ')}`);
  for (const bad of [NaN, Infinity]) {
    assert.ok(caughtBy(mutant((gaps) => { gaps[100] = bad; }), /non-finite gap at 100/), `${bad} at one cell`);
    const twin = Float64Array.from(build[0]);
    twin[3] = bad;
    assert.throws(() => buildSpread([...build, twin], sides, REFERENCE), /not finite/, `${bad} in a twin`);
  }
  assert.throws(() => buildSpread(build, sides, { max: NaN, rms: 1 }), /reference/);
});

test('twins that are all zero (M = 0) bound a group at the f32 floors', () => {
  const zeros = Array.from({ length: 20 }, () => new Float64Array(n)), floor = buildSpread(zeros, sides, REFERENCE), mean = roundUp(F32_RMS * REFERENCE.rms);
  assert.deepEqual(floor.counts, []);
  assert.equal(floor.largest, roundUp(F32_CELL * REFERENCE.max));
  for (const side of ['north', 'south']) assert.deepEqual([floor.hemisphere[side], floor.bulk[side]], [mean, { trimmed: mean, signed: mean }]);
  assert.deepEqual(failures(new Float64Array(n), floor), []);
  const cell = new Float64Array(n);
  cell[7] = 2 * F32_CELL * REFERENCE.max;
  assert.ok(caughtBy(failures(cell, floor), /^largest/));
  const zero = { max: 0, rms: 0 }, exact = buildSpread(zeros, sides, zero);
  assert.deepEqual(distributionCheck(new Float64Array(n), sides, exact, zero, { require: S9_TIERS }).failed, []);
  cell[7] = 1e-30;
  assert.ok(caughtBy(distributionCheck(cell, sides, exact, zero, { require: S9_TIERS }).failed, /^largest/));
});

test('an empty group, a where of the wrong length, a non-finite bound, a missing tier and a hemisphere the group lacks fail or throw', () => {
  assert.ok(caughtBy(distributionCheck([], new Int8Array(0), spread, REFERENCE).failed, /an empty group/));
  assert.throws(() => buildSpread([new Float64Array(0)], new Int8Array(0), REFERENCE), /an empty group/);
  assert.throws(() => buildSpread([], sides, REFERENCE), /no twins/);
  assert.throws(() => buildSpread([build[0], build[1].subarray(1)], sides, REFERENCE), /different numbers/);
  assert.throws(() => failures(fresh[0].subarray(1)), /entries for/);
  assert.throws(() => distributionCheck(fresh[0], sides, spread, REFERENCE, { cells: [1, 2] }), /cell ids/);
  assert.throws(() => spreadStatistics([1, 1], [0.5, NaN]), /where holds NaN/);
  const broken = {
    largest: { ...spread, largest: Infinity },
    rung: { ...spread, counts: [[Infinity, 3], ...spread.counts.slice(1)] },
    hemisphere: { ...spread, hemisphere: { ...spread.hemisphere, north: Infinity } },
    bulk: { ...spread, bulk: { ...spread.bulk, south: { ...spread.bulk.south, signed: Infinity } } },
  };
  for (const [name, bounds] of Object.entries(broken)) assert.ok(caughtBy(failures(fresh[0], bounds), /bound is Infinity/), `an infinite ${name} bound`);
  for (const tier of S9_TIERS) {
    const without = { ...spread };
    delete without[tier];
    assert.ok(caughtBy(failures(fresh[0], without), new RegExp(`the spread has no ${tier} tier`)), `a spread without ${tier}`);
  }
  assert.deepEqual(distributionCheck(fresh[0], sides, { largest: spread.largest, counts: spread.counts, hemisphere: spread.hemisphere }, REFERENCE).failed, []);
  assert.ok(caughtBy(failures(fresh[0], { ...spread, bulk: { north: spread.bulk.north } }), /no southern bulk bounds/));
  assert.ok(caughtBy(distributionCheck(fresh[0], new Int8Array(n).fill(1), spread, REFERENCE).failed, /southern bounds for a group with no southern cells/));
});

test('groupLayout puts each cell in one surface and each edge between land and sea on the coast, and bandLayers splits at 700 and 100 hPa', () => {
  const land = Uint8Array.from(mesh.latCell, (lat) => (lat > 0.5 ? 1 : 0)), cells = groupLayout(mesh, { land }), edges = groupLayout(mesh, { land, edges: true });
  assert.deepEqual(Object.keys(cells.members).sort(), ['land', 'sea']);
  assert.equal(cells.members.land.length + cells.members.sea.length, n);
  assert.ok(cells.members.land.every((i) => land[i]) && cells.members.sea.every((i) => !land[i]));
  assert.deepEqual(Object.keys(edges.members).sort(), ['coast', 'land', 'sea']);
  assert.equal(Object.values(edges.members).reduce((sum, list) => sum + list.length, 0), mesh.nEdges);
  assert.ok(edges.members.coast.every((e) => land[mesh.cellsOnEdge[2 * e]] !== land[mesh.cellsOnEdge[2 * e + 1]]));
  assert.deepEqual(bandLayers([0.9, 0.7, 0.5, 0.1, 0.05], 5), { low: [0, 1], mid: [2, 3], high: [4] });
  assert.deepEqual(bandLayers(null, 1), { all: [0] });
});

test('groupChecks fails the band whose layers carry a bias and no other, and a group without bounds or a bulk tier', () => {
  const layers = 3, bands = bandLayers([0.85, 0.3, 0.05], layers), layout = groupLayout(mesh), { uniform } = generator(99);
  const cpu = Float64Array.from({ length: layers * n }, () => 20 * (2 * uniform() - 1)), engine = (gaps) => Float64Array.from(cpu, (v, x) => v + gaps[x]);
  const field = twins(4242, 30, layers * n), runs = field.map((gaps) => spreadGroups(cpu, engine(gaps), layout, bands)), names = Object.keys(runs[0]).sort();
  assert.deepEqual(names, ['sea high', 'sea low', 'sea mid']);
  const spreads = Object.fromEntries(names.map((name) => [name, buildSpread(runs.slice(0, 20).map((groups) => groups[name].gaps), runs[0][name].where, runs[0][name].reference)]));
  runs.forEach((groups, t) => assert.deepEqual(groupChecks(groups, spreads).failed, [], `twin ${t}`));
  const biased = Float64Array.from(field[20]);
  for (let x = n; x < 2 * n; x++) biased[x] += 3e-4;
  const failed = groupChecks(spreadGroups(cpu, engine(biased), layout, bands), spreads).failed;
  assert.ok(failed.length > 0 && failed.every((text) => text.startsWith('sea mid:')), failed.join('; '));
  const partial = { ...spreads };
  delete partial['sea high'];
  assert.ok(groupChecks(runs[20], partial).failed.includes('sea high: no bounds for this group'));
  assert.ok(groupChecks(runs[20], { ...spreads, 'land low': spreads['sea low'] }).failed.includes('land low: bounds for a group the run does not have'));
  const noBulk = { ...spreads['sea low'] };
  delete noBulk.bulk;
  assert.ok(groupChecks(runs[20], { ...spreads, 'sea low': noBulk }).failed.includes('sea low: the spread has no bulk tier'));
});
