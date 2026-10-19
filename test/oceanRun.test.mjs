import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
import { HYBRID_LAYERS, MODE_DEPTH as MODEL_MODE_DEPTH } from '../js/ocean/column.module.js';
import { readFrames } from '../explain/frames.module.js';
import { heatTransport } from '../explain/sunlitRun.module.js';
import { storage, waterTransport, recordedYear as wetYear } from '../explain/wetRun.module.js';
import {
  DATA_URL, PART7_TRANSPORT, PART6_TRANSPORT, GROUPS, MODE_DEPTH, calendarDate, peakLatitude, recordedYear, globalMean, integrate, latentNet,
  heatBooks, cellHeat, transportSplit, layerGroups, waterCycle, yearlyMeans, zonalMeans, layerSection, layerOverturning, groupOf, modeWeight, modeName,
} from '../explain/oceanRun.module.js';

const open = (name) => { const file = new URL(`../explain/data/${name}`, import.meta.url), bytes = readFileSync(file); return readFrames(file.href, bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength)); };
const run = open('ocean_N16.bin'), part7 = open('water_N16.bin');
const { header } = run, DAY = 86400, YEAR = 365, D = Math.PI / 180, L = header.layers, F = header.fixedLayers;
const cells = [...new Grid(header.N)].map((c) => ({ lat: Math.asin(c.centerVertex.z) / D, area: c.area }));
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const date = (day) => { const c = calendarDate(day); return `${c.month} ${c.date}`; };
const zonal = (values) => { const sum = new Float64Array(36), area = new Float64Array(36); cells.forEach((c, i) => { const b = Math.min(35, Math.floor((c.lat + 90) / 5)); sum[b] += values[i] * c.area; area[b] += c.area; }); return Array.from(sum, (v, b) => v / area[b]); };
const band = (values, lo, hi) => { let sum = 0, area = 0; cells.forEach((c, i) => { if (c.lat >= lo && c.lat < hi) { sum += values[i] * c.area; area += c.area; } }); return sum / area; };
const field = (h, name) => h.fields.findIndex((f) => f.name === name);
const sinAreas = (radius) => Array.from({ length: 36 }, (_, b) => 2 * Math.PI * radius * radius * (Math.sin((-85 + 5 * b) * D) - Math.sin((-90 + 5 * b) * D)));
const area = header.zonalBudget.binArea, sphere = 4 * Math.PI * header.constants.radius ** 2;
const books = heatBooks(header), split = transportSplit(header), water = waterCycle(header), year = recordedYear(header), part7Year = wetYear(part7.header);
const lat = (b) => -87.5 + 5 * b, bin = (l) => (l + 87.5) / 5;
const month = (name, part = 'zonalBudget') => header[part].months.find((m) => m.name === name);
const capMean = (values, from, to) => { let sum = 0, a = 0; for (let b = from; b < to; b++) { sum += values[b] * area[b]; a += area[b]; } return sum / a; };
const lastDay = (name, over = 0) => { let n = -1; header.daily[name].forEach((v, d) => { if (v > over) n = d; }); return n + 1; };

test('the run is Part 7’s wet planet with the hybrid sea let go: ten years of spin-up and a recorded year from the March equinox', () => {
  assert.equal(new URL(DATA_URL).pathname.endsWith('/explain/data/ocean_N16.bin'), true);
  assert.deepEqual([header.case, header.engine, header.N, header.K, header.layers, header.dt, header.grid, header.spinupDays, header.recordDays], ['ocean', 'gpu', 16, 36, 58, 450, 'bl36', 3650, 365]);
  assert.equal(header.commit, '40fdb03');
  assert.equal(header.date.slice(0, 10), '2026-10-19');
  assert.equal(header.options.ocean.everySteps, 4);
  assert.deepEqual(header.ice.options, part7.header.ice.options);
  assert.deepEqual(header.floor, { kind: 'flat', depth: 4000, minimum: 4000, maximum: 4000, land: 'none: every cell and edge is sea' });
  assert.equal(header.recordStart, 3650 * DAY);
  assert.deepEqual(header.sequences.map((s) => [s.name, s.times.length]), [['year', 74], ['spinup', 61]]);
  assert.equal(header.daily.surfaceTemperature.length, header.spinupDays + header.recordDays);
  for (const counts of [header.nonFinite, header.clipped]) assert.ok(Object.values(counts).every((n) => n === 0), 'every frame value finite and in range');
  assert.ok(header.resumedFrom === null ? header.resumeCheck === null && header.segments.length === 1 : !header.resumeCheck.forced && header.resumeCheck.sameEngine && Object.values(header.resumeCheck.stores).every((v) => v === 0), 'one leg, or resumed legs that are one run');
  assert.equal(Math.max(...header.daily.clampedEdges), 0);
});

test('the layers are the model’s hybrid table: twelve fixed-depth surface layers, labels and minimum thicknesses as the code has them', () => {
  assert.deepEqual(header.classes, Array.from(HYBRID_LAYERS.densities));
  assert.deepEqual(header.minimumThickness, Array.from(HYBRID_LAYERS.minimum));
  assert.equal(F, HYBRID_LAYERS.fixed); assert.equal(F, 12); assert.equal(L, 58);
  assert.equal(MODE_DEPTH, MODEL_MODE_DEPTH); assert.equal(MODE_DEPTH, header.oceanOptions.modeDepth);
  assert.equal(header.oceanOptions.minimumLayer, 0.5); assert.equal(header.oceanOptions.generator, 'hybrid');
  const depths = HYBRID_LAYERS.depths;
  close(depths[12], 56.9, 0.05, 'the twelve surface layers'); close(depths[21], 144, 0.5, 'the warm-pool layers end'); close(depths[44], 374, 0.5, 'the thermocline layers end'); close(depths[58], 3019, 0.5, 'the minima');
  assert.deepEqual([header.classes[0], header.classes[11], header.classes[12], header.classes[20], header.classes[21], header.classes[43], header.classes[44], header.classes[57]], [1016.5, 1019.25, 1019.5, 1021.875, 1022, 1026.6, 1026.65, 1027.5]);
  assert.deepEqual([GROUPS.length, groupOf(0, 1016.5, F), groupOf(11, 1019.25, F), groupOf(12, 1019.5, F), groupOf(38, 1026.1, F), groupOf(49, 1026.8, F), groupOf(57, 1027.5, F)], [4, 0, 0, 1, 2, 3, 3]);
  const { groups } = layerGroups(header);
  assert.deepEqual(groups.map((g) => g.name), GROUPS);
  assert.deepEqual([modeWeight(3, 3), modeWeight(5.5, 3), modeWeight(9, 3), modeName(0), modeName(0.4), modeName(1)], [0, 0.5, 1, 'fixed-depth', 'between', 'isopycnal']);
  const section = layerSection(header);
  section.columns.forEach((col, b) => {
    for (let k = 0; k < F; k++) close(col[k].thickness, header.minimumThickness[k], 0.2, `layer ${k} at its minimum at ${lat(b)}`);
    assert.ok(col.slice(0, F).every((c) => c.mode === 'fixed-depth'));
  });
  for (const l of [-2.5, 2.5]) {
    const col = section.columns[bin(l)], held = col.filter((c) => c.weight < 0.02), first = col.find((c) => c.weight > 0.5);
    assert.equal(first.layer, 35); assert.equal(first.label, 1025.5); close(first.top, 284.3, 0.5, `the first isopycnal layer starts at ${l}`);
    assert.ok([...Array(34).keys()].every((k) => col[k].weight < 0.02), `layers to 33 held at ${l}`); assert.ok(col[34].weight < 0.3, 'and the next within 1.5 m of its minimum');
    close(held.length, 43, 0, 'layers held at their minimum, the nine sheets on the floor with them');
    assert.equal(col.filter((c) => c.weight > 0.5).length, 14);
  }
  const first = (b) => section.columns[b].find((c) => c.layer >= F && c.weight > 0.5);
  close(first(bin(17.5)).top, 233.7, 0.5, 'isopycnal layers begin at 234 m at 17.5° N'); close(first(bin(62.5)).top, 364, 0.5, 'and at 364 m at 62.5° N');
});

test('Part 7’s calendar carries over: time 0 is noon on March 20, and the months are the header’s', () => {
  assert.deepEqual(calendarDate(0), { month: 'March', date: 20, hour: 12 });
  assert.equal(date(run.times('spinup')[0] / DAY), 'March 20');
  assert.equal(date(run.times('spinup').at(-1) / DAY), 'May 19');
  for (const m of header.zonalBudget.months) {
    assert.equal(calendarDate(m.ranges[0][0]).month, m.name, `${m.name} starts its range`);
    const end = m.ranges.at(-1)[1];
    assert.notEqual(calendarDate(end).month, m.name, `${m.name} ends its range`);
    assert.equal(calendarDate(end - 1).month, m.name, `${m.name}’s last day`);
  }
  assert.deepEqual(header.heatBooks.boundaries.map((b) => b.time / DAY), [0, 12, 42, 73, 103, 134, 165, 195, 226, 256, 287, 318, 346, 365]);
});

test('the bins tile the sphere, and the books’ bands add up to their globe', () => {
  close(area.reduce((s, a) => s + a, 0) / sphere, 1, 1e-9, 'binArea');
  const own = new Float64Array(36);
  cells.forEach((c) => { own[Math.min(35, Math.floor((c.lat + 90) / 5))] += c.area; });
  const unit = own.reduce((s, a) => s + a, 0);
  own.forEach((a, b) => close(a / unit, area[b] / sphere, 1e-9, `bin ${b} holds the cells centered in it`));
  const first = header.heatBooks.boundaries[0], last = header.heatBooks.boundaries.at(-1);
  assert.equal(books.seconds, YEAR * DAY);
  for (const name of Object.keys(books.stores)) {
    books.stores[name].forEach((v, b) => close(v * books.seconds, last.bands[name][b] - first.bands[name][b], 1e-6 * Math.max(1, Math.abs(last.bands[name][b])), `${name} ${b}`));
    close(books.global[name], (last.global[name] - first.global[name]) / books.seconds, 1e-6, `${name} global`);
  }
  for (const name of ['ocean', 'ice', 'airEnthalpy', 'airKinetic', 'latent']) close(books.global[name], header.books.record.storageChange[name], 1e-6, name);
  for (const name of ['toaNet', 'surfaceNet', 'snowFusion']) close(books.global[name], header.books.record[name], 1e-5, name);
  for (const name of ['residual', 'surfaceResidual', 'atmosphereResidual']) close(books.global[name], header.books.record[name], 1e-5, name);
  close(books.global.surfaceResidual + books.global.atmosphereResidual, books.global.residual, 1e-12, 'the two parts');
  const weighted = header.heatBooks.periods.reduce((sum, p) => sum.map((v, b) => v + p.storageChange.ocean[b] * (p.to - p.from) / books.seconds), new Array(36).fill(0));
  weighted.forEach((v, b) => close(v, books.stores.ocean[b], 1e-4, `the months add up to the year at ${lat(b)}`));
  books.latent.forEach((v, b) => close(v, header.constants.latentHeat * (books.stores.vapor[b] + books.stores.condensate[b]), 1e-9, `latent ${b}`));
  books.stores.latent.forEach((v, b) => close(v, header.constants.latentHeat * books.stores.vapor[b], 1e-9, `the books’ latent store is L times the vapor at ${b}`));
});

test('each transport is zero at both poles once its own imbalance is spread, and the parts add up', () => {
  const { total, sea, air, latent, dry } = split;
  for (const t of [total, sea, air, latent, dry]) { assert.equal(t.transport[0], 0); close(t.transport.at(-1), 0, 1e-12, 'north pole'); }
  air.transport.forEach((v, e) => close(v, total.transport[e] - sea.transport[e], 1e-12, `air at ${split.edges[e]}`));
  dry.transport.forEach((v, e) => close(v + latent.transport[e], air.transport[e], 1e-12, `dry at ${split.edges[e]}`));
  close(total.imbalance, -header.books.record.residual, 1e-5, 'the total spreads the books’ residual');
  close(sea.imbalance, -header.books.record.surfaceResidual, 1e-5, 'the sea its surface side');
  close(air.imbalance, -header.books.record.atmosphereResidual, 1e-5, 'the air its own');
  assert.ok(Math.abs(water.waterImbalance) < 0.001, `water is kept: ${water.waterImbalance} mm a day`);
  for (const [name, t] of [['total', total], ['sea', sea], ['latent', latent]]) {
    const unspread = t.net.reduce((sum, v, b) => sum + v * area[b], 0) / 1e15;
    close(unspread, t.imbalance * sphere / 1e15, 1e-9, `${name}: without the spread the north pole would be off by the imbalance`);
  }
  for (const l of [-63, -20, 0, 7, 42]) { const at = split.at(l); close(at.air - at.latent, at.dry, 1e-12, `at ${l}`); close(at.total - at.sea, at.air, 1e-12, `at ${l}`); }
  const corrected = integrate(books.toaNet, area).transport.map((v, e) => v + books.corrections.total[e]);
  corrected.forEach((v, e) => close(v, total.transport[e], 1e-12, `storage correction at ${split.edges[e]}`));
  const seaCorrected = split.bare.sea.map((v, e) => v + books.corrections.sea[e]);
  seaCorrected.forEach((v, e) => close(v, sea.transport[e], 1e-12, `sea correction at ${split.edges[e]}`));
});

test('the latent part is wetRun’s method, and on Part 7’s own file it gives Part 7’s numbers', () => {
  const h = part7.header, kept = storage(h), L7 = h.constants.latentHeat, areas = sinAreas(h.constants.radius), wet = waterTransport(h);
  const latent = integrate(latentNet(h.zonalBudget, kept.vapor.map((v, b) => (v + kept.cloud[b]) / L7), L7), areas);
  latent.transport.forEach((v, e) => close(v, wet.latent.transport[e], 1e-12, `Part 7 latent at ${wet.edges[e]}`));
  const total = integrate(h.zonalBudget.absorbedSolar.map((v, b) => v - h.zonalBudget.outgoingLongwave[b] - kept.total[b]), areas);
  total.transport.forEach((v, e) => close(v, wet.total.transport[e], 1e-12, `Part 7 total at ${wet.edges[e]}`));
  const part6 = heatTransport({ bins: 36, absorbedSolar: h.zonalBudget.absorbedSolar, outgoingLongwave: h.zonalBudget.outgoingLongwave }, kept.total, h.constants.radius);
  close(part6.imbalance, total.imbalance, 1e-12, 'the sin-band integral is heatTransport’s');
  const own = new Float64Array(36);
  cells.forEach((c) => { own[Math.min(35, Math.floor((c.lat + 90) / 5))] += c.area; });
  const binned = Array.from(own, (a) => a * sphere / own.reduce((s, v) => s + v, 0));
  const recut = integrate(h.zonalBudget.absorbedSolar.map((v, b) => v - h.zonalBudget.outgoingLongwave[b] - kept.total[b]), binned);
  const shift = Math.max(...recut.transport.map((v, e) => Math.abs(v - wet.total.transport[e])));
  assert.ok(shift < 0.11, `Part 7 recut with the cells’ own areas moves by at most ${shift} PW`);
  close(recut.transport[26], 4.25, 0.005, 'and keeps its 4.2 PW at 40° N'); close(recut.transport[10], -5.03, 0.005, 'and its 5.0 PW at 40° S');
});

test('Part 7’s air and Part 6’s dry sphere are their own runs’, and give the page’s numbers', () => {
  const wet = waterTransport(part7.header);
  wet.total.transport.forEach((v, e) => close(PART7_TRANSPORT.air[e], v, 0.0005, `air at ${wet.edges[e]}`));
  wet.latent.transport.forEach((v, e) => close(PART7_TRANSPORT.latent[e], v, 0.0005, `latent at ${wet.edges[e]}`));
  assert.deepEqual([split.peak('part7', 1).latitude, split.peak('part7', -1).latitude], [40, -40]);
  close(split.peak('part7', 1).transport, 4.2, 0.05, 'Part 7 north'); close(split.peak('part7', -1).transport, -5.0, 0.05, 'Part 7 south');
  close(split.at(0).part7, -0.6, 0.05, 'Part 7 across the equator');
  close(split.at(35).part7Latent, 2.1, 0.05, 'Part 7 vapor at 35° N'); close(split.at(-35).part7Latent, -2.1, 0.05, 'at 35° S');
  assert.deepEqual(split.part6, PART6_TRANSPORT);
  close(split.peak('part6', 1).transport, 1.1, 0.005, 'Part 6');
});

test('the sea’s flow by layer: its water adds up to the free surface’s, its eddies move none, its groups are the file’s, and its cells’ heat is what the page draws', () => {
  const { direct, sea, edges } = split, ov = header.overturning;
  assert.equal(direct.transport[0], 0); assert.equal(direct.transport.at(-1), 0);
  ov.caps.forEach((l, c) => close(direct.mass[c + 1], ov.etaNet[c] + ov.imbalance[c], 0.05, `the layer sum is etaNet plus the imbalance at ${l}`));
  ov.caps.forEach((l, c) => { let eddy = 0; for (let k = 0; k < L; k++) eddy += ov.eddyMass[c * L + k]; close(eddy, 0, 0.3, `eddies move almost no water across ${l}`); });
  const { groups } = layerGroups(header);
  for (const flow of ov.returnFlow) {
    const c = ov.caps.indexOf(flow.latitude);
    assert.deepEqual(flow.groups.map((g) => g.name), ['fixed-depth layers', 'classes lighter than 1026', 'classes 1026 to 1026.8', 'classes from 1026.8']);
    flow.groups.forEach((g, n) => { close(groups[n].mass[c], g.mass, 0.02, `${GROUPS[n]} at ${flow.latitude}`); close(groups[n].eddyMass[c], g.eddyMass, 0.02, 'its eddies'); });
  }
  const top = ov.returnFlow.find((r) => r.latitude === -10).groups[0];
  close(top.southward.mass, -271.78, 0.01, 'the surface layers’ flow south across 10° S, at its edges and snapshots'); close(top.northward.mass, 234.5, 0.01, 'and north'); close(top.mass, -37.29, 0.01, 'a net south of');
  assert.deepEqual(split.peak('direct', 1).latitude, -5);
  const gap = (l) => Math.abs(direct.transport[edges.indexOf(l)] - sea.transport[edges.indexOf(l)]);
  assert.ok(gap(-5) > 5 && gap(5) > 5, `the snapshots’ heat sum is far from the books near the equator: ${gap(-5)}, ${gap(5)}`);
  const cellsHeat = cellHeat(header);
  assert.equal(cellsHeat.length, ov.caps.length);
  ov.caps.forEach((l, c) => close(direct.cells[c + 1], cellsHeat[c], 1e-12, `the page’s cells at ${l}`));
  assert.equal(direct.cells[0], 0); assert.equal(direct.cells.at(-1), 0);
  for (const l of [-85, 85, -80, 80]) close(direct.cells[edges.indexOf(l)], 0, 0.05, `no cells at ${l}`);
  for (const l of [-20, 20]) close(direct.cells[edges.indexOf(l)], sea.transport[edges.indexOf(l)], 0.1, `the cells carry what the books say at ${l}`);
  close(direct.cells[edges.indexOf(-20)], -2.19, 0.005, 'the cells across 20° S'); close(sea.transport[edges.indexOf(-20)], -2.12, 0.005, 'against the books');
});

test('the overturning by layer closes at the floor, and its depths are the section’s', () => {
  const closed = layerOverturning(header), open = layerOverturning(header, null, { closed: false }), ov = header.overturning, M = L + 1;
  ov.caps.forEach((l, c) => {
    assert.equal(closed.psi[c * M], 0);
    let sum = 0;
    for (let k = 0; k < L; k++) sum += ov.mass[c * L + k] + ov.eddyMass[c * L + k];
    close(open.bottom[c], sum, 1e-9, `open at ${l}`); close(closed.bottom[c], sum - ov.imbalance[c], 1e-9, `closed at ${l}`);
    close(closed.bottom[c], ov.etaNet[c], 0.05, `closed to the sea level’s own transport at ${l}`);
    close(closed.depth[c * M + L], 4000, 6, `the floor at ${l}`);
    close(closed.depth[c * M + F], HYBRID_LAYERS.depths[F], 0.1, `the base of the twelve surface layers at ${l}`);
  });
  for (const m of ['March', 'September']) { const flow = layerOverturning(header, m), mo = month(m, 'overturning'); ov.caps.forEach((l, c) => close(flow.bottom[c], mo.etaNet[c], 0.05, `${m} at ${l}`)); }
  const section = layerSection(header), sep = layerSection(header, 'September');
  section.columns.forEach((col, b) => {
    assert.equal(col[0].layer, 0); assert.equal(col[0].top, 0);
    for (let n = 1; n < col.length; n++) assert.ok(col[n].top >= col[n - 1].base - 1e-9, `layers stack at ${lat(b)}`);
    close(header.section.base[b * L + L - 1], 4000, 6, `the column reaches the floor at ${lat(b)}`);
  });
  assert.equal(sep.v, null); assert.ok(sep.columns.every((col) => col.every((c) => c.v === null)));
  const zm = zonalMeans(header), june = zonalMeans(header, 'June');
  assert.equal(zm.name, 'year'); assert.equal(june.samples, 120);
  zm.thermoclineTop.forEach((v, b) => assert.ok(v >= 0 && v < 400, `the 1024.0 isopycnal at ${lat(b)}`));
});

test('the prose numbers: the planet’s year and its books', () => {
  close(year.surfaceTemperature, 288.46, 0.005, 'surface');
  close(year.surfaceTemperature - part7Year.surfaceTemperature, 3.07, 0.005, 'warmer than Part 7');
  close(year.absorbedSolar, 226.2, 0.05, 'sunlight absorbed'); close(year.outgoingLongwave, 217.8, 0.05, 'infrared out'); close(part7Year.absorbedSolar, 216.4, 0.05, 'Part 7’s sunlight'); close(part7Year.outgoingLongwave, 217.6, 0.05, 'Part 7’s infrared');
  close(year.planetaryAlbedo, 0.336, 0.0005, 'reflected'); close(part7Year.planetaryAlbedo, 0.364, 0.0005, 'Part 7 reflected');
  close(1 - year.clearAbsorbedSolar / 340.5, 0.107, 0.0005, 'clear sky reflected'); close(1 - part7Year.clearAbsorbedSolar / 340.5, 0.136, 0.0005, 'Part 7’s');
  close(year.cloudSunlight, -78.0, 0.05, 'clouds in sunlight'); close(year.cloudInfrared, 34.2, 0.05, 'clouds in infrared');
  close(year.balance, 8.4, 0.05, 'gained at the top');
  close(year.seaSurfaceTemperature, 15.3, 0.05, 'the sea’s mean surface');
  const sst = header.zonalBudget.seaSurfaceTemperature.map((t) => t - 273.15);
  close(globalMean(sst, area), year.seaSurfaceTemperature, 0.01, 'the band mean is the daily mean');
  close(Math.max(...sst), 23.1, 0.05, 'warmest band'); assert.equal(lat(sst.indexOf(Math.max(...sst))), 7.5); close(sst[0], 0.9, 0.05, 'south pole'); close(sst[35], 1.7, 0.05, 'north pole');
  assert.deepEqual(sst.map((t, b) => (Math.max(...sst) - t < 0.4 ? lat(b) : null)).filter((v) => v !== null), [-12.5, -7.5, -2.5, 2.5, 7.5, 12.5]);
  close(books.global.total, 8.44, 0.005, 'stored'); close(books.global.ocean, 8.49, 0.005, 'by the sea'); close(books.global.air, -0.06, 0.005, 'by the air');
  close(books.global.toaNet, 8.43, 0.005, 'from the top'); close(books.global.seaInput, 8.51, 0.005, 'into the sea'); close(books.global.surfaceNet, 8.74, 0.005, 'through the surface'); close(books.global.snowFusion, 0.24, 0.005, 'the snow melts');
  close(books.global.residual, 0.003, 0.0005, 'unplaced'); close(books.global.surfaceResidual, -0.015, 0.0005, 'on the sea’s side'); close(books.global.atmosphereResidual, 0.018, 0.0005, 'on the air’s');
  assert.ok(Math.abs(books.global.residual) + Math.abs(books.global.surfaceResidual) + Math.abs(books.global.atmosphereResidual) < 0.04, 'the books close to 0.04 W/m²');
  const sigmaT4 = header.zonalBudget.surfaceTemperature.reduce((sum, t, b) => sum + 5.670374419e-8 * t ** 4 * area[b], 0) / sphere;
  close(1 - globalMean(header.zonalBudget.outgoingLongwave, area) / sigmaT4, 0.4475, 0.0005, 'infrared held back');
  close(year.precipitableWater, 22.4, 0.05, 'vapor'); close(part7Year.precipitableWater, 23.0, 0.05, 'Part 7’s vapor');
  close(year.precipitableWater / part7Year.precipitableWater - 1, -0.030, 0.0005, '3 percent less vapor');
  const interior = header.daily.interiorTemperature;
  close(interior[3649], 3.9986, 0.0005, 'the interior at the start of the record'); close(interior[4014], 4.0151, 0.0005, 'and at its end'); close(interior[4014] - interior[3649], 0.0165, 0.0005, 'the rise over the record');
  close(books.global.ocean * YEAR * DAY / (header.constants.rhoCp * 3700), 0.0165, 0.003, 'the books’ 8.5 W/m² would warm 3,700 m of sea by about that');
});

test('the prose numbers: why it is warmer than Part 7', () => {
  const z8 = header.zonalBudget, z7 = part7.header.zonalBudget, sin = Array.from(area, (a) => a / sphere);
  const part = (values, lo, hi) => values.reduce((sum, v, b) => sum + (Math.abs(lat(b)) > lo && Math.abs(lat(b)) < hi ? v * sin[b] : 0), 0);
  const all = z8.absorbedSolar.map((v, b) => v - z7.absorbedSolar[b]), clear = z8.clearAbsorbedSolar.map((v, b) => v - z7.clearAbsorbedSolar[b]);
  close(part(all, 0, 90), 9.8, 0.05, 'more sunlight absorbed'); close(part(clear, 0, 90), 9.8, 0.05, 'under clear skies');
  close(part(clear, 50, 90), 8.3, 0.05, 'clear-sky poleward of 50°'); close(part(all, 50, 90), 4.3, 0.05, 'all-sky poleward of 50°'); close(part(all, 0, 30), 4.4, 0.05, 'all-sky within 30°');
  close(year.outgoingLongwave - part7Year.outgoingLongwave, 0.14, 0.01, 'more infrared out');
  const warmer = z8.surfaceTemperature.map((t, b) => t - z7.surfaceTemperature[b]);
  close(warmer[0], 18.1, 0.05, 'the south pole'); close(warmer[35], 8.8, 0.05, 'the north pole');
  close(capMean(warmer, 0, 7), 14.2, 0.05, 'poleward of 55° S'); close(capMean(warmer, 29, 36), 5.4, 0.05, 'poleward of 55° N');
  const near = warmer.filter((_, b) => Math.abs(lat(b)) < 10);
  assert.equal(near.length, 4); close(Math.max(...near.map(Math.abs)), 0.6, 0.05, 'within 10° of the equator');
  close(Math.min(...warmer), -1.2, 0.05, 'least'); assert.equal(lat(warmer.indexOf(Math.min(...warmer))), 42.5);
  close(part7Year.ice, 11.9, 0.05, 'Part 7’s ice, percent of the globe'); close(year.ice, 0.057, 0.0005, 'this run’s');
});

test('the prose numbers: the drift through the spin-up and the ice that went', () => {
  const sst = yearlyMeans(header, 'seaSurfaceTemperature'), toa = yearlyMeans(header, 'toaNet'), interior = yearlyMeans(header, 'interiorTemperature');
  assert.equal(sst.length, 11);
  close(sst[0], 16.33, 0.005, 'the first year'); close(sst[10], 15.315, 0.005, 'the recorded year'); close(sst[10] - sst[7], -0.04, 0.005, 'the last three years');
  const d = header.daily;
  close(d.seaSurfaceTemperature[364], 16.18, 0.005, 'day 365'); close(d.seaSurfaceTemperature[4014], 15.22, 0.005, 'day 4015'); close(d.interiorTemperature[364], 3.874, 0.0005, 'interior day 365'); close(d.interiorTemperature[4014], 4.015, 0.0005, 'day 4015');
  close(toa[0], 3.5, 0.05, 'the top in year 1'); close(toa[5], 7.6, 0.05, 'year 6'); close(toa[10], 8.4, 0.05, 'the record');
  assert.ok(toa.every((v, y) => y === 0 || v > toa[y - 1]), 'the planet takes in more each year');
  close(interior[10] - interior[0], 0.136, 0.0005, 'the interior’s rise');
  const surface = yearlyMeans(header, 'surfaceNet');
  close(surface[10], 8.74, 0.005, 'the sea takes in through its surface in the record'); close(surface[5], 7.9, 0.05, 'in year 6');
  const melt = yearlyMeans(header, 'snowFusion');
  close(surface[1] - melt[1], 4.85, 0.005, 'the sea takes in, net of the snow it melts, in the second year'); close(surface[10] - melt[10], 8.51, 0.005, 'and in the record');
  assert.ok(surface.slice(1).every((v) => v > 5), 'the sea took in more than 5 W/m² every year after the first');
  const fastest = yearlyMeans(header, 'fastestMixedLayerCurrent'), fl = d.fastestCurrentLayer.slice(365), layer0 = fl.filter((v) => v === 0).length;
  close(layer0 / fl.length, 0.968, 0.001, 'days the fastest current was in the top layer');
  close(Math.min(...d.fastestMixedLayerCurrent.slice(3650)), 1.19, 0.005, 'the top layer’s fastest current, least over the record'); close(Math.max(...d.fastestMixedLayerCurrent.slice(3650)), 1.84, 0.005, 'greatest');
  close(Math.min(...fastest.slice(1)), 1.17, 0.005, 'in the year means'); close(Math.max(...fastest), 1.53, 0.005, 'and at most');
  const lastIce = (name, over = 0) => { const l = lastDay(name, over); return { day: l, date: date(d.from / DAY + l), year: Math.floor((l - 1) / YEAR) + 1 }; };
  assert.deepEqual(lastIce('iceNorth'), { day: 2595, date: 'April 29', year: 8 });
  assert.deepEqual(lastIce('iceSouth'), { day: 3894, date: 'November 19', year: 11 });
  close(d.ice[0], 6.57, 0.005, 'ice at the start'); close(d.ice[59], 4.63, 0.005, 'day 60'); close(part7.header.daily.ice[59], 8.95, 0.005, 'Part 7’s day 60');
  const peak = (name, y) => Math.max(...d[name].slice(y * YEAR, y * YEAR + YEAR));
  close(peak('ice', 1), 1.52, 0.005, 'the most ice in year 2'); close(peak('ice', 10), 0.25, 0.005, 'in the record'); close(peak('iceNorth', 7), 0.07, 0.005, 'north in year 8');
  close(peak('iceSouth', 10), 0.50, 0.005, 'the southern patch, percent of its hemisphere'); assert.equal(peak('iceNorth', 10), 0);
  const record = d.iceSouth.slice(3650), at = record.indexOf(Math.max(...record));
  assert.equal(date(at + 0.5), 'October 10');
  assert.deepEqual([record.findIndex((v) => v > 0.001), record.findLastIndex((v) => v > 0.001)].map((n) => date(n + 0.5)), ['July 22', 'November 19']);
  const zm = zonalMeans(header);
  close(zm.iceCover[0], 0.174, 0.0005, 'cover at 87.5° S'); close(zm.iceCover[1], 0.057, 0.0005, 'at 82.5° S'); assert.ok(zm.iceCover.slice(2).every((v) => v === 0));
  const cold = header.zonalBudget.seaSurfaceTemperature.map((t) => t - 273.15);
  assert.ok(Math.min(...cold) > 0.9);
});

test('the prose numbers: the sea’s surface through the seasons', () => {
  const zm = zonalMeans(header), eq = (v) => (v[17] + v[18]) / 2;
  close(eq(zm.mixedLayerDepth), 18.4, 0.05, 'equatorial mixed layer'); close(Math.min(...zm.mixedLayerDepth), 18.2, 0.05, 'the shallowest band'); assert.equal(lat(zm.mixedLayerDepth.indexOf(Math.min(...zm.mixedLayerDepth))), -2.5);
  close(Math.max(...zm.mixedLayerDepth), 157, 0.5, 'deepest band in the year’s mean'); assert.equal(lat(zm.mixedLayerDepth.indexOf(Math.max(...zm.mixedLayerDepth))), 87.5);
  close(zm.mixedLayerDepth[bin(-37.5)], 124, 0.5, 'at 37.5° S'); close(zm.mixedLayerDepth[bin(42.5)], 130, 0.5, 'at 42.5° N'); close(zm.mixedLayerDepth[bin(2.5)], 18.6, 0.05, 'at 2.5° N');
  assert.deepEqual(zm.mixedLayerDepth.map((v, b) => (v < HYBRID_LAYERS.depths[F] ? lat(b) : null)).filter((v) => v !== null), [-12.5, -7.5, -2.5, 2.5, 7.5, 12.5, 17.5]);
  close(zm.thermoclineTop[bin(-17.5)], 150, 0.5, '1024.0 at 17.5° S'); close(zm.thermoclineTop[bin(17.5)], 159, 0.5, 'at 17.5° N'); close(zm.thermoclineTop[bin(-2.5)], 67, 0.5, 'at 2.5° S'); close(zm.thermoclineTop[bin(2.5)], 69, 0.5, 'at 2.5° N');
  assert.deepEqual(zm.thermoclineTop.map((v, b) => (v < 5 ? lat(b) : null)).filter((v) => v !== null && Math.abs(v) < 40), [-37.5, -32.5, 32.5, 37.5]);
  const hump = (v, lo, hi) => Math.max(...v.slice(bin(lo), bin(hi) + 1)), ssh = zm.seaSurfaceHeight.map((v) => v * 100);
  close(ssh[bin(-22.5)], 90, 0.5, 'hump at 22.5° S'); close(ssh[bin(22.5)], 80, 0.5, 'at 22.5° N'); close(eq(ssh), -32, 0.5, 'dip at the equator');
  assert.equal(lat(ssh.indexOf(hump(ssh, -45, 0))), -22.5); assert.equal(lat(ssh.indexOf(hump(ssh, 0, 45))), 22.5);
  close(ssh[bin(-57.5)], -104, 0.5, 'low at 57.5° S'); close(ssh[bin(62.5)], -101, 0.5, 'at 62.5° N'); close(Math.max(...ssh) - Math.min(...ssh), 194, 0.5, 'hump to trough');
  const g = 9.81, f = (l) => 2 * 7.2921e-5 * Math.sin(l * D), span = 5 * D * header.constants.radius;
  const balance = (a, b) => -g * (zm.seaSurfaceHeight[bin(b)] - zm.seaSurfaceHeight[bin(a)]) / (f((a + b) / 2) * span) * 100;
  close(ssh[bin(32.5)] - ssh[bin(37.5)], 38, 0.5, 'the fall between 32.5° and 37.5° N'); close(balance(32.5, 37.5), 7.9, 0.05, 'balances a current of cm/s east');
  close((zm.uSurface[bin(32.5)] + zm.uSurface[bin(37.5)]) / 2 * 100, 8.7, 0.05, 'the top layer’s there');
  close(zm.uSurface[bin(-2.5)] * 100, -90, 0.5, 'westward at 2.5° S'); close(zm.uSurface[bin(2.5)] * 100, -86, 0.5, 'at 2.5° N'); close(zm.uSurface[bin(-7.5)] * 100, -49, 0.5, 'at 7.5° S'); close(zm.uSurface[bin(7.5)] * 100, -45, 0.5, 'at 7.5° N');
  for (const l of [-42.5, -37.5, 37.5, 42.5]) assert.ok(zm.uSurface[bin(l)] * 100 > 9.5 && zm.uSurface[bin(l)] * 100 < 12.5, `east at ${l}: ${zm.uSurface[bin(l)] * 100}`);
  close(zm.upwelling[bin(-2.5)], 0.41, 0.005, 'rising at 2.5° S'); close(zm.upwelling[bin(2.5)], 0.32, 0.005, 'at 2.5° N'); close(zm.upwelling[bin(-7.5)], -0.31, 0.005, 'sinking at 7.5° S'); close(zm.upwelling[bin(7.5)], -0.31, 0.005, 'at 7.5° N');
  for (const l of [-37.5, -32.5, -27.5, -22.5, 22.5, 27.5, 32.5, 37.5]) assert.ok(zm.upwelling[bin(l)] <= -0.055 && zm.upwelling[bin(l)] >= -0.085, `sinking at ${l}: ${zm.upwelling[bin(l)]}`);
  const deepest = (m) => { const v = month(m, 'zonalMean').mixedLayerDepth, d = Math.max(...v); return [Math.round(d), lat(v.indexOf(d))]; };
  assert.deepEqual(deepest('March'), [380, 87.5]); assert.deepEqual(deepest('September'), [356, -77.5]); assert.deepEqual(deepest('July'), [265, -77.5]); assert.deepEqual(deepest('January'), [337, 87.5]);
  const times = run.times('year'), k = field(header, 'seaSurfaceTemperature'), at45 = [], at45s = [], tropics = [];
  for (let n = 0; n + 1 < times.length; n++) { const v = run.values('year', n, k); at45.push(band(v, 42.5, 47.5)); at45s.push(band(v, -47.5, -42.5)); tropics.push(band(v, -5, 5)); }
  const swing = (s) => Math.max(...s) - Math.min(...s);
  close(swing(at45), 3.2, 0.05, 'the sea’s swing at 45° N'); assert.equal(date(times[at45.indexOf(Math.max(...at45))] / DAY), 'September 1'); assert.equal(date(times[at45.indexOf(Math.min(...at45))] / DAY), 'March 30');
  close(Math.min(...at45), 8.06, 0.005, 'its coolest'); close(Math.max(...at45), 11.27, 0.005, 'its warmest'); close(swing(at45s), 2.8, 0.05, 'and at 45° S'); close(swing(tropics), 1.4, 0.05, 'within 5° of the equator');
  const part7k = field(part7.header, 'surfaceTemperature'), p7 = [], p7s = [];
  for (let n = 0; n + 1 < part7.times('year').length; n++) { const v = part7.values('year', n, part7k); p7.push(band(v, 42.5, 47.5)); p7s.push(band(v, -47.5, -42.5)); }
  close(swing(p7), 25.6, 0.05, 'Part 7’s surface at 45° N'); close(swing(p7s), 29.1, 0.05, 'and at 45° S');
  const sss = zm.sss;
  close(Math.max(...sss.slice(0, 18)), 35.30, 0.005, 'saltiest south'); assert.equal(lat(sss.indexOf(Math.max(...sss.slice(0, 18)))), -27.5);
  close(Math.max(...sss.slice(18)), 35.27, 0.005, 'saltiest north'); assert.equal(lat(18 + sss.slice(18).indexOf(Math.max(...sss.slice(18)))), 27.5);
  const gap = header.zonalBudget.evaporation.map((v, b) => v - header.zonalBudget.precipitation[b]);
  close(gap[bin(-22.5)], 1.68, 0.005, 'evaporation over rain at 22.5° S'); close(gap[bin(22.5)], 1.60, 0.005, 'at 22.5° N');
  assert.deepEqual([gap.indexOf(Math.max(...gap.slice(0, 18))), 18 + gap.slice(18).indexOf(Math.max(...gap.slice(18)))], [bin(-22.5), bin(22.5)]);
  close(Math.min(...sss), 33.98, 0.005, 'the freshest band');
  let largest = 0; const km = field(header, 'mixedLayerDepth');
  times.forEach((_, n) => run.values('year', n, km).forEach((v) => { largest = Math.max(largest, v); }));
  close(largest, 485, 1, 'no cell goes deeper');
});

test('the prose numbers: the rain band stays near the equator', () => {
  const wander = (r) => {
    const times = r.times('year'), k = field(r.header, 'rain'), out = [];
    for (let n = 1; n < times.length; n++) out.push({ day: times[n] / DAY - 2.5, ...peakLatitude(zonal(r.values('year', n, k)), 12, 24) });
    return [out.reduce((a, b) => (b.latitude > a.latitude ? b : a)), out.reduce((a, b) => (b.latitude < a.latitude ? b : a))];
  };
  const [north, south] = wander(run), [north7, south7] = wander(part7);
  close(north.latitude, 15.9, 0.05, 'farthest north'); assert.equal(date(north.day), 'September 24');
  close(south.latitude, -17.8, 0.05, 'farthest south'); assert.equal(date(south.day), 'January 2');
  close(north7.latitude, 23.3, 0.05, 'Part 7 north'); close(south7.latitude, -18.4, 0.05, 'Part 7 south');
  const monthly = header.zonalBudget.months.map((m) => peakLatitude(m.precipitation, 12, 24).latitude);
  close(Math.min(...monthly), -12.9, 0.05, 'monthly south'); close(Math.max(...monthly), 12.6, 0.05, 'monthly north');
  const annual = peakLatitude(header.zonalBudget.precipitation);
  close(annual.latitude, 11.7, 0.05, 'the year’s band'); close(annual.value, 4.58, 0.005, 'its rain');
  close(peakLatitude(header.zonalBudget.precipitation, 24, 36).latitude, 44.1, 0.05, 'northern storm track'); close(peakLatitude(header.zonalBudget.precipitation, 0, 12).latitude, -45.7, 0.05, 'southern');
});

test('the prose numbers: the currents and the overturning cells', () => {
  const zm = zonalMeans(header), runs = [];
  zm.zonalTransport.forEach((v, b) => { const s = Math.sign(v); if (runs.length && runs.at(-1).sign === s) runs.at(-1).sum += v; else runs.push({ sign: s, from: lat(b) - 2.5, sum: v }); });
  assert.deepEqual(runs.map((r) => [r.from, Math.round(r.sum) + 0]), [[-90, -490], [-55, 320], [-30, -2412], [30, 275], [60, -180]]);
  close(zm.uDepthMean[17] * 100, -18.7, 0.05, 'depth-mean at 2.5° S'); close(zm.uDepthMean[18] * 100, -18.2, 0.05, 'at 2.5° N');
  const section = layerSection(header), col = section.columns[bin(2.5)], at = (k) => col.find((c) => c.layer === k);
  close(at(0).u, -86, 0.5, 'top layer'); close(at(11).u, -70, 0.5, 'the twelfth layer, to 57 m'); close(at(20).u, -61, 0.5, 'to 144 m'); close(at(29).u, -56, 0.5, 'to 234 m'); close(at(38).u, -58, 0.5, 'to 428 m'); close(at(42).u, -27, 0.5, 'to 748 m');
  close(at(44).u, -14, 0.5, 'to 1,141 m'); close(at(48).u, -9, 0.5, 'the deepest layer');
  close(at(48).top, 1991, 1, 'its top'); close(at(48).base, 3995, 1, 'its base');
  for (const b of [bin(-2.5), bin(2.5)]) assert.ok(section.columns[b].every((c) => c.u < -9 && c.u > -91), `every layer on the equator runs west at ${lat(b)}`);
  assert.deepEqual(col.filter((c) => c.thickness < 1).map((c) => c.label), [1026.8, 1026.83, 1026.86, 1026.89, 1026.92, 1026.95, 1026.98, 1027.1, 1027.5]);
  const flow = layerOverturning(header), M = L + 1, under = (l) => flow.psi[flow.caps.indexOf(l) * M + F];
  close(under(-5), -107.3, 0.05, 'south across 5° S'); close(under(5), 102.5, 0.05, 'north across 5° N'); close(under(0), 13.7, 0.05, 'across the equator');
  close(under(-10), -37.3, 0.05, 'south across 10° S'); close(under(10), 32.1, 0.05, 'north across 10° N');
  close(under(-45), 31.4, 0.05, 'north across 45° S'); close(under(45), -26.9, 0.05, 'south across 45° N'); close(under(-35), 14.2, 0.05, 'north across 35° S'); close(under(35), -12.3, 0.05, 'south across 35° N');
  close(flow.depth[flow.caps.indexOf(-5) * M + F], 56.9, 0.05, 'the base of the twelve layers');
  const psiAt = (f, l, z) => { const c = f.caps.indexOf(l); let j = 0; while (j < L - 1 && f.depth[c * M + j + 1] < z) j++; const d0 = f.depth[c * M + j], d1 = f.depth[c * M + j + 1], t = Math.min(1, Math.max(0, (z - d0) / (d1 - d0))); return f.psi[c * M + j] + t * (f.psi[c * M + j + 1] - f.psi[c * M + j]); };
  for (const [l, expected] of [[-5, [-107.2, -96.2, -96.2, -6.9, -26.0, -14.1, -6.7]], [5, [102.5, 92.9, 91.5, -7.1, 39.0, 15.3, 5.8]]]) [57, 100, 300, 500, 1000, 2000, 3000].forEach((z, n) => close(psiAt(flow, l, z), expected[n], 0.1, `${l}° at ${z} m`));
  const peak = (sign, from, to, deep = 700) => { let best = { psi: 0 }; flow.caps.forEach((l, c) => { if (l >= from && l <= to) for (let k = 0; k <= L; k++) if (flow.depth[c * M + k] <= deep && sign * flow.psi[c * M + k] > sign * best.psi) best = { psi: flow.psi[c * M + k], l, depth: flow.depth[c * M + k] }; }); return best; };
  const w1 = peak(1, -60, -30), w2 = peak(-1, 30, 60);
  assert.deepEqual([w1.l, Math.round(w1.depth)], [-45, 351]); close(w1.psi, 40.0, 0.05, 'the southern westerly cell'); assert.deepEqual([w2.l, Math.round(w2.depth)], [40, 289]); close(w2.psi, -36.6, 0.05, 'the northern');
  close(psiAt(flow, -45, 500), 14.7, 0.05, 'the southern westerly cell at 500 m'); close(psiAt(flow, 45, 500), -12.5, 0.05, 'the northern');
  close(psiAt(flow, -45, 1000), 1.6, 0.05, 'at 1,000 m'); close(psiAt(flow, 45, 1000), -3.7, 0.05, 'and the north');
  const deep = (sign, from, to) => { let best = { psi: 0 }; flow.caps.forEach((l, c) => { if (l >= from && l <= to) for (let k = 0; k <= L; k++) if (flow.depth[c * M + k] >= 1000 && sign * flow.psi[c * M + k] > sign * best.psi) best = { psi: flow.psi[c * M + k], l, depth: flow.depth[c * M + k] }; }); return best; };
  close(deep(-1, -90, 0).psi, -33.8, 0.05, 'the southern deep cell'); assert.equal(deep(-1, -90, 0).l, -10); close(deep(1, 0, 90).psi, 37.6, 0.05, 'the northern'); assert.equal(deep(1, 0, 90).l, 5);
  let across = 0, weak = 0;
  for (const z of Array.from({ length: 400 }, (_, i) => 10 + 10 * i)) { across = Math.max(across, Math.abs(psiAt(flow, 0, z))); for (const l of flow.caps) if (Math.abs(l) >= 60) weak = Math.max(weak, Math.abs(psiAt(flow, l, z))); }
  close(across, 14.7, 0.1, 'the most the year’s flow carries across the equator at any depth'); close(weak, 6.7, 0.1, 'the most any cell poleward of 60° carries');
  const sep = layerOverturning(header, 'September'), mar = layerOverturning(header, 'March'), sepUnder = (l) => sep.psi[sep.caps.indexOf(l) * M + F], marUnder = (l) => mar.psi[mar.caps.indexOf(l) * M + F];
  close(sepUnder(-5), -237.3, 0.05, 'September, 5° S'); close(sepUnder(5), 46.0, 0.05, 'September, 5° N'); close(marUnder(5), 230.7, 0.05, 'March, 5° N'); close(marUnder(-5), -31.5, 0.05, 'March, 5° S');
  close(sepUnder(0), -34.0, 0.05, 'September across the equator, top 57 m'); close(marUnder(0), 88.7, 0.05, 'March');
  close(psiAt(sep, 0, 300), -90.1, 0.05, 'September across the equator, above 300 m'); close(psiAt(mar, 0, 100), 110.5, 0.05, 'March, above 100 m');
});

test('the prose numbers: the heat carried, split between the air and the sea', () => {
  const s = split, p = (which, sign) => s.peak(which, sign);
  assert.deepEqual([p('total', 1).latitude, p('total', -1).latitude], [30, -30]);
  close(p('total', 1).transport, 4.80, 0.005, 'in all, north'); close(p('total', -1).transport, -4.86, 0.005, 'south');
  assert.deepEqual([p('air', 1).latitude, p('air', -1).latitude], [35, -35]);
  close(p('air', 1).transport, 3.89, 0.005, 'the air north'); close(p('air', -1).transport, -3.93, 0.005, 'south');
  close(s.at(40).air, 3.85, 0.005, 'the air at 40° N'); close(s.at(-40).air, -3.77, 0.005, 'at 40° S');
  close(s.at(40).sea, 0.40, 0.005, 'the sea at 40° N'); close(s.at(-40).sea, -0.43, 0.005, 'at 40° S');
  close(s.at(40).total, 4.24, 0.005, 'in all at 40° N'); close(s.at(-40).total, -4.20, 0.005, 'at 40° S');
  close(s.at(35).sea, 0.64, 0.005, 'the sea at 35° N'); close(s.at(-35).sea, -0.67, 0.005, 'at 35° S'); close(s.at(35).total, 4.53, 0.005, 'in all at 35° N'); close(s.at(-35).total, -4.61, 0.005, 'at 35° S');
  close(s.at(35).sea / s.at(35).total, 0.141, 0.0005, 'the sea’s share across 35° N'); close(s.at(-35).sea / s.at(-35).total, 0.146, 0.0005, 'across 35° S');
  close(s.at(30).sea / s.at(30).total, 0.26, 0.005, 'across 30° N'); close(s.at(-30).sea / s.at(-30).total, 0.25, 0.005, 'across 30° S');
  close(s.at(-10).sea / s.at(-10).total, 0.72, 0.005, 'across 10° S'); close(s.at(10).sea / s.at(10).total, 0.74, 0.005, 'across 10° N');
  close(s.at(-20).sea / s.at(-20).total, 0.47, 0.005, 'across 20° S'); close(s.at(20).sea / s.at(20).total, 0.48, 0.005, 'across 20° N');
  close(s.at(-40).sea / s.at(-40).total, 0.10, 0.005, 'across 40° S'); close(s.at(40).sea / s.at(40).total, 0.09, 0.005, 'across 40° N');
  close(p('sea', 1).transport, 2.42, 0.005, 'the sea north'); close(p('sea', -1).transport, -2.26, 0.005, 'south'); assert.deepEqual([p('sea', 1).latitude, p('sea', -1).latitude], [15, -15]);
  close(s.at(-10).sea, -2.21, 0.005, 'the sea at 10° S'); close(s.at(10).sea, 2.19, 0.005, 'at 10° N'); close(s.at(20).sea, 2.15, 0.005, 'at 20° N'); close(s.at(-20).sea, -2.12, 0.005, 'at 20° S');
  close(s.at(0).total, -0.01, 0.005, 'in all across the equator'); close(s.at(0).sea, 0.03, 0.005, 'the sea'); close(s.at(0).air, -0.04, 0.005, 'the air');
  const crossings = []; for (let e = 1; e < s.edges.length - 1; e++) { const d0 = Math.abs(s.sea.transport[e]) - Math.abs(s.air.transport[e]), d1 = Math.abs(s.sea.transport[e + 1]) - Math.abs(s.air.transport[e + 1]); if (d0 * d1 < 0) crossings.push(s.edges[e] + 5 * d0 / (d0 - d1)); }
  assert.equal(crossings.length, 4); close(crossings[0], -18.7, 0.05, 'the sea carries more than the air from'); close(crossings[3], 19.1, 0.05, 'to');
  close(p('latent', 1).transport, 1.87, 0.005, 'vapor north'); close(p('latent', -1).transport, -1.97, 0.005, 'vapor south'); assert.deepEqual([p('latent', 1).latitude, p('latent', -1).latitude], [35, -35]);
  close(p('dry', 1).transport, 2.49, 0.005, 'dry air north'); close(p('dry', -1).transport, -2.49, 0.005, 'south'); assert.deepEqual([p('dry', 1).latitude, p('dry', -1).latitude], [20, -20]);
  close(p('bareTotal', 1).transport, 4.67, 0.005, 'from the top of the air alone, north'); close(p('bareTotal', -1).transport, -4.58, 0.005, 'south'); assert.deepEqual([p('bareTotal', 1).latitude, p('bareTotal', -1).latitude], [35, -35]);
  close(s.at(0).bareTotal, 0.08, 0.005, 'and across the equator');
  const input = books.seaInput, stored = books.sea;
  close(capMean(input, bin(-2.5), bin(2.5) + 1), 50.2, 0.05, 'the sea takes in within 5° of the equator, W/m²'); close(input[bin(-2.5)], 50.6, 0.05, 'at 2.5° S'); close(input[bin(2.5)], 49.7, 0.05, 'at 2.5° N');
  close(capMean(input, bin(-7.5), bin(7.5) + 1) * (area.slice(bin(-7.5), bin(7.5) + 1).reduce((a, b) => a + b, 0)) / 1e15, 3.6, 0.05, 'PW between 10° S and 10° N');
  close(capMean(stored, bin(-2.5), bin(2.5) + 1), -17.4, 0.05, 'the equatorial sea stores, W/m²');
  close(capMean(input, bin(-2.5), bin(2.5) + 1) - capMean(stored, bin(-2.5), bin(2.5) + 1), 67.6, 0.1, 'and the currents carry away');
  close(capMean(input, bin(-2.5), bin(2.5) + 1) * area.slice(bin(-2.5), bin(2.5) + 1).reduce((a, b) => a + b, 0) / 1e15, 2.43, 0.005, 'PW into the sea within 5° of the equator'); close(capMean(stored, bin(27.5), bin(37.5) + 1), 22.5, 0.05, 'the northern subtropical sea 25° to 40° N'); close(capMean(stored, bin(-37.5), bin(-27.5) + 1), 25.0, 0.05, 'and the southern');
  close(capMean(input, bin(27.5), bin(37.5) + 1), -1.1, 0.05, 'while the northern subtropics take in'); close(capMean(input, bin(-37.5), bin(-27.5) + 1), 2.4, 0.05, 'and the southern');
  close(capMean(stored, 30, 36), -0.3, 0.05, 'north of 60° N'); close(capMean(stored, 0, 6), -0.5, 0.05, 'south of 60° S'); close(capMean(input, 30, 36), -9.0, 0.05, 'taking in'); close(capMean(input, 0, 6), -10.5, 0.05, 'south');
  close(capMean(stored, 30, 36) - capMean(input, 30, 36), 8.7, 0.05, 'the currents bring the northern polar sea'); close(capMean(stored, 0, 6) - capMean(input, 0, 6), 10.0, 0.05, 'the southern');
  const bz = header.zonalBudget, per = (name, l) => bz[name][bin(l)];
  close(per('surfaceShortwave', 2.5), 191.8, 0.05, 'sunlight at 2.5° N'); close(per('surfaceShortwave', 17.5), 174.2, 0.05, 'at 17.5° N'); close(per('latentHeat', 2.5), 88.1, 0.05, 'latent heat at 2.5° N'); close(per('latentHeat', 17.5), 114.6, 0.05, 'at 17.5° N'); close(per('latentHeat', -17.5), 115.3, 0.05, 'at 17.5° S');
  close(per('evaporation', 2.5), 3.04, 0.005, 'evaporation at 2.5° N'); close(per('evaporation', 17.5), 3.96, 0.005, 'at 17.5° N'); close(per('evaporation', -17.5), 3.98, 0.005, 'at 17.5° S'); close(per('evaporation', -2.5), 3.07, 0.005, 'at 2.5° S');
  const ir = (l) => per('surfaceShortwave', l) - per('sensibleHeat', l) - per('latentHeat', l) - per('netSurfaceHeat', l);
  close(ir(2.5), 45.4, 0.05, 'infrared at 2.5° N'); close(ir(17.5), 45.3, 0.05, 'at 17.5° N');
  close(books.global.residual * sphere / 2 / 1e15, 0.0007, 0.0003, 'the residual put in one hemisphere moves the equator by, PW');
});

test('the prose numbers: the water cycle', () => {
  const g = water.global;
  close(g.evaporation, 2.67, 0.005, 'evaporation'); close(g.precipitation, 2.67, 0.005, 'rain');
  close(g.evaporation / part7Year.evaporation - 1, 0.003, 0.003, 'the same as Part 7’s 2.66');
  close(g.latentHeat, 77.3, 0.05, 'latent heat'); close(g.sensibleHeat, 12.7, 0.05, 'sensible heat'); close(g.latentHeat / g.sensibleHeat, 6.1, 0.05, 'times');
  close(g.convectiveRain / g.precipitation, 0.64, 0.005, 'from the plume');
  close(year.residenceDays, 8.4, 0.05, 'days of rain in the air');
  assert.deepEqual(water.gapCrossings.map((l) => +l.toFixed(0)), [-37, -15, -5, 5, 16, 37]);
});

test('the prose numbers: the spin-up’s first current and the polar sea', () => {
  const st = run.times('spinup'), east = field(header, 'currentEast');
  const equator = (day) => band(run.values('spinup', st.findIndex((t) => Math.round(t / DAY + header.spinupDays) === day), east), -5, 5) * 100;
  close(equator(1), 5.8, 0.05, 'eastward on day 1'); close(equator(3), 20.2, 0.05, 'at most on day 3'); close(equator(10), 16.4, 0.05, 'on day 10'); close(equator(25), 6.6, 0.05, 'on day 25');
  assert.ok(equator(31) > 0 && equator(32) < 0, 'westward from day 32'); close(equator(60), -16.6, 0.05, 'westward on day 60');
  for (let d = 1; d <= 31; d++) assert.ok(equator(d) > 0, `east on day ${d}`);
  const times = run.times('year'), sst = field(header, 'seaSurfaceTemperature'), sss = field(header, 'surfaceSalinity');
  let coldest = Infinity, densest = 0;
  times.forEach((_, k) => { const t = run.values('year', k, sst), s = run.values('year', k, sss); t.forEach((v, i) => { coldest = Math.min(coldest, v); densest = Math.max(densest, seawaterDensity(v + 273.15, s[i])); }); });
  close(coldest, -1.86, 0.005, 'the coldest cell, °C'); close(densest, 1026.60, 0.005, 'the densest surface water');
  assert.ok(densest < header.classes[44], 'lighter than the 1026.65 layer');
  close(seawaterDensity(275.15, 34.7), 1026.77, 0.005, 'the start’s deep water'); assert.equal(header.classes[48], 1026.77);
  const tops = layerSection(header).columns.map((c) => c.find((x) => x.layer === 48));
  tops.forEach((x, b) => { assert.ok(x.T > 1.96 && x.T < 2.1 && x.S > 34.68 && x.S < 34.705, `the deepest layer at ${lat(b)}`); if (Math.abs(lat(b)) < 30) assert.ok(x.top > 1900 && x.top < 2340, `its top at ${lat(b)}`); if (Math.abs(lat(b)) > 60) assert.ok(x.top > 700 && x.top < 1000, `its top at ${lat(b)}`); });
  close(tops[0].top, 712, 1, 'the deepest layer’s top at the South Pole'); close(tops[35].top, 713, 1, 'and the North Pole'); close(Math.min(...tops.filter((_, b) => Math.abs(lat(b)) < 30).map((x) => x.top)), 1914, 1, 'its shallowest within 30°'); close(Math.max(...tops.map((x) => x.top)), 2339, 1, 'and deepest');
  const deep = layerSection(header).columns.map((c) => c.at(-1));
  assert.ok(deep.every((c) => c.layer === 57 && c.T > 1.9 && c.T < 2.05), 'the sheets on the floor hold the start’s deep water');
});

test('the prose numbers: the section’s layers, its salt, and the cross-equatorial flow of September', () => {
  const section = layerSection(header), column = (l) => section.columns[section.centers.indexOf(l)], at = (l, k) => column(l).find((c) => c.layer === k);
  for (const [l, depth] of [[-27.5, 784], [-22.5, 774], [22.5, 735], [27.5, 765]]) {
    const salty = column(l).filter((x) => x.S > 35);
    close(Math.max(...salty.map((x) => x.base)), depth, 5, `saltier than 35 psu down to ${depth} m at ${l}`);
    assert.ok(salty.every((x) => x.top < depth), `and no deeper at ${l}`);
  }
  close(column(-22.5)[0].S, 35.24, 0.005, 'surface salinity at 22.5° S'); close(column(22.5)[0].S, 35.19, 0.005, 'at 22.5° N');
  close(column(-87.5)[0].S, 33.98, 0.005, 'the south polar top layer'); close(column(87.5)[0].S, 34.11, 0.005, 'the north polar'); close(column(-87.5).at(-1).S, 34.70, 0.005, 'over the deep water'); close(column(87.5).at(-1).S, 34.70, 0.005, 'north');
  const probe = (k) => at(2.5, k);
  close(probe(0).T, 22.8, 0.05, 'the top layer at 2.5° N'); close(probe(29).T, 13.2, 0.05, 'the layer to 234 m'); close(probe(48).T, 2.1, 0.05, 'the deepest');
  close(probe(29).top, 224, 0.5, 'its top'); close(probe(29).base, 234, 0.5, 'its base');
  const water29 = seawaterDensity(probe(29).T + 273.15, probe(29).S);
  close(water29, 1025.4, 0.05, 'the water in the layer labeled 1024.0'); assert.ok(water29 - 1024 > 1.3, 'denser than its label');
  close(probe(7).T, 22.1, 0.05, 'the probe’s opening layer'); close(probe(7).u, -77, 0.5, 'its current'); close(probe(7).top, 27, 0.5, 'its top'); close(probe(7).base, 32, 0.5, 'its base');
  assert.deepEqual([probe(7).label, probe(7).mode], [1018.25, 'fixed-depth']);
  close(probe(0).u, -86, 0.5, 'the surface current at 2.5° N');
  const iso = zonalMeans(header).thermoclineTop[bin(2.5)];
  close(iso, 69, 0.5, 'the 1024.0 isopycnal at 2.5° N');
  const sep = layerOverturning(header, 'September'), M = L + 1, c = sep.caps.indexOf(0);
  const psi = (z) => { let j = 0; while (j < L - 1 && sep.depth[c * M + j + 1] < z) j++; const d0 = sep.depth[c * M + j], d1 = sep.depth[c * M + j + 1], t = Math.min(1, Math.max(0, (z - d0) / (d1 - d0))); return sep.psi[c * M + j] + t * (sep.psi[c * M + j + 1] - sep.psi[c * M + j]); };
  close(psi(100), -72.2, 0.05, 'south above 100 m'); close(psi(300), -90.1, 0.05, 'above 300 m'); close(psi(500), -42.5, 0.05, 'above 500 m'); close(psi(1000), -67.5, 0.05, 'above 1,000 m'); close(psi(3000), -37.2, 0.05, 'above 3,000 m');
  close(sep.imbalance[c], -1.8, 0.05, 'what the samples miss there');
  const gap = (m) => Math.max(...Array.from(layerOverturning(header, m).imbalance, Math.abs));
  close(gap(), 0.51, 0.005, 'the samples miss up to 0.5 Sv over the year'); close(Math.max(...header.overturning.months.map((m) => gap(m.name))), 2.61, 0.005, 'and up to 2.6 Sv in a month');
});

test('the prose numbers: the sea’s share, against Part 7 and the Earth, and how its cells carry it', () => {
  const s = split, at = (l) => s.at(l);
  close(at(35).bareTotal, 4.67, 0.005, 'the top alone at 35° N'); close(at(-35).bareTotal, -4.58, 0.005, 'at 35° S');
  const z7 = part7.header.zonalBudget, top7 = integrate(z7.absorbedSolar.map((v, b) => v - z7.outgoingLongwave[b]), area).transport, e = s.edges;
  close(top7[e.indexOf(35)], 4.6, 0.05, 'Part 7 from the top alone at 35° N'); close(top7[e.indexOf(-35)], -4.64, 0.005, 'and at 35° S');
  close(at(35).bareTotal - top7[e.indexOf(35)], 0.06, 0.005, 'more across 35° N'); close(at(-35).bareTotal - top7[e.indexOf(-35)], 0.06, 0.005, 'across 35° S');
  const ov = header.overturning, rc = header.constants.rhoCp;
  const poleward = (l) => {
    const c = ov.caps.indexOf(l), sign = Math.sign(layerOverturning(header).psi[c * (L + 1) + F] || l), sec = header.section;
    let out = 0, back = 0, ho = 0, hb = 0;
    for (let k = 0; k < L; k++) {
      const m = ov.mass[c * L + k] + ov.eddyMass[c * L + k], th = (b) => sec.base[b * L + k] - (k ? sec.base[b * L + k - 1] : 0), T0 = sec.T[c * L + k], T1 = sec.T[(c + 1) * L + k];
      const T = T0 === null ? T1 ?? 0 : T1 === null ? T0 : (T0 * th(c) + T1 * th(c + 1)) / (th(c) + th(c + 1));
      if (Math.sign(m) === sign) { out += m; ho += m * T; } else { back += m; hb += m * T; }
    }
    return { out, back, T: ho / out, Tback: hb / back, heat: rc * (ho + hb) * 1e-9 };
  };
  const ten = poleward(-10);
  close(ten.out, -97, 0.5, 'water going south across 10° S, Sv'); close(ten.back, 97, 0.5, 'coming back'); close(ten.T, 17.2, 0.05, 'at °C'); close(ten.Tback, 10.2, 0.05, 'returning at °C'); close(ten.heat, -2.83, 0.005, 'the cells carry');
  close(poleward(10).T, 17.3, 0.05, 'north across 10° N'); close(poleward(10).Tback, 8.1, 0.05, 'returning at'); close(poleward(10).heat, 3.13, 0.005, 'carrying');
  close(poleward(20).T, 19.9, 0.05, 'across 20° N, °C'); close(poleward(20).Tback, 8.1, 0.05, 'and back at'); close(poleward(20).out, 46, 0.5, 'Sv');
  const c = (l) => at(l).cells;
  close(c(-35), 0.30, 0.005, 'the cells carry toward the equator across 35° S'); close(c(35), -0.29, 0.005, 'across 35° N'); close(c(-40), 0.40, 0.005, 'across 40° S'); close(c(40), -0.48, 0.005, 'across 40° N');
  close(c(30), 0.54, 0.005, 'across 30° N, poleward'); close(c(-30), -0.50, 0.005, 'across 30° S');
  close(c(-5), -4.82, 0.005, 'the cells across 5° S'); close(at(-5).sea, -1.65, 0.005, 'the sea'); close(c(5), 4.65, 0.005, 'across 5° N'); close(at(5).sea, 1.62, 0.005, 'the sea');
  close(c(-10) - at(-10).sea, -0.62, 0.005, 'within the layers across 10° S'); close(at(35).sea - c(35), 0.92, 0.005, 'the sea’s poleward heat within layers across 35° N'); close(at(40).sea - c(40), 0.87, 0.005, 'across 40° N');
  close(at(5).sea - c(5), -3.02, 0.005, 'the sea less its cells at 5° N');
});
