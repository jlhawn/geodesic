import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
import { readFrames } from '../explain/frames.module.js';
import { heatTransport } from '../explain/sunlitRun.module.js';
import { storage, waterTransport, recordedYear as wetYear } from '../explain/wetRun.module.js';
import {
  DATA_URL, PART7_TRANSPORT, PART6_TRANSPORT, GROUPS, OUTCROP, calendarDate, peakLatitude, recordedYear, globalMean, integrate, latentNet,
  heatBooks, transportSplit, classGroups, waterCycle, yearlyMeans, zonalMeans, classSection, classOverturning,
} from '../explain/oceanRun.module.js';

const open = (name) => { const file = new URL(`../explain/data/${name}`, import.meta.url), bytes = readFileSync(file); return readFrames(file.href, bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength)); };
const run = open('ocean_N16.bin'), part7 = open('water_N16.bin');
const { header } = run, DAY = 86400, YEAR = 365, D = Math.PI / 180, L = header.layers;
const cells = [...new Grid(header.N)].map((c) => ({ lat: Math.asin(c.centerVertex.z) / D, area: c.area }));
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const date = (day) => { const c = calendarDate(day); return `${c.month} ${c.date}`; };
const zonal = (values) => { const sum = new Float64Array(36), area = new Float64Array(36); cells.forEach((c, i) => { const b = Math.min(35, Math.floor((c.lat + 90) / 5)); sum[b] += values[i] * c.area; area[b] += c.area; }); return Array.from(sum, (v, b) => v / area[b]); };
const band = (values, lo, hi) => { let sum = 0, area = 0; cells.forEach((c, i) => { if (c.lat >= lo && c.lat < hi) { sum += values[i] * c.area; area += c.area; } }); return sum / area; };
const field = (h, name) => h.fields.findIndex((f) => f.name === name);
const sinAreas = (radius) => Array.from({ length: 36 }, (_, b) => 2 * Math.PI * radius * radius * (Math.sin((-85 + 5 * b) * D) - Math.sin((-90 + 5 * b) * D)));
const area = header.zonalBudget.binArea, sphere = 4 * Math.PI * header.constants.radius ** 2;
const books = heatBooks(header), split = transportSplit(header), water = waterCycle(header), year = recordedYear(header), part7Year = wetYear(part7.header);
const lat = (b) => -87.5 + 5 * b;
const month = (name, part = 'zonalBudget') => header[part].months.find((m) => m.name === name);

test('the run is Part 7’s wet planet with its sea let go: ten years of spin-up and a recorded year from the March equinox', () => {
  assert.equal(new URL(DATA_URL).pathname.endsWith('/explain/data/ocean_N16.bin'), true);
  assert.deepEqual([header.case, header.engine, header.N, header.K, header.layers, header.dt, header.grid, header.spinupDays, header.recordDays], ['ocean', 'gpu', 16, 36, 45, 450, 'bl36', 3650, 365]);
  assert.equal(header.options.ocean.everySteps, 4);
  assert.deepEqual(header.ice.options, part7.header.ice.options);
  assert.deepEqual(header.floor, { kind: 'flat', depth: 4000, minimum: 4000, maximum: 4000, land: 'none: every cell and edge is sea' });
  assert.equal(header.recordStart, 3650 * DAY);
  assert.deepEqual(header.sequences.map((s) => [s.name, s.times.length]), [['year', 74], ['spinup', 61]]);
  assert.equal(header.daily.surfaceTemperature.length, header.spinupDays + header.recordDays);
  for (const counts of [header.nonFinite, header.clipped]) assert.ok(Object.values(counts).every((n) => n === 0), 'every frame value finite and in range');
  assert.ok(!header.resumeCheck.forced && header.resumeCheck.sameEngine && Object.values(header.resumeCheck.stores).every((v) => v === 0), 'the resumed legs are one run');
  assert.equal(header.classes.length, L);
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
  for (const name of ['toaNet', 'surfaceNet', 'snowFusion']) close(books.global[name], header.books.record[name], 5e-6, name);
  for (const name of ['residual', 'surfaceResidual', 'atmosphereResidual']) close(books.global[name], header.books.record[name], 5e-6, name);
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
  close(total.imbalance, -header.books.record.residual, 5e-6, 'the total spreads the books’ residual');
  close(sea.imbalance, -header.books.record.surfaceResidual, 5e-6, 'the sea its surface side');
  close(air.imbalance, -header.books.record.atmosphereResidual, 5e-6, 'the air its own');
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

test('the sea’s own transport, added up by class, agrees with the books', () => {
  const { direct, sea, edges } = split, ov = header.overturning;
  assert.equal(direct.transport[0], 0); assert.equal(direct.transport.at(-1), 0);
  let worst = 0;
  edges.forEach((l, e) => {
    const gap = Math.abs(direct.transport[e] - sea.transport[e]);
    worst = Math.max(worst, gap);
    assert.ok(gap <= 0.25, `within 0.25 PW at ${l}: ${direct.transport[e]} vs ${sea.transport[e]}`);
    if (Math.abs(l) >= 40) assert.ok(gap <= 0.05, `within 0.05 PW poleward of 40° at ${l}: ${gap}`);
    if (Math.abs(sea.transport[e]) > 0.3) assert.equal(Math.sign(direct.transport[e]), Math.sign(sea.transport[e]), `same way at ${l}`);
  });
  close(worst, 0.24, 0.01, 'the largest gap, at 10° S');
  close(direct.perKelvin[edges.indexOf(-10)], 0.034, 0.0005, 'a kelvin of reference moves 10° S by');
  ov.caps.forEach((l, c) => close(direct.mass[c + 1], ov.etaNet[c] + ov.imbalance[c], 0.02, `the class sum is etaNet plus the imbalance at ${l}`));
  ov.caps.forEach((l, c) => { let eddy = 0; for (let k = 0; k < L; k++) eddy += ov.eddyMass[c * L + k]; close(eddy, 0, 0.005, `eddies move no water across ${l}`); });
  assert.deepEqual(split.peak('direct', 1).latitude, split.peak('sea', 1).latitude);
  assert.deepEqual(split.peak('direct', -1).latitude, split.peak('sea', -1).latitude);
  const { groups } = classGroups(header);
  assert.deepEqual(groups.map((g) => g.name), GROUPS);
  for (const flow of ov.returnFlow) {
    const c = ov.caps.indexOf(flow.latitude);
    flow.groups.forEach((g, n) => { assert.equal(g.name, GROUPS[n]); close(groups[n].mass[c], g.mass, 0.02, `${g.name} at ${flow.latitude}`); close(groups[n].eddyMass[c], g.eddyMass, 0.02, `its eddies`); });
    close(groups.reduce((s, g) => s + g.heat[c] + g.eddyHeat[c], 0), flow.heat, 0.002, `heat at ${flow.latitude}`);
  }
});

test('the overturning by class closes at the floor, and its depths are the section’s', () => {
  const closed = classOverturning(header), open = classOverturning(header, null, { closed: false }), ov = header.overturning, M = L + 1;
  ov.caps.forEach((l, c) => {
    assert.equal(closed.psi[c * M], 0);
    let sum = 0;
    for (let k = 0; k < L; k++) sum += ov.mass[c * L + k] + ov.eddyMass[c * L + k];
    close(open.bottom[c], sum, 1e-9, `open at ${l}`); close(closed.bottom[c], sum - ov.imbalance[c], 1e-9, `closed at ${l}`);
    close(closed.bottom[c], ov.etaNet[c], 0.02, `closed to the sea level’s own transport at ${l}`);
    close(closed.depth[c * M + L], 4000, 2, `the floor at ${l}`);
    close(closed.depth[c * M + 1], 0.5 * (header.section.base[c * L] + header.section.base[(c + 1) * L]), 1e-9, `the mixed layer’s base at ${l}`);
  });
  for (const m of ['March', 'September']) { const flow = classOverturning(header, m), mo = month(m, 'overturning'); ov.caps.forEach((l, c) => close(flow.bottom[c], mo.etaNet[c], 0.02, `${m} at ${l}`)); }
  const section = classSection(header), sep = classSection(header, 'September');
  section.columns.forEach((col, b) => {
    assert.equal(col[0].class, 0); assert.equal(col[0].top, 0);
    for (let n = 1; n < col.length; n++) assert.ok(col[n].top >= col[n - 1].base - 1e-9, `classes stack at ${lat(b)}`);
    close(header.section.base[b * L + L - 1], 4000, 2, `the column reaches the floor at ${lat(b)}`);
  });
  assert.equal(sep.v, null); assert.ok(sep.columns.every((col) => col.every((c) => c.v === null)));
  const zm = zonalMeans(header), june = zonalMeans(header, 'June');
  assert.equal(zm.name, 'year'); assert.equal(june.samples, 120);
  zm.outcropped.forEach((o, b) => assert.equal(o, zm.thermoclineTop[b] - zm.mixedLayerDepth[b] < OUTCROP.zonal));
});

test('the prose numbers: the planet’s year and its books', () => {
  close(year.surfaceTemperature, 292.5, 0.05, 'surface');
  close(year.surfaceTemperature - part7Year.surfaceTemperature, 7.1, 0.05, 'warmer than Part 7');
  close(year.absorbedSolar, 229.3, 0.05, 'sunlight absorbed'); close(year.outgoingLongwave, 228.9, 0.05, 'infrared out');
  close(year.planetaryAlbedo, 0.327, 0.0005, 'reflected'); close(part7Year.planetaryAlbedo, 0.364, 0.0005, 'Part 7 reflected');
  close(1 - year.clearAbsorbedSolar / 340.5, 0.105, 0.0005, 'clear sky reflected'); close(1 - part7Year.clearAbsorbedSolar / 340.5, 0.136, 0.0005, 'Part 7’s');
  close(year.cloudSunlight, -75.4, 0.05, 'clouds in sunlight'); close(year.cloudInfrared, 33.0, 0.05, 'clouds in infrared');
  close(year.balance, 0.4, 0.05, 'gained at the top');
  close(year.seaSurfaceTemperature, 19.4, 0.05, 'the sea’s mean surface');
  const sst = header.zonalBudget.seaSurfaceTemperature.map((t) => t - 273.15);
  close(globalMean(sst, area), year.seaSurfaceTemperature, 0.01, 'the band mean is the daily mean');
  close(Math.max(...sst), 29.9, 0.05, 'warmest band'); close(sst[0], 3.3, 0.05, 'south pole'); close(sst[35], 3.3, 0.05, 'north pole');
  close(books.global.total, 1.95, 0.005, 'stored'); close(books.global.ocean, 1.87, 0.005, 'by the sea'); close(books.global.air, 0.08, 0.005, 'by the air');
  close(books.global.toaNet, 0.42, 0.005, 'from the top'); close(books.global.seaInput, 1.25, 0.005, 'into the sea');
  close(books.global.residual, 1.54, 0.005, 'unplaced'); close(books.global.surfaceResidual, 0.63, 0.005, 'on the sea’s side'); close(books.global.atmosphereResidual, 0.91, 0.005, 'on the air’s');
  const sigmaT4 = header.zonalBudget.surfaceTemperature.reduce((sum, t, b) => sum + 5.670374419e-8 * t ** 4 * area[b], 0) / sphere;
  close(1 - globalMean(header.zonalBudget.outgoingLongwave, area) / sigmaT4, 0.452, 0.0005, 'infrared held back');
  close(year.precipitableWater, 33.4, 0.05, 'vapor'); close(part7Year.precipitableWater, 23.0, 0.05, 'Part 7’s vapor');
});

test('the prose numbers: why it is warmer than Part 7', () => {
  const z8 = header.zonalBudget, z7 = part7.header.zonalBudget, sin = sinAreas(1).map((a) => a / (4 * Math.PI));
  const part = (values, lo, hi) => values.reduce((sum, v, b) => sum + (Math.abs(lat(b)) > lo && Math.abs(lat(b)) < hi ? v * sin[b] : 0), 0);
  const all = z8.absorbedSolar.map((v, b) => v - z7.absorbedSolar[b]), clear = z8.clearAbsorbedSolar.map((v, b) => v - z7.clearAbsorbedSolar[b]);
  close(part(all, 0, 90), 12.9, 0.05, 'more sunlight absorbed'); close(part(clear, 0, 90), 10.4, 0.05, 'under clear skies');
  close(part(clear, 50, 90), 8.3, 0.05, 'clear-sky poleward of 50°'); close(part(all, 0, 90) - part(clear, 0, 90), 2.5, 0.05, 'clouds');
  close(year.outgoingLongwave - part7Year.outgoingLongwave, 11.3, 0.05, 'more infrared out');
  close((year.outgoingLongwave - part7Year.outgoingLongwave) / (year.surfaceTemperature - part7Year.surfaceTemperature), 1.58, 0.005, 'W/m² per kelvin');
  const warmer = z8.surfaceTemperature.map((t, b) => t - z7.surfaceTemperature[b]);
  close(warmer[0], 22.0, 0.05, 'the south pole'); close(warmer[35], 10.5, 0.05, 'the north pole');
  close(Math.min(...warmer), 1.7, 0.05, 'least'); assert.equal(lat(warmer.indexOf(Math.min(...warmer))), 32.5); close(warmer[11], 6.2, 0.05, 'at 32.5° S');
  close(part7Year.ice, 11.9, 0.05, 'Part 7’s ice, percent of the globe');
  assert.equal(year.ice, 0);
});

test('the prose numbers: the drift through the spin-up and the ice that went', () => {
  const sst = yearlyMeans(header, 'seaSurfaceTemperature'), ts = yearlyMeans(header, 'surfaceTemperature'), toa = yearlyMeans(header, 'toaNet');
  assert.equal(sst.length, 11);
  close(sst[0], 17.08, 0.005, 'the first year'); close(sst[10], 19.38, 0.005, 'the recorded year');
  close((sst[10] - sst[8]) / 2, 0.12, 0.005, 'K a year over the last two years'); close(sst[9] - sst[8], 0.135, 0.0005, 'year 10'); close(sst[10] - sst[9], 0.096, 0.0005, 'the record');
  close(ts[10] - ts[9], 0.096, 0.0005, 'the surface too');
  close(toa[9], 0.40, 0.005, 'the top in year 10');
  const daily = header.daily, residualYear = (y) => (daily.surfaceResidual[y * YEAR - 1] * y - daily.surfaceResidual[(y - 1) * YEAR - 1] * (y - 1));
  for (let y = 3; y <= 11; y++) { const r = residualYear(y); assert.ok(r > 0.58 && r < 0.73, `the sea’s unplaced heat in year ${y}: ${r}`); }
  close(0.63 * YEAR * DAY / (header.constants.rhoCp * year.mixedLayerDepth), 0.072, 0.0005, 'K a year over the mean mixed layer');
  const lastDay = (name) => { let n = -1; daily[name].forEach((v, d) => { if (v > 0) n = d; }); return n + 1; };
  assert.equal(lastDay('iceNorth'), 423); assert.equal(date(daily.from / DAY + 423), 'May 17');
  assert.equal(lastDay('ice'), 966); assert.equal(date(daily.from / DAY + 966), 'November 11');
  close(daily.ice[0], 6.58, 0.005, 'ice at the start');
  close(daily.ice[59], 4.96, 0.005, 'day 60'); close(part7.header.daily.ice[59], 8.95, 0.005, 'Part 7’s day 60');
  const peaks = [1, 2].map((y) => Math.max(...daily.iceSouth.slice(y * YEAR, y * YEAR + YEAR)));
  close(peaks[0], 2.91, 0.005, 'southern peak in year 2'); close(peaks[1], 0.36, 0.005, 'in year 3');
});

test('the prose numbers: the sea’s surface through the seasons', () => {
  const zm = zonalMeans(header), eq = (v) => (v[17] + v[18]) / 2;
  close(eq(zm.mixedLayerDepth), 51, 0.5, 'equatorial mixed layer'); close(Math.max(...zm.mixedLayerDepth), 82, 0.5, 'deepest band in the year’s mean');
  close(zm.thermoclineTop[15], 290.7, 0.05, '1024 class at 12.5° S'); close(zm.thermoclineTop[21], 291.4, 0.05, 'at 17.5° N'); close(eq(zm.thermoclineTop), 239, 0.5, 'shallower over the equator');
  assert.deepEqual(zm.outcropped.map((o, b) => (o ? lat(b) : null)).filter((v) => v !== null && Math.abs(v) < 40), [-37.5, 37.5]);
  close(zm.seaSurfaceHeight[14] * 100, 86, 0.5, 'hump at 17.5° S'); close(zm.seaSurfaceHeight[21] * 100, 83, 0.5, 'at 17.5° N'); close(eq(zm.seaSurfaceHeight) * 100, 50, 0.5, 'dip at the equator');
  close(zm.seaSurfaceHeight[7] * 100, -111, 0.5, 'low at 52.5° S'); close(zm.seaSurfaceHeight[28] * 100, -109, 0.5, 'at 52.5° N');
  close(zm.uSurface[18] * 100, -31.8, 0.05, 'westward at 2.5° N'); close(zm.uSurface[10] * 100, 10.2, 0.05, 'eastward at 37.5° S'); close(zm.uSurface[24] * 100, 10.1, 0.05, 'at 32.5° N');
  close(zm.upwelling[17], 0.15, 0.005, 'rising at 2.5° S'); close(Math.min(...zm.upwelling.slice(12, 24)), -0.09, 0.005, 'sinking in the subtropics');
  const deepest = (m) => { const v = month(m, 'zonalMean').mixedLayerDepth, d = Math.max(...v); return [Math.round(d), lat(v.indexOf(d))]; };
  assert.deepEqual(deepest('March'), [156, 77.5]); assert.deepEqual(deepest('July'), [139, -37.5]); assert.deepEqual(deepest('September'), [127, -77.5]);
  const times = run.times('year'), k = field(header, 'seaSurfaceTemperature'), at45 = [], north = [];
  for (let n = 0; n + 1 < times.length; n++) { const v = run.values('year', n, k); at45.push(band(v, 42.5, 47.5)); north.push(band(v, 80, 90.1)); }
  const swing = (s) => Math.max(...s) - Math.min(...s);
  close(swing(at45), 3.8, 0.05, 'the sea’s swing at 45° N'); assert.equal(date(times[at45.indexOf(Math.max(...at45))] / DAY), 'August 27');
  close(swing(north), 5.0, 0.05, 'beyond 80° N');
});

test('the prose numbers: the rain band stays near the equator', () => {
  const wander = (r) => {
    const times = r.times('year'), k = field(r.header, 'rain'), out = [];
    for (let n = 1; n < times.length; n++) out.push({ day: times[n] / DAY - 2.5, ...peakLatitude(zonal(r.values('year', n, k)), 12, 24) });
    return [out.reduce((a, b) => (b.latitude > a.latitude ? b : a)), out.reduce((a, b) => (b.latitude < a.latitude ? b : a))];
  };
  const [north, south] = wander(run), [north7, south7] = wander(part7);
  close(north.latitude, 12.6, 0.05, 'farthest north'); assert.equal(date(north.day), 'September 19');
  close(south.latitude, -14.2, 0.05, 'farthest south'); assert.equal(date(south.day), 'June 6');
  close(north7.latitude, 23.3, 0.05, 'Part 7 north'); close(south7.latitude, -18.4, 0.05, 'Part 7 south');
  const monthly = header.zonalBudget.months.map((m) => peakLatitude(m.precipitation, 12, 24).latitude);
  close(Math.min(...monthly), -8.4, 0.05, 'monthly south'); close(Math.max(...monthly), 9.5, 0.05, 'monthly north');
  const annual = peakLatitude(header.zonalBudget.precipitation);
  close(annual.latitude, -3.5, 0.05, 'the year’s band'); close(annual.value, 5.0, 0.05, 'its rain');
  close(peakLatitude(header.zonalBudget.precipitation, 24, 36).latitude, 41.0, 0.05, 'northern storm track'); close(peakLatitude(header.zonalBudget.precipitation, 0, 12).latitude, -41.5, 0.05, 'southern');
});

test('the prose numbers: the currents and the overturning', () => {
  const zm = zonalMeans(header), runs = [];
  zm.zonalTransport.forEach((v, b) => { const s = Math.sign(v); if (runs.length && runs.at(-1).sign === s) runs.at(-1).sum += v; else runs.push({ sign: s, from: lat(b) - 2.5, sum: v }); });
  assert.deepEqual(runs.map((r) => [r.from, Math.round(r.sum)]), [[-90, -186], [-50, 196], [-25, -831], [25, 175], [45, -206]]);
  close(zm.uDepthMean[17] * 100, -5.8, 0.05, 'depth-mean at 2.5° S');
  const flow = classOverturning(header), M = L + 1, under = (l) => flow.psi[flow.caps.indexOf(l) * M + 1];
  close(under(-10), -58.4, 0.05, 'south across 10° S'); close(under(10), 52.4, 0.05, 'north across 10° N');
  close(under(-35), 35.1, 0.05, 'north across 35° S'); close(under(35), -34.6, 0.05, 'south across 35° N');
  const deep = (sign) => { let best = { psi: 0 }; flow.caps.forEach((l, c) => { for (let k = 0; k <= L; k++) if (flow.depth[c * M + k] > 500 && sign * flow.psi[c * M + k] > sign * best.psi) best = { psi: flow.psi[c * M + k], l, depth: flow.depth[c * M + k] }; }); return best; };
  const north = deep(1), south = deep(-1);
  assert.deepEqual([north.l, south.l], [55, -55]); close(north.psi, 32.7, 0.05, 'the northern deep cell'); close(south.psi, -26.2, 0.05, 'the southern');
  close(classOverturning(header, 'September').psi[flow.caps.indexOf(-10) * M + 3], -90.3, 0.05, 'the southern cell in September');
  const rf = (l) => header.overturning.returnFlow.find((r) => r.latitude === l);
  close(rf(-10).groups[0].temperature, 30.0, 0.05, 'mixed layer going south'); close(rf(10).groups[0].temperature, 30.1, 0.05, 'going north');
  const { groups } = classGroups(header), at = (l) => header.overturning.caps.indexOf(l), back = (l) => (groups[1].heat[at(l)] + groups[1].eddyHeat[at(l)]) / ((groups[1].mass[at(l)] + groups[1].eddyMass[at(l)]) * header.constants.rhoCp * 1e-9);
  close(back(-10), 24.7, 0.05, 'the thermocline’s return at 10° S'); close(back(10), 26.2, 0.05, 'at 10° N');
});

test('the prose numbers: the heat carried, split between the air and the sea', () => {
  const s = split, p = (which, sign) => s.peak(which, sign);
  assert.deepEqual([p('total', 1).latitude, p('total', -1).latitude], [50, -50]);
  close(p('total', 1).transport, 6.1, 0.05, 'in all, north'); close(p('total', -1).transport, -5.7, 0.05, 'south');
  assert.deepEqual([p('air', 1).latitude, p('air', -1).latitude], [40, -35]);
  close(p('air', 1).transport, 4.4, 0.05, 'the air north'); close(p('air', -1).transport, -4.5, 0.05, 'south');
  close(s.at(40).air, 4.4, 0.05, 'the air at 40° N'); close(s.at(-40).air, -4.4, 0.05, 'at 40° S');
  close(s.at(40).sea, 1.0, 0.05, 'the sea at 40° N'); close(s.at(-40).sea, -1.1, 0.05, 'at 40° S');
  close(s.at(40).total, 5.4, 0.05, 'in all at 40° N'); close(s.at(-40).total, -5.6, 0.05, 'at 40° S');
  close(p('sea', 1).transport, 2.2, 0.05, 'the sea north, at 50°'); close(p('sea', -1).transport, -1.8, 0.05, 'south');
  close(s.at(-10).sea, -1.4, 0.05, 'the sea at 10° S'); close(s.at(10).sea, 0.5, 0.05, 'at 10° N'); close(s.at(20).sea, 0.7, 0.05, 'at 20° N');
  close(s.at(0).total, -0.4, 0.05, 'in all across the equator'); close(s.at(0).sea, -0.6, 0.05, 'the sea'); close(s.at(0).air, 0.2, 0.05, 'the air');
  close(s.at(0).latent, -0.2, 0.05, 'vapor across the equator'); close(s.at(0).dry, 0.3, 0.05, 'dry air');
  close(p('latent', 1).transport, 3.0, 0.05, 'vapor north'); close(p('latent', -1).transport, -3.1, 0.05, 'vapor south');
  assert.deepEqual([p('latent', 1).latitude, p('latent', -1).latitude], [35, -35]);
  assert.deepEqual(s.zeros.latent.map((l) => +l.toFixed(1)), [-15.8, -1.0, 14.4]);
  close(p('bareTotal', 1).transport, 5.1, 0.05, 'from the top of the air alone, north'); close(p('bareTotal', -1).transport, -5.1, 0.05, 'south');
  close(s.at(0).bareTotal, 0.0, 0.05, 'and across the equator');
  close(p('bareSea', 1).transport, 1.0, 0.05, 'the sea without its storage, north'); close(p('bareSea', -1).transport, -1.1, 0.05, 'south');
  assert.deepEqual([p('bareSea', 1).latitude, p('bareSea', -1).latitude], [10, -10]);
  close(s.at(50).sea - s.at(50).bareSea, 1.87, 0.005, 'what the sea stores poleward of 50° N');
  const cap = (from, to) => { let sum = 0, a = 0; for (let b = from; b < to; b++) { sum += books.stores.ocean[b] * area[b]; a += area[b]; } return sum / a; };
  close(cap(29, 36), 41.5, 0.05, 'the northern polar sea warms, W/m²'); close(cap(0, 7), 37.5, 0.05, 'the southern'); close(cap(26, 28), -45.7, 0.05, 'while 40–50° N cools');
  close(books.global.residual * sphere / 2 / 1e15, 0.39, 0.005, 'the residual put in one hemisphere moves the equator by');
});

test('the prose numbers: the water cycle', () => {
  const g = water.global;
  close(g.evaporation, 3.08, 0.005, 'evaporation'); close(g.precipitation, 3.08, 0.005, 'rain');
  close(g.evaporation / part7Year.evaporation - 1, 0.16, 0.005, 'more than Part 7’s 2.66');
  close(g.latentHeat, 89.1, 0.05, 'latent heat'); close(g.sensibleHeat, 11.2, 0.05, 'sensible heat'); close(g.latentHeat / g.sensibleHeat, 7.9, 0.05, 'times');
  close(g.convectiveRain / g.precipitation, 0.52, 0.005, 'from the plume');
  close(year.residenceDays, 10.9, 0.05, 'days of rain in the air');
  assert.deepEqual(water.gapCrossings.map((l) => +l.toFixed(0)), [-34, -8, 8, 33]);
  close(Math.min(...water.vaporSverdrups), -1.25, 0.005, 'fresh water carried south by the air, Sv'); close(Math.max(...water.vaporSverdrups), 1.19, 0.005, 'north');
});

test('the prose numbers: the spin-up’s first current, the polar sea that stays open and the deep cells that drain', () => {
  const st = run.times('spinup'), east = field(header, 'currentEast');
  const equator = (day) => band(run.values('spinup', st.findIndex((t) => Math.round(t / DAY + header.spinupDays) === day), east), -5, 5) * 100;
  close(equator(10), 4.6, 0.05, 'eastward on day 10'); assert.ok(equator(20) < 0, 'westward by day 20'); close(equator(60), -18.5, 0.05, 'westward on day 60');
  const lastDay = (name) => { let n = -1; header.daily[name].forEach((v, d) => { if (v > 0) n = d; }); return n + 1; };
  assert.equal(lastDay('iceSouth'), 965); assert.equal(date(header.daily.from / DAY + 965), 'November 10');
  const times = run.times('year'), sst = field(header, 'seaSurfaceTemperature'), sss = field(header, 'surfaceSalinity');
  let coldest = Infinity, densest = 0;
  times.forEach((_, k) => { const t = run.values('year', k, sst), s = run.values('year', k, sss); t.forEach((v, i) => { coldest = Math.min(coldest, v); densest = Math.max(densest, seawaterDensity(v + 273.15, s[i])); }); });
  close(coldest, 0.86, 0.005, 'the coldest cell, °C'); close(coldest + 1.8, 2.7, 0.05, 'above freezing');
  close(densest, 1026.02, 0.005, 'the densest surface water');
  const cap = (values, from, to) => { let sum = 0, a = 0; for (let b = from; b < to; b++) { sum += values[b] * area[b]; a += area[b]; } return sum / a; };
  close(cap(books.sea, 29, 36) - cap(books.seaInput, 29, 36), 46.1, 0.05, 'the currents bring the northern polar sea, W/m²');
  close(cap(books.sea, 0, 7) - cap(books.seaInput, 0, 7), 39.9, 0.05, 'the southern');
  const k = header.classes.indexOf(1026.95), months = header.section.months, Sv = 1e6;
  const volume = (m, from, to) => { let v = 0; for (let b = from; b < to; b++) v += (m.base[b * L + k] - m.base[b * L + k - 1]) * area[b]; return v; };
  months.forEach((m) => close(volume(m, 0, 36) / volume(months[0], 0, 36), 1, 1e-4, `the class’s whole volume in ${m.name}`));
  const april = months.find((m) => m.name === 'April'), february = months.find((m) => m.name === 'February');
  const seconds = ((february.ranges[0][0] + february.ranges[0][1]) / 2 - (april.ranges[0][0] + april.ranges[0][1]) / 2) * DAY;
  const across = (l) => { const c = header.overturning.caps.indexOf(l); return header.overturning.mass[c * L + k] + header.overturning.eddyMass[c * L + k]; };
  const north = (volume(april, 29, 36) - volume(february, 29, 36)) / seconds / Sv, south = (volume(april, 0, 7) - volume(february, 0, 7)) / seconds / Sv;
  close(north, 32.5, 0.05, 'the north loses, Sv'); close(across(55), -32.8, 0.05, 'and as much flows south across 55° N');
  close(south, 25.9, 0.05, 'the south loses, Sv'); close(across(-55), 26.2, 0.05, 'and as much flows north across 55° S');
  const yearly = (from, to) => (volume(april, from, to) - volume(february, from, to)) / volume(april, from, to) * YEAR * DAY / seconds;
  close(yearly(29, 36), 0.022, 0.0005, 'a year’s loss in the north'); close(yearly(0, 7), 0.017, 0.0005, 'in the south');
});

test('the prose numbers: the cells of the section, its salt, and the September return', () => {
  const flow = classOverturning(header), M = L + 1, under = (f, l) => f.psi[f.caps.indexOf(l) * M + 1];
  close(under(flow, -55), -10.8, 0.05, 'the southern polar cell'); close(under(flow, 55), 12.2, 0.05, 'the northern');
  const sep = classOverturning(header, 'September'), c = sep.caps.indexOf(0);
  const psi = (z) => { let j = 0; while (j < L - 1 && sep.depth[c * M + j + 1] < z) j++; const d0 = sep.depth[c * M + j], d1 = sep.depth[c * M + j + 1], t = Math.min(1, Math.max(0, (z - d0) / (d1 - d0))); return sep.psi[c * M + j] + t * (sep.psi[c * M + j + 1] - sep.psi[c * M + j]); };
  close(psi(100), -71.3, 0.05, 'south above 100 m'); close(psi(400), -62.1, 0.05, 'above 400 m'); close(psi(1000), -32.8, 0.05, 'above 1,000 m'); close(psi(3000), -17.8, 0.05, 'above 3,000 m');
  close(sep.imbalance[c], 0.5, 0.05, 'what the samples miss there');
  const section = classSection(header), column = (l) => section.columns[section.centers.indexOf(l)];
  for (const l of [-22.5, 22.5]) {
    const salty = column(l).filter((x) => x.S > 35);
    close(Math.max(...salty.map((x) => x.base)), 600, 10, `saltier than 35 psu down to 600 m at ${l}`);
    assert.ok(salty.every((x) => x.top < 610), `and no deeper at ${l}`);
  }
  for (const b of [0, 35]) {
    const col = section.columns[b];
    close(col[0].S, 33.6, 0.05, `the polar mixed layer at ${lat(b)}`); close(col.at(-1).S, 34.65, 0.06, `over the deep water at ${lat(b)}`);
  }
});

test('the prose numbers: the sea’s share without its storage, against Part 7 and the Earth', () => {
  const s = split, at = (l) => s.at(l);
  close(at(40).sea - at(40).bareSea, 0.58, 0.005, 'storage in the sea’s share at 40° N'); close(at(-40).sea - at(-40).bareSea, -0.71, 0.005, 'at 40° S');
  close(at(40).bareSea, 0.41, 0.005, 'the sea without its storage at 40° N'); close(at(-40).bareSea, -0.40, 0.005, 'at 40° S');
  assert.ok(at(40).bareSea / at(40).bareTotal < 0.1 && at(-40).bareSea / at(-40).bareTotal < 0.1, 'under a tenth');
  close(at(40).bareTotal - at(40).bareSea, at(40).air, 0.02, 'the air just as much at 40° N'); close(at(-40).bareTotal - at(-40).bareSea, at(-40).air, 0.02, 'at 40° S');
  close(at(35).bareSea / at(35).bareTotal, 0.135, 0.0005, 'the sea’s share across 35° N'); close(at(-35).bareSea / at(-35).bareTotal, 0.114, 0.0005, 'across 35° S');
  close(at(10).sea, 0.46, 0.005, 'the sea across 10° N');
  const z7 = part7.header.zonalBudget, top7 = integrate(z7.absorbedSolar.map((v, b) => v - z7.outgoingLongwave[b]), area).transport, e = s.edges;
  close(top7[e.indexOf(35)], 4.6, 0.05, 'Part 7 from the top alone at 35° N'); close(Math.min(...top7), -4.7, 0.05, 'and its southern peak');
  close(at(35).bareTotal - top7[e.indexOf(35)], 0.45, 0.005, 'more across 35° N'); close(at(-35).bareTotal - top7[e.indexOf(-35)], -0.41, 0.005, 'across 35° S');
  close(at(35).bareTotal / top7[e.indexOf(35)] - 1, 0.097, 0.0005, 'a tenth more'); close(at(-35).bareTotal / top7[e.indexOf(-35)] - 1, 0.088, 0.0005, 'nearly a tenth');
  const z8 = header.zonalBudget, change = (name) => z8[name].map((v, b) => v - z7[name][b]), cap = (values, from, to) => { let sum = 0, a = 0; for (let b = from; b < to; b++) { sum += values[b] * area[b]; a += area[b]; } return sum / a; };
  close(cap(change('absorbedSolar'), 0, 36), 12.9, 0.05, 'more sunlight'); close(cap(change('absorbedSolar'), 25, 36), 4.9, 0.05, 'poleward of 35° N');
  close(cap(change('outgoingLongwave'), 0, 36), 11.3, 0.05, 'more infrared'); close(cap(change('outgoingLongwave'), 0, 11), 17.3, 0.05, 'poleward of 35° S');
  close(cap(change('surfaceTemperature'), 0, 11), 12.8, 0.05, 'warmer poleward of 35° S');
  const anomaly = (z, name, from, to) => { const mean = globalMean(z[name], area); let sum = 0; for (let b = from; b < to; b++) sum += (z[name][b] - mean) * area[b] / 1e15; return sum; };
  const extra = (name, from, to) => anomaly(z8, name, from, to) - anomaly(z7, name, from, to);
  close(-extra('absorbedSolar', 25, 36), 0.885, 0.005, 'the sunlight adds across 35° N, PW'); close(extra('outgoingLongwave', 25, 36), -0.438, 0.005, 'the infrared takes back');
  close(extra('absorbedSolar', 0, 11), 0.253, 0.005, 'the sunlight alone would carry less south across 35° S'); close(-extra('outgoingLongwave', 0, 11), -0.662, 0.005, 'the infrared carries more south');
  const ov = header.overturning, sec = header.section, rc = header.constants.rhoCp;
  const parts = (l) => {
    const c = ov.caps.indexOf(l);
    let total = 0, mean = 0;
    for (let k = 0; k < L; k++) {
      const i = c * L + k, m = ov.mass[i] + ov.eddyMass[i], h = (b) => sec.base[b * L + k] - (k ? sec.base[b * L + k - 1] : 0), T0 = sec.T[c * L + k], T1 = sec.T[(c + 1) * L + k];
      const T = T0 === null ? T1 ?? 0 : T1 === null ? T0 : (T0 * h(c) + T1 * h(c + 1)) / (h(c) + h(c + 1));
      total += ov.heat[i] + ov.eddyHeat[i]; mean += rc * m * 1e6 * T / 1e15;
    }
    return { total, cells: mean, within: total - mean };
  };
  for (const l of [-40, -35, -30, 30, 35, 40]) {
    const p = parts(l);
    assert.ok(Math.sign(l) * p.cells < 0.05, `the cells carry next to none of it poleward across ${l}: ${p.cells}`);
    assert.ok(p.within / p.total > 0.93, `nearly all of it moves within classes across ${l}: ${p.within} of ${p.total}`);
  }
  close(parts(50).within, 1.47, 0.005, 'within classes across 50° N'); close(parts(50).total, 2.14, 0.005, 'of');
  close(parts(-50).within, -1.33, 0.005, 'across 50° S'); close(parts(-50).total, -1.76, 0.005, 'of');
  const tropics = parts(-10);
  close(tropics.cells / tropics.total, 0.83, 0.005, 'across 10° S the cells carry most of it');
});
