import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { readFrames } from '../explain/frames.module.js';
import { heatTransport, storedHeat } from '../explain/sunlitRun.module.js';
import { calendarDate, storage, recordedYear, waterTransport, overturning, peakLatitude, PART6_TRANSPORT } from '../explain/wetRun.module.js';

const open = (name) => { const file = new URL(`../explain/data/${name}`, import.meta.url), bytes = readFileSync(file); return readFrames(file.href, bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength)); };
const run = open('water_N16.bin'), part6 = open('radiation_N16.bin');
const { header } = run, DAY = 86400, YEAR = 365, K = header.K, D = Math.PI / 180;
const cells = [...new Grid(header.N)].map((c) => ({ lat: Math.asin(c.centerVertex.z) / D, area: c.area }));
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const field = (name) => header.fields.findIndex((f) => f.name === name);
const lat = (b) => -90 + (b + 0.5) * 5;
const bandArea = (b) => (Math.sin((lat(b) + 2.5) * D) - Math.sin((lat(b) - 2.5) * D)) / 2;
const hemisphere = (values, north) => values.reduce((sum, v, b) => sum + ((b >= 18) === north ? 2 * v * bandArea(b) : 0), 0);
const month = (name, part = 'zonalBudget') => header[part].months.find((m) => m.name === name);
const band = (values, lo, hi) => { let sum = 0, area = 0; cells.forEach((c, i) => { if (c.lat >= lo && c.lat < hi) { sum += values[i] * c.area; area += c.area; } }); return sum / area; };
const zonal = (values) => { const sum = new Float64Array(36), area = new Float64Array(36); cells.forEach((c, i) => { const b = Math.min(35, Math.floor((c.lat + 90) / 5)); sum[b] += values[i] * c.area; area[b] += c.area; }); return Array.from(sum, (v, b) => v / area[b]); };
const date = (day) => { const c = calendarDate(day); return `${c.month} ${c.date}`; };
const year = recordedYear(header), water = waterTransport(header), kept = storage(header);
const sigmaNear = (s) => header.zonalMean.sigma.reduce((best, v, k) => (Math.abs(v - s) < Math.abs(header.zonalMean.sigma[best] - s) ? k : best), 0);

test('the run is the climate model at N = 16 over a slab sea, a year of spin-up and a recorded year from the March equinox', () => {
  assert.deepEqual([header.case, header.N, K, header.dt, header.grid, header.spinupDays, header.recordDays], ['wet', 16, 36, 450, 'bl36', 365, 365]);
  assert.equal(header.options.ocean, false);
  assert.deepEqual(header.ice.options, { iceAlbedo: 0.07, meltingIceAlbedo: 0.07, iceSnowAlbedo: 0.07, slabHeatCapacity: 4.2e7 });
  assert.deepEqual(header.ice.diffuseAlbedo, { water: 0.06, coldIce: 0.07, meltingIce: 0.07, freshSnowOnIce: 0.85 });
  close(header.ice.slabHeatCapacity / 4.2e7, 1, 1e-12, 'ten meters of water');
  assert.equal(run.times('year').length, 74);
  assert.equal(run.times('weather').length, 121);
  assert.equal(date(run.times('weather')[0] / DAY), 'January 1');
  assert.equal(date(run.times('weather').at(-1) / DAY), 'January 31');
  assert.equal(header.daily.surfaceTemperature.length, header.spinupDays + header.recordDays);
});

test('the stored heat and water change by what the stored states say', () => {
  const { start, end } = header.storage;
  assert.equal(kept.seconds, YEAR * DAY);
  kept.surface.forEach((v, b) => close(v * kept.seconds, end.surfaceEnergy[b] - start.surfaceEnergy[b], 1e-3, `surface ${b}`));
  kept.air.forEach((v, b) => close(v * kept.seconds, end.airEnergy[b] - start.airEnergy[b], 1e-3, `air ${b}`));
  kept.vapor.forEach((v, b) => close(v * kept.seconds / header.constants.latentHeat, end.precipitableWater[b] - start.precipitableWater[b], 1e-9, `vapor ${b}`));
  kept.total.forEach((v, b) => close(v, kept.surface[b] + kept.air[b] + kept.vapor[b], 1e-12, `total ${b}`));
  const global = (values) => values.reduce((sum, v, b) => sum + v * bandArea(b), 0);
  close(global(kept.surface), -0.338, 0.0005, 'the sea and its ice');
  close(global(kept.air), 0.050, 0.0005, 'the air');
  close(global(kept.vapor), 0.079, 0.0005, 'the vapor');
  close(global(kept.total), -0.209, 0.0005, 'all of it');
});

test('the heat carried is zero at both poles, and its latent and dry parts add up to it', () => {
  const { total, latent, dry } = water;
  for (const t of [total, latent]) { assert.equal(t.transport[0], 0); close(t.transport.at(-1), 0, 1e-12, 'north pole'); }
  dry.forEach((v, b) => close(v + latent.transport[b], total.transport[b], 1e-12, `split at ${water.edges[b]}`));
  close(total.imbalance, total.absorbedMean - total.outgoingMean - total.storedMean, 1e-12, 'imbalance');
  close(total.imbalance, -0.958, 0.0005, 'what the run leaves unaccounted, spread evenly');
  assert.ok(Math.abs(water.waterImbalance) < 0.002, `water is kept: ${water.waterImbalance} mm a day`);
  for (const l of [-60, -20, 0, 35]) { const at = water.at(l); close(at.total - at.latent, at.dry, 1e-9, `at ${l}`); }
});

test('Part 6’s dry transport is its own run’s', () => {
  const t = heatTransport(part6.header.zonalBudget, storedHeat(part6, cells).total);
  t.transport.forEach((v, b) => close(PART6_TRANSPORT[b], v, 0.0005, `edge ${b}`));
});

test('the overturning closes at the ground and turns the way the rain band says', () => {
  const flow = overturning(header);
  for (let b = 0; b < 36; b++) { assert.equal(flow.psi[b * (K + 1)], 0); close(flow.psi[b * (K + 1) + K] / 1e9, 0, 1e-3, `ground at ${lat(b)}`); }
  const north = flow.cell(1, -30, 40), south = flow.cell(-1, -40, 30);
  assert.ok(north.latitude > 0 && south.latitude < 0, 'a cell on each side');
  const k500 = sigmaNear(0.5);
  assert.ok(flow.omega[18 * (K + 1) + k500] < 0, 'rising near the equator');
  assert.ok(flow.omega[13 * (K + 1) + k500] > 0, 'sinking at 22.5° S');
  const september = overturning(header, 'September'), across = september.cell(-1, -40, 40);
  assert.ok(Math.abs(across.latitude) < 5 && Math.abs(across.psi) > 4 * Math.abs(north.psi), 'the winter cell crosses the equator in September');
});

test('the prose numbers: the planet’s year', () => {
  close(year.surfaceTemperature, 285.4, 0.05, 'surface');
  close(year.absorbedSolar, 216.4, 0.05, 'sunlight absorbed');
  close(year.outgoingLongwave, 217.6, 0.05, 'infrared out');
  close(year.planetaryAlbedo, 0.364, 0.0005, 'reflected');
  close(1 - year.clearAbsorbedSolar / 340.5, 0.136, 0.0005, 'clear sky reflected');
  close(year.cloudSunlight, -77.9, 0.05, 'clouds in sunlight');
  close(year.cloudInfrared, 31.4, 0.05, 'clouds in infrared');
  close(year.surfaceTemperature - (year.absorbedSolar / 5.670374419e-8) ** 0.25, 36.8, 0.05, 'above the effective temperature');
  const from = -part6.header.daily.from / DAY, dry = part6.header.daily.surfaceTemperature.slice(from, from + YEAR).reduce((s, v) => s + v, 0) / YEAR;
  close(year.surfaceTemperature - dry, 6.1, 0.05, 'warmer than Part 6’s sphere');
  close(year.balance, -1.2, 0.05, 'lost at the top');
  close(year.balance - kept.total.reduce((sum, v, b) => sum + v * bandArea(b), 0), -1.0, 0.05, 'the gap');
  close(-kept.total.reduce((sum, v, b) => sum + v * bandArea(b), 0), 0.2, 0.05, 'what the sea, ice and air gave up');
  const sigmaT4 = header.zonalBudget.surfaceTemperature.reduce((sum, t, b) => sum + 5.670374419e-8 * t ** 4 * bandArea(b), 0);
  const olr = header.zonalBudget.outgoingLongwave.reduce((sum, v, b) => sum + v * bandArea(b), 0);
  assert.ok(1 - olr / sigmaT4 > 0.4 && 1 - olr / sigmaT4 < 0.45, `held back ${1 - olr / sigmaT4}`);
});

test('the prose numbers: the drift over the recorded year', () => {
  const { start, end } = header.storage;
  assert.equal(start.time, 0); assert.equal(end.time, YEAR * DAY);
  close(hemisphere(end.surfaceTemperature, true) - hemisphere(start.surfaceTemperature, true), -2.1, 0.05, 'the north cooled');
  close(hemisphere(end.surfaceTemperature, false) - hemisphere(start.surfaceTemperature, false), 1.1, 0.05, 'the south warmed');
  const growth = hemisphere(end.iceVolume, true) / hemisphere(start.iceVolume, true) - 1;
  assert.ok(growth > 0.4 && growth < 0.5, `northern ice grew by nearly half: ${growth}`);
});

test('the prose numbers: the water cycle', () => {
  close(year.evaporation, 2.66, 0.005, 'evaporation');
  close(year.precipitation, 2.66, 0.005, 'rain');
  close(year.evaporation * header.constants.latentHeat / DAY, 77, 0.5, 'latent heat');
  close(year.evaporation * header.constants.latentHeat / DAY / year.sensibleHeat, 5.5, 0.05, 'against the sensible heat');
  close(year.sensibleHeat, 14.0, 0.05, 'sensible heat');
  close(year.convectiveRain / year.precipitation, 0.53, 0.005, 'from the plume');
  close(year.precipitableWater, 23.0, 0.05, 'vapor');
  close(year.residenceDays, 8.7, 0.05, 'days of rain');
  const { evaporation, precipitation } = header.zonalBudget, tropics = evaporation.slice(13, 23);
  close(Math.min(...tropics), 3.7, 0.05, 'slowest between 25° S and 25° N'); close(Math.max(...tropics), 4.0, 0.05, 'fastest');
  const peak = peakLatitude(precipitation), north = peakLatitude(precipitation, 24, 36), south = peakLatitude(precipitation, 0, 12);
  close(peak.latitude, 15.7, 0.05, 'rain band'); close(peak.value, 4.2, 0.05, 'its rain');
  close(north.latitude, 43.0, 0.05, 'northern storm track'); close(north.value, 3.0, 0.05, 'its rain');
  close(south.latitude, -41.6, 0.05, 'southern storm track'); close(south.value, 2.5, 0.05, 'its rain');
  assert.deepEqual(water.gapCrossings.map((l) => +l.toFixed(0)), [-37, 6, 22, 37]);
  const sources = [peakLatitude(water.gap, 0, 18), peakLatitude(water.gap, 18, 36)];
  close(sources[0].latitude, -25.5, 0.05, 'southern source'); close(sources[0].value, 1.7, 0.05, 'its excess');
  close(sources[1].latitude, 29.4, 0.05, 'northern source'); close(sources[1].value, 1.1, 0.05, 'its excess');
});

test('the prose numbers: the heat carried', () => {
  assert.deepEqual([water.peak('total', 1).latitude, water.peak('total', -1).latitude], [40, -40]);
  close(water.peak('total', 1).transport, 4.2, 0.05, 'north'); close(water.peak('total', -1).transport, -5.0, 0.05, 'south');
  close(water.peak('part6', 1).transport, 1.1, 0.005, 'Part 6');
  close(water.at(0).total, -0.6, 0.05, 'across the equator');
  close(water.at(35).latent, 2.1, 0.05, 'vapor at 35° N'); close(water.at(-35).latent, -2.1, 0.05, 'vapor at 35° S');
  close(water.at(35).latent / water.at(35).total, 0.49, 0.005, 'half'); close(water.at(-35).latent / water.at(-35).total, 0.42, 0.005, 'two-fifths');
  assert.equal(water.zeros.latent.length, 1); close(water.zeros.latent[0], -21.2, 0.05, 'the vapor turns north');
  for (let l = -21; l <= 90; l++) assert.ok(water.at(l).latent > 0 || l === 90, `vapor carried north at ${l}`);
  close(water.zeros.dry[0], 10.8, 0.05, 'the dry part turns north');
  close(peakLatitude(header.zonalBudget.precipitation).latitude, 15.7, 0.05, 'the band on the year’s average');
  close(water.at(0).latent, 1.3, 0.05, 'vapor across the equator'); close(water.at(0).dry, -1.9, 0.05, 'dry air across the equator');
  close(water.at(-20).dry, -4.2, 0.05, 'dry air at 20° S'); assert.equal(water.peak('dry', -1).latitude, -20);
  close(water.at(-15).total, -3.3, 0.05, 'all at 15° S'); close(water.at(-15).dry, -3.9, 0.05, 'dry air at 15° S');
  const bare = heatTransport(header.zonalBudget, null, header.constants.radius);
  close(bare.at(0).transport, -0.1, 0.05, 'from the top of the air alone'); close(bare.at(40).transport, 4.6, 0.05, 'north'); close(bare.at(-40).transport, -4.7, 0.05, 'south');
  const sphere = 4 * Math.PI * header.constants.radius ** 2, net = header.zonalBudget.absorbedSolar.map((v, b) => v - header.zonalBudget.outgoingLongwave[b] - kept.total[b]);
  const crossing = (north) => net.slice(0, 18).reduce((sum, v, b) => sum + (v - (north ? 0 : 2 * water.imbalance)) * bandArea(b) * sphere / 1e15, 0);
  close(crossing(true), -0.8, 0.05, 'the gap all in the north'); close(crossing(false), -0.3, 0.05, 'all in the south');
  const drop = (values) => (values[17] + values[18]) / 2 - (values[0] + values[35]) / 2, dry = part6.header.zonalBudget;
  close(drop(header.zonalBudget.absorbedSolar), 189, 0.5, 'sunlight absorbed'); close(drop(header.zonalBudget.outgoingLongwave), 47, 0.5, 'infrared');
  close(drop(dry.absorbedSolar), 218, 0.5, 'Part 6’s sunlight'); close(drop(dry.outgoingLongwave), 179, 0.5, 'Part 6’s infrared');
});

test('the prose numbers: the rain band, the dry subtropics and the vapor', () => {
  const times = run.times('year'), rain = field('rain'), band5 = [];
  for (let n = 1; n < times.length; n++) band5.push({ day: times[n] / DAY - 2.5, ...peakLatitude(zonal(run.values('year', n, rain)), 6, 30) });
  const north = band5.reduce((a, b) => (b.latitude > a.latitude ? b : a)), south = band5.reduce((a, b) => (b.latitude < a.latitude ? b : a));
  close(north.latitude, 23.3, 0.05, 'farthest north'); assert.equal(date(north.day), 'September 14');
  close(south.latitude, -18.4, 0.05, 'back south'); assert.equal(date(south.day), 'March 18');
  const opening = peakLatitude(zonal(run.values('year', 0, rain)), 6, 30);
  assert.ok(opening.latitude > -15 && opening.latitude < -10, `the band on the opening day: ${opening.latitude}`);
  const olr = zonal(run.values('year', 0, field('outgoingLongwave'))), ts = zonal(run.values('year', 0, field('surfaceTemperature')));
  assert.equal(Math.round(olr[15]), 171); assert.equal(Math.round(ts[15]), 298);
  const frozenNorth = olr.slice(28); assert.ok(Math.min(...frozenNorth) < 171 && Math.max(...frozenNorth) > 171, 'about as little as the frozen north');
  assert.equal(Math.round(Math.max(...olr.slice(10, 26))), 267);
  const september = month('September');
  for (let b = 16; b >= 12; b--) assert.ok(september.precipitation[b] < 0.7, `dry at ${lat(b)}: ${september.precipitation[b]}`);
  const vapor = september.precipitableWater;
  assert.deepEqual([Math.min(...vapor.slice(12, 17)), Math.max(...vapor.slice(12, 17))].map(Math.round), [10, 25]);
  assert.deepEqual([Math.min(...vapor.slice(19, 23)), Math.max(...vapor.slice(19, 23))].map(Math.round), [45, 53]);
  const zm = header.zonalMean, low = Array.from({ length: 36 }, (_, b) => zm.RH[b * K + K - 1]);
  assert.deepEqual([Math.min(...low), Math.max(...low)].map((v) => Math.round(100 * v)), [74, 92]);
  const k500 = sigmaNear(0.5), rh = (b, m = zm) => Math.round(100 * m.RH[b * K + k500]);
  assert.deepEqual([rh(13), rh(22)], [28, 37]);
  const middle = [6, 7, 8, 9, 26, 27, 28, 29].map((b) => rh(b));
  assert.deepEqual([Math.min(...middle), Math.max(...middle)], [52, 56]);
  const k950 = sigmaNear(0.95);
  close(zm.q[13 * K + k950], 8.9, 0.05, 'vapor near the sea at 22.5° S'); close(zm.q[17 * K + k950], 12.5, 0.05, 'at 2.5° S');
  const sep = month('September', 'zonalMean');
  close(sep.v[18 * K + k950], 8.6, 0.05, 'September trades across the equator');
  assert.equal(Math.min(...Array.from({ length: 12 }, (_, n) => rh(n + 6, sep))), 9);
});

test('the prose numbers: the seasons at the sea and the ice', () => {
  const times = run.times('year'), ts = field('surfaceTemperature'), last = times.length - 1, at45 = [], north = [], south = [];
  for (let n = 0; n < last; n++) { const v = run.values('year', n, ts); at45.push(band(v, 42.5, 47.5)); north.push(band(v, 80, 90.1)); south.push(band(v, -90.1, -80)); }
  const half = times.findIndex((t) => t >= 182.5 * DAY), lowest = (values, from, to) => Math.min(...values.slice(from, to)), highest = (values, from, to) => Math.max(...values.slice(from, to));
  const warm = Math.max(...at45), firstCold = lowest(at45, 0, half), lastCold = lowest(at45, half, last);
  close(warm, 294.4, 0.05, 'warmest at 45° N'); assert.equal(date(times[at45.indexOf(warm)] / DAY), 'September 16');
  close(firstCold, 274.2, 0.05, 'coldest at 45° N in the first March'); assert.equal(calendarDate(times[at45.indexOf(firstCold)] / DAY).month, 'March');
  close(lastCold, 268.8, 0.05, 'and in the last'); assert.equal(date(times[at45.indexOf(lastCold)] / DAY), 'March 15');
  assert.deepEqual([warm - firstCold, warm - lastCold].map(Math.round), [20, 26]);
  const swings = [highest(north, 0, last) - lowest(north, 0, half), highest(north, 0, last) - lowest(north, half, last), highest(south, half, last) - lowest(south, 0, last), highest(south, 0, half) - lowest(south, 0, last)];
  assert.deepEqual([Math.min(...swings), Math.max(...swings)].map(Math.round), [33, 45]);
  const dryTs = part6.header.fields.findIndex((f) => f.name === 'surfaceTemperature'), dryTimes = part6.times('year').slice(0, -1);
  const cap = (lo, hi) => dryTimes.map((_, n) => band(part6.values('year', n, dryTs), lo, hi)), range = (v) => Math.max(...v) - Math.min(...v);
  assert.deepEqual([range(cap(80, 90.1)), range(cap(-90.1, -80))].map(Math.round), [89, 87]);
  close(hemisphere(month('September').ice, false), 0.328, 0.0005, 'southern ice in September');
  assert.equal(hemisphere(month('September').ice, true), 0);
  const ice = header.daily.ice, spinupSouth = 2 * Math.max(...ice.slice(150, 200)), recordSouth = 2 * Math.max(...ice.slice(YEAR + 150, YEAR + 200));
  close(spinupSouth, 37.8, 0.2, 'the spin-up’s southern peak, with the north ice-free'); close(recordSouth, 33.0, 0.15, 'the record’s');
});

test('the prose numbers: the jets and the tropopause', () => {
  const zm = header.zonalMean, sigma = zm.sigma;
  const jet = (u, from, to) => { let best = { u: -Infinity }; for (let b = from; b < to; b++) for (let k = 0; k < K; k++) if (sigma[k] > 0.1 && u[b * K + k] > best.u) best = { u: u[b * K + k], lat: lat(b), p: 1000 * sigma[k] }; return best; };
  const south = jet(zm.u, 0, 18), north = jet(zm.u, 18, 36), september = jet(month('September', 'zonalMean').u, 0, 36), march = jet(month('March', 'zonalMean').u, 0, 36);
  assert.deepEqual([Math.round(south.u), south.lat, Math.round(north.u), north.lat], [42, -32.5, 32, 37.5]);
  assert.deepEqual([Math.round(september.u), september.lat, Math.round(march.u), march.lat], [61, -37.5, 57, 32.5]);
  for (const j of [south, north, september, march]) assert.ok(j.p > 160 && j.p < 200, `core near ${j.p} hPa`);
  const sep = month('September', 'zonalMean'), k300 = sigmaNear(0.314), k37 = sigmaNear(0.037), k85 = sigmaNear(0.085);
  assert.deepEqual([sep.u[7 * K + k300], sep.u[7 * K + k37]].map(Math.round), [21, 46]);
  assert.deepEqual([sep.T[1 * K + k85], sep.T[7 * K + k85]].map(Math.round), [196, 213]);
  const equator = (k) => (zm.T[17 * K + k] + zm.T[18 * K + k]) / 2;
  let coldest = 0; for (let k = 0; k < K; k++) if (equator(k) < equator(coldest)) coldest = k;
  close(equator(coldest), 193.3, 0.05, 'coldest over the equator'); close(1000 * sigma[coldest], 100, 1, 'near 100 hPa');
  close(zm.T[30 * K + coldest], 211.6, 0.05, '60° N at that height');
  const k164 = sigmaNear(0.164), k193 = sigmaNear(0.193);
  assert.ok(equator(k164) < zm.T[30 * K + k164] && equator(k193) > zm.T[30 * K + k193], 'the tropics turn colder between 193 and 164 hPa');
});
