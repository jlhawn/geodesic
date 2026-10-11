import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { readFrames } from '../explain/frames.module.js';
import { calendarDate, season, heatTransport, storedHeat, RADIUS } from '../explain/sunlitRun.module.js';

const open = (name) => { const file = new URL(`../explain/data/${name}`, import.meta.url), bytes = readFileSync(file); return readFrames(file.href, bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength)); };
const run = open('radiation_N16.bin'), part5 = open('heldsuarez_N16.bin');
const { header } = run, DAY = 86400, YEAR = 365;
const cells = [...new Grid(header.N)].map((c) => ({ lat: Math.asin(c.centerVertex.z) * 180 / Math.PI, area: c.area }));
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const mean = (values) => values.reduce((sum, v) => sum + v, 0) / values.length;
const field = (name) => header.fields.findIndex((f) => f.name === name);
const stored = storedHeat(run, cells), corrected = heatTransport(header.zonalBudget, stored.total), uncorrected = heatTransport(header.zonalBudget);
const zonal = (name, h = header) => h.zonalMean.months.find((m) => m.name === name);
const band = (values, lo, hi) => {
  let sum = 0, area = 0;
  cells.forEach((c, i) => { if (c.lat >= lo && c.lat < hi) { sum += values[i] * c.area; area += c.area; } });
  return sum / area;
};

test('the run’s calendar starts at noon on March 20 and its January is days 287 to 317', () => {
  assert.deepEqual(calendarDate(0), { month: 'March', date: 20, hour: 12 });
  assert.deepEqual(calendarDate(287), { month: 'January', date: 1, hour: 12 });
  assert.deepEqual(calendarDate(287.5), { month: 'January', date: 2, hour: 0 });
  assert.deepEqual(calendarDate(145), { month: 'August', date: 12, hour: 12 });
  assert.deepEqual(calendarDate(364.75), { month: 'March', date: 20, hour: 6 });
  assert.deepEqual([91.25, 182.5, 273.75].map((d) => { const c = calendarDate(d); return `${c.month} ${c.date}`; }), ['June 19', 'September 19', 'December 19']);
  const weather = run.times('weather');
  assert.equal(weather.length, 121);
  assert.equal(calendarDate(weather[0] / DAY).month, 'January');
  assert.equal(calendarDate(weather.at(-1) / DAY).date, 31);
});

test('seasons turn at the equinoxes and solstices, opposite in the two hemispheres', () => {
  assert.deepEqual(season(0), { north: 'spring', south: 'fall' });
  assert.deepEqual(season(91), { north: 'summer', south: 'winter' });
  assert.deepEqual(season(200), { north: 'fall', south: 'spring' });
  assert.deepEqual(season(300), { north: 'winter', south: 'summer' });
  assert.deepEqual(season(365), season(0));
});

test('the heat carried is zero at both poles and peaks where the net, less the heat stored, changes sign', () => {
  const t = corrected;
  assert.equal(t.transport[0], 0);
  close(t.transport.at(-1), 0, 1e-12, 'north pole');
  close(t.imbalance, t.absorbedMean - t.outgoingMean - t.storedMean, 1e-12, 'imbalance');
  const area = 4 * Math.PI * RADIUS ** 2, raw = t.net.reduce((sum, v, b) => sum + v * 2 * Math.PI * RADIUS ** 2 * (Math.sin(t.edges[b + 1] * Math.PI / 180) - Math.sin(t.edges[b] * Math.PI / 180)), 0);
  close(raw / area, t.imbalance, 1e-12, 'the bands tile the sphere');
  assert.equal(t.crossings.length, 2);
  for (const lat of t.crossings) close(t.at(lat).net, 0, 1e-9, `net at ${lat}`);
  const north = t.peak(1), south = t.peak(-1);
  assert.ok(Math.abs(north.latitude - t.crossings[1]) < 5 && Math.abs(south.latitude - t.crossings[0]) < 5, 'the peaks sit at the crossings');
  for (const lat of [-60, 0, 40]) close(t.at(lat).absorbed - t.at(lat).outgoing - t.at(lat).stored, t.at(lat).net, 1e-12, `net at ${lat}`);
});

test('the heat the ground and air stored over the year, and what is left over', () => {
  const area = (b) => Math.sin(corrected.edges[b + 1] * Math.PI / 180) - Math.sin(corrected.edges[b] * Math.PI / 180);
  const global = (values) => values.reduce((sum, v, b) => sum + v * area(b), 0) / 2;
  close(global(stored.ground), -0.565, 0.0005, 'ground');
  close(global(stored.air), -0.193, 0.0005, 'air');
  close(corrected.storedMean, -0.758, 0.0005, 'ground and air');
  close(uncorrected.imbalance, -0.087, 0.0005, 'top of the air');
  close(corrected.imbalance, 0.67, 0.005, 'left over and spread evenly');
  close(stored.ground[0] * YEAR * DAY / header.ground.heatCapacity, -5.70, 0.005, 'south polar ground');
  const share = [0, 35].map((b) => stored.total[b] / uncorrected.net[b]);
  close(share[0], 0.20, 0.005, 'a fifth in the south');
  close(share[1], 0.31, 0.005, 'nearly a third in the north');
});

test('the prose numbers: balance, the gap by latitude and the heat carried', () => {
  const t = corrected;
  close(t.absorbedMean, 295.9, 0.05, 'sunlight absorbed');
  close(t.outgoingMean, 296.0, 0.05, 'infrared out');
  const equator = t.at(0), poles = [t.absorbed[0], t.absorbed.at(-1), t.outgoing[0], t.outgoing.at(-1)];
  close(equator.absorbed, 363.6, 0.05, 'equator sunlight');
  close(equator.outgoing, 351.0, 0.05, 'equator infrared');
  close(equator.absorbed - equator.outgoing, 12.6, 0.05, 'equator surplus');
  close(poles[0], 145.4, 0.05, 'south pole sunlight');
  close(poles[1], 145.4, 0.05, 'north pole sunlight');
  close((poles[2] + poles[3]) / 2, 172, 0.5, 'polar infrared');
  close(poles[2] - poles[0], 26.9, 0.05, 'south polar loss');
  close(poles[3] - poles[1], 26.3, 0.05, 'north polar loss');
  close(t.crossings[0], -40.8, 0.05, 'southern crossing');
  close(t.crossings[1], 39.7, 0.05, 'northern crossing');
  assert.deepEqual([t.peak(1).latitude, t.peak(-1).latitude], [40, -40]);
  close(t.peak(1).transport, 1.10, 0.005, 'northward peak');
  close(t.peak(-1).transport, -1.08, 0.005, 'southward peak');
  close(t.at(0).transport, 0.02, 0.005, 'across the equator');
  close(uncorrected.peak(-1).transport, -1.27, 0.005, 'southward peak without the correction');
  close(equator.outgoing - (poles[2] + poles[3]) / 2, 179, 0.5, 'infrared drop');
  close(equator.absorbed - (poles[0] + poles[1]) / 2, 218, 0.5, 'sunlight drop');
});

test('the prose numbers: the year at the ground', () => {
  const from = -header.daily.from / DAY, year = (key) => mean(header.daily[key].slice(from, from + YEAR));
  close(year('surfaceTemperature'), 279.3, 0.05, 'global ground');
  close(year('absorbedSolar'), 295.9, 0.1, 'daily sunlight absorbed');
  close(year('outgoingLongwave'), 296.0, 0.1, 'daily infrared out');
  close(1 - year('absorbedSolar') / 340.5, 0.131, 0.0005, 'planetary albedo');
  close(year('surfaceTemperature') - (year('absorbedSolar') / 5.670374419e-8) ** 0.25, 10.5, 0.05, 'above the effective temperature');
  const Ts = header.zonalBudget.surfaceTemperature, equator = (Ts[17] + Ts[18]) / 2, pole = (Ts[0] + Ts[35]) / 2;
  close(equator, 295.8, 0.05, 'equator');
  close(pole, 236.7, 0.05, 'poles');
  close(equator - pole, 59, 0.5, 'equator to pole');
  const ts = field('surfaceTemperature'), times = run.times('year'), last = times.length - 1, at45 = [], north = [], south = [], equatorial = [];
  for (let n = 0; n <= last; n++) {
    const values = run.values('year', n, ts);
    north.push(band(values, 80, 90.1)); south.push(band(values, -90.1, -80));
    if (n < last) { at45.push(band(values, 42.5, 47.5)); equatorial.push(band(values, -5, 5)); }
  }
  close(Math.min(...equatorial), 293.9, 0.05, 'equator coolest');
  close(Math.max(...equatorial), 297.9, 0.05, 'equator warmest');
  const warm = Math.max(...at45), cold = Math.min(...at45);
  close(warm, 298.8, 0.05, 'warmest at 45° N');
  close(cold, 244.6, 0.05, 'coldest at 45° N');
  assert.deepEqual([calendarDate(times[at45.indexOf(warm)] / DAY), calendarDate(times[at45.indexOf(cold)] / DAY)].map((c) => `${c.month} ${c.date}`), ['August 12', 'February 8']);
  const date = (series, value) => { const c = calendarDate(times[series.indexOf(value)] / DAY); return `${c.month} ${c.date}`; };
  const summer = Math.max(...north), floor = Math.min(...north), southFloor = Math.min(...south);
  close(summer, 285.7, 0.05, 'northern polar summer');
  assert.equal(date(north, summer), 'August 2');
  close(floor, 196.5, 0.05, 'end of the northern polar night');
  assert.equal(north.indexOf(floor), last);
  close(southFloor, 197.3, 0.05, 'end of the southern polar night');
  assert.equal(date(south, southFloor), 'September 26');
  close(summer - floor, 89, 0.5, 'northern polar swing');
  close(north[0], 201.9, 0.05, 'the northern night a year earlier');
});

test('the prose numbers: the planet was still cooling, most of it in the south, the poles fastest at first', () => {
  const ts = field('surfaceTemperature'), times = run.times('year'), last = times.length - 1;
  const start = run.values('year', 0, ts), end = run.values('year', last, ts), change = end.map((v, i) => v - start[i]);
  close(band(change, -91, 91), -0.85, 0.005, 'globe');
  close(band(change, -91, 0), -1.51, 0.005, 'south');
  close(band(change, 0, 91), -0.18, 0.005, 'north');
  close(header.daily.surfaceTemperature[-header.daily.from / DAY - 1], 279.40, 0.005, 'ground at the start of the record');
  close(header.daily.surfaceTemperature.at(-1), 278.55, 0.005, 'ground at the end of the record');
  const cap = (n, north) => band(run.values('year', n, ts), north ? 80 : -90.1, north ? 90.1 : -80);
  const halfYearLater = (n, north) => (cap(n + 36, north) + cap(n + 37, north)) / 2;
  const gap = (n) => [cap(n, false) - halfYearLater(n, true), cap(n, true) - halfYearLater(n, false)];
  for (const g of gap(0)) assert.ok(g > 4, `each pole 4 K warmer than the other half a year later at first: ${g}`);
  for (const g of gap(30)) assert.ok(g < 1.5, `and about 1 K by day 150: ${g}`);
});

test('the prose numbers: the inversion and the infrared over the poles', () => {
  const weather = run.times('weather'), ts = field('surfaceTemperature'), t850 = field('temperature850');
  let ground = 0, air = 0;
  for (let n = 0; n < weather.length - 1; n++) { ground += band(run.values('weather', n, ts), 80, 90.1); air += band(run.values('weather', n, t850), 80, 90.1); }
  ground /= weather.length - 1; air /= weather.length - 1;
  close(ground, 210.6, 0.05, 'January ground near the North Pole');
  close(air - ground, 17.1, 0.05, 'the air 1.5 km up is warmer');
  const times = run.times('year'), n = times.findIndex((t) => calendarDate(t / DAY).month === 'August' && calendarDate(t / DAY).date === 12);
  const olr = run.values('year', n, field('outgoingLongwave'));
  let hi = 0, lo = 0;
  olr.forEach((v, i) => { if (v > olr[hi]) hi = i; if (v < olr[lo]) lo = i; });
  close(olr[hi], 388.6, 0.05, 'brightest');
  close(cells[hi].lat, 21.8, 0.05, 'near 22° north');
  assert.ok(olr[lo] < 100 && cells[lo].lat < -80, 'dark over the winter pole');
});

test('the prose numbers: the air near the ground against Held and Suarez and Part 5', () => {
  const sigma = header.zonalMean.sigma, K = header.K, k = sigma.findIndex((s) => Math.abs(s - 0.93) < 0.005);
  assert.equal(k, 24);
  close(287 * 280 / 9.80616 * Math.log(1 / sigma[k]), 600, 10, 'about 600 m up');
  const at = (T, layer, lat) => { const x = (lat + 90) / 5 - 0.5, b = Math.floor(x), t = x - b; return T[b * K + layer] * (1 - t) + T[(b + 1) * K + layer] * t; };
  const contrast = (T, layer) => at(T, layer, 0) - (at(T, layer, 60) + at(T, layer, -60)) / 2;
  close(contrast(header.zonalMean.T, k), 32.2, 0.05, 'equator to 60° about 600 m up');
  close(contrast(part5.header.zonalMean.T, k), 37.0, 0.05, 'Part 5 at the same height');
  const T = header.zonalMean.T, low = K - 1;
  close(at(T, low, 0) - (T[low] + T[35 * K + low]) / 2, 56.3, 0.05, 'the lowest layer from the equator to the poles');
});

test('the prose numbers: the westerlies grow to the model top, weak at 250 hPa, and the winds near the ground', () => {
  const { zonalMean } = header, K = header.K, sigma = zonalMean.sigma, k = sigma.findIndex((s) => s > 0.25), lat = (b) => -90 + (b + 0.5) * 5;
  const at250 = (u, b) => { const w = Math.log(0.25 / sigma[k - 1]) / Math.log(sigma[k] / sigma[k - 1]); return u[b * K + k - 1] * (1 - w) + u[b * K + k] * w; };
  const jet = (u, north) => {
    let b = north ? 18 : 0;
    for (let c = b; c < (north ? 36 : 18); c++) if (at250(u, c) > at250(u, b)) b = c;
    const lo = at250(u, b - 1), mid = at250(u, b), hi = at250(u, b + 1);
    return { u: mid, lat: lat(b + 0.5 * (lo - hi) / (lo - 2 * mid + hi)) };
  };
  close(jet(zonalMean.u, true).u, 19.2, 0.05, 'northern annual mean');
  close(jet(zonalMean.u, false).u, 17.9, 0.05, 'southern annual mean');
  close(jet(zonal('April').u, true).u, 34.6, 0.05, 'northern westerlies in April');
  close(jet(zonal('April').u, true).lat, 43.6, 0.05, 'near 44° north');
  close(jet(zonal('September').u, false).u, 29.9, 0.05, 'southern westerlies in September');
  close(jet(part5.header.zonalMean.u, true).u, 32.2, 0.05, 'Part 5 at 250 hPa');
  for (const [month, b] of [['January', 27], ['July', 9]]) {
    const u = zonal(month).u;
    for (let layer = 0; layer < k; layer++) assert.ok(u[b * K + layer] > u[b * K + layer + 1], `${month} at ${lat(b)}: stronger with height above 250 hPa`);
  }
  const top = Math.max(...zonalMean.months.flatMap((m) => Array.from({ length: 36 }, (_, b) => m.u[b * K])));
  close(top, 150, 0.5, 'top layer');
  close(sigma[0] * 1000, 1.1, 0.05, 'top layer near 1 hPa');
  const ground = (u) => Array.from({ length: 36 }, (_, b) => u[b * K + K - 1]);
  close(Math.max(...ground(zonalMean.u)), 2.3, 0.05, 'westerlies near the ground over the year');
  close(Math.max(...ground(part5.header.zonalMean.u)), 9.3, 0.05, 'Part 5 near the ground');
  for (const [month, b] of [['January', 20], ['July', 15]]) {
    const u = ground(zonal(month).u), trades = u.map((v, c) => (Math.abs(lat(c)) < 35 ? v : 0));
    assert.equal(trades.indexOf(Math.min(...trades)), b, `${month}: strongest trades at ${lat(b)}`);
    close(-u[b], 8.7, 0.15, `${month}: trades about 9 m/s`);
  }
});

test('the prose numbers: January storms against Part 5', () => {
  const lats = cells.map((c) => c.lat);
  const spread = (values, lo, hi) => {
    const rings = new Map();
    lats.forEach((l, i) => { if (l >= lo && l < hi) { const key = Math.round(l / 2.5); if (!rings.has(key)) rings.set(key, []); rings.get(key).push(values[i]); } });
    let sum = 0, count = 0;
    for (const ring of rings.values()) { const m = mean(ring); for (const v of ring) { sum += (v - m) ** 2; count++; } }
    return Math.sqrt(sum / count);
  };
  const storms = (frames, sequence, lo, hi) => { const k = frames.header.fields.findIndex((f) => f.name === 'temperature850'); return mean(frames.times(sequence).map((_, n) => spread(frames.values(sequence, n, k), lo, hi))); };
  close(storms(run, 'weather', 35, 65), 4.7, 0.05, 'northern winter');
  close(storms(run, 'weather', -65, -35), 1.5, 0.05, 'southern summer');
  close(storms(part5, 'settled', 35, 65), 4.0, 0.05, 'Part 5');
});
