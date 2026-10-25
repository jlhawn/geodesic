import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore } from '../js/dynamics/sigmaCore.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { createRadiation, sunDirection, sunDistance, seasonTimes, orbitOf, DEFAULT_ORBIT, SOLAR_CONSTANT, AXIAL_TILT, ECCENTRICITY, PERIHELION, DAY, YEAR } from '../js/physics/radiation.module.js';
import { seasonPhrase } from '../js/levels.module.js';

const D = Math.PI / 180;
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);

function areaLawOrbit(e = ECCENTRICITY, varpi = PERIHELION, step = 3600) {
  const k = 2 * Math.PI / YEAR / (1 - e * e) ** 1.5, rate = (lambda) => k * (1 + e * Math.cos(lambda - varpi)) ** 2;
  const advance = (lambda, h) => {
    const a = rate(lambda), b = rate(lambda + 0.5 * h * a), c = rate(lambda + 0.5 * h * b), d = rate(lambda + h * c);
    return lambda + h / 6 * (a + 2 * b + 2 * c + d);
  };
  const steps = Math.round(YEAR / step), samples = new Float64Array(steps + 1);
  for (let s = 0; s < steps; s++) samples[s + 1] = advance(samples[s], step);
  const longitude = (t) => { const s = Math.min(steps, Math.floor(t / step)); return advance(samples[s], t - s * step); };
  const distance = (lambda) => ((1 + e * Math.cos(lambda - varpi)) / (1 - e * e)) ** 2;
  const timeOf = (lambda) => {
    let lo = 0, hi = YEAR;
    for (let n = 0; n < 60; n++) { const mid = 0.5 * (lo + hi); if (longitude(mid) < lambda) lo = mid; else hi = mid; }
    return 0.5 * (lo + hi);
  };
  return { longitude, distance, timeOf, end: samples[steps] };
}

function textbookDailyMean(lat, declination, distance) {
  const c = -Math.tan(lat) * Math.tan(declination), h0 = c >= 1 ? 0 : c <= -1 ? Math.PI : Math.acos(c);
  return SOLAR_CONSTANT / Math.PI * distance * (h0 * Math.sin(lat) * Math.sin(declination) + Math.cos(lat) * Math.cos(declination) * Math.sin(h0));
}

const law = areaLawOrbit();
const declinationAt = (t) => Math.asin(Math.sin(AXIAL_TILT) * Math.sin(law.longitude(t)));
const textbookAt = (lat, t) => textbookDailyMean(lat, declinationAt(t), law.distance(law.longitude(t)));

function dailyMean(lat, t, orbit, samples = 288) {
  const x = Math.cos(lat), z = Math.sin(lat), sun = new Float64Array(3);
  let total = 0;
  for (let k = 0; k < samples; k++) {
    const at = t - DAY / 2 + (k + 0.5) * DAY / samples;
    sunDirection(at, sun, orbit);
    total += SOLAR_CONSTANT * sunDistance(at, orbit) * Math.max(0, x * sun[0] + z * sun[2]);
  }
  return total / samples;
}

test('the integral of Kepler\'s second law takes the year to turn the sun once round', () => {
  close(law.end, 2 * Math.PI, 1e-9, 'longitude after YEAR');
  close(law.longitude(0), 0, 0, 'longitude at the spring equinox');
});

test('the kepler sun\'s declination and distance are those of the area law, with the equinox at t = 0 and the hour angle of the circular orbit', () => {
  const out = new Float64Array(3), circle = new Float64Array(3);
  for (let n = 0; n < 150; n++) {
    const t = n * 2.43 * DAY + n * 3671;
    sunDirection(t, out, 'kepler'); sunDirection(t, circle, 'circular');
    const lambda = law.longitude(t % YEAR);
    close(out[2], Math.sin(AXIAL_TILT) * Math.sin(lambda), 1e-9, `sin δ, t = ${t}`);
    close(Math.hypot(out[0], out[1], out[2]), 1, 1e-14, `unit vector, t = ${t}`);
    close(sunDistance(t, 'kepler') / law.distance(lambda), 1, 1e-9, `(a/r)², t = ${t}`);
    close(Math.atan2(out[1], out[0]), Math.atan2(circle[1], circle[0]), 1e-12, `hour angle, t = ${t}`);
  }
  close(sunDirection(0, undefined, 'kepler')[2], 0, 1e-12, 'spring equinox');
  assert.equal(sunDistance(123456, 'circular'), 1);
});

test('the kepler orbit is the default: of the sun functions, of the radiation and of the page\'s options', () => {
  assert.equal(DEFAULT_ORBIT, 'kepler');
  assert.equal(orbitOf(undefined), 'kepler');
  assert.equal(orbitOf({}), 'kepler');
  assert.equal(orbitOf({ radiation: {} }), 'kepler');
  assert.equal(orbitOf({ radiation: { orbit: 'circular' } }), 'circular');
  assert.equal(orbitOf({ radiation: { orbit: 'kepler' }, moist: {} }), 'kepler');
  const defaulted = new Float64Array(3), kepler = new Float64Array(3);
  for (const t of [0, 30 * DAY + 77, 100.3 * DAY, 250 * DAY]) {
    sunDirection(t, defaulted); sunDirection(t, kepler, 'kepler');
    assert.deepEqual(defaulted, kepler, `sun at ${t}`);
    assert.equal(sunDistance(t), sunDistance(t, 'kepler'));
  }
  const mesh = buildMesh(new Grid(16)), core = createSigmaCore(mesh), plain = createRadiation(mesh, core), explicit = createRadiation(mesh, core, { orbit: 'kepler' }), undecided = createRadiation(mesh, core, { orbit: undefined });
  for (const t of [0, 40.5 * DAY, 200.2 * DAY]) {
    plain.setTime(t); explicit.setTime(t); undecided.setTime(t);
    for (let i = 0; i < mesh.nCells; i += 7) { assert.equal(plain.insolation(i), explicit.insolation(i)); assert.equal(undecided.insolation(i), explicit.insolation(i)); }
  }
  const model = createModel(new Grid(16), { ocean: false, surface: { exchange: 'roughness' } }), C = model.mesh.nCells;
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  model.time = 150 * DAY;
  model.step(900);
  model.radiation.setTime(150 * DAY);
  const f = sunDistance(150 * DAY, 'kepler');
  let checked = 0;
  for (let i = 0; i < C; i++) {
    if (!(model.radiation.cosZenith(i) > 0)) continue;
    close(model.radiation.summed.insolation[i], SOLAR_CONSTANT * f * model.radiation.cosZenith(i), 1e-9, `cell ${i}`);
    checked++;
  }
  assert.ok(checked > C / 3, `${checked} cells lit`);
});

test('the June solstice falls on day 92.697, the September equinox on 186.282 and the December solstice on 276.064 of the 365-day year, perihelion on day 288.715 and aphelion on 106.215', () => {
  const days = { june: law.timeOf(90 * D) / DAY, september: law.timeOf(180 * D) / DAY, december: law.timeOf(270 * D) / DAY, perihelion: law.timeOf(PERIHELION) / DAY, aphelion: law.timeOf(PERIHELION - Math.PI) / DAY };
  console.log(`kepler orbit, e = ${ECCENTRICITY}, ϖ = ${(PERIHELION / D).toFixed(1)}°, tilt ${(AXIAL_TILT / D).toFixed(2)}°: June solstice day ${days.june.toFixed(3)}, September equinox ${days.september.toFixed(3)}, December solstice ${days.december.toFixed(3)}, perihelion ${days.perihelion.toFixed(3)}, aphelion ${days.aphelion.toFixed(3)}; (a/r)² ${sunDistance(days.june * DAY, 'kepler').toFixed(5)}, ${sunDistance(days.december * DAY, 'kepler').toFixed(5)}`);
  close(days.june, 92.697, 0.002, 'June solstice');
  close(days.september, 186.282, 0.002, 'September equinox');
  close(days.december, 276.064, 0.002, 'December solstice');
  close(days.perihelion, 288.715, 0.002, 'perihelion');
  close(days.aphelion, 106.215, 0.002, 'aphelion');
  seasonTimes('kepler').forEach((at, k) => close(at / DAY, [0, days.june, days.september, days.december][k], 1e-4, `seasonTimes ${k}`));
  assert.deepEqual(seasonTimes('circular'), [0, YEAR / 4, YEAR / 2, 3 * YEAR / 4]);
  const sun = new Float64Array(3);
  sunDirection(days.june * DAY, sun, 'kepler');
  close(Math.asin(sun[2]), AXIAL_TILT, 1e-9, 'declination at the June solstice');
  sunDirection(days.december * DAY, sun, 'kepler');
  close(Math.asin(sun[2]), -AXIAL_TILT, 1e-9, 'declination at the December solstice');
  sunDirection(days.september * DAY, sun, 'kepler');
  close(sun[2], 0, 1e-9, 'declination at the September equinox');
  const near = (day) => { sunDirection(day * DAY, sun, 'kepler'); return sun[2]; };
  assert.ok(near(days.june) > near(days.june - 1) && near(days.june) > near(days.june + 1), 'the declination peaks at the solstice');
  close(sunDistance(days.perihelion * DAY, 'kepler'), 1 / (1 - ECCENTRICITY) ** 2, 1e-9, '(a/r)² at perihelion');
  close(sunDistance(days.aphelion * DAY, 'kepler'), 1 / (1 + ECCENTRICITY) ** 2, 1e-9, '(a/r)² at aphelion');
});

test('daily means at the model\'s cells match S₀/π (a/r)² (h₀ sin φ sin δ + cos φ cos δ sin h₀) at the area law\'s declination and distance of the local noon', () => {
  const mesh = buildMesh(new Grid(16)), radiation = createRadiation(mesh, createSigmaCore(mesh), { orbit: 'kepler' });
  const samples = 288;
  let worst = 0, at = '';
  for (const day of [0, 40, 92, 140, 186, 230, 276, 330]) {
    for (let i = 0; i < mesh.nCells; i++) {
      const lat = mesh.latCell[i];
      if (Math.abs(lat) > 62 * D) continue;
      const turn = -Math.atan2(mesh.xCell[3 * i + 1], mesh.xCell[3 * i]) / (2 * Math.PI), noon = (day + turn - Math.floor(turn)) * DAY;
      let mean = 0;
      for (let k = 0; k < samples; k++) {
        radiation.setTime(noon - DAY / 2 + (k + 0.5) * DAY / samples);
        mean += radiation.insolation(i) / samples;
      }
      const error = Math.abs(mean - textbookAt(lat, noon));
      if (error > worst) { worst = error; at = `day ${day}, ${(lat / D).toFixed(1)}°`; }
    }
  }
  console.log(`daily means of the kepler orbit on the N = 16 mesh against the textbook form: worst difference ${worst.toFixed(3)} W/m² (${at})`);
  assert.ok(worst < 0.1, `worst difference ${worst} W/m² at ${at}`);
});

test('50°N at the June solstice gets 3.2 % less than on the circular orbit, 482 against 498 W/m²', () => {
  const lat = 50 * D, june = law.timeOf(90 * D), circular = dailyMean(lat, YEAR / 4, 'circular'), kepler = dailyMean(lat, june, 'kepler');
  console.log(`50°N at the June solstice: ${circular.toFixed(2)} W/m² circular, ${kepler.toFixed(2)} kepler (${(100 * (kepler / circular - 1)).toFixed(2)} %)`);
  close(circular, textbookDailyMean(lat, AXIAL_TILT, 1), 0.3, 'circular against the textbook form');
  close(kepler, textbookDailyMean(lat, AXIAL_TILT, sunDistance(june, 'kepler')), 0.3, 'kepler against the textbook form');
  close(circular, 498, 0.5, 'circular');
  close(kepler, 482, 0.5, 'kepler');
  close(kepler / circular - 1, -0.0318, 0.0005, 'ratio');
});

test('the annual global mean rises by (1 − e²)^(−1/2) and no more', () => {
  const steps = 365 * 48;
  let sum = 0;
  for (let n = 0; n < steps; n++) sum += sunDistance((n + 0.5) / steps * YEAR, 'kepler');
  close(sum / steps, 1 / Math.sqrt(1 - ECCENTRICITY * ECCENTRICITY), 1e-9, 'time mean of (a/r)²');
  const mesh = buildMesh(new Grid(16)), kepler = createRadiation(mesh, createSigmaCore(mesh), { orbit: 'kepler' }), circular = createRadiation(mesh, createSigmaCore(mesh), { orbit: 'circular' });
  let total = { kepler: 0, circular: 0 }, area = 0;
  for (let i = 0; i < mesh.nCells; i++) area += mesh.areaCell[i];
  const perDay = 6, count = 365 * perDay;
  for (let n = 0; n < count; n++) {
    const t = (n + 0.5) / count * YEAR;
    kepler.setTime(t); circular.setTime(t);
    for (let i = 0; i < mesh.nCells; i++) { total.kepler += mesh.areaCell[i] * kepler.insolation(i); total.circular += mesh.areaCell[i] * circular.insolation(i); }
  }
  console.log(`annual global mean at N = 16: ${(total.circular / count / area).toFixed(3)} W/m² circular, ${(total.kepler / count / area).toFixed(3)} kepler; S₀/4 ${(SOLAR_CONSTANT / 4).toFixed(3)}, ratio ${(total.kepler / total.circular).toFixed(6)} against ${(1 / Math.sqrt(1 - ECCENTRICITY * ECCENTRICITY)).toFixed(6)}`);
  close(total.kepler / total.circular, 1 / Math.sqrt(1 - ECCENTRICITY * ECCENTRICITY), 5e-5, 'ratio of the annual means on the mesh');
});

test('the half-years between the equinoxes receive equal totals, Kepler\'s second law, where the distance factor alone would take 2 % from the northern summer', () => {
  const september = law.timeOf(180 * D), n = 40000;
  const total = (from, to, f) => { let sum = 0; for (let k = 0; k < n; k++) sum += f(from + (k + 0.5) * (to - from) / n); return sum * (to - from) / n; };
  const spring = total(0, september, (t) => sunDistance(t, 'kepler')), autumn = total(september, YEAR, (t) => sunDistance(t, 'kepler')), factorOnly = total(0, YEAR / 2, (t) => sunDistance(t, 'kepler'));
  console.log(`global insolation in S₀/4 from the March to the September equinox ${(spring / (YEAR / 2)).toFixed(5)} against ${(autumn / (YEAR / 2)).toFixed(5)} from the September to the March (each in units of half a year); with the distance factor alone over the circular half-year ${(factorOnly / (YEAR / 2)).toFixed(5)}`);
  close(spring / autumn, 1, 2e-6, 'ratio of the half-year totals');
  close(spring / (YEAR / 2), 1 / Math.sqrt(1 - ECCENTRICITY * ECCENTRICITY), 1e-8, 'each half-year against the mean of (a/r)²');
  assert.ok(factorOnly / (YEAR / 2) < 0.985 && factorOnly / (YEAR / 2) > 0.975, `the distance factor alone changes the northern summer by ${factorOnly / (YEAR / 2)}`);
});

test('at every latitude the annual mean on the kepler orbit is the circular timing\'s times (1 − e²)^(−1/2): the timing and the distance together', () => {
  const flat = areaLawOrbit(0), sun = new Float64Array(3), perDay = 360, steps = 365 * perDay, step = YEAR / steps;
  for (const degrees of [0, 30, 50, 80]) {
    const x = Math.cos(degrees * D), z = Math.sin(degrees * D);
    let kepler = 0, uniform = 0;
    for (let n = 0; n < steps; n++) {
      const t = (n + 0.5) * step;
      sunDirection(t, sun, 'kepler');
      kepler += SOLAR_CONSTANT * sunDistance(t, 'kepler') * Math.max(0, x * sun[0] + z * sun[2]);
      const lambda = flat.longitude(t), declination = Math.asin(Math.sin(AXIAL_TILT) * Math.sin(lambda)), spin = -(t % DAY) * 2 * Math.PI / DAY;
      uniform += SOLAR_CONSTANT * Math.max(0, x * Math.cos(declination) * Math.cos(spin) + z * Math.sin(declination));
    }
    console.log(`${degrees}°N annual mean: ${(kepler / steps).toFixed(3)} W/m² on the kepler orbit, ${(uniform / steps).toFixed(3)} on the circular timing with the same declination formula; ratio ${(kepler / uniform).toFixed(6)} against ${(1 / Math.sqrt(1 - ECCENTRICITY * ECCENTRICITY)).toFixed(6)}`);
    close(kepler / uniform, 1 / Math.sqrt(1 - ECCENTRICITY * ECCENTRICITY), 2e-5, `${degrees}°N`);
  }
});

test('every cell gets the model\'s own sunlight scaled by the distance factor, and the circular orbit still gets exactly S₀ times the cosine', () => {
  const mesh = buildMesh(new Grid(16)), core = createSigmaCore(mesh), kepler = createRadiation(mesh, core), circular = createRadiation(mesh, core, { orbit: 'circular' }), sun = new Float64Array(3);
  for (const t of [0, 0.3 * DAY, 91.25 * DAY + 5000, 200.7 * DAY, 300 * DAY + 77]) {
    kepler.setTime(t); circular.setTime(t);
    sunDirection(t, sun, 'kepler');
    const f = sunDistance(t, 'kepler');
    let lit = 0;
    for (let i = 0; i < mesh.nCells; i++) {
      const cosine = Math.max(0, mesh.xCell[3 * i] * sun[0] + mesh.xCell[3 * i + 1] * sun[1] + mesh.xCell[3 * i + 2] * sun[2]);
      close(kepler.insolation(i), SOLAR_CONSTANT * f * cosine, 1e-9, `kepler cell ${i} at t = ${t}`);
      close(kepler.cosZenith(i), cosine, 1e-14, `kepler cosine of cell ${i}`);
      if (cosine > 0) lit++;
      const cosineCircular = circular.cosZenith(i);
      assert.equal(circular.insolation(i), SOLAR_CONSTANT * cosineCircular, `circular cell ${i} at t = ${t}`);
    }
    assert.ok(lit > mesh.nCells / 3, `${lit} cells lit`);
  }
});

test('the page names the season by the orbit\'s own equinoxes and solstices', () => {
  const [, june, september, december] = seasonTimes('kepler');
  assert.equal(seasonPhrase(0, 'kepler'), 'the Northern spring equinox');
  assert.equal(seasonPhrase(june, 'kepler'), 'the Northern summer solstice');
  assert.equal(seasonPhrase(september + 2 * DAY, 'kepler'), '2 days past the Northern autumn equinox');
  assert.equal(seasonPhrase(december - 3 * DAY, 'kepler'), '3 days before the Northern winter solstice');
  assert.equal(seasonPhrase(YEAR - DAY, 'kepler'), '1 day before the Northern spring equinox');
  assert.equal(seasonPhrase(182.5 * DAY, 'kepler'), '4 days before the Northern autumn equinox');
  assert.equal(seasonPhrase(91.25 * DAY, 'circular'), 'the Northern summer solstice');
  assert.equal(seasonPhrase(185.5 * DAY, 'circular'), '3 days past the Northern autumn equinox');
  assert.equal(seasonPhrase(60 * DAY), seasonPhrase(60 * DAY, DEFAULT_ORBIT));
});

test('the orbit is circular or kepler', () => {
  const mesh = buildMesh(new Grid(16)), core = createSigmaCore(mesh);
  for (const orbit of ['ellipse', 'Kepler', null, 1, '']) assert.throws(() => createRadiation(mesh, core, { orbit }), /orbit must be/, String(orbit));
});

function heldModel(radiationEvery, orbit) {
  const model = createModel(new Grid(16), { ocean: false, radiation: { stratus: true, clearSkyPass: true, radiationEvery, orbit }, surface: { exchange: 'roughness' } });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  return model;
}

test('with the kepler orbit what the radiation absorbs and reflects makes the step\'s sunlight, every step, held or not', () => {
  for (const every of [1, 4]) {
    const model = heldModel(every, 'kepler'), radiation = model.radiation, C = model.mesh.nCells, dt = 900, start = 150 * DAY + 2 * dt;
    model.time = start;
    const absorbed = new Float64Array(C), reflected = new Float64Array(C), insolation = new Float64Array(C);
    let worst = 0, lit = 0, shone = 0;
    for (let n = 0; n < 9; n++) {
      const time = model.time;
      model.step(dt);
      radiation.setTime(time);
      const f = sunDistance(time, 'kepler');
      for (let i = 0; i < C; i++) {
        const beam = radiation.summed.insolation[i] - insolation[i], top = radiation.summed.absorbedSolar[i] - absorbed[i], lost = radiation.summed.reflectedSolar[i] - reflected[i];
        close(beam, SOLAR_CONSTANT * f * radiation.cosZenith(i), 1e-9, `step ${n} cell ${i}: insolation`);
        worst = Math.max(worst, Math.abs(top + lost - beam));
        if (beam > 0) lit++;
        if (top > 0) shone++;
      }
      absorbed.set(radiation.summed.absorbedSolar); reflected.set(radiation.summed.reflectedSolar); insolation.set(radiation.summed.insolation);
      radiation.setTime(model.time);
    }
    console.log(`radiationEvery ${every} on the kepler orbit: absorbed plus reflected differs from the insolation by at most ${worst.toExponential(1)} W/m² over ${lit} lit cell-steps`);
    assert.ok(lit > 3 * C && shone > 3 * C, `${lit} lit and ${shone} absorbing cell-steps`);
    assert.ok(worst < 1e-8, `absorbed plus reflected misses the insolation by ${worst} W/m²`);
    if (every > 1) {
      let held = 0;
      for (let i = 0; i < C; i++) if (radiation.held.mu[i] > 0) held++;
      assert.ok(held > C / 4, `${held} columns hold a lit call`);
    }
  }
});
