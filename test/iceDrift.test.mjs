import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createOcean } from '../js/ocean/layered.module.js';
import { ICE_DRIFT, iceDriftOptions } from '../js/ocean/iceDrift.module.js';
import { createSeaIce, FREEZING_POINT } from '../js/physics/ice.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { DEG, RHO_I, meshOf, rig, packRig, packSolve, packDifference, windOn, rotationOn, streamOn, nearestCell, localFrame, edgesNear, closedForm, unit, dot, cross } from './helpers/iceDrift.mjs';

const near = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const bits = (a, b) => a.length === b.length && a.every((v, i) => Object.is(v, b[i]));
const withinCap = (m, lat, lon, degrees) => {
  const c = [Math.cos(lat * DEG) * Math.cos(lon * DEG), Math.cos(lat * DEG) * Math.sin(lon * DEG), Math.sin(lat * DEG)];
  return Array.from({ length: m.nCells }, (_, i) => i).filter((i) => Math.acos(Math.min(1, dot(unit(m.xCell, i), c))) < degrees * DEG);
};

function freeDrift(N, latDeg, { theta = 0, h = 0.5, A = 0.8, wind = 0.1, current = 0, tilt = 0 } = {}) {
  const m = meshOf(N), r = rig(m, { rheology: 'free', waterTurn: theta, freePasses: 12 });
  const c0 = nearestCell(m, latDeg * DEG, 0), frame = localFrame(m, c0);
  r.thickness.fill(h); r.concentration.fill(A);
  const tau3 = frame.east.map((v) => wind * v), current3 = frame.east.map((v) => current * v);
  r.stress.set(windOn(m, tau3));
  if (current) r.waterU.set(windOn(m, current3));
  if (tilt) for (let i = 0; i < r.C; i++) r.eta[i] = tilt * m.xCell[3 * i + 2];
  let change = 1, calls = 0;
  while (change > 1e-10 && calls < 300) {
    const before = Float64Array.from(r.drift.velocity);
    r.drift.solve(r.env());
    let moved = 0, top = 0;
    for (let e = 0; e < r.E; e++) { moved = Math.max(moved, Math.abs(r.drift.velocity[e] - before[e])); top = Math.max(top, Math.abs(r.drift.velocity[e])); }
    change = moved / (top + 1e-12); calls++;
  }
  const mu = RHO_I * h;
  let worst = 0, scale = 0;
  for (const e of edgesNear(m, c0, 3.5)) {
    const n = unit(m.nEdge, e), x = unit(m.xEdge, e), t = cross(x, n), f = m.fEdge[e];
    let Fn = dot(n, tau3), Ft = dot(t, tau3);
    const Un = dot(n, current3), Ut = dot(t, current3);
    Fn += f * mu * Ut; Ft -= f * mu * Un;
    if (tilt) {
      const k = dot([0, 0, 1], x), p = [-k * x[0], -k * x[1], 1 - k * x[2]];
      Fn -= mu * 9.81 * tilt * dot(n, p) / m.radius; Ft -= mu * 9.81 * tilt * dot(t, p) / m.radius;
    }
    const w = closedForm(Fn, Ft, f, mu, theta);
    worst = Math.max(worst, Math.abs(r.drift.velocity[e] - (Un + w.n))); scale = Math.max(scale, Math.hypot(w.n, w.t));
  }
  const v = cellVector(m, r.drift.velocity), vc = [v[3 * c0], v[3 * c0 + 1], v[3 * c0 + 2]];
  return { error: worst / scale, speed: scale, angle: Math.atan2(dot(vc, frame.north), dot(vc, frame.east)) / DEG };
}

test('free drift is the closed form of McPhee and Leppäranta in both hemispheres, with turning, a thin ice, a current and a sea-surface tilt', () => {
  const cases = [
    ['16 S', 16, -60, {}], ['16 S 20°', 16, -60, { theta: 20 }], ['16 N', 16, 60, {}], ['16 N 20°', 16, 60, { theta: 20 }],
    ['16 S thin', 16, -60, { h: 0.05 }], ['16 S current', 16, -60, { theta: 20, current: 0.2 }], ['16 S tilt', 16, -60, { wind: 0, tilt: 20 }],
    ['32 S 20°', 32, -60, { theta: 20 }], ['32 N', 32, 60, {}],
  ];
  const got = {};
  for (const [name, N, lat, options] of cases) {
    got[name] = freeDrift(N, lat, options);
    console.log(`free drift ${name}: ${(100 * got[name].error).toFixed(3)} % of ${got[name].speed.toFixed(3)} m/s, turned ${got[name].angle.toFixed(1)}° from the wind`);
    assert.ok(got[name].error < 0.01, `${name}: ${got[name].error}`);
  }
  assert.ok(got['16 N'].angle < -3 && got['16 N'].angle > -6 && got['16 S'].angle > 3 && got['16 S'].angle < 6, 'mass alone turns the ice right in the north and left in the south');
  near(got['16 N 20°'].angle, -24, 1.5, 'north with 20° of turning'); near(got['16 S 20°'].angle, 24, 1.5, 'south with 20° of turning');
  assert.ok(got['32 S 20°'].error < got['16 S 20°'].error, 'the error falls with the resolution');
});

test('the Coriolis force turns the ice by atan(fΔt) a step, to the right in the north and the left in the south, and does no work', () => {
  for (const latDeg of [60, -60]) for (const dt of [1350, 5400]) {
    const m = meshOf(16), r = rig(m, { rheology: 'free', waterDrag: 0, freePasses: 400 });
    r.dt = dt; r.thickness.fill(0.5); r.concentration.fill(0.8);
    const c0 = nearestCell(m, latDeg * DEG, 0), frame = localFrame(m, c0);
    r.drift.velocity.set(windOn(m, frame.east));
    const u0 = Float64Array.from(r.drift.velocity);
    r.drift.solve(r.env());
    const u1 = r.drift.velocity;
    let e0 = 0, e1 = 0, x01 = 0;
    for (let e = 0; e < r.E; e++) { const w = m.dcEdge[e] * m.dvEdge[e]; e0 += w * u0[e] * u0[e]; e1 += w * u1[e] * u1[e]; x01 += w * u0[e] * u1[e]; }
    near((e1 - x01) / e0, 0, 1e-9, `${latDeg}° ${dt} s: the new kinetic energy equals the overlap of old and new`);
    const a = cellVector(m, u0), b = cellVector(m, u1), pick = (v) => [dot([v[3 * c0], v[3 * c0 + 1], v[3 * c0 + 2]], frame.east), dot([v[3 * c0], v[3 * c0 + 1], v[3 * c0 + 2]], frame.north)];
    const [ae, an] = pick(a), [be, bn] = pick(b), turn = Math.atan2(bn * ae - be * an, be * ae + bn * an), f = m.fCell[c0];
    console.log(`Coriolis ${latDeg}° ${dt} s: turned ${(turn / DEG).toFixed(3)}° against ${(-Math.atan(f * dt) / DEG).toFixed(3)}°, speed ratio ${(Math.hypot(be, bn) / Math.hypot(ae, an)).toFixed(4)}`);
    assert.ok(turn * f < 0, 'right in the north, left in the south');
    near(turn / -Math.atan(f * dt), 1, 0.06, `${latDeg}° ${dt} s: the turning`);
    near(Math.hypot(be, bn) / Math.hypot(ae, an), 1 / Math.sqrt(1 + (f * dt) ** 2), 0.02, `${latDeg}° ${dt} s: the decay of the speed`);
  }
});

test('cavitating ice against a coast holds the wind with the gradient of its strength: the thickness ramp has the slope τ/P*', () => {
  const N = 16, tau = 0.05, P = ICE_DRIFT.strength, r = packRig(N, { rheology: 'cavitating' }, { wind: 0, dt: 3600 }), m = r.m;
  const sea = Array.from({ length: r.C }, (_, i) => i).filter((i) => r.cellOpen[i]), X = (i) => m.radius * m.xCell[3 * i];
  for (const i of sea) { r.thickness[i] = 4.5 + tau * X(i) / P; r.concentration[i] = 1; r.drift.sigma[i] = -2 * P * r.thickness[i]; }
  for (let e = 0; e < r.E; e++) if (r.edgeOcean[e]) r.stress[e] = tau * (X(m.cellsOnEdge[2 * e + 1]) - X(m.cellsOnEdge[2 * e])) / m.dcEdge[e];
  const start = Float64Array.from(r.thickness);
  let peak = 0, top = 0;
  for (let n = 0; n < 24; n++) { r.drift.advance(r.env()); top = r.drift.count.topSpeed; peak = Math.max(peak, top); }
  let change = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (const i of sea) { change = Math.max(change, Math.abs(r.thickness[i] - start[i]) / start[i]); sx += X(i); sy += r.thickness[i]; sxx += X(i) * X(i); sxy += X(i) * r.thickness[i]; }
  const slope = (sea.length * sxy - sx * sy) / (sea.length * sxx - sx * sx), free = Math.sqrt(tau / (1025 * 5.5e-3));
  const creep = change / (24 * 3600), regularization = ICE_DRIFT.minimumDeformation;
  console.log(`ramp held for 24 calls: top speed ${top.toFixed(4)} m/s at the end after ${peak.toFixed(4)} m/s at most, against ${free.toFixed(3)} m/s of free drift, thickness changed by ${(100 * change).toFixed(2)} %, ${(creep / regularization).toFixed(1)} Δ_min, slope ${slope.toExponential(3)} against ${(tau / P).toExponential(3)}`);
  assert.ok(creep < 10 * regularization && peak < free, `the ramp holds the wind, creeping at the regularization rate: ${creep} /s`);
  near(slope / (tau / P), 1, 0.05, 'slope of the ramp');

  const piled = packRig(N, { rheology: 'cavitating' }, { wind: 0.1, dt: 5400 }), cells = Array.from({ length: piled.C }, (_, i) => i).filter((i) => piled.cellOpen[i]);
  peak = 0;
  const x = cells.map((i) => m.radius * m.xCell[3 * i]), lo = Math.min(...x), hi = Math.max(...x);
  const mean = (pick) => { const set = cells.filter((_, j) => pick(x[j])); return set.reduce((s, i) => s + piled.thickness[i] * piled.concentration[i], 0) / set.length; };
  const report = [];
  for (let n = 1; n <= 3000; n++) {
    piled.drift.advance(piled.env()); peak = Math.max(peak, piled.drift.count.topSpeed);
    if (n % 750 === 0) report.push(`day ${(n * 5400 / 86400).toFixed(0)}: ${mean((v) => v > lo + 0.75 * (hi - lo)).toFixed(2)} m downwind, ${mean((v) => v < lo + 0.25 * (hi - lo)).toFixed(2)} m upwind, ${piled.drift.count.topSpeed.toFixed(3)} m/s`);
  }
  const down = mean((v) => v > lo + 0.75 * (hi - lo)), up = mean((v) => v < lo + 0.25 * (hi - lo)), thickest = Math.max(...cells.map((i) => piled.thickness[i]));
  console.log(`0.1 N/m² on a uniform meter of ice: ${report.join('; ')}; thickest ${thickest.toFixed(2)} m, ${piled.drift.count.closures} closures`);
  assert.ok(cells.every((i) => piled.concentration[i] <= 1 && piled.thickness[i] >= 0 && Number.isFinite(piled.thickness[i])));
  assert.ok(down > 1 && up < 0.5 && thickest > 1.5 && piled.drift.count.closures > 100, `the ice leaves the upwind coast and piles against the downwind one: ${down} ${up} ${thickest}`);
});

function bump(N) {
  const m = meshOf(N), r = rig(m, { rheology: 'frozen' }), steps = 1000 * N / 16, c3 = [Math.cos(20 * DEG) * Math.cos(30 * DEG), Math.cos(20 * DEG) * Math.sin(30 * DEG), Math.sin(20 * DEG)];
  r.dt = 3600 * 16 / N;
  for (let i = 0; i < r.C; i++) {
    const g = Math.exp(-((Math.acos(Math.min(1, dot(unit(m.xCell, i), c3))) * m.radius / 1.2e6) ** 2));
    r.concentration[i] = 0.25 + 0.6 * g; r.thickness[i] = 0.4 + 1.2 * g; r.snow[i] = 25 * g; r.surfaceT[i] = FREEZING_POINT - 8 * g;
  }
  r.drift.velocity.set(streamOn(m, 2 * Math.PI / (steps * r.dt)));
  const start ={ A: Float64Array.from(r.concentration), h: Float64Array.from(r.thickness), s: Float64Array.from(r.snow) };
  const total = () => {
    const out = { V: 0, S: 0, E: 0 };
    for (let i = 0; i < r.C; i++) {
      const a = r.concentration[i], Om = m.areaCell[i];
      out.V += Om * a * r.thickness[i]; out.S += Om * a * r.snow[i]; out.E += Om * a * r.seaIce.skinHeatCapacity * (r.surfaceT[i] - FREEZING_POINT);
    }
    return out;
  };
  const before = total();
  for (let n = 0; n < steps; n++) r.drift.transport(r.env());
  const after = total();
  let error = 0, area = 0, hLow = Infinity, hHigh = -Infinity, sHigh = -Infinity;
  for (let i = 0; i < r.C; i++) { error += m.areaCell[i] * (r.concentration[i] - start.A[i]) ** 2; area += m.areaCell[i]; hLow = Math.min(hLow, r.thickness[i]); hHigh = Math.max(hHigh, r.thickness[i]); sHigh = Math.max(sHigh, r.snow[i]); }
  return { r, before, after, rms: Math.sqrt(error / area), hLow, hHigh, sHigh, start };
}

test('transport carries area, thickness, snow and skin enthalpy once around the sphere, conserving volume, snow and enthalpy and creating no extrema', () => {
  const coarse = bump(16), fine = bump(32);
  for (const run of [coarse, fine]) {
    for (const key of ['V', 'S']) near(run.after[key] / run.before[key], 1, 1e-12, `integral ${key}`);
    near(run.after.E / run.before.E, 1, 1e-9, 'integral E, to the rounding of the skin temperature near 271 K');
    assert.ok(run.r.concentration.every((a) => a >= 0 && a <= 1), 'the cover stays between 0 and 1');
    assert.ok(run.hLow >= 0.4 - 0.02 && run.hHigh <= 1.6 + 0.02 && run.sHigh <= 25 + 0.5, `extrema: ${run.hLow} ${run.hHigh} ${run.sHigh}`);
  }
  console.log(`one revolution of a bump at a Courant number of 0.08: area error ${coarse.rms.toExponential(2)} at N = 16, ${fine.rms.toExponential(2)} at N = 32; thickness ${coarse.hLow.toFixed(3)}–${coarse.hHigh.toFixed(3)} and ${fine.hLow.toFixed(3)}–${fine.hHigh.toFixed(3)} m from 0.4–1.6`);
  assert.ok(fine.rms < 0.75 * coarse.rms, 'the error falls with the resolution');
});

test('cells filled past full cover keep their volume, snow mass and skin enthalpy, and the excess area becomes thickness', () => {
  const m = meshOf(16), c0 = nearestCell(m, -50 * DEG, 20 * DEG);
  const run = (maximum) => {
    const r = rig(m, { rheology: 'frozen', maximumConcentration: maximum });
    r.dt = 3600; r.thickness.fill(0.8); r.concentration.fill(0.9); r.snow.fill(12); r.surfaceT.fill(FREEZING_POINT - 6);
    for (let k = 0; k < m.nEdgesOnCell[c0]; k++) r.drift.velocity[m.edgesOnCell[m.maxEdges * c0 + k]] = -m.edgeSignOnCell[m.maxEdges * c0 + k] * 6;
    r.drift.transport(r.env());
    return r;
  };
  const closed = run(1), open = run(10);
  assert.ok(open.concentration[c0] > 1.05, `without the closure the cell holds ${open.concentration[c0]}`);
  assert.equal(closed.concentration[c0], 1);
  const volume = (r) => r.concentration[c0] * r.thickness[c0], snow = (r) => r.concentration[c0] * r.snow[c0], skin = (r) => r.concentration[c0] * (r.surfaceT[c0] - FREEZING_POINT);
  near(volume(closed) / volume(open), 1, 1e-13, 'volume'); near(snow(closed) / snow(open), 1, 1e-13, 'snow mass'); near(skin(closed) / skin(open), 1, 1e-13, 'skin enthalpy');
  assert.ok(closed.thickness[c0] > open.thickness[c0] && closed.drift.count.closures >= 1 && closed.drift.count.closedArea > 0);
  console.log(`a cell filled to ${open.concentration[c0].toFixed(3)} of its area closes to 1 and goes from ${open.thickness[c0].toFixed(3)} to ${closed.thickness[c0].toFixed(3)} m thick`);
});

function oceanRig(iceDrift, setup = () => {}, options = {}) {
  const m = meshOf(16), C = m.nCells, E = m.nEdges;
  const ocean = createOcean(m, { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34.5, iceDrift, ...options });
  const seaIce = createSeaIce(m, { heatCapacity: ocean.capacity });
  ocean.useSeaIce(seaIce.concentration, seaIce);
  const surfaceT = new Float64Array(C).fill(FREEZING_POINT + 2), ice = new Float64Array(C), o = { m, C, E, ocean, seaIce, surfaceT, ice, concentration: seaIce.concentration, calm: new Float64Array(E) };
  setup(o);
  ocean.initialize(surfaceT, ice, { concentration: seaIce.concentration });
  for (let i = 0; i < C; i++) if (!(ice[i] > 0)) surfaceT[i] = FREEZING_POINT + 2;
  o.env = (dt = 1350) => ({ surfaceT, ice, concentration: seaIce.concentration, stress: o.calm, waterU: ocean.u, eta: ocean.eta, dt });
  o.step = (stress = o.calm) => ocean.advance(surfaceT, ice, seaIce.oceanFlux, stress, 1350, seaIce.concentration);
  return o;
}
const icePatch = (lat, lon, degrees, { A = 0.8, h = 0.6, s = 15, T = FREEZING_POINT - 3 } = {}) => (o) => {
  for (const i of withinCap(o.m, lat, lon, degrees)) { o.ice[i] = h; o.concentration[i] = A; o.seaIce.snow[i] = s; o.surfaceT[i] = T; }
};

test('transport shifts the ocean\'s stored ice volume with it, so only freezing and melting reach the salt', () => {
  const setup = icePatch(-68, 30, 22), X = oceanRig({ rheology: 'frozen' }, setup), Y = oceanRig(false, setup), m = X.m, drift = X.ocean.iceDrift;
  for (const o of [X, Y]) for (let i = 0; i < o.C; i++) if (!(o.ice[i] > 0)) o.surfaceT[i] = FREEZING_POINT;
  X.step(); Y.step();
  assert.ok(bits(X.ocean.W, Y.ocean.W) && bits(X.ocean.Q, Y.ocean.Q), 'the same ocean without motion');
  const volume = (o) => Float64Array.from(o.ice, (h, i) => h * o.concentration[i]), before = volume(X);
  drift.velocity.set(rotationOn(m, 4 / m.radius));
  drift.transport(X.env());
  const after = volume(X);
  const shifted = X.ocean.restartArrays().previousIceVolume;
  let moved = 0, shift = 0, melted = 0, total = [0, 0], emptied = 0;
  for (let i = 0; i < X.C; i++) {
    moved = Math.max(moved, Math.abs(after[i] - before[i])); total[0] += m.areaCell[i] * before[i]; total[1] += m.areaCell[i] * after[i];
    if (after[i] > 0) shift = Math.max(shift, Math.abs(shifted[i] - after[i])); else { emptied++; melted = Math.min(melted, shifted[i] - after[i]); total[1] += m.areaCell[i] * (shifted[i] - after[i]); }
  }
  console.log(`transport moved up to ${moved.toFixed(4)} m of ice volume in a cell; the stored volume differs from the measured by ${shift.toExponential(2)} where ice remains, and ${emptied} cells hold none`);
  assert.ok(moved > 1e-3, 'ice moved');
  assert.ok(shift < 1e-12 && melted >= 0, 'previousVolume follows the transport, leaving what melts away to the salt');
  near(total[1] / total[0], 1, 1e-13, 'the ice volume of the globe, less what melted away');

  drift.velocity.fill(0);
  Y.ice.set(X.ice); Y.concentration.set(X.concentration); Y.seaIce.snow.set(X.seaIce.snow); Y.surfaceT.set(X.surfaceT);
  X.step(); Y.step();
  let worst = 0, scale = 0;
  for (let i = 0; i < X.C; i++) {
    if (!X.ocean.cellOcean[i]) continue;
    const c = 917 / 1000, h = X.ocean.h[i], gX = after[i] - shifted[i], pre = (X.ocean.W[i] + 5 * c * gX) / (1 + c * gX / h), expected = (pre / h - 5) * c * (shifted[i] - before[i]);
    worst = Math.max(worst, Math.abs((Y.ocean.W[i] - X.ocean.W[i]) - expected)); scale = Math.max(scale, Math.abs(expected));
  }
  console.log(`without the shift the salt of a cell would change by up to ${scale.toExponential(2)} (psu·m); the two oceans differ from that by ${worst.toExponential(2)}`);
  assert.ok(scale > 1e-2 && worst < 1e-9 * scale, `${worst} of ${scale}`);
  assert.ok(Y.ocean.W.some((w, i) => w !== X.ocean.W[i]), 'the control feels the spurious flux');
});

test('freezing in one cell, drift into another and melting there move salt from the one to the other by exactly the brine of the ice', () => {
  const m = meshOf(16), A = nearestCell(m, -62 * DEG, 0), B = m.cellsOnCell[m.maxEdges * A], edge = m.edgesOnCell[m.maxEdges * A];
  const o = oceanRig({ rheology: 'frozen' }, icePatch(-62, 0, 7, { A: 0.9, h: 1 })), { ocean } = o, drift = ocean.iceDrift;
  assert.ok(o.ice[A] > 0 && o.ice[B] > 0);
  o.step();
  const salt = () => Float64Array.from(ocean.W), W0 = salt(), brine = (i, dV) => (ocean.W[i] / ocean.h[i] - 5) * 917 / 1000 * dV;
  const grow = 0.2 * o.concentration[A]; o.ice[A] += 0.2;
  let W1 = salt(), expectedA = brine(A, grow);
  o.step();
  near(ocean.W[A] - W1[A], expectedA, 1e-9 * Math.abs(expectedA), 'brine of the freezing');
  drift.velocity[edge] = (m.cellsOnEdge[2 * edge] === A ? 1 : -1) * 3;
  W1 = salt(); const stored = ocean.restartArrays().previousIceVolume[B];
  o.step();
  drift.velocity.fill(0);
  assert.ok(Math.abs(ocean.restartArrays().previousIceVolume[B] - o.ice[B] * o.concentration[B]) < 1e-12 && Math.abs(stored - o.ice[B] * o.concentration[B]) > 1e-4, 'ice moved into B');
  const meltedAt = o.concentration[B] * 0.2; o.ice[B] -= 0.2;
  const W2 = salt(), expectedB = brine(B, -meltedAt);
  o.step();
  near(ocean.W[B] - W2[B], expectedB, 1e-5 * Math.abs(expectedB), 'freshwater of the melting');
  console.log(`brine ${expectedA.toExponential(3)} psu·m at A, freshening ${expectedB.toExponential(3)} at B`);
  assert.ok(W0.length > 0 && expectedA > 0 && expectedB < 0);
});

test('ice that drifts into water above freezing melts exactly the heat that the top layer holds, and the surface energy is conserved', () => {
  const m = meshOf(16), A = nearestCell(m, -60 * DEG, 0), B = m.cellsOnCell[m.maxEdges * A], edge = m.edgesOnCell[m.maxEdges * A];
  const o = oceanRig({ rheology: 'frozen' }, (r) => { r.ice[A] = 1; r.concentration[A] = 0.9; r.seaIce.snow[A] = 10; r.surfaceT[A] = FREEZING_POINT - 5; }), { ocean, seaIce } = o, drift = ocean.iceDrift;
  o.step();
  const energy = () => [A, B].reduce((s, i) => s + m.areaCell[i] * seaIce.energy(o.surfaceT[i], o.ice[i], seaIce.snow[i], o.concentration[i]), 0);
  assert.equal(o.ice[B], 0);
  const water = o.surfaceT[B], before = energy(), volume = m.areaCell[A] * o.ice[A] * o.concentration[A];
  drift.velocity[edge] = (m.cellsOnEdge[2 * edge] === A ? 1 : -1) * 8;
  drift.transport(o.env(3600));
  assert.equal(drift.arrived[B], 1); assert.equal(drift.waterT[B], water); assert.ok(o.concentration[B] > 0.01);
  near((energy() + m.areaCell[B] * seaIce.slabHeatCapacity * (drift.waterT[B] - FREEZING_POINT)) / before, 1, 1e-12, 'ice energy plus the water\'s held heat');
  near((m.areaCell[A] * o.ice[A] * o.concentration[A] + m.areaCell[B] * o.ice[B] * o.concentration[B]) / volume, 1, 1e-12, 'volume of the pair');
  ocean.readSurface(o.surfaceT, o.ice, o.concentration);
  assert.equal(ocean.T0[B], water); assert.equal(ocean.Q[B], ocean.h[B] * water);
  ocean.writeSurface(o.surfaceT, seaIce.oceanFlux, 3600);
  const handed = seaIce.oceanFlux[B] * 3600 / (ocean.rhoCp * ocean.h[B]);
  console.log(`B arrives at ${(water - FREEZING_POINT).toFixed(1)} K above freezing and hands ${handed.toFixed(6)} K of its top layer to the ice base`);
  near(handed, water - FREEZING_POINT, 1e-9, 'the heat the ice base receives');
});

test('the surface energy changes by exactly the surface flux while the ice drifts, closes up and melts away', () => {
  const m = meshOf(16), r = rig(m, { rheology: 'frozen' }), dt = 3600, patch = withinCap(m, -65, 40, 20);
  r.dt = dt;
  for (const i of patch) { r.thickness[i] = 1.2; r.concentration[i] = 0.7; r.snow[i] = 40; r.surfaceT[i] = FREEZING_POINT - 10; }
  const tiny = nearestCell(m, 40 * DEG, 100 * DEG);
  r.thickness[tiny] = 0.01; r.concentration[tiny] = 0.02; r.surfaceT[tiny] = FREEZING_POINT - 2;
  r.drift.velocity.set(rotationOn(m, 6 / m.radius));
  for (let k = 0; k < m.nEdgesOnCell[tiny]; k++) r.drift.velocity[m.edgesOnCell[m.maxEdges * tiny + k]] = m.edgeSignOnCell[m.maxEdges * tiny + k] * 8;
  const flux = new Float64Array(r.C), energy = () => { let s = 0; for (let i = 0; i < r.C; i++) s += m.areaCell[i] * r.seaIce.energy(r.surfaceT[i], r.thickness[i], r.snow[i], r.concentration[i]); return s; };
  let expected = energy(), cleared = 0, worst = 0, held = 0;
  const scale = patch.reduce((s, i) => s + m.areaCell[i] * r.seaIce.latent * 1.2, 0);
  for (let n = 0; n < 60; n++) {
    flux.fill(n < 20 ? -30 : 80);
    for (let i = 0; i < r.C; i++) r.seaIce.update(r.surfaceT, r.thickness, flux, i, dt);
    expected += dt * m.areaCell.reduce((s, a, i) => s + a * flux[i], 0);
    const had = r.thickness.map((h) => h > 0);
    r.drift.arrived.fill(0);
    r.drift.transport(r.env());
    for (let i = 0; i < r.C; i++) {
      if (had[i] && !(r.thickness[i] > 0)) cleared++;
      if (r.drift.arrived[i]) held += m.areaCell[i] * r.seaIce.slabHeatCapacity * (r.drift.waterT[i] - FREEZING_POINT);
    }
    worst = Math.max(worst, Math.abs(energy() + held - expected) / scale);
  }
  console.log(`60 steps with ${cleared} cells emptied by transport: the energy follows the flux to ${worst.toExponential(2)} of the patch's latent heat`);
  assert.ok(cleared >= 2 && worst < 1e-11, `${cleared} ${worst}`);
});

function history(iceDrift) {
  const o = oceanRig(iceDrift, (r) => { icePatch(-70, 0, 25)(r); icePatch(75, 120, 15, { A: 0.6, h: 0.9 })(r); }), flux = new Float64Array(o.C).fill(-30), wind = windOn(o.m, [0.08, 0, 0]), copy = (a) => Float64Array.from(a);
  for (let n = 0; n < 20; n++) {
    for (let i = 0; i < o.C; i++) o.seaIce.update(o.surfaceT, o.ice, flux, i, 1350);
    o.step(wind);
  }
  const saved = o.ocean.restartArrays();
  return { stress: copy(o.ocean.stress), Q: copy(o.ocean.Q), W: copy(o.ocean.W), h: copy(o.ocean.h), u: copy(o.ocean.u), previousIceVolume: saved.previousIceVolume, oceanFlux: copy(o.seaIce.oceanFlux), surfaceT: copy(o.surfaceT), ice: copy(o.ice), concentration: copy(o.concentration), snow: copy(o.seaIce.snow) };
}

test('iceDrift false, absent and frozen reproduce the ocean\'s stress, ice volume, heat, salt and ocean flux bit for bit over 20 calls', () => {
  const base = history(false), absent = history(undefined), frozen = history({ rheology: 'frozen' });
  assert.ok(base.ice.some((h) => h > 0) && base.stress.some((s) => s !== 0));
  for (const key of Object.keys(base)) { assert.ok(bits(base[key], absent[key]), `absent: ${key}`); assert.ok(bits(base[key], frozen[key]), `frozen: ${key}`); }
});

test('the option validates, a running ocean keeps A at most 1 and V non-negative, and the ice velocity rides the restart arrays', () => {
  assert.throws(() => iceDriftOptions({ rheology: 'vp' }), /rheology/);
  assert.throws(() => iceDriftOptions({ hemisphere: 'equator' }), /hemisphere/);
  assert.throws(() => iceDriftOptions({ strenth: 1 }), /no option/);
  assert.throws(() => iceDriftOptions({ iterations: 0 }), /iterations/);
  assert.deepEqual(iceDriftOptions(true), ICE_DRIFT);
  const setup = icePatch(-62, 10, 25, { A: 0.9, h: 0.5 }), wind = windOn(meshOf(16), [0.1, 0, 0]);
  const run = (hemisphere) => { const o = oceanRig({ rheology: 'cavitating', hemisphere }, setup); for (let n = 0; n < 6; n++) o.step(wind); return o; };
  const south = run('south'), north = run('north');
  assert.ok(south.ocean.iceDrift.velocity.some((u) => u !== 0) && north.ocean.iceDrift.velocity.every((u) => u === 0), 'the hemisphere switch keeps the other hemisphere\'s ice still');
  assert.ok(south.concentration.every((a) => a >= 0 && a <= 1) && south.ice.every((h) => h >= 0 && Number.isFinite(h)));
  const d = south.ocean.diagnostics();
  console.log(`south: top speed ${d.iceSpeed.toFixed(4)} m/s, ${d.iceSpeedClamps} clamps, ${d.iceLimiterFallbacks} fallbacks, ${d.iceClosures} closures, ice transport across 60S ${d.iceTransport60S.toExponential(2)} m³/s`);
  assert.ok(d.iceSpeed > 0.01 && d.iceSpeed <= 2);
  const saved = { ...south.ocean.serialize(), ...south.ocean.restartArrays(), flux: new Float64Array(south.C) };
  assert.ok(saved.iceVelocity && saved.iceVelocity.some((u) => u !== 0));
  const again = oceanRig({ rheology: 'cavitating' }, setup);
  again.ocean.load(saved, again.surfaceT, again.ice, again.concentration);
  assert.ok(bits(again.ocean.iceDrift.velocity, saved.iceVelocity), 'the velocity loads');
  const { iceVelocity, ...older } = saved;
  const plain = oceanRig({ rheology: 'cavitating' }, setup);
  plain.ocean.load(older, plain.surfaceT, plain.ice, plain.concentration);
  assert.ok(plain.ocean.iceDrift.velocity.every((u) => u === 0), 'a state without it loads unchanged');
  const off = oceanRig(false, setup);
  off.ocean.load(saved, off.surfaceT, off.ice, off.concentration);
  assert.equal('iceVelocity' in off.ocean.restartArrays(), false);
});

test('a cavitating drift restarted from its arrays steps on bit for bit: velocity, the internal stress σ1, the ice and the ocean', () => {
  const setup = icePatch(-62, 10, 25, { A: 0.9, h: 0.5 }), wind = windOn(meshOf(16), [0.1, 0, 0]), drift = { rheology: 'cavitating', hemisphere: 'both' };
  const whole = oceanRig(drift, setup), first = oceanRig(drift, setup);
  for (let n = 0; n < 12; n++) whole.step(wind);
  for (let n = 0; n < 6; n++) first.step(wind);
  const saved = { ...first.ocean.serialize(), ...first.ocean.restartArrays(), flux: Float64Array.from(first.seaIce.oceanFlux) };
  assert.ok(saved.iceSigma && saved.iceSigma.some((x) => x !== 0), 'the pack carries an internal stress');
  const resume = (state) => {
    const o = oceanRig(drift, setup);
    for (const [name, values] of Object.entries(first.seaIce)) if (values instanceof Float64Array && values.length === o.C) o.seaIce[name].set(values);
    o.surfaceT.set(first.surfaceT); o.ice.set(first.ice);
    o.ocean.load(state, o.surfaceT, o.ice, o.concentration);
    return o;
  };
  const again = resume(saved);
  for (let n = 0; n < 6; n++) again.step(wind);
  for (const [name, a, b] of [['velocity', whole.ocean.iceDrift.velocity, again.ocean.iceDrift.velocity], ['σ1', whole.ocean.iceDrift.sigma, again.ocean.iceDrift.sigma], ['ice', whole.ice, again.ice], ['concentration', whole.concentration, again.concentration], ['Q', whole.ocean.Q, again.ocean.Q], ['W', whole.ocean.W, again.ocean.W], ['u', whole.ocean.u, again.ocean.u]]) assert.ok(bits(a, b), name);
  const { iceSigma, ...older } = saved, without = resume(older);
  assert.ok(without.ocean.iceDrift.sigma.every((x) => x === 0), 'a state without σ1 starts it at zero');
  for (let n = 0; n < 6; n++) without.step(wind);
  const apart = (a, b) => a.reduce((m, x, n) => Math.max(m, Math.abs(x - b[n])), 0);
  console.log(`restarted without σ1 (up to ${apart(saved.iceSigma, new Float64Array(saved.iceSigma.length)).toExponential(2)} N/m), six steps on: ice velocity ${apart(whole.ocean.iceDrift.velocity, without.ocean.iceDrift.velocity).toExponential(2)} m/s apart, thickness ${apart(whole.ice, without.ice).toExponential(2)} m, concentration ${apart(whole.concentration, without.concentration).toExponential(2)}`);
});

test('with its defaults the solver is within 3 percent of the converged solution on a compact pack from rest, and below its stability limit it is not', () => {
  for (const [N, dt] of [[16, 5400], [16, 1350], [32, 2700]]) {
    const reference = packSolve(N, dt, { alpha: 40, beta: 40, iterations: 3000 }), run = packSolve(N, dt, {});
    const error = packDifference(run, reference), { alpha, beta, passes } = run.r.drift.count;
    console.log(`pack N = ${N}, Δt = ${dt} s: α = ${alpha}, β = ${beta}, ${passes} iterations, ${(100 * error).toFixed(2)} % from the converged speeds`);
    assert.ok(error < 0.03, `${error}`);
  }
  const unstable = packDifference(packSolve(32, 2700, { alpha: 10, beta: 10, iterations: 400 }), packSolve(32, 2700, { alpha: 40, beta: 40, iterations: 3000 }));
  console.log(`α = β = 10 at N = 32, Δt = 2700 s: ${(100 * unstable).toFixed(2)} % after 400 iterations`);
  assert.ok(unstable > 0.01, 'α β of 100 is under the limit 3 P* Δt / (ρ_i Δ_min d²) = 420 there');
});
