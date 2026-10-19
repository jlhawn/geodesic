import { test } from 'node:test';
import assert from 'node:assert/strict';
import * as original from './helpers/layered.mjs';
import {
  caseMesh, zonalWindOnEdges, ekman, lensCase, bumpCase, runCase, internalWaveSpeed, ekmanLayer, meanWithin, EKMAN_REACH, PUMP_DEPTH, RHO, G, DEG, DAY,
} from '../explain/oceanCases.module.js';
import { ekmanSpiral, fitViscosity, coriolis, columnDays } from '../explain/ekmanTheory.module.js';
import { HYBRID_LAYERS } from '../js/ocean/column.module.js';

const mesh = caseMesh(8), C = mesh.nCells;
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const between = (v, lo, hi, what) => assert.ok(v >= lo && v <= hi, `${what}: ${v} not within ${lo} to ${hi}`);
const sameArray = (a, b, what) => { assert.equal(a.length, b.length, `${what} length`); for (let n = 0; n < a.length; n++) assert.equal(a[n], b[n], `${what}[${n}]`); };
const angles = (profile, upTo) => Array.from({ length: upTo }, (_, k) => Math.atan2(-profile.north[k], profile.east[k]) / DEG);

test('the case mesh and wind helper give what test/helpers/layered.mjs gives', () => {
  const speed = (lat) => 8 * Math.cos(3 * lat);
  sameArray(zonalWindOnEdges(mesh, speed), original.zonalWindOnEdges(original.mesh, speed), 'zonal wind');
  for (const name of ['xCell', 'xEdge', 'areaCell', 'dcEdge', 'fEdge', 'latCell']) sameArray(mesh[name], original.mesh[name], name);
});

test('the Ekman layer is read against the water below it, and the stack’s gravity-wave speed follows its stratification', () => {
  const tau = 0.1, f = coriolis(45), spiral = ekmanSpiral({ tau, f, viscosity: 0.02 });
  const thick = Float64Array.from({ length: 80 }, () => 5), east = new Float64Array(80), north = new Float64Array(80);
  for (let k = 0; k < 80; k++) { const [u, v] = spiral.at(-(k + 0.5) * 5); east[k] = 0.004 + u; north[k] = -0.002 + v; }
  const layer = ekmanLayer(thick, east, north, 40);
  close(layer.depth, EKMAN_REACH * 40, 1e-12, 'depth the transport is summed over');
  close(layer.east, 0.004, 5e-4, 'the reference current east'); close(layer.north, -0.002, 5e-4, 'the reference current north');
  close(layer.transportNorth / spiral.transport[1], 1, 0.06, 'the truncated spiral carries nearly all of τ/(ρ₀f)');
  close(layer.transportEast, 0, 0.05 * Math.abs(spiral.transport[1]), 'and no transport along the wind');
  close(meanWithin([10, 10, 10], [1, 2, 3], 15), 4 / 3, 1e-12, 'mean within 15 m of three 10 m layers');

  const c2 = internalWaveSpeed([100, 300], [1024, 1026]);
  close(c2, Math.sqrt(G * (2 / RHO) * 100 * 300 / 400), 1e-9, 'two layers: √(g′ h₁h₂/(h₁+h₂))');
  const n = 400, uniform = internalWaveSpeed(new Array(n).fill(10), Array.from({ length: n }, (_, k) => 1024 + 1e-3 * (k + 0.5) * 10));
  close(uniform, Math.sqrt(G / RHO * 1e-3) * 4000 / Math.PI, 0.01 * uniform, 'uniform stratification: N D/π');
  close(internalWaveSpeed([100, 300, 0.5, 0.5], [1024, 1026, 1026.00001, 1026.00002]), c2, 2e-3 * c2, 'a meter of vanished layers of the same water joins the layer above it');
});

test('a band of westerlies at 45°N drives a spiral of currents whose water carries τ/(ρ₀f) at right angles to the wind, to the right, and to the left at 45°S', () => {
  const north = ekman({ latitude: 45, days: columnDays(45) }), south = ekman({ latitude: -45, days: columnDays(-45) });
  const show = (r) => `${r.transport.magnitude.toFixed(4)} m²/s at ${r.transport.angle.toFixed(1)}° (${r.transport.ratio.toFixed(3)} of τ/(ρ₀f) = ${r.theory.toFixed(4)}); top layer ${(100 * r.top.magnitude).toFixed(2)} cm/s at ${r.top.angle.toFixed(1)}°; mixed layer ${(100 * r.mixed.magnitude).toFixed(2)} cm/s at ${r.mixed.angle.toFixed(1)}°, ${r.mixedDepth.toFixed(1)} m deep, boundary layer ${r.boundaryDepth.toFixed(1)} m, summed over ${r.reach.toFixed(1)} m; mean viscosity ${r.meanViscosity.toFixed(4)} m²/s`;
  console.log(`45°N: ${show(north)}; inertial period ${(north.inertialPeriod / 3600).toFixed(2)} h, ${north.ringCells} cells`);
  console.log(`45°S: ${show(south)}`);
  close(north.inertialPeriod / 3600, 16.92, 0.01, 'inertial period at 45° in hours');
  close(north.theory, 0.1 / (RHO * 2 * mesh.omega * Math.sin(45 * DEG)), 0.05, 'the ring’s τ/(ρ₀f) is near the one at 45°');
  assert.ok(north.settled && north.averagedOver >= 2 * north.inertialPeriod - 1350 && north.averagedOver <= 2 * north.inertialPeriod + 1350);
  close(north.transport.ratio, 1, 0.05, 'Ekman transport against τ/(ρ₀f) at 45°N');
  between(north.transport.angle, 85, 92, 'the transport’s angle to the right of the wind');
  close(south.transport.ratio, north.transport.ratio, 0.01, 'the same transport at 45°S');
  close(south.transport.angle, -north.transport.angle, 0.5, 'to the left of the wind at 45°S');
  between(north.top.angle, 45, 60, 'the top layer drifts about 45° to the right'); between(100 * north.top.magnitude, 5, 6.2, 'top layer cm/s');
  between(north.mixed.angle, 80, 92, 'the mixed layer’s drift is at right angles'); between(100 * north.mixed.magnitude, 2.3, 2.9, 'mixed layer cm/s');
  close(north.mixed.magnitude * north.mixedDepth, north.transport.magnitude, 0.1 * north.transport.magnitude, 'the mixed layer carries the transport');
  between(north.mixedDepth, 32, 38, 'mixed layer depth'); between(north.boundaryDepth, 35, 41, 'boundary layer depth');
  assert.ok(north.mixedDepth < north.boundaryDepth && north.boundaryDepth * EKMAN_REACH === north.reach);
  between(north.meanViscosity, 0.011, 0.017, 'viscosity averaged over the boundary layer');

  const p = north.profile, a = angles(p, 12);
  assert.equal(p.center.length, HYBRID_LAYERS.K);
  between(a[0], 45, 60, 'the top layer’s angle'); assert.ok(a[7] > a[0] + 90, `the current turns more than 90° to the right between the top layer and ${p.center[7].toFixed(0)} m: ${a[0]} to ${a[7]}`);
  for (let k = 1; k <= 7; k++) assert.ok(a[k] > a[k - 1], `the current turns further right at layer ${k}`);
  assert.ok(Math.hypot(p.east[9], p.north[9]) < 0.05 * Math.hypot(p.east[0], p.north[0]), 'the current has died away by 40 m');
  const nu = Array.from(p.viscosity);
  between(Math.max(...nu.slice(0, 9)), 0.02, 0.03, 'the largest viscosity within the boundary layer'); assert.ok(nu[12] < 5e-4, `viscosity below the boundary layer: ${nu[12]}`);

  const n = north.series.time.length;
  assert.equal(n, columnDays(45) * 64 + 1);
  for (const name of Object.keys(north.series)) assert.equal(north.series[name].length, n, name);
  let fastest = 0;
  for (let s = 0; s < n; s++) fastest = Math.max(fastest, Math.hypot(north.series.mixedEast[s] - north.series.deepEast[s], north.series.mixedNorth[s] - north.series.deepNorth[s]));
  const deepFit = (r) => { const q = fitViscosity({ tau: r.stressMean, f: r.fMean, center: r.profile.center, thickness: r.profile.thickness, east: r.profile.east, north: r.profile.north, reach: r.reach }); return { ...q, surface: Math.abs(ekmanSpiral({ tau: r.stressMean, f: r.fMean, viscosity: q.viscosity }).surface) }; };
  const fit45 = deepFit(north);
  console.log(`transport ${north.transport.magnitude.toFixed(5)} m²/s against ${north.theory.toFixed(5)}; fitted viscosity ${fit45.viscosity.toFixed(4)} m²/s, misfit ${(100 * fit45.rms).toFixed(2)} cm/s of ${(100 * fit45.surface).toFixed(2)}, ${(100 * fit45.rms / fit45.surface).toFixed(1)} percent`);
  between(fit45.rms / fit45.surface, 0.04, 0.09, 'the spiral fits the 45°N layers to within a few percent of its surface current');
  let loops = 0;
  for (let s = 0; s < north.series.time.length; s++) if (north.series.time[s] > north.series.time[north.series.time.length - 1] - 17 * 3600) loops = Math.max(loops, Math.hypot(north.series.mixedEast[s] - north.series.deepEast[s] - north.mixed.east, north.series.mixedNorth[s] - north.series.deepNorth[s] - north.mixed.north));
  console.log(`on day 4 the mixed layer’s loops reach ${(100 * loops).toFixed(2)} cm/s around its drift of ${(100 * north.mixed.magnitude).toFixed(2)} cm/s`);
  between(100 * loops, 0.8, 1.3, 'the inertial loops on day 4, cm/s');
  console.log(`the mixed layer’s speed against the water below peaks at ${(100 * fastest).toFixed(2)} cm/s while it circles and drifts at ${(100 * north.mixed.magnitude).toFixed(2)} cm/s`);
  assert.ok(fastest > 2 * north.mixed.magnitude, 'the mixed layer overshoots in an inertial circle before it settles');
  console.log(`the mixed layer deepens from ${north.series.mixedDepth[0].toFixed(1)} m to ${north.series.mixedDepth[n - 1].toFixed(1)} m and the boundary layer from ${north.series.boundaryDepth[0].toFixed(1)} to ${north.series.boundaryDepth[n - 1].toFixed(1)} m; layer angles ${a.slice(0, 9).map((v) => v.toFixed(0)).join(', ')}° at ${Array.from(p.center.slice(0, 9), (z) => z.toFixed(0)).join(', ')} m; speeds ${Array.from({ length: 9 }, (_, k) => (100 * Math.hypot(p.east[k], p.north[k])).toFixed(2)).join(', ')} cm/s; viscosity ${nu.slice(0, 10).map((v) => v.toFixed(4)).join(', ')} m²/s`);
  between(north.series.mixedDepth[0], 4.5, 6, 'the analytic start’s mixed layer'); between(north.series.mixedDepth[n - 1], 34, 38, 'the mixed layer after four days');
  between(north.series.boundaryDepth[0], 0, 0.01, 'no boundary layer before the first step'); between(north.series.boundaryDepth[n - 1], 37, 41, 'the boundary layer after four days');
  between(a[4], 88, 100, 'the current at 16 m points at right angles to the wind'); between(Math.hypot(p.east[7], p.north[7]) / Math.hypot(p.east[0], p.north[0]), 0.2, 0.35, 'a quarter of the top layer’s speed remains at 29 m');
});

test('toward the equator the Ekman balance weakens: the transport falls short of τ/(ρ₀f) and turns less than a right angle', () => {
  const mid = ekman({ latitude: 15, days: columnDays(15) }), low = ekman({ latitude: 5, days: columnDays(5) });
  for (const [name, r] of [['15°N', mid], ['5°N', low]]) console.log(`${name}: ${r.transport.magnitude.toFixed(3)} m²/s at ${r.transport.angle.toFixed(1)}° (${r.transport.ratio.toFixed(3)} of τ/(ρ₀f) = ${r.theory.toFixed(3)}); top layer ${r.top.angle.toFixed(1)}°; mixed ${r.mixedDepth.toFixed(1)} m, boundary ${r.boundaryDepth.toFixed(1)} m; cells ${r.ringCells}`);
  const share = (r) => { const q = fitViscosity({ tau: r.stressMean, f: r.fMean, center: r.profile.center, thickness: r.profile.thickness, east: r.profile.east, north: r.profile.north, reach: r.reach }); return q.rms / Math.abs(ekmanSpiral({ tau: r.stressMean, f: r.fMean, viscosity: q.viscosity }).surface); };
  console.log(`the best constant viscosity misses the layers by ${(100 * share(mid)).toFixed(1)} percent of its surface current at 15°N and ${(100 * share(low)).toFixed(1)} percent at 5°N; boundary layers ${mid.boundaryDepth.toFixed(2)} and ${low.boundaryDepth.toFixed(2)} m`);
  between(share(mid), 0.1, 0.18, 'misfit at 15°N'); between(share(low), 0.15, 0.25, 'misfit at 5°N');
  between(mid.boundaryDepth, 38, 43, 'boundary layer at 15°N'); between(low.boundaryDepth, 48, 53, 'boundary layer at 5°N');
  between(mid.transport.ratio, 0.95, 1.2, 'Ekman transport against τ/(ρ₀f) at 15°N'); between(mid.transport.angle, 85, 105, 'angle at 15°N');
  between(low.transport.ratio, 0.75, 0.95, 'Ekman transport against τ/(ρ₀f) at 5°N'); between(low.transport.angle, 55, 78, 'angle at 5°N');
  const lats = Array.from(mesh.latCell, (v) => v / DEG), ring = [];
  for (let i = 0; i < C; i++) if (Math.abs(lats[i] - 5) < 0.5 * (mesh.dcEdge.reduce((s, d) => s + d, 0) / mesh.nEdges) / mesh.radius / DEG) ring.push(lats[i]);
  between(Math.min(...ring), 0, 3, 'the 5°N ring reaches nearly to the equator'); between(Math.max(...ring), 7, 10, 'and up to 8°N');
});

test('over the planet the trades lift water at the equator and sink it near 15° and 30°, and the westerlies’ Ekman transport is τ/(ρ₀f)', () => {
  const trades = ekman({ view: 'planet', wind: 'trades' }), westerlies = ekman({ view: 'planet', wind: 'westerlies' });
  const bin = (r, lat) => r.zonal.latitude.findIndex((l) => l === lat);
  const w = (r, lat) => r.zonal.upwelling[bin(r, lat)] * DAY;
  const band = (r, lat) => { let sum = 0, n = 0; for (let i = 0; i < C; i++) if (Math.abs(mesh.latCell[i] / DEG - lat) < 2.5) { sum += r.upwelling[i]; n++; } return sum / n * DAY; };
  console.log(`trades, as the readout averages them: ${band(trades, 15).toFixed(3)} m/day at 15°, ${band(trades, 25).toFixed(3)} at 25°, ${band(trades, 30).toFixed(3)} at 30°; Ekman transport against τ/(ρ₀f) at 42.5° and 47.5°: ${[42.5, 47.5].map((lat) => (trades.zonal.ekmanNorth[bin(trades, lat)] / trades.zonal.theory[bin(trades, lat)]).toFixed(3)).join(', ')}`);
  between(band(trades, 15), -0.4, -0.3, 'sinking at 15°N, m/day'); between(band(trades, 25), 0.05, 0.13, 'rising at 25°N, m/day'); between(band(trades, 30), -0.2, -0.12, 'sinking at 30°N, m/day');
  for (const lat of [42.5, 47.5]) between(trades.zonal.ekmanNorth[bin(trades, lat)] / trades.zonal.theory[bin(trades, lat)], 0.9, 1.0, `Ekman transport under the trades and westerlies at ${lat}°`);
  console.log(`trades: ${PUMP_DEPTH} m down the equator rises ${(trades.equatorUpwelling * DAY).toFixed(3)} m/day; 12.5° ${w(trades, 12.5).toFixed(3)}, 17.5° ${w(trades, 17.5).toFixed(3)}, 22.5° ${w(trades, 22.5).toFixed(3)}, 27.5° ${w(trades, 27.5).toFixed(3)}, 32.5° ${w(trades, 32.5).toFixed(3)}; westerlies only: ${(westerlies.equatorUpwelling * DAY).toFixed(3)} m/day`);
  for (const lat of [42.5, 47.5, -42.5, -47.5]) {
    const b = bin(westerlies, lat), got = westerlies.zonal.ekmanNorth[b], want = westerlies.zonal.theory[b];
    console.log(`westerlies at ${lat}°: ${got.toFixed(3)} m²/s north against ${want.toFixed(3)}`);
    close(got / want, 1, 0.06, `zonal-mean Ekman transport at ${lat}°`);
  }
  between(trades.equatorUpwelling * DAY, 0.7, 1.1, 'equatorial upwelling under the trades, m/day');
  assert.ok(w(trades, 17.5) < -0.15 && w(trades, 12.5) < -0.25, 'the trades’ drift piles up and sinks near 15°');
  assert.ok(w(trades, 22.5) > 0.05, 'and the water rises a little between 20° and 25°');
  assert.ok(w(trades, 27.5) < -0.05 && w(trades, 32.5) < -0.05, 'and sinks again near 30°');
  assert.ok(Math.abs(westerlies.equatorUpwelling * DAY) < 0.02, 'no wind at the equator, no upwelling');
  close(trades.averagedOver, 2 * DAY, 1350, 'the window');
  assert.ok(Math.max(...trades.boundaryDepth) < PUMP_DEPTH / 2, 'the boundary layer is far above the depth the pumping is read at');
  for (const name of ['ekmanEast', 'ekmanNorth', 'upwelling', 'eta', 'mixedDepth', 'boundaryDepth', 'stress']) assert.equal(trades[name].length, C, name);
});

test('a lens of one layer 60 m thick pushes the layers beneath it down and raises the sea surface by Δρ/ρ₀ of that, less the bottom pressure the barotropic waves leave, with the surface current geostrophic and the deep water slow', () => {
  const sim = lensCase({ days: 30 }), layout = sim.layout(), { ocean } = sim, Cn = sim.mesh.nCells;
  assert.equal(layout.densities[layout.lens], 1025.25); assert.equal(layout.densities[layout.donor], 1026.77); assert.equal(layout.fixedDepth, 33);
  let top = 0;
  for (let k = 0; k < layout.fixedDepth; k++) { close(ocean.h[k * Cn + sim.centre], HYBRID_LAYERS.minimum[k], 1e-9, `fixed-depth layer ${k} at the lens`); top += ocean.h[k * Cn + sim.centre]; }
  close(top, 263.9, 0.1, 'the fixed-depth layers reach 264 m');
  const far = layout.slice.cells[0], at = (k, i) => ocean.h[k * Cn + i];
  const shift = (k, i) => { let a = 0, b = 0; for (let j = 0; j <= k; j++) { a += at(j, i); b += at(j, far); } return a - b; };
  const start = sim.frame();
  console.log(`start: lens layer ${at(layout.lens, sim.centre).toFixed(1)} m thick at the center against ${at(layout.lens, far).toFixed(1)} m, base pushed ${start.depression.toFixed(2)} m; bottom water ${at(layout.donor, sim.centre).toFixed(1)} against ${at(layout.donor, far).toFixed(1)}; layers pushed down: ${[layout.lens + 1, 40, layout.donor].map((k) => shift(k - 1, sim.centre).toFixed(2)).join(', ')} m`);
  close(start.depression, 60, 1, 'the lens is 60 m thick at the center'); close(at(layout.lens, sim.centre) - at(layout.lens, far), start.depression, 0.2, 'the lens layer is that much thicker');
  for (let k = layout.lens + 1; k < layout.donor; k++) close(shift(k, sim.centre), start.depression, 0.2, `interface ${k + 1} is pushed down by the lens`);
  close(shift(layout.donor, sim.centre), 0, 1e-6, 'and the column holds the same water');
  for (let k = 0; k < layout.fixedDepth; k++) close(shift(k, sim.centre), 0, 1e-9, `fixed-depth layer ${k} keeps its depth`);
  close(start.rise, 0, 1e-9, 'the sea surface starts flat');

  const flat = lensCase({ depth: 0 }), settled = Float64Array.from(flat.ocean.h);
  flat.ocean.regrid(Infinity);
  let moved = 0;
  for (let n = 0; n < settled.length; n++) moved = Math.max(moved, Math.abs(flat.ocean.h[n] - settled[n]));
  console.log(`a further regrid of the start moves no layer by more than ${(100 * moved).toFixed(1)} cm`);
  assert.ok(moved < 0.1, 'the start is already on the layers’ targets');

  const speed = sim.waveSpeed, radius = layout.deformationRadius;
  console.log(`first internal gravity wave ${speed.toFixed(3)} m/s, deformation radius ${(radius / 1e3).toFixed(1)} km at 45°`);
  between(speed, 2.6, 2.75, 'wave speed along the stratification'); close(radius, speed / (2 * mesh.omega * Math.sin(45 * DEG)), 1, 'the deformation radius is c/f at 45°');
  assert.ok(radius < 0.03 * 1.5e6, 'far smaller than the lens');
  const f30 = 2 * mesh.omega * Math.sin(30 * DEG), beta30 = 2 * mesh.omega * Math.cos(30 * DEG) / mesh.radius, years = 80 * DEG * mesh.radius * Math.cos(30 * DEG) / (beta30 * (speed / f30) ** 2) / (365 * DAY);
  console.log(`a first-mode Rossby wave at 30°N, β(c/f)², takes ${years.toFixed(1)} years to cross 80° of longitude`);
  between(years, 8.5, 10, 'years for a Rossby wave to cross the basin');

  const frames = {};
  for (const day of [1, 5, 10, 20, 30]) { runCase(sim, day); frames[day] = sim.frame(); }
  for (const f of Object.values(frames)) console.log(`day ${f.day}: surface ${(100 * f.rise).toFixed(2)} cm up over a lens base ${f.depression.toFixed(2)} m down (expected ${(100 * f.expected).toFixed(2)} cm, bottom pressure ${(100 * f.bottomPressure).toFixed(2)} cm); f·u against g∂η/∂n ${f.geostrophic.ratio.toFixed(3)}, correlation ${f.geostrophic.correlation.toFixed(4)}; speeds top ${(100 * f.layerSpeeds[0]).toFixed(3)}, lens layer ${(100 * f.layerSpeeds[layout.lens]).toFixed(3)}, bottom water ${(100 * f.layerSpeeds[layout.donor]).toFixed(3)} cm/s`);
  for (const day of [5, 10, 20, 30]) {
    const f = frames[day];
    between(f.rise / f.expected, 0.78, 0.95, `day ${day}: the rise against Δρ/ρ₀ times the lens`);
    close(f.rise, f.expected + f.bottomPressure, 0.004, `day ${day}: the rise is the expectation plus the bottom pressure`);
    between(f.bottomPressure / -f.expected, 0.05, 0.25, `day ${day}: the bottom pressure keeps a small part of the mass the lens removes`);
    close(f.geostrophic.ratio, 1, 0.05, `day ${day}: geostrophic ratio`);
    assert.ok(f.geostrophic.correlation > 0.995, `day ${day}: correlation ${f.geostrophic.correlation}`);
    assert.ok(f.layerSpeeds[layout.donor] < 0.25 * f.layerSpeeds[0], `day ${day}: the bottom water is slow`);
    assert.ok(f.layerSpeeds[layout.lens + 2] < f.layerSpeeds[layout.lens], `day ${day}: the layers beneath the lens move more slowly than it does`);
  }
  const speeds10 = Array.from(frames[10].layerSpeeds).slice(0, layout.fixedDepth);
  assert.ok(Math.max(...speeds10) < 1.02 * Math.min(...speeds10), 'the fixed-depth layers move alike');
  console.log(`lens layer density ${frames[10].lensDensity.toFixed(3)}, bottom water ${frames[10].donorDensity.toFixed(3)}, difference ${(frames[10].donorDensity - frames[10].lensDensity).toFixed(3)} kg/m³; bottom pressure ${(100 * -frames[10].bottomPressure / frames[10].expected).toFixed(1)} percent of the mirror on day 10; layer speeds beneath the lens ${Array.from(frames[10].layerSpeeds).slice(layout.lens, layout.lens + 6).map((v) => (100 * v).toFixed(3)).join(', ')} cm/s`);
  assert.ok(frames[10].depression > 58.5 && frames[10].depression < 60, `the lens is still ${frames[10].depression} m deep on day 10`);
  const lifted = start.depression - frames[30].depression;
  console.log(`the eddy transport lifts the lens base by ${lifted.toFixed(2)} m in 30 days`);
  between(lifted, 0.1, 1, 'the lens survives 30 days, lifted only slowly');
  between(frames[30].rise / frames[30].expected, 0.78, 0.95, 'and the surface stays up');
  assert.equal(frames[10].slice.distance.length, layout.slice.cells.length);
  assert.equal(frames[10].slice.interfaces.length, layout.slice.cells.length * layout.layers);
  assert.equal(ocean.diagnostics().oceanLimited, 0);
});

test('a 1 m bump of sea surface races away at about √(gD) and leaves a geostrophic eddy that every layer circles alike, keeping the share of its mass that the lens’s bottom pressure keeps', () => {
  const sim = bumpCase({ days: 1 }), depth = 4000, layout = sim.layout();
  const Cn = sim.mesh.nCells, farCell = layout.slice.cells[0];
  assert.equal(layout.bumped, 48);
  close(sim.ocean.h[layout.bumped * Cn + sim.centre] - sim.ocean.h[layout.bumped * Cn + farCell], 1, 0.02, 'the deepest layer holding water is 1 m thicker at the center');
  runCase(sim, 8 / 24);
  const wave = sim.frame();
  console.log(`after 8 h the edge of the ring is ${(wave.front / 1e3).toFixed(0)} km out, moving at ${wave.frontSpeed.toFixed(1)} m/s against √(gD) ${wave.gravitySpeed.toFixed(1)} m/s`);
  close(wave.gravitySpeed, Math.sqrt(G * depth), 1e-9, '√(gD)');
  between(wave.frontSpeed / wave.gravitySpeed, 0.9, 1.02, 'front speed against √(gD)');
  runCase(sim);
  const eddy = sim.frame(), speeds = Array.from(eddy.layerSpeeds).filter(Number.isFinite);
  console.log(`after a day the center stands ${(100 * eddy.centreHeight).toFixed(1)} cm up, ${(100 * eddy.kept).toFixed(0)} % of the bump's volume within twice its width; f·u against g∂η/∂n ${eddy.geostrophic.ratio.toFixed(3)}; speeds ${(100 * Math.min(...speeds)).toFixed(3)} to ${(100 * Math.max(...speeds)).toFixed(3)} cm/s over ${speeds.length} layers`);
  between(eddy.centreHeight, 0.12, 0.2, 'the eddy at the center');
  between(eddy.kept, 0.25, 0.4, 'the share of the raised water still within twice its width');
  close(eddy.geostrophic.ratio, 1, 0.05, 'geostrophic eddy');
  assert.ok(Math.max(...speeds) < 1.02 * Math.min(...speeds), 'the same speed in every layer');
  close(eddy.level * 4 * Math.PI * mesh.radius ** 2, Math.PI * 1.5e6 ** 2, 0.05 * Math.PI * 1.5e6 ** 2, 'the volume of the bump spread over the sea');
  assert.equal(sim.ocean.diagnostics().oceanLimited, 0);
});

test('the worker streams a layout and frames whose transferred buffers are fresh, and reports an unknown case', async () => {
  const messages = [];
  globalThis.self = { postMessage: (message, transfer = []) => messages.push({ message, transfer }) };
  await import('../explain/ocean.worker.js');
  const send = (data) => self.onmessage({ data });
  send({ type: 'start', run: 7, kind: 'ekman', options: { days: 0.25 } });
  assert.deepEqual(messages.map((m) => m.message.type), ['layout', 'frame']);
  const [layout, first] = messages.map((m) => m.message);
  assert.equal(layout.run, 7); assert.equal(layout.nCells, C); assert.equal(layout.latitude.length, C); assert.equal(layout.layers, HYBRID_LAYERS.K);
  assert.equal(first.time, 0);
  for (let n = 0; n < 60 && !messages[messages.length - 1].message.done; n++) send({ type: 'advance', seconds: 3 * 3600 });
  const last = messages[messages.length - 1].message;
  assert.equal(last.type, 'frame'); assert.ok(last.done); close(last.time, 0.25 * DAY, 1350, 'stops at the end');
  for (const { transfer } of messages) {
    assert.equal(new Set(transfer).size, transfer.length, 'each buffer once');
    for (const buffer of transfer) assert.ok(buffer instanceof ArrayBuffer && !(buffer instanceof SharedArrayBuffer), 'only plain buffers move');
  }
  send({ type: 'start', run: 8, kind: 'gyre' });
  const error = messages[messages.length - 1].message;
  assert.equal(error.type, 'error'); assert.equal(error.run, 8); assert.match(error.message, /no ocean case called gyre/);
  delete globalThis.self;
});
