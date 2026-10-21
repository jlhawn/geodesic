import { test } from 'node:test';
import assert from 'node:assert/strict';
import { caseMesh, lensCase, runCase, DEG, DAY } from '../explain/oceanCases.module.js';
import { HYBRID_LAYERS } from '../js/ocean/column.module.js';

const N = 16;
const mesh = caseMesh(N);
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const between = (v, lo, hi, what) => assert.ok(v >= lo && v <= hi, `${what}: ${v} not within ${lo} to ${hi}`);

test('a lens of one layer 60 m thick pushes the layers beneath it down and raises the sea surface by Δρ/ρ₀ of that, less the bottom pressure the barotropic waves leave, with the surface current geostrophic and the deep water slow', () => {
  const sim = lensCase({ N, days: 30 }), layout = sim.layout(), { ocean } = sim, Cn = sim.mesh.nCells;
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

  const flat = lensCase({ N, depth: 0 }), settled = Float64Array.from(flat.ocean.h);
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
