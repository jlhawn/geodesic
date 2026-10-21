import { test } from 'node:test';
import assert from 'node:assert/strict';
import { caseMesh, bumpCase, runCase, G } from '../explain/oceanCases.module.js';

const N = 16;
const mesh = caseMesh(N);
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const between = (v, lo, hi, what) => assert.ok(v >= lo && v <= hi, `${what}: ${v} not within ${lo} to ${hi}`);

test('a 1 m bump of sea surface races away at about √(gD) and leaves a geostrophic eddy that every layer circles alike, keeping the share of its mass that the lens’s bottom pressure keeps', () => {
  const sim = bumpCase({ N, days: 1 }), depth = 4000, layout = sim.layout();
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
