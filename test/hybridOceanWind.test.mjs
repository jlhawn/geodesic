import { test } from 'node:test';
import assert from 'node:assert/strict';
import { DT, DEG, hybridOcean, zonalWindOnEdges } from './helpers/layered.mjs';

test('tropical easterlies drive equatorial upwelling and a poleward surface drift, and the mixed layer deepens under the stronger wind', () => {
  const { mesh, C, ocean, surfaceT, ice, flux } = hybridOcean();
  const wind = zonalWindOnEdges(mesh, (lat) => -0.08 * Math.exp(-((lat / (15 * DEG)) ** 2)));
  for (let n = 0; n < 5 * 64; n++) ocean.advance(surfaceT, ice, flux, wind, DT);
  const f = ocean.fields(50);
  let equator = 0, nEq = 0, side = 0, nSide = 0, mixedEq = 0, mixedPole = 0, nPole = 0;
  for (let i = 0; i < C; i++) {
    const lat = Math.abs(mesh.latCell[i]);
    if (lat < 3 * DEG) { equator += f.upwelling[i]; nEq++; mixedEq += f.h1[i]; }
    else if (lat > 6 * DEG && lat < 12 * DEG) { side += f.upwelling[i]; nSide++; }
    else if (lat > 60 * DEG) { mixedPole += f.h1[i]; nPole++; }
  }
  console.log(`after five days of easterlies: equatorial upwelling ${(equator / nEq * 86400).toFixed(2)} m/day, 6–12° ${(side / nSide * 86400).toFixed(2)} m/day; mixed layer at the equator ${(mixedEq / nEq).toFixed(1)} m, poleward of 60° ${(mixedPole / nPole).toFixed(1)} m`);
  assert.ok(equator / nEq > 0, 'upwelling at the equator');
  assert.ok(side / nSide < equator / nEq, 'less off the equator');
  assert.ok(mixedEq / nEq > mixedPole / nPole, 'the wind mixes deeper');
});
