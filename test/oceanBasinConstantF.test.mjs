import { test } from 'node:test';
import assert from 'node:assert/strict';
import { basin, basinMesh, caseMesh } from '../explain/oceanCases.module.js';

test('with the same f everywhere the basin has no Sverdrup interior and its western current is matched by an eastern one', () => {
  const planet = caseMesh(16), flat = basinMesh(16, false), f30 = 2 * planet.omega * Math.sin(Math.PI / 6);
  assert.ok(flat.fEdge.every((f) => f === f30) && flat.fVertex.every((f) => f === f30));
  assert.notEqual(planet.fEdge[0], f30);
  const f = basin({ beta: false });
  console.log(`day ${f.day}: western band ${f.west.toFixed(2)} m²/s north, interior ${f.interior.toFixed(2)}, eastern band ${f.east.toFixed(2)}; ${f.northward.toFixed(1)} Sv north and ${f.southward.toFixed(1)} Sv south across 30°N; strongest 4° band at ${f.strongest}°; fastest ${f.fastest.toFixed(3)} m/s`);
  assert.equal(f.sverdrup, null);
  assert.ok(f.east < 0 && Math.abs(f.west + f.east) < 0.25 * f.west, `western band ${f.west} against the eastern band ${f.east}`);
  assert.ok(f.west < 2 * Math.abs(f.east), 'no western intensification');
  assert.ok(Math.abs(f.interior) < 0.1 * f.west, `interior ${f.interior}`);
  assert.equal(f.limited, 0);
});
