import { test } from 'node:test';
import assert from 'node:assert/strict';
import { basin, basinCase, basinMesh, caseMesh, runCase } from '../explain/oceanCases.module.js';
import { gyreReadings } from '../explain/gyreReadings.module.js';

const between = (v, lo, hi, what) => assert.ok(v >= lo && v <= hi, `${what}: ${v} not within ${lo} to ${hi}`);

test('with the same f everywhere the basin has no Sverdrup interior and its western current is matched by an eastern one that spins up', () => {
  const planet = caseMesh(16), flat = basinMesh(16, false), f30 = 2 * planet.omega * Math.sin(Math.PI / 6);
  assert.ok(flat.fEdge.every((f) => f === f30) && flat.fVertex.every((f) => f === f30));
  assert.notEqual(planet.fEdge[0], f30);
  const sim = basinCase({ beta: false }), seen = {};
  for (const day of [20, 100]) {
    const f = runCase(sim, day);
    seen[day] = f;
    console.log(`day ${f.day}: western band ${(100 * f.west / 4000).toFixed(3)} cm/s north, interior ${(100 * f.interior / 4000).toFixed(4)}, eastern band ${(100 * f.east / 4000).toFixed(3)}; ${f.northward.toFixed(1)} Sv north and ${f.southward.toFixed(1)} Sv south across 30°N; strongest 4° band at ${f.strongest}°; fastest ${f.fastest.toFixed(3)} m/s`);
    assert.equal(f.sverdrup, null);
    assert.ok(f.east < 0 && Math.abs(f.west + f.east) < 0.4 * f.west, `western band ${f.west} against the eastern band ${f.east}`);
    assert.ok(f.west < 2 * Math.abs(f.east), 'no western intensification');
    assert.ok(Math.abs(f.interior) < 0.1 * f.west, `interior ${f.interior}`);
    assert.ok(Math.abs(f.northward + f.southward) < 0.05 * f.northward, 'what goes north comes back south');
    assert.equal(f.limited, 0);
  }
  between(seen[20].northward, 70, 80, 'Sv each way across 30°N on day 20'); between(seen[100].northward, 245, 270, 'Sv each way across 30°N on day 100');
  between(100 * seen[100].west / 4000, 1.7, 1.95, 'western band cm/s on day 100'); between(-100 * seen[100].east / 4000, 1.25, 1.45, 'eastern band cm/s on day 100');
  assert.ok(seen[100].northward > 3 * seen[20].northward, 'the gyre keeps spinning up');
});

test('doubling the wind doubles the speeds and the transports on day 20', () => {
  const one = gyreReadings(basin({ days: 20 })), two = gyreReadings(basin({ days: 20, tau0: 0.2 }));
  const ratios = [two.west / one.west, two.interior / one.interior, two.current.transport / one.current.transport, two.across.north / one.across.north];
  console.log(`doubling the wind: western band ${ratios[0].toFixed(3)}, interior ${ratios[1].toFixed(3)}, western current ${ratios[2].toFixed(3)}, north across 30°N ${ratios[3].toFixed(3)} times`);
  for (const r of ratios) between(r, 1.9, 2.1, 'ratio of speeds and transports');
  between(two.current.transport, 70, 80, 'Sv in the western current at 0.2 N/m²');
});
