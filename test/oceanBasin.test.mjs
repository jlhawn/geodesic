import { test } from 'node:test';
import assert from 'node:assert/strict';
import { basinCase, runCase } from '../explain/oceanCases.module.js';
import { gyreReadings, DEPTH } from '../explain/gyreReadings.module.js';
import { threaded } from './helpers/oceanThreads.mjs';

const between = (v, lo, hi, what) => assert.ok(v >= lo && v <= hi, `${what}: ${v} not within ${lo} to ${hi}`);
const describe = (f, r) => `day ${f.day}: western band ${(100 * r.west).toFixed(3)} cm/s north, interior ${(100 * r.interior).toFixed(4)}, eastern band ${(100 * r.east).toFixed(3)}; current ${(r.current.width / 1e3).toFixed(0)} km wide at ${(100 * r.current.peak).toFixed(3)} cm/s carrying ${r.current.transport.toFixed(1)} Sv; ${f.northward.toFixed(1)} Sv north and ${f.southward.toFixed(1)} Sv south across 30°N; ratio ${r.ratio.toFixed(1)}; strongest 4° band at ${f.strongest}°; fastest ${f.fastest.toFixed(3)} m/s`;

test('the subtropical wind over a walled basin on the sphere drives a western boundary current that carries about what Sverdrup’s balance asks, with an interior 13 to 15 times slower', async () => {
  const sim = await threaded(basinCase());
  for (const day of [20, 100]) {
    const f = runCase(sim, day), r = gyreReadings(f);
    console.log(`${describe(f, r)}; Sverdrup ${(100 * r.sverdrup.speed).toFixed(4)} cm/s, ${r.sverdrup.transport.toFixed(1)} Sv across the basin`);
    assert.ok(f.interior < 0, `the interior returns water south: ${f.interior}`);
    assert.ok(f.west >= 10 * Math.abs(f.interior), `western band ${f.west} against the interior ${f.interior}`);
    assert.ok(f.west >= 5 * Math.abs(f.east), `western band ${f.west} against the eastern band ${f.east}`);
    assert.ok(f.strongest === -38 || f.strongest === -34, `strongest band at ${f.strongest}`);
    assert.ok(Math.abs(f.northward + f.southward) < 0.05 * f.northward, 'what goes north across 30°N comes back south');
    assert.ok(Math.abs(f.northward / -f.sverdrupTotal - 1) < 0.2, `${f.northward} Sv against Sverdrup's ${-f.sverdrupTotal}`);
    between(r.ratio, 12, 16, 'the western band against the interior'); between(r.current.width / 1e3, 800, 1050, 'width of the western current in km');
    assert.equal(f.limited, 0);
    if (day === 20) {
      between(100 * r.west, 1.1, 1.3, 'western band cm/s on day 20'); between(r.current.transport, 35, 41, 'Sv in the western current on day 20'); between(f.northward, 38, 43, 'Sv north across 30°N on day 20');
      between(-100 * r.interior, 0.07, 0.09, 'interior cm/s on day 20'); between(f.interior / f.sverdrup, 0.6, 0.75, 'the interior is still catching up with Sverdrup’s balance');
    } else {
      between(100 * r.west, 1.15, 1.35, 'western band cm/s on day 100'); between(r.current.transport, 37, 42, 'Sv in the western current on day 100'); between(-100 * r.interior, 0.09, 0.105, 'interior cm/s on day 100'); between(f.interior / f.sverdrup, 0.7, 1.0, 'the interior against Sverdrup’s balance');
      between(r.across.north, 38, 43, 'Sv north across 30°N on day 100');
    }
  }
  await sim.close();
  const r = gyreReadings(sim.frame());
  assert.ok(Math.abs(r.sverdrup.speed * DEPTH + 4.63) < 0.01 && Math.abs(r.sverdrup.transport + 35.7) < 0.1, 'Sverdrup’s balance at 30°N for 0.1 N/m²: −4.63 m²/s, −35.7 Sv');
});
