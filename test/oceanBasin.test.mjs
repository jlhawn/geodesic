import { test } from 'node:test';
import assert from 'node:assert/strict';
import { basin } from '../explain/oceanCases.module.js';

const describe = (f) => `day ${f.day}: western band ${f.west.toFixed(2)} m²/s north, interior ${f.interior.toFixed(2)}, eastern band ${f.east.toFixed(2)}; ${f.northward.toFixed(1)} Sv north and ${f.southward.toFixed(1)} Sv south across 30°N; strongest 4° band at ${f.strongest}°; fastest ${f.fastest.toFixed(3)} m/s`;

test('the subtropical wind over a walled basin on the sphere drives a western boundary current carrying over 10 times the interior\'s return after 100 days', () => {
  const f = basin();
  console.log(`${describe(f)}; Sverdrup ${f.sverdrup.toFixed(2)} m²/s, ${f.sverdrupTotal.toFixed(1)} Sv across the basin`);
  assert.ok(f.interior < 0, `the interior returns water south: ${f.interior}`);
  assert.ok(f.west >= 10 * Math.abs(f.interior), `western band ${f.west} against the interior ${f.interior}`);
  assert.ok(f.west >= 5 * Math.abs(f.east), `western band ${f.west} against the eastern band ${f.east}`);
  assert.ok(f.strongest === -38 || f.strongest === -34, `strongest band at ${f.strongest}`);
  assert.ok(Math.abs(f.northward + f.southward) < 0.05 * f.northward, 'what goes north across 30°N comes back south');
  assert.ok(Math.abs(f.interior / f.sverdrup - 1) < 0.3, `interior ${f.interior} against Sverdrup ${f.sverdrup}`);
  assert.ok(Math.abs(f.northward / -f.sverdrupTotal - 1) < 0.3, `${f.northward} Sv against Sverdrup's ${-f.sverdrupTotal}`);
  assert.equal(f.limited, 0);
});
