import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel, dragConstants } from '../js/model.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { orographicColumn, namedConstants, OROGRAPHY_DEFAULTS, OROGRAPHY_CONSTANTS, LOTT_MILLER } from '../js/physics/orography.module.js';
import { orographyConstants } from '../js/gpu/orography.gpu.js';

const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);

function ridgeColumn() {
  const K = 30, z = Float64Array.from({ length: K }, (_, k) => 150 * (K - k) - 75), theta = Float64Array.from(z, (h) => 290 * Math.exp(1e-4 * h / 9.80616));
  const p = Float64Array.from(z, (h) => 1e5 * Math.exp(-h / 8000)), edge = (h) => 1e5 * Math.exp(-h / 8000);
  return { z, p, rho: Float64Array.from(p, () => 1.2), theta, east: new Float64Array(K).fill(10), north: new Float64Array(K), pTop: Float64Array.from(z, (h) => edge(h + 75)), pBottom: Float64Array.from(z, (h) => edge(Math.max(0, h - 75))) };
}

test('the orography constants \'cy41r2\' take the IFS\'s operational G, 0.3, over the defaults; constants given by name win, and the raster\'s fields take the set over LOTT_MILLER', () => {
  assert.deepEqual(namedConstants({}), {});
  assert.deepEqual(namedConstants({ constants: 'cy41r2' }), { ...OROGRAPHY_DEFAULTS, waveDrag: 0.3 });
  assert.equal(namedConstants({ constants: 'cy41r2', waveDrag: 0.5 }).waveDrag, 0.5);
  assert.throws(() => namedConstants({ constants: 'cy47r3' }), /cy41r2/);
  assert.deepEqual(dragConstants({ raster: true }, { constants: 'cy41r2' }), OROGRAPHY_CONSTANTS.cy41r2);
  assert.deepEqual(dragConstants({ raster: true }, {}), LOTT_MILLER);
  assert.deepEqual(dragConstants({}, { waveDrag: 2 }), { waveDrag: 2 });
  assert.equal(createModel(new Grid(16), { topography, orography: { constants: 'cy41r2' } }).orography.options.waveDrag, 0.3);
  assert.match(orographyConstants({ orography: dragConstants({}, { constants: 'cy41r2' }) }), /const ORO_G: f32 = 0\.3;/);
});

test('with \'cy41r2\' an unblocked ridge launches 0.3/1.23 of the default\'s wave stress and the same blocking', () => {
  const sub = { deviation: 200, anisotropy: 1, orientation: Math.PI / 2, slope: 0.02 };
  const base = orographicColumn(sub, ridgeColumn(), OROGRAPHY_DEFAULTS), ifs = orographicColumn(sub, ridgeColumn(), OROGRAPHY_CONSTANTS.cy41r2);
  assert.ok(base.launch > 0, `the default launches ${base.launch} N/m²`);
  assert.ok(Math.abs(ifs.launch / base.launch - 0.3 / 1.23) < 1e-12, `launch ratio ${ifs.launch / base.launch}`);
  assert.equal(ifs.blocking, base.blocking);
  assert.deepEqual(Array.from(ifs.beta), Array.from(base.beta));
});
