// The top-of-atmosphere cloud radiative effect of saved states by latitude band, from the saves' 30-day means:
//   node scripts/cloudEffectBands.mjs STATE.bin...
// Each save's meanShortwaveCloudEffect and meanLongwaveCloudEffect (absorbed minus clear-sky absorbed sunlight, and the clear-sky minus
// all-sky outgoing longwave) are weighted by mesh.areaCell. The bands go by cell-center latitude: global; 30S-30N is |lat| < 30;
// 30-60 is 30 <= |lat| < 60, both hemispheres pooled; 60-90 is |lat| >= 60. One row a save, then the mean of the saves.
import { readFileSync } from 'node:fs';
import { decodeState } from '../js/stateFile.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { savedGeography, checkGeography, readTopography } from '../js/geography.module.js';

const files = process.argv.slice(2), DEG = 180 / Math.PI;
const BANDS = { global: () => true, '30S-30N': (la) => Math.abs(la) < 30, '30-60': (la) => Math.abs(la) >= 30 && Math.abs(la) < 60, '60-90': (la) => Math.abs(la) >= 60 };
const rows = [];
for (const file of files) {
  const saved = await decodeState(new Uint8Array(readFileSync(file)));
  const mesh = stateGrid(saved, readTopography(checkGeography({ ...savedGeography(saved) }).topography)).mesh;
  const row = {};
  for (const [name, test] of Object.entries(BANDS)) {
    let weight = 0, shortwave = 0, longwave = 0;
    for (let i = 0; i < mesh.nCells; i++) {
      if (!test(mesh.latCell[i] * DEG)) continue;
      const a = mesh.areaCell[i];
      weight += a; shortwave += a * saved.meanShortwaveCloudEffect[i]; longwave += a * saved.meanLongwaveCloudEffect[i];
    }
    row[name] = [shortwave / weight, longwave / weight];
  }
  rows.push(row);
  console.log(`${file.split('/').slice(-2).join('/')}: ${Object.entries(row).map(([name, [sw, lw]]) => `${name} SW ${sw.toFixed(1)} LW ${lw.toFixed(1)}`).join(' | ')}`);
}
if (rows.length > 1) console.log(`mean of ${rows.length} saves: ${Object.keys(BANDS).map((name) => `${name} SW ${(rows.reduce((s, r) => s + r[name][0], 0) / rows.length).toFixed(1)} LW ${(rows.reduce((s, r) => s + r[name][1], 0) / rows.length).toFixed(1)}`).join(' | ')}`);
