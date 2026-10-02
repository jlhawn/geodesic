// A saved state carried onto another sigma grid:
//   LEVELS=bl36 node scripts/remapState.mjs <in.bin> <out.bin>
// The atmosphere (theta, u, q, qc) goes through remapLevels
// (js/physics/regrid.module.js), each layer the σ-weighted mean of the
// layers it overlaps; everything else, the clock included, is kept as
// saved, so a spin-up continues from the result on the new grid.
import { readFileSync, writeFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { decodeState, encodeState, savedLevels } from '../js/stateFile.module.js';
import { sigmaInterfaces, sigmaGridName } from '../js/dynamics/sigmaCore.module.js';
import { remapLevels } from '../js/physics/regrid.module.js';

const [input, output] = process.argv.slice(2);
if (!input || !output || !process.env.LEVELS) { console.error('usage: LEVELS=<grid> node scripts/remapState.mjs <in.bin> <out.bin>'); process.exit(1); }
const saved = await decodeState(new Uint8Array(readFileSync(input)));
const from = savedLevels(saved), to = sigmaInterfaces(process.env.LEVELS);
const remapped = remapLevels(from, to, saved, buildMesh(new Grid(saved.N)));
const f64 = [];
for (const [key, value] of Object.entries(saved)) {
  if (value instanceof Float64Array) f64.push(key);
  else if ((key === 'ocean' || key === 'land') && value) for (const [inner, values] of Object.entries(value)) if (values instanceof Float64Array) f64.push(`${key}.${inner}`);
}
const state = { ...saved, K: to.length - 1, levels: to, theta: remapped.theta, u: remapped.u, ...(saved.q ? { q: remapped.q } : {}), ...(saved.qc ? { qc: remapped.qc } : {}) };
writeFileSync(output, encodeState(state, { f64 }));
console.log(`${input} (${sigmaGridName(from) ?? 'a saved grid'}, ${from.length - 1} layers, day ${saved.day}) → ${output} (${process.env.LEVELS}, ${to.length - 1} layers)`);
