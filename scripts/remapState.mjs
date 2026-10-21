// A saved state carried onto another sigma grid:
//   LEVELS=bl42 node scripts/remapState.mjs <in.bin> <out.bin>
// The atmosphere (theta, u, q, qc) goes through remapLevels
// (js/physics/regrid.module.js); everything else but the SST ring, the
// clock included, is kept as saved, so a spin-up continues from the result
// on the new grid.
import { readFileSync, writeFileSync } from 'node:fs';
import { decodeState, encodeState, savedLevels, f64Names } from '../js/stateFile.module.js';
import { withoutSstRing } from '../js/sstMeans.module.js';
import { sigmaInterfaces, sigmaGridName } from '../js/dynamics/sigmaCore.module.js';
import { remapLevels } from '../js/physics/regrid.module.js';
import { stateGrid } from '../js/stateGrid.module.js';

const [input, output] = process.argv.slice(2);
if (!input || !output || !process.env.LEVELS) { console.error('usage: LEVELS=<grid> node scripts/remapState.mjs <in.bin> <out.bin>'); process.exit(1); }
const saved = await decodeState(new Uint8Array(readFileSync(input)));
const from = savedLevels(saved), to = sigmaInterfaces(process.env.LEVELS);
const { mesh, relax } = stateGrid(saved);
const remapped = remapLevels(from, to, saved, mesh);
const state = { ...withoutSstRing(saved), grid: { relax }, K: to.length - 1, levels: to, theta: remapped.theta, u: remapped.u, ...(saved.q ? { q: remapped.q } : {}), ...(saved.qc ? { qc: remapped.qc } : {}) };
writeFileSync(output, encodeState(state, { f64: f64Names(saved) }));
console.log(`${input} (${sigmaGridName(from) ?? 'a saved grid'}, ${from.length - 1} layers, day ${saved.day}) → ${output} (${process.env.LEVELS}, ${to.length - 1} layers)`);
