// Converts a saved state's slab ocean into the hybrid ocean's layers, as
// the page does when it loads a slab state with vertical: 'hybrid':
//   node scripts/convertOcean.mjs <in.bin|in.bin.gz|in.parts.json> <out.bin[.gz]>
// Each column's temperature and salinity are remapped onto the depth table
// at its floor, each edge's velocity onto the mean of its cells' new
// thicknesses, then one regrid to the generator's targets (GENERATOR=hybrid,
// the default, isopycnal or zstar), on the CPU ocean in double precision
// (js/ocean/layered.module.js convertSlab). Everything else in the state
// but its SST ring, the clock included, is kept as saved, and the ocean's
// arrays are written in the precision they arrived in; the output records
// vertical and generator in its header, so the page runs it on the hybrid
// ocean without being asked. A state already on the hybrid layers is
// refused. Split the result for the page with scripts/splitState.mjs.
import { readFileSync, writeFileSync } from 'node:fs';
import { gzipSync } from 'node:zlib';
import { decodeState, encodeState, f64Names } from '../js/stateFile.module.js';
import { withoutSstRing } from '../js/sstMeans.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { createGeography, topographyFromInt16 } from '../js/geography.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { GENERATORS } from '../js/ocean/column.module.js';

const [input, output] = process.argv.slice(2), generator = process.env.GENERATOR ?? 'hybrid';
if (!input || !output) { console.error('usage: [GENERATOR=hybrid|isopycnal|zstar] node scripts/convertOcean.mjs <in.bin|in.parts.json> <out.bin[.gz]>'); process.exit(1); }
if (!GENERATORS.includes(generator)) { console.error(`the generator is one of ${GENERATORS.join(', ')}, not ${generator}`); process.exit(1); }

function readState(file) {
  if (!file.endsWith('.parts.json')) return readFileSync(file);
  const { parts } = JSON.parse(readFileSync(file, 'utf8')), dir = file.replace(/[^/]*$/, '');
  return Buffer.concat(parts.map((part) => readFileSync(dir + part.file)));
}
const saved = await decodeState(new Uint8Array(readState(input)));
if (!saved.ocean || !saved.ocean.h) { console.error(`${input} holds no ocean`); process.exit(1); }
if (saved.ocean.vertical === 'hybrid') { console.error(`${input} is already on the hybrid layers`); process.exit(1); }
const { mesh, relax } = stateGrid(saved), C = mesh.nCells;
const topographyBytes = readFileSync(new URL('../data/topography_0p25.bin', import.meta.url));
const geography = createGeography(mesh, topographyFromInt16(topographyBytes.buffer.slice(topographyBytes.byteOffset, topographyBytes.byteOffset + topographyBytes.byteLength)), {});
const ocean = createOcean(mesh, { geography, vertical: 'hybrid', generator });
const surfaceT = Float64Array.from(saved.surfaceT), ice = Float64Array.from(saved.ice);
if (saved.concentration) ocean.useSeaIce(Float64Array.from(saved.concentration));
const slabLayers = saved.ocean.h.length / C;
const heatOf = (h, T, L) => { let sum = 0; for (let i = 0; i < C; i++) if (ocean.cellOcean[i]) for (let k = 0; k < L; k++) sum += mesh.areaCell[i] * h[k * C + i] * T[k * C + i]; return sum; };
const before = heatOf(saved.ocean.h, saved.ocean.T, slabLayers);
const started = Date.now();
ocean.load(saved.ocean, surfaceT, ice);
const converted = ocean.serialize();
const after = heatOf(converted.h, converted.T, ocean.layers);
const bytes = encodeState({ ...withoutSstRing(saved), grid: { relax }, ocean: converted }, { f64: f64Names(saved) });
writeFileSync(output, output.endsWith('.gz') ? gzipSync(bytes, { level: 9 }) : bytes);
console.log(`${input} (N=${saved.N}, day ${saved.day}, ${slabLayers} slab layers) → ${output} (${ocean.layers} hybrid layers, ${generator} generator, ${((Date.now() - started) / 1000).toFixed(0)} s); the ocean's heat changed by ${((after - before) / Math.abs(before)).toExponential(2)} of itself in the conversion and regrid`);
