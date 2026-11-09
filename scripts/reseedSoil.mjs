#!/usr/bin/env node
/*
 * Re-seeds a saved state's land water from a soil-water file (scripts/soilWater.mjs's format, the committed
 * data/soilwater_march.bin or a seed scaled with --reference --out) as an observed fresh start seeds it: on open land
 * (land off the ice sheets) where the file has a value, the root zone holds the fill times its capacity, the surface
 * layer the fill times its own, and the deep store, when the state carries one, the fill times DEEP_CAPACITY. Everything
 * else in the state is kept but the header's soilSeed, which names the file and its sha256. The capacities are those the
 * state's header records among its land options (given, as scripts/stateOptions.mjs writes them, then options), else
 * the defaults, and the mesh is the state's own (stateGrid). The page builds a state's land from its header's
 * options.land and honors a capacity only where that carries one. A spin-up's save does not: headerOptions writes
 * options.land through soilOptions, which drops the capacities, and records them under given. splitState passes an input
 * header's options.land through as it is and filters only the LAND it is given. So a capacity recorded only under given
 * reaches a spin-up given the same LAND but not the page. The command line warns when the state records a capacity its
 * land uses (the root zone or, without vegetation, the bucket; the surface layer; the deep store) at other than the
 * default and its options.land does not carry it.
 *
 *   node scripts/reseedSoil.mjs state.bin seed.bin out.bin
 */
import { readFileSync, writeFileSync, realpathSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { decodeState, encodeState, f64Names } from '../js/stateFile.module.js';
import { createGeography, decodeSoilWater, geographyOptions, readTopography, savedGeography } from '../js/geography.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { soilSeedRecord } from './soilWater.mjs';
import { DEEP_CAPACITY } from '../js/physics/land.module.js';

export const ROOT_ZONE = 300, SURFACE = 15, BUCKET = 150;

export function reseedLand(saved, soilWater, soilSeed = null) {
  const record = savedGeography(saved), topography = readTopography(record.topography), { mesh } = stateGrid(saved, topography);
  const geography = createGeography(mesh, topography, { ...geographyOptions(record).geography, landCover: null, soilWater });
  const option = (key) => saved.given?.land?.[key] ?? saved.options?.land?.[key], vegetated = option('vegetation') !== false;
  const rootZone = vegetated ? option('rootZoneCapacity') ?? ROOT_ZONE : option('bucketCapacity') ?? BUCKET, surface = option('surfaceCapacity') ?? SURFACE, deep = option('deepCapacity') ?? (vegetated ? DEEP_CAPACITY : 0);
  const { land } = saved, before = Float64Array.from(land.soil);
  let cells = 0;
  for (let i = 0; i < mesh.nCells; i++) {
    if (!geography.land[i] || (geography.iceSheet && geography.iceSheet[i]) || !Number.isFinite(geography.wetness[i])) continue;
    const fill = Math.min(1, Math.max(0, geography.wetness[i]));
    land.soil[i] = rootZone * fill;
    if (land.surface) land.surface[i] = surface * fill;
    if (land.deep) land.deep[i] = deep * fill;
    cells++;
  }
  if (soilSeed) saved.soilSeed = soilSeed;
  else delete saved.soilSeed;
  const mean = (soil) => { let a = 0, s = 0; for (let i = 0; i < mesh.nCells; i++) if (geography.land[i] && !(geography.iceSheet && geography.iceSheet[i])) { a += mesh.areaCell[i]; s += mesh.areaCell[i] * soil[i] / rootZone; } return s / a; };
  const defaults = { ...(vegetated ? { rootZoneCapacity: ROOT_ZONE } : { bucketCapacity: BUCKET }), surfaceCapacity: SURFACE, deepCapacity: vegetated ? DEEP_CAPACITY : 0 };
  const nonDefault = Object.entries(defaults).filter(([key, value]) => option(key) !== undefined && option(key) !== value && saved.options?.land?.[key] !== option(key)).map(([key]) => `${key} ${option(key)}`);
  return { cells, before: mean(before), after: mean(land.soil), deep: !!land.deep, nonDefault };
}

if (process.argv[1] && realpathSync(process.argv[1]) === realpathSync(fileURLToPath(import.meta.url))) {
  const [input, seed, output] = process.argv.slice(2);
  if (!output) throw new Error('usage: node scripts/reseedSoil.mjs state.bin seed.bin out.bin');
  const saved = await decodeState(new Uint8Array(readFileSync(input)));
  const { cells, before, after, deep, nonDefault } = reseedLand(saved, decodeSoilWater(readFileSync(seed)), soilSeedRecord(seed));
  if (nonDefault.length) console.warn(`warning: ${input} records ${nonDefault.join(', ')}; its options.land does not carry them, so the page loads its land at the defaults and only a spin-up given the same LAND keeps this fill`);
  writeFileSync(output, encodeState(saved, { f64: f64Names(saved) }));
  console.log(`N=${saved.N}: re-seeded ${cells} open land cells from ${seed}; mean root-zone fill ${before.toFixed(3)} → ${after.toFixed(3)}${deep ? ', the deep store with it' : '; no deep store saved, so a load starts it from the root zone\'s fill'}; wrote ${output}`);
}
