#!/usr/bin/env node
/*
 * Re-seeds a saved state's land water from a soil-water file (scripts/soilWater.mjs's format, the committed
 * data/soilwater_march.bin or a seed scaled with --reference --out) as an observed fresh start seeds it: on open land
 * (land off the ice sheets) where the file has a value, the root zone holds the fill times its capacity, the surface
 * layer the fill times its own, and the deep store, when the state carries one, the fill times DEEP_CAPACITY. Everything
 * else in the state is kept.
 *
 *   node scripts/reseedSoil.mjs state.bin seed.bin out.bin
 */
import { readFileSync, writeFileSync } from 'node:fs';
import { pathToFileURL } from 'node:url';
import { decodeState, encodeState, f64Names } from '../js/stateFile.module.js';
import { createGeography, decodeSoilWater, geographyOptions, readTopography, savedGeography } from '../js/geography.module.js';
import { Grid } from '../js/grid.module.js';
import { gridMesh, savedRelax } from '../js/stateGrid.module.js';
import { DEEP_CAPACITY } from '../js/physics/land.module.js';

export const ROOT_ZONE = 300, SURFACE = 15;

export function reseedLand(saved, soilWater) {
  const record = savedGeography(saved), mesh = gridMesh(new Grid(saved.N, { relax: savedRelax(saved) ?? undefined }));
  const geography = createGeography(mesh, readTopography(record.topography), { ...geographyOptions(record).geography, landCover: null, soilWater });
  const { land } = saved, before = Float64Array.from(land.soil);
  let cells = 0;
  for (let i = 0; i < mesh.nCells; i++) {
    if (!geography.land[i] || (geography.iceSheet && geography.iceSheet[i]) || !Number.isFinite(geography.wetness[i])) continue;
    const fill = Math.min(1, Math.max(0, geography.wetness[i]));
    land.soil[i] = ROOT_ZONE * fill;
    if (land.surface) land.surface[i] = SURFACE * fill;
    if (land.deep) land.deep[i] = DEEP_CAPACITY * fill;
    cells++;
  }
  const mean = (soil) => { let a = 0, s = 0; for (let i = 0; i < mesh.nCells; i++) if (geography.land[i] && !(geography.iceSheet && geography.iceSheet[i])) { a += mesh.areaCell[i]; s += mesh.areaCell[i] * soil[i] / ROOT_ZONE; } return s / a; };
  return { cells, before: mean(before), after: mean(land.soil), deep: !!land.deep };
}

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
  const [input, seed, output] = process.argv.slice(2);
  if (!output) throw new Error('usage: node scripts/reseedSoil.mjs state.bin seed.bin out.bin');
  const saved = await decodeState(new Uint8Array(readFileSync(input)));
  const { cells, before, after, deep } = reseedLand(saved, decodeSoilWater(readFileSync(seed)));
  writeFileSync(output, encodeState(saved, { f64: f64Names(saved) }));
  console.log(`N=${saved.N}: re-seeded ${cells} open land cells from ${seed}; mean root-zone fill ${before.toFixed(3)} → ${after.toFixed(3)}${deep ? ', the deep store with it' : '; no deep store saved, so a load starts it from the root zone\'s fill'}; wrote ${output}`);
}
