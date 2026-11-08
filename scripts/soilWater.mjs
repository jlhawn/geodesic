#!/usr/bin/env node
/*
 * Writes data/soilwater_march.bin (encodeSoilWater in js/geography.module.js), the fill of the soil's plant-available water at
 * March 20, the model's day 0, on a 0.5° grid: 720 columns from −180° and 360 rows from −90°, cell centers at −179.75 + 0.5 i
 * and −89.75 + 0.5 j, each a uint8 of the fill times 250 and 255 where the source has no value. The source is CPC Soil
 * Moisture, NOAA PSL; Fan and van den Dool 2004, J. Geophys. Res. 109, D10102: the monthly means of 1991–2020 of the water
 * its 760 mm column holds (soilw, mm) on a 0.5° grid, rows from the north and columns from 0.25°E. Fetch it as the binary DAP
 * response of PSL's OPeNDAP server, over HTTP/1.1 (the .ascii form truncates, and the netCDF file is HDF5):
 *
 *   curl --http1.1 -o soilw.dods 'https://psl.noaa.gov/thredds/dodsC/Datasets/cpcsoil/soilw.mon.ltm.1991-2020.v2.nc.dods?soilw'
 *
 * The fill is soilw / 760 interpolated linearly between the mid-month means that bracket day 0 (March and April), clamped to
 * 0–1; a point holds a value where both months hold one of at least 0.
 *
 * CPC's 760 mm column is not the model's root zone: at its own March 20 the model's fill (root-zone water over
 * ROOT_ZONE) matches soilw / 760 in the tropics but runs 2 times it at 40–60N and 3 to 4 times it poleward of 60N. With
 * --reference, the model states listed (March 20 saves of one grid) set the scale: the fill is soilw / 760 times a factor
 * linear in latitude between 10° band centers, fitted so that each band's mean fill over the states' open land matches
 * the states' mean, clamped to 0–1, the factors held to SCALE_RANGE; a band with less than MINIMUM_BAND of the open land
 * takes its factor from the fitted bands either side. --out writes elsewhere than data/soilwater_march.bin.
 *
 *   node scripts/soilWater.mjs soilw.dods
 *   node scripts/soilWater.mjs soilw.dods --reference=a.bin,b.bin --out=soilwater_march_model.bin
 */
import { readFileSync, writeFileSync } from 'node:fs';
import { pathToFileURL } from 'node:url';
import { encodeSoilWater, decodeSoilWater, soilWaterUrl, SOIL_SCALE, SOIL_NONE, createGeography, readTopography, geographyOptions, savedGeography } from '../js/geography.module.js';
import { Grid } from '../js/grid.module.js';
import { gridMesh, savedRelax } from '../js/stateGrid.module.js';
import { decodeState } from '../js/stateFile.module.js';

const COLS = 720, ROWS = 360, MONTHS = 12, CAPACITY = 760, ROOT_ZONE = 300, BAND = 10, SCALE_RANGE = [0.5, 5], MINIMUM_BAND = 0.005;
const MID = [15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349];

function monthWeights(day) {
  let d = (day + 79) % 365;
  if (d < MID[0]) d += 365;
  for (let m = 0; m < 12; m++) {
    const a = MID[m], b = m === 11 ? MID[0] + 365 : MID[m + 1];
    if (d >= a && d < b) return [m, (m + 1) % 12, (d - a) / (b - a)];
  }
  throw new Error(`no month for day ${day}`);
}

function readDods(path) {
  const buf = readFileSync(path), marker = buf.indexOf('Data:\n');
  if (marker < 0) throw new Error(`${path}: no Data: marker`);
  let p = marker + 6;
  const take = (size, read) => {
    const n = buf.readUInt32BE(p);
    if (buf.readUInt32BE(p + 4) !== n) throw new Error(`${path}: count mismatch at byte ${p}`);
    p += 8;
    const values = new Float64Array(n);
    for (let k = 0; k < n; k++) values[k] = read.call(buf, p + size * k);
    p += size * n;
    return values;
  };
  const soilw = take(4, buf.readFloatBE), time = take(8, buf.readDoubleBE), lat = take(4, buf.readFloatBE), lon = take(4, buf.readFloatBE);
  if (p !== buf.length) throw new Error(`${path}: ${buf.length - p} bytes after the maps`);
  return { soilw, time, lat, lon };
}

export async function referenceScale(raw, paths) {
  let mesh = null, geography = null, N = null, relax = null, ground = null;
  const target = [];
  for (const path of paths) {
    const saved = await decodeState(new Uint8Array(readFileSync(path)));
    if (!mesh) {
      const record = savedGeography(saved);
      N = saved.N; relax = savedRelax(saved); ground = JSON.stringify(record); mesh = gridMesh(new Grid(N, { relax: savedRelax(saved) ?? undefined }));
      geography = createGeography(mesh, readTopography(record.topography), { ...geographyOptions(record).geography, landCover: null, soilWater: raw });
    } else if (saved.N !== N || savedRelax(saved) !== relax || JSON.stringify(savedGeography(saved)) !== ground) throw new Error(`${path}: N = ${saved.N}, relax ${savedRelax(saved)}, ${JSON.stringify(savedGeography(saved))}; the first reference is N = ${N}, relax ${relax}, ${ground}`);
    saved.land.soil.forEach((v, i) => { target[i] = (target[i] ?? 0) + v / ROOT_ZONE / paths.length; });
  }
  const { land, iceSheet, wetness } = geography, cells = [];
  for (let i = 0; i < mesh.nCells; i++) if (land[i] && !(iceSheet && iceSheet[i]) && Number.isFinite(wetness[i])) cells.push(i);
  const bands = 180 / BAND, knots = new Float64Array(bands).fill(1), latitude = (i) => mesh.latCell[i] * 180 / Math.PI;
  const factor = (la) => {
    const x = Math.min(bands - 1, Math.max(0, (la + 90) / BAND - 0.5)), b = Math.min(bands - 2, Math.floor(x));
    return knots[b] + (x - b) * (knots[b + 1] - knots[b]);
  };
  const band = (i) => Math.min(bands - 1, Math.floor((latitude(i) + 90) / BAND));
  for (let pass = 0; pass < 40; pass++) {
    const area = new Float64Array(bands), model = new Float64Array(bands), seed = new Float64Array(bands);
    for (const i of cells) { const b = band(i), a = mesh.areaCell[i]; area[b] += a; model[b] += a * Math.min(1, Math.max(0, target[i])); seed[b] += a * Math.min(1, factor(latitude(i)) * wetness[i]); }
    const total = area.reduce((x, y) => x + y, 0), fitted = (b) => area[b] >= MINIMUM_BAND * total && seed[b] > 0;
    for (let b = 0; b < bands; b++) if (fitted(b)) knots[b] = Math.min(SCALE_RANGE[1], Math.max(SCALE_RANGE[0], knots[b] * model[b] / seed[b]));
    for (let b = 0; b < bands; b++) {
      if (fitted(b)) continue;
      let below = b - 1, above = b + 1;
      while (below >= 0 && !fitted(below)) below--;
      while (above < bands && !fitted(above)) above++;
      knots[b] = below < 0 ? knots[above] : above >= bands ? knots[below] : knots[below] + (b - below) / (above - below) * (knots[above] - knots[below]);
    }
  }
  let area = 0, model = 0, before = 0, after = 0;
  for (const i of cells) { const a = mesh.areaCell[i]; area += a; model += a * Math.min(1, Math.max(0, target[i])); before += a * wetness[i]; after += a * Math.min(1, factor(latitude(i)) * wetness[i]); }
  console.log(`N=${N}, ${cells.length} open land cells: mean fill ${(before / area).toFixed(3)} at soilw / ${CAPACITY}, ${(after / area).toFixed(3)} scaled, ${(model / area).toFixed(3)} in the ${paths.length} reference states`);
  console.log(`factors by band center: ${Array.from(knots, (k, b) => `${-90 + BAND * (b + 0.5)}° ${k.toFixed(2)}`).join(', ')}`);
  return factor;
}

export function scaleFill(fill, factor) {
  for (let r = 0; r < ROWS; r++) {
    const scale = factor(-89.75 + 0.5 * r);
    for (let c = 0; c < COLS; c++) if (fill[r * COLS + c] !== SOIL_NONE) fill[r * COLS + c] = Math.round(Math.min(1, scale * fill[r * COLS + c] / SOIL_SCALE) * SOIL_SCALE);
  }
  return fill;
}

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
  const args = process.argv.slice(2), given = (name) => { const arg = args.find((a) => a.startsWith(`--${name}=`)); return arg ? arg.slice(name.length + 3) : null; };
  const inputs = args.filter((a) => !a.startsWith('--'));
  if (given('reference') && !given('out')) throw new Error('a seed scaled to reference states goes to its own file (--out); data/soilwater_march.bin stays the observed fill');
  if (inputs.length !== 1 || args.some((a) => a.startsWith('--') && !/^--(reference|out)=/.test(a))) throw new Error('usage: node scripts/soilWater.mjs soilw.dods [--reference=a.bin,b.bin] [--out=file.bin]');
  const references = given('reference') ? given('reference').split(',') : [], outUrl = given('out') ? pathToFileURL(given('out')) : soilWaterUrl;
  const { soilw, time, lat, lon } = readDods(inputs[0]);
  if (time.length !== MONTHS || lat.length !== ROWS || lon.length !== COLS || soilw.length !== MONTHS * ROWS * COLS) throw new Error(`shape: time ${time.length}, lat ${lat.length}, lon ${lon.length}, soilw ${soilw.length}`);
  lat.forEach((v, j) => { if (Math.abs(v - (89.75 - 0.5 * j)) > 1e-4) throw new Error(`lat[${j}] = ${v}`); });
  lon.forEach((v, k) => { if (Math.abs(v - (0.25 + 0.5 * k)) > 1e-4) throw new Error(`lon[${k}] = ${v}`); });

  const [m0, m1, f] = monthWeights(0), valid = (v) => v >= 0 && v < 1e30;
  const fill = new Uint8Array(ROWS * COLS).fill(SOIL_NONE);
  let points = 0, partial = 0, clamped = 0, largest = 0;
  for (let r = 0; r < ROWS; r++) {
    const j = ROWS - 1 - r;
    for (let c = 0; c < COLS; c++) {
      const k = (c + COLS / 2) % COLS, at = (m) => soilw[(m * ROWS + j) * COLS + k];
      let months = 0;
      for (let m = 0; m < MONTHS; m++) if (valid(at(m))) { months++; largest = Math.max(largest, at(m)); }
      if (months > 0 && months < MONTHS) partial++;
      if (!valid(at(m0)) || !valid(at(m1))) continue;
      const value = ((1 - f) * at(m0) + f * at(m1)) / CAPACITY;
      if (value < 0 || value > 1) clamped++;
      fill[r * COLS + c] = Math.round(Math.min(1, Math.max(0, value)) * SOIL_SCALE);
      points++;
    }
  }
  if (references.length) scaleFill(fill, await referenceScale({ cols: COLS, rows: ROWS, fill: Uint8Array.from(fill) }, references));
  const bytes = new Uint8Array(encodeSoilWater({ cols: COLS, rows: ROWS, fill }));
  writeFileSync(outUrl, bytes);
  const back = decodeSoilWater(readFileSync(outUrl));
  const sample = (la, lo) => { const v = back.fill[Math.round((la + 89.75) / 0.5) * COLS + Math.round((lo + 179.75) / 0.5)]; return v === SOIL_NONE ? 'none' : (v / SOIL_SCALE).toFixed(3); };
  console.log(`months ${m0 + 1} and ${m1 + 1} weighted ${(1 - f).toFixed(4)}/${f.toFixed(4)}; ${points} points with a value, ${partial} with some months missing, ${clamped} clamped, largest soilw ${largest.toFixed(1)} mm`);
  console.log(Object.entries({ Sahara: [23, 10], Amazon: [-3, -60], Europe: [50, 10], 'eastern US': [38, -80], India: [22, 78], Pacific: [0, -140], Greenland: [72, -40], Antarctica: [-80, 0] }).map(([name, [la, lo]]) => `${name} ${sample(la, lo)}`).join(', '));
  console.log(`wrote ${outUrl.pathname} (${bytes.length} bytes)`);
}
