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
 *   node scripts/soilWater.mjs soilw.dods
 */
import { readFileSync, writeFileSync } from 'node:fs';
import { encodeSoilWater, decodeSoilWater, soilWaterUrl, SOIL_SCALE, SOIL_NONE } from '../js/geography.module.js';

const COLS = 720, ROWS = 360, MONTHS = 12, CAPACITY = 760;
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

if (process.argv.length !== 3) throw new Error('usage: node scripts/soilWater.mjs soilw.dods');
const { soilw, time, lat, lon } = readDods(process.argv[2]);
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
const bytes = new Uint8Array(encodeSoilWater({ cols: COLS, rows: ROWS, fill }));
writeFileSync(soilWaterUrl, bytes);
const back = decodeSoilWater(readFileSync(soilWaterUrl));
const sample = (la, lo) => { const v = back.fill[Math.round((la + 89.75) / 0.5) * COLS + Math.round((lo + 179.75) / 0.5)]; return v === SOIL_NONE ? 'none' : (v / SOIL_SCALE).toFixed(3); };
console.log(`months ${m0 + 1} and ${m1 + 1} weighted ${(1 - f).toFixed(4)}/${f.toFixed(4)}; ${points} points with a value, ${partial} with some months missing, ${clamped} clamped, largest soilw ${largest.toFixed(1)} mm`);
console.log(Object.entries({ Sahara: [23, 10], Amazon: [-3, -60], Europe: [50, 10], 'eastern US': [38, -80], India: [22, 78], Pacific: [0, -140], Greenland: [72, -40], Antarctica: [-80, 0] }).map(([name, [la, lo]]) => `${name} ${sample(la, lo)}`).join(', '));
console.log(`wrote ${soilWaterUrl.pathname} (${bytes.length} bytes)`);
