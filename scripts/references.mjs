/*
 * Reads the scorecards' references in data/references/ (README.md there gives each one's source and definition). Each
 * loader checks a binary file against the sha256 its .json records and every file against the shape its writer gives
 * it, and throws on a mismatch, so a score never reads a truncated or replaced file.
 *
 * Gridded files come back as { south, west, step, rows, columns, months, values } with values [month][row][column],
 * rows northward from `south` and columns eastward from `west`, all in degrees. sampleGrid interpolates one bilinearly at
 * a point (a global grid holds its edge rows to the poles; a band throws for a point more than half a step outside
 * it, and any grid for a month index outside 0–11); given ERA5's land-sea mask it weighs the four points by lsm for a land cell and by 1 − lsm for a sea cell, so a
 * coastal cell takes the points of its own kind. boxMean averages the samples over a set of cells by area.
 */
import { readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { decodeMonthlySst } from '../js/ocean/climatology.module.js';

export const REFERENCES = new URL('../data/references/', import.meta.url).pathname;
export const MONTH_DAYS = [31, 28 + 8 / 30, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
export const TAO_STATIONS = ['0n165e', '0n170w', '0n140w', '0n110w'];

const fail = (name, message) => { throw new Error(`data/references/${name}: ${message}`); };
const json = (name, dir) => JSON.parse(readFileSync(`${dir}${name}.json`, 'utf8'));
const twelve = (name, xs, what) => { if (!Array.isArray(xs) || xs.length !== 12 || !xs.every(Number.isFinite)) fail(name, `${what} is not twelve numbers`); return xs; };

export function checkSha256(name, meta, bytes) {
  const got = createHash('sha256').update(bytes).digest('hex');
  if (got !== meta.sha256) fail(`${name}.bin`, `sha256 ${got}, the .json records ${meta.sha256}`);
  return bytes;
}

function binary(name, dir) {
  const meta = json(name, dir), bytes = checkSha256(name, meta, readFileSync(`${dir}${name}.bin`));
  return { meta, bytes, floats: new Float32Array(bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.length)) };
}

const grid = (values, { south, west, step, rows, columns, months }) => ({ south, west, step, rows, columns, months, values });

export function gpcp(dir = REFERENCES) {
  const name = 'gpcp-precip-1991-2020', { meta, floats } = binary(name, dir);
  if (floats.length !== 12 * 72 * 144) fail(name, `${floats.length} values, not 12 × 72 × 144`);
  return { meta, grid: { ...grid(floats, { south: -88.75, west: 1.25, step: 2.5, rows: 72, columns: 144, months: 12 }), global: true } };
}

export function era5Fields(dir = REFERENCES) {
  const name = 'era5-fields-1991-2020', { meta, floats } = binary(name, dir), fields = {};
  const names = meta.fields.map((f) => f.name);
  if (names.join() !== 't2m,msl,metss,lsm') fail(name, `fields ${names}`);
  let at = 0;
  for (const f of meta.fields) {
    if (f.rows !== f.north - f.south + 1 || f.columns !== 360) fail(name, `${f.name} is ${f.rows} × ${f.columns} over ${f.south}–${f.north}`);
    const n = f.months * f.rows * f.columns;
    fields[f.name] = { ...grid(floats.subarray(at, at + n), { south: f.south, west: 0, step: 1, rows: f.rows, columns: f.columns, months: f.months }), units: f.units };
    at += n;
  }
  if (at !== floats.length) fail(name, `${floats.length} values, the fields hold ${at}`);
  return { meta, fields };
}

export function woa23Sst(dir = REFERENCES) {
  const name = 'woa23-sst-1991-2020', { meta, bytes } = binary(name, dir);
  return { meta, sst: decodeMonthlySst(bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.length)) };
}

export function globalTs(dir = REFERENCES) {
  const name = 'global-ts-era5', meta = json(name, dir);
  for (const key of ['monthly', 'air', 'gap']) twelve(name, meta[key], key);
  const annual = (xs) => xs.reduce((a, x, m) => a + x * MONTH_DAYS[m], 0) / MONTH_DAYS.reduce((a, b) => a + b, 0);
  if (Math.abs(annual(meta.monthly) - meta.annual) > 0.001 || Math.abs(annual(meta.air) - meta.airAnnual) > 0.001) fail(name, 'annual is not the month-length-weighted mean of monthly');
  return meta;
}

export function seaIce(dir = REFERENCES) {
  const name = 'nsidc-sea-ice-monthly', meta = json(name, dir);
  for (const key of ['extent', 'area', 'areaPublished']) for (const h of ['N', 'S']) twelve(name, meta[key]?.[h], `${key}.${h}`);
  return meta;
}

export function tao(dir = REFERENCES) {
  const name = 'tao-adcp-equatorial-undercurrent', meta = json(name, dir);
  for (const s of TAO_STATIONS) {
    const st = meta.stations[s];
    if (!st || !Number.isFinite(st.core?.speed) || !Number.isFinite(st.core?.depth) || !(st.profile?.length >= 5)) fail(name, `${s} has no core or profile`);
    if (!st.profile.every((p, k) => Number.isFinite(p.u) && (k === 0 || p.depth > st.profile[k - 1].depth))) fail(name, `${s}: the profile's depths do not increase`);
  }
  return meta;
}

export function ceres(dir = REFERENCES) {
  const name = 'ceres-ebaf-global-monthly', meta = json(name, dir);
  for (const [key, xs] of Object.entries(meta.fields ?? {})) twelve(name, xs, `fields.${key}`);
  return meta;
}

export function sampleGrid(g, month, lat, lon, { lsm = null, land = null } = {}) {
  if (lsm && typeof land !== 'boolean') throw new Error('sampling with lsm needs land: true or false');
  if (!(Number.isInteger(month) && month >= 0 && month < 12 && (g.months === 1 || month < g.months))) throw new Error(`month index ${month} is not 0–11`);
  const north = g.south + (g.rows - 1) * g.step;
  if (!g.global && (lat < g.south - g.step / 2 || lat > north + g.step / 2)) throw new Error(`latitude ${lat} lies outside the grid's ${g.south}–${north}°`);
  const m = g.months === 1 ? 0 : month, y = (lat - g.south) / g.step, x = ((((lon - g.west) % 360) + 360) % 360) / g.step;
  const j0 = Math.min(Math.max(Math.floor(y), 0), g.rows - 2), fy = Math.min(Math.max(y - j0, 0), 1);
  const i0 = Math.floor(x) % g.columns, i1 = (i0 + 1) % g.columns, fx = x - Math.floor(x);
  const points = [[j0, i0, (1 - fy) * (1 - fx)], [j0, i1, (1 - fy) * fx], [j0 + 1, i0, fy * (1 - fx)], [j0 + 1, i1, fy * fx]];
  let sum = 0, weight = 0, plain = 0;
  for (const [j, i, w] of points) {
    const v = g.values[(m * g.rows + j) * g.columns + i];
    plain += w * v;
    if (!lsm) continue;
    const r = Math.round(g.south + j * g.step - lsm.south), c = Math.round((((g.west + i * g.step - lsm.west) % 360) + 360) % 360);
    if (r < 0 || r >= lsm.rows) throw new Error(`lsm covers ${lsm.south}–${lsm.south + lsm.rows - 1}°, not ${g.south + j * g.step}°`);
    const f = lsm.values[r * lsm.columns + c], k = w * (land ? f : 1 - f);
    sum += k * v; weight += k;
  }
  return lsm && weight > 1e-9 ? sum / weight : plain;
}

export function boxMean(mesh, g, month, cells, { lsm = null, land = null } = {}) {
  const deg = 180 / Math.PI;
  let sum = 0, area = 0;
  for (const c of cells) {
    const v = sampleGrid(g, month, mesh.latCell[c] * deg, mesh.lonCell[c] * deg, { lsm, land: land ? !!land[c] : null });
    sum += mesh.areaCell[c] * v; area += mesh.areaCell[c];
  }
  return sum / area;
}

export function gridMean(g, month = null, { south = -90, north = 90 } = {}) {
  let sum = 0, weight = 0;
  for (let m = 0; m < g.months; m++) {
    if (month !== null && m !== month) continue;
    for (let j = 0; j < g.rows; j++) {
      const lat = g.south + j * g.step;
      if (lat < south || lat > north) continue;
      const w = Math.cos(lat * Math.PI / 180);
      for (let i = 0; i < g.columns; i++) { sum += w * g.values[(m * g.rows + j) * g.columns + i]; weight += w; }
    }
  }
  return sum / weight;
}
