#!/usr/bin/env node
/*
 * Writes data/lakedepth.bin (encodeLakeDepths in js/physics/lake.module.js) from HydroLAKES v1.0 (Messager et al. 2016,
 * Nat. Commun. 7, 13603; CC BY 4.0): every lake of at least MINIMUM_AREA km² with its pour point (longitude, latitude),
 * area (km²), mean depth (m, total volume over area) and the source of that volume (1 reported, 2 a reservoir's from
 * GRanD, 3 the geostatistical estimate). Fetch the pour points' attributes as
 *
 *   curl -L -O https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_points_v10_shp.zip
 *   unzip -j HydroLAKES_points_v10_shp.zip HydroLAKES_points_v10_shp/HydroLAKES_points_v10.dbf
 *
 *   node scripts/lakeDepths.mjs HydroLAKES_points_v10.dbf
 */
import { openSync, readSync, closeSync, writeFileSync } from 'node:fs';
import { encodeLakeDepths, lakeDepthsUrl } from '../js/physics/lake.module.js';

const MINIMUM_AREA = 100;

function readDbf(path, wanted) {
  const fd = openSync(path, 'r'), head = Buffer.alloc(32);
  readSync(fd, head, 0, 32, 0);
  const count = head.readUInt32LE(4), headerLength = head.readUInt16LE(8), recordLength = head.readUInt16LE(10);
  const fields = Buffer.alloc(headerLength);
  readSync(fd, fields, 0, headerLength, 0);
  const columns = {};
  for (let at = 32, offset = 1; fields[at] !== 0x0d; at += 32) {
    const name = fields.toString('latin1', at, at + 11).replace(/\0.*$/, ''), length = fields[at + 16];
    columns[name] = { offset, length };
    offset += length;
  }
  for (const name of wanted) if (!columns[name]) throw new Error(`${path}: no field ${name}`);
  const rows = [], chunk = 65536, buffer = Buffer.alloc(chunk * recordLength);
  for (let first = 0; first < count; first += chunk) {
    const n = Math.min(chunk, count - first);
    readSync(fd, buffer, 0, n * recordLength, headerLength + first * recordLength);
    for (let r = 0; r < n; r++) {
      const base = r * recordLength;
      if (buffer[base] === 0x2a) continue;
      rows.push(Object.fromEntries(wanted.map((name) => { const { offset, length } = columns[name], text = buffer.toString('latin1', base + offset, base + offset + length).trim(); return [name, name === 'Lake_name' ? text : Number(text)]; })));
    }
  }
  closeSync(fd);
  return rows;
}

if (process.argv.length !== 3) throw new Error('usage: node scripts/lakeDepths.mjs HydroLAKES_points_v10.dbf');
const rows = readDbf(process.argv[2], ['Hylak_id', 'Lake_name', 'Lake_area', 'Depth_avg', 'Vol_src', 'Pour_long', 'Pour_lat']);
const kept = rows.filter((row) => row.Lake_area >= MINIMUM_AREA && row.Depth_avg > 0).sort((a, b) => b.Lake_area - a.Lake_area);
const bytes = new Uint8Array(encodeLakeDepths(kept.map((row) => ({ lon: row.Pour_long, lat: row.Pour_lat, area: row.Lake_area, depth: row.Depth_avg, source: row.Vol_src }))));
writeFileSync(lakeDepthsUrl, bytes);
const bySource = [1, 2, 3].map((s) => kept.filter((row) => row.Vol_src === s).length);
console.log(`${rows.length} lakes, ${kept.length} of at least ${MINIMUM_AREA} km² (volume reported ${bySource[0]}, reservoir ${bySource[1]}, estimated ${bySource[2]})`);
console.log(kept.slice(0, 12).map((row) => `${row.Lake_name || row.Hylak_id} ${row.Lake_area.toFixed(0)} km² ${row.Depth_avg} m`).join(', '));
console.log(`wrote ${lakeDepthsUrl.pathname} (${bytes.length} bytes)`);
