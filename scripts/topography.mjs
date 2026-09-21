#!/usr/bin/env node
// Writes data/topography_0p25.bin: global ice-surface elevation (ETOPO1, i.e.
// the top of the Antarctic and Greenland ice sheets, not bedrock) on a
// regular 0.25 degree grid, int16 little-endian metres, 1440 columns x 720
// rows, row-major. Row 0 is centred at +89.875N, column 0 at -179.875W;
// rows run north to south, columns west to east. Each value is the area
// mean of the source 5 arc-minute data over the 0.25 degree cell (a 3x3
// block average), not a point sample.

import { writeFile } from 'node:fs/promises';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const SOURCE_URL = 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.nc'
  + '?altitude%5B(-90.0):5:(90.0)%5D%5B(-180.0):5:(180.0)%5D';

const OUT_COLS = 1440;
const OUT_ROWS = 720;
const OUT_PATH = join(dirname(fileURLToPath(import.meta.url)), '..', 'data', 'topography_0p25.bin');

const NC_TYPE_SIZE = { 1: 1, 2: 1, 3: 2, 4: 4, 5: 4, 6: 8 };

function readName(view, off) {
  const n = view.getInt32(off);
  off += 4;
  const bytes = new Uint8Array(view.buffer, view.byteOffset + off, n);
  const name = Buffer.from(bytes).toString('ascii');
  off += n + ((4 - (n % 4)) % 4);
  return [name, off];
}

function skipAttrs(view, off) {
  const count = view.getInt32(off + 4);
  off += 8;
  for (let i = 0; i < count; i++) {
    [, off] = readName(view, off);
    const type = view.getInt32(off);
    const n = view.getInt32(off + 4);
    off += 8;
    const size = NC_TYPE_SIZE[type] * n;
    off += size + ((4 - (size % 4)) % 4);
  }
  return off;
}

function parseNetCDF3(buffer) {
  const view = new DataView(buffer.buffer, buffer.byteOffset, buffer.byteLength);
  if (view.getUint8(0) !== 0x43 || view.getUint8(1) !== 0x44 || view.getUint8(2) !== 0x46) {
    throw new Error('not a NetCDF3 classic file');
  }
  let off = 8;
  const dims = [];
  const dimCount = view.getInt32(off + 4);
  off += 8;
  for (let i = 0; i < dimCount; i++) {
    let name;
    [name, off] = readName(view, off);
    const length = view.getInt32(off);
    off += 4;
    dims.push({ name, length });
  }
  off = skipAttrs(view, off);
  const vars = [];
  const varCount = view.getInt32(off + 4);
  off += 8;
  for (let i = 0; i < varCount; i++) {
    let name;
    [name, off] = readName(view, off);
    const ndims = view.getInt32(off);
    off += 4;
    const dimids = [];
    for (let d = 0; d < ndims; d++) {
      dimids.push(view.getInt32(off));
      off += 4;
    }
    off = skipAttrs(view, off);
    const type = view.getInt32(off);
    off += 8;
    const begin = view.getInt32(off);
    off += 4;
    vars.push({ name, dimids, type, begin });
  }
  return { dims, vars, view };
}

function readVar(nc, name) {
  const v = nc.vars.find((x) => x.name === name);
  if (!v) throw new Error(`variable ${name} not found in source file`);
  const shape = v.dimids.map((d) => nc.dims[d].length);
  const count = shape.reduce((a, b) => a * b, 1);
  const { view } = nc;
  let values;
  if (v.type === 3) {
    values = new Float64Array(count);
    for (let i = 0; i < count; i++) values[i] = view.getInt16(v.begin + 2 * i);
  } else if (v.type === 5) {
    values = new Float64Array(count);
    for (let i = 0; i < count; i++) values[i] = view.getFloat32(v.begin + 4 * i);
  } else if (v.type === 6) {
    values = new Float64Array(count);
    for (let i = 0; i < count; i++) values[i] = view.getFloat64(v.begin + 8 * i);
  } else {
    throw new Error(`unsupported NetCDF type ${v.type} for variable ${name}`);
  }
  return { shape, values };
}

async function download(url) {
  process.stderr.write(`fetching ${url}\n`);
  const t0 = Date.now();
  const res = await fetch(url);
  if (!res.ok) throw new Error(`download failed: ${res.status} ${res.statusText}`);
  const buffer = Buffer.from(await res.arrayBuffer());
  process.stderr.write(`downloaded ${buffer.length} bytes in ${((Date.now() - t0) / 1000).toFixed(1)}s\n`);
  return buffer;
}

function areaAverage(altitude, nLatSrc, nLonSrc) {
  const latTrim = nLatSrc - (nLatSrc % 3);
  const lonTrim = nLonSrc - (nLonSrc % 3);
  const outRows = latTrim / 3;
  const outCols = lonTrim / 3;
  if (outRows !== OUT_ROWS || outCols !== OUT_COLS) {
    throw new Error(`unexpected source grid ${nLatSrc}x${nLonSrc}, expected a 3x downsample to ${OUT_ROWS}x${OUT_COLS}`);
  }
  const out = new Uint8Array(OUT_ROWS * OUT_COLS * 2);
  const outView = new DataView(out.buffer);
  for (let sj = 0; sj < outRows; sj++) {
    const r = outRows - 1 - sj; // source rows run south to north; output rows run north to south
    for (let sk = 0; sk < outCols; sk++) {
      let sum = 0;
      for (let di = 0; di < 3; di++) {
        for (let dj = 0; dj < 3; dj++) {
          sum += altitude[(sj * 3 + di) * nLonSrc + (sk * 3 + dj)];
        }
      }
      outView.setInt16((r * outCols + sk) * 2, Math.round(sum / 9), true);
    }
  }
  return out;
}

async function main() {
  const buffer = await download(SOURCE_URL);
  const nc = parseNetCDF3(buffer);
  const { shape: latShape, values: lat } = readVar(nc, 'latitude');
  const { shape: lonShape, values: lon } = readVar(nc, 'longitude');
  const { shape, values: altitude } = readVar(nc, 'altitude');
  const [nLatSrc, nLonSrc] = shape;
  process.stderr.write(`source grid ${nLatSrc} x ${nLonSrc}, lat ${lat[0]}..${lat[latShape[0] - 1]}, lon ${lon[0]}..${lon[lonShape[0] - 1]}\n`);

  const out = areaAverage(altitude, nLatSrc, nLonSrc);
  await writeFile(OUT_PATH, out);
  process.stderr.write(`wrote ${OUT_PATH} (${out.length} bytes)\n`);
}

main().catch((err) => {
  process.stderr.write(`${err.stack || err}\n`);
  process.exitCode = 1;
});
