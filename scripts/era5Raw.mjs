// Range reads of ARCO-ERA5's raw pressure-level netCDF3 files
// (gs://gcp-public-data-arco-era5/raw/date-variable-pressure_level/Y/M/D/<var>/<hPa>.nc,
// one variable, one level, 24 hours of 721x1440 int16 each, public).
// fieldRows(date, variable, level, hour, rowFrom, rowTo) -> Float32Array of
// (rowTo - rowFrom) * 1440 values, rows from 90N southward in 0.25 deg.
const BASE = 'https://storage.googleapis.com/gcp-public-data-arco-era5/raw/date-variable-pressure_level';
const NX = 1440, NY = 721;

async function get(url, from, to, tries = 10) {
  for (let t = 0; ; t++) {
    try {
      const r = await fetch(url, { headers: { Range: `bytes=${from}-${to}` } });
      if (r.status !== 206 && r.status !== 200) throw new Error(`${r.status} ${url}`);
      return new DataView(await r.arrayBuffer());
    } catch (e) { if (t + 1 >= tries) throw e; await new Promise((s) => setTimeout(s, 3000 * (t + 1))); }
  }
}

function parseHeader(v) {
  let o = 0;
  const u32 = () => { const x = v.getUint32(o); o += 4; return x; };
  const name = () => { const n = u32(); const s = new TextDecoder().decode(new Uint8Array(v.buffer, v.byteOffset + o, n)); o += Math.ceil(n / 4) * 4; return s; };
  const magic = u32(), version = magic & 0xff;
  if ((magic >>> 8) !== 0x434446) throw new Error('not netCDF');
  u32();
  const dims = [];
  if (u32() === 0x0a) { const n = u32(); for (let i = 0; i < n; i++) dims.push([name(), u32()]); } else u32();
  const size = [0, 1, 1, 2, 4, 4, 8];
  const attrs = () => {
    const out = {};
    if (u32() !== 0x0c) { u32(); return out; }
    const n = u32();
    for (let i = 0; i < n; i++) {
      const key = name(), type = u32(), count = u32(), bytes = count * size[type];
      let value;
      if (type === 2) value = new TextDecoder().decode(new Uint8Array(v.buffer, v.byteOffset + o, count));
      else if (type === 6) value = v.getFloat64(o);
      else if (type === 5) value = v.getFloat32(o);
      else if (type === 3) value = v.getInt16(o);
      else if (type === 4) value = v.getInt32(o);
      o += Math.ceil(bytes / 4) * 4;
      out[key] = value;
    }
    return out;
  };
  const global = attrs();
  const vars = {};
  if (u32() === 0x0b) {
    const n = u32();
    for (let i = 0; i < n; i++) {
      const vname = name(), nd = u32(), ids = [];
      for (let d = 0; d < nd; d++) ids.push(u32());
      const at = attrs(), type = u32(), vsize = u32();
      const begin = version === 2 ? Number(v.getBigUint64((o += 8) - 8)) : u32();
      vars[vname] = { dims: ids.map((d) => dims[d][0]), at, type, vsize, begin };
    }
  }
  return { dims, global, vars };
}

const headers = new Map();
export async function fieldRows(date, variable, level, hour, rowFrom = 0, rowTo = NY) {
  const url = `${BASE}/${date}/${variable}/${level}.nc`;
  let h = headers.get(url);
  if (!h) { h = parseHeader(await get(url, 0, 16383)); headers.set(url, h); }
  const name = Object.keys(h.vars).find((k) => h.vars[k].dims.length === 3);
  const meta = h.vars[name], scale = meta.at.scale_factor, offset = meta.at.add_offset, fill = meta.at._FillValue, missing = meta.at.missing_value;
  const from = meta.begin + ((hour * NY + rowFrom) * NX) * 2, to = meta.begin + ((hour * NY + rowTo) * NX) * 2 - 1;
  const d = await get(url, from, to), out = new Float32Array((rowTo - rowFrom) * NX);
  for (let m = 0; m < out.length; m++) { const raw = d.getInt16(2 * m); out[m] = raw === fill || raw === missing ? NaN : raw * scale + offset; }
  return out;
}

export async function pool(tasks, width = 8) {
  const results = new Array(tasks.length);
  let next = 0;
  await Promise.all(Array.from({ length: width }, async () => { while (next < tasks.length) { const n = next++; results[n] = await tasks[n](); } }));
  return results;
}

export const ERA5 = { NX, NY, latOfRow: (r) => 90 - 0.25 * r };
