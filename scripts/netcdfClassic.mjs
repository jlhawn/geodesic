import { openSync, readSync, closeSync, fstatSync } from 'node:fs';

/*
 * The variables of a netCDF-3 classic file (CDF-1 or CDF-2) of float or double
 * data, which is all that scripts/aerosolFit.mjs reads of MACv2's files:
 * netcdfClassic(path) gives { dims, vars, numrecs, read(name, record = 0) },
 * read returning one record of a record variable, or a fixed variable,
 * as a Float64Array.
 */
const TYPES = { 1: 1, 2: 1, 3: 2, 4: 4, 5: 4, 6: 8 };

export function netcdfClassic(path) {
  const fd = openSync(path, 'r');
  const head = Buffer.alloc(Math.min(1 << 20, fstatSync(fd).size));
  readSync(fd, head, 0, head.length, 0);
  if (head.toString('latin1', 0, 3) !== 'CDF' || ![1, 2].includes(head[3])) throw new Error(`${path}: not a netCDF-3 classic file`);
  const version = head[3];
  let at = 4;
  const u32 = () => { const v = head.readUInt32BE(at); at += 4; return v; };
  const u64 = () => { const v = Number(head.readBigUInt64BE(at)); at += 8; return v; };
  const text = () => { const n = u32(), s = head.toString('latin1', at, at + n); at += (n + 3) & ~3; return s; };
  const attributes = () => { u32(); for (let n = u32(); n > 0; n--) { text(); const size = TYPES[u32()], count = u32(); at += (count * size + 3) & ~3; } };
  const numrecs = u32();
  u32();
  const dims = Array.from({ length: u32() }, () => [text(), u32()]);
  attributes();
  u32();
  const vars = {};
  for (let n = u32(); n > 0; n--) {
    const name = text(), ids = Array.from({ length: u32() }, u32);
    attributes();
    const type = u32(), vsize = u32(), begin = version === 1 ? u32() : u64();
    const record = ids.length > 0 && dims[ids[0]][1] === 0, shape = ids.map((id) => dims[id][1] || numrecs);
    vars[name] = { type, begin, record, shape, vsize, width: shape.slice(record ? 1 : 0).reduce((a, b) => a * b, TYPES[type]) };
  }
  const records = Object.values(vars).filter((v) => v.record);
  const recordSize = records.length === 1 ? records[0].width : records.reduce((s, v) => s + v.vsize, 0);
  function read(name, record = 0) {
    const v = vars[name], size = TYPES[v.type], count = v.width / size, buffer = Buffer.alloc(v.width);
    if (readSync(fd, buffer, 0, v.width, v.begin + (v.record ? record * recordSize : 0)) !== v.width) throw new Error(`${path}: ${name} is cut short`);
    return Float64Array.from({ length: count }, (_, i) => (v.type === 5 ? buffer.readFloatBE(4 * i) : v.type === 6 ? buffer.readDoubleBE(8 * i) : v.type === 4 ? buffer.readInt32BE(4 * i) : buffer.readInt16BE(2 * i)));
  }
  return { dims: Object.fromEntries(dims), vars, numrecs, read, close: () => closeSync(fd) };
}
