/*
 * A netCDF-3 classic (CDF-1) file of float variables, for the tests of scripts/packAerosol.py and
 * scripts/netcdfClassic.mjs: writeClassic({ dims: [[name, length]] (length 0 for the record dimension,
 * first), numrecs, vars: [{ name, dims: [dimension names], values }] }) gives its bytes. A record
 * variable's values hold its records one after the other.
 */
export function writeClassic({ dims, numrecs = 0, vars }) {
  const pad = (n) => (4 - (n % 4)) % 4, ids = new Map(dims.map(([name], n) => [name, n]));
  const size = (v) => 4 * v.dims.reduce((a, name) => a * (dims[ids.get(name)][1] || 1), 1);
  const text = (s) => { const b = Buffer.alloc(4 + s.length + pad(s.length)); b.writeUInt32BE(s.length, 0); b.write(s, 4, 'latin1'); return b; };
  const u32 = (...values) => { const b = Buffer.alloc(4 * values.length); values.forEach((v, n) => b.writeUInt32BE(v, 4 * n)); return b; };
  const record = (v) => v.dims.length > 0 && dims[ids.get(v.dims[0])][1] === 0;
  const head = (offsets) => Buffer.concat([
    Buffer.from('CDF\x01', 'latin1'), u32(numrecs),
    u32(0x0A, dims.length), ...dims.map(([name, length]) => Buffer.concat([text(name), u32(length)])),
    u32(0, 0),
    u32(0x0B, vars.length), ...vars.map((v, n) => Buffer.concat([text(v.name), u32(v.dims.length), u32(...v.dims.map((d) => ids.get(d))), u32(0, 0), u32(5, size(v)), u32(offsets[n])])),
  ]);
  const fixed = vars.filter((v) => !record(v)), recorded = vars.filter(record);
  const start = head(vars.map(() => 0)).length, offsets = new Array(vars.length);
  let at = start;
  for (const v of fixed) { offsets[vars.indexOf(v)] = at; at += size(v); }
  const first = at;
  let within = 0;
  for (const v of recorded) { offsets[vars.indexOf(v)] = first + within; within += size(v); }
  const perRecord = within;
  const out = Buffer.alloc(first + perRecord * numrecs);
  head(offsets).copy(out, 0);
  const put = (position, values) => values.forEach((value, n) => out.writeFloatBE(value, position + 4 * n));
  for (const v of fixed) put(offsets[vars.indexOf(v)], v.values);
  for (const v of recorded) {
    const width = v.values.length / numrecs;
    for (let r = 0; r < numrecs; r++) put(offsets[vars.indexOf(v)] + r * perRecord, v.values.subarray(r * width, (r + 1) * width));
  }
  return out;
}
