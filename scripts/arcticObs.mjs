// Arctic (60-90N) monthly climatologies for the surface sunlight, nothing kept but OUTDIR's binaries and JSON:
//   ERA5 (NCAR RDA ds633.1 monthly means, OPeNDAP, 1-degree subsample, 2001-2020) and CERES-EBAF Ed4.2.1 surface (obs4MIPs, ESGF ORNL, 2001-2024).
//   node scripts/arcticObs.mjs OUTDIR
// The files are era5_NAME.f32 (12 months x 31 rows from 90N to 60N x 360) and ceres_NAME.f32 (12 x 30 rows from 60.5N x 360), float32.
import { writeFileSync, existsSync, readFileSync, mkdirSync } from 'node:fs';
const OUT = process.argv[2];
mkdirSync(OUT, { recursive: true });
const RDA = 'https://thredds.rda.ucar.edu/thredds/dodsC/files/g/d633001';
const ERA = [
  ['msdwswrf', 'e5.moda.fc.sfc.meanflux', '235_035'], ['msdwswrfcs', 'e5.moda.fc.sfc.meanflux', '235_068'], ['mtdwswrf', 'e5.moda.fc.sfc.meanflux', '235_053'],
  ['msnswrf', 'e5.moda.fc.sfc.meanflux', '235_037'], ['msnswrfcs', 'e5.moda.fc.sfc.meanflux', '235_051'], ['mtnswrf', 'e5.moda.fc.sfc.meanflux', '235_039'], ['mtnswrfcs', 'e5.moda.fc.sfc.meanflux', '235_049'],
  ['tcc', 'e5.moda.an.sfc', '128_164'], ['lcc', 'e5.moda.an.sfc', '128_186'], ['tclw', 'e5.moda.an.sfc', '128_078'], ['tciw', 'e5.moda.an.sfc', '128_079'],
  ['tcwv', 'e5.moda.an.sfc', '128_137'], ['tco3', 'e5.moda.an.sfc', '128_206'], ['ci', 'e5.moda.an.sfc', '128_031'], ['fal', 'e5.moda.an.sfc', '128_243'],
];
const YEARS = Array.from({ length: 20 }, (_, k) => 2001 + k), NX = 360, ROWS = 31;
let bytes = 0;
const get = async (url) => { for (let a = 1; ; a++) { try { const r = await fetch(url); if (!r.ok) throw new Error(`${r.status}`); return r; } catch (e) { if (a === 4) throw new Error(`${url}: ${e.message}`); await new Promise((s) => setTimeout(s, 5000 * a)); } } };
const names = {};
async function era(name, path, code, year) {
  const base = `${RDA}/${path}/${year}/${path}.${code}_${name}.ll025sc.${year}010100_${year}120100.nc`;
  if (!names[name]) { const dds = await (await get(`${base}.dds`)).text(); const m = [...dds.matchAll(/Float32 (\w+)\[time/g)].map((x) => x[1]); if (m.length !== 1) throw new Error(`${name}: ${dds.slice(0, 400)}`); names[name] = m[0]; }
  const v = names[name];
  const buf = Buffer.from(await (await get(`${base}.dods?${v}.${v}%5B0:1:11%5D%5B0:4:120%5D%5B0:4:1439%5D`)).arrayBuffer());
  const at = buf.indexOf('\nData:\n') + 7, n = buf.readUInt32BE(at); if (at < 7 || n !== 12 * ROWS * NX) throw new Error(`${name} ${year}: ${n}`);
  bytes += buf.length; const out = new Float64Array(n);
  for (let k = 0; k < n; k++) { const x = buf.readFloatBE(at + 8 + 4 * k); out[k] = Number.isFinite(x) && Math.abs(x) < 1e19 ? x : NaN; }
  return out;
}
const meta = existsSync(`${OUT}/obs.json`) ? JSON.parse(readFileSync(`${OUT}/obs.json`, 'utf8')) : { era5: { source: 'ERA5 monthly means of daily means, NCAR RDA ds633.1 (e5.moda.an.sfc, e5.moda.fc.sfc.meanflux), OPeNDAP 1-degree subsample', years: [2001, 2020], rows: 'lat 90..60 step -1', cols: 'lon 0..359E', fields: {} }, ceres: { source: 'CERES-EBAF Ed4.2.1 surface, obs4MIPs NASA-GISS gn 200003-202509, ESGF ORNL OPeNDAP', years: [2001, 2024], rows: 'lat 60.5..89.5', cols: 'lon 0.5..359.5E', fields: {} } };
for (const [name, path, code] of ERA) {
  if (existsSync(`${OUT}/era5_${name}.f32`)) continue;
  const sum = new Float64Array(12 * ROWS * NX), cnt = new Float64Array(12 * ROWS * NX);
  for (const y of YEARS) { const x = await era(name, path, code, y); for (let k = 0; k < x.length; k++) if (Number.isFinite(x[k])) { sum[k] += x[k]; cnt[k]++; } }
  writeFileSync(`${OUT}/era5_${name}.f32`, Buffer.from(Float32Array.from(sum, (s, k) => (cnt[k] ? s / cnt[k] : NaN)).buffer));
  meta.era5.fields[name] = `${path}.${code}_${name} (${names[name]})`;
  writeFileSync(`${OUT}/obs.json`, JSON.stringify(meta, null, 1));
  console.log(name, `${(bytes / 1e6).toFixed(1)} MB so far`);
}
const U = 'https://esgf-node.ornl.gov/thredds/dodsC/user_pub_work/obs4MIPs/NASA-LaRC/CERES-EBAF-4-2-1/mon';
const file = (v) => `${U}/${v}/gn/v20260220/${v}_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_200003-202509.nc`;
const axis = (text, name) => { const l = text.split('\n'), at = l.findIndex((s) => s.startsWith(`${name}[`)); return l[at + 1].split(',').map(Number); };
const dateOf = (days) => new Date(Date.UTC(2000, 2, 1) + days * 86400e3);
for (const v of ['rsds', 'rsdscs', 'rsdt', 'rsus', 'rsuscs', 'rsut', 'rsutcs']) {
  if (existsSync(`${OUT}/ceres_${v}.f32`)) continue;
  let ax;
  try { ax = await (await get(`${file(v)}.ascii?time,lat,lon`)).text(); } catch (e) { console.log(v, 'unavailable', e.message.slice(0, 200)); continue; }
  const time = axis(ax, 'time'), lat = axis(ax, 'lat'), lon = axis(ax, 'lon');
  if (lat.length !== 180 || lat[150] !== 60.5 || lon[0] !== 0.5) throw new Error(`${v}: grid`);
  const NY = 24, R = 30, out = new Float32Array(12 * R * NX);
  for (let m = 1; m <= 12; m++) {
    const A = 9 + m, Z = A + 12 * (NY - 1);
    for (let t = A; t <= Z; t += 12) { const d = dateOf(time[t]); if (d.getUTCMonth() !== m - 1 || d.getUTCFullYear() !== 2001 + (t - A) / 12) throw new Error(`${v}: time ${t}`); }
    const buf = Buffer.from(await (await get(`${file(v)}.dods?${v}.${v}%5B${A}:12:${Z}%5D%5B150:1:179%5D%5B0:1:359%5D`)).arrayBuffer());
    const at = buf.indexOf('\nData:\n') + 7, n = buf.readUInt32BE(at); if (n !== NY * R * NX) throw new Error(`${v} ${m}: ${n}`);
    bytes += buf.length;
    for (let p = 0; p < R * NX; p++) { let s = 0, c = 0; for (let y = 0; y < NY; y++) { const x = buf.readFloatBE(at + 8 + 4 * (y * R * NX + p)); if (Number.isFinite(x) && Math.abs(x) < 1e19) { s += x; c++; } } out[(m - 1) * R * NX + p] = c ? s / c : NaN; }
  }
  writeFileSync(`${OUT}/ceres_${v}.f32`, Buffer.from(out.buffer));
  meta.ceres.fields[v] = file(v).split('/').pop();
  writeFileSync(`${OUT}/obs.json`, JSON.stringify(meta, null, 1));
  console.log(v, `${(bytes / 1e6).toFixed(1)} MB so far`);
}
console.log('done', `${(bytes / 1e6).toFixed(1)} MB`);
