// The upward mass flux through each pressure over |φ| < EDGE from an arm's
// TEM means (fixedSst.mjs TEM_HOURS, OUT/TAG_tem.json), against the part
// the gravity-wave drag drives by downward control (Haynes et al. 1991):
// the poleward mass flux above p across ±EDGE, (2π a cos φ / g) ∫ -X/f̂ dp,
// with f̂ = f - (u cos φ)_φ / (a cos φ) from the TEM's zonal wind (and with
// f alone), the drag split into the lid layers' deposit and the breaking
// below them. 10^8 kg/s.
//   node scripts/downwardControl.mjs ARM_tem.json [EDGE]
import { readFileSync } from 'node:fs';
const t = JSON.parse(readFileSync(process.argv[2], 'utf8')), EDGE = Number(process.argv[3] ?? 15);
const g = 9.806, a = 6371220, OMEGA = 7.292e-5, P0 = 101325, B = t.lat.length, J = t.pressures.length, K = t.layers.length;
const at = (field, j, b) => field[j * B + b];
const bandAt = (lat) => { const x = (lat + 90) / 2.5 - 0.5, b = Math.floor(x); return [b, x - b]; };
const interpLat = (get, lat) => { const [b, w] = bandAt(lat); return (1 - w) * get(b) + w * get(b + 1); };
const uAt = (p, b) => {
  let j = 0;
  while (j < J - 2 && t.pressures[j + 1] < p) j++;
  const w = Math.min(1, Math.max(0, Math.log(p / t.pressures[j]) / Math.log(t.pressures[j + 1] / t.pressures[j])));
  return (1 - w) * at(t.u, j, b) + w * at(t.u, j + 1, b);
};
const fHat = (p, lat) => {
  const [b] = bandAt(lat), d = 2.5 * Math.PI / 180, phi = lat * Math.PI / 180;
  const uc = (bb) => uAt(p, bb) * Math.cos((-90 + (bb + 0.5) * 2.5) * Math.PI / 180);
  return 2 * OMEGA * Math.sin(phi) - (uc(b + 1) - uc(b)) / d / (a * Math.cos(phi));
};
function polewardAbove(p, lat, from, to, hat) {
  let sum = 0;
  for (let k = from; k < to; k++) {
    const top = t.sigmaUpper[k] * P0, bottom = Math.min(p, t.sigmaLower[k] * P0);
    if (bottom <= top) break;
    const X = interpLat((b) => t.waveEast[k * B + b], lat), f = hat ? fHat(t.layers[k], lat) : 2 * OMEGA * Math.sin(lat * Math.PI / 180);
    sum += -X / f * (bottom - top);
  }
  return 2 * Math.PI * a * Math.cos(lat * Math.PI / 180) / g * sum;
}
let belt = 0;
for (let b = 0; b < B; b++) if (Math.abs(t.lat[b]) < EDGE) belt += t.area[b];
console.log(`${process.argv[2].replace(/.*\//, '')}: ${t.count ?? '?'} samples, lid ${t.lid} layers (above ${(t.sigmaLower[t.lid - 1] * P0 / 100).toPrecision(3)} hPa), |φ| < ${EDGE}°; 1e8 kg/s, with f̂ (with f)`);
console.log('    hPa   TEM ascent   waves: lid    below lid    all       w* mm/s  (waves alone, with f̂)');
for (let j = 0; j < J; j++) {
  const p = t.pressures[j];
  if (p > 2000) break;
  let M = 0, rho = 0;
  for (let b = 0; b < B; b++) if (Math.abs(t.lat[b]) < EDGE) { M += -at(t.omegaStar, j, b) * t.area[b] / g; }
  const parts = [true, false].map((hat) => [[0, t.lid], [t.lid, K]].map(([from, to]) => polewardAbove(p, EDGE, from, to, hat) - polewardAbove(p, -EDGE, from, to, hat)));
  const [[lid, below], [lidF, belowF]] = parts, all = lid + below;
  const wOf = (m) => m * g * 7000 / (p * belt) * 1e3;
  console.log(`  ${(p / 100).toPrecision(3).padStart(5)}  ${(M / 1e8).toFixed(2).padStart(8)}   ${(lid / 1e8).toFixed(2).padStart(6)} (${(lidF / 1e8).toFixed(2)})  ${(below / 1e8).toFixed(2).padStart(6)} (${(belowF / 1e8).toFixed(2)})  ${(all / 1e8).toFixed(2).padStart(6)}   ${wOf(M).toFixed(2).padStart(6)}  (${wOf(all).toFixed(2)})`);
}
console.log('zonal-mean gravity-wave force by layer, m/s/day, at 0, ±7.5, ±15, ±30 and ±45 degrees, and the absolute flux leaving the layer upward at the equator, mPa');
for (let k = 0; k < Math.min(K, 9); k++) {
  const row = [0, 7.5, -7.5, 15, -15, 30, -30, 45, -45].map((lat) => (86400 * interpLat((b) => t.waveEast[k * B + b], lat)).toFixed(2).padStart(7));
  console.log(`  ${(t.layers[k] / 100).toPrecision(3).padStart(6)} ${row.join('')}   ${(1e3 * interpLat((b) => t.waveFlux[k * B + b], 0)).toFixed(2)}`);
}
