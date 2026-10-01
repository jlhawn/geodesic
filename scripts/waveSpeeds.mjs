// The gravity-wave speeds of the mean equatorial column (2S-2N, 170W-110W) of
// saved states, from the layered model's linear modes (the mixed layer at its
// own density, the classes at their labels), with the crossing times of the
// Pacific from 150E to 80W by the first baroclinic Kelvin wave, its first
// meridional Rossby wave at a third of its speed, and the second Kelvin wave.
//   node scripts/waveSpeeds.mjs <state.bin>...
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { decodeState } from '../js/stateFile.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
const meshes = {};
function jacobi(A) {
  const n = A.length, a = A.map((r) => r.slice());
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += a[p][q] ** 2;
    if (off < 1e-22) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(a[p][q]) < 1e-30) continue;
      const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]), t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { const akp = a[k][p], akq = a[k][q]; a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq; }
      for (let k = 0; k < n; k++) { const apk = a[p][k], aqk = a[q][k]; a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk; }
    }
  }
  return a.map((r, i) => r[i]).sort((x, y) => y - x);
}
for (const file of process.argv.slice(2)) {
  const s = await decodeState(new Uint8Array(readFileSync(file)));
  const mesh = meshes[s.N] ?? (meshes[s.N] = buildMesh(new Grid(s.N)));
  const C = mesh.nCells, L = s.ocean.h.length / C, deg = 180 / Math.PI, g = 9.81, rho0 = 1025;
  const H = new Float64Array(L), R = new Float64Array(L); let n = 0;
  for (let i = 0; i < C; i++) {
    const lon = ((mesh.lonCell[i] * deg) % 360 + 360) % 360;
    if (Math.abs(mesh.latCell[i] * deg) > 2 || lon < 190 || lon > 250 || !(s.ocean.h[i] > 1)) continue;
    n++;
    for (let k = 0; k < L; k++) H[k] += s.ocean.h[k * C + i];
    R[0] += seawaterDensity(s.ocean.T[i], s.ocean.S[i]);
  }
  for (let k = 0; k < L; k++) H[k] /= n;
  R[0] /= n; for (let k = 1; k < L; k++) R[k] = s.ocean.densities[k - 1];
  const keep = [...H.keys()].filter((k) => H[k] > 1);
  const A = keep.map((k) => keep.map((j) => Math.sqrt(H[k] * H[j]) * g * R[Math.min(j, k)] / rho0));
  const c = jacobi(A).map((x) => Math.sqrt(Math.max(0, x)));
  const L1 = 135 * 111e3;
  console.log(`${file.split('/').pop()}: ${n} cells, depth ${H.reduce((a, b) => a + b, 0).toFixed(0)} m; speeds ${c.slice(0, 4).map((x) => x.toFixed(2)).join(', ')} m/s; first baroclinic Kelvin over 150E-80W (${(L1 / 1e3).toFixed(0)} km) ${(L1 / c[1] / 86400).toFixed(0)} d, first Rossby back ${(3 * L1 / c[1] / 86400).toFixed(0)} d, second baroclinic Kelvin ${(L1 / c[2] / 86400).toFixed(0)} d`);
}
