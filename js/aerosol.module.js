import { AEROSOL_SPECTRA } from './physics/aerosolSpectra.module.js';
import { monthWeights } from './ocean/climatology.module.js';

/*
 * The MACv2 aerosol climatology (Kinne 2019, Tellus B 71, 1623639; reference year 2005) that
 * radiation option aerosol 'mac' puts in the shortwave, in data/aerosol_mac_2deg.bin
 * (scripts/packAerosol.py, which documents the layout): for each month and cell of a 2 degree
 * grid, eight bytes, the 550 nm optical depth, 1 - single scattering albedo and asymmetry of the
 * fine mode and the coarse mode and the mean height above the ground of each mode's aerosol.
 * aerosolFields blends two months of them onto a mesh's cells, 8 Float32 per cell
 * (AEROSOL_FIELD, [field][cell]); both engines are fed by it, and what they do with a cell's
 * eight numbers is the arithmetic below.
 *
 * Each mode's optics at the 12 wavelengths of AEROSOL_SPECTRA follow from its 550 nm triple
 * (tau, x = 1 - omega, g) by modeBand's constants, and a cell's two modes make
 *  - the sequential scheme's two slabs, one per stream: streamOptics gives the mixture's depth,
 *    albedo and asymmetry in the visible (0.297-0.7 um) and in the rest of the beam, and the
 *    share of the absorption that is the fine mode's, which weighs the two profiles;
 *  - the layered scheme's layers: a slice's aerosol in a layer is the sum over the modes of its
 *    depth, absorption depth, scattering depth and asymmetry-weighted scattering depth times
 *    the mode's share of its column in the layer.
 * A mode's shares come from profileFractions: a layer of constant density from the ground up to
 * the depth `MIXED_DEPTH` (at most MIXED_CAP times the mean height), an exponential above it,
 * the exponential's scale height such that the mean height is the cell's, on the layers' own
 * heights above the ground (their interfaces the midpoints of the mid-layer heights, the top
 * one infinitely high, the lowest the ground).
 * aerosolScale of the radiation multiplies the depths.
 */
export const DEFAULT_AEROSOL = 'mac';
export const climatologyOn = (radiation) => (radiation?.aerosol ?? DEFAULT_AEROSOL) === 'mac';
export const AEROSOL_FILE = 'data/aerosol_mac_2deg.bin';
export const AEROSOL_FIELD = { fineDepth: 0, coarseDepth: 1, fineAbsorption: 2, coarseAbsorption: 3, fineAsymmetry: 4, coarseAsymmetry: 5, fineHeight: 6, coarseHeight: 7 };
export const AEROSOL_FIELDS = 8;
export const WAVELENGTHS = AEROSOL_SPECTRA.wavelengths.length;
export const MIXED_DEPTH = [1.5, 1.0];
export const HEIGHT_RANGE = [0.3, 6];
export const MIXED_CAP = 1.8;
export const ASYMMETRY_CEILING = 0.99;
export const AEROSOL_RANGES = [[0, 1.2, true], [0, 1.2, true], [0, 0.3, false], [0, 0.3, false], [0.45, 0.85, false], [0.45, 0.85, false], [0, 8, false], [0, 8, false]];
const MAGIC = 'AERM', DAY = 86400, EQUINOX_DAY_OF_YEAR = 78;

export function decodeAerosol(bytes) {
  const view = ArrayBuffer.isView(bytes) ? new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength) : new DataView(bytes);
  const magic = String.fromCharCode(...[0, 1, 2, 3].map((n) => view.getUint8(n)));
  if (magic !== MAGIC) throw new Error(`not an aerosol climatology (${magic})`);
  const offset = view.getUint32(4, true), nLon = view.getUint16(8, true), nLat = view.getUint16(10, true), nMonth = view.getUint16(12, true), nField = view.getUint16(14, true);
  if (nField !== AEROSOL_FIELDS) throw new Error(`an aerosol climatology of ${nField} fields, not ${AEROSOL_FIELDS}`);
  const [lon0, dLon, lat0, dLat] = [0, 1, 2, 3].map((n) => view.getFloat32(16 + 4 * n, true));
  const centers = Float64Array.from({ length: nMonth }, (_, m) => view.getFloat32(32 + 4 * m, true));
  let at = 32 + 4 * nMonth;
  const ranges = Array.from({ length: nField }, () => { const range = { low: view.getFloat32(at, true), high: view.getFloat32(at + 4, true), square: view.getUint32(at + 8, true) === 1 }; at += 12; return range; });
  const textLength = view.getUint16(at, true);
  const source = String.fromCharCode(...Array.from({ length: textLength }, (_, n) => view.getUint8(at + 2 + n)));
  const count = nField * nMonth * nLat * nLon;
  if (view.byteLength !== offset + count) throw new Error(`an aerosol climatology of ${view.byteLength} bytes, not ${offset + count}`);
  const values = new Uint8Array(view.buffer, view.byteOffset + offset, count);
  const lookup = new Float64Array(nField * 256);
  ranges.forEach(({ low, high, square }, f) => { for (let q = 0; q < 256; q++) lookup[f * 256 + q] = low + (high - low) * (square ? (q / 255) ** 2 : q / 255); });
  return { nLon, nLat, nMonth, nField, lon0, dLon, lat0, dLat, centers, ranges, values, lookup, source };
}

/*
 * The file for fields given as [field][month][lat][lon] flattened real values (rows south to
 * north, columns west to east from lon0) on a grid of nLon x nLat cells and the month centers
 * in days of a 365-day year from January 1.
 */
export function encodeAerosol({ nLon, nLat, lon0, dLon, lat0, dLat, centers, fields, source = '', ranges = AEROSOL_RANGES }) {
  const nMonth = centers.length, nField = AEROSOL_FIELDS, plane = nMonth * nLat * nLon;
  let head = 32 + 4 * nMonth + 12 * nField + 2 + source.length;
  head += (8 - (head % 8)) % 8;
  const buffer = new ArrayBuffer(head + nField * plane), view = new DataView(buffer);
  for (let n = 0; n < 4; n++) view.setUint8(n, MAGIC.charCodeAt(n));
  view.setUint32(4, head, true); view.setUint16(8, nLon, true); view.setUint16(10, nLat, true); view.setUint16(12, nMonth, true); view.setUint16(14, nField, true);
  [lon0, dLon, lat0, dLat].forEach((value, n) => view.setFloat32(16 + 4 * n, value, true));
  centers.forEach((c, m) => view.setFloat32(32 + 4 * m, c, true));
  let at = 32 + 4 * nMonth;
  ranges.forEach(([low, high, square]) => { view.setFloat32(at, low, true); view.setFloat32(at + 4, high, true); view.setUint32(at + 8, square ? 1 : 0, true); at += 12; });
  view.setUint16(at, source.length, true);
  for (let n = 0; n < source.length; n++) view.setUint8(at + 2 + n, source.charCodeAt(n));
  ranges.forEach(([low, high, square], f) => {
    for (let n = 0; n < plane; n++) {
      const s = Math.min(1, Math.max(0, (fields[f][n] - low) / (high - low)));
      view.setUint8(head + f * plane + n, Math.round(255 * (square ? Math.sqrt(s) : s)));
    }
  });
  return new Uint8Array(buffer);
}

export function loadAerosol(data = null) {
  if (data && data.lookup) return data;
  if (data) return decodeAerosol(data);
  const fs = globalThis.process?.getBuiltinModule?.('node:fs');
  const url = new URL(`../${AEROSOL_FILE}`, import.meta.url);
  if (!fs) throw new Error(`aerosol 'mac' needs the climatology's bytes: pass them as aerosolData`);
  if (!fs.existsSync(url)) throw new Error(`aerosol 'mac' needs ${AEROSOL_FILE}, which is not in the repository: build it with scripts/packAerosol.py (docs/aerosol-climatology.md)`);
  return decodeAerosol(fs.readFileSync(url));
}

export const aerosolDay = (time) => Math.floor(time / DAY);
export const aerosolDayOfYear = (day) => (((day + 0.5 + EQUINOX_DAY_OF_YEAR) % 365) + 365) % 365;

/*
 * The eight fields of the climatology on a mesh's cells on a day of the year (days of a
 * 365-day year from January 1): bilinear in the grid, linear between the two months whose
 * centers bracket the day, in the stored variables. Rows beyond the grid's first and last
 * take theirs, and columns wrap.
 */
export function aerosolFields(mesh, data, dayOfYear, out = new Float32Array(AEROSOL_FIELDS * mesh.nCells)) {
  const { nLon, nLat, nMonth, lon0, dLon, lat0, dLat, values, lookup } = data, C = mesh.nCells, deg = 180 / Math.PI, plane = nMonth * nLat * nLon;
  const months = monthWeights(data, dayOfYear);
  for (let i = 0; i < C; i++) {
    const y = Math.min(nLat - 1, Math.max(0, (mesh.latCell[i] * deg - lat0) / dLat)), r0 = Math.min(nLat - 2, Math.floor(y)), fy = y - r0;
    const x = (mesh.lonCell[i] * deg - lon0) / dLon, fx = x - Math.floor(x), c0 = ((Math.floor(x) % nLon) + nLon) % nLon, c1 = (c0 + 1) % nLon;
    const w00 = (1 - fy) * (1 - fx), w01 = (1 - fy) * fx, w10 = fy * (1 - fx), w11 = fy * fx;
    const t00 = r0 * nLon + c0, t01 = r0 * nLon + c1, t10 = (r0 + 1) * nLon + c0, t11 = (r0 + 1) * nLon + c1;
    for (let f = 0; f < AEROSOL_FIELDS; f++) {
      let sum = 0;
      for (const [m, weight] of months) {
        const base = f * plane + m * nLat * nLon, table = f * 256;
        sum += weight * (w00 * lookup[table + values[base + t00]] + w01 * lookup[table + values[base + t01]] + w10 * lookup[table + values[base + t10]] + w11 * lookup[table + values[base + t11]]);
      }
      out[f * C + i] = sum;
    }
  }
  return out;
}

/*
 * A sampler of the climatology on a mesh: `fields` is refreshed by update(time) when the model
 * day (the day since the March equinox, 78 days into the year) differs from that of the last
 * update, which returns whether it did.
 */
export function createAerosolSampler(mesh, data) {
  const fields = new Float32Array(AEROSOL_FIELDS * mesh.nCells);
  let day = null;
  return {
    fields, data,
    get day() { return day; },
    update(time) {
      const now = aerosolDay(time);
      if (now === day) return false;
      aerosolFields(mesh, data, aerosolDayOfYear(now), fields);
      day = now;
      return true;
    },
  };
}

/*
 * The optics of one mode in wavelength b of AEROSOL_SPECTRA from its 550 nm depth, x = 1 - omega
 * and asymmetry: out.tau, out.x and out.g.
 */
const BAND = { tau: 0, x: 0, g: 0 };
export function modeBand(mode, b, tau, x, g, out) {
  const s = mode === 0 ? AEROSOL_SPECTRA.fine : AEROSOL_SPECTRA.coarse;
  out.tau = s.r[b] * tau;
  out.x = Math.min(1, Math.max(0, s.a[b] + s.c[b] * x));
  out.g = Math.min(ASYMMETRY_CEILING, Math.max(0, s.d[b] + s.e[b] * g));
}

/*
 * A cell's aerosol: from the two modes' 550 nm depths, 1 - omega and asymmetries (tau, x, g,
 * each [fine, coarse]) fills out.tau, out.x and out.g [mode * WAVELENGTHS + b] and, for the
 * sequential scheme's visible stream and its rest, out.tauV, out.omegaV, out.gV, out.fineV and
 * out.tauN, out.omegaN, out.gN, out.fineN: the Planck-weighted depth, its single scattering
 * albedo and asymmetry (1 and 0 without depth) and the fine mode's share of the absorption
 * depth (of the depth where nothing absorbs).
 */
export function streamOptics(tau, x, g, out) {
  const one = BAND, { visible, nearInfrared } = AEROSOL_SPECTRA;
  let tV = 0, sV = 0, pV = 0, tN = 0, sN = 0, pN = 0, aV0 = 0, aV1 = 0, dV0 = 0, dV1 = 0, aN0 = 0, aN1 = 0, dN0 = 0, dN1 = 0;
  for (let m = 0; m < 2; m++) {
    for (let b = 0; b < WAVELENGTHS; b++) {
      modeBand(m, b, tau[m], x[m], g[m], one);
      const n = m * WAVELENGTHS + b, scattering = one.tau * (1 - one.x);
      out.tau[n] = one.tau; out.x[n] = one.x; out.g[n] = one.g;
      const wv = visible[b], wn = nearInfrared[b];
      if (wv > 0) {
        tV += wv * one.tau; sV += wv * scattering; pV += wv * scattering * one.g;
        if (m === 0) { aV0 += wv * one.tau * one.x; dV0 += wv * one.tau; } else { aV1 += wv * one.tau * one.x; dV1 += wv * one.tau; }
      }
      if (wn > 0) {
        tN += wn * one.tau; sN += wn * scattering; pN += wn * scattering * one.g;
        if (m === 0) { aN0 += wn * one.tau * one.x; dN0 += wn * one.tau; } else { aN1 += wn * one.tau * one.x; dN1 += wn * one.tau; }
      }
    }
  }
  out.tauV = tV; out.omegaV = tV > 0 ? sV / tV : 1; out.gV = sV > 0 ? pV / sV : 0;
  out.tauN = tN; out.omegaN = tN > 0 ? sN / tN : 1; out.gN = sN > 0 ? pN / sN : 0;
  const absorbingV = aV0 + aV1 > 0, absorbingN = aN0 + aN1 > 0;
  const fineV = absorbingV ? aV0 : dV0, totalV = absorbingV ? aV0 + aV1 : dV0 + dV1, fineN = absorbingN ? aN0 : dN0, totalN = absorbingN ? aN0 + aN1 : dN0 + dN1;
  out.fineV = totalV > 0 ? fineV / totalV : 0;
  out.fineN = totalN > 0 ? fineN / totalN : 0;
  return out;
}

export const aerosolTable = () => ({ tau: new Float64Array(2 * WAVELENGTHS), x: new Float64Array(2 * WAVELENGTHS), g: new Float64Array(2 * WAVELENGTHS), tauV: 0, omegaV: 1, gV: 0, fineV: 0, tauN: 0, omegaN: 1, gN: 0, fineN: 0 });

/*
 * A mode's share of its column in each layer, out[k] for the K layers from the top: `mean` its
 * mean height above the ground (km, within HEIGHT_RANGE), `depth` the depth of its mixed
 * layer (km, at most MIXED_CAP times the mean height, so that the exponential above it has a scale
 * height) and zMid[k] the layers' heights above the ground (km). The shares sum to 1.
 */
export function profileFractions(mean, depth, zMid, out) {
  const K = zMid.length, mu = Math.min(HEIGHT_RANGE[1], Math.max(HEIGHT_RANGE[0], mean)), D = Math.min(depth, MIXED_CAP * mu);
  const H = 0.5 * ((mu - D) + Math.sqrt(mu * mu + 2 * mu * D - D * D)), total = D + H;
  let upper = 0;
  for (let k = 0; k < K; k++) {
    const height = k === K - 1 ? 0 : 0.5 * (zMid[k] + zMid[k + 1]);
    const lower = height < D ? 1 - height / total : H / total * Math.exp(-(height - D) / H);
    out[k] = lower - upper;
    upper = lower;
  }
  return out;
}
