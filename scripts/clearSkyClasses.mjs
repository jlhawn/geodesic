// The clear-sky shortwave of a loaded CPU model by surface class, for
// scripts/clearSkyBudget.mjs and scripts/verticalAudit.mjs. Each sunlit
// column with its cloud, deck and cumulus taken away is lit at `times`
// instants spread over day `day` (the day ending `day` days after the
// equinox), with the state's vapour, sea ice, snow and soil held fixed. A
// sea cell with ice is lit twice more, over its open water alone and over
// its ice alone, so that each class's numbers are those of columns over that
// surface; the cell's own column (its area-weighted albedos) gives the
// global and banded rows. Per class: the area share, the surface albedo (the
// sunlight the surface reflects over the sunlight reaching it, direct and
// diffuse, summed over the day), the clear-sky albedo at the top, the
// atmosphere's own share (the same column over a black surface), each beside
// a reference range (REFERENCES) and a verdict. Open water also by the
// sun's cosine, its direct-beam albedo beside the curves of Taylor et al.
// (1996) and of Fresnel reflection from flat water, the range they span
// widened by 0.01 on each side.
import { DAY } from '../js/physics/radiation.module.js';

export const TERMS = ['insolation', 'reflected', 'atmosphere', 'ozone', 'vapour', 'aerosol', 'absorbedSurface', 'down', 'directAlbedo'];
export const TYPES = ['open sea', 'sea ice', 'land', 'ice sheets'];
export const MU_BINS = [0, 0.1, 0.2, 0.4, 0.7, 1];
const SNOWY = 10, DRY_FILL = 0.35, FOREST = 0.5;

// [name, surface range, top-of-atmosphere range, what the surface model lacks when outside]
export const REFERENCES = [
  ['open sea 0-30', null, [0.08, 0.10], ''],
  ['open sea 30-50', null, [0.10, 0.13], ''],
  ['open sea 50-70', null, [0.13, 0.20], ''],
  ['open sea 70-90', null, null, ''],
  ['bare dry soil (v < 0.2, fill < 0.35)', [0.30, 0.40], null, 'one dry soil albedo 0.30, no sand or soil colour'],
  ['bare wet soil (v < 0.2, fill >= 0.35)', [0.10, 0.20], null, ''],
  ['partly vegetated (0.2-0.7)', [0.18, 0.25], null, 'the cover blends bare soil into forest; no grass or crop albedo'],
  ['dense vegetation (v > 0.7)', [0.12, 0.15], null, ''],
  ['thin snow on land (< 10 kg/m2)', null, null, ''],
  ['snow on open land (v < 0.5)', [0.60, 0.85], null, 'one snow albedo 0.55 with no ageing'],
  ['snow under forest (v >= 0.5)', [0.20, 0.35], null, 'snow albedo 0.55 whatever the cover: no masking by forest'],
  ['thin sea ice (< 0.5 m, snow < 10 kg/m2)', [0.20, 0.50], null, 'ice albedo ramps from the water\'s at 0 m to 0.5 at 0.5 m'],
  ['bare sea ice (>= 0.5 m, snow < 10 kg/m2)', [0.50, 0.60], null, 'bare ice 0.5, no ponds'],
  ['snow-covered sea ice (snow >= 10 kg/m2)', [0.80, 0.85], null, 'snow on ice 0.75 with no ageing'],
  ['ice sheets', [0.80, 0.85], null, 'one ice-sheet albedo 0.8'],
];
export const REFERENCE_SOURCES = 'open sea at the top: CERES EBAF clear-sky ocean (approximate); surfaces: textbook ranges (approximate), thin ice from memory';

export function taylorAlbedo(mu) { return 0.037 / (1.1 * mu ** 1.4 + 0.15); }
export function fresnelAlbedo(mu, n = 1.333) {
  if (!(mu > 0)) return 1;
  const i = Math.acos(Math.min(1, mu)), t = Math.asin(Math.sin(i) / n);
  if (i < 1e-6) return ((n - 1) / (n + 1)) ** 2;
  return 0.5 * ((Math.sin(i - t) / Math.sin(i + t)) ** 2 + (Math.tan(i - t) / Math.tan(i + t)) ** 2);
}

export function verdict(value, range) {
  if (!range || !Number.isFinite(value)) return 'n/a';
  const [lo, hi] = range, v = Math.round(value * 1000) / 1000;
  if (v >= lo && v <= hi) return 'matches';
  return value < lo ? `low by ${(lo - value).toFixed(3)}` : `high by ${(value - hi).toFixed(3)}`;
}

export function landClass(land, i, iceSheet) {
  if (iceSheet) return 'ice sheets';
  const v = land.vegetation[i], snow = land.snow[i];
  if (snow >= SNOWY) return v < FOREST ? 'snow on open land (v < 0.5)' : 'snow under forest (v >= 0.5)';
  if (snow > 0) return 'thin snow on land (< 10 kg/m2)';
  if (v < 0.2) return land.soil[i] / land.capacity(i) < DRY_FILL ? 'bare dry soil (v < 0.2, fill < 0.35)' : 'bare wet soil (v < 0.2, fill >= 0.35)';
  return v <= 0.7 ? 'partly vegetated (0.2-0.7)' : 'dense vegetation (v > 0.7)';
}
export function iceClass(h, snow) {
  if (snow >= SNOWY) return 'snow-covered sea ice (snow >= 10 kg/m2)';
  return h < 0.5 ? 'thin sea ice (< 0.5 m, snow < 10 kg/m2)' : 'bare sea ice (>= 0.5 m, snow < 10 kg/m2)';
}
export function seaBand(latDegrees) {
  const a = Math.abs(latDegrees);
  return a < 30 ? 'open sea 0-30' : a < 50 ? 'open sea 30-50' : a < 70 ? 'open sea 50-70' : 'open sea 70-90';
}

export function clearSkyClasses(model, { day, times = 48, ozone = 0.03 }) {
  const { mesh, core, state, radiation, seaIce, land, geography } = model;
  const { K, C } = core.diagnostics;
  const [pi, theta, , surfaceT, q, , ice] = state;
  const bottom = (K - 1) * C, deg = 180 / Math.PI;
  radiation.useCumulus(null, null);
  for (let i = 0; i < C; i++) core.diagnoseColumn(i, pi, theta, q, null);
  const zero = () => Object.fromEntries(TERMS.map((t) => [t, 0]));
  const cells = Array.from({ length: C }, zero);
  const classes = new Map(), types = new Map(), muBins = MU_BINS.slice(1).map(() => ({ ...zero(), mu: 0 }));
  const areaOf = new Map(), typeArea = new Map();
  const accumulate = (target, key, weight, values) => {
    if (!target.has(key)) target.set(key, zero());
    const r = target.get(key);
    for (const t of TERMS) r[t] += weight * values[t];
  };
  const run = (i, beam, adir, adif) => {
    radiation.column(i, pi[i], theta, surfaceT[i], 5, undefined, beam, q[bottom + i], q, null, adir, adif, 1, undefined, 0, 0, 0, 0);
    return { ...radiation.budget };
  };
  const terms = (beam, all, black, adir) => ({
    insolation: beam, reflected: all.reflectedSolar, atmosphere: black.reflectedSolar, ozone: beam * ozone, aerosol: all.aerosolSolar,
    vapour: all.atmosphereSolar - beam * ozone - all.aerosolSolar, absorbedSurface: all.absorbedSolar - all.atmosphereSolar, down: all.surfaceShortwave, directAlbedo: beam * adir,
  });
  const parts = (i) => {
    if (geography.land[i]) {
      const iceSheet = !!(geography.iceSheet && geography.iceSheet[i]);
      return [{ share: 1, name: landClass(land, i, iceSheet), type: iceSheet ? 'ice sheets' : 'land', albedo: () => [land.albedo(i), land.albedo(i)] }];
    }
    const h = ice[i], area = seaIce.cover(i, h), snow = seaIce.snow[i], out = [];
    if (area < 1) out.push({ share: 1 - area, name: seaBand(mesh.latCell[i] * deg), type: 'open sea', water: true, albedo: (mu) => [seaIce.albedo(0, mu), seaIce.albedo(0, null)] });
    if (area > 0) out.push({ share: area, name: iceClass(h, snow), type: 'sea ice', albedo: (mu) => [seaIce.albedo(h, mu, snow, 1), seaIce.albedo(h, null, snow, 1)] });
    return out;
  };
  const cellParts = Array.from({ length: C }, (_, i) => parts(i));
  for (let i = 0; i < C; i++) for (const p of cellParts[i]) {
    areaOf.set(p.name, (areaOf.get(p.name) ?? 0) + mesh.areaCell[i] * p.share);
    typeArea.set(p.type, (typeArea.get(p.type) ?? 0) + mesh.areaCell[i] * p.share);
  }
  for (let n = 0; n < times; n++) {
    radiation.setTime((day - 1) * DAY + (n + 0.5) * DAY / times);
    for (let i = 0; i < C; i++) {
      const mu = radiation.cosZenith(i), beam = radiation.insolation(i);
      if (!(beam > 0)) continue;
      let adir, adif;
      if (geography.land[i]) adir = adif = land.albedo(i);
      else { const h = ice[i], area = seaIce.cover(i, h); adir = seaIce.albedo(h, mu, seaIce.snow[i], area); adif = seaIce.albedo(h, null, seaIce.snow[i], area); }
      const black = run(i, beam, 0, 0), all = run(i, beam, adir, adif), cell = terms(beam, all, black, adir);
      for (const t of TERMS) cells[i][t] += cell[t];
      const list = cellParts[i];
      for (const p of list) {
        const [pd, pf] = p.albedo(mu);
        const values = list.length === 1 ? cell : terms(beam, run(i, beam, pd, pf), black, pd);
        const w = mesh.areaCell[i] * p.share;
        accumulate(classes, p.name, w, values);
        accumulate(types, p.type, w, values);
        if (p.water) {
          const b = muBins[Math.min(MU_BINS.length - 2, MU_BINS.findIndex((edge, k) => mu >= edge && mu < (MU_BINS[k + 1] ?? Infinity)))];
          for (const t of TERMS) b[t] += w * values[t];
          b.mu += w * beam * mu;
        }
      }
    }
  }
  let globe = 0;
  for (let i = 0; i < C; i++) globe += mesh.areaCell[i];
  const describe = (r, area) => ({
    areaShare: area / globe, insolation: r.insolation / times / area,
    surfaceAlbedo: 1 - r.absorbedSurface / r.down, toaAlbedo: r.reflected / r.insolation, atmosphereShare: r.atmosphere / r.insolation,
    directAlbedo: r.directAlbedo / r.insolation, raw: Object.fromEntries(TERMS.map((t) => [t, r[t] / times / area])),
  });
  const classRows = REFERENCES.filter(([name]) => classes.has(name)).map(([name, surfaceRange, toaRange, lacks]) => {
    const d = describe(classes.get(name), areaOf.get(name));
    const surfaceVerdict = verdict(d.surfaceAlbedo, surfaceRange), toaVerdict = verdict(d.toaAlbedo, toaRange);
    const off = (surfaceRange && surfaceVerdict !== 'matches') || (toaRange && toaVerdict !== 'matches');
    return { name, ...d, surfaceRange, toaRange, surfaceVerdict, toaVerdict, note: off && lacks ? lacks : '' };
  });
  const typeRows = TYPES.filter((t) => types.has(t)).map((t) => ({ name: t, ...describe(types.get(t), typeArea.get(t)) }));
  const muRows = muBins.map((b, k) => {
    const mean = b.mu / b.insolation;
    const curves = [taylorAlbedo(mean), fresnelAlbedo(mean)];
    const range = [Math.min(...curves) - 0.01, Math.max(...curves) + 0.01];
    const direct = b.directAlbedo / b.insolation;
    return { name: `open sea, mu ${MU_BINS[k].toFixed(1)}-${MU_BINS[k + 1].toFixed(1)}`, mu: mean, insolationShare: b.insolation, directAlbedo: direct, surfaceAlbedo: 1 - b.absorbedSurface / b.down, toaAlbedo: b.reflected / b.insolation, atmosphereShare: b.atmosphere / b.insolation, taylor: curves[0], fresnel: curves[1], directVerdict: verdict(direct, range) };
  });
  const seaInsolation = muBins.reduce((s, b) => s + b.insolation, 0);
  for (const r of muRows) r.insolationShare /= seaInsolation;
  return { cells, classRows, typeRows, muRows, times, globe };
}

export function classTable(result, f = (x, d) => (Number.isFinite(x) ? x.toFixed(d) : 'n/a')) {
  const lines = [];
  const range = (r) => (r ? `${r[0].toFixed(2)}-${r[1].toFixed(2)}` : '');
  lines.push('class                                      area   surfAlb  ref        verdict          TOA    ref        verdict          atmos');
  for (const r of result.classRows) lines.push(`${r.name.padEnd(42)} ${f(r.areaShare, 3)}  ${f(r.surfaceAlbedo, 3).padStart(6)}  ${range(r.surfaceRange).padEnd(9)}  ${(r.surfaceRange ? r.surfaceVerdict : '').padEnd(15)}  ${f(r.toaAlbedo, 3)}  ${range(r.toaRange).padEnd(9)}  ${(r.toaRange ? r.toaVerdict : '').padEnd(15)}  ${f(r.atmosphereShare, 3)}${r.note ? `  [${r.note}]` : ''}`);
  lines.push('open sea by the sun\'s cosine     mean mu  sun share  direct  Taylor  Fresnel  verdict          surfAlb  TOA    atmos');
  for (const r of result.muRows) lines.push(`${r.name.padEnd(32)} ${f(r.mu, 3)}    ${f(r.insolationShare, 3)}      ${f(r.directAlbedo, 3)}   ${f(r.taylor, 3)}   ${f(r.fresnel, 3)}    ${r.directVerdict.padEnd(15)}  ${f(r.surfaceAlbedo, 3)}    ${f(r.toaAlbedo, 3)}  ${f(r.atmosphereShare, 3)}`);
  return lines;
}
