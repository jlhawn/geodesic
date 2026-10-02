// The clear-sky shortwave budget of a saved state on the CPU engine, term by
// term, by 10-degree band and by surface type:
//   node scripts/clearSkyBudget.mjs <state.bin>
// Each column with its cloud, deck and cumulus taken away is lit at TIMES
// (48) instants spread over day DAY (the day ending DAY days after the
// equinox, as the spin-up's day lines count; the state's day by default),
// with the state's vapour, sea ice and snow held fixed, and its sunlight
// split into what leaves the top (atmos, the same column's reflection over
// a black surface, and surf, the rest: the surface's reflection as seen
// from the top), what ozone, vapour and aerosol absorb, and what the
// surface absorbs (sfcabs). Fluxes are day means over the row's area
// (W/m2); albedo is the clear-sky reflected over the incoming, surfAlb the
// direct-beam surface albedo weighted by the incoming, and the last column
// the surface's reflection seen from the top over surfAlb times the
// incoming. RADIATION (JSON) passes options to the radiation.
// Surface types: open sea, sea ice (the ice-covered share of sea cells),
// land, land ice (ice sheets). A sea cell's surface reflection is
// attributed to its open water and its ice in proportion to area times
// albedo, its atmosphere's to their areas.
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { DAY } from '../js/physics/radiation.module.js';

const FILE = process.argv[2];
if (!FILE) { console.error('usage: node scripts/clearSkyBudget.mjs <state.bin>'); process.exit(1); }
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const TIMES = Number(process.env.TIMES ?? 48), RADIATION = JSON.parse(process.env.RADIATION ?? '{}'), AT = Number(process.env.DAY ?? saved.day);
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(saved.N), { topography, levels: savedLevels(saved), ocean: false, radiation: { clearSkyPass: true, ...RADIATION } });
const { mesh, core, state, radiation, seaIce, land, geography } = model;
const { K, C } = core.diagnostics;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
land.load({ soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow), ...(saved.land.vegetation ? { vegetation: Float64Array.from(saved.land.vegetation) } : {}), ...(saved.land.surface ? { surface: Float64Array.from(saved.land.surface) } : {}) }, state[6]);
radiation.useCumulus(null, null);
const [pi, theta, , surfaceT, q, , ice] = state;
const bottom = (K - 1) * C, deg = 180 / Math.PI;
for (let i = 0; i < C; i++) core.diagnoseColumn(i, pi, theta, q, null);

const TERMS = ['insolation', 'reflected', 'atmosphere', 'surface', 'ozone', 'vapour', 'aerosol', 'absorbedSurface', 'albedoWeight'];
const TYPES = ['open sea', 'sea ice', 'land', 'land ice'];
const rows = new Map();
const add = (key, area, values) => {
  if (!rows.has(key)) rows.set(key, { area: 0, ...Object.fromEntries(TERMS.map((t) => [t, 0])) });
  const r = rows.get(key);
  r.area += area;
  for (const t of TERMS) r[t] += area * values[t];
};
const sums = Array.from({ length: C }, () => Object.fromEntries(TERMS.map((t) => [t, 0])));
const shares = Array.from({ length: C }, () => new Float64Array(4)), reflShares = Array.from({ length: C }, () => new Float64Array(4));
const OZONE = RADIATION.ozoneAbsorption ?? 0.03;
function run(i, beam, adir, adif) {
  radiation.column(i, pi[i], theta, surfaceT[i], 5, undefined, beam, q[bottom + i], q, null, adir, adif, 1, undefined, 0, 0, 0, 0);
  return { ...radiation.budget };
}
for (let n = 0; n < TIMES; n++) {
  radiation.setTime((AT - 1) * DAY + (n + 0.5) * DAY / TIMES);
  for (let i = 0; i < C; i++) {
    const mu = radiation.cosZenith(i), beam = radiation.insolation(i);
    if (!(beam > 0)) continue;
    let adir, adif, iceShare = 0, landType = -1, waterDir = 0, iceDir = 0;
    if (geography.land[i]) { adir = adif = land.albedo(i); landType = geography.iceSheet[i] ? 3 : 2; }
    else {
      const h = ice[i], area = seaIce.cover(i, h);
      adir = seaIce.albedo(h, mu, seaIce.snow[i], area); adif = seaIce.albedo(h, null, seaIce.snow[i], area);
      iceShare = area; waterDir = seaIce.albedo(0, mu); iceDir = waterDir + seaIce.albedoContrast(h, mu, seaIce.snow[i]);
    }
    const all = run(i, beam, adir, adif), black = run(i, beam, 0, 0);
    const s = sums[i];
    s.insolation += beam;
    s.reflected += all.reflectedSolar;
    s.atmosphere += black.reflectedSolar;
    s.surface += all.reflectedSolar - black.reflectedSolar;
    s.ozone += beam * OZONE;
    s.aerosol += all.aerosolSolar ?? 0;
    s.vapour += all.atmosphereSolar - beam * OZONE - (all.aerosolSolar ?? 0);
    s.absorbedSurface += all.absorbedSolar - all.atmosphereSolar;
    s.albedoWeight += beam * adir;
    if (landType >= 0) { shares[i][landType] += beam; reflShares[i][landType] += all.reflectedSolar - black.reflectedSolar; }
    else {
      shares[i][0] += beam * (1 - iceShare); shares[i][1] += beam * iceShare;
      const w = (1 - iceShare) * waterDir, v = iceShare * iceDir, surf = all.reflectedSolar - black.reflectedSolar;
      reflShares[i][0] += w + v > 0 ? surf * w / (w + v) : surf; reflShares[i][1] += w + v > 0 ? surf * v / (w + v) : 0;
    }
  }
}
for (let i = 0; i < C; i++) {
  const a = mesh.areaCell[i], s = Object.fromEntries(TERMS.map((t) => [t, sums[i][t] / TIMES]));
  add('global', a, s);
  const lat = mesh.latCell[i] * deg, b = Math.min(80, Math.floor(lat / 10) * 10);
  add(`${b >= 0 ? `${b}N` : `${-b}S`}..${b + 10 >= 0 ? `${b + 10}N` : `${-(b + 10)}S`}`.replace('0N..', '0..').replace(/^0\.\./, 'EQ..'), a, s);
  if (lat >= -30 && lat <= 30) add('30S-30N', a, s);
  const total = shares[i].reduce((x, y) => x + y, 0);
  if (!(total > 0)) {
    const t = geography.land[i] ? (geography.iceSheet[i] ? 3 : 2) : ice[i] > 0 ? 1 : 0;
    add(TYPES[t], a, s); continue;
  }
  for (let t = 0; t < 4; t++) {
    const f = shares[i][t] / total;
    if (!(f > 0)) continue;
    const part = { ...Object.fromEntries(TERMS.map((x) => [x, f * s[x]])), surface: reflShares[i][t] / TIMES, reflected: f * s.atmosphere + reflShares[i][t] / TIMES };
    add(TYPES[t], a * f, Object.fromEntries(TERMS.map((x) => [x, part[x] / f])));
  }
}
const f = (x, d = 1) => (Number.isFinite(x) ? x.toFixed(d) : 'n/a');
console.log(`clear-sky shortwave budget of ${FILE.split('/').pop()} (N=${saved.N}, K=${K}) lit over day ${AT} at ${TIMES} instants; RADIATION ${JSON.stringify(RADIATION)}`);
console.log('row          insol  refl  atmos  surf  ozone vapour aeros sfcabs  albedo  atmos/insol  surfAlb  surf/(surfAlb*insol)');
const order = ['global', '30S-30N', ...[...rows.keys()].filter((k) => k.includes('..')).sort((x, y) => parseFloat(x) * (x.includes('S') ? -1 : 1) - parseFloat(y) * (y.includes('S') ? -1 : 1)), ...TYPES];
for (const key of order) {
  const r = rows.get(key);
  if (!r) continue;
  const m = (t) => r[t] / r.area, alb = m('albedoWeight') / m('insolation');
  console.log(`${key.padEnd(11)} ${f(m('insolation')).padStart(6)} ${f(m('reflected')).padStart(5)} ${f(m('atmosphere')).padStart(6)} ${f(m('surface')).padStart(5)} ${f(m('ozone')).padStart(6)} ${f(m('vapour')).padStart(6)} ${f(m('aerosol')).padStart(5)} ${f(m('absorbedSurface')).padStart(6)}  ${f(m('reflected') / m('insolation'), 3)}   ${f(m('atmosphere') / m('insolation'), 3)}       ${f(alb, 3)}    ${f(m('surface') / (alb * m('insolation')), 3)}${TYPES.includes(key) ? `   (area share ${f(r.area / rows.get('global').area, 3)})` : ''}`);
}
