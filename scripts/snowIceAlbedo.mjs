// The surface albedo of snow, sea ice and the ice sheets in a saved state,
// on the CPU engine, by the state that sets it:
//   node scripts/snowIceAlbedo.mjs <state.bin>
// Each albedo is the surface's own (land.albedo, and for sea ice the ice
// part of its cell alone, under the direct beam) weighted by the sunlight
// reaching the top of each column at TIMES (48) instants spread over day
// DAY (the state's day by default; the day ending DAY days after the
// equinox). Rows: snow-covered land (at least SNOWY, 10 kg/m2, off the ice
// sheets) by latitude band and by the standing cover the masking reads
// (the vegetation cover on code without one), with the cover, the standing
// cover and the per-cell snow albedo where the land has them, the mean
// pace of the cold snow's ageing, and the share of its snow within
// 2 K of melting; the land's snow cover by 5-degree band (the share of the
// land, off the ice sheets, under at least 1 kg/m2) and the snow line, the
// lowest band at least half covered; the cover of the boreal belt (50-70N);
// sea ice by hemisphere, snow load and skin temperature, thin ice (< 0.5 m)
// apart; the ice sheets by hemisphere and skin temperature. Each snow row
// also gives the state's last-day precipitation where the skin is below
// freezing (mm/d, a proxy for snowfall) and the albedo at which the row's
// mean snowfall, refreshing the snow, balances the cold snow's ageing
// (agedSnowAlbedo and refreshedSnowAlbedo in js/physics/ice.module.js),
// 0.85 - 0.008 p / (fall / 10 kg/m2) and not below the floor (0.5 on land,
// 0.7 on ice), once with p the row's mean slowing of the ageing in the
// cold (eq) and once at the plain 0.008 a day (eqPlain).
// LAND (JSON) passes options to the land, ICE (JSON) to the sea ice.
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { DAY } from '../js/physics/radiation.module.js';
import { SNOW_AGEING } from '../js/physics/ice.module.js';

const FILE = process.argv[2];
if (!FILE) { console.error('usage: node scripts/snowIceAlbedo.mjs <state.bin>'); process.exit(1); }
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const TIMES = Number(process.env.TIMES ?? 48), LAND = JSON.parse(process.env.LAND ?? '{}'), ICE = JSON.parse(process.env.ICE ?? '{}'), AT = Number(process.env.DAY ?? saved.day);
const SNOWY = 10, MELTING = 273.15;
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(saved.N), { topography, levels: savedLevels(saved), ocean: false, land: LAND, ice: ICE });
const { mesh, state, seaIce, land, geography, radiation } = model;
const C = mesh.nCells, deg = 180 / Math.PI;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
land.load(saved.land, state[6]);
const [, , , surfaceT, , , ice] = state;
const fall = Float64Array.from({ length: C }, (_, i) => (surfaceT[i] < MELTING && saved.convectiveRain ? saved.convectiveRain[i] + saved.largeScaleRain[i] : 0));
const iceSheet = (i) => !!(geography.iceSheet && geography.iceSheet[i]);
const standing = land.canopy ?? land.vegetation, snowAlbedo = land.snowAlbedo ?? null;
const pace = (i) => Math.min(1, Math.exp(SNOW_AGEING.ageingActivation * (surfaceT[i] - MELTING) / (MELTING * surfaceT[i])));
const balance = (fallRate, p, floor) => (fallRate > 0 ? Math.max(floor, SNOW_AGEING.freshSnowAlbedo - SNOW_AGEING.coldSnowAgeing * p / (fallRate / SNOW_AGEING.refreshSnowfall)) : floor);
const iceAlbedo = (i, mu) => seaIce.albedo(ice[i], mu, seaIce.snow[i], 1, surfaceT[i], seaIce.snowAlbedo ? seaIce.snowAlbedo[i] : undefined);

const sun = new Float64Array(C), reflectedLand = new Float64Array(C), reflectedIce = new Float64Array(C);
for (let n = 0; n < TIMES; n++) {
  radiation.setTime((AT - 1) * DAY + (n + 0.5) * DAY / TIMES);
  for (let i = 0; i < C; i++) {
    const beam = radiation.insolation(i);
    if (!(beam > 0)) continue;
    sun[i] += beam / TIMES;
    if (geography.land[i]) reflectedLand[i] += beam * land.albedo(i) / TIMES;
    else if (ice[i] > 0) reflectedIce[i] += beam * iceAlbedo(i, radiation.cosZenith(i)) / TIMES;
  }
}

const rows = new Map();
const add = (key, i, reflected, extra = {}) => {
  if (!rows.has(key)) rows.set(key, { area: 0, sun: 0, reflected: 0, fall: 0, n: {}, sums: {} });
  const r = rows.get(key), a = mesh.areaCell[i] * (extra.share ?? 1);
  r.area += a; r.sun += a * sun[i]; r.reflected += a * reflected; r.fall += a * fall[i];
  for (const [k, v] of Object.entries(extra)) if (k !== 'share' && Number.isFinite(v)) { r.sums[k] = (r.sums[k] ?? 0) + a * v; r.n[k] = (r.n[k] ?? 0) + a; }
};
let globe = 0;
for (let i = 0; i < C; i++) globe += mesh.areaCell[i];
const band = (lat, edges) => { const a = Math.abs(lat), k = edges.findIndex((e, j) => a >= e && a < (edges[j + 1] ?? 91)); return k < 0 ? null : `${lat >= 0 ? 'N' : 'S'} ${edges[k]}-${edges[k + 1] ?? 90}`; };
const coverBins = [0, 0.1, 0.3, 0.5, 0.7, 1.0001];
const coverBin = (v) => { const k = coverBins.findIndex((e, j) => v >= e && v < coverBins[j + 1]); return `v ${coverBins[k].toFixed(1)}-${Math.min(1, coverBins[k + 1]).toFixed(1)}`; };
const snowBin = (s) => (s < 1 ? 'snow < 1' : s < 10 ? 'snow 1-10' : s < 20 ? 'snow 10-20' : 'snow >= 20');
const tBin = (t) => (t < MELTING - 10 ? 'T < -10' : t < MELTING - 2 ? 'T -10..-2' : t < MELTING - 1 ? 'T -2..-1' : 'T >= -1');
const snowLine = new Map(), boreal = { area: 0, v: 0, forest: 0, snow: 0 };

for (let i = 0; i < C; i++) {
  const lat = mesh.latCell[i] * deg;
  if (geography.land[i]) {
    if (iceSheet(i)) { add(`ice sheet ${lat >= 0 ? 'N' : 'S'} ${tBin(surfaceT[i]).replace('T -2..-1', 'T >= -2').replace('T >= -1', 'T >= -2')}`, i, reflectedLand[i], { snow: land.snow[i], snowAlbedo: snowAlbedo ? snowAlbedo[i] : NaN, pace: pace(i) }); continue; }
    const b5 = band(lat, [30, 35, 40, 45, 50, 55, 60, 65, 70, 75]);
    if (b5) { if (!snowLine.has(b5)) snowLine.set(b5, { area: 0, covered: 0 }); const s = snowLine.get(b5); s.area += mesh.areaCell[i]; if (land.snow[i] >= 1) s.covered += mesh.areaCell[i]; }
    if (lat >= 50 && lat < 70) { const a = mesh.areaCell[i]; boreal.area += a; boreal.v += a * land.vegetation[i]; boreal.forest += land.vegetation[i] >= 0.5 ? a : 0; boreal.snow += land.snow[i] > 0 ? a : 0; }
    if (land.snow[i] < SNOWY) continue;
    const extra = { v: land.vegetation[i], standing: standing[i], snowAlbedo: snowAlbedo ? snowAlbedo[i] : NaN, melting: surfaceT[i] >= MELTING - 2 ? 1 : 0, pace: pace(i) };
    const b = band(lat, [0, 30, 40, 50, 60, 70]);
    add(`land snow ${b}`, i, reflectedLand[i], extra);
    add(`land snow ${coverBin(standing[i])}`, i, reflectedLand[i], extra);
    add('land snow, all', i, reflectedLand[i], extra);
  } else if (ice[i] > 0) {
    const share = seaIce.cover(i, ice[i]), hemi = lat >= 0 ? 'N' : 'S', extra = { share, snow: seaIce.snow[i], snowAlbedo: seaIce.snowAlbedo ? seaIce.snowAlbedo[i] : NaN, h: ice[i], pace: pace(i) };
    const key = ice[i] < 0.5 ? `sea ice ${hemi} thin (< 0.5 m)` : `sea ice ${hemi} ${snowBin(seaIce.snow[i]).padEnd(10)} ${tBin(surfaceT[i])}`;
    add(key, i, reflectedIce[i], extra);
    add(`sea ice ${hemi}, all`, i, reflectedIce[i], extra);
  }
}

const f = (x, d = 3) => (Number.isFinite(x) ? x.toFixed(d) : 'n/a');
console.log(`snow and ice albedo of ${FILE.split('/').pop()} (N=${saved.N}, day ${saved.day}) lit over day ${AT} at ${TIMES} instants; LAND ${JSON.stringify(LAND)}; ICE ${JSON.stringify(ICE)}; the state ${land.snowAlbedo && saved.land.snowAlbedo ? 'carries' : 'does not carry'} a snow albedo${land.canopy ? `, ${saved.land.canopy ? 'carries' : 'does not carry'} a standing cover` : ''}`);
console.log('row                                        area     sun W/m2  albedo  fall mm/d  means');
const order = [...rows.keys()].sort();
for (const key of order) {
  const r = rows.get(key), fallRate = r.fall / r.area;
  const floor = key.startsWith('sea ice') ? 0.7 : 0.5, p = r.n.pace ? r.sums.pace / r.n.pace : NaN;
  const means = Object.keys(r.sums).map((k) => `${k} ${f(r.sums[k] / r.n[k], k === 'snow' ? 0 : 2)}`).join(', ') + (r.n.pace ? `, eq ${f(balance(fallRate, p, floor), 2)}, eqPlain ${f(balance(fallRate, 1, floor), 2)}` : '');
  console.log(`${key.padEnd(42)} ${f(r.area / globe, 4)}  ${f(r.sun / r.area, 1).padStart(7)}  ${r.sun > 0 ? f(r.reflected / r.sun) : ' n/a '}  ${f(fallRate, 2).padStart(9)}  ${means}`);
}
console.log('land snow cover by band (share of the land off the ice sheets under at least 1 kg/m2):');
for (const hemi of ['N', 'S']) {
  const bands = [...snowLine.keys()].filter((k) => k.startsWith(hemi)).sort((a, b) => parseFloat(a.slice(2)) - parseFloat(b.slice(2)));
  const line = bands.find((k) => snowLine.get(k).covered / snowLine.get(k).area >= 0.5);
  console.log(`  ${hemi}: ${bands.map((k) => `${k.slice(2)} ${f(snowLine.get(k).covered / snowLine.get(k).area, 2)}`).join(', ')}; snow line ${line ? line.slice(2) : 'none'}`);
}
console.log(`boreal belt 50-70N land: mean cover ${f(boreal.v / boreal.area)}, share with cover >= 0.5 ${f(boreal.forest / boreal.area)}, share under snow ${f(boreal.snow / boreal.area)}`);
