/*
 * One day of the polar surface and cloud numbers that spinup.mjs writes with
 * POLAR_LINES, as a JSON object. `values` are the day's means from the
 * interval-means recorder (concentration, ice, surfaceT, shortwave, cloud,
 * t2m, sensibleHeat, latentHeat); snow and snowAlbedo the day-end snow on the
 * ice (kg/m² of ice) and its albedo; absorbed, atmosphere and clear the day's
 * sums over `steps` of the sunlight absorbed by the column, by the atmosphere
 * and by the clear-sky column.
 *
 * Zones: the Arctic's pack (day-mean concentration at least 0.8) and marginal
 * ice (0.15 to 0.8) north of 50N, the central Arctic (sea north of 80N) and
 * the Beaufort Sea box of SHEBA's drift (72–80N, 180–130W). In each: the cell
 * area (Mkm²), the area-weighted concentration and downward sunlight, the
 * ice-area-weighted thickness, snow, skin temperature, 2 m temperature, the
 * share of the ice area whose skin is within 2 K of melting (wet), within
 * 0.1 K (melting) and under 1 kg/m² of snow (bare), the snow-weighted snow
 * albedo, and the effective surface albedo, 1 − absorbed / downward. The
 * shortwave cloud effect at the top by 10° band (all sky less clear sky
 * absorbed), the column cloud water of 70–90N and the downward sunlight over
 * the iced cells (concentration at least 0.15) north of 50N.
 */
import { MELTING_POINT } from '../js/physics/ice.module.js';

const deg = 180 / Math.PI;
const BANDS = [[40, 50], [50, 60], [60, 70], [70, 80], [80, 90], [-70, -60], [-90, -70]];
const round = (x, digits) => (Number.isFinite(x) ? Number(x.toFixed(digits)) : null);

export function polarLine({ mesh, land, day, steps, values, snow, snowAlbedo, absorbed, atmosphere, clear }) {
  const C = mesh.nCells, lat = (i) => mesh.latCell[i] * deg, lon = (i) => ((mesh.lonCell[i] * deg + 540) % 360) - 180;
  const { concentration, ice, surfaceT, shortwave, cloud, t2m } = values;
  const zones = {
    pack: (i) => lat(i) > 50 && concentration[i] >= 0.8,
    marginal: (i) => lat(i) > 50 && concentration[i] >= 0.15 && concentration[i] < 0.8,
    central: (i) => lat(i) >= 80,
    beaufort: (i) => lat(i) >= 72 && lat(i) < 80 && lon(i) >= -180 && lon(i) <= -130,
  };
  const line = { day };
  for (const [name, inZone] of Object.entries(zones)) {
    const s = { a: 0, ac: 0, h: 0, snow: 0, sw: 0, abs: 0, skin: 0, t2m: 0, wet: 0, melting: 0, bare: 0, snowMass: 0, snowAlb: 0 };
    for (let i = 0; i < C; i++) {
      if (land[i] || !inZone(i)) continue;
      const a = mesh.areaCell[i], c = concentration[i] > 0 ? concentration[i] : 0, sw = shortwave[i];
      s.a += a; s.sw += a * sw; s.abs += a * (absorbed[i] - atmosphere[i]) / steps;
      if (!(c > 0)) continue;
      const ac = a * c, T = surfaceT[i];
      s.ac += ac; s.h += ac * ice[i]; s.snow += ac * snow[i]; s.skin += ac * T; s.t2m += ac * (t2m ? t2m[i] : NaN);
      if (T >= MELTING_POINT - 2) s.wet += ac;
      if (T >= MELTING_POINT - 0.1) s.melting += ac;
      if (snow[i] < 1) s.bare += ac;
      if (snow[i] > 0) { s.snowMass += ac * snow[i]; s.snowAlb += ac * snow[i] * snowAlbedo[i]; }
    }
    line[name] = {
      area: round(s.a / 1e12, 3), conc: round(s.ac / s.a, 3), h: round(s.h / s.ac, 3), snow: round(s.snow / s.ac, 2), snowAlbedo: round(s.snowAlb / s.snowMass, 3),
      skin: round(s.skin / s.ac - MELTING_POINT, 2), t2m: round(s.t2m / s.ac - MELTING_POINT, 2), wet: round(s.wet / s.ac, 3), melting: round(s.melting / s.ac, 3), bare: round(s.bare / s.ac, 3),
      swdn: round(s.sw / s.a, 1), absorbed: round(s.abs / s.a, 1), albedo: round(1 - s.abs / s.sw, 3),
    };
  }
  const swcre = {};
  for (const [from, to] of BANDS) {
    let a = 0, x = 0;
    for (let i = 0; i < C; i++) if (lat(i) >= from && lat(i) < to) { a += mesh.areaCell[i]; x += mesh.areaCell[i] * (absorbed[i] - clear[i]) / steps; }
    swcre[`${from}..${to}`] = round(x / a, 2);
  }
  let ca = 0, cw = 0, ia = 0, isw = 0;
  for (let i = 0; i < C; i++) {
    if (lat(i) >= 70) { ca += mesh.areaCell[i]; cw += mesh.areaCell[i] * cloud[i]; }
    if (lat(i) > 50 && !land[i] && concentration[i] >= 0.15) { ia += mesh.areaCell[i]; isw += mesh.areaCell[i] * shortwave[i]; }
  }
  line.swcre = swcre;
  line.cloudWater70 = round(cw / ca, 4);
  line.iceSwdn = round(isw / ia, 1);
  return line;
}
