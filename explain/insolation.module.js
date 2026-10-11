/*
 * Daily-mean sunlight at the top of the air from the model's own sun
 * (sunDirection: a circular 365-day orbit, no refraction). Day d is the 24
 * hours centered on noon at longitude 0 on model day d, so day 0 is the
 * March equinox; each mean samples the sun perDay times at longitude 0.
 */
import { sunDirection, SOLAR_CONSTANT, DAY } from '../js/physics/radiation.module.js';

export const DAYS = 365;
export const MARKERS = [[0, 'March equinox'], [91, 'June solstice'], [183, 'September equinox'], [274, 'December solstice'], [365, 'March equinox']];
const D = Math.PI / 180;

export function sunTrack(day, perDay, out = { x: new Float64Array(perDay), z: new Float64Array(perDay) }) {
  const sun = new Float64Array(3);
  for (let k = 0; k < perDay; k++) {
    sunDirection((day - 0.5 + (k + 0.5) / perDay) * DAY, sun);
    out.x[k] = sun[0]; out.z[k] = sun[2];
  }
  return out;
}

export function insolationMap({ rows = 181, perDay = 288 } = {}) {
  const latitude = Float64Array.from({ length: rows }, (_, i) => -90 + 180 * i / (rows - 1));
  const daily = new Float32Array(rows * DAYS), annual = new Float64Array(rows), x = new Float64Array(DAYS * perDay), z = new Float64Array(DAYS * perDay);
  for (let d = 0; d < DAYS; d++) sunTrack(d, perDay, { x: x.subarray(d * perDay, (d + 1) * perDay), z: z.subarray(d * perDay, (d + 1) * perDay) });
  for (let i = 0; i < rows; i++) {
    const c = Math.cos(latitude[i] * D), s = Math.sin(latitude[i] * D);
    let year = 0;
    for (let d = 0, k = 0; d < DAYS; d++) {
      let sum = 0;
      for (const end = k + perDay; k < end; k++) { const mu = c * x[k] + s * z[k]; if (mu > 0) sum += mu; }
      const mean = SOLAR_CONSTANT * sum / perDay;
      daily[i * DAYS + d] = mean;
      year += mean;
    }
    annual[i] = year / DAYS;
  }
  return { rows, latitude, daily, annual };
}

export function dailyMean(latitude, day, perDay = 288) {
  const { x, z } = sunTrack(day, perDay), c = Math.cos(latitude * D), s = Math.sin(latitude * D);
  let sum = 0;
  for (let k = 0; k < perDay; k++) { const mu = c * x[k] + s * z[k]; if (mu > 0) sum += mu; }
  return SOLAR_CONSTANT * sum / perDay;
}

export function daylightHours(latitude, day, perDay = 1440) {
  const { x, z } = sunTrack(day, perDay), c = Math.cos(latitude * D), s = Math.sin(latitude * D);
  let lit = 0;
  for (let k = 0; k < perDay; k++) if (c * x[k] + s * z[k] > 0) lit++;
  return 24 * lit / perDay;
}

export const subsolarLatitude = (day) => Math.asin(sunDirection(day * DAY)[2]) / D;

export const noonElevation = (latitude, day) => 90 - Math.abs(latitude - subsolarLatitude(day));

export function seasonPhrase(day) {
  const d = ((Math.round(day) % DAYS) + DAYS) % DAYS;
  const [at, name] = MARKERS.reduce((best, m) => (Math.abs(d - m[0]) < Math.abs(d - best[0]) ? m : best));
  const n = Math.abs(d - at);
  return n === 0 ? `the ${name}` : `${n} day${n === 1 ? '' : 's'} ${d > at ? 'after' : 'before'} the ${name}`;
}
