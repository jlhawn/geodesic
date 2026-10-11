/*
 * The model's sun on the model's cells. Day d at hour h, solar time at
 * Greenwich, is model time t = d·DAY + (h − 12)·3600, so day 0 at noon is
 * t = 0: the March equinox with the sun over longitude 0.
 */
import { SOLAR_CONSTANT, DAY } from '../js/physics/radiation.module.js';

export { SOLAR_CONSTANT };

export const timeOf = (day, hour) => day * DAY + (hour - 12) * 3600;

export const subsolarPoint = (sun) => ({ lat: Math.asin(Math.max(-1, Math.min(1, sun[2]))), lon: Math.atan2(sun[1], sun[0]) });

export const sunlightAt = (x, y, z, sun) => SOLAR_CONSTANT * Math.max(0, x * sun[0] + y * sun[1] + z * sun[2]);

export function sunlightOnCells(xCell, sun, out = new Float64Array(xCell.length / 3)) {
  for (let i = 0; i < out.length; i++) out[i] = sunlightAt(xCell[3 * i], xCell[3 * i + 1], xCell[3 * i + 2], sun);
  return out;
}

export function areaMean(values, areaCell) {
  let sum = 0, area = 0;
  for (let i = 0; i < values.length; i++) { sum += values[i] * areaCell[i]; area += areaCell[i]; }
  return sum / area;
}
