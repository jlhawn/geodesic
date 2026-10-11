/*
 * Columns of radiativeColumn.module.js from the equator to the pole, each
 * under its own annual-mean sunlight and ozone, none trading heat with
 * another. Hemisphere means weight them by area, each quantity taken as
 * linear in latitude between neighboring columns.
 */
import { makeColumn, meanOzone, K } from './radiativeColumn.module.js';
import { STEFAN_BOLTZMANN } from '../js/physics/radiation.module.js';

export const DATA_URL = new URL('./data/latitudes.json', import.meta.url).href;
export const STEP = 3 * 3600;
export const STEPS_PER_DAY = 86400 / STEP;
const BALANCED = 0.05;
export const CALM_DAYS = 10;
const D = Math.PI / 180;

export const heldSuarezSurface = (latitude) => 315 - 60 * Math.sin(latitude * D) ** 2;
export const effectiveTemperature = (flux) => (Math.max(0, flux) / STEFAN_BOLTZMANN) ** 0.25;

export function areaWeights(latitudes) {
  const weights = new Float64Array(latitudes.length);
  for (let n = 0; n + 1 < latitudes.length; n++) {
    const a = latitudes[n] * D, b = latitudes[n + 1] * D, chord = (Math.cos(a) - Math.cos(b)) / (b - a);
    weights[n] += chord - Math.sin(a);
    weights[n + 1] += Math.sin(b) - chord;
  }
  const total = Math.sin(latitudes[latitudes.length - 1] * D) - Math.sin(latitudes[0] * D);
  return weights.map((w) => w / total);
}

export const areaMean = (weights, values) => values.reduce((sum, v, n) => sum + weights[n] * v, 0);

export function temperatureAt({ pressure, T }, target) {
  for (let k = 1; k < K; k++) if (pressure[k] >= target) {
    const f = Math.log(target / pressure[k - 1]) / Math.log(pressure[k] / pressure[k - 1]);
    return T[k - 1] + f * (T[k] - T[k - 1]);
  }
  return T[K - 1];
}

export function readLatitudes(data) {
  return { ...data, temperature: data.temperature.map((T) => Float64Array.from(T)) };
}

let loading = null;
export function loadLatitudes(url = DATA_URL) {
  loading ??= fetch(url).then((response) => {
    if (!response.ok) throw new Error(`${url}: ${response.status}`);
    return response.json();
  }).then(readLatitudes);
  return loading;
}

export function createLatitudes(data, settings = {}) {
  const { latitudes, insolation, mu } = data, count = latitudes.length, weights = areaWeights(latitudes);
  const columns = latitudes.map((latitude, n) => makeColumn({
    beam: insolation[n], mu: mu[n], ozoneProfile: meanOzone(latitude),
    temperature: data.temperature?.[n] ?? 288, groundTemperature: data.groundTemperature?.[n] ?? 288, ...settings,
  }));
  const calm = new Int32Array(count), steps = new Int32Array(count), needed = CALM_DAYS * STEPS_PER_DAY;

  function changed() { calm.fill(0); steps.fill(0); }

  return {
    latitudes, columns, weights,
    set(changes) { for (const column of columns) column.set(changes); changed(); },
    restart(temperature = 288) { for (const column of columns) column.set({ temperature, groundTemperature: temperature }); changed(); },
    load(states, changes = {}) { columns.forEach((column, n) => column.set({ ...changes, temperature: states.temperature[n], groundTemperature: states.groundTemperature[n] })); changed(); },
    advance(n) {
      let active = 0;
      columns.forEach((column, i) => {
        for (let s = 0; s < n && calm[i] < needed; s++) {
          column.step(STEP);
          steps[i]++;
          const { imbalance, groundNet } = column.budget;
          calm[i] = Math.abs(imbalance) < BALANCED && Math.abs(groundNet) < BALANCED ? calm[i] + 1 : 0;
        }
        if (calm[i] < needed) active++;
      });
      return active;
    },
    get settled() { return calm.every((c) => c >= needed); },
    get days() { return Math.max(...steps) / STEPS_PER_DAY; },
    summary() {
      const rows = columns.map((column, n) => ({
        latitude: latitudes[n], Ts: column.Ts, T500: temperatureAt(column, 50000), target: heldSuarezSurface(latitudes[n]),
        insolation: column.budget.insolation, absorbed: column.budget.absorbed, olr: column.budget.olr, reflected: column.budget.reflected,
        albedo: column.budget.albedo, imbalance: column.budget.imbalance, coldPoint: column.coldPoint(),
      }));
      const mean = (key) => areaMean(weights, rows.map((r) => r[key]));
      const absorbed = mean('absorbed'), insolationMean = mean('insolation');
      return {
        rows,
        mean: { Ts: mean('Ts'), T500: mean('T500'), absorbed, olr: mean('olr'), insolation: insolationMean, reflected: mean('reflected'), albedo: mean('reflected') / insolationMean, Te: effectiveTemperature(absorbed) },
        contrast: rows[0].Ts - rows[count - 1].Ts,
        targetContrast: rows[0].target - rows[count - 1].target,
      };
    },
  };
}
