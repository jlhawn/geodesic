/*
 * A resolved column of Price, Weller and Pinkel (1986) under a daily
 * cycle, the reference for the mixed layer's retreat: `dz` levels to
 * `depth` from a 20 °C layer 40 m deep over 0.05 K/m, a constant
 * linear expansion `alpha`, absorbed sunlight `sun` (W/m², a function of
 * the seconds from day 0) spread by depth over a Jerlov `type`, a
 * constant dark loss `loss` (W/m²) from the top level and a constant
 * `wind` (m/s, stress 1.2 · 1.2e-3 U²) on what the surface mixes, the
 * inertial current damped over two days. Each step mixes downward while
 * the column is statically unstable, then to a bulk Richardson number of
 * 0.65, then relaxes a gradient Richardson number under 0.25 over four
 * passes. The dawn depth is where the profile first lies 0.02 K below the
 * surface at 05:52, the median over the last ten days, and `range` the
 * median daily range of the top level.
 */
import { JERLOV } from '../../js/ocean/layered.module.js';

const GRAVITY = 9.81, RHO0 = 1025, RHO_CP = 1025 * 3985;

export function pwpColumn({ wind, sun, loss, latitude = 45, days = 30, dz = 0.5, depth = 200, dt = 900, type = 'IB', alpha = 2.5e-4 }) {
  const [R, z1, z2] = JERLOV[type], passing = (z) => R * Math.exp(-z / z1) + (1 - R) * Math.exp(-z / z2);
  const nz = Math.round(depth / dz), f = 2 * 7.292e-5 * Math.sin(latitude * Math.PI / 180);
  const t = new Float64Array(nz), u = new Float64Array(nz), v = new Float64Array(nz), share = new Float64Array(nz);
  for (let k = 0; k < nz; k++) {
    const z = (k + 0.5) * dz;
    t[k] = z < 40 ? 20 : 20 - 0.05 * (z - 40);
    share[k] = passing(k * dz) - passing((k + 1) * dz);
  }
  const tau = 1.2 * 1.2e-3 * wind * wind, damp = 1 / (2 * 86400), c = Math.cos(f * dt), s = Math.sin(f * dt);
  const mixTo = (n) => {
    let st = 0, su = 0, sv = 0;
    for (let k = 0; k <= n; k++) { st += t[k]; su += u[k]; sv += v[k]; }
    for (let k = 0; k <= n; k++) { t[k] = st / (n + 1); u[k] = su / (n + 1); v[k] = sv / (n + 1); }
  };
  const steps = Math.round(days * 86400 / dt), dawn = [], tops = new Map();
  for (let step = 0; step < steps; step++) {
    const middle = (step + 0.5) * dt, hour = (middle / 3600) % 24, light = sun(middle);
    for (let k = 0; k < nz; k++) t[k] += light * share[k] * dt / (RHO_CP * dz);
    t[0] -= loss * dt / (RHO_CP * dz);
    let n = 0;
    while (n + 1 < nz && t[n + 1] >= t[0] - 1e-12 / alpha) { n++; mixTo(n); }
    let mixed = (n + 1) * dz;
    for (let k = 0; k < nz; k++) {
      const uu = u[k], vv = v[k];
      u[k] = (uu * c + vv * s) * (1 - damp * dt); v[k] = (vv * c - uu * s) * (1 - damp * dt);
    }
    for (let k = 0; k <= n; k++) u[k] += tau / (RHO0 * mixed) * dt;
    while (n + 1 < nz) {
      const db = GRAVITY * alpha * (t[0] - t[n + 1]), du2 = (u[0] - u[n + 1]) ** 2 + (v[0] - v[n + 1]) ** 2;
      if (db * mixed / Math.max(du2, 1e-12) >= 0.65) break;
      n++; mixTo(n); mixed = (n + 1) * dz;
    }
    for (let pass = 0; pass < 4; pass++) for (let k = n; k + 1 < nz && k < n + 40; k++) {
      const db = GRAVITY * alpha * (t[k] - t[k + 1]), du2 = (u[k] - u[k + 1]) ** 2 + (v[k] - v[k + 1]) ** 2, rg = db * dz / Math.max(du2, 1e-12);
      if (rg < 0.25) {
        const r = (1 - rg / 0.25) / 2, dtk = (t[k] - t[k + 1]) * r, duk = (u[k] - u[k + 1]) * r, dvk = (v[k] - v[k + 1]) * r;
        t[k] -= dtk; t[k + 1] += dtk; u[k] -= duk; u[k + 1] += duk; v[k] -= dvk; v[k + 1] += dvk;
      }
    }
    const day = Math.floor(middle / 86400);
    if (day < days - 10) continue;
    if (Math.abs(hour - 5.875) < 0.2) {
      let k = 0;
      while (k + 1 < nz && t[0] - t[k + 1] < 0.02) k++;
      dawn.push((k + 1) * dz);
    }
    const [low, high] = tops.get(day) ?? [Infinity, -Infinity];
    tops.set(day, [Math.min(low, t[0]), Math.max(high, t[0])]);
  }
  const median = (a) => [...a].sort((x, y) => x - y)[Math.floor(a.length / 2)];
  return { dawn: median(dawn), range: median([...tops.values()].map(([low, high]) => high - low)) };
}
