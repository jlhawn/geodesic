import { cellVector } from '../dynamics/operators.module.js';

/*
 * Non-orographic gravity-wave drag after Alexander & Dunkerton (1999), as
 * GFDL's cg_drag runs it: each column launches, in each of two horizontal
 * directions (east and north), a spectrum of waves of one horizontal
 * wavelength at ground-relative phase speeds c = u₀ ± j Δc, j = 1 …
 * maxSpeed/Δc, where u₀ is the source layer's wind along the direction,
 * carrying the upward flux of momentum B(c) = sgn(c − u₀) B_m
 * exp(−ln 2 ((c − u₀)/halfWidth)²), B_m such that the absolute fluxes of
 * each direction sum to the column's flux F(φ) (Pa). F(φ) is Garfinkel et
 * al.'s (2022) eq. A3: `equatorialFlux` within `width` degrees of the
 * equator, rising linearly to `flux` at `edge` degrees, and poleward of
 * `edge` `flux` + ½ `northFlux` (1 + tanh((φ − edge)/width)) +
 * ½ `southFlux` (1 + tanh((φ + edge)/−width)). The source layer is the one
 * nearest σ = (`sourcePressure`/p₀)^cos φ (with `sourceDescent`; the level
 * index (K + 1) − (K + 1 − k₀) cos φ of cg_drag read as log-pressure
 * height), never the lowest. The spectrum is antisymmetric about u₀, so
 * a column launches no net momentum. Each wave rises unchanged until the
 * first layer below the lid layers where it meets its critical level
 * (c − u changes sign) or breaks, and leaves all of its momentum there. It breaks where its flux
 * exceeds the saturation flux ρ k |c − u|³ / (2N); with
 * `breakingAmplitude` (B_w, m²/s²) the flux tested is that of the wave
 * where it is present, ρ₀ B_w exp(−ln 2 ((c − u₀)/halfWidth)²) with ρ₀ the
 * source layer's density, of which the grid-box mean flux is the
 * intermittent fraction (AD99's ε); with `breakingAmplitude` null it is
 * the grid-box mean flux of the wave itself. What rises into the lid
 * layers, those whose midpoints lie above `lidPressure` (at least the top
 * layer), is spread over them at one acceleration (as Garfinkel et al.
 * 2022 spread what leaves the top), so each column's momentum is
 * conserved (Shepherd & Shaw 2004); with `lidTests` the waves are tested
 * in the lid layers too and only what rises out of the top layer is
 * spread, as MiMA's cg_drag does. The layer's
 * acceleration is g times its deposited flux over its mass per area. N²
 * is (g/θ) ∂θ/∂z between the neighbouring layers, z from the hypsometric
 * equation, floored at `minimumFrequency`². `compute` fills the cells'
 * east and north accelerations (`east`, `north`, K·C, m/s²) from the
 * state; `applyEdges` adds dt times their mean over an edge's two cells
 * to its normal velocity and, given `dissipation`, adds the kinetic
 * energy each edge layer loses (twice, per unit mass, as the ∇⁴ closure
 * does) for the model to return as heat. The cell on a pole, whose east
 * is undefined, launches nothing. With `diagnose`, `absoluteFlux` (K·C,
 * Pa) holds the absolute flux of the waves of both directions that rise
 * through each layer's upper interface.
 */
export const GRAVITY_WAVES = {
  flux: 4.3e-3, equatorialFlux: 4.3e-3, northFlux: 0, southFlux: 0, edge: 15, width: 10, sourcePressure: 31500, sourceDescent: true,
  halfWidth: 40, maxSpeed: 100, speedStep: 2, wavelength: 300e3, minimumFrequency: 0.005, breakingAmplitude: 0.4, lidPressure: 85, lidTests: false,
};

export function gravityWaveSpectrum({ flux = 1, halfWidth, maxSpeed, speedStep }) {
  const J = Math.floor(maxSpeed / speedStep), shape = Float64Array.from({ length: J }, (_, j) => Math.exp(-Math.LN2 * (((j + 1) * speedStep) / halfWidth) ** 2));
  const total = 2 * shape.reduce((s, x) => s + x, 0);
  return shape.map((x) => (flux * x) / total);
}

/*
 * The amplitude each wave's breaking test takes, times the column's
 * `gravityWaveScale`: B_w times the spectrum's shape with
 * `breakingAmplitude`, the unit spectrum without it.
 */
export function gravityWaveBreaking({ halfWidth, maxSpeed, speedStep, breakingAmplitude }) {
  if (!breakingAmplitude) return gravityWaveSpectrum({ halfWidth, maxSpeed, speedStep });
  return Float64Array.from({ length: Math.floor(maxSpeed / speedStep) }, (_, j) => breakingAmplitude * Math.exp(-Math.LN2 * (((j + 1) * speedStep) / halfWidth) ** 2));
}

/*
 * The spectrum's partial sums, cumulative[n] = Σ_{j<n} amplitude[j]:
 * the waves of one side leave a layer as a contiguous run of phase
 * speeds nearest the source wind, since both the critical level and the
 * saturation flux take the slowest waves first.
 */
export function gravityWaveSums(amplitude) {
  const sums = new Float64Array(amplitude.length + 1);
  for (let j = 0; j < amplitude.length; j++) sums[j + 1] = sums[j] + amplitude[j];
  return sums;
}

export function gravityWaveFlux(latitude, { flux, equatorialFlux, northFlux, southFlux, edge, width }) {
  const phi = latitude * 180 / Math.PI;
  if (Math.abs(phi) <= width) return equatorialFlux;
  if (Math.abs(phi) <= edge) return flux + (equatorialFlux - flux) * (edge - Math.abs(phi)) / (edge - width);
  return flux + 0.5 * northFlux * (1 + Math.tanh((phi - edge) / width)) + 0.5 * southFlux * (1 + Math.tanh((phi + edge) / -width));
}

export function gravityWaveSource(sigmaMid, sourcePressure, p0, latitude = 0, sourceDescent = false) {
  const target = (sourcePressure / p0) ** (sourceDescent ? Math.cos(latitude) : 1);
  let best = 0;
  for (let k = 1; k < sigmaMid.length - 1; k++) if (Math.abs(sigmaMid[k] - target) < Math.abs(sigmaMid[best] - target)) best = k;
  return best;
}

export function gravityWaveLid(sigmaMid, lidPressure, p0) {
  let layers = 1;
  while (layers < sigmaMid.length && sigmaMid[layers] * p0 < lidPressure) layers++;
  return layers;
}

/*
 * Per cell: the source layer, the flux F(φ) and whether it launches (not
 * on a pole); `deepest` is the lowest source layer.
 */
export function gravityWaveColumns(mesh, sigmaMid, p0, options) {
  const o = { ...GRAVITY_WAVES, ...options }, C = mesh.nCells;
  const source = new Int32Array(C), scale = new Float64Array(C);
  let deepest = 0;
  for (let i = 0; i < C; i++) {
    const lat = mesh.latCell[i];
    source[i] = gravityWaveSource(sigmaMid, o.sourcePressure, p0, lat, o.sourceDescent);
    scale[i] = Math.abs(lat) > Math.PI / 2 - 1e-9 ? 0 : gravityWaveFlux(lat, o);
    deepest = Math.max(deepest, source[i]);
  }
  return { source, scale, deepest };
}

export function createGravityWaveDrag(mesh, core, options = {}) {
  const o = { ...GRAVITY_WAVES, ...options };
  const { halfWidth, maxSpeed, speedStep, wavelength, minimumFrequency, breakingAmplitude, lidPressure, lidTests, buffers = null, diagnose = false } = o;
  const { K, C, E, sigmaMid, dSigma, R, g, kappa, p0 } = core.diagnostics;
  const { cellsOnEdge, nEdge, latCell, lonCell } = mesh;
  const amplitude = gravityWaveSpectrum({ halfWidth, maxSpeed, speedStep }), J = amplitude.length, cumulative = gravityWaveSums(amplitude);
  const breaking = gravityWaveBreaking({ halfWidth, maxSpeed, speedStep, breakingAmplitude });
  const columns = gravityWaveColumns(mesh, sigmaMid, p0, o), source = columns.deepest, wavenumber = 2 * Math.PI / wavelength, floor = minimumFrequency * minimumFrequency;
  const lid = gravityWaveLid(sigmaMid, lidPressure, p0), lidShare = new Float64Array(lid), tested = lidTests ? 0 : lid;
  for (let k = 0, total = dSigma.subarray(0, lid).reduce((s, x) => s + x, 0); k < lid; k++) lidShare[k] = dSigma[k] / total;
  const shared = { east: buffers ? buffers.east : new SharedArrayBuffer(8 * K * C), north: buffers ? buffers.north : new SharedArrayBuffer(8 * K * C) };
  const east = new Float64Array(shared.east), north = new Float64Array(shared.north);
  const projection = new Float64Array(4 * E);
  for (let e = 0; e < E; e++) {
    for (let s = 0; s < 2; s++) {
      const i = cellsOnEdge[2 * e + s], lat = latCell[i], lon = lonCell[i], n = nEdge.subarray(3 * e, 3 * e + 3);
      if (columns.scale[i] === 0) continue;
      projection[4 * e + 2 * s] = -Math.sin(lon) * n[0] + Math.cos(lon) * n[1];
      projection[4 * e + 2 * s + 1] = -Math.sin(lat) * Math.cos(lon) * n[0] - Math.sin(lat) * Math.sin(lon) * n[1] + Math.cos(lat) * n[2];
    }
  }
  const vector = new Float64Array(3 * C), windEast = new Float64Array(K * C), windNorth = new Float64Array(K * C);
  const absoluteFlux = diagnose ? new Float64Array(K * C) : null, passing = new Float64Array(K);
  const density = new Float64Array(K), frequency = new Float64Array(K), temperature = new Float64Array(K), pressure = new Float64Array(K), deposit = new Float64Array(K);

  function launch(wind, i, out, columnPressure) {
    deposit.fill(0);
    if (absoluteFlux) passing.fill(0);
    const top = columns.source[i], scale = columns.scale[i];
    const u0 = wind[top * C + i], present = breakingAmplitude ? density[top] : scale;
    for (let side = -1; side <= 1; side += 2) {
      let gone = 0;
      for (let k = top - 1; k >= tested && gone < J; k--) {
        const ahead = side * (u0 - wind[k * C + i]), saturation = density[k] * wavenumber / (2 * frequency[k]);
        let reached = Math.max(gone, Math.min(J, Math.floor(-ahead / speedStep)));
        while (reached < J) {
          const relative = ahead + (reached + 1) * speedStep;
          if (present * breaking[reached] < saturation * relative * relative * relative) break;
          reached++;
        }
        const left = scale * (cumulative[reached] - cumulative[gone]);
        deposit[k] += side * left;
        if (absoluteFlux) passing[k] += left;
        gone = reached;
      }
      const escaping = scale * (cumulative[J] - cumulative[gone]);
      for (let k = 0; k < lid; k++) {
        deposit[k] += side * escaping * lidShare[k];
        if (absoluteFlux) passing[k] += escaping * lidShare[k];
      }
    }
    for (let k = 0; k < top; k++) out[k * C + i] = g * deposit[k] / (columnPressure * dSigma[k]);
    if (absoluteFlux) for (let k = 0, above = 0; k < top; k++) { absoluteFlux[k * C + i] += above; above += passing[k]; }
  }

  function windsOf(u, k, iFrom, iTo) {
    cellVector(mesh, u.subarray(k * E, (k + 1) * E), vector, iFrom, iTo);
    for (let i = iFrom; i < iTo; i++) {
      const lat = latCell[i], lon = lonCell[i], x = vector[3 * i], y = vector[3 * i + 1], z = vector[3 * i + 2];
      windEast[k * C + i] = -Math.sin(lon) * x + Math.cos(lon) * y;
      windNorth[k * C + i] = -Math.sin(lat) * Math.cos(lon) * x - Math.sin(lat) * Math.sin(lon) * y + Math.cos(lat) * z;
    }
  }

  function compute(state, iFrom = 0, iTo = C) {
    const [pi, theta, u] = state;
    for (let k = 0; k <= source; k++) windsOf(u, k, iFrom, iTo);
    for (let i = iFrom; i < iTo; i++) {
      for (let k = 0; k < K; k++) { east[k * C + i] = 0; north[k * C + i] = 0; if (absoluteFlux) absoluteFlux[k * C + i] = 0; }
      if (columns.scale[i] === 0) continue;
      const top = columns.source[i];
      for (let k = 0; k <= top; k++) {
        pressure[k] = sigmaMid[k] * pi[i];
        temperature[k] = theta[k * C + i] * Math.pow(pressure[k] / p0, kappa);
        density[k] = pressure[k] / (R * temperature[k]);
      }
      for (let k = 0; k < top; k++) {
        const above = Math.max(0, k - 1), below = k + 1;
        const depth = R / g * 0.5 * (temperature[above] + temperature[below]) * Math.log(pressure[below] / pressure[above]);
        const n2 = g / theta[k * C + i] * (theta[above * C + i] - theta[below * C + i]) / depth;
        frequency[k] = Math.sqrt(Math.max(n2, floor));
      }
      launch(windEast, i, east, pi[i]);
      launch(windNorth, i, north, pi[i]);
    }
  }

  function applyEdges(u, eFrom, eTo, dt, dissipation = null) {
    for (let e = eFrom; e < eTo; e++) {
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      for (let k = 0; k < source; k++) {
        const change = 0.5 * dt * (east[k * C + a] * projection[4 * e] + north[k * C + a] * projection[4 * e + 1] + east[k * C + b] * projection[4 * e + 2] + north[k * C + b] * projection[4 * e + 3]);
        if (change === 0) continue;
        const n = k * E + e, before = u[n];
        u[n] = before + change;
        if (dissipation) dissipation[n] += before * before - u[n] * u[n];
      }
    }
  }

  return { compute, applyEdges, east, north, absoluteFlux, source, columns, lid, amplitude, shared };
}
