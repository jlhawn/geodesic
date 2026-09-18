import { cellVector } from '../dynamics/operators.module.js';

const REFERENCE_SURFACE_T = 305.086;
const AVERAGE_SURFACE_T = 288;
const EQUATOR_POLE_CONTRAST = 45;

/*
 * Radiative–convective equilibrium θ(σ) of a single column over a surface
 * held at surfaceT: the model's own gray radiation and sensible heat
 * flux, with convective adjustment after every step, integrated until
 * the column stops changing. Computed on column 0 of the given model.
 */
export function equilibriumProfile(model, { surfaceT = REFERENCE_SURFACE_T, days = 600, dt = 900, windSpeed = 3, p0 = 101325 } = {}) {
  const { core, radiation, surface } = model;
  const { K, C, dSigma, cp, g, sigmaMid, exnerLayer } = core.diagnostics;
  const pi = new Float64Array(C).fill(p0);
  const theta = new Float64Array(K * C);
  core.diagnoseColumn(0, pi, theta);
  for (let k = 0; k < K; k++) theta[k * C] = surfaceT * Math.pow(sigmaMid[k], 0.19) / exnerLayer[k * C];
  const steps = Math.round(days * 86400 / dt);
  for (let n = 0; n < steps; n++) {
    core.diagnoseColumn(0, pi, theta);
    radiation.column(0, p0, theta, surfaceT, windSpeed);
    for (let k = 0; k < K; k++) {
      const massPerArea = p0 * dSigma[k] / g;
      theta[k * C] += dt * radiation.layerFlux[k] / (cp * massPerArea) / exnerLayer[k * C];
    }
    surface.convectiveAdjustColumn(0, pi, theta);
  }
  return Float64Array.from({ length: K }, (_, k) => theta[k * C]);
}

export function surfaceTemperature(lat) {
  const s = Math.sin(lat);
  return AVERAGE_SURFACE_T + EQUATOR_POLE_CONTRAST * (1 / 3 - s * s);
}

/*
 * Geopotential height of a pressure level in column i, interpolated in
 * ln p between the layer midpoints that bracket it. Needs the column's
 * diagnostics to be current.
 */
function geopotentialHeightAt(core, i, pi, pressure) {
  const { K, C, sigmaMid, geopotential, g } = core.diagnostics;
  let upper = 0;
  while (upper < K - 2 && pi * sigmaMid[upper + 1] < pressure) upper++;
  const lower = upper + 1;
  const pUpper = pi * sigmaMid[upper], pLower = pi * sigmaMid[lower];
  const t = (Math.log(pressure) - Math.log(pUpper)) / (Math.log(pLower) - Math.log(pUpper));
  return (geopotential[upper * C + i] + t * (geopotential[lower * C + i] - geopotential[upper * C + i])) / g;
}

/*
 * Initial state: the equilibrium profile shifted by each column's
 * surface-temperature offset (tapered by σ), a
 * wavenumber-5 θ seed at ±45°, surface pressure set by bisection so the
 * 500 hPa surface is level, optionally the hand-drawn subtropical-high /
 * subpolar-low pressure bands (off: they are not balanced by anything and
 * ring at 10 hPa), and winds in geostrophic balance with the model's own
 * pressure gradient force, tapered to zero inside ±15°.
 */
export function initializeState(model, {
  p0 = 101325, seedAmplitude = 2, seedWavenumber = 5, seedLatitude = 45, seedWidth = 15,
  bands = false, geostrophic = true, referencePressure = 50000, taperLatitude = 15, profile = null,
} = {}) {
  const { mesh, core } = model;
  const { K, C, E, sigmaMid, cp, exnerLayer, dExnerDpi, geopotential } = core.diagnostics;
  const reference = profile ?? equilibriumProfile(model, { p0 });
  const { nCells, latCell, lonCell, areaCell, cellsOnEdge, dcEdge, nEdge, xCell, fCell } = mesh;
  const deg = Math.PI / 180;
  const pi = new Float64Array(C).fill(p0);
  const theta = new Float64Array(K * C);
  const u = new Float64Array(K * E);
  const surfaceT = Float64Array.from(latCell, surfaceTemperature);

  for (let i = 0; i < C; i++) {
    const offset = surfaceT[i] - REFERENCE_SURFACE_T;
    const lat = latCell[i] / deg;
    const north = (lat - seedLatitude) / seedWidth, south = (lat + seedLatitude) / seedWidth;
    const envelope = Math.exp(-north * north) + Math.exp(-south * south);
    const wave = Math.cos(seedWavenumber * lonCell[i]);
    for (let k = 0; k < K; k++) {
      theta[k * C + i] = reference[k] + offset * sigmaMid[k] + seedAmplitude * envelope * wave * sigmaMid[k];
    }
  }

  let warmest = 0;
  for (let i = 1; i < C; i++) if (surfaceT[i] > surfaceT[warmest]) warmest = i;
  core.diagnoseColumn(warmest, pi, theta);
  const targetHeight = geopotentialHeightAt(core, warmest, pi[warmest], referencePressure);
  for (let i = 0; i < C; i++) {
    let lo = 0.7 * p0, hi = 1.3 * p0;
    for (let iteration = 0; iteration < 48; iteration++) {
      pi[i] = 0.5 * (lo + hi);
      core.diagnoseColumn(i, pi, theta);
      if (geopotentialHeightAt(core, i, pi[i], referencePressure) > targetHeight) hi = pi[i]; else lo = pi[i];
    }
  }
  let area = 0, piSum = 0;
  for (let i = 0; i < C; i++) { area += areaCell[i]; piSum += areaCell[i] * pi[i]; }
  const scale = p0 / (piSum / area);
  for (let i = 0; i < C; i++) pi[i] *= scale;

  if (bands) {
    const shape = (lat, center, width) => Math.exp(-(((lat - center) / width) ** 2));
    const taper = (lat) => {
      const a = Math.abs(lat);
      if (a <= 68) return 1;
      if (a >= 80) return 0;
      return Math.cos((a - 68) / 12 * Math.PI / 2) ** 2;
    };
    const added = new Float64Array(C);
    let addedSum = 0;
    for (let i = 0; i < C; i++) {
      const lat = latCell[i] / deg;
      const wobble = 5 * Math.sin(4 * lonCell[i]);
      let dp = 600 * (shape(lat, 30 + wobble, 12) + shape(lat, -30 - wobble, 12));
      dp += -2200 * (shape(lat, 60 - wobble, 12) + shape(lat, -60 + wobble, 12));
      added[i] = dp * taper(lat);
      addedSum += added[i] * areaCell[i];
    }
    const mean = addedSum / area;
    for (let i = 0; i < C; i++) pi[i] += added[i] - mean;
  }

  if (geostrophic) {
    core.diagnose(pi, theta);
    const force = new Float64Array(E);
    const vector = new Float64Array(3 * C);
    const wind = new Float64Array(3 * C);
    const sinTaper = Math.sin(taperLatitude * deg) ** 2;
    for (let k = 0; k < K; k++) {
      const off = k * C;
      for (let e = 0; e < E; e++) {
        const i = cellsOnEdge[2 * e], j = cellsOnEdge[2 * e + 1];
        const gradPhi = (geopotential[off + j] - geopotential[off + i]) / dcEdge[e];
        const pgfPi = cp * 0.5 * (theta[off + i] * dExnerDpi[off + i] + theta[off + j] * dExnerDpi[off + j]) * (pi[j] - pi[i]) / dcEdge[e];
        force[e] = gradPhi + pgfPi;
      }
      cellVector(mesh, force, vector);
      for (let i = 0; i < C; i++) {
        const s = Math.sin(latCell[i]);
        const taper = Math.min(1, s * s / sinTaper);
        const f = fCell[i];
        if (f === 0 || taper === 0) { wind[3 * i] = wind[3 * i + 1] = wind[3 * i + 2] = 0; continue; }
        const kx = xCell[3 * i], ky = xCell[3 * i + 1], kz = xCell[3 * i + 2];
        const gx = vector[3 * i], gy = vector[3 * i + 1], gz = vector[3 * i + 2];
        wind[3 * i] = taper * (ky * gz - kz * gy) / f;
        wind[3 * i + 1] = taper * (kz * gx - kx * gz) / f;
        wind[3 * i + 2] = taper * (kx * gy - ky * gx) / f;
      }
      for (let e = 0; e < E; e++) {
        const i = cellsOnEdge[2 * e], j = cellsOnEdge[2 * e + 1];
        u[k * E + e] = 0.5 * ((wind[3 * i] + wind[3 * j]) * nEdge[3 * e] + (wind[3 * i + 1] + wind[3 * j + 1]) * nEdge[3 * e + 1] + (wind[3 * i + 2] + wind[3 * j + 2]) * nEdge[3 * e + 2]);
      }
    }
  }

  return [pi, theta, u, surfaceT];
}
