export const R = 287.06, CP = 1003.5, KAPPA = R / CP, G = 9.806, P0 = 1e5;
export const T0 = 288.15, LAPSE = 6.5e-3, TROPOPAUSE = 11000, T_STRATOSPHERE = 216.65;

export const exner = (p) => (p / P0) ** KAPPA;

export function pressureAt(z, surfaceT = T0, lapse = LAPSE) {
  if (lapse < 1e-6) return P0 * Math.exp(-G * z / (R * surfaceT));
  const zt = Math.min(z, TROPOPAUSE), tt = surfaceT - lapse * zt;
  const pt = P0 * (tt / surfaceT) ** (G / (R * lapse));
  return z <= TROPOPAUSE ? pt : pt * Math.exp(-G * (z - TROPOPAUSE) / (R * tt));
}

export function heightOf(p, surfaceT = T0, lapse = LAPSE) {
  if (lapse < 1e-6) return -R * surfaceT * Math.log(p / P0) / G;
  const tt = surfaceT - lapse * TROPOPAUSE, pt = P0 * (tt / surfaceT) ** (G / (R * lapse));
  if (p >= pt) return (surfaceT / lapse) * (1 - (p / P0) ** (R * lapse / G));
  return TROPOPAUSE - R * tt * Math.log(p / pt) / G;
}

export const temperatureAt = (z, surfaceT = T0, lapse = LAPSE) => Math.max(surfaceT - lapse * Math.min(z, TROPOPAUSE), z > TROPOPAUSE ? surfaceT - lapse * TROPOPAUSE : 0);

export const thetaAt = (z, surfaceT = T0, lapse = LAPSE) => temperatureAt(z, surfaceT, lapse) / exner(pressureAt(z, surfaceT, lapse));
