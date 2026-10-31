import { encodeAerosol, decodeAerosol, AEROSOL_FIELDS } from '../../js/aerosol.module.js';

/*
 * A synthetic aerosol climatology in the file's own form, on a grid of nLon x nLat cells, with a month
 * center every 365 / nMonth days: optical depths 0.05 to 0.9 that fall from the equator to the poles and
 * move around the year, 1 - omega and asymmetries within the stored ranges, mean heights of 0.6 to 4.5 km,
 * and no aerosol at the south pole. `sample` overrides it: (field, lat, lon, month) in degrees -> value.
 */
export function syntheticAerosol({ nLon = 36, nLat = 18, nMonth = 12, sample = null } = {}) {
  const dLon = 360 / nLon, dLat = 180 / nLat, lon0 = -180 + 0.5 * dLon, lat0 = -90 + 0.5 * dLat, plane = nMonth * nLat * nLon, deg = Math.PI / 180;
  const centers = Array.from({ length: nMonth }, (_, m) => (m + 0.5) * 365 / nMonth);
  const fields = Array.from({ length: AEROSOL_FIELDS }, () => new Float64Array(plane));
  const natural = (f, lat, lon, month) => {
    const season = Math.sin(2 * Math.PI * (month + 0.5) / nMonth), north = Math.sin(lat * deg);
    const fine = 0.06 + 0.5 * Math.cos(lat * deg) ** 2 * (1 + 0.5 * Math.sin(lon * deg + 1)) * (1 + 0.4 * season * (north > 0 ? 1 : -1));
    const coarse = 0.03 + 0.4 * Math.exp(-(((lat - 20) / 25) ** 2)) * (1 + 0.5 * Math.cos(2 * lon * deg)) * (1 - 0.3 * season);
    switch (f) {
      case 0: return lat < -88 ? 0 : fine;
      case 1: return lat < -88 ? 0 : coarse;
      case 2: return 0.05 + 0.05 * Math.sin(lon * deg) ** 2 + 0.02 * season;
      case 3: return 0.015 + 0.025 * Math.cos(lat * deg) ** 2;
      case 4: return 0.6 + 0.06 * Math.cos(lon * deg);
      case 5: return 0.75 + 0.04 * Math.sin(lat * deg);
      case 6: return 1 + 3.5 * Math.cos(lat * deg) ** 2 * (0.6 + 0.4 * Math.sin(2 * lon * deg));
      default: return 0.6 + 2.4 * Math.exp(-(((lat - 15) / 30) ** 2)) * (0.7 + 0.3 * Math.cos(lon * deg));
    }
  };
  for (let f = 0; f < AEROSOL_FIELDS; f++) {
    for (let m = 0; m < nMonth; m++) for (let r = 0; r < nLat; r++) for (let c = 0; c < nLon; c++) {
      const lat = lat0 + r * dLat, lon = lon0 + c * dLon;
      fields[f][(m * nLat + r) * nLon + c] = sample ? sample(f, lat, lon, m) : natural(f, lat, lon, m);
    }
  }
  const bytes = encodeAerosol({ nLon, nLat, lon0, dLon, lat0, dLat, centers, fields, source: 'synthetic' });
  return { bytes, data: decodeAerosol(bytes) };
}

/*
 * One value in every cell: the eight fields [fine depth, coarse depth, fine 1 - omega, coarse 1 - omega, fine
 * asymmetry, coarse asymmetry, fine height, coarse height] of a uniform climatology.
 */
export const uniformAerosol = (values, options = {}) => syntheticAerosol({ ...options, sample: (f) => values[f] });

export function fillAerosol(fields, values) {
  const C = fields.length / AEROSOL_FIELDS;
  values.forEach((value, f) => fields.fill(value, f * C, (f + 1) * C));
  return fields;
}
