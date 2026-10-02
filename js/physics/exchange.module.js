import { VIRTUAL_FACTOR } from '../dynamics/sigmaCore.module.js';
import { saturationHumidity } from './moist.module.js';
import { SEA_DRAG, LAND_DRAG } from './surface.module.js';

/*
 * The surface layer: per cell the transfer coefficients C_D (momentum)
 * and C_H (heat and vapour) between the surface and the lowest layer,
 * at that layer's own height z above the ground, from the roughness
 * lengths of what the cell carries and the stability of its bulk
 * Richardson number.
 *
 * Roughness (z0m for momentum, z0h for heat and vapour, m):
 *  - land: forest, grass and bare soil at the IFS's values (IFS Cy47r3
 *    Part IV, Table 8.3, calibrated so that the 10 m wind's error against
 *    SYNOP vanishes per vegetation type: evergreen and deciduous trees
 *    2.0 / 2.0, short grass 0.1 / 0.001, desert 0.013 / 0.00013). Snow
 *    covers grass and bare soil over the share min(1, S / snowCover) of
 *    their area (the IFS's c_sn with its 0.1 m depth at 300 kg/m³) and
 *    leaves the trees standing; snow and the ice sheets take the IFS's
 *    ice caps' z0m 1.3e-3 m. With the land's vegetation off the land is
 *    grass.
 *  - open sea: COARE 3.5 (Edson et al. 2013), z0m = α u*²/g + 0.11 ν/u*
 *    with Charnock's α = 0.0017 U10N − 0.005 (U10N at most 19 m/s, α at
 *    least 0) and ν the air's viscosity at its temperature, iterated
 *    from the neutral u*; z0h = min(1.6e-4, 5.8e-5 Rr^−0.72), Rr = z0m u* / ν.
 *  - sea ice of concentration A: the IFS's (Cy47r3 eq. 3.30, after
 *    Andreas et al. 2010) z0m = max(1e-3, 0.93e-3 (1 − A) + 6.05e-3
 *    exp(−17 (A − 0.5)²)).
 *  - the scalar roughness over snow, ice sheets and sea ice: Andreas
 *    (1987, Table 2 of Andreas 2002), ln(z0h/z0m) = b0 + b1 ln R* +
 *    b2 (ln R*)², R* = u* z0m / ν with the tile's neutral u*.
 * A cell's tiles are blended as the IFS aggregates its tiles' roughness:
 * the area-weighted neutral coefficients at blendingHeight (10 m) give
 * the cell's z0m and, with it, z0h.
 *
 * Stability: Monin–Obukhov with the IFS's surface-layer functions
 * (Cy47r3 eqs. 3.16–3.26): Dyer–Hicks integrated by Paulson when
 * unstable, Holtslag and De Bruin (1988) with a = 1, b = 2/3, c = 5,
 * d = 0.35 when stable, ζ = z/L found from
 * Ri_b = g z (θv − θv_s) / (θ̄v U²) by `iterations` steps (fixed point
 * when unstable, Newton in ln ζ when stable). U is the lowest wind with
 * the gustiness floor, θv_s the surface's virtual potential temperature,
 * with the saturation humidity at the skin over sea and sea ice and dry
 * over land.
 *
 * `reference` is the FAO-56 reference grass's neutral C_H at the same
 * height (Allen et al. 1998 eq. 4: crop height 0.12 m, d = 2/3 h,
 * z0m = 0.123 h, z0h = 0.1 z0m, κ 0.41), which the radiation's reference
 * evapotranspiration takes.
 *
 * exchange 'fixed' gives constant coefficients: the sea's dragCoefficient
 * (SEA_DRAG) and the land's (LAND_DRAG) for momentum, heat and vapour
 * alike. An option set that names either coefficient and no `exchange` is
 * fixed.
 */
export const KARMAN = 0.4;
export const FAO_KARMAN = 0.41;
export const ROUGHNESS = { forest: [2.0, 2.0], grass: [0.1, 1e-3], bare: [0.013, 1.3e-4], snow: 1.3e-3 };
export const EXCHANGE_DEFAULTS = { roughness: ROUGHNESS, snowCover: 30, blendingHeight: 10, charnock: [0.0017, -0.005, 19], smoothFlow: 0.11, iterations: 5, referenceCrop: 0.12 };
export const ANDREAS = [[0.135, 1.25, 0, 0], [2.5, 0.149, -0.55, 0], [Infinity, 0.317, -0.565, -0.183]];
const GRAVITY = 9.81;
const STABLE = { a: 1, b: 2 / 3, c: 5, d: 0.35 };

export function exchangeMode(surface = {}, land = {}) {
  const named = surface.dragCoefficient !== undefined || land.dragCoefficient !== undefined;
  const mode = surface.exchange ?? (named ? 'fixed' : 'roughness');
  if (mode !== 'roughness' && mode !== 'fixed') throw new Error(`exchange must be 'roughness' or 'fixed', not ${mode}`);
  if (mode === 'roughness' && named) throw new Error('a dragCoefficient belongs to exchange \'fixed\'');
  return mode;
}

export function exchangeOptions(surface = {}) {
  const o = { ...EXCHANGE_DEFAULTS, ...Object.fromEntries(Object.keys(EXCHANGE_DEFAULTS).filter((key) => surface[key] !== undefined).map((key) => [key, surface[key]])) };
  return { ...o, roughness: { ...ROUGHNESS, ...o.roughness } };
}

export const airViscosity = (celsius) => 1.326e-5 * (1 + 6.542e-3 * celsius + 8.301e-6 * celsius * celsius - 4.84e-9 * celsius * celsius * celsius);

export function psiMomentum(zeta) {
  if (zeta < 0) {
    const x = Math.pow(1 - 16 * zeta, 0.25);
    return Math.PI / 2 - 2 * Math.atan(x) + Math.log((1 + x) * (1 + x) * (1 + x * x) / 8);
  }
  const { a, b, c, d } = STABLE;
  return -b * (zeta - c / d) * Math.exp(-d * zeta) - a * zeta - b * c / d;
}

export function psiHeat(zeta) {
  if (zeta < 0) return 2 * Math.log((1 + Math.sqrt(1 - 16 * zeta)) / 2);
  const { a, b, c, d } = STABLE;
  return -b * (zeta - c / d) * Math.exp(-d * zeta) - Math.pow(1 + 2 * a * zeta / 3, 1.5) - b * c / d + 1;
}

const slopeMomentum = (zeta) => { const { a, b, c, d } = STABLE; return -b * Math.exp(-d * zeta) * (1 + c - d * zeta) - a; };
const slopeHeat = (zeta) => { const { a, b, c, d } = STABLE; return -b * Math.exp(-d * zeta) * (1 + c - d * zeta) - a * Math.sqrt(1 + 2 * a * zeta / 3); };

/*
 * The denominators of C_D = κ²/F_m² and C_H = κ²/(F_m F_h) at ζ = z/L.
 */
function profiles(zeta, z, z0m, z0h, out) {
  const top = zeta * (1 + z0m / z);
  out.m = Math.log((z + z0m) / z0m) - psiMomentum(top) + psiMomentum(zeta * z0m / z);
  out.h = Math.log((z + z0m) / z0h) - psiHeat(top) + psiHeat(zeta * z0h / z);
  return out;
}

const scratch = { m: 0, h: 0 };
export function transfer(ri, z, z0m, z0h, iterations = EXCHANGE_DEFAULTS.iterations, out = { drag: 0, heat: 0, zeta: 0 }) {
  let f = profiles(0, z, z0m, z0h, scratch);
  let zeta = ri * f.m * f.m / f.h;
  if (ri > 0) {
    const target = Math.log(ri);
    let s = Math.log(zeta);
    for (let n = 0; n < iterations; n++) {
      zeta = Math.exp(s);
      f = profiles(zeta, z, z0m, z0h, scratch);
      const top = zeta * (1 + z0m / z), lowM = zeta * z0m / z, lowH = zeta * z0h / z;
      const dm = -slopeMomentum(top) * top + slopeMomentum(lowM) * lowM, dh = -slopeHeat(top) * top + slopeHeat(lowH) * lowH;
      s -= (s - target - Math.log(f.m * f.m / f.h)) / (1 - 2 * dm / f.m + dh / f.h);
    }
    zeta = Math.exp(s);
  } else {
    for (let n = 0; n < iterations; n++) { f = profiles(zeta, z, z0m, z0h, scratch); zeta = ri * f.m * f.m / f.h; }
  }
  f = profiles(zeta, z, z0m, z0h, scratch);
  out.drag = KARMAN * KARMAN / (f.m * f.m);
  out.heat = KARMAN * KARMAN / (f.m * f.h);
  out.zeta = zeta;
  return out;
}

export function andreasScalar(z0, ustar, viscosity) {
  const reynolds = ustar * z0 / viscosity, l = Math.log(Math.min(1000, reynolds));
  const [, b0, b1, b2] = ANDREAS.find(([limit]) => reynolds <= limit) ?? ANDREAS[2];
  return z0 * Math.exp(b0 + b1 * l + b2 * l * l);
}

export const seaIceRoughness = (concentration) => Math.max(1e-3, 0.93e-3 * (1 - concentration) + 6.05e-3 * Math.exp(-17 * (concentration - 0.5) ** 2));

export function charnockRoughness(wind, z, viscosity, { charnock: [slope, offset, ceiling], smoothFlow } = EXCHANGE_DEFAULTS, iterations = 4, out = { momentum: 0, heat: 0 }) {
  let z0 = 1e-4, ustar = 0;
  for (let n = 0; n < iterations; n++) {
    ustar = KARMAN * wind / Math.log((z + z0) / z0);
    const wind10 = ustar / KARMAN * Math.log((10 + z0) / z0);
    const alpha = Math.max(0, slope * Math.min(wind10, ceiling) + offset);
    z0 = alpha * ustar * ustar / GRAVITY + smoothFlow * viscosity / ustar;
  }
  ustar = KARMAN * wind / Math.log((z + z0) / z0);
  out.momentum = z0;
  out.heat = Math.min(1.6e-4, 5.8e-5 / Math.pow(z0 * ustar / viscosity, 0.72));
  return out;
}

/*
 * Accumulates tiles' neutral coefficients at the blending height and
 * backs out the blend's roughness lengths.
 */
export function createBlend(height) {
  let drag = 0, heat = 0;
  const blend = {
    momentum: 0, heat: 0,
    clear() { drag = 0; heat = 0; return blend; },
    add(share, z0m, z0h) {
      if (!(share > 0)) return blend;
      const lm = Math.log(height / z0m);
      drag += share * KARMAN * KARMAN / (lm * lm);
      heat += share * KARMAN * KARMAN / (lm * Math.log(height / z0h));
      return blend;
    },
    finish() {
      const lm = KARMAN / Math.sqrt(drag);
      blend.momentum = height * Math.exp(-lm);
      blend.heat = height * Math.exp(-KARMAN * KARMAN / (heat * lm));
      return blend;
    },
  };
  return blend;
}

export function referenceCoefficient(z, crop = EXCHANGE_DEFAULTS.referenceCrop) {
  const d = 2 / 3 * crop, z0m = 0.123 * crop, z0h = 0.1 * z0m;
  return FAO_KARMAN * FAO_KARMAN / (Math.log((z - d) / z0m) * Math.log((z - d) / z0h));
}

export function createSurfaceExchange(mesh, core, { geography = null, vegetated = true, mode = 'roughness', seaDrag = SEA_DRAG, landDrag = LAND_DRAG, gustiness = 3, buffers = null, ...options } = {}) {
  const C = mesh.nCells;
  const { K, cp, g, exnerLayer, exnerLower } = core.diagnostics;
  const thetaV = core.arrays.thetaV;
  const landMask = geography ? geography.land : null, iceSheet = geography ? geography.iceSheet : null;
  const shared = (name) => (buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * C));
  const dragBuffer = shared('drag');
  const drag = new Float64Array(dragBuffer);
  const fixed = mode === 'fixed';
  const heatBuffer = fixed ? dragBuffer : shared('heat'), referenceBuffer = fixed ? dragBuffer : shared('reference');
  const heat = fixed ? drag : new Float64Array(heatBuffer), reference = fixed ? drag : new Float64Array(referenceBuffer);
  const o = exchangeOptions(options);
  const { roughness: { forest, grass, bare, snow: snowRoughness }, snowCover, iterations } = o;
  if (!(buffers && buffers.drag)) {
    for (let i = 0; i < C; i++) {
      drag[i] = landMask && landMask[i] ? landDrag : seaDrag;
      if (!fixed) { heat[i] = drag[i]; reference[i] = drag[i]; }
    }
  }
  const blend = createBlend(o.blendingHeight), sea = { momentum: 0, heat: 0 }, coefficients = { drag: 0, heat: 0, zeta: 0 }, last = { z: 0, ri: 0, momentum: 0, heat: 0 };
  const bottom = K - 1;

  /*
   * Cell i from the state, its skin temperature `skin`, its lowest wind
   * `wind`, its sea-ice concentration and, on land, its snow (kg/m²),
   * cover and trees; returns the height, bulk Richardson number and
   * roughness lengths it used.
   */
  function cell(i, pi, theta, q, qc, skin, wind, concentration = 0, snow = 0, cover = 0, trees = 0) {
    const b = bottom * C + i;
    const z = cp * thetaV[b] * (exnerLower[b] - exnerLayer[b]) / g;
    const speed = Math.max(wind, gustiness);
    const viscosity = airViscosity(theta[b] * exnerLayer[b] - 273.15);
    const onLand = landMask && landMask[i];
    blend.clear();
    const snowTile = (share) => {
      const ustar = KARMAN * speed / Math.log((z + snowRoughness) / snowRoughness);
      blend.add(share, snowRoughness, andreasScalar(snowRoughness, ustar, viscosity));
    };
    if (onLand && iceSheet && iceSheet[i]) snowTile(1);
    else if (onLand) {
      const t = vegetated ? Math.min(1, trees) : 0, gr = vegetated ? Math.max(0, cover - t) : 1, bareShare = Math.max(0, 1 - t - gr);
      const covered = Math.min(1, snow / snowCover);
      blend.add(t, forest[0], forest[1]).add(gr * (1 - covered), grass[0], grass[1]).add(bareShare * (1 - covered), bare[0], bare[1]);
      snowTile((gr + bareShare) * covered);
    } else {
      charnockRoughness(speed, z, viscosity, o, 4, sea);
      blend.add(1 - concentration, sea.momentum, sea.heat);
      if (concentration > 0) {
        const z0 = seaIceRoughness(concentration), ustar = KARMAN * speed / Math.log((z + z0) / z0);
        blend.add(concentration, z0, andreasScalar(z0, ustar, viscosity));
      }
    }
    blend.finish();
    const exS = exnerLower[b];
    const airV = theta[b] * (1 + VIRTUAL_FACTOR * (q ? q[b] : 0) - (qc ? qc[b] : 0));
    const surfaceQ = onLand ? (q ? q[b] : 0) : (q ? saturationHumidity(skin, pi[i]) : 0);
    const surfaceV = skin / exS * (1 + VIRTUAL_FACTOR * surfaceQ);
    const ri = g * z * (airV - surfaceV) / (0.5 * (airV + surfaceV) * speed * speed);
    transfer(ri, z, blend.momentum, blend.heat, iterations, coefficients);
    drag[i] = coefficients.drag;
    heat[i] = coefficients.heat;
    reference[i] = referenceCoefficient(z, o.referenceCrop);
    last.z = z; last.ri = ri; last.momentum = blend.momentum; last.heat = blend.heat;
    return last;
  }

  return { mode, fixed, drag, heat, reference, cell, options: o, shared: { drag: dragBuffer, heat: heatBuffer, reference: referenceBuffer } };
}
