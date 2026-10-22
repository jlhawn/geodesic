import { SEAWATER, seawaterDensity, thermalExpansion, halineContraction } from './seawater.module.js';

/*
 * The column library of the ALE ocean (docs/ocean-hybrid-coordinate.md):
 * pure functions on one column of K layers, each with a thickness h_k, a
 * temperature and a salinity, interface k being the top of layer k and
 * depth measured downward from the free surface. The ocean's step moves
 * the layers with the flow, then regrids each column to a target made by
 * one of three generators and remaps the water onto it, then runs the
 * column physics on the new layers.
 *
 * Layers (hybridLayerTable): label densities ρ_k and a minimum
 * thickness Δz_k per layer whose running sum is
 * the depth table z_k, the z grid a column takes where every layer is
 * held at its minimum. Twelve surface layers have labels lighter than any
 * seawater, so they are fixed-depth everywhere, 3 m at the top growing
 * 8 % a layer to 57 m; the nine warm-pool labels follow at 8 to 10 m,
 * to 144 m; the 23 thermocline and intermediate labels keep 10 m each,
 * less than their isopycnal spacing wherever the thermocline is
 * stratified, so there they are isopycnal; and the 14 deep labels grow
 * from 40 to 500 m so that a weakly stratified column is a z grid to
 * about 3 km.
 *
 * Targets (targetInterfaces): the isopycnal position d_k of interface k is
 * where the column's potential density, its layer means held from
 * decreasing with depth and reconstructed by the remap's limited
 * parabolas, crosses the boundary half-way between labels k−1 and k. A
 * layer whose two interfaces are both isopycnal is a uniform class, so
 * its reconstruction is flattened by the smaller of their mode weights,
 * given as `previousWeight` or read from the layers as they stand (the
 * excess of each layer over its minimum, or 1 for the isopycnal
 * generator), and a boundary that falls in the
 * jump at an interface puts d_k at that interface: a column of uniform
 * layers at their labels is already on target. Water lighter than the
 * surface value puts the interface at the surface and water denser than
 * the bottom value at the floor. The
 * isopycnal generator takes d_k; the z* generator the depth table
 * stretched to the column's height; the hybrid generator the deeper of
 * d_k and ẑ_{k−1} + Δz_k, so a layer whose isopycnal spacing is thinner
 * than its minimum is held at its minimum thickness and holds whatever
 * water lies there, and one whose isopycnals are further apart follows
 * them (Bleck 2002). An isopycnal interface moves a fraction `relax` of
 * the way from where it is to its target in a step, a held one all the
 * way. Every layer then keeps at least hMin of water, the interfaces being
 * spaced from the top and clipped from the floor, so the layers whose
 * minimum depth lies below a shallow floor collapse to hMin at the
 * bottom. The mode weight w_k = clamp((ẑ_k − ẑ_{k−1} − Δz_k)/modeDepth,
 * 0, 1) says how isopycnal an interface is.
 *
 * Remap (remapColumn): a conservative piecewise-parabolic reconstruction
 * (Colella and Woodward 1984) of a layer-mean quantity on the K old
 * layers, integrated over the new ones (`layers` of them, K unless given). Interior interface values are the slope
 * at the interface of the quartic through the five running integrals of
 * the four layers around it, the next-to-boundary ones the thickness
 * weighted mean of the two layers, the surface and floor values the
 * boundary layers' own means; each is held between the two means it
 * separates and each layer's parabola is limited to lie between its
 * interface values, so no new extremum appears and water, heat and salt
 * are conserved to rounding. Order 2 is the piecewise-linear form with
 * the monotonized central slope, for velocity.
 *
 * Column physics (mixColumn), per step and in this order: the surface
 * heat and salt into layer 0 and the sunlight over the layers by the
 * Jerlov transmission across their interfaces (sunlightShares, the light
 * passing the floor staying in the bottom layer); convective adjustment
 * (convectiveAdjustment), which merges statically unstable neighbors into
 * uniform slabs and measures the potential energy each slab released; the
 * energetic boundary layer (boundaryLayer): the depth h at which the work
 * of mixing the column over the step with the diffusivity
 * κ v* z (1 − z/h)², v*³ = u*³ + E_c/(ρ₀ dt), equals the energy
 * ρ₀ m* u*³ dt + n* E_c, E_c the energy the surface slab's convection
 * released, m* = mstar Λ / (1 + (h |f| / (ekman u*))²) with Λ the
 * Langmuir factor of the wind over the depth h. The work is the potential
 * energy the implicit mixing raises, interface by interface from the heat
 * and salt it moves and the thermal expansion and haline contraction
 * there (Reichl and Hallberg 2018 in its budget, Kraus and Turner in its
 * sources), each trial solving only the layers the mixing reaches; h is
 * the root in log depth of (work − energy)/(work + energy), a bounded
 * function that is near ±1 away from the root, found by the
 * Anderson–Björck secant to a part in 10⁶ in a bracket around the
 * previous step's depth, or between the floor and the whole column, the
 * side where the work is within the energy taken; below and beside it the Pacanowski and Philander
 * (1981) interior viscosity ν = ν₀/(1 + 5Ri)² + ν_b and diffusivity
 * ν/(1 + 5Ri) + κ_b on the gradient Richardson number of the resolved
 * shear and stratification, a convective diffusivity where the column is
 * still unstable; and one implicit solve (solveColumn) of the column's
 * heat and salt with the summed diffusivities, which conserves both. The
 * viscosity profile is published for the edges' solve of their normal
 * velocity, which the same solver does; the wind's stress enters the top
 * layer's tendency in the dynamics, not this solve. The mixed-layer depth (mixedLayerDepth) is where the
 * density first exceeds the surface layer's by an excess.
 */
/*
 * The Langmuir enhancement of the wind's mixing: the theory wave of a fully
 * developed sea under a 10 m wind U10 = u* / FRICTION_PER_WIND (air 1.2
 * kg/m³, drag 1.2e-3), its Stokes drift averaged over the top fifth of a
 * layer `depth` thick in the form of CVMix's cvmix_kpp_ustokes_SL_model
 * (Li, Fox-Kemper, Breivik and Webb 2017), and the factor
 * √(1 + La⁻²/1.5² + La⁻⁴/5.4⁴) of Li et al. (2016) for its surface-layer
 * Langmuir number La = √(u* / u_SL), La no less than 0.1, the factor no more
 * than `cap`. erfc is Abramowitz and Stegun's 7.1.26, as on the GPU.
 */
export const FRICTION_PER_WIND = 1.19e-3, CALM_WIND = 1e-3, RETREAT_EKMAN = 0.7;
export function erfc(x) {
  const t = 1 / (1 + 0.3275911 * x);
  return t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * Math.exp(-x * x);
}
export function stokesSurfaceLayer(u10, depth) {
  const drift = 0.0162 * u10, height = 0.0246 * u10 * u10, mean = 1.296 * 0.877 * 9.81 / (2 * Math.PI * 1.075 * u10);
  const k = 0.176 * drift / (0.125 * Math.PI * 0.667 * mean * height * height), ks = 2.56 * k, z = 0.2 * depth;
  const r1 = (0.151 / (k * z) - 0.84) * (1 - Math.exp(-2 * k * z));
  const r2 = -(0.84 + 0.0591 / (k * z)) * Math.sqrt(2 * Math.PI * k * z) * erfc(Math.sqrt(2 * k * z));
  const r3 = (0.0632 / (ks * z) + 0.125) * (1 - Math.exp(-2 * ks * z));
  const r4 = (0.125 + 0.0946 / (ks * z)) * Math.sqrt(2 * Math.PI * ks * z) * erfc(Math.sqrt(2 * ks * z));
  return drift * (0.715 + r1 + r2 + r3 + r4);
}
export function langmuirFactor(ustar, depth, cap) {
  const u10 = ustar / FRICTION_PER_WIND;
  if (!(u10 > CALM_WIND)) return 1;
  const inverse = Math.min(100, Math.max(0, stokesSurfaceLayer(u10, depth) / ustar));
  return Math.min(cap, Math.sqrt(1 + inverse / 2.25 + inverse * inverse / 850.3056));
}

/*
 * The share of the absorbed surface shortwave that passes depth z, two
 * bands R e^(−z/ζ₁) + (1 − R) e^(−z/ζ₂) with [R, ζ₁, ζ₂] by Jerlov water
 * type (Paulson and Simpson 1977, table 2), less its value at
 * PENETRATION_DEPTH so that none passes below that.
 */
export const JERLOV = { I: [0.58, 0.35, 23], IA: [0.62, 0.6, 20], IB: [0.67, 1.0, 17], II: [0.77, 1.5, 14], III: [0.78, 1.4, 7.9] };
export const PENETRATION_DEPTH = 300;
const jerlovShare = ([R, z1, z2], z) => R * Math.exp(-z / z1) + (1 - R) * Math.exp(-z / z2);
export function transmitted(z, type = 'IB') {
  const bands = JERLOV[type];
  return jerlovShare(bands, Math.min(z, PENETRATION_DEPTH)) - jerlovShare(bands, PENETRATION_DEPTH);
}
export function deposited(z, depth, type = 'IB') {
  const [R, z1, z2] = JERLOV[type];
  return z < Math.min(depth, PENETRATION_DEPTH) ? R / z1 * Math.exp(-z / z1) + (1 - R) / z2 * Math.exp(-z / z2) : 0;
}

export const GRAVITY = 9.81, H_MIN = 0.5, MODE_DEPTH = 5, DZ_FLOOR = 1e-3, H_NEGLECT = 1e-6;
export const GENERATORS = ['isopycnal', 'zstar', 'hybrid'];

export function hybridLayerTable({ top = 3, growth = 1.08, thermocline = 10, deepTop = 40, deepBottom = 500 } = {}) {
  const densities = [], minimum = [];
  for (let k = 0; k < 12; k++) { densities.push(1016.5 + 0.25 * k); minimum.push(top * growth ** k); }
  const warm = [1019.5, 1020.0, 1020.5, 1021.0, 1021.25, 1021.5, 1021.625, 1021.75, 1021.875];
  for (let j = 0; j < warm.length; j++) { densities.push(warm[j]); minimum.push(Math.min(thermocline, 8 + j)); }
  for (let j = 0; j <= 16; j++) { densities.push(1022.0 + 0.25 * j); minimum.push(thermocline); }
  for (let j = 1; j <= 6; j++) { densities.push(Number((1026.0 + 0.1 * j).toFixed(1))); minimum.push(thermocline); }
  const deep = [1026.65, 1026.68, 1026.71, 1026.74, 1026.77, 1026.8, 1026.83, 1026.86, 1026.89, 1026.92, 1026.95, 1026.98, 1027.1, 1027.5];
  for (let j = 0; j < deep.length; j++) { densities.push(deep[j]); minimum.push(deepTop * (deepBottom / deepTop) ** (j / (deep.length - 1))); }
  const K = densities.length, depths = new Float64Array(K + 1);
  for (let k = 0; k < K; k++) depths[k + 1] = depths[k] + minimum[k];
  return { densities: Float64Array.from(densities), minimum: Float64Array.from(minimum), depths, K, fixed: 12, spacing: thermocline };
}
export const HYBRID_LAYERS = hybridLayerTable();

export function labelBoundaries(densities) {
  const K = densities.length, out = new Float64Array(K + 1);
  out[0] = -Infinity;
  out[K] = Infinity;
  for (let k = 1; k < K; k++) out[k] = 0.5 * (densities[k - 1] + densities[k]);
  return out;
}

export function columnScratch(K) {
  const array = (n) => new Float64Array(n);
  return {
    K, centers: array(K), sorted: array(K), flat: array(K), interfaces: array(K + 1), target: array(K + 1), weight: array(K + 1),
    edges: array(K + 1), left: array(K), slope: array(K), curve: array(K), zOld: array(K + 1), zNew: array(K + 1), knots: array(5), primitive: array(5),
    gain: array(K), trialT: array(K), trialS: array(K), alpha: array(K), beta: array(K), rho: array(K),
    boundaryKappa: array(K), interiorKappa: array(K), interiorNu: array(K), start: new Int32Array(K), thick: array(K), heat: array(K), salt: array(K), slabRho: array(K), shares: array(K),
  };
}

export function interfacesOf(h, K, out) {
  out[0] = 0;
  for (let k = 0; k < K; k++) out[k + 1] = out[k] + h[k];
  return out[K];
}

export function layerDensities(T, S, K, out) {
  for (let k = 0; k < K; k++) out[k] = seawaterDensity(T[k], S[k]);
  return out;
}

export function isopycnalInterfaces(h, rho, K, boundaries, out, scratch, flatness = null) {
  const { sorted, zOld, left, slope, curve } = scratch;
  for (let k = 0; k < K; k++) sorted[k] = k > 0 ? Math.max(sorted[k - 1], rho[k]) : rho[k];
  const H = interfacesOf(h, K, zOld);
  reconstruct(h, sorted, K, 3, scratch);
  if (flatness) for (let k = 0; k < K; k++) {
    const keep = 1 - flatness[k];
    left[k] = sorted[k] + keep * (left[k] - sorted[k]);
    slope[k] *= keep;
    curve[k] *= keep;
  }
  out[0] = 0;
  out[K] = H;
  let j = 0;
  for (let k = 1; k < K; k++) {
    const b = boundaries[k];
    while (j < K && left[j] + slope[j] < b) j++;
    if (j >= K) { out[k] = H; continue; }
    if (b <= left[j]) { out[k] = zOld[j]; continue; }
    out[k] = zOld[j] + h[j] * crossing(left[j], slope[j], curve[j], b);
  }
  return out;
}

function crossing(qL, d, c, b) {
  const target = b - qL;
  if (Math.abs(c) < 1e-12 * Math.max(Math.abs(d), 1e-300)) return d > 0 ? Math.min(1, Math.max(0, target / d)) : 0;
  const p = d + c, discriminant = p * p - 4 * c * target;
  let xi = discriminant >= 0 ? (p - Math.sqrt(discriminant)) / (2 * c) : 0.5;
  if (!(xi >= 0 && xi <= 1)) {
    let low = 0, high = 1;
    for (let i = 0; i < 40; i++) { xi = 0.5 * (low + high); if (xi * (d + c * (1 - xi)) < target) low = xi; else high = xi; }
  }
  return xi;
}

export function spaceInterfaces(z, K, H, hMin) {
  const spacing = Math.min(hMin, H / K);
  z[0] = 0;
  z[K] = H;
  for (let k = 1; k < K; k++) z[k] = Math.max(z[k], z[k - 1] + spacing);
  for (let k = K - 1; k >= 1; k--) z[k] = Math.min(z[k], H - (K - k) * spacing);
  return z;
}

export function targetInterfaces(kind, { h, rho, K, depth, table, boundaries, hMin = H_MIN, relax = 1, modeDepth = MODE_DEPTH, previousWeight = null }, out, weight, scratch) {
  if (!GENERATORS.includes(kind)) throw new Error(`the generator is one of ${GENERATORS.join(', ')}, not ${kind}`);
  const H = interfacesOf(h, K, scratch.interfaces), current = scratch.interfaces;
  weight[0] = 0;
  weight[K] = 0;
  if (kind === 'zstar') {
    for (let k = 1; k < K; k++) { out[k] = Math.min(table.depths[k], depth) * H / depth; weight[k] = 0; }
    return spaceInterfaces(out, K, H, hMin);
  }
  const flatness = scratch.flat, mode = (k) => (k <= 0 || k >= K ? 0 : previousWeight ? previousWeight[k] : kind === 'isopycnal' ? 1 : Math.min(1, Math.max(0, (h[k - 1] - table.minimum[k - 1]) / modeDepth)));
  for (let k = 0; k < K; k++) flatness[k] = Math.min(mode(k), mode(k + 1));
  isopycnalInterfaces(h, rho, K, boundaries, out, scratch, flatness);
  let above = 0;
  for (let k = 1; k < K; k++) {
    const floor = kind === 'hybrid' ? above + table.minimum[k - 1] : 0;
    let target = Math.max(out[k], floor);
    if (kind === 'hybrid' && k <= (table.fixed ?? 0)) target = floor;
    else if (target > floor) target = Math.max(floor, current[k] + relax * (target - current[k]));
    out[k] = target;
    weight[k] = kind === 'hybrid' ? Math.min(1, Math.max(0, (target - floor) / modeDepth)) : 1;
    above = target;
  }
  return spaceInterfaces(out, K, H, hMin);
}

function lagrangeSlope(x, P, n, m) {
  let sum = 0;
  for (let i = 0; i < n; i++) {
    if (i === m) {
      let s = 0;
      for (let j = 0; j < n; j++) if (j !== m) s += 1 / (x[m] - x[j]);
      sum += P[m] * s;
    } else {
      let numerator = 1, denominator = 1;
      for (let j = 0; j < n; j++) {
        if (j === i) continue;
        denominator *= x[i] - x[j];
        if (j !== m) numerator *= x[m] - x[j];
      }
      sum += P[i] * numerator / denominator;
    }
  }
  return sum;
}

function reconstruct(hOld, qOld, K, order, scratch) {
  const { edges, left, slope, curve, zOld, knots, primitive } = scratch;
  if (order === 2) {
    for (let k = 0; k < K; k++) {
      let s = 0;
      if (k > 0 && k < K - 1) {
        const below = qOld[k + 1] - qOld[k], above = qOld[k] - qOld[k - 1];
        if (below * above > 0) {
          const central = hOld[k] * (qOld[k + 1] - qOld[k - 1]) / Math.max(0.5 * (hOld[k - 1] + 2 * hOld[k] + hOld[k + 1]), DZ_FLOOR);
          s = Math.sign(central) * Math.min(Math.abs(central), 2 * Math.abs(below), 2 * Math.abs(above));
        }
      }
      left[k] = qOld[k] - 0.5 * s;
      slope[k] = s;
      curve[k] = 0;
    }
    return;
  }
  edges[0] = qOld[0];
  edges[K] = qOld[K - 1];
  for (let k = 1; k < K; k++) {
    const a = k - 1, b = k, low = Math.min(qOld[a], qOld[b]), high = Math.max(qOld[a], qOld[b]);
    let value;
    const wide = k >= 2 && k <= K - 2 && hOld[k - 2] > H_NEGLECT && hOld[k - 1] > H_NEGLECT && hOld[k] > H_NEGLECT && hOld[k + 1] > H_NEGLECT;
    if (wide) {
      knots[0] = zOld[k - 2] - zOld[k];
      knots[1] = zOld[k - 1] - zOld[k];
      knots[2] = 0;
      knots[3] = zOld[k + 1] - zOld[k];
      knots[4] = zOld[k + 2] - zOld[k];
      primitive[0] = 0;
      primitive[1] = hOld[k - 2] * qOld[k - 2];
      primitive[2] = primitive[1] + hOld[k - 1] * qOld[k - 1];
      primitive[3] = primitive[2] + hOld[k] * qOld[k];
      primitive[4] = primitive[3] + hOld[k + 1] * qOld[k + 1];
      value = lagrangeSlope(knots, primitive, 5, 2);
    } else {
      value = (hOld[b] * qOld[a] + hOld[a] * qOld[b]) / Math.max(hOld[a] + hOld[b], H_NEGLECT);
    }
    edges[k] = Math.min(high, Math.max(low, value));
  }
  for (let k = 0; k < K; k++) {
    let qL = edges[k], qR = edges[k + 1];
    const mean = qOld[k];
    if ((qR - mean) * (mean - qL) <= 0) { qL = mean; qR = mean; } else {
      const d = qR - qL, excess = d * (mean - 0.5 * (qL + qR));
      if (excess > d * d / 6) qL = 3 * mean - 2 * qR;
      else if (excess < -d * d / 6) qR = 3 * mean - 2 * qL;
    }
    left[k] = qL;
    slope[k] = qR - qL;
    curve[k] = 6 * (mean - 0.5 * (qL + qR));
  }
}

export function remapColumn(hOld, qOld, K, hNew, amount, scratch, { order = 3, layers = K } = {}) {
  const { zOld, zNew, left, slope, curve } = scratch;
  const H = interfacesOf(hOld, K, zOld);
  interfacesOf(hNew, layers, zNew);
  zNew[layers] = H;
  reconstruct(hOld, qOld, K, order, scratch);
  const piece = (j, lo, hi) => {
    const hj = hOld[j];
    if (!(hj > 0)) return 0;
    const a = Math.min(1, Math.max(0, (lo - zOld[j]) / hj)), b = Math.min(1, Math.max(0, (hi - zOld[j]) / hj));
    const d1 = b - a, d2 = 0.5 * (b * b - a * a), d3 = (b * b * b - a * a * a) / 3;
    return hj * (left[j] * d1 + slope[j] * d2 + curve[j] * (d2 - d3));
  };
  let j = 0;
  for (let i = 0; i < layers; i++) {
    const A = zNew[i], B = zNew[i + 1];
    while (j < K - 1 && zOld[j + 1] <= A) j++;
    let total = 0;
    for (let jj = j; jj < K; jj++) {
      const lo = Math.max(A, zOld[jj]), hi = Math.min(B, zOld[jj + 1]);
      if (hi > lo) total += piece(jj, lo, hi);
      if (zOld[jj + 1] >= B) break;
    }
    amount[i] = total;
  }
  return amount;
}

export function remapMeans(hOld, qOld, K, hNew, out, scratch, options = {}) {
  const layers = options.layers ?? K;
  remapColumn(hOld, qOld, K, hNew, out, scratch, options);
  for (let i = 0; i < layers; i++) out[i] = hNew[i] > 0 ? out[i] / hNew[i] : qOld[Math.min(i, K - 1)];
  return out;
}

/*
 * The implicit diffusion of the first n layers of a column with kappa at
 * their interfaces, kappa[k] between layers k and k + 1 and the interface
 * below layer n − 1 taken as closed: the elimination down the layers with
 * the coefficient and gain of the layer above carried in scalars and the
 * eliminated right-hand sides written to the outputs, then the
 * substitution back up. Two tracers share the matrix, so a pair is solved
 * in one sweep; an output may be its own input.
 */
function solveLayers(h, n, kappa, dt, gain, q1, out1, q2 = null, out2 = null) {
  let wPrev = 0, gPrev = 0, r1 = 0, r2 = 0;
  for (let k = 0; k < n; k++) {
    const w = k < n - 1 ? dt * kappa[k] / Math.max(0.5 * (h[k] + h[k + 1]), DZ_FLOOR) : 0;
    const m = h[k] + wPrev + w + wPrev * gPrev, g = -w / m;
    r1 = (h[k] * q1[k] + wPrev * r1) / m;
    out1[k] = r1;
    if (q2) { r2 = (h[k] * q2[k] + wPrev * r2) / m; out2[k] = r2; }
    gain[k] = g;
    wPrev = w;
    gPrev = g;
  }
  let n1 = out1[n - 1], n2 = q2 ? out2[n - 1] : 0;
  for (let k = n - 2; k >= 0; k--) {
    n1 = out1[k] - gain[k] * n1;
    out1[k] = n1;
    if (q2) { n2 = out2[k] - gain[k] * n2; out2[k] = n2; }
  }
}

const anyMixing = (kappa, K) => { for (let k = 0; k < K - 1; k++) if (kappa[k] > 0) return true; return false; };

export function solveColumn(h, q, K, kappa, dt, scratch) {
  if (!anyMixing(kappa, K)) return false;
  solveLayers(h, K, kappa, dt, scratch.gain, q, q);
  return true;
}

export function interfaceFluxes(h, q, K, kappa, dt, out) {
  for (let k = 0; k < K - 1; k++) out[k] = dt * kappa[k] * (q[k] - q[k + 1]) / Math.max(0.5 * (h[k] + h[k + 1]), DZ_FLOOR);
  return out;
}

export function sunlightShares(z, K, type, out) {
  let sum = 0;
  for (let k = 0; k < K; k++) { out[k] = transmitted(z[k], type) - transmitted(z[k + 1], type); sum += out[k]; }
  out[K - 1] += 1 - sum;
  return out;
}

export function mixedLayerDepth(h, rho, K, excess = 0.03) {
  let z = 0.5 * h[0], previous = rho[0];
  const limit = rho[0] + excess;
  for (let k = 1; k < K; k++) {
    const center = z + 0.5 * (h[k - 1] + h[k]);
    if (rho[k] > limit) return z + (limit - previous) / (rho[k] - previous) * (center - z);
    z = center;
    previous = rho[k];
  }
  let H = 0;
  for (let k = 0; k < K; k++) H += h[k];
  return H;
}

export function convectiveAdjustment(h, T, S, K, scratch, { gravity = GRAVITY } = {}) {
  const { start, thick, heat, salt, slabRho, rho, centers } = scratch;
  let z = 0;
  for (let k = 0; k < K; k++) { centers[k] = z + 0.5 * h[k]; z += h[k]; rho[k] = seawaterDensity(T[k], S[k]); }
  let n = 0;
  for (let k = 0; k < K; k++) {
    start[n] = k; thick[n] = h[k]; heat[n] = h[k] * T[k]; salt[n] = h[k] * S[k]; slabRho[n] = rho[k];
    n++;
    while (n >= 2 && slabRho[n - 2] > slabRho[n - 1]) {
      thick[n - 2] += thick[n - 1];
      heat[n - 2] += heat[n - 1];
      salt[n - 2] += salt[n - 1];
      slabRho[n - 2] = thick[n - 2] > 0 ? seawaterDensity(heat[n - 2] / thick[n - 2], salt[n - 2] / thick[n - 2]) : slabRho[n - 2];
      n--;
    }
  }
  let released = 0, interior = 0, depth = h[0];
  for (let s = 0; s < n; s++) {
    const from = start[s], to = s + 1 < n ? start[s + 1] : K;
    if (to - from < 2) continue;
    const t = heat[s] / thick[s], w = salt[s] / thick[s], mixed = seawaterDensity(t, w);
    let drop = 0;
    for (let k = from; k < to; k++) { drop += gravity * centers[k] * h[k] * (mixed - rho[k]); T[k] = t; S[k] = w; }
    if (from === 0) { released = drop; depth = thick[s]; } else interior += drop;
  }
  return { released, interior, depth };
}

export function boundaryLayer({ h, T, S, K, ustar, f, dt, released = 0, guess = 0, rho0 = SEAWATER.rho0, gravity = GRAVITY, mstar = 1.25, ekman = 0.5, nstar = 0.2, langmuir = 2, vonKarman = 0.4, iterations = 24, tolerance = 1e-6, floor = 0.1 }, kappa, scratch) {
  const { alpha, beta, trialT, trialS, gain } = scratch;
  let H = 0;
  for (let k = 0; k < K; k++) H += h[k];
  for (let k = 0; k < K - 1; k++) {
    alpha[k] = 0.5 * (thermalExpansion(T[k], S[k]) + thermalExpansion(T[k + 1], S[k + 1]));
    beta[k] = 0.5 * (halineContraction(T[k], S[k]) + halineContraction(T[k + 1], S[k + 1]));
  }
  const convective = Math.max(0, released), vstar = Math.cbrt(ustar * ustar * ustar + convective / (rho0 * dt));
  for (let k = 0; k < K; k++) kappa[k] = 0;
  if (!(vstar > 0)) return { depth: 0, energy: 0, work: 0, velocity: 0 };
  const energy = (depth) => {
    const rotation = ustar > 0 ? 1 + (depth * Math.abs(f) / (ekman * ustar)) ** 2 : Infinity;
    const factor = ustar > 0 && langmuir > 0 ? langmuirFactor(ustar, Math.max(depth, 1), langmuir) : 1;
    return rho0 * dt * mstar * factor * ustar * ustar * ustar / rotation + nstar * convective;
  };
  /*
   * The work of mixing to `depth`: the diffusivity at each interface the
   * depth reaches is set as the elimination descends the layers it
   * couples, which stops at the first interface it does not, and the
   * potential energy the trial profile raises is summed as the
   * substitution climbs back. kappa[k] is left for k < solved, the layers
   * the last trial solved.
   */
  let solved = 0;
  const work = (depth) => {
    let z = 0, wPrev = 0, gPrev = 0, rT = 0, rS = 0;
    solved = 0;
    for (let k = 0; k < K; k++) {
      z += h[k];
      const inside = k < K - 1 && z < depth, kap = inside ? vonKarman * vstar * z * (1 - z / depth) ** 2 : 0;
      kappa[k] = kap;
      const w = inside ? dt * kap / Math.max(0.5 * (h[k] + h[k + 1]), DZ_FLOOR) : 0;
      const m = h[k] + wPrev + w + wPrev * gPrev, g = -w / m;
      rT = (h[k] * T[k] + wPrev * rT) / m;
      rS = (h[k] * S[k] + wPrev * rS) / m;
      trialT[k] = rT;
      trialS[k] = rS;
      gain[k] = g;
      wPrev = w;
      gPrev = g;
      solved = k + 1;
      if (!inside) break;
    }
    let raised = 0, nT = trialT[solved - 1], nS = trialS[solved - 1];
    for (let k = solved - 2; k >= 0; k--) {
      const tBelow = nT, sBelow = nS;
      nT = trialT[k] - gain[k] * nT;
      nS = trialS[k] - gain[k] * nS;
      const lift = 0.5 * (h[k] + h[k + 1]), dz = Math.max(lift, DZ_FLOOR);
      const heatDown = dt * kappa[k] * (nT - tBelow) / dz, saltDown = dt * kappa[k] * (nS - sBelow) / dz;
      raised += gravity * rho0 * lift * (alpha[k] * heatDown - beta[k] * saltDown);
    }
    return raised;
  };
  let done = 0, available = 0;
  const excess = (depth) => { done = work(depth); available = energy(depth); const w = Math.max(0, done); return (w - available) / (w + available); };
  /*
   * The search evaluates the excess at the whole column (phase 0) and, when
   * the work exceeds the energy there, at the floor (1), then at the
   * Anderson–Björck secant points of the bracket (6) until it is
   * `tolerance` wide in log depth or `iterations` points have been tried,
   * then once more at the depth settled on (7) so that kappa holds its
   * diffusivity. A guess, the previous step's depth, is tried first (2),
   * then half (3) or twice (4) of it, which brackets the root in three
   * evaluations when it has moved by less than a factor two, else the
   * whole column (5) or the floor completes the bracket.
   */
  const lowest = Math.min(floor, H), half = Math.max(lowest, 0.5 * guess), twice = Math.min(H, 2 * guess);
  let phase = guess > lowest && guess < H ? 2 : 0, a = Math.log(lowest), b = Math.log(H), fa = 0, fb = 0, x = 0, side = 0, count = 0, depth = H;
  for (;;) {
    if (phase === 1) depth = lowest; else if (phase === 2) depth = guess; else if (phase === 3) depth = half; else if (phase === 4) depth = twice; else if (phase === 0 || phase === 5) depth = H;
    else if (phase === 6) {
      x = (a * fb - b * fa) / (fb - fa);
      if (!(x > a && x < b)) x = 0.5 * (a + b);
      if (!(x > a && x < b)) phase = 7; else depth = Math.exp(x);
    }
    if (phase === 7) depth = Math.exp(a);
    const fx = excess(depth);
    if (phase === 0 || phase === 5) { if (fx <= 0) break; fb = fx; phase = phase === 0 ? 1 : 6; }
    else if (phase === 1) { if (fx >= 0) { depth = lowest; break; } fa = fx; phase = 6; }
    else if (phase === 2) { if (fx > 0) { b = Math.log(guess); fb = fx; phase = 3; } else { a = Math.log(guess); fa = fx; phase = 4; } }
    else if (phase === 3) { if (fx < 0) { a = Math.log(half); fa = fx; phase = 6; } else { b = Math.log(half); fb = fx; if (half <= lowest) { depth = lowest; break; } phase = 1; } }
    else if (phase === 4) { if (fx > 0) { b = Math.log(twice); fb = fx; phase = 6; } else { a = Math.log(twice); fa = fx; if (twice >= H) { depth = H; break; } phase = 5; } }
    else if (phase === 6) {
      if (fx > 0) { if (side > 0) { const m = 1 - fx / fb; fa *= m > 0 ? m : 0.5; } b = x; fb = fx; side = 1; }
      else { if (side < 0) { const m = 1 - fx / fa; fb *= m > 0 ? m : 0.5; } a = x; fa = fx; side = -1; }
      count++;
      if (b - a <= tolerance || count >= iterations) phase = 7;
    } else break;
  }
  if (depth <= floor) { depth = 0; solved = 0; done = 0; }
  for (let k = solved; k < K; k++) kappa[k] = 0;
  return { depth, energy: available, work: done, velocity: vstar };
}

export function interiorMixing({ h, T, S, u = null, v = null, shear = null, K, rho0 = SEAWATER.rho0, gravity = GRAVITY, shearViscosity = 1e-2, backgroundViscosity = 1e-4, backgroundDiffusivity = 1e-5, convectiveDiffusivity = 0.1 }, kappa, nu, scratch) {
  const { rho } = scratch;
  layerDensities(T, S, K, rho);
  for (let k = 0; k < K - 1; k++) {
    const dz = Math.max(0.5 * (h[k] + h[k + 1]), DZ_FLOOR), buoyancy = gravity * (rho[k + 1] - rho[k]) / rho0;
    if (buoyancy < 0) { kappa[k] = convectiveDiffusivity; nu[k] = convectiveDiffusivity; continue; }
    const du = u ? u[k] - u[k + 1] : 0, dv = v ? v[k] - v[k + 1] : 0, squared = shear ? shear[k] : du * du + dv * dv;
    const richardson = buoyancy * dz / (squared + 1e-12);
    const viscosity = shearViscosity / (1 + 5 * richardson) ** 2 + backgroundViscosity;
    nu[k] = viscosity;
    kappa[k] = viscosity / (1 + 5 * richardson) + backgroundDiffusivity;
  }
  kappa[K - 1] = 0;
  nu[K - 1] = 0;
}

export function mixColumn({ h, T, S, u = null, v = null, shear = null, K, heat = 0, salt = 0, light = 0, type = 'IB', ustar, f, dt, capped = false, options = {} }, { kappa, nu }, scratch) {
  const { interfaces, shares, boundaryKappa, interiorKappa, interiorNu } = scratch;
  interfacesOf(h, K, interfaces);
  if (light !== 0) {
    sunlightShares(interfaces, K, type, shares);
    for (let k = 0; k < K; k++) T[k] += light * shares[k] / h[k];
  }
  T[0] += heat / h[0];
  S[0] += salt / h[0];
  const convection = convectiveAdjustment(h, T, S, K, scratch, options);
  const layer = boundaryLayer({ h, T, S, K, ustar, f, dt, released: capped ? 0 : convection.released, ...options }, boundaryKappa, scratch);
  interiorMixing({ h, T, S, u, v, shear, K, ...options }, interiorKappa, interiorNu, scratch);
  for (let k = 0; k < K; k++) { kappa[k] = boundaryKappa[k] + interiorKappa[k]; nu[k] = boundaryKappa[k] + interiorNu[k]; }
  if (anyMixing(kappa, K)) solveLayers(h, K, kappa, dt, scratch.gain, T, T, S, S);
  return { depth: layer.depth, energy: layer.energy, work: layer.work, released: convection.released, convectiveDepth: convection.depth };
}
