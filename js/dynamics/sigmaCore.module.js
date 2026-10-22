import { divergence, gradient, curl, kineticEnergy, laplacianVelocity, laplacianScalar } from './operators.module.js';
import { spongeGeometry, dampEddies } from './sponge.module.js';

export const R_DRY = 287.06;
export const VIRTUAL_FACTOR = 0.608;
export const CP_DRY = 1003.5;
export const P0 = 101325;
export const GRAVITY = 9.806;
export const DIVERGENCE_DAMPING = 0.03;

const CAM_L26_HYAI = [
  0.00219406700000001, 0.00489520900000001, 0.009882418, 0.01805201,
  0.02983724, 0.0446233400000002, 0.0616058700000002, 0.0785124300000004,
  0.0773127100000002, 0.0759013100000003, 0.0742408600000002,
  0.0722874400000002, 0.0699893299999998, 0.06728574, 0.06410509,
  0.0603632200000002, 0.0559611100000001, 0.0507822500000001,
  0.0446896000000001, 0.0375219099999999, 0.0290894900000001, 0.02084739,
  0.01334443, 0.00708499000000001, 0.00252136, 0, 0,
];
const CAM_L26_HYBI = [
  0, 0, 0, 0, 0, 0, 0, 0, 0.01505309, 0.03276228, 0.05359622,
  0.0781062700000006, 0.1069411, 0.140863700000001, 0.180772, 0.227722,
  0.282956200000001, 0.347936400000002, 0.4243822, 0.514316800000003,
  0.620120200000002, 0.723535500000004, 0.817676800000001,
  0.896215300000001, 0.953476100000003, 0.9851122, 1,
];

export const DIVERGENCE_SURFACES = ['sigma', 'hybrid'];

/*
 * The models' `dynamics` options with their defaults: divergenceSurfaces,
 * the surfaces the divergence damping acts along, 'sigma' or 'hybrid'.
 */
export function dynamicsOptions({ divergenceSurfaces = 'sigma', ...rest } = {}) {
  const unknown = Object.keys(rest);
  if (unknown.length) throw new Error(`the dynamics options are divergenceSurfaces, not ${unknown.join(', ')}`);
  if (!DIVERGENCE_SURFACES.includes(divergenceSurfaces)) throw new Error(`divergenceSurfaces is ${DIVERGENCE_SURFACES.join(' or ')}, not ${divergenceSurfaces}`);
  return { divergenceSurfaces };
}

/*
 * The weight w = 1 − B/η of the hybrid σ–p surface through each σ, with
 * B(η) CAM L26's HYBI, linear in η = HYAI + HYBI (Collins et al. 2004),
 * at η = σ: 1 where CAM's coordinate is pure pressure, 0 at the ground.
 * Along such a surface d ln p = (B/η) d ln π, so d ln σ = −w d ln π.
 */
export function hybridWeights(sigma) {
  const eta = CAM_L26_HYAI.map((a, k) => a + CAM_L26_HYBI[k]);
  return Float64Array.from(sigma, (s) => {
    if (s <= eta[0]) return 1;
    let k = 0;
    while (k < eta.length - 2 && eta[k + 1] < s) k++;
    const t = (s - eta[k]) / (eta[k + 1] - eta[k]);
    return 1 - (CAM_L26_HYBI[k] + t * (CAM_L26_HYBI[k + 1] - CAM_L26_HYBI[k])) / s;
  });
}

/*
 * 1/(ln σ_below − ln σ_above) for the derivative in ln σ at each layer's
 * middle, centered over the layers on either side and one-sided at the
 * top and bottom layers.
 */
export function logSigmaSteps(sigmaMid) {
  const K = sigmaMid.length;
  return Float64Array.from(sigmaMid, (_, k) => 1 / Math.log(sigmaMid[Math.min(k + 1, K - 1)] / sigmaMid[Math.max(k - 1, 0)]));
}

const BOUNDARY_LAYER_HEIGHTS = [1698, 1213, 882, 649, 482, 359, 265, 192, 133, 83, 40];

/*
 * σ at a geopotential height (m) in the US Standard Atmosphere 1976, and
 * the height at a σ. Its layers from the ground up, each as its base
 * height, temperature, lapse rate (K/m, positive where it cools upward)
 * and σ, the base temperature and σ carried up from the layer below.
 */
const STANDARD = { T0: 288.15, lapse: 0.0065, g: 9.80665, R: 287.053 };
const STANDARD_LAYERS = [[11000, 0], [20000, -0.001], [32000, -0.0028], [47000, 0], [51000, 0.0028], [71000, 0.002], [84852, 0]].reduce((layers, [z, lapse]) => {
  const [z0, T0, lapse0] = layers.at(-1);
  return [...layers, [z, T0 - lapse0 * (z - z0), lapse, sigmaInLayer(layers.at(-1), z)]];
}, [[0, STANDARD.T0, STANDARD.lapse, 1]]);
function sigmaInLayer([z0, T0, lapse, sigma0], z) {
  return lapse === 0 ? sigma0 * Math.exp(-STANDARD.g * (z - z0) / (STANDARD.R * T0)) : sigma0 * Math.pow(1 - lapse * (z - z0) / T0, STANDARD.g / (STANDARD.R * lapse));
}
export function standardSigma(z) {
  return sigmaInLayer(STANDARD_LAYERS.findLast((layer) => layer[0] <= z) ?? STANDARD_LAYERS[0], z);
}
export function standardHeight(sigma) {
  const [z0, T0, lapse, sigma0] = STANDARD_LAYERS.findLast((layer) => layer[3] >= sigma) ?? STANDARD_LAYERS[0];
  return lapse === 0 ? z0 + STANDARD.R * T0 / STANDARD.g * Math.log(sigma0 / sigma) : z0 + T0 / lapse * (1 - Math.pow(sigma / sigma0, STANDARD.R * lapse / STANDARD.g));
}

export const SIGMA_GRIDS = ['cam26', 'bl34', 'bl36', 'bl42', 'bl48', 'bl60'];
const STRATOPAUSE_INTERFACES = [30, 100];
const STRATOSPHERE = { base: 0.09, layers: 27 };
const INVERSION_HEIGHTS = [4100, 3498, 3160, 2830, 2520, 2230, 1960, 1698, 1450, ...BOUNDARY_LAYER_HEIGHTS.slice(1)];
const LADDER = { bottom: 40, ratio: 2 ** (1 / 5), troposphere: [25, 7], stratosphere: [5, 9] };

/*
 * The interfaces from the ground up, as σ, of `growing` layers in the
 * standard atmosphere from σ `from`, the first `first` m thick and each
 * after it LADDER.ratio times the one below, then `equal` layers of one
 * thickness up to σ `to`; and that thickness.
 */
function ladder(from, to, first, growing, equal) {
  const heights = [standardHeight(from)];
  for (let j = 0; j < growing; j++) heights.push(heights[j] + first * LADDER.ratio ** j);
  const step = (standardHeight(to) - heights[growing]) / equal;
  for (let j = 1; j < equal; j++) heights.push(heights[growing] + j * step);
  return { sigma: [...heights.slice(1).map(standardSigma), to], step };
}

/*
 * Interface sigma values, top (0) to ground (1), of the grid by that
 * name. 'cam26' is the CAM 26-level hybrid grid that Jablonowski &
 * Williamson 2006 ran on, read as sigma with p_s = p0 (hyai + hybi from
 * HOMME's cami-26.ascii), plus the cap above CAM's 2.19 hPa lid that a
 * sigma coordinate needs: 27 layers, the lowest three within about 850 m
 * of the ground. 'bl34' is cam26 down to its σ 0.7444 interface (2.4 km)
 * and below that the interfaces at BOUNDARY_LAYER_HEIGHTS in the
 * standard atmosphere: a 40 m lowest layer, each layer 7 % to 49 %
 * thicker than the one below it up into the first cam26 layer, ten
 * layers below 1.25 km, 34 in all. 'bl36' is bl34 with its 0-2.19 hPa
 * top layer split at STRATOPAUSE_INTERFACES (Pa over p0): 0-0.3,
 * 0.3-1 and 1-2.19 hPa, so the stratopause near 1 hPa lies between
 * layer midpoints and the mesosphere above it has a layer of its own.
 * 'bl42' is bl36 with its 1.2-4.7 km span cut at INVERSION_HEIGHTS:
 * layers of 237 to 338 m up to 3.5 km and two of about 600 m to the
 * cam26 interface at 4.7 km (standard-atmosphere heights, about 8 %
 * higher in the tropics), so that a trade inversion falls within one
 * thin layer instead of across a 500-1100 m one. 'bl60' is bl42 with
 * its 1-93.6 hPa span cut into STRATOSPHERE.layers layers equal in ln σ,
 * 1.18 km each in log-pressure height (H = 7 km): the 1.2 km of MiMA's
 * lower stratosphere (Garfinkel et al. 2022) carried up to the
 * stratopause, where bl42's are 1.1-5.7 km. 'bl48' keeps bl60's lid and
 * 93.6 hPa interface and lays out the rest in standard-atmosphere
 * heights by LADDER: from a 40 m lowest layer, 25 layers each 2^(1/5)
 * times as thick as the one below, to 8.34 km, then 7 of 1.18 km to
 * 93.6 hPa; above it 5 more that keep growing by 2^(1/5) from 1.18 km,
 * to 25.7 km, then 9 of 2.46 km to 1 hPa. No layer below 1 hPa is more
 * than 1.15 times as thick as the one below it, and those of 1.4-2.5 km
 * in the stratosphere are finer than GFDL AM3's 3-4 km (Donner et al.
 * 2011): 48 layers.
 */
export function sigmaInterfaces(name = 'cam26') {
  const levels = [0];
  for (let k = 0; k < CAM_L26_HYAI.length; k++) levels.push(CAM_L26_HYAI[k] + CAM_L26_HYBI[k]);
  levels[levels.length - 1] = 1;
  if (name === 'cam26') return Float64Array.from(levels);
  if (name === 'bl48') {
    const top = STRATOPAUSE_INTERFACES.map((p) => p / P0), base = levels.find((sigma) => sigma > STRATOSPHERE.base);
    const troposphere = ladder(1, base, LADDER.bottom, ...LADDER.troposphere), stratosphere = ladder(base, top.at(-1), troposphere.step * LADDER.ratio, ...LADDER.stratosphere);
    return Float64Array.from([0, ...top, ...stratosphere.sigma.slice(0, -1).reverse(), ...troposphere.sigma.reverse(), 1]);
  }
  if (name === 'bl34' || name === 'bl36' || name === 'bl42' || name === 'bl60') {
    const refined = (name === 'bl42' || name === 'bl60' ? INVERSION_HEIGHTS : BOUNDARY_LAYER_HEIGHTS).map(standardSigma), top = name === 'bl34' ? [] : STRATOPAUSE_INTERFACES.map((p) => p / P0);
    let free = levels.filter((sigma) => sigma > 0 && sigma < refined[0]);
    if (name === 'bl60') {
      const from = top.at(-1), to = free.find((sigma) => sigma > STRATOSPHERE.base), n = STRATOSPHERE.layers;
      free = [...Array.from({ length: n - 1 }, (_, j) => from * (to / from) ** ((j + 1) / n)), ...free.filter((sigma) => sigma >= to)];
    }
    return Float64Array.from([0, ...top, ...free, ...refined, 1]);
  }
  throw new Error(`no sigma grid is named ${name}; the grids are ${SIGMA_GRIDS.join(', ')}`);
}

/*
 * The name of the grid with these interfaces, to within single-precision
 * rounding, or null.
 */
export function sigmaGridName(levels) {
  return SIGMA_GRIDS.find((name) => {
    const grid = sigmaInterfaces(name);
    return grid.length === levels.length && grid.every((sigma, k) => Math.abs(sigma - levels[k]) <= 1e-6);
  }) ?? null;
}

/*
 * Hydrostatic primitive equations in sigma coordinates on the C-grid.
 * State: pi[C] surface pressure, theta[K*C] layer potential temperature,
 * u[K*E] layer normal velocity (layer k occupies [k*C, (k+1)*C) and
 * [k*E, (k+1)*E)), and optionally state[4] = q[K*C], specific humidity,
 * and state[5] = qc[K*C], cloud condensate, both carried with the same
 * flux-form transport as theta. Interface
 * quantities use K+1 slots per column, top first. Steps follow the
 * A-grid column: mass flux and divergence per layer, dpi/dt,
 * pi*sigma-dot telescoping to zero at the ground, Exner ratios from the
 * exact layer integral of sigma^kappa, interface theta interpolated in
 * Exner, geopotential integrated upward with each layer's virtual theta
 * over its own Exner span, then flux-form theta and q transport and the
 * vector-invariant momentum equation with the RTSK PV flux. The
 * pressure-gradient force uses virtual theta (with condensate loading;
 * with no q it equals theta), and theta takes the pressure work that
 * adds (virtualWork), so the core keeps c_p T + K.
 */
export function createSigmaCore(mesh, options = {}) {
  const {
    levels = sigmaInterfaces(), g = GRAVITY, cp = CP_DRY, R = R_DRY, p0 = P0,
    nu4 = 0, nu4Theta = 0, divergenceDamping = 0, divergenceSurfaces = 'sigma', surfaceGeopotential = null, buffers = null, splitClosure = false, spongeRates = null, spongeMeanRates = null,
  } = options;
  if (!DIVERGENCE_SURFACES.includes(divergenceSurfaces)) throw new Error(`divergenceSurfaces is ${DIVERGENCE_SURFACES.join(' or ')}, not ${divergenceSurfaces}`);
  const {
    nCells: C, nEdges: E, nVertices: V, maxEdgesOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge,
    cellsOnEdge, verticesOnEdge, cellsOnVertex, kiteAreasOnVertex, areaTriangle, dcEdge, dvEdge, fVertex,
  } = mesh;
  const K = levels.length - 1;
  const kappa = R / cp;
  if (divergenceDamping > 0 && !splitClosure) throw new Error('divergence damping needs a core with splitClosure');
  const sponged = [spongeRates, spongeMeanRates].some((rates) => rates && rates.some((rate) => rate > 0));
  if (sponged && !splitClosure) throw new Error('the sponge needs a core with splitClosure');
  const sigmaUpper = levels.subarray(0, K);
  const sigmaLower = levels.subarray(1, K + 1);
  const dSigma = Float64Array.from(sigmaLower, (s, k) => s - sigmaUpper[k]);
  const sigmaMid = Float64Array.from(sigmaLower, (s, k) => 0.5 * (s + sigmaUpper[k]));

  const shared = {};
  const sharedArray = (name, n) => {
    const buffer = buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * n);
    shared[name] = buffer;
    return new Float64Array(buffer);
  };
  const flux = sharedArray('flux', K * E);
  const dissipation = sharedArray('dissipation', K * E);
  const uBefore = new Float64Array(E);
  const divFlux = sharedArray('divFlux', K * C);
  const piSigmaDot = sharedArray('piSigmaDot', (K + 1) * C);
  const exnerLower = sharedArray('exnerLower', K * C);
  const exnerLayer = sharedArray('exnerLayer', K * C);
  const dExnerDpi = sharedArray('dExnerDpi', K * C);
  const thetaLower = sharedArray('thetaLower', K * C);
  const qLower = sharedArray('qLower', K * C);
  const qcLower = sharedArray('qcLower', K * C);
  const thetaV = sharedArray('thetaV', K * C);
  const geopotential = sharedArray('geopotential', K * C);
  const piVertex = sharedArray('piVertex', V);
  const piEdge = new Float64Array(E);
  const thetaFlux = new Float64Array(E);
  const divThetaFlux = new Float64Array(C);
  const zeta = new Float64Array(V);
  const qVertex = new Float64Array(V);
  const qEdge = new Float64Array(E);
  const kinetic = new Float64Array(C);
  const phi = new Float64Array(C);
  const gradPhi = new Float64Array(E);
  const lnPi = new Float64Array(C), gradLnPi = new Float64Array(E);
  const lap = new Float64Array(E);
  const lap2 = new Float64Array(E);
  const divScratch = new Float64Array(C);
  const curlScratch = new Float64Array(V);
  const lapTheta = new Float64Array(C);
  const lapTheta2 = new Float64Array(C);
  const qFlux = new Float64Array(E);
  const divQFlux = new Float64Array(C);
  const qcFlux = new Float64Array(E);
  const divQcFlux = new Float64Array(C);
  const massField = new Float64Array(C);
  const pvWeights = Float64Array.from(weightsOnEdge, (w, slot) => w * dvEdge[edgesOnEdge[slot]]);
  let spacing = 0;
  for (let e = 0; e < E; e++) spacing += dcEdge[e];
  spacing /= E;
  const divergenceStep = divergenceDamping * spacing * spacing;
  const gradDiv = new Float64Array(E);
  const coupled = divergenceSurfaces === 'hybrid' && divergenceStep > 0;
  const surfaceWeight = coupled ? hybridWeights(sigmaMid) : null, lnSigmaStep = coupled ? logSigmaSteps(sigmaMid) : null;
  const hybridDivergence = coupled ? new Float64Array(K * C) : null, halfLnRatio = coupled ? new Float64Array(E) : null;
  const sponge = sponged ? spongeGeometry(mesh) : null, bandMeans = sponge ? new Float64Array(2 * sponge.bands) : null;

  function diagnoseColumn(i, pi, theta, q = null, qc = null) {
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      const exLower = Math.pow(pi[i] * sigmaLower[k] / p0, kappa);
      const exUpper = k === 0 ? 0 : exnerLower[idx - C];
      const span = exLower * sigmaLower[k] - exUpper * sigmaUpper[k];
      exnerLower[idx] = exLower;
      exnerLayer[idx] = span / ((1 + kappa) * dSigma[k]);
      dExnerDpi[idx] = (kappa / (1 + kappa)) * span / (pi[i] * dSigma[k]);
    }
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      thetaV[idx] = q ? theta[idx] * (1 + VIRTUAL_FACTOR * q[idx] - (qc ? qc[idx] : 0)) : theta[idx];
    }
    for (let k = 0; k < K - 1; k++) {
      const idx = k * C + i;
      const t = (exnerLower[idx] - exnerLayer[idx]) / (exnerLayer[idx + C] - exnerLayer[idx]);
      thetaLower[idx] = theta[idx] + t * (theta[idx + C] - theta[idx]);
      qLower[idx] = q ? q[idx] + t * (q[idx + C] - q[idx]) : 0;
      qcLower[idx] = qc ? qc[idx] + t * (qc[idx + C] - qc[idx]) : 0;
    }
    const bottom = (K - 1) * C + i;
    geopotential[bottom] = (surfaceGeopotential ? surfaceGeopotential[i] : 0) + cp * thetaV[bottom] * (exnerLower[bottom] - exnerLayer[bottom]);
    for (let k = K - 2; k >= 0; k--) {
      const idx = k * C + i;
      const below = idx + C;
      geopotential[idx] = geopotential[below]
        + cp * thetaV[below] * (exnerLayer[below] - exnerLower[idx])
        + cp * thetaV[idx] * (exnerLower[idx] - exnerLayer[idx]);
    }
  }

  function diagnose(pi, theta, q = null, qc = null) {
    for (let i = 0; i < C; i++) diagnoseColumn(i, pi, theta, q, qc);
  }

  function edgePi(pi) {
    for (let e = 0; e < E; e++) piEdge[e] = 0.5 * (pi[cellsOnEdge[2 * e]] + pi[cellsOnEdge[2 * e + 1]]);
  }

  function phaseFlux(state, kFrom, kTo) {
    const [pi, , u] = state;
    edgePi(pi);
    for (let k = kFrom; k < kTo; k++) {
      const fk = flux.subarray(k * E, (k + 1) * E);
      const uk = u.subarray(k * E, (k + 1) * E);
      for (let e = 0; e < E; e++) fk[e] = piEdge[e] * uk[e];
      divergence(mesh, fk, divFlux.subarray(k * C, (k + 1) * C));
    }
  }

  function phaseColumn(state, out, iFrom, iTo) {
    const [pi, theta] = state;
    const q = state[4] ?? null, qc = state[5] ?? null;
    const [dPi] = out;
    for (let i = iFrom; i < iTo; i++) {
      let sum = 0;
      for (let k = 0; k < K; k++) sum += divFlux[k * C + i] * dSigma[k];
      dPi[i] = -sum;
      let cumulative = 0;
      piSigmaDot[i] = 0;
      for (let k = 0; k < K; k++) {
        cumulative += divFlux[k * C + i] * dSigma[k];
        piSigmaDot[(k + 1) * C + i] = -cumulative - sigmaLower[k] * dPi[i];
      }
      piSigmaDot[K * C + i] = 0;
      diagnoseColumn(i, pi, theta, q, qc);
    }
  }

  function phaseVertex(state, vFrom, vTo) {
    const [pi] = state;
    for (let v = vFrom; v < vTo; v++) {
      let sum = 0;
      for (let m = 0; m < 3; m++) sum += kiteAreasOnVertex[3 * v + m] * pi[cellsOnVertex[3 * v + m]];
      piVertex[v] = sum / areaTriangle[v];
    }
  }

  /*
   * The ∇⁴ closure of one layer of a scalar field as a tendency (or, with
   * a time step, as an explicit increment of the field). `conservative`
   * applies it to the mass-weighted field so the column's total is
   * preserved exactly under it (used for water). Theta (`exner` given)
   * takes ∇·(π_e ∇∇²T)/π on its temperature T = θΠ, whose fluxes between
   * cells are weighted by the edge's π, so the closure keeps the layer's
   * c_p T Δm exactly.
   */
  function scalarClosure(k, pi, field, target, coefficient, conservative, exner = null) {
    const off = k * C;
    if (exner) {
      for (let i = 0; i < C; i++) massField[i] = field[off + i] * exner[off + i];
      laplacianScalar(mesh, massField, lapTheta);
      laplacianScalar(mesh, lapTheta, lapTheta2, pi);
      for (let i = 0; i < C; i++) target[off + i] -= coefficient * lapTheta2[i] / (pi[i] * exner[off + i]);
    } else if (conservative) {
      for (let i = 0; i < C; i++) massField[i] = pi[i] * field[off + i];
      laplacianScalar(mesh, massField, lapTheta);
      laplacianScalar(mesh, lapTheta, lapTheta2);
      for (let i = 0; i < C; i++) target[off + i] -= coefficient * lapTheta2[i] / pi[i];
    } else {
      laplacianScalar(mesh, field.subarray(off, off + C), lapTheta);
      laplacianScalar(mesh, lapTheta, lapTheta2);
      for (let i = 0; i < C; i++) target[off + i] -= coefficient * lapTheta2[i];
    }
  }

  function momentumClosure(k, u, target, coefficient) {
    const uk = u.subarray(k * E, (k + 1) * E);
    laplacianVelocity(mesh, uk, lap, divScratch, curlScratch);
    laplacianVelocity(mesh, lap, lap2, divScratch, curlScratch);
    const off = k * E;
    for (let e = 0; e < E; e++) target[off + e] -= coefficient * lap2[e];
  }

  /*
   * Divergence damping of one layer over one step: u += c d² ∇(∇·u) with
   * d the mean distance between cell centres, the tendency ν_d ∇δ with
   * ν_d = c d²/dt. The gradient of a divergence has no curl on the C-grid,
   * so the vorticity is untouched.
   */
  function divergenceClosure(k, u) {
    const off = k * E, uk = u.subarray(off, off + E);
    divergence(mesh, uk, divScratch);
    gradient(mesh, divScratch, gradDiv);
    for (let e = 0; e < E; e++) u[off + e] += divergenceStep * gradDiv[e];
  }

  /*
   * Divergence damping of every layer over one step along the hybrid
   * σ–p surfaces of hybridWeights (Simmons & Burridge 1981). The σ
   * surfaces slope over terrain as much in the stratosphere as at the
   * ground, and damping divergence along them mixes momentum vertically
   * (Mellor & Blumberg 1985), K_z ≈ H² ν_d |∇ln π|². A cell's
   * divergence is taken on the surface through its center, each edge's
   * wind moved onto it by w (∂u/∂ln σ)(ln π_i − ln π_e), so
   * D_η = D_σ − w (∂v/∂ln σ)·∇ln π, and its gradient along the surface
   * through the edge, ∇_σ D_η − w (∂D_η/∂ln σ) ∇ln π. The derivatives
   * couple the layers, so all of them are damped at once.
   */
  function hybridDivergenceClosure(pi, u) {
    const { maxEdges, nEdgesOnCell, edgesOnCell, edgeSignOnCell, areaCell } = mesh;
    for (let e = 0; e < E; e++) halfLnRatio[e] = 0.5 * Math.log(pi[cellsOnEdge[2 * e + 1]] / pi[cellsOnEdge[2 * e]]);
    for (let k = 0; k < K; k++) {
      const off = k * E, above = Math.max(k - 1, 0) * E, below = Math.min(k + 1, K - 1) * E, tilt = surfaceWeight[k] * lnSigmaStep[k];
      for (let i = 0; i < C; i++) {
        let flow = 0, tilted = 0;
        for (let m = 0; m < nEdgesOnCell[i]; m++) {
          const slot = maxEdges * i + m, e = edgesOnCell[slot], side = edgeSignOnCell[slot] * dvEdge[e];
          flow += side * u[off + e];
          tilted += side * (cellsOnEdge[2 * e] === i ? -halfLnRatio[e] : halfLnRatio[e]) * (u[below + e] - u[above + e]);
        }
        hybridDivergence[k * C + i] = (flow + tilt * tilted) / areaCell[i];
      }
    }
    for (let k = 0; k < K; k++) {
      const off = k * E, here = k * C, above = Math.max(k - 1, 0) * C, below = Math.min(k + 1, K - 1) * C, w = surfaceWeight[k], step = lnSigmaStep[k];
      for (let e = 0; e < E; e++) {
        const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
        const slope = 0.5 * step * (hybridDivergence[below + a] - hybridDivergence[above + a] + hybridDivergence[below + b] - hybridDivergence[above + b]);
        u[off + e] += divergenceStep * (hybridDivergence[here + b] - hybridDivergence[here + a] - w * slope * 2 * halfLnRatio[e]) / dcEdge[e];
      }
    }
  }

  /*
   * Flux-form transport of a layer scalar by the mass fluxes, with the
   * ∇⁴ closure as part of the tendency unless the closure is split off.
   */
  function transportLayer(k, pi, field, fieldLower, edgeFlux, divField, dField, conservative = false) {
    const off = k * C;
    const fk = flux.subarray(k * E, (k + 1) * E);
    for (let e = 0; e < E; e++) {
      edgeFlux[e] = fk[e] * 0.5 * (field[off + cellsOnEdge[2 * e]] + field[off + cellsOnEdge[2 * e + 1]]);
    }
    divergence(mesh, edgeFlux, divField);
    for (let i = 0; i < C; i++) {
      const idx = off + i;
      const lowerFlow = piSigmaDot[(k + 1) * C + i];
      const upperFlow = piSigmaDot[k * C + i];
      const lower = k === K - 1 ? 0 : fieldLower[idx];
      const upper = k === 0 ? 0 : fieldLower[idx - C];
      const vertical = (lowerFlow * lower - upperFlow * upper - field[idx] * (lowerFlow - upperFlow)) / (pi[i] * dSigma[k]);
      dField[idx] = -(divField[i] - field[idx] * divFlux[idx]) / pi[i] - vertical;
    }
    if (nu4Theta > 0 && !splitClosure) scalarClosure(k, pi, field, dField, nu4Theta, conservative, conservative ? null : exnerLayer);
  }

  /*
   * The pressure work of moist air beyond the dry air's, added to one
   * layer's θ tendency: the pressure-gradient force and the geopotential
   * take θ_v, so the energy equation must take α_v ω where the dry θ
   * equation takes α ω, dθ/dt += κ θ (0.608 q − q_c) ω / p, with
   * ω = σ(∂π/∂t + v·∇π) + πσ̇ at the layer's middle.
   */
  function virtualWork(k, pi, u, theta, q, qc, dPi, dTheta) {
    divergence(mesh, u.subarray(k * E, (k + 1) * E), divScratch);
    for (let i = 0; i < C; i++) {
      const idx = k * C + i;
      const omega = sigmaMid[k] * (dPi[i] + divFlux[idx] - pi[i] * divScratch[i]) + 0.5 * (piSigmaDot[idx] + piSigmaDot[idx + C]);
      dTheta[idx] += kappa * theta[idx] * (VIRTUAL_FACTOR * q[idx] - (qc ? qc[idx] : 0)) * omega / (pi[i] * sigmaMid[k]);
    }
  }

  /*
   * The ∇⁴ closures, and the divergence damping and the top sponge
   * (`spongeRates` on the eddies and `spongeMeanRates` on the zonal
   * mean, 1/s per layer, see sponge.module.js) after the
   * momentum one, applied to the state itself over one time step, for
   * cores built with splitClosure: the model runs this once per step
   * after the RK4 dynamics instead of inside every stage. With
   * divergenceSurfaces 'hybrid' its momentum part takes every layer at once.
   */
  function phaseClosure(state, kFrom, kTo, dt, part = 'all') {
    const [pi, theta, u] = state;
    const q = state[4] ?? null, qc = state[5] ?? null;
    const together = coupled && part !== 'tracers';
    if (together && (kFrom !== 0 || kTo !== K)) throw new Error('the hybrid divergence damping couples the layers: close all of them at once');
    for (let k = kFrom; k < kTo; k++) {
      if (part !== 'momentum' && nu4Theta > 0) {
        scalarClosure(k, pi, theta, theta, dt * nu4Theta, false, exnerLayer);
        if (q) scalarClosure(k, pi, q, q, dt * nu4Theta, true);
        if (qc) scalarClosure(k, pi, qc, qc, dt * nu4Theta, true);
      }
      if (part === 'tracers') continue;
      const off = k * E;
      if (together) {
        for (let e = 0; e < E; e++) dissipation[off + e] = u[off + e] * u[off + e];
        if (nu4 > 0) momentumClosure(k, u, u, dt * nu4);
        continue;
      }
      const spongeRate = sponge && spongeRates ? spongeRates[k] : 0, meanRate = sponge && spongeMeanRates ? spongeMeanRates[k] : 0;
      if (nu4 > 0 || divergenceStep > 0 || spongeRate > 0 || meanRate > 0) {
        for (let e = 0; e < E; e++) uBefore[e] = u[off + e];
        if (nu4 > 0) momentumClosure(k, u, u, dt * nu4);
        if (divergenceStep > 0) divergenceClosure(k, u);
        if (spongeRate > 0 || meanRate > 0) dampEddies(mesh, sponge, u.subarray(off, off + E), spongeRate, dt, bandMeans, meanRate);
        for (let e = 0; e < E; e++) dissipation[off + e] = uBefore[e] * uBefore[e] - u[off + e] * u[off + e];
      } else dissipation.fill(0, off, off + E);
    }
    if (!together) return;
    hybridDivergenceClosure(pi, u);
    for (let k = 0; k < K; k++) {
      const off = k * E;
      const spongeRate = sponge && spongeRates ? spongeRates[k] : 0, meanRate = sponge && spongeMeanRates ? spongeMeanRates[k] : 0;
      if (spongeRate > 0 || meanRate > 0) dampEddies(mesh, sponge, u.subarray(off, off + E), spongeRate, dt, bandMeans, meanRate);
      for (let e = 0; e < E; e++) dissipation[off + e] -= u[off + e] * u[off + e];
    }
  }

  /*
   * The layer phase in two independent parts, so a parallel run can
   * treat a layer's tracer transport and its momentum tendency as
   * separate work units: `part` is 'all', 'tracers' or 'momentum'.
   */
  function phaseLayer(state, out, kFrom, kTo, part = 'all') {
    const [pi, theta, u] = state;
    const [, dTheta, dU] = out;
    const q = state[4] ?? null, dQ = out[4] ?? null, qc = state[5] ?? null, dQc = out[5] ?? null;
    if (part !== 'momentum') {
      for (let k = kFrom; k < kTo; k++) {
        transportLayer(k, pi, theta, thetaLower, thetaFlux, divThetaFlux, dTheta);
        if (q) virtualWork(k, pi, u, theta, q, qc, out[0], dTheta);
        if (q && dQ) transportLayer(k, pi, q, qLower, qFlux, divQFlux, dQ, true);
        if (qc && dQc) transportLayer(k, pi, qc, qcLower, qcFlux, divQcFlux, dQc, true);
      }
    }
    if (part === 'tracers') return;
    edgePi(pi);
    for (let i = 0; i < C; i++) lnPi[i] = Math.log(pi[i]);
    gradient(mesh, lnPi, gradLnPi);
    for (let k = kFrom; k < kTo; k++) {
      const off = k * C;
      const fk = flux.subarray(k * E, (k + 1) * E);
      const uk = u.subarray(k * E, (k + 1) * E);
      const dUk = dU.subarray(k * E, (k + 1) * E);
      curl(mesh, uk, zeta);
      for (let v = 0; v < V; v++) qVertex[v] = (zeta[v] + fVertex[v]) / piVertex[v];
      for (let e = 0; e < E; e++) qEdge[e] = 0.5 * (qVertex[verticesOnEdge[2 * e]] + qVertex[verticesOnEdge[2 * e + 1]]);
      kineticEnergy(mesh, uk, kinetic);
      for (let i = 0; i < C; i++) phi[i] = geopotential[off + i] + kinetic[i];
      gradient(mesh, phi, gradPhi);
      for (let e = 0; e < E; e++) {
        const i = cellsOnEdge[2 * e], j = cellsOnEdge[2 * e + 1];
        let pv = 0;
        const qHere = 0.5 * qEdge[e];
        for (let s = 0; s < nEdgesOnEdge[e]; s++) {
          const slot = maxEdgesOnEdge * e + s, other = edgesOnEdge[slot];
          pv += pvWeights[slot] * fk[other] * (qHere + 0.5 * qEdge[other]);
        }
        const pgfPi = R * 0.5 * (thetaV[off + i] * exnerLayer[off + i] + thetaV[off + j] * exnerLayer[off + j]) * gradLnPi[e];
        const lowerFlow = 0.5 * (piSigmaDot[(k + 1) * C + i] + piSigmaDot[(k + 1) * C + j]);
        const upperFlow = 0.5 * (piSigmaDot[k * C + i] + piSigmaDot[k * C + j]);
        const lowerU = k === K - 1 ? 0 : 0.5 * (uk[e] + u[(k + 1) * E + e]);
        const upperU = k === 0 ? 0 : 0.5 * (uk[e] + u[(k - 1) * E + e]);
        const vertical = (lowerFlow * lowerU - upperFlow * upperU - uk[e] * (lowerFlow - upperFlow)) / (piEdge[e] * dSigma[k]);
        dUk[e] = pv / dcEdge[e] - gradPhi[e] - pgfPi - vertical;
      }
      if (nu4 > 0 && !splitClosure) momentumClosure(k, u, dU, nu4);
    }
  }

  function tendency(state, out) {
    phaseFlux(state, 0, K);
    phaseColumn(state, out, 0, C);
    phaseVertex(state, 0, V);
    phaseLayer(state, out, 0, K);
  }

  const diagnostics = { K, C, E, V, levels, sigmaMid, sigmaLower, sigmaUpper, dSigma, kappa, cp, R, g, p0, exnerLayer, exnerLower, dExnerDpi, geopotential, piSigmaDot, diagnose, diagnoseColumn };

  function mass(pi) {
    let m = 0;
    for (let i = 0; i < C; i++) m += mesh.areaCell[i] * pi[i];
    return m / g;
  }

  return { K, levels, sigmaMid, nu4, nu4Theta, divergenceDamping, divergenceSurfaces, coupledClosure: coupled, spongeRates, spongeMeanRates, spacing, tendency, phaseFlux, phaseColumn, phaseVertex, phaseLayer, phaseClosure, splitClosure, diagnose, diagnoseColumn, diagnostics, mass, shared, arrays: { exnerLayer, exnerLower, dExnerDpi, geopotential, piSigmaDot, thetaLower, qLower, qcLower, thetaV, dissipation } };
}
