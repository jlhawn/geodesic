import { LATENT_HEAT, EPSILON, saturationHumidity, liftingCondensationLevel } from './moist.module.js';
import { createMixedLayer, dycomsLongwave } from './mixedLayer.module.js';
export const STEFAN_BOLTZMANN = 5.670374419e-8;
export const SOLAR_CONSTANT = 1362;
export const AXIAL_TILT = 23.44 * Math.PI / 180;
export const DAY = 86400;
export const YEAR = 365 * DAY;

/*
 * Unit vector toward the sun at model time t. t = 0 is the spring equinox
 * with the sun over the +x meridian; the subsolar latitude follows
 * AXIAL_TILT·sin(2πt/YEAR) and the sun circles westward once per day.
 */
export function sunDirection(t, out = new Float64Array(3)) {
  const tilt = -AXIAL_TILT * Math.sin((t % YEAR) * 2 * Math.PI / YEAR);
  const x = Math.cos(tilt), z = -Math.sin(tilt);
  const spin = -((t % DAY) * 2 * Math.PI / DAY);
  out[0] = x * Math.cos(spin);
  out[1] = x * Math.sin(spin);
  out[2] = z;
  return out;
}

/*
 * Three-band gray longwave column with a slab-ocean surface and a bulk
 * sensible heat flux. A window band carrying the fraction `window` of
 * blackbody emission is transparent: the surface radiates it straight
 * to space. A vapour band whose optical depth follows Frierson et al.
 * (2006) — tau0(lat) = tauEquator + (tauPole − tauEquator) sin²lat,
 * distributed in the vertical as tau0 (f·σ + (1 − f)·σ⁴) — concentrates
 * near the surface like water vapour. A well-mixed-gas band carrying
 * `gasFraction` of the emission has the optical depth gasOpticalDepth
 * spread uniformly per unit mass, so thin high layers keep an
 * emissivity they can cool with, as CO₂'s 15 µm band lets the
 * stratosphere do. In each absorbing band every emission is either
 * absorbed on its way or leaves through the top or reaches the surface,
 * so the layer and surface energy fluxes sum exactly to absorbed solar
 * minus outgoing longwave. With vaporCoupling > 0 (m²/kg) and a
 * humidity field, the vapour band's optical depth is instead
 * vaporCoupling times each layer's water mass, so the greenhouse
 * follows the model's own humidity.
 *
 * Clouds: each layer's cloud water path gives it a gray emissivity
 * 1 − exp(−cloudAbsorption × path) that joins every longwave band —
 * including the window, which is transparent only where there is no
 * cloud. In the shortwave the column's cloud optical depth
 * cloudScattering × path reflects the beam by the two-stream
 * reflectance τ / (τ + 2μ). What reaches the surface is direct beam,
 * exp(−τ/μ) of it less the clear-sky `skylight` fraction, and diffuse
 * light, the rest; the surface reflects each with its own albedo, and
 * the multiple reflections between surface and cloud base (diffuse,
 * at the mean cosine DIFFUSE_MU) are summed. The two albedos are given
 * per cell (open water or sea ice); `albedo` is the default for both.
 * Cloud water also absorbs: a cloud of path W passes on only the
 * fraction exp(−cloudSolarAbsorption × W) of what it reflects or
 * transmits, both of the beam from above and of the light the surface
 * sends back up through it — with the default 0.4 m²/kg it absorbs
 * 3.9 % at 100 g/m² and 15 % at 400 g/m² (Stephens 1978) — and what it
 * absorbs heats the cloudy layers in proportion to their water, the
 * deck's layer counting the deck's water in the overcast column.
 *
 * With cloudCover 'pdf' (the default) a layer's resolved cloud covers
 * the fraction f of the layer that a uniform distribution of its total
 * water q + qc, of half-width (1 − RHc) qs about the mean, holds above
 * saturation, clipped to [coverFloor, 1]: RHc is
 * `boundaryCriticalHumidity` (0.85) for a layer whose midpoint lies below
 * the boundary-layer top and `criticalHumidity` (0.8) above it, so a
 * just-saturated layer is half covered and one holding (1 − RHc) qs of
 * cloud water overcast. Cloud under a strong inversion is stratiform and
 * uniform (the EIS cover of Wood & Bretherton reaches 1 near 11 K), so
 * where the column's estimated inversion strength, the deck's EIS of the
 * lowest layer's air, rises through `overcastInversion` (8 to 12 K) f
 * blends linearly into the f of a distribution whose half-width is also
 * at most the layer's cloud water but not below `overcastWater`
 * (5·10⁻⁵ kg/kg; null: no bound): there a saturated layer holding more
 * cloud water than that is overcast. Each cell keeps that share of the
 * ramp in `stratiform`, whatever the cover, for the boundary layer. Its emissivity is
 * f (1 − exp(−cloudAbsorption × path / f)), and the shortwave is the
 * blend, at the column's cover f̄, of the clear column and the column
 * whose cloud path lies in f̄, as the deck below blends its two columns.
 * Each layer is seen through its visibility 1 − exp(−path /
 * VISIBLE_PATH), 1 g/m². With `cloudOverlap` 'maximumRandom' (the
 * default) the layers of each run of adjacent cloudy layers overlap
 * maximally and the runs randomly: f̄ is 1 − Π(1 − f_run), f_run the
 * largest of its layers' f times their visibility; 'maximum' overlaps
 * every layer maximally, f̄ the largest over the column.
 * 'overcast' gives every cloudy layer the whole cell. With `cumulusCloud`
 * and the moist physics' shallow cumulus (`useCumulus`) a plume layer adds
 * its cumulus fraction times the plume's condensate to its water and
 * covers the larger of that fraction and its resolved cloud's f.
 *
 * Marine stratocumulus: over the part of a cell that is ice-free sea
 * (`openSea`, the per-cell fraction the caller passes; no deck without it)
 * a diagnostic deck covers the fraction f of the column. By default f
 * and the deck's water come from the mixed-layer model (below); with
 * mixedLayerDeck: false from an empirical fit, in which f is
 * 0.19 + 0.08 (EIS − 1) clamped to [0, 1] — 0.2 at the warm pool's EIS
 * of about 1 K, 0.67 at the south-east Pacific deck's 7 K, the 6–8 %
 * per K of Wood & Bretherton (2006) — times a ramp from 0 at a 5 °C
 * surface to 1 at 10 °C that keeps the deck off polar seas, times
 * `openSea`. EIS is their estimated inversion strength
 * LTS − Γ_θ (z_700 − z_LCL): the lower-tropospheric stability LTS is the
 * potential temperature of the layer nearest σ = 0.7 less the lowest
 * layer's, z_700 that layer's height above the lowest layer, z_LCL the
 * height of the lifting condensation level of the lowest layer's air
 * (Bolton's, as the moist physics finds it, reached along the dry
 * adiabat), and Γ_θ = g/c_p − Γ_m the potential-temperature gradient of
 * the moist adiabat at 850 hPa and the mean temperature of the two
 * layers, so EIS counts only the θ at σ = 0.7 beyond what a moist
 * adiabat from cloud base reaches. With stratusIndex 'ectei' the fit
 * takes instead the estimated cloud-top entrainment index of Kawai,
 * Koshiro & Webb (2017), ECTEI = EIS − 0.23 (L/c_p)(q_lowest − q_700)
 * with q_700 the humidity of the σ = 0.7 layer, which lowers the cover
 * where the air the deck entrains is dry. The deck fills the boundary
 * layer from z_LCL to the boundary-layer top, `mixedDepth` metres above
 * the lowest layer, and its water path is stratusScale × ½ Γ_l Δz² for
 * that thickness Δz, at most stratusWaterMax, with Γ_l the adiabatic
 * liquid-water lapse rate at cloud base. The water sits in the layer
 * nearest σ = stratusSigma and never enters qc. The deck and the clear
 * part of the column are two independent columns: every shortwave
 * quantity is the f-weighted mean of the column with and without the
 * deck's water, and the deck layer's cloud emissivity is the f-weighted
 * mean of its emissivity with and without it. `stratus: false` removes
 * the deck.
 *
 * The mixed-layer deck (mixedLayerDeck, the default) takes its cover and
 * water path from the mixed-layer model (mixedLayer.module.js, options
 * `mixedLayer`), advanced one physics step in each column from this
 * start: h the inversion height carried per cell in mlmHeight
 * (prognosticHeight, the default), held by the model's `bound` to at
 * least the boundary-layer top above the surface and at most
 * maximumHeight and the column's inversion ceiling (below), or that top
 * itself where mlmHeight is unset (0) or with prognosticHeight: false;
 * θ_l and q_t the dσ-weighted means of the layers whose midpoints lie
 * below h, the free troposphere the first layer above it, the
 * subsidence w_s = −πσ̇/(ρ g) at h of the last dynamics stage, πσ̇
 * averaged with equal weights over the cell and its neighbours
 * subsidenceSmoothing times over (2; 0 to 2) at the two interfaces
 * bracketing h and interpolated between them, the surface fluxes this
 * column's bulk sensible heat and evaporation, and the longwave the
 * DYCOMS-II form (dycomsLongwave) driven by the mixed layer's own liquid
 * water. The step's h, bounded the same way,
 * goes back into mlmHeight and the cover and water path are diagnosed
 * there; where the deck does not run, mlmHeight relaxes
 * toward the boundary-layer top with the model's `relax`
 * (heightMemory, 1 day). With deckRest 'inversion' that resting height,
 * and the start of an unset one, is the inversion ceiling (below) where
 * the column has one: under cumulus the Richardson depth is the top of
 * the subcloud layer, well below the inversion that caps the cloud
 * layer, and a deck started there finds no jump. Only h is carried: θ_l and q_t are the
 * column's again at each step, and the deck acts on the column through
 * its radiation and through the boundary layer's mixing depth, mlmTop —
 * the deck's h in the boundary layer's height coordinate (that of its
 * `depth`) where the deck ran with a carried height, 0 elsewhere (see
 * boundaryLayer.module.js).
 * The inversion ceiling keeps the deck under the column's own inversion,
 * the lowest interface whose upper layer's midpoint lies above the
 * boundary-layer top (and whose lower one's below maximumHeight) across
 * which θ_v rises by minimumInversion: the ceiling is 1 m below the
 * midpoint of the layer above that interface, so that layer stays the
 * free troposphere the deck entrains. Above the lowest kilometre the
 * layers are thick enough that the free troposphere's own
 * stratification across one of them passes the 2 K test below, and a
 * deck entraining under such a weak jump would deepen into it. A column
 * with no such interface has no ceiling but maximumHeight.
 * The regime test is the capping inversion: Δθ_v ≥ minimumInversion
 * (2 K, a capping inversion rather than the top of a subcloud layer
 * under cumulus) at the start's h. A stratocumulus-topped layer also
 * needs large-scale subsidence (DYCOMS-II has 3 mm/s at 840 m), but the
 * model resolves it as a residual of a few mm/s: in the SE Pacific and
 * Peru boxes the running mean w̄_s below sinks at 1.8 and 2.1 mm/s with
 * a spread over the cells as large, so a floor on the sink turns columns
 * of the regime away on synoptic swings, while the inversion is the
 * resolved record of the subsidence that built it. The subsidence test
 * therefore only vetoes large-scale ascent, under which an inversion at
 * the boundary-layer top is transient: it passes where
 * w̄_s ≤ −stratusSubsidence, and stratusSubsidence is −1 mm/s, about
 * twice the grid-scale residual of w̄_s, so the deck is refused where
 * the mean rises faster than 1 mm/s.
 * The divergent computational mode of the hexagonal C-grid puts most of
 * the variance of one stage's πσ̇ at the neighbouring-cell scale: on the
 * day-183 N=128 state the sink at h spreads over the SE Pacific's cells
 * by 30 mm/s about a mean of 1.5, 97 % of it at that scale, and the two
 * ring passes leave 8.6 mm/s about 1.6, 12 % of it there. w̄_s is the
 * running mean of w_s, w̄_s ← w̄_s e + w_s (1 − e) with
 * e = exp(−dt/subsidenceMemory) (2 days), because the large-scale
 * subsidence that defines the regime is still a small residual of the
 * synoptic swings in w_s; it is kept per cell in mlmSubsidence, starts
 * at 0, and is saved with the state (the key `mlmSubsidence`; a state
 * saved without it starts from 0). The model's own dh/dt keeps the
 * step's w_s.
 * The deck runs where the gates have mostly passed of late: mlmGate
 * holds the running mean of the pass indicator P (1 or 0),
 * G ← G + (P − G)(1 − exp(−dt/gateMemory)) (1 day), and the deck runs
 * while G > 0.5, or G = 0.5 and P = 1; elsewhere the column has no
 * deck. A standing deck so outlives its gates by ln 2 × gateMemory
 * (17 h), a new one waits as long for them, and a column whose Δθ_v
 * hovers about 2 K keeps the state its recent majority gives it. G
 * starts at 0.5, undecided, where the first step's gates decide;
 * gateMemory: 0 makes G the instantaneous P. mlmHeight and mlmGate are
 * saved with the state (keys `mlmHeight`, `mlmGate`; a state saved
 * without them starts from 0 and 0.5) and, like mlmSubsidence, keep
 * their values in cells where no deck is diagnosed. With
 * prognosticHeight: false and gateMemory: 0 the deck is re-diagnosed at
 * each step from the boundary-layer top behind the instantaneous gates.
 * With stratusSolar (the default) the mixed layer is also heated by the
 * sunlight its cloud absorbs: its forcing `absorbedSolar(W)` is what this
 * column absorbs in the deck's layer per unit deck area under a deck of
 * water path W — the overcast column's absorption in that layer less
 * the clear column's — so the mixed layer feels exactly the power the
 * column puts into that layer. Without it the mixed layer sees no
 * sunlight; the column's radiation is the same either way.
 * The deck covers the mixed layer's cover times `openSea`, with its water
 * path (at most stratusWaterMax) in the same layer and the same
 * two-column blend; the EIS is still diagnosed. The mixed layer's cover,
 * water path and entrainment rate are kept per cell in mlmCover,
 * mlmWater and mlmEntrainment (0 where it does not run).
 *
 * Γ_l: a saturated parcel conserves q_s + q_l, so it condenses −dq_s/dz
 * per metre of ascent. With d ln q_s/dT = L/(R_v T²), d ln q_s/d ln p
 * = −1, dp/dz = −p g/(R T) and the moist lapse rate
 * Γ_m = (g/c_p)(1 + L q_s/(R T))/(1 + L² q_s/(c_p R_v T²)),
 * −dq_s/dz = q_s (L Γ_m/(R_v T²) − g/(R T)), and Γ_l is that times the
 * air density p/(R T): 2.44e-6 kg/m³ per m at 290 K and 950 hPa.
 *
 * Shortwave: the fraction `ozoneAbsorption` of the incoming beam is
 * absorbed aloft. The ozone column follows Lacis & Hansen (1974)
 * (centred at ozoneHeight with width ozoneWidth, heights from σ with the
 * scale height) and the absorbing part of the beam decays through it
 * with the optical depth ozoneOpacity, so the heating peaks above the
 * ozone maximum as it does at the stratopause. Water vapour absorbs the
 * beam below it by the Lacis & Hansen (1974) absorptivity of the water
 * path the beam has crossed — pressure-scaled by √σ and lengthened by
 * their magnification 35/√(1224μ² + 1) — times `vaporAbsorption`, each
 * layer taking what its own vapour adds to the path above it; what is
 * left goes on to the clouds and the surface. Dry air absorbs nothing.
 */
export function waterVaporAbsorptivity(path) {
  return 2.9 * path / (Math.pow(1 + 141.5 * path, 0.635) + 5.925 * path);
}

const DIFFUSE_MU = 0.6;
export const STABILITY_SIGMA = 0.7;
export const DECK_CLOUD_LEVELS = 8;
export const UNDECIDED = 0.5;
export const VISIBLE_PATH = 1e-3;

/*
 * values[offset + i] averaged with equal weights over cell i and its
 * neighbours, `passes` times over: the value smoothCells in
 * levels.module.js gives cell i after that many passes over the field.
 */
export function ringMean(mesh, values, offset, i, passes) {
  if (passes <= 0) return values[offset + i];
  const { maxEdges, nEdgesOnCell, cellsOnCell } = mesh;
  let sum = ringMean(mesh, values, offset, i, passes - 1);
  for (let m = 0; m < nEdgesOnCell[i]; m++) sum += ringMean(mesh, values, offset, cellsOnCell[maxEdges * i + m], passes - 1);
  return sum / (nEdgesOnCell[i] + 1);
}

export function nearestLayer(sigmaMid, sigma) {
  let best = 0;
  for (let k = 1; k < sigmaMid.length; k++) if (Math.abs(sigmaMid[k] - sigma) < Math.abs(sigmaMid[best] - sigma)) best = k;
  return best;
}

export function stratusFraction(index, surfaceT) {
  return Math.min(1, Math.max(0, 0.19 + 0.08 * (index - 1))) * Math.min(1, Math.max(0, (surfaceT - 278.15) / 5));
}

function moistLapse(T, qs, cp, R, g, latentHeat) {
  return g / cp * (1 + latentHeat * qs / (R * T)) / (1 + latentHeat * latentHeat * qs / (cp * (R / EPSILON) * T * T));
}

export function inversionStrength(stability, lowerT, upperT, depth, cp, R, g, latentHeat = LATENT_HEAT) {
  const T = 0.5 * (lowerT + upperT);
  return stability - (g / cp - moistLapse(T, saturationHumidity(T, 85000), cp, R, g, latentHeat)) * depth;
}

export function entrainmentIndex(inversion, lowerQ, upperQ, cp, latentHeat = LATENT_HEAT) {
  return inversion - 0.23 * latentHeat / cp * (lowerQ - upperQ);
}

export function adiabaticWaterLapse(T, p, cp, R, g, latentHeat = LATENT_HEAT) {
  const qs = saturationHumidity(T, p), vaporR = R / EPSILON;
  const moist = moistLapse(T, qs, cp, R, g, latentHeat);
  return p / (R * T) * qs * (latentHeat * moist / (vaporR * T * T) - g / (R * T));
}

export function createRadiation(mesh, core, {
  solarConstant = SOLAR_CONSTANT, albedo = 0.07, cloudAbsorption = 130, cloudScattering = 95, stratus = true, stratusIndex = 'eis', stratusScale = 0.15, stratusWaterMax = 0.15, stratusSigma = 0.92,
  mixedLayerDeck = true, mixedLayer: mixedLayerOptions = {}, stratusSubsidence = -1e-3, minimumInversion = 2, subsidenceMemory = 2 * DAY, stratusSolar = true, cloudSolarAbsorption = 0.4,
  prognosticHeight = true, deckRest = 'depth', gateMemory = DAY, subsidenceSmoothing = 2, cloudCover = 'pdf', criticalHumidity = 0.8, boundaryCriticalHumidity = 0.85, coverFloor = 0.01, overcastWater = 5e-5, overcastInversion = [8, 12], cloudOverlap = 'maximumRandom',
  cumulusCloud = true, window = 0.25, tauEquator = 5.3, tauPole = 1.325, linearFraction = 0.1, gasFraction = 0.2, gasOpticalDepth = 7,
  ozoneAbsorption = 0.03, ozoneHeight = 25e3, ozoneWidth = 5e3, ozoneOpacity = 4, scaleHeight = 7e3, vaporAbsorption = 1,
  exchangeCoefficient = 1.5e-3, exchangeCoefficients = null, gustiness = 3, latentHeat = LATENT_HEAT, vaporCoupling = 0.55, skylight = 0.15, buffers = null,
} = {}) {
  const { K, C, dSigma, sigmaMid, cp, R, g, kappa, exnerLayer, exnerLower, geopotential, piSigmaDot, p0 } = core.diagnostics;
  const { thetaV } = core.arrays;
  const levels = core.levels;
  const vaporFraction = 1 - window - gasFraction;
  const opticalDepth = (lat) => tauEquator + (tauPole - tauEquator) * Math.sin(lat) ** 2;
  const tauCell = Float64Array.from({ length: C }, (_, i) => opticalDepth(mesh.latCell[i]));
  const shape = Float64Array.from({ length: K }, (_, k) => linearFraction * (levels[k + 1] - levels[k]) + (1 - linearFraction) * (levels[k + 1] ** 4 - levels[k] ** 4));
  const ozoneAbove = (sigma) => (sigma <= 0 ? 0 : (1 + Math.exp(-ozoneHeight / ozoneWidth)) / (1 + Math.exp((-scaleHeight * Math.log(sigma) - ozoneHeight) / ozoneWidth)));
  const beamLeft = (sigma) => Math.exp(-ozoneOpacity * ozoneAbove(sigma));
  const ozoneFraction = Float64Array.from({ length: K }, (_, k) => (beamLeft(levels[k]) - beamLeft(levels[k + 1])) / (1 - Math.exp(-ozoneOpacity)));
  const emissivity = new Float64Array(K);
  const cloudEmissivity = new Float64Array(K);
  const layerCover = new Float64Array(K).fill(1);
  let cumulusCover = null, cumulusWater = null;
  const clearSky = { absorbed: 0, down: 0, direct: 0, reflectance: 0, cloud: 0 };
  const vaporEmissivity = new Float64Array(K);
  const mixedEmissivity = new Float64Array(K);
  const surfaceFlux = new Float64Array(C), surfaceDirect = new Float64Array(C);
  const outgoingBuffer = buffers && buffers.outgoing ? buffers.outgoing : new SharedArrayBuffer(8 * C);
  const shortwaveBuffer = buffers && buffers.surfaceShortwave ? buffers.surfaceShortwave : new SharedArrayBuffer(8 * C);
  const outgoing = new Float64Array(outgoingBuffer), surfaceShortwave = new Float64Array(shortwaveBuffer);
  const evaporationBuffer = buffers && buffers.evaporation ? buffers.evaporation : new SharedArrayBuffer(8 * C);
  const evaporation = new Float64Array(evaporationBuffer);
  const stratusBuffer = buffers && buffers.stratus ? buffers.stratus : new SharedArrayBuffer(8 * C);
  const stratusPath = new Float64Array(stratusBuffer);
  const coverBuffer = buffers && buffers.stratusFraction ? buffers.stratusFraction : new SharedArrayBuffer(8 * C);
  const stratusCover = new Float64Array(coverBuffer);
  const indexBuffer = buffers && buffers.stabilityIndex ? buffers.stabilityIndex : new SharedArrayBuffer(8 * C);
  const stabilityIndex = new Float64Array(indexBuffer);
  const mlmCoverBuffer = buffers && buffers.mlmCover ? buffers.mlmCover : new SharedArrayBuffer(8 * C);
  const mlmWaterBuffer = buffers && buffers.mlmWater ? buffers.mlmWater : new SharedArrayBuffer(8 * C);
  const mlmEntrainmentBuffer = buffers && buffers.mlmEntrainment ? buffers.mlmEntrainment : new SharedArrayBuffer(8 * C);
  const mlmSubsidenceBuffer = buffers && buffers.mlmSubsidence ? buffers.mlmSubsidence : new SharedArrayBuffer(8 * C);
  const mlmHeightBuffer = buffers && buffers.mlmHeight ? buffers.mlmHeight : new SharedArrayBuffer(8 * C);
  const mlmGateBuffer = buffers && buffers.mlmGate ? buffers.mlmGate : new SharedArrayBuffer(8 * C);
  const mlmTopBuffer = buffers && buffers.mlmTop ? buffers.mlmTop : new SharedArrayBuffer(8 * C);
  const mlmCover = new Float64Array(mlmCoverBuffer), mlmWater = new Float64Array(mlmWaterBuffer), mlmEntrainment = new Float64Array(mlmEntrainmentBuffer), mlmSubsidence = new Float64Array(mlmSubsidenceBuffer);
  const mlmHeight = new Float64Array(mlmHeightBuffer), mlmGate = new Float64Array(mlmGateBuffer), mlmTop = new Float64Array(mlmTopBuffer);
  const stratiformBuffer = buffers && buffers.stratiform ? buffers.stratiform : new SharedArrayBuffer(8 * C);
  const stratiformShare = new Float64Array(stratiformBuffer);
  if (!(buffers && buffers.mlmGate)) mlmGate.fill(UNDECIDED);
  const shadow = mixedLayerDeck ? createMixedLayer({ cp, R, g, latentHeat, referencePressure: p0, cloudLevels: DECK_CLOUD_LEVELS, ...mixedLayerOptions }) : null;
  const shadowLongwave = dycomsLongwave();
  if (stratusIndex !== 'eis' && stratusIndex !== 'ectei') throw new Error(`stratusIndex must be 'eis' or 'ectei', not ${stratusIndex}`);
  if (cloudOverlap !== 'maximum' && cloudOverlap !== 'maximumRandom') throw new Error(`cloudOverlap must be 'maximum' or 'maximumRandom', not ${cloudOverlap}`);
  if (!(overcastInversion?.[1] > overcastInversion?.[0])) throw new Error(`overcastInversion must rise from its first to its second EIS, not ${overcastInversion}`);
  if (![0, 1, 2].includes(subsidenceSmoothing)) throw new Error(`subsidenceSmoothing must be 0, 1 or 2, not ${subsidenceSmoothing}`);
  if (deckRest !== 'depth' && deckRest !== 'inversion') throw new Error(`deckRest must be 'depth' or 'inversion', not ${deckRest}`);
  const entraining = stratusIndex === 'ectei';
  const stratusLayer = nearestLayer(sigmaMid, stratusSigma), stabilityLayer = nearestLayer(sigmaMid, STABILITY_SIGMA);
  const gasEmissivity = Float64Array.from({ length: K }, (_, k) => 1 - Math.exp(-gasOpticalDepth * (levels[k + 1] - levels[k])));
  const temperature = new Float64Array(K);
  const cloudWater = new Float64Array(K);
  const vaporTaken = new Float64Array(K);
  const emitted = new Float64Array(K);
  const netFlux = new Float64Array(K);
  const sun = new Float64Array([1, 0, 0]);
  const budget = { absorbedSolar: 0, atmosphereSolar: 0, outgoingLongwave: 0, sensibleHeat: 0, evaporation: 0, surfaceFlux: 0, insolation: 0, reflectedSolar: 0, cloudReflectance: 0, cloudCover: 0, cloudSolar: 0, stratus: 0, stratusFraction: 0, stabilityIndex: NaN, mlmCover: 0, mlmWater: 0, mlmEntrainment: 0, mlmSolar: 0, mlmTop: 0, stratiform: 0 };
  const sky = { absorbed: 0, down: 0, direct: 0, reflectance: 0, cloud: 0 }, decked = { absorbed: 0, down: 0, direct: 0, reflectance: 0, cloud: 0 }, probe = { absorbed: 0, down: 0, direct: 0, reflectance: 0, cloud: 0 };
  const deckLight = { incident: 0, mu: 0, direct: 0, diffuse: 0, path: 0, layer: 0, clear: 0 };

  function setTime(t) {
    sunDirection(t, sun);
  }

  function cosZenith(i) {
    return Math.max(0, mesh.xCell[3 * i] * sun[0] + mesh.xCell[3 * i + 1] * sun[1] + mesh.xCell[3 * i + 2] * sun[2]);
  }

  function insolation(i) {
    return solarConstant * cosZenith(i);
  }

  function band(fraction, eps, surfaceEmission) {
    for (let k = 0; k < K; k++) emitted[k] = fraction * eps[k] * STEFAN_BOLTZMANN * temperature[k] ** 4;
    let down = 0;
    for (let k = 0; k < K; k++) {
      netFlux[k] += eps[k] * down - 2 * emitted[k];
      down = down * (1 - eps[k]) + emitted[k];
    }
    let up = fraction * surfaceEmission;
    for (let k = K - 1; k >= 0; k--) {
      netFlux[k] += eps[k] * up;
      up = up * (1 - eps[k]) + emitted[k];
    }
    return [up, down];
  }

  function deckWater(lcl, base, mixedDepth) {
    const thickness = mixedDepth - base;
    if (thickness <= 0) return 0;
    return Math.min(stratusWaterMax, stratusScale * 0.5 * adiabaticWaterLapse(lcl.temperature, lcl.pressure, cp, R, g, latentHeat) * thickness * thickness);
  }

  function inversionCeiling(i, surface, floor) {
    for (let k = K - 2; k >= 1; k--) {
      const upper = (geopotential[k * C + i] - surface) / g;
      if ((geopotential[(k + 1) * C + i] - surface) / g >= shadow.maximumHeight) break;
      if (upper > floor && thetaV[k * C + i] - thetaV[(k + 1) * C + i] >= minimumInversion) return upper - 1;
    }
    return Infinity;
  }

  function shadowDeck(i, pi, theta, q, qc, mixedDepth, sensible, evaporation, dt, absorbedSolar) {
    const bottom = (K - 1) * C + i;
    const surface = geopotential[bottom] - cp * thetaV[bottom] * (exnerLower[bottom] - exnerLayer[bottom]);
    const depth = mixedDepth + (geopotential[bottom] - surface) / g;
    const ceiling = prognosticHeight ? inversionCeiling(i, surface, depth) : Infinity;
    const resting = deckRest === 'inversion' && ceiling < shadow.maximumHeight ? ceiling : depth;
    const h = prognosticHeight && mlmHeight[i] > 0 ? shadow.bound(mlmHeight[i], depth, ceiling) : prognosticHeight ? resting : depth;
    const rest = () => { mlmHeight[i] = shadow.relax(mlmHeight[i], resting, dt); return false; };
    let weight = 0, heat = 0, water = 0, k = K - 1;
    for (; k >= 0 && geopotential[k * C + i] - surface < g * h; k--) {
      const idx = k * C + i, cloud = qc ? Math.max(0, qc[idx]) : 0;
      heat += dSigma[k] * (theta[idx] - latentHeat * cloud / (cp * exnerLayer[idx]));
      water += dSigma[k] * (Math.max(0, q[idx]) + cloud);
      weight += dSigma[k];
    }
    if (k < 1) return rest();
    const above = k * C + i, aboveCloud = qc ? Math.max(0, qc[above]) : 0;
    const interfaceHeight = (m) => (geopotential[m * C + i] + cp * thetaV[m * C + i] * (exnerLayer[m * C + i] - exnerLower[(m - 1) * C + i]) - surface) / g;
    let lowerHeight = 0, lower = K, m = K - 1;
    for (; m > k && interfaceHeight(m) < h; m--) { lowerHeight = interfaceHeight(m); lower = m; }
    const upperHeight = interfaceHeight(m);
    const lowerFlow = lower < K ? ringMean(mesh, piSigmaDot, lower * C, i, subsidenceSmoothing) : 0;
    const flow = lowerFlow + (ringMean(mesh, piSigmaDot, m * C, i, subsidenceSmoothing) - lowerFlow) * (h - lowerHeight) / (upperHeight - lowerHeight);
    const density = pi * sigmaMid[m] / (R * thetaV[m * C + i] * exnerLayer[m * C + i]);
    const subsidence = -flow / (density * g), keep = Math.exp(-dt / subsidenceMemory);
    mlmSubsidence[i] = mlmSubsidence[i] * keep + subsidence * (1 - keep);
    const sinking = !(mlmSubsidence[i] > -stratusSubsidence);
    const forcing = {
      surfacePressure: pi, sensibleHeat: sensible, evaporation, radiation: shadowLongwave, subsidence: () => subsidence, absorbedSolar,
      thetaLAbove: theta[above] - latentHeat * aboveCloud / (cp * exnerLayer[above]), qtAbove: Math.max(0, q[above]) + aboveCloud,
    };
    const start = { h, thetaL: heat / weight, qt: water / weight };
    let now = sinking ? shadow.diagnose(start, forcing) : null;
    const pass = sinking && now.virtualJump >= minimumInversion ? 1 : 0;
    const gate = gateMemory > 0 ? mlmGate[i] - (pass - mlmGate[i]) * Math.expm1(-dt / gateMemory) : pass;
    mlmGate[i] = gate;
    if (!(gate > UNDECIDED || (gate === UNDECIDED && pass === 1))) return rest();
    now ??= shadow.diagnose(start, forcing);
    const stepped = dt > 0 ? shadow.step(start, forcing, dt, now) : start;
    const top = prognosticHeight ? shadow.bound(stepped.h, depth, ceiling) : stepped.h;
    const next = dt > 0 ? shadow.diagnose(top === stepped.h ? stepped : { h: top, thetaL: stepped.thetaL, qt: stepped.qt }, forcing) : now;
    if (!(Number.isFinite(next.liquidWaterPath) && Number.isFinite(next.cover) && Number.isFinite(now.entrainment))) return rest();
    mlmHeight[i] = top;
    budget.mlmCover = next.cover;
    budget.mlmWater = next.liquidWaterPath;
    budget.mlmEntrainment = now.entrainment;
    budget.mlmSolar = next.absorbedShortwave;
    budget.mlmTop = prognosticHeight ? top + surface / g : 0;
    return true;
  }

  function shortwave(out, cloudDepth, keep, mu, surfaceAlbedo, diffuseAlbedo) {
    const reflectance = mu > 0 && cloudDepth > 0 ? cloudDepth / (cloudDepth + 2 * mu) : 0;
    let direct = (1 - skylight) * (cloudDepth > 0 && mu > 0 ? Math.exp(-cloudDepth / mu) : 1);
    let diffuse = 1 - reflectance - direct;
    let returned = cloudDepth > 0 ? cloudDepth / (cloudDepth + 2 * DIFFUSE_MU) : 0;
    direct *= keep; diffuse *= keep; returned *= keep;
    const upward = surfaceAlbedo * direct + diffuseAlbedo * diffuse;
    const reflections = returned * upward / (1 - diffuseAlbedo * returned);
    out.absorbed = (1 - surfaceAlbedo) * direct + (1 - diffuseAlbedo) * (diffuse + reflections);
    out.down = direct + diffuse + reflections;
    out.direct = direct;
    out.reflectance = keep * reflectance;
    out.cloud = (1 - keep) * (1 + upward / (1 - diffuseAlbedo * returned));
  }

  function deckAbsorption(liquidWaterPath) {
    const water = Math.min(stratusWaterMax, liquidWaterPath);
    if (!(water > 0) || !(deckLight.incident > 0)) return 0;
    const total = deckLight.path + water;
    shortwave(probe, cloudScattering * total, Math.exp(-cloudSolarAbsorption * total), deckLight.mu, deckLight.direct, deckLight.diffuse);
    return deckLight.incident * probe.cloud / total * (deckLight.layer + water) - deckLight.clear;
  }

  function column(i, pi, theta, surfaceT, windSpeed, tau0 = tauCell[i], beam = insolation(i), qAir = null, q = null, qc = null, surfaceAlbedo = albedo, diffuseAlbedo = surfaceAlbedo, wetness = 1, exchangeCoefficientAt = exchangeCoefficient, openSea = 0, mixedDepth = 0, dt = 0) {
    const ozoneHeating = beam * ozoneAbsorption;
    const surfaceEmission = STEFAN_BOLTZMANN * surfaceT * surfaceT * surfaceT * surfaceT;
    const coupled = vaporCoupling > 0 && q !== null;
    const bottom = K - 1;
    const airTemperature = theta[bottom * C + i] * exnerLayer[bottom * C + i];
    const airDensity = pi * sigmaMid[bottom] / (R * airTemperature);
    const exchange = airDensity * exchangeCoefficientAt * Math.max(windSpeed, gustiness);
    const sensible = exchange * cp * (surfaceT - airTemperature);
    const evaporation = qAir === null ? 0 : wetness * Math.max(0, exchange * (saturationHumidity(surfaceT, pi) - qAir));
    const mu = beam / solarConstant, lit = vaporAbsorption > 0 && q !== null && mu > 0;
    if (lit) {
      const magnification = 35 / Math.sqrt(1224 * mu * mu + 1);
      let path = 0;
      for (let k = 0; k < K; k++) {
        path += Math.max(0, q[k * C + i]) * pi * dSigma[k] / g * Math.sqrt(sigmaMid[k]) * 0.1 * magnification;
        vaporTaken[k] = vaporAbsorption * waterVaporAbsorptivity(path);
      }
    }
    let cloudPath = 0, columnCover = 0, block = 0, clearColumn = 1, inversionShare = 0;
    if (qAir !== null) {
      const lcl = liftingCondensationLevel(airTemperature, qAir, pi * sigmaMid[bottom], kappa);
      if (lcl) {
        const lower = bottom * C + i, upper = stabilityLayer * C + i, base = Math.max(0, cp * (airTemperature - lcl.temperature) / g);
        const inversion = inversionStrength(theta[upper] - theta[lower], airTemperature, theta[upper] * exnerLayer[upper], (geopotential[upper] - geopotential[lower]) / g - base, cp, R, g, latentHeat);
        inversionShare = Math.min(1, Math.max(0, (inversion - overcastInversion[0]) / (overcastInversion[1] - overcastInversion[0])));
      }
    }
    budget.stratiform = inversionShare;
    const stratiform = cloudCover === 'pdf' && overcastWater !== null ? inversionShare : 0;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      cloudWater[k] = qc ? Math.max(0, qc[idx]) * (pi * dSigma[k] / g) : 0;
      cloudPath += cloudWater[k];
      layerCover[k] = 1;
      if (cloudCover === 'pdf' && q !== null && cloudWater[k] > 0) {
        const inside = (geopotential[idx] - geopotential[bottom * C + i]) / g < mixedDepth;
        const qs = saturationHumidity(theta[idx] * exnerLayer[idx], pi * sigmaMid[k]), water = Math.max(0, qc[idx]), excess = Math.max(0, q[idx]) + water - qs;
        const width = (1 - (inside ? boundaryCriticalHumidity : criticalHumidity)) * qs;
        layerCover[k] = Math.min(1, Math.max(coverFloor, (excess + width) / (2 * width)));
        if (stratiform > 0) {
          const bound = Math.min(width, Math.max(water, overcastWater));
          layerCover[k] = (1 - stratiform) * layerCover[k] + stratiform * Math.min(1, Math.max(coverFloor, (excess + bound) / (2 * bound)));
        }
      }
      const cumulus = cumulusCover !== null ? cumulusCover[idx] * cumulusWater[idx] * (pi * dSigma[k] / g) : 0;
      if (cumulus > 0) {
        layerCover[k] = cloudWater[k] > 0 ? Math.max(layerCover[k], cumulusCover[idx]) : cumulusCover[idx];
        cloudWater[k] += cumulus;
        cloudPath += cumulus;
      }
      const seen = cloudCover === 'pdf' && cloudWater[k] > 0 ? layerCover[k] * -Math.expm1(-cloudWater[k] / VISIBLE_PATH) : 0;
      if (seen > 0) block = Math.max(block, seen);
      if (block > 0 && (!(seen > 0) || k === K - 1)) { clearColumn *= 1 - block; columnCover = Math.max(columnCover, block); block = 0; }
    }
    if (cloudOverlap === 'maximumRandom' && cloudCover === 'pdf') columnCover = 1 - clearColumn;
    if (!(columnCover > 0)) columnCover = 1;
    const sunlit = beam - ozoneHeating, vaporHeating = lit ? sunlit * vaporTaken[K - 1] : 0, incident = sunlit - vaporHeating;
    const inCloud = cloudPath / columnCover;
    shortwave(sky, cloudScattering * inCloud, Math.exp(-cloudSolarAbsorption * inCloud), mu, surfaceAlbedo, diffuseAlbedo);
    if (columnCover < 1) {
      shortwave(clearSky, 0, 1, mu, surfaceAlbedo, diffuseAlbedo);
      for (const key of ['absorbed', 'down', 'direct', 'reflectance', 'cloud']) sky[key] = columnCover * sky[key] + (1 - columnCover) * clearSky[key];
    }
    let clearShare = cloudPath > 0 ? incident * sky.cloud / cloudPath : 0, deckShare = 0;
    let fraction = 0, deck = 0, index = NaN;
    budget.mlmCover = 0; budget.mlmWater = 0; budget.mlmEntrainment = 0; budget.mlmSolar = 0; budget.mlmTop = 0;
    if (stratus && openSea > 0 && qAir !== null && mixedDepth > 0) {
      const lower = bottom * C + i, upper = stabilityLayer * C + i, lowerT = theta[lower] * exnerLayer[lower];
      const lcl = liftingCondensationLevel(lowerT, qAir, pi * sigmaMid[bottom], kappa);
      if (lcl) {
        const base = Math.max(0, cp * (lowerT - lcl.temperature) / g);
        const inversion = inversionStrength(theta[upper] - theta[lower], lowerT, theta[upper] * exnerLayer[upper], (geopotential[upper] - geopotential[lower]) / g - base, cp, R, g, latentHeat);
        index = entraining ? entrainmentIndex(inversion, qAir, q ? q[upper] : qAir, cp, latentHeat) : inversion;
        if (!shadow) {
          fraction = stratusFraction(index, surfaceT) * openSea;
          if (fraction > 0) deck = deckWater(lcl, base, mixedDepth);
        }
      }
      if (shadow && q) {
        deckLight.incident = incident; deckLight.mu = mu; deckLight.direct = surfaceAlbedo; deckLight.diffuse = diffuseAlbedo;
        deckLight.path = cloudPath; deckLight.layer = cloudWater[stratusLayer]; deckLight.clear = clearShare * cloudWater[stratusLayer];
        if (shadowDeck(i, pi, theta, q, qc, mixedDepth, sensible, evaporation, dt, stratusSolar ? deckAbsorption : null)) {
          fraction = budget.mlmCover * openSea;
          if (fraction > 0) deck = Math.min(stratusWaterMax, budget.mlmWater);
        }
      }
      if (deck <= 0) fraction = 0;
    }
    for (let k = 0; k < K; k++) {
      const mass = pi * dSigma[k] / g;
      emissivity[k] = 1 - Math.exp(coupled ? -vaporCoupling * Math.max(0, q[k * C + i]) * mass : -tau0 * shape[k]);
      const water = cloudWater[k];
      cloudEmissivity[k] = water > 0 ? layerCover[k] * (1 - Math.exp(-cloudAbsorption * water / layerCover[k])) : 0;
      if (deck > 0 && k === stratusLayer) cloudEmissivity[k] = fraction * (1 - Math.exp(-cloudAbsorption * (water + deck))) + (1 - fraction) * cloudEmissivity[k];
      const clear = 1 - cloudEmissivity[k];
      vaporEmissivity[k] = 1 - (1 - emissivity[k]) * clear;
      mixedEmissivity[k] = 1 - (1 - gasEmissivity[k]) * clear;
      temperature[k] = theta[k * C + i] * exnerLayer[k * C + i];
      netFlux[k] = ozoneHeating * ozoneFraction[k];
    }
    if (lit) {
      let taken = 0;
      for (let k = 0; k < K; k++) {
        const through = vaporTaken[k];
        netFlux[k] += sunlit * (through - taken);
        taken = through;
      }
    }
    if (deck > 0) {
      const total = cloudPath + deck;
      shortwave(decked, cloudScattering * total, Math.exp(-cloudSolarAbsorption * total), mu, surfaceAlbedo, diffuseAlbedo);
      deckShare = fraction * incident * decked.cloud / total;
      clearShare *= 1 - fraction;
      sky.absorbed = fraction * decked.absorbed + (1 - fraction) * sky.absorbed;
      sky.down = fraction * decked.down + (1 - fraction) * sky.down;
      sky.direct = fraction * decked.direct + (1 - fraction) * sky.direct;
      sky.reflectance = fraction * decked.reflectance + (1 - fraction) * sky.reflectance;
      sky.cloud = fraction * decked.cloud + (1 - fraction) * sky.cloud;
    }
    let cloudHeating = 0;
    if (incident > 0 && sky.cloud > 0) {
      for (let k = 0; k < K; k++) {
        const share = clearShare * cloudWater[k] + deckShare * (k === stratusLayer ? cloudWater[k] + deck : cloudWater[k]);
        netFlux[k] += share;
        cloudHeating += share;
      }
    }
    const absorbedSolar = incident * sky.absorbed;
    const [outVapor, backVapor] = band(vaporFraction, vaporEmissivity, surfaceEmission);
    const [outGas, backGas] = band(gasFraction, mixedEmissivity, surfaceEmission);
    const [outWindow, backWindow] = band(window, cloudEmissivity, surfaceEmission);
    const outgoing = outVapor + outGas + outWindow;
    const back = backVapor + backGas + backWindow;
    netFlux[bottom] += sensible;
    const net = absorbedSolar - surfaceEmission + back - sensible - latentHeat * evaporation;
    budget.absorbedSolar = absorbedSolar + ozoneHeating + vaporHeating + cloudHeating;
    budget.atmosphereSolar = ozoneHeating + vaporHeating + cloudHeating;
    budget.cloudSolar = cloudHeating;
    budget.outgoingLongwave = outgoing;
    budget.sensibleHeat = sensible;
    budget.evaporation = evaporation;
    budget.surfaceFlux = net;
    budget.insolation = beam;
    budget.reflectedSolar = incident - absorbedSolar - cloudHeating;
    budget.surfaceShortwave = incident * sky.down;
    budget.surfaceDirect = incident * sky.direct;
    budget.cloudReflectance = sky.reflectance;
    budget.cloudCover = cloudPath > 0 ? columnCover : 0;
    budget.stratus = deck;
    budget.stratusFraction = fraction;
    budget.stabilityIndex = index;
    return net;
  }

  /*
   * Heating tendencies of the layers and the evaporation tendency of the
   * lowest layer for the cells in range; the net surface flux of each
   * cell is left in `surfaceFlux` for the surface model to apply, with
   * the sunlight reaching the surface in `surfaceShortwave`, of which
   * `surfaceDirect` is the direct beam, and the stratocumulus deck's
   * water path and cover in `stratus` and `stratusFraction`, with the
   * EIS or ECTEI its cover follows in `stabilityIndex` (NaN where the
   * deck is not diagnosed: over land or full ice, without humidity or a
   * boundary layer, or with `stratus: false`). `depth` is the height of
   * each cell's boundary-layer top; `dt` is the physics step the
   * mixed-layer deck advances by.
   */
  function apply(state, out, windSpeed, totals, iFrom = 0, iTo = C, surfaceAlbedo = null, diffuseAlbedo = null, wetness = null, openSea = null, depth = null, dt = 0) {
    const [pi, theta, , surfaceT] = state;
    const [, dTheta] = out;
    const q = state[4] ?? null, dQ = out[4] ?? null, qc = state[5] ?? null;
    const bottom = (K - 1) * C;
    if (totals) for (const name of ['absorbedSolar', 'atmosphereSolar', 'outgoingLongwave', 'sensibleHeat', 'evaporation', 'insolation', 'reflectedSolar']) totals[name] = 0;
    for (let i = iFrom; i < iTo; i++) {
      surfaceFlux[i] = column(i, pi[i], theta, surfaceT[i], windSpeed[i], tauCell[i], insolation(i), q && dQ ? q[bottom + i] : null, q && dQ ? q : null, q && dQ ? qc : null, surfaceAlbedo ? surfaceAlbedo[i] : albedo, diffuseAlbedo ? diffuseAlbedo[i] : surfaceAlbedo ? surfaceAlbedo[i] : albedo, wetness ? wetness[i] : 1, exchangeCoefficients ? exchangeCoefficients[i] : exchangeCoefficient, openSea ? openSea[i] : 0, depth ? depth[i] - geopotential[bottom + i] / g : 0, dt);
      outgoing[i] = budget.outgoingLongwave;
      stratusPath[i] = budget.stratus;
      stratusCover[i] = budget.stratusFraction;
      stabilityIndex[i] = budget.stabilityIndex;
      mlmCover[i] = budget.mlmCover;
      mlmWater[i] = budget.mlmWater;
      mlmEntrainment[i] = budget.mlmEntrainment;
      mlmTop[i] = budget.mlmTop;
      stratiformShare[i] = budget.stratiform;
      evaporation[i] = budget.evaporation;
      surfaceShortwave[i] = budget.surfaceShortwave;
      surfaceDirect[i] = budget.surfaceDirect;
      for (let k = 0; k < K; k++) {
        const massPerArea = pi[i] * dSigma[k] / g;
        dTheta[k * C + i] += netFlux[k] / (cp * massPerArea) / exnerLayer[k * C + i];
      }
      if (q && dQ) dQ[bottom + i] += budget.evaporation * g / (pi[i] * dSigma[K - 1]);
      if (totals) {
        const a = mesh.areaCell[i];
        totals.absorbedSolar += a * budget.absorbedSolar;
        totals.atmosphereSolar += a * budget.atmosphereSolar;
        totals.outgoingLongwave += a * budget.outgoingLongwave;
        totals.sensibleHeat += a * budget.sensibleHeat;
        totals.evaporation += a * budget.evaporation;
        totals.insolation += a * budget.insolation;
        totals.reflectedSolar += a * budget.reflectedSolar;
      }
    }
  }

  function useCumulus(cover, water) {
    const on = cumulusCloud && cloudCover === 'pdf' && cover && water;
    cumulusCover = on ? cover : null;
    cumulusWater = on ? water : null;
  }

  const deckGates = { subsidenceSmoothing, subsidenceMemory, stratusSubsidence, minimumInversion, gateMemory, deckRest };
  return { setTime, sun, cosZenith, insolation, column, apply, useCumulus, deckGates, layerFlux: netFlux, surfaceFlux, outgoing, surfaceShortwave, surfaceDirect, evaporation, stratus: stratusPath, stratusFraction: stratusCover, stabilityIndex, mlmCover, mlmWater, mlmEntrainment, mlmSubsidence, mlmHeight, mlmGate, mlmTop, stratiform: stratiformShare, stratusLayer, stabilityLayer, budget, emissivity, opticalDepth, ozoneFraction, shared: { outgoing: outgoingBuffer, surfaceShortwave: shortwaveBuffer, evaporation: evaporationBuffer, stratus: stratusBuffer, stratusFraction: coverBuffer, stabilityIndex: indexBuffer, mlmCover: mlmCoverBuffer, mlmWater: mlmWaterBuffer, mlmEntrainment: mlmEntrainmentBuffer, mlmSubsidence: mlmSubsidenceBuffer, mlmHeight: mlmHeightBuffer, mlmGate: mlmGateBuffer, mlmTop: mlmTopBuffer, stratiform: stratiformBuffer } };
}
