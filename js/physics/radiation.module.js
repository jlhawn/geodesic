import { LATENT_HEAT, EPSILON, R_VAPOR, saturationHumidity, liftingCondensationLevel, cloudSaturation, criticalHumidityAt, uniformCover } from './moist.module.js';
import { createMixedLayer, dycomsLongwave } from './mixedLayer.module.js';
import { REGIME } from './boundaryLayer.module.js';
import { SEA_DRAG } from './surface.module.js';
import { LONGWAVE_TABLE, LONGWAVE_CONSTANTS, GAS_MOLAR, layerPaths, planckShare } from './longwave.module.js';
import { ozoneWeights, ozoneAbove as climatologyAbove } from './ozone.module.js';
import { ozoneAbsorptivity, visibleVaporAbsorptivity, nearInfraredVaporAbsorptivity, oxygenAbsorptivity, carbonDioxideAbsorptivity, pressureScaling, vaporScaling, OXYGEN, OZONE_CM_ATM, STP_DEPTH, OZONE_SHARES, OZONE_COEFFICIENTS, VAPOR_STRENGTH } from './shortwaveGases.module.js';
export const STEFAN_BOLTZMANN = 5.670374419e-8;
export const SOLAR_CONSTANT = 1362;
export const AXIAL_TILT = 23.44 * Math.PI / 180;
export const DAY = 86400;
export const YEAR = 365 * DAY;
export const REFERENCE_RESISTANCE = 70;

/*
 * The FAO-56 Penman–Monteith reference evapotranspiration (Allen et al.
 * 1998) in kg/m²/s from the surface's net radiation (W/m²), the lowest
 * air's temperature, humidity and density at surface pressure p and the
 * aerodynamic conductance (m/s), through REFERENCE_RESISTANCE s/m.
 */
export function referenceEvaporation(netRadiation, airTemperature, q, p, density, conductance, cp, latentHeat) {
  const qs = saturationHumidity(airTemperature, p), slope = qs * 4302.645 / ((airTemperature - 29.65) * (airTemperature - 29.65));
  return (slope * netRadiation + density * cp * conductance * (qs - q)) / (latentHeat * slope + cp * (1 + REFERENCE_RESISTANCE * conductance));
}

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
 * Longwave column with a slab-ocean surface and a bulk sensible heat flux.
 * With longwaveScheme 'correlated' (the default) and a humidity field the
 * clear-sky gases are the g-points of js/physics/longwave.module.js
 * (scripts/longwaveFit.mjs): water vapour lines and self continuum, CO₂
 * (`carbonDioxide`), ozone, methane (`methane`) and nitrous oxide
 * (`nitrousOxide`), volume mixing ratios of dry air (by default
 * GREENHOUSE_GASES, NOAA GML's global annual means for 2010, the model's
 * present day), the ozone the column's (below); each g-point is a gray
 * band whose layer emissivity 1 − exp(−τ) joins the cloud's,
 * 1 − (1 − ε_gas)(1 − ε_cloud), emitting
 * its share of σT⁴. Without humidity, or with 'gray', the gases are the
 * three-band gray column: a window band carrying the fraction `window` of
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
 * Clouds: each layer's condensate path W gives it a gray emissivity
 * 1 − exp(−κ W) that joins every longwave band — including the window,
 * which is transparent only where there is no cloud. In the shortwave
 * the column's cloud depth τ', the sum of each layer's, reflects the beam
 * by the two-stream reflectance τ' / (τ' + 2μ), Coakley and Chýlek's
 * (1975) for conservative scattering with the upscatter fraction
 * (1 − g)/2, so that τ' is (1 − g) times the extinction optical depth τ.
 * The optics follow the phase (cloudOptics): a layer at temperature T
 * holds the liquid share (T − iceTemperature)/(liquidTemperature −
 * iceTemperature), clipped to [0, 1], all liquid above 273.15 K and all
 * ice below 235.15 K, the homogeneous freezing point: half of it is liquid
 * at −19 °C, near the −20 °C at which CALIPSO finds half of the cloud
 * tops supercooled liquid (Hu et al. 2010). Liquid: τ = 3 W/(2 ρ_w r_e) with the
 * effective radius r_e `seaDropletRadius` (11.8 µm) over sea and ice
 * sheets and `landDropletRadius` (8.5 µm) over land, the ISCCP maritime
 * and continental means of Han et al. (1994); g = 0.829 + 2.482·10⁻³ r_e
 * (Slingo 1989, 0.25–0.69 µm); κ = `diffusivity` (1.66) times
 * `liquidInfrared` (0.090361 m²/g, CAM3). Ice: r_e from Ou and Liou's
 * (1995) effective size D_e = 326.3 + 12.42 T_c + 0.197 T_c² +
 * 0.0012 T_c³ µm at the layer's T_c in °C held to [`iceFitColdest`,
 * `iceFitWarmest`] (−60 to −20 °C, their fit's range), r_e = D_e/2
 * (15.55 to 73.55 µm); τ = IWP (3.448·10⁻³ + 2.431/r_e) and
 * g = 0.7661 + 5.851·10⁻⁴ r_e (Ebert and Curry 1992, IWP in g/m², r_e in
 * µm, 0.25–0.69 µm); κ = 1.66 (0.005 + 1/r_e) m²/g (Ebert and Curry). A
 * mixed layer's τ' and κ are the shares' sums. The deck's water and the
 * cumulus are condensate of their layer like any other. Set, the gray
 * optics `cloudScattering` (m²/kg) give every cloud τ' = cloudScattering
 * × W and `cloudAbsorption` (m²/kg) κ = cloudAbsorption in their place.
 * What reaches the surface is direct beam,
 * exp(−τ'/μ) of it less the `skylight` fraction, and diffuse light, the
 * rest; the surface reflects each with its own albedo, and
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
 * cloud water overcast. Under the moist physics' uniform condensation
 * (useCondensation, its `condensation` settings) such a layer instead
 * covers sqrt(q_c/b) up to 1, the cover of the condensate the moist
 * physics' distribution of half-width b = a (1 − RH_c) q_s(T_l) holds at
 * the layer's θ_l, its saturation and critical humidity (see
 * moist.module.js), and the variance cover below takes its saturation over
 * ice too where the moist physics does. Cloud under a strong inversion is stratiform and
 * uniform (the EIS cover of Wood & Bretherton reaches 1 near 11 K), so
 * where the column's estimated inversion strength, the deck's EIS of the
 * lowest layer's air, rises through `overcastInversion` (8 to 12 K) f
 * blends linearly into the f of a distribution whose half-width is also
 * at most the layer's cloud water but not below `overcastWater`
 * (5·10⁻⁵ kg/kg; null: no bound): there a saturated layer holding more
 * cloud water than that is overcast; under the uniform condensation that
 * distribution is of the saturation excess a (q_t − q_s(T_l)) at the
 * cover's saturation, its half-width at most a (1 − RH_c) q_s(T_l), RH_c
 * the moist physics' profile above the mixing top and
 * `boundaryCriticalHumidity` below. Each cell with open sea keeps that
 * share of the ramp in `stratiform`, whatever the cover, for the boundary
 * layer (0 over land and full ice, where no deck forms). Its emissivity is
 * f (1 − exp(−κ × path / f)), and the shortwave is the
 * blend, at the column's cover f̄, of the clear column and the column
 * whose cloud path lies in f̄, as the deck below blends its two columns.
 * Each layer is seen through its visibility 1 − exp(−path /
 * VISIBLE_PATH), 1 g/m². With `cloudOverlap` 'maximumRandom' the
 * layers of each run of adjacent cloudy layers overlap
 * maximally and the runs randomly: f̄ is 1 − Π(1 − f_run), f_run the
 * largest of its layers' f times their visibility; 'maximum' overlaps
 * every layer maximally, f̄ the largest over the column. With
 * 'exponentialRandom' (the default) adjacent layers overlap with
 * α = exp(−Δz/z₀) between maximum and random (Hogan and Illingworth
 * 2000), Δz the distance between their midpoints and the decorrelation
 * length z₀ = `decorrelationLength` − `decorrelationSlope` |latitude in
 * degrees| (2899 m and 27.59 m, Shonk et al. 2010, from CloudSat and
 * CALIPSO), layers apart from each other by clear layers randomly.
 * 'overcast' gives every cloudy layer the whole cell. With `cumulusCloud`
 * and the moist physics' shallow cumulus (`useCumulus`) a plume layer adds
 * its cumulus fraction times the plume's condensate to its water and
 * covers the larger of that fraction and its resolved cloud's f; it is
 * seen at least as its plume, the fraction times the visibility of the
 * plume's own path, so that its in-cloud path is the plume's however
 * small the fraction.
 *
 * In the correlated longwave, `longwaveOverlap` 'exponentialRandom' (the
 * default) overlaps the cloudy layers as the shortwave's cover does: in
 * each g-point every layer is two regions, its cloud's f with emissivity
 * 1 − (1 − ε_gas)(1 − exp(−κ W/f)) and the clear rest with ε_gas, and the
 * flux leaving a region of one layer enters the regions of the next in
 * proportion to the pair's joint areas over the source region's, both
 * cloudy f_a + f_b − C with C = α max(f_a, f_b) + (1 − α)(f_a + f_b −
 * f_a f_b) and α the shortwave's (Hogan and Illingworth 2000; the
 * two-region transfer of Shonk and Hogan 2008); each layer's heating is
 * the divergence of the summed net flux. The deck's layer is one region
 * with its blended emissivity. 'random' takes each layer's mean emissivity
 * f (1 − exp(−κ W/f)), the expectation under random overlap (α = 0).
 *
 * With `boundaryCover` 'variance' (the default; 'pdf' keeps the cover
 * above) a cloudy layer whose midpoint lies below the moist boundary
 * layer's mixing top (`useBoundaryLayer`) covers the Gaussian share of
 * Sommeria and Deardorff (1977) and Mellor (1977) above saturation,
 * f = ½ [1 + erf(Q₁/√2)], Q₁ = s/σ_s, s = a_l (q_t − q_sl(T_l)) the
 * saturation deficit of the layer's θ_l and q_t, a_l = 1/(1 + (L/c_p)
 * dq_s/dT), and σ_s = max(`varianceFloor` q_sl, `varianceScale` l a_l
 * |∂q_t/∂z − Π (dq_s/dT) ∂θ_l/∂z|), the turbulence's spread of s from
 * its mixing length l = κ z/(1 + κ z/`mixingLength`) (z above the
 * surface; 300 m, or `stableMixingLength` 30 m over a surface whose
 * buoyancy flux, the boundary layer's of the step before, is not
 * positive), `varianceScale` 5 (above Mellor and Yamada's B2^½ of
 * 3.2; see M22) and the mean gradients to the neighbouring layers
 * inside the mixed layer; erf is Abramowitz and Stegun's 7.1.26. Such a
 * layer is half covered at s = 0 and overcast once s exceeds 2–3 σ_s,
 * and blends on the overcast ramp of the estimated inversion strength
 * into the bounded cover above as the humidity PDF does. The layers
 * above keep the cover above. The longwave heating of each layer
 * (W/m²) is kept in `longwave` for the
 * boundary layer's cloud-top cooling. Per cell, for the audit (this
 * engine only): `lowCover`, the overlapped cover of the layers below
 * LOW_CLOUD_PRESSURE (680 hPa) combined at random with the deck's, and
 * `lowWater`, their cloud water with the deck's share (kg/m²).
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
 * the subcloud layer, well below the inversion that caps the cloud layer,
 * and a deck started there finds no jump; 'depth' rests it at the
 * boundary-layer top. 'regime' (the default) reads the moist boundary
 * layer's regime of the step before (`useBoundaryLayer`): a coupled
 * column rests at the boundary-layer top, which is then its cloud top; a
 * surface-driven or decoupled column whose ceiling's θ_v jump lies above
 * `cumulusCeiling` (2000 m, the interface's height above the surface), or
 * that has no ceiling, is a cumulus layer whose inversion the plumes have
 * lifted: its gate shuts (G = 0), it runs no deck and rests at the
 * boundary-layer top; any other column rests as under 'inversion', as
 * does every column without the moist boundary layer. Only h is carried: θ_l and q_t are the
 * column's again at each step, and the deck acts on the column through
 * its radiation and through the boundary layer's mixing depth, mlmTop —
 * the deck's h in the boundary layer's height coordinate (that of its
 * `depth`) where the deck ran with a carried height, 0 elsewhere (see
 * boundaryLayer.module.js).
 * The inversion ceiling keeps the deck under the column's own inversion,
 * the lowest interface whose upper layer's midpoint lies above the
 * boundary-layer top (and whose lower one's below maximumHeight) across
 * which θ_v rises by `ceilingInversion` (null, the default:
 * minimumInversion, so that a weaker jump below the regime test's does
 * not hold the deck): the ceiling is 1 m below the
 * midpoint of the layer above that interface, so that layer stays the
 * free troposphere the deck entrains. Above the lowest kilometre the
 * layers are thick enough that the free troposphere's own
 * stratification across one of them passes a 2 K test, and a deck
 * entraining under such a weak jump would deepen into it. A column with
 * no such interface has no ceiling but maximumHeight.
 * With `deckRegime` 'boundaryLayer' the regime test is instead the moist
 * boundary layer's diagnosis of the step before: a coupled
 * stratocumulus-topped layer (REGIME.COUPLED) passes, and with
 * `deckBypass` a passing column runs no mixed-layer model, so that its
 * deck is the resolved cloud at the layer's top under the cover above
 * and the gate only vetoes convection. 'inversion' (the default):
 * The regime test is the capping inversion: Δθ_v ≥ minimumInversion
 * (4 K) at the start's h. Resting at the ceiling under a 2 K test, the
 * deck's gate stood open over 0.28 of the globe and 0.21 of 10S–10N
 * three days on from eight64_day0183, trade cumulus included; under 4 K
 * over 0.15 and 0.07. A stratocumulus-topped layer also
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
 * Shortwave gases: with solarGases 'clirad' (the default) and a humidity
 * field, ozone, water vapour, O₂ and CO₂ absorb the beam by
 * js/physics/shortwaveGases.module.js (after CLIRAD-SW) along their paths
 * from the top lengthened by the magnification 35/√(1224μ² + 1) of
 * Lacis & Hansen (1974), each layer taking what its own gas adds; ozone
 * and the visible band's vapour come out of the visible part of the
 * beam, the rest out of the near infrared. With `ozone` 'afgl' (the
 * default) each layer's ozone is that of the climatology of
 * js/physics/ozone.module.js (the AFGL atmospheres by latitude and season)
 * between its interfaces' pressures; with 'idealized' the column is
 * ozoneColumn[0] + (ozoneColumn[1] − ozoneColumn[0]) sin²lat (cm-atm)
 * distributed in the vertical as the column above σ below; an
 * `ozoneProfile` of per-layer amounts replaces either, for one-column
 * use. With upwardAbsorption the visible light the surface sends out of
 * the column, and the visible light cloud and clear air reflect, also
 * cross the ozone column at 5/3, losing
 * 1 − exp(−0.0542 Ω 5/3) (the ozone coefficients of the 0.32-0.7 µm bands
 * weighted by their share), and the near-infrared light the gases'
 * absorption of the path down plus 5/3 of the column's beyond the first.
 * Without humidity, or with 'lacisHansen', the fraction `ozoneAbsorption`
 * of the incoming beam is absorbed aloft. The ozone column follows Lacis & Hansen (1974)
 * (centred at ozoneHeight with width ozoneWidth, heights from σ with the
 * scale height) and the absorbing part of the beam decays through it
 * with the optical depth ozoneOpacity, so the heating peaks above the
 * ozone maximum as it does at the stratopause. Water vapour absorbs the
 * beam below it by the Lacis & Hansen (1974) absorptivity of the water
 * path the beam has crossed — pressure-scaled by √σ and lengthened by
 * their magnification 35/√(1224μ² + 1) — times `vaporAbsorption`, each
 * layer taking what its own vapour adds to the path above it; what is
 * left goes on to the clouds and the surface. Dry air absorbs nothing.
 *
 * Clear air scatters: the part `visibleFraction` (VISIBLE_FRACTION, the
 * share 0.4707 of CLIRAD-SW's bands 1-8 below 0.7 µm, Chou & Suarez 1999,
 * Table 3) of the beam, less what ozone takes, is split into the
 * sub-bands `rayleighBands`, [weight, depth] pairs, each with the
 * scattering depth depth × p_s / REFERENCE_PRESSURE + (1 − aerosolAsymmetry)
 * aerosolAlbedo τ_a added to the cloud's in its own two-stream,
 * conservative Rayleigh scattering with asymmetry 0 and the aerosol's
 * forward peak counted as transmitted; the rest of the beam meets the
 * cloud with the scattering depth nearInfraredRayleigh × p_s /
 * REFERENCE_PRESSURE added. The two default sub-bands, weights 0.7049 and
 * 0.2951 at depths 0.0957 and 0.5806 for a full atmosphere, follow within
 * 0.1 % over μ 0.1–1 the band-mean reflectance over a black surface that
 * scripts/rayleighReference.mjs finds for the molecular atmosphere with
 * the same two-stream at 40 wavelengths (the depths of Hansen & Travis
 * 1974, a 5778 K spectrum from 0.297 to 0.683 µm, the share 0.4407 of the
 * beam this band carries); one grey depth misses it by 10 % at both ends.
 * The near infrared's NEAR_INFRARED_RAYLEIGH, 0.0114, gives the same
 * reference's reflection of the rest of the spectrum to 0.1 % (2.0 W/m²
 * of the 25.2 of the whole spectrum over a black surface). A number for
 * rayleighDepth replaces the sub-bands by one band of that depth.
 * τ_a, the aerosol's mid-visible depth, is landAerosol (0.12) over land
 * and seaAerosol (0.07) over sea and ice sheets; before the two-stream
 * it absorbs 1 − exp(−(1 − aerosolAlbedo) τ_a m) of that part of the beam,
 * m the vapour's magnification, and heats each layer by its share of an
 * aerosol whose density falls off over aerosolHeight (2 km), the fraction
 * σ^(scaleHeight/aerosolHeight) of it lying above σ. With
 * upwardAbsorption (the default) the light the surface reflects out of
 * the column is absorbed on its way up: in the rest of the beam the
 * vapour takes what the path the beam crossed coming down plus 5/3 of the
 * column's pressure-scaled path absorbs beyond the first, as a share of
 * that part of the beam left after the descent, each layer what its own
 * vapour adds going up; in the visible part the aerosol takes
 * 1 − exp(−(1 − aerosolAlbedo) τ_a 5/3), heating the layers as above.
 * The light cloud and clear air reflect leaves unabsorbed by vapour and
 * aerosol. The clear-sky
 * pass has the same scattering and absorption. `skylight` (0) sends that
 * fraction of the beam reaching the surface as diffuse light in every band.
 */
export function erf(x) {
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const y = 1 - t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * Math.exp(-x * x);
  return x < 0 ? -y : y;
}

export const DECORRELATION_LENGTH = 2899, DECORRELATION_SLOPE = 27.59;

/*
 * The cover of the layers down to this one under exponential-random
 * overlap, from that of the layers above (`cumulative`), the layer above's
 * cover `previous`, this layer's `cover` and the overlap parameter alpha
 * between the two (Hogan and Illingworth 2000; the recursion of ecRad).
 */
export function overlapped(cumulative, previous, cover, alpha) {
  if (!(previous < 1)) return 1;
  const pair = alpha * Math.max(previous, cover) + (1 - alpha) * (previous + cover - previous * cover);
  return 1 - (1 - cumulative) * (1 - pair) / (1 - previous);
}

export function varianceCover(deficit, spread) {
  if (!(spread > 0)) return deficit > 0 ? 1 : deficit < 0 ? 0 : 0.5;
  return 0.5 * (1 + erf(deficit / (Math.SQRT2 * spread)));
}

export function waterVaporAbsorptivity(path) {
  return 2.9 * path / (Math.pow(1 + 141.5 * path, 0.635) + 5.925 * path);
}

const DIFFUSE_MU = 0.6;
export const LOW_CLOUD_PRESSURE = 680e2;
export const STABILITY_SIGMA = 0.7;
export const DECK_CLOUD_LEVELS = 8;
export const UNDECIDED = 0.5;
export const VISIBLE_PATH = 1e-3;
export const REFERENCE_PRESSURE = 101325;
export const RAYLEIGH_BANDS = [[0.7049, 0.0957], [0.2951, 0.5806]], NEAR_INFRARED_RAYLEIGH = 0.0114, VISIBLE_FRACTION = 0.4707, LAND_AEROSOL = 0.12, SEA_AEROSOL = 0.07;
const DIFFUSE_PATH = 5 / 3;
export const GREENHOUSE_GASES = { carbonDioxide: 388.75e-6, methane: 1798.93e-9, nitrousOxide: 323.18e-9 };
export const OZONE_COLUMN = [0.26, 0.35];
const VISIBLE_OZONE = (OZONE_SHARES[6] * OZONE_COEFFICIENTS[6] + OZONE_SHARES[7] * OZONE_COEFFICIENTS[7]) / (OZONE_SHARES[6] + OZONE_SHARES[7]);
const SHORTWAVE_KEYS = ['absorbed', 'down', 'direct', 'reflectance', 'cloud', 'visibleEscape', 'restEscape', 'visibleReflectance'];
export const SUMMED = ['absorbedSolar', 'atmosphereSolar', 'outgoingLongwave', 'insolation', 'reflectedSolar'];
export const CLEAR_SUMMED = ['clearAbsorbedSolar', 'clearOutgoingLongwave'];

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

export const CLOUD_OPTICS = {
  liquidTemperature: 273.15, iceTemperature: 235.15, seaDropletRadius: 11.8, landDropletRadius: 8.5, liquidInfrared: 0.090361, diffusivity: 1.66,
  iceFitWarmest: -20, iceFitColdest: -60,
};

export function liquidShare(T, { liquidTemperature, iceTemperature } = CLOUD_OPTICS) {
  return Math.min(1, Math.max(0, (T - iceTemperature) / (liquidTemperature - iceTemperature)));
}

export function iceRadius(T, { iceFitWarmest, iceFitColdest } = CLOUD_OPTICS) {
  const c = Math.min(iceFitWarmest, Math.max(iceFitColdest, T - 273.15));
  return 0.5 * (326.3 + c * (12.42 + c * (0.197 + c * 0.0012)));
}

/*
 * The optics of a kilogram of cloud condensate per m² at temperature T
 * (cloudOptics(T, continental, options) with the options of CLOUD_OPTICS):
 * `liquid` the liquid share, `visible` its mid-visible extinction optical
 * depth, `solar` the two-stream's depth (1 − g) of it, `infrared` its
 * longwave absorption with the diffusivity factor (all m²/kg).
 */
export function cloudOptics(T, continental, options = CLOUD_OPTICS, out = { liquid: 0, visible: 0, solar: 0, infrared: 0 }) {
  const liquid = liquidShare(T, options), droplet = continental ? options.landDropletRadius : options.seaDropletRadius;
  const crystal = iceRadius(T, options);
  const liquidVisible = 1500 / droplet, iceVisible = 1000 * (3.448e-3 + 2.431 / crystal);
  out.liquid = liquid;
  out.visible = liquid * liquidVisible + (1 - liquid) * iceVisible;
  out.solar = liquid * (1 - (0.829 + 2.482e-3 * droplet)) * liquidVisible + (1 - liquid) * (1 - (0.7661 + 5.851e-4 * crystal)) * iceVisible;
  out.infrared = options.diffusivity * 1000 * (liquid * options.liquidInfrared + (1 - liquid) * (0.005 + 1 / crystal));
  return out;
}

export function createRadiation(mesh, core, {
  solarConstant = SOLAR_CONSTANT, albedo = 0.07, cloudAbsorption = null, cloudScattering = null, stratus = true, stratusIndex = 'eis', stratusScale = 0.15, stratusWaterMax = 0.15, stratusSigma = 0.92,
  mixedLayerDeck = true, mixedLayer: mixedLayerOptions = {}, stratusSubsidence = -1e-3, minimumInversion = 4, ceilingInversion = null, subsidenceMemory = 2 * DAY, stratusSolar = true, cloudSolarAbsorption = 0.4,
  prognosticHeight = true, deckRest = 'regime', cumulusCeiling = 2000, gateMemory = DAY, subsidenceSmoothing = 2, cloudCover = 'pdf', criticalHumidity = 0.8, boundaryCriticalHumidity = 0.85, coverFloor = 0.01, overcastWater = 5e-5, overcastInversion = [8, 12], cloudOverlap = 'exponentialRandom', decorrelationLength = DECORRELATION_LENGTH, decorrelationSlope = DECORRELATION_SLOPE,
  boundaryCover = 'variance', varianceFloor = 0.002, varianceScale = 5, mixingLength = 300, stableMixingLength = 30, deckRegime = 'inversion', deckBypass = false,
  cumulusCloud = true, window = 0.25, tauEquator = 5.3, tauPole = 1.325, linearFraction = 0.1, gasFraction = 0.2, gasOpticalDepth = 7,
  ozoneAbsorption = 0.03, ozoneHeight = 25e3, ozoneWidth = 5e3, ozoneOpacity = 4, scaleHeight = 7e3, vaporAbsorption = 1,
  exchangeCoefficient = SEA_DRAG, exchangeCoefficients = null, referenceCoefficients = null, surfaceLayer = false, gustiness = 3, latentHeat = LATENT_HEAT, vaporCoupling = 0.55, skylight = 0, clearSkyPass = false, buffers = null,
  longwaveScheme = 'correlated', longwaveOverlap = 'exponentialRandom', solarGases = 'clirad', carbonDioxide = GREENHOUSE_GASES.carbonDioxide, methane = GREENHOUSE_GASES.methane, nitrousOxide = GREENHOUSE_GASES.nitrousOxide, ozone = 'afgl', ozoneColumn = OZONE_COLUMN, ozoneProfile = null, vaporStrength = VAPOR_STRENGTH,
  rayleighBands = RAYLEIGH_BANDS, rayleighDepth = null, nearInfraredRayleigh = NEAR_INFRARED_RAYLEIGH, upwardAbsorption = true, visibleFraction = VISIBLE_FRACTION, landAerosol = LAND_AEROSOL, seaAerosol = SEA_AEROSOL, aerosolAlbedo = 0.95, aerosolAsymmetry = 0.7, aerosolHeight = 2000, land = null, iceSheet = null,
  liquidTemperature = CLOUD_OPTICS.liquidTemperature, iceTemperature = CLOUD_OPTICS.iceTemperature, seaDropletRadius = CLOUD_OPTICS.seaDropletRadius, landDropletRadius = CLOUD_OPTICS.landDropletRadius,
  liquidInfrared = CLOUD_OPTICS.liquidInfrared, diffusivity = CLOUD_OPTICS.diffusivity, iceFitWarmest = CLOUD_OPTICS.iceFitWarmest, iceFitColdest = CLOUD_OPTICS.iceFitColdest,
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
  if (longwaveScheme !== 'correlated' && longwaveScheme !== 'gray') throw new Error(`longwaveScheme must be 'correlated' or 'gray', not ${longwaveScheme}`);
  if (longwaveOverlap !== 'exponentialRandom' && longwaveOverlap !== 'random') throw new Error(`longwaveOverlap must be 'exponentialRandom' or 'random', not ${longwaveOverlap}`);
  if (solarGases !== 'clirad' && solarGases !== 'lacisHansen') throw new Error(`solarGases must be 'clirad' or 'lacisHansen', not ${solarGases}`);
  if (ozone !== 'afgl' && ozone !== 'idealized') throw new Error(`ozone must be 'afgl' or 'idealized', not ${ozone}`);
  if (ozoneProfile && ozoneProfile.length !== K) throw new Error(`ozoneProfile must give one ozone amount per layer, ${K}, not ${ozoneProfile.length}`);
  const ozoneShare = Float64Array.from({ length: K }, (_, k) => ozoneAbove(levels[k + 1]) - ozoneAbove(levels[k]));
  const ozoneCell = Float64Array.from({ length: C }, (_, i) => ozoneColumn[0] + (ozoneColumn[1] - ozoneColumn[0]) * Math.sin(mesh.latCell[i]) ** 2);
  const wellMixed = [carbonDioxide * GAS_MOLAR.co2 / GAS_MOLAR.air, methane * GAS_MOLAR.ch4 / GAS_MOLAR.air, nitrousOxide * GAS_MOLAR.n2o / GAS_MOLAR.air];
  const layerOzone = new Float64Array(K), paths = Array.from({ length: 6 }, () => new Float64Array(K)), pathRow = new Float64Array(6);
  const gasEmissivityG = new Float64Array(K), totalEmissivity = new Float64Array(K), planck = new Float64Array(K);
  const chained = longwaveOverlap === 'exponentialRandom';
  const chainCover = new Float64Array(K), chainCloud = new Float64Array(K), chainOutside = new Float64Array(K), clearShareInverse = new Float64Array(K), cloudShareInverse = new Float64Array(K);
  const joint = Array.from({ length: 4 }, () => new Float64Array(K)), leavingUp = new Float64Array(K), leavingDown = new Float64Array(K);
  const ozoneTaken = new Float64Array(K), visibleVaporTaken = new Float64Array(K), nearInfraredTaken = new Float64Array(K), upwardOzone = new Float64Array(K);
  const gasSplit = { vapour: 0, oxygen: 0, co2: 0 };
  const rayleigh = rayleighDepth !== null ? [[1, rayleighDepth]] : rayleighBands;
  if (!(rayleigh.length >= 1 && rayleigh.length <= 3 && Math.abs(rayleigh.reduce((s, [w]) => s + w, 0) - 1) < 1e-9 && rayleigh.every(([w, tau]) => w > 0 && tau >= 0))) throw new Error(`rayleighBands must be one to three [weight, depth] pairs whose weights sum to 1, not ${JSON.stringify(rayleigh)}`);
  const scatters = rayleigh.some(([, tau]) => tau > 0) || nearInfraredRayleigh > 0 || landAerosol > 0 || seaAerosol > 0;
  const aerosolCell = Float64Array.from({ length: C }, (_, i) => (land && land[i] && !(iceSheet && iceSheet[i]) ? landAerosol : seaAerosol));
  const aerosolFraction = Float64Array.from({ length: K }, (_, k) => levels[k + 1] ** (scaleHeight / aerosolHeight) - levels[k] ** (scaleHeight / aerosolHeight));
  const light = { share: 0, pressure: 0, aerosol: 0 }, visibleLight = { absorbed: 0, down: 0, direct: 0, reflectance: 0, cloud: 0 };
  const upwardVapor = new Float64Array(K);
  const emissivity = new Float64Array(K);
  const cloudEmissivity = new Float64Array(K);
  const layerCover = new Float64Array(K).fill(1);
  let cumulusCover = null, cumulusWater = null, condensing = null;
  const saturated = { qs: 0, slope: 0, liquid: 1 };
  const clearSky = { absorbed: 0, down: 0, direct: 0, reflectance: 0, cloud: 0, visibleEscape: 0, restEscape: 0, visibleReflectance: 0 };
  const vaporEmissivity = new Float64Array(K);
  const mixedEmissivity = new Float64Array(K);
  const surfaceFlux = new Float64Array(C), surfaceDirect = new Float64Array(C);
  const outgoingBuffer = buffers && buffers.outgoing ? buffers.outgoing : new SharedArrayBuffer(8 * C);
  const shortwaveBuffer = buffers && buffers.surfaceShortwave ? buffers.surfaceShortwave : new SharedArrayBuffer(8 * C);
  const outgoing = new Float64Array(outgoingBuffer), surfaceShortwave = new Float64Array(shortwaveBuffer);
  const evaporationBuffer = buffers && buffers.evaporation ? buffers.evaporation : new SharedArrayBuffer(8 * C);
  const evaporation = new Float64Array(evaporationBuffer), potentialEvaporation = new Float64Array(C);
  const surfaceBuoyancyBuffer = buffers && buffers.surfaceBuoyancy ? buffers.surfaceBuoyancy : new SharedArrayBuffer(8 * C);
  const surfaceBuoyancy = new Float64Array(surfaceBuoyancyBuffer);
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
  const summedBuffers = Object.fromEntries([...SUMMED, ...CLEAR_SUMMED].map((name) => [name, buffers && buffers.summed && buffers.summed[name] ? buffers.summed[name] : new SharedArrayBuffer(8 * C)]));
  const summed = Object.fromEntries([...SUMMED, ...CLEAR_SUMMED].map((name) => [name, new Float64Array(summedBuffers[name])]));
  const meanAbsorbedSolar = new Float64Array(C), meanOutgoingLongwave = new Float64Array(C), meanPlanetaryAlbedo = new Float64Array(C);
  const meanShortwaveCloudEffect = new Float64Array(C), meanLongwaveCloudEffect = new Float64Array(C);
  if (!(buffers && buffers.mlmGate)) mlmGate.fill(UNDECIDED);
  const shadow = mixedLayerDeck ? createMixedLayer({ cp, R, g, latentHeat, referencePressure: p0, cloudLevels: DECK_CLOUD_LEVELS, ...mixedLayerOptions }) : null;
  const shadowLongwave = dycomsLongwave();
  if (stratusIndex !== 'eis' && stratusIndex !== 'ectei') throw new Error(`stratusIndex must be 'eis' or 'ectei', not ${stratusIndex}`);
  if (cloudOverlap !== 'maximum' && cloudOverlap !== 'maximumRandom' && cloudOverlap !== 'exponentialRandom') throw new Error(`cloudOverlap must be 'maximum', 'maximumRandom' or 'exponentialRandom', not ${cloudOverlap}`);
  const exponential = cloudOverlap === 'exponentialRandom';
  const decorrelation = Float64Array.from({ length: C }, (_, i) => decorrelationLength - decorrelationSlope * Math.abs(mesh.latCell[i]) * 180 / Math.PI);
  if (!(overcastInversion?.[1] > overcastInversion?.[0])) throw new Error(`overcastInversion must rise from its first to its second EIS, not ${overcastInversion}`);
  if (![0, 1, 2].includes(subsidenceSmoothing)) throw new Error(`subsidenceSmoothing must be 0, 1 or 2, not ${subsidenceSmoothing}`);
  if (deckRest !== 'depth' && deckRest !== 'inversion' && deckRest !== 'regime') throw new Error(`deckRest must be 'depth', 'inversion' or 'regime', not ${deckRest}`);
  if (boundaryCover !== 'variance' && boundaryCover !== 'pdf') throw new Error(`boundaryCover must be 'variance' or 'pdf', not ${boundaryCover}`);
  if (deckRegime !== 'inversion' && deckRegime !== 'boundaryLayer') throw new Error(`deckRegime must be 'inversion' or 'boundaryLayer', not ${deckRegime}`);
  const longwaveBuffer = buffers && buffers.longwave ? buffers.longwave : new SharedArrayBuffer(8 * K * C);
  const longwave = new Float64Array(longwaveBuffer), beforeBands = new Float64Array(K);
  let boundaryRegime = null, boundaryTop = null, boundaryBuoyancy = null;
  const lowCover = new Float64Array(C), lowWater = new Float64Array(C);
  const ceilingJump = ceilingInversion ?? minimumInversion;
  const entraining = stratusIndex === 'ectei';
  const stratusLayer = nearestLayer(sigmaMid, stratusSigma), stabilityLayer = nearestLayer(sigmaMid, STABILITY_SIGMA);
  const gasEmissivity = Float64Array.from({ length: K }, (_, k) => 1 - Math.exp(-gasOpticalDepth * (levels[k + 1] - levels[k])));
  const temperature = new Float64Array(K);
  const optics = { liquidTemperature, iceTemperature, seaDropletRadius, landDropletRadius, liquidInfrared, diffusivity, iceFitWarmest, iceFitColdest };
  const graySolar = cloudScattering !== null, grayInfrared = cloudAbsorption !== null;
  const continental = Uint8Array.from({ length: C }, (_, i) => (land && land[i] && !(iceSheet && iceSheet[i]) ? 1 : 0));
  const solarDepth = new Float64Array(K), infrared = new Float64Array(K), layerOptics = { liquid: 0, visible: 0, solar: 0, infrared: 0 };
  let cloudMask = null;
  const cloudWater = new Float64Array(K);
  const vaporTaken = new Float64Array(K);
  const emitted = new Float64Array(K);
  const netFlux = new Float64Array(K);
  const sun = new Float64Array([1, 0, 0]), ozoneWeight = new Float64Array(5);
  let yearFraction = 0;
  const budget = { absorbedSolar: 0, atmosphereSolar: 0, aerosolSolar: 0, outgoingLongwave: 0, clearAbsorbedSolar: 0, clearOutgoingLongwave: 0, sensibleHeat: 0, evaporation: 0, surfaceBuoyancy: 0, surfaceFlux: 0, insolation: 0, reflectedSolar: 0, cloudReflectance: 0, cloudCover: 0, cloudSolar: 0, stratus: 0, stratusFraction: 0, stabilityIndex: NaN, mlmCover: 0, mlmWater: 0, mlmEntrainment: 0, mlmSolar: 0, mlmTop: 0, stratiform: 0, ozoneSolar: 0, vaporSolar: 0, oxygenSolar: 0, carbonDioxideSolar: 0, upwardGasSolar: 0, downwardLongwave: 0 };
  const sky = { absorbed: 0, down: 0, direct: 0, reflectance: 0, cloud: 0, visibleEscape: 0, restEscape: 0, visibleReflectance: 0 }, decked = { ...sky }, probe = { ...sky };
  const deckLight = { incident: 0, mu: 0, direct: 0, diffuse: 0, path: 0, depth: 0, layer: 0, clear: 0 };

  function setTime(t) {
    sunDirection(t, sun);
    yearFraction = (t % YEAR) / YEAR;
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

  function upward(fraction, eps, surfaceEmission) {
    let up = fraction * surfaceEmission;
    for (let k = K - 1; k >= 0; k--) up = up * (1 - eps[k]) + fraction * eps[k] * STEFAN_BOLTZMANN * temperature[k] ** 4;
    return up;
  }

  function bandPlanck(source, eps, surfaceUp) {
    let down = 0;
    for (let k = 0; k < K; k++) {
      netFlux[k] += eps[k] * (down - 2 * source[k]);
      down = down * (1 - eps[k]) + eps[k] * source[k];
    }
    let up = surfaceUp;
    for (let k = K - 1; k >= 0; k--) {
      netFlux[k] += eps[k] * up;
      up = up * (1 - eps[k]) + eps[k] * source[k];
    }
    return [up, down];
  }

  function upwardPlanck(source, eps, surfaceUp) {
    let up = surfaceUp;
    for (let k = K - 1; k >= 0; k--) up = up * (1 - eps[k]) + eps[k] * source[k];
    return up;
  }

  function chainPlanck(source, gas, surfaceUp) {
    const [p00, p01, p10, p11] = joint;
    let up0 = 0, up1 = 0;
    for (let k = K - 1; k >= 0; k--) {
      const f = chainCover[k], clear = 1 - (1 - gas[k]) * (1 - chainOutside[k]), cloudy = 1 - (1 - gas[k]) * (1 - chainCloud[k]);
      let in0 = (1 - f) * surfaceUp, in1 = f * surfaceUp;
      if (k < K - 1) {
        const from0 = up0 * clearShareInverse[k + 1], from1 = up1 * cloudShareInverse[k + 1];
        in0 = from0 * p00[k + 1] + from1 * p01[k + 1]; in1 = from0 * p10[k + 1] + from1 * p11[k + 1];
      }
      up0 = in0 * (1 - clear) + clear * source[k] * (1 - f);
      up1 = in1 * (1 - cloudy) + cloudy * source[k] * f;
      leavingUp[k] = up0 + up1;
    }
    let down0 = 0, down1 = 0;
    for (let k = 0; k < K; k++) {
      const f = chainCover[k], clear = 1 - (1 - gas[k]) * (1 - chainOutside[k]), cloudy = 1 - (1 - gas[k]) * (1 - chainCloud[k]);
      let in0 = 0, in1 = 0;
      if (k > 0) {
        const from0 = down0 * clearShareInverse[k - 1], from1 = down1 * cloudShareInverse[k - 1];
        in0 = from0 * p00[k] + from1 * p10[k]; in1 = from0 * p01[k] + from1 * p11[k];
      }
      down0 = in0 * (1 - clear) + clear * source[k] * (1 - f);
      down1 = in1 * (1 - cloudy) + cloudy * source[k] * f;
      leavingDown[k] = down0 + down1;
    }
    for (let k = 0; k < K; k++) netFlux[k] += (k > 0 ? leavingDown[k - 1] : 0) - leavingUp[k] - leavingDown[k] + (k < K - 1 ? leavingUp[k + 1] : surfaceUp);
    return [leavingUp[0], leavingDown[K - 1]];
  }

  function solarGasPaths(i, pi, theta, q, mu) {
    const magnification = 35 / Math.sqrt(1224 * mu * mu + 1);
    let ozone = 0, vapour = 0, oxygen = 0, co2 = 0;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i, p = pi * sigmaMid[k], mass = pi * dSigma[k] / g, dry = mass * Math.max(0, 1 - Math.max(0, q[idx])), scaling = pressureScaling(p) * magnification;
      ozone += layerOzone[k] * magnification;
      vapour += Math.max(0, q[idx]) * mass * 0.1 * vaporScaling(p, theta[idx] * exnerLayer[idx]) * magnification;
      oxygen += OXYGEN.mixingRatio * STP_DEPTH * dry * scaling;
      co2 += carbonDioxide * STP_DEPTH * dry * scaling;
      ozoneTaken[k] = ozoneAbsorptivity(ozone);
      visibleVaporTaken[k] = vaporAbsorption * visibleVaporAbsorptivity(vapour, vaporStrength);
      nearInfraredTaken[k] = vaporAbsorption * nearInfraredVaporAbsorptivity(vapour, vaporStrength) + oxygenAbsorptivity(oxygen) + carbonDioxideAbsorptivity(co2);
    }
    gasSplit.vapour = vaporAbsorption * (visibleVaporAbsorptivity(vapour, vaporStrength) + nearInfraredVaporAbsorptivity(vapour, vaporStrength));
    gasSplit.oxygen = oxygenAbsorptivity(oxygen);
    gasSplit.co2 = carbonDioxideAbsorptivity(co2);
    return { vapour, oxygen, co2 };
  }

  function nearInfraredUpward(i, pi, theta, q, down) {
    let vapour = down.vapour, oxygen = down.oxygen, co2 = down.co2;
    const total = () => vaporAbsorption * nearInfraredVaporAbsorptivity(vapour, vaporStrength) + oxygenAbsorptivity(oxygen) + carbonDioxideAbsorptivity(co2);
    let before = total(), loss = 0;
    for (let k = K - 1; k >= 0; k--) {
      const idx = k * C + i, p = pi * sigmaMid[k], mass = pi * dSigma[k] / g, dry = mass * Math.max(0, 1 - Math.max(0, q[idx])), scaling = pressureScaling(p) * DIFFUSE_PATH;
      vapour += Math.max(0, q[idx]) * mass * 0.1 * vaporScaling(p, theta[idx] * exnerLayer[idx]) * DIFFUSE_PATH;
      oxygen += OXYGEN.mixingRatio * STP_DEPTH * dry * scaling;
      co2 += carbonDioxide * STP_DEPTH * dry * scaling;
      const through = total();
      upwardVapor[k] = through - before;
      loss += upwardVapor[k];
      before = through;
    }
    return loss;
  }

  function deckWater(lcl, base, mixedDepth) {
    const thickness = mixedDepth - base;
    if (thickness <= 0) return 0;
    return Math.min(stratusWaterMax, stratusScale * 0.5 * adiabaticWaterLapse(lcl.temperature, lcl.pressure, cp, R, g, latentHeat) * thickness * thickness);
  }

  function inversionLayer(i, surface, floor) {
    for (let k = K - 2; k >= 1; k--) {
      if ((geopotential[(k + 1) * C + i] - surface) / g >= shadow.maximumHeight) break;
      if ((geopotential[k * C + i] - surface) / g > floor && thetaV[k * C + i] - thetaV[(k + 1) * C + i] >= ceilingJump) return k;
    }
    return -1;
  }

  function shadowDeck(i, pi, theta, q, qc, mixedDepth, sensible, evaporation, dt, absorbedSolar) {
    const bottom = (K - 1) * C + i;
    const surface = geopotential[bottom] - cp * thetaV[bottom] * (exnerLower[bottom] - exnerLayer[bottom]);
    const depth = mixedDepth + (geopotential[bottom] - surface) / g;
    const interfaceHeight = (m) => (geopotential[m * C + i] + cp * thetaV[m * C + i] * (exnerLayer[m * C + i] - exnerLower[(m - 1) * C + i]) - surface) / g;
    const regime = deckRest === 'regime' && boundaryRegime !== null ? boundaryRegime[i] : -1;
    const lifted = regime === REGIME.SURFACE || regime === REGIME.DECOUPLED;
    const capping = prognosticHeight || lifted ? inversionLayer(i, surface, depth) : -1;
    const ceiling = prognosticHeight && capping >= 0 ? (geopotential[capping * C + i] - surface) / g - 1 : Infinity;
    const standDown = lifted && !(capping >= 0 && interfaceHeight(capping + 1) <= cumulusCeiling);
    const resting = deckRest !== 'depth' && regime !== REGIME.COUPLED && !standDown && ceiling < shadow.maximumHeight ? ceiling : depth;
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
    let now = sinking && !standDown ? shadow.diagnose(start, forcing) : null;
    const pass = sinking && !standDown && (deckRegime === 'boundaryLayer' ? boundaryRegime !== null && boundaryRegime[i] === REGIME.COUPLED : now.virtualJump >= minimumInversion) ? 1 : 0;
    const gate = standDown ? 0 : gateMemory > 0 ? mlmGate[i] - (pass - mlmGate[i]) * Math.expm1(-dt / gateMemory) : pass;
    mlmGate[i] = gate;
    if (!(gate > UNDECIDED || (gate === UNDECIDED && pass === 1)) || deckBypass) return rest();
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
    stream(out, cloudDepth + nearInfraredRayleigh * light.pressure, keep, mu, surfaceAlbedo, diffuseAlbedo);
    const escape = 1 - out.absorbed - out.reflectance - out.cloud;
    if (!scatters) { out.visibleEscape = light.share * escape; out.restEscape = (1 - light.share) * escape; out.visibleReflectance = light.share * out.reflectance; return; }
    const rest = 1 - light.share;
    for (const key of ['absorbed', 'down', 'direct', 'reflectance', 'cloud']) out[key] *= rest;
    out.restEscape = rest * escape; out.visibleEscape = 0; out.visibleReflectance = 0;
    for (const [weight, depth] of rayleigh) {
      stream(visibleLight, cloudDepth + depth * light.pressure + light.aerosol, keep, mu, surfaceAlbedo, diffuseAlbedo);
      const share = light.share * weight;
      for (const key of ['absorbed', 'down', 'direct', 'reflectance', 'cloud']) out[key] += share * visibleLight[key];
      out.visibleEscape += share * (1 - visibleLight.absorbed - visibleLight.reflectance - visibleLight.cloud);
      out.visibleReflectance += share * visibleLight.reflectance;
    }
  }

  function stream(out, cloudDepth, keep, mu, surfaceAlbedo, diffuseAlbedo) {
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
    shortwave(probe, graySolar ? cloudScattering * total : deckLight.depth + solarDepth[stratusLayer] * water, Math.exp(-cloudSolarAbsorption * total), deckLight.mu, deckLight.direct, deckLight.diffuse);
    return deckLight.incident * probe.cloud / total * (deckLight.layer + water) - deckLight.clear;
  }

  function turbulentCover(i, k, pi, theta, q, qc, mixingDepth) {
    const bottom = (K - 1) * C + i, idx = k * C + i, ex = exnerLayer[idx], water = Math.max(0, qc[idx]), total = Math.max(0, q[idx]) + water;
    const level = theta[idx] - latentHeat * water / (cp * ex), liquidT = level * ex;
    const iced = condensing !== null && condensing.iceSaturation;
    const qs = iced ? cloudSaturation(liquidT, pi * sigmaMid[k], true, condensing.liquidTemperature, condensing.iceTemperature, saturated).qs : saturationHumidity(liquidT, pi * sigmaMid[k]);
    const slope = iced ? saturated.slope : qs * latentHeat / (R_VAPOR * liquidT * liquidT);
    const a = 1 / (1 + latentHeat / cp * slope);
    let gradientQ = 0, gradientL = 0, n = 0;
    for (const j of [k - 1, k + 1]) {
      if (j < 0 || j > K - 1 || (j < k && (geopotential[j * C + i] - geopotential[bottom]) / g >= mixingDepth)) continue;
      const jdx = j * C + i, jWater = Math.max(0, qc[jdx]), dz = (geopotential[jdx] - geopotential[idx]) / g;
      gradientQ += (Math.max(0, q[jdx]) + jWater - total) / dz;
      gradientL += (theta[jdx] - latentHeat * jWater / (cp * exnerLayer[jdx]) - level) / dz;
      n++;
    }
    if (n > 1) { gradientQ /= n; gradientL /= n; }
    const surface = geopotential[bottom] - cp * thetaV[bottom] * (exnerLower[bottom] - exnerLayer[bottom]), z = (geopotential[idx] - surface) / g;
    const asymptote = boundaryBuoyancy !== null && !(boundaryBuoyancy[i] > 0) ? stableMixingLength : mixingLength;
    const length = 0.4 * z / (1 + 0.4 * z / asymptote);
    const spread = Math.max(varianceFloor * qs, varianceScale * length * a * Math.abs(gradientQ - ex * slope * gradientL));
    return varianceCover(a * (total - qs), spread);
  }

  function overcastCover(idx, k, pi, theta, q, qc, rhc) {
    const water = Math.max(0, qc[idx]), pressure = pi * sigmaMid[k];
    if (condensing !== null && condensing.uniform) {
      cloudSaturation(theta[idx] * exnerLayer[idx] - latentHeat * water / cp, pressure, condensing.iceSaturation, condensing.liquidTemperature, condensing.iceTemperature, saturated);
      const a = 1 / (1 + latentHeat * saturated.slope / cp), bound = Math.min(a * (1 - rhc) * saturated.qs, Math.max(water, overcastWater));
      return Math.min(1, Math.max(coverFloor, (a * (Math.max(0, q[idx]) + water - saturated.qs) + bound) / (2 * bound)));
    }
    const qs = saturationHumidity(theta[idx] * exnerLayer[idx], pressure), excess = Math.max(0, q[idx]) + water - qs;
    const bound = Math.min((1 - rhc) * qs, Math.max(water, overcastWater));
    return Math.min(1, Math.max(coverFloor, (excess + bound) / (2 * bound)));
  }

  function column(i, pi, theta, surfaceT, windSpeed, tau0 = tauCell[i], beam = insolation(i), qAir = null, q = null, qc = null, surfaceAlbedo = albedo, diffuseAlbedo = surfaceAlbedo, wetness = 1, exchangeCoefficientAt = exchangeCoefficient, openSea = 0, mixedDepth = 0, dt = 0, mixingDepth = 0) {
    const clirad = solarGases === 'clirad' && q !== null, correlated = longwaveScheme === 'correlated' && q !== null;
    let ozoneColumnAmount = 0;
    if ((clirad || correlated) && !ozoneProfile && ozone === 'afgl') {
      ozoneWeights(mesh.latCell[i], yearFraction, ozoneWeight);
      let above = 0;
      for (let k = 0; k < K; k++) { const below = climatologyAbove(pi * levels[k + 1], ozoneWeight); layerOzone[k] = below - above; above = below; }
    }
    if (clirad || correlated) for (let k = 0; k < K; k++) { if (ozoneProfile || ozone !== 'afgl') layerOzone[k] = ozoneProfile ? ozoneProfile[k] : ozoneCell[i] * ozoneShare[k]; ozoneColumnAmount += layerOzone[k]; }
    const solarPaths = clirad && beam > 0 ? solarGasPaths(i, pi, theta, q, beam / solarConstant) : null;
    const ozoneHeating = clirad ? (solarPaths ? beam * ozoneTaken[K - 1] : 0) : beam * ozoneAbsorption;
    const surfaceEmission = STEFAN_BOLTZMANN * surfaceT * surfaceT * surfaceT * surfaceT;
    const coupled = vaporCoupling > 0 && q !== null;
    const bottom = K - 1;
    const airTemperature = theta[bottom * C + i] * exnerLayer[bottom * C + i];
    const airDensity = pi * sigmaMid[bottom] / (R * airTemperature);
    const exchange = airDensity * exchangeCoefficientAt * Math.max(windSpeed, gustiness);
    const sensible = surfaceLayer ? exchange * (cp * (surfaceT - airTemperature) - cp * thetaV[bottom * C + i] * (exnerLower[bottom * C + i] - exnerLayer[bottom * C + i])) : exchange * cp * (surfaceT - airTemperature);
    const evaporation = qAir === null ? 0 : wetness * Math.max(0, exchange * (saturationHumidity(surfaceT, pi) - qAir));
    budget.surfaceBuoyancy = g / theta[bottom * C + i] * (sensible / (airDensity * cp * exnerLayer[bottom * C + i]) + 0.61 * theta[bottom * C + i] * evaporation / airDensity);
    const mu = beam / solarConstant, lit = !clirad && vaporAbsorption > 0 && q !== null && mu > 0;
    let downPath = 0;
    if (lit) {
      const magnification = 35 / Math.sqrt(1224 * mu * mu + 1);
      let path = 0;
      for (let k = 0; k < K; k++) {
        path += Math.max(0, q[k * C + i]) * pi * dSigma[k] / g * Math.sqrt(sigmaMid[k]) * 0.1 * magnification;
        vaporTaken[k] = vaporAbsorption * waterVaporAbsorptivity(path);
      }
      downPath = path;
    }
    let cloudPath = 0, cloudDepth = 0, columnCover = 0, block = 0, clearColumn = 1, inversionShare = 0, lowBlock = 0, lowClear = 1, lowMaximum = 0, lowPath = 0;
    let previous = 0, cumulative = 0, lowPrevious = 0, lowCumulative = 0;
    if (qAir !== null) {
      const lcl = liftingCondensationLevel(airTemperature, qAir, pi * sigmaMid[bottom], kappa);
      if (lcl) {
        const lower = bottom * C + i, upper = stabilityLayer * C + i, base = Math.max(0, cp * (airTemperature - lcl.temperature) / g);
        const inversion = inversionStrength(theta[upper] - theta[lower], airTemperature, theta[upper] * exnerLayer[upper], (geopotential[upper] - geopotential[lower]) / g - base, cp, R, g, latentHeat);
        inversionShare = Math.min(1, Math.max(0, (inversion - overcastInversion[0]) / (overcastInversion[1] - overcastInversion[0])));
      }
    }
    budget.stratiform = openSea > 0 ? inversionShare : 0;
    const stratiform = cloudCover === 'pdf' && overcastWater !== null ? inversionShare : 0;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      cloudWater[k] = qc ? Math.max(0, qc[idx]) * (pi * dSigma[k] / g) : 0;
      if (cloudMask) cloudWater[k] *= cloudMask.resolved[k];
      cloudPath += cloudWater[k];
      if (graySolar && grayInfrared) { solarDepth[k] = cloudScattering; infrared[k] = cloudAbsorption; } else {
        cloudOptics(theta[idx] * exnerLayer[idx], continental[i] === 1, optics, layerOptics);
        solarDepth[k] = graySolar ? cloudScattering : layerOptics.solar;
        infrared[k] = grayInfrared ? cloudAbsorption : layerOptics.infrared;
      }
      layerCover[k] = 1;
      if (cloudCover === 'pdf' && q !== null && cloudWater[k] > 0 && boundaryCover === 'variance' && (geopotential[idx] - geopotential[bottom * C + i]) / g < mixingDepth) {
        layerCover[k] = Math.min(1, Math.max(coverFloor, turbulentCover(i, k, pi, theta, q, qc, mixingDepth)));
        if (stratiform > 0) layerCover[k] = (1 - stratiform) * layerCover[k] + stratiform * overcastCover(idx, k, pi, theta, q, qc, boundaryCriticalHumidity);
      } else if (cloudCover === 'pdf' && q !== null && cloudWater[k] > 0) {
        const inside = (geopotential[idx] - geopotential[bottom * C + i]) / g < mixedDepth;
        const qs = saturationHumidity(theta[idx] * exnerLayer[idx], pi * sigmaMid[k]), water = Math.max(0, qc[idx]), excess = Math.max(0, q[idx]) + water - qs;
        let rhc = inside ? boundaryCriticalHumidity : criticalHumidity;
        const width = (1 - rhc) * qs;
        if (condensing !== null && condensing.uniform) {
          const pressure = pi * sigmaMid[k], liquidT = theta[idx] * exnerLayer[idx] - latentHeat * water / cp;
          cloudSaturation(liquidT, pressure, condensing.iceSaturation, condensing.liquidTemperature, condensing.iceTemperature, saturated);
          rhc = criticalHumidityAt(pressure, pi, condensing.surfaceCriticalHumidity, condensing.topCriticalHumidity, condensing.criticalExponent);
          layerCover[k] = Math.min(1, Math.max(coverFloor, uniformCover(water, (1 - rhc) * saturated.qs / (1 + latentHeat * saturated.slope / cp))));
        } else layerCover[k] = Math.min(1, Math.max(coverFloor, (excess + width) / (2 * width)));
        if (stratiform > 0) layerCover[k] = (1 - stratiform) * layerCover[k] + stratiform * overcastCover(idx, k, pi, theta, q, qc, rhc);
      }
      const cumulus = cumulusCover !== null ? cumulusCover[idx] * cumulusWater[idx] * (pi * dSigma[k] / g) * (cloudMask ? cloudMask.cumulus[k] : 1) : 0;
      if (cumulus > 0) {
        layerCover[k] = cloudWater[k] > 0 ? Math.max(layerCover[k], cumulusCover[idx]) : cumulusCover[idx];
        cloudWater[k] += cumulus;
        cloudPath += cumulus;
      }
      cloudDepth += solarDepth[k] * cloudWater[k];
      let seen = cloudCover === 'pdf' && cloudWater[k] > 0 ? layerCover[k] * -Math.expm1(-cloudWater[k] / VISIBLE_PATH) : 0;
      if (cumulus > 0 && cloudCover === 'pdf') seen = Math.max(seen, cumulusCover[idx] * -Math.expm1(-cumulus / (cumulusCover[idx] * VISIBLE_PATH)));
      if (seen > 0) block = Math.max(block, seen);
      if (block > 0 && (!(seen > 0) || k === K - 1)) { clearColumn *= 1 - block; columnCover = Math.max(columnCover, block); block = 0; }
      if (exponential) {
        const alpha = k > 0 ? Math.exp(-(geopotential[(k - 1) * C + i] - geopotential[idx]) / (g * decorrelation[i])) : 0;
        cumulative = overlapped(cumulative, previous, seen, alpha);
        previous = seen;
        if (pi * sigmaMid[k] > LOW_CLOUD_PRESSURE) { lowCumulative = overlapped(lowCumulative, lowPrevious, seen, alpha); lowPrevious = seen; }
      }
      if (pi * sigmaMid[k] > LOW_CLOUD_PRESSURE) {
        lowPath += cloudWater[k];
        if (seen > 0) lowBlock = Math.max(lowBlock, seen);
        if (lowBlock > 0 && (!(seen > 0) || k === K - 1)) { lowClear *= 1 - lowBlock; lowMaximum = Math.max(lowMaximum, lowBlock); lowBlock = 0; }
      }
    }
    if (cloudOverlap === 'maximumRandom' && cloudCover === 'pdf') columnCover = 1 - clearColumn;
    if (exponential && cloudCover === 'pdf') columnCover = cumulative;
    if (!(columnCover > 0)) columnCover = 1;
    const sunlit = beam - ozoneHeating, vaporHeating = solarPaths ? beam * (visibleVaporTaken[K - 1] + nearInfraredTaken[K - 1]) : lit ? sunlit * vaporTaken[K - 1] : 0;
    let incident = sunlit - vaporHeating, aerosolHeating = 0, restLoss = 0, visibleLoss = 0, aerosolLoss = 0, ozoneLoss = 0;
    light.share = 0;
    if (scatters || upwardAbsorption) {
      const visible = Math.max(0, visibleFraction * beam - ozoneHeating - (solarPaths ? beam * visibleVaporTaken[K - 1] : 0)), aerosol = aerosolCell[i];
      if (mu > 0 && aerosol > 0) aerosolHeating = -visible * Math.expm1(-(1 - aerosolAlbedo) * aerosol * 35 / Math.sqrt(1224 * mu * mu + 1));
      incident -= aerosolHeating;
      light.share = incident > 0 ? Math.min(1, (visible - aerosolHeating) / incident) : 0;
      light.pressure = pi / REFERENCE_PRESSURE;
      light.aerosol = (1 - aerosolAsymmetry) * aerosolAlbedo * aerosol;
      if (upwardAbsorption && mu > 0) {
        aerosolLoss = aerosol > 0 ? -Math.expm1(-(1 - aerosolAlbedo) * aerosol * DIFFUSE_PATH) : 0;
        ozoneLoss = clirad ? -Math.expm1(-VISIBLE_OZONE * ozoneColumnAmount * DIFFUSE_PATH) : 0;
        visibleLoss = ozoneLoss > 0 ? 1 - (1 - aerosolLoss) * (1 - ozoneLoss) : aerosolLoss;
        const restAfter = sunlit - visible - vaporHeating;
        upwardVapor.fill(0);
        if (solarPaths && restAfter > 0) {
          restLoss = beam * nearInfraredUpward(i, pi, theta, q, solarPaths) / restAfter;
          for (let k = 0; k < K; k++) upwardVapor[k] *= beam / restAfter;
          if (restLoss > 1) { for (let k = 0; k < K; k++) upwardVapor[k] /= restLoss; restLoss = 1; }
        } else if (lit && restAfter > 0) {
          let path = downPath, before = waterVaporAbsorptivity(downPath);
          for (let k = K - 1; k >= 0; k--) {
            path += Math.max(0, q[k * C + i]) * pi * dSigma[k] / g * Math.sqrt(sigmaMid[k]) * 0.1 * DIFFUSE_PATH;
            const through = waterVaporAbsorptivity(path);
            upwardVapor[k] = sunlit * vaporAbsorption * (through - before) / restAfter;
            restLoss += upwardVapor[k];
            before = through;
          }
          if (restLoss > 1) { for (let k = 0; k < K; k++) upwardVapor[k] /= restLoss; restLoss = 1; }
        }
      }
    }
    const inCloud = cloudPath / columnCover;
    shortwave(sky, graySolar ? cloudScattering * inCloud : cloudDepth / columnCover, Math.exp(-cloudSolarAbsorption * inCloud), mu, surfaceAlbedo, diffuseAlbedo);
    if (columnCover < 1) {
      shortwave(clearSky, 0, 1, mu, surfaceAlbedo, diffuseAlbedo);
      for (const key of SHORTWAVE_KEYS) sky[key] = columnCover * sky[key] + (1 - columnCover) * clearSky[key];
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
        deckLight.path = cloudPath; deckLight.depth = cloudDepth; deckLight.layer = cloudWater[stratusLayer]; deckLight.clear = clearShare * cloudWater[stratusLayer];
        if (shadowDeck(i, pi, theta, q, qc, mixedDepth, sensible, evaporation, dt, stratusSolar ? deckAbsorption : null)) {
          fraction = budget.mlmCover * openSea;
          if (fraction > 0) deck = Math.min(stratusWaterMax, budget.mlmWater);
        }
      }
      if (deck <= 0) fraction = 0;
    }
    if (cloudMask) { fraction = cloudMask.fraction; deck = cloudMask.deck; }
    for (let k = 0; k < K; k++) {
      const mass = pi * dSigma[k] / g;
      emissivity[k] = 1 - Math.exp(coupled ? -vaporCoupling * Math.max(0, q[k * C + i]) * mass : -tau0 * shape[k]);
      const water = cloudWater[k];
      cloudEmissivity[k] = water > 0 ? layerCover[k] * (1 - Math.exp(-infrared[k] * water / layerCover[k])) : 0;
      if (deck > 0 && k === stratusLayer) cloudEmissivity[k] = fraction * (1 - Math.exp(-infrared[k] * (water + deck))) + (1 - fraction) * cloudEmissivity[k];
      const clear = 1 - cloudEmissivity[k];
      vaporEmissivity[k] = 1 - (1 - emissivity[k]) * clear;
      mixedEmissivity[k] = 1 - (1 - gasEmissivity[k]) * clear;
      temperature[k] = theta[k * C + i] * exnerLayer[k * C + i];
      netFlux[k] = solarPaths ? beam * (ozoneTaken[k] - (k > 0 ? ozoneTaken[k - 1] : 0)) : ozoneHeating * ozoneFraction[k];
      if (aerosolHeating > 0) netFlux[k] += aerosolHeating * aerosolFraction[k];
    }
    if (lit) {
      let taken = 0;
      for (let k = 0; k < K; k++) {
        const through = vaporTaken[k];
        netFlux[k] += sunlit * (through - taken);
        taken = through;
      }
    } else if (solarPaths) {
      let taken = 0;
      for (let k = 0; k < K; k++) {
        const through = visibleVaporTaken[k] + nearInfraredTaken[k];
        netFlux[k] += beam * (through - taken);
        taken = through;
      }
    }
    if (deck > 0) {
      const total = cloudPath + deck;
      shortwave(decked, graySolar ? cloudScattering * total : cloudDepth + solarDepth[stratusLayer] * deck, Math.exp(-cloudSolarAbsorption * total), mu, surfaceAlbedo, diffuseAlbedo);
      deckShare = fraction * incident * decked.cloud / total;
      clearShare *= 1 - fraction;
      sky.absorbed = fraction * decked.absorbed + (1 - fraction) * sky.absorbed;
      sky.down = fraction * decked.down + (1 - fraction) * sky.down;
      sky.direct = fraction * decked.direct + (1 - fraction) * sky.direct;
      sky.reflectance = fraction * decked.reflectance + (1 - fraction) * sky.reflectance;
      sky.cloud = fraction * decked.cloud + (1 - fraction) * sky.cloud;
      sky.visibleEscape = fraction * decked.visibleEscape + (1 - fraction) * sky.visibleEscape;
      sky.restEscape = fraction * decked.restEscape + (1 - fraction) * sky.restEscape;
      sky.visibleReflectance = fraction * decked.visibleReflectance + (1 - fraction) * sky.visibleReflectance;
    }
    let cloudHeating = 0;
    if (incident > 0 && sky.cloud > 0) {
      for (let k = 0; k < K; k++) {
        const share = clearShare * cloudWater[k] + deckShare * (k === stratusLayer ? cloudWater[k] + deck : cloudWater[k]);
        netFlux[k] += share;
        cloudHeating += share;
      }
    }
    let upwardHeating = 0, upwardAerosol = 0, upwardOzoneHeating = 0;
    if (restLoss > 0 || visibleLoss > 0) {
      const rest = incident * sky.restEscape;
      for (let k = 0; k < K; k++) netFlux[k] += rest * upwardVapor[k];
      upwardAerosol = incident * sky.visibleEscape * aerosolLoss;
      if (upwardAerosol > 0) for (let k = 0; k < K; k++) netFlux[k] += upwardAerosol * aerosolFraction[k];
      upwardOzoneHeating = incident * (sky.visibleEscape * (1 - aerosolLoss) + sky.visibleReflectance) * ozoneLoss;
      if (upwardOzoneHeating > 0) for (let k = 0; k < K; k++) netFlux[k] += upwardOzoneHeating * layerOzone[k] / ozoneColumnAmount;
      upwardHeating = rest * restLoss + upwardAerosol + upwardOzoneHeating;
    }
    const absorbedSolar = incident * sky.absorbed;
    beforeBands.set(netFlux);
    let outgoing = 0, back = 0, clearOutgoing = 0;
    if (correlated) {
      for (let k = 0; k < K; k++) {
        const idx = k * C + i, dry = Math.max(0, 1 - Math.max(0, q[idx]));
        layerPaths(pathRow, pi * sigmaMid[k], pi * dSigma[k] / g, temperature[k], q[idx], layerOzone[k] * OZONE_CM_ATM, wellMixed[0] * dry, wellMixed[1] * dry, wellMixed[2] * dry);
        for (let j = 0; j < 6; j++) paths[j][k] = pathRow[j];
      }
      if (chained) {
        for (let k = 0; k < K; k++) {
          const water = cloudWater[k], decked = deck > 0 && k === stratusLayer;
          chainCover[k] = water > 0 && !decked ? layerCover[k] : 0;
          chainCloud[k] = chainCover[k] > 0 ? -Math.expm1(-infrared[k] * water / chainCover[k]) : 0;
          chainOutside[k] = decked ? cloudEmissivity[k] : 0;
          clearShareInverse[k] = chainCover[k] < 1 ? 1 / (1 - chainCover[k]) : 0;
          cloudShareInverse[k] = chainCover[k] > 0 ? 1 / chainCover[k] : 0;
          if (k === 0) continue;
          const a = chainCover[k - 1], b = chainCover[k], alpha = Math.exp(-(geopotential[(k - 1) * C + i] - geopotential[k * C + i]) / (g * decorrelation[i]));
          const both = a + b - (alpha * Math.max(a, b) + (1 - alpha) * (a + b - a * b));
          joint[0][k] = 1 - a - b + both; joint[1][k] = b - both; joint[2][k] = a - both; joint[3][k] = both;
        }
      }
      for (const row of LONGWAVE_TABLE.points) {
        for (let k = 0; k < K; k++) {
          const tau = LONGWAVE_CONSTANTS.diffusivity * (row[0] * paths[0][k] + row[1] * paths[1][k] + row[2] * paths[2][k] + row[3] * paths[3][k] + row[4] * paths[4][k] + row[5] * paths[5][k]);
          gasEmissivityG[k] = -Math.expm1(-tau);
          totalEmissivity[k] = 1 - (1 - gasEmissivityG[k]) * (1 - cloudEmissivity[k]);
          planck[k] = planckShare(row, temperature[k]) * STEFAN_BOLTZMANN * temperature[k] ** 4;
        }
        const surfaceUp = planckShare(row, surfaceT) * surfaceEmission;
        const [up, down] = chained ? chainPlanck(planck, gasEmissivityG, surfaceUp) : bandPlanck(planck, totalEmissivity, surfaceUp);
        outgoing += up; back += down;
        if (clearSkyPass) clearOutgoing += upwardPlanck(planck, gasEmissivityG, surfaceUp);
      }
    } else {
      const [outVapor, backVapor] = band(vaporFraction, vaporEmissivity, surfaceEmission);
      const [outGas, backGas] = band(gasFraction, mixedEmissivity, surfaceEmission);
      const [outWindow, backWindow] = band(window, cloudEmissivity, surfaceEmission);
      outgoing = outVapor + outGas + outWindow;
      back = backVapor + backGas + backWindow;
      if (clearSkyPass) clearOutgoing = upward(vaporFraction, emissivity, surfaceEmission) + upward(gasFraction, gasEmissivity, surfaceEmission) + window * surfaceEmission;
    }
    for (let k = 0; k < K; k++) longwave[k * C + i] = netFlux[k] - beforeBands[k];
    if (clearSkyPass) {
      shortwave(clearSky, 0, 1, mu, surfaceAlbedo, diffuseAlbedo);
      const clearUpward = restLoss > 0 || visibleLoss > 0 ? incident * clearSky.restEscape * restLoss + incident * clearSky.visibleEscape * aerosolLoss + (ozoneLoss > 0 ? incident * (clearSky.visibleEscape * (1 - aerosolLoss) + clearSky.visibleReflectance) * ozoneLoss : 0) : 0;
      budget.clearAbsorbedSolar = incident * clearSky.absorbed + ozoneHeating + vaporHeating + aerosolHeating + clearUpward;
      budget.clearOutgoingLongwave = clearOutgoing;
    }
    netFlux[bottom] += sensible;
    const net = absorbedSolar - surfaceEmission + back - sensible - latentHeat * evaporation;
    budget.absorbedSolar = absorbedSolar + ozoneHeating + vaporHeating + aerosolHeating + cloudHeating + upwardHeating;
    budget.atmosphereSolar = ozoneHeating + vaporHeating + aerosolHeating + cloudHeating + upwardHeating;
    budget.aerosolSolar = aerosolHeating + upwardAerosol;
    budget.ozoneSolar = ozoneHeating + upwardOzoneHeating;
    budget.vaporSolar = solarPaths ? beam * gasSplit.vapour : vaporHeating;
    budget.oxygenSolar = solarPaths ? beam * gasSplit.oxygen : 0;
    budget.carbonDioxideSolar = solarPaths ? beam * gasSplit.co2 : 0;
    budget.upwardGasSolar = (restLoss > 0 ? incident * sky.restEscape * restLoss : 0) + upwardOzoneHeating;
    budget.downwardLongwave = back;
    budget.cloudSolar = cloudHeating;
    budget.outgoingLongwave = outgoing;
    budget.sensibleHeat = sensible;
    budget.evaporation = evaporation;
    budget.potentialEvaporation = qAir === null ? 0 : referenceEvaporation(absorbedSolar - surfaceEmission + back, airTemperature, qAir, pi, airDensity, (referenceCoefficients ? referenceCoefficients[i] : exchangeCoefficientAt) * Math.max(windSpeed, gustiness), cp, latentHeat);
    budget.surfaceFlux = net;
    budget.insolation = beam;
    budget.reflectedSolar = incident - absorbedSolar - cloudHeating - upwardHeating;
    budget.surfaceShortwave = incident * sky.down;
    budget.surfaceDirect = incident * sky.direct;
    budget.cloudReflectance = sky.reflectance;
    budget.cloudCover = cloudPath > 0 ? columnCover : 0;
    const lowResolved = cloudCover !== 'pdf' ? (lowPath > 0 ? 1 : 0) : exponential ? lowCumulative : cloudOverlap === 'maximumRandom' ? 1 - lowClear : lowMaximum;
    budget.lowCover = 1 - (1 - lowResolved) * (1 - fraction);
    budget.lowWater = lowPath + fraction * deck;
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
   * mixed-layer deck advances by. Each cell's absorbed, atmosphere-absorbed,
   * incoming and reflected sunlight and its outgoing longwave add to
   * `summed` (shared, one array per name in SUMMED) until restartSums;
   * readMeans(steps) turns the sums of that many steps into the per-cell
   * means meanAbsorbedSolar and meanOutgoingLongwave and the albedo
   * meanPlanetaryAlbedo, reflected over incoming summed (0 where no sun
   * rose). With clearSkyPass each column also finds its clear-sky
   * absorbed sunlight and outgoing longwave, the same column with no
   * resolved cloud, deck or cumulus (the clear two-stream; the vapour and
   * gas bands' upward pass without cloud emissivity and the open window),
   * summed in clearAbsorbedSolar and clearOutgoingLongwave, and readMeans
   * also gives meanShortwaveCloudEffect (absorbed less clear-sky absorbed)
   * and meanLongwaveCloudEffect (clear-sky less all-sky outgoing).
   * useCloudMask(mask) (this engine only; null clears it) weighs each
   * column's resolved water by mask.resolved[k] and its cumulus by
   * mask.cumulus[k] and gives it the deck mask.fraction and mask.deck in
   * place of the deck it diagnoses, for scripts/cloudClasses.mjs.
   * With surfaceLayer the sensible heat flux is that of the dry static
   * energy, ρ C_H U (c_p T_s − c_p T − g z) with z the lowest layer's
   * height (IFS Cy47r3 Part IV eq. 8.6), and the reference evaporation
   * takes referenceCoefficients. Each cell's surface buoyancy flux from
   * the sensible and latent fluxes, g/θ (H/(ρ c_p Π) + 0.61 θ E/ρ), goes
   * to `surfaceBuoyancy`.
   */
  function apply(state, out, windSpeed, totals, iFrom = 0, iTo = C, surfaceAlbedo = null, diffuseAlbedo = null, wetness = null, openSea = null, depth = null, dt = 0) {
    const [pi, theta, , surfaceT] = state;
    const [, dTheta] = out;
    const q = state[4] ?? null, dQ = out[4] ?? null, qc = state[5] ?? null;
    const bottom = (K - 1) * C;
    if (totals) for (const name of ['absorbedSolar', 'atmosphereSolar', 'outgoingLongwave', 'sensibleHeat', 'evaporation', 'insolation', 'reflectedSolar']) totals[name] = 0;
    for (let i = iFrom; i < iTo; i++) {
      surfaceFlux[i] = column(i, pi[i], theta, surfaceT[i], windSpeed[i], tauCell[i], insolation(i), q && dQ ? q[bottom + i] : null, q && dQ ? q : null, q && dQ ? qc : null, surfaceAlbedo ? surfaceAlbedo[i] : albedo, diffuseAlbedo ? diffuseAlbedo[i] : surfaceAlbedo ? surfaceAlbedo[i] : albedo, wetness ? wetness[i] : 1, exchangeCoefficients ? exchangeCoefficients[i] : exchangeCoefficient, openSea ? openSea[i] : 0, depth ? depth[i] - geopotential[bottom + i] / g : 0, dt, boundaryTop ? boundaryTop[i] - geopotential[bottom + i] / g : 0);
      outgoing[i] = budget.outgoingLongwave;
      stratusPath[i] = budget.stratus;
      stratusCover[i] = budget.stratusFraction;
      stabilityIndex[i] = budget.stabilityIndex;
      mlmCover[i] = budget.mlmCover;
      mlmWater[i] = budget.mlmWater;
      mlmEntrainment[i] = budget.mlmEntrainment;
      mlmTop[i] = budget.mlmTop;
      stratiformShare[i] = budget.stratiform;
      lowCover[i] = budget.lowCover;
      lowWater[i] = budget.lowWater;
      evaporation[i] = budget.evaporation;
      surfaceBuoyancy[i] = budget.surfaceBuoyancy;
      potentialEvaporation[i] = budget.potentialEvaporation;
      surfaceShortwave[i] = budget.surfaceShortwave;
      surfaceDirect[i] = budget.surfaceDirect;
      for (let k = 0; k < K; k++) {
        const massPerArea = pi[i] * dSigma[k] / g;
        dTheta[k * C + i] += netFlux[k] / (cp * massPerArea) / exnerLayer[k * C + i];
      }
      if (q && dQ) dQ[bottom + i] += budget.evaporation * g / (pi[i] * dSigma[K - 1]);
      summed.absorbedSolar[i] += budget.absorbedSolar;
      summed.atmosphereSolar[i] += budget.atmosphereSolar;
      summed.outgoingLongwave[i] += budget.outgoingLongwave;
      summed.insolation[i] += budget.insolation;
      summed.reflectedSolar[i] += budget.reflectedSolar;
      if (clearSkyPass) {
        summed.clearAbsorbedSolar[i] += budget.clearAbsorbedSolar;
        summed.clearOutgoingLongwave[i] += budget.clearOutgoingLongwave;
      }
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

  function readMeans(steps) {
    for (let i = 0; i < C; i++) {
      meanAbsorbedSolar[i] = summed.absorbedSolar[i] / steps;
      meanOutgoingLongwave[i] = summed.outgoingLongwave[i] / steps;
      meanPlanetaryAlbedo[i] = summed.insolation[i] > 0 ? summed.reflectedSolar[i] / summed.insolation[i] : 0;
      if (clearSkyPass) {
        meanShortwaveCloudEffect[i] = (summed.absorbedSolar[i] - summed.clearAbsorbedSolar[i]) / steps;
        meanLongwaveCloudEffect[i] = (summed.clearOutgoingLongwave[i] - summed.outgoingLongwave[i]) / steps;
      }
    }
  }

  function restartSums() {
    for (const name of [...SUMMED, ...CLEAR_SUMMED]) summed[name].fill(0);
  }

  function useBoundaryLayer(regime, mixingTop, buoyancyFlux = null) {
    boundaryRegime = regime;
    boundaryTop = mixingTop;
    boundaryBuoyancy = buoyancyFlux;
  }

  function useCloudMask(mask) {
    cloudMask = mask;
  }

  function useCondensation(settings) {
    condensing = settings;
  }

  function useCumulus(cover, water) {
    const on = cumulusCloud && cloudCover === 'pdf' && cover && water;
    cumulusCover = on ? cover : null;
    cumulusWater = on ? water : null;
  }

  const deckGates = { subsidenceSmoothing, subsidenceMemory, stratusSubsidence, minimumInversion, ceilingInversion: ceilingJump, gateMemory, deckRest, cumulusCeiling, deckRegime, deckBypass, stratusWaterMax };
  return { setTime, sun, cosZenith, insolation, column, apply, readMeans, restartSums, summed, clearSkyPass, useCloudMask, useCondensation, solarDepth, infrared, meanAbsorbedSolar, meanOutgoingLongwave, meanPlanetaryAlbedo, meanShortwaveCloudEffect, meanLongwaveCloudEffect, useCumulus, useBoundaryLayer, longwave, layerCover, lowCover, lowWater, deckGates, layerFlux: netFlux, surfaceFlux, outgoing, surfaceShortwave, surfaceDirect, evaporation, potentialEvaporation, surfaceBuoyancy, stratus: stratusPath, stratusFraction: stratusCover, stabilityIndex, mlmCover, mlmWater, mlmEntrainment, mlmSubsidence, mlmHeight, mlmGate, mlmTop, stratiform: stratiformShare, stratusLayer, stabilityLayer, budget, emissivity, opticalDepth, ozoneFraction, shared: { outgoing: outgoingBuffer, surfaceShortwave: shortwaveBuffer, evaporation: evaporationBuffer, surfaceBuoyancy: surfaceBuoyancyBuffer, stratus: stratusBuffer, stratusFraction: coverBuffer, stabilityIndex: indexBuffer, mlmCover: mlmCoverBuffer, mlmWater: mlmWaterBuffer, mlmEntrainment: mlmEntrainmentBuffer, mlmSubsidence: mlmSubsidenceBuffer, mlmHeight: mlmHeightBuffer, mlmGate: mlmGateBuffer, mlmTop: mlmTopBuffer, stratiform: stratiformBuffer, longwave: longwaveBuffer, summed: summedBuffers } };
}
