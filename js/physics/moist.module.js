import { R_DRY, VIRTUAL_FACTOR } from '../dynamics/sigmaCore.module.js';
import { cellVector } from '../dynamics/operators.module.js';

export const LATENT_HEAT = 2.5e6;
export const EPSILON = 0.622;
export const R_VAPOR = R_DRY / EPSILON;
export const CLEAR_AIR = 1e-7;
export const DECK_OPEN = 0.5;
export const DECK_CLOSED = 0.6;
export const CUMULUS_FLOOR = 1e-6;
export const MAXIMUM_SURFACE_PRESSURE = 110000;
export const DEEP_REFERENCE = 1e5;
export const COUPLED_REGIME = 3;
export const BECHTOLD = { resolution: 1.66, reference: 125e3, shortest: 720, longest: 10800, boundaryWind: 2, temperatureScale: 1 };
export const SUBCLOUD_LAYERS = 12;
export const DEEP_CLOUD_DEPTH = 200e2;
export const IFS_ENTRAINMENT = { entrainment: 1.75e-3, humidity: 1.3, detrainment: 0.75e-4, detrainmentHumidity: 1.6, scaleExponent: 3, drag: 1 + 1.875 * 0.506 };
export const TEST_PARCEL = { entrainment: 0.4, removal: 0.5 };
export const IFS_PRECIPITATION = { conversion: 1.4e-3, liquidFactor: 1.3, critical: 5e-4, seaThreshold: 3e-4, landThreshold: 5e-4, bergeron: 268.16, ice: 250.16, speed: 10, velocityScale: 0.75 };
export const SOURCE_EXCESS = { coefficient: 1.5, temperature: 3, humidity: 2e-3, scale: 1.2, layer: 1.5, karman: 0.4, friction: 0.1, virtual: 0.61 };

/*
 * The convective-scale velocity of the IFS's test parcel at the lowest
 * model level (Cy43r1 eq. 6.20), w* = 1.2 (u*³ + 1.5 g z κ / T (J_s / (ρ c_p)
 * + 0.61 T J_q / ρ))^⅓ with u* 0.1 m/s, from the surface sensible heat
 * `sensible` (W/m²) and evaporation (kg/m²/s), the lowest layer's density,
 * temperature and height above the ground; 0 unless the surface buoyancy
 * flux is upward, where the IFS gives the parcel no excess.
 */
export function surfaceLayerVelocity(sensible, evaporation, density, temperature, height, cp, g) {
  const flux = (sensible / cp + SOURCE_EXCESS.virtual * temperature * evaporation) / density;
  if (!(flux > 0)) return 0;
  return SOURCE_EXCESS.scale * Math.cbrt(SOURCE_EXCESS.friction ** 3 + SOURCE_EXCESS.layer * g * height * SOURCE_EXCESS.karman / temperature * flux);
}

export function saturationVaporPressure(T) {
  return 611.2 * Math.exp(17.67 * (T - 273.15) / (T - 29.65));
}

export function saturationHumidity(T, p) {
  const es = saturationVaporPressure(T);
  const dry = p - (1 - EPSILON) * es;
  return dry > 0 ? EPSILON * es / dry : 1;
}

export const FUSION_HEAT = 3.34e5;
export const LIQUID_TEMPERATURE = 273.15, ICE_TEMPERATURE = 235.15;

export function iceVaporPressure(T) {
  return 611.21 * Math.exp(22.587 * (T - 273.16) / (T + 0.7));
}

export function liquidFraction(T, liquidTemperature = LIQUID_TEMPERATURE, iceTemperature = ICE_TEMPERATURE) {
  return Math.min(1, Math.max(0, (T - iceTemperature) / (liquidTemperature - iceTemperature)));
}

/*
 * The saturation humidity of cloud at temperature T and pressure p into
 * out.qs, with out.slope its Clausius–Clapeyron derivative q_s L/(R_v T²)
 * and out.liquid the liquid share. Over liquid water (Bolton's) unless
 * `ice`; with it the vapour pressure is the liquid share's mix
 * α e_w + (1 − α) e_i, α linear in T from 0 at iceTemperature to 1 at
 * liquidTemperature (the radiation's phase ramp), e_i the IFS form, and
 * the slope's latent heat L + (1 − α) FUSION_HEAT.
 */
export function cloudSaturation(T, p, ice, liquidTemperature, iceTemperature, out) {
  const alpha = ice ? liquidFraction(T, liquidTemperature, iceTemperature) : 1;
  const es = alpha < 1 ? alpha * saturationVaporPressure(T) + (1 - alpha) * iceVaporPressure(T) : saturationVaporPressure(T);
  const dry = p - (1 - EPSILON) * es;
  out.qs = dry > 0 ? EPSILON * es / dry : 1;
  out.slope = out.qs * (LATENT_HEAT + (1 - alpha) * FUSION_HEAT) / (R_VAPOR * T * T);
  out.liquid = alpha;
  return out;
}

/*
 * The critical relative humidity at pressure p under surface pressure ps,
 * top + (surface − top) exp(1 − (ps/p)^exponent) (ECHAM's profile).
 */
export function criticalHumidityAt(p, ps, surface, top, exponent) {
  return top + (surface - top) * Math.exp(1 - Math.pow(ps / p, exponent));
}

/*
 * The uniform total-water distribution of half-width b (in saturation
 * deficit, b > 0) about a mean deficit Q: its condensate, and the cover
 * that a condensate qc implies, sqrt(qc / b) up to 1.
 */
export function uniformCondensate(Q, b) {
  return Q >= b ? Q : Q > -b ? (Q + b) * (Q + b) / (4 * b) : 0;
}

export function uniformCover(qc, b) {
  return qc >= b ? 1 : qc > 0 ? Math.sqrt(qc / b) : 0;
}

/*
 * Lifting condensation level of a parcel (T, q, p) by Bolton (1980):
 * the dew point from the vapour pressure, the LCL temperature from his
 * eq. 15, and the pressure along the dry adiabat. Returns null when the
 * parcel is already saturated or has no vapour.
 */
export function liftingCondensationLevel(T, q, p, kappa) {
  if (q <= 0) return null;
  const e = q * p / (EPSILON + (1 - EPSILON) * q);
  const y = Math.log(e / 611.2);
  const dewPoint = (273.15 * 17.67 - 29.65 * y) / (17.67 - y);
  if (dewPoint >= T) return { temperature: T, pressure: p };
  const temperature = 1 / (1 / (dewPoint - 56) + Math.log(T / dewPoint) / 800) + 56;
  return { temperature, pressure: p * Math.pow(temperature / T, 1 / kappa) };
}

/*
 * Moist physics for the sigma core, applied to the state after each
 * step in the order: saturation adjustment, convection, rain, filler.
 *
 * Saturation adjustment: supersaturated vapour condenses into cloud
 * water and cloud water evaporates into subsaturated air, latent heat to
 * the layer. With `iceSaturation` (true) saturation is cloudSaturation's
 * mix over water and ice along the phase ramp from `iceTemperature`
 * (235.15 K) to `liquidTemperature` (273.15 K), the radiation's; the
 * model's condensate carries no enthalpy, every phase change takes the
 * latent heat L of vaporisation, and the fusion heat of what falls as snow
 * is released at the surface. With `condensation` 'uniform' (the
 * default; 'saturation' adjusts every layer to its own saturation) a layer
 * above the moist boundary layer's mixing top (every layer without it)
 * holds the condensate of a uniform distribution of total water about its
 * mean (LeTreut and Li 1991; with fixed width the cover of Sundqvist et
 * al. 1989): in the saturation deficit of the layer's θ_l and q_t,
 * Q = a (q_t − q_s(T_l)) with a = 1/(1 + (L/c_p) dq_s/dT), half-width
 * b = a (1 − RH_c) q_s(T_l), the condensate is Q above b, (Q + b)²/(4b)
 * between −b and b and none below, the cover (Q + b)/(2b) = sqrt(q_c/b),
 * in one linearised step about T_l. RH_c follows ECHAM's profile
 * criticalHumidityAt with `surfaceCriticalHumidity` 0.975,
 * `topCriticalHumidity` 0.75 and `criticalExponent` 2 (ECHAM6 at T63,
 * mo_cloud's crs, crt and nex). The layers below the mixing top, whose
 * cover is the boundary layer's variance cover, adjust to saturation.
 * With `iceNucleation` (false; under iceSaturation) a cloud-free layer (at
 * most CLEAR_AIR of condensate) colder than iceTemperature forms cloud only
 * where its distribution exceeds q_ref = min(q_sw, RH_homo q_si),
 * RH_homo = 2.583 − T/207.8 the homogeneous nucleation threshold of Kärcher
 * and Lohmann (2002), so clear air may be supersaturated over ice; the part
 * above q_ref deposits to ice saturation, and from then on the layer holds
 * the distribution's condensate about q_si (the scheme of Tompkins et al.
 * 2007 in the IFS, Cy43r1 documentation §7.2.4c).
 *
 * Convection is one entraining mass-flux plume with a shallow and a
 * deep branch. Both run in proportion to the deck's opening, 1 where
 * the mixed-layer deck's gate `deckGate` (radiation.mlmGate) is at most
 * DECK_OPEN, falling linearly to 0 at DECK_CLOSED (`deckVeto` false:
 * 1 everywhere), and with `coupledVeto` not at all where the moist
 * boundary layer (`boundaryRegime`) is a coupled stratocumulus-topped
 * layer, COUPLED_REGIME.
 *
 * The shallow cumulus mass flux, a bulk entraining plume after
 * Bretherton, McCaa and Grenier (2004), simplified. Its source is the
 * mass-weighted mean liquid-water static energy s_l = c_p T + g z − L q_c
 * and total water q_t of the layers whose lower interfaces lie below
 * the boundary-layer top, or whose midpoints lie within
 * `cumulusSourceDepth` of the surface, and below `shallowTop`. The plume
 * leaves the source's top interface with its air, rises without mixing
 * to the Bolton LCL of that air and above it entrains the layer's s_l
 * and q_t at `cumulusEntrainment` and detrains at `cumulusDetrainment`
 * (per metre), its mass flux growing by exp((ε − δ) Δz) through each
 * layer. At each layer's midpoint the plume's temperature and condensate
 * follow from its s_l and q_t, and its buoyancy is its virtual
 * temperature against the air's, condensate loading included in both.
 * Its inhibition is the negative buoyant energy from the source's top to
 * the first layer where it is saturated; from there on the plume ends in
 * the first layer where it is not buoyant, or in the last layer below
 * `shallowTop`, and detrains all it carries there, `cumulusOvershoot` of
 * the mass flux entering that layer and the rest in the layer below. The
 * base mass flux is ρ_LCL c w exp(−inhibition / w²), c the
 * `cumulusClosure`, w the larger of (B₀ h)^⅓ and `cumulusFriction` u*,
 * B₀ the boundary layer's surface buoyancy flux and h its depth above the
 * lowest layer, times the deck's opening, and zero where B₀ ≤ 0, below
 * CUMULUS_FLOOR (kg/m²/s), where
 * the LCL lies above `shallowTop` or where the plume does not saturate
 * below it; at most `cumulusBoundaryLoss` of the source
 * layers' mass leaves in a step. Inside the source the flux draws on
 * each layer in proportion to its mass. The layers change in flux form:
 * the flux M (X_u − X_above) of s_l and q_t through each interface, X_u
 * the plume's and X_above the layer above's, and each layer's X by
 * g Δt/Δp times the flux through its lower interface less that through
 * its upper one, the mass flux scaled down where needed so that
 * M g Δt/Δp ≤ 1 in every layer; the change in s_l goes to the
 * temperature, that in q_t to the vapour, and a column the plume moved
 * is saturation-adjusted again, so detrained condensate evaporates into
 * dry layers and stays as cloud water in saturated ones. Column enthalpy
 * and water are exact. With `cumulusRain` (kg/kg; null: none) the
 * plume's condensate above it rains at each interface, convective rain.
 * `cumulusCover` holds each plume layer's cumulus fraction, its mean
 * mass flux over ρ `cumulusUpdraft`, and `cumulusWater` the plume's
 * condensate at its midpoint, for the radiation; `cumulusBaseFlux` the
 * base mass flux (kg/m²/s) and `cumulusTop` the pressure of the plume's
 * top interface.
 *
 * The deep plume leaves the layers whose midpoints lie within
 * `cumulusSourceDepth` of the surface with their mean s_l and q_t
 * (`plumeSourceDepth` 'surface50', the IFS's 50 hPa test parcel, Cy43r1
 * §6.5; 'boundaryLayer': the shallow plume's source layers; `plumeSource`
 * 'lowest': the lowest layer's air), with 'surface50' raised by the IFS
 * excess (eq. 6.19) ΔT = min(3 K, 1.5 J_s / (ρ c_p w*)) and
 * Δq = min(2 g/kg, 1.5 J_q / (ρ L w*)), J_s and J_q the cell's surface
 * sensible and latent fluxes (`surfaceSensible`, `surfaceEvaporation`; none
 * without them), ρ the lowest layer's density and w* with `excessVelocity`
 * 'surfaceLayer' (the default) the IFS's own (eq. 6.20,
 * surfaceLayerVelocity), each part at least 0 and none under a downward
 * buoyancy flux as in the IFS's test parcel, or with 'convective' the
 * shallow closure's velocity, and rises unmixed to their LCL. From the first layer above it, its vertical
 * velocity follows d(w²)/dz = 2 a B − 2 b ε w² across each layer exactly
 * for constant B and ε, from `plumeVelocity` w0 at the layer's base, a
 * `plumeAcceleration`, b `plumeDrag`, B its buoyancy at the layer's
 * midpoint (virtual temperature with condensate loading, as above, g ΔTv
 * / Tv), and with `plumeEntrainmentLaw` 'gregory' its fractional
 * entrainment is ε = max(`plumeEntrainmentFloor`, c_ε B / w²) with c_ε
 * `plumeEntrainment`, B the layer below's buoyancy and w² at the layer's
 * base (Gregory 2001), mixing s_l and q_t toward the layer's at exp(−ε Δz).
 * With 'ifs' (the default; IFS Cy43r1 eqs 6.7, 6.8, 6.10 and 6.12, Bechtold
 * et al. 2008, the IFS_ENTRAINMENT constants) ε = 1.75e-3 /m (1.3 − RH)
 * (q_s(T)/q_s(T_base))³ where the layer below is buoyant and 0 elsewhere,
 * RH the layer's relative humidity (at most 1) and q_s its saturation as
 * the condensation takes it, T_base the first cloud layer's temperature;
 * the turbulent detrainment δ = 0.75e-4 /m (1.6 − RH); in the w² equation b
 * is 1 + βC_d = 1 + 1.875 · 0.506 and the mixing rate ε, or δ where ε is 0;
 * and the mass flux grows by exp((ε − δ) Δz) through each cloud layer whose
 * plume is buoyant, and through the others falls by exp(−δ Δz) times the
 * organised detrainment's min(1, (1.6 − RH) √(w²_above / w²_below)) to zero
 * at the top interface. At each layer's upper interface the condensate
 * above `plumeRainThreshold` rains at the fraction 1 − exp(−c0 Δz), c0
 * `plumeRainRate` (Zhang and McFarlane 1995), q_t falling and s_l rising
 * by L times it; with `plumeConversion` 'sundqvist' (the default;
 * 'zhangMcFarlane' the fraction above) the IFS's updraught conversion
 * instead (Cy43r1 eqs 6.38–6.40 after Sundqvist 1978, IFS_PRECIPITATION):
 * where the plume's condensate l at the interface exceeds 0.3 g/kg over sea
 * and 0.5 g/kg over `land`, l (1 − exp(−a Δz)) rains with
 * a = c0 / (0.75 w) (1 − exp(−(l/l_crit)²)), c0 = 1.4e-3 /s (1.3 α + 1 − α),
 * l_crit 0.5 g/kg, w the plume's speed there within plumeVelocity and 10
 * m/s, and below 268.16 K c0 times and l_crit over 1 + 0.5 √min(268.16 − T,
 * 18) (the Bergeron–Findeisen factor). With `plumePhase` 'mixed' (the default; 'liquid' saturates
 * the plumes over liquid water) both plumes and the test parcel below carry
 * s_li = c_p T + g z − L q_l − (L + L_f) q_i, saturated over cloudSaturation's
 * mix with the condensate's ice share 1 − α(T) on the radiation's phase ramp
 * (IFS Cy43r1 §6.6.2, with the model's linear ramp in place of the IFS's
 * quadratic one); a rain's ice share takes L_f more to s_li, falls as its own
 * stream and melts in the first layer below it at liquidTemperature or
 * warmer, L_f times it taken from that layer (§6.6.6, there a relaxation over
 * a few layers); below cloud base the evaporation of the stream's ice share
 * takes L + L_f, and what reaches the ground frozen is `convectiveSnow`, to
 * which the surface adds no fusion heat (and from which it takes L_f where the
 * air at the ground is not freezing); the downdraft evaporates only the
 * liquid and melted rain. The flux form carries s_li, so condensate frozen in
 * the plume that detrains returns to the environment's L-only convention in
 * the detraining layer, and column c_p T + L q stays exact. The plume ends in the layer where w² falls to zero, or in
 * the top layer but one. With `convectionType` 'testParcel' (the default)
 * the IFS's first-guess deep updraught types the column first (Cy43r1 §6.4,
 * eqs 6.18–6.21): a test parcel of the deep source's s_l and q_t leaves the
 * source's top interface at `plumeVelocity` and mixes toward each layer's
 * at ε = 0.4 · 1.75e-3 /m (q_s(T)/q_s(T_lowest))³ (TEST_PARCEL, q_s as the
 * condensation takes it, T_lowest the lowest layer's), below its cloud as
 * in it, its w² following the IFS form above with the mixing rate ε, and at
 * each layer's upper interface half its condensate is removed (L times it
 * to s_l); its cloud runs from the lower interface of its first cloudy
 * layer to the height inside a layer where its w² vanishes, exact for that
 * layer's buoyancy and mixing, at the pressure interpolated in ln p. A
 * column whose test cloud is deeper than DEEP_CLOUD_DEPTH (200 hPa) is deep
 * and runs the deep plume; else, or where the deep plume is not cloudy, the
 * shallow plume alone. With 'cloudDepth' the plume's own cloud, from its
 * base interface to its top interface, deeper than DEEP_CLOUD_DEPTH makes
 * the column deep. Under either a column convects as one type only: where
 * the deep plume runs the shallow one does not, whatever `plumeClosure`
 * says, and where its closure gives no flux the shallow plume runs alone
 * (IFS Cy43r1 §6.4 and §6.4.2: a column is deep convective if its test
 * parcel's cloud is thicker than 200 hPa, otherwise shallow, and only one
 * cloud type can exist); with 'top' a plume whose top interface
 * lies above σ = `shallowTop` / DEEP_REFERENCE is deep. Any plume that is
 * not deep is handed to the shallow cumulus mass flux above, unchanged. Its CAPE is the positive work of its
 * cloudy layers (`plumeCapeParcel` 'plume'; 'undilute': of the source air
 * lifted without mixing or rain, with the plume's top), its inhibition the
 * negative work below its first cloudy layer. Its mass flux per unit base
 * flux: 1 at the source's top, drawn from the source layers in proportion
 * to their mass, growing by exp((ε − δ) Δz) with δ = max(0, ε −
 * `plumeMassGrowth`) up to the height of neutral buoyancy, interpolated
 * linearly in B between the midpoints of the highest buoyant layer and the
 * one above, and from there falling linearly in height to zero at the top
 * interface. A saturated downdraft starts at the lower interface of the
 * layer of least moist static energy between the plume's top and its
 * cloud-base layer with that layer's air saturated by evaporating rain,
 * its flux `downdraftShare` α of the base flux, growing by entrainment of
 * each layer's air at `downdraftEntrainment` down to the cloud-base layer,
 * then detraining into the layers below in proportion to their mass, moist
 * static energy conserved and saturated at each interface by evaporating
 * rain; α is lowered where the rain made above a level would not cover what
 * the downdraft evaporates down to it. s_l and q_t move in the shallow
 * plume's flux form, the downdraft's flux at each interface α M_d
 * (X_d − X_below), the rain made and evaporated as sources in the layers
 * (−made + evaporated for q_t, L times the opposite for s_l), so column
 * enthalpy and water are exact. With `capeClosure` 'bechtold' (the
 * default) the deep base flux removes the plume's PCAPE less its
 * boundary-layer part over the convective turnover time (Bechtold et al.
 * 2014; IFS Cy43r1 eqs 6.22–6.29): M_b = max(0, PCAPE − PCAPE_bl) / (τ F_P)
 * with PCAPE = Σ (T_v,u − T_v)/T_v Δp over the layers whose work the CAPE
 * counts (condensate loading in both), F_P = Σ (dT_v/dt)/T_v Δp per unit
 * base flux over the layers F counts below, τ = α_x H / w̄ within
 * BECHTOLD.shortest–longest (720–10800 s), H from the plume's base
 * interface to its top, w̄ the depth-weighted mean of the plume's speed
 * over its cloud layers, α_x = 1 + 1.66 dx / 125 km with dx = √(cell area),
 * and PCAPE_bl = τ_bl / T* Σ dT_v/dt|nc Δp over the layers below the
 * plume's base (at most the lowest SUBCLOUD_LAYERS), T* = 1 K, the
 * non-convective tendency the change of each layer's T_v since the end of
 * the previous step's adjustment (`subcloudVirtual`, 0 until one has run:
 * no tendency), τ_bl the base's height over the mass-mean wind speed of
 * those layers (at least 2 m/s) over sea and sea ice, and H / w̄ over
 * `land`; with `pcapeBoundary` 'positive' (the default) a cooling subcloud
 * layer gives PCAPE_bl 0, the boundary-layer production that shallow
 * convection takes up (Bechtold et al. 2014, §2b; max(0, zcape2) in WRF's
 * IFS-derived module_cu_ntiedtke.F), 'signed' keeps its sign. 'threshold' relaxes the CAPE toward `plumeCape` over
 * `plumeRelaxation`: the deep base flux is (CAPE − CAPE0) / (τ F), F the
 * change of the plume's net work over the layers between its source and
 * its top per second and unit base flux, the plume held fixed and the
 * layers moved by the tendencies above (dTv from the change in s_l over
 * c_p and in q_t as vapour). Either is times the deck's opening and an
 * inhibition ramp, 1 up to half `inhibitionThreshold` and 0 from one and a
 * half times it. With `plumeClosure` 'separate' (under convectionType
 * 'top'; 'testParcel' and 'cloudDepth' run it as 'cape') the shallow
 * plume then runs in the same column on its own closure, 'cape' leaves it
 * out where the deep plume runs, 'maximum' gives the deep plume the larger
 * of that base flux and the shallow plume's closure and no shallow plume;
 * where the deep base flux is below CUMULUS_FLOOR the shallow plume runs
 * alone ('separate', 'cape') or nothing does ('maximum'). The base flux is
 * at most `cumulusBoundaryLoss` of the source's mass a step and scaled so
 * that (M + α |M_d|) g Δt / Δp ≤ 1 in every layer. The rain left after the
 * downdraft falls from the layer that made it; below the cloud-base layer
 * it evaporates into each cloud-free subsaturated layer at the fraction
 * 1 − exp(−`plumeRainEvaporation` (1 − q / q_s) Δz) of what falls, at
 * most `rainEvaporation` of what would saturate the layer and never below
 * what the downdraft needs further down; what reaches the ground is
 * convective. The cumulus cover of each layer at or below the shallow
 * top at the highest surface pressure (`cumulusK0`) is the larger of the
 * shallow plume's and the deep plume's M / (ρ w_u), w_u at least w0.
 * With `plumeMomentum` each cell keeps the deep plume's mass fluxes and
 * mixing factors, and `transportMomentum` moves the edges' normal
 * velocity by them (see there).
 *
 * Rain: with `iceFall` (m/s; null: none) the ice share 1 − α of each
 * layer's cloud falls at v = iceFall (ρ q_i / f)^`iceFallExponent`, ρ q_i/f
 * the in-cloud ice content in kg/m³ over the uniform distribution's cover
 * f (Heymsfield and Donner 1990: 3.29 and 0.16; ECHAM6 takes 2.5 at T63
 * and 3.0 at other resolutions), implicitly in flux form
 * from the top down within the step: a layer keeps 1/(1 + v Δt/Δz) of the
 * ice it holds with what fell into it and passes the rest to the layer
 * below, which takes its ice share 1 − α as cloud ice and its liquid share
 * α as precipitation falling on with the autoconversion's (all of it
 * above liquidTemperature, none below iceTemperature), and what
 * leaves the lowest layer is large-scale precipitation, snow or rain by
 * the surface's rule. The column is then adjusted again, so that ice
 * falling into subsaturated layers sublimates, moistening and cooling
 * them. Only the liquid share α of the cloud converts to rain as below.
 * Kessler autoconversion of cloud water above the threshold at
 * autoconversionRate, and of all cloud water over cloudLifetime
 * (`upperCloudLifetime` where the layer's pressure is below `shallowTop`,
 * the anvils' layers; null: cloudLifetime throughout). Stratiform cloud
 * under an inversion lives longer: the lifetime moves to
 * `stratiformLifetime` (3 h; null: no such cloud) by the layer's share
 * s, 0 in the layers at and below the top of a plume that ran in the
 * column this step (cumulusTop), whose cloud the plume detrained, and
 * elsewhere the larger of the cell's sea-ice cover and, with the moist
 * boundary layer (`boundaryTop`, its mixing top, and `boundaryRegime`),
 * 1 below the mixing top of a coupled column, 0 below that of a
 * surface-driven, decoupled or stable one and the radiation's EIS share
 * (`stratiform`) above it. Every layer converts, except
 * in the lowest two layers (`autoconversionFloor` 'lowest') or in the
 * layers wholly below the boundary-layer top ('boundaryLayer'; the
 * lowest two without a boundary layer). The rain falls through the layers below within the
 * step and evaporates into each cloud-free (at most CLEAR_AIR of cloud
 * water) subsaturated one up to `rainEvaporation` of what would saturate
 * it (with `iceSaturation`, of cloudSaturation's mix), latent cooling
 * included. A filler removes negative humidity by
 * borrowing from the layer below.
 *
 * Precipitation accumulates per cell (kg/m²), and so do its two parts:
 * convective, the plumes' rain that reaches the ground; large-scale, the
 * autoconversion rain less what evaporates on the way down. `readRain`
 * turns the parts' sums into their means over the interval they cover,
 * in mm/d. The budget sums are area-weighted masses (kg). `trace`, when
 * its arrays (K·C) are set, receives each layer's temperature change (K)
 * from convection, the evaporation of the deep plume's rain below cloud
 * base included, and from the large-scale condensation, autoconversion
 * and rain evaporation.
 *
 * Defaults: inhibitionThreshold 50 J/kg, shallowTop 700 hPa,
 * autoconversionThreshold 2e-4, autoconversionRate 1e-3 /s,
 * cloudLifetime 1 h, no upperCloudLifetime, stratiformLifetime 3 h, autoconversionFloor 'lowest', rainEvaporation 1,
 * cumulusClosure 0.03 (Grant 2001), cumulusEntrainment 2.5e-3 /m, cumulusDetrainment
 * 3e-3 /m, cumulusSourceDepth 50 hPa, cumulusBoundaryLoss 0.1,
 * cumulusFriction 1, cumulusOvershoot 1, cumulusUpdraft 1 m/s, no
 * cumulusRain, cumulusSource 'mean' (or 'lowest': the plume leaves with
 * the lowest layer's air), plumeClosure 'separate',
 * plumeCapeParcel 'plume', plumeSource 'mean', plumeVelocity 1 m/s,
 * plumeAcceleration 1/3, plumeDrag 1, plumeEntrainmentLaw 'ifs', plumeEntrainment 0.1,
 * plumeEntrainmentFloor 1e-4 /m, plumeMassGrowth 0, plumeConsumption 'all'
 * (or 'buoyant': F counts only the layers whose work the CAPE counts),
 * plumeRainRate 3e-3 /m,
 * plumeRainThreshold 0, plumeRainEvaporation 1e-3 /m, downdraftShare 0.3, convectionType 'testParcel',
 * downdraftEntrainment 1e-4 /m, capeClosure 'bechtold' with pcapeBoundary
 * 'positive' (with 'threshold' plumeCape 120 J/kg and plumeRelaxation 1 h), no
 * plumeMomentum, condensation 'uniform', iceSaturation true, no
 * iceNucleation, iceFall 2.5 m/s, iceFallExponent 0.16.
 */
export const MOIST_DEFAULTS = {
  latentHeat: LATENT_HEAT, inhibitionThreshold: 50, shallowTop: 700e2,
  autoconversionThreshold: 2e-4, autoconversionRate: 1e-3, cloudLifetime: 3600, upperCloudLifetime: null, stratiformLifetime: 3 * 3600, rainEvaporation: 1, autoconversionFloor: 'lowest',
  deckVeto: true, coupledVeto: false, evaporationInCloud: false, virtualBuoyancy: true,
  cumulusClosure: 0.03, cumulusEntrainment: 2.5e-3, cumulusDetrainment: 3e-3, cumulusSourceDepth: 50e2, cumulusBoundaryLoss: 0.1,
  cumulusFriction: 1, cumulusOvershoot: 1, cumulusUpdraft: 1, cumulusRain: null, cumulusSource: 'mean',
  plumeClosure: 'separate', plumeCapeParcel: 'plume', plumeSource: 'mean', plumeSourceDepth: 'surface50', excessVelocity: 'surfaceLayer', plumeVelocity: 1, plumeAcceleration: 1 / 3, plumeDrag: 1, plumeEntrainmentLaw: 'ifs', plumeEntrainment: 0.1, plumeEntrainmentFloor: 1e-4, plumeMassGrowth: 0,
  plumeRainRate: 3e-3, plumeRainThreshold: 0, plumeRainEvaporation: 1e-3, plumePhase: 'mixed', plumeConversion: 'sundqvist', convectionType: 'testParcel', downdraftShare: 0.3, downdraftEntrainment: 1e-4, capeClosure: 'bechtold', pcapeBoundary: 'positive', plumeCape: 120, plumeRelaxation: 3600, plumeMomentum: false, plumeConsumption: 'all',
  condensation: 'uniform', iceSaturation: true, iceNucleation: false, surfaceCriticalHumidity: 0.975, topCriticalHumidity: 0.75, criticalExponent: 2, iceFall: 2.5, iceFallExponent: 0.16,
  liquidTemperature: LIQUID_TEMPERATURE, iceTemperature: ICE_TEMPERATURE,
};
export const RETIRED_OPTIONS = ['convection', 'shallowScheme', 'cumulusWithDeep', 'relaxationTime', 'referenceHumidity', 'parcelDepth', 'entrainmentRate', 'capeThreshold', 'activityMemory', 'detrainment', 'anvilDepth',
  'downdraftEvaporation', 'downdraftSpread', 'shallowHumidity', 'shallowCape', 'shallowInhibition', 'shallowStability', 'shallowReference', 'shallowRain', 'boundaryParcel', 'adjustFrom'];

export function createMoistPhysics(mesh, core, { boundaryDepth = null, boundaryRegime = null, boundaryTop = null, stratiform = null, deckGate = null, surfaceBuoyancy = null, frictionVelocity = null, land = null, surfaceSensible = null, surfaceEvaporation = null, buffers = null, ...options } = {}) {
  for (const retired of RETIRED_OPTIONS) if (retired in options) throw new Error(`${retired} belongs to the retired Betts–Miller convection; the plume is the only scheme`);
  const {
    latentHeat, inhibitionThreshold, shallowTop, autoconversionThreshold, autoconversionRate, cloudLifetime, upperCloudLifetime, stratiformLifetime, rainEvaporation, autoconversionFloor,
    deckVeto, coupledVeto, evaporationInCloud, virtualBuoyancy,
    cumulusClosure, cumulusEntrainment, cumulusDetrainment, cumulusSourceDepth, cumulusBoundaryLoss, cumulusFriction, cumulusOvershoot, cumulusUpdraft, cumulusRain, cumulusSource,
    plumeClosure, plumeCapeParcel, plumeSource, plumeSourceDepth, excessVelocity, plumePhase, plumeConversion, convectionType, plumeVelocity, plumeAcceleration, plumeDrag, plumeEntrainmentLaw, plumeEntrainment, plumeEntrainmentFloor, plumeMassGrowth, plumeRainRate, plumeRainThreshold, plumeRainEvaporation,
    downdraftShare, downdraftEntrainment, capeClosure, pcapeBoundary, plumeCape, plumeRelaxation, plumeMomentum, plumeConsumption,
    condensation, iceSaturation, iceNucleation, surfaceCriticalHumidity, topCriticalHumidity, criticalExponent, iceFall, iceFallExponent, liquidTemperature, iceTemperature,
  } = { ...MOIST_DEFAULTS, ...options };
  if (condensation !== 'uniform' && condensation !== 'saturation') throw new Error(`condensation must be 'uniform' or 'saturation', not ${condensation}`);
  const uniform = condensation === 'uniform', nucleating = iceNucleation && iceSaturation;
  const saturated = { qs: 0, slope: 0, liquid: 1 };
  const saturation = (T, p) => cloudSaturation(T, p, iceSaturation, liquidTemperature, iceTemperature, saturated);
  const halfWidth = (p, ps) => (1 - criticalHumidityAt(p, ps, surfaceCriticalHumidity, topCriticalHumidity, criticalExponent)) * saturated.qs / (1 + latentHeat * saturated.slope / cp);
  if (plumeConsumption !== 'all' && plumeConsumption !== 'buoyant') throw new Error(`plumeConsumption must be 'all' or 'buoyant', not ${plumeConsumption}`);
  const buoyantConsumption = plumeConsumption === 'buoyant';
  if (capeClosure !== 'bechtold' && capeClosure !== 'threshold') throw new Error(`capeClosure must be 'bechtold' or 'threshold', not ${capeClosure}`);
  const bechtold = capeClosure === 'bechtold';
  if (pcapeBoundary !== 'positive' && pcapeBoundary !== 'signed') throw new Error(`pcapeBoundary must be 'positive' or 'signed', not ${pcapeBoundary}`);
  const positiveBoundary = pcapeBoundary === 'positive';
  if (plumeSource !== 'mean' && plumeSource !== 'lowest') throw new Error(`plumeSource must be 'mean' or 'lowest', not ${plumeSource}`);
  const deepLowest = plumeSource === 'lowest';
  if (plumeSourceDepth !== 'surface50' && plumeSourceDepth !== 'boundaryLayer') throw new Error(`plumeSourceDepth must be 'surface50' or 'boundaryLayer', not ${plumeSourceDepth}`);
  const surfaceSource = plumeSourceDepth === 'surface50';
  if (excessVelocity !== 'surfaceLayer' && excessVelocity !== 'convective') throw new Error(`excessVelocity must be 'surfaceLayer' or 'convective', not ${excessVelocity}`);
  const layerExcess = excessVelocity === 'surfaceLayer';
  if (plumeClosure !== 'maximum' && plumeClosure !== 'separate' && plumeClosure !== 'cape') throw new Error(`plumeClosure must be 'maximum', 'separate' or 'cape', not ${plumeClosure}`);
  if (convectionType !== 'testParcel' && convectionType !== 'cloudDepth' && convectionType !== 'top') throw new Error(`convectionType must be 'testParcel', 'cloudDepth' or 'top', not ${convectionType}`);
  const byDepth = convectionType === 'cloudDepth', byParcel = convectionType === 'testParcel';
  if (plumePhase !== 'mixed' && plumePhase !== 'liquid') throw new Error(`plumePhase must be 'mixed' or 'liquid', not ${plumePhase}`);
  const mixedPlume = plumePhase === 'mixed';
  if (plumeConversion !== 'sundqvist' && plumeConversion !== 'zhangMcFarlane') throw new Error(`plumeConversion must be 'sundqvist' or 'zhangMcFarlane', not ${plumeConversion}`);
  const sundqvist = plumeConversion === 'sundqvist';
  if (plumeEntrainmentLaw !== 'ifs' && plumeEntrainmentLaw !== 'gregory') throw new Error(`plumeEntrainmentLaw must be 'ifs' or 'gregory', not ${plumeEntrainmentLaw}`);
  const ifsEntrainment = plumeEntrainmentLaw === 'ifs';
  const separate = plumeClosure === 'separate' && convectionType === 'top', relaxedOnly = plumeClosure !== 'maximum';
  if (plumeCapeParcel !== 'plume' && plumeCapeParcel !== 'undilute') throw new Error(`plumeCapeParcel must be 'plume' or 'undilute', not ${plumeCapeParcel}`);
  const undilute = plumeCapeParcel === 'undilute';
  if (cumulusSource !== 'mean' && cumulusSource !== 'lowest') throw new Error(`cumulusSource must be 'mean' or 'lowest', not ${cumulusSource}`);
  if (autoconversionFloor !== 'lowest' && autoconversionFloor !== 'boundaryLayer' && autoconversionFloor !== 'none') throw new Error(`autoconversionFloor must be 'lowest' or 'boundaryLayer', not ${autoconversionFloor}`);
  const { K, C, levels, dSigma, sigmaMid, cp, R, g, kappa, exnerLayer, exnerLower, geopotential } = core.diagnostics;
  const thetaV = core.arrays.thetaV;
  const KL = Math.min(SUBCLOUD_LAYERS, K), subcloudBuffer = buffers && buffers.subcloudVirtual ? buffers.subcloudVirtual : new SharedArrayBuffer(8 * KL * C);
  const subcloudVirtual = new Float64Array(subcloudBuffer), windVector = new Float64Array(3 * C);
  const resolutionScale = Float64Array.from(mesh.areaCell, (a) => 1 + BECHTOLD.resolution * Math.sqrt(a) / BECHTOLD.reference);
  const upperInterface = (i, k) => (geopotential[k * C + i] + cp * thetaV[k * C + i] * (exnerLayer[k * C + i] - exnerLower[(k - 1) * C + i])) / g;
  const shared = (name, n) => (buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * n));
  const precipBuffer = shared('precipitation', C), rainBuffer = shared('rain', C), convectiveBuffer = shared('convectivePrecipitation', C);
  const largeScaleBuffer = shared('largeScalePrecipitation', C);
  const cumulusCoverBuffer = shared('cumulusCover', K * C), cumulusWaterBuffer = shared('cumulusWater', K * C), baseFluxBuffer = shared('cumulusBaseFlux', C), cumulusTopBuffer = shared('cumulusTop', C);
  const momentumLayers = plumeMomentum ? K : 0;
  const momentumBuffers = { up: shared('momentumUp', (momentumLayers + 1) * C), upKeep: shared('momentumUpKeep', momentumLayers * C), down: shared('momentumDown', (momentumLayers + 1) * C), downKeep: shared('momentumDownKeep', momentumLayers * C), source: shared('momentumSource', C) };
  const momentumUp = new Float64Array(momentumBuffers.up), momentumUpKeep = new Float64Array(momentumBuffers.upKeep), momentumDown = new Float64Array(momentumBuffers.down), momentumDownKeep = new Float64Array(momentumBuffers.downKeep), momentumSource = new Float64Array(momentumBuffers.source);
  const edgeUp = new Float64Array(K + 1), edgeDown = new Float64Array(K + 1), edgeKeep = new Float64Array(K), edgeDownKeep = new Float64Array(K), edgeFlux = new Float64Array(K + 1), edgeBefore = new Float64Array(K);
  const cumulusCover = new Float64Array(cumulusCoverBuffer), cumulusWater = new Float64Array(cumulusWaterBuffer), cumulusBaseFlux = new Float64Array(baseFluxBuffer), cumulusTop = new Float64Array(cumulusTopBuffer);
  const precipitation = new Float64Array(precipBuffer), rain = new Float64Array(rainBuffer);
  const convectivePrecipitation = new Float64Array(convectiveBuffer), largeScalePrecipitation = new Float64Array(largeScaleBuffer);
  const convectiveRain = new Float64Array(C), largeScaleRain = new Float64Array(C);
  const T = new Float64Array(K), p = new Float64Array(K), dp = new Float64Array(K), z = new Float64Array(K);
  const downdraftCooling = new Float64Array(K);
  const plumeS = new Float64Array(K), plumeQ = new Float64Array(K), plumeLiquid = new Float64Array(K), plumeGrowth = new Float64Array(K), plumeRain = new Float64Array(K);
  const envS = new Float64Array(K), envQ = new Float64Array(K), cumulusFlux = new Float64Array(K + 1), fluxS = new Float64Array(K + 1), fluxQ = new Float64Array(K + 1);
  const envHumidity = new Float64Array(K), plumeDetrained = new Float64Array(K), envVirtual = new Float64Array(K), plumeCounted = new Uint8Array(K), plumeSpeed = new Float64Array(K + 1), plumeWork = new Float64Array(K), plumeBuoyancy = new Float64Array(K), plumeEntrained = new Float64Array(K), plumeDepth = new Float64Array(K);
  const draftFlux = new Float64Array(K + 1), draftS = new Float64Array(K + 1), draftQ = new Float64Array(K + 1), draftEvaporation = new Float64Array(K);
  const deepCover = new Float64Array(K), deepWater = new Float64Array(K), tendencyS = new Float64Array(K), tendencyQ = new Float64Array(K), convectiveFall = new Float64Array(K), convectiveReserve = new Float64Array(K);
  const plumeCarried = new Float64Array(K), plumeInterfaceT = new Float64Array(K), plumeFrozen = new Float64Array(K), convectiveMelted = new Float64Array(K), convectiveFrozen = new Float64Array(K), convectiveSnow = new Float64Array(C);
  const plume = { T: 0, liquid: 0, ice: 0 }, draft = { q: 0, s: 0 }, plumeSaturated = { qs: 0, slope: 0, liquid: 1 };
  const deep = { deep: false, shallowRain: 0, shallowSnow: 0, top: -1, base: K, cape: 0, consumption: 0, inhibition: 0, start: -1, baseFlux: 0, downdraft: 0, pcape: 0, pcapeBoundary: 0, consumptionP: 0, tau: 0, speed: 0, depth: 0, boundaryWind: 0, boundaryTime: 0, excessT: 0, excessQ: 0, sourceS: 0, sourceQ: 0, sourceMass: 0, parcelDeep: false, parcelBase: 0, parcelTop: 0, conversionMargin: Infinity, typeMargin: Infinity };
  let cumulusK0 = K;
  while (cumulusK0 > 0 && 0.5 * (levels[cumulusK0 - 1] + levels[cumulusK0]) * MAXIMUM_SURFACE_PRESSURE > shallowTop) cumulusK0--;
  const cumulus = { top: -1, source: K - 1, inhibition: 0, lclPressure: 0, velocity: 0, baseFlux: 0, snow: 0 };
  const falling = { evaporated: 0, convective: 0, frozen: 0, ice: 0, moved: false };
  const budget = { condensation: 0, convection: 0, lost: 0 };
  const trace = { convection: null, largeScale: null };
  const marked = new Float64Array(K);
  const ramp = (x) => Math.min(1, Math.max(0, x));
  function mark(i, theta) { for (let k = 0; k < K; k++) marked[k] = theta[k * C + i]; }
  function charge(into, i, theta) {
    for (let k = 0; k < K; k++) { const idx = k * C + i; if (into) into[idx] += (theta[idx] - marked[k]) * exnerLayer[idx]; marked[k] = theta[idx]; }
  }

  /*
   * The temperature at which saturated air at `pressure` holds `energy`
   * as cp T + L q_s(T, p), by four Newton steps from `guess`.
   */
  function saturatedTemperature(energy, pressure, guess) {
    let t = guess;
    for (let n = 0; n < 4; n++) {
      const qs = saturationHumidity(t, pressure);
      t -= (cp * t + latentHeat * qs - energy) / (cp + latentHeat * latentHeat * qs / (R_VAPOR * t * t));
    }
    return t;
  }

  /*
   * Saturation adjustment: supersaturated vapour condenses into cloud
   * water and cloud water evaporates into subsaturated air, each with
   * one implicit step, so afterwards a layer is either saturated or
   * cloud-free. Returns the condensate formed (kg/m², negative when
   * cloud evaporated); nothing rains here.
   */
  function condenseColumn(i, pi, theta, q, qc) {
    let formed = 0;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      const ex = exnerLayer[idx];
      const temperature = theta[idx] * ex;
      const pressure = pi[i] * sigmaMid[k];
      let change;
      if (uniform && !(boundaryTop !== null && geopotential[idx] / g < boundaryTop[i])) {
        const liquidT = temperature - latentHeat * qc[idx] / cp;
        saturation(liquidT, pressure);
        const a = 1 / (1 + latentHeat * saturated.slope / cp), b = halfWidth(pressure, pi[i]), total = q[idx] + qc[idx], Q = a * (total - saturated.qs);
        if (nucleating && !(qc[idx] > CLEAR_AIR) && liquidT < iceTemperature) {
          const reference = Math.min(saturationHumidity(liquidT, pressure), (2.583 - liquidT / 207.8) * saturated.qs), width = b / (a * saturated.qs) * reference;
          const lowest = Math.max(reference, total - width);
          change = total + width > reference ? Math.max(0, a * (total + width - lowest) / (2 * width) * (0.5 * (lowest + total + width) - saturated.qs)) - qc[idx] : -qc[idx];
        } else change = (b > 0 ? uniformCondensate(Q, b) : Math.max(0, Q)) - qc[idx];
      } else {
        const qs = iceSaturation ? saturation(temperature, pressure).qs : saturationHumidity(temperature, pressure);
        const slope = iceSaturation ? saturated.slope : qs * latentHeat / (R_VAPOR * temperature * temperature);
        change = (q[idx] - qs) / (1 + latentHeat * slope / cp);
      }
      if (change < 0) change = Math.max(change, -qc[idx]);
      if (change === 0) continue;
      q[idx] -= change;
      qc[idx] += change;
      theta[idx] += latentHeat * change / (cp * ex);
      formed += pi[i] * dSigma[k] / g * change;
    }
    return formed;
  }

  /*
   * Autoconversion and the rain's fall (see the header). `stream` holds
   * each layer's share of the deep plume's rain (plumeColumn's
   * convectiveFall); falling.evaporated receives what of it evaporated
   * below cloud base and falling.convective what reached the ground.
   * Returns the autoconversion rain that reaches the ground (kg/m²).
   */
  function longCloudShare(i, k, pi, iced) {
    if (cumulusBaseFlux[i] > 0 && pi[i] * sigmaMid[k] >= cumulusTop[i]) return 0;
    let share = 0;
    if (boundaryTop !== null) {
      if (!(geopotential[k * C + i] / g < boundaryTop[i])) share = stratiform !== null ? stratiform[i] : 0;
      else if (boundaryRegime !== null && boundaryRegime[i] === COUPLED_REGIME) share = 1;
    }
    return Math.max(share, iced);
  }

  function autoconvertColumn(i, pi, theta, q, qc, dt, stream = null, iced = 0, frozenStream = null) {
    let rain = 0, convective = 0, streamed = 0, descending = 0, frozen = 0;
    falling.moved = false;
    const floor = autoconversionFloor === 'boundaryLayer' && boundaryDepth ? boundaryDepth[i] : null;
    if (trace.convection) downdraftCooling.fill(0);
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      if (rain > 0 && (evaporationInCloud || !(qc[idx] > CLEAR_AIR)) && rainEvaporation > 0) {
        const ex = exnerLayer[idx], mass = pi[i] * dSigma[k] / g;
        const temperature = theta[idx] * ex;
        const qs = iceSaturation ? saturation(temperature, pi[i] * sigmaMid[k]).qs : saturationHumidity(temperature, pi[i] * sigmaMid[k]);
        const slope = iceSaturation ? saturated.slope : qs * latentHeat / (R_VAPOR * temperature * temperature);
        const deficit = Math.max(0, (qs - q[idx]) / (1 + latentHeat * slope / cp)) * mass;
        const evaporated = Math.min(rain, rainEvaporation * deficit);
        if (evaporated > 0) {
          rain = Math.max(0, rain - evaporated);
          q[idx] += evaporated / mass;
          theta[idx] -= latentHeat * evaporated / (mass * cp * ex);
        }
      }
      if (stream) {
        convective = Math.max(0, convective + stream[k]);
        if (frozenStream) frozen = Math.min(convective, Math.max(0, frozen + frozenStream[k]));
        const spare = convective - convectiveReserve[k];
        if (spare > 0 && k > deep.base && plumeRainEvaporation > 0 && rainEvaporation > 0 && (evaporationInCloud || !(qc[idx] > CLEAR_AIR))) {
          const ex = exnerLayer[idx], mass = pi[i] * dSigma[k] / g;
          const temperature = theta[idx] * ex;
          const qs = saturationHumidity(temperature, pi[i] * sigmaMid[k]);
          const slope = qs * latentHeat / (R_VAPOR * temperature * temperature);
          const airborne = -convective * Math.expm1(-plumeRainEvaporation * Math.max(0, 1 - q[idx] / qs) * R * temperature * dSigma[k] / (sigmaMid[k] * g));
          const evaporated = Math.min(spare, airborne, rainEvaporation * Math.max(0, (qs - q[idx]) / (1 + latentHeat * slope / cp)) * mass);
          if (evaporated > 0) {
            const sublimated = frozen > 0 ? evaporated * frozen / convective : 0;
            convective -= evaporated;
            streamed += evaporated;
            q[idx] += evaporated / mass;
            if (frozenStream) {
              frozen = Math.max(0, frozen - sublimated);
              theta[idx] -= (latentHeat * evaporated + FUSION_HEAT * sublimated) / (mass * cp * ex);
              if (trace.convection) downdraftCooling[k] += (latentHeat * evaporated + FUSION_HEAT * sublimated) / (mass * cp);
            } else {
              theta[idx] -= latentHeat * evaporated / (mass * cp * ex);
              if (trace.convection) downdraftCooling[k] += latentHeat * evaporated / (mass * cp);
            }
          }
        }
      }
      let liquid = qc[idx];
      if (iceFall !== null) {
        const mass = pi[i] * dSigma[k] / g;
        if (descending > 0) {
          const melted = descending * liquidFraction(theta[idx] * exnerLayer[idx], liquidTemperature, iceTemperature);
          qc[idx] += (descending - melted) / mass;
          rain += melted;
          descending = 0;
        }
        liquid = qc[idx];
        if (qc[idx] > 0) {
          const temperature = theta[idx] * exnerLayer[idx], pressure = pi[i] * sigmaMid[k];
          const share = liquidFraction(temperature, liquidTemperature, iceTemperature), ice = (1 - share) * qc[idx];
          liquid = share * qc[idx];
          if (ice > 0) {
            saturation(temperature, pressure);
            const cover = Math.max(CLEAR_AIR, uniformCover(qc[idx], halfWidth(pressure, pi[i])));
            const speed = iceFall * Math.pow(pressure / (R * temperature) * ice / cover, iceFallExponent);
            const courant = speed * dt * sigmaMid[k] * g / (R * temperature * dSigma[k]), leaving = ice * courant / (1 + courant);
            qc[idx] -= leaving;
            descending = leaving * mass;
            falling.moved = true;
          }
        }
      }
      if (!(qc[idx] > 0)) continue;
      if (autoconversionFloor !== 'none' && (floor === null ? k >= K - 2 : k > 0 && upperInterface(i, k) < floor)) continue;
      const excess = Math.max(0, liquid - autoconversionThreshold);
      let lifetime = upperCloudLifetime !== null && pi[i] * sigmaMid[k] < shallowTop ? upperCloudLifetime : cloudLifetime;
      if (stratiformLifetime !== null) lifetime += longCloudShare(i, k, pi, iced) * (stratiformLifetime - lifetime);
      const converted = Math.min(liquid, excess * (1 - Math.exp(-autoconversionRate * dt)) + liquid * (1 - Math.exp(-dt / lifetime)));
      qc[idx] -= converted;
      rain += pi[i] * dSigma[k] / g * converted;
    }
    falling.evaporated = streamed;
    falling.convective = convective;
    falling.frozen = Math.min(frozen, convective);
    falling.ice = descending;
    return rain + descending;
  }

  function clearCumulus(i) {
    cumulusBaseFlux[i] = 0; cumulusTop[i] = 0;
    for (let k = 0; k < K; k++) { cumulusCover[k * C + i] = 0; cumulusWater[k * C + i] = 0; }
  }

  /*
   * The temperature and condensate of plume air of liquid-water static
   * energy `energy` and total water `water` at `height` and `pressure`.
   */
  function plumeState(energy, water, height, pressure, guess) {
    if (mixedPlume) { frozenPlumeState(energy, water, height, pressure, guess); return; }
    const dry = (energy - g * height) / cp;
    if (!(water > saturationHumidity(dry, pressure))) { plume.T = dry; plume.liquid = 0; return; }
    const t = saturatedTemperature(energy - g * height + latentHeat * water, pressure, Math.max(dry, guess));
    plume.T = t;
    plume.liquid = Math.max(0, water - saturationHumidity(t, pressure));
  }

  /*
   * With `plumePhase` 'mixed': `energy` is s_li = c_p T + g z − L q_l − (L +
   * L_f) q_i, the condensate l = q_t − q_s,mix(T) of cloudSaturation's mix,
   * q_i = (1 − α(T)) l; T solves c_p T + L q_s − L_f (1 − α)(q_t − q_s) =
   * s_li − g z + L q_t by four Newton steps; plume.ice is q_i.
   */
  function frozenPlumeState(energy, water, height, pressure, guess) {
    const dry = (energy - g * height) / cp;
    plume.ice = 0;
    if (!(water > cloudSaturation(dry, pressure, true, liquidTemperature, iceTemperature, plumeSaturated).qs)) { plume.T = dry; plume.liquid = 0; return; }
    const target = energy - g * height + latentHeat * water, span = 1 / (liquidTemperature - iceTemperature);
    let t = Math.max(dry, guess);
    for (let n = 0; n < 4; n++) {
      const { qs, slope, liquid } = cloudSaturation(t, pressure, true, liquidTemperature, iceTemperature, plumeSaturated), held = water - qs, frozen = 1 - liquid;
      t -= (cp * t + latentHeat * qs - FUSION_HEAT * frozen * held - target) / (cp + (latentHeat + FUSION_HEAT * frozen) * slope + (liquid > 0 && liquid < 1 ? FUSION_HEAT * span * held : 0));
    }
    cloudSaturation(t, pressure, true, liquidTemperature, iceTemperature, plumeSaturated);
    plume.T = t;
    plume.liquid = Math.max(0, water - plumeSaturated.qs);
    plume.ice = (1 - plumeSaturated.liquid) * plume.liquid;
  }

  function fillEnvironment(i, pi, theta, q, qc) {
    for (let k = 0; k < K; k++) {
      const idx = k * C + i, cloud = qc ? Math.max(0, qc[idx]) : 0;
      T[k] = theta[idx] * exnerLayer[idx];
      p[k] = pi[i] * sigmaMid[k];
      dp[k] = pi[i] * dSigma[k];
      z[k] = geopotential[idx] / g;
      envS[k] = cp * T[k] + g * z[k] - latentHeat * cloud;
      envQ[k] = Math.max(0, q[idx]) + cloud;
    }
  }

  function sourceLayers(i, pi, lowest, surface = false) {
    const bottom = K - 1, depth = boundaryDepth[i];
    let mass = 0, energy = 0, water = 0, source = bottom;
    for (let k = bottom; k >= 0; k--) {
      if (k < bottom && (((surface || !(upperInterface(i, k + 1) < depth)) && !(p[k] >= pi[i] - cumulusSourceDepth)) || !(p[k] > shallowTop))) break;
      mass += dp[k]; energy += dp[k] * envS[k]; water += dp[k] * envQ[k];
      source = k;
    }
    return { source, mass, sourceS: lowest ? envS[bottom] : energy / mass, sourceQ: lowest ? envQ[bottom] : water / mass };
  }

  /*
   * The shallow cumulus mass flux of column i over dt (see the header).
   * Returns its rain (kg/m²).
   */
  function cumulusColumn(i, pi, theta, q, qc, dt) {
    const bottom = K - 1;
    cumulus.top = -1; cumulus.inhibition = 0; cumulus.baseFlux = 0; cumulus.velocity = 0; cumulus.snow = 0;
    clearCumulus(i);
    const open = (deckVeto && deckGate !== null ? ramp((DECK_CLOSED - deckGate[i]) / (DECK_CLOSED - DECK_OPEN)) : 1) * (coupledVeto && boundaryRegime !== null && boundaryRegime[i] === COUPLED_REGIME ? 0 : 1);
    const buoyancy = surfaceBuoyancy ? surfaceBuoyancy[i] : 0;
    if (!(open > 0) || !(buoyancy > 0) || !boundaryDepth) return 0;
    fillEnvironment(i, pi, theta, q, qc);
    const { source, mass, sourceS, sourceQ } = sourceLayers(i, pi, cumulusSource === 'lowest');
    const depth = boundaryDepth[i], lowest = cumulusSource === 'lowest';
    const lcl = liftingCondensationLevel((sourceS - g * z[bottom]) / cp, sourceQ, p[bottom], kappa);
    if (!lcl || !(lcl.pressure > shallowTop)) return 0;
    cumulus.source = source; cumulus.lclPressure = lcl.pressure;
    const virtual = virtualBuoyancy ? VIRTUAL_FACTOR : 0, loading = virtualBuoyancy ? 1 : 0;
    let s = sourceS, w = sourceQ, inhibition = 0, cloudy = false, top = -1, guess = 0;
    plumeS[source] = s; plumeQ[source] = w;
    for (let k = source - 1; k >= 0; k--) {
      if (!(p[k] > shallowTop)) { if (cloudy) top = k + 1; break; }
      const below = upperInterface(i, k + 1), above = upperInterface(i, k);
      const mixes = pi[i] * levels[k + 1] <= lcl.pressure, epsilon = mixes ? cumulusEntrainment : 0;
      const half = Math.exp(-epsilon * (z[k] - below));
      const midS = envS[k] + (s - envS[k]) * half, midQ = envQ[k] + (w - envQ[k]) * half;
      plumeState(midS, midQ, z[k], p[k], guess);
      guess = plume.T;
      plumeLiquid[k] = plume.liquid;
      const air = Math.max(0, q[k * C + i]), cloud = qc ? Math.max(0, qc[k * C + i]) : 0;
      const work = R * (plume.T * (1 + virtual * (midQ - plume.liquid) - loading * plume.liquid) - T[k] * (1 + virtual * air - loading * cloud)) * dp[k] / p[k];
      if (plume.liquid > 0) cloudy = true;
      if (!cloudy) { if (work < 0) inhibition -= work; } else if (!(work > 0)) { top = k; break; }
      const full = Math.exp(-epsilon * (above - below));
      s = envS[k] + (s - envS[k]) * full; w = envQ[k] + (w - envQ[k]) * full;
      plumeRain[k] = 0; plumeFrozen[k] = 0;
      if (cumulusRain !== null) {
        plumeState(s, w, above, pi[i] * levels[k], guess);
        const excess = plume.liquid - cumulusRain;
        if (excess > 0) { w -= excess; s += latentHeat * excess; plumeRain[k] = excess; if (mixedPlume) { plumeFrozen[k] = excess * plume.ice / plume.liquid; s += FUSION_HEAT * plumeFrozen[k]; } }
      }
      plumeS[k] = s; plumeQ[k] = w;
      plumeGrowth[k] = Math.exp((epsilon - (mixes ? cumulusDetrainment : 0)) * (above - below));
    }
    cumulus.inhibition = inhibition;
    if (top < 0) return 0;
    const convective = Math.cbrt(buoyancy * Math.max(0, depth - z[bottom]));
    const velocity = Math.max(convective, cumulusFriction * (frictionVelocity ? frictionVelocity[i] : 0));
    cumulus.velocity = velocity;
    if (!(velocity > 0)) return 0;
    let base = open * cumulusClosure * lcl.pressure / (R * lcl.temperature) * velocity * Math.exp(-inhibition / (velocity * velocity));
    base = Math.min(base, cumulusBoundaryLoss * mass / (g * dt));
    if (!(base > CUMULUS_FLOOR)) return 0;
    cumulusFlux.fill(0);
    let sourceBelow = 0;
    for (let j = bottom; j > source; j--) { sourceBelow += dp[j]; cumulusFlux[j] = base * sourceBelow / mass; }
    cumulusFlux[source] = base;
    for (let k = source - 1; k > top; k--) cumulusFlux[k] = cumulusFlux[k + 1] * plumeGrowth[k];
    cumulusFlux[top + 1] *= cumulusOvershoot;
    let scale = 1;
    for (let k = top; k <= bottom; k++) {
      const courant = Math.max(cumulusFlux[k], cumulusFlux[k + 1]) * g * dt / dp[k];
      if (courant * scale > 1) scale = 1 / courant;
    }
    for (let j = top + 1; j <= bottom; j++) cumulusFlux[j] *= scale;
    let belowMass = 0, belowS = 0, belowQ = 0;
    fluxS[K] = 0; fluxQ[K] = 0; fluxS[top] = 0; fluxQ[top] = 0;
    for (let j = bottom; j > top; j--) {
      let upS = plumeS[j], upQ = plumeQ[j];
      if (j > source) {
        belowMass += dp[j]; belowS += dp[j] * envS[j]; belowQ += dp[j] * envQ[j];
        upS = lowest ? sourceS : belowS / belowMass; upQ = lowest ? sourceQ : belowQ / belowMass;
      }
      fluxS[j] = cumulusFlux[j] * (upS - envS[j - 1]);
      fluxQ[j] = cumulusFlux[j] * (upQ - envQ[j - 1]);
    }
    let rain = 0;
    for (let k = top; k <= bottom; k++) {
      const idx = k * C + i, per = g * dt / dp[k];
      let dS = (fluxS[k + 1] - fluxS[k]) * per, dQ = (fluxQ[k + 1] - fluxQ[k]) * per;
      if (cumulusRain !== null && k > top && k < source && plumeRain[k] > 0) {
        const fallen = cumulusFlux[k] * plumeRain[k] * dt;
        rain += fallen; dQ -= fallen * g / dp[k]; dS += latentHeat * fallen * g / dp[k];
        if (mixedPlume && plumeFrozen[k] > 0) { const frozen = cumulusFlux[k] * plumeFrozen[k] * dt; cumulus.snow += frozen; dS += FUSION_HEAT * frozen * g / dp[k]; }
      }
      theta[idx] += dS / (cp * exnerLayer[idx]);
      q[idx] += dQ;
    }
    for (let k = top; k < source; k++) {
      if (!(plumeLiquid[k] > 0)) continue;
      const idx = k * C + i;
      cumulusCover[idx] = Math.min(1, 0.5 * (cumulusFlux[k] + cumulusFlux[k + 1]) * R * T[k] / (p[k] * cumulusUpdraft));
      cumulusWater[idx] = plumeLiquid[k];
    }
    cumulus.top = top; cumulus.baseFlux = base * scale;
    cumulusBaseFlux[i] = base * scale; cumulusTop[i] = pi[i] * levels[top];
    return rain;
  }

  /*
   * The temperature, vapour and liquid-water static energy of saturated
   * downdraft air of moist static energy `energy` at `height` and
   * `pressure` into `draft`.
   */
  function saturatedDraft(energy, height, pressure, guess) {
    const t = saturatedTemperature(energy - g * height, pressure, guess);
    draft.q = saturationHumidity(t, pressure);
    draft.s = energy - latentHeat * draft.q;
  }

  /*
   * The IFS's first-guess deep updraught (see the header) from the deep
   * source's s_l and q_t, on the environment fillEnvironment holds: true
   * where its cloud is deeper than DEEP_CLOUD_DEPTH, with deep.parcelBase and
   * deep.parcelTop its cloud's base and top pressures (the top at least the
   * interface it reached once the cloud is deep).
   */
  function parcelDeep(i, pi, q, qc, source, s, w) {
    const bottom = K - 1, virtual = virtualBuoyancy ? VIRTUAL_FACTOR : 0, loading = virtualBuoyancy ? 1 : 0;
    const surfaceSaturation = saturation(T[bottom], p[bottom]).qs;
    let w2 = plumeVelocity * plumeVelocity, guess = 0, base = 0;
    deep.typeMargin = Infinity;
    for (let k = source - 1; k > 0; k--) {
      const lower = upperInterface(i, k + 1), upper = upperInterface(i, k), depth = upper - lower;
      const epsilon = TEST_PARCEL.entrainment * IFS_ENTRAINMENT.entrainment * (saturation(T[k], p[k]).qs / surfaceSaturation) ** IFS_ENTRAINMENT.scaleExponent;
      const half = Math.exp(-epsilon * (z[k] - lower));
      const midS = envS[k] + (s - envS[k]) * half, midQ = envQ[k] + (w - envQ[k]) * half;
      plumeState(midS, midQ, z[k], p[k], guess);
      guess = plume.T;
      if (!(base > 0) && plume.liquid > 0) base = pi[i] * levels[k + 1];
      const idx = k * C + i, environment = T[k] * (1 + virtual * Math.max(0, q[idx]) - loading * (qc ? Math.max(0, qc[idx]) : 0));
      const buoyancy = g * (plume.T * (1 + virtual * (midQ - plume.liquid) - loading * plume.liquid) - environment) / environment;
      const mixing = IFS_ENTRAINMENT.drag * epsilon, x = 2 * mixing * depth;
      const next = w2 * Math.exp(-x) + 2 * plumeAcceleration * buoyancy * depth * (x > 0 ? -Math.expm1(-x) / x : 1);
      if (!(next > 0)) {
        if (!(base > 0)) return false;
        const lift = -plumeAcceleration * buoyancy, y = mixing * w2 / lift, reach = mixing > 0 ? (y > 1e-3 ? Math.log1p(y) : y * (1 - 0.5 * y)) / (2 * mixing) : w2 / (2 * lift);
        deep.parcelBase = base;
        deep.parcelTop = pi[i] * levels[k + 1] * Math.pow(levels[k] / levels[k + 1], Math.min(1, reach / depth));
        deep.parcelDeep = base - deep.parcelTop > DEEP_CLOUD_DEPTH;
        deep.typeMargin = Math.min(deep.typeMargin, Math.abs(base - deep.parcelTop - DEEP_CLOUD_DEPTH) / DEEP_CLOUD_DEPTH);
        return deep.parcelDeep;
      }
      w2 = next;
      deep.typeMargin = Math.min(deep.typeMargin, w2 / (plumeVelocity * plumeVelocity));
      if (base > 0 && base - pi[i] * levels[k] > DEEP_CLOUD_DEPTH) { deep.parcelBase = base; deep.parcelTop = pi[i] * levels[k]; deep.parcelDeep = true; return true; }
      const full = Math.exp(-epsilon * depth);
      s = envS[k] + (s - envS[k]) * full; w = envQ[k] + (w - envQ[k]) * full;
      plumeState(s, w, upper, pi[i] * levels[k], guess);
      w -= TEST_PARCEL.removal * plume.liquid; s += latentHeat * TEST_PARCEL.removal * plume.liquid;
      if (mixedPlume) s += FUSION_HEAT * TEST_PARCEL.removal * plume.ice;
    }
    if (!(base > 0)) return false;
    deep.parcelBase = base; deep.parcelTop = pi[i] * levels[1]; deep.parcelDeep = base - deep.parcelTop > DEEP_CLOUD_DEPTH;
    return deep.parcelDeep;
  }

  /*
   * The convective mass flux of column i over dt (see the header).
   * Returns the rain it leaves falling (kg/m²);
   * convectiveFall holds each layer's share of it and convectiveReserve
   * what must still fall past each layer for the downdraft below it.
   */
  function plumeColumn(i, pi, theta, q, qc, dt, u = null) {
    const bottom = K - 1;
    convectiveFall.fill(0); convectiveReserve.fill(0);
    if (mixedPlume) { convectiveFrozen.fill(0); convectiveMelted.fill(0); }
    deep.deep = false; deep.top = -1; deep.cape = 0; deep.consumption = 0; deep.start = -1; deep.downdraft = 0; deep.baseFlux = 0; deep.inhibition = 0; deep.shallowRain = 0; deep.shallowSnow = 0;
    deep.pcape = 0; deep.pcapeBoundary = 0; deep.consumptionP = 0; deep.tau = 0; deep.speed = 0; deep.depth = 0; deep.boundaryWind = 0; deep.boundaryTime = 0; deep.excessT = 0; deep.excessQ = 0;
    deep.parcelDeep = false; deep.parcelBase = 0; deep.parcelTop = 0; deep.conversionMargin = Infinity; deep.typeMargin = Infinity; cumulus.snow = 0;
    if (momentumLayers) {
      momentumSource[i] = K;
      for (let k = 0; k <= K; k++) { momentumUp[k * C + i] = 0; momentumDown[k * C + i] = 0; }
      for (let k = 0; k < K; k++) { momentumUpKeep[k * C + i] = 1; momentumDownKeep[k * C + i] = 1; }
    }
    const open = (deckVeto && deckGate !== null ? ramp((DECK_CLOSED - deckGate[i]) / (DECK_CLOSED - DECK_OPEN)) : 1) * (coupledVeto && boundaryRegime !== null && boundaryRegime[i] === COUPLED_REGIME ? 0 : 1);
    if (!(open > 0) || !boundaryDepth) return cumulusColumn(i, pi, theta, q, qc, dt);
    fillEnvironment(i, pi, theta, q, qc);
    const { source, mass, sourceS: meanS, sourceQ: meanQ } = sourceLayers(i, pi, false, surfaceSource);
    let sourceS = deepLowest ? envS[bottom] : meanS, sourceQ = deepLowest ? envQ[bottom] : meanQ;
    if (surfaceSource && surfaceSensible && surfaceEvaporation) {
      const buoyancy = surfaceBuoyancy ? surfaceBuoyancy[i] : 0, density = p[bottom] / (R * T[bottom]), b = bottom * C + i;
      const velocity = layerExcess ? surfaceLayerVelocity(surfaceSensible[i], surfaceEvaporation[i], density, T[bottom], cp * thetaV[b] * (exnerLower[b] - exnerLayer[b]) / g, cp, g)
        : Math.max(buoyancy > 0 ? Math.cbrt(buoyancy * Math.max(0, boundaryDepth[i] - z[bottom])) : 0, frictionVelocity ? frictionVelocity[i] : 0);
      if (velocity > 0) {
        deep.excessT = Math.min(SOURCE_EXCESS.temperature, SOURCE_EXCESS.coefficient * surfaceSensible[i] / (density * cp * velocity));
        deep.excessQ = Math.min(SOURCE_EXCESS.humidity, SOURCE_EXCESS.coefficient * surfaceEvaporation[i] / (density * velocity));
        if (layerExcess) { deep.excessT = Math.max(0, deep.excessT); deep.excessQ = Math.max(0, deep.excessQ); }
        sourceS += cp * deep.excessT; sourceQ += deep.excessQ;
      }
    }
    deep.sourceS = sourceS; deep.sourceQ = sourceQ; deep.sourceMass = mass;
    const lcl = liftingCondensationLevel((sourceS - g * z[bottom]) / cp, sourceQ, p[bottom], kappa);
    if (!lcl || !(lcl.pressure > pi[i] * levels[1])) return cumulusColumn(i, pi, theta, q, qc, dt);
    if (byParcel && !parcelDeep(i, pi, q, qc, source, sourceS, sourceQ)) return cumulusColumn(i, pi, theta, q, qc, dt);
    const virtual = virtualBuoyancy ? VIRTUAL_FACTOR : 0, loading = virtualBuoyancy ? 1 : 0;
    let s = sourceS, w = sourceQ, w2 = 0, below = 0, inhibition = 0, cloudy = false, started = false, top = -1, guess = 0, cape = 0, base = -1, pcape = 0, baseSaturation = 0;
    plumeS[source] = s; plumeQ[source] = w;
    plumeSpeed.fill(0); plumeCounted.fill(0);
    for (let k = source - 1; k >= 0; k--) {
      const lower = upperInterface(i, k + 1), upper = k > 0 ? upperInterface(i, k) : Infinity, depth = upper - lower;
      const mixes = pi[i] * levels[k + 1] <= lcl.pressure;
      if (mixes && !started) { started = true; base = k + 1; w2 = plumeVelocity * plumeVelocity; plumeSpeed[k + 1] = w2; if (ifsEntrainment) baseSaturation = saturation(T[k], p[k]).qs; }
      let epsilon = 0, mixing = 0;
      if (ifsEntrainment) {
        const humidity = Math.min(1, Math.max(0, q[k * C + i]) / saturation(T[k], p[k]).qs);
        envHumidity[k] = humidity;
        plumeDetrained[k] = mixes ? IFS_ENTRAINMENT.detrainment * (IFS_ENTRAINMENT.detrainmentHumidity - humidity) : 0;
        epsilon = mixes && below > 0 ? IFS_ENTRAINMENT.entrainment * (IFS_ENTRAINMENT.humidity - humidity) * (saturated.qs / baseSaturation) ** IFS_ENTRAINMENT.scaleExponent : 0;
        mixing = IFS_ENTRAINMENT.drag * (epsilon > 0 ? epsilon : plumeDetrained[k]);
      } else {
        epsilon = mixes ? Math.max(plumeEntrainmentFloor, plumeEntrainment * Math.max(0, below) / w2) : 0;
        mixing = plumeDrag * epsilon;
      }
      const half = Math.exp(-epsilon * (z[k] - lower));
      const midS = envS[k] + (s - envS[k]) * half, midQ = envQ[k] + (w - envQ[k]) * half;
      plumeState(midS, midQ, z[k], p[k], guess);
      guess = plume.T;
      plumeLiquid[k] = plume.liquid;
      const idx = k * C + i, air = Math.max(0, q[idx]), cloud = qc ? Math.max(0, qc[idx]) : 0;
      const environment = T[k] * (1 + virtual * air - loading * cloud), rising = plume.T * (1 + virtual * (midQ - plume.liquid) - loading * plume.liquid);
      const work = R * (rising - environment) * dp[k] / p[k], buoyancy = g * (rising - environment) / environment;
      if (plume.liquid > 0) cloudy = true;
      if (!cloudy && work < 0) inhibition -= work;
      plumeWork[k] = work; plumeBuoyancy[k] = buoyancy; plumeEntrained[k] = epsilon; plumeDepth[k] = depth; envVirtual[k] = environment;
      if (undilute) {
        plumeState(sourceS, sourceQ, z[k], p[k], plume.T);
        plumeWork[k] = R * (plume.T * (1 + virtual * (sourceQ - plume.liquid)) - environment) * dp[k] / p[k];
      }
      if (mixes) {
        if (!(depth < Infinity)) { top = k; break; }
        const x = 2 * mixing * depth, decay = Math.exp(-x);
        w2 = w2 * decay + 2 * plumeAcceleration * buoyancy * depth * (x > 0 ? -Math.expm1(-x) / x : 1);
        if (!(w2 > 0)) { top = k; break; }
      }
      plumeSpeed[k] = mixes ? w2 : 0;
      if (cloudy && plumeWork[k] > 0) { cape += plumeWork[k]; plumeCounted[k] = 1; if (bechtold) pcape += plumeWork[k] * p[k] / (R * environment); }
      below = buoyancy;
      const full = Math.exp(-epsilon * depth);
      s = envS[k] + (s - envS[k]) * full; w = envQ[k] + (w - envQ[k]) * full;
      plumeRain[k] = 0; plumeFrozen[k] = 0; plumeCarried[k] = 0;
      if (mixes && sundqvist) {
        plumeState(s, w, upper, pi[i] * levels[k], guess);
        const threshold = land && land[i] ? IFS_PRECIPITATION.landThreshold : IFS_PRECIPITATION.seaThreshold;
        plumeInterfaceT[k] = plume.T;
        if (plume.liquid > 0) deep.conversionMargin = Math.min(deep.conversionMargin, Math.abs(plume.liquid - threshold) / threshold);
        if (plume.liquid > threshold) {
          const alpha = liquidFraction(plume.T, liquidTemperature, iceTemperature), P = IFS_PRECIPITATION;
          const bergeron = plume.T < P.bergeron ? 1 + 0.5 * Math.sqrt(Math.min(P.bergeron - plume.T, P.bergeron - P.ice)) : 1;
          const speed = Math.min(P.speed, Math.max(plumeVelocity, Math.sqrt(w2))), critical = P.critical / bergeron;
          const rate = P.conversion * (P.liquidFactor * alpha + 1 - alpha) * bergeron / (P.velocityScale * speed) * -Math.expm1(-((plume.liquid / critical) ** 2));
          const fallen = -plume.liquid * Math.expm1(-rate * depth); w -= fallen; s += latentHeat * fallen; plumeRain[k] = fallen;
          if (mixedPlume && plume.ice > 0) { plumeFrozen[k] = fallen * plume.ice / plume.liquid; s += FUSION_HEAT * plumeFrozen[k]; }
        }
        plumeCarried[k] = plume.liquid - plumeRain[k];
      } else if (mixes) {
        plumeState(s, w, upper, pi[i] * levels[k], guess);
        const excess = plume.liquid - plumeRainThreshold;
        if (excess > 0) {
          const fallen = -excess * Math.expm1(-plumeRainRate * depth); w -= fallen; s += latentHeat * fallen; plumeRain[k] = fallen;
          if (mixedPlume && plume.ice > 0) { plumeFrozen[k] = fallen * plume.ice / plume.liquid; s += FUSION_HEAT * plumeFrozen[k]; }
        }
        plumeCarried[k] = plume.liquid - plumeRain[k];
      }
      plumeS[k] = s; plumeQ[k] = w;
    }
    if (top === 0) top = 1;
    if (!cloudy || top < 1 || !(byParcel || (byDepth ? pi[i] * (levels[base] - levels[top]) > DEEP_CLOUD_DEPTH : levels[top] * DEEP_REFERENCE < shallowTop))) return cumulusColumn(i, pi, theta, q, qc, dt);
    clearCumulus(i);
    deep.deep = true; deep.top = top; deep.cape = cape; deep.inhibition = inhibition; deep.base = base;
    let neutral = -1;
    for (let k = top + 1; k < source; k++) if (plumeBuoyancy[k] > 0) { neutral = k; break; }
    const topHeight = upperInterface(i, top);
    let neutralHeight = upperInterface(i, top + 1);
    if (neutral > 0 && !(plumeBuoyancy[neutral - 1] > 0)) neutralHeight = Math.min(neutralHeight, z[neutral] + (z[neutral - 1] - z[neutral]) * plumeBuoyancy[neutral] / (plumeBuoyancy[neutral] - plumeBuoyancy[neutral - 1]));
    cumulusFlux.fill(0);
    let sourceBelow = 0;
    for (let j = bottom; j > source; j--) { sourceBelow += dp[j]; cumulusFlux[j] = sourceBelow / mass; }
    cumulusFlux[source] = 1;
    if (ifsEntrainment) {
      for (let k = source - 1; k > top; k--) {
        if (k >= base) cumulusFlux[k] = cumulusFlux[k + 1];
        else if (plumeBuoyancy[k] > 0) cumulusFlux[k] = cumulusFlux[k + 1] * Math.exp((plumeEntrained[k] - plumeDetrained[k]) * plumeDepth[k]);
        else cumulusFlux[k] = cumulusFlux[k + 1] * Math.exp(-plumeDetrained[k] * plumeDepth[k]) * Math.min(1, (IFS_ENTRAINMENT.detrainmentHumidity - envHumidity[k]) * Math.sqrt(plumeSpeed[k] / plumeSpeed[k + 1]));
      }
    } else {
      let k = source - 1;
      for (; k > top && !(upperInterface(i, k) > neutralHeight); k--) {
        const epsilon = plumeEntrained[k];
        cumulusFlux[k] = cumulusFlux[k + 1] * Math.exp((epsilon - Math.max(0, epsilon - plumeMassGrowth)) * plumeDepth[k]);
      }
      for (const anchor = cumulusFlux[k + 1]; k > top; k--) cumulusFlux[k] = anchor * (topHeight - upperInterface(i, k)) / (topHeight - neutralHeight);
    }
    let rainAbove = 0;
    for (let k = top + 1; k < source; k++) rainAbove += cumulusFlux[k] * plumeRain[k];
    let flyingIce = 0;
    if (mixedPlume) {
      for (let k = top; k <= bottom; k++) {
        convectiveMelted[k] = 0;
        if (flyingIce > 0 && T[k] >= liquidTemperature) { convectiveMelted[k] = flyingIce; flyingIce = 0; }
        convectiveFrozen[k] = -convectiveMelted[k];
        if (k > top && k < source && plumeFrozen[k] > 0) { flyingIce += cumulusFlux[k] * plumeFrozen[k]; convectiveFrozen[k] += cumulusFlux[k] * plumeFrozen[k]; }
      }
    }
    draftFlux.fill(0); draftEvaporation.fill(0);
    let start = -1, share = 0;
    if (downdraftShare > 0 && rainAbove > 0) {
      for (let k = top + 1; k < base; k++) if (start < 0 || envS[k] + latentHeat * envQ[k] < envS[start] + latentHeat * envQ[start]) start = k;
    }
    if (start >= 0) {
      let energy = envS[start] + latentHeat * envQ[start], flux = 1, subcloud = 0;
      for (let k = base + 1; k <= bottom; k++) subcloud += dp[k];
      saturatedDraft(energy, upperInterface(i, start + 1), pi[i] * levels[start + 1], T[start]);
      draftFlux[start + 1] = -1; draftS[start + 1] = draft.s; draftQ[start + 1] = draft.q;
      draftEvaporation[start] = draft.q - envQ[start];
      let water = draft.q, left = subcloud, atBase = 1;
      for (let k = start + 1; k < bottom; k++) {
        let mixedQ = water;
        if (k <= base) {
          const keep = Math.exp(-downdraftEntrainment * (upperInterface(i, k) - upperInterface(i, k + 1))), air = envS[k] + latentHeat * envQ[k];
          flux /= keep;
          energy = air + (energy - air) * keep;
          mixedQ = envQ[k] + (water - envQ[k]) * keep;
          atBase = flux;
        } else {
          left -= dp[k];
          flux = atBase * left / subcloud;
        }
        saturatedDraft(energy, upperInterface(i, k + 1), pi[i] * levels[k + 1], T[k]);
        draftFlux[k + 1] = -flux; draftS[k + 1] = draft.s; draftQ[k + 1] = draft.q;
        draftEvaporation[k] = (draft.q - mixedQ) * flux;
        water = draft.q;
      }
      share = downdraftShare;
      let produced = 0, taken = 0;
      for (let k = 0; k < bottom; k++) {
        if (k > top && k < source) produced += mixedPlume ? cumulusFlux[k] * (plumeRain[k] - plumeFrozen[k]) : cumulusFlux[k] * plumeRain[k];
        if (mixedPlume && k >= top) produced += convectiveMelted[k];
        taken += draftEvaporation[k];
        if (taken > 0 && share * taken > produced) share = produced / taken;
      }
      if (!(share > 0)) { share = 0; draftFlux.fill(0); draftEvaporation.fill(0); }
      deep.start = start;
    }
    let belowMass = 0, belowS = 0, belowQ = 0;
    fluxS.fill(0); fluxQ.fill(0);
    for (let j = bottom; j > top; j--) {
      let upS = plumeS[j], upQ = plumeQ[j];
      if (j > source) {
        belowMass += dp[j]; belowS += dp[j] * envS[j]; belowQ += dp[j] * envQ[j];
        upS = deepLowest ? sourceS : belowS / belowMass; upQ = deepLowest ? sourceQ : belowQ / belowMass;
      }
      fluxS[j] = cumulusFlux[j] * (upS - envS[j - 1]) + share * draftFlux[j] * (draftS[j] - envS[j]);
      fluxQ[j] = cumulusFlux[j] * (upQ - envQ[j - 1]) + share * draftFlux[j] * (draftQ[j] - envQ[j]);
    }
    let consumption = 0, consumptionP = 0;
    for (let k = top; k <= bottom; k++) {
      const per = g / dp[k], made = k > top && k < source ? cumulusFlux[k] * plumeRain[k] : 0, evaporated = share * draftEvaporation[k];
      tendencyS[k] = (mixedPlume ? fluxS[k + 1] - fluxS[k] + latentHeat * (made - evaporated) + FUSION_HEAT * convectiveFrozen[k] : fluxS[k + 1] - fluxS[k] + latentHeat * (made - evaporated)) * per;
      tendencyQ[k] = (fluxQ[k + 1] - fluxQ[k] - made + evaporated) * per;
      convectiveFall[k] = made - evaporated;
      if (k < source && k > top && (!buoyantConsumption || plumeCounted[k])) {
        const idx = k * C + i, air = Math.max(0, q[idx]), cloud = qc ? Math.max(0, qc[idx]) : 0;
        const warming = tendencyS[k] / cp * (1 + virtual * air - loading * cloud) + virtual * T[k] * tendencyQ[k];
        consumption += R * warming * dp[k] / p[k];
        if (bechtold) consumptionP += warming * dp[k] / envVirtual[k];
      }
    }
    deep.consumption = consumption;
    let relaxed = 0;
    if (bechtold) {
      const baseHeight = upperInterface(i, base), cloudDepth = topHeight - baseHeight;
      let weighted = 0, thickness = 0;
      for (let k = top; k < base; k++) { weighted += plumeDepth[k] * Math.sqrt(Math.max(0, 0.5 * (plumeSpeed[k] + plumeSpeed[k + 1]))); thickness += plumeDepth[k]; }
      const speed = weighted / thickness, turnover = cloudDepth / speed;
      const tau = Math.min(BECHTOLD.longest, Math.max(BECHTOLD.shortest, resolutionScale[i] * turnover));
      const b = bottom * C + i, ground = (geopotential[b] - cp * thetaV[b] * (exnerLower[b] - exnerLayer[b])) / g;
      let forcing = 0, mass = 0, wind = 0;
      for (let k = Math.max(base, K - KL); k <= bottom; k++) {
        const idx = k * C + i, saved = subcloudVirtual[(k - K + KL) * C + i];
        if (saved > 0) forcing += (T[k] * (1 + virtual * Math.max(0, q[idx]) - loading * (qc ? Math.max(0, qc[idx]) : 0)) - saved) / dt * dp[k];
        if (u) { cellVector(mesh, u.subarray(k * mesh.nEdges, (k + 1) * mesh.nEdges), windVector, i, i + 1); wind += dp[k] * Math.hypot(windVector[3 * i], windVector[3 * i + 1], windVector[3 * i + 2]); }
        mass += dp[k];
      }
      const boundaryWind = Math.max(BECHTOLD.boundaryWind, wind / mass), boundaryTime = land && land[i] ? turnover : (baseHeight - ground) / boundaryWind;
      const pcapeBoundary = boundaryTime / BECHTOLD.temperatureScale * (positiveBoundary ? Math.max(0, forcing) : forcing);
      relaxed = consumptionP > 0 ? Math.max(0, pcape - pcapeBoundary) / (tau * consumptionP) : 0;
      deep.pcape = pcape; deep.pcapeBoundary = pcapeBoundary; deep.consumptionP = consumptionP; deep.tau = tau; deep.speed = speed; deep.depth = cloudDepth; deep.boundaryWind = boundaryWind; deep.boundaryTime = boundaryTime;
    } else relaxed = consumption > 0 && cape > plumeCape ? (cape - plumeCape) / (plumeRelaxation * consumption) : 0;
    const gate = ramp(0.5 + (inhibitionThreshold - inhibition) / Math.max(1, inhibitionThreshold));
    const buoyancy = surfaceBuoyancy ? surfaceBuoyancy[i] : 0;
    let shallowBase = 0;
    if (buoyancy > 0 && !relaxedOnly) {
      const convective = Math.cbrt(buoyancy * Math.max(0, boundaryDepth[i] - z[bottom]));
      const velocity = Math.max(convective, cumulusFriction * (frictionVelocity ? frictionVelocity[i] : 0));
      if (velocity > 0) shallowBase = cumulusClosure * lcl.pressure / (R * lcl.temperature) * velocity * Math.exp(-inhibition / (velocity * velocity));
    }
    let baseFlux = Math.min(open * Math.max(shallowBase, gate * relaxed), cumulusBoundaryLoss * mass / (g * dt));
    for (let k = top; k <= bottom; k++) {
      const courant = baseFlux * (Math.max(cumulusFlux[k], cumulusFlux[k + 1]) + share * Math.max(-draftFlux[k], -draftFlux[k + 1])) * g * dt / dp[k];
      if (courant > 1) baseFlux /= courant;
    }
    if (!(baseFlux > CUMULUS_FLOOR)) { deep.deep = false; return relaxedOnly ? cumulusColumn(i, pi, theta, q, qc, dt) : 0; }
    let fallen = 0;
    for (let k = top; k <= bottom; k++) {
      const idx = k * C + i;
      theta[idx] += baseFlux * dt * tendencyS[k] / (cp * exnerLayer[idx]);
      q[idx] += baseFlux * dt * tendencyQ[k];
      convectiveFall[k] *= baseFlux * dt;
      if (mixedPlume) convectiveFrozen[k] *= baseFlux * dt;
      fallen += convectiveFall[k];
    }
    for (let k = bottom - 1; k >= 0; k--) convectiveReserve[k] = Math.max(0, convectiveReserve[k + 1] - convectiveFall[k + 1]);
    for (let k = Math.max(top, cumulusK0); k < source; k++) {
      deepCover[k] = 0; deepWater[k] = plumeLiquid[k];
      if (!(plumeLiquid[k] > 0)) continue;
      const speed = Math.max(plumeVelocity, Math.sqrt(Math.max(0, 0.5 * (plumeSpeed[k] + plumeSpeed[k + 1]))));
      deepCover[k] = Math.min(1, 0.5 * (cumulusFlux[k] + cumulusFlux[k + 1]) * baseFlux * R * T[k] / (p[k] * speed));
    }
    if (momentumLayers) {
      for (let k = 0; k <= K; k++) { momentumUp[k * C + i] = baseFlux * cumulusFlux[k]; momentumDown[k * C + i] = baseFlux * share * draftFlux[k]; }
      for (let k = 0; k < K; k++) {
        momentumUpKeep[k * C + i] = k > top && k < source ? Math.exp(-plumeEntrained[k] * plumeDepth[k]) : 1;
        momentumDownKeep[k * C + i] = start < 0 ? 1 : k === start ? 0 : k > start && k <= base ? Math.exp(-downdraftEntrainment * (upperInterface(i, k) - upperInterface(i, k + 1))) : 1;
      }
      momentumSource[i] = source;
    }
    const shallowRain = separate ? cumulusColumn(i, pi, theta, q, qc, dt) : 0;
    for (let k = Math.max(top, cumulusK0); k < source; k++) {
      const idx = k * C + i;
      if (deepCover[k] > cumulusCover[idx]) { cumulusCover[idx] = deepCover[k]; cumulusWater[idx] = deepWater[k]; }
    }
    deep.baseFlux = baseFlux; deep.downdraft = share * baseFlux; deep.shallowRain = shallowRain; deep.shallowSnow = separate ? cumulus.snow : 0;
    cumulus.top = top; cumulus.baseFlux = baseFlux + cumulusBaseFlux[i];
    cumulusBaseFlux[i] += baseFlux; cumulusTop[i] = pi[i] * levels[top];
    return fallen;
  }

  /*
   * The deep plume's transport of the normal velocity of edges eFrom to
   * eTo over dt with `plumeMomentum`: on each edge the mean of its two
   * cells' updraft and downdraft mass fluxes, mixing factors and the
   * shallower source, the updraft leaving with the mass-weighted velocity
   * of the layers below each source interface and mixing toward each
   * layer's velocity above, the downdraft starting with its first layer's
   * velocity; the same flux form as s_l and q_t with the edge's layer
   * masses, so each edge's column momentum is exact. The kinetic energy
   * removed goes to `dissipation` as the boundary layer's does.
   */
  function transportMomentum(pi, u, eFrom, eTo, dt, dissipation = null) {
    if (!momentumLayers) return;
    const { cellsOnEdge } = mesh, E = mesh.nEdges, bottom = K - 1;
    for (let e = eFrom; e < eTo; e++) {
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      let moving = false;
      for (let k = 1; k < K; k++) {
        edgeUp[k] = 0.5 * (momentumUp[k * C + a] + momentumUp[k * C + b]);
        edgeDown[k] = 0.5 * (momentumDown[k * C + a] + momentumDown[k * C + b]);
        if (edgeUp[k] > 0 || edgeDown[k] < 0) moving = true;
      }
      if (!moving) continue;
      for (let k = 0; k < K; k++) {
        edgeKeep[k] = 0.5 * (momentumUpKeep[k * C + a] + momentumUpKeep[k * C + b]);
        edgeDownKeep[k] = 0.5 * (momentumDownKeep[k * C + a] + momentumDownKeep[k * C + b]);
        edgeBefore[k] = u[k * E + e];
      }
      const source = Math.min(momentumSource[a], momentumSource[b]), columnMass = 0.5 * (pi[a] + pi[b]);
      edgeFlux.fill(0);
      let rising = 0, below = 0, weight = 0;
      for (let j = bottom; j >= 1; j--) {
        if (j >= source) { below += dSigma[j] * edgeBefore[j]; weight += dSigma[j]; rising = below / weight; }
        else rising = edgeBefore[j] + (rising - edgeBefore[j]) * edgeKeep[j];
        edgeFlux[j] = edgeUp[j] * (rising - edgeBefore[j - 1]);
      }
      let sinking = 0;
      for (let j = 1; j < K; j++) {
        sinking = edgeBefore[j - 1] + (sinking - edgeBefore[j - 1]) * edgeDownKeep[j - 1];
        edgeFlux[j] += edgeDown[j] * (sinking - edgeBefore[j]);
      }
      let loss = 0, total = 0;
      for (let k = 0; k < K; k++) {
        const mass = columnMass * dSigma[k] / g, now = edgeBefore[k] + (edgeFlux[k + 1] - edgeFlux[k]) * dt / mass;
        u[k * E + e] = now;
        loss += mass * (edgeBefore[k] * edgeBefore[k] - now * now);
        edgeUp[k] = mass * (now - edgeBefore[k]) * (now - edgeBefore[k]);
        total += edgeUp[k];
      }
      if (!dissipation || !(loss > 0) || !(total > 0)) continue;
      for (let k = 0; k < K; k++) dissipation[k * E + e] += loss * edgeUp[k] / (total * columnMass * dSigma[k] / g);
    }
  }

  function fillColumn(i, pi, q) {
    if (!q) return;
    for (let k = 0; k < K - 1; k++) {
      const idx = k * C + i;
      if (q[idx] < 0) {
        q[idx + C] += q[idx] * dSigma[k] / dSigma[k + 1];
        q[idx] = 0;
      }
    }
    const bottom = (K - 1) * C + i;
    if (q[bottom] < 0) {
      budget.lost -= mesh.areaCell[i] * pi[i] * dSigma[K - 1] / g * q[bottom];
      q[bottom] = 0;
    }
  }

  function saveSubcloud(i, theta, q, qc) {
    const virtual = virtualBuoyancy ? VIRTUAL_FACTOR : 0, loading = virtualBuoyancy ? 1 : 0;
    for (let k = K - KL; k < K; k++) {
      const idx = k * C + i;
      subcloudVirtual[(k - K + KL) * C + i] = theta[idx] * exnerLayer[idx] * (1 + virtual * Math.max(0, q[idx]) - loading * (qc ? Math.max(0, qc[idx]) : 0));
    }
  }

  let iceConcentration = null;
  function useSeaIce(concentration) { iceConcentration = concentration; }

  function adjust(state, iFrom, iTo, dt) {
    const [pi, theta, u, , q, qc, ice = null] = state;
    for (let i = iFrom; i < iTo; i++) {
      core.diagnoseColumn(i, pi, theta, q, qc);
      const traced = trace.convection || trace.largeScale;
      if (traced) mark(i, theta);
      condenseColumn(i, pi, theta, q, qc);
      if (traced) charge(trace.largeScale, i, theta);
      const produced = plumeColumn(i, pi, theta, q, qc, dt, u);
      if (traced) charge(trace.convection, i, theta);
      if (cumulusBaseFlux[i] > 0) condenseColumn(i, pi, theta, q, qc);
      const rained = autoconvertColumn(i, pi, theta, q, qc, dt, deep.deep ? convectiveFall : null, ice && ice[i] > 0 ? (iceConcentration !== null && iceConcentration[i] > 0 ? iceConcentration[i] : 1) : 0, deep.deep && mixedPlume ? convectiveFrozen : null);
      const convected = deep.deep ? falling.convective + deep.shallowRain : produced;
      if (mixedPlume) convectiveSnow[i] = deep.deep ? falling.frozen + deep.shallowSnow : cumulus.snow;
      if (falling.moved) condenseColumn(i, pi, theta, q, qc);
      if (traced) {
        charge(trace.largeScale, i, theta);
        if (trace.convection) {
          for (let k = 0; k < K; k++) {
            const idx = k * C + i;
            trace.convection[idx] -= downdraftCooling[k];
            if (trace.largeScale) trace.largeScale[idx] += downdraftCooling[k];
          }
        }
      }
      fillColumn(i, pi, q);
      fillColumn(i, pi, qc);
      if (bechtold) saveSubcloud(i, theta, q, qc);
      precipitation[i] += rained + convected;
      convectivePrecipitation[i] += convected;
      largeScalePrecipitation[i] += rained;
      rain[i] = rained + convected;
      budget.condensation += mesh.areaCell[i] * rained;
      budget.convection += mesh.areaCell[i] * convected;
    }
  }

  function readRain(interval) {
    const scale = 86400 / interval;
    for (let i = 0; i < C; i++) {
      convectiveRain[i] = scale * convectivePrecipitation[i];
      largeScaleRain[i] = scale * largeScalePrecipitation[i];
    }
  }

  function columnWater(pi, q, i) {
    let water = 0;
    for (let k = 0; k < K; k++) water += pi[i] * dSigma[k] / g * q[k * C + i];
    return water;
  }

  const settings = { uniform, iceSaturation, surfaceCriticalHumidity, topCriticalHumidity, criticalExponent, liquidTemperature, iceTemperature };
  return {
    condensation: settings, adjust, useSeaIce, condenseColumn, autoconvertColumn, cumulusColumn, plumeColumn, transportMomentum, fillColumn, columnWater, readRain,
    plumeAir(energy, water, height, pressure, guess) { plumeState(energy, water, height, pressure, guess); return { ...plume }; },
    precipitation, rain, convectivePrecipitation, largeScalePrecipitation, convectiveRain, largeScaleRain, convectiveSnow, budget, latentHeat, trace, falling,
    cumulus, cumulusCover, cumulusWater, cumulusBaseFlux, cumulusTop, cumulusFlux, deep, subcloudVirtual, subcloudLayers: KL, saveSubcloud, deepSigma: shallowTop / DEEP_REFERENCE, cumulusK0, convectiveFall, convectiveFrozen, convectiveMelted, plumeFrozen, plumeCarried, plumeInterfaceT, draftFlux, plumeSpeed, plumeRain, plumeBuoyancy, plumeEntrained, plumeDetrained, envHumidity, plumeCounted, resolutionScale,
    momentum: { up: momentumUp, upKeep: momentumUpKeep, down: momentumDown, downKeep: momentumDownKeep, source: momentumSource },
    shared: { momentumUp: momentumBuffers.up, momentumUpKeep: momentumBuffers.upKeep, momentumDown: momentumBuffers.down, momentumDownKeep: momentumBuffers.downKeep, momentumSource: momentumBuffers.source, precipitation: precipBuffer, rain: rainBuffer, convectivePrecipitation: convectiveBuffer, largeScalePrecipitation: largeScaleBuffer, cumulusCover: cumulusCoverBuffer, cumulusWater: cumulusWaterBuffer, cumulusBaseFlux: baseFluxBuffer, cumulusTop: cumulusTopBuffer, subcloudVirtual: subcloudBuffer },
  };
}
