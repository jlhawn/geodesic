export const FREEZING_POINT = 271.35;
export const MELTING_POINT = 273.15;
const DEVICE_MELTING_POINT = Math.fround(MELTING_POINT);
export const MINIMUM_CONCENTRATION = 0.01;
export const MINIMUM_VOLUME = 1e-4;
export const NEGLIGIBLE_CONCENTRATION = 1e-6;
export const NEGLIGIBLE_VOLUME = 1e-6;
export const COVER_GRAIN = 2 ** -24;

/*
 * Open water forms ice only from MINIMUM_CONCENTRATION of its cell and
 * MINIMUM_VOLUME of volume up. Ice that exists melts away, its latent heat
 * drawn from the water, below the same area or volume; with `keepSmallIce`
 * (a drift carries ice into open cells a few thousandths of a cell at a time)
 * only below NEGLIGIBLE_CONCENTRATION or NEGLIGIBLE_VOLUME, what single-
 * precision arithmetic cannot tell from zero: floats near one are 2^-23
 * apart, and a cover or a volume is a sum of several terms of order one.
 * CICE removes ice by the same rule, at `puny` = 1e-11 of the cell in double
 * precision (Hunke et al., the CICE Consortium's CICE documentation, index
 * of primary variables; Icepack's icepack_itd.F90 zap_small_areas and
 * icepack_parameters.F90, itd_area_min = 1e-11), and sends the enthalpy,
 * fresh water and salt of what it removes to the ocean as meltAway does.
 */
export const seaIceRemoval = (keepSmallIce) => (keepSmallIce ? { concentration: NEGLIGIBLE_CONCENTRATION, volume: NEGLIGIBLE_VOLUME } : { concentration: MINIMUM_CONCENTRATION, volume: MINIMUM_VOLUME });

/*
 * How the ice and the water under it trade heat, `iceExchange`. 'instant'
 * gives the ice all the heat the ocean's top layer holds above freezing,
 * through `oceanFlux`, and the ocean holds that layer at the freezing point
 * under ice; the heat the leads gain melts ice. 'turbulent' lets the water
 * under the ice keep the heat it holds above freezing, `waterHeat` in J/m² of
 * cell, and the ice take it at the rates real ice does. The turbulent
 * boundary layer under the ice carries heat to its base at
 *   F_bot = ρ_w c_p c_h u* (T_w − T_f)
 * per unit ice area (McPhee 1992; Maykut and McPhee 1995; the form of CICE and
 * Icepack, Hunke et al., the Icepack documentation, "Basal heat flux"), with
 * c_h = 0.006 and u* = max(u*_min, √C_w |u_ice − u_ocean|), the friction
 * velocity of the ice–ocean stress. The ocean hands the ice ρ_w c_p c_h u*
 * as `conductance` (W m⁻² K⁻¹) each time it steps, and u*_min is 5·10⁻⁴ m/s,
 * the value of CICE's namelist and the one Icepack's documentation
 * recommends (Icepack's source default is 5·10⁻³ and the documentation
 * allows 0). A floor of 2·10⁻² m/s stands in for stirring the model does not
 * resolve (tides, internal waves, eddies): √C_w |u_ice − u_ocean| passes it
 * only above about 0.27 m/s, so the drift then barely sets the conductance,
 * about 490 W m⁻² K⁻¹. Floes also melt at their sides, at the rate
 *   w_lat = m1 (T_w − T_f)^m2,  m1 = 1.6·10⁻⁶ m s⁻¹ K^−m2,  m2 = 1.36
 * (Maykut and Perovich 1987), which takes the fraction w_lat π/(α D) of their
 * area per second for floes of diameter D = 300 m and shape α = 0.66
 * (Steele 1992; Rothrock 1984), and the same fraction of their volume, so
 * that the thickness stays. T_w − T_f is waterHeat over the top layer's heat
 * capacity. Bottom and lateral melt are scaled together where they would take
 * more than waterHeat in the step, as CICE scales them against the heat the
 * ocean can give. The constants are ICE_EXCHANGE; `iceExchange` is the mode
 * name or an object of them with a `mode`.
 */
export const ICE_EXCHANGES = ['instant', 'turbulent'];
export const DEFAULT_ICE_EXCHANGE = 'instant';
export const ICE_EXCHANGE = { heatTransfer: 0.006, frictionMinimum: 5e-4, lateralMelt: 1.6e-6, lateralExponent: 1.36, floeDiameter: 300, floeShape: 0.66 };
export function iceExchangeOptions(option = DEFAULT_ICE_EXCHANGE) {
  const given = typeof option === 'string' ? { mode: option } : option;
  if (given === null || typeof given !== 'object') throw new Error(`iceExchange is ${ICE_EXCHANGES.join(' or ')} or an object of ${Object.keys(ICE_EXCHANGE).join(', ')} with a mode, not ${option}`);
  const { mode, turbulent, lateralRate, ...rest } = given;
  if (mode === undefined) throw new Error(`iceExchange given as an object names its mode, ${ICE_EXCHANGES.join(' or ')}: ${JSON.stringify(given)}`);
  if (!ICE_EXCHANGES.includes(mode)) throw new Error(`iceExchange is one of ${ICE_EXCHANGES.join(', ')}, not ${mode}`);
  for (const key of Object.keys(rest)) if (!(key in ICE_EXCHANGE)) throw new Error(`iceExchange has no option ${key}; its options are ${Object.keys(ICE_EXCHANGE).join(', ')}`);
  const set = { ...ICE_EXCHANGE, ...rest };
  for (const key of ['heatTransfer', 'lateralExponent', 'floeDiameter', 'floeShape']) if (!(set[key] > 0)) throw new Error(`iceExchange.${key} is positive, not ${set[key]}`);
  for (const key of ['frictionMinimum', 'lateralMelt']) if (!(set[key] >= 0)) throw new Error(`iceExchange.${key} is not negative, not ${set[key]}`);
  const resolved = { mode, turbulent: mode === 'turbulent', ...set, lateralRate: set.lateralMelt * Math.PI / (set.floeShape * set.floeDiameter) };
  if (turbulent !== undefined && turbulent !== resolved.turbulent) throw new Error(`iceExchange.turbulent follows its mode ${mode}, not ${turbulent}`);
  if (lateralRate !== undefined && lateralRate !== resolved.lateralRate) throw new Error(`iceExchange.lateralRate follows lateralMelt, floeShape and floeDiameter (${resolved.lateralRate}), not ${lateralRate}`);
  return resolved;
}

// A header that records no exchange is an instant save, whatever DEFAULT_ICE_EXCHANGE is; scripts/savedLand.mjs resumedIceExchange refuses one that holds ocean.iceWater.
export const savedIceExchange = (saved) => saved?.given?.ocean?.iceExchange ?? saved?.options?.ocean?.iceExchange ?? 'instant';

/*
 * The slab ocean and its sea ice, a zero-layer thermodynamic model: the
 * E1, E3 and E4 simplifications of Fichefet and Morales Maqueda (1997, JGR
 * 102, 12609; their Table 2: no heat stored in the ice, no sunlight into
 * it, one thickness) to Semtner's (1976) three-layer model, by default with
 * no brine reservoir; iceLight lifts E3, brineReservoir the latent half of
 * E1. Open water is a mixed layer of heat capacity slabHeatCapacity; when
 * it cools to the freezing point the deficit freezes into ice of that
 * latent heat. Ice has a skin of small heat capacity whose temperature
 * answers the surface flux and the conduction through the ice from its
 * base at the freezing point; the conducted heat freezes water onto the
 * base, and a skin that would
 * pass the melting point melts the ice from the top instead. Ice that
 * melts away returns its leftover energy to the mixed layer. The surface
 * energy — mixed-layer heat over the freezing point, skin heat, minus
 * the ice's latent heat — changes by exactly the surface flux plus
 * `oceanFlux`, the heat above freezing that the ocean's mixed layer
 * hands to the base of the ice. surfaceT is the skin temperature of
 * the ice, or of the water where there is none. Open water reflects
 * the direct beam with the zenith-angle albedo of Briegleb et al.
 * (1986), 0.02 under a high sun and 0.3 near the horizon, and diffuse
 * light (`albedo` without a zenith cosine) with diffuseWaterAlbedo,
 * unless oceanAlbedo fixes both. With `heatCapacity` given, that per-cell
 * array (a dynamic ocean's upper layer) replaces slabHeatCapacity in the
 * cell update; with `transmission` and `penetrated` too, open water keeps
 * the share transmission[i] of the absorbed sunlight `solar` out of its
 * temperature and adds it to penetrated[i] (J/m²) for the ocean to take
 * below its mixed layer.
 *
 * Ice covers the fraction `concentration` of its cell (A), with the
 * thickness `ice` over that part, so the volume per cell area is A·h;
 * the open part, the leads, is water held at the freezing point. The
 * atmosphere sees the area-weighted albedo and, for its fluxes, the
 * area-weighted surface temperature. The cell's net flux is split
 * between the parts by the sunlight the leads absorb beyond the ice
 * (`contrast` in update) and by the heat the leads at the freezing
 * point lose beyond the ice, leadExchange per kelvin of the ice skin's
 * difference from freezing, so that the ice's flux per unit ice area
 * and the water's per unit water area average back to it. The ice part
 * evolves as above; the heat the leads gain or lose, with the ocean's
 * flux under them, melts ice or freezes new ice. Melting takes area as
 * Hibler (1979) does, half the relative loss of volume from the area,
 * times areaLossFactor;
 * ice frozen in the leads closes them as (1 − A) of new ice leadClosing
 * thick, so the open area closes as its square, or with leadClosure
 * 'hibler' as all of it, Hibler's (1979) h₀ law with h₀ = leadClosing
 * (Hibler and FMM97's Table 1 take 0.5 m). Open
 * water that cools below freezing forms ice leadClosing thick over the
 * area its volume covers. Ice below the floors of seaIceRemoval melts
 * away, and open water forms ice only from MINIMUM_CONCENTRATION and
 * MINIMUM_VOLUME up, so water cools a little below freezing before its
 * first ice forms. When the area shrinks the lost
 * area's snow melts, its latent heat drawn from the water at freezing,
 * which freezes the same mass onto the ice, and its skin heat goes to
 * the ice volume; new area joins at the freezing point without snow.
 * A cell whose ice carries no concentration is fully covered.
 *
 * Snow lies on the ice: deposit() adds precipitation that falls on an
 * iced cell from air below the melting point, in water equivalent per
 * unit ice area, to `snow` — the ocean cells of the array the land
 * surface keeps its own snow in, when its buffer is shared. Snow
 * brightens the surface toward iceSnowAlbedo over iceFullSnow kg/m²,
 * conducts in series with the ice (snowConductivity over its depth at
 * snowDensity), melts before the ice does, and goes into the water when
 * the ice is gone, its latent heat drawn from the mixed layer. Snow
 * that falls on open water melts at once, cooling the water by that
 * latent heat; what falls on the leads of an iced cell freezes the same
 * mass of water onto the ice. The ocean's freshwater counts
 * precipitation when it falls, snow or not. Snow heavier than the ice's
 * freeboard (waterDensity − iceDensity per metre of ice) floods and
 * freezes into snow-ice: the surplus mass leaves `snow` and joins the
 * ice at iceDensity, which conserves both mass and energy.
 *
 * The bare ice's albedo depends on its skin temperature as CCSM3's sea
 * ice does (Briegleb et al. 2004): iceAlbedo while cold, falling linearly
 * over the last iceMeltingRange kelvin below the melting point to
 * meltingIceAlbedo, which stands for melting ice with its ponds; ice
 * thinner than fullAlbedoThickness blends toward the water's albedo. With
 * snowAgeing the snow on the ice has its own albedo in `snowAlbedo` (the
 * ocean cells of the land's array when its buffer is shared), aged by
 * agedSnowAlbedo toward iceSnowFloor and refreshed by snowfall
 * (refreshedSnowAlbedo); without it the snow's albedo is iceSnowAlbedo.
 *
 * With iceAlbedoPhases the ice follows the phases Perovich et al. (2002,
 * JGR 107(C10) 8044; the sequence of Perovich et al. 2011, Ann. Glaciol.
 * 52(57) 192) observed on SHEBA's multiyear ice: dry snow keeps its
 * albedo until it melts, ageing in the cold toward iceDrySnowFloor, and
 * ages wet toward iceMeltingSnowFloor only in a step in which its surface
 * melts (phasedSnowAlbedo); neither ageing brightens snow. Melting bare ice
 * is white ice with melt ponds over pondFraction of it (pondedIceAlbedo),
 * in place of meltingIceAlbedo: white ice 0.64 (MOSAiC, Light et al. 2022),
 * ponds 0.2 (0.1-0.3 observed), the fraction 0.25 the MODIS peak in July
 * (Ding et al. 2019; Rösel et al. 2012); dry snow 0.84 and melting snow
 * 0.77 as the HIRHAM scheme of Køltzow (2007) has them.
 *
 * With `iceLight` the share i0 = 0.18 (1 - c) + 0.35 c (FMM97, with c the
 * diffuse share of the sunlight reaching the surface, `diffuse` in update)
 * of the sunlight the ice absorbs passes its surface, none under snow and
 * blended over iceFullSnow as the albedo is, and leaves the skin's balance.
 * Of it, exp(-1.5 (h - 0.1)) reaches the water: the share transmission[i]
 * goes to penetrated[i] and the rest melts the base. What the ice holds
 * back returns to the surface, or with `brineReservoir`, once the skin is
 * at the melting point and while any is held, fills `brine` (J/m² of ice),
 * Semtner's reservoir, which holds at most half the latent heat of the ice
 * at the end of a step, the rest melting ice. It stands in for FMM97's
 * eq. 9: while it lasts and the skin is below the melting point, an upper
 * half-layer of the ice (its midpoint at h/4) is held at the melting
 * point, the skin taking (T_m - T) / (h / 4k + snow's resistance) and the base
 * (T_f - T_m) / (3h / 4k), and the reservoir pays the difference. With no
 * heat in the ice to hold, this only approximates eq. 9.
 *
 * With iceExchange 'turbulent' (ICE_EXCHANGE above) the ocean does not hand
 * the ice the heat of its top layer above freezing. The layer keeps it as
 * `waterHeat`, which the ocean writes after its step, the ice changes in its
 * steps and the ocean reads back before the next, and the ice draws on it at
 * the bottom and lateral rates of conductance[i] and w_lat. The leads' net
 * heat gain, sunlight included, goes into waterHeat and melts no ice; a
 * net loss first takes the heat above freezing and freezes new ice for what
 * remains, as the open water of CICE does; the atmosphere sees the leads at
 * the freezing point, an approximation where the water under them is warmer.
 * Lateral melt takes the fraction r of the area and, with it, of the volume
 * and snow, the lost area's snow and skin heat going where area loss sends
 * them, so that the thickness stays; the heat of what melts, the volume's
 * latent heat, comes from waterHeat. Ice that forms in open water leaves the
 * water at freezing and ice that melts away leaves it with waterHeat. Water
 * the ocean finds supercooled under ice freezes onto the base through
 * `oceanFlux`, as with the instant exchange, which then holds no heat for the
 * ice to melt. The ice's surface energy plus waterHeat changes by exactly the
 * surface flux plus `oceanFlux`. Over a cell with no ocean column under it
 * (`oceanCells` 0: a lake, or any cell of a model without an ocean) nothing is held and the ice takes the instant exchange.
 */
export const ICE_LIGHT = { clear: 0.18, cloudy: 0.35, extinction: 1.5, surfaceLayer: 0.1, reservoirShare: 0.5 };

export const LEAD_CLOSURES = ['squared', 'hibler'];

export function openWaterAlbedo(mu) {
  return 0.026 / (Math.pow(mu, 1.7) + 0.065) + 0.15 * (mu - 0.1) * (mu - 0.5) * (mu - 1);
}

/*
 * The snow-albedo ageing of Douville et al. (1995) as the ECMWF land
 * scheme carries it (Dutra et al. 2010, appendix eq. A7 and eq. 9): snow
 * whose surface lies within wetSnowRange kelvin of the melting point
 * relaxes toward `floor` at meltingSnowAgeing per day; colder snow loses
 * coldSnowAgeing per day, down to `floor`, scaled with ageingActivation
 * by exp(ageingActivation (1/T_melt − 1/T)), the temperature dependence of
 * grain growth in BATS (Dickinson et al. 1993), so that cold dry snow
 * keeps its brightness. A snowfall of `fall` kg/m² moves the albedo
 * min(1, fall / refreshSnowfall) of the way to freshSnowAlbedo.
 */
export const SNOW_AGEING = { freshSnowAlbedo: 0.85, coldSnowAgeing: 0.008, meltingSnowAgeing: 0.24, refreshSnowfall: 10, wetSnowRange: 2, ageingActivation: 5000 };
export function agedSnowAlbedo(albedo, temperature, dt, floor, { coldSnowAgeing, meltingSnowAgeing, wetSnowRange, ageingActivation } = SNOW_AGEING) {
  const days = dt / 86400;
  if (temperature >= MELTING_POINT - wetSnowRange) return floor + (albedo - floor) * Math.exp(-meltingSnowAgeing * days);
  const pace = ageingActivation > 0 ? Math.min(1, Math.exp(ageingActivation * (temperature - MELTING_POINT) / (MELTING_POINT * temperature))) : 1;
  return Math.max(floor, albedo - coldSnowAgeing * pace * days);
}
export function refreshedSnowAlbedo(albedo, fall, { freshSnowAlbedo, refreshSnowfall } = SNOW_AGEING) {
  return albedo + Math.min(1, fall / refreshSnowfall) * (freshSnowAlbedo - albedo);
}

export const ALBEDO_PHASES = { iceDrySnowFloor: 0.84, iceMeltingSnowFloor: 0.77, whiteIceAlbedo: 0.64, pondAlbedo: 0.2, pondFraction: 0.25 };
export function phasedSnowAlbedo(albedo, temperature, melting, dt, dryFloor, wetFloor, { coldSnowAgeing, meltingSnowAgeing, ageingActivation } = SNOW_AGEING) {
  const days = dt / 86400;
  if (melting) return Math.min(albedo, wetFloor + (albedo - wetFloor) * Math.exp(-meltingSnowAgeing * days));
  const pace = ageingActivation > 0 ? Math.min(1, Math.exp(ageingActivation * (temperature - MELTING_POINT) / (MELTING_POINT * temperature))) : 1;
  return Math.min(albedo, Math.max(dryFloor, albedo - coldSnowAgeing * pace * days));
}
export const pondedIceAlbedo = ({ whiteIceAlbedo, pondAlbedo, pondFraction } = ALBEDO_PHASES) => (1 - pondFraction) * whiteIceAlbedo + pondFraction * pondAlbedo;

export function createSeaIce(mesh, options = {}) {
  const {
  slabHeatCapacity = 2.1e7, skinHeatCapacity = 2e5, conductivity = 2.0, minimumThickness = 0.1,
  iceDensity = 917, latentHeatFusion = 3.34e5, oceanAlbedo = null, diffuseWaterAlbedo = 0.06, iceAlbedo = 0.62, meltingIceAlbedo = 0.48, iceMeltingRange = 1, fullAlbedoThickness = 0.5,
  iceSnowAlbedo = 0.75, iceFullSnow = 20, snowConductivity = 0.31, snowDensity = 300, waterDensity = 1026, leadClosing = 0.3, leadExchange = 10, areaLossFactor = 1, iceExchange = DEFAULT_ICE_EXCHANGE,
  leadClosure = 'squared', snowAgeing = true, iceSnowFloor = 0.70, freshSnowAlbedo = SNOW_AGEING.freshSnowAlbedo, coldSnowAgeing = SNOW_AGEING.coldSnowAgeing, meltingSnowAgeing = SNOW_AGEING.meltingSnowAgeing,
  refreshSnowfall = SNOW_AGEING.refreshSnowfall, wetSnowRange = SNOW_AGEING.wetSnowRange, ageingActivation = SNOW_AGEING.ageingActivation,
  heatCapacity = null, oceanCells = null, transmission = null, penetrated = null, buffers = null, keepSmallIce = false,
  iceLight = false, brineReservoir = false, iceAlbedoPhases = false, iceDrySnowFloor = ALBEDO_PHASES.iceDrySnowFloor, iceMeltingSnowFloor = ALBEDO_PHASES.iceMeltingSnowFloor,
  whiteIceAlbedo = ALBEDO_PHASES.whiteIceAlbedo, pondAlbedo = ALBEDO_PHASES.pondAlbedo, pondFraction = ALBEDO_PHASES.pondFraction,
  } = options;
  const unread = iceAlbedoPhases ? ['meltingIceAlbedo', 'iceSnowFloor', 'wetSnowRange'].filter((name) => options[name] !== undefined) : [];
  if (unread.length) console.warn(`iceAlbedoPhases leaves ${unread.join(', ')} unread: the phases set the melting ice's albedo and age the snow on the ice their own way`);
  if (brineReservoir && !iceLight) throw new Error('brineReservoir holds the light iceLight lets through the surface');
  if (!LEAD_CLOSURES.includes(leadClosure)) throw new Error(`leadClosure is ${LEAD_CLOSURES.join(' or ')}, not ${leadClosure}`);
  const C = mesh.nCells;
  const latent = iceDensity * latentHeatFusion, openShare = leadClosure === 'hibler' ? () => 1 : (A) => 1 - A;
  const removal = seaIceRemoval(keepSmallIce);
  const transfer = iceExchangeOptions(iceExchange), turbulent = transfer.turbulent, holds = (i) => turbulent && (!oceanCells || oceanCells[i] > 0);
  const budget = { frozen: 0, melted: 0, snowfall: 0, snowMelted: 0, snowIce: 0, leadFrozen: 0, lateralMelted: 0 };
  const oceanFlux = new Float64Array(buffers && buffers.oceanFlux ? buffers.oceanFlux : new SharedArrayBuffer(8 * C));
  const waterHeat = new Float64Array(buffers && buffers.waterHeat ? buffers.waterHeat : new SharedArrayBuffer(8 * C));
  const conductance = new Float64Array(buffers && buffers.conductance ? buffers.conductance : new SharedArrayBuffer(8 * C));
  const snow = new Float64Array(buffers && buffers.snow ? buffers.snow : new SharedArrayBuffer(8 * C));
  const concentration = new Float64Array(buffers && buffers.concentration ? buffers.concentration : new SharedArrayBuffer(8 * C));
  const snowAlbedo = new Float64Array(buffers && buffers.snowAlbedo ? buffers.snowAlbedo : new SharedArrayBuffer(8 * C));
  if (!(buffers && buffers.snowAlbedo)) snowAlbedo.fill(freshSnowAlbedo);
  const brine = brineReservoir ? new Float64Array(buffers && buffers.brine ? buffers.brine : new SharedArrayBuffer(8 * C)) : null;
  const ageing = { freshSnowAlbedo, coldSnowAgeing, meltingSnowAgeing, refreshSnowfall, wetSnowRange, ageingActivation };

  const waterAlbedo = (mu) => oceanAlbedo ?? (mu === null ? diffuseWaterAlbedo : openWaterAlbedo(mu));
  const meltingEnd = iceAlbedoPhases ? pondedIceAlbedo({ whiteIceAlbedo, pondAlbedo, pondFraction }) : meltingIceAlbedo;
  const bareIceAlbedo = (temperature) => iceAlbedo + (meltingEnd - iceAlbedo) * Math.min(1, Math.max(0, (temperature - MELTING_POINT + iceMeltingRange) / iceMeltingRange));
  function coverAlbedo(thickness, water, snowCover, temperature, snowy) {
    if (thickness <= 0) return water;
    const bare = water + (bareIceAlbedo(temperature) - water) * Math.min(1, thickness / fullAlbedoThickness);
    return bare + ((snowAgeing ? snowy : iceSnowAlbedo) - bare) * Math.min(1, snowCover / iceFullSnow);
  }

  /*
   * temperature is the ice's skin temperature (cold when not given) and
   * snowy the albedo of the snow on it (fresh when not given).
   */
  function albedo(thickness, mu = null, snowCover = 0, fraction = thickness > 0 ? 1 : 0, temperature = -Infinity, snowy = freshSnowAlbedo) {
    const water = waterAlbedo(mu);
    return fraction * coverAlbedo(thickness, water, snowCover, temperature, snowy) + (1 - fraction) * water;
  }

  function albedoContrast(thickness, mu = null, snowCover = 0, temperature = -Infinity, snowy = freshSnowAlbedo) {
    const water = waterAlbedo(mu);
    return coverAlbedo(thickness, water, snowCover, temperature, snowy) - water;
  }

  function cover(i, thickness) {
    return thickness > 0 ? (concentration[i] > 0 ? concentration[i] : 1) : 0;
  }

  function load(ice, saved = null, savedBrine = null) {
    for (let i = 0; i < C; i++) concentration[i] = ice[i] > 0 ? (saved && saved[i] > 0 ? Math.min(1, saved[i]) : 1) : 0;
    if (brine) for (let i = 0; i < C; i++) brine[i] = ice[i] > 0 && savedBrine ? savedBrine[i] : 0;
    else if (savedBrine) {
      let heat = 0, area = 0;
      for (let i = 0; i < C; i++) { area += mesh.areaCell[i]; if (ice[i] > 0) heat += mesh.areaCell[i] * concentration[i] * savedBrine[i]; }
      if (heat > 0) console.warn(`the saved sea ice holds ${(heat / area).toExponential(2)} J/m² of brine (a global mean), dropped: this run has no ice.brineReservoir`);
    }
  }

  function energy(surfaceT, ice, snowCover = 0, area = ice > 0 ? 1 : 0, reservoir = 0) {
    return ice > 0 ? area * (skinHeatCapacity * (surfaceT - FREEZING_POINT) - latent * ice - latentHeatFusion * snowCover + reservoir) : slabHeatCapacity * (surfaceT - FREEZING_POINT);
  }

  /*
   * Ice that melts away, whether by the thermodynamics or because transport
   * left too little of it: the cell is open water again, its temperature
   * `base` (freezing for a cell that held ice) plus `energyLeft`, the skin
   * heat and brine less the latent heat of the ice and snow, over the water's heat
   * capacity; `volume` and `snowMass` (m³ and kg) go to the budget.
   */
  function meltAway(surfaceT, ice, i, energyLeft, volume, snowMass, base = FREEZING_POINT) {
    budget.melted += volume;
    budget.snowMelted += snowMass;
    surfaceT[i] = base + energyLeft / (heatCapacity ? heatCapacity[i] : slabHeatCapacity);
    ice[i] = 0;
    snow[i] = 0;
    snowAlbedo[i] = freshSnowAlbedo;
    concentration[i] = 0;
    if (brine) brine[i] = 0;
    if (turbulent) waterHeat[i] = 0;
  }

  const capacityOf = (i) => (heatCapacity ? heatCapacity[i] : slabHeatCapacity);
  /*
   * The temperature of the water under the ice of cell i: the freezing point,
   * and above it by waterHeat over the layer's heat capacity when the ice
   * exchanges heat with the water as turbulence does.
   */
  const waterTemperature = (i) => (holds(i) ? FREEZING_POINT + waterHeat[i] / capacityOf(i) : FREEZING_POINT);
  const holdWater = (i, temperature) => { if (holds(i)) waterHeat[i] = capacityOf(i) * (temperature - FREEZING_POINT); };

  function deposit(i, amount, airTemperature, ice, surfaceT) {
    if (airTemperature >= MELTING_POINT || amount <= 0) return false;
    budget.snowfall += mesh.areaCell[i] * amount;
    if (ice[i] > 0) {
      const area = cover(i, ice[i]), leads = (1 - area) * amount;
      snow[i] += amount;
      ice[i] += leads / (iceDensity * area);
      snowAlbedo[i] = refreshedSnowAlbedo(snowAlbedo[i], amount, ageing);
      budget.snowMelted += mesh.areaCell[i] * leads;
      budget.frozen += mesh.areaCell[i] * leads / iceDensity;
    } else { surfaceT[i] -= latentHeatFusion * amount / (heatCapacity ? heatCapacity[i] : slabHeatCapacity); budget.snowMelted += mesh.areaCell[i] * amount; }
    return true;
  }

  function update(surfaceT, ice, flux, i, dt, contrast = 0, solar = 0, diffuse = 0) {
    const ocean = oceanFlux[i];
    const capacity = heatCapacity ? heatCapacity[i] : slabHeatCapacity;
    const cellArea = mesh.areaCell[i];
    if (ice[i] <= 0) {
      concentration[i] = 0;
      if (snow[i] > 0) { surfaceT[i] -= latentHeatFusion * snow[i] / capacity; budget.snowMelted += cellArea * snow[i]; snow[i] = 0; }
      snowAlbedo[i] = freshSnowAlbedo;
      const pen = transmission ? solar * transmission[i] : 0;
      surfaceT[i] += dt * (flux[i] + ocean - pen) / capacity;
      if (transmission) penetrated[i] += dt * pen;
      if (surfaceT[i] < FREEZING_POINT) {
        const volume = (FREEZING_POINT - surfaceT[i]) * capacity / latent, area = Math.min(1, volume / leadClosing);
        if (area >= MINIMUM_CONCENTRATION && volume >= MINIMUM_VOLUME) {
          ice[i] = volume / area;
          concentration[i] = area;
          surfaceT[i] = FREEZING_POINT;
          if (brine) brine[i] = 0;
          if (turbulent) waterHeat[i] = 0;
          budget.frozen += cellArea * volume;
        }
      }
      return;
    }
    const h = ice[i], A = cover(i, h);
    let T = surfaceT[i], s = snow[i];
    const split = contrast - leadExchange * (FREEZING_POINT - T);
    let iceFlux = flux[i] - (1 - A) * split;
    const waterFlux = flux[i] + A * split;
    let through = 0, held = 0;
    if (iceLight) {
      const passed = Math.max(0, solar - (1 - A) * contrast) * (ICE_LIGHT.clear + (ICE_LIGHT.cloudy - ICE_LIGHT.clear) * diffuse) * Math.max(0, 1 - s / iceFullSnow);
      through = passed * Math.exp(-ICE_LIGHT.extinction * Math.max(0, h - ICE_LIGHT.surfaceLayer));
      held = brine && (T >= DEVICE_MELTING_POINT || brine[i] > 0) ? passed - through : 0;
      iceFlux -= through + held;
    }
    const depth = Math.max(h, minimumThickness), insulation = s / (snowDensity * snowConductivity);
    let conduction = (FREEZING_POINT - T) / (depth / conductivity + insulation), base = conduction;
    let Q = brine ? brine[i] + dt * held : 0;
    if (Q > 0 && T < DEVICE_MELTING_POINT) {
      const upper = (MELTING_POINT - T) / (depth / (4 * conductivity) + insulation), lower = (FREEZING_POINT - MELTING_POINT) / (3 * depth / (4 * conductivity));
      const cost = dt * (upper - lower), part = Math.min(1, Q / cost);
      conduction += part * (upper - conduction);
      base += part * (lower - base);
      Q = part < 1 ? 0 : Q - cost;
    }
    T += dt * (iceFlux + conduction) / skinHeatCapacity;
    let thickness = h + dt * (base - ocean) / latent;
    if (iceLight) {
      const share = transmission ? transmission[i] : 0;
      thickness -= dt * through * (1 - share) / latent;
      if (transmission) penetrated[i] += dt * A * through * share;
    }
    const melting = T > MELTING_POINT;
    if (melting) {
      let excess = (T - MELTING_POINT) * skinHeatCapacity;
      T = MELTING_POINT;
      const fromSnow = Math.min(s, excess / latentHeatFusion);
      s -= fromSnow;
      budget.snowMelted += cellArea * A * fromSnow;
      excess -= fromSnow * latentHeatFusion;
      thickness -= excess / latent;
    }
    let leadHeat = (1 - A) * (waterFlux + ocean) * dt, excess = 0, side = 0, lost = 0;
    if (holds(i)) {
      const stored = waterHeat[i] + leadHeat, rise = Math.max(0, stored) / capacity;
      leadHeat = Math.min(0, stored);
      excess = Math.max(0, stored);
      const wantBasal = conductance[i] * rise * dt, lateral = Math.min(1, transfer.lateralRate * Math.pow(rise, transfer.lateralExponent) * dt), wantSide = lateral >= COVER_GRAIN ? lateral : 0;
      const wanted = A * wantBasal + latent * wantSide * A * Math.max(0, thickness);
      const share = wanted > excess ? excess / wanted : 1;
      side = share * wantSide;
      thickness -= share * wantBasal / latent;
      lost = side * A * Math.max(0, thickness);
      excess = Math.max(0, excess - A * share * wantBasal - latent * lost);
    }
    const grown = A * (thickness - h), leadMelt = Math.max(0, leadHeat / latent), leadIce = Math.max(0, -leadHeat / latent);
    const melted = Math.max(0, -grown) + leadMelt;
    budget.leadFrozen += cellArea * leadIce;
    budget.lateralMelted += cellArea * (leadMelt + lost);
    let volume = A * thickness - lost - leadHeat / latent;
    const energyLeft = A * skinHeatCapacity * (T - FREEZING_POINT) - latent * volume - latentHeatFusion * A * s + A * Q;
    const area = Math.min(1, (A - areaLossFactor * melted / (2 * h)) * (1 - side) + openShare(A) * leadIce / leadClosing);
    let skin = T, snowOnIce = s, kept = Q;
    if (area < A) volume += (A - area) * (s / iceDensity - (skinHeatCapacity * (T - FREEZING_POINT) + Q) / latent);
    else if (area > A) { skin = FREEZING_POINT + A * (T - FREEZING_POINT) / area; snowOnIce = A * s / area; kept = A * Q / area; }
    if (brine && area > 0) {
      const over = Math.min(kept, Math.max(0, kept - ICE_LIGHT.reservoirShare * latent * Math.max(0, volume) / area) / (1 - ICE_LIGHT.reservoirShare));
      kept -= over;
      volume -= area * over / latent;
    }
    if (area < removal.concentration || volume < removal.volume) {
      meltAway(surfaceT, ice, i, energyLeft, cellArea * A * h, cellArea * A * s, turbulent ? FREEZING_POINT + excess / capacity : FREEZING_POINT);
      return;
    }
    if (area < A) budget.snowMelted += cellArea * (A - area) * s;
    budget.frozen += cellArea * (volume - A * h);
    const newThickness = volume / area;
    const flooded = Math.max(0, snowOnIce - (waterDensity - iceDensity) * newThickness) * iceDensity / waterDensity;
    if (turbulent) waterHeat[i] = excess;
    surfaceT[i] = skin;
    snow[i] = snowOnIce - flooded;
    snowAlbedo[i] = snow[i] > 0 ? (iceAlbedoPhases ? phasedSnowAlbedo(snowAlbedo[i], skin, melting, dt, iceDrySnowFloor, iceMeltingSnowFloor, ageing) : agedSnowAlbedo(snowAlbedo[i], skin, dt, iceSnowFloor, ageing)) : freshSnowAlbedo;
    budget.snowIce += cellArea * area * flooded;
    ice[i] = newThickness + flooded / iceDensity;
    concentration[i] = area;
    if (brine) brine[i] = kept;
  }

  return {
    albedo, albedoContrast, cover, load, energy, update, deposit, meltAway, budget, removal, slabHeatCapacity, skinHeatCapacity, freshSnowAlbedo, latent, latentHeatFusion, leadClosing, leadExchange, oceanFlux, snow, snowAlbedo, concentration, brine, iceLight,
    iceExchange: transfer, turbulent, waterHeat, conductance, waterTemperature, holdWater,
    shared: { oceanFlux: oceanFlux.buffer, snow: snow.buffer, snowAlbedo: snowAlbedo.buffer, concentration: concentration.buffer, waterHeat: waterHeat.buffer, conductance: conductance.buffer, ...(brine ? { brine: brine.buffer } : {}) },
  };
}
