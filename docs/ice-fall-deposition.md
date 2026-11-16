# Deposition-limited ice fall (`moist.iceDeposition`)

Under the ramp (the default), a cloud layer's ice share is 1 − α(T), linear from 0 at liquidTemperature (273.15 K) to 1 at iceTemperature (235.15 K). That share falls at the Heymsfield–Donner speed, and the liquid share α converts to rain.

With `iceDeposition` (false by default; it needs `iceFall`), a layer between those temperatures handles its cloud differently:
- What falls is the layer's cloud ice q_i, which only deposition, the ice falling in and the ice share of post-plume condensate grow, following the Wegener–Bergeron–Findeisen process as the IFS models it.
- The liquid q_l = q_c − q_i converts to rain at the plain liquid rates.

Saturation, latent heat and the radiation keep the ramp. The option changes how a layer's condensate is split between the sink that falls and the sink that rains out.

## Equations

These follow IFS Documentation Cy41r2 Part IV §7.2.4(d), eqs 7.33–7.43 (after Rotstayn et al. 2000 and Meyers et al. 1992). Units are SI throughout: T in K, p and e in Pa, ρ in kg/m³, q in kg/kg, N_i in m⁻³, Δt in s, and ρ_i^(1/3) = 8.88 kg^(1/3)/m. The model's R is 287.06 J/(kg K) and R_v = R/0.622 = 461.5 J/(kg K).

Each step, in every layer with iceTemperature < T < liquidTemperature and q_l > 0, after the ice falling in has been added:

- Vapor pressures: e_sl is Bolton's (`saturationVaporPressure`) and e_si the IFS form (`iceVaporPressure`). The air density is ρ = p/(R T).
- Active nuclei: N_i = 1000 exp[12.96 (e_sl − e_si)/e_si − 0.639] (7.40). This is the IFS text's and Meyers's form, and with it f_nuc reaches 1 at −23 °C as the text says.
  - CLOUDSC divides by e_sl instead. The e_si form has 1.35× the nuclei at 258 K, 1.7× at 253 K and 2.25× at 248 K, so it glaciates 1.2–1.7× faster than the operational code at 248–258 K and lifts F_dep earlier.
  - Bolton's e_sl over the IFS e_si stands in for CLOUDSC's min(Koop, e_l/e_i), which is e_l/e_i above about 238 K, to about 2 % in the supersaturation.
- Conduction and diffusion terms:
  - A″ = (L_s/(K_a T)) (L_s/(R_v T) − 1)
  - B″ = R_v T p/(2.21 e_si)
  - L_s = L + L_f and K_a = 2.4e-2 W m⁻¹ K⁻¹ (7.34, 7.35).
- Growth coefficient: c_svd = 7.8 (N_i/ρ)^(2/3) (e_sl − e_si)/(ρ_i^(1/3) (A″ + B″) e_si), with ρ_i = 700 kg/m³ (7.38).
- In-cloud ice: a is the layer's cover, the uniform distribution's cover that the fall speed also uses. The in-cloud ice, with the seed as its floor, is q₀ = max(q_i/a, m₀ N_i/ρ), where m₀ = 1e-12 kg (7.41).
- New in-cloud ice: q₁ = ((2/3) c_svd Δt + q₀^(2/3))^(3/2). This is exact over the step for a constant c_svd (7.39).
- Deposition: D = F_dep · min(max(0, a (q₁ − q₀)), q_l). D moves from q_l to q_i; q_c is unchanged.
- Cloud-top reduction:
  - F_dep = min[f_nuc + (1 − f_nuc)(0.1 + Δz_top/500 m), 1], with f_nuc = min(N_i/15000, 1) (7.42, 7.43).
  - The tracker's cover is the uniform distribution's cover of the layer's q_c at the start of its visit in the rain step, before the ice falling in.
  - Δz_top is 0 in a cloud-top layer: one whose cover is at least 0.01 under a layer with less. The column's top layer is a top whenever its cover is at least 0.01.
  - Below a top, Δz_top grows by each layer's own thickness R T Δσ/(g σ), as CLOUDSC's ZCLDTOPDIST does.

Time to grow 0.1 g/kg of in-cloud ice from the seed (a = 1, F_dep = 1, 600 hPa):

| T | 268 K | 263 K | 258 K | 253 K | 248 K | 240 K (400 hPa) |
|---|---|---|---|---|---|---|
| N_i (m⁻³) | 1050 | 2110 | 4380 | 9410 | 20900 | 79300 |
| Time | 171 min | 68 min | 37 min | 23 min | 15.6 min | 5.2 min |
| F_dep in the top layer | 0.16 | 0.23 | 0.36 | 0.66 | 1 | 1 |

- The times grow about as p^1.5: at 248 K they are 8.7, 11.9, 15.6 and 24.0 min at 400, 500, 600 and 800 hPa.
- They grow as q^(2/3) of the target: a thin cloud of 0.02 g/kg glaciates 2.9× faster than the table, and one of 0.01 g/kg 4.6× faster, often within one step.
- With F_dep < 1 they lengthen by about 1/F_dep.

## The memory

`moist.cloudIceShare` is each layer's glaciated fraction φ = q_i/q_c, K·C values in [0, 1]. It is a Float64Array on the CPU (in a shared buffer for the workers) and PH region `CICE` (f32) on the GPU. It is allocated only with the option, so the option-off layout and saves are unchanged.

Between steps the memory holds φ; within the moist step the same array holds the ice mass q_i:

1. **Start.** Before the first condensation pass, q_i = φ q_c. On the CPU this is the top of `adjust`; on the GPU it is `adjustMix`, after the mixing and before its saturation pass, the same point.
2. **Condensation passes.** New large-scale condensate is liquid between 235.15 and 273.15 K and ice at or below 235.15 K, through the limits below. Condensate formed in the post-plume pass enters with the ramp's ice share (see Detrainment). A loss of cloud takes its liquid first: q_i is capped at q_c when it is next read.
3. **Rain step**, top down, at each layer's visit:
   - the limits: q_i = q_c at or below iceTemperature, q_i = 0 at or above liquidTemperature, and q_i ≤ q_c between;
   - the ice falling in (see Arriving ice), deposition, the fall, and the liquid's conversion.
4. **End.** After the third condensation pass and the filler, a single cap and the store: φ = 0 where q_c ≤ CLEAR_AIR (1e-7 kg/kg), 1 at or below iceTemperature, 0 at or above liquidTemperature, and min(q_i, q_c)/q_c between.

Within the step only these change q_i: deposition, the ice kept from what falls in, the ice share of post-plume condensate, the fall, and the caps and limits.

The memory is not advected. Cloud that moves keeps its share, because each step starts from q_i = φ q_c of the cell the cloud is in.
- A memory of ice mass would err both ways: an arrival or growth of q_c would be all liquid, and a thinning would keep the ice until the cap.
- With φ two errors remain. Old cloud that arrives in a clear cell (φ = 0) is liquid and restarts its glaciation, and φ has gradients across cell boundaries.
- Neither touches a stationary Arctic deck. For storm-track and frontal cloud at 258–273 K, N=64 arms bound the error from below.
- A fourth advected tracer would still need the cap, since independent limiters break q_i ≤ q_c, and would add layout, mixing, plume and conservation work on both engines.

## Detrainment

The model has no detrained-condensate stream: the plume writes only θ and q, and the saturation pass after it (`condenseColumn`'s second call on the CPU, `adjustPlume`'s on the GPU, in columns where a plume ran) makes the condensate.

- `detrainedPhase: 'ramp'` (the default) gives that pass's positive change the ramp's ice share: q_i = min(max(0, q_i), q_c) + (1 − α(T)) Δq_c, with q_c and T the layer's cloud and temperature after the plume. The cap matters where the first pass evaporated cloud and left q_i above q_c; without it, the regrown cloud would all count as ice. The plume's own ice is the same ramp share at its own temperature.
- `detrainedPhase: 'liquid'` makes it liquid, for an attribution arm. It needs `iceDeposition`. In both settings the pass first caps a layer's stored ice at the cloud the first pass left, so a regrowth never counts as the stale ice; the ramp then adds its share of the regrowth and 'liquid' adds none.
- `iceDeposition` runs only with `rainEvaporationLaw` 'deficit' and no `rainEvaporationCeiling`, and both engines refuse the others. Under 'kessler' with `rainCover` 'overlap', at N=16, the engines part by 1e-2 K where ice falls into rain below an ice cloud.
- The model's ramp is linear over 235.15–273.15 K. The IFS's (eq 7.6) is quadratic over 250.16–273.16 K, so this entry is less icy than the IFS's at 250–268 K. Deposition makes up the difference within the times above.

## Arriving ice

Ice falling into a layer below liquidTemperature joins its q_i, less what passes on through the implicit fall. It melts, into the precipitation, only at or above liquidTemperature, without heat in the column: the surface books L_f on the precipitation, as in the default. The ice leaving the lowest layer is large-scale precipitation.

## The liquid a layer holds back

Held liquid has three sinks:
- deposition onto the ice, which then falls;
- Kessler autoconversion and the lifetime at the plain rates: threshold 2e-4, 1e-3 /s, 1 h, or 3 h for stratiform cloud;
- evaporation.

There is no Bergeron factor under the option:
- Cy41r2's large-scale liquid autoconversion has none (7.49–7.52); the factor survives only in its convection (6.36–6.38, where 6.38 defines it). Its default autoconversion is Khairoutdinov–Kogan; the model's Kessler stands in for the Sundqvist form.
- Deposition is the Bergeron–Findeisen process here, so keeping `bergeronAutoconversion` would count it twice at 235–255 K.
- `bergeronAutoconversion` therefore has no effect anywhere under the option, on both engines: the CPU gates it on `bergeronAutoconversion && !iceDeposition`, and the GPU kernel leaves out its patch.
- From 258 to 263 K the factor (2.6 to 2.1) converts faster than deposition glaciates (23–28 min against 37–68 min), so the liquid there lives 2.1–2.6× longer than in the default. At 268 K the factor is 1.2: conversion takes 50 min against 171 for glaciation, and the liquid lives 1.2× longer. That is intended: it is the Arctic gain, and also the SWCRE risk.
- The model has no riming sink (IFS eq 7.56, S_rime in eqs 7.7 and 7.10). If the liquid overshoots, at 250–450 hPa or in the Arctic SWCRE, riming is the next lever, not the factor.

coldfall64 (`iceFallColdOnly`) converted its held cloud with the factor: 1 + 0.5 √18 = 3.1 at and below 250 K, its cap, so the whole cloud decayed in about a third of the lifetime. That this caused its LWCRE loss is a lead, which the attribution arms bfoff64 and coldbf64 test.

What I expect under the option:
- **235–255 K:** deposition glaciates 0.1 g/kg within 5–30 min, so the ice falls much as in the default. Tropical mid-level cloud should stay near ctl64, which is the LWCRE guard. The 30S–30N 250–450 hPa liquid is the risk: 12.3 g/m² without the factor and 6.7 with it in six-hour uncoupled pairs, against ERA5's about 1.
- **255–273 K:** the liquid lives longer than in the default. Boundary-layer mixed-phase cloud should hold more supercooled liquid, which is the Arctic rows.
- **Risk:** global SWCRE gets stronger. It is already −51 W/m² against CERES's −45.

## Thermodynamics and radiation

- The heat released by condensation or deposition is L, whatever the phase. L + (1 − α) L_f is the latent heat in the saturation slope only; see "Where the heat of fusion is" in docs/cloud-gas-overlap.md.
- Deposition's L_f is not booked: about 0.1 K/day in a cloud layer and 2–4 W/m² over the column. The memory carries no enthalpy, so energy closes by construction.
- Saturation stays on the ramp. A layer of mostly liquid cloud is held 1–10 % below water saturation, which inflates its held liquid at 253–263 K by roughly 0.05–0.1 g/kg. If the Arctic rows fall short, a saturation on the memory's share is the next lead.
- The radiation's cloud phase stays on the ramp unless `radiation.cloudPhase` is 'memory' (see "Arm B: radiation on the memory").
- The plumes keep their own ice.

## CPU (f64) and WGSL (f32)

- Both engines use the same expressions, with max(0, ·) on the growth. q₀ ≥ the seed > 0, so f32 `pow` is defined. The subtraction q₁ − q₀ loses about 6e-8 of q_i.
- The GPU applies three string patches, only with the option. Each anchor must appear exactly once:
  - `cloudIceRain` on `adjustRain`, in place of `bergeronRain` and before `deepStoreRain`;
  - `cloudIceMix` on `adjustMix`;
  - `cloudIcePlume` on `adjustPlume`, under 'ramp'.
- With the option off the kernel texts and the PH layout are unchanged. `CICE` sits immediately before `T2M`, and every field from `T2M` on, the turbulent exchange's among them, moves by K × C.
- The GPU runs one thread per column, top down. It carries the cover above and Δz_top as scalars, and reads and writes the memory in the PH buffer, so there are no private arrays or pointers (CLAUDE.md's M5 and M1 rules).
- Four decisions can split the engines at a few cells: the cloud-top test (cover 0.01), the two temperature limits, min(·, q_l), and the CLEAR_AIR floor.

## Saves and restores

- `placeCloudIce(model, saved)` places a saved `cloudIceShare` of the model's length. Otherwise, as for a fresh start, it lays the ramp's share 1 − α(T) of the model's state, 0 where q_c ≤ CLEAR_AIR. It runs after the state is in place.
- Callers: scripts/spinup.mjs (resume and fresh start), the page's `placeDeck`, explainOcean, fixedSst, figureState, deckGate, shortwaveBudget, explainSurfaceBudget, adjustStages, tropicalHeating and compareEngines.
- Writers: spinup, explainOcean and fixedSst write `cloudIceShare` with the other moist memories, and the page saves it in `savedCumulus`.
- remapState drops it, so a remapped state starts from the ramp.
- The GPU model mirrors `CICE` into `moist.cloudIceShare` when it reads back, and uploads it with the retained physics on load. compareEngines maps `CICE` to it.

## Tests

- `test/iceDeposition.test.mjs`:
  - the defaults and the option's checks;
  - columns at the limits bit for bit;
  - a hand-computed mixed-phase cloud-top layer;
  - the cover step across 0.01 and a q_l-limited layer;
  - the glaciation times;
  - doubling q_c between steps doubles q_i;
  - the Bergeron factor has no effect under the option;
  - the post-plume cap by hand: a cloud all ice that the first pass halves and the second regrows by 0.3 at α = 0.5 holds 0.8125 of its water as ice;
  - the same cap on the GPU: after adjustMix and adjustPlume, with every mixed-phase cloud all ice and `cumulusMemory` 0, the ice matches the CPU's stages where the post-plume pass regrows cloud by more than 1e-7 kg/kg under stale ice, and the hand example set above a deep plume gives 0.8125 on both engines;
  - `placeCloudIce`, and a GPU model's ramp start from a fresh start, a save without the fraction and one of another length;
  - the option-off kernel text and layout;
  - the source lint over restorers and writers;
  - GPU parity at N=16 on bl34 and bl48, at two time steps, with `detrainedPhase` 'ramp' and 'liquid' and with `bergeronAutoconversion`;
  - the refusal of other rain-evaporation laws and of a ceiling, on each engine;
  - a GPU core's refusal of a glaciated fraction whose length is not K × C.
- `test/gpuReload.test.mjs`: under the option a load restores every device buffer word for word.
- `test/pageRestore.test.mjs`: under the option a page save restores and steps on as the run that did not stop.
- `test/spinupRestart.test.mjs`: under the option a split spin-up ends bit for bit on the whole one.
- `test/atmosphereConservation.test.mjs`: the moist row "over a slab with iceDeposition".
- CLAUDE.md's M5 checks: `test/boundaryLayer.test.mjs`, `test/gpuModel.test.mjs`, and compareEngines at the page state with and without the option.

## Open items

- `rainEvaporationLaw` 'kessler' with `rainCover` 'overlap': the engines part where ice falls into rain below an ice cloud, in the `iceAbove` rule for the rain's cover. The case is `test/iceDeposition.test.mjs`'s random state (seed 4242) at N=16 on bl34 with a 1350 s step. There the θ gap moves by 1.1e-2 K beside the engines' gap without the option, at column 1017, layers 19–20, below the ice cloud at k17. Both engines refuse the combination until the rule matches.

## Arm B: radiation on the memory

`radiation.cloudPhase` is 'ramp' (the default) or 'memory'; 'memory' needs `moist.iceDeposition` and throws without it, in `createModel`, `createGpuCore` and at the first radiation column. It also needs `cloudOptics`'s optics, and both engines refuse it with a gray `cloudScattering` or `cloudAbsorption`.

The rule: a layer's resolved cloud water q_c has liquid share 1 − φ between iceTemperature and liquidTemperature, and the ramp's share outside them, where the moist physics ignores φ (ice at or below iceTemperature, liquid at or above liquidTemperature). Its cumulus water keeps the ramp's α(T). The optics per kilogram are linear in the liquid share, so the layer's combined water has the water-weighted mean share (q_c (1 − φ) + w_cu α) / (q_c + w_cu), with the masses and the clear-sky mask applied as the radiation applies them. The deck's water (the mixed-layer or stratus deck, its home layer's unit) stays on α(T).

Sites that decide the radiation's cloud phase:

| Site | CPU (`js/physics/radiation.module.js`) | GPU |
|---|---|---|
| Layer optics, solar depth and longwave absorption | `cloudOptics(T, continental, options, out, share)`, called a second time after the cumulus water is added; the ramp call stays | `cloudOpticsAt(T, continental, liquid)`, with `cloudOptics` its ramp wrapper; radiation kernel |
| Layer optics in the deck's shading loop | same arrays | physics kernel, `MEMORY_SHADE` |
| Layered shortwave's liquid share | `layerLiquid[k]` set from the mixed share | `liquid` in `layeredKernel`, `MEMORY_LAYERED` |
| Deck unit at its home layer | `rampSolar`, in `shadowDeck` | `light.unit` already ramp; `unit` in the shading loop and `deckUnit` in the radiation kernel read the ramp optics |
| Deck water over the layer's longwave and shortwave | `rampInfrared`, `rampSolar` at `deckK` | `rampOptics.y`, `deckUnit` |

Under 'ramp' every kernel text, the shared functions and the PH layout are those of the default. Under 'memory' only the physics, radiation and layered kernels change, by string patches that each anchor exactly once, and the shared `cloudOptics` is split into `cloudOpticsAt` and its wrapper. The kernels read `PH_CICE`, which exists under `iceDeposition`.

Tests are in `test/cloudPhaseOptics.test.mjs`: the defaults and throws (including a model without moisture), the CPU optics at 255 K with φ = 0.2, with cumulus water and outside the ramp's limits, the mixed-layer deck on both `deckCloud` modes (the deck takes the ramp, the resolved water above it 1 − φ), kernel identity, and N=16 CPU against GPU rows on bl34 and bl48 with cumulus water in layers that also hold resolved water, and on bl48 with the mixed-layer deck on (a deck on most columns on both engines). The parity rows bound the engines' gap on the change the memory makes at 5 % of that change. The deck rows leave out columns where the engines took different deck decisions or had parted states before the last radiation call, and bound the gap layer by layer. Resolved water of a trace in the deck's own layer, with limits set around its temperature, must leave the deck's water and sunlight those of the ramp, on both `deckCloud` modes: that row fails with the shadow's or the final fluxes' deck unit on the memory-mixed optics. `memory` is refused when moist's `liquidTemperature` or `iceTemperature` differs from radiation's, and without moisture.

What I expect from arm B, against dep64 (iceDeposition on the ramp) from the same fresh start and the same tree, written before the run. Both arms run from one export of ecf628e4 and differ only in `RADIATION`.

Measured inputs, not predictions:

- φ on the first day of icedep/f7, a CPU N=64 run of dep64's configuration (iceDeposition, radiation on the ramp) from ctl64's day-30 save, exceeds the ramp's ice share 1 − α in every latitude band and temperature bin. At 253–263 K it is 0.66–0.94 against 0.38–0.39; at 263–273 K, 0.42–0.76 against 0.12–0.15. So resolved cloud between the limits reads icier under 'memory' than under the ramp.
- Per kilogram, ice has 0.4–0.7 of liquid's solar depth and 0.2–0.3 of its longwave absorption.
- ctl64 (seaice/code at 44021ee3, without iceDeposition), from its monthly saves:
  - 70–90N ice, surface shortwave cloud effect (all-sky minus clear-sky sunlight at the surface): −72.7 W/m² in the day-61 save (Apr 21–May 20, model days 32–61) and −58.1 in the day-426 save. These are snapshots of each save's state, not 30-day means: `scripts/arcticSunlight.mjs` builds the CPU model from the save, holds the state as saved, lights each column north of 60N at the sun of the save's window (8 instants a day on the days 25, 15 and 5 before the save), and runs the radiation cell by cell, once as is and once with the resolved cloud, cumulus and deck removed. The ice70 cells are the sea cells north of 70N whose saved concentration exceeds 0.15, weighted by area. Day 61 is 187.6 all-sky minus 260.3 cloud-free of offline downward surface sunlight, day 426 is 209.0 minus 267.1.
  - Days 395–730: global SWCRE −48.1 and LWCRE 23.9 W/m²; 30–60° in both hemispheres −51.2 and 23.6; 30S–30N −52.2 and 25.8; net TOA +3.6 against CERES's +1.1, day means at the save days. The cloud effects are the means over ctl64's 12 saves in that range of `scripts/cloudEffectBands.mjs`, which area-weights each save's meanShortwaveCloudEffect and meanLongwaveCloudEffect (absorbed minus clear-sky absorbed sunlight, and clear-sky minus all-sky outgoing longwave). Those are day means, of the day each save closes, not 30-day means (the radiation sums reset at every daily diagnostics pass). The registered rows are the mean over days 395–1095 for every arm and for ctl64 alike, estimated as the mean of the saves' day means, one day every 30. The bands: global; 30S–30N is |lat| < 30; 30–60 is 30 ≤ |lat| < 60 with both hemispheres pooled; 60–90 is |lat| ≥ 60, by cell-center latitude. The net TOA is the mean of the saves' day means of absorbed sunlight minus outgoing longwave (meanAbsorbedSolar and meanOutgoingLongwave, also day means), from icedep/tools/guards.mjs (scorecard64.mjs reads the same quantity per day).
  - Arctic September extent (the day-183 save): 0.49 million km² against NSIDC's 5.91.
- Tropical φ at 263–273 K, about 0.75 against a ramp of 0.14, is set by ice that falls in from above and by the memory, which is not advected; local deposition does not set it, since there the caps remove about as much ice as deposition adds (1.82, 1.46, 1.35, 1.27, 1.20 and 0.98 times as much over f7's days 1–6). This bin, not the Arctic, carries the global SWCRE prediction, so read the tropical and storm-track rows separately.

The predictions as registered at 3049ac24:

- **255–273 K cloud:** its SWCRE weakens. Global SWCRE weakens by 1 to 4 W/m², toward CERES.
- **Thin 235–255 K cloud:** its LWCRE weakens. Global LWCRE weakens by 0.3 to 1.5 W/m².
- **Arctic:** the spring surface shortwave cloud effect weakens by 2 to 6 W/m², the wrong way for the melt-out rows. September extent ends below dep64's.
- **Net TOA:** 1 to 3 W/m² more absorbed than in dep64.

The rows these are scored on, added on Nov 15 before the run. Δ is arm B minus dep64. The global and band rows are means over days 395–1095, each the mean of the saves' day means at the save days: a save's meanShortwaveCloudEffect, meanLongwaveCloudEffect, meanAbsorbedSolar, meanOutgoingLongwave and meanPlanetaryAlbedo are day means, and the rows are estimated the same way for every arm and for ctl64 (the window quantity of the Arctic reading, from means.shortwave and means.concentration, stays a true 30-day mean), and the Arctic row is the mean of the day-61, 426 and 791 saves. dep64's expected value is ctl64's, since iceDeposition's own shift is not yet known. Every W/m² range is my conversion of the φ gap through the optics ratios above, a judgment that no radiation run has checked.

| Row | dep64 expected | Target | Predicted Δ | Falsified if |
|---|---|---|---|---|
| 70–90N ice surface SW cloud effect, Apr 21–May 20 (model days 32–61) | −65 ± 8 (ctl64 −72.7, −58.1) | ERA5 −61 | +2 to +6: toward ERA5 from below, and more sunlight on the ice for the melt-out rows | Δ ≤ 0 |
| Global SWCRE | −48 ± 2 | CERES −45.1 | +1 to +4, toward CERES | Δ < 0 |
| Global LWCRE | 24 ± 1 | CERES 25.6 | −0.3 to −1.5, away from CERES | Δ > 0 |
| 30S–30N SWCRE, LWCRE | −52, 26 | none on disk (the CERES file holds global means only) | +1 to +5; −0.3 to −1.5 | ΔSW < 0 or ΔLW > 0 |
| 30–60° SWCRE, LWCRE (storm tracks) | −51, 24 | none on disk | +1 to +5; −0.2 to −1.5 | ΔSW < 0 or ΔLW > 0 |
| Arctic September extent | 0.5 to 2.5 million km² (ctl64 0.49; iceFallColdOnly 2.55 in year 1) | NSIDC 5.91 | −0.2 to −1.0, and never below 0 | Δ > +0.2 |
| Net TOA | +3.6 ± 1.5 | CERES +1.1 | +1 to +3 as registered; the SW and LW ranges together allow −0.5 to +3.7 | Δ ≤ 0 |
| Inert | | | | \|ΔSWCRE\| < 0.3 and \|ΔLWCRE\| < 0.1 W/m² globally: the coupling did nothing, and the pair says nothing about it |

Amendment, Nov 15, written at 06:42 (dacee882) and restated in c21f2603 at 07:28, before the run, not a measurement. It comes from icedep/f7 (log ~/.cache/webgcm-scratch/icedep/f7/run.log): a CPU probe of dep64's configuration, not arm B's, started from ctl64's day-30 save. Its days 1–6 therefore fall in the Apr 19–May 19 window that the Arctic row scores. Its 60–90N φ at 253–263 K went 0.696, 0.606, 0.518, 0.511, 0.477 and 0.420 over days 1–6, against a ramp of 0.38, and was still falling. From that gap I would put the Arctic row at +1 to +4 W/m², but that conversion is a judgment, and a gap that keeps closing would lower it further. The registered +2 to +6 stays the prediction to score.

The 06:42 revision (dacee882) set the Arctic row at +1 to +4 W/m²; the 07:28 table scores the registered +2 to +6 on purpose, because the falsifier (Δ ≤ 0) is the same in both. The outcome will be reported against both ranges. The day-6 φ of 0.420 is newer than the 06:42 revision.

Note, Nov 15 13:30, after the launch and before any row is scored: the Arctic row's quantity and its observation.
- The row scores on the snapshot, as registered: the surface shortwave cloud effect of the day-61, 426 and 791 saves, computed as above. Arm B minus dep64 of that quantity is Δ, and the falsifier Δ ≤ 0 stands as written. dep64's R1, its Arctic row, scores on the window mean as registered; the registration is in dep64's launch brief, which this repository does not hold.
- The ERA5 value behind "−61" is a window mean, and the snapshot cannot be compared with it. The clause "toward ERA5 from below" therefore cannot be scored and is marked unscorable; it stays in the row. What can be compared like for like is the model's window mean, the save's 30-day mean all-sky surface sunlight minus the same offline cloud-free sunlight: ctl64's is −35.3 W/m² at day 61, above ERA5's −60.7. On the window mean a positive Δ moves away from ERA5, unless dep64 already sits below it.
- The verdict will report a like-for-like reading beside the registered one, labeled unregistered: each arm's window mean against ERA5 and CERES, sampled at that arm's own ice70 cells, and |arm − obs| for arm B and dep64. Any observational claim in the verdict rests on that reading, not on the snapshot.
- Sources. ERA5: NCAR RDA ds633.1, ERA5 monthly means, e5.moda.fc.sfc.meanflux; msdwswrf (235.035) minus msdwswrfcs (235.068); 2001–2020, read on a 1° subsample of the 0.25° grid and sampled bilinearly at each model ice70 cell; the months weighted by the window's days (Apr 21–May 20, model days 32–61, at day 61). At ctl64's day-61 ice70 cells, the cells of the mean concentration above 0.15, that is −60.7 W/m². CERES-EBAF Ed4.2.1 surface (rsds minus rsdscs), 2001–2024, the nearest cell: −60.1. `scripts/arcticObs.mjs` fetches both into float32 files; `scripts/arcticSunlight.mjs` with that directory samples them.
- Noise floor for this row: the spread across ctl64's three springs, days 61, 426 and 791. Snapshot (the cloud effect of the snapshot set): −72.6, −58.1 and −67.5 W/m²; the window mean: −35.3, −45.3 and −37.7. The snapshots span 14.5 W/m² (the first two, as the table has them, differ by 14.6), several times the predicted +2 to +6; the window means span 10.0. The snapshot sign test is likely to be dominated by weather.

The committed tools reproduce the measured inputs from ctl64's saves (seaice/ctl64) and the cached observations: `scripts/arcticSunlight.mjs` gives the snapshot −72.64 at day 61 and −58.10 at day 426 (−67.49 at day 791), the window mean −35.31 at day 61 (−45.32 and −37.65 at days 426 and 791), and, by the mean-concentration cells at day 61, ERA5 −60.69 and CERES −60.12; `scripts/cloudEffectBands.mjs` over the 12 saves of days 395–730 prints global SWCRE −48.2 and LWCRE 23.9, 30–60° −51.2 and 23.6, 30S–30N −52.2 and 25.8 (the table's −48.1 is that mean, −48.15, cut toward zero).

Reading rules, Nov 15 13:30, before any scored row of dep64 or arm B has been read. Every falsifier stays as written. Each Δ is also given against the noise floor of its row, labeled as an unregistered reading beside the registered verdict. An extent row with both arms at 0 is reported as not testable.

Early look, Nov 15, before the verdict: a reviewer's probe read the Arctic row's day-61 snapshot of arm B and dep64 (review-armb/work2/tools/runSun1.log, run with arcticSunlight.probe.mjs). Arm B read −98.89 W/m² with the tool as committed and −93.76 with the cloud-ice share placed; dep64 read −82.73. The registered rules and falsifiers are frozen, and none changed after this look. The day-1095 verdict will repeat it.
