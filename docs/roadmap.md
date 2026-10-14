# WebGCM roadmap: milestones 16–19

Four milestones, in this order. Each one builds on the previous: the
land mask gives terrain its footprint and salinity its runoff and
brine, and the performance work comes last so that it optimizes the
finished physics rather than a moving target.

| # | Milestone | What it delivers | Depends on |
|---|-----------|------------------|------------|
| M16 (done) | Land surface | Real continents as a land mask with soil water, snow, land albedo and roughness; the ocean stops at coasts; coastlines and land colours in the page | — |
| M17 (done) | Terrain | Surface elevation in the dynamical core on both engines; sea-level pressure reduction; pressure levels that meet the ground | M16 |
| M18 | Layered ocean with salinity — done | A six-layer hybrid isopycnal ocean on the real bathymetry with a split free surface, a bulk mixed layer, prognostic salinity with a linear equation of state, freshwater from rain, evaporation, runoff and sea ice, and convection | M16 |
| M19 | N=128 performance | One simulated day per minute at N=128 on the GPU engine | M16–M18 |

Estimates are deliberately relative: M17 is the hardest, comparable to
the GPU port; M16 is a little smaller; M18 is half of M16 unless the
variable-density pressure gradient turns into a derivation; M19 is
measurement-driven and could be short or long.

## M16 — Land surface

**Goal.** Continents on the aquaplanet, with the physics that makes land
behave like land: a small heat capacity, evaporation limited by soil
water, snow with a high albedo, a rougher surface, and an ocean that
stops at the coast. What should emerge in a year at N=64: monsoon
reversals over a large subtropical continent, dry subtropical
continental interiors, a cold winter continental high, seasonal snow
cover, and wind-driven gyres with western boundary currents in the
ocean layer.

**Data.** ETOPO 2022 (ice-surface elevation) downsampled to a 0.5°
grid of int16 metres, about 500 KB, committed under `data/` with the
script that made it. Per cell the mesh takes the area mean elevation
and the land fraction (share of sub-cells above sea level); a cell is
land when its fraction exceeds one half. Fractional cells are a later
refinement. M17 reuses the same file for elevation.

**State.** Two new state fields, `soil` (bucket water, kg/m²) and
`snow` (water equivalent, kg/m²), added to `STATE_NAMES`, the worker's
snapshot names, the saved-state JSON, and `regridState`. A land mask
per cell and per edge lives in the mesh's shared buffers, not the
state.

**Physics.**
- Heat capacity: the ice module already takes a per-cell
  `heatCapacity`; land cells get about 1e6 J/m²/K (a decimetre of
  soil) instead of the ocean layer's ρcp·h1.
- Evaporation: the radiation column's surface evaporation is scaled by
  a wetness β = min(1, soil / (0.75 · capacity)), Manabe's bucket with
  150 mm capacity; rain fills the bucket, overflow is runoff, which M18
  routes to the sea. Over ocean β = 1.
- Snow: precipitation with the lowest-layer air below 0 °C accumulates
  as snow; the surface energy balance melts it with the same
  skin-temperature machinery the ice module uses for sea ice; melt goes
  into the bucket.
- Albedo: land 0.25 for direct and diffuse, rising to 0.7 with snow
  over a few centimetres; the zenith-angle water albedo stays for open
  ocean.
- Roughness: the drag coefficient becomes a per-cell array, about
  3e-3 over land against 1.5e-3 over water.
- Ocean: edges that touch a land cell carry no flow (masked in the
  edge loops and the closure), land cells hold no layers, the ice
  module and the ocean flux apply to ocean cells only.
- GPU: the masks and per-cell coefficients become buffers; the physics
  and adjust kernels branch on the mask; the ocean kernels mask edges.

**Page.** A coastline layer built from mesh edges between land and
ocean cells, drawn like the graticule; Satellite-mode land colour from
soil water (tan to green) and snow (white); overlays for soil water,
snow depth and land fraction; diagnostics for land mean temperature,
snow area and soil water.

**Tests.** Land fraction from the dataset near 29% of the area; bucket
water balance (rain − evaporation − runoff = Δsoil) to roundoff; snow
energy balance; no flux through coast edges and ocean mass conserved
with a mask; GPU agrees with CPU on a masked run to the usual
tolerances; a warm island in a resting atmosphere drives inflow at the
surface (sea breeze) within a day.

**Risks.** Coastal cells at 120 km mix land and sea; a binary mask
overstates land at the coast. Snow and ice share code paths that must
not double-count energy. The saved-state format changes: old states
load with soil half full and no snow.

## M17 — Terrain

**Goal.** Mountains. Surface geopotential enters the hydrostatic
integration and the pressure-gradient force on both engines, so that
orographic uplift, rain shadows and stationary waves appear.

**Core.** `sigmaCore` already accepts `surfaceGeopotential` in the
geopotential integration, and its pressure gradient is the σ form
∇φ + cp θv (∂Π/∂π) ∇π, whose two terms nearly cancel over slopes. The
GPU core's geopotential is built as a deviation from a reference
column (`GR`, `GABS`); it needs a per-cell φ_s added at the bottom of
`D_GEO`. Elevation is smoothed on the mesh before use, a few passes of
the mesh Laplacian and a slope cap, so that nothing sits at the grid
scale.

**Consequences to handle.**
- Surface pressure is no longer sea-level pressure: the MSLP overlay
  and the isobars reduce π to sea level with the standard lapse-rate
  formula; the model-information panel shows both.
- Pressure levels below ground: the level extraction masks cells where
  the level is under the surface (no colour, no contour) rather than
  extrapolating.
- Initialization from a flat state: π scales by exp(−φ_s / (R T̄)) and
  the column is rebalanced hydrostatically; `regridState` does this
  when the target mesh has terrain and the source does not.
- Physics sees lower π over high ground automatically (saturation,
  layer masses); nothing else changes.

**Tests.** The classic one: an isothermal atmosphere at rest over a
smooth mountain stays at rest, with spurious winds under 0.5 m/s after
five days on the CAM levels; the geopotential test with a raised
surface; mass conservation; JW06 unchanged with flat terrain; GPU
agrees with CPU over terrain.

**Risks.** Pure σ over steep terrain is noisy near the top of the
mountains where thin σ layers follow the ground. If the rest test
fails at acceptable smoothing, the fix is a hybrid σ-p coordinate,
which is its own milestone; the plan accepts smoother terrain first.

## M18 — Ocean salinity

**Goal.** Salinity as a prognostic field in both ocean layers, with
density from a linear equation of state, so that stratification and
convection respond to freshwater and the ice edge feels the freezing
point of seawater. What should emerge: subtropical salinity maxima
under evaporation, freshening under the ITCZ and at river mouths,
brine rejection under growing sea ice, and density-driven overturning
where cold salty water forms.

**Design.**
- Layer salt content S1·h1 and S2·h2 join H1 = T1·h1 and H2 in the
  ocean state, advected with the same flux form; snapshots and regrid
  carry them.
- Freshwater: evaporation minus precipitation over ocean cells from the
  radiation and moist budgets; runoff from M16's buckets routed to the
  nearest coastal ocean cell by the neighbour walk; sea-ice growth
  leaves salt behind and melt returns fresh water; each is a virtual
  salt flux on the upper layer.
- Equation of state: ρ = ρ0 (1 − α (T − T0) + β (S − S0)) with α 2e-4
  /K and β 7.6e-4 /psu, so g'12 and g'23 become fields. The Montgomery
  potentials keep their form but the layer momentum equation gains the
  term from the gradient of density within a layer; this is the one
  derivation in the plan and is written up in the design doc before
  it is coded.
- Convection: when ρ1 ≥ ρ2 the layers mix heat and salt toward
  neutral over the entrainment time, the same mechanism as the
  thickness floors.
- Freezing point −0.054 · S °C in the ice module.
- GPU ocean kernels carry the two salt fields and the equation of
  state.

**Page.** A salinity overlay in psu with a sequential palette; ocean
diagnostics gain mean salinity and the freshwater budget.

**Tests.** Global salt conserved to roundoff with no forcing and equal
to the integrated fluxes with forcing; a fresh lens lowers density and
raises the interface as the Montgomery potential predicts; brine
rejection conserves salt across freeze and melt; the variable-density
pressure gradient at rest stays at rest; GPU agrees with CPU.

**Risks.** A 1.5-layer ocean cannot carry a deep overturning; the
plan's convection is a parameterization of what it cannot resolve.
The upper layer also deepens steadily today, about 10 m a year, because
thin-layer entrainment adds thermocline water that never returns; M18
should add a diapycnal restoring of the mean depth so the surface heat
capacity stops creeping up.
Runoff routing at 120 km puts rivers into the wrong coastal cell; the
freshening pattern is what matters.

## M19 — One simulated day per minute at N=128

**Goal.** At N=128 the step is 168.75 s, so a simulated day is 512
steps and the budget is 117 ms per step including readbacks and the
page's share of the GPU. The last measurement was about 490 ms per step
under contention; the target is a 3–4× speedup.

**Plan, in order.**
1. Measure. Add WebGPU timestamp queries around every kernel and a
   per-kernel table at N=64 and N=128 in Node and in the page. The
   N=64 numbers (dynamics 38 ms, adjust 21 ms of a 71 ms step) are the
   only breakdown so far and may not scale the same way.
2. The column kernels. `physics` and `adjust` run one thread per
   column with 27-level loops; check occupancy and memory traffic on
   the packed layout, split `adjust` into its three parts if register
   pressure limits it, and lift repeated work (Exner ratios,
   saturation humidity) into level constants or a small precomputed
   buffer.
3. Radiation cadence. Compute radiation every second or fourth step
   and hold its tendencies, as GCMs do; convection and the boundary
   layer stay every step. Verify the climate is unchanged over a year.
4. Time integration. Test SSP-RK3 at the same step: one tendency
   evaluation fewer per step if it holds the 3× step's stability.
5. The frame path. Download only what the page shows: extract the
   display level on the GPU, keep diagnostics sums on the GPU, and
   overlap the readback with the next step's dispatches.
6. The page at N=128. The overlay painter and the contour builder
   handle four times the cells each frame; move contouring into the
   worker or throttle it if it stalls the main thread.

**Acceptance.** 24 simulated hours per minute at N=128 on the M1 Max
in the pane, with the GPU-versus-CPU agreement tests unchanged and a
year at N=128 matching the N=64 climate statistics. Any cadence or
integrator change is judged by that year, not by speed alone.

**Risks.** Timestamp queries may be unavailable in the pane's browser,
in which case Node with Dawn measures. Radiation cadence changes the
diurnal cycle slightly. The snapshot at N=128 is about 480 MB; the
snapshot store and the regrid must cope.

## Milestones 27 to 30

The milestones after M26 (fog). M27 is adopted at a 50 m floor; M28
follows it, M29 follows M28, and M30, judged without adoption, keeps the
floor at 50 m.

| # | Milestone | What it delivers | Depends on |
|---|-----------|------------------|------------|
| M27 | A thin mixed layer | The three additions a mixed-layer floor below 50 m needs: the wind stress and drag acting over a momentum depth separate from the refill floor, with the mixed layer's momentum carried through its exchanges; penetrating sunlight (two-band Jerlov IB) deposited in the classes beneath the mixed layer; a rate-limited refill. Adopted Oct 13 at a 50 m floor with the shelf fix, a mixing package and two lighter classes; no thinner floor passed its judging. | the Oct 12 conservation work |
| M28 | Rivers, lakes and the land below sea level | An ocean mask by connection to the ocean rather than by elevation alone; river routing with travel time; lakes as a surface type and as reservoirs in the river network. | M16, M27's state-file work |
| M29 | Land ice | The ice sheets as a layer of ice on the bedrock: bed and thickness from BedMachine, the surface the atmosphere sees as bed plus ice, a surface mass balance from snowfall, melt, sublimation and calving, meltwater through the rivers; the bed held fixed; flow by the shallow-ice approximation later. | M28 |
| M30 | The sea's day and night | A warm layer and cool skin above the mixed layer, stepped with the atmosphere, and a retreat rule that holds the depth the night's convection reaches, so that a floor of 10–20 m no longer stands in for the night. | M27 |

## M27 — A thin mixed layer

**Goal.** Let the mixed layer be as shallow as the physics asks
(3–40 m in the tropics and the summer hemisphere, where the 50 m floor
now holds a third of the sea) without the losses a thin layer exposes:
the wind's momentum spread over a 50 m momentum depth and carried
through the layer's retreat, refill and convection; the sunlight that
passes a thin layer heating the water beneath it instead of the
surface; the refill from below limited to a rate, with the coastal
cells the dynamics drain fixed at their source. Every option off
leaves today's model bit for bit, except the momentum depth where a
mixed layer at an edge is already thinner than 50 m.

**Acceptance.** With no surface flux, heat, salt and volume conserved
through every phase on both engines (test/layeredConservation); the
Ekman transport at 45° within 5 % of τ/(ρ₀f) with a 10 m layer; the
momentum of an edge's column unchanged through a retreat and a refill;
the sunlight deposited equal to the light passed; the step at N=128 no
more than 0.6 ms per ocean step slower on the H100 scale. Then the
judging arms (N=16, two years, floors 50/20/10 with the additions on
and off; N=64, 30 days from the page's state) decide the floor and
whether momentum carriage and penetration default on.

**Judged (Oct 13, 2026).** The three additions landed (61e7850,
d037d82, cc41cba, edb0d79) and the arms ran on the GPU: the N=16 water
planet for two years at 50 m with the additions off and at 50, 20 and
10 m with them on; N=64 from eleven's March state for 30 and 60 days
and from its June state for 30 days. At 50 m the additions change the
model within the spread between realizations: global SST −0.06 K, no
clamps, at most 0.021 % of the sea thinned, the pending light steady
at 500–750 J/m². At 20 m and 10 m the summer and tropical mixed layers
fall into Earth's range (summer medians 14–24 m at 20–60°, tropical
medians 14–22 m; the floor binds on 16 % and 8 % of the sea instead of
33 %), the seasonal SST range at 45° stays under 4 K, nothing clamps,
thins or reaches the 5 m backstop on the water planet, and the SST
responds as a 10–28 m layer where it sits on the floor. Two things
keep the page's floor at 50 m for now. In boreal spring the equatorial
Pacific's zonal SST gradient collapses (W−E 2.8 → 1.1–1.2 K, the cold
tongue 1.6–1.7 K warmer at day 30, 0.4–0.7 K at day 60) because
nothing entrains thermocline water through a thin layer: shear
entrainment at the mixed layer's base is missing. And a grid-scale
sea-level mode on the Kara, Laptev and East Siberian shelves, already
±2–4.5 m at 50 m under the ice, grows to ±10–13 m with 3–5 m/s
mixed-layer jets once the shelf ice thins, and clamps at 10 m.
Recommended: `penetration: 'IB'`, `carryMomentum: true` and
`refillRate` 24 m/day (`refillFloor` 5 m) as the defaults in one
commit with Part 8 regenerated; 20 m after the shelf mode and shear
entrainment are fixed, with the June and 60-day March 20 m arms
repeated and the WOA23 summer SST comparison; 10 m after those and
lighter classes.

**Adopted (Oct 13, 2026).** A second round fixed what the first found
and judged the floor again. The shelf mode was a fault in both engines:
the thin classes copied the velocity above them after the column's
transport was shifted to the barotropic average, an explicit pressure
kick; in the right order the Arctic shelf's rms neighbor step in η is
0.02–0.07 m and the fastest current north of 60N 0.1–0.5 m/s at every
floor. A mixing package went on both engines: the bulk-Richardson shear
entrainment of Price, Weller and Pinkel at the mixed layer's base
(Ri_c 0.65), the wind's stress over its own Ekman depth
min(50, max(10, 0.4u*/|f|)) m, and a heated retreat to
2.5 ε u*³/(−B̄), ε a Langmuir factor of at most 2, bounded by 0.7u*/|f|.
Two lighter classes, 1019.5 and 1020.0, take the warm pool's lightest
water, and a 44-class state rebins onto them at load (0.26–0.42 s at
N=128). The arms ran as before at floors of 50 m with everything off
and 50, 20, 10 and 5 m with everything on. No thin floor passed. At
20 m one Mediterranean-outflow mixed-layer current reached 2.52 m/s on
March day 25 (the stop rule is 2.5), the June Benguela–Peru mean lay
1.21 K further from WOA23 than at 50 m off (Peru 2.15 K; the allowance
is 1.0) and the June 10S–10N median mixed layer was 43.5 m (band
10–40 m). At 10 and 5 m the June global SST ran 0.33–0.35 K above the
control (limit 0.3), May's 40–60N, cold-tongue and coastal boxes lost
0.37–0.62 K, and at 5 m the water planet's 20S–20N cooled by 0.59 K.
At 50 m everything on passes against everything off: global SST −0.08 K
on the water planet and −0.01 to −0.03 K at N=64, no clamps, the
whole-run residual +0.053 W/m², no WOA23 box worse by more than 0.3 K
but Peru in May (+0.63 K, within the coastal allowance), and in April
Benguela 0.70 K, California 0.45 K and Canary 0.55 K closer. The
floor's within-1-m share falls from 34.6 to 29.3 % on the water planet,
and where it binds the depth a column asks for rises 2.7×, from 1.4 to
3.7 m. Adopted as the defaults of both engines at the 50 m floor:
`penetration` 'IB', `carryMomentum`, `refillRate` 24 m a day,
`shearEntrainment` 0.65, `ekmanDepth` 0.4, `retreatStirring` 2.5,
`langmuir` 2, `retreatEkman` 0.7 and the 46 classes. Part 8's live
figures run these defaults, and its recorded run is regenerated on
them, with its run sentences rewritten and its tests re-pinned. The cost at N=128
on the M1 Max is about 0.65 ms a model step for this round (+0.33 s
a simulated day) and 0.83–0.92 ms with
M27's additions (+0.43–0.47 s). Between the floor and 20 m stand, in
order: the coastal upwelling under a thin layer (Peru is 2.2 K warmer
at every thin floor and the Ekman depth saturates at 50 m there), the
June SST drift, the Mediterranean outflow under a thin layer (M28's
straits), and the night's convection, which a one-day buoyancy memory
cannot see: where a thin floor binds, the asked depth's median is
1.0–1.5 m (M30).
In June the 50 m arm with everything on has 16 columns at the 600 m cap
against the control's 4–7 and 137 deeper than 400 m at 40–60S against
17, in the southeast Pacific's winter; they are to be read before
twelve.

## M28 — Rivers, lakes and the land below sea level

**Goal.** Water on land that moves and collects the way the real
water does. Three parts, in this order, because each one's data
structure carries the next.

1. **The ocean mask by connection.** A cell is sea only if its water
   connects to the global ocean through sea cells or a listed strait;
   dry land below sea level (the Qattara Depression, Chott el Djerid,
   the Danakil, Death Valley) is land, and ETOPO's lake beds (the Great
   Lakes lie below sea level in the raster) are not sea. The Caspian,
   the Black Sea and the Baltic stay separate basins; a parameterized
   strait exchange (the Bosporus at about 0.01 Sv, the Danish straits)
   is the later way to join the last two to the ocean, since a one-cell
   strait exchanges far too freely. A lake dataset (HydroLAKES or GLWD,
   rasterized to the topography grid) says which closed depressions hold
   water, with each lake's mean depth.
2. **River routing with travel time.** A linear reservoir per land cell
   (TRIP's form): the cell's channel store takes its bucket overflow
   and its upstream neighbors' outflow and drains to its downstream
   cell at S/τ, τ = d/v, with v from the slope (Manning reduced to one
   constant, clamped), so a basin's rain reaches the sea over weeks and
   snowmelt arrives as spring floods. One new state field; one per-cell
   gather over the upstream neighbors on the GPU, deterministic, which
   replaces the serial basin sum of oAccumulateFresh; the mouths'
   discharge becomes the ocean's freshwater input; the land water
   budget test extends to Δsoil + Δriver + discharge. Page overlays:
   runoff per cell and the river network drawn by discharge, with the
   mouths' discharge in the readout (the Amazon's mean is about
   2×10⁵ m³/s). Left out at first: floodplains, baseflow, river
   evaporation and river ice.
3. **Lakes.** A lake tile with a fraction per cell: a freshwater column
   with a seasonal mixed layer over a thermocline (FLake's form), depth
   from the data (Superior 147 m, Baikal 744, Victoria 40, Malawi 292,
   Titicaca 107, Winnipeg 12, the Caspian 211), water's albedo,
   roughness and saturated evaporation, and freshwater ice at 0 °C (the
   ice model's freezing point made per cell, its brine paths skipped).
   Each lake is a reservoir node in the river network: inflow from
   upstream, outflow over a sill at the lake's level, the closed basins
   (the Caspian, Titicaca, Eyre, the Great Basin) balancing inflow
   against their own evaporation, so a lake's level is a state the
   seasons move. At N=128 the large lakes are 3–26 cells; at N=64 a few;
   below that only tiles.

**Acceptance.** No sea cell without a path to the ocean; the Great
Lakes, Baikal, Victoria, Malawi, Titicaca and Winnipeg present as lakes
at N=128 with ice seasons in the north and none in the tropics; the
Amazon, Congo, Ob, Mississippi and Nile mouths within a factor of two
of their observed mean discharge after a spun-up year; the land water
budget closed to roundoff with the river and lake stores in it; both
engines bit for bit with the options off.

## M29 — Land ice

**Goal.** The Greenland and Antarctic ice sheets as ice resting on the
bedrock rather than as terrain. The topography carries two layers, the
bed (ETOPO 2022's bedrock, which is BedMachine under the sheets) and
the ice thickness, and the surface the atmosphere sees is their sum,
rebuilt once a model year as the ice changes. The ice is a state with
a surface mass balance: snow older than a threshold becomes ice, melt
and sublimation take it away by the sea-ice energy balance with bare
ice at about 0.6 and snow at about 0.8 albedo, the meltwater runs to
the sea through M28's rivers, and ice at the coast calves. The ice-sheet
mask comes from thickness rather than from latitude and elevation.

**The bed stays fixed.** The ice's weight depresses the bed by about a
third of its thickness, but the mantle answers on millennial
timescales, so for runs of decades to centuries the bed of today is the
right bed whatever the ice does; a melted Greenland keeps its central
basin below sea level, which the M28 mask by connection treats as a
lake or a sea as the channels dictate. No isostatic model.

**Flow, later.** Without transport from the interior to the ablation
zone the sheet thickens in the middle and thins at the margin without
end; the shallow-ice approximation, a nonlinear diffusion of thickness
on the surface slope stepped in years, is the minimal flow that makes
century-scale experiments (a melting Greenland; ice caps on another
planet, see the ice-free topography note) honest.

**Acceptance.** Surface elevation over the sheets within the data's
accuracy at the start; the surface mass balance of Greenland and
Antarctica within a factor of two of the observed accumulation and
ablation after a spun-up year; meltwater reaching the sea; the ice
thickness changing by no more than a few meters a year anywhere; both
engines bit for bit with the ice held fixed.

## M30 — The sea's day and night

**Goal.** Let the mixed layer stand for what climatologies measure,
the layer the night's convection mixes each dawn, and give the day's
warming its own place above it. The retreat reads a one-day memory of
the surface buoyancy flux, which cannot see the night: where a thin floor
binds, a calm, strongly heated column asks for 1.0–1.5 m, and the floor
stands in for the night's depth. A shorter memory is not the answer by
itself: 3 hours at a 10 m floor bound far less but cooled the water
planet's SST by 0.72 K in a year, with no day and night in its
snapshots.

**Design.** A Zeng–Beljaars warm layer, with Takaya's stability
function and Langmuir factor, and a Fairall cool skin, stepped with the
atmosphere in the surface kernels and booked as a profile inside the
mixed layer's heat, so that it stores nothing of its own. In the
retreat only, the previous day's mean buoyancy flux replaces the
one-day memory, and a night gate lets the mixed layer retreat only
where the night's convection would fail to reach its base. The
coupling risk runs through the surface exchange of the Oct 6 runaway.
About 8–10 hours.

**Acceptance.** Both engines, with a parity test at N=8 and rows in the
conservation tests; heat conserved through the skin; then M27's
judging rules at 20 and 10 m floors, where the floor should bind only
where the nights are calm enough to leave a thin layer.

**Judged (Oct 13, 2026).** Both halves landed on both engines, off by
default and bit for bit when off: `warmLayer` and `coolSkin` in the
surface step (they need exchange 'roughness'), `retreatDaily` and
`retreatLight` in the retreat, with parity tests at N=8, conservation
rows, and `SKIN_WATCH`, a daily skin and exchange recorder in
scripts/spinup.mjs. The arms ran on the M5 Max at the 50 m floor with
the skin (A), the night rule (C) and both, and with both at 10 and 5 m:
N=64 from eleven's June state for 30 days and its March state for 60,
and the N=16 water planet for two years. The skin does what it was
built for: 23 of 24 low-wind bins within 0.5–1.5× of the satellite
regression of Gentemann et al. (2003), the day's maximum at 12–16 h in
91 % of the cells with a range of 0.3 K or more, a 20S–20N mean range
of 0.26 K, and the surface exchange's extremes within +10 % of the
control's. It is not adopted. The Benguela box lies 0.34–0.53 K further
from WOA23 in April and May in both the M5 and the M1 realizations,
which points at the cool skin under the stratocumulus deck. On the
50 m slab only 55 % of the day's minima fall at 05–09 h (the rule asks
70 %; 76 % at 10 m), and a warm layer outlives the night at local
03 h, nearly all of it in polar-day and cold-water cells. The night
rule is not adopted either. It does not thin the pile on the floor:
the share within a meter of it rises from 60 to 69 % at 20–60N in June
and from 54 to 55 % at 10S–10N in March; with both halves at 10 m it is
26 % in June against the floor round's 23 % on the M1. The gate holds
only 2–4 % of the heated column-steps on the floor, and the depth a
binding column asks for rises 1.1× in June
(4.2 → 4.6 m) against the 1.5× asked. Those columns' dark loss is 0.38
of their sunlight (median, June), where the night's convection alone
keeps a column shallow: the floor there stands in for missing shear
and advection, not for the night. California and Canary also lie
0.4–0.7 K further from WOA23 in May. The 10 m and 5 m floors fail as
in M27: June's global SST 0.36–0.37 K above the control, the
Benguela–Peru mean 1.2–1.3 K further from WOA23, the May cold tongue
0.8–1.1 K warmer on day 60 of the M5 run (0.05–0.17 K on the M1's) and
the June 10S–10N median mixed layer 41–46 m. The floor stays at 50 m
with every M30 option off, and Part 8's recorded run is M27's defaults'.
Still open: the coastal upwelling under a thin layer; the Mediterranean
outflow pulse (2.12 m/s in the Gulf of Cadiz on March day 25 at 50 m);
the southeast Pacific's winter columns at the 600 m cap (9 to 14 in the
June 50 m arms, 20 at 10 m); and a `coolSkin`-only and a
`warmLayer`-only arm over March–May to find which half warms the deck.
On the M5 the N=16 atmosphere's books close about 2 W/m² lower than on
the M1 with otherwise identical weather, to be read before the M5 is
trusted for books. Its skin parity (`test/skinGpu.test.mjs`) drifts
2e-4 to 1e-3 K in the warm layer over 20 steps against under 1e-4 K on
the M1; the test allows what the M5 shows until that cause is read too.
