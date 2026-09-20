# Roadmap: milestones 16–19

Four milestones, in this order. Each one builds on the previous: the
land mask gives terrain its footprint and salinity its runoff and
brine, and the performance work comes last so that it optimizes the
finished physics rather than a moving target.

| # | Milestone | What it delivers | Depends on |
|---|-----------|------------------|------------|
| M16 | Land surface | Real continents as a land mask with soil water, snow, land albedo and roughness; the ocean stops at coasts; coastlines and land colours in the page | — |
| M17 | Terrain | Surface elevation in the dynamical core on both engines; sea-level pressure reduction; pressure levels that meet the ground | M16 |
| M18 | Ocean salinity | Prognostic salinity in both ocean layers with a linear equation of state; freshwater from rain, evaporation, runoff and sea ice; density-driven convection | M16 |
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
