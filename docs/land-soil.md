# Land soil: conduction under the skin slab

Oct 27, 2026. Branch `claude/land-soil`. `land.soil` is `'layers'` by default on both engines: three conducting soil layers under the land's skin slab, which is layer 1. `soil: 'slab'` is the skin alone with `heatCapacity` (1e6 J/m²/K); its land update, saved land, WGSL text and PH layout are bit-identical to 5738db9's, the last tree with the slab as the default (`PREVIOUS_LAND` in `test/helpers/previousDefaults.mjs`).

## Options

| option | default | meaning |
|---|---|---|
| `soil` | `'layers'` | `'layers'` the layered soil; `'slab'` the skin slab with `heatCapacity`. With the layers, `heatCapacity` is an error: the skin has `soilHeatCapacity × soilDepths[0]` |
| `soilDepths` | `[0.07, 0.21, 0.72, 1.89]` | layer thicknesses, m (IFS Cy41r2 Part IV ch. 8, Table 8.7); any count of one or more, the WGSL is generated for it |
| `soilHeatCapacity` | `2.19e6` | J/m³/K, constant (IFS §8.5.1) |
| `soilConductivity` | `0.7` | W/m/K, constant: Peters-Lidard et al. (1998) for a dry to moderately wet topsoil |
| `firnDensity` | `400` | kg/m³; ice-sheet cells take (ρC, λ) = (2.05e6 ρ/917, 2.2 (ρ/917)^1.88) = (8.94e5, 0.46) (IFS eqs. 8.29 and 8.44, Table 8.5) |

Per class the layers' capacities are C_k = ρC D_k (soil 1.533e5, 4.599e5, 1.577e6, 4.139e6 J/m²/K; firn 6.26e4, 1.88e5, 6.44e5, 1.69e6) and the conductances between layers are g_k = λ / (½ (D_k + D_k+1)) (soil 5.00, 1.51, 0.536 W/m²/K; firn 3.30, 0.99, 0.35). There is no flux through the bottom, at 2.89 m (IFS eq. 8.60).

## The step

Layer 1 is the skin the atmosphere sees: the radiation column, the exchange and the evaporation are as they are with the slab, and the net surface flux Q enters the skin explicitly. Per land cell and step the layers' changes δ_k solve

    (C_1/Δt + g_1) δ_1 − g_1 δ_2 = Q − F_1
    −g_k−1 δ_k−1 + (C_k/Δt + g_k−1 + g_k) δ_k − g_k δ_k+1 = F_k−1 − F_k
    −g_n−1 δ_n−1 + (C_n/Δt + g_n−1) δ_n = F_n−1,     F_k = g_k (T_k − T_k+1),

by elimination downward and substitution upward. The rows sum to Σ C_k δ_k = Q Δt: the layers' heat content changes by exactly Q Δt, and the conduction is stable at any step. The solve is for the changes from the old temperatures' fluxes, not for the new temperatures, so that float32 rounding stays at the size of the changes. The sublimation heat, the snow melt and the glacier melt use the skin's capacity of the cell's class, so Σ C_k ΔT_k = Q Δt + L_f Δ(snow) − L_f Δ(ice lost).

The explicit coupling is stable while Δt Λ / C_1 < 2 for the surface's sensitivity Λ = 4σT³ + c_p ρC_H U + L β ρC_H U dq_s/dT; with the exchange capped at the lowest layer's mass per step the ratio is at most 1.1 at N = 16 and 0.6 at N = 128 over the August N = 128 probe's land cells.

Ice-sheet cells (`geography.iceSheet`) have the same layers with firn properties: an annual wave has a damping depth of about 2.3 m in firn, so the upper 2.9 m must be firn for the column to store the winter's heat; soil or ice there puts three times the flux through the surface. Snow stays as it was: it has no thermal mass and gives no insulation.

## State

The layers below the skin are the arrays `soilT2`, `soilT3`, … (K, zero on non-land cells): shared buffers on the CPU, slots appended to the PH layout on the GPU (a model with `soil: 'slab'` has the layout without them), fields `land.soilT2` … of a saved state in float32. `regridLand` carries them (the target tile's surface temperature where a cell has no land tile). A state header records `land.soil` and any non-default constants in its `options`, which the worker merges group by group; `scripts/spinup.mjs` and `scripts/splitState.mjs` stamp the soil options the run was given (`LAND`) over `DEFAULT_SOIL` (`scripts/stateOptions.mjs`), so that a state made now starts on the page as it was made when the default moves again. The page's default states record no soil option, so they start on the layers, with the deep layers on the anchored start below. A slab run that reads a layered save ignores the arrays.

### The start of the deep layers

A state saved without the arrays (every state made with the slab) starts them on the annual temperature wave of the column: the periodic solution of the heat equation in a homogeneous soil of depth H, the layers' total (2.89 m), with no flux through its bottom (Carslaw and Jaeger 1959, §2.6; Hillel 1982),

    T(z, t) = T₀ + A F(z, t),     F(z, t) = Re[ e^{iω(t − t_peak)} cosh(k (H − z)) / cosh(k H) ],     k = (1 + i)/d,     d = √(2κ/ω),     κ = λ/ρC,

taken at each layer's mid-depth z_k, with ω = 2π per year, κ of the cell's class (soil or firn: d = 1.79 m and 2.29 m) and t the model time of the state (day 0 is the March equinox). A is the amplitude of `airCycle(lat, elevation)` (`js/physics/land.module.js`). T₀ puts the wave's surface value on the state's skin temperature: T_k = T_skin + A [F(z_k, t) − F(0, t)], with F(0, t) = cos(ω (t − t_peak)). Layer 1 keeps the saved skin, and the mean of the wave is the skin's, whatever `airCycle` makes of the annual mean there. A fresh start (`initialize`), and a cell whose skin is missing or not above zero, have no skin to anchor to and take T₀ = `airCycle`'s mean. A save that carries the arrays uses them as saved. The skin is the instantaneous saved value, so its diurnal phase enters the start: the 1e6 slab's skin ranges 5 to 7 K in a day in August at 10 to 60N.

`t_peak` is the surface temperature's seasonal maximum: the June solstice (day 91.25 of `insolationCycle`'s year) in the northern hemisphere and the December solstice in the southern, plus `SOIL_LAG`. Nothing in `insolationCycle` or `airCycle` records a lag (the fit behind `airCycle` kept the mean and the amplitude, not the phase), and no saved field of the model gives one, so `SOIL_LAG` is 30 days, the usual three to five weeks by which the observed continental maximum of surface air temperature follows the solstice. Ten days of error in it moves a layer's start by at most 0.17 A.

What the start leaves: the column's four layers are thick against the damping depth (layer 4 is 1.89 m, d is 1.79 m), so the stepped column's own periodic state differs from F at the layers' mid-depths by 0.03 to 0.05 of the surface amplitude. Forced at its top by the wave for a year (`test/landSoil.test.mjs`), a column started this way keeps its layer means within 0.21 K (north), 0.13 K (south) and 0.13 K (firn) of the mean, and its deepest layer ends the year 0.55, 0.36 and 0.53 K from where it began; a start at the skin temperature gives 0.90, 0.58 and 1.15 K and 2.41, 1.56 and 4.82 K. Started on the stepped column's own periodic state (the continued fraction of its layers' admittances) the same year drifts 0.014, 0.009 and 0.005 K in the means and 0.08, 0.05 and 0.05 K at its end.

Layers 2, 3 and 4 minus the saved skin at the start, on the N = 128 state of Aug 19 (`s6klm/arm-s0_day0882`, area-weighted means; the departures are the wave's own, A [F(z_k) − F(0)], since the skin cancels):

| region | layer 2 | layer 3 | layer 4 |
|---|---|---|---|
| Antarctica | +0.3 K | +1.5 K | +4.8 K |
| Greenland | −0.3 K | −1.5 K | −4.8 K |
| India, 10–28N 70–90E | −0.3 K | −1.3 K | −3.9 K |
| 40–60N land | −0.6 K | −2.6 K | −8.0 K |

`scripts/savedLand.mjs` loads a save's land onto a model for the spin-up, the fixed-SST run, the ocean spin-up and the figure scripts: the model's own sea ice keeps the snow the save carries on it, and the deep layers start from the state's skin at its time.

## The GPU

`js/gpu/soil.gpu.js` generates the constants of both classes and the step as straight-line WGSL on scalars (`colT2`, `colF1`, `colB2`, …) with no arrays, loops, pointers or module-scope private storage, which avoids the Metal compiler failures listed in `CLAUDE.md` by construction. `physicsSurfaceKernel(held, seaSkin, soil)` takes the step from its third argument; with `soil` null the kernel text is the slab's. The layers take `soilDeep` PH slots (`layoutFor`), travel through `uploadPhysics` and `uploadLand`, and come back with `pullLand` and `land.serialize()`.

## Tests

- `test/landSoil.test.mjs` (CPU): the column's constants and errors; the land energy row (soil and firn columns under a constant flux, a flux that changes sign with steps of one to five times 1350 s, melting snow with sublimation, and glacier melt); a uniform column unchanged without a flux; the stepped column against the backward-step solution of its own equations, step by step over a period, at 12 hours to a year; a one-layer column against the slab; the start against the formula (the mean over a year is the skin and the amplitude by depth is A |c_k − 1|, in both hemispheres and in firn; airCycle's mean without a skin; a missing, zero or NaN skin; saved arrays used as saved; a cell with an empty layer); the drift over a year; the slab option's digest over 200 steps and its serialized keys, and the layers as the default; the state-file round trip and `regridLand`.
- `test/landLoadPath.test.mjs` (GPU file): the scripts' land load from a save with snow on sea ice. The CPU land after the load and a CPU model one step later reproduce the digests 5738db9 gave; the GPU model keeps the sea-ice snow with the layers and with the slab, where a null ice argument drops it; and a slab GPU model loaded that way steps once to the same state, bit for bit, as one loaded with the call 5738db9 made.
- `test/landSoilGpu.test.mjs` (GPU file): eight coupled steps at N = 16 against the CPU model beside a slab of the skin's capacity; one step with the GPU's own fluxes fed to the CPU land; the energy row on the GPU's pulled fields; the WGSL for one, two and five layers against the CPU land fed with the GPU's fluxes; and the slab's kernel text, constants and PH layout.
- The suite with the layers as the default passes (883 CPU tests, 407 in 54 GPU files). The files whose expectations are the slab's name it with `PREVIOUS_LAND`: `land`, `landIce`, `snowAlbedo`, `frameGpu` and the load digests of `landLoadPath`; `docs/cloud-gas-overlap.md` lists them with the other defaults.

## Measured

- Energy row: the CPU column's heat content changes by Q Δt + L_f Δ(snow) − L_f Δ(ice lost) to 1e-16 of the column's heat content at every step in all four cases; the GPU's pulled fields, with the net flux of the step, to 2.5e-8 (float32), over 1193 land cells of which 141 are on ice sheets.
- The CPU column follows the backward-step solution of its own equations to 1e-10 K at periods from 12 hours to a year in soil and firn; the skin's response to a sinusoidal flux is 0.98 (12 hours) and 0.9998 (10 days) of the continuous column's.
- Engines at N = 16, eight coupled steps: the skin differs by 5.1e-4 K (2.9e-4 K for a slab of the skin's capacity), the layers by 7.6e-5, 9.2e-5 and 1.3e-4 K at every land cell; one step with the GPU's own fluxes fed to the CPU land: 3.6e-5 K in the skin and 1.5e-5 K in the layers; the GPU energy row 2.4e-8 of the column's heat content over 1199 land cells.
- Engines at N = 128 on the page's default state (`fourteen128_day1825`, 60 layers, the options its header records; `node scripts/compareEngines.mjs` at 1, 4 and 16 steps): the skin differs at all cells by a relative rms of 4e-7, 2e-6 and 6e-6 (at most 0.05, 0.18 and 0.45 K), soilT2 by at most 1e-4, 9e-4 and 8.6e-3 K, soilT3 and soilT4 by at most 3e-5, 8e-5 and 2.6e-4 K. The tool's other regions are over its 1e-3 tolerance in the same 24 families as with `soil: 'slab'` (31 failures against 34); the lighter skin lets a flux difference between the engines move the surface about four times as far as the slab's 1e6 J/m²/K does (skin rms 5.7e-6 against 1.4e-6 after 16 steps), so OLR separates earlier (1.4e-3 against 2.6e-4 after 4 steps).
- Pace, N = 128 on the M5, 512 steps queued 8 to a submission, slab and layers alternating over three rounds: 39.4 and 39.3 s per model day in the first, 39.4 and 39.3 in the second, 42.9 and 41.6 in the third (a busier GPU); no cost can be told from the noise, under 0.3 percent in the quiet rounds.

## Not included

Snow insulation (IFS §8.4.3), soil freezing (§8.5.2), a conductivity and capacity that follow the soil's water (`soilConductivity: 'wetness'` from the bucket), the canopy's skin coupling (IFS Λ_sk) and a graded firn-to-ice column deeper than 2.9 m are separate changes.
