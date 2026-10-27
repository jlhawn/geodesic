# Land soil: conduction under the skin slab

Oct 26, 2026. Branch `claude/land-soil`. `land: { soil: 'layers' }` puts three conducting soil layers under the land's skin slab on both engines. The default `'slab'` is the skin alone with `heatCapacity` (1e6 J/m²/K) and runs the code it always has: the land update, the saved land, the WGSL text and the PH layout are bit-identical to what they were.

## Options

| option | default | meaning |
|---|---|---|
| `soil` | `'slab'` | `'slab'` the skin slab with `heatCapacity`; `'layers'` the layered soil. With `'layers'`, `heatCapacity` is an error: the skin has `soilHeatCapacity × soilDepths[0]` |
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

The layers below the skin are the arrays `soilT2`, `soilT3`, … (K, zero on non-land cells): shared buffers on the CPU, slots appended to the PH layout on the GPU (the default layout is unchanged), fields `land.soilT2` … of a saved state in float32. `regridLand` carries them (the target tile's surface temperature where a cell has no land tile). A state header records `land.soil` and any non-default constants in its `options`, which the worker merges group by group; `scripts/stateOptions.mjs` pins `land: { soil: 'slab' }` as a default so that states made before the default moves keep starting as they were made. A slab run that reads a layered save ignores the arrays.

### The start of the deep layers

A state saved without the arrays (every state made with the slab) starts them on the annual temperature wave of a homogeneous half-space, the periodic solution of the heat equation under a sinusoidal surface temperature (Carslaw and Jaeger 1959; Hillel 1982):

    T_k = T̄ + A exp(−z_k/d) cos(ω (t − t_peak) − z_k/d),     d = √(2κ/ω),     κ = λ/ρC

with z_k the layer's mid-depth, ω = 2π per year, κ the cell's class (soil or firn: d = 1.79 m and 2.29 m), T̄ and A the mean and amplitude of `airCycle(lat, elevation)` (`js/physics/land.module.js`) and t the model time of the state (day 0 is the March equinox). Layer 1 keeps the saved skin. `t_peak` is the surface temperature's seasonal maximum: the June solstice (day 91.25 of `insolationCycle`'s year) in the northern hemisphere and the December solstice in the southern, plus `SOIL_LAG`. Nothing in `insolationCycle` or `airCycle` records a lag (the fit behind `airCycle` kept the mean and the amplitude, not the phase), and no saved field of the model gives one, so `SOIL_LAG` is 30 days, the usual three to five weeks by which the observed continental maximum of surface air temperature follows the solstice. A deep layer's error from 10 days of lag is below 0.1 A exp(−z/d).

What the start leaves, measured with the layered column forced at its top by the same wave for a year (`test/landSoil.test.mjs`): the layers' means over the year depart from T̄ by at most 0.24 K (north), 0.16 K (south) and 0.55 K (firn), and the deepest layer ends the year 0.65, 0.42 and 2.39 K from where it started. A start at the skin temperature gives 0.79, 0.51 and 1.09 K, and 2.26, 1.46 and 4.56 K. What remains is the half-space's own misfit: the column's zero-flux bottom reflects the wave, so the deepest layer's amplitude is 27 percent (soil) to 42 percent (firn) larger than the half-space's, and the periodic solution of the finite column, T_k = T̄ + A Re[cosh(κ'(L − z))/cosh(κ' L) e^{iω(t − t_peak)}], would remove it.

Over the ice sheets `airCycle` is outside the data it was fitted to (cells off the ice sheets): on the N = 128 probe state of Aug 19 (`s6klm/arm-s0_day0882`) its annual mean is 256 K over Antarctica against a saved skin of 235 K, so the deep firn starts about 10 to 15 K warmer than the skin there; Greenland, the northern land and the southern land start within 2 to 11 K of the saved August skin in the directions the seasonal cycle has.

## The GPU

`js/gpu/soil.gpu.js` generates the constants of both classes and the step as straight-line WGSL on scalars (`colT2`, `colF1`, `colB2`, …) with no arrays, loops, pointers or module-scope private storage, which avoids the Metal compiler failures listed in `CLAUDE.md` by construction. `physicsSurfaceKernel(held, seaSkin, soil)` takes the step from its third argument; with `soil` null the kernel text is the slab's. The layers take `soilDeep` PH slots (`layoutFor`), travel through `uploadPhysics` and `uploadLand`, and come back with `pullLand` and `land.serialize()`.

## Tests

- `test/landSoil.test.mjs` (CPU): the column's constants and errors; the land energy row (soil and firn columns under a constant flux, a flux that changes sign with steps of one to five times 1350 s, melting snow with sublimation, and glacier melt); a uniform column unchanged without a flux; the stepped column against the backward-step solution of its own equations, step by step over a period, at 12 hours to a year; a one-layer column against the slab; the start against the formula (mean, amplitude and phase by depth, both hemispheres, firn); saved arrays used as saved; the drift over a year; the default slab's digest over 200 steps and its serialized keys; the state-file round trip and `regridLand`.
- `test/landSoilGpu.test.mjs` (GPU file): eight coupled steps at N = 16 against the CPU model beside a slab of the skin's capacity; one step with the GPU's own fluxes fed to the CPU land; the energy row on the GPU's pulled fields; and the default kernel text, constants and PH layout.
- Unchanged and passing: `land`, `landModel`, `landGpu`, `landIce`, `landStart`, `layeredDefault`, `regrid`, `stateFile`, `splitState` (its expected options gain `land: { soil: 'slab' }`), the four digest files, `boundaryLayer`, `gpuModel`.

## Measured

- Energy row: the CPU column's heat content changes by Q Δt + L_f Δ(snow) − L_f Δ(ice lost) to 1e-16 of the column's heat content at every step in all four cases; the GPU's pulled fields, with the net flux of the step, to 2.5e-8 (float32), over 1193 land cells of which 141 are on ice sheets.
- The CPU column follows the backward-step solution of its own equations to 1e-10 K at periods from 12 hours to a year in soil and firn; the skin's response to a sinusoidal flux is 0.98 (12 hours) and 0.9998 (10 days) of the continuous column's.
- Engines at N = 16, eight coupled steps: the skin differs by 7.9e-4 K (2.9e-4 K for a slab of the skin's capacity), the layers by 7.8e-5, 1.1e-4 and 1.3e-4 K at every land cell; one step with the GPU's own fluxes fed to the CPU land: 3.2e-5 K in the skin and 1.5e-5 K in the layers.
- Engines at N = 128 (`LAND='{"soil":"layers"}' node scripts/compareEngines.mjs` on `s6klm/arm-s0_day0882`, 1, 4 and 16 steps): the skin differs at all cells by a relative rms of 8e-8, 6e-7 and 7e-6 (maximum 0.001, 0.03 and 0.63 K); soilT2 by at most 3e-5, 1e-4 and 2.7e-3 K, soilT3 and soilT4 by at most 3e-5, 8e-5 and 2.5e-4 K. The tool's other regions are over its 1e-3 tolerance in the same families as with the slab (44 against 39 failures); the lighter skin lets a flux difference between the engines move the surface 6.5 times as far as the slab's 1e6 J/m²/K does, so SFLUX, sensible heat and OLR separate earlier (OLR at 4 steps 1.2e-3 against 1.9e-4 of the slab).
- Pace, N = 128 on the M5, 512 steps queued 8 to a submission, slab and layers alternating over three rounds: 39.4 and 39.3 s per model day in the first, 39.4 and 39.3 in the second, 42.9 and 41.6 in the third (a busier GPU); no cost can be told from the noise, under 0.3 percent in the quiet rounds.

## Not included

Snow insulation (IFS §8.4.3), soil freezing (§8.5.2), a conductivity and capacity that follow the soil's water (`soilConductivity: 'wetness'` from the bucket), the canopy's skin coupling (IFS Λ_sk) and a graded firn-to-ice column deeper than 2.9 m are separate changes.
