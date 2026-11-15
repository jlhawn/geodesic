# Deposition-limited ice fall (`moist.iceDeposition`)

Under the ramp (the default), a cloud layer's ice share is 1 − α(T), linear from 0 at liquidTemperature (273.15 K) to 1 at iceTemperature (235.15 K). That share falls at the Heymsfield–Donner speed, and the liquid share α converts to rain.

With `iceDeposition` (false by default; it needs `iceFall`), a layer between those temperatures handles its cloud differently:
- What falls is the layer's cloud ice q_i. This is a per-layer memory that only deposition grows, following the Wegener–Bergeron–Findeisen process as the IFS models it.
- The liquid q_l = q_c − q_i converts to rain at the plain liquid rates.

Saturation, latent heat and the radiation keep the ramp. The option changes only how a layer's condensate is split between the sink that falls and the sink that rains out.

## Equations

These follow IFS Documentation Cy41r2 Part IV §7.2.4(d), eqs 7.33–7.43 (after Rotstayn et al. 2000 and Meyers et al. 1992). Each step in every layer with iceTemperature < T < liquidTemperature and q_l > 0, after the ice falling in has been added:

- Vapor pressures: e_sl is Bolton's (`saturationVaporPressure`) and e_si the IFS form (`iceVaporPressure`). The air density is ρ = p/(R T).
- Active nuclei: N_i = 1000 exp[12.96 (e_sl − e_si)/e_si − 0.639] m⁻³ (7.40). The CLOUDSC code divides by e_sl. Eq. 7.40 and Meyers divide by e_si, which is what makes f_nuc reach 1 at −23 °C as the text says.
- Conduction and diffusion terms:
  - A″ = (L_s/(K_a T)) (L_s/(R_v T) − 1)
  - B″ = R_v T p/(2.21 e_si)
  - L_s = L + L_f and K_a = 2.4e-2 W m⁻¹ K⁻¹ (7.34, 7.35).
- Growth coefficient: c_svd = 7.8 (N_i/ρ)^(2/3) (e_sl − e_si)/(ρ_i^(1/3) (A″ + B″) e_si), with ρ_i = 700 kg m⁻³ (7.38).
- In-cloud ice: a is the layer's cover, the uniform distribution's cover that the fall speed already uses. The in-cloud ice, with the seed as its floor, is q₀ = max(q_i/a, m₀ N_i/ρ), where m₀ = 1e-12 kg (7.41).
- New in-cloud ice: q₁ = ((2/3) c_svd Δt + q₀^(2/3))^(3/2). This is exact over the step for a constant c_svd (7.39).
- Deposition: D = F_dep · min(a (q₁ − q₀), q_l). D moves from q_l to q_i; q_c is unchanged.
- Cloud-top reduction:
  - F_dep = min[f_nuc + (1 − f_nuc)(0.1 + Δz_top/500 m), 1], with f_nuc = min(N_i/15000, 1) (7.42, 7.43).
  - Δz_top is 0 in a cloud-top layer: one whose cover is at least 0.01 under a layer with less.
  - Below that, Δz_top grows by each lower layer's thickness, as CLOUDSC's ZCLDTOPDIST does.

Time to grow 0.1 g/kg from the seed (a = 1, F_dep = 1, 600 hPa):

| T | 268 K | 263 K | 258 K | 253 K | 248 K | 240 K (400 hPa) |
|---|---|---|---|---|---|---|
| N_i (m⁻³) | 1050 | 2110 | 4380 | 9410 | 20900 | 79300 |
| Time | 2.8 h | 69 min | 37 min | 23 min | 16 min | 5 min |
| F_dep in the top layer | 0.16 | 0.23 | 0.36 | 0.66 | 1 | 1 |

## The memory

`moist.cloudIce` is the cloud ice of each layer, K·C values in kg/kg. It is a Float64Array on the CPU (in a shared buffer for the workers) and PH region `CICE` (f32) on the GPU. It is allocated only with the option, so the option-off layout and saves are unchanged.

Invariants:
- 0 ≤ q_i ≤ q_c.
- q_i = q_c at or below iceTemperature (homogeneous freezing).
- q_i = 0 at or above liquidTemperature.

Only three things change q_i:
- deposition;
- the share of falling ice that a layer keeps as cloud;
- the fall itself.

Every other process that changes q_c changes its liquid: condensation, evaporation, mixing, transport, cumulus detrainment, conversion and the filler. q_i is capped at q_c when the rain step reads it, and again after the column's second adjustment and filler.

The memory is not advected:
- In a steady cloud field each cell's q_c holds, so the local memory tracks that cloud's glaciation as if the cloud stood still.
- Transients at cloud edges err toward liquid. Condensate that arrives is liquid, and ice left behind is capped away.
- Advecting it would add a fourth moist tracer to both dynamical cores and their conservation rows.

New condensate between the two temperatures is liquid, as in the IFS, and so is detrained cumulus condensate. The IFS detrains with its ramp. Here detrained condensate glaciates at the deposition rate: 16–37 min between 248 and 258 K.

## Saves and restores

- `scripts/spinup.mjs` writes `cloudIce` when the model has it and restores it on resume with the other moist memories.
- The page saves it with the cumulus memories (`savedCumulus` in `js/model.worker.js`), and every restore places it through `placeCumulus`.
- The GPU model mirrors `CICE` into `moist.cloudIce` when it reads back for a save, and uploads it with the retained physics on load.
- A save without `cloudIce` restored under the option starts the memory at 0. All mixed-phase cloud then begins as liquid, as a fresh start does.
- compareEngines maps `CICE` to `moist.cloudIce`.

## CPU (f64) and WGSL (f32)

- Both engines use the same expressions. q₀ > 0 always, so f32 `pow` is defined. The subtraction q₁ − q₀ loses about 1e-7 of q₀, under 1e-11 kg/kg.
- The GPU runs one thread per column with a top-down loop. It reads and writes the memory in the PH buffer at each layer and carries the cover above and Δz_top as scalars. There are no private arrays and no pointers, per CLAUDE.md's M5 and M1 rules.
- Three decisions can split the engines at a few cells, so compareEngines will judge them where they show:
  - the cloud-top test (cover 0.01);
  - the two temperature limits;
  - min(·, q_l).

## Where it runs in the moist step

The memory acts in `autoconvertColumn` (CPU) and the `adjustRain` kernel (GPU). It runs top down, per layer, inside the iceFall block, after the layer's rain evaporation and the plume's stream:

1. **Ice falling in.** The receiving layer's share α melts into the rain, as now. The rest joins q_c and q_i.
2. **Cloud top.** The tracker updates Δz_top for every layer, clear ones included.
3. **Deposition.** The growth above.
4. **The fall.** q_i falls at v = iceFall (ρ q_i/a)^iceFallExponent, implicit and in flux form, as now. The layer keeps 1/(1 + v Δt/Δz) of it and passes the rest down. The ice leaving the lowest layer is large-scale precipitation.
5. **Conversion.** q_l converts by Kessler and the lifetime, with no Bergeron factor. In floored layers the fog scheme acts on q_l.

Then come the second adjustment (the fallen ice sublimates in subsaturated layers), the filler and the cap.

## The liquid a layer holds back

Held liquid has three sinks:
- deposition onto the ice, which then falls;
- Kessler autoconversion and the lifetime at the plain rates: threshold 2e-4, 1e-3 /s, 1 h, or 3 h for stratiform cloud;
- evaporation.

There is no Bergeron factor:
- Cy41r2's large-scale liquid autoconversion has none (7.49–7.52). The factor survives only in its convection (6.37), because the deposition now is the Bergeron–Findeisen process.
- Keeping `bergeronAutoconversion` under the option would count that process twice. The option therefore sets the factor to 1 in the layers it splits. The factor matters only between 235.15 and 268.16 K, so this is exactly those layers.

coldfall64 (`iceFallColdOnly`) converted its held cloud with the factor:
- Below 250 K the factor is 1 + 0.5 √18 = 3.1, so the 1 h lifetime became 19 min and the 3 h stratiform lifetime 58 min.
- That fast conversion, rather than the missing fall, is the likely cause of its LWCRE loss. It is untested.

What I expect under the option:
- **235–255 K:** deposition glaciates 0.1 g/kg within 5–30 min, so the ice falls much as in the default. Tropical mid-level cloud should stay near ctl64, which is the LWCRE guard.
- **255–273 K:** the liquid lives longer than in the default, with no factor and slower glaciation. Boundary-layer mixed-phase cloud should hold more supercooled liquid, which is the Arctic rows.
- **Risk:** global SWCRE gets stronger. It is already −51 W/m² against CERES's −45.

An attribution arm, `bergeronAutoconversion: false` alone under the same protocol, would separate the factor from the deposition.

## Unchanged

- The radiation's cloud phase stays on the ramp. Coupling it to the memory would change both shortwave and longwave and is a separate step.
- The thermodynamics: the saturation mix and the latent heat L + (1 − α) L_f.
- The plumes' own ice.

## Tests

- `test/cloudIce.test.mjs`:
  - The option is off by default and needs iceFall. A column whose cloud is all colder than 235.15 K or all warmer than 273.15 K steps under it as without it, bit for bit.
  - A mixed-phase cloud-top layer's deposition, fall and liquid conversion match a hand computation. The Bergeron factor changes nothing under the option.
  - The glaciation times of the table above.
  - `placeCumulus` places a saved cloud ice and starts it at zero without one.
  - A GPU parity row at N=16 at the file's tolerances; the engines' cloud ice agrees to 1e-7 kg/kg.
- `test/gpuReload.test.mjs`: under the option, a load restores every device buffer word for word, the cloud ice included.
- `test/atmosphereConservation.test.mjs`: the moist row "over a slab with iceDeposition".
- CLAUDE.md's M5 checks: `test/boundaryLayer.test.mjs`, `test/gpuModel.test.mjs`, and compareEngines at the page state (combined128 day 1460) with and without the option.
