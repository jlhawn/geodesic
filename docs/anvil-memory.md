# Anvil memory: `moist.anvilLifetime`

The tropics of seventeen128 have almost no anvil. This note records what the
option `moist.anvilLifetime` is, what the first probes of it do at N = 64, and
what limits it.

## The gap

`scripts/cloudRegimes.mjs` on seventeen128 day 1095 (year 3, March), against
the values it prints beside each regime:

| Region | High cover | Observed | LWCRE (W/m²) | Observed LWCRE |
|---|---|---|---|---|
| Pacific ITCZ (5–12N, 160E–100W) | 0.13 | 0.45–0.6 | 13 | 40+ |
| Warm pool | 0.36 | 0.55–0.7 | 36 | |
| Global | 0.25 | 0.2–0.3 | | |

Total cover is 0.51 against 0.65, and the relative humidity over ice at
150–350 hPa is 0.41–0.53 in the convective regimes. The deep plumes do reach
anvil height: a one-step read over the warm pool and the ITCZ finds a plume in
98–100 % of the columns, top pressures at the quartiles near 200–280, 450–510
and 700 hPa, and 35–48 % of the columns topping above 400 hPa.

The cause is structural. The plume's flux form sends the change in q_t to the
vapor, so the condensate it detrains evaporates into a dry layer, and the
uniform condensation (LeTreut and Li; RH_c from 0.975 at the surface to 0.75
at the top) then re-diagnoses cloud from the layer's mean q_t at every step.
Ice detrained into a layer whose mean humidity is under RH_c is gone at once.
Nothing gives an anvil a lifetime.

## The scheme

With `anvilLifetime` τ (s; null, the default, leaves the model unchanged) each
layer above `cumulusK0` (above `shallowTop` at the highest surface pressure:
29 of the 48 bl48 layers) keeps an anvil cover a and a grid-mean anvil ice
q_a. Neither is advected. The header of `js/physics/moist.module.js` has the
full statement; in short:

- Source, in the deep plume. Where its mass flux M falls upward through a
  layer, the detrained flux D = M(lower interface) − M(upper interface)
  carries the ice share 1 − α(T) of the condensate the plume holds at the
  lower interface after the rain made there. g Δt D q_i,u / Δp, at most the
  layer's vapor, leaves the vapor increment of the flux form for q_a, the
  layer warms by L / c_p times it, and the cover rises by (1 − a) g Δt D / Δp
  (Tiedtke 1993), times `anvilSpread` (1 by default). The column's c_p T + L
  q_v and its q_v + q_c + q_a stay exact.
- Sinks, in the rain's pass, top down. Ice falling from the layer above joins
  a layer's, with the larger of the two covers. Cover and ice then decay by
  exp(−Δt/τ), the ice lost sublimating into the layer's vapor, and the share
  v Δt / (Δz + v Δt) of the ice leaves downward, v = `anvilFall`
  (ρ q_a / a)^`iceFallExponent` with `anvilFall` defaulting to `iceFall`. The
  cover stays, so a cloud thins as its ice falls out and falls more slowly as
  it does. A layer at or above `liquidTemperature` passes all its ice down;
  ice under 1e-8 kg/kg or cover under 1e-6 sublimates whole. What leaves the
  lowest anvil layer joins the large-scale ice stream, which melts,
  sublimates or snows it.
- Radiation. The cumulus arrays hold a and q_a / a above `cumulusK0`, so the
  radiation takes the anvil as it takes the plumes' cumulus below: the larger
  of the resolved and the anvil cover, the waters added. On the GPU the same
  cloud is read from the anvil's own buffers by `cumulusCloud`, in the
  physics, radiation and layered kernels, so no cumulus buffer grows.
- State. `anvilCover` and `anvilIce` (K0 × C) are saved and restored by
  `scripts/spinup.mjs`, `fixedSst.mjs`, `explainOcean.mjs` and the page's
  worker (`placeAnvil`); a state without them starts the anvil at zero, and a
  state that records `anvilLifetime` in its header options starts the page
  with it. `scripts/cloudRegimes.mjs`, `compareEngines.mjs` and the figure
  scripts that rebuild a model from a state (`gpuModelFrom`) place it too; a
  script that rebuilds a CPU model by its own code starts the anvil at zero.
- Memory. On the GPU the anvil is two K0 × C float fields: 2 × 29 × 163,842 × 4
  B = 38 MB at N = 128 on bl48.

`anvilFall` and `anvilSpread` move the two quantities the probes below found
to decide the result, the fall speed of the ice and the cover the detrained
air makes; at their defaults (`iceFall`, 1) the scheme is the one above.

## The reference lifetime

The decay is the e-folding of the cover that is left once the ice has fallen.
In trajectories of air detrained from convective peaks in storm-resolving
E3SM over the western Pacific an anvil lasts about 15 h (Gasparini et al.
2021, JGR Atmospheres 126, e2020JD033487), and satellite tracking of tropical
cirrus systems gives 19 to 30 h (Luo and Rossow 2004, J. Climate 17, 4541).
An e-folding of 6 h leaves a tenth of the cover after those 15 h, and it is
the reference here; 3 h (the "few hours" of a mesoscale system's stratiform
anvil) and 12 h bracket it.

## The probes at N = 64

Each run is `scripts/spinup.mjs` at N = 64 on bl48 (`STRATOSPHERE=1
LEVELS=bl48`) from seventeen64 day 365 (1 March) for two days on the GPU. The
control has no `MOIST`, the others the options named. `scripts/cloudRegimes.mjs`
then reads each day-367 state after one CPU step with the anvil restored from
the state and counted as its own class: "high + anvil" is the cover of the
resolved cloud with a top above 440 hPa and the anvil above 440 hPa together,
and the cloud effects (W/m²) are the state's day means. The warm pool is
10S–10N, 120–170E over sea; the ITCZ is 5–12N, 160E–100W.

| Run | Warm pool high + anvil | LWCRE / SWCRE | ITCZ high + anvil | LWCRE / SWCRE | Global high + anvil | LWCRE / SWCRE | OLR / ASR |
|---|---|---|---|---|---|---|---|
| control (no anvil) | 0.304 | 34.8 / -59.0 | 0.042 | 5.5 / -44.3 | 0.252 | 24.0 / -50.3 | 235.1 / 241.1 |
| τ 3 h | 0.310 (0.273 + anvil 0.092) | 35.8 / -61.2 | 0.049 (0.043 + anvil 0.010) | 6.3 / -44.9 | 0.239 (0.233 + anvil 0.029) | 23.3 / -50.6 | 235.8 / 240.8 |
| τ 6 h (reference) | 0.333 (0.266 + anvil 0.159) | 37.4 / -61.9 | 0.053 (0.045 + anvil 0.017) | 6.7 / -45.0 | 0.238 (0.228 + anvil 0.044) | 23.4 / -50.7 | 235.8 / 240.7 |
| τ 12 h | 0.364 (0.258 + anvil 0.241) | 40.1 / -63.6 | 0.058 (0.045 + anvil 0.026) | 7.3 / -45.1 | 0.239 (0.223 + anvil 0.058) | 23.6 / -50.9 | 235.7 / 240.5 |
| τ 6 h, anvilFall 0.5 | 0.356 (0.293 + anvil 0.174) | 41.8 / -71.3 | 0.052 (0.044 + anvil 0.021) | 7.0 / -45.5 | 0.249 (0.239 + anvil 0.048) | 25.1 / -53.0 | 234.1 / 238.4 |
| τ 6 h, anvilFall 0.1 | 0.374 (0.318 + anvil 0.176) | 43.8 / -75.3 | 0.055 (0.047 + anvil 0.022) | 7.3 / -45.9 | 0.258 (0.248 + anvil 0.049) | 26.4 / -54.3 | 232.7 / 237.1 |
| τ 6 h, anvilSpread 5 | 0.423 (0.298 + anvil 0.328) | 46.2 / -65.5 | 0.068 (0.053 + anvil 0.035) | 8.2 / -44.2 | 0.249 (0.227 + anvil 0.084) | 24.8 / -51.4 | 234.5 / 240.0 |
| τ 6 h, anvilSpread 5, anvilFall 0.5 | 0.463 (0.343 + anvil 0.378) | 54.1 / -77.9 | 0.071 (0.058 + anvil 0.046) | 9.3 / -45.5 | 0.263 (0.242 + anvil 0.096) | 27.1 / -54.1 | 232.2 / 237.3 |

Day-mean cloud effects of one day on a state two days from the control carry
the weather of that day: regions smaller than the warm pool (the South
Pacific convergence zone, for one) move by 20 W/m² with an anvil path under
1 g/m² to explain it. The trends down a column of the table are the signal.

The second table takes the weather out. It is the radiation of each state with
the anvil in and with it hidden, the longwave and the sunlight absorbed, the
mean of four times of day, in W/m²: the anvil's direct effect on the
planet, before the atmosphere answers.

| Run | Global LW / SW (net) | Warm pool LW / SW (net) | ITCZ LW / SW (net) |
|---|---|---|---|
| τ 3 h | +1.6 / -1.7 (-0.1) | +6.5 / -5.2 (+1.4) | +0.9 / -0.3 (+0.6) |
| τ 6 h (reference) | +2.2 / -2.3 (-0.0) | +10.3 / -7.6 (+2.7) | +1.3 / -0.4 (+0.9) |
| τ 12 h | +3.1 / -2.9 (+0.1) | +15.3 / -10.4 (+4.8) | +1.8 / -0.6 (+1.2) |
| τ 6 h, anvilFall 0.5 | +2.6 / -3.4 (-0.9) | +10.4 / -11.8 (-1.4) | +1.6 / -0.7 (+0.9) |
| τ 6 h, anvilFall 0.1 | +2.7 / -3.9 (-1.2) | +10.0 / -13.4 (-3.5) | +1.7 / -0.9 (+0.8) |
| τ 6 h, anvilSpread 5 | +3.9 / -3.2 (+0.7) | +19.1 / -10.9 (+8.1) | +1.9 / -0.5 (+1.4) |
| τ 6 h, anvilSpread 5, anvilFall 0.5 | +4.4 / -4.4 (-0.1) | +20.2 / -16.6 (+3.6) | +2.4 / -0.9 (+1.5) |

The plumes of the control's step: in the warm pool 886 of 895 columns fire,
518 of them with a top above 400 hPa and a mean base flux of 1.1e-2 kg/m²/s;
in the ITCZ box 687 of 687 fire, 53 of them (8 %) above 400 hPa, at 1.7e-2.

## What the probes say

- At the reference 6 h the anvil adds cloud where deep plumes top high. The
  warm pool's high cover goes from 0.304 to 0.333 (the anvil 0.159 of it, the
  resolved cloud losing 0.038), its longwave cloud effect from 34.8 to 37.4,
  and the anvil's direct effect there is +10.3 W/m² of longwave and −7.6 of
  sunlight. The ITCZ box gains 0.011 of cover and 1.2 W/m²: 92 % of its
  plumes top below 400 hPa, so there is little ice to detrain at anvil height,
  and no lifetime changes that.
- Globally nothing is gained. The anvil's direct effect is +2.2 W/m² of
  longwave and −2.3 of sunlight, but the high cover goes from 0.252 to 0.238
  and the longwave cloud effect from 24.0 to 23.4, so the resolved high cloud
  has lost 0.024 of cover and about 2.8 W/m² of longwave (the −0.6 less the
  direct +2.2). That fits the ice the flux form used to detrain into a layer
  near saturation, which the uniform condensation kept as resolved cloud
  with its own cover (sqrt(q_c / b), a third of the cell or more), going to an
  anvil of a few per cent cover instead, where the same ice is weaker
  radiatively; the probes did not isolate it. Over land and the middle
  latitudes the loss exceeds the gain (land 60S–60N: high + anvil 0.244 to
  0.217).
- The lifetime τ moves the cover more than the ice. By the large-scale ice's
  law the anvil's ice falls at about 0.5 m/s at an in-cloud content of
  0.1 g/kg, a layer of 500 m in a quarter of an hour, so the ice is gone from
  the anvil in an hour or two and the cover outlasts it as thin cloud for τ:
  the warm pool's anvil cover is 0.092, 0.159 and 0.241 for 3, 6 and 12 h, its
  anvil path 8, 12 and 16 g/m².
- The fall speed sets the ice path and with it the sunlight. With `anvilFall`
  0.5 and 0.1 the global anvil path doubles (5 to 10 and 12 g/m²), the direct
  longwave effect hardly changes and the sunlight effect grows from −2.3 to
  −3.4 and −3.9, so the net turns negative.
- The spread sets the cover and with it the longwave. With `anvilSpread` 5 the
  warm pool's anvil cover doubles to 0.328 and the direct effect there is
  +19.1 W/m² of longwave against −10.9 of sunlight, net +8.1; with the fall at
  0.5 as well the warm pool's high + anvil cover is 0.463 and its longwave
  cloud effect 54.1, against 0.304 and 34.8 in the control, with the sunlight
  at −77.9 against −59.0. The ITCZ box reaches only 0.071.

So the option at its defaults gives the anvil a lifetime and a place to be
radiatively, but the amount of cloud it can make is bounded by two things it
does not control: how much ice the plumes carry to anvil height, which in the
ITCZ box is little at N = 64, and how far that ice spreads (a cover of a
detrained updraft's own area is a few per cent), which `anvilSpread` stands
in for. Neither the fall speed nor the spread has a measured value here; the
runs bracket them.

## GPU cost

`scripts/paceGpu.mjs` at N = 64 on bl48, 512 steps in batches of 8, the GPU
to itself, two pairs alternated, from the control's and the 6 h run's day-367
states: 22.4 and 22.5 ms a step without the anvil, 22.6 and 22.6 ms with it
(5.7–5.8 s a model day for both), about 0.1 ms or 0.5 % more. The anvil's two
fields are 2 × 29 × C floats of the physics buffer, 38 MB at N = 128 on bl48.

## Tests

- `test/anvil.test.mjs` (CPU, 13 tests): the option and its refusals; the
  plume's source against a twin run without the anvil (the vapor gives up the
  ice exactly, the layer warms by L / c_p times it, water and c_p T + L q
  equal to 1e-13); the cover as g Δt D / Δp and `anvilSpread`; decay, trace,
  fall, the large-scale stream and melting each with water and enthalpy exact
  to 1e-12; `anvilFall` leaving the large-scale ice alone; the radiation
  seeing the anvil; `placeAnvil`.
- `test/atmosphereConservation.test.mjs`: over a spun-up slab with the anvil
  the adjustment changes c_p T + L q_v + S by 2.9e-8 W/m² at most and the water
  of vapor, cloud, anvil and rain by 1.5e-14 of itself.
- `test/gpuModel.test.mjs`: two parity cases at N = 16, cam26 and bl48, cell by
  cell outside the columns where a discrete decision parted (cover rms
  2.8e-4 and 1.8e-4 of itself, ice 6.9e-4 and 2.1e-4, OLR within 3.5e-2 W/m²).
- `scripts/compareEngines.mjs` on the 6 h probe's day-367 state (N = 64, bl48),
  one and two steps: the GPU's anvil cover and ice differ from the CPU's by
  1.1e-3 and 2.7e-3 (step 1) and 5.1e-4 and 1.0e-3 (step 2) of themselves in
  rms, and fewer than 6e-5 of the cells part by more than the tolerance; the
  rows it flags besides (the lake region unmapped, the large-scale rain, the
  sea ice concentration, the discharge, the mixed-layer top, the stratiform
  share) are not the anvil's.
- `test/anvilGpu.test.mjs`: the PH layout, the constants, the anchors, the
  page's cloud totals, and the kernels' text: with the option off the
  adjustment, radiation, physics and layered kernels hash as they did before
  the option.
- `test/gpuReload.test.mjs`, `test/spinupRestart.test.mjs`,
  `test/cumulusLoader.test.mjs`, `test/parallel.test.mjs`: the device buffers
  after a load word for word, a two-day run against a day and a day bit for
  bit, the page's worker, and three worker threads against one thread bit for
  bit, each with the anvil on.
- Off, the model is the parent's: the CPU digests of five bl48 variants after
  40 steps at N = 16, the pinned digest tests, and the WGSL of seven GPU
  configurations at N = 16 are unchanged from the commit before the option.

## Not done

- The anvil is not advected and has no deposition: it keeps no ice against
  the fall except by the plume's next detrainment, and the upper troposphere's
  humidity over ice stays where the rest of the model puts it.
- No N = 128 run and no run longer than two days has used it; the two-day
  probes start from a state spun up without the anvil, which the anvil's own
  feedbacks (heating, the plumes' response) have not had time to change.
- `anvilLifetime`, `anvilFall` and `anvilSpread` have no tuned values.
