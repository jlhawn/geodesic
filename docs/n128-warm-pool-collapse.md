# The N=128 warm-pool collapse

Why spin-up thirteen's N=128 run lost the west Pacific warm pool from
day ~1316 while N=64 kept it, measured on 2026-10-20. The cause is
numerical: the ocean step counts an edge's transport two ways, and
every step's ALE cycle ratchets the difference into the barotropic
transport in the direction of the flow. A jet that is fast enough feeds
the ratchet faster than the ocean dissipates it.

## What the run did

Thirteen ran N=128 on a Verda A100 (bl42, STRATOSPHERE=0, ocean every
8 steps, a 1350 s ocean step) from WOA's annual ocean, continuously to
about day 436 and in five-minute segments after that. Spot preemptions
replaced the instance at about days 365, 770 and 1560. N=64 ran on the
Mac. The daily log line's "currents ≤ … m/s, transport … Sv" is the
fastest edge velocity and the largest single-edge column transport:

- Through day 1250, currents stay at 2.0–2.6 m/s and the transport at
  32–45 Sv, as at N=64 (≤ 2.3 m/s and mostly 35–57 Sv all run).
- From day 1253 the transport grows: 61 Sv on day 1258, 77 on day
  1261, 100 on day 1267, and 81–116 to day 1305. Currents rise from
  2.2 to 3.5 m/s by day 1297.
- Days 1306–1308 reach 4.1–4.5 m/s. One edge is held at the 5 m/s
  limit on day 1307 and ten on day 1308. From day 1309 the fastest edge
  sits at the limit, and by the 1320s hundreds are held.

The flow that grows is a broad eastward jet at 6–8°N, 128–142°E:
1.5–2.5 m/s at the surface on day 1261 and 3–3.5 m/s on day 1297, with
a depth-mean of 1.8 m/s at 7–8°N, 137–140°E. It is not grid noise. In
the box 0–12°N, 125–150°E the barotropic kinetic energy is 74 PJ on day
1242 and 468 PJ on day 1297.

The warm pool goes after the jet, starting where the jet does:

| day | global Ts (°C) | warm-pool SST, equator line (°C) |
|---|---|---|
| 1297 | 14.41 | 25.6 |
| 1316 | 13.48 | 25.5 |
| 1335 | 12.14 | 23.2 |
| 1354 | 11.42 | 19.8 |
| 1369 | 10.88 | 16.4 |
| 1407 | 10.32 | 12.1 |
| 1460 | 10.73 | 8.3 |

N=64 holds 24.8–25.7 °C on the same line from day 1278 to day 1825.

Between days 1297 and 1316 the 0–300 m mean cools by up to 14 K at
6.5°N, 140–145°E, and by 8–9 K at 1.5°N/1°S, 140°E. That is over the
ridges and trenches around Palau and Yap, not in the Sulu, Celebes or
Indonesian seas, which stay at 26–27 °C. Cell 100762 (6.92°N, 142.54°E,
3276 m deep) holds a 0–300 m mean of 19.8 °C on day 1297. By day 1316
the whole column is 1.4–1.55 °C water, η is −3.77 m, and |u| is 4.8 m/s
at every depth: a barotropic vortex at the speed limit that has stirred
the column to the bottom. The 0–300 m mean temperature over 5°S–5°N,
130–150°E falls from 21.9 °C (day 1278) through 12.7 (1335) to 9.1
(1354), with the salinity mixed to 34.73 psu down to 300 m.

## Not the segments, not the A100

`runs/n128-collapse/c128.log` reruns the unchanged model on the M5 Max
from `thirteen128_day1297.bin` without segments. It does the same thing
on the same days. The first edges are held on day 1308 (A100: day
1307), and the fastest edge reaches the limit on day 1309, as on the
A100:

| day | currents (m/s) | transport (Sv) | edges at the limit | warm pool (°C) |
|---|---|---|---|---|
| 1299 | 3.59 | 85 | 0 | 25.5 |
| 1305 | 3.73 | 106 | 0 | 25.8 |
| 1308 | 4.36 | 166 | 4 | 25.8 |
| 1309 | 5.00 | 176 | 14 | 25.9 |
| 1317 | 5.00 | 405 | 104 | 25.5 |
| 1323 | 5.00 | 514 | 215 | 25.2 |
| 1330 | 5.00 | 547 | 382 | 24.3 |

The equator line has the warm pool at 24.3 °C on day 1330 against the
A100's 23.2 on day 1335. The ocean alone does it too: stepped from day
1297 on the CPU with no wind and no surface flux, the jet box's kinetic
energy grows 4.7% in a day (below).

## Where the energy comes from

`scripts/oceanPhaseEnergy.mjs` steps the ocean alone from a saved state
(no stress, no surface flux, CPU, f64) and charges each phase of the
step with the kinetic energy it changes. The kinetic energy is
Σ ½ρ h̄ u² dc dv, with h̄ the mean of the two cells' thicknesses. From
day 1297, over 24 hours (PJ):

| phase | jet box | globe |
|---|---|---|
| layer step (dynamics) | −53.4 | −106.8 |
| rescale | +0.01 | +0.02 |
| eddy transport | +0.42 | +0.46 |
| barotropic shift | **+118.0** | **+132.3** |
| regrid | −5.9 | −7.5 |
| column physics | −6.5 | −29.1 |

The box goes from 1119.8 to 1172.5 PJ. Everything but the shift takes
energy out. The barotropic Coriolis term does −0.065 TW of work in the
box, so it is not the source.

The shift adds a depth-uniform velocity that puts each edge's column
transport on the barotropic solver's average. Its correction splits
into two parts:

- the barotropic tendency less the layers' own tendency;
- the change in the column transport from the edge thicknesses alone
  over the layer step, Σ (h_e,after − h_e,before) u.

As power against each edge's depth-mean flow, the shift does 0.97 TW in
the box over 11 steps from day 1297. The thickness part is 0.85 TW and
steady from step to step. The tendency part swings from −1.3 to
+0.36 TW between steps and averages −0.08 to +0.12 TW, depending on the
window. Split by phase, the thickness part is all from the layer step
(+0.90 TW; rescale −0.004, eddy transport −0.05). Split by layer, it
is the isopycnal classes' (0.84 TW; the twelve fixed-depth layers
0.010).

## The mechanism

The ocean counts the water at an edge in two ways.

- **The minimum rule (`edgeColumn`).** The top layer takes the mean of
  its two cells; each layer below takes the smaller of its two cells'
  thicknesses; the column is cut at the edge's sill. The continuity's
  volume fluxes, the barotropic solver's transport and forcing, and the
  shift all use this rule.
- **The cell mean.** The regrid's velocity remap (`regridEdges`) and the
  vertical viscosity (`columnEdges`) weight each layer by the mean of its
  two cells' thicknesses, so they conserve Σ h̄ u, not Σ h_min u.

Σ h_min u is concave in the interfaces. Raise an interface on one side
of an edge and lower it on the other: the mean thicknesses do not
change, but both layers' minimum thicknesses fall, and the column
transport by the minimum rule falls with them. Every ocean step then
runs a cycle:

1. The layer step moves the interfaces. Where they move differently on
   the two sides of an edge, the minimum-rule transport falls short of
   the barotropic average by some W.
2. The shift adds W/Σh_min to every layer, which puts the transport back
   on the barotropic average.
3. The regrid returns the interfaces to its targets. The minimum
   thicknesses recover. The remap carries the velocities, refund
   included, at constant Σ h̄ u, so the minimum-rule transport now
   exceeds the barotropic average by about W.
4. The next barotropic step starts from that transport.

A worked example, both cells [50, 100, 850] m and u = [1, 1, 0] m/s:

| | minimum-rule transport |
|---|---|
| start | 150 m²/s |
| layer step: the interface between the second and third layers 20 m deeper in one cell, 20 m shallower in the other; minimum thicknesses [50, 80, 830] | 130 |
| shift: +20/960 = +0.0208 m/s in every layer | 150 |
| regrid back to [50, 100, 850]: the mean thicknesses never moved, so the remap changes nothing | 170.8 |

No water moved, yet the edge carries 14% more and its deep layer is
now flowing. The refund always points along the flow, and it grows with
the flow: a faster jet moves its interfaces further each step. In the
unforced ocean the shift injects about 1% of the jet box's kinetic
energy per 12 hours at N=64 on day 1278 (1.2 of 127 PJ) and at N=128 on
day 1242 (2.6 of 316 PJ). There the box's energy holds about steady
(−0.8 and +0.5 PJ). By day 1297 the injection is 5% per 12 hours
(60 of 1120 PJ) and the box gains 30 PJ. The N=128 box held 316 PJ on
day 1242 against N=64's 127 PJ on day 1278. This round did not
separate the grid's own share in that from the run's history.

## The fix: `edgeTransport: 'kept'`

A change of vertical coordinate moves no water across an edge, and an
internal stress moves no momentum out of the column. With
`edgeTransport: 'kept'` (CPU `js/ocean/layered.module.js`, WGSL
`js/gpu/layeredOcean.gpu.js`), the regrid's velocity remap and the
vertical viscosity end with a depth-uniform velocity. That velocity
returns each edge's transport, counted by the minimum rule with its
sill cut, to its value before the operation. The barotropic solver then
starts every step from the transport it ended the last one with.

The remap and the solve are otherwise unchanged and still conserve
Σ h̄ u to within that one uniform velocity. `kept` is the default;
`edgeTransport: 'remapped'` runs the model as thirteen ran it.

Two other ways were tried and set aside:

- **Barotropic averaging over 2Δt**, centered at the new time, as in
  MPAS-O's split-explicit scheme. It cuts the shift's injection in the box from 59.7 to
  42.1 PJ in 12 hours, but the box still gains 14 PJ.
- **One partition everywhere.** The cell-mean thicknesses, cut at the
  sill (HYCOM's u-point thickness), in the barotropic solver, the shift,
  the remap and the solve. It stops the ratchet: the shift does −3.5 PJ
  in 12 hours and the box loses 16 PJ. But the continuity still moves
  water by the minimum rule, so each step the barotropic free surface
  and the layers' column sums disagree by 2.8 m rms (46 m at worst) in
  the box, against 0.3–0.5 m before. The rescale spreads that water
  through each column in proportion to the layers' thicknesses.

## The fix in the unforced ocean

From day 1297, in the jet box (PJ):

| | remapped, 12 h | kept, 12 h | kept, 36 h |
|---|---|---|---|
| box energy | 1119.8 → 1149.9 | 1119.8 → 1117.8 | 1119.8 → 1087.8 |
| its barotropic part | 467.5 → 491.7 | 467.5 → 460.6 | 467.5 → 438.0 |
| globe | 2443.2 → 2437.6 | 2443.2 → 2403.7 | 2443.2 → 2325.0 |
| columns against the free surface before the rescale, rms | 0.27–0.52 m | 0.035–0.056 m | 0.029 m |

N=64 from day 1278 over 24 hours (ocean step 2700 s), the same box:

| | remapped | kept |
|---|---|---|
| box energy | 127.1 → 125.4 | 127.1 → 124.4 |
| its barotropic part | 33.2 → 33.7 | 33.2 → 32.7 |
| globe | 1304.0 → 1216.1 | 1304.0 → 1211.5 |
| columns against the free surface before the rescale, rms | 0.064–0.090 m | 0.008–0.009 m |

The ratchet runs at N=64 too, but there it is small against the
dissipation. With `kept` the globe loses 0.35% a day more of its
kinetic energy than it does without.

At N=128 the shift still adds energy, 56.7 PJ in the box in 12 hours,
while the regrid's uniform correction takes 32.6 PJ out. The refund of step 2 and
its removal at step 3 nearly cancel. What remains is the barotropic
tendency the shift is there to carry. The disagreement between the
barotropic free surface and the layers falls tenfold, because the
barotropic solver and the continuity start each step from the same
transport.

## The fix in the coupled model

`runs/n128-collapse/k128.log`: the same rerun from day 1297 on the M5
Max with `OCEAN={"everySteps":8,"edgeTransport":"kept"}`, 163 days in
93 minutes. That is 34 s a day, the same as without the option.

| day | currents (m/s) | transport (Sv) | edges at the limit | warm pool (°C) | jet box (PJ) |
|---|---|---|---|---|---|
| 1297 | 3.53 | 94 | 0 | 25.6 | 1120 |
| 1315 | 2.09 | 49 | 0 | 25.9 | 728 |
| 1333 | 2.07 | 50 | 0 | 26.1 | 533 |
| 1369 | 1.74 | 46 | 0 | 25.8 | 384 |
| 1405 | 1.41 | 46 | 0 | 24.8 | 246 |
| 1460 | 1.36 | 38 | 0 | 25.4 | 201 |

The jet runs down within three weeks and keeps decaying. By day 1405
its box holds less than the 316 PJ of day 1242. No edge reaches the
limit, and the warm pool stays within N=64's range.

| on the same days | A100, unchanged model | N=64 |
|---|---|---|
| transport | 404–729 Sv | — |
| edges at the limit | 107 to 4992 | none |
| warm pool | 25.5, 23.2, 16.4, 12.1, 8.3 °C | 25.7 °C on day 1369, 24.8 on day 1460 |
| jet box | 85,400 PJ by day 1369 | 95 PJ on day 1369 |

On day 1460 the global surface temperature is 12.39 °C, against N=64's
12.33 and the A100's 10.73. Across 65°W the Drake Passage carries
87–113 Sv over the rerun, against N=64's 84 Sv on day 1369 and 81 on
day 1825.

N=64 barely notices the fix. The pair below runs 90 days from
`thirteen64_day1825.bin` (`runs/n128-collapse/n64r.log`, `n64k.log`),
with the ranges taken over the 15-day log lines and the equator lines:

| | remapped | kept |
|---|---|---|
| currents | 0.96–1.31 m/s | 1.10–1.43 m/s |
| transport | 39–49 Sv | 38–47 Sv |
| warm pool, equator line | 25.1–25.6 °C | 25.3–25.8 °C |
| undercurrent | 0.36–0.41 m/s | 0.33–0.39 m/s |
| Drake | 65–84 Sv | 72–80 Sv |
| jet box, day 1915 | 59 PJ | 54 PJ |
| global Ts, day 1915 | 16.00 °C | 16.04 °C |

## Tests

- `test/layeredShift.test.mjs`: with `kept`, the regrid and the vertical
  viscosity leave every edge's minimum-rule transport where the shift put
  it, to 2·10⁻¹⁵ of the largest. With `remapped` they move it by 6·10⁻⁵
  (regrid) and 9·10⁻⁷ (viscosity).
- `test/hybridOcean.test.mjs`: water, heat and salt conserved to
  rounding over 24 steps under both settings.
- `test/hybridOceanGpu.test.mjs`: at N=16 with `kept`, the GPU tracks the
  CPU through the regrid, the column physics and sixteen more steps.
  Its regrid and column physics hold every edge's transport to 1.0·10⁻⁶
  of the largest, float32 rounding.

## Open

- The barotropic average over [t, t+Δt] is installed as the state at
  t+Δt, half a step off center. That is a separate item; centering it
  alone does not stop the ratchet (above).
- The rescale scales each column's heat with its water, so the column
  mismatch it absorbs carries the column's own temperature profile. It
  is ten times smaller with `kept`.
