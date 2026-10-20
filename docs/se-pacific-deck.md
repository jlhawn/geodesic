# The SE Pacific deck

Why the mixed-layer stratocumulus deck runs on so little of the SE Pacific,
measured on spin-up thirteen's N=64 states and on the fixed-SST arms of
`docs/trade-cumulus-cover.md` (thirteen64_day1825 over WOA23, bl42).

Spin-up thirteen was later found to carry the hybrid ocean's edge-transport
ratchet (`docs/n128-warm-pool-collapse.md`, on
claude/priceless-dijkstra-3ce393), and nothing new is to be based on its
states. Rerun the N=64 numbers here from a clean state before relying on
them. The N=128 arms start from twelve.

## Answer

The deck's gate fails because the resolved boundary layer is 0.5–1 km too
deep and decoupled, and the inversion that caps it is held up by the shallow
plume. The model has two branches, and the gate picks between them:

- **Deck branch** (September, thirteen's own day 1643): the gate is open on
  0.73 of the Peruvian box (Klein and Hartmann's 10–20S 80–90W). The open gate
  vetoes the plumes, and the boundary layer's K-profile, laid over the deck's
  depth, supplies the moisture at the layer's top. The inversion is at 1.6 km
  with a 13.5 K jump, and the low cover is 0.94.
- **Plume branch** (March and April, every arm): the gate is shut. The shallow
  plume supplies 80–85 percent of the physics' moistening at the moist layer's
  top. The inversion rises to 1.8–2.3 km in the Peruvian box and to 2.7 km at
  20S 85–75W. There, the deck's starting height sits about 1 km below the jump
  in coupled columns, and decoupled columns stand down under a ceiling above
  2 km.

Neither branch is the observed one. The SE Pacific's boundary layer is
deepest in March and April. Even then its inversion lies between 925 and
850 hPa, and the deck covers 0.57 of 10–30S 75–95W.

Once the gate shuts, nothing brings the column back. Under the same WOA23
SST, the arms move further onto the plume branch over ten days. The plume
fires at the full trade-cumulus closure under an inversion of
stratocumulus strength, which is what operational models switch off. The
IFS does it with an EIS threshold of 7 K (below), and the leading change
follows it. In the arms, the switch holds a deck whose EIS stays above 7 K.
It does not reopen the MAM deck, whose EIS has already fallen to the
threshold (The switch in the arms).

## Instrument

`scripts/deckGate.mjs STATE` (CPU, ocean off, `STEPS` steps, 4 by default;
`RADIATION`, `MOIST`, `BOUNDARY_LAYER` and `SURFACE` work as elsewhere) wraps the
physics phase. Before each step it replays the deck's gate per column
exactly as `shadowDeck` takes it: the ceiling, the stand-down, the resting
and carried height, the fractional slab, the smoothed and running-mean
subsidence, the free troposphere, `mixedLayer.diagnose` and the gate's
memory. Over 3,388 column-steps its decisions differ from the model's in
none.

For the boxes (Peruvian 10–20S 90–80W; SE Pacific 10–30S 110–80W; Peru
5–20S 90–75W; 20S 85–75W, VOCALS-REx's line; Namibia 10–20S 0–10E;
California 20–30N 130–120W; ice-free sea cells only) it prints:

- each column-step's class and the shares that fail each test alone;
- the regimes;
- the resolved inversion (the interface of largest θ_v jump between 300 m
  and 3 km), the ceiling and the lowest cloudy run;
- the boundary layer's depth, mixing top, decoupling height, entrainment and
  cloud-top cooling;
- the plume's share, base flux and top;
- the subsidence at the inversion and at 700 hPa;
- the deck's start, its jump and its diagnosis;
- EIS, LTS, the radiation's low cover and water, the surface fluxes and SST;
- the box's mean profile below 4 km.

`scripts/tropicalHeating.mjs` takes its boxes from `BOXES` (JSON) for the
process budgets below. The ERA5 numbers are the WeatherBench 2 1990–2019
climatology (128×64, 13 levels), averaged over days of year 85–99 to match
the arms, and over October. They are interpolated in ln p like the model's
layers.

## The gate in the boxes

Peruvian box, shares of column-steps (four CPU steps from each state):

| state | runs | stood down | off: jump | resolved inversion, km | EIS, K | plume on | low cover |
|---|---|---|---|---|---|---|---|
| thirteen64 day 1643 (Sept, own SST) | 0.73 | 0.00 | – | 1.60 | 11.3 | 0.28 | 0.94 |
| thirteen64 day 1734 (Dec) | 0.57 | | | 1.86 | 7.6 | | 0.83 |
| thirteen64 day 1825 (March) | 0.00 | | | 2.60 | 2.9 | 0.75 | 0.43 |
| defaults arm, day 1835 | 0.28 | 0.14 | 0.53 | 1.83 | 6.1 | 0.56 | 0.75 |
| defaults arm, day 1840 | 0.05 | 0.33 | 0.55 | 2.01 | 5.8 | 0.72 | 0.68 |
| defaults arm, day 1845 | 0.06 | 0.55 | 0.38 | 2.32 | 4.0 | 0.71 | 0.61 |
| variance pair, day 1835 | 0.62 | | | 1.87 | 7.0 | | 0.71 |
| variance pair, day 1840 | 0.49 | | | 1.84 | 5.5 | | 0.53 |
| variance pair, day 1845 | 0.18 | 0.42 | | 1.95 | 5.6 | 0.64 | 0.26 |

The plume share counts all four steps. The plume's fields are not in the
saved states, so the first step reads none: over the three later steps the
share on day 1845 is about 0.94.

**The subsidence test never fails.** The ring-meaned sink at the resolved
inversion is 3.1–7.6 mm/s across the arms and days. The deck's running mean
at its start is 3.7 mm/s on day 1845.

**Why the gate fails on day 1845 of the defaults arm:**

- Coupled columns (0.34) start the deck at or near their boundary-layer
  top. The box's mean start is 1.0 km, the mean jump there 0.12 K, and the
  4 K jump about 1 km higher.
- Of the decoupled (0.57) and surface-driven (0.09) columns, 0.46 of the box
  stands down under a ceiling above 2 km and 0.09 has no ceiling.

The other boxes on day 1845:

| box | runs | resolved inversion, km |
|---|---|---|
| SE Pacific | 0.01 | 2.37 |
| Peru | 0.07 | 2.08 |
| 20S 85–75W | 0.00 | 2.73 (sink 1.2 mm/s that day) |
| Namibia | 0.23 | 1.26 (EIS 2.4) |
| California | 0.78 | 1.56 (EIS 10.2) |

Where the Peruvian deck runs, it holds 300 g/m² of water against Earth's
40–150 g/m². Its h is 1.67 km, its slab's cloud base 1.13 km and its
entrainment 9.6 mm/s. The radiation takes at most 150 g/m².

## Against observation

| | model, Peruvian box | Earth |
|---|---|---|
| inversion, MAM | 1.8–2.3 km (arms); 2.6 km (own, March) | between 925 and 850 hPa (ERA5, days 85–99); the year's deepest (COSMIC, Chan and Wood 2013) |
| inversion, SON | 1.60 km | 1.14 ± 0.10 km at 80–90W 15–25S (MODIS, Wood and Bretherton 2004); 1.0, 1.2 and 1.6 km at the coast, 75W and 85W along 20S (Bretherton et al. 2010) |
| inversion at 20S 85–75W, MAM | 2.73 km (defaults, day 1845) | 1.2–1.6 km in Oct–Nov; 1.2–1.7 km at 85W in Dec 2004 |
| jump | 6.2 K (MAM), 13.5 K (SON) | 10–12 K in Oct–Nov (Bretherton et al. 2010) |
| EIS | 4.0–7.0 K (MAM arms), 11.3 (SON) | 7.3 K (ERA5, days 85–99), 11.9 (ERA5, Oct); 6.9 annual (Wood and Bretherton 2006) |
| subsidence at the inversion | 3.1–7.6 mm/s | 2.6 ± 0.9 mm/s (Wood and Bretherton 2004); 2–4 mm/s (Wood 2012) |
| entrainment | equals the subsidence at the inversion in steady state | 2.7 ± 1.2 mm/s (Wood and Bretherton 2004); 4 ± 1 mm/s at 20S 85W (Caldwell et al. 2005) |
| cloud-top longwave cooling | 24 W/m² (MAM, the boundary layer's ΔF), 52 (SON) | 70 W/m² nearshore, 100 offshore (Bretherton et al. 2010); 50–90 (Wood 2012) |
| liquid water path | 45 g/m² of low water (MAM); 300 in the deck where it runs | 40–150 g/m² (Wood 2012); 67 at the coast, 166 offshore (Zuidema et al. 2012) |
| decoupling | 0.57 of columns; z_i − LCL about 1.8 km | about half of the clouds, with z_i − LCL 0.3 km at the coast to 1 km at 85W (Jones et al. 2011) |
| low cover | 0.61–0.75 (defaults, days 1835–1845), 0.26–0.71 (pair) | MAM 0.57, SON 0.70 (CloudSat–CALIPSO 2006–11, 10–30S 75–95W, Muhlbauer et al. 2014) |

ERA5's Peruvian column on days 85–99:

| level, hPa | θ, K | RH |
|---|---|---|
| 1000 | 294.7 | 0.73 |
| 925 | 295.3 | 0.88 |
| 850 | 301.0 | 0.52 |
| 700 | 315.2 | 0.16 |

The arms against it:

- At 850 hPa they are 2–2.8 K colder and 2.2–2.9 g/kg moister (RH
  0.84–0.87). They lie inside the boundary layer.
- At 1000 and 925 hPa they are 1.1–1.3 K and 1 K warmer. That is the
  deepening–warming of a decoupled layer (Bretherton and Wyant 1997), and it
  is what lowers their EIS.
- At 700 hPa they match ERA5 within synoptic swings (315.1 to 313.5 K over
  the ten days).

Over the whole eastern Pacific (30S–0, 120–70W), the model's 700 hPa stays
within 0.4 K of ERA5. Its 850 hPa is 1–2 K colder, with RH 0.75–0.80 against
0.63, so the deep layer is regional. In October, ERA5's 850 hPa θ is 302.0 K
(RH 0.18), and thirteen's day 1643 is 295.7 K: the September inversion is
too high as well.

The 0.6–0.8 in the task is Klein and Hartmann's for SON or the annual mean.
For the arms' late March and early April, the measure is the CloudSat–CALIPSO
MAM 0.57. The defaults meet it (0.58 in `tradeCover.mjs`'s SE Pacific box
over days 1835–1845; 0.61–0.75 in the Peruvian box) with the wrong cloud:
resolved cumulus layers covered at about one half under a 2 km inversion,
with 45–98 g/m² of low water.

## The inversion's budget

Moistening is (L/c_p) dq_t/dt in K/day, a day's mean over 256 CPU steps
(`tropicalHeating.mjs`, `BOXES` the Peruvian and 20S boxes). The rows are the
two layers at the top of the moist layer, below the resolved inversion.
"dyn" is the dynamics' tendency (subsidence and advection), "BL" the
boundary layer's mixing, "shallow" and "deep" the plume's two branches; rain
is the box's surface rain.

| state, box | layers | dyn | BL | shallow | deep | rain, mm/d |
|---|---|---|---|---|---|---|
| day 1643 (SON), Peruvian | 841, 866 hPa | −9.8, −7.6 | +9.9, +11.6 | +2.0, +3.5 | 0 | 0.00 |
| day 1643, 20S | 842, 867 hPa | −6.3, −8.2 | +7.9, +9.5 | +0.5, +0.7 | 0 | 0.00 |
| defaults day 1835, Peruvian | 838, 863 hPa | −12.0, −5.0 | +7.9, +1.5 | +7.0, +8.7 | −0.1, −0.3 | 0.06 |
| defaults day 1835, 20S | 840, 865 hPa | −12.9, −4.6 | +16.1, +10.3 | +0.8, +0.6 | 0 | 0.00 |
| defaults day 1845, Peruvian | 811, 837 hPa | −6.7, −5.1 | +1.5, +1.3 | +7.9, +7.6 | +0.1, −0.9 | 0.64 |
| defaults day 1845, 20S | 838, 863 hPa | −0.5, −0.3 | −0.5, −0.1 | +3.7, +4.7 | −2.1, −3.3 | 1.08 |

**Where the deck runs, the boundary layer holds the top.** In September
its mixing over the deck's depth puts back 100–150 percent of the dynamics'
drying, and the plume adds a tenth to a third as much again.

**Where the gate has shut, the shallow plume holds it.** On day 1845, the
plume puts back 115–150 percent of the drying 1.6–1.9 km up, and the
boundary layer a quarter. The box's longwave cooling over its three top
layers falls from about 85 W/m² in September to 45, and the boundary layer's
own cloud-top ΔF from 52 to 24 W/m².

**The deep branch rains** 0.64 and 1.08 mm/d from tops near 3 km. M22 asked
for SE Pacific rain below 0.3 mm/d.

The transition is visible within the arm. At 20S, day 1835's top is still
the boundary layer's (16 K/day). By day 1845 the subsidence there has
weakened for a few days, and the plumes have taken the top to 2.7 km.

## Attribution

**The plume holds the inversion up.**

- On the gate-shut branch, the shallow plume is the only process that
  carries the boundary layer's water to 1.6–2.3 km against 5–7 K/day of
  drying.
- It fires at the full Grant (2001) closure, 0.03 w*, giving 21–28 g/m²/s at
  base, with almost no inhibition. The resolved subcloud layer is within
  0.4 K of well-mixed in θ_v, while the cloud layer above it is a
  trade-cumulus layer: θ_l rises about 3 K and q_t falls about 3.5 g/kg per
  kilometer.
- Its tops sit at the resolved inversion (2.31 km against 2.32 on day 1845).

**The bistable switch keeps it there.** The deck's veto is the only thing that
stops the plume, and the gate tests the jump at a height (the boundary-layer
top) that the plume branch has left 1 km below the inversion. A column that
has lost its deck has no way back. The ten days of the defaults arm move the
Peruvian box from 0.28 of its columns running to 0.06, and the inversion
from 1.83 to 2.32 km.

**Not the cause:**

- **Subsidence.** At the model's inversion it is at or above Earth's at
  Earth's inversion. The steady inversion is 0.5–1 km higher because the
  plume's entrainment matches the subsidence there.
- **SST.** The arms hold WOA23. In the coupled run, a free troposphere 6 K
  cold at 700 hPa (the tropical cold bias, `docs/toa-radiation-bias.md`)
  lowers the March EIS to 2.9 K. That deepens the coupled run's collapse but
  not the arms'.
- **The gate's thresholds.** M22 found that a deck run at the plume-lifted
  inversion grows 800–900 m thick with 400–900 g/m² of water, which is why the
  cumulus ceiling stands it down. The thresholds see the plume's inversion
  correctly.

## Changes, ranked

**1. The IFS's stratocumulus switch on the plumes.**

The IFS separates a well-mixed, stratocumulus-topped boundary layer (eddy
diffusivity and mass flux, no shallow convection) from a decoupled cumulus
layer (shallow convection allowed) with "an EIS threshold of 7 Kelvin"
(Bechtold et al. 2014, ECMWF Technical Memorandum 725, §2.1, after Köhler et
al. 2011 and Wood and Bretherton 2006). Here that is the opening the deck's
veto already applies: 0 for both branches where the column's EIS (the
radiation's 'eis': the σ = 0.7 layer against the lowest, moist-adiabatic from
the lowest layer's LCL) exceeds 7 K. It is an option, off by default, on both
engines. The GPU already computes the EIS inline (`stratiformShare`), and the
CPU has `inversionStrength`.

Expected effect:

- **The SON state** (EIS 11.3): it removes the 2.0–3.5 K/day of shallow
  moistening in the 0.27 of columns whose gate is shut.
- **The MAM arms** (medians 4.0–7.4 K; ERA5 7.3 K): in the columns above
  7 K, the boundary layer's 1.3–1.5 K/day cannot hold a top at 1.6–1.9 km
  against 5–7 K/day of drying. The top falls about a layer a day until the
  jump at h passes 4 K, the gate opens on its one-day memory, and the deck's
  veto holds it, as in September. The near-surface air then cools, which
  raises the EIS and spreads the switch.
- **Elsewhere:** the trade boxes (EIS below 3 K) are untouched; Namibia
  (2.4 K in the arm) moves little; California (10.2 K) and the subtropical
  decks gain; the cold tongue may cross 7 K in SON, so the N=128 stress
  protocol applies.

Caveats: the IFS's authors call the fixed threshold unsatisfactory, because
EIS differs between regions (the NE Atlantic's decks form at lower EIS than
the SE Pacific's) and because the switch is not smooth. The threshold is
theirs, not a value fitted here.

**2. The deck's water at its top.**

The deck's water sits in the layer nearest σ = 0.92 (about 0.7 km) and never
enters q_c. Its longwave cooling therefore lands 0.5–1 km below the deck's
top, outside the cloudy run that the boundary layer's cloud-top closure
sums. Earth's 70–100 W/m² of divergence sits at cloud top (Bretherton et al.
2010; Wood 2012).

Putting the water in the layer that holds h would sharpen the jump the gate
tests where the deck runs, and would give the closure its cooling. It acts
only where the deck runs, so it lengthens the deck branch but does not
reopen it.

**3. A decoupled slab.**

The deck averages θ_l and q_t over everything below h. A decoupled column's
moist subcloud air therefore lowers its cloud base to 1.13 km and gives it
300 g/m² of water. Earth's decoupled decks keep their cloud above the
decoupling (Wood and Bretherton 2004; Jones et al. 2011).

A slab from the decoupling height z_b to h would bring the water toward
40–150 g/m². The radiation's 150 g/m² cap limits the shortwave change, so
this is consistency more than cover.

## The switch in the arms

Change 1 is the moist option `stratocumulusInversion` (K; null, the default,
is none) on both engines. The arms set it to the IFS's 7 K. They are N=64
fixed-SST runs over WOA23 on bl42, 20 days each, with states saved on days
10, 15 and 20, on the protocol of `docs/trade-cumulus-cover.md`:

| arm | start | options |
|---|---|---|
| defaults (control) | thirteen64 day 1825 | – |
| switch | day 1825 | `MOIST='{"stratocumulusInversion":7}'` |
| pair (control) | day 1825 | `RADIATION='{"varianceCover":"water","buoyancyLength":0.53}'` |
| switch with the pair | day 1825 | both |
| September (control) | thirteen64 day 1643 | – |
| September with the switch | day 1643 | `MOIST` as above |

Checks:

- **Parity.** The plume parity harness (N=6, 362 random columns) matches
  the engines with the switch at 2 K. Few of its columns reach 7 K. The test
  files `test/convection.test.mjs` (47 tests), `test/boundaryLayer.test.mjs`
  (16) and `test/gpuModel.test.mjs` (32) pass on the GPU, one at a time, and
  the CPU files pass, pins included.
- **Engine comparison.**
  `MOIST='{"stratocumulusInversion":7}' scripts/compareEngines.mjs` on
  eleven128_day1825 (N=128) fails on the same 62 regions as at the defaults,
  every one of them over 1e-3 there already. The cumulus base flux differs between the engines by 3.6e-2 relative rms after
  one step (3.3e-2 at the defaults) and by 4.8e-2 after two (4.5e-2).

**The switch holds a column on the deck side while its EIS stays above 7 K.
It does not bring back one that has fallen below.**

Peruvian box (`deckGate.mjs`, four CPU steps from each state):

| arm, day | runs | resolved inversion, km | EIS, K | above 7 K | plume on | low cover |
|---|---|---|---|---|---|---|
| defaults, 1835 | 0.28 | 1.83 | 6.1 | | 0.56 | 0.75 |
| switch, 1835 | 0.36 | 1.77 | 6.1 | 0.45 | 0.43 | 0.79 |
| defaults, 1840 | 0.05 | 2.01 | 5.8 | | 0.72 | 0.68 |
| switch, 1840 | 0.05 | 2.21 | 4.5 | 0.11 | 0.70 | 0.58 |
| defaults, 1845 | 0.06 | 2.32 | 4.0 | | 0.71 | 0.61 |
| switch, 1845 | 0.18 | 1.95 | 6.4 | 0.43 | 0.36 | 0.76 |
| pair, 1835 | 0.61 | 1.87 | 7.0 | | 0.30 | 0.71 |
| switch with the pair, 1835 | 0.55 | 1.73 | 6.3 | 0.39 | 0.27 | 0.81 |
| pair, 1840 | 0.49 | 1.84 | 5.5 | | 0.40 | 0.53 |
| switch with the pair, 1840 | 0.20 | 1.91 | 5.0 | 0.19 | 0.53 | 0.72 |
| pair, 1845 | 0.18 | 1.95 | 5.6 | | 0.64 | 0.26 |
| switch with the pair, 1845 | 0.20 | 2.11 | 5.4 | 0.05 | 0.60 | 0.52 |

"Above 7 K" is the share of column-steps whose radiation EIS exceeds the
switch's threshold. The controls' reports predate that row.

On the 20S 85–75W line, the EIS stays above the threshold. By day 1845 the
defaults arm has gone over to the plume branch, and the switch keeps it
from doing so:

| arm, day 1845 | runs | resolved inversion, km | EIS, K | plume on | low cover |
|---|---|---|---|---|---|
| defaults | 0.00 | 2.73 | 5.1 | 0.75 | 0.44 |
| switch | 0.44 | 1.98 | 8.8 | 0.02 | 0.84 |
| pair | 0.53 | 1.92 | 8.0 | 0.39 | 0.58 |
| switch with the pair | 0.13 | 2.40 | 6.7 | 0.59 | 0.43 |

In the Peruvian box, the mechanism expected under change 1 runs only part of
the way:

- On day 1845 the inversion is 0.4 km lower than in the defaults arm.
- The EIS is 6.4 K against 4.0, and half as many columns run the plume.
- The deck runs on 0.18 of the box against 0.06. The boundary layer holds the
  cloud as resolved cloud, but the gate stays shut on the jump at h.

The box's EIS sits around the threshold (means 4.5–6.4 K; ERA5 7.3 K), so a
synoptic dip starts the plume again:

- On day 1840 the box's EIS falls to 4.5 K, and only 0.11 of the box lies
  above 7 K.
- On the 20S line, with the pair, the EIS falls to 6.7 K by day 1845. The
  run with the switch has gone further over to the plume branch than the run
  without it (inversion 2.40 km against 1.92). By then the two runs have
  parted synoptically.

**SE Pacific cover and the global cost.** The table gives low cover from
`tradeCover.mjs`, one CPU step from each state, averaged over the three
states. The radiative effects and rain are global means over days 6–20.
Fluxes are in W/m² and rain in mm/d.

| arm | SE Pacific (cover/water, g/m²) | N Pacific trades | S Pacific trades | Atlantic trades | sea 20–40S | sea 20–40N | SWCRE | LWCRE | net CRE | rain |
|---|---|---|---|---|---|---|---|---|---|---|
| defaults | 0.583/50 | 0.515 | 0.522 | 0.513 | 0.549 | 0.527 | −60.5 | 33.7 | −26.7 | 3.12 |
| switch | 0.577/52 | 0.524 | 0.537 | 0.550 | 0.553 | 0.549 | −61.7 | 33.6 | −28.2 | 3.05 |
| pair | 0.302/44 | 0.269 | 0.349 | 0.223 | 0.340 | 0.354 | −50.6 | 32.1 | −18.6 | 3.02 |
| switch with the pair | 0.498/65 | 0.438 | 0.317 | 0.453 | 0.439 | 0.492 | −56.8 | 32.8 | −24.0 | 3.06 |
| September | 0.665/77 | 0.516 | 0.473 | 0.423 | 0.592 | 0.461 | −61.1 | 31.3 | −29.9 | 3.02 |
| September with the switch | 0.724/100 | 0.539 | 0.492 | 0.432 | 0.601 | 0.461 | −61.9 | 31.0 | −30.9 | 3.00 |

The observed values are:

- SE Pacific: MAM 0.57 and SON 0.70 (CloudSat–CALIPSO, Muhlbauer et al.
  2014, over 10–30S 75–95W).
- Trades: low cover 0.2–0.4 (MODIS, CALIPSO-GOCCP).

Against those:

- **Defaults.** With the switch, the SE Pacific box keeps its MAM cover
  (0.577 against 0.57) and its water, at −1.2 W/m² of SWCRE.
- **Pair.** The switch brings the box from 0.30 to 0.50, at −6.2 W/m² of
  SWCRE.

The switch acts on every column above 7 K, not on the decks alone:

- Between 60S and 60N, 0.16–0.22 of the sea lies above 7 K. In the
  controls, plumes run over 0.07–0.09 of the sea's area within it, and with
  the switch over about 0.01.
- With the pair, the switch raises the low cover:
  - 20–40N: 0.35 → 0.49;
  - N Pacific trades: 0.27 → 0.44;
  - Atlantic trades: 0.22 → 0.45.

  These rises undo the pair's trade cover and pass the trades' 0.4.
- The N Pacific box has 0.23–0.27 of its columns above 7 K, in its
  northeastern corner toward the California deck.
- The Atlantic box has 0.01, so its cover arrives from upstream. In the
  Canary deck's region (20–30N 30–15W), the switch carries the EIS from 6.1
  to 8.6 K by day 1840. There 0.93 of the region lies above 7 K against 0.21,
  and the plume runs on 0.09 against 0.79: the mechanism expected under
  change 1 runs its full course.

**September.** Over WOA23, the September control holds the deck branch on
its own:

- 0.89–0.90 of the Peruvian box runs on days 1653–1663.
- The EIS is 9.7–12.0 K, the inversion 1.25–1.53 km and the low cover
  0.80–0.96.

The September start therefore does not test whether the switch holds the
deck. What the switch does there:

- It removes the plume left in the Peruvian box (0.10–0.16 → 0.00–0.01).
- It raises the SE Pacific box's cover from 0.665 to 0.724 (SON 0.70) and
  its water from 77 to 100 g/m².
- The cost is −0.8 W/m² of SWCRE.

**No direct reach into the cold tongue.** At most 0.01 of 5S–5N 150–80W
lies above 7 K in any arm (the September control's last state; 0.003 in
twelve_day1597 at N=128), and with the switch at most 0.001. In the September pair, the time-mean equatorial stress
(`armMetrics.mjs`, three states) still moves by −5, −10 and +2 mN/m² in
160–140W, 140–120W and 120–100W (−112, −81, −49 → −117, −91, −47).
**N=128.** The stress protocol runs twelve_day1597 over WOA23 on bl42 and
takes the means of days 1609–1617. Each arm ran three times: unperturbed,
and as two replicates with every θ scaled by 1 + 10⁻⁴u, u uniform in ±1/2
(`PERTURB=1e-4`, `SEED` 1 and 2). The bins give the stress in mN/m². SLP
is 180–160W less 120–100W, and SWCRE is the global mean over days 6–20.

| arm | 160–140W | 140–120W | 120–100W | SLP, hPa | SWCRE, W/m² |
|---|---|---|---|---|---|
| defaults, unperturbed | −142.7 | −100.2 | −63.1 | −5.34 | −60.3 |
| defaults, seed 1 | −135.5 | −93.6 | −61.0 | −4.94 | −60.4 |
| defaults, seed 2 | −135.7 | −93.5 | −60.6 | −5.07 | −60.5 |
| switch, unperturbed | −128.7 | −90.0 | −61.0 | −4.98 | −61.6 |
| switch, seed 1 | −138.1 | −95.8 | −61.3 | −5.08 | −61.5 |
| switch, seed 2 | −121.6 | −90.0 | −64.5 | −4.50 | −60.9 |

- **The defaults are reproducible.**
  - The two seeds agree within 0.4 mN/m².
  - The unperturbed run at the current commit reproduces the one made before
    it to 0.1 mN/m².
  - The unperturbed run sits 7, 7 and 2 mN/m² stronger than the seeds, so the
    perturbation itself shifts the means. The seeds measure the spread among
    replicates, not their distance from the unperturbed run.
- **The switch makes the stress depend on the realization.**
  - Its two seeds differ by 16.5 mN/m² in 160–140W.
  - By day 1617 it shuts the plume on 0.55–0.65 of the SE Pacific box, and
    the seeds part there: EIS 9.5 against 8.1 K, plume on 0.21 against 0.32
    of the box.
- **Its equatorial effect is not established.** Paired realization by
  realization, the switch moves:
  - 160–140W by +14.0, −2.6 and +14.1 mN/m² (mean +8.5, a weakening of 6
    percent);
  - 140–120W by +10.2, −2.2 and +3.5;
  - 120–100W by +2.1, −0.3 and −3.9.

  The weakening seen without seeds is the top of that spread. The mean is
  about the size of the switch's own spread, so three pairs do not separate
  it from zero.
- **No direct reach.** No column of 5S–5N 150–80W with EIS above 7 K runs a
  plume in any switch state, so any effect there is remote.
- **Global SWCRE** moves by −1.3, −1.1 and −0.4 W/m².

**What the arms say about change 1.** At the IFS's 7 K, the switch is a
hold, not a restore:

- It keeps a stratocumulus-topped column on the deck side for as long as the
  column's EIS stays above 7 K. That is the Peruvian box in September, and
  in MAM the 20S line, California (low cover 0.95 against 0.66 on day 1845)
  and the Canary deck.
- It does not reopen the MAM deck. That column's EIS has already fallen to
  the threshold: the boundary layer is too deep and decoupled at 850 hPa
  (above), and the plume branch's own warming of the lower layer keeps it
  there.
- With the pair, its cover gain in the SE Pacific comes with a larger one in
  the subtropics of both hemispheres.

## What remains

- **The MAM column's EIS.** Earth's 7.3 K against the arms' 4.0–6.4 K is
  the deep, warm lower layer of "Against observation", and the switch cannot
  supply it. Changes 2 and 3 act where the deck runs: they lengthen its
  branch but do not reopen it.
- **A threshold or a ramp fitted here.** A lower threshold, or a smooth
  ramp over EIS, would act on more of the MAM box. It would also act on more
  of the trades and the 20–40° seas, which the pair arm shows are sensitive.
  The IFS's authors call the fixed threshold unsatisfactory. Such a value
  would be fitted here and is not in these arms.
- **More replicates, if the equatorial effect matters.** With the switch,
  the N=128 realizations of 160–140W span 16.5 mN/m². Separating a 6
  percent mean weakening from zero at that spread takes about eight
  replicates of each arm.
- **A coupled judgment.** The arms are fixed-SST runs of 20 days. Before
  any default changes, a coupled run has to judge the switch's global SWCRE
  cost: −1.2 W/m² at N=64 and −0.4 to −1.3 at N=128 under the defaults, −6.2
  under the pair.
