# The top-of-atmosphere radiation bias

Why spin-up thirteen's planet reflects too much sunlight and its clear sky
emits too little longwave, measured on 2026-10-19 on the N=64 run's
year 4–5 states, the N=128 run's first two years and N=64 fixed-SST arms.

## The question

Day means from the spin-up logs (every step of the day), against CERES
EBAF Ed4.0 (Loeb et al. 2018, Table 5: clear-sky fluxes from cloud-free
scenes) and the same with the clear sky of the whole region, the model's
definition (Loeb et al. 2020: clear-sky OLR −2.2, SW +0.5 W/m²):

| | Ts °C | ASR | OLR | albedo | SWCRE | LWCRE | clear-sky OLR |
|---|---|---|---|---|---|---|---|
| thirteen64 year 4 | 13.68 | 230.6 | 228.1 | 0.323 | −61.9 | +28.6 | 256.7 |
| thirteen64 year 5 | 13.55 | 231.0 | 227.9 | 0.321 | −61.3 | +28.5 | 256.4 |
| thirteen128 year 1 | 14.59 | 232.7 | 232.0 | 0.317 | −59.1 | +28.4 | 260.3 |
| thirteen128 year 2 | 14.16 | 233.2 | 231.1 | 0.315 | −59.6 | +27.8 | 258.9 |
| eleven128 year 1 | 14.86 | 237.7 | 234.1 | 0.302 | −53.9 | +27.0 | 261.1 |
| eleven128 year 5 | 15.86 | 237.9 | 236.8 | 0.301 | −54.1 | +26.6 | 263.4 |
| EBAF Ed4.0, clear scenes | 14.85 | 240.9 | 240.1 | 0.291 | −45.8 | +27.9 | 268.1 |
| EBAF, whole-region clear sky | | | | | about −45.3 | about +25.8 | about 265.9 |

Against the whole-region clear sky the cloud reflects 16 W/m² too much
and the clear sky emits 9.5 W/m² too little. The clear sky's reflectance,
0.148, is a little darker than EBAF's 0.157 (about 3 W/m²), and the
atmosphere absorbs 80 W/m² against about 79. The errors cancel in the
net; the surface, short of sunlight (about 150 W/m² absorbed against
163.8, Kato et al. 2018), runs 1.3 K cold at N=64, its open sea 2.3 K
below WOA23 on day 1825.

eleven (bl36, 8-pass mesh, the isopycnal ocean, the defaults before
Oct 19) was already 8 W/m² too bright in SWCRE; thirteen adds about 5.

## Instruments

- `scripts/independentColumns.mjs STATE`: each column of a
  saved state after one CPU step, its cloud read as the radiation reads
  it (layer covers from `radiation.layerCover`, the visibility rule, the
  plumes' cumulus, the deck), and its day-mean reflected sunlight at 24
  hour angles in one band with the radiation's own two-stream. Two
  treatments of the same cloud: the model's, every layer's water pooled
  over the overlapped total cover f̄, and 128 independent columns per
  cell drawn from the same exponential-random overlap, a Markov chain on
  each layer's cloud indicator whose pair covers are those of
  `overlapped()`, each layer's water in its own cover; `FSD` adds gamma
  distributed in-cloud water in each cloudy layer of each column.
- `scripts/dayRadiation.mjs STATE` with `RADIATION`: the model's own CPU
  radiation on one state at 24 hours of its day, day means of ASR, OLR
  and the cloud effects.
- `scripts/fixedSst.mjs`, which takes a diagnostics frame every day,
  prints its day means (Ts, ASR, OLR, albedo, SWCRE, LWCRE, clear-sky
  OLR, rain) and saves the cloud effects: N=64 arms from
  thirteen64_day1825 over WOA23's sea surface, 20 days, judged on days
  6–20. PERTURB and SEED make replicates of an arm.
- `scripts/radiationBenchmark.mjs` with `LEVELS=bl42`: the clear-sky gas
  optics on the standard atmospheres against RRTMG and the ICRCCM
  line-by-line fluxes.
- `scripts/cloudRegimes.mjs` in the quarterly reports, and the saved day
  means of the year-5 quarterly states by latitude band.
- `scripts/armMetrics.mjs` (`docs/equatorial-wind-bias.md`): an arm's
  equatorial stress and sea-level pressure by 20° bin over its saved
  states; `scripts/tradeCover.mjs` (`docs/trade-cumulus-cover.md`): low
  cover, low water, the mixing top and the deck's gate in the regime
  boxes.

## Findings

### 1. The shortwave pools every layer's water over the column's total cover

Joshua's hypothesis, refined. The radiation does not spread a cell's
cloud over the whole cell: the shortwave blends a clear column with a
cloudy one at the overlapped cover f̄ (radiation.module.js, `column()`).
But the cloudy column holds the sum of every layer's water spread
uniformly over f̄, so a layer whose own cover f_k is smaller than f̄ has
its water spread over f̄ instead of f_k. A low cloud of in-cloud optical
depth 20 covering 0.3 of a cell under a cirrus of optical depth 1
covering 0.8 (random overlap, f̄ = 0.86) becomes one cloud of optical
depth 7.9 over 0.86 of the cell: at μ = 0.5 over a black surface it
reflects 0.45 of the beam, where the two layers as they lie reflect
0.29. Reflectance is concave in optical depth, so spreading a given water
path over more area brightens it. The longwave does not make this
error: it carries each layer's cloud in its own cover (the two-region
transfer under the same overlap). The plumes' cumulus water, about
2.5 g/m² of the grid mean but over 1000 g/m² inside its 0.2–0.4 percent
of the cell, is pooled the same way.

Independent columns of the same overlap, each layer's water in its own
cover, reflect 7.6 to 7.9 percent less (one band; day 1825 83.7 against
77.1 W/m², day 1643 93.0 against 85.9). Columns with cloud both below
and above 680 hPa (46 percent of the globe) carry 70 percent of it,
low-only columns 26 percent, and every latitude band shares it. In the
model's own radiation on day 1825 (24 hourly calls) the change is
**+5.0 W/m² of absorbed sunlight**, all of it in SWCRE (−56.0 to −51.0),
the longwave untouched.

### 2. The water inside each cloud is uniform

Real cloud water varies inside a cloud, and the independent-column albedo
of a distribution of optical depths lies below the albedo of its mean.
MODIS gives the reduction χ = exp⟨ln τ⟩/⟨τ⟩ at 0.65–0.8 in hemispheric
and global monthly means (Oreopoulos and Cahalan 2005; Cahalan et al.
1994 found about 0.7 in marine stratocumulus), which for a gamma
distribution is a fractional standard deviation (FSD) of the water of
0.6 to 0.85. Shonk et al. (2010) give FSD 0.75 ± 0.18 from their review
of the observations, Hill et al. (2012) 0.60 as CloudSat's mean for ice
water at its scale; the IFS used 1. At FSD 0.75 the model's own
radiation on day 1825 gains **+2.0 W/m² more** absorbed sunlight and its
cloud longwave effect falls by **1.0 W/m²** (OLR +1.0), from the thinner
emissivity of an uneven cloud at the same mean water; the net is +6.0
W/m² with finding 1.

| | ASR | SWCRE | LWCRE | net change |
|---|---|---|---|---|
| pooled (today) | 231.5 | −56.0 | +28.6 | |
| independent columns | 236.5 | −51.0 | +28.6 | +5.0 |
| independent columns, FSD 0.75 | 238.5 | −49.0 | +27.6 | +6.0 |

(model CPU radiation, thirteen64 day 1825 held fixed, day means of 24
calls)

### 3. Where the excess sits: the subtropics and tropics

Year-5 SWCRE by latitude (four saved day means, one per season), with
CERES EBAF's annual zonal means as recalled (approximate, to be replaced
by the EBAF file):

| band | area | SWCRE | EBAF (approx.) | excess × area |
|---|---|---|---|---|
| 60–90S | 0.067 | −27.2 | about −20 | −0.5 |
| 40–60S | 0.111 | −71.7 | about −60 | −1.3 |
| 20–40S | 0.150 | −67.1 | about −40 | −4.0 |
| 0–20S | 0.172 | −71.6 | about −45 | −4.6 |
| 0–20N | 0.171 | −73.1 | about −55 | −3.1 |
| 20–40N | 0.150 | −56.0 | about −40 | −2.4 |
| 40–60N | 0.111 | −55.0 | about −55 | 0 |
| 60–90N | 0.067 | −37.7 | about −35 | −0.2 |

The excess lies between 40S and 40N, most of it in the southern tropics
and subtropics, where the regime tables (day 730) show the trade-cumulus
regimes at twice Earth's cover:

| regime (day 730, N=64) | total cover | low cover | low in-cloud τ | Earth |
|---|---|---|---|---|
| S Pacific trades 10–20S 160–120W | 0.58 | 0.51 | 7.1 | total 0.35–0.55, low 0.2–0.4 |
| N Pacific trades 15–25N 170–130W | 0.90 | 0.57 | 14.9 | as above |
| SE Pacific 10–30S 110–80W | 0.65 | 0.56 | 8.3 | low 0.6–0.8 |

Globally the total cover matches ISCCP (0.65) only because the middle
cloud is a third of Earth's (0.066 against about 0.2) while the low cloud
is generous (0.40). The high class is thick: in-cloud 224 g/m², τ 24,
0.25 of its cover above τ 23 and 0.39 below 3.6 against Earth's thin
share of about 0.6, and only 7 of its 64 g/m² colder than 235 K; by the
optics' phase ramp the rest is partly supercooled liquid, which reflects
2.6 times as much per gram as ice at −20 °C.

### 4. eleven to thirteen: the kinetic overshoot

N=64 arms from thirteen64_day1825 over WOA23's sea surface, means of days
6–20 (true OLR, the arm's clear-sky OLR less its LWCRE):

| arm | Ts °C | ASR | OLR | SWCRE | LWCRE | clear-sky OLR | rain mm/d |
|---|---|---|---|---|---|---|---|
| defaults | 15.11 | 229.7 | 226.1 | −60.5 | +33.7 | 259.8 | 3.12 |
| cloudTop.cumulusCapped false | 15.09 | 228.5 | 225.8 | −61.6 | +33.9 | 259.8 | 3.12 |
| **cumulusOvershoot 1** (the default before Oct 19) | 15.12 | **235.5** | 226.6 | **−54.7** | +33.0 | 259.7 | 3.08 |
| plumeMomentum false | 15.10 | 229.7 | 226.2 | −60.4 | +33.6 | 259.9 | 3.09 |
| the run's own sea surface | 13.09 | 231.5 | 226.1 | −58.3 | +29.1 | 255.1 | 2.62 |

cumulusOvershoot 'kinetic' brightens the cloud by **5.8 W/m²**, the whole
of the difference from eleven; the capped cumulus and the plume momentum
do not (the capped columns are 1.1 W/m² darker). The extra reflection is
subtropical: 9.5 W/m² at 20–40S, 8.1 at 20–40N, 5–6 in the tropics and
the southern midlatitudes. The kinetic share detrains most of a shallow
plume's condensate in the layer below its first non-buoyant layer (shares
0.2–0.77 east of the warm pool, 7cb2e0c) instead of in that layer, the
inversion's: the water stays in the moist top of the cumulus layer, where
the variance cover spreads it, rather than in the dry inversion layer,
where it evaporates. On the arms' last day (1845) the low class holds
25.6 g/m² of grid-mean water at cover 0.39 under the kinetic share and
20.1 g/m² at 0.35 under the full-layer deposit (single days are too
noisy to split the trade boxes). The share follows the plume's physics;
what has no observational support is the cover the detrained water then
takes (the trades' low cover at twice Earth's).

### 5. The changes in the arm protocol

The same N=64 arms (thirteen64_day1825, WOA23's sea surface, days 6–20)
with findings 1 and 2 on both engines, and the equatorial stress (2S–2N
sea cells within 3°, `scripts/armMetrics.mjs` over the arm's eight saved
states; a mean of eight states carries about ±12 percent of day-to-day
noise):

| arm | ASR | OLR | SWCRE | LWCRE | net | rain | τ_x 160–140W / 140–120W / 120–100W mN/m² |
|---|---|---|---|---|---|---|---|
| defaults | 229.7 | 226.1 | −60.5 | +33.7 | +3.6 | 3.12 | −66 / −56 / −45 |
| cloudColumns | 235.6 | 226.3 | −54.6 | +33.6 | +9.3 | 3.11 | −66 / −54 / −44 |
| cloudColumns, cloudSpread 0.75 | 238.4 | 229.1 | −51.8 | +30.9 | +9.2 | 3.18 | −75 / −59 / −44 |
| the same and cumulusOvershoot 1 | 244.4 | 229.8 | −45.9 | +30.0 | +14.6 | 3.12 | −58 / −44 / −43 |
| cumulusOvershoot 1 alone | 235.5 | 226.7 | −54.7 | +33.0 | +8.8 | 3.08 | −71 / −51 / −50 |

The independent columns add 5.9 W/m² of absorbed sunlight with the
cloud free to adjust (5.0 on the fixed state), the spread 2.8 more and
2.8 less cloud longwave effect, and the three together put SWCRE on
EBAF's −45.8. The rain moves by at most 0.06 mm/d and the stress by no
more than the noise. The cloud itself barely moves: on day 1845 the
regime table gives the low class cover 0.375 against 0.391 and the high
0.276 against 0.271, with the same water and in-cloud τ; what changes
is how the same cloud reflects. Over a sea surface held at WOA23 the net top of the
atmosphere is then +14.6 W/m²: with the clear sky 6 W/m² short and the
cloud longwave effect 4 above the whole-region 25.8, the coupled planet
would warm until its emission caught up, as the question's arithmetic
said. The shortwave changes need the longwave's column (finding 7) to
balance.

Coupled, the same N=64 state continued 90 days (March to June) with its
own ocean and sea ice (`scripts/spinup.mjs` with FROM and ICE_FROM=1), by
30-day means:

| | days | Ts °C | ASR | OLR | net | albedo | SWCRE | LWCRE | clear-sky OLR | rain |
|---|---|---|---|---|---|---|---|---|---|---|
| defaults | 1–30 | 13.17 | 231.0 | 225.5 | +5.6 | 0.322 | −58.8 | +29.5 | 254.9 | 2.59 |
| | 31–60 | 14.48 | 230.2 | 227.3 | +2.9 | 0.324 | −59.3 | +30.5 | 257.8 | 2.66 |
| | 61–90 | 15.69 | 230.6 | 230.5 | 0.0 | 0.323 | −61.3 | +30.2 | 260.8 | 2.69 |
| cloudColumns, cloudSpread 0.75 | 1–30 | 13.29 | 239.6 | 227.9 | +11.7 | 0.296 | −50.2 | +27.3 | 255.2 | 2.62 |
| | 31–60 | 14.72 | 239.0 | 230.3 | +8.7 | 0.298 | −50.4 | +28.3 | 258.5 | 2.72 |
| | 61–90 | 16.00 | 238.8 | 233.8 | +5.1 | 0.299 | −53.4 | +28.1 | 261.9 | 2.78 |

The albedo comes to 0.296–0.299 against EBAF's 0.291 and holds; the
planet's net gain (+5.1 to +11.7 W/m²) runs 5 to 6 W/m² above the
defaults' and its surface is 0.3 K warmer by the third month and still
warming. The
equator line's stress at day 90 (one ocean step, 160E–100W) reads −0.071
against −0.055 N/m²; single readings move by 30–50 percent from day to
day, so the stress needs the 20-day-mean protocol at N=128 before either
option becomes a default.

### 6. Ruled out or not a lever

- **Droplet radius**: the sea's 11.8 µm (Han et al. 1994, AVHRR at
  3.7 µm, a cloud-top retrieval) agrees with the in situ radii MODIS's
  3.7 µm band matches best, while its 2.1 µm band reads 13 percent above
  them (King et al. 2013). In an adiabatic cloud the vertically uniform
  radius that gives the right optical depth for a water path is 5/6 of
  the top's (τ = 9 W/(5 ρ r_top), Szczodrak et al. 2001), so 11.8 µm is
  not too small.
- **Ice radius**: Ou and Liou's fit gives 16 to 74 µm between −60 and
  −20 °C, within the range of the cloud-top retrievals.
- **The capped cumulus and the plume momentum** (finding 4).
- **The gas optics** (finding 7).

### 7. The clear-sky longwave gap

From 11.7 W/m² (year 5 at N=64, 256.4 against EBAF's 268.1):

| part | W/m² | how |
|---|---|---|
| definition | 2.2 | EBAF's clear sky is cloud-free scenes, drier than the cloudy columns the model's clear-sky pass also counts (Loeb et al. 2020, global mean) |
| the cold sea surface | 3.5 to 4.5 | the arm over the run's own sea surface (2.3 K below WOA23) emits 255.1, over WOA23 259.8 |
| gas optics | about 0 (±1) | bl42 against RRTMG: tropical +0.3, MLS +0.4, MLW −1.0, SAW −1.4 W/m²; against ICRCCM line-by-line +0.2, +0.4, −0.1, −1.8, −1.3 |
| the column | about 5 to 6 | over WOA23's sea surface (Ts 15.1 °C) the clear sky still emits 259.8 against the whole-region 265.9; ERA5 attributes it below |

The column's 5 to 6 W/m² is the largest part.
The surface part goes once the shortwave warms the surface.

A first look at the tropics: the mean column of the ice-free sea within
20° of the equator in the defaults arm (eight states, WOA23's sea
surface at 27.6 °C) through the scheme's own clear-sky longwave, against
RRTMG's tropical atmosphere (TROP, 26.9 °C) with TROP's trace gases and
ozone on both:

| column | clear-sky OLR, black surface |
|---|---|
| the model's | 283.9 |
| the model's temperature with TROP's relative humidity | 287.8 |
| TROP's temperature with the model's relative humidity | 284.2 |
| TROP | 287.8 |

Against TROP the whole tropical deficit, about 5 W/m² at equal surface
temperature, is the humidity: above 330 hPa the model's relative
humidity stays at 0.28–0.29 (about 0.45 over ice) up to 110 hPa, where
TROP falls from 0.24 to 0.06, so the model carries two to three times
TROP's vapor at 150–300 hPa; below 850 hPa it is also moister (0.78–0.84
against 0.71) and 1–1.5 K colder, and at 300–600 hPa 2–4 K warmer.

TROP is the wrong reference aloft. Against ERA5 at the arm's dates
(`docs/clear-sky-olr-deficit.md`) the model's tropical upper troposphere
is 0.7–0.8 times as moist at 150–250 hPa, the flat profile is ERA5's
own, and moving the humidity to ERA5's would lower the clear-sky OLR by
0.9 W/m². The column is cold, not moist. Of the 5.6 W/m² between the
model's eight states (260.3) and EBAF's whole-region 265.9, the arm's
dates take 0.6, the troposphere's temperature 1.3 (the trade boundary
layer through 850 hPa 1.6 K cold, the upper troposphere and the
extratropical tropopause, under a lowermost stratosphere 1.3–2.2 times
as moist as ERA5's), the stratosphere above 50 hPa 1.1 (3–11 K cold, its
shortwave heating 10–35 percent below RRTMG's at 14–170 hPa), the
relative humidity 0.2 and the gas optics about 0.3; 2.1 is unexplained.

### 8. The three together

Findings 1 and 2 with the variance pair of `docs/trade-cumulus-cover.md`
(`varianceCover` 'water' with `buoyancyLength` 0.53, change 3 below) in
the arms of finding 5; the pair's row is that doc's:

| arm | ASR | OLR | SWCRE | LWCRE | clear-sky OLR | net | rain | τ_x 160–140W / 140–120W / 120–100W mN/m² |
|---|---|---|---|---|---|---|---|---|
| defaults | 229.7 | 226.1 | −60.5 | +33.7 | 259.8 | +3.6 | 3.12 | −66 / −56 / −45 |
| cloudColumns, cloudSpread 0.75 | 238.4 | 229.1 | −51.8 | +30.9 | 260.0 | +9.2 | 3.18 | −75 / −59 / −44 |
| the variance pair | 234.2 | 227.4 | −55.9 | +32.5 | 260.0 | +6.8 | 3.02 | |
| all three | 243.0 | 230.5 | −47.2 | +29.7 | 260.2 | +12.5 | 3.05 | −68 / −49 / −38 |

The changes add: the two arms' departures from the defaults predict an
ASR of 242.9, SWCRE −47.2, LWCRE +29.7 and a net +12.4. With
cumulusOvershoot 'kinetic' kept, SWCRE comes within 1.4 W/m² of EBAF's
−45.8. The cloud longwave effect stays 3.9 above the whole-region 25.8
and the clear sky 5.7 below 265.9, so over WOA23's sea the planet gains
12.5 W/m². The equatorial stress (−51.5 mN/m² over the three bins
against −55.4) and the sea-level pressure of 180–160W less 120W–100W's
(−1.46 hPa against −1.29) move within the eight states' noise.

Low cover in the regime boxes (`scripts/tradeCover.mjs`, means of days
1835–1845), with the SE Pacific's low water and the share of it where
the deck's gate is open; the pair's row is the trade doc's, whose
defaults read within 0.02 of these:

| arm | N Pacific trades | S Pacific trades | Atlantic trades | SE Pacific | sea 20–40S | sea 20–40N | SE Pacific water g/m² | deck open |
|---|---|---|---|---|---|---|---|---|
| defaults | 0.533 | 0.531 | 0.526 | 0.569 | 0.552 | 0.523 | 49.2 | 0.06 |
| cloudColumns, cloudSpread 0.75 | 0.556 | 0.555 | 0.501 | 0.601 | 0.554 | 0.538 | 67.5 | 0.19 |
| the variance pair | 0.424 | 0.317 | 0.441 | 0.472 | 0.448 | 0.463 | | |
| all three | 0.381 | 0.297 | 0.409 | 0.377 | 0.435 | 0.463 | 41.7 | 0.07 |
| Earth | 0.2–0.4 | 0.2–0.4 | 0.2–0.4 | 0.6–0.8 | | | | |

In the trades the three put the low cover in Earth's range. In the SE
Pacific they do more than the sum of their parts: columns and spread
raise the box's low water from 49 to 68 g/m² and open the deck's gate on
0.19 of it; with the pair as well the water falls to 42 g/m², the gate
to 0.07 and the cover to 0.377, 0.19 below the defaults' where the pair
alone took 0.11 off. Klein and Hartmann's 0.6–0.8 is the SON or annual
cover; for late March and early April, the arms' season, the SE Pacific
deck study (`docs/se-pacific-deck.md` on claude/cranky-noether-1eea9a)
reads about 0.57 from CloudSat–CALIPSO, which the defaults match with the
wrong cloud: cumulus under an inversion 0.5–1 km too high, where the
shallow plume holds it once the deck's gate has shut. The pair thins
that cumulus; the deck that should replace it is that study's change
(the plumes off where EIS exceeds 7 K, as in the IFS), not these.

At N=128, by the protocol of `docs/equatorial-wind-bias.md`
(twelve_day1597 over WOA23 on bl42, the mean of the states of days
1609–1617, sea cells within 3° of the equator; on twelve's 8-pass mesh
the sub-grid orography comes from the raster, without form drag), with
the sea-level pressure of 180–160W less 120W–100W's. Every arm but
columns and spread is the mean of three realizations: the unperturbed
run (the trade doc's for the defaults and the pair) and two with every θ
of the start scaled by 1 + 10⁻⁴·u (`fixedSst.mjs` PERTURB, seeds 1 and
2):

| τ_x, mN/m² | 140E–160E | 160E–180 | 180–160W | 160W–140W | 140W–120W | 120W–100W | 100W–80W | central three | SLP difference, hPa |
|---|---|---|---|---|---|---|---|---|---|
| observed, about | −15 | −20 | −40 | −50 | −55 | −50 | −30 | −52 | |
| defaults | −6.9 | −38.5 | −113.6 | −138.0 | −95.8 | −61.6 | −14.5 | −98.5 (−102.0, −96.7, −96.6) | −5.12 |
| cloudColumns, cloudSpread 0.75, one run | −11.3 | −25.3 | −102.7 | −129.1 | −99.6 | −62.4 | −15.9 | −97.0 | −4.72 |
| the variance pair | −11.7 | −29.0 | −89.8 | −107.6 | −87.0 | −65.0 | −16.0 | −86.5 (−79.6, −91.7, −88.3) | −3.98 |
| all three | −9.2 | −24.1 | −83.4 | −112.4 | −85.1 | −58.9 | −13.5 | −85.5 (−87.1, −87.4, −82.0) | −4.13 |

Day means, days 6–20, of the unperturbed runs (the realizations agree to
0.2 W/m²):

| arm | ASR | OLR | SWCRE | LWCRE | clear-sky OLR | net | rain |
|---|---|---|---|---|---|---|---|
| defaults | 235.8 | 239.1 | −60.3 | +29.0 | 268.1 | −3.3 | 3.08 |
| cloudColumns, cloudSpread 0.75 | 245.0 | 241.8 | −51.0 | +26.6 | 268.4 | +3.2 | 3.09 |
| the variance pair | 239.6 | 239.6 | −56.4 | +28.5 | 268.1 | 0.0 | 3.01 |
| all three | 248.2 | 242.4 | −47.8 | +26.2 | 268.6 | +5.7 | 3.08 |

The realizations of one arm differ by a standard deviation of 4.4 mN/m²
in the central three bins' mean (pooled over six degrees of freedom) and
0.14 hPa in the gradient; within one run the five states differ by 6–13
mN/m². Against the defaults the pair takes 11.9 ± 3.6 mN/m² off the
stress, 12 percent, and 1.13 ± 0.11 hPa off the gradient, 22 percent;
the three together take 13.0 ± 3.6 and 0.99 ± 0.11. The optics take
nothing back: all three less the pair is +1.0 ± 3.6 mN/m² and −0.14 ±
0.11 hPa. Columns and spread alone leave the stress where the defaults
have it (+1.4 ± 5.1, one run) and may ease the gradient by 0.40 ± 0.16
hPa. Single runs cannot resolve these differences: the trade doc's
defaults and pair fall at opposite ends of their spreads, which made
the pair's gain look like 22 percent. With the pair, with or without
the optics, 160W–140W stays 2.2 times the observed and 120W–100W does
not move. The shortwave changes add at N=128 as at N=64: the parts
predict SWCRE −47.1 and the arm gives −47.8.

Coupled, the continuation of finding 5 with the three, by 30-day means:

| | days | Ts °C | ASR | OLR | net | albedo | SWCRE | LWCRE | clear-sky OLR | rain |
|---|---|---|---|---|---|---|---|---|---|---|
| all three | 1–30 | 13.32 | 243.2 | 229.4 | +13.8 | 0.286 | −46.8 | +26.2 | 255.6 | 2.59 |
| | 31–60 | 14.72 | 241.4 | 231.2 | +10.2 | 0.291 | −48.1 | +27.4 | 258.6 | 2.66 |
| | 61–90 | 16.07 | 241.7 | 235.0 | +6.7 | 0.290 | −50.5 | +27.4 | 262.5 | 2.76 |

After the first month the albedo is EBAF's 0.291 and SWCRE within 2–5
W/m² of its −45.8; over the model's own sea the cloud longwave effect is 26.2–27.4,
0.4–1.6 above the whole-region 25.8. The planet gains 6.7 to 8.2 W/m²
more than with the defaults (finding 5's table) and its surface is
0.4 K warmer by the third month and still warming: the clear-sky
column (4) is the larger part of what is left to balance it. The equator line's stress at
day 90 reads −0.049 N/m² against the defaults' −0.055, a single reading.

## Ranked changes

ΔASR and ΔOLR from the radiation on the fixed day-1825 state (findings 1
and 2) and, in brackets, from the 20-day arms over WOA23 (findings 5 and
8):

| | change | reference | ΔASR | ΔOLR | Δclear-sky OLR | status |
|---|---|---|---|---|---|---|
| 1 | Independent columns in the shortwave (`cloudColumns`): each layer's water in its own cover under the overlap the cover and the longwave already use; the column's τ' lognormal with the overlap's exact first two moments, at four Gauss-Hermite values | Cahalan et al. 1994 (independent columns); Pincus, Barker and Morcrette 2003; Räisänen et al. 2004 (the overlap's columns); Barker 1996 (a distribution of τ in the two-stream) | +5.0 (+5.9) | 0 (+0.2) | 0 | both engines (8bcecac); parity test, trace-cloud case included; off by default |
| 2 | Water spread inside each cloud (`cloudSpread` 0.75), in both bands | Shonk et al. 2010; Oreopoulos and Cahalan 2005; Hill et al. 2012 | +2.0 (+2.8) | +1.0 (+2.8) | 0 | both engines (8bcecac), with 1; off by default |
| 3 | The cover of the shallow plumes' water: the variance cover's spread (1.0–2.5 g/kg, four to eight times the condensation's half-width) covers every cloudy layer of the cumulus layer by about one half whatever its water (BOMEX: 0.43–0.46 against the LES's 0.03). The variance pair (`docs/trade-cumulus-cover.md`): `varianceCover` 'water', the Gaussian that holds the layer's water, with `buoyancyLength` 0.53, the length bounded in stable layers | Sommeria and Deardorff 1977; Bougeault 1981; Galperin et al. 1988; Siebesma et al. 2003 (BOMEX) | (+4.5) | (+1.3) | (+0.2) | both engines (6254bce); off by default |
| — | The buoyant shallow plume (`cumulusEntrainmentHeight` 1, `cumulusAscent` 'kinetic'): BOMEX's plume top and the same 4.6 W/m², but the trade stress unmoved (alone and with the pair, single runs within one realization of the defaults' and the pair's means), the SE Pacific's cover at 0.35 | Siebesma et al. 2003 | (+4.6) | (+0.5) | (0) | both engines (7f950f6); an option |
| 4 | The column's clear-sky emission against ERA5 (`docs/clear-sky-olr-deficit.md`, finding 7): cold, not moist. The stratosphere's near-infrared heating, k-distributions fitted to RRTMG_SW (`nearInfrared` 'bands'), +4.4 K at 50 hPa at equilibrium; the extratropical tropopause's excess vapor; the trade boundary layer's 1.6 K | Chou and Suarez 1999; Iacono et al. 2008; Mlawer et al. 1997 | | | +1.0 to 1.1; up to +1.0; up to +1.2 | the refit on both engines (87b5a21 on claude/silly-chaum-3dfd05), off by default; the others to trace |
| 5 | The clear sky's reflectance, 0.148 against 0.157: the aerosol depths (sea 0.07, land 0.12) against MODIS and MISR | | −3 | | | examine |
| 6 | Judge the clear sky against EBAF's whole-region fluxes (265.9, LWCRE 25.8), not the clear scenes' 268.1 | Loeb et al. 2020 | | | 2.2 by definition | done here |

Changes 1 and 2 remove 7 to 9 of the 16 W/m² cloud excess and add one
observed parameter, the FSD; the overlap and its decorrelation length
stay those the cover and the longwave already use. They change the
cloud's optics, not its amount, and bring the shortwave and the longwave
to the same independent-column view of the cloud. Change 3 gives the
cumulus layer BOMEX's cover and the trades Earth's, and takes 12
percent off the N=128 central-Pacific stress and 22 percent off its
pressure gradient (three realizations each). The three add (finding 8):
SWCRE −47.2 in the arm and the albedo 0.290 coupled, against EBAF's
−45.8 and 0.291, with the kinetic overshoot kept. They should not go
into a running spin-up, nor into a new one ahead of the longwave side:
over the observed sea surface the planet then gains 12.5 W/m², and
coupled 6.7 to 8.2 more than with the defaults, because the clear sky's
column (4) is not yet fixed and the cloud longwave effect stands 3.9
above the whole-region 25.8 over that sea surface (0.4–1.6 over the
model's own).

## Tests

- `test/gpuModel.test.mjs` (33 of 33), with a parity test of 1 and 2
  at N=6: the engines' layer heating agrees to 2.4·10⁻⁴ K/day of 26.5,
  the absorbed sunlight to 4·10⁻³ W/m² and the OLR to 5·10⁻⁴ W/m², the
  columns absorb 2.6 W/m² more than the pooled column and the spread
  6.5 more; its trace-cloud case (columns whose only cloud is 10⁻³⁰ kg/kg
  in one layer) stays finite on the GPU, where the moments formed as
  cover · Σ square / depth² gave NaN. The variance pair's parity test is
  `docs/trade-cumulus-cover.md`'s.
- `test/boundaryLayer.test.mjs` (16 of 16), `test/convection.test.mjs`
  (45 of 45) and `test/cloudEffect.test.mjs` (3 of 3) on the GPU; the
  CPU's radiation, column and moist files (cloudEffects, cloudHeights,
  column, infraredSlices, insolation, latitudeColumns, longwaveOverlap,
  radiationEvery, radiativeColumn, wetColumn, plumeColumn, plumePath,
  moist, moistConvective, moistParcel, radiativeConvective, defaultRun:
  130 of 130). `test/physics.test.mjs` (37 of 37) pins the CPU's
  twelve-step digests of the default radiation, which hold with every
  option off.
- `scripts/compareEngines.mjs` on thirteen64_day1825 (`RADIATION`, and
  `MOIST` for the buoyant plume's options). With every option off it
  flags 56 regions: 54 over the tool's 10⁻³ tolerance (mixing, rain,
  convection, snow, the held clear-sky shortwave), the sea ice's ocean
  flux, which the CPU leaves at zero, and a momentum scratch input it
  reads as non-finite. With 1 and 2 every radiation region differs
  between the engines by what it differs by with the options off
  (reflected 1.2·10⁻² against 1.1·10⁻², absorbed 4.8·10⁻³ both, OLR
  3·10⁻⁵ both). With 1, 2 and the variance pair it flags 59: the
  radiation as before (absorbed 4.7·10⁻³, OLR 1.3–2.0·10⁻⁴ against
  0.3–1.5·10⁻⁴); the boundary layer's cover-dependent regions (regime,
  radiative velocity, mixing) differ at the first step by 2.4–4.0·10⁻³,
  where the defaults differ by under 2·10⁻⁴, and at the second by
  0.8–1.0·10⁻², inside the defaults' 1.2–1.4·10⁻²; the orographic drag's
  blocking rate and wave (1.4–1.6·10⁻³, from 1.4–3.7·10⁻⁴) and at the
  second step the surface drag and heat exchange (1.1–1.2·10⁻³, from
  4.5–4.7·10⁻⁴) cross the tolerance, and the surface shortwave and the
  stress fall under it. The buoyant plume's options flag the same 59.

## What remains

- The longwave side, before 1–3 go into a spin-up: the clear-sky
  column's 5.7 W/m² (4: the near-infrared refit's 1.0–1.1 is on
  claude/silly-chaum-3dfd05, off by default; the extratropical
  tropopause's vapor and the trade boundary layer's 1.6 K are to trace;
  2.1 is unexplained) and the cloud longwave effect, 29.7 against the
  whole-region 25.8 with the three on.
- The SE Pacific's deck (`docs/se-pacific-deck.md`): the plumes off
  where EIS exceeds 7 K, on both engines, then these arms again.
- Replicates for every N=128 stress change under about 9 mN/m² (two
  realization standard deviations): columns and spread's, the buoyant
  plume's and the single-run arms of `docs/equatorial-wind-bias.md`.
- The buoyant plume's stress and SE Pacific effects
  (`docs/trade-cumulus-cover.md`) before it is more than an option.
- The clear sky's reflectance (5) against MODIS and MISR.
- Why compareEngines flags 54 regions at thirteen64_day1825 with the
  defaults, and the boundary layer's departures at the first step with
  the variance pair (2.4–4.0·10⁻³).
