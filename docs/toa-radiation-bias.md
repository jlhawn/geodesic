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
  6–20.
- `scripts/radiationBenchmark.mjs` with `LEVELS=bl42`: the clear-sky gas
  optics on the standard atmospheres against RRTMG and the ICRCCM
  line-by-line fluxes.
- `scripts/cloudRegimes.mjs` in the quarterly reports, and the saved day
  means of the year-5 quarterly states by latitude band.

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
| the column | about 5 to 6 | over WOA23's sea surface (Ts 15.1 °C) the clear sky still emits 259.8 against the whole-region 265.9: a moist or cold upper troposphere, unmeasured (no ERA5 or AIRS climatology in the repository yet) |

The column's 5 to 6 W/m² is the largest part and the one left to find.
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

The whole tropical deficit, about 5 W/m² at equal surface temperature,
is the humidity. Above 330 hPa the model's relative humidity stays at
0.28–0.29 (about 0.45 over ice) up to 110 hPa, where TROP falls from 0.24
to 0.06, so the model carries two to three times TROP's vapor at
150–300 hPa; below 850 hPa it is also moister (0.78–0.84 against 0.71)
and 1–1.5 K colder, and at 300–600 hPa 2–4 K warmer. TROP's upper
troposphere is drier than the AIRS climatology, so the size of the
model's excess needs the AIRS or ERA5 tropical mean. The flat 0.28 from
330 hPa to the tropopause, below the condensation's critical humidity
there (0.75 over ice), is the profile to explain: the deep plumes'
detrainment, the ice cloud's evaporation and fall, the resolved ascent.

## Ranked changes

ΔASR and ΔOLR from the radiation on the fixed day-1825 state (finding 1
and 2) and, in brackets, from the 20-day arms over WOA23 (finding 5):

| | change | reference | ΔASR | ΔOLR | Δclear-sky OLR | status |
|---|---|---|---|---|---|---|
| 1 | Independent columns in the shortwave (`cloudColumns`): each layer's water in its own cover under the overlap the cover and the longwave already use; the column's τ' lognormal with the overlap's exact first two moments, at four Gauss-Hermite values | Cahalan et al. 1994 (independent columns); Pincus, Barker and Morcrette 2003; Räisänen et al. 2004 (the overlap's columns); Barker 1996 (a distribution of τ in the two-stream) | +5.0 (+5.9) | 0 (+0.2) | 0 | both engines; parity test, trace-cloud case included; off by default |
| 2 | Water spread inside each cloud (`cloudSpread` 0.75), in both bands | Shonk et al. 2010; Oreopoulos and Cahalan 2005; Hill et al. 2012 | +2.0 (+2.8) | +1.0 (+2.8) | 0 | both engines, with 1; off by default |
| 3 | The cover of the water the shallow plumes detrain under the inversion, where cumulusOvershoot 'kinetic' puts it: the variance cover's scale (5, against Mellor and Yamada's B2^½ of 3.2) and the detrained cloud's lifetime there, against the trades' observed 0.2–0.4 low cover | Mellor and Yamada 1982; Siebesma et al. 2003 (BOMEX) | up to +5.8 (the overshoot's whole effect) | | | arms to run |
| 4 | The column's clear-sky emission: the upper troposphere's temperature and humidity against ERA5 and AIRS | | | | +5 to 6 | data to fetch |
| 5 | The clear sky's reflectance, 0.148 against 0.157: the aerosol depths (sea 0.07, land 0.12) against MODIS and MISR | | −3 | | | examine |
| 6 | Judge the clear sky against EBAF's whole-region fluxes (265.9, LWCRE 25.8), not the clear scenes' 268.1 | Loeb et al. 2020 | | | 2.2 by definition | done here |

Changes 1 and 2 remove 7 to 9 of the 16 W/m² cloud excess and add one
observed parameter, the FSD; the overlap and its decorrelation length
stay those the cover and the longwave already use. With 3 the arms put
SWCRE on EBAF's (−45.9). They should not go into a running spin-up
alone: over the observed sea surface the planet then gains 14.6 W/m²,
because the clear sky's column (4) and the cloud longwave effect, 30 to
34 against 25.8 over that sea surface, are not yet fixed. 1 and 2 are
the least-regret pair: they change the cloud's optics, not its amount,
and bring the shortwave and the longwave to the same independent-column
view of the cloud.

## Tests

- `test/gpuModel.test.mjs` (32 of 32), with a new parity test of 1 and 2
  at N=6: the engines' layer heating agrees to 2.4·10⁻⁴ K/day of 26.5,
  the absorbed sunlight to 4·10⁻³ W/m² and the OLR to 5·10⁻⁴ W/m², the
  columns absorb 2.6 W/m² more than the pooled column and the spread
  6.5 more; its trace-cloud case (columns whose only cloud is 10⁻³⁰ kg/kg
  in one layer) stays finite on the GPU, where the moments formed as
  cover · Σ square / depth² gave NaN.
- `test/boundaryLayer.test.mjs` (16 of 16) and the CPU radiation files
  (cloudEffects, column, longwaveOverlap, radiationEvery, cloudHeights,
  infraredSlices, insolation: 62 of 62). `test/physics.test.mjs` (37 of
  37) pins the CPU's twelve-step digests of the default radiation, which
  hold with both options off.
- `scripts/compareEngines.mjs` on thirteen64_day1825 with 1 and 2
  (`RADIATION`): every radiation region differs between the engines by
  what it differs by with the options off (reflected 1.2·10⁻² against
  1.1·10⁻², absorbed 4.8·10⁻³ both, OLR 3·10⁻⁵ both). The 54 regions over
  the tool's 10⁻³ tolerance at this state (mixing, rain, convection,
  snow, the held clear-sky shortwave) are flagged with the options off as
  well.

## What remains

- The central-Pacific stress under 1 and 2 by the 20-day-mean fixed-SST
  protocol at N=128 (twelve_day1597, as in equatorial-wind-bias.md).
- Arms for 3; the ERA5/AIRS profile comparison for 4, starting from the
  tropics' flat upper-tropospheric humidity.
- Why compareEngines flags 54 regions at thirteen64_day1825 with the
  defaults.
