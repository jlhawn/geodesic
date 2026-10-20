# The clear-sky OLR deficit at the observed sea surface

Where spin-up thirteen's clear-sky outgoing longwave falls short of CERES,
measured on 2026-10-19 on the N=64 fixed-SST arm of master fb9fff5's
defaults over WOA23's sea surface (thirteen64_day1825, days 6–20, the
eight saved states of days 1831–1845, late March to early April), against
ERA5's climatology for the same days of the year.

## The question

The arm's clear sky emits 259.8 W/m² (day means of every step; 260.3 on
the eight end-of-day states used below) against CERES EBAF's whole-region
clear sky of about 265.9 W/m² (Loeb et al. 2018, 2020: 268.1 for
cloud-free scenes less 2.2). The gas optics on bl42 match RRTMG within
+0.3 to +0.4 W/m² on the warm standard atmospheres and −1.0 to −1.4 on
the winter ones, about −0.3 globally at the equinox. Against RRTMG's
tropical atmosphere (TROP) the model's tropical column is 3.9 W/m² short,
all of it humidity: its upper troposphere holds two to three times TROP's
vapor. TROP is dry aloft, so the question was whether the model's upper
troposphere is too moist against an observed climatology.

## Instruments

- ERA5 (Hersbach et al. 2020), the WeatherBench 2 climatology (Rasp et
  al. 2024) `1990-2019_6h_128x64_equiangular_with_poles_conservative.zarr`
  on Google's public bucket: temperature and specific humidity on 13
  levels from 1000 to 50 hPa, surface pressure, 2 m temperature and column
  vapor by day of year and hour, reduced to the arm's days (85–99) and to
  monthly means. 925 MB, Blosc/LZ4 Zarr, read by a 60-line decoder.
- Profiles: each state's T and q interpolated in ln p to ERA5's levels
  cell by cell, ERA5 bilinearly to the same cells, both masked below
  either surface; band means of T and q (arithmetic), relative humidity
  from the band-mean T and q by the model's own saturation formulas.
- Swaps: the model's own CPU radiation (`radiation.apply`, clear-sky
  pass) on the eight states with each cell's time-mean T shifted by
  ERA5 − model and its q scaled by ERA5/model (arithmetic means over the
  eight states) at ERA5's levels, interpolated in ln p onto the model's
  layers between 50 hPa and the surface. Each state keeps its own day-to-day
  departures, so the swap moves the mean state only.
- The scripts are archived outside the repository with the side sessions' notes (`swapOlr.mjs`,
  `profiles.mjs`, `stratProfile.mjs`, `zarr.mjs`, `reduceEra5.mjs`; for finding 6
  `heatBudget.mjs`, `swOzoneBands.mjs`, `rrtmgRun.mjs`, `lwTop.mjs`, `lwDecomp.mjs`, `lwRows.mjs`).

## Findings

### 1. The tropical upper troposphere is drier than ERA5, not moister

Band means at ERA5's levels (RH over ice below 0 °C, from band-mean T
and q):

| 0–20N | T model | T ERA5 | q model/ERA5 | RH_ice model | RH_ice ERA5 |
|---|---|---|---|---|---|
| 100 hPa | 191.1 | 192.8 | 0.68 | 0.63 | 0.70 |
| 150 | 204.6 | 205.8 | 0.72 | 0.52 | 0.60 |
| 200 | 218.7 | 219.8 | 0.71 | 0.45 | 0.55 |
| 250 | 230.9 | 231.8 | 0.73 | 0.38 | 0.47 |
| 300 | 240.9 | 241.7 | 0.82 | 0.34 | 0.39 |
| 500 | 266.7 | 267.5 | 0.89 | 0.33 | 0.35 |
| 700 | 281.6 | 282.9 | 1.03 | 0.52 | 0.46 |
| 850 | 289.9 | 291.6 | 1.09 | 0.71 | 0.59 |
| 1000 | 298.0 | 298.2 | 1.00 | 0.80 | 0.79 |

20S–0 is similar (q ratio 0.72–0.83 at 150–250 hPa, RH_ice 0.49–0.57
against 0.54–0.67). The flat RH of 0.28 over liquid from 330 to 110 hPa is
ERA5's own shape (0.28–0.38 over liquid); TROP, not the model, is the
outlier. Column vapor matches (39.2 against 38.0 kg/m² at 0–20N, 43.4
against 44.0 at 20S–0). The model is moister than ERA5 at 850 hPa
everywhere equatorward of 60° (RH +0.07 to +0.16) and through the
troposphere at 40–60° and over the Antarctic (q 1.06–1.23 times ERA5 at
400–850 hPa; the Arctic is drier, 0.78–0.94, and 2 K colder), and 1.1–2.2
times ERA5 at 100–200 hPa poleward of 40° (2.0–2.2 at 150–200 hPa).

### 2. Attribution: temperature, not humidity

Clear-sky OLR of the eight states, W/m², global and by band, with ERA5's
mean T and/or q swapped in:

| | global | 90–60S | 60–40S | 40–20S | 20S–0 | 0–20N | 20–40N | 40–60N | 60–90N |
|---|---|---|---|---|---|---|---|---|---|
| model | 260.26 | 199.8 | 245.8 | 277.6 | 281.6 | 287.1 | 273.2 | 234.3 | 197.6 |
| T | +2.45 | −0.2 | +0.4 | +0.9 | +3.7 | +4.1 | +2.0 | +2.0 | +6.6 |
| T at fixed RH | +1.27 | 0.0 | +0.6 | +0.7 | +1.4 | +1.5 | +1.1 | +1.4 | +4.2 |
| q | −0.91 | +0.8 | +1.4 | +0.5 | −3.3 | −2.5 | −1.3 | +0.7 | −1.3 |
| T and q | +1.51 | +0.6 | +1.9 | +1.3 | +0.3 | +1.4 | +0.6 | +2.7 | +5.6 |
| T, 600 hPa–surface | +1.20 | −1.0 | −0.1 | +1.2 | +2.1 | +2.3 | +1.4 | +0.3 | +1.5 |
| T, 300–600 hPa | +0.24 | −0.5 | −0.6 | −0.7 | +1.0 | +1.1 | −0.1 | −0.1 | +1.8 |
| T, 50–300 hPa | +1.01 | +1.4 | +1.0 | +0.5 | +0.6 | +0.7 | +0.7 | +1.7 | +3.4 |
| q, 50–300 hPa | −0.52 | +0.2 | +0.2 | −0.4 | −1.4 | −1.6 | 0.0 | +0.2 | 0.0 |
| q, 300–600 hPa | −0.81 | +0.3 | +0.7 | −0.2 | −1.9 | −1.9 | −1.1 | +0.2 | −1.1 |
| q, 600 hPa–surface | +0.19 | +0.3 | +0.5 | +0.8 | −0.5 | +0.6 | −0.3 | +0.2 | −0.2 |
| T +5 K above 50 hPa | +1.11 | (+1.1 in every band) | | | | | | | |
| T and q, the 50 hPa T bias carried to the top | +2.63 | +0.7 | +2.9 | +2.8 | +1.4 | +2.6 | +2.0 | +3.4 | +6.9 |

ERA5's humidity lowers the model's clear-sky OLR: the drier tropical
upper and middle troposphere outweighs the moister extratropics and
850 hPa. ERA5's temperature raises it, by 2.45 W/m² at fixed q or 1.27
at fixed relative humidity, the RH then adding +0.24. The temperature
biases behind it:

- **The tropical and subtropical lower troposphere** is 0.7–1.7 K cold
  at 600–925 hPa with the sea at WOA23's surface temperature, coldest and
  moistest at 850 hPa (−1.6 K, RH +0.10 to +0.12): the trade boundary
  layer reaching through 850 hPa, where ERA5 is above the inversion. It
  is the cold-tongue inversion "too high and smeared" of
  docs/equatorial-wind-bias.md and the low cloud at twice Earth's cover
  of docs/toa-radiation-bias.md (finding 3), and it is worth +1.2 W/m²
  (+2.1 and +2.3 in the two tropical bands).
- **The polar and midlatitude tropopause** is 4–14 K cold at
  100–200 hPa (60–90N −10.8, −13.8 and −13.3 K at 100, 150 and 200 hPa;
  40–60N −7.8 K at 150 hPa), under a lowermost stratosphere 1.1–2.2
  times as moist as ERA5's, which is itself moist there (Krüger et al.
  2022). A moist lowermost stratosphere cools itself in the infrared
  (Charlesworth et al. 2023). Worth +1.0 W/m² (50–300 hPa), most of it
  poleward of 40°.
- **The stratosphere** is cold: 3–7 K at 50 hPa against ERA5 at every
  latitude but the Antarctic's (−0.5 K), and against the AFGL profiles 5–11 K from 14 to 100 hPa
  outside the poles and 14–18 K at the tropical stratopause (finding 6). At
  0.22 W/m² per kelvin above 50 hPa, a 5 K bias is worth +1.1 W/m².

### 3. The stratosphere's missing shortwave heating

Between 5 and 150 hPa the radiation's shortwave heats AER's four
standard atmospheres 4–36 percent less than RRTMG for the same sun
(`scripts/radiationBenchmark.mjs` columns on bl42, surface albedo 0.2,
overhead sun, K/day, CLIRAD-SW's near infrared / the bands of finding 5 /
RRTMG):

| hPa | TROP | MLS |
|---|---|---|
| 1.6 | 29.03 / 30.01 / 28.88 | 26.66 / 27.65 / 26.15 |
| 7.5 | 8.86 / 9.43 / 9.32 | 9.19 / 9.76 / 9.57 |
| 14 | 4.60 / 5.06 / 5.19 | 4.79 / 5.24 / 5.34 |
| 24 | 2.83 / 3.21 / 3.29 | 2.93 / 3.30 / 3.42 |
| 38 | 1.72 / 2.02 / 2.10 | 1.80 / 2.08 / 2.17 |
| 54 | 1.13 / 1.39 / 1.44 | 1.17 / 1.41 / 1.48 |
| 71 | 0.70 / 0.94 / 1.00 | 0.78 / 0.99 / 1.04 |
| 102 | 0.48 / 0.67 / 0.70 | 0.44 / 0.62 / 0.67 |
| 141 | 0.41 / 0.56 / 0.57 | 0.38 / 0.54 / 0.57 |

RRTMG_SW (AER's column model, built here and run by
`scripts/nearInfraredFit.mjs`) puts the deficit in the near infrared. Run
with each gas alone, its direct beam heats TROP's stratosphere at 14, 54
and 102 hPa (overhead sun) by: CO₂ 0.44, 0.25 and 0.21 K/day (its 4.3,
2.7 and 2.0 µm bands); water vapor 0.19, 0.15 and 0.14; CH₄ 0.03, 0.02
and 0.03; O₂ 0.08, 0.08 and 0.07, about the model's. The four together
heat 0.67, 0.44 and 0.40 K/day; the model's near infrared heats 0.13,
0.16 and 0.19, its CO₂ 0.05 at each.

CLIRAD-SW's near infrared (`js/physics/shortwaveGases.module.js`)
starves the stratosphere three ways: the paths of vapor, O₂ and CO₂ scale
as (p/300 hPa)^0.8, a strong-line law with no Doppler floor, which divides
a layer's path by 5 at 40 hPa and by 12 at 14 hPa; CO₂ is a√u, calibrated
to what it adds to the vapor's absorption over a whole column; CH₄ is
absent. `VAPOR_STRENGTH` (1.48), fitted to RRTMG's column absorption,
gives the troposphere what the stratosphere misses. A Doppler floor on
the scaling alone recovers 10–30 percent of the deficit from 100 to
14 hPa (most of the vapor's, little of the CO₂'s) and has O₂ heat two to
three times RRTMG's at 7–14 hPa, so the remedy is a refit of the near
infrared (finding 5). With radiative damping times of weeks to months in
the lower stratosphere (Fels 1982), a deficit of 0.2–0.6 K/day under an
overhead sun is worth several kelvin; finding 5's arm measures it.

### 4. The humidity's spread

Each cell's q replaced by its own eight-state mean lowers the model's
clear-sky OLR by 1.96 W/m² (T as well: 2.03), evenly across the bands.
OLR goes roughly as −ln q, so at a fixed mean a broader humidity
distribution emits more (Pierrehumbert 1995). The swaps attach the
model's 14 days of variability to ERA5's 30-year mean, which also holds
interannual variability and ERA5's sub-cell variance. They therefore
read low against a real atmosphere of the same mean by an amount of this
order. This term is not measured here.

### 5. The near infrared in bands

The radiation option `nearInfrared: 'bands'` (`nearInfraredBands` in
`js/physics/shortwaveGases.module.js`, its table
`js/physics/nearInfraredTable.module.js` written by
`scripts/nearInfraredFit.mjs`) replaces CLIRAD-SW's near infrared with a
k-distribution in each of RRTMG_SW's bands 16–24 and 29
(820–16000 cm⁻¹). Each gas (water vapor, CO₂, CH₄, O₂) has four terms,
each coefficient scaled by the Voigt pressure √(p² + p_D²) to a fitted
power and by exp(t (T − 250 K)), with vapor's self continuum; a band's
gases overlap at random. The beam's magnification 35/√(1224μ² + 1) and
the reflected beam's 5/3 stay. Ozone and the visible band's vapor stay
CLIRAD-SW's, the visible vapor's strength refitted to RRTMG's bands 25–28
(0.68 against 1.48), since the bands take all of band 24 (625–778 nm).

The fit is to RRTMG_SW's direct beam over 25 columns: AER's four
standard atmospheres, the five AFGL atmospheres, and sixteen columns of
thirteen's day 1825, from the warm pool to the Antarctic plateau. Each
column is run at sun cosines 1, 0.5 and 0.25, with each gas alone and
with the four together. Each gas is fitted alone, then a band's gases
together against the four-gas run: multiplying the gases' own
transmissions from RRTMG takes 1.3–4.5 W/m² too much from a column
(2.5 on average). Against the four-gas run the heating misses by
0.027 K/day rms, both between 5 and 200 hPa and below, and the column's
absorption by +0.19 W/m² on average (rms 0.28, −0.36 to +0.68).
Overhead, the fit stays 2–13 percent short at 10–30 hPa (11 at the
median); weighting the layers above 100 hPa three times over moves that
by under a percent.

On the benchmark's four atmospheres under an overhead sun (table in
finding 3) the atmosphere absorbs 0.1–0.5 percent less than RRTMG
(CLIRAD-SW 0.1–1.3 percent more), and the layers between 5 and 150 hPa
heat within −9 and +3 percent of RRTMG (CLIRAD-SW −36 to −4). At
1.6 hPa the heating is 1.1–1.5 K/day above RRTMG's: that is the ozone's
excess there, which CLIRAD-SW's missing near infrared had hidden.

What CO₂ adds is still too much. On MLS, removing it takes 2.6 W/m² from
the atmosphere's absorption under a 60° sun (RRTMG 1.5, CLIRAD-SW 3.3,
calibrated to Chou and Suarez's line-by-line 3.3), and doubling it adds
0.8 (RRTMG 0.6, CLIRAD-SW 1.4); overhead, doubling adds 1.5 (RRTMG 1.1,
CLIRAD-SW 2.2). A joint fit under random overlap gets a band's sum
right, not each gas's margin. The budget's terms split a band's
absorption among its gases by what each absorbs alone (MLS overhead:
CO₂ 11.4 W/m², O₂ 7.0, CH₄ 1.8), so they are not what each gas adds.

Both engines carry it: `test/gpuModel.test.mjs` holds an N=8 parity test
(the layers above 100 hPa warm 0.357 K/day on each, every layer within
1.7e-4 K/day), and `test/gasRadiation.test.mjs` the benchmark's numbers.
`scripts/compareEngines.mjs` (env RADIATION) on eleven128_day1825 flags
the same fields with the option on as with it off. On the M5 Max the
passes that hold the radiation take 1.80 → 2.22 ms (0.8 a step) and
2.60 → 3.34 ms (0.3 a step) at N=128 (`scripts/profileGpu.mjs` on
eleven128_day1825, alternated twice): 0.56 ms a step, 1.3 percent of
42 ms. The N=64 fixed-SST arm ran 4.1 s/day against 3.7.

The arm: N=64 from thirteen64_day1825 over WOA23's sea surface for 30
days, against the same arm with CLIRAD-SW's near infrared
(`scripts/fixedSst.mjs`, states every 3 days). The global-mean warming,
with the approach dT = A (1 − e^(−t/τ)) fitted at each level:

| hPa | 1.6 | 3.5 | 7.3 | 24 | 37 | 52 | 69 | 99 | 161 |
|---|---|---|---|---|---|---|---|---|---|
| day 30 (K) | 2.7 | 2.5 | 2.6 | 2.4 | 2.3 | 2.1 | 1.9 | 1.1 | 0.3 |
| A (K) | 2.7 | 2.6 | 3.0 | 3.4 | 3.8 | 4.4 | 4.5 | 2.7 | |
| τ (days) | 8 | 9 | 14 | 24 | 34 | 46 | 56 | 48 | |

At equilibrium the warming at 50 hPa, about 4.4 K, is the size of the
stratosphere's cold bias against ERA5 there (3–7 K, finding 2). Added
above 100–150 hPa to the control's day-30 state, it raises the clear-sky
OLR by 1.0–1.1 W/m² (clear-sky longwave of every eighth column). At day
30 the stratosphere's part is +0.7 W/m² and the pair's net −0.03: its
troposphere is 0.06–0.25 K cooler, within the ±0.3 K the pair's
difference wanders over the 30 days.

### 6. The upper stratosphere, 0.1–3 hPa

The control arm's last state (day 1855, the sun at 11.6° N) is colder
than AFGL's tropical atmosphere in every layer above 10 hPa: 14 K in the
top layer, 17 K at 0.64 hPa and 7–8 K at 1.6–7 hPa. Its clear-sky
heating over 15° S–15° N, the mean of 24 hours on every sixteenth column
with the model's ozone, at its own temperatures and with AFGL's above
8 hPa, in K/day, with `nearInfrared: 'clirad'` (with 'bands' in
parentheses; the equilibrium is where the net heating vanishes,
interpolated between the two temperatures, extrapolated at 1.6 hPa):

| hPa | T, K | AFGL | shortwave | longwave | net | net at AFGL's T | equilibrium, K |
|---|---|---|---|---|---|---|---|
| 0.15 | 227.4 | 241.2 | 4.21 (5.26) | −4.42 | −0.22 (+0.83) | −2.25 (−1.22) | 226 (233) |
| 0.64 | 250.1 | 266.9 | 9.40 (9.95) | −8.03 | +1.37 (+1.92) | −3.12 (−2.55) | 255 (257) |
| 1.59 | 256.6 | 264.7 | 10.37 (10.80) | −8.50 | +1.88 (+2.30) | +0.92 (+1.35) | about 272 (276) |
| 3.54 | 244.5 | 251.5 | 6.46 (6.78) | −5.47 | +1.00 (+1.31) | +0.08 (+0.40) | 252 (255) |

In a steady state a layer's net radiative heating is what the
circulation removes. At 1.6–3.5 hPa the radiation alone would hold
AFGL's temperature or more, and the tropical ascent's 1.0–1.9 K/day of
adiabatic cooling holds the layers 7–8 K below it. At 0.64 hPa the two
share the bias: the radiation's own equilibrium is 12 K below AFGL's (10
with the bands), and the ascent holds the layer 5 K below that. In
30–60° S, in autumn, the circulation instead warms 0.15–3.6 hPa by
1.0–2.4 K/day, and in 30–60° N, in spring, it cools them by 0.4–0.5;
AFGL has no spring or autumn profile to set against them.

The radiation's errors, each against RRTMG on the same column, in K/day:

| | 0.15 | 0.65 | 1.6 | 3.6 hPa |
|---|---|---|---|---|
| CLIRAD-SW's near infrared, the tropical mean above (what the bands add back) | −1.05 | −0.55 | −0.43 | −0.32 |
| the longwave table against RRTMG_LW, AER's tropical atmosphere | +0.37 | −1.12 | +0.62 | +0.03 |
| the same, midlatitude summer, midlatitude winter, subarctic winter | +1.71, −0.86, −2.16 | +0.28, −1.33, −0.79 | +0.21, +0.51, +0.65 | −0.07, −0.03, +0.14 |
| ozone's ultraviolet and visible on the same ozone, overhead sun, tropical and midlatitude summer (percent) | −0.2, +13 | +0.9, −5.2 | +4.7, +6.6 | +3.0, +3.3 |
| O₂'s ultraviolet (RRTMG's bands 27–28), which the model lacks, under any sun | −0.23 | −0.23 | −0.23 | −0.22 |

The near infrared's part in the top layer includes O₂'s A and B bands
(RRTMG's band 24, 0.55 K/day under an overhead sun). The CO₂ amount is no error
but counts against AFGL's profile: 388.75 ppmv against the 354.8 of
AER's tropical column cools the four layers by 0.32, 0.52, 0.43 and
0.32 K/day (the table). At 0.64 hPa, at the 0.27 K/day per kelvin
between the two temperatures, the 16.8 K divides:

| part | K |
|---|---|
| the ascent | −5.1 |
| the longwave table | −4.2 |
| CLIRAD-SW's near infrared (the bands take it back) | −2.1 |
| CO₂ at 388.75 against 354.8 ppmv | −1.9 |
| the model's colder layers at 8–30 hPa | −0.7 |
| O₂'s ultraviolet (0.11 K/day in the tropical mean, the sun up half the day) | −0.4 |
| ozone's optics (about +1.5 percent between overhead and 60°) | about +0.5 |
| the rest: RRTMG's own equilibrium against AFGL's measured profile, and the layers below 30 hPa | −2.9 |

Two more 30-day arms from the same state (SST=woa) each take away one
driver at the top. Tropical means at day 30, temperature (K) and net
radiative heating (K/day):

| arm | 0.15 | 0.64 | 1.59 | 3.54 hPa |
|---|---|---|---|---|
| control | 227.4, −0.22 | 250.1, +1.37 | 256.6, +1.88 | 244.5, +1.00 |
| no lid friction (`SURFACE='{"lidFriction":null}'`) | 224.6, +0.51 | 249.5, +1.42 | 256.2, +1.98 | 244.9, +0.89 |
| no gravity-wave drag (`GRAVITY_WAVES=false`) | 224.9, +0.70 | 251.4, +1.19 | 260.0, +0.96 | 246.4, +0.58 |

The non-orographic gravity-wave drag (`js/physics/gravityWaves.module.js`,
which spreads what reaches the lid over the layers above 85 Pa) drives
half of the ascent's cooling at 1.6 hPa and 40 percent at 3.5 hPa:
without it 1.6 hPa warms 3.4 K in 30 days. The lid friction moves no
tropical layer below the top by more than 0.6 K. Each keeps the top
layer warmer, by 2.8 K (the friction) and 2.5 K (the waves).

No option comes of this finding. After the bands, the radiation's
largest error at 0.1–3 hPa is the longwave table's. It alternates
between layers: it over-cools 0.65 hPa and under-cools 1.6 hPa, and the
top layer goes either way by atmosphere. The table's CO₂ g-points reach unit
optical depth from the top at 0.011, 0.035 and 0.44 hPa, then not until
about 10 hPa; its opaque bins are 90-fold wide (`THICK_BIN` 4.5 in
`scripts/longwaveFit.mjs`). Narrower opaque bins with the same spectral
model (37, 41 and 46 g-points) put the tropical miss at 0.65 hPa at
−1.78, −0.91 and −1.57 K/day. Meanwhile the midlatitude atmospheres'
3–30 hPa rms grows from 0.06–0.07 to 0.09–0.20 K/day. A 300-step refit
with ten times the upper layers' weight moves the tropical miss to
−0.81, and the midlatitude summer's top layer from +1.71 to +1.94.

The likely cause is the spectral model itself, which the bins only
sample: above 7 hPa (p_D 727 Pa) its CO₂ lines keep a Lorentz profile at
the Doppler width, with far wings a Doppler line does not have. The
change to try is a line shape for the Doppler regime, refitted to
RRTMG_LW's top layers. Two items don't earn a pass of their own:

- ozone's optics, which stay within the spread between RRTMG_SW and
  CLIRAD-SW;
- O₂'s ultraviolet, worth 0.4 K at 0.64 hPa.

The ozone climatology is not a cause: it matches AER's tropical profile
within 3 percent at 0.3–5 hPa (12 percent less in the top layer). The
ascent is the circulation's question. Is 1.0–1.9 K/day at 1.6–3.5 hPa
too much against a reanalysis's residual circulation? And how much of
the gravity-wave drag's share is what it deposits at the lid?

## The ledger

| part | W/m² | |
|---|---|---|
| EBAF whole-region clear sky, annual | 265.9 | Loeb et al. 2018, 2020 |
| the arm's dates against the annual mean | −0.6 | ERA5's annual T, q and surface against its days 85–99 |
| the model's tropospheric temperature at fixed RH (50 hPa–surface) | −1.3 | lower troposphere −1.2, middle −0.2, upper and tropopause −1.0 (by layer at fixed q, which sums to 2.45) |
| the model's relative humidity | −0.2 | the drier tropical upper troposphere +, the moist 850 hPa and extratropics − |
| the stratosphere above 50 hPa | −1.1 | if its 50 hPa bias holds to the top (5 K); the bands' equilibrium warming gives +1.0 to +1.1 (finding 5) |
| gas optics | about −0.3 | bl42 against RRTMG by atmosphere, equinox weights |
| unexplained | −2.1 | the humidity's spread (finding 4), ERA5, EBAF |
| the model, eight states | 260.3 | |

## Ranked changes

| | change | references | clear-sky OLR | status |
|---|---|---|---|---|
| 1 | The stratosphere's shortwave heating: the near infrared's gases in RRTMG_SW's bands, fitted to its transmissions by level (finding 5), with nothing tuned to the model's climate. | Chou and Suarez 1999; Iacono et al. 2008; Mlawer et al. 1997 | +1.0 to +1.1 once the stratosphere adjusts (τ up to 56 days) | done: `nearInfrared: 'bands'`, the default |
| 2 | The extratropical tropopause: the lowermost stratosphere's vapor at 1.3–2.2 times ERA5's, and its temperature 5–14 K low. Trace the cross-tropopause vapor transport (vertical advection, diffusion, the ice cloud's sublimation) and the radiative cooling of the excess vapor. 1 adds part of the heating. | Krüger et al. 2022; Charlesworth et al. 2023 | up to +1.0 | to trace |
| 3 | The trade boundary layer through 850 hPa (−1.6 K, RH +0.1): the inversion's height and sharpness, in the cumulus and boundary-layer work already open (equatorial-wind-bias.md; toa-radiation-bias.md change 3) | Siebesma et al. 2003; Mellor and Yamada 1982 | up to +1.2, and less low cloud | with that work |
| 4 | Measure the unexplained 2.1 W/m²: EBAF Ed4.2's whole-region clear sky by latitude and month (toa_lw_clr_t; NASA Earthdata login), and ERA5's 6-hourly humidity for the spread | Loeb et al. 2020 | | data to fetch |
| — | Moisten or dry the upper troposphere: not a lever. Toward ERA5 it lowers the clear-sky OLR by 0.9 W/m². | | | ruled out |

The humidity's flat upper-tropospheric profile needs no trace: it is
ERA5's. The sources the question named (the deep plumes' detrainment, the
ice cloud's fall and sublimation, the resolved ascent) leave the tropical
upper troposphere 0.7–0.8 times as moist as ERA5 at 150–250 hPa, within
the spread of reanalyses and AIRS there (AIRS, which sees clear and
partly cloudy scenes, is drier than ERA5 in the convective tropics).

Change 1 is `nearInfrared: 'bands'` (finding 5), the default. Its effect
on the OLR takes the stratosphere's one to two months to show: a 30-day
arm reaches about half of it at 50–70 hPa.

## References

- Charlesworth, E., F. Plöger, T. Birner, et al., 2023: Stratospheric
  water vapor affecting atmospheric circulation. Nat. Commun., 14, 3925.
- Chou, M.-D., and M. J. Suarez, 1999: A solar radiation parameterization
  for atmospheric studies. NASA TM-1999-104606, vol. 15.
- Fels, S. B., 1982: A parameterization of scale-dependent radiative
  damping rates in the middle atmosphere. J. Atmos. Sci., 39, 1141–1152.
- Freidenreich, S. M., and V. Ramaswamy, 1999: Solar radiation absorption
  by carbon dioxide, overlap with water, and a parameterization for
  general circulation models. J. Geophys. Res., 104, 31389–31409.
- Hersbach, H., et al., 2020: The ERA5 global reanalysis. Q. J. R.
  Meteorol. Soc., 146, 1999–2049.
- Iacono, M. J., et al., 2008: Radiative forcing by long-lived greenhouse
  gases: calculations with the AER radiative transfer models. J. Geophys.
  Res., 113, D13103.
- Krüger, K., A. Schäfler, M. Wirth, M. Weissmann, and G. C. Craig, 2022:
  Vertical structure of the lower-stratospheric moist bias in the ERA5
  reanalysis and its connection to mixing processes. Atmos. Chem. Phys.,
  22, 15559–15577.
- Loeb, N. G., et al., 2018: Clouds and the Earth's Radiant Energy System
  (CERES) Energy Balanced and Filled (EBAF) top-of-atmosphere (TOA)
  Edition-4.0 data product. J. Climate, 31, 895–918.
- Loeb, N. G., et al., 2020: Toward a consistent definition between
  satellite and model clear-sky radiative fluxes. J. Climate, 33, 61–75.
- Mellor, G. L., and T. Yamada, 1982: Development of a turbulence closure
  model for geophysical fluid problems. Rev. Geophys., 20, 851–875.
- Mlawer, E. J., et al., 1997: Radiative transfer for inhomogeneous
  atmospheres: RRTM. J. Geophys. Res., 102, 16663–16682.
- Pierrehumbert, R. T., 1995: Thermostats, radiator fins, and the local
  runaway greenhouse. J. Atmos. Sci., 52, 1784–1806.
- Rasp, S., et al., 2024: WeatherBench 2: a benchmark for the next
  generation of data-driven global weather models. J. Adv. Model. Earth
  Syst., 16, e2023MS004019.
- Siebesma, A. P., et al., 2003: A large eddy simulation intercomparison
  study of shallow cumulus convection. J. Atmos. Sci., 60, 1201–1219.
