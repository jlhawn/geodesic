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
  `profiles.mjs`, `stratProfile.mjs`, `zarr.mjs`, `reduceEra5.mjs`).

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
  outside the poles and 14–18 K at the tropical stratopause. At
  0.22 W/m² per kelvin above 50 hPa, a 5 K bias is worth +1.1 W/m².

### 3. The stratosphere's missing shortwave heating

The radiation's shortwave heats the lower stratosphere 10–35 percent less
than RRTMG for the same atmosphere and sun (`scripts/radiationBenchmark.mjs`
columns, bl42, overhead sun, K/day, model / RRTMG / model without ozone):

| hPa | TROP | MLS | SAW |
|---|---|---|---|
| 7.5 | 8.86 / 9.32 / 0.15 | 9.19 / 9.57 | 7.18 / 7.52 |
| 14 | 4.60 / 5.19 / 0.15 | 4.79 / 5.34 | 4.13 / 4.58 |
| 24 | 2.83 / 3.29 / 0.15 | 2.93 / 3.42 | 2.94 / 3.36 |
| 38 | 1.72 / 2.10 / 0.15 | 1.80 / 2.17 | 2.27 / 2.63 |
| 54 | 1.13 / 1.44 / 0.16 | 1.17 / 1.48 | 1.74 / 2.03 |
| 71 | 0.70 / 1.00 / 0.17 | 0.78 / 1.04 | 1.13 / 1.45 |
| 102 | 0.48 / 0.70 / 0.19 | 0.44 / 0.67 | 0.76 / 0.98 |
| 141 | 0.41 / 0.57 / 0.21 | 0.38 / 0.57 | 0.49 / 0.68 |

The upper stratosphere (1.6–3.6 hPa) matches. The deficit, 0.15–0.6
K/day, peaks with the ozone at 14–24 hPa and grows in relative terms
downward, where the Chappuis band carries the heating; at 100–170 hPa,
where the model's ozone heating is 0.2–0.3 K/day, part of it must be
the near infrared. Two pieces of `js/physics/shortwaveGases.module.js`
fit that pattern:

- ozone in CLIRAD-SW's eight bands of one absorption coefficient each
  (Chou and Suarez 1999, Table 3). Inside a band the cross-section spans
  one to two orders of magnitude (the Huggins bands at 310–340 nm), and a
  single mean coefficient spends the band's sunlight too high: the weakly
  absorbed part that reaches the lower stratosphere decays as exp(−k̄O)
  instead of as the slower mixture of its own coefficients. The column's
  ozone absorption is RRTMG's within a few W/m² (38.0 W/m², MLS
  overhead, about 23 of it below 0.4 µm, against RRTMG's 24.7 in the
  ultraviolet), so the error is where, not how much;
- vapor, O₂ and CO₂ paths scaled by (p/300 hPa)^0.8 without a Doppler
  limit. Near-infrared lines turn Doppler-broadened at about 50 hPa
  (CO₂ near 4000 cm⁻¹), 80 hPa (vapor near 6000 cm⁻¹) and 260 hPa (the
  O₂ A band), so the scaling starves the stratosphere of their
  absorption. The longwave already carries this limit (`voigtScale`).

The column totals agree with RRTMG because `VAPOR_STRENGTH` (1.48) is
fitted to the column's absorption, so what the stratosphere misses the
troposphere's vapor takes. With radiative damping times of weeks to
months in the lower stratosphere (Fels 1982), a 0.2–0.3 K/day (overhead)
deficit is worth several kelvin, the size of the bias.

### 4. The humidity's spread

Each cell's q replaced by its own eight-state mean lowers the model's
clear-sky OLR by 1.96 W/m² (T as well: 2.03), evenly across the bands.
OLR goes roughly as −ln q, so at a fixed mean a broader humidity
distribution emits more (Pierrehumbert 1995). The swaps attach the
model's 14 days of variability to ERA5's 30-year mean, which also holds
interannual variability and ERA5's sub-cell variance. They therefore
read low against a real atmosphere of the same mean by an amount of this
order. This term is not measured here.

## The ledger

| part | W/m² | |
|---|---|---|
| EBAF whole-region clear sky, annual | 265.9 | Loeb et al. 2018, 2020 |
| the arm's dates against the annual mean | −0.6 | ERA5's annual T, q and surface against its days 85–99 |
| the model's tropospheric temperature at fixed RH (50 hPa–surface) | −1.3 | lower troposphere −1.2, middle −0.2, upper and tropopause −1.0 (by layer at fixed q, which sums to 2.45) |
| the model's relative humidity | −0.2 | the drier tropical upper troposphere +, the moist 850 hPa and extratropics − |
| the stratosphere above 50 hPa | −1.1 | if its 50 hPa bias holds to the top (5 K) |
| gas optics | about −0.3 | bl42 against RRTMG by atmosphere, equinox weights |
| unexplained | −2.1 | the humidity's spread (finding 4), ERA5, EBAF |
| the model, eight states | 260.3 | |

## Ranked changes

| | change | references | clear-sky OLR | status |
|---|---|---|---|---|
| 1 | The stratosphere's shortwave heating: split the ozone bands whose coefficient varies most (6 and 7, the Huggins bands) into sub-bands of their own coefficients, give the near-infrared path scaling the Doppler limit the longwave has, then refit VAPOR_STRENGTH, all against RRTMG's heating profiles by level and band. No free parameters: the targets are the reference's heating. | Chou and Suarez 1999; Iacono et al. 2008; Mlawer et al. 1997; Freidenreich and Ramaswamy 1999 | about +1 to +1.5 once the stratosphere adjusts (months) | needs RRTMG's band-by-level heating (the repo holds totals by level and by band, not both) |
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

No model change was made: the lead, 1, is a refit of the shortwave gas
scheme against a reference the repository does not yet hold, and its
effect on the OLR takes the stratosphere's months to show, beyond a 20-day
fixed-SST arm.

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
