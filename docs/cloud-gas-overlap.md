# Gas absorption under clouds: vaporShielding and the layered shortwave

Why the model's clouds absorb 9.7 W/m² of sunlight above what a clear sky
does where observations and the reference models say 2 to 4, and what two
structural fixes do to that and to the planet's reflection, and what the
brighter planet asks of the clouds. Measured on 2026-10-25 on the CPU engine at
fourteen64's days 1460, 1551, 1643 and 1734 (March, June, September and
December of year 5, eight times of day each) and r18flux's days 91 and 183, on
the single columns of `test/vaporShielding.test.mjs` and
`test/layeredShortwave.test.mjs`, on CIRC Phase I's cases 6 and 7, and against
ESA Cloud_cci retrievals of MODIS-Aqua and AVHRR NOAA-18.

## The question

The `sequential` shortwave sends the beam through every layer's gases
first (CLIRAD-SW's ozone and water vapor, the fitted near-infrared
k-distributions, O2 and CO2) and only then lets the cloud act on what is
left: reflect `d / (d + 2 mu)` of it, with `d` the column's delta-scaled
cloud depth, and keep `exp(-cloudSolarAbsorption W)` of what it reflects or
transmits, with one coefficient (0.4 m²/kg) for all spectra and every
cloud. Three things follow. A cloud shields no gas below it, so a thick
high cloud leaves the troposphere's near-infrared absorption untouched
where the light that reaches the low levels is a fraction of the clear
sky's. The light a cloud reflects is cut in the visible as much as in the
near infrared, though only the near infrared is absorbed. And the cloud's own
absorption is one number for every phase and size.

## Three schemes

| `radiation.shortwave` / `vaporShielding` | what the cloud does | reference |
|---|---|---|
| `sequential` (default) | the gray stream above | |
| `vaporShielding: 'cloudTop'` | one effective cloud top per column; the near-infrared stream passes the gases around it as in a column with a reflecting layer: gases above the top see the beam and the light the cloud reflects, gases below see the light it transmits and the light the surface reflects | Lacis and Hansen 1974, with the overlap of Hogan and Illingworth 2000 |
| `shortwave: 'layered'` | every layer of the cloudy region holds the gases, Rayleigh scattering, aerosol and cloud of the slice (32 near-infrared g-points and 2 visible slices, the ultraviolet ozone bands left to the clear column); delta-Eddington two-stream, layers added over the surface; cloud depth at the model's four lognormal nodes | Joseph, Wiscombe and Weinman 1976; Meador and Weaver 1980; Slingo 1989 (liquid) and Ebert and Curry 1992 (ice), cloud optics by band |

Both options leave the clear-sky pass, the clear columns and the mixed-layer
deck's stream (2.6 % of the area, 4 % of the ocean's, in the days above) as
they were, and both keep absorbed plus reflected sunlight equal to the
insolation to rounding, held radiation included. `layered` adds to the clear
column of the sequential scheme the cover times what the cloudy region's
layered solution adds to the clear region's, so a column without cloud is
the sequential scheme's bit for bit; `cloudSolarAbsorption` does not act on
it.

## The folded k-distribution

The near-infrared table (`js/physics/nearInfraredTable.module.js`, fitted to
RRTMG_SW) has 230 combinations of one term per gas in its ten bands.
`scripts/layeredFit.mjs` merges them to 32 g-points by taking, pair by
pair, the two groups of a band whose merger changes least the direct and the
surface-reflected transmission to each interface on the four standard
atmospheres at humidities 0.25, 1 and 2 and magnifications 1, 2 and 4, with
a scale on each group's depth that takes back the absorption the mean of
unlike depths adds. Against the full table the clear sky's direct beam
absorbs 1.2 W/m² of 1360 differently on average (worst 4.1; the worst layer
3.2) over the benchmark columns, their humidities and three magnifications.
In cloudy columns (24 single-column cases) the atmosphere's
absorption from 32 g-points differs from that of the 210 unfolded
combinations by 0.07 % of the beam on average (worst 0.20 %), the
reflection by at most 0.13 %.

## The clear region and the cloud's difference

The layered column is the sequential scheme's clear column plus the cover times
the layered cloudy region's difference from the layered clear region. A column
without cloud is then the sequential scheme's bit for bit, and a thin cloud
changes it by its own effect; but the sum is the cloud's effect only where the
layered clear region is the sequential one, because under an opaque cloud it
becomes the sequential clear column less the layered clear region plus what the
cloud leaves. Two differences between the clear regions are set right in the
visible slices:

- CLIRAD-SW's ozone bands 1 to 6 (the ultraviolet, 2.1 % of the beam) are
  absorbed by the ozone before any cloud, 1.6 % of the beam at an overhead sun
  in the tropics and all of them toward the horizon; in the layered slices they
  were light, reflected by the cloud and unabsorbed. The slices now leave them
  to the clear column, in the column's own ozone as `ozoneAbsorptivity` finds
  it, and weigh what remains (0.3636 of the beam at an overhead sun), with the
  visible ozone and vapor coefficients raised to take what the sequential scheme
  takes per unit path in bands 7 and 8.
- The visible slices hold 0.279-0.607 um where the sequential scheme's Rayleigh
  bands hold 0.297-0.683 um, and the mean Rayleigh depth of the first is 1.329
  times that of the second; the slices' depths are those of the sequential scheme
  times 1.329 (`scripts/rayleighReference.mjs` with `OZONE=0.0208 VISIBLE=0.3844`).

The layered clear region less the sequential clear column, in % of the beam,
over the tropical atmosphere (the midlatitude one is within 0.1 of it), with the
ultraviolet bands in the slices and the sequential scheme's Rayleigh depths
(before) and as the slices are now:

| zenith | surface before | reflected before | atmosphere before | surface | reflected | atmosphere |
|---|---|---|---|---|---|---|
| 1 | +2.05 | -0.33 | -1.72 | -0.30 | +0.20 | +0.11 |
| 30 | +2.15 | -0.38 | -1.77 | -0.31 | +0.20 | +0.11 |
| 60 | +2.64 | -0.59 | -2.05 | -0.16 | +0.11 | +0.05 |
| 80 | +4.88 | -1.16 | -3.72 | +1.83 | -0.63 | -1.21 |

(the atmosphere is the layered atmosphere plus the ultraviolet share the slices
leave out). Above a zenith angle of 70 degrees the two clear regions differ
more, the sequential scheme scattering with the cosine and the layered one with
the cosine of the airmass that the gas paths take.

### The floor under a thick cloud near the horizon

Where the two clear regions differ, the cloud's difference can ask more of the
surface than the clear column gives it: under a thick cloud near the horizon the
sum of the sequential clear column's surface absorption and the layered
difference would be below nothing. The column's surface absorption,
downward flux and direct flux are each kept at no less than nothing
(`max(difference, -clear)`). What the floor does not take from the surface goes
nowhere but to the reflection: the heating profile is the layered solution's, and
the top-of-atmosphere reflection is the beam less what the atmosphere and the
surface absorb, so absorbed plus reflected equals the beam. Over the 163,847
cell-times of the four seasonal days that closure is at most 4.9e-16 of the
beam in every record, with the floor holding the surface's absorption in 1.24 %
of the sunlit samples (any of the three fluxes in 7.61 %). The sunlight the floor
sends from the surface to the reflection is 0.009 W/m² in the global mean:

| zenith angle | below 70 | 70-80 | 80-85 | 85-90 |
|---|---|---|---|---|
| global-mean sunlight moved from the surface to the reflection, W/m² | 0.000 | 0.000 | 0.002 | 0.007 |
| sunlit samples where the surface's absorption is floored | 0.0 % | 0.2 % | 2.0 % | 11.8 % |

(the downward and direct fluxes, floored the same way, gain 0.013 and 0.074
W/m²). The floor is small against the two-stream method's own error at low sun:
against Monte Carlo, in domain means of A-Train clouds, Barker et al. (2015) find
the delta-Eddington top-of-atmosphere reflected flux off by about 10 W/m² at a
cosine of the zenith angle above 0.9 and by about 3 W/m² the other way below 0.2
(for the approximations' accuracy against exact solutions see also King and
Harshvardhan 1986).

## Results, W/m²

Global means of the four seasonal days, offline on the saved states with
their own clouds and decks (`scripts/shortwaveBudget.mjs`); ERA5 2020's 20th
of the same months; CERES EBAF Ed4.0 as Loeb et al. (2018, Table 5, July
2005 to June 2015) give it: reflected 99.1 of 340.0, so ASR 240.9, clear-sky
286.7 and SWCRE -45.8 (the 2001-2024 monthly means in the EBAF file give
-45.1).

| | sequential | cloudTop | cloudTop, kappa 0.55 | layered | ERA5 | CERES |
|---|---|---|---|---|---|---|
| ASR | 242.8 | 237.8 | 240.3 | 235.0 | 243.7 | 240.9 |
| atmosphere | 82.8 | 78.1 | 81.2 | 76.9 | 79.6 | |
| cloud-induced atmosphere | 9.7 | 5.0 | 8.1 | 3.8 | 2.4 | |
| surface | 160.0 | 159.8 | 159.1 | 158.1 | 164.1 | |
| SWCRE | -50.0 | -55.1 | -52.5 | -57.9 | -46.7 | -45.8 |
| June 0-20N ASR | 312.6 | 302.2 | 307.6 | 292.6 | 326.5 | |
| June 0-20N SWCRE | -77.7 | -88.1 | -82.7 | -97.7 | -44.4 | -52.1 |
| June 40-60N ASR | 334.8 | 327.5 | 332.1 | 320.3 | 316.7 | |
| June 40-60N SWCRE | -94.1 | -101.5 | -96.9 | -108.7 | -99.3 | -90.4 |

The clear-sky ASR is 292.9 in every column. CERES's 286.7 is that of cloud-free
scenes, which differ from the clear sky of cloudy columns by 0.24 W/m² in the
global shortwave (Kato et al. 2013, as quoted by Loeb et al.), so the model's
clear sky absorbs 6.2 W/m² more (ERA5: 290.4, 2.5 less): the surface 219.7
against ERA5's 213.2, the atmosphere 73.1 against 77.2. The Kepler orbit
changes none of the differences between the columns by more than 0.05 W/m²
(annual means; June, near aphelion, is 3 % dimmer).

Shielding takes 5.9 W/m² of gas absorption, 0.8 of it returns in the cloud,
and the planet brightens by 5.0 W/m² (4.7 to 5.3 by season). The layered
scheme brightens it by 7.9 (7.4 to 9.0 by season): its cloud-induced
absorption, 3.8, is the closest of the three to ERA5's 2.4, and its reflection
is what a delta-Eddington cloud of the model's optical depths gives (the
single columns below). The tropics carry it: June 0-20N loses 20 W/m² of ASR
(its SWCRE is -98 against CERES's -52), the clouds of that state hold more
water than the sky does.

## Where the layered scheme's brightening sits

Layered less sequential, the same clouds and decks, over the 163,847 cell-times
of the four seasonal days, global means in W/m²: reflected +7.85 (97.68 to
105.53), absorbed -7.85, of which the atmosphere -5.93 and the surface -1.92.
By latitude band:

| band | area | reflected, band mean | atmosphere | surface | share of the global change |
|---|---|---|---|---|---|
| 90S-60S | 6.7 % | +2.74 | -2.39 | -0.35 | 2 % |
| 60S-30S | 18.3 % | +7.53 | -5.64 | -1.89 | 18 % |
| 30S-0 | 25.0 % | +7.67 | -5.31 | -2.36 | 24 % |
| 0-30N | 25.0 % | +10.57 | -8.00 | -2.57 | 34 % |
| 30N-60N | 18.3 % | +7.76 | -6.07 | -1.69 | 18 % |
| 60N-90N | 6.7 % | +4.62 | -4.40 | -0.22 | 4 % |

By cloud class, ISCCP-style: the pressure of the first cloud the light meets
(the level by optical depth) below 440 hPa, from 440 to 680 and above 680, and the
column's in-cloud visible optical depth (its grid-mean depth over its cover) from
0.3 to 3.6, from 3.6 to 23 and above 23. Shares of the sunlit samples, the
class's reflectance (reflected over incident sunlight) in the two schemes, and
the class's part of the global change:

| class | sunlit samples | reflectance, sequential / layered | change, W/m² | share |
|---|---|---|---|---|
| high, depth > 23 | 8.4 % | 0.481 / 0.649 | +4.78 | 61 % |
| high, 3.6-23 | 10.3 % | 0.397 / 0.431 | +1.18 | 15 % |
| high, 0.3-3.6 | 7.4 % | 0.191 / 0.196 | +0.13 | 2 % |
| mid, depth > 23 | 2.8 % | 0.466 / 0.533 | +0.66 | 8 % |
| mid, 3.6-23 | 9.1 % | 0.333 / 0.347 | +0.42 | 5 % |
| mid, 0.3-3.6 | 2.9 % | 0.331 / 0.333 | +0.01 | 0 % |
| low, depth > 23 | 1.8 % | 0.371 / 0.393 | +0.14 | 2 % |
| low, 3.6-23 | 35.4 % | 0.265 / 0.269 | +0.49 | 6 % |
| low, 0.3-3.6 | 6.4 % | 0.211 / 0.213 | +0.03 | 0 % |
| clear or depth < 0.3 | 15.3 % | 0.165 / 0.165 | 0.00 | 0 % |

The tropics carry 58 % of it and the towers 61 %: columns topped above 440 hPa
with a depth above 23, whose mean depth is 80 with 466 g/m² of liquid and 217 g/m²
of ice, reflect 0.48 of the sunlight in the sequential scheme and 0.65 in the
layered one, and 4.61 of their 4.78 W/m² comes out of the atmosphere's
absorption: the sequential scheme lets the gases below such a cloud take
sunlight that the cloud then does not reflect. The low clouds, 44 % of the sunlit
samples, change it by 0.7 W/m².

## The model's clouds against satellite retrievals

The sequential scheme's absorbed sunlight is near CERES's (242.8 against 240.9)
with clouds that are optically thick. The tropical means (30S-30N) of the model's
four seasonal days against the ESA Cloud_cci L3C v2.0 monthly products of
MODIS-Aqua and AVHRR NOAA-18 (Stengel et al. 2017):

| 30S-30N | model | MODIS-Aqua | AVHRR NOAA-18 |
|---|---|---|---|
| cloud cover | 57.8 % | 58.0 % | 54.9 % |
| cover by top: high, mid, low | 17.9, 7.8, 22.2 % (9.9 % below depth 0.3) | 22.6, 10.0, 24.5 % | 24.4, 8.2, 21.8 % |
| area topped above 440 hPa with depth > 23 | 7.3 % | 1.5 % | |
| in-cloud optical depth: mean (geometric mean) | 24.1 (11.5) | 6.9 (3.6) | 8.8 (4.2) |
| over ocean | 19.5 (10.2) | 6.4 (3.2) | 8.1 (3.8) |
| over land | 38.1 (16.6) | 8.8 (5.1) | 11.2 (5.5) |
| all-sky liquid water path, g/m² | 69.9 | 25.1 | 21.8 |
| all-sky ice water path, g/m² | 22.3 | 45.0 | 87.1 |
| liquid droplet radius over ocean, land, um | 11.8, 8.5 | 16.3, 11.7 | 13.9, 10.5 |
| global: in-cloud optical depth, mean (geometric mean) | 21.2 (10.4) | 10.6 (5.4) | 11.7 (5.8) |
| global: area topped above 440 hPa with depth > 23 | 8.3 % | 2.6 % | |

The model: cover, top pressure and in-cloud depth as above, weighted by area and
cover over the columns with cover above 0.05; water paths are the column's
(liquid and ice by the radiation's temperature rule, resolved plus cumulus
water, below cloud tops too); the droplet radii are the options' (11.8 over sea,
8.5 over land). The retrievals see the topmost layer in daytime at the
overpass (13:30 for Aqua, about 14:00 for NOAA-18): cloud fractions and
all-sky paths are means of March, June, September and December 2007 to 2009, the
optical depths and radii (weighted by area and cloud fraction, by liquid or ice
share for the phase variables) of the same months of 2008, and the depth-pressure
frequencies the daytime joint histogram of 2008 on every sixth 0.5-degree cell,
area-weighted. What the retrievals do not see makes the model's liquid path an
upper bound of the comparison, and the overpass time and the
sensors' sampling matter at the tens of percent, not at factors of three.

Sources (ESA Cloud_cci, anonymous OPeNDAP at CEDA, `dap.ceda.ac.uk/thredds/dodsC/
neodc/esacci/cloud/data/phase-2/L3C/`; no login), each set in a directory of its own:

| set | source path | files | size |
|---|---|---|---|
| cloud fraction (total, high, mid, low), all-sky liquid and ice water path, 2007-2009, four months | `MODIS-AQUA/v2.0/<year>/<yyyymm>-ESACCI-L3C_CLOUD-CLD_PRODUCTS-MODIS_AQUA-fv2.0.nc` | 73 | 74.9 MB |
| the same, nine months (2007, 2008 and March 2009: the product has no June 2009) | `AVHRR-PM/v2.0/AVHRR_NOAA-18/<year>/<yyyymm>-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-18-fv2.0.nc` | 55 | 56.2 MB |
| cot, cot_log, cer_liq, cer_ice, cot_liq, cot_ice, cph, 2008, four months | MODIS-AQUA | 29 | 29.1 MB |
| the same | AVHRR NOAA-18 | 29 | 29.1 MB |
| hist2d_cot_ctp and nobs_day, every sixth 0.5-degree cell, 2008, four months | MODIS-AQUA | 9 | 45.1 MB |

The model's clouds hold the observed cover and 2.7 to 3.5 times the observed
optical depth in the tropics (2.0 and 1.8 times in the global mean). The tropical
in-cloud depth is 24 where the satellites give 7 to 9, the liquid water path about
three times theirs (the model's counts the water below the tops, which the
satellites do not see), the droplets smaller over ocean by 28 % than MODIS's and
15 % than AVHRR's (the 16.3 um of MODIS over the tropical ocean is itself about
4 um above the MISR-referenced values of Fu et al. 2019, which the model's
11.8 um matches), and the ice path half of MODIS's and a quarter of AVHRR's. The
excess is in the towers: columns topped above 440 hPa with a depth above 23 cover
7.3 % of the tropics where MODIS gives 1.5 %. The sequential scheme meets CERES
with them in part because it reflects too little of the thick ones: in CIRC case
6 (depth about 57) its top-of-atmosphere reflection is 25.8 W/m² (4.1 %) below
CHARTS's 628.8, the layered scheme's 7.5 W/m² (1.2 %) above.

## Single columns, tropical atmosphere, % of the beam

Atmosphere's absorption, reflection and surface's absorption, 30 degrees (60
degrees in the next lines of the test's output); stratocumulus 100 g/m² (10
um), mixed-phase 300 g/m² at 450 hPa, ice 300 and 100 g/m² at 250 hPa, two
layers (100 g/m² at 250 and 850 hPa) and cirrus of 5 g/m².

| | clear | sequential | cloudTop, kappa 0.55 | layered |
|---|---|---|---|---|
| stratocumulus | 21.1 / 8.8 / 70.1 | 24.6 / 43.6 / 31.8 | 27.3 / 41.7 / 31.0 | 27.1 / 44.0 / 28.9 |
| mixed 300 | | 30.6 / 51.0 / 18.4 | 26.4 / 56.6 / 17.1 | 22.2 / 61.6 / 16.3 |
| ice 300 | | 30.6 / 51.6 / 17.8 | 24.1 / 59.4 / 16.4 | 22.7 / 60.9 / 16.3 |
| ice 100 | | 24.5 / 38.6 / 36.9 | 19.1 / 45.3 / 35.5 | 22.1 / 42.5 / 35.3 |
| two layers | | 27.7 / 50.4 / 21.9 | 23.6 / 55.7 / 20.7 | 24.8 / 56.1 / 19.1 |
| cirrus 5 | | 21.3 / 12.4 / 66.3 | 20.5 / 13.4 / 66.1 | 21.2 / 12.1 / 66.6 |

## CIRC Phase I, the two cloudy cases (line by line: CHARTS)

Overcast liquid cloud, the case's profile and cloud water laid on the model's
48 layers, surface albedo the CHARTS surface reflection over its downward
flux (0.136 and 0.164), the case's incoming irradiance (967.8 and 994.5
W/m²). Errors against CHARTS in W/m².

| | TOA up | surface down | surface net | atmosphere absorbs |
|---|---|---|---|---|
| case 6, tau about 57: CHARTS | 628.8 | 92.1 | 79.6 | 259.4 |
| sequential | -25.8 | +14.9 | +12.9 | +13.0 |
| cloudTop, kappa 0.55 | -49.0 | +9.0 | +7.8 | +41.2 |
| layered | +7.5 | -3.8 | -3.3 | -4.2 |
| case 7, tau about 7: CHARTS | 356.4 | 473.7 | 396.0 | 242.1 |
| sequential | +11.9 | +17.3 | +14.4 | -26.4 |
| cloudTop, kappa 0.55 | -3.5 | +13.1 | +10.9 | -7.4 |
| layered | +7.6 | -7.5 | -6.3 | -1.3 |

Their clear-sky subcases (6c, 7b): TOA up +2.7 and +6.7, surface down +3.8 and
-7.2, atmosphere -6.0 and -0.8 (3 % and 0.4 %), the same in all schemes. The
layered scheme has the atmosphere's absorption right to 2 % and the other
fluxes to under 1 % of the incoming sunlight. Under the cloud of depth 7 the
sequential scheme's surface receives 181 W/m² as a direct beam and the layered
scheme's 47 (the delta-scaling counts the forward peak as direct; CHARTS's
unscattered beam is 0.06). Slingo's liquid asymmetry parameter (0.846 at 6.9 um
where Mie gives about 0.86) is the largest sensitivity left: raised by 0.02 it
moves the layered TOA by -10.9 and -18.0 W/m² in the two cases.

## The two engines

The CPU solver, `js/physics/layeredShortwave.module.js`, is the reference.
On the GPU the layered scheme is one pass, `layered`, between the radiation
and longwave kernels, present only with the option (the default passes and
their shader text are those of the sequential scheme). The radiation kernel
solves the clear sky as for the sequential scheme and, for each cell with
cloud outside the deck, leaves nine values in `PH_LYR` (whether the pass
applies, the sun's cosine, the beam, the two surface albedos, the cover, the
relative variance of the cloud's depth, the deck's fraction, the aerosol's
depth). The pass solves the clear region and each node of the cloud in each
of the 34 slices, adds the cloud's difference from the clear region,
weighted by the cover, to the stashed sunlight (the layers' heating, the
surface's absorption and fluxes, the top's absorption and reflection) and
returns at its first line for every other cell. Its 17 thread-private arrays
of K floats are local to the entry point, every loop that reads them back
counts up, and the previous layer's albedo and source are carried in scalars,
as the rules for the M5 and M1 Metal compilers in CLAUDE.md require. Held
radiation keeps the pass's results in the sequential scheme's buffers.

In single precision the layer's coefficients need care. The exponential form
of the delta-Eddington solution (`E = exp(-k tau')`, `1 - E^2`, `1 - k^2 mu^2`
in a denominator) cancels terms of order 1 in a thin layer that absorbs
little, and is wrong by 1e-4 of the beam in a layer's reflection and
transmission. Rounding every operation of the CPU solver to single precision
(`Math.fround`) reproduced the GPU's difference from the double-precision
solver, absorbed 2.1e-4, atmosphere 3.5e-4 to 4.4e-4 and reflected 3.1e-4
(99th percentile of the cells over the field's rms, N = 16, 4 steps); the
same rounding of the adding method alone left under 3e-6, so the differencing
of the fluxes into a layer's heating was not the cause. The layer is therefore
written with `Q = tau' (1 - exp(-2 k tau')) / (2 k tau')` and the factor
`1 - k mu0` divided out of the beam's terms: regular at `k = 0` and at
`k mu0 = 1`, with no guard at the singular cosine and no floor on the
co-albedo, equal to the exponential form to 2e-13 in double precision, and
good to 2e-7 in single precision over 4,000 layers from `tau'` of 1e-6 to
1000, co-albedos of 1e-9 to 1 and every cosine the model takes. The emulated
single-precision difference falls to 5e-7 to 2.4e-6, and the GPU matches
the CPU model as the sequential scheme does:

| N = 16, 4 steps, 99th percentile of the cells over the field's rms | absorbed | atmosphere | reflected | clear-sky absorbed | OLR | theta |
|---|---|---|---|---|---|---|
| sequential | 1.7e-5 | 7.1e-6 | 3.0e-5 | 2.4e-7 | 1.8e-4 | 1.4e-6 |
| layered | 2.6e-5 | 1.2e-5 | 3.9e-5 | 2.6e-7 | 2.0e-4 | 1.4e-6 |
| sequential, radiation every 4th step | 2.4e-5 | 1.0e-5 | 3.9e-5 | 2.8e-7 | 3.0e-4 | 2.4e-6 |
| layered, radiation every 4th step | 3.6e-5 | 1.6e-5 | 5.0e-5 | 3.0e-7 | 2.6e-4 | 2.5e-6 |

The area means move alike: from the sequential to the layered scheme the
absorbed sunlight (W/m² per step sum) goes from 866.74 to 819.47 on the CPU and
from 866.86 to 819.47 on the GPU, the atmosphere's from 309.78 to 268.44 and
from 309.75 to 268.44. Under cloud ten times as thick the layered cells differ
from the CPU's by 2.9e-5 (absorbed) and 1.5e-5 (atmosphere), and the least
downward shortwave at the surface is 0.000 W/m².

GPU memory the options add: `cloudTop` enlarges the longwave stash buffer K1
(3 K C floats) by 78.0 MB at N = 128 with 48 layers (93.7 MB with 60) and by
19.5 MB at N = 64 (23.4 MB); `layered` adds the nine `PH_LYR` slots, 5.9 MB at
N = 128 and 1.5 MB at N = 64, and its private arrays are thread memory, not
buffers.

## The held radiation between calls

With `radiationEvery` above 1 (the page's cadence: 11.25 minutes of model time,
every 4th step at N = 128, every 2nd at N = 64) both engines run the full
radiation only at a step whose number is a multiple of `radiationEvery`, and at
a column's first step, under the mean cosine of the zenith angle over the steps
the call covers (`callCosZenith`: the positive cosines weighted by themselves,
so that a cell whose sun rises between calls still has sunlight to scale). The
call's layer heating, surface fluxes and top-of-atmosphere fluxes are kept, and
every step between calls the shortwave is those values times the ratio of the
step's cosine, times the sun's distance factor, to the call's: the heating
profile keeps the call's shape and the cloud and gas amounts of the call, while
the surface's emission change since the call is absorbed by the layers in the
call's shares (Hogan and Bozzo 2015) and the rest of the longwave stays the
call's. The layered scheme changes nothing in this: its heating, surface and top
terms enter the call's values, so the held scaling and the closure of absorbed
plus reflected sunlight to the insolation apply to them as to the sequential
scheme's (the CPU and GPU tests hold both). What the hold cannot follow is the
change of the heating's shape with the zenith angle and the change of the
clouds and the humidity between calls.

## What would close the gap to CERES (none applied)

The reflection is the cloud's: its water (above), its optical depths (droplet
radius 11.8 and 8.5 um, Han et al. 1994; MODIS 2.1 um's ocean mean 16.3 to 16.6 um
falls to 11.2 to 12.1 once its bias against MISR is removed, Fu et al. 2019), its
inhomogeneity (cloudSpread 0.75, 0.75 +- 0.18, Shonk et al. 2010) and its cover.
Absorbed sunlight and the shortwave cloud effect (SWCRE) of the four seasonal
days with the clouds or their optics changed offline, in W/m² (CERES: 240.9 and
-45.8; the clear-sky absorbed sunlight is 292.9 against CERES's 286.7 in every
row):

| change | sequential | cloudTop | layered |
|---|---|---|---|
| none | 242.8 / -50.0 | 237.8 / -55.1 | 235.0 / -57.9 |
| cloud water above 440 hPa x 0.7 | 243.9 / -48.9 | | 236.8 / -56.0 |
| cloud water above 440 hPa x 0.5 | 244.9 / -48.0 | | 238.3 / -54.5 |
| cloud water above 440 hPa x 0.3 | 246.1 / -46.8 | | 240.2 / -52.7 |
| cloud water above 440 hPa x 0.15 | 247.3 / -45.6 | | 241.9 / -51.0 |
| cloud water at every level x 0.7 | 250.9 / -42.0 | | 245.0 / -47.9 |
| cloud water at every level x 0.5 | 257.8 / -35.1 | | 253.3 / -39.6 |
| droplet radii +20 % (sea 14.2, land 10.2 um) | | 241.8 / -51.0 | 239.3 / -53.5 |
| cloudSpread 1.0 | | 239.3 / -53.6 | 236.4 / -56.5 |
| both | | 243.2 / -49.7 | 240.6 / -52.3 |
| liquid asymmetry +0.02 (Slingo's 0.846 at 6.9 um; Mie about 0.86) | | | 237.5 / -55.4 |
| cloudSolarAbsorption 0.55 | 245.1 / -47.8 | 240.3 / -52.5 | |

The layered scheme reaches CERES's absorbed sunlight with the cloud water above
440 hPa at about 0.24 of its present amount, with the water at every level at
about 0.8 of it, or with the droplet radii and the inhomogeneity both at the
edge of their observed ranges (240.6); its SWCRE of -45.8 needs the water at
every level at about 0.65, with the clear sky's 6.2 W/m² excess as it is. The
same changes of the water put the sequential scheme 3 to 17 W/m² above CERES.
Cutting the water by a fifth would leave the tropical in-cloud depth near 19
(the depth scales with the water), two to three times the satellites', so
closing the balance does not by itself bring the clouds to the observed
thickness; that would also take the clear sky's excess absorption and the cover
into account.

## Cheaper layered variants (offline, not options)

The layered work is the number of slices (the g-points and the two visible
slices) times the regions each solves (the clear one and one per lognormal node)
times the layers. The same four days with the table folded to fewer g-points
(`POINTS=24 WRITE=1 node scripts/layeredFit.mjs` in a scratch tree) and with
two nodes in place of four:

| g-points | nodes | work (34 x 5 = 1) | direct beam absorption miss, W/m² of 1360 (mean, worst) | global ASR | SWCRE |
|---|---|---|---|---|---|
| 32 | 4 | 1.00 | 1.2, 4.1 | 235.0 | -57.9 |
| 24 | 4 | 0.76 | 2.0, 7.3 | 234.9 | -58.0 |
| 16 | 4 | 0.53 | 5.5, 29.9 | 234.8 | -58.1 |
| 32 | 2 | 0.60 | 1.2, 4.1 | 234.5 | -58.3 |
| 24 | 2 | 0.46 | 2.0, 7.3 | 234.5 | -58.4 |
| 16 | 2 | 0.32 | 5.5, 29.9 | 234.3 | -58.5 |

The global balance hardly notices either (0.1 to 0.6 W/m²), the clear-sky
absorption of the folded table does at 16 g-points.

## Limitations

The deck keeps the gray stream. The cloudy region is one region of cover c
with the cloud water of each layer over c, as the column's total depth has
always been taken; Tripleclouds (Shonk and Hogan 2008) is the next step if
the overlap needs more than its four nodes. The two-stream is
delta-Eddington (with the Eddington coefficient gamma2 held at or above 0,
which keeps diffuse reflection positive and makes the layers add exactly):
against Monte Carlo it misses a thick layer's reflection by 0.02 rms.
Diffuse light goes through an absorber at 1.75 rather than 5/3. The
visible slices carry the ozone of CLIRAD-SW's bands 7 and 8 and the vapor of
0.4-0.7 um with one coefficient each, and the ultraviolet bands 1 to 6 stay in
the clear column, so the ozone does not take cloud-reflected ultraviolet a
second time. Above a zenith angle of 70 degrees the layered and the sequential
clear regions differ by more than 0.3 % of the beam (1.8 % at 80 degrees), and
the surface's sunlight under a thick cloud is held at no less than nothing.
