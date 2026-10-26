# Gas absorption under clouds: vaporShielding and the layered shortwave

Why the model's clouds absorb 10.3 W/m² of sunlight above what a clear sky
does where observations and the reference models say 2 to 4, and what two
structural fixes do to that and to the planet's reflection, and what the
brighter planet asks of the clouds. Measured on 2026-10-25 on the CPU engine at
fourteen64's days 1460, 1551, 1643 and 1734 (March, June, September and
December of year 5, eight times of day each) and r18flux's days 91 and 183, on
the single columns of `test/vaporShielding.test.mjs` and
`test/layeredShortwave.test.mjs`, on CIRC Phase I's cases 6 and 7, and against
ESA Cloud_cci retrievals of MODIS-Aqua and AVHRR NOAA-18 and ERA5's cloud water.

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
cell-times of the four seasonal days that closure is at most 5.4e-16 of the
beam in every record, with the floor holding the surface's absorption in 1.26 %
of the sunlit samples (any of the three fluxes in 7.54 %). The sunlight the floor
sends from the surface to the reflection is 0.009 W/m² in the global mean:

| zenith angle | below 70 | 70-80 | 80-85 | 85-90 |
|---|---|---|---|---|
| global-mean sunlight moved from the surface to the reflection, W/m² | 0.000 | 0.000 | 0.002 | 0.007 |
| sunlit samples where the surface's absorption is floored | 0.0 % | 0.2 % | 2.1 % | 12.0 % |

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
| ASR | 242.4 | 237.4 | 240.1 | 234.0 | 243.7 | 240.9 |
| atmosphere | 83.4 | 78.7 | 82.0 | 77.0 | 79.6 | |
| cloud-induced atmosphere | 10.3 | 5.6 | 8.9 | 3.9 | 2.4 | |
| surface | 158.9 | 158.7 | 158.1 | 157.0 | 164.1 | |
| SWCRE | -50.5 | -55.5 | -52.8 | -58.9 | -46.7 | -45.8 |
| June 0-20N ASR | 312.3 | 302.1 | 307.8 | 291.2 | 326.5 | |
| June 0-20N SWCRE | -78.0 | -88.3 | -82.5 | -99.1 | -44.4 | -52.1 |
| June 40-60N ASR | 333.8 | 326.5 | 331.3 | 318.9 | 316.7 | |
| June 40-60N SWCRE | -95.1 | -102.5 | -97.7 | -110.1 | -99.3 | -90.4 |

The clear-sky ASR is 292.9 in every column. CERES's 286.7 is that of cloud-free
scenes, which differ from the clear sky of cloudy columns by 0.24 W/m² in the
global shortwave (Kato et al. 2013, as quoted by Loeb et al.), so the model's
clear sky absorbs 6.2 W/m² more (ERA5: 290.4, 2.5 less): the surface 219.7
against ERA5's 213.2, the atmosphere 73.1 against 77.2. The Kepler orbit
changes none of the differences between the columns by more than 0.05 W/m²
(annual means; June, near aphelion, is 3 % dimmer).

Shielding takes 5.9 W/m² of gas absorption, the clouds' own absorption rises
by 1.4, and the planet brightens by 5.0 W/m² (4.7 to 5.2 by season). The layered
scheme brightens it by 8.4 (7.9 to 9.5 by season): its cloud-induced
absorption, 3.9, is the closest of the three to ERA5's 2.4, and its reflection
is what a delta-Eddington cloud of the model's optical depths gives (the
single columns below). The tropics carry it: June 0-20N loses 21 W/m² of ASR
(its SWCRE is -99 against CERES's -52), the clouds of that state hold more
water than the sky does.

The shallow-cumulus cloud that the saved states hold for the radiation (the
lowest 18 layers' plume cover and water) is part of these numbers: left out it
would read 0.5 W/m² more absorbed sunlight in the sequential scheme (242.8) and
1.0 more in the layered one (235.0).

## Where the layered scheme's brightening sits

Layered less sequential, the same clouds and decks, over the 163,847 cell-times
of the four seasonal days, global means in W/m²: reflected +8.39 (98.14 to
106.54), absorbed -8.39, of which the atmosphere -6.45 and the surface -1.94.
By latitude band:

| band | area | reflected, band mean | atmosphere | surface | share of the global change |
|---|---|---|---|---|---|
| 90S-60S | 6.7 % | +2.75 | -2.42 | -0.32 | 2 % |
| 60S-30S | 18.3 % | +7.77 | -5.88 | -1.89 | 17 % |
| 30S-0 | 25.0 % | +8.47 | -6.06 | -2.41 | 25 % |
| 0-30N | 25.0 % | +11.56 | -8.97 | -2.60 | 34 % |
| 30N-60N | 18.3 % | +8.01 | -6.33 | -1.68 | 17 % |
| 60N-90N | 6.7 % | +4.67 | -4.46 | -0.21 | 4 % |

By cloud class, ISCCP-style: the pressure of the first cloud the light meets
(the level by optical depth) below 440 hPa, from 440 to 680 and above 680, and the
column's in-cloud visible optical depth (its grid-mean depth over its cover) from
0.3 to 3.6, from 3.6 to 23 and above 23. Shares of the sunlit samples, the
class's reflectance (reflected over incident sunlight) in the two schemes, and
the class's part of the global change:

| class | sunlit samples | reflectance, sequential / layered | change, W/m² | share |
|---|---|---|---|---|
| high, depth > 23 | 8.8 % | 0.479 / 0.643 | +4.93 | 59 % |
| high, 3.6-23 | 10.4 % | 0.381 / 0.414 | +1.14 | 14 % |
| high, 0.3-3.6 | 4.9 % | 0.224 / 0.230 | +0.09 | 1 % |
| mid, depth > 23 | 3.3 % | 0.448 / 0.510 | +0.73 | 9 % |
| mid, 3.6-23 | 10.4 % | 0.309 / 0.323 | +0.48 | 6 % |
| mid, 0.3-3.6 | 2.8 % | 0.332 / 0.334 | +0.01 | 0 % |
| low, depth > 23 | 8.1 % | 0.195 / 0.204 | +0.29 | 4 % |
| low, 3.6-23 | 36.3 % | 0.257 / 0.263 | +0.70 | 8 % |
| low, 0.3-3.6 | 4.1 % | 0.244 / 0.246 | +0.02 | 0 % |
| clear or depth < 0.3 | 10.9 % | 0.193 / 0.193 | 0.00 | 0 % |

The tropics carry 59 % of it and the towers 59 %: columns topped above 440 hPa
with a depth above 23, whose mean depth is 80 with 464 g/m² of liquid and 211 g/m²
of ice, reflect 0.48 of the sunlight in the sequential scheme and 0.64 in the
layered one, and 4.74 of their 4.93 W/m² comes out of the atmosphere's
absorption: the sequential scheme lets the gases below such a cloud take
sunlight that the cloud then does not reflect. The low clouds, 49 % of the sunlit
samples, change it by 1.0 W/m².

## The model's clouds against satellite retrievals and ERA5

The sequential scheme's absorbed sunlight is near CERES's (242.4 against 240.9)
with clouds that are optically thick. The tropical means (30S-30N) of the model's
four seasonal days against the ESA Cloud_cci L3C v2.0 monthly products of
MODIS-Aqua and AVHRR NOAA-18 (Stengel et al. 2017) and against ERA5 (the cloud
liquid and ice water content of ARCO-ERA5 at 12 UTC on 20 March, June, September
and December 2020, 0.25 degrees, with the model's optics laid on them):

| 30S-30N | model | ERA5 | MODIS-Aqua | AVHRR NOAA-18 |
|---|---|---|---|---|
| cloud cover | 53.1 % | | 58.0 % | 54.9 % |
| cover by top: high, mid, low | 17.6, 8.9, 21.8 % (5.0 % below depth 0.3) | | 22.6, 10.0, 24.5 % | 24.4, 8.2, 21.8 % |
| condensate path, g/m²: liquid + ice | 80 + 22 = 102 | 55 + 16 = 71 | 25 + 45 = 70 | 22 + 87 = 109 |
| the same in the retrievals' convention (a column topped by ice has its whole path in ice) | 39 + 63 | 36 + 36 | as retrieved | as retrieved |
| grid-mean optical depth | 12.8 | 8.6 | 4.0 | 4.8 |
| in-cloud optical depth: mean (geometric mean) | 26.5 (14.2) | | 6.9 (3.6) | 8.8 (4.2) |
| over ocean | 22.1 (13.1) | | 6.4 (3.2) | 8.1 (3.8) |
| over land | 40.0 (18.3) | | 8.8 (5.1) | 11.2 (5.5) |
| share of the depth in liquid | 88 % | 87 % | | |
| liquid droplet radius over ocean, land, um | 11.8, 8.5 | | 16.3, 11.7 | 13.9, 10.5 |
| area topped above 440 hPa with depth > 23 (ERA5: grid-mean depth > 23) | 7.8 % | 7.9 % | 1.5 % | |
| global: in-cloud optical depth, mean (geometric mean) | 22.6 (11.7) | | 10.6 (5.4) | 11.7 (5.8) |
| global: area topped above 440 hPa with depth > 23 | 8.6 % | | 2.6 % | |

The model: cover, top pressure and in-cloud depth as above, weighted by area and
cover over the columns with cover above 0.05; water paths are the column's
(liquid and ice by the radiation's temperature rule, resolved plus the cumulus
plume's water, below cloud tops too); the droplet radii are the options' (11.8
over sea, 8.5 over land). The retrievals see the topmost layer in daytime at the
overpass (13:30 for Aqua, about 14:00 for NOAA-18): cloud fractions and
all-sky paths are means of March, June, September and December 2007 to 2009, the
optical depths and radii (weighted by area and cloud fraction, by liquid or ice
share for the phase variables) of the same months of 2008, and the depth-pressure
frequencies the daytime joint histogram of 2008 on every sixth 0.5-degree cell,
area-weighted. Grid-mean depth is the in-cloud mean times the cover. The
retrievals' paths are their liquid and ice retrievals' all-sky means; the model's
and ERA5's are the water of the whole column, which the retrievals do not see
below the top of a thick cloud.

Sources (ESA Cloud_cci, anonymous OPeNDAP at CEDA, `dap.ceda.ac.uk/thredds/dodsC/
neodc/esacci/cloud/data/phase-2/L3C/`; ARCO-ERA5 on Google Cloud Storage,
`gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3`; no login),
each set in a directory of its own:

| set | source path | files | size |
|---|---|---|---|
| cloud fraction (total, high, mid, low), all-sky liquid and ice water path, 2007-2009, four months | `MODIS-AQUA/v2.0/<year>/<yyyymm>-ESACCI-L3C_CLOUD-CLD_PRODUCTS-MODIS_AQUA-fv2.0.nc` | 73 | 74.9 MB |
| the same, nine months (2007, 2008 and March 2009: the product has no June 2009) | `AVHRR-PM/v2.0/AVHRR_NOAA-18/<year>/<yyyymm>-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-18-fv2.0.nc` | 55 | 56.2 MB |
| cot, cot_log, cer_liq, cer_ice, cot_liq, cot_ice, cph, 2008, four months | MODIS-AQUA | 29 | 29.1 MB |
| the same | AVHRR NOAA-18 | 29 | 29.1 MB |
| hist2d_cot_ctp and nobs_day, every sixth 0.5-degree cell, 2008, four months | MODIS-AQUA | 9 | 45.1 MB |
| specific_cloud_liquid_water_content and specific_cloud_ice_water_content on 37 levels and total_column_cloud_liquid_water and total_column_cloud_ice_water, 30S-30N, four dates of 2020 (the Blosc blocks of the band, read by byte ranges from the global chunks) | `specific_cloud_liquid_water_content/<hours since 1900>.0.0.0` and the others | 320 | 79.6 MB |

## Why the tropical clouds are three times too thick

The model's tropical grid-mean optical depth, 12.8, is 3.2 times MODIS's 4.0
(2.7 times AVHRR's 4.8) where the cover is within 5 points of theirs. Three
factors make the 3.2:

- 1.5 is the model's condensate against ERA5's: 102 g/m² against 71, and 12.8
  against 8.6 with the same optics.
- 1.4 is the droplet size: with the retrievals' radii (16.3 um over ocean, 11.7
  over land, where the model has 11.8 and 8.5; the ocean's MODIS value is itself
  about 4 um above the MISR-referenced 11.2 to 12.1 of Fu et al. 2019) ERA5's depth
  falls to 6.3.
- 1.6 is what is left of ERA5 against MODIS with the same optics: passive
  retrievals see the top of a thick cloud and give an ice-topped column's whole
  depth to ice (Platnick et al. 2017), the overpass is the afternoon, and ERA5's
  own error.

The phase is not the cause. Liquid carries 88 % of the model's depth and 87 % of
ERA5's, and in the retrievals' convention the model's paths (39 liquid, 63 ice)
lie between AVHRR's and MODIS's ice and 1.5 to 1.8 times their liquid. The
model's excess over ERA5 is in two places. By layer, in g/m² of the tropical
mean:

| layer | model | ERA5 | model / ERA5 |
|---|---|---|---|
| above 250 hPa | 5.6 | 4.7 | 1.2 |
| 250-450 hPa (-38 to -15 C) | 25.6 (liquid 10.2) | 9.8 (liquid 0.9) | 2.6 |
| 450-700 hPa | 24.2 | 22.4 | 1.1 |
| 700-1000 hPa | 46.6 (the plume's 8.6) | 34.8 | 1.3 |
| 250-1000 hPa, resolved cloud only (the plume's 9.7 left out) | 86.8 | 67.0 | 1.3 |

Below 700 hPa the resolved cloud (38.0 g/m² of liquid at 700-1000 hPa) is ERA5's
34.9 within 9 %, and the cumulus plume's water (9.7 g/m², 10 % of the in-cloud
depth), which the radiation adds to it, is the excess. In the mixed-phase layer
the model holds 2.6 times ERA5's condensate and eleven times its liquid. By
temperature band the in-cloud depth of 26.5 is 66 % warm liquid, 20 % in
the band 0 to -20 C, 8 % in -20 to -40 C and 5 % in ice colder than -40 C; the
condensate there is 60, 22, 14 and 6 g/m² (ERA5: 45, 14, 7 and 5). The towers,
columns topped above 440 hPa with a depth above 23, are about 8 % of the area with
a mean depth of 91, 83 % of it liquid and 4 % the plume's, and 744 g/m² of water.
The liquid share of the condensate mass by layer temperature, against the
liquid share of the cloud tops the retrievals see:

| temperature, C | -36 | -31 | -26 | -21 | -16 | -11 | -6 | -1 |
|---|---|---|---|---|---|---|---|---|
| model, the radiation's ramp (273.15 to 235.15 K) | 6 % | 20 % | 31 % | 47 % | 58 % | 71 % | 87 % | 98 % |
| ERA5, share of the mass | 0 % | 2 % | | 8 % | 30 % | 52 % | 81 % | 94 % |
| MODIS-Aqua, share of the cloud tops (daytime, 30S-30N, 2008) | 6 % (-38) | | 33 % (-27) | | 64 % | | 84 % (-8) | 100 % (+1) |

The ramp is the retrievals' cloud-top phase within 6 points, and it gives
the layer 250-450 hPa eleven times ERA5's liquid mass; the two measure different
things, the top layer and the mass, and ERA5's partition is itself a model's.

The processes that would have to change, each with the sunlight the layered
scheme gives (annual global absorbed sunlight and SWCRE in W/m²; CERES 240.9 and
-45.8) when the tropical cloud (30S-30N) or the rule is taken to ERA5's, nothing
else altered:

| process | evidence | reference | layered | sequential |
|---|---|---|---|---|
| the model as is | | | 234.0 / -58.9 | 242.4 / -50.5 |
| mixed-phase condensate, 250-450 hPa: glaciation (Bergeron-Findeisen, riming), conversion to precipitation, ice fall speed | 2.6 times ERA5's mass | Korolev et al. 2017; Tiedtke 1993; Heymsfield and Donner 1990 | 235.2 / -57.6 | 242.9 / -50.0 |
| anvil condensate above 250 hPa: detrainment, lifetime | 1.2 times | Tiedtke 1993; Li et al. 2012; Waliser et al. 2009 | 234.3 / -58.5 | 242.6 / -50.2 |
| shallow-cumulus plume water added to the resolved cloud, 9.7 g/m² | 10 % of the depth | Siebesma et al. 2003 | 235.0 / -57.9 | 242.8 / -50.0 |
| the cloud below 450 hPa as a whole (the plume's water, the resolved excess at 650-750 and 950-1000 hPa) | 1.24 times | Khairoutdinov and Kogan 2000; Siebesma et al. 2003 | 237.3 / -55.5 | 245.1 / -47.7 |
| phase partition: the ramp taken to ERA5's mass shares (271.15 to 251.15 K, everywhere) | 47 % against 8 % liquid at -21 C | Korolev et al. 2017; Hu et al. 2010; ECMWF IFS Documentation Part IV | 234.9 / -57.9 | 242.9 / -50.0 |
| the tropical profile at ERA5's, bin by bin, reductions only | 71 g/m² against 102 | | 239.1 / -53.7 | 246.1 / -46.8 |
| the same with ERA5's phase partition | | | 240.0 / -52.9 | 246.6 / -46.3 |
| autoconversion of the warm cloud (Kessler, 1 h) | the resolved liquid at 700-1000 hPa is ERA5's within 9 % | Kessler 1969; Khairoutdinov and Kogan 2000 | not implicated | |
| overlap, the two bounds (random, maximum) | cover 53 % against 55 to 58 % | Hogan and Illingworth 2000; Shonk et al. 2010 | 229.1 / -63.8; 235.7 / -57.2 | 237.6 / -55.3; 244.0 / -48.9 |

The tropical condensate at ERA5's takes the layered scheme to within 1.8 W/m² of
CERES's absorbed sunlight (0.9 with ERA5's phase partition), and its cloud
effect from -58.9 to -53.7 (-52.9). The layer 250-450 hPa holds half of the
excess condensate and gives a quarter of the sunlight, because three quarters of
a column's depth lies below it; the cloud below 450 hPa gives two thirds. What
remains of the cloud effect, 7 W/m² too strong, is outside what this changes: the
clouds outside 30S-30N (June 40-60N is -110 against CERES's -90), the optics and
the cover. The clear sky's 6.2 W/m² excess absorption is separate; with the
cloud effect at CERES's -45.8 the layered scheme would absorb 247.

## What else would close the gap to CERES (none applied)

The reflection is the cloud's: its water (above), its optical depths (droplet
radius 11.8 and 8.5 um, Han et al. 1994; MODIS 2.1 um's ocean mean 16.3 to 16.6 um
falls to 11.2 to 12.1 once its bias against MISR is removed, Fu et al. 2019), its
inhomogeneity (cloudSpread 0.75, 0.75 +- 0.18, Shonk et al. 2010) and its cover.
Absorbed sunlight and SWCRE of the four seasonal days with the clouds or their
optics changed offline, in W/m²:

| change | sequential | cloudTop | layered |
|---|---|---|---|
| none | 242.4 / -50.5 | 237.4 / -55.5 | 234.0 / -58.9 |
| cloud water above 440 hPa x 0.7 | 243.4 / -49.4 | | 235.8 / -57.1 |
| cloud water above 440 hPa x 0.5 | 244.3 / -48.6 | | 237.2 / -55.7 |
| cloud water above 440 hPa x 0.3 | 245.5 / -47.4 | | 239.0 / -53.9 |
| cloud water above 440 hPa x 0.15 | 246.6 / -46.3 | | 240.6 / -52.2 |
| cloud water at every level x 0.7 | 250.1 / -42.8 | | 243.8 / -49.1 |
| cloud water at every level x 0.5 | 256.8 / -36.0 | | 252.1 / -40.8 |
| droplet radii +20 % (sea 14.2, land 10.2 um) | | 241.3 / -51.6 | 238.2 / -54.6 |
| cloudSpread 1.0 | | 238.7 / -54.1 | 235.2 / -57.6 |
| both | | 242.5 / -50.4 | 239.4 / -53.5 |
| liquid asymmetry +0.02 (Slingo's 0.846 at 6.9 um; Mie about 0.86) | | | 236.5 / -56.4 |
| cloudSolarAbsorption 0.55 | 244.8 / -48.1 | 240.1 / -52.8 | |

The layered scheme reaches CERES's absorbed sunlight with the cloud water above
440 hPa at about 0.12 of its present amount, with the water at every level at
about 0.8, or with the droplet radii and the inhomogeneity both at the edge of
their observed ranges (239.4); its SWCRE of -45.8 needs the water at every level
at about 0.6, with the clear sky's 6.2 W/m² excess as it is. The same changes of
the water put the sequential scheme 3 to 16 W/m² above CERES.

## Mixed-phase liquid at 250-450 hPa: `moist.bergeronAutoconversion`

The layer 250-450 hPa (-38 to -15 C) holds 25.6 g/m² of condensate in the
tropical mean against ERA5's 9.8, and 10.2 g/m² of it as liquid against about 1.
ERA5's cloud ice excludes snow, which the model lumps with its ice, so the
liquid is the like-for-like number. Four things keep it there.

- One phase rule. `liquidFraction` (linear in temperature from 273.15 K to
  235.15 K, 44 % liquid at -21 C) splits the condensate for the saturation mix,
  the plume's liquid-ice static energy, the Kessler conversion, the melting of
  falling ice and the radiation; the moist module takes the radiation's two
  temperatures unless its own options name others. Nothing in the resolved cloud
  turns liquid into ice faster than the temperature does.
- No Bergeron-Findeisen or riming term in the resolved cloud's conversion. The
  factor 1 + 0.5 sqrt(268.16 - T) appears only in the plume's Sundqvist
  conversion; the liquid share of the resolved cloud converts at Kessler's rate
  above 0.2 g/kg (1e-3 per s) and over the cloud lifetime (1 h, 3 h under an
  inversion) at -20 C as at +20 C.
- One fall speed for the ice share, Heymsfield and Donner's form, 2.5 m/s times the
  in-cloud ice water content (kg/m³) to the power 0.16 (0.4 to 0.8 m/s), for the
  cloud ice and for the snow the model lumps with it.
- The plume's detrained water takes the same ramp.

One pass of the moist physics on the tropical mean of fourteen64's June state:
the layer holds 27.4 g/m², 12.1 of it liquid, and the conversion and fall stage
takes 24.4 g/m² per hour out of it, half Kessler's conversion of the liquid and
half the ice's settling, a residence of 1.1 to 1.3 h against 0.7 h at 450-700 hPa
and 1.4 h above 250 hPa.

Six hours of the CPU model from that state (N = 64, 64 steps), the tropical
(30S-30N) mean in g/m² of condensate by layer, the liquid of 250-450 hPa in
parentheses, and the rain, with the moist options named changed:

| options | 250-450 hPa | above 250 | 450-700 | 700-1000 | column | rain, mm/day |
|---|---|---|---|---|---|---|
| none (the start: 27.4, 5.5, 26.2, 43.5, 102.7) | 27.7 (12.3) | 5.0 | 26.7 | 42.8 | 102.3 | 3.29 |
| `upperCloudLifetime` 1200 s | 18.2 (7.5) | 4.9 | 16.5 | 43.6 | 83.1 | 3.54 |
| `iceFall` 5 | 20.0 (9.7) | 2.0 | 26.3 | 43.0 | 91.3 | 3.36 |
| `cloudLifetime` 1800 s | 21.8 (9.3) | 4.9 | 20.2 | 32.0 | 78.9 | 3.50 |
| the moist ramp 271.15 K to 251.15 K (the radiation's unchanged) | 34.4 (7.8) | 5.3 | 26.9 | 42.7 | 109.3 | 3.39 |
| the three together | 17.5 (4.1) | 2.1 | 16.4 | 43.3 | 79.3 | 3.68 |
| `bergeronAutoconversion` | 16.4 (6.7) | 4.8 | 26.0 | 43.0 | 90.2 | 3.40 |
| `bergeronAutoconversion` and the glaciated ramp | 26.8 (4.8) | 5.2 | 26.0 | 42.7 | 100.7 | 3.45 |
| `bergeronAutoconversion` and `iceFall` 4 | 13.9 (6.1) | 2.6 | 25.8 | 43.1 | 85.4 | 3.43 |

A shorter lifetime or a faster fall acts on every layer that has the water: the
lifetime on 450-700 hPa and the fall on the anvils, whose 5.0 g/m² is ERA5's 4.5.
A glaciated ramp lowers the liquid and raises the condensate, because more of the
vapor condenses over ice. The Bergeron-Findeisen factor is the one change that
acts only where the temperature is below 268 K (450-700 hPa loses 3 %, 700-1000 hPa
nothing), and it is a process the cloud has in the IFS and the plume here has.

`bergeronAutoconversion` (false) multiplies the liquid share's two conversion
rates and divides its threshold, below 268.16 K, by F = 1 + 0.5 sqrt(min(268.16 -
T, 18)): F is 1.9 at 265 K, 2.4 at 260 K, 3.0 at 252 K (-21 C) and 3.1 from 250.16
K to the colder end of the liquid ramp. The IFS's large-scale precipitation does
this to the two parameters of Sundqvist's conversion, the rate c0 and the
critical water content l_crit, as c0 = c0* F1 F2 and l_crit = l_crit* / (F1 F2)
with F2 = 1 + b2 sqrt(T_BF - T) for T_ice < T < T_BF and Sundqvist's constants,
T_BF = 268 K and b2 = 0.5 (IFS Documentation CY25r1 Part IV chapter 6, equations
6.49 to 6.53; F1 is the collection by falling rain, which this model's conversion
does not have). Here Kessler's rate and threshold take the places of c0 and
l_crit; the factor of the plume's conversion, which carries the later cycles'
268.16 K and 250.16 K, is the one used, and below 250.16 K the factor stays at its
value there, because the liquid ramp reaches down to 235.15 K where the IFS's
cloud is ice. The option changes the liquid share's conversion and nothing else:
the phase ramp, the radiation, the ice's fall and the plume are as they were, the
shader text is unchanged with the option off, and the default path's digests hold.
The converted liquid leaves as rain without freezing aloft or melting at 0 C; at
+0.1 mm/day in the tropical mean (the table's 3.29 to 3.40) that is 0.4 W/m² of
fusion heat (3.34e5 J/kg times 0.1 kg/m² per day over 86400 s), which a rule that
froze the water where it forms and melted the snow would move. It is not built.

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

## Cheaper layered variants (offline, not options)

The layered work is the number of slices (the g-points and the two visible
slices) times the regions each solves (the clear one and one per lognormal node)
times the layers. The same four days with the table folded to fewer g-points
(`POINTS=24 WRITE=1 node scripts/layeredFit.mjs` in a scratch tree) and with
two nodes in place of four:

| g-points | nodes | work (34 x 5 = 1) | direct beam absorption miss, W/m² of 1360 (mean, worst) | global ASR | SWCRE |
|---|---|---|---|---|---|
| 32 | 4 | 1.00 | 1.2, 4.1 | 234.0 | -58.9 |
| 24 | 4 | 0.76 | 2.0, 7.3 | 233.9 | -59.0 |
| 16 | 4 | 0.53 | 5.5, 29.9 | 233.8 | -59.1 |
| 32 | 2 | 0.60 | 1.2, 4.1 | 232.7 | -60.2 |
| 24 | 2 | 0.46 | 2.0, 7.3 | 232.6 | -60.3 |
| 16 | 2 | 0.32 | 5.5, 29.9 | 232.4 | -60.4 |

The global balance hardly notices the g-points (0.1 to 0.2 W/m²) and takes
1.3 W/m² less sunlight with two nodes in place of four; the clear-sky
absorption of the folded table does notice them at 16 g-points.

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

## References for the cloud and low-sun sections

- Barker, H. W., J. N. S. Cole, J. Li, B. Yi and P. Yang, 2015: Estimation of
  errors in two-stream approximations of the solar radiative transfer equation
  for cloudy-sky conditions. J. Atmos. Sci., 72, 4053-4074.
- ECMWF, IFS Documentation, Part IV: Physical Processes, chapter 7, Clouds and
  large-scale precipitation.
- ECMWF, 2003: IFS Documentation Cycle CY25r1, Part IV: Physical processes,
  chapter 6, Clouds and large-scale precipitation (equations 6.49 to 6.53).
- Heymsfield, A. J., and L. J. Donner, 1990: A scheme for parameterizing ice-cloud
  water content in general circulation models. J. Atmos. Sci., 47, 1865-1877.
- Hogan, R. J., and A. J. Illingworth, 2000: Deriving cloud overlap statistics
  from radar. Q. J. R. Meteorol. Soc., 126, 2903-2909.
- Hu, Y., et al., 2010: Occurrence, liquid water content, and fraction of
  supercooled water clouds from combined CALIOP/IIR/MODIS measurements. J.
  Geophys. Res., 115, D00H34.
- Kessler, E., 1969: On the distribution and continuity of water substance in
  atmospheric circulations. Meteor. Monogr., 10 (32).
- Khairoutdinov, M., and Y. Kogan, 2000: A new cloud physics parameterization in
  a large-eddy simulation model of marine stratocumulus. Mon. Wea. Rev., 128,
  229-243.
- King, M. D., and Harshvardhan, 1986: Comparative accuracy of selected multiple
  scattering approximations. J. Atmos. Sci., 43, 784-801.
- Korolev, A., et al., 2017: Mixed-phase clouds: Progress and challenges. Meteor.
  Monogr., 58, 5.1-5.50.
- Li, J.-L. F., et al., 2012: An observationally based evaluation of cloud ice
  water in CMIP3 and CMIP5 GCMs and contemporary reanalyses using contemporary
  satellite data. J. Geophys. Res., 117, D16105.
- Platnick, S., et al., 2017: The MODIS cloud optical and microphysical products:
  Collection 6 updates and examples from Terra and Aqua. IEEE Trans. Geosci.
  Remote Sens., 55, 502-525.
- Shonk, J. K. P., R. J. Hogan, J. M. Edwards and G. G. Mace, 2010: Effect of
  improving representation of horizontal and vertical cloud structure on the
  Earth's global radiation budget. Part I: Review and parametrization. Q. J. R.
  Meteorol. Soc., 136, 1191-1204.
- Siebesma, A. P., et al., 2003: A large eddy simulation intercomparison study of
  shallow cumulus convection. J. Atmos. Sci., 60, 1201-1219.
- Stengel, M., et al., 2017: Cloud property datasets retrieved from AVHRR, MODIS,
  AATSR and MERIS in the framework of the Cloud_cci project. Earth Syst. Sci.
  Data, 9, 881-904.
- Sundqvist, H., 1978: A parameterization scheme for non-convective condensation
  including prediction of cloud water content. Q. J. R. Meteorol. Soc., 104,
  677-690.
- Tiedtke, M., 1993: Representation of clouds in large-scale models. Mon. Wea.
  Rev., 121, 3040-3061.
- Waliser, D. E., et al., 2009: Cloud ice: A climate model challenge with signs
  and expectations of progress. J. Geophys. Res., 114, D00A21.
