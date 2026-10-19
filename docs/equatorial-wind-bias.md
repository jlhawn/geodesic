# The equatorial Pacific wind bias

Why the coupled model's trade winds over the equatorial Pacific are 1.5 to
3 times too strong, measured on 2026-10-18 and 19 after the hybrid ocean's
gate campaign (ocean-hybrid-coordinate.md, row 5) had shown the ocean's
equatorial state to be the forced response to those winds.

## The question

The twelve-local recording (N=128, year 5, days 1461–1825) carries a
2S–2N zonal stress of −0.06 to −0.09 N/m² over 160E–100W from June to
October and −0.08 to −0.12 over the warm pool, against about −0.04 to
−0.055 and −0.02 observed (QuikSCAT and TAO climatologies, approximate).
Replaying the recording on the ocean at 0.6 of its stress puts the
tropics, the 140W thermocline and the undercurrent on the World Ocean
Atlas; nothing in the ocean moves them. The stress is the atmosphere's,
and the question is which part of the atmosphere.

## Instruments

- `scripts/equatorialMomentum.mjs STATE`: the boundary-layer zonal
  momentum budget on the equator by 20° of longitude from one saved state
  (CPU, one step, ocean off), with the pressure-gradient force by layer,
  the layer temperatures along the equator and the sea-level pressure
  difference from 120W–100W split by height band (the hydrostatic weight
  of each band, so that the gradient can be charged to the layers that
  carry it).
- `scripts/fixedSst.mjs STATE OUT TAG`: the atmosphere alone on the GPU
  from a saved state over a sea surface held fixed (the ocean off, a slab
  of 10¹² J/m²/K): `SST=own` keeps the state's sea surface, `SST=woa`
  sets every open-sea cell to WOA23's monthly SST on the state's date.
  `MOIST`, `BOUNDARY_LAYER`, `SURFACE` and `RADIATION` (JSON) pass
  options; `DAYS`, `SAVE_FROM` and `SAVE_DAYS` choose the saves. The AMIP
  question in twenty days at 21 s a day at N=128 on the M5 Max.
- `scripts/armMetrics.mjs LABEL STATE...`: the same equatorial bins
  averaged over several saved states (sea-level pressure, lowest-layer
  wind and stress, layer temperatures, relative humidity, the virtual
  temperature and the hydrostatic weight of each height band against
  120W–100W, the deck's height, the cumulus top, and the mean column of
  180–160W and 120W–100W).

All runs start from `runs/twelve-local/twelve_day1597.bin` (August 4 of
year 5, when the bias is largest). Single states vary by 30 to 50 percent
in the stress from one day to the next (easterly waves and the diurnal
cycle), so the arms of section 4 run twenty days and are judged on the
mean of five states two days apart, days 1609 to 1617 (August 16 to 24).

## Findings

### 1. Over the warm pool the excess is the sea surface, not the atmosphere

Twelve days from the same state, the sea surface held at the state's own
(the coupled run's, 3 K cold over the warm pool) and at WOA23's:

| 2S–2N, day +12 | 140E–160E | 160E–180 | 180–160W | 160W–140W | 140W–120W | 120W–100W | 100W–80W |
|---|---|---|---|---|---|---|---|
| SST held at the state's (°C) | 25.7 | 24.6 | 22.1 | 21.5 | 21.2 | 20.6 | 21.6 |
| stress τx (N/m²) | −0.084 | −0.134 | −0.098 | −0.093 | −0.103 | −0.080 | −0.019 |
| u lowest layer (m/s) | −8.1 | −9.8 | −8.6 | −8.2 | −7.9 | −6.2 | −1.8 |
| SST held at WOA23 (°C) | 29.3 | 29.3 | 28.4 | 27.2 | 25.2 | 23.7 | 23.4 |
| stress τx (N/m²) | −0.028 | −0.007 | −0.058 | −0.089 | −0.103 | −0.091 | −0.020 |
| u lowest layer (m/s) | −3.0 | −1.9 | −6.5 | −8.3 | −8.5 | −7.1 | −1.8 |
| observed stress, about | −0.015 | −0.02 | −0.04 | −0.05 | −0.055 | −0.05 | −0.03 |
| observed 10 m wind, about | −2 | −3 | −5 | −6 | −6 | −5 | −3 |

West of the date line the coupled bias is the Bjerknes loop on the
ocean's cold warm pool: the atmosphere over the observed sea surface
makes the observed winds there within six days.

### 2. Over the central and eastern Pacific the atmosphere alone is 1.3 to 2 times too strong

The time mean over WOA23's sea surface (the base of section 4, days
1609–1617, `plumeMomentum` on):

| 3S–3N sea, mean of days 1609–1617 | 140E–160E | 160E–180 | 180–160W | 160W–140W | 140W–120W | 120W–100W | 100W–80W |
|---|---|---|---|---|---|---|---|
| sea-level pressure (hPa) | 1007.1 | 1007.0 | 1007.4 | 1008.8 | 1010.5 | 1011.9 | 1011.9 |
| u lowest layer (m/s) | −2.0 | −4.7 | −8.3 | −8.8 | −7.2 | −5.9 | −1.5 |
| stress τx (mN/m²) | −18 | −42 | −100 | −113 | −83 | −66 | −16 |
| observed stress, about (mN/m²) | −15 | −20 | −40 | −50 | −55 | −50 | −30 |
| observed sea-level pressure, about (hPa) | 1008.5 | 1009 | 1010 | 1011 | 1011.5 | 1012 | 1012.5 |

The sea-level pressure rises by 4.5 hPa from 180–160W to 120W–100W where
about 2 hPa is observed in August (ICOADS and ERA climatologies,
approximate: 1010 at 170W, 1011.3 at 140W, 1012 at 110W). The
pressure-gradient force across the central Pacific is about twice the
observed, the wind 1.2 to 1.45 times (u² follows the gradient), the
stress 1.3 to 2.2 times, largest at 160W–140W. The model's trough is too
flat over the warm pool and the date line (0.3 hPa from 150E to 170W
where 1.5 is observed) and too steep east of it.

Where the gradient sits, hydrostatically. The sea-level pressure of each
bin less that of 120W–100W, split by the height band that carries it
(the integral of g ρ over the band, interpolated in ln ρ; the remainder is
the air above 18 km and the interpolation):

| band | 140E–160E | 160E–180 | 180–160W | 160W–140W | 140W–120W |
|---|---|---|---|---|---|
| 0–1 km | −2.48 | −2.56 | −2.22 | −1.52 | −0.71 |
| 1–2 km | −2.06 | −2.05 | −1.66 | −1.07 | −0.46 |
| 2–3 km | −0.82 | −0.77 | −0.55 | −0.25 | −0.06 |
| 3–5 km | −0.44 | −0.44 | −0.40 | −0.17 | −0.03 |
| 5–8 km | −0.24 | −0.22 | −0.25 | −0.12 | 0.00 |
| 8–12 km | −0.99 | −0.74 | −0.58 | −0.38 | −0.03 |
| 12–18 km | +1.61 | +1.33 | +0.82 | +0.16 | −0.12 |
| sea-level pressure difference | −4.8 | −4.9 | −4.5 | −3.1 | −1.4 |

Eighty-seven percent of the model's gradient between 180–160W and
120W–100W is the lowest two kilometers; 3 to 12 km add 1.2 hPa more,
and the tropopause layer takes 0.8 back. The virtual temperature of
each bin less 120W–100W's at fixed heights, with the observed contrast
(ERA climatology along the equator, approximate) for 180–160W:

| height | 140E–160E | 180–160W | 160W–140W | 140W–120W | observed 180–160W, about |
|---|---|---|---|---|---|
| 500 m | 5.6 | 4.9 | 3.4 | 1.6 | 4 |
| 1000 m | 6.0 | 4.8 | 3.2 | 1.5 | 4 |
| 1500 m | 6.2 | 4.9 | 3.2 | 1.3 | 3 |
| 2000 m | 4.0 | 2.6 | 1.3 | 0.5 | 2 |
| 3000 m | 1.5 | 1.0 | 0.4 | −0.2 | 1 |
| 4000 m | 1.1 | 0.8 | 0.3 | 0.1 | ≤ 0.5 |
| 9000 m | 3.0 | 2.0 | 1.1 | −0.1 | ≤ 0.5 |
| 12000 m | 2.5 | 1.0 | 0.2 | 0.0 | ≤ 0.5 |

The excess has two parts of similar size:

1. The lowest 2.5 km over the cold tongue are 2 to 3 K too cold at 1 to
   2 km (about 1 hPa). The eastern column in the time mean is a dry
   adiabatic mixed layer to 0.8 km (23.4 °C at 21 m, 16.7 at 791 m), a
   saturated cloud layer to 1.7 km (14.7 °C at 1.08 km, 12.5 at 1.50 km,
   relative humidity 0.92 to 0.95), then an isothermal layer to 2.1 km
   (12.4 °C, relative humidity 0.67) and a stable layer to 3 km (10.3 °C
   at 3.07 km, relative humidity 0.40). The observed cold-tongue column
   has its inversion at about 850 hPa, 1.5 km (Yin and Albrecht 2000,
   the FGGE soundings), with the free troposphere above it at 15 to 16 °C
   at 1.5 km and 14 to 15 at 2 km. The model's inversion lies 200 to 400 m
   higher and is smeared over the 1.8 to 2.5 km layer (the layers there
   are 450 to 600 m thick), so the 1.3 to 2.5 km column is cloud-layer or
   inversion air, 2 to 3 K colder than the free troposphere that should
   be there. The mixed-layer deck's time-mean height is 1.65 to 1.76 km
   over 140W–100W (its ceiling rule lets it stand 1 m below the midpoint
   of the layer above the inversion interface, 2.1 km) and the boundary
   layer mixes momentum to 1.5 to 1.9 km.
2. The free troposphere over the central Pacific is 1 to 2 K warmer than
   over the east at 3 to 12 km (about 1.2 hPa). The weak-temperature-
   gradient balance holds the real tropical free troposphere within about
   0.5 K along the equator; the model's 9 km level is 2 K warmer at
   180–160W and 3 K at 140E–160E than at 120W–100W, and the surplus
   reaches down to 3 km (1.0 and 1.5 K). It is not local convection (the
   3S–3N rain east of the date line is under 0.1 mm/d); it is the
   model's Walker and Hadley structure, the warm column of the convective
   regions extending too far east along the equator, and the twelve
   recording's rain (the coupled run) peaks at 140W–120W in the 5N–12N
   band (11.4 mm/d against 6 to 7 observed) and is weak at 120W–100W
   (3.2 against 7.5).

The budget closes the same way in every arm: in the boundary layer the
surface drag (−τx over the layer's mass, −7 to −8 × 10⁻⁵ m/s²) balances
the pressure-gradient force (−2 to −7), the relative-vorticity flux of
the cross-equatorial southerlies (−1 to −5, largest in the east where
v = 4 to 7 m/s carries the stronger south-east trades across the
equator) and the mixing-down of the 10 to 15 m/s easterly jet that sits
between 1 and 2 km.

What cools the eastern column between 1.5 and 3 km (the process budget
of `scripts/tropicalHeating.mjs` on the 3S–3N 180–100W box, K/day, one
state): longwave −2.3 to −3.0, the shallow plume's detrainment −0.9 to
−1.8, rain evaporation −0.4 to −0.6, the boundary layer's mixing −1.5 at
1.5 km and −0.35 at 2.1 km, the plume's downdraft −0.4 to −0.6, against
subsidence warming of +0.9 to +1.9 and the condensation of the deep
plume, which fires as a congestus (top 600–700 hPa, CAPE 42 J/kg) on
half the column-steps there and rains 2 mm/d. The longwave cooling of
2.3 to 3.0 K/day at 850 to 700 hPa is the emission of the moist layer's
top (relative humidity 0.67 at 2.1 km where about 0.4 is observed).

### 3. Ruled out

- The drag law: C_D is neutral COARE over the cold tongue (Ri_b ≈ 0), the
  stable functions engage only below 17.5 °C (earlier finding).
- The ocean: twelve ocean-only arms over every column, coordinate, eddy,
  tracer and closure knob leave the equatorial state unchanged
  (ocean-hybrid-coordinate.md, row 5).
- The depth of the boundary layer and of the deck as the lever: holding
  the deck at the Richardson depth (0.75 km instead of 1.76 over the cold
  tongue, the momentum mixing at 0.85 km instead of 1.8) leaves the stress
  and the pressure gradient unchanged (`noprog` below). The wind is not
  √(pgf · h / C_D) with h the mixing depth; the gradient is.
- The shallow plume's reach: its detrainment rate (4.5 and 6 × 10⁻³ /m
  against 3), its cap (800 and 850 hPa against 700) and the deck's
  entrainment efficiency leave the eastern 1.5 to 2.1 km layer at 12.4 to
  13.5 °C and the stress within 10 percent.
- The plume's momentum transport as the cause: `plumeMomentum` (on since
  this diagnosis) takes 18 percent off the time-mean stress over
  160W–100W (−106 to −87 mN/m², the last row of the table below) and is
  kept on as a sink the model lacked, not as the cause.

### 4. Knob arms over WOA23's sea surface

Twenty days from `twelve_day1597.bin` over WOA23's August sea surface,
judged on the mean of days 1609, 1611, 1613, 1615 and 1617; the base is
the model's defaults (`plumeMomentum` on). Observed, about: the
sea-level pressure difference −2 hPa, the stress −50, −55, −50 mN/m²,
the wind −6, −6, −5.5 m/s at 20 m.

| arm | SLP 180–160W less 120W–100W (hPa) | τx 160W–140W | 140W–120W | 120W–100W (mN/m²) | u 160W–140W | 140W–120W | 120W–100W (m/s) | deck h 120W–100W (m) | T 1.5 km east (°C) | T 2.1 km east |
|---|---|---|---|---|---|---|---|---|---|---|
| base | −4.49 | −112 | −83 | −66 | −8.8 | −7.2 | −5.9 | 1764 | 12.5 | 12.4 |
| cumulusDetrainment 4.5e-3 | −4.58 | −115 | −96 | −70 | −8.9 | −7.9 | −6.2 | 1850 | 12.7 | 12.5 |
| cumulusDetrainment 6e-3 | −4.58 | −112 | −88 | −69 | −8.8 | −7.6 | −6.1 | 1676 | 12.6 | 12.7 |
| cumulusEntrainment 5e-3 | −4.00 | −85 | −84 | −65 | −7.6 | −7.3 | −5.7 | 1642 | 12.4 | 12.7 |
| shallowTop 800 hPa | −3.88 | −88 | −88 | −73 | −7.8 | −7.5 | −6.2 | 1680 | 12.5 | 12.7 |
| shallowTop 850 hPa | −4.15 | −92 | −90 | −65 | −8.0 | −7.6 | −5.8 | 1678 | 12.6 | 12.4 |
| deck entrainmentEfficiency 0.1 | −4.47 | −112 | −91 | −69 | −8.8 | −7.6 | −6.1 | 1630 | 12.9 | 13.5 |
| deck prognosticHeight off | −4.62 | −107 | −83 | −63 | −8.6 | −7.2 | −5.6 | 752 | 13.4 | 13.2 |
| cloudTop.maximumHeight 1500 m | −4.45 | −120 | −95 | −70 | −9.1 | −7.7 | −5.9 | 1822 | 12.7 | 12.9 |
| turbulence 'dry' | −3.95 | −91 | −68 | −57 | −7.9 | −6.4 | −5.4 | 1912 | 12.6 | 12.7 |
| cloudTop.maximumHeight 1000 m | −3.98 | −93 | −77 | −62 | −8.0 | −6.9 | −5.7 | 1802 | 12.8 | 13.0 |
| plumeMomentum off | −5.01 | −130 | −108 | −80 | −9.5 | −8.4 | −6.5 | 1695 | 13.0 | 12.8 |

Single states (day 1609 of twelve-day runs) also showed the Gregory
entrainment law (`plumeEntrainmentLaw: 'gregory'`) at −124 and −130 mN/m²
over 180–140W with a 5.5 hPa gradient, worse, and `convectionType:
'cloudDepth'` unchanged.

No knob of the boundary layer, the deck or the shallow plume closes more
than 0.6 hPa of the gradient or 20 percent of the stress, and
none warms the eastern 1.5 to 2.5 km layer by more than 1 K. The three
that help, a shallow plume that entrains twice as fast (−10 percent), a
cloud-top mixing ceiling at 1 km (−11 percent) and the dry turbulence
closure (−18 percent), are mitigations: if their effects multiplied,
160W–140W would still carry 1.4 times the observed stress, and none
touches the free-tropospheric part.

## What remains

The cause is identified to two structures and not yet fixed:

1. The free-tropospheric warm column over the central Pacific (finding
   2, part 2): a weak-temperature-gradient violation of 1 to 2 K at 3 to
   12 km along the equator. Diagnose the heating that holds it (the
   ITCZ and SPCZ rain over the fixed sea surface, read back from the GPU
   in `fixedSst.mjs`; the heating budget of `tropicalHeating.mjs` on the
   5N–12N band) and the Walker cell's upper branch; the convective
   closure's timescale (`BECHTOLD`, α_x = 1.8 at N=128) and the deep
   plume's entrainment are the candidate levers, and the Gregory arm shows
   the stress follows the deep convection's vigor.
2. The inversion over the cold tongue (finding 2, part 1): 200 to 400 m
   too high and smeared over a 600 m layer. The deck's ceiling rule
   (radiation.module.js `inversionLayer`: the ceiling 1 m below the
   midpoint of the layer above the inversion interface) lets the deck
   and the boundary layer's mixing stand half a layer above the
   inversion; a ceiling at the interface itself, and a finer vertical
   grid between 1.5 and 3 km, are the candidates. The deck's height alone
   is not the lever (`noprog`): the layer stays cold with the deck at
   0.75 km, so the smearing is the moist physics' (the shallow plume,
   the resolved cloud and its longwave cooling) on a 600 m layer.
3. The warm pool's part of the coupled bias is the ocean's cold warm
   pool under the recorded winds, closed in the coupled run by the
   Bjerknes loop; it needs no atmospheric change and will follow the
   central Pacific's.

Until 1 and 2 are fixed, spin-ups that start the ocean from WOA23
(`OCEAN_START=climatology`) avoid the warm-pool Bjerknes start but will
still cool the cold tongue under the central Pacific's 1.5 to 2 times
stress.
