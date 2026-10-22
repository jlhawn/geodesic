# The SST dipole in fourteen and fifteen

## The question

Spin-ups fourteen and fifteen run at N=128 on bl60 with the hybrid ocean;
fifteen adds the ETOPO 2022 geography. Both keep the same SST bias against
WOA23 in every year. Fifteen's quarterly reports give, in K (fourteen's
are close):

- **Cold tropics:** 20° S–20° N is 1.5–2.0 K cold, and the cold tongue
  (5° S–5° N) 2.7 K. The equatorial west (140–170° E) cooled from −1.8 K
  in fifteen's year 1 to −3.5 K in year 2; fourteen held it at −2.7 to
  −3.4 K through years 3–5.
- **Warm north:** 40–60° N is 0.6–1.6 K warm, and north of 66° N is
  +1.0 K.
- **Split eastern boundaries:** California −3.0, Canary +2.2 and
  Peru −1.3 K.

Gemini 4 Argon read fifteen's day-365 bias panel and offered four groups
of causes. This doc asks:

- which of them apply to this model;
- what ranks first;
- what follows for the level set and the stratocumulus deck.

The panel averages days 336–365. Day 0 is the March equinox, so that is
18 February–20 March, not December as the reading assumed.

## Instruments

- **Saved states and logs:**
  - fifteen (`runs/verda-fifteen/fifteen128_day*.bin`), fourteen
    (`runs/verda-fourteen/fourteen128_day*.bin`) and fourteen64
    (`runs/fourteen64/`, days 1643, 1734 and 1825), decoded with
    `js/stateFile.module.js`;
  - both spin-up logs (`fifteen128.log`, `fourteen128.log`), for the
    ocean's interior mean and the Antarctic extent, and the quarterly
    reports under `runs/fifteen/reports/` and `runs/fourteen/reports/`.
- **Means:** zonal and box means over four states a year (June,
  September, December, March). Mixed-layer depth uses a 0.03 kg/m³
  density step.
- **`scripts/deckTops.mjs`:** where the deck's top sits against the
  grid.
  - Deck columns are those with mlmGate > 0.5.
  - It reads each stratocumulus box, plus the 4 K and the largest θ_v
    jumps at the VOCALS 20° S sites and DYCOMS-II.
- **`scripts/levelSets.mjs`:** named level sets side by side in the US
  Standard Atmosphere 1976.
- **Earlier work:** the fixed-SST arms and budgets of
  `docs/equatorial-wind-bias.md` and `docs/toa-radiation-bias.md`.
- **The throughflow:** fifteen's own strait transports and island-rule
  diagnostics.

## Findings

### 1. Gemini's points

| # | Claim | What the model and data show | Verdict |
|---|---|---|---|
| 1 | The Gulf Stream and Kuroshio overshoot because harmonic viscosity widens the boundary layer; use ∇⁴ and lateral slip | The closure is already ∇⁴: ν₄ ≈ 2.5 × 10¹³ m⁴/s at N=128, a Munk width (ν₄/β)^(1/5) of about 67 km, one cell. The +4–5 K spots are real, and the North Atlantic warms 0.9–1.6 K at 400–1000 m in five years. There is no lateral-slip option. | Mechanism does not apply. At 60 km, separation is limited by resolution (Chassignet and Xu 2017). The 40–60° N bias peaks in June–September (+1.1 to +2.6 K), not in winter. |
| 2a | Tropical instability waves are damped | Resolved at N=128 since the third-order tracer flux (10–12°, 40 cm/s) | Does not apply |
| 2b | GM is untapered at the equator | κ = 1000/(1 + (L_d/Δx)²) is 76 m²/s on the equator and 185 at 5°, with mode weights and a 200 m taper (`js/ocean/layered.module.js`) | Does not apply |
| 2c | No shortwave penetration, or κ_v too large | Two-band Jerlov IB light to 300 m (`js/ocean/column.module.js`); κ_b = 10⁻⁵ m²/s | Does not apply. Second order: one water type everywhere. Type II in the upwelling zones would keep about 9 W/m² more in a 25 m layer. |
| 2d | The trades are too strong, amplified by the Bjerknes feedback | Over observed SST the trades are 1.3–2.2 times observed (`docs/equatorial-wind-bias.md`, findings 1–2). At 0.6 times the stress, the ocean-only hybrid matches WOA. | Already diagnosed, and the leading cause (finding 3) |
| 3a | December mixed layers are over-mixed at 15–45° S | December–March mixed layers are 24–29 m at 20–30° S and 25–54 m at 30–40° S, and the band is cold in every season | Refuted |
| 3b | The warm ring around Antarctica comes from missing supercooled cloud or early melt | SWCRE is −49 W/m² at 60–65° S, the ring persists through winter, and September mixed layers are 370–980 m deep | Does not apply: open-ocean convection (finding 5) |
| 3c | The −1.56 K global mean implies a 15–25 W/m² surface heat loss | Net TOA is +7 W/m², and the ocean gains 9 W/m² of ocean area (finding 2) | Refuted |
| 3d | The subtropical decks are too bright, or the boundary layer too dry | SWCRE at 20–40° S is about −48 W/m² against about −40 (CERES EBAF), down from thirteen's 27 W/m² excess (`docs/toa-radiation-bias.md`, finding 3) | Already diagnosed. Reduced, not closed; the latent-heat part is untested. |
| 4 | The eastern-boundary coasts are warm because coastal upwelling is unresolved | Only the Canary coast is warm (+1.5 to +3 K in all five years). California (−1.8 to −3.1 K in December–March) and Peru are cold at the same spacing. | Does not fit. The Canary warmth needs a local cause. |

### 2. The ocean gains heat and stores it at 100–700 m

- **Interior mean (fourteen's log):** 3.64 °C at day 1 and 3.63 at
  days 61–183. It then reads 3.65, 3.67, 3.69, 3.71 and 3.73 °C at the
  ends of years 1–5.
  - That is +0.02 K a year: about +6 W/m² of Earth, or +9 W/m² of
    ocean.
  - It matches the net TOA of +7 W/m².
- **Day 1 of both runs:** Ts is 12.1 °C, OLR 190 W/m² and TOA
  +40 W/m², from the cold atmosphere of a start from nothing.
  60° S–60° N SST is already 1.1 K cold by day 91.

Fourteen, day 1825 minus day 13 (the same season): area-mean ΔT (K) by
depth (m), and the column's heat gain in W/m² of that region's ocean.

| Region | 0–20 | 20–50 | 50–100 | 100–200 | 200–400 | 400–700 | 700–1000 | 1000–2000 | W/m² |
|---|---|---|---|---|---|---|---|---|---|
| 90–60° S | +3.63 | +3.63 | +3.39 | +0.95 | −0.44 | −0.36 | −0.16 | −0.05 | 5.0 |
| 60–40° S | +2.10 | +2.06 | +1.94 | +1.17 | +1.04 | +0.66 | +0.15 | −0.05 | 17.3 |
| 40–20° S | +1.98 | +1.82 | −0.05 | −0.15 | +1.09 | +0.22 | 0.00 | −0.02 | 9.5 |
| 20° S–20° N | +0.09 | −0.35 | −0.43 | +0.97 | +0.92 | −0.02 | −0.03 | 0.00 | 6.2 |
| 20–40° N | −1.60 | −1.66 | −1.39 | +0.66 | +1.17 | −0.09 | −0.15 | +0.09 | 6.9 |
| 40–60° N | −0.42 | −0.39 | +0.06 | +0.91 | +1.49 | +0.85 | +0.41 | +0.16 | 20.9 |
| Equatorial Pacific, 5° S–5° N, 160° E–90° W | −1.06 | −1.59 | −1.56 | +1.49 | +1.46 | 0.00 | −0.02 | −0.05 | 5.3 |
| North Atlantic, 40–60° N, 70–10° W | −0.67 | −0.63 | −0.41 | −0.09 | +0.74 | +1.64 | +0.92 | +0.31 | 30.9 |
| Global | +0.66 | +0.45 | +0.14 | +0.70 | +0.92 | +0.16 | +0.02 | +0.01 | 9.3 |

- **The surface bins** carry the start's offset: a WOA annual mean
  against a March start, so the south starts warm and the north cold.
  The signal is in the subsurface bins and the five-year rates.
- **The equatorial Pacific** cools 1.1–1.6 K at 0–100 m and warms 1.5 K
  at 100–400 m. That is Ekman and subtropical-cell pumping under strong
  trades.

### 3. The causes, ranked

1. **The atmosphere's excess trades, and the ocean's response to them
   (documented).** They explain:
   - the cold 20° S–20° N band and the cold tongue;
   - the equatorial west's cooling, through the Bjerknes feedback;
   - the thermocline's heat uptake (finding 2).

   Finding 4 attributes the excess.
2. **Southern Ocean open-ocean convection (finding 5).** It explains
   the warm ring and the missing Antarctic ice. By ice-albedo
   arithmetic it is also about 1 W/m² of the TOA excess.
3. **Excess subtropical cloud (documented).**
   - The 8–12 W/m² at 20–40° S is worth about 0.5 K of the cold southern
     subtropics, at a surface response of about 20 W/m² per K.
   - The colder SST then raises the inversion strength and the cloud.
4. **Summer warmth at 40–60° N (unresolved).**
   - It is largest in the marginal seas and the NE Pacific, and the
     Arctic is ice-free every September.
   - The zonal SWCRE is −45 W/m² against about −55, though the June
     ocean boxes below reach −109 and −111 W/m².
   - The boundary-current overshoot adds to it, locally and at depth.
5. **The Indonesian throughflow excess (secondary, on the ocean side;
   finding 6).** It is a candidate for part of the cold warm pool only.

Fourteen's year-5 boxes, as SWCRE day mean (W/m²) / wind speed (m/s):

| Box | March | June | September | December |
|---|---|---|---|---|
| Canary, 15–30° N, 25–10° W | −29 / 7.8 | −40 / 7.9 | −40 / 5.3 | −62 / 9.5 |
| California, 30–42° N, 130–115° W | −72 / 6.4 | −33 / 7.4 | −27 / 6.5 | −17 / 9.4 |
| Peru, 20–5° S, 85–70° W | −10 / 7.0 | −94 / 9.1 | −73 / 5.6 | −86 / 6.5 |
| NE Pacific, 40–60° N, 180–130° W | −65 | −111 | −44 | −19 |
| North Atlantic, 40–60° N | −71 | −109 | −53 | −18 |

The ocean is far from equilibrium: the 60° S–60° N cold bias shrinks by
about 0.1 K a year. Arms are judged on the pattern, under fixed-SST or
ocean-only protocols, not on the coupled global mean.

### 4. Why the trades are too strong

**The excess pressure gradient.**
- Over observed SST, sea-level pressure rises about 4.5 hPa from
  180–160° W to 120–100° W, against about 2 hPa observed. 87 % of it
  sits below 2 km (`docs/equatorial-wind-bias.md`, finding 2), so the
  boundary layer's temperature sets most of it: the mechanism Lindzen and
  Nigam (1987) proposed for the tropics' surface winds.
- With fifteen's cloud options (`docs/toa-radiation-bias.md`,
  finding 8), the rise is 4.0–4.1 hPa.
- The central-three stress is −85 to −91 mN/m² against about −52
  observed: 1.7 times, and 2.2 times at 160–140° W.

**Its split,** measured on bl36 and not re-measured on fifteen's options:
- about 1 hPa comes from the east's 1–2 km being 2–3 K cold;
- about 1.2 hPa comes from the central Pacific's 3–12 km being 1–2 K
  warm relative to the east, of which the tropopause layer takes back
  0.8.

**The SST start is not the cause.**
- The fixed-SST arms hold WOA's SST, so no ocean feedback runs, and the
  atmosphere alone builds the excess within days.
- The coupled loop (Bjerknes) then amplifies it west of the dateline.

**The inversion is no longer smeared, but it is still too high.**
- On bl36, the layers between 1.2 and 4.7 km were 450–600 m thick. The
  jump became a ramp from 1.5 to 3 km, with an in-between layer at
  relative humidity 0.67.
- On bl42 and bl60 it is a 2.5 K step between 1.62 and 2.17 km, with
  humidity 0.87 below it and 0.43 above.
- It still sits at about 1.7 km, against about 1.5 km observed (Yin and
  Albrecht 2000). The trades did not move when the smearing went.

**Why the inversion sits too high.** Over the cold tongue the deck runs
in 1 % of columns (fifteen, year 3), so this is the trade-cumulus
regime.
- **Shallow plumes cool the inversion layers.** They detrain there at
  1.7–2.7 K/day, the largest term after radiation. The budgets, in
  K/day:

  | Term | At 1.6 km | At 1.9 km |
  |---|---|---|
  | Dynamics | +5.6 | +6.0 |
  | Longwave | −2.7 | −3.3 |
  | Plumes | −2.7 | −1.7 |
  | Boundary layer | −0.5 | −0.7 |

- **The plumes feed a loop.** Grant's (2001) closure sizes them from the
  surface buoyancy flux, which here is mostly evaporation under 8–9 m/s
  winds: 15–22 g/m²/s at cloud base.
  - The loop runs from strong trades, through evaporation and plumes, to
    a high, cold cloud layer.
  - That cold layer feeds back into a larger pressure gradient.
- **Longwave cooling** at the moist cloud top is 2.7–3.3 K/day.
- **Arms have ruled out:**
  - the deck's height;
  - the boundary layer's exchange;
  - entrainment efficiency: maximumEfficiency 0.5 lowered the inversion
    250 m but cooled the air above it, and the stress rose;
  - the grid alone.
- **The best arm** in the wind doc's table is cumulusOvershoot 0 with
  cloudTop.cumulusCapped.
  - It gives −108, −79 and −44 mN/m² at 160–140, 140–120 and
    120–100° W, against the base's −112, −83 and −66.
  - It saturates the 1.64 km layer.

**Why the air above it is too cold.**
- **At 1.3–2.2 km,** the cause is the inversion's height. Cloud-layer air
  sits there at 12.5 °C at 1.5 km, where 15–16 °C belongs.
- **At 3–4 km,** 120–100° W is 1.2–1.35 K colder than the central
  Pacific, which is worth about 1 hPa.
  - Longwave cooling there is 2.8–2.9 K/day, against 1.8–2.0 K/day of
    subsidence warming.
  - The layer is too moist from plume detrainment: relative humidity
    0.67 at 2.1 km, against about 0.4.
  - The column above an overcast deck cools 0.8 K/day faster than the
    central Pacific's.
  - Excess radiative cooling in the subsiding region forces faster
    descent, which means a stronger Walker cell (Fermepin and Bony
    2014).
- **Numerics (minor):** the centered interface value passes only about
  half of the subsiding air's excess heat to the layer under the jump,
  about 11 W/m².

**Why the central Pacific is warm relative to the east.**
- **Only the east–west difference is wrong.** Against ERA5 the tropical
  free troposphere is about 1 K cold (0–20° N: −0.8 K at 500 hPa,
  −1.3 K at 700 hPa). The east–west difference has not yet been
  compared with ERA5 longitude by longitude.
- **Deep convection heats too far east.** Coupled rain at 5–12° N is
  11.4 mm/day at 140–120° W (observed 6–7) and 3.2 at 120–100° W
  (observed 7.5).
  - The stress follows convective vigor: the Gregory entrainment arm
    made it worse.
  - The levers the wind doc names are the closure's timescale and the
    deep plume's entrainment.
- **Lower ranked:**
  - Longwave heating under thick anvils, which strengthens the Walker
    circulation (Rädel et al. 2016). In-cloud τ is about 24, and only
    0.39 of the high cover is thin, against about 0.6 observed.
  - The moist turbulence closure, which keeps the central mixed layer
    about 1 K warm. The dry closure is the best single arm (−22 %) and
    still leaves 1.4 times observed.

**Ruled out** in the wind doc (findings 3–7):
- the drag law;
- the ocean;
- the boundary-layer and deck depth;
- the plumes' reach;
- plume momentum (on by default, −18 %);
- the boundary layer's, the deck's and the shallow plume's knobs, none
  of which closes more than 0.6 hPa of the gradient or 20 % of the
  stress;
- the resolved inversion;
- the cloud-variance pair (−12 %);
- the buoyant plume.

### 5. Southern Ocean open-ocean convection

Sea cells at 60–75° S deeper than 500 m, in fifteen's first winter:

| Day | SST (°C) | SSS | 150–400 m T (°C) / S | Mixed layer (m) | Boundary layer (m) | Share with mixed layer > 300 m | Iced |
|---|---|---|---|---|---|---|---|
| 5 | −0.86 | 34.040 | 0.94 / 34.544 | 72 | 66 | 0.00 | 0.25 |
| 56 | −0.84 | 34.207 | 0.71 / 34.509 | 214 | 168 | 0.08 | 0.24 |
| 91 | −0.85 | 34.257 | 0.47 / 34.485 | 305 | 215 | 0.24 | 0.24 |
| 141 | −0.79 | 34.308 | 0.18 / 34.465 | 537 | 304 | 0.44 | 0.22 |
| 183 | −0.72 | 34.333 | 0.08 / 34.462 | 684 | 253 | 0.53 | 0.22 |
| 232 | −0.34 | 34.323 | 0.12 / 34.463 | 344 | 56 | 0.29 | 0.16 |

- **The surface never nears freezing** (about −1.85 °C).
  - The ePBL boundary layer is 168 m deep by mid-May.
  - The deep water gives up its heat and salt upward.
  - Once the halocline is gone, the convection sustains itself
    (Martinson 1990).
- **September mixed layers**, fifteen / fourteen:
  - 979 / 814 m at 65–75° S;
  - 389 / 370 m at 60–65° S;
  - 427 / 383 m at 50–60° S.
  - About 100 m is observed under ice (Pellichero et al. 2017).
- **Ice:**
  - The iced share at 60–65° S is 0.00 every September.
  - Fourteen's September Antarctic extent falls 5.9, 3.8, 3.0, 2.1 and
    1.9 Mkm² over years 1–5. Fifteen's is 5.5 and 4.1.
- **Fourteen at day 1643:** south of 60° S, the top 200 m is homogeneous
  at +0.15 °C and 34.36. The start had −0.95 °C and 34.13 over +0.92 °C
  and 34.57 at 200–400 m.
- **fourteen64 at day 1643** is in the same regime: SST +0.29 °C,
  mixed layer 446 m, boundary layer 261 m, iced 0.06.
- **Cloud does not explain it.** SWCRE at 60–65° S is −49 W/m², and the
  ring persists through the dark months.

### 6. The Indonesian throughflow excess

- **Size:** both runs carry 2–2.5 times the observed throughflow. That
  is 30 Sv (fourteen) and 42 Sv (fifteen) through all exits, against
  about 15 Sv (INSTANT; Gordon et al. 2010).
- **Not wind-driven:**
  - The island rule (Godfrey 1989) on each run's own stress gives 18
    and 21 Sv, against about 16 from observed winds.
  - So 12–16 Sv flows beyond what the winds require, independent of the
    trades.
- **Fifteen's extra 8–12 Sv:**
  - 3–5 Sv comes from its stronger winds.
  - The rest goes through two coastline gaps near Bali: +5.2 Sv at
    Lombok's second edge and +3.4 Sv through a Java–Bali gap.
- **Where the extra water goes:**
  - 80 % of fifteen's extra flow is above 300 m, and fifteen carries
    0.41 ± 0.23 PW more heat west across 120° E.
  - Its Leeuwin box (22–35° S, 105–115° E) is 0.7 K warmer than
    fourteen's from year 2.
  - The outflow boxes are still 1.8–2.1 K below WOA23 in both runs.
- **A candidate for part of the cold warm pool:**
  - In coupled models a stronger throughflow cools the west Pacific
    warm pool (Schneider 1998; Song, Vecchi and Rosati 2007).
  - Fifteen's extra flow did not make its warm pool (10° S–10° N,
    120–160° E) robustly colder than fourteen's. The two runs' gaps
    against WOA23 differ by −0.01 and −0.06 K in years 1–2; in year 3
    fifteen is 0.5 K colder at two of the four quarters.
  - It ranks behind the trades and their Bjerknes feedback.
  - It cannot explain the cold tongue, the southern subtropics,
    40–60° N or the Southern Ocean.
- **The figures are instants, not means.** Every state is saved at
  00 UTC, and a fast barotropic exchange locked to that time leaves
  budgets out of balance: 12.9 Sv in fourteen's Indonesian seas.

### 7. The deck's top locks to the grid

**How the ceiling works** (`js/physics/radiation.module.js`,
`inversionLayer()` and `shadowDeck()`):
- The deck's prognostic top h moves by entrainment against subsidence.
- It is capped 1 m below the midpoint of the layer above the lowest
  interface with a θ_v jump of at least 4 K.
- When the jump drops a layer, h drops the whole step at once.
- Where the deck runs, the boundary layer's K-profile spans h
  (`deckTop` in `js/physics/boundaryLayer.module.js`). K at the jump,
  κu* z (1 − z/h)², therefore grows with the square of half the
  thickness of the layer above.

**In fifteen's year 3** (days 779–1095, eight states):
- In the three layers between 1.2 and 2.0 km, 45–69 % of the deck tops
  each layer holds sit in its middle fifth, in every deck box. A top that
  moved freely would put 20 % there.
- 12–32 % of a box's deck tops lie within 5 m of the ceiling.
- In the SE Pacific (25–10° S, 90–72° W), the three 100 m bins that hold
  the ceilings (1.3, 1.5 and 1.8 km) hold 62 % of deck tops. The bins
  between them hold 4–8 % each.
- Over the cold tongue (3° S–3° N, 140–90° W) the deck runs in 1 % of
  columns.

**Against VOCALS and DYCOMS-II:**
- On the VOCALS line in September–December of years 2–3 (days 548, 599,
  639, 913, 962 and 1004), the 4 K jump sits at:
  - 1.31 km at the coast (73–70° W), against 1.0 observed;
  - 1.44 km at 77–73° W, against 1.2;
  - 1.60 km at 87–83° W, against 1.6 (Bretherton et al. 2010).
- The largest θ_v jumps there are at 1.46, 1.48 and 1.70 km.
- At the coast in year 3's eight states, the 4 K jump is the interface
  near 1.2 km in 62 % of deck columns and the one near 0.9 km in only
  10 %.
- In the DYCOMS-II box (30–33° N, 124–119° W) in June–August (days 456,
  506, 821 and 872), the largest jump is at 0.76 km and the 4 K jump at
  0.72 km, against 0.84 km observed.

**Each step changes the water twofold.** The deck's water path goes as
the square of its thickness, so one ceiling step of 250–290 m changes
it twofold or more.

**The grid sets part of the entrainment.** The index z(1 − z/h)² (per
κu*, in m), with the deck's top half a layer above the jump:

| Jump at (m) | 482 | 649 | 826 | 882 | 1014 | 1213 | 1450 | 1698 |
|---|---|---|---|---|---|---|---|---|
| bl60 | 10.5 | 15.0 | | 22.0 | | 9.6 | 9.0 | 8.7 |
| bl60 with 882 m split at 826 and 1014 m | 10.5 | 9.3 | 8.6 | | 8.1 | 9.6 | 9.0 | 8.7 |

- A deck capped at 882 m is mixed through its inversion 2.3 times as
  hard as one at 1.2 km, which pushes it up exactly where the coast is
  too high.
- Lock (2001) describes this grid-dependent entrainment.

**The deck's radiation sits in a fixed layer.**
- The deck's water sits in the one layer nearest σ = stratusSigma =
  0.92, wherever h is (`stratusLayer` in `js/physics/radiation.module.js`).
- The WGSL takes the same constant from `js/gpu/core.gpu.js` and
  `js/gpu/forcing.gpu.js`.
- On bl60 that layer is 649–882 m. With tops at 1.3–1.8 km in the
  SE Pacific, the deck's longwave emission and cooling come from 0.4–1.0
  km below its top.

### 8. The level set: what bl60's layers bought, and bl48

**bl60's 60 layers,** in standard heights:
- **Top:** two layers above 1 hPa, 0–0.3 hPa and 0.3–1 hPa (9,334 m).
- **Stratosphere:** 27 layers equal in ln σ from 1 to 93.6 hPa,
  1,066–1,330 m each.
- **Free troposphere:** layers of 1.0–1.2 km.
- **Inversion zone:** bl42's cuts, 237–338 m from 1.2 to 3.5 km, with
  two layers of about 600 m above them. The cold-tongue inversion sits in
  layers 46–48, at 1.45–2.23 km.
- **Boundary layer:** thinning to a 40 m lowest layer.
- **Layer 50 (882–1,213 m, 331 m)** lies between layers of 233 and
  237 m.
  - It is the top rung of the bl34 ladder, each layer about 1.4 times
    the one below.
  - bl42 re-cut only 1.2–4.7 km, so the rung stayed.

**What the 18 layers bl60 adds to bl42 bought.** All of them sit between
1 and 93.6 hPa, and they were added for the QBO (`3de3e6c`;
`docs/tropical-upper-stratosphere.md`, finding 5).
- **The radiation gain came from the table, not the layers.** bl60's
  longwave table was fitted on its own columns and is within
  0.1–0.3 K/day of RRTMG at 1–2.2 hPa. bl36's table on bl60's layers
  misses by 0.6–0.9 K/day.
- **Temperatures barely changed.** Over 30 days, bl42 (stratospheric
  layers of 1.1–6.2 km) and bl60 differ by at most 1.4 K at 3–50 hPa.
- **The exception is under the lid.** bl42's 6.2 km layer there ascends
  too fast: w* is 3.88 mm/s at 1 hPa, against 0.66 on bl60, and the
  layer is 4.6 K colder.
- **No QBO.** Fourteen's oscillation still lacks its easterly phase.
- **The cost:** 1.40 times bl42 per simulated day at N=128, and 1.48
  times at N=64. Radiation is about 3 % of an N=128 step (longwave 6.7
  and shortwave 2.1 ms every fourth step). The atmosphere's cost scales
  with the layer count.

**Geometric ladders.** A ladder that grows by a fixed ratio r from 40 m
gives each layer a thickness of about (r − 1)(z + 40 m/(r − 1)), where z
is the height of its base.
- **The ratio fixes the thickness at every height,** so it cannot add a
  layer to one band alone. At 2^(1/5), every doubling of z + 269 m holds
  five layers.
- **∛2 for 54 layers runs out of atmosphere.** It passes the stratopause
  at layer 25, and its 54th layer would end 40,000 km up.
- **A single ratio from 40 m to 16.6 km ends too thick.** Its last
  layers are 2,171 m (2^(1/5), 59 layers in all), 1,858 m (2^(1/6), 63),
  1,421 m (2^(1/8), 71) or 1,163 m (2^(1/10), 78). The stratosphere's
  first layer is 1,066 m.
- **So the ladders grow to a cap and hold it.** With a 1.17 km cap and
  bl60's 29 layers above 93.6 hPa:

| Rule | Layers | At 0.9 km | At 1.5 km | At 2 km | First capped layer |
|---|---|---|---|---|---|
| ∛2 | 54 | 254 m | 403 m | 508 m | 16, at 4.8 km |
| 2^(1/4) | 58 | 190 m | 320 m | 381 m | 21, at 6.6 km |
| 2^(1/5) | 61 | 160 m | 243 m | 320 m | 26, at 8.3 km |
| bl60 | 60 | 331 m | 248 m | 270 m | |

- **The four-section sets:** bl63 and bl64 put bl60's stratosphere over a
  2^(1/5) troposphere and grow ×2 per layer above 1 hPa (four layers
  there, against bl60's two).
  - bl64 matches its cap to the stratosphere's first layer.
  - Neither was built.
  - A K of 64 buys nothing on the GPU: every kernel dispatches a fixed
    workgroup, and the buffers are K × C and K × E.

**bl48, the set Joshua chose under a 48-layer limit:**
- **Troposphere:** ×2^(1/5) per layer from 40 m to 8.34 km (25 layers),
  then 7 layers of 1,180 m to bl60's 93.59 hPa interface
  (σ = 0.0923658).
- **Stratosphere:** the ratio continues from the cap through 1,356,
  1,557, 1,789, 2,055 and 2,361 m to 25.7 km, then 9 layers of 2,456 m
  run to 1 hPa (14 layers).
- **Top:** bl60's 0.3 hPa interface and the open top (2 layers).
- **Ratios:** every layer below 1 hPa is 1.00–1.15 times the one below
  it. The layer above 1 hPa is 3.8 times the last stratospheric one; on
  bl60 it is 7 times.
- **The conservative alternative R1** keeps bl60's troposphere (31
  layers). Above it come 5 layers growing ×2^(1/5) from 1,031 m, then
  10 of 2,326 m to 1 hPa. It keeps bl60's thickness jumps: ×0.72 at
  1.2 km and ×1.88 at 4.7 km.

Thickness (m) of the layer that holds each height, and layers per band
by middle height:

| Set | Layers | 0.5 km | 1 | 2 | 2.5 | 3 | 5 | 12 | 18 | 25 | 45 | 0–1 / 1–2.5 / 2.5–5 / 5–16.6 km / 16.6 km–1 hPa / above |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| bl42 | 42 | 167 | 331 | 270 | 290 | 330 | 1,205 | 1,031 | 1,538 | 3,267 | 6,175 | 9 / 6 / 5 / 11 / 9 / 2 |
| bl60 | 60 | 167 | 331 | 270 | 290 | 330 | 1,205 | 1,031 | 1,066 | 1,089 | 1,298 | 9 / 6 / 5 / 11 / 27 / 2 |
| bl48 | 48 | 106 | 184 | 320 | 368 | 485 | 735 | 1,180 | 1,557 | 2,360 | 2,456 | 11 / 6 / 4 / 11 / 14 / 2 |
| R1 | 48 | 167 | 331 | 270 | 290 | 330 | 1,205 | 1,031 | 1,360 | 2,326 | 2,326 | 9 / 6 / 5 / 11 / 15 / 2 |

- **Finer and coarser than bl60:** bl48 is finer below 1.2 km and at
  4.7–7 km, about equal at 1.2–2 km, and 1.2–1.5 times coarser at
  2–3.5 km. The deck's mixing index for jumps between 1.9 and 3.6 km is
  13–19 on bl48. On bl60 it is 8–9, except 22 at its 3.5 km interface.
- **Estimated pace,** linear in the layer count from the timings above:
  31.8 s per simulated day at N=128 (bl42 28.0, bl60 39.3), and 4.6 at
  N=64 (3.95 and 5.85).
- **Precedents:**
  - MPI-ESM-LR runs 47 levels to 0.01 hPa without a QBO (Stevens et al.
    2013).
  - GFDL AM3's layers are 1–1.5 km near the tropopause and 3–4 km in
    much of the stratosphere (Donner et al. 2011).
  - GISS-E2.1 has 40 levels to 0.1 hPa (Kelley et al. 2020).
- **The σ interfaces,** top to ground:

```
[0, 0.0002960770, 0.0009869233, 0.001347988, 0.001855411, 0.002575878,
 0.003608678, 0.005104253, 0.007293249, 0.01052185, 0.01525079,
 0.02219543, 0.03195852, 0.04403563, 0.05834921, 0.07458910, 0.09236580,
 0.1112576, 0.1340133, 0.1614233, 0.1944395, 0.2341587, 0.2806937,
 0.3344493, 0.3925796, 0.4496087, 0.5045532, 0.5566991, 0.6055710,
 0.6508927, 0.6925483, 0.7305454, 0.7649822, 0.7960205, 0.8238641,
 0.8487406, 0.8708888, 0.8905489, 0.9079552, 0.9233316, 0.9368886,
 0.9488215, 0.9593096, 0.9685162, 0.9765892, 0.9836615, 0.9898519,
 0.9952667, 1]
```

**What bl48 needs to be built:**
- a rule in `sigmaInterfaces`;
- bl60's sponge and lid friction;
- its own longwave table, fitted by `scripts/longwaveFit.mjs` and
  checked against RRTMG as bl60's was;
- a CPU/GPU parity test at N=16;
- a conservative remap of bl60 states.

### 9. The deck fix

Two options, both off by default, on both engines, each with a CPU/GPU
parity test at N=16:

1. **Sub-grid inversion.**
   - Reconstruct the jump's height z_i inside the layer above the capping
     interface. Treat that layer's mean θ_l and q_t as a mass-weighted
     mix of the mixed-layer values below and the free-tropospheric values
     above (Grenier and Bretherton 2001).
   - Use z_i for the deck's ceiling and as the top of the K-profile, so
     that only the explicit entrainment crosses the jump (Lock 2001).
2. **Radiative top.** Put the deck's water, column by column, in the
   layer that holds its top (z_i with option 1, h without it), instead
   of the fixed σ = 0.92 layer. A real deck's longwave cooling sits in
   its top few tens of meters (Wood 2012).

The options are level-agnostic. With them the deck's height and water
no longer depend on where the interfaces fall, and the level set can be
chosen for everything else.

## What this means

The dipole's main parts have separate causes:
- the atmosphere's trades, through the ocean's response;
- convection in the Southern Ocean;
- subtropical cloud;
- an unresolved northern summer.

The throughflow is a smaller, ocean-side part of the cold warm pool.

On the deck: its height, its water path and the height its cloud
radiates from are all set by the grid's interfaces. The sub-grid
inversion and the radiative top make them independent of the level set.

That leaves the level set free to follow the cost. bl60's extra
stratospheric layers bought no QBO. The radiation accuracy that came
with them came from the per-grid longwave table, which bl48 gets too.

Proposed and not run; every GPU run needs Joshua's go:
- **E1, the Southern Ocean.** Three N=64 starts from nothing on
  fifteen's settings, run to day 183, about 13 minutes of GPU each.
  - Arms: control, `langmuir: 0`, and a reduced `nstar`.
  - Langmuir enhancement needs surface waves, which ice damps (Squire
    2007).
  - Confirmation: an arm whose surface reaches freezing by day 91, keeps
    the boundary layer under 120 m and grows more than 12 Mkm² of ice
    names the source.
  - If no arm holds the halocline, add the missing Antarctic meltwater,
    about 0.08 Sv (Rignot et al. 2013; Bronselaer et al. 2018), as a
    change on both engines.
- **E2, wind attribution.** A coupled N=64 year from fourteen64's day
  1825, the control against ocean stress ×0.6 over 20° S–20° N. It
  needs a diagnostic stress-scale option on both engines, and about 52
  minutes of GPU for the pair. It shows how much of the southern
  subtropics, the northern warmth and the heat uptake share the wind's
  cause.
- **E3, the northern summer budget.**
  - First, seasonal SWCRE from fifteen's states against CERES EBAF. The
    EBAF file is not in `data/`.
  - Then a 30-day June surface-flux record at N=64, about 3 minutes of
    GPU.
- **Trade-wind arms.** 20-day fixed-SST runs at N=128, three runs per
  arm because one run's noise is ±4.4 mN/m². That is about 20 minutes
  per arm, plus a control set at fifteen's settings.
  - Halve the shallow plumes' cloud-base mass flux east of 140° W: the
    1.5–2.5 km east column should warm and dry, and the 1–2 km share of
    the gradient should fall.
  - Thin the high cloud's water.
  - First, on the CPU: `scripts/tropicalHeating.mjs` on the 5–12° N band
    of fixed-SST states.
- **The throughflow (CPU).** Run the island-rule and strait-transport
  diagnostics on the end states of the Oct 18 ocean-only ladder's 1.0×
  and 0.6× arms.
  - If the 0.6× arm, whose tropics match WOA, still carries 12–16 Sv
    beyond its own winds, the excess is not what holds the warm pool
    cold.
  - If the excess shrinks with the stress, a 1.0× arm with the Bali
    gaps closed measures the warm pool's change per Sv.
- **The Canary coast (CPU, about 10 minutes):** its alongshore stress
  against California's, from existing states.
- **bl48:**
  - finding 5's 30-day arm at N=64, bl60 against bl48 (w* and T at
    1–50 hPa, the cold point);
  - then a 20-day fixed-SST pair at N=128 (deck heights, the
    central-three stress, the 2–3 km inversion at 140–160° W, and
    pace).
- **The deck fix:** a 20-day fixed-SST set at N=128 with both options
  off, option 1 alone, and options 1 and 2 together. It reads the 20° S
  tops, Namibia, Canary, DYCOMS-II, the deck's shortwave and longwave
  cloud effects, and the central-three stress as a no-harm check.
- **Two rows of the quarterly report look stale for the hybrid ocean:**
  - "Mixed layer, ocean mean 3 m", where the ocean panels show about
    70 m;
  - a flat "Thermocline 233/237/234 m", where the section's 20 °C
    isotherm runs from about 150 m in the west to about 40 m in the
    east.

## References

- Bretherton, C. S., R. Wood, R. C. George, D. Leon, G. Allen, and
  X. Zheng, 2010: Southeast Pacific stratocumulus clouds, precipitation
  and boundary layer structure sampled along 20° S during VOCALS-REx.
  Atmos. Chem. Phys., 10, 10639–10654.
- Bronselaer, B., M. Winton, S. M. Griffies, W. J. Hurlin, K. B. Rodgers,
  O. V. Sergienko, R. J. Stouffer, and J. L. Russell, 2018: Change in
  future climate due to Antarctic meltwater. Nature, 564, 53–58.
- Chassignet, E. P., and X. Xu, 2017: Impact of horizontal resolution
  (1/12° to 1/50°) on Gulf Stream separation, penetration, and
  variability. J. Phys. Oceanogr., 47, 1999–2021.
- Donner, L. J., et al., 2011: The dynamical core, physical
  parameterizations, and basic simulation characteristics of the
  atmospheric component AM3 of the GFDL global coupled model CM3.
  J. Climate, 24, 3484–3519.
- Fermepin, S., and S. Bony, 2014: Influence of low-cloud radiative
  effects on tropical circulation and precipitation. J. Adv. Model.
  Earth Syst., 6, 513–526.
- Godfrey, J. S., 1989: A Sverdrup model of the depth-integrated flow for
  the world ocean allowing for island circulations. Geophys. Astrophys.
  Fluid Dyn., 45, 89–112.
- Gordon, A. L., et al., 2010: The Indonesian throughflow during
  2004–2006 as observed by the INSTANT program. Dyn. Atmos. Oceans, 50,
  115–128.
- Grant, A. L. M., 2001: Cloud-base fluxes in the cumulus-capped boundary
  layer. Quart. J. Roy. Meteor. Soc., 127, 407–421.
- Grenier, H., and C. S. Bretherton, 2001: A moist PBL parameterization
  for large-scale models and its application to subtropical
  cloud-topped marine boundary layers. Mon. Wea. Rev., 129, 357–377.
- Kelley, M., et al., 2020: GISS-E2.1: Configurations and climatology.
  J. Adv. Model. Earth Syst., 12, e2019MS002025.
- Lindzen, R. S., and S. Nigam, 1987: On the role of sea surface
  temperature gradients in forcing low-level winds and convergence in
  the tropics. J. Atmos. Sci., 44, 2418–2436.
- Lock, A. P., 2001: The numerical representation of entrainment in
  parameterizations of boundary layer turbulent mixing. Mon. Wea. Rev.,
  129, 1148–1163.
- Martinson, D. G., 1990: Evolution of the Southern Ocean winter mixed
  layer and sea ice: open ocean deepwater formation and ventilation.
  J. Geophys. Res., 95, 11641–11654.
- Pellichero, V., J.-B. Sallée, S. Schmidtko, F. Roquet, and
  J.-B. Charrassin, 2017: The ocean mixed layer under Southern Ocean
  sea-ice: seasonal cycle and forcing. J. Geophys. Res. Oceans, 122,
  1608–1633.
- Rädel, G., T. Mauritsen, B. Stevens, D. Dommenget, D. Matei,
  K. Bellomo, and A. Clement, 2016: Amplification of El Niño by cloud
  longwave coupling to atmospheric circulation. Nature Geosci., 9,
  106–110.
- Rignot, E., S. Jacobs, J. Mouginot, and B. Scheuchl, 2013: Ice-shelf
  melting around Antarctica. Science, 341, 266–270.
- Schneider, N., 1998: The Indonesian throughflow and the global climate
  system. J. Climate, 11, 676–689.
- Song, Q., G. A. Vecchi, and A. J. Rosati, 2007: The role of the
  Indonesian throughflow in the Indo-Pacific climate variability in the
  GFDL coupled climate model. J. Climate, 20, 2434–2451.
- Squire, V. A., 2007: Of ocean waves and sea-ice revisited. Cold Reg.
  Sci. Technol., 49, 110–133.
- Stevens, B., et al., 2013: Atmospheric component of the MPI-M Earth
  System Model: ECHAM6. J. Adv. Model. Earth Syst., 5, 146–172.
- Wood, R., 2012: Stratocumulus clouds. Mon. Wea. Rev., 140, 2373–2423.
- Yin, B., and B. A. Albrecht, 2000: Spatial variability of atmospheric
  boundary layer structure over the eastern equatorial Pacific.
  J. Climate, 13, 1574–1592.
