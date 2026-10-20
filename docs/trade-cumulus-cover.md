# The trade cumulus's cover

Why the shallow plumes' detrained water brightens the subtropics under
`cumulusOvershoot` 'kinetic', measured on 2026-10-19 on N=64 fixed-SST
arms from `runs/thirteen/thirteen64_day1825.bin` over WOA23's sea surface
(days 6–20) and the N=128 stress protocol of `equatorial-wind-bias.md`.

## The question

The kinetic share (7cb2e0c) detrains most of a shallow plume's condensate
in the layer below the inversion instead of in the dry inversion layer. In
the arms it makes SWCRE −60.5 W/m² against −54.7 under the old share 1, the
difference subtropical (9.5 W/m² at 20–40S, 8.1 at 20–40N;
`toa-radiation-bias.md`, finding 4). The share follows the plume's physics;
the question is the cover the detrained water then takes.

## Instrument

`scripts/tradeCover.mjs STATE` (CPU, one step, ocean off; `RADIATION`,
`MOIST`, `BOUNDARY_LAYER`, `SURFACE` as elsewhere): for the trade boxes of
`cloudRegimes.mjs`, the SE Pacific and the open sea of 20–40S and 20–40N,
the radiation's low cover and water, the stratiform ramp's share and the
deck, and per layer below 3.5 km the radiated cover, the variance cover's
σ_s beside the moist physics' uniform half-width b, the deficit s and Q₁,
the covers of the same layer under the alternatives below, the plume tops,
and the cloudy mixed layers split by the part of the boundary layer that
mixes them.

## Findings

### 1. The kinetic share adds cloudy layers at 0.8–1.9 km

Means of the arms' states on days 1835, 1840 and 1845:

| | low cover 20–40S | low cover 20–40N | low water 20–40S / 20–40N, g/m² |
|---|---|---|---|
| kinetic (the default) | 0.549 | 0.527 | 69.0 / 57.5 |
| share 1 | 0.456 | 0.427 | 57.9 / 46.2 |

Between 0.8 and 1.9 km the share of columns whose layer holds cloud rises
from 0.27–0.39 to 0.32–0.49 (20–40S) and the layer's water by 20–70
percent; each cloudy layer is covered at 0.45–0.52 under either share.

### 2. Every cloudy mixed layer is covered by about one half, whatever its water

Below the mixing top the radiation covers a cloudy layer by the variance
cover Φ(s/σ_s). In the trade cumulus layer its σ_s (medians 1.0–2.5 g/kg)
is four to eight times the half-width b (0.2–0.5 g/kg) of the uniform
distribution with which the moist physics condensed the water. A layer
holds water only where s > −b, so wherever it holds any, |s| ≪ σ_s,
Q₁ = −0.06 to −0.15 and the cover is 0.44–0.50: a layer of 0.001 g/kg and
one of 0.05 g/kg both cover half the cell. This holds in every part of the
boundary layer (kinetic arm, day 1845, 20–40S at 0.8–1.4 km: the
surface-driven layer, the cumulus layer that the decoupled column leaves to
the plume, the cloud-top layer and coupled columns all at 0.43–0.50). The
layers above the mixing top cover 0.35–0.48 by sqrt(q_c/b). The stratiform
ramp is not involved: its share is 0.00 in the three trade boxes (EIS below
8 K) and 0.02–0.12 over the bands.

The model's trade cloud is a veil: at 1–1.6 km a cloudy layer holds about
0.06 g/kg over 0.45–0.5 of the cell, in-cloud about 0.12 g/kg (τ about 5
per 300 m layer), where BOMEX's clouds hold 0.3–1 g/kg.

### 3. Against BOMEX

Siebesma et al. (2003), the LES intercomparison of BOMEX's trade cumulus:

- the layer cover peaks at about 0.06 at cloud base and falls to about
  0.015 at 1.5 km (their Fig. 6); the total cover is 0.13 (0.08–0.17
  across the models, none beyond about 0.25 under doubled fluxes, no
  subsidence or a moister subcloud layer);
- the domain-mean liquid water peaks at about 0.006 g/kg (Fig. 3d), the
  in-cloud mean is 0.3–0.5 g/kg at 1 km and about 1 g/kg at 1.5 km
  (Fig. 7d);
- σ_s "has a typical value in the cloud layer of 0.2 g kg⁻¹", and the
  Gaussian cover with it is close to the LES's, a little low below 10
  percent from the cumulus field's skewness (§3b.3);
- the cloud layer's TKE is 0.10–0.13 m²/s² (Fig. 5a).

The variance cover's own σ_s on BOMEX's mean profile at 1 km (∂q_t/∂z
−5.8 g/kg/km, ∂θ_l/∂z 3.85 K/km, a = 0.29, l = κz/(1 + κz/300 m) = 171 m,
scale 5) is 2.4 g/kg, twelve times the LES's, and with BOMEX's mean deficit
it covers 0.43–0.46 of the layer against the LES's 0.03. The formula is the
local-gradient variance of Mellor and Yamada's level 2 with a neutral
length; the cumulus layer's environment is stably stratified (θ_v rises
2.5–3.5 K/km, N ≈ 0.01 s⁻¹), where the turbulence length is bounded by
the stratification: l ≤ 0.53 q/N (Galperin et al. 1988), or 0.76 e^½/N
(Deardorff 1980), about 25–30 m at BOMEX's cloud-layer TKE. With that length the
formula gives 0.46 g/kg at the scale 5 (0.30 at Mellor and Yamada's
B2^½ = 3.2), within a factor of two of the LES.

A shorter length alone does not darken the model's trade layers: on the
kinetic arm's states their medians fall to 0.3–0.5 g/kg and their covers
only to 0.42–0.48, because the layers that hold water sit at s ≈ 0, where
any Gaussian covers about half. What the model lacks against BOMEX is the
relation between the cover and the water: BOMEX holds a fifth of the
model's grid-mean water per layer (0.005–0.006 against 0.02–0.03 g/kg) in
a fourth to a tenth of its cover (0.015–0.06 against 0.19–0.24).

### 4. The cover of the water the layer holds

The layers above the mixing top already cover the water their distribution
holds (sqrt(q_c/b) is the uniform distribution's cover of its condensate
q_c). The variance cover can do the same: place its Gaussian, of the same
σ_s, where it holds the layer's water, σ_s G(Q) = q_c with G(Q) = QΦ(Q) +
φ(Q) the Gaussian's condensate (Sommeria and Deardorff 1977; Bougeault
1981), and cover Φ(Q). On BOMEX it reproduces the LES: σ_s 0.2 g/kg and
q_l 0.005 g/kg give 0.058, the length-bounded 0.46 g/kg gives 0.028,
against the LES's 0.03–0.06. In a well-mixed stratocumulus layer, where
σ_s sits at its floor, it is overcast as before.

## The change

Two options of the radiation, both engines, off by default (the defaults
are bitwise unchanged):

- `varianceCover` 'water': the cover of section 4, by eight Newton steps
  from Q = q_c/σ_s (G is convex, so they descend onto the root; the cover
  floor takes the far tail).
- `buoyancyLength` c: in a layer whose θ_v rises to its neighbors under a
  buoyant surface, the variance cover's length is at most c w*/N, w* =
  (B₀ h)^⅓ the boundary layer's velocity standing for q, N the layer's
  buoyancy frequency; c = 0.53 (Galperin et al. 1988). BOMEX: w* 0.62
  m/s, l ≈ 34 m.

Used together; 'water' alone covers with the unbounded σ_s, twelve times
BOMEX's.

On the kinetic arm's states held fixed (one step), low cover, means of days
1835–1845:

| | N Pacific trades | S Pacific trades | Atlantic trades | SE Pacific | sea 20–40S | sea 20–40N |
|---|---|---|---|---|---|---|
| defaults | 0.515 | 0.522 | 0.513 | 0.583 | 0.549 | 0.527 |
| 'water' | 0.248 | 0.197 | 0.280 | 0.350 | 0.357 | 0.379 |
| `buoyancyLength` 0.53 | 0.509 | 0.519 | 0.496 | 0.599 | 0.554 | 0.541 |
| both | 0.364 | 0.302 | 0.388 | 0.473 | 0.459 | 0.477 |

## Arms

N=64, thirteen64_day1825 over WOA23, means of days 6–20 (`fixedSst.mjs`
with a diagnostics frame each day; "both" is `varianceCover` 'water' with
`buoyancyLength` 0.53):

| arm | ASR | OLR | SWCRE | LWCRE | clear-sky OLR | net | rain mm/d |
|---|---|---|---|---|---|---|---|
| defaults (kinetic) | 229.7 | 226.1 | −60.5 | +33.7 | 259.8 | +3.6 | 3.12 |
| cumulusOvershoot 1 | 235.5 | 226.7 | −54.7 | +33.0 | 259.7 | +8.8 | 3.08 |
| **both** | 234.2 | 227.4 | **−55.9** | +32.5 | 260.0 | +6.8 | 3.02 |
| both, cumulusOvershoot 1 | 239.6 | 228.0 | −50.6 | +32.1 | 260.0 | +11.6 | 3.02 |
| 'water' alone | 240.4 | 228.9 | −49.7 | +31.4 | 260.3 | +11.4 | 2.95 |
| cloudLifetime 0.5 h | 233.9 | 227.5 | −56.3 | +32.3 | 259.8 | +6.5 | 3.09 |

Low cover in the arms, means of days 1835–1845:

| arm | N Pacific trades | S Pacific trades | Atlantic trades | SE Pacific | sea 20–40S | sea 20–40N |
|---|---|---|---|---|---|---|
| defaults | 0.515 | 0.522 | 0.513 | 0.583 | 0.549 | 0.527 |
| cumulusOvershoot 1 | 0.430 | 0.554 | 0.395 | 0.438 | 0.456 | 0.427 |
| both | 0.424 | 0.317 | 0.441 | 0.472 | 0.448 | 0.463 |
| both, cumulusOvershoot 1 | 0.269 | 0.349 | 0.223 | 0.302 | 0.340 | 0.354 |
| Earth | 0.2–0.4 | 0.2–0.4 | 0.2–0.4 | 0.6–0.8 | | |

(Earth as `cloudRegimes.mjs` gives it: MODIS and CALIPSO-GOCCP for the
trades, Klein and Hartmann 1993 for the SE Pacific.)

The pair removes 4.6 of the overshoot's 5.8 W/m² and brings the subtropical
bands' low cover to share 1's, the trades into Earth's range or near it,
with the kinetic share kept. Under it the kinetic share still brightens by
5.3 W/m² against share 1: it adds about a fifth to the low water, and the
cover now grows with the water. The cost: the SE Pacific's resolved low
cover falls from 0.58 to 0.47 against Earth's 0.6–0.8 (the deck's gate is
open on 0.01–0.10 of it in these states; the resolved cloud there is the
decoupled cumulus under a high inversion that M22 found), LWCRE by 1.2
W/m², and over WOA23's sea the planet gains 3.2 W/m² more, the clear-sky
column's deficit still uncorrected (`toa-radiation-bias.md`, finding 7).

The lifetime arm brightens less by turning cloud water into rain
everywhere; the detrained remnants of real trade cumulus evaporate within
minutes (Sakradzija, Seifert and Heus 2015), which in the model is the
saturation adjustment's part, so it is not the lever.

### The trade stress at N=128

The protocol of `equatorial-wind-bias.md`: N=128, `twelve_day1597` over
WOA23 on bl42, means of the states of days 1609–1617 by `armMetrics.mjs`,
sea cells within 3° of the equator.

| τ_x, mN/m² | 140E–160E | 160E–180 | 180–160W | 160W–140W | 140W–120W | 120W–100W | 100W–80W |
|---|---|---|---|---|---|---|---|
| observed, about | −15 | −20 | −40 | −50 | −55 | −50 | −30 |
| defaults | −7.3 | −45.0 | −116.2 | −142.7 | −100.2 | −63.1 | −15.9 |
| both | −9.0 | −30.4 | −85.5 | −99.7 | −81.4 | −57.7 | −18.0 |

Every bin moves toward the observed stress: 160W–140W from 2.9 to 2.0
times it, the mean of the three central bins from −102.0 to −79.6. The
sea-level pressure of 180–160W less 120W–100W's goes from −5.34 to −3.85
hPa. Of the 1.49 hPa, the layers below 2 km give 0.41 and those of 2–5 km
0.99: the 120W–100W column is 0.6–0.9 K warmer from 1.9 to 4 km and drier
(relative humidity at 2.8–4 km 0.20–0.39 against 0.45–0.49), while the
180–160W column there stays within 0.4 K. Over 180–160W the cumulus layer's
top (the highest layer covered above 0.01) falls from 1.91 to 1.69 km under
a sharper inversion: at 1.92 km 16.7 °C and relative humidity 0.49 against
15.1 °C and 0.72. How the pair warms the east column aloft is not traced
here.

Day means, days 6–20:

| arm | ASR | OLR | SWCRE | LWCRE | clear-sky OLR | net | rain mm/d |
|---|---|---|---|---|---|---|---|
| defaults | 235.8 | 239.1 | −60.3 | +29.0 | 268.1 | −3.3 | 3.08 |
| both | 239.6 | 239.6 | −56.4 | +28.5 | 268.1 | 0.0 | 3.01 |

## Tests

- `test/gpuModel.test.mjs`: a parity test of 'water', `buoyancyLength`
  0.53 and both at N=8 beside the variance cover's: the options move the
  layer heating by 3.5–6.5 K/day and the engines agree to 2.5–7.8·10⁻⁴
  K/day against a largest 26.4.
- CPU radiation and column files (cloudEffects, cloudHeights, column,
  infraredSlices, insolation, latitudeColumns, longwaveOverlap,
  radiationEvery, radiativeColumn, wetColumn, plumeColumn): 93 of 93.
- `test/physics.test.mjs`, whose pinned CPU digests hold with both options
  off: 37 of 37. `test/gpuModel.test.mjs` 32 of 32 and
  `test/boundaryLayer.test.mjs` 16 of 16, with no command-buffer aborts
  during the arms or the tests.

## What remains

- The scale 5 against Mellor and Yamada's 3.2 (M22's choice): with the
  bounded length BOMEX's σ_s is 0.46 g/kg at 5 and 0.30 at 3.2.
- The SE Pacific deck's cover, which the resolved cloud does not stand for.
- The regime covers at N=128 (measured here at N=64 only), and the path
  by which the pair warms the 120W–100W column at 2–4 km.
- A coupled arm: the stress and the cover both act on the ocean.
- With `cloudColumns` and `cloudSpread` (`toa-radiation-bias.md`) the pair
  is the shortwave's remaining cloud-amount change; the net balance needs
  the clear-sky column first.
