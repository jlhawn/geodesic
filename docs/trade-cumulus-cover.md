# The trade cumulus's cover

Why the shallow plumes' detrained water brightens the subtropics under
`cumulusOvershoot` 'kinetic', measured on 2026-10-19 and 20 on N=64 fixed-SST
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

## The plume's reach

Under 'kinetic' the share is the reach of w², from `plumeVelocity` 1 m/s
at the first cloudy layer by dw²/dz = 2aB − 2bεw² (a = 1/3, b = 1, ε the
shallow plume's 2.5·10⁻³/m), into the first non-buoyant layer.
`scripts/plumeReach.mjs` repeats `cumulusColumn`'s ascent column by
column. It agrees with the model's top layer and share in all 7,911–7,943
plume columns of each state.

### On BOMEX's profile

The same ascent, source and w equation through Table B1's profile, with
the source the mean of the lowest 500 m. Siebesma et al.'s Figs. 7e and 10
give the LES's values: cloud base 0.5 km, core w 0.6 m/s there (the
cloud ensemble's too, about w*), core about 2.7 m/s at 1.5 km, core top
1.75–1.8 km and cloud top about 2.0 km over an inversion base at 1.48 km.

| plume | top, 20 m layers (km) | w at 1.5 km (m/s) | implied top on the model's levels (km) |
|---|---|---|---|
| the model's | 0.95 | 0 | 0.68 |
| w0 0.6 | 0.81 | 0 | 0.68 |
| a 1/3, b 2 (Gregory 2001), w0 0.6 | 0.75 | 0 | 0.68 |
| a 2/3, b 2 (Simpson and Wiggert 1969), w0 0.6 | 0.68 | 0 | 0.60 |
| ε = 1/z | 1.75 | 1.73 | 0.68 |
| the adiabat's source | 1.20 | 0 | 1.27 |
| ε = 1/z and the adiabat's source | 1.83 | 2.30 | 1.78 |
| Siebesma et al.'s eq. (15) (a 1/3, b 2, ε = 1/z, w0 0.6), the adiabat's source | 1.76 | 1.84 | 1.78 |

The adiabat is that of Siebesma et al.'s Fig. 7 (θ_l 298.8 K, q_t 17.25
g/kg), 0.6 g/kg moister than the source mean. Their eq. (15) solution
reaches about 1.75 km at about 2 m/s (Fig. 10), as the last row does, so
the setup holds. The model's plume is negatively buoyant from cloud base in
BOMEX: w falls from 1 to 0.28 m/s by 0.9 km. On the model's levels it ends
in its first cloudy layer.

The speed is not what is short. The model's w0 of 1 m/s is above the LES's
0.6 m/s at cloud base. It is also above EUREC4A's: Vogel et al. (2022)
give a cloud-base mass flux of 17.4 ± 7.5 mm/s and a cloud-base cover of
5.4 ± 3.1 percent, with C ≈ 2 C_core ∝ 2M/w*, so the core rises at about
w*. The model's b = 1 is half the b = 2 of Gregory (2001) and of Siebesma
and Holtslag. The LES fits make the reach shorter still. Tian et al.
(2019) find a ≈ 0.6–0.7 in their shallow-cumulus LES, with b varying by
more than a factor of three across datasets. de Roode et al. (2012) find
the pressure term the dominant sink, with a and b linearly related.

What is short is the buoyancy. Entraining at 2.5·10⁻³/m from the LCL and
leaving with the mean air below 500 m, the plume is never buoyant in
BOMEX. With either ε ∝ 1/z (the rate Siebesma et al. 2003 use, falling from
2·10⁻³/m at cloud base) or the updraft's own cloud-base air alone, it falls
short. With both, the same code gives the LES's core top and speeds.

ε = 1/z alone reaches 1.75 km on 20 m layers, but on the model's levels it
still ends in its first cloudy layer. The ascent stops at the first
non-buoyant cloudy layer, and at that layer's midpoint, 11 m above the
LCL, the plume has hardly condensed. On 20 m layers the same plume slows
to 0.78 m/s through the transition layer and then accelerates.

### In the trade boxes

The kinetic arm's states on days 1835, 1840 and 1845; box medians, the
range over the three days:

| | N Pacific trades | S Pacific trades | Atlantic trades | SE Pacific | sea 20–40S | sea 20–40N |
|---|---|---|---|---|---|---|
| plumes with no buoyant cloudy layer | 0.41–0.58 | 0.23–0.25 | 0.12–0.24 | 0.18–0.24 | 0.30–0.36 | 0.35–0.40 |
| largest speed, m/s | 1.00 | 1.25–1.54 | 1.08–1.39 | 1.31–1.36 | 1.00–1.04 | 1.00 |
| share | 0.19–0.31 | 0.42–0.44 | 0.09–0.26 | 0.16–0.26 | 0.23–0.25 | 0.17–0.23 |
| implied top, km | 1.31–1.51 | 1.32–1.50 | 1.52–1.53 | 1.51–1.77 | 1.45–1.53 | 1.45–1.56 |
| implied top less the inversion, km | −0.2 to −1.7 | −0.1 to −0.2 | 0 to −1.1 | −0.2 to −1.1 | −0.24 to −0.38 | −0.21 to −0.25 |
| ε = 1/z: no buoyant cloudy layer | 0.26–0.40 | 0.16–0.19 | 0.11–0.14 | 0.12–0.18 | 0.24–0.29 | 0.30–0.33 |
| ε = 1/z: implied top, km | 1.48–1.61 | 1.52–1.84 | 1.57–1.75 | 1.75–1.91 | 1.63–1.76 | 1.51–1.73 |
| ε = 1/z: less the inversion, km | 0 to −1.6 | +0.25 to +0.26 | +0.02 to −0.9 | 0 to −0.85 | −0.11 to −0.25 | −0.04 to −0.21 |

The inversion here is the interface of strongest θ_l stability between the
LCL and 4 km, at 1.5–2.9 km. In 12–58 percent of the plume columns the
plume saturates but is not buoyant in its first cloudy layer. There its
share is only the 1 m/s start's reach into that layer. Elsewhere it is
buoyant for one or two layers, at 1.0–1.5 m/s (upper quartiles 1.3–1.8).
w* is 0.70–0.85 m/s.

The literature's w equations lower the shares on day 1845, from 0.09–0.42
to:
- 0.05–0.24 with a 1/3, b 2;
- 0.05–0.17 with a 2/3, b 2 and w0 = w*.

ε = 1/z raises the tops by 0.05–0.37 km. The S Pacific trades' clouds then
reach 0.25 km into the inversion, as BOMEX's do.

So the reach is too short, and the water's placement is part of the
excess. A third to a half of the trade plumes detrain at cloud base
because they are not buoyant there, not because they rise too slowly.
Faster plumes would not move the water up. A return to share 1 would not
either: it detrains the same plumes one layer up. Three changes would carry
the water up toward the inversion:
- an entrainment rate that falls with height;
- a source that is the updraft's air;
- an ascent that w² carries through a non-buoyant layer, as the w equation
  itself does.

The next section tries them in the model. How large the source excess
should be is open: the IFS test parcel's excess (Cy43r1 eq. 6.19) gives
about 0.15 g/kg on BOMEX's fluxes, while the LES adiabat holds 0.6.

## The plume's buoyancy

Two options, on both engines and off by default (the defaults' fields are
bitwise unchanged):

- `cumulusEntrainmentHeight` c: above the LCL the shallow plume entrains
  at ε = c/z, z the layer midpoint's height above the surface, in place of
  `cumulusEntrainment`. Siebesma et al. (2003, section 6) find the LES's
  entrainment rate about 2·10⁻³/m near cloud base and falling with height,
  and with a fit ε ≈ 1/z their plume model reproduces the cores' q_t and
  θ_l. The plume detrains at ε + `cumulusDetrainment` −
  `cumulusEntrainment`, so its mass flux falls with height as it does with
  the constant rates. Siebesma et al. diagnose δ as the residual that the
  core mass flux's decrease requires, so the option changes the dilution
  alone.
- `cumulusAscent` 'kinetic' (with `cumulusOvershoot` 'kinetic'): the plume
  crosses a cloudy layer where it is not buoyant when its w² carries it
  through the whole layer (share 1). It ends in the first layer where w²
  does not, and detrains there by the share. This is the w equation's
  cloud top, "the height where the vertical velocity vanishes" (Siebesma
  et al. 2003), at the resolution of the layers.
- The updraft's air as the source is the existing `cumulusSource`
  'lowest', the lowest layer's air in place of the source layers' mean.

### On BOMEX's profile, in the model

`test/convection.test.mjs` places Table B1's profile over 1015 hPa. The
boundary layer is 500 m deep, the surface buoyancy flux 5.7·10⁻⁴ m²/s³
(BOMEX's 8·10⁻³ K m/s and 5.2·10⁻⁵ m/s) and u* 0.28 m/s. The test runs
`cumulusColumn` once and takes the top layer's base plus the share of its
depth (km):

| plume | bl34 | bl36 | bl42 |
|---|---|---|---|
| defaults | 0.92 | 0.92 | 0.92 |
| `cumulusAscent` 'kinetic' | 0.93 | 0.93 | 0.93 |
| `cumulusEntrainmentHeight` 1 | 1.80 | 1.80 | 1.69 |
| both | 1.80 | 1.80 | 1.69 |
| both, `cumulusSource` 'lowest' | 1.84 | 1.84 | 1.80 |
| LES | | core top 1.75–1.8, cloud top about 2.0 | |

The entrainment profile alone carries the plume to the LES's core top.
On these levels a plume entraining at 1/z is buoyant in its first cloudy
layer, so the carried ascent adds nothing on BOMEX. The test holds the
defaults below 1.1 km and the three 1/z plumes within 1.6–2.1 km.

### In the trade boxes

`plumeReach.mjs` follows both options. Under each set it agrees with the
model's top and share in every plume column (7,911–8,677 per state).
Median implied tops (km) on the defaults' day-1845 state, after one step:

| plume | N Pacific trades | S Pacific trades | Atlantic trades | SE Pacific | sea 20–40S | sea 20–40N |
|---|---|---|---|---|---|---|
| defaults | 1.31 | 1.32 | 1.52 | 1.51 | 1.45 | 1.45 |
| 'kinetic' ascent | 1.31 | 1.32 | 1.52 | 1.51 | 1.47 | 1.46 |
| 1/z | 1.48 | 1.52 | 1.57 | 1.75 | 1.62 | 1.53 |
| both | 1.48 | 1.71 | 1.57 | 1.83 | 1.68 | 1.57 |
| both, 'lowest' | 1.68 | 1.99 | 1.61 | 2.02 | 1.84 | 1.71 |

With 1/z the share of plumes that are not buoyant in their first cloudy
layer falls from 0.12–0.51 to 0.13–0.39 (0.01–0.31 with 'lowest'). In the
trades the carried ascent adds up to 0.19 km, where the 1/z plume meets a
non-buoyant layer it can cross. In each arm's own day-1845 state (the
arms below), the 1/z plume with the carried ascent tops at 1.69, 1.60,
1.26, 1.74, 1.54 and 1.57 km, against the defaults' 1.31–1.52.

### Arms at N=64

The protocol of the arms above. "1/z + ascent" is the two options
together, "both" the variance pair:

| arm | ASR | OLR | SWCRE | LWCRE | clear-sky OLR | net | rain mm/d |
|---|---|---|---|---|---|---|---|
| defaults (kinetic) | 229.7 | 226.1 | −60.5 | +33.7 | 259.8 | +3.6 | 3.12 |
| both | 234.2 | 227.4 | −55.9 | +32.5 | 260.0 | +6.8 | 3.02 |
| **1/z + ascent** | 234.3 | 226.6 | **−55.9** | +33.1 | 259.8 | +7.7 | 3.10 |
| 1/z + ascent, 'lowest' | 232.3 | 227.2 | −57.9 | +33.0 | 260.2 | +5.1 | 3.15 |
| 1/z + ascent, both | 237.4 | 226.9 | −52.9 | +32.9 | 259.9 | +10.5 | 3.03 |
| 1/z + ascent, 'lowest', both | 235.3 | 227.9 | −54.9 | +32.5 | 260.4 | +7.4 | 3.13 |

| arm | N Pacific trades | S Pacific trades | Atlantic trades | SE Pacific | sea 20–40S | sea 20–40N |
|---|---|---|---|---|---|---|
| defaults | 0.515 | 0.522 | 0.513 | 0.583 | 0.549 | 0.527 |
| both | 0.424 | 0.317 | 0.441 | 0.472 | 0.448 | 0.463 |
| 1/z + ascent | 0.466 | 0.498 | 0.295 | 0.350 | 0.406 | 0.403 |
| 1/z + ascent, 'lowest' | 0.493 | 0.498 | 0.299 | 0.326 | 0.411 | 0.380 |
| 1/z + ascent, both | 0.400 | 0.351 | 0.232 | 0.291 | 0.347 | 0.354 |
| 1/z + ascent, 'lowest', both | 0.489 | 0.454 | 0.262 | 0.189 | 0.353 | 0.373 |
| Earth | 0.2–0.4 | 0.2–0.4 | 0.2–0.4 | 0.6–0.8 | | |

The buoyant plume removes 4.6 of the kinetic share's 5.8 W/m², as much as
the variance pair, with LWCRE 0.6 W/m² lower. With the pair it removes 7.6,
1.8 more than share 1. The lowest-layer source puts back 2.0 W/m²
(−57.9) and raises BOMEX's top by only 0.04–0.11 km, so it is not taken.
Low cover falls in the Atlantic trades to 0.30 (0.23 with the pair) and
in the subtropical bands to 0.40–0.41 (0.35). The Pacific trade boxes
move less: 0.47 and 0.50, or 0.40 and 0.35 with the pair. The cost is the
SE Pacific box, from 0.58 to 0.35 (0.29 with the pair) against Klein and
Hartmann's 0.6–0.8 for the deck's core (their Peruvian region, 10–20S
80–90W). There the layers from 0.4 to 1.6 km lose most of their cover
(0.06–0.35 per layer to 0.00–0.16 over the three days) under plumes that
top higher (1.74 km against 1.51 in each arm's own day-1845 state). The
deck's gate is open on few of the box's columns, so the defaults' 0.58
there was itself resolved cumulus layers covered at about one half
(finding 2).

### The trade stress at N=128

The protocol above, with the 1/z plume and the carried ascent alone and
with the pair:

| τ_x, mN/m² | 140E–160E | 160E–180 | 180–160W | 160W–140W | 140W–120W | 120W–100W | 100W–80W |
|---|---|---|---|---|---|---|---|
| observed, about | −15 | −20 | −40 | −50 | −55 | −50 | −30 |
| defaults | −7.3 | −45.0 | −116.2 | −142.7 | −100.2 | −63.1 | −15.9 |
| both | −9.0 | −30.4 | −85.5 | −99.7 | −81.4 | −57.7 | −18.0 |
| 1/z + ascent | −10.2 | −32.6 | −110.5 | −132.4 | −99.8 | −74.0 | −21.4 |
| 1/z + ascent, both | −17.8 | −20.5 | −94.5 | −116.4 | −87.9 | −67.8 | −20.5 |

The buoyant plume leaves the mean of the three central bins where the
defaults have it, −102.1 against −102.0. It weakens 180–140W by 5–7
percent and strengthens 120W–100W by 17 percent. With the pair it gives
back half of the pair's gain: −90.7 against −79.6. The sea-level pressure
of 180–160W less 120W–100W's goes from −5.34 to −4.88 hPa (−4.44 with the
pair, against the pair's −3.85). Like the pair, the buoyant plume warms
the 120W–100W column by 0.7–0.9 K at 2.8–4 km, and the layers of 2–5 km
give back 0.77 hPa of the gradient (the pair 0.99). But that column stays
moist there (relative humidity 0.33–0.50, the pair's 0.20–0.39, the
defaults' 0.45–0.49), the layers below 2 km give nothing back (the pair
0.41), and those of 5–12 km take 0.44 hPa more (the pair gives back
0.33). Each arm is one realization, so how much of the 5–12 km part is
noise is not measured. Over 180–160W the cumulus layer's top is 1.78 km
(defaults 1.91, the pair 1.69).

Day means, days 6–20:

| arm | ASR | OLR | SWCRE | LWCRE | clear-sky OLR | net | rain mm/d |
|---|---|---|---|---|---|---|---|
| defaults | 235.8 | 239.1 | −60.3 | +29.0 | 268.1 | −3.3 | 3.08 |
| both | 239.6 | 239.6 | −56.4 | +28.5 | 268.1 | 0.0 | 3.01 |
| 1/z + ascent | 240.4 | 239.6 | −55.7 | +28.3 | 267.9 | +0.8 | 3.06 |
| 1/z + ascent, both | 243.4 | 240.1 | −52.6 | +28.2 | 268.3 | +3.3 | 3.00 |

### What the arms say

Entraining at 1/z, the model's shallow plume is BOMEX's plume: it tops
where the LES's cores do. In the arms it removes the kinetic share's
brightening as fully as the variance pair does, 4.6 W/m² at N=64 and at
N=128, by thinning the cloudy layers below the inversion. On its own it
is not the trade cumulus's fix. The trade stress does not move beyond
one realization's spread, alone or with the pair, and the SE Pacific box
loses cover it should keep. The variance pair remains the change that
both the cover and the stress support. Neither is made a default here.

Replicated (`docs/toa-radiation-bias.md`, finding 8: three realizations
each of the defaults, the pair and the pair with cloudColumns and
cloudSpread), one run of the protocol carries a standard deviation of
4.4 mN/m² in the central three bins' mean and 0.14 hPa in the gradient.
The pair takes 11.9 ± 3.6 mN/m² off the stress (−98.5 to −86.5), 12
percent rather than the single runs' 22, and 1.13 ± 0.11 hPa off the
gradient. The buoyant plume's single runs, −102.1 alone and −90.7 with
the pair, sit within one realization of those means.

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
- `test/convection.test.mjs` for the buoyant plume. The engine parity test
  also runs the carried ascent with ε = 1/z, and the carried ascent from
  the lowest layer with rain. Under them 169 and 138 of 362 random columns
  convect shallow (112 under the defaults). The engines agree to 1.7·10⁻⁴
  K and in every plume's top. A refusal test covers both options, and
  BOMEX's test is the one above. With the kinetic share's test, 6 of 6.
- The moist, plume and pinned CPU files (physics, moist, moistConvective,
  moistParcel, plumeColumn, plumePath, column, radiativeConvective,
  cloudEffects, cloudHeights, defaultRun): 115 of 115, the pinned digests
  holding with both options off. After the kernel change,
  `test/gpuModel.test.mjs` 32 of 32 and `test/boundaryLayer.test.mjs` 16 of
  16, with no command-buffer aborts.

## What remains

- The scale 5 against Mellor and Yamada's 3.2 (M22's choice): with the
  bounded length BOMEX's σ_s is 0.46 g/kg at 5 and 0.30 at 3.2.
- The SE Pacific deck's cover, which the resolved cloud does not stand for.
- The regime covers at N=128 (measured here at N=64 only), and the path
  by which the pair warms the 120W–100W column at 2–4 km.
- A coupled arm: the stress and the cover both act on the ocean.
- The buoyant plume's 17 percent stronger stress at 120W–100W, its moist
  east column at 2.8–4 km and its SE Pacific cover; whether the 5–12 km
  part of its pressure gradient is beyond the protocol's noise.
- A CPU/GPU field comparison at the page's configuration with the options
  of either change on.
- With `cloudColumns` and `cloudSpread` (`toa-radiation-bias.md`) the pair
  is the shortwave's remaining cloud-amount change; the net balance needs
  the clear-sky column first.
