# The tropical upper stratosphere's ascent

## The question

Finding 6 of `docs/clear-sky-olr-deficit.md` left the circulation's
question open. At day 1855 of the 30-day fixed-SST arm (N=64, from
`runs/thirteen/thirteen64_day1825.bin`, SST=woa) 15° S–15° N is 7–8 K
colder than AFGL's tropical atmosphere at 1.6–7 hPa. The clear-sky
radiative heating at the model's own temperatures, +1.88 K/day at
1.59 hPa and +1.00 at 3.54, says the ascent's adiabatic cooling holds
the layers there, and without the non-orographic gravity-wave drag
1.6 hPa warms 3.4 K in 30 days. Is that ascent too strong against a
reanalysis's residual circulation? If it is, what drives it: what the
drag deposits at the lid, or what breaks below it?

## Instruments

- `scripts/residualCirculation.mjs`: the transformed Eulerian mean (TEM)
  of saved states on pressure surfaces (Andrews, Holton and Leovy 1987),
  each column interpolated in ln p from the layer midpoints σπ, ω from
  the interfaces' πσ̇ and σ(π_t + v·∇π), zonal means and eddy
  covariances over 2.5° bands, w* = −Hω*/p with H = 7 km and the
  residual ascent's adiabatic cooling ω*θ̄_p T̄/θ̄.
- `scripts/fixedSst.mjs` with `TEM_HOURS`: the same TEM sampled every
  three model hours during an arm, with the gravity-wave drag's zonal
  force and absolute flux by layer from the CPU engine on the same state.
- ERA5 (Hersbach et al. 2020) from ARCO-ERA5's raw pressure-level files
  (public; one variable, one level and 24 hours of 0.25° int16 per
  file), range reads of the rows within 20° of the equator at 1, 2, 3, 5,
  7, 10, 20, 30, 50 and 70 hPa, and the same TEM on each sample over
  1440 longitudes, the eddy term's φ-derivative over ±1°.
- The equatorial wind's annual cycle from `runs/verda-eleven/eleven64.log`,
  whose daily `upper winds` line has the 5° S–5° N zonal-mean wind of
  every layer above 200 hPa for five years (bl36, whose layers above
  100 hPa are bl42's; the gravity waves, sponge and lid friction as now).

Both the model and ERA5 carry a zonally symmetric diurnal and semidiurnal
signal in the tropical ω̄ at 1–5 hPa, so a mean from one hour of the day
is biased. ERA5's means over 2007–2019 at 00 and at 12 UTC alone differ
by 3.2 mm/s at 1 hPa and 2.1 at 3 hPa, and the pair 00 and 12 UTC still
sits 0.5–0.65 below the four-hour mean at 1–3 hPa. The model's states,
saved at 00 UTC, put w* 0.3–0.5 mm/s below its three-hourly mean at
1.5–5 hPa over the same 30 days. The means below take four (ERA5) or
eight (the model) samples a day.

## Findings

### 1. In the arm's weeks the ascent at 1–3 hPa is 1–4 mm/s too strong

15° S–15° N; the model is the control arm (days 3–30, 224 samples), and
ERA5 is 2007–2019 on March 23 and 29 and April 4, 10 and 16 at 00, 06,
12 and 18 UTC (260 samples), the same weeks after the equinox. The wind
is 5° S–5° N.

| hPa | w*, mm/s | | cooling, K/day | | T, K | | wind, m/s | |
|---|---|---|---|---|---|---|---|---|
| | model | ERA5 | model | ERA5 | model | ERA5 | model | ERA5 |
| 1 | 3.88 | 0.05 | 3.15 | 0.02 | 253.2 | 265.6 | −20.0 | +30.5 |
| 1.5 | 2.30 | | 2.12 | | 256.4 | | +16.3 | |
| 2 | 1.60 | −0.15 | 1.60 | −0.18 | 253.5 | 264.5 | +34.3 | +13.3 |
| 3 | 0.89 | −0.06 | 0.97 | −0.10 | 247.8 | 256.9 | +56.4 | −6.8 |
| 5 | 0.46 | 0.28 | 0.48 | 0.32 | 239.2 | 244.6 | +56.5 | −17.8 |
| 7 | 0.18 | 0.38 | 0.19 | 0.37 | 233.2 | 237.9 | +47.9 | −15.8 |
| 10 | −0.01 | 0.43 | −0.01 | 0.40 | 227.2 | 232.5 | +27.8 | −15.7 |
| 20 | 0.19 | 0.32 | 0.18 | 0.33 | 215.9 | 222.2 | +3.2 | −8.6 |
| 50 | 0.24 | 0.26 | 0.23 | 0.27 | 200.5 | 205.0 | +11.1 | +0.3 |

The heat budget closes: the ascent's adiabatic cooling, about 2.0 K/day
at 1.59 hPa and 0.75 at 3.54, against finding 6's radiative heating of
+1.88 and +1.00 at the same temperatures.

In these weeks ERA5 has no residual ascent at 1–3 hPa. The
semiannual oscillation's westerly phase (+31 m/s at 1 hPa, +13 at 2)
lies over easterlies of 7–18 m/s at 3–10 hPa, and its westerly shear
holds a warm stratopause with the ascent cancelled. The model's weeks
have the opposite structure. Westerlies of 48–62 m/s at 3–7 hPa lie
under easterlies above 1 hPa (−20 m/s at 1 hPa, −53 at 0.5), so the
easterly shear at 1–2 hPa carries the column's strongest ascent. The
ascent is 1.7–3.8 mm/s too strong at 1–2 hPa and 0.95 too strong at
3 hPa; at 7–10 hPa it is 0.2–0.4 too weak.

### 2. What drives it: the drag that the 3–7 hPa westerlies filter

The westerlies at 3–7 hPa take the eastward waves slower than 50–60 m/s
at their critical levels, so what rises past 3 hPa is mostly westward.
The arm's zonal-mean drag by layer, m/s/day:

| layer, hPa | equator | 7.5° N, S | 15° N, S | 30° N, S |
|---|---|---|---|---|
| 0.15 and 0.65 (the lid layers, one acceleration) | −7.00 | −5.52, −4.56 | −0.71, −0.06 | +0.18, −1.74 |
| 1.61 | −1.44 | −2.20, −2.27 | −2.77, −2.56 | +1.24, +2.38 |
| 3.59 | +0.29 | +0.37, +0.42 | +0.35, +0.51 | +0.32, +0.77 |
| 7.49 | +0.86 | +0.91, +0.98 | +0.65, +0.76 | +0.09, +0.08 |

At the equator 4.4 mPa of absolute flux (both directions) leaves the
54 hPa layer upward, 2.8 mPa leaves the 7.5 hPa layer, and 1.77 mPa
reaches 1 hPa, where all of it is spread over the two lid layers.

By downward control (Haynes et al. 1991), the upward mass flux through p
over |φ| < 15° is the poleward residual mass flux above p across ±15°.
In a steady state that flux is −X/f̂ for each zonal force X, with
f̂ = f − (ū cos φ)_φ/(a cos φ). The drag's share, 10⁸ kg/s (with f
alone in parentheses):

| hPa | TEM ascent | the lid layers' deposit | the breaking below the lid | the drag in all |
|---|---|---|---|---|
| 1 | 7.47 | 5.65 (0.93) | 0 | 5.65 |
| 1.5 | 6.63 | 5.65 (0.93) | 2.93 (3.22) | 8.58 |
| 2 | 6.17 | 5.65 (0.93) | 5.87 (6.44) | 11.52 |
| 3 | 5.13 | 5.65 (0.93) | 6.67 (7.07) | 12.32 |
| 5 | 4.40 | 5.65 (0.93) | 5.35 (4.95) | 11.00 |
| 7 | 2.44 | 5.65 (0.93) | 3.00 (1.53) | 8.65 |
| 10 | −0.16 | 5.65 (0.93) | −0.53 (−3.59) | 5.11 |

The breaking at 1.6 hPa drives as much of the ascent at 2–5 hPa as the
whole TEM ascent; the resolved waves and the lid friction take some of it
back. The lid's share hangs on f̂. The lid layers' easterlies (−53 m/s at
0.5 hPa) shear so strongly toward 15° that f̂ there is a fraction of f,
and the steady balance is ill-conditioned. The arms decide it.

Three more 30-day arms from the same state take the drag away or move
its lid. Means over days 3–30, w* (mm/s) and T (K) over 15° S–15° N:

| hPa | control | no drag (`GRAVITY_WAVES=false`) | `lidTests: true` | one lid layer (`lidPressure: 0`) | ERA5, the same weeks |
|---|---|---|---|---|---|
| 1 | 3.88, 253.2 | 2.09, 255.0 | 5.05, 252.6 | 5.26, 252.6 | 0.05, 265.6 |
| 1.5 | 2.30, 256.4 | 1.44, 258.4 | 2.80, 255.5 | 2.83, 255.6 | |
| 2 | 1.60, 253.5 | 1.05, 255.3 | 1.84, 252.8 | 1.87, 252.8 | −0.15, 264.5 |
| 3 | 0.89, 247.8 | 0.66, 249.0 | 0.94, 247.2 | 0.94, 247.0 | −0.06, 256.9 |
| 5 | 0.46, 239.2 | 0.36, 239.5 | 0.49, 239.0 | 0.51, 238.9 | 0.28, 244.6 |
| 7 | 0.18, 233.2 | 0.26, 232.8 | 0.20, 233.2 | 0.21, 233.2 | 0.38, 237.9 |
| 10 | −0.01, 227.2 | 0.20, 226.6 | 0.01, 227.2 | 0.02, 227.3 | 0.43, 232.5 |

The drag drives 46 percent of the ascent at 1 hPa, where only the lid
layers' deposit acts above, and a third of it at 1.5–2 hPa. Without it
the ascent through 1 hPa falls by 3.5 × 10⁸ kg/s, between the f and f̂
estimates of the lid's share. MiMA's own algorithm, the waves tested in
the lid layers too, and a lid of the top layer alone both make the
ascent stronger: 1.2–1.4 mm/s more at 1 hPa and 0.2–0.5 more at
1.5–2 hPa, with those layers 0.6–0.9 K colder. Tested in the 0.65 hPa
layer, the westward waves stop at their critical levels there. Its
easterlies grow from −52 to −59 m/s, and the wind at 1–1.5 hPa below it
turns 7 m/s easterly with them. Even
with no drag, the ascent at 1–2 hPa is 1.1–2.1 mm/s in these weeks. That
ascent rides on the easterly shear above a westerly jet of 55 m/s at
3–4 hPa, which 30 days without the drag do not remove.

### 3. The semiannual oscillation is close; a two-year mode sits too high

5° S–5° N, the model from eleven's daily winds over days 366–1825, ERA5
over 2007–2019 (four dates a month at 00 UTC) interpolated in ln p to the
model's layers. The semiannual harmonic of the monthly composite (m/s,
and the model day of its westerly maximum, day 0 the March equinox),
the standard deviation of the monthly anomalies from that composite,
and the mean:

| layer, hPa | semiannual, model | ERA5 | anomalies' sd, model | ERA5 | mean, model | ERA5 |
|---|---|---|---|---|---|---|
| 0.642 | 22.7, day 6 | | 26.0 | | −8.1 | |
| 1.59 | 20.7, day 43 | 24.5, day 37 | 36.7 | 13.5 | +0.8 | −5.6 |
| 3.54 | 10.9, day 64 | 14.9, day 59 | 33.5 | 14.2 | +6.9 | −10.8 |
| 7.39 | 2.7 | 4.3, day 67 | 32.1 | 16.5 | +20.2 | −11.6 |
| 14.0 | 2.2 | 1.7 | 21.9 | 18.8 | +33.2 | −10.2 |
| 23.9 | 0.5 | 1.2 | 13.8 | 18.5 | +36.1 | −7.6 |
| 37.2 | 0.9 | 1.0 | 9.9 | 14.5 | +23.3 | −3.5 |
| 53.1 | 0.6 | 0.8 | 6.8 | 10.9 | +7.8 | −0.2 |

The semiannual oscillation has 73–85 percent of ERA5's amplitude at
1.6–3.5 hPa and its westerly maxima come 5–6 days later. ERA5's own SAO is
stronger than SABER's and MLS's (Kawatani et al. 2020; the temperature
SAO at 2 hPa is 4.3 K in SABER, 4.2 in MLS and 4.7 in ERA5). The SAO is
not the error.

The error is a two-year mode, the periodogram's strongest peak at
1.6–53 hPa (740–785 days, 44–50 m/s at 1.6–7 hPa). Its anomalies' standard
deviation is 1.9–2.7 times ERA5's at 1.6–7 hPa and 0.6–0.75 times it at
24–53 hPa, where ERA5's quasi-biennial oscillation peaks (its 966-day
peak: 21 m/s at 20 hPa). Its mean is 6–44 m/s more westerly than ERA5's
at 1.6–37 hPa, where ERA5's is easterly. The arm's weeks fall in its
westerly phase at 3–7 hPa. MiMA, whose cg_drag settings the model takes
(B_eq = B₀ = 4.3 mPa, what reaches the lid spread over the layers above
0.85 hPa), keeps a QBO with 1.2 km layers in its lower stratosphere (40
levels). Its QBO's period and its amplitude in the lowermost stratosphere
grow steadily from 30 to 120 levels, and only the 120-level run (400 m)
has the observed asymmetry between the westerly and easterly regimes
(Garfinkel et al. 2022). CAM5 and the GISS model need about 500 m for a
realistic QBO (Richter et al. 2014; Geller et al. 2016). bl42's layers
are 1.1–1.7 km thick at 62–110 hPa and 2.3–5.7 km from 62 up to 1 hPa.

### 4. Over the mode's phases only the layer under the lid ascends too fast

Four more arms of 10 days each start from thirteen's states at days
1460, 1551, 1643 and 1734 (the equinoxes and solstices of year 5). Each
samples the TEM every three hours from day 1. With the control they give
w* (mm/s) and the 5° S–5° N wind (m/s), and their mean against ERA5's
2009–2018 annual mean (960 samples, two dates a month at four hours):

| hPa | day 1460 | 1551 | 1643 | 1734 | 1825 (control) | five-phase mean, w*, T | ERA5, w*, T |
|---|---|---|---|---|---|---|---|
| 1 | 0.35, +46 | 3.17, −20 | 1.23, +49 | 3.93, +6 | 3.88, −20 | 2.51, 255.1 | 0.46, 263.9 |
| 1.5 | −0.52, +18 | 1.87, +1 | −0.04, +49 | 2.24, +43 | 2.30, +16 | 1.17, 258.9 | |
| 2 | −1.02, 0 | 1.28, +5 | −0.70, +27 | 1.41, +50 | 1.60, +34 | 0.52, 255.8 | 0.63, 260.7 |
| 3 | −0.25, −26 | 0.67, +5 | −0.39, −13 | 0.66, +52 | 0.89, +56 | 0.32, 249.4 | 0.59, 252.5 |
| 5 | 0.76, −19 | 0.24, 0 | 0.48, −33 | −0.02, +26 | 0.46, +56 | 0.38, 239.0 | 0.55, 241.4 |
| 7 | 0.61, −2 | 0.33, −5 | 0.65, −37 | 0.20, +1 | 0.18, +48 | 0.39, 231.5 | 0.50, 235.6 |
| 10 | 0.46, +24 | 0.39, +4 | 0.87, −16 | 0.32, −8 | −0.01, +28 | 0.41, 225.2 | 0.47, 230.7 |
| 20 | 0.35, +47 | 0.10, +32 | 0.50, +28 | 0.51, +6 | 0.19, +3 | 0.33, 214.8 | 0.37, 221.5 |
| 50 | 0.22, +8 | 0.08, +21 | 0.28, +24 | 0.25, +19 | 0.24, +11 | 0.21, 201.9 | 0.23, 206.2 |

The ascent at 1–3 hPa swings from −1.0 to +3.9 mm/s with the
equatorial wind. Where westerlies lie above easterlies, as at days 1460
and 1643,
the layers descend and 1.5–2 hPa reaches 260–264 K, ERA5's April. Over
the five phases the ascent is too strong only at 1–1.5 hPa, the layer
under the lid. There 4.8 × 10⁸ kg/s rises into the lid layers over
|φ| < 15°, against ERA5's 0.9 × 10⁸ through 1 hPa. At 2–10 hPa the
ascent is 0.5–0.9 times ERA5's, yet the layers are 2.4–5.5 K colder, so
there the cold bias is not the ascent's.

## What this means

- Finding 6's ascent is too strong against ERA5's residual circulation
  for the same weeks, by 1.7–3.8 mm/s at 1–2 hPa and 0.95 at 3 hPa. At
  1.5–3 hPa most of that, and at 1 hPa half of it, is the phase the arm
  started in: the westerly phase at 3–7 hPa of a two-year mode whose
  counterpart in ERA5, the QBO, peaks 2–3 scale heights lower at half
  the amplitude. The drag that the mode's westerlies let through drives
  a third to a half of the ascent at 1–2 hPa.
- Over the mode's phases only the layer under the lid ascends too fast:
  2.51 mm/s at 1 hPa against ERA5's 0.46. In the arm's weeks the drag
  drives half of the ascent there, but MiMA's own treatment (the waves
  tested in the lid layers) and a thinner lid make it stronger. The
  waves that pass 1 hPa belong to a mesosphere that bl42 does not have,
  and no reference says where else in the column their momentum should
  go.
- At 2–10 hPa the mean ascent is 0.5–0.9 times ERA5's, while the layers
  are 2.4–5.5 K colder than ERA5's annual mean. That cold bias lies with
  the radiation and the other terms of the heat budget, not with the
  residual ascent.
- The SAO is close: 73–85 percent of ERA5's amplitude at 1.6–3.5 hPa and
  5–6 days late, and ERA5's own SAO is stronger than the satellites'. A
  QBO driven by parameterized waves needs layers of 0.4–1.2 km in the
  lower stratosphere (Garfinkel et al. 2022; Richter et al. 2014; Geller
  et al. 2016), and bl42's above 62 hPa are 2–5 times as thick as MiMA's.

No option comes of this finding. The change the references point to is
the stratosphere's vertical resolution. A level set with layers near
1.2 km from 100 to 1 hPa would have to be judged over several years
against ERA5's QBO (the anomalies' standard deviation and the mean by
level), its SAO, and the five-phase w* here. That is more than a 30-day
arm can show. Before the lid is touched, the absolute momentum flux
that reaches 1 hPa (1.77 mPa at the equator in the arm's weeks) should
be measured against the satellite climatologies (Ern et al. 2018).

## References

- Alexander, M. J., and T. J. Dunkerton, 1999: A spectral
  parameterization of mean-flow forcing due to breaking gravity waves.
  J. Atmos. Sci., 56, 4167–4182.
- Andrews, D. G., J. R. Holton, and C. B. Leovy, 1987: Middle Atmosphere
  Dynamics. Academic Press, 489 pp.
- Ern, M., Q. T. Trinh, P. Preusse, J. C. Gille, M. G. Mlynczak, J. M.
  Russell III, and M. Riese, 2018: GRACILE: a comprehensive climatology
  of atmospheric gravity wave parameters based on satellite limb
  soundings. Earth Syst. Sci. Data, 10, 857–892.
- Garfinkel, C. I., E. P. Gerber, O. Shamir, J. Rao, M. Jucker, I. White,
  and N. Paldor, 2022: A QBO cookbook: sensitivity of the quasi-biennial
  oscillation to resolution, resolved waves, and parameterized gravity
  waves. J. Adv. Model. Earth Syst., 14, e2021MS002568.
- Haynes, P. H., C. J. Marks, M. E. McIntyre, T. G. Shepherd, and K. P.
  Shine, 1991: On the "downward control" of extratropical diabatic
  circulations by eddy-induced mean zonal forces. J. Atmos. Sci., 48,
  651–678.
- Geller, M. A., et al., 2016: Modeling the QBO—Improvements resulting
  from higher-model vertical resolution. J. Adv. Model. Earth Syst., 8,
  1092–1105.
- Hersbach, H., et al., 2020: The ERA5 global reanalysis. Q. J. R.
  Meteorol. Soc., 146, 1999–2049.
- Kawatani, Y., T. Hirooka, K. Hamilton, A. K. Smith, and M. Fujiwara,
  2020: Representation of the equatorial stratopause semiannual
  oscillation in global atmospheric reanalyses. Atmos. Chem. Phys., 20,
  9115–9133.
- Richter, J. H., A. Solomon, and J. T. Bacmeister, 2014: On the
  simulation of the quasi-biennial oscillation in the Community
  Atmosphere Model, version 5. J. Geophys. Res. Atmos., 119, 3045–3062.
- Shepherd, T. G., and T. A. Shaw, 2004: The angular momentum constraint
  on climate sensitivity and downward influence in the middle atmosphere.
  J. Atmos. Sci., 61, 2899–2908.
