# Observational references

Each reference the scorecards compare the model with is a file here that records where it comes from: the data set with its DOI or URL and license, the period, the box, the mask, the definition (including the model quantity that matches it), and the script that wrote it. A row of a scorecard takes its target from one of these files, never from a number typed into a script.

| file | quantity | source | period | script |
|---|---|---|---|---|
| `ceres-ebaf-global-monthly.json` | global TOA incoming solar, all-sky and clear-sky outgoing fluxes, cloud radiative effects | CERES-EBAF Ed4.2.1 (obs4MIPs), doi:10.5067/TERRA-AQUA-NOAA20/CERES/EBAF_L3B004.2.1, CC BY 4.0 | 2001–2024 | `scripts/referenceCeres.mjs` |
| `global-ts-berkeley.json` | global mean surface temperature by month | Berkeley Earth Land/Ocean, doi:10.5194/essd-12-3469-2020, CC BY-NC 4.0 | 1991–2020 | `scripts/referenceTs.mjs` |
| `nsidc-sea-ice-monthly.json` | Arctic and Antarctic sea-ice extent and area by month | NSIDC Sea Ice Index v4 (G02135), doi:10.7265/a98x-0f50, no restrictions | 1991–2020 | `scripts/referenceSeaIce.mjs` |
| `woa23-sst-1991-2020.bin` (`.json`) | monthly SST on 1° | World Ocean Atlas 2023, 1991–2020 normal (decav91C0), doi:10.25923/54bh-1613, public domain | 1991–2020 | `scripts/packWoaMonthlySst.py … decav91C0` |

Definitions that the earlier scorecards got wrong, as the Nov 7 audit found them:

- **Clear sky.** EBAF's `rsutcs` and `rlutcs` are the clear-sky fluxes of cloud-free areas. The model's clear-sky pass removes the clouds from every column, which is EBAF's total-region clear sky; Loeb et al. (2020) put that 2.2 W/m² lower in OLR and 0.5 W/m² higher in reflected sunlight. The CERES file's `rstcre` and `rltcre` are EBAF's effects against the total-region clear sky (annual −45.13 and 25.70 W/m², against −45.61 and 27.80 from the cloud-free fluxes); score the model against those.
- **SST period.** WOA23's `decav` (`data/woa_sst_monthly_1deg.bin`, which the ocean's initial state still uses) averages the atlas's decades and is centered near 1988. Score SST against the 1991–2020 normal here, 0.20 K warmer over 30S–30N.
- **Global surface temperature.** The 14.85 °C target had no source; the Berkeley Earth 1991–2020 climatology here is 14.76 °C for the year, with a seasonal cycle of 3.6 K.
- **Sea-ice extent.** The passive-microwave index counts no lake ice. The model's extent must count sea cells only; its FLake ice set the Arctic's March bias to +0.42 Mkm² where the sea ice alone gives +0.07.

Attribution: CERES data courtesy of NASA Langley Research Center's CERES team; Berkeley Earth data under CC BY-NC 4.0, attributed to Berkeley Earth (Rohde and Hausfather 2020); NSIDC's Sea Ice Index cited as Fetterer et al.; WOA23 from NOAA NCEI.
