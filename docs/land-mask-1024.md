# Land mask from an N = 1024 resampling: phase 1

Oct 24, 2026. Branch `claude/land-mask-1024`. Model defaults, sub-grid files, saved states and page defaults are untouched.

## What was built

- `js/fineGeography.module.js`
  - `relaxedMesh(N, relax)`: Grid's seeds, topology and Lloyd passes on typed arrays. Its centers, areas and neighbors match `new Grid(N, { relax })` bit for bit (checked at N = 16 and 64; `test/fineGeography.test.mjs`).
  - The file format of `data/geography1024.bin`.
  - `coarseGeography(mesh, fine)`: area-weighted fractions and means on a model mesh, plus each model edge's sill from the fine edges.
  - `connectedGeography(mesh, coarse, fine)`: the proposed connectivity rule.
- `scripts/fineGeography.mjs`: builds the N = 1024 dataset from ETOPO 2022's 60″ surface and bed GeoTIFFs (the cache that `scripts/landIce.mjs` reads).
- `scripts/landMaskCompare.mjs`: compares the new masks with `createGeography` on the 0.25° ETOPO 2022 raster, the fresh-start geography.

## Cost

| run | wall time | peak RSS | output |
|---|---|---|---|
| Grid(256, relax 30), object mesh | 5.2 s | 2.0 GB | |
| Grid(512, relax 30), object mesh | 25.5 s | 6.6 GB (6.4 GB heap) | |
| dataset N = 256 | 75 s | 3.1 GB | 62 MiB |
| dataset N = 512 | 75 s | 4.8 GB | 248 MiB |
| dataset N = 1024 | 108 s (mesh 44.5, 233M points 14, fine-edge sills 39) | 12.3 GB | 990 MiB |
| comparison, N = 16 to 192 | 129 s | 4.8 GB | |

At N = 1024 the object Grid would need about 25 GB of heap. The typed-array mesh peaks at 8.9 GB, most of it the unrelaxed Grid used as the seed.

## The dataset

`data/geography1024.bin` (gitignored, 990 MiB) holds 10,485,762 cells and 31,457,280 edges. Each 60″ point goes to its nearest cell, weighted by its area.

- Per cell:
  - center, area and neighbors;
  - land fraction;
  - mean surface, mean land surface and mean sea surface;
  - bed and deepest point;
  - ice fraction and thickness;
  - floating fraction, thickness, surface and bed, by `scripts/landIce.mjs`'s rules (60″ points, 10 m minimum, 15 m firn).
- Per edge: the sill, the depth at which the deepest 60″ points of the two cells join through their own points (eight-connected).

Every fine cell holds at least one point; no point lacks data. The total land area is 148.36 million km².

## Per-N comparison

Three masks:
- **today**: `createGeography` with ETOPO 2022 + observed ice. Land means more than half of the 0.25° points; LAND_BRIDGES and SEA_STRAITS apply.
- **area**: area-weighted land fraction > 0.5 from N = 1024.
- **rule**: area plus the connectivity rule.

Box columns give the land cells of today / area / rule, out of the box's cells.

| N | land cells, today / area / rule | land, million km² | today→area: to sea, to land | rule flips: to sea, to land | sea cells with mean z > 0, today / area | elevation RMS vs today, all / land (m) | Darién 7–10N 80–76.5W | Cuba box | SF box |
|---|---|---|---|---|---|---|---|---|---|
| 16 | 735 / 726 / 725 | 146.30 / 144.48 / 144.28 | 10, 1 | 7, 6 | 23 / 26 | 57 / 13 | 1/0/1 of 1 | 0/0/0 of 2 | – |
| 32 | 2939 / 2942 / 2944 | 146.39 / 146.54 / 146.63 | 14, 17 | 8, 10 | 75 / 78 | 23 / 13 | 2/1/2 of 3 | 1/2/3 of 9 | – |
| 64 | 11849 / 11841 / 11848 | 147.56 / 147.46 / 147.55 | 71, 63 | 8, 15 | 185 / 220 | 28 / 19 | 5/5/5 of 12 | 7/8/8 of 35 | – |
| 128 | 47258 / 47489 / 47493 | 147.13 / 147.84 / 147.86 | 237, 468 | 18, 22 | 656 / 553 | 60 / 39 | 14/16/16 of 41 | 25/34/34 of 144 | 2/2/2 of 3 |
| 192 | 106056 / 107074 / 107052 | 146.74 / 148.15 / 148.12 | 489, 1507 | 49, 27 | 1430 / 866 | 79 / 48 | 31/34/34 of 91 | 59/79/78 of 318 | 6/6/6 of 7 |

Today's rule loses land as N grows: 1.6 million km² short of the 60″ total at N = 192. A cell at exactly half its points counts as sea, and at N ≥ 128 a cell has only 3 to 6 points. The area rule stays within 0.2–0.5 million km² of the total at N ≥ 128.

**Panama, N = 128.** Two of the three Darién sea cells become land:
- 8.86N 77.80W: half of the points today → 0.73 of the area;
- 8.45N 78.08W: half → 0.86.

The third, 9.33N 78.09W (0.40 → 0.33), is the San Blas coast, mostly sea, and stays sea. The ocean is closed across the isthmus at every N in all three masks, and at N = 1024.

**Cuba, N = 128.** 9 cells become land (25 → 34). Sea cells with a positive mean elevation in the box fall from 7 to 2.

**San Francisco, N = 128.** No change. The Golden Gate cell is 0.39 land and stays sea; the Bay cell is 0.79 and stays land.

## Strait connectivity

The table gives the widest ocean path's sill in meters: the deepest level at which water still connects the two basins.
- 0: closed.
- "one cell": both basins fall in the same model cell.
- Coarse paths are confined to a box around the strait, widened by 1.5 cell spacings.
- Today's values include SEA_STRAITS, which forces 300 m at Gibraltar, 150 at Bab-el-Mandeb and 90 at Hormuz.

| strait (N = 1024) | N=16 today / area / rule | N=32 | N=64 | N=128 | N=192 |
|---|---|---|---|---|---|
| Gibraltar (341) | 0 / 0 / one cell | 300 / 341 / 341 | 300 / 341 / 341 | 300 / 318 / 318 | 300 / 341 / 341 |
| Bab-el-Mandeb (137) | 0 / 0 / 137 | 150 / 0 / 137 | 150 / 0 / 169 | 104 / 59 / 59 | 99 / 137 / 137 |
| Hormuz (83) | 0 / 0 / 83 | 90 / 83 / 83 | 90 / 83 / 83 | 69 / 83 / 83 | 47 / 83 / 83 |
| Bering (44) | 0 / 39 / 39 | 44 / 39 / 39 | 44 / 30 / 30 | 44 / 44 / 44 | 44 / 44 / 44 |
| Lombok (216) | one cell everywhere | one cell | 336 / 216 / 216 | 336 / 216 / 216 | 336 / 216 / 216 |
| Panama (0) | 0 / – / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |

The N = 1024 sills agree with the literature at Bab-el-Mandeb (Hanish sill, 137 m), Hormuz and Bering. Gibraltar reads 341 m against the Camarinal sill's ~290 m: the 60″ raster is coarser than the sill.

## The proposed rule (`connectedGeography`)

1. **Class.** A cell is land where its area-weighted land fraction exceeds 0.5.
2. **Edge sills.** A model edge's sill is the depth at which the deepest fine cells of its two cells join through fine edges between those cells. An ocean edge is open only where that sill is above 0, so water crosses at fine scale.
3. **Straits.**
   - Resolved water is fine cells with a land fraction ≤ 0.5, joined by fine edges with sill > 0.
   - Where resolved water joins model ocean components that the open edges keep apart, the model land cells on the joining fine paths become sea.
   - Paths come from a maximum spanning tree over the fine edges, with each model ocean component as one node, pruned to the branches that join components. Weight: the sill capped at `openSill` = 150 m, then the least land in the model cells the edge lies in.
   - A flipped cell is at least as deep as the sill through it, and none of its edges is deeper than that sill.
4. **Bridges.** The same construction applies to land:
   - resolved land is fine cells with a land fraction > 0.5 that touch;
   - model sea cells on the paths become land, preferring the model cells with the most land;
   - a cell that a strait needs stays sea.

This replaces LAND_BRIDGES and SEA_STRAITS.

### Cases it changes

**N = 64:**
- To sea: Bab-el-Mandeb (12.7N 43.0E); the Gulf of Bothnia and Gulf of Finland throats; the White Sea; Dease and Queen Maud straits.
- To land: Sulawesi's arms; Jutland and Fehmarn; Novaya Zemlya; the Gulf of California head; Cape Breton; Auckland; Panama at 8.2N 80.1W.

**N = 128:**
- To sea:
  - Europe and Russia: the Great Belt (55.4N 11.4E), the White Sea, the Ob and Taz gulfs, Khatanga.
  - Alaska and Canada: Cook Inlet, Chatham Strait, Juan de Fuca, Hudson Strait and Frobisher Bay (one cell at 0.91 land), Canadian Arctic channels.
  - Río de la Plata.
- To land:
  - Asia and the western Pacific: Kanmon (Kyushu–Honshu), Hokkaido, Palawan and Sulawesi pieces, Bali.
  - The Alaska Peninsula.
  - Elsewhere: Calabria, New Caledonia, Auckland, Valdés.

**N = 192:** 76 flips. Besides those above:
- To sea:
  - Europe and Russia: the Gulf of Corinth; the Volga delta and the north Caspian; the Yenisei.
  - Asia: Tatar Strait; the Seto Inland Sea; the Gulf of Khambhat; Halmahera.
  - The Middle East: the Gulf of Aqaba.
  - The Americas: the St. Lawrence; the Amazon mouth; Lagoa dos Patos; Chiloé.
  - The Arctic: a Svalbard fjord.
- To land: the Mississippi delta, Andros, the Solomons.
- Three conflicts, kept as sea: north Caspian 44.5N 51.1E, Chatham Strait, Halmahera.

## Open points for Joshua

1. **Narrow features.** The rule carves every feature resolved at 7 km into the model mesh, including estuaries and lagoons that leave a one- to four-cell coarse sea component. This conflicts with the unresolved-bathymetry ruling. A minimum size for the smaller component, or turning such components to land or lakes, is his call.
2. **Edge-sill limit.** The deepest-to-deepest edge sill can miss the crossing inside one cell that holds both basins:
   - Bab-el-Mandeb at N = 64 reads 169 m against 137 m.
   - It can also read a shelf shallower than the channel: Bab-el-Mandeb at N = 128 gives 59 m, Bering at N = 64 gives 30 m.

   Both limits apply to today's rule too.
3. **Porous barriers.** Adcroft (2013) open widths are not built. The fine edges give the data for them.
4. **Not yet carried to the model.** The ice fields and `iceShelves` 'sea' / 'cavity' (the floating ice's bed) are not passed through `createGeography`. Sub-grid files, sills, cavities and saved states would all need rebuilding at every N.
