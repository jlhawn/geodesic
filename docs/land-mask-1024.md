# Land mask from an N = 1024 resampling: phase 1

Oct 24, 2026. Branch `claude/land-mask-1024`. Model defaults, sub-grid files, saved states and page defaults are untouched.

## What was built

- `js/fineGeography.module.js`
  - `relaxedMesh(N, relax)`: Grid's seeds, topology and Lloyd passes on typed arrays. Its centers, areas and neighbors match `new Grid(N, { relax })` bit for bit (checked at N = 16 and 64; `test/fineGeography.test.mjs`).
  - The file format of `data/geography1024.bin`.
  - `oceanDepth(fine)`: how deep each fine cell's water stays joined to the world ocean.
  - `coarseGeography(mesh, fine)`: area-weighted fractions and means on a model mesh, each cell's representative fine water, and each model edge's sill.
  - `connectedGeography(mesh, coarse, fine)`: the connectivity rule.
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
| comparison, N = 16 to 192 | 135 s | 4.7 GB | |

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

## The rule (`coarseGeography` + `connectedGeography`)

1. **Class.** A cell is land where its area-weighted land fraction exceeds 0.5.

2. **Edge depth.** Exchange through a strait is controlled at its sill (Armi & Farmer 1986, J. Fluid Mech. 164; Bryden & Kinder 1991, Deep-Sea Res. 38). Adcroft (2013, Ocean Modelling 67, 13–27) derives model face depths from high-resolution topography the same way. Here:
   - Each model cell's water is represented by one fine cell: the one whose water stays joined to the world ocean deepest (`oceanDepth`, a widest path from the deepest fine cell). Lakes have no such join, so for them the deepest fine cell is used.
   - A model edge's `edgeSill` is the depth to which the representatives of its two cells stay joined through fine edges.
     - The search covers the fine cells of the two cells and of the two cells at the edge's ends, so a channel that leaves through a third cell counts.
     - Where land parts them at the fine scale, the edge is 0.
   - Where the edge comes out shallower than both representatives' depth to the ocean, the search widens to every neighbor of the two cells. This can only deepen the edge.
   - Since every edge equals a fine bottleneck between representatives, no model path can be deeper than the fine water between its ends (the ultrametric inequality). Over-reading is ruled out by construction.
   - An ocean edge is open where `edgeSill` > 0.

3. **Minimum size.** Model ocean components are sea cells joined by open edges.
   - A component of fewer than 3 cells that resolved fine water joins to a component of 3 or more is a feature the mesh does not resolve. It turns to land, under the Oct 18 ruling that features the mesh cannot resolve are ignored. Resolved fine water means fine cells ≤ 0.5 land, joined by fine edges with sill > 0.
   - A component the fine water joins to no other is left as it is (see Lakes and depressions).

4. **Straits.** Where resolved fine water joins components of 3 or more cells that the model keeps apart, the model land cells on the joining paths turn to sea.
   - The paths come from a maximum spanning tree over the fine edges, with each component as one node, pruned to the branches that join components.
   - Weight: the sill capped at 150 m, then the least land in the model cells.
   - A flipped cell is at least as deep as its deepest edge.

5. **Bridges.** These keep no minimum size.
   - Where resolved fine land (fine cells > 0.5 land that touch) joins model land components kept apart by sea cells, the sea cells on the paths through the most land turn to land.
   - A cell that a strait needs stays sea. There are no conflicts at any N now.

The rule replaces LAND_BRIDGES and SEA_STRAITS. `test/fineGeography.test.mjs` checks it at N = 16 and 64 on synthetic geographies:
- a two-cell estuary that a channel joins to the ocean turns to land;
- a three-cell sea is joined through its 100 m channel;
- a 50 m sill inside one cell controls the model path;
- a 200 m channel through a third, land, cell joins two cells that a wall parts;
- the wall alone keeps them apart.

## Per-N comparison

Three masks:
- **today**: `createGeography` with ETOPO 2022 + observed ice. Land means more than half of the 0.25° points; LAND_BRIDGES and SEA_STRAITS apply.
- **area**: area land fraction > 0.5 from N = 1024.
- **rule**: area plus the rule above.

Column notes:
- "to land" counts the cells of the small components that (a) removes as well as the bridge cells.
- Box columns give the land cells of today / area / rule, out of the box's cells.

| N | land cells, today / area / rule | land, million km² | today→area: to sea, to land | rule: to sea, to land | small components to land (cells) | isolated components | sea cells with mean z > 0, today / area / rule | Darién | Cuba | SF |
|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 735 / 726 / 731 | 146.30 / 144.48 / 145.49 | 10, 1 | 3, 8 | 2 (2) | 2 | 23 / 26 / 24 | 1/0/1 of 1 | 0/0/0 of 2 | – |
| 32 | 2939 / 2942 / 2953 | 146.39 / 146.54 / 147.08 | 14, 17 | 4, 15 | 4 (5) | 1 | 75 / 78 / 75 | 2/1/2 of 3 | 1/2/3 of 9 | – |
| 64 | 11849 / 11841 / 11861 | 147.56 / 147.46 / 147.71 | 71, 63 | 2, 22 | 6 (7) | 4 | 185 / 220 / 214 | 5/5/5 of 12 | 7/8/8 of 35 | – |
| 128 | 47258 / 47489 / 47518 | 147.13 / 147.84 / 147.94 | 237, 468 | 7, 36 | 11 (14) | 18 | 656 / 553 / 537 | 14/16/16 of 41 | 25/34/34 of 144 | 2/2/2 of 3 |
| 192 | 106056 / 107074 / 107126 | 146.74 / 148.15 / 148.23 | 489, 1507 | 5, 57 | 26 (30) | 25 | 1430 / 866 / 841 | 31/34/34 of 91 | 59/79/79 of 318 | 6/6/6 of 7 |

The 60″ land total is 148.36 million km². Today's rule falls 1.6 million km² short of it at N = 192, because a cell at exactly half its points counts as sea.

**Panama, N = 128.** Two of the three Darién sea cells become land:
- 8.86N 77.80W: half of the points today → 0.73 of the area;
- 8.45N 78.08W: half → 0.86.

The third, 9.33N 78.09W (0.33), is the San Blas coast and stays sea. **Panama stays closed at every N** in all three masks, and at N = 1024.

**Cuba, N = 128:** 25 → 34 land cells. **San Francisco:** no change.

## Controlling depths

The table compares the model's controlling depth with the fine one, in meters.
- **Model depth:** the widest model ocean path's sill between the sea cells nearest two points inside the basins. Paths are confined to a box around the passage, widened by half a cell spacing.
- **N = 1024 depth:** the depth at which the two basins' fine water first joins, from the deepest fine cell within 0.3° of each point.
- **Cell format:** today / area / rule, then in brackets the fine depth between the representatives of the same two model cells, within the same box.
- "1c" means both points fall in one model cell.
- Today's figures include the forced SEA_STRAITS sills.

| passage (N = 1024) | N = 16 | N = 32 | N = 64 | N = 128 | N = 192 |
|---|---|---|---|---|---|
| Gibraltar (341) | 0/0/341 [341] | 0/0/341 [341] | 300/341/341 [341] | 300/341/341 [341] | 300/341/341 [341] |
| Bab-el-Mandeb (137) | 0/0/137 [0] | 150/0/137 [137] | 150/0/137 [137] | 104/137/137 [137] | 99/137/137 [137] |
| Hormuz (83) | 1c | 90/86/86 [86] | 90/83/83 [83] | 69/83/83 [83] | 47/83/83 [83] |
| Bering (44) | 0/44/44 [44] | 44/44/44 [44] | 44/44/44 [44] | 44/44/44 [44] | 44/44/44 [44] |
| Lombok (216) | 1c | 1c | 336/216/216 [216] | 336/216/216 [216] | 298/216/216 [216] |
| Danish straits (18, Darss sill) | 1c | 1/24/24 [24] | 0/18/18 [18] | 0/0/18 [18] | 8/18/18 [18] |
| Bosporus (0) | 0/0/0 [0] | 1c/1c/0 [0] | 1983/0/1982 [1982] | 0/0/0 [0] | 0/0/0 [0] |
| Hudson Strait (163) | 0/0/195 [195] | 195/195/195 [195] | 176/176/176 [176] | 143/165/165 [165] | 148/163/163 [163] |
| Tsushima (122) | 109/122/122 [0] | 109/122/122 [0] | 103/122/122 [122] | 109/122/122 [122] | 109/122/122 [122] |
| Tatar Strait (5) | 943/1311/1311 [0] | 0/0/0 [5] | 0/0/0 [5] | 0/5/5 [5] | 0/0/0 [5] |
| Florida Straits (737) | 576/737/737 [737] | 827/737/737 [660] | 632/737/737 [737] | 641/737/737 [737] | 634/737/737 [737] |
| Denmark Strait (609) | 563/609/609 [609] | 563/609/609 [609] | 563/609/609 [609] | 563/609/609 [609] | 516/609/609 [609] |
| Faroe Bank Channel (832) | 922/832/832 [832] | 922/832/832 [832] | 899/832/832 [832] | 922/832/832 [832] | 922/832/832 [832] |
| Panama (0) | 0/1c/0 [0] | 0/0/0 [0] | 0/0/0 [0] | 0/0/0 [0] | 0/0/0 [0] |

At N ≥ 64 the rule's value matches the fine value between the same cells everywhere except Tatar. The N = 1024 sills agree with the literature: Hanish 137 m, Darss 18 m, Denmark Strait ~620 m, Faroe Bank Channel ~840 m, Bering ~50 m. Gibraltar reads 341 m against the Camarinal sill's ~290 m, because the 60″ raster is coarser than the sill.

What remains:
- **Coarse cells straddling a passage.** At N ≤ 64 one cell can hold both basins; its representative then sits on the ocean side. Hudson at N = 16–64 reads 176–195 between the representatives, and the rule matches that. The bracketed 0s at N = 16–32 and the 1311/1982 readings come from the box: at those N the nearest model sea cell to a point lies outside the passage. They are not readings of the passage.
- **Bosporus at N = 64.** The Sea of Marmara (one cell) turns to land under (a), and a Black Sea cell is the nearest to its point. At N ≥ 128 the Black Sea (152 / 330 cells) and Marmara are lakes: the Bosporus and Dardanelles, 0.7–1.2 km wide, are unresolved even at 60″. Today's mesh also keeps the Black Sea closed at N = 128 and 192.
- **Tatar at N = 32, 64 and 192.** The model reads 0 against the fine 5 m. The Japan Sea and the Okhotsk are already one component through La Pérouse, so the rule cuts no second passage. Only joining separate components is in the rule.

## What (a) removes

These are the small ocean components that turn to land, and the strait flips they made under the earlier version of the rule.

**Small components to land:**
- N = 16: the Baltic and the Persian Gulf, one cell each.
- N = 32: the White Sea, the Sea of Marmara, a Gulf of Boothia pocket, the Gulf of California head.
- N = 64: Bothnian Bay, the Gulf of Finland, Marmara, a Volga pocket by the Caspian, two Canadian Arctic inlets.
- N = 128 (11, 14 cells):
  - Alaska: Cook Inlet, Chatham Strait.
  - The Canadian Arctic: Eureka Sound, Prince Regent Inlet, Roes Welcome, Frobisher Bay, Kennedy Channel.
  - Elsewhere: Chesapeake Bay, a Scoresby Sund pocket, Río de la Plata, Lagoa dos Patos.
- N = 192 (26, 30 cells): Svalbard's Isfjorden; the Gulf of Corinth; the Volga delta; the Gulf of Aqaba; the Mangyshlak gulf; the Ob and Yenisei pockets; the Seto Inland Sea; Kolyuchin Bay; Alaska's inside passage (3); Eureka Sound; Queen Charlotte Strait; Long Island Sound; the St. Lawrence estuary (2); Greenland fjords (2); the Amazon mouth; the Banc d'Arguin; the Gulf of Khambhat; Halmahera; Shark Bay; Chiloé; Lagoa dos Patos.

**Strait flips no longer made:**
- N = 128: Cook Inlet; Chatham Strait (2); the Canadian Arctic channels (4); Hudson and Frobisher (3, one cell at 0.91 land); Río de la Plata.
- N = 192: 44, including Corinth, Aqaba, the Volga, the Ob and Yenisei, Tatar, Kanmon and Seto, Alaska's inside passage, the St. Lawrence, the Amazon, Khambhat, Halmahera, Chiloé.

**Strait flips still made:**
- N = 128: the Great Belt, the White Sea throat, the Ob and Taz gulfs (3), Khatanga, Juan de Fuca.
- N = 192: the Mangyshlak gulf of the Caspian (2), Khatanga, Juan de Fuca, an inlet of north Baffin Island.

## Lakes and depressions

These are the components the fine water joins to nothing. Today's geography counts every point with surface ≤ 0 m as sea. A closed basin whose cells are mostly below sea level therefore becomes an isolated ocean basin in both today's mask and the new one, and the rule leaves it there.

That covers:
- the Caspian (218 cells at N = 128), the Black Sea and Marmara (see above);
- dry depressions below sea level: Qattara (6 cells at N = 128), the Chott Melrhir, Turpan, Danakil, the Dead Sea, the Salton Sink, Lake Eyre;
- lakes whose ETOPO 2022 bed lies below sea level: Superior and Michigan (6), Ontario (3), Baikal (4 + 1), Lake Melville;
- Lake Maracaibo (5);
- Patagonian sounds.

At N = 128, 18 such components; at N = 192, 25. Nearly all of their cells are sea in today's geography too: 405 of 414 at N = 128.

## Open points

1. **Porous barriers.** Adcroft (2013) open widths wait. The fine edges give the data for them.
2. **Phase 2 is not started, by Joshua's choice.** Not done: ice and cavity pass-through; sub-grid files and sills for the model; an option. Nothing in the model reads this code.
3. **Lakes and depressions** stay sea, as today. Turning them to land, or giving lakes their own surface level, is a separate decision.
4. **Passages between joined basins** (Tatar) are not cut. Cutting them would mean a second pass that matches every pair of basins, not only separate components.
