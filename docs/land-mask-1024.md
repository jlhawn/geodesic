# Land mask from an N = 1024 resampling: phases 1 and 2

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
- Per edge: the sill, the depth at which the deepest 60″ points of the two cells join through their own points, four-connected. A diagonal step slips past a crest one point wide (see Strait depths against the literature).

Every fine cell holds at least one point; no point lacks data. The total land area is 148.36 million km².

## The rule (`coarseGeography` + `connectedGeography`)

1. **Class.** A cell is land where its area-weighted land fraction exceeds 0.5.

2. **Edge depth.** Exchange through a strait is controlled at its sill (Armi & Farmer 1986, J. Fluid Mech. 164; Bryden & Kinder 1991, Deep-Sea Res. 38). Adcroft (2013, Ocean Modelling 67, 13–27) derives model face depths from high-resolution topography the same way. Here:
   - Each model cell's water is represented by one fine cell: the one whose water stays joined to the world ocean deepest (`oceanDepth`, a widest path from the deepest fine cell). Lakes have no such join, so for them the deepest fine cell is used.
   - A model edge's `edgeSill` is the depth to which the representatives of its two cells stay joined through fine edges.
     - The search covers the fine cells of the two cells and of the two cells at the edge's ends, so a channel that leaves through a third cell counts.
  - It widens to every neighbor of the two cells where the edge comes out shallower than both representatives' depth to the ocean, and the cells are joined at all or fine water crosses straight between them. Without the direct-crossing case, the Danish straits at N = 64 were cut: the water crosses the Darss sill straight between two cells, but the small region misses the representatives' path.
     - Where land parts them at the fine scale, the edge is 0.
   - Widening can only deepen an edge.
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
   - A cell that a strait needs stays sea, and so does any cell whose flip would split an ocean component. There are no conflicts at any N.

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
| 32 | 2939 / 2942 / 2950 | 146.39 / 146.54 / 146.93 | 14, 17 | 5, 13 | 3 (3) | 1 | 75 / 78 / 75 | 2/1/2 of 3 | 1/2/3 of 9 | – |
| 64 | 11849 / 11841 / 11860 | 147.56 / 147.46 / 147.70 | 71, 63 | 3, 22 | 6 (7) | 4 | 185 / 220 / 214 | 5/5/5 of 12 | 7/8/8 of 35 | – |
| 128 | 47258 / 47489 / 47516 | 147.13 / 147.84 / 147.93 | 237, 468 | 7, 34 | 9 (12) | 18 | 656 / 553 / 538 | 14/16/16 of 41 | 25/34/34 of 144 | 2/2/2 of 3 |
| 192 | 106056 / 107074 / 107124 | 146.74 / 148.15 / 148.22 | 489, 1507 | 5, 55 | 24 (28) | 25 | 1430 / 866 / 841 | 31/34/34 of 91 | 59/79/79 of 318 | 6/6/6 of 7 |

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
| Gibraltar (283) | 0/0/283 [283] | 0/0/283 [283] | 300/283/283 [283] | 300/283/283 [283] | 300/283/283 [283] |
| Bab-el-Mandeb (134) | 0/0/134 [0] | 150/0/134 [134] | 150/0/134 [134] | 104/134/134 [134] | 99/134/134 [134] |
| Hormuz (83) | 1c | 90/85/85 [85] | 90/83/83 [83] | 69/83/83 [83] | 47/83/83 [83] |
| Bering (44) | 0/44/44 [44] | 44/44/44 [44] | 44/44/44 [44] | 44/44/44 [44] | 44/44/44 [44] |
| Lombok (215) | 1c | 1c | 336/215/215 [215] | 336/215/215 [215] | 298/215/215 [215] |
| Danish straits (17) | 1c | 1/24/24 [24] | 0/17/17 [17] | 0/0/17 [17] | 8/17/17 [17] |
| Bosporus (0) | 0/0/0 [0] | 1c/1c/0 [0] | 1983/0/1982 [1982] | 0/0/0 [0] | 0/0/0 [0] |
| Hudson Strait (163) | 0/0/195 [195] | 195/195/195 [195] | 176/176/176 [176] | 143/165/165 [165] | 148/163/163 [163] |
| Tsushima (116) | 109/116/116 [0] | 109/116/116 [0] | 103/116/116 [116] | 109/116/116 [116] | 109/116/116 [116] |
| Tatar Strait (4) | 943/1310/1310 [0] | 0/0/0 [4] | 0/0/0 [4] | 0/4/4 [4] | 0/0/0 [4] |
| Florida Straits (735) | 576/735/735 [735] | 827/735/735 [660] | 632/735/735 [735] | 641/735/735 [735] | 634/735/735 [735] |
| Denmark Strait (609) | 563/609/609 [609] | 563/609/609 [609] | 563/609/609 [609] | 563/609/609 [609] | 516/609/609 [609] |
| Faroe Bank Channel (827) | 922/827/827 [827] | 922/827/827 [827] | 899/827/827 [827] | 922/827/827 [827] | 922/827/827 [827] |
| Panama (0) | 0/1c/0 [0] | 0/0/0 [0] | 0/0/0 [0] | 0/0/0 [0] | 0/0/0 [0] |

At N ≥ 64 the rule's value matches the fine value between the same cells everywhere except Tatar.

## Strait depths against the literature

`scripts/straitSills.mjs` reads each passage's controlling depth straight from the 60″ surface: the widest path between the basins' deepest nearby points through the points of the passage's box. It steps either to a point's eight neighbors or to its four.

| passage | 60″ points, 8-connected | 60″ points, 4-connected | N = 1024 mesh, 4-connected | published | source |
|---|---|---|---|---|---|
| Gibraltar | 341 | 283 | 283 | ~284 (Camarinal Sill) | Armi & Farmer 1988 |
| Bab-el-Mandeb | 137 | 134 | 134 | ~137 (Hanish Sill) | Werner & Lange 1975 |
| Hormuz | 83 | 83 | 83 | ~70–100 | approximate, source not checked |
| Bering | 44 | 44 | 44 | ~50 | Woodgate et al. 2005 |
| Lombok | 216 | 215 | 215 | ~300 | Murray & Arief 1988; Sprintall et al. 2009 |
| Danish straits | 18 | 18 | 17 | 18 (Darss Sill) | approximate, source not checked |
| Bosporus | 0 | 0 | 0 | ~35, 0.7 km wide | Gregg & Özsoy 2002 |
| Hudson Strait | 164 | 164 | 163 | – | not checked |
| Tsushima | 122 | 118 | 116 | ~130 | approximate, source not checked |
| Tatar Strait | 5 | 5 | 4 | ~7–10 | approximate, source not checked |
| Florida Straits | 737 | 735 | 735 | ~750–800 | approximate, source not checked |
| Denmark Strait | 609 | 609 | 609 | ~620 | approximate, source not checked |
| Faroe Bank Channel | 832 | 829 | 827 | ~840 | Hansen & Østerhus 2000 |

- **Gibraltar.** The eight-connected step slips diagonally past the Camarinal crest. Four-connected reads 283 m, the published sill. The fine edges are four-connected.
- **Lombok.** 215–216 m in either connectivity, so the shallow reading is in ETOPO 2022's 60″ data. It is left as it is, under the ruling against hand-made depths.
- **Bosporus.** Closed at 60″ in either connectivity: the strait is narrower than the grid.
- **Other passages.** Within 0–6 m between the two connectivities. Open or closed changes nowhere.

The N = 1024 column comes from each fine edge's sill, the join of the deepest points of two ~7 km cells. It can read a few meters shallower than the point-level path (Tsushima 116 against 118).

What remains:
- **Coarse cells straddling a passage.** At N ≤ 64 one cell can hold both basins; its representative then sits on the ocean side. Hudson at N = 16–64 reads 176–195 between the representatives, and the rule matches that. The bracketed 0s at N = 16–32 and the 1311/1982 readings come from the box: at those N the nearest model sea cell to a point lies outside the passage. They are not readings of the passage.
- **Bosporus at N = 64.** The Sea of Marmara (one cell) turns to land under (a), and a Black Sea cell is the nearest to its point. At N ≥ 128 the Black Sea (152 / 330 cells) and Marmara are lakes: the Bosporus and Dardanelles, 0.7–1.2 km wide, are unresolved even at 60″. Today's mesh also keeps the Black Sea closed at N = 128 and 192.
- **Tatar at N = 32, 64 and 192.** The model reads 0 against the fine 5 m. The Japan Sea and the Okhotsk are already one component through La Pérouse, so the rule cuts no second passage. Only joining separate components is in the rule.

## What (a) removes

These are the small ocean components that turn to land, and the strait flips they made under the earlier version of the rule.

**Small components to land:**
- N = 16: the Baltic and the Persian Gulf, one cell each.
- N = 32: the White Sea, the Sea of Marmara, the Gulf of California head.
- N = 64: Bothnian Bay, the Gulf of Finland, Marmara, a Volga pocket by the Caspian, two Canadian Arctic inlets.
- N = 128 (9, 12 cells):
  - Alaska: Cook Inlet, Chatham Strait.
  - The Canadian Arctic: Eureka Sound, Prince Regent Inlet, Roes Welcome, Frobisher Bay, Kennedy Channel.
  - Elsewhere: Río de la Plata, Lagoa dos Patos.
- N = 192 (24, 28 cells): Svalbard's Isfjorden; the Gulf of Corinth; the Volga delta; the Gulf of Aqaba; the Mangyshlak gulf; the Ob and Yenisei pockets; the Seto Inland Sea; Kolyuchin Bay; Alaska's inside passage (3); Eureka Sound; Queen Charlotte Strait; the St. Lawrence estuary (2); Greenland fjords (2); the Amazon mouth; the Banc d'Arguin; the Gulf of Khambhat; Halmahera; Shark Bay; Chiloé.

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

## Phase 2: the `n1024` geography option

Geography records gain a `mask` key: `'raster'` or `'n1024'`. A raster record has no key, so every existing record, state and digest is unchanged. An n1024 record:
- needs `topography: 'etopo2022'`, observed ice sheets, and ice shelves `'land'` or `'cavity'`;
- reads `data/geography_n1024_N<N>.bin` (`maskUrl`, `encodeMask`, `decodeMask`) and `data/subgrid_etopo2022_n1024_N<N>.bin`.

Build:
```
LAKES=~/.cache/webgcm-lakes/ne_10m_lakes.geojson LAKE_REPORT=lakes.json node --max-old-space-size=20000 scripts/fineGeography.mjs
LAKE_REPORT=lakes.json REPORT=mask.json node --max-old-space-size=16000 scripts/geographyMask.mjs
MASK=n1024 TOPOGRAPHY=etopo2022 node scripts/subgridTerrain.mjs CACHE 16 32 64 128 192
```
The fine build takes 5.4 min (18.5 GB peak); the five masks take 5.4 min (5.7 GB).

| N | land cells | land, 10⁶ km² | lakes | lake cells | depressions (cells) | seawater lakes (cells) | closed seas | file |
|---|---|---|---|---|---|---|---|---|
| 16 | 731 | 145.49 | 0 | 0 | 0 | 0 | Black Sea, Caspian | 35 KB |
| 32 | 2948 | 146.83 | 2 | 2 | 0 | 0 | Black Sea, Caspian | 137 KB |
| 64 | 11832 | 147.35 | 14 | 33 | 2 (3) | 1 (1) | Black Sea, Caspian | 544 KB |
| 128 | 47359 | 147.44 | 47 | 196 | 6 (15) | 3 (7) | + Marmara, a Magellan fjord | 2.1 MB |
| 192 | 106670 | 147.59 | 106 | 556 | 11 (39) | 6 (17) | + Marmara | 4.8 MB |

**Ice.** The fine dataset also carries a sea mode, in which floating ice is water to its bed. Each mask file stores the shelf mode in full and the cells and edges where the sea mode differs. `'cavity'` builds on the sea mode.

**Lakes.** The lakes come from Natural Earth `ne_10m_lakes` (1355 features), rasterized at each 60″ point.
- A coarse cell is a lake where lakes cover more than half of its area. It takes the id of the lake that covers most of it.
- Level: the mode of ETOPO's surface over the lake when at least 30 % of its points share it (the raster holds the water); otherwise the 10th percentile of the shore ring.
- Depth: level − bed, averaged over the cell.
- A lake cell is water at elevation −depth with its own surface (`surface`, used by `surfaceGeopotential`). Its edges open only to cells of the same lake.
- Natural Earth's "Alkaline Lake" class marks a lake as saline.

Natural Earth lists neither the Caspian nor Maracaibo. The dataset has no intermittent class: salt pans such as Eyre, Frome and Torrens are "Alkaline Lake", as are the Dead Sea and the Great Salt Lake. They first take cells at N = 128 (10 lakes).

**Depressions and cut-off seawater.** Every sea component but the largest is a closed basin.
- A closed basin deeper than 500 m below sea level at its deepest point stays a closed sea: the Black Sea, the Caspian, Marmara.
- A shallower one the lakes don't mark is seawater the mesh cuts off if its fine water (60″ points below 0 m outside lakes) joins the world ocean's through neighboring fine cells that both hold some. It becomes saline lake cells at level 0 (`seawater`).
- Otherwise it is a dry depression and becomes land at its own elevation (`depression`).

The dry depressions:
- N = 64: Qattara, Eyre.
- N = 128 adds: the Chott Melrhir, Danakil, Turpan, Salton.
- N = 192 adds: the Dead Sea, a Libyan sabkha (30.2° N 21.1° E), more of Eyre, and one Skyring cell.

The seawater lakes:
- N = 64: Maracaibo.
- N = 128 adds: Lake Melville, Chesterfield Inlet.
- N = 192 adds: Lagoa dos Patos, Skyring and Otway, and the Sacramento delta.

The delta is wrong. ETOPO counts every point at or below 0 m as water, so its subsided islands (to −7 m) are water. The same holds for the polders and for the Caspian lowland between −28 and 0 m: the model's Caspian is 605,000–751,000 km² against about 371,000. An elevation-independent land–water source would fix these.

**Lake model.** `lake: { model: 'flake' }` (off by default; `js/physics/lake.module.js`) runs FLake (Mironov 2008; Mironov et al. 2010) on lake cells. It follows the IFS's choices (Balsamo et al. 2012; Choulga et al. 2019):
- depth bounded to 2–50 m;
- no snow layer and no bottom sediments;
- one shortwave band at 3 m⁻¹, with the ice opaque.

Where ETOPO holds no bed below a lake's level (Victoria, Tanganyika, Malawi, Great Bear), the depth is GLDB v3's default: 10 m fresh, 5 m saline.

Coupling:
- The atmosphere sees open sea water; the ice has FLake's albedo.
- Snow's heat of fusion is taken from the lake.
- FLake's own ice does not conserve heat: the residual is 3.7 % of the flux under the ice in the 20 m column test, −87 MJ/m² over three years. The model takes the ice thickness from the energy instead (`lakeSettle`): FLake's equations set the water and the ice's temperature, and the ice is as thick as the remaining heat requires.
  - The 20 m column then drifts 4e-4 J/m² in three years.
  - Against FLake's own ice: peak 1.04 vs 1.07 m, break-up day 122 vs 123, freeze-up day 345 for both.
- The shortwave that reaches the bottom stays in the water (no sediments).
- Snowfall's heat of fusion is taken from the lake at once, through the same closure.
- A lake passes its whole water gain (P − E plus inflowing runoff) to the sea cell its outlet reaches, negative where evaporation wins. Its depth stays fixed, and global water is conserved.

## Open points

1. **Porous barriers.** Adcroft (2013) open widths wait. The fine edges give the data for them.
2. **Salt pans** count as lakes where Natural Earth's polygons cover most of a cell: Eyre, Frome, Torrens and Mackay at N = 128. Natural Earth's playas file could separate them.
3. **Passages between joined basins** (Tatar) are not cut. Cutting them would mean a second pass that matches every pair of basins, not only separate components.
