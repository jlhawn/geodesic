# The MACv2 aerosol climatology: `radiation.aerosol: 'mac'`

Oct 27, 2026. Branch `claude/aerosol-mac`. `radiation: { aerosol: 'mac' }`, the default on both engines (CPU f64, WGSL f32), replaces the flat visible aerosol of the shortwave (0.12 over land, 0.07 over sea, one albedo, one asymmetry, one 2 km profile) with the MACv2 climatology of Kinne (2019): for each month and cell the 550 nm optical depth, single scattering albedo and asymmetry of a fine and a coarse mode and the mean height of each. A table of 120 constants gives each mode's optics at the twelve solar wavelengths of the data, which are the centers of the shortwave bands of RRTMG_SW, so the aerosol enters the near infrared for the first time. `aerosol: 'flat'` is the previous default: asked for it, a model's CPU results, WGSL text, PH layout and saved states are bit-identical to what they were before the climatology became the default.

## Sources

- Kinne, S. (2019). The MACv2 aerosol climatology. Tellus B 71, 1623639, https://doi.org/10.1080/16000889.2019.1623639. The climatology: how the fine and coarse modes are built from sun- and sky-photometer statistics with ECHAM-HAM as context, their optical properties, the reference year 2005.
- Kinne, S. (2019). Aerosol radiative effects with MACv2. Atmospheric Chemistry and Physics 19, 10919 to 10959. The effects the climatology has on the radiation, which the model's change is read against: for the year, all sky, the total aerosol changes the absorbed sunlight by -1.8 W/m² at the top, +3.7 in the atmosphere and -5.5 at the surface.

## Options

| option | default | meaning |
|---|---|---|
| `aerosol` | `'mac'` | `'mac'` the climatology (`DEFAULT_AEROSOL` in `js/aerosol.module.js`); `'flat'` the previous default, the visible aerosol of `landAerosol`, `seaAerosol`, `aerosolAlbedo`, `aerosolAsymmetry` and `aerosolHeight`, which play no part under `'mac'` |
| `aerosolScale` | `1` | factor on both modes' depths, 0 or more. 0 is the model without aerosol, bit for bit |
| `aerosolData` | the file | the climatology's bytes, or the decoded climatology (`decodeAerosol`). Under node it defaults to `data/aerosol_mac_2deg.bin`, and a missing file is an error that names the builder. The page's worker fetches the file and passes the bytes |

`'mac'` is refused, with the same words by the CPU radiation and the GPU model, with `solarGases` other than `'clirad'` and with `vaporShielding: 'cloudTop'`; the message says to give `aerosol: 'flat'`, as the layered soil's says to give `soil: 'slab'` for a `heatCapacity`. It runs with both shortwave schemes (`'layered'`, the default, and `'sequential'`), with either near infrared (`'bands'`, the default, and `'clirad'`), with held radiation (`radiationEvery`) and with the clear-sky pass. A saved state's header records `radiation.aerosol` with the other options the run was given and not the data, and `scripts/spinup.mjs` and `scripts/splitState.mjs` stamp the default (`recordedOptions` in `scripts/stateOptions.mjs`) where the run gave none; the page's default states record no aerosol option, so they start with the climatology.

The page's worker fetches `data/aerosol_mac_2deg.bin` (1.56 MB) whenever the run's resolved option is `'mac'` (`climatologyOn`): a start that names no aerosol, a state that records none and a speed probe included; a start from a state that records `'flat'`, or `options: { radiation: { aerosol: 'flat' } }` in the worker's start message, fetches nothing. The scripts read the file under node, and take `RADIATION='{"aerosol":"flat"}'` for the previous default. A model built without physics (`physics: false`, the worker's regridding source) has no sampler, and `js/charts.module.js`, which builds a model for its core alone, gives `aerosol: 'flat'`.

## The default

The project owner decided on October 27, 2026 that the climatology would be the default if the results were good on net. The evidence is the fixed-SST A/B of the bl48 control, the climatology with the layered soil against the same without it, three seeds over model days 882 to 892 (August 17 to 27), which ran on the earlier file, the one whose heights were August's for every month.

For:

- Net top of atmosphere +0.06 W/m². The effective forcing is near zero, so a coupled global mean should hardly move.
- 10S-10N clear atmosphere 94.2 to 99.0 W/m² against CERES 100.0; clear surface 276.5 to 268.5 against 265.3.
- Equatorial stress rms (140E to 100W) 27.0 to 23.6 in every seed; the excess at 180 to 160W fell 16 percent.
- SWCRE -51.2 to -48.4 W/m² against -43.6.
- India rain 0.69 to 0.77 times GPCP; the Indian Ocean's 2.02 to 1.71 times.

Against:

- Open warm-pool rain 1.00 to 0.92 times GPCP in every seed.
- Daily temperature range lower by 0.1 to 0.5 K outside the deserts.
- The global clear atmosphere overshoots by about 1.7 W/m².

It was judged good on net. For the days of the A/B the full-year file differs from the earlier one only in referring each mode's heights to the dataset's surface altitude, which lowers the mean heights by 0.1 to 0.3 km (the table of heights below), and the A/B was not repeated on it.

In a run the default adds a fetch of 1.56 MB (the page's worker, once), eight floats a cell of memory (1.3 MB at N = 64, 5.2 MB at N = 128, on the host and on the device), 2.5 percent of the GPU's pace at N = 128 and 8 percent of a sequential CPU column (below), and moves the sunlight the model absorbs as in the table under Measured. The page's default states and every state that records no aerosol option start with the climatology. Anything that needs `solarGases` other than `'clirad'` or `vaporShielding: 'cloudTop'`, and every digest, figure and bound made under the previous default, gives `aerosol: 'flat'` (`PREVIOUS_AEROSOL` in `test/helpers/previousDefaults.mjs`, which `PREVIOUS_RADIATION` includes; `PUBLISHED_RADIATION` in `explain/published.module.js` for the explainers, whose configuration and figures are unchanged).

## The data

The files are the Max Planck Institute for Meteorology's, on its anonymous FTP server, ftp://ftp-projects.mpimet.mpg.de/aerocom/climatology/MACv2_2018/ (the README is `for_ICON/README_MACv2`). The builder reads four of them, 145 MB:

| file under `for_ICON/` | size | content |
|---|---|---|
| `spec_2005/by_wave/gt_f_00550nm.nc` | 10.1 MB | fine mode at 550 nm: `aod`, `ssa`, `asy` for the 12 months on a 1 degree grid |
| `spec_2005/by_wave/gt_c_00550nm.nc` | 10.1 MB | the coarse mode |
| `altitude/aeropt_fin_alt_km20.nc` | 62.5 MB | the fine mode's share of the column depth in 20 layers above sea level, by month, and `asl`, the surface altitude |
| `altitude/aeropt_coa_alt_km20.nc` | 62.5 MB | the coarse mode's |

### License and attribution

The server, the README and the two papers state no terms of use for these files. The derived file, `data/aerosol_mac_2deg.bin`, is in the repository and the page serves it, with the attribution in `docs/aerosol-data-notice.md`: the source, the files that were read and when, what the derivation did to them (the 2 degree resampling, the eight one-byte fields, the heights above the dataset's surface altitude) and the license status, which is that none is stated.

### Building the file

The file is committed; to make it again from the sources:

    BASE=ftp://ftp-projects.mpimet.mpg.de/aerocom/climatology/MACv2_2018/for_ICON
    mkdir -p RAW_DIR && cd RAW_DIR
    curl -O $BASE/spec_2005/by_wave/gt_f_00550nm.nc
    curl -O $BASE/spec_2005/by_wave/gt_c_00550nm.nc
    curl -O $BASE/altitude/aeropt_fin_alt_km20.nc
    curl -O $BASE/altitude/aeropt_coa_alt_km20.nc
    python3 -I scripts/packAerosol.py RAW_DIR data/aerosol_mac_2deg.bin

`scripts/packAerosol.py` uses the standard library only and takes a third argument, FACTOR (2 by default), the number of source cells per side that make one cell of the output. It coarsens the 1 degree data to 2 degrees: the depth by area, 1 - SSA by the depth, the asymmetry by depth times SSA and the height by the depth. An altitude file may instead be a one-month record, `aeropt_fin_alt_km20.record` or `aeropt_coa_alt_km20.record`, cut from the file by byte range (a big-endian float32 `time` and the 20 layers); every month then has that month's heights, and the ground of a cell is the lower edge of its lowest layer that holds aerosol, there being no surface altitude in a record. A complete file gives each month its own heights over its `asl`, as the committed file has.

The file is little-endian, in the layout of `data/woa_sst_monthly_1deg.bin`:

| bytes | content |
|---|---|
| 0 | `'AERM'` |
| 4 | uint32 byte offset of the data |
| 8, 10, 12, 14 | uint16 nLon, nLat, nMonth, nField (8) |
| 16 | float32 lon0, dLon, lat0, dLat (degrees; the center of column 0 and row 0) |
| 32 | float32 center[nMonth], the day of a 365-day year from January 1 |
|  | per field float32 low, float32 high, uint32 scale (0 linear, 1 square root) |
|  | uint16 n, then n bytes of ASCII provenance, zero padding to a multiple of 8 |
| data | uint8 [field][month][lat][lon], rows south to north, columns west to east from lon0 |

The eight fields, one byte each, as low + (high - low) s with s = q / 255 (squared for the depths): the 550 nm depth of the fine and the coarse mode (0 to 1.2, a relative step of 2.7 percent at 0.1), their 1 - SSA (0 to 0.3), their asymmetry (0.45 to 0.85) and the mean height of each mode's aerosol above the ground (0 to 8 km). At 2 degrees, 12 months and 8 fields the file is 1,555,624 bytes (0.7 MB deflated).

### The spectral constants

`js/physics/aerosolSpectra.module.js` is generated by `scripts/aerosolFit.mjs` from the band files `spec_2005/by_wave/gt_{f,c}_WWWWWnm.nc` of the twelve wavelengths (or the August record of each, cut by byte range, 0.78 MB):

    WRITE=1 node scripts/aerosolFit.mjs DIR [DIR...]

For each mode and wavelength b it fits, by least squares over the cells of August weighted by each cell's 550 nm depth and area, the five constants that take the mode's 550 nm triple (τ, x = 1 - ω, g) to the wavelength's:

    τ_b = r τ,   x_b = a + c x,   g_b = d + e g

(r is the weighted mean ratio of the depths), 120 constants in all. Beyond 2.5 micrometer (bands 17, 16 and 29, 3 percent of the beam) the albedo and asymmetry do not follow the 550 nm values and the band means stand in (c = e = 0). The script also writes the Planck weights (5778 K) of each wavelength's range in the sequential scheme's two streams and the layered scheme's visible slices, and prints the fit's misses against the band data in six regions. Without `WRITE=1` it only prints.

| wavelength (nm) | 400 | 550 | 700 | 1000 | 1270 | 1460 | 1780 | 2050 | 2320 | 2790 | 3470 | 8000 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| RRTMG_SW band | 26 | 25 | 24 | 23 | 22 | 21 | 20 | 19 | 18 | 17 | 16 | 29 |

The ultraviolet bands of the data (230 and 300 nm) are not used: the ultraviolet is ozone's.

## In the model

### Time

The file holds monthly means at the data's mid-month days. The model day d (the days since the March equinox, which is day 78 of the year) maps to the day of the year (d + 0.5 + 78) mod 365; the cell's values are bilinear in the grid and linear between the two months whose centers bracket the day (`monthWeights` of `js/ocean/climatology.module.js`), in the stored variables (depth, 1 - SSA, asymmetry, height) before they become optics. The year is the 2005 of the data and does not change. `createAerosolSampler(mesh, data)` of `js/aerosol.module.js` holds the eight fields of a mesh's cells ([field][cell], Float32) and refreshes them when the model day changes; `createModel` calls it from the wrapper around `radiation.setTime`.

### Each cell's optics

From a cell's two 550 nm triples `streamOptics` gives each mode's depth, 1 - ω and asymmetry in the twelve wavelengths by the constants above, held within [0, 1] and [0, 0.99], and the sequential scheme's mixtures: the Planck-weighted depth of the visible stream (0.297 to 0.7 µm) and of the rest of the beam, the depth-weighted single scattering albedo and the scattering-weighted asymmetry of each, and the fine mode's share of the absorption depth of each stream, which weighs the two profiles.

### The profile

One mean height μ above the ground (km, within [0.3, 6]) per mode and cell, placed on the model's own layers: constant density from the ground to the depth D = min(1.5 km fine, 1.0 km coarse; 1.8 μ), an exponential above it whose scale height H = ½((μ - D) + √(μ² + 2μD - D²)) gives the mean height μ exactly. The heights above the ground are those of the layers' mid-levels in the model's own geopotential (the ground's from the lowest layer's, less the thickness down to the surface of the column's virtual potential temperature), the layers' interfaces are midpoints of those heights, the lowest is the ground and the top is infinitely high, so a layer's share is a difference of the cumulative fraction and the shares sum to 1. Over high terrain the aerosol sits by its height above the ground. The well-mixed layer is what the data's own profiles look like (flat through the lowest kilometer or two, which an exponential from the ground cannot do); the cap at 1.8 μ keeps the scale height at 0.18 μ or more (at 2 μ the mixed layer alone has the mean height and the exponential has none).

### The sequential scheme

The visible slab takes the visible mixture's depth, albedo and asymmetry in place of the flat aerosol's. The second slab is new: the rest of the beam (the near infrared) meets the near-infrared mixture after the gases above the aerosol (the profile-weighted cumulative absorption of the near infrared gases, `nearInfraredTaken`), absorbs 1 - exp(-(1 - ω) τ m) of what is left there, adds the delta-scaled scattering depth (1 - g) ω τ to the rest stream's depth, absorbs 1 - exp(-(1 - ω) τ 5/3) of the light the surface and the clear air send back, and heats the layers by the modes' profile weighted by their absorption depths in the stream. The clear column and the clear-sky pass are this scheme's.

### The layered scheme

The module takes the aerosol of each slice as arrays by layer (`aerosolLayers`): for the two visible slices and each of RRTMG_SW's near infrared bands, the depth, the absorption depth, the scattering depth and the asymmetry-weighted scattering depth, each the sum over the modes of the mode's value times its share of the layer. For `'flat'` the arrays are the old expressions to the last bit, and the near infrared has none. The layered scheme adds the cloud's difference from its own clear region, so its aerosol is in both.

### The GPU

`js/gpu/aerosol.gpu.js` generates the same arithmetic as WGSL, only when the option is `'mac'`; without it no text of the aerosol is generated and the default shaders hash alike before and after. The eight fields are the slots `PH_AER + field * C + i` of the physics buffer, appended after the engine's own regions, written once a model day through a staging buffer and a copy at the step's place in the queue (a batch submits a copy still waiting before the staging buffer is written again). The 12 wavelengths' constants are `const` arrays. The kernels fill the two profile arrays `macF` and `macG` counting up, in the entry point, read them in place and never pass them by pointer, which avoids the Metal failures listed in `CLAUDE.md`. The `light` of a column becomes a `vec4` whose fourth number is the near-infrared slab's scattering depth. The cumulative absorption of the near infrared gases is kept in the stash of the radiation kernel for the slab's beam. Memory: 8 floats a cell for the fields (1.3 MB at N = 64, 5.2 MB at N = 128, on the host and on the device).

### Cold starts

`initializeState` builds the initial profile from the radiative-convective equilibrium of one column of the model's own radiation (`equilibriumProfile` in `js/physics/init.module.js`), a mean column of the sun at a quarter of the solar constant that describes no place. That column meets the flat aerosol of cell 0 whatever the model's option (`column` takes `aerosol: 'flat'` in its `sunlight`, which changes nothing under `'flat'`), so the cold start of the default model is the previous default's, bit for bit, and the same state on the CPU engine, on the GPU engine (whose own radiation object holds no climatology) and in a model without physics. A column that met the climatology would meet that of cell 0 on the model's day in the CPU model and nothing in the other two, and a CPU and a GPU model started from `initializeState` would part by up to 3 K in the top layer and 6 Pa at the surface before their first step.

## Tests

- `test/aerosolClimatology.test.mjs` (CPU): the file's encode and decode round trip within the byte's step; the month blend and the wrap from December to January; the sampler's bilinear taps at a hundred cells against a direct evaluation; the model day's mapping to the year and the daily refresh; the loading from bytes, from the decoded form and from the file, and the missing-file error; the spectral constants (twelve wavelengths, the 550 nm triple taken to itself, the weights of the streams and slices summing to 1); `netcdfClassic`; `scripts/packAerosol.py` on a synthetic file set made in the test (with the dataset's surface altitude, and from one-month records) and its refusals; and, when `data/aerosol_mac_2deg.bin` exists, the shipped file's global means, byte ranges and stream depths against a fixture of the band data (skipped, saying why, when the file is absent).
- `test/aerosolOptics.test.mjs` (CPU): a mode's optics at 550 nm and the limits of the others; the mixing rules of `streamOptics`; the profile shares (sum to 1, the mean height, the constant density to the mixed depth, thinning above it); the shares over 4.5 km of terrain; the option checks; the flat aerosol of the cold start's column, which gives a model with the climatology the flat model's column and profile bit for bit.
- `test/layeredConservation.test.mjs` (CPU): absorbed plus reflected sunlight is the beam, and the layers' shortwave heating the atmosphere's absorption, in columns of the standard atmospheres, clear and cloudy, in both schemes, over sea and land, at a high and a low sun; the aerosol's absorption is positive and its heating lies where its mean height says; no depth and `aerosolScale: 0` are the model without aerosol bit for bit and the absorption scales with the depth; the closure of the clear region (below); and the N = 16 model's steps, held or not, over an aquaplanet and the fresh geography, across a change of the model day.
- `test/aerosolGpu.test.mjs` (GPU file): no WGSL of the aerosol and the engine's own PH layout and stash without `'mac'`, and the added region and array with it; the GPU's refusals are the CPU's; a cold start is the same state on both engines and with or without physics; the climatology of the model day lands in PH as the CPU sampler has it and the N = 16 model matches the CPU's, layered, sequential, held and with the CLIRAD near infrared, as it does without the aerosol; twelve steps across a change of day and of the month's weights; the change `'mac'` makes to the heating of each layer of the sunlit columns agrees between the engines; and over a 4 km mountain the engines agree on the layers more than 5 km above the ground.
- `test/parallel.test.mjs` (CPU): the worker-thread step under `'mac'` reproduces the single-thread step bit for bit across a change of the model day.
- Unchanged and passing: the shader-hash pins of `test/landSoilGpu.test.mjs`, `test/boundaryLayer.test.mjs`, `test/gpuModel.test.mjs` and every other test file that passed before the option.
- `scripts/compareEngines.mjs` maps the PH region to `radiation.aerosolFields`.

## Measured

### Explicit `'flat'` is the previous default

Against 10aa5d7, the last commit before the climatology became the default:

- The CPU model's arrays after 24 steps at N = 16 in ten configurations (plain, held radiation, sequential, clear-sky pass, land, no flat aerosol, no mixed-layer deck, the layered scheme held with its clear region, and the page's options with the fresh geography and the default ocean; the eleventh, the shielded one, is refused without `'flat'`) have the same SHA-256 under `aerosol: 'flat'`. A model day (64 steps of 1350 s) from the cold start at N = 16, of an aquaplanet and of the fresh geography with the default ocean, with and without the page's options, leaves the state, the ocean and the land bit-identical on the CPU and on the GPU.
- The 661 WGSL modules of sixteen GPU configurations at N = 16 (the ten above, the gray gases, the shielded sequential scheme, the CLIRAD near infrared and the default model with the default ocean, with and without the page's options) hash alike under `aerosol: 'flat'`; the PH layout and the stash are the previous ones.
- On an August N = 128 state of the KLM control (163,842 columns, eight instants of the day, the layered shortwave with the saved clouds and the clear-sky pass), 15 day-mean radiation fields under `aerosol: 'flat'` are identical to the bit to those of the model before the option existed.

The default's WGSL differs from the previous one's by the climatology and nothing else. The 85 shaders of the default model (N = 16, the fresh geography, the default ocean) hash `3fbb9c19778b60e1`, against `b41944bfabf0ff2d` with `aerosol: 'flat'` and at 10aa5d7; with the page's options, `2588fd1219dbd54f` against `265d9092be2dcae5`. Each hash is the SHA-256, to 16 digits, of the 85 module hashes in the order the model compiles them.

### What the climatology does to the sunlight

An August N = 128 state of the KLM control (model day 880, 17 August), the model's CPU shortwave with the saved clouds, deck and cumulus, day means of eight instants, `'mac'` minus `'flat'`, W/m² (negative is less absorbed sunlight):

| region | clear top | clear atmosphere | clear surface | all-sky top | all-sky atmosphere | all-sky surface |
|---|---|---|---|---|---|---|
| global | -2.31 | +3.44 | -5.75 | -0.54 | +3.30 | -3.84 |
| 10S-10N sea | -3.40 | +3.26 | -6.66 | -1.50 | +3.12 | -4.62 |
| 40-60N land | -2.01 | +4.35 | -6.37 | -0.42 | +4.03 | -4.46 |
| India | -4.64 | +9.91 | -14.55 | +0.66 | +7.41 | -6.76 |
| Sahara and Arabia | -1.67 | +21.02 | -22.69 | +2.28 | +20.03 | -17.75 |
| East Asia | -7.61 | +10.84 | -18.45 | -1.60 | +9.80 | -11.40 |

The aerosol's own absorption rises by 3.48 W/m² globally. The design study's prototype of the same configuration, run on the same state, gives -2.31, +3.45 and -5.76 W/m² for the global clear sky and -0.54, +3.32 and -3.85 for the all sky; each entry of the table agrees with its prototype's to 0.4 W/m² (the largest miss is the clear atmosphere of the Sahara and Arabia, 21.02 against 21.40) and most to 0.05. Kinne (2019) gives -1.8, +3.7 and -5.5 W/m² for the year, all sky; the model's annual figures were -0.98, +3.78 and -4.76 in the study, with the climatology's August larger than its year by about a fifth.

### The heights

Each month has its own mean heights, above the ground, which the earlier file did not: it had August's for every month. AOD-weighted regional means of the file, km (the range the earlier file had in every month):

| coarse mode | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec | earlier |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Sahara and Arabia | 1.0 | 1.1 | 1.4 | 1.7 | 1.9 | 2.1 | 2.2 | 2.3 | 2.1 | 1.7 | 1.2 | 0.9 | 2.5 to 2.6 |
| East Asia | 1.9 | 2.1 | 2.0 | 1.8 | 1.9 | 2.0 | 2.2 | 2.3 | 1.9 | 1.5 | 1.2 | 1.5 | 2.4 to 2.5 |
| global | 1.2 | 1.3 | 1.4 | 1.5 | 1.6 | 1.6 | 1.7 | 1.6 | 1.5 | 1.3 | 1.2 | 1.2 | 1.5 to 1.7 |

| fine mode | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec | earlier |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Sahel | 2.2 | 2.2 | 2.4 | 2.6 | 2.7 | 2.5 | 2.0 | 1.9 | 2.1 | 2.2 | 2.2 | 2.2 | 1.9 to 2.2 |
| southern Africa | 1.8 | 1.8 | 2.2 | 2.5 | 2.3 | 1.8 | 1.7 | 1.9 | 2.3 | 2.4 | 2.4 | 2.2 | 2.1 |
| Amazon | 2.2 | 2.1 | 2.0 | 2.1 | 2.1 | 1.8 | 1.9 | 1.9 | 2.1 | 2.5 | 2.6 | 2.3 | 2.1 to 2.2 |
| India | 1.3 | 1.5 | 1.8 | 2.2 | 2.3 | 2.4 | 2.0 | 1.8 | 1.9 | 1.6 | 1.4 | 1.3 | 1.9 to 2.0 |
| global | 1.8 | 1.8 | 1.9 | 2.0 | 2.0 | 2.0 | 1.9 | 1.9 | 1.9 | 1.9 | 1.9 | 1.8 | 2.0 to 2.1 |

Saharan dust lies 1.5 km lower in winter than the earlier file had it, 0.9 to 1.2 km above the ground from November to February where it was 2.5, and rises through the spring to 2.3 km in August, 0.2 to 0.3 km below the earlier value: over land the earlier heights counted from the lower edge of the lowest layer that held aerosol, which is below the surface, and now count from the surface altitude (over the sea the two agree). The smoke of the burning seasons is higher in them: the fine mode over southern Africa is 2.3 to 2.4 km from September to November (2.1 before) and 1.7 to 1.8 in June and July; over the Amazon 2.5 to 2.6 km in October and November (2.1 to 2.2 before); over the Sahel 2.2 km in the December to February burning (1.9 to 2.0 before) and 2.6 to 2.7 in April and May. The winter haze of India is 1.3 km (2.0 before) and its pre-monsoon aerosol 2.3 to 2.4.

On the model's seasonal N = 64 states, the day means of the shortwave with the climatology on the full-year file minus on the earlier one, all-sky top of atmosphere, W/m² (the clear-sky atmosphere in brackets):

| region | Jun 20 | Sep 20 | Dec 18 | Feb 17 |
|---|---|---|---|---|
| global | -0.02 (-0.01) | -0.05 (-0.01) | -0.06 (-0.01) | -0.07 (-0.02) |
| Sahara and Arabia | -0.36 (-0.14) | -0.49 (-0.13) | -0.26 (-0.11) | -0.26 (-0.14) |
| Sahel | +0.31 (+0.04) | -0.35 (-0.04) | -0.14 (-0.08) | +0.23 (-0.10) |
| southern Africa | -0.09 (-0.05) | +0.11 (+0.03) | -0.08 (-0.00) | -0.18 (-0.02) |
| Amazon | -0.13 (-0.04) | +0.03 (+0.00) | -0.05 (-0.01) | -0.20 (-0.04) |
| India | +0.82 (+0.08) | -0.26 (-0.03) | -0.27 (-0.09) | -0.05 (-0.05) |
| East Asia | -0.58 (-0.08) | -0.56 (-0.10) | -0.93 (-0.12) | -1.49 (-0.14) |

The heights move the sunlight by less than 0.1 W/m² in the global mean. Regionally the largest changes are the all-sky ones of the cloudy regions: East Asia's fine mode in February, 0.5 km lower, leaves 1.5 W/m² less sunlight absorbed, ten times its change in the clear sky. In the sequential slab the near infrared that meets the aerosol is what the gases above it have left, so lower aerosol takes less of it. The Sahara and Arabia on September 20 absorb 0.5 W/m² less in the all sky against the 17 W/m² of clear-sky atmospheric absorption that their aerosol adds.

### The two shortwave schemes and the clear sky

The clear column and the clear-sky pass are the sequential scheme's, and the layered scheme adds the cloud's difference from its own clear region. `test/layeredConservation.test.mjs` runs both on the clear columns of two standard atmospheres with the climatology's global mean modes and records the increments of the absorbed sunlight (the tropical column, W/m², the sequential column's increment, the layered clear region's and their ratio; the midlatitude column's ratios differ from these by 0.02 at most):

| albedo | zenith | top | atmosphere | surface |
|---|---|---|---|---|
| 0.06 | 2 | -12.57, -5.95, 0.47 | 8.92, 10.25, 1.15 | -21.50, -16.20, 0.75 |
| 0.06 | 30 | -12.33, -8.66, 0.70 | 8.78, 9.57, 1.09 | -21.11, -18.23, 0.86 |
| 0.06 | 45 | -11.86, -11.52, 0.97 | 8.59, 8.59, 1.00 | -20.45, -20.11, 0.98 |
| 0.06 | 60 | -10.78, -14.26, 1.32 | 8.25, 6.90, 0.84 | -19.03, -21.16, 1.11 |
| 0.06 | 70 | -9.26, -14.97, 1.62 | 7.87, 5.09, 0.65 | -17.12, -20.06, 1.17 |
| 0.2 | 2 | -5.83, -0.04, 0.01 | 9.93, 12.09, 1.22 | -15.76, -12.13, 0.77 |
| 0.2 | 30 | -6.21, -3.05, 0.49 | 9.63, 11.09, 1.15 | -15.84, -14.14, 0.89 |
| 0.2 | 45 | -6.51, -6.33, 0.97 | 9.25, 9.73, 1.05 | -15.76, -16.06, 1.02 |
| 0.2 | 60 | -6.49, -9.76, 1.50 | 8.67, 7.60, 0.88 | -15.16, -17.35, 1.14 |
| 0.2 | 70 | -5.86, -11.20, 1.91 | 8.11, 5.49, 0.68 | -13.97, -16.69, 1.19 |

The two agree with the sun 45 degrees from the zenith (the test asserts the top and the surface within 10 percent there), and in the atmosphere's absorption within 0.75 to 1.35 up to 60 degrees (0.84 to 1.22 measured). At the top the layered increment is smaller than the sequential one with the sun high and larger with the sun low, over bright surfaces most. The clear sky of every sunlit column is the sequential scheme's, so the clear-sky change at the top is uncertain by this much; the all-sky change, where the layered scheme's clear region carries the cloudy columns, is less so. The disagreement exists without the climatology: for the flat aerosol the layered clear region's increment is 0.85 of the sequential column's at the top and 0.71 in the atmosphere (the design study's three August N = 128 states), and its single-column test put the difference in the sequential visible slab at low sun and over bright surfaces. Carrying the layered clear region for every sunlit column would remove it and is a separate change.

### The engines

At N = 16 on a synthetic climatology (`test/aerosolGpu.test.mjs`, GPU against CPU after three steps from day 200, 99th percentile of the cells of the relative difference of the day sums, the rms in parentheses), the engines agree under `'mac'` as they do without it:

| configuration | absorbed sunlight | atmosphere | theta |
|---|---|---|---|
| layered, flat | 1.3e-5 (4.7e-6) | 6.4e-6 (4.1e-6) | 1.0e-6 (1.5e-6) |
| layered, `'mac'` | 9.0e-6 (2.7e-6) | 4.9e-6 (1.9e-6) | 1.1e-6 (9.8e-7) |
| sequential, flat | 1.4e-5 (6.3e-4) | 4.9e-6 (2.2e-4) | 1.1e-6 (2.3e-6) |
| sequential, `'mac'` | 7.6e-6 (2.7e-6) | 3.7e-6 (1.8e-6) | 1.0e-6 (2.5e-7) |
| layered, radiation held, flat | 2.1e-5 (2.2e-5) | 9.9e-6 (7.3e-6) | 2.3e-6 (4.4e-6) |
| layered, radiation held, `'mac'` | 1.8e-5 (1.1e-5) | 8.3e-6 (3.1e-6) | 2.4e-6 (5.6e-6) |
| sequential, CLIRAD near infrared, flat | 1.2e-5 (2.7e-6) | 4.6e-6 (1.1e-6) | 1.0e-6 (1.1e-6) |
| sequential, CLIRAD near infrared, `'mac'` | 7.4e-6 (2.3e-6) | 3.7e-6 (1.0e-6) | 1.0e-6 (1.3e-6) |
| twelve steps across a change of day, flat | 3.6e-4 (3.7e-4) | 1.7e-4 (2.6e-4) | 1.7e-5 (1.0e-5) |
| twelve steps across a change of day, `'mac'` | 2.9e-4 (4.6e-4) | 1.7e-4 (3.1e-4) | 2.2e-5 (1.6e-5) |

The uploaded fields equal the CPU sampler's exactly, and the area means of the absorbed sunlight agree to the printed digits. On a single full call of the radiation at N = 16 over the 1400 sunlit columns where the engines' decisions agree without the aerosol, the aerosol changes a layer's shortwave heating by up to 1.5 K/day (sequential) and 2.4 K/day (layered), and the engines' change differs by at most 7.0e-4 and 2.2e-3 K/day, 1.9e-4 and 3.0e-4 below 100 hPa, and 1.7e-4 and 3.2e-4 W/m² as a flux in any layer. Over a 4 km mountain the aerosol adds 1813.6 W/m² to the layers' heating on both engines, 20.0 percent of it more than 5 km above the ground; the engines' change differs by 8.6e-5 W/m² in a layer and 4.7e-4 W/m² in the absorbed sunlight.

At the page's state (`runs/fourteen128_day1825`, N = 128, dt 168.75 s, with the options its header records), `scripts/compareEngines.mjs` after one and two steps, with the default and with `aerosol: 'flat'`: the shortwave rows (surface and held shortwave, the sums of the absorbed sunlight, the top) differ by 3e-4 to 6e-4 in relative rms in both, PH_AER by 0, and the only rows over the tool's 1e-3 are the same in both, the longwave (1.1e-3) and the held direct beam (1.2e-3), at each of the two steps.

### Cost

On the CPU one column costs the same in the layered scheme (371 to 376 µs flat, 368 to 386 µs `'mac'`, N = 16 grid with 48 levels, two cloud layers) and 8 percent more in the sequential scheme (156 to 160 µs flat, 172 µs `'mac'`).

On the GPU at N = 128 on the M5 Max, from the page's state with its options (`runs/fourteen128_day1825`), alone on the device at level 2 of the queue, `'flat'` and `'mac'` alternating over three rounds (measured before the climatology was the default and had its heights per month, neither of which changes the kernels): one model day with the steps queued 8 to a submission (`scripts/paceGpu.mjs`) takes 47.3, 47.5 and 47.4 s flat and 48.5, 48.7 and 48.6 s with `'mac'`, 2.5 percent more. `scripts/profileGpu.mjs` gives a step 92.0 ms of GPU time flat and 94.5 ms with `'mac'` (69.2 and 70.6 ms awaited one at a time); the difference is in the radiation passes, 2.4 ms a step amortised over the radiation's cadence.


## What the option leaves out

- The aerosol is tropospheric and of one year. No stratospheric sulfate (ERA5's cycle adds 0.045 of it), no trend and no volcanic years.
- No infrared. The coarse mode's absorption at 10 µm (dust) is a later stage; it matters for the surface in the Sahara and Arabia and for a coupled model's top-of-atmosphere balance by about 0.7 W/m².
- The spectral constants are fitted to the August maps. For the fine mode's 1000 to 550 nm depth ratio and its 1 - SSA they hold in the other months to 5 percent; the other bands were not tested against other months.
- The vertical distribution is one number per mode, from ECHAM-HAM.
- The default moves the calibrated state: in August the all-sky top of atmosphere absorbs 0.5 W/m² less, the atmosphere 3.3 W/m² more and the surface 3.8 W/m² less (above), and the page's default states were spun up with the flat aerosol, so a run from them adjusts to the climatology. The A/B behind the default is at fixed sea-surface temperature, on the earlier file; a coupled run with the climatology has not been done.
