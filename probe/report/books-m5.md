# The M5's 2 W/m² in the air's books: the upward ozone heating read as zeros

Branch `m5-books` (off `c-grid-mesh-prep` at e1909fb), one commit, not pushed:

- 3ee0c63 Read the upward ozone heating's layers from the stash so the Apple M5 deposits it

Only `js/gpu/physics.gpu.js` changes: one loop in the radiation kernel's upward block, plus a one-line comment.

## Which term moved

On the M5, the energy that went missing was the air's own heating. The fluxes the recorder reads (ASR, OLR, the surface net, the snow) match the M1's. Each step's change in c_p T + L q splits into three parts (one N = 16 day, M5, before the fix):

| part of the step | Δ(c_p T + L q) | Δ KE | expected |
|---|---|---|---|
| dynamics (RK) | −0.994 | +1.002 | the sum ≈ 0 |
| physics passes (physics, radiation, longwave, physicsSurface, pblDiagnose) | 342.590 | 0 | TOA net − surface net + snow = 344.802, **short by 2.212** |
| closure, adjust, mixMomentum, dissipation | +0.605 | −0.520 | the sum ≈ 0 |
| total against the fluxes | −2.602 | | the recorder's "air −2.60" |

Inside the physics passes, with the radiation held over two steps (radiationEvery 2 at 450 s), the held profiles were checked against their boundary fluxes:

- Longwave: Σ_k LWH equals −OLR + emission − back-radiation to 1e-3 W/m², and Σ_k RADDF equals 1 − T.
- Shortwave: **Σ_k RADSW is 0.84–0.86 W/m² short of RADATM**, the absorbed sunlight the air reports. RADATM equals RADTOA − RADABS.

Stashing Σ netFlux after each stage of the shortwave puts all of the gap in the upward beam's ozone term:

- ozone, gas, aerosol, cloud and near-infrared upward: 0.000
- upward ozone: −0.838 of the upward heating's 1.31 W/m²

## The kernel and the loop

`radiationKernel` (held and unheld share the text), upward block:

```wgsl
var layerOzone: array<f32, K>;                                   // K = 36
for (var k = 0; k < K; k++) { layerOzone[k] = OUT[LWS_OZONE + k * C + i]; }
...
if (SOLAR_CLIRAD && mu > 0.0) { solar = solarGases(i, pi, mu, &layerOzone, &ozoneTaken, &gasTaken); }   // reads it through a pointer: correct
...
if (ozoneLoss > 0.0) {
  for (var k = 0; k < K; k++) { netFlux[k] += upwardOzone * layerOzone[k] / ozoneColumn; }   // M5: every layerOzone[k] reads 0.0
  upwardHeating += upwardOzone;
}
```

Probes on the M5:

- A stash added after this loop read Σ layerOzone / ozoneColumn = 0.0000 in all 1292 sunlit columns, with all 36 elements equal to 0.0.
- The stash the array is copied from is right. Σ_k OUT[LWS_OZONE] equals stash 24 (ozoneColumn) to 1e-8, and no layer is zero.
- `solarGases` sees the right values through the pointer, because the downward ozone heating and the ASR match the M1's.
- The upward ozone absorption entered RADATM, RADTOA and ASR but never heated the air.

### Why the toggle and the down-counting scan missed it

- With `disable_robustness` the books are the same (−2.60). This loop's bounds are the constant K, so Tint emits no clamp either way. The generated MSL is a plain up-counting `while` loop reading `layerOzone[k]`.
- It is not the down-counting pattern. `nearInfraredUpward`'s down-counting loop, rewritten to count up, changed nothing.
- A standalone probe (`books/repro_ptr.mjs`: an array filled in main, read by a pointer-taking function, then read again in main) compiles correctly for K = 32, 33 and 36. The trigger needs more of the real kernel, which has twelve K-sized private arrays.

## The form that survives

The loop now reads the stash the array was filled from, with the same arithmetic:

```wgsl
for (var k = 0; k < K; k++) { netFlux[k] += upwardOzone * OUT[LWS_OZONE + k * C + i] / ozoneColumn; }
```

`OUT[LWS_OZONE]` is written only by the physics kernel, so the values are bit-for-bit those the array was filled with.

Bits on the M5:

- **Page state, N = 64, 32 steps (S, PH, D).** The fix and the old code with `layerOzone` refilled from the stash just before the loop give the same hashes: 0b72836af058 / e2c5f944a053 / fe6d1f13a41f. The old code gives 2c6a75aeb7ca / 2ea9c2a5c4be / 59d17fa78e34.
- **gpuIdentity, base and fix the same:**
  - coupled N = 6, stepBatch(33): S 3358de065d0f, PH 10c437f9a086, D 8d6f745c501b
  - ocean N = 8: OS 2e144e1be10f
  - These configurations take the radiation every step. That kernel either compiled correctly on the M5 or yields the same bits.

## Before and after

**N = 16 water planet, wctl options** (`books/wctl.json`), `scripts/explainOcean.mjs` from day 0. The residual is the whole run's, as the log line prints it.

| | day 1: residual / surface / air (W/m²) | day 1: ASR / OLR | day 365: residual / air |
|---|---|---|---|
| M1 a50all (d27c5dc, `diurnal/gpu/a50all.log`) | +0.58 / 1.06 / −0.48 | 205.6 / 191.6 | +0.31 / −0.03 |
| M5 wctl (2391a70, `diurnal/gpu/wctl.log`) | −1.55 / 1.05 / −2.60 | 205.6 / 191.5 | −1.60 / −1.95 |
| M5 e1909fb, robustness on | −1.54 / 1.07 / −2.60 | 205.7 / 191.5 | |
| M5 e1909fb, robustness off | −1.54 / 1.07 / −2.60 | 205.6 / 191.5 | |
| M5 d27c5dc (the M1's commit), robustness off | −1.54 / 1.06 / −2.60 | 205.7 / 191.5 | |
| M5 m5-books, robustness on | +0.58 / 1.07 / −0.49 | 205.6 / 191.6 | +0.28 / −0.03 (`books/fix365.log`) |
| M5 m5-books, robustness off | +0.58 / 1.06 / −0.49 | 205.6 / 191.5 | |

- Over the year, m5-books follows a50all day by day: the mean of the logged residual / air is 0.111 / −0.167 against 0.127 / −0.161, and the old M5 gave −1.786 / −2.133.
- d27c5dc on the M5 with robustness on blows up (−50 W/m²): that is the Thomas bug fixed in 8a3f033.
- Under the fix, robustness on and off differ in the last bits (relaxed math) and agree in the books to 0.01 W/m².

**The page is affected too.** N = 64 also holds its radiation over two steps. After 2 steps from `runs/eleven64_day1825.bin`, Σ RADSW − RADATM is −1.748 W/m² on e1909fb and 0.000 on m5-books (`books/pagesw.mjs`). In Chrome on the M5 the page's air has been losing about 1–2 W/m² of upward ozone heating, mostly in the stratosphere. ASR does not show it.

**One-day check** (`node scripts/paceGpu.mjs runs/eleven64_day1825.bin`):

| | ASR | OLR | Ts |
|---|---|---|---|
| e1909fb (M5) | 227.96 | 229.29 | 14.436 °C |
| m5-books (M5) | **227.94** | 229.33 | 14.436 °C |

The 0.02 drop in ASR is the stratosphere receiving its upward ozone heat. The M1 gives 228.0.

## Tests on the M5 Max (node 18, webgpu 0.6.1), each GPU file alone, after the fix

TAP files are in `books/tests/`.

| file | result |
|---|---|
| gpuModel | 30/31. Test 17 (drags laid every fourth step over terrain) fails, as it does on the M1 and before the fix (fix.md). |
| atmosphereConservation | 6/6 |
| physics | 33/37. Tests 34–37 are CPU-only twelve-step digests and fail identically on e1909fb (`physics.base.tap`). They are an M5 libm/f64 matter, not this change. |
| frameGpu | 4/4 |
| boundaryLayer | 12/12 |

`/usr/bin/log show` found no command-buffer aborts.

## Verifying on the M1 Max (bit identity)

The old code is wrong on the M5 with robustness on or off, so the bit check of old against new has to run on the M1.

1. Get `m5-books` onto the M1. It is not pushed: push it from the M5, or fetch it.
2. Copy these scripts from the M5 scratchpad. They take the repo path as an argument and load webgpu from `/Users/jlhawn/git_repos/jlhawn/geodesic/node_modules`:
   - `m5probe/pagecmp.mjs`
   - `momwork/gpuIdentity.mjs`
   - `m5probe/books/eo.mjs`, `wctl.json`, `hashState.mjs` and `pagesw.mjs`

```sh
cd /Users/jlhawn/git_repos/jlhawn/geodesic
git fetch origin m5-books   # or however the branch arrives
B=$TMPDIR/booksbase; F=$TMPDIR/booksfix
git worktree add --detach $B e1909fb
git worktree add --detach $F m5-books
ln -s $PWD/node_modules $B/node_modules; ln -s $PWD/node_modules $F/node_modules
# one GPU job at a time; check port 8000 and other node jobs first
node pagecmp.mjs $B $PWD/runs/eleven64_day1825.bin 32    # page state, held radiation: S, PH, D hashes
node pagecmp.mjs $F $PWD/runs/eleven64_day1825.bin 32    # must match the line above (the M5 fix gives S 0b72836af058)
node gpuIdentity.mjs $B; node gpuIdentity.mjs $F         # must match line for line (M5: S 3358de065d0f, OS 2e144e1be10f)
export OCEAN_OPTS="$(cat wctl.json)"
CHECKPOINT=$TMPDIR/bb_{day}.bin node eo.mjs $B 3 0 $TMPDIR/bb.bin 2>&1 | grep 'ocean day'
CHECKPOINT=$TMPDIR/bf_{day}.bin node eo.mjs $F 3 0 $TMPDIR/bf.bin 2>&1 | grep 'ocean day'   # same lines; air ≈ −0.49, −0.46, −0.38
node hashState.mjs $TMPDIR/bb_00003.bin $TMPDIR/bf_00003.bin                              # must match
node pagesw.mjs $B $PWD/runs/eleven64_day1825.bin; node pagesw.mjs $F $PWD/runs/eleven64_day1825.bin   # both 0.000 on the M1
/usr/bin/log show --last 15m --predicate 'eventMessage CONTAINS "aborted"' --style compact | tail -3
git worktree remove --force $B; git worktree remove --force $F; rm -f $TMPDIR/bb* $TMPDIR/bf*
```

The `pagecmp` pair is the check that matters, because the page holds its radiation as the N = 16 run does. If it differs on the M1, relaxed math reshaped the loop there. In that case compare the two S buffers with `fielddiff.mjs` before trusting either. The `eo.mjs` lines should also reproduce a50all's day 1 (air −0.48).

## How it was found (to reuse)

- **Books per stage.** `books/seg.inc.mjs` is spliced into a scratch copy of `scripts/explainOcean.mjs`, after the `const global = …` line, with `awk`. It wraps `model.step` and the core's `beforePhysics` and `afterPhysics` hooks. Each `seg day` line splits Δ(c_p T + L q) and ΔKE into dynamics, physics passes and the rest, and sets them against TOA − surface + snow.
- **Held-profile checks.** `held.inc.mjs` checks the held profiles against their boundary fluxes: SW Σ RADSW against RADATM, LW Σ LWH against −OLR + emission − back-radiation, and Σ RADDF against 1 − T.
- **Shortwave stages.** `dbg.inc.mjs` and the kernel stashes in `physics.dbg.js` (slots 26–36 of the K1 stash) give Σ netFlux after each shortwave stage.
- **The array the loop reads.** `oz.inc.mjs` checks the stashed ozone. `oz2.inc.mjs` reads what the upward loop saw (slots 37–38).
- **Runner.** `eo.mjs` runs any repo's (or worktree's) explainOcean with Dawn toggles from FLAGS. SCRIPT picks the instrumented copy.
- **Comparison.** `cmp.mjs` sets run logs side by side, and `engineOf.mjs` prints a checkpoint's engine commit.
- **Shader dump.** `n16.dump` holds the WGSL and MSL of the N = 16 configuration, with `dump_shaders,disable_symbol_renaming`.

## Still open

The trigger is not reduced to a standalone shader. A private K-sized array that is written in main, read through a pointer by a helper, and read again in main also occurs in these places, each of which this round did not cover:

- The orography kernel (`z`, `th` through `oroN2`, with `z[k]` read afterward in main). It is off on the water planet. On the page its robustness check cannot catch this bug class, since the toggle does not matter here. A CPU/GPU comparison of the orography drag at N = 64 on the M5 would.
- `mixShare(&rhs)` in mixMomentum.
- The gray longwave `band()` path, which the page does not use.
