# The Apple M5 cloud bug: a Metal miscompile under Tint's index clamps

Branch `m5-fix` (off `c-grid-mesh-prep` at 421662b), two commits:

- 8a3f033 Count the column Thomas solves' back-substitution up so the Apple M5 compiles it (`thomas`, `thomasRate`)
- dadcbb4 Count the edge wind solve's back-substitution up so the Apple M5 compiles it (`mixEdge`)

Only `js/gpu/physics.gpu.js` changes: three loops plus two one-line comments.

## The construct

The miscompile needs all of the following:

- a function-scope (private) array of **more than 32** `f32` (K = 33 or more; at K = 32 the code is right),
- written by a forward **recurrence** loop (`gain[j]` from `gain[j - 1]`),
- then read by a **down-counting** loop (`for (var j = m - 1; j >= 0; j--)`),
- with loop bounds that become **compile-time constants only after the Metal compiler inlines** them (a function parameter, or a `let` of a constant expression). Tint cannot see them as constant, so it keeps its robustness clamps `arr[min(uint(j), K - 1u)]`.

On the M5 Max the down-counting loop then reads back the array's zero initializer for every element except the last one written. If any bound is a runtime value, or if Dawn's `disable_robustness` toggle removes the clamps, the code is right. It is also right when the same indices are visited by an up-counting loop. When Tint can prove the bounds constant (literal or `const` bounds), it drops the clamps itself, which is why most of the model's down-counting loops are unaffected.

In the model this is the back-substitution of the tridiagonal solves: `thomas()` (adjustMix: `mixConserved` and `mixField`, so θ, q and qc) and `mixEdge()` (mixMomentum: the wind). Both have n = K − KTOP < K whenever the boundary layer has a top (`kTop` > 0, as on the page). The skipped back-substitution left q about half of the CPU's value through wet mixed layers. The deck closure then found nothing to hold, so clouds and the deck overlays emptied.

## The rewrite

```wgsl
for (var j = n - 2; j >= 0; j--) { (*rhs)[j] -= gain[j] * (*rhs)[j + 1]; }
// became
for (var m = 0; m < n - 1; m++) { let j = n - 2 - m; (*rhs)[j] -= gain[j] * (*rhs)[j + 1]; }
```

The rewrite keeps the same index sequence and the same operations in the same order. On the M5 under `disable_robustness`, `c-grid-mesh-prep` and `m5-fix` step to identical bits:

| run | `c-grid-mesh-prep` | `m5-fix` |
|---|---|---|
| page N = 64 state, 32 steps (S / PH / D) | 5f26e5995e37 / 3961493bde68 / dec655853b43 | the same hashes |
| gpuIdentity coupled N = 6, stepBatch(33) (S) | 6ae82179d5b6 | 6ae82179d5b6 |
| gpuIdentity ocean N = 8 (OS) | 2e144e1be10f | 2e144e1be10f |

Other loop forms tested on the real kernel or the probes:

- u32 counters: still wrong.
- `j > 0` with `j - 1` indexing: still wrong.
- Arrays sized to n: right for n ≤ 32 but wrong at n = 33, so it is not a fix.
- Flat, inlined `const` bounds: right only because Tint drops the clamps.

The up-counting loop is right for every K in {33, 34, 36, 48, 64} and n from 8 to 34 (`repro.mjs`, BACK=up).

## Tests on the M5 Max (node 18, webgpu 0.6.1)

"Before" is `c-grid-mesh-prep` from the probe session. "After" is `m5-fix`, each GPU file run alone.

| file | before | after |
|---|---|---|
| boundaryLayer | 11/12 (test 11: θ 0.37 K, q 7.5e-3, wind 2.2 m/s) | 12/12 (θ 6.4e-5 K, q 4.9e-8, wind 2.1e-4 m/s) |
| gpuModel | 23/31 (21–25, 27, 31 and 17) | 30/31 (17 only, fails on the M1 too) |
| atmosphereConservation | 5/6 (6) | 6/6 |
| skinGpu | 2/6 | 3/6, see below |
| convection | 41/41 | 41/41 |
| frameGpu | 4/4 | 4/4 |
| cloudHeights | 6/6 | 6/6 |
| layeredGpu | 21/21 | 21/21 |
| layeredConservation | 15/15 + 1 todo | 15/15 + 1 todo |
| stepBatch | 3/3 | 3/3 |
| layeredSunlight | 9/9 | 9/9 |

- **skinGpu 1, 3 and 4 are a separate issue.** They fail with identical numbers with robustness off: the warm layer drifts 2.0e-4, 1.0e-3 and 2.2e-4 K over 20 steps against a 1e-4 K tolerance, and the bulk temperatures drift 4e-3 to 4e-2 K apart. The first step passes. This is an M5-versus-M1 numeric difference in the coupled drift, not this bug, and is still open.
- **One-day check:** `node scripts/paceGpu.mjs runs/eleven64_day1825.bin` gives a day-mean ASR of **227.96 W/m²**, OLR 229.29 and Ts 14.436 °C (3.4 s per model day at N = 64). It printed 259.6 with the bug; the M1 gives 228.0.
- `log show` showed no command-buffer aborts.

## No other kernel is affected at the page configuration

- **Robustness on against off, page N = 64 state, one step.**
  - Before the fix: q differs by more than 1e-2 of its field scale in 301,732 values, θ by more than 1e-4 in 52,163, and the wind by more than 1e-2 in 570,810.
  - After the fix: nothing beyond 1e-3 of scale in θ, q or the wind. One column of CUMF and two values of CUFALL differ by more than 1e-2.
  - Dawn compiles the MSL under `#pragma METAL fp math_mode(relaxed)`, so the code shape changes the rounding. These leftovers are that rounding flipping a few marginal convection triggers. Tools: `pagecmp.mjs`, `dist.mjs` and `fielddiff.mjs`.
- **A scan of the page's generated MSL** (`downloops.py page.dump`) for down-counting loops that read clamped private arrays found three more with bounds known at compile time:
  - the `cuReserve` loop in diagnosedPlume,
  - `transportEdge`,
  - four loops in the orography kernel (three of them exit through `break`).
- Rewriting those loops the same way left the robustness-on result bit-identical, so they compile correctly now. Under `disable_robustness` the orography rewrite changed OBETA at 1e-7 relative (relaxed math). To keep the M1's bits, those rewrites were reverted and are not on the branch. The other matches have a runtime bound (`source`, `top`, `src`), which compiles correctly.

## Verifying on the M1 Max (bit identity)

1. Get `m5-fix` onto the M1. It is not pushed yet: push it from the M5, or fetch it.
2. Copy these scripts from the M5 scratchpad. They take the repo path as their first argument and load webgpu from `/Users/jlhawn/git_repos/jlhawn/geodesic/node_modules`:
   - `m5probe/pagecmp.mjs`
   - `momwork/gpuIdentity.mjs`
   - `momwork/cpu20.mjs`

```sh
cd /Users/jlhawn/git_repos/jlhawn/geodesic
git fetch origin m5-fix   # or however the branch arrives
B=$TMPDIR/m5base; F=$TMPDIR/m5fix
git worktree add --detach $B c-grid-mesh-prep
git worktree add --detach $F m5-fix
# one GPU job at a time; check port 8000 and other node jobs first
node gpuIdentity.mjs $B; node gpuIdentity.mjs $F          # N = 6 stepBatch(33) and N = 8 ocean hashes: must match line for line
node cpu20.mjs $B; node cpu20.mjs $F                      # CPU ocean after 20 steps at N = 8: must match (the CPU is untouched)
node pagecmp.mjs $B $PWD/runs/eleven64_day1825.bin 32     # page state, 32 steps: S, PH, D hashes
node pagecmp.mjs $F $PWD/runs/eleven64_day1825.bin 32     # must match the line above
log show --last 15m --predicate 'eventMessage CONTAINS "aborted"' --style compact | tail -3
git worktree remove $B; git worktree remove $F
```

The `pagecmp` pair is the check that matters. At N = 6 the boundary layer has no top (n = K) and the old loop compiled correctly even on the M5. On the M5 with robustness on, `gpuIdentity` gives base = fix (S 3358de065d0f), while the page state gives fix S 2c6a75aeb7ca and the base differs because it hits the bug. On the M1 all three pairs should be equal. If the page pair differs, relaxed math reshaped something on the M1; compare with `fielddiff.mjs` before trusting either.

## Bug report (Apple Feedback, and Dawn/Tint for awareness)

**Title:** Apple M5 Metal compiler: a down-counting loop over a private array of more than 32 floats reads back the zero initializer when its indices are clamped with `min(uint(j), N - 1u)` and the bounds fold to constants after inlining.

**Environment:**

- Apple M5 Max, macOS 26.7.1.
- Chrome, and Dawn via the npm `webgpu` 0.6.1 package, adapter "apple metal-3".
- The same code is right on an M1 Max. The page is fine in Safari on the M5.

**Reproduction:** `node m5probe/repro_min.mjs`, a self-contained WebGPU compute shader:

```wgsl
fn solve(i: i32, m: i32) {
  var gain: array<f32, 34>;
  gain[0] = L[i];
  for (var j = 1; j < m; j++) { gain[j] = L[j * N + i] * gain[j - 1]; }
  for (var j = m - 1; j >= 0; j--) { OUT[j * N + i] = gain[j]; }
}
// main: solve(i, 15)
```

**Expected:** `OUT` holds the running product.

**Actual (M5, robustness on):**

- 57,344 of 61,440 values are wrong. Rows j = 0 to 13 read back as 0; only j = 14 is right.
- With `enable-dawn-features=disable_robustness`: 0 wrong.
- With the read loop counting up (`FORM=up`): 0 wrong.
- With an array of 32 (`K=32`): 0 wrong.
- With any runtime loop bound: 0 wrong.

**Attachments:**

- `m5probe/min_robust_on.msl`, the failing MSL. It has the clamps `gain[min(uint(j), 33u)]`, a `uint2 tint_loop_idx` guard on the forward loop and `#pragma METAL fp math_mode(relaxed)`.
- `m5probe/min_robust_off.msl`, the same shader without clamps, which compiles correctly.

The bug is in Apple's compiler, since the MSL is valid and in-bounds. A Dawn/Tint report could ask whether the clamps could be elided when the bound is a `let` of a constant expression or a parameter passed one. That would hide this pattern in Chrome.

**Further probes:**

- `repro.mjs`: the model's Thomas solve against a float32 JS twin. Forms: ptr, sized, flat, noback; BACK=down|up|unsigned|above; K and NN.
- `onoff.mjs` with `forms/*.wgsl`: robustness on against off on any kernel body.
- `downloop.mjs`: a plain down loop with no recurrence, which compiles correctly.
