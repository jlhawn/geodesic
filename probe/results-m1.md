# M1 Max results

Machine: MacBook Pro, Apple M1 Max, 32 GB, macOS 26.7.1, Node v26.7.0, webgpu 0.6.1 (Dawn, Metal).
Branch m5-probe at 90de321 (= c-grid-mesh-prep 421662b + probe/). Each job ran alone on the GPU
(port 8000 free, no other node job); `log show --last 15m … "aborted"` found nothing afterward.
Raw outputs: probe/out/ (not committed).

## Step 1 (2026-10-14 03:12–03:14 UTC)

| Job | Result |
|---|---|
| `node --test test/boundaryLayer.test.mjs` | 12/12 |
| `node --test test/gpuModel.test.mjs` | 30/31; the one failure is test 17 (drags every fourth step over terrain), the known one |
| `node probe/bisect.mjs` | theta: 0 of 362 columns differ, worst 1.01e-4 K at column 20 (per-layer worst 1e-5 to 1e-4 K through the mixed layers, 0 above); q: 0 of 362 differ, worst 1.96e-7 at column 20; qc: 0 differ, worst 5.48e-13; wind: worst 2.09e-4 m/s at k 27 edge 533, edges differing > 1e-3: 0 of 36720. Column 20's mixing coefficients agree layer by layer (mix cpu/gpu identical to the printed digits; w_e 1.16e-3 on both) |
| `node scripts/paceGpu.mjs runs/eleven64_day1825.bin` | N = 64, 256 steps (1 d) awaited one at a time: 6.7 s per model day, 26.2 ms per step; **day mean ASR 228.08**, OLR 229.36 W/m², Ts 14.442 °C |

So the M1 computes the boundary-layer mixing as the CPU does (the README's 0.37 K θ gap on the M5
is 1e-4 K here) and gives the right ASR (228.08 against the M5's 259.6).

## Step 2: m5-fix at 2391a70 (2026-10-14 03:45–03:52 UTC)

m5-fix = c-grid-mesh-prep 421662b + three commits: 8a3f033 (thomas/thomasRate back-substitution
counts up), dadcbb4 (mixEdge's wind solve likewise), 2391a70 (CLAUDE.md: /usr/bin/log). The whole
engine diff is 8 lines of js/gpu/physics.gpu.js: each `for (var j = n - 2; j >= 0; j--)` becomes
`for (var m = 0; m < n - 1; m++) { let j = n - 2 - m; … }` with the same body. Each GPU job ran
alone; `/usr/bin/log show` found no aborted command buffer (the only "aborted" hit is the log
command's own line).

| Job | m5-fix on the M1 Max | m5-probe (step 1) |
|---|---|---|
| `test/boundaryLayer.test.mjs` | 12/12 | 12/12 |
| `test/gpuModel.test.mjs` | 30/31, test 17 as before | 30/31 |
| `test/atmosphereConservation.test.mjs` | 6/6 | – |
| `test/skinGpu.test.mjs` | 6/6 | – |
| `test/layeredGpu.test.mjs` | 21/21 | – |
| `probe/bisect.mjs` | theta 0 of 362 differ, worst 1.01e-4 at column 20; q worst 1.96e-7; qc worst 5.48e-13; wind worst 2.09e-4 at k 27 edge 533, 0 of 36720 edges > 1e-3 | the same digits |
| `paceGpu` eleven64_day1825 | day mean ASR 228.08, OLR 229.36 W/m², Ts 14.442 °C; 6.6 s per model day | ASR 228.08, OLR 229.36, Ts 14.442 |

### verify: bit identity, m5-fix against c-grid-mesh-prep (421662b)

Reference 421662b in a detached worktree; the fix in the main checkout. Scripts: the session
scratchpad's momwork/cpu20.mjs (the CPU ocean after 20 steps at N = 8, synthetic and real
topography) and momwork/gpuIdentity.mjs (the coupled GPU model at N = 6 after stepBatch(33) as
test/stepBatch.test.mjs prepares it; the N = 8 GPU ocean of layeredGpu's scenario after 1 + 20
advances), each a SHA-256 of the state buffers.

| Hash | 421662b | 2391a70 (m5-fix) |
|---|---|---|
| CPU 20 steps N = 8, synthetic | d6e57038e6701b8f | d6e57038e6701b8f |
| CPU 20 steps N = 8, real topography | 6ce79c14425c57a6 | 6ce79c14425c57a6 |
| GPU coupled N = 6 stepBatch(33): S / PH / D / OS / OD | ef02e2aee9ca / 674f74199bfb / 2432f16e14b9 / 0451a6f1246d / 1479957dc347 | identical |
| GPU ocean N = 8, 21 advances: OS / OD / PH | e6ffe7267fc9 / a8b028262b02 / 2b991f687529 | identical |

The check the M5 session's fix-m5.md calls the one that matters, `probe/pagecmp.mjs` on the page's
own state (eleven64_day1825, N = 64, 32 steps; at N = 6 the boundary layer has no top, so only
this path runs the changed loops with n < K):

| Hash after 32 steps | 421662b | 2391a70 (m5-fix) |
|---|---|---|
| S / PH / D | 338567a3639b / b7e7737821f2 / 144631da4b8b | identical |

**The fix changes nothing on the M1 Max**: every state hash is identical, every test count is
the same, the bisect harness prints the same digits, and the day-mean ASR is 228.08 in both.
m5-fix can be merged.

## m5-fix, the extra items the README asks for (2026-10-14 17:00–17:08 UTC, m5-fix worktree at 2391a70)

| Job | M1 Max |
|---|---|
| `test/skinGpu.test.mjs` | 6/6. Parity after one step, the four cases (slab from no skin; slab from a 1 K warm layer; layered ocean with 2 m mixed layers under IB light; layered ocean under type III light): ΔT_w 1.4e-6 / 1.0e-5 / 1.3e-5 / 7.1e-6 K, ΔT_c 3.7e-6 / 4.2e-6 / 4.3e-6 / 4.9e-6 K, T_s 2.9e-5 / 3.2e-5 / 3.3e-5 / 3.3e-5 K. Worst over 20 steps: ΔT_w 9.7e-6 / 6.9e-5 / 1.4e-5 / 1.5e-5 K (the M5's 2e-4 to 1e-3 K drift does not occur here), ΔT_c 5.7e-3 / 1.1e-2 / 5.7e-3 / 1.0e-2 K, T_s 5.9e-3 / 4.0e-2 / 1.3e-2 / 1.2e-2 K, with the bulk temperatures drifting 5.3e-4 / 4.0e-2 / 1.3e-2 / 1.2e-2 K apart; net 2.3e-2 / 1.9 / 1.0e-1 / 5.5e-1 and sensible 1.2e-3 / 1.5e-2 / 2.0e-3 / 4.8e-3 relative |
| `paceGpu` on eleven128_day1825 (N = 128, 512 steps awaited one at a time) | 45.4 s per model day, 88.7 ms per step; **day mean ASR 237.00**, OLR 234.60 W/m², Ts 14.766 °C (the M5 gives 236.82) |

## Step 3: m5-books at ad60082 (2026-10-14 17:00–17:04 UTC)

Base b3 = origin/c-grid-mesh-prep e1909fb (m5-fix merged, Part 8's run recorded); fix f3 =
origin/m5-books ad60082 (94d4285: the upward ozone heating reads its layers from the OUT stash;
ad60082: the Node 22 guard, no kernel change). Detached worktrees, each GPU job alone, `/usr/bin/log`
found no aborted command buffer (the only hit is a Wi-Fi auto-join line).

| Hash | e1909fb | ad60082 (m5-books) |
|---|---|---|
| GPU coupled N = 6 stepBatch(33): S / PH / D / OS / OD | ef02e2aee9ca / 674f74199bfb / 2432f16e14b9 / 0451a6f1246d / 1479957dc347 | identical |
| GPU ocean N = 8, 21 advances: OS / OD / PH | e6ffe7267fc9 / a8b028262b02 / 2b991f687529 | identical |
| `pagecmp` eleven64_day1825, 32 steps: S / PH / D | 338567a3639b / b7e7737821f2 / 144631da4b8b | identical |

All three sets of hashes are also identical to those of 421662b and 2391a70 above: neither fix
changes a bit on the M1.

Tests in the f3 worktree, one file at a time:

| File | M1 Max |
|---|---|
| `test/gpuModel.test.mjs` | 30/31. Test 17 fails at its last assertion, "the limiter over 16 steps binds" (AssertionError, expected true, actual false, at test/gpuModel.test.mjs:535), the same place the fix fails on the M5; the earlier CPU/GPU comparisons of the test (beta and the drags) pass here |
| `test/atmosphereConservation.test.mjs` | 6/6 |
| `test/physics.test.mjs` | 33/37; 34–37 fail (the known CPU digest cases: the mixed-layer deck's bit-identity, the moist defaults' digest, the fog scheme's digest, the deck's sunlight) |

**m5-books changes nothing on the M1 Max** and can be merged.
