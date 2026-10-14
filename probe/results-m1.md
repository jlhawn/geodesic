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

## Step 2 (m5-fix)

Pending: the branch m5-fix has not appeared on origin yet. The "verify" section below is to be
filled in by the M5 session; the M1 will run it when both are there.

### verify

(to be added by the M5 session)
