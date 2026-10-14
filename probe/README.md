# M5 miscompile probe

The physics kernels compute the boundary-layer mixing wrong on an Apple M5 Max
(gpuModel 23/31, boundaryLayer parity off by 0.37 K in θ, q halved in wet
columns, wind 2.2 m/s) and right on the M1 Max. With Dawn's toggle
`enable-dawn-features=disable_robustness` the M5 matches the CPU to rounding,
so Tint's robustness clamping miscompiles on the M5 for some access pattern.

Scripts (run from the repo root, nothing else on the GPU, `log show` afterward):

    node probe/bisect.mjs                      # the parity harness: per-layer and per-column differences
    FLAGS="enable-dawn-features=disable_robustness" node probe/withflags.mjs   # the same under a Dawn toggle
    FORM=ptr node probe/thomas.mjs             # the kernel's Thomas solve alone (forms ptr, inline, fwd)

Results from another machine go in probe/results-<machine>.md on this branch.

More probes: `COEF=1 node probe/bisect.mjs` prints the mixing coefficients the GPU kernel
holds after its solve against the expected values for two columns (needs the
instrumented mixField described in the session notes; without it the qc column is
the plain result); `node probe/lengths.mjs` and `node probe/arrlen.mjs` check
robustness clamping of buffer reads in helper functions and across eight bindings.
On the M5 Max all three small probes pass; only the full kernels fail.
