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

## Verifying the fix (branch m5-fix) on the M1 Max: bit identity

The fix only reverses the direction of three loops; the arithmetic is the same, so the
M1's results must not change at all. From the repo root, nothing else on the GPU:

    git fetch origin m5-fix m5-probe && git checkout m5-probe
    B=$TMPDIR/m5base; F=$TMPDIR/m5fix
    git worktree add --detach $B origin/c-grid-mesh-prep
    git worktree add --detach $F origin/m5-fix
    node probe/gpuIdentity.mjs $B; node probe/gpuIdentity.mjs $F     # must match line for line
    node probe/cpu20.mjs $B; node probe/cpu20.mjs $F                 # must match
    node probe/pagecmp.mjs $B $PWD/runs/eleven64_day1825.bin 32      # the check that matters: S, PH, D hashes after 32 page steps
    node probe/pagecmp.mjs $F $PWD/runs/eleven64_day1825.bin 32      # must match the line above
    /usr/bin/log show --last 15m --predicate 'eventMessage CONTAINS "aborted"' --style compact | tail -3
    git worktree remove $B; git worktree remove $F

Also run, in the m5fix worktree, one file at a time: test/boundaryLayer.test.mjs (12/12),
test/gpuModel.test.mjs (30/31), test/skinGpu.test.mjs (report every number: tests 1, 3 and 4
fail on the M5 by 2e-4 to 1e-3 K of warm-layer drift over 20 steps, a separate question), and
node scripts/paceGpu.mjs runs/eleven128_day1825.bin (report the day-mean ASR; the M5 gives 236.82).
Write everything to probe/results-m1.md, commit on m5-probe, push. Use /usr/bin/log, not log: in
zsh `log` is a builtin and the abort check silently prints nothing.

## Step 3: verify branch m5-books on the M1 Max (bit identity again)

A second M5-only compiler bug: in the radiation kernel the private array layerOzone, after
being passed by pointer to solarGases, reads back as zeros in the upward-ozone heating loop,
so the M5 lost about 2 W/m² of stratospheric heating (N = 16 books air −2.60 against the
M1's −0.48; the page at N = 64 lost 1.75 W/m²). The fix (one loop reads the same values from
the OUT stash instead) is on origin/m5-books, one commit on e1909fb. On the M1 the old code is
right, so the fix must be bit-identical there:

    git fetch origin m5-books m5-probe && git checkout m5-probe
    B=$TMPDIR/b3; F=$TMPDIR/f3
    git worktree add --detach $B origin/c-grid-mesh-prep      # e1909fb
    git worktree add --detach $F origin/m5-books
    node probe/gpuIdentity.mjs $B; node probe/gpuIdentity.mjs $F          # must match
    node probe/pagecmp.mjs $B $PWD/runs/eleven64_day1825.bin 32            # S, PH, D hashes
    node probe/pagecmp.mjs $F $PWD/runs/eleven64_day1825.bin 32            # must match the line above
    /usr/bin/log show --last 15m --predicate 'eventMessage CONTAINS "aborted"' --style compact | tail -3
    git worktree remove $B; git worktree remove $F

In the f3 worktree, one at a time: test/gpuModel.test.mjs (record WHICH assertion of test 17
fails on the M1: on the M5 the base fails its first CPU/GPU comparison, beta off by 2.2e-2,
and the fix fails later at "the limiter over 16 steps binds"), test/atmosphereConservation,
test/physics (34–37 are known CPU digest failures). Append to probe/results-m1.md, commit, push.

Step 3 note: origin/m5-books now has a second commit (the Node 22 guard, no kernel change); verify at
its tip. On the M5 with Node 22 and the fix, paceGpu on eleven64_day1825 gives ASR 228.10, OLR 229.36,
Ts 14.442 against the M1's 228.08 / 229.36 / 14.442, so the two machines now agree.
