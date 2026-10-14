# Working in this repository

WebGCM: a browser climate model (JavaScript; CPU and WebGPU engines under `js/`), its page (`climate.html`, `js/climate.module.js`, `js/model.worker.js`), the explainer series under `explain/`, spin-up and analysis scripts under `scripts/`, tests under `test/` (`node --test`).

## Tools and tokens

- Put several steps in one Bash call; one command per turn is the expensive habit. Never poll in a loop of separate calls: use one `until … ; do sleep 20; done`, or run the job in the background and read its log once.
- Keep tool output short: `head`, `tail`, `grep`, `wc`, `cut -c1-200`. Send long output to a file and read the part you need. Read large files in ranges (`sed -n a,bp`, or Read with offset and limit); `docs/c-grid-dynamical-core.md` is 13,000 lines and `js/gpu/physics.gpu.js` 2,600.
- Do not `cat` binary files or state files (`runs/*.bin`, `data/*.bin`, `explain/data/*.bin`); read them through the modules that decode them.
- Report with file paths and a few lines, not pasted dumps. Finish in under about 80 tool calls; if the job needs more, stop and report what remains and why.

## The GPU

- One GPU model at a time on this Mac: macOS aborts command buffers under contention, silently. Before any GPU run (a test file that steps a GPU model, `scripts/spinup.mjs`, `oceanSpinup.mjs`, `explainOcean.mjs`, a profile) check that no page server listens on port 8000 (`lsof -nP -iTCP:8000 -sTCP:LISTEN`) and no other node job is running (`ps -axo pid,command | grep -E 'node (scripts/|--test)'`); if one is, wait, polling every 30 s. Afterward: `/usr/bin/log show --last 15m --predicate 'eventMessage CONTAINS "aborted"' --style compact | tail -3`.
- Timings taken under contention are worthless; alternate base and candidate runs and repeat them. The page loop on the M1 Max is already about 96 percent GPU-busy.
- Run GPU test files one at a time; CPU test files together as concurrent runners.

## Conventions

- Never commit `socol.txt`, anything under `runs/` beyond the page's default state parts, or scratch files.
- American spelling in prose, UI text, messages and comments. Code comments default to none and never narrate a change.
- The static server on http://localhost:8001 serves the repo root; do not start another, and do not open a GUI browser from an agent.
- Both engines must stay in step: a model change lands on the CPU (f64) and the WGSL (f32) together, with a parity test at N = 8 and, where it touches storage or transport, a row in `test/layeredConservation.test.mjs` or `test/atmosphereConservation.test.mjs`.
- Saved states record their grid (`grid: { relax }`); loaders use `js/stateGrid.module.js`.
