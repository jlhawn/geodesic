import { initializeState } from '../../js/physics/init.module.js';

export const SHARE = (i) => 0.1 + 0.8 * ((i * 0.6180339887) % 1);

/*
 * A layer mixed to σ 0.85 under a free troposphere warmed by `jump`, with
 * the first layer above the mixed layer a mass-weighted mix holding
 * `share(i)` of free-tropospheric air over the mixed layer's; the deck is
 * seeded sinking at `seed` and at `height`, above any ceiling.
 */
export function deckProfile(model, { jump = 10, share = SHARE, height = 5000, seed = -1e-3 } = {}) {
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const { K, sigmaMid } = model.core, C = model.mesh.nCells, theta = model.state[1], q = model.state[4];
  let mixed = K - 1;
  while (mixed > 0 && sigmaMid[mixed - 1] >= 0.85) mixed--;
  const ambiguous = mixed - 1;
  for (let i = 0; i < C; i++) {
    const bottom = (K - 1) * C + i;
    for (let k = 0; k < K - 1; k++) {
      const idx = k * C + i;
      if (k >= mixed) { theta[idx] = theta[bottom]; q[idx] = q[bottom]; } else theta[idx] += jump;
    }
    const idx = ambiguous * C + i, f = share(i);
    theta[idx] = (1 - f) * theta[bottom] + f * theta[idx];
    q[idx] = (1 - f) * q[bottom] + f * q[idx];
  }
  model.radiation.mlmSubsidence.fill(seed);
  model.radiation.mlmHeight.fill(height);
  return { ambiguous, mixed };
}
