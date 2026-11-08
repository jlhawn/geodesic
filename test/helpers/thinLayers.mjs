import { EPS } from '../../js/ocean/layered.module.js';

export function thinTopLayers(ocean, ice, cuts) {
  const served = ocean.serialize(), C = served.eta.length;
  const h = Float64Array.from(served.h), T = Float64Array.from(served.T), S = Float64Array.from(served.S), Q = Float64Array.from(ocean.Q), W = Float64Array.from(ocean.W);
  for (const [i, top] of cuts) {
    h[C + i] += h[i] - top;
    h[i] = top;
    for (const n of [i, C + i]) {
      const span = Math.max(EPS, h[n]);
      [T[n], S[n]] = [(h[n] / span) * served.T[n], (h[n] / span) * served.S[n]];
      [Q[n], W[n]] = [span * T[n], span * S[n]];
    }
  }
  const state = { ...served, h, T, S, previousT0: ocean.T0, previousIceVolume: ice, fresh: new Float64Array(C), buoyancy: new Float64Array(C), flux: new Float64Array(C), capacity: ocean.capacity };
  return { Q, W, state };
}
