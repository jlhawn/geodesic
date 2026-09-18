export function createRK4(nCells, nEdges) {
  const stages = [0, 1, 2, 3].map(() => ({ h: new Float64Array(nCells), u: new Float64Array(nEdges) }));
  const trial = { h: new Float64Array(nCells), u: new Float64Array(nEdges) };

  function advance(h, u, stage, dt) {
    for (let i = 0; i < nCells; i++) trial.h[i] = h[i] + dt * stage.h[i];
    for (let e = 0; e < nEdges; e++) trial.u[e] = u[e] + dt * stage.u[e];
  }

  return function step(tendency, h, u, dt) {
    const [k1, k2, k3, k4] = stages;
    tendency(h, u, k1.h, k1.u);
    advance(h, u, k1, dt / 2);
    tendency(trial.h, trial.u, k2.h, k2.u);
    advance(h, u, k2, dt / 2);
    tendency(trial.h, trial.u, k3.h, k3.u);
    advance(h, u, k3, dt);
    tendency(trial.h, trial.u, k4.h, k4.u);
    const w = dt / 6;
    for (let i = 0; i < nCells; i++) h[i] += w * (k1.h[i] + 2 * k2.h[i] + 2 * k3.h[i] + k4.h[i]);
    for (let e = 0; e < nEdges; e++) u[e] += w * (k1.u[e] + 2 * k2.u[e] + 2 * k3.u[e] + k4.u[e]);
  };
}

export function createRK4Arrays(lengths) {
  const allocate = () => lengths.map((n) => new Float64Array(n));
  const stages = [allocate(), allocate(), allocate(), allocate()];
  const trial = allocate();

  function advance(state, stage, dt) {
    for (let a = 0; a < state.length; a++) {
      const s = state[a], k = stage[a], t = trial[a];
      for (let i = 0; i < s.length; i++) t[i] = s[i] + dt * k[i];
    }
  }

  return function step(tendency, state, dt) {
    const [k1, k2, k3, k4] = stages;
    tendency(state, k1);
    advance(state, k1, dt / 2);
    tendency(trial, k2);
    advance(state, k2, dt / 2);
    tendency(trial, k3);
    advance(state, k3, dt);
    tendency(trial, k4);
    const w = dt / 6;
    for (let a = 0; a < state.length; a++) {
      const s = state[a];
      for (let i = 0; i < s.length; i++) s[i] += w * (k1[a][i] + 2 * k2[a][i] + 2 * k3[a][i] + k4[a][i]);
    }
  };
}
