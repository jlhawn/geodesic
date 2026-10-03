import { CLOUD_TOP_DEFAULTS } from '../../js/physics/boundaryLayer.module.js';

/*
 * The discrete decisions of the boundary layer, the plume and the
 * condensation, read from either engine after a step, and the columns where
 * two runs took one of them apart: the regime, the plume firing, its base
 * flux by `flux` of itself (near its onset), its top by a pascal, a layer's
 * cloud water against the cloud-top threshold, and a merge of the dry
 * adjustment (two adjacent layers' θ equal to 1e-6 of itself). The same
 * rule serves the CPU against the GPU in the parity tests and the CPU
 * against itself under θ noise in scripts/perturbedRain.mjs.
 */
export const cpuDecisions = (model) => ({ regime: model.boundaryLayer.regime, flux: model.moist.cumulusBaseFlux, top: model.moist.cumulusTop, theta: model.state[1], qc: model.state[5] });
export const gpuDecisions = (state, physics) => ({ regime: physics.REGIME, flux: physics.CUMF, top: physics.CUTOP, theta: state[1], qc: state[5] });

export const DECISION_KINDS = ['regime', 'fires', 'flux', 'top', 'cloud', 'merge'];

export function decisionTracker(K, C, { flux: fluxShare = 1e-2, threshold = CLOUD_TOP_DEFAULTS.threshold } = {}) {
  const kinds = Object.fromEntries(DECISION_KINDS.map((kind) => [kind, new Set()])), parted = new Set(), first = new Map();
  let step = 0;
  const note = (kind, i, what) => { kinds[kind].add(i); parted.add(i); if (!first.has(i)) first.set(i, `step ${step}: ${what}`); };
  const merges = (theta, x) => Math.abs(theta[x] - theta[x + C]) <= 1e-6 * theta[x + C];
  return {
    kinds, parted, first,
    check(a, b) {
      step++;
      for (let i = 0; i < C; i++) {
        if (a.regime[i] !== b.regime[i]) note('regime', i, `regime ${a.regime[i]}/${b.regime[i]}`);
        const fa = a.flux[i], fb = b.flux[i];
        if ((fa > 0) !== (fb > 0)) note('fires', i, `the plume fires ${fa > 0}/${fb > 0}`);
        else if (fa > 0 && Math.abs(fa - fb) > fluxShare * Math.max(fa, fb)) note('flux', i, `plume base flux ${fa.toExponential(3)}/${fb.toExponential(3)} kg/m²/s`);
        else if (fa > 0 && Math.abs(a.top[i] - b.top[i]) > 1) note('top', i, `plume top ${(a.top[i] / 100).toFixed(0)}/${(b.top[i] / 100).toFixed(0)} hPa`);
        for (let k = 0; k < K; k++) { const x = k * C + i; if ((a.qc[x] > threshold) !== (b.qc[x] > threshold)) { note('cloud', i, `cloud in layer ${k} ${a.qc[x].toExponential(2)}/${b.qc[x].toExponential(2)}`); break; } }
        for (let k = 0; k < K - 1; k++) { const x = k * C + i; if (merges(a.theta, x) !== merges(b.theta, x)) { note('merge', i, `dry adjustment merges layer ${k}`); break; } }
      }
    },
  };
}

export function neighbourhood(mesh, cells, rings = 2) {
  const { cellsOnCell, nEdgesOnCell, maxEdges } = mesh, near = new Set(cells);
  let edge = [...cells];
  for (let r = 0; r < rings; r++) {
    const next = [];
    for (const i of edge) for (let j = 0; j < nEdgesOnCell[i]; j++) { const n = cellsOnCell[i * maxEdges + j]; if (!near.has(n)) { near.add(n); next.push(n); } }
    edge = next;
  }
  return near;
}
