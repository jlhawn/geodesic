import { CLOUD_TOP_DEFAULTS } from '../../js/physics/boundaryLayer.module.js';

/*
 * The discrete decisions of the boundary layer, the plume and the
 * condensation, read from either engine after a step, and the columns where
 * two runs took one of them apart: the regime, the plume firing, its base
 * flux by `flux` of itself (near its onset), its top by a pascal, a layer's
 * cloud water against the cloud-top threshold, and a merge of the dry
 * adjustment (two adjacent layers' θ and q each equal to 1e-6 of
 * themselves: the adjustment writes both, and only θ changes after it). The same
 * rule serves the CPU against the GPU in the parity tests and the CPU
 * against itself under θ noise in scripts/perturbedRain.mjs. A run watched
 * by radiationTraces also notes the plume layers whose cloud water, as the
 * radiation is called, lies within `share` of the largest around it: any
 * cloud water lifts a layer's cover to the floor, and a remainder that small
 * is the transport's round-off, whose sign the engines need not share.
 */
export const cpuDecisions = (model, traces = null) => ({ regime: model.boundaryLayer.regime, flux: model.moist.cumulusBaseFlux, top: model.moist.cumulusTop, theta: model.state[1], q: model.state[4], qc: model.state[5], traces });
export const gpuDecisions = (state, physics) => ({ regime: physics.REGIME, flux: physics.CUMF, top: physics.CUTOP, theta: state[1], q: state[4], qc: state[5] });

export const DECISION_KINDS = ['regime', 'fires', 'flux', 'top', 'cloud', 'merge', 'trace'];

export function radiationTraces(model, share = 1e-5) {
  const { mesh: { cellsOnCell, nEdgesOnCell, maxEdges, nCells: C }, core: { K }, radiation, moist: { cumulusCover, cumulusWater } } = model, apply = radiation.apply, traces = new Map();
  radiation.apply = (state, out, windSpeed, totals, iFrom = 0, iTo = C, ...rest) => {
    const qc = state[5];
    if (qc) for (let i = iFrom; i < iTo; i++) for (let k = 0; k < K; k++) {
      const x = k * C + i;
      if (qc[x] === 0 || !(cumulusCover[x] * cumulusWater[x] > 0)) continue;
      let largest = Math.max(k > 0 ? Math.abs(qc[x - C]) : 0, k < K - 1 ? Math.abs(qc[x + C]) : 0);
      for (let j = 0; j < nEdgesOnCell[i]; j++) largest = Math.max(largest, Math.abs(qc[k * C + cellsOnCell[i * maxEdges + j]]));
      if (Math.abs(qc[x]) <= share * largest) { traces.set(i, `plume layer ${k} holds ${qc[x].toExponential(2)} kg/kg of cloud water as the radiation is called`); break; }
    }
    return apply(state, out, windSpeed, totals, iFrom, iTo, ...rest);
  };
  return traces;
}

export function decisionTracker(K, C, { flux: fluxShare = 1e-2, threshold = CLOUD_TOP_DEFAULTS.threshold } = {}) {
  const kinds = Object.fromEntries(DECISION_KINDS.map((kind) => [kind, new Set()])), parted = new Set(), first = new Map();
  let step = 0;
  const note = (kind, i, what) => { kinds[kind].add(i); parted.add(i); if (!first.has(i)) first.set(i, `step ${step}: ${what}`); };
  const merges = (theta, q, x) => Math.abs(theta[x] - theta[x + C]) <= 1e-6 * theta[x + C] && (!q || Math.abs(q[x] - q[x + C]) <= 1e-6 * Math.max(Math.abs(q[x]), Math.abs(q[x + C])));
  return {
    kinds, parted, first,
    check(a, b) {
      step++;
      for (const run of [a, b]) if (run.traces) { for (const [i, what] of run.traces) note('trace', i, what); run.traces.clear(); }
      for (let i = 0; i < C; i++) {
        if (a.regime[i] !== b.regime[i]) note('regime', i, `regime ${a.regime[i]}/${b.regime[i]}`);
        const fa = a.flux[i], fb = b.flux[i];
        if ((fa > 0) !== (fb > 0)) note('fires', i, `the plume fires ${fa > 0}/${fb > 0}`);
        else if (fa > 0 && Math.abs(fa - fb) > fluxShare * Math.max(fa, fb)) note('flux', i, `plume base flux ${fa.toExponential(3)}/${fb.toExponential(3)} kg/m²/s`);
        else if (fa > 0 && Math.abs(a.top[i] - b.top[i]) > 1) note('top', i, `plume top ${(a.top[i] / 100).toFixed(0)}/${(b.top[i] / 100).toFixed(0)} hPa`);
        for (let k = 0; k < K; k++) { const x = k * C + i; if ((a.qc[x] > threshold) !== (b.qc[x] > threshold)) { note('cloud', i, `cloud in layer ${k} ${a.qc[x].toExponential(2)}/${b.qc[x].toExponential(2)}`); break; } }
        for (let k = 0; k < K - 1; k++) { const x = k * C + i; if (merges(a.theta, a.q, x) !== merges(b.theta, b.q, x)) { note('merge', i, `dry adjustment merges layer ${k}`); break; } }
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
