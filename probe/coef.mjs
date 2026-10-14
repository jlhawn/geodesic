import { sigmaInterfaces } from '/Users/jlhawn/git_repos/jlhawn/geodesic/js/dynamics/sigmaCore.module.js';
globalThis.__coefHook = (on) => {
  const { K: nK, C: nC, kTop } = on; const dt = 900, GRAV = 9.806;
  const iface = sigmaInterfaces('bl34'); const ds = Array.from({ length: nK }, (_, k) => iface[k + 1] - iface[k]);
  for (const i of [110, 114]) {
    const pi = on.before[0][i]; let worstU = 0, worstL = 0;
    const rows = [];
    for (let k = kTop; k < nK; k++) {
      const j = k - kTop, mass = pi * ds[k] / GRAV;
      const expU = j > 0 ? dt * on.ph.MIX[(k - 1) * nC + i] / mass : 0, expL = k < nK - 1 ? dt * on.ph.MIX[k * nC + i] / mass : 0;
      const got = on.after[5][k * nC + i]; const gotL = Math.round(got / 1000 * 1e9) / 1e9; // crude split: upper is < 1000·lower scale? print raw instead
      rows.push(`k${k} raw ${got.toExponential(4)} expected upper ${expU.toExponential(4)} lower ${expL.toExponential(4)} (upper+1000·lower = ${(expU + 1000 * expL).toExponential(4)})`);
    }
    console.log(`column ${i} (pi ${pi.toFixed(0)}):\n  ` + rows.join('\n  '));
  }
};
