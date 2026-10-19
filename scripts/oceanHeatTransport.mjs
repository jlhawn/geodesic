// The resolved meridional ocean heat transport of saved states across
// latitude circles (PW, northward +): over the edges whose two cells
// straddle each latitude, Σ sign · dv · h_e · u · (T_e − T̄) ρ c_p, the
// edge column cut at the sill as the model's edgeColumn cuts it and T̄
// the section's thickness-weighted mean temperature, so a net volume
// flux across the section adds nothing. The global and the Atlantic
// sections of each latitude, one table per state and, with several
// states, their mean with the standard deviation across them beside the
// observed values (Trenberth and Caron 2001, Ganachaud and Wunsch 2000,
// RAPID at 26.5N). A single day's value carries the day's Ekman and eddy
// cells, several petawatts at the equator; read the mean over a year's
// saves. SILLS=1 applies the geography's edge sills, for states run
// with them.
//   [SILLS=1] [QUIET=1] node scripts/oceanHeatTransport.mjs STATE.bin [STATE.bin ...]
import { readFileSync } from 'node:fs';
import { decodeState } from '../js/stateFile.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { createGeography, topographyFromInt16 } from '../js/geography.module.js';
import { edgeSills } from '../js/ocean/layered.module.js';
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const KELVIN = 273.15, deg = 180 / Math.PI, RHO_CP = 1025 * 3990;
const LATS = [-60, -50, -40, -30, -20, -10, 0, 10, 20, 26.5, 30, 40, 50, 60];
const rows = new Map(); const push = (la, g, atl) => { if (!rows.has(la)) rows.set(la, []); rows.get(la).push([g.V, g.H, atl ? atl.H : NaN]); };
const f = (x, d = 2, w = 7) => (Number.isFinite(x) ? x.toFixed(d) : '—').padStart(w);
for (const path of process.argv.slice(2)) {
  const st = await decodeState(new Uint8Array(readFileSync(path)));
  const hybrid = st.ocean.vertical === 'hybrid', mesh = stateGrid(st, topography).mesh, C = mesh.nCells, E = mesh.nEdges;
  const geo = createGeography(mesh, topography), land = geo.land, oc = st.ocean, L = oc.densities.length, { h, u, T } = oc, at = (k, i) => k * C + i;
  const sea = new Uint8Array(C); for (let i = 0; i < C; i++) { let col = 0; for (let k = 0; k < L; k++) col += h[at(k, i)]; sea[i] = !land[i] && col > 1 ? 1 : 0; }
  const lat = (i) => mesh.latCell[i] * deg, lon = (i) => mesh.lonCell[i] * deg, { cellsOnEdge, dvEdge } = mesh;
  const D = new Float64Array(C); for (let i = 0; i < C; i++) { let col = 0; for (let k = 0; k < L; k++) col += h[at(k, i)]; D[i] = col - (oc.eta ? oc.eta[i] : 0); }
  const edgeOcean = new Uint8Array(E); for (let e = 0; e < E; e++) edgeOcean[e] = sea[cellsOnEdge[2 * e]] && sea[cellsOnEdge[2 * e + 1]] ? 1 : 0;
  const sillEdge = hybrid ? edgeSills(mesh, D, edgeOcean, process.env.SILLS ? geo.edgeSill : null, 50) : null;
  const column = new Float64Array(L);
  const edgeColumn = (e, a, b) => {
    column[0] = 0.5 * (h[a] + h[b]); for (let k = 1; k < L; k++) column[k] = Math.min(h[at(k, a)], h[at(k, b)]);
    if (!hybrid) return column;
    const sill = Math.max(1e-6, sillEdge[e] + 0.5 * ((oc.eta ? oc.eta[a] : 0) + (oc.eta ? oc.eta[b] : 0)));
    let above = 0; for (let k = 0; k < L; k++) { const hk = column[k]; column[k] = Math.min(hk, Math.max(0, sill - above)); above += hk; }
    return column;
  };
  const inLon = (l, a, b) => (a <= b ? l >= a && l <= b : l >= a || l <= b);
  const section = (lat0, west, east) => {
    let V = 0, H = 0, A = 0, TA = 0;
    for (let e = 0; e < E; e++) {
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1]; if (!sea[a] || !sea[b]) continue;
      if (!inLon(0.5 * (lon(a) + lon(b)), west, east)) continue;
      const la = lat(a), lb = lat(b); if ((la < lat0) === (lb < lat0)) continue;
      const sign = lb > la ? 1 : -1, col = edgeColumn(e, a, b);
      for (let k = 0; k < L; k++) {
        const he = col[k]; if (he <= 1e-6) continue;
        const flux = sign * dvEdge[e] * he * u[k * E + e], Te = 0.5 * (T[at(k, a)] + T[at(k, b)]) - KELVIN;
        V += flux; H += flux * Te; A += dvEdge[e] * he; TA += dvEdge[e] * he * Te;
      }
    }
    const Tbar = TA / A; return { V: V / 1e6, H0: RHO_CP * H / 1e15, H: RHO_CP * (H - Tbar * V) / 1e15, Tbar };
  };
  console.log(`\n=== ${path.split('/').pop()}: ${hybrid ? 'hybrid' : 'slab'} day ${st.day}`);
  console.log('  lat      net V (Sv)   heat (PW, T rel. 0 °C)   heat (PW, rel. section mean T)   Atlantic heat (PW)');
  for (const la of LATS) {
    const g = section(la, -180, 180), atl = la > -35 ? section(la, la > 10 ? -100 : -70, la > 10 ? 0 : 20) : null;
    push(la, g, atl);
    if (process.env.QUIET) continue;
    console.log(`  ${String(la).padStart(5)}N ${f(g.V, 1, 9)} ${f(g.H0, 2, 12)} ${f(g.H, 2, 22)} ${atl ? f(atl.H, 2, 20) : ''}`);
  }
}

const stat = (xs) => { const ok = xs.filter(Number.isFinite), m = ok.reduce((a, b) => a + b, 0) / ok.length; return [m, Math.sqrt(ok.reduce((a, b) => a + (b - m) ** 2, 0) / ok.length)]; };
if (process.argv.length > 3) {
console.log(`\n=== mean over ${process.argv.length - 2} states (± standard deviation across them)`);
console.log('  lat      net V (Sv)      heat (PW, rel. section mean T)   Atlantic heat (PW)   observed global / Atlantic (PW)');
const OBS = { '-60': '−0.2', '-50': '−0.4', '-40': '−0.8', '-30': '−1.0 / +0.3', '-20': '−1.3 / +0.4', '-10': '−1.0 / +0.6', '0': '+0.2 / +0.8', '10': '+1.3 / +1.0', '20': '+1.7 / +1.2', '26.5': '+1.6 / +1.25', '30': '+1.4 / +1.1', '40': '+0.9 / +0.7', '50': '+0.5 / +0.5', '60': '+0.2 / +0.3' };
for (const [la, xs] of rows) { const [v, sv] = stat(xs.map((r) => r[0])), [hh, sh] = stat(xs.map((r) => r[1])), [a, sa] = stat(xs.map((r) => r[2])); console.log(`  ${String(la).padStart(5)}N ${f(v, 1, 7)} ± ${f(sv, 1, 4)} ${f(hh, 2, 12)} ± ${f(sh, 2, 4)} ${Number.isFinite(a) ? `${f(a, 2, 16)} ± ${f(sa, 2, 4)}` : ''.padStart(23)}   ${OBS[String(la)]}`); }
}
