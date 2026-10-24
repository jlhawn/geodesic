import { Grid } from '../../js/grid.module.js';
import { buildMesh } from '../../js/mesh.module.js';
import { beyondCells, limitedEdgeValue } from '../../js/ocean/layered.module.js';
import { createIceDrift } from '../../js/ocean/iceDrift.module.js';
import { createSeaIce, FREEZING_POINT } from '../../js/physics/ice.module.js';
import { createGeography, syntheticTopography } from '../../js/geography.module.js';

export const DEG = Math.PI / 180, RHO_W = 1025, RHO_I = 917;

const meshes = new Map();
export const meshOf = (N) => { if (!meshes.has(N)) meshes.set(N, buildMesh(new Grid(N))); return meshes.get(N); };

export const unit = (a, i) => [a[3 * i], a[3 * i + 1], a[3 * i + 2]];
export const dot = (p, q) => p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
export const cross = (p, q) => [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]];

export function nearestCell(m, lat, lon) {
  const target = [Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat)];
  let best = -2, at = 0;
  for (let i = 0; i < m.nCells; i++) { const d = dot(unit(m.xCell, i), target); if (d > best) { best = d; at = i; } }
  return at;
}
export function localFrame(m, i) {
  const x = unit(m.xCell, i), r = Math.hypot(x[0], x[1]), east = [-x[1] / r, x[0] / r, 0];
  return { east, north: cross(x, east) };
}
export const windOn = (m, v) => Float64Array.from({ length: m.nEdges }, (_, e) => dot(unit(m.nEdge, e), v));
export function rotationOn(m, omega) {
  return Float64Array.from({ length: m.nEdges }, (_, e) => { const x = unit(m.xEdge, e), n = unit(m.nEdge, e); return omega * m.radius * (-x[1] * n[0] + x[0] * n[1]); });
}
/*
 * The same solid-body rotation from a streamfunction at the vertices, whose
 * discrete divergence is zero to roundoff on any mesh.
 */
export function streamOn(m, omega) {
  const psi = (v) => -omega * m.radius * m.radius * m.xVertex[3 * v + 2];
  const u = Float64Array.from({ length: m.nEdges }, (_, e) => (psi(m.verticesOnEdge[2 * e + 1]) - psi(m.verticesOnEdge[2 * e])) / m.dvEdge[e]);
  const solid = rotationOn(m, omega);
  let agree = 0;
  for (let e = 0; e < m.nEdges; e++) agree += u[e] * solid[e];
  if (agree < 0) for (let e = 0; e < m.nEdges; e++) u[e] = -u[e];
  return u;
}
/*
 * The same rotation about any axis, from a streamfunction at the vertices.
 */
export function streamAbout(m, axis, omega) {
  const norm = Math.hypot(...axis), a = axis.map((v) => v / norm);
  const psi = (v) => -omega * m.radius * m.radius * (a[0] * m.xVertex[3 * v] + a[1] * m.xVertex[3 * v + 1] + a[2] * m.xVertex[3 * v + 2]);
  const u = Float64Array.from({ length: m.nEdges }, (_, e) => (psi(m.verticesOnEdge[2 * e + 1]) - psi(m.verticesOnEdge[2 * e])) / m.dvEdge[e]);
  const solid = Float64Array.from({ length: m.nEdges }, (_, e) => {
    const x = unit(m.xEdge, e), n = unit(m.nEdge, e), v = cross(a, x);
    return omega * m.radius * dot(v, n);
  });
  let agree = 0;
  for (let e = 0; e < m.nEdges; e++) agree += u[e] * solid[e];
  if (agree < 0) for (let e = 0; e < m.nEdges; e++) u[e] = -u[e];
  return u;
}
/*
 * The component of the local northward unit vector along each edge's
 * normal, where `keep(latitude)` and zero elsewhere: times a speed it is a
 * drift toward the equator in the Southern Hemisphere.
 */
export function northwardOn(m, keep = () => true) {
  return Float64Array.from({ length: m.nEdges }, (_, e) => {
    const x = unit(m.xEdge, e), r = Math.hypot(x[0], x[1]);
    if (!keep(m.latEdge[e]) || r === 0) return 0;
    return dot(unit(m.nEdge, e), [-x[2] * x[0] / r, -x[2] * x[1] / r, r]);
  });
}
export function edgesNear(m, c0, spacings) {
  const x0 = unit(m.xCell, c0), spacing = m.dcEdge.reduce((s, d) => s + d, 0) / m.nEdges, out = [];
  for (let e = 0; e < m.nEdges; e++) if (Math.acos(Math.min(1, dot(unit(m.xEdge, e), x0))) * m.radius < spacings * spacing) out.push(e);
  return out;
}

/*
 * The steady solution w = u − U_w of 0 = F − c|w| e^{iθe} w − i μ f w, F the
 * force on the ice other than the water drag and the Coriolis force, in the
 * edge's (n, t) frame with t = k × n; θe = sign(f) θ turns the water stress
 * counterclockwise in the Northern Hemisphere.
 */
export function closedForm(Fn, Ft, f, mu, theta, c = RHO_W * 5.5e-3) {
  const te = Math.sign(f) * theta * DEG, ce = Math.cos(te), se = Math.sin(te), F2 = Fn * Fn + Ft * Ft;
  const g = (R) => R * R * ((c * R * ce) ** 2 + (c * R * se + mu * f) ** 2) - F2;
  let lo = 0, hi = 50;
  for (let k = 0; k < 200; k++) { const mid = 0.5 * (lo + hi); if (g(mid) > 0) hi = mid; else lo = mid; }
  const R = 0.5 * (lo + hi), dr = c * R * ce, di = c * R * se + mu * f, den = dr * dr + di * di;
  return { n: (Fn * dr + Ft * di) / den, t: (Ft * dr - Fn * di) / den };
}

/*
 * The drift module on its own: a sea-ice object of the thermodynamics, open
 * cells everywhere or where a geography has sea, and the arrays the
 * ocean would hand it.
 */
export function rig(m, options, geography = null, extra = {}, iceOptions = {}) {
  const C = m.nCells, E = m.nEdges;
  const cellOpen = geography ? Uint8Array.from(geography.land, (l) => (l ? 0 : 1)) : new Uint8Array(C).fill(1);
  const edgeOcean = geography ? geography.edgeOcean : new Uint8Array(E).fill(1);
  const previousVolume = new Float64Array(C), seaIce = createSeaIce(m, iceOptions);
  const drift = createIceDrift(m, options, { cellOpen, edgeOcean, waterDensity: RHO_W, iceDensity: RHO_I, gravity: 9.81, previousVolume, beyond: beyondCells(m), limited: limitedEdgeValue, ...extra });
  drift.attach(seaIce);
  const r = {
    m, C, E, drift, seaIce, cellOpen, edgeOcean, previousVolume, dt: 1350,
    surfaceT: new Float64Array(C).fill(FREEZING_POINT), thickness: new Float64Array(C), concentration: seaIce.concentration, snow: seaIce.snow,
    stress: new Float64Array(E), waterU: new Float64Array(E), eta: new Float64Array(C),
  };
  r.env = () => ({ surfaceT: r.surfaceT, ice: r.thickness, concentration: r.concentration, stress: r.stress, waterU: r.waterU, eta: r.eta, dt: r.dt });
  return r;
}

/*
 * The compact pack: sea poleward of lat0 closed by land, ice everywhere in
 * it, a uniform wind along the x axis.
 */
export function packRig(N, options, { lat0 = 70, h = 1, wind = 0.1, dt = 5400 } = {}) {
  const m = meshOf(N), geography = createGeography(m, syntheticTopography(180, 360, (lat) => (lat > lat0 * DEG ? -4000 : 500)));
  const r = rig(m, options, geography);
  r.dt = dt;
  for (let i = 0; i < r.C; i++) if (r.cellOpen[i]) { r.thickness[i] = h; r.concentration[i] = 1; }
  r.stress.set(windOn(m, [wind, 0, 0]));
  return r;
}
export function packSolve(N, dt, options, calls = 1) {
  const r = packRig(N, options, { dt });
  for (let n = 0; n < calls; n++) r.drift.advance(r.env());
  return { velocity: Float64Array.from(r.drift.velocity), valid: Uint8Array.from(r.drift.valid), volume: Float64Array.from(r.thickness, (h, i) => h * r.concentration[i]), r };
}
export function packDifference(a, b) {
  let num = 0, den = 0;
  for (let e = 0; e < a.velocity.length; e++) if (b.valid[e] || a.valid[e]) { num += (a.velocity[e] - b.velocity[e]) ** 2; den += b.velocity[e] ** 2; }
  return Math.sqrt(num / den);
}
