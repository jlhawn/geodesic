import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean, THIN } from '../js/ocean/layered.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';

const DT = 1800, DEG = Math.PI / 180;
const mesh = buildMesh(new Grid(16)), C = mesh.nCells, E = mesh.nEdges;

/*
 * An aquaplanet's sea a few steps into its life under trades and
 * westerlies, with tokens at the outcrops and under the thermocline.
 */
function spunUp(steps = 8, options = {}) {
  const ocean = createOcean(mesh, { everySteps: 1, ...options });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => {
    const polar = Math.max(0, Math.min(1, (Math.abs(lat) / DEG - 50) / 20));
    return Math.max(FREEZING_POINT + 0.3, 273.15 + 30 * Math.cos(lat) ** 2 - 2 - 4 * polar);
  });
  const ice = new Float64Array(C), flux = new Float64Array(C), stress = new Float64Array(E);
  for (let e = 0; e < E; e++) {
    const x = mesh.xEdge[3 * e], y = mesh.xEdge[3 * e + 1], r = Math.hypot(x, y) || 1e-12;
    stress[e] = -0.1 * Math.cos(3 * mesh.latEdge[e]) * ((-y / r) * mesh.nEdge[3 * e] + (x / r) * mesh.nEdge[3 * e + 1]);
  }
  ocean.initialize(surfaceT, ice);
  for (let n = 0; n < steps; n++) ocean.advance(surfaceT, ice, flux, stress, DT);
  return ocean;
}

/*
 * The shift lands each edge's column transport on the barotropic
 * solver's, the classes thinner than THIN carrying the velocity of the
 * layer above, so a second shift finds nothing to move.
 */
test('the velocity shift leaves the column transport on the barotropic transport, tokens included', () => {
  const ocean = spunUp();
  ocean.phases.dynamics(DT); ocean.phases.rescale(); ocean.phases.eddy(DT); ocean.phases.shift();
  const once = Float64Array.from(ocean.u);
  ocean.phases.shift();
  const L = ocean.layers, hEdge = new Float64Array(ocean.shared.hEdge);
  let thinEdges = 0, worst = 0;
  for (let e = 0; e < E; e++) {
    if (!ocean.edgeOcean[e]) continue;
    for (let k = 1; k < L; k++) if (hEdge[k * E + e] < THIN) { thinEdges++; break; }
    for (let k = 0; k < L; k++) worst = Math.max(worst, Math.abs(ocean.u[k * E + e] - once[k * E + e]));
  }
  assert.ok(thinEdges > 50, `only ${thinEdges} edges hold a class thinner than THIN`);
  assert.ok(worst < 1e-10, `a second shift moved a velocity by ${worst.toExponential(2)} m/s`);
});

/*
 * Each edge's transport as the continuity, the barotropic solver and the
 * shift count it: the top layer's centered thickness, the smaller of the
 * two below, the column cut at the edge's sill.
 */
function edgeTransports(ocean) {
  const L = ocean.layers, { h, u, eta, sillEdge } = ocean, out = new Float64Array(E);
  for (let e = 0; e < E; e++) {
    if (!ocean.edgeOcean[e]) continue;
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1], sill = Math.max(1e-9, sillEdge[e] + 0.5 * (eta[a] + eta[b]));
    let above = 0;
    for (let k = 0; k < L; k++) {
      const hk = k === 0 ? 0.5 * (h[a] + h[b]) : Math.min(h[k * C + a], h[k * C + b]);
      out[e] += Math.min(hk, Math.max(0, sill - above)) * u[k * E + e]; above += hk;
    }
  }
  return out;
}
const largestChange = (from, to) => { let worst = 0, scale = 0; for (let e = 0; e < E; e++) { worst = Math.max(worst, Math.abs(to[e] - from[e])); scale = Math.max(scale, Math.abs(from[e])); } return worst / scale; };

test('with edgeTransport kept the regrid and the vertical viscosity leave every edge\'s transport where the shift put it; remapped, they move it', () => {
  const changes = {};
  for (const edgeTransport of ['kept', 'remapped']) {
    const ocean = spunUp(8, { edgeTransport });
    ocean.phases.dynamics(DT); ocean.phases.rescale(); ocean.phases.eddy(DT); ocean.phases.shift();
    const shifted = edgeTransports(ocean);
    ocean.phases.regrid(DT);
    const regridded = edgeTransports(ocean);
    ocean.phases.column(DT);
    changes[edgeTransport] = { regrid: largestChange(shifted, regridded), column: largestChange(regridded, edgeTransports(ocean)) };
  }
  for (const phase of ['regrid', 'column']) {
    assert.ok(changes.kept[phase] < 1e-12, `kept: the ${phase} moved an edge's transport by ${changes.kept[phase].toExponential(2)} of the largest`);
    assert.ok(changes.remapped[phase] > 1e-7, `remapped: the ${phase} moved no edge's transport (${changes.remapped[phase].toExponential(2)})`);
  }
});
