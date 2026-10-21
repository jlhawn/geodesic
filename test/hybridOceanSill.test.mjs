import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createOcean } from '../js/ocean/layered.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';
import { HYBRID_LAYERS, H_MIN } from '../js/ocean/column.module.js';
import { DT, DEG, FLAT, mesh, C, E, maxSpeed } from './helpers/layered.mjs';

test('a class held at its floor on one side of an edge, with no depth in common with the layer across it, follows the layer above while the front between two stratifications adjusts', () => {
  const deep = (i) => Math.cos(mesh.lonCell[i]) > 0;
  const bathymetry = Float64Array.from({ length: C }, (_, i) => (deep(i) ? 4000 : 1500));
  const ocean = createOcean(mesh, { ...FLAT, bathymetry });
  const surfaceT = Float64Array.from({ length: C }, (_, i) => (deep(i) ? 300 : 280)), ice = new Float64Array(C), flux = new Float64Array(C), calm = new Float64Array(E);
  ocean.initialize(surfaceT, ice);
  const L = ocean.layers, { cellsOnEdge } = mesh, tokenEdge = Array.from(HYBRID_LAYERS.minimum, (m) => 2 * Math.min(m, HYBRID_LAYERS.spacing));
  const tokens = () => {
    const Ha = new Float64Array(E), Hb = new Float64Array(E), list = [];
    for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) {
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], ha = ocean.h[k * C + a], hb = ocean.h[k * C + b], he = Math.min(ha, hb);
      if (k >= HYBRID_LAYERS.fixed && ocean.edgeOcean[e] && he > 2 * H_MIN && he <= tokenEdge[k] && Math.min(Ha[e] + ha, Hb[e] + hb) <= Math.max(Ha[e], Hb[e])) list.push([k, e]);
      Ha[e] += ha; Hb[e] += hb;
    }
    return list;
  };
  let shear = 0, count = Infinity;
  for (let n = 0; n < 64; n++) {
    ocean.advance(surfaceT, ice, flux, calm, DT);
    const list = tokens();
    count = Math.min(count, list.length);
    for (const [k, e] of list) shear = Math.max(shear, Math.abs(ocean.u[k * E + e] - ocean.u[(k - 1) * E + e]));
  }
  const fastest = maxSpeed(ocean);
  assert.ok(count > 1000, `${count} token edge-layers`);
  assert.ok(fastest > 0.5, `the front adjusts: ${fastest} m/s`);
  assert.ok(shear < 0.012 * fastest, `a token's shear to the layer above reached ${shear} m/s against a fastest current of ${fastest} m/s`);
  assert.equal(ocean.diagnostics().oceanLimited, 0);
});

test('an edge takes the deepest channel its two cells hold between them as its sill and the ocean cuts the edge column there so the classes below carry nothing across', () => {
  const { cellsOnEdge, lonCell, latCell } = mesh;
  const wrap = (x) => Math.atan2(Math.sin(x), Math.cos(x));
  const ridge = (lon) => Math.abs(lon) < 1.5 * DEG, strip = (lon) => Math.abs(wrap(lon - Math.PI)) < 1.5 * DEG, within = (lat) => Math.abs(lat) < 40 * DEG;
  const topography = syntheticTopography(180, 360, (lat, lon) => (within(lat) && ridge(lon) ? -300 : within(lat) && strip(lon) ? 100 : -4000));
  const geography = createGeography(mesh, topography);
  const spacing = mesh.dcEdge.reduce((sum, d) => sum + d, 0) / E / mesh.radius, near = 0.38 * spacing, far = 1.74 * spacing;
  const across = (e, lon0) => { const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], la = wrap(lonCell[a] - lon0), lb = wrap(lonCell[b] - lon0); return geography.edgeOcean[e] && la * lb < 0 && Math.abs(la) > near && Math.abs(lb) > near && Math.abs(la) < far && Math.abs(lb) < far && Math.abs(latCell[a]) < 30 * DEG && Math.abs(latCell[b]) < 30 * DEG; };
  const ridgeEdges = [], stripEdges = [], openEdges = [];
  for (let e = 0; e < E; e++) {
    if (across(e, 0)) ridgeEdges.push(e);
    else if (across(e, Math.PI)) stripEdges.push(e);
    else if (geography.edgeOcean[e] && Math.abs(lonCell[cellsOnEdge[2 * e]]) > 20 * DEG && Math.abs(lonCell[cellsOnEdge[2 * e]]) < 70 * DEG) openEdges.push(e);
  }
  assert.ok(ridgeEdges.length >= 8 && stripEdges.length >= 8 && openEdges.length > 0.1 * E, `${ridgeEdges.length} ridge, ${stripEdges.length} strip, ${openEdges.length} open edges`);
  for (const e of ridgeEdges) assert.equal(geography.edgeSill[e], 300, `ridge edge ${e}`);
  for (const e of stripEdges) assert.equal(geography.edgeSill[e], 0, `land-strip edge ${e}`);
  for (const e of openEdges) assert.equal(geography.edgeSill[e], 4000, `open edge ${e}`);
  const hybrid = createOcean(mesh, { ...FLAT, geography });
  const floor = (e) => Math.min(hybrid.D[cellsOnEdge[2 * e]], hybrid.D[cellsOnEdge[2 * e + 1]]);
  for (const e of ridgeEdges) { assert.equal(hybrid.sillEdge[e], 300); assert.ok(floor(e) > 1000, `the cells stay deep at ${e}: ${floor(e)}`); }
  for (const e of stripEdges) assert.equal(hybrid.sillEdge[e], 50);
  for (const e of openEdges) assert.equal(hybrid.sillEdge[e], floor(e));
  const front = Float64Array.from({ length: C }, (_, i) => (lonCell[i] < 0 ? 280 : 300)), ice = new Float64Array(C), flux = new Float64Array(C), calm = new Float64Array(E);
  const L = hybrid.layers;
  const transport = (ocean) => { let through = 0, below = 0; for (const e of ridgeEdges) { let z = 0; for (let k = 0; k < L; k++) { const t = Math.abs(ocean.hEdge[k * E + e] * ocean.u[k * E + e]) * mesh.dvEdge[e]; through += t; if (z >= 300) below += t; z += ocean.hEdge[k * E + e]; } } return { through, below }; };
  const open = createOcean(mesh, { ...FLAT, geography, edgeSill: null });
  hybrid.initialize(front, ice); open.initialize(front, ice);
  let sillThrough = 0, sillBelow = 0, openThrough = 0;
  for (let n = 0; n < 48; n++) {
    hybrid.advance(front, ice, flux, calm, DT); open.advance(front, ice, flux, calm, DT);
    const a = transport(hybrid), b = transport(open);
    sillThrough += a.through; sillBelow += a.below; openThrough += b.through;
  }
  for (const e of ridgeEdges) {
    let z = 0, sum = 0;
    for (let k = 0; k < L; k++) { const he = hybrid.hEdge[k * E + e]; if (z >= 300) assert.equal(he, 0, `layer ${k} of ridge edge ${e} at ${z} m`); z += Math.min(hybrid.h[k * C + cellsOnEdge[2 * e]], hybrid.h[k * C + cellsOnEdge[2 * e + 1]]); sum += he; }
    assert.ok(Math.abs(sum - 300) < 1, `ridge edge ${e} column ${sum} m`);
  }
  assert.ok(maxSpeed(hybrid) > 0.1 && maxSpeed(open) > 0.1, 'the front adjusts');
  assert.equal(sillBelow, 0);
  assert.ok(openThrough > 0 && sillThrough < 0.5 * openThrough, `transport through the ridge edges over 48 steps: ${(sillThrough / 48 / 1e6).toFixed(1)} Sv with the sill, ${(openThrough / 48 / 1e6).toFixed(1)} Sv without`);
  assert.equal(hybrid.diagnostics().oceanLimited, 0);
});
