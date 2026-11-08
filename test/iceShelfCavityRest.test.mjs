import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createGeography, ICE_DENSITY } from '../js/geography.module.js';
import { basalMelt, baseFreezing, SHELF_MELT } from '../js/ocean/iceShelfMelt.module.js';
import { N, shelfBands, warmOcean, restingShelf, shelfModel, prepareShelfModel, shelfTops } from './helpers/iceShelf.mjs';
import { seawaterDensity } from '../js/ocean/seawater.module.js';

test(`on ISOMIP+'s warm profile the shelf melts, fastest near its grounding line, at meters of ice a year at N=${N}, and with no discharge to take it from the sheet's store pays for all of it`, () => {
  const { mesh, raster, options } = shelfBands(250), geography = createGeography(mesh, raster, { ...options, iceShelves: 'cavity' }), C = mesh.nCells, E = mesh.nEdges, days = 2;
  const run = (tidalSpeed) => {
    const { ocean, surfaceT, ice } = warmOcean(mesh, geography, { tidalSpeed }), flux = new Float64Array(C), calm = new Float64Array(E);
    for (let n = 0; n < days * 64; n++) ocean.advance(surfaceT, ice, flux, calm, 1350);
    return ocean;
  };
  const tide = run(SHELF_MELT.tidalSpeed), still = run(0), year = 365.25 * 86400, iceRate = (o, i) => o.shelfMelted[i] / (days * 86400) * year / ICE_DENSITY;
  const cells = [...geography.cavity.keys()].filter((i) => geography.cavity[i]), drafts = cells.map((i) => geography.iceDraft[i]), deep = Math.max(...drafts), shallow = Math.min(...drafts);
  const mean = (o, where) => { let a = 0, m = 0; for (const i of cells) if (where(i)) { a += mesh.areaCell[i]; m += mesh.areaCell[i] * iceRate(o, i); } return m / a; };
  const near = mean(tide, (i) => geography.iceDraft[i] > deep - 1), front = mean(tide, (i) => geography.iceDraft[i] < shallow + 1), all = [mean(tide, () => true), mean(still, () => true)];
  console.log(`warm profile, N=${N}, ${days} days: ${near.toFixed(2)} m/yr of ice at the ${deep.toFixed(0)} m draft by the grounding line, ${front.toFixed(2)} at ${shallow.toFixed(0)} m by the front; mean ${all[0].toFixed(2)} with u_tidal ${SHELF_MELT.tidalSpeed} m/s, ${all[1].toFixed(2)} without`);
  assert.ok(near > front && front > 0 && near > 0.3 && near < 30, `${near} m/yr by the grounding line, ${front} by the front`);
  assert.ok(all[0] > all[1], `the tide's stirring adds melt: ${all}`);
  const { discharge, melt, drawn } = tide.shelfBudget;
  let draw = 0, melted = 0;
  for (const i of cells) { draw += mesh.areaCell[i] * tide.shelfDraw[i]; melted += mesh.areaCell[i] * tide.shelfMelted[i]; }
  assert.ok(discharge[0] === 0 && drawn[0] === melt[0] && Math.abs(draw - melted) < 1e-12 * melted, `with no discharge: ${draw} kg drawn from the store, ${melted} kg melted`);
});

test('with a 20 m melt layer, as in Losch (2008), the melt reads T_w and S_w as the mean of the 20 m below the ice base, spreads its heat and salt over them, closes both in each column, and melts more by the grounding line than the 3 m top layer does', () => {
  const { mesh, raster, options } = shelfBands(250), geography = createGeography(mesh, raster, { ...options, iceShelves: 'cavity' }), C = mesh.nCells, E = mesh.nEdges, DT = 3600;
  const cells = [...geography.cavity.keys()].filter((i) => geography.cavity[i]), { ocean } = warmOcean(mesh, geography, { shelfMeltLayer: 20 }), L = ocean.layers, load = ocean.serialize().shelfLoad;
  const k = { rho0: 1025, rhoCp: ocean.rhoCp, latentHeat: 3.34e5 };
  ocean.u.fill(0);
  const before = cells.map((i) => {
    let slab = 0, t = 0, s = 0, q = 0, w = 0;
    for (let n = 0; n < L; n++) {
      const j = n * C + i, part = Math.min(ocean.h[j], 20 - slab);
      if (part > 0) { t += part * ocean.Q[j] / ocean.h[j]; s += part * ocean.W[j] / ocean.h[j]; slab += part; }
      q += ocean.Q[j]; w += ocean.W[j];
    }
    return { i, T: t / slab, S: s / slab, q, w };
  });
  ocean.phases.column(DT);
  for (const { i, T, S, q, w } of before) {
    const { melt, Tb } = basalMelt(T, S, baseFreezing(9.81 * load[i]), Math.sqrt(2.5e-3 * SHELF_MELT.tidalSpeed ** 2), k), m = melt * DT;
    let q1 = 0, w1 = 0;
    for (let n = 0; n < L; n++) { q1 += ocean.Q[n * C + i]; w1 += ocean.W[n * C + i]; }
    assert.ok(Math.abs(ocean.melted[i] - m) < 1e-12 * Math.abs(m) && Math.abs(ocean.meltBase[i] - Tb) < 1e-12 * Tb, `cell ${i}: melt ${ocean.melted[i]} against ${m} from the 20 m means ${T} K, ${S} psu`);
    assert.ok(Math.abs(q1 - q + 3.34e5 * m / ocean.rhoCp + m / 1000 * (T - Tb)) < 1e-12 * q && Math.abs(w1 - w + S * m / 1000) < 1e-12 * w, `cell ${i}: heat ${q1 - q}, salt ${w1 - w}`);
  }
  const run = (shelfMeltLayer) => {
    const { ocean: o, surfaceT, ice } = warmOcean(mesh, geography, { shelfMeltLayer }), flux = new Float64Array(C), calm = new Float64Array(E);
    for (let n = 0; n < 2 * 64; n++) o.advance(surfaceT, ice, flux, calm, 1350);
    return o;
  };
  const top = run(0), losch = run(20), deep = Math.max(...cells.map((i) => geography.iceDraft[i]));
  const mean = (o, where) => { let a = 0, m = 0; for (const i of cells) if (where(i)) { a += mesh.areaCell[i]; m += mesh.areaCell[i] * o.shelfMelted[i]; } return m / a / (2 * 86400) * 365.25 * 86400 / ICE_DENSITY; };
  const near = [mean(top, (i) => geography.iceDraft[i] > deep - 1), mean(losch, (i) => geography.iceDraft[i] > deep - 1)], all = [mean(top, () => true), mean(losch, () => true)];
  console.log(`warm profile, N=${N}, 2 days, m/yr of ice: by the grounding line ${near.map((x) => x.toFixed(2)).join(' (3 m top layer) against ')} (20 m layer); mean ${all.map((x) => x.toFixed(2)).join(' against ')}`);
  assert.ok(near[1] > near[0] && all[1] > all[0], `the 20 m layer's melt ${near[1]}, ${all[1]} against the top layer's ${near[0]}, ${all[0]}`);
});

test('with the melt off, a homogeneous ocean under a sloping shelf with a front stays as still for a month as it does with the shelves as land, on the CPU', () => restingShelf(['cpu']));

test('a save without the shelves\' load gives the CPU each shelf\'s load from the first layer thicker than 1e-6 m under it, at that layer\'s own Q/h and W/h: a top layer of 5 mm or 2e-6 m, the layer below a top layer of 1e-7 m, 1e-6 m or none, and the rebuilt top layer of an emptied column', async () => {
  const source = await prepareShelfModel(shelfModel());
  const { C, thin, empty, dry, submicron, micron, twoMicron, Q, W, saved } = shelfTops(source.ocean, source.state[6]);
  const model = await prepareShelfModel(shelfModel()), { iceDraft } = model.geography, { ocean } = model;
  ocean.load(saved, Float64Array.from(model.state[3]), Float64Array.from(model.state[6]));
  const { shelfLoad } = ocean.serialize();
  assert.ok(ocean.h[dry] > 1e-6, `the load left the emptied column's top layer at ${ocean.h[dry]} m`);
  for (const [i, n, [q, w, h]] of [[thin, thin, [Q, W, saved.h]], [twoMicron, twoMicron, [Q, W, saved.h]], [empty, C + empty, [Q, W, saved.h]], [submicron, C + submicron, [Q, W, saved.h]], [micron, C + micron, [Q, W, saved.h]], [dry, dry, [ocean.Q, ocean.W, ocean.h]]]) {
    const weight = iceDraft[i] * seawaterDensity(q[n] / h[n], w[n] / h[n]), layer = (m) => iceDraft[i] * seawaterDensity(Q[m] / saved.h[m], W[m] / saved.h[m]);
    if ([submicron, micron, twoMicron].includes(i)) assert.ok(Math.abs(layer(i) - layer(C + i)) > 1e-9 * weight, `under the ${saved.h[i]} m top layer at ${i}, layers 0 and 1 load the shelf within 1e-9 of each other: ${layer(i)} and ${layer(C + i)}`);
    assert.ok(Math.abs(shelfLoad[i] - weight) <= 1e-12 * weight, `the shelf at ${i}, over a ${saved.h[i]} m top layer, loads at ${shelfLoad[i]}; the water under it gives ${weight}`);
  }
});
