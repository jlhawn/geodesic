import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createSeaIce, openWaterAlbedo, FREEZING_POINT, MELTING_POINT } from '../js/physics/ice.module.js';
import { initializeState } from '../js/physics/init.module.js';

const model = createModel(new Grid(2));
const seaIce = createSeaIce(model.mesh, { oceanHeatFlux: 0, oceanDiffusivity: 0 });

test('the surface energy changes by exactly the surface flux through freezing, growth, melting and thaw', () => {
  const surfaceT = new Float64Array([FREEZING_POINT + 2]), ice = new Float64Array([0]);
  const flux = new Float64Array(1);
  const dt = 900;
  let energy = seaIce.energy(surfaceT[0], ice[0]);
  const scale = seaIce.slabHeatCapacity * 2;
  let froze = false, melted = false;
  for (let n = 0; n < 4 * 96 * 60; n++) {
    flux[0] = n < 96 * 90 ? -150 : 250;
    seaIce.update(surfaceT, ice, flux, 0, dt);
    energy += dt * flux[0];
    assert.ok(Math.abs(seaIce.energy(surfaceT[0], ice[0]) - energy) < 1e-9 * scale, `step ${n}: energy ${seaIce.energy(surfaceT[0], ice[0])} vs ${energy}`);
    if (ice[0] > 0) { froze = true; assert.ok(surfaceT[0] <= MELTING_POINT + 1e-9 && surfaceT[0] < FREEZING_POINT + 20); }
    if (froze && ice[0] === 0) melted = true;
  }
  assert.ok(froze && melted, 'expected the cell to freeze over and later thaw');
  assert.ok(surfaceT[0] > FREEZING_POINT, 'open water ends above freezing');
});

test('a cold skin over thick ice conducts heat up from the base and grows the ice', () => {
  const surfaceT = new Float64Array([FREEZING_POINT - 20]), ice = new Float64Array([1]);
  const flux = new Float64Array([0]);
  seaIce.update(surfaceT, ice, flux, 0, 900);
  assert.ok(ice[0] > 1, 'ice grew');
  assert.ok(surfaceT[0] > FREEZING_POINT - 20, 'the skin warmed');
});

test('the ocean heat convergence has zero global mean and warms the poles at the tropics\' expense', () => {
  const withFlux = createSeaIce(model.mesh, { oceanHeatFlux: 25 });
  let sum = 0, area = 0, pole = 0, equator = 0;
  for (let i = 0; i < model.mesh.nCells; i++) {
    sum += model.mesh.areaCell[i] * withFlux.convergence[i]; area += model.mesh.areaCell[i];
    if (Math.abs(model.mesh.latCell[i]) > 1.3) pole = withFlux.convergence[i];
    if (Math.abs(model.mesh.latCell[i]) < 0.2) equator = withFlux.convergence[i];
  }
  assert.ok(Math.abs(sum / area) < 0.5, `global mean ${sum / area} W/m²`);
  assert.ok(pole > 30 && equator < -20);
});

test('mixed-layer diffusion conserves energy and carries heat from open water into the ice edge', () => {
  const m = createModel(new Grid(4));
  const { mesh } = m;
  const withDiffusion = createSeaIce(mesh, { oceanHeatFlux: 0, oceanDiffusivity: 0.3 });
  const surfaceT = new Float64Array(mesh.nCells), ice = new Float64Array(mesh.nCells);
  for (let i = 0; i < mesh.nCells; i++) {
    const polar = Math.abs(mesh.latCell[i]) > Math.PI / 3;
    surfaceT[i] = polar ? FREEZING_POINT - 10 : 300 - 25 * Math.sin(mesh.latCell[i]) ** 2;
    ice[i] = polar ? 1 : 0;
  }
  withDiffusion.prepare(surfaceT, ice);
  let total = 0, magnitude = 0, edgeIce = 0, interiorIce = 0, edgeWater = 0;
  for (let i = 0; i < mesh.nCells; i++) {
    const a = mesh.areaCell[i], f = withDiffusion.oceanFlux[i];
    total += a * f; magnitude += a * Math.abs(f);
    let openNeighbour = false, icyNeighbour = false;
    for (let k = 0; k < mesh.nEdgesOnCell[i]; k++) {
      const j = mesh.cellsOnCell[mesh.maxEdges * i + k];
      if (ice[j] > 0) icyNeighbour = true; else openNeighbour = true;
    }
    if (ice[i] > 0 && openNeighbour) edgeIce = Math.max(edgeIce, f);
    if (ice[i] > 0 && !openNeighbour) interiorIce = Math.max(interiorIce, Math.abs(f));
    if (ice[i] === 0 && icyNeighbour) edgeWater = Math.min(edgeWater, f);
  }
  assert.ok(Math.abs(total) < 1e-9 * magnitude, `net convergence ${total / magnitude} of the gross`);
  assert.ok(edgeIce > 10, `ice at the edge receives ${edgeIce} W/m²`);
  assert.ok(interiorIce < 1e-9, `ice interior receives ${interiorIce} W/m²`);
  assert.ok(edgeWater < -10, `open water at the edge loses ${edgeWater} W/m²`);
});

test('open water is dark under a high sun and bright near the horizon', () => {
  assert.ok(openWaterAlbedo(1) < 0.03);
  assert.ok(Math.abs(openWaterAlbedo(0.5) - 0.07) < 0.01);
  assert.ok(openWaterAlbedo(0.1) > 0.25 && openWaterAlbedo(0.1) < 0.4);
  for (let mu = 0.05; mu < 0.8; mu += 0.05) assert.ok(openWaterAlbedo(mu) > openWaterAlbedo(mu + 0.05));
  assert.ok(openWaterAlbedo(1) < openWaterAlbedo(0.5));
  assert.ok(seaIce.albedo(0, 0.1) > seaIce.albedo(0, 0.9));
  assert.ok(Math.abs(seaIce.albedo(0) - 0.06) < 1e-12, 'diffuse light sees the diffuse albedo');
  assert.equal(seaIce.albedo(1, 0.1), seaIce.albedo(1, 0.9));
});

test('albedo rises from open water to thick ice', () => {
  assert.ok(seaIce.albedo(0) < 0.1);
  assert.ok(seaIce.albedo(0.25) > seaIce.albedo(0) && seaIce.albedo(0.25) < seaIce.albedo(1));
  assert.equal(seaIce.albedo(1), seaIce.albedo(3));
});

test('the initial state has ice at the poles and the model steps with it', () => {
  const m = createModel(new Grid(4));
  const init = initializeState(m, {});
  for (let a = 0; a < init.length; a++) m.state[a].set(init[a]);
  const d0 = m.diagnostics();
  assert.ok(d0.iceFraction > 0.02 && d0.iceFraction < 0.3, `initial ice fraction ${d0.iceFraction}`);
  for (let n = 0; n < 24; n++) m.step(1800);
  const d = m.diagnostics();
  assert.ok(Number.isFinite(d.maxWind) && d.maxWind < 80);
  assert.ok(d.planetaryAlbedo > 0.05 && d.planetaryAlbedo < 0.6, `planetary albedo ${d.planetaryAlbedo}`);
  assert.ok(d.iceFraction > 0);
  console.log(`N=4 after 12 h: ice fraction ${d.iceFraction.toFixed(3)}, mean thickness ${d.iceThickness.toFixed(2)} m, surface albedo ${d.surfaceAlbedo.toFixed(3)}, planetary albedo ${d.planetaryAlbedo.toFixed(3)}, TCW ${(1000 * d.columnCloud).toFixed(0)} g/m²`);
});
