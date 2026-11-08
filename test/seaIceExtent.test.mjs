import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { readTopography, geographyOptions, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { cpuFrameFields } from '../js/cpuFrame.module.js';
import { decodeSeaIceStart, SEA_ICE_START_FILE } from '../js/physics/seaIceStart.module.js';
import { seaCells, seaIceExtent, polarIcedCells } from '../js/seaIceExtent.module.js';

const LAKES = [731, 3622], SUPERIOR = 3622, deg = 180 / Math.PI;

const observedStart = (lake) => {
  const start = decodeSeaIceStart(readFileSync(new URL(`../${SEA_ICE_START_FILE}`, import.meta.url)));
  const topography = readTopography(FRESH_GEOGRAPHY.topography), sheets = geographyOptions(FRESH_GEOGRAPHY);
  const model = createModel(new Grid(32), { topography, geography: sheets.geography, ocean: false, land: { ...sheets.land }, levels: sigmaInterfaces('bl48'), ...(lake ? { lake } : {}) });
  initializeState(model, { geostrophic: !model.surfaceGeopotential, seaIce: start }).forEach((values, a) => model.state[a].set(values));
  model.land.initialize();
  return model;
};

const extentOf = (mesh, land, concentration, skip) => {
  const out = [0, 0];
  for (let i = 0; i < mesh.nCells; i++) if (!land[i] && !skip.includes(i) && concentration[i] >= 0.15) out[mesh.latCell[i] > 0 ? 0 : 1] += mesh.areaCell[i];
  return out.map((a) => a / 1e12);
};
const icedOf = (mesh, land, ice, skip) => [...ice.keys()].filter((i) => !land[i] && !skip.includes(i) && ice[i] > 0 && Math.abs(mesh.latCell[i] * deg) >= 60);

test("at N=32 on the n1024 mask at the observed fresh start, spinup.mjs's sea-ice extents and iced polar cells leave out the frozen lake the lake model holds, and a lake poleward of 60°, and with the lake model off count both as sea", () => {
  for (const lake of [null, { model: 'off' }]) {
    const model = observedStart(lake), { mesh } = model, land = model.geography.land, name = model.lake ? 'FLake' : 'the lake model off';
    assert.equal(model.lake === null, lake !== null, 'the lake model runs unless it is off');
    assert.deepEqual([...model.geography.lake.keys()].filter((i) => model.geography.lake[i]), LAKES, 'the mask holds Lake Victoria and Lake Superior');
    const held = model.lake ? LAKES : [];
    const { concentration, ice } = cpuFrameFields(model, { fields: ['concentration', 'ice'] });
    assert.ok(concentration[SUPERIOR] >= 0.15 && ice[SUPERIOR] > 0 && mesh.latCell[SUPERIOR] > 0, `${name}: Lake Superior starts frozen, cover ${concentration[SUPERIOR]}`);

    const sea = seaCells(model);
    assert.deepEqual(Array.from(sea), Array.from(land, (l, i) => (l || held.includes(i) ? 0 : 1)), `${name}: the sea cells are the cells off land${held.length ? ' but the two lakes' : ''}`);
    const extent = seaIceExtent(mesh, sea, concentration), all = extentOf(mesh, land, concentration, []);
    assert.deepEqual(extent, extentOf(mesh, land, concentration, held), `${name}: the extents count the sea cells at least 15 % covered`);
    console.log(`${name}: extent N ${extent[0].toFixed(3)} S ${extent[1].toFixed(3)} Mkm², every cell off land N ${all[0].toFixed(3)} S ${all[1].toFixed(3)}`);
    if (held.length) assert.ok(all[0] - extent[0] > 0.04, `FLake: Lake Superior's cell leaves the northern extent, ${all[0] - extent[0]} Mkm²`);
    else assert.deepEqual(extent, all, 'with the lake model off the frozen lake counts');

    const iced = polarIcedCells(mesh, sea, ice);
    assert.deepEqual(iced, icedOf(mesh, land, ice, held), `${name}: the iced cells poleward of 60°`);
    assert.ok(iced.some((i) => mesh.latCell[i] > 0) && iced.some((i) => mesh.latCell[i] < 0), `${name}: ice lies poleward of 60° in each hemisphere`);

    const polar = iced.find((i) => mesh.latCell[i] > 0 && concentration[i] >= 0.15), lakes = Uint8Array.from(model.geography.lake);
    lakes[polar] = 1;
    const withPolarLake = seaCells({ geography: { land, lake: lakes }, lake: model.lake ? { isLake: lakes } : null }), dropped = model.lake ? [...LAKES, polar] : [];
    const where = `a lake at cell ${polar}, ${(mesh.latCell[polar] * deg).toFixed(1)}°N, under ice ${ice[polar].toFixed(2)} m thick`;
    assert.deepEqual(polarIcedCells(mesh, withPolarLake, ice), icedOf(mesh, land, ice, dropped), `${name}: ${where}, and the iced cells`);
    assert.deepEqual(seaIceExtent(mesh, withPolarLake, concentration), extentOf(mesh, land, concentration, dropped), `${name}: ${where}, and the extents`);
    assert.equal(polarIcedCells(mesh, withPolarLake, ice).includes(polar), !model.lake, `${name}: ${where} is ${model.lake ? 'left out' : 'counted'}`);
  }
});
