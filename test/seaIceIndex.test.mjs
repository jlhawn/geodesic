import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { inIndexGrid, indexCells, seaIceIndexCut } from '../js/seaIceIndex.module.js';

const at = (lat, lon) => inIndexGrid(lat * Math.PI / 180, lon * Math.PI / 180);

test('the index\x27s grids hold the frozen seas and leave the closed basins and the far south out', () => {
  for (const [name, lat, lon] of [['pole', 90, 0], ['Bohai', 39, 120], ['Baltic', 55, 20], ['Okhotsk', 55, 150], ['Hudson Bay', 60, -85], ['St Lawrence', 48, -62], ['Japan Sea', 42, 135], ['Weddell', -70, -40], ['Ross', -75, 180]]) assert.ok(at(lat, lon), name);
  for (const [name, lat, lon] of [['Caspian', 42, 50.5], ['north Caspian', 46.5, 51], ['Black Sea', 44, 35], ['Azov', 46, 36.5], ['50S 0E', -50, 0], ['30S', -30, 20], ['equator', 0, -140]]) assert.ok(!at(lat, lon), name);
});

test('the cut counts cells of the index\x27s grids at 15 % or more, extent by area and area by cover', () => {
  const mesh = gridMesh(new Grid(16)), C = mesh.nCells, deg = 180 / Math.PI;
  const land = Uint8Array.from({ length: C }, (_, i) => (i % 5 === 0 ? 1 : 0)), lake = Uint8Array.from({ length: C }, (_, i) => (i % 7 === 0 ? 1 : 0));
  const cells = indexCells(mesh, { land, lake }), cover = Float64Array.from({ length: C }, (_, i) => (i % 3) * 0.1);
  let n = 0, a = 0;
  for (let i = 0; i < C; i++) {
    assert.equal(cells[i], !land[i] && !lake[i] && at(mesh.latCell[i] * deg, mesh.lonCell[i] * deg) ? 1 : 0);
    if (cells[i] && cover[i] >= 0.15 && mesh.latCell[i] > 0) { n += mesh.areaCell[i] / 1e12; a += cover[i] * mesh.areaCell[i] / 1e12; }
  }
  const { N } = seaIceIndexCut(mesh, cells, cover);
  assert.ok(n > 0 && Math.abs(N.extent - n) < 1e-9 && Math.abs(N.area - a) < 1e-9);
  assert.ok(cells.some((x, i) => !x && !land[i] && !lake[i] && mesh.latCell[i] > 0), 'some northern sea cells lie off the grid');
});
