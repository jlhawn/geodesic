import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid, GridVertex } from '../js/grid.module.js';

const SIZES = (process.env.GRID_TEST_N ?? '2,3,5,8,16').split(',').map(Number);
const RELAX = (process.env.GRID_TEST_RELAX ?? '0,10').split(',').map(Number);
const EPS = 1e-12;

for (const relax of RELAX) for (const N of SIZES) {
  const grid = new Grid(N, { relax });
  const cells = [...grid];
  const C = 10 * N * N + 2;

  test(`N=${N} relax=${relax}: cell count and index follow iteration order`, () => {
    assert.equal(grid.size, C);
    assert.equal(cells.length, C);
    cells.forEach((cell, i) => assert.equal(cell.index, i));
    assert.equal(cells[0], grid.northPole);
    assert.equal(cells[C - 1], grid.southPole);
    for (const cell of cells) assert.ok(Math.abs(cell.centerVertex.length() - 1) < EPS);
  });

  test(`N=${N} relax=${relax}: neighbors are symmetric and distinct, 5 on the 12 pentagons and 6 elsewhere`, () => {
    let pentagons = 0;
    for (const cell of cells) {
      const nb = cell.neighbors;
      assert.equal(nb.length, cell.isPentagon ? 5 : 6);
      assert.equal(new Set(nb).size, nb.length);
      for (const other of nb) assert.ok(other.neighbors.includes(cell));
      if (cell.isPentagon) pentagons++;
    }
    assert.equal(pentagons, 12);
  });

  test(`N=${N} relax=${relax}: neighbors and vertices are counter-clockwise viewed from outside`, () => {
    for (const cell of cells) {
      const c = cell.centerVertex;
      const nb = cell.neighbors;
      const v = cell.vertices;
      const n = nb.length;
      for (let k = 0; k < n; k++) {
        const nbTurn = nb[k].centerVertex.clone().sub(c).cross(nb[(k + 1) % n].centerVertex.clone().sub(c)).dot(c);
        assert.ok(nbTurn > 0);
        const vTurn = v[k].clone().sub(c).cross(v[(k + 1) % n].clone().sub(c)).dot(c);
        assert.ok(vTurn > 0);
      }
    }
  });

  test(`N=${N} relax=${relax}: vertices are unit circumcenters, each shared by exactly three cells`, () => {
    assert.equal(grid.vertices.length, 2 * C - 4);
    grid.vertices.forEach((v, i) => {
      assert.ok(v instanceof GridVertex);
      assert.equal(v.index, i);
      assert.ok(Math.abs(v.length() - 1) < EPS);
    });
    const uses = new Array(grid.vertices.length).fill(0);
    for (const cell of cells) {
      const nb = cell.neighbors;
      const n = nb.length;
      for (let k = 0; k < n; k++) {
        const a = nb[k];
        const b = nb[(k + 1) % n];
        const v = cell.vertices[k];
        uses[v.index]++;
        const d = [cell, a, b].map((x) => v.angleTo(x.centerVertex));
        assert.ok(Math.max(...d) - Math.min(...d) < EPS);
        assert.equal(a.vertices[a.neighbors.indexOf(b)], v);
        assert.equal(b.vertices[b.neighbors.indexOf(cell)], v);
      }
    }
    assert.ok(uses.every((u) => u === 3));
  });

  test(`N=${N} relax=${relax}: the edge to neighbors[k] runs from vertices[k-1] to vertices[k], perpendicular to the dual edge`, () => {
    let edges = 0;
    for (const cell of cells) {
      const c = cell.centerVertex;
      const nb = cell.neighbors;
      const v = cell.vertices;
      const n = nb.length;
      for (let k = 0; k < n; k++) {
        const j = nb[k].centerVertex;
        const v1 = v[(k - 1 + n) % n];
        const v2 = v[k];
        for (const p of [v1, v2]) assert.ok(Math.abs(p.angleTo(c) - p.angleTo(j)) < EPS);
        const dual = j.clone().sub(c);
        const primal = v2.clone().sub(v1);
        assert.ok(Math.abs(dual.dot(primal)) / (dual.length() * primal.length()) < EPS);
        if (cell.index < nb[k].index) edges++;
      }
    }
    assert.equal(edges, 3 * C - 6);
  });

  test(`N=${N} relax=${relax}: spherical cell areas are positive and sum to 4π`, () => {
    let sum = 0;
    for (const cell of cells) {
      assert.ok(cell.area > 0);
      sum += cell.area;
    }
    assert.ok(Math.abs(sum - 4 * Math.PI) < EPS * 4 * Math.PI);
  });
}
