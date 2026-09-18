import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { divergence, gradient, curl, tangential, kineticEnergy, cellVector } from '../js/dynamics/operators.module.js';

const SIZES = (process.env.MESH_TEST_N ?? '4,8,16').split(',').map(Number);
const RELAX = (process.env.MESH_TEST_RELAX ?? '0,10').split(',').map(Number);
const EPS = 1e-12;
const OMEGA = 1e-5;

function random(seed) {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296 - 0.5;
  };
}

function randomArray(n, seed) {
  const next = random(seed);
  return Float64Array.from({ length: n }, next);
}

function relative(residual, scale) {
  return scale === 0 ? Math.abs(residual) : Math.abs(residual) / scale;
}

function rotationComponent(mesh, direction) {
  const { nEdges, xEdge, radius } = mesh;
  const out = new Float64Array(nEdges);
  for (let e = 0; e < nEdges; e++) {
    const mx = xEdge[3 * e], my = xEdge[3 * e + 1];
    const dx = direction[3 * e], dy = direction[3 * e + 1];
    out[e] = radius * OMEGA * (-my * dx + mx * dy);
  }
  return out;
}

const convergence = {};
const reported = {};
function record(key, N, value) {
  (convergence[key] ??= []).push({ N, value });
}
function report(key, N, value) {
  (reported[key] ??= []).push({ N, value });
}
function tag(key, relax) {
  return `${key} (relax=${relax})`;
}

for (const relax of RELAX) for (const N of SIZES) {
  const grid = new Grid(N, { relax });
  const mesh = buildMesh(grid);
  const { nCells, nEdges, nVertices, maxEdges, maxEdgesOnEdge, radius } = mesh;
  const sphere = 4 * Math.PI * radius * radius;

  test(`N=${N} relax=${relax}: counts and connectivity`, () => {
    assert.equal(nCells, 10 * N * N + 2);
    assert.equal(nEdges, 3 * nCells - 6);
    assert.equal(nVertices, 2 * nCells - 4);
    const edgeUses = new Int32Array(nEdges);
    for (let i = 0; i < nCells; i++) {
      const n = mesh.nEdgesOnCell[i];
      assert.ok(n === 5 || n === 6);
      for (let k = 0; k < n; k++) {
        const e = mesh.edgesOnCell[maxEdges * i + k];
        edgeUses[e]++;
        const sign = mesh.edgeSignOnCell[maxEdges * i + k];
        assert.equal(mesh.cellsOnEdge[2 * e + (sign === 1 ? 0 : 1)], i);
        assert.equal(mesh.cellsOnEdge[2 * e + (sign === 1 ? 1 : 0)], mesh.cellsOnCell[maxEdges * i + k]);
      }
    }
    assert.ok(edgeUses.every((u) => u === 2));
    for (let e = 0; e < nEdges; e++) {
      const ni = mesh.nEdgesOnCell[mesh.cellsOnEdge[2 * e]];
      const nj = mesh.nEdgesOnCell[mesh.cellsOnEdge[2 * e + 1]];
      assert.equal(mesh.nEdgesOnEdge[e], ni + nj - 2);
      assert.notEqual(mesh.verticesOnEdge[2 * e], mesh.verticesOnEdge[2 * e + 1]);
    }
    for (let v = 0; v < nVertices; v++) {
      const cells = [0, 1, 2].map((m) => mesh.cellsOnVertex[3 * v + m]);
      assert.equal(new Set(cells).size, 3);
      for (let m = 0; m < 3; m++) {
        const e = mesh.edgesOnVertex[3 * v + m];
        assert.ok(e >= 0);
        const pair = new Set([mesh.cellsOnEdge[2 * e], mesh.cellsOnEdge[2 * e + 1]]);
        assert.ok(pair.has(cells[m]) && pair.has(cells[(m + 1) % 3]));
        assert.ok(mesh.verticesOnEdge[2 * e] === v || mesh.verticesOnEdge[2 * e + 1] === v);
      }
    }
  });

  test(`N=${N} relax=${relax}: spherical areas and kites partition exactly`, () => {
    let cellSum = 0, triangleSum = 0;
    const cells = [...grid];
    for (let i = 0; i < nCells; i++) {
      cellSum += mesh.areaCell[i];
      assert.ok(relative(mesh.areaCell[i] - cells[i].area * radius * radius, mesh.areaCell[i]) < EPS);
      let kites = 0;
      for (let k = 0; k < mesh.nEdgesOnCell[i]; k++) {
        const kite = mesh.kiteAreasOnCell[maxEdges * i + k];
        assert.ok(kite > 0);
        kites += kite;
      }
      assert.ok(relative(kites - mesh.areaCell[i], mesh.areaCell[i]) < EPS);
    }
    for (let v = 0; v < nVertices; v++) {
      triangleSum += mesh.areaTriangle[v];
      let kites = 0;
      for (let m = 0; m < 3; m++) kites += mesh.kiteAreasOnVertex[3 * v + m];
      assert.ok(relative(kites - mesh.areaTriangle[v], mesh.areaTriangle[v]) < EPS);
    }
    assert.ok(relative(cellSum - sphere, sphere) < EPS);
    assert.ok(relative(triangleSum - sphere, sphere) < EPS);
  });

  test(`N=${N} relax=${relax}: edge frames are orthogonal and oriented as documented`, () => {
    for (let e = 0; e < nEdges; e++) {
      const i = mesh.cellsOnEdge[2 * e], j = mesh.cellsOnEdge[2 * e + 1];
      const a = mesh.verticesOnEdge[2 * e], b = mesh.verticesOnEdge[2 * e + 1];
      const m = [0, 1, 2].map((c) => mesh.xEdge[3 * e + c]);
      const n = [0, 1, 2].map((c) => mesh.nEdge[3 * e + c]);
      const t = [0, 1, 2].map((c) => mesh.tEdge[3 * e + c]);
      const chord = [0, 1, 2].map((c) => mesh.xVertex[3 * b + c] - mesh.xVertex[3 * a + c]);
      const chordLength = Math.hypot(...chord);
      const dot = (p, q) => p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
      assert.ok(Math.abs(dot(n, chord)) / chordLength < 1e-10);
      assert.ok(dot(t, chord) > 0);
      assert.ok(Math.abs(dot(m, n)) < EPS && Math.abs(dot(m, t)) < EPS && Math.abs(dot(n, t)) < EPS);
      const dij = [0, 1, 2].map((c) => mesh.xCell[3 * i + c] - mesh.xCell[3 * j + c]);
      assert.ok(Math.abs(dot(m, dij)) < EPS);
      assert.ok(Math.abs(mesh.dvEdge[e] - radius * 2 * Math.asin(chordLength / 2)) < 1e-6);
    }
    for (let v = 0; v < nVertices; v++) {
      const xv = [0, 1, 2].map((c) => mesh.xVertex[3 * v + c]);
      for (let m = 0; m < 3; m++) {
        const e = mesh.edgesOnVertex[3 * v + m];
        const i = mesh.cellsOnEdge[2 * e], j = mesh.cellsOnEdge[2 * e + 1];
        const pi = [0, 1, 2].map((c) => mesh.xCell[3 * i + c] - xv[c]);
        const pj = [0, 1, 2].map((c) => mesh.xCell[3 * j + c] - xv[c]);
        const turn = xv[0] * (pi[1] * pj[2] - pi[2] * pj[1]) + xv[1] * (pi[2] * pj[0] - pi[0] * pj[2]) + xv[2] * (pi[0] * pj[1] - pi[1] * pj[0]);
        assert.equal(Math.sign(turn), mesh.edgeSignOnVertex[3 * v + m]);
      }
    }
  });

  test(`N=${N} relax=${relax}: divergence conserves mass and gradient is its negative adjoint`, () => {
    const F = randomArray(nEdges, 11);
    const phi = randomArray(nCells, 23);
    const D = divergence(mesh, F);
    let total = 0, scale = 0;
    for (let i = 0; i < nCells; i++) {
      total += mesh.areaCell[i] * D[i];
      scale += mesh.areaCell[i] * Math.abs(D[i]);
    }
    assert.ok(relative(total, scale) < EPS);
    const G = gradient(mesh, phi);
    let edges = 0, cells = 0, edgeScale = 0;
    for (let e = 0; e < nEdges; e++) {
      edges += mesh.dcEdge[e] * mesh.dvEdge[e] * G[e] * F[e];
      edgeScale += Math.abs(mesh.dcEdge[e] * mesh.dvEdge[e] * G[e] * F[e]);
    }
    for (let i = 0; i < nCells; i++) cells += mesh.areaCell[i] * phi[i] * D[i];
    assert.ok(relative(edges + cells, edgeScale) < EPS);
    const rotation = rotationComponent(mesh, mesh.nEdge);
    const Drot = divergence(mesh, rotation);
    report(tag('divergence of solid rotation, max / omega', relax), N, Math.max(...Drot.map(Math.abs)) / OMEGA);
  });

  test(`N=${N} relax=${relax}: curl of solid-body rotation is 2ω sin(lat), and integrates to zero`, () => {
    const rotation = rotationComponent(mesh, mesh.nEdge);
    const zeta = curl(mesh, rotation);
    let worst = 0;
    for (let v = 0; v < nVertices; v++) {
      worst = Math.max(worst, Math.abs(zeta[v] - 2 * OMEGA * Math.sin(mesh.latVertex[v])) / (2 * OMEGA));
    }
    report(tag('curl of solid rotation, max relative error', relax), N, worst);
    assert.ok(worst < 0.05);
    const zetaRandom = curl(mesh, randomArray(nEdges, 37));
    let total = 0, scale = 0;
    for (let v = 0; v < nVertices; v++) {
      total += mesh.areaTriangle[v] * zetaRandom[v];
      scale += mesh.areaTriangle[v] * Math.abs(zetaRandom[v]);
    }
    assert.ok(relative(total, scale) < EPS);
  });

  test(`N=${N} relax=${relax}: T1 — TRiSK weights are antisymmetric and the Coriolis term does no work`, () => {
    let worst = 0, scale = 0;
    for (let e = 0; e < nEdges; e++) {
      for (let s = 0; s < mesh.nEdgesOnEdge[e]; s++) {
        const other = mesh.edgesOnEdge[maxEdgesOnEdge * e + s];
        const w = mesh.weightsOnEdge[maxEdgesOnEdge * e + s];
        scale = Math.max(scale, Math.abs(w));
        let back = null;
        for (let r = 0; r < mesh.nEdgesOnEdge[other]; r++) {
          if (mesh.edgesOnEdge[maxEdgesOnEdge * other + r] === e) back = mesh.weightsOnEdge[maxEdgesOnEdge * other + r];
        }
        assert.notEqual(back, null);
        worst = Math.max(worst, Math.abs(w + back));
      }
    }
    assert.ok(worst / scale < 1e-13);
    const u = randomArray(nEdges, 41);
    const uPerp = tangential(mesh, u);
    let work = 0, energy = 0;
    for (let e = 0; e < nEdges; e++) {
      work += mesh.dcEdge[e] * mesh.dvEdge[e] * u[e] * uPerp[e];
      energy += mesh.dcEdge[e] * mesh.dvEdge[e] * u[e] * u[e];
    }
    assert.ok(relative(work, energy) < EPS);
  });

  test(`N=${N} relax=${relax}: T2 — tangential reconstruction of solid-body rotation converges`, () => {
    const u = rotationComponent(mesh, mesh.nEdge);
    const exact = rotationComponent(mesh, mesh.tEdge);
    const uPerp = tangential(mesh, u);
    let worst = 0, sumError = 0, sumExact = 0;
    for (let e = 0; e < nEdges; e++) {
      worst = Math.max(worst, Math.abs(uPerp[e] - exact[e]));
      sumError += Math.abs(uPerp[e] - exact[e]);
      sumExact += Math.abs(exact[e]);
    }
    const speed = radius * OMEGA;
    report(tag('T2 max tangential error / speed', relax), N, worst / speed);
    record(tag('T2 mean tangential error / mean speed', relax), N, sumError / sumExact);
    assert.ok(worst / speed < 0.15);
  });

  test(`N=${N} relax=${relax}: T3 — flux out of every dual triangle equals its kite-weighted divergence`, () => {
    const u = randomArray(nEdges, 53);
    const uPerp = tangential(mesh, u);
    const D = divergence(mesh, u);
    let worst = 0;
    for (let v = 0; v < nVertices; v++) {
      let flux = 0, scale = 0;
      for (let m = 0; m < 3; m++) {
        const e = mesh.edgesOnVertex[3 * v + m];
        const term = mesh.edgeSignOnVertex[3 * v + m] * uPerp[e] * mesh.dcEdge[e];
        const kite = mesh.kiteAreasOnVertex[3 * v + m] * D[mesh.cellsOnVertex[3 * v + m]];
        flux += term + kite;
        scale += Math.abs(term) + Math.abs(kite);
      }
      worst = Math.max(worst, relative(flux, scale));
    }
    assert.ok(worst < EPS);
  });

  test(`N=${N} relax=${relax}: kinetic energy and cell-vector reconstruction of solid-body rotation`, () => {
    const u = rotationComponent(mesh, mesh.nEdge);
    const K = kineticEnergy(mesh, u);
    const V = cellVector(mesh, u);
    let worstK = 0, worstV = 0;
    for (let i = 0; i < nCells; i++) {
      const x = mesh.xCell[3 * i], y = mesh.xCell[3 * i + 1];
      const vx = -radius * OMEGA * y, vy = radius * OMEGA * x;
      const speed2 = vx * vx + vy * vy;
      if (speed2 < 0.01 * (radius * OMEGA) ** 2) continue;
      worstK = Math.max(worstK, Math.abs(K[i] - 0.5 * speed2) / (0.5 * speed2));
      worstV = Math.max(worstV, Math.hypot(V[3 * i] - vx, V[3 * i + 1] - vy, V[3 * i + 2]) / Math.sqrt(speed2));
    }
    report(tag('KE max relative error', relax), N, worstK);
    report(tag('cell vector max relative error', relax), N, worstV);
    assert.ok(worstK < 0.2 && worstV < 0.2);
  });
}

test('mean truncation errors decrease with resolution', () => {
  for (const [key, series] of Object.entries(convergence)) {
    for (let k = 1; k < series.length; k++) {
      assert.ok(series[k].value < series[k - 1].value, `${key}: N=${series[k - 1].N} ${series[k - 1].value} → N=${series[k].N} ${series[k].value}`);
    }
    console.log(key, series.map((s) => `N=${s.N}: ${s.value.toExponential(2)}`).join(', '));
  }
  for (const [key, series] of Object.entries(reported)) {
    console.log(key, series.map((s) => `N=${s.N}: ${s.value.toExponential(2)}`).join(', '));
  }
});
