export function divergence(mesh, flux, out = new Float64Array(mesh.nCells)) {
  const { nCells, maxEdges, nEdgesOnCell, edgesOnCell, edgeSignOnCell, dvEdge, areaCell } = mesh;
  for (let i = 0; i < nCells; i++) {
    let sum = 0;
    for (let k = 0; k < nEdgesOnCell[i]; k++) {
      const e = edgesOnCell[maxEdges * i + k];
      sum += edgeSignOnCell[maxEdges * i + k] * flux[e] * dvEdge[e];
    }
    out[i] = sum / areaCell[i];
  }
  return out;
}

export function gradient(mesh, phi, out = new Float64Array(mesh.nEdges)) {
  const { nEdges, cellsOnEdge, dcEdge } = mesh;
  for (let e = 0; e < nEdges; e++) {
    out[e] = (phi[cellsOnEdge[2 * e + 1]] - phi[cellsOnEdge[2 * e]]) / dcEdge[e];
  }
  return out;
}

export function curl(mesh, u, out = new Float64Array(mesh.nVertices)) {
  const { nVertices, edgesOnVertex, edgeSignOnVertex, dcEdge, areaTriangle } = mesh;
  for (let v = 0; v < nVertices; v++) {
    let sum = 0;
    for (let m = 0; m < 3; m++) {
      const e = edgesOnVertex[3 * v + m];
      sum += edgeSignOnVertex[3 * v + m] * u[e] * dcEdge[e];
    }
    out[v] = sum / areaTriangle[v];
  }
  return out;
}

export function tangential(mesh, u, out = new Float64Array(mesh.nEdges)) {
  const { nEdges, maxEdgesOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge, dvEdge, dcEdge } = mesh;
  for (let e = 0; e < nEdges; e++) {
    let sum = 0;
    for (let s = 0; s < nEdgesOnEdge[e]; s++) {
      const other = edgesOnEdge[maxEdgesOnEdge * e + s];
      sum += weightsOnEdge[maxEdgesOnEdge * e + s] * dvEdge[other] * u[other];
    }
    out[e] = sum / dcEdge[e];
  }
  return out;
}

export function kineticEnergy(mesh, u, out = new Float64Array(mesh.nCells)) {
  const { nCells, maxEdges, nEdgesOnCell, edgesOnCell, dcEdge, dvEdge, areaCell } = mesh;
  for (let i = 0; i < nCells; i++) {
    let sum = 0;
    for (let k = 0; k < nEdgesOnCell[i]; k++) {
      const e = edgesOnCell[maxEdges * i + k];
      sum += 0.25 * dcEdge[e] * dvEdge[e] * u[e] * u[e];
    }
    out[i] = sum / areaCell[i];
  }
  return out;
}

export function cellVector(mesh, u, out = new Float64Array(3 * mesh.nCells)) {
  const { nCells, maxEdges, nEdgesOnCell, edgesOnCell, dcEdge, dvEdge, nEdge, areaCell } = mesh;
  for (let i = 0; i < nCells; i++) {
    let x = 0, y = 0, z = 0;
    for (let k = 0; k < nEdgesOnCell[i]; k++) {
      const e = edgesOnCell[maxEdges * i + k];
      const w = 0.5 * dcEdge[e] * dvEdge[e] * u[e];
      x += w * nEdge[3 * e];
      y += w * nEdge[3 * e + 1];
      z += w * nEdge[3 * e + 2];
    }
    out[3 * i] = x / areaCell[i];
    out[3 * i + 1] = y / areaCell[i];
    out[3 * i + 2] = z / areaCell[i];
  }
  return out;
}

export function laplacianVelocity(mesh, u, out = new Float64Array(mesh.nEdges), div = new Float64Array(mesh.nCells), vort = new Float64Array(mesh.nVertices)) {
  const { nEdges, cellsOnEdge, verticesOnEdge, dcEdge, dvEdge } = mesh;
  divergence(mesh, u, div);
  curl(mesh, u, vort);
  for (let e = 0; e < nEdges; e++) {
    out[e] = (div[cellsOnEdge[2 * e + 1]] - div[cellsOnEdge[2 * e]]) / dcEdge[e]
      - (vort[verticesOnEdge[2 * e + 1]] - vort[verticesOnEdge[2 * e]]) / dvEdge[e];
  }
  return out;
}
