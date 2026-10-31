export const MAX_EDGES = 6, MAX_EDGES_ON_EDGE = 10;

/*
 * The mesh as the kernels read it: MI the connectivity and MF the geometry,
 * each a list of [name, length] sections over C cells, E edges and V
 * vertices. The unused slots of a cell's edges and of an edge's neighbors
 * hold a valid index and a zero sign or weight, so kernels loop over the
 * full MAXE and MAXEE slots. fillMeshTables writes the sections a layout
 * holds, so the atmosphere's tables and the ocean's, which cover part of the
 * mesh, come from one builder.
 */
export const meshSections = (C, E, V) => ({
  MI: [['COE', 2 * E], ['VOE', 2 * E], ['EOC', MAX_EDGES * C], ['ESC', MAX_EDGES * C], ['COC', MAX_EDGES * C], ['NEC', C], ['COV', 3 * V], ['EOV', 3 * V], ['ESV', 3 * V], ['EOE', MAX_EDGES_ON_EDGE * E], ['NEE', E]],
  MF: [['AREA', C], ['ATRI', V], ['DC', E], ['DV', E], ['FV', V], ['KAV', 3 * V], ['PVW', MAX_EDGES_ON_EDGE * E], ['NEDGE', 3 * E], ['LAT', C], ['XC', 3 * C]],
});

export function fillMeshTables(mesh, MI, MF, mi, mf) {
  const C = mesh.nCells, E = mesh.nEdges;
  const padded = (source, fill) => Int32Array.from({ length: MAX_EDGES * C }, (_, slot) => {
    const i = Math.floor(slot / MAX_EDGES), m = slot % MAX_EDGES;
    return m < mesh.nEdgesOnCell[i] ? source[mesh.maxEdges * i + m] : fill(i);
  });
  const edgesOnEdge = () => Int32Array.from({ length: MAX_EDGES_ON_EDGE * E }, (_, slot) => (slot % MAX_EDGES_ON_EDGE < mesh.nEdgesOnEdge[Math.floor(slot / MAX_EDGES_ON_EDGE)] ? mesh.edgesOnEdge[mesh.maxEdgesOnEdge * Math.floor(slot / MAX_EDGES_ON_EDGE) + slot % MAX_EDGES_ON_EDGE] : Math.floor(slot / MAX_EDGES_ON_EDGE)));
  const integers = {
    COE: () => mesh.cellsOnEdge, VOE: () => mesh.verticesOnEdge,
    EOC: () => padded(mesh.edgesOnCell, (i) => mesh.edgesOnCell[mesh.maxEdges * i]), ESC: () => padded(mesh.edgeSignOnCell, () => 0), COC: () => padded(mesh.cellsOnCell, (i) => i),
    NEC: () => mesh.nEdgesOnCell, COV: () => mesh.cellsOnVertex, EOV: () => mesh.edgesOnVertex, ESV: () => mesh.edgeSignOnVertex, EOE: edgesOnEdge, NEE: () => mesh.nEdgesOnEdge,
  };
  const reals = {
    AREA: () => mesh.areaCell, ATRI: () => mesh.areaTriangle, DC: () => mesh.dcEdge, DV: () => mesh.dvEdge, FV: () => mesh.fVertex, KAV: () => mesh.kiteAreasOnVertex,
    PVW: () => Float64Array.from({ length: MAX_EDGES_ON_EDGE * E }, (_, slot) => {
      const e = Math.floor(slot / MAX_EDGES_ON_EDGE), s = slot % MAX_EDGES_ON_EDGE;
      if (s >= mesh.nEdgesOnEdge[e]) return 0;
      const source = mesh.maxEdgesOnEdge * e + s;
      return mesh.weightsOnEdge[source] * mesh.dvEdge[mesh.edgesOnEdge[source]];
    }),
    NEDGE: () => mesh.nEdge, LAT: () => mesh.latCell, XC: () => mesh.xCell,
  };
  for (const [name, make] of Object.entries(integers)) if (MI[name] !== undefined) mi.set(make(), MI[name]);
  for (const [name, make] of Object.entries(reals)) if (MF[name] !== undefined) mf.set(make(), MF[name]);
}
