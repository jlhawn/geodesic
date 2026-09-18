import * as THREE from "./three.module.js";
import * as ISEA from "./isea.module.js";

const NORTH_POLE_QUAD_ID = 'NP';
const SOUTH_POLE_QUAD_ID = 'SP';

const _ab = new THREE.Vector3();
const _ac = new THREE.Vector3();
const _bc = new THREE.Vector3();

/**
 * Circumcenter of a triangle whose vertices lie on the unit sphere, written
 * to `target` as a unit vector on the triangle's side of the sphere: the
 * point equidistant (in arc length) from a, b and c. Degenerate (colinear or
 * coincident) inputs fall back to the normalized centroid.
 */
function circumcenter(a, b, c, target = new THREE.Vector3()) {
  target.crossVectors(_ab.subVectors(b, a), _ac.subVectors(c, a));
  if (target.lengthSq() === 0) {
    return target.copy(centroid(a, b, c)).normalize();
  }
  target.normalize();
  if (target.dot(a) < 0) target.negate();
  return target;
}

function centroid(a, b, c) {
  return a.clone().add(b).add(c).multiplyScalar(1 / 3);
}

function sphericalExcess(a, b, c) {
  return 2 * Math.atan2(_bc.crossVectors(b, c).dot(a), 1 + a.dot(b) + b.dot(c) + c.dot(a));
}

class GridVertex extends THREE.Vector3 {
  index = -1;
}

/*
 * Invariants every consumer may rely on:
 *  - `neighbors` and `vertices` are parallel arrays, counter-clockwise when
 *    viewed from outside the sphere; both have 5 entries on pentagons.
 *  - `vertices[k]` is the unit-sphere circumcenter of (centerVertex,
 *    neighbors[k], neighbors[k+1]) and is the same GridVertex object in the
 *    three cells that share it.
 *  - The edge shared with `neighbors[k]` runs from `vertices[k-1]` to
 *    `vertices[k]`.
 *  - `index` is the cell's position in Grid iteration order.
 */
class GridCell {
  constructor(N, quadId, x, y, center) {
    this.coords = [quadId, N, x, y];
    this.index = -1;
    this.isNorthPole = quadId === NORTH_POLE_QUAD_ID;
    this.isSouthPole = quadId === SOUTH_POLE_QUAD_ID;
    this.isPole = this.isNorthPole || this.isSouthPole;
    this.centerVertex = center;
    this.isPentagon = this.isPole || (x === N - 1 && y === 0);
    this.isAlongIcosahedronEdge = this.isPentagon || (x === N - 1) || (y === 0) || (x + y === N - 1);
    this.neighbors = null;
    this.vertices = null;
    this.area = 0;
  }

  get id() {
    const [quadId, N, x, y] = this.coords;
    return `${quadId}-${N}-${x}-${y}`;
  }

  calculateNeighbors(grid) {
    this.neighbors = this.#lookupNeighbors(grid);
  }

  #lookupNeighbors({ quadCells, northPole, southPole }) {
    const [quadId, N, x, y] = this.coords;

    if (quadId >= 0 && (0 < x && x < N - 1) && (0 < y && y < N - 1)) {
      // Common case: all neighbors are on the same quad.
      const quad = quadCells[quadId];
      return [
        quad[N * (x - 1) + y],
        quad[N * (x - 1) + y + 1],
        quad[N * x + y + 1],
        quad[N * (x + 1) + y],
        quad[N * (x + 1) + y - 1],
        quad[N * x + y - 1],
      ];
    }

    if (0 <= quadId && quadId < 5) {
      // In an edge along one of the upper quads.
      if (x === 0) {
        // Along top-left edge.
        if (y === 0) {
          // top corner.
          return [
            northPole,
            quadCells[(quadId + 4) % 5][0],
            quadCells[quadId][1],
            quadCells[quadId][N],
            quadCells[(quadId + 1) % 5][1],
            quadCells[(quadId + 1) % 5][0],
          ];
        }

        if (y < N - 1) {
          // Inner edge.
          return [
            quadCells[(quadId + 4) % 5][N * (y - 1)],
            quadCells[(quadId + 4) % 5][N * y],
            quadCells[quadId][y + 1],
            quadCells[quadId][N + y],
            quadCells[quadId][N + y - 1],
            quadCells[quadId][y - 1],
          ];
        }

        // Must be left corner. y === N-1
        return [
          quadCells[(quadId + 4) % 5][N * (N - 2)],
          quadCells[(quadId + 4) % 5][N * (N - 1)],
          quadCells[5 + ((quadId + 4) % 5)][0],
          quadCells[quadId][2 * N - 1],
          quadCells[quadId][2 * N - 2],
          quadCells[quadId][N - 2],
        ]
      }

      if (x < N - 1) {
        if (y === 0) {
          // Along top-right edge.
          return [
            quadCells[quadId][N * (x - 1)],
            quadCells[quadId][N * (x - 1) + 1],
            quadCells[quadId][N * x + 1],
            quadCells[quadId][N * (x + 1)],
            quadCells[(quadId + 1) % 5][x + 1],
            quadCells[(quadId + 1) % 5][x],
          ];
        }

        // Must be bottom-left edge. y === N-1
        return [
          quadCells[quadId][N * x - 1],
          quadCells[5 + ((quadId + 4) % 5)][N * (x - 1)],
          quadCells[5 + ((quadId + 4) % 5)][N * x],
          quadCells[quadId][N * (x + 2) - 1],
          quadCells[quadId][N * (x + 2) - 2],
          quadCells[quadId][N * (x + 1) - 2],
        ];
      }

      // Bottom right edge. x === N-1
      if (y === 0) {
        // right-corner
        return [
          quadCells[quadId][N * (N - 2)],
          quadCells[quadId][N * (N - 2) + 1],
          quadCells[quadId][N * (N - 1) + 1],
          quadCells[quadId + 5][0],
          quadCells[(quadId + 1) % 5][N - 1],
        ];
      }

      if (y < N - 1) {
        // inner right edge.
        return [
          quadCells[quadId][N * (N - 2) + y],
          quadCells[quadId][N * (N - 2) + y + 1],
          quadCells[quadId][N * (N - 1) + y + 1],
          quadCells[quadId + 5][y],
          quadCells[quadId + 5][y - 1],
          quadCells[quadId][N * (N - 1) + y - 1],
        ];
      }

      // bottom corner. y === N-1
      return [
        quadCells[quadId][N * (N - 1) - 1],
        quadCells[5 + ((quadId + 4) % 5)][N * (N - 2)],
        quadCells[5 + ((quadId + 4) % 5)][N * (N - 1)],
        quadCells[quadId + 5][N - 1],
        quadCells[quadId + 5][N - 2],
        quadCells[quadId][N * N - 2],
      ];
    }

    if (5 <= quadId && quadId < 10) {
      // In an edge along one of the lower quads.
      if (x === 0) {
        // Along top-left edge.
        if (y === 0) {
          // top corner.
          return [
            quadCells[quadId - 5][N * (N - 1)],
            quadCells[quadId - 5][N * (N - 1) + 1],
            quadCells[quadId][1],
            quadCells[quadId][N],
            quadCells[(quadId - 4) % 5][2 * N - 1],
            quadCells[(quadId - 4) % 5][N - 1],
          ];
        }

        if (y < N - 1) {
          // Inner edge.
          return [
            quadCells[quadId - 5][N * (N - 1) + y],
            quadCells[quadId - 5][N * (N - 1) + y + 1],
            quadCells[quadId][y + 1],
            quadCells[quadId][N + y],
            quadCells[quadId][N + y - 1],
            quadCells[quadId][y - 1],
          ];
        }

        // Must be left corner. y === N-1
        return [
          quadCells[quadId - 5][N * N - 1],
          quadCells[5 + ((quadId - 1) % 5)][N * (N - 1)],
          quadCells[5 + ((quadId - 1) % 5)][N * (N - 1) + 1],
          quadCells[quadId][2 * N - 1],
          quadCells[quadId][2 * N - 2],
          quadCells[quadId][N - 2],
        ]
      }

      if (x < N - 1) {
        if (y === 0) {
          // Along top-right edge.
          return [
            quadCells[quadId][N * (x - 1)],
            quadCells[quadId][N * (x - 1) + 1],
            quadCells[quadId][N * x + 1],
            quadCells[quadId][N * (x + 1)],
            quadCells[(quadId - 4) % 5][N * (x + 2) - 1],
            quadCells[(quadId - 4) % 5][N * (x + 1) - 1],
          ];
        }

        // Must be bottom-left edge. y === N-1
        return [
          quadCells[quadId][N * x - 1],
          quadCells[5 + ((quadId - 1) % 5)][N * (N - 1) + x],
          quadCells[5 + ((quadId - 1) % 5)][N * (N - 1) + x + 1],
          quadCells[quadId][N * (x + 2) - 1],
          quadCells[quadId][N * (x + 2) - 2],
          quadCells[quadId][N * (x + 1) - 2],
        ];
      }

      // Bottom right edge. x === N-1
      if (y === 0) {
        // right-corner
        return [
          quadCells[quadId][N * (N - 2)],
          quadCells[quadId][N * (N - 2) + 1],
          quadCells[quadId][N * (N - 1) + 1],
          quadCells[5 + ((quadId + 1) % 5)][N - 1],
          quadCells[(quadId - 4) % 5][N * N - 1],
        ];
      }

      if (y < N - 1) {
        // inner right edge.
        return [
          quadCells[quadId][N * (N - 2) + y],
          quadCells[quadId][N * (N - 2) + y + 1],
          quadCells[quadId][N * (N - 1) + y + 1],
          quadCells[5 + ((quadId + 1) % 5)][N * (y + 1) - 1],
          quadCells[5 + ((quadId + 1) % 5)][N * y - 1],
          quadCells[quadId][N * (N - 1) + y - 1],
        ];
      }

      // bottom corner. y === N-1
      return [
        quadCells[quadId][N * (N - 1) - 1],
        quadCells[5 + ((quadId - 1) % 5)][N * N - 1],
        southPole,
        quadCells[5 + ((quadId + 1) % 5)][N * N - 1],
        quadCells[5 + ((quadId + 1) % 5)][N * (N - 1) - 1],
        quadCells[quadId][N * N - 2],
      ];
    }

    if (this.isNorthPole) {
      return [
        quadCells[0][0],
        quadCells[1][0],
        quadCells[2][0],
        quadCells[3][0],
        quadCells[4][0],
      ];
    }

    // Must be south pole.
    return [
      quadCells[9][N * N - 1],
      quadCells[8][N * N - 1],
      quadCells[7][N * N - 1],
      quadCells[6][N * N - 1],
      quadCells[5][N * N - 1],
    ];
  }

  calculateVertices(grid) {
    const neighbors = this.neighbors;
    const n = neighbors.length;
    this.vertices = new Array(n);
    for (let k = 0; k < n; k++) {
      const a = neighbors[k];
      const b = neighbors[(k + 1) % n];
      let vertex;
      if (a.index < this.index) {
        vertex = a.vertexBetween(b, this);
      } else if (b.index < this.index) {
        vertex = b.vertexBetween(this, a);
      } else {
        vertex = circumcenter(this.centerVertex, a.centerVertex, b.centerVertex, new GridVertex());
        vertex.index = grid.vertices.length;
        grid.vertices.push(vertex);
      }
      this.vertices[k] = vertex;
    }
  }

  vertexBetween(first, second) {
    const neighbors = this.neighbors;
    const n = neighbors.length;
    for (let k = 0; k < n; k++) {
      if (neighbors[k] === first && neighbors[(k + 1) % n] === second) {
        return this.vertices[k];
      }
    }
    throw new Error(`cell ${this.id} has no vertex between ${first.id} and ${second.id}`);
  }

  calculateArea() {
    const c = this.centerVertex;
    const v = this.vertices;
    const n = v.length;
    let area = 0;
    for (let k = 0; k < n; k++) {
      area += sphericalExcess(c, v[k], v[(k + 1) % n]);
    }
    this.area = area;
  }

  calculateCentroid() {
    const c = this.centerVertex;
    const v = this.vertices;
    const n = v.length;
    const centroid = new THREE.Vector3();
    for (let k = 0; k < n; k++) {
      const a = v[k];
      const b = v[(k + 1) % n];
      _ab.copy(c).add(a).add(b).normalize();
      centroid.addScaledVector(_ab, sphericalExcess(c, a, b));
    }
    return centroid.normalize();
  }
}

class Grid {
  constructor(N, { relax = 8 } = {}) {
    const { quadCells, northPole, southPole } = Grid.make(N);

    this.N = N;
    this.quadCells = quadCells;
    this.northPole = northPole;
    this.southPole = southPole;
    this.size = 10 * N * N + 2;
    this.vertices = [];

    let index = 0;
    for (const cell of this) {
      cell.index = index++;
      cell.calculateNeighbors(this);
    }
    for (const cell of this) {
      cell.calculateVertices(this);
    }
    for (let i = 0; i < relax; i++) {
      this.relax();
    }
    for (const cell of this) {
      cell.calculateArea();
    }
  }

  /*
   * One Lloyd iteration toward a centroidal Voronoi tessellation: every
   * center moves to the centroid of its cell, then the shared vertices are
   * recomputed as circumcenters. Topology and ordering are unchanged; the
   * ISEA points become the seed rather than the result.
   */
  relax() {
    const centroids = new Array(this.size);
    for (const cell of this) {
      centroids[cell.index] = cell.calculateCentroid();
    }
    for (const cell of this) {
      cell.centerVertex.copy(centroids[cell.index]);
    }
    for (const cell of this) {
      const neighbors = cell.neighbors;
      const n = neighbors.length;
      for (let k = 0; k < n; k++) {
        if (neighbors[k].index > cell.index && neighbors[(k + 1) % n].index > cell.index) {
          circumcenter(cell.centerVertex, neighbors[k].centerVertex, neighbors[(k + 1) % n].centerVertex, cell.vertices[k]);
        }
      }
    }
  }

  *[Symbol.iterator]() {
    yield this.northPole;
    for (const quad of this.quadCells) {
      for (const cell of quad) {
        yield cell;
      }
    }
    yield this.southPole;
  }

  static make(N) {
    const ico = new ISEA.Icosahedron();
    const refGrid = Grid.makeRefQuad(ico.quads[0], N);

    const quadCells = new Array(10);
    // Then do the inverse projections for each cell center.
    for (let i = 0; i < 10; i++) {
      const upFrame = ico.quads[i].up.frame;
      const downFrame = ico.quads[i].down.frame;

      const cells = new Array(N * N);
      // The refGrid is (N+1)^2 in order to span all edges
      // of the triangle faces. The points we care about
      // range from x=1 to x=N and from y=0 to Y=N-1.
      for (let x = 1; x <= N; x++) {
        for (let y = 0; y < N; y++) {
          let refPoint, refFrame;
          if (x + y <= N) {
            // Point is in UP triangle.
            refPoint = refGrid[x][y].up;
            refFrame = upFrame;
          } else {
            // Point is in DOWN triangle.
            refPoint = refGrid[x][y].down;
            if (refPoint === undefined) {
              console.log(`ERROR!! refPoint (${x}, ${y}) is undefined.`);
            }
            refFrame = downFrame;
          }
          const center = ISEA.unprojectFromFace(refPoint, refFrame);
          cells[N * (x - 1) + y] = new GridCell(N, i, x - 1, y, center);
        }
      }
      quadCells[i] = cells;
    }

    const northPoleCenter = ISEA.unprojectFromFace(refGrid[0][0].up, ico.quads[0].up.frame);
    const southPoleCenter = ISEA.unprojectFromFace(refGrid[N][N].down, ico.quads[5].down.frame);

    const northPole = new GridCell(N, NORTH_POLE_QUAD_ID, 0, 0, northPoleCenter);
    const southPole = new GridCell(N, SOUTH_POLE_QUAD_ID, 0, 0, southPoleCenter);

    return { northPole, southPole, quadCells };
  }

  // makeRefQuad returns a grid of centerpoints for a hexagon grid
  // on the given triangle face with N hexagons along its edges along
  // with the coordinates of the vertices of those hexagons being
  // the centroids of the triangular grid.
  static makeRefQuad({ up: { tri, frame } }, N) {
    // A "Quad" (quadrilateral) is two icosahedral face triangles, one
    // pointing North and the other pointing South.
    // 
    //        A
    //       / \
    //      /   \               x→
    //     /     \            A - - - - - - C
    //    /       \         y |         o / |
    //   /    0    \        ↓ |         /   |
    //  /           \         |   *   /     |
    // B - - - - - - C   ->   |     /   *   |
    //  \           /         |   /         |
    //   \         /          | / o         |
    //    \   1   /           B - - - - - - D
    //     \     /
    //      \   /
    //       \ /
    //        D
    //
    // Because the icosahedron is symmetrical, we only need to subdivide a
    // quad once and that subdivision can be applied to all 10 quads on the
    // icosahedron. Here we'll be making an (N+1)x(N+1) square grid which 
    // contains the (u, v) coordinates for points on the flat triangular face
    // given the ISEA projection.

    // We'll then get the projected coordinates in the frame for UP triangle.
    const { a: A, b: B, c: C } = tri;
    const { x: Ax, y: Ay } = ISEA.projectVectorToFace(A, frame);
    const { x: Bx, y: By } = ISEA.projectVectorToFace(B, frame);
    const { x: Cx, y: Cy } = ISEA.projectVectorToFace(C, frame);

    // Compute vectors in our local X (A->C) and Y (A->B) directions.
    const A0 = new THREE.Vector2(Ax, Ay);
    const x0 = (new THREE.Vector2(Cx, Cy)).sub(A0).divideScalar(N);
    const y0 = (new THREE.Vector2(Bx, By)).sub(A0).divideScalar(N);

    const grid = new Array(N + 1);
    for (let col = 0; col <= N; col++) {
      grid[col] = new Array(N + 1);
      for (let row = 0; row <= N; row++) {
        grid[col][row] = {};
      }
    }

    for (let x = 0; x <= N; x++) {
      // Only traverse up to x+y <= N.
      // This keeps the loop within the face of triangle ABC.
      for (let y = 0; x + y <= N; y++) {
        const p = A0.clone();
        if (x > 0) p.add(x0.clone().multiplyScalar(x));
        if (y > 0) p.add(y0.clone().multiplyScalar(y));
        // Triangle ABC is referred to as the "up" triangle
        // while DCB is referred to as teh "down" triangle.
        grid[x][y].up = p;
        // Due to the symmetry of the triangle faces, the frame
        // coordinates at (x, y) are also the same at (N-x, N-y)
        // which is in triangle DCB. Along the diangonal BC, this
        // is equivalent to swapping the row and column.
        grid[N - x][N - y].down = p;
      }
    }

    return grid;
  }
}

export {
  Grid, GridCell, GridVertex, circumcenter, centroid
};
