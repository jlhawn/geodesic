import * as THREE from "./three.module.js";
import * as ISEA from "./isea.module.js";

/**
 * Computes the circumcenter of a triangle on (or near) the unit sphere.
 *
 * The circumcenter is the unique point equidistant from the triangle’s
 * three vertices a, b, and c. In Euclidean space, it is the intersection
 * of the triangle’s perpendicular bisectors; for small spherical
 * triangles this Cartesian version is an excellent approximation
 * before normalizing back to the sphere.
 *
 * Implementation notes:
 *  - Uses the vector cross-product formulation found in
 *    *Shewchuk, J. R., “Lecture Notes on Geometric Robustness,”*
 *    and in *Eberly, D. H., “Intersection of Three Planes,” Geometric Tools, 2008*.
 *  - The formula avoids explicit matrix inversion and is robust for
 *    acute and obtuse triangles alike.
 *  - Falls back to the centroid if the points are degenerate (colinear
 *    or nearly coincident).
 *
 * Returns: THREE.Vector3 — the unnormalized circumcenter (normalize it
 * externally if you need it on the unit sphere).
 */
function circumcenter(a, b, c) {
  const tmp1 = new THREE.Vector3();
  const tmp2 = new THREE.Vector3();
  const tmp3 = new THREE.Vector3();

  // ba = b - a, ca = c - a
  const ba = tmp1.copy(b).sub(a);
  const ca = tmp2.copy(c).sub(a);

  // n = ba x ca
  const n = tmp3.copy(ba).cross(ca);
  const denom = 2 * n.lengthSq();
  if (denom === 0) {
    // Degenerate: fall back to centroid (then normalize by caller)
    return centroid(a, b, c);
  }

  // weights from |v|^2 and cross products; see Eberly/Shewchuk-style formula
  const ba2 = ba.lengthSq();
  const ca2 = ca.lengthSq();

  // u = (ba2 * (ca x n) + ca2 * (n x ba)) / denom
  const u = new THREE.Vector3()
    .copy(ca).cross(n).multiplyScalar(ba2)
    .add( new THREE.Vector3().copy(n).cross(ba).multiplyScalar(ca2) )
    .multiplyScalar(1/denom);

  return new THREE.Vector3().copy(a).add(u);
}

function centroid(a, b, c) {
  return a.clone().add(b).add(c).multiplyScalar(1/3);
}

class GridCell {
  constructor(N, quadId, x, y, center) {
    this.coords = [quadId, N, x, y];
    this.isNorthPole = quadId === -1;
    this.isSouthPole = quadId === -2;
    this.id = `${quadId}-${N}-${x}-${y}`;
    this.centerVertex = center;
    this.isPentagon = quadId < 0 || (x === N-1 && y === 0);
    this.isAlongIcosahedronEdge = this.isPentagon || (x === N-1) || (y === 0) || (x+y === N-1);
    this.vertices = null;
    this.faceTriangles = null;
    this.area = 0;
  }

  neighbors({ quadCells, northPole, southPole }) {
    const [quadId, N, x, y] = this.coords;

    if (quadId >= 0 && (0 < x && x < N-1) && (0 < y && y < N-1)) {
      // Common case: all neighbors are on the same quad.
      const quad = quadCells[quadId];
      return [
        quad[N*(x-1) + y],
        quad[N*(x-1) + y+1],
        quad[N*x + y+1],
        quad[N*(x+1) + y],
        quad[N*(x+1) + y-1],
        quad[N*x + y-1],
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
            quadCells[(quadId+4)%5][0],
            quadCells[quadId][1],
            quadCells[quadId][N],
            quadCells[(quadId+1)%5][1],
            quadCells[(quadId+1)%5][0],
          ];
        }

        if (y < N-1) {
          // Inner edge.
          return [
            quadCells[(quadId+4)%5][N*(y-1)],
            quadCells[(quadId+4)%5][N*y],
            quadCells[quadId][y+1],
            quadCells[quadId][N + y],
            quadCells[quadId][N + y-1],
            quadCells[quadId][y-1],
          ];
        }

        // Must be left corner. y === N-1
        return [
          quadCells[(quadId+4)%5][N*(N-2)],
          quadCells[(quadId+4)%5][N*(N-1)],
          quadCells[5+((quadId+4)%5)][0],
          quadCells[quadId][2*N - 1],
          quadCells[quadId][2*N - 2],
          quadCells[quadId][N-2],
        ]
      }

      if (x < N-1) {
        if (y === 0) {
          // Along top-right edge.
          return [
            quadCells[quadId][N*(x-1)],
            quadCells[quadId][N*(x-1) + 1],
            quadCells[quadId][N*x + 1],
            quadCells[quadId][N*(x+1)],
            quadCells[(quadId+1)%5][x+1],
            quadCells[(quadId+1)%5][x],
          ];
        }

        // Must be bottom-left edge. y === N-1
        return [
            quadCells[quadId][N*x - 1],
            quadCells[5+((quadId+4)%5)][N*(x-1)],
            quadCells[5+((quadId+4)%5)][N*x],
            quadCells[quadId][N*(x+2) - 1],
            quadCells[quadId][N*(x+2) - 2],
            quadCells[quadId][N*(x+1) - 2],
          ];
      }

      // Bottom right edge. x === N-1
      if (y === 0) {
        // right-corner
        return [
          quadCells[quadId][N*(N-2)],
          quadCells[quadId][N*(N-2)+1],
          quadCells[quadId][N*(N-1)+1],
          quadCells[quadId+5][0],
          quadCells[(quadId+1)%5][N-1],
        ];
      }

      if (y < N-1) {
        // inner right edge.
        return [
          quadCells[quadId][N*(N-2) + y],
          quadCells[quadId][N*(N-2) + y+1],
          quadCells[quadId][N*(N-1) + y+1],
          quadCells[quadId+5][y],
          quadCells[quadId+5][y-1],
          quadCells[quadId][N*(N-1) + y-1],
        ];
      }

      // bottom corner. y === N-1
      return [
        quadCells[quadId][N*(N-1) - 1],
        quadCells[5+((quadId+4)%5)][N*(N-2)],
        quadCells[5+((quadId+4)%5)][N*(N-1)],
        quadCells[quadId+5][N-1],
        quadCells[quadId+5][N-2],
        quadCells[quadId][N*N - 2],
      ];
    }

    if (5 <= quadId && quadId < 10) {
      // In an edge along one of the lower quads.
      if (x === 0) {
        // Along top-left edge.
        if (y === 0) {
          // top corner.
          return [
            quadCells[quadId-5][N*(N-1)],
            quadCells[quadId-5][N*(N-1) + 1],
            quadCells[quadId][1],
            quadCells[quadId][N],
            quadCells[(quadId-4)%5][2*N - 1],
            quadCells[(quadId-4)%5][N - 1],
          ];
        }

        if (y < N-1) {
          // Inner edge.
          return [
            quadCells[quadId-5][N*(N-1) + y],
            quadCells[quadId-5][N*(N-1) + y+1],
            quadCells[quadId][y+1],
            quadCells[quadId][N + y],
            quadCells[quadId][N + y-1],
            quadCells[quadId][y-1],
          ];
        }

        // Must be left corner. y === N-1
        return [
          quadCells[quadId-5][N*N - 1],
          quadCells[5+((quadId-1)%5)][N*(N-1)],
          quadCells[5+((quadId-1)%5)][N*(N-1) + 1],
          quadCells[quadId][2*N - 1],
          quadCells[quadId][2*N - 2],
          quadCells[quadId][N - 2],
        ]
      }

      if (x < N-1) {
        if (y === 0) {
          // Along top-right edge.
          return [
            quadCells[quadId][N*(x-1)],
            quadCells[quadId][N*(x-1) + 1],
            quadCells[quadId][N*x + 1],
            quadCells[quadId][N*(x+1)],
            quadCells[(quadId-4)%5][N*(x+2) - 1],
            quadCells[(quadId-4)%5][N*(x+1) - 1],
          ];
        }

        // Must be bottom-left edge. y === N-1
        return [
            quadCells[quadId][N*x - 1],
            quadCells[5+((quadId-1)%5)][N*(N-1) + x],
            quadCells[5+((quadId-1)%5)][N*(N-1) + x+1],
            quadCells[quadId][N*(x+2) - 1],
            quadCells[quadId][N*(x+2) - 2],
            quadCells[quadId][N*(x+1) - 2],
          ];
      }

      // Bottom right edge. x === N-1
      if (y === 0) {
        // right-corner
        return [
          quadCells[quadId][N*(N-2)],
          quadCells[quadId][N*(N-2) + 1],
          quadCells[quadId][N*(N-1) + 1],
          quadCells[5+((quadId+1)%5)][N-1],
          quadCells[(quadId-4)%5][N*N - 1],
        ];
      }

      if (y < N-1) {
        // inner right edge.
        return [
          quadCells[quadId][N*(N-2) + y],
          quadCells[quadId][N*(N-2) + y+1],
          quadCells[quadId][N*(N-1) + y+1],
          quadCells[5+((quadId+1)%5)][N*(y+1) - 1],
          quadCells[5+((quadId+1)%5)][N*y - 1],
          quadCells[quadId][N*(N-1) + y-1],
        ];
      }

      // bottom corner. y === N-1
      return [
        quadCells[quadId][N*(N-1) - 1],
        quadCells[5+((quadId-1)%5)][N*N - 1],
        southPole,
        quadCells[5+((quadId+1)%5)][N*N - 1],
        quadCells[5+((quadId+1)%5)][N*(N-1) - 1],
        quadCells[quadId][N*N - 2],
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
      quadCells[9][N*N - 1],
      quadCells[8][N*N - 1],
      quadCells[7][N*N - 1],
      quadCells[6][N*N - 1],
      quadCells[5][N*N - 1],
    ];
  }

  calculateVertices(grid) {
    const neighborCells = this.neighbors(grid);
    this.vertices = new Array(neighborCells.length);
    for (let i = 0; i < neighborCells.length; i++) {
      const neighborA = neighborCells[i];
      const neighborB = neighborCells[(i+1)%neighborCells.length];
      const vertex = circumcenter(this.centerVertex, neighborA.centerVertex, neighborB.centerVertex);
      this.vertices[i] = vertex;
    }
    this.calculateFaceTriangles();
  }

  calculateFaceTriangles() {
    // NOTE: This only works correctly if the vertices array
    // is in a counter-clockwise order around this cell.
    this.faceTriangles = new Array(this.vertices.length-2);
    let tri = 0;
    this.faceTriangles[tri++] = new THREE.Triangle(
      this.vertices[0],
      this.vertices[1],
      this.vertices[2],
    );
    for(let i = 2; i+1 < this.vertices.length; i++) {
      this.faceTriangles[tri++] = new THREE.Triangle(
        this.vertices[0],
        this.vertices[i],
        this.vertices[i+1],
      );
    }
    for (let triangle of this.faceTriangles) {
      this.area += triangle.getArea();
    }
  }

  /**
   * Returns the average angle (in radians) between this cell's center vertex
   * and all of its neighbors' center vertices.
   * Assumes vertices are unit vectors on a sphere.
   */
  averageNeighborAngle() {
    const p = this.centerVertex;
    let total = 0;

    for (const n of this.neighbors) {
      // Clamp dot product to avoid numerical errors outside [-1, 1]
      const dot = Math.max(-1, Math.min(1, p.dot(n.centerVertex)));
      total += Math.acos(dot);
    }

    return total / this.neighbors.length;
  }
}

class Grid {
  constructor(N) {
    let start = performance.now()
    const { quadCells, northPole, southPole } = Grid.make(N);
    console.log(`Grid make completed in ${performance.now()-start}ms`)

    this.quadCells = quadCells;
    this.northPole = northPole;
    this.southPole = southPole;
    this.size = 10*N*N + 2

    for (const cell of this) {
      cell.calculateVertices(this);
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

      const cells = new Array(N*N);
      // The refGrid is (N+1)^2 in order to span all edges
      // of the triangle faces. The points we care about
      // range from x=1 to x=N and from y=0 to Y=N-1.
      for (let x = 1; x <= N; x++) {
        for (let y = 0; y < N; y++) {
          let refPoint, refFrame;
          if (x+y <= N) {
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
          cells[N*(x-1)+y] = new GridCell(N, i, x-1, y, center);
        }
      }
      quadCells[i] = cells;
    }

    const northPoleCenter = ISEA.unprojectFromFace(refGrid[0][0].up, ico.quads[0].up.frame);
    const southPoleCenter = ISEA.unprojectFromFace(refGrid[N][N].down, ico.quads[5].down.frame);

    const northPole = new GridCell(N, -1, 0, 0, northPoleCenter);
    const southPole = new GridCell(N, -2, 0, 0, southPoleCenter);

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

    const grid = new Array(N+1);
    for (let col = 0; col <= N; col++) {
      grid[col] = new Array(N+1);
      for (let row =0; row <= N; row++) {
        grid[col][row] = {};
      }
    }

    for (let x = 0; x <= N; x++) {
      // Only traverse up to x+y <= N.
      // This keeps the loop within the face of triangle ABC.
      for (let y = 0; x+y <= N; y++) {
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
        grid[N-x][N-y].down = p;
      }
    }

    return grid;
  }
}

export {
  Grid, GridCell, circumcenter, centroid
};
