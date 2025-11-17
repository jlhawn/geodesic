import * as THREE from "./three.module.js";

export class Icosahedron {
  // Unfolded Icosahedron Diagram:
  //  Vertices labeled A through L
  //  Faced numbered 0 through 19
  //
  //        A             A             A             A             A
  //       / \           / \           / \           / \           / \
  //      /   \         /   \         /   \         /   \         /   \
  //     /     \       /     \       /     \       /     \       /     \
  //    /       \     /       \     /       \     /       \     /       \
  //   /    0    \   /    1    \   /    2    \   /    3    \   /    4    \
  //  /           \ /           \ /           \ /           \ /           \
  // B - - - - - - C - - - - - - D - - - - - - E - - - - - - F - - - - - - B
  //  \           / \           / \           / \           / \           / \
  //   \         /   \         /   \         /   \         /   \         /   \
  //    \   5   /     \   6   /     \   7   /     \   8   /     \   9   /     \
  //     \     /       \     /       \     /       \     /       \     /       \
  //      \   /   10    \   /   11    \   /   12    \   /   13    \   /   14    \
  //       \ /           \ /           \ /           \ /           \ /           \
  //        G - - - - - - H - - - - - - I - - - - - - J - - - - - - K - - - - - - G
  //         \           / \           / \           / \           / \           /
  //          \         /   \         /   \         /   \         /   \         /
  //           \  15   /     \  16   /     \  17   /     \  18   /     \  19   /
  //            \     /       \     /       \     /       \     /       \     /
  //             \   /         \   /         \   /         \   /         \   /
  //              \ /           \ /           \ /           \ /           \ /
  //               L             L             L             L             L
  //
  constructor() {
    // The simplist costruction is with a = 1 and b = (1+√5)/2
    // however we want all vectors to begin normalized so we
    // scale those values so that the magnitude of <a,b> is 1.
    const a = 1 / Math.sqrt((5 + Math.sqrt(5))/2)
    const b = a * (1 + Math.sqrt(5))/2;

    // 12 Vertices of an Icosahedron.

    const A = new THREE.Vector3(-a,  0,  b);
    const B = new THREE.Vector3( a,  0,  b);
    const C = new THREE.Vector3( 0,  b,  a);
    const D = new THREE.Vector3(-b,  a,  0);
    const E = new THREE.Vector3(-b, -a,  0);
    const F = new THREE.Vector3( 0, -b,  a);
    const G = new THREE.Vector3( b,  a,  0);
    const H = new THREE.Vector3( 0,  b, -a);
    const I = new THREE.Vector3(-a,  0, -b);
    const J = new THREE.Vector3( 0, -b, -a);
    const K = new THREE.Vector3( b, -a,  0);
    const L = new THREE.Vector3( a,  0, -b);
    

    this.vertices = [A, B, C, D, E, F, G, H, I, J, K, L];

    /*
    These vertices start out pretty symmetrical, each being within a plane
    of two of the x, y, or z axes. We want to rotate it about the y axis
    such that A becomes a unit vector in the direction of the z axis (and
    L becomes a unit vector in the negatize-z direction). The sine of that
    rotation angle happens to be 'a' and the cosine happens to be 'b'. The
    transformation matrix for this rotation is given by:

      [ b   0  a ]
      [ 0   1  0 ]
      [ -a  0  b ]

    After applying this rotation, A is the "north pole" and L is the "south
    pole". This is a bit more straightforward than using a quaternion to
    rotate the vectors.

    TODO: support options for having the pole centered on an edge or face
    instead of a vertex.
    */
    const m = new THREE.Matrix3(
      b, 0, a,
      0, 1, 0,
      -a, 0, b,
    );
    this.vertices.forEach(function(vec) {
      vec.applyMatrix3(m);
    });

    // Each triangle needs to specify its vertices in counter-clockwise
    // (CCW) order. We also use the convention that if the triangle
    // points up (in the Unfolded Icosahedron Diagram above) then we
    // start with the vertex at the top, while if it points down then we
    // start with the vertex at the bottom.
    this.faces = [
      new THREE.Triangle(A, B, C),
      new THREE.Triangle(A, C, D),
      new THREE.Triangle(A, D, E),
      new THREE.Triangle(A, E, F),
      new THREE.Triangle(A, F, B),

      new THREE.Triangle(G, C, B),
      new THREE.Triangle(H, D, C),
      new THREE.Triangle(I, E, D),
      new THREE.Triangle(J, F, E),
      new THREE.Triangle(K, B, F),

      new THREE.Triangle(C, G, H),
      new THREE.Triangle(D, H, I),
      new THREE.Triangle(E, I, J),
      new THREE.Triangle(F, J, K),
      new THREE.Triangle(B, K, G),

      new THREE.Triangle(L, H, G),
      new THREE.Triangle(L, I, H),
      new THREE.Triangle(L, J, I),
      new THREE.Triangle(L, K, J),
      new THREE.Triangle(L, G, K),
    ];

    this.quads = [
      Icosahedron.makeQuad(this.faces[0], this.faces[5]),
      Icosahedron.makeQuad(this.faces[1], this.faces[6]),
      Icosahedron.makeQuad(this.faces[2], this.faces[7]),
      Icosahedron.makeQuad(this.faces[3], this.faces[8]),
      Icosahedron.makeQuad(this.faces[4], this.faces[9]),
      Icosahedron.makeQuad(this.faces[10], this.faces[15]),
      Icosahedron.makeQuad(this.faces[11], this.faces[16]),
      Icosahedron.makeQuad(this.faces[12], this.faces[17]),
      Icosahedron.makeQuad(this.faces[13], this.faces[18]),
      Icosahedron.makeQuad(this.faces[14], this.faces[19]),
    ]
  }

  /**
   *            A       
   *           / \      
   *          /   \     
   *         /     \    
   *        /   .   \   
   *       /         \  
   *  v↑  /           \ 
   *     B - - - - - - C
   *  n⊙   u→
   * Build a per-face frame:
   *  - n: unit face-center direction (spherical centroid)
   *  - u: x-axis in tangent plane (projected b→c edge)
   *  - v: y-axis in tangent plane, v = n × u (right-handed)
   */
  static triangleFrame({a, b, c}) {
    // 1) face-center direction (spherical centroid)
    const n = a.clone().add(b).add(c).normalize();
    // 2) base direction: b -> c, then Gram–Schmidt into tangent plane
    const ref = c.clone().sub(b);
    // (edge is already parallel to tangent; this projection zeros any 1e-15 drift)
    const u = ref.clone().sub(n.clone().multiplyScalar(ref.dot(n))).normalize();
    // 3) complete right-handed basis
    const v = new THREE.Vector3().crossVectors(n, u).normalize();

    return { n, u, v };
  }

  static makeQuad(upTri, downTri) {
    const upFrame = Icosahedron.triangleFrame(upTri);
    const downFrame = Icosahedron.triangleFrame(downTri);

    return {
      up: { tri: upTri, frame: upFrame },
      down: { tri: downTri, frame: downFrame },
    }
  }
};
