// icosahedral-snyder-equal-area.js
// ES module implementing Snyder's equal-area projection on an icosahedron.
// Now builds frames from your custom Icosahedron class (not THREE.IcosahedronGeometry).

import * as THREE from "./three.module.js";

/**
 * Snyder equal-area constants for the *icosahedron* (radians)
 */
const DEG = Math.PI / 180;
const RAD = 1 / DEG;
const g  = 37.37736814 * DEG;
const G  = 36.0        * DEG;
const th = 30.0        * DEG;

const sin_g = Math.sin(g);
const cos_g = Math.cos(g);
const tan_g = Math.tan(g);
const tan2_g = tan_g ** 2;

const sin_G = Math.sin(G);
const cos_G = Math.cos(G);

const sin_th = Math.sin(th);
const cos_th = Math.cos(th);
const tan_th = Math.tan(th);
const cot_th = 1 / tan_th;


const SECTOR = Math.PI * 2/3; // 120°

const EPS        = 1e-9;
const NEWTON_MAX = 8;

// Roughly 0.91038
// (G-th)*R^2 = R'^2 tan(g)^2 sin(th)cos(th)/2
// R'^2 = (G-th)*R^2 / (tan(g)^2 sin(th)cos(th)/2)
// R' = R √((G-th) / (tan(g)^2 sin(th)cos(th)/2))
const RpOverR = Math.sqrt(2 * (G-th) / (tan2_g*sin_th*cos_th))
const RpOverR_squared = RpOverR ** 2;

function roundIfNearZero(x) {
  return Math.abs(x) < EPS ? 0 : x;
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
export function triangleFrame({a, b, c}) {
  // 1) face-center direction (spherical centroid)
  const n = a.clone().add(b).add(c).normalize();
  // 2) base direction: b -> c, then Gram–Schmidt into tangent plane
  const ref = c.clone().sub(b);
  // (edge is already parallel to tangent; this projection zeros any 1e-15 drift)
  const u = ref.clone().sub(n.clone().multiplyScalar(ref.dot(n))).normalize();
  // 3) complete right-handed basis
  const v = new THREE.Vector3().crossVectors(n, u).normalize();

  return { n, u, v };
};

function foldAzimuth(Az_) {
  let k = 0;
  let Az = Az_;
  while (Az > SECTOR) {
    Az -= SECTOR; k++;
  }
  return [k, Az];
}

function divmod(x, y) {
  if (y === 0) {
    throw new Error("Division by zero is not allowed.");
  }
  const quotient = Math.floor(x / y);
  const remainder = x % y;
  return [quotient, remainder];
};

/** Convert 3D unit vector -> (lon, lat) in radians. */
function xyzToLonLat(p) {
  const lon = Math.atan2(p.y, p.x);
  const rxy = Math.hypot(p.x, p.y);
  const lat = Math.atan2(p.z, rxy);
  return { lon, lat };
}

/** Convert (lon, lat) in radians -> 3D unit vector. */
function lonLatToXyz(lon, lat) {
  const cl = Math.cos(lat);
  return new THREE.Vector3(
    Math.cos(lon) * cl,
    Math.sin(lon) * cl,
    Math.sin(lat)
  );
}

/** azimuthFromApex
 ** Returns the angle in the plane (in radians) between the positive
 ** y-axis and the ray from (0, 0) to the point (x, y) ranging between
 ** 0 and 2π in the clockwise direction.
 **/
function azimuthFromApex(x, y) {
  return Math.PI - Math.atan2(x, -y);
}

function q_of_Az(Az, cos_Az, sin_Az) {
  cos_Az = cos_Az === undefined ? Math.cos(Az) : cos_Az;
  sin_Az = sin_Az === undefined ? Math.sin(Az) : sin_Az;

  return Math.atan2(tan_g, cos_Az + sin_Az*cot_th);
}

function H_of_Az(Az, cos_Az, sin_Az) {
  cos_Az = cos_Az === undefined ? Math.cos(Az) : cos_Az;
  sin_Az = sin_Az === undefined ? Math.sin(Az) : sin_Az;

  const val = sin_Az*sin_G*cos_g - cos_Az*cos_G;
  return Math.acos(THREE.MathUtils.clamp(val, -1, 1));
}

function F_of_Az_AG_H(Az, AG, H) {
  return Math.PI + AG - G - H - Az;
}

function Fprime_of_AZ_H(Az, H, cos_Az, sin_Az) {
  cos_Az = cos_Az === undefined ? Math.cos(Az) : cos_Az;
  sin_Az = sin_Az === undefined ? Math.sin(Az) : sin_Az;

  return ((cos_Az*sin_G*cos_g + sin_Az*cos_G) / Math.sin(H) ) - 1;
}

function delta_Az(Az, AG, cos_Az, sin_Az) {
  cos_Az = cos_Az === undefined ? Math.cos(Az) : cos_Az;
  sin_Az = sin_Az === undefined ? Math.sin(Az) : sin_Az;

  const H = H_of_Az(Az, cos_Az, sin_Az);

  return - F_of_Az_AG_H(Az, AG, H) / Fprime_of_AZ_H(Az, H, cos_Az, sin_Az);
}

function dp_of_Azp(Azp) {
  return RpOverR*tan_g / (Math.cos(Azp) + Math.sin(Azp)*cot_th);
}

function f_of_dp_q(dp, q) {
  return dp / (2 * RpOverR * Math.sin(q/2));
}

function rho_of_f_z(f, z) {
  return 2 * RpOverR * f * Math.sin(z/2);
}

function z_of_rho_f(rho, f) {
  const arg = rho / (2 * RpOverR * f);
  return 2 * Math.asin(THREE.MathUtils.clamp(arg, -1, 1));
}

/**
 * Forward: unit Vector3 P → { x, y }
 */
export function projectVectorToFace(P, frame) {
  const { n, u, v } = frame;

  // z is the spherical angle beween P and n.
  const dot = THREE.MathUtils.clamp(n.dot(P), -1, 1);
  const z   = Math.acos(dot);

  // If z is greater than g, then P is not in the given frame as it is
  // further away in angle than any vertex on this face.
  if (z > (g+EPS)) {
    console.log(`Warning - z (${z * RAD} degrees) exceeds g (${g * RAD} degrees) by ${(z-g) * RAD} degrees.`);
  }

  // t is the projection of P into the normal plane.
  // We use this gnomonic projection to find spherical azimuth,
  // Az by breaking t into its components of u and v.
  const t   = P.clone().sub(n.clone().multiplyScalar(dot)).normalize();
  const Az_ = azimuthFromApex(t.dot(u), t.dot(v));

  // The following code assumes an azimuth between 0 and SECTOR so we fold
  // Az_ into that range but keep track of the index so we can rotate back
  // later.
  const [k, Az] = foldAzimuth(Az_);
  const cos_Az = Math.cos(Az);
  const sin_Az = Math.sin(Az);

  // q is the sphereical angle between P and a point on the edge of the face
  // along the azimuth.
  const q = q_of_Az(Az, cos_Az, sin_Az);

  // If z is greater than q, then P is not in the given frame as it is beyond
  // the edge of the triangle face (though not as far as a vertex is, g).
  if (z > (q+EPS)) {
    console.log(`Warning - z (${z * RAD} degrees) exceeds q (${q * RAD} degrees) by ${(z-q) * RAD} degrees.`);
  }

  // H is the internal angle on the spherical triangle A'B'D' along the edge at D'.
  const H = H_of_Az(Az, cos_Az, sin_Az);

  // AG is the Area of the triangle A'B'D'.
  const AG = Az + G + H - Math.PI;
  // Azp is the azimuth on our flat projection.
  const Azp = Math.atan2(2*AG, RpOverR_squared*tan2_g - 2*AG*cot_th);
  // dp is the point along the flat edge BC that is Azp from A.
  const dp = dp_of_Azp(Azp);
  // f is a scaling factor which varies with azimuth to maintain equal area.
  const f = f_of_dp_q(dp, q);
  // rho is the length of the flattened arc to point P in our projection.
  const rho = 2 * RpOverR * f * Math.sin(z/2);

  // Now rotate the azimuth back to the correct sector.
  const Azp_ = Azp + k * SECTOR;

  const x = rho * Math.sin(Azp_);
  const y = rho * Math.cos(Azp_);

  // console.log({
  //   z: z*RAD,
  //   Az_: Az_*RAD,
  //   Az: Az*RAD,
  //   k,
  //   q: q*RAD,
  //   H,
  //   AG,
  //   Azp: Azp*RAD,
  //   dp,
  //   f,
  //   rho,
  //   Azp_: Azp_*RAD,
  //   x,
  //   y
  // });

  // return {
  //   x: roundIfNearZero(x),
  //   y: roundIfNearZero(y)
  // };
  return { x, y };
}

/**
 * Inverse: (x, y) → P (Vector3)
 */
export function unprojectFromFace({x, y}, frame, trace) {
  const { n, u, v } = frame;

  // x = roundIfNearZero(x);
  // y = roundIfNearZero(y);

  const Azp_ = azimuthFromApex(x, y);
  const rho = Math.hypot(x, y);

  // The following code assumes an azimuth between 0 and SECTOR so we fold
  // Azp_ into that range but keep track of the index so we can rotate back
  // later.
  const [k, Azp] = foldAzimuth(Azp_);

  if (trace) {
    console.log({ n, u, v, x, y, Azp_: Azp_*RAD, rho, k, Azp: Azp*RAD});
  }

  // Calculate the area AG using Azp.
  const cos_AZp = Math.cos(Azp);
  const sin_Azp = Math.sin(Azp);
  const AG_denom = 2 * (cos_AZp + cot_th*sin_Azp); // Never near zero in our range of Azp.
  const AG = RpOverR_squared*tan2_g*sin_Azp/AG_denom;

  let Az = Azp;
  let cos_Az = cos_AZp;
  let sin_Az = sin_Azp;

  let iter = 0;
  while (iter++ < NEWTON_MAX) {
    const dAz = delta_Az(Az, AG, cos_Az, sin_Az);

    Az += dAz;
    cos_Az = Math.cos(Az);
    sin_Az = Math.sin(Az);

    if (Math.abs(dAz) < EPS) break;
  }
  if (iter == NEWTON_MAX) {
    console.log(`Azp=${Azp*RAD} converged to Az=${Az*RAD} in ${iter} iterations.`);
  }

  const q = q_of_Az(Az, cos_Az, sin_Az);
  const dp = dp_of_Azp(Azp);
  const f = f_of_dp_q(dp, q);
  const z = z_of_rho_f(rho, f);

  // Unfold azimuth.
  const Az_ = Az + k * SECTOR;

  // Reconstruct the 3D point.
  const t = u.clone().multiplyScalar(Math.sin(Az_))
    .add(v.clone().multiplyScalar(Math.cos(Az_)))
    .normalize();
  const P = n.clone().multiplyScalar(Math.cos(z))
    .add(t.multiplyScalar(Math.sin(z)))
    .normalize();

  // P.x = roundIfNearZero(P.x);
  // P.y = roundIfNearZero(P.y);
  // P.z = roundIfNearZero(P.z);

  return P;
}


/** Convenience wrappers (lon/lat in radians). */
export function projectLonLat(lon, lat, frame) {
  const P = lonLatToXyz(lon, lat);
  return projectVectorToFace(P, frame);
}

export function unprojectToLonLat(x, y, frame) {
  const P = unprojectFromFace(x, y, frame);
  return xyzToLonLat(P);
}

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

  static makeQuad(upTri, downTri) {
    const upFrame = triangleFrame(upTri);
    const downFrame = triangleFrame(downTri);

    return {
      up: { tri: upTri, frame: upFrame },
      down: { tri: downTri, frame: downFrame },
    };
  }
};
