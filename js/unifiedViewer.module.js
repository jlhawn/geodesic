import * as THREE from "./three.module.js";
import { twoFingerMotion } from "./gestures.module.js";

// ----------------------------------------------------------------------------
// SHADER DEFINITIONS
// ----------------------------------------------------------------------------

const vertexHead = `
uniform mat4 uModelRotation; 
uniform sampler2D uCenterTexture;
uniform vec2 uTexSize; 
uniform float uBlend; // 0.0 = Sphere, 1.0 = Map
uniform float uMapLift;

attribute float cellIndex; 

// Projection Constants
const float A1 = 1.340264;
const float A2 = -0.081106;
const float A3 = 0.000893;
const float A4 = 0.003796;
const float SQRT3 = 1.73205080757;
const float M_SQRT3_2 = 0.86602540378;
const float PI = 3.14159265359;

vec2 projectEqualEarth(float lat, float lon) {
  float sinPhi = sin(lat);
  float theta = asin(M_SQRT3_2 * sinPhi);
  
  float theta2 = theta * theta;
  float theta6 = theta2 * theta2 * theta2;
  float theta7 = theta6 * theta;
  float theta8 = theta6 * theta2;
  float theta9 = theta8 * theta;

  float cosTheta = cos(theta);
  float denom = 3.0 * (9.0 * A4 * theta8 + 7.0 * A3 * theta6 + 3.0 * A2 * theta2 + A1);
  
  float projX = (2.0 * SQRT3 * lon * cosTheta) / denom;
  float projY = A4 * theta9 + A3 * theta7 + A2 * theta2 * theta + A1 * theta;

  return vec2(projX, projY);
}
`;

const vertexLogic = `
  // 1. Data Lookup
  float col = mod(cellIndex, uTexSize.x);
  float row = floor(cellIndex / uTexSize.x);
  vec2 uv = (vec2(col, row) + 0.5) / uTexSize; 
  vec3 cellCenter = texture2D(uCenterTexture, uv).xyz;
  vec3 sourcePos = position;
  // [vertex source]

  // 2. Apply Rotation
  vec4 rotatedPos = uModelRotation * vec4(sourcePos, 1.0);
  vec4 rotatedCenter = uModelRotation * vec4(cellCenter, 1.0);
  
  // 3. SPHERE POSITION
  vec3 posSphere = rotatedPos.xyz;

  // 4. MAP POSITION
  // We assume the sphere is rotated such that Y is Up (North) and Z is Front.
  // Lat: comes from Y height (-1 to 1)
  // Lon: comes from angle in XZ plane (atan(x, z))
  
  float lat = asin(clamp(rotatedPos.y, -1.0, 1.0));
  float lon = atan(rotatedPos.x, rotatedPos.z);
  float centerLon = atan(rotatedCenter.x, rotatedCenter.z);

  // Fix Tearing
  float delta = lon - centerLon;
  if (delta > PI) lon -= 2.0 * PI;
  if (delta < -PI) lon += 2.0 * PI;

  // Project
  vec2 proj = projectEqualEarth(lat, lon);
  vec3 posMap = vec3(proj, uMapLift);

  // 5. MORPH
  vec3 finalPos = mix(posSphere, posMap, uBlend);
  vec3 transformed = finalPos;
`;

const EARTH_RADIUS = 6.371e6, VERTICAL_EXAGGERATION = 1, METRE = VERTICAL_EXAGGERATION / EARTH_RADIUS;
const SCALE_HEIGHT = 8000 * METRE, AEROSOL_HEIGHT = 1500 * METRE, OZONE_LOW = 1 + 15000 * METRE, OZONE_HIGH = 1 + 35000 * METRE;
const RAYLEIGH = [0.06, 0.12, 0.29], AEROSOL = [0.05, 0.06, 0.07], OZONE = [0.035, 0.025, 0.003];
const CLOUD_ROUGHNESS = 0.2, CLOUD_TINT = 2.5, CLOUD_WARMTH_RISE = 1, DAY_START = 0.05, DAY_FULL = 0.175, NIGHT_DEEP = -0.35, CAMERA_TOE = 0.08, TOE_SATURATION = 0.5, CLOUD_RELIEF = 20;
const glsl = (x) => { const s = Number(x).toPrecision(9); return /[.e]/.test(s) ? s : `${s}.0`; };
const glslVec3 = (v) => `vec3(${v.map(glsl).join(', ')})`;

/*
 * Sunlight through the air, in the globe's units (radius 1). Every height
 * the Satellite view draws is VERTICAL_EXAGGERATION (two) times its real
 * value, the air's scale height included, so cloud relief, cloud shadows,
 * how far past the terminator a high top stays lit and the layers of the
 * limb follow one exaggeration (METRE converts real metres to globe radii).
 * The depths are sea-level vertical optical depths in linear R, G, B of
 * Rayleigh air, a boundary-layer aerosol, and the ozone layer's Chappuis
 * band as a shell between 15 and 35 km. Light is in units of a white
 * surface under the overhead sun.
 */
const atmosphere = `
const float SCALE_HEIGHT = ${glsl(SCALE_HEIGHT)};
const float AEROSOL_HEIGHT = ${glsl(AEROSOL_HEIGHT)};
const float OZONE_LOW = ${glsl(OZONE_LOW)};
const float OZONE_HIGH = ${glsl(OZONE_HIGH)};
const float METRE = ${glsl(METRE)};
const float CLOUD_ROUGHNESS = ${glsl(CLOUD_ROUGHNESS)};
const float CLOUD_TINT = ${glsl(CLOUD_TINT)};
const float CLOUD_WARMTH_RISE = ${glsl(CLOUD_WARMTH_RISE)};
const float CAMERA_TOE = ${glsl(CAMERA_TOE)};
const float TOE_SATURATION = ${glsl(TOE_SATURATION)};
const float AIR_NODES = 6.0;
const vec3 RAYLEIGH = ${glslVec3(RAYLEIGH)};
const vec3 AEROSOL = ${glslVec3(AEROSOL)};
const vec3 OZONE = ${glslVec3(OZONE)};
const vec3 BRIGHTNESS = vec3(0.2126, 0.7152, 0.0722);
const float AEROSOL_ALBEDO = 0.9, AEROSOL_FORWARD = 0.65;

// Chapman's grazing-incidence function at x = r / H: the slant column toward cosine mu >= 0 over the vertical one.
float chapman(float x, float mu) { return 1.0 / (0.65 * mu + sqrt(0.1225 * mu * mu + 0.6366 / x)); }

/*
 * Column of a species of scale height h, in its sea-level vertical columns,
 * from radius r toward cosine mu, past the horizon over the path that dips
 * to its tangent radius; -1 when the ray meets the globe. The grazing
 * function takes the globe's radius as 'curve' times its drawn ratio to h.
 */
float grazingColumn(float r, float mu, float h, float curve) {
  if (mu >= 0.0) return exp((1.0 - r) / h) * chapman(curve * r / h, mu);
  float rt = r * sqrt(1.0 - mu * mu);
  if (rt < 1.0) return -1.0;
  return 2.0 * exp((1.0 - rt) / h) * chapman(curve * rt / h, 0.0) - exp((1.0 - r) / h) * chapman(curve * r / h, -mu);
}
float slantColumn(float r, float mu, float h) { return grazingColumn(r, mu, h, 1.0); }

// Path through the ozone shell from radius r toward cosine mu, over the shell's thickness.
float ozonePath(float r, float mu) {
  float p2 = r * r * (1.0 - mu * mu);
  float top = sqrt(max(OZONE_HIGH * OZONE_HIGH - p2, 0.0)), bottom = sqrt(max(OZONE_LOW * OZONE_LOW - p2, 0.0));
  float s = r > OZONE_HIGH ? (mu < 0.0 ? 2.0 * (top - bottom) : 0.0)
    : r < OZONE_LOW ? top - bottom
    : top - r * mu - (mu < 0.0 ? 2.0 * bottom : 0.0);
  return max(s, 0.0) / (OZONE_HIGH - OZONE_LOW);
}

// Optical depth of the sunbeam reaching radius r where the sun's cosine is mu; the globe's shadow is opaque.
vec3 pathDepth(float r, float mu, float curve) {
  float air = grazingColumn(r, mu, SCALE_HEIGHT, curve);
  if (air < 0.0) return vec3(1.0e4);
  return RAYLEIGH * air + AEROSOL * max(grazingColumn(r, mu, AEROSOL_HEIGHT, curve), 0.0) + OZONE * ozonePath(r, mu);
}
vec3 sunDepth(float r, float mu) { return pathDepth(r, mu, 1.0); }

/*
 * The same for the beam that lights a surface or a cloud top, its Rayleigh
 * and aerosol columns grazing at the real ratio of the globe's radius to the
 * scale height: about sqrt(pi R / 2 H) = 35 vertical columns at sea level,
 * which the drawn air, thicker on the same globe, would cut by the square
 * root of the exaggeration. Where the globe's shadow falls and how far below
 * a top the beam dips stay in the drawn heights.
 */
vec3 beamDepth(float r, float mu) { return pathDepth(r, mu, ${glsl(VERTICAL_EXAGGERATION)}); }

vec3 zenithDepth(float r) {
  return RAYLEIGH * exp((1.0 - r) / SCALE_HEIGHT) + AEROSOL * exp((1.0 - r) / AEROSOL_HEIGHT) + OZONE * ozonePath(r, 1.0);
}

// The beam's colour on a surface at radius r, white where the sun stands overhead.
vec3 sunColour(float r, float mu) { return exp(zenithDepth(r) - beamDepth(r, mu)); }

// Lambert's cosine averaged over facets whose slope toward the sun spreads by sigma radians.
float lambertSoft(float c, float sigma) {
  float s = 0.576 * sigma;
  return s > 0.0 ? max(c, 0.0) + s * log(1.0 + exp(-abs(c) / s)) : max(c, 0.0);
}

/*
 * Sky light on a level surface at radius r with the sun's cosine mu, in
 * units of the sun's flux: the sunbeam scattered once by the air above
 * that the globe's shadow leaves lit, taken at four heights up through it,
 * each slab lit by its own beam (deep red where that beam grazed the
 * ground, bluer where it only crossed the ozone layer) and sending half of
 * what it scatters down through the Rayleigh air and aerosol below at the
 * diffuse air mass 1.66. After sunset the shadow climbs, the lit air thins
 * and sees a redder beam, so the light goes orange, pink and purple as it
 * fades, and is gone once the shadow stands twelve scale heights up.
 */
vec3 skyLight(float r, float mu) {
  float rb = mu < 0.0 ? max(r, inversesqrt(max(1.0 - mu * mu, 1.0e-6))) : r;
  float hb = (rb - 1.0) / SCALE_HEIGHT;
  if (hb > 12.0) return vec3(0.0);
  float hr = (r - 1.0) / SCALE_HEIGHT, ar = (r - 1.0) / AEROSOL_HEIGHT;
  const float NODE_H[4] = float[4](0.2, 0.7, 1.6, 3.5);
  const float NODE_W[4] = float[4](0.4, 0.6, 1.2, 2.6);
  vec3 sum = vec3(0.0);
  for (int k = 0; k < 4; k++) {
    float hk = hb + NODE_H[k], rk = 1.0 + hk * SCALE_HEIGHT;
    vec3 beam = exp(-beamDepth(rk, mu));
    vec3 slab = RAYLEIGH * exp(-hk) * NODE_W[k];
    vec3 below = RAYLEIGH * max(exp(-hr) - exp(-hk), 0.0) + AEROSOL * max(exp(-ar) - exp(-hk * SCALE_HEIGHT / AEROSOL_HEIGHT), 0.0);
    sum += beam * slab * exp(-1.66 * below);
  }
  return 0.5 * sum;
}

/*
 * Sun and sky on an element at radius r with facet cosine 'facet', relative
 * to the same element under the overhead sun; 'unblocked' is the share of
 * the beam that nothing sunward casts off. The beam's hue is the colour the
 * beam has with the sun 'rise' times as high above the horizon (its
 * brightness stays the beam's own), raised to the power 'tint' at the same
 * brightness, which deepens a reddened beam and leaves a white one white.
 */
vec3 illumination(float r, float mu, float facet, float roughness, float unblocked, float tint, float rise) {
  vec3 zenith = zenithDepth(r), overhead = exp(-zenith);
  vec3 colour = exp(zenith - beamDepth(r, mu));
  vec3 hue = rise == 1.0 ? colour : exp(zenith - beamDepth(r, mu > 0.0 ? mu * rise : mu));
  vec3 deeper = pow(hue, vec3(tint));
  colour = deeper * dot(colour, BRIGHTNESS) / max(dot(deeper, BRIGHTNESS), 1.0e-12);
  vec3 beam = overhead * colour * lambertSoft(facet, roughness) * unblocked;
  return (beam + skyLight(r, mu)) / (overhead + skyLight(r, 1.0));
}

/*
 * Sunlight scattered once toward the camera by the air along a ray from o
 * along the unit d, out to space: the ray's Rayleigh column cut into 'nodes'
 * equal parts, each lit as at its middle (the sunbeam's own depth there, the
 * globe's shadow included) and dimmed by the Rayleigh air between it and the
 * camera, integrated exactly over the part, and by the aerosol and ozone,
 * plus the aerosol of each part (its density against the air's at the
 * part's middle) scattering the beam forward in a Henyey-Greenstein lobe,
 * which is the glow around a setting sun. The camera is beyond the far end,
 * or, with 'back', behind o (the far half of a limb ray). 'depth' is the
 * whole ray's depth.
 */
vec3 airAlong(vec3 o, vec3 d, vec3 sun, float nodes, bool back, out vec3 depth) {
  float r0 = length(o), mud = max(dot(o, d) / r0, 0.0), b = r0 * mud;
  float column = exp((1.0 - r0) / SCALE_HEIGHT) * chapman(r0 / SCALE_HEIGHT, mud);
  float aerosol0 = max(slantColumn(r0, mud, AEROSOL_HEIGHT), 0.0), ozone0 = ozonePath(r0, mud);
  vec3 total = RAYLEIGH * column;
  depth = total + AEROSOL * aerosol0 + OZONE * ozone0;
  float cosT = dot(d, sun);
  vec3 light = vec3(0.0), glow = vec3(0.0);
  for (int k = 0; k < 8; k++) {
    if (float(k) >= nodes) break;
    float rk = r0 - SCALE_HEIGHT * log(1.0 - (float(k) + 0.5) / nodes);
    vec3 p = o + (sqrt(max(b * b + rk * rk - r0 * r0, 0.0)) - b) * d;
    float r = length(p);
    vec3 up = p / r;
    float mk = dot(up, d);
    float a = max(slantColumn(r, mk, AEROSOL_HEIGHT), 0.0), oz = ozonePath(r, mk);
    vec3 view = back ? AEROSOL * max(aerosol0 - a, 0.0) + OZONE * max(ozone0 - oz, 0.0) : AEROSOL * a + OZONE * oz;
    vec3 lo = total * (float(k) / nodes), hi = total * ((float(k) + 1.0) / nodes);
    vec3 share = back ? exp(-lo) - exp(-hi) : exp(hi - total) - exp(lo - total);
    vec3 lit = exp(-sunDepth(r, dot(up, sun)) - view) * share;
    light += lit;
    glow += lit * (AEROSOL / AEROSOL_HEIGHT * exp((1.0 - r) / AEROSOL_HEIGHT)) / (RAYLEIGH / SCALE_HEIGHT * exp((1.0 - r) / SCALE_HEIGHT));
  }
  float g = AEROSOL_FORWARD, lobe = 0.25 * (1.0 - g * g) / pow(1.0 + g * g - 2.0 * g * cosT, 1.5);
  return 0.1875 * (1.0 + cosT * cosT) * light + AEROSOL_ALBEDO * lobe * glow;
}

// The air along a ray past the globe through tangent point t, direction v away from the camera: the near half, and the far half seen through it.
vec3 limbLight(vec3 t, vec3 v, vec3 sun) {
  vec3 nearDepth, farDepth;
  vec3 near = airAlong(t, -v, sun, 4.0, false, nearDepth);
  vec3 far = airAlong(t, v, sun, 4.0, true, farDepth);
  return near + exp(-nearDepth) * far;
}

// The air light already in the colours the page paints: looking straight down with the sun overhead.
vec3 bakedAir(float r) { return 0.1875 * (1.0 - exp(-2.0 * RAYLEIGH * exp((1.0 - r) / SCALE_HEIGHT))); }

// What the slant path to the camera takes from a surface's light beyond the straight-down path: aerosol and ozone absorb; of the light Rayleigh
// scattering takes out of the view, half comes back from the surroundings, taken to have the surface's own colour.
vec3 viewLoss(float r, float muv) {
  float air = exp((1.0 - r) / SCALE_HEIGHT) * (chapman(r / SCALE_HEIGHT, muv) - chapman(r / SCALE_HEIGHT, 1.0));
  float aerosol = exp((1.0 - r) / AEROSOL_HEIGHT) * (chapman(r / AEROSOL_HEIGHT, muv) - chapman(r / AEROSOL_HEIGHT, 1.0));
  return 0.5 * (1.0 + exp(-RAYLEIGH * air)) * exp(-AEROSOL * aerosol - OZONE * (ozonePath(r, muv) - ozonePath(r, 1.0)));
}

/*
 * The day side's own look, which the twilight above hands over to as the
 * sun climbs from DAY_START to DAY_FULL in cosine: the sunbeam reddened over
 * the air mass of a sea-level path, a flat blue sky light, blue air seen at
 * a slant toward the limb (DAY_AIR) and the blue rim past it (DAY_RIM),
 * which fades out through the terminator (dayRim).
 */
const vec3 DAY_SKY = vec3(0.1, 0.14, 0.2);
const vec3 DAY_AIR = vec3(0.45, 0.65, 1.0);
const vec3 DAY_RIM = vec3(0.1, 0.3, 1.0);
float daylight(float mu) { return smoothstep(${glsl(DAY_START)}, ${glsl(DAY_FULL)}, mu); }
float dayRim(float mu) { return smoothstep(-0.08, 0.025, mu); }
vec3 dayBeam(float mu) {
  float sunHeight = max(mu, 0.0);
  float low = min(sunHeight, 0.5);
  float airMass = 1.0 / (sunHeight + low * (1.0 - 2.0 * low) * (1.0 - 2.0 * low) + 0.025);
  return exp(-RAYLEIGH * (airMass - 1.0));
}
`;

// ----------------------------------------------------------------------------
// CPU MATH HELPER: INVERSE PROJECTION
// ----------------------------------------------------------------------------

const A1 = 1.340264;
const A2 = -0.081106;
const A3 = 0.000893;
const A4 = 0.003796;
const SQRT3 = 1.73205080757;
const M_SQRT3_2 = 0.86602540378;
const PI = 3.14159265359;

function equalEarth(lat, lon, out) {
  const theta = Math.asin(M_SQRT3_2 * Math.sin(lat));
  const theta2 = theta * theta;
  const theta6 = theta2 * theta2 * theta2;
  const denom = 3 * (9 * A4 * theta6 * theta2 + 7 * A3 * theta6 + 3 * A2 * theta2 + A1);
  out[0] = (2 * SQRT3 * lon * Math.cos(theta)) / denom;
  out[1] = A4 * theta6 * theta2 * theta + A3 * theta6 * theta + A2 * theta2 * theta + A1 * theta;
  return out;
}

function getThetaFromY(y) {
  let theta = y / A1; 
  const tolerance = 1e-6;
  const maxIter = 50;
  for (let i = 0; i < maxIter; i++) {
    const theta2 = theta * theta;
    const theta6 = theta2 * theta2 * theta2;
    const f = A4 * theta6 * theta2 * theta + A3 * theta6 * theta + A2 * theta2 * theta + A1 * theta - y;
    const fp = 9 * A4 * theta6 * theta2 + 7 * A3 * theta6 + 3 * A2 * theta2 + A1;
    const delta = f / fp;
    theta -= delta;
    if (Math.abs(delta) < tolerance) break;
  }
  return theta;
}

function inverseEqualEarthToVector(x, y) {
  if (Math.abs(y) > 1.35) return null; 
  const theta = getThetaFromY(y);
  const sinTheta = Math.sin(theta);
  const val = (2 / SQRT3) * sinTheta;
  if (Math.abs(val) > 1) return null;
  const lat = Math.asin(val);
  const theta2 = theta * theta;
  const theta6 = theta2 * theta2 * theta2;
  const theta8 = theta6 * theta2;
  const denom = 3 * (9 * A4 * theta8 + 7 * A3 * theta6 + 3 * A2 * theta2 + A1);
  const cosTheta = Math.cos(theta);
  if (Math.abs(cosTheta) < 1e-6) return null;
  const lon = (x * denom) / (2 * SQRT3 * cosTheta);
  if (Math.abs(lon) > PI + 0.1) return null;
  
  const cosLat = Math.cos(lat);
  
  // FIX: Match the shader's Y-Up orientation.
  // Shader: lat=asin(y), lon=atan(x, z)
  // Therefore: x = cosLat * sin(lon), y = sin(lat), z = cosLat * cos(lon)
  return new THREE.Vector3(
    cosLat * Math.sin(lon), // X
    Math.sin(lat),          // Y (Up)
    cosLat * Math.cos(lon)  // Z (Front)
  );
}

// ----------------------------------------------------------------------------
// MAIN VIEWER
// ----------------------------------------------------------------------------

export function initUnifiedViewer(container, grid, config = {}) {
  const {
    backgroundColor = 0x111111,
    getColor = (cell) => cell.color,
    controls = true,
  } = config;

  // --- 1. Geometry Generation ---
  const pData = [], kData = [], iData = [], idxData = [], texDataArray = [], centerData = [];
  const cellVertexStart = [], cellVertexCount = [], cellRadius = [];
  const cellsOnVertex = [];
  const colorHelper = new THREE.Color();
  let vertexCounter = 0, cellCounter = 0;

  for (const cell of grid) {
    const cv = cell.centerVertex; 
    const verts = cell.vertices || []; 
    const currentCellID = cellCounter++;
    cellVertexStart.push(vertexCounter);
    cellVertexCount.push(verts.length);

    texDataArray.push(cv.x, cv.y, cv.z, 1.0);
    centerData.push(cv.x, cv.y, cv.z);
    cellRadius.push(verts.reduce((sum, v) => sum + Math.hypot(v.x - cv.x, v.y - cv.y, v.z - cv.z), 0) / Math.max(1, verts.length));
    for (const v of verts) if (v.index !== undefined) (cellsOnVertex[v.index] ??= []).push(currentCellID);

    const cVal = getColor(cell);
    if (cVal && typeof cVal === 'object' && 'r' in cVal) colorHelper.setRGB(cVal.r, cVal.g, cVal.b);
    else colorHelper.set(cVal);
    
    const r8 = Math.floor(colorHelper.r * 255);
    const g8 = Math.floor(colorHelper.g * 255);
    const b8 = Math.floor(colorHelper.b * 255);

    const baseIndex = vertexCounter;
    for (let i = 0; i < verts.length; i++) {
      const v = verts[i];
      pData.push(v.x, v.y, v.z);
      idxData.push(currentCellID); 
      kData.push(r8, g8, b8);
      vertexCounter++;
    }

    const numTriangles = verts.length - 2;
    for (let k = 1; k <= numTriangles; k++) {
      iData.push(baseIndex, baseIndex + k, baseIndex + k + 1);
    }
  }

  // --- 2. Data Texture ---
  const width = Math.ceil(Math.sqrt(cellCounter));
  const height = Math.ceil(cellCounter / width);
  const floatBuffer = new Float32Array(width * height * 4);
  floatBuffer.set(texDataArray);
  
  const centerTexture = new THREE.DataTexture(floatBuffer, width, height, THREE.RGBAFormat, THREE.FloatType);
  centerTexture.minFilter = THREE.NearestFilter;
  centerTexture.magFilter = THREE.NearestFilter;
  centerTexture.needsUpdate = true;

  /*
   * Per-cell data the page updates every frame, on the centre texture's
   * layout so that a vertex finds its cell's texel from cellIndex: colours
   * as linear RGB bytes, the lighting surface terms, and one value per
   * cell that the vertex shader colours through the colour map.
   */
  const cellTexture = (array, format, type) => {
    const texture = new THREE.DataTexture(array, width, height, format, type);
    texture.minFilter = THREE.NearestFilter;
    texture.magFilter = THREE.NearestFilter;
    texture.needsUpdate = true;
    return texture;
  };
  const cellColors = cellTexture(new Uint8Array(4 * width * height), THREE.RGBAFormat, THREE.UnsignedByteType);
  const cellSurface = cellTexture(new Float32Array(4 * width * height), THREE.RGBAFormat, THREE.FloatType);
  const cellValues = cellTexture(new Float32Array(width * height), THREE.RedFormat, THREE.FloatType);
  const cellTerrain = cellTexture(new Float32Array(width * height), THREE.RedFormat, THREE.FloatType);
  const MISSING = 1e30, MAX_STOPS = 16;
  const colorMap = {
    uColorMode: { value: 0 }, uCellColors: { value: cellColors }, uCellSurface: { value: cellSurface }, uCellValues: { value: cellValues }, uCellTerrain: { value: cellTerrain }, uTerrain: { value: 0 },
    uStops: { value: Array.from({ length: MAX_STOPS }, () => new THREE.Vector3()) }, uStopCount: { value: 2 }, uValueMap: { value: new THREE.Vector2(1, 0) },
    uMissing: { value: new THREE.Vector3() }, uFlat: { value: new THREE.Vector3() }, uCoverBase: { value: new THREE.Vector3() },
  };

  // --- 3. Scene & Buffers ---
  function disposeArray() { this.array = null; }
  const geometry = new THREE.BufferGeometry();
  geometry.setIndex(new THREE.BufferAttribute(new Uint32Array(iData), 1).onUpload(disposeArray));
  geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(pData), 3).onUpload(disposeArray));
  geometry.setAttribute('cellIndex', new THREE.BufferAttribute(new Float32Array(idxData), 1).onUpload(disposeArray));
  geometry.setAttribute('color', new THREE.BufferAttribute(new Uint8Array(kData), 3, true).onUpload(disposeArray));
  const cornerData = new Float32Array(2 * vertexCounter);
  {
    let k = 0, c = 0;
    for (const cell of grid) {
      for (const v of cell.vertices || []) {
        const others = (v.index === undefined ? [] : cellsOnVertex[v.index]).filter((other) => other !== c);
        cornerData[k++] = others[0] ?? c;
        cornerData[k++] = others[1] ?? c;
      }
      c++;
    }
  }
  geometry.setAttribute('corner', new THREE.BufferAttribute(cornerData, 2).onUpload(disposeArray));
  const slopeAttribute = new THREE.BufferAttribute(Float32Array.from(pData, (v, k) => v / Math.hypot(pData[k - (k % 3)], pData[k - (k % 3) + 1], pData[k - (k % 3) + 2])), 3);
  geometry.setAttribute('slope', slopeAttribute);
  geometry.computeBoundingSphere();

  function fillPerCell(attribute, values, size) {
    const array = attribute.array;
    for (let c = 0; c < cellCounter; c++) {
      let v = size * cellVertexStart[c];
      for (let k = 0; k < cellVertexCount[c]; k++) for (let j = 0; j < size; j++) array[v++] = values[size * c + j];
    }
    attribute.needsUpdate = true;
  }
  const updateSlopes = (normals) => fillPerCell(slopeAttribute, normals, 3);

  function updateColors(rgb) {
    const array = cellColors.image.data;
    for (let c = 0; c < cellCounter; c++) {
      array[4 * c] = rgb[3 * c]; array[4 * c + 1] = rgb[3 * c + 1]; array[4 * c + 2] = rgb[3 * c + 2];
    }
    cellColors.needsUpdate = true;
    colorMap.uColorMode.value = 1;
  }
  /*
   * The Satellite view's surface per cell: the sea's glint weight, the
   * cloud's opacity and the heights above sea level, in metres, where light
   * from above and from below first meets the cloud; 'highest' is the
   * highest of those tops, which bounds the light pass's march.
   */
  const lightPass = { dirty: true, last: -Infinity, sun: new THREE.Vector3(2, 0, 0), cloud: 0, land: 0 };
  function updateSurface(values, highest = 0) {
    cellSurface.image.data.set(values.subarray(0, 4 * cellCounter));
    cellSurface.needsUpdate = true;
    lightPass.cloud = highest;
    lightPass.dirty = true;
  }
  function updateValues(values) {
    const array = cellValues.image.data;
    for (let c = 0; c < cellCounter; c++) { const v = values[c]; array[c] = v === v ? v : MISSING; }
    cellValues.needsUpdate = true;
  }
  function updateTerrain(elevation) {
    const array = cellTerrain.image.data;
    let land = 0;
    for (let c = 0; c < cellCounter; c++) { array[c] = elevation ? elevation[c] : 0; land = Math.max(land, array[c]); }
    cellTerrain.needsUpdate = true;
    lightPass.land = land;
    lightPass.dirty = true;
  }

  /*
   * How the shader colours the cell values: 'palette' maps
   * t = value·a + b onto sRGB stops interpolated as the page's legend
   * does, with missing values (NaN) in the linear colour `missing`, or
   * with `terrain` as the grey relief of the elevation updateTerrain
   * gave, dark for the sea floor and light for the land;
   * 'cover' composites white over the sRGB `base` with opacity
   * 1 − exp(−value·a); 'flat' paints every cell the linear `color`.
   * updateColors switches back to per-cell colours.
   */
  function setColorMap({ kind, stops = null, a = 1, b = 0, missing = null, color = null, base = null, terrain = false }) {
    colorMap.uColorMode.value = { palette: 2, cover: 3, flat: 4 }[kind];
    colorMap.uTerrain.value = terrain ? 1 : 0;
    if (stops) {
      if (stops.length > MAX_STOPS) throw new Error(`at most ${MAX_STOPS} palette stops`);
      stops.forEach((stop, k) => colorMap.uStops.value[k].set(stop[0], stop[1], stop[2]));
      colorMap.uStopCount.value = stops.length;
    }
    colorMap.uValueMap.value.set(a, b);
    if (missing) colorMap.uMissing.value.set(missing[0], missing[1], missing[2]);
    if (color) colorMap.uFlat.value.set(color[0], color[1], color[2]);
    if (base) colorMap.uCoverBase.value.set(base[0], base[1], base[2]);
  }

  const scene = new THREE.Scene();

  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.autoClear = false;

  /*
   * The view from space: a star sphere and the sun, drawn behind the
   * globe by a perspective camera so that directions are correct, and
   * sunlight on the cells with a dark ambient. The stars sit in an
   * inertial frame that setSpace() turns about the pole by the sidereal
   * angle; the drag rotation applies on top of both. The camera keeps
   * station in a frame between the planet's (space.frame 0: it hangs over
   * one point, the sun and stars wheel past) and the stars' (1: it holds
   * still against the stars while the planet turns beneath it once per
   * sidereal day and the sun creeps along the ecliptic): the globe turns
   * about the pole by that share of every advance of the sidereal angle,
   * summed in space.turn, so moving the slider changes the rate and never
   * the orientation.
   */
  const SKY_RADIUS = 100;
  const skyScene = new THREE.Scene();
  const skyCamera = new THREE.PerspectiveCamera(60, 1, 1, 10 * SKY_RADIUS);
  const stars = buildStars();
  skyScene.add(stars);
  const sun = buildSun();
  skyScene.add(sun);
  const space = { enabled: false, sun: new THREE.Vector3(1, 0, 0), sidereal: new THREE.Quaternion(), spin: null, frame: 0, turn: 0 };
  const POLE = new THREE.Vector3(0, 0, 1), frameTurn = new THREE.Quaternion(), frameQ = new THREE.Quaternion();
  function frameQuaternion() {
    return frameQ.copy(sphereQuaternion).multiply(frameTurn.setFromAxisAngle(POLE, space.enabled ? space.turn : 0));
  }
  const updateRotation = () => rotationMatrix.makeRotationFromQuaternion(frameQuaternion());

  function buildStars() {
    const random = mulberry32(7);
    const gaussian = () => Math.sqrt(-2 * Math.log(1 - random())) * Math.cos(2 * Math.PI * random());
    const positions = [], colors = [], sizes = [];
    const put = (x, y, z, brightness, warmth, size) => {
      positions.push(SKY_RADIUS * x, SKY_RADIUS * y, SKY_RADIUS * z);
      colors.push(brightness * (1 + 0.25 * warmth), brightness * (1 + 0.05 * warmth), brightness * (1 - 0.3 * warmth));
      sizes.push(size);
    };
    for (let n = 0; n < 5000; n++) {
      const z = 2 * random() - 1, phi = 2 * Math.PI * random(), r = Math.sqrt(1 - z * z);
      const brightness = 0.25 + 0.75 * random() ** 3;
      put(r * Math.cos(phi), r * Math.sin(phi), z, brightness, 2 * random() - 1, 1 + 3 * brightness ** 2);
    }
    const tilt = 62 * Math.PI / 180, cosT = Math.cos(tilt), sinT = Math.sin(tilt);
    for (let n = 0; n < 24000; n++) {
      const along = 2 * Math.PI * random(), off = 0.14 * gaussian();
      const x = Math.cos(off) * Math.cos(along), y = Math.cos(off) * Math.sin(along), z = Math.sin(off);
      put(x, y * cosT - z * sinT, y * sinT + z * cosT, 0.12 + 0.18 * random(), 0.5 * random(), 1 + random());
    }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    geometry.setAttribute('size', new THREE.Float32BufferAttribute(sizes, 1));
    const starMaterial = new THREE.ShaderMaterial({
      vertexColors: true, transparent: true, depthWrite: false, blending: THREE.AdditiveBlending,
      vertexShader: `
attribute float size;
varying vec3 vColor;
void main() {
  vColor = color;
  gl_PointSize = size;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}`,
      fragmentShader: `
varying vec3 vColor;
void main() {
  float r = 2.0 * length(gl_PointCoord - 0.5);
  float a = smoothstep(1.0, 0.2, r);
  gl_FragColor = vec4(vColor * a, a);
}`,
    });
    const points = new THREE.Points(geometry, starMaterial);
    points.frustumCulled = false;
    return points;
  }

  function buildSun() {
    const size = 256;
    const canvas = document.createElement('canvas');
    canvas.width = canvas.height = size;
    const context = canvas.getContext('2d');
    const gradient = context.createRadialGradient(size / 2, size / 2, 0, size / 2, size / 2, size / 2);
    gradient.addColorStop(0, 'rgba(255, 255, 255, 1)');
    gradient.addColorStop(0.16, 'rgba(255, 252, 240, 1)');
    gradient.addColorStop(0.24, 'rgba(255, 236, 190, 0.55)');
    gradient.addColorStop(0.5, 'rgba(255, 220, 160, 0.12)');
    gradient.addColorStop(1, 'rgba(255, 200, 120, 0)');
    context.fillStyle = gradient;
    context.fillRect(0, 0, size, size);
    const texture = new THREE.CanvasTexture(canvas);
    const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: texture, blending: THREE.AdditiveBlending, depthTest: false, depthWrite: false, transparent: true }));
    sprite.scale.setScalar(0.14 * SKY_RADIUS);
    sprite.renderOrder = 1;
    return sprite;
  }

  /*
   * The camera's flare from the sun, drawn over the whole picture in
   * screen space with half the frame height as the unit: a halo, a
   * horizontal streak and a starburst on the sun, and ghosts along the
   * line from the sun through the centre of the view. Each part is
   * [texture, width, height, tint, place, air], the place along that line
   * with the sun at 0 and the centre at 1. All are dimmed as the visible
   * part of the sun's disc is; 'air' parts take the colour of the sunlight
   * along each fragment's own line of sight, the others that disc's colour.
   */
  const SUN_ANGULAR_RADIUS = 0.012;
  const flareScene = new THREE.Scene();
  const flareCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, -1, 1);
  const flare = buildFlare();
  flareScene.add(flare.group);
  const flareSun = new THREE.Vector3();

  function buildFlare() {
    const texture = (width, height, value) => {
      const canvas = document.createElement('canvas');
      canvas.width = width; canvas.height = height;
      const context = canvas.getContext('2d');
      const image = context.createImageData(width, height);
      for (let j = 0; j < height; j++) {
        for (let i = 0; i < width; i++) {
          const linear = Math.min(1, Math.max(0, value(2 * (i + 0.5) / width - 1, 2 * (j + 0.5) / height - 1)));
          const v = 255 * (linear <= 0.0031308 ? 12.92 * linear : 1.055 * linear ** (1 / 2.4) - 0.055);
          const k = 4 * (j * width + i);
          image.data[k] = image.data[k + 1] = image.data[k + 2] = v;
          image.data[k + 3] = 255;
        }
      }
      context.putImageData(image, 0, 0);
      const map = new THREE.CanvasTexture(canvas);
      map.colorSpace = THREE.SRGBColorSpace;
      return map;
    };
    const edge = (r, inner) => r >= 1 ? 0 : r <= inner ? 1 : 1 - THREE.MathUtils.smoothstep(r, inner, 1);
    const halo = texture(256, 256, (x, y) => { const r = Math.hypot(x, y); return (0.55 * Math.exp(-r * r / 0.01) + 0.45 * (1 - r) ** 4) * edge(r, 0.9); });
    const streak = texture(512, 32, (x, y) => Math.exp(-((y / 0.3) ** 2)) * (1 - Math.abs(x)) ** 2 * (0.35 + 0.65 * Math.exp(-6 * Math.abs(x))));
    const random = mulberry32(11);
    const rays = Array.from({ length: 18 }, (_, n) => [Math.PI * n / 9 + 0.08 * (random() - 0.5), 0.5 + 0.5 * random()]);
    const burst = texture(256, 256, (x, y) => {
      const r = Math.hypot(x, y), theta = Math.atan2(y, x);
      let spikes = 0;
      for (const [angle, length] of rays) {
        const off = Math.abs(Math.atan2(Math.sin(theta - angle), Math.cos(theta - angle)));
        spikes += Math.exp(-((off * Math.max(r, 0.02) / 0.012) ** 2)) * Math.max(0, 1 - r / length) ** 2;
      }
      return spikes / (1 + 6 * r) * edge(r, 0.8);
    });
    const ring = texture(128, 128, (x, y) => { const r = Math.hypot(x, y); return (0.3 + 0.7 * THREE.MathUtils.smoothstep(r, 0.7, 0.96)) * edge(r, 0.93); });
    const blob = texture(128, 128, (x, y) => { const r = Math.hypot(x, y); return Math.exp(-r * r / 0.18) * edge(r, 0.6); });
    const parts = [
      [halo, 0.9, 0.9, [0.24, 0.225, 0.2], 0, true],
      [streak, 2.4, 0.05, [0.06, 0.085, 0.13], 0],
      [burst, 0.6, 0.6, [0.09, 0.087, 0.078], 0, true],
      [ring, 0.1, 0.1, [0.045, 0.032, 0.018], 1.35],
      [blob, 0.06, 0.06, [0.03, 0.07, 0.04], 1.7],
      [ring, 0.24, 0.24, [0.012, 0.017, 0.03], 2.0],
      [ring, 0.13, 0.13, [0.025, 0.015, 0.032], 2.45],
      [blob, 0.36, 0.36, [0.01, 0.02, 0.028], 2.9],
    ].map(([map, width, height, tint, place, air = false]) => {
      const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map, blending: THREE.AdditiveBlending, depthTest: false, depthWrite: false, transparent: true }));
      return { sprite, width, height, tint: new THREE.Color(...tint), place, air };
    });
    const group = new THREE.Group();
    for (const part of parts) group.add(part.sprite);
    group.visible = false;
    return { group, parts, textures: [halo, streak, burst, ring, blob] };
  }

  /*
   * The colour the air gives the sun's light seen past the globe: the
   * transmission of a line of sight whose closest height to the globe is h
   * (globe radii, rescaled to the air's own scale height as the limb shell
   * draws it), through the Rayleigh air and the aerosol on both sides of
   * the tangent point and the ozone shell. The flare is dimmed by its
   * brightness, and its streak and ghosts take its colour, at the middle of
   * the part of the sun's disc that clears the limb.
   */
  const limbTint = new THREE.Color(1, 1, 1), toSun = new THREE.Vector3(), nearest = new THREE.Vector3(), sunTop = { value: 0 };
  function sunThroughLimb(h, out) {
    const b = 1 + Math.max(h, 0), chord = (o) => Math.sqrt(Math.max(o * o - b * b, 0));
    const air = 2 * Math.exp((1 - b) / SCALE_HEIGHT) * Math.sqrt(Math.PI * b / (2 * SCALE_HEIGHT));
    const aerosol = 2 * Math.exp((1 - b) / AEROSOL_HEIGHT) * Math.sqrt(Math.PI * b / (2 * AEROSOL_HEIGHT));
    const ozone = 2 * (chord(OZONE_HIGH) - chord(OZONE_LOW)) / (OZONE_HIGH - OZONE_LOW);
    const t = (k) => Math.exp(-(RAYLEIGH[k] * air + AEROSOL[k] * aerosol + OZONE[k] * ozone));
    return out.setRGB(t(0), t(1), t(2));
  }

  // The closest height above the globe of the camera's line of sight to the top of the sun's disc.
  function sunDiscTop() {
    toSun.copy(sun.position).normalize();
    const eye = lighting.uCameraPosition.value, along = Math.max(-eye.dot(toSun), 0);
    return nearest.copy(eye).addScaledVector(toSun, along).length() - 1 + SUN_ANGULAR_RADIUS * along;
  }

  // The main camera sees the sun along the sky camera's line of sight, so its position decides the limb test.
  function updateFlare(aspect) {
    flare.group.visible = false;
    const strength = lighting.uSun.value * glow.uFade.value;
    if (!state.perspective || strength <= 0 || sun.position.z >= 0) return;
    flareSun.copy(sun.position).project(skyCamera);
    const inFrame = 1 - THREE.MathUtils.smoothstep(Math.max(Math.abs(flareSun.x), Math.abs(flareSun.y)), 0.75, 1);
    const distance = camera.position.length();
    const fromCentre = Math.acos(THREE.MathUtils.clamp(-camera.position.dot(sun.position) / (distance * SKY_RADIUS), -1, 1));
    const limb = Math.asin(Math.min(1 / distance, 1));
    const unhidden = THREE.MathUtils.smoothstep(fromCentre, limb - SUN_ANGULAR_RADIUS, limb + SUN_ANGULAR_RADIUS);
    const brightness = strength * inFrame * unhidden;
    if (brightness <= 0) return;
    flare.group.visible = true;
    toSun.copy(sun.position).normalize();
    const along = Math.max(-camera.position.dot(toSun), 0);
    const past = nearest.copy(camera.position).addScaledVector(toSun, along).length() - 1, disc = SUN_ANGULAR_RADIUS * along;
    sunThroughLimb(0.5 * (past + disc + Math.max(past - disc, 0)) * SCALE_HEIGHT / glow.uScaleHeight.value, limbTint);
    const dimmed = 0.2126 * limbTint.r + 0.7152 * limbTint.g + 0.0722 * limbTint.b;
    flareCamera.left = -aspect; flareCamera.right = aspect;
    flareCamera.updateProjectionMatrix();
    const size = Math.sqrt(lighting.uSun.value);
    const x = flareSun.x * aspect, y = flareSun.y, centreX = 2 * inset.shiftX / container.clientHeight, centreY = 2 * inset.shift / container.clientHeight;
    for (let k = 0; k < flare.parts.length; k++) {
      const { sprite, width, height, tint, place, air } = flare.parts[k];
      sprite.position.set(x + place * (centreX - x), y + place * (centreY - y), 0);
      sprite.scale.set(width * size, height * size, 1);
      if (air) sprite.material.color.copy(tint).multiplyScalar(brightness * dimmed);
      else sprite.material.color.copy(tint).multiplyScalar(brightness).multiply(limbTint);
    }
  }

  function mulberry32(seed) {
    let a = seed >>> 0;
    return () => {
      a = (a + 0x6D2B79F5) >>> 0;
      let t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  renderer.domElement.style.width = "100%";
  renderer.domElement.style.height = "100%";
  renderer.domElement.style.display = "block";
  container.appendChild(renderer.domElement);

  /*
   * Two cameras on the view axis. The orthographic one shows 500 / zoom
   * globe radii of height; the perspective one sits above the front of
   * the globe (or the map) at the height that shows the same span there,
   * so the zoom keeps its meaning and the camera flies in as it grows.
   */
  const FIELD_OF_VIEW = 40;
  const TAN_HALF = Math.tan(THREE.MathUtils.degToRad(FIELD_OF_VIEW / 2));
  const orthographic = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.1, 1000);
  orthographic.position.set(0, 0, 10);
  orthographic.lookAt(0, 0, 0);
  const perspective = new THREE.PerspectiveCamera(FIELD_OF_VIEW, 1, 0.02, 200);
  let camera = orthographic;

  // --- 4. Material & Rotation ---
  
  // FIX: Rotate the sphere -90 degrees around X initially.
  // This maps the Source Z (North) to World Y (Up).
  // This maps the Source X (Prime Meridian) to World X.
  // The Seam (Date Line) ends up at World -X or -Z depending on View.
  // We actually want X to point towards camera (Z) so Prime Meridian is center.
  // Rotate -90 X: (0,0,1)->(0,1,0) [N->Up]. (1,0,0)->(1,0,0).
  // Then Rotate -90 Y: (1,0,0)->(0,0,1) [Prime->Front].
  // Combined Euler: (-PI/2, -PI/2, 0).
  const rotationMatrix = new THREE.Matrix4(); 
  const sphereQuaternion = new THREE.Quaternion(); 
  
  // Initialize rotation: North Up, Prime Meridian Front
  const initialEuler = new THREE.Euler(-Math.PI/2, -Math.PI/2, 0, 'YXZ');
  sphereQuaternion.setFromEuler(initialEuler);
  rotationMatrix.makeRotationFromQuaternion(sphereQuaternion);

  const viewState = {
    blend: 0.0,
    targetBlend: 0.0,
    version: 0,
  };
  let disposed = false;
  const aborter = new AbortController();
  const on = (target, type, handler, options = {}) => target.addEventListener(type, handler, { ...options, signal: aborter.signal });

  const projectedMaterials = [];
  /*
   * Projects a material's vertices through the globe rotation and the
   * Equal Earth morph. `head` declares extra attributes and uniforms,
   * and `source` may reassign sourcePos from them before projection;
   * `uniforms` are shared objects whose values the caller can change.
   * Three keys its program cache by the onBeforeCompile source text,
   * which is the same for every material here, so the key must carry
   * the variant.
   */
  function projectMaterial(material, mapLift = 0.0, { head = '', source = '', uniforms = {} } = {}) {
    material.customProgramCacheKey = () => `projected:${head}${source}`;
    material.onBeforeCompile = (shader) => {
      shader.uniforms.uModelRotation = { value: rotationMatrix };
      shader.uniforms.uCenterTexture = { value: centerTexture };
      shader.uniforms.uTexSize = { value: new THREE.Vector2(width, height) };
      shader.uniforms.uBlend = { value: viewState.blend };
      shader.uniforms.uMapLift = { value: mapLift };
      Object.assign(shader.uniforms, uniforms);
      material.userData.shader = shader;
      shader.vertexShader = vertexHead + head + shader.vertexShader;
      shader.vertexShader = shader.vertexShader.replace('#include <begin_vertex>', vertexLogic.replace('// [vertex source]', source));
    };
    projectedMaterials.push(material);
  }

  const lighting = { uSunDirection: { value: new THREE.Vector3(1, 0, 0) }, uCameraPosition: { value: new THREE.Vector3(0, 0, 1e5) }, uLighting: { value: 0 }, uAmbient: { value: 0.004 }, uSun: { value: 1 } };

  /*
   * The light pass: one texel per cell, on the cell textures' layout,
   * drawn by a full-screen triangle into an RGBA8 target whenever a frame
   * lands or the sun turns a tenth of a degree, at most every LIGHT_INTERVAL ms.
   * r and g are the square roots of the share of the sunbeam reaching the
   * cell's ground and its cloud top past what lies sunward (occlusion());
   * b carries the cloud top's slope toward the sun as sign·√|slope| about
   * texel 128, and a the cosine of the top's tilt, the slope exaggerated
   * CLOUD_RELIEF (twenty) times for the shading as the terrain's is.
   */
  const LIGHT_INTERVAL = 66, LIGHT_TURN = Math.cos(THREE.MathUtils.degToRad(0.1)), MAX_STEPS = 48;
  const lightTarget = new THREE.WebGLRenderTarget(width, height, { type: THREE.UnsignedByteType, format: THREE.RGBAFormat, minFilter: THREE.NearestFilter, magFilter: THREE.NearestFilter, depthBuffer: false, generateMipmaps: false });
  const lightUniforms = {
    uCenterTexture: { value: centerTexture }, uNeighboursA: { value: null }, uNeighboursB: { value: null }, uCellSurface: { value: cellSurface }, uCellTerrain: { value: cellTerrain },
    uSunDirection: lighting.uSunDirection, uTexWidth: { value: width }, uCellCount: { value: cellCounter }, uSpacing: { value: 0 }, uHighest: { value: 0 },
  };
  const lightCamera = new THREE.Camera();
  let lightScene = null;

  function buildLightPass() {
    const index = new Map();
    let n = 0;
    for (const cell of grid) index.set(cell, n++);
    const first = new Float32Array(4 * width * height), second = new Float32Array(2 * width * height);
    const unit = (c) => { const x = centerData[3 * c], y = centerData[3 * c + 1], z = centerData[3 * c + 2], r = Math.hypot(x, y, z); return [x / r, y / r, z / r]; };
    let arc = 0, pairs = 0;
    n = 0;
    for (const cell of grid) {
      const around = (cell.neighbors || []).slice(0, 6).map((other) => index.get(other));
      const [cx, cy, cz] = unit(n);
      for (let k = 0; k < 6; k++) {
        const j = around[k] ?? n;
        if (k < 4) first[4 * n + k] = j; else second[2 * n + k - 4] = j;
        if (k < around.length) { const [x, y, z] = unit(j); arc += Math.acos(Math.min(1, cx * x + cy * y + cz * z)); pairs++; }
      }
      n++;
    }
    lightUniforms.uNeighboursA.value = cellTexture(first, THREE.RGBAFormat, THREE.FloatType);
    lightUniforms.uNeighboursB.value = cellTexture(second, THREE.RGFormat, THREE.FloatType);
    lightUniforms.uSpacing.value = arc / Math.max(1, pairs);
    const lightMaterial = new THREE.RawShaderMaterial({
      glslVersion: THREE.GLSL3, uniforms: lightUniforms, depthTest: false, depthWrite: false,
      vertexShader: `in vec3 position;
void main() { gl_Position = vec4(position.xy, 0.0, 1.0); }`,
      fragmentShader: `precision highp float; precision highp int; precision highp sampler2D;
uniform sampler2D uCenterTexture;
uniform sampler2D uNeighboursA;
uniform sampler2D uNeighboursB;
uniform sampler2D uCellSurface;
uniform sampler2D uCellTerrain;
uniform vec3 uSunDirection;
uniform int uTexWidth;
uniform int uCellCount;
uniform float uSpacing;
uniform float uHighest;
out vec4 cellLight;
const float METRE = ${glsl(METRE)};
const int MAX_STEPS = ${MAX_STEPS};
ivec2 texel(int i) { int row = int((float(i) + 0.5) / float(uTexWidth)); return ivec2(i - row * uTexWidth, row); }
vec3 centreOf(int i) { return normalize(texelFetch(uCenterTexture, texel(i), 0).xyz); }
void neighboursOf(int i, out int js[6]) {
  vec4 a = texelFetch(uNeighboursA, texel(i), 0);
  vec2 b = texelFetch(uNeighboursB, texel(i), 0).xy;
  js[0] = int(a.x); js[1] = int(a.y); js[2] = int(a.z); js[3] = int(a.w); js[4] = int(b.x); js[5] = int(b.y);
}
// The cell nearest p, walking from cur to whichever neighbour lies nearer until none does.
int walk(int cur, vec3 p) {
  float near = dot(centreOf(cur), p);
  for (int it = 0; it < 4; it++) {
    int js[6];
    neighboursOf(cur, js);
    int best = cur;
    for (int m = 0; m < 6; m++) { float d = dot(centreOf(js[m]), p); if (d > near) { near = d; best = js[m]; } }
    if (best == cur) break;
    cur = best;
  }
  return cur;
}
/*
 * The share of the sunbeam reaching a point h0 above cell i's centre c
 * (heights in globe radii, as drawn): marching sunward
 * along the great circle in steps of half the cell spacing, or longer when
 * MAX_STEPS would not take the ray above the highest top, where the ray at
 * arc a stands at (1 + h0) cos e / cos(a + e) - 1 for the sun's elevation
 * e. Ground above sea level blocks the part of a step's height range below
 * it (the sea never does: the globe's shadow is the beam's own); each
 * cell's cloud is a slab from its base to its top holding the optical
 * depth 7.7 A / (1 - A) of its albedo A evenly, which the step's path crosses
 * over the part of its height range inside the slab. The cell's own cloud
 * is left out, being drawn over its ground already, and so is the globe's
 * shadow, which the light through the air (sunDepth) holds.
 */
float occlusion(int i, vec3 c, float h0) {
  float mu = clamp(dot(c, uSunDirection), -1.0, 1.0);
  float e = asin(mu);
  if (e < -acos(1.0 / (1.0 + h0))) return 1.0;
  vec3 t = uSunDirection - mu * c;
  float tl = length(t);
  if (tl < 1.0e-6) return 1.0;
  t /= tl;
  float k0 = (1.0 + h0) * sqrt(max(0.0, 1.0 - mu * mu));
  float reach = acos(min(1.0, k0 / (1.0 + uHighest))) - e;
  if (reach <= 0.0) return 1.0;
  float ds = max(0.5 * uSpacing, reach / float(MAX_STEPS));
  int steps = min(MAX_STEPS, int(ceil(reach / ds)));
  float depth = 0.0, clear = 1.0, hPrev = h0;
  int cur = i;
  for (int k = 0; k < MAX_STEPS; k++) {
    if (k >= steps) break;
    float a0 = float(k) * ds, a1 = a0 + ds, a = a0 + 0.5 * ds;
    cur = walk(cur, cos(a) * c + sin(a) * t);
    float hNext = k0 / cos(min(a1 + e, 1.5)) - 1.0;
    float lo = min(hPrev, hNext), hi = max(hPrev, hNext);
    if (-e > a0 && -e < a1) lo = k0 - 1.0;
    float range = max(hi - lo, 1.0e-7);
    hPrev = hNext;
    if (cur == i) continue;
    float g = METRE * texelFetch(uCellTerrain, texel(cur), 0).r;
    if (g > 0.0 && g > lo) clear = min(clear, clamp((hi - g) / range, 0.0, 1.0));
    vec4 s = texelFetch(uCellSurface, texel(cur), 0);
    if (s.y > 0.0) {
      float tau = min(8.0, 7.7 * s.y / max(1.0 - s.y, 1.0e-3));
      float top = METRE * s.z, base = METRE * s.w;
      float inside = max(0.0, min(hi, top) - max(lo, base));
      depth += tau / max(top - base, 1.0e-5) * length(vec2(ds, hi - lo)) * inside / range;
    }
    if (clear * exp(-depth) < 1.0e-3) return 0.0;
  }
  return clear * exp(-depth);
}
void main() {
  int i = int(gl_FragCoord.y) * uTexWidth + int(gl_FragCoord.x);
  if (i >= uCellCount) { cellLight = vec4(1.0, 1.0, 128.0 / 255.0, 1.0); return; }
  vec3 c = centreOf(i);
  vec4 s = texelFetch(uCellSurface, texel(i), 0);
  float ground = METRE * max(0.0, texelFetch(uCellTerrain, texel(i), 0).r);
  float toGround = occlusion(i, c, ground), toTop = 1.0, lean = 0.0, scale = 1.0;
  if (s.y > 0.01) {
    float top = METRE * s.z;
    toTop = occlusion(i, c, max(top, ground));
    int js[6];
    neighboursOf(i, js);
    vec3 gradient = vec3(0.0);
    float n = 0.0;
    for (int m = 0; m < 6; m++) {
      if (js[m] == i) continue;
      vec3 d = centreOf(js[m]) - c;
      vec4 o = texelFetch(uCellSurface, texel(js[m]), 0);
      gradient += o.y * (METRE * o.z - top) / dot(d, d) * d;
      n += 1.0;
    }
    gradient *= ${glsl(CLOUD_RELIEF)} * 2.0 / n;
    float tilt = clamp(-dot(gradient, uSunDirection), -1.0, 1.0);
    lean = sign(tilt) * sqrt(abs(tilt));
    scale = inversesqrt(1.0 + dot(gradient, gradient));
  }
  cellLight = vec4(sqrt(toGround), sqrt(toTop), (128.0 + 127.0 * lean) / 255.0, scale);
}`,
    });
    const triangle = new THREE.BufferGeometry();
    triangle.setAttribute('position', new THREE.BufferAttribute(new Float32Array([-1, -1, 0, 3, -1, 0, -1, 3, 0]), 3));
    const quad = new THREE.Mesh(triangle, lightMaterial);
    quad.frustumCulled = false;
    lightScene = new THREE.Scene();
    lightScene.add(quad);
  }

  function runLightPass() {
    const now = performance.now();
    if (!lightPass.dirty && lightPass.sun.dot(space.sun) >= LIGHT_TURN) return;
    if (now - lightPass.last < LIGHT_INTERVAL) return;
    if (!lightScene) buildLightPass();
    lightUniforms.uHighest.value = METRE * Math.max(lightPass.cloud, lightPass.land) + 1e-6;
    renderer.setRenderTarget(lightTarget);
    renderer.render(lightScene, lightCamera);
    renderer.setRenderTarget(null);
    lightPass.dirty = false;
    lightPass.last = now;
    lightPass.sun.copy(space.sun);
  }

  const material = new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide });
  projectMaterial(material, 0.0, {
    uniforms: { ...lighting, ...colorMap, uCellLight: { value: lightTarget.texture } },
    head: `
uniform vec3 uSunDirection;
uniform vec3 uCameraPosition;
uniform float uLighting;
uniform float uAmbient;
uniform float uSun;
uniform float uColorMode;
uniform sampler2D uCellColors;
uniform sampler2D uCellSurface;
uniform sampler2D uCellValues;
uniform vec3 uStops[${MAX_STOPS}];
uniform float uStopCount;
uniform vec2 uValueMap;
uniform vec3 uMissing;
uniform vec3 uFlat;
uniform vec3 uCoverBase;
uniform sampler2D uCellTerrain;
uniform float uTerrain;
uniform sampler2D uCellLight;
attribute vec3 slope;
attribute vec2 corner;
vec2 cellUv(float i) { return (vec2(mod(i, uTexSize.x), floor(i / uTexSize.x)) + 0.5) / uTexSize; }
vec3 terrainGrey(float z) {
  float g = z < 0.0 ? 0.03 + 0.12 * clamp((z + 6000.0) / 6000.0, 0.0, 1.0) : 0.22 + 0.6 * clamp(z / 5000.0, 0.0, 1.0);
  return vec3(g);
}
${atmosphere}
vec3 srgbToLinear(vec3 c) { return mix(c / 12.92, pow((c + 0.055) / 1.055, vec3(2.4)), step(vec3(0.04045), c)); }
vec3 paletteColor(float t) {
  float x = clamp(t, 0.0, 1.0) * (uStopCount - 1.0);
  int k = int(min(uStopCount - 2.0, floor(x)));
  return mix(uStops[k], uStops[k + 1], x - float(k));
}
`,
    source: `
  vec4 surface = texture2D(uCellSurface, uv);
  if (uColorMode > 0.5) {
    float value = texture2D(uCellValues, uv).r;
    if (uColorMode < 1.5) vColor.rgb = texture2D(uCellColors, uv).rgb;
    else if (uColorMode > 3.5) vColor.rgb = uFlat;
    else if (value > 1.0e29) vColor.rgb = uTerrain > 0.5 ? terrainGrey(texture2D(uCellTerrain, uv).r) : uMissing;
    else if (uColorMode < 2.5) vColor.rgb = srgbToLinear(paletteColor(value * uValueMap.x + uValueMap.y));
    else vColor.rgb = srgbToLinear(uCoverBase + (1.0 - exp(-max(0.0, value) * uValueMap.x)) * (1.0 - uCoverBase));
  }
  /*
   * The Satellite view's light on the cell: the sunbeam through the air
   * down to the ground or the cloud top, on the terrain's or the cloud's
   * facet, less what lies sunward casts off it, plus sky light; the air
   * between the surface and the camera; the sea's glint. The cloud top the
   * beam reaches is, at each corner, the mean of the tops of the cells
   * meeting there ('corner' holds the other two) weighted by their cloud's
   * opacity, so the colour and the reach of the sunset light on cloud vary
   * smoothly across the deck; what lies sunward casts off the cell's own
   * top. Where the sun stands high the day side's own look takes over
   * (daylight()), the cloud tops' relief kept; below a sun cosine of NIGHT_DEEP
   * even the air a grazing view crosses lies in the globe's shadow, leaving
   * the ambient alone. The sunlit part goes through a camera's toe,
   * L² / (L + CAMERA_TOE), which leaves the day side alone and crushes the
   * dim light of the penumbra as a photograph does; applied to each channel
   * it also deepens the colour of dark things, applied to the brightness
   * alone it keeps their hue, and TOE_SATURATION sets the share of the
   * former.
   */
  if (uLighting > 0.0) {
    vec3 n = normalize(position);
    float mu = dot(n, uSunDirection);
    vec4 cell = texture2D(uCellLight, uv);
    float lean = (cell.b * 255.0 - 128.0) / 127.0;
    float facet = dot(normalize(slope), uSunDirection);
    float groundUnblocked = cell.r * cell.r;
    float ground = 1.0 + max(texture2D(uCellTerrain, uv).r, 0.0) * METRE;
    float day = daylight(mu);
    float cloudFacet = (mu + lean * abs(lean)) * cell.a;
    float cloudUnblocked = cell.g * cell.g;
    mat3 spin = mat3(uModelRotation);
    vec3 nView = spin * n;
    vec3 toCamera = normalize(uCameraPosition - nView);
    vec3 halfway = normalize(spin * uSunDirection + toCamera);
    float glint = pow(max(0.0, dot(nView, halfway)), 90.0) * surface.x * (1.0 - surface.y) * smoothstep(0.0, 0.025, mu) * groundUnblocked;
    vec3 lit = vColor.rgb * uAmbient;
    if (day < 1.0 && mu > ${glsl(NIGHT_DEEP)}) {
      vec4 next = texture2D(uCellSurface, cellUv(corner.x)), last = texture2D(uCellSurface, cellUv(corner.y));
      float cover = surface.y + next.y + last.y;
      float top = max(ground, 1.0 + METRE * (surface.y * surface.z + next.y * next.z + last.y * last.z) / max(cover, 1.0e-6));
      vec3 light = mix(illumination(ground, mu, facet, 0.0, groundUnblocked, 1.0, 1.0),
                       illumination(top, mu, cloudFacet, CLOUD_ROUGHNESS, cloudUnblocked, CLOUD_TINT, CLOUD_WARMTH_RISE), surface.y);
      float visible = mix(ground, top, surface.y);
      vec3 view = normalize(mix(toCamera * spin, n, uBlend));
      vec3 airDepth;
      vec3 haze = max(airAlong(visible * n, view, uSunDirection, AIR_NODES, false, airDepth) - light * bakedAir(visible), 0.0);
      lit = vColor.rgb * (uAmbient + uSun * light * viewLoss(visible, max(dot(n, view), 0.0)))
        + uSun * (0.9 * glint * sunColour(ground, mu) + haze);
    }
    if (day > 0.0) {
      vec3 beam = dayBeam(mu);
      float diffuse = mix(max(0.0, facet) * groundUnblocked, max(0.0, cloudFacet) * cloudUnblocked, surface.y);
      float slant = pow(1.0 - max(0.0, dot(nView, toCamera)), 2.0) * (1.0 - uBlend);
      vec3 dayLit = vColor.rgb * (uAmbient + uSun * (diffuse * beam + DAY_SKY)) + uSun * (0.9 * glint * beam + 0.45 * slant * DAY_AIR);
      lit = mix(lit, dayLit, day);
    }
    vec3 sunlit = max(lit - vColor.rgb * uAmbient, 0.0);
    float luma = dot(sunlit, BRIGHTNESS);
    vec3 toed = mix(sunlit * (luma / (luma + CAMERA_TOE)), sunlit * sunlit / (sunlit + CAMERA_TOE), TOE_SATURATION);
    lit = vColor.rgb * uAmbient + toed;
    vColor.rgb = mix(vColor.rgb, lit, uLighting);
  }
`,
  });

  const mesh = new THREE.Mesh(geometry, material);
  mesh.frustumCulled = false;
  scene.add(mesh);

  /*
   * The sunlit air beyond the limb, on the back faces of a shell that the
   * globe hides inside its silhouette. Each ray is taken at its closest
   * height h, measured in a drawn scale height H that stays at least
   * GLOW_PIXELS of a pixel so the rim is resolved when far out, and
   * rescaled to the air's own scale height there: sunlight scattered once
   * along the ray's near and far halves, the far half seen through the
   * near, with the globe's shadow on both (limbLight), so that with the
   * sun behind the limb the band is red low down, yellow-white above and
   * blue higher up. That band shows as the sun nears the horizon at the
   * closest point and stands near the ray's own direction; elsewhere the day
   * side's blue rim takes over, its column exp(-h / H) for a drawn scale
   * height of DAY_RIM_HEIGHT under the same pixel floor. When the floor
   * stretches the drawn air well past the real scale height, a pixel holds
   * the whole band and shows its light summed over the band's real heights,
   * which is the thin red-orange ring of a sunset seen from far away.
   */
  const GLOW_PIXELS = 0.6, GLOW_LARGEST = 1 / 60, GLOW_CUT = 8, DAY_RIM_HEIGHT = 0.005, DAY_RIM_CUT = 6, GLOW_SHELL = 1 + GLOW_CUT * GLOW_LARGEST;
  const glow = {
    uModelRotation: { value: rotationMatrix },
    uSunDirection: lighting.uSunDirection,
    uCameraPosition: lighting.uCameraPosition,
    uSun: lighting.uSun,
    uScaleHeight: { value: SCALE_HEIGHT },
    uRimHeight: { value: DAY_RIM_HEIGHT },
    uFade: { value: 0 },
  };
  const glowMaterial = new THREE.ShaderMaterial({
    uniforms: glow, side: THREE.BackSide, transparent: true, depthWrite: false, blending: THREE.AdditiveBlending,
    vertexShader: `
varying vec3 vWorld;
void main() {
  vWorld = (modelMatrix * vec4(position, 1.0)).xyz;
  gl_Position = projectionMatrix * viewMatrix * vec4(vWorld, 1.0);
}`,
    fragmentShader: `
uniform mat4 uModelRotation;
uniform vec3 uSunDirection;
uniform vec3 uCameraPosition;
uniform float uSun;
uniform float uScaleHeight;
uniform float uRimHeight;
uniform float uFade;
varying vec3 vWorld;
${atmosphere}
void main() {
  vec3 ray = normalize(vWorld - uCameraPosition);
  vec3 closest = vWorld - dot(vWorld, ray) * ray;
  if (dot(closest - uCameraPosition, ray) < 0.0) closest = uCameraPosition;
  float b = length(closest);
  float h = max(b - 1.0, 0.0) / uScaleHeight, hRim = max(b - 1.0, 0.0) / uRimHeight;
  if (h > ${glsl(GLOW_CUT)} && hRim > ${glsl(DAY_RIM_CUT)}) discard;
  vec3 sun = mat3(uModelRotation) * uSunDirection;
  float mu = dot(closest / b, sun);
  float twilight = (1.0 - daylight(mu)) * smoothstep(-0.2, 0.6, dot(ray, sun));
  vec3 light = vec3(0.0);
  if (twilight > 0.0 && h < ${glsl(GLOW_CUT)}) {
    vec3 up = closest / b;
    light = limbLight(up * (1.0 + h * SCALE_HEIGHT), ray, sun);
    float far = smoothstep(1.0, 3.0, uScaleHeight / SCALE_HEIGHT);
    if (far > 0.0) {
      vec3 whole = (0.3 * limbLight(up * (1.0 + 0.15 * SCALE_HEIGHT), ray, sun) + 0.5 * limbLight(up * (1.0 + 0.55 * SCALE_HEIGHT), ray, sun)
        + 1.2 * limbLight(up * (1.0 + 1.4 * SCALE_HEIGHT), ray, sun) + 3.0 * limbLight(up * (1.0 + 3.5 * SCALE_HEIGHT), ray, sun)) / ${glsl(GLOW_CUT)};
      light = mix(light, whole * ${glsl(GLOW_CUT / 1.5)} * exp(-h / 1.5) / (1.0 - exp(${glsl(-GLOW_CUT / 1.5)})), far);
    }
    light *= 1.0 - smoothstep(${glsl(GLOW_CUT - 2)}, ${glsl(GLOW_CUT)}, h);
  }
  if (twilight < 1.0) light = mix(1.2 * max(exp(-hRim) - exp(-${glsl(DAY_RIM_CUT)}), 0.0) / (1.0 - exp(-${glsl(DAY_RIM_CUT)})) * dayRim(mu) * DAY_RIM, light, twilight);
  gl_FragColor = vec4(uFade * uSun * light, 1.0);
  #include <colorspace_fragment>
}`,
  });
  const glowShell = new THREE.Mesh(new THREE.SphereGeometry(GLOW_SHELL, 64, 32), glowMaterial);
  glowShell.frustumCulled = false;
  glowShell.renderOrder = 1;
  glowShell.visible = false;
  scene.add(glowShell);

  /*
   * The sun's disc and halo, and the flare's halo and rays on the sun, each
   * fragment dimmed and reddened by the air its line of sight from the
   * camera crosses past the globe (sunThroughAir), as sunThroughLimb() does
   * for one height, so that where the sun stands in the band at the limb its
   * glare is reddened with it. Glare beyond the top of the disc is that
   * top's light, so the height is taken no higher than the line of sight to
   * the disc's top (uSunTop, from sunDiscTop()). The sky camera sits at the
   * origin, so a fragment's direction is the main camera's ray there; the
   * flare, drawn in screen space, finds it through the sky camera's inverse
   * projection.
   */
  const airOnSun = { uCameraPosition: lighting.uCameraPosition, uScaleHeight: glow.uScaleHeight, uSunTop: sunTop, uSkyInverse: { value: skyCamera.projectionMatrixInverse } };
  const sunThroughAir = `uniform vec3 uCameraPosition;
uniform float uScaleHeight;
uniform float uSunTop;
${atmosphere}
vec3 sunThroughAir(vec3 toSun) {
  float past = min(length(uCameraPosition + max(-dot(uCameraPosition, toSun), 0.0) * toSun) - 1.0, uSunTop);
  return exp(-2.0 * sunDepth(1.0 + max(past, 0.0) * SCALE_HEIGHT / uScaleHeight, 0.0));
}`;
  sun.material.customProgramCacheKey = () => 'sun through the air';
  sun.material.onBeforeCompile = (shader) => {
    Object.assign(shader.uniforms, airOnSun);
    shader.vertexShader = `varying vec3 vSkyDirection;
${shader.vertexShader.replace('#include <fog_vertex>', `#include <fog_vertex>
  vSkyDirection = (vec4(mvPosition.xyz, 0.0) * viewMatrix).xyz;`)}`;
    shader.fragmentShader = `varying vec3 vSkyDirection;
${sunThroughAir}
${shader.fragmentShader.replace('#include <map_fragment>', `#include <map_fragment>
  diffuseColor.rgb *= sunThroughAir(normalize(vSkyDirection));`)}`;
  };
  for (const { sprite } of flare.parts.filter((part) => part.air)) {
    sprite.material.customProgramCacheKey = () => 'flare through the air';
    sprite.material.onBeforeCompile = (shader) => {
      Object.assign(shader.uniforms, airOnSun);
      shader.vertexShader = `varying vec2 vScreen;
${shader.vertexShader.replace('#include <fog_vertex>', `#include <fog_vertex>
  vScreen = gl_Position.xy / gl_Position.w;`)}`;
      shader.fragmentShader = `uniform mat4 uSkyInverse;
varying vec2 vScreen;
${sunThroughAir}
${shader.fragmentShader.replace('#include <map_fragment>', `#include <map_fragment>
  vec4 far = uSkyInverse * vec4(vScreen, 1.0, 1.0);
  diffuseColor.rgb *= sunThroughAir(normalize(far.xyz / far.w));`)}`;
    };
  }
  

  // Create GUI / Buttons
  const toggleMode = function() {
    if (viewState.targetBlend === 1.0) {
      viewState.targetBlend = 0.0;
    } else {
      viewState.targetBlend = 1.0;
    }
  }

  const ui = document.createElement('div');
  ui.style.position = 'absolute';
  ui.style.top = '20px';
  ui.style.right = '20px';
  ui.style.zIndex = '999';
  ui.style.display = 'flex';
  
  const btnToggleMode = document.createElement('button');
  btnToggleMode.innerText = "Toggle";
  btnToggleMode.style.padding = "8px 16px";
  btnToggleMode.style.cursor = "pointer";
  btnToggleMode.onclick = () => toggleMode('sphere');

  ui.appendChild(btnToggleMode);
  if (controls) container.appendChild(ui);

  // --- 5. Animation Loop ---
  const state = { isDragging: false, lastX: 0, lastY: 0, zoom: 150, pan: new THREE.Vector3(0, 0, 0), lastVector: null, perspective: false };
  const viewHeight = () => 500 / state.zoom;
  const cameraDistance = () => 1 - viewState.blend + viewHeight() / (2 * TAN_HALF);

  /*
   * setInsets() centres the view in the part of the canvas that page
   * elements leave uncovered: the picture moves up by `inset.shift` and
   * right by `inset.shiftX` pixels, easing towards the targets over the
   * 0.3 s in which the page's panel slides, so the two move together,
   * through a translation of the projection itself, so the zoom and the pan keep
   * their meaning and the raycaster (which reads the projection) stays exact.
   */
  const inset = { target: 0, shift: 0, targetX: 0, shiftX: 0, from: 0, fromX: 0, since: 0, matrix: new THREE.Matrix4() };
  const INSET_SECONDS = 0.3;
  function applyInset(target) {
    if (!inset.shift && !inset.shiftX) return;
    target.projectionMatrix.premultiply(inset.matrix.makeTranslation(2 * inset.shiftX / container.clientWidth, 2 * inset.shift / container.clientHeight, 0));
    target.projectionMatrixInverse.copy(target.projectionMatrix).invert();
  }

  let busy = 0;
  function render() {
    if (disposed) return;
    requestAnimationFrame(render);
    const started = performance.now();
    try { draw(); } finally { busy += performance.now() - started; }
  }
  function draw() {
    if (container.clientHeight === 0) return;

    if (viewState.targetBlend === 1.0) {
      btnToggleMode.innerText = "To Spherical";
    } else {
      btnToggleMode.innerText = "To Equal Earth";
    }

    if (Math.abs(viewState.blend - viewState.targetBlend) > 0.001) {
      viewState.blend += (viewState.targetBlend - viewState.blend) * 0.05;
      viewState.version++;
    } else {
       viewState.blend = viewState.targetBlend; 
    }
    for (const projected of projectedMaterials) {
      if (projected.userData.shader) projected.userData.shader.uniforms.uBlend.value = viewState.blend;
    }
    if (inset.shift !== inset.target || inset.shiftX !== inset.targetX) {
      const part = Math.min(1, (performance.now() - inset.since) / (1000 * INSET_SECONDS)), eased = 1 - (1 - part) ** 3;
      inset.shift = part < 1 ? inset.from + (inset.target - inset.from) * eased : inset.target;
      inset.shiftX = part < 1 ? inset.fromX + (inset.targetX - inset.fromX) * eased : inset.targetX;
      viewState.version++;
    }

    const aspect = container.clientWidth / container.clientHeight;
    if (state.perspective) {
      camera = perspective;
      camera.aspect = aspect;
      camera.position.set(state.pan.x, state.pan.y, cameraDistance());
      lighting.uCameraPosition.value.copy(camera.position);
    } else {
      camera = orthographic;
      const frustumSize = viewHeight();
      camera.left = -frustumSize * aspect / 2 + state.pan.x;
      camera.right = frustumSize * aspect / 2 + state.pan.x;
      camera.top = frustumSize / 2 + state.pan.y;
      camera.bottom = -frustumSize / 2 + state.pan.y;
      lighting.uCameraPosition.value.set(0, 0, 1e5);
    }
    camera.updateProjectionMatrix();
    applyInset(camera);

    glow.uFade.value = space.enabled ? 1 - THREE.MathUtils.smoothstep(viewState.blend, 0, 0.15) : 0;
    glowShell.visible = glow.uFade.value > 0 && lighting.uSun.value > 0;
    if (glowShell.visible) {
      const limb = state.perspective ? 2 * TAN_HALF * Math.sqrt(Math.max(camera.position.lengthSq() - 1, 0)) : viewHeight();
      const pixelFloor = GLOW_PIXELS * limb / container.clientHeight;
      glow.uScaleHeight.value = Math.min(Math.max(SCALE_HEIGHT, pixelFloor), GLOW_LARGEST);
      glow.uRimHeight.value = Math.min(Math.max(DAY_RIM_HEIGHT, pixelFloor), GLOW_LARGEST);
    }

    updateRotation();
    if (space.enabled && lighting.uLighting.value > 0) runLightPass();
    renderer.setClearColor(space.enabled ? 0x000000 : backgroundColor);
    renderer.clear();
    if (space.enabled) {
      skyCamera.aspect = aspect;
      skyCamera.fov = state.perspective ? FIELD_OF_VIEW : 60;
      skyCamera.updateProjectionMatrix();
      applyInset(skyCamera);
      stars.quaternion.copy(frameQuaternion()).multiply(space.sidereal);
      sun.position.copy(space.sun).applyQuaternion(frameQ).multiplyScalar(SKY_RADIUS);
      sunTop.value = sunDiscTop();
      renderer.render(skyScene, skyCamera);
      updateFlare(aspect);
    } else flare.group.visible = false;
    renderer.render(scene, camera);
    if (flare.group.visible) renderer.render(flareScene, flareCamera);
  }
  requestAnimationFrame(render);

  // --- 6. Interaction ---
  const raycaster = new THREE.Raycaster();
  const planeZ0 = new THREE.Plane(new THREE.Vector3(0, 0, 1), 0);
  const sphereOrigin = new THREE.Sphere(new THREE.Vector3(0,0,0), 1.0);
  const intersectPoint = new THREE.Vector3();

  function getCursorOnWorld(clientX, clientY) {
    const rect = container.getBoundingClientRect();
    const x = ((clientX - rect.left) / rect.width) * 2 - 1;
    const y = -((clientY - rect.top) / rect.height) * 2 + 1;
    raycaster.setFromCamera({ x, y }, camera);

    if (viewState.targetBlend < 0.5) {
      if (raycaster.ray.intersectSphere(sphereOrigin, intersectPoint)) {
        return intersectPoint.clone().normalize(); 
      }
    } else {
      if (raycaster.ray.intersectPlane(planeZ0, intersectPoint)) {
        return inverseEqualEarthToVector(intersectPoint.x, intersectPoint.y);
      }
    }
    return null;
  }

  const canvas = renderer.domElement;
  on(canvas, 'contextmenu', (e) => e.preventDefault());

  on(canvas, 'wheel', (e) => {
    e.preventDefault();
    const zoomSpeed = 0.001;
    state.zoom += -e.deltaY * zoomSpeed * state.zoom; 
    state.zoom = Math.max(10, Math.min(state.zoom, 10000)); 
    viewState.version++;
  }, { passive: false });

  function panBy(dx, dy) {
    if (container.clientHeight <= 0) return;
    const pxToWorld = viewHeight() / container.clientHeight;
    state.pan.x -= dx * pxToWorld;
    state.pan.y += dy * pxToWorld;
    viewState.version++;
  }

  on(window, 'keydown', (e) => {
    if (e.altKey || e.ctrlKey || e.metaKey || ['INPUT', 'SELECT', 'TEXTAREA'].includes(e.target.tagName)) return;
    const step = e.shiftKey ? 120 : 30;
    const move = { ArrowLeft: [step, 0], ArrowRight: [-step, 0], ArrowUp: [0, step], ArrowDown: [0, -step] }[e.key];
    if (!move) return;
    e.preventDefault();
    panBy(move[0], move[1]);
  });

  const VIEW_AXIS = new THREE.Vector3(0, 0, 1);
  function rollBy(angle) {
    sphereQuaternion.premultiply(new THREE.Quaternion().setFromAxisAngle(VIEW_AXIS, angle)).normalize();
    updateRotation();
  }

  /*
   * Turns the globe so that the point grabbed at state.lastVector
   * follows the pointer.
   */
  function dragGlobeTo(clientX, clientY) {
    const currentVector = getCursorOnWorld(clientX, clientY);
    if (currentVector && state.lastVector) {
      const axis = new THREE.Vector3().crossVectors(state.lastVector, currentVector);
      const angle = Math.acos(Math.max(-1, Math.min(1, state.lastVector.dot(currentVector))));
      if (angle > 0.0001) {
        sphereQuaternion.premultiply(new THREE.Quaternion().setFromAxisAngle(axis.normalize(), angle)).normalize();
        updateRotation();
      }
    }
    state.lastVector = currentVector;
  }

  /*
   * Zooms by `factor` while the point under (clientX, clientY) stays put:
   * a screen offset from the view's centre spans viewHeight() /
   * clientHeight world units per pixel in the orthographic view and at
   * the front of the globe in perspective.
   */
  function zoomAbout(factor, clientX, clientY) {
    if (container.clientHeight <= 0) return;
    const before = viewHeight();
    state.zoom = Math.max(10, Math.min(state.zoom * factor, 10000));
    const shrink = (before - viewHeight()) / container.clientHeight, rect = container.getBoundingClientRect();
    state.pan.x += (clientX - rect.left - rect.width / 2 - inset.shiftX) * shrink;
    state.pan.y -= (clientY - rect.top - rect.height / 2 + inset.shift) * shrink;
  }

  /*
   * A mouse or pen drags the globe with the main button, rolls it about
   * the line of sight with alt or meta held, and pans with the right
   * button. On a touch screen one finger drags the globe; two pan with
   * their midpoint and zoom with their spread about it, and roll with
   * their turn once it passes TWIST_START, so that a pinch doesn't wobble.
   */
  const TWIST_START = 0.2;
  const touches = new Map();
  let pair = null;
  canvas.style.touchAction = 'none';
  canvas.style.webkitUserSelect = canvas.style.userSelect = 'none';
  canvas.style.webkitTapHighlightColor = 'transparent';

  function touchesChanged() {
    const [first, second] = touches.values();
    state.lastVector = first && !second ? getCursorOnWorld(first.x, first.y) : null;
    pair = second ? { a: { ...first }, b: { ...second }, twist: 0, twisting: false } : null;
  }

  function pairMoved() {
    const [a, b] = touches.values();
    const motion = twoFingerMotion(pair.a, pair.b, a, b);
    panBy(motion.dx, motion.dy);
    zoomAbout(motion.scale, motion.x, motion.y);
    if (pair.twisting) rollBy(-motion.turn);
    else pair.twisting = Math.abs(pair.twist += motion.turn) > TWIST_START;
    pair.a = { ...a };
    pair.b = { ...b };
  }

  on(canvas, 'pointerdown', (e) => {
    canvas.setPointerCapture(e.pointerId);
    if (e.pointerType === 'touch') {
      touches.set(e.pointerId, { x: e.clientX, y: e.clientY });
      touchesChanged();
      return;
    }
    state.isDragging = true;
    state.lastX = e.clientX;
    state.lastY = e.clientY;
    state.lastVector = getCursorOnWorld(e.clientX, e.clientY);
    if (e.buttons === 1) canvas.style.cursor = 'grabbing';
    else if (e.buttons === 2) canvas.style.cursor = 'move';
  });

  on(canvas, 'pointermove', (e) => {
    if (e.pointerType === 'touch') {
      const touch = touches.get(e.pointerId);
      if (!touch) return;
      touch.x = e.clientX;
      touch.y = e.clientY;
      if (pair) pairMoved();
      else dragGlobeTo(e.clientX, e.clientY);
      viewState.version++;
      return;
    }
    if (!state.isDragging) return;
    const dx = e.clientX - state.lastX;
    const dy = e.clientY - state.lastY;
    viewState.version++;
    if (e.buttons === 2) panBy(dx, dy);
    else if (e.buttons === 1 && (e.altKey || e.metaKey)) rollBy((dx + dy) * 0.01);
    else if (e.buttons === 1) dragGlobeTo(e.clientX, e.clientY);
    state.lastX = e.clientX;
    state.lastY = e.clientY;
  });

  function release(e) {
    if (e.pointerType === 'touch') {
      if (touches.delete(e.pointerId)) touchesChanged();
      return;
    }
    state.isDragging = false;
    canvas.style.cursor = 'default';
  }
  on(canvas, 'pointerup', release);
  on(canvas, 'pointercancel', release);
  
  const resizeObserver = new ResizeObserver(() => {
    if (container.clientWidth > 0) renderer.setSize(container.clientWidth, container.clientHeight, false);
  });
  resizeObserver.observe(container);

  /*
   * Screen position of a point given in the grid's coordinates, following
   * the same rotation and sphere-to-map morph as the shader. out[2] is
   * positive where the point is on the visible side of the globe.
   */
  const mapPoint = [0, 0];
  function projectPoint(x, y, z, out) {
    const e = rotationMatrix.elements;
    const rx = e[0] * x + e[4] * y + e[8] * z;
    const ry = e[1] * x + e[5] * y + e[9] * z;
    const rz = e[2] * x + e[6] * y + e[10] * z;
    const blend = viewState.blend;
    let fx = rx, fy = ry;
    if (blend > 0) {
      equalEarth(Math.asin(Math.max(-1, Math.min(1, ry))), Math.atan2(rx, rz), mapPoint);
      fx += (mapPoint[0] - rx) * blend;
      fy += (mapPoint[1] - ry) * blend;
    }
    if (state.perspective) {
      const distance = cameraDistance();
      const half = (distance - rz * (1 - blend)) * TAN_HALF;
      const aspect = container.clientWidth / container.clientHeight;
      out[0] = ((fx - state.pan.x) / (half * aspect) + 1) / 2 * container.clientWidth + inset.shiftX;
      out[1] = (1 - (fy - state.pan.y) / half) / 2 * container.clientHeight - inset.shift;
      out[2] = blend > 0.5 ? 1 : rx * state.pan.x + ry * state.pan.y + rz * distance - 1;
    } else {
      out[0] = (fx - camera.left) / (camera.right - camera.left) * container.clientWidth + inset.shiftX;
      out[1] = (camera.top - fy) / (camera.top - camera.bottom) * container.clientHeight - inset.shift;
      out[2] = blend > 0.5 ? 1 : rz;
    }
    return out;
  }

  /*
   * Contour lines of a cell field on the globe, drawn on the GPU: the
   * Delaunay triangles (the three cells around each vertex) carry the
   * field linearly between cell centres, and the fragment shader draws a
   * line wherever it crosses a multiple of the interval, `width` pixels
   * wide whatever the zoom. Where lines come closer than a couple of
   * pixels it fades to their mean coverage. Each corner is projected
   * against its triangle's first cell so nothing tears at the map seam.
   * Triangles with a missing corner, or flat on a level, draw nothing.
   * update() takes the field and the contour interval.
   */
  function addContourLayer({ color = 0xffffff, opacity = 0.7, width: lineWidth = 1 } = {}) {
    const values = cellTexture(new Float32Array(width * height), THREE.RedFormat, THREE.FloatType);
    const uniforms = { uContourValues: { value: values }, uContourStep: { value: 1 }, uContourWidth: { value: lineWidth } };
    const contourMaterial = new THREE.MeshBasicMaterial({ color, transparent: true, opacity, depthWrite: false, side: THREE.DoubleSide });
    projectMaterial(contourMaterial, 0.008, {
      head: `
uniform sampler2D uContourValues;
uniform float uContourStep;
attribute float valueCell;
varying float vContour;
varying float vContourMissing;
`,
      source: `
  float contourValue = texture2D(uContourValues, (vec2(mod(valueCell, uTexSize.x), floor(valueCell / uTexSize.x)) + 0.5) / uTexSize).r;
  vContourMissing = contourValue > 1.0e29 ? 1.0 : 0.0;
  vContour = contourValue / uContourStep;
`,
      uniforms,
    });
    const project = contourMaterial.onBeforeCompile;
    contourMaterial.onBeforeCompile = (shader) => {
      project(shader);
      shader.fragmentShader = `
uniform float uContourWidth;
varying float vContour;
varying float vContourMissing;
` + shader.fragmentShader.replace('#include <color_fragment>', `#include <color_fragment>
  float spacing = length(vec2(dFdx(vContour), dFdy(vContour)));
  if (vContourMissing > 0.0 || spacing < 1.0e-5) discard;
  float line = clamp(0.5 * uContourWidth + 0.5 - abs(fract(vContour - 0.5) - 0.5) / spacing, 0.0, 1.0);
  diffuseColor.a *= mix(line, min(1.0, uContourWidth * spacing), smoothstep(0.5, 1.0, spacing));
  if (diffuseColor.a < 0.002) discard;`);
    };
    const contourGeometry = new THREE.BufferGeometry();
    const mesh = new THREE.Mesh(contourGeometry, contourMaterial);
    mesh.frustumCulled = false;
    mesh.visible = false;
    scene.add(mesh);
    let shown = false, ready = false;

    function build() {
      const corners = cellsOnVertex.filter((cells) => cells && cells.length === 3).flat();
      const positions = new Float32Array(3 * corners.length), reference = new Float32Array(corners.length), own = new Float32Array(corners.length);
      const lift = 1.003;
      corners.forEach((cell, k) => {
        const x = centerData[3 * cell], y = centerData[3 * cell + 1], z = centerData[3 * cell + 2], scale = lift / Math.hypot(x, y, z);
        positions[3 * k] = x * scale; positions[3 * k + 1] = y * scale; positions[3 * k + 2] = z * scale;
        reference[k] = corners[k - (k % 3)];
        own[k] = cell;
      });
      contourGeometry.setAttribute('position', new THREE.BufferAttribute(positions, 3).onUpload(disposeArray));
      contourGeometry.setAttribute('cellIndex', new THREE.BufferAttribute(reference, 1).onUpload(disposeArray));
      contourGeometry.setAttribute('valueCell', new THREE.BufferAttribute(own, 1).onUpload(disposeArray));
      contourGeometry.computeBoundingSphere();
      ready = true;
    }

    function update(field, step) {
      if (!ready) build();
      const array = values.image.data;
      for (let c = 0; c < cellCounter; c++) { const v = field[c]; array[c] = v === v ? v : MISSING; }
      values.needsUpdate = true;
      uniforms.uContourStep.value = step;
      mesh.visible = shown;
    }

    return {
      update,
      setColor(value) { contourMaterial.color.set(value); },
      setVisible(visible) { shown = visible; mesh.visible = shown && ready; },
      dispose() { scene.remove(mesh); contourGeometry.dispose(); contourMaterial.dispose(); values.dispose(); },
    };
  }

  /*
   * A graticule: parallels and meridians every `step` degrees as line
   * segments on the sphere, projected with the cell shader. Each segment
   * references a cell of similar longitude so both of its ends unwrap
   * against the same centre and nothing tears at the map seam; the
   * meridians run only between the outermost parallels.
   */
  function addGraticuleLayer({ color = 0xffffff, opacity = 0.4 } = {}) {
    const capacity = 2 * (36 * 180 + 72 * 180);
    const positions = new Float32Array(3 * capacity);
    const cellIndex = new Float32Array(capacity);
    const graticuleGeometry = new THREE.BufferGeometry();
    const positionAttribute = new THREE.BufferAttribute(positions, 3).setUsage(THREE.DynamicDrawUsage);
    const cellAttribute = new THREE.BufferAttribute(cellIndex, 1).setUsage(THREE.DynamicDrawUsage);
    graticuleGeometry.setAttribute('position', positionAttribute);
    graticuleGeometry.setAttribute('cellIndex', cellAttribute);
    graticuleGeometry.setDrawRange(0, 0);
    const graticuleMaterial = new THREE.LineBasicMaterial({ color, transparent: true, opacity });
    projectMaterial(graticuleMaterial, 0.006);
    const lines = new THREE.LineSegments(graticuleGeometry, graticuleMaterial);
    lines.frustumCulled = false;
    lines.visible = false;
    scene.add(lines);
    const lift = 1.002;
    const buckets = new Map();
    const bucketOf = (lat, lon) => `${Math.floor((lat + Math.PI / 2) / (Math.PI / 36))},${(Math.floor((lon + Math.PI) / (Math.PI / 36)) + 72) % 72}`;
    for (let i = 0; i < centerData.length / 3; i++) {
      const x = centerData[3 * i], y = centerData[3 * i + 1], z = centerData[3 * i + 2];
      const key = bucketOf(Math.atan2(z, Math.hypot(x, y)), Math.atan2(y, x));
      (buckets.get(key) ?? buckets.set(key, []).get(key)).push(i);
    }
    const referenceFor = (lat, lon) => {
      const px = Math.cos(lat) * Math.cos(lon), py = Math.cos(lat) * Math.sin(lon), pz = Math.sin(lat);
      const latBand = Math.floor((lat + Math.PI / 2) / (Math.PI / 36)), lonBand = Math.floor((lon + Math.PI) / (Math.PI / 36));
      let best = 0, bestDot = -2;
      for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) {
        const candidates = buckets.get(`${latBand + a},${(lonBand + b + 72) % 72}`);
        if (!candidates) continue;
        for (const i of candidates) {
          const dot = px * centerData[3 * i] + py * centerData[3 * i + 1] + pz * centerData[3 * i + 2];
          if (dot > bestDot) { bestDot = dot; best = i; }
        }
      }
      return best;
    };
    let built = 0;

    function update(step) {
      if (step === built) return;
      built = step;
      let n = 0;
      const put = (lat, lon, reference) => {
        positions[3 * n] = lift * Math.cos(lat) * Math.cos(lon);
        positions[3 * n + 1] = lift * Math.cos(lat) * Math.sin(lon);
        positions[3 * n + 2] = lift * Math.sin(lat);
        cellIndex[n] = reference;
        n++;
      };
      const segment = (lat0, lon0, lat1, lon1) => {
        if (n + 2 > capacity) return;
        const reference = referenceFor(0.5 * (lat0 + lat1), 0.5 * (lon0 + lon1));
        put(lat0, lon0, reference);
        put(lat1, lon1, reference);
      };
      const rad = Math.PI / 180, piece = 2 * rad;
      const outermost = step * Math.floor(89 / step);
      for (let lat = -outermost; lat <= outermost + 1e-9; lat += step) {
        for (let lon = -180; lon < 180; lon += 2) segment(lat * rad, lon * rad, lat * rad, (lon + 2) * rad);
      }
      for (let lon = -180; lon < 180; lon += step) {
        for (let lat = -outermost; lat < outermost - 1e-9; lat += 2) segment(lat * rad, lon * rad, Math.min(lat + 2, outermost) * rad, lon * rad);
      }
      graticuleGeometry.setDrawRange(0, n);
      positionAttribute.needsUpdate = true;
      cellAttribute.needsUpdate = true;
    }

    return {
      update,
      setColor(value) { graticuleMaterial.color.set(value); },
      setVisible(visible) { lines.visible = visible; },
      dispose() { scene.remove(lines); graticuleGeometry.dispose(); graticuleMaterial.dispose(); },
    };
  }

  /*
   * A layer of line segments given as unit-sphere positions, two per
   * segment, each segment with a reference cell for the projection's
   * seam handling; set() replaces them all.
   */
  function addSegmentLayer({ color = 0xffffff, opacity = 0.6 } = {}) {
    const segmentGeometry = new THREE.BufferGeometry();
    segmentGeometry.setDrawRange(0, 0);
    const segmentMaterial = new THREE.LineBasicMaterial({ color, transparent: true, opacity });
    projectMaterial(segmentMaterial, 0.007);
    const lines = new THREE.LineSegments(segmentGeometry, segmentMaterial);
    lines.frustumCulled = false;
    lines.visible = false;
    scene.add(lines);
    const lift = 1.0025;

    function set(positions, cells) {
      const n = positions.length / 3;
      const lifted = Float32Array.from(positions, (v) => lift * v);
      const cellIndex = new Float32Array(n);
      for (let k = 0; k < n; k++) cellIndex[k] = cells[k >> 1];
      segmentGeometry.setAttribute('position', new THREE.BufferAttribute(lifted, 3));
      segmentGeometry.setAttribute('cellIndex', new THREE.BufferAttribute(cellIndex, 1));
      segmentGeometry.setDrawRange(0, n);
    }

    return {
      set,
      setColor(value) { segmentMaterial.color.set(value); },
      setVisible(visible) { lines.visible = visible; },
      dispose() { scene.remove(lines); segmentGeometry.dispose(); segmentMaterial.dispose(); },
    };
  }

  /*
   * A layer of arrows, one per cell, drawn in the cell's tangent plane
   * by the vertex shader from a wind texture: the geometry is static
   * (six vertices per cell with a role) and update() only refreshes the
   * texture. Each arrow is centred on its cell and spans up to half of
   * the cell's diameter at referenceSpeed.
   */
  function addArrowLayer({ color = 0xffffff, opacity = 0.8 } = {}) {
    const positions = new Float32Array(3 * 6 * cellCounter);
    const cellIndex = new Float32Array(6 * cellCounter);
    const role = new Float32Array(6 * cellCounter);
    const radius = new Float32Array(6 * cellCounter);
    for (let c = 0; c < cellCounter; c++) {
      for (let k = 0; k < 6; k++) {
        positions.set(centerData.slice(3 * c, 3 * c + 3), 18 * c + 3 * k);
        cellIndex[6 * c + k] = c;
        role[6 * c + k] = k;
        radius[6 * c + k] = cellRadius[c];
      }
    }
    const arrowGeometry = new THREE.BufferGeometry();
    arrowGeometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    arrowGeometry.setAttribute('cellIndex', new THREE.BufferAttribute(cellIndex, 1));
    arrowGeometry.setAttribute('role', new THREE.BufferAttribute(role, 1));
    arrowGeometry.setAttribute('cellRadius', new THREE.BufferAttribute(radius, 1));
    const windBuffer = new Float32Array(width * height * 4);
    const windTexture = new THREE.DataTexture(windBuffer, width, height, THREE.RGBAFormat, THREE.FloatType);
    windTexture.minFilter = THREE.NearestFilter;
    windTexture.magFilter = THREE.NearestFilter;
    const uniforms = { uWindTexture: { value: windTexture }, uReferenceSpeed: { value: 20 } };
    const arrowMaterial = new THREE.LineBasicMaterial({ color, transparent: opacity < 1, opacity });
    projectMaterial(arrowMaterial, 0.01, {
      uniforms,
      head: `
attribute float role;
attribute float cellRadius;
uniform sampler2D uWindTexture;
uniform float uReferenceSpeed;
`,
      source: `
  vec3 wind = texture2D(uWindTexture, uv).xyz;
  vec3 tangentWind = wind - dot(wind, cellCenter) * cellCenter;
  float speed = length(tangentWind);
  vec3 dir = speed > 1e-6 ? tangentWind / speed : vec3(0.0);
  vec3 side = cross(cellCenter, dir);
  float len = cellRadius * min(1.0, speed / uReferenceSpeed);
  float head = 0.35 * len;
  vec3 tip = 1.004 * cellCenter + 0.5 * len * dir;
  sourcePos = tip;
  if (role < 0.5) sourcePos = 1.004 * cellCenter - 0.5 * len * dir;
  else if (role > 2.5 && role < 3.5) sourcePos = tip - 0.866 * head * dir + 0.5 * head * side;
  else if (role > 4.5) sourcePos = tip - 0.866 * head * dir - 0.5 * head * side;
`,
    });
    const lines = new THREE.LineSegments(arrowGeometry, arrowMaterial);
    lines.frustumCulled = false;
    lines.visible = false;
    scene.add(lines);

    function update(vectors, { referenceSpeed = 20 } = {}) {
      for (let c = 0; c < cellCounter; c++) {
        windBuffer[4 * c] = vectors[3 * c];
        windBuffer[4 * c + 1] = vectors[3 * c + 1];
        windBuffer[4 * c + 2] = vectors[3 * c + 2];
      }
      windTexture.needsUpdate = true;
      uniforms.uReferenceSpeed.value = referenceSpeed;
    }

    return {
      update,
      setVisible(visible) { lines.visible = visible; },
      dispose() { scene.remove(lines); arrowGeometry.dispose(); arrowMaterial.dispose(); windTexture.dispose(); },
    };
  }

  function canonicalQuaternion(lat, lon) {
    const tilt = new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0, 1, 0), THREE.MathUtils.degToRad(lat));
    const spin = new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0, 0, 1), -THREE.MathUtils.degToRad(lon));
    return new THREE.Quaternion().setFromEuler(initialEuler).multiply(tilt).multiply(spin).normalize();
  }

  return {
    updateColors, updateSurface, updateSlopes, updateValues, updateTerrain, setColorMap,
    setSpace({ enabled, sun: direction = null, sidereal = null, frame = null, ambient = 0.004, intensity = 1, perspective = enabled } = {}) {
      space.enabled = enabled;
      if (frame !== null) space.frame = Math.min(Math.max(frame, 0), 1);
      if (state.perspective !== perspective) { state.perspective = perspective; viewState.version++; }
      lighting.uLighting.value = enabled ? 1 : 0;
      lighting.uAmbient.value = ambient;
      lighting.uSun.value = intensity;
      if (direction) { space.sun.set(direction[0], direction[1], direction[2]); lighting.uSunDirection.value.copy(space.sun); }
      if (sidereal !== null) {
        if (space.spin !== null) space.turn = (space.turn + space.frame * (sidereal - space.spin)) % (2 * Math.PI);
        space.spin = sidereal;
        space.sidereal.setFromAxisAngle(POLE, -sidereal);
      }
    },
    addArrowLayer,
    addContourLayer,
    addGraticuleLayer,
    addSegmentLayer,
    setProjection(mode) { viewState.targetBlend = mode === 'map' ? 1.0 : 0.0; },
    projection: () => (viewState.targetBlend === 1.0 ? 'map' : 'sphere'),
    /*
     * The view is the point facing the camera and the roll about the
     * line of sight: the orientation with that point in front and north
     * up is the canonical one, and whatever rotation about the view axis
     * separates the actual orientation from it is the roll.
     */
    view() {
      const whole = frameQuaternion().clone();
      const front = new THREE.Vector3(0, 0, 1).applyQuaternion(whole.clone().invert());
      const lat = THREE.MathUtils.radToDeg(Math.atan2(front.z, Math.hypot(front.x, front.y))), lon = THREE.MathUtils.radToDeg(Math.atan2(front.y, front.x));
      const residual = whole.multiply(canonicalQuaternion(lat, lon).invert());
      let roll = THREE.MathUtils.radToDeg(2 * Math.atan2(residual.z, residual.w));
      if (roll > 180) roll -= 360; else if (roll <= -180) roll += 360;
      return { lat, lon, zoom: state.zoom, roll, x: state.pan.x, y: state.pan.y };
    },
    setView({ lat = null, lon = null, zoom = null, roll = null, x = null, y = null } = {}) {
      if (lat !== null || lon !== null || roll !== null) {
        const current = this.view();
        sphereQuaternion.copy(canonicalQuaternion(lat ?? current.lat, lon ?? current.lon));
        const turn = roll ?? current.roll;
        if (turn) sphereQuaternion.premultiply(new THREE.Quaternion().setFromAxisAngle(POLE, THREE.MathUtils.degToRad(turn))).normalize();
        if (space.enabled && space.turn) sphereQuaternion.multiply(new THREE.Quaternion().setFromAxisAngle(POLE, -space.turn)).normalize();
        updateRotation();
      }
      if (zoom !== null) state.zoom = Math.max(10, Math.min(zoom, 10000));
      if (x !== null) state.pan.x = x;
      if (y !== null) state.pan.y = y;
      viewState.version++;
    },
    projectPoint,
    unprojectPoint(px, py, out) {
      raycaster.setFromCamera({ x: (px / container.clientWidth) * 2 - 1, y: -(py / container.clientHeight) * 2 + 1 }, camera);
      let hit = null;
      if (viewState.targetBlend < 0.5) { if (raycaster.ray.intersectSphere(sphereOrigin, intersectPoint)) hit = intersectPoint.normalize(); }
      else if (raycaster.ray.intersectPlane(planeZ0, intersectPoint)) hit = inverseEqualEarthToVector(intersectPoint.x, intersectPoint.y);
      if (!hit) return null;
      const e = rotationMatrix.elements;
      out[0] = e[0] * hit.x + e[1] * hit.y + e[2] * hit.z;
      out[1] = e[4] * hit.x + e[5] * hit.y + e[6] * hit.z;
      out[2] = e[8] * hit.x + e[9] * hit.y + e[10] * hit.z;
      return out;
    },
    pixelsPerUnit: () => container.clientHeight / viewHeight(),
    setInsets({ top = 0, bottom = 0, left = 0, right = 0 } = {}) {
      const target = (bottom - top) / 2, targetX = (left - right) / 2;
      if (target === inset.target && targetX === inset.targetX) return;
      inset.from = inset.shift; inset.fromX = inset.shiftX; inset.since = performance.now();
      inset.target = target; inset.targetX = targetX;
    },
    viewVersion: () => viewState.version,
    takeBusyTime() { const ms = busy; busy = 0; return ms; },
    dispose: () => {
      disposed = true;
      aborter.abort();
      resizeObserver.disconnect();
      renderer.dispose();
      geometry.dispose();
      material.dispose();
      glowShell.geometry.dispose(); glowMaterial.dispose();
      lightTarget.dispose();
      if (lightScene) { for (const quad of lightScene.children) { quad.geometry.dispose(); quad.material.dispose(); } lightUniforms.uNeighboursA.value.dispose(); lightUniforms.uNeighboursB.value.dispose(); }
      for (const part of flare.parts) part.sprite.material.dispose();
      for (const map of flare.textures) map.dispose();
      centerTexture.dispose();
      cellColors.dispose(); cellSurface.dispose(); cellValues.dispose();
      renderer.domElement.remove();
      ui.remove();
    }
  };
}
