import * as THREE from "./three.module.js";
import { SetupCameraAndRenderer } from "./equalEarthCamera.module.js";
import { Stats } from "./stats.module.js";

// --- Shader Injection Blocks ---
const vertexHead = `
uniform mat4 uModelRotation; 

// -- Data Texture Lookups --
uniform sampler2D uCenterTexture;
uniform vec2 uTexSize; 
attribute float cellIndex; 

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
  // -- 1. Fetch Cell Center from Texture --
  float col = mod(cellIndex, uTexSize.x);
  float row = floor(cellIndex / uTexSize.x);
  vec2 uv = (vec2(col, row) + 0.5) / uTexSize; 
  vec3 cellCenter = texture2D(uCenterTexture, uv).xyz;

  // -- 2. Rotate --
  vec4 rotatedPos = uModelRotation * vec4(position, 1.0);
  vec4 rotatedCenter = uModelRotation * vec4(cellCenter, 1.0);
  
  // -- 3. Cartesian -> Spherical --
  float lon = atan(rotatedPos.y, rotatedPos.x);
  float centerLon = atan(rotatedCenter.y, rotatedCenter.x);

  // -- 4. Fix Tearing --
  float delta = lon - centerLon;
  if (delta > PI) lon -= 2.0 * PI;
  if (delta < -PI) lon += 2.0 * PI;

  float lat = asin(clamp(rotatedPos.z, -1.0, 1.0));

  // -- 5. Project --
  vec2 proj = projectEqualEarth(lat, lon);
  vec3 transformed = vec3(proj, 0.0);
`;

// --- Main Module ---

export function initEqualEarthViewer(container, grid, config = {}) {
  const {
    backgroundColor = 0x111111,
    getColor = (cell) => cell.color 
  } = config;

  // 1. PREPARE STORAGE
  const pData = []; // Positions
  const kData = []; // Colors
  const iData = []; // Indices
  const idxData = []; // Cell IDs
  const texDataArray = []; // Texture Data (Centers)
  
  const colorHelper = new THREE.Color();
  let vertexCounter = 0;
  let cellCounter = 0;

  for (const cell of grid) {
    const cv = cell.centerVertex; 
    const verts = cell.vertices || []; 
    const currentCellID = cellCounter++;

    // A. Texture Data
    texDataArray.push(cv.x, cv.y, cv.z, 1.0);

    // B. Colors
    const cVal = getColor(cell);
    if (cVal && typeof cVal === 'object' && 'r' in cVal) {
      colorHelper.setRGB(cVal.r, cVal.g, cVal.b);
    } else {
      colorHelper.set(cVal);
    }
    const r8 = Math.floor(colorHelper.r * 255);
    const g8 = Math.floor(colorHelper.g * 255);
    const b8 = Math.floor(colorHelper.b * 255);

    // C. Geometry (Perimeter Only)
    const baseIndex = vertexCounter;

    for (let i = 0; i < verts.length; i++) {
      const v = verts[i];
      pData.push(v.x, v.y, v.z);
      idxData.push(currentCellID); 
      kData.push(r8, g8, b8);
      vertexCounter++;
    }

    // D. Triangulation (Corner Fan: 0-1-2, 0-2-3...)
    const numTriangles = verts.length - 2;
    for (let k = 1; k <= numTriangles; k++) {
      iData.push(baseIndex, baseIndex + k, baseIndex + k + 1);
    }
  }

  // 2. CREATE DATA TEXTURE
  const width = Math.ceil(Math.sqrt(cellCounter));
  const height = Math.ceil(cellCounter / width);
  const floatBuffer = new Float32Array(width * height * 4);
  floatBuffer.set(texDataArray);
  
  const centerTexture = new THREE.DataTexture(floatBuffer, width, height, THREE.RGBAFormat, THREE.FloatType);
  centerTexture.minFilter = THREE.NearestFilter;
  centerTexture.magFilter = THREE.NearestFilter;
  centerTexture.needsUpdate = true;

  // 3. CREATE GEOMETRY
  function disposeArray() { this.array = null; }
  const geometry = new THREE.BufferGeometry();
  geometry.setIndex(new THREE.BufferAttribute(new Uint32Array(iData), 1).onUpload(disposeArray));
  geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(pData), 3).onUpload(disposeArray));
  geometry.setAttribute('cellIndex', new THREE.BufferAttribute(new Float32Array(idxData), 1).onUpload(disposeArray));
  geometry.setAttribute('color', new THREE.BufferAttribute(new Uint8Array(kData), 3, true).onUpload(disposeArray));
  geometry.computeBoundingSphere();

  const rotationMatrix = new THREE.Matrix4(); 
  const material = new THREE.MeshBasicMaterial({
    vertexColors: true,
    side: THREE.DoubleSide,
  });

  material.onBeforeCompile = (shader) => {
    shader.uniforms.uModelRotation = { value: rotationMatrix };
    shader.uniforms.uCenterTexture = { value: centerTexture };
    shader.uniforms.uTexSize = { value: new THREE.Vector2(width, height) };

    shader.vertexShader = vertexHead + shader.vertexShader;
    shader.vertexShader = shader.vertexShader.replace('#include <begin_vertex>', vertexLogic);
  };

  // 4. Camera, Renderer, and Scene.

  const { scene, render } = SetupCameraAndRenderer(container, rotationMatrix);
  scene.background = new THREE.Color(backgroundColor);

  const mesh = new THREE.Mesh(geometry, material);
  mesh.frustumCulled = false; 
  scene.add(mesh);

  const stats = new Stats();
  document.body.appendChild( stats.dom );

  function animate(t) {
    render();
    stats.update();

    requestAnimationFrame( animate );
  }
  requestAnimationFrame( animate );
}