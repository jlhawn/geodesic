import * as THREE from "./three.module.js";
import { Stats } from "./stats.module.js";

// ----------------------------------------------------------------------------
// SHADER DEFINITIONS
// ----------------------------------------------------------------------------

const vertexHead = `
uniform mat4 uModelRotation; 
uniform sampler2D uCenterTexture;
uniform vec2 uTexSize; 
uniform float uBlend; // 0.0 = Sphere, 1.0 = Map

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

  // 2. Apply Rotation
  vec4 rotatedPos = uModelRotation * vec4(position, 1.0);
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
  vec3 posMap = vec3(proj, 0.0);

  // 5. MORPH
  vec3 finalPos = mix(posSphere, posMap, uBlend);
  vec3 transformed = finalPos;
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
    getColor = (cell) => cell.color 
  } = config;

  // --- 1. Geometry Generation ---
  const pData = [], kData = [], iData = [], idxData = [], texDataArray = [];
  const colorHelper = new THREE.Color();
  let vertexCounter = 0, cellCounter = 0;

  for (const cell of grid) {
    const cv = cell.centerVertex; 
    const verts = cell.vertices || []; 
    const currentCellID = cellCounter++;

    texDataArray.push(cv.x, cv.y, cv.z, 1.0);

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

  // --- 3. Scene & Buffers ---
  function disposeArray() { this.array = null; }
  const geometry = new THREE.BufferGeometry();
  geometry.setIndex(new THREE.BufferAttribute(new Uint32Array(iData), 1).onUpload(disposeArray));
  geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(pData), 3).onUpload(disposeArray));
  geometry.setAttribute('cellIndex', new THREE.BufferAttribute(new Float32Array(idxData), 1).onUpload(disposeArray));
  geometry.setAttribute('color', new THREE.BufferAttribute(new Uint8Array(kData), 3, true).onUpload(disposeArray));
  geometry.computeBoundingSphere();

  const scene = new THREE.Scene();
  scene.background = new THREE.Color(backgroundColor);

  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.domElement.style.width = "100%";
  renderer.domElement.style.height = "100%";
  renderer.domElement.style.display = "block";
  container.appendChild(renderer.domElement);

  const camera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.1, 1000);
  camera.position.set(0, 0, 10);
  camera.lookAt(0, 0, 0);

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
  const initialEuler = new THREE.Euler(-Math.PI/2, -Math.PI/2, 0, 'XYZ');
  sphereQuaternion.setFromEuler(initialEuler);
  rotationMatrix.makeRotationFromQuaternion(sphereQuaternion);

  const viewState = {
    blend: 0.0,
    targetBlend: 0.0
  };

  const material = new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide });
  material.onBeforeCompile = (shader) => {
    shader.uniforms.uModelRotation = { value: rotationMatrix };
    shader.uniforms.uCenterTexture = { value: centerTexture };
    shader.uniforms.uTexSize = { value: new THREE.Vector2(width, height) };
    shader.uniforms.uBlend = { value: 0.0 };
    material.userData.shader = shader;
    shader.vertexShader = vertexHead + shader.vertexShader;
    shader.vertexShader = shader.vertexShader.replace('#include <begin_vertex>', vertexLogic);
  };

  const mesh = new THREE.Mesh(geometry, material);
  mesh.frustumCulled = false; 
  scene.add(mesh);
  
  const stats = new Stats();
  document.body.appendChild(stats.dom);

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
  container.appendChild(ui);

  // --- 5. Animation Loop ---
  const state = { isDragging: false, lastX: 0, lastY: 0, zoom: 150, pan: new THREE.Vector3(0, 0, 0), lastVector: null };

  function render() {
    if (container.clientHeight === 0) return;

    if (viewState.targetBlend === 1.0) {
      btnToggleMode.innerText = "To Spherical";
    } else {
      btnToggleMode.innerText = "To Equal Earth";
    }

    if (Math.abs(viewState.blend - viewState.targetBlend) > 0.001) {
      viewState.blend += (viewState.targetBlend - viewState.blend) * 0.05;
      if (material.userData.shader) {
        material.userData.shader.uniforms.uBlend.value = viewState.blend;
      }
    } else {
       viewState.blend = viewState.targetBlend; 
    }

    const aspect = container.clientWidth / container.clientHeight;
    const frustumSize = 500 / state.zoom; 
    
    camera.left = -frustumSize * aspect / 2 + state.pan.x;
    camera.right = frustumSize * aspect / 2 + state.pan.x;
    camera.top = frustumSize / 2 + state.pan.y;
    camera.bottom = -frustumSize / 2 + state.pan.y;
    camera.updateProjectionMatrix();
    
    renderer.render(scene, camera);
    stats.update();
    requestAnimationFrame(render);
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
  canvas.addEventListener('contextmenu', e => e.preventDefault());

  canvas.addEventListener('wheel', (e) => {
    e.preventDefault();
    const zoomSpeed = 0.001;
    state.zoom += -e.deltaY * zoomSpeed * state.zoom; 
    state.zoom = Math.max(10, Math.min(state.zoom, 10000)); 
  }, { passive: false });

  canvas.addEventListener('mousedown', (e) => {
    state.isDragging = true;
    state.lastX = e.clientX;
    state.lastY = e.clientY;
    state.lastVector = getCursorOnWorld(e.clientX, e.clientY);
    if (e.buttons == 1) {
      canvas.style.cursor = 'grabbing';
    } else if (e.buttons == 2) {
      canvas.style.cursor = 'move';
    }
  });

  window.addEventListener('mousemove', (e) => {
    if (!state.isDragging) return;
    const dx = e.clientX - state.lastX;
    const dy = e.clientY - state.lastY;

    if (e.buttons === 2) {
       if (container.clientHeight > 0) {
        const worldHeight = (camera.top - camera.bottom);
        const pxToWorld = worldHeight / container.clientHeight;
        state.pan.x -= dx * pxToWorld;
        state.pan.y += dy * pxToWorld;
      }
    }
    else if (e.buttons === 1 && (e.altKey || e.metaKey)) {
       // Roll Axis: In View Space, Roll is Z.
       // We want to rotate around the view vector (Camera Z).
       const rollAxis = new THREE.Vector3(0, 0, 1); 
       const angle = (dx + dy) * 0.01;
       const qRot = new THREE.Quaternion().setFromAxisAngle(rollAxis, angle);
       sphereQuaternion.premultiply(qRot).normalize();
       rotationMatrix.makeRotationFromQuaternion(sphereQuaternion);
    }
    else if (e.buttons === 1) {
      const currentVector = getCursorOnWorld(e.clientX, e.clientY);
      if (currentVector && state.lastVector) {
        const axis = new THREE.Vector3().crossVectors(state.lastVector, currentVector);
        const dot = Math.max(-1, Math.min(1, state.lastVector.dot(currentVector)));
        const angle = Math.acos(dot);
        if (angle > 0.0001) {
          axis.normalize();
          const qRot = new THREE.Quaternion().setFromAxisAngle(axis, angle);
          sphereQuaternion.premultiply(qRot).normalize();
          rotationMatrix.makeRotationFromQuaternion(sphereQuaternion);
        }
      }
      state.lastVector = currentVector;
    }
    state.lastX = e.clientX;
    state.lastY = e.clientY;
  });

  window.addEventListener('mouseup', () => {
    state.isDragging = false;
    canvas.style.cursor = 'default';
  });
  
  const resizeObserver = new ResizeObserver(() => {
    if (container.clientWidth > 0) renderer.setSize(container.clientWidth, container.clientHeight, false);
  });
  resizeObserver.observe(container);

  return {
    dispose: () => {
      resizeObserver.disconnect();
      renderer.dispose();
      geometry.dispose();
      material.dispose();
      centerTexture.dispose();
    }
  };
}
