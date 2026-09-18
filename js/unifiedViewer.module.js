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
  vec3 posMap = vec3(proj, uMapLift);

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
    dynamicColors = false,
  } = config;

  // --- 1. Geometry Generation ---
  const pData = [], kData = [], iData = [], idxData = [], texDataArray = [], centerData = [];
  const cellVertexStart = [], cellVertexCount = [];
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
  const colorAttribute = new THREE.BufferAttribute(new Uint8Array(kData), 3, true);
  if (!dynamicColors) colorAttribute.onUpload(disposeArray);
  geometry.setAttribute('color', colorAttribute);
  geometry.computeBoundingSphere();

  function updateColors(rgb) {
    const array = colorAttribute.array;
    for (let c = 0; c < cellCounter; c++) {
      const r = rgb[3 * c], g = rgb[3 * c + 1], b = rgb[3 * c + 2];
      let v = 3 * cellVertexStart[c];
      for (let k = 0; k < cellVertexCount[c]; k++) {
        array[v++] = r;
        array[v++] = g;
        array[v++] = b;
      }
    }
    colorAttribute.needsUpdate = true;
  }

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
  const initialEuler = new THREE.Euler(-Math.PI/2, -Math.PI/2, 0, 'YXZ');
  sphereQuaternion.setFromEuler(initialEuler);
  rotationMatrix.makeRotationFromQuaternion(sphereQuaternion);

  const viewState = {
    blend: 0.0,
    targetBlend: 0.0,
    version: 0,
  };

  const projectedMaterials = [];
  function projectMaterial(material, mapLift = 0.0) {
    material.onBeforeCompile = (shader) => {
      shader.uniforms.uModelRotation = { value: rotationMatrix };
      shader.uniforms.uCenterTexture = { value: centerTexture };
      shader.uniforms.uTexSize = { value: new THREE.Vector2(width, height) };
      shader.uniforms.uBlend = { value: viewState.blend };
      shader.uniforms.uMapLift = { value: mapLift };
      material.userData.shader = shader;
      shader.vertexShader = vertexHead + shader.vertexShader;
      shader.vertexShader = shader.vertexShader.replace('#include <begin_vertex>', vertexLogic);
    };
    projectedMaterials.push(material);
  }

  const material = new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide });
  projectMaterial(material);

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
      viewState.version++;
    } else {
       viewState.blend = viewState.targetBlend; 
    }
    for (const projected of projectedMaterials) {
      if (projected.userData.shader) projected.userData.shader.uniforms.uBlend.value = viewState.blend;
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
    viewState.version++;
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
    viewState.version++;

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
    out[0] = (fx - camera.left) / (camera.right - camera.left) * container.clientWidth;
    out[1] = (camera.top - fy) / (camera.top - camera.bottom) * container.clientHeight;
    out[2] = blend > 0.5 ? 1 : rz;
    return out;
  }

  /*
   * A layer of arrows, one per cell, drawn in the cell's tangent plane
   * and projected with the same shader as the cells so they follow the
   * globe in both views. update() takes a vector per cell in the grid's
   * coordinates (3 components per cell) and scales the arrow with the
   * speed up to referenceSpeed, at which it spans about a cell.
   */
  function addArrowLayer({ color = 0xffffff, opacity = 0.8 } = {}) {
    const centers = Float32Array.from(centerData);
    const positions = new Float32Array(3 * 6 * cellCounter);
    const cellIndex = new Float32Array(6 * cellCounter);
    for (let c = 0; c < cellCounter; c++) cellIndex.fill(c, 6 * c, 6 * c + 6);
    const arrowGeometry = new THREE.BufferGeometry();
    const positionAttribute = new THREE.BufferAttribute(positions, 3).setUsage(THREE.DynamicDrawUsage);
    arrowGeometry.setAttribute('position', positionAttribute);
    arrowGeometry.setAttribute('cellIndex', new THREE.BufferAttribute(cellIndex, 1));
    const arrowMaterial = new THREE.LineBasicMaterial({ color, transparent: opacity < 1, opacity });
    projectMaterial(arrowMaterial, 0.01);
    const lines = new THREE.LineSegments(arrowGeometry, arrowMaterial);
    lines.frustumCulled = false;
    lines.visible = false;
    scene.add(lines);
    const cellSpan = 0.8 * Math.sqrt(4 * Math.PI / cellCounter);
    const lift = 1.004;

    function update(vectors, { referenceSpeed = 20, stride = 1 } = {}) {
      positions.fill(0);
      for (let c = 0; c < cellCounter; c += stride) {
        const cx = centers[3 * c], cy = centers[3 * c + 1], cz = centers[3 * c + 2];
        let dx = vectors[3 * c], dy = vectors[3 * c + 1], dz = vectors[3 * c + 2];
        const radial = dx * cx + dy * cy + dz * cz;
        dx -= radial * cx; dy -= radial * cy; dz -= radial * cz;
        const speed = Math.hypot(dx, dy, dz);
        if (!(speed > 1e-6)) continue;
        const length = cellSpan * Math.min(1, speed / referenceSpeed);
        dx /= speed; dy /= speed; dz /= speed;
        const tx = cy * dz - cz * dy, ty = cz * dx - cx * dz, tz = cx * dy - cy * dx;
        const half = 0.5 * length, head = 0.35 * length, along = 0.866 * head, across = 0.5 * head;
        const tipX = lift * cx + half * dx, tipY = lift * cy + half * dy, tipZ = lift * cz + half * dz;
        let at = 18 * c;
        positions[at++] = lift * cx - half * dx; positions[at++] = lift * cy - half * dy; positions[at++] = lift * cz - half * dz;
        positions[at++] = tipX; positions[at++] = tipY; positions[at++] = tipZ;
        positions[at++] = tipX; positions[at++] = tipY; positions[at++] = tipZ;
        positions[at++] = tipX - along * dx + across * tx; positions[at++] = tipY - along * dy + across * ty; positions[at++] = tipZ - along * dz + across * tz;
        positions[at++] = tipX; positions[at++] = tipY; positions[at++] = tipZ;
        positions[at++] = tipX - along * dx - across * tx; positions[at++] = tipY - along * dy - across * ty; positions[at++] = tipZ - along * dz - across * tz;
      }
      positionAttribute.needsUpdate = true;
    }

    return {
      update,
      setVisible(visible) { lines.visible = visible; },
      dispose() { scene.remove(lines); arrowGeometry.dispose(); arrowMaterial.dispose(); },
    };
  }

  return {
    updateColors: dynamicColors ? updateColors : null,
    addArrowLayer,
    projectPoint,
    pixelsPerUnit: () => container.clientHeight / (camera.top - camera.bottom),
    viewVersion: () => viewState.version,
    dispose: () => {
      resizeObserver.disconnect();
      renderer.dispose();
      geometry.dispose();
      material.dispose();
      centerTexture.dispose();
    }
  };
}
