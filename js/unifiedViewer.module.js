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

  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.autoClear = false;

  /*
   * The view from space: a star sphere and the sun, drawn behind the
   * globe by a perspective camera so that directions are correct, and
   * sunlight on the cells with a dark ambient. The stars sit in an
   * inertial frame that setSpace() turns about the pole by the sidereal
   * angle; the drag rotation applies on top of both.
   */
  const SKY_RADIUS = 100;
  const skyScene = new THREE.Scene();
  const skyCamera = new THREE.PerspectiveCamera(60, 1, 1, 10 * SKY_RADIUS);
  const stars = buildStars();
  skyScene.add(stars);
  const sun = buildSun();
  skyScene.add(sun);
  const space = { enabled: false, sun: new THREE.Vector3(1, 0, 0), sidereal: new THREE.Quaternion() };

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

  const lighting = { uSunDirection: { value: new THREE.Vector3(1, 0, 0) }, uLighting: { value: 0 }, uAmbient: { value: 0.015 } };
  const material = new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide });
  projectMaterial(material, 0.0, {
    uniforms: lighting,
    head: `
uniform vec3 uSunDirection;
uniform float uLighting;
uniform float uAmbient;
`,
    source: `
  float daylight = uAmbient + (1.0 - uAmbient) * max(0.0, dot(normalize(position), uSunDirection));
  vColor.rgb *= mix(1.0, daylight, uLighting);
`,
  });

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
  if (controls) container.appendChild(ui);

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

    renderer.setClearColor(space.enabled ? 0x000000 : backgroundColor);
    renderer.clear();
    if (space.enabled) {
      skyCamera.aspect = aspect;
      skyCamera.updateProjectionMatrix();
      stars.quaternion.copy(sphereQuaternion).multiply(space.sidereal);
      sun.position.copy(space.sun).applyQuaternion(sphereQuaternion).multiplyScalar(SKY_RADIUS);
      renderer.render(skyScene, skyCamera);
    }
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

  function panBy(dx, dy) {
    if (container.clientHeight <= 0) return;
    const pxToWorld = (camera.top - camera.bottom) / container.clientHeight;
    state.pan.x -= dx * pxToWorld;
    state.pan.y += dy * pxToWorld;
    viewState.version++;
  }

  window.addEventListener('keydown', (e) => {
    if (e.altKey || e.ctrlKey || e.metaKey || ['INPUT', 'SELECT', 'TEXTAREA'].includes(e.target.tagName)) return;
    const step = e.shiftKey ? 120 : 30;
    const move = { ArrowLeft: [-step, 0], ArrowRight: [step, 0], ArrowUp: [0, -step], ArrowDown: [0, step] }[e.key];
    if (!move) return;
    e.preventDefault();
    panBy(move[0], move[1]);
  });

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
      panBy(dx, dy);
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
   * Contour lines of a cell field on the globe: marching triangles over
   * the Delaunay triangles (the three cells around each vertex), drawn
   * through the same projection as the cells. update() takes the field
   * and the contour interval.
   */
  function addContourLayer({ color = 0xffffff, opacity = 0.7 } = {}) {
    const triangles = cellsOnVertex.filter((cells) => cells && cells.length === 3);
    const capacity = 6 * triangles.length;
    const positions = new Float32Array(3 * capacity);
    const cellIndex = new Float32Array(capacity);
    const contourGeometry = new THREE.BufferGeometry();
    const positionAttribute = new THREE.BufferAttribute(positions, 3).setUsage(THREE.DynamicDrawUsage);
    const cellAttribute = new THREE.BufferAttribute(cellIndex, 1).setUsage(THREE.DynamicDrawUsage);
    contourGeometry.setAttribute('position', positionAttribute);
    contourGeometry.setAttribute('cellIndex', cellAttribute);
    contourGeometry.setDrawRange(0, 0);
    const contourMaterial = new THREE.LineBasicMaterial({ color, transparent: true, opacity });
    projectMaterial(contourMaterial, 0.008);
    const lines = new THREE.LineSegments(contourGeometry, contourMaterial);
    lines.frustumCulled = false;
    lines.visible = false;
    scene.add(lines);
    const lift = 1.003;

    function update(field, step) {
      let fmin = Infinity, fmax = -Infinity;
      for (const value of field) { if (value < fmin) fmin = value; if (value > fmax) fmax = value; }
      let n = 0;
      const crossing = (a, b, level, reference) => {
        const t = (level - field[a]) / (field[b] - field[a]);
        const x = centerData[3 * a] + t * (centerData[3 * b] - centerData[3 * a]);
        const y = centerData[3 * a + 1] + t * (centerData[3 * b + 1] - centerData[3 * a + 1]);
        const z = centerData[3 * a + 2] + t * (centerData[3 * b + 2] - centerData[3 * a + 2]);
        const norm = Math.hypot(x, y, z) / lift;
        positions[3 * n] = x / norm; positions[3 * n + 1] = y / norm; positions[3 * n + 2] = z / norm;
        cellIndex[n] = reference;
        n++;
      };
      for (let level = Math.ceil(fmin / step) * step; level < fmax && n + 2 <= capacity; level += step) {
        for (const [a, b, c] of triangles) {
          if (n + 2 > capacity) break;
          const sa = field[a] - level, sb = field[b] - level, sc = field[c] - level;
          const before = n;
          if (sa * sb < 0) crossing(a, b, level, a);
          if (sb * sc < 0) crossing(b, c, level, a);
          if (sc * sa < 0) crossing(c, a, level, a);
          if (n - before !== 2) n = before;
        }
      }
      contourGeometry.setDrawRange(0, n);
      positionAttribute.clearUpdateRanges(); positionAttribute.addUpdateRange(0, 3 * n); positionAttribute.needsUpdate = true;
      cellAttribute.clearUpdateRanges(); cellAttribute.addUpdateRange(0, n); cellAttribute.needsUpdate = true;
    }

    return {
      update,
      setColor(value) { contourMaterial.color.set(value); },
      setVisible(visible) { lines.visible = visible; },
      dispose() { scene.remove(lines); contourGeometry.dispose(); contourMaterial.dispose(); },
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
   * A layer of arrows, one per cell, drawn in the cell's tangent plane
   * by the vertex shader from a wind texture: the geometry is static
   * (six vertices per cell with a role) and update() only refreshes the
   * texture. Each arrow is centred on its cell and spans up to 80% of
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
  float len = 1.6 * cellRadius * min(1.0, speed / uReferenceSpeed);
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

  return {
    updateColors: dynamicColors ? updateColors : null,
    setSpace({ enabled, sun: direction = null, sidereal = 0, ambient = 0.015 } = {}) {
      space.enabled = enabled;
      lighting.uLighting.value = enabled ? 1 : 0;
      lighting.uAmbient.value = ambient;
      if (direction) { space.sun.set(direction[0], direction[1], direction[2]); lighting.uSunDirection.value.copy(space.sun); }
      space.sidereal.setFromAxisAngle(new THREE.Vector3(0, 0, 1), -sidereal);
    },
    addArrowLayer,
    addContourLayer,
    addGraticuleLayer,
    setProjection(mode) { viewState.targetBlend = mode === 'map' ? 1.0 : 0.0; },
    projection: () => (viewState.targetBlend === 1.0 ? 'map' : 'sphere'),
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
