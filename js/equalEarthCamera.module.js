import * as THREE from "./three.module.js";

const FIELD_OF_VIEW = 10;

const viewSize = new THREE.Vector2();

// --- Inverse Projection ---
// --- Constants ---
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
    return new THREE.Vector3(cosLat * Math.cos(lon), cosLat * Math.sin(lon), Math.sin(lat));
}

function SetupCameraAndRenderer(container, rotationMatrix) {
  viewSize.set(window.innerWidth, window.innerHeight);

  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.domElement.style.width = "100%";
  renderer.domElement.style.height = "100%";
  renderer.domElement.style.display = "block";
  container.appendChild(renderer.domElement);

  const camera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.1, 1000);
  camera.position.set(0, 0, 1);

  const scene = new THREE.Scene();

  const render = SetupCameraControls(camera, renderer, scene, rotationMatrix);

  return { scene, render };
};

function SetupCameraControls(camera, renderer, scene, rotationMatrix) {
  const domElement = renderer.domElement;
  const state = { isDragging: false, lastX: 0, lastY: 0, zoom: 150, pan: new THREE.Vector3(0, 0, 0), lastVector: null };
  
  if (domElement.clientHeight > 0) {
      renderer.setSize(domElement.clientWidth, domElement.clientHeight, false);
  }

  const raycaster = new THREE.Raycaster();
  const planeZ0 = new THREE.Plane(new THREE.Vector3(0, 0, 1), 0);
  const intersectPoint = new THREE.Vector3();

  function getCursorOnMap(clientX, clientY) {
      const rect = domElement.getBoundingClientRect();
      const x = ((clientX - rect.left) / rect.width) * 2 - 1;
      const y = -((clientY - rect.top) / rect.height) * 2 + 1;
      raycaster.setFromCamera({ x, y }, camera);
      const hit = raycaster.ray.intersectPlane(planeZ0, intersectPoint);
      return hit ? { x: hit.x, y: hit.y } : null;
  }

  domElement.addEventListener('contextmenu', e => e.preventDefault());

  domElement.addEventListener('wheel', (e) => {
    e.preventDefault();
    const zoomSpeed = 0.001;
    state.zoom += -e.deltaY * zoomSpeed * state.zoom; 
    state.zoom = Math.max(10, Math.min(state.zoom, 10000));
  }, { passive: false });

  domElement.addEventListener('mousedown', (e) => {
    state.isDragging = true;
    state.lastX = e.clientX;
    state.lastY = e.clientY;
    const mapPoint = getCursorOnMap(e.clientX, e.clientY);
    state.lastVector = mapPoint ? inverseEqualEarthToVector(mapPoint.x, mapPoint.y) : null;
    domElement.style.cursor = 'grabbing';
  });

  const sphereQuaternion = new THREE.Quaternion();

  window.addEventListener('mousemove', (e) => {
    if (!state.isDragging) return;
    const dx = e.clientX - state.lastX;
    const dy = e.clientY - state.lastY;
    
    if (e.buttons === 1 && !(e.altKey || e.metaKey)) {
      const mapPoint = getCursorOnMap(e.clientX, e.clientY);
      if (mapPoint) {
        const currentVector = inverseEqualEarthToVector(mapPoint.x, mapPoint.y);
        const prevVector = state.lastVector;
        if (currentVector && prevVector) {
          const axis = new THREE.Vector3().crossVectors(prevVector, currentVector);
          const dot = Math.max(-1, Math.min(1, prevVector.dot(currentVector)));
          const angle = Math.acos(dot);
          if (angle > 0.0001) {
            axis.normalize();
            const qRot = new THREE.Quaternion();
            qRot.setFromAxisAngle(axis, angle);
            sphereQuaternion.premultiply(qRot);
            sphereQuaternion.normalize();
            rotationMatrix.makeRotationFromQuaternion(sphereQuaternion);
          }
        }
        state.lastVector = currentVector;
      }
    } else if (e.buttons === 1 && (e.altKey || e.metaKey)) {
      const rollAxis = new THREE.Vector3(1, 0, 0); 
      const sensitivity = 0.01;
      const angle = (dx + dy) * sensitivity;
      const qRot = new THREE.Quaternion();
      qRot.setFromAxisAngle(rollAxis, angle);
      sphereQuaternion.premultiply(qRot);
      sphereQuaternion.normalize();
      rotationMatrix.makeRotationFromQuaternion(sphereQuaternion);
    } else if (e.buttons === 2) {
      if (domElement.clientHeight > 0) {
        const worldHeight = (camera.top - camera.bottom);
        const pxToWorld = worldHeight / domElement.clientHeight;
        state.pan.x -= dx * pxToWorld;
        state.pan.y += dy * pxToWorld;
      }
    }
    state.lastX = e.clientX;
    state.lastY = e.clientY
  });

  window.addEventListener('mouseup', () => {
    state.isDragging = false;
    domElement.style.cursor = 'default';
  });

  const resizeObserver = new ResizeObserver(() => {
    if (domElement.clientWidth > 0) {
      renderer.setSize(domElement.clientWidth, domElement.clientHeight, false);
    }
  });
  resizeObserver.observe(domElement);

  const render = function() {
      if (domElement.clientHeight === 0) return;
      const aspect = domElement.clientWidth / domElement.clientHeight;
      const frustumSize = 500 / state.zoom; 
      camera.left = -frustumSize * aspect / 2 + state.pan.x;
      camera.right = frustumSize * aspect / 2 + state.pan.x;
      camera.top = frustumSize / 2 + state.pan.y;
      camera.bottom = -frustumSize / 2 + state.pan.y;
      camera.updateProjectionMatrix();
      renderer.render(scene, camera);
  };

  return render;
};

export { SetupCameraAndRenderer, SetupCameraControls };
