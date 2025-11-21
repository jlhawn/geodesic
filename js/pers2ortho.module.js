import * as THREE from './three.module.js';
import { OrbitControls } from './OrbitControls.module.js';

let scene, renderer;
let perspCam, orthoCam, activeCamera, transitionCamera;
let controls;
let clock;
let isTransitioning = false;
let transitionStart = 0;
const transitionDuration = 2;

// Now used generically for "from" and "to"
const fromProj = new THREE.Matrix4();
const toProj = new THREE.Matrix4();

// Track which cameras we’re interpolating between
let transitionFromCam = null;
let transitionToCam = null;

function init() {
  renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setSize(window.innerWidth, window.innerHeight);
  document.body.appendChild(renderer.domElement);

  scene = new THREE.Scene();
  scene.background = new THREE.Color(0x222222);

  const grid = new THREE.GridHelper(100, 100);
  scene.add(grid);

  const geom = new THREE.BoxGeometry(1, 1, 1);
  const mat = new THREE.MeshNormalMaterial();

  const count = 200;
  const range = 100;

  for (let i = 0; i < count; i++) {
    const m = new THREE.Mesh(geom, mat);
    m.position.set(
      THREE.MathUtils.randFloatSpread(range),
      THREE.MathUtils.randFloatSpread(range * 0.2),
      THREE.MathUtils.randFloatSpread(range)
    );
    m.scale.setScalar(THREE.MathUtils.randFloat(2, 4));
    scene.add(m);
  }

  const aspect = window.innerWidth / window.innerHeight;

  perspCam = new THREE.PerspectiveCamera(50, aspect, 0.1, 1000);
  perspCam.position.set(10, 10, 20);
  perspCam.lookAt(0, 0, 0);
  perspCam.updateProjectionMatrix();

  orthoCam = new THREE.OrthographicCamera(-10, 10, 10, -10, 0.1, 1000);
  orthoCam.position.copy(perspCam.position);
  orthoCam.quaternion.copy(perspCam.quaternion);
  matchOrthoToPerspective(orthoCam, perspCam, 20);

  transitionCamera = new THREE.Camera();

  activeCamera = perspCam;

  controls = new OrbitControls(activeCamera, renderer.domElement);
  controls.enableDamping = true;
  controls.target.set(0, 0, 0);

  clock = new THREE.Clock();

  window.addEventListener('resize', onWindowResize);

  const btn = document.getElementById('btn');
  btn.addEventListener('click', () => {
    if (isTransitioning) return;

    if (activeCamera === perspCam) {
      startPerspectiveToOrtho();
      // optional: btn.textContent = 'Switch to Perspective';
    } else {
      startOrthoToPerspective();
      // optional: btn.textContent = 'Switch to Orthographic';
    }
  });
}

function matchOrthoToPerspective(ortho, persp, referenceDistance) {
  const aspect = renderer.domElement.clientWidth / renderer.domElement.clientHeight;
  const fovRad = THREE.MathUtils.degToRad(persp.fov);
  const h = 2 * referenceDistance * Math.tan(fovRad / 2);
  const w = h * aspect;

  ortho.left = -w / 2;
  ortho.right = w / 2;
  ortho.top = h / 2;
  ortho.bottom = -h / 2;
  ortho.near = persp.near;
  ortho.far = persp.far;
  ortho.updateProjectionMatrix();
}

function lerpMatrix4(out, m1, m2, t) {
  const e = out.elements;
  const a = m1.elements;
  const b = m2.elements;
  for (let i = 0; i < 16; i++) e[i] = THREE.MathUtils.lerp(a[i], b[i], t);
}

// perspective -> orthographic
function startPerspectiveToOrtho() {
  orthoCam.position.copy(perspCam.position);
  orthoCam.quaternion.copy(perspCam.quaternion);
  matchOrthoToPerspective(orthoCam, perspCam, 20);

  transitionFromCam = perspCam;
  transitionToCam = orthoCam;

  fromProj.copy(perspCam.projectionMatrix);
  toProj.copy(orthoCam.projectionMatrix);

  isTransitioning = true;
  transitionStart = clock.getElapsedTime();
}

// orthographic -> perspective
function startOrthoToPerspective() {
  // Align perspective camera to current ortho pose
  perspCam.position.copy(orthoCam.position);
  perspCam.quaternion.copy(orthoCam.quaternion);
  perspCam.updateMatrixWorld(true);
  // aspect is already updated in resize handler

  transitionFromCam = orthoCam;
  transitionToCam = perspCam;

  fromProj.copy(orthoCam.projectionMatrix);
  toProj.copy(perspCam.projectionMatrix);

  isTransitioning = true;
  transitionStart = clock.getElapsedTime();
}

function onWindowResize() {
  const aspect = window.innerWidth / window.innerHeight;

  perspCam.aspect = aspect;
  perspCam.updateProjectionMatrix();

  matchOrthoToPerspective(orthoCam, perspCam, 20);

  renderer.setSize(window.innerWidth, window.innerHeight);
}

function switchControlsTo(camera) {
  const t = controls.target.clone();
  controls.dispose();
  controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.target.copy(t);
}

function animate() {
  requestAnimationFrame(animate);

  const tNow = clock.getElapsedTime();

  if (isTransitioning) {
    const t = THREE.MathUtils.clamp(
      (tNow - transitionStart) / transitionDuration,
      0,
      1
    );

    // Interpolate between whichever cameras are active for this transition
    transitionCamera.position.lerpVectors(
      transitionFromCam.position,
      transitionToCam.position,
      t
    );
    transitionCamera.quaternion.slerpQuaternions(
      transitionFromCam.quaternion,
      transitionToCam.quaternion,
      t
    );
    transitionCamera.updateMatrixWorld(true);

    // And between the stored projection matrices
    lerpMatrix4(transitionCamera.projectionMatrix, fromProj, toProj, t);
    transitionCamera.projectionMatrixInverse
      .copy(transitionCamera.projectionMatrix)
      .invert();

    renderer.render(scene, transitionCamera);

    if (t >= 1.0) {
      isTransitioning = false;
      activeCamera = transitionToCam;
      switchControlsTo(activeCamera);
    }
  } else {
    controls.update();
    renderer.render(scene, activeCamera);
  }
}

export function run() {
  init();
  animate();
}
