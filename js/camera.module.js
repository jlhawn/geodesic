import * as THREE from "./three.module.js";

const RADIUS = 1;
const KEY_DELTA = 4;    // pixel-equivalent per frame

function SetupCameraControls(camera, domElement) {
  const raycaster = new THREE.Raycaster();
  const q = new THREE.Quaternion();

  let viewWidth = domElement.clientWidth;
  let viewHeight = domElement.clientHeight;
  const pointerInCamera = new THREE.Vector2();

  domElement.addEventListener("contextmenu", function(e) {
      e.preventDefault();
  });

  window.addEventListener("resize", function() {
    viewWidth = domElement.clientWidth;
    viewHeight = domElement.clientHeight;
  });

  // Convert pixel coordinates (x, y) into
  // Normalized Device Coordinates (NDC) where:
  //   NDC.x = -1 (left)  to +1 (right)
  //   NDC.y = -1 (bottom) to +1 (top)
  //
  // This format is required by THREE.Raycaster.setFromCamera()
  // and matches clip space.
  const normalizeToCamera = function(x, y) {
    // Convert pixel to NDC:
    //   x_ndc =  (pixelX / w) * 2 - 1
    //   y_ndc = -(pixelY / h) * 2 + 1
    //
    // Notes:
    // - The minus in the Y formula flips the browser’s top-origin
    //   coordinate system into OpenGL-style bottom-origin coordinates.
    const nx = ( x / viewWidth )*2 - 1;
    const ny = (-y / viewHeight)*2 + 1;
    pointerInCamera.set(nx, ny);
    return pointerInCamera;
  }

  const _origin = new THREE.Vector3();
  const _unitSphere = { center: _origin, radius: 1 };
  const anchorPoint = new THREE.Vector3();
  const pointerOnSphere = new THREE.Vector3();

  // Temps for rotation:
  const _from = new THREE.Vector3();
  const _to   = new THREE.Vector3();
  const _axis = new THREE.Vector3();
  const _epsilon = 1e-6;

  const calculatePointerOnSphere = function(e) {
    // Convert from pixel coords to normalized device coordinates (NDC).
    const ndc = normalizeToCamera(e.clientX, e.clientY);

    // We want to find where on the unit sphere the cursor was pressed down.
    // To do this, we use the raycaster to find where the ray from the
    // camera through the normalized device coordinates of the click passes
    // through a unit sphere in the scene. If it doesn't intersect the unit
    // sphere then the pointerdown event can be ignored.

    // Build a ray from the camera through the NDC point.
    raycaster.setFromCamera(ndc, camera);

    // Test the ray against the unit sphere centered at the origin.
    return raycaster.ray.intersectSphere(_unitSphere, pointerOnSphere) !== null;
  };

  const _inertia = { active: false };
  const _minInertialAngle = 0.1 * Math.PI/180;
  const _inertiaDecay = 0.96;

  const inertialTranslate = function() {
    q.setFromAxisAngle(_inertia.axis, -_inertia.angle);

    camera.applyQuaternion(q);
    camera.position.applyQuaternion(q);

    _inertia.angle = _inertiaDecay * _inertia.angle;
    if (_inertia.angle >= _minInertialAngle) {
      _inertia.timeoutID = setTimeout(inertialTranslate, _inertia.dt);
    } else {
      cancelInertia();
    }
  };

  const activateInertia = function(angleScalar=1) {
    if (!(_inertia.active && _inertia.timeoutID === undefined)) {
      return // Never active or already running.
    }
    _inertia.dt = performance.now() - _inertia.lastTranslateAt;
    _inertia.angle *= angleScalar;
    _inertia.timeoutID = setTimeout(inertialTranslate, _inertia.dt);
  };

  const cancelInertia = function() {
    if (!_inertia.active) {
      return; // Already cancelled.
    }
    if (_inertia.timeoutID !== undefined) {
      clearTimeout(_inertia.timeoutID);
    }
    _inertia.active = false;
    _inertia.timeoutID = undefined;
  };

  const _sphereEdgeInertiaAngleScalar = 0.25;

  const onMoveTranslate = function(e) {
    if (pointerOnSphere.equals(_origin)) {
      // The pointer is not on the sphere. Zero out the anchor point.
      anchorPoint.copy(_origin);
      activateInertia(_sphereEdgeInertiaAngleScalar);
      return;
    }

    cancelInertia();

    if (anchorPoint.equals(_origin)) {
      // The current point on the sphere should be our new anchor point.
      // We need one more 'pointermove' event to decide how to translate.
      anchorPoint.copy(pointerOnSphere);
      return;
    }

    // We want to move the camera so that, from its perspective,
    // the anchorPoint is now at pointerOnSphere.

    // Normalize the two points to ensure they lie on the unit sphere.
    _from.copy(anchorPoint).normalize();
    _to.copy(pointerOnSphere).normalize();

    // Compute the angle between them, clamped to [-1, 1] for numerical safety.
    const dot = THREE.MathUtils.clamp(_from.dot(_to), -1, 1);
    const angle = Math.acos(dot);

    // If the angle is extremely small, skip to avoid jitter.
    if (angle < _epsilon) {
      return;
    }

    // Rotation axis is from × to.
    _axis.crossVectors(_from, _to);

    // If axis is near zero, the points are (almost) colinear; skip.
    if (_axis.lengthSq() < _epsilon) {
      return;
    }
    _axis.normalize();

    // Build quaternion that rotates _from -> _to in world space.
    // To get the same visual effect by moving the camera, apply the inverse
    // (i.e. rotate camera by -angle around the same axis).
    q.setFromAxisAngle(_axis, -angle);

    camera.applyQuaternion(q);
    camera.position.applyQuaternion(q);

    _inertia.active = true;
    _inertia.lastTranslateAt = performance.now(),
    _inertia.axis = _axis;
    _inertia.angle = angle;
  };

  let prevViewX;
  let prevViewY;

  const onMoveRotate = function(e) {
    // TODO: IGNORE FOR NOW.
    // Convert pointerInCamera to x,y coordinates about the center of the screen.
    const x = viewWidth * pointerInCamera.x;
    const y = viewHeight * pointerInCamera.y;

    if (prevViewX === undefined) {
      prevViewX = x;
    }
  };

  const onPointerMove = function(e) {
    if (!calculatePointerOnSphere(e)) {
      // The pointer has moved off the sphere.
      // Set it to the origin.
      pointerOnSphere.copy(_origin);
    }

    if (e.altKey || e.metaKey) {
      cancelInertia();
      onMoveRotate(e);
      // Set the anchor point to the point under the cursor (if on the sphere)
      // so that when the key is released that can be the new anchor point.
      anchorPoint.copy(pointerOnSphere);
    } else {
      onMoveTranslate(e);
    }
  };

  domElement.addEventListener('pointerdown', function(e) {
    cancelInertia();
    domElement.addEventListener('pointermove', onPointerMove);
    domElement.setPointerCapture(e.pointerId);

    domElement.addEventListener('pointerup', function(e) {
      anchorPoint.copy(_origin);
      activateInertia();
      domElement.removeEventListener('pointermove', onPointerMove);
      domElement.releasePointerCapture(e.pointerId);
    });
  });

  domElement.addEventListener('wheel', function(e) {
    e.preventDefault();
    cancelInertia();

    let factor = 1 + e.deltaY*0.002;

    camera.position.multiplyScalar(factor);

    let distanceSq = camera.position.lengthSq();
    if (distanceSq < 1.44) {
      camera.position.normalize();
      camera.position.multiplyScalar(1.2);
    } else if (distanceSq > 400) { // 20 squared.
      camera.position.normalize();
      camera.position.multiplyScalar(20);
    }
  });

  // For mobile web, need to prevent swiping on the screen from scrolling
  // around the page.
  const preventDefault = function(e) { e.preventDefault(); };
  ["touchstart", "touchend", "touchmove", "touchcancel"].forEach(function(eventType) {
    domElement.addEventListener(eventType, preventDefault);
  });
};

export { SetupCameraControls };
