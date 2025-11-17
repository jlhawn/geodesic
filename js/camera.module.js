import * as THREE from "./three.module.js";

const RADIUS = 1;
const KEY_DELTA = 4;    // pixel-equivalent per frame
const FIELD_OF_VIEW = 50;

const viewSize = new THREE.Vector2();

function SetupCameraAndRenderer() {
  viewSize.set(window.innerWidth, window.innerHeight);

  const camera = new THREE.PerspectiveCamera( FIELD_OF_VIEW, viewSize.x/viewSize.y, 0.1, 1000 );

  const renderer = new THREE.WebGLRenderer();
  renderer.setSize( viewSize.x, viewSize.y );
  document.body.appendChild( renderer.domElement );

  const stats = new Stats();
  document.body.appendChild( stats.dom );

  const setWindowSize = function() {
    viewSize.set(window.innerWidth, window.innerHeight);

    camera.aspect = viewSize.x / viewSize.y;
    camera.updateProjectionMatrix();
    renderer.setSize(viewSize.x, viewSize.y);
  };
  window.addEventListener("resize", setWindowSize, false);

  SetupCameraControls(camera, renderer.domElement);

  return { camera, renderer };
};

function SetupCameraControls(camera, domElement) {
  const raycaster = new THREE.Raycaster();
  const q = new THREE.Quaternion();

  // Normalized from [-1, 1] for camera field of view.
  const pointerInCamera = new THREE.Vector2();
  // Same as pointerInCamera but scaled by (width,height).
  const pointerInView = new THREE.Vector2();
  const prevPointerinView = new THREE.Vector2();

  domElement.addEventListener("contextmenu", function(e) {
      e.preventDefault();
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
    const nx = ( x / viewSize.x )*2 - 1;
    const ny = (-y / viewSize.y)*2 + 1;
    pointerInCamera.set(nx, ny);
    prevPointerinView.copy(pointerInView);
    pointerInView.set(nx*viewSize.x, ny*viewSize.y);
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

  // Because a small movement at the edge can translate to a big angle,
  // we want to dramatically reduce that angle so that the inertial
  // translation doesn't seem too fast.
  const _sphereEdgeInertiaAngleScalar = 0.25;
  const _epsilon = 1e-6;

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
      // This happens when the pointer is held down, moves off the sphere,
      // and then back onto the sphere.
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

  const _ndcCenter = new THREE.Vector2(0, 0); // center of the view in NDC
  const _centerOnSphere = new THREE.Vector3(); // intersection at center ray

  // degrees per full-view-width drag; tweak to taste
  const ROLL_SPEED = 180 * Math.PI/180;
  const MIN_TILT = 0;
  const MAX_TILT = 60 * Math.PI/180;
  const _forward = new THREE.Vector3();

  const onMoveRotate = function(e) {
    // TODO: IGNORE FOR NOW.
    const dx = pointerInView.x - prevPointerinView.x;
    const dy = pointerInView.y - prevPointerinView.y;

    // Need to find the point on the sphere at the center of the view
    // using the raycaster and rotate around the axis that goes from
    // the center of the sphere (origin) to that point.
    // Cast a ray through the center of the view (0,0 in NDC).
    // This gives us the point on the unit sphere that lies under
    // the exact center of the camera's view.
    raycaster.setFromCamera(_ndcCenter, camera);
    if (!raycaster.ray.intersectSphere(_unitSphere, _centerOnSphere)) {
      // If the center ray doesn't hit the sphere (e.g. zoomed far out
      // or looking away), there's no meaningful roll axis.
      return;
    }

    if (Math.abs(dx) > 0) {
      // Rotation axis: from origin to the intersection point at the center.
      _axis.copy(_centerOnSphere).normalize();

      // Map horizontal motion to a small roll angle. dx is in "view units"
      // (scaled NDC), so dividing by viewSize.x approximates "fraction
      // of screen width" and ROLL_SPEED scales that to radians.
      const angle = (dx / viewSize.x) * ROLL_SPEED;

      // For a right-handed system, positive rotation is CCW when looking
      // along the axis. If you want drag-right to feel like clockwise roll
      // from the viewer's perspective, you may need to flip the sign here.
      q.setFromAxisAngle(_axis, -angle);

      camera.applyQuaternion(q);
      camera.position.applyQuaternion(q);
    }

    if (Math.abs(dy) > 0) {
      // Our current tilt can be determined by finding the angle between
      // the vector from the origin to the _centerOnSphere (done) and the
      // negative vector in the direction the camera is looking.
      camera.getWorldDirection(_forward);
      let currentTilt = Math.acos(_centerOnSphere.dot(_forward.negate()));
      if (isNaN(currentTilt)) {
        currentTilt = 0;
      }

      // Need to tilt the camera so that it continues to face the point on the
      // sphere at the center of the screen but its position should be rotated
      // up or down around a horizontal (from the camera's perspective) axis
      // that is tangent to that point on the sphere.
      // (1, 0, 0) is right from the camera's local perspective. Use its
      // perspective quaternion to move it to the actual rightward direction.
      _axis.set(1, 0, 0).applyQuaternion(camera.quaternion).normalize();

      let angle = (dy / viewSize.y) * ROLL_SPEED;
      if (currentTilt + angle > MAX_TILT) {
        angle = MAX_TILT - currentTilt;
      } else if (currentTilt + angle < MIN_TILT) {
        angle = MIN_TILT - currentTilt;
      }

      if (Math.abs(angle) < _epsilon) {
        // Tilt should have no effect.
        return;
      }

      q.setFromAxisAngle(_axis, angle);

      camera.position.sub(_centerOnSphere); // move pivot to origin.
      camera.position.applyQuaternion(q);   // rotate around right axis through origin.
      camera.position.add(_centerOnSphere); // move back.

      camera.applyQuaternion(q); // rotate camera to face pivot point.
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
    if (!calculatePointerOnSphere(e)) {
      // The pointer has moved off the sphere.
      // Set it to the origin.
      pointerOnSphere.copy(_origin);
    }

    domElement.addEventListener('pointermove', onPointerMove);
    domElement.setPointerCapture(e.pointerId);

    domElement.addEventListener('pointerup', function(e) {
      anchorPoint.copy(_origin);
      activateInertia();
      domElement.removeEventListener('pointermove', onPointerMove);
      domElement.releasePointerCapture(e.pointerId);
    });
  });

  const _pivotToCamera = new THREE.Vector3();
  const MIN_DISTANCE = 0.2;
  const MAX_DISTANCE = 5;

  domElement.addEventListener('wheel', function(e) {
    e.preventDefault();
    cancelInertia();

    const factor = 1 + e.deltaY*0.002;

    raycaster.setFromCamera(_ndcCenter, camera);
    if (!raycaster.ray.intersectSphere(_unitSphere, _centerOnSphere)) {
      // If the center ray doesn't hit the sphere (e.g. zoomed far out
      // or looking away), there's no meaningful roll axis.
      return;
    }

    _pivotToCamera.subVectors(camera.position, _centerOnSphere);
    const dist = _pivotToCamera.length();
    const newDist = THREE.MathUtils.clamp(dist*factor, MIN_DISTANCE, MAX_DISTANCE);

    _pivotToCamera.normalize();
    camera.position.copy(_centerOnSphere).addScaledVector(_pivotToCamera, newDist);
  });

  // For mobile web, need to prevent swiping on the screen from scrolling
  // around the page.
  const preventDefault = function(e) { e.preventDefault(); };
  ["touchstart", "touchend", "touchmove", "touchcancel"].forEach(function(eventType) {
    domElement.addEventListener(eventType, preventDefault);
  });
};

export { SetupCameraAndRenderer, SetupCameraControls };
