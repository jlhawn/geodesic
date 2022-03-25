const RADIUS = 1;

// Length of simulated second in real milliseconds.
const SECOND = 0.5;
const MINUTE = 60*SECOND;
const HOUR = 60*MINUTE;
const DAY = 24*HOUR;
const YEAR = 365*DAY;
const AXIAL_TILT = 23.4*Math.PI/180;

// The equation for the number of cells is 5*2^(2*n+1)+2
//  0 ->         12
//  1 ->         42
//  2 ->        162
//  3 ->        642
//  4 ->      2,562
//  5 ->     10,242
//  6 ->     40,962
//  7 ->    163,842
//  8 ->    655,362
//  9 ->  2,621,442
// 10 -> 10,485,762
const SUBDIVISIONS = 6;

if (SUBDIVISIONS < 1) {
    throw "SUBDIVISIONS must be greater than zero.";
}

const NUM_CELLS = 5 * (1 << (2*SUBDIVISIONS + 1)) + 2;
const EARTH_SURFACE_AREA_SQKM = 5.1007e8;
const CELL_AREA_SQKM = EARTH_SURFACE_AREA_SQKM / NUM_CELLS;
// Assuming each cell is a regular hexagon, the equation for the area of a
// hexagon is A = 3*√3/2 * a^2 where a is a side length. Solving for a side
// length gives: a = √(2√3*A/9). The diameter is double a side length.
const CELL_AVG_DIAMETER_KM = 2*Math.sqrt(2*Math.sqrt(3)*CELL_AREA_SQKM/9);

const simplex = new SimplexNoise(Math.random());

function lerp(aVal, aMin, aMax, bMin, bMax) {
    var aRange = aMax - aMin;
    var valRatio = (aVal - aMin)/aRange;

    var bRange = bMax - bMin;
    var bVal = valRatio*bRange + bMin;
    return bVal;
}

function terrainSimplex(x, y, z) {
    // Use 4 different harmonics.
    const e0 = simplex.noise3D(x/2, y/2, z/2)*2;
    const e1 = simplex.noise3D(x, y, z);
    const e2 = simplex.noise3D(2*x, 2*y, 2*z)/2;
    const e3 = simplex.noise3D(3*x, 3*y, 3*z)/3;
    const e4 = simplex.noise3D(4*x, 4*y, 4*z)/4;

    return (e0+e1+e2+e3+e4)/(2 + 1 + 1/2 + 1/3 + 1/4);
}

function simplexColor(x, y, z) {
    const height = terrainSimplex(x, y, z);
    let r = height;
    let g = height;
    let b = height;

    if (height < 0) {
        r = lerp(height, -1, 0, 0.1, 0.25);
        b = lerp(height, -1, 0, 0.25, 1);
        g = lerp(height, -1, 0, 0.1, 0.25);
    } else {
        r = lerp(height, 0, 1, 0.2, 0.5);
        g = lerp(height, 0, 1, 0.5, 1);
        b = lerp(height, 0, 1, 0.2, 0.5);
    }

    return {r, g, b};
}

function runApp() {
    let start = performance.now();

    const grid = new Grid(SUBDIVISIONS);

    console.log(`completed grid construction in ${performance.now()-start}ms`);
    start = performance.now();

    const numPentagons = 12;
    const numHexagons = NUM_CELLS - numPentagons;
    const numTriangles = 4*numPentagons + 5*numHexagons;

    const vertices = new Float32Array(9*numTriangles);
    const normals = new Float32Array(9*numTriangles);
    const colors = new Float32Array(12*numTriangles);

    const color = new THREE.Color();
    const n = new THREE.Vector3();

    i = 0; // Triangle counter.

    for (const gridCell of grid) {
        {
            const {x, y, z} = gridCell.centerVertex;
            let {r, g, b} = simplexColor(x, y, z);
            if (gridCell.isAlongIcosahedronEdge) {
                r = lerp(r, 0, 1, 0.1, 1);
                g = lerp(g, 0, 1, 0.1, 1);
                b = lerp(b, 0, 1, 0.1, 1);
            }
            color.setRGB(r, g, b);
        }

        gridCell.faceTriangles.forEach(function(triangle) {
            const { a, b, c } = triangle;

            // Set vertex positions.
            vertices.set([a.x, a.y, a.z], i*9);
            vertices.set([b.x, b.y, b.z], i*9+3);
            vertices.set([c.x, c.y, c.z], i*9+6);

            // Get normal vector.
            triangle.getNormal(n);

            // One for each vertex.
            normals.set([n.x, n.y, n.z], i*9);
            normals.set([n.x, n.y, n.z], i*9+3);
            normals.set([n.x, n.y, n.z], i*9+6);

            // Set color for each vertex to be the same.
            colors.set([color.r, color.g, color.b, 1], i*12);
            colors.set([color.r, color.g, color.b, 1], i*12+4);
            colors.set([color.r, color.g, color.b, 1], i*12+8);

            i++; // IMPORTANT! Increment triangle counter.
        });
    }

    console.log(`completed vertex positions/normals/colors in ${performance.now()-start}ms`);

    function disposeArray() { this.array = null; };

    const geometry = new THREE.BufferGeometry();

    geometry.setAttribute('position', new THREE.BufferAttribute(vertices, 3).onUpload(disposeArray));
    geometry.setAttribute('normal', new THREE.BufferAttribute(normals, 3).onUpload(disposeArray));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 4).onUpload(disposeArray));

    geometry.computeBoundingSphere();

    const material = new THREE.MeshLambertMaterial({vertexColors: true});

    const mesh = new THREE.Mesh(geometry, material);

    const scene = new THREE.Scene();
    scene.add( new THREE.AmbientLight(0x404040) );
    
    const { camera, renderer, stats } = SetupCameraAndRenderer();

    scene.add(mesh);

    const X_AXIS = new THREE.Vector3(1, 0, 0);
    const Y_AXIS = new THREE.Vector3(0, 1, 0);
    const Z_AXIS = new THREE.Vector3(0, 0, 1);
    var q = new THREE.Quaternion();

    const SUN = new THREE.DirectionalLight(0xffffff);
    SUN.position.copy(X_AXIS);
    scene.add( SUN );

    // Move the camera away from the origin and rotate it to point back
    // to the origin (it starts out pointing in the negative-z) direction but
    // we want it to point in the positive-y direction.
    camera.position.y = -10;
    q.setFromAxisAngle(X_AXIS, Math.PI/2);
    camera.applyQuaternion(q);

    // Next we want to rotate the camera pi/4 radians around the z-axis to
    // have a more interesting perspective.
    q.setFromAxisAngle(Z_AXIS, Math.PI/4);
    camera.position.applyQuaternion(q);
    camera.applyQuaternion(q);

    SetupCameraControls(camera, renderer.domElement);

    const SUN_START_POS = X_AXIS.clone();
    const TILT_AXIS = Z_AXIS.clone();

    const timer = document.createElement("div");
    timer.style.position = "fixed";
    timer.style.top = "50px";
    timer.style.left = "5px";
    timer.style.color = "cyan";
    timer.style.fontFamily = "sans-serif";
    document.body.appendChild(timer);

    var t0 = 0;
    function animate(t) {
        requestAnimationFrame( animate );

        const years = Math.floor(t / YEAR);
        const days = Math.floor((t%YEAR) / DAY);
        const hours = Math.floor((t%DAY) / HOUR);
        const minutes = Math.floor((t%HOUR) / MINUTE);

        timer.innerText = `${years} years, ${days} days, ${hours} hours, ${minutes} minutes`;

        // Simulate the axial tilt of the earth.
        // To do this, we tilt the sun towards or away from the Z_AXIS by a
        // tiny bit with each time step. Assuming t=0 is the spring equinox,
        // the angle between the sun's position vector and the equator follows
        // a sine function which ossilates between 23.4 degrees at its max at
        // the summer solstice (tilted towards the Z_AXIS) to negative 23.4 at
        // its minimum (tilted towards the negative Z_AXIS)
        // 
        // Equation for axial tilt (from perspective of the sun):
        //
        //     AXIAL_TILT*sin(2*PI*t/YEAR)
        //
        let theta = AXIAL_TILT*Math.sin(2*Math.PI*(t % YEAR)/YEAR);
        q.setFromAxisAngle(Y_AXIS, -theta); // Clockwise (negative) rotation.
        SUN.position.copy(SUN_START_POS).applyQuaternion(q);

        // Simulate earth rotating around its axis.
        // To avoid modifying the sphere mesh, we will rotate the sun
        // around the Z_AXIS.
        // NOTE: while the earth actually rotates counter-clockwise around its
        // axis, the sun appears to rotate clockwise so we negate this angle.
        theta = -(t % DAY)*2*Math.PI/DAY;
        q.setFromAxisAngle(Z_AXIS, theta);
        // Note that we're changing the *position* of the light source and not
        // the orientation which would have no effect.
        SUN.position.applyQuaternion(q);

        renderer.render( scene, camera );
        stats.update();
    }
    requestAnimationFrame( animate );
}

function SetupCameraAndRenderer() {
    const camera = new THREE.PerspectiveCamera( 15, window.innerWidth / window.innerHeight, 0.1, 1000 );

    const renderer = new THREE.WebGLRenderer();
    renderer.setSize( window.innerWidth, window.innerHeight );
    document.body.appendChild( renderer.domElement );

    const stats = new Stats();
    document.body.appendChild( stats.dom );

    const setWindowSize = function() {
        camera.aspect = window.innerWidth / window.innerHeight;
        camera.updateProjectionMatrix();

        renderer.setSize(window.innerWidth, window.innerHeight);
    };
    window.addEventListener('resize', setWindowSize, false);

    return {
        camera: camera,
        renderer: renderer,
        stats: stats,
    }
}

function SetupCameraControls(camera, domElement) {
    // The camera is always looking towards its local negative-z axis.
    // Up is its local positive-y axis.
    // Right is its local positive-x axis.
    const forward = new THREE.Vector3(0, 0, -1);
    const      up = new THREE.Vector3(0, 1,  0);
    const   right = new THREE.Vector3(1, 0,  0);

    // We can apply the camera's current orientation quaternion to each of
    // these vectors to determine their real directions. If it's still the
    // default quaternion (w=1, x=0, y=0, z=0) then this has no effect.
    [forward, up, right].forEach((v) => v.applyQuaternion(camera.quaternion));

    // Now that we know which way is right and which way is up, we can tilt the
    // camera up and down around the origin when the mouse is dragged up or
    // down by rotating its position about the right vector axis and rotating
    // the camera orientation in the same way. Likewise, when the mouse is
    // dragged left or right, we can rotate the camera about the up vector axis
    // and rotate teh camera orientation in the same way.

    const raycaster = new THREE.Raycaster();
    const pointer = new THREE.Vector2();

    const target = new THREE.Vector2();
    const tangentPlane = new THREE.Plane(forward.clone().negate(), -RADIUS);

    const q = new THREE.Quaternion();

    const move = function(e) {
        let dx = e.clientX - pointer.x; // Horizontal delta.
        let dy = pointer.y - e.clientY; // Vertical delta.

        pointer.x = e.clientX;
        pointer.y = e.clientY;

        // Normalize the deltas to a coordinate within the camera's field of
        // view (between -1 to 1) relative to the center.
        dx = (dx / domElement.clientWidth) * 2;
        dy = (dy / domElement.clientHeight) * 2;

        target.set(dx, dy);

        // Now that we have the delta between where the pointer was and where
        // it is now, we can cast a ray from the camera to a point on a plane
        // which is tangent to the surface of the sphere and facing the camera.
        // The means that the plane's normal vector will need to be the
        // negative of the camera's forward vector.
        tangentPlane.normal.copy(forward).negate();
        raycaster.setFromCamera(target, camera);
        let intersection = new THREE.Vector3();
        raycaster.ray.intersectPlane(tangentPlane, intersection);
        intersection.normalize();


        // Next we find the horizontal and vertical angles between these two
        // points. The horizontal angle we will denote as theta. The vertical
        // angle we will denote as phi.
        let upPlane = new THREE.Plane(up);
        let rightPlane = new THREE.Plane(right);
        let proj = new THREE.Vector3();
        upPlane.projectPoint(intersection, proj);
        let theta = tangentPlane.normal.angleTo(proj);
        theta = dx > 0 ? -theta : theta;
        rightPlane.projectPoint(intersection, proj);
        let phi = tangentPlane.normal.angleTo(proj);
        phi = dy > 0 ? phi : -phi;

        q.setFromAxisAngle(up, theta);
        camera.applyQuaternion(q);
        camera.position.applyQuaternion(q);

        // Up vector has not changed but forward and right have.
        forward.applyQuaternion(q);
        right.applyQuaternion(q);

        q.setFromAxisAngle(right, phi);
        camera.applyQuaternion(q);
        camera.position.applyQuaternion(q);

        // Right vector has not changed but forward and up have.
        forward.applyQuaternion(q);
        up.applyQuaternion(q);
    };

    domElement.addEventListener('pointerdown', function(e) {
        pointer.x = e.clientX;
        pointer.y = e.clientY;

        domElement.addEventListener('pointermove', move);
        domElement.setPointerCapture(e.pointerId);

        domElement.addEventListener('pointerup', function(e) {
            domElement.removeEventListener('pointermove', move);
            domElement.releasePointerCapture(e.pointerId);
        });
    });

    domElement.addEventListener('wheel', function(e) {
        e.preventDefault();

        let factor = 1 + e.deltaY*0.002;

        camera.position.multiplyScalar(factor);

        let distanceSq = camera.position.lengthSq();
        if (distanceSq < 4) { // 2 squared.
            camera.position.normalize();
            camera.position.multiplyScalar(2);
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
}
