const RADIUS = 1;
const HOUR = 1000; // Length of simulated hour in real milliseconds.
const SUBDIVISIONS = 4;

function projectOntoUnitSphere(triangle) {
    let { a, b, c } = triangle;
    a.normalize();
    b.normalize();
    c.normalize();
}

function circumcenter(triangle) {
    const { a, b, c } = triangle;

    const p1p2 = a.clone().sub(b);
    const p2p1 = b.clone().sub(a);
    const p1p3 = a.clone().sub(c);
    const p3p1 = c.clone().sub(a);
    const p2p3 = b.clone().sub(c);
    const p3p2 = c.clone().sub(b);

    const p1p2Xp2p3 = p1p2.clone().cross(p2p3);

    const alpha = p2p3.lengthSq() * p1p2.dot(p1p3) / (2 * p1p2Xp2p3.lengthSq());
    const beta = p1p3.lengthSq() * p2p1.dot(p2p3) / (2 * p1p2Xp2p3.lengthSq());
    const gamma = p1p2.lengthSq() * p3p1.dot(p3p2) / (2 * p1p2Xp2p3.lengthSq());

    return a.clone().multiplyScalar(alpha).add(b.clone().multiplyScalar(beta)).add(c.clone().multiplyScalar(gamma));
}

function subdivideTriangle(triangle, N, callback) {
    if (N == 0) {
        callback(triangle.clone());
        return;
    }

    const { a, b, c } = triangle;

    // It will be helpful to orient our thinking such that a is at the 'top' of
    // the triangle, b is the bottom-left, and c is the bottom-right.
    //
    //         a
    //          /\
    //         /  \
    //      d /____\ f 
    //       /\    /\
    //      /  \  /  \
    //   b /____\/____\ c
    //          e
    //

    const d = a.clone().add(b).normalize();
    const e = b.clone().add(c).normalize();
    const f = c.clone().add(a).normalize();

    let t = new THREE.Triangle();
    t = t.set(a, d, f);
    subdivideTriangle(t, N-1, callback);
    t = t.set(d, b, e);
    subdivideTriangle(t, N-1, callback);
    t = t.set(e, f, d);
    subdivideTriangle(t, N-1, callback);
    t = t.set(f, e, c);
    subdivideTriangle(t, N-1, callback);
}

// subdivideTriangleOLD divides the given triangle into multiple smaller triangles
// with 'N' triangles along each side. The number of new smaller
// triangles grows as the square of 'N':
//     2->4, 3->9, 4->16, etc.
//
// NOTE: I stopped using this function because the size of the triangles varied
// too much. The recursive solution creates much better triangles because it
// projects the vertexes to the sphere with each iteration.
function subdivideTriangleOLD(triangle, N) {
    // NOTE: Points a, b, c make a counter-clockwise turn in the front-facing
    // direction. i.e., by right hand rule, the normal vector at b points in
    // the direction that the triangle's "front" side faces. This is important
    // to know because the rederer performs back-face culling by default.
    const { a, b, c } = triangle;

    // It will be helpful to orient our thinking such that a is at the 'top' of
    // the triangle, b is the bottom-left, and c is the bottom-right.

    const a0 = a.clone().divideScalar(N);
    const b0 = b.clone().divideScalar(N);
    const c0 = c.clone().divideScalar(N);

    const triangles = [];

    // Each vertex in our subdivided triangles will be some linear combination
    // of a0, b0, and c0 such that for each each vertex:
    //     a0*i + b0*j + c0*k = vertex
    // where i + j + k = N.

    // We'll be creating the subdivided triangles by going row-by row from top
    // to bottom of the triangle.

    // Each row begins with an upward-pointing triangle followed by a
    // downward-pointing triangle. This repeats 0 or more times and always ends
    // with an upward-pointing triangle.

    // For the pair of upward and downward facing triangles, there will be 4
    // vertices with 2 shared by both triangles:
    //     v1  ____  v4
    //       /\    /
    //      /  \  /
    //  v2 /____\/ v3
    //

    for (var i = N; i > 0; i--) {
        var j = N-i;

        var v1 = a0.clone().multiplyScalar(i).add(b0.clone().multiplyScalar(j));
        var v2 = v1.clone().sub(a0).add(b0);
        var v3 = v1.clone().sub(a0).add(c0);

        while (j > 0) {
            var v4 = v1.clone().sub(b0).add(c0);
            
            triangles.push(new THREE.Triangle(v1, v2, v3));
            triangles.push(new THREE.Triangle(v3.clone(), v4, v1.clone()));

            v1 = v4.clone();
            v2 = v3.clone();
            v3 = v1.clone().sub(a0).add(c0);

            j--;
        }

        triangles.push(new THREE.Triangle(v1, v2, v3));
    }

    return triangles;
}

function GeodesicTriangles(N, callback) {
    const a = 0.525731112119133606;
    const b = 0.850650808352039932;

    // 12 Vertices of an Icosahedron.

    const vA = new THREE.Vector3(-a, 0, b);
    const vB = new THREE.Vector3( a, 0, b);
    const vC = new THREE.Vector3(-a, 0, -b);
    const vD = new THREE.Vector3( a, 0, -b);
    const vE = new THREE.Vector3( 0, b, a);
    const vF = new THREE.Vector3( 0, b, -a);
    const vG = new THREE.Vector3( 0, -b, a);
    const vH = new THREE.Vector3( 0, -b, -a);
    const vI = new THREE.Vector3( b, a, 0);
    const vJ = new THREE.Vector3(-b, a, 0);
    const vK = new THREE.Vector3( b, -a, 0);
    const vL = new THREE.Vector3(-b, -a, 0);

    const vertices = [vA, vB, vC, vD, vE, vF, vG, vH, vI, vJ, vK, vL];

    /*

    These vertices start out pretty symmetrical, each being within a plane
    of two of the x, y, or z axes. We want to rotate it about the y axis such
    that vA becomes a unit vector in the direction of the z axis. The sine of
    that rotation angle happens to be 'a' and the cosine happens to be 'b'. The
    transformation matrix for this rotation is given by:

        [ b   0  a ]
        [ 0   1  0 ]
        [ -a  0  b ]

    After applying this rotation, vA is the "north pole" and vD is the
    "south pole".

    */
    const m = new THREE.Matrix3();
    m.set( b, 0, a,
           0, 1, 0,
          -a, 0, b );

    vertices.forEach(function(vec) {
        vec.applyMatrix3(m);
    });

    // Next, we'll make the 20 triangle faces of the icosahedron.
    // The order of these vertices matters in each triangle to establish which
    // side of the triangle is the "front" side according to the right hand
    // rule. If the "back" side faces the camera then it will not be rendered
    // by default, so all of the triangles here are created such that the
    // "front" sides face out from the base icosahedron.
    const triangles = [
        new THREE.Triangle( vA, vB, vE ),
        new THREE.Triangle( vE, vB, vI ),
        new THREE.Triangle( vE, vI, vF ),
        new THREE.Triangle( vF, vI, vD ),
        new THREE.Triangle( vD, vC, vF ),
        new THREE.Triangle( vF, vC, vJ ),
        new THREE.Triangle( vF, vJ, vE ),
        new THREE.Triangle( vE, vJ, vA ),
        new THREE.Triangle( vA, vJ, vL ),
        new THREE.Triangle( vL, vJ, vC ),
        new THREE.Triangle( vL, vC, vH ),
        new THREE.Triangle( vH, vC, vD ),
        new THREE.Triangle( vD, vK, vH ),
        new THREE.Triangle( vH, vK, vG ),
        new THREE.Triangle( vH, vG, vL ),
        new THREE.Triangle( vL, vG, vA ),
        new THREE.Triangle( vA, vG, vB ),
        new THREE.Triangle( vB, vG, vK ),
        new THREE.Triangle( vB, vK, vI ),
        new THREE.Triangle( vI, vK, vD ),
    ];

    triangles.forEach(function(triangle) {
        subdivideTriangle(triangle, N, callback);
    });
}

function vectorKey(v) {
    return [v.x.toFixed(5), v.y.toFixed(5), v.z.toFixed(5)].toString();
}

function vertexKeys(triangle) {
    const { a, b, c } = triangle;

    return [
        vectorKey(a),
        vectorKey(b),
        vectorKey(c),
    ];
}

function canonicalizeVertexOrder(firstVertexKey, triangle) {
    const {a, b, c} = triangle;
    if (vectorKey(a) == firstVertexKey) {
        return;
    }
    if (vectorKey(b) == firstVertexKey) {
        triangle.a = b;
        triangle.b = c;
        triangle.c = a;
        return;
    }
    if (vectorKey(c) != firstVertexKey) {
        throw `expected vertex c of triangle ${triangle.toString()} to be ${firstVertexKey} but it was ${vectorKey(c)}`;
    }
    triangle.a = c;
    triangle.b = a;
    triangle.c = b;
}

function findTriangleWithSecondVertexKey(triangles, vertexKey) {
    for (const triangle of triangles) {
        let secondKey = vectorKey(triangle.b);
        if (secondKey == vertexKey) {
            return triangle;
        }
    }
    throw `could not find triangle with 2nd vertex with key ${vertexKey} in ${triangles.toString()}`;
}

function veronoi(commonVertexKey, triangles) {
    // The given triangles all meet at a common vertex, making either
    // a pentagon (5 triangles) or a hexagon (6 triangles);
    // First, we canonicalize each of the triangles such that their first
    // vertex is this shared vertex.
    for (const triangle of triangles) {
        canonicalizeVertexOrder(commonVertexKey, triangle);
    }

    // Next, we'll want to arrange the ordering of these triangles such that
    // we go around the common vertex in a counter-clockwise direction. If we
    // chose the first triangle to be {a, b, c}, then there should be another
    // triangle whose vertices are {a, c, d}, then the next would have vertices
    // {a, d, e}, etc.
    //
    //            b _____ g
    //             /\    /\
    //            /  \  /  \
    //         c /___a\/____\ 
    //           \    /\    / f
    //            \  /  \  /
    //             \/____\/
    //              d     e

    let triangle = triangles[0];
    const orderedTriangles = [triangle];
    while (triangles.length > orderedTriangles.length) {
        triangle = findTriangleWithSecondVertexKey(triangles, vectorKey(triangle.c));
        orderedTriangles.push(triangle);
    }

    // Now that we have the ordered triangles, we can calculate the
    // circumcenter of each triangle.
    const circumcenters = orderedTriangles.map(function(triangle) {
        return circumcenter(triangle).normalize();
    });

    const finalTriangles = [new THREE.Triangle(
        circumcenters[0].clone(),
        circumcenters[1].clone(),
        circumcenters[2].clone(),
    )];
    for(let i = 2; i+1 < circumcenters.length; i++) {
        finalTriangles.push(new THREE.Triangle(
            circumcenters[0].clone(),
            circumcenters[i].clone(),
            circumcenters[i+1].clone(),
        ));
    }

    return finalTriangles;
}

function runApp() {
    const positions = [];
    const normals = [];
    const colors = [];

    const color = new THREE.Color();
    const n = new THREE.Vector3();

    const vertexMap = new Map();

    GeodesicTriangles(SUBDIVISIONS, function(triangle) {
        for (const key of vertexKeys(triangle)) {
            const existing = vertexMap.get(key);
            if (existing) {
                existing.push(triangle);
            } else {
                vertexMap.set(key, [triangle]);
            }
        }
    });

    for (let [vertex, triangles] of vertexMap.entries()) {
        // Make pentagon or Hexagon.
        if (triangles.length !== 5 && triangles.length !== 6) {
            throw `vertex ${vertex} has ${triangles.length} triangles`;
        }

        const veronoiTriangles = veronoi(vertex, triangles);

        {
            // Generate a random color.
            const r = Math.abs(triangles[0].a.x);
            const g = Math.abs(triangles[0].a.y);
            const b = Math.abs(triangles[0].a.z);
            color.setRGB(r, g, b);
        }

        veronoiTriangles.forEach(function(triangle) {
            const { a, b, c } = triangle;

            // Set vertex positions.
            positions.push(a.x, a.y, a.z);
            positions.push(b.x, b.y, b.z);
            positions.push(c.x, c.y, c.z);

            // Get normal vector.
            triangle.getNormal(n);

            // One for each vertex.
            normals.push(n.x, n.y, n.z);
            normals.push(n.x, n.y, n.z);
            normals.push(n.x, n.y, n.z);

            // Set color for each vertex to be the same.
            colors.push(color.r, color.g, color.b, 1);
            colors.push(color.r, color.g, color.b, 1);
            colors.push(color.r, color.g, color.b, 1);
        });
    }

    function disposeArray() { this.array = null; };

    const geometry = new THREE.BufferGeometry();

    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3).onUpload(disposeArray));
    geometry.setAttribute('normal', new THREE.Float32BufferAttribute(normals, 3).onUpload(disposeArray));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 4).onUpload(disposeArray));

    geometry.computeBoundingSphere();

    const material = new THREE.MeshLambertMaterial({vertexColors: true});

    const mesh = new THREE.Mesh(geometry, material);

    const scene = new THREE.Scene();
    scene.add( new THREE.AmbientLight(0x404040) );

    const SUN = new THREE.DirectionalLight(0xffffff);
    SUN.position.set(1, 0, 0);
    scene.add( SUN );
    
    const { camera, renderer, stats } = SetupCameraAndRenderer();

    scene.add(mesh);

    const X_AXIS = new THREE.Vector3(1, 0, 0);
    const Y_AXIS = new THREE.Vector3(0, 1, 0);
    const Z_AXIS = new THREE.Vector3(0, 0, 1);
    var q = new THREE.Quaternion();

    // Move the camera away from the origin and rotate it to point back
    // to the origin (it starts out pointing in the negative-z) direction but
    // we want it to point in the positive-y direction.
    camera.position.y = -10;
    q.setFromAxisAngle(X_AXIS, Math.PI/2);
    camera.applyQuaternion(q);

    SetupCameraControls(camera, renderer.domElement);

    // We want to tilt the sphere mesh to be 23.4 degrees away from the z axis.
    // We can rotate it about any axis in the xy-plane, but we'll choose the
    // y-axis for simplicity. This will rotate it clockwise from the
    // perspective of the camera.
    q.setFromAxisAngle(Y_AXIS, 23.4*Math.PI/180);
    mesh.applyQuaternion(q);

    const EARTH_AXIS = Z_AXIS.clone().applyQuaternion(q);

    var t0 = 0;
    function animate(elapsedMilliseconds) {
        var dt = elapsedMilliseconds - t0;
        t0 = elapsedMilliseconds;

        requestAnimationFrame( animate );

        // Rotate the earth around the sun.
        // While the earth actually rotates around the sun, it's easier to
        // model this as the sun rotating around the earth that way we only
        // need to move the position of the sun with each frame instead of the
        // position of the earth and the camera.
        const orbitAngularVelocity = 2*Math.PI/(365*24*HOUR);
        var theta = orbitAngularVelocity * dt % (2*Math.PI);
        q.setFromAxisAngle(Z_AXIS, theta);
        // Note that we're changing the *position* of the light source and not
        // the orientation which would have no effect.
        SUN.position.applyQuaternion(q);


        // Rotate the earth around its axis.
        const rotationAngularVelocity = 2*Math.PI/(24*HOUR); 
        theta = rotationAngularVelocity*dt % (2*Math.PI);
        q.setFromAxisAngle(EARTH_AXIS, theta);
        mesh.applyQuaternion(q);

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
