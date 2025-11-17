import SimplexNoise from "./simplex-noise.module.js";
import * as THREE from "./three.module.js";
import { Grid, circumcenter } from "./grid.module.js";
import { SetupCameraControls } from "./camera.module.js";

const DRAW_CENTER_VERTICES = false;
const DRAW_AREA_DIFF = false;
const DRAW_CELL_EDGES = false;

const HIGHLIGHT_ICOSAHEDRON_EDGE = false;

const AMBIENT_LIGHT = 0x555555;

// Length of simulated second in real milliseconds.
const SECOND = 0.5;
const MINUTE = 60*SECOND;
const HOUR = 60*MINUTE;
const DAY = 24*HOUR;
const YEAR = 365*DAY;
const AXIAL_TILT = 23.4*Math.PI/180;

// The equation for the number of cells is 10*N^2+2
//  10 ->     1,002
//  20 ->     4,002
//  50 ->    25,002
// 100 ->   100,002
// 200 ->   400,002
// 320 -> 1,024,002
// 400 -> 1,600,002
// 500 -> 2,500,002
const SUBDIVISIONS = 100;

if (SUBDIVISIONS < 1) {
    throw "SUBDIVISIONS must be greater than zero.";
}

const NUM_CELLS = 10*SUBDIVISIONS*SUBDIVISIONS + 2;
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

    // --- NEW: precompute per-cell average neighbor angles and min/max ---
    let minArea = Infinity;
    let maxArea = -Infinity;
    let avgArea = 0;
    const areas = [];
    for (const cell of grid) {
        if (cell.isPentagon) continue;
        const a = cell.area;
        areas.push(a);
        if (a < minArea) minArea = a;
        if (a > maxArea) maxArea = a;
        avgArea += a;
    }
    avgArea /= NUM_CELLS;

    // Calculate percentiles
    areas.sort((a, b) => a - b);
    const p1Index = Math.floor(NUM_CELLS * 0.01);
    const p50Index = Math.floor(NUM_CELLS * 0.50);
    const p99Index = Math.floor(NUM_CELLS * 0.99);
    const p1Area = areas[p1Index];
    const p50Area = areas[p50Index];
    const p99Area = areas[p99Index];

    const areaRange = maxArea - minArea;
    const areaRangeOverAverage = areaRange/avgArea;
    const areaRangeOverp50 = areaRange/p50Area;
    const invRange = 1 / Math.max(1e-12, areaRange); // guard against division by zero
    console.log({avgArea, minArea, maxArea, areaRange, areaRangeOverAverage, areaRangeOverp50, p1Area, p50Area, p99Area});

    // console.log(`completed grid construction in ${performance.now()-start}ms`);
    start = performance.now();

    const numPentagons = 12;
    const numHexagons = NUM_CELLS - numPentagons;
    const numTriangles = 3*numPentagons + 4*numHexagons;

    const vertices = new Float32Array(9*numTriangles);
    const normals = new Float32Array(9*numTriangles);
    const colors = new Float32Array(12*numTriangles);

    const color = new THREE.Color();
    const n = new THREE.Vector3();

    let tri = 0; // Triangle counter.

    for (const gridCell of grid) {
        if (DRAW_AREA_DIFF) {
            let r, g, b;
            if (gridCell.isPentagon) {
                r = 0; g = 1; b = 0;
            } else {
                // Map angle to [0,1]
                const t = (gridCell.area - minArea) * invRange;

                // Red (low) → Blue (high)
                r = 1 - t;
                g = 0;
                b = t;

                // if (gridCell.isAlongIcosahedronEdge) {
                //     r = 0.5, g = 0.5, b = 0.5;
                // }

                // if (gridCell.area > p99Area) {
                //     console.log(`grid cell ${gridCell.id} is in the top 1% of largest cells with an area of ${gridCell.area}`);
                //     r = 0.9, g = 0.9, b = 0.9;
                // }
            }

            color.setRGB(r, g, b);
        } else {
            const {x, y, z} = gridCell.centerVertex;
            let {r, g, b} = simplexColor(x, y, z);
            if (HIGHLIGHT_ICOSAHEDRON_EDGE && gridCell.isAlongIcosahedronEdge) {
                r = lerp(r, 0, 1, 0.1, 1);
                g = lerp(g, 0, 1, 0.1, 1);
                b = lerp(b, 0, 1, 0.1, 1);
                // [r, g, b] = [0.9, 0.9, 0.9];
            }
            color.setRGB(r, g, b);
        }

        gridCell.faceTriangles.forEach(function(triangle) {
            const { a, b, c } = triangle;

            // Set vertex positions.
            vertices.set([a.x, a.y, a.z], tri*9);
            vertices.set([b.x, b.y, b.z], tri*9+3);
            vertices.set([c.x, c.y, c.z], tri*9+6);

            // Get normal vector.
            triangle.getNormal(n);

            // One for each vertex.
            normals.set([n.x, n.y, n.z], tri*9);
            normals.set([n.x, n.y, n.z], tri*9+3);
            normals.set([n.x, n.y, n.z], tri*9+6);

            // Set color for each vertex to be the same.
            colors.set([color.r, color.g, color.b, 1], tri*12);
            colors.set([color.r, color.g, color.b, 1], tri*12+4);
            colors.set([color.r, color.g, color.b, 1], tri*12+8);

            tri++; // IMPORTANT! Increment triangle counter.
        });
    }

    console.log(`completed vertex positions/normals/colors in ${performance.now()-start}ms`);

    function disposeArray() { this.array = null; };

    const geometry = new THREE.BufferGeometry();

    geometry.setAttribute('position', new THREE.BufferAttribute(vertices, 3).onUpload(disposeArray));
    geometry.setAttribute('normal', new THREE.BufferAttribute(normals, 3).onUpload(disposeArray));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 4).onUpload(disposeArray));

    geometry.computeBoundingSphere();

    const material = new THREE.MeshLambertMaterial({
        vertexColors: true,
        transparent: true,
        opacity: 1.0
    });

    const mesh = new THREE.Mesh(geometry, material);

    const scene = new THREE.Scene();
    scene.add( new THREE.AmbientLight(AMBIENT_LIGHT) );
    
    const { camera, renderer, stats } = SetupCameraAndRenderer();

    scene.add(mesh);

    // Create lines for cell edges
    if (DRAW_CELL_EDGES) {
        const edgePositions = [];
        for (const gridCell of grid) {
            const circumcenters = gridCell.vertices;

            // Scale slightly above surface
            const hoverScale = 1.0001;

            // Create line segments connecting adjacent vertices
            for (let i = 0; i < circumcenters.length; i++) {
                const v1 = circumcenters[i];
                const v2 = circumcenters[(i + 1) % circumcenters.length];

                edgePositions.push(
                    v1.x * hoverScale, v1.y * hoverScale, v1.z * hoverScale,
                    v2.x * hoverScale, v2.y * hoverScale, v2.z * hoverScale
                );
            }
        }

        const edgesGeometry = new THREE.BufferGeometry();
        edgesGeometry.setAttribute('position', new THREE.Float32BufferAttribute(edgePositions, 3));

        const edgesMaterial = new THREE.LineBasicMaterial({
            color: 0xffffff,
            transparent: true,
            opacity: 0.05
        });

        const edges = new THREE.LineSegments(edgesGeometry, edgesMaterial);
        scene.add(edges);
    }

    // Create points for cell centers (scaled to float above surface)
    if (DRAW_CENTER_VERTICES) {
        const pointPositions = new Float32Array(NUM_CELLS * 3);
        let pointIndex = 0;
        for (const gridCell of grid) {
            const {x, y, z} = gridCell.centerVertex;
            const hoverScale = 1.0001;
            pointPositions[pointIndex++] = x * hoverScale;
            pointPositions[pointIndex++] = y * hoverScale;
            pointPositions[pointIndex++] = z * hoverScale;
        }

        const pointsGeometry = new THREE.BufferGeometry();
        pointsGeometry.setAttribute('position', new THREE.BufferAttribute(pointPositions, 3).onUpload(disposeArray));

        const pointsMaterial = new THREE.PointsMaterial({
            color: 0xffffff,
            size: 0.005,
            transparent: true,
            opacity: 0.5,
        });

        const points = new THREE.Points(pointsGeometry, pointsMaterial);
        scene.add(points);
    }

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

export default runApp;
