import * as THREE from "./three.module.js";
import { SetupCameraAndRenderer } from "./sphericalCamera.module.js";
import { Stats } from "./stats.module.js";

export function initSphericalViewer(container, grid, config = {}) {
  const numPentagons = 12;
  const numHexagons = grid.size - numPentagons;
  const numTriangles = 3*numPentagons + 4*numHexagons;

  const vertices = new Float32Array(9*numTriangles);
  const normals = new Float32Array(9*numTriangles);
  const colors = new Float32Array(12*numTriangles);

  const color = new THREE.Color();
  const n = new THREE.Vector3();

  let tri = 0; // Triangle counter.

  for (const gridCell of grid) {
    const {x, y, z} = gridCell.centerVertex;
    let {r, g, b} = config.getColor(gridCell);
    color.setRGB(r, g, b);

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

  function disposeArray() { this.array = null; };

  const geometry = new THREE.BufferGeometry();

  geometry.setAttribute('position', new THREE.BufferAttribute(vertices, 3).onUpload(disposeArray));
  geometry.setAttribute('normal', new THREE.BufferAttribute(normals, 3).onUpload(disposeArray));
  geometry.setAttribute('color', new THREE.BufferAttribute(colors, 4).onUpload(disposeArray));

  geometry.computeBoundingSphere();

  const material = new THREE.MeshBasicMaterial({vertexColors: true});

  const mesh = new THREE.Mesh(geometry, material);
  
  const { scene, render } = SetupCameraAndRenderer(container);
  scene.background = new THREE.Color(config.backgroundColor);

  scene.add(mesh);

  const stats = new Stats();
  document.body.appendChild( stats.dom );

  function animate(t) {
    render();
    stats.update();

    requestAnimationFrame( animate );
  }
  requestAnimationFrame( animate );
}
