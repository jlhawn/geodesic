import * as THREE from '../js/three.module.js';
import { Figure } from './runtime.module.js';

const NORTH_UP = new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 0, 1), new THREE.Vector3(0, 1, 0));

export const rgb = (css) => css.match(/\d+(\.\d+)?/g).slice(0, 3).map((v) => Number(v) / 255);

export function arc(a, b, steps = 12) {
  const points = [];
  for (let s = 0; s <= steps; s++) points.push(a.clone().lerp(b, s / steps).normalize());
  return points;
}

export class Globe {
  constructor(root, { height = 440, zoom = 1.08, lookAt = null } = {}) {
    root.classList.add('drag');
    this.fig = new Figure(root, { height, minHeight: 300, context: 'webgl', draw: () => this.render() });
    this.renderer = new THREE.WebGLRenderer({ canvas: this.fig.canvas, antialias: true, alpha: true });
    this.scene = new THREE.Scene();
    this.camera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.1, 10);
    this.camera.position.set(0, 0, 4);
    this.camera.lookAt(0, 0, 0);
    this.group = new THREE.Group();
    this.scene.add(this.group);
    const light = new THREE.DirectionalLight(0xffffff, 2.8);
    light.position.set(-1, 1.2, 2.5);
    this.scene.add(light, new THREE.AmbientLight(0xffffff, 0.7));
    this.zoom = zoom;
    this.base = lookAt ? new THREE.Quaternion().setFromUnitVectors(lookAt.clone().normalize(), new THREE.Vector3(0, 0, 1)) : NORTH_UP.clone();
    this.spin = new THREE.Quaternion();
    this.layers = [];
    let last = null;
    this.fig.render();
    this.fig.pointer({
      down: (p) => { last = p; },
      move: (p) => {
        const dx = p.x - last.x, dy = p.y - last.y;
        last = p;
        this.spin.premultiply(new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0, 1, 0), dx * 0.008)).premultiply(new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1, 0, 0), dy * 0.008));
        this.fig.render();
      },
    });
  }

  clear() {
    for (const layer of this.layers) { this.group.remove(layer); layer.geometry.dispose(); layer.material.dispose(); }
    this.layers = [];
  }

  cells(polygons, colorOf) {
    const positions = [], colors = [];
    polygons.forEach((cell, i) => {
      const [r, g, b] = colorOf(cell, i), n = cell.vertices.length;
      for (let k = 0; k < n; k++) {
        const a = cell.vertices[k], c = cell.vertices[(k + 1) % n];
        positions.push(cell.center.x, cell.center.y, cell.center.z, a.x, a.y, a.z, c.x, c.y, c.z);
        colors.push(r, g, b, r, g, b, r, g, b);
      }
    });
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('normal', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    const mesh = new THREE.Mesh(geometry, new THREE.MeshLambertMaterial({ vertexColors: true }));
    this.group.add(mesh);
    this.layers.push(mesh);
    return mesh;
  }

  lines(paths, { color = 0xffffff, opacity = 0.35, lift = 1.002 } = {}) {
    const positions = [];
    for (const path of paths) for (let k = 0; k + 1 < path.length; k++) { const a = path[k], b = path[k + 1]; positions.push(a.x * lift, a.y * lift, a.z * lift, b.x * lift, b.y * lift, b.z * lift); }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    const lines = new THREE.LineSegments(geometry, new THREE.LineBasicMaterial({ color, transparent: opacity < 1, opacity }));
    this.group.add(lines);
    this.layers.push(lines);
    return lines;
  }

  coloredLines(segments, { lift = 1.003 } = {}) {
    const positions = [], colors = [];
    for (const { a, b, color } of segments) { positions.push(a.x * lift, a.y * lift, a.z * lift, b.x * lift, b.y * lift, b.z * lift); colors.push(...color, ...color); }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    const lines = new THREE.LineSegments(geometry, new THREE.LineBasicMaterial({ vertexColors: true }));
    this.group.add(lines);
    this.layers.push(lines);
    return lines;
  }

  render() {
    if (!this.renderer) return;
    const { width: w, h, dpr } = this.fig;
    if (!w) return;
    this.renderer.setPixelRatio(dpr);
    this.renderer.setSize(w, h, false);
    const aspect = w / h;
    Object.assign(this.camera, { left: -this.zoom * aspect, right: this.zoom * aspect, top: this.zoom, bottom: -this.zoom });
    this.camera.updateProjectionMatrix();
    this.group.quaternion.copy(this.spin).multiply(this.base);
    this.renderer.render(this.scene, this.camera);
  }
}
