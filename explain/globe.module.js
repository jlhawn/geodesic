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
  constructor(root, { height = 440, zoom = 1.08, lookAt = null, center = null, lighting = { directional: 2.8, ambient: 0.7 } } = {}) {
    root.classList.add('drag');
    this.fig = new Figure(root, { height, minHeight: 300, context: 'webgl', draw: () => this.render() });
    this.renderer = new THREE.WebGLRenderer({ canvas: this.fig.canvas, antialias: true, alpha: true });
    this.scene = new THREE.Scene();
    this.camera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.1, 10);
    this.camera.position.set(0, 0, 4);
    this.camera.lookAt(0, 0, 0);
    this.group = new THREE.Group();
    this.scene.add(this.group);
    const light = new THREE.DirectionalLight(0xffffff, lighting.directional);
    light.position.set(-1, 1.2, 2.5);
    this.scene.add(light, new THREE.AmbientLight(0xffffff, lighting.ambient));
    this.zoom = zoom;
    this.base = lookAt ? new THREE.Quaternion().setFromUnitVectors(lookAt.clone().normalize(), new THREE.Vector3(0, 0, 1)) : NORTH_UP.clone();
    if (center) {
      const spin = new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0, 1, 0), -Math.PI / 2 - center.lon);
      const tilt = new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1, 0, 0), center.lat);
      this.base = tilt.multiply(spin).multiply(NORTH_UP.clone());
    }
    this.spin = new THREE.Quaternion();
    this.layers = [];
    this.fig.render();
    let grabbed = null;
    this.fig.pointer({
      down: (p) => { grabbed = this.surfacePoint(p).applyQuaternion(this.spin.clone().invert()); },
      move: (p) => {
        if (!grabbed) return;
        this.spin.premultiply(new THREE.Quaternion().setFromUnitVectors(grabbed.clone().applyQuaternion(this.spin), this.surfacePoint(p)));
        this.fig.render();
      },
      up: () => { grabbed = null; },
    });
  }

  surfacePoint({ x, y }) {
    const { width: w, h } = this.fig, aspect = w / h;
    const sx = (x / w - 0.5) * 2 * this.zoom * aspect, sy = (0.5 - y / h) * 2 * this.zoom, r2 = sx * sx + sy * sy;
    return r2 < 1 ? new THREE.Vector3(sx, sy, Math.sqrt(1 - r2)) : new THREE.Vector3(sx, sy, 0).normalize();
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

  dynamicCells(polygons) {
    const positions = [], spans = [];
    for (const cell of polygons) {
      const n = cell.vertices.length, first = positions.length / 3;
      for (let k = 0; k < n; k++) {
        const a = cell.vertices[k], c = cell.vertices[(k + 1) % n];
        positions.push(cell.center.x, cell.center.y, cell.center.z, a.x, a.y, a.z, c.x, c.y, c.z);
      }
      spans.push(first, 3 * n);
    }
    const geometry = new THREE.BufferGeometry(), colors = new Float32Array(positions.length);
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('normal', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    const mesh = new THREE.Mesh(geometry, new THREE.MeshLambertMaterial({ vertexColors: true }));
    this.group.add(mesh);
    this.layers.push(mesh);
    const rgb = [0, 0, 0];
    return (colorOf) => {
      for (let i = 0; i < polygons.length; i++) {
        colorOf(i, rgb);
        for (let v = spans[2 * i], end = v + spans[2 * i + 1]; v < end; v++) { colors[3 * v] = rgb[0]; colors[3 * v + 1] = rgb[1]; colors[3 * v + 2] = rgb[2]; }
      }
      geometry.attributes.color.needsUpdate = true;
    };
  }

  dynamicArrows(count, { color = 0xffffff, opacity = 0.9, lift = 1.006 } = {}) {
    const positions = new Float32Array(count * 6 * 3), geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    const lines = new THREE.LineSegments(geometry, new THREE.LineBasicMaterial({ color, transparent: opacity < 1, opacity }));
    this.group.add(lines);
    this.layers.push(lines);
    return (centers, vectors, scale) => {
      for (let k = 0; k < count; k++) {
        const cx = centers[3 * k], cy = centers[3 * k + 1], cz = centers[3 * k + 2];
        const vx = vectors[3 * k] * scale, vy = vectors[3 * k + 1] * scale, vz = vectors[3 * k + 2] * scale, len = Math.hypot(vx, vy, vz);
        const tx = cx + vx, ty = cy + vy, tz = cz + vz;
        let sx = cy * vz - cz * vy, sy = cz * vx - cx * vz, sz = cx * vy - cy * vx;
        const sl = Math.hypot(sx, sy, sz) || 1, head = Math.min(0.35 * len, 0.02);
        sx /= sl; sy /= sl; sz /= sl;
        const bx = len ? vx / len : 0, by = len ? vy / len : 0, bz = len ? vz / len : 0;
        const pts = [cx - vx / 2, cy - vy / 2, cz - vz / 2, tx - vx / 2, ty - vy / 2, tz - vz / 2];
        const hx = tx - vx / 2, hy = ty - vy / 2, hz = tz - vz / 2;
        pts.push(hx, hy, hz, hx - head * (bx + 0.6 * sx), hy - head * (by + 0.6 * sy), hz - head * (bz + 0.6 * sz));
        pts.push(hx, hy, hz, hx - head * (bx - 0.6 * sx), hy - head * (by - 0.6 * sy), hz - head * (bz - 0.6 * sz));
        for (let j = 0; j < 18; j++) positions[18 * k + j] = pts[j] * lift;
      }
      geometry.attributes.position.needsUpdate = true;
    };
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
