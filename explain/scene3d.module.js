import * as THREE from '../js/three.module.js';
import { Figure, clamp } from './runtime.module.js';

export const linear = (c) => (c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);
export const cssColor = (css) => { const m = css[0] === '#' ? [1, 3, 5].map((i) => parseInt(css.slice(i, i + 2), 16)) : css.match(/\d+(\.\d+)?/g).map(Number); return [m[0] / 255, m[1] / 255, m[2] / 255]; };

export class Scene3D {
  constructor(root, { height = 440, minHeight = 300, distance = 16, target = [0, 2, 0], yaw = -0.6, pitch = 0.45, fov = 32, step = null, draw = null, update = null } = {}) {
    root.classList.add('drag');
    this.drawOverlay = draw;
    this.update = update;
    this.fig = new Figure(root, { height, minHeight, context: 'webgl', step, draw: () => this.render() });
    this.overlay = document.createElement('canvas');
    Object.assign(this.overlay.style, { position: 'absolute', left: '0', top: '0', width: '100%', pointerEvents: 'none' });
    this.fig.stage.style.position = 'relative';
    this.fig.stage.append(this.overlay);
    this.renderer = new THREE.WebGLRenderer({ canvas: this.fig.canvas, antialias: true, alpha: true });
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(fov, 1, 0.1, 400);
    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x404050, 2.2));
    const sun = new THREE.DirectionalLight(0xffffff, 1.6);
    sun.position.set(-4, 10, 6);
    this.scene.add(sun);
    this.group = new THREE.Group();
    this.scene.add(this.group);
    Object.assign(this, { distance, target: new THREE.Vector3(...target), yaw, pitch });
    let last = null;
    this.fig.pointer({
      down: (p) => { last = p; },
      move: (p) => { this.yaw -= (p.x - last.x) * 0.008; this.pitch = clamp(this.pitch + (p.y - last.y) * 0.006, 0.05, 1.45); last = p; this.fig.render(); },
    });
    queueMicrotask(() => { this.started = true; this.fig.render(); });
  }

  render() {
    const { width: w, h, dpr } = this.fig;
    if (!this.started || !w) return;
    this.update?.();
    const size = `${w}x${h}@${dpr}`;
    if (size !== this.size) {
      this.size = size;
      this.renderer.setPixelRatio(dpr);
      this.renderer.setSize(w, h, false);
      this.camera.aspect = w / h;
      this.camera.updateProjectionMatrix();
      this.overlay.width = Math.round(w * dpr);
      this.overlay.height = Math.round(h * dpr);
      this.overlay.style.height = `${h}px`;
    }
    const c = Math.cos(this.pitch), d = this.distance * Math.max(1, 1.45 * h / w);
    this.camera.position.set(this.target.x + d * c * Math.sin(this.yaw), this.target.y + d * Math.sin(this.pitch), this.target.z + d * c * Math.cos(this.yaw));
    this.camera.lookAt(this.target);
    this.renderer.render(this.scene, this.camera);
    const ctx = this.overlay.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    this.drawOverlay?.(ctx, w, h);
  }

  project(x, y, z) {
    const v = new THREE.Vector3(x, y, z).project(this.camera);
    return { x: (v.x + 1) / 2 * this.fig.width, y: (1 - v.y) / 2 * this.fig.h, visible: v.z < 1 };
  }

  add(object) { this.group.add(object); return object; }
}

export function hexPatch(cols, rows, radius = 1) {
  const dx = Math.sqrt(3) * radius, dz = 1.5 * radius;
  const offset = [((cols - 1) + 0.5) * dx / 2, (rows - 1) * dz / 2];
  const index = (c, r) => r * cols + c;
  const cells = [];
  for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) {
    const x = c * dx + (r & 1) * dx / 2 - offset[0], z = r * dz - offset[1];
    const corners = Array.from({ length: 6 }, (_, k) => { const a = Math.PI / 6 + k * Math.PI / 3; return [x + radius * Math.cos(a), z + radius * Math.sin(a)]; });
    cells.push({ c, r, x, z, corners });
  }
  const around = (c, r) => ((r & 1) ? [[1, 0], [0, 1], [1, 1], [-1, 0], [0, -1], [1, -1]] : [[1, 0], [-1, 1], [0, 1], [-1, 0], [-1, -1], [0, -1]]).map(([a, b]) => [c + a, r + b]);
  const edges = [];
  for (const cell of cells) for (const [c, r] of around(cell.c, cell.r)) {
    if (c < 0 || c >= cols || r < 0 || r >= rows) continue;
    const j = index(c, r), i = index(cell.c, cell.r);
    if (j < i) continue;
    const other = cells[j], len = Math.hypot(other.x - cell.x, other.z - cell.z);
    edges.push({ i, j, x: (cell.x + other.x) / 2, z: (cell.z + other.z) / 2, nx: (other.x - cell.x) / len, nz: (other.z - cell.z) / len, d: len });
  }
  for (const e of edges) { (cells[e.i].edges ??= []).push({ edge: e, sign: 1 }); (cells[e.j].edges ??= []).push({ edge: e, sign: -1 }); }
  return { cells, edges, index, dx, dz };
}

const UP = new THREE.Vector3(0, 1, 0);

export class Arrows {
  constructor(scene, capacity, { radius = 0.035, head = 0.16, headRadius = 0.09, dashed = false } = {}) {
    const material = new THREE.MeshLambertMaterial({ color: 0xffffff });
    this.dashed = dashed;
    this.segments = dashed ? 3 : 1;
    this.shafts = new THREE.InstancedMesh(new THREE.CylinderGeometry(radius, radius, 1, 8, 1, false).translate(0, 0.5, 0), material, capacity * this.segments);
    this.heads = new THREE.InstancedMesh(dashed ? new THREE.ConeGeometry(headRadius, head, 12, 1, true).translate(0, -head / 2, 0) : new THREE.ConeGeometry(headRadius, head, 12).translate(0, -head / 2, 0), dashed ? new THREE.MeshLambertMaterial({ color: 0xffffff, side: THREE.DoubleSide }) : material, capacity);
    for (const mesh of [this.shafts, this.heads]) { mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage); mesh.frustumCulled = false; scene.add(mesh); }
    Object.assign(this, { capacity, head, n: 0 });
    this.m = new THREE.Matrix4(); this.q = new THREE.Quaternion(); this.v = new THREE.Vector3(); this.s = new THREE.Vector3(); this.p = new THREE.Vector3(); this.color = new THREE.Color();
  }

  begin() { this.n = 0; }

  push(x, y, z, vx, vy, vz, rgb) {
    const len = Math.hypot(vx, vy, vz);
    if (len < 1e-4 || this.n >= this.capacity) return;
    this.v.set(vx / len, vy / len, vz / len);
    this.q.setFromUnitVectors(UP, this.v);
    this.color.setRGB(linear(rgb[0]), linear(rgb[1]), linear(rgb[2]));
    const shaft = Math.max(0, len - this.head);
    for (let k = 0; k < this.segments; k++) {
      const piece = this.dashed ? shaft / (2 * this.segments - 1) : shaft, start = this.dashed ? 2 * k * piece : 0;
      this.p.set(x + this.v.x * start, y + this.v.y * start, z + this.v.z * start);
      this.m.compose(this.p, this.q, this.s.set(1, Math.max(piece, 1e-4), 1));
      this.shafts.setMatrixAt(this.n * this.segments + k, this.m);
      this.shafts.setColorAt(this.n * this.segments + k, this.color);
    }
    this.p.set(x + vx, y + vy, z + vz);
    this.m.compose(this.p, this.q, this.s.set(1, 1, 1));
    this.heads.setMatrixAt(this.n, this.m);
    this.heads.setColorAt(this.n, this.color);
    this.n++;
  }

  end() {
    this.shafts.count = this.n * this.segments; this.heads.count = this.n;
    for (const mesh of [this.shafts, this.heads]) { mesh.instanceMatrix.needsUpdate = true; if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true; }
  }
}

export function hexTiles(scene, { thickness = 0.06, opacity = 0.85 } = {}) {
  const shape = new THREE.Shape();
  for (let k = 0; k < 6; k++) { const a = Math.PI / 6 + k * Math.PI / 3; const x = 0.96 * Math.cos(a), z = 0.96 * Math.sin(a); if (k) shape.lineTo(x, z); else shape.moveTo(x, z); }
  const geometry = new THREE.ExtrudeGeometry(shape, { depth: thickness, bevelEnabled: false }).rotateX(Math.PI / 2);
  const material = new THREE.MeshLambertMaterial({ color: 0xffffff, transparent: opacity < 1, opacity });
  return { geometry, material, mesh: (count) => { const mesh = new THREE.InstancedMesh(geometry, material, count); mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage); mesh.frustumCulled = false; scene.add(mesh); return mesh; } };
}
