import * as THREE from '../../js/three.module.js';
import { slider, choice, legend, readout, text, MUTED, rampRGB, paletteVersion } from '../runtime.module.js';
import { Scene3D, Arrows, linear } from '../scene3d.module.js';
import { R, G } from '../physics.module.js';

const LEVELS = [850, 700, 500, 300, 200], SIZE = 12, SEGMENTS = 40, METERS = 1800, EXAG = 15, WIDTH = 2.6;

export function mountThickness(root) {
  const controls = root.querySelector('.controls');
  let warming = 6, depth = 'column', painted = -1;
  const temperature = (p) => 288.15 * (p / 1013.25) ** (R * 0.0065 / G);
  const anomaly = (x, z) => warming * Math.exp(-(x * x + z * z) / (2 * WIDTH * WIDTH));
  function height(p, x, z) {
    const steps = 24, a = anomaly(x, z), top = depth === 'column' ? 0 : 700;
    let total = 0;
    for (let s = 0; s < steps; s++) { const p0 = 1000 - (1000 - p) * s / steps, p1 = 1000 - (1000 - p) * (s + 1) / steps, mid = (p0 + p1) / 2; total += R * (temperature(mid) + (mid >= top ? a : 0)) / G * Math.log(p0 / p1); }
    return total;
  }
  const scene = new Scene3D(root, { height: 460, distance: 30, target: [0, 3.6, 0], yaw: -0.5, pitch: 0.32, update, draw: labels });
  const ground = scene.add(new THREE.Mesh(new THREE.PlaneGeometry(SIZE, SIZE).rotateX(-Math.PI / 2), new THREE.MeshLambertMaterial({ color: new THREE.Color(0.05, 0.045, 0.035) })));
  ground.position.y = -0.01;
  const sheets = LEVELS.map(() => {
    const geometry = new THREE.PlaneGeometry(SIZE, SIZE, SEGMENTS, SEGMENTS).rotateX(-Math.PI / 2);
    geometry.setAttribute('color', new THREE.BufferAttribute(new Float32Array(geometry.attributes.position.count * 3), 3));
    return scene.add(new THREE.Mesh(geometry, new THREE.MeshLambertMaterial({ vertexColors: true, transparent: true, opacity: 0.7, side: THREE.DoubleSide, depthWrite: false })));
  });
  const forces = new Arrows(scene.group, 64, { radius: 0.04, head: 0.22, headRadius: 0.11, dashed: true });
  legend(root, [['ramp', 'how much each surface of equal pressure has risen, drawn 15 times larger', 'neutral', 'warm'], ['force', 'the pressure gradient force on the top surface', 'warm']]);
  slider(controls, { label: 'Warm the middle by', min: 0, max: 10, step: 0.5, value: warming, format: (v) => `${v} °C`, onInput: (v) => { warming = v; painted = -1; scene.fig.render(); } });
  choice(controls, { label: 'Warm', options: [['the whole column', 'column'], ['only below 700 hPa', 'low']], value: depth, onChange: (v) => { depth = v; painted = -1; scene.fig.render(); }, span: true });
  const out = readout(controls);

  const rgb = [0, 0, 0];
  function update() {
    if (painted === paletteVersion) return;
    painted = paletteVersion;
    let topRise = 0;
    LEVELS.forEach((p, s) => {
      const geometry = sheets[s].geometry, pos = geometry.attributes.position, col = geometry.attributes.color;
      const base = height(p, 1e9, 1e9);
      for (let v = 0; v < pos.count; v++) {
        const x = pos.getX(v), z = pos.getZ(v), rise = height(p, x, z) - base;
        pos.setY(v, base / METERS + EXAG * rise / METERS);
        rampRGB(0.5 + Math.min(1, rise / 150) / 2, rgb);
        col.setXYZ(v, linear(rgb[0]), linear(rgb[1]), linear(rgb[2]));
        if (s === LEVELS.length - 1) topRise = Math.max(topRise, rise);
      }
      pos.needsUpdate = true; col.needsUpdate = true; geometry.computeVertexNormals();
    });
    forces.begin();
    const warm = rampRGB(1, [0, 0, 0]), top = LEVELS[LEVELS.length - 1], base = height(top, 1e9, 1e9);
    for (let ring = 1; ring <= 2; ring++) for (let k = 0; k < 8; k++) {
      const a = k * Math.PI / 4 + ring * 0.3, r = ring * 1.6, x = r * Math.cos(a), z = r * Math.sin(a), e = 0.05;
      const gx = -G * (height(top, x + e, z) - height(top, x - e, z)) / (2 * e), gz = -G * (height(top, x, z + e) - height(top, x, z - e)) / (2 * e), scale = 0.012;
      const y = base / METERS + EXAG * (height(top, x, z) - base) / METERS + 0.15;
      forces.push(x, y, z, gx * scale, 0, gz * scale, warm);
    }
    forces.end();
    const lowRise = height(850, 0, 0) - height(850, 1e9, 1e9);
    out.set([['the 200 hPa surface rises', `${topRise.toFixed(0)} m over the warm spot`], ['the 850 hPa surface rises', `${lowRise.toFixed(0)} m`], ['the ground', 'does not move: no air has gone anywhere yet']]);
  }

  function labels(ctx) {
    LEVELS.forEach((p, s) => { const y = height(p, 1e9, 1e9) / METERS, at = scene.project(-SIZE / 2, y, -SIZE / 2); if (at.visible) text(ctx, `${p} hPa`, at.x - 6, at.y, { align: 'right', color: MUTED, size: 10 }); });
  }
}
