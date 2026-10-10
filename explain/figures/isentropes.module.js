import * as THREE from '../../js/three.module.js';
import { slider, buttons, legend, readout, text, rampRGB, INK, ACCENT, paletteVersion } from '../runtime.module.js';
import { Scene3D, linear } from '../scene3d.module.js';
import { jw06Point } from '../layerCases.module.js';
import { heightOf } from '../physics.module.js';

const SURFACES = [295, 300, 305, 310, 320], X = 6, Z = 5, LAT_SOUTH = 25, LAT_NORTH = 65, METERS = 1000, NORTHWARD = 8, SPEED = 3600 * 5, TRAIL = 2000, EVERY = 900, HALO = 'rgba(20, 20, 22, 0.85)';

export function mountIsentropes(root) {
  const controls = root.querySelector('.controls');
  let heating = 0, painted = -1, parcel = null;
  const latOf = (z) => (LAT_SOUTH + (Z - z) / (2 * Z) * (LAT_NORTH - LAT_SOUTH)) * Math.PI / 180;
  const zOf = (lat) => Z - (lat * 180 / Math.PI - LAT_SOUTH) / (LAT_NORTH - LAT_SOUTH) * 2 * Z;
  function pressureOn(theta, lat) {
    if (jw06Point(lat, 1000).theta >= theta) return 1000;
    let lo = 150, hi = 1000;
    for (let i = 0; i < 40; i++) { const mid = (lo + hi) / 2; if (jw06Point(lat, mid).theta > theta) lo = mid; else hi = mid; }
    return (lo + hi) / 2;
  }
  const yOf = (p) => heightOf(p * 100) / METERS;

  const scene = new Scene3D(root, { height: 460, distance: 25, target: [0, 3.6, 0], yaw: -0.85, pitch: 0.24, update, draw: labels, step });
  const ground = scene.add(new THREE.Mesh(new THREE.PlaneGeometry(2 * X, 2 * Z).rotateX(-Math.PI / 2), new THREE.MeshLambertMaterial({ color: new THREE.Color(0.05, 0.045, 0.035) })));
  ground.position.y = -0.01;
  const sheets = SURFACES.map((theta) => {
    const geometry = new THREE.PlaneGeometry(2 * X, 2 * Z, 4, 40).rotateX(-Math.PI / 2), pos = geometry.attributes.position;
    for (let v = 0; v < pos.count; v++) pos.setY(v, yOf(pressureOn(theta, latOf(pos.getZ(v)))));
    geometry.computeVertexNormals();
    return scene.add(new THREE.Mesh(geometry, new THREE.MeshLambertMaterial({ color: 0xffffff, transparent: true, opacity: 0.5, side: THREE.DoubleSide, depthWrite: false })));
  });
  const trail = scene.add(new THREE.Line(new THREE.BufferGeometry().setAttribute('position', new THREE.BufferAttribute(new Float32Array(3 * TRAIL), 3)), new THREE.LineBasicMaterial({ color: new THREE.Color(linear(1), linear(0.91), linear(0.63)) })));
  const ball = scene.add(new THREE.Mesh(new THREE.SphereGeometry(0.18, 16, 12), new THREE.MeshLambertMaterial({ color: 0xffffff })));
  legend(root, [['ramp', 'surfaces of equal potential temperature, 295 to 320 K', 'cool', 'warm', 'neutral'], ['line', 'the parcel’s path', ACCENT]]);
  slider(controls, { label: 'Heat the parcel by', min: 0, max: 8, step: 0.5, value: heating, format: (v) => v === 0 ? 'nothing: it keeps its θ' : `${v} K a day`, onInput: (v) => { heating = v; } });
  buttons(controls, [['Launch it again', launch]]);
  const out = readout(controls);

  function launch() { parcel = { lat: 34 * Math.PI / 180, x: -X + 0.5, theta: 300, t: 0, points: 0, next: 0 }; trail.geometry.setDrawRange(0, 0); scene.fig.play(true); }

  function step(dt) {
    if (!parcel) return;
    const seconds = dt * SPEED;
    if (parcel.lat < (LAT_NORTH - 2) * Math.PI / 180 && parcel.x < X - 0.3) {
      const p = pressureOn(parcel.theta, parcel.lat), u = jw06Point(parcel.lat, p).u;
      parcel.lat += NORTHWARD * seconds / 6.371e6;
      parcel.x += u * seconds / 2.6e6;
      parcel.theta += heating * seconds / 86400;
      parcel.t += seconds;
      const arr = trail.geometry.attributes.position.array, n = parcel.points;
      if (parcel.t >= parcel.next && n < TRAIL) { arr[3 * n] = parcel.x; arr[3 * n + 1] = yOf(pressureOn(parcel.theta, parcel.lat)); arr[3 * n + 2] = zOf(parcel.lat); parcel.points++; parcel.next += EVERY; trail.geometry.attributes.position.needsUpdate = true; trail.geometry.setDrawRange(0, parcel.points); }
    } else scene.fig.play(false);
    const p = pressureOn(parcel.theta, parcel.lat);
    ball.position.set(parcel.x, yOf(p), zOf(parcel.lat));
    out.set([['day', (parcel.t / 86400).toFixed(1)], ['latitude', `${(parcel.lat * 180 / Math.PI).toFixed(0)}° north`], ['the parcel’s θ', `${parcel.theta.toFixed(1)} K`], ['its pressure', `${p.toFixed(0)} hPa, about ${(heightOf(p * 100) / 1000).toFixed(1)} km up`]], 4);
  }

  const rgb = [0, 0, 0];
  function update() {
    if (painted === paletteVersion) return;
    painted = paletteVersion;
    SURFACES.forEach((theta, s) => { rampRGB((theta - 295) / 25, rgb); sheets[s].material.color.setRGB(linear(rgb[0]), linear(rgb[1]), linear(rgb[2])); });
  }

  function labels(ctx) {
    const south = scene.project(0, 0.05, Z + 0.6), north = scene.project(0, 0.05, -Z - 0.6);
    if (south.visible) text(ctx, `${LAT_SOUTH}° north`, south.x, south.y, { align: 'center', color: INK, size: 11, halo: HALO });
    if (north.visible) text(ctx, `${LAT_NORTH}° north, toward the pole`, north.x, north.y, { align: 'center', color: INK, size: 11, halo: HALO });
    const middle = 45 * Math.PI / 180;
    SURFACES.forEach((theta) => { const at = scene.project(X, yOf(pressureOn(theta, middle)), zOf(middle)); if (at.visible) text(ctx, `${theta} K`, at.x + 6, at.y, { color: INK, size: 10, halo: HALO }); });
  }

  launch();
}
