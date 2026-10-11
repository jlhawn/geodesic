import * as THREE from '../../js/three.module.js';
import { slider, buttons, legend, caption, readout, text, termColor, anomalyColor, ACCENT, INK, MUTED, LINE } from '../runtime.module.js';
import { G } from '../physics.module.js';
import { Scene3D, Arrows, linear, cssColor } from '../scene3d.module.js';
import { OMEGA, RADIUS, launch, step as advance, parts, latitude } from '../coriolisSphere.module.js';

const D = Math.PI / 180, SPEED = 5400, H = 60, TRAIL = 3 * 86400, CAPACITY = Math.ceil(TRAIL / H) + 64, LIFT = 1.004, WIDTH = 0.0075, OFFSET = 55 * D, MARK = -90 * D, FAR = 0.4, SOUTH = -Math.sin(0.1 * D);
const AXIS = [-1.3, 1.36], HALO = 'rgba(20, 20, 22, 0.85)';
const COMPASS = ['north', 'north-northeast', 'northeast', 'east-northeast', 'east', 'east-southeast', 'southeast', 'south-southeast', 'south', 'south-southwest', 'southwest', 'west-southwest', 'west', 'west-northwest', 'northwest', 'north-northwest'];
const PALE = cssColor(INK);
const place = (x) => [x[0] / RADIUS, x[2] / RADIUS, -x[1] / RADIUS];
const turn = (v) => [v[0], v[2], -v[1]];
const tilt = (lat) => 0.65 + 0.25 * Math.max(0, Math.abs(lat) - 60) / 30;
const latitudeName = (v) => v === 0 ? 'equator' : `${Math.abs(v)}° ${v > 0 ? 'north' : 'south'}`;
const rounded = (x) => { const p = 10 ** (Math.floor(Math.log10(x)) - 1); return (Math.round(x / p) * p).toLocaleString('en-US'); };
const layout = (w) => { const r = Math.round(Math.max(60, Math.min(95, 0.21 * w))), globe = Math.round(Math.min(450, Math.max(320, 0.7 * w))); return { r, globe, height: globe + 2 * r + 62 }; };
const enlarge = (w) => Math.max(1, Math.min(1.25, 1 + 0.25 * (640 - w) / 300));

export function mountCoriolis(root) {
  const controls = root.querySelector('.controls'), SURFACE = termColor('c'), VERTICAL = termColor('d');
  let lat = 45, heading = 90, speed = 100, state = null, trail = [], pending = 0;

  const scene = new Scene3D(root, { height: (w) => layout(w).height, distance: 4.8, target: [0, 0.1, 0], yaw: 0, pitch: tilt(45), fit: 1, pitchRange: [-1.45, 1.45], view: (w) => [0, layout(w).globe], step, update, draw: labels });
  caption(root, 'The globe is drawn turning with the planet; beneath it, the same parcel seen from above the north pole. Time runs 1.5 hours per second. The force arrows have a scale of their own, and where the force and one of its parts coincide, the force is drawn over the part.');
  legend(root, [['arrow', 'velocity', INK], ['force', 'Coriolis force, −2Ω × v', 'cool'], ['force', 'its part along the ground', SURFACE], ['force', 'its part along the vertical', VERTICAL], ['faint', 'how the force splits into the two', 'rgba(255,255,255,0.7)'], ['faint', 'the local vertical', 'rgba(255,255,255,0.4)'], ['faint', 'the plane through the parcel at right angles to the axis', 'cool'], ['line', 'the path over the last three days', ACCENT]]);

  const globe = new THREE.SphereGeometry(1, 64, 48), shell = { color: new THREE.Color(linear(0.2), linear(0.24), linear(0.32)), transparent: true, opacity: 0.25, depthWrite: false };
  scene.add(new THREE.Mesh(globe, new THREE.MeshLambertMaterial({ ...shell, side: THREE.BackSide }))).renderOrder = 0;
  scene.add(new THREE.Mesh(globe, new THREE.MeshLambertMaterial({ ...shell, side: THREE.FrontSide }))).renderOrder = 2;
  const circle = (points) => new THREE.BufferGeometry().setFromPoints(points.map((p) => new THREE.Vector3(...p)));
  const parallel = (phi) => Array.from({ length: 129 }, (_, k) => { const a = 2 * Math.PI * k / 128; return [1.001 * Math.cos(phi) * Math.cos(a), 1.001 * Math.sin(phi), 1.001 * Math.cos(phi) * Math.sin(a)]; });
  const meridian = (a) => Array.from({ length: 65 }, (_, k) => { const phi = -Math.PI / 2 + Math.PI * k / 64; return [1.001 * Math.cos(phi) * Math.cos(a), 1.001 * Math.sin(phi), 1.001 * Math.cos(phi) * Math.sin(a)]; });
  const faint = new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.13 });
  for (const phi of [-60, -30, 30, 60]) scene.add(new THREE.Line(circle(parallel(phi * D)), faint)).renderOrder = 1;
  for (let k = 0; k < 12; k++) scene.add(new THREE.Line(circle(meridian(k * 30 * D)), faint)).renderOrder = 1;
  scene.add(new THREE.Line(circle(parallel(0)), new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.42 }))).renderOrder = 1;

  const axis = new Arrows(scene.group, 1, { radius: 0.007, head: 0.09, headRadius: 0.032 });
  axis.begin();
  axis.push(0, AXIS[0], 0, 0, AXIS[1] - AXIS[0], 0, [0.75, 0.75, 0.78]);
  axis.end();
  const cameraLongitude = () => Math.atan2(-Math.cos(scene.yaw), Math.sin(scene.yaw));

  const parcel = scene.add(new THREE.Mesh(new THREE.SphereGeometry(0.03, 20, 14), new THREE.MeshBasicMaterial({ color: 0xffffff })));
  const segments = (count, opacity) => new THREE.LineSegments(new THREE.BufferGeometry().setAttribute('position', new THREE.BufferAttribute(new Float32Array(6 * count), 3)), new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity, depthWrite: false }));
  const plumb = scene.add(segments(1, 0.3));
  const pin = scene.add(new THREE.Group());
  const plane = new THREE.Mesh(new THREE.CircleGeometry(1, 48).rotateX(-Math.PI / 2), new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.07, depthWrite: false, side: THREE.DoubleSide }));
  const rim = new THREE.LineLoop(circle(Array.from({ length: 64 }, (_, k) => [Math.cos(2 * Math.PI * k / 64), 0, Math.sin(2 * Math.PI * k / 64)])), new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.45 }));
  const guides = segments(2, 0.6);
  pin.add(plane, rim, guides);
  plumb.renderOrder = plane.renderOrder = rim.renderOrder = guides.renderOrder = 3;
  plumb.frustumCulled = guides.frustumCulled = false;
  const velocity = new Arrows(pin, 1, { radius: 0.009, head: 0.07, headRadius: 0.028, headShare: 0.4 });
  const whole = new Arrows(pin, 1, { radius: 0.0075, head: 0.07, headRadius: 0.026, dashed: true, headShare: 0.4 });
  const split = new Arrows(pin, 2, { radius: 0.0135, head: 0.07, headRadius: 0.036, dashed: true, headShare: 0.4 });

  const ribbon = new Float32Array(CAPACITY * 6), index = new Uint32Array((CAPACITY - 1) * 6);
  for (let i = 0; i < CAPACITY - 1; i++) index.set([2 * i, 2 * i + 1, 2 * i + 2, 2 * i + 1, 2 * i + 3, 2 * i + 2], 6 * i);
  const path = new THREE.BufferGeometry().setAttribute('position', new THREE.BufferAttribute(ribbon, 3).setUsage(THREE.DynamicDrawUsage)).setIndex(new THREE.BufferAttribute(index, 1));
  const [ar, ag, ab] = cssColor(ACCENT);
  scene.add(new THREE.Mesh(path, new THREE.MeshBasicMaterial({ color: new THREE.Color(linear(ar), linear(ag), linear(ab)), side: THREE.DoubleSide }))).frustumCulled = false;

  slider(controls, { label: 'Latitude', min: -90, max: 90, step: 1, value: lat, format: latitudeName, onInput: (v) => { lat = v; scene.pitch = (v < 0 || (!v && scene.pitch < 0) ? -1 : 1) * tilt(v); course.value = heading; relaunch(); } });
  const course = slider(controls, { label: 'Heading', min: 0, max: 359, step: 1, value: heading, format: (v) => Math.abs(lat) === 90 ? `${v}°` : `${v}°, ${COMPASS[Math.round(v / 22.5) % 16]}`, onInput: (v) => { heading = v; relaunch(); } });
  slider(controls, { label: 'Speed', min: 10, max: 300, step: 10, value: speed, format: (v) => `${v} m/s`, onInput: (v) => { speed = v; relaunch(); } });
  const [, toggle] = buttons(controls, [['Launch again', () => { relaunch(); scene.fig.play(true); toggle.textContent = 'Pause'; }], ['Pause', (b) => { scene.fig.play(!scene.fig.running); b.textContent = scene.fig.running ? 'Pause' : 'Run'; }]]);
  const out = readout(controls);

  function relaunch() {
    state = launch(lat * D, cameraLongitude() + OFFSET, heading * D, speed);
    trail = [[0, ...place(state.x)]];
    pending = 0;
    show(1);
    scene.fig.render();
  }

  function step(dt) {
    for (pending += dt * SPEED; pending >= H; pending -= H) {
      advance(state, H);
      trail.push([state.t, ...place(state.x)]);
    }
    const old = Math.max(trail.findIndex((p) => p[0] >= state.t - TRAIL), trail.length - CAPACITY);
    if (old > 0) trail.splice(0, old);
    show(6);
  }

  function show(every) {
    const phi = latitude(state.x), { surface, vertical } = parts(state.x, state.v), f = 2 * OMEGA * Math.sin(phi), deg = phi / D;
    const days = Math.floor(state.t / 86400), hours = Math.floor(state.t % 86400 / 3600);
    const ground = Math.hypot(...surface) * 3600, lift = Math.abs(vertical) * 3600, wide = Math.abs(deg) < 3 || speed / Math.abs(f) > RADIUS / 2;
    out.set([
      ['latitude now', latitudeName(Math.round(deg))],
      ['along the ground', ground < 0.05 ? 'none' : `${ground.toFixed(1)} m/s an hour, to the ${phi > 0 ? 'right' : 'left'}`],
      ['along the vertical', lift < 0.05 ? 'none' : `${lift.toFixed(1)} m/s an hour ${vertical > 0 ? 'upward' : 'downward'}, 1 part in ${rounded(G / Math.abs(vertical))} of gravity`],
      ['inertial period here', Math.abs(deg) < 3 ? 'none this close to the equator' : `${(2 * Math.PI / Math.abs(f) / 3600).toLocaleString('en-US', { maximumFractionDigits: 1, minimumFractionDigits: 1 })} h`],
      ['turning radius here', Math.abs(deg) < 3 ? 'no circle this close to the equator' : wide ? 'more than half the planet' : `${Math.round(speed / Math.abs(f) / 1000).toLocaleString('en-US')} km`],
      ['time', days ? `${days} d ${hours} h` : `${hours} h`],
    ], every);
  }

  const cool = () => cssColor(anomalyColor(-1)), tips = [];
  function update() {
    if (!state) return;
    const at = place(state.x), { up, full, vertical, surface } = parts(state.x, state.v), k = enlarge(scene.fig.width), length = 0.34 + 0.22 * (speed - 10) / 290, scale = length / (2 * OMEGA * speed);
    const ring = cool(), c = Math.cos(scene.pitch), toward = [c * Math.sin(scene.yaw), Math.sin(scene.pitch), c * Math.cos(scene.yaw)].map((x) => 0.012 * x);
    parcel.position.set(...at);
    pin.position.set(...at); pin.scale.setScalar(k);
    plane.scale.setScalar(length + 0.06); rim.scale.setScalar(length + 0.06);
    plane.material.color.setRGB(linear(ring[0]), linear(ring[1]), linear(ring[2])); rim.material.color.copy(plane.material.color);
    const v = turn(state.v).map((x) => x * length / speed), a = turn(full).map((x) => x * scale), s = turn(surface).map((x) => x * scale), u = turn(up).map((x) => x * vertical * scale);
    const ground = Math.hypot(...s) > 0.04, lift = Math.abs(vertical) * scale > 0.04;
    velocity.begin(); velocity.push(0, 0, 0, ...v, PALE); velocity.end();
    whole.begin(); whole.push(...toward, ...a, ring); whole.end();
    split.begin();
    if (ground) split.push(0, 0, 0, ...s, cssColor(SURFACE));
    if (lift) split.push(0, 0, 0, ...u, cssColor(VERTICAL));
    split.end();
    const g = guides.geometry.attributes.position, lines = [ground && s, lift && u].filter(Boolean);
    lines.forEach((tip, i) => g.array.set([...a, ...tip], 6 * i));
    g.needsUpdate = true;
    guides.geometry.setDrawRange(0, 2 * lines.length);
    const reach = 1 + k * length + 0.12, p = plumb.geometry.attributes.position;
    p.array.set([0, 0, 0, ...at.map((x) => x * reach)]);
    p.needsUpdate = true;
    tips.length = 0; tips.push(...[v, a, s, u].map((tip) => tip.map((x) => x * k)));
    const n = trail.length;
    for (let i = 0; i < n; i++) {
      const q = trail[i], prev = trail[Math.max(0, i - 1)], next = trail[Math.min(n - 1, i + 1)];
      const tx = next[1] - prev[1], ty = next[2] - prev[2], tz = next[3] - prev[3];
      let sx = q[2] * tz - q[3] * ty, sy = q[3] * tx - q[1] * tz, sz = q[1] * ty - q[2] * tx;
      const len = Math.hypot(sx, sy, sz) || 1;
      sx *= WIDTH / len; sy *= WIDTH / len; sz *= WIDTH / len;
      ribbon.set([q[1] * LIFT + sx, q[2] * LIFT + sy, q[3] * LIFT + sz, q[1] * LIFT - sx, q[2] * LIFT - sy, q[3] * LIFT - sz], 6 * i);
    }
    path.attributes.position.needsUpdate = true;
    path.setDrawRange(0, Math.max(0, n - 1) * 6);
  }

  function labels(ctx, w) {
    const { r, globe } = layout(w), style = { color: MUTED, size: 11, halo: HALO };
    ctx.save();
    ctx.beginPath(); ctx.rect(0, 0, w, globe); ctx.clip();
    const north = scene.project(0, 1, 0), south = scene.project(0, -1, 0), top = scene.project(0, AXIS[1] + 0.06, 0);
    const e = cameraLongitude() - 0.75, edge = scene.project(Math.cos(e), 0, -Math.sin(e));
    const near = state ? latitude(state.x) / D : 0;
    let clear = () => true;
    if (state) {
      const at = place(state.x), p = scene.project(...at), name = latitudeName(Math.round(near));
      let dx = 0, dy = 0;
      for (const tip of tips) { const q = scene.project(...at.map((x, i) => x + tip[i])), len = Math.hypot(q.x - p.x, q.y - p.y) || 1; dx += (q.x - p.x) / len; dy += (q.y - p.y) / len; }
      ctx.font = '400 11px system-ui, sans-serif';
      const half = ctx.measureText(name).width / 2, len = Math.hypot(dx, dy), ux = len > 0.3 ? -dx / len : -0.7, uy = len > 0.3 ? -dy / len : 0.7, x = Math.min(w - half - 4, Math.max(half + 4, p.x + ux * (half + 10))), y = p.y + uy * 14;
      text(ctx, name, x, y, { ...style, color: '#fff', align: 'center' });
      clear = (cx, cy, size) => Math.abs(cx - x) > half + size || Math.abs(cy - y) > 14;
    }
    if (near < 80 && clear(north.x + 9, north.y - 7, 8)) text(ctx, 'N', north.x + 9, north.y - 7, { ...style, color: INK, weight: 600 });
    if (near > -80 && clear(south.x + 9, south.y + 8, 8)) text(ctx, 'S', south.x + 9, south.y + 8, style);
    text(ctx, 'Ω', top.x, Math.max(10, top.y - 9), { ...style, color: INK, size: 14, align: 'center', weight: 600 });
    const ex = Math.max(30, edge.x), ey = edge.y + (scene.pitch >= 0 ? 11 : -11);
    if (clear(ex, ey, 24)) text(ctx, 'equator', ex, ey, { ...style, align: 'center' });
    ctx.restore();
    const yc = globe + 10 + r, base = -cameraLongitude() - Math.PI / 2;
    polar(ctx, w * 0.25, yc, r, base, true, 'seen from above the pole, standing still');
    polar(ctx, w * 0.75, yc, r, base, false, 'seen from above the pole, turning with the planet');
  }

  function polar(ctx, xc, yc, r, base, fixed, label) {
    const now = state ? state.t : 0, spin = base + (fixed ? OMEGA * now : 0);
    const at = (x, y, a) => [xc + r * (x * Math.cos(a) - y * Math.sin(a)), yc - r * (x * Math.sin(a) + y * Math.cos(a))];
    ctx.save();
    ctx.fillStyle = 'rgba(255,255,255,0.06)'; ctx.beginPath(); ctx.arc(xc, yc, r, 0, 2 * Math.PI); ctx.fill();
    ctx.strokeStyle = 'rgba(255,255,255,0.15)'; ctx.lineWidth = 1;
    for (const phi of [30, 60]) { ctx.beginPath(); ctx.arc(xc, yc, r * Math.cos(phi * D), 0, 2 * Math.PI); ctx.stroke(); }
    for (let k = 1; k < 12; k++) { const [x, y] = at(Math.cos(MARK + k * Math.PI / 6), Math.sin(MARK + k * Math.PI / 6), spin); ctx.beginPath(); ctx.moveTo(xc, yc); ctx.lineTo(x, y); ctx.stroke(); }
    const [mx, my] = at(Math.cos(MARK), Math.sin(MARK), spin);
    ctx.strokeStyle = 'rgba(255,255,255,0.5)'; ctx.beginPath(); ctx.moveTo(xc, yc); ctx.lineTo(mx, my); ctx.stroke();
    ctx.fillStyle = 'rgba(255,255,255,0.5)'; ctx.beginPath(); ctx.arc(xc + 0.92 * (mx - xc), yc + 0.92 * (my - yc), 3.5, 0, 2 * Math.PI); ctx.fill();
    ctx.strokeStyle = LINE; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(xc, yc, r, 0, 2 * Math.PI); ctx.stroke();
    ctx.strokeStyle = ACCENT; ctx.lineWidth = 1.75; ctx.lineJoin = 'round'; ctx.lineCap = 'round';
    let south = null;
    trail.forEach(([t, x, z, y], i) => {
      const p = at(x, -y, base + (fixed ? OMEGA * t : 0)), below = z < SOUTH;
      if (i && below === south) { ctx.lineTo(...p); return; }
      if (i) { ctx.lineTo(...p); ctx.stroke(); }
      ctx.globalAlpha = below ? FAR : 1;
      ctx.beginPath(); ctx.moveTo(...p);
      south = below;
    });
    if (trail.length) ctx.stroke();
    ctx.globalAlpha = 1;
    ctx.fillStyle = INK; ctx.beginPath(); ctx.arc(xc, yc, 2.5, 0, 2 * Math.PI); ctx.fill();
    const here = state ? at(...place(state.x).filter((_, i) => i !== 1).map((c, i) => (i ? -c : c)), spin) : [xc + 1, yc - 1];
    const away = Math.atan2(here[1] - yc, here[0] - xc) + Math.PI;
    text(ctx, 'N', xc + 9 * Math.cos(away), yc + 9 * Math.sin(away), { color: INK, size: 11, weight: 600, halo: HALO, align: 'center' });
    if (state) {
      const [x, z, y] = place(state.x), p = at(x, -y, spin);
      ctx.globalAlpha = z < SOUTH ? FAR : 1;
      ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(...p, 4.5, 0, 2 * Math.PI); ctx.fill();
    }
    ctx.restore();
    ctx.font = '400 11px system-ui, sans-serif';
    label.split(', ').forEach((line, i) => text(ctx, line, xc, yc + r + 16 + 13 * i, { align: 'center', color: MUTED, size: 11 }));
  }

  relaunch();
  scene.fig.play(true);
}
