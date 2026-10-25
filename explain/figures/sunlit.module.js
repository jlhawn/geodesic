import * as THREE from '../../js/three.module.js';
import { Grid } from '../../js/grid.module.js';
import { buildMesh } from '../../js/mesh.module.js';
import { sunDirection, AXIAL_TILT } from '../../js/physics/radiation.module.js';
import { Globe } from '../globe.module.js';
import { PUBLISHED_RADIATION } from '../published.module.js';
import { slider, buttons, legend, caption, readout, text, sequentialRGB, paletteVersion, ACCENT } from '../runtime.module.js';
import { SOLAR_CONSTANT, timeOf, subsolarPoint, sunlightAt, sunlightOnCells, areaMean } from '../subsolar.module.js';

const D = Math.PI / 180, N = 32, HOURS_PER_SECOND = 3, VIEW = { lat: 25 * D, lon: -50 * D };
const SEASONS = { 0: 'March equinox', 91: 'June solstice', 183: 'September equinox', 274: 'December solstice' };
const HALO = 'rgba(20, 20, 22, 0.85)', BRIGHT = 'rgba(255, 255, 255, 0.6)';
const NORTH_UP = new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 0, 1), new THREE.Vector3(0, 1, 0));

const point = (lat, lon) => new THREE.Vector3(Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat));
const parallel = (lat, from = -180, to = 180) => Array.from({ length: Math.round((to - from) / 2) + 1 }, (_, k) => point(lat, (from + 2 * k) * D));
const meridian = (lon) => Array.from({ length: 91 }, (_, k) => point((-90 + 2 * k) * D, lon));
const dashed = (lat) => Array.from({ length: 45 }, (_, k) => parallel(lat, -180 + 8 * k, -176 + 8 * k));

const clock = (hour) => {
  const minutes = Math.round(hour * 60), name = minutes === 720 ? ', noon' : minutes % 1440 === 0 ? ', midnight' : '';
  return `${String(Math.floor(minutes / 60)).padStart(2, '0')}:${String(minutes % 60).padStart(2, '0')}${name}`;
};
const angle = (radians, digits, plus, minus) => {
  const v = Number((radians / D).toFixed(digits)), size = Math.abs(v).toFixed(digits);
  return v === 0 || Math.abs(v) === 180 ? `${size}°` : `${size}° ${v > 0 ? plus : minus}`;
};
const place = (lat, lon, digits) => `${angle(lat, digits, 'N', 'S')}, ${angle(lon, digits, 'E', 'W')}`;
const watts = (v) => `${v.toLocaleString('en-US', { minimumFractionDigits: 1, maximumFractionDigits: 1 })} W/m²`;

export function mountSunlit(root) {
  const controls = root.querySelector('.controls');
  const globe = new Globe(root, { height: 440, center: VIEW, lighting: { directional: 0, ambient: Math.PI } });
  const grid = new Grid(N), mesh = buildMesh(grid);
  const paint = globe.dynamicCells([...grid].map((c) => ({ center: c.centerVertex, vertices: c.vertices })));
  globe.lines([...[-60, -30, 30, 60].map((lat) => parallel(lat * D)), ...Array.from({ length: 11 }, (_, k) => meridian((k + 1) * 30 * D))], { opacity: 0.12 });
  globe.lines([parallel(0), meridian(0)], { opacity: 0.6 });
  globe.lines([...dashed(AXIAL_TILT), ...dashed(-AXIAL_TILT)], { opacity: 0.6 });

  caption(root, `The model’s sun on a coarser copy of its grid, ${mesh.nCells.toLocaleString('en-US')} cells, each colored by the sunlight arriving at the top of the atmosphere over it. The dot marks the subsolar point, where the sun is straight overhead; the readout reports the sunlight under the cross at the center of the globe. Drag to turn the globe.`);
  legend(root, [['ramp', `sunlight arriving, 0 to ${SOLAR_CONSTANT.toLocaleString('en-US')} W/m²`, 'rgb(52, 55, 62)', 'warm'], ['faint', 'the terminator, where the sun is on the horizon', ACCENT], ['faint', 'the equator and the Greenwich meridian', BRIGHT], ['dash', `the tropics, ${(AXIAL_TILT / D).toFixed(2)}° north and south`, BRIGHT]]);

  const overlay = document.createElement('canvas');
  Object.assign(overlay.style, { position: 'absolute', left: '0', top: '0', width: '100%', pointerEvents: 'none' });
  globe.fig.stage.style.position = 'relative';
  globe.fig.stage.append(overlay);

  let day = 0, hour = 12, running = false, painted = null;
  const sun = new Float64Array(3), values = new Float64Array(mesh.nCells), view = new THREE.Quaternion(), probe = new THREE.Vector3();
  let mean = 0;

  function refresh() {
    const t = timeOf(day, hour), key = `${t}/${paletteVersion}`;
    if (key === painted) return;
    sunDirection(t, sun, PUBLISHED_RADIATION.orbit);
    sunlightOnCells(mesh.xCell, sun, values);
    mean = areaMean(values, mesh.areaCell);
    paint((i, rgb) => sequentialRGB(values[i] / SOLAR_CONSTANT, rgb));
    painted = key;
  }

  function project(x, y, z, w, h) {
    probe.set(x, y, z).applyQuaternion(view);
    const scale = h / (2 * globe.zoom);
    return [w / 2 + probe.x * scale, h / 2 - probe.y * scale, probe.z];
  }

  function drawOverlay() {
    const { width: w, h, dpr } = globe.fig;
    if (overlay.width !== Math.round(w * dpr) || overlay.height !== Math.round(h * dpr)) {
      overlay.width = Math.round(w * dpr);
      overlay.height = Math.round(h * dpr);
      overlay.style.height = `${h}px`;
    }
    const ctx = overlay.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    const s = new THREE.Vector3(...sun), e1 = new THREE.Vector3(0, 0, 1).cross(s);
    if (e1.lengthSq() < 1e-12) e1.set(1, 0, 0); else e1.normalize();
    const e2 = s.clone().cross(e1);
    ctx.save();
    ctx.lineCap = 'round';
    for (const [color, width] of [[HALO, 4], [ACCENT, 2]]) {
      ctx.strokeStyle = color; ctx.lineWidth = width;
      ctx.beginPath();
      let last = null;
      for (let k = 0; k <= 256; k++) {
        const a = 2 * Math.PI * k / 256, c = Math.cos(a), n = Math.sin(a);
        const p = project(c * e1.x + n * e2.x, c * e1.y + n * e2.y, c * e1.z + n * e2.z, w, h);
        if (p[2] > -1e-9) { if (last) ctx.lineTo(p[0], p[1]); else ctx.moveTo(p[0], p[1]); last = p; } else last = null;
      }
      ctx.stroke();
    }
    const [sx, sy, front] = project(sun[0], sun[1], sun[2], w, h);
    if (front > 0) {
      ctx.fillStyle = ACCENT; ctx.strokeStyle = HALO; ctx.lineWidth = 2.5;
      ctx.beginPath(); ctx.arc(sx, sy, 5, 0, 2 * Math.PI); ctx.stroke(); ctx.fill();
      ctx.font = '400 11px system-ui, sans-serif';
      const left = sx + 9 + ctx.measureText('sun overhead').width > w - 4;
      text(ctx, 'sun overhead', left ? sx - 9 : sx + 9, sy < 18 ? sy + 11 : sy - 9, { color: ACCENT, size: 11, halo: HALO, align: left ? 'right' : 'left' });
    }
    for (const [color, width] of [[HALO, 4], ['#fff', 1.5]]) {
      ctx.strokeStyle = color; ctx.lineWidth = width;
      ctx.beginPath();
      for (const [dx, dy] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) { ctx.moveTo(w / 2 + 3 * dx, h / 2 + 3 * dy); ctx.lineTo(w / 2 + 9 * dx, h / 2 + 9 * dy); }
      ctx.stroke();
    }
    ctx.restore();
  }

  function report() {
    const center = new THREE.Vector3(0, 0, 1).applyQuaternion(view.clone().invert());
    const lat = Math.asin(Math.max(-1, Math.min(1, center.z))), lon = Math.atan2(center.y, center.x);
    const there = sunlightAt(center.x, center.y, center.z, sun), height = Math.asin(there / SOLAR_CONSTANT) / D, below = there <= 0;
    const overhead = subsolarPoint(sun);
    out.set([
      ['sun overhead at', place(overhead.lat, overhead.lon, 2)],
      [`under the cross, ${place(lat, lon, 1)}:`, watts(there)],
      ['sun there', below ? 'below the horizon' : `${height.toFixed(1)}° above the horizon`],
      ['mean over the globe', watts(mean)],
    ]);
  }

  const render = globe.render.bind(globe);
  globe.render = () => {
    if (!globe.fig.width) return;
    refresh();
    render();
    view.copy(globe.group.quaternion);
    drawOverlay();
    report();
  };

  globe.fig.step = (dt) => {
    if (!running) return false;
    hour += dt * HOURS_PER_SECOND;
    if (hour >= 24) { hour -= 24; day = (day + 1) % 365; days.value = day; }
    hours.value = hour;
  };

  const days = slider(controls, { label: 'Day of the year', min: 0, max: 364, step: 1, value: day, format: (v) => SEASONS[v] ? `${v}, ${SEASONS[v]}` : `${v}`, onInput: (v) => { day = v; globe.fig.render(); } });
  const hours = slider(controls, { label: 'Hour at Greenwich', min: 0, max: 24, step: 0.25, value: hour, format: clock, onInput: (v) => { hour = v; globe.fig.render(); } });
  for (const output of controls.querySelectorAll('.control output')) output.style.minWidth = '11.5em';
  buttons(controls, [
    ['Run', (b) => { running = !running; b.textContent = running ? 'Pause' : 'Run'; globe.fig.play(running); }],
    ['Look from the sun', () => {
      const { lat, lon } = subsolarPoint(sun);
      const target = new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1, 0, 0), lat)
        .multiply(new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0, 1, 0), -Math.PI / 2 - lon)).multiply(NORTH_UP);
      globe.spin.copy(target.multiply(globe.base.clone().invert()));
      globe.fig.render();
    }],
  ]);
  const out = readout(controls);
}
