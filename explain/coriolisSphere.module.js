/*
 * A parcel sliding without friction on a turning sphere, seen in the turning
 * frame, with no pressure force. Earth-centered coordinates in meters with z
 * along the axis toward the north pole and x through longitude 0. Effective
 * gravity and the ground's push are both along the local vertical, so only
 * the surface part of the Coriolis acceleration −2Ω × v turns the parcel.
 */
export const OMEGA = 7.2921e-5, RADIUS = 6371e3;

const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const scale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const add = (a, b, s = 1) => [a[0] + b[0] * s, a[1] + b[1] * s, a[2] + b[2] * s];
const unit = (a) => scale(a, 1 / Math.hypot(...a));

export const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];

export function frame(lat, lon) {
  const cp = Math.cos(lat), sp = Math.sin(lat), cl = Math.cos(lon), sl = Math.sin(lon);
  return { up: [cp * cl, cp * sl, sp], east: [-sl, cl, 0], north: [-sp * cl, -sp * sl, cp] };
}

export const latitude = (x) => Math.asin(Math.max(-1, Math.min(1, x[2] / Math.hypot(...x))));
export const longitude = (x) => Math.atan2(x[1], x[0]);
export const coriolis = (v, omega = OMEGA) => [2 * omega * v[1], -2 * omega * v[0], 0];

export function parts(x, v, omega = OMEGA) {
  const up = unit(x), full = coriolis(v, omega), vertical = dot(full, up);
  return { up, full, vertical, surface: add(full, up, -vertical) };
}

export function launch(lat, lon, heading, speed, omega = OMEGA, radius = RADIUS) {
  const { up, east, north } = frame(lat, lon);
  return { x: scale(up, radius), v: add(scale(north, speed * Math.cos(heading)), east, speed * Math.sin(heading)), t: 0, omega };
}

function rate(x, v, omega) {
  const up = unit(x), a = coriolis(v, omega);
  return add(add(a, up, -dot(a, up)), x, -dot(v, v) / dot(x, x));
}

export function step(state, dt, omega = state.omega ?? OMEGA, radius = RADIUS) {
  const { x, v } = state;
  const k1x = v, k1v = rate(x, v, omega);
  const x2 = add(x, k1x, dt / 2), v2 = add(v, k1v, dt / 2), k2v = rate(x2, v2, omega);
  const x3 = add(x, v2, dt / 2), v3 = add(v, k2v, dt / 2), k3v = rate(x3, v3, omega);
  const x4 = add(x, v3, dt), v4 = add(v, k3v, dt), k4v = rate(x4, v4, omega);
  const nx = scale(unit(add(add(add(add(x, k1x, dt / 6), v2, dt / 3), v3, dt / 3), v4, dt / 6)), radius);
  const nv = add(add(add(add(v, k1v, dt / 6), k2v, dt / 3), k3v, dt / 3), k4v, dt / 6), up = scale(nx, 1 / radius);
  state.x = nx;
  state.v = add(nv, up, -dot(nv, up));
  state.t += dt;
  return state;
}
