import { slider, buttons, legend, caption, readout, rampRGB, clamp, INK, DARK_NEUTRAL } from '../runtime.module.js';
import { simulationGlobe } from '../swGlobe.module.js';
import { EARTH } from '../swCases.module.js';

const N = 24, DEPTH = 4000, HEIGHT = 100, WIDTH = 1.5e6, SPEED = 6 * 3600, SCALE = 25;

export function mountBump(root) {
  const controls = root.querySelector('.controls');
  let spin = 1;
  const sim = simulationGlobe(root, {
    N, speed: SPEED, center: { lat: Math.PI / 4, lon: 0 }, arrowScale: 0.012,
    colorOf: (value, rgb) => rampRGB(0.5 + clamp(value / SCALE, -1, 1) / 2, rgb, DARK_NEUTRAL),
    onFrame: (frame) => {
      let top = -Infinity;
      for (const v of frame.field) top = Math.max(top, v);
      const hours = Math.round(frame.time / 3600), f = 2 * EARTH.omega * spin * Math.SQRT1_2, c = Math.sqrt(EARTH.g * DEPTH);
      out.set([['time', `${Math.floor(hours / 24)} d ${String(hours % 24).padStart(2, '0')} h`], ['highest point of the surface', `${top.toFixed(0)} m above the mean`], ['gravity waves travel at', `${c.toFixed(0)} m/s`], ['Rossby radius at 45°', spin ? `${(c / f / 1000).toFixed(0)} km` : 'unlimited: no spin']], 4);
    },
  });
  caption(root, `A ${DEPTH / 1000} km deep layer of water at rest covers the globe, with a hump ${HEIGHT} m high and about ${(2 * WIDTH / 1e6).toFixed(0)},000 km across dropped onto it at 45° north. Time runs ${SPEED / 3600} hours per second. Drag to turn the globe.`);
  legend(root, [['ramp', `surface ${SCALE} m below to ${SCALE} m above the mean`, 'cool', 'warm', 'rgb(52, 55, 62)'], ['arrow', 'the water’s motion', INK]]);
  slider(controls, { label: 'Spin of the planet', min: 0, max: 4, step: 0.5, value: spin, format: (v) => v === 0 ? 'none' : v === 1 ? 'the Earth’s' : `${v}× the Earth’s`, onInput: (v) => { spin = v; start(); } });
  buttons(controls, [['Drop it again', () => start()]]);
  const out = readout(controls);
  const start = () => sim.start({ kind: 'bump', options: { depth: DEPTH, height: HEIGHT, width: WIDTH }, spin, closureHours: 12, coarse: 8 });
  start();
}
