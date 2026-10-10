import { buttons, legend, caption, readout, rampRGB, clamp, INK, DARK_NEUTRAL } from '../runtime.module.js';
import { simulationGlobe } from '../swGlobe.module.js';

const N = 24, SPEED = 12 * 3600, LOW = 0, HIGH = 2600, FLOW = 7.848e-6 * 86400 * 180 / Math.PI;

export function mountHaurwitz(root) {
  const controls = root.querySelector('.controls');
  let last = null, travelled = 0;
  const sim = simulationGlobe(root, {
    N, speed: SPEED, center: { lat: 35 * Math.PI / 180, lon: 0 }, arrowScale: 0.0007, arrowOpacity: 0.6,
    colorOf: (value, rgb) => rampRGB(clamp((value - LOW) / (HIGH - LOW), 0, 1), rgb, DARK_NEUTRAL),
    onFrame: (frame) => {
      if (frame.time === 0) { last = frame.phase; travelled = 0; }
      else if (last !== null) { let d = frame.phase - last; if (d > Math.PI / 4) d -= Math.PI / 2; if (d < -Math.PI / 4) d += Math.PI / 2; travelled += d; last = frame.phase; }
      const days = frame.time / 86400, moved = travelled * 180 / Math.PI;
      out.set([['day', days.toFixed(1)], ['the pattern has moved', `${moved.toFixed(0)}° east, ${days > 0.5 ? (moved / days).toFixed(1) : '…'}° a day`], ['the wind carries the fluid east at', `${FLOW.toFixed(0)}° a day`]], 4);
    },
  });
  caption(root, `A layer of water 8 km deep at its shallowest, with four highs and four lows around each hemisphere and the wind they balance. The grid is N = ${N}. Time runs ${SPEED / 3600} hours per second.`);
  legend(root, [['ramp', `surface from 8 km to ${(8 + HIGH / 1000).toFixed(1)} km deep`, 'cool', 'warm', 'rgb(52, 55, 62)'], ['arrow', 'the wind', INK]]);
  buttons(controls, [['Restart', () => start()]]);
  const out = readout(controls);
  const start = () => sim.start({ kind: 'haurwitz', spin: 1, closureHours: 0, coarse: 8 });
  start();
}
