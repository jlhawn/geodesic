import { choice, buttons, legend, caption, readout, rampRGB, clamp, INK, DARK_NEUTRAL } from '../runtime.module.js';
import { simulationGlobe } from '../swGlobe.module.js';

const N = 32, SPEED = 6 * 3600, SCALE = 1.0e-4;

export function mountGalewsky(root) {
  const controls = root.querySelector('.controls');
  let perturbed = true;
  const sim = simulationGlobe(root, {
    N, speed: SPEED, center: { lat: 50 * Math.PI / 180, lon: 0 }, arrowScale: 0.0009, arrowOpacity: 0.6,
    colorOf: (value, rgb) => rampRGB(0.5 + clamp(value / SCALE, -1, 1) / 2, rgb, DARK_NEUTRAL),
    onFrame: (frame) => {
      let most = 0;
      for (const v of frame.field) most = Math.max(most, Math.abs(v));
      out.set([['day', `${(frame.time / 86400).toFixed(1)}`], ['strongest spin', `${(most * 1e5).toFixed(1)} × 10⁻⁵ per second`], ['water gained or lost', Math.abs(frame.massDrift) < 1e-13 ? 'none, to 13 decimal places' : `${(frame.massDrift * 100).toExponential(0)} %`]], 4);
    },
  });
  caption(root, `A jet blowing east at up to 80 m/s around 45° north, on a 10 km deep layer of water, balanced exactly by the slope of the surface. The grid is N = ${N}, 10,242 cells about 220 km across. Time runs ${SPEED / 3600} hours per second.`);
  legend(root, [['ramp', 'spinning clockwise to counterclockwise', 'cool', 'warm', 'rgb(52, 55, 62)'], ['arrow', 'the wind', INK]]);
  choice(controls, { label: 'The jet starts', options: [['with a small bump in it', 'bump'], ['perfectly smooth', 'smooth']], value: 'bump', onChange: (v) => { perturbed = v === 'bump'; start(); }, span: true });
  buttons(controls, [['Restart', () => start()]]);
  const out = readout(controls);
  const start = () => sim.start({ kind: 'galewsky', options: { perturbed }, spin: 1, closureHours: 3, coarse: 10 });
  start();
}
