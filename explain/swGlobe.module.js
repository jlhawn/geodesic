import { Grid } from '../js/grid.module.js';
import { Globe } from './globe.module.js';
import { paletteVersion } from './runtime.module.js';

export function simulationGlobe(root, { N, height = 460, center = null, zoom = 1.08, speed, colorOf, arrowScale = 0.01, arrowOpacity = 0.85, onFrame = () => {} }) {
  const globe = new Globe(root, { height, center, zoom, lighting: { directional: 1.5, ambient: 1.9 } });
  const polygons = [...new Grid(N)].map((c) => ({ center: c.centerVertex, vertices: c.vertices }));
  const paint = globe.dynamicCells(polygons);
  const worker = new Worker(new URL('./sw.worker.js', import.meta.url), { type: 'module' });
  let latest = null, pending = false, drawArrows = null, centers = null, config = null, run = 0, owed = 0, painted = null, paintedPalette = -1;

  const render = globe.render.bind(globe);
  globe.render = () => {
    if (latest && (latest !== painted || paintedPalette !== paletteVersion)) {
      paint((i, rgb) => colorOf(latest.field[i], rgb));
      if (drawArrows) drawArrows(centers, latest.arrows, arrowScale);
      painted = latest; paintedPalette = paletteVersion;
    }
    render();
  };

  worker.onmessage = ({ data }) => {
    if (data.run !== run) return;
    pending = false;
    if (!drawArrows) {
      drawArrows = globe.dynamicArrows(data.picks.length, { opacity: arrowOpacity });
      centers = new Float32Array(3 * data.picks.length);
      data.picks.forEach((i, k) => { const c = polygons[i].center; centers[3 * k] = c.x; centers[3 * k + 1] = c.y; centers[3 * k + 2] = c.z; });
    }
    latest = data;
    onFrame(data);
    if (!globe.fig.frame) globe.fig.render();
  };

  globe.fig.step = (dt) => {
    owed += dt;
    if (!pending && latest) { pending = true; worker.postMessage({ type: 'advance', seconds: owed * speed }); owed = 0; }
  };
  globe.fig.play(true);

  return {
    globe,
    start(next) { config = { N, ...next }; pending = true; latest = null; owed = 0; worker.postMessage({ type: 'start', run: ++run, ...config }); },
    restart() { if (config) this.start(config); },
    get frame() { return latest; },
  };
}
