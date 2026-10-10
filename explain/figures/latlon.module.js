import * as THREE from '../../js/three.module.js';
import { Globe, rgb } from '../globe.module.js';
import { slider, legend, readout, rampColor } from '../runtime.module.js';

const SIZES = [30, 20, 15, 12, 10, 6, 5, 4, 3], RADIUS = 6371;

export function mountLatLon(root) {
  const controls = root.querySelector('.controls');
  const globe = new Globe(root, { height: 440 });
  legend(root, [['ramp', 'cells smaller to larger than the average, on a scale of four times either way', 'cool', 'warm', 'neutral']]);
  let index = 4;
  slider(controls, { label: 'Cell size', min: 0, max: SIZES.length - 1, step: 1, value: index, format: (v) => `${SIZES[v]}° × ${SIZES[v]}°`, onInput: (v) => { index = v; build(); } });
  const out = readout(controls);

  function build() {
    const d = SIZES[index], rows = Math.round(180 / d), cols = Math.round(360 / d), toRad = Math.PI / 180;
    const point = (lat, lon) => new THREE.Vector3(Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat));
    const polygons = [], areas = [];
    for (let j = 0; j < rows; j++) for (let i = 0; i < cols; i++) {
      const lat0 = (-90 + j * d) * toRad, lat1 = lat0 + d * toRad, lon0 = i * d * toRad, lon1 = lon0 + d * toRad;
      const corners = [point(lat0, lon0), point(lat0, lon1), point(lat1, lon1), point(lat1, lon0)];
      polygons.push({ center: point((lat0 + lat1) / 2, (lon0 + lon1) / 2), vertices: corners.filter((v, k) => k === 0 || v.distanceTo(corners[k - 1]) > 1e-9) });
      areas.push((lon1 - lon0) * (Math.sin(lat1) - Math.sin(lat0)));
    }
    const mean = areas.reduce((a, b) => a + b) / areas.length;
    globe.clear();
    globe.cells(polygons, (cell, k) => rgb(rampColor(0.5 + Math.log2(areas[k] / mean) / 4)));
    const paths = [];
    for (let j = 1; j < rows; j++) { const lat = (-90 + j * d) * toRad, ring = []; for (let s = 0; s <= 90; s++) ring.push(point(lat, (s / 90) * 2 * Math.PI)); paths.push(ring); }
    for (let i = 0; i < cols; i++) { const lon = i * d * toRad, meridian = []; for (let s = 0; s <= 36; s++) meridian.push(point(-Math.PI / 2 + (s / 36) * Math.PI, lon)); paths.push(meridian); }
    globe.lines(paths, { color: 0x000000, opacity: 0.4 });
    const widest = RADIUS * d * toRad, narrowest = widest * Math.cos((90 - d) * toRad);
    out.set([['cells', `${polygons.length}`], ['largest cell', `${(Math.max(...areas) / Math.min(...areas)).toFixed(0)}× the smallest`], ['cell width', `${widest.toFixed(0)} km at the equator, ${narrowest.toFixed(0)} km beside the pole`]]);
    globe.fig.render();
  }

  build();
}
