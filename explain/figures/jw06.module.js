import { playbackGlobe } from '../playbackGlobe.module.js';
import { LAYER_FIELDS } from '../layerCases.module.js';

export function mountJw06(root) {
  playbackGlobe(root, {
    url: new URL('../data/jw06_N16.bin', import.meta.url).href, N: 16, center: { lat: 50 * Math.PI / 180, lon: 160 * Math.PI / 180 },
    sequences: [{ name: 'wave', label: 'twelve days' }], fields: LAYER_FIELDS, framesPerSecond: 6,
    note: 'Computed in advance by the model’s layered core, 27 layers on the N = 16 grid, with a frame every three hours. Drag the globe to turn it, or the time slider to step through the days.',
  });
}
