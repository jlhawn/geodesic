import { playbackGlobe } from '../playbackGlobe.module.js';
import { LAYER_FIELDS } from '../layerCases.module.js';

export function mountHeldSuarez(root) {
  playbackGlobe(root, {
    url: new URL('../data/heldsuarez_N16.bin', import.meta.url).href, N: 16, center: { lat: 35 * Math.PI / 180, lon: 0 },
    sequences: [{ name: 'spinup', label: 'the first 120 days, a frame a day' }, { name: 'settled', label: 'a month of the settled climate, four frames a day' }],
    fields: LAYER_FIELDS, framesPerSecond: 8,
    note: 'Computed in advance by the model’s layered core under the Held–Suarez forcing, 27 layers on the N = 16 grid, starting from rest. Drag the globe to turn it, or the time slider to step through the days.',
  });
}
