import { mountMoistlift } from './figures/moistlift.module.js';
import { mountVaporslices } from './figures/vaporslices.module.js';
import { mountEvaporation } from './figures/evaporation.module.js';
import { mountWetcolumn } from './figures/wetcolumn.module.js';
import { mountPlume } from './figures/plume.module.js';
import { mountCloudlayer } from './figures/cloudlayer.module.js';
import { paletteControl } from './runtime.module.js';
import { previews } from './previews.module.js';

const mounts = { moistlift: mountMoistlift, vaporslices: mountVaporslices, evaporation: mountEvaporation, wetcolumn: mountWetcolumn, plume: mountPlume, cloudlayer: mountCloudlayer };
for (const figure of document.querySelectorAll('figure.fig[data-figure]')) mounts[figure.dataset.figure]?.(figure);
paletteControl(document.querySelector('nav.parts .palette'));
previews(document);
