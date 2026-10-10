import { mountBump } from './figures/bump.module.js';
import { mountGalewsky } from './figures/galewsky.module.js';
import { paletteControl } from './runtime.module.js';
import { previews } from './previews.module.js';

const mounts = { bump: mountBump, galewsky: mountGalewsky };
for (const figure of document.querySelectorAll('figure.fig[data-figure]')) mounts[figure.dataset.figure]?.(figure);
paletteControl(document.querySelector('nav.parts .palette'));
previews(document);
