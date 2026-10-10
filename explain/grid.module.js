import { mountLatLon } from './figures/latlon.module.js';
import { mountIcosa } from './figures/icosa.module.js';
import { mountRelax } from './figures/relax.module.js';
import { paletteControl } from './runtime.module.js';
import { previews } from './previews.module.js';

const mounts = { latlon: mountLatLon, icosa: mountIcosa, relax: mountRelax };
for (const figure of document.querySelectorAll('figure.fig[data-figure]')) mounts[figure.dataset.figure]?.(figure);
paletteControl(document.querySelector('nav.parts .palette'));
previews(document);
