import { mountLatLon } from './figures/latlon.module.js';
import { mountIcosa } from './figures/icosa.module.js';
import { mountRelax } from './figures/relax.module.js';
import { mountStagger } from './figures/stagger.module.js';
import { mountOperators } from './figures/operators.module.js';
import { mountTrisk } from './figures/trisk.module.js';
import { mountThreeColor } from './figures/threecolor.module.js';
import { mountSlices } from './figures/slices.module.js';
import { paletteControl } from './runtime.module.js';
import { previews } from './previews.module.js';

const mounts = { latlon: mountLatLon, icosa: mountIcosa, relax: mountRelax, stagger: mountStagger, operators: mountOperators, trisk: mountTrisk, threecolor: mountThreeColor, slices: mountSlices };
for (const figure of document.querySelectorAll('figure.fig[data-figure]')) mounts[figure.dataset.figure]?.(figure);
paletteControl(document.querySelector('nav.parts .palette'));
previews(document);
