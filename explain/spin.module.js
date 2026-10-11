import { mountTurntable } from './figures/turntable.module.js';
import { mountCoriolis } from './figures/coriolis.module.js';
import { mountInertial } from './figures/inertial.module.js';
import { paletteControl } from './runtime.module.js';
import { previews } from './previews.module.js';

const mounts = { turntable: mountTurntable, coriolis: mountCoriolis, inertial: mountInertial };
for (const figure of document.querySelectorAll('figure.fig[data-figure]')) mounts[figure.dataset.figure]?.(figure);
paletteControl(document.querySelector('nav.parts .palette'));
previews(document);
