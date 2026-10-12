import { mountSeastore } from './figures/seastore.module.js';
import { mountSeawater } from './figures/seawater.module.js';
import { mountEkman } from './figures/ekman.module.js';
import { mountSeasurface } from './figures/seasurface.module.js';
import { mountGyre } from './figures/gyre.module.js';
import { mountSeaice } from './figures/seaice.module.js';
import { paletteControl } from './runtime.module.js';
import { previews } from './previews.module.js';

const mounts = { seastore: mountSeastore, seawater: mountSeawater, ekman: mountEkman, seasurface: mountSeasurface, gyre: mountGyre, seaice: mountSeaice };
for (const figure of document.querySelectorAll('figure.fig[data-figure]')) mounts[figure.dataset.figure]?.(figure);
paletteControl(document.querySelector('nav.parts .palette'));
previews(document);
