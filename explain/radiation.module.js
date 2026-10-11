import { mountSunlit } from './figures/sunlit.module.js';
import { mountSeasons } from './figures/seasons.module.js';
import { mountBeam } from './figures/beam.module.js';
import { mountInfrared } from './figures/infrared.module.js';
import { mountEquilibrium } from './figures/equilibrium.module.js';
import { mountLatitudes } from './figures/latitudes.module.js';
import { mountSunlitsphere } from './figures/sunlitsphere.module.js';
import { mountTransport } from './figures/transport.module.js';
import { paletteControl } from './runtime.module.js';
import { previews } from './previews.module.js';

const mounts = { sunlit: mountSunlit, seasons: mountSeasons, beam: mountBeam, infrared: mountInfrared, equilibrium: mountEquilibrium, latitudes: mountLatitudes, sunlitsphere: mountSunlitsphere, transport: mountTransport };
for (const figure of document.querySelectorAll('figure.fig[data-figure]')) mounts[figure.dataset.figure]?.(figure);
paletteControl(document.querySelector('nav.parts .palette'));
previews(document);
