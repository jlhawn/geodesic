import { mountBump } from './figures/bump.module.js';
import { mountGalewsky } from './figures/galewsky.module.js';
import { mountHaurwitz } from './figures/haurwitz.module.js';
import { mountJw06 } from './figures/jw06.module.js';
import { mountJwSection, mountHsSection } from './figures/sections.module.js';
import { mountHeldSuarez } from './figures/heldsuarez.module.js';
import { paletteControl } from './runtime.module.js';
import { previews } from './previews.module.js';

const mounts = { bump: mountBump, galewsky: mountGalewsky, haurwitz: mountHaurwitz, jwsection: mountJwSection, jw06: mountJw06, hssection: mountHsSection, heldsuarez: mountHeldSuarez };
for (const figure of document.querySelectorAll('figure.fig[data-figure]')) mounts[figure.dataset.figure]?.(figure);
paletteControl(document.querySelector('nav.parts .palette'));
previews(document);
