import { mountStack } from './figures/stack.module.js';
import { mountContinuity } from './figures/continuity.module.js';
import { mountThickness } from './figures/thickness.module.js';
import { mountMountain } from './figures/mountain.module.js';
import { mountInvariant } from './figures/invariant.module.js';
import { mountIsentropes } from './figures/isentropes.module.js';
import { mountStorm3d } from './figures/storm3d.module.js';
import { mountEnergy } from './figures/energy.module.js';
import { paletteControl } from './runtime.module.js';
import { previews } from './previews.module.js';

const mounts = { stack: mountStack, continuity: mountContinuity, thickness: mountThickness, mountain: mountMountain, invariant: mountInvariant, isentropes: mountIsentropes, storm3d: mountStorm3d, energy: mountEnergy };
for (const figure of document.querySelectorAll('figure.fig[data-figure]')) mounts[figure.dataset.figure]?.(figure);
paletteControl(document.querySelector('nav.parts .palette'));
previews(document);
