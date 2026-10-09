import { mountPiston } from './figures/piston.module.js';
import { mountLift } from './figures/lift.module.js';
import { mountColumn } from './figures/column.module.js';
import { mountBob } from './figures/bob.module.js';
import { mountCells } from './figures/cells.module.js';
import { paletteControl } from './runtime.module.js';
import { previews } from './previews.module.js';

const mounts = { piston: mountPiston, lift: mountLift, column: mountColumn, bob: mountBob, cells: mountCells };
for (const figure of document.querySelectorAll('figure.fig[data-figure]')) mounts[figure.dataset.figure]?.(figure);
paletteControl(document.querySelector('nav.parts .palette'));
previews(document);
