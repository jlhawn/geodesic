import { test } from 'node:test';
import { steppedDraft } from './helpers/iceShelf.mjs';

test('under a shelf whose draft steps 300 m from cell to cell over 200 m of water, where neighbors share no depth, every edge passes the sill floor, and with the melt off the ocean stays as still for a month as with the shelves as land, homogeneous or on ISOMIP+\'s warm profile, on the CPU', () => steppedDraft(['cpu']));
