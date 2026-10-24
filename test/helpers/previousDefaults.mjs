import { sigmaInterfaces } from '../../js/dynamics/sigmaCore.module.js';

/*
 * The settings the pinned digests and figures were made with, where the
 * model's defaults differ.
 */
export const PREVIOUS_LEVELS = sigmaInterfaces('cam26');
export const PREVIOUS_RADIATION = { cloudColumns: false, cloudSpread: 0, varianceCover: 'deficit', buoyancyLength: null, deckResolved: 'added' };
export const PREVIOUS_MOIST = { stratocumulusInversion: null, cumulusEntrainmentHeight: null, cumulusAscent: 'buoyant' };
export const PREVIOUS_OCEAN = { rotationLimit: 'algebraic', langmuirIce: 'none', shelfMeltLayer: 0, iceDrift: false };
