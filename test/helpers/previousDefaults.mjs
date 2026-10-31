import { sigmaInterfaces } from '../../js/dynamics/sigmaCore.module.js';

/*
 * The settings the pinned digests and figures were made with, where the
 * model's defaults differ.
 */
export const PREVIOUS_LEVELS = sigmaInterfaces('cam26');
export const PREVIOUS_SHORTWAVE = { shortwave: 'sequential' };
export const PREVIOUS_BERGERON = { bergeronAutoconversion: false };
export const PREVIOUS_AEROSOL = { aerosol: 'flat' };
export const PREVIOUS_RADIATION = { cloudColumns: false, cloudSpread: 0, varianceCover: 'deficit', buoyancyLength: null, deckResolved: 'added', orbit: 'circular', ...PREVIOUS_SHORTWAVE, ...PREVIOUS_AEROSOL };
export const PREVIOUS_MOIST = { stratocumulusInversion: null, cumulusEntrainmentHeight: null, cumulusAscent: 'buoyant', ...PREVIOUS_BERGERON };
export const PREVIOUS_LAND = { soil: 'slab' };
export const PREVIOUS_SURFACE = { seaHumidity: 1 };
export const PREVIOUS_OCEAN = { rotationLimit: 'algebraic', langmuirIce: 'none', shelfMeltLayer: 0, iceDrift: false, mstarScheme: 'fixed', langmuirScheme: 'factor', mixedLayerRestratification: false };

const lowestLayer = (levels) => levels[levels.length - 1] - levels[levels.length - 2];
export const ROWS = [['cam26', PREVIOUS_LEVELS], ['bl48', sigmaInterfaces('bl48')]];
export const lowestLayerRatio = (levels) => lowestLayer(PREVIOUS_LEVELS) / lowestLayer(levels);
