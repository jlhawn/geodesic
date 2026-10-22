import { test } from 'node:test';
import assert from 'node:assert/strict';
import { modelDigest, OVERCAST, UNSCATTERED, GRAY_GASES, GREY_ICE, RAYLEIGH_TOP } from './helpers/modelDigest.mjs';

test('with mixedLayerDeck: false and the purely scattering clouds of cloudSolarAbsorption: 0, cloudScattering: 55 the model is bit-identical to the engine before the mixed-layer deck; by default the deck follows the mixed-layer model', () => {
  const before = 'afbccef040c18dad3d4a15487681dd04';
  assert.equal(modelDigest({ mixedLayerDeck: false, cloudSolarAbsorption: 0, cloudScattering: 55, cloudAbsorption: 130, ...OVERCAST, ...UNSCATTERED, ...GRAY_GASES }).digest, before);
  assert.equal(modelDigest({ mixedLayerDeck: false, cloudSolarAbsorption: 0, cloudScattering: 55, cloudAbsorption: 130, ...OVERCAST, ...UNSCATTERED, ...GRAY_GASES }, {}, {}, RAYLEIGH_TOP).digest, '3864f85972740bb3640b3a3c111ff01e', 'with the Rayleigh top and neither sponge nor gravity waves, the engine before the model top');
  assert.equal(modelDigest({ mixedLayerDeck: false, cloudSolarAbsorption: 0, cloudScattering: 55, cloudAbsorption: 130, ...OVERCAST, ...UNSCATTERED, ...GRAY_GASES }, {}, GREY_ICE, RAYLEIGH_TOP).digest, 'bb924e00d741ac658132742934d24a61', 'on grey ice with unaged snow and the Rayleigh top, the engine before the ice and snow albedo depended on temperature and age');
  assert.equal(modelDigest({ mixedLayerDeck: false, cloudSolarAbsorption: 0, cloudScattering: 55, cloudAbsorption: 130, ...OVERCAST, ...UNSCATTERED, ...GRAY_GASES }, {}, GREY_ICE).digest, 'cf3cefdfcbd177c312e275780a78002a', 'on grey ice with unaged snow, the model-top parent\'s digest');
  assert.equal(modelDigest({ mixedLayerDeck: false, cloudSolarAbsorption: 0, cloudScattering: 55, cloudAbsorption: 130, ...OVERCAST, ...UNSCATTERED, longwaveOverlap: 'random' }, {}, GREY_ICE, RAYLEIGH_TOP).digest, 'ba8fb08dab8472ebc790cbd18d958c6f', 'on grey ice under the spectral gases with the longwave\'s random overlap and the Rayleigh top, the gas parent\'s digest before the model top');
  assert.equal(modelDigest({ mixedLayerDeck: false, cloudSolarAbsorption: 0, cloudScattering: 55, cloudAbsorption: 130, ...OVERCAST, ...UNSCATTERED }, {}, GREY_ICE, RAYLEIGH_TOP).digest, '5163aa93e5ee2977c1f3de5867984e0f', 'on grey ice under the spectral gases with the longwave\'s exponential-random overlap and the Rayleigh top');
  assert.notEqual(modelDigest({ mixedLayerDeck: false }).digest, before, 'by default cloud water absorbs sunlight');
  const fresh = modelDigest();
  assert.equal(fresh.digest, modelDigest({ mixedLayerDeck: true }).digest);
  assert.notEqual(fresh.digest, before);
  assert.ok(fresh.radiation.stratusFraction.every((f) => f === 0) && fresh.radiation.mlmCover.every((f) => f === 0), 'twelve steps from rest build no 2 K inversion, and no empirical deck stands in');
  const shadow = modelDigest({ mixedLayerDeck: true, stratusSubsidence: 0, minimumInversion: 0 });
  const { mlmCover, mlmWater, mlmEntrainment, stratusFraction } = shadow.radiation;
  let covered = 0;
  for (let i = 0; i < mlmCover.length; i++) {
    assert.ok(mlmCover[i] >= 0 && mlmCover[i] <= 1 && mlmWater[i] >= 0 && mlmEntrainment[i] >= 0 && Number.isFinite(mlmWater[i] + mlmEntrainment[i]));
    assert.ok(stratusFraction[i] <= mlmCover[i]);
    if (stratusFraction[i] > 0) covered++;
  }
  assert.ok(covered > 0, 'some cells carry a mixed-layer deck');
});
