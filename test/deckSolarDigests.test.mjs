import { test } from 'node:test';
import assert from 'node:assert/strict';
import { DAY } from '../js/physics/radiation.module.js';
import { modelDigest, OVERCAST, UNSCATTERED, GRAY_GASES, GREY_ICE, REDIAGNOSED } from './helpers/modelDigest.mjs';

test('with the purely scattering clouds of cloudSolarAbsorption: 0, cloudScattering: 55 the engine is bit-identical to the deck before it absorbed sunlight whether or not stratusSolar is on, while absorbing clouds, the carried height and the gate\'s memory each change the digest', () => {
  const forced = { stratusSubsidence: 0, minimumInversion: 0, subsidenceSmoothing: 0, subsidenceMemory: 10 * DAY }, scatteringOnly = { cloudSolarAbsorption: 0, cloudScattering: 55, cloudAbsorption: 130, ...OVERCAST, ...UNSCATTERED, ...GRAY_GASES };
  assert.equal(modelDigest({ stratusSolar: false, ...scatteringOnly }).digest, '03cebf2835d18c5346aa6f8ff7575c47');
  assert.equal(modelDigest({ stratusSolar: false, ...scatteringOnly }, {}, GREY_ICE).digest, 'b9e7763656da872aed702e3b2d50681b');
  assert.equal(modelDigest({ ...forced, ...REDIAGNOSED, deckSlab: 'midpoint', deckReference: 'layer', ...scatteringOnly }, {}, GREY_ICE).digest, '9e4cf3aff63eddaa8f6065d17771b2ac');
  assert.equal(modelDigest({ ...forced, ...REDIAGNOSED, deckSlab: 'midpoint', deckReference: 'layer', stratusSolar: false, ...scatteringOnly }).digest, '1c0726f9f9f190a8f99f1248afba0b43');
  assert.equal(modelDigest({ ...forced, ...REDIAGNOSED, deckSlab: 'midpoint', deckReference: 'layer', ...scatteringOnly }).digest, '1c0726f9f9f190a8f99f1248afba0b43');
  assert.notEqual(modelDigest({ ...forced, ...REDIAGNOSED }).digest, '1c0726f9f9f190a8f99f1248afba0b43');
  assert.notEqual(modelDigest({ ...forced, ...scatteringOnly }).digest, '1c0726f9f9f190a8f99f1248afba0b43', 'the carried height and the gate\'s memory change the deck');
});
