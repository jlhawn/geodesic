import { test } from 'node:test';
import assert from 'node:assert/strict';
import { MOIST_DEFAULTS } from '../js/physics/moist.module.js';
import { modelDigest, GRAY_GASES, RAYLEIGH_TOP } from './helpers/modelDigest.mjs';

test('the uniform condensation, saturation over ice and falling ice of the moist defaults give their own twelve-step digest, which the overlap moves (the deck digests run with all three off and are the parent\'s), and the cumulus memory its own again', () => {
  const defaults = { condensation: MOIST_DEFAULTS.condensation, iceSaturation: MOIST_DEFAULTS.iceSaturation, iceFall: MOIST_DEFAULTS.iceFall, cumulusMemory: 0 };
  const now = modelDigest({}, defaults).digest, random = modelDigest({ cloudOverlap: 'maximumRandom' }, defaults).digest;
  const remembered = { ...defaults, cumulusMemory: MOIST_DEFAULTS.cumulusMemory }, memory = modelDigest({}, remembered).digest, memoryRandom = modelDigest({ cloudOverlap: 'maximumRandom' }, remembered).digest;
  console.log(`12 steps at N=16: ${now} under the defaults without the cumulus memory, ${random} with maximum-random overlap; ${memory} and ${memoryRandom} with it`);
  assert.equal(now, 'aeeea325273097cc8876b2210863e26a');
  assert.equal(memory, '54e47e59d5a21845e79d0df05a2f7855');
  assert.equal(memoryRandom, '052577499e36a13b4f0b0b8ce58bec68');
  assert.equal(modelDigest({ longwaveOverlap: 'random' }, defaults).digest, 'a84c7875159a514907ab857d19bcf4ec', 'with the longwave\'s random overlap, the integration parent\'s digest');
  assert.equal(modelDigest({ longwaveOverlap: 'random' }, defaults, {}, RAYLEIGH_TOP).digest, 'a1265117281b3c22f1df4b021bb61789', 'with the longwave\'s random overlap and the Rayleigh top, the engine before the longwave overlapped as the shortwave does');
  assert.equal(modelDigest({}, defaults, {}, RAYLEIGH_TOP).digest, '89b063a335bb006d3b99854d7726b111', 'with the Rayleigh top, the convection parent\'s digest');
  assert.equal(modelDigest({ ...GRAY_GASES, visibleFraction: 0.5, rayleighBands: [[0.712, 0.0874], [0.288, 0.5687]], nearInfraredRayleigh: 0 }, defaults, {}, RAYLEIGH_TOP).digest, '8050c05f0cfc8a8fc3b108045d0734d8', 'under the gray gases, the visible split and Rayleigh bands before the gas branch\'s and the Rayleigh top, the cloud parent\'s digest');
  assert.notEqual(random, now);
  assert.notEqual(modelDigest({ cloudOverlap: 'maximumRandom' }).digest, random);
});

test('the fog scheme of the moist defaults leaves the twelve-step digest as it is: no lowest-layer cloud forms in twelve steps from rest at N=16', () => {
  const fog = modelDigest({}, { fogDroplets: MOIST_DEFAULTS.fogDroplets, fogDeposition: MOIST_DEFAULTS.fogDeposition });
  assert.equal(fog.digest, 'e1c1475216c90ea3ac50fa8052bc7459');
  assert.equal(fog.digest, modelDigest({}).digest);
});

test('the Doppler-core lines of longwaveLines \'voigt\', the default, give their own twelve-step digest (the other digests run with the scaled lines and are the parent\'s)', () => {
  assert.equal(modelDigest({ longwaveLines: 'voigt' }).digest, 'a0edcffa8b8b2a1c874b497ac066e0aa');
});
