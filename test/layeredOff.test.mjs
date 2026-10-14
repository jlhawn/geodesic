import { test } from 'node:test';
import assert from 'node:assert/strict';
import { oceanShaderDigest, stateDigest, SEA_CONFIGS } from './helpers/offIdentity.mjs';

const PACKAGE = { everySteps: 1, penetration: 'IB', carryMomentum: true, momentumDepth: 50, refillRate: 1 / 3600, refillFloor: 5, minimumThickness: 1, shallowestMixedDepth: 1, shearEntrainment: 0.65, ekmanDepth: 0.4, momentumFloor: 10, stressThickness: 20, retreatStirring: 2.5, langmuir: 2, retreatEkman: 0.7 };
const SHADERS = {
  layered: {},
  packaged: PACKAGE,
  daily: { ...PACKAGE, retreatDaily: true, retreatLight: 1 },
  rings: { ...PACKAGE, closureFill: 1 },
  shearMix: { ...PACKAGE, interiorShearMixing: true },
};
const SHADER_PINS = { layered: '42:acc5dab238e1ece3', packaged: '42:8a1be7753ab9db69', daily: '42:f5390582335fd054', rings: '42:599097d39e56fecd', shearMix: '42:ae5b07c9a880c062' };
const STATE_PINS = ['6c8721dcf7cb8b3d', 'c09cf1bf7adb3fa9', '73ebe79a3a9f7785'];
const OFF = { ekmanBlock: false, upwellingSource: 0, blockDrag: false };

test('with the Ekman block, the upwelling source and the block drag absent or off, the ocean shaders are the parent\'s, source for source', async () => {
  for (const [name, options] of Object.entries(SHADERS)) {
    const absent = await oceanShaderDigest(options), off = await oceanShaderDigest({ ...options, ...OFF });
    assert.equal(absent, SHADER_PINS[name], `${name}: the default ocean shaders changed; if that is the commit's purpose, set the pin to ${absent}`);
    assert.equal(off, absent, `${name}: the options present but off compile other shaders than absent`);
  }
});

test('with the options absent, the CPU ocean after 36 hours on a sea with a continent is the parent\'s, bit for bit', async () => {
  for (let which = 0; which < SEA_CONFIGS; which++) {
    const absent = await stateDigest(which);
    assert.equal(absent, STATE_PINS[which], `configuration ${which}: the default ocean changed; if that is the commit's purpose, set the pin to ${absent}`);
  }
});
