import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore, sigmaInterfaces, GRAVITY } from '../js/dynamics/sigmaCore.module.js';
import { createModel } from '../js/model.module.js';
import { createRadiation } from '../js/physics/radiation.module.js';
import { AEROSOL_SPECTRA } from '../js/physics/aerosolSpectra.module.js';
import { WAVELENGTHS, MIXED_DEPTH, MIXED_CAP, modeBand, streamOptics, aerosolTable, profileFractions } from '../js/aerosol.module.js';
import { BENCHMARK, modelColumn } from '../scripts/standardAtmospheres.mjs';
import { syntheticAerosol, fillAerosol } from './helpers/aerosol.mjs';

const levels = sigmaInterfaces('bl48'), mesh = buildMesh(new Grid(16)), C = mesh.nCells, band = (nm) => AEROSOL_SPECTRA.wavelengths.indexOf(nm);
const reference = (tau, x, g) => {
  const out = { visible: { depth: 0, scatter: 0, tilt: 0, absorbing: [0, 0], depths: [0, 0] }, rest: { depth: 0, scatter: 0, tilt: 0, absorbing: [0, 0], depths: [0, 0] } };
  for (let m = 0; m < 2; m++) {
    const set = m === 0 ? AEROSOL_SPECTRA.fine : AEROSOL_SPECTRA.coarse;
    AEROSOL_SPECTRA.wavelengths.forEach((nm, b) => {
      const t = set.r[b] * tau[m], xb = Math.min(1, Math.max(0, set.a[b] + set.c[b] * x[m])), gb = Math.min(0.99, Math.max(0, set.d[b] + set.e[b] * g[m]));
      for (const [kind, weights] of [['visible', AEROSOL_SPECTRA.visible], ['rest', AEROSOL_SPECTRA.nearInfrared]]) {
        const w = weights[b], s = out[kind];
        s.depth += w * t; s.scatter += w * t * (1 - xb); s.tilt += w * t * (1 - xb) * gb; s.absorbing[m] += w * t * xb; s.depths[m] += w * t;
      }
    });
  }
  return out;
};

test('a mode\'s optics at 550 nm are its own, at the others its depth times the wavelength\'s ratio and its albedo and asymmetry from the fitted lines, held within their limits', () => {
  const out = { tau: 0, x: 0, g: 0 };
  for (const mode of [0, 1]) {
    modeBand(mode, band(550), 0.37, 0.081, 0.66, out);
    assert.deepEqual([out.tau, out.x, out.g], [0.37, 0.081, 0.66]);
    const set = mode === 0 ? AEROSOL_SPECTRA.fine : AEROSOL_SPECTRA.coarse;
    for (let b = 0; b < WAVELENGTHS; b++) {
      modeBand(mode, b, 0.2, 0.05, 0.7, out);
      assert.ok(Math.abs(out.tau - 0.2 * set.r[b]) < 1e-15 && out.x >= 0 && out.x <= 1 && out.g >= 0 && out.g <= 0.99, `${AEROSOL_SPECTRA.wavelengths[b]} nm: ${JSON.stringify(out)}`);
      if (b < 9) assert.ok(Math.abs(out.x - Math.min(1, Math.max(0, set.a[b] + set.c[b] * 0.05))) < 1e-15 && Math.abs(out.g - Math.min(0.99, Math.max(0, set.d[b] + set.e[b] * 0.7))) < 1e-15);
    }
    modeBand(mode, band(1780), 0, 0, 0, out);
    assert.equal(out.tau, 0);
    modeBand(mode, band(400), 1, 1e3, -1e3, out);
    assert.deepEqual([out.x, out.g], [1, 0]);
    modeBand(mode, band(2320), 1, 0, 1e3, out);
    assert.equal(out.g, 0.99);
  }
  modeBand(0, band(1000), 1, 0, 0, out);
  assert.equal(out.x, 0, 'the fine mode\'s fitted 1 - SSA at 1 um is negative for a non-absorbing aerosol and is held at 0');
});

test('the streams\' mixture of the modes adds their depths, weighs albedo by depth and asymmetry by scattering, and gives the profile to the mode that absorbs', () => {
  const table = aerosolTable(), single = aerosolTable(), second = aerosolTable(), both = [[0.083, 0.059], [0.0703, 0.0268], [0.644, 0.76]];
  streamOptics(both[0], both[1], both[2], table);
  const wanted = reference(both[0], both[1], both[2]);
  assert.ok(Math.abs(table.tauV - wanted.visible.depth) < 1e-14 && Math.abs(table.tauN - wanted.rest.depth) < 1e-14);
  assert.ok(Math.abs(table.omegaV - wanted.visible.scatter / wanted.visible.depth) < 1e-14 && Math.abs(table.gV - wanted.visible.tilt / wanted.visible.scatter) < 1e-14);
  assert.ok(Math.abs(table.omegaN - wanted.rest.scatter / wanted.rest.depth) < 1e-14 && Math.abs(table.gN - wanted.rest.tilt / wanted.rest.scatter) < 1e-14);
  assert.ok(Math.abs(table.fineV - wanted.visible.absorbing[0] / (wanted.visible.absorbing[0] + wanted.visible.absorbing[1])) < 1e-14);
  assert.ok(Math.abs(table.fineN - wanted.rest.absorbing[0] / (wanted.rest.absorbing[0] + wanted.rest.absorbing[1])) < 1e-14);
  console.log(`global mean modes of August: visible stream tau ${table.tauV.toFixed(3)} omega ${table.omegaV.toFixed(3)} g ${table.gV.toFixed(3)} (the band data, weighted by depth: 0.152, 0.928, 0.701); rest of the beam ${table.tauN.toFixed(3)} ${table.omegaN.toFixed(3)} ${table.gN.toFixed(3)} (0.087, 0.940, 0.691; the albedo of mean absorption is higher than the depth-weighted mean's); the fine mode takes ${(100 * table.fineV).toFixed(0)} % of the visible absorption`);
  assert.ok(Math.abs(table.tauV - 0.152) < 0.008 && table.omegaV > 0.92 && table.omegaV < 0.95 && Math.abs(table.gV - 0.701) < 0.015, 'visible stream');
  assert.ok(Math.abs(table.tauN - 0.087) < 0.006 && table.omegaN > 0.93 && table.omegaN < 0.955 && Math.abs(table.gN - 0.691) < 0.015, 'rest of the beam');
  streamOptics([both[0][0], 0], [both[1][0], 0], [both[2][0], 0], single);
  streamOptics([0, both[0][1]], [0, both[1][1]], [0, both[2][1]], second);
  assert.ok(Math.abs(single.tauV + second.tauV - table.tauV) < 1e-15 && Math.abs(single.tauN + second.tauN - table.tauN) < 1e-15, 'the depths add');
  assert.equal(single.fineV, 1);
  assert.equal(second.fineV, 0);
  for (const key of ['omegaV', 'gV', 'omegaN', 'gN']) assert.ok(table[key] >= Math.min(single[key], second[key]) && table[key] <= Math.max(single[key], second[key]), `${key} lies between the modes'`);
  const doubled = aerosolTable();
  streamOptics(both[0].map((t) => 2 * t), both[1], both[2], doubled);
  for (const key of ['omegaV', 'gV', 'omegaN', 'gN', 'fineV', 'fineN']) assert.ok(Math.abs(doubled[key] - table[key]) < 1e-12, `${key} does not change with the depth`);
  assert.ok(Math.abs(doubled.tauV / table.tauV - 2) < 1e-12);
  const none = aerosolTable();
  streamOptics([0, 0], [0.05, 0.02], [0.7, 0.7], none);
  assert.deepEqual([none.tauV, none.omegaV, none.gV, none.fineV, none.tauN, none.omegaN, none.gN, none.fineN], [0, 1, 0, 0, 0, 1, 0, 0]);
  const clear = aerosolTable();
  streamOptics([0.2, 0.1], [0, 0], [0.7, 0.7], clear);
  assert.ok(clear.omegaV > 0.97 && clear.omegaV < 1, 'a non-absorbing visible aerosol is nearly conservative');
  const dark = aerosolTable();
  streamOptics([0.2, 0.1], [0.2, 0.2], [0.7, 0.7], dark);
  assert.ok(dark.omegaV < clear.omegaV && dark.omegaN < clear.omegaN, 'more absorption at 550 nm lowers the albedo of both streams');
});

test('a mode\'s shares of the layers sum to 1, put its mean height where it is asked to, hold the density constant to the mixed depth and thin out above it', () => {
  const fine = 3000, top = 60, dz = top / fine, heights = Float64Array.from({ length: fine }, (_, k) => top - (k + 0.5) * dz), shares = new Float64Array(fine);
  for (const [mean, depth] of [[0.3, 1], [0.45, 1.5], [0.7, 1.5], [1, 1.5], [2.06, 1.5], [1.7, 1], [4, 1.5], [6, 1]]) {
    profileFractions(mean, depth, heights, shares);
    assert.ok(Math.abs(shares.reduce((a, b) => a + b, 0) - 1) < 1e-12 && shares.every((f) => f >= 0), `shares at ${mean} km`);
    let centroid = 0;
    for (let k = 0; k < fine; k++) centroid += shares[k] * heights[k];
    assert.ok(Math.abs(centroid - mean) < 0.01 * mean + 0.5 * dz, `mean height ${centroid} against ${mean} km (mixed depth ${depth})`);
    const mixed = Math.min(depth, MIXED_CAP * mean), rate = (k) => shares[k] / dz, low = fine - 1, high = fine - 1 - Math.floor(0.5 * mixed / dz);
    assert.ok(Math.abs(rate(low) - rate(high)) < 1e-9 * rate(low), `constant density below ${mixed} km`);
    assert.ok(rate(fine - 1 - Math.floor(2 * mixed / dz)) < rate(low), 'and less above');
  }
  profileFractions(2.06, 1.5, heights, shares);
  let below = 0;
  for (let k = 0; k < fine; k++) if (heights[k] < 1) below += shares[k];
  const scale = 0.5 * ((2.06 - 1.5) + Math.sqrt(2.06 ** 2 + 2 * 2.06 * 1.5 - 1.5 ** 2));
  assert.ok(Math.abs(below - 1 / (1.5 + scale)) < 2e-3, `${below} of the aerosol below 1 km against ${1 / (1.5 + scale)}`);
  const model = [0.02, 0.1, 0.3, 0.7, 1.5, 3, 6, 12, 25].reverse(), onModel = new Float64Array(model.length);
  profileFractions(1.5, 1, model, onModel);
  assert.ok(Math.abs(onModel.reduce((a, b) => a + b, 0) - 1) < 1e-15 && onModel[0] > 0, 'the top layer holds what lies above the layers');
  profileFractions(-3, 1, model, onModel);
  profileFractions(60, 1, model, onModel);
  assert.ok(Math.abs(onModel.reduce((a, b) => a + b, 0) - 1) < 1e-15, 'a mean height beyond its range is held to it');
  assert.deepEqual(MIXED_DEPTH, [1.5, 1.0]);
});

function column(phis) {
  const core = createSigmaCore(mesh, { levels, surfaceGeopotential: phis === 0 ? null : new Float64Array(C).fill(phis) }), { K } = core.diagnostics, atmosphere = modelColumn(BENCHMARK.atmospheres.TROP, levels);
  const pi = new Float64Array(C).fill(atmosphere.ps), theta = new Float64Array(K * C), q = new Float64Array(K * C), qc = new Float64Array(K * C);
  core.diagnose(pi, theta);
  const { exnerLayer } = core.diagnostics;
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) { theta[k * C + i] = atmosphere.T[k] / exnerLayer[k * C + i]; q[k * C + i] = atmosphere.q[k]; }
  core.diagnose(pi, theta, q, qc);
  return { core, atmosphere, pi, theta, q, qc, K };
}

test('over 4.5 km of terrain the aerosol sits by its height above the ground, as it does over the sea', () => {
  const sea = column(0), plateau = column(4500 * GRAVITY), profiles = [];
  assert.ok(Math.abs(plateau.core.diagnostics.geopotential[(plateau.K - 1) * C] - sea.core.diagnostics.geopotential[(plateau.K - 1) * C] - 4500 * GRAVITY) < 1e-6, 'the plateau lifts the geopotential by the terrain');
  for (const { core, atmosphere, pi, theta, q, qc, K } of [sea, plateau]) {
    const fields = fillAerosol(new Float32Array(8 * C), [0.12, 0.08, 0.06, 0.02, 0.65, 0.76, 2.2, 1.4]);
    const radiation = createRadiation(mesh, core, { aerosol: 'mac', aerosolFields: fields, landAerosol: 0, seaAerosol: 0, clearSkyPass: true });
    radiation.setTime(0);
    const i = Array.from({ length: C }, (_, n) => n).reduce((best, n) => (radiation.cosZenith(n) > radiation.cosZenith(best) ? n : best), 0);
    radiation.column(i, atmosphere.ps, theta, atmosphere.Ts, 5, undefined, undefined, q[(K - 1) * C + i], q, qc, 0.1, 0.1, 0, 0, 0, 0, 0, 0);
    assert.ok(radiation.cosZenith(i) > 0.5 && radiation.budget.aerosolSolar > 0);
    profiles.push({ heights: Float64Array.from(radiation.aerosolProfiles.heights), fine: Float64Array.from(radiation.aerosolProfiles.fine), visible: Float64Array.from(radiation.aerosolProfiles.visible), aerosol: radiation.budget.aerosolSolar, heating: Float64Array.from(radiation.layerFlux) });
  }
  for (const key of ['fine', 'visible', 'heating']) for (let k = 0; k < sea.K; k++) assert.ok(Math.abs(profiles[0][key][k] - profiles[1][key][k]) < 1e-9 * Math.max(1, Math.abs(profiles[0][key][k])), `${key} of layer ${k}`);
  assert.ok(Math.abs(profiles[0].aerosol - profiles[1].aerosol) < 1e-9 * profiles[0].aerosol, 'and so does its absorption');
  const { heights, fine } = { heights: profiles[1].heights, fine: profiles[1].fine }, D = 1.5, scale = 0.5 * ((2.2 - D) + Math.sqrt(2.2 ** 2 + 2 * 2.2 * D - D ** 2));
  for (const j of [sea.K - 12, sea.K - 6, sea.K - 3, sea.K - 1]) {
    const height = 0.5 * (heights[j - 1] + heights[j]);
    let below = 0;
    for (let k = j; k < sea.K; k++) below += fine[k];
    const wanted = 1 - (height < D ? 1 - height / (D + scale) : scale / (D + scale) * Math.exp(-(height - D) / scale));
    assert.ok(Math.abs(below - wanted) < 1e-6, `the fine mode holds ${below} of its column below the interface ${height} km above the ground, the profile ${wanted}`);
  }
  assert.ok(heights[sea.K - 1] > 0 && heights[sea.K - 1] < 0.1 && heights[0] > 20, `layer heights above the ground from ${heights[sea.K - 1]} to ${heights[0]} km`);
});

test('aerosol takes \'flat\' or \'mac\', mac needs the CLIRAD gases and no cloud-top shielding, and aerosolScale a factor of 0 or more', () => {
  const { core } = column(0), make = (extra) => createRadiation(mesh, core, extra);
  assert.doesNotThrow(() => make({}));
  assert.doesNotThrow(() => make({ aerosol: 'flat', aerosolScale: 3 }));
  for (const extra of [{}, { shortwave: 'sequential' }, { nearInfrared: 'clirad', shortwave: 'sequential' }, { aerosolScale: 0 }, { aerosolScale: 2.5 }]) assert.doesNotThrow(() => make({ aerosol: 'mac', ...extra }), JSON.stringify(extra));
  for (const [bad, message] of [[{ aerosol: 'dust' }, /aerosol must be 'flat' or 'mac', not dust/], [{ aerosol: 'mac', solarGases: 'lacisHansen' }, /aerosol 'mac' needs solarGases 'clirad'/],
    [{ aerosol: 'mac', shortwave: 'sequential', vaporShielding: 'cloudTop' }, /aerosol 'mac' does not run with vaporShielding 'cloudTop'/], [{ aerosolScale: -1 }, /aerosolScale/], [{ aerosolScale: NaN }, /aerosolScale/], [{ aerosolScale: Infinity }, /aerosolScale/],
    [{ aerosol: 'mac', aerosolFields: new Float32Array(5) }, /aerosolFields must hold 8 fields/]]) assert.throws(() => make(bad), message, JSON.stringify(bad));
  const { bytes } = syntheticAerosol();
  assert.throws(() => createModel(new Grid(16), { ocean: false, radiation: { aerosol: 'mac', aerosolData: new Uint8Array(100) } }), /not an aerosol climatology|bytes, not/);
  const model = createModel(new Grid(16), { ocean: false, radiation: { aerosol: 'mac', aerosolData: bytes } });
  assert.equal(model.radiation.aerosolFields.length, 8 * C);
  assert.equal(model.aerosol.fields, model.radiation.aerosolFields);
  assert.equal(createModel(new Grid(16), { ocean: false }).aerosol, null);
});
