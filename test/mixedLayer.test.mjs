import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createMixedLayer, dycomsLongwave } from '../js/physics/mixedLayer.module.js';
import { adiabaticWaterLapse } from '../js/physics/radiation.module.js';

const CONSTANTS = { cp: 1015, R: 287, latentHeat: 2.47e6 };
const G = 9.806;
const ZI = 840;
const INITIAL = { h: ZI, thetaL: 289.0, qt: 9.0e-3 };
const DT = 5;

/*
 * DYCOMS-II RF01 (Stevens et al. 2005): the layer, its free troposphere,
 * subsidence, prescribed surface fluxes and idealised longwave, with the
 * case's physical constants. The Boussinesq density is the initial
 * layer's mean density, held for the whole run.
 */
function rf01(mlm, overrides = {}) {
  const forcing = {
    surfacePressure: 101780, sensibleHeat: 15, evaporation: 115 / CONSTANTS.latentHeat,
    thetaLAbove: (z) => 297.5 + Math.cbrt(Math.max(0, z - ZI)), qtAbove: 1.5e-3,
    divergence: 3.75e-6, radiation: dycomsLongwave(), ...overrides,
  };
  forcing.density ??= mlm.diagnose(INITIAL, forcing).density;
  return forcing;
}

function run(mlm, forcing, onStep = null, hours = 4) {
  let state = { ...INITIAL };
  const hourly = [mlm.diagnose(state, forcing)];
  for (let n = 1; n <= hours * 3600 / DT; n++) {
    const d = mlm.diagnose(state, forcing);
    if (onStep) onStep(state, d);
    state = mlm.step(state, forcing, DT, d);
    if ((n * DT) % 3600 === 0) hourly.push(mlm.diagnose(state, forcing));
  }
  return { state, hourly, final: hourly[hourly.length - 1] };
}

const row = (d, t) => `${t} h: z_i ${d.h.toFixed(1)} m, LWP ${(1000 * d.liquidWaterPath).toFixed(1)} g/m², z_b ${d.cloudBase.toFixed(1)} m, w_e ${(1000 * d.entrainment).toFixed(2)} mm/s, A ${d.efficiency.toFixed(3)}, ΔF ${d.radiativeDivergence.toFixed(1)} W/m², Δθ_v ${d.virtualJump.toFixed(2)} K, BIR ${d.buoyancyIntegralRatio.toFixed(3)}, cover ${d.cover.toFixed(2)}`;

test('the RF01 initial layer has its cloud base within 10 m of 600 m, about 0.45 g/kg of cloud-top liquid, and the adiabatic water path of that thickness', () => {
  const mlm = createMixedLayer(CONSTANTS);
  const d = mlm.diagnose(INITIAL, rf01(mlm));
  const lapse = adiabaticWaterLapse(INITIAL.thetaL * d.cloudBaseExner, 1e5 * Math.pow(d.cloudBaseExner, CONSTANTS.cp / CONSTANTS.R), CONSTANTS.cp, CONSTANTS.R, G, CONSTANTS.latentHeat);
  const adiabatic = 0.5 * lapse * (ZI - d.cloudBase) ** 2;
  console.log(`RF01 at t = 0: cloud base ${d.cloudBase.toFixed(1)} m, cloud-top liquid ${(1000 * d.topLiquid).toFixed(3)} g/kg, LWP ${(1000 * d.liquidWaterPath).toFixed(1)} g/m² against ½ Γ_l Δz² = ${(1000 * adiabatic).toFixed(1)}, mean density ${d.density.toFixed(3)} kg/m³`);
  assert.ok(Math.abs(d.cloudBase - 600) < 10, `cloud base ${d.cloudBase}`);
  assert.ok(d.topLiquid > 0.43e-3 && d.topLiquid < 0.49e-3, `cloud-top liquid ${d.topLiquid}`);
  assert.ok(Math.abs(d.liquidWaterPath - adiabatic) < 0.05 * adiabatic, `LWP ${d.liquidWaterPath} against ${adiabatic}`);
  assert.ok(d.cover === 1 && d.buoyancyIntegralRatio < 0.15);
});

test('DYCOMS-II RF01 after 4 h: inversion 800–900 m, LWP 40–110 g/m², cloud base 500–700 m, entrainment 2–6 mm/s, under either closure', () => {
  for (const closure of ['radiative', 'buoyancy']) {
    const mlm = createMixedLayer({ ...CONSTANTS, closure });
    const { hourly, final } = run(mlm, rf01(mlm));
    console.log(`RF01, ${closure} closure (a_1 0.2, a_2 25):\n  ${hourly.map(row).join('\n  ')}`);
    assert.ok(final.h > 800 && final.h < 900, `${closure}: z_i ${final.h}`);
    assert.ok(final.liquidWaterPath > 0.040 && final.liquidWaterPath < 0.110, `${closure}: LWP ${final.liquidWaterPath}`);
    assert.ok(final.cloudBase > 500 && final.cloudBase < 700, `${closure}: z_b ${final.cloudBase}`);
    assert.ok(final.entrainment > 2e-3 && final.entrainment < 6e-3, `${closure}: w_e ${final.entrainment}`);
    for (const d of hourly) assert.ok(d.cover === 1 && d.buoyancyIntegralRatio < 0.15, `${closure}: the layer stays coupled`);
  }
  const plain = createMixedLayer({ ...CONSTANTS, evaporativeEnhancement: 0 });
  console.log(`without the evaporative enhancement (A = 0.2): ${row(run(plain, rf01(plain)).final, 4)}`);
});

test('without radiative cooling the layer stops entraining: it sinks with the subsidence and warms and moistens by the surface fluxes alone', () => {
  const mlm = createMixedLayer(CONSTANTS);
  const forcing = rf01(mlm, { radiation: () => 0 });
  let heat = INITIAL.thetaL, water = INITIAL.qt, h = INITIAL.h, steps = 0;
  const { state, final } = run(mlm, forcing, (s, d) => {
    assert.equal(d.entrainment, 0);
    assert.equal(d.radiativeDivergence, 0);
    h *= 1 - forcing.divergence * DT;
    heat += DT * forcing.sensibleHeat / (forcing.density * CONSTANTS.cp * h);
    water += DT * forcing.evaporation / (forcing.density * h);
    steps++;
  });
  console.log(`no longwave: after 4 h z_i ${state.h.toFixed(1)} m, θ_l ${state.thetaL.toFixed(3)} K, q_t ${(1000 * state.qt).toFixed(3)} g/kg, LWP ${(1000 * final.liquidWaterPath).toFixed(1)} g/m²`);
  assert.equal(steps, 4 * 3600 / DT);
  assert.ok(Math.abs(state.h - h) < 1e-9 * h, `z_i ${state.h} against ${h}`);
  assert.ok(Math.abs(state.thetaL - heat) < 1e-12 * heat, `θ_l ${state.thetaL} against ${heat}`);
  assert.ok(Math.abs(state.qt - water) < 1e-12 * water, `q_t ${state.qt} against ${water}`);
  assert.ok(state.h < INITIAL.h && state.thetaL > INITIAL.thetaL && state.qt > INITIAL.qt);
});

test('doubling the free-tropospheric moisture jump from 3.75 to 7.5 g/kg thins the cloud by entraining drier air', () => {
  const mlm = createMixedLayer(CONSTANTS);
  const drying = (qtAbove) => {
    let total = 0;
    const result = run(mlm, rf01(mlm, { qtAbove }), (s, d) => { total += DT * d.entrainment * -d.qtJump; });
    return { ...result, drying: total };
  };
  const control = drying(1.5e-3), moist = drying(5.25e-3);
  console.log(`after 4 h, Δq_t −7.5 g/kg: ${row(control.final, 4)}; entrained drying ${(1e3 * control.drying).toFixed(3)} m·g/kg\n            Δq_t −3.75 g/kg: ${row(moist.final, 4)}; entrained drying ${(1e3 * moist.drying).toFixed(3)} m·g/kg`);
  assert.ok(control.final.liquidWaterPath < moist.final.liquidWaterPath - 0.010, `LWP ${control.final.liquidWaterPath} against ${moist.final.liquidWaterPath}`);
  assert.ok(control.final.cloudBase > moist.final.cloudBase, `cloud base ${control.final.cloudBase} against ${moist.final.cloudBase}`);
  assert.ok(control.drying > 1.5 * moist.drying, `entrainment drying ${control.drying} against ${moist.drying}`);
});

test('over 4 h of RF01 with drizzle the column\'s water and heat change by exactly the entrainment, surface, radiative, drizzle and subsidence fluxes', () => {
  const mlm = createMixedLayer({ ...CONSTANTS, drizzle: true, dropletNumber: 100 });
  const forcing = rf01(mlm);
  const radiation = dycomsLongwave();
  const rho = forcing.density, cp = CONSTANTS.cp, L = CONSTANTS.latentHeat;
  const water = { entrainment: 0, subsidence: 0, surface: 0, drizzle: 0, gross: 0 };
  const heat = { entrainment: 0, subsidence: 0, surface: 0, radiation: 0, drizzle: 0, gross: 0 };
  const { state } = run(mlm, forcing, (s, d) => {
    assert.ok(Math.abs(d.radiativeDivergence - (radiation(d.liquidWaterPath, 0) - radiation(0, d.liquidWaterPath))) < 1e-12);
    assert.ok(Math.abs(d.entrainment - d.efficiency * d.radiativeDivergence / (rho * cp * d.virtualJump)) < 1e-15);
    const w = [rho * d.entrainment * d.qtAbove, rho * d.subsidence * s.qt, d.evaporation, -d.drizzle];
    const e = [rho * cp * d.entrainment * d.thetaLAbove, rho * cp * d.subsidence * s.thetaL, d.sensibleHeat, -d.radiativeDivergence, L * d.drizzle / d.cloudBaseExner];
    ['entrainment', 'subsidence', 'surface', 'drizzle'].forEach((name, i) => { water[name] += DT * w[i]; water.gross += DT * Math.abs(w[i]); });
    ['entrainment', 'subsidence', 'surface', 'radiation', 'drizzle'].forEach((name, i) => { heat[name] += DT * e[i]; heat.gross += DT * Math.abs(e[i]); });
  });
  const waterChange = rho * (state.h * state.qt - INITIAL.h * INITIAL.qt);
  const heatChange = rho * cp * (state.h * state.thetaL - INITIAL.h * INITIAL.thetaL);
  const waterFlux = water.entrainment + water.subsidence + water.surface + water.drizzle;
  const heatFlux = heat.entrainment + heat.subsidence + heat.surface + heat.radiation + heat.drizzle;
  console.log(`4 h water budget (kg/m²): change ${waterChange.toFixed(6)} = entrainment ${water.entrainment.toFixed(4)} + subsidence ${water.subsidence.toFixed(4)} + evaporation ${water.surface.toFixed(4)} + drizzle ${water.drizzle.toFixed(4)}, residual ${((waterChange - waterFlux) / water.gross).toExponential(1)} of the gross flux`);
  console.log(`4 h heat budget (MJ/m²): change ${(heatChange / 1e6).toFixed(6)} = entrainment ${(heat.entrainment / 1e6).toFixed(3)} + subsidence ${(heat.subsidence / 1e6).toFixed(3)} + sensible ${(heat.surface / 1e6).toFixed(4)} + longwave ${(heat.radiation / 1e6).toFixed(4)} + drizzle ${(heat.drizzle / 1e6).toFixed(5)}, residual ${((heatChange - heatFlux) / heat.gross).toExponential(1)} of the gross flux`);
  assert.ok(water.drizzle < 0 && heat.drizzle > 0, 'the drizzle takes water and leaves its latent heat');
  assert.ok(Math.abs(heat.surface - 15 * 4 * 3600) < 1e-6 && Math.abs(water.surface - 115 / L * 4 * 3600) < 1e-12);
  assert.ok(Math.abs(waterChange - waterFlux) < 1e-9 * water.gross, `water residual ${waterChange - waterFlux}`);
  assert.ok(Math.abs(heatChange - heatFlux) < 1e-9 * heat.gross, `heat residual ${heatChange - heatFlux}`);
});

test('the cover is 1 while the layer is coupled and falls to the trade-cumulus 0.3 as the buoyancy integral ratio passes from 0.15 to 0.4; a clear or uncapped layer has none', () => {
  const mlm = createMixedLayer(CONSTANTS);
  const cover = (sensibleHeat, latentHeat) => mlm.diagnose(INITIAL, rf01(mlm, { sensibleHeat, evaporation: latentHeat / CONSTANTS.latentHeat }));
  const coupled = cover(15, 115), partial = cover(-10, 60), decoupled = cover(-20, 30);
  console.log(`BIR and cover of the RF01 layer: surface fluxes 15/115 W/m² ${coupled.buoyancyIntegralRatio.toFixed(3)} → ${coupled.cover.toFixed(2)}, −10/60 ${partial.buoyancyIntegralRatio.toFixed(3)} → ${partial.cover.toFixed(2)}, −20/30 ${decoupled.buoyancyIntegralRatio.toFixed(3)} → ${decoupled.cover.toFixed(2)}`);
  assert.ok(coupled.buoyancyIntegralRatio < 0.15 && coupled.cover === 1);
  assert.ok(partial.buoyancyIntegralRatio > 0.15 && partial.buoyancyIntegralRatio < 0.4);
  assert.ok(Math.abs(partial.cover - (1 - 0.7 * (partial.buoyancyIntegralRatio - 0.15) / 0.25)) < 1e-12, `cover ${partial.cover}`);
  assert.ok(decoupled.buoyancyIntegralRatio > 0.4 && decoupled.cover === 0.3);
  const clear = mlm.diagnose({ ...INITIAL, h: 500 }, rf01(mlm));
  assert.ok(clear.liquidWaterPath === 0 && clear.cover === 0 && clear.entrainment === 0 && clear.cloudBase === 500);
  const uncapped = mlm.diagnose(INITIAL, rf01(mlm, { thetaLAbove: 289 }));
  assert.ok(uncapped.virtualJump < 0 && uncapped.cover === 0 && uncapped.entrainment === 0);
});

test('with the sea-surface temperature instead of fluxes the bulk formulas give about the case\'s 15 W/m² sensible and 115 W/m² latent heat', () => {
  const mlm = createMixedLayer(CONSTANTS);
  const d = mlm.diagnose(INITIAL, rf01(mlm, { sensibleHeat: undefined, evaporation: undefined, seaSurfaceTemperature: 292.5, transferVelocity: 0.0011 * Math.hypot(6, 4.25) }));
  const latent = CONSTANTS.latentHeat * d.evaporation;
  console.log(`SST 292.5 K, C_T = 0.0011, 7.35 m/s: sensible ${d.sensibleHeat.toFixed(1)} W/m², latent ${latent.toFixed(1)} W/m²`);
  assert.ok(d.sensibleHeat > 10 && d.sensibleHeat < 25, `sensible ${d.sensibleHeat}`);
  assert.ok(latent > 100 && latent < 130, `latent ${latent}`);
  assert.throws(() => createMixedLayer({ closure: 'lilly' }));
});
