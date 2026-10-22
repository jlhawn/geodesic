import { test } from 'node:test';
import assert from 'node:assert/strict';
import { availableParallelism } from 'node:os';
import { Grid } from '../js/grid.module.js';
import { createParallelModel } from '../js/parallel.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { withCadence } from '../js/cadence.module.js';
import { transfer } from '../js/physics/exchange.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { SEA_RHO_CP, warmLayerStep } from '../js/physics/skin.module.js';
import { skinRecords, SEEN, STEP } from './helpers/skinProbe.mjs';

/*
 * The exchange guard: N = 16 over the layered ocean with Jerlov IB light,
 * warmLayer and coolSkin on, from initializeState at the March equinox,
 * with half the longitudes land below 69°, a polar cap of sea ice that
 * keeps its skin cleared, and the sea 6 K cooler than the start's between
 * 15° and 45°, so that warm air crosses cooler water by day. At each of
 * the physics' exchange calls over open water the cell is also asked at
 * T_slab, to find where the skin alone puts the
 * surface layer on its unstable branch, Ri_b(T_s) < 0 ≤ Ri_b(T_slab), and
 * each such cell-step must be real instability. A cell's branch turns
 * only where Ri_b crosses out of ±0.01, the near-neutral band in which the
 * two branches' coefficients differ by a few percent. The slab's daily
 * change is bounded from day 5, after the start's cooled water has
 * adjusted. SKIN_EXCHANGE_DAYS sets the length.
 */
const N = 16, DT = 1350, DAYS = Number(process.env.SKIN_EXCHANGE_DAYS ?? 12), NEUTRAL = 0.01, PER_DAY = 86400 / DT;

test(`over ${DAYS} days at N = 16 every cell-step the skin turns unstable is open, sunlit water below the cap, on the warm-layer equation, without a limit cycle or the Oct 6 runaway`, async (t) => {
  const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));
  const grid = new Grid(N), records = skinRecords(grid.size);
  const model = await createParallelModel(grid, { topography, levels: sigmaInterfaces('bl36'), ...withCadence({ ocean: { penetration: 'IB' } }, DT, N), surface: { warmLayer: true, coolSkin: true } }, Math.min(8, availableParallelism()), { module: new URL('./helpers/skinProbe.mjs', import.meta.url).href, data: records.data });
  t.after(() => model.close());
  initializeState(model, {}).forEach((values, s) => model.state[s].set(values));
  const concentration = new Float64Array(model.mesh.nCells);
  for (let i = 0; i < model.mesh.nCells; i++) if (!model.geography.land[i] && Math.abs(model.mesh.latCell[i]) > 1.2) { model.state[6][i] = 1.5; concentration[i] = 0.9; model.state[3][i] = 262; }
  for (let i = 0; i < model.mesh.nCells; i++) { const lat = Math.abs(model.mesh.latCell[i]) * 180 / Math.PI; if (!model.geography.land[i] && lat > 15 && lat < 45) model.state[3][i] -= 6; }
  model.seaIce.load(model.state[6], concentration);
  model.ocean.initialize(model.state[3], model.state[6]);
  const { exchange, radiation, skin, state, core, ocean } = model, C = model.mesh.nCells, K = core.K, bottom = (K - 1) * C;
  const { g, cp, dSigma } = core.diagnostics, cap = skin.options.warmLayerCap;

  const { seen, seenSet, stepped, steppedSet, flux } = records;

  const lastSun = new Float64Array(C).fill(-Infinity), branch = new Int8Array(C), slabBranch = new Int8Array(C), flips = new Int32Array(C), rawFlips = new Int32Array(C), dayStart = Float64Array.from(state[3]), icedToday = new Uint8Array(C);
  const worst = { dayCell: '', raw: 0, ratio: 0, sensible: 0, evaporation: 0, surfaceTheta: 0, theta: 0, alternation: 0, dayChange: 0, equation: 0 };
  let nightWarm = 0, warmestNight = 0, open = 0, sunlit = 0, flipped = 0, flippedSunlit = 0, capped = 0, iced = 0, warmest = 0;
  for (let n = 0; n < DAYS * PER_DAY; n++) {
    const lastSunBefore = Float64Array.from(lastSun), thetaBefore = Float64Array.from(state[1].subarray(bottom, bottom + C)), warmBefore = Float64Array.from(skin.warm), iceBefore = Float64Array.from(state[6]);
    steppedSet.fill(0);
    model.step(DT);
    for (const array of state) assert.ok(array.every(Number.isFinite), `finite state at step ${n}`);
    assert.ok(skin.warm.every(Number.isFinite) && skin.cool.every(Number.isFinite), `finite skin at step ${n}`);
    for (let i = 0; i < C; i++) {
      if (model.geography.land[i]) continue;
      if (radiation.surfaceSolar[i] > 0) lastSun[i] = model.time;
      if (iceBefore[i] > 0 || state[6][i] > 0) { iced++; icedToday[i] = 1; assert.equal(skin.warm[i], 0); assert.equal(skin.cool[i], 0); continue; }
      open++;
      const lit = radiation.surfaceSolar[i] > 0;
      if (lit) sunlit++;
      if (skin.warm[i] >= cap) capped++;
      warmest = Math.max(warmest, skin.warm[i]);
      const dark = model.time - lastSun[i] > 43200 + DT;
      if (dark && skin.warm[i] > 0) { nightWarm++; warmestNight = Math.max(warmestNight, skin.warm[i]); }
      assert.ok(steppedSet[i], `cell ${i} stepped`);
      const call = { before: stepped[STEP * i], args: Array.from(stepped.subarray(STEP * i + 1, STEP * (i + 1))) };
      const [nonSolar, solar, ustarAir, , slabT, capacity, dt] = call.args;
      assert.equal(nonSolar, flux[i] - radiation.surfaceSolar[i]);
      assert.equal(solar, radiation.surfaceSolar[i]);
      assert.equal(ustarAir, Math.sqrt(exchange.drag[i]) * exchange.wind[i]);
      assert.ok(Math.abs(slabT - state[3][i]) < 0.1, `cell ${i}: the skin stepped at the slab's temperature after its update`);
      assert.equal(dt, DT);
      const again = warmLayerStep(call.before, ...call.args, skin.options, skin.light);
      worst.equation = Math.max(worst.equation, Math.abs(again - skin.warm[i]));
      assert.ok(Math.abs(again - skin.warm[i]) <= 1e-9, `(d) cell ${i}: ${skin.warm[i]} against ${again}`);
      if (!seenSet[i]) continue;
      const s = { slab: seen[SEEN * i], ri: seen[SEEN * i + 1], z: seen[SEEN * i + 2], z0m: seen[SEEN * i + 3], z0h: seen[SEEN * i + 4] };
      const now = s.ri < -NEUTRAL ? -1 : s.ri > NEUTRAL ? 1 : branch[i], slabNow = s.slab < -NEUTRAL ? -1 : s.slab > NEUTRAL ? 1 : slabBranch[i], turned = branch[i] && now && now !== branch[i];
      rawFlips[i] = turned ? rawFlips[i] + 1 : 0;
      flips[i] = turned && slabNow === slabBranch[i] ? flips[i] + 1 : 0;
      branch[i] = now; slabBranch[i] = slabNow;
      worst.alternation = Math.max(worst.alternation, flips[i]); worst.raw = Math.max(worst.raw, rawFlips[i]);
      if (!(s.ri < 0 && s.slab >= 0)) continue;
      flipped++;
      if (lit) flippedSunlit++;
      assert.ok(capacity / SEA_RHO_CP >= 1 - 1e-9, `(a) cell ${i}: h₀ ${capacity / SEA_RHO_CP}`);
      assert.ok(!(warmBefore[i] > 0 && model.time - DT - lastSunBefore[i] > 43200 + DT), `(b) cell ${i}: ΔT_w ${warmBefore[i]} with no sun for 12 h`);
      assert.ok(warmBefore[i] < cap, `(c) cell ${i} at the cap`);
      const ratio = exchange.heat[i] / transfer(0, s.z, s.z0m, s.z0h).heat, surfaceTheta = Math.abs(DT * radiation.sensibleHeat[i] * g / (cp * state[0][i] * dSigma[K - 1]));
      worst.ratio = Math.max(worst.ratio, ratio); worst.sensible = Math.max(worst.sensible, Math.abs(radiation.sensibleHeat[i])); worst.evaporation = Math.max(worst.evaporation, radiation.evaporation[i] * 86400);
      worst.surfaceTheta = Math.max(worst.surfaceTheta, surfaceTheta); worst.theta = Math.max(worst.theta, Math.abs(state[1][bottom + i] - thetaBefore[i]));
      assert.ok(ratio <= 5, `(f) cell ${i}: C_H ${ratio.toFixed(2)}× neutral`);
      assert.ok(Math.abs(radiation.sensibleHeat[i]) <= 300, `(f) cell ${i}: sensible ${radiation.sensibleHeat[i]} W/m²`);
      assert.ok(radiation.evaporation[i] * 86400 <= 30, `(f) cell ${i}: evaporation ${radiation.evaporation[i] * 86400} mm/day`);
      assert.ok(surfaceTheta <= 1, `(f) cell ${i}: the surface flux moves the lowest layer's θ ${surfaceTheta} K in a step`);
    }
    if ((n + 1) % PER_DAY === 0) {
      if (n + 1 >= 5 * PER_DAY) for (let i = 0; i < C; i++) if (!model.geography.land[i] && !icedToday[i] && Math.abs(state[3][i] - dayStart[i]) > worst.dayChange) { worst.dayChange = Math.abs(state[3][i] - dayStart[i]); worst.dayCell = `cell ${i} at ${(model.mesh.latCell[i] * 180 / Math.PI).toFixed(0)}° on day ${(n + 1) / PER_DAY}`; }
      dayStart.set(state[3]); icedToday.fill(0);
    }
  }
  console.log(`skin exchange guard, ${DAYS} days: ${open} open-water cell-steps, ${flipped} flipped (${(100 * flippedSunlit / Math.max(1, sunlit)).toFixed(2)} % of sunlit ones, ${flipped - flippedSunlit} unlit), cap ${capped} (${(100 * capped / open).toFixed(4)} %), warmest ΔT_w ${warmest.toFixed(3)} K, ${iced} iced cell-steps held at 0, ${nightWarm} unflipped cell-steps warm after 12 h dark (warm air over cool water, at most ${warmestNight.toFixed(4)} K); flipped worst: C_H ${worst.ratio.toFixed(2)}× neutral, |sensible| ${worst.sensible.toFixed(1)} W/m², evaporation ${worst.evaporation.toFixed(2)} mm/day, surface θ step ${worst.surfaceTheta.toFixed(3)} K (whole step ${worst.theta.toFixed(3)} K); longest alternation the skin makes ${worst.alternation} steps (of the branch at T_s ${worst.raw}); largest daily slab change from day 5 ${worst.dayChange.toFixed(3)} K (${worst.dayCell}); equation ${worst.equation.toExponential(1)} K; ocean h₀ min ${Math.min(...ocean.capacity) / SEA_RHO_CP} m`);
  assert.ok(flipped > 0 && iced > 0, 'the run reaches the unstable branch and the ice');
  assert.ok(worst.alternation <= 2, `(e) a cell alternates branch on ${worst.alternation} consecutive steps`);
  assert.ok(capped < 1e-3 * open, `the cap binds in ${capped} of ${open} open-water cell-steps`);
  assert.ok(worst.dayChange <= 2, `the slab changes by ${worst.dayChange} K in a day`);
});
