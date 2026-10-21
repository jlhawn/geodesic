import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { withCadence } from '../js/cadence.module.js';
import { createSkin, SEA_RHO_CP, SKIN_DEFAULTS, absorbedAbove, surfaceShare } from '../js/physics/skin.module.js';
import { JERLOV } from '../js/ocean/layered.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';

/*
 * The skin driven as the survey's column drives its schemes: a 12 h
 * half-sine sun of daily mean S̄, a dark loss Q_ns = −(60 + 18U), the
 * air's u* = √(1.2e-3) U at 1.2 kg/m³, the foundation at 28 °C under a
 * 50 m slab, local time = UTC.
 */
const sun = (S, hour) => (hour > 6 && hour < 18 ? Math.PI * S * Math.sin(Math.PI * (hour - 6) / 12) : 0);
const ON = { warmLayer: true, coolSkin: true };

function column(U, S, { dt = 337.5, days = 4, depth = 50, options = ON, slab = null } = {}) {
  const skin = createSkin(1, options), perDay = Math.round(86400 / dt), ua = Math.sqrt(1.2e-3) * U;
  let capacity = depth * SEA_RHO_CP, T = 301.15;
  const rows = [];
  for (let n = 0; n < days * perDay; n++) {
    const hour = ((n + 0.5) * dt / 3600) % 24, solar = sun(S, hour);
    let nonSolar = -(60 + 18 * U);
    if (slab) {
      const surface = T + skin.offset(0, capacity);
      nonSolar -= slab.feedback * (surface - 301.15);
      const kept = solar * absorbedAbove(JERLOV.IB, depth), net = nonSolar + solar;
      T += dt * (net + slab.base - (solar - kept)) / capacity;
    }
    skin.step(0, nonSolar, solar, ua, 1.2, T, capacity, dt);
    rows.push({ day: Math.floor(n / perDay), hour: ((n + 1) * dt / 3600) % 24, offset: skin.offset(0, capacity), warm: skin.warm[0], cool: skin.cool[0], T, surface: T + skin.offset(0, capacity) });
  }
  return rows;
}
const ofDay = (rows, day) => rows.filter((r) => r.day === day);
const range = (rows, key = 'offset') => Math.max(...rows.map((r) => r[key])) - Math.min(...rows.map((r) => r[key]));
const at = (rows, hour) => rows.reduce((best, r) => (Math.abs(r.hour - hour) < Math.abs(best.hour - hour) ? r : best));

/*
 * The survey's tables C (Takaya with La_t = 0.3) and D (COARE 3.6) at
 * S̄ = 200, 250 and 300 W/m², and Gentemann et al. (2003)'s TMI fit at
 * Q_TOA = S̄/0.7 (survey §6a–6b).
 */
const WINDS = [1, 2, 4, 6, 8];
const TAKAYA = { 200: [1.84, 0.94, 0.21, 0.08, 0.04], 250: [2.45, 1.34, 0.32, 0.12, 0.06], 300: [3.07, 1.76, 0.45, 0.17, 0.09] };
const COARE = { 200: [2.02, 0.89, 0.27, 0.13, 0.06], 250: [3.06, 1.37, 0.42, 0.21, 0.11], 300: [4.20, 1.90, 0.60, 0.30, 0.17] };
const gentemannShape = (t) => { const w = 0.2668 * t; return 0.001 * (6.814 - 6.837 * Math.cos(w) - 8.427 * Math.sin(w) + 1.447 * Math.cos(2 * w) + 4.274 * Math.sin(2 * w) - 0.407 * Math.cos(3 * w) - 0.851 * Math.sin(3 * w) + 0.457 * Math.cos(4 * w) - 0.555 * Math.sin(4 * w) - 0.101 * Math.cos(5 * w) + 0.375 * Math.sin(5 * w)); };
const shapePeak = Math.max(...Array.from({ length: 96 }, (_, n) => gentemannShape(n / 4)));
const gentemann = (Q, U) => (Q < 132 ? 0 : shapePeak * ((Q - 132) - 9.632e-4 * (Q - 132) ** 2) * Math.exp(-0.53 * U));
const band = (S, w) => {
  const values = [TAKAYA[S][w], COARE[S][w], gentemann(S / 0.7, WINDS[w])];
  return [0.8 * Math.min(...values), 1.25 * Math.max(...values)];
};
const DESIGN_BAND = [[1.8, 3.8], [1.05, 1.7], [0.26, 0.56], [0.10, 0.26], [0.04, 0.14]];

test('the day-4 range of the atmosphere-facing skin lies within the surveyed schemes and the TMI fit at every wind', () => {
  const table = [];
  for (const S of [200, 250, 300]) {
    const row = WINDS.map((U, w) => {
      const r = range(ofDay(column(U, S), 3)), [lo, hi] = S === 250 ? DESIGN_BAND[w] : band(S, w);
      assert.ok(r >= lo && r <= hi, `S̄ ${S}, U ${U}: range ${r.toFixed(3)} K outside ${lo.toFixed(2)}–${hi.toFixed(2)}`);
      return r.toFixed(2);
    });
    table.push(`S̄ ${S}: ${row.join(' / ')}`);
  }
  console.log(`skin range (K) at U = ${WINDS.join(', ')} m/s — ${table.join('; ')}`);
});

test('the skin peaks in the early afternoon where the warm layer matters and at the latest by mid-afternoon at any wind, and reaches its minimum by dawn', () => {
  for (const U of WINDS) {
    const day = ofDay(column(U, 250), 3), top = day.reduce((a, r) => (r.offset > a.offset ? r : a)), bottom = Math.min(...day.map((r) => r.offset)), span = range(day);
    assert.ok(top.hour >= 12 && top.hour <= 16, `U ${U}: maximum at ${top.hour.toFixed(2)} h`);
    if (span >= 1) assert.ok(top.hour >= 13 && top.hour <= 16, `U ${U}: a warm layer of ${span.toFixed(2)} K peaks at ${top.hour.toFixed(2)} h`);
    assert.ok(at(day, 6).offset - bottom <= 0.01 * span, `U ${U}: 06 h is the day's minimum`);
    assert.ok(at(day, 9).offset - bottom > 0.01 * span, `U ${U}: warming has begun by 09 h`);
  }
});

test('the warm layer is gone by dawn and repeats from day to day', () => {
  for (const U of WINDS) {
    const rows = column(U, 250), dawn = at(ofDay(rows, 3), 6).warm;
    if (U >= 2) assert.ok(dawn < 0.05, `U ${U}: ΔT_w ${dawn} at 06 h`);
    const peak3 = Math.max(...ofDay(rows, 2).map((r) => r.warm)), peak4 = Math.max(...ofDay(rows, 3).map((r) => r.warm));
    assert.ok(Math.abs(peak4 - peak3) <= 0.02 * peak3, `U ${U}: day-4 peak ${peak4} against day 3's ${peak3}`);
  }
});

test('the implicit step gives the same day at N = 128\'s, N = 64\'s and N = 16\'s atmosphere steps', () => {
  for (const U of WINDS) {
    const ranges = [168.75, 337.5, 1350].map((dt) => range(ofDay(column(U, 250, { dt }), 3)));
    assert.ok(Math.max(...ranges) <= 1.1 * Math.min(...ranges), `U ${U}: ranges ${ranges.map((r) => r.toFixed(3)).join(', ')}`);
  }
});

test('the warm layer stops at its cap, never goes negative, and the cool skin stays within [0, 1] K', () => {
  const calm = column(0.2, 350);
  assert.ok(Math.max(...calm.map((r) => r.warm)) === SKIN_DEFAULTS.warmLayerCap, 'the cap binds at 0.2 m/s under 350 W/m²');
  for (const rows of [calm, ...WINDS.map((U) => column(U, 300))]) for (const r of rows) {
    assert.ok(r.warm >= 0 && r.warm <= SKIN_DEFAULTS.warmLayerCap, `ΔT_w ${r.warm}`);
    assert.ok(r.cool >= 0 && r.cool <= 1, `ΔT_c ${r.cool}`);
  }
});

test('the cool skin is 0.1–0.4 K at night for winds of 2–10 m/s', () => {
  const night = [2, 4, 6, 8, 10].map((U) => at(ofDay(column(U, 250), 3), 3).cool);
  for (const c of night) assert.ok(c >= 0.1 && c <= 0.4, `ΔT_c ${c}`);
  console.log(`cool skin at 03 h, U = 2, 4, 6, 8, 10 m/s: ${night.map((c) => c.toFixed(3)).join(' / ')} K`);
});

test('over a 20 m slab the skin adds its own range to the slab\'s', () => {
  for (const U of [1, 2, 4]) {
    const S = 250, slab = { feedback: 20, base: -(S * absorbedAbove(JERLOV.IB, 20) - (60 + 18 * U)) };
    const rows = column(U, S, { depth: 20, slab }), alone = range(ofDay(column(U, S, { depth: 20 }), 3));
    const day = ofDay(rows, 3), surface = range(day, 'surface'), slabRange = range(day, 'T');
    assert.ok(surface <= alone + slabRange + 0.05, `U ${U}: T_s range ${surface.toFixed(3)} against ${alone.toFixed(3)} + ${slabRange.toFixed(3)}`);
    assert.ok(slabRange >= 0.05 && slabRange <= 0.2, `U ${U}: the slab's own range ${slabRange.toFixed(3)} K`);
  }
  assert.ok(Math.abs(surfaceShare(5 * SEA_RHO_CP, SKIN_DEFAULTS) - 0.8615) < 1e-3 && Math.abs(surfaceShare(50 * SEA_RHO_CP, SKIN_DEFAULTS) - 0.9862) < 1e-3, 's = 1 − ν d/((ν+1) h₀)');
});

test('with both options off the skin adds nothing and steps nothing, and the model builds none', () => {
  const skin = createSkin(1, {});
  skin.step(0, -100, 800, 0.1, 1.2, 300, 50 * SEA_RHO_CP, 600);
  assert.equal(skin.offset(0, 50 * SEA_RHO_CP), 0);
  assert.equal(skin.warm[0], 0);
  assert.equal(skin.cool[0], 0);
  assert.equal(createModel(new Grid(16), {}).skin, null);
  assert.throws(() => createModel(new Grid(16), { surface: { exchange: 'fixed', warmLayer: true } }), /exchange 'roughness'/);
  assert.throws(() => createModel(new Grid(16), { surface: { dragCoefficient: 1.2e-3, coolSkin: true } }), /exchange 'roughness'/);
});

const N = 16, DT = 1350, BASE = { levels: sigmaInterfaces('bl36'), ...withCadence({ ocean: { penetration: 'IB' } }, DT, N) };
function started(surface) {
  const model = createModel(new Grid(N), { ...BASE, surface });
  initializeState(model).forEach((values, s) => model.state[s].set(values));
  model.ocean.initialize(model.state[3], model.state[6]);
  return model;
}

test('options absent and present but off give byte-equal states over twelve hours at N = 16', () => {
  const absent = started({}), off = started({ ...SKIN_DEFAULTS });
  for (let n = 0; n < 12 * 3600 / DT; n++) { absent.step(DT); off.step(DT); }
  absent.state.forEach((values, s) => assert.ok(Buffer.from(values.buffer).equals(Buffer.from(off.state[s].buffer)), `state ${s}`));
  const a = absent.ocean.serialize(), b = off.ocean.serialize();
  for (const name of Object.keys(a)) assert.ok(Buffer.from(Float64Array.from(a[name]).buffer).equals(Buffer.from(Float64Array.from(b[name]).buffer)), `ocean ${name}`);
});

test('the skin goes into a saved state and comes back from it, and a state without it starts the skin at 0', async () => {
  const model = started(ON);
  for (let n = 0; n < 16; n++) model.step(DT);
  const { warm, cool } = model.skin;
  assert.ok(warm.some((x) => x > 0) && cool.some((x) => x > 0), 'a lit sea has both');
  const saved = await decodeState(encodeState({ N, K: model.core.K, day: 0, time: model.time, ...model.skin.serialize() }, { f64: ['warmLayer', 'coolSkin'] }));
  const other = createSkin(warm.length, ON);
  other.load(saved);
  assert.deepEqual(Array.from(other.warm), Array.from(warm));
  assert.deepEqual(Array.from(other.cool), Array.from(cool));
  other.load({ N });
  assert.ok(other.warm.every((x) => x === 0) && other.cool.every((x) => x === 0));
});
