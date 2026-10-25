import { test } from 'node:test';
import assert from 'node:assert/strict';
import { encodeForcing, decodeForcing, FORCING_FIELDS, forcingName, forcingDay } from '../js/forcing.module.js';
import { encodeState } from '../js/stateFile.module.js';

const N = 16, C = 10 * N * N + 2, E = 30 * N * N;
const header = { N, grid: { relax: 30 }, day: 1827, time: 1827 * 86400, seconds: 86400, steps: 64, oceanSteps: 16 };
const day = () => Object.fromEntries(FORCING_FIELDS.map(([name, where], f) => [name, Float32Array.from({ length: where === 'edges' ? E : C }, (_, i) => Math.sin(0.37 * i + f) * (f + 1) * 100)]));

test('a day of forcing comes back from its bytes as it went in', async () => {
  const fields = day();
  const back = await decodeForcing(encodeForcing(header, fields));
  assert.deepEqual({ ...back, fields: undefined }, { ...header, fields: undefined });
  for (const [name] of FORCING_FIELDS) assert.deepEqual(back.fields[name], fields[name], name);
});

test('a day from a drifting run carries the air stress on the edges, and a version-1 day without it still decodes', async () => {
  const airStress = Float32Array.from({ length: E }, (_, e) => Math.cos(0.11 * e));
  const back = await decodeForcing(encodeForcing(header, { ...day(), airStress }));
  assert.deepEqual(back.fields.airStress, airStress);
  assert.throws(() => encodeForcing(header, { ...day(), airStress: new Float32Array(C) }), /airStress has 2562 values, not the 7680 edges/);
  const older = await decodeForcing(encodeState({ kind: 'forcing', version: 1, ...header, ...day() }));
  assert.equal(older.fields.airStress, undefined);
  assert.deepEqual(older.fields.stress, day().stress);
  await assert.rejects(decodeForcing(encodeState({ kind: 'forcing', version: 3, ...header, ...day() })), /forcing version 3, not 1 or 2/);
});

test('forcing files are named by the day they end', () => {
  assert.equal(forcingName(7), 'forcing-0007.bin');
  assert.equal(forcingDay(forcingName(1827)), 1827);
  assert.equal(forcingDay('spin64_day1827.bin'), null);
});

test('a state file, a missing field or a field of the wrong size is refused', async () => {
  await assert.rejects(decodeForcing(encodeState({ N, day: 1, surfaceT: new Float32Array(C) })), /not a forcing file/);
  const { rain, ...partial } = day();
  assert.throws(() => encodeForcing(header, partial), /lacks rain/);
  assert.throws(() => encodeForcing(header, { ...day(), stress: new Float32Array(C) }), /stress has 2562 values/);
});
