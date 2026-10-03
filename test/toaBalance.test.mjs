import { test } from 'node:test';
import assert from 'node:assert/strict';
import { writeFileSync, mkdtempSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { readDays, window, summary, report } from '../scripts/toaBalance.mjs';

const dayLine = (day, asr, olr, extra = '') => `day ${day} (0.1 min): Ts ${(15 + day / 1000).toFixed(2)} °C, ASR ${asr.toFixed(1)} (atmosphere 80.0) OLR ${olr.toFixed(1)} W/m², ps 540–1030 hPa, max wind 80.0 m/s, precip 2.50 mm/d, ice 2.5% (N 1.0 S 1.5 Mkm²), albedo 0.300, SWCRE -50.0 LWCRE 25.0${extra}`;

function writeLog(dir, name, lines) {
  const file = join(dir, name);
  writeFileSync(file, lines.join('\n') + '\n');
  return file;
}

test('the windows are the last 1, 7, 30 and 365 logged days and a whole year is 365 of them', () => {
  const dir = mkdtempSync(join(tmpdir(), 'toa-'));
  const lines = ["ASR, atmosphere, OLR and albedo below are day means over the day's 256 steps"];
  for (let d = 1; d <= 400; d++) lines.push(dayLine(d, 240 + 10 * Math.sin((2 * Math.PI * d) / 365), 236));
  const days = readDays([writeLog(dir, 'a.log', lines)]);
  assert.equal(days.length, 400);
  assert.equal(days.dayMeans, true);
  assert.equal(window(days, 400, 7).length, 7);
  assert.equal(window(days, 400, 7)[0].day, 394);
  const year = summary(window(days, 365, 365));
  assert.equal(year.n, 365);
  assert.ok(Math.abs(year.net - 4) < 1e-6, `a whole year's net is the constant part, got ${year.net}`);
  const out = report(days);
  assert.match(out, /^last day \(to day 400\)/m);
  assert.match(out, /^year 1 \(days 1–365\) +net \+4\.00 W\/m²/m);
  assert.doesNotMatch(out, /aliased/);
  assert.match(out, /last 365 days \(to day 400\)[^\n]*\n/);
  assert.doesNotMatch(out, /year 2/);
});

test('a day logged again takes its last line, a missing day is counted, and a sampled log is flagged', () => {
  const dir = mkdtempSync(join(tmpdir(), 'toa-'));
  const a = writeLog(dir, 'a.log', ["ASR, atmosphere, OLR and albedo below are day means over the day's 256 steps", dayLine(1, 240, 230), dayLine(2, 240, 230), dayLine(3, 999, 230)]);
  const b = writeLog(dir, 'b.log', [dayLine(3, 240, 230), dayLine(5, 240, 230)]);
  const days = readDays([a, b]);
  assert.deepEqual(days.map((d) => d.day), [1, 2, 3, 5]);
  assert.equal(days[2].asr, 240);
  assert.equal(window(days, 5, 7).length, 4);
  assert.match(report(days), /last 7 days \(to day 5\)[^\n]*\[4 of 7 days\]/);
  const old = readDays([writeLog(dir, 'c.log', ['day 1 (0.1 min): Ts 15.00 °C, ASR 240.0 OLR 230.0 W/m², ps 540–1030 hPa, max wind 80.0 m/s, precip 2.50 mm/d, ice 2.5% (N 1.0 S 1.5 Mkm²), albedo 0.300'])]);
  assert.equal(old.length, 1);
  assert.ok(Number.isNaN(old[0].atmosphere));
  assert.match(report(old), /^the daily lines are samples at one time of day/);
});
