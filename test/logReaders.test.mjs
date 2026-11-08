import { test } from 'node:test';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdtempSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = fileURLToPath(new URL('..', import.meta.url));

const dayLine = (day, { fog = true, label = false, cloud = true, current = 1.2, clamped = 0 } = {}) => [
  `day ${day} (${(0.5 * day).toFixed(1)} min): Ts ${(14 + day / 100).toFixed(2)} °C, ASR ${(240 + day / 10).toFixed(1)} (atmosphere ${(75 + day / 10).toFixed(1)}) OLR ${(238.3 + day / 10).toFixed(1)} W/m², ps 505–1041 hPa, max wind ${(80 + day).toFixed(1)} m/s, precip ${(2.9 + day / 100).toFixed(2)} mm/d, `,
  fog ? 'fog 0.012 mm/d global (land 0.010, sea 0.013; 41 cells over 0.1 mm/d, mean there 0.180 mm/d), ' : '',
  `ice ${(3 + day / 10).toFixed(1)}% (N ${(15 + day / 10).toFixed(1)} S ${(5 + day / 10).toFixed(1)} Mkm²${label ? ', lakes left out' : ''}), albedo ${(0.3 + day / 1000).toFixed(3)}, `,
  cloud ? `SWCRE ${(-50 - day / 10).toFixed(1)} LWCRE ${(28 + day / 10).toFixed(1)}, clear-sky reflectance 0.1500, ` : '',
  `sea surface shortwave 170.2 net longwave -55.1 W/m² (iced cells poleward of 60°: N 20.1 S 80.2), ocean h1 ${30 + day} m, interior ${(3.1 + day / 100).toFixed(2)} °C, currents ≤ ${current.toFixed(2)} m/s, transport ${150 + day} Sv, clamped ${clamped}; balance day +${(1 + day / 10).toFixed(1)} (${day} of 365 d) W/m²`,
].join('');

const SEPTEMBER = 'day 1 (0.5 min): Ts 14.01 °C, ASR 240.1 OLR 238.4 W/m², ps 505–1041 hPa, max wind 81.0 m/s, precip 2.91 mm/d, ice 3.1% (N 15.1 S 5.1 Mkm²), albedo 0.301, ocean h1 31 m, interior 3.11 °C, currents ≤ 1.20 m/s, transport 151 Sv, clamped 0';
const LOGS = {
  fx64: [SEPTEMBER, ...[2, 3].map((day) => dayLine(day, { fog: false }))],
  fx128: [dayLine(1), dayLine(2), dayLine(3, { label: true, current: 5.2, clamped: 1 })],
  ns64: [1, 2, 3].map((day) => dayLine(day, { label: true, cloud: false })),
};
const QUARTER = 'Q1 days 1-3 (Mar 21–Mar 23): Ts 14.02 °C, TOA +1.7 W/m², albedo 0.302, rain 2.92 mm/d, ice N 15.1–15.3 S 5.1–5.3 Mkm² (mean N 15.2 S 5.2';

const withLogs = (body) => {
  const dir = mkdtempSync(join(tmpdir(), 'log-readers-'));
  try {
    for (const [name, lines] of Object.entries(LOGS)) writeFileSync(join(dir, `${name}.log`), lines.join('\n') + '\n');
    body(dir);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
};

test("rounds.py takes each field of spinup.mjs's daily lines from before the atmosphere's shortwave and the fog text, without the lakes' label and with it, and marks the quarters whose extents count frozen lakes", () => withLogs((dir) => {
  const rounds = (name) => execFileSync('python3', [join(root, 'scripts/rounds.py'), join(dir, `${name}.log`)], { encoding: 'utf8', env: { ...process.env, PER_YEAR: '4' } }).trim();
  assert.equal(rounds('fx64'), `fx64 ${QUARTER}; frozen lakes counted)`, 'a line from before the atmosphere\'s shortwave, and lines from before the fog text and the label');
  assert.equal(rounds('fx128'), `fx128 ${QUARTER}; frozen lakes counted on 2 of 3 days)`, 'two lines without the label and one with it');
  assert.equal(rounds('ns64'), `ns64 ${QUARTER})`, 'labeled lines without the cloud effects');
}));

test("elevenReport.mjs takes each field of spinup.mjs's newest daily line, from before the fog text, with the lakes' label and without the cloud effects, and flags clamped ocean edges and currents at the 5 m/s cap", () => withLogs((dir) => {
  const report = (prefix) => execFileSync(process.execPath, [join(root, 'scripts/elevenReport.mjs'), dir, prefix, `log-readers-test-${process.pid}`], { cwd: root, encoding: 'utf8', env: { ...process.env, PER_YEAR: '4' } }).trim().split('\n');

  const fx = report('fx');
  assert.match(fx[0], /^fx at \d\d:\d\d\s[AP]M: pulled N=128 day 3, N=64 day 3 \(Mar 23, year 1\); rounds done to day —; instance unknown, 0 recreations$/);
  assert.deepEqual(fx.slice(1), [
    'N=128 day 3: Ts 14.03 °C, ASR 240.3 OLR 238.6; balance day +1.3 (3 of 365 d); albedo 0.303, SWCRE -50.3 LWCRE 28.3; rain 2.93 mm/d; ice N 15.3 S 5.3 Mkm²; max wind 83 m/s; ocean h1 33 m, interior 3.13 °C, currents ≤ 5.20 m/s, transport 153 Sv, clamped 1',
    'N=64 day 3: Ts 14.03 °C, ASR 240.3 OLR 238.6; balance day +1.3 (3 of 365 d); albedo 0.303, SWCRE -50.3 LWCRE 28.3; rain 2.93 mm/d; ice N 15.3 S 5.3 Mkm² (frozen lakes counted); max wind 83 m/s; ocean h1 33 m, interior 3.13 °C, currents ≤ 1.20 m/s, transport 153 Sv, clamped 0',
    `fx128 ${QUARTER}; frozen lakes counted on 2 of 3 days)`,
    'FLAGS: N=128: 1 clamped ocean edges on day 3 | N=128: ocean currents at the 5 m/s cap on day 3',
  ], 'a labeled N=128 line with a clamped edge and a current over 5 m/s, and N=64 lines from before the fog text and the label');

  const ns = report('ns');
  assert.match(ns[0], /^ns at \d\d:\d\d\s[AP]M: pulled N=128 day —, N=64 day 3 \(Mar 23, year 1\); rounds done to day —; instance unknown, 0 recreations$/);
  assert.deepEqual(ns.slice(1), [
    'N=128: no day yet',
    'N=64 day 3: Ts 14.03 °C, ASR 240.3 OLR 238.6; balance day +1.3 (3 of 365 d); albedo 0.303; rain 2.93 mm/d; ice N 15.3 S 5.3 Mkm²; max wind 83 m/s; ocean h1 33 m, interior 3.13 °C, currents ≤ 1.20 m/s, transport 153 Sv, clamped 0',
    `ns64 ${QUARTER})`,
    'no flags: no NaN, no stop, no clamped edges, winds and Courant bounded, Ts drift within limits, both ice caps present',
  ], 'labeled N=64 lines without the cloud effects, and no N=128 log');
}));
