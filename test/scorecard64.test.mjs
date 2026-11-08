import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { ROWS, BOXES, SIGMA_FILE, logDays, undercurrentCore, dayWeighted, monthAt, windowDays, windowMeans, pairSigma, compareScores } from '../scripts/scorecard64.mjs';

test('the log parser reads the day-mean TOA line and the daily sea-ice cut', () => {
  const text = [
    'day 410 (1.4 min): Ts 16.14 °C, ASR 231.5 (atmosphere 77.0) OLR 234.5 W/m², ps 539–1035 hPa, albedo 0.309, SWCRE -50.7 LWCRE 26.2, clear-sky reflectance 0.1573',
    'sea ice on day 410\'s mean concentration, sea cells, 15 %: N extent 12.14 area 10.73, S extent 1.59 area 0.89 Mkm²',
    'regions after 2 days',
  ].join('\n');
  const { days, ice } = logDays(text);
  assert.deepEqual(days.get(410), { asr: 231.5, olr: 234.5, swcre: -50.7, lwcre: 26.2 });
  assert.deepEqual(ice.get(410), { nExtent: 12.14, nArea: 10.73, sExtent: 1.59, sArea: 0.89 });
  assert.equal(days.size + ice.size, 2);
});

test('the undercurrent core is the profile\x27s eastward maximum refined by a parabola', () => {
  const profile = Array.from({ length: 60 }, (_, j) => { const z = 20 + 5 * j; return { z, u: 0.9 - ((z - 107.5) / 40) ** 2 }; });
  const { speed, depth } = undercurrentCore(profile);
  assert.ok(Math.abs(speed - 0.9) < 1e-12 && Math.abs(depth - 107.5) < 1e-9, `${speed} m/s at ${depth} m`);
  assert.deepEqual(undercurrentCore([{ z: 20, u: 0.5 }, { z: 25, u: 0.2 }, { z: 30, u: NaN }]), { speed: 0.5, depth: 20 });
});

test('window weights: day 0 is the March equinox and a window\x27s reference weighs each month by its hours', () => {
  assert.deepEqual([monthAt(0), monthAt(12), monthAt(365), monthAt(-1)], [3, 4, 3, 3]);
  assert.deepEqual(windowDays(425).slice(0, 2), [396, 397]);
  assert.equal(dayWeighted(395, () => 7), 7);
  const april = dayWeighted(395, (m) => (m === 4 ? 1 : 0));
  assert.ok(april > 0.6 && april < 0.75, `${april} of days 366–395 fall in April`);
});

test('a row moves only past twice its pair sigma; proxies stand in where a twin has no value', () => {
  const rows = (model, proxy = NaN) => Object.fromEntries(ROWS.map(([key]) => [key, { model, ref: 0, proxy }]));
  const a = [{ day: 395, rows: rows(1) }, { day: 425, rows: rows(3) }], b = [{ day: 395, rows: rows(NaN, 0) }, { day: 425, rows: rows(NaN, 2) }];
  const sigma = pairSigma(a, b);
  assert.ok(ROWS.every(([key]) => sigma[key].sigma === 1 && sigma[key].proxy === false));
  assert.ok(Number.isNaN(pairSigma(a, b, { proxy: false }).sst30.sigma));
  const base = windowMeans(a), moved = windowMeans([{ day: 395, rows: rows(4.1) }]), inside = windowMeans([{ day: 395, rows: rows(3.9) }]);
  assert.equal(base.sst30.model, 2);
  assert.ok(compareScores(base, moved, sigma).every((c) => c.verdict === 'moved'));
  assert.ok(compareScores(base, inside, sigma).every((c) => c.verdict === 'unresolved'));
  assert.equal(compareScores(base, windowMeans([{ day: 395, rows: rows(NaN) }]), sigma)[0].verdict, 'no value');
});

test('the committed sigma covers every row with its twins, windows and expiry, and the rows and boxes are named once', () => {
  const record = JSON.parse(readFileSync(SIGMA_FILE, 'utf8'));
  assert.ok(ROWS.every(([key]) => Number.isFinite(record.rows[key]?.sigma) && record.rows[key].sigma > 0), ROWS.filter(([key]) => !(record.rows[key]?.sigma > 0)).map(([key]) => key).join(', '));
  assert.ok(record.twins.length === 2 && record.n === record.days.length && record.expires && record.measured && record.definition);
  assert.equal(new Set(ROWS.map(([key]) => key)).size, ROWS.length);
  assert.ok(Object.values(BOXES).every(([kind, s, n, w, e]) => ['land', 'sea', 'all'].includes(kind) && s < n && w < e));
});
