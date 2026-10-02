// The second sweep's candidates (runs/sweep2/candidates.json, from fit.py)
// on the longer tests, each with two replicates whose sea drag is 1 + 1e-4
// (<name>p) and 1 - 1e-4 (<name>n) of the default:
//   node scripts/sweep/candidates2.mjs
// For each member: thirty N=64 days from a fresh atlas start on bl34 (TAG
// c2<member>f; the means of days 25-30 take the fresh terms of TERMS2),
// ten N=64 days from eight64_day0183 (c2<member>t; the balance's mean over
// days 188-193, the clear-sky albedo of day 186, the day-193 audit over
// three windows two hours apart and the equatorial stress), five N=128
// days from eight128_day0183 (c2<member>h; the balance's mean over days
// 186-188, the extra term hBalance) and three N=64 days from
// nine64_day0091 (c2<member>a; the Arctic loss). WORKERS (2)
// members run the N=64 tests at once, beside one running the N=128 test.
// One row per member goes to runs/sweep2/candidates.csv and the means and
// standard deviations over each candidate's three members, with the
// decision, to candidates.txt. ONLY (comma-separated candidate names)
// limits the runs to those candidates.
import { existsSync, readFileSync, writeFileSync } from 'node:fs';
import { PARAMETERS2, SWEEP2, STATES, spinup, auditWindows, arcticStart } from './runs.mjs';
import { TERMS2, score, readLog, dayMeans, arcticVolume } from './score.mjs';

export const TERMS2C = [...TERMS2, { key: 'hBalance', label: 'N=128 from eight128: ASR - OLR, days 186-188 (W/m2)', target: 0, tolerance: 3, weight: 2 }];
const WORKERS = Number(process.env.WORKERS ?? 2);
const ONLY = process.env.ONLY ? process.env.ONLY.split(',') : null;
const candidates = JSON.parse(readFileSync(`${SWEEP2}/candidates.json`, 'utf8')).filter((c) => !ONLY || ONLY.includes(c.name));
const members = candidates.flatMap((c) => [{ ...c, of: c.name }, { ...c, name: `${c.name}p`, of: c.name, dragScale: 1 + 1e-4 }, { ...c, name: `${c.name}n`, of: c.name, dragScale: 1 - 1e-4 }]);
const run = (member, kind, spec) => spinup({ tag: `c2${member.name}${kind}`, options: member, parameters: PARAMETERS2, ...spec });
const clampedOf = (log) => log.days.reduce((s, d) => s + d.clamped, 0);

async function lowRes(member) {
  const fresh = await run(member, 'f', { n: 64, days: 30, env: { LEVELS: 'bl34' } });
  const ten = await run(member, 't', { n: 64, from: `${STATES}/eight64_day0183.bin`, fromDay: 183, days: 10 });
  const arctic = await run(member, 'a', { n: 64, from: `${STATES}/nine64_day0091.bin`, fromDay: 91, days: 3 });
  const flog = readLog(fresh.log), f = dayMeans(flog, 25, 30), tlog = readLog(ten.log), t = dayMeans(tlog, 188, 193), alog = readLog(arctic.log), firstDay = tlog.days.find((d) => d.day === 186);
  const values = {
    fBalance: f.balance, fAlbedo: f.albedo, fOlr: f.olr, fSwcre: f.swcre, fLwcre: f.lwcre, fRain: f.rain, fAsr: f.asr, fNan: flog.nan ? 1 : 0, fClamped: clampedOf(flog),
    balance: t.balance, clearAlbedo: firstDay && !tlog.nan ? firstDay.clearAlbedo : NaN, tAlbedo: t.albedo, tOlr: t.olr, tSwcre: t.swcre, tLwcre: t.lwcre, tDayRain: t.rain, stress: tlog.stress ?? NaN, tNan: tlog.nan ? 1 : 0, tClamped: clampedOf(tlog),
  };
  const audit = existsSync(ten.state) && !tlog.nan ? auditWindows(ten.state, member, `${SWEEP2}/c2${member.name}t`, 64) : Promise.resolve({});
  const iceStart = await arcticStart();
  if (existsSync(arctic.state) && !alog.nan) { const iceEnd = await arcticVolume(arctic.state); Object.assign(values, { iceStart, iceEnd, arctic: (iceStart - iceEnd) / 3 }); }
  return { values, audit };
}

async function highRes(member) {
  const h = await run(member, 'h', { n: 128, from: `${STATES}/eight128_day0183.bin`, fromDay: 183, days: 5 });
  const log = readLog(h.log), m = dayMeans(log, 186, 188);
  return { hBalance: m.balance, hAlbedo: m.albedo, hOlr: m.olr, hSwcre: m.swcre, hLwcre: m.lwcre, hDayRain: m.rain, hStress: log.stress ?? NaN, hNan: log.nan ? 1 : 0, hClamped: clampedOf(log) };
}

const EXTRA = ['fAsr', 'fNan', 'fClamped', 'tAlbedo', 'tOlr', 'tSwcre', 'tLwcre', 'tDayRain', 'tNan', 'tClamped', 'hAlbedo', 'hOlr', 'hSwcre', 'hLwcre', 'hDayRain', 'hStress', 'hNan', 'hClamped',
  'rain', 'evaporation', 'peruRain', 'sepLwp', 'peruLwp', 'sepRuns', 'peruRuns', 'namibiaLow', 'sepInversion', 'zonalPeak', 'iceStart', 'iceEnd'];

if (import.meta.url === `file://${process.argv[1]}`) {
  const header = ['name', 'of', 'dragScale', ...PARAMETERS2.map((p) => p.key), ...TERMS2C.map((t) => t.key), ...TERMS2C.map((t) => `w_${t.key}`), 'score', 'scoreN64', ...EXTRA];
  const file = `${SWEEP2}/candidates.csv`;
  const written = existsSync(file) ? readFileSync(file, 'utf8').split('\n')[0] : header.join(',');
  if (written !== header.join(',')) throw new Error(`${file} has the columns ${written}, not ${header.join(',')}: move it aside before adding members`);
  const results = {}, audits = [];
  const queue64 = [...members], queue128 = [...members];
  const low = async () => {
    for (let m = queue64.shift(); m; m = queue64.shift()) {
      const { values, audit } = await lowRes(m);
      results[m.name] = { ...results[m.name], ...values };
      audits.push(audit.then((a) => { results[m.name] = { ...results[m.name], ...a }; }));
      console.log(`${m.name}: N=64 tests done`);
    }
  };
  const high = async () => { for (let m = queue128.shift(); m; m = queue128.shift()) { results[m.name] = { ...results[m.name], ...(await highRes(m)) }; console.log(`${m.name}: N=128 test done`); } };
  await Promise.all([high(), ...Array.from({ length: WORKERS }, low)]);
  await Promise.all(audits);

  const cell = (x) => (typeof x === 'number' ? (Number.isFinite(x) ? Number(x.toPrecision(6)) : '') : x ?? '');
  const old = existsSync(file) ? readFileSync(file, 'utf8').trim().split('\n').slice(1).filter((l) => !members.some((m) => l.startsWith(`${m.name},`))) : [];
  const lines = members.map((m) => {
    const v = results[m.name], s = score(v, TERMS2C), s64 = score(v, TERMS2);
    return [m.name, m.of, m.dragScale ?? 1, ...PARAMETERS2.map((p) => m[p.key]), ...TERMS2C.map((t) => v[t.key]), ...TERMS2C.map((t) => s.parts[t.key]), s.total, s64.total, ...EXTRA.map((k) => v[k])].map(cell).join(',');
  });
  writeFileSync(file, [header.join(','), ...old, ...lines].join('\n') + '\n');

  const rows = readFileSync(file, 'utf8').trim().split('\n').slice(1).map((l) => Object.fromEntries(l.split(',').map((x, k) => [header[k], x])));
  const stats = (xs) => { const m = xs.reduce((s, x) => s + x, 0) / xs.length; return { mean: m, sd: xs.length > 1 ? Math.sqrt(xs.reduce((s, x) => s + (x - m) ** 2, 0) / (xs.length - 1)) : NaN }; };
  const groups = [...new Set(rows.map((r) => r.of))].map((of) => {
    const g = rows.filter((r) => r.of === of);
    return { of, n: g.length, ...stats(g.map((r) => Number(r.score))), n64: stats(g.map((r) => Number(r.scoreN64))) };
  }).sort((a, b) => a.mean - b.mean);
  const out = groups.map((g) => `${g.of}: score ${g.mean.toFixed(1)} +- ${g.sd.toFixed(1)} over ${g.n} members (without the N=128 balance ${g.n64.mean.toFixed(1)} +- ${g.n64.sd.toFixed(1)})`);
  const best = groups[0], base = groups.find((g) => g.of === 'base');
  if (base) {
    const margin = Math.max(best.sd, base.sd);
    out.push(best.of === 'base' ? 'the defaults score lowest: they stay'
      : base.mean - best.mean > margin ? `winner ${best.of}: beats the defaults by ${(base.mean - best.mean).toFixed(1)}, more than the larger standard deviation ${margin.toFixed(1)}`
        : `no winner: ${best.of} beats the defaults by ${(base.mean - best.mean).toFixed(1)}, inside the larger standard deviation ${margin.toFixed(1)}; the defaults stay`);
  }
  writeFileSync(`${SWEEP2}/candidates.txt`, out.join('\n') + '\n');
  console.log(out.join('\n'));
}
