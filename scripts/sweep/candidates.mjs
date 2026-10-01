// The sweep's candidates (runs/sweep/candidates.json, from fit.py) on the
// longer tests, scored with the screens' formula:
//   node scripts/sweep/candidates.mjs
// For each candidate: ten N=64 days from eight64_day0183 (TAG cx<name>t;
// the albedo's day means over days 186-193, ASR - OLR over 188-193, the
// audit of day 193 and its equatorial stress), five N=128 days from
// eight128_day0183 (cx<name>h; albedo over 186-188, ASR - OLR over
// 186-188, the equatorial stress) and thirty N=64 days from a fresh atlas
// start on bl34 (cx<name>f; day 30's albedo and ASR - OLR, NaN, clamps),
// with the three-day screens for a candidate the design did not run. The
// N=128 runs go beside the N=64 ones on the GPU. The score is the sum of
// the screens' formula over the day-193 numbers (with the screens' Arctic
// loss), its balance, albedo and stress terms over the day-188 numbers
// and its balance and albedo terms over day 30's; one row per candidate
// goes to runs/sweep/candidates.csv. ONLY (comma-separated names) limits
// the runs to those candidates.
import { existsSync, readFileSync, writeFileSync } from 'node:fs';
import { PARAMETERS, OUT, SWEEP, STATES, spinup, audit, screenValues } from './runs.mjs';
import { TERMS, score, readLog, readAudit } from './score.mjs';

const candidates = JSON.parse(readFileSync(`${SWEEP}/candidates.json`, 'utf8'));
const ONLY = process.env.ONLY ? process.env.ONLY.split(',') : null;
const chosen = candidates.filter((c) => !ONLY || ONLY.includes(c.name));
const mean = (xs) => xs.reduce((s, x) => s + x, 0) / xs.length;
const between = (days, a, b) => days.filter((d) => d.day >= a && d.day <= b);
const pick = (keys) => TERMS.filter((t) => keys.includes(t.key));

async function longN128(c) {
  const run = await spinup({ tag: `cx${c.name}h`, n: 128, from: `${STATES}/eight128_day0183.bin`, fromDay: 183, days: 5, options: c });
  const log = readLog(run.log), window = between(log.days, 186, 188);
  return {
    albedos: log.days.map((d) => d.meanAlbedo), balances: log.days.map((d) => d.meanAsr - d.meanOlr),
    albedo: mean(window.map((d) => d.meanAlbedo)), balance: mean(window.map((d) => d.meanAsr - d.meanOlr)), stress: log.stress, nan: log.nan, clamped: log.days.reduce((s, d) => s + d.clamped, 0),
  };
}

async function longN64(c) {
  const ten = await spinup({ tag: `cx${c.name}t`, n: 64, from: `${STATES}/eight64_day0183.bin`, fromDay: 183, days: 10, options: c });
  const log = readLog(ten.log);
  const values = existsSync(ten.state) ? readAudit(await audit(ten.state, c, `${OUT}/cx${c.name}t.audit`)) : {};
  Object.assign(values, {
    albedos: log.days.map((d) => d.meanAlbedo), balances: log.days.map((d) => d.meanAsr - d.meanOlr),
    albedo: mean(between(log.days, 186, 193).map((d) => d.meanAlbedo)), balance: mean(between(log.days, 188, 193).map((d) => d.meanAsr - d.meanOlr)), dayRain: mean(between(log.days, 188, 193).map((d) => d.meanPrecip)),
    stress: log.stress, nan: log.nan, clamped: log.days.reduce((s, d) => s + d.clamped, 0),
  });
  const fresh = await spinup({ tag: `cx${c.name}f`, n: 64, days: 30, options: c, env: { LEVELS: 'bl34' } });
  const flog = readLog(fresh.log), day30 = flog.days.find((d) => d.day === 30);
  const freshValues = { albedo: day30?.meanAlbedo ?? NaN, balance: day30 ? day30.meanAsr - day30.meanOlr : NaN, nan: flog.nan || !day30, clamped: flog.days.reduce((s, d) => s + d.clamped, 0), stress: flog.stress, albedos: flog.days.slice(-6).map((d) => d.meanAlbedo), balances: flog.days.slice(-6).map((d) => d.meanAsr - d.meanOlr) };
  let screen = null;
  if (c.point === null) {
    const m = await spinup({ tag: `cx${c.name}m`, n: 64, from: `${STATES}/eight64_day0183.bin`, fromDay: 183, days: 3, options: c });
    const a = await spinup({ tag: `cx${c.name}a`, n: 64, from: `${STATES}/nine64_day0091.bin`, fromDay: 91, days: 3, options: c });
    screen = await screenValues(`cx${c.name}m`, c, m.log, m.state, a.state);
  } else {
    const row = readFileSync(`${SWEEP}/results.csv`, 'utf8').trim().split('\n').map((l) => l.split(','));
    const header = row[0], mine = row.find((r) => Number(r[0]) === c.point);
    screen = Object.fromEntries(header.map((h, k) => [h, mine[k] === '' ? NaN : Number(mine[k])]));
  }
  return { ten: values, fresh: freshValues, screen };
}

const results = {};
const n128 = (async () => { for (const c of chosen) { const h = await longN128(c); results[c.name] = { ...results[c.name], h }; } })();
const n64 = (async () => { for (const c of chosen) { const long = await longN64(c); results[c.name] = { ...results[c.name], ...long }; } })();
await Promise.all([n128, n64]);

const lines = [];
for (const c of chosen) {
  const { ten, h, fresh, screen } = results[c.name];
  const s193 = score({ ...ten, arctic: screen.arctic });
  const s188 = score(h, pick(['balance', 'albedo', 'stress']));
  const s30 = score(fresh, pick(['balance', 'albedo']));
  const total = s193.total + s188.total + s30.total;
  results[c.name].scores = { s193, s188, s30, total, screen: score(screen).total };
  lines.push({ name: c.name, ...Object.fromEntries(PARAMETERS.map((p) => [p.key, c[p.key]])), screenScore: score(screen).total, s193: s193.total, s188: s188.total, s30: s30.total, total,
    ...Object.fromEntries(TERMS.map((t) => [t.key, t.key === 'arctic' ? screen.arctic : ten[t.key]])), evaporation: ten.evaporation, dayRain: ten.dayRain, namibiaLow: ten.namibiaLow, sepRuns: ten.sepRuns, sepInversion: ten.sepInversion,
    h128Albedo: h.albedo, h128Balance: h.balance, h128Stress: h.stress, h128Clamped: h.clamped, h128Nan: h.nan ? 1 : 0, freshAlbedo: fresh.albedo, freshBalance: fresh.balance, freshStress: fresh.stress, freshClamped: fresh.clamped, freshNan: fresh.nan ? 1 : 0,
    albedos193: ten.albedos.map((x) => x.toFixed(3)).join(' '), balances193: ten.balances.map((x) => x.toFixed(1)).join(' '), albedos188: h.albedos.map((x) => x.toFixed(3)).join(' '), balances188: h.balances.map((x) => x.toFixed(1)).join(' '), albedos30: fresh.albedos.map((x) => x.toFixed(3)).join(' '), balances30: fresh.balances.map((x) => x.toFixed(1)).join(' ') });
}
const file = `${SWEEP}/candidates.csv`, keys = Object.keys(lines[0]);
const old = existsSync(file) ? readFileSync(file, 'utf8').trim().split('\n').slice(1).filter((l) => !chosen.some((c) => l.startsWith(`${c.name},`))) : [];
writeFileSync(file, [keys.join(','), ...old, ...lines.map((l) => keys.map((k) => (typeof l[k] === 'number' ? Number(l[k].toPrecision(6)) : l[k])).join(','))].join('\n') + '\n');
const scoresFile = `${SWEEP}/candidates.scores.json`, kept = existsSync(scoresFile) ? JSON.parse(readFileSync(scoresFile, 'utf8')) : {};
writeFileSync(scoresFile, JSON.stringify({ ...kept, ...Object.fromEntries(Object.entries(results).map(([k, v]) => [k, v.scores])) }, null, 1));
for (const l of lines) console.log(`${l.name}: total ${l.total.toFixed(2)} (day 193 ${l.s193.toFixed(2)}, N=128 ${l.s188.toFixed(2)}, fresh ${l.s30.toFixed(2)}; screen ${l.screenScore.toFixed(2)})`);
