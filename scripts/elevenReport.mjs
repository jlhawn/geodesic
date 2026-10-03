// A brief report on a paired spin-up from the files its puller has
// brought back (scripts/verdaPull.sh): the newest day of each resolution
// and its daily line, the energy balance windows, the change over the
// last 30 days, the stratosphere's top, the latest comparison round, the
// land jumps, the pace from the driver's log, the instance's status from
// the relauncher's log, and the signs of trouble it finds: NaN, a stop,
// clamped ocean edges, winds past WIND m/s or a Courant number past
// COURANT in the top layers, currents at the ocean's 5 m/s cap, a surface
// temperature moving more than DRIFT K in 30 days after day 120 (the
// fresh start warms through its first months), sea ice
// gone or past 25 Mkm² in a hemisphere, a driver silent for STALL minutes
// while the instance runs (the driver's log is stamped in UTC), and the quarter, year and rolling-year table of
// scripts/rounds.py for the leading resolution. Nothing here touches the run.
//   node scripts/elevenReport.mjs [runs/verda-eleven] [eleven] [gcm-eleven]
import { readFileSync, existsSync, statSync } from 'node:fs';
import { join } from 'node:path';
import { homedir } from 'node:os';
import { spawnSync } from 'node:child_process';

const DIR = process.argv[2] ?? 'runs/verda-eleven', PREFIX = process.argv[3] ?? 'eleven', NAME = process.argv[4] ?? 'gcm-eleven';
const PRICE = Number(process.env.PRICE ?? 1.8911), UNTIL = Number(process.env.UNTIL ?? 1095);
const WIND = 180, COURANT = 0.7, DRIFT = 1.5, STALL = 20;
const text = (file) => (existsSync(file) ? readFileSync(file, 'utf8') : '');
const DAY = /^day (\d+) \(([\d.]+) min\): Ts ([-\d.]+) °C, ASR ([-\d.]+) \(atmosphere ([-\d.]+)\) OLR ([-\d.]+) W\/m².*?max wind ([-\d.]+) m\/s, precip ([-\d.]+) mm\/d, ice ([-\d.]+)% \(N ([-\d.]+) S ([-\d.]+) Mkm²\), albedo ([-\d.]+)(?:, SWCRE ([-\d.]+) LWCRE ([-\d.]+))?.*?(?:ocean h1 (\d+) m, interior ([-\d.]+) °C, currents ≤ ([-\d.]+) m\/s, transport ([-\d.]+) Sv, clamped (\d+))?/;

function readLog(file) {
  const days = new Map(), lines = text(file).split('\n');
  let nan = null, jumps = [], upper = null, stratosphere = null;
  for (const line of lines) {
    const m = line.match(DAY);
    const b = m && line.match(/; balance (.*?) W\/m²/);
    if (m) days.set(+m[1], { day: +m[1], minutes: +m[2], ts: +m[3], asr: +m[4], olr: +m[6], wind: +m[7], precip: +m[8], ice: +m[9], iceN: +m[10], iceS: +m[11], albedo: +m[12], swcre: m[13] === undefined ? NaN : +m[13], lwcre: m[14] === undefined ? NaN : +m[14], h1: m[15] === undefined ? NaN : +m[15], interior: m[16] === undefined ? NaN : +m[16], current: m[17] === undefined ? NaN : +m[17], transport: m[18] === undefined ? NaN : +m[18], clamped: m[19] === undefined ? NaN : +m[19], balance: b ? b[1] : null });
    if (/NaN on day/.test(line) && !nan) nan = line;
    if (/^land jump at the end of day/.test(line)) jumps.push(line);
    if (/^upper winds day/.test(line)) upper = line;
    if (/^stratosphere day/.test(line)) stratosphere = line;
  }
  return { days: [...days.values()].sort((a, b) => a.day - b.day), nan, jumps, upper, stratosphere };
}

// The top layer of the upper-winds line: its max wind and Courant numbers.
function top(upper) {
  if (!upper) return null;
  const m = upper.match(/: ([\d.]+) hPa [^;]*? ([\d.]+) ([\d.]+) ([\d.]+) ([\d.]+) ([\d.]+)\/([\d.]+)/);
  if (!m) return null;
  let worst = 0, worstC = 0;
  for (const part of upper.split(';')) { const p = part.match(/([\d.]+) hPa .*? ([\d.]+) [\d.]+ [\d.]+ [\d.]+ ([\d.]+)\/([\d.]+)\s*$/); if (p) { worst = Math.max(worst, +p[2]); worstC = Math.max(worstC, +p[3], +p[4]); } }
  return { pressure: +m[1], wind: +m[2], courant: `${m[6]}/${m[7]}`, worstWind: worst, worstCourant: worstC, day: +upper.match(/day (\d+)/)[1] };
}

const calendar = (day) => new Date(Date.UTC(2001, 2, 20) + ((day % 365) * 86400e3)).toUTCString().slice(8, 11) + ' ' + new Date(Date.UTC(2001, 2, 20) + ((day % 365) * 86400e3)).getUTCDate();
const signed = (x, digits = 1) => `${x >= 0 ? '+' : ''}${x.toFixed(digits)}`;
const stamp = (s) => new Date(s.replace(' ', 'T') + (s.endsWith('Z') ? '' : '')).getTime();

export function report(now = new Date()) {
  const out = [], flags = [];
  const driver = text(join(DIR, `${PREFIX}.log`)).split('\n');
  const done = driver.map((l) => l.match(/^(\S+ \S+) day (\d+) done/)).filter(Boolean).map((m) => ({ at: new Date(m[1].replace(' ', 'T') + 'Z').getTime(), day: +m[2] }));
  const lastStop = driver.findLastIndex((l) => /stopped \(/.test(l)), lastStart = driver.findLastIndex((l) => /^\S+ \S+ paired spin-up of/.test(l));
  const stopped = lastStop > lastStart ? driver[lastStop] : null;
  const pulled = text(join(DIR, 'pull.log')).split('\n').filter((l) => /states (there|wanted)/.test(l)).pop();
  const relaunch = text(join(homedir(), `verda-relaunch-${NAME}.log`)).split('\n').filter(Boolean);
  const statusLines = relaunch.filter((l) => / is /.test(l)), evictions = relaunch.filter((l) => /creat/.test(l)).length;
  const started = text(join(DIR, `STARTED_${PREFIX}`)).trim();
  const logs = Object.fromEntries(['64', '128'].map((n) => [n, readLog(join(DIR, `${PREFIX}${n}.log`))]));
  const newest = Object.fromEntries(Object.entries(logs).map(([n, l]) => [n, l.days[l.days.length - 1] ?? null]));
  const main = newest['128'] ?? newest['64'], lead = newest['128'] ? '128' : '64';
  // pace from the rounds done in the last hour, else all of them
  let pace = null;
  if (done.length >= 2) {
    const recent = done.filter((d) => d.at >= done[done.length - 1].at - 3600e3), use = recent.length >= 2 ? recent : done;
    pace = (use[use.length - 1].at - use[0].at) / 1000 / (use[use.length - 1].day - use[0].day);
  }
  const startedAt = started ? new Date(started.split(' ')[0]).getTime() : null;
  const hours = startedAt ? (now.getTime() - startedAt) / 3600e3 : null;
  const lastDone = done.length ? done[done.length - 1] : null;
  const silent = lastDone ? (now.getTime() - lastDone.at) / 60e3 : null;
  const season = main ? calendar(main.day) : '';
  const eta = pace && lastDone ? new Date(lastDone.at + pace * (UNTIL - lastDone.day) * 1000) : null;
  out.push(`${PREFIX} at ${now.toLocaleTimeString('en-US', { hour: '2-digit', minute: '2-digit' })}: pulled N=128 day ${newest['128']?.day ?? '—'}, N=64 day ${newest['64']?.day ?? '—'}${main ? ` (${season}, year ${Math.floor((main.day - 1) / 365) + 1})` : ''}; rounds done to day ${lastDone?.day ?? '—'}${pace ? `, ${pace.toFixed(1)} s per paired day` : ''}${eta ? `, day ${UNTIL} at ${eta.toLocaleString('en-US', { weekday: 'short', hour: '2-digit', minute: '2-digit' })}` : ''}${hours !== null ? `; ${hours.toFixed(1)} h since the start, ~$${(hours * PRICE).toFixed(2)}` : ''}; instance ${statusLines.length ? statusLines[statusLines.length - 1].replace(/^\S+ \S+ /, '') : 'unknown'}, ${evictions} recreation${evictions === 1 ? '' : 's'}${stopped ? `; DRIVER ${stopped.replace(/^\S+ \S+ /, '')}` : ''}`);
  for (const n of ['128', '64']) {
    const l = logs[n], d = newest[n];
    if (!d) { out.push(`N=${n}: no day yet`); continue; }
    const back = l.days.filter((x) => x.day <= d.day - 30).pop();
    const drift = back ? d.ts - back.ts : null;
    const t = top(l.upper);
    out.push(`N=${n} day ${d.day}: Ts ${d.ts.toFixed(2)} °C${drift !== null ? ` (${signed(drift, 2)} K over ${d.day - back.day} d)` : ''}, ASR ${d.asr.toFixed(1)} OLR ${d.olr.toFixed(1)}${d.balance ? `; balance ${d.balance}` : ''}; albedo ${d.albedo.toFixed(3)}${Number.isFinite(d.swcre) ? `, SWCRE ${d.swcre.toFixed(1)} LWCRE ${d.lwcre.toFixed(1)}` : ''}; rain ${d.precip.toFixed(2)} mm/d; ice N ${d.iceN.toFixed(1)} S ${d.iceS.toFixed(1)} Mkm²; max wind ${d.wind.toFixed(0)} m/s${Number.isFinite(d.current) ? `; ocean h1 ${d.h1} m, interior ${d.interior.toFixed(2)} °C, currents ≤ ${d.current.toFixed(2)} m/s, transport ${d.transport.toFixed(0)} Sv, clamped ${d.clamped}` : ''}${t ? `; top ${t.pressure} hPa wind ${t.wind.toFixed(0)} m/s Courant ${t.courant} (worst in the top layers ${t.worstWind.toFixed(0)} m/s, ${t.worstCourant.toFixed(2)}, day ${t.day})` : ''}`);
    if (l.nan) flags.push(`N=${n}: ${l.nan.trim()}`);
    if (Number.isFinite(d.clamped) && d.clamped > 0) flags.push(`N=${n}: ${d.clamped} clamped ocean edges on day ${d.day}`);
    if (d.wind > WIND) flags.push(`N=${n}: max wind ${d.wind.toFixed(0)} m/s on day ${d.day}`);
    if (t && t.worstCourant > COURANT) flags.push(`N=${n}: Courant ${t.worstCourant.toFixed(2)} in the top layers on day ${t.day}`);
    if (Number.isFinite(d.current) && d.current >= 4.99) flags.push(`N=${n}: ocean currents at the 5 m/s cap on day ${d.day}`);
    if (drift !== null && d.day > 120 && Math.abs(drift) > DRIFT) flags.push(`N=${n}: Ts moved ${signed(drift, 2)} K over the last ${d.day - back.day} days`);
    if (d.day > 30 && (d.iceN < 0.5 || d.iceS < 0.5)) flags.push(`N=${n}: sea ice gone in a hemisphere (N ${d.iceN.toFixed(1)} S ${d.iceS.toFixed(1)} Mkm²)`);
    if (d.iceN > 25 || d.iceS > 25) flags.push(`N=${n}: sea ice past 25 Mkm² (N ${d.iceN.toFixed(1)} S ${d.iceS.toFixed(1)})`);
    for (const j of l.jumps.slice(-2)) out.push(`N=${n} ${j.split('; by band')[0]}`);
  }
  const compare = text(join(DIR, `${PREFIX}_compare.md`)).split(/^### /m).filter((c) => /Surface temperature \(°C\) \|[^|]+\|[^|]*\d[^|]*\|/.test(c)).pop();
  if (compare) {
    const want = ['Surface temperature (°C)', 'SST, west Pacific warm pool', 'SST, east Pacific cold tongue', 'Thermocline, equator', 'Surface wind, trades', 'Surface wind, strongest westerlies'];
    const rows = compare.split('\n').filter((r) => want.some((w) => r.includes(w))).map((r) => r.replace(/\*\*/g, '').split('|').map((c) => c.trim()).filter(Boolean)).map((c) => `${c[0].replace(/ \(.*?\)$/, '')}: ${c[1]} | ${c[2]}`);
    out.push(`compare ${compare.split('\n')[0].trim()} (N=64 | N=128): ${rows.join('; ')}`);
  }
  if (silent !== null && silent > STALL && !stopped && /running/.test(statusLines[statusLines.length - 1] ?? '')) flags.push(`the driver has not finished a round for ${silent.toFixed(0)} min while the instance runs`);
  if (stopped && !/day \d+ reached/.test(stopped)) flags.push(`the driver stopped: ${stopped}`);
  if (pulled) out.push(`puller: ${pulled.replace(/^\S+ \S+ /, '').replace(/; newest:.*/, '')}, last round ${pulled.slice(0, 19)}`);
  const rounds = spawnSync('python3', [join(import.meta.dirname, 'rounds.py'), join(DIR, `${PREFIX}${lead}.log`)], { encoding: 'utf8', env: { ...process.env, SEGMENTS: '0' } });
  if (rounds.status === 0 && rounds.stdout.trim()) out.push(rounds.stdout.trim());
  out.push(flags.length ? `FLAGS: ${flags.join(' | ')}` : 'no flags: no NaN, no stop, no clamped edges, winds and Courant bounded, Ts drift within limits, both ice caps present');
  return out.join('\n');
}

if (process.argv[1] && import.meta.url.endsWith(process.argv[1].split('/').pop())) console.log(report());
