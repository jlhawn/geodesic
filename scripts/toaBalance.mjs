// The top-of-atmosphere energy balance of a spin-up log over its last day,
// week, month and year, and over every whole model year, from the daily
// lines (each a mean over the day's steps). The windows are the last 1, 7,
// 30 and 365 logged days; only the 365-day window is a balance, the
// shorter ones carry the seasonal cycle (Earth's global net swings about
// ±10 W/m² through the year). A window with days missing from the log is
// reported with the days it has. Several logs are read in order, a day
// logged again (a run continued from an earlier snapshot) taking its last
// line. Logs written before the daily line became a day mean (its header
// line says so) hold one sample a day at a fixed time, and the report says
// that their means are aliased.
//   node scripts/toaBalance.mjs runs/eleven128.log
//   node scripts/toaBalance.mjs runs/verda-eleven/eleven64.log runs/verda-eleven/eleven64_b.log
import { readFileSync } from 'node:fs';

const DAY = /^day (\d+) .*Ts ([-\d.]+) °C, ASR ([-\d.]+) (?:\(atmosphere ([-\d.]+)\) )?OLR ([-\d.]+) W\/m².*precip ([-\d.]+) mm\/d.*albedo ([-\d.]+)/;

export function readDays(files) {
  const byDay = new Map();
  let dayMeans = false;
  for (const file of files) for (const line of readFileSync(file, 'utf8').split('\n')) {
    if (/ are day means over the day's /.test(line)) dayMeans = true;
    const m = line.match(DAY);
    if (!m) continue;
    const row = { day: +m[1], ts: +m[2], asr: +m[3], atmosphere: m[4] === undefined ? NaN : +m[4], olr: +m[5], precip: +m[6], albedo: +m[7] };
    const c = line.match(/SWCRE ([-\d.]+) LWCRE ([-\d.]+)/);
    if (c) Object.assign(row, { swcre: +c[1], lwcre: +c[2] });
    byDay.set(row.day, row);
  }
  const days = [...byDay.values()].sort((a, b) => a.day - b.day);
  days.dayMeans = dayMeans;
  return days;
}

const mean = (rows, f) => rows.reduce((s, r) => s + f(r), 0) / rows.length;

export function summary(rows) {
  const r = { last: rows[0], first: rows[0], n: rows.length, net: mean(rows, (d) => d.asr - d.olr), asr: mean(rows, (d) => d.asr), olr: mean(rows, (d) => d.olr), albedo: mean(rows, (d) => d.albedo), precip: mean(rows, (d) => d.precip), ts: mean(rows, (d) => d.ts) };
  r.last = rows[rows.length - 1];
  if (rows.every((d) => Number.isFinite(d.swcre))) Object.assign(r, { swcre: mean(rows, (d) => d.swcre), lwcre: mean(rows, (d) => d.lwcre) });
  return r;
}

// The last `width` days before and including `end`, as logged.
export function window(days, end, width) {
  return days.filter((d) => d.day > end - width && d.day <= end);
}

export function report(days) {
  if (!days.length) return 'no daily lines';
  const end = days[days.length - 1].day, lines = days.dayMeans ? [] : ['the daily lines are samples at one time of day, not day means: the means below are aliased by the sampling time'];
  const line = (label, rows, width) => {
    const s = summary(rows), cre = Number.isFinite(s.swcre) ? `, SWCRE ${s.swcre.toFixed(1)} LWCRE ${s.lwcre.toFixed(1)}` : '';
    lines.push(`${label.padEnd(26)} net ${s.net >= 0 ? '+' : ''}${s.net.toFixed(2)} W/m² (ASR ${s.asr.toFixed(1)} OLR ${s.olr.toFixed(1)}), albedo ${s.albedo.toFixed(3)}${cre}, rain ${s.precip.toFixed(2)} mm/d, Ts ${s.ts.toFixed(2)} °C${rows.length < width ? ` [${rows.length} of ${width} days]` : ''}`);
  };
  for (const [label, width] of [['last day', 1], ['last 7 days', 7], ['last 30 days', 30], ['last 365 days', 365]]) {
    const rows = window(days, end, width);
    if (rows.length) line(`${label} (to day ${end})`, rows, width);
  }
  for (let y = 1; y * 365 <= end; y++) {
    const rows = window(days, y * 365, 365);
    if (rows.length) {
      line(`year ${y} (days ${(y - 1) * 365 + 1}–${y * 365})`, rows, 365);
      const first = rows[0], last = rows[rows.length - 1];
      if (rows.length === 365) lines[lines.length - 1] += `, Ts day ${first.day} → ${last.day} ${first.ts.toFixed(2)} → ${last.ts.toFixed(2)} °C`;
    }
  }
  return lines.join('\n');
}

if (process.argv[1] && import.meta.url.endsWith(process.argv[1].split('/').pop())) {
  const files = process.argv.slice(2);
  if (!files.length) { console.error('usage: node scripts/toaBalance.mjs <spin-up log> [more logs]'); process.exit(2); }
  console.log(report(readDays(files)));
}
