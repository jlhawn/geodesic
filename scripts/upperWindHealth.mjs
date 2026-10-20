// One line per simulated year of a fixedSst WINDS=1 log: the day reached, the pace at each save,
// the largest wind and Courant numbers of any upper layer, and the 5S-5N wind's range at 3/10/30 hPa.
//   node scripts/upperWindHealth.mjs RUN.log
import { readFileSync } from 'node:fs';
const text = readFileSync(process.argv[2], 'utf8').split('\n');
const years = new Map(), paces = [];
let first = null;
for (const line of text) {
  const save = line.match(/^day \+(\d+): .*?([\d.]+) s\/day/);
  if (save) paces.push(`${save[1]}:${save[2]}`);
  if (/NaN|non-finite|Error/.test(line)) console.log('PROBLEM', line.slice(0, 200));
  const m = line.match(/^upper winds day (\d+) \([^)]*\): (.*)$/);
  if (!m) continue;
  const day = Number(m[1]); first ??= day;
  const y = Math.floor((day - first) / 365), e = years.get(y) ?? { from: day, to: day, wind: 0, ch: 0, cv: 0, eq: { 3: [1e9, -1e9], 10: [1e9, -1e9], 30: [1e9, -1e9] } };
  e.to = day;
  const layers = m[2].split('; ').map((s) => s.split(' ')).map(([p, , w, max, , , , co]) => ({ p: Number(p), eq: Number(w.split('/')[2]), max: Number(max), ch: Number(co.split('/')[0]), cv: Number(co.split('/')[1]) }));
  for (const l of layers) { e.wind = Math.max(e.wind, l.max); e.ch = Math.max(e.ch, l.ch); e.cv = Math.max(e.cv, l.cv); }
  for (const p of [3, 10, 30]) {
    const l = layers.reduce((a, b) => (Math.abs(Math.log(b.p / p)) < Math.abs(Math.log(a.p / p)) ? b : a));
    e.eq[p] = [Math.min(e.eq[p][0], l.eq), Math.max(e.eq[p][1], l.eq)];
  }
  years.set(y, e);
}
for (const [y, e] of years) console.log(`year ${y + 1} (days ${e.from}-${e.to}): max wind ${e.wind} m/s, Courant ${e.ch}/${e.cv}; 5S-5N u range 3 hPa ${e.eq[3].join('..')}, 10 hPa ${e.eq[10].join('..')}, 30 hPa ${e.eq[30].join('..')}`);
console.log(`pace at saves (day:s/day) ${paces.join(' ')}`);
