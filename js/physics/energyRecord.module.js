/*
 * The planet's energy balance over the last model year: a ring of day
 * means of the absorbed sunlight, the outgoing longwave and the surface
 * temperature, one entry a model day, which the drivers add to at the end
 * of each whole day and which travels with the saved state. The balance
 * over the last day, week, month and year is read from it; only the year
 * cancels the seasonal cycle (Earth's global net swings about ±10 W/m²
 * through the year, absorbed sunlight peaking at perihelion, outgoing
 * longwave in the northern summer), the shorter windows show the trend.
 * The days held are consecutive and end at the newest; a day added out
 * of sequence (a run continued from an older snapshot, or days skipped)
 * restarts the sequence, and a day added again replaces its entry. A day
 * whose means do not cover it whole (a run continued from inside the
 * day) is added as a gap, NaN, which the windows leave out and count.
 *
 * Layout of `values` (float64): [count, newest day, then DAYS rows of
 * FIELDS], a day in the row day mod DAYS.
 */
export const RECORD_DAYS = 365;
export const FIELDS = ['asr', 'olr', 'ts'];
export const WINDOWS = { day: 1, week: 7, month: 30, year: RECORD_DAYS };

export function createEnergyRecord(days = RECORD_DAYS) {
  const F = FIELDS.length, values = new Float64Array(2 + days * F);
  const slot = (day) => 2 + (((day % days) + days) % days) * F;
  const record = {
    values,
    days,
    get count() { return values[0]; },
    get newest() { return values[1]; },
    add(day, entry) {
      if (!(Number.isInteger(day) && day > 0)) throw new Error(`a model day, not ${day}`);
      const { asr, olr, ts } = entry ?? { asr: NaN, olr: NaN, ts: NaN };
      const count = values[0], newest = values[1];
      values[0] = count === 0 || day === newest + 1 ? Math.min(count + 1, days) : day > newest || day <= newest - count ? 1 : count - (newest - day);
      values[1] = day;
      values.set([asr, olr, ts], slot(day));
    },
    // The entries of the last `width` days, oldest first, without the gaps.
    last(width) {
      const n = Math.min(width, values[0]), out = [];
      for (let d = values[1] - n + 1; d <= values[1]; d++) { const s = slot(d); if (Number.isFinite(values[s])) out.push({ day: d, asr: values[s], olr: values[s + 1], ts: values[s + 2] }); }
      return out;
    },
    // Means over each window of WINDOWS: { net, asr, olr, ts, n, width, from, to },
    // n the days the window holds, null for a window with none.
    windows() {
      const out = {};
      for (const [name, width] of Object.entries(WINDOWS)) {
        const rows = record.last(width);
        if (!rows.length) { out[name] = null; continue; }
        const mean = (f) => rows.reduce((s, r) => s + f(r), 0) / rows.length;
        out[name] = { net: mean((r) => r.asr - r.olr), asr: mean((r) => r.asr), olr: mean((r) => r.olr), ts: mean((r) => r.ts), n: rows.length, width, from: rows[0].day, to: rows[rows.length - 1].day };
      }
      return out;
    },
    load(saved) {
      if (saved.length !== values.length) throw new Error(`an energy record of ${values.length} values, not ${saved.length}`);
      values.set(saved);
    },
    clear() { values.fill(0); },
  };
  return record;
}

// One line for a log: the balance over each window, with the days a
// window short of its width holds.
export function balanceLine(record) {
  const w = record.windows();
  return Object.entries(WINDOWS).map(([name, width]) => {
    const x = w[name];
    if (!x) return `${name} —`;
    return `${name} ${x.net >= 0 ? '+' : ''}${x.net.toFixed(1)}${x.n < width ? ` (${x.n} of ${width} d)` : ''}`;
  }).join(', ');
}
