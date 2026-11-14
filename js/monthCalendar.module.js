/*
 * The twelve mid-month days on the model's axis. Day 0 is the March equinox,
 * day of year 79, and the year is 365 days. monthWeights(day) gives the two
 * months whose mid-month values bracket a model time and the weight of the
 * second, cyclic across December and January.
 */
export const YEAR = 365, DAY_OF_YEAR_AT_ZERO = 79;
export const MID_MONTH_DAY_OF_YEAR = [15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349];
export const MID_MONTH_DAYS = MID_MONTH_DAY_OF_YEAR.map((d) => d - DAY_OF_YEAR_AT_ZERO);

export function monthWeights(day) {
  const MID = MID_MONTH_DAY_OF_YEAR;
  let d = (((day + DAY_OF_YEAR_AT_ZERO) % YEAR) + YEAR) % YEAR;
  if (d < MID[0]) d += YEAR;
  for (let m = 0; m < 12; m++) {
    const a = MID[m], b = m === 11 ? MID[0] + YEAR : MID[m + 1];
    if (d >= a && d < b) return [m, (m + 1) % 12, (d - a) / (b - a)];
  }
  throw new Error(`no month for day ${day}`);
}
