export const F32_RMS = 16 * 2 ** -24, F32_CELL = 64 * 2 ** -24;

export function distributionCheck(gaps, lat, { largest, counts, hemisphere }, { max, rms }) {
  const top = Math.max(largest, F32_CELL * max), floor = Math.max(hemisphere, F32_RMS * rms);
  let worst = 0, at = -1;
  gaps.forEach((gap, x) => { if (gap > worst) { worst = gap; at = x; } });
  const over = counts.map(([threshold]) => gaps.reduce((n, gap) => n + (gap >= threshold), 0));
  const means = [true, false].map((north) => { let sum = 0, n = 0; gaps.forEach((gap, x) => { if ((lat[x] >= 0) === north) { sum += gap; n++; } }); return sum / n; });
  const failed = [];
  if (!(worst <= top)) failed.push(`largest ${worst.toExponential(3)} at ${at}, past ${top}`);
  counts.forEach(([threshold, limit], j) => { if (over[j] > limit) failed.push(`${over[j]} at or over ${threshold}, past ${limit}`); });
  means.forEach((mean, j) => { if (!(mean <= floor)) failed.push(`${j ? 'southern' : 'northern'} mean ${mean.toExponential(3)}, past ${floor}`); });
  return { failed, worst, over, means, text: `largest ${worst.toExponential(2)} at ${at} of ${top.toExponential(2)}, ${counts.map(([threshold, limit], j) => `${over[j]} of ${limit} at or over ${threshold.toExponential(1)}`).join(', ')}, mean ${means.map((mean) => mean.toExponential(2)).join('/')} north/south of ${floor.toExponential(2)}` };
}
