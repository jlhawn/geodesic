/*
 * The saved runs the page starts from when its address names none, by
 * resolution. A chosen or given resolution starts from defaultRunFor(N),
 * regridded when no run is saved at that resolution.
 */
export const DEFAULT_RUNS = { 128: 'runs/eleven128_day1825.parts.json', 64: 'runs/eleven64_day1825.parts.json' };
export const DEFAULT_RUN = DEFAULT_RUNS[128];

/*
 * The saved run for a resolution: the one at the lowest resolution that is
 * at least N, so a coarser model regrids a finer state down without
 * fetching more than it needs, and above the highest the highest.
 */
export function defaultRunFor(N, runs = DEFAULT_RUNS) {
  const resolutions = Object.keys(runs).map(Number).sort((a, b) => a - b);
  return runs[resolutions.find((n) => n >= N) ?? resolutions[resolutions.length - 1]];
}
