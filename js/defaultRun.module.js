/*
 * The saved runs the page starts from when its address names none, by
 * resolution. The device test's choice starts from the run at its
 * resolution when there is one, and otherwise from DEFAULT_RUN, regridded.
 */
export const DEFAULT_RUNS = { 128: 'runs/eleven128_day1825.parts.json', 64: 'runs/eleven64_day1825.parts.json' };
export const DEFAULT_RUN = DEFAULT_RUNS[128];
