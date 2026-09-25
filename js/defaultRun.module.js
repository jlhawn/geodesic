/*
 * The saved runs the page starts from when its address names none, by
 * resolution. The device test's choice starts from the run at its
 * resolution when there is one, and otherwise from DEFAULT_RUN, regridded.
 */
export const DEFAULT_RUNS = { 128: 'runs/spin128c_day0678.parts.json' };
export const DEFAULT_RUN = DEFAULT_RUNS[128];
