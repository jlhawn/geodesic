/*
 * The distribution of one group's engine gaps against an ensemble of CPU
 * twins (S7, S9). A spread holds a group's bounds, each tier applying only
 * when its constants are given:
 * - largest: the largest |gap|;
 * - counts: [threshold, limit] rungs, the cells at or past each threshold;
 * - hemisphere: the mean |gap| per hemisphere, one bound or { north, south };
 * - bulk: { north, south } of { trimmed, signed }, the mean |gap| and the
 *   absolute mean signed gap over the cells left once a run drops its own
 *   largest TRIM share (rounded up) by |gap|.
 * buildSpread computes every one of them from the twins' gap arrays, at k
 * times the twins' worst (counts plus 2), rounded up to three figures and
 * floored at the f32 floors: F32_CELL of the reference's largest value for
 * a cell, F32_RMS of its rms for a mean. Its ladder halves the threshold
 * from the twins' largest gap M down to the cell floor, M/1024, or the last
 * rung whose limit stays under half the group's cells. A non-finite gap or
 * reference fails every check.
 *
 * Hemispheres: a cell whose center lies within EQUATOR of the equatorial
 * plane (|z| of the unit center) counts in both; a latitude array in place
 * of sides puts a cell north when its latitude is at least 0.
 */
export const F32_RMS = 16 * 2 ** -24, F32_CELL = 64 * 2 ** -24, TRIM = 0.05, LADDER_DEPTH = 10, EQUATOR = 1e-9;

export function magnitude(values) {
  let max = 0, sum = 0, n = 0;
  for (const v of values) { if (!Number.isFinite(v)) return { max: NaN, rms: NaN }; max = Math.max(max, Math.abs(v)); sum += v * v; n++; }
  return { max, rms: Math.sqrt(sum / Math.max(n, 1)) };
}

export function hemispheres(mesh, cells = null) {
  const { xCell } = mesh, list = cells ?? Array.from({ length: mesh.nCells }, (_, i) => i);
  return Int8Array.from(list, (i) => { const x = xCell[3 * i], y = xCell[3 * i + 1], z = xCell[3 * i + 2] / Math.hypot(x, y, xCell[3 * i + 2]); return z > EQUATOR ? 1 : z < -EQUATOR ? -1 : 0; });
}

export function columnLargest(gaps, C, layers) {
  return Float64Array.from({ length: C }, (_, i) => { let worst = 0; for (const k of layers) { const g = gaps[k * C + i]; if (!(Math.abs(g) <= Math.abs(worst))) worst = g; } return worst; });
}

export const roundUp = (x) => { if (!(x > 0) || !Number.isFinite(x)) return x; const unit = 10 ** (Math.floor(Math.log10(x)) - 2), up = +(Math.ceil(x / unit) * unit).toPrecision(3); return up >= x ? up : +(up + unit).toPrecision(3); };

const sidesOf = (where, n) => (where instanceof Int8Array ? where : Int8Array.from({ length: n }, (_, x) => (where[x] >= 0 ? 1 : -1)));
const inHemisphere = (side, north) => side === 0 || (side > 0) === north;

export function spreadStatistics(gaps, where, thresholds = [], trim = TRIM) {
  const n = gaps.length, sides = sidesOf(where, n);
  let finite = true, largest = 0, at = -1;
  for (let x = 0; x < n; x++) { const a = Math.abs(gaps[x]); if (!Number.isFinite(a)) { if (finite) { finite = false; largest = NaN; at = x; } } else if (finite && a > largest) { largest = a; at = x; } }
  const over = thresholds.map((threshold) => { let count = 0; for (let x = 0; x < n; x++) if (!(Math.abs(gaps[x]) < threshold)) count++; return count; });
  const hemisphere = {};
  for (const [name, north] of [['north', true], ['south', false]]) {
    const cells = [];
    for (let x = 0; x < n; x++) if (inHemisphere(sides[x], north)) cells.push(x);
    if (!cells.length) continue;
    const size = (x) => (Number.isFinite(gaps[x]) ? Math.abs(gaps[x]) : Infinity), sorted = cells.slice().sort((a, b) => size(b) - size(a)), kept = sorted.slice(Math.ceil(trim * cells.length));
    const mean = (list, f) => list.reduce((sum, x) => sum + f(gaps[x]), 0) / Math.max(list.length, 1);
    hemisphere[name] = { cells: cells.length, mean: mean(cells, Math.abs), signed: mean(cells, (g) => g), trimmed: mean(kept, Math.abs), trimmedSigned: Math.abs(mean(kept, (g) => g)), top: sorted.slice(0, 5) };
  }
  return { finite, largest, at, over, hemisphere };
}

export function buildSpread(runs, where, reference, { k = 2, trim = TRIM } = {}) {
  const n = runs[0].length, cellFloor = F32_CELL * reference.max, meanFloor = F32_RMS * reference.rms;
  if (!Number.isFinite(cellFloor) || !Number.isFinite(meanFloor)) throw new Error('the reference field is not finite');
  const first = runs.map((gaps) => spreadStatistics(gaps, where, [], trim));
  if (first.some((s) => !s.finite)) throw new Error('a twin gap is not finite');
  const M = Math.max(...first.map((s) => s.largest)), counts = [];
  for (let j = 0; j <= LADDER_DEPTH; j++) {
    const threshold = +(M * 2 ** -j).toPrecision(3);
    if (!(threshold >= cellFloor) || threshold <= 0) break;
    const worst = Math.max(...runs.map((gaps) => spreadStatistics(gaps, where, [threshold], trim).over[0])), limit = Math.ceil(k * worst) + 2;
    if (limit >= n / 2) break;
    counts.push([threshold, limit]);
  }
  const bound = (values, floor) => roundUp(Math.max(k * Math.max(...values), floor));
  const hemisphere = {}, bulk = {};
  for (const side of ['north', 'south']) {
    if (!first[0].hemisphere[side]) continue;
    hemisphere[side] = bound(first.map((s) => s.hemisphere[side].mean), meanFloor);
    bulk[side] = { trimmed: bound(first.map((s) => s.hemisphere[side].trimmed), meanFloor), signed: bound(first.map((s) => s.hemisphere[side].trimmedSigned), meanFloor) };
  }
  return { largest: bound([M], cellFloor), counts, hemisphere, bulk, twins: { largest: M, over: counts.map(([threshold]) => Math.max(...runs.map((gaps) => spreadStatistics(gaps, where, [threshold], trim).over[0]))) } };
}

export function distributionCheck(gaps, where, { largest, counts = [], hemisphere, bulk }, { max, rms }, { cells = null, trim = TRIM } = {}) {
  const s = spreadStatistics(gaps, where, counts.map(([threshold]) => threshold), trim), failed = [], id = (x) => (cells ? cells[x] : x);
  if (!s.finite) failed.push(`a non-finite gap at ${id(s.at)}`);
  if (!Number.isFinite(max) || !Number.isFinite(rms)) failed.push(`a non-finite reference: largest ${max}, rms ${rms}`);
  const cellFloor = F32_CELL * max, meanFloor = F32_RMS * rms;
  if (largest !== undefined) { const top = Math.max(largest, cellFloor); if (!(s.largest <= top)) failed.push(`largest ${s.largest} at ${id(s.at)}, past ${top}`); }
  counts.forEach(([threshold, limit], j) => {
    if (!(limit < gaps.length)) failed.push(`the count at or over ${threshold} allows ${limit} of the ${gaps.length} entries, so it bounds nothing`);
    else if (!(s.over[j] <= limit)) failed.push(`${s.over[j]} at or over ${threshold}, past ${limit}`);
  });
  for (const [name, h] of Object.entries(s.hemisphere)) {
    if (hemisphere !== undefined) { const b = Math.max(typeof hemisphere === 'number' ? hemisphere : hemisphere[name], meanFloor); if (!(h.mean <= b)) failed.push(`${name}ern mean ${h.mean}, past ${b}`); }
    if (bulk?.[name]) {
      const t = Math.max(bulk[name].trimmed, meanFloor), g = Math.max(bulk[name].signed, meanFloor);
      if (!(h.trimmed <= t)) failed.push(`${name}ern trimmed mean ${h.trimmed}, past ${t}`);
      if (!(h.trimmedSigned <= g)) failed.push(`${name}ern trimmed signed mean ${h.trimmedSigned}, past ${g}`);
    }
  }
  const e = (v) => (Number.isFinite(v) ? v.toExponential(2) : String(v));
  const text = [`largest ${e(s.largest)} at ${id(s.at)}${largest !== undefined ? ` of ${e(Math.max(largest, cellFloor))}` : ''}`, counts.length ? `ladder ${counts.map(([threshold, limit], j) => `${s.over[j]}/${limit} ≥ ${threshold.toExponential(2)}`).join(' ')}` : '',
    ...Object.entries(s.hemisphere).map(([name, h]) => `${name}: signed ${e(h.signed)}, mean ${e(h.mean)}, trimmed ${e(h.trimmed)}, trimmed signed ${e(h.trimmedSigned)}, top ${h.top.map((x) => `${id(x)} ${e(gaps[x])}`).join(' ')}`)].filter(Boolean).join('; ');
  return { failed, worst: s.largest, over: s.over, means: ['north', 'south'].map((name) => s.hemisphere[name]?.mean), statistics: s, text };
}
