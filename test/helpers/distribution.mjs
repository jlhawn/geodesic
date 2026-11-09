/*
 * The distribution of one group's engine gaps against an ensemble of CPU
 * twins (S7, S9). A spread holds a group's bounds, each tier applying when
 * its constants are given:
 * - largest: the largest |gap|;
 * - counts: [threshold, limit] rungs, the cells at or past each threshold;
 * - emptyLadder: why counts is empty, { why: 'floor', threshold } when the
 *   top rung would sit under the cell floor or at 0, or { why: 'half',
 *   threshold, limit } when its limit would reach half the group;
 * - hemisphere: the mean |gap| per hemisphere, one bound or { north, south };
 * - hemisphereDropped: { north?, south? } of { largest, median, loo }, a side
 *   whose mean bound buildSpread leaves out under dropEventHemispheres,
 *   because single columns' events set the twins' means there: their largest
 *   is over EVENT_SPREAD times their median, or a twin's mean passes the
 *   bound the other twins set (loo, the largest share of that bound, over 1);
 * - bulk: { north, south } of { trimmed, signed }, the mean |gap| and the
 *   absolute mean signed gap over the cells left once a run drops its own
 *   largest TRIM share (rounded up) by |gap|.
 * distributionCheck fails a tier its require names that the spread lacks;
 * an S9 group requires S9_TIERS. buildSpread computes every tier from the
 * twins' gap arrays, at k times the twins' worst (counts plus 2), rounded up
 * to three figures and floored at the f32 floors: F32_CELL of the
 * reference's largest value for a cell, F32_RMS of its rms for a mean. Its
 * ladder halves the threshold from the twins' largest gap M down to the cell
 * floor, M/1024, or the last rung whose limit stays under half the group's
 * cells. Thresholds round down to three figures, so the top rung counts the
 * twin that set M. A check fails an empty group, a non-finite gap, reference
 * or bound, a hemisphere's bound that is missing for cells the group has
 * there or given where it has none, an empty ladder that a required counts
 * tier gives no reason for, a reason the group's floor or size contradicts
 * or that the largest bound does not tie to (the bound at least the
 * reason's threshold and at most K_MOST, S9(f)'s ceiling on k, times the
 * larger of that threshold and the cell floor, past the two roundings; a
 * half reason's limit at most K_MOST times the group plus 2),
 * and a dropped hemisphere that has a bound too or whose reason falls short;
 * buildSpread throws on an empty group or a non-finite gap.
 *
 * Hemispheres: `where` holds one value per gap, a side as hemispheres()
 * gives it or a latitude: positive is north, negative south, and 0 counts in
 * both. hemispheres() gives 0 to a cell whose center lies within EQUATOR of
 * the equatorial plane (|z| of the unit center).
 */
export const F32_RMS = 16 * 2 ** -24, F32_CELL = 64 * 2 ** -24, TRIM = 0.05, LADDER_DEPTH = 10, EQUATOR = 1e-9, EVENT_SPREAD = 3, K_MOST = 3;
export const S9_TIERS = ['largest', 'counts', 'hemisphere', 'bulk'];

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
  return Float64Array.from({ length: C }, (_, i) => { let worst = 0; for (const k of layers) { const g = gaps[k * C + i]; if (Number.isNaN(worst)) break; if (!(Math.abs(g) <= Math.abs(worst))) worst = g; } return worst; });
}

const nextFigure = (x, by) => {
  const [mantissa, exponent] = x.toExponential(2).split('e'), digits = Math.round(mantissa * 100) + by;
  return digits < 100 ? +`999e${exponent - 3}` : +`${digits}e${exponent - 2}`;
};
export const roundUp = (x) => { if (!(x > 0) || !Number.isFinite(x)) return x; const near = +x.toPrecision(3); return near >= x ? near : nextFigure(near, 1); };
export const roundDown = (x) => { if (!(x > 0) || !Number.isFinite(x)) return x; const near = +x.toPrecision(3); return near <= x ? near : nextFigure(near, -1); };

const sidesOf = (where, n) => {
  if (where.length !== n) throw new Error(`where holds ${where.length} entries for ${n} gaps`);
  if (where instanceof Int8Array) return where;
  return Int8Array.from({ length: n }, (_, x) => { const v = where[x]; if (!Number.isFinite(v)) throw new Error(`where holds ${v} at ${x}`); return v > 0 ? 1 : v < 0 ? -1 : 0; });
};
const inHemisphere = (side, north) => side === 0 || (side > 0) === north;
const middle = (values) => { const s = [...values].sort((a, b) => a - b), m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };
const eventSet = ({ largest, median, loo } = {}) => largest > EVENT_SPREAD * median || loo > 1;
const ROUNDING = 1.01 / 0.99;
const reasonHolds = (reason, largest, cellFloor, n) => {
  const { why, threshold, limit } = reason ?? {}, tied = threshold <= largest && largest <= K_MOST * ROUNDING * Math.max(threshold, cellFloor);
  return why === 'floor' ? tied && (threshold === 0 || (threshold > 0 && threshold < cellFloor)) : why === 'half' && tied && threshold > 0 && limit >= n / 2 && limit <= K_MOST * n + 2;
};

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

export function buildSpread(runs, where, reference, { k = 2, trim = TRIM, dropEventHemispheres = false } = {}) {
  if (!runs.length) throw new Error('no twins to build a spread from');
  const n = runs[0].length, cellFloor = F32_CELL * reference.max, meanFloor = F32_RMS * reference.rms;
  if (!n) throw new Error('an empty group');
  if (runs.some((gaps) => gaps.length !== n)) throw new Error('the twins hold different numbers of gaps');
  if (!(k > 0) || !Number.isFinite(k)) throw new Error(`k is ${k}`);
  if (!Number.isFinite(cellFloor) || !Number.isFinite(meanFloor)) throw new Error('the reference field is not finite');
  const first = runs.map((gaps) => spreadStatistics(gaps, where, [], trim));
  if (first.some((s) => !s.finite)) throw new Error('a twin gap is not finite');
  const M = Math.max(...first.map((s) => s.largest)), counts = [];
  let stop = null;
  for (let j = 0; j <= LADDER_DEPTH; j++) {
    const threshold = roundDown(M * 2 ** -j);
    if (!(threshold >= cellFloor) || threshold <= 0) { stop = { why: 'floor', threshold }; break; }
    const worst = Math.max(...runs.map((gaps) => spreadStatistics(gaps, where, [threshold], trim).over[0])), limit = Math.ceil(k * worst) + 2;
    if (limit >= n / 2) { stop = { why: 'half', threshold, limit }; break; }
    counts.push([threshold, limit]);
  }
  const bound = (values, floor) => roundUp(Math.max(k * Math.max(...values), floor));
  const hemisphere = {}, hemisphereDropped = {}, bulk = {};
  for (const side of ['north', 'south']) {
    if (!first[0].hemisphere[side]) continue;
    const means = first.map((s) => s.hemisphere[side].mean), loo = means.length > 1 ? Math.max(...means.map((v, t) => v / bound(means.filter((_, s) => s !== t), meanFloor))) : 0;
    const reason = { largest: roundUp(Math.max(...means)), median: roundDown(middle(means)), loo: roundUp(loo) };
    if (dropEventHemispheres && eventSet(reason)) hemisphereDropped[side] = reason;
    else hemisphere[side] = bound(means, meanFloor);
    bulk[side] = { trimmed: bound(first.map((s) => s.hemisphere[side].trimmed), meanFloor), signed: bound(first.map((s) => s.hemisphere[side].trimmedSigned), meanFloor) };
  }
  return { largest: bound([M], cellFloor), counts, ...(counts.length ? {} : { emptyLadder: stop }), hemisphere, ...(Object.keys(hemisphereDropped).length ? { hemisphereDropped } : {}), bulk, twins: { largest: M, over: counts.map(([threshold]) => Math.max(...runs.map((gaps) => spreadStatistics(gaps, where, [threshold], trim).over[0]))) } };
}

export function distributionCheck(gaps, where, spread, { max, rms }, { cells = null, trim = TRIM, require = [] } = {}) {
  const { largest, counts = [], emptyLadder, hemisphere, hemisphereDropped = {}, bulk } = spread, failed = [], id = (x) => (cells ? cells[x] : x);
  if (cells && cells.length !== gaps.length) throw new Error(`${cells.length} cell ids for ${gaps.length} gaps`);
  const s = spreadStatistics(gaps, where, counts.map(([threshold]) => threshold), trim);
  const finite = (name, bound) => Number.isFinite(bound) || (failed.push(`the ${name} bound is ${bound}`), false);
  for (const tier of require) if (spread[tier] === undefined) failed.push(`the spread has no ${tier} tier`);
  if (!gaps.length) failed.push('an empty group');
  if (!s.finite) failed.push(`a non-finite gap at ${id(s.at)}`);
  if (!Number.isFinite(max) || !Number.isFinite(rms)) failed.push(`a non-finite reference: largest ${max}, rms ${rms}`);
  const cellFloor = F32_CELL * max, meanFloor = F32_RMS * rms;
  if (emptyLadder !== undefined) {
    if (counts.length) failed.push('a ladder with rungs and a reason for having none');
    else if (!reasonHolds(emptyLadder, largest, cellFloor, gaps.length)) failed.push(`no ladder, for a reason the group contradicts: ${JSON.stringify(emptyLadder)} against a cell floor of ${cellFloor}, a largest bound of ${largest} and ${gaps.length} entries`);
  } else if (require.includes('counts') && spread.counts?.length === 0) failed.push('an empty ladder with no reason given');
  if (largest !== undefined && finite('largest', largest)) { const top = Math.max(largest, cellFloor); if (!(s.largest <= top)) failed.push(`largest ${s.largest} at ${id(s.at)}, past ${top}`); }
  counts.forEach(([threshold, limit], j) => {
    if (!finite(`rung ${j} threshold`, threshold)) return;
    if (!(limit < gaps.length)) failed.push(`the count at or over ${threshold} allows ${limit} of the ${gaps.length} entries, so it bounds nothing`);
    else if (!(s.over[j] <= limit)) failed.push(`${s.over[j]} at or over ${threshold}, past ${limit}`);
  });
  const sidesBounded = (tier) => (tier === undefined || typeof tier === 'number' ? [] : Object.keys(tier));
  for (const name of new Set([...sidesBounded(hemisphere), ...sidesBounded(hemisphereDropped ?? {}), ...sidesBounded(bulk)])) if (!s.hemisphere[name]) failed.push(`${name}ern bounds for a group with no ${name}ern cells`);
  for (const [name, h] of Object.entries(s.hemisphere)) {
    if (hemisphere !== undefined) {
      const b = typeof hemisphere === 'number' ? hemisphere : hemisphere[name], dropped = hemisphereDropped?.[name];
      if (dropped !== undefined) {
        if (b !== undefined) failed.push(`a ${name}ern mean bound and a reason for having none`);
        else if (!eventSet(dropped ?? {})) failed.push(`no ${name}ern mean bound, for a reason short of single columns' events: ${JSON.stringify(dropped)}`);
      } else if (b === undefined) failed.push(`no ${name}ern mean bound`);
      else if (finite(`${name}ern mean`, b)) { const top = Math.max(b, meanFloor); if (!(h.mean <= top)) failed.push(`${name}ern mean ${h.mean}, past ${top}`); }
    }
    if (bulk === undefined) continue;
    if (!bulk[name]) { failed.push(`no ${name}ern bulk bounds`); continue; }
    if (finite(`${name}ern trimmed mean`, bulk[name].trimmed)) { const t = Math.max(bulk[name].trimmed, meanFloor); if (!(h.trimmed <= t)) failed.push(`${name}ern trimmed mean ${h.trimmed}, past ${t}`); }
    if (finite(`${name}ern trimmed signed mean`, bulk[name].signed)) { const g = Math.max(bulk[name].signed, meanFloor); if (!(h.trimmedSigned <= g)) failed.push(`${name}ern trimmed signed mean ${h.trimmedSigned}, past ${g}`); }
  }
  const e = (v) => (Number.isFinite(v) ? v.toExponential(2) : String(v));
  const text = [`largest ${e(s.largest)} at ${id(s.at)}${largest !== undefined ? ` of ${e(Math.max(largest, cellFloor))}` : ''}`, counts.length ? `ladder ${counts.map(([threshold, limit], j) => `${s.over[j]}/${limit} ≥ ${threshold.toExponential(2)}`).join(' ')}` : emptyLadder ? `no ladder: ${emptyLadder.why === 'floor' ? `the twins' largest gap rounds to ${e(emptyLadder.threshold)}, under the cell floor ${e(cellFloor)}` : `its top rung, ${e(emptyLadder.threshold)}, would allow ${emptyLadder.limit} of ${gaps.length}`}` : '',
    ...Object.entries(s.hemisphere).map(([name, h]) => `${name}: signed ${e(h.signed)}, mean ${e(h.mean)}${hemisphereDropped?.[name] ? ` unbounded (the twins' means: largest ${e(hemisphereDropped[name].largest)}, median ${e(hemisphereDropped[name].median)}, a twin at ${hemisphereDropped[name].loo} of the others' bound)` : ''}, trimmed ${e(h.trimmed)}, trimmed signed ${e(h.trimmedSigned)}, top ${h.top.map((x) => `${id(x)} ${e(gaps[x])}`).join(' ')}`)].filter(Boolean).join('; ');
  return { failed, worst: s.largest, over: s.over, means: ['north', 'south'].map((name) => s.hemisphere[name]?.mean), statistics: s, text };
}
