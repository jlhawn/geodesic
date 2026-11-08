import { FIELDS } from './frames.module.js';
import { cpuFrameFields } from './cpuFrame.module.js';

/*
 * Means over the last MEANS_WINDOW days of per-cell fields of
 * frames.module.js, at the lowest layer and the ocean's top, so that a
 * saved state can be set against day-mean climatologies and not only at
 * the one UTC hour it was saved.
 *
 * A recorder (createIntervalMeans) samples its fields every `every` model
 * steps, counted from `offset`, through the model's sampleFrame on the GPU
 * engine (model.gpu.js; it folds and clears nothing a frame would) or
 * cpuFrameFields on the CPU engines, and adds each sample into a float64
 * sum. On the GPU a sample's read is queued behind the steps before it
 * and summed once it lands, so the steps after it are not held. `step`
 * goes after every model step, or after a batch of them with the batch's
 * count, which samples once at the batch's end if the batch passed a
 * multiple of `every`. `day` goes at each day's end: the day's mean joins
 * a ring of daily means, float32 and `depth` days deep, which the oldest
 * day leaves, unless the recorder began partway through the day without a
 * checkpoint. `checkpoint` resolves, once the queued reads have landed, to
 * { fields, every, days, daySamples, sums }: the ring's days oldest first,
 * each { samples, values } with the values by field, and the samples and
 * float64 sums by field of the day under way; a recorder created with one
 * as `from` carries it on. Its fields must be among those the model's
 * frames give, which availableMeans keeps of a list: the GPU engine's
 * frames give the surface stress only with the implicit drag.
 *
 * A state carries the means as the group `means`: the header's
 * { fields, every, days, samples }, the days the ring held and the
 * samples they took, and `means.<field>` the mean of the ring's daily
 * means as float32, NaN where any sample was (the ocean's fields over
 * land). A state saved for a run to carry the ring on adds the header's
 * daySamples, `means.ringSamples` (float32, each day's samples),
 * `means.<field>Ring` (float32, the daily means oldest first) and, when
 * the day under way has samples, `means.<field>Sum` (float64).
 */
export const MEANS = ['mslp', 'ps', 'speed', 'temperature', 'humidity', 'sst', 'surfaceT', 'precipitation', 'cloud', 'albedo', 'shortwave', 'longwave', 'ice', 'concentration', 'zonalStress', 'meridionalStress', 'latentHeat', 'sensibleHeat', 't2m'];
export const MEANS_PER_DAY = 64, MEANS_WINDOW = 30;
export const MEMORY_FIELDS = ['rain', 'vertical', 'upwelling'];

export const meansEvery = (dt) => Math.max(1, Math.round(86400 / dt / MEANS_PER_DAY));
export const availableMeans = (model, fields) => (model.frameFields ? fields.filter((name) => model.frameFields.has(name)) : fields.slice());

/*
 * The fields a MEANS setting names: the default list when unset, none for
 * 'none', else a JSON list of frame fields other than the memories.
 */
export function meansFields(setting) {
  if (setting === undefined || setting === '') return MEANS.slice();
  if (setting === 'none') return [];
  const fields = JSON.parse(setting);
  if (!Array.isArray(fields) || fields.some((name) => !(name in FIELDS) || MEMORY_FIELDS.includes(name))) throw new Error(`MEANS is 'none' or a JSON list of the fields of js/frames.module.js other than ${MEMORY_FIELDS.join(', ')}, not ${setting}`);
  return [...new Set(fields)];
}

export function createIntervalMeans(model, { fields = MEANS, every = 1, depth = MEANS_WINDOW, step = 0, offset = 0, from = null } = {}) {
  if (!(Number.isInteger(every) && every >= 1)) throw new Error(`every is a whole number of steps, not ${every}`);
  if (!(Number.isInteger(depth) && depth >= 1)) throw new Error(`depth is a whole number of days, not ${depth}`);
  const missing = fields.filter((name) => !availableMeans(model, [name]).length);
  if (missing.length) throw new Error(`the model's frames give no ${missing.join(', ')}`);
  const sampled = !!model.sampleFrame;
  let sums = {}, daySamples = 0, days = [], partial = offset > 0 && !from, counted = offset, landed = Promise.resolve();
  const add = (got) => {
    for (const name of fields) {
      const values = got[name];
      if (!values) throw new Error(`the frame has no ${name}`);
      const sum = sums[name] ??= new Float64Array(values.length);
      for (let x = 0; x < values.length; x++) sum[x] += values[x];
    }
    daySamples++;
  };
  const close = () => {
    if (daySamples && !partial) {
      days.push({ samples: daySamples, values: Object.fromEntries(fields.map((name) => [name, Float32Array.from(sums[name], (s) => s / daySamples)])) });
      if (days.length > depth) days.shift();
    }
    sums = {};
    daySamples = 0;
    partial = false;
  };
  if (from) {
    if (from.every !== every || from.fields.join() !== fields.join()) throw new Error(`the checkpoint samples ${from.fields.join(', ')} every ${from.every} steps, not ${fields.join(', ')} every ${every}`);
    days = from.days.slice(-depth).map(({ samples, values }) => ({ samples, values: Object.fromEntries(fields.map((name) => [name, Float32Array.from(values[name])])) }));
    daySamples = from.daySamples;
    if (daySamples) sums = Object.fromEntries(fields.map((name) => [name, Float64Array.from(from.sums[name])]));
  }
  const sample = () => {
    if (!sampled) { add(cpuFrameFields(model, { fields }, { step })); return; }
    const read = model.sampleFrame({ fields });
    landed = Promise.all([landed, read]).then(([, got]) => add(got.fields));
  };
  return {
    fields, every, depth,
    step(count = 1) {
      const before = counted;
      counted += count;
      if (Math.floor(counted / every) > Math.floor(before / every)) sample();
    },
    day() {
      if (sampled) landed = landed.then(close);
      else close();
    },
    async checkpoint() {
      await landed;
      return { fields: fields.slice(), every, days: days.slice(), daySamples, sums: Object.fromEntries(fields.map((name) => [name, Float64Array.from(sums[name] ?? [])])) };
    },
  };
}

/*
 * The `means` group of a state from a recorder's checkpoint, with the ring
 * and the day under way when a run is to carry them on (`inside`); {} for
 * a checkpoint with nothing to save.
 */
export function meansGroup(part, { inside = false } = {}) {
  const n = part ? part.days.length : 0;
  if (!n && !(inside && part && part.daySamples)) return {};
  const group = { fields: part.fields, every: part.every, days: n, samples: part.days.reduce((sum, day) => sum + day.samples, 0) };
  for (const name of n ? part.fields : []) {
    const mean = new Float64Array(part.days[0].values[name].length);
    for (const { values } of part.days) {
      const v = values[name];
      for (let x = 0; x < mean.length; x++) mean[x] += v[x];
    }
    group[name] = Float32Array.from(mean, (s) => s / n);
  }
  if (inside) {
    group.daySamples = part.daySamples;
    if (n) group.ringSamples = Float32Array.from(part.days, (day) => day.samples);
    for (const name of part.fields) {
      if (n) {
        const length = part.days[0].values[name].length, ring = new Float32Array(n * length);
        part.days.forEach(({ values }, d) => ring.set(values[name], d * length));
        group[`${name}Ring`] = ring;
      }
      if (part.daySamples) group[`${name}Sum`] = Float64Array.from(part.sums[name]);
    }
  }
  return { means: group };
}

/*
 * A recorder's checkpoint from a saved state's `means` group, when the
 * state was saved for a run to carry the ring on; otherwise null.
 */
export function savedMeansPart(saved) {
  const group = saved && saved.means;
  if (!group || !Array.isArray(group.fields) || !Number.isInteger(group.days) || !Number.isInteger(group.daySamples)) return null;
  const n = group.days, counts = group.ringSamples, rings = group.fields.map((name) => group[`${name}Ring`]), sums = group.fields.map((name) => group[`${name}Sum`]);
  if (n && !(counts && counts.length === n && rings.every((ring) => ring instanceof Float32Array && ring.length % n === 0))) return null;
  if (group.daySamples && !sums.every((sum) => sum instanceof Float64Array)) return null;
  const days = Array.from({ length: n }, (_, d) => ({ samples: counts[d], values: Object.fromEntries(group.fields.map((name, f) => [name, rings[f].subarray(d * rings[f].length / n, (d + 1) * rings[f].length / n)])) }));
  return { fields: group.fields.slice(), every: group.every, days, daySamples: group.daySamples, sums: Object.fromEntries(group.fields.map((name, f) => [name, group.daySamples ? sums[f] : new Float64Array(0)])) };
}

export function withoutMeans(saved) {
  const { means, ...kept } = saved;
  return kept;
}
