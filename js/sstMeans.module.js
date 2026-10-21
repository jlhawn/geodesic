/*
 * Daily means of the sea surface temperature, and the ring of the last
 * SST_WINDOW of them that the spin-ups keep in their states, so that a
 * state's SST can be set against a climatology over a window of days and
 * not only at the one UTC hour it was saved, when each longitude stands
 * at its own local time and that hour's eddies and instability waves
 * count as bias.
 *
 * A recorder (createSstRecorder for a model of model.module.js,
 * createGpuSstRecorder in gpu/sstMeans.gpu.js for one of model.gpu.js)
 * adds up, after every ocean step, each ocean cell's top-layer
 * temperature less SST_REFERENCE; it reads the ocean and writes nothing
 * of the model's. It counts the model's steps from its own creation to
 * know which were ocean steps, so it must be created after the ocean is
 * loaded and before the model's first step, and `step` called after
 * every step. `day` resolves to the day's { mean, samples }: the mean in
 * K of the samples since the last `day` (0 off the ocean's cells), as
 * float32, and starts the next day; `checkpoint` resolves to the part of
 * a day summed so far, { sum, samples }, which a recorder created with it
 * as `from` carries on.
 *
 * The ring (createSstRing) holds whole consecutive days, oldest first; a
 * day pushed that does not follow the last starts it again. A state
 * carries it as `sstMeans`, the daily means one day after another, and
 * `sstMeanDays`, the model day each ends: the mean of day d covers model
 * time d − 1 to d days. A state saved inside a day adds its recorder's
 * checkpoint as `sstSum` and `sstSamples`.
 */
export const SST_WINDOW = 30, SST_REFERENCE = 273.15;

export const sstMeanOf = (sum, samples, ocean) => Float32Array.from(sum, (s, i) => (ocean[i] && samples > 0 ? SST_REFERENCE + s / samples : 0));

export function createSstRecorder(model, from = null) {
  const { ocean } = model, C = model.mesh.nCells, cells = ocean.cellOcean, sum = new Float64Array(C);
  if (from) sum.set(from.sum);
  let counted = 0, samples = from ? from.samples : 0;
  return {
    step() {
      if (++counted % ocean.everySteps !== 0) return;
      const T0 = ocean.T0;
      for (let i = 0; i < C; i++) if (cells[i]) sum[i] += T0[i] - SST_REFERENCE;
      samples++;
    },
    day() {
      const out = { mean: sstMeanOf(sum, samples, cells), samples };
      sum.fill(0);
      samples = 0;
      return Promise.resolve(out);
    },
    checkpoint() { return Promise.resolve({ sum: Float64Array.from(sum), samples }); },
  };
}

export function createSstRing(C, saved = null, length = SST_WINDOW) {
  const days = [], means = [];
  function push(day, mean) {
    if (days.length && day !== days[days.length - 1] + 1) { days.length = 0; means.length = 0; }
    days.push(day);
    means.push(Float32Array.from(mean));
    while (days.length > length) { days.shift(); means.shift(); }
  }
  const savedDays = saved && saved.sstMeanDays, savedMeans = saved && saved.sstMeans;
  if (savedDays && savedMeans && savedMeans.length === savedDays.length * C) for (let d = 0; d < savedDays.length; d++) push(savedDays[d], savedMeans.subarray(d * C, (d + 1) * C));
  return {
    length, push,
    get days() { return days.slice(); },
    serialize() {
      if (!days.length) return {};
      const sstMeans = new Float32Array(days.length * C);
      means.forEach((mean, d) => sstMeans.set(mean, d * C));
      return { sstMeans, sstMeanDays: Float64Array.from(days) };
    },
  };
}

/*
 * The mean over a saved state's ring when its last `length` days are
 * consecutive and end at the state's own day: { days, mean }, the days
 * oldest first and the mean in K, one value a cell; otherwise null.
 */
export function sstWindowOf(saved, length = SST_WINDOW) {
  const days = saved && saved.sstMeanDays ? Array.from(saved.sstMeanDays) : [], means = saved && saved.sstMeans, n = days.length;
  if (n < length || !means || means.length % n || days[n - 1] !== saved.day) return null;
  for (let d = n - length + 1; d < n; d++) if (days[d] !== days[d - 1] + 1) return null;
  const C = means.length / n, mean = new Float64Array(C);
  for (let d = n - length; d < n; d++) for (let i = 0; i < C; i++) mean[i] += means[d * C + i];
  for (let i = 0; i < C; i++) mean[i] /= length;
  return { days: days.slice(n - length), mean };
}
