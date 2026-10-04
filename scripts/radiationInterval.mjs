// What holding the radiation between full calls (the radiation's
// radiationEvery) does to the climate of a saved state, on the GPU.
//   node scripts/radiationInterval.mjs run <state.bin> <out.json>
// steps DAYS (3) model days from the state with the options of
// scripts/figures/figureState.mjs (RADIATION carries radiationEvery),
// with every θ scaled by 1 + PERTURB·u (u uniform in ±1/2, from a fixed
// seed) when PERTURB is set, for a replicate, and writes out.json: each
// day's global means (ASR, OLR, albedo, SWCRE, LWCRE, the surface's
// absorbed shortwave and net longwave, rain, Ts), the diurnal cycle of
// the surface temperature over land and over sea (each cell's departure
// from its own mean, binned by local solar hour and averaged by area), and
// beside it out.json.bin, the per-cell day-mean ASR and OLR of each day and
// the surface net radiation (absorbed shortwave plus net longwave) of
// every step of the first day, float32.
//   node scripts/radiationInterval.mjs compare <reference.json> <run.json>...
// prints each run's day means less the reference's, the amplitude and the
// local hour of the peak of its diurnal cycles, the area-weighted rms over
// cells of its day-mean ASR and OLR less the reference's, and the rms over
// cells of its surface net radiation less the reference's, averaged over
// the first day's steps.
import { readFileSync, writeFileSync } from 'node:fs';

const [mode, ...files] = process.argv.slice(2);
const HOURS = 24;

if (mode === 'run') {
  const { gpuModelFrom, readState } = await import('./figures/figureState.mjs');
  const { readRanges } = await import('../js/gpu/device.module.js');
  const { sunDirection } = await import('../js/physics/radiation.module.js');
  const [stateFile, out] = files;
  const DAYS = Number(process.env.DAYS ?? 3), PERTURB = Number(process.env.PERTURB ?? 0);
  const saved = await readState(stateFile);
  if (PERTURB) {
    let seed = 12345;
    const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
    saved.theta = Float64Array.from(saved.theta, (x) => x * (1 + PERTURB * (random() - 0.5)));
  }
  const model = await gpuModelFrom(saved);
  const { mesh } = model, C = mesh.nCells, dt = 1350 * 16 / saved.N, perDay = Math.round(86400 / dt);
  const { PH, S } = model.gpu.layout, land = model.geography.land, lheat = model.moist.latentHeat;
  const sums = new Float64Array(HOURS * C), counts = new Float64Array(HOURS * C), sun = new Float64Array(3);
  const cells = new Float32Array(2 * DAYS * C + perDay * C), days = [];
  for (let day = 0; day < DAYS; day++) {
    for (let step = 0; step < perDay; step++) {
      const time = model.time;
      await model.step(dt);
      const [ts, flux, sensible, evaporation] = await Promise.all([
        readRanges(model.gpu.device, model.gpu.buffers.S, [{ offset: S.TS, length: C }]).then(([x]) => x),
        ...['SFLUX', 'SH', 'EVAP'].map((name) => readRanges(model.gpu.device, model.gpu.buffers.PH, [{ offset: PH[name], length: C }]).then(([x]) => x)),
      ]);
      sunDirection(time, sun);
      const sunLongitude = Math.atan2(sun[1], sun[0]);
      for (let i = 0; i < C; i++) {
        const hour = ((12 + (mesh.lonCell[i] - sunLongitude) * 12 / Math.PI) % 24 + 24) % 24, bin = Math.floor(hour) % HOURS;
        sums[bin * C + i] += ts[i]; counts[bin * C + i]++;
        if (day === 0) cells[2 * DAYS * C + step * C + i] = flux[i] + sensible[i] + lheat * evaporation[i];
      }
    }
    const [absorbedSum, atmosphereSum, longwaveSum] = await readRanges(model.gpu.device, model.gpu.buffers.PH, [{ offset: PH.ABSSUM, length: C }, { offset: PH.ATMSUM, length: C }, { offset: PH.LWSFCSUM, length: C }]);
    let area = 0, surfaceSolar = 0, surfaceLongwave = 0;
    for (let i = 0; i < C; i++) { area += mesh.areaCell[i]; surfaceSolar += mesh.areaCell[i] * (absorbedSum[i] - atmosphereSum[i]) / perDay; surfaceLongwave += mesh.areaCell[i] * longwaveSum[i] / perDay; }
    const d = await model.diagnostics();
    await model.sync();
    cells.set(model.radiation.meanAbsorbedSolar, 2 * day * C);
    cells.set(model.radiation.meanOutgoingLongwave, (2 * day + 1) * C);
    days.push({ asr: d.absorbedSolar, olr: d.outgoingLongwave, albedo: d.planetaryAlbedo, swcre: d.shortwaveCloudEffect, lwcre: d.longwaveCloudEffect, surfaceSolar: surfaceSolar / area, surfaceLongwave: surfaceLongwave / area, rain: 86400 * d.precipitation, ts: d.meanSurfaceT - 273.15 });
    console.log(`day ${day + 1}: ${JSON.stringify(days[day])}`);
  }
  const diurnal = { land: new Float64Array(HOURS), sea: new Float64Array(HOURS) }, weight = { land: 0, sea: 0 };
  for (let i = 0; i < C; i++) {
    let mean = 0, n = 0;
    for (let b = 0; b < HOURS; b++) if (counts[b * C + i]) { mean += sums[b * C + i] / counts[b * C + i]; n++; }
    if (n < HOURS) continue;
    mean /= n;
    const kind = land[i] ? 'land' : 'sea', a = mesh.areaCell[i];
    for (let b = 0; b < HOURS; b++) diurnal[kind][b] += a * (sums[b * C + i] / counts[b * C + i] - mean);
    weight[kind] += a;
  }
  for (const kind of ['land', 'sea']) for (let b = 0; b < HOURS; b++) diurnal[kind][b] /= weight[kind];
  writeFileSync(out, JSON.stringify({ state: stateFile, N: saved.N, C, perDay, days: DAYS, perturb: PERTURB, radiation: process.env.RADIATION ?? '{}', ocean: process.env.OCEAN ?? '{}', means: days, diurnal: { land: Array.from(diurnal.land), sea: Array.from(diurnal.sea) } }, null, 1));
  writeFileSync(`${out}.bin`, new Uint8Array(cells.buffer));
  process.exit(0);
}

if (mode === 'compare') {
  const { Grid } = await import('../js/grid.module.js');
  const { buildMesh } = await import('../js/mesh.module.js');
  const load = (file) => ({ ...JSON.parse(readFileSync(file, 'utf8')), cells: new Float32Array(readFileSync(`${file}.bin`).buffer.slice(0)), file });
  const [reference, ...runs] = files.map(load);
  const mesh = buildMesh(new Grid(reference.N)), C = reference.C, area = mesh.areaCell, total = area.reduce((s, a) => s + a, 0);
  const rms = (a, aOff, b, bOff) => { let s = 0; for (let i = 0; i < C; i++) { const x = a[aOff + i] - b[bOff + i]; s += area[i] * x * x; } return Math.sqrt(s / total); };
  const cycle = (curve) => {
    const top = curve.indexOf(Math.max(...curve)), before = curve[(top + HOURS - 1) % HOURS], after = curve[(top + 1) % HOURS];
    const shift = before - 2 * curve[top] + after < 0 ? 0.5 * (before - after) / (before - 2 * curve[top] + after) : 0;
    return { amplitude: Math.max(...curve) - Math.min(...curve), peak: (top + 0.5 + shift + HOURS) % HOURS };
  };
  const keys = Object.keys(reference.means[0]);
  const show = (run) => {
    const lines = [`${run.file.replace(/.*\//, '')} (RADIATION ${run.radiation}, OCEAN ${run.ocean}${run.perturb ? `, θ perturbed ${run.perturb}` : ''})`];
    for (let d = 0; d < run.days; d++) lines.push(`  day ${d + 1}: ${keys.map((k) => `${k} ${run.means[d][k].toFixed(3)} (${(run.means[d][k] - reference.means[d][k] >= 0 ? '+' : '')}${(run.means[d][k] - reference.means[d][k]).toFixed(3)})`).join(', ')}`);
    const mean = (k) => run.means.reduce((s, m) => s + m[k], 0) / run.days - reference.means.reduce((s, m) => s + m[k], 0) / reference.days;
    lines.push(`  ${run.days}-day means less the reference's: ${keys.map((k) => `${k} ${mean(k) >= 0 ? '+' : ''}${mean(k).toFixed(3)}`).join(', ')}`);
    for (const kind of ['land', 'sea']) { const c = cycle(run.diurnal[kind]), r = cycle(reference.diurnal[kind]); lines.push(`  diurnal Ts over ${kind}: amplitude ${c.amplitude.toFixed(3)} K (reference ${r.amplitude.toFixed(3)}, ${(100 * (c.amplitude / r.amplitude - 1)).toFixed(2)} %), peak at ${c.peak.toFixed(2)} h local (reference ${r.peak.toFixed(2)})`); }
    lines.push(`  per-cell rms of the day-mean ASR / OLR less the reference's: ${Array.from({ length: run.days }, (_, d) => `${rms(run.cells, 2 * d * C, reference.cells, 2 * d * C).toFixed(2)} / ${rms(run.cells, (2 * d + 1) * C, reference.cells, (2 * d + 1) * C).toFixed(2)}`).join(', ')} W/m² by day`);
    const steps = run.perDay, base = 2 * run.days * C;
    let first = 0, all = 0;
    for (let s = 0; s < steps; s++) { const x = rms(run.cells, base + s * C, reference.cells, base + s * C); all += x; if (s < 8) first += x; }
    lines.push(`  per-step rms of the surface net radiation less the reference's: ${(first / 8).toFixed(2)} W/m² over the first 8 steps, ${(all / steps).toFixed(2)} over the first day`);
    return lines.join('\n');
  };
  console.log(`reference ${reference.file.replace(/.*\//, '')}: ${reference.days} days from ${reference.state.replace(/.*\//, '')}; diurnal Ts land ${JSON.stringify(cycle(reference.diurnal.land))}, sea ${JSON.stringify(cycle(reference.diurnal.sea))}`);
  for (const run of runs) console.log(show(run));
  process.exit(0);
}

console.log('node scripts/radiationInterval.mjs run <state.bin> <out.json> | compare <reference.json> <run.json>...');
