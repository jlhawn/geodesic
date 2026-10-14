// The sea's day as the atmosphere steps see it, read after each physics
// pass of a GPU model and summed over each UTC day; it reads only, so a
// run with it steps as the run without. watchSkin(model, { prefix, sea,
// surface, fieldDays }) appends one JSON line a day to <prefix>.skin.jsonl
// over the open water (sea, ice-free cells): the open cell-steps, those
// with the warm layer at its cap, the warmest skin outside the Persian
// Gulf, the Gulf of Oman and the Red Sea, the warm layer at local 03 h,
// the largest C_H over its neutral value, the same warm layer at 03 h in the
// cells that are dark from local 00 h and lose heat over those hours (dark,
// darkWarm, darkOver01, darkOver1, darkWarmest, darkAt; scripts/nightCensus.mjs),
// and sparse histograms of C_H,
// the sensible plus latent heat and the physics' step of the lowest
// layer's θ. On the days in fieldDays it writes <prefix>.skin_dayD.f32,
// the per-cell fields FIELDS one after another as float32 (C each).
// Before the day `from` it reads nothing. Local hours are solar: model
// time 0 is noon at longitude 0 (radiation.module.js's sunDirection).
import { appendFileSync, writeFileSync } from 'node:fs';
import { readRanges } from '../js/gpu/device.module.js';
import { skinOptions, SEA_RHO_CP } from '../js/physics/skin.module.js';
import { LATENT_HEAT } from '../js/physics/moist.module.js';
import { createNightCensus } from './nightCensus.mjs';

export const FIELDS = ['skinMax', 'skinMin', 'skinHourMax', 'skinHourMin', 'warmMax', 'coolMean', 'dawnH0', 'wind', 'insolation', 'slabMean', 'open'];
export const HISTOGRAMS = { exchange: { from: 0, width: 1e-5 }, ratio: { from: 0, width: 0.01 }, heat: { from: -1000, width: 1 }, theta: { from: 0, width: 1e-3 } };
const DAY = 86400, deg = 180 / Math.PI;
const GULFS = [[23.5, 30.5, 47.5, 56.5], [22, 26.8, 56.5, 61.5], [12, 30, 32, 43.6]];

export function watchSkin(model, { prefix, sea, surface = {}, fieldDays = [], from = 0 }) {
  const { mesh, gpu } = model, C = mesh.nCells, K = model.core.K, { S, PH, D } = gpu.layout, device = gpu.device;
  const o = skinOptions(surface), skinned = PH.WARM !== undefined;
  const lat = Float64Array.from(mesh.latCell, (x) => x * deg), lon = Float64Array.from(mesh.lonCell, (x) => ((x * deg + 540) % 360) - 180);
  const gulf = Uint8Array.from({ length: C }, (_, i) => (GULFS.some(([s, n, w, e]) => lat[i] >= s && lat[i] <= n && lon[i] >= w && lon[i] <= e) ? 1 : 0));
  const want = new Set(fieldDays.map(Number));
  const f = Object.fromEntries(FIELDS.map((name) => [name, new Float64Array(C)]));
  const census = createNightCensus(C);
  let steps = 0, totals = null, before = null;
  const histogram = () => new Map();
  const add = (map, name, x) => { const { from, width } = HISTOGRAMS[name], b = Math.max(0, Math.floor((x - from) / width)); map.set(b, (map.get(b) ?? 0) + 1); };
  function reset() {
    f.skinMax.fill(-Infinity); f.skinMin.fill(Infinity); f.skinHourMax.fill(NaN); f.skinHourMin.fill(NaN); f.warmMax.fill(0); f.coolMean.fill(0);
    f.dawnH0.fill(NaN); f.wind.fill(0); f.insolation.fill(0); f.slabMean.fill(0); f.open.fill(0);
    steps = 0;
    totals = { open: 0, capped: 0, hottest: -Infinity, hottestAt: null, night: 0, nightWarm: 0, nightWarmest: 0, nightAt: null, nightOver01: 0, ratioMax: 0, ratioAt: null, hist: { exchange: histogram(), ratio: histogram(), heat: histogram(), theta: histogram() } };
  }
  reset();
  const ocean = gpu.hooks.beforePhysics;
  const idle = () => model.time < from * DAY - 1;
  gpu.hooks.beforePhysics = async (dt, n) => {
    if (ocean) await ocean(dt, n);
    if (idle()) return;
    [before] = await readRanges(device, gpu.buffers.S, [{ offset: S.TH + (K - 1) * C, length: C }]);
  };
  const later = gpu.hooks.afterPhysics;
  gpu.hooks.afterPhysics = async (dt, n) => {
    if (later) await later(dt, n);
    if (idle()) return;
    const time = model.time;
    const [[ts, ice, theta], ph, [wind]] = await Promise.all([
      readRanges(device, gpu.buffers.S, [{ offset: S.TS, length: C }, { offset: S.ICE, length: C }, { offset: S.TH + (K - 1) * C, length: C }]),
      readRanges(device, gpu.buffers.PH, [PH.CAP, PH.HEATX, PH.REFX, PH.SH, PH.EVAP, PH.INS, ...(skinned ? [PH.WARM, PH.COOL] : []), PH.SFLUX].map((offset) => ({ offset, length: C }))),
      readRanges(device, gpu.buffers.D, [{ offset: D.WIND, length: C }]),
    ]);
    const [cap, heatx, refx, sensible, evap, ins] = ph, warm = skinned ? ph[6] : null, cool = skinned ? ph[7] : null, sflux = ph[ph.length - 1];
    const solarHour = 12 + ((time + 0.5 * dt) % DAY) / 3600;
    steps++;
    for (let i = 0; i < C; i++) {
      if (!sea(i)) continue;
      f.wind[i] += wind[i]; f.insolation[i] += ins[i];
      if (ice[i] > 0) continue;
      const hour = (((solarHour + lon[i] / 15) % 24) + 24) % 24, h0 = cap[i] / SEA_RHO_CP;
      const w = warm ? warm[i] : 0, c = cool ? cool[i] : 0, d = Math.min(o.warmLayerDepth, h0);
      const skin = ts[i] + (o.warmLayer ? (1 - o.warmLayerProfile * d / ((o.warmLayerProfile + 1) * h0)) * w : 0) - c - 273.15;
      f.open[i]++; f.slabMean[i] += ts[i] - 273.15; f.coolMean[i] += c;
      census.observe(i, hour, ins[i], sflux[i], dt, w);
      if (skin > f.skinMax[i]) { f.skinMax[i] = skin; f.skinHourMax[i] = hour; }
      if (skin < f.skinMin[i]) { f.skinMin[i] = skin; f.skinHourMin[i] = hour; }
      if (w > f.warmMax[i]) f.warmMax[i] = w;
      if (hour >= 5.5 && hour < 6.5) f.dawnH0[i] = h0;
      totals.open++;
      if (o.warmLayer && w >= o.warmLayerCap - 1e-4) totals.capped++;
      if (!gulf[i] && skin > totals.hottest) { totals.hottest = skin; totals.hottestAt = [+lat[i].toFixed(1), +lon[i].toFixed(1)]; }
      if (hour >= 2.5 && hour < 3.5) {
        totals.night++;
        if (w > 0) { totals.nightWarm++; if (w > 0.1) totals.nightOver01++; if (w > totals.nightWarmest) { totals.nightWarmest = w; totals.nightAt = [+lat[i].toFixed(1), +lon[i].toFixed(1), +sensible[i].toFixed(1)]; } }
      }
      const ratio = refx[i] > 0 ? heatx[i] / refx[i] : 0;
      if (ratio > totals.ratioMax) { totals.ratioMax = ratio; totals.ratioAt = [+lat[i].toFixed(1), +lon[i].toFixed(1), +skin.toFixed(2), +hour.toFixed(1)]; }
      add(totals.hist.exchange, 'exchange', heatx[i]); add(totals.hist.ratio, 'ratio', ratio);
      add(totals.hist.heat, 'heat', sensible[i] + LATENT_HEAT * evap[i]); add(totals.hist.theta, 'theta', Math.abs(theta[i] - before[i]));
    }
    const ending = Math.round((time + dt) / DAY);
    if (Math.abs(time + dt - ending * DAY) < 0.5 * dt && steps > 0) close(ending);
  };
  function close(day) {
    for (let i = 0; i < C; i++) {
      f.wind[i] /= steps; f.insolation[i] /= steps;
      if (f.open[i] > 0) { f.slabMean[i] /= f.open[i]; f.coolMean[i] /= f.open[i]; } else { f.skinMax[i] = NaN; f.skinMin[i] = NaN; f.slabMean[i] = NaN; f.coolMean[i] = NaN; }
      f.open[i] /= steps;
    }
    const sparse = (map) => [...map.entries()].sort((a, b) => a[0] - b[0]);
    const { hist, ...rest } = totals, dark = census.take(), { darkCell, ...darkCounts } = dark;
    appendFileSync(`${prefix}.skin.jsonl`, JSON.stringify({ day, steps, ...rest, ...darkCounts, darkAt: darkCell >= 0 ? [+lat[darkCell].toFixed(1), +lon[darkCell].toFixed(1)] : null, capShare: totals.open ? totals.capped / totals.open : 0, hist: Object.fromEntries(Object.entries(hist).map(([k, m]) => [k, sparse(m)])) }) + '\n');
    if (want.has(day)) {
      const out = new Float32Array(FIELDS.length * C);
      FIELDS.forEach((name, k) => out.set(f[name], k * C));
      writeFileSync(`${prefix}.skin_day${day}.f32`, out);
    }
    reset();
  }
}
