import { test } from 'node:test';
import assert from 'node:assert/strict';
import { plumeColumn, SPACINGS } from '../explain/plumeColumn.module.js';
import { readPlume, heightAt } from '../explain/plumePath.module.js';
import { LATENT_HEAT } from '../js/physics/moist.module.js';
import { GRAVITY } from '../js/dynamics/sigmaCore.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const CASES = [{}, { humidity: 0.4 }, { humidity: 0.6 }, { humidity: 0.9 }, { humidity: 1.1 }, { warming: 0.5 }, { warming: 1 }, { ice: false }, { spacing: SPACINGS[1] }, { ice: false, humidity: 0.7, warming: 0.3 }];

test('the shaded excess of the plume over the air adds up to the model\'s CAPE', () => {
  for (const c of CASES) {
    const r = plumeColumn(c), read = readPlume(r), what = JSON.stringify(c);
    const condensed = r.heating.largeScale.some((h) => h !== 0);
    close(read.cape, r.cape, (condensed ? 1e-4 : 1e-9) * r.cape, `${what} CAPE`);
    for (let k = 0; k < r.pressure.length; k++) assert.equal(Number.isFinite(read.excess[k]), Number.isFinite(r.plume.buoyancy[k]), `${what} layer ${k}`);
  }
  const r = plumeColumn(), { excess } = readPlume(r);
  const most = Math.max(...excess.filter(Number.isFinite)), least = Math.min(...excess.filter(Number.isFinite));
  console.log(`Jordan: the plume is up to ${most.toFixed(2)} K warmer than the air and, overshooting at its top, ${(-least).toFixed(2)} K colder`);
  close(most, 1.14, 0.01, 'warmest');
});

test('the source air leaves the lowest 50 hPa dry and saturates between the plume\'s base interface and the one below it', () => {
  for (const c of CASES) {
    const r = plumeColumn(c), { path, source } = readPlume(r), what = JSON.stringify(c), K = r.pressure.length;
    const base = r.interfaces.findIndex((p) => Math.abs(p - r.base) < 1e-6);
    assert.ok(base > 0 && base < K, `${what} base interface`);
    assert.ok(source.lcl.pressure >= r.base && source.lcl.pressure < r.interfaces[base + 1], `${what} LCL ${source.lcl.pressure} between ${r.base} and ${r.interfaces[base + 1]}`);
    close(path[0].pressure, r.pressure[K - 1], 0, `${what} the path starts in the lowest layer`);
    for (let n = 1; n < path.length; n++) assert.ok(path[n].pressure < path[n - 1].pressure, `${what} the path rises`);
    const atLcl = path.findIndex((point) => point.pressure === source.lcl.pressure);
    close(path[atLcl].temperature, source.lcl.temperature, 1e-9, `${what} the dry leg ends at the LCL's temperature`);
    const highest = r.interfaces.findIndex((p, k) => Number.isFinite(r.plume.interfaceTemperature[k]));
    close(path.at(-1).pressure, r.interfaces[highest], 0, `${what} the path ends at the plume's highest interface`);
    const top = r.pressure.findIndex((p) => p >= r.surfacePressure - 50e2);
    close(source.mass * GRAVITY, r.surfacePressure - r.interfaces[top], 1e-9 * r.surfacePressure, `${what} the source is every layer centered within 50 hPa of the ground`);
  }
  const { source } = readPlume(plumeColumn());
  console.log(`Jordan's source air: ${source.temperature.toFixed(2)} K and ${(source.vapor * 1e3).toFixed(2)} g/kg at ${(source.pressure / 100).toFixed(1)} hPa, saturating at ${(source.lcl.pressure / 100).toFixed(1)} hPa and ${source.lcl.temperature.toFixed(2)} K`);
});

test('the heating and the drying each add up through the column to the latent heat of the rain and of the cloud water left behind', () => {
  for (const c of CASES) {
    const r = plumeColumn(c), { column } = readPlume(r), what = JSON.stringify(c);
    let cloud = 0;
    for (let k = 0; k < r.pressure.length; k++) cloud += r.mass[k] * r.cloudWater[k];
    const left = LATENT_HEAT * cloud / r.dt;
    close(column.heating, column.drying, 1e-9 * column.heating, `${what} heating against drying`);
    close(column.heating, column.rain + left, 1e-9 * column.heating, `${what} heating against rain`);
  }
  const { column, peak } = readPlume(plumeColumn());
  console.log(`Jordan: the column gains ${column.heating.toFixed(1)} W/m², the latent heat of ${(column.rain / LATENT_HEAT * 86400).toFixed(2)} mm/day; heating peaks at ${peak.rate.toFixed(2)} K/day at ${(peak.pressure / 100).toFixed(0)} hPa`);
  close(column.rain, 1159.2, 0.1, 'the latent heat of 40.06 mm/day');
  close(peak.rate, 19.28, 0.01, 'peak heating');
  close(peak.pressure / 100, 440.4, 0.1, 'at 440 hPa');
});

test('the base mass flux stops at a tenth of the source air per step', () => {
  for (const c of [{}, { humidity: 1.2 }, { warming: 3 }, { humidity: 0.4, spacing: SPACINGS[1] }]) {
    const r = plumeColumn(c), { cap } = readPlume(r), what = JSON.stringify(c);
    assert.ok(r.baseFlux <= cap * (1 + 1e-12), `${what} ${r.baseFlux} within ${cap}`);
    if (c.humidity === 1.2 || c.warming === 3) close(r.baseFlux, cap, 1e-12 * cap, `${what} at the cap`);
  }
});

test('heights interpolate between the layers in the logarithm of pressure', () => {
  const r = plumeColumn(), K = r.pressure.length;
  for (let k = 0; k < K; k++) close(heightAt(r, r.pressure[k]), r.height[k], 1e-9, `layer ${k}`);
  close(heightAt(r, r.surfacePressure), 0, 1e-9, 'the ground');
  for (let p = r.surfacePressure; p > 2e3; p *= 0.97) assert.ok(heightAt(r, p * 0.97) > heightAt(r, p), `rises above ${p}`);
  const { topHeight } = readPlume(r);
  console.log(`Jordan's plume top: ${(r.top / 100).toFixed(0)} hPa, ${(topHeight / 1000).toFixed(2)} km`);
  assert.ok(topHeight > r.height[r.pressure.findIndex((p) => p > r.top)] && topHeight < r.height[r.pressure.findLastIndex((p) => p < r.top)], 'the top lies between its neighbors');
});
