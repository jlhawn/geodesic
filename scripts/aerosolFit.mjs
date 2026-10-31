// The spectral constants of radiation option aerosol 'mac' (js/physics/aerosolSpectra.module.js) from MACv2's band files:
//   node scripts/aerosolFit.mjs DIR [DIR...]            the fit and its misses against the band data
//   WRITE=1 node scripts/aerosolFit.mjs DIR [DIR...]    also rewrites js/physics/aerosolSpectra.module.js
//
// MACv2 (Kinne 2019, Tellus B 71, 1623639) gives each mode's AOD, single scattering albedo and asymmetry at 550 nm and at
// the 14 solar wavelengths of ICON's RRTM_SW bands, the centers of RRTMG_SW's bands 16-26 and 29 (ftp://ftp-projects.mpimet.mpg.de/
// aerocom/climatology/MACv2_2018/for_ICON/spec_2005/by_wave/gt_{f,c}_WWWWWnm.nc, the fine and coarse mode at WWWWW nm). The page
// carries the 550 nm triple of each mode and each cell, and this script the five constants r, a, c, d, e per mode and
// wavelength with which the triple (tau, x = 1 - omega, g) becomes the wavelength's:
//     tau_b = r tau,   x_b = a + c x,   g_b = d + e g
// fitted to August (the data of DIR: a file gt_{f,c}_WWWWWnm.nc whole, or its August record gt_{f,c}_WWWWWnm.aug.bin cut from
// the file by byte range, in DIR or in a directory in it) by least squares weighted by each cell's 550 nm AOD times its area; r is
// the weighted mean of the ratio of AOD. Beyond 2.5 micrometer (bands 17, 16 and 29, 3 % of the beam) the albedo and asymmetry do
// not follow the 550 nm values and the band means stand in (c = e = 0).
// The Planck weights (5778 K) of the wavelengths' ranges in the two streams of the sequential scheme (0.297-0.7 um, the rest) and
// in the layered scheme's visible slices (0.297-0.625 um) are written with them.
import { readFileSync, readdirSync, existsSync, statSync, writeFileSync } from 'node:fs';
import { netcdfClassic } from './netcdfClassic.mjs';

const WAVELENGTHS = [400, 550, 700, 1000, 1270, 1460, 1780, 2050, 2320, 2790, 3470, 8000];
const BANDS = [26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 29], MEANS_FROM = 9;
const VISIBLE_RANGES = { 400: [0.297, 0.441], 550: [0.441, 0.625], 700: [0.625, 0.7] };
const NEAR_INFRARED_RANGES = { 700: [0.7, 0.778], 1000: [0.778, 1.24], 1270: [1.24, 1.30], 1460: [1.30, 1.63], 1780: [1.63, 1.94], 2050: [1.94, 2.15], 2320: [2.15, 2.50], 2790: [2.50, 3.08], 3470: [3.08, 3.85], 8000: [3.85, 1000] };
const SLICE_RANGES = { 400: [0.297, 0.441], 550: [0.441, 0.625] };
const AUGUST = 7, ROWS = 180, COLS = 360, CELLS = ROWS * COLS, TEMPERATURE = 5778;
const REGIONS = { global: [-90, 90, -180, 180], 'tropics 10S-10N': [-10, 10, -180, 180], 'midlatitudes 40-60N': [40, 60, -180, 180], India: [8, 28, 70, 90], 'Sahara and Arabia': [15, 32, -17, 60], 'East Asia': [22, 45, 100, 125] };

const pad = (nm) => String(nm).padStart(5, '0');
function locate(dirs, mode, nm) {
  for (const dir of dirs) {
    const places = [dir, ...readdirSync(dir).map((name) => `${dir}/${name}`).filter((path) => statSync(path).isDirectory())];
    for (const place of places) for (const suffix of ['.nc', '.aug.bin']) {
      const path = `${place}/gt_${mode}_${pad(nm)}nm${suffix}`;
      if (existsSync(path)) return path;
    }
  }
  throw new Error(`gt_${mode}_${pad(nm)}nm.nc or .aug.bin is in none of ${dirs.join(', ')}`);
}

function load(path) {
  if (path.endsWith('.nc')) {
    const file = netcdfClassic(path), out = ['aod', 'ssa', 'asy'].map((name) => file.read(name, AUGUST));
    file.close();
    return out;
  }
  const bytes = readFileSync(path);
  if (bytes.length !== 4 + 12 * CELLS) throw new Error(`${path}: ${bytes.length} bytes, not the 4 + 3 x ${CELLS} floats of an August record`);
  return [0, 1, 2].map((v) => Float64Array.from({ length: CELLS }, (_, c) => bytes.readFloatBE(4 + 4 * (v * CELLS + c))));
}

const latitude = (row) => 89.5 - row, longitude = (col) => -179.5 + col;
const area = Float64Array.from({ length: CELLS }, (_, c) => Math.cos(latitude(Math.floor(c / COLS)) * Math.PI / 180));

function weightedFit(x, y, weight) {
  let sw = 0, mx = 0, my = 0;
  for (let c = 0; c < CELLS; c++) { sw += weight[c]; mx += weight[c] * x[c]; my += weight[c] * y[c]; }
  mx /= sw; my /= sw;
  let sxx = 0, sxy = 0;
  for (let c = 0; c < CELLS; c++) { sxx += weight[c] * (x[c] - mx) ** 2; sxy += weight[c] * (x[c] - mx) * (y[c] - my); }
  const slope = sxx > 0 ? sxy / sxx : 0, offset = my - slope * mx;
  let rss = 0;
  for (let c = 0; c < CELLS; c++) rss += weight[c] * (y[c] - offset - slope * x[c]) ** 2;
  return { offset, slope, mean: my, residual: Math.sqrt(rss / sw) };
}

function planck(wavelength) {
  const h = 6.62607015e-34, c = 299792458, k = 1.380649e-23, l = wavelength * 1e-6;
  return 1 / (l ** 5 * (Math.exp(h * c / (l * k * TEMPERATURE)) - 1));
}
function integral(a, b, n = 20000) {
  let s = 0;
  for (let i = 0; i < n; i++) { const x0 = a + (b - a) * i / n, x1 = a + (b - a) * (i + 1) / n; s += (planck(x0) + 4 * planck(0.5 * (x0 + x1)) + planck(x1)) * (x1 - x0) / 6; }
  return s;
}
function weights(ranges) {
  const parts = WAVELENGTHS.map((nm) => (ranges[nm] ? integral(...ranges[nm]) : 0)), total = parts.reduce((a, b) => a + b, 0);
  return parts.map((p) => p / total);
}

const six = (x) => Number(x.toPrecision(6));
const dirs = process.argv.slice(2);
if (!dirs.length) throw new Error('usage: node scripts/aerosolFit.mjs DIR [DIR...]');
const data = { f: {}, c: {} };
for (const mode of ['f', 'c']) for (const nm of WAVELENGTHS) data[mode][nm] = load(locate(dirs, mode, nm));

const constants = {};
for (const mode of ['f', 'c']) {
  const [tau550, ssa550, asy550] = data[mode][550], weight = Float64Array.from(tau550, (t, c) => t * area[c]);
  const x550 = Float64Array.from(ssa550, (s) => 1 - s);
  const total = weight.reduce((a, b) => a + b, 0);
  constants[mode] = { r: [], a: [], c: [], d: [], e: [] };
  console.log(`${mode === 'f' ? 'fine' : 'coarse'} mode (weighted by 550 nm AOD times area)`);
  WAVELENGTHS.forEach((nm, b) => {
    const [tau, ssa, asy] = data[mode][nm];
    let mean = 0, spread = 0;
    for (let c = 0; c < CELLS; c++) mean += weight[c] * tau[c] / Math.max(tau550[c], 1e-30);
    mean /= total;
    for (let c = 0; c < CELLS; c++) spread += weight[c] * (tau[c] / Math.max(tau550[c], 1e-30) - mean) ** 2;
    const absorption = weightedFit(x550, Float64Array.from(ssa, (s) => 1 - s), weight), asymmetry = weightedFit(asy550, asy, weight);
    const means = b >= MEANS_FROM, fit = nm === 550 ? [0, 1, 0, 1] : means ? [absorption.mean, 0, asymmetry.mean, 0] : [absorption.offset, absorption.slope, asymmetry.offset, asymmetry.slope];
    constants[mode].r.push(nm === 550 ? 1 : six(mean));
    ['a', 'c', 'd', 'e'].forEach((key, n) => constants[mode][key].push(six(fit[n])));
    console.log(`  ${String(nm).padStart(4)} nm (band ${BANDS[b]}): AOD ratio ${mean.toFixed(3)} (sd ${Math.sqrt(spread / total).toFixed(3)}); 1-SSA = ${fit[0].toFixed(4)} + ${fit[1].toFixed(3)} x (miss ${means ? 'band mean' : absorption.residual.toFixed(4)}); ASY = ${fit[2].toFixed(3)} + ${fit[3].toFixed(3)} g (miss ${means ? 'band mean' : asymmetry.residual.toFixed(4)})`);
  });
}

const spectra = {
  wavelengths: WAVELENGTHS, bands: BANDS, fine: constants.f, coarse: constants.c,
  visible: weights(VISIBLE_RANGES).map(six), nearInfrared: weights(NEAR_INFRARED_RANGES).map(six), slices: weights(SLICE_RANGES).map(six),
};

function modeBand(mode, b, tau, x, g) {
  const s = mode === 0 ? spectra.fine : spectra.coarse;
  return [s.r[b] * tau, Math.min(1, Math.max(0, s.a[b] + s.c[b] * x)), Math.min(0.99, Math.max(0, s.d[b] + s.e[b] * g))];
}
function regionMeans() {
  const rows = [];
  for (const [name, [south, north, west, east]] of Object.entries(REGIONS)) {
    const sums = { band: [0, 0, 0, 0], fit: [0, 0, 0, 0] };
    let total = 0;
    for (let c = 0; c < CELLS; c++) {
      const lat = latitude(Math.floor(c / COLS)), lon = longitude(c % COLS);
      if (lat < south || lat > north || lon < west || lon > east) continue;
      total += area[c];
      for (const kind of ['band', 'fit']) {
        const out = sums[kind];
        for (let mode = 0; mode < 2; mode++) {
          const m = mode === 0 ? 'f' : 'c', [tau550, ssa550, asy550] = data[m][550];
          WAVELENGTHS.forEach((nm, b) => {
            let t, x, g;
            if (kind === 'band') { t = data[m][nm][0][c]; x = 1 - data[m][nm][1][c]; g = data[m][nm][2][c]; } else [t, x, g] = modeBand(mode, b, tau550[c], 1 - ssa550[c], asy550[c]);
            const wv = spectra.visible[b], wn = spectra.nearInfrared[b];
            out[0] += area[c] * wv * t; out[1] += area[c] * wv * t * x; out[2] += area[c] * wn * t; out[3] += area[c] * wn * t * x;
          });
        }
      }
    }
    rows.push({ name, band: sums.band.map((v) => v / total), fit: sums.fit.map((v) => v / total) });
  }
  return rows;
}
console.log('\nstreams\' effective depths over each region, area means, from the band data and from the constants:');
console.log('region'.padEnd(22), 'visible tau', 'visible abs', 'near-infrared tau', 'near-infrared abs');
const rows = regionMeans();
for (const { name, band, fit } of rows) console.log(name.padEnd(22), band.map((v, n) => `${v.toFixed(4)} (${(100 * (fit[n] / v - 1)).toFixed(1)}%)`.padEnd(n < 2 ? 18 : 22)).join(''));
console.log('fixture of the band data:', JSON.stringify(Object.fromEntries(rows.map(({ name, band }) => [name, band.map((v) => Number(v.toPrecision(4)))]))));
console.log('Planck weights, visible:', spectra.visible.join(' '), '\n near infrared:', spectra.nearInfrared.join(' '), '\n slices:', spectra.slices.join(' '));

if (process.env.WRITE) {
  const header = `// Generated by scripts/aerosolFit.mjs from MACv2's fine and coarse mode (Kinne 2019, Tellus B 71, 1623639) at its 12 solar wavelengths (nm), the centers of RRTMG_SW's
// bands (\`bands\`): for each mode and wavelength the constants r, a, c, d, e that take the mode's 550 nm AOD tau, 1 - SSA x and asymmetry g to the wavelength's, tau_b = r tau,
// x_b = a + c x, g_b = d + e g, and the Planck weights (5778 K) of the wavelengths in the sequential scheme's visible stream (0.297-0.7 um), its near-infrared stream (0.7 um on)
// and the layered scheme's visible slices (0.297-0.625 um).
`;
  const line = (key, value) => `  ${key}: ${JSON.stringify(value).replaceAll(',', ', ')},\n`;
  const mode = (m) => `{ ${['r', 'a', 'c', 'd', 'e'].map((k) => `${k}: ${JSON.stringify(m[k]).replaceAll(',', ', ')}`).join(', ')} }`;
  const text = `${header}export const AEROSOL_SPECTRA = {\n${line('wavelengths', spectra.wavelengths)}${line('bands', spectra.bands)}  fine: ${mode(spectra.fine)},\n  coarse: ${mode(spectra.coarse)},\n${line('visible', spectra.visible)}${line('nearInfrared', spectra.nearInfrared)}${line('slices', spectra.slices)}};\n`;
  writeFileSync(new URL('../js/physics/aerosolSpectra.module.js', import.meta.url), text);
  console.log('wrote js/physics/aerosolSpectra.module.js');
}
