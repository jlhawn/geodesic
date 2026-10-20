// The near-infrared gas absorption of js/physics/shortwaveGases.module.js's 'bands' scheme, fitted
// to RRTMG_SW, and the visible water vapor strength that goes with it.
//   RRTMG_SW=<rrtmg_sw> RRTMG_INPUTS=<dir> STATE=<state.bin> node scripts/nearInfraredFit.mjs reference <ref.json>
//   node scripts/nearInfraredFit.mjs fit <ref.json> [<fit.json>]     WRITE=1 also rewrites js/physics/nearInfraredTable.module.js
//   node scripts/nearInfraredFit.mjs report <ref.json> [<fit.json>]  the table (or a fit) against the reference
//
// The reference: AER's RRTMG_SW column model (github.com/AER-RC/RRTMG_SW, its rrtmg_sw.1col.f90 built
// with gfortran and the two flux tables' write statements given the format
// '(1x,i3,1x,es14.7,6(1x,es17.10))' for their precision), run on AER's four standard atmospheres
// (RRTMG_INPUTS: the repository's run_examples_std_atm), the five AFGL atmospheres of
// data/radiationBenchmark.json and the sixteen TARGETS columns of a saved model state (STATE; the
// table's was spin-up thirteen's day 1825 at N=64), each at cosines MUS (1, 0.5, 0.25) of the sun's
// zenith angle with no gas, with water vapor, CO2, CH4 and O2 each alone, and with the four together.
// Every band's direct-beam transmission to each level, the run's direct flux over the gas-free run's
// (which takes the Rayleigh extinction out exactly: it is one coefficient per band), is what is fitted.
//
// The scheme: in each of RRTMG_SW's bands 16-24 and 29 (820-16000 cm-1), the share of the beam the band
// carries times a transmission that is the product of its gases' (random overlap), each gas's
//   T = w0 exp(-m s) + sum_j w_j exp(-m (tau_j + s)),
// tau_j summed down the layers above, a layer adding k_j u (p_V/p_ref)^n_j exp(t_j (T - 250)) for its
// mass u (kg/m2) of the gas, p_V = sqrt(p^2 + p_D^2) the Voigt half-width as a pressure (p_D: the
// pressure at which the gas's Lorentz half-width at the band's line center, CENTERS, equals its Doppler
// half-width at 250 K), and s, water vapor's self continuum, summed down as self u e/p_ref (296/T)^selfT
// with e the vapor pressure; m the path's magnification. Each gas is fitted alone (TERMS terms,
// Levenberg-Marquardt from six starts), then all of a band's gases together against the four-gas run
// with the lone runs kept in at SOLO_WEIGHT: random overlap of RRTMG_SW's lone transmissions takes
// 1.3-4.5 W/m2 too much from a column. The misses are the layers' absorption, each over the
// layer's share of the four-gas heating plus 0.3 % of the column's absorption per unit of the layer's
// mass share, and the transmission to the surface over 0.2 % of the column's absorption.
//
// The visible water vapor (CLIRAD-SW's single coefficient over 0.4-0.7 um) keeps its coefficient,
// times the strength fitted here to RRTMG_SW's water vapor alone in bands 25-28, since the scheme takes
// all of band 24 (12850-16000 cm-1), the 625-700 nm it shares with the visible band included.
import { readFileSync, writeFileSync, mkdirSync, rmSync, mkdtempSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';

const REPO = fileURLToPath(new URL('..', import.meta.url));
const MOL = ['h2o', 'co2', 'o3', 'n2o', 'co', 'ch4', 'o2'];
const MASS = { h2o: 18.015, co2: 44.01, o3: 47.997, n2o: 44.013, co: 28.01, ch4: 16.04, o2: 31.999, air: 28.964 };
const AVOGADRO = 6.02214e23, G = 9.806, CP = 1004.64, P_REF = 50000, EPS = 1e-6;
const GAMMA = { h2o: 0.08, co2: 0.07, ch4: 0.06, o2: 0.045 };
const CENTERS = {
  16: { h2o: 3100, ch4: 3020 }, 17: { h2o: 3700, co2: 3650 }, 18: { h2o: 4300, ch4: 4300 }, 19: { h2o: 4900, co2: 4980 }, 20: { h2o: 5330, ch4: 5900 },
  21: { h2o: 7200, co2: 6350 }, 22: { h2o: 7900, o2: 7880 }, 23: { h2o: 10600 }, 24: { h2o: 13900, o2: 13120 }, 29: { h2o: 1595, co2: 2350 },
};
const BAND_OF_BLOCK = [null, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29];
const NIR_BLOCKS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 14], VISIBLE_BLOCKS = [10, 11, 12, 13];
const CONFIGS = { none: [], h2o: ['h2o'], co2: ['co2'], ch4: ['ch4'], o2: ['o2'], nir: ['h2o', 'co2', 'ch4', 'o2'] };
const TARGETS = [['warm pool', 0, 160], ['cold tongue', 0, -100], ['SE Pacific', -20, -90], ['Sahara', 25, 5], ['Amazon', -3, -60], ['Tibet', 32, 88], ['N Pacific', 45, -160],
  ['Southern Ocean', -55, 0], ['Greenland', 74, -40], ['Antarctic plateau', -80, 90], ['Arctic', 85, 0], ['Great Plains', 40, -100], ['Siberia', 60, 100], ['W Pacific 15N', 15, 140], ['S Atlantic', -35, -10], ['Arabian Sea', 12, 65]];

export function dopplerPressure(gas, nu, T = 250) {
  const vth = Math.sqrt(2 * Math.LN2 * 1.380649e-23 * T / (MASS[gas] * 1.66054e-27));
  return 101325 * nu * vth / 2.99792458e8 / (GAMMA[gas] * Math.pow(296 / T, 0.75));
}

// ---- the reference ----
const dryColumn = (dp, r) => dp / (G * (MASS.air + r * MASS.h2o) / AVOGADRO / 1000) / 1e4;
const pf = (p) => { const h = p / 100; return (h >= 1000 ? h.toFixed(3) : h >= 100 ? h.toFixed(4) : h >= 10 ? h.toFixed(5) : h.toFixed(6)).padStart(8); };
const e15 = (x) => x.toExponential(7).toUpperCase().padStart(15);

// Layers bottom first: { pBot, pTop, p (Pa), T, vmr per dry air, coldry (molecules/cm2) }
function parseInput(text) {
  const lines = text.split('\n'), start = lines.findIndex((l) => l.startsWith('$'));
  const n = Number(lines[start + 4].slice(2, 5)), layers = [];
  let pz = Number(lines[start + 5].slice(48, 56));
  for (let l = 0; l < n; l++) {
    const a = lines[start + 5 + 2 * l], b = lines[start + 6 + 2 * l];
    const top = Number(a.slice(70, 78)), w = Array.from({ length: 8 }, (_, f) => Number(b.slice(15 * f, 15 * f + 15)));
    layers.push({ pave: Number(a.slice(0, 15)), tave: Number(a.slice(15, 25)), pBot: pz * 100, pTop: top * 100, w });
    pz = top;
  }
  const mixing = layers[0].w.slice(0, 7).every((x) => x <= 1);
  return layers.map(({ pave, tave, pBot, pTop, w }) => {
    const others = w.slice(1, 7).reduce((s, x) => s + x, 0), coldry = mixing ? w[7] / (1 - others) : w[7] + others;
    return { pBot, pTop, p: pave * 100, T: tave, vmr: Object.fromEntries(MOL.map((m, j) => [m, mixing ? w[j] : w[j] / coldry])), coldry };
  });
}

function writeInput(layers, sza, keep) {
  const out = ['$ fitted column', ' HI=1 F4=1 CN=1 AE 0 EM=0 SC=0 FI=0 PL=0 TS=0 AM=0 MG=1 LA=0    1        00   00  1 0   98   00   00',
    `${' '.repeat(12)}${'0'.padStart(3)}${' '.repeat(3)}${sza.toFixed(4).padStart(7)}`, `${' '.repeat(11)}1  0${(0.8).toFixed(3).padStart(5)}`,
    ` 1${String(layers.length).padStart(3)}${'7'.padStart(5)}  1.000000`];
  layers.forEach((L, l) => {
    const head = `${e15(L.p / 100)}${L.T.toFixed(4).padStart(10)}${' '.repeat(15)}`;
    out.push(l === 0 ? `${head} ${' '.repeat(7)}${pf(L.pBot)}${L.T.toFixed(2).padStart(7)}${' '.repeat(7)}${pf(L.pTop)}${L.T.toFixed(2).padStart(7)}` : `${head}${' '.repeat(30)}${pf(L.pTop)}${L.T.toFixed(2).padStart(7)}`);
    const v = MOL.map((m) => (keep.includes(m) ? L.vmr[m] : 0)), others = v.slice(1).reduce((s, x) => s + x, 0);
    out.push([...v, L.coldry * (1 - others)].map(e15).join(''));
  });
  out.push('%%%%%', '');
  return out.join('\n');
}

function parseOutput(text) {
  const bands = [];
  let band = null;
  for (const line of text.split('\n')) {
    const m = line.match(/Wavenumbers:\s*([\d.]+)\s*-\s*([\d.]+)\s*cm-1/);
    if (m) { band = { range: [Number(m[1]), Number(m[2])], levels: [] }; bands.push(band); continue; }
    if (!band) continue;
    const v = line.trim().split(/\s+/).map((t) => Number(t.replace(/(\d)([+-]\d{3})$/, '$1E$2')));
    if (v.length === 8 && v.every(Number.isFinite)) band.levels.push({ p: v[1] * 100, dirdown: v[4], down: v[5] });
  }
  for (const b of bands) b.levels.sort((a, c) => a.p - c.p);
  return bands;
}

function runRrtmg(text) {
  const dir = mkdtempSync(join(tmpdir(), 'rrtmg-'));
  writeFileSync(join(dir, 'INPUT_RRTM'), text);
  execFileSync(process.env.RRTMG_SW, { cwd: dir, stdio: 'ignore' });
  const bands = parseOutput(readFileSync(join(dir, 'OUTPUT_RRTM'), 'utf8'));
  rmSync(dir, { recursive: true });
  return bands;
}

function standardColumns() {
  const cols = ['TROP', 'MLS', 'MLW', 'SAW'].map((a) => ({ name: a, layers: parseInput(readFileSync(join(process.env.RRTMG_INPUTS, `input_rrtm_${a}-clr`), 'utf8')) }));
  const bench = JSON.parse(readFileSync(join(REPO, 'data/radiationBenchmark.json'), 'utf8')).atmospheres;
  for (const a of ['AFGL_TR', 'AFGL_MS', 'AFGL_SS', 'AFGL_MW', 'AFGL_SW']) {
    cols.push({ name: a, layers: bench[a].layers.filter((l) => l.pBottom > 1).map((l) => {
      const pTop = Math.max(l.pTop, 0.5);
      return { pBot: l.pBottom, pTop, p: 0.5 * (l.pBottom + pTop), T: 0.5 * (l.tBottom + l.tTop), vmr: { h2o: l.h2o, co2: l.co2, o3: l.o3, n2o: l.n2o, co: 0, ch4: l.ch4, o2: 0.2095 }, coldry: dryColumn(l.pBottom - l.pTop, l.h2o) };
    }) });
  }
  return cols;
}

async function modelColumns(file) {
  const { topographyFromInt16 } = await import(join(REPO, 'js/geography.module.js'));
  const { decodeState, savedLevels } = await import(join(REPO, 'js/stateFile.module.js'));
  const { stateGrid } = await import(join(REPO, 'js/stateGrid.module.js'));
  const { ozoneWeights, ozoneAbove } = await import(join(REPO, 'js/physics/ozone.module.js'));
  const { GREENHOUSE_GASES } = await import(join(REPO, 'js/physics/radiation.module.js'));
  const topography = topographyFromInt16(readFileSync(join(REPO, 'data/topography_0p25.bin')).buffer);
  const s = await decodeState(new Uint8Array(readFileSync(file)));
  const mesh = stateGrid(s, topography).mesh, levels = savedLevels(s), K = levels.length - 1, C = mesh.latCell.length;
  const kappa = 287.06 / 1003.5, p0 = 101325, f = ((s.day ?? 1825) % 365) / 365;
  return TARGETS.map(([name, lat, lon]) => {
    const rad = Math.PI / 180, x = Math.cos(lat * rad) * Math.cos(lon * rad), y = Math.cos(lat * rad) * Math.sin(lon * rad), z = Math.sin(lat * rad);
    let i = 0, bd = Infinity;
    for (let c = 0; c < C; c++) { const a = mesh.latCell[c], b = mesh.lonCell[c], d = (Math.cos(a) * Math.cos(b) - x) ** 2 + (Math.cos(a) * Math.sin(b) - y) ** 2 + (Math.sin(a) - z) ** 2; if (d < bd) { bd = d; i = c; } }
    const ps = s.pi[i], w = ozoneWeights(mesh.latCell[i], f), layers = [];
    let exUpper = 0;
    for (let k = 0; k < K; k++) {
      const exLower = Math.pow(ps * levels[k + 1] / p0, kappa), ex = (exLower * levels[k + 1] - exUpper * levels[k]) / ((1 + kappa) * (levels[k + 1] - levels[k]));
      exUpper = exLower;
      const q = Math.max(1e-9, s.q[k * C + i]), r = q / (1 - q) * MASS.air / MASS.h2o, pBot = levels[k + 1] * ps;
      const coldry = dryColumn(pBot - levels[k] * ps, r), o3 = (ozoneAbove(pBot, w) - ozoneAbove(levels[k] * ps, w)) * 2.6868e19 / coldry;
      layers.push({ pBot, pTop: Math.max(levels[k] * ps, 0.5), p: 0.5 * (levels[k] + levels[k + 1]) * ps, T: s.theta[k * C + i] * ex, vmr: { h2o: r, co2: GREENHOUSE_GASES.carbonDioxide, o3, n2o: GREENHOUSE_GASES.nitrousOxide, co: 0, ch4: GREENHOUSE_GASES.methane, o2: 0.2095 }, coldry });
    }
    return { name: `${name} (${(mesh.latCell[i] / rad).toFixed(0)}, ${(mesh.lonCell[i] / rad).toFixed(0)})`, layers: layers.reverse() };
  });
}

async function reference(out) {
  const mus = (process.env.MUS ?? '1,0.5,0.25').split(',').map(Number);
  const cols = standardColumns().concat(process.env.STATE ? await modelColumns(process.env.STATE) : []);
  const result = { mus, bands: null, columns: [] };
  for (const c of cols) {
    const entry = { name: c.name, layers: c.layers.slice().reverse().map((L) => {
      const u = (m) => L.coldry * 1e4 * L.vmr[m] * MASS[m] / AVOGADRO / 1000;
      return { pTop: L.pTop, pBot: L.pBot, p: L.p, T: L.T, h2o: u('h2o'), co2: u('co2'), ch4: u('ch4'), o2: u('o2'), e: L.p * L.vmr.h2o / (1 + L.vmr.h2o) };
    }), runs: {} };
    for (const mu of mus) {
      const r = {};
      for (const [key, keep] of Object.entries(CONFIGS)) {
        const bands = runRrtmg(writeInput(c.layers, Math.acos(mu) * 180 / Math.PI, keep));
        if (!result.bands && mu === 1) result.bands = bands.map((b) => ({ range: b.range, toa: b.levels[0].down }));
        r[key] = bands.map((b) => b.levels.map((x) => x.dirdown));
      }
      entry.runs[mu] = r;
    }
    result.columns.push(entry);
  }
  writeFileSync(out, JSON.stringify(result));
  console.log(`${result.columns.length} columns, cosines ${mus} -> ${out}`);
}

// ---- the scheme ----
let SELF_BUFFER = new Float64Array(256);
export function transmission(par, layers, gas, mu, out) {
  const J = par.k.length, L = layers.length, self = par.self ?? 0, selfT = par.selfT ?? 0;
  if (SELF_BUFFER.length <= L) SELF_BUFFER = new Float64Array(L + 1);
  const ts = SELF_BUFFER;
  let w0 = 1;
  for (let j = 0; j < J; j++) w0 -= par.w[j];
  ts[0] = 0;
  for (let l = 0; l < L; l++) { const y = layers[l]; ts[l + 1] = ts[l] + (self > 0 ? self * y[gas] * y.e / P_REF * Math.pow(296 / y.T, selfT) : 0); }
  for (let l = 0; l <= L; l++) out[l] = w0 * Math.exp(-ts[l] / mu);
  for (let j = 0; j < J; j++) {
    let tau = 0;
    out[0] += par.w[j];
    for (let l = 0; l < L; l++) {
      const y = layers[l];
      tau += par.k[j] * y[gas] * Math.pow(Math.sqrt(y.p * y.p + par.pD * par.pD) / P_REF, par.n[j]) * Math.exp(par.t[j] * (y.T - 250));
      out[l + 1] += par.w[j] * Math.exp(-(tau + ts[l + 1]) / mu);
    }
  }
  return out;
}

function unpack(x, J, withSelf) {
  const e = x.slice(0, J).map(Math.exp), s = 1 + e.reduce((a, b) => a + b, 0);
  const par = { w: e.map((v) => v / s), k: x.slice(J, 2 * J).map(Math.exp), n: x.slice(2 * J, 3 * J).map((v) => 1.5 * Math.tanh(v)), t: x.slice(3 * J, 4 * J).map((v) => 0.02 * Math.tanh(v)) };
  if (withSelf) { par.self = Math.exp(x[4 * J]); par.selfT = 4 + 4 * Math.tanh(x[4 * J + 1]); }
  return par;
}

function pack(par) {
  const w0 = Math.max(1e-9, 1 - par.w.reduce((a, b) => a + b, 0)), clip = (v) => Math.max(-0.999999, Math.min(0.999999, v));
  const x = [...par.w.map((v) => Math.log(v / w0)), ...par.k.map(Math.log), ...par.n.map((v) => Math.atanh(clip(v / 1.5))), ...par.t.map((v) => Math.atanh(clip(v / 0.02)))];
  if (par.self !== undefined) x.push(Math.log(par.self), Math.atanh(clip((par.selfT - 4) / 4)));
  return x;
}

function solve(A, b) {
  const n = b.length, M = A.map((row, i) => [...row, b[i]]);
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let i = c + 1; i < n; i++) if (Math.abs(M[i][c]) > Math.abs(M[p][c])) p = i;
    [M[c], M[p]] = [M[p], M[c]];
    if (Math.abs(M[c][c]) < 1e-300) return null;
    for (let i = c + 1; i < n; i++) { const f = M[i][c] / M[c][c]; for (let j = c; j <= n; j++) M[i][j] -= f * M[c][j]; }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) { let s = M[i][n]; for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j]; x[i] = s / M[i][i]; }
  return x;
}

function levenbergMarquardt(f, x0, iterations) {
  let x = x0.slice(), r = f(x), cost = r.reduce((s, v) => s + v * v, 0), lambda = 1e-2;
  for (let it = 0; it < iterations; it++) {
    const n = x.length, Jac = [];
    for (let i = 0; i < n; i++) {
      const h = 1e-6 * Math.max(1, Math.abs(x[i])), xp = x.slice();
      xp[i] += h;
      const rp = f(xp);
      Jac.push(rp.map((v, m) => { const d = (v - r[m]) / h; return Number.isFinite(d) ? d : 0; }));
    }
    const JtJ = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => { let s = 0; for (let m = 0; m < r.length; m++) s += Jac[i][m] * Jac[j][m]; return s; }));
    const Jtr = Array.from({ length: n }, (_, i) => { let s = 0; for (let m = 0; m < r.length; m++) s += Jac[i][m] * r[m]; return -s; });
    const floor = 1e-4 * JtJ.reduce((m, row, i) => Math.max(m, row[i]), 0);
    let improved = false;
    for (let tries = 0; tries < 10; tries++) {
      const d = solve(JtJ.map((row, i) => row.map((v, j) => (i === j ? v + lambda * (v + floor) + 1e-30 : v))), Jtr);
      if (d) {
        const xn = x.map((v, i) => v + d[i]), rn = f(xn), cn = rn.reduce((s, v) => s + v * v, 0);
        if (Number.isFinite(cn) && cn < cost) {
          const rel = (cost - cn) / cost;
          x = xn; r = rn; cost = cn; lambda = Math.max(1e-9, lambda / 3); improved = true;
          if (rel < 1e-9) return { x, cost };
          break;
        }
      }
      lambda *= 4;
    }
    if (!improved) break;
  }
  return { x, cost };
}

function casesFor(REF, block, config) {
  const cases = [];
  for (const col of REF.columns) {
    const layers = col.layers, L = layers.length, ps = layers[L - 1].pBot;
    for (const mu of REF.mus) {
      const run = col.runs[mu], heat = new Float64Array(L);
      let column = 0;
      for (const b of NIR_BLOCKS) {
        const F = REF.bands[b].toa * mu, t = run.nir[b].map((v, l) => v / run.none[b][l]);
        for (let l = 0; l < L; l++) heat[l] += F * (t[l] - t[l + 1]);
        column += F * (1 - t[L]);
      }
      cases.push({ mu, layers, ref: run[config][block].map((v, l) => Math.min(1, v / run.none[block][l])), buf: new Float64Array(L + 1), flux: REF.bands[block].toa * mu,
        scale: Float64Array.from(heat, (h, l) => Math.abs(h) + 0.003 * column * (layers[l].pBot - layers[l].pTop) / ps), columnScale: 0.002 * column });
    }
  }
  return cases;
}

function soloResiduals(par, cases, gas, r = [], weight = 1) {
  for (const c of cases) {
    const T = transmission(par, c.layers, gas, c.mu, c.buf), L = c.layers.length;
    for (let l = 0; l < L; l++) r.push(weight * c.flux * ((T[l] - T[l + 1]) - (c.ref[l] - c.ref[l + 1])) / c.scale[l]);
    r.push(weight * c.flux * (T[L] - c.ref[L]) / c.columnScale);
  }
  return r;
}

function fitSolo(REF, block, gas, J) {
  const cases = casesFor(REF, block, gas), withSelf = gas === 'h2o';
  const pD = dopplerPressure(gas, CENTERS[BAND_OF_BLOCK[block]][gas]);
  const column = Math.max(...cases.map((c) => c.layers.reduce((s, y) => s + y[gas], 0)));
  let best = null;
  for (let start = 0; start < 6; start++) {
    const spread = 1 + 0.5 * (start % 3), shift = start < 3 ? 0 : 1.5;
    const x0 = [...new Array(J).fill(-1.5), ...Array.from({ length: J }, (_, j) => Math.log(10 ** (spread * (j - 1) + shift) / column)),
      ...Array.from({ length: J }, (_, j) => (start < 3 ? 0.3 : 0.8 - 0.6 * j)), ...new Array(J).fill(0)];
    if (withSelf) x0.push(Math.log(1e-3 / column * 10 ** (start % 3)), 0);
    const fit = levenbergMarquardt((x) => soloResiduals({ pD, ...unpack(x, J, withSelf) }, cases, gas), x0, 200);
    if (!best || fit.cost < best.cost) best = fit;
  }
  return { pD, ...unpack(best.x, J, withSelf) };
}

function fitJoint(REF, block, start) {
  const gases = Object.keys(start), sizes = gases.map((g) => pack(start[g]).length), weight = Number(process.env.SOLO_WEIGHT ?? 0.3);
  const combined = casesFor(REF, block, 'nir'), solo = Object.fromEntries(gases.map((g) => [g, casesFor(REF, block, g)]));
  const split = (x) => { let o = 0; return gases.map((g, i) => { const p = { pD: start[g].pD, ...unpack(x.slice(o, o + sizes[i]), start[g].w.length, start[g].self !== undefined) }; o += sizes[i]; return p; }); };
  const f = (x) => {
    const pars = split(x), r = [];
    for (const c of combined) {
      const L = c.layers.length, T = new Float64Array(L + 1).fill(1);
      gases.forEach((g, i) => { const t = transmission(pars[i], c.layers, g, c.mu, c.buf); for (let l = 0; l <= L; l++) T[l] *= t[l]; });
      for (let l = 0; l < L; l++) r.push(c.flux * ((T[l] - T[l + 1]) - (c.ref[l] - c.ref[l + 1])) / c.scale[l]);
      r.push(c.flux * (T[L] - c.ref[L]) / c.columnScale);
    }
    gases.forEach((g, i) => soloResiduals(pars[i], solo[g], g, r, weight));
    return r;
  };
  const fit = levenbergMarquardt(f, gases.flatMap((g) => pack(start[g])), 300);
  return Object.fromEntries(split(fit.x).map((p, i) => [gases[i], p]));
}

function visibleStrength(REF) {
  const { VISIBLE_VAPOR, vaporScaling } = SHORTWAVE;
  let num = 0, den = 0;
  for (const col of REF.columns) {
    let w = 0;
    for (const y of col.layers) w += y.h2o * 0.1 * vaporScaling(y.p, y.T);
    for (const mu of REF.mus) {
      const run = col.runs[mu], L = col.layers.length, m = 35 / Math.sqrt(1224 * mu * mu + 1);
      let ref = 0;
      for (const b of VISIBLE_BLOCKS) ref += REF.bands[b].toa * mu * (1 - Math.min(1, run.h2o[b][L] / run.none[b][L]));
      const slope = REF.bands[0].toa * mu * VISIBLE_VAPOR.share * VISIBLE_VAPOR.coefficient * w * m;
      num += slope * ref; den += slope * slope;
    }
  }
  return num / den;
}

function round(x) { return Number(x.toPrecision(6)); }
function tableOf(REF, fits, strength) {
  return {
    pRef: P_REF, visibleVaporStrength: round(strength),
    bands: fits.map(({ block, gases }) => ({ band: BAND_OF_BLOCK[block], range: REF.bands[block].range, share: round(REF.bands[block].toa / REF.bands[0].toa),
      gases: Object.fromEntries(Object.entries(gases).map(([g, p]) => [g, { pD: round(p.pD), w: p.w.map(round), k: p.k.map(round), n: p.n.map(round), t: p.t.map(round), ...(p.self !== undefined ? { self: round(p.self), selfT: round(p.selfT) } : {}) }])) })),
  };
}

// The table's heating (K/day) and absorption against the reference's four-gas run.
function report(REF, table) {
  const P = [1, 3, 7, 14, 27, 54, 100, 200, 300, 500, 700, 850, 950];
  let sq = 0, n = 0, sqS = 0, nS = 0;
  const misses = [];
  for (const col of REF.columns) {
    for (const mu of REF.mus) {
      const run = col.runs[mu], L = col.layers.length, ref = new Float64Array(L), mod = new Float64Array(L), buf = new Float64Array(L + 1);
      let refColumn = 0, modColumn = 0;
      for (const entry of table.bands) {
        const b = BAND_OF_BLOCK.indexOf(entry.band), F = REF.bands[b].toa * mu, tr = run.nir[b].map((v, l) => Math.min(1, v / run.none[b][l])), tm = new Float64Array(L + 1).fill(1);
        for (const [gas, par] of Object.entries(entry.gases)) { transmission(par, col.layers, gas, mu, buf); for (let l = 0; l <= L; l++) tm[l] *= buf[l]; }
        for (let l = 0; l < L; l++) { ref[l] += F * (tr[l] - tr[l + 1]); mod[l] += F * (tm[l] - tm[l + 1]); }
        refColumn += F * (1 - tr[L]); modColumn += F * (1 - tm[L]);
      }
      const rate = (f, l) => f * G / (CP * (col.layers[l].pBot - col.layers[l].pTop)) * 86400;
      col.layers.forEach((y, l) => { const d = rate(mod[l] - ref[l], l); if (y.p > 20000) { sq += d * d; n++; } else if (y.p > 500) { sqS += d * d; nS++; } });
      misses.push(modColumn - refColumn);
      if (mu === 1 && (col.name.length <= 4 || col.name.startsWith('warm') || col.name.startsWith('Antarctic'))) {
        const rows = P.map((p) => col.layers.reduce((b, y, l) => (Math.abs(Math.log(y.p / (p * 100))) < Math.abs(Math.log(col.layers[b].p / (p * 100))) ? l : b), 0));
        console.log(`${col.name}, overhead sun: absorption ${modColumn.toFixed(1)} against RRTMG ${refColumn.toFixed(1)} W/m2; heating K/day at hPa`);
        console.log('  ' + rows.map((l) => `${(col.layers[l].p / 100).toPrecision(3)}: ${rate(mod[l], l).toFixed(3)}/${rate(ref[l], l).toFixed(3)}`).join('  '));
      }
    }
  }
  const mean = misses.reduce((a, b) => a + b, 0) / misses.length, rms = Math.sqrt(misses.reduce((a, b) => a + b * b, 0) / misses.length);
  console.log(`over ${misses.length} columns and angles: heating rms ${Math.sqrt(sqS / nS).toFixed(3)} K/day at 5-200 hPa, ${Math.sqrt(sq / n).toFixed(3)} below; absorption miss mean ${mean.toFixed(2)}, rms ${rms.toFixed(2)}, ${Math.min(...misses).toFixed(2)} to ${Math.max(...misses).toFixed(2)} W/m2`);
}

const SHORTWAVE = await import(join(REPO, 'js/physics/shortwaveGases.module.js'));
const [mode, refFile, fitFile] = process.argv.slice(2);
if (mode === 'reference') await reference(refFile);
else if (mode === 'fit') {
  const REF = JSON.parse(readFileSync(refFile, 'utf8')), J = Number(process.env.TERMS ?? 4), fits = [];
  for (const block of NIR_BLOCKS) {
    const band = BAND_OF_BLOCK[block], solo = {};
    for (const gas of Object.keys(CENTERS[band])) solo[gas] = fitSolo(REF, block, gas, J);
    const gases = Object.keys(solo).length > 1 ? fitJoint(REF, block, solo) : solo;
    fits.push({ block, gases });
    console.log(`band ${band}: ${Object.entries(gases).map(([g, p]) => `${g} w ${p.w.map((v) => v.toFixed(3))} k ${p.k.map((v) => v.toExponential(2))} n ${p.n.map((v) => v.toFixed(2))}`).join('; ')}`);
  }
  const table = tableOf(REF, fits, visibleStrength(REF));
  console.log(`visible water vapor strength ${table.visibleVaporStrength}`);
  if (fitFile) writeFileSync(fitFile, JSON.stringify(table));
  report(REF, table);
  if (process.env.WRITE) {
    writeFileSync(join(REPO, 'js/physics/nearInfraredTable.module.js'), `// Generated by scripts/nearInfraredFit.mjs: for each of RRTMG_SW's bands 16-24 and 29, its share of the
// solar beam and, for each gas, the Doppler pressure p_D (Pa) and the k-distribution's weights w,
// coefficients k (m2/kg at p_ref), pressure exponents n and temperature coefficients t (1/K), and water
// vapor's self continuum (m2/kg per p_ref of vapor pressure at 296 K) and its temperature exponent.
export const NEAR_INFRARED_TABLE = ${JSON.stringify(table)};
`);
    console.log('wrote js/physics/nearInfraredTable.module.js');
  }
} else if (mode === 'report') {
  const REF = JSON.parse(readFileSync(refFile, 'utf8'));
  report(REF, fitFile ? JSON.parse(readFileSync(fitFile, 'utf8')) : (await import(join(REPO, 'js/physics/nearInfraredTable.module.js'))).NEAR_INFRARED_TABLE);
} else {
  console.error('usage: node scripts/nearInfraredFit.mjs reference|fit|report <ref.json> [<fit.json>]');
  process.exit(2);
}
