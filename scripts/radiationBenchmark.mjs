// The gas radiation on fixed profiles against its references, each term by
// itself (data/radiationBenchmark.json, sources in the file):
//   node scripts/radiationBenchmark.mjs
// One column of the CPU radiation (js/physics/radiation.module.js) on the
// bl34 levels for each standard atmosphere: clear sky, a black surface in the
// longwave, the reference's own temperature, vapour, ozone and well-mixed
// gases (the layer means of scripts/standardAtmospheres.mjs), no aerosol.
// Two schemes side by side: BEFORE (RADIATION_BEFORE, JSON; by default the
// three-band grey longwave and the Lacis-Hansen vapour with a fixed ozone
// share) and AFTER (RADIATION, JSON; by default the radiation's defaults).
// Longwave: OLR, surface downward flux, net flux at 200 hPa and the
// cooling-rate profile against RRTMG (tropical, midlatitude summer and
// winter, subarctic winter) and against the ICRCCM line-by-line fluxes of
// Feigelson et al. (1991) for the five AFGL atmospheres (CO2 300 ppmv, no
// methane or nitrous oxide); doubled CO2 and vapour x1.2 on the midlatitude
// summer profile against LBLRTM (Iacono et al. 2008), and methane and nitrous
// oxide from none to their 1860 amounts on the same profile, with the
// stratosphere-adjusted forcing (fixed dynamical heating); the OLR's slope
// with surface temperature at fixed relative humidity. Shortwave: the
// atmosphere's absorption, the surface's downward flux and the heating
// profile against RRTMG (overhead sun and, midlatitude summer, 65 degrees,
// surface albedo 0.2, Rayleigh scattering) and each gas against the
// line-by-line values of Chou & Suarez (1999) on the midlatitude summer
// profile at 60 degrees without scattering.
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore, sigmaInterfaces, GRAVITY, CP_DRY } from '../js/dynamics/sigmaCore.module.js';
import { createRadiation } from '../js/physics/radiation.module.js';
import { OZONE_CM_ATM } from '../js/physics/shortwaveGases.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { BENCHMARK, MOLAR, modelColumn, meanMixingRatio, referenceAt, referenceHeating, layerHeating, interfaceAt } from './standardAtmospheres.mjs';

const BEFORE = { longwaveScheme: 'gray', solarGases: 'lacisHansen', ...JSON.parse(process.env.RADIATION_BEFORE ?? '{}') };
const AFTER = JSON.parse(process.env.RADIATION ?? '{}');
const levels = sigmaInterfaces('bl34');
const mesh = buildMesh(new Grid(2));
const core = createSigmaCore(mesh, { levels });
const { K, C } = core.diagnostics;
const CLEAN = { landAerosol: 0, seaAerosol: 0, clearSkyPass: true };

function gasesOf(column) {
  const ozone = Array.from(column.o3, (x, k) => x * (levels[k + 1] - levels[k]) * column.ps / GRAVITY / OZONE_CM_ATM);
  return { carbonDioxide: meanMixingRatio(column, 'co2'), methane: meanMixingRatio(column, 'ch4'), nitrousOxide: meanMixingRatio(column, 'n2o'), ozoneProfile: Float64Array.from(ozone) };
}

// One column: longwave by night (beam 0) and shortwave at the given beam.
export function runColumn(options, column, { beam = 0, albedo = 0.2, solarConstant = 1360.85 } = {}) {
  const radiation = createRadiation(mesh, core, { ...CLEAN, solarConstant, ...gasesOf(column), ...options });
  const pi = new Float64Array(C).fill(column.ps), theta = new Float64Array(K * C).fill(280), q = new Float64Array(K * C), qc = new Float64Array(K * C);
  core.diagnose(pi, theta);
  const { exnerLayer } = core.diagnostics;
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) { theta[k * C + i] = column.T[k] / exnerLayer[k * C + i]; q[k * C + i] = column.q[k]; }
  core.diagnose(pi, theta, q, qc);
  radiation.column(0, column.ps, theta, column.Ts, 5, undefined, beam, q[(K - 1) * C], q, qc, albedo, albedo, 0, 0, 0, 0, 0, 0);
  const b = { ...radiation.budget };
  const lw = Array.from({ length: K }, (_, k) => radiation.longwave[k * C]);
  const sw = Array.from({ length: K }, (_, k) => radiation.layerFlux[k] - lw[k]);
  const net = [b.outgoingLongwave];
  for (let k = 0; k < K; k++) net.push(net[k] + lw[k]);
  return { budget: b, olr: b.outgoingLongwave, dlr: b.downwardLongwave, net, lw, sw };
}

const f = (x, d = 1) => (Number.isFinite(x) ? x.toFixed(d) : 'n/a');
const pct = (x, ref) => `${x >= 0 ? '+' : ''}${f(x, 1)} (${x >= 0 ? '+' : ''}${f(100 * x / ref, 1)} %)`;
const rms = (a, b, keep) => { let s = 0, n = 0; a.forEach((x, k) => { if (keep(k)) { s += (x - b[k]) ** 2; n++; } }); return Math.sqrt(s / n); };
const columns = Object.fromEntries(Object.keys(BENCHMARK.atmospheres).map((a) => [a, modelColumn(BENCHMARK.atmospheres[a], levels)]));
const midPressure = (column, k) => 0.5 * (levels[k] + levels[k + 1]) * column.ps;
const withGas = (column, gas, vmr) => ({ ...column, [gas]: column[gas].map((_, k) => vmr * MOLAR[gas] / MOLAR.air * (1 - column.q[k])) });

function longwaveTable() {
  console.log('LONGWAVE, clear sky, black surface; misses against the reference, W/m2 (per cent); cooling rms K/day for p > 200 hPa and 3-200 hPa');
  console.log('atmosphere  ref         | before: OLR            DLR            net 200 hPa    cooling        | after: OLR             DLR            net 200 hPa    cooling');
  for (const a of ['TROP', 'MLS', 'MLW', 'SAW']) {
    const column = columns[a], ref = BENCHMARK.rrtmgLongwave[a].levels, top = ref[ref.length - 1], sfc = ref[0];
    const refHeat = referenceHeating(column, ref), ref200 = referenceAt(ref, 20000);
    const cells = [BEFORE, AFTER].map((o) => {
      const r = runColumn(o, column), heat = layerHeating(column, r.net);
      return `${pct(r.olr - top.up, top.up).padEnd(15)}${pct(r.dlr - sfc.down, sfc.down).padEnd(15)}${pct(interfaceAt(column, r.net, 20000) - ref200, ref200).padEnd(15)}${f(rms(heat, refHeat, (k) => midPressure(column, k) > 20000), 2)}/${f(rms(heat, refHeat, (k) => midPressure(column, k) <= 20000 && midPressure(column, k) > 300), 2)}`.padEnd(64);
    });
    console.log(`${a.padEnd(11)} RRTMG ${f(top.up)}/${f(sfc.down)} | ${cells[0]}| ${cells[1]}`);
  }
  for (const a of ['AFGL_TR', 'AFGL_MS', 'AFGL_SS', 'AFGL_MW', 'AFGL_SW']) {
    const column = columns[a], [, down, olr] = BENCHMARK.iae[a];
    const cells = [BEFORE, AFTER].map((o) => { const r = runColumn(o, column); return `${pct(r.olr - olr, olr).padEnd(15)}${pct(r.dlr - down, down).padEnd(15)}`.padEnd(64); });
    console.log(`${a.padEnd(11)} IAE ${f(olr)}/${f(down)}   | ${cells[0]}| ${cells[1]}`);
  }
  const mls = columns.MLS, m = BENCHMARK.mlawerMls;
  const r = runColumn(AFTER, mls);
  console.log(`MLS against LBLRTM (Mlawer et al. 1997, CKD_2.0): OLR ${f(r.olr)} (${m.toaUp}), net at 179 hPa ${f(interfaceAt(mls, r.net, m.tropopause.p))} (${m.tropopause.net}), DLR ${f(r.dlr)} (${m.surface.down})`);
  console.log('cooling of the layers above 30 hPa, K/day (before / after / RRTMG):');
  for (const a of ['TROP', 'MLS', 'MLW', 'SAW']) {
    const column = columns[a], refHeat = referenceHeating(column, BENCHMARK.rrtmgLongwave[a].levels);
    const heats = [BEFORE, AFTER].map((o) => layerHeating(column, runColumn(o, column).net));
    console.log(`  ${a.padEnd(4)} ` + Array.from({ length: K }, (_, k) => k).filter((k) => midPressure(column, k) < 3000).map((k) => `${f(midPressure(column, k) / 100, 1)} hPa ${f(heats[0][k], 2)}/${f(heats[1][k], 2)}/${f(refHeat[k], 2)}`).join('; '));
  }
  for (const a of ['MLS', 'TROP']) {
    const column = columns[a], ref = BENCHMARK.rrtmgLongwave[a].levels, refHeat = referenceHeating(column, ref);
    const heats = [BEFORE, AFTER].map((o) => layerHeating(column, runColumn(o, column).net));
    console.log(`${a} cooling profile, K/day (before / after / RRTMG):`);
    console.log('  ' + Array.from({ length: K }, (_, k) => `${f(midPressure(column, k) / 100, 0)} hPa ${f(heats[0][k], 2)}/${f(heats[1][k], 2)}/${f(refHeat[k], 2)}`).join('; '));
  }
}

// Stratospheric temperatures adjusted (fixed dynamical heating) above the
// tropopause until each layer's longwave heating is its base value again.
function adjusted(options, base, perturbed, tropopause) {
  const target = runColumn(options, base).lw;
  const column = { ...perturbed, T: Float64Array.from(perturbed.T) };
  const strat = Array.from({ length: K }, (_, k) => midPressure(column, k) < tropopause);
  for (let it = 0; it < 60; it++) {
    const r = runColumn(options, column);
    const dT = 0.5;
    const warmer = { ...column, T: column.T.map((t, k) => (strat[k] ? t + dT : t)) };
    const rw = runColumn(options, warmer);
    for (let k = 0; k < K; k++) if (strat[k]) {
      const slope = (rw.lw[k] - r.lw[k]) / dT;
      if (slope < 0) column.T[k] -= Math.max(-5, Math.min(5, (r.lw[k] - target[k]) / slope)) * 0.7;
    }
  }
  return column;
}

function sensitivityTable() {
  const mls = columns.MLS;
  const rtmip = (co2) => withGas(withGas(withGas(mls, 'co2', co2), 'ch4', 806e-9), 'n2o', 275e-9);
  const ref = BENCHMARK.iacono.co2Doubling.longwave, vap = BENCHMARK.iacono.vapour120.longwave;
  console.log('\nSENSITIVITIES, midlatitude summer, LBLRTM (Iacono et al. 2008)');
  for (const [name, o] of [['before', BEFORE], ['after', AFTER]]) {
    const one = rtmip(287e-6), two = rtmip(574e-6), wet = { ...two, q: two.q.map((x) => x * 1.2) };
    const r1 = runColumn(o, one), r2 = runColumn(o, two), r3 = runColumn(o, wet);
    const at = (r, p) => interfaceAt(mls, r.net, p);
    const forcing = (a, b) => `TOA ${f(a.olr - b.olr, 2)}  200 hPa ${f(at(a, 20000) - at(b, 20000), 2)}  surface ${f(b.dlr - a.dlr, 2)}`;
    console.log(`  ${name.padEnd(6)} CO2 287 -> 574 ppmv: ${forcing(r1, r2)}   (LBLRTM ${ref.toa} / ${ref.p20000} / ${ref.surface})`);
    console.log(`  ${name.padEnd(6)} vapour x1.2 at 574:  ${forcing(r2, r3)}   (LBLRTM ${vap.toa} / ${vap.p20000} / ${vap.surface})`);
    const minor = BENCHMARK.iacono.minorGases.longwave, none = runColumn(o, withGas(withGas(one, 'ch4', 0), 'n2o', 0));
    console.log(`  ${name.padEnd(6)} CH4 0 -> 806 ppbv and N2O 0 -> 275 ppbv at 287 ppmv CO2: ${forcing(none, r1)}   (LBLRTM ${minor.toa} / ${minor.p20000} / ${minor.surface})`);
    if (o === AFTER) {
      const tropopause = 17900, adjustedColumn = adjusted(o, one, two, tropopause), ra = runColumn(o, adjustedColumn);
      const dT = adjustedColumn.T.map((t, k) => t - two.T[k]);
      console.log(`  ${name.padEnd(6)} stratosphere-adjusted (fixed dynamical heating above 179 hPa): TOA ${f(r1.olr - ra.olr, 2)}  179 hPa ${f(interfaceAt(mls, r1.net, tropopause) - interfaceAt(mls, ra.net, tropopause), 2)} W/m2; stratospheric cooling ${f(Math.min(...dT), 1)} K at most`);
    }
    for (const a of ['MLS', 'TROP']) {
      const c = columns[a], slope = (delta) => {
        const warm = { ...c, Ts: c.Ts + delta, T: c.T.map((t, k) => (midPressure(c, k) > 15000 ? t + delta : t)) };
        warm.q = c.q.map((x, k) => (midPressure(c, k) > 15000 ? x * saturationHumidity(warm.T[k], midPressure(c, k)) / saturationHumidity(c.T[k], midPressure(c, k)) : x));
        return runColumn(o, warm).olr;
      };
      console.log(`  ${name.padEnd(6)} ${a} OLR slope with surface and tropospheric temperature at fixed relative humidity: ${f((slope(1) - slope(-1)) / 2, 2)} W/m2/K (about 2; Koll & Cronin 2018, from memory)`);
    }
  }
}

function shortwaveTable() {
  console.log('\nSHORTWAVE, clear sky, no aerosol; RRTMG with Rayleigh scattering, surface albedo 0.2 (65 degrees: its spectral albedo, 0.213 in effect)');
  console.log('atmosphere  mu    RRTMG atm / sfc down | before: atmosphere      surface down    heating rms | after: atmosphere       surface down    heating rms');
  for (const [a, key, mu, albedo] of [['TROP', 'TROP', 1, 0.2], ['MLS', 'MLS', 1, 0.2], ['MLW', 'MLW', 1, 0.2], ['SAW', 'SAW', 1, 0.2], ['MLS', 'MLS65', Math.cos(65 * Math.PI / 180), 0.2127]]) {
    const column = columns[a], ref = BENCHMARK.rrtmgShortwave[key].levels, top = ref[ref.length - 1], sfc = ref[0];
    const atm = top.net - sfc.net, refHeat = referenceHeating(column, ref).map((x) => -x);
    const cells = [BEFORE, AFTER].map((o) => {
      const r = runColumn(o, column, { beam: top.down, albedo, solarConstant: top.down / mu });
      const heat = layerHeating(column, (() => { const n = [0]; for (let k = 0; k < K; k++) n.push(n[k] + r.sw[k]); return n; })());
      return `${pct(r.budget.atmosphereSolar - atm, atm).padEnd(16)}${pct(r.budget.surfaceShortwave - sfc.down, sfc.down).padEnd(16)}${f(rms(heat, refHeat, (k) => midPressure(column, k) > 20000), 2)}/${f(rms(heat, refHeat, (k) => midPressure(column, k) <= 20000 && midPressure(column, k) > 100), 2)}`.padEnd(50);
    });
    console.log(`${a.padEnd(11)} ${f(mu, 3)} ${f(atm)} / ${f(sfc.down)}    | ${cells[0]}| ${cells[1]}`);
  }
  for (const a of ['MLS']) {
    const bands = BENCHMARK.rrtmgShortwave.MLS.bands;
    const group = (lo, hi) => bands.filter((b) => { const [x, y] = b.range.match(/[\d.]+/g).map(Number); return x >= lo && y <= hi; }).reduce((s, b) => s + b.toaNet - b.surfaceNet, 0);
    const column = columns[a];
    for (const o of [BEFORE, AFTER]) {
      const r = runColumn(o, column, { beam: 1360.85, albedo: 0.2 });
      console.log(`  ${o === BEFORE ? 'before' : 'after '} MLS overhead sun by term: ozone ${f(r.budget.ozoneSolar)} vapour ${f(r.budget.vaporSolar)} O2 ${f(r.budget.oxygenSolar)} CO2 ${f(r.budget.carbonDioxideSolar)} reflected light on its way up ${f(r.budget.upwardGasSolar)} W/m2`);
    }
    console.log(`  RRTMG MLS overhead sun by band group: ultraviolet and visible above 12850 cm-1 (O3, O2 A and B, weak H2O) ${f(group(12850, 50000))}, near infrared (H2O, CO2, CH4, O2 1.27 um) ${f(group(820, 12850))} W/m2`);
  }
  const mls = columns.MLS, cs = BENCHMARK.chouShortwave, beam = 0.5 * 1365;
  console.log(`  Chou & Suarez (1999) line-by-line, MLS, 60 degrees, no scattering (insolation taken as ${beam}, solar constant 1365):`);
  for (const o of [BEFORE, AFTER]) {
    const r = runColumn({ ...o, rayleighDepth: 0, skylight: 0, upwardAbsorption: true }, withGas(mls, 'co2', 350e-6), { beam, albedo: 0.2, solarConstant: 1365 });
    const r0 = runColumn({ ...o, rayleighDepth: 0, skylight: 0 }, withGas(mls, 'co2', 350e-6), { beam, albedo: 0, solarConstant: 1365 });
    console.log(`    ${o === BEFORE ? 'before' : 'after '} atmosphere ${f(r.budget.atmosphereSolar)} (${cs.noScattering.total[2]}), of the net at the top ${f(100 * r.budget.atmosphereSolar / (beam - r.budget.reflectedSolar))} % (${f(100 * cs.noScattering.total[2] / cs.noScattering.total[0])} %); albedo 0: O2 ${f(r0.budget.oxygenSolar, 2)} (${-cs.surfaceReduction.o2}) CO2 ${f(r0.budget.carbonDioxideSolar, 2)} (${-cs.surfaceReduction.co2}) ozone ${f(r0.budget.ozoneSolar, 2)} (ultraviolet ${cs.noScattering.uv[2]} of the band 1-7 absorption with albedo 0.2) vapour ${f(r0.budget.vaporSolar, 2)}`);
  }
}

// The vapour strength that minimises the squared misses of the atmosphere's
// absorption against RRTMG over the five shortwave cases.
function fitVaporStrength() {
  const cases = [['TROP', 'TROP', 1, 0.2], ['MLS', 'MLS', 1, 0.2], ['MLW', 'MLW', 1, 0.2], ['SAW', 'SAW', 1, 0.2], ['MLS', 'MLS65', Math.cos(65 * Math.PI / 180), 0.2127]];
  const misses = (strength) => cases.map(([a, key, mu, albedo]) => {
    const ref = BENCHMARK.rrtmgShortwave[key].levels, top = ref[ref.length - 1], sfc = ref[0];
    return runColumn({ ...AFTER, vaporStrength: strength }, columns[a], { beam: top.down, albedo, solarConstant: top.down / mu }).budget.atmosphereSolar - (top.net - sfc.net);
  });
  let lo = 1, hi = 2;
  for (let it = 0; it < 30; it++) { const mid = 0.5 * (lo + hi); if (misses(mid).reduce((a, b) => a + b, 0) < 0) lo = mid; else hi = mid; }
  console.log(`vapour strength ${f(0.5 * (lo + hi), 3)}: misses ${misses(0.5 * (lo + hi)).map((x) => f(x, 1)).join(', ')} W/m2 (TROP, MLS, MLW, SAW overhead, MLS 65 degrees); at 1: ${misses(1).map((x) => f(x, 1)).join(', ')}`);
}

if (import.meta.url === `file://${process.argv[1]}` && process.env.VAPOR_FIT) fitVaporStrength();
else if (import.meta.url === `file://${process.argv[1]}`) {
  console.log(`before ${JSON.stringify(BEFORE)}; after ${JSON.stringify(AFTER)} (the radiation's defaults otherwise); bl34, ${K} layers`);
  longwaveTable();
  sensitivityTable();
  shortwaveTable();
}
