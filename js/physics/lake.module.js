/*
 * FLake, the lake model of Mironov (2008, COSMO Technical Report 11) and
 * Mironov et al. (2010, Boreal Environment Research 15, 218–230), as the
 * IFS runs it (Balsamo et al. 2012, Tellus A 64; Choulga et al. 2019,
 * HESS 23, 4051): a mixed layer of temperature T_ML and depth h over a
 * thermocline down to the bottom at depth D whose temperature falls to T_b
 * along the self-similar shape Φ_T with shape factor C_T, so the column's
 * mean is T̄ = T_ML − C_T (1 − h/D)(T_ML − T_b); ice above it with a
 * linear-to-quadratic temperature profile (Mironov and Ritter 2004), no
 * snow layer and no bottom sediments. The column follows the reference
 * code's flake_driver step for step (FLake, the release at
 * http://www.flake.igb-berlin.de/source/src_flake_sfcflx.tar.gz), with its
 * constants, flake_parameters, flake_configure, flake_albedo_ref and
 * flake_paramoptic_ref. Radiation enters as the shortwave the surface absorbs, taken
 * down through the water with one band of extinction 3 m⁻¹; the ice is
 * opaque. The lake's depth lies between 2 and 50 m, the IFS's bounds: a
 * deeper lake has a false bottom at 50 m.
 *
 * FLake is distributed under the MIT license
 * (http://www.flake.igb-berlin.de/site/download):
 *
 * Copyright (c) the FLake developers (D. Mironov et al., Deutscher
 * Wetterdienst)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
export const FLAKE = Object.freeze({
  cbl1: 0.17, cbl2: 1, zmN: 0.5, zmS: 10, zmI: 20, relaxH: 0.03, relaxC: 0.003,
  shapeMin: 0.5, shapeMax: 0.8, phi1: 40 / 3, phi2: 20 / 3, ctt1: 11 / 18, ctt2: 7 / 45,
  iceShape: 0.5, iceSlope: 1, iceAst: 2, iceMR: 1 / 12, iceMax: 3,
  iceMin: 1e-9, mixedMin: 0.01, mixedMax: 1000, ustarMin: 1e-6, small: 1e-10,
  g: 9.81, maximumDensityT: 277.13, freezing: 273.15, expansion: 1.6509e-5, rhoWater: 1000, rhoIce: 910,
  fusion: 3.3e5, cWater: 4200, cIce: 2100, kWater: 0.546, kIce: 2.29,
  extinction: 3, minimumDepth: 2, maximumDepth: 50,
});

const buoyancy = (T) => FLAKE.g * FLAKE.expansion * (T - FLAKE.maximumDensityT);

/*
 * The albedo of the lake's ice at surface temperature T (the reference
 * interface's): white ice well below freezing, blue ice at melt.
 */
export const lakeIceAlbedo = (T) => { const blue = Math.exp(-95.6 * (FLAKE.freezing - T) / FLAKE.freezing); return 0.6 * (1 - blue) + 0.1 * blue; };

export const lakeDepth = (depth) => Math.min(FLAKE.maximumDepth, Math.max(FLAKE.minimumDepth, depth));

/*
 * A new lake column at temperature T (K): mixed to the bottom, no ice.
 */
export function lakeColumn(T, depth) {
  const t = Math.max(FLAKE.freezing, T);
  return { mixed: t, mean: t, bottom: t, shape: FLAKE.shapeMin, mixedDepth: depth, ice: 0, iceT: FLAKE.freezing };
}

/*
 * One step of dt seconds of column `c` (changed in place) of depth D under
 * `absorbed`, the shortwave the surface takes in (W/m², into the ice where
 * there is ice, which passes none on), `flux`, the rest of the surface's
 * heat gain (W/m², downward), `ustar`, the water's friction velocity
 * (m/s), and the Coriolis parameter f. Returns the surface temperature.
 * With `closed` false the ice is FLake's own, which does not conserve heat.
 */
export function lakeStep(c, D, absorbed, flux, ustar, f, dt, closed = true) {
  const k = FLAKE, rc = k.rhoWater * k.cWater, target = lakeEnthalpy(c, D) + (absorbed + flux) * dt;
  const Tml = c.mixed, Tmn = c.mean, Tb = c.bottom, CT = c.shape, h = c.mixedDepth, hi = c.ice, Ti = c.iceT;
  const iced = hi >= k.iceMin;
  const Iw = iced ? 0 : absorbed, Iice = iced ? absorbed : 0;
  const decay = (z) => Math.exp(-k.extinction * z);
  const Ih = h >= k.mixedMin ? Iw * decay(h) : Iw, Ibot = Iw * decay(D);
  const I0h = h >= k.mixedMin ? Iw * (1 - decay(h)) / k.extinction / h : Ih;
  const IhD = h <= D - k.mixedMin ? Iw * (decay(h) - decay(D)) / k.extinction / (D - h) : Ih;
  let Qw = iced ? 0 : flux;
  const Qice = iced ? flux : 0;
  if (iced) Qw = h <= k.mixedMin ? -k.kWater * (Tb - Tml) / D * Math.max(k.phi1 * CT - k.phi2, 1) : 0;
  const Qstar = Qw + Iw + Ih - 2 * I0h;
  let hiN = hi, TiN = Ti, hN = h, CTN = CT, TbN = Tb, TmlN = Tml, create = false;
  if (!iced) {
    create = Tml <= k.freezing + k.small && Qw < 0;
    if (create) {
      hiN = hi - Qw / k.rhoIce / k.fusion * dt;
      TiN = k.freezing + hiN * Qw / k.kIce / k.iceSlope;
    }
  } else {
    let melting = false;
    if (Ti >= k.freezing - k.small) {
      const gain = Qice + Iice - Iw - Qw;
      if (gain >= 0) { melting = true; hiN = hi - gain / k.fusion / k.rhoIce * dt; TiN = k.freezing; }
    }
    if (!melting) {
      const ratio = hi / k.iceMax, CI = k.iceShape - k.iceMR * (1 + k.iceAst) * ratio, slope = k.iceSlope - ratio;
      let threshold = Math.max(1, 2 * CI * k.cIce * (k.freezing - Ti) / k.fusion);
      threshold = Math.sqrt(slope / CI * k.kIce / k.rhoIce / k.cIce * threshold * dt);
      threshold = Math.min(0.9 * k.iceMax, Math.max(threshold, k.iceMin));
      if (hi < threshold) {
        const top = Qice + Iice - Iw;
        hiN = hi - (top - Qw) / k.fusion / k.rhoIce * dt;
        TiN = k.freezing + hiN * top / k.kIce / slope;
      } else {
        const conducted = k.kIce * (k.freezing - Ti) / hi * slope;
        hiN = hi + (Qw + conducted) / k.fusion / k.rhoIce * dt;
        const star = (1 - CI) * k.cIce * (k.freezing - Ti) / k.fusion;
        const rate = (conducted * (1 - star) - star * Qw + (Qice + Iice - Iw)) / k.rhoIce / k.cIce / (CI * hi);
        TiN = Ti + rate * dt;
      }
    }
  }
  hiN = Math.min(hiN, k.iceMax);
  TiN = Math.max(73.15, Math.min(TiN, k.freezing));
  if (hiN < k.iceMin) { hiN = 0; TiN = k.freezing; create = false; }
  if (create) Qw = 0;
  let TmnN = Math.max(k.freezing, Tmn + (Qw + Iw) / rc / D * dt);
  if (hiN >= k.iceMin) {
    TmnN = Math.min(TmnN, k.maximumDensityT);
    TmlN = k.freezing;
    if (create) {
      if (h >= D - k.mixedMin) { hN = 0; CTN = k.shapeMin; } else { hN = h; CTN = CT; }
      TbN = TmlN - (TmlN - TmnN) / CTN / (1 - hN / D);
    } else if (Tb < k.maximumDensityT) {
      TbN = TmlN - (TmlN - TmnN) / CTN / (1 - hN / D);
    } else {
      TbN = k.maximumDensityT;
      if (h >= k.small) hN = Math.max(0, D * (1 - (TmlN - TmnN) / (TmlN - TbN) / CT));
      else CTN = Math.min(k.shapeMax, Math.max((TmlN - TmnN) / (TmlN - TbN), k.shapeMin));
    }
    TbN = Math.min(TbN, k.maximumDensityT);
  } else {
    const B = buoyancy(Tml) * Qstar / rc, wstar = B < 0 ? Math.cbrt(-B * h) : 0;
    let equilibrium = -Qw / Math.max(Iw, k.small);
    if (equilibrium > 0 && equilibrium < 1 && Tml > k.maximumDensityT) equilibrium = Math.min(D, (Math.sqrt(6 * equilibrium) + 2 * equilibrium / (1 - equilibrium)) / k.extinction);
    else equilibrium = 0;
    const N = h <= D - k.mixedMin ? Math.sqrt(Math.max(0, buoyancy(0.5 * (Tml + Tb)) * (Tml - Tb)) / (D - h)) : 0;
    let dCT = (k.shapeMax - k.shapeMin) / Math.max(N * (D - h) ** 2 / k.relaxC / Math.max(wstar, ustar, k.ustarMin) ** 2, k.small);
    const CTT = k.ctt1 * CT - k.ctt2, CQ = 2 * CTT / CT;
    let dh;
    if (B < 0) {
      CTN = Math.min(k.shapeMax, Math.max(CT + dCT * dt, k.shapeMin));
      dCT = (CTN - CT) / dt;
      if (h <= D - k.mixedMin) {
        if (h <= k.mixedMin) dh = k.cbl1 / k.cbl2 * Math.max(wstar, k.small);
        else {
          const RH = D / h, Rr = RH - 1, RT = CT / CTT, Rs = (RT / 2 - 1) * Rr + 1;
          dh = -Qstar * (Rs * (1 + k.cbl1) - 1) / rc;
          dh += (D - h) * (Tml - Tb) * k.ctt2 / CTT * dCT;
          dh += (Ibot + (RT - 1) * Ih - RT * IhD + (RT - 2) * Rr * (Ih - I0h)) / rc;
          dh /= -k.cbl2 * Rs * Qstar / rc / Math.max(wstar, k.small) + CT * (Tml - Tb);
        }
        dh = Math.max(dh, k.small);
        hN = Math.max(k.mixedMin, Math.min(h + dh * dt, D));
      } else hN = D;
    } else {
      const u = Math.max(ustar, k.ustarMin);
      let scale = Math.max((Math.abs(f) / k.zmN + N / k.zmI) * u * u + B / k.zmS, k.small);
      scale = Math.max(Math.max(k.mixedMin, Math.min(u ** 3 / scale, k.mixedMax)), equilibrium);
      hN = Math.max(k.mixedMin, Math.min(scale - (scale - h) * Math.exp(-k.relaxH * u / scale * dt), D));
      if (hN <= h) dCT = -dCT;
      CTN = Math.min(k.shapeMax, Math.max(CT + dCT * dt, k.shapeMin));
      dCT = (CTN - CT) / dt;
    }
    dh = (hN - h) / dt;
    if (hN <= D - k.mixedMin) {
      let dTb = 0;
      if (hN > h) {
        const RH = h / D, Rr = 1 - RH, R0 = 0.5 * CT * Rr + CTT * (2 * RH - 1), Rs = (0.5 + CTT - CQ) / R0, RT = (1 - CT * Rr) / R0;
        dTb = ((Qw + Iw - Ibot) / rc - CT * (Tml - Tb) * dh) * Rs / D;
        dTb += (IhD - (1 - CQ) * Ih - CQ * Ibot) * RT / (D - h) / rc;
        dTb += (1 - k.ctt2 * RT) / CT * (Tml - Tb) * dCT;
      }
      TbN = Math.max(k.freezing, Tb + dTb * dt);
      if ((TbN - k.maximumDensityT) * buoyancy(TmnN) < 0) TbN = k.maximumDensityT;
      const w = CTN * (1 - hN / D);
      TmlN = Math.max(k.freezing, (TmnN - TbN * w) / (1 - w));
    } else { hN = D; TmlN = TmnN; TbN = TmnN; CTN = k.shapeMin; }
  }
  if ((TmlN - TbN) * buoyancy(TmnN) < 0) { hN = D; TmlN = TmnN; TbN = TmnN; CTN = k.shapeMin; }
  c.mixed = TmlN; c.mean = TmnN; c.bottom = TbN; c.shape = CTN; c.mixedDepth = hN; c.ice = hiN; c.iceT = TiN;
  if (closed) lakeSettle(c, D, target);
  return c.ice > 0 ? c.iceT : c.mixed;
}

/*
 * Brings column `c` to enthalpy `target` (J/m², lakeEnthalpy's) through
 * its ice: the water keeps what the step gave it, and the ice is as thick
 * as the rest requires at its temperature, from E = −ρ_I h (L_f + c_I C_I(h)
 * (T_f − T_I)) with C_I(h) linear in h, the root of the quadratic that
 * goes to the latent-heat thickness as T_I nears T_f. With less ice than
 * iceMin the rest warms or cools the whole water; at iceMax the ice's
 * temperature takes it. Open water does not freeze on a deficit within
 * SETTLE_SLACK kelvin of the whole water, the round-off of f32's enthalpy.
 */
export const SETTLE_SLACK = 1e-4;
export function lakeSettle(c, D, target) {
  const k = FLAKE, rc = k.rhoWater * k.cWater;
  const left = target - rc * D * (c.mean - k.freezing), need = -left / k.rhoIce;
  if (need <= k.fusion * k.iceMin + (c.ice > 0 ? 0 : rc * D * SETTLE_SLACK / k.rhoIce)) {
    const shift = left / rc / D;
    c.mixed += shift; c.mean += shift; c.bottom += shift; c.ice = 0; c.iceT = k.freezing;
    return;
  }
  const deficit = c.ice > 0 ? k.freezing - c.iceT : 0, slope = k.iceMR * (1 + k.iceAst) / k.iceMax;
  const a = k.fusion + k.cIce * deficit * k.iceShape, b = k.cIce * deficit * slope, root = a * a - 4 * b * need;
  let h = root > 0 ? 2 * need / (a + Math.sqrt(root)) : k.iceMax + 1;
  if (h > k.iceMax) {
    h = k.iceMax;
    c.iceT = k.freezing - (need / h - k.fusion) / (k.cIce * (k.iceShape - slope * h));
  } else if (c.ice === 0) c.iceT = k.freezing;
  c.ice = h;
}

/*
 * The heat a column holds above the freezing point less the ice's latent
 * and sensible deficit, J/m².
 */
export function lakeEnthalpy(c, D) {
  const k = FLAKE, CI = k.iceShape - k.iceMR * (1 + k.iceAst) * Math.min(1, c.ice / k.iceMax);
  return k.rhoWater * k.cWater * D * (c.mean - k.freezing) - k.rhoIce * c.ice * (k.fusion + k.cIce * CI * (k.freezing - c.iceT));
}

export const LAKE_DEFAULTS = Object.freeze({ model: 'flake', freshDepth: 10, salineDepth: 5 });
export const lakeModelOf = (geography, { model } = {}) => model ?? (geography?.lake ? LAKE_DEFAULTS.model : 'off');
export const lakeStateName = (name) => `lake${name[0].toUpperCase()}${name.slice(1)}`;
export const lakeArrays = (serialized) => Object.fromEntries(Object.entries(serialized).map(([name, values]) => [lakeStateName(name), Float64Array.from(values)]));
export const savedLake = (saved) => (saved?.lakeMean ? Object.fromEntries([...LAKE_FIELDS, 'passed'].map((name) => [name, saved[lakeStateName(name)]]).filter(([, values]) => values)) : null);
export const LAKE_FIELDS = ['mixed', 'mean', 'bottom', 'shape', 'mixedDepth', 'ice', 'iceT'];

/*
 * The lakes of a geography (its `lake` cells, js/fineGeography.module.js)
 * under FLake: each its own column, as deep as the geography's lake depth
 * within the IFS's bounds, or, where ETOPO 2022 holds no bed below the
 * lake's level, the GLDB v3 default the IFS takes (Choulga et al. 2019):
 * freshDepth, or salineDepth for Natural Earth's alkaline lakes. The water
 * a lake gains (rain and snow and inflowing runoff, less evaporation) it
 * passes on whole to the sea cell its outlet reaches, so its depth stays
 * fixed and a closed basin's loss is drawn from that sea; `passed` holds
 * the running total per lake cell (kg/m², negative for a net loss). Snow
 * on a lake joins its water and takes its heat of fusion from it. FLake's
 * equations set the water; the ice holds what remains of the heat the
 * fluxes give (lakeSettle), so enthalpy() changes by exactly those.
 */
export function createLake(mesh, geography, { model = LAKE_DEFAULTS.model, freshDepth = LAKE_DEFAULTS.freshDepth, salineDepth = LAKE_DEFAULTS.salineDepth, fusion = 3.34e5, melting = 273.15 } = {}) {
  if (model !== 'flake') throw new Error(`lake model is 'flake' or 'off', not ${model}`);
  if (!geography || !geography.lake) throw new Error('the lake model needs a geography with lakes (mask n1024)');
  const C = mesh.nCells, cells = [];
  for (let i = 0; i < C; i++) if (geography.lake[i]) cells.push(i);
  const depth = new Float64Array(C), coriolis = new Float64Array(C), passed = new Float64Array(C);
  const fields = Object.fromEntries(LAKE_FIELDS.map((name) => [name, new Float64Array(C)]));
  for (const i of cells) {
    const bed = -geography.elevation[i];
    depth[i] = lakeDepth(bed >= 1 ? bed : geography.saline?.[i] ? salineDepth : freshDepth);
    coriolis[i] = 2 * 7.2921e-5 * Math.sin(mesh.latCell[i]);
  }
  const column = { mixed: 0, mean: 0, bottom: 0, shape: 0, mixedDepth: 0, ice: 0, iceT: 0 };
  const read = (i) => { for (const name of LAKE_FIELDS) column[name] = fields[name][i]; return column; };
  const write = (i) => { for (const name of LAKE_FIELDS) fields[name][i] = column[name]; };
  return {
    cells, depth, coriolis, ...fields, passed,
    isLake: geography.lake,
    initialize(surfaceT) {
      for (const i of cells) { Object.assign(column, lakeColumn(surfaceT[i], depth[i])); write(i); surfaceT[i] = column.mixed; }
    },
    step(i, surfaceT, ice, concentration, absorbed, flux, ustar, dt) {
      if (fields.mean[i] === 0) { Object.assign(column, lakeColumn(surfaceT[i], depth[i])); write(i); }
      const T = lakeStep(read(i), depth[i], absorbed, flux, ustar, coriolis[i], dt);
      write(i);
      surfaceT[i] = T; ice[i] = column.ice; concentration[i] = column.ice > 0 ? 1 : 0;
    },
    albedo(ice, surfaceT) { return ice > 0 ? lakeIceAlbedo(surfaceT) : null; },
    deposit(i, amount, airTemperature) {
      passed[i] += amount;
      const frozen = airTemperature < melting;
      if (frozen && fields.mean[i] !== 0) { const c = read(i); lakeSettle(c, depth[i], lakeEnthalpy(c, depth[i]) - fusion * amount); write(i); }
      return frozen;
    },
    evaporate(i, evaporation, dt) { passed[i] -= evaporation * dt; },
    enthalpy(i) { return lakeEnthalpy(read(i), depth[i]); },
    serialize() { return { ...Object.fromEntries(LAKE_FIELDS.map((name) => [name, Float64Array.from(fields[name])])), passed: Float64Array.from(passed) }; },
    load(saved) {
      for (const name of [...LAKE_FIELDS, 'passed']) if (saved[name]) (name === 'passed' ? passed : fields[name]).set(saved[name]);
    },
  };
}

/*
 * The geography the ocean takes when FLake holds the lakes: the lake cells
 * dry in its water (`ocean`, with the cavities' where there are any) and
 * land to its runoff, at their level, so what runs into a lake runs on
 * through it to the sea.
 */
export function lakeWater(mesh, geography) {
  const { nEdges: E, cellsOnEdge } = mesh, water = geography.ocean ?? geography, lake = geography.lake;
  const land = Uint8Array.from(water.land, (l, i) => (l || lake[i] ? 1 : 0)), edgeOcean = Uint8Array.from(water.edgeOcean), coast = [];
  for (let e = 0; e < E; e++) {
    const a = land[cellsOnEdge[2 * e]], b = land[cellsOnEdge[2 * e + 1]];
    if (a || b) edgeOcean[e] = 0;
    if (a !== b) coast.push(e);
  }
  const runoffLand = Uint8Array.from(geography.land, (l, i) => (l || lake[i] ? 1 : 0));
  const runoffElevation = Float64Array.from(geography.elevation, (h, i) => (lake[i] ? geography.lakeLevel[i] : h));
  return { ...geography, ocean: { ...water, land, edgeOcean, coastEdges: Int32Array.from(coast) }, runoffLand, runoffElevation };
}
