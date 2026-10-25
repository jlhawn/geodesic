import { FREEZING_POINT, MELTING_POINT, MINIMUM_CONCENTRATION } from '../physics/ice.module.js';

/*
 * The drift of sea ice on the ocean's own C-grid: one velocity component
 * along the normal of every ocean edge, stepped once per ocean call with
 * the wind stress the atmosphere hands to the ocean, and the ice's area A,
 * volume per cell area V = A h, snow mass, skin enthalpy and snow albedo
 * carried by it. The thermodynamics (physics/ice.module.js) never moves
 * ice; this module does, and keeps the ocean's salt, heat and stress
 * consistent with the move.
 *
 * Momentum, per unit ice area, for the mass μ = ρ_i h + s of ice and snow:
 *   μ du/dt = τ_a + τ_w − μ f k×u − μ g ∇η + (1/A) ∇·σ,
 * with τ_a the air stress (the edge stress the atmosphere gives the ocean),
 * τ_w = ρ_w C_w |U_w − u| R(θ)(U_w − u) the water stress on the ice
 * (McPhee 1980; ρ_w C_w with C_w = 5.5e-3 is MITgcm's default and about twice
 * the 10 m air drag), R(θ) a turning of the relative velocity by θ counter-
 * clockwise in the Northern Hemisphere and clockwise in the Southern (the
 * form of CICE's ice-ocean stress; θ = 0 leaves the Ekman turning to the
 * ocean's 3 m top layer, whose velocity is U_w), η the ocean's free surface,
 * and σ the vertically integrated internal stress. The advection of ice
 * momentum is dropped (|u| is of order 0.1 m/s). The ice reaches the ocean as
 * τ_io = −τ_w, so that the stress the ocean feels is (1 − Ā) τ_a + Ā τ_io
 * for the mean cover Ā of an edge's two cells.
 *
 * Rheologies. 'free' is free drift: σ = 0, which with the water drag and the
 * Coriolis force gives the steady solution of McPhee (2008) and Leppäranta
 * (2011). 'cavitating' is the cavitating fluid of Flato and Hibler (1992):
 * the ice resists convergence only, with the strength of Hibler (1979),
 * P = P* V exp(−C (1 − A)) (P* = 2.75e4 N/m², C = 20), carried by one
 * isotropic stress σ1 = σ11 + σ22 per cell whose equilibrium is
 * 2 P min(D, 0) / max(|D|, Δ_min) for the divergence D of the ice velocity,
 * the replacement-pressure form of Hibler's law with the elliptic yield curve
 * opened to a half line; the force on an edge is the TRiSK gradient of σ1/2
 * (Thuburn et al. 2009, Ringler et al. 2010) over the edge's mean cover, which
 * needs no tangential reconstruction and no strain tensor. 'frozen' keeps the
 * velocity it holds (zero unless a field is prescribed), takes no stress from
 * it and leaves the ocean's stress to iceStressTransmission, so that the
 * ocean runs as it does with the drift off.
 *
 * Solver. The modified EVP of Bouillon et al. (2013, eqs. 11–23), read as a
 * pseudo-time iteration p towards the viscous-plastic solution
 * (Kimmritz et al. 2015; Danilov et al. 2015, eqs. 13–16): each iteration
 * relaxes σ1 towards its equilibrium at the latest velocity by 1/α, then
 *   β (u^{p+1} − u^p) + (u^{p+1} − u^n) = (Δt/μ) [ τ_a + τ_w + μ f u_t − μ g ∂η/∂n + F ]
 * with u_t the TRiSK tangential component of u^p (the C-grid form of
 * Bouillon et al. 2013 eqs. 22–23, which also lags the other component), the
 * water drag implicit in u^{p+1} and its magnitude |U_w − u^p| lagged. The
 * Coriolis force is the TRiSK sum with ½(f_e + f_e') of the ocean's
 * edgeCoriolis, so that it does no work. The pseudo-time damping is raised to
 * the drag number a = (Δt/μ) ρ_w C_w |U_w − u^p| cos θ wherever that is larger
 * than β, which leaves the fixed point alone and makes the iteration of a thin
 * or free-drifting edge contract (cf. the adaptive EVP of Kimmritz et al.
 * 2016); with β = 0 and no stress it converges in a few passes. Stability
 * wants α β of at least 3 P Δt / (Δ_min m d_c²) for the stiffest ice, the
 * compact pack, where P/m = P* / ρ_i (Bouillon et al. 2013 §4; Kimmritz et al.
 * 2015; Danilov et al. 2015: the same ratio with a larger constant), and the
 * residual falls as exp(−p/β) where the stress is small and exp(−p/2α) where
 * it is large. With alpha, beta and iterations left null, α = β is 1.5 times
 * the limit sqrt(3 P* Δt / (ρ_i Δ_min d_c²)) of the mean spacing d_c and the
 * iteration count 3 α, which brings a pack's first call from rest within 2
 * percent of the converged solution (test/iceDrift.test.mjs). σ1 carries over
 * from call to call. An edge whose mean cover is below dynamicConcentration
 * rides the water.
 *
 * Transport is flux form with the ocean's Koren-limited third-order edge
 * values (limitedEdgeValue, upwind where the upwind-biased line meets land or
 * a cavity or an empty cell): the area flux u A* carries the thickness, the
 * snow per ice area, the skin deficit C_s (T_skin − T_f) per ice area and the
 * snow albedo, so that those stay inside the range of their neighbors and the
 * volume, snow mass and skin enthalpy are conserved. Forward Euler is
 * subcycled when the Courant number passes courantLimit.
 *
 * Closure. A cell whose area passes maximumConcentration is set to it at
 * constant volume, snow mass and skin enthalpy, so the excess area becomes
 * thickness (the two-category closure of Hibler 1979). Ice left below the
 * sea ice's removal floors (seaIce.removal, the rule of ice.update) melts
 * away through seaIce.meltAway, which draws its latent heat from the water.
 * With keepSmallIce the floors are NEGLIGIBLE_CONCENTRATION and
 * NEGLIGIBLE_VOLUME, so that the few thousandths of a cell that a Courant
 * number of 0.005 carries into open water survive to freeze on, or to melt
 * in the water's heat. The limited edge values lie between the donor's and
 * the acceptor's, so the thickness, snow per ice area and skin temperature
 * of a nearly emptied cell, quotients of small numbers, are not bounded by
 * its neighbors'; a cell under MINIMUM_CONCENTRATION whose quotients have
 * left ICE_REMNANT melts away as well, which keeps them finite in single
 * precision.
 *
 * Coupling with the ocean's bookkeeping. The ocean's virtual salt flux follows
 * the change of each cell's ice volume, so `previousVolume` is shifted by the
 * volume that transport moved, leaving the freezing and melting for salt().
 * A cell the ice reaches while open keeps its water's temperature in `waterT`
 * (flagged in `arrived`); the ocean's readSurface holds the top layer of such
 * a cell at that temperature for the call, where it holds every other iced
 * cell at the freezing point, so that the water's heat melts the arriving ice
 * through the ocean flux. Ice edges are open only where both cells are open
 * sea: not land, not an ice-shelf cavity.
 */
export const ICE_RHEOLOGIES = ['frozen', 'free', 'cavitating'], ICE_HEMISPHERES = ['both', 'south', 'north'];
export const COUPLED_ICE_DRIFT = { rheology: 'cavitating', hemisphere: 'both' };
export const ICE_DRIFT = {
  rheology: 'cavitating', hemisphere: 'both',
  strength: 27.5e3, strengthDecay: 20, waterDrag: 5.5e-3, waterTurn: 0,
  maximumConcentration: 1, minimumDeformation: 1e-8,
  alpha: null, beta: null, iterations: null, freePasses: 6, dynamicsEvery: 1,
  speedLimit: 2, dynamicConcentration: MINIMUM_CONCENTRATION, minimumRelativeSpeed: 1e-3, courantLimit: 0.2,
};

export const ICE_REMNANT = { thickness: 50, snow: 5e3, skinLow: 150, skinHigh: MELTING_POINT + 1 };
export const iceDriftMoves = (option) => !!option && (option === true || (option.rheology ?? ICE_DRIFT.rheology) !== 'frozen');

export const coverOf = (ice, concentration, i) => (ice[i] > 0 ? (concentration && concentration[i] > 0 ? concentration[i] : 1) : 0);

export function iceDriftOptions(options) {
  if (options !== true && (options === null || typeof options !== 'object')) throw new Error(`iceDrift is false or an object of options, not ${options}`);
  const o = { ...ICE_DRIFT, ...(options === true ? {} : options) };
  for (const key of Object.keys(o)) if (!(key in ICE_DRIFT)) throw new Error(`iceDrift has no option ${key}; its options are ${Object.keys(ICE_DRIFT).join(', ')}`);
  if (!ICE_RHEOLOGIES.includes(o.rheology)) throw new Error(`iceDrift.rheology is one of ${ICE_RHEOLOGIES.join(', ')}, not ${o.rheology}`);
  if (!ICE_HEMISPHERES.includes(o.hemisphere)) throw new Error(`iceDrift.hemisphere is one of ${ICE_HEMISPHERES.join(', ')}, not ${o.hemisphere}`);
  for (const key of ['strength', 'minimumDeformation', 'speedLimit', 'courantLimit', 'minimumRelativeSpeed', 'maximumConcentration']) if (!(o[key] > 0)) throw new Error(`iceDrift.${key} is positive, not ${o[key]}`);
  for (const key of ['strengthDecay', 'waterDrag', 'dynamicConcentration']) if (!(o[key] >= 0)) throw new Error(`iceDrift.${key} is not negative, not ${o[key]}`);
  if (o.alpha !== null && !(o.alpha > 0)) throw new Error(`iceDrift.alpha is positive or null, not ${o.alpha}`);
  if (o.beta !== null && !(o.beta >= 0)) throw new Error(`iceDrift.beta is not negative or null, not ${o.beta}`);
  if (o.iterations !== null && !(Number.isInteger(o.iterations) && o.iterations >= 1)) throw new Error(`iceDrift.iterations is a whole number of at least 1 or null, not ${o.iterations}`);
  for (const key of ['freePasses', 'dynamicsEvery']) if (!(Number.isInteger(o[key]) && o[key] >= 1)) throw new Error(`iceDrift.${key} is a whole number of at least 1, not ${o[key]}`);
  if (!(Math.abs(o.waterTurn) < 90)) throw new Error(`iceDrift.waterTurn is a turning angle in degrees under 90, not ${o.waterTurn}`);
  return o;
}

export function iceDriftPasses(o, dt, iceDensity, spacing) {
  const stable = Math.sqrt(3 * o.strength * dt / (iceDensity * o.minimumDeformation * spacing * spacing));
  const alpha = o.alpha ?? Math.max(4, Math.ceil(1.5 * stable)), cavitating = o.rheology === 'cavitating';
  return { alpha, beta: cavitating ? (o.beta ?? alpha) : 0, passes: cavitating ? (o.iterations ?? Math.ceil(3 * alpha)) : o.freePasses };
}

export function createIceDrift(mesh, options, { cellOpen, edgeOcean, waterDensity, iceDensity, gravity, previousVolume, beyond = null, limited = null }) {
  const o = iceDriftOptions(options);
  const {
    nCells: C, nEdges: E, maxEdges, nEdgesOnCell, edgesOnCell, edgeSignOnCell, cellsOnEdge, areaCell, dcEdge, dvEdge, fEdge, latEdge, latCell,
    maxEdgesOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge,
  } = mesh;
  if (beyond && !limited) throw new Error('createIceDrift needs the limited edge value that goes with the cells beyond');
  const cavitating = o.rheology === 'cavitating', moving = o.rheology !== 'frozen';
  const cosTurn = Math.cos(o.waterTurn * Math.PI / 180), sinTurn = Math.sin(o.waterTurn * Math.PI / 180), drag = waterDensity * o.waterDrag;
  let spacing = 0;
  for (let e = 0; e < E; e++) spacing += dcEdge[e] / E;
  const dynamic = new Uint8Array(E);
  for (let e = 0; e < E; e++) {
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
    const side = o.hemisphere === 'both' || (o.hemisphere === 'south' ? latEdge[e] < 0 : latEdge[e] > 0);
    dynamic[e] = edgeOcean[e] && cellOpen[a] && cellOpen[b] && side ? 1 : 0;
  }

  const velocity = new Float64Array(E), spare = new Float64Array(E), start = new Float64Array(E), sigma = new Float64Array(C);
  const iceStress = new Float64Array(E), volumeFlux = new Float64Array(E), valid = new Uint8Array(E);
  const waterT = new Float64Array(C), arrived = new Uint8Array(C);
  const A = new Float64Array(C), V = new Float64Array(C), Sn = new Float64Array(C), En = new Float64Array(C), Sa = new Float64Array(C), mass = new Float64Array(C), pressure = new Float64Array(C);
  const hi = new Float64Array(C), si = new Float64Array(C), di = new Float64Array(C), ai = new Float64Array(C);
  const dA = new Float64Array(C), dV = new Float64Array(C), dS = new Float64Array(C), dE = new Float64Array(C), dSa = new Float64Array(C);
  const Vstart = new Float64Array(C), Astart = new Float64Array(C), touched = new Uint8Array(C), touchedList = new Int32Array(C), cellList = new Int32Array(C);
  const Ae = new Float64Array(E), mu = new Float64Array(E), uw = new Float64Array(E), uwt = new Float64Array(E), tilt = new Float64Array(E);
  const act = new Int32Array(E), upd = new Int32Array(E);
  let nAct = 0, nUpd = 0, nCellsActive = 0, nTouched = 0, calls = 0, ice = null, removal = null, enthalpyLow = 0, enthalpyHigh = 0;
  const count = { clamps: 0, fallbacks: 0, closures: 0, closedArea: 0, courant: 0, thickest: 0, topSpeed: 0, removed: 0 };

  function attach(seaIce) {
    if (!seaIce || !seaIce.snow || !seaIce.meltAway || !seaIce.removal) throw new Error('the ice drift needs the sea ice of physics/ice.module.js: ocean.useSeaIce(concentration, seaIce)');
    ice = seaIce;
    removal = seaIce.removal;
    enthalpyLow = seaIce.skinHeatCapacity * (ICE_REMNANT.skinLow - FREEZING_POINT);
    enthalpyHigh = seaIce.skinHeatCapacity * (ICE_REMNANT.skinHigh - FREEZING_POINT);
  }
  const remnant = (a, v, s, e) => v > ICE_REMNANT.thickness * a || s > ICE_REMNANT.snow * a || e < enthalpyLow * a || e > enthalpyHigh * a;
  function reset(saved = null, savedSigma = null) {
    velocity.fill(0); spare.fill(0); sigma.fill(0); iceStress.fill(0); volumeFlux.fill(0); valid.fill(0); arrived.fill(0); calls = 0;
    if (saved && saved.length === E) velocity.set(saved);
    if (savedSigma && savedSigma.length === C) sigma.set(savedSigma);
  }

  function prepare({ surfaceT, ice: thickness, concentration }) {
    if (!ice) throw new Error('the ice drift has no sea ice: call ocean.useSeaIce(concentration, seaIce) first');
    if (!concentration) throw new Error('the ice drift moves the sea ice concentration, which the ocean was not given');
    const { snow, snowAlbedo, skinHeatCapacity, freshSnowAlbedo } = ice;
    nCellsActive = 0;
    for (let i = 0; i < C; i++) {
      const a = cellOpen[i] ? coverOf(thickness, concentration, i) : 0;
      if (a > 0) {
        A[i] = a; V[i] = a * thickness[i]; Sn[i] = a * snow[i]; En[i] = a * skinHeatCapacity * (surfaceT[i] - FREEZING_POINT); Sa[i] = Sn[i] * snowAlbedo[i];
        mass[i] = iceDensity * V[i] + Sn[i];
        cellList[nCellsActive++] = i;
      } else { A[i] = 0; V[i] = 0; Sn[i] = 0; En[i] = 0; Sa[i] = 0; mass[i] = 0; sigma[i] = 0; }
    }
    nAct = 0;
    for (let e = 0; e < E; e++) {
      if (!dynamic[e]) continue;
      if (A[cellsOnEdge[2 * e]] + A[cellsOnEdge[2 * e + 1]] > 0) act[nAct++] = e;
      else if (moving) velocity[e] = 0;
    }
  }

  const tangentialOf = (field, e) => {
    let sum = 0;
    const base = maxEdgesOnEdge * e;
    for (let s = 0; s < nEdgesOnEdge[e]; s++) { const other = edgesOnEdge[base + s]; sum += weightsOnEdge[base + s] * dvEdge[other] * field[other]; }
    return sum / dcEdge[e];
  };

  function iterate({ stress, waterU, eta }, dt) {
    valid.fill(0);
    if (!moving) return;
    const { clamps } = count;
    nUpd = 0;
    for (let n = 0; n < nAct; n++) {
      const e = act[n], a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], area = A[a] + A[b];
      Ae[e] = 0.5 * area; mu[e] = (mass[a] + mass[b]) / area;
      uw[e] = waterU[e]; uwt[e] = tangentialOf(waterU, e); tilt[e] = gravity * (eta[b] - eta[a]) / dcEdge[e];
      if (Ae[e] < o.dynamicConcentration) velocity[e] = uw[e]; else { upd[nUpd++] = e; start[e] = velocity[e]; }
    }
    spare.set(velocity);
    const { alpha, beta: betaBase, passes } = iceDriftPasses(o, dt, iceDensity, spacing);
    const Dmin = o.minimumDeformation, limit = o.speedLimit, floor = o.minimumRelativeSpeed;
    count.alpha = alpha; count.beta = betaBase; count.passes = passes;
    if (cavitating) for (let n = 0; n < nCellsActive; n++) { const i = cellList[n]; pressure[i] = o.strength * V[i] * Math.exp(-o.strengthDecay * (1 - A[i])); }
    let cur = velocity, nxt = spare, clamped = 0;
    for (let p = 0; p < passes; p++) {
      if (cavitating) {
        for (let n = 0; n < nCellsActive; n++) {
          const i = cellList[n], base = maxEdges * i;
          let div = 0;
          for (let k = 0; k < nEdgesOnCell[i]; k++) { const e = edgesOnCell[base + k]; div += edgeSignOnCell[base + k] * dvEdge[e] * cur[e]; }
          div /= areaCell[i];
          sigma[i] += ((div < 0 ? 2 * pressure[i] * div / Math.max(-div, Dmin) : 0) - sigma[i]) / alpha;
        }
      }
      const last = p === passes - 1;
      for (let n = 0; n < nUpd; n++) {
        const e = upd[n], a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], fe = fEdge[e], base = maxEdgesOnEdge * e;
        let tang = 0, cor = 0;
        for (let s = 0; s < nEdgesOnEdge[e]; s++) {
          const other = edgesOnEdge[base + s], w = weightsOnEdge[base + s] * dvEdge[other] * cur[other];
          tang += w; cor += w * (fe + fEdge[other]);
        }
        tang /= dcEdge[e]; cor *= 0.5 / dcEdge[e];
        const dn = uw[e] - cur[e], dtn = uwt[e] - tang;
        const dr = drag * Math.max(Math.sqrt(dn * dn + dtn * dtn), floor), sinT = fe > 0 ? sinTurn : fe < 0 ? -sinTurn : 0;
        const k = dt / mu[e], implicit = k * dr * cosTurn, bet = Math.max(betaBase, implicit);
        const force = cavitating ? (sigma[b] - sigma[a]) / (2 * dcEdge[e] * Ae[e]) : 0;
        const rhs = stress[e] + dr * (uw[e] * cosTurn - dtn * sinT) + mu[e] * (cor - tilt[e]) + force;
        let un = (bet * cur[e] + start[e] + k * rhs) / (1 + bet + implicit);
        if (un > limit) { un = limit; if (last) clamped++; } else if (un < -limit) { un = -limit; if (last) clamped++; }
        nxt[e] = un;
      }
      const swap = cur; cur = nxt; nxt = swap;
    }
    if (cur !== velocity) velocity.set(cur);
    count.clamps = clamps + clamped;
    let top = 0;
    for (let n = 0; n < nAct; n++) {
      const e = act[n], fe = fEdge[e];
      const dn = uw[e] - velocity[e], dtn = uwt[e] - tangentialOf(velocity, e), dr = drag * Math.max(Math.sqrt(dn * dn + dtn * dtn), floor);
      iceStress[e] = -dr * (dn * cosTurn - dtn * (fe > 0 ? sinTurn : fe < 0 ? -sinTurn : 0));
      valid[e] = 1;
      top = Math.max(top, Math.abs(velocity[e]));
    }
    count.topSpeed = top;
  }

  function carry({ surfaceT, ice: thickness, concentration }, dt) {
    const { skinHeatCapacity, freshSnowAlbedo, snow, snowAlbedo, latent, latentHeatFusion, meltAway } = ice;
    let courant = 0;
    for (let n = 0; n < nAct; n++) { const e = act[n]; courant = Math.max(courant, Math.abs(velocity[e]) * dt / dcEdge[e]); }
    count.courant = courant; count.removed = 0;
    volumeFlux.fill(0);
    if (!(courant > 0)) return;
    const sub = Math.ceil(courant / o.courantLimit), step = dt / sub;
    nTouched = 0;
    for (let s = 0; s < sub; s++) {
      for (let i = 0; i < C; i++) {
        const a = A[i];
        if (a > 0) { hi[i] = V[i] / a; si[i] = Sn[i] / a; di[i] = En[i] / a; ai[i] = Sn[i] > 0 ? Sa[i] / Sn[i] : freshSnowAlbedo; } else { hi[i] = 0; si[i] = 0; di[i] = 0; ai[i] = freshSnowAlbedo; }
      }
      for (let n = 0; n < nAct; n++) {
        const e = act[n], u = velocity[e];
        if (u === 0) continue;
        const side = u > 0 ? 0 : 1, a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], up = side === 0 ? a : b, down = side === 0 ? b : a;
        if (!(A[up] > 0)) continue;
        let Astar = A[up], hs = hi[up], ss = si[up], ds = di[up], as = ai[up];
        if (beyond) {
          const uu = beyond.cell[2 * e + side];
          if (cellOpen[uu]) {
            const dUp = dcEdge[beyond.edge[2 * e + side]], dCentral = dcEdge[e];
            Astar = Math.max(0, limited(A, up, down, uu, dUp, dCentral));
            if (A[down] > 0 && A[uu] > 0) {
              hs = limited(hi, up, down, uu, dUp, dCentral); ss = limited(si, up, down, uu, dUp, dCentral);
              ds = limited(di, up, down, uu, dUp, dCentral); as = limited(ai, up, down, uu, dUp, dCentral);
            }
          } else count.fallbacks++;
        }
        const fA = u * Astar, fV = fA * hs, fS = fA * ss, fE = fA * ds, fSa = fS * as, scale = step * dvEdge[e], ia = scale / areaCell[a], ib = scale / areaCell[b];
        volumeFlux[e] += fV / sub;
        dA[a] -= fA * ia; dV[a] -= fV * ia; dS[a] -= fS * ia; dE[a] -= fE * ia; dSa[a] -= fSa * ia;
        dA[b] += fA * ib; dV[b] += fV * ib; dS[b] += fS * ib; dE[b] += fE * ib; dSa[b] += fSa * ib;
        if (!touched[a]) { touched[a] = 1; touchedList[nTouched++] = a; Vstart[a] = V[a]; Astart[a] = A[a]; }
        if (!touched[b]) { touched[b] = 1; touchedList[nTouched++] = b; Vstart[b] = V[b]; Astart[b] = A[b]; }
      }
      for (let n = 0; n < nTouched; n++) {
        const i = touchedList[n];
        A[i] = Math.max(0, A[i] + dA[i]); V[i] = Math.max(0, V[i] + dV[i]); Sn[i] = Math.max(0, Sn[i] + dS[i]); En[i] += dE[i]; Sa[i] = Math.max(0, Sa[i] + dSa[i]);
        dA[i] = 0; dV[i] = 0; dS[i] = 0; dE[i] = 0; dSa[i] = 0;
      }
    }
    for (let n = 0; n < nTouched; n++) {
      const i = touchedList[n];
      touched[i] = 0;
      let A1 = A[i];
      const V1 = V[i], S1 = Sn[i], E1 = En[i], Sa1 = Sa[i];
      previousVolume[i] += V1 - Vstart[i];
      if (A1 > o.maximumConcentration) { count.closedArea += (A1 - o.maximumConcentration) * areaCell[i]; count.closures++; A1 = o.maximumConcentration; }
      if (A1 < removal.concentration || V1 < removal.volume || (A1 < MINIMUM_CONCENTRATION && remnant(A1, V1, S1, E1))) {
        count.removed += areaCell[i] * V1;
        meltAway(surfaceT, thickness, i, E1 - latent * V1 - latentHeatFusion * S1, areaCell[i] * V1, areaCell[i] * S1, Astart[i] > 0 ? FREEZING_POINT : surfaceT[i]);
        continue;
      }
      if (Astart[i] === 0) { waterT[i] = surfaceT[i]; arrived[i] = 1; }
      thickness[i] = V1 / A1; concentration[i] = A1; snow[i] = S1 / A1; surfaceT[i] = FREEZING_POINT + E1 / (A1 * skinHeatCapacity); snowAlbedo[i] = S1 > 0 ? Math.min(1, Sa1 / S1) : freshSnowAlbedo;
      count.thickest = Math.max(count.thickest, thickness[i]);
    }
  }

  function solve(env) { prepare(env); iterate(env, env.dt); }
  function transport(env) { prepare(env); carry(env, env.dt); }
  function advance(env) {
    arrived.fill(0);
    if (++calls % o.dynamicsEvery !== 0) return false;
    const dt = env.dt * o.dynamicsEvery;
    prepare(env);
    iterate(env, dt);
    carry(env, dt);
    return true;
  }

  function transportAcross(latitude) {
    let sum = 0;
    for (let n = 0; n < nAct; n++) {
      const e = act[n], la = latCell[cellsOnEdge[2 * e]], lb = latCell[cellsOnEdge[2 * e + 1]];
      if ((la < latitude) === (lb < latitude)) continue;
      sum += volumeFlux[e] * dvEdge[e] * (lb > la ? 1 : -1);
    }
    return sum;
  }
  function diagnostics() {
    const circle = 60 * Math.PI / 180;
    return {
      iceSpeed: count.topSpeed, iceSpeedClamps: count.clamps, iceLimiterFallbacks: count.fallbacks, iceClosures: count.closures, iceClosedArea: count.closedArea, iceThickest: count.thickest, iceRemoved: count.removed,
      iceTransport60S: transportAcross(-circle), iceTransport60N: transportAcross(circle),
    };
  }

  return { options: o, rheology: o.rheology, velocity, sigma, waterT, arrived, valid, iceStress, volumeFlux, dynamic, count, attach, reset, solve, transport, advance, transportAcross, diagnostics };
}
