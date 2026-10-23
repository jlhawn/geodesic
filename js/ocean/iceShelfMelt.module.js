/*
 * Melt at the base of a floating ice shelf by the three-equation boundary
 * condition of Holland & Jenkins (1999), with the constants of ISOMIP+
 * (Asay-Davis et al. 2016, eqs. 24–27 and Table 4):
 *
 *   m L = ρ0 c γ_T (T_w − T_b)        heat
 *   m S_b = ρ0 γ_S (S_w − S_b)        salt
 *   T_b = λ1 S_b + λ2 + λ3 p_b        the interface at its freezing point
 *
 * m is the melt in kg/m²/s (negative when water freezes onto the base),
 * γ = Γ u* with u*² = C_d (|u_w|² + u_tidal²), and T_w, S_w, u_w are the
 * water's at the ice base. Temperatures are potential: the model's are,
 * and ISOMIP+'s liquidus coefficients are fitted to the potential
 * freezing point. The ice conducts no heat (κ_i = 0, as in ISOMIP+).
 *
 * In δ = S_w − S_b the balances are a δ² + b δ + c = 0, with
 *   a = −A λ1,  b = A (λ1 S_w − T*) − B,  c = A T* S_w,
 * A = ρ0 c Γ_T u* / L, B = ρ0 Γ_S u* and T* = T_w − (λ1 S_w + λ2 + λ3 p_b).
 * The root that vanishes with T* is c / q, q = −(b + sgn(b) √(b² − 4ac)) / 2,
 * and m = B δ / S_b: neither subtracts near-equal numbers as the melt goes
 * to zero, so f32 keeps the root there.
 */
export const SHELF_MELT = Object.freeze({
  liquidusSlope: -0.0573, liquidusConstant: 0.0832, liquidusPressure: -7.53e-8,
  heatTransfer: 0.011, saltTransfer: 3.1e-4, tidalSpeed: 0.01, freshDensity: 1000,
});

export const baseFreezing = (pressure, zero = 273.15) => zero + SHELF_MELT.liquidusConstant + SHELF_MELT.liquidusPressure * pressure;

export function basalMelt(T, S, freezing, ustar, { rho0, rhoCp, latentHeat, heatTransfer = SHELF_MELT.heatTransfer, saltTransfer = SHELF_MELT.saltTransfer }) {
  const slope = SHELF_MELT.liquidusSlope, A = rhoCp * heatTransfer * ustar / latentHeat, B = rho0 * saltTransfer * ustar;
  const driving = T - (slope * S + freezing);
  const a = -A * slope, b = A * (slope * S - driving) - B, c = A * driving * S;
  const q = -0.5 * (b + (b > 0 ? 1 : -1) * Math.sqrt(Math.max(0, b * b - 4 * a * c)));
  const fresher = q !== 0 ? c / q : 0, Sb = S - fresher;
  return { melt: B * fresher / Sb, Tb: slope * Sb + freezing, Sb };
}

export const basalMeltWgsl = ({ rho0, rhoCp, latentHeat, heatTransfer = SHELF_MELT.heatTransfer, saltTransfer = SHELF_MELT.saltTransfer }) => `fn basalMelt(t: f32, s: f32, freezing: f32, ustar: f32) -> vec3<f32> {
  let slope = f32(${SHELF_MELT.liquidusSlope}); let A = f32(${rhoCp * heatTransfer / latentHeat}) * ustar; let B = f32(${rho0 * saltTransfer}) * ustar;
  let driving = t - (slope * s + freezing);
  let a = -A * slope; let b = A * (slope * s - driving) - B; let c = A * driving * s;
  let q = -0.5 * (b + select(-1.0, 1.0, b > 0.0) * sqrt(max(0.0, b * b - 4.0 * a * c)));
  let fresher = select(0.0, c / q, q != 0.0); let sb = s - fresher;
  return vec3<f32>(B * fresher / sb, slope * sb + freezing, sb);
}
`;
