/*
 * Seawater density from the simplified equation of state of Roquet et al.
 * (2015), with NEMO's nn_eos = 1 coefficients, at the surface (potential
 * density):
 *   ρ = ρ₀ − a₀ (1 + ½ λ₁ Tₐ) Tₐ + b₀ (1 − ½ λ₂ Sₐ) Sₐ − ν Tₐ Sₐ,
 *   Tₐ = T − 10 °C, Sₐ = S − 35.
 * The quadratic temperature term makes cold water far less sensitive to
 * temperature than warm water (α ≈ 6×10⁻⁵ /K at 0 °C against 3×10⁻⁴ at
 * 25 °C), so fresh water at the freezing point floats on warmer, saltier
 * deep water as it does in the polar oceans.
 */
export const SEAWATER = { rho0: 1026, a0: 0.1655, b0: 0.76554, lambda1: 0.05952, lambda2: 7.4914e-4, nu: 2.4341e-3, t0: 283.15, s0: 35 };

const { rho0, a0, b0, lambda1, lambda2, nu, t0, s0 } = SEAWATER;

export function seawaterDensity(t, s) {
  const ta = t - t0, sa = s - s0;
  return rho0 - a0 * (1 + 0.5 * lambda1 * ta) * ta + b0 * (1 - 0.5 * lambda2 * sa) * sa - nu * ta * sa;
}

export function thermalExpansion(t, s) {
  return (a0 * (1 + lambda1 * (t - t0)) + nu * (s - s0)) / rho0;
}

/*
 * The temperature at which water of salinity s has density rho, on the
 * warm side of the density maximum; a density beyond the maximum gives
 * the temperature of the maximum.
 */
export function labelTemperature(rho, s = s0) {
  const sa = s - s0;
  const A = 0.5 * a0 * lambda1, B = a0 + nu * sa, Cc = rho - rho0 - b0 * (1 - 0.5 * lambda2 * sa) * sa;
  const discriminant = B * B - 4 * A * Cc;
  return t0 + (discriminant > 0 ? (-B + Math.sqrt(discriminant)) / (2 * A) : -B / (2 * A));
}

export const SEAWATER_WGSL = `
fn eosAnomaly(t: f32, s: f32) -> f32 {
  let ta = t - ${t0}; let sa = s - ${s0.toFixed(1)};
  return -${a0} * (1.0 + ${0.5 * lambda1} * ta) * ta + ${b0} * (1.0 - ${0.5 * lambda2} * sa) * sa - ${nu} * ta * sa;
}
fn eos(t: f32, s: f32) -> f32 { return ${rho0.toFixed(1)} + eosAnomaly(t, s); }
fn alphaT(t: f32, s: f32) -> f32 { return (${a0} * (1.0 + ${lambda1} * (t - ${t0})) + ${nu} * (s - ${s0.toFixed(1)})) / ${rho0.toFixed(1)}; }
`;
