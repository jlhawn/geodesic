/*
 * Turbulent orographic form drag (Beljaars, Brown and Wood 2004) as the
 * IFS documents it (Cy47r3 Part IV §3.4, eqs. 3.55–3.57): the drag of
 * hills below 5 km, spread with height,
 *   ∂U/∂t = −C_tofd(z) |U| U,
 *   C_tofd = α β C_md C_corr 2.109 e^(−(z/1500)^1.5) a₂ z^(−1.2),
 *   a₂ = a₁ k₁^(n₁ − n₂),  a₁ = σ_flt² (I_H k_flt^n₁)⁻¹,
 * z the height above the model's ground in metres and σ_flt the standard
 * deviation of the orography in the 3–22 km band (§11.3.3). The boundary
 * layer's edge solve takes it implicitly with |U| from the step's start.
 */
export const FORM_DRAG_DEFAULTS = { alpha: 35, beta: 1, mountainDrag: 0.005, correction: 0.6, n1: -1.9, n2: -2.8, k1: 0.003, kFilter: 0.00035, filterIntegral: 0.00102, decayHeight: 1500 };

export function formDragScale(o = FORM_DRAG_DEFAULTS) {
  const c = { ...FORM_DRAG_DEFAULTS, ...o };
  return c.alpha * c.beta * c.mountainDrag * c.correction * 2.109 * Math.pow(c.k1, c.n1 - c.n2) / (c.filterIntegral * Math.pow(c.kFilter, c.n1));
}

export function formDragCoefficient(sigmaFiltered, z, o = FORM_DRAG_DEFAULTS, scale = formDragScale(o)) {
  if (!(sigmaFiltered > 0) || !(z > 0)) return 0;
  const decay = o.decayHeight ?? FORM_DRAG_DEFAULTS.decayHeight;
  return scale * sigmaFiltered * sigmaFiltered * Math.exp(-Math.pow(z / decay, 1.5)) * Math.pow(z, -1.2);
}
