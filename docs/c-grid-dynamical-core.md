# WebGCM: C-Grid Dynamical Core on the ISEA Icosahedral Mesh

Design document for replacing the climate model's A-grid horizontal
dynamical core with a C-grid (TRiSK) formulation built on
`js/grid.module.js`. Status: **M0 complete** (mesh, operators, TRiSK
weights, all cross-checked against MPAS and the papers); M1 next. Every
operator carries a test that catches sign and convention errors
independently of anyone's memory of the formula.

References:

- Thuburn, Ringler, Skamarock & Klemp (2009), *Numerical representation of
  geostrophic modes on arbitrarily structured C-grids*, J. Comput. Phys. 228.
- Ringler, Thuburn, Klemp & Skamarock (2010), *A unified approach to energy
  conservation and potential vorticity dynamics for arbitrarily-structured
  C-grids*, J. Comput. Phys. 229. ("RTSK" below; the MPAS formulation.)
- Williamson, Drake, Hack, Jakob & Swarztrauber (1992), the standard
  shallow-water test suite.
- Galewsky, Scott & Polvani (2004), barotropic-instability test case.
- MPAS mesh specification (variable names below follow it where possible,
  so the MPAS documentation and source can be used as a cross-check).

---

## 1. Why

The A-grid core (all variables at cell centers) has been made stable, but
every mechanism that stabilizes it also suppresses the weather it is
supposed to produce. Measured on the current code:

| mechanism | why it exists on the A-grid | cost to eddies |
|---|---|---|
| adjoint (transposed-gradient) divergence | energy-neutral pairing kills the 2Δx checkerboard null mode | rough pointwise output |
| divergence-field Jacobi smoothing | hides the adjoint operator's roughness | removes ~10% of eddy-scale divergence per step |
| divergence damping (ν=1e7) | residual checkerboard, gravity-wave noise | ~11 h e-fold on eddy secondary circulations |
| strong ∇⁴ hyperdiffusion (1–3 h at 2Δx) | operator-generated grid noise | measured to be the dominant eddy suppressor: EKE doubling 18 d with it, 1.5 d without it (= the Eady-predicted rate) — but without it the grid dies in 10 days |

The root cause is structural: with pressure and velocity collocated, the
pressure gradient cannot see a 2Δx pressure checkerboard, so the
gravity-wave subsystem has a null mode that must be damped from outside.
On a C-grid the normal velocity lives on cell edges; the pressure gradient
across an edge is a two-point difference between the adjacent cells, and
the divergence of a cell is the sum of fluxes through its own edges. A
checkerboard produces the *maximum* pressure force on every edge; the
null mode does not exist. The stabilization stack above is deleted and
only a closure-strength ∇⁴ (the sub-grid turbulence closure every model
carries) remains.

Everything vertical and everything physical is untouched (Section 8).

### Non-goals

Moisture, non-hydrostatic dynamics, topography (the C-grid makes it
easier later, but the aqua-planet stays flat for now), higher-order
transport schemes (2nd-order centered first; upgrade path noted), and a
new time integrator (AB4 stays; RK3 is an optional later swap).

---

## 2. Mesh

### 2.1 What `grid.module.js` already provides

- Cell centers `centerVertex` (unit vectors) from the Snyder equal-area
  (ISEA) projection; `10N²+2` cells with `N` cells per icosahedron edge.
- `cell.index`: position in `Grid` iteration order (north pole 0, quads
  0..9 row-major, south pole last); `grid.size = 10N²+2`.
- `cell.neighbors`: the 5 or 6 neighbors in counter-clockwise order (viewed
  from outside), every seam case handled, computed once at construction;
  symmetric (A lists B ⇔ B lists A) and duplicate-free.
- `cell.vertices[k]`: the unit-sphere circumcenter of `(center,
  neighbors[k], neighbors[k+1])` — the Voronoi vertex, in CCW order — as a
  `GridVertex` object shared by the three cells around it; `grid.vertices`
  lists all `2C−4` of them by `vertex.index`. Because they are true
  circumcenters of the Delaunay triangles, the mesh is exactly orthogonal:
  every primal edge is perpendicular to its dual edge to roundoff. This is
  the property TRiSK requires.
- `cell.area`: exact spherical polygon area (sums to 4π to roundoff);
  `cell.isPentagon`, `cell.isPole`.
- `new Grid(N)` applies 10 Lloyd iterations by default (Section 8 item
  2); `new Grid(N, { relax: 0 })` is the raw ISEA tessellation.
- `test/grid.test.mjs` certifies all of the above at N = 2, 3, 5, 8, 16
  (`node --test`).

Measured at N=32 (10,242 cells, ≈240 km spacing): 30,720 edges (= 3C−6),
12 pentagons, primal/dual edge-length ratio `l_e/d_e` in [0.40, 0.89]
(regular hexagon: 0.577), max/min cell area 1.24. The distortion is
concentrated around the 12 pentagons and is intrinsic to icosahedral
meshes; it degrades operator accuracy locally but not the scheme's
conservation properties.

### 2.2 The three staggered locations

```
        primal cell i (hexagon/pentagon): scalars   — mass, θ, π, Φ, K, divergence
        primal edge e (between cells i,j): u_e      — normal velocity component
        dual vertex v (triangle center):   ζ_v, q_v — vorticity, PV
```

Each edge `e` separates cells `i,j` and joins vertices `v1,v2`. Two
lengths: `d_e` = arc distance between the cell centers (the dual edge),
`l_e` = arc distance between the vertices (the primal edge). The rhombus
`(i, v1, j, v2)` has area `½ d_e l_e`; the primal edge splits it into two
triangles of area `¼ d_e l_e`, one in each cell, and the dual edge splits
it into two triangles of area `¼ d_e l_e`, one in each dual triangle.
These are planar identities: on the sphere with arc lengths the kite sum
over all cells falls short of 4π by 4e-5 at N=32 (1e-5 at N=64), so they
hold only to O(Δx²). The scheme needs an exact partition of area, not
planar areas, so all areas are spherical:

```
A_i     = cell.area from Grid                       spherical polygon
A_v     = spherical excess of triangle (i, j, k)    dual-triangle area
R_{i,v} = area of spherical quadrilateral (x_i, m_e, x_v, m_e')   kite,
          m_e = normalize(x_i + x_j) the point where edge e crosses its dual edge
```

With these, `Σ_v R_{i,v} = A_i` and `Σ_i R_{i,v} = A_v` hold to roundoff
(measured 6e-15 at N=8, 2e-12 at N=128). The planar `¼ d_e l_e` appears
only as the edge weight in the operators of Section 3, where the
identities that matter (3.2, T1, T3) need only internal consistency.

`EC(i)`: edges of cell i. `EV(v)`: the 3 edges meeting at vertex v.
`ECP(e)`: the edges of the two cells sharing e, excluding e (10 for two
hexagons).

### 2.3 Orientation conventions (fixed; every operator below uses them)

- `n_e`: unit normal to the primal edge at `m_e = normalize(x_i + x_j)`, the
  point where the primal and dual edges cross (the primal edge's own
  midpoint differs from it by up to 0.2 l_e near pentagons), pointing from
  `cellsOnEdge[e][0]` to `cellsOnEdge[e][1]`, tangent to the sphere.
- `t_e = k_e × n_e`, where `k_e = m_e` is the local vertical. `t_e` is `n_e` rotated 90° counter-clockwise when
  viewed from outside the sphere.
- `verticesOnEdge[e] = [v1, v2]` ordered so that `(x_v2 − x_v1) · t_e > 0`.
- `n_{e,i} = +1` if `n_e` points out of cell i (i.e. `i = cellsOnEdge[e][0]`), else −1.
- `t_{e,v} = +1` if traversing the dual edge in the `n_e` direction circulates
  counter-clockwise around vertex v, else −1.
- `u_e = u · n_e` is the prognostic normal velocity; `u⊥_e = u · t_e` is
  the diagnosed tangential velocity.

### 2.4 Data layout (structure-of-arrays, typed arrays)

Cells (C = 10N²+2):

| array | length | contents |
|---|---|---|
| `xCell` | 3C | unit position vectors |
| `latCell`, `lonCell` | C | for Coriolis at cells (init, radiation, diagnostics) |
| `areaCell` | C | `A_i` on the unit sphere × `a²` |
| `nEdgesOnCell` | C | 5 or 6 |
| `edgesOnCell`, `cellsOnCell`, `verticesOnCell` | 6C (padded with −1) | CCW order; `verticesOnCell[i][k]` is between `edgesOnCell[i][k]` and `[k+1]` |
| `edgeSignOnCell` | 6C | `n_{e,i}` |
| `kiteAreasOnCell` | 6C | `R_{i,v}` aligned with `verticesOnCell` |

Edges (E = 3C−6):

| array | length | contents |
|---|---|---|
| `cellsOnEdge`, `verticesOnEdge` | 2E | oriented per 2.3 |
| `dcEdge`, `dvEdge` | E | `d_e`, `l_e` (meters) |
| `xEdge`, `nEdge`, `tEdge` | 3E | midpoint, normal, tangent |
| `nEdgesOnEdge` | E | size of `ECP(e)` (≤ 10) |
| `edgesOnEdge`, `weightsOnEdge` | 10E (padded) | `ECP(e)` and TRiSK weights `w_{e,e'}` |
| `fEdge` | E | Coriolis parameter at edge midpoint |

Vertices (V = 2C−4):

| array | length | contents |
|---|---|---|
| `xVertex` | 3V | circumcenters (deduplicated across the 3 cells that share each) |
| `areaTriangle` | V | `A_v` |
| `cellsOnVertex`, `edgesOnVertex` | 3V | CCW |
| `edgeSignOnVertex` | 3V | `t_{e,v}` |
| `kiteAreasOnVertex` | 3V | `R_{i,v}` |
| `fVertex` | V | Coriolis parameter at vertex |

Prognostic state, per layer k (K=20 layers, unchanged):
`u[k*E+e]`, `theta[k*C+i]`; per column: `pi[i]`, `surfaceT[i]`.
Diagnostic per layer: `massFlux[E]`, `uPerp[E]`, `divergence[C]`,
`vorticity[V]`, `ke[C]`, plus the column diagnostics already present
(Exner, geopotential, `pi*dSigma/dt`).

Sizes at N=32, K=20: `u` 614k doubles, `theta` 205k. Fine.

### 2.5 Mesh builder (`js/mesh.module.js`, pure geometry, testable standalone)

1. Cell indices are `cell.index`; neighbors are `cell.neighbors`.
2. Enumerate edges: for cell i and neighbor index k, the edge to
   `cellsOnCell[i][k]` is created once (when `i < j`) and lies between
   `verticesOnCell[i][k−1]` and `[k]`.
3. Vertices are `grid.vertices`, already deduplicated and indexed;
   `cellsOnVertex[v]` is `(i, nb[k], nb[k+1])` for any cell `i` whose
   `vertices[k]` is `v` (CCW by construction).
4. Compute lengths, areas, kites, normals/tangents, orientation signs.
5. Compute TRiSK weights (Section 3.5).
6. Compute `f` at vertices, edges, cells from latitude
   (`f = 2Ω sin φ`, sidereal Ω).
7. Emit typed arrays; the worker consumes only these (no `Grid` objects).

Geometry tests (M0 acceptance):

- Counts: `E = 3C−6`, `V = 2C−4`; every edge has 2 cells and 2 vertices;
  every vertex has exactly 3 cells and 3 edges.
- `Σ_v R_{i,v} = A_i` and `Σ_i R_{i,v} = A_v` to 1e-12 relative (spherical
  kites); `Σ_e ¼ d_e l_e` matches `A_i` to O(Δx²) and no better.
- Orthogonality: `|n_e · (x_v2 − x_v1)| / l_e < 1e-10` for every edge.
- Orientation: `edgeSignOnCell` sums to zero flux for the constant vector
  field (discrete divergence of a uniform tangent field is O(Δx²), not
  O(1)); `t_{e,v}` gives positive circulation for a solid-body rotation.
- Total area `4πa²` to 1e-12 relative.

---

## 3. Discrete operators

All per layer; `i,j` cells, `e` edges, `v` vertices.

### 3.1 Divergence (edges → cells)

```
D_i = (1/A_i) Σ_{e ∈ EC(i)} n_{e,i} F_e l_e
```

with `F_e` any edge flux (e.g. `m_e u_e`). Mass-conserving by
construction: each edge's flux enters one cell and leaves the other with
opposite sign, so `Σ_i A_i D_i = 0` to roundoff for any `F`.

Test: uniform tangent field → `D` is O(Δx²) and its global area integral
is < 1e-12 relative; solid-body rotation → `D` ≈ 0 pointwise.

### 3.2 Gradient (cells → edge normal component)

```
(∇φ)_e = (φ_j − φ_i) / d_e        j = cellsOnEdge[e][1], i = cellsOnEdge[e][0]
```

This is the *only* horizontal derivative the momentum equation needs, and
it is the one that sees the checkerboard. `Σ_e d_e l_e (∇φ)_e F_e =
−Σ_i A_i φ_i D_i(F)` holds exactly (discrete integration by parts; RTSK
Appendix A.1 — the weight is the full `d_e l_e`, twice the rhombus area
`½ d_e l_e` that an edge owns geometrically). This is the
energy-consistency property the A-grid needed the adjoint construction to
fake.

Test: the identity above for random `φ`, `F`, to 1e-12 relative
(`test/mesh.test.mjs`, passes at 1e-15).

### 3.3 Curl (edges → vertices)

```
ζ_v = (1/A_v) Σ_{e ∈ EV(v)} t_{e,v} u_e d_e
```

(Circulation around the dual triangle: `u_e` is the component along the
dual edge of length `d_e`.) Absolute vorticity `η_v = ζ_v + f_v`.

Test: solid-body rotation with angular velocity `ω` about the polar axis
gives `ζ_v = 2ω sin φ_v` to O(Δx²); `Σ_v A_v ζ_v = 0` for any `u`.

### 3.4 Kinetic energy (edges → cells)

```
K_i = (1/A_i) Σ_{e ∈ EC(i)} ¼ d_e l_e u_e²
```

Test: uniform field `U`: `K_i = ½|U|²` exactly on a regular hexagon, to
O(Δx) elsewhere.

### 3.5 Tangential velocity reconstruction (TRiSK)

```
u⊥_e = (1/d_e) Σ_{e' ∈ ECP(e)} w_{e,e'} l_{e'} u_{e'}

w_{e,e'} = n_{e,i} n_{e',i} ( ½ − Σ_{v ∈ V(e→e', i)} R_{i,v} / A_i )
```

where `i` is the cell containing both `e` and `e'` and `V(e→e', i)` is
the set of vertices of cell `i` strictly between `e` and `e'` walking
counter-clockwise around `i`. This is Thuburn et al. 2009 eq. 33–34 with
the energy-conserving split constant `a = ½`, and it follows from one
physical statement: the mass leaving each kite `R_{i,v}` through its two
primal half-edges (half of each edge's flux) and its two dual half-edges
balances the cell divergence apportioned by kite area. Implemented in
`mesh.module.js` and checked three independent ways: the property tests
below; a literal port of MPAS's `buildEdgesOnEdgeArrays` (which folds
`l_{e'}/d_e` into `weightsOnEdge`) agrees bit for bit at N=8; and the
closed form was re-derived from both papers.

- **T1 (energy):** `w_{e,e'} = −w_{e',e}` for every pair, so the Coriolis
  term does no work: `Σ_e d_e l_e u_e u⊥_e = 0` for any `u`. Measured
  1e-13 relative.
- **T2 (consistency):** the reconstruction is exact for a uniform field
  only when the dual edge bisects the primal edge. On the raw ISEA grid
  the crossing point `m_e` sits up to 0.19 `l_e` from the primal midpoint —
  a property of the projection, not of the pentagons — and the worst-edge
  error for solid-body rotation is 12% at every N (mean 1.5% at N=32),
  tracking that offset edge for edge. Eight Lloyd iterations (cell centers
  moved to their spherical-polygon centroids, vertices recomputed) reduce
  it to 0.9% worst / 0.2% mean with the area ratio unchanged; see Section
  8 item 2.
- **T3 (mass consistency on the dual mesh):** for any `u` and every vertex
  `v`, `Σ_{e ∈ EV(v)} t_{e,v} u⊥_e d_e = −Σ_{i ∈ CV(v)} R_{i,v} D_i`: the
  flux of `u` out of each dual triangle equals its kite-weighted
  divergence (RTSK eq. 25 with the dual divergence of eq. 28; T09 eq. 12
  states it with the opposite circulation sign). It is what makes the
  dual-cell mass — and hence PV — evolve consistently and keeps
  geostrophic modes stationary. Measured 1e-15.

### 3.6 Laplacian of velocity (for the ∇⁴ closure)

```
(∇²u)_e = (D_j − D_i)/d_e  −  (ζ_v2 − ζ_v1)/l_e
(∇⁴u)_e = ∇²(∇²u)_e
```

(vector Laplacian = grad div − curl curl, evaluated with the operators
above). Scale-selective by construction — no first-derivative
contamination, unlike the hex ring-average Laplacian on the A-grid — so
the (2Δx/L)⁴ selectivity estimate should actually hold. Test: on a
uniform field both terms vanish; on a checkerboard `u` the result has the
expected sign and magnitude.

### 3.7 Cell-center velocity reconstruction (edges → cell vectors)

```
v_i = (1/A_i) Σ_{e ∈ EC(i)} ½ d_e l_e u_e n_e
```

Exact for a uniform field on a regular hexagon (`Σ_e ½ d l (U·n_e) n_e =
½ d l · 3U = A_i U`). Used only by physics that need a speed (surface
drag, sensible heat flux), by the top sponge's zonal mean and the
gravity-wave drag (M21 item 1, The model top), and by the viewer.
Pole cells need no special basis: `v_i` is a 3-vector in the tangent
plane; zonal/meridional components are computed for display only.

---

## 4. Governing equations on the C-grid

Vertical coordinate `σ = p/π`, K layers, unchanged. Per layer `k` the
layer mass per area is `m_i = π_i Δσ_k / g`; edge mass `m_e = ½(m_i + m_j)`
(centered; an upwinded variant is a later option). `Δσ_k` is fixed so
`m_e u_e Δ`-bookkeeping reduces to `π_e u_e Δσ_k`.

### 4.1 Continuity and vertical velocity (Steps 1–3, structure preserved)

```
F_{e,k}      = π_e u_{e,k}                                    (per unit Δσ)
D_{i,k}      = (1/A_i) Σ_e n_{e,i} F_{e,k} l_e
∂π_i/∂t      = −Σ_k D_{i,k} Δσ_k
(π σ̇)_{i,k+½} = −Σ_{k'≤k} D_{i,k'} Δσ_{k'}  −  σ_{k+½} ∂π_i/∂t
```

The last line is exactly the current `calculate_lower_pi_dSigma_dt`; it
telescopes to zero at the ground because `∂π/∂t` is exactly the column
sum (no smoothing is folded in — the π smoothing filter is deleted along
with the reason for it).

### 4.2 Momentum (vector-invariant form, on edges)

```
∂u_e/∂t = + η_e u⊥_e
          − (K_j − K_i)/d_e
          − (Φ_j − Φ_i)/d_e  −  c_p θ_e (∂Π/∂π)_e (π_j − π_i)/d_e
          − [ (π σ̇ u)_{e,k+½} − (π σ̇ u)_{e,k−½} − u_{e,k} ((π σ̇)_{e,k+½} − (π σ̇)_{e,k−½}) ] / (π_e Δσ_k)
          + F_e
```

- `η_e = ½(η_v1 + η_v2)`: absolute vorticity averaged to the edge. RTSK's
  fully PV-conserving variant uses the PV-weighted mass flux `q_e F⊥_e`
  with `q_v = η_v/m_v`, `m_v = (1/A_v) Σ_i R_{i,v} m_i`; adopt that form
  in M1 if the simpler one shows PV drift in the Rossby–Haurwitz test.
  With variable `q` energy conservation needs the PV averaged *inside* the
  reconstruction sum, `Q⊥_e = (1/d_e) Σ_{e'} w_{e,e'} l_{e'} F_{e'}
  ½(q̃_e + q̃_{e'})` with `q̃_e = ½(q_v1 + q_v2)` (RTSK eq. 49–50), not
  multiplied onto `F⊥_e` afterwards.
- **Sign of the Coriolis/vorticity term**, derived from the conventions in
  2.3: `−η (k×u) · n_e = +η u⊥_e`. Sanity check: NH, `Φ` decreasing
  poleward, edge with `n_e` northward ⇒ `t_e` westward, steady state gives
  `u = −(∂Φ/∂y)/f > 0`, westerly. The geostrophic-init sign error that
  produced the "very wobbly" start on the A-grid came from deriving this
  by hand from the textbook convention instead of from the code's; here
  the convention is fixed above and **Williamson TC2 (Section 6) is the
  test that certifies the sign** — a wrong sign fails it within hours.
- Pressure-gradient force: the same two-term σ-coordinate PGF as now
  (`−∇Φ|_σ − c_p θ ∂Π/∂π ∇π`), with both gradients now honest two-point
  edge differences; `θ_e`, `(∂Π/∂π)_e` are cell averages to the edge.
- Vertical advection: the current energy-consistent flux-difference form,
  with `(π σ̇)_e` averaged from the two cells and `u_{e,k±½}` averaged
  between layers.
- `F_e`: bulk surface drag on the lowest layer using `|v|_e = ½(|v_i| +
  |v_j|)` from 3.7, the ∇⁴ closure `−K₄ (∇⁴u)_e`, and at the top the
  sponge on the zonally asymmetric wind and the gravity-wave drag (M21
  item 1, The model top).

### 4.3 Thermodynamics (cells)

Flux form, consistent with the mass flux used in continuity so that
`∫θ dm` is conserved exactly by horizontal transport:

```
∂(π θ)_i/∂t Δσ_k = −(1/A_i) Σ_e n_{e,i} F_{e,k} θ_e l_e Δσ_k
                   − [ (π σ̇ θ)_{i,k+½} − (π σ̇ θ)_{i,k−½} ]
                   + π_i Δσ_k Q_i / (c_p Π_i)
```

with `θ_e = ½(θ_i + θ_j)` (second-order centered) initially. Upgrade
path: 3rd/4th-order upwind-biased flux (Skamarock & Gassmann 2011) when
the centered scheme's dispersive ripples at fronts become the limiting
noise source. The interface `θ_{k±½}` interpolation and the radiation
heating `Q_i` are unchanged from the current code. The thermal ∇²
diffusion on θ is dropped in favour of a ∇⁴ on θ (`nu4Theta`), which
the thin top layers require (Section 6, M2).

### 4.4 Time integration

M1 uses classical RK4 for clean convergence and conservation studies.
For the full model, Adams–Bashforth (order ramp 1→4) as now; the state to step is the flat
`u`, `theta`, `pi`, `surfaceT` arrays, so the stepper history becomes a
few large typed arrays instead of an object per variable. `dt` from the
gravity-wave CFL using the *minimum* `d_e` (pentagon neighborhoods have
the shortest dual edges), with the same ~2× margin found empirically on
the A-grid (`c·dt/d_min ≲ 0.3–0.4` with `c ≈ 300 m/s`). At N=32 expect
`dt ≈ 300 s` again. RK3 is the optional later swap if the AB4
imaginary-axis limit bites at finer N.

---

## 5. What is deleted, what is kept

Deleted from `sim.js` (A-grid life support, never ported):

- `GradientHelper` (lerp stencils, 2-neighbor selections), `ColumnNeighbor`
  / `LayerNeighbor` and all velocity rotation between local bases.
- Adjoint divergence: `extractStencilWeights`, `divIncoming`,
  `divergenceOfPiV`/`divergenceOfV`.
- Divergence-field Jacobi smoothing (`smooth_divergence_layers`,
  `DIVERGENCE_SMOOTHING_*`).
- Divergence damping (`gradDivergence`, `DIVERGENCE_DAMPING_COEFFICIENT`).
  Keep a *small* optional `ν ∇_e D` term available: hexagonal C-grids do
  carry an extra branch of velocity degrees of freedom (E = 3C vs. the 2C
  a vector field needs) and MPAS retains weak divergence damping for it.
- ∇² eddy viscosity and the hex ring-average `laplacian()`; the π
  smoothing filter (`pi_smoothing`); the surface-pressure diffusion.
- `RayleighSponge` (zonal-mean relaxation): replace with the existing
  top-of-model Rayleigh friction toward zero (`σ < 0.05`), which is the
  wave absorber the sponge was meant to be. The ablation fleet showed the
  sponge neither helps stability nor limits eddies at current settings.
- Per-column `Layer`/`Column` object graph for the horizontal state
  (`layer.v`, `lap2v`, `divPiV`, …) → flat arrays.

Kept, ported to operate on flat arrays (no algorithmic change):

- σ levels, `calculateSigma`, `initialThetas`, the Exner steps 4–6, the
  interface-θ hydrostatic integration (step 8), the σ̇ telescoping (steps
  1–3, now fed by the edge divergence).
- Gray longwave radiation column (`RadiationColumn`), solar geometry, slab
  ocean heat capacity, surface albedo.
- Bulk surface drag and sensible heat flux (speed from 3.7), gustiness
  floor, dry convective adjustment, PBL and top Rayleigh drag.
- Initialization: latitude-dependent surface temperature and σ-tapered θ
  shift, balanced surface pressure (bisection on a level 500 hPa
  surface), the wavenumber-5 eddy seed, the ±30°/±60° pressure bands.
  Geostrophic wind initialization changes representation: compute the
  cell-center geostrophic vector `v_i` from the cell-center PGF, then
  project to edges, `u_e = ½(v_i + v_j) · n_e`; TC2 certifies balance.
- Positivity floors, sidereal day, axial tilt, dt scaling.
- The synoptic-chart viewer path: cell-centered scalars (π, θ, T, z500)
  are unchanged; wind vectors come from 3.7; the contour layer's
  triangulation is the dual mesh, which the mesh builder now provides
  directly.

---

## 6. Validation protocol and milestones

Each milestone has acceptance tests that gate the next. All tests are Node
scripts under `test/` (the same harness style used throughout the A-grid
work; `three.module.js` imports cleanly in Node).

### M0 — Mesh (`js/mesh.module.js`) — done

`buildMesh(grid, {radius, omega})` emits the arrays of Section 2.4;
`js/dynamics/operators.module.js` implements 3.1–3.5 and 3.7.
`test/mesh.test.mjs` (`node --test`) certifies at N = 4, 8, 16: the
counts and connectivity of 2.5; kites partitioning every cell and dual
triangle and summing to `4πa²` (1e-12); edge frames orthogonal and
oriented as documented; divergence conserving mass and the gradient its
negative adjoint (1e-12); curl of solid-body rotation converging to
`2ω sin φ` and integrating to zero; T1 and T3 at roundoff; T2 and the
cell reconstructions converging in the mean. The velocity Laplacian (3.6)
is built and tested in M1 with the closure.

### M1 — Shallow-water core (`js/dynamics/shallowWater.module.js`) — done

The vector-invariant equations of Section 4.2 with `h` in place of `π`,
the energy-conserving PV flux of RTSK eq. 49, the ∇⁴ closure of 3.6, and
RK4 time stepping (`integrators.module.js`; AB4 is the M2 swap when
cost matters). Shared initial-condition and norm helpers live in
`test/helpers/sphere.mjs`; every test builds `new Grid(N)` (relaxed).

- **Williamson TC2** (`test/sw_tc2.test.mjs`, 5 days): `h` l₂ error
  2.2e-4 at N=16, 9.7e-5 at N=32; `u` l₂ 5.1e-3 and 1.3e-3; mass to
  1e-15, energy to 1e-10. This certified the sign conventions: the
  momentum tendency is `+Q⊥_e − (Φ_j − Φ_i)/d_e` with `u⊥_e = u·t_e`,
  `t_e = m_e × n_e`. On the raw ISEA grid the N=32 error is 4.8e-4 — 5×
  worse — which settled the relaxation default.
- **Williamson TC6** (`test/sw_tc6.test.mjs`, 14 days, N=32): mass
  exact, energy drift 1.5e-9, potential enstrophy drift 7.6e-4 (bounded;
  the energy-conserving PV average does not also conserve enstrophy),
  wavenumber-4 phase speed 11.3°/day against the nondivergent analytic
  12.2°/day, the usual shallow-water lag. No breakup.
- **Galewsky et al. 2004** (`test/sw_galewsky.test.mjs`, N=32, ∇⁴ at a
  3 h timescale on the 2Δx mode): the unperturbed balanced jet holds at
  2–5e-4 for three days and then goes unstable on its own, growing ×3.6
  per day (e-folding ≈ 0.8 d) from grid-scale truncation noise — the
  behaviour the paper describes for under-resolved grids, so the test
  asserts the pre-onset window. At N=64 the pre-onset error is 4.4×
  smaller (4.8e-5, second order) and the same growth starts about half a
  day later; the seed shrinks with resolution, the growth rate is the
  jet's own. With the 120 m perturbation the eddy
  kinetic energy grows from 2e-4 to 0.9 of the total over six days,
  e-folding ≈ 0.6 d through days 2–4, and the jet has rolled up by day
  6. Mass exact.

All three run in about a minute at N=32; `SW_TEST_N` selects other
resolutions.

### M2 — Multi-layer σ-coordinate dynamics (`js/dynamics/sigmaCore.module.js`) — done

`createSigmaCore(mesh, {levels, g, cp, R, p0, nu4, nu4Theta, forcing,
surfaceGeopotential})` carries the A-grid column onto flat arrays: state
`pi[C]`, `theta[K·C]`, `u[K·E]`; per step the layer mass fluxes and
divergences, `dπ/dt`, `πσ̇` telescoping to zero at the ground, Exner
ratios from the exact layer integral of `σ^κ`, interface θ interpolated
in Exner, geopotential integrated upward with each layer's θ over its
own Exner span, flux-form θ transport, and the vector-invariant momentum
equation of 4.2 with the RTSK PV flux, the two-term PGF, and the
flux-difference vertical advection. `createHeldSuarez` supplies the
Held & Suarez 1994 forcing through the `forcing` hook. RK4 over the
three arrays (`createRK4Arrays`).

- **Rest state** (`test/sigma_rest.test.mjs`): uniform `π` with
  horizontally uniform θ(σ) gives identically zero tendencies and a
  bit-identical state after 50 steps; an isentropic column reproduces
  the analytic hydrostatic geopotential and the exact layer-mean Exner
  to 1e-12.
- **Jablonowski & Williamson 2006** (`test/sigma_jw06.test.mjs`, N=16,
  θ initialized so the model's own hydrostatic integration reproduces
  the analytic geopotential, `π` uniform, the JW06 surface
  geopotential): the balanced base state holds for five days with
  surface pressure within 999.5–1000.8 hPa and winds within 0.9 m/s of
  the initial field, mass to 1e-15. The 1 m/s perturbation grows into
  the baroclinic wave on schedule — surface pressure 961/1023 hPa at
  day 10 at this 480 km resolution, deepening fastest through days
  7–10, against the paper's ~940 hPa at high resolution.
- **Levels.** The model's default is the CAM 26-level grid the paper was
  run on (`sigmaInterfaces()`: HOMME's `cami-26.ascii` read as σ with
  `p_s = p0`, plus the cap above CAM's 2.19 hPa lid, 27 layers). The
  A-grid's 20 levels — eleven of them above 135 hPa and a single 800 m
  boundary layer — are gone. A 22-level set with a 10 hPa top and five
  30 hPa boundary layers was tried on the way: the sets agree on the
  JW06 wave to 3 hPa at day 10, and that one was stable without the θ
  closure, but the CAM grid is the standard and is what stays.
  `sigmaInterfaces('bl34')` is the same grid above its 2.4 km
  interface (σ 0.7444) with the five layers below it replaced by
  twelve: interfaces at 40, 83, 133, 192, 265, 359, 482,
  649, 882, 1213 and 1698 m in the standard atmosphere, each layer up
  to 1.5 times as thick as the one below, 34 layers in all. A saved
  state records its interfaces (`levels`), and every loader builds the
  model on them; a state without them is on cam26. `remapLevels`
  (`js/physics/regrid.module.js`) carries a state onto another grid
  conservatively in σ, each layer the σ-weighted mean of the layers it
  overlaps (θ as its temperature, keeping enthalpy, and stratified where
  a layer is split: M21), so cam26 and bl34 exchange their layers above
  2.4 km unchanged; a spin-up seeded from a state on the other grid
  carries its atmosphere across with it. The page starts a run from
  nothing (`climate.html?from=none`) on bl36; a run started from a saved
  state is on that state's grid.
- **A required closure on θ.** Without any θ dissipation the A-grid
  set's thin top layers (Δσ ≈ 0.0008 above ~2 hPa) went unstable from
  day 3 in the steady run — θ departures of hundreds of kelvin with no surface
  signal — independently of the time step. Uniform layers avoid it,
  top-of-model Rayleigh drag makes it worse (it unbalances the jet),
  and a scale-selective ∇⁴ on θ at the same 3 h timescale as the
  momentum closure removes it without touching the troposphere. Both
  closures are on in the tests; Section 4.3's "only if fronts demand
  it" has been answered by the stratosphere instead.

The 200-day Held–Suarez climatology runs as an experiment, not a test.

### M3 — Physics hookup

Port radiation, surface fluxes, convective adjustment, drags, closure,
initialization onto the `forcing` hook. The A-grid model is not used as
a baseline — its vertical grid and horizontal operators differ too much
for a profile comparison to mean anything. Acceptance:

- Column-local pieces conserve exactly: the radiation column's layer and
  surface fluxes sum to absorbed solar minus outgoing longwave; the
  convective adjustment conserves column enthalpy and leaves a
  statically stable column; drag only removes kinetic energy.
- The balanced initialization rings at ≤ 1 hPa over the first days.
- 90-day stability with hyperdiffusion at closure strength (2Δx timescale
  ≥ 3 h) and **no** divergence smoothing/damping. **Met** at N=16 from
  the balanced initialization with the full physics: mass to roundoff,
  surface temperature steady at 287.5–287.9 K, absorbed solar 238 and
  OLR 238–240 W/m² after day 30, and the ∇⁴ closures as the only
  dissipation. Eddy kinetic energy at 250 hPa grew from 4 to 184 m²/s²
  (5-day doubling early; the A-grid reached 4 in 90 days at an 18-day
  doubling), the jets to 23 m/s in the winter hemisphere, and by day 90
  the winter hemisphere's zonal-mean surface pressure had a subtropical
  high (1025 hPa at 35°S) above a subpolar minimum (1016 hPa at 65°S)
  with surface westerlies beginning at 50–60°S — the structure the
  A-grid never produced. The maximum wind anywhere climbed 1.5 m/s per
  day to 136 m/s, presumably the stratospheric winter jet with nothing
  above the closure to bound a zonally symmetric flow; the M5 runs
  locate it and test a top-of-model drag.

Status: `js/physics/radiation.module.js` (two-band gray column, ozone
shortwave absorption, insolation with tilt, slab ocean, sensible heat
flux — see the radiation note below), `js/physics/surface.module.js`
(bulk drag, boundary-layer drag, convective adjustment) and
`js/physics/init.module.js` are ported and `js/model.module.js`
assembles them on the `forcing` hook with state `[π, θ, u, T_s]`. The
reference θ(σ) is not inherited from the old model: `equilibriumProfile`
integrates a single column of this model's own radiation and convective
adjustment over a 305 K surface to radiative–convective equilibrium
(converged to roundoff in 1200 days of column time, a few seconds of
compute; θ ≈ 299/315/371/532 K at 850/500/200/50 hPa, surface air 298 K).
`test/physics.test.mjs` certifies the column budgets to roundoff;
`test/init.test.mjs` shows the balanced state ringing at 1.3 hPa (N=8)
and 1.2 hPa (N=16) over three days with the physics off. The
subtropical/subpolar pressure bands of the A-grid initialization are off
by default: nothing balances them and they ring at 10 hPa, and the
Held–Suarez run produces the surface wind structure they imitated on its
own. Held–Suarez at N=16 for 200 days: 213 hPa jets 28–32 m/s at
±40–50°, surface westerlies +9 m/s at ±50°, tropical easterlies −7,
upper-level eddy kinetic energy ≈ 210–240 m²/s², surface pressure
975–1035 hPa, statistically steady from day 100.
- EKE doubling time ≤ 3 days (Eady prediction for the current base state:
  1.5 days; A-grid delivered 14–19 days).
- The emergence experiment: subpolar surface lows at ±60° appearing in the
  zonal-mean profile within ~60 simulated days at N=32.

**Radiation (Sept 19, 2026).** The column is a three-band gray scheme.
A window band carrying 25% of blackbody emission is transparent, so the
surface radiates it straight to space. A vapour band (55%) has the
optical depth of Frierson et al. (2006): τ₀(φ) = τ_e + (τ_p − τ_e)
sin²φ with τ_e = 5.3, τ_p = 1.325, distributed as τ₀(0.1 σ + 0.9 σ⁴) so
it sits near the surface like water vapour (a prescribed stand-in for
it — fixed in time, no feedback). A well-mixed-gas band (20%, the
15 µm band's share) has optical depth 5 spread uniformly per unit mass,
so thin high layers keep an emissivity they can cool with, as CO₂ lets
the stratosphere do. Shortwave: 3% of the beam is absorbed aloft,
spread with a Lacis–Hansen ozone column (25 km ± 5 km) through which
the absorbing part of the beam decays with optical depth 4, so the
heating peaks near the stratopause; the surface takes (1 − albedo) of
the rest. Every flux still closes exactly (`test/physics.test.mjs`).

Tuned on 500-day N=3 runs to a 288 K annual mean (τ_e 5 → 287.4 K,
6 → 289.5 K; tropics ~300 K, poles ~253 K, OLR ≈ 240 W/m² against 241
absorbed). The single-column equilibrium over a 300 K surface has a
~205 K tropopause near 70–100 hPa and a stratopause of 265–277 K at
1–5 hPa, where the original gray column was isothermal at 207 K above
200 hPa. Without the gas band, ozone heating concentrated at the top
drove the thinnest layers to 300–400 K, because a mass-proportional
emissivity gives them no way to radiate. The equilibrium profile is
built at the area-mean optical depth (sin²φ = 1/3) under the global-mean
beam S/4 and needs 1200 days of column time to converge.

Effect on the circulation (120-day N=16 spin-ups): the tropospheric
lapse rate went from 6.2 to 8.9 K/km in the tropics — nearly
dry-adiabatic, the correct answer for a dry model — and the strongest
tropospheric wind is a ~75 m/s subtropical jet core near 160 hPa.
The top-of-model sponge stays, in the role of gravity-wave drag: with
the gas band cooling the dark winter pole, the polar-night jet at
σ ≈ 0.001 grows without bound when undamped (250 m/s by day 100) and
the original σ < 0.005 sponge only slows it (134 m/s at day 120), while
a sponge over σ < 0.02 (the top four layers, above ~20 hPa) with a
5-day timescale holds it to 44/59/71/84 m/s at days 30/60/90/120 —
the range of the real polar-night jet — at no measurable tropospheric
cost (EKE 152 vs 157–161, jets identical). Radiation sets up the
pole-to-pole contrast that drives that jet; only a momentum sink
bounds it, which is why removing the drag was not an option. On the gas
radiation of M21 the drag on the mean wind warms the winter pole by
10-37 K; the sink is now a gravity-wave drag and the sponge acts on the
eddies alone (M21 item 1, The model top).

### M4 — Worker and viewer — done

`js/model.worker.js` builds the grid and the model in a module worker,
integrates it, and posts surface pressure, surface temperature and
surface wind speed as transferable arrays every few simulated hours.
`unifiedViewer.module.js` gains a `dynamicColors` mode that keeps the
vertex color array and exposes `updateColors(rgbPerCell)`;
`climate.html` colors the globe by the chosen field with the global
diagnostics in a readout. The synoptic contour layer of the A-grid
viewer is not ported.

### M5 — Dissipation diet and emergence — done

Runs from the balanced initialization with the full physics, tilt on
from the spring equinox, N=16 unless noted; `emergence.mjs` logs the
zonal means every 5 days and writes surface snapshots every 10 days for
`climate.html?snapshot=`. Any panel setting can be named in the query string (`?view=space`, `?overlay=wind&level=250&animate=arrows`, `?projection=map&isobars=on`, `?palette=`, `?panel=closed`), and `?lat=&lon=&zoom=&roll=` set the globe's orientation, so a link opens an exact view. The trough metric is the zonal-mean surface
pressure at 35° minus that at 65° in each hemisphere; it starts at
−34 hPa because the balanced initialization has polar highs.

- **Emergent subpolar lows.** In the winter (southern) hemisphere the
  trough metric crosses zero at day 60–80 and reaches +18 to +23 hPa by
  day 120–180, the zonal-mean surface wind at 50–60°S turns from −4.6 to
  +3 m/s westerly, the 250 hPa jet grows to 49–53 m/s, and eddy kinetic
  energy at 250 hPa saturates around 320–340 m²/s². The summer
  hemisphere stays quiet (jet 7–9 m/s, no trough), as its base state
  provides little baroclinicity. This is the structure the A-grid never
  reached.
- **Dissipation diet.** Relaxing the ∇⁴ closures from 3 h to 10 h at
  the 2Δx mode changes nothing that matters (EKE 326 vs 335, jet 49 vs
  51, trough 17 vs 20 hPa at day 180): the closure is not what limits
  the eddies here, unlike the A-grid where it was the dominant sink.
- **The cap layer.** With nothing above the closure to bound a zonally
  symmetric flow, the winter jet in the σ < 0.002 cap above CAM's lid
  grows about 1.3 m/s per day to 230 m/s by day 180. A Rayleigh drag
  above σ = 0.05 with a 10-day timescale at the top holds it near
  60 m/s and moves the model's maximum wind to the real subtropical jet
  near 50–100 hPa at 80–90 m/s, at a small cost to the troposphere
  (EKE 293 vs 318, trough 15 vs 19 hPa at day 150). Whether to make it
  the default, or to replace the cap with a rigid lid at CAM's top, is
  open.
- **Energy balance.** Absorbed solar 238 W/m² against OLR falling from
  253 to 240: the initialization is warmer than this gray atmosphere's
  equilibrium, and the slab ocean cools at 0.02 K/day, from 288.9 K at
  day 40 to 286.2 K at day 180, slowing as OLR approaches 238.

- **A full seasonal cycle** (base, 360 days, 97 minutes at N=16). The
  storm track follows the winter hemisphere: the southern trough peaks
  at +24 hPa around day 210 and decays to zero by day 300 as the sun
  moves north, while the northern one goes from −9 hPa at day 230
  through zero at day 280 to +24 hPa at day 360, with surface
  westerlies of +3 m/s at 50–60°N and a 43 m/s jet by the end. Eddy
  kinetic energy at 250 hPa cycles from 335 at the southern solstice
  through 125 near the equinox to 238 approaching the northern
  solstice. The energy budget closes from about day 250 (OLR 235–240
  against 238 absorbed) with the slab ocean settling near 286.5–287.2 K.

- **Top drag over a full year.** The 10-day drag above σ = 0.05 keeps
  the model's maximum wind at 60–100 m/s (in the real subtropical jet
  near 50–100 hPa rather than the cap) but costs the troposphere: eddy
  kinetic energy 7–10% lower through the year, jets 1–5 m/s weaker, and
  the winter surface trough weaker or later (11 vs 24 hPa at day 360,
  −2 vs +23 hPa at day 120). Part of that is the chaotic timing of a
  single realization, part is the ramp reaching the lower stratosphere.
  A drag confined to the cap layer (σ < 0.005) separates the two: over
  the same year its winter-mean eddy kinetic energy (306/196 vs
  308/199), jets (49.8/38.2 vs 50.7/37.4 m/s), troughs and westerlies
  are indistinguishable from the undamped run, while the maximum wind
  stays at 90–140 m/s in the upper CAM layers and follows the season
  instead of growing. That drag is now the model default
  (`createModel`: `topSigma: 0.005, topDragDays: 10`).

- **Resolution.** N=32 (240 km, 120 days, 4.6 hours) tracks N=16
  closely: eddy kinetic energy 235 vs 255, jet 40 vs 46 m/s, the
  winter trough +8.5 vs +22.6 hPa at day 120 with the same sign
  reversal near day 50–70 — the emergence is not a resolution artifact,
  and the timing of the surface trough is the chaotic part of a single
  realization at either resolution.

The A-grid project's stated goal — surface pressure cells maintained by
the model's own eddies rather than imposed — is met in the winter
hemisphere at both resolutions. Open next steps: ensembles or longer
runs to average the trough statistics; the gray radiation's climate
(the slab settles near 286.5 K, cooler than the 288 K it was aimed at);
and the summer-hemisphere weakness, which follows from the base state's
lack of baroclinicity there.

### M6 — Worker-thread time step (`js/parallel.module.js`) — done

The dependency chain inside one RK4 stage rules out pipelining across
workers (every phase needs the previous phase's output for the whole
sphere), so the parallel engine partitions each phase instead. The
σ-core tendency is split into phases that partition cleanly:

| phase  | partitioned over | work                                                  |
|--------|------------------|-------------------------------------------------------|
| flux   | layers           | edge mass fluxes and their divergence                  |
| column | cells, vertices  | dπ/dt, σ̇, Exner/geopotential column, vertex π, surface wind |
| layer  | layers           | θ tendency, momentum (PV flux, PGF, vertical advection), ∇⁴, drag |
| cell   | cells            | radiation and surface fluxes, slab ocean               |

Every array that crosses a phase boundary — the state, the four RK4
stages, the trial state and the core's intermediate arrays — lives once
in `SharedArrayBuffer`s (`createModel(..., { buffers })` adopts them;
`shareMesh`/`meshFromShared` pass the mesh the same way). The main
thread drives the phases through an `Atomics` generation counter, so a
step costs 16 phase barriers plus advance, combine and adjust. Within a
phase the workers claim work units from a shared chunk counter — a
layer's momentum tendency or its tracer transport, a block of 512–2048
cells, a block of vertices, a 64k-element slice of an array — so the
split adapts to the speed of each core (Apple silicon mixes performance
and efficiency cores) and no worker waits on a straggler. Every array
element is computed by exactly one worker with the single-thread
arithmetic, so the state is bit-identical to `createModel`'s
(`test/parallel.test.mjs` asserts element-wise equality of every state
array after four steps); only the radiation totals are summed in a
different order.

Measured on the 10-core laptop (8 performance + 2 efficiency cores),
moist model: N=64 2558 ms/step serial → 365 ms on 10 workers (7.0×),
against 425 ms (6.0×) with static per-worker partitions. The layer
phase, 60% of the step, scales to 6.5× and is most likely memory-bound.
The emergence driver takes `WORKERS=<n>` in the environment; the default
is every core.

**Step cost (Sept 19, 2026).** Three changes took the N=64 step from
359 to 253 ms on 10 workers: the ∇⁴ closures are applied once per step
after the RK4 dynamics instead of inside all four stages (with a 3 h
timescale the explicit step is far inside stability; JW06 unchanged);
radiation, surface fluxes and evaporation likewise run once per step as
an explicit increment of the state, which also makes the water budget
close to roundoff because evaporation is applied exactly once; and the
TRiSK weights carry the neighbour's `dvEdge` from mesh build. Fusing the
three tracer transports into one sweep gained nothing and was dropped.
The time step is now 900·16/N s (2× the original; RK4's gravity-wave
limit is ~3×, where a 3-day N=16 run differs by 0.1 hPa RMS and 0.36 m/s
from the original step), so N=64 runs a simulated day in 1.65 min of
wall time, 1.07 min at 3× (`DT_FACTOR=3` for the driver). The 3× step
(337.5 s at N=64, 1350 s at N=16) became the default for the page and
the drivers after M14: 100 days at N=16 from the tuned state track the
2× run within 0.2 K and the same ice cover, and 5 days at N=64 from
the page snapshot run clean at 1.24 min per simulated day, 19 simulated
hours per minute on 10 workers.

The same two files run in the browser: `threads.module.js` provides
the spawn/receive primitives from `worker_threads` or Web Workers, and
the page's model worker acts as coordinator — a dedicated worker may
block in `Atomics.wait`, and it spawns the phase workers as nested
workers. `SharedArrayBuffer` needs the page cross-origin isolated,
which `httpd.py` provides (COOP/COEP). `climate.html?workers=<n>`
chooses the count (cores minus two by default).

---

### M7 — Moisture — done (first tuning)

The moist gray-radiation aquaplanet of Frierson, Held & Zurita-Gotor
(2006), built on the dry model without changing its results when the
sources are off (`moist: false` carries q but never sources it).

- **Tracer.** `q` (specific humidity) is state[4], transported by the
  same flux-form scheme and mass fluxes as θ, with interface values
  interpolated in Exner like θ's. The ∇⁴ closure acts on the
  mass-weighted field π q, so the column's water is preserved exactly
  under transport (checked to roundoff in `test/moist.test.mjs`); θ's
  closure is unchanged. Virtual potential temperature θ(1 + 0.608 q)
  enters the hydrostatic integration and the pressure-gradient force.
  A filler removes negative q by borrowing from the layer below (any
  residue at the ground is counted as lost).
- **Evaporation.** The bulk formula of the sensible-heat flux applied
  to moisture, E = ρ C |v| (q_sat(T_s) − q_air), in the radiation
  column; the slab loses L·E and the lowest layer gains E.
- **Large-scale condensation.** Supersaturation is removed with one
  implicit step, warming the layer by L Δq / c_p and raining out at
  once; moist enthalpy c_p T + L q is conserved exactly.
- **Convection.** The simplified Betts–Miller scheme of Frierson
  (2007): a parcel from the lowest layer rises dry to its LCL and
  moist-adiabatically above; the column up to its level of zero
  buoyancy relaxes over 2 h toward that profile and a reference RH of
  60%. When the implied rain is positive the reference temperature is
  shifted so the enthalpy change equals L times the rain; otherwise
  both references are shifted so nothing is gained or lost. Dry
  convective adjustment follows, mixing q with θ.
- **Saturation** by Bolton's formula; the initial humidity is a
  relative humidity of 0.7 σ² times saturation.
- **Diagnostics.** Precipitation accumulates per cell (mm) and is
  reported as a rate; evaporation, latent and sensible heat, and total
  precipitable water join the budget line. The page overlays RH (at the
  selected height), precipitation and TPW.
- **Radiation coupling.** By default the vapour band's optical depth
  follows the model's own humidity: `vaporCoupling` m²/kg times each
  layer's water mass (0.55 with clouds in radiation; 0 restores the
  prescribed Frierson profile, which the single-column initial profile
  still uses). This is the water-vapour feedback the prescribed profile
  lacked.

Tuning (500-day N=3 runs, annual means). With the dry model's τ_e the
moist model settles ~6 K colder: evaporation cools the surface and
moist convection warms the upper troposphere, so the column radiates
more for the same surface temperature. Prescribed τ_e would have to
rise to ~12.5 for 288 K (contrast falling to 29 K at this resolution);
coupling reaches it at ~2.1 m²/kg with a larger contrast (35 K), and
the N=16 spin-up keeps 43 K at day 30 against the dry model's 52 K
(tropics cooled by evaporation, poles unchanged). Chosen: coupling 2.

Order of operations each step: RK4 dynamics with evaporation as a
tendency, then condensation, Betts–Miller, dry adjustment and the
filler as adjustments. One day at N=4 from the dry equilibrium closes
the water budget to 3% (the residual is the last-stage evaporation
estimate in the diagnostics, not a loss).

First moist climate (`moist2`, N=16, 120 days from the dry equilibrium
profile with 70% humidity, 8 workers, 14.5 min): latent heat 90–100
W/m², sensible 19–22, rain 3.1–3.4 mm/day, precipitable water rising
15 → 22 kg/m²; the tropical 850–500 hPa lapse rate is 6.3 K/km (dry
model 8.4) — the moist adiabat — and mid-latitudes 4.9. Eddies are far
stronger than in the dry model: EKE(250) 222/227/290/469 at days
30/60/90/120 against 63/89/117/157, the winter jet 56 m/s against 35,
and the winter subpolar trough appears by day 120. Global Ts 283–285 K
over these 120 days with the budget still +6 W/m² in the atmosphere's
favour, so the slab has not settled; the tropics are cool (292 K) and
the poles moist (18 kg/m²) for Earth, which the next tuning pass
should address (reference humidity, relaxation time, coupling).

### M8 — Cloud water — done

Cloud condensate `qc` is state[5], transported like `q` (conservative
closure, interface values in Exner) and loading the density
temperature, θ_v = θ(1 + 0.608 q − qc). The saturation adjustment
replaces instant rain-out: supersaturated vapour condenses into cloud
water and cloud water evaporates into subsaturated air, each with one
implicit step, so a layer is afterwards either saturated or cloud-free
and c_p T + L q is conserved exactly. Kessler autoconversion turns
cloud water above 0.2 g/kg into rain at 10⁻³ s⁻¹, and all cloud water
decays over 3 h; both are exact exponential decays per step. The rain
falls through the layers below within the step and evaporates into
each subsaturated one, up to `rainEvaporation` (default 1) of what
would saturate it with its latent cooling counted, so mid-level cloud
over dry air gives virga rather than surface rain; at N=128 this cuts
the Sahara's June large-scale rain from 1.44 to 0.13 mm/day, and over
90 days at N=64 it raises the global cloud water from 59 to 81 g/m²
and the planetary albedo from 0.279 to 0.303, for 0.2 % of the step's
GPU time. Betts–Miller still rains directly. Dry
adjustment mixes `qc` with `q`. The Betts–Miller reference profile uses
Bolton's closed-form LCL and a two-substep moist ascent that stops once
the parcel is 10 K colder than the air (87 → 5 ms per N=16 step).
Column cloud water (TCW) is a diagnostic and a page overlay; the step
costs 126 ms serial at N=16 against 98 ms dry. Cloud–radiation
coupling (albedo and longwave emissivity of cloud) is the next step.

### M9 — Sea ice and surface albedo — done

`js/physics/ice.module.js` is a zero-layer thermodynamic sea-ice model
in the manner of Semtner (1976) on the slab ocean, with ice thickness
as the seventh state array. Open water is the mixed layer (2.1×10⁷
J/m²/K); when it cools to the seawater freezing point (271.35 K) the
deficit freezes into ice of latent heat ρ_i L_f. Ice has a skin of
small heat capacity (2×10⁵ J/m²/K) whose temperature answers the
surface flux and the conduction k (T_f − T_skin)/h from the base (k = 2
W/m/K, h floored at 0.1 m for stability); the conducted heat freezes
water onto the base, a skin that would pass 273.15 K melts the ice from
the top instead, and ice that melts away returns its leftover energy to
the mixed layer. The surface energy — mixed-layer heat over the freezing
point, skin heat, minus the ice's latent heat — changes by exactly the
surface flux through every transition (`test/ice.test.mjs`).
`surfaceT` is the skin temperature the atmosphere sees in both states.
Albedo is 0.07 for open water, rising linearly to 0.6 at 0.5 m of ice.
Snow lies on the ice: precipitation that falls on an iced cell from air
below the melting point accumulates in water equivalent in the ocean
cells of the land surface's snow array (one field of snow on the
ground, land or ice, in the state file and on the page), brightens the
surface toward 0.75 over 20 kg/m², conducts in series with the ice
(0.31 W/m/K over its depth at 300 kg/m³, so a winter snow cover slows
the growth beneath it), melts before the ice, and goes into the water
when the ice is gone, its latent heat drawn from the mixed layer. Snow
that falls on open water melts at once and cools the water by its
latent heat, so cold water under snowfall freezes over; the ocean's
freshwater counts the precipitation when it falls. Snow heavier than
the ice's freeboard (109 kg/m² per metre of ice at 1026 and 917 kg/m³)
floods and freezes into snow-ice: the surplus leaves the snow and joins
the ice at ice density, conserving mass and energy, so a deep Antarctic
snow load thickens the pack instead of insulating it indefinitely.
The initial state carries 1.5 m of ice on the Arctic Ocean poleward
of 72°N and 0.7 m around Antarctica poleward of 68°S, over a zonal
climatological SST that reaches the freezing point near 70°. A prescribed ocean heat transport
(`oceanHeatFlux`, Q₀ = 20 W/m²) converges Q₀(3 sin²φ − 1) into the
mixed layer — zero in the global mean, cooling the tropics, warming
the poles by up to 2 Q₀ — and melts ice from below where it is
covered. It stands in for the poleward heat carried by ocean currents,
without which the ice–albedo feedback runs the ice edge to ~45°.

Ice covers a fraction A of its cell (`concentration`, a per-cell array
beside the snow, in the state file and on the page), with the thickness
h over that part, so the volume per cell area is A h; the rest, the
leads, is water held at the freezing point. The atmosphere sees the
area-weighted albedo and, in its fluxes, the area-weighted surface
temperature A T_skin + (1 − A) T_f. The cell's net flux is split
between the parts by the sunlight the leads absorb beyond the ice — the
direct and diffuse beams at the surface times the albedo contrast of
each — and by the longwave and turbulent heat the leads at the freezing
point lose beyond the colder ice, linearised as λ (T_f − T_skin) with
λ = 10 W/m²/K (`leadExchange`): the ice gains λ (1 − A)(T_f − T_skin)
per unit ice area and the water loses λ A (T_f − T_skin) per unit water
area, so the two still average back to the cell's flux. The ice part
grows and melts as above under its share; the leads' heat, with the
ocean's flux beneath them, melts ice or freezes new ice. Melting takes
area as Hibler (1979) does, half the relative loss of volume from the
area; ice frozen in the leads closes them as new ice 0.3 m thick
(`leadClosing`), and open water that cools below freezing forms ice
0.3 m thick over the area its volume covers. Ice below 1 % of its cell
or 10⁻⁴ m of volume melts away. The snow on lost area melts into the
water, its latent heat freezing the same mass onto the ice, and new
area joins at the freezing point without snow; the surface energy, A
times the ice's terms, still changes by exactly the surface and ocean
fluxes (`test/ice.test.mjs`). The dynamic ocean treats any cell with
ice as covered, holding its mixed layer at the freezing point. A state
saved without a concentration loads fully covered wherever it has ice.
The spin-up logs each hemisphere's extent, the area of the cells at
least 15 % covered, the measure behind Earth's 15 and 6 Mkm² (Arctic)
and 18 and 3 Mkm² (Antarctic) of ice at the seasonal maximum and
minimum.

### M10 — Clouds in radiation — done

Each layer's cloud water path gives it a gray emissivity
1 − exp(−130 m²/kg × path) that joins every longwave band: the vapour
and gas bands as 1 − (1 − ε_gas)(1 − ε_cloud), and the window, which is
now an exchange band of its own that is transparent only where there is
no cloud. In the shortwave the column's cloud optical depth
(`cloudScattering` = 95 m²/kg × path) reflects the beam
with the two-stream reflectance τ/(τ + 2μ); what passes is absorbed by
the surface with its per-cell albedo, with the multiple reflections
between surface and cloud summed. Cloud water also absorbs: a cloud of
water path W passes on only exp(−`cloudSolarAbsorption` × W) of what it
reflects or transmits, of the beam from above and of the light the
surface returns through it, and the rest heats its layers in proportion
to their water. The default 0.4 m²/kg absorbs 3.9 % at 100 g/m² and
15 % at 400 g/m² (Stephens 1978), and the scattering was re-balanced
from 55 to 95 m²/kg so that the planetary albedo is unchanged by the
absorption (two-day mean absorbed solar from `five64_day1734` within
0.3 W/m² of the purely scattering clouds'); `cloudSolarAbsorption` 0
with `cloudScattering` 55 gives back the purely scattering engine bit
for bit in both engines. The fixed planetary albedo of 0.3 is
gone: it is now produced by clouds and ice, and diagnosed. Radiation's
closure still holds exactly (with the latent heat of evaporation
counted as leaving the surface). The cloud–albedo, cloud–longwave and
ice–albedo feedbacks are all live from here.

Marine stratocumulus is diagnosed, not condensed. By default
(`mixedLayerDeck`) each ice-free sea column starts the bulk mixed-layer
model of `js/physics/mixedLayer.module.js` (Lilly 1968; Bretherton &
Wyant 1997; Stevens 2002) from its own state — h the deck's inversion
height carried from the last step (`prognosticHeight`), θ_l and q_t
the means of the layers below it, the first layer above h as the free
troposphere, the subsidence −πσ̇/(ρg) at h (πσ̇ averaged with equal
weights over the cell and its neighbours, twice over,
`subsidenceSmoothing`, at the two interfaces bracketing h), the
column's bulk surface fluxes and the DYCOMS-II longwave on the model's own liquid water —
advances it one physics step, keeps the new h, and takes its cover (1
when coupled, down to 0.3 as the buoyancy-integral ratio decouples it)
times 1 − A and its water path, at most `stratusWaterMax`, as the
deck's. The carried h never falls below the boundary-layer scheme's
Richardson depth, never rises above 3 km nor past the midpoint of the
layer above the lowest interface (above that depth) where θ_v jumps by
2 K — on the 26-layer grid the free troposphere's own stratification
across one of the thick layers above a kilometre passes the 2 K test,
and an unbounded deck deepens into it — and relaxes back to the
Richardson depth over a day (`mixedLayer.heightMemory`) while no deck
runs; unset (0) it starts from that depth. Where the deck runs, the
boundary layer's K-profile spans max(Richardson depth, h), so the
column is mixed through the deck's layer. The deck runs only in the
stratocumulus regime, and the regime test is the capping jump
Δθ_v ≥ 2 K at h. The subsidence only vetoes: the deck is refused
where the running-mean sink at h below is negative by more than
1 mm/s (`stratusSubsidence` −1 mm/s), large-scale ascent under which
an inversion at the boundary-layer top is transient. The model resolves
the deck regions' subsidence as a residual of a few mm/s (1.8 and
2.1 mm/s in the SE Pacific and Peru boxes, with a spread over the
cells as large), so a floor on the sink would turn columns of the regime
away on synoptic swings, while the inversion is the resolved record of
the subsidence that built it. One stage's πσ̇ carries the
divergent computational mode of the hexagonal C-grid at the
neighbouring-cell scale: on the day-183 N=128 state the sink at h
spreads over the SE Pacific's cells by 30 mm/s about a mean of
1.5 mm/s, 97 % of its variance at that scale, and the two ring passes
leave 8.6 mm/s about 1.6 mm/s with 12 % there. The subsidence test
reads a two-day running mean of the smoothed sink
(`subsidenceMemory`), whose grid-scale residual is about 0.5 mm/s, and
the deck follows a one-day
running mean G of the two tests' pass (`gateMemory`), running while
G > 0.5: a standing deck outlives failing tests by 17 h and a new one
waits as long, where the instantaneous Δθ_v test, hovering about 2 K,
would switch it in most columns from one step to the next. The three
running quantities are saved with the state as `mlmSubsidence`,
`mlmHeight` and `mlmGate`, and a state saved without them starts from
0, 0 (unset) and 0.5 (undecided); `prognosticHeight: false` with
`gateMemory: 0` and `subsidenceSmoothing: 0` gives back the deck
restarted each step from the Richardson depth behind the instantaneous
tests, bit for bit in both engines. By day the mixed layer is
heated by exactly the sunlight the column's radiation absorbs in the
deck's layer, per unit deck area (the overcast column's absorption
there less the clear column's; `stratusSolar`, false leaving the mixed
layer unlit while the column still absorbs), which weakens its
buoyancy flux, so a deck thins and decouples by day (RF01 run
standalone under a July sun at 30°N, absorbing 4 % per 100 g/m² by the
mixed-layer model's own formula: 0.56–0.59 of its night-time water by
mid-afternoon). Since θ_l and q_t restart from the column every step,
the deck feels one physics step of it, about 1.2 % less water after
900 s under a high sun. With
`mixedLayerDeck: false` the deck comes instead from an empirical fit:
over the ice-free part of each sea cell it covers the fraction f of
the column, whose predictor is the estimated inversion strength of
Wood & Bretherton (2006), EIS = LTS − Γ_θ (z_700 − z_LCL): the lower-tropospheric
stability LTS = θ(σ ≈ 0.7) − θ(lowest layer) less the rise in θ along
a moist adiabat from the lifting condensation level z_LCL to the
height z_700 of the σ ≈ 0.7 layer (both above the lowest layer), with
Γ_θ = g/c_p − Γ_m the moist adiabat's potential-temperature gradient
at 850 hPa and the two layers' mean temperature. f is
0.19 + 0.08 (EIS − 1), clamped to [0, 1] — 0.2 over the warm pool's
EIS of about 1 K, 0.67 over the south-east Pacific deck's 7 K, their
6–8 % per K — times a ramp from 0 over a 5 °C sea to 1 over 10 °C
(the fit is for subtropical and mid-latitude decks; polar fog is
another regime and would only add longwave warming over the ice),
times 1 − A for the ice concentration A; land carries none, nor does
the model with `moist: false`. With `stratusIndex: 'ectei'` the fit
takes instead the estimated cloud-top entrainment index of Kawai,
Koshiro & Webb (2017),
ECTEI = EIS − 0.23 (L/c_p)(q(lowest layer) − q(σ ≈ 0.7)), which
lowers the cover where the air the deck entrains is dry. The deck
fills the boundary layer above its condensation level: its base is
the lifting condensation level of the lowest layer's air (Bolton's,
as the convection finds it, reached along
the dry adiabat), its top the boundary-layer top the previous step
diagnosed, and for that thickness Δz it holds `stratusScale` × ½ Γ_l Δz²
of water, at most `stratusWaterMax` = 0.15 kg/m². Γ_l is the adiabatic
liquid-water lapse rate at cloud base, the density times the vapour a
saturated parcel condenses per metre of moist-adiabatic ascent:
2.4 × 10⁻⁶ kg/m³ per metre at 290 K and 950 hPa, 110 g/m² for a 300 m
deck at the adiabatic limit. The water sits in the layer nearest
σ = `stratusSigma` = 0.92, the top of a 1 km marine boundary layer,
where it scatters and emits as condensed water does but never enters
qc. The deck and the rest of the cell are independent columns: the
shortwave is computed with and without the deck's water and the two
weighted by f, and the deck layer's longwave emissivity is the
f-weighted mean of its emissivity with and without it. The page shows
the deck too: the `cloud` frame field, which the satellite view and the
cloud overlays draw, adds f times the deck's water to the column's
condensate, so a patchy deck shows as a fainter sheet. Without the deck
the eastern subtropical oceans, where Earth keeps its decks, were
nearly cloud-free on a year-two state (planetary albedo 0.11–0.14
against Earth's ~0.38), absorbed about 100 W/m² too much sunlight and
stayed too warm for a Pacific cold tongue or coastal upwelling. The
default `stratusScale` is 0.15, far below adiabatic, because on a flat
stability field the fit engages over every warm sea alike; it can be
raised once the eastern oceans carry an inversion for it to answer to.
A run's first step, before any boundary layer has been diagnosed,
carries no deck. `stratus: false` removes it, and both engines are
then bit-identical to the model without it, as they are wherever f or
the deck's water is zero; the column tests are in
`test/physics.test.mjs` and the engines' parity with the deck engaged
under an imposed inversion in `test/gpuModel.test.mjs`. Both engines
diagnose the boundary layer before the ∇⁴ closure, so the deck rests
on the same depth in each.

Tuning. 400-day N=8 runs (annual means) put the three knobs — the
cloud optical scale, the vapour coupling and the ocean heat transport
— at 60 m²/kg, 0.55 and Q₀ = 20 W/m²: 288.4 K, planetary albedo
0.304, ice on 21 % of the area, tropics 302 K, poles 256 K. **That
climate does not survive at N=16.** With eddies resolved the same
settings freeze: ice 25 % at day 100, 43 % and 274.5 K at day 400 and
still cooling, with 100 % ice cover to 30°S and a 100–170 g/m² cloud
band sitting on both ice edges. Cold air off the ice evaporates hard
over the open water, the cloud band reflects the sunlight at the edge,
the mixed layer under it freezes and the edge advances; the prescribed
ocean transport crosses zero at 35° and delivers nothing there. N=8
has no cold-air outbreaks, so it never sees this. The opposite knob
setting (scale 30, coupling 1.0) runs the other way at N=16: the ice is
gone by day 300 and the planet is at 298 K and warming at day 400.
Started from the warm moist day-100 state the defaults still freeze,
one hemisphere at a time (the south 100 % iced to 20°S at 2 m while the
north is ice-free at 301 K).

The instability is the ice–albedo feedback at Earth-unlike strength:
with ice on a fifth of the sphere and an albedo contrast of 0.53 a 10 K
warming frees ~25 W/m², eight times Earth's. Tuning cannot fix that;
it needs (a) ocean heat transport that responds to the ice edge
(diffusion of the mixed-layer temperature, the standard slab-aquaplanet
device, ~0.3 W/m²/K in energy-balance units ≈ 2 PW), (b) a smaller
albedo contrast at the latitudes that matter (zenith-angle-dependent
open-water albedo, ice nearer 0.5), and (c) tuning at N=16, where the
climate is the model's own. Tuning runs at N=8 are not representative.

### M11 — A stable ice edge: ocean diffusion and zenith albedo — done

Three changes to `js/physics/ice.module.js` and one to radiation,
aimed at the two feedbacks that ran away in M10:

- **Diffusive ocean heat transport.** Once per step, before the
  physics phase, `seaIce.prepare` forms the mixed-layer temperature
  (the surface temperature over open water, the freezing point under
  ice), takes its mesh Laplacian and stores the convergence
  `oceanDiffusivity` R² ∇²T plus the fixed Q-flux in a shared
  `oceanFlux` array that the cell update then consumes. Heat flows
  down the gradient, so an advancing ice edge pulls heat toward itself
  from the water beside it and an ice interior at uniform freezing
  point receives nothing. The coefficient is the energy-balance
  diffusivity in W/m²/K; 0.3 carries about 2 PW poleward at 35°, the
  real ocean's share. The Laplacian is the divergence of edge fluxes, so
  the convergence sums to zero over the sphere to roundoff
  (`test/ice.test.mjs`). With the layered ocean carrying the heat in
  its own currents the coefficient is the sub-grid eddy residual,
  0.01 W/m²/K (about 2000 m²/s over a 55 m mixed layer); at 0.3 it
  fed the freezing water at the ice edge 160 W/m² from the warmer
  water beside it and no winter ice could form there. In the layered
  ocean it is part of the mixed layer's tendency, which the workers
  compute with the other layers' (`parallel.module.js`), each layer by
  one thread, so the parallel step stays bit-identical to the serial
  one. The prescribed profile stays available but defaults to zero.
- **Direct and diffuse light at the surface.** Open water reflects the
  direct beam with the zenith-angle albedo of Briegleb et al. (1986),
  0.02 under a high sun, 0.07 at 60° zenith, 0.3 near the horizon, and
  diffuse light with 0.06. Radiation splits what reaches the surface:
  the direct beam is exp(−τ/μ) of the incident less a clear-sky
  skylight fraction of 0.15 (Rayleigh scattering the model does not
  otherwise have); the rest of what the cloud passes is diffuse; the
  multiple reflections between surface and cloud base are diffuse at a
  mean cosine of 0.6. Under thick cloud the surface therefore sees the
  diffuse albedo whatever the sun's height, so a low sun over open
  water near the ice edge is bright only in clear sky. Radiation
  exposes `cosZenith(i)` and the physics phase evaluates both albedos
  per cell.
- **Ice albedo 0.5**, the value of bare summer sea ice, instead of 0.6.
- **Tuning at N=16 only.** 400 days cost 9 minutes on 10 workers;
  `scratchpad/trisk/annual.sh <log>` prints the annual means,
  `zonal.mjs <state.json>` the 10° bands of temperature, ice and cloud,
  and `eke.mjs <state.json>` the eddy kinetic energy by band and
  hemisphere; the run log carries the per-hemisphere eddy energy.

Tuning at N=16 (400–500-day runs, means of the last year; cloud
optical scale c in m²/kg, ocean diffusivity D in W/m²/K, vapour
coupling 0.55 throughout):

| c | D | Ts (K) | solar − OLR | albedo | ice, annual (range) | note |
|---|---|---|---|---|---|---|
| 60 | 0.3 | 290.2 | +12 | 0.26 | 13 % (6–19) | still warming toward ~296; no runaway either way |
| 100 | 0.3 | 286.9 | +5.5 | 0.30 | 17 % (10–20) | |
| 150 | 0.3 | 284.0 | +3 | 0.32 | 20 % (13–24) | |
| 100 | 0.6 | 291.0 | +13 | 0.25 | 9 % (3–16) | summer cap melts away, warming |
| 120 | 0.6 | 289.9 | +11 | 0.26 | 10 % (4–16) | same |
| 115 | 0.45 | 289.6 | +6 | 0.28 | 11 % (7–18) | equilibrated ~290 |
| 115 | 0.45, slab ×2 (10 m) | 289.4 | +9 | 0.28 | 10 % (6–16) | summer storm track unchanged |
| **120** | **0.45** | 289.5 | +7, within 2 by day 500 | 0.28 | 11 % (6–18) | the defaults; 500 days, 12 min on 10 workers |

The diffusive ocean removed the runaway: at every setting above the
ice edge finds a seasonal equilibrium, with 2–3 m of ice surviving the
summer at the pole when D ≤ 0.45. The direct/diffuse split of the
surface reflection changed the climate by under 0.2 K; the warmth of
the c = 60 climate is the thin tropical cloud (6–20 g/m² under
Betts–Miller convection that rains without detraining condensate,
against 100–250 g/m² in the storm tracks), so the planetary albedo is
the lever and the cloud scale is what sets it. Doubling D melts the
summer cap and re-arms the ice–albedo feedback. Doubling the slab's
heat capacity does not change the seasonal cycle of the storm tracks:
the summer hemisphere's eddy kinetic energy at 250 hPa is a third to a
half of the winter's either way, which is Earth's northern-hemisphere
seasonality. Note that the slab's 2.1×10⁷ J/m²/K is only 5 m of
seawater (ρc_p = 4.1×10⁶ J/m³/K), a tenth of a real mixed layer, so
that test compared 5 m with 10 m; a 50 m layer (2.1×10⁸ J/m²/K) has a
600-day adjustment time and is the M13 ocean's upper layer.
Chosen defaults: c = 120, D = 0.45. Their climate at day 500 (northern
midsummer): tropics 299–300 K, summer subtropics 304 K, the winter
hemisphere iced to 55° with 3.5 m at the pole, the summer hemisphere
ice-free with its pole at 280 K, storm-track cloud 130–270 g/m² and
tropical cloud 4–40 g/m². Eddy kinetic energy at 250 hPa runs 500–680
m²/s² in the winter hemisphere and 210–290 in the summer one.


### M12 — Convective detrainment — done

Betts–Miller convection produced rain and nothing else, so the deep
tropics of the M11 climate carried 5–20 g/m² of cloud against 100–250
in the storm tracks, and the cloud optical scale had to be raised to
120 m²/kg to reach an Earth-like planetary albedo from the extratropics
alone. `convectColumn` now keeps the fraction `detrainment` (0.25) of
the condensate it produces in the column as cloud water, spread by mass
through the anvil — the layers from the level of zero buoyancy down
through `anvilDepth` (150 hPa) of pressure — and rains the rest. The
column's enthalpy change still equals the latent heat of all the
condensate, and vapour, cloud and rain sum to what was there
(`test/moist.test.mjs`). The anvil then lives under the existing cloud
physics: saturation adjustment evaporates it into subsaturated air and
autoconversion rains it out over its 3-hour lifetime, so a persistent
anvil needs the outflow layers near saturation, which the reference
profile's 70 % relative humidity does not guarantee.

Tuning (500-day N=16 runs, slab ocean, means of the last year): a
quarter of the condensate freezes the planet at every cloud scale
tried (60/80/100 → planetary albedo 0.44/0.48/0.51, 277/269/264 K and
falling, ice 24–42 %), because with no cloud fraction the anvil
overcasts its whole 965 km cell and a quarter of 5 mm/d with a 3-hour
lifetime is a 150 g/m² tropical overcast. A tenth, the all-sky
equivalent of a quarter under 40 % cloud cover, with cloud scale 40/50/
60 gives 293.8/290.9/288.5 K, albedo 0.27/0.30/0.32, ice 7/9/11 %.
Defaults: detrainment 0.1, cloud scale 60 (down from the 120 that the
cloud-free tropics had demanded).

### M13 — A shallow dynamic ocean — done, replaced by M18 (its module `js/ocean/reducedGravity.module.js` and the GPU port `js/gpu/ocean.gpu.js` were deleted once the layered ocean ran on both engines)

Designed against a panel of three (z-level free-surface, reduced-
gravity isopycnal, diagnostic Ekman): the z-level ocean's free-surface
wave (77 m/s at 600 m) sits near the explicit stability limit at the
atmosphere's step and needs seven new engine phases; the diagnostic
scheme has no momentum and is a dead end. The reduced-gravity model
has prognostic momentum, no barotropic mode, and fits the existing
main-thread ocean hook.

Two active layers over a motionless abyss: an upper layer (initially
50 m, the mixed layer whose temperature is the SST) and a thermocline
layer (350 m), with reduced gravities 0.02 and 0.01 m/s² at the two
interfaces. Each layer is the TRiSK shallow-water layer of M1 in the
vector-invariant form, driven by the Montgomery potential of a resting
abyss, M₁ = g′₁₂ h₁ + g′₂₃ (h₁ + h₂) and M₂ = g′₂₃ (h₁ + h₂), so the
fastest wave is internal (√(g′h) ≈ 2–3 m/s) and the ocean takes four
atmosphere steps at once with the RK4 of `integrators.module.js`.
Momentum forcing: wind stress ρ_a C_D max(|U|, gust) u on the upper
layer from the atmosphere's lowest-layer wind (zero under ice),
interfacial and bottom drag as linear stresses ρ r Δu (r = 2×10⁻⁴
m/s), and a ∇⁴ closure with its own 12-hour grid-scale e-folding —
the atmosphere's 3-hour coefficient at the ocean's 4× step is past
the RK4 stability limit at coarse resolution and was the first bug.
Heat: each layer carries h·T in flux form with centred edge
temperatures; the upper layer keeps the M11 diffusion (D R² ∇²T);
a layer thinner than 10 m entrains from below over an hour (a day at
first; with coastlines the upper layer of the Alboran Sea at N=64,
fed only through the Gibraltar cell, was drained by the wind faster
than a day's relaxation refilled it and the skin temperature over it
ran away), the thermocline layer from an abyss at 275 K (the one
exchange the ocean's heat budget does not close); the heat capacity the
surface sees is ρcp times at least a metre. Coupling, in `phases.ocean` on the main
thread: the ocean reads the SST from `surfaceT` over open water (the
freezing point under ice), steps, writes the SST back, publishes the
upper layer's heat capacity ρc_p h₁ per cell (a shared array that
replaces the slab's constant in the sea-ice cell update, so the
slab's 5 m becomes a real 50 m), and hands the heat converged under
ice to `oceanFlux` for the ice base over the steps until the next
ocean step. Workers never step the ocean; they read only the shared
capacity, so the parallel engine stays bit-identical. Snapshots carry
{h₁, h₂, u₁, u₂, T₂} and regrid across resolutions. Tests
(`test/ocean.test.mjs`): rest stays at rest bit for bit; westerlies at
45° give an equatorward transport within 9 % of τ/ρf in both
hemispheres; wind-driven flow conserves heat to 10⁻¹¹; heat converged
under ice reaches the ice base with the water at the freezing point;
the coupled model stays bounded. 168 tests pass.

Cost: at N=64 the ocean step takes 70–110 ms once every four
atmosphere steps on the main thread, 7–11 % of the 253 ms step. The
CPU engine's ocean also steps on the main thread but hands each
tendency evaluation to the workers a layer at a time, momentum and
tracers apart, from the 2D fields the main thread prepares (the free
surface and density gradients, the edge thicknesses and the layer
pressure potentials as column sums); with the ocean on, the parallel
step is about 3× the serial one at N=32 where it had been 1.6×.

First spin-up (N=16, 1500 days, the M12 defaults with the slab's
diffusivity 0.45 kept in the upper layer): the ice is gone by day
1050 and the planet sits at 292.8 K still warming, poles at 275 K,
currents 0.1–0.2 m/s, the upper layer pumped from 22 m under the
westerlies to 87 m in the subtropics, thermocline stable. The Ekman
cells themselves carry only ~0.1 PW (`scratchpad/trisk/oht.mjs`);
the warming is the mixed layer: 50 m of water holds its summer heat
through the polar winter where 5 m froze, and the diffusion tuned for
the slab now duplicates the transport the currents do explicitly.
Retuning the ocean's diffusivity with the ocean on (0, 0.1, 0.2).

**The wind stress was a third to an eighth of what the atmosphere
loses.** The zonal-mean surface winds are 2–3 m/s in every run against
Earth's 5–8, so the aerodynamic stress ρ_a C_D |U| U is 0.01–0.05
N/m². But the atmosphere's total momentum sink at the surface is
Earth-like, 0.1–0.27 N/m² in the zonal mean (`scratchpad/trisk/
stress.mjs`), because the Held–Suarez Rayleigh damping through the
lowest 30 % of the column — boundary-layer friction by another name —
takes three to eight times more momentum than the aerodynamic drag on
the thin lowest layer, and the ocean never saw it. The ocean now
receives the whole sink, `surface.stress`: the aerodynamic stress plus
Σ_k r_k u_k Δm_k over the damped layers, which is the momentum-
conserving coupling. The Ekman transport rises by the same factor
(the N=4 test's currents 0.026 → 0.146 m/s), and the ocean's heat
transport should approach the ~1 PW of Earth's wind-driven cells,
which is what would let the diffusion go.


### M14 — A diffusive boundary layer (`js/physics/boundaryLayer.module.js`) — done (tuning)

The last Held–Suarez placeholder was the Rayleigh damping of momentum
through the lowest 30 % of the column. It set the surface winds at 2–3
m/s in every run against Earth's 5–8: the circulation fixes how much
momentum must leave through the ground, and friction spread over that
much air lets the lowest layer give up its share without blowing. The
bulk surface fluxes then ran on the 3 m/s gustiness floor rather than
on the wind, the ocean received an eighth of the atmosphere's momentum
sink until M13 summed it explicitly, and heat and moisture from the
surface reached only the lowest 120 m layer until dry convection
carried them up.

The replacement is the diffusive boundary layer of Troen and Mahrt
(1986) that the simple moist GCMs use. Per cell, once per step in the
physics phase, `diagnose` finds the boundary-layer top as the height
where the bulk Richardson number of the lowest layer's virtual
potential temperature and wind, with the convective floor 100 u*² in
the shear, first exceeds 0.5 (interpolated between layers; u* is
√C_D times the lowest layer's wind with the gustiness floor), and lays
the K-profile κ u* z (1 − z/h)² over the layer interfaces below it.
Where the surface is warmer than the lowest layer, u* in that profile
becomes the unstable velocity scale of Holtslag and Boville (1993),
u* (1 − 15 ζ)^¼ with ζ = 0.1 h/L from the bulk surface buoyancy flux
(virtual, with a sea surface's saturation humidity; dry over land;
floored at −2), so a convective marine boundary layer mixes momentum
down to the surface instead of leaving the lowest layer to the drag;
stable columns keep the neutral profile.
The interface coefficients ρK/Δz live in a shared array. In the adjust
phase the cell units mix θ, q and qc down each column and new edge
units mix the normal velocity down each edge, both by implicit Euler on
the same tridiagonal system, which conserves each column's mass-
weighted total exactly and is unconditionally stable. Nothing mixes
above the top, and the search stops at σ = 0.5. The surface fluxes and
the aerodynamic drag stay explicit sources on the lowest layer, which
the diffusion spreads upward; the Rayleigh damping is off
(`pblRate` 0) whenever the boundary layer is on. Tests
(`test/boundaryLayer.test.mjs`): an unstable sheared column gets a
1.7 km boundary layer that mixes θ toward uniform and conserves θ and
q; a strongly stable column stays unmixed; momentum mixing brings wind
down to the surface layer and conserves each edge column's momentum;
serial and parallel engines stay bit-identical.

First N=16 run (500 days, with the ocean, the M12/M13 defaults):
the surface winds are Earth's — trades −6 m/s, westerlies +7 to +9 —
and the aerodynamic stress 0.07–0.19 N/m² is the whole momentum sink.
The ocean answers: currents to 0.8 m/s, the upper layer pumped to the
10 m floor under the westerlies and at the equator and to 120 m in
the subtropics, Ekman transports of 2 m²/s, and an upper-layer heat
transport of 1–2 PW poleward across the tropics with the mid-latitude
cells carrying ~1 PW equatorward. Climate: 287.3 K over the last year
but warming at +10–17 W/m², planetary albedo 0.33 with cloud water
50–77 g/m² (the mixing moistens the lower troposphere), ice 5 % of the
area (3–9 % with the seasons), latent heat 110 W/m².

Tuning with the boundary layer and the ocean (800-day N=16 runs, last-
year means). Removing the ocean's diffusion does not work: with D = 0
the planet cools to 280.5 K (cloud 60) or 278.5 K (cloud 80) with ice
on 28–31 % of the area and still cooling, and D = 0.15 with cloud 70
gives 282.2 K and 19 % ice, also cooling. The wind-driven cells carry
heat within the tropics and back toward the equator in mid-latitudes;
nothing dynamic reaches past 50°, which on Earth is the buoyancy-
driven overturning this two-layer ocean does not have, so the
diffusion keeps standing in for it. The other change is the cloud: the
mixed boundary layer keeps the lower troposphere moist and the cloud
water sits at 60–85 g/m² against 45 before, which with cloud scale
60–80 gives planetary albedos of 0.39–0.43. The retune therefore keeps
D at 0.3–0.45 and lowers the cloud scale to 45–55.

| D | cloud | Ts (K) | budget | albedo | ice, annual (range) |
|---|---|---|---|---|---|
| 0 | 60 | 280.5 | −3, cooling | 0.40 | 28 % (26–31) |
| 0 | 80 | 278.5 | −12 | 0.43 | 31 % (29–34) |
| 0.15 | 70 | 282.2 | −0.5, cooling | 0.39 | 19 % (16–22) |
| 0.45 | 60 (to day 1000) | 288.1 | +15 | 0.31 | 2 % (0.2–4) |
| 0.45 | 45 | 289.5 | +22 | 0.29 | 1.8 % (0–5) |
| **0.3** | **55** | 286.5 | +13 | 0.33 | 7.4 % (5–11) |

Defaults: D = 0.3, cloud scale 55. Its surface climate held at 286.5–
287 K over days 400–800 with Earth's sea-ice fraction and a 350/280
m²/s² winter/summer storm-track contrast; the +13 W/m² is going into
the 400 m ocean, which warms by ~0.25 K a year at that rate, so the
surface will drift up toward 288 K over decades. The warmer settings
lose their ice within a few years.

The page default is the end of a cascade from that run's day-800
state (N=32 for 50 days, N=64 for 25 days, 53 minutes on 10
workers): `runs/pbl64_state_day875.json`, 285.3 K, sea ice on 11 % of
the area at 0.7 m, planetary albedo 0.30, +8 W/m², currents to 2 m/s
and the upper layer entrained to a 70 m mean at N=64.

### M15 — The model on the GPU (`js/gpu/`) — done

The whole step runs on the GPU through WebGPU: in Node through Google's
Dawn (the `webgpu` package, a dev dependency, so the GPU kernels sit in
the same test suite as the CPU engine and never need a browser to be
checked), in the page through `navigator.gpu`. Everything on the device
is single precision.

- `core.gpu.js` packs the mesh into one integer and one float buffer,
  the level constants into a third, the state into one f32 buffer with
  the CPU layout inside it (the RK4 stages, the trial state and the
  tendency share it), the diagnostics into one scratch buffer and the
  physics arrays into another; every kernel binds the same eight
  buffers so a bind group is a choice of input and output. A tendency
  is eight dispatches — edge mass fluxes, layer divergences, the column
  (dπ/dt, πσ̇, the Exner and geopotential diagnosis, the drag rate),
  kite-weighted π on vertices, PV on vertices and on edges, the cell
  tendencies (kinetic energy, θ/q/qc transport) and the momentum
  tendency — then a fused advance; the ∇⁴ closures are two Laplacian
  passes per field, and dispatches beyond 65535 workgroups go through a
  second dimension.
- `physics.gpu.js` is the column physics one thread per column: the
  three-band radiation with clouds and the direct/diffuse surface
  reflection, the bulk fluxes, the sea ice, the boundary-layer
  diagnosis, and the adjustment (boundary-layer mixing by the
  tridiagonal solve, saturation adjustment, Betts–Miller with
  detrainment, autoconversion, the filler, the dry adjustment), with
  momentum mixing one thread per edge. `layeredOcean.gpu.js` is the
  ocean with its own state, stages and a binding that adds the
  atmosphere's state and diagnostics for the coupling kernels.
- `model.gpu.js` presents the CPU model's interface: double-precision
  mirrors refreshed by `sync`, an asynchronous `step`, `diagnostics`
  from the per-cell energy terms read back and summed on the CPU, and
  the ocean's initialize/load/serialize. The run driver takes
  `ENGINE=gpu`; the page takes `?engine=gpu`, the default when WebGPU
  is available, and says so on the status line.

Precision. The σ-coordinate column is the one place single precision
bit: the Exner span of a bottom layer is the difference of two numbers
near 1, and that cost 1e-5 relative in the layer Exner, 2.4 J/kg in the
geopotential and 3.6×10⁻⁶ m/s² in the pressure gradient. Two exact
factorings remove it: all σ-dependence of the Exner integrals is
precomputed in double precision as level constants (the layer Exner,
its π-derivative and the interface interpolation weights become one
well-conditioned power per column times a constant), and the
geopotential is carried as a deviation from a reference column built
from the initial mean θ profile, so gradients come from small numbers.
What remains, ~1e-7 m/s² in the momentum tendency at the top of the
atmosphere, is the f32 rounding of θ itself near 2000 K.

Verification (`test/gpu.test.mjs`, `test/gpuModel.test.mjs`): the
tendency matches the CPU core to 3×10⁻⁶ in dπ/dt and 4×10⁻⁵ in dθ/dt; a
rest state stays at rest; twenty dynamics steps agree to a surface-
pressure RMS of 0.008 Pa and 4×10⁻⁴ m/s in wind; one full step with
physics agrees to 1.8×10⁻⁷ in θ and 3×10⁻⁵ K in surface temperature;
twelve steps give the same mean temperature, absorbed solar and OLR to
the last printed digit; eight coupled steps with the ocean agree to
10⁻⁴ K and 10⁻⁶ m/s. A 100-day N=16 run from the tuned state tracks
the CPU run's climate (temperature within 0.2 K, the same ice cover,
fluxes within weather noise).

Cost at N=64: dynamics and closures 38 ms per step, the whole step 71
ms, against 291 ms on ten CPU cores — 79 simulated hours per minute at
the 3× step in Node, 40 in the pane's browser while sharing the GPU.
The gain is smaller at N=16 (1.4 against 1.9 minutes per 100 days),
where the grid is too small to fill the GPU.

Batched stepping. `model.stepBatch(count, dt)` records `count` steps
into one command encoder, submits it once and waits once. The eight
parameters every kernel reads change 15 times a step: the RK factors,
the sun's direction at each step's model time, the closure factors, the
ocean's. Each set is staged on the host and copied into the parameter
buffer by a copy command at its place in the encoder, and the staged
sets reach a ring buffer in one queue write just before the
submission. The ring holds 16384 sets, about 1090 steps; a longer batch
submits each time it fills. The device ends bit for bit where the same
steps taken one at a time leave it (`test/stepBatch.test.mjs`, at N=6
with land and the ocean, with batches that start off the ocean's
every-fourth-step cadence), and a day of `scripts/spinup.mjs` from a
saved N=64 state writes the same snapshot byte for byte either way.
The spin-up takes a day per batch (`BATCH`; 1 steps one at a time,
waiting every eighth step). With `RECORD` the forcing recorder's
per-step kernel is recorded into the batch after its step through
stepBatch's per-step callback, and writes the same day file byte for
byte. At N=64 batching gained nothing measurable.
The step is bound by the GPU: 32 ms of kernel time in pass
timestamps, against 0.3 ms to record it and 0.02 ms for an empty round
trip. Stepping one at a time already kept the device busy. (An N=128
spin-up shared the GPU during these measurements.) The page still steps
one at a time. Its worker keeps two steps in flight so that the globe's
frames never wait behind the model's work, and the pacer's pauses act
between steps (see the page as a client of the model worker). A batch
is a long queue by design, and it would buy the page no time.

### M16 — Land surface (`js/geography.module.js`, `js/physics/land.module.js`) — done

Continents come from `data/topography_0p25.bin`: ETOPO1 ice-surface
elevation (the top of the Antarctic and Greenland ice sheets) from
NOAA CoastWatch's ERDDAP, subsampled to 5 arc-minutes and block-averaged
to a 0.25° grid of int16 metres, 2 MB, 29.2% land by area;
`scripts/topography.mjs` regenerates it and `test/topography.test.mjs`
checks it. `createGeography` assigns every raster point to its nearest
cell by walking the neighbour graph from the previous point's cell and
takes each cell's mean elevation and land fraction; a cell is land when
more than half of its points are above sea level. Edges between two
ocean cells are the ocean's, edges between land and ocean are the
coast, which the page draws through a segment layer.

Land cells carry a skin of heat capacity 1e6 J/m²/K in place of the
ocean's upper layer, a Manabe bucket of 150 kg/m² of soil water whose
wetness β = min(1, soil/(0.75·150)) scales evaporation and which spills
what it cannot hold into runoff, and a snow cover in water equivalent
that precipitation builds when the lowest air is below 0 °C and the
surface energy melts, holding the skin at the melting point while it
does. Evaporation draws on the snow first, then the soil. Land albedo is
0.2 for both the direct and diffuse beams, rising to 0.55 over 20 kg/m²
of snow (the tuning below settled both; 0.25 and 0.7 held a snowy,
cold climate); the drag and exchange coefficients are 1.5e-3 over land
and water alike, as per-cell arrays that the surface drag, the
boundary layer's friction velocity and the radiation column's bulk
exchange all read. The ocean carries no flux, stress or diffusion
through edges that touch land, and its diagnostics average over ocean
cells. The land update sits in the physics phase; the deposit of the
step's rain as snow or soil water sits in the adjustment phase, where
the moist physics leaves each column's rain (`moist.rain`).

On the GPU the physics buffer carries the land mask, per-cell drag,
soil, snow and runoff; the physics kernel branches on the mask, the
adjustment kernel deposits rain, and the ocean kernels carry edge and
cell masks. Eight steps over a synthetic continent agree with the CPU
engine to 2e-4 K in surface temperature and 1e-4 kg/m² in soil and
snow. The parallel engine reproduces the serial model bit for bit.

Snapshots, saved states and the regrid carry `land: {soil, snow}`;
states without it start with half-full buckets and no snow, and sea
ice is zeroed on land when an aquaplanet state is placed on continents.
The page draws coastlines in the Atmosphere and Ocean modes, colours land in Satellite mode
from soil water (dry tan to wet green) with snow whitening it, and
adds Soil water, Snow and Elevation overlays; `?land=off` keeps the
aquaplanet and `?topography=<url>` takes another raster.
Beyond the globe's limb the Satellite mode draws the sunlit air as a
thin blue rim, an exponential column of about 0.5 % of the radius in
scale height lit with the same sun-elevation ramps as the surface's
twilight, so it reddens and ends where the globe's terminator does.
With the sun in the frame a camera flare is drawn over the picture in
screen space, a halo, a horizontal streak and a starburst on the sun and
five ghost discs on the line through the view centre, scaled by the
sunlight slider and faded out toward the frame edge and over the sun's
own width as it passes behind the limb.

The first 400-day N=16 run with continents (from the aquaplanet
`pbl16b` state) was stable but cold: planetary albedo 0.37 against the
aquaplanet's 0.30, Ts falling from 285 to 281 K into the northern
winter, snow on half the land. Four sweeps of 365-day N=16 runs with
terrain re-tuned it. Two starts bracket each configuration's
equilibrium: the cold `terr16` state (280.8 K) and the warm aquaplanet
state placed on the continents (286.7 K); the drift is Ts at the same
season a year later, since the top-of-atmosphere imbalance stays near
+15 W/m² regardless while the ocean's upper layer deepens (section
M18's note). Land albedo 0.2 and snow albedo 0.55 held in every run
after the first pair.

| run | cloud scattering | gas optical depth | start | annual Ts | albedo | drift K/yr | land T | snow on land |
|---|---|---|---|---|---|---|---|---|
| terr16 (0.25 / 0.7) | 55 | 5 | cold | 280.8 | 0.391 | | 263 | 53% |
| c45 (0.25 / 0.7) | 45 | 5 | cold | 279.7 | 0.378 | | 271 | 47% |
| a20s55 | 55 | 5 | cold | 280.5 | 0.360 | | 273 | 42% |
| c40a20s55 | 40 | 5 | cold | 282.3 | 0.334 | +1.5 | 276 | 37% |
| w_c40a20s55 | 40 | 5 | warm | 287.4 | 0.314 | −1.8 | 284 | 24% |
| w_c40g6 | 40 | 6 | warm | 287.7 | 0.310 | −1.0 | 285 | 23% |
| w_c40g7 | 40 | 7 | warm | 287.9 | 0.309 | −0.7 | 285 | 22% |
| w_c35 | 35 | 5 | warm | 287.9 | 0.302 | −0.8 | 285 | 23% |
| **w_c35g7** | **35** | **7** | warm | **288.6** | **0.296** | **+0.3** | 286 | 21% |
| w_c35g8 | 35 | 8 | warm | 288.7 | 0.299 | +0.3 | 286 | 21% |
| k_c35g7 | 35 | 7 | cold | 283.4 | 0.317 | +2.1 | | |

With an Earth-like albedo already reached at cloud scattering 40, the
remaining cold bias was longwave, and the well-mixed gas band's optical
depth — the model's CO₂ knob — closed it. The defaults are now cloud
scattering 35 and gas optical depth 7 (from 55 and 5), land albedo 0.2
and snow albedo 0.55: annual mean 288.6 K, planetary albedo 0.30, sea
ice 1% (0.1–3% over the year), land 286 K with snow on a fifth of it,
and the warm start drifts by +0.3 K a year toward an equilibrium at or
just above it. The tuned N=16 state cascaded to N=32 for 50 days and
N=64 for 25 days on the GPU; `runs/cont64_state_day1440.json` is the
page's default, and the GPU engine runs it at 72 simulated hours per
minute, the aquaplanet's pace.

### M17 — Terrain — done

The surface geopotential is the land elevation clamped at sea level,
relaxed twice toward the neighbour mean with ocean cells fixed at zero
so that coasts slope down to the sea over two cells, times g. The CPU
core already integrated the geopotential from it; the GPU core keeps
its geopotential as the deviation from a reference column and adds the
per-edge gradient of φ_s in the momentum kernel from a mesh-float
buffer filled in double precision.

The gate was the classic test: an isothermal atmosphere at rest over a
3 km Gaussian mountain at N=16 must stay at rest. Two changes to the
core were needed, and with them it stays at rest exactly, 0.000 m/s
over five days, with the GPU within 2e-4 m/s of it:

1. The pressure-gradient force's second term is R T̄ ∇ln π with T̄ the
   edge mean of the layer temperature θ_v Π_layer, in place of
   cp θ_v (∂Π/∂π) ∇π. The discretely isothermal state has
   θ_k = T0 / Π_layer,k, for which the geopotential of every σ layer is
   φ_s plus a column-independent constant, so its horizontal gradient
   is exactly ∇φ_s, and the ln π form of the second term is exactly
   −∇φ_s; the earlier form left the second-order cancellation error of
   Δπ/π across an edge, tens of Pa/m of it over a mountain.
2. The ∇⁴ closure diffuses the temperature θ Π_layer rather than θ.
   Along a σ surface over terrain θ varies as π^−κ, 25 K across a 3 km
   mountain even in an isothermal atmosphere, and diffusing that
   variation was the entire spurious wind: 3.1 m/s after five days with
   the θ closure, exactly zero without it. Diffusing the temperature
   leaves the isothermal state alone and differs from diffusing θ by
   about a percent of the closure elsewhere, where π varies by a few
   percent.

Surface pressure is no longer sea-level pressure. Frames carry
`mslp`, π reduced to sea level through a column at the lowest layer's
temperature plus half a standard lapse rate over the terrain height,
and the page's pressure overlay and isobars use it. Where a pressure
level lies below the surface, wind, temperature and humidity show the
column's lowest layer, its real near-ground values, and only the height
is extrapolated hydrostatically, the same operation as the sea-level
reduction, so that its contours remain a smooth pressure field rather
than a map of the terrain; at 1000 hPa that is most of the land. A state placed on a different terrain is rebalanced: π scales
by exp(−Δφ_s / (R T_bottom)), where Δφ_s is the target minus the source
terrain regridded to the target mesh, so a flat state loads onto
mountains without a shock; saved states and snapshots carry a
`terrain` flag. `?terrain=off` keeps flat continents.

### M18 — A layered ocean (`js/ocean/layered.module.js`) — done (spinning up)

The two-layer ocean of M13 could hold an Ekman layer and a gyre but had
nothing below to return the flow, so it grew no western boundary current
and its upper layer deepened ten metres a year. It is replaced by a
hybrid isopycnal ocean in the MICOM design: a bulk mixed layer with its
own temperature and salinity over 44 interior classes of fixed
reference density (`LAYER_DENSITIES`, **Classes** below), on the real
bathymetry `D` (the cell-mean ETOPO depth, at least 50 m, zero on land).
Every layer is a TRiSK shallow-water layer in the vector-invariant form
carrying thickness, edge velocity, heat h·T and salt h·S. Density is the simplified equation of state of Roquet et
al. (2015) with NEMO's nn_eos = 1 coefficients at the surface
(`js/ocean/seawater.module.js`):
ρ = 1026 − a₀(1 + ½λ₁Tₐ)Tₐ + b₀(1 − ½λ₂Sₐ)Sₐ − νTₐSₐ, Tₐ = T − 10 °C,
Sₐ = S − 35, a₀ = 0.1655, b₀ = 0.76554, λ₁ = 0.05952, λ₂ = 7.4914×10⁻⁴,
ν = 2.4341×10⁻³. The first version had six layers (σ = 24.0 to 27.7) and
a linear equation of state with α = 2×10⁻⁴ K⁻¹ everywhere; under it,
water at the freezing point was denser than warmer deep water, so polar
columns were only stable if the interior started at the freezing point,
the deep ocean sat at −0.7 and −1.8 °C, and every class warmer than
15 °C was missing from the tropical thermocline.

**Classes.** The classes lie where the World Ocean Atlas 2023 (**Start
from the World Ocean Atlas**) holds its water in this equation of state,
1.38×10¹⁸ m³ on its 1° grid, half of it between 1026.72 and 1026.91:

| Classes (kg/m³) | Spacing | Water | Atlas volume per class (10¹⁵ m³) |
|---|---|---|---|
| 1020.5, 1021.0, 1021.5 | 0.5 | the warm pool (1020.7) and tropical surface water | 1.6, 1.5, 2.1 |
| 1022.0–1026.0 | 0.25 | the thermocline, the 20 °C isotherm near the top of 1024.0 | 1.3–44 |
| 1026.1–1026.6 | 0.1 | mode and intermediate water, Antarctic Intermediate Water at 1026.27 | 30–77 |
| 1026.65–1026.98 | 0.03 | the Southern Ocean's winter water (−1.8 °C, 34.4 psu: 1026.79) and Labrador Sea water (1026.80) in 1026.80, Circumpolar Deep Water (+1 °C, 34.7: 1026.85) in 1026.86, North Atlantic Deep Water (1026.89) in 1026.89, Antarctic Bottom Water (1026.93) in 1026.92, Weddell shelf water (1026.95) in 1026.95 | 1.0–155 |
| 1027.05, 1027.15, 1027.3 | 0.07–0.15 | the Arctic's intermediate and deep water (−0.7 °C, 34.93) and brine-enriched shelf water | 2.9, 12, 0.2 |
| 1027.5, 1027.75, 1028.0 | 0.25 | the Red Sea's deep water (1027.54) and the Mediterranean's intermediate and deep water (1028.03) | 0.15, 0.13, 3.5 |

Between 50 and 70°S the atlas's water spans 1025.7–1026.95, 70% of it in
1026.83–1026.92; north of 50°N in the Atlantic sector it runs from 1025.8
through the dense band to the Nordic Seas' 1027.15. The label salinities
(`LAYER_SALINITIES`) are the tropical atlas's 34.0, 34.5 and 34.85 psu in
the warm-pool classes, 35 through 1025.0 falling to 34.8 at 1026.95,
34.82–35.05 in the classes from 1026.98 to 1027.3 and 39.9, 39.0 and
38.65 in the Red Sea's and Mediterranean's, so that every label
temperature lies above the freezing point: 28.6 °C at 1020.5, 19.4 °C at
1024.0, 2.8 °C at 1026.80, −1.7 °C at 1027.3 and 14.1 °C at 1028.0. At
65°S 0°E the N=64 atlas start holds the −1.4 °C winter water in
1026.71–1026.77 at 50–83 m and the +0.7 °C deep water in 1026.86 at
150–440 m, with 1026.80 (−1.0 °C) and 1026.83 (−0.1 °C) between them.

**Pressure force.** The pressure in interior layer k at height z is
P_{k−1} + ρ_k g (z_{k−1} − z), with P_{k−1} the pressure at the top of
the layer and z_{k−1} = η − H_{k−1} that top; its horizontal gradient is
the gradient of the cell potential

    Φ_k = g η + (g/ρ₀) [ (ρ_ml − ρ_k) h₀ + Σ_{j=1}^{k−1} (ρ_j − ρ_k) h_j ],

exact on the mesh whatever the bathymetry, and telescoping through a
layer of token thickness so the layers it separates feel the density
step across it. In a column where layer k lies below the bottom the
potential is the bottom pressure minus ρ_k g D, which is the Montgomery
potential MICOM assigns to a massless layer. The mixed layer's own
density varies, so its force is g∇η plus the depth mean of its density
gradient, −(g/ρ₀)(h₀/2)∇ρ_ml, applied at the edges.

**Edge thicknesses.** The mixed layer, present everywhere, takes the
centred thickness at an edge; an interior layer takes the smaller of the
two cells', so it never flows into a cell where it has no water, whether
it has outcropped there or the bottom lies above it (the fictitious
potential of a layer below the bottom would otherwise fling water off
every shelf). A column flows through an edge only within the water that
exists on both sides: the thicknesses at an edge are scaled so their sum
is at most min(D_a, D_b) + η. The mixed layer's mass and tracer flux
through an edge uses at most the thickness of the cell the water leaves,
so a deep convective column beside a thin mixed layer cannot drain its
neighbour below zero within a step (the centred thickness stays in the
momentum equation, and the Coriolis term uses the unlimited flux h_e u:
with the limited flux and a potential vorticity built on the centred
thickness, a 20 m mixed layer beside a 200 m one felt a fifth of the
Coriolis force and ran down the sea-level gradient at the Falkland
Plateau edge unbalanced, to the 5 m/s clamp). Heat and salt ride the mass flux with the donor
cell's temperature and salinity: a centred edge value let a step that
moves a large share of a thin cell's water leave the remainder at a
temperature outside anything in the neighbourhood (370 K off Chile on
day 34 of the first N=64 run). The edge potential vorticity is
(ζ̄ + f̄)/max(h, 20 m) with h the same edge thickness, except that an
interior layer's is at least `vorticityCentring` = ½ of its two cells'
mean thickness ½(h_a + h_b); the floor keeps the term bounded where a
layer thins to nothing. On the smaller thickness alone, a layer left a
few metres thick under a deep mixed layer beside the full layer next
door assembles a Coriolis force from its thick neighbours' fluxes many
times what its pressure gradient balances: a 12 m remnant under a
600 m mixed layer at 45°S ran to 1.3 m/s in five days at N=32 (0.13 m/s
with the centring), and the year-6 N=64 state had 36 interior edges
faster than 1 m/s, up to 4.4 m/s in a 9.5 m layer at 53°N. The full
centred thickness (`vorticityCentring: 1`) does as well there but moves
the western boundary current of the 80° test basin at N=16 two cells
off the wall; half of it leaves unchanged every edge whose thinner side
holds at least a third of the thicker's water. `vorticityCentring: 0` builds it
on the edge thickness alone. A layer thinner than 5 m at
an edge follows the velocity of the layer above, relaxing at the lesser
of 1/hour and 1/step.

**Free surface.** η = Σ h − D moves at √(gD) ≈ 240 m/s, far too fast for
the ocean's step (four atmosphere steps, 1350 s at N=64), so the depth
mean is split off. Each ocean step: η is frozen at its starting value
through one RK4 step of the layers; the barotropic transport U = Σ h_e u
and η take RK4 sub-steps forced by the depth integral of the slow
tendencies (every layer tendency with its g∇η removed, minus the
Coriolis force on U), with U ← −g H_e ∇η + f×U + slow and η ← −∇·U,
H_e the same edge thickness sum, at a wave Courant number of 0.35 on the
shortest edge over the deepest cell (eleven sub-steps at N=64); the
layers are then rescaled to the sub-step-averaged η and shifted so
their transport equals the averaged U. Before the rescale an interior
class left below its 0.01 m token is made up to it at its label from the
mixed layer, which gives up that water with its heat and salt, so the
tokens cost the column nothing. Velocities are finally clamped
to 5 m/s and the count of clamps reported (`oceanLimited`, zero in every
run so far).

**Mixed layer.** After each step, per column, driven by the surface
buoyancy loss of the step, B = g(α ΔT h + β ΔW)/Δt, with ΔT the surface
update's cooling of the mixed layer at the thickness h its capacity
used and ΔW the salt (psu·m) that evaporation, rain, runoff and ice
growth add, remembered over `buoyancyMemory` = 1 day as
B̄ ← B̄ + (B − B̄) Δt/day, so a day's sunshine does not undo a winter's
convection. Each interior class stands for water spanning half-way to
its neighbours' labels, from ρ_k − δ⁻ to ρ_k + δ⁺. The mixed layer
(1) swallows a class it is denser than all of (ρ_ml ≥ ρ_k + δ⁺, static
instability), at up to `convectiveRate` = 100 m a day;
(2) erodes the first class beneath it while it is within that class's
span and B̄ > 0, at the rate of convection into stratified water,
dh₀/dt = B̄/(h₀ N²) with N² = g (ρ_k + δ⁺ − ρ_ml)/(ρ₀ h_k), the class's
remaining span of density over its thickness — `convectiveErosion:
false` swallows every class no denser than the mixed layer instead;
(3) entrains the first layer below at the Kraus–Turner wind-stirring
rate w = 2 m u*³/(h₀ Δb) with m = 0.8 exp(−h₀/100 m), so the wind's
stirring fades below the depth it can reach, and Δb the buoyancy step
to that layer, at least 10⁻³ m/s² (a 0.1 kg/m³ step; the KT rate is
unbounded as Δb → 0); all three within `maximumMixedDepth` = 600 m (the
Southern Ocean's mode-water mixed layers reach 500–700 m; the deepest
class, one label for everything below about 1600 m, would otherwise let
a polar mixed layer that reaches its density erode it to the bottom at
its thickness's weak N², and the subpolar North Atlantic's 1000–2000 m
Labrador Sea convection is cut to the cap) and,
when `mixedNeighbourRatio` is set, within that multiple of the mean
mixed-layer depth of its sea neighbours. Otherwise it holds its depth,
convectively neutral or not. It detrains at once anything beyond
600 m, and, over a day and never above 50 m (`shallowestMixedDepth`,
also the minimum thickness, refilled from below), anything beyond its
neighbours' reach and anything beyond the Monin–Obukhov depth
2 m u*³/(−B̄) while B̄ < −10⁻⁹ m²/s³, which is how a winter's mixed layer
retreats in spring and leaves its water behind as mode water.
Detrained water goes to the interior layer whose density is nearest its
own, so water swallowed from a layer returns to that layer (the first
rule, the first layer at least as dense, put 900 m of 1027.24 water into
the 1027.7 layer and made a 0.5 m sea-level step with 2.5 m/s currents
around it; Bleck's mass-conserving split between the two bracketing
layers ratcheted a tenth of a metre of every swallow-and-return cycle
into the denser class); water lighter than the first interior layer
goes there whole (the tropical and summer mixed layers, a known slow
drift). Every warming day pulls the tropical mixed layer to the 50 m
floor, where 43–75 % of the equatorial cells east of 120W sit (N=128,
day 274). With a 20 m floor an earlier set of rules lost 2 K of SST a
month; with the present ones 60 coupled N=64 days of boreal spring keep
the equatorial mixed layer at 21–25 m, but warm the 20S–20N sea surface
by 0.3–0.5 K and, unless the South Equatorial Current is at its
strongest, the cold tongue by 0.4–1.2 K against the 50 m floor on the
same days, the thin layer resting on the warm water above an eastern
thermocline still 70–80 m deep (M21 item 4). `mixedNeighbourRatio` (off by default) caps the depth at that
multiple of the neighbours' mean and detrains the excess over a day; at
3 it moved the Southern Ocean's depths by less than 10 m in a 60-day
N=64 test, the pressure force and the centring of the potential
vorticity keeping even an isolated 600 m column balanced.

Earth's winter mixed layer is 100–400 m deep in the Southern Ocean,
500–700 m in the mode-water regions of the south-east Pacific and
Indian sectors, 300–1000 m in the subpolar North Atlantic, 100–150 m in
the subtropics and 20–60 m in the tropics (de Boyer Montégut et al.
2004). Earlier rules held it far shallower: the mixed layer returned
everything below 50 m to the interior whenever it was as dense as the
water beneath it (`neutralSnap`), kept at most 200 m, and swallowed any
class no denser than itself at 100 m a day. Those rules were added when
a kilometre-deep mixed layer beside 40 m ones opened a 5.7 m sea-level
hole off Cape Farewell and 500–900 m ones ran the Drake Passage at
5 m/s, but both happened while the mixed layer's Coriolis term used its
donor-limited flux, so that where a thin mixed layer fed a deep one it
felt a fraction of its Coriolis force (the fault behind the Falkland
Plateau jet that the centred flux h_e u removed); the pressure force,
(g/ρ₀)(h₀/2)∇ρ_ml, is the exact layer mean for any thickness, and a
600 m mixed layer beside 50 m ones moves the free surface by its steric
deficit (7 cm) and the water by 2 cm/s at N=16. The snap pinned the
June Southern Ocean of year 6 at exactly 50 m in 70 % of its cells,
over water within 0.03 kg/m³ of their own density down to a median
480 m; the snap and return each step churned that water through the
mixed layer. Letting a neutral layer hold, but still swallowing
whole classes, turned the pinned layer into zonal stripes of 600 and
100 m, one per class outcrop, as each 700–1000 m class at 50–55°S went
whole where the surface water just passed its label; eroding the class
at B̄/(h₀N²) instead gives, 60 days on from that June at N=64, 180 m
(90–360 m, 10th–90th percentile) at 40–65°S, smooth on the large scale
and deepest, to the 600 m cap, in the south-east Pacific west of Chile
where Subantarctic Mode Water forms, 100–150 m under the Antarctic ice
and 50 m in the northern summer; a tenth of the Southern
Ocean's columns still differ by more than 25 m from all their
neighbours, most of them in those deepest regions, where the depth
turns on how near the surface water is to the density of the whole
class beneath it.

**Salinity.** Evaporation minus rain is accumulated per cell over the
atmosphere steps between ocean steps, and each land cell's runoff flows
down the terrain to the sea: to the lowest of its neighbours that is
lower than it, or, from a pit, to the lowest lower cell two steps away,
then three, and so on (`runoffOutlets`, built once per model; at N=32
the largest outlets are the Ob, the Amazon, the Río de la Plata, the
Congo and the Nile), all applied as virtual salt fluxes to the mixed
layer; ice growth
rejects brine and melt freshens with an ice salinity of 5. The freezing
point is still the fixed 271.35 K of M9.

**Start.** From rest: the mixed layer 60 m deep with the atmosphere's
initial surface temperature and a salinity 34 + 2 exp(−((|φ|−25°)/20°)²);
each interior layer starts at its label salinity (`LAYER_SALINITIES`,
**Classes**) and the temperature that gives its label density there (the
1026.95 class at about 0.8 °C), blending poleward of 45° over 15° of
latitude toward −1 °C at the salinity that keeps the density (the 1026.6 class at
34.21 psu, the deepest at 34.65), as polar oceans hold cold, fresh water
on the density surfaces of the warm subtropical thermocline; so its
density is its label; interior layer bases (`LAYER_BOTTOMS`) in the
subtropics from 65 m under 1020.5 through 300 m under 1024.0 and 1061 m
under 1026.6 to 3550 m under 1026.92, scaled by 0.7 + 0.6 cos²φ toward
the poles; every layer lighter than the local surface water outcropped;
the 1026.95 class, whose base, like those of the denser classes, lies
below any sea floor, filling the column to the bottom, the denser classes
tokens beneath it and the deepest layer holding water giving up their
metres. The free surface starts level in pressure at 3500 m across the
open ocean's abyss, the cells at least that deep joined through such
cells over the largest area (`abyssalCells`), and every other column
takes the sea level around it by relaxation, so a basin behind a sill
starts at the ocean's sea level: levelled at 3500 m in their own deep
cells, the Mediterranean with its 1028.0 water behind Gibraltar would
start 5 m low and the Arctic behind Fram Strait 2 m low (N=64 and 128).
Polar surface water at the freezing point and 34 psu (1026.47) floats on
the 1026.95 class; brine that raises it past about 34.6 psu sinks into
it, as bottom water forms. With
35 psu in every class and at the poles, a first spin-up year lost all its
sea ice: polar mixed layers denser than the deepest class convected
2.4 °C water up all winter. A
saved ocean on other classes is carried onto these (**Loading across
class lists**), and one whose classes are unknown loads as this
climatology.

**Start from the World Ocean Atlas.** Fresh starts of `scripts/spinup.mjs`
and the page take their ocean instead from the World Ocean Atlas 2023
(NOAA National Centers for Environmental Information, public domain;
Locarnini et al. 2023, Reagan et al. 2023): the annual-mean objectively
analysed temperature and salinity of all decades on the 1° grid, at 25 of
its standard depths (0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300,
400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000, 3000, 4000 and
5000 m). `scripts/packWoa.py` writes them from the atlas's NetCDF files
to `data/woa_annual_1deg.bin` (6.5 MB): a header with the grid, the
depths, the scales and a provenance line, then potential temperature (the
atlas's in-situ temperature brought to the surface by the UNESCO
algorithm of Fofonoff and Millard 1983 at the pressure of Saunders 1981)
and practical salinity as int16 in steps of 0.001 °C and 0.001, −32768
where the atlas has no water. `js/ocean/climatology.module.js` reads it
in Node and the browser; `columnAt(lat, lon)` interpolates each level
bilinearly between the four surrounding grid points when all four are
wet, takes the nearest wet one otherwise, ends the column at the first
level none of them reaches, and over atlas land takes the nearest wet
column within 5°. Each sea cell's column (`atlasColumns` in
`js/ocean/layered.module.js`) is sampled down to its bottom; where the
atlas column ends a standard level or more above it, as off coasts that
the cell's mean depth runs past, the levels below come from the nearest
atlas grid point that holds them (`columnReaching`, within 30°), and the
atlas's deepest values are carried on below the deepest level. Held from
the cell's own last level instead, 165 columns within 25° of the equator
at N=64 started more than 70 % lighter than 1025 over more than 300 m,
full-depth plugs of surface water such as 1020.5 from 50 to 1194 m at
1S 100E, where the ocean reached its speed limit; with the deeper column
13 remain, 300–430 m columns of the South China, Timor and Coral Seas whose
stratification the atlas itself holds; its potential density, held from
decreasing with depth, places each interface where it crosses half-way
between the neighbouring labels, since a class holds the water between
those midpoints: the lightest class takes all lighter water, the densest
all denser, a class with no water keeps its token and the densest class
present reaches the bottom. The mixed layer reaches to where the density
first exceeds the surface's by 0.03 kg/m³ (the criterion of de Boyer
Montégut et al. 2004, from the surface), within 50–600 m, and every layer
takes the atlas's mean temperature and salinity over its depths. So
averaged, the water lies off its label: at N=32 73% of the classes
holding water are within the 0.01 kg/m³ restoring tolerance (97% within
0.05); the lightest class is a median 0.47 light, up to 1.9 in the
tropical open ocean and 15 in the brackish Baltic, the densest at most
0.17 dense in the Mediterranean, and the rest a median 0.002 dense from
the curvature of the equation of state. A class further off than the
tolerance therefore keeps the atlas temperature and takes the salinity
of its label, as `interiorWater` does (below the lightest class and above
the densest, 99% of the classes move by less than 0.15 psu). Under the
initial ice the mixed layer is at the freezing point and elsewhere never
colder, and over open water its temperature becomes the atmosphere's
initial sea surface: `initialize` writes it into the surface temperature
it is given and returns how many sea cells took the atlas and how many
the analytic start (at N=64, Lakes Superior and Eyre, more than 5° from
any atlas water), and the GPU model sends that surface to the device.
Poleward of 72°N the initial ice lies over the atlas's warm Atlantic
water in the Barents and Greenland Seas, where the salty mixed layer at
the freezing point is denser than the class beneath and convects, as the
Mediterranean's does. The start puts the top of the 1024 class (about
the 20 °C isotherm) at 173 m at 0°N 160°E and 50 m at 0°N 100°W, where the
atlas crosses at 171 and 42 m, and outcrops it poleward of 35°; the 1025
class top lies at 210 and 203 m at 30°S 90°W and 30°N 150°W, where the
atlas crosses. The 60–70°S column averages −0.5, 0.0, 1.2, 1.1 and
0.5 °C over 0–60, 60–200, 200–500, 500–1000 and 1000–3000 m against the
atlas's −0.5, −0.1, 1.2, 1.2 and 0.5 at the same cells, where the
analytic start has 0.8 °C in its mixed layer and −1 °C at every depth
below; the interior's mean temperature is 3.5 °C. Five coupled days
from it at N=64 on bl34 (the ocean every eighth step) kept currents under
1.3 m/s with no clamped edge and ended with the warm pool at 27.9 °C,
the cold tongue at 21.8 °C and the 1024 class top at 172 m in the west
Pacific and 63 m in the east. `spinup.mjs` reads the file for a fresh
start without FROM, `CLIMATOLOGY=<file>` another, and `CLIMATOLOGY=none` starts
from the analytic climatology; the page fetches it for `from=none`,
`?climatology=off` keeps the analytic start and `?climatology=<url>`
takes another file. The option is `ocean: { climatology }`, a decoded
climatology (`loadClimatology`), or for the GPU and parallel engines a
path or URL; null, the default, is the analytic start byte for byte.
`test/climatology.test.mjs` checks the packing, the interpolation, the
N=32 columns and these values, and the GPU start against the CPU's.

**Interfaces.** `advance(surfaceT, ice, oceanFlux, stress, dt, concentration)` as in M13, the wind stress reaching the water under ice scaled by 1 − A(1 − `iceStressTransmission`) with the transmission 0.8, since drifting floes pass most of the air's stress to the ocean and the Antarctic Divergence upwells its deep water under a pack that covers it most of the year,
plus `accumulate(evaporation, rain, dt, runoff)` each atmosphere step;
`fields()` gives the frame its mixed-layer depth, SST, SSS, surface
velocity, thermocline depth (the top of the 1024.0 class, labelled
19.4 °C under the 20.4 °C of 1023.75, close to the 20 °C isotherm) and η;
`serialize()` is {h, u, T, S, eta} flattened layer-major with the class
list `densities`, and
`regridOcean` carries each layer's cell fields by masked interpolation
over sea tiles and its velocities by vector interpolation. Diagnostics
add the thermocline depth, mean salinity, the largest |η| and the
strongest depth-integrated transport through an edge in Sverdrups. The
page shows Mixed layer depth, Thermocline depth, Salinity and Sea
surface height under the ocean overlays.

**Checks.** A flat uniform column over a bumpy bottom on the real
coastline stays at rest to 10⁻⁶ m/s; a 1 m surface bump spreads as a
gravity wave; heat and salt are conserved to the last digit without
surface fluxes; in a closed basin 80° wide under a subtropical wind the
northward return flow sits in the two westernmost bands after 40 days
with a weak southward interior, the Stommel picture. Instabilities met
and removed on the way: PV at a vertex divided by a vanishing
thickness, a 2000 m column pouring into a 30 m shelf cell, layers pushed
onto a shelf by the potential of a layer below the bottom, forward-Euler
Coriolis in the sub-steps, and a token-layer relaxation faster than the
step at N=8. Cost on one CPU thread: 93 ms per model step at N=16
against about 30 before, the ocean now the larger part.

**GPU.** `js/gpu/layeredOcean.gpu.js` runs the same step on WebGPU: the
baroclinic RK4 over all layers, the barotropic sub-steps batched into one
submission, the mixed-layer exchanges, salt and the surface write-back,
with initialisation, loading and serialisation kept in double precision
on the CPU side. One ocean step at N=64 takes 16 ms on the GPU against
389 ms on a CPU thread, and `test/layeredGpu.test.mjs` holds the two
engines to single precision over one and twenty steps. Its freshwater
kernel takes the rain accumulator's increment since its last call, so
the precipitation diagnostic keeps working.

**Loading across resolutions.** A saved ocean carried to another mesh
is fitted column by column (`fitColumns`): the interpolated layers keep
their interface depths from the top down and are cut or extended at the
bottom to the new bathymetry plus the interpolated sea level (held
within ±5 m); a sea cell that arrives without usable water (a coast the
finer grid resolves, a stencil with no sea source) takes the climatology
column; the mixed layer keeps its floor; a column that already fits is
left exactly as it is. Without this a cell off Sydney held 4038 m of
water over a 900 m bottom and new sea cells had a 0 K surface, and the
first N=128 step from the N=64 state went non-finite. Basins the finer
grid connects through a strait the coarser one closed (the Red Sea at
Bab-el-Mandeb, the Gulf at Hormuz, the Gulf of Finland) arrive with
their own sea levels and surge through the new strait at the clamp for
about two days, then settle (2.3 m/s at Bab-el-Mandeb after 2.5 days
from the day-930 state).

**Loading across class lists.** A saved ocean carries its class list
(`densities`, float64 in the binary state); one without it is on the
seven classes of the page's saved runs or the 23 that followed them,
told apart by its layer count (`UNLISTED_LAYER_DENSITIES`). One on other
classes is carried onto these column by column before it is fitted
(`rebinOcean`, in `load` on both engines, so `spinup.mjs` continuing or
seeding FROM, `oceanSpinup.mjs`, an ocean handed over by
`oceanHandOff.module.js` and the page all take it). The mixed layer
copies through. An old class stands for water spread over its span,
half-way to its neighbours' labels, the lightest and densest over the
same width centred on their own water's density since they hold all
lighter and all denser water; each new class takes the share of that span
within its own, tilted linearly so that the shares' labels average to
the water's density. Each share keeps the old class's temperature and
salinity moved by the least change, counted in 4 K and 0.5 psu, that
takes the density from the shares' mean label to its own: warmer or
colder in the thermocline, saltier or fresher in cold water. A class
left under its token is made up to it from the column's thickest class,
and each new class moves with the old ones its label span overlaps. Every
column keeps its water, heat and salt to rounding; from a 24-layer
ocean at N=8 the median class lies within 10⁻⁴ kg/m³ of its label and the
worst 0.002 off, from an 8-layer one 0.0015 and 0.03, and the restoring
takes it from there (`test/layeredState.test.mjs`). An ocean already on
these classes loads as it did.

**Drag.** Interfacial drag is the linear stress ρ₀ r Δu between each
layer holding water at an edge and the nearest ones above and below that
do, with r = `interfacialDrag` = 2×10⁻⁴ m/s. Each layer divides the
stress by its own edge thickness, no less than THIN for an interior class
and no less than `minimumThickness` for the mixed layer, which divides its
wind stress and bottom drag the same way, so the momentum one layer loses
the other gains wherever the mixed layer is that thick at the edge; a
divisor of `minimumThickness` on every layer made an eastward column
source of 0.9–1.3×10⁻⁵ m²/s² on the equator at 180–100W (N=64, day 365),
a quarter to three quarters of the wind's. Under the option
`shearMixing`, r is instead the Pacanowski and Philander (1981)
viscosity ν = ν₀/(1 + 5 Ri)² + ν_b, ν₀ = 10⁻² m²/s (`shearViscosity`)
and ν_b = 10⁻⁴ m²/s (`backgroundViscosity`), over the distance Δz between
the two layers' middles at the edge, and no less than `interfacialDrag`;
Ri = Δb Δz/|Δu|² takes the buoyancy step between them (the mixed layer's
own density at the edge, the labels of the classes) and the difference
of their full velocities, the tangential part reconstructed as the
Coriolis term's, and r is held to half of what the thinner of the two
can take in one ocean step (M21 item 4). Under the option
`interiorShearMixing` (off) the same viscosity applies between two
interior classes with no `interfacialDrag` floor, so a quiet interface
takes ν_b/Δz, 2–4×10⁻⁶ m/s over 25–50 m, where the constant 2×10⁻⁴ m/s
is a viscosity r Δz of 5×10⁻³ to 10⁻² m²/s through the thermocline
against Earth's 10⁻⁴ to 10⁻³ below the undercurrent core; the mixed
layer's base keeps `interfacialDrag`. In 60 coupled N=64 days it let
thermocline classes run at 1.5–3 m/s beside the coasts near the equator
(M23). Bottom drag is quadratic,
C_D |u| u with C_D = 3×10⁻³, applied to the deepest layer with water at
the edge. A linear bottom drag of
2×10⁻⁴ m/s let a 28 m bottom layer slide down the flank of a seamount
at 4.5 m/s (60 times weaker than the quadratic law at that speed); the
quadratic law holds such layers near the gravity-current speed
√(g′h) ≈ 0.4 m/s.

**Volume and energy.** The GPU rescale adds a small correction,
h ← h + h·(η̄ − (Σh − D))/Σh, rather than multiplying each column by
(D + η̄)/Σh: in single precision that product was biased upward by about
5×10⁻⁵ m per column per ocean step, a mean sea-level rise of 1.2 m a
year that also carried heat in at each layer's temperature. The
correction keeps the drift near 2×10⁻⁷ m per step with no bias, and
`test/layeredGpu.test.mjs` holds the volume over 400 wind-driven steps.
The atmosphere returns the kinetic energy it dissipates as heat, in the
cells whose edges lost it: the surface and top drags inside the RK4
tendency, and the ∇⁴ momentum closure and the boundary-layer momentum
mixing through a per-edge record of each layer's loss of u² that a
cell phase gathers at the end of the step (the mixing's loss is shared
between layers in proportion to the shear dissipation at each interface
and each layer's own increment, so the heating is never negative). At
day 930 that is 2.66 W/m²: 1.57 from the drags, 0.57 from the closure
and 0.52 from the mixing. Snowfall releases its heat of fusion into the
lowest layer and sublimation takes it from the ground. From the day-930
state the whole-model budget then closes to −0.2 W/m² over a day, of a
time-mean top-of-atmosphere surplus of 13.0 W/m², 12.8 of it into the
ocean; what remains is the adiabatic dynamics gaining 0.29 W/m², the
temperature closure losing 0.21 (it is area-weighted rather than
mass-weighted over terrain) and the boundary layer's mixing of θ rather
than enthalpy gaining 0.11. Before these, 2.57 W/m² simply vanished.
The ocean uptake is spin-up: the volume-mean
ocean temperature is 1.0 °C against Earth's 3.5 °C, the two deepest
layers sit at −0.7 and −1.8 °C, and because the interior layers push
on the flow only through their label densities their temperature and
salinity are passive, so the heat they absorb has no dynamical brake
and the interior takes decades to centuries to come into balance.

**Ocean-only spin-up.** Those decades are cheaper without the
atmosphere. `RECORD=<dir>` on `scripts/spinup.mjs` writes one file a
model day, `forcing-DDDD.bin` (`js/forcing.module.js`; 2.3 MB at N=64,
so 0.8 GB a year), of the daily means of what the ocean and sea ice
received: the stress on every edge as the ocean steps used it (the
ocean's own stress buffer, after the ice's transmission), the net
surface heat flux with its absorbed and downward shortwave and the
sensible heat, evaporation, precipitation, the snow falling on the sea,
each land cell's runoff, and the surface temperature, SST (freezing
under ice), ice thickness and concentration. A kernel sums them after
every step (`createForcingRecorder`, `js/gpu/forcing.gpu.js`), rain
among them from each step's own total (PH STEPRAIN); runoff is the
change of the land's runoff tally over the day.

`scripts/oceanSpinup.mjs` loads a coupled snapshot (STATE) into the GPU
model and loops the first DAYS_PER_YEAR (365) recorded days over it for
YEARS years, never stepping the atmosphere (`createForcedOcean`). At
every atmosphere step of dt = 1350·16/N s a kernel runs the physics
kernel's own sea-cell surface update, spliced in as the same WGSL text
(`SEA_SURFACE_WGSL`, and `snowOnSea` for the snowfall), with the
recorded net flux in place of the computed one; the lead/ice split
takes the recorded downward shortwave times the ice's diffuse albedo
contrast. Every fourth step takes the recorded freshwater through the
ocean's accumulation and runoff routing and steps the ocean as
`advanceCoupled` does, with the recorded stress copied into its stress
buffer. Surface temperature, ice and concentration stay in the
atmosphere's buffers, as in the coupled model. A prescribed flux does
not answer the surface temperature, so the open water is restored to
the day's recorded SST by −RESTORE·(SST − SST_rec), RESTORE = 30 W/m²/K
by default (about 95 days on a 60 m mixed layer); a partly iced cell's
leads take it per unit lead area, and the ice itself only the recorded
flux. Replaying two recorded days at N=6 from the same state ends
within 0.2 K of the coupled run's SST everywhere
(`test/oceanSpinup.test.mjs`). A daily mean has no diurnal cycle, so
the mixed layer's buoyancy loss, and the exchanges it drives, differ:
after those two days the mixed layer is up to 6 m off and a few thin
interior layers the coupled run filled stay empty.

At each year's end the driver saves `<TAG>_yearYYYY.bin`, a full state:
the atmosphere and land of STATE as loaded (the land keeping the snow
on the sea ice), the spun-up ocean, surface temperature, ice and
concentration, and the time advanced by the years run, so that with a
365-day cycle recorded from STATE's day the season still matches the
atmosphere. It logs the ocean line of spinup.mjs, the 60–70S mean
temperature over 0–60, 60–200, 200–500 and 500–1000 m, the SST and
ice extent against the recorded last day, the global mean SST with its
drift since STATE and over the last ten years, the ocean's interior
temperature, fastest current, largest transport and clamped count, and
the wall time, and exits with 2 on NaN. A coupled run
starts from it as a snapshot: copied into the coupled run's OUT as
`<tag>_dayDDDD.bin` (DDDD its `day`), `scripts/spinup.mjs` continues
from it with the carried atmosphere, which then adjusts to the new
surface; LAND_FROM is not needed, since the land is already in the
file. At N=64, with the N=128 spin-up sharing the GPU, a looped day took
3.1–4.2 s against 9 s for a coupled day under the same load (4.5 s
alone); nearly all of it is the ocean step itself.

**Asynchronous coupling.** The deep ocean needs centuries and the
atmosphere does not, so `scripts/asyncSpinup.sh` alternates the two
kinds of run at one resolution N. Each of CYCLES cycles (9) is a
coupled phase of COUPLED_YEARS (10) model years, run as
`scripts/spinup.mjs` segments that end on the PER_YEAR (4) snapshot
days of the model year as `pairedSpinup.sh` places them, recording the
forcing over its last RECORD_YEARS (1) years, and then an ocean-only
phase of OCEAN_YEARS (100) years looping that record. A second
invocation with the other N runs the other resolution. With the coupled
run starting at day D₀ (its newest `<PREFIX><N>_dayDDDD.bin` in OUT, or
a fresh start at 0), cycle c's coupled phase ends at
E_c = D₀ + 365·c·COUPLED_YEARS and records days E_c − 365·RECORD_YEARS + 1
to E_c into `<PREFIX><N>_cCC_forcing/`. The ocean-only phase
(`oceanSpinup.mjs`, TAG `<PREFIX><N>_cCC`) starts from the snapshot at
E_c, whose season matches the record's first day, and its last year's
file hands the ocean back: the next coupled segment continues from the
snapshot at E_c with `OCEAN_FROM` set to it, which replaces the ocean,
sea ice and sea surface (surface temperature, ice, concentration and
the snow on the ice over the sea cells) and keeps the atmosphere with
the deck's carried state (its running-mean subsidence, inversion height
and gate), the land cells, the day and the time of the coupled snapshot
(`js/oceanHandOff.module.js`). The coupled calendar
therefore never counts the ocean-only years: the coupled run resumes at
the day it stopped, and its days stay consistent with the season, since
the ocean-only phase covers whole years of the same cycle. The
ocean-only files' own `day` runs on from E_c by the days they ran and
is outside that calendar. The snapshots carry `oceanYears`, the years
their ocean has spent alone in all, and `oceanFrom`, the file whose
ocean they took; a segment whose snapshot already names its OCEAN_FROM
continues without replacing it again.

The script keeps no state of its own beyond `<PREFIX><N>_cycles.txt`
(the start day D₀, the cycles done and a convergence stop) and reads
where it stands from the files in OUT, so it can be started again at
any point. After each ocean-only phase it writes a summary line to
`<PREFIX><N>_cycles.log`, from the last year's lines of the phase's log:
the cycle, the coupled and ocean-only years done, the ocean line (warm
pool, cold tongue, the 1024 class top in the west and east Pacific),
the Southern Ocean's 60–70S column by depth band, the global SST and
its drift over the phase and over its last ten years, and the ice
extent against the record. It then deletes the record (unless
KEEP_FORCING=1) and all but the last of the phase's files. With
CONVERGED set to a drift in K, an ocean-only phase whose global SST
moved less than that over its last ten years is the last: no further
cycle starts. FINAL_YEARS (0) coupled years can follow the last
ocean-only phase, handed over as before. The other variables are
KEEP (4 coupled snapshots) and OCEAN_KEEP (2 year files per phase),
SNAPSHOT_DAYS, RESTORE (the ocean-only restoring, 30 W/m²/K), OCEAN
and RADIATION for both phases, OCEAN_ONLY for the ocean-only phase's
ocean options alone (for example '{"everySteps":8}'), BATCH and
LAND_FROM passed to spinup.mjs, and YEAR_DAYS (365), which only the
tests shorten.

On the M1 Max at N=64 a coupled day takes 4.9 s and a looped
ocean-only day 1.57 s with the ocean stepping every fourth atmosphere
step (1350 s). Stepping it every eighth (OCEAN '{"everySteps":8}',
2700 s) makes them 4.05 and 0.85 s; over 30 coupled days from the
paired run's day 2190 the two cadences stayed as close as a run
restarted once stays to its uninterrupted twin (daily global Ts within
0.3 K against 0.22 K, ASR within 5 W/m² against 5.9), with currents up
to 0.88 m/s and one clamped edge on two days in both, and over three
looped five-day years alone they ended within 0.01 K in the SST and
0.01 K in the Southern Ocean column. A default cycle of 10 coupled and
100 ocean-only years then takes 21 hours at N=64, 13.6 with
OCEAN_ONLY='{"everySteps":8}', so the nine cycles take 7.8 or 5.1 days.

**Interruptions.** Both spin-up scripts turn the first SIGTERM or
SIGINT into a stop after the ocean step in progress (a later one is
logged and ignored, since the driver forwards the signal its process
group may already have had): they save where they stand and exit 0 with
a `stopped by SIGTERM` line, well inside the 30 s a preempted cloud
machine gets: on the M1 Max the stop takes 0.3 s at N=64 and 2 s at
N=128, most of it writing the checkpoint (60 and 242 MB). Inside a day
`spinup.mjs` saves `<TAG>_dayDDDD_stepSSSS.bin` (DDDD days and SSSS
steps done) with the forcing recorder's part of the day when RECORD is
set, so that the day's forcing file still covers the whole day, and
`oceanSpinup.mjs` saves `<TAG>_yearYYYY_dayDDD_stepSSSS.bin`; the next
run continues from these as from any snapshot. The ocean-only phase
also saves `<TAG>_yearYYYY_dayDDD.bin` every SNAPSHOT_DAYS (30) days of
the year. Only the newest of these in-day or in-year files is kept, and
KEEP counts only whole days or whole years; the shell drivers' day
patterns see only whole days. Every ocean-only file and in-day
checkpoint also carries the ocean's restart arrays
(`restartArrays` in `js/gpu/layeredOcean.gpu.js`: the heat and salt
contents as stored, the mixed layer's previous temperature and ice,
the freshwater not yet taken and the heat flux and capacity the surface
update is using), which the ocean's upload restores as they are rather
than rebuilding them from temperature and salinity and fitting the
columns to the bathymetry, which puts the ice skin 0.08 K and the open
water 3×10⁻⁴ K off the uninterrupted run within two N=6 days. With them
an ocean-only run stopped inside a day or at any
checkpoint and continued ends byte for byte where the uninterrupted run
does (`test/asyncSpinup.test.mjs`). A coupled run continues from the
exact step, but not bit for bit: the atmosphere's step-to-step
diagnostics (the surface wind the next step's drag uses among them)
start again from zero, and a recorded day spanning an interruption
differs from the uninterrupted one by about 0.5% in its mean fluxes.

`SYNC_CMD` is a shell command both scripts run through `/bin/sh` after
every snapshot, recorded forcing day and log update, with the file's
path as `$1`, one at a time and in their own process group so that the
stop signal does not cut an upload short (three tries 5 s apart; a file
pruned before its turn is skipped); the scripts wait for the queue
before they exit. asyncSpinup.sh runs
`RESTORE_CMD` once before anything else and does not start if it fails,
and it supervises its phases itself: a node process that exits with an
error is started again from its last file after RETRY_WAIT (30) s, up
to five failures in a row, while NaN (exit 2) stops it. It stops at
`STOP_<PREFIX>` or `STOP_<PREFIX><N>` in OUT before its next segment or
phase, and on SIGTERM or SIGINT, which it passes on to the running
phase, once that has saved. On a preemptible machine with an
S3-compatible bucket, rclone configured from the environment (no config
file; Cloudflare R2 here):

```
export RCLONE_CONFIG_STORE_TYPE=s3 RCLONE_CONFIG_STORE_PROVIDER=Cloudflare \
  RCLONE_CONFIG_STORE_ENDPOINT=https://<account>.r2.cloudflarestorage.com \
  RCLONE_CONFIG_STORE_ACCESS_KEY_ID=<key> RCLONE_CONFIG_STORE_SECRET_ACCESS_KEY=<secret>
export OUT=$HOME/runs/async64
export RESTORE_CMD='rclone copy store:gcm-runs/async64 "$OUT"'
export SYNC_CMD='rclone sync "$OUT" store:gcm-runs/async64 --exclude "*.partial" --exclude "*.lock/**"'
until N=64 scripts/asyncSpinup.sh || [ $? -eq 2 ]; do sleep 60; done >> "$OUT.driver.out" 2>&1
```

This SYNC_CMD ignores `$1` and mirrors OUT, deletions included, so the
bucket holds what the pruning leaves (at N=64 a coupled snapshot is
49 MB, an ocean-only file 58 MB and a recorded year 0.8 GB) and a new
machine's RESTORE_CMD fetches only that; a mirror must never run
against an OUT that was not restored, which the driver's order
guarantees. `rclone copyto "$1" "store:gcm-runs/async64/${1#$OUT/}"`
copies each file as it is saved instead, forcing days into their
directory, but never deletes. The last line, run from the
machine's boot script (a cloud startup script or `@reboot` in cron),
restarts the driver after a failure and ends once it finishes, is
stopped or meets NaN (exit 2); each boot after a preemption starts it
again.

**Verda spot instances.** `scripts/verdaRelaunch.sh` keeps one spot
instance alive from an always-on machine with the `verda` CLI (1.8.2)
logged in. Every POLL (120) s it lists the instances and looks at the
one called NAME: a running one is left alone and its OS volume
remembered in STATE_FILE (`$HOME/.verda-relaunch-NAME`); one starting
or stopping, or of unknown status, is waited for; an offline one is
started; and when there is
none, or it was discontinued, it creates a spot instance again with
`verda --agent vm create --kind gpu --instance-type INSTANCE_TYPE
--location LOCATION --is-spot --os <OS volume or image>
--os-volume-size OS_VOLUME_SIZE --os-volume-on-spot-discontinue
keep_detached --ssh-key SSH_KEY --hostname NAME --startup-script
STARTUP_SCRIPT --wait -o json`, on OS_VOLUME when it is set and
otherwise on the remembered OS volume, once `verda volume list` shows it
detached (at the size that list gives it), and on the image OS only
while no volume is known; without OS it creates nothing until it knows
a volume. A remembered volume that is no longer
listed stops it from creating anything until OS_VOLUME is set or the
state file removed. Every action and every change of the instance's
status goes to LOG (`$HOME/verda-relaunch-NAME.log`) with a timestamp,
and it stops at STOP_FILE (`$HOME/STOP_verda-relaunch-NAME`). The OS
volume carries the repository, node_modules and the run, so after an
eviction the new instance boots the same disk, and the startup script
only restarts the driver, which continues from the newest files (a
second driver on the same OUT finds the lock and exits).
`test/verdaRelaunch.test.mjs` drives the loop through a stand-in for the
CLI that answers with the
documented JSON fields (`id`, `hostname`, `status`, `os_volume_id` or
`volumes[].is_os_volume`); it has not yet run against the service.

Once: `verda auth login`; register the key with
`verda --agent ssh-key add --name gcm --public-key "$(cat ~/.ssh/id_ed25519.pub)" -o json`
and keep its `id`; write the startup script and register it with
`verda --agent startup-script add --name gcm-async --file startup.sh -o json`,
keeping that `id` too:

```
#!/bin/bash
export REPO=/root/geodesic OUT=/root/runs/async64 N=64 PATH=/usr/local/bin:$PATH
[ -d "$REPO" ] || exit 0
mkdir -p "$OUT" && cd "$REPO" && nohup sh -c 'until scripts/asyncSpinup.sh || [ $? -eq 2 ]; do sleep 60; done' >> "$OUT.driver.out" 2>&1 &
```

(the paths are wherever the repository and the run live on the OS
volume, and PATH must reach node; Verda runs the script as root when
it creates the instance.)

Pick INSTANCE_TYPE and LOCATION from `verda --agent vm availability --spot -o json`
and the image from `verda --agent images -o json`, then start the relauncher:

```
NAME=gcm64 INSTANCE_TYPE=<type> LOCATION=FIN-01 OS=<image> SSH_KEY=<key id> \
  STARTUP_SCRIPT=<script id> nohup scripts/verdaRelaunch.sh > /dev/null 2>&1 &
```

The first instance comes from the image and the startup script finds no
repository there; log in (`verda ssh gcm64`), install node, clone the
repository to REPO, `npm install`, put the starting state in OUT and run
the startup script by hand. From then on the relauncher knows the OS
volume. Should the state file be lost, the
detached OS volume's ID is in `verda --agent volume list --status detached -o json`
(it must be in LOCATION); pass it as OS_VOLUME.

**The spin-up eleven on Verda.** Three model years of the paired run
(N=64 and N=128, fresh from the atlas on bl36, every feature at its
default) on one H100 spot instance (1H100.80S.30V, FIN-02 on Oct 2),
with the Mac pulling every state back. `scripts/verdaEleven.sh` is the
instance side: `setup` installs what the CUDA image lacks (node 22 or
later, rsync, the Vulkan loader and tools, node_modules), checks that
nvidia-smi, vulkaninfo and node's webgpu adapter all see the NVIDIA GPU,
and enables the systemd unit `gcm-eleven-resume.service`, which runs
`resume` at every boot of the OS volume;
`suite` runs every test file as a concurrent runner
(`scripts/suiteReport.sh`) and writes each failing test's name and
output to `eleven_suite.txt`; `bench` (`scripts/verdaBenchmark.mjs`)
times fresh starts of two days at N=128 and five at N=64 from the
arrival of their log lines and writes `eleven_benchmark.txt`: each N's
setup, first day, steady seconds per model day, finish and state size,
and the run's projection, per N segments × (setup + first day − steady
+ finish) + days × steady, plus one `compareStates.mjs` per round, the
hours times PRICE (the spot price, 1.85 $/h by default) and the GB of
states at the end (the setup is a fresh start's, which on the Mac at
N=64 took 7 s against 2 s for a segment continuing from a state, so the
projection errs high); `bootstrap` does the three; `run` starts
`NS="64 128" PREFIX=eleven LEVELS=bl36 PER_YEAR=36 KEEP=1000
OCEAN='{"everySteps":8}' STRATOSPHERE=1 UNTIL=1095 SYNC_CMD='sync "$1"'
scripts/pairedSpinup.sh` with OUT `/root/runs/eleven` on the OS volume,
after checking the volume has room for every state still to come, and
writes `ENDED_eleven` only when pairedSpinup logs its own stop, so a
run killed by an eviction's shutdown stays resumable; the SYNC_CMD
fsyncs every saved state and checkpoint, which spinup.mjs renames into
place without; `resume` is what the startup script
`scripts/verdaElevenStartup.sh` and the unit run at boot (the lock lets
one through), and continues only a run that `run` started and that has
neither ended (`ENDED_eleven`) nor been stopped (`STOP_eleven`); `status`
says where it stands. PER_YEAR 36 puts pairedSpinup's targets
int(k·365/36 + 0.5) 10 or 11 days apart, 108 segments a resolution, with
days 365 and 730, where a fresh land jumps, and 1095 among them: an
eviction loses at most one segment, and UNTIL ends the run there. KEEP
1000 keeps every state on the volume: 108 at each N, 105.7 MB each at
N=64 on bl36 (measured) and about four times that at N=128, some 57 GB
in all (a 150 GB OS volume). `scripts/verdaPush.sh` sends the checkout's
commit (a partial, sparse clone of HEAD without runs/: 53 MB, the data
files and a .git recording the commit included) and checks it there,
refusing while a run started there has not ended;
`scripts/verdaPull.sh` copies the states, logs and reports into
`runs/verda-eleven` of the main checkout (PULLED below) every ten
minutes, never deleting anything and never replacing a state it already
has, and with `--verify` compares the SHA-256 of every file on both
sides, pulling a differing one again and keeping the copy it replaces
as `<file>.replaced-<time>`. On the Mac, from the checkout of the
commit to run (a clean tree; `IP` is the instance's current address):

```
verda --agent instance-types -o json | node -e 'let s="";process.stdin.on("data",d=>s+=d).on("end",()=>console.log(JSON.parse(s).find(t=>t.instance_type==="1H100.80S.30V").spot_price))'
verda --agent cost balance -o json
verda --agent startup-script add --name gcm-eleven --file scripts/verdaElevenStartup.sh -o json
verda --agent vm create --kind gpu --instance-type 1H100.80S.30V --location FIN-02 --is-spot \
  --os 24.04.cuda13.0 --os-volume-size 150 --os-volume-on-spot-discontinue keep_detached \
  --ssh-key <key id> --hostname gcm-eleven --startup-script <script id> --wait --wait-timeout 15m -o json
NAME=gcm-eleven INSTANCE_TYPE=1H100.80S.30V LOCATION=FIN-02 OS_VOLUME_SIZE=150 SSH_KEY=<key id> \
  STARTUP_SCRIPT=<script id> nohup caffeinate -is scripts/verdaRelaunch.sh > /dev/null 2>&1 &
scripts/verdaPush.sh
IP=$(verda --agent vm list -o json | node scripts/verdaInstances.mjs address gcm-eleven | cut -d' ' -f3)
ssh -o StrictHostKeyChecking=accept-new root@$IP 'cd /root/geodesic && PRICE=<price> scripts/verdaEleven.sh bootstrap --detach'
nohup caffeinate -is scripts/verdaPull.sh > /dev/null 2>&1 &
PULLED=$(dirname "$(git rev-parse --path-format=absolute --git-common-dir)")/runs/verda-eleven
ssh root@$IP 'tail -n 5 /root/runs/eleven/eleven.stages.log'
cat $PULLED/eleven_gpu.txt $PULLED/eleven_suite.txt $PULLED/eleven_benchmark.txt
ssh root@$IP 'cd /root/geodesic && scripts/verdaEleven.sh run --detach'
ssh root@$IP 'cd /root/geodesic && scripts/verdaEleven.sh status'
scripts/verdaPull.sh --verify
touch ~/STOP_verda-relaunch-gcm-eleven $PULLED/STOP_pull
while pgrep -f scripts/verdaRelaunch.sh > /dev/null; do sleep 10; done
verda --agent vm delete <instance id> --with-volumes --yes -o json
verda --agent vm list -o json; verda --agent volume list -o json; verda --agent volume trash -o json
```

If the create's wait runs out, the instance may still be coming: look
in `vm list` before creating again. The relauncher starts without OS,
so it only ever recreates the instance on its kept volume, and it
remembers that volume within a poll of the instance running. A
recreated instance has new host keys and
perhaps a new address: take IP again and `ssh-keygen -R $IP` before the
next ssh (the puller keeps a known-hosts file per instance id and needs
nothing). The run is started (`run`) only once the
benchmark's projection has been weighed against the balance. The end is
`ENDED_eleven` in the run's directory (and pulled into
PULLED); `--verify` must then exit 0 with no problems in
`$PULLED/verify.txt`, and only after that are the relauncher
and the puller stopped and the instance deleted with its volume, once
the relauncher has exited, since it would otherwise recreate the
instance. An eviction near the end leaves the OS volume detached until
the relauncher, still running, recreates the instance on it, whose
resume then finds `ENDED_eleven` and does nothing; the verification
waits for that instance. A volume left detached with no instance (the
relauncher stopped first) is deleted, once verified, by `verda --agent
volume delete <volume id> --yes -o json`. A deleted volume waits 96
hours in `volume trash` before it is gone. The
instance lives for the bootstrap (about half an hour), the wait for the
go, the run (the benchmark's hours; some 4–5 h by the Sept 28 H100
timings scaled to bl36) and the final pull and verification (under half
an hour, the puller having kept up): its cost is about (PRICE + 0.04
$/h for the 150 GB volume) × that lifetime, plus the same rate over
any eviction gap's lost segment and restart, and a volume left detached
by an eviction bills on until it is reused or deleted.

**Report figures.** `scripts/figures/` draws four figures of a pulled
state, each a node dump (`<figure>.mjs <state.bin> <out.json>`) and a
python plot (`<figure>.py <json> <png> [title]`, matplotlib), titled
`<tag> day N (<season>)` on the model calendar: `stateMaps`, the mixed
layer's temperature, the sea ice's thickness (iced cells and their mean
snow), the land's vegetation cover v (0 bare, 1 closed by trees or
grass; the land means of v and of the trees in the title) and the
surface temperature;
`eqsection`, the equatorial Pacific's 2S–2N temperature to 300 m on the
ocean's layers; `eqpanels`, four GPU steps from the state and then the
lowest layer's wind, air temperature and sea-level pressure over the mixed
layer's current, temperature and the 1024.0 class top, 90E–70W within
15°; `mlmdeck`, the deck after four GPU steps with a host port of its
column checked against the GPU on the night side (where the GPU's
absorbed sunlight is zero, as the host's) and the deck's area and box
statistics. `snapshot.sh` runs the four under the GPU lock, writes
`<tag>_<figure>_dayNNNN.png` with its JSON and log, prints the paths, the
wall times and the deck's lines, and carries on past a failing figure
(exit 1); at N=128 the four take about 35 s on the Mac. The test
(`test/figures.test.mjs`) draws from the newest pulled eleven64 state
or FIGURES_STATE and skips without one.

```
scripts/figures/snapshot.sh $PULLED/eleven128_day0132.bin <outdir>
```

**Density-consistent interior.** After the mixed-layer exchanges, an
interior layer more than 0.01 kg/m³ from its label mixes in water from
the nearest layer lying clearly (by more than 0.01) on the other side of
the label, a fraction dt/2 days of the full correction per step. The move
is the same conservative transfer of mass, heat and salt as detrainment;
the correction is gradual and dead-banded because the curvature of the
equation of state makes a mixture denser than the linear estimate, and
correcting in full every step overshot and pulled water back the other
way. Two cases have no donor and are left alone: the shallowest water
under the mixed layer when too dense, and the deepest class, the
Mediterranean's 1028.0, when too light. Deep water formed like NADW
(1026.89) and AABW (1026.93) has classes of its own among those 0.03
apart and is restored like the rest. The GPU ocean
carries its free surface from the barotropic solve rather than
re-summing the layers each step: in single precision the sum rounded
about 4×10⁻⁶ m low in the same way every step, a steady loss of volume.

**Closure on fine meshes.** The ∇⁴ closure takes 12 hours to damp the
grid-scale wave at 120 km spacing and coarser (N ≤ 64). Below that, its
coefficient falls only in proportion to the spacing, so at N=128 the
grid-scale decay time is 1.5 hours. Both weaker choices ran away in the
first N=128 spin-ups:
- plain Δ⁴ scaling let a mixed-layer jet over the Vitória–Trindade
  seamounts reach the 5 m/s speed limit by day 23;
- Δ³ scaling let a barotropic meander, 120–180 km long, grow along f/H
  over the Southeast Indian Ridge from day 220, carrying up to 430 Sv
  across one edge.

A run with the present strength held to day 692. The barotropic
sub-steps themselves have no damping and a plain average, which is the
root cause still open below. The clamped-edge count in the diagnostics
is the number of edges where any layer sat at the speed limit after the
last ocean step, in both engines.

**Closure at token edges.** A class's ∇⁴ closure reads its token edges,
which carry the velocity of the layer above, the mixed layer's where the
class has pinched out. The patchy thermocline classes 1022.25–1024.0 hold
water in 12–43 % of the equatorial cells, and through their tokens the
closure tied them to the mixed layer at 0.5–3.1×10⁻⁵ s⁻¹ at N=128 (the
interfacial drag's r/h is 4×10⁻⁶), a source of momentum, since the
closure on a token edge itself is discarded: on the nine128 states of
days 91–365 its column total was 1.1–5.3 times the stress over
140–120W and 0.07–1.1 times over 120–100W, with the mixed layer's sign
(M23). By
default (`closureFill` 1, `closureTokens` 'interior') a token edge whose
two cells both hold a class thicker than THIN denser than the class
takes the least-squares uniform flow of the class's own neighbouring
edges (`closureVelocity`), and the token edges next to those are fitted
from the first ring; a class cut off by the sea floor keeps the velocity
above. What the closure then gives the fitted edges goes back to the
edges they were fitted from through the transpose of the fit under the
dc·dv inner product (`closureAdjoint`), so its work on the class is the
∇⁴ energy of the filled flow. Read on the thick edges alone, the closure
of the filled flow is not dissipative: on nine128_day0091 one ocean step
of it grows the fastest mode by 1.12–1.30 in every class from 1022.0 to
1026.8 (1.18–1.29 filling every token edge beside the class, M21's
`closureFill`), against 0.99 with the tokens as they are and with the
transpose (`scripts/closureStability.mjs`), and five N=128 days from that
state clamped 497–585 edges a day, 1963 layer edges, most of them in
the thermocline classes on the equator, at the limit after one. On the nine128 states of days 183,
274 and 365 the closure on the classes 1022.0–1025.5 at 180–110W follows
the mixed layer at −0.9 to 2.5×10⁻⁶ s⁻¹ (1021.5 on day 274: 19.9)
against 0.0–47.6×10⁻⁶ through the tokens, and over 140–120W and 120–100W
it takes +15.6/+3.5, −11.5/−7.1 and −26.1/−2.7×10⁻⁵ m²/s² from the
column, against the slab's flow, where through the tokens it gave
−7.8/−0.4, +4.8/+2.7 and +11.0/+0.5 with it (stress −7.1/−5.4,
−3.0/−3.0 and −2.1/−2.1; `scripts/closureCoupling.mjs`, 'extended').
`closureTokens` 'beside' is M21's fill without the transpose. On the
same states the transposed fill damps each thermocline class's own flow
near its tokens, 1022.25–1024.5 over 2S–2N 180–100W, at 6.6–95×10⁻⁶
s⁻¹, where the tokens give −4.4 to +0.2×10⁻⁶ and the first ring alone
with its transpose −15.8 to +10.6×10⁻⁶; on nine64_day0365 1.0–7.8×10⁻⁶
(`scripts/closureDrag.mjs`, M23).

**Friction across empty classes.** Interfacial drag acts between each
layer that holds water and the nearest layers holding water above and
below it; an empty class in between only relaxes to the velocity above.
With drag only between neighbouring classes, a mixed layer resting on
the empty tropical classes felt none, since each empty class moved with
it. In the N=64 spin-up the trades then drove the South Equatorial
Current to 3.4 m/s by day 540 and jets to the 5 m/s limit in the
Maritime Continent's seas from day 573. With the drag reaching the water
below, the same run's fastest current fell to 1.3 m/s, with no clamps.

**Eddy transport.** The layers carry the eddy-induced transport of Gent
and McWilliams (1990) in the interface-height form of MICOM and HYCOM,
one explicit step after the dynamics of every ocean step and before the
mixed-layer exchanges (`eddyTransport`; the `oEddyFlux` and `oEddyApply`
kernels). The water above each interior interface crosses an edge at
−κτ Δz dv/dc, Δz the change across the edge of the interface's depth
below the free surface, and each class carries the difference of the
fluxes at its top and its base, so the interfaces diffuse while every
class keeps its volume, heat and salt and every column its sum; the
mixed-layer base and the floor carry none, heat and salt go with the
donor cell's water, and the velocities are untouched. κ is
`eddyDiffusivity` (1000 m²/s) times 1/(1 + (L_d/dc)²), L_d = c/√(f² + 2βc)
with c = 2 m/s (Hallberg 2013's resolution function): 983 m²/s at 60°S,
902 at 20°S and 253 at the equator at N=64, and 935, 699 and 78 at N=128,
where κΔt/dc² is 9×10⁻⁵ and 2×10⁻⁴. The taper τ is linear over the top
200 m (`eddyTaperDepth`, on the shallower cell), over 100 m above the sill
and of the water below the interface that both columns hold, and over
5 m of interior water above it in the fuller column, so an outcropped
class can still spread under the mixed layer from the side where it
exists. A class thinner than 5 m on both sides carries that fraction of
its flux and the rest goes to the classes both columns hold, by
thickness, so a token class keeps its token exactly. A class is held to
1/nEdges of its water above the token per edge and step, the difference
spread the same way; on the day-2281 N=64 state one step holds a class at
1.7% of the edges, 0.02% of the flux, and no edge needs the last-resort
scaling of its whole column. With `eddyDiffusivity: 0` both engines are
byte-identical to the ocean without it. On the GPU it adds 1.2 ms to the
24 ms ocean step at N=64 and 4 ms to 96 ms at N=128. From that state
(June), 60 coupled days with κ = 1000 and with κ = 0 at N=64: the rms
slope of the 1026.6 and 1026.85 interfaces over 35–70S fell from 0.57 and
0.48 m/km to 0.48 and 0.38 by day 30 and 0.43 and 0.36 by day 60 against
0.58 and 0.52 without it, the grid-scale roughness going first (a 240 km
wave decays in 17 days, the 2000 km tilt across the ACC in about three
years, and the 45S–65S depth differences were 2–5 m smaller); the 60–70S
column stayed within 0.03 K of the run without it at every depth, and
nothing else moved outside the two runs' day-to-day spread (clamped
edges on 4 days against 3, currents ≤ 0.8 m/s against 1.0, Drake
Passage 61 Sv against 51–66).

**Open.** The freezing point ignores salinity; the deepest class has no restoring when light; the Kraus–Turner constants
and the 50 m minimum depth are first guesses; the barotropic mode has
no explicit filter beyond the sub-step average.

### M19 — Emergent vegetation (`js/physics/land.module.js`) — done (first tuning)

The uniform land of M16 (albedo 0.2, a 150 mm bucket everywhere) left
tropical land in a self-sustaining dry state: in every saved season of
the day-810 spin-up the Congo held 9–34 of 150 kg/m² at 34–40 °C and
India was 43–45 °C and dry in June and September, with the monsoon
rain falling offshore (Bay of Bengal 10–12 mm/day, India 0.7–1.8).
Over land the Betts–Miller scheme triggered in 30–50% of columns but
rained in 2–6%, because the whole column above the LCL was 4–6 kg/m²
short of its 70% reference humidity. Filling the tropical buckets gave
realistic monsoons at 27–28 °C for a month (Sahel 7.4, India
4.5 mm/day) before the 150 mm buckets drained back to the dry state;
filled 500 mm buckets carried the monsoons through the season.

So the land now grows a vegetation cover v from 0 (bare) to 1 (dense
forest) per cell, and the land's properties follow it rather than a
map: the bare-ground albedo runs from `bareAlbedo` 0.30 to
`vegetatedAlbedo` 0.13, while the bucket holds a fixed
`rootZoneCapacity` of 300 kg/m² whatever the cover: a soil keeps its
water capacity when its plants die, so a region that browns in a run
of dry years can regreen when the rain returns instead of shrinking
its own bucket and locking the desert in (the Sahel did exactly that
over three spin-up years when the capacity followed the cover).
Snow-free, v relaxes toward a goal set by how full the bucket
is — 0 below `dryWetness` 0.1 of the capacity, 1 above `wetWetness`
0.6 — over `growthTime` (180 days) when rising and `declineTime`
(365 days) when falling; under snow it fades toward 0 over
`snowDeclineTime` (720 days), so ice sheets go
bare while a boreal forest survives its winters. Water above a
shrinking bucket runs off. The soil holds two stores: rain fills a
15 kg/m² surface layer first, what it cannot hold infiltrates the
bucket with a share (soil/capacity)⁴ running off, and the surface
layer seeps into the bucket over a day. Bare ground evaporates from
the surface layer alone, so a desert dries within a day of rain; the
cover transpires from the bucket through stomata, at the aerodynamic
rate times 1/(1 + r_s g_a) with r_s = 70 s/m divided by a warmth
rising from 0 at 5 °C to 1 at 15 °C, so cold or dry roots close them
(a well-watered canopy transpires at about 60% of the potential rate
instead of the bucket's 100%); growth toward its goal needs the same
warmth, decline does not. The ice sheets — Antarctica's land and
Greenland's interior above 800 m, the topography carrying no ice mask
— grow nothing and keep an albedo of 0.8 whatever lies on them. Both engines carry v (the GPU in the PH
buffer's VEG range), snapshots save it under `land.vegetation`,
regridding samples it by tile, and the page shows it as the VEG
overlay and uses it for the land colour of the Satellite view.

A fresh land surface starts at half cover with half-full buckets, so
forests and deserts both have to emerge and neither transient is
large; `LAND_FROM=<state.bin>` on `scripts/spinup.mjs` seeds the
land from a saved state instead, regridded if its N differs. A saved
state without vegetation loads fully vegetated with full buckets where
it is free of snow. From a green start at the day-810
N=64 state, one year gives forest over the Congo (v 0.94, 6.1 mm/day
over the year), the Amazon (0.99, 7.1), Borneo, Europe, the eastern
United States and Siberia, savanna in the Sahel (0.49) and India
(0.75), and bare ground drying out over the Sahara (0.45 and falling),
Arabia, the Horn of Africa, the Kalahari, central Australia, central
Asia and the dry northeast of Brazil. From v = 0.5 everywhere the same
year sends tropical Africa and most of the Amazon to desert: the dry
state is still an attractor, so the start matters.

### M20 — Water vapour absorbs sunlight (`js/physics/radiation.module.js`) — done (tuning)

With emergent vegetation the tropical land still dried out in the
second year of a green start: over the Sahel, India and the Congo
evaporation plus runoff exceeded rain in wet years and dry ones alike,
so the land never imported ocean moisture (all land converged
+0.1 mm/day against about +0.75 on the Earth). In the dry state a
1015 hPa high sat over the Sahara where the Earth has a heat low near
1006–1008 hPa; its northeasterlies closed the Sahel to the monsoon and
the Somali jet weakened. The Sahara's column was 10–12 K too cold
(850 hPa 15 °C against about 30) although its surface budget (266
W/m² absorbed, 91 W/m² sensible) and its boundary layer (4.3 km
deep by day) were realistic, and the global rain was 3.5 mm/day
against 2.7 observed. Both point at the shortwave: only ozone absorbed
the beam aloft (3%), where the Earth's atmosphere takes 70–80 W/m²,
most of it in water vapour, so the surface received that energy and
evaporated it.

Water vapour now absorbs the beam by the Lacis & Hansen (1974)
absorptivity A(y) = 2.9y / ((1 + 141.5y)^0.635 + 5.925y) of the water
path y (cm) the beam has crossed, pressure-scaled by √σ and lengthened
by their magnification 35/√(1224μ² + 1); each layer is heated by what
its own vapour adds to A above it, and what is left goes on to the
clouds and the surface. A tropical column with 40–50 kg/m² of vapour
takes 10–20% of an overhead beam, most of it in its lower half; dry air
takes nothing. `vaporAbsorption` scales A (0 turns it off). Both
engines carry it line by line, and the absorbed-solar diagnostics
count it.

Tuned on two-year N=64 runs from a green start (`runs/twin64_day0810.bin`,
early June). Vapour absorption alone brought the global rain from 3.5
to 2.7–2.9 mm/day and warmed the Saharan column by 6 K at 850 hPa, but
the second year still dried the Sahel (0.3 mm/day) and the Congo (3.4),
and with the surface receiving less sunlight the planetary albedo fell
from 0.34 to 0.28 and the surface warmed. The defaults are therefore
`cloudScattering` 55 (planetary albedo 0.30; 95 once cloud water
absorbs, M10), a convective reference
humidity of 0.6 (the Congo's lever: 2.9 → 5.7 mm/day on its own),
`declineTime` 365 days with `bareAlbedo` 0.30 (one bad season no
longer tips a region into the bare state), and a land drag and
exchange coefficient of 1.5e-3, as over water (India's lever: 1.3 →
3.5 mm/day). Together the second year holds the Congo at 5.1 mm/day
(v 0.95), India 3.4 (0.82), the Sahel 1.4 (0.39) and the Amazon 4.3
(0.93) with the Sahara at 0.3 (0.25); global rain 2.66 mm/day,
planetary albedo 0.30, global surface temperature 14.6 °C and still
rising about 1 K/yr from the green start, though the June net
top-of-atmosphere flux (+3.6 W/m²) is nearer balance than the previous
calibration's (+5.0). The Sahara at 850 hPa is 21 °C against the
previous 15 (observed about 30), and over July–August a heat low of
1008 hPa sits on it where the dry state had a 1015 hPa high, so the
Sahel's surface wind is a southwesterly monsoon inflow (u +1.3, v
+1.5 m/s) where it was northeasterly (−3, −2), and India's trough is
1004 hPa under southwesterlies.

### The page as a client of the model worker

The page subscribes to exactly what it draws — the level, the depth,
the fields of the current overlay, animation and contours, and the
global diagnostics only while the Model dialog is open
(`js/frames.module.js`) — and every frame the worker posts carries
exactly that. The depth chooses, per cell, the isopycnal layer holding
it for the ocean's temperature, current and animation, and the
vertical velocity through it as the divergence of the transport above
(`depthFields` in `js/ocean/layered.module.js`); the atmosphere's
vertical velocity at the level is πσ̇ from the layers' mass-flux
divergences plus σ dπ/dt, converted with the level's density
(`verticalVelocity` in `js/levels.module.js`), and it agrees with the
core's own πσ̇ to rounding; shown averaged with the neighbouring cells
and over a two-hour memory, as forecasters look at ω, so that the
gravity waves and the cell-by-cell convection do not speckle it. Both engines compute them where the state
lives, the GPU only when subscribed. The level and the depth are
a Surface button beside a slider that gives every model layer the
same width: the level runs over the σ interfaces from 1000 to 10 hPa,
in ln p within a layer, and the depth over the ocean layers' nominal
bases from 1 to 5500 m; any value is accepted in
the address, and the colour ranges interpolate between the standard
levels. Where the ground rises above the level (its surface pressure
below the level's) the page masks the field and keeps the particles
out, and the viewer draws every masked cell — under the ground, below
the sea floor, or land in Ocean mode — as the grey relief of its
elevation, so the terrain shows through the slice. On the GPU engine the
frame is computed where the state lives: `frameFields` interpolates the
level fields (with the surface geopotential in the heights), the comfort
measures, the column water, cloud and sea-level pressure; `frameRain`
keeps the three-hour rain memory and the runoff tally; the diagnostics
are a workgroup reduction finished in double precision on the host
(`reductionKernel`). Only the subscribed ranges are copied back — a few
arrays of one value per cell instead of the whole state, about 140 MB
per frame at N=128. `test/frameGpu.test.mjs` checks every field and
diagnostic against double-precision sums over the downloaded state.

The page colours an overlay on the GPU: it uploads one float per cell to
a texture on the centre texture's layout and sets a colour map (palette
stops, cloud cover or flat); the vertex shader maps the value and
linearizes it. A frame costs the page 0.02 ms and 160 KB at N=64
(0.14 ms and 640 KB at N=128) against 0.9 ms and 720 KB (3.7 ms and
2.9 MB) for colouring every vertex on the main thread.

The cloud overlays have a row of their own: All clouds, the `cloud`
field the Satellite view draws (white over grey at opacity
1 − exp(−g / 40 g/m²), its legend to 100 g/m²), and each type on its
own from its own frame field, Low, Mid and High cloud (`cloudLow`,
`cloudMid`, `cloudHigh`: the resolved condensate in the layers whose
midpoint pressure is above 800 hPa, between 800 and 500, 500 or less),
Cumulus (`cloudCumulus`, the plumes' cover × condensate) and
Stratocumulus deck (`cloudDeck`, the mixed-layer deck's cover × water
path); `cloud` is their sum, and `test/frameGpu.test.mjs` checks each
type against the state and against the CPU engine. Each type paints the
same curve stretched to its own range (`CLOUD_RANGES` in
`js/frames.module.js`), the opacity 1 − exp(−g / (0.4 × range)), so all
the legends share their stops; the ranges sit near the 95th percentile
of the cells holding each type on eleven64 and eleven128 at day 1825:
low 200 g/m² (175 and 221), mid 500 (486, 589), high 400 (379, 358),
cumulus 40 (35, 39) and deck 150 (150 at both, the deck's water-path
cap). Total cloud water stays a palette overlay in the rain row.

The browser draws the globe on the same GPU, so the worker must not let
its queue run ahead of the page's frames. With no wait on the device, a
frame that asked for nothing let the worker queue steps far faster than
they ran: 3.8 frames per second with nothing drawn, 44 with the wind
shown. The worker keeps two steps in flight: after queuing one, it waits
for the one before to finish. That holds 60 frames per second at the
throughput of an unbounded queue. In Chrome on the M1 Max that is 391
simulated hours a minute at N=64 and 46 at N=128. Waiting after every
step gave 367 and 45.

The page counts frames that come late against the display's own frame
interval, leaving out delays that its own work explains: frame messages,
its loop, and the globe's and the particles' drawing. It reports them
each second, and the pacer (`js/pace.module.js`) tries idle pauses of 2,
4 and 8 ms between steps, keeping one only while it cuts the late frames
by a third. On an iPhone 17 Pro the late frames came as often without a
pause as with one; the earlier controller had pinned 12 ms there at 10%
of the rate. `?pace=off` stops the reports.

The model dialog's GPU profile (`js/gpu/profile.module.js`) runs 16
steps one at a time. Where the browser offers `timestamp-query`, it puts
timestamps around every compute pass and charges each pass to the
kernels it dispatched. It also measures an empty round trip and one
frame, and gives the result as text to copy from another device. On the
M1 Max at N=64 a step takes 12.7 ms: the atmosphere's RK stages 7.4 ms,
the boundary layer 1.7, the ocean about 2, and a frame 3.6 ms.

When the address names no resolution, engine or saved run, a first visit
tests the device. It times an N=64 model on the GPU, or one CPU thread
at N=16 when there is no usable GPU. It then picks the highest of
N=128, 64 and 32 (for the CPU: 64, 32 and 16), capped at N=128 on
desktop browsers and N=64 on phones and tablets, that is projected to
clear 30 simulated hours a minute (`js/deviceChoice.module.js`). The
chosen resolution starts from its own default run when there is one
(`js/defaultRun.module.js`): the paired spin-up's day-810 states at
N=128 and N=64. The choice is
kept for the same browser and GPU, and "Test again" in the model dialog
repeats the test. The M1 Max measures 14.2 ms a step at N=64 and runs
N=128.

Below 600 px wide, or on a short touch screen, the page takes a phone
layout:
- the settings panel is a bottom sheet, and the globe lifts into the gap
  above it (`viewer.setInsets`, a translation of the projection);
- overlays are listed by name with their notes, in place of hover tips;
- the legend spans the width.

On any touch screen, one finger turns the globe; two pan, pinch to zoom
about their midpoint, and roll once they twist past about 11° (the
maths is in `js/gestures.module.js`).

### M21 — Measurable subsidence and convection — done (tuning; audit of Sept 29 2026)

An audit of the day-183 state of the fresh-atlas paired run (N=128,
bl34, 44-class ocean) measured what the atmosphere does with vertical
motion and convection, against reanalysis values for September. The
mean circulation is right: the southern Hadley cell peaks at
144·10⁹ kg/s, the Walker cell rises over the warm pool and sinks over
the east Pacific at Earth's rates, and the SE Pacific box sinks at
0.050 Pa/s at 700 hPa. Per cell the diagnosed vertical motion is not
usable: 55–68 % of its variance in the lower troposphere sits at the
neighbouring-cell scale (0.2 % for π and θ), the lowest-layer
divergence is white noise, and the deck's ten-day mean subsidence is
still 85 % grid-scale, so the 30 % of SE Pacific columns that fail
the 0.3 mm/s gate fail on noise around an adequate 1.2 mm/s mean. The
noise is the divergent computational branch of the hexagonal C-grid
that §8 anticipates. Convection is the simplified Betts–Miller
relaxation: it lifts the 40 m surface layer with no inhibition, no
entrainment and no downdraft, fires shallow tops in the deck regions
(SE Pacific 1.47 mm/d against Earth's 0.1–0.3, 7 % of the box firing
per step, low cloud 0.08 against 0.6–0.7), heats and dries the
subcloud layer by 34 K/d and 45 g/kg/d in firing ITCZ columns, and
keeps its convective/large-scale split only as unreported global
sums. The ITCZ sits on the equator (zonal-mean rain peak 4.9 mm/d at
2.5S; Earth 6–7 at 8N) and the central Pacific ITCZ is missing, with
descent at 500 hPa where Earth's strongest rain band rises; the
equatorial SST is uniform at 26 °C, the eastern thermocline at 90 m,
and the surface current east of 140W 0.06 m/s against Earth's 0.2.

The work, in order:

1. Diagnostics. Per-cell daily convective and large-scale rain,
   saved in the state, printed in the spin-up log and mapped in the
   quarterly report; the audit scripts promoted into `scripts/` as a
   snapshot tool whose headline numbers (deck-box rain, firing
   fraction, inversion jump, stability, ten-day sink; ITCZ ascent and
   rain) the report prints; an equatorial ocean line (surface and
   thermocline-class zonal current, stress, mixed-layer depth,
   thermocline tilt).

   Built. Both engines sum each cell's convective rain (the
   Betts–Miller rain less its detrained share) and large-scale rain
   (autoconversion less the rain evaporated on the way down, which
   counts the detrained anvil water once it rains out), clear the sums
   with the precipitation, and at each diagnostics turn them into means
   over the interval in mm/d, `moist.convectiveRain` and
   `moist.largeScaleRain`: daily means in a spin-up, mirrored from the
   GPU on sync, saved in states as `convectiveRain` and `largeScaleRain`
   (zero in older ones) and regridded with them. Each spin-up segment
   ends with two lines from `js/audit.module.js`: `convection after N
   days` (the convective share globally and 15S–15N; the SE Pacific
   box's rain, its convective part and the fraction of its column-days
   with any convective rain; the Pacific ITCZ box's rain) and `equator
   after N days` (2S–2N: the mixed layer's zonal current over 160E–100W
   and 140W–100W, the fastest thickness-weighted eastward current of the
   classes to 1026.0 over 180–100W and its depth, the zonal stress over
   160E–100W, the mixed-layer depth over 140W–100W, the 1024 class top at
   150E–180 and 120W–90W). `node scripts/verticalAudit.mjs <state.bin>`
   prints the audit's headline numbers from one state on the CPU with
   their Earth references and verdicts in 84 s at N=128, the window
   numbers over 8 steps; its low cloud (cloud water below 680 hPa, or
   the deck's cover) is 0.035 where the 0.08 above counted cloud in any
   layer. The baselines, day 183 of the paired run:

   | | N=128 | N=64 | Earth |
   |---|---|---|---|
   | SE Pacific rain, mm/d | 1.47 | 1.06 | 0.1–0.3 |
   | its convective share | 1.00 | 1.00 | ≤ 0.1 |
   | columns firing a step | 0.067 | 0.047 | < 0.01 |
   | deck's virtual jump, K | 1.49 | 1.90 | 6–12 |
   | estimated inversion strength, K | 2.96 | 2.72 | 5–8 |
   | saved ten-day sink, mm/s (at h, m) | 1.17 (683) | 2.52 (857) | 3–5·10⁻³ h |
   | ω700, Pa/s | 0.050 | 0.056 | 0.03–0.05 |
   | low cloud | 0.035 | 0.098 | 0.6–0.7 |
   | deck height / resolved inversion, m | 735 / 1328 | 872 / 1378 | 1000–1500 |
   | Peru rain, mm/d (firing) | 0.59 (0.015) | 0.81 (0.025) | 0.1–0.3 |
   | Pacific ITCZ 5–12N rain, mm/d | 2.73 | 4.99 | 6–9 |
   | its ω500, Pa/s | +0.013 | −0.018 | −0.05 to −0.10 |
   | zonal-mean rain peak, mm/d (lat) | 4.94 (2.5S) | 5.52 (4.5N) | 6–7 (8N) |
   | ω700 grid-scale share | 0.665 | 0.548 | < 0.065 |
   | Hadley peaks, 10⁹ kg/s S / N | −144 / 36 | −152 / 34 | 100–200 / 10–50 |

   The spin-up's daily ASR, OLR and albedo are day means. Both engines
   sum each cell's absorbed, atmosphere-absorbed, incoming and reflected
   sunlight and its outgoing longwave over the steps since the last
   read-out and clear the sums with the rain's; the diagnostics give the
   area-weighted sums over the step count and the albedo as the summed
   reflected over the summed incoming sunlight, the last step's values
   under `instantaneous`, and the per-cell means
   (`radiation.meanAbsorbedSolar`, `meanOutgoingLongwave`,
   `meanPlanetaryAlbedo`) are mirrored from the GPU on sync and saved in
   states (zero in older ones; `test/dayMeans.test.mjs`). The last step
   alone sees the sun at the same UTC hour every day: from the N=64
   day-183 state the global albedo runs 0.220–0.296 through day 184,
   whose last step reads 0.287 against the day's 0.251 (day 185: 0.312
   against 0.280), so runs with different cloud geography under that
   hour differ by up to 0.04 of albedo and 15 W/m² in a single step's
   values.

   The cloud-radiative effects. With the radiation's `clearSkyPass`
   (off by default; on in `scripts/spinup.mjs` and
   `scripts/verticalAudit.mjs`) each column also finds the clear-sky
   absorbed sunlight and outgoing longwave of the same column with no
   resolved cloud, deck or cumulus: the clear two-stream with the ozone
   and vapour absorption kept, and the vapour and gas bands' upward pass
   with their clear emissivities plus the open window, which is the
   column's own OLR with the cloud terms zeroed. Both engines sum them
   per cell beside the day-mean sums, and the diagnostics give
   `clearAbsorbedSolar`, `clearOutgoingLongwave`, the shortwave effect
   `shortwaveCloudEffect` (ASR less clear-sky ASR) and the longwave
   effect `longwaveCloudEffect` (clear-sky OLR less OLR) as means over
   the same steps, the per-cell effects in
   `radiation.meanShortwaveCloudEffect` and `meanLongwaveCloudEffect`
   (mirrored from the GPU, saved in states). The spin-up's daily line
   gives them after the albedo (`SWCRE x LWCRE y`); the audit gives
   them globally (Earth −47 ± 4 and +26 ± 3 W/m², CERES EBAF) and over
   30S–30N, as the saved day means of the state's last day or else over
   its window. The pass runs every step: in `js/gpu/profile.module.js`
   over 128 steps from ten64_day0183 and ten128_day0183, two runs each,
   the physics pass took 2.24 and 2.23 ms without it and 2.28 and
   2.28 ms with it at N=64 (+0.04 ms of a 26 ms GPU step, 0.15 %), and
   9.52 and 9.65 against 9.70 and 9.69 ms at N=128 (+0.11 ms of
   115 ms, 0.1 %), the step medians unchanged within their run-to-run
   spread (21.0–21.7 and 90.3–92.0 ms). A column's clear-sky fluxes equal
   those of the same column without its cloud exactly, and over 48 steps
   at N=6 the engines' per-cell clear-sky sums agree to rms 1.1·10⁻⁵ and
   their read-outs to 0.002 W/m² (`test/cloudEffect.test.mjs`).

   Clear-sky scattering. With the pass on, the three-day N=64 run from
   eight64_day0183 gave a day-186 clear-sky albedo of 0.094 (CERES EBAF
   about 0.155, 53 W/m² reflected): the clear column reflected only at the
   surface. The clear two-stream now has Rayleigh and aerosol scattering in
   the visible half of the beam (`rayleighDepth`, `landAerosol`,
   `seaAerosol`; the header of `js/physics/radiation.module.js`), the clear-sky
   pass with it, and `skylight`, the diffuse fraction that stood in for
   that scattering, is 0. The budget of the day-186 sun on the day-183
   state, 48 instants, W/m² of 340.5 incoming (`scripts/clearSkyBudget.mjs`):

   | | before | after | Earth |
   |---|---|---|---|
   | atmosphere over a black surface | 0.0 | 23.7 (Rayleigh alone 21.9, aerosol alone 3.1) | Rayleigh 15–20, aerosol 3–5 |
   | surface, seen from the top | 32.1 | 26.5 | about 30 |
   | clear-sky reflected | 32.1 | 50.3 | 53 |
   | ozone / vapour / aerosol absorbed | 10.2 / 46.5 / 0 | 10.2 / 46.5 / 1.1 | |
   | absorbed at the surface | 251.6 | 232.4 | |
   | clear-sky albedo, global / 30S–30N | 0.094 / 0.074 | 0.148 / 0.122 | 0.15 ± 0.01 |
   | clear-sky albedo over open sea / land / sea ice / ice sheets | 0.047 / 0.174 / 0.640 / 0.709 | 0.105 / 0.219 / 0.651 / 0.723 | |

   The surface albedos, insolation-weighted (open sea 0.059, land 0.208,
   sea ice 0.622, ice sheets 0.800), are unchanged. Three-day N=64 GPU
   runs, day 186: `rayleighDepth` 0.16, 0.18, 0.20 give 0.143, 0.147 and
   0.152; aerosol depths 0.12 and 0.08 at 0.18 give 0.148. Defaults:
   `rayleighDepth` 0.18, `visibleFraction` 0.5, `landAerosol` 0.11,
   `seaAerosol` 0.06, `aerosolAlbedo` 0.95, `aerosolAsymmetry` 0.7,
   `aerosolHeight` 2000 m, `skylight` 0. Day 186 before and after: albedo
   0.279 and 0.312, ASR 245.6 and 234.3, OLR 242.3 and 242.1, SWCRE −62.9
   and −56.0, LWCRE 17.8 and 17.6 W/m², rain 1.74 and 1.70 mm/d (the audit's
   2.06 and 2.04), sunlight absorbed at the sea's surface 192.3 and
   178.9 W/m². A fresh atlas start on bl34, days 6–10: albedo 0.438 and
   0.459, ASR 191.1 and 184.2, OLR 220.2 and 220.5, SWCRE −115.4 and
   −104.3, LWCRE 37.0 and 36.4 W/m², rain 4.56 and 4.53 mm/d, clear-sky
   albedo 0.100 and 0.153, the sea's surface 116.8 and 108.8 W/m². With
   `rayleighDepth` 0, both aerosol depths 0 and `skylight` 0.15 the CPU
   engine keeps its previous digests bit for bit and the GPU's day-186
   line from eight64_day0183 is the previous one to every printed digit.
   Over a black surface the clear column reflects τ/(τ + 2μ) of the
   visible beam less ozone and aerosol absorption to 9.1e-15 at three
   zenith angles; in the 721 sunlit columns of an N=12 real-geography
   state absorbed plus reflected equals the beam to 4.0e-16 and the
   layers take the atmosphere's absorption to 4.1e-15; on a random set of
   sunlit columns the engines' fluxes differ by 1.1e-5 of the beam
   (`test/clearScattering.test.mjs`), and physics alone at N=6 the
   scattering changes a layer's heating by up to 0.053 K/day with the
   engines' change apart by 5.4e-6 K/day (`test/gpuModel.test.mjs`). On
   eight64_day0183 at N=64 every sunlit column closes to 4.1e-16 of the
   beam on the CPU at four instants and to 2.5e-7 (single precision) on
   the GPU over 64 steps, no dark column carries any shortwave on either,
   and the scattering leaves the longwave unchanged. On nine64_day0091 and
   nine64_day0365 (Oct 1 review, the final defaults) the columns close
   alike (CPU 4.1·10⁻¹⁶ and 4.2·10⁻¹⁶, GPU 2.4·10⁻⁷ and 2.4·10⁻⁷ over 64
   steps), the surface absorbs (1 − albedo) of its direct and diffuse
   light to 6.5·10⁻¹⁶ of the beam on the CPU and over land to 1.9·10⁻⁷ on
   the GPU, and the clear-sky OLR is unchanged, but the all-sky OLR of
   934 and 861 sunlit deck columns moves by up to 0.087 W/m² on both
   engines: the mixed-layer deck's cover follows the sunlight its layer
   absorbs, which the scattering above it changes. The daily line gives the clear-sky
   reflectance and the sea's surface sunlight; the audit and the second
   sweep's `clearAlbedo` term (target 0.15, tolerance 0.01, weight 3) read
   the clear-sky albedo.

   Acceptance by surface class (Oct 1). The global clear-sky albedo is an
   outcome, reported without a verdict: the planet's mix of surfaces
   emerges, so each class is checked against its own reference. The
   second sweep's `clearAlbedo` term in `scripts/sweep/score.mjs` still
   scores the global value against 0.15. Four
   changes, each on both engines with parity tests.

   - The Rayleigh band from the spectrum (`scripts/rayleighReference.mjs`):
     the visible band as 0.297–0.711 µm of a 5778 K Planck spectrum (the
     0.5 of the beam below 0.711 µm less ozone's 0.03 below 0.297 µm), 40
     wavelengths, τ_R(λ) of Hansen & Travis (1974) at 1013.25 hPa (0.0973
     at 0.55 µm, band mean 0.229), each through the model's τ/(τ + 2μ);
     beside it scalar doubling-adding with the azimuth-averaged Rayleigh
     phase function. The grey 0.18 departs from the reference by −11.8 % at
     μ = 1 and +10.4 % at μ = 0.1 (minimax grey 0.182: 11.0 %); two
     sub-bands, weights 0.712 / 0.288 at depths 0.0874 / 0.5687
     (`rayleighBands`), follow it within 0.1 % over μ 0.1–1 and 0.8 % at
     0.05. The grazing-sun excess that remains is the two-stream's own:
     against doubling-adding it is +17.7 % at μ = 0.05, +9.6 % at 0.1,
     +3.8 % at 0.2 and −0.9 % at 1.

     | μ | 0.05 | 0.1 | 0.2 | 0.5 | 1 |
     |---|---|---|---|---|---|
     | spectral two-stream (the reference) | 0.573 | 0.429 | 0.297 | 0.162 | 0.094 |
     | doubling-adding | 0.487 | 0.392 | 0.286 | 0.162 | 0.095 |
     | grey 0.18 | 0.643 | 0.474 | 0.310 | 0.153 | 0.083 |
     | two sub-bands | 0.577 | 0.430 | 0.297 | 0.162 | 0.094 |

     Rayleigh-only reflection of a full atmosphere over a black surface,
     global mean: 23.5 W/m² by the reference, 23.4 by doubling-adding,
     22.3 under the grey 0.18, 23.6 under the two sub-bands; on
     eight64_day0183 lit over day 186 (the terrain's lower pressures) 21.9
     under the grey depth and 23.1 under the sub-bands.
   - Aerosol from natural backgrounds, mid-visible: `seaAerosol` 0.07, the
     remote ocean's 0.06–0.07 at 500 nm (Smirnov et al. 2009, the Maritime
     Aerosol Network; tropical Pacific mean 0.07, mode 0.06, Smirnov et al.
     2003); `landAerosol` 0.12, a remote continental background of about
     0.05 (AERONET's clean and aged-background sites 0.04–0.1, Eck et al.
     2009) plus the land-mean dust, 0.068 by CALIOP and 0.103 by MODIS
     (Song et al. 2021; global dust 0.030 ± 0.005, Ridley et al. 2016),
     spread over all land. Area-weighted they give about 0.085, against the
     natural 0.09 of MAC-v1 (0.13 in all, 0.037 of it anthropogenic; Kinne
     et al. 2013). Single-scattering albedo 0.95 and asymmetry 0.7 kept
     (from memory: dust 0.93–0.97, sea salt near 0.99). Aerosol alone over
     a black surface: 3.1 → 3.5 W/m².
   - The surface's reflection absorbed on its way up (`upwardAbsorption`):
     the vapour's share by the Lacis–Hansen absorptivity of the path down
     plus 5/3 of the column's, the aerosol's 1 − exp(−(1 − ω) τ_a 5/3)
     (`test/clearScattering.test.mjs`: to 1.5·10⁻¹⁴ and 1.6·10⁻¹⁵ against
     the formulas, engines' change apart by 5.1·10⁻⁶ of the beam).
   - Wet bare soil (`soilDarkening`, land): the bare soil's albedo falls
     linearly with the bucket's fill from `bareAlbedo` 0.30 at a fill of
     0.2, below which the cover settles under 0.2, to `wetSoilAlbedo` 0.15
     at 0.5 (`darkeningWetness` [0.2, 0.5]); the engines agree to 4.4·10⁻⁸
     (`test/landGpu.test.mjs`). The snow-free albedo by bucket fill and
     vegetation cover v:

     | fill | 0 | 0.2 | 0.3 | 0.4 | 0.5 | 1 |
     |---|---|---|---|---|---|---|
     | v = 0 | 0.300 | 0.300 | 0.250 | 0.200 | 0.150 | 0.150 |
     | v = 0.5 | 0.215 | 0.215 | 0.190 | 0.165 | 0.140 | 0.140 |
     | v = 1 | 0.130 | 0.130 | 0.130 | 0.130 | 0.130 | 0.130 |

   The global clear-sky albedo of eight64_day0183 lit over day 186
   (`scripts/clearSkyBudget.mjs`, 48 instants) moves 0.148 (the grey build)
   → 0.152 (sub-bands) → 0.152 (aerosol) → 0.149 (upward absorption) →
   0.145 (wet soil). Three N=64 GPU days from a copy of eight64_day0183
   with the final defaults, day 186, before → after: albedo 0.312 → 0.310,
   ASR 234.3 → 234.9, OLR 242.1 → 242.2, SWCRE −56.0 → −56.2, LWCRE 17.6 →
   17.6 W/m², clear-sky albedo 0.147 → 0.145, rain 1.70 → 1.71 mm/d (the
   audit's 2.04 → 2.05), the sea's surface sunlight 178.9 → 177.7 W/m²; the
   60–90N ice lost from nine64_day0091 0.166·10³ km³/day (9.191 → 8.692;
   0.174 before the scattering). The physics pass's added sub-band and
   escape streams left the three days at 0.6–0.7 wall minutes, as before.
   The day-186 state by class (`scripts/clearSkyBudget.mjs` and the audit's
   `clear sky,` rows, `scripts/clearSkyClasses.mjs`):

   | class | area | surface albedo | at the top | atmosphere's own | reference | verdict |
   |---|---|---|---|---|---|---|
   | open sea 0–30° | 0.372 | 0.044 | 0.095 | 0.066 | top 0.08–0.10 | matches |
   | open sea 30–50° | 0.194 | 0.058 | 0.117 | 0.077 | top 0.10–0.13 | matches |
   | open sea 50–70° | 0.098 | 0.090 | 0.159 | 0.097 | top 0.13–0.20 | matches |
   | open sea 70–90° | 0.019 | 0.184 | 0.275 | 0.161 | | |
   | partly vegetated (v 0.2–0.7) | 0.191 | 0.199 | 0.208 | 0.074 | 0.18–0.25 | matches |
   | dense vegetation (v > 0.7) | 0.056 | 0.137 | 0.166 | 0.075 | 0.12–0.15 | matches |
   | thin snow on land (< 10 kg/m²) | 0.009 | 0.186 | 0.241 | 0.119 | | |
   | snow on open land (v < 0.5) | 0.001 | 0.542 | 0.499 | 0.147 | 0.60–0.85 | low by 0.058 |
   | snow under forest (v ≥ 0.5) | 0.003 | 0.472 | 0.437 | 0.132 | 0.20–0.35 | high by 0.122 |
   | thin sea ice (< 0.5 m) | 0.000 | 0.446 | 0.427 | 0.114 | 0.20–0.50 | matches |
   | bare sea ice | 0.000 | 0.568 | 0.520 | 0.174 | 0.50–0.60 | matches |
   | snow-covered sea ice (≥ 10 kg/m²) | 0.026 | 0.749 | 0.649 | 0.118 | 0.80–0.85 | low by 0.051 |
   | ice sheets | 0.030 | 0.800 | 0.715 | 0.132 | 0.80–0.85 | matches |

   | open sea by μ | 0–0.1 | 0.1–0.2 | 0.2–0.4 | 0.4–0.7 | 0.7–1 |
   |---|---|---|---|---|---|
   | mean μ; share of the sea's sunlight | 0.067; 0.008 | 0.156; 0.025 | 0.312; 0.103 | 0.566; 0.331 | 0.860; 0.532 |
   | direct-beam albedo (Briegleb et al. 1986) | 0.343 | 0.246 | 0.136 | 0.059 | 0.027 |
   | Taylor et al. (1996); Fresnel, flat water | 0.212; 0.660 | 0.160; 0.387 | 0.101; 0.158 | 0.057; 0.045 | 0.036; 0.022 |
   | surface albedo, direct and diffuse | 0.262 | 0.206 | 0.124 | 0.058 | 0.029 |
   | at the top | 0.373 | 0.290 | 0.194 | 0.115 | 0.074 |

   Global and by type: clear-sky albedo 0.145 (open sea 0.108 at the top
   over a surface of 0.053, land 0.201 over 0.187, sea ice 0.646 over
   0.744, ice sheets 0.715 over 0.800); the atmosphere over a black surface
   reflects 25.3 W/m², 0.074 of the incoming. The eight64 state has no bare
   land: of its ice-free land 0.782 carries a cover of 0.2–0.7 and 0.218
   more, after 186 days from the 0.5 every cell starts with. The year-six
   state five64_day2190 (0.121 of its ice-free land under 0.2) lit over its
   own day: bare dry soil (fill < 0.35, 0.029 of the globe) 0.275, low by
   0.025 against 0.30–0.40; bare wet soil 0.191 (0.10–0.20, matches, a
   trace of area); partly vegetated 0.193 and dense 0.138 (both match);
   snow on open land 0.547 and under forest 0.546 (0.025 and 0.033 of the
   globe); its global clear-sky albedo 0.156 (0.159 without the darkening).
   Where a class misses, the surface model lacks it: one snow albedo 0.55
   with no ageing and no masking by forest (snow on open land low by
   0.05–0.06, under forest high by 0.12–0.20); snow on sea ice 0.75 with
   no ageing (low by 0.051); one dry soil albedo 0.30, no sand or soil
   colour, and a cover that blends bare soil into forest with no grass
   (bare dry soil low by 0.025). References: the open sea at the top,
   CERES EBAF clear-sky ocean (approximate); the surfaces, textbook ranges
   (approximate; thin ice from memory); the sea's curves as named.
   Defaults: `rayleighBands` [[0.712, 0.0874], [0.288, 0.5687]],
   `rayleighDepth` null, `upwardAbsorption` true, `landAerosol` 0.12,
   `seaAerosol` 0.07, `aerosolAlbedo` 0.95, `aerosolAsymmetry` 0.7,
   `visibleFraction` 0.5, `aerosolHeight` 2000 m, `skylight` 0; land
   `soilDarkening` true, `wetSoilAlbedo` 0.15, `darkeningWetness`
   [0.2, 0.5], `bareAlbedo` 0.30, `vegetatedAlbedo` 0.13.
   Snow and sea ice by surface (Oct 1). The model's clock starts at the
   March equinox, so the states measured here are nine64_day0091 at the
   June solstice (the Arctic pack bare and melting, no snow on it),
   eight64_day0183 at the September equinox (the Arctic minimum, the
   southern pack at its maximum under a returning sun) and nine64_day0365
   at the March equinox (the northern snow at its extent, the Arctic pack
   under 90 kg/m² of snow). Before the change land snow was 0.55 at
   20 kg/m² whatever its age or the trees above it, snow on sea ice 0.75
   and bare sea ice 0.5 whatever the temperature
   (`scripts/snowIceAlbedo.mjs`, each state lit over its own day,
   sunlight-weighted):

   - Snow-covered land on nine64_day0365: 0.066 of the globe, 0.535–0.549
     in every band from 30N to 90N; the cover under it 0.44–0.62 by band,
     0.031 of the globe at 0.3–0.5 and 0.035 at 0.5–0.7, 0.001 below 0.3
     (the atlas start's 0.5 a year on); 0.06 of it within 2 K of melting.
   - The cover through a winter: under snow it decays over
     `snowDeclineTime` 720 days, keeping 0.81–0.76 of itself over 150–200
     days of snow. The boreal belt 50–70N carries 0.648 at the September
     equinox (0.138 under snow) and 0.506 at the March one (0.861 under
     snow); in the year-six state five64_day2190 (also a March equinox)
     0.564, with 0.60 under the snow at 50–60N, 0.46 at 60–70N and 0.30 at
     70–90N.
   - Sea ice: the March Arctic pack 0.744 (snow ≥ 20 kg/m² below −10 °C,
     0.019 of the globe, 0.750), the southern pack in September 0.746, the
     June Arctic pack 0.500 (bare, its skin within 1 K of melting).
   - The snowfall against the ageing: on the March Arctic pack below
     −10 °C 0.66 mm/d (skin below freezing, the state's last day) balances
     the cold ageing at 0.82 when it slows in the cold and 0.73 at the
     plain 0.008 a day; the Antarctic pack 0.40 mm/d, 0.79 and 0.71; land
     at 60–70N 1.63 mm/d, 0.83 and 0.80.

   The scheme, on both engines (`js/physics/ice.module.js`,
   `js/physics/land.module.js`, `js/gpu/physics.gpu.js`):

   - Snow ageing (`agedSnowAlbedo`, `refreshedSnowAlbedo`): the Douville et
     al. (1995) scheme of the ECMWF land surface (Dutra et al. 2010,
     appendix eq. A7 and eq. 9): fresh 0.85; snow whose skin is within
     2 K of melting (Dutra's revised test) relaxes toward its floor at
     0.24 a day; colder snow loses 0.008 a day down to the floor; a fall of
     F kg/m² moves it min(1, F/10) of the way back to 0.85. The cold
     ageing is slowed by the temperature dependence of grain growth in
     BATS (Dickinson et al. 1993), exp(5000 (1/273.15 − 1/T)): 0.50 at
     −10 °C, 0.24 at −20 °C (from memory). BATS sums that grain-growth
     term with a melt term (its tenth power) and a dirt term 0.3, 2.3 at
     melting (as CLM carries it, from memory); taken whole and scaled to 1
     at melting it would age cold snow at 0.35, 0.23 and 0.18 of the full
     pace at −10, −20 and −30 °C where the grain-growth term alone gives
     0.50, 0.24 and 0.10. Without it the snowfall above
     holds cold snow on the Arctic pack at 0.73 against 0.80–0.85. Each
     cell's albedo is in the land's `snowAlbedo` (its ocean cells the sea
     ice's); ground and ice without snow hold 0.85 so the next snow
     starts fresh. Floors: 0.50 on land (Dutra's), 0.70 on sea ice
     (melting snow on Arctic ice about 0.7, Perovich et al. 2002, and
     CCSM3's melting snow about 0.72 in the broadband: visible and
     near-infrared cold snow 0.98 and 0.70, Briegleb et al. 2004, its
     melting reductions from memory).
   - Masking by trees: a full snow cover's albedo falls linearly from the
     snow's own to `forestSnowAlbedo` 0.27 as the standing cover rises to
     `closedCanopy` 0.7, constant above; 0.27 is MODIS's snow-covered
     evergreen needleleaf forest (Moody et al. 2007, tabulated in Dutra et
     al. 2010: deciduous needleleaf 0.33, deciduous broadleaf 0.31, mixed
     0.29), the knee the MODIS finding that snow changes the albedo little
     above about 70 % tree cover (from memory, no paper named). The
     cover ramp `fullSnow` 20 kg/m² is
     kept. The standing cover `canopy` rises with the cover at once and
     falls toward it over `canopyMemory` 365 days, so it keeps 0.94 of an
     autumn cover through 200 days of snow while the cover keeps 0.76, and
     follows a cell under snow for years down at the cover's own pace.
   - Bare sea ice: 0.62 while its skin is colder than −1 °C, falling
     linearly to 0.48 at the melting point, as CCSM3 lowers its bare ice
     over the last kelvin (Briegleb et al. 2004: cold visible and
     near-infrared 0.78 and 0.36, 0.57 broadband; melting about 0.50,
     from memory) and between SHEBA's cold bare ice 0.60–0.65 and its
     ponded July ice 0.45–0.55 (Perovich et al. 2002, as given); the
     ramp from the water's albedo below 0.5 m kept.
   - Ice sheets keep 0.80. The ageing would hold their dry interiors at
     0.80–0.85 (Antarctica below −10 °C, 0.58 mm/d, balances at 0.83;
     0.71 without the cold slowing), but at the June solstice 0.0015 of
     the globe, 44 % of the northern sheet, lies within 2 K of melting,
     where the field falls from 0.85 to 0.73 in ten days toward 0.50, and
     a margin cell bare of snow (3 kg/m² at the September equinox) would
     need a glacier-ice albedo the model does not have.

   Saved states and the page's snapshots carry `snowAlbedo` and `canopy`
   (both engines; on the GPU at the end of PH); an older state loads its
   snow at 0.85 and its standing cover at its cover. Unit tests against
   hand-computed values (`test/snowAlbedo.test.mjs`): a day of snow at
   −20 °C from 0.85 to 0.84812 (0.842 without the slowing), a wet day
   0.85 → 0.77532 (0.81799 on ice), 5 kg/m² on 0.6 → 0.725, the cover
   0.8 → 0.69626 and the standing cover 0.8 → 0.77514 after 100 days
   under snow. At N=6 on random snow, snow albedo, standing cover and
   skin temperature the engines' albedos agree to 4.9·10⁻⁸ over land and
   2.7·10⁻⁶ over 99 iced cells (22 within 1 K of melting, the single
   precision of the skin); over 48 steps with snow falling on 42 of 57
   snow cells, wet and cold, the snow albedo agrees to rms 4.8·10⁻⁶ (max
   1.9·10⁻⁵) and the standing cover to 2.8·10⁻⁷. Digests re-pinned to
   b892e42f, e228ab4c and d8b73e96; the scattering's layer-heating parity
   is held on the grey ice of 0.5 it was measured with.

   Short N=64 GPU runs from copies of the three states with the final
   defaults (OCEAN `{"everySteps":8}`, bl34), before → after. By class
   (`scripts/clearSkyBudget.mjs`, surface albedo; area share of the
   globe):

   | class | June solstice + 3 d | September equinox + 3 d | March equinox + 10 d | reference |
   |---|---|---|---|---|
   | snow among sparse trees (standing 0.2–0.7) | trace: 0.528 → 0.433 | 0.003: 0.482 → 0.335 | 0.063: 0.546 → 0.416 | |
   | snow under forest (standing ≥ 0.7) | | trace: 0.461 → 0.205 | trace: 0.487 → 0.245 (3 d) | 0.20–0.35 |
   | thin snow on land (< 10 kg/m²) | trace: 0.162 → 0.154 | 0.009: 0.186 → 0.160 | 0.004: 0.232 → 0.196 | |
   | bare sea ice, melting (T ≥ −1 °C) | 0.012: 0.500 → 0.487 | | trace: 0.500 → 0.480 | 0.45–0.55 |
   | bare sea ice, cold | trace: 0.500 → 0.620 | | trace: 0.500 → 0.620 | 0.60–0.65 |
   | sea ice under 1–10 kg/m² of snow | 0.001: 0.548 → 0.602 | 0.001: 0.491 → 0.572 | 0.001: 0.533 → 0.616 | |
   | snow-covered sea ice, cold (T < −2 °C) | 0.023: 0.738 → 0.827 | 0.024: 0.749 → 0.836 | 0.032: 0.748 → 0.828 | 0.80–0.85 |
   | snow-covered sea ice, wet (T ≥ −2 °C) | trace: 0.747 → 0.837 | 0.002: 0.747 → 0.818 | 0.002: 0.749 → 0.811 | 0.65–0.75 |
   | ice sheets | 0.030: 0.800 | 0.030: 0.800 | 0.030: 0.800 | 0.80–0.85 |

   The open sea and the snow-free land change by at most 0.001.
   No state has snow on open land (standing below 0.2): the atlas start's
   0.5 stands everywhere. The year-six state lit over its own day with
   its snow fresh: cold snow on open land 0.544 → 0.719 (0.002 of the
   globe; a standing cover of 0.1–0.2 already masks it by 0.08–0.17),
   under sparse trees 0.547 → 0.410, under forest 0.539 → 0.266, cold
   snow on sea ice 0.749 → 0.849.

   - The northern snow over ten days from the March equinox (day 375,
     sunlight-weighted): before 0.538–0.550 in every band; after 0.387
     (30–40N), 0.417 (40–50N), 0.411 (50–60N), 0.421 (60–70N) and 0.471
     (70–90N), the snow's own albedo 0.76–0.83; the snow line (the lowest
     5° band at least half covered) 50–55N both, the snow-covered land
     0.0629 → 0.0617 of the globe (45–50N 0.39 → 0.37 covered, 40–45N
     0.23 → 0.20). Ten days from the June solstice leave snow only at
     70–90N (0.04–0.05 of that land): 0.550 → 0.382, the snow at 0.57.
   - The Arctic from the June solstice: the 60–90N ice lost 0.166 → 0.168
     ·10³ km³ a day over three days (9.191 → 8.692 and 8.687), 0.192 →
     0.191 over ten (7.268 and 7.277); the pack's ice albedo 0.500 →
     0.491 (day 94) and 0.500 → 0.484 (day 101); the sunlight absorbed at
     the surface of the iced cells poleward of 60N 62.6 → 64.0 and 72.1 →
     76.8 W/m². The September Arctic pack 0.705 → 0.780.
   - The southern pack at the September equinox (day 186): 0.746 → 0.831,
     its iced cells' surface sunlight 12.9 → 10.9 W/m². The March Arctic
     pack 0.743 → 0.836 (day 368) and 0.823 (day 375), its iced cells'
     surface sunlight 5.9 → 5.3 W/m² on day 375.
   - Day means (outcomes), before → after: day 94 albedo 0.322 → 0.321,
     ASR 230.8 → 231.1, OLR 240.7 → 240.7, SWCRE −60.0 → −59.8, LWCRE
     20.2 → 20.2 W/m², clear-sky reflectance 0.1459 → 0.1457; day 186
     0.310 → 0.311, 234.9 → 234.8, 242.2 → 242.2, −56.2 → −56.1, 17.6 →
     17.6, 0.1451 → 0.1458; day 375 0.342 → 0.340, 224.1 → 224.6, 233.3 →
     233.4, −62.6 → −63.6, 21.8 → 21.7, 0.1579 → 0.1536. The clear-sky
     albedo of the day-375 state lit over its day 0.158 → 0.153 (land
     0.244 → 0.225 at the top, sea ice 0.647 → 0.706). A three-day run
     took 0.3 wall minutes and a ten-day one 0.8, before and after.

   Review (Oct 1). A cold week at −15 °C, a melting week at −0.5 °C and
   4 kg/m² of snow, in closed form, take fresh snow to 0.830669, 0.561628
   and 0.676977 on land and 0.830669, 0.724353 and 0.774612 on sea ice;
   stepped every half hour the land module gives the land's to 2·10⁻¹⁴
   and the GPU's own WGSL ageing and refresh functions both to 5.2·10⁻⁶
   (single precision), and the sea-ice module, which ages at the skin it
   ends each step with, falls 5.2·10⁻⁵ short where conduction warms that
   skin. A cell whose snow melts away and is snowed on starts at 0.85. In
   the day-94 and day-186 states of reruns of the two three-day runs
   (which reproduce the day means above and the Arctic loss of 0.168)
   every land snow albedo lies in 0.670–0.850 and every sea-ice one in
   0.783–0.850, every snow-free cell holds 0.85, and the standing cover is
   never below the cover. The ice albedo falls monotonically with the
   skin temperature, by at most 1.4·10⁻⁵ over 10⁻⁴ K through the melting
   ramp, and the masking monotonically with the standing cover. With
   grey ice (`iceAlbedo` and `meltingIceAlbedo` 0.5, `snowAgeing` false)
   the three physics digests are those pinned before the change
   (3ca002d1, 3d0c610f, da3ea94c).

   What still misses:

   - One cover cannot tell a forest from tundra or grassland. The masking
     reads the cover as tree cover, so the atlas start's 0.5 masks all
     northern snow 0.71 of the way to 0.27: 70–90N reads 0.471 against
     tundra's 0.60–0.80 (MODIS grassland and tundra, approximate), and
     50–60N reads 0.411 where Earth's snow-covered boreal forest reads
     0.27–0.33 and its steppe about 0.6–0.7 (approximate). Without the
     masking the forest snow read 0.46–0.55 against 0.20–0.35.
   - The standing cover's memory needs a winter to show: runs from older
     states start it at the cover. From the boreal belt's cover of 0.648
     in September it would hold about 0.61 by March, where the cover is
     0.506: a masked albedo of about 0.34 instead of 0.41.
   - Wet snow on sea ice reads 0.81–0.84 against melting snow's
     0.65–0.75: these cells (0.002 of the globe) are mostly at the ice
     edge, where 2–6 mm/d of snow refreshes it, and the runs start their
     snow fresh.
   - The June Arctic pack carries no snow from the solstice on, so it
     goes from bare to ponded (0.48) without SHEBA's month of melting
     snow near 0.7.
   - The ice sheets keep one albedo, 0.80.
   - Wet soil at the June solstice. Three days from nine64_day0091 (day
     94, review of Oct 1) partly vegetated land reads 0.167 against
     0.18–0.25 (0.259 of the globe; 0.205 with `soilDarkening` false),
     its snow-free buckets 0.49 full on the mean, where the September
     state's are 0.31 (0.199, inside) and the year-six state's 0.36. The
     darkening reads the root zone's fill, which stays high for months
     after rain, where BATS and CLM read the top layer's water (from
     memory). The open
     sea 0–30° reads 0.102 at the top on that day against 0.08–0.10.

   Defaults: sea ice `iceAlbedo` 0.62 (cold bare ice), `meltingIceAlbedo`
   0.48, `iceMeltingRange` 1 K, `fullAlbedoThickness` 0.5 m,
   `iceFullSnow` 20 kg/m², `snowAgeing` true, `iceSnowFloor` 0.70
   (`iceSnowAlbedo` 0.75 only without the ageing); the ageing shared by
   land and ice (`SNOW_AGEING` and `snowAgeing`; the GPU has one set for
   both and refuses a model whose land and ice options differ, a default
   counting as given) `freshSnowAlbedo` 0.85,
   `coldSnowAgeing` 0.008 /day, `meltingSnowAgeing` 0.24 /day,
   `refreshSnowfall` 10 kg/m², `wetSnowRange` 2 K, `ageingActivation`
   5000 K; land `snowAgeing` true, `oldSnowAlbedo` 0.50 (`snowAlbedo`
   0.55 only without the ageing), `fullSnow` 20 kg/m², `snowMasking`
   true, `forestSnowAlbedo` 0.27, `closedCanopy` 0.7, `canopyMemory` 365
   days, `iceSheetAlbedo` 0.80; an older state's snow albedo 0.85 and
   standing cover its cover.
   Treeline and top-soil wetness (Oct 1). Two changes on both engines
   (`js/physics/land.module.js`, `js/gpu/physics.gpu.js`), each with
   parity and hand-computed unit tests (`test/treeline.test.mjs`,
   `test/land.test.mjs`, `test/landGpu.test.mjs`).

   - A tree cover the season admits (`treeline`). The masking's `canopy`
     is a tree cover. Each land cell keeps two running means of its
     lowest layer's air temperature over `seasonMemory` (3 years):
     `seasonLength`, the share of the time at or above `seasonThreshold`
     0.9 °C, and `seasonWarmth`, the mean excess over it (K). Their ratio
     plus 0.9 °C is the growing season's mean temperature, taken over at
     least `minimumSeason` 94 days (a shorter season counts as cooler; no
     season gives 0.9 °C). The limits are TREELIM's (Paulsen & Körner
     2014, Alpine Botany 124: 1–12, 376 treelines from satellite images):
     a season of days whose mean air is at least 0.9 °C, at least 94 days
     long, at least 6.4 °C on the mean; Körner & Paulsen (2004, J.
     Biogeogr. 31: 713–732) measured 6.7 ± 0.8 °C in the root zone at 46
     treelines from 68N to 42S. The treeline factor f rises linearly from
     0 at 6.4 °C to 1 at 8.0 °C (`treelineWarmth`; the width, twice the
     2004 sites' spread, is not itself sourced). Köppen's 10 °C warmest
     month gives the same line: a sine year of −9 °C ± 19 K has a 10 °C
     warmest month and a 119-day season at 6.9 °C, a maritime 2 °C ± 8 K
     a 199-day season at 6.6 °C, where growing degree days above 5 °C
     read 283 and 448 K·d against LPJ's 350 for boreal summergreen and
     600 for boreal evergreen trees (Sitch et al. 2003, Global Change
     Biol.), which is why the season's mean and not degree days is used.
     Snow-free, the tree cover relaxes toward f v over `treeGrowthTime` 10
     years rising and `treeDeclineTime` 3 years falling; under snow it
     holds, falling toward f over 3 years only where it stands above f (a
     failed warmth, or summers that stay under snow). Times: jack pine
     reaches crown closure 20–21 years after fire with its seed in place
     (Porter et al. 2023, Sci. Rep., a chronosequence); ten years of
     snow-free time closes 90 % of the gap in 23 years where it never
     snows and about 46 where the ground is bare half the year. The
     decline time is not sourced: in the 2000–03 drought of the US
     Southwest piñon mortality reached 90 % at some sites and about 25 %
     region-wide, after 15 consecutive months of dry soil (Breshears et
     al. 2005, PNAS 102: 15144–15148), faster than three years, which
     remove 90 % in 7; the masking needs the time the dead stems stop
     shading the snow, for which no source was checked. The air is
     sampled every step, where TREELIM gates daily means.
   - Starts. A state without the season means (every state saved before
     this, and a fresh start) starts them from a sine year with mean
     −31.8 + 0.148 Q̄ − 6.5 K/km × the ground's height (°C) and amplitude
     min(0.069 ΔQ, 19.7 − 0.0205 ΔQ) (K), Q̄ and ΔQ the annual mean and
     annual harmonic of the top-of-atmosphere daily insolation at the
     cell's latitude (W/m²; `SEASON_ESTIMATE`, `seasonEstimate`), and its
     trees at the cover times f. The constants are least-squares fits to
     the model's own first year: the lowest air of nine64 and eight64 at
     days 91, 183, 274 and 365, harmonically fitted per land cell off the
     ice sheets (the amplitude's second branch over 40–85N). The atlas
     start's land temperature, 288 + 45(1/3 − sin²φ) K with no season and
     no lapse, is 3–7 K colder at 45–75N than that first year and is not
     used. Land added by regridding takes the estimate and its trees the
     guessed cover (0.5 snow-free) times the estimate's f. The estimate
     against the first years' own season (season mean °C, mean f, days):

     | band | estimate | nine64 | eight64 | nine64 / eight64 warmest month |
     |---|---|---|---|---|
     | 55–60N | 11.4, 0.99, 196 | 11.2, 0.90, 206 | 10.3, 0.90, 217 | 17.4 / 16.0 |
     | 60–65N | 9.1, 0.86, 167 | 8.8, 0.69, 160 | 9.2, 0.78, 183 | 13.0 / 13.9 |
     | 65–70N | 7.6, 0.67, 148 | 8.0, 0.66, 145 | 7.8, 0.58, 163 | 11.8 / 11.6 |
     | 70–75N | 7.0, 0.44, 142 | 6.5, 0.48, 131 | 6.9, 0.46, 145 | 9.4 / 10.0 |
     | 75–80N | 5.6, 0.03, 121 | 2.9, 0.01, 75 | 5.5, 0.15, 118 | 3.9 / 7.8 |
     | 80–85N | 4.2, 0.00, 101 | 1.9, 0.00, 40 | 5.5, 0.34, 110 | 1.9 / 8.0 |

     In the model's own climate the season's mean crosses 6.4 °C at about
     73N (nine64) and 74N (eight64) and the warmest month 10 °C at about
     71N and 72.5N, the band means' factor falling through a half between
     65–70N and 70–75N. Earth's northern treeline lies at 56N in
     Labrador (53N in parts), 61N by Hudson Bay, 68N in Alaska, 69N in
     the Northwest Territories, 70N in Norway and 73N on the Central
     Siberian Plateau (Wikipedia, "Tree line", its table of Arctic tree
     lines). The zonal band means lie at the poleward end of that range,
     and the estimate, being zonal, puts the line at 67.5–70N in
     Labrador and by Hudson Bay (2.5° bands, mean factor under a half),
     9–11° poleward of Earth's there, and at 65N in Taimyr and eastern
     Siberia, where it subtracts the plateaus' height. The model's own
     first year by sector, from per-cell fits to the four instants, is
     too noisy to place the line: nine64 and eight64 differ by up to 10°
     in a sector (Labrador 67.5N and 75N, Hudson Bay 72.5N and 65N).
   - Top-soil wetness (`soilDarkening` 'surface', the default). The bare
     soil darkens linearly with the 15 kg/m² surface layer's fill
     (`darkeningWetness` [0, 1]), from
     `bareAlbedo` 0.30 dry to `wetSoilAlbedo` 0.15 full, whatever the root
     zone holds: Idso et al. (1975, J. Appl. Meteor. 14: 109–113) found
     a loam's albedo linear in its top layer's water, 0.30 dry to 0.14 wet (a
     ratio of 0.47; 0.15 / 0.30 is 0.50), linear in a layer under 0.2 cm
     thick (as cited by later studies; the paper itself was not read).
     The 15 kg/m² store is far deeper than that, so a few millimetres of
     rain darken the model's soil by a fraction where they would darken
     the top 0.2 cm fully. The layer seeps into the root
     zone over a day, so a wetted bare soil brightens 0.150, 0.209, 0.244,
     0.265, 0.278 every 12 hours after the rain stops (with the cover
     growing over a full bucket). 'rootZone' keeps the bucket's ramp
     ([0.2, 0.5]); without the vegetated land there is no surface layer
     and no darkening. The snow-free albedo by the layer's fill and the
     cover v:

     | fill | 0 | 0.25 | 0.5 | 0.75 | 1 |
     |---|---|---|---|---|---|
     | v = 0 | 0.300 | 0.263 | 0.225 | 0.188 | 0.150 |
     | v = 0.5 | 0.215 | 0.196 | 0.177 | 0.159 | 0.140 |
     | v = 1 | 0.130 | 0.130 | 0.130 | 0.130 | 0.130 |

   Engines: at N=6 over 48 steps (season memory 6 h, trees 3 h up and 2
   h down, f over 6–22 °C; 151 land cells, 14 under snow, 113 in season,
   trees grew on 13 and died back on 128) the season length agrees to
   4.1·10⁻⁷, the season warmth to 3.1·10⁻³ K under lowest-air
   differences of 1.4·10⁻² K, the tree cover to 3.8·10⁻⁵; the albedo of
   173 land cells with a random surface layer to 7.6·10⁻⁸ (it darkens
   136, by up to 0.132; by the root zone 1.1·10⁻⁷). Hand values: a sine
   year −4 °C ± 15 K over 0.9 °C, share 0.394075 and excess 2.581744 K;
   a 146-day season at 7.9 °C, f 0.9375; a 73-day one with 1.5 K of
   excess, f 0.202793; trees 0.2 → 0.705696 in 10 years at f = 1, 0.6
   → 0.220728 in 3 at f = 0, 0.8 → 0.610364 in 3 under snow at f = 0.5.
   States save `seasonLength` and `seasonWarmth` (on the GPU at the end
   of PH, read back by `land.serialize`, which the page's snapshot uses);
   land regridding samples them by tile. The digests are unchanged (their
   model has no land).

   Runs: N=64 GPU, bl34, OCEAN `{"everySteps":8}`, from copies of the
   states, before (`treeline` false, `soilDarkening` 'rootZone', the
   previous rule) → after; the five64_day2190 state lit over its own day
   without running. Snow-covered land (`scripts/snowIceAlbedo.mjs`,
   sunlight-weighted surface albedo, area share of the globe):

   | row | June solstice + 3 d | September equinox + 3 d | March equinox + 3 d | March equinox + 10 d | five64_day2190 | reference |
   |---|---|---|---|---|---|---|
   | 70–90N (tundra) | trace: 0.450 → 0.758 | 0.0018: 0.370 → 0.714 | 0.0056: 0.473 → 0.707 | 0.0056: 0.471 → 0.704 | 0.0056: 0.589 → 0.744 | 0.60–0.80 |
   | 50–70N (boreal) | | 0.0012: 0.298 → 0.596 | 0.044: 0.424 → 0.470 | 0.044: 0.416 → 0.464 | 0.041: 0.409 → 0.449 | 0.27–0.45 |
   | 50–60N | | | 0.019: 0.423 → 0.429 | 0.018: 0.411 → 0.417 | 0.016: 0.352 → 0.359 | |
   | 60–70N | | 0.0012: 0.298 → 0.596 | 0.026: 0.425 → 0.511 | 0.025: 0.421 → 0.505 | 0.024: 0.461 → 0.531 | |
   | 30–40N | | | 0.005: 0.392 → 0.728 | 0.005: 0.387 → 0.741 | 0.004: 0.366 → 0.785 | |
   | trees < 0.1 (open; none before) | | 0.0027: 0.698 | 0.011: 0.818 | 0.011: 0.816 | 0.012: 0.839 | 0.80–0.85 cold |
   | trees 0.1–0.3 | | | 0.006: 0.645 | 0.006: 0.631 | 0.007: 0.684 | |
   | all land snow | 0.0002: 0.433 → 0.724 | 0.003: 0.333 → 0.648 | 0.066: 0.427 → 0.517 | 0.062: 0.416 → 0.509 | 0.058: 0.405 → 0.503 | |

   The tree cover by band on day 375 (cover, trees, f; share of the
   land with trees ≥ 0.2): 55–60N 0.51, 0.51, 0.99 (0.99); 60–65N 0.51,
   0.45, 0.86 (0.88); 65–70N 0.49, 0.34, 0.67 (0.74); 70–75N 0.46,
   0.21, 0.44 (0.56); 75–80N 0.40, 0.01, 0.03 (0); 80–85N 0.37, 0, 0
   (0); south of 55N and in the south the trees are within 0.02 of the
   cover. The 30–40N snow is the high plateau's, whose season's mean is
   below 6.4 °C at its height. Cold snow
   on open land (trees < 0.2) reads 0.806 and 0.802 on days 368 and 375
   (0.80–0.85, matches; 0.013 of the globe), 0.814 on the year-six
   state. Bare and vegetated land (`scripts/clearSkyBudget.mjs`, surface
   albedo; with the darkening off in brackets):

   | class | June solstice + 3 d | September equinox + 3 d | five64_day2190 | reference |
   |---|---|---|---|---|
   | partly vegetated (v 0.2–0.7) | 0.259: 0.167 → 0.200 (0.205) | 0.191: 0.199 → 0.213 (0.217) | 0.084: 0.193 → 0.222 (0.227) | 0.18–0.25 |
   | dense vegetation (v > 0.7) | trace: 0.136 → 0.178 | 0.056: 0.137 → 0.165 (0.171) | 0.083: 0.138 → 0.147 (0.149) | 0.12–0.15 |
   | bare dry soil (v < 0.2, layer < half full) | | | 0.030: 0.274 → 0.278 (0.279) | 0.30–0.40 |
   | bare wet soil (layer ≥ half full) | | | trace: 0.271 → 0.172 | 0.10–0.20 |

   The partly vegetated land's surface layer is 0.09 full on day 94
   (its root zone 0.49) and 0.06 on day 186 (0.31); the dense
   vegetation's 0.18 (0.81); the year-six deserts' 0.017 (0.13).
   Outcomes, before → after: clear-sky albedo of the state lit over its
   day 0.146 → 0.153 (day 94), 0.146 → 0.149 (186), 0.153 → 0.160 (368),
   0.153 → 0.161 (375), 0.153 → 0.159 (year six); day means (spin-up
   log) day 94 albedo 0.321 → 0.326, ASR 231.1 → 229.5, OLR 240.7 →
   240.6, SWCRE −59.8 → −58.7, clear-sky reflectance 0.1457 → 0.1535,
   Ts 16.66 → 16.57 °C; day 186 0.311 → 0.312, 234.8 → 234.2, 242.2 →
   242.2, −56.1 → −55.5, 0.1458 → 0.1493, 16.71 → 16.65; day 368 0.325
   → 0.328, 229.8 → 228.8, 233.3 → 233.3, −58.6 → −57.2, 0.1532 →
   0.1602, 13.78 → 13.69; day 375 0.340 → 0.345, 224.6 → 223.2, 233.4
   → 233.1, −63.6 → −62.4, 0.1536 → 0.1613, 13.86 → 13.70. Three days
   took 0.3 wall minutes and ten 0.8 → 1.0 (CPU work alongside).

   What still misses:

   - The rule cannot tell grassland from forest where the season is
     warm: there f is 1 and the trees are the cover. On day 375 the snow
     at 40–60N under f = 1 (0.025 of the globe) reads 0.416 under trees
     of 0.49; its part with a cover under 0.5 (0.012 of the globe, steppe
     and prairie by their moisture) reads 0.471 under trees of 0.41,
     against MODIS snow-covered grassland and cropland (Moody et al.
     2007 tabulate them by IGBP class; the values were not checked
     here). A moisture split (the cover's goal already keeps a
     dry steppe sparse) or fire would be needed.
   - The boreal belt 50–70N reads 0.464 against 0.27–0.45: its cover is
     still the atlas start's 0.5 a year on (0.50), so its trees are at
     most 0.5, masking 0.71 of the way, where MODIS snow-covered
     needleleaf and mixed forest reads 0.27–0.33 (Moody et al. 2007 as
     tabulated by Dutra et al. 2010).
   - The zonal treeline (70–75N) sits at the poleward end of Earth's
     56–73N: the model's 70–75N land has a 9.4–10.0 °C warmest month,
     and the zonal start estimate cannot hold Labrador's and Hudson
     Bay's line far south until three years of the model's own seasons
     replace it.
   - Dense vegetation in September reads 0.165 against 0.12–0.15: the
     cover blends the dry soil (0.30 at a layer 0.18 full) at 1 − v =
     0.25 into the forest, where a forest floor is shaded and littered;
     the root zone's darkening had hidden it at 0.137.
   - Bare dry soil still reads 0.278 against 0.30–0.40 (one soil colour).
   - Ten days cannot show the rule's own times: the season means move
     by about 1 % in ten days, the trees by under 1 % of their gap.

   Defaults: land `treeline` true, `seasonThreshold` 0.9 °C,
   `minimumSeason` 94 days, `treelineWarmth` [6.4, 8.0] °C,
   `seasonMemory` 3 years, `treeGrowthTime` 10 years (snow-free),
   `treeDeclineTime` 3 years, `SEASON_ESTIMATE` mean [−31.8, 0.148],
   amplitude [0.069, 19.7, −0.0205], lapse 0.0065 K/m; `soilDarkening`
   'surface', `darkeningWetness` [0, 1] ('rootZone' [0.2, 0.5]),
   `wetSoilAlbedo` 0.15, `bareAlbedo` 0.30, `surfaceCapacity` 15 kg/m²;
   the GPU takes these and `percolationTime`, `stomatalResistance`,
   `growthColdest`, `growthWarmest` and `iceSheetAlbedo` from the land's
   options (it had taken the last five from its own defaults whatever
   the options said); `canopyMemory` 365 days only without the
   treeline. Older states start their season means from the estimate
   and their trees at the cover times f.

   Trees and grass (Oct 1). Where the season admits trees, water now
   decides between forest and grass, and the cover not under trees is
   grass, on both engines (`js/physics/land.module.js`,
   `js/physics/radiation.module.js`, `js/gpu/physics.gpu.js`), with
   parity and hand-computed unit tests (`test/grassland.test.mjs`).

   - The moisture gate (`treeMoisture`). Each land cell keeps running
     means over `moistureMemory` (3 years, the season means' memory; no
     source was found for another) of its rain and snow (`rainMean`,
     mm/d) and of its evaporative demand (`demandMean`, mm/d). The
     demand is the FAO-56 Penman–Monteith reference evapotranspiration
     (Allen et al. 1998, FAO Irrigation and Drainage Paper 56) from the
     model's own fluxes: the surface's net radiation, the lowest air's
     temperature, humidity and density at surface pressure, and the
     land's aerodynamic conductance C|V| (the exchange the land's
     evaporation already uses), through a surface resistance of 70 s/m
     (`REFERENCE_RESISTANCE`, FAO-56's reference grass; the land's
     `stomatalResistance` is also 70). The treeline factor f is
     multiplied by a moisture factor m that rises linearly with the
     aridity index P/PET from 0 at 0.2 to 1 at 1.0 (`forestAridity`),
     and the trees relax toward f m v with the treeline's times (10
     years rising snow-free, 3 falling). Grass is the rest of the cover,
     v − trees, and keeps the cover's own times (180 days rising, 365
     falling, 720 under snow).
   - Why P/PET. The bulk potential evaporation the land's exchange
     implies, E with wetness 1 from the skin, reads the skin's heat: a
     dry surface runs hot, so the drier the land the larger its
     "demand". Four ten-day means of the model's first year (nine64 and
     eight64 copied and run from days 91, 183, 274 and 365, the means
     kept with a 100-day memory and the start's weight removed) by band,
     mm/yr:

     | band | 60–70N | 50–60N | 40–50N | 30–40N | 20–30N | 0–10N | 10–20S | 20–30S |
     |---|---|---|---|---|---|---|---|---|
     | rain | 718 | 744 | 735 | 702 | 665 | 1487 | 838 | 649 |
     | bulk potential (skin) | 748 | 1356 | 2355 | 3151 | 4803 | 3364 | 4928 | 4346 |
     | FAO-56 reference | 288 | 490 | 709 | 854 | 1122 | 1132 | 1211 | 1131 |
     | P/PET (reference) | 2.50 | 1.52 | 1.04 | 0.82 | 0.59 | 1.31 | 0.69 | 0.57 |

     The reference by region: Europe 627, Siberia 452, Congo 1163,
     Amazon 1339, Sahara 1355 (the hot dry skin cutting the net
     radiation); the bulk potential reads 1855 in Europe and 6739 in the
     Sahara. No regional values of Earth's reference evaporation were
     read (Zomer et al. 2022, Sci. Data, doi:10.1038/s41597-022-01493-1,
     the FAO-56 global database, gives them only as maps). What can be checked is the
     aridity classes the index draws: by the first year's own P/PET
     (rain and reference demand of the four ten-day means, per cell off
     the ice sheets) the land is 4.3 % hyper-arid, 12.9 % arid, 17.6 %
     semi-arid, 8.0 % dry subhumid and 57.3 % humid, against UNEP's
     6.6, 10.6, 15.3, 9.0 and 58.5 % (Bastin et al. 2017's 978, 1566,
     2263 and 1326 Mha of drylands making 41.5 % of the land).
     Rain alone against P/PET and the root zone's fill, the model's
     first year (rain, reference demand, P/PET, mean fill of the four
     states) with the Earth biome of each region of the spin-up's log:

     | region | P mm/yr | PET mm/yr | P/PET | fill | Earth |
     |---|---|---|---|---|---|
     | Sahara | 331 | 1355 | 0.24 | 0.20 | desert |
     | Gobi | 348 | 703 | 0.50 | 0.28 | desert steppe |
     | Kalahari | 171 | 1356 | 0.13 | 0.27 | semi-desert, savanna |
     | Sahel | 1327 | 1298 | 1.02 | 0.47 | savanna, steppe |
     | Cerrado | 1170 | 1232 | 0.95 | 0.44 | savanna |
     | Congo | 1437 | 1163 | 1.24 | 0.51 | rainforest |
     | Amazon | 1057 | 1339 | 0.79 | 0.38 | rainforest |
     | Europe | 1186 | 627 | 1.89 | 0.63 | temperate forest |
     | Siberia | 524 | 452 | 1.16 | 0.74 | taiga |

     Rain alone puts Siberia's 524 mm below Sankaran's 650 mm ceiling
     and Staver's 1000 mm, so the taiga would be capped as an African
     savanna; its P/PET of 1.16 is humid. The fill is the model's own
     store: it compresses the range (R² 0.51 against P/PET over the
     first year's land), reads the Congo drier than Siberia, swings with
     the season (the Sahel's fills 0.04 at the year-six March equinox
     after 657 mm), and the cover v already follows it (its goal ramps
     with the fill from 0.1 to 0.6), so gating trees by it would count
     the bucket twice. P/PET is the index Earth's classes are drawn in.
   - The limits. UNEP's classes (as Bastin et al. 2017, Science 356:
     635–638, give them): hyper-arid below 0.05, arid 0.05–0.2,
     semi-arid 0.2–0.5, dry subhumid 0.5–0.65, humid above. Bastin's
     Table 1 (2015, photo-interpreted plots, Mha): of 1566 Mha arid land
     103 carry at least 10 % tree cover and 28 at least 40 %; semi-arid
     2263: 559 and 276; dry subhumid 1326: 652 and 469. Taking 0.6 for
     the closed class, 0.22 for the open one and 0–0.02 for the rest (my
     assumption, not Bastin's), the mean tree cover is 0.02–0.04 arid,
     0.10–0.12 semi-arid and 0.24–0.25 dry subhumid. The dry limit 0.2 is
     the arid border, below which 1.8 % of the land carries a closed
     canopy; the wet limit 1.0 makes the year-one land (v about 0.5)
     carry those means: 0.00 arid, 0.09 semi-arid, 0.24 dry subhumid
     (0.75 instead gives 0.13 and 0.35). In rainfall: maximum woody cover
     in African savannas rises linearly with mean annual rain up to about
     650 mm and is held below closure above it by fire and herbivory
     (Sankaran et al. 2005, Nature 438: 846–849, 854 sites); between
     1000 and 2500 mm with a dry season under 7 months tree cover is
     bimodal and only fire tells savanna from forest (Staver et al. 2011,
     Science 334: 230–232); Hirota et al. (2011, Science 334: 232–235)
     find forest, savanna and treeless states, bimodal at about
     1000–2000 mm (as Baudena et al. 2015, Biogeosciences 12: 1833–1848,
     summarise it; Hirota's own numbers not read). With a tropical
     reference demand of about 1300 mm, 650 mm is P/PET 0.5 and 2500 mm
     1.9; the ramp, with no fire in the model, stands for the mixture.
   - Grass's albedo (`grassland`). The vegetated albedo runs from
     `grassAlbedo` 0.20 to `forestAlbedo` 0.13 with the trees' share of
     the cover (trees above v count as all forest), the ground showing
     through at 1 − v as before. Crops and natural vegetation under 1 m
     with a full green cover 0.18–0.25 (Oke 1987, Boundary Layer
     Climates, p. 132); forest 0.12–0.15 is the class table's textbook
     range (Oke's Table 1.1 is an image in the copy found, not read). Snow buries grass: a snow
     cover's own albedo less `grassSnowDarkening` 0.06 times the grass's
     share of the cell, then masked by the trees as before. MODIS
     broadband (0.3–5.0 µm) white-sky albedo in the presence of snow,
     Northern Hemisphere 2000–04 (Moody et al. 2007, Remote Sens.
     Environ. 111: 337–345, Table 1; standard deviations 0.10–0.15):
     grassland 0.59, cropland 0.58, open shrubland 0.56, savanna 0.47,
     barren 0.65, permanent snow 0.74, evergreen needleleaf forest 0.27,
     deciduous needleleaf 0.33, deciduous broadleaf 0.31, mixed 0.29.
     These are means over snow-flagged pixels, partial and patchy cover
     included, which the model's own snow fraction carries; 0.06 is
     grassland's (and cropland's 0.07) difference from barren land under
     the same sampling, the grass's own effect. Snow-covered grassland in
     the model reads its snow's albedo less 0.06 g: 0.79 fresh, 0.44 at
     the 0.50 floor, at full grass.
   - The cover v keeps its dynamics. In the year-six state the cover
     follows the last year's rain of the log (the four 91-day means to
     day 2190): Sahara 0.20 at 265 mm, Gobi 0.29 at 301, ausNorth 0.38 at
     456, Sahel 0.48 at 657, Europe 0.84 at 812, Congo 0.96 at 1387,
     Amazon 0.90 at 1360; Borneo 0.47 at 739 and New Guinea 0.32 at 557
     where the model's rain is a third of Earth's. Nothing there calls
     for a change: the cover is low where the model is dry.
   - Starts. A state without the moisture means (every state before
     this, and a fresh start) takes the demand −2.14 + 0.0134 Q̄ mm/d
     (Q̄ the annual mean insolation at the top, W/m²; R² 0.72, rms 0.57
     mm/d over the first year's land) and the rain that demand times
     0.01 + 0.79 × fill + 0.63 × v (`MOISTURE_ESTIMATE`, the same on both
     engines through the land's load), and its trees start at f m v (a
     state with season means and trees keeps those below f m v). The
     index's fit is to the 21 regions of five64's log (last year's rain
     over the estimated demand against the year-six state's fill and
     cover, R² 0.64); the fill alone fits the first year's cells at
     0.27 + 1.69 × fill (R² 0.35). Mean |m − m_true|: over the year-six
     regions 0.111 with the fill and cover, 0.255 with the fill alone;
     over the first year's cells (v still near the atlas's 0.5) 0.251
     and 0.271. Its rain against the log's at the year-six start (mm/yr):
     Sahara 335 / 265, Sahel 416 / 657, Congo 1366 / 1387, Amazon 1573 /
     1360, Europe 790 / 812, Siberia 443 / 757, Gobi 344 / 301, Cerrado
     1084 / 694. A fresh start (fill and cover 0.5) reads P/PET 0.72 and
     m 0.65. Land added by regridding takes the estimate at the guessed
     fill and cover (0.5 snow-free); saved states, the GPU (at the end of
     PH), `land.serialize` (the page's snapshot) and land regridding
     carry both means.
   What else differs between forest and grass, ranked by its effect on
   the model's surface fluxes (the first built Oct 2: M22, the surface
   layer by roughness):

   1. Roughness. The land's exchange coefficient is one value, 1.5·10⁻³
      for momentum, heat and vapour, which at the bl34 grid's lowest
      level (20 m) is the neutral drag of a roughness length of 0.65 mm.
      Neutral at 20 m: grass (z₀ 0.03 m) C_D 3.8·10⁻³ and, with the heat
      roughness a tenth, C_H 2.8·10⁻³; forest (z₀ 1 m, 0.1 m) C_D
      1.8·10⁻² and C_H 1.0·10⁻². The IFS (Cy43r1 Part IV, Table 8.3)
      calibrates z₀ for momentum against SYNOP winds at 0.2 m for short
      grass, 0.25 m crops, 0.034 m tundra, 0.013 m desert and 2.0 m for
      all forests (heat 0.002 m over grass, 2.0 m over forest). One step
      on eight64_day0183 at its instant over snow-free land: with 1.5·10⁻³
      H 22.1, LE 42.6 W/m²; grass-like 2.8·10⁻³ H 41.3, LE 62.9; forest-
      like 1.0·10⁻² H 148.9, LE 124.6. Three N=64 days from the same state
      with the whole land at 3.8·10⁻³: on day 186 the snow-free land's
      instantaneous H 24.8 → 24.3 W/m² (its skin−air difference 2.57 →
      1.20 K), LE 38.4 → 46.7 W/m², the lowest wind 6.8 → 5.5 m/s, the
      land's skin 23.22 → 21.86 °C, the global surface 16.65 → 16.34 °C.
      Building it means a drag that follows the trees on both engines in
      the momentum, the boundary layer and the surface fluxes.
   2. Interception. The model has no canopy store: rain falls into the
      15 kg/m² surface layer, of which only the bare share evaporates.
      The IFS holds 0.2 mm per unit leaf area (0.2 mm on bare ground,
      Cy43r1 eq. 8.2), so a forest of leaf area 5 holds about 1 mm and
      grass of 2 about 0.4; the intercepted water evaporates at the
      potential rate. Forest interception loss is 10–50 % of the gross
      annual rain (Carlyle-Moses and Gash 2011, Forest Hydrology and
      Biogeochemistry, Springer, as later papers quote it; the chapter
      not read): on 2 mm/d over forest 0.2–1.0 mm/d, 6–29 W/m² of LE
      where it rains, nothing on dry days.
   3. Rooting depth. One bucket of 300 kg/m² under every cover. Maximum
      rooting depth 7.0 ± 1.2 m for trees, 5.1 ± 0.8 for shrubs, 2.6 ±
      0.1 for herbaceous plants, 2.0 ± 0.3 m in boreal forest (Canadell
      et al. 1996, Oecologia 108: 583–595). The IFS root profiles
      (Table 8.1, after Zeng et al. 1998) put 90 % of the roots above
      0.62 m for short grass, 0.67 crops, 0.75 evergreen needleleaf, 0.84
      deciduous broadleaf, 1.24 evergreen broadleaf. A one-step or
      three-day effect is nil; through a dry season a forest bucket half
      again as deep would carry 150 mm more, 50 days at 3 mm/d.
   Engines: at N=6 on random cover, trees, snow, snow albedo and moisture
   means the land albedos agree to 4.7·10⁻⁸ over 151 land cells (58 with
   grass above 0.2, 18 of them under snow; the grass moves the albedo
   by up to 0.059); over 48 steps (memories 6 h, trees 3 h up and 2 h
   down) the rain means agree to 6.1·10⁻⁴ mm/d (up to 3.1), the demand
   means to 2.1·10⁻² mm/d (up to 4.9; the engines' lowest air and net
   radiation drift apart) and the trees, thinned by the moisture on 88
   cells, to 3.4·10⁻⁵. Hand values: P/PET 0.4 gives m 0.25; a half-life
   of 2 mm/d of rain moves a mean of 1 to 1.5, snowfall counts; trees
   0.1 → 0.5 − 0.4/e in ten years at P/PET 0.6, 0.8 → 0.8/e in three
   at P/PET 0.1, under snow 0.8 → 0.125 + 0.675/e at P/PET 0.3; a cover
   of 0.6 a third under trees on dry soil 0.30 + (0.1767 − 0.30) × 0.6;
   snow of 0.85 on it 0.826 masked 0.2/0.7 of the way to 0.27; the
   estimate at the equator 3.43 mm/d of demand and 2.47 of rain. An
   older state's trees restart at f m v, or keep their own where
   those are lower; saved states keep the means
   through a state file; regridding carries them. The digests are
   unchanged (their model has no land).

   Runs: N=64 GPU, OCEAN `{"everySteps":8}`, from copies of the states
   (bl34; five64_day2190 on its own cam26 grid), before (`treeMoisture`
   and `grassland` false, the previous rule; the moisture means kept as a
   diagnostic with a 100-day memory) → after (the defaults, means and
   trees from the estimate). Surface classes, insolation-weighted
   surface albedo (`land.albedo` by the day's mean insolation at the
   cell's latitude), area share of the globe:

   | class | June solstice + 3 d | September equinox + 3 d | December solstice + 3 d | March equinox + 10 d | year six + 10 d | reference |
   |---|---|---|---|---|---|---|
   | dense forest (v > 0.7, trees ≥ half of it) | trace | 0.052: 0.165 → 0.167 | 0.053: 0.159 → 0.162 | 0.056: 0.157 → 0.160 | 0.077: 0.148 → 0.150 | 0.12–0.15 |
   | dense vegetation (v > 0.7) | trace | 0.056: 0.165 → 0.170 | 0.054: 0.159 → 0.163 | 0.058: 0.157 → 0.162 | 0.084: 0.148 → 0.154 | |
   | grassland (v ≥ 0.4, trees < 0.3 of it) | 0.014 → 0.050: 0.199 → 0.241 | 0.009 → 0.047: 0.180 → 0.243 | 0.002 → 0.010: 0.180 → 0.239 | 0.002 → 0.010: 0.179 → 0.239 | 0.002 → 0.009: 0.153 → 0.232 | 0.18–0.25 |
   | partly vegetated (v 0.2–0.7) | 0.259: 0.200 → 0.215 | 0.191: 0.213 → 0.232 | 0.122: 0.217 → 0.234 | 0.134: 0.216 → 0.231 | 0.088: 0.219 → 0.236 | 0.18–0.25 |
   | bare dry soil (v < 0.2) | | | | | 0.029: 0.278 → 0.286 | 0.30–0.40 |
   | snow under forest (trees ≥ 0.5, ≥ 20 kg/m²) | | trace | 0.030 → 0.023: 0.333 → 0.319 | 0.024 → 0.020: 0.377 → 0.376 | 0.024 → 0.022: 0.343 → 0.342 | 0.27–0.33 (MODIS) |
   | snow on grassland (trees < 0.2, v ≥ 0.3) | | trace | 0.003 → 0.007: 0.735 → 0.732 | 0.003 → 0.010: 0.725 → 0.692 | 0.004 → 0.005: 0.725 → 0.708 | 0.59 ± 0.14 (MODIS, patchy included) |
   | tundra snow (f < 0.1) | trace | 0.002: 0.802 → 0.773 | 0.007: 0.845 → 0.808 | 0.010: 0.830 → 0.799 | 0.010: 0.830 → 0.802 | 0.60–0.80 |
   | snow 40–60N under a cover < 0.5 | | | 0.015: 0.385 → 0.498 | 0.013: 0.457 → 0.595 | 0.003: 0.527 → 0.629 | grassland 0.59, cropland 0.58 |
   | snow 50–70N (boreal), all | | | 0.049: 0.348 → 0.398 | 0.045: 0.459 → 0.501 | 0.040: 0.442 → 0.448 | |

   The trees and grass by band after ten days (cover v, trees, grass;
   the bands' P/PET from their own means), year one (from the March
   equinox) and year six:

   | band | year one | year six |
   |---|---|---|
   | 70–80N | 0.45, 0.15, 0.30 | 0.31, 0.11, 0.20 |
   | 60–70N | 0.50, 0.36, 0.14 | 0.47, 0.39, 0.09 |
   | 50–60N | 0.50, 0.45, 0.06 | 0.64, 0.62, 0.02 |
   | 40–50N | 0.47, 0.35, 0.12 | 0.55, 0.46, 0.10 |
   | 30–40N | 0.53, 0.34, 0.19 | 0.58, 0.39, 0.19 |
   | 20–30N | 0.52, 0.32, 0.20 | 0.45, 0.26, 0.18 |
   | 10–20N | 0.57, 0.30, 0.27 | 0.39, 0.14, 0.25 |
   | 0–10N | 0.68, 0.55, 0.13 | 0.69, 0.48, 0.21 |
   | 0–10S | 0.66, 0.55, 0.11 | 0.74, 0.66, 0.08 |
   | 10–20S | 0.54, 0.38, 0.16 | 0.53, 0.37, 0.16 |
   | 20–30S | 0.51, 0.28, 0.22 | 0.58, 0.40, 0.18 |
   | 30–40S | 0.67, 0.46, 0.20 | 0.75, 0.56, 0.18 |

   Before, the trees were f v (year six: 0.86 under the dense
   vegetation's cover of 0.88, 0.43 under the partly vegetated land's
   0.43, 0.12 under the bare land's 0.12). By the aridity class of the
   state's means, year six after ten days (still the estimate's: ten
   days move them by 1 %, so the classes and the trees both come from
   it; the check against the model's own P/PET follows) (area share of the globe; cover, trees, grass; share of the land with
   trees ≥ 0.1 and ≥ 0.4) against Bastin's land with tree cover ≥ 10 %
   and ≥ 40 %: arid (0.019) 0.12, 0.00, 0.12; 0.00, 0.00 against 0.066,
   0.018; semi-arid (0.054) 0.32, 0.07, 0.25; 0.25, 0.00 against 0.247,
   0.122; dry subhumid (0.026) 0.54, 0.25, 0.29; 0.98, 0.09 against
   0.492, 0.354; humid 0.65–1 (0.047) 0.61, 0.45, 0.16; humid 1–2
   (0.113) 0.74, 0.68, 0.06. The zone means match; the spread does not:
   the model's trees are a smooth fraction of every cell where Earth's
   are bimodal (closed or open), so the dry subhumid land is almost all
   lightly wooded and rarely closed. By rain (mm/yr): 250–500 0.38,
   0.21, 0.17; 500–650 0.57, 0.39, 0.18; 650–1000 0.72, 0.57, 0.15;
   1000–1500 0.85, 0.75, 0.10; 1500–2500 0.94, 0.89, 0.05. By region (cover, trees,
   grass, P/PET of the means): Sahara 0.21, 0.05, 0.16, 0.30; Sahel
   0.46, 0.13, 0.33, 0.34; Kalahari 0.61, 0.33, 0.28, 0.58; Gobi 0.29,
   0.10, 0.19, 0.44; Congo 0.96, 0.89, 0.07, 1.09; Amazon 0.90, 0.86,
   0.04, 1.28; Europe 0.85, 0.83, 0.02, 1.27; Siberia 0.59, 0.59, 0.00,
   1.17; Cerrado 0.58, 0.48, 0.10, 0.91.

   Against the model's own P/PET (the first year's four ten-day means
   per cell, v from nine64's day 365, f from the season estimate), the
   trees at f m v by class: hyper-arid and arid 0.000, semi-arid 0.090
   (trees ≥ 0.1 on 0.38 of it, ≥ 0.4 on none), dry subhumid 0.243 (0.97,
   0.07), humid 0.65–1 0.420, 1–2 0.573, above 2 0.515: the limits put
   the class means on Bastin's with the model's own rain and demand, not
   only with the estimate's. The estimate puts them on the wrong cells:
   over the land whose season factor is at least 0.5, the estimate's m
   by the own P/PET class is 0.34 where the own m is 0.00 (own arid,
   0.18 of that land), 0.52 against 0.18 (semi-arid), 0.62 against 0.47
   (dry subhumid), 0.71 against 0.77, 0.83 against 1.00 and 0.93 against
   1.00 (humid 0.65–1, 1–2, above 2); the two agree on which side of m
   = 0.5 a cell lies on 0.74 of that land, mean |Δm| 0.245.

   Review (Oct 1). The CPU land stepped through three wet years, three
   dry, two of regrowth and one under snow (6 h steps) matches a float64
   reference of the moisture means, cover, season means and trees to
   5·10⁻¹⁴; the GPU over 48 steps at N=6 matches an f32 hand formula of
   the rain mean to 3.2·10⁻⁷ (relative) and of the trees to 6.0·10⁻⁸
   over 4315 cell-steps (672 under snow), and the potential its first
   step implies matches the CPU's to 2.2·10⁻⁵ mm/d. Land surfaces absorb
   (1 − α) of the shortwave reaching them on both engines (one albedo
   for the direct and diffuse beams). Under snow the cover decays over
   720 days while the trees hold, so trees can stand above the cover;
   the grass is then 0 and the cover reads as forest. A state with
   season means and trees but no moisture means keeps trees below f m v
   instead of raising them to it. Three N=64 days from nine64_day0365
   and five64_day2190 (tags bmr365, bmr2190), the class table with
   `grassland` false → true on the same states: dense vegetation 0.156
   → 0.162 and 0.147 → 0.155 (0.12–0.15), partly vegetated 0.218 →
   0.232 and 0.221 → 0.236 (0.18–0.25), bare dry soil 0.277 → 0.285
   (year six), snow on open land, cold (standing < 0.2) 0.773 → 0.753
   and 0.806 → 0.783 against the class's 0.80–0.85: the grass darkens
   deep snow as much as patchy snow, where Moody's 0.06 is a mean over
   both. Snow at 40–60N under a cover below 0.5 reads 0.604 and 0.661
   after three days.

   Day means (outcomes), before → after: day 94 albedo 0.326 → 0.328,
   ASR 229.6 → 228.8, OLR 240.6 → 240.5, SWCRE −58.7 → −58.3, clear-sky
   reflectance 0.1535 → 0.1569, Ts 16.57 → 16.52 °C; day 101 0.336 →
   0.337, 225.9 → 225.8, 240.2 → 240.1, 0.1525 → 0.1561, 16.42 → 16.33;
   day 186 0.312 → 0.314, 234.2 → 233.4, 242.2 → 242.1, 0.1493 → 0.1523,
   16.65 → 16.60; day 193 0.326 → 0.325, 229.6 → 229.7, 0.1509 → 0.1538,
   16.07 → 15.99; day 277 0.349 → 0.350, 221.7 → 221.2, 0.1650 → 0.1671,
   13.56 → 13.53; day 284 0.359 → 0.360, 218.1 → 218.0, 0.1636 → 0.1657,
   13.11 → 13.07; day 368 0.328 → 0.330, 228.7 → 228.1, 0.1602 →
   0.1636, 13.69 → 13.64; day 375 0.345 → 0.347, 223.0 → 222.3, 233.2 →
   233.3, −62.5 → −62.1, 0.1613 → 0.1649, 13.73 → 13.67; day 2193 0.302
   → 0.304, 237.6 → 237.1, 0.1588 → 0.1607, 15.46 → 15.43; day 2200
   0.310 → 0.312, 235.0 → 234.3, 240.4 → 240.3, −51.0 → −50.9, 0.1602 →
   0.1624, 15.73 → 15.69. Rain 2.19–2.73 mm/d, within 0.07 of before.
   The runs shared the GPU with other jobs; their wall times say nothing.

   What it cannot do:

   - No fire and no competition for light: the moisture factor is a
     ramp fitted to the zone means, so where Earth splits into closed
     forest and open savanna (Staver's 1000–2500 mm) every cell carries
     the mixture's fraction; no herbivory, no succession, one tree type
     and one grass (no evergreen against deciduous, no C3 against C4).
   - The model's own climate decides: its Sahara rains 265–340 mm/yr
     (P/PET 0.24–0.30, Earth's below 0.05), so it is a grassy semi-desert,
     and its mid-latitude continents are humid (40–50N P/PET 1.04 in the
     first year), so the 40–60N snow under a sparse cover is still masked
     by trees of 0.2 (0.595 against MODIS grassland's 0.59, but under a
     cover the masking assigns to trees).
   - Dense forest still reads 0.148–0.167 against 0.12–0.15: the ground
     showing through at 1 − v (0.25 at v 0.75) is the dry soil.
   - Snow on grassland reads 0.69–0.73 on full snow; MODIS's 0.59 mixes
     in patchy snow, which the model's own snow fraction carries in the
     thin-snow class.
   - Ten days move the means by 1 % of their gap and the trees by under
     1 %: what the runs show is the start, the estimate's trees.
   - The estimate's error (mean |m − m_true| 0.11 on the year-six regions,
     0.25 on the first year's cells) stays in the trees for about three
     years, the means' memory.
   - The reference evapotranspiration uses the model's own albedo and
     skin in the net radiation, not FAO-56's 0.23 reference surface; the
     hot dry skin lowers it in deserts (Sahara 1355 mm/yr against
     the bulk potential's 6739).
   - Interception and rooting depth do not yet differ between forest and
     grass (above); roughness does since Oct 2 (M22, the surface layer by
     roughness).

   Defaults: land `treeMoisture` true, `moistureMemory` 3 years,
   `forestAridity` [0.2, 1.0] (P/PET), `grassland` true, `forestAlbedo`
   0.13, `grassAlbedo` 0.20, `grassSnowDarkening` 0.06,
   `MOISTURE_ESTIMATE` demand [−2.14, 0.0134] (mm/d against W/m²),
   aridity [0.01, 0.79, 0.63] (constant, fill, cover); radiation
   `REFERENCE_RESISTANCE` 70 s/m; `vegetatedAlbedo` 0.13 only without
   grassland. The spin-up takes land options through `LAND`. Older
   states start their moisture means from the estimate and their trees
   at f m v, or at their own where those are lower.
   Soil colour from the soil's own carbon (Oct 2). The dry bare soil's
   albedo follows a topsoil organic-carbon store that each cell's own
   cover and climate build, on both engines (`js/physics/land.module.js`,
   `js/gpu/physics.gpu.js`), with parity and hand-computed unit tests
   (`test/soilCarbon.test.mjs`). No map of soils or soil colour enters.

   - The store (`soilCarbon`, kg C/m² in the top 0.2 m, Jobbágy & Jackson
     2000's top interval; one pool): dS/dt = A (I − k S), stepped exactly
     for each step's I and k. Litter I = `litterInput` × (`treeLitter` ×
     trees + `grassLitter` × grass) × min(1, fill / 0.75) × m(T), with
     m(T) = 1 / (1 + e^(1.315 − 0.119 T)) in the lowest air (°C) while it
     is at least 0.9 °C and 0 otherwise: the temperature curve of Lieth's
     (1975) Miami model of net primary production (there in the annual
     mean temperature; here step by step inside TREELIM's growing
     season), and min(1, fill / 0.75) the root-zone factor the
     transpiration already uses. Decomposition k = r(T) M(fill) /
     `soilTurnover`, r = exp(308.56 (1/56.02 − 1/(T − 227.13))) of the
     lowest air in K (Lloyd & Taylor 1994, Funct. Ecol. 8: 315–323): 1 at
     10 °C, 0.30 at 0 °C, 0.047 at −10 °C, 0.002 at −20 °C, 2.30 at
     20 °C, 4.26 at 30 °C. M is TRIFFID's moisture factor (Cox 2001,
     Hadley Centre Technical Note 24, as given by Clark et al. 2011,
     Geosci. Model Dev. 4: 701–722, eq. 67, read: 0.2 up to s_w, linear
     to 1 at s_o = 0.5 (1 + s_w), then 1 − 0.8 (s − s_o); Cox's original
     lower bound s_w, which JULES later raised to 1.7 s_w) in the
     bucket's fill: 0.2 below the wilting fill 0.1
     (`decompositionWilting`), 1 at 0.55, 0.64 full; 0.2 whenever the air
     is below freezing, as JULES's s counts only the unfrozen water. There
     s is the top layer's water as a fraction of saturation; here it is
     the fill, the share of the root zone's capacity, so a full bucket
     decomposes at 0.64 where a loam at field capacity would sit near
     s_o and decompose at about 1 (recalled loam values). Mapping the fill linearly onto a saturation
     fraction from wilting to field capacity (loam, clay, sandy loam; the
     hydraulic values recalled, not read) fits the seven climates below
     worse, rms log ratio to J&J 0.67–0.73 against 0.58, the desert at
     0.4 kg/m² against 2.05, and leaves the boreal forest at 10.7–11.4.
   - Trees against grass. Jobbágy & Jackson (2000, Ecol. Appl. 10:
     423–436, Table 4): temperate grasslands keep 0.21 of their biomass
     above ground and 70 % of their roots in the top 20 cm, forests
     0.75–0.85 above ground; root-derived carbon stays in soil 2.4 times
     as long as shoot-derived (Rasse et al. 2005, Plant Soil 269:
     341–356). Per unit of production grass therefore feeds the mineral
     topsoil more; per unit of cover forest produces more (from memory,
     about 1.5–2 times grassland's NPP in one climate, not checked). The
     stocks they leave in 0–20 cm (J&J Tables 3 and 4, the first metre
     times its top-20-cm share): temperate grassland 4.8 kg/m², temperate
     deciduous and evergreen forest 9.0 and 6.8, under wetter climates.
     The weights are 1 and 1. What makes grassland soils dark in the
     field is depth: SOC under grass lies deeper (41 % of the first metre
     in the top 20 cm against temperate forests' 47–52 %), the thick
     mollic A horizon,
     which one 0.2 m pool does not resolve.
   - Calibration. The equilibrium S* = ⟨I⟩/⟨k⟩ depends on `litterInput` ×
     `soilTurnover` alone, 35 kg/m², the geometric-mean fit (34.8, rms
     log ratio 0.58) to J&J's 0–20 cm stocks, Table 3's first metre times
     Table 4's top-20-cm share, for seven sine-year climates (mean ±
     amplitude °C, fill, cover → model S*, J&J, carbon %, dry albedo;
     `carbonEquilibrium` integrates the year by Simpson's rule between
     the rates' jumps at 0.9 °C and freezing, converged to 10⁻⁷): tropical
     evergreen forest 26 ± 1.5, 0.55, 1.0 → 6.45, 8.2, 2.5 %, 0.141;
     tropical savanna 25 ± 3, 0.3, 0.6 → 3.94, 4.8, 1.5 %, 0.175;
     temperate deciduous forest 10 ± 11, 0.6, 0.95 → 10.5, 9.0, 4.0 %,
     0.124; temperate grassland 7 ± 14, 0.3, 0.6 → 5.58, 4.8, 2.1 %,
     0.149; boreal forest −3 ± 18, 0.7, 0.7 → 11.3, 4.7, 4.4 %, 0.123;
     tundra −11 ± 16, 0.75, 0.4 → 8.57, 5.7, 3.3 %, 0.129; desert 22 ± 8,
     0.12, 0.1 → 0.66, 2.0, 0.25 %, 0.314. Without the water factor on the
     litter the desert reads 2.2 kg/m² and 0.227; with the growth warmth
     (5–15 °C) in place of Lieth's curve the tundra gets no litter.
     J&J's deserts include shrub deserts; their top-20-cm SOC
     (0.8 % at 260 kg/m²) would read 0.23 by the curve below.
   - Time. `soilTurnover` 70 years at 10 °C and optimal moisture with
     `litterInput` 0.5 kg C/m² a year at full cover, warmth and water
     (about a third of Lieth's ceiling of 3000 g dry matter, 1.4 kg C,
     /m²/yr): the year's mean decomposition gives real turnover times of
     21 years in the tropical forest climate, 63 temperate, 108 desert,
     129 steppe, 171 boreal, 628 tundra. Carvalhais et al. (2014, Nature
     514: 213–217) give whole-ecosystem carbon 23 (+7/−4) years globally,
     15 near the equator, 255 north of 75N; TRIFFID's κ_s 0.5·10⁻⁸ s⁻¹ at
     25 °C is 6.3 years there, 20 at 10 °C under Lloyd–Taylor (recalled).
     `carbonAcceleration` A = 100 multiplies input and decomposition
     alike, leaving S* where it is under steady forcing: e-folding times
     0.2 years tropical, 0.6 temperate, 1.1 desert, 1.3 steppe, 1.7
     boreal, 6.3 tundra. Under the seasons the store's swing correlates
     with the decomposition, so the year's mean moves with A: stepped
     daily to its periodic state, the seven climates' annual means at
     A = 100 lie within 0.7 % of A = 1's (temperate forest +0.6 %), at
     A = 1000 up to +9 %; a monsoon year (27 ± 3 °C, fill 0.35 ± 0.30 in
     phase, grass 0.4) −0.7 % at A = 100. What
     it loses is the soil's lag: a cell whose cover dies keeps its dark
     soil for decades to centuries on Earth (a 0.2 m store at 63 years in
     a temperate climate) and for months to a year here, so the soil
     brightens with the cover's own 365-day decline instead of after it,
     and the seasonal swing of a store follows the seasons more than a
     real one.
   - Colour. Dry albedo = `mineralAlbedo` − (`mineralAlbedo` −
     `humusAlbedo`) (1 − exp(−c / `organicScale`)), c = 100 S / `topsoilMass`
     (% organic carbon; 260 kg/m², 0.2 m at an assumed 1300 kg/m³). The
     shape: soil reflectance falls with organic matter in all visible and
     near-infrared bands, other constituents mask organic matter below
     about 2 % (1.2 % C; Baumgardner et al. 1985, Adv. Agron. 38, as
     quoted by later studies) and above about 5 % (2.9 % C) more changes
     little (Page 1974, as quoted); dry Munsell value falls
     logarithmically with organic carbon (Konen et al. 2003, SSSAJ 67:
     1823–1830, r² 0.74 air-dry, 130 Iowa Ap horizons on one parent
     material; steep below 1 % C in a multi-state set, as summarised);
     broadband albedo = 0.069 × Munsell value − 0.114 (Post et al.
     2000, SSSAJ 64: 1027–1034, 26 smoothed US soils air-dry and wet, the
     value of each as measured, r² 0.93, as the abstract gives it):
     value 7 0.369, 6
     0.300, 5 0.231, 4 0.162, 3 0.093. `organicScale` 1 % C puts 70 % of
     the darkening below 1.2 % C and 95 % below 3 %. `mineralAlbedo` 0.37
     is Post's value 7, inside desert sand and rock's 0.30–0.40;
     `humusAlbedo` 0.12 is the darkest dry CLM soil colour, (vis 0.08 +
     nir 0.16)/2 of class 20, whose brightest class 1 is (0.36 + 0.61)/2
     = 0.485 (CLM5 Technical Note, Table 2.3.3; the colours fitted to
     MODIS by Lawrence & Chase 2007). Carbon 0.25, 0.6, 1, 2.5, 5 % →
     0.315, 0.257, 0.212, 0.141, 0.122. The wet darkening takes the dry
     value d to d − (d − 0.15 d/0.30) × the surface layer's fill, a
     saturated/dry ratio of 0.5 (Idso et al. 1975's 0.47; CLM's 0.50 for
     class 20 to 0.69 visible, 0.82 near-infrared for class 1). With no
     carbon over a 0.30 mineral, or `soilCarbon` false, the land albedo
     is the fixed soil's to the bit (on five64_day2190's 11 882 land cells
     against the previous commit, and in the test).
   - Starts and plumbing. A fresh start and a state without the field
     (every state before this) start each cell at S* for its own cover,
     trees and present bucket fill under SEASON_ESTIMATE's sine year of
     the lowest air (`airCycle`), on both engines through the land's load;
     a state with the field keeps it. Ice sheets hold 0. Saved states,
     the GPU (at the end of PH, `SOILC`), `land.serialize` (the page's
     snapshot) and land regridding (S* at the guessed cover 0.5 for new
     land) carry it. Every option of `createLandSurface` reaches the GPU
     through `VEGETATION_OPTIONS` or the named heat, bucket, albedo, snow
     and ageing options (the test reads the signature); a land
     `latentHeatFusion` that differs from the sea ice's is refused, the
     engines sharing one.

   Tests. Hand values: a forest cell at 26 °C, fill 0.55, cover 1 from 0
   reaches 4.0140 kg/m² in 0.2 years and S* 6.4594 (dry albedo 0.14084);
   a desert at 30 °C, fill 0.12, grass 0.1: 0.50518 (0.32585); a tundra
   summer at 8 °C, fill 0.75, grass 0.4: 8.3839 (0.12994); a tundra winter
   at −5 °C keeps 9.8112 of 10 over half a year (no litter, frozen
   decay); a cell losing its cover at 20 °C, fill 0.55: 6.4594 → 1.2466 in
   half a year, the dry soil 0.1408 → 0.2748. A = 1 over 100 years equals
   A = 100 over one; three sine-year climates stepped 1460 times a year
   hold S* within 0.1 % at A = 1 and within 1 % after 30 years at A = 100.
   Over 2000 random cells the dry albedo stays within 0.12–0.37 and the
   land albedo never rises with the carbon under any cover, trees,
   surface layer or snow. Stepped by hand (review, Oct 2): a warm
   vegetated year, a cold year, a dry year and a lost cover, 1460 steps
   at A = 1, 100 and 3·10⁵, agree with the exact solution
   I/k + (S − I/k) e^(−A k Δt) in float64 to 3·10⁻¹⁴ kg/m²; one GPU
   physics pass over the N=6 cells against the same step in float32
   arithmetic agrees to the bit at A = 100 and to 1.1·10⁻⁶ kg/m² at
   A = 3·10⁵ (steps up to 1.35 kg/m²). A state's carbon through the
   state file reloads to 2.4·10⁻⁷ kg/m² (float32) on the CPU, to the bit
   on the GPU. At N=6 with random cover, trees, fill, surface layer
   and carbon the engines' land albedos agree to 3.6·10⁻⁸ (the carbon moves
   them by up to 0.135); over 24 steps at A = 3·10⁵ the carbon rose on 88
   cells and fell on 31, by up to 18.3 kg/m², engines apart by
   5.3·10⁻⁴ kg/m². The treeline parity test runs with `soilCarbon` false:
   darker soils move its N=6 trajectory (lowest air apart by 9.6·10⁻³ K
   against 2.9·10⁻³; its tree cover 1.3·10⁻⁴ apart against the test's
   10⁻⁴). Digests unchanged (their model has no land).

   The year-six state five64_day2190 at the start (S* from its cover and
   fill, carbon %, dry albedo; class by its own season means and P/PET,
   area share of the globe):

   | class | globe | cover, trees | fill | P/PET | carbon kg/m² (%) | dry albedo | reference |
   |---|---|---|---|---|---|---|---|
   | desert (P/PET < 0.2) | 0.020 | 0.12, 0.00 | 0.07 | 0.13 | 0.40 (0.15) | 0.336 | sand and rock 0.30–0.40; J&J desert 2.0 kg/m² |
   | steppe and savanna (0.2–0.65) | 0.079 | 0.39, 0.13 | 0.21 | 0.42 | 2.24 (0.86) | 0.238 | J&J temperate grassland 4.8, savanna 4.8 |
   | temperate forest (humid) | 0.057 | 0.76, 0.71 | 0.78 | 1.10 | 11.5 (4.4) | 0.133 | J&J 6.8–9.0 |
   | boreal (humid, season < 0.6 yr) | 0.049 | 0.56, 0.51 | 0.94 | 1.10 | 12.4 (4.8) | 0.125 | J&J 4.7 |
   | tundra (treeline factor < 0.1) | 0.013 | 0.51, 0.00 | 0.92 | 1.06 | 9.4 (3.6) | 0.185 | J&J 5.7 |
   | tropical forest (humid, < 25°, no cold season) | 0.042 | 0.84, 0.78 | 0.71 | 1.11 | 9.2 (3.5) | 0.136 | J&J 8.2 |

   The class means of the dry albedo average the exponential over cells
   (a tundra mean of 3.6 % reads 0.185 because its sparse cells are
   pale). No measured broadband albedo by soil group was read beyond the
   desert range; by Post's regression the chernic horizon's moist value
   ≤ 3 (WRB 2022, as summarised, not read) would be dark, about 0.09–0.16
   dry. By band (cover, carbon kg/m², dry albedo): 0–10N 0.70, 4.8,
   0.203; 10–20N 0.40, 1.7, 0.285; 20–30N 0.45, 4.4, 0.220; 30–40N 0.58,
   8.2, 0.185; 40–50N 0.55, 8.8, 0.167; 50–60N 0.64, 12.8, 0.128; 60–70N
   0.48, 11.2, 0.125; 70–80N 0.31, 7.1, 0.143; 80–90N 0.18, 4.0, 0.180;
   0–10S 0.73, 8.3, 0.172; 10–20S 0.52, 6.3, 0.196; 20–30S 0.58, 5.8,
   0.186; 30–40S 0.75, 7.4, 0.156; 40–50S 0.95, 12.9, 0.126.

   Runs: N=64 GPU, OCEAN `{"everySteps":8}`, three days from copies of
   nine64_day0091, eight64_day0183, nine64_day0365 (bl34) and
   five64_day2190 (cam26), before (the previous commit) → after (the
   defaults, carbon from S*). Snow-free land, insolation-weighted surface
   albedo over the end day (`land.albedo`), area share of the globe:

   | class | June solstice + 3 d | September equinox + 3 d | March equinox + 3 d | year six + 3 d | reference |
   |---|---|---|---|---|---|
   | bare dry soil (v < 0.2, layer < half) | | | | 0.029: 0.285 → 0.300 | 0.30–0.40 |
   | grassland (v ≥ 0.4, trees < 0.3 v) | 0.050: 0.241 → 0.198 | 0.047: 0.243 → 0.214 | 0.010: 0.239 → 0.206 | 0.009: 0.233 → 0.218 | 0.18–0.25 |
   | partly vegetated (v 0.2–0.7) | 0.259: 0.215 → 0.159 | 0.191: 0.232 → 0.199 | 0.128: 0.232 → 0.186 | 0.086: 0.236 → 0.194 | 0.18–0.25 |
   | dense vegetation (v > 0.7) | trace: 0.178 → 0.128 | 0.056: 0.170 → 0.133 | 0.057: 0.161 → 0.135 | 0.083: 0.155 → 0.138 | 0.12–0.15 |
   | dense forest (v > 0.7, trees ≥ v/2) | trace: 0.178 → 0.128 | 0.052: 0.167 → 0.129 | 0.054: 0.159 → 0.133 | 0.075: 0.149 → 0.133 | 0.12–0.15 |

   Their soils' dry albedo after (carbon kg/m²): bare dry soil 0.316
   (0.66); grassland 0.210–0.244 (3.4–5.6); partly vegetated 0.164–0.220
   (3.4–6.0); dense forest 0.126–0.132 (9.9–11.9). Dense forest comes
   inside 0.12–0.15 on every state: the forest floor's own soil (0.13)
   shows at 1 − v in place of 0.30. The partly vegetated land of the June
   state reads 0.159, low by 0.021: the atlas start's cover of 0.5 a year
   on, over humid mid-latitude land whose estimated store is 6.0 kg/m².
   The class table (`scripts/clearSkyBudget.mjs`, surface albedo; LAND
   `{"soilCarbon":false}` before): day 186 partly vegetated 0.233 →
   0.200, dense vegetation 0.171 → 0.133 (high by 0.021 → matches);
   year six lit over its day bare dry soil 0.286 → 0.301 (low by 0.014 →
   matches), partly vegetated 0.237 → 0.195, dense vegetation 0.155 →
   0.139 (high by 0.005 → matches). Outcomes: clear-sky albedo of the
   state lit 0.152 → 0.146 (day 186), 0.161 → 0.157 (year six); land at
   the top 0.225 → 0.201 and 0.252 → 0.239. Day means (spin-up log),
   before → after: day 94 albedo 0.328 → 0.320, ASR 228.8 → 231.7, OLR
   240.5 → 240.7, SWCRE −58.3 → −59.7, LWCRE 20.1 → 20.2, clear-sky
   reflectance 0.1569 → 0.1441, Ts 16.52 → 16.69 °C, rain 2.19 → 2.21
   mm/d; day 186 0.314 → 0.311, 233.4 → 234.7, 242.1 → 242.2, −55.2 →
   −56.1, 17.6 → 17.6, 0.1523 → 0.1462, 16.60 → 16.70, 1.70 → 1.72; day
   368 0.330 → 0.326, 228.1 → 229.4, 233.2 → 233.3, −56.7 → −57.5, 20.2 →
   20.3, 0.1636 → 0.1574, 13.64 → 13.73, 2.33 → 2.34; day 2193 0.304 →
   0.301, 237.1 → 237.9, 240.6 → 240.7, −48.7 → −49.1, 17.2 → 17.3,
   0.1607 → 0.1571, 15.43 → 15.48, 1.88 → 1.89. Three days took 0.5–1.1
   wall minutes on a shared GPU, before and after. Repeated from
   nine64_day0365 and five64_day2190 in the review (Oct 2), the day lines
   and the class table of the year-six state lit are the same to the
   digits above. The review's converged equilibrium moves single cells'
   start carbon by up to 1.9 kg/m² (26 %) and their land albedo by up to
   0.008, the land means unchanged (6.974 and 7.520 kg/m² on those two
   states), and leaves the class tables of both start states the same to
   three decimals.

   The feedback (year-six state, snow-free cells, surface sunlight taken
   as 0.75 of the annual-mean insolation at the top; d(absorbed)/d(cover),
   W/m² per unit of cover; with the soil held, and at equilibrium, where
   S* falls in proportion to the litter):

   | land | cover | sunlight W/m² | fixed 0.30 soil | carbon, soil held | carbon at S* | all cover lost: Δalbedo, Δabsorbed W/m² (fixed → carbon) |
   |---|---|---|---|---|---|---|
   | semi-arid and dry subhumid (P/PET 0.2–0.65) | 0.39 | 284 | 33.8 | 16.7 | 65.3 | 0.049 → 0.151, 13.9 → 42.8 |
   | Sahel box 10–20N, 20W–40E | 0.32 | 302 | 32.5 | 35.4 | 77.1 | 0.037 → 0.094, 11.2 → 28.4 |
   | all snow-free land | 0.58 | 277 | 38.7 | 9.5 | 37.4 | 0.088 → 0.183, 24.4 → 50.7 |

   The sign is Charney's: losing cover brightens the ground and cuts the
   sunlight it absorbs, now twice as much per unit of cover in the
   drylands once the soil has followed. Charney (1975, Q. J. R. Meteorol.
   Soc. 101: 193–202) raised the Sahel's albedo from 0.14 to 0.35 (from
   memory, not re-read), 0.21; a dryland that loses all its cover here
   gains 0.09–0.15 (0.04–0.05 with the fixed soil). With A = 100 the
   soil's part arrives within about a year of the cover's loss.

   How bare soil's albedo divides between organic matter and moisture
   on one side and parent material and iron on the other. The evidence
   read: (1) organic matter dominates the spectrum above about 2 % organic
   matter (1.2 % C) and is masked by mineral constituents below it
   (Baumgardner et al. 1985 as quoted); (2) over Northern Africa and
   Arabia, where it is below that, MODIS diffuse shortwave albedo varies
   by a factor of about 2.5 from the darkest volcanic terrain to the
   brightest sand sheets, with the soil and rock type (Tsvetsinskaya et
   al. 2002, GRL 29(9), abstract): parent material alone spans about
   0.15–0.40 there, as wide as the organic range 0.12–0.37; (3) the top
   20 cm's carbon correlates with precipitation (r² 0.33) more than
   with sand (−0.28), temperature (−0.17) or clay (0.07), one variable
   at a time (J&J Table 2): precipitation leads, texture comes second,
   ahead of temperature; (4) on one
   parent material carbon explains 74 % of the dry Munsell value (Konen
   et al. 2003); (5) the first principal component of 12 509 soil
   spectra (55 % of the variance) is iron oxides and kaolinite (Viscarra
   Rossel et al. 2016, Earth-Sci. Rev. 155: 198–230), but the spectra
   were continuum-removed, which takes out the overall brightness, so it
   says nothing of the albedo's split; (6) water halves any soil's
   albedo (Idso; CLM 0.5–0.8). No global map of bare-soil albedo (CLM's
   colour classes, MODIS) was read, so no global share was computed.
   Plainly: on humid vegetated land, roughly the 58 % outside UNEP's
   drylands (Bastin et al. 2017), organic matter and water set the soil's
   darkness and parent material matters little; on drylands, 42 %,
   parent material and iron set it, with a spread as wide as the whole
   organic range; the single largest contrast, pale deserts against dark
   vegetated soils, follows climate. Neither "the climate-driven part is
   the larger" nor "parent material dominates" holds over all land;
   each holds where its regime applies.

   What it cannot do:

   - One mineral colour: no pale quartz sand against dark basalt or
     red iron-rich soils, no salt pans, no carbonate crusts; the model's
     deserts read 0.30–0.34 dry where Earth's span 0.15–0.45.
   - One 0.2 m pool: no litter layer, no peat, no depth of the dark
     horizon; boreal soils read the darkest (12.4 kg/m² against J&J's
     4.7, whose profiles leave out the organic layer), tundra 9.4 against
     5.7.
   - No fire, erosion, photodegradation of dry litter or grazing; the
     water factor alone keeps deserts pale.
   - The start reads the present fill and cover, so a seasonally dry cell
     started in its dry season (the Sahel's fill 0.04 in the year-six
     March) starts low; A = 100 corrects it within a year or two, the
     tundra within about six.
   - The lag of soil behind vegetation is lost (above).

   Defaults: land `soilCarbon` true, `mineralAlbedo` 0.37, `humusAlbedo`
   0.12, `organicScale` 1 % C, `topsoilMass` 260 kg/m², `soilTurnover` 70
   years, `litterInput` 0.5 kg C/m²/yr, `treeLitter` 1, `grassLitter` 1,
   `decompositionWilting` 0.1, `carbonAcceleration` 100; `LLOYD_TAYLOR`
   308.56 K, 56.02 K, 227.13 K; `MIAMI` [1.315, 0.119]; the wet
   darkening's ratio `wetSoilAlbedo` / `bareAlbedo` 0.5; `bareAlbedo`
   0.30 only without the carbon. Older states start at S*.
   An unbiased fresh start (Oct 2). A start from the atlas (no FROM
   state) had put the cover at 0.5 everywhere and filled the slow means
   from the latitude and bucket estimates, whose errors stay for the
   means' three-year memory (with the half-full bucket and half cover
   the estimate's moisture factor is 0.65 on every one of the 10640
   land cells off the ice sheets at N=64, whatever its climate).
   The two-stage start replaces that on both engines
   (`js/physics/land.module.js`, `js/gpu/physics.gpu.js`, `core.gpu.js`,
   `model.gpu.js`), with hand-computed and parity tests
   (`test/landStart.test.mjs`).

   - The record. The six slow means (`seasonLength`, `seasonWarmth`,
     `rainMean`, `demandMean`, and two new ones the soil's equilibrium
     needs: `litterMean`, the mean of min(1, fill/0.75) × m(T), and
     `decayMean`, the mean of r(T) M(fill), both over `seasonMemory`)
     share one record [age, hold, start] (`record`, float64 in state
     files). While age + Δt ≤ the memory a mean takes the weight
     Δt/(age + Δt), the plain average of every step so far; after it the
     exponential weight 1 − e^(−Δt/memory) as before. Both memories are 3
     years, so the record is a plain average to day 1095 and a running
     mean from then on. A state saved without a record (every state
     before this) has age −1 and the exponential weight at once (the CPU
     weight is the old one to the bit); the estimates stay for those
     states and for land added by regridding, which carries the record.
     The GPU takes the two weights and the hold as step parameters P[5],
     P[6], P[7], computed on the host in float64.
   - Float32. The GPU's mean is m + (x − m) w in float32. Over three
     years of steps (Math.fround in every operation, against a
     compensated float64 sum over the step count) the worst error is 1.1·10⁻⁵
     (N=64 steps, 337.5 s) and 4.3·10⁻⁵ (N=128, 168.75 s) of the mean
     for a boreal season length, 5.8·10⁻⁶ and 2.3·10⁻⁵ for its warmth,
     2.7·10⁻⁶ and 5.8·10⁻⁶ for a tundra decomposition, 6.3·10⁻⁵ and
     1.6·10⁻⁴ for rain in showers (3 % of the steps); the exponential
     mean the model already ran in float32 is off by 9.2·10⁻⁴ and
     2.1·10⁻³ on the same showers in years 3–6. At 1.6·10⁻⁴ of P/PET
     the moisture factor moves by 2·10⁻⁴, and 2.3·10⁻⁵ of the season's
     warmth moves the treeline factor by under 10⁻⁴: no compensated or
     count-based form is needed.
   - The placeholder year. During a fresh start's first year the trees
     and the topsoil carbon are held and the cover runs free. The trees
     (10 years rising, 3 falling) close 9.5 % of a rising gap in a year
     and the carbon (e-folding 0.2–6.3 years at A = 100) follows the
     trees and cover, so both carry their start for years; they are held
     at values that favour no region (`startPlaceholders`): trees at
     half the cover, the grass–forest split's middle (vegetated albedo
     0.165, half way from grass 0.20 to forest 0.13; snow masked by
     0.25/0.70 of the way at v = 0.5), and carbon 1.802 kg/m² (c = ln 2 ×
     `organicScale`, 0.69 % C), where the dry soil reads 0.245, midway
     between mineral 0.37 and humus 0.12. The cover runs free from 0.5:
     in a year it closes 63 % of a falling gap (365 days), 87 % of a
     rising one at full warmth (180 days) and 40 % of its decay under
     snow (720 days), so it follows the year; holding it would need its
     equilibrium over the seasons, which depends on which side of its
     goal it stands at each moment (the rates differ rising and falling
     and the growth reads the warmth), a record of the goal's spread
     through the year rather than a mean. What it leaves: a cell whose
     goal is 0 all year keeps 0.5 e⁻² = 0.068 of cover after two years.
     The bucket starts half full and fills or drains within months.
   - The jump (`land.jump()`). At the end of year one and again at the
     end of year two every land cell's trees are set to f m v (the
     record's treeline factor times its moisture factor times the cover
     at that moment, or under snow `snowFreeCover`, the cover at the
     cell's last snow-free step: the trees do not follow the cover under
     snow and no litter falls there, and the March equinox finds the
     northern boreal land under snow with its cover decayed, five64's
     Siberia box 0.60 against 0.75 at the September equinox) and its
     topsoil carbon to `litterInput` ×
     `soilTurnover` × litter × `litterMean` / `decayMean` (35 kg/m² ×
     the type-weighted cover × the ratio), the store's equilibrium for
     its cover, trees and record; ice sheets 0. The first jump ends the
     hold. The cover and the record are not touched: the record goes on
     as a plain average, so the second jump reads years one and two
     with equal weights. The day-365 record covers exactly one year,
     equinox to equinox. On the GPU the jump reads the land back, applies
     the same function with each value rounded to float32 as stored and
     writes the trees and carbon back, so a repeat at once changes
     nothing (to the bit, tested). The equilibrium is the ratio of the
     year's means: the store's seasonal swing at A = 100 is not set,
     and the annual mean of the periodic store lies within 0.7 % of the
     ratio (the review's measurement above).
   - Starts (`start`, LAND `{"start": ...}`): 'neutral' (the default)
     as above; 'bare': no cover, no trees, mineral soil (0.37) held;
     'green': cover 1, the trees held at f m v of the record as it
     builds (none until the record admits them, since a part-year record
     reads the season so far) and carbon 13 kg/m² (c = 5 `organicScale`,
     dry soil 0.122) held. The cover runs free in all three, so a bare
     and a green twin keep their own covers through the jumps.
   - Why not a random start. A start drawn per cell from 0–1 has an rms
     error of 0.29 (a cell whose own value is 0.5) to 0.58 (0 or 1)
     against the cell's own value; the trees keep 82 % of
     it after two years where they rise and 51 % where they fall, and
     the noise sits at the grid scale in the albedo and evaporation. The
     placeholders have no error that varies by region, and the jump
     removes what error they have wholesale once a year of the cell's
     own climate is known.
   - Where it runs. scripts/spinup.mjs takes `LAND_JUMPS`: 'fresh' (the
     default) jumps when a record that started fresh passes 365 and 730
     days, at the end of the day that reaches them in whichever segment
     holds it; 'none' never; a list of model days at their ends. A FROM
     run of a state saved before the record never jumps unless asked. Each
     jump is logged with the land means of trees, carbon, dry and full
     land albedo, globally and by 10° band. The page jumps a fresh land
     at the same record ages. A saved negative demand mean (FAO-56's
     reference can be negative under dew) now reloads as saved rather
     than clamped to 0, so a state saved after the jump reloads to the
     bit on the GPU and with its record exact on the CPU.

   Tests. CPU: ten steps crossing 0.9 °C and freezing match the
   arithmetic mean of the steps' own values to 10⁻¹⁴ while within the
   memory and the exponential recursion after; the jump gives 0.1875 of
   trees and 14 kg/m² for a 183-day season at 7.4 °C, P/PET 0.5, cover
   0.8 and record ratio 0.5, and from a sine year's record lands on
   `carbonEquilibrium` to 10⁻¹³. GPU, N=6: a twin whose record is reset
   every step gives each step's own values, and the record matches their
   plain average to 4.2·10⁻⁸ (season length) to 1.5·10⁻⁶ (demand,
   scale 13.5 mm/d) over ten steps, six within the memory; the held
   placeholders of the three starts, the jump against the formula on
   the read-back fields, its repeat and a state file reloaded to the
   bit; after the jump the trees move.

   Run (mechanics, not a result): `jump64`, N=64 GPU, bl34, OCEAN
   `{"everySteps":8}`, ten days from the atlas with `LAND_JUMPS=10`,
   4.0 wall minutes on a shared GPU. A twin run from the same start with
   its record reset every step (a scratch driver; the covers of the two
   stayed identical to the bit, as the held land feeds nothing back)
   gives the model's own per-step values. Land means off the ice sheets
   (area-weighted), the record against the plain average of the
   per-step values and, in brackets, the day's own mean:

   | day | season length | season warmth K | rain mm/d | demand mm/d | litter factor | decay factor |
   |---|---|---|---|---|---|---|
   | 1 | 0.5442 = 0.5442 (0.5442) | 4.4032 = 4.4032 (4.4032) | 0.8824 = 0.8824 (0.8824) | 1.4148 = 1.4148 (1.4148) | 0.1593 = 0.1593 (0.1593) | 0.4781 = 0.4781 (0.4781) |
   | 2 | 0.5399 = 0.5399 (0.5357) | 4.3978 = 4.3978 (4.3924) | 0.9473 = 0.9473 (1.0123) | 1.5173 = 1.5173 (1.6198) | 0.1584 = 0.1584 (0.1574) | 0.4772 = 0.4772 (0.4764) |
   | 5 | 0.5324 = 0.5324 (0.5302) | 4.5579 = 4.5579 (4.9630) | 0.9919 = 0.9919 (1.2893) | 1.7140 = 1.7140 (1.8753) | 0.1598 = 0.1598 (0.1662) | 0.4897 = 0.4897 (0.5227) |
   | 10 | 0.5456 = 0.5456 (0.5819) | 5.4591 = 5.4591 (7.2178) | 1.2191 = 1.2191 (1.3556) | 1.7766 = 1.7766 (1.8797) | 0.1770 = 0.1770 (0.2113) | 0.5723 = 0.5723 (0.7411) |

   The worst cell differs from its own average by 1.3·10⁻⁶ to 2.4·10⁻⁶
   of the field's land maximum on day 10. An exponential mean from an
   empty start would stand at 10/1095 = 0.9 % of these after ten days.
   The jump at the end of day 10 (a ten-day record, March to mid-March:
   the northern land has no season yet), land means including the ice
   sheets: trees 0.226 → 0.037, topsoil carbon 1.61 → 3.30 kg/m², dry
   soil albedo 0.258 → 0.243, land albedo 0.300 → 0.316; by band, trees
   and carbon: 70–60N 0.24 → 0.00 and 1.7 → 0.0, 50–40N 0.25 → 0.00 and
   1.8 → 1.4, 30–20N 0.25 → 0.07 and 1.8 → 6.1, 0–10S 0.26 → 0.07 and
   1.8 → 5.2, 30–40S 0.25 → 0.13 and 1.8 → 6.6 (runs/jump64.log). Repeated at once: trees
   0.0374 → 0.0374, carbon 3.2947 → 3.2947. The saved day-10 state holds
   the record [864000 s, 0, 1].

   What a paired spin-up from the atlas now does (scripts/pairedSpinup.sh,
   the defaults): day 0, the cover 0.5 and free, trees held at half the
   cover, carbon held at 1.802 kg/m², the six means empty at age 0; to
   day 365 the means are the plain average of every step; at the end of
   day 365 (a snapshot day for PER_YEAR 4, saved after the jump) the
   trees jump to f m v and the carbon to its equilibrium from the
   one-year record, and both run free; at the end of day 730 they jump
   again from the two-year record; at day 1095 the record reaches its
   memory and the means become 3-year running means; no later jump.

   What still misses:

   - The year-one record is of a planet whose trees and soils are the
     placeholders, and the second jump weighs that year as much as the
     second; the third year's running mean dilutes it.
   - The cover's own start (0.5) is left to its own times (0.068 of
     cover after two years where the goal is 0).
   - A snow-free cell's trees jump to f m v of its cover at the March
     equinox, wherever its wet and dry seasons leave the cover then
     (five64's Sahel box 0.48 at day 2190 against 0.54, 0.58 and 0.38 at
     days 2008, 2099 and 2281), and the carbon to the year's ratio, not
     its seasonal phase.
   - Land added by regridding during the record takes the estimates,
     counted as a record as old as the rest.
   - The twins ('bare', 'green') are built, not run.

   Review (Oct 2), each check run independently of the work's own
   tests:

   - The CPU record of one cell driven through `update` for three years
     and 30 days matches a compensated float64 average of the cell's
     per-step values, and the exponential recursion after day 1095, to
     4.1·10⁻¹⁴ of the mean at 337.5 s and 168.75 s.
   - The GPU's update m + (x − m) w in a WGSL kernel with the host's
     float64 weights rounded to float32, over seven signals (boreal
     season length and warmth, tropical warmth, rain in showers, demand
     with dew, the litter and tundra decay factors): worst in year one
     2.3·10⁻⁶ to 7.1·10⁻⁵ of the mean's scale; over the three plain
     years 6.2·10⁻⁶ to 1.0·10⁻⁴ at 337.5 s and 1.2·10⁻⁵ to 2.7·10⁻⁴ at
     168.75 s, the showers the worst at both; the kernel equals a
     Math.fround emulation's worst to two digits. The step at the
     hand-over equals its neighbours' (no jump); in the 60 days after it
     the float32 exponential mean errs by up to 4.6·10⁻⁴ (showers,
     337.5 s). The engine itself (N=6, 1200 steps with a 600-step
     memory, against a twin whose record restarts every step): within
     1.0·10⁻⁶ of the scale before the hand-over and 1.3·10⁻⁶ after,
     the placeholders held to the bit on every step.
   - From nine64_day0365 (saved without the means) the CPU engine is
     bit-identical to the work's parent c753e21 over three steps in
     every state and land array; the state loads with the record
     [−1, 0, 0] and the estimates exactly on both engines.
   - N=64 starts on both engines: the same uniform initial fields
     (neutral 0.5, 0.25, 1.802 kg/m²; bare 0, 0, 0; green 1, 0,
     13 kg/m²; bucket 150 mm), the six means 0; over 16 GPU steps the
     carbon is held, the bare trees 0 and the green trees f m v of the
     record to 1.4·10⁻⁷.
   - Runs, N=64 bl34 from the atlas: rvA five days without a jump, rvB
     a jump at the end of day 5 and a second segment to day 10, rvC the
     same in one segment, rvE and rvD five days without and with the
     jump on the snow-free cover. At day 5 all 10640 land cells off the
     ice sheets hold trees 0.5 v and carbon 1.8022 kg/m² to the bit, the
     jump changes only the trees, the carbon and the hold, and both
     engines' jumps of rvE's state equal a hand computation of
     the formula (trees to 1.1·10⁻¹⁶, carbon exactly) and rvD's saved
     state; a repeat changes nothing on either engine; rvD reloads to
     the bit on both engines and through the state file; after it no
     value is negative or non-finite, the carbon is at most 11.38 kg/m²
     and the dry soil spans 0.123–0.370. Four cells of rvE's five-day
     record: warm and wet (20.1°N 75.6°W, P/PET 3.6) trees 0.50747,
     carbon 4.7268 kg/m²; beyond the treeline (76.0°N) 0 and 0; arid
     (11.9°N 125.1°E, P/PET 0.11) 0 and 4.5142; under 29.9 mm of snow
     (23.8°N 105.8°W) 0 and 6.4215. The jump at day 5: trees 0.225 →
     0.025, carbon 1.61 → 3.15 kg/m², dry soil 0.258 → 0.253, land
     albedo 0.279 → 0.296. On day 5 1889 cells with cover lie under
     snow, their snow-free cover above their cover by 1.3·10⁻³ on
     average (3.3·10⁻³ at most). Between days 5 and 10 the trees move
     in 1720 cells by at most 1.4·10⁻³ and the carbon in 6863 cells by
     −0.095 to +0.107 kg/m². rvB and rvC part from day 6 on as rvA's
     second segment and jump64 do without any jump: a whole-day segment
     break is not bit-exact in this spin-up, while the land's own fields
     reload to the bit.

   Defaults: land `start` 'neutral' (`startPlaceholders`: neutral
   cover 0.5, share 0.5, carbon ln 2 × `topsoilMass` × `organicScale`/100
   = 1.802 kg/m²; bare 0, 0, 0; green 1, f m v, 13 kg/m²);
   `FRESH_JUMPS` 365 and 730 days of the record; the spin-up's
   `LAND_JUMPS` 'fresh'; `seasonMemory` and `moistureMemory` 3 years (the
   plain average's length).

   The gases on fixed profiles (Oct 1). `scripts/radiationBenchmark.mjs`
   runs one CPU column on bl34 for each standard atmosphere with the
   reference's own temperature, vapour, ozone and well-mixed gases, clear
   sky, no aerosol, a black surface in the longwave, against: RRTMG's
   tropical, midlatitude summer and winter and subarctic winter examples
   (AER's `run_examples_std_atm`; RRTMG is within 1.5 W/m² of LBLRTM,
   Iacono et al. 2008), the ICRCCM line-by-line fluxes of the five AFGL
   atmospheres (Feigelson et al. 1991, Table 6: CO₂ 300 ppmv, no CH₄ or
   N₂O), LBLRTM's doubled-CO₂ and vapour × 1.2 forcings on the
   midlatitude summer profile (Iacono et al. 2008, RTMIP cases of Collins
   et al. 2006; Mlawer et al. 1997, Table 6, by band) and the CLIRAD
   line-by-line shortwave terms (Chou & Suarez 1999, Tables 7-8, 60°, no
   scattering); profiles and references in `data/radiationBenchmark.json`.
   Misses, W/m² (per cent), before → after:

   | | reference | before | after |
   |---|---|---|---|
   | OLR, TROP / MLS / MLW / SAW (RRTMG) | 287.6 / 281.5 / 230.6 / 199.5 | +5.3 / +3.0 / +9.0 / +12.0 | +0.4 / +0.2 / −0.2 / −0.3 |
   | surface downward longwave (RRTMG) | 398.1 / 348.5 / 224.0 / 172.4 | −65.4 / −39.8 / −2.0 / −2.4 | −1.0 / −0.7 / +2.2 / +2.2 |
   | net longwave at 200 hPa (RRTMG) | | +14.5 / +10.7 / +14.0 / +16.3 | +0.3 / +0.2 / +0.1 / +0.2 |
   | cooling rms below 200 hPa / 3-200 hPa, K/day | | 1.65 / 1.33 / 0.72 / 0.75; 0.82 / 0.86 / 0.45 / 0.34 | 0.28 / 0.24 / 0.12 / 0.09; 0.23 / 0.29 / 0.31 / 0.37 |
   | OLR, AFGL TR / MS / SS / MW / SW (ICRCCM) | 294.0 / 286.6 / 268.1 / 235.9 / 202.6 | −1.8 / −2.0 / −1.3 / +3.6 / +8.8 | −2.0 / −1.2 / −1.2 / −2.3 / −1.2 |
   | surface downward longwave (ICRCCM) | 397.8 / 348.9 / 298.2 / 218.8 / 166.7 | −66.6 / −40.1 / −20.4 / +3.0 / +3.1 | −4.4 / −2.0 / +1.7 / +4.0 / +4.5 |
   | CO₂ 287 → 574 ppmv, TOA / 200 hPa / surface (LBLRTM) | 2.84 / 5.54 / 1.68 | 0 / 0 / 0 | 2.81 / 5.53 / 1.70 |
   | vapour × 1.2 at 574 ppmv (LBLRTM) | 3.79 / 4.52 / 11.55 | 4.10 / 4.12 / 0.38 | 4.07 / 4.90 / 12.55 |
   | OLR slope at fixed RH, MLS / TROP, W/m²/K | about 2 (Koll & Cronin 2018, from memory) | 1.93 / 2.12 | 2.27 / 2.06 |
   | atmosphere's shortwave absorption, μ = 1, albedo 0.2 (RRTMG) | 287.1 / 265.5 / 204.4 / 174.0 | −16.6 / −16.5 / −17.3 / −17.6 % | −0.4 / −0.1 / +0.5 / +0.1 % |
   | the same, MLS at μ = 0.42 (RRTMG, its spectral albedo, 0.213) | 140.8 | −21.8 % | +0.4 % |
   | surface downward shortwave, μ = 1 (RRTMG) | 1053.2 / 1074.4 / 1133.6 / 1162.8 | +3.3 / +2.8 / +1.9 / +1.5 % | −0.7 / −0.8 / −0.8 / −0.7 % |
   | solar heating rms below 200 hPa / 1-200 hPa, K/day | | 0.64 / 0.59 / 0.39 / 0.32; 3.54 / 3.27 / 4.35 / 4.42 | 0.15 / 0.16 / 0.10 / 0.09; 0.47 / 0.48 / 0.42 / 0.39 |
   | MLS, 60°, no scattering: atmosphere; O₂; CO₂ (CLIRAD line-by-line) | 148.2; 4.29; 3.30 | 126.5; 0; 0 | 160.1; 4.17; 3.32 |

   The changes, on both engines (`test/gasRadiation.test.mjs`):
   - Longwave (`longwaveScheme` 'correlated'; 'gray' keeps the three-band
     column): 22 g-points of the simple spectral model of Jeevanjee &
     Fueglistaler (2020) and Williams et al. (2025), exponential envelopes
     of the water-vapour rotation and vibration-rotation bands and CO₂'s
     15 µm band (with its laser and 4.3 µm bands), lines split ±spread
     about the envelope, the self continuum in the spectral shape of
     Roberts et al. (1976), ozone's 9.6 µm band, methane and nitrous oxide,
     line strengths scaled by p/500 hPa, at the diffusivity 1.66;
     `scripts/longwaveFit.mjs` fits its 19 coefficients through the
     g-point reduction to the RRTMG fluxes and cooling rates of the four
     atmospheres, LBLRTM's doubled-CO₂ forcing (total and by band) and the
     minor gases' effects of Chou et al. (2001, Table 16, computed with
     their parameterization at CH₄ 1.75 and N₂O 0.28 ppmv), and writes
     `js/physics/longwaveTable.module.js`; each g-point emits the share of
     σT⁴ a quartic in T gives, normalised to sum to 1. CO₂ 390 ppmv, CH₄
     1.80 ppmv, N₂O 0.323 ppmv (`carbonDioxide`, `methane`,
     `nitrousOxide`; NOAA GML's global annual means for 2010 are 388.8 ppm,
     1798.9 ppb and 323.2 ppb), ozone as below. The cloud's emissivity
     joins every g-point.
   - Shortwave (`solarGases` 'clirad'; 'lacisHansen' keeps the Lacis &
     Hansen vapour and the fixed ozone share), after CLIRAD-SW (Chou &
     Suarez 1999): ozone in the eight bands of their Table 3, vapour by the
     visible band's coefficient and the ten-term near-infrared
     k-distribution of their Table 2 on the path scaled by (p/300 hPa)^0.8
     (1 + 0.00135 (T − 240 K)), every vapour coefficient times 1.48
     (`vaporStrength`, fitted to RRTMG's absorption in the five cases:
     misses at 1 were −23.6, −21.2, −15.5, −13.7 and −9.8 W/m²; the tables
     come from HITRAN-96 lines cut off 10 cm⁻¹ out, with no continuum), O₂
     by their eq. 3.16 over 6.33 % of the beam, CO₂ by a sqrt(u) of its
     scaled path with a giving their Table 7's 3.30 W/m². The light the
     surface reflects crosses the ozone column at 5/3 in the visible and
     the gases' path down plus 5/3 of the column in the near infrared.
   - Ozone (`ozoneColumn` [0.26, 0.35] cm-atm, equator to pole in
     sin²lat; the zonal annual totals from memory) in the vertical
     distribution the old heating used; `ozoneProfile` (CPU, one column)
     replaces it in the benchmark.

   What still misses: the dry atmospheres' surface downward longwave
   (+2.2 against RRTMG, +4.0 and +4.5 against ICRCCM; the tropical −4.4
   against ICRCCM, whose line-by-line codes differ among themselves by up
   to 8 W/m² there, Feigelson et al. 1991, Table 8); the top layer's
   cooling (0-2.2 hPa, MLS −3.7 K/day against −10.3 over the same
   layer, TROP −3.2 against −9.0: the k-distribution with linear pressure
   scaling is not accurate between 0.01 and 10 hPa, where Doppler
   broadening matters, Chou et al. 2001, section 4.2); the vapour × 1.2
   forcing at the surface and at 200 hPa (+8.7 % and +8.4 %); methane's
   and nitrous oxide's own effects on the OLR (MLS 1.65 and 1.18 W/m²
   against the 2.22 and 1.83 of Chou's parameterization), the total being
   in the RRTMG fit; the surface's downward shortwave (−0.8 % at μ = 1
   while the atmosphere's absorption matches); no methane in the
   shortwave (in the fitted vapour strength). The stratosphere-adjusted doubled-CO₂ forcing of the MLS
   column (fixed dynamical heating above 179 hPa, Newton's method to a
   layer heating left of 1e-14 W/m²) is 4.92 W/m² at the top and at
   179 hPa, the top layer cooling 9.5 K.

   Measured with the final defaults, before → after. The clear-sky
   shortwave budget of eight64_day0183 lit over day 186
   (`scripts/clearSkyBudget.mjs`, global, W/m²): reflected 49.5 → 46.9,
   ozone 10.2 → 10.7, vapour 47.3 → 56.7, O₂ and CO₂ 0 → 3.5, aerosol 1.4
   → 1.3 (the atmosphere 58.9 → 72.2; Wild et al. 2019: 73), at the
   surface 232.1 → 221.3 (214); ten64_day0183: reflected 49.3 → 46.8,
   atmosphere 58.7 → 72.0, surface 232.5 → 221.6. Three N=64 GPU days,
   day 186, eight64_day0183 / ten64_day0183: ASR 235.0 → 235.1 / 224.3 →
   224.9, of it in the atmosphere 66.7 → 80.3 / 67.9 → 81.4, OLR 242.2 →
   236.9 / 239.4 → 235.4, clear-sky OLR 259.8 → 256.9 / 258.4 → 255.6,
   clear-sky ASR 291.1 → 293.8 / 291.0 → 293.6, SWCRE −56.1 → −58.7 /
   −66.7 → −68.7, LWCRE 17.6 → 20.0 / 19.0 → 20.2, clear-sky reflectance
   0.1451 → 0.1373 / 0.1454 → 0.1378, rain 1.71 → 1.46 / 2.68 → 2.45
   mm/d, at the sea's surface shortwave 177.8 → 166.6 / 164.2 → 150.3 and
   net longwave −78.2 → −47.6 / −71.8 → −46.5; the day-186 states' clear
   sky (`scripts/clearSkyBudget.mjs`) reflects 49.4 → 46.7 / 49.5 → 46.9,
   absorbs 60.1 → 74.2 / 58.4 → 72.1 in the atmosphere and 231.0 → 219.5 /
   232.6 → 221.4 at the surface. A fresh atlas start on bl34, means of
   days 6-10: Ts 11.24 → 12.50 °C, ASR 187.4 → 198.0, of it in the
   atmosphere 70.1 → 79.6, OLR 221.0 → 209.9, albedo 0.450 → 0.418, SWCRE
   −105.0 → −96.7, LWCRE 36.4 → 40.8, clear-sky reflectance 0.141 →
   0.134, rain 4.52 → 4.39 mm/d, at the sea's surface shortwave 109.0 →
   110.8 and net longwave −75.9 → −63.5. Cost, `js/gpu/profile.module.js`
   over 128 steps from ten64_day0183 at N=64 under the exclusive lock, two
   runs each: the physics pass 2.64 and 2.64 → 3.61 and 3.59 ms, the step
   median 20.79 and 20.79 → 21.75 and 21.81 ms (+4.8 %). The spin-up's
   daily line adds the sea's net surface longwave.

   On the model's own states (one N=64 CPU step, every column alone, and
   the same step on the GPU). eight64_day0183 and nine64_day0091 close in
   every column: absorbed plus reflected is the incoming sunlight to
   5e-13 W/m², the layers' shortwave heating sums to the atmosphere's
   absorption and their longwave heating to σTs⁴ less the downward
   longwave less the OLR to 5e-13 (GPU 5e-4), and dark columns get no
   shortwave. With 'gray' and 'lacisHansen' the CPU gives every column's
   fluxes and heating bit for bit as 5dde347 did. CPU against GPU per
   cell, eight64: clear-sky ASR rms 7.5e-8 and clear-sky OLR 5.8e-7 of
   the field, the layers' longwave heating 1.3e-5; the all-sky misses
   (one cell at 171 W/m² in the shortwave on eight64, cloudy cells in the
   longwave on nine64) are as large with the gray gases. Instantaneous
   CO₂ 390 → 780 ppmv on eight64: OLR −1.08, clear-sky OLR −1.78,
   surface downward +1.89 W/m². The state's stratosphere is warm at the
   top (global means 276, 275 and 264 K at 1.1, 3.6 and 7.5 hPa): the MLS
   column with those temperatures above 110 hPa gives 1.91 W/m² at the
   top instead of 2.74, with 5.60 at 200 hPa in both. OLR slope on
   eight64 with Ts and every layer below 150 hPa ±1 K at fixed relative
   humidity: 2.17 W/m²/K, clear sky 1.89. In the three days from
   eight64_day0183 the global mean of the top layer (1.1 hPa) warms
   275.6 → 290.7 K with the spectral gases (275.6 with the gray ones),
   the 3.6 and 7.5 hPa layers cool 275.3 → 260.2 and 264.3 → 248.3 K:
   the top layer's shortwave heating is RRTMG's (MLS overhead 22.9
   against 23.9 K/day) while its cooling is a third of it.

   The upper layers, the stratosphere and the sourced amounts (Oct 2).
   The top layer's cooling. RRTMG's cooling over the 0-2.2 hPa layer is
   −9.04 / −10.25 / −8.74 / −7.13 K/day (TROP / MLS / MLW / SAW); line by
   line, the largest stratospheric cooling of the midlatitude summer and
   subarctic winter columns is about −12 °C/day (Chou et al. 2001,
   section 10), and RRTMG's stratospheric heating is generally within
   0.4 K/day of LBLRTM's (Iacono et al. 2008). By RRTMG's sixteen bands (MLS, the spectral model above 200 hPa)
   the deficit was CO₂'s 630-700 cm⁻¹ centre (top layer −1.1 against
   −5.7 K/day), its wings (500-630 and 700-820 cm⁻¹: −0.05 against −1.6)
   and the vapour's rotation band (10-500 cm⁻¹: 0 against −0.8), while
   ozone's 980-1080 cm⁻¹ band overcooled the 3.6 and 7.5 hPa layers
   (−3.4 and −1.4 against −1.35 and −0.63, its paths scaled by
   (p/500 hPa)^0.17). The changes (`longwave.module.js`,
   `scripts/longwaveFit.mjs`, both engines): the vapour lines', CO₂'s and
   ozone's paths scale with the Voigt pressure √(p² + p_D²), p_D the
   pressure at which the gas's Lorentz half-width equals its Doppler
   half-width at 250 K (air-broadened half-widths 0.09, 0.07 and
   0.075 cm⁻¹/atm at 296 K, typical of HITRAN's lines at 300, 667 and
   1042 cm⁻¹, from memory: 400, 727 and 1012 Pa, not fitted); the
   spectral model carries each absorber's line cores, the k(g) ∝ g⁻²
   tail of a Lorentz line's k-distribution from the upper half's k at g₀
   down to g₀ e^(−depth) in four nodes (fitted: CO₂ g₀ 0.152, depth
   2.37, on an envelope of e-folding 16.5 cm⁻¹; ozone 0.504, 16.9; vapour
   below 500 cm⁻¹ 0.012, 8.2), ozone's pressure exponent refits to 0.94;
   the reduction bins the cores into g-points up to a
   midlatitude-summer column depth of 10⁹ (22 → 34 g-points), each
   g-point's Planck share following its members' widths; the fit scores
   the 3-30 hPa layers and the top layer of the four atmospheres, each
   RRTMG band's cooling above 200 hPa, the dry atmospheres' surface
   downward flux (×3), LBLRTM's doubled CO₂ (×150) and methane and
   nitrous oxide from none to their 1860 amounts (Iacono et al. 2008,
   Table 2, case 3a-1a: 3.60 / 3.45 / 1.08 W/m² at the top, 200 hPa and
   the surface). A Doppler width alone (refit, no cores) brought the top
   layer to −7.3 / −8.1 / −7.1 / −5.4 K/day but the doubled-CO₂ forcing
   at the top to 2.43 W/m² (−14 %), as a floor had.
   Cooling above 30 hPa, K/day, before (1f61d5c) → after, RRTMG:

   | | 1.1 hPa | 3.6 hPa | 7.5 hPa | 14 hPa | 24 hPa |
   |---|---|---|---|---|---|
   | TROP | −3.20 → −9.25 (−9.04) | −5.97 → −6.13 (−6.00) | −4.11 → −3.81 (−3.92) | −2.22 → −2.62 (−2.49) | −1.07 → −1.67 (−1.56) |
   | MLS | −3.68 → −10.12 (−10.25) | −6.92 → −6.67 (−6.42) | −4.49 → −3.91 (−4.06) | −2.17 → −2.53 (−2.44) | −1.03 → −1.60 (−1.53) |
   | MLW | −3.22 → −9.17 (−8.74) | −3.17 → −3.85 (−3.86) | −1.58 → −1.83 (−2.03) | −1.09 → −1.43 (−1.49) | −0.93 → −1.26 (−1.28) |
   | SAW | −2.29 → −7.24 (−7.13) | −1.96 → −2.68 (−2.90) | −1.35 → −1.70 (−1.94) | −1.11 → −1.43 (−1.53) | −0.82 → −1.06 (−1.14) |

   The top layer is within 5 % in all four, every 3-30 hPa layer within
   0.25 K/day (0.94 before). By band (MLS, cooling rms above 200 hPa)
   10-350 / 500-630 / 630-700 / 700-820 / 980-1080 cm⁻¹: 0.26 / 0.21 /
   1.25 / 0.32 / 0.61 → 0.26 / 0.09 / 0.32 / 0.17 / 0.19 K/day. The rest
   of the benchmark (`scripts/radiationBenchmark.mjs`), before → after:
   OLR +0.4 → +0.8 / +0.2 → +0.8 / −0.2 → −0.8 / −0.3 → −1.2 W/m²
   against RRTMG; surface downward −1.0 → −1.1 / −0.7 → −0.7 / +2.2 →
   +1.4 / +2.2 → +2.2; net at 200 hPa +0.3 → +0.8 / +0.2 → +1.1 / +0.1
   → +0.1 / +0.2 → +0.2; cooling rms below 200 hPa 0.28 → 0.31 / 0.24 →
   0.27 / 0.12 → 0.12 / 0.09 → 0.09, 3-200 hPa 0.23 → 0.08 / 0.29 →
   0.10 / 0.31 → 0.09 / 0.37 → 0.12 K/day; ICRCCM OLR −2.0 → +0.8 / −1.2
   → +0.9 / −1.2 → +0.4 / −2.3 → −1.5 / −1.2 → −1.0, surface downward
   −4.4 → −5.3 / −2.0 → −2.6 / +1.7 → +0.2 / +4.0 → +2.0 / +4.5 → +3.1;
   LBLRTM MLS OLR 281.8 → 282.3 (283.3); doubled CO₂ 2.81 → 2.63 /
   5.53 → 5.72 / 1.70 → 1.65 (2.84 / 5.54 / 1.68); vapour × 1.2 4.07 →
   3.97 / 4.90 → 4.67 / 12.55 → 11.99 (3.79 / 4.52 / 11.55); methane and
   nitrous oxide from none 2.02 → 3.19 / 1.91 → 2.99 / 0.76 → 1.26
   (3.60 / 3.45 / 1.08); without methane, MLS, OLR / surface 1.65 → 2.37
   / −0.55 → −0.79 (Chou's 2.22 / −0.86), without nitrous oxide 1.18 →
   1.91 / −0.51 → −0.87 (1.83 / −0.58); the stratosphere-adjusted doubled
   CO₂ 4.92 → 5.49 W/m² at the top and at 179 hPa (the top layer cools
   9.5 → 19.3 K, the 3.6 and 7.5 hPa layers 8.6 → 14.7 and 7.6 → 11.4 K;
   the earlier solver stopped with up to 1.5 W/m² of layer heating left
   and gave 4.22 → 5.07 at the top); OLR
   slope 2.27 → 2.36 (MLS), 2.06 → 2.14 W/m²/K (TROP). The GPU longwave
   heating of the layers above 30 hPa matches the CPU's to 1e-6 rms
   (`test/gasRadiation.test.mjs`).

   What one 0-2.2 hPa layer can and cannot hold. It holds 0.2 % of the
   atmosphere's mass, everything above about 42 km, and its temperature
   is that mass's mean: in RRTMG's tropical column the air warms from
   263 K at 2.1 hPa to 269 K near 0.9 hPa (the stratopause) and cools to
   229 K at 0.1 hPa, half the layer's mass lies below 1.1 hPa, and its
   cooling runs from −9.0 through −11.3 (1.07 hPa) to −3.6 K/day at
   0.14 hPa, its solar heating (MLS, overhead) from 28.5 to 10 K/day.
   The fitted layer gives the mass-weighted mean of both at the layer's
   mean temperature; it cannot hold the stratopause's maximum, the
   mesosphere's fall of temperature above it, or the difference between
   the two. Two measurable costs: on doubling CO₂ the band centre's
   emission rises within one isothermal layer instead of into the colder
   mesosphere, so the 630-700 cm⁻¹ forcing at the top is −0.74 W/m²
   against LBLRTM's −0.57 (the total 2.63 against 2.84), and the
   stratosphere-adjusted forcing lets the whole layer cool 19.3 K. Report
   only: two more interfaces, at about 0.3 and 1 hPa, would put the
   stratopause between layer midpoints and give the band centre a cold
   mesospheric layer to emit from, at 36 layers (+6 % of the column
   physics); the sponge (`topSigma` 0.02) would then cover them.

   Sourced amounts. The model's present day is 2010: CO₂ 388.75 ppm,
   CH₄ 1798.93 ppb, N₂O 323.18 ppb, NOAA GML's global annual means
   (`co2_annmean_gl.txt`, `ch4_annmean_gl.txt`, `n2o_annmean_gl.txt`,
   gml.noaa.gov) for 2010, in `GREENHOUSE_GASES`. Ozone (`ozone` 'afgl',
   `js/physics/ozone.module.js`, table by `scripts/ozoneTable.mjs`): the
   five AFGL atmospheres' profiles (Anderson et al. 1986, AFGL-TR-86-0110,
   as libRadtran distributes them, the same in
   `data/radiationBenchmark.json`) as the column above each pressure,
   tropical equatorward of 15°, midlatitude at 45°, subarctic poleward of
   60°, linear in latitude between, each hemisphere's summer profile on
   15 July and winter one on 15 January with the cosine of the time of
   year between; each layer takes the column between its interfaces'
   pressures, so the column over high ground is shorter. The columns are
   282 DU (tropical), 334 / 378 (45° summer / winter), 348 / 376 (60°),
   320.6 DU in the global mean over a 1013 hPa surface at any date (the
   sin² column of [0.26, 0.35] gave 290). Against the idealised shape the
   tropics hold 30 % more ozone at 3-25 hPa (5.9 / 17.2 / 38.3 / 56.0 /
   55.6 DU in the five layers above 30 hPa, 7.2 / 13.7 / 28.3 / 42.4 /
   44.8 before) and a third as much below 85 hPa; at 60° less above
   30 hPa and twice as much below 100 hPa. The AFGL profiles are of the
   1960s-70s: WMO (2010, chapter 2) puts the 2006-2009 global total
   ozone 3.5 % below its 1964-1980 mean, and the Antarctic spring
   depletion is absent; no zonal-mean climatology of the 2010s was
   reachable from here (McPeters & Labow 2012 and the SBUV merged set
   need a login), so the AFGL set stays, with a seasonal phase at the
   solstices' month where the observed columns peak in spring.
   CLIRAD-SW's near-infrared ozone (band 9) is already in: Chou & Suarez
   (1999, section 3.6, eqs. 3.14-3.15) fold it into band 8's coefficient,
   0.0572 including Δk 0.0032-0.0033, and the model uses that
   coefficient over band 8's share. CLIRAD-SW has no methane, so the
   shortwave methane stays inside the vapour strength fitted to RRTMG.
   Methane's and nitrous oxide's longwave are now fitted to LBLRTM's
   case 3a-1a (above); the split between them follows Chou et al.'s
   (2001) Table 16 (their parameterisation, not line by line).
   The visible light that cloud and clear air reflect now crosses the
   ozone column on its way out at 5/3, as the surface's did, losing
   1 − exp(−0.0542 Ω 5/3) to the layers by their ozone (both engines);
   the light cloud and clear air reflect had left unabsorbed. On the
   RRTMG cases at μ = 1 the overhead solar heating at 14 / 24 / 38 /
   54 hPa rises 4.55 → 4.79 / 2.75 → 2.93 / 1.68 → 1.80 / 1.09 → 1.17
   K/day (MLS; RRTMG 5.34 / 3.42 / 2.17 / 1.48) and the atmosphere's
   absorption −0.7 → +1.0 W/m² (+0.4 %).

   The visible / near-infrared split. `visibleFraction` (VISIBLE_FRACTION)
   0.5 → 0.4707, the share of CLIRAD-SW's bands 1-8 (Chou & Suarez 1999,
   Table 3), so that the near infrared the gases absorb from is the
   0.5293 their k-distribution's weights are shares of. The Rayleigh
   sub-bands refit (`node scripts/rayleighReference.mjs`, whose VISIBLE
   is now VISIBLE_FRACTION) to the 0.297-0.683 µm band of a 5778 K
   spectrum this leaves: [[0.7049, 0.0957], [0.2951, 0.5806]], within
   0.1 % of the band's reference reflectance over μ 0.1-1 and 23.21 W/m²
   over a black surface against its 23.20 (the 0.5 band's 23.54 counted
   the 0.683-0.711 µm sliver that now falls in the near infrared). The
   near infrared scatters too (`nearInfraredRayleigh` 0.0114:
   `NEAR_INFRARED=1` at 400 wavelengths, where the 40 of the visible band
   give 1.95 W/m² and 0.0112, fits one grey depth to the rest of the
   spectrum, 2.00 W/m² against the reference's 2.00, within 2.2 % over
   μ 0.1-1; CLIRAD's bands 9-10 give 0.0101 weighted by their shares):
   the black-surface reflection is 25.2 W/m², the reference's whole
   spectrum (25.20 two-stream, 25.02 doubling-adding),
   where the 0.5 band alone gave 23.5 and nothing beyond 0.711 µm. On the
   day-186 state of three days from eight64_day0183 the open sea's
   clear-sky albedo at the top by class (`scripts/clearSkyBudget.mjs`)
   is 0.091 / 0.110 / 0.148 for 0-30 / 30-50 / 50-70° (0.091 / 0.111 /
   0.150 with 1f61d5c on its own day 186; ranges 0.08-0.10 / 0.10-0.13 /
   0.13-0.20).

   The stratosphere's own balance. Thirty N=64 GPU days with the final
   defaults from copies of eight64_day0183 (to day 213, September-October)
   and nine64_day0091 (to day 121, June-July), the spin-up logging every
   layer above 200 hPa by zone (`STRATOSPHERE=1`: the layer-mean
   temperature θ times the layer's own Exner function, the temperature
   the radiation sees, 5 % below θ (p_mid/p₀)^κ in the 0-2.2 hPa layer);
   heating from one CPU step at six times of the day on the day-30 state.
   Layer-mean temperature (K) on day 213 / day 121, global, 20S-20N,
   35-55N, 35-55S, 70-90N, 70-90S, against the AFGL atmospheres' layer
   means (tropical; midlatitude summer, winter; subarctic summer, winter):

   | layer | 1f61d5c, day 213 global | global | tropics | 35-55N | 35-55S | 70-90N | 70-90S | AFGL TR; MS, MW; SS, SW |
   |---|---|---|---|---|---|---|---|---|
   | 1.1 hPa | 283 / 282 | 252 / 252 | 254 / 254 | 248 / 259 | 258 / 245 | 224 / 268 | 258 / 215 | 261; 266, 256; 269, 247 |
   | 3.5 hPa | 245 / 245 | 243 / 244 | 245 / 245 | 239 / 250 | 246 / 237 | 227 / 255 | 251 / 227 | 252; 257, 236; 261, 227 |
   | 7.4 hPa | 231 / 231 | 230 / 231 | 231 / 231 | 226 / 233 | 231 / 225 | 222 / 241 | 242 / 223 | 240; 243, 222; 245, 218 |
   | 14 hPa | 223 / 224 | 220 / 221 | 220 / 221 | 217 / 221 | 221 / 217 | 219 / 231 | 235 / 219 | 230; 233, 216; 235, 214 |
   | 24 hPa | 218 / 218 | 214 / 215 | 213 / 214 | 212 / 215 | 216 / 213 | 216 / 228 | 231 / 213 | 222; 226, 215; 229, 212 |
   | 37 hPa | 211 / 212 | 210 / 211 | 207 / 208 | 210 / 211 | 214 / 212 | 217 / 227 | 227 / 208 | 216; 223, 215; 226, 213 |
   | 53 hPa | 203 / 206 | 207 / 208 | 202 / 204 | 208 / 208 | 211 / 212 | 217 / 225 | 223 / 206 | 208; 220, 215; 225, 214 |
   | 70 hPa | 200 / 201 | 205 / 205 | 199 / 200 | 206 / 207 | 209 / 212 | 216 / 223 | 219 / 205 | 201; 218, 216; 225, 216 |
   | 85 hPa | 199 / 200 | 203 / 204 | 198 / 198 | 205 / 206 | 207 / 211 | 215 / 220 | 215 / 204 | 197; 216, 216; 225, 216 |
   | 101 hPa | 200 / 200 | 203 / 203 | 198 / 197 | 205 / 206 | 207 / 209 | 213 / 217 | 212 / 203 | 196; 216, 217; 225, 217 |

   Global net radiative heating on the day-30 states is −0.01 to −0.17
   K/day in the layers from 3.5 to 101 hPa, −0.17 / −0.34 in the top
   one (solar / longwave, day 213: 8.24 / −8.41 at 1.1 hPa, 5.35 /
   −5.48, 2.95 / −3.03, 1.80 / −1.92, 1.15 / −1.25, 0.70 / −0.80 and
   0.50 / −0.59 K/day down to 53 hPa), and the
   3.5-85 hPa layers still cool 0.03-0.18 K/day in the global mean and
   the tropical 70-101 hPa layers 0.15-0.30 K/day (days 21-30): they
   settle colder still. Against the climatology (the targets of the review, from the US
   Standard Atmosphere 1976, CIRA-86 and Seidel et al. 2001, not
   re-checked here): the tropical cold point is 197-198 K at 85-101 hPa
   (190-195 K near 90-100 hPa, warm by 3-7 K); midlatitudes at 53 hPa
   208-212 K (215-220, cold by 5-10); about 10 hPa, between the 7.4 and
   14 hPa layers, about 222-227 K (228-235, cold by about 6); the 0-2.2 hPa
   layer 245-259 K at midlatitudes (the stratopause near 1 hPa 265-270 K;
   AFGL's layer means 256-266 K, cold by about 8); the winter pole
   (70-90S, July) 203-227 K at 1-101 hPa with no cold vortex, where the
   vortex is 185-195 K. Which terms:
   - The winter and summer poles are the sponge's. With the Rayleigh drag
     on σ < 0.02 off (`SURFACE='{"topDragDays":0}'`, nine64_day0091 to
     day 121, otherwise the same) the 70-90S column is 181 / 190 / 187 /
     187 / 187 / 188 / 190 / 192 / 193 / 193 K from 1.1 to 101 hPa, a
     vortex, its radiative cooling down to −0.3 to −1.2 K/day from −0.4 to
     −3.8 with the sponge: the drag on the polar-night jet drives a
     circulation that warms the winter pole by 10-37 K and supplies up to
     3.8 K/day there. The summer pole warms 7-10 K above 7 hPa without
     it (275 / 265 / 250 K at 1.1 / 3.5 / 7.4 hPa against 268 / 255 / 241);
     the global and tropical means move by 3 K or less, the 30-day means
     of balance, OLR and rain by under 0.3. No change made: the sponge is
     the dynamics', for the core's own review (The model top, below).
   - The cold middle stratosphere is in part a radiation term: on the
     standard columns at an overhead sun the solar heating is below
     RRTMG's (MLS: −3.6 % at 1.1 hPa, −4.0 / −10.3 / −14.3 / −17.1 /
     −20.9 / −25 % at 7.5 / 14 / 24 / 38 / 54 / 71 hPa, after the
     reflected light's ozone above, which lifted 14-54 hPa by 5-7 %),
     while the longwave matches; at each layer's global longwave and the
     cooling-to-space sensitivity c₂ν/T² of the 667 cm⁻¹ band, that
     deficit is worth about 2 / 2 / 5 / 7 / 8 K at 1.1 / 7.4 / 14 / 24 /
     37 hPa. RRTMG's shortwave reference here has no band or gas split
     to say which term misses (CLIRAD's single ozone coefficient per
     ultraviolet band, or the near-infrared CO₂, vapour and O₂ heating,
     0.15 K/day together in the model at 7-100 hPa); not corrected
     further.
   - The ozone amount and shape: on the day-213 state the AFGL ozone
     heats the layers 8.24 / 5.35 / 2.95 / 1.80 K/day in the global mean
     at 1.1 / 3.5 / 7.4 / 14 hPa where the idealized shape of 1f61d5c
     gave 9.70 / 4.54 / 2.78 / 1.81: the top layer holds less ozone (5.9
     DU at the equator against 7.2), worth about 11 K at its longwave
     sensitivity, as much as its departure from the AFGL layer means. The
     top layer's thickness: one 0-2.2 hPa layer cannot be the stratopause
     (above). Before (1f61d5c) the top layer settled at 282-283 K, 30 K
     above the AFGL layer mean, its cooling a third of RRTMG's.
   Instantaneous CO₂ ×2 (388.75 → 777.5 ppmv) on the day-30 states: OLR
   −2.04 / −2.07, clear-sky −2.73 / −2.77 W/m² (the starting states,
   whose stratosphere 1f61d5c warmed: −0.93 / −0.90, clear −1.68 /
   −1.68; 1f61d5c on its own day-30 states, 390 → 780: −1.84 / −1.86,
   clear −2.55 / −2.57).

   Measured with the final defaults, before (1f61d5c) → after. The
   shortwave benchmark: the atmosphere's absorption at μ = 1 −1.2 → +0.1
   / −0.4 → +1.0 / +1.0 → +2.6 / +0.1 → +1.8 W/m² (TROP / MLS / MLW /
   SAW; +0.5 → +1.7 at μ = 0.42), the surface's downward −0.7 → −0.8 /
   −0.8 → −0.9 / −0.8 → −0.9 / −0.7 → −0.8 / −0.7 → −1.1 %, solar heating
   rms below 200 hPa 0.15 → 0.16 / 0.16 → 0.17 / 0.10 → 0.10 / 0.09 →
   0.09 / 0.09 → 0.09 and 1-200 hPa 0.47 → 0.37 / 0.48 → 0.37 / 0.42 →
   0.32 / 0.39 → 0.29 / 0.27 → 0.18 K/day; MLS overhead ozone 36.5 →
   38.0, the reflected light's absorption on its way up 14.8 → 16.1 W/m².
   Three N=64 GPU days from a copy of eight64_day0183, day 186: ASR 235.1
   → 236.8, of it in the atmosphere 80.3 → 82.3, OLR 236.9 → 237.5 W/m²,
   albedo 0.310 → 0.304, SWCRE −58.7 → −57.2, LWCRE 20.0 → 20.1, clear-sky
   reflectance 0.1373 → 0.1365, rain 1.46 → 1.46 mm/d, at the sea's
   surface shortwave 166.6 → 166.5 and net longwave −47.6 → −47.9 W/m²;
   the top layer 273.9 → 259.2 K in the global mean (282.7 → 252.1 by
   day 213). The day-186 states' clear sky (`scripts/clearSkyBudget.mjs`,
   global): reflected 46.7 → 46.5, absorbed in the atmosphere 74.2 →
   75.4 (ozone 10.7 → 12.0, vapour 58.7, O₂ and CO₂ 3.5, aerosol 1.3 →
   1.2), at the surface 219.5 → 218.7 W/m² (Wild et al. 2019: 53 / 73 /
   214); the black-surface reflection 24.6 → 24.5. Thirty days, means of
   days 1-10 and 21-30, from eight64_day0183 (1f61d5c → final): balance
   ASR − OLR −1.3 → −0.2 and −9.3 → −5.5, OLR 233.5 → 234.2 and 228.1 →
   228.5 W/m², rain 1.99 → 1.98 and 2.68 → 2.60 mm/d, the sea's surface
   shortwave 163.5 → 163.1 and 150.0 → 153.9, net longwave −47.8 → −48.0
   and −44.5 → −45.7 W/m², albedo 0.318 → 0.313 and 0.358 → 0.345; from
   nine64_day0091: balance −7.4 → −6.2 and −13.7 → −11.7, OLR 233.6 →
   234.2 and 232.7 → 232.5, rain 2.17 → 2.16 and 2.62 → 2.70, sea
   shortwave 138.0 → 137.9 and 126.7 → 125.8, net longwave −46.4 → −46.5
   and −40.0 → −40.4, albedo 0.336 → 0.330 and 0.357 → 0.351. Both
   versions drift the same way over the month (albedo +0.02-0.04 and
   rain +0.4-0.7 mm/d from days 1-10 to 21-30, with or without these
   changes); the radiation moves the balance by +1 to +4 W/m², mostly
   through the shortwave (ASR +1.7 to +4.2 W/m² by days 21-30). Cost,
   `js/gpu/profile.module.js` over 128 steps under the exclusive lock, two
   runs each: N=64 (ten64_day0183) step median 21.75 and 21.71 → 22.47
   and 22.48 ms (+3.4 %), the physics pass 3.58 → 4.34 ms; N=128
   (eight128_day0183) 92.88 and 92.79 → 95.91 and 95.92 ms (+3.3 %), the
   physics pass 15.11 → 18.38 ms, 49 s a model day.

   Final defaults: `longwaveScheme` 'correlated' with the 34 g-points of
   `js/physics/longwaveTable.module.js` (p_D 400 / 727 / 1012 Pa),
   `carbonDioxide` 388.75e-6, `methane` 1798.93e-9, `nitrousOxide`
   323.18e-9, `ozone` 'afgl', `visibleFraction` 0.4707, `rayleighBands`
   [[0.7049, 0.0957], [0.2951, 0.5806]], `nearInfraredRayleigh` 0.0114,
   the rest as before. What still misses: the water vapour's rotation
   band cools the stratosphere too little by itself (10-350 cm⁻¹, top
   layer −0.05 against −0.63 K/day, −0.03 against −0.3 to −0.4 at
   3-40 hPa in the spectral model), the CO₂ and ozone g-points making up
   the total; by band the spectral model's top layer cools −6.9 against
   −10.25 K/day, the g-points' correlated overlap the rest; the stratosphere's
   solar heating below RRTMG's at 7-100 hPa (−4 to −25 %); the doubled-CO₂
   forcing at the top −7 % (2.63 against 2.84; its 630-700 cm⁻¹ band
   −0.74 against −0.57, the single top layer), the stratosphere-adjusted
   one 5.49 W/m² with a 19.3 K cooling of the top layer (9.5 K with
   1f61d5c); methane's and
   nitrous oxide's own longwave −11 % at the top and −13 % at 200 hPa
   against LBLRTM, +17 % at the surface, and without either gas the
   subarctic winter's OLR change 1.82 / 1.27 against Chou's 1.00 / 1.10;
   the vapour × 1.2 forcing +5 / +3 / +4 %; the ICRCCM tropical surface
   downward longwave −5.3 W/m²; the AFGL ozone a 1960s-70s climatology,
   3.5 % above the 2006-2009 global total (WMO 2010), with no Antarctic
   spring depletion and its seasons at the solstices' months; no
   shortwave methane or CFCs; the 23.5 W/m² visible-band Rayleigh
   reflection of the review is now 25.2 for the whole spectrum (the
   visible band's alone 23.2, −1.4 %); the sponge's warming of both polar
   stratospheres (the dynamics); three GPU parity tests now fail at one
   column each where an f32 rounding tips a threshold: the cloud effects
   after 16 steps (`test/cloudEffect.test.mjs`; at step 14 one column's
   boundary layer mixes layers 22-25 on one engine only, 0.3-0.7 g/kg of
   vapour apart, its longwave cloud effect then 1.9 W/m² apart against
   0.5 allowed), the cloudy columns' layer heating after a 10-day
   mixed-layer step (`test/gpuModel.test.mjs`, 1.3e-4 against 1e-4 K/day
   at a deck's top, its cover 1.6e-6 apart; 9.1e-5 with 1f61d5c) and the
   regime-gated deck (one decoupled column's water 0.73 against
   0.61 kg/m², a layer at the 1.9 km deck top counted on one engine);
   the per-layer longwave heating of a step agrees as before (1e-5 rms),
   and the tolerances are unchanged.

   Review (Oct 2), checked independently of the scripts above.
   - RRTMG's own layer heating averaged by mass over the model's layers
     (not the scripts' interpolation of its fluxes): the 0-2.2 hPa layer
     −8.97 / −10.54 / −8.74 / −7.42 K/day over the part the reference
     covers (its MLS and SAW columns stop at 6.7 and 10 Pa, 97.0 and
     95.5 % of the layer's mass, which `referenceHeating` counts as not
     cooling: −10.25 and −7.13); the model's −9.25 / −10.12 / −9.17 /
     −7.24 are within 5 % of both; the 3-30 hPa layers within 0.25 K/day
     of both.
   - Cool-to-space of the 0-2.2 hPa layer by line-by-line CO₂ 15 µm
     lines built from band constants (the 626 fundamental and the
     02201, 10002 and 10001 hot bands and the 636 fundamental, strengths
     at 296 K 7.97e-18, 6.2e-19, 1.4e-19, 1.9e-19 and 8.3e-20 cm/molecule
     from memory; rigid rotor, Hönl-London factors, Voigt lines with an
     air-broadened half-width of 0.07 cm⁻¹/atm at the Curtis-Godson
     1.1 hPa, the layer's mean temperature, diffusivity 1.66): −4.3 /
     −4.7 / −4.0 / −3.3 K/day (TROP / MLS / MLW / SAW), −3.2 / −3.5 /
     −3.0 / −2.5 of it at 630-700 cm⁻¹, and −1.8 / −2.0 / −1.6 / −1.4
     from the Doppler cores alone. Same sign as RRTMG's CO₂ bands (MLS
     −5.7 in the centre, −1.6 in the wings) and two thirds of their size,
     the rest within the bands left out and the exchange with the colder
     layers below; at 1.1 hPa the Lorentz wings of the Voigt lines carry
     more than half of the layer's CO₂ cooling.
   - The spectral model by itself (1 cm⁻¹) misses RRTMG's OLR by +8.5 /
     +8.1 / +4.8 / +3.0 W/m² (1f61d5c: +2.6 / +2.2 / +0.9 / +0.2), the
     net flux at 200 hPa by +10.4 / +10.5 / +7.4 / +6.1, the 3-30 hPa
     cooling by 0.99 / 1.02 / 0.70 / 0.68 K/day rms, cools the top layer
     −6.3 / −6.9 / −6.3 / −4.9 K/day, and gives methane and nitrous oxide
     from none 5.82 / 5.34 / 3.35 W/m² (LBLRTM 3.60 / 3.45 / 1.08): the
     34 g-points meet RRTMG through their reduction (binning by the
     midlatitude summer column's depth, each g-point's absorbers
     correlated), which supplies about 7 W/m² of the OLR and 3 K/day of
     the top layer. Two refits that also score the spectral model's own
     misses (weights 1 and 0.3 of the g-points', 3000 iterations each)
     bring its OLR to +2.0 / +2.8 (MLS) and its top layer within 10 %, but
     move the g-points' top layer to −11.5 / −11.1 K/day (MLS), methane
     and nitrous oxide to 1.98 / 2.58 W/m² at the top, doubled CO₂ at the
     surface to 1.33 / 1.44 (1.68) and the ICRCCM subarctic winter's
     surface downward longwave to +6.3 / +5.2: not adopted.
   - Closure on eight64_day0183, the CPU every column alone at four
     times of the day: absorbed plus reflected less incoming at most
     4.5e-13 W/m², the layers' shortwave less the atmosphere's 3.4e-13,
     the layers' longwave less σTs⁴ − DLR − OLR 9.1e-13, the dark
     columns' shortwave 2.8e-14; one GPU step (also nine64_day0091):
     1.8e-4, longwave 2.0e-4 (2.3e-4) W/m², dark columns 0. CPU against
     GPU after one step, longwave heating of the 1.1 / 3.5 / 7.4 hPa
     layers: eight64 rms 1.1e-5 / 1.5e-5 / 1.1e-5 K/day, at most 4.7e-5
     / 5.9e-5 / 4.4e-5 (1f61d5c 0.8e-5 / 2.8e-5 / 2.3e-5, at most 3.7e-5
     / 1.1e-4 / 9.4e-5); nine64 rms 2.1e-5 / 4.7e-5 / 6.8e-5, at most
     7.7e-4 / 1.9e-3 / 3.0e-3 in one column (1f61d5c 1.6e-4 / 3.0e-4 /
     4.0e-4, at most 7.6e-3 / 1.5e-2 / 2.0e-2). The all-sky per-cell
     shortwave misses (one cell 168 W/m², rms 4.6, eight64) are
     1f61d5c's (171, 4.7).
   - With the gas options at 1f61d5c's values (`visibleFraction` 0.5,
     its `rayleighBands`, `nearInfraredRayleigh` 0, `ozone` 'idealized',
     CO₂ 390, CH₄ 1.8, N₂O 0.323 ppmv) and 'gray' with 'lacisHansen', or
     'gray' with 'clirad' and no upward absorption, every column's
     budget and layer heating on eight64_day0183 is 1f61d5c's bit for
     bit on the CPU; on the GPU the OLR and the deck are, and the
     shortwave fields differ by one or two ulps in 2 % of the cells.
   - The failing parity tests: the cloudy columns' heating passes with
     `ozone` 'idealized' and fails with 'afgl' whatever the visible split
     and near-infrared scattering, its per-layer rms 7.6e-6 against
     8.1e-6 K/day and its worst cell a deck top in both (6.9e-5 against
     1.28e-4); the regime-gated deck fails at `visibleFraction` 0.4707 and
     passes at 0.4708 (and at `nearInfraredRayleigh` 0): one-column
     threshold flips.
   - The 30-day logs re-read (`STRATOSPHERE=1` lines): the table above
     holds. Trends over days 21-30 (least squares): the global layers at
     3.5-85 hPa −0.11 to −0.18 K/day from eight64_day0183 and −0.03 to
     −0.14 from nine64_day0091; the tropical 70-101 hPa layers −0.15 to
     −0.30 K/day in both, so the cold point of 197-198 K is still
     falling; at the poles the trends follow the season (70-90S +0.4 to
     +0.6 K/day at 24-85 hPa in September-October), and without the
     sponge the 70-90S column still cools 0.6-1.2 K/day at day 121.
   - Sources: NOAA GML's `co2_annmean_gl.txt`, `ch4_annmean_gl.txt` and
     `n2o_annmean_gl.txt` (read Oct 2) give 388.75 ppm, 1798.93 ppb and
     323.18 ppb for 2010; the electronic Tables 1a-1e of AFGL-TR-86-0110
     (github.com/rayference/afgl1986) integrate to ozone columns of 282.0
     / 334.4 / 378.3 / 347.8 / 375.5 DU, the table's 281.9 / 334.2 /
     377.9 / 347.6 / 375.8; p_D from the stated half-widths is 395 / 726
     / 1012 Pa.
   - Reproduced: the benchmark's numbers above; three days from a copy of
     eight64_day0183 give the day-186 line of the run above bit for bit
     and the open sea's classes 0.091 / 0.110 / 0.148; the step medians
     21.75 / 21.74 → 22.45 / 22.42 ms at N=64 and 92.61 → 95.95 ms at
     N=128 (49.1 s a model day).

   The model top (Oct 2). Thirty N=64 GPU days from copies of
   nine64_day0091 (to day 121, July) and nine64_day0274 (to day 304,
   January), one at a time, `STRATOSPHERE=1` (with the second daily line
   of the upper layers' winds) and `scripts/upperAtmosphere.mjs` on the
   day-30 states.
   - What the Rayleigh sponge prevents (σ < 0.02, 5 days at the top:
     rates 0.95 / 0.82 / 0.63 / 0.30 of 1/5 days at 1.1 / 3.5 / 7.4 /
     14 hPa; off with `topDragDays` 0). Nothing breaks without it: the
     step stays stable, the largest wind's horizontal Courant number in
     the top layer reaches 0.54 (July) / 0.62 (January) against 0.24 /
     0.22, the vertical one stays at 0.08 or less above 100 hPa. The top layer
     has no momentum sink: its winter jet grows from 66 to 131 m/s at 54S
     (still 1.2 m/s/day at the end) and to 122 m/s at 59N, the largest
     edge wind to 167 / 195 m/s, its eddy kinetic energy to 197 / 250
     m²/s² (27 / 43 with the sponge), twice the 3.5 hPa layer's, its
     eddy temperature to 4.4 / 4.3 K rms (2.1 / 2.5) and its rms
     divergence to 8.9 / 10.6·10⁻⁶ /s (4.7 / 6.0): the waves that reach
     the lid stay under it. The winter pole cools without bound: 70-90S
     181 K at 1.1 hPa and 187 K at 7-24 hPa on day 121, still falling
     0.7 and 1.0-1.2 K/day.
   - What it does where it should not. Its zonal force on day 121
     (m/s/day): −11.3 / −11.2 / −9.8 at 65 / 55 / 45S and +4.5 at 25-35N
     (on the summer easterlies) in the 1.1 hPa layer, −7.0 / −6.7 at
     3.5 hPa, −5.1 / −5.0 at 7.4 hPa, −2.4 / −2.5 at 14 hPa (65 / 55S);
     on day 304 −10.1 / −9.0 at 65 / 55N at 1.1 hPa, −5.6 / −4.7 at
     3.5 hPa. By downward control (f v* = −F, Haynes et al. 1991) it
     drives 80 kg/m/s poleward across 50S below 24 hPa (57 across 50N in
     January), and the descent over the cap at the 70-90° cap's static
     stability warms it 2.3 / 2.0 / 1.8 / 1.4 / 0.9 / 0.5 K/day at 1.1 /
     3.5 / 7.4 / 14 / 24 / 37 hPa; the ascent over the summer cap cools it
     1.4 / 1.3 / 0.9 / 0.6 K/day at 1.1-14 hPa. The polar-night jet it
     leaves is 43-47 m/s at 3.5-37 hPa near 60S in July and 33-35 m/s near
     64N in January, against Earth's 60-80 m/s near 60° at 1-10 hPa in
     winter (CIRA-86, Fleming et al. 1990; the SPARC climatologies,
     Randel et al. 2004; from memory, not re-read here), and the winter
     pole 205-227 K at 1-53 hPa where the vortex is 185-195 K at 30-50 hPa.
   - What published models do at a low lid (from memory of their
     documentation, not re-read here): CAM-SE adds ∇² viscosity on the
     wind and temperature in its top three layers (`nu_top`, Lauritzen et
     al. 2018, JAMES 10, 1537); CAM-FV raises its divergence damping
     there; ECHAM5 lowers the order of its hyperdiffusion in the
     uppermost layers (Roeckner et al. 2003, MPI Report 349); the IFS
     raises its diffusion in the top levels and replaced its Rayleigh
     friction above 10 hPa by Scinocca's (2003) non-orographic
     gravity-wave drag in 2009 (Orr et al. 2010, J. Climate 23, 5905);
     CMAM's sponge acts on the departures from the zonal mean alone,
     after Shepherd, Semeniuk & Koshyk (1996, JGR 101, 23447), who showed
     that a Rayleigh drag on the zonal-mean wind drives a spurious
     downward-control circulation below it. None of these five, as
     remembered, keeps a drag on the mean wind at its top.
   - Built, both engines (`js/dynamics/sponge.module.js`,
     `js/physics/gravityWaves.module.js`). The sponge: in the layers above
     `spongeSigma` 0.005 the departure of each edge's normal velocity
     from its zonal mean relaxes at up to 1/`spongeDays` (1 day), falling
     linearly in σ (0.78 and 0.29 /day at 1.1 and 3.5 hPa on bl34),
     implicitly and once per step with the closures; the zonal mean is
     the band mean of Section 3.7's cell winds over bands one mesh
     spacing wide, interpolated to the edges as angular velocity, so a
     zonal flow is 0.42 / 0.11 m/s rms (1.7 / 0.4 %) from it at N=16 / 32
     and the sponge's torque is 1.4·10⁻³ / 2.9·10⁻⁴ of a Rayleigh drag's
     at the same rate (`test/sponge.test.mjs`); its zonal force on the
     day-121 state is under 0.08 m/s/day. The gravity-wave drag:
     Alexander & Dunkerton (1999, JAS 56, 4167), from the layer nearest
     315 hPa (316 hPa), in two directions (east, north), phase speeds
     u₀ ± j·4 m/s to ±100 m/s about the source wind, a Gaussian of
     half-width 40 m/s antisymmetric about u₀ carrying 4.3 mPa of
     absolute flux in each direction, one wavelength of 300 km (values of
     the order used with this scheme in GFDL's AM3 and in MiMA, Jucker &
     Gerber 2017, J. Climate 30, 7339; from memory, unchecked); each wave
     leaves its whole flux at its critical level or where it first
     exceeds the saturation flux ρk|c − u|³/(2N) (its grid-box mean flux
     in the runs of this section; its amplitude where present since the
     review below), and what reaches the
     top layer is left there (Shepherd & Shaw 2004, JAS 61, 2899), so
     each column launches and keeps zero net momentum (to 1·10⁻¹⁸ Pa
     against 4·10⁻³ Pa of deposit, `test/gravityWaves.test.mjs`); computed per cell after the
     physics, applied to the edges with the momentum mixing; the kinetic
     energy both change returns as heat (total energy to 1e-9 of the
     change). CPU and GPU: the accelerations after one step at N=8 differ
     by 4.4·10⁻⁷ m/s/day rms of up to 7.8 m/s/day, the sponge's wind
     after 20 steps by 1.6·10⁻³ m/s of the 4.2 m/s it moves. The Rayleigh drag
     stays as an option (`topDragDays`, default 0).
   - What the drag does. Its absolute flux rising through 53 hPa
     (~20 km) is 6.4-8.1 mPa in July and 6.8-8.1 in January at every
     latitude (Earth: 1-5 mPa at 20-25 km outside the southern winter's
     5-10, Ern et al. 2004, JGR 109, D20103, Geller et al. 2013, J.
     Climate 26, 6383; from memory), 5.2-7.3 mPa through 3.5 hPa, nearly
     all of it left in the top layer. Its force on day 121: −6.6 / −6.9 /
     −5.1 m/s/day at 65 / 55 / 45S and +3.7 to +5.0 at 25-55N in the
     1.1 hPa layer; +1.0 to +2.0 at 25-45S at 3.5 hPa where the eastward
     waves meet their critical levels in the jet, under 0.6 in magnitude
     below; on day 304 −3.9 to −4.7 at 55-85N and +3.1 to +4.9 at 15-55S.
     Its descent warms the 70-90S cap 1.3 / 0.8 / 0.4 / 0.2 / 0.1 K/day
     at 1.1 / 3.5 / 7.4 / 14 / 24 hPa (13 kg/m/s across 50S at 3.5 hPa,
     7-10 below 24 hPa: an eighth of the Rayleigh drag's there).
   - Before → after, the strongest zonal-mean westerly of the winter
     hemisphere (m/s at latitude) and the winter cap's layer-mean
     temperature (K) on the last day: Rayleigh (as it was) / none /
     sponge and drag (the defaults) / drag alone / sponge and drag at
     8.6 mPa (`GRAVITY_WAVES='{"flux":8.6e-3}'`), July 70-90S and
     January 70-90N:

     | layer | July jet | July 70-90S | January jet | January 70-90N |
     |---|---|---|---|---|
     | 1.1 hPa | 61@61S / 131@54S / 112@39S / 112@44S / 110@46S | 215 / 181 / 204 / 199 / 214 | 54@66N / 122@59N / 92@41N / – / 86@31N | 217 / 185 / 212 / – / 214 |
     | 3.5 hPa | 44 / 92 / 67 / 81 / 63 | 227 / 190 / 212 / 205 / 222 | 35 / 76 / 62 / – / 62 | 225 / 192 / 216 / – / 218 |
     | 7.4 hPa | 43 / 77 / 65 / 74 / 57 | 223 / 187 / 203 / 196 / 213 | 33 / 64 / 51 / – / 47 | 221 / 192 / 208 / – / 211 |
     | 14 hPa | 44 / 66 / 63 / 68 / 58 | 219 / 187 / 197 / 192 / 203 | 33 / 57 / 51 / – / 50 | 218 / 193 / 202 / – / 204 |
     | 24 hPa | 46 / 58 / 59 / 62 / 57 | 213 / 187 / 193 / 190 / 197 | 34 / 51 / 49 / – / 51 | 213 / 194 / 199 / – / 200 |
     | 37 hPa | 47 / 52 / 55 / 58 / 54 | 208 / 188 / 192 / 190 / 195 | 34 / 46 / 47 / – / 49 | 209 / 196 / 199 / – / 199 |
     | 53 hPa | 48 / 50 / 53 / 56 / 51 | 206 / 190 / 193 / 192 / 196 | 33 / 43 / 44 / – / 47 | 207 / 197 / 199 / – / 199 |

     The jets at 7.4-53 hPa sit at 49-56S and 61-64N in every run. The
     zonal-mean wind in the 55° and 65° bands (m/s, 1.1 / 3.5 / 7.4 / 14
     / 24 hPa): July 55S 59 / 41 / 39 / 42 / 44 → 84 / 67 / 63 / 60 / 55
     (none: 129 / 90 / 75 / 64 / 56), 65S 60 / 43 / 41 / 39 / 37 → 69 /
     53 / 49 / 43 / 36; January 55N 47 / 28 / 25 / 26 / 29 → 74 / 50 /
     47 / 44 / 42 (none: 120 / 74 / 55 / 45 / 40), 65N 53 / 34 / 32 / 33
     / 34 → 57 / 47 / 49 / 50 / 47. The summer cap at 1.1 / 3.5 / 7.4
     hPa: July 70-90N 268 / 255 / 240 → 269 / 259 / 247 (none 275 / 265 /
     250), January 70-90S 269 / 262 / 248 → 271 / 265 / 253 (274 / 270 /
     255). Trends over the last ten days at 70-90S in July (K/day, 1.1 /
     3.5 / 7.4 / 14 / 24 / 37 / 53 hPa): Rayleigh +0.10 / −0.21 / −0.20 /
     −0.22 / −0.17 / −0.14 / −0.12, defaults −0.64 / −0.71 / −0.76 /
     −0.70 / −0.58 / −0.45 / −0.35, none −0.73 / −0.94 / −1.15 / −1.12 /
     −0.99 / −0.85 / −0.71; at 70-90N in January −0.33 / −0.25 / −0.34 /
     −0.33 / −0.41 / −0.46 / −0.49 → −0.38 / −0.36 / −0.43 / −0.51 /
     −0.52 / −0.48 / −0.43. The global layer means move by 1.3 K or
     less; the tropics (20S-20N) warm 1.7 / 3.0 / 3.0 / 2.2 / 1.4 K at
     1.1 / 3.5 / 7.4 / 14 / 24 hPa in July and 2.9 / 3.3 / 2.1 / 1.4 /
     0.8 K in January (to 256 / 248 / 234 / 223 / 215 K, AFGL's tropical
     layer means 261 / 252 / 240 / 230 / 222), as the Rayleigh drag's
     ascent over them goes; the tropical cold point (85-101 hPa) is
     197-198 K and falls 0.2-0.3 K/day in every run. Without the
     sponge (drag alone) the top two layers carry 4 times the eddy
     kinetic energy (98 / 113 against 24 / 27 m²/s² at 1.1 / 3.5 hPa on
     day 121) and the largest winds 139 / 120 against 124 / 85 m/s, and
     the 70-90S cap is 5-7 K colder at 1.1-7.4 hPa, inside what a
     round-off twin moves it (review, below): the sponge keeps the top
     quiet. The troposphere does
     not move: means of days 1-10 and 21-30, Rayleigh → defaults, July
     balance −6.16 → −6.01 and −12.24 → −12.25, OLR 234.23 → 234.21 and
     232.36 → 232.48 W/m², rain 2.163 → 2.164 and 2.697 → 2.710 mm/d;
     January −4.91 → −4.71 and −8.88 → −9.01, 225.21 → 225.18 and
     224.22 → 224.24, 2.121 → 2.117 and 2.650 → 2.691 (the runs without a
     sponge or at 8.6 mPa differ from the defaults by up to 2.0 W/m² and
     0.15 mm/d over days 21-30; a round-off twin of the defaults moves
     them by up to 1.4 W/m² and 0.07 mm/d, review below). Three days
     with the defaults gave the explicit run's lines bit for bit.
   - Whether a gravity-wave drag is needed at this lid: yes. Without any
     momentum sink the top layer's jet and its cooling do not stop
     (above); the drag closes the vortex where Earth's is closed, by
     waves breaking above the lid, and gives what reaches the lid to the
     top layer, as momentum conservation asks. Its cost: 0.8 ms of the
     23 ms N=64 step and 2.4 ms of the 98 ms N=128 step, the sponge's
     included (below). Doubling its flux warms the 1.1-7.4 hPa winter
     cap by 10 K (July) and 2-3 K (January) and weakens the July jet at
     3.5-24 hPa by 2-8 m/s; not adopted: the flux it would send through
     20 km, 13-16 mPa, is beyond what is observed.
   - bl36 (optional; `sigmaInterfaces('bl36')`, `LEVELS=bl36`): bl34
     with its 0-2.19 hPa layer split at 0.3 and 1 hPa (36 layers,
     midpoints 0.15 / 0.64 / 1.6 hPa); `scripts/remapState.mjs` carries a
     saved state onto it (each new layer takes the old top layer's θ and
     wind; back onto bl34 to 1e-12, `test/levels.test.mjs`). The
     radiation benchmark on it (`LEVELS=bl36 node
     scripts/radiationBenchmark.mjs`, the g-points as fitted on bl34),
     cooling K/day of the 0.15 / 0.64 / 1.6 / 3.6 hPa layers against
     RRTMG: TROP −3.94 / −11.09 / −9.90 / −6.03 (−5.40 / −9.92 / −9.42 /
     −6.00), MLS −3.53 / −11.98 / −11.23 / −6.53 (−6.14 / −12.19 /
     −10.16 / −6.42), MLW −6.75 / −12.12 / −8.17 / −3.87 (−7.69 / −10.67
     / −7.90 / −3.86), SAW −9.08 / −9.88 / −5.32 / −2.76 (−9.89 / −9.02 /
     −5.37 / −2.90; the MLS and SAW references stop at 6.7 and 10 Pa, 78
     and 67 % of the top layer's mass, and count the rest as not
     cooling); the cooling rms at 3-200 hPa 0.08 / 0.10 / 0.09 / 0.12 →
     0.07 / 0.09 / 0.09 / 0.11 K/day, OLR and surface fluxes within 0.03
     W/m² of bl34's, the shortwave heating rms at 1-200 hPa 0.37 / 0.37 /
     0.32 / 0.29 / 0.18 → 0.31 / 0.33 / 0.29 / 0.28 / 0.13 K/day; doubled
     CO₂ at the top 2.63 → 2.65 W/m² (LBLRTM 2.84; its 630-700 cm⁻¹ band
     −0.74 → −0.72 against −0.57, `LEVELS=bl36 node
     scripts/longwaveFit.mjs`), the stratosphere-adjusted forcing 5.49 →
     5.49 W/m² with the layers above 10 hPa cooling 19.3 / 14.7 / 11.4 K
     (1.1 / 3.5 / 7.4 hPa) → 16.8 / 20.6 / 19.0 / 14.6 / 11.3 K (0.15 /
     0.64 / 1.6 / 3.5 / 7.4 hPa). The g-points' cooling rms above 200 hPa
     by band (MLS) rises on bl36, 630-700 cm⁻¹ 0.32 → 0.59, 700-820 0.17
     → 0.25, 10-350 0.26 → 0.32 K/day: the reduction was binned on bl34's
     columns, and a bl36 run would want its own (`LEVELS=bl36 WRITE=1`).
   - Final defaults: `topDragDays` 0 (was 5; `topSigma` 0.02 unused),
     `spongeSigma` 0.005 and `spongeDays` 1 (`SPONGE`, the sponge on the
     eddies; `surface` options), `gravityWaves` on with `flux` 4.3·10⁻³
     Pa, `sourcePressure` 31500 Pa, `halfWidth` 40 m/s, `maxSpeed` 100
     m/s, `speedStep` 4 m/s, `wavelength` 300 km, `minimumFrequency`
     0.005 /s and, since the review below, `breakingAmplitude` 0.4 m²/s²
     (`GRAVITY_WAVES`; `gravityWaves: false` for none); the
     levels stay bl34. What still misses: the winter cap's top layer,
     204 K in July and 212 K in January against AFGL's 247 K subarctic
     winter layer mean and still cooling 0.4-0.6 K/day, and the 1.1 hPa
     jet growing 1.9 m/s/day at 39S in July and 0.4 at 41N in January
     over the last ten days (one 0-2.2 hPa layer
     holds the whole mesosphere's drag and the polar stratopause it
     drives); the July jet at 3.5-37 hPa 55-67 m/s against 60-80; the
     January vortex at 14-53 hPa 199-202 K; the waves' flux at 20 km
     4.6-5.7 mPa at every latitude with `breakingAmplitude` (6.2-8.1
     without it), at the top of the observed 1-5 outside the southern
     winter, from a source uniform in latitude and season; the tropical
     upper stratosphere 5-7 K below AFGL's means; and
     the 30-day runs start from states whose stratosphere the earlier
     radiation warmed, so the caps are still adjusting.
     Thirty N=64 days on bl36 from nine64_day0091 remapped, the defaults
     (the sponge then covers the 0.15 / 0.64 / 1.6 / 3.5 hPa layers at
     0.97 / 0.87 / 0.68 / 0.29 /day, and the waves' flux that reaches the
     lid is left in the 0-0.3 hPa layer): the winter cap (70-90S, day 121)
     247 / 257 / 236 / 219 / 207 / 198 / 193 / 192 / 193 K at 0.15 / 0.64
     / 1.6 / 3.5 / 7.4 / 14 / 24 / 37 / 53 hPa, a polar winter
     stratopause, 244 K in the mass mean of 0-2.2 hPa where bl34's one
     layer is 204 K (AFGL subarctic winter 247); the global means 234 /
     253 / 254 / 243 K at 0.15 / 0.64 / 1.6 / 3.5 hPa, the tropical
     stratopause 261 K at 1.6 hPa, the summer cap 216 / 253 / 266 / 258 K;
     the winter jet 62 / 63 / 61 / 56 / 54 m/s at 3.5 / 7.4 / 14 / 24 /
     37 hPa near 51-54S (bl34 67 / 65 / 63 / 59 / 55). But the 0-0.3 hPa
     layer runs away at the equator: a zonal-mean westerly of 221 m/s at
     5S (202 at 5N) by day 121, the largest edge wind 349 m/s on day 113
     (horizontal Courant number 1.12; 258 m/s, 0.83 on day 121). The
     sponge drives it: without it (`spongeDays` 0, otherwise the same)
     the layer's equatorial westerly is 47 m/s at 5S on day 121, the
     sponge having taken the eastward equatorial waves' momentum into the
     zonal mean of a 30 Pa layer, but its eddy kinetic energy reaches
     1080 m²/s² and its largest wind 275 m/s (Courant 0.88). On bl34 the
     1.1 hPa layer's equatorial wind is −19 / −34 m/s at 5S / 5N with the
     sponge and −24 / −39 without. bl36's 0-0.3 hPa layer wants a
     treatment of its own before a season on it. Days
     21-30 of the troposphere: balance −11.03, OLR 233.37 W/m², rain
     2.604 mm/d (bl34 −12.25, 232.48, 2.710, within the realizations'
     spread). Three N=128 days from eight128_day0183 remapped: largest
     wind 94-103 m/s, Courant numbers 0.33 / 0.09 at most, day 186 ASR
     241.8, OLR 239.9 W/m², rain 1.85 mm/d (bl34 241.8, 239.8, 1.85).
   - Cost, `js/gpu/profile.module.js` over 128 steps under the exclusive
     lock, two runs each, the step median: N=64 (nine64_day0091) 22.12 /
     22.17 ms with the Rayleigh drag → 23.22 / 22.91 ms (+3.6 %; the
     physics pass, where the drag's column kernel runs, 4.31 → 4.78 ms,
     the sponge 0.08 ms), on bl36 25.36 / 25.54 ms (+10 % over bl34);
     N=128 (eight128_day0183) 95.91 / 95.88 → 98.20 / 98.37 ms (+2.5 %;
     physics pass 18.40 → 19.84, sponge 0.27), 50.3 s a model day where
     it was 49.1, on bl36 107.49 / 107.33 ms (+9 %), 55.0 s a model day.
     The tables' runs used a search up the column for each wave; the
     committed form leaves each side's waves as contiguous runs of phase
     speed and agrees with it to 5·10⁻¹⁴ m/s/day on their states.
   - Review (Oct 2), N=64 and N=128 GPU runs one at a time, the CPU
     engine's operators on the saved states.
     - Round-off twins: the defaults above from nine64_day0091 and
       nine64_day0274 with the lowest layer's θ perturbed by ±1·10⁻⁴ K.
       July day 121, 70-90S at 1.1 / 3.5 / 7.4 / 14 / 24 / 37 / 53 hPa
       208 / 217 / 209 / 201 / 196 / 195 / 196 K against 204 / 212 / 203
       / 197 / 193 / 192 / 193, the winter jet 107 / 62 / 60 / 58 / 55 /
       51 / 50 against 112 / 67 / 65 / 63 / 59 / 55 / 53 m/s; January day
       304, 70-90N 212 / 213 / 208 / 205 / 203 / 204 / 204 against 212 /
       216 / 208 / 202 / 199 / 199 / 199 K, the jet at 7.4-53 hPa 42 /
       33 / 33 / 34 / 34 against 51 / 51 / 49 / 47 / 44 m/s; days 21-30
       balance / OLR / rain July −11.67 / 233.41 / 2.640 against −12.25 /
       232.48 / 2.710, January −7.63 / 224.17 / 2.642 against −9.01 /
       224.24 / 2.691. One 30-day realization moves the winter cap by up
       to 6 K and the January vortex by up to 18 m/s. What stands outside
       that: the Rayleigh drag's warming of the July cap at 1.1-53 hPa,
       10-22 K over the defaults and 16-37 K over no sink at all. What
       does not: the sponge's 5-7 K, the doubled flux's 0-10 K, and every January
       difference of the jet at 7.4-53 hPa in the table above. The
       troposphere's day means under the Rayleigh drag and the defaults
       agree within the twins.
     - The sponge on real states (nine64_day0091, nine64_day0274,
       eight128_day0183; one step, zonal means in 10° bands by Section
       3.7's reconstruction): its zonal force is 0.053 m/s/day at most at
       1.1 hPa and 0.021 at 3.5 hPa, where a Rayleigh drag at the same
       rates gives up to 50 and 14; its change of the layer's axial
       angular momentum is 0.3-3 % of that drag's, of either sign.
     - CPU and GPU, the whole model from nine64_day0091 (ocean off), the
       top six layers: the wind 1.0 / 2.7 / 6.4·10⁻⁴ m/s rms apart after
       1 / 4 / 16 steps (1.2·10⁻² at most) of the 0.09 / 0.33 / 1.3 m/s
       it moves, θ 1.5·10⁻³ K rms after 16 steps of 3.9 K; with the
       Rayleigh drag in their place 7.1·10⁻⁴ (4.6·10⁻² at most).
     - The gravity waves' breaking. With each wave's grid-box mean flux
       tested against the saturation flux, saturation takes 0.3 mPa of
       the 8.4 that rise from the source before the lid on the nine64
       states (6.5 against 6.8 mPa through 3.5 hPa with breaking left
       out): the scheme acts as a critical-level filter. Alexander &
       Dunkerton (1999) test each wave's amplitude where it is present
       and take the grid-box mean flux as its intermittent fraction ε;
       `breakingAmplitude` B_w tests ρ₀ B_w times the spectrum's shape,
       ρ₀ the source layer's density, with B_w 0.4 m²/s² (GFDL cg_drag's
       wide-spectrum amplitude, from memory, not re-read). The flux
       rising through 53 hPa is then 4.6-5.7 mPa against 6.2-8.1, through
       3.5 hPa 2.4-3.3 against 5.6-7.3, on the day-121 and day-304
       states. Thirty days with it from the same states: July 70-90S 206
       / 213 / 205 / 198 / 194 / 193 / 194 K at 1.1-53 hPa, jet 110@41S /
       65 / 61 / 58 / 54 / 52 / 50 m/s; January 70-90N 210 / 215 / 209 /
       203 / 199 / 199 / 199 K, jet 88@44N / 61 / 45 / 44 / 42 / 38 / 35
       m/s; days 21-30 July −12.29 / 233.11 / 2.654, January −8.11 /
       223.99 / 2.705: within the twins. It is the default
       (`GRAVITY_WAVES`; `breakingAmplitude: null` tests the grid-box
       mean as before); the GPU's accelerations match the CPU's to
       8.7·10⁻⁸ m/s/day rms at N=8 (`test/gravityWaves.test.mjs`).
     - Stability with the final defaults: the July and January 30-day
       runs reach a horizontal Courant number of 0.40 (124 m/s) in the
       top six layers, the vertical 0.09. The strongest polar-night jet:
       30 days from nine64_day0091 with no sink at all (the top layer's
       jet 131 m/s at 54S, Courant 0.56; the earlier run bit for bit),
       then 10 days with the defaults: Courant 0.54 → 0.45, the top
       layer's jet 129 → 127 m/s, the jet at 3.5-37 hPa still
       strengthening to 98 / 86 / 77 / 68 / 61 m/s and the cap at 7.4-37
       hPa 184-186 K on day 131: the drag does not undo a vortex that
       strong in ten days. Three N=128 days from eight128_day0183: largest
       wind 85.5 m/s, Courant 0.23 at most in the top six layers, day 186
       ASR 241.8, OLR 239.8 W/m², rain 1.85 mm/d, albedo 0.290.
     - bl36: nine64_day0091 carried onto it and back keeps the mass, the
       enthalpy (the layer-integral Exner function of each new layer
       sums to the old one's) and the kinetic energy to 1·10⁻¹⁴, θ to
       2·10⁻¹³ K; the radiation benchmark on it gives the numbers above.
       Ten days on it with the final defaults: largest wind 129 → 126 m/s
       (Courant 0.41 → 0.40), the 0-0.3 hPa layer's zonal mean −56 / −48
       m/s at 5S / 5N on day 101 and 1.0-1.9 mPa reaching it; the earlier
       run's runaway began after day 103, so whether the intermittent
       breaking prevents it is open. bl34 stays the default everywhere.
     - Cost under the exclusive lock, 128 steps, two runs each, the step
       median: N=64 21.94 / 21.90 ms with the Rayleigh drag, 22.63 /
       22.63 with the grid-box breaking, 22.71 / 22.70 with the defaults
       (+3.6 %); N=128 95.82 / 95.82, 98.07 / 98.07, 98.21 / 98.07 ms
       (+2.4 %), 50.2 s a model day.
   The gravity waves' sources and bl36's top (Oct 2, second round).
   N=64 GPU runs one at a time from copies of nine64_day0091 (to day
   121) and nine64_day0274 (to day 304), `STRATOSPHERE=1`, the upper
   winds line now with the 5S-5N zonal-mean wind; budgets from the CPU
   engine's operators on the saved states.
   - The parameters against their sources. Read here: MiMA's
     `cg_drag.f90` and `input/input.nml` (github mjucker/MiMA, master);
     GFDL's `atmos_param/cg_drag/cg_drag.F90` (NOAA-GFDL/atmos_phys,
     main) and AM4's `run/input.nml` (NOAA-GFDL/AM4); Garfinkel et al.
     (2022, JAMES 14, e2021MS002568) Appendix A, MiMA's control; Hertzog
     et al. (2008, JAS 65, 3056), Vorcore superpressure balloons;
     Corcos et al. (2021, JGR 126, e2021JD035165, abstract), Strateole-2
     balloons; Geller et al. (2013, J. Climate 26, 6383, abstract) and its
     Table 1 as Holt et al. (2017, QJRMS 143, 2481) quote it. Alexander &
     Dunkerton (1999) itself was not reachable (the AMS server refused);
     its values below are the two codes' comments on it.

     | parameter | GFDL cg_drag / AM4 | MiMA code / input.nml / Garfinkel 2022 | observed | before | set |
     |---|---|---|---|---|---|
     | spectrum | Gaussian, peak at ground-relative c = 0 (`flag` 1, not in AM4's namelist) | peak at c − u₀ = 0 (`flag` 0 in input.nml; Garfinkel: "symmetric about the zonal wind at the source level") | meridional flux symmetric about intrinsic c = 0 (Hertzog 2008) | about u₀ | about u₀ |
     | half-width c_w | 40 m/s (earlier 50, 25) | 35 (input.nml, Garfinkel) | intrinsic spread broader than AD99's broad spectrum, exponential scale ≤ 92 m/s, an upper bound (Hertzog 2008) | 40 | 40 |
     | B_w | 0.4 m²/s² | 0.4 | – | 0.4 | 0.4 |
     | flux | `Bt_0` 0.005, `Bt_nh` 0.002, `Bt_sh` −0.00025 m²/s² × 1.5 ρ₀, tanh at ±30° over 5° (AM4) | code 4 mPa + 1 / −1 at ±30°; input.nml 4.3 mPa everywhere (`Bt_eq` = `Bt_0`, `Bt_nh` = `Bt_sh` = 0); Garfinkel 4.3 + 3.5 mPa beyond ±15° (tanh over 10°), "to keep the polar vortex from becoming too strong" | a non-zero background away from convection, no value in the abstract (Corcos 2021); 2.5 raw, 3.2 corrected, ≤ 6.4 with the unresolved high frequencies over 50-75S in spring (Hertzog 2008); HIRDLS global means 1.8-4.1 mPa at 20 km (Geller 2013 Table 1), perhaps 2-4 times low | 4.3 mPa per direction, uniform | 4.3, uniform (eq. A3 built: `equatorialFlux`, `northFlux`, `southFlux`, `edge`, `width`) |
     | launch | 315 hPa at the equator, level index (K+1) − (K+1−k₀) cos φ | the same code; "315 hPa in the tropics" after Donner et al. (2011) | – | 316 hPa everywhere | nearest σ = 0.315^cos φ (the index rule read as log-pressure height), never the lowest layer (`sourceDescent`) |
     | phase speeds | ±99.6 m/s by 1.2 (code), 2.4 (AM4) | ±99.6 by 1.2 (code), 2 (Garfinkel) | – | u₀ ± 4 j to 100 | u₀ ± 2 j to 100 |
     | wavelength | 300 km (`nk` 1) | 300 km | – | 300 km | 300 km |
     | breaking | B₀/(c−u)³ ≥ ½ (ρ/ρ₀) k/N on the intermittent amplitude, ε = F/(ρ₀ ΣB₀) (GFDL with a 1.5 "unexplained", in m²/s²) | the same, the 1.5 dropped and F in Pa (Cohen et al. 2013, appendix) | – | the same | the same |
     | at the lid | `dump_flux`: what remains into the top level | "deposited evenly in the levels above 0.85 hPa" (`damp_level_pressure` 85 Pa; Cohen 2013: 0.5 hPa); no sponge | – | the top layer | spread over the layers whose midpoints lie above 85 Pa, untested (below) |

     The flux rising through 53 hPa, both directions, mPa, at 85S … 85N
     on the bl34 states nine64_day0091 / nine64_day0274: before 5.0-5.7
     (July); the defaults 4.5-5.4 / 4.5-5.4; with Garfinkel's 3.5 mPa
     increment 8.1-9.5 poleward of 30° (4.9-5.1 at the equator); with
     c_w 35 as well 7.4-8.9 there; with the increment and no descent
     9.0-10.1.
     The increment sends 2.5 times the balloons' 3.2 mPa through 20 km
     and beyond their 6.4 upper bound: not adopted. On the bl36 day-30
     states with the final defaults: 4.4-5.3 (July) and 4.4-5.3 (January)
     through 53 hPa, 1.9-3.0 / 1.6-2.9 through 3.5 hPa, 1.2-2.6 / 1.1-2.3
     into the lid layers at 1 hPa, 0.35-0.79 / 0.33-0.69 through 0.3 hPa.
     Both engines; the accelerations after one N=8 step agree to 2.6e-7
     m/s/day rms on bl36, 1.1e-7 on cam26 (`test/gravityWaves.test.mjs`).
   - What drives bl36's top. Thirty July days with the parameters
     above, waves tested in every layer and the lid's flux spread over
     the 0-1 hPa layers (MiMA's rule), the sponge from σ 0.005: the 0-0.3
     hPa layer's 5S-5N wind −45 → −115 m/s, the 0.3-1 hPa layer's
     −4 → +73, the 1-2.2 hPa layer's −8 → +41 (the westerly descending
     from day 96), largest wind 188 m/s. On the day-121 state the waves
     push the top layer's equatorial band at −13.9 m/s/day (−4.1 with that
     layer untested), against +16.9 from the resolved dynamics: the
     westward waves meet their critical levels in its easterly, which
     drives it further east-to-west, a wave-mean-flow feedback in a 30 Pa
     layer whose mean wind stands for the whole mesosphere, where the
     real winds turn. Below it, the 0.3-1 hPa westerly is the sponge's:
     with every lid layer untested it persists (+72), and it goes when the
     sponge leaves the layer. July day 121, the variants (5S-5N wind at
     0.15 / 0.64 / 1.6 hPa, m/s; 70-90S at the same layers, K; the
     strongest westerly at 0.15 hPa; the run's largest wind):

     | top treatment | equator | winter cap | jet | largest |
     |---|---|---|---|---|
     | all tested, lid 0-1 hPa, sponge from σ 0.005 (bl34's table) | −115 / +73 / +41 | 229 / 238 / 223 | 113@74S | 188 |
     | top layer untested, lid 0-1 hPa (bl34's table) | −56 / +76 / −2 | 206 / 222 / 217 | 111@59S | 143 |
     | top layer untested, lid the top layer | −39 / +100 / +1 | 222 / 234 / 231 | 96@56S | 155 |
     | lid layers untested, lid 0-1 hPa | −4 / +72 / −33 | 203 / 220 / 218 | 117@59S | 145 |
     | lid layers untested, lid 0-2.2 hPa | −69 / −16 / +1 | 192 / 209 / 213 | 143@54S | 174 |
     | lid layers untested, lid 0-1 hPa, sponge from 0.77 hPa | −22 / −4 / −21 | 209 / 225 / 218 | 118@61S | 157 |
     | the same without a sponge | −26 / −8 / −11 | 199 / 216 / 211 | 135@61S | 251 |
     | top layer untested, lid the top layer, sponge from 0.77 hPa | −52 / +58 / +9 | 232 / 240 / 224 | 92@46S | 163 |
     | the defaults (sponge from 78 Pa) / its round-off twin | −15 / +2 / −20; −26 / −4 / −13 | 213 / 229 / 221; 206 / 222 / 216 | 109@64S; 125@59S | 194; 170 |

     The earlier bl36 run's +221 m/s westerly in the top layer came with
     the grid-box breaking test, which lets about twice the flux reach
     the lid. At solstice the equatorial stratopause is in the SAO's
     easterly phase (Fleming et al. 1990, CIRA-86, section B: "easterlies
     during the solstices near the stratopause maximum"); every variant
     with a lid layer tested, or with the sponge reaching 0.64 hPa, builds
     a 58-100 m/s westerly there.
   - Built (both engines). The waves are tested for breaking and
     critical levels only below the lid layers; what rises into those
     (midpoints above `lidPressure` 85 Pa: 0.15 and 0.64 hPa on bl36,
     the top layer on bl34 and cam26, where nothing changes) is spread
     over them at one acceleration, conserving each column's momentum
     (`test/gravityWaves.test.mjs`). The sponge on the eddies begins on
     bl36 at 78 Pa, the IFS's onset (IFS CY48r1 Part III, 2.2.11: from
     0.78 hPa up, the zonal mean undamped after Shepherd et al. 1996),
     falling linearly in σ from 1/day at the top (0.81 / 0.17 /day at
     0.15 / 0.64 hPa, nothing below); bl34 keeps σ 0.005 (0.78 / 0.29 at
     1.1 / 3.5 hPa; `spongeSigmaFor`). No mean-flow drag: the variants
     above settle without one, and the published mean-flow drags at a low
     lid (the IFS's Rayleigh friction above the stratopause before
     Cy35r3, Orr et al. 2010) were dropped for gravity-wave drag.
   - The longwave g-points on bl36 (`longwaveTableFor`, table
     `js/physics/longwaveTableBl36.module.js`, `LEVELS=bl36
     UPPER_WEIGHT=100 FIT=4000 WRITE=1 node scripts/longwaveFit.mjs`;
     the fit now scores every layer above 3 hPa that RRTMG's column
     covers to 90 % of its mass, so bl34's score is unchanged). Cooling
     K/day at 0.15 / 0.64 / 1.6 hPa against RRTMG, bl34's table → bl36's:
     TROP −3.92 / −11.09 / −9.90 → −5.05 / −11.08 / −8.83 (−5.40 / −9.92 /
     −9.42); MLS −3.53 / −11.98 / −11.23 → −4.43 / −11.91 / −9.95 (−6.14
     over the 78 % of the top layer RRTMG covers / −12.19 / −10.16); MLW
     −6.75 / −12.11 / −8.16 → −8.56 / −12.01 / −7.40 (−7.69 / −10.67 /
     −7.90); SAW −9.08 / −9.88 / −5.32 → −12.05 / −9.81 / −4.72 (−9.89
     over 67 % / −9.02 / −5.37). Where RRTMG covers the layer the misses
     are −6 / +12 / −6 % (TROP), −2 / −2 % (MLS), +11 / +13 / −6 % (MLW),
     +9 / −12 % (SAW), against −27 / +12 / +5, −2 / +11, −12 / +14 / +3,
     +10 / −1 % with bl34's table; UPPER_WEIGHT 30 and 300 and a
     continued fits gave 24 / 13.5 / 25 / 17 % at worst. The 3-30 hPa layers
     cool within 0.18 K/day of RRTMG (rms 0.05 / 0.06 / 0.07 / 0.12 K/day,
     bl34's table 0.11 / 0.11 / 0.11 / 0.15); OLR +0.3 / +0.5 / −1.0 /
     −1.3 W/m², surface downward −1.1 / −0.7 / +1.7 / +2.5, at 200 hPa
     +0.8 / +1.1 / −0.1 / −0.1; LBLRTM's doubled CO₂ 2.66 / 5.71 / 1.66
     W/m² (2.84 / 5.54 / 1.68; bl34's table on bl34 2.63 / 5.72 / 1.65),
     its stratosphere-adjusted forcing 5.48 with the layers above 10 hPa
     cooling 16.4 / 18.7 / 18.7 / 16.1 / 11.9 K; vapour × 1.2 4.04 / 4.76
     / 12.18 (3.79 / 4.52 / 11.55); methane and nitrous oxide from none
     3.13 / 2.94 / 1.15 (3.60 / 3.45 / 1.08); ICRCCM tropical surface
     downward −5.2 W/m². The 0.64 hPa layer still overcools by 9-13 % in
     three of four columns: the shortfall of the reduction, not refitted
     further. Both engines carry the table alike on bl36
     (`test/gasRadiation.test.mjs`, the layers above 30 hPa 1e-6 to 2e-5
     in relative rms).
   - Fresh starts (no FROM; the atlas ocean). The initial stratosphere is
     the radiative-convective column of `equilibriumProfile` (the model's
     radiation, on the polar cell 0 at the global-mean sunlight, AFGL
     ozone) laid on every column: on bl36 226 / 247 / 260 / 265 / 256 K at
     0.15 / 0.64 / 1.6 / 3.6 / 7.5 hPa, where AFGL's tropical and
     midlatitude-summer layer means are 238 / 265 / 264 / 252 / 240 and
     239 / 270 / 271 / 257 / 243 and subarctic winter's 253 / 256 / 241 /
     227 / 218: cold by 10-23 K above 1 hPa, warm by 8-16 K at 3.6-7.5 hPa,
     the stratopause at 3.6 hPa instead of near 1 hPa, and no pole-to-pole
     contrast. It was read at each layer's σ, so over high ground the
     stratosphere was sea level's lifted by up to a scale height; on bl36
     at N=128 the top layers then reached 292 m/s in eight steps and NaN
     within the first day (bl34: 205 m/s after eight steps, 215 after 64,
     finite). Above σ 0.1 (blended
     over σ 0.3-0.1) it is now read at the layer's pressure,
     `js/physics/init.module.js`: 41 m/s after eight N=128 steps;
     at N=16 on terrain 65 m/s after 32 steps from rest, 175 before
     (`test/init.test.mjs`). Ten N=64 days on bl36: layer means global
     226 / 248 / 260 / 261 / 252 K at 0.15 / 0.64 / 1.6 / 3.5 / 7.4 hPa on
     day 1 and 229 / 254 / 260 / 248 / 235 on day 10 (the July runs' 230 /
     253 / 257 / 244 / 229), largest wind 95 m/s (Courant 0.28); the
     troposphere goes as on bl34 (albedo 0.29 → 0.43 on day 2 and 0.41 on
     day 10, ASR 193-203, OLR 208-213 W/m², rain 3.3-5.0 mm/d; bl34 the
     same within 0.004 in albedo and 1.2 W/m²): the fresh cloud is the
     troposphere's, not the levels'. Two N=128 days: largest wind 79 /
     94 m/s, Courant 0.25 / 0.30, day 2 ASR 201.4, OLR 213.6, albedo 0.408.
   - Defaults for new runs. `scripts/spinup.mjs` starts a fresh run on
     bl36 (`LEVELS` was cam26 by default) and the page's fresh start
     (`from=none`) on bl36; a run continues on its snapshot's grid, the
     page loads whatever levels a saved state carries, bl34 stays
     selectable (`LEVELS=bl34`) and its states load and remap as before
     (`test/levels.test.mjs`). A paired spin-up runs on bl36 when it
     starts fresh: `scripts/pairedSpinup.sh` calls spinup.mjs with the
     environment it is given, so with no `<PREFIX><N>_day*.bin` in OUT
     every N starts from the atlas on bl36; passing `LEVELS=bl36`
     explicitly makes a leftover bl34 snapshot under the prefix an error
     ("… is on bl34, not bl36") instead of a silent bl34 continuation.
   - Thirty days, the defaults, bl36 (and its round-off twin, the lowest
     layer's θ ± 1·10⁻⁴ K) against bl34 with the same gravity waves.
     Layer-mean temperature on the last day (K), global / 20S-20N /
     35-55N / 35-55S / 70-90N / 70-90S:

     | layer | July bl36 | July twin | January bl36 | January twin |
     |---|---|---|---|---|
     | 0.15 hPa | 230/228/232/231/245/213 | 230/227/232/235/244/206 | 229/227/228/234/208/243 | 229/227/230/234/210/242 |
     | 0.64 hPa | 253/255/257/249/270/229 | 254/254/258/255/270/222 | 253/254/250/259/220/270 | 253/253/251/259/222/270 |
     | 1.6 hPa | 257/261/265/246/274/221 | 257/261/266/247/274/216 | 257/262/244/266/212/278 | 257/261/245/266/214/278 |
     | 3.5 hPa | 244/245/253/233/260/211 | 244/247/254/229/260/207 | 244/247/230/252/205/269 | 244/247/230/252/207/268 |
     | 7.4 hPa | 229/230/235/222/246/201 | 229/232/237/218/245/199 | 229/231/219/235/199/253 | 229/231/218/236/201/253 |
     | 14 hPa | 221/222/225/217/238/195 | 221/223/226/213/237/193 | 220/222/212/225/196/243 | 220/222/212/225/197/243 |
     | 24 hPa | 216/215/219/215/234/193 | 216/217/220/211/234/191 | 214/214/210/218/196/237 | 214/214/210/218/197/238 |
     | 37 hPa | 211/209/214/214/231/193 | 211/210/214/211/231/192 | 209/206/210/213/197/233 | 209/206/208/212/198/235 |
     | 53 hPa | 208/204/211/214/228/194 | 208/205/211/211/228/193 | 205/200/209/209/198/230 | 205/200/207/209/199/231 |
     | 85 hPa | 204/198/207/213/222/196 | 204/198/207/210/224/195 | 203/197/207/207/197/222 | 203/197/206/207/200/224 |

     bl34, the same days: July 1.1 hPa 252/255/260/239/271/206, 3.5 hPa
     243/247/252/226/260/214, 14 hPa 220/223/225/211/237/200, 53 hPa
     208/204/211/211/229/196; January 1.1 hPa 251/255/240/261/211/271,
     3.5 hPa 243/246/228/252/217/266, 14 hPa 220/222/209/224/203/242, 53
     hPa 205/200/208/210/198/230. The winter cap over 0-2.2 hPa by mass:
     bl36 July 222 K (twin 216), January 214 (216); bl34 206 / 211; AFGL
     subarctic winter 247. Trends over days 21-30 (K/day): the global
     layers −0.19 to 0.00, the July 70-90S cap +1.06 / +0.79 / +0.49 at
     0.15-1.6 hPa and −0.12 to −0.83 at 3.5-101 hPa (twin +0.18 / +0.12 /
     −0.06 and −0.15 to −0.63), the January 70-90N cap +0.56 / +0.51 / −0.12
     and −0.44 to −0.68 at 3.5-53 hPa (twin +0.18 / +0.31 / −0.25, −0.27 to
     −0.75); bl34's caps −0.71 to −0.13 (July) and −0.30 to −1.00
     (January). The strongest zonal-mean westerly of the winter hemisphere
     (m/s at latitude), July / January, 0.15 / 0.64 / 1.6 / 3.5 / 7.4 / 14
     / 24 / 37 / 53 / 85 hPa: bl36 109@64S / 101@59S / 91@56S / 82@59S /
     73@59S / 63@56S / 55@56S / 50@54S / 47@49S / 49@46S, twin 125 / 106 /
     89 / 78 / 71 / 66 / 61 / 57 / 56 / 54; bl34 101@41S (1.1 hPa) / 66 /
     64 / 61 / 58 / 54 / 51 / 48; January bl36 108@61N / 97@54N / 82@51N /
     65@56N / 57@61N / 52@61N / 48@61N / 45@61N / 42@61N / 39@59N, twin
     104 / 93 / 77 / 62 / 55 / 51 / 49 / 45 / 44 / 43; bl34 88@46N / 56 /
     48 / 49 / 48 / 45 / 43 / 41. The 5S-5N wind at 0.15 / 0.64 / 1.6 /
     3.5 hPa on the last day: July −15 / +2 / −20 / −26 (twin −26 / −4 /
     −13 / −24), January −43 / −19 / +19 / −17 (−42 / −21 / +23 / −15); bl34
     at 1.1 / 3.5 hPa −31 / −22 and −16 / −27. Over the 30 days the top six
     layers' largest wind and horizontal / vertical Courant number: bl36
     July 194 m/s, 0.62 / 0.10 (twin 170, 0.54 / 0.10), January 176, 0.56 /
     0.13 (175, 0.56 / 0.13); bl34 117, 0.37 / 0.08 and 123, 0.39 / 0.08;
     by day (m/s / horizontal Courant, days 92-121): 122/0.39 117/0.37
     121/0.39 127/0.41 136/0.44 137/0.44 137/0.44 140/0.45 140/0.45
     138/0.44 134/0.43 131/0.42 166/0.53 131/0.42 142/0.45 130/0.42
     135/0.43 139/0.44 144/0.46 152/0.49 160/0.51 152/0.49 152/0.49
     194/0.62 151/0.48 145/0.46 152/0.49 149/0.48 139/0.45 136/0.44
     (July) and, days 275-304, 148/0.48 138/0.44 125/0.40 126/0.40
     113/0.36 109/0.35 122/0.39 130/0.42 137/0.44 139/0.45 138/0.44
     136/0.44 136/0.44 141/0.45 139/0.45 152/0.49 164/0.53 176/0.56
     174/0.56 167/0.54 167/0.54 162/0.52 162/0.52 149/0.48 138/0.44
     135/0.43 136/0.43 134/0.43 130/0.42 127/0.41 (January);
     their eddy kinetic energy at 1.6-7.4 hPa on the last day is 77-137
     m²/s² on bl36 in July (twin 58-102), 88-113 in January (70-111),
     where bl34's sponge holds 27-36 (July) and 16-27 (January) at
     3.5-7.4 hPa. The troposphere,
     means of days 1-10 and 21-30, balance / OLR (W/m²) / rain (mm/d):
     July bl36 −5.73 / 233.82 / 2.160 and −11.69 / 232.31 / 2.636, twin
     −5.66 / 233.81 / 2.160 and −12.17 / 232.81 / 2.677, bl34 −6.09 /
     234.15 / 2.166 and −11.06 / 232.77 / 2.639; January bl36 −4.32 /
     224.95 / 2.125 and −5.95 / 224.68 / 2.563, twin −4.37 / 224.90 / 2.129
     and −6.66 / 224.61 / 2.609, bl34 −4.68 / 225.18 / 2.127 and −7.06 /
     224.35 / 2.573: bl36 and bl34 differ by up to 1.1 W/m² and 0.01 mm/d,
     the twins by up to 0.7 and 0.05 (the earlier bl34 twins by 1.4 and
     0.07). Three N=128 days from eight128_day0183 remapped: largest wind
     114 / 87 / 86 m/s, Courant 0.37 / 0.28 / 0.28 at most in the top six
     layers, day 186 ASR 241.5, OLR 239.5 W/m², rain 1.85 mm/d, albedo
     0.291 (bl34 earlier 241.8 / 239.8 / 1.85 / 0.290).
   - Cost, `js/gpu/profile.module.js`, 128 steps under the exclusive lock,
     two runs each, the step median: N=64 bl34 24.69 / 24.73 ms, bl36
     26.96 / 27.14 ms (+9.5 %); N=128 bl34 (eight128_day0183) 99.85 /
     100.26 ms, bl36 112.78 / 113.09 ms (+12.9 %), 57.8 s of GPU steps a
     model day; the three N=128 days above took 72 s each with the
     `STRATOSPHERE` diagnostics. The physics pass with the gravity waves
     5.2 → 5.5 ms (N=64) and 20.5 → 21.9 ms (N=128); the 2 m/s phase
     speeds and the descending source moved bl34's N=128 step from 98.1
     to 100.1 ms.
   - Final defaults: `GRAVITY_WAVES` `flux` and `equatorialFlux` 4.3e-3
     Pa, `northFlux` and `southFlux` 0, `edge` 15°, `width` 10°,
     `sourcePressure` 31500 Pa with `sourceDescent`, `halfWidth` 40 m/s,
     `maxSpeed` 100 m/s, `speedStep` 2 m/s, `wavelength` 300 km,
     `minimumFrequency` 0.005 /s, `breakingAmplitude` 0.4 m²/s²,
     `lidPressure` 85 Pa; the sponge 1 day at the top from σ 0.005, on
     bl36 from 78 Pa; `topDragDays` 0; the longwave table by level set;
     the fresh stratosphere read by pressure above σ 0.1; bl36 for fresh
     spin-ups and the page's fresh start. What still misses: the winter
     polar stratopause, the caps 206-213 / 220-229 / 212-221 K at 0.15 /
     0.64 / 1.6 hPa against AFGL subarctic winter's 253 / 256 / 241 (bl36
     resolves the shape, not the warmth, with the intermittent breaking:
     the 247 / 257 / 236 K measured on bl36 earlier came with the grid-box
     test's doubled flux and its runaway top); the winter mesospheric jet
     at 0.15 hPa 104-125 m/s near 59-64° (no climatology read for it);
     the 0.64 hPa cooling 9-13 % above RRTMG's; the larger eddies and edge winds at 1.6-7.4 hPa with the sponge
     raised (Courant 0.62 at most at N=64); the fresh start's stratopause
     at 3.6 hPa and its first-week cloud (as on bl34); the flux through 20
     km 4.4-5.3 mPa, at the top of the observed 2-6, uniform in latitude
     and season. Whether bl36 can take a multi-year run: stable in both
     solstice months at N=64, from a fresh atlas start at N=64 and N=128
     and from a remapped N=128 state, with equatorial winds of the
     solstice's easterly sign on day 30 (their trend: the review below) and
     a troposphere within the twins of bl34's, at +9 to +13 % a step; its
     polar winter stratopause and mesospheric jets are not yet Earth's.
   The second round's review (Oct 2). Read here: MiMA's `cg_drag.f90`
   and `input/input.nml` (mjucker/MiMA, master), GFDL's `cg_drag.F90`
   (NOAA-GFDL/atmos_phys, main) and AM4's `run/input.nml`, Garfinkel et
   al. (2022) Appendix A (PMC9286580), Hertzog et al. (2008) in full, the
   abstracts of Corcos et al. (2021) and Geller et al. (2013), IFS CY48r1
   Part III section 2.2.11(b). Alexander & Dunkerton (1999) was again
   unreachable (403), as was Holt et al. (2017).
   - Confirmed in the sources: B_w 0.4 m²/s² (both codes, Garfinkel);
     c_w 40 m/s (GFDL code) and 35 (MiMA input.nml, Garfinkel); 300 km
     (`nk` 1 in both); c_max 99.6 m/s with Δc 1.2 (both codes), 2.4
     (AM4), 2 (Garfinkel: "the spectral resolution for the phase speed
     bins is 2 m/s"); 315 hPa at the equator descending as (K + 1) − (K +
     1 − k₀) cos φ, at most the second-lowest level (both codes); MiMA's
     4.3 mPa everywhere in input.nml and 4.3 + 3.5 mPa poleward of 15° in
     Garfinkel's CONTROL, with eq. A3 as `gravityWaveFlux` writes it; the
     breaking test and ε (MiMA divides by ρ₀ and drops GFDL's 1.5); MiMA
     tests the waves in every model level and spreads only what would
     leave the top evenly over the levels above 85 Pa
     (`damp_level_pressure`; Garfinkel: "deposited evenly in the levels
     above 0.85 hPa"), with no sponge; the IFS sponge from 0.78 hPa with
     the zonal mean undamped (Shepherd et al. 1996). Hertzog (2008): 2.5
     mPa raw, 3.2 corrected, about 6.4 an upper bound; ĉp 92 m/s an upper
     bound. Unconfirmed: AD99's own values, the HIRDLS 1.8-4.1 mPa
     (Geller's Table 1 through Holt), the Fleming et al. (1990) quotation.
     Corcos (2021)'s abstract gives no background value. Not as the
     sources: the source level is the nearest σ to 0.315^cos φ (the codes
     take the level above the first one below 315 hPa and truncate
     the index); the spectrum moves with u₀ (the codes keep a fixed
     ground-relative grid); no reflection test (MiMA drops waves with |c −
     u| k at or above N k / (k² + 1/(4H²))^½); the lid layers untested
     (MiMA tests them).
   - On real bl36 states (CPU, N=64: nine64_day0091 remapped, and the July
     run's day 121): each column's momentum is kept to 3e-16 of its
     absolute deposit; the kinetic energy one step takes (2.3 and 6.0
     mW/m²) returns as heat to 1.5e-8 of itself. The flux through the 45
     and 62 hPa interfaces is 4.5-5.6 mPa in every 10° band, the summer
     pole included (Geller 2013: the observed fluxes "are very small at
     summer high latitudes"). Source layers by latitude: 314 hPa near the
     equator and 15°, 314-369 near 30°, 369-510 near 45°, 510-601 near
     60°, 601-839 near 75°, 882-993 near the poles. Zonal-mean torque on
     day 121 (July, m/s/day): the two lid layers alike, −15 to −12 at
     50-80S, +8 to +9 at 10-40N, +3 to +6 within 10° of the equator;
     below them at most 1.9 (1.6 hPa), 1.4 (3.6 hPa) and 0.5 (7.5 hPa and
     down). With `lidTests` (MiMA's rule, added on both engines, 3.9e-7
     m/s/day rms apart at N=8) the top layer takes −24 to −25 at 50-70S
     and +15 at 20-40N.
   - The remap. `remapLevels` gave the three layers that split bl34's
     0-2.2 hPa layer its θ: N² = 0 across them and, on nine64_day0091,
     global means of 155 / 236 / 306 K at 0.15 / 0.65 / 1.6 hPa where the
     source layer holds 275 K; every remapped bl36 run above began so
     (day 92 of the July runs 162 / 242 / 289 K, the 70-90S cap 137 K at
     0.15 hPa). θ now remaps as θ·Π (the core's layer Exner function), so
     each column's enthalpy is kept to round-off, and a split layer's
     temperature is log-linear in σ with the minmod of its slopes to its
     neighbours, at most isentropic: 244 / 258 / 267 K above 272 K at 3.6
     hPa (day 274: 243 / 257 / 266). The column θ integral moves by 6.6e-4
     of its global value, all in the split layer; the round trip to bl34
     returns every field to 1.1e-12; mass, q, qc and u keep to 1e-14.
     Fresh starts do not remap: ten N=64 days log line for line as above.
   - Thirty July days (N=64) from nine64_day0091 so remapped, the
     defaults, and a round-off twin; the July run above reproduced log
     line for line. Top six layers' largest wind and horizontal / vertical Courant
     number 152 m/s, 0.49 / 0.10 (twin 185, 0.59 / 0.13, on day 120); the
     70-90S cap at 0.15 / 0.64 / 1.6 hPa 216 / 232 / 225 K (twin 210 / 226
     / 221), so the cold cap is not the remap's; days 21-30 balance −11.97
     / −10.27 W/m², OLR 233.47 / 233.03, rain 2.649 / 2.618 mm/d. The 5S-5N
     wind at 0.64 hPa went from −11 m/s to +17 (twin +40) on day 121,
     gaining 2.1 (twin 4.0) m/s/day over the last ten days; at 0.15 / 1.6
     hPa +8 / −31 (twin +5 / −29). The waves push both lid layers there at
     +4.1 to +5.5 m/s/day (the easterlies below filter the westward half),
     and nothing in an untested lid answers its own wind. The build's
     January runs (old remap) do the same at 1.6 hPa: −6 → +19 (twin +23) m/s, 1.8
     (2.1) m/s/day over days 295-304. With `lidTests` (one July run): the
     top layer's 5S-5N wind −16 → −93 m/s, falling 11 m/s/day at the end,
     0.64 hPa +53, 1.6 hPa +30, largest wind 176 m/s, as the build found
     for that rule.
   - Three N=128 days from eight128_day0183 so remapped: top six layers'
     largest wind 86 / 82 / 82 m/s, Courant 0.27 at most; day 186 ASR
     241.9, OLR 239.5 W/m², rain 1.85 mm/d, albedo 0.290 (bl34 241.8 /
     239.8 / 1.85 / 0.290).
   - Engines from the review's July day-121 state (bl36, N=64, no ocean,
     full physics), the top eight layers: after 1 / 4 / 16 steps the
     temperature 2.8e-3 / 3.6e-2 / 0.37 K apart at most (rms 5.5e-5 /
     2.7e-3 / 5.1e-2), the wind 1.4e-3 / 3.1e-2 / 0.44 m/s (rms 8.1e-5 /
     6.0e-4 / 2.7e-2), most in the top layer; below them 4.2 K and 2.8
     m/s at most after 16 steps (rms 1.3e-2 K, 9.0e-3 m/s).
   - bl34: the CPU model three steps from nine64_day0091 with `speedStep`
     4 and no `sourceDescent` against 9b7b476: 18645 of 8.5 million values
     apart, by 7e-11 Pa at most (the flux now scales a unit spectrum); the
     bl34 benchmark prints as 9b7b476's, the bl36 one as recorded above.
   - Cost, 128 steps, exclusive lock, step median: N=64 bl34 24.22 / 24.39
     ms, bl36 27.28 / 27.52 (+12.7 %); N=128 bl34 99.16 / 99.35, bl36
     111.79 / 112.33 (+12.9 %), 57.4 s of steps a model day.
   - Whether bl36 can take a multi-year run: thirty days hold (Courant
     0.59 at most), but the equatorial lid has not settled: under the
     default its 0.64 hPa wind was still gaining 2-4 m/s/day on day 30 in
     July (1.6 hPa in January, 2 m/s/day), and under MiMA's rule the top
     layer's easterly runs away instead. The winter cap stays 30-40 K
     below AFGL subarctic winter above 2 hPa from either start.
   The lid's budget and the lid friction (Oct 2, third round). N=64 GPU
   runs on bl36 one at a time, `STRATOSPHERE=1`, `TOP_BUDGET=8` (the top
   six layers' zonal-mean budgets sampled eight times a day:
   `createTopBudget` in `scripts/upperAtmosphere.mjs`, the gravity waves'
   force as the step applies it, the radiative heating as the physics
   pass applies it, the sponge's and the friction's by the CPU operators,
   the resolved terms of the Eulerian mean in σ; the residual is the
   total less the parameterized forces, so it holds the resolved
   dynamics, the closures and the dissipation heat). Thirty fresh days
   with it repeat `runs/igpre64.log` line for line, and days 31-60 from
   its day-30 snapshot repeat the rest (180 of 180 daily lines).
   - The lid's momentum from the fresh start, m/s/day at 0.15 / 0.64 hPa
     (the gravity-wave drag is the lid deposit, one acceleration over
     both layers):

     | band, days | total | gravity waves | resolved | Coriolis (f − ζ̄) v̄ | eddy u′v′ | vertical advection | eddy u′ω′ |
     |---|---|---|---|---|---|---|---|
     | 50-70S, 21-30 | +1.64 / +1.01 | −5.92 | +7.57 / +6.93 | +6.12 / +4.67 | +1.93 / +1.33 | +0.15 / +0.43 | −0.31 / −0.16 |
     | 50-70S, 51-60 | +3.50 / +2.12 | −11.84 | +15.40 / +13.97 | +16.05 / +11.84 | +1.77 / +1.09 | +0.40 / +0.86 | +0.87 / +0.40 |
     | 5S-5N, 21-30 | −1.17 / −0.20 | −0.01 | −1.13 / −0.19 | −0.54 / −0.33 | −0.42 / −0.36 | +0.01 / +0.16 | −0.07 / +0.08 |
     | 5S-5N, 41-50 | −4.85 / −1.85 | −3.93 | −0.99 / +2.08 | −2.19 / −1.23 | +0.32 / +0.22 | +0.44 / +1.92 | +0.04 / +0.72 |
     | 5S-5N, 51-60 | −1.36 / −1.03 | −5.24 | +3.82 / +4.20 | +1.70 / −1.62 | +0.22 / +1.32 | +0.88 / +3.04 | +0.33 / +0.49 |

     The sponge's zonal force is 0.07 m/s/day at most. The
     winter jet at the lid is driven by the Coriolis torque on the
     poleward mean flow (the radiatively driven summer-to-winter drift),
     which outgrows the lid deposit as the jet strengthens (63 → 127 m/s
     at 0.15 hPa, 50-70S, days 30 → 60). The top layer's equatorial
     easterly is driven by the waves: over days 41-50, −3.93 of its −4.85
     m/s/day, the 1.6 hPa westerly (+47 m/s, the waves' eastward deposit
     there +2.10 m/s/day) filtering the eastward half so that the westward
     half reaches the lid. Temperature, K/day at 0.15 / 0.64 / 1.6 hPa,
     days 51-60: the 70-90S cap radiates −4.33 / −5.74 / −3.61 and the
     dynamics returns +3.96 / +5.37 / +3.12 (cooling 0.38 / 0.37 / 0.49);
     20S-20N over days 41-50 +0.60 / +1.08 / +0.95 against −0.83 / −1.26 /
     −1.03 of ascent. The descent is too weak by the cap's 0.4-0.5 K/day
     because the drag that drives it is weak: the winter lid's force is
     −5.9 to −11.8 m/s/day, where the atmosphere's gravity-wave forcing
     in the winter reaches 100-160 m/s/day near 60° in the upper
     mesosphere (Sato et al. 2018, JAS 75, the abstract as a search quoted
     it, the paper not read), and the
     drag that closed the jet in models of the time was a 2 ± 1 day decay
     at 65 km, 1 at 70 and ½ at 75 km in winter (the GISS 21-layer model,
     Rind, Suozzo, Lacis, Russell & Hansen 1984, NASA TM-86183, read),
     with Holton & Wehrbein (1980) at 5 to 2 days over the same heights
     (as Rind et al. quote them; HW80 itself, Lindzen 1981, Holton 1982,
     1983 and Garcia & Solomon 1985 were not reachable). At the
     equatorial stratopause the observed gravity-wave forcing of the SAO
     peaks at about 5-7 m/s/day eastward and about 2 westward in
     reanalyses, the westward up to half the eastward on average in the
     middle mesosphere (Ern et al. 2021, ACP 21, 13763); the SAO's
     amplitude near the stratopause is over 30 m/s (Hirota 1978, through
     Kawatani et al. 2020, ACP 20, 9115), with westerlies all year at 0.1
     hPa in SABER and MLS winds (Kawatani 2020; both read through a
     summary of the paper), and easterlies at the solstices and
     westerlies at the equinoxes near the stratopause (Fleming et al.
     1990, CIRA-86, read). The lid's westward
     −5.2 m/s/day at the equator (days 51-60) is about two and a half
     times the westward forcing observed at the stratopause. CIRA-86 puts the
     winter jet's maximum in the midlatitude mesosphere ("Maximum
     velocities generally occur in the midlatitude mesosphere"), so a jet
     that does not close between 1.6 and 0.15 hPa is not by itself a
     misfit; its speed is: 146 m/s at 60S on day 60.
   - Built (both engines): a Rayleigh friction on the zonal-mean wind,
     the sponge's band mean, u ← ū/(1 + r̄ dt) + (u − ū)/(1 + r dt), the
     kinetic energy it removes returned as heat through the closure's
     dissipation (`dampEddies`, the core's `spongeMeanRates`, the GPU's
     `spongeApply`); its rate a profile of decay time against log-pressure
     height (H = 7 km) averaged over each layer's mass (`LID_FRICTION`,
     `lidFrictionRates`, `surface.lidFriction`): `holtonWehrbein` 5 days
     at 65 km to 2 at 75 km, `rind` 2 / 1 / 0.5 days at 65 / 70 / 75 km,
     none below 65 km, the last held above. On bl36 both reach the 0-0.3
     hPa layer alone (9.30 and 2.71 days; nothing at 0.64 hPa or below),
     on bl34 its 0-2.19 hPa layer (68.7 and 20.0 days). The friction does
     not conserve the column's momentum: on the July day-121 state it
     removes the top layer's axial angular momentum at 1/(2.69 days)
     (2.2·10¹⁷ of 5.1·10²² kg m²/s a second), the layers below untouched,
     and returns 8.0 mW/m² as heat (each column's, measured from its θ
     change, equal to its kinetic energy lost to 2.5·10⁻¹² of the largest
     column's loss: the review below); one step's change of the top layer's
     wind (up to 0.10 m/s) agrees between the engines to 8.3·10⁻⁶ m/s
     (1.4·10⁻⁶ rms); `test/sponge.test.mjs`: torque 1.0001 of a Rayleigh
     drag on a zonal flow with the wave moved by 1.3·10⁻³ m/s rms, heat to
     10⁻⁹ of the kinetic energy removed with the sponge, 20 N=8 steps
     1.7·10⁻³ m/s apart of the 6.58 m/s the treatment moves.
   - The candidates, thirty days each from nine64_day0091 and
     nine64_day0274 remapped (days 92-121 and 275-304) and days 31-60
     from the fresh start's day-30 snapshot, defaults (no friction) →
     `holtonWehrbein` → `rind`, round-off twins (the lowest layer's θ ±
     10⁻⁴ K) in brackets; 0.15 / 0.64 / 1.6 hPa unless stated:

     | | July day 121 | January day 304 | fresh day 60 |
     |---|---|---|---|
     | 5S-5N u | −2/+15/−35 [−4/+5/−24] → 0/+30/−31 → +2/+53/−30 [+1/+55/−25] | −18/−16/0 [−18/−15/−9] → −13/−14/+3 → −8/−14/+8 [−8/−16/+15] | −100/−27/+42 → −53/−32/+42 → −23/−30/+43 |
     | winter 60° u | 105/87/75 [97/84/73] → 93/83/74 → 67/74/72 [67/68/65] | 57/30/18 [50/34/24] → 59/34/20 → 58/38/23 [44/35/24] | 146/98/70 → 91/76/60 → 68/62/53 |
     | winter jet | 106@59S/91@49S/81@36S → 95@64S/85@54S/78@56S → 70@64S/76@59S/75@56S | 67@34N/81@34N/72@31N → 69@51N/67@34N/65@31N → 58@61N/61@31N/60@29N | 146@59S/102@56S/76@54S → 91@59S/85@49S/73@49S → 68@59S/72@39S/65@39S |
     | winter cap, 70-90° | 212/226/222 [214/228/224] → 219/234/228 → 234/249/241 [233/248/239] | 218/227/228 [223/232/232] → 218/226/228 → 227/233/232 [232/238/234] | 210/223/218 → 220/232/224 → 233/245/236 |
     | winter cap, 14 / 24 / 37 / 53 hPa | 196/193/193/194 [197/195/195/196] → 198/195/194/195 → 202/198/196/196 [200/196/195/195] | 215/212/210/208 [213/209/207/205] → 216/214/213/212 → 217/213/210/208 [215/210/207/206] | 197 (14 hPa) → 198 → 202 |
     | summer cap, 70-90° | 245/270/274 [245/271/274] → 241/267/272 → 235/262/270 [235/262/270] | 242/270/278 [242/270/278] → 238/266/275 → 233/262/273 [232/261/272] | 244/270/270 → 242/268/269 → 238/264/267 |
     | 20S-20N T | 229/256/262 → 228/256/262 → 227/255/261 | 226/252/260 → 226/252/260 → 226/251/260 | 220/247/257 → 226/250/257 → 226/249/258 |
     | global T, 0.15-14 hPa | 230.0/253.6/257.0/243.6/229.2/220.6 → 230.2/253.4/257.3/243.8/229.3/220.7 → 230.0/253.4/257.8/244.1/229.4/220.8 | 229.2/252.6/257.0/243.9/229.0/220.1 → same → 229.2/252.5/257.5/244.0/229.1/220.2 | 227.7/252.8/256.6/242.5/227.4/219.0 → 229.1/252.6/256.8/242.5/227.4/219.0 → 229.0/252.7/257.4/242.7/227.6/219.1 |
     | top six layers' largest wind, Courant, last day | 132, 0.42 [133, 0.42] → 155, 0.50 → 144, 0.46 [135, 0.43] | 118, 0.38 [103, 0.33] → 119, 0.38 → 95, 0.30 [95, 0.30] | 156, 0.50 → 123, 0.39 → 95, 0.30 |
     | the same, largest of the run | 140, 0.45 [142, 0.45] → 155, 0.50 → 161, 0.52 [146, 0.47] | 155, 0.50 [148, 0.48] → 132, 0.42 → 126, 0.40 [127, 0.41] | 159, 0.51 → 124, 0.40 → 112, 0.36 |

     AFGL's layer means at 0.15 / 0.64 / 1.6 / 3.5 hPa: subarctic winter
     250 / 256 / 241 / 227, subarctic summer 233 / 274 / 275 / 261,
     tropical 235 / 266 / 265 / 252 (`scripts/standardAtmospheres.mjs`
     on bl36). Trends over the last ten days (m/s/day; K/day): July
     `rind` 5S-5N +0.48 / +3.52 / +0.16 [+0.18 / +3.65 / +0.32], 60S +1.56
     / +2.15 / +2.05 [−0.33 / +0.42 / +0.47], the 70-90S cap +0.02 / −0.02
     / −0.18 [+0.05 / +0.06 / −0.20] (defaults −0.31 / −0.26 / −0.39);
     January `rind` 60N +4.38 / +3.66 / +3.03 [+2.92 / +3.83 / +4.09]
     (defaults +5.79 / +3.80 / +3.44), the 70-90N cap −0.64 / −1.07 /
     −1.28 [+0.01 / −0.21 / −0.92] (defaults −1.20 / −1.55 / −1.46); fresh
     days 51-60 `rind` 5S-5N −0.29 / −1.88 / −0.73, 60S +0.75 / +0.52 /
     +0.74, the 70-90S cap +0.20 / +0.08 / −0.15 (defaults −0.95 / −0.81 /
     −0.60, +5.21 / +2.67 / +1.38, −0.41 / −0.39 / −0.53); global layer
     means −0.18 to +0.04 K/day in every run. Thirty fresh days (to day
     30): 5S-5N 0/+10/+34 → −14/+3/+36 → −6/+8/+30, 60S 67/49/34 →
     48/35/23 → 43/42/35, the 70-90S cap 213/231/230 → 221/238/237 →
     229/244/238, the top six layers' largest wind and Courant number
     over the run 99 m/s, 0.32 → 92, 0.30 → 79, 0.25.
   - The lid's budget under `rind`, m/s/day at 0.15 / 0.64 hPa: fresh days
     51-60, 50-70S friction −23.13 / 0, waves −8.85, Coriolis +33.79 /
     +8.46, total +0.67 / +0.53; July days 112-121 friction −22.16, waves
     −10.25, Coriolis +30.01 / +6.19; the 70-90S cap's descent +8.02 /
     +9.55 / +5.32 K/day against radiation −7.94 / −9.53 / −5.48 (fresh
     days 51-60). At the equator the friction holds the top layer
     (+8.03 m/s/day at −23 m/s against the waves' −5.26, fresh days
     51-60), and in July it does not reach the 0.64 hPa layer, whose
     westerly the waves push at +4.28 / +4.47 m/s/day [twin] with the
     Coriolis torque −1.91 / −1.97 against the defaults' −3.06 / −3.26
     (days 92-121): +53 / +55 m/s on day 121, where the solstice's
     stratopause is easterly.
   - Troposphere, means of the first and the last ten days, balance /
     OLR (W/m²) / rain (mm/d): July defaults +1.98 / 234.26 / 2.362 and
     −6.07 / 238.46 / 2.695 [+1.97 / 234.30 / 2.359 and −6.23 / 238.23 /
     2.774], `rind` +2.03 / 234.27 / 2.363 and −6.07 / 238.81 / 2.714
     [+2.01 / 234.26 / 2.365 and −6.03 / 238.71 / 2.693]; January
     defaults +8.69 / 223.53 / 2.166 and +3.12 / 229.33 / 2.672 [+8.67 /
     223.53 / 2.163 and +2.90 / 229.63 / 2.751], `rind` +8.62 / 223.56 /
     2.161 and +2.13 / 229.79 / 2.710 [+8.57 / 223.57 / 2.163 and +3.40 /
     229.79 / 2.640]; fresh days 21-30 defaults +2.13 / 229.92 / 2.875,
     `rind` +3.46 / 229.23 / 2.796. `rind` and the defaults differ by up to
     1.0 W/m² in balance, 0.5 in OLR and 0.04 mm/d over the last ten days,
     the twins by up to 1.3, 0.3 and 0.08.
   - Cost, `js/gpu/profile.module.js`, 128 steps after 64 under the
     exclusive lock, twice each, the step median: N=64 (nine64_day0091 on
     bl36) 27.7 / 27.9 ms without the friction, 27.8 / 27.8 with it;
     N=128 (eight128_day0183 on bl36) 116.0 / 116.6 → 116.3 / 117.3 ms
     (+0.4 %), 59.8 s of steps a model day. The friction runs in the
     sponge's pass, on layers it already covers.
   - Final defaults: `lidFriction` `rind` on bl36 (`lidFrictionFor`), none
     on bl34 and cam26; the sponge, the gravity waves and their lid (85
     Pa) as above. The validation the readiness rests on, a fresh atlas
     start to day 150 at N=64 with the defaults (`TAG=<new> N=64 DAYS=150
     LEVELS=bl36 STRATOSPHERE=1 TOP_BUDGET=8 OCEAN='{"everySteps":8}'
     node scripts/spinup.mjs`), was not run: its launch was refused by the
     session's permission check, though the user's approval of N=64 runs
     of up to 150 days for the model top stands. It would test whether the
     top six layers' winds and temperatures level or oscillate within
     bounds over days 100-150 (through the June solstice and a month past
     it): the top layer's equatorial wind against the SAO's range, the
     0.64 hPa equatorial westerly that the July runs build to +53 m/s, the
     winter jet at 60S (68 m/s and +0.75 m/s/day on day 60), the winter
     cap's approach to AFGL, and the Courant number's margin (0.30 on day
     60, 0.52 at most in the July runs). Verdict: bl36 is not shown ready for a
     multi-year run from a fresh start. The friction removes the fresh
     start's runaway at the lid (day 60: top-layer equatorial wind −23 m/s
     against −100, the 60S jet 68 against 146 m/s, Courant 0.30 against
     0.50, the winter cap rising 0.08-0.20 K/day at 0.15-0.64 hPa against
     falling 0.4) and brings the July winter cap to AFGL at 1.6 hPa, but
     the 0.64 hPa equatorial wind in July is +53 / +55 m/s and still
     gaining 3.5 m/s/day, the winter vortex at 14 hPa is 4-6 K warmer in
     July (202 / 200 K against 196 / 197), and nothing beyond day 60 of a
     fresh start has been seen.

   The third round's review (Oct 2). Read here: Rind, Suozzo, Lacis,
   Russell & Hansen (1984, NASA TM-86183), the model description and its
   Fig. 1.
   - Confirmed in the source: the drag in layers 19-21 ("approximately
     65-75km", mean pressures 0.09 / 0.05 / 0.03 mb in Fig. 1, so H = 7 km
     puts the profile's 65 / 70 / 75 km at 0.094 / 0.047 / 0.023 hPa, on
     those layers); the winter decay times 2 ± 1, 1 ± ½ and ½ ± ¼ day,
     "varying with wind speed, and thus latitude"; Holton & Wehrbein's 5
     to 2 days over the same heights. Not as the source: GISS's drag has
     the surface drag's form with a stability-dependent coefficient on the
     whole wind, where `rind` is a linear decay of the zonal mean at those
     times, the eddies left to the sponge; GISS's final version discards
     the energy the drag removes (returned as heat, its 70N winter
     mesosphere ran 10-30 K warmer, at variance with observations), where
     the model returns it. On the July day-121 states below that heat is
     0.45-0.82 K/day in the top layer at 70-90S and 1.0-1.5 at 50-70S (two
     members), where the cap's top layer radiates −7.9 (fresh days 51-60). The rate above 75 km,
     held at ½ day, is 43 % of the 0-0.3 hPa layer's (its 0-2.4 Pa); the
     rate is set per σ layer at 1013 hPa, so over the Antarctic plateau
     (680 hPa, the layer at 0-0.20 hPa) the profile by pressure would give
     1.82 days for the 2.71.
   - Profile: 2.709 days on bl36's 0-0.3 hPa layer (an independent
     10⁶-point mean 2.7085), `holtonWehrbein` 9.30, zero in every other
     layer of bl36, bl34 and cam26; the sponge 1.24 and 6.00 days at 0.15
     and 0.65 hPa.
   - Accounting on nine64_day0091 remapped to bl36 (N=64), from the state
     changes: the CPU closure and dissipation with the friction alone take
     14.155 mW/m² of kinetic energy and return 14.155 as heat, each
     column's heat (cp Π Δθ times its mass) within 2.5·10⁻¹² of the
     largest column's loss (22.6 J/m²) of its kinetic energy from the
     winds before and after; with the sponge 16.628 and 5.1·10⁻¹². The top
     layer's axial angular momentum 3.09·10²² kg m²/s falls by 1.33·10¹⁷
     a second (2.69 days), the layers below untouched. GPU, one full step
     with the friction against one without: 14.226 mW/m² both ways, each
     column to 1.2·10⁻³, the angular momentum −1.335·10¹⁷ a second, every
     other layer equal. The friction's own change of the wind, 9.3·10⁻²
     m/s at most, GPU less CPU on the same pre-closure wind 8.3·10⁻⁶ at
     most (1.2·10⁻⁶ rms).
   - Engines from the same state (full physics, no ocean), the top eight
     layers after 1 / 4 / 16 steps: T 1.5·10⁻⁴ / 3.6·10⁻² / 0.37 K apart
     at most (rms 2.7·10⁻⁵ / 2.7·10⁻³ / 4.9·10⁻²), u 1.6·10⁻³ / 0.15 /
     0.60 m/s (rms 8.5·10⁻⁵ / 9.9·10⁻⁴ / 2.7·10⁻²), as in the second
     round; the friction itself moves the top layer's wind by 9.3·10⁻² /
     0.37 / 1.9 m/s (rms 3.4·10⁻² / 0.14 / 0.56). bl34: the CPU model three
     steps from nine64_day0091 equals 027d2c6's byte for byte.
   - Tests: 63 of 63 files, 580 tests. dayMeans' two columns: 324's plume
     runs on one engine only at steps 17 and 18; 135's keeps its top, its
     base flux 4.287·10⁻² against 4.504·10⁻² at step 7 with the boundary
     layer's depth (2575.12 / 2575.07 m) and buoyancy flux agreeing; on the
     CPU ×(1 + 10⁻⁵) on the column's vapour moves that flux by 1.3·10⁻⁴ of
     itself and ×(1 + 10⁻⁴) by 5.2 % (to 4.510·10⁻²), a jump. gpuModel: the
     CPU's ±1 ulp response 2.43·10⁻⁴ K/day again in 20 and 200 draws, the
     limit 4.0·10⁻⁴ (1.5 × 1e-5 × 26.6).
   - Logs: ten-day linear trends reproduce the tables above; `top3fb`'s
     180 daily lines are `igpre64`'s. A July run with the defaults
     (`top3rvjul`) repeats `top3rijul` line for line, so neither
     `TOP_BUDGET` nor naming the friction changes the run; a third member
     (`top3rvjulT`, θ ± 10⁻⁴ K, another draw): day 121 5S-5N +2 / +47 /
     −32, 60S 54 / 61 / 63, the 70-90S cap 236 / 252 / 245 and 202 K at 14
     hPa, the summer cap 234 / 262 / 269, the largest wind and Courant
     number over the run 165 m/s, 0.53, the 0.64 hPa equatorial wind
     gaining 3.93 m/s/day over the last ten days. Against the members'
     spread, `rind` against the defaults: in July the 0.64 hPa equatorial
     westerly (+47 to +55 against +5 / +15), the 60S wind at 0.15 hPa (54
     to 67 against 97 / 105), the winter cap (+15 to +26 K) and the summer
     cap at 0.15-0.64 hPa (−10 / −8 K) differ, and 14 hPa by 3-6 K; the
     July Courant number does not fall (0.47-0.53 against 0.45 / 0.45), the
     largest wind sitting at 1.6 hPa (130-144 m/s on day 121). In January
     the 60N wind lies within the twins at 0.15 and 1.6 hPa and 4 m/s
     above them at 0.64, the Courant number falls (0.40 / 0.41 against
     0.50 / 0.48), the winter cap rises 4-14 K at 0.15 hPa against a twin
     spread of 5. The fresh "day 60" runs take the
     defaults to day 30 and the candidate after; with `rind` from the start
     30 days exist (`top3rif`: 5S-5N −6 / +8 / +30, 60S 43 / 42 / 35, the
     cap 229 / 244 / 238, Courant 0.25 at most). On day 60 the 1.6 hPa
     equatorial westerly is +42 / +42 / +43 in all three (the defaults'
     +42 to +51 since day 35), below the friction's reach.
   - The 150-day fresh validation was refused again by the session's
     permission check and not retried. Verdict: bl36 is not ready for a
     multi-year run from a fresh start. No run with the friction lasts more
     than 30 days; the July 0.64 hPa
     equatorial westerly grows in every member, the friction's July
     Courant margin is no better than the defaults', and the friction's
     treatment of the removed energy is the opposite of its source's.
2. The deck gate. The vertical mass flux smoothed over neighbouring
   cells before it is interpolated to the deck height (the page's
   overlay already does this), the memory shortened from ten days to
   about two, the inversion jump made the primary regime test, the
   height bound checked (the SE Pacific deck sits at 730 m under a
   resolved inversion at 1330 m). Separately, an ablation of weak
   divergence damping in the core against its eddy cost.

   Built. Both engines average πσ̇ with equal weights over the cell and
   its neighbours, twice over (`subsidenceSmoothing` 2), at the two
   interfaces bracketing h, and the running mean remembers two days
   (`subsidenceMemory`); the smoothed value feeds the mixed layer's
   dh/dt as well. On the day-183 N=128 state the SE Pacific's
   deck-height sink keeps its box mean (1.52 → 1.57 mm/s) while its
   per-cell spread falls from 29.7 to 8.6 mm/s (× 3.45) and the
   grid-scale share of its variance from 0.97 to 0.12; Peru's spread
   falls from 17.5 to 6.7 mm/s about 0.93 → 1.12 mm/s, inside the raw
   box mean's grid-noise uncertainty of 1.2 mm/s (N=64: 13.1 → 5.0 mm/s
   about 2.86 → 2.48). After five days the saved running mean is
   3–8 % grid-scale in the two boxes, where the ten-day unsmoothed one
   was 55–85 % on day 183 and 69–99 % after the same five days.
   `scripts/verticalAudit.mjs` prints the deck-height sink as the
   dynamics leaves it and as the deck reads it, replicates the gate from
   the physics phase's own inputs (exactly, in running mean, gate and
   decision) and attributes each column-step: the deck runs, or is off
   on the subsidence, on the jump, or on the gate's memory. On the
   day-183 state, with its saved ten-day mean, the SE Pacific deck ran
   on 0.089 of the column-steps and 0.303 failed the 0.3 mm/s floor;
   after five days of the smoothed two-day mean (run `p2a128`, the floor
   kept) the mean sink at h is 1.76 mm/s in the SE Pacific and
   2.12 mm/s in Peru with a spread over the cells as large, the floor
   still refuses 0.19 and 0.22 of the column-steps, and the jump test
   fails on 0.82 and 0.50. The regime test is therefore the jump alone
   and the subsidence only vetoes ascent faster than 1 mm/s
   (`stratusSubsidence` −1 mm/s, about twice the mean's grid-scale
   residual), which on day 188 of that run refuses 6 % of the SE Pacific
   and 5 % of Peru against 51 % of the warm pool (10S–10N, 120–170E) and
   44 % of the east Pacific ITCZ (5–12N, 140–90W), by area of the saved
   running mean (6, 6, 53 and 42 % in the package gate's own run).
   Five-day N=128 continuations from day 183, audited on day 188 (SE
   Pacific; the saved running mean's box mean, spread and grid-scale
   share):

   | day 188, SE Pacific | ten-day, unsmoothed, 0.3 mm/s | package gate | package gate, c = 0.03 |
   |---|---|---|---|
   | saved sink, mm/s (spread, grid share) | 1.47 (2.46, 0.99) | 1.76 (1.71, 0.08) | 2.17 (1.57, 0.04) |
   | deck runs, share of column-steps | 0.056 | 0.074 | 0.160 |
   | failing the subsidence test | 0.276 | 0.059 | 0.026 |
   | failing the jump test | 0.810 | 0.811 | 0.641 |
   | deck's virtual jump at its start, K | 1.42 | 1.42 | 1.78 |
   | deck height where it runs, m | 790 | 866 | 821 |
   | resolved inversion, m (θv jump, K) | 1035 (2.75) | 1054 (2.78) | 1158 (3.37) |
   | low cloud | 0.050 | 0.051 | 0.083 |
   | rain, mm/d | 0.64 | 0.52 | 0.57 |
   | Peru: deck runs (low cloud) | 0.272 (0.086) | 0.372 (0.072) | 0.510 (0.137) |
   | ω700 grid-scale share, global | 0.642 | 0.632 | 0.181 |

   On the day-183 state itself the start height (735 m), the height where
   the deck runs (714 m under the old gate, 713 m under the new), the
   resolved inversion (1328 m) and its jump (3.44 K) are the same under
   both gates: an eight-step window cannot move the one-day gate memory.

   The height. Where the SE Pacific deck is off (0.91 of the
   column-steps on day 183) its start height is the Richardson depth by
   construction: it relaxes there over `heightMemory`. Where it runs it
   stands at 714 m, within 20 m of that floor on 0.54 of the
   column-steps and at the ceiling on 0.03, never above the resolved
   inversion: its mixed layer is cloud-free on 0.93 of them (cloud base
   706 m for h 713 m, 2.1 g/m² of water), so the radiative closure
   entrains 0.48 mm/s against 1.48 mm/s of subsidence and the layer
   sinks to the floor; Peru is cloud-free on 0.98. The ceiling
   (1063 m where the deck runs, lowest single-interface 2 K jump) lies
   above the midpoint of the resolved inversion's upper layer (1006 m)
   in some columns,
   where a smeared inversion splits its jump over two interfaces, but it
   binds on 3–6 % of the running column-steps, and replacing it by the
   largest-gradient interface would favour the thin lowest layers. So
   neither `bound`, the ceiling nor the height memory holds the deck
   down: the entrainment–subsidence balance of a layer whose condensation
   level sits at its top does, and the dry subcloud layer is the shallow
   Betts–Miller firing's (item 3). Nothing in the bound was changed.
   After five days the running SE Pacific decks are still cloud-free on
   0.80 of their column-steps (0.90 with the damping) and within 20 m of
   the floor on 0.44 (0.72), at the ceiling on 0.06 (0.02) and never
   above the resolved inversion.

   Divergence damping (`divergenceDamping`, c): each step the closure
   adds c d² ∇(∇·u) to the edge velocities, the tendency ν_d ∇δ with
   ν_d = c d²/dt and d the mean distance between cell centres, after
   the ∇⁴ momentum closure in both engines; the kinetic energy it removes
   is returned as heat with the closure's, it leaves the mass and the
   vorticity untouched, and the engines agree to 3·10⁻⁴ m/s over twenty
   N=8 steps. Ten days at N=64 from day 183 (day 193):

   | c | ω700 grid-scale share | global mean KE, J/kg | EKE, J/kg |
   |---|---|---|---|
   | 0 | 0.557 | 150.5 | 49.2 |
   | 0.01 | 0.320 | 149.0 (−1.0 %) | 50.6 (+2.7 %) |
   | 0.03 | 0.183 | 147.7 (−1.9 %) | 50.4 (+2.4 %) |
   | 0.05 | 0.143 | 146.5 (−2.7 %) | 49.3 (+0.1 %) |

   (EKE from the cell-reconstructed winds less their 2° zonal means,
   mass-weighted.) The eddy cost is nil within the runs' spread, so
   `createModel` and `createGpuModel` default to c = 0.03
   (`DIVERGENCE_DAMPING`); the SE Pacific's ten-day rain was 1.01, 0.99,
   1.19 and 1.34 mm/d across the four runs, which item 3 has to watch.
   At N=128 (the last column of the table above) it takes the global
   ω700 grid-scale share from 0.63 to 0.18 in five days, the deck-height
   sink's spread as the dynamics leaves it from 18.5 to 4.6 mm/s, and it
   sharpens the SE Pacific inversion (3.37 against 2.78 K, EIS 3.11
   against 2.65 K), so that the deck runs twice as often.
3. Convection. A boundary-layer-mean parcel with virtual temperature
   and an entraining ascent; a trigger on dilute CAPE and inhibition
   with a short memory; adjustment from cloud base up only, with a
   share of the rain evaporating in the subcloud layer; a shallow,
   non-precipitating mixing-line branch for tops below about 700 hPa
   and no firing under an active deck; autoconversion kept out of the
   lowest layers and evaporation below cloud base, to stop the
   ±10 K/d condensation–evaporation churn in the lowest 100 m.
   Sounding tests (a stratocumulus profile must not fire, a deep
   tropical one must, with peak heating at 400–500 hPa and none below
   cloud base), conservation, CPU–GPU parity, and regression on the
   audited state: SE Pacific rain under 0.3 mm/d, firing under 1 % a
   step, low cloud above 0.4.

   Built, in both engines (`js/physics/moist.module.js`, the adjust
   kernel of `js/gpu/physics.gpu.js`). The parcel is the lowest layer's
   air (`boundaryParcel` false, `parcelDepth` 0; with `boundaryParcel` the
   mass-weighted mean θ and q of the layers below the boundary layer's
   Richardson depth or of the lowest `parcelDepth`); it rises dry to its
   Bolton LCL, then saturated with its moist static energy relaxed toward
   the air's at `entrainmentRate` 5·10⁻⁵ m⁻¹, buoyant in virtual
   temperature. Its inhibition is the negative buoyant energy from the
   top of its source layers to the first buoyant layer above the LCL,
   its CAPE the positive energy above that, its top the highest buoyant
   layer; tops below `shallowTop` 700 hPa are shallow. The deep pass is
   the product of two ramps, 0 to 1 as the CAPE rises through
   `capeThreshold` 100 J/kg ± half of it and as the inhibition falls
   through `inhibitionThreshold` 50 J/kg ± half, and of the deck's
   opening, 1 at a mixed-layer gate of one half or less and 0 at 0.6; the
   per-cell `convectiveActivity` relaxes toward it over `activityMemory`
   2 h (saved with the state, one half in older ones) and the column
   fires while it is above one half. A shallow top also vents, without
   memory, at the same product of ramps about `shallowCape` 10 J/kg and
   `shallowInhibition` 15 J/kg times the deck's opening (and not where
   the estimated inversion strength exceeds `shallowStability`, off by
   default), relaxing over the relaxation time divided by its vent. Only
   the layers from the LCL's layer to the top relax, toward the parcel's
   temperature and 60 % of its saturation, shallow tops too
   (`shallowReference` 'parcel', `shallowRain`; 'mixingLine' without
   `shallowRain` is the non-raining Betts mixing line of the first
   build, its water at most `shallowHumidity` 0.8 of saturation); the
   anvil keeps 10 % of the rain and `downdraftEvaporation` 0.01 of the
   rest may evaporate into the subcloud layers, offered to each in
   proportion to its mass (`downdraftSpread` 'mass'). Autoconversion
   stays out of the lowest two layers (`autoconversionFloor` 'lowest';
   'boundaryLayer' keeps it out of every layer wholly inside the
   boundary layer), and rain evaporates only into cloud-free layers (at
   most 10⁻⁷ kg/kg of cloud water). Jordan's (1958) mean
   hurricane-season sounding heats most at 440 hPa with its top at
   195 hPa. The
   radiation gives resolved cloud a cover (`cloudCover` 'pdf', both
   engines): a layer covers the part of a uniform total-water
   distribution of half-width (1 − RHc) qs above saturation, RHc 0.85
   inside the boundary layer and 0.8 above, with emissivity
   f (1 − exp(−a W / f)), and the shortwave blends the clear column with
   the column whose path lies in the column's cover, the layers' largest
   f times their visibility 1 − exp(−W / 1 g/m²) (overlapped maximum-random
   and bounded by the condensate under strong inversions since the Arctic
   cover below). A Sundqvist cover from
   the vapour's humidity would be overcast: the radiation sees the state
   after the dynamics and before the saturation adjustment, and on the
   defaults' day-186 N=64 state 94 % of the cloud water below σ 0.9 lies
   in layers at 99 % humidity or more, a water-weighted Sundqvist cover
   of 0.94. The tests (`test/convection.test.mjs`): a stratocumulus
   column over a 26 °C sea under a 1.5 K θv inversion at 1.3 km with dry
   air above never convects; the Jordan sounding convects from cloud
   base with no change below it and its heating peak at 440 hPa; a trade-
   cumulus column vents at an activity of 0, at half the rate at its
   CAPE threshold, and not without the trigger, below half its threshold,
   above its stability bound or under a deck; column enthalpy and water
   close to 10⁻¹⁶ through the rain, anvil and downdraft, and heat and
   water to 10⁻¹⁶ in the mixing-line branch; the activity's switching
   times, both autoconversion floors, and the engines agree on 362
   random columns to 10⁻⁴ K and 5·10⁻⁸ kg/kg under the defaults, the
   first build's options and the stability veto. `test/physics.test.mjs`
   checks the cover's blend against the overcast column with W / f.
   `test/gpuModel.test.mjs` compares the engines' layer heating under
   resolved cloud of partial cover inside and above the boundary layer
   (1.1·10⁻⁴ K/day against 23 K/day) and the buoyancy-closure deck under
   the default convection, venting included (cover 7·10⁻⁵); the random
   columns include gates on the deck's opening ramp. Both engines take
   the deck's clear-column sunlight from the column blended at the
   resolved cloud's cover; the GPU first took it from the overcast
   column, and over a column that had vented the deck's cover parted by
   3.6·10⁻³.
   `scripts/verticalAudit.mjs` adds, from the moist physics' `trace`,
   the pressure of the maximum of the Pacific ITCZ firing columns'
   convective heating and its mean over the lowest 100 m, the box's
   large-scale heating below 1 km, the global rain, the convective share
   of 15S–15N and the fraction of columns convection changes at all;
   `MOIST` passes moist options to it and to `scripts/spinup.mjs`. The
   GPU model mirrors the boundary layer's depth and spin-up states carry
   it (`boundaryDepth`): the physics reads the depth of the step before,
   and a continued run's first step read 0 (no deck, no boundary-layer
   humidity for the cover); a coupled run stopped inside a day now
   continues within 2.2 % of the uninterrupted run's daily forcing over
   two days, where it was within 6.5 %.

   As first built, with the boundary layer's dilute parcel, the
   trigger and a shallow branch that never rained, convection stopped at
   once on the audited states (8-step windows before and after; the
   dilute parcel has no CAPE in a troposphere adjusted to the 40 m
   layer's undilute adiabat):

   | | eight128 d183 | eight64 d183 | eight128 d274 | eight64 d365 |
   |---|---|---|---|---|
   | SE Pacific rain, mm/d | 1.36 → 0.00 | 0.92 → 0.00 | 0.28 → 0.04 | 1.13 → 0.15 |
   | its columns firing a step | 0.065 → 0 | 0.042 → 0 | 0.007 → 0 | 0.036 → 0 |
   | Pacific ITCZ rain, mm/d | 2.69 → 0.42 | 4.86 → 0.57 | 2.55 → 0.41 | 1.89 → 0.28 |
   | global rain, mm/d | 2.36 → 0.48 | 2.61 → 0.55 | 2.77 → 0.59 | 2.66 → 0.50 |
   | ITCZ firing columns' heating peak, hPa | 973 → none | 438 → none | 439 → none | 439 → none |
   | ITCZ large-scale heating below 1 km, max \|K/d\| | 12.3 → 13.9 | 13.1 → 14.1 | 11.2 → 12.5 | 10.7 → 11.9 |

   Five days of that build at N=128 from eight128_day0183 (run
   `m21d128`) failed: the planetary albedo rose from 0.315 to 0.574 (ASR
   234 → 145 W/m²), the cloud water path from 78 to 177 g/m², and the SE
   Pacific drizzled 0.70 mm/d from a layer 0.975 overcast. Forty-nine
   three-day N=64 runs from eight64_day0183, audited on day 186,
   attributed it: the pre-package scheme gives 0.323 and 15.4 g/m² of
   cloud water below σ 0.9, the build 0.568 and 57.3 g/m². The
   boundary-layer parcel, the trigger and the non-raining shallow
   branch each remove the old scheme's continuous venting of the
   boundary layer on their own (added one at a time to the old scheme
   0.466, 0.503, 0.479; removed one at a time from the build 0.525,
   0.555, 0.564); every other element moves the albedo by less than
   0.04. Moving water inside the column (a non-raining shallow branch
   that fires at once, from the surface or not, entrainment across the
   boundary-layer top at 5–15 mm/s, a mass-spread downdraft) leaves the
   humidity below σ 0.9 in 30S–30N at 0.80–0.83 against the old scheme's
   0.70; only raining it out, or the cover, moves the albedo, and the
   cover alone takes 0.135 off the build and 0.057 off the old scheme.
   Hence the defaults above: the lowest layer's parcel, shallow tops
   venting on their own trigger and raining, the deep trigger kept, and
   the cover. Three days at N=64 (day-186 albedo): the build 0.568, the
   decision 0.319. Variants, each against a run of the code it was
   measured with: a shallow CAPE threshold of 5 J/kg 0.313 against 0.327
   at 10, one of 20 J/kg 0.348 against 0.331; a vent ramp from the
   threshold to twice it instead of about it 0.330 against 0.319, with
   no fewer SE Pacific firings after the three days; a stability veto on
   the vent at an EIS of 1, 2 or 3 K 0.378, 0.365, 0.358 against 0.327
   (at 2 K it vetoes a third of the venting ocean columns worldwide and
   only half of the SE Pacific's, whose venting columns have a median EIS
   of 2.2 K), so it stays an option; shallow tops raining along the
   mixing line instead of toward the parcel 0.333 against 0.331, with SE
   Pacific firing 0.105 against 0.059. In the audit window on the
   decided configuration's day-186 N=64 state a downdraft share of 0.1,
   0.05, 0.02 or 0.01 changes the ITCZ firing columns' lowest 100 m by
   −48, −32, −15 or −6 K/d (they rain tens of mm/d into a subcloud layer
   a few hundred metres deep), 0 by +4.6, hence 0.01.

   Five days at N=128 on the GPU from eight128_day0183 (run `r128c`,
   audited on day 188), against the old scheme (package 2's gate with
   c = 0.03, run `p2c128`) and the first build (`m21d128`):

   | day 188, N=128 | old convection | first build | the decision |
   |---|---|---|---|
   | planetary albedo, days 184–188 | 0.323, 0.329, 0.337, 0.338, 0.339 | 0.481, 0.526, 0.554, 0.571, 0.574 | 0.296, 0.311, 0.320, 0.328, 0.331 |
   | ASR, W/m², day 184 → 188 | 230.6 → 225.1 | → 145 | 239.7 → 227.7 |
   | cloud water path, g/m² (below σ 0.9) | 78 (15.5), day 183 | 177 | 113 (36.8) |
   | SE Pacific rain, mm/d (convective share) | 0.57 | 0.70 (0.00) | 0.81 (0.81) |
   | its columns firing a step (adjusting at all) | | 0.000 (0.006) | 0.100 (0.142) |
   | its low cloud | 0.083 | 0.975 | 0.493 |
   | its deck runs, share of column-steps | 0.160 | 0.049 | 0.181 |
   | its resolved inversion, m (θv jump, K); EIS, K | 1158 (3.37); 3.11 | 1186 (6.13); 5.03 | 994 (3.31); 3.59 |
   | Pacific ITCZ rain, mm/d (convective share); ω500, Pa/s | | 5.68 (0.78); −0.018 | 5.64 (0.80); −0.009 |
   | ITCZ firing columns' heating peak, hPa (K/d); lowest 100 m, K/d | | 515 (46.4); −71.3 | 438 (29.0); −0.3 |
   | ITCZ large-scale heating below 1 km, K/d | | +18.1 at 793 m, −31.8 at 21 m | +15.5 at 437 m, −34.0 at 21 m |
   | global rain, mm/d (convective share); 15S–15N share | | 2.53 (0.49); 0.79 | 2.70 (0.71); 0.86 |
   | zonal-mean rain peak, mm/d (latitude) | | 4.58 (7.5N) | 5.54 (9.5N) |

   With the GPU deck's sunlight blended at the cover (run `rv128f`, the
   same five days) the decision's numbers move by little: albedo 0.296,
   0.312, 0.320, 0.328, 0.331; ASR 239.8 → 227.9 W/m²; SE Pacific rain
   0.81 mm/d (0.83), firing 0.102 (0.146), low cloud 0.488, deck runs
   0.186; ITCZ peak 438 hPa (29.0 K/d), lowest 100 m −0.45 K/d,
   large-scale +15.2 at 437 m and −32.9 at 21 m; global rain 2.69 (0.71).

   Ten days at N=64 from the atlas (bl34): the old scheme's albedo peaks
   at 0.397 on day 7 and ends at 0.339, the first build's climbs to
   0.591 and ends at 0.580, the decision's peaks at 0.401 and ends at
   0.335. Against the acceptance: the daily albedo stays within
   0.30–0.34 but for day 184 (0.296), and the ASR falls 12.0 W/m² from
   day 184 where 10 was asked (the old scheme 5.5, from a first day
   0.027 brighter); the ITCZ heats aloft and its lowest 100 m by −0.3 K/d;
   the SE Pacific's low cloud (0.49) and deck runs (0.18) pass; the global
   rain passes; the fresh start passes. The SE Pacific's rain
   (0.81 mm/d) and firing (0.100 a step) fail, by more than under the
   old scheme (0.57 mm/d; 0.065 a step on day 183): the venting that keeps the boundary layer from saturating
   worldwide vents the SE Pacific wherever the deck's gate is shut, and
   the deck runs on 0.18 of its column-steps because its virtual jump at
   h is 1.7 K against the resolved inversion's 3.3 K. The ITCZ's
   large-scale dipole (15.5 and −34.0 K/d) fails too: splitting it on the
   decided configuration's day-186 N=64 state, the saturation adjustment
   condenses 8–20 K/d between 240 and 590 m and
   evaporates cloud at −8 to −11 K/d between 60 and 170 m, and the rain
   evaporates at −24 K/d into the lowest layer; neither the rain
   evaporation (0.1–1 of the deficit), autoconversion (no floor, no
   lifetime), the adjustment's base nor the boundary layer's mixing of
   cloud water changes the condensation peak.

   The Arctic cover (Sept 30). The paired spin-up "nine" melts its
   northern summer ice under the resolved cloud's cover: on
   nine128_day0091 the 70–90N ice's stratus holds about 0.18 g/kg of
   cloud water in saturated layers against a half-width of 0.67 g/kg, so
   each layer covers about 0.64 and the column 0.69 by maximum overlap,
   and over a day of sun positions on the frozen state the ice absorbs
   97.9 W/m² of sunlight against 39.8 under overcast; three days at N=64
   from nine64_day0091 lose 0.227·10³ km³ of northern (60–90N) ice a day,
   0.143 under overcast. Both engines now overlap the resolved cloud
   maximally within each run of adjacent cloudy layers and randomly
   between runs (`cloudOverlap` 'maximumRandom'); where the column's EIS
   (the deck's, from the lowest layer's air) rises through
   `overcastInversion` 8–12 K they blend a layer's cover into that of a
   distribution whose half-width is also at most the layer's cloud water
   and not below `overcastWater` 5·10⁻⁵ kg/kg, so that under a strong
   inversion a saturated layer holding more is overcast. A shallow CAPE
   threshold of 5 J/kg was tried with them and dropped: a ten-day N=64
   run from eight64_day0183 gave the same albedo and balance as the
   threshold of 10 and twice the SE Pacific firing. The bound without the EIS gate
   overcasts nearly all resolved cloud: 92.7 % of its water lies in
   layers at 99.9 % humidity or more, and the Arctic's water-weighted
   qc/qs (0.077) is below the globe's (0.131); gating it by temperature
   instead (layers below 0 °C) reaches the Arctic (0.941) at a larger
   global cost (0.501). On nine128_day0091 after one CPU step at the
   defaults (column cover by area; the ice's absorbed sunlight over 24
   hourly sun positions; the frozen state's planetary albedo, 1 − ASR /
   insolation summed over four sun positions 6 h apart from the step):

   | | maximum (package 3) | maximum-random | + bound, no gate | + bound, EIS 8–12 K | overcast |
   |---|---|---|---|---|---|
   | 70–90N ice columns | 0.692 | 0.777 | 0.957 | 0.945 | 0.961 |
   | global | 0.450 | 0.469 | 0.589 | 0.479 | 0.699 |
   | 15S–15N | 0.469 | 0.490 | 0.649 | 0.490 | 0.745 |
   | SE Pacific (where the deck is off) | 0.334 (0.345) | 0.339 (0.350) | 0.524 (0.539) | 0.339 (0.350) | 0.643 (0.657) |
   | ice's absorbed sunlight, W/m² | 98.0 | 78.1 | 40.1 | 42.3 | 39.8 |
   | global albedo, frozen state | 0.393 | 0.401 | 0.459 | 0.406 | 0.469 |

   Confining the gate to the layers below σ 0.7 changes the 70–90N ice
   columns' cover from 0.945 to 0.938 and the global from 0.479 to 0.478.

   Three days at N=64 on the CPU from nine64_day0091 (`arc3d64`) lose
   0.173·10³ km³ of northern ice a day (volume 9.19 → 8.67), above the
   0.14–0.16 asked; the 70–90N ice's surface takes 64.9, 71.7, 69.3 W/m²
   net in the daily means (package 3 91.8, 93.1, 88.6; overcast 53.6,
   58.4, 57.6). The N=64 state's Arctic inversions are weaker: on
   nine64_day0091 the ice columns cover 0.883 (maximum 0.720, overcast
   0.963) and the ice absorbs 60.8 W/m² (92.7, 49.5); a ramp of 4–8 K
   would take it to 52.5 at N=64 but the N=128 global cover to 0.490
   (6–10 K: 55.5, 0.483), and the bound without the gate reaches only
   51.0. Ten days at N=64 on the GPU from eight64_day0183 (`arc10d64`):
   planetary albedo 0.287, 0.307, 0.318, 0.323, 0.334, 0.340, 0.340,
   0.341, 0.343, 0.338 on days 184–193 (the last physics step's, as the
   log prints it), the first three within 0.005 of package 3's three-day
   run from the same state (`j0`: 0.288, 0.312, 0.319), then above the
   0.30–0.31 asked; ASR − OLR +1.9 on day 184, −13.2 W/m² at its lowest
   on day 192 and −11.0 on day 193; global rain 2.04, 2.30, then
   2.31–2.58 mm/d. No ten-day run of package 3 from that state exists to
   set the later days against. `scripts/verticalAudit.mjs` on day 193, against
   package 2's ten-day run with c = 0.03 (`dd364`, day 193) and package
   3's three-day run (`j0`, day 186):

   | | this, day 193 | package 2, day 193 | package 3, day 186 |
   |---|---|---|---|
   | SE Pacific rain, mm/d (convective share) | 1.48 (0.96) | 0.30 (1.00) | 0.98 (0.95) |
   | its columns firing a step (convecting) | 0.142 (0.229) | 0.053 (0.085) | 0.079 (0.131) |
   | its low cloud | 0.390 | 0.060 | 0.356 |
   | its deck runs, share of column-steps | 0.219 | 0.320 | 0.365 |
   | its EIS, K; deck's virtual jump, K | 2.07; 1.52 | 2.43; 1.94 | 3.44; 2.30 |
   | Pacific ITCZ rain, mm/d (convective share) | 6.96 (0.87) | 4.48 (0.98) | 4.20 (0.72) |
   | ITCZ firing columns' heating peak, hPa; lowest 100 m, K/d | 439; −1.39 | 961; +17.7 | 439; +0.17 |
   | global rain, mm/d | 2.59 | 1.94 | 2.37 |

   The SE Pacific's columns fire on 0.142 of the column-steps on day 193
   against package 3's 0.079 on day 186; the two differ in day and in all
   three changes, so the rise is not attributed to the vent threshold.
   Five days at N=128 on the GPU from
   nine128_day0183 (`arc5d128`): albedo 0.296, 0.296, 0.299, 0.296, 0.304
   on days 184–188 against the spin-up's own 0.301, 0.302, 0.304, 0.303,
   0.310; ASR − OLR +1.8, +1.1, +0.9, +0.8, −2.2 W/m² (mean +0.5; the
   spin-up's +0.1, −0.7, −0.9, −1.0, −3.9); the northern ice extent
   0.072 → 0.090 Mkm² (cells at least 15 % covered).

   Tests: a synthetic column's single run of cloudy layers covers as its
   largest layer, two runs two clear layers apart 1 − (1 − f1)(1 − f2),
   and maximum overlap the largest; a saturated layer holding 0.02 qs
   reads the uniform cover below EIS 8 K, the linear blend on the ramp
   and cover 1 above 12 K; the GPU parity test of partial cover adds
   maximum overlap, the unbounded half-width, and a bound of 5·10⁻⁴ kg/kg
   on a −40 to 40 K ramp and over every column, on columns with separate
   runs (1.5·10⁻⁴ against 24 K/day). At the default 5·10⁻⁵ kg/kg the
   bounded half-width amplifies the engines' f32 and f64 difference in
   q + qc − qs: over every column the engines then part by 5.0·10⁻⁴
   K/day, on the ramp 3.3·10⁻⁴. The scattering-only overcast
   engine's digests are re-pinned for the vent threshold alone.
4. The ocean's wind response, diagnosed on the day-274 state before
   it is changed: stress against 0.05 N/m² on the equator, mixed-layer
   depth against 30–50 m in the east, whether an undercurrent exists
   in the thermocline classes, the momentum sinks. Acceptance: a
   0.2 m/s westward surface current east of 140W, an undercurrent of
   0.5–1 m/s near 100 m, an eastern thermocline at 40–60 m, and 2 K
   between the warm pool and the cold tongue held through a year.
5. One-year N=64 runs from the atlas start, audited at day 183 with
   the tools of item 1, before a new paired spin-up. A mass-flux
   convection scheme stays deferred until the central Pacific ITCZ and
   the Walker cell still fail with all of the above in place.
6. Boundary-layer entrainment (Sept 30). The K-profile vanishes at h,
   so outside the deck nothing entrained. Both engines give the first
   interface above h the coefficient ρ w_e,
   w_e = min(cap, (A B0 + A_s u*³/h) / max(Δb, b_min)), with B0 the bulk
   surface buoyancy flux, Δb = g Δθv/θv between the layer above h and the
   boundary layer's mass mean, A 0.2 and A_s 5 (the buoyancy and
   friction-velocity sources of Tennekes 1973 with the constants of
   Driedonks 1982), b_min 0.015 m/s² and cap 0.05 m/s (`entrainment`:
   `efficiency`, `shear`, `jumpFloor`, `cap`), and zero where B0 ≤ 0 or
   the deck's gate exceeds one half. Over the sea at 30S–30N where the
   deck is off, after one CPU step: nine128_day0183 7.74 mm/s area mean
   (median 4.28, 90th percentile 20.1), 2.24 mm/d of water carried out of
   the boundary layer; eight64_day0183 2.94 mm/s, 1.14 mm/d. A floor of
   3·10⁻³ m/s² gave 11.6 mm/s on nine128 (0.070 of the columns at the cap,
   Δθv 0.2–0.3 K in the tail) and at N=64 an albedo of 0.274, 0.291,
   0.296 on days 184–186 against 0.274, 0.292, 0.297 with 0.015.
   Relative humidity below σ 0.9 at 30S–30N, mass-weighted, sea
   (all surfaces): nine128_day0183 0.808 (0.747), eight64_day0183 0.738
   (0.685). Three days at N=64 on the GPU from eight64_day0183, day 186:

   | | off | A 0.15 | A 0.2 | A 0.3 |
   |---|---|---|---|---|
   | w_e, mm/s; water out of the layer, mm/d | 0 | 4.49; 1.21 | 5.19; 1.37 | 6.53; 1.67 |
   | humidity below σ 0.9, sea (all) | 0.809 (0.741) | 0.807 (0.736) | 0.807 (0.736) | 0.807 (0.736) |
   | cloud water below σ 0.9, sea 30S–30N, g/m² | 46.8 | | 37.5 | 36.8 |
   | albedo, days 184–186 | 0.293, 0.319, 0.328 | 0.275, 0.293, 0.299 | 0.274, 0.292, 0.297 | 0.272, 0.288, 0.293 |

   The layers inside the boundary layer (σ 0.95–1) dry from 0.86–0.88 to
   0.84–0.86 and those at σ 0.88–0.91 moisten by 0.02–0.03, so the mean
   below σ 0.9 stays at 0.807 for every A in 0.15–0.3. Ten days at N=64
   on the GPU from eight64_day0183 (days 184–193):

   | | off | A 0.2 (default) | A 0.3 |
   |---|---|---|---|
   | albedo | 0.293, 0.319, 0.328, 0.331, 0.343, 0.348, 0.349, 0.348, 0.351, 0.348 | 0.274, 0.292, 0.297, 0.300, 0.311, 0.319, 0.322, 0.324, 0.326, 0.319 | 0.272, 0.288, 0.293, 0.297, 0.306, 0.312, 0.315, 0.320, 0.318, 0.315 |
   | ASR − OLR, W/m², days 188–193 (mean) | −13.6 to −15.7 (−14.9) | −3.7 to −8.5 (−6.8) | −2.2 to −6.6 (−4.8) |
   | global rain, mm/d | 2.03–2.56 | 2.00–2.59 | 2.00–2.58 |
   | humidity below σ 0.9, sea (all), day 193 | 0.821 (0.762) | 0.830 (0.765) | |

   `scripts/verticalAudit.mjs` on day 193 (off; A 0.2): SE Pacific rain
   1.72; 1.44 mm/d, firing 0.176; 0.137, low cloud 0.402; 0.223, deck runs
   0.232; 0.190, deck's jump 1.59; 1.29 K, EIS 2.13; 2.04 K, resolved
   inversion's jump 3.78; 4.06 K; Pacific ITCZ rain 6.43; 6.10 mm/d, ω500
   −0.023; −0.026 Pa/s, heating peak 439 hPa in both (24.8; 24.1 K/d,
   lowest 100 m −1.53; −1.77 K/d), large-scale +15.4 at 436 m and −23.0
   at 21 m; +13.4 and −19.4; global rain 2.57; 2.53 mm/d; zonal-mean peak
   5.46 at 1.5S; 4.91 at 0.5S. The spin-up log's ten-day SE Pacific rain
   0.93; 0.80 mm/d (`base10d64`, before the Arctic cover, 0.89), Pacific
   ITCZ 3.95; 4.12 (4.27). Five days at N=128 on the GPU from
   nine128_day0183 at A 0.2: albedo 0.278, 0.272, 0.272, 0.268, 0.269
   (`arc5d128` 0.296–0.304), ASR − OLR +7.5, +8.7, +9.1, +9.5, +8.8 W/m²
   (mean +8.7; `arc5d128` +0.5), northern ice extent 0.072 → 0.097 Mkm²,
   w_e 6.72 mm/s and humidity below σ 0.9 0.808 (0.744) on day 188.
   Three days at N=64 from nine64_day0091, 60–90N ice volume loss per
   day: GPU 0.172 (off 0.168)·10³ km³, CPU 0.173 (off 0.167; volume
   9.19 → 8.67 and 8.69).
   In the first CPU step from eight64_day0183, over the sea at 30S–30N
   where the deck is off, the shear term gives 1.62 and the buoyancy term
   1.33 of the 2.94 mm/s (each uncapped) and is the larger on 0.58 of the
   entraining columns; from day 193 of the A 0.2 run, 2.23 and 3.46 of
   5.69 mm/s. In that first step 1113 columns have a buoyancy term below
   5% of the shear term; their w_e, which drops to zero where B0 changes
   sign, has a median of 3.8 and a 90th percentile of 11.1 mm/s. The step
   keeps each column's mass-weighted θ and water and each edge's momentum
   to 1·10⁻¹⁵ relative; because θ rather than cp T is mixed, the
   global-mean enthalpy gains 0.18 W/m² with entrainment off and 0.24
   with it on. Continued from day 193 in one-day segments (`rv4on`,
   `rv4off`), days 194–199, A 0.2: albedo 0.314, 0.319, 0.317, 0.313,
   0.314, 0.316, ASR − OLR −4.0, −5.1, −3.6, −2.3, −1.8, −2.2 W/m²; off:
   0.347, 0.354, 0.351, 0.339, 0.346, 0.350 and −13.9, −16.0, −14.2,
   −10.6, −11.9, −12.1. The ocean reaches its 5 m/s speed limit at
   1S 99–100E on days 193 (4 edges) and 194 (3) with entrainment and on
   days 195 (4) and 196 (3) without it.

**Item 4, diagnosed and tried (Sept 30).** `scripts/equatorialOcean.mjs`
takes a saved state apart on the CPU by 20° of longitude along the
equator: wind and stress with the implied drag coefficient, mixed-layer
depth and how often it sits on its floor, the 20 °C isotherm and the
class tops, the zonal current by class and by depth, a meridional
section at 140–110W, the momentum budget term by term (the terms sum to
the ocean's own tendency within 10⁻¹⁹ m/s²) for the mixed layer, the
water above the 1024 class and the whole column, and the closure's pull
near the thermocline classes' token edges. On the N=128 day-183 and
day-274 states and the N=64 day-365 one:

- The stress reaches the ocean intact (the ice factor is 1 on the
  equator; sampling it once an ocean step rather than averaging the
  eight atmosphere steps changes it by under 3 %), with an effective
  10 m drag coefficient of 1.8×10⁻³ against Large and Pond's 1.2×10⁻³,
  so the drag law is not the cause. The trades are weak: −3.2 to
  −3.7 m/s at 10 m over 180–120W and −2.1 m/s at 120–100W, the stress
  falling from −0.048 N/m² at 180–160W to −0.021 at 120–100W (day 183),
  about 60 % of Earth's, and the meridional stress at 140–110W northerly
  or nil against Earth's southerly 0.02–0.04 N/m².
- Above the 1024 class the pressure force balances 0.4–1.6 of the
  stress: the tilt (20 °C from 172 m to 90 m, sea level down 26 cm from
  150E to 110W) is in balance with the weak stress it is given. East of
  120W the stress falls to a third while the sea-level slope goes on, and
  the mixed layer flows east at 7–22 cm/s.
- The mixed layer sat on its 50 m floor in 43–75 % of the cells at
  120–80W, spreading the stress over 50–67 m.
- The interfacial drag capped the shear under the mixed layer at
  τ/(ρ₀ r), 0.15 m/s at 0.03 N/m² (measured 4–6 cm/s), and, divided by
  at least 50 m in every layer, made the eastward column source that
  **Drag** (M18) describes.
- The thermocline classes 1022.5–1024.75 hold more than 5 m in only
  20–60 % of the equatorial cells and have tokens on 55–91 % of the
  edges. Their ∇⁴ closure reads the tokens' velocity, the westward mixed
  layer's, and removed −1 to −3.8×10⁻⁷ m/s² per class at N=128 against an
  eastward pressure force of +1 to +3.4×10⁻⁷ (−0.3 to +0.6×10⁻⁷ away
  from the tokens, where 1024.00 flows east at 9 cm/s). The eastward flow
  that exists is 10–17 cm/s, broad (3S–7N) and at 100–200 m.

Changed in both engines, with CPU–GPU parity tests: the drag divisor
(**Drag**, the only change to the defaults); the option `shearMixing`
(**Drag**); and the option `closureFill`, the share of the least-squares
uniform flow through a class's own neighbouring edges that its token
edges beside it take in the closure's input (`closureVelocity`). The
floors (`minimumThickness`, `shallowestMixedDepth`) and `interfacialDrag`
were already options. 60 coupled days at N=64 from `eight64_day0365`
(March 20 to May 19), 2S–2N at day 60, the second row of a pair a
second realization (the surface diffusivity 1 % higher); W−E and the
cold tongue as the spin-up log prints them, clamps as edge-days (days):

| Configuration | u(5 m) 140–100W | u(5 m) 160E–100W | strongest eastward, 160–120W | h₀ 120–100W | 1024 top 120–100W | W−E | cold tongue | clamps |
|---|---|---|---|---|---|---|---|---|
| before | −0.05, −0.05 m/s | −0.13, −0.14 | 0.04 m/s at 195 m, 0.06 at 185 | 52 m | 83, 83 m | 1.9, 2.2 K | 25.7, 25.5 °C | 7 (3), 2 (1) |
| the drag divisor (the defaults now) | −0.13, −0.14 | −0.12, −0.18 | 0.03 at 195, 0.02 at 265 | 51 | 76, 71 | 2.0, 2.4 | 25.4, 25.6 | 2 (2), 0 (0) |
| 20 m floors, r = 2×10⁻⁴ | −0.04 | −0.11 | 0.02 at 190 | 22 | 80 | 2.0 | 26.3 | 3 (2) |
| 20 m floors, r = 10⁻⁴ | −0.25 | −0.48 | 0.01 at 285 | 22 | 77 | 1.1 | 26.5 | 6 (4) |
| 20 m floors, r = 5×10⁻⁵ | −0.77 | −0.97 | 0.03 at 245 | 24 | 60 | 2.0 | 25.4 | 5 (2) |
| the same, closureHours 48 | −0.74 | −0.88 | 0.06 at 195 | 24 | 67 | 2.1 | 25.7 | 7 (4) |
| 20 m floors, shearMixing, r ≥ 0 | −0.29 | −0.67 | 0.10 at 235 | 22 | 71 | 0.6 | 26.9 | 4 (4) |
| 20 m floors, shearMixing, r ≥ 5×10⁻⁵ | −0.18, −0.38 | −0.58, −0.52 | 0.02 at 220, 0.07 at 205 | 22 | 76, 68 | 1.3, 1.6 | 26.6, 26.6 | 17 (6), 18 (6) |
| 20 m floors, shearMixing, r ≥ 0, closureFill 0.5 | −0.36 | −0.64 | 0.15 at 195 | 21 | 68 | 1.6 | 26.1 | 60 (13) |

The drag divisor alone moves the surface current east of 140W from
−0.05 to −0.13 m/s, lifts the eastern class top by 7–12 m and keeps W−E
and the cold tongue, with fewer clamped edges. No configuration makes an
undercurrent: the strongest eastward flow stays 1–15 cm/s at 185–285 m.
The weaker the friction the stronger the South Equatorial Current,
which with a constant 5×10⁻⁵ m/s was −0.24 m/s at
140–100W at day 30 and −0.77 at day 60, still accelerating, and with
the Richardson form −0.13 to −0.21 at day 30 and −0.18 to −0.38 at day
60. The stronger current pulls the patchy thermocline classes west
through the closure: in the Richardson run at day 60 the closure beside
the tokens is −0.6 to −5.2×10⁻⁷ m/s² in 1022.0–1025.0 against −0.1 to
−1.4×10⁻⁷ before, and the whole column loses 1–3×10⁻⁵ m²/s² to it at
160E–140W, as much as the stress gives; the interfacial drag's column
sum, 0.3–1.1×10⁻⁵ eastward at 160E–80W before, is now zero. The 20 m
floor puts the eastern mixed layer at 21–24 m and lifts the eastern
class top by 3–23 m, but the thin layer rests on warm water above a
thermocline that is still 70–80 m deep, and the cold tongue warms by
0.4–1.2 K except in the two runs with the strongest current. The two
realizations of a configuration differ by up to 0.4 K in W−E and 0.2 m/s
in the surface current, and 60 days of boreal spring, when Earth's cold
tongue is weakest, cannot show whether W−E holds through a year. Some
edges clamp in nearly every run, the one before included; those examined lie in
the Maritime Continent's seas, where by day 60 the run before had a
12.7 °C mixed layer at 1S 131E.

Filling the token edges entirely (`closureFill` 1, with a second ring
filled from the first) removed the closure's pull but also its hold on a
class's edges beside its tokens, where the Coriolis force, assembled
from the neighbours' fluxes, is weak: within a week 588 layer edges ran
above 1.2 m/s and edges clamped at 5 m/s, with column transports of up
to 165 Sv against 25; filling only the token edges with three or more
class neighbours did the same. Half the fit held for 60 days with the
strongest undercurrent of the set, within about 2S–3N, but clamped
edges on 13 days.

So the defaults keep the 50 m floor and r = 2×10⁻⁴ m/s. Before the
ocean can take the 20 m floor and the weaker friction, the atmosphere
must give it Earth's trades in the east (items 2–3) so that the eastern
thermocline rises under the thin layer, and the thermocline classes need
the westward pull of their token edges removed without losing the hold,
for example by filling the patchy 1022.5–1024.75 classes with a weak
diapycnal exchange between adjacent interior classes.

**Item 6, a shallow cumulus mass flux (Sept 30).** Both engines
(`cumulusColumn` in `js/physics/moist.module.js`, the adjust kernel of
`js/gpu/physics.gpu.js`) replace the shallow Betts–Miller vent
(`shallowScheme` 'massFlux'; 'bettsMiller' keeps the vent) with a bulk
plume after Bretherton, McCaa and Grenier (2004). Source: the
mass-weighted s_l = c_p T + g z − L q_c and q_t of the boundary layer, or
of the lowest `cumulusSourceDepth` 50 hPa where that is deeper, lifted
unmixed to its Bolton LCL (`cumulusSource` 'lowest': the lowest layer's
air). Base mass flux ρ_LCL c w exp(−CIN / w²), c = `cumulusClosure`
0.06, w = max((B₀ h)^⅓, `cumulusFriction` 1 × u*) from the boundary
layer's surface buoyancy flux and depth (kept per cell in both engines,
`buoyancyFlux`, `friction`), CIN the negative buoyant energy from the
source's top to the first saturated layer; zero where B₀ ≤ 0, where the
LCL lies above `shallowTop` 700 hPa, under the deck's gate (the deep
branch's ramp from 0.5 to 0.6), below 10⁻⁶ kg/m²/s, and, without
`cumulusWithDeep`, where the deep branch fires; at most
`cumulusBoundaryLoss` 0.1 of the source's mass a step. Above the LCL
the plume entrains at `cumulusEntrainment` 2.5·10⁻³ m⁻¹ and detrains at
`cumulusDetrainment` 3·10⁻³ m⁻¹, buoyant in virtual temperature with
condensate loading, and ends in its first non-buoyant cloudy layer or
the last below 700 hPa, `cumulusOvershoot` 1 of the flux entering that
layer detraining there. s_l and q_t move in flux form, M (X_u − X_above)
at each interface, M g Δt/Δp ≤ 1 in every layer, followed by a
saturation adjustment; no rain unless `cumulusRain` (kg/kg of plume
condensate). The deep branch acts on tops above 700 hPa only. The
radiation (`cumulusCloud`, both engines) gives each plume layer the
cover max(f, M/(ρ `cumulusUpdraft`)), `cumulusUpdraft` 1 m/s, and adds
that fraction of the plume's condensate. Tests
(`test/convection.test.mjs`): a trade-wind column (26 °C sea, 1.5 K θv
inversion at 1.3 km) lifts 0.042 kg/m²/s with no inhibition into the
1285–1814 m layer, which moistens by 55 g/kg/d, drying the boundary
layer by 7.9 kg/m²/d, with enthalpy and water to 2·10⁻¹⁶ and nothing
raining; a stable surface, no flux, a deck and an LCL above the shallow
top give none, half the flux on the deck's ramp; the Jordan column's
deep branch is bit-identical under either scheme; on 362 random columns
(98 plumes) the engines agree on which columns lift and where the plumes
end, on the base flux to 5·10⁻⁵ of 0.059 kg/m²/s, θ to 1.8·10⁻⁴ K and q
to 3·10⁻⁷; `test/gpuModel.test.mjs`: the cumulus cover moves the layer
heating by up to 6.0 K/day and the engines agree to 9.8·10⁻⁵ K/day.
The moist adjustment of every column after one CPU step from
eight64_day0183 (8153 plumes) and twin64_day0540 (9930) keeps each
column's c_p T + L q to 1.1·10⁻¹⁵ and its water with the rain to
7.7·10⁻¹⁶ of themselves, the winds and surface pressure untouched.

On nine128_day0183 after one CPU step, over ice-free sea in 30S–30N
where the gate is at most 0.5: mean base flux 0.0077 kg/m²/s, 0.174 of
the columns above 10⁻³ with a mean of 0.0446 there (Earth 0.02–0.05);
plume tops at the 950, 900, 850 and 800 hPa interfaces on 0.30, 0.52,
0.16 and 0.015 of them; SE Pacific 0.0155 (0.32 of its columns). With
c 0.04, ε 1.5·10⁻³, δ 2·10⁻³ and plumes beside deep convection 0.0239,
0.77, 0.031; with those parameters kept out of deep columns and a
buoyant saturated layer required below 700 hPa (CIN to the LFC) 0.0013,
with 0.043 of the columns above zero and 0.72 never buoyant.

Ten days at N=64 on the GPU from eight64_day0183 (planetary albedo on
days 184 / 187 / 193, ASR − OLR on day 193), plumes beside deep
convection: c 0.04, ε 1.5·10⁻³, δ 2·10⁻³ 0.266 / 0.281 / 0.279 (+8.7);
c 0.02 0.271 / 0.285 / 0.294 (+4.1); c 0.03, overshoot 0.5 0.273 /
0.287 / 0.295 (+4.1); `cumulusUpdraft` 0.3 0.270 / 0.285 / 0.284
(+7.1). Kept out of deep columns: the same 0.271 / 0.291 / 0.310
(−1.1); c 0.03 0.272 / 0.298 / 0.322 (−5.0); c 0.06 0.269 / 0.291 /
0.303 (+1.0); overshoot 0.5 0.274 / 0.299 / 0.330 (−7.5); overshoot
0.75 0.271 / 0.293 / 0.315 (−2.9); ε 10⁻³, δ 1.5·10⁻³ 0.263 / 0.289 /
0.306 (−0.1); ε 2·10⁻³, δ 2.5·10⁻³ 0.280 / 0.296 / 0.323 (−5.4);
ε 2·10⁻³, δ 3·10⁻³ 0.281 / 0.300 / 0.324 (−5.8); ε 2.5·10⁻³,
δ 3·10⁻³ at c 0.04 0.290 / 0.311 / 0.331 (−7.4) and at c 0.06 (the
defaults) 0.285 / 0.304 / 0.313 (−2.5); the lowest layer's air 0.267 /
0.289 / 0.291 (+4.9).

The defaults (`mfv64`, the same run): albedo 0.285, 0.299, 0.304,
0.304, 0.312, 0.314, 0.310, 0.310, 0.310, 0.313 on days 184–193
(package 3 with the Arctic cover, `arc10d64`: 0.287 … 0.338; the log's
M21 baseline `base10d64`: 0.288 … 0.337); ASR − OLR +3.6, −0.7, −1.5,
−0.1, −1.8, −2.4, −1.4, −1.0, −1.0, −2.5 W/m² (mean over days 188–193
−1.7; `base10d64` −11.5); global rain 1.76, 2.19, 2.23, 2.37, 2.41,
2.45, 2.59, 2.52, 2.51, 2.48 mm/d. Day 186 of the same code (`mfv3`)
against package 3 (`j0`): humidity below σ 0.9 in 30S–30N (the mean of
q/q_s by σ thickness and area) 0.712 against 0.737 (eight64_day0183
0.682; the old scheme's day 186 0.695); cloud water in the layers below
σ 0.9 20.7 against 30.6 g/m² globally, 13.4 against 37.1 in 30S–30N;
SE Pacific rain 0.36 against 0.56 mm/d in the day's mean, convective
rain on 0.09 against 0.51 of its column-days. `scripts/verticalAudit.mjs`
on day 193 against `arc10d64`'s day 193:

| day 193, N=64 | cumulus mass flux | package 3 + Arctic cover |
|---|---|---|
| humidity below σ 0.9, 30S–30N | 0.737 | 0.757 |
| cloud water below σ 0.9, g/m² (30S–30N) | 23.7 (16.8) | 32.1 (42.7) |
| SE Pacific rain, mm/d (convective share) | 0.80 (0.96) | 1.48 (0.96) |
| its columns firing a step (convecting at all) | 0.018 (0.605) | 0.142 (0.229) |
| its low cloud; deck runs | 0.137; 0.014 | 0.390; 0.219 |
| its EIS; deck's virtual jump, K | 3.69; 0.76 | 2.07; 1.52 |
| its resolved inversion, m (θv jump, K) | 1614 (5.11) | 1410 (3.85) |
| Pacific ITCZ rain, mm/d; ω500, Pa/s | 4.75; −0.032 | 6.96; −0.032 |
| ITCZ heating peak, hPa (K/d); lowest 100 m, K/d | 438 (20.0); −1.35 | 439 (24.9); −1.39 |
| ITCZ large-scale heating below 1 km, largest \|K/d\| | 1.28 | 20.03 |
| global rain, mm/d (convective share) | 2.37 (0.64) | 2.59 (0.75) |
| zonal-mean rain peak, mm/d (latitude) | 4.95 (0.5S) | 5.63 (1.5S) |

Five days at N=128 on the GPU from nine128_day0183 (`mfv128`): albedo
0.286, 0.278, 0.274, 0.275, 0.269 (`arc5d128` 0.296, 0.296, 0.299,
0.296, 0.304); ASR − OLR +5.9, +8.6, +10.3, +9.4, +10.8 W/m² (mean
+9.0; `arc5d128` +0.5); global rain 2.40–2.66 mm/d; cloud water below
σ 0.9 in 30S–30N 39.6 → 10.6 g/m² (`arc5d128` 39.9 on day 188), globally
32.7 → 19.5; the northern ice extent 0.072 → 0.093 Mkm² (`arc5d128`
0.090); the audit of day 188: SE Pacific rain 0.81 mm/d, firing 0.016
a step, low cloud 0.234, deck runs 0.159; Pacific ITCZ 7.57 mm/d; global
2.55 mm/d. Three days with c 0.04 end at +8.4 W/m² and with overshoot 0.5
at +6.9 (both −7.4 and −7.0 on day 193 at N=64). Pace: 0.95 min a model
day over days 185–188 (57 s, the log's 0.1 min resolution; two other
N=128 spin-ups shared the GPU at the start), the b4cc733 code's day 185
run next 1.0 min (60 s). `js/gpu/profile.module.js` over 128 steps
from nine128_day0183 after 64, alternating with the b4cc733 code twice
with nothing else on the GPU: a step's median 85.8 and 85.7 ms against
83.9 and 83.6 (+2.5 %, 1.1 s a model day), the adjust kernel's pass
16.2 against 14.7 ms and the physics kernel's 7.55 against 7.20, of
which the second saturation adjustment is 0.3 ms; three-day runs back
to back took 62.5 s a model day on days 185–186 against 59.8 for
b4cc733, and two-day runs earlier the same evening 77.5 and 74.9 s on
day 185 against 78.3 and 66.1. Three days at N=64 from nine64_day0091 lose
0.178·10³ km³ of northern ice a day on the GPU (9.191 → 8.657) and
0.177 on the CPU (9.19 → 8.66; package 3 with the Arctic cover 0.173,
9.19 → 8.67), the 70–90N ice's surface taking 67.1, 74.0, 70.4 W/m² net
(64.9, 71.7, 69.3).

**Item 7, one mass-flux scheme for all convection (Sept 30).** Both
engines (`plumeColumn` and `transportMomentum` in
`js/physics/moist.module.js`, the adjust and mixMomentum kernels of
`js/gpu/physics.gpu.js`) take all convection with one bulk plume
(`convection` 'plume', the default until item 8; 'bettsMiller' keeps the Betts–Miller
relaxation, its trigger and its 2 h activity, bit for bit, for side-by-side
runs; the activity is carried and saved but unused under 'plume'). The
plume leaves the shallow plume's source with its mean s_l and q_t
(`plumeSource` 'mean'; 'lowest'), rises unmixed to the LCL and from there
integrates d(w²)/dz = 2aB − 2bεw² exactly across each layer for constant
B and ε, from `plumeVelocity` 1 m/s, a = `plumeAcceleration` 1/3,
b = `plumeDrag` 1, B the virtual buoyancy with condensate loading;
ε = max(`plumeEntrainmentFloor` 10⁻⁴ m⁻¹, c_ε B/w²), c_ε =
`plumeEntrainment` 0.1, from the layer below's B and w² at the layer's
base. The net acceleration is 2B(a − b c_ε): with b = 2 and c_ε = 0.5
(the values first proposed) no plume accelerates, so b = 1 and c_ε = 0.1.
Condensate above `plumeRainThreshold` 0 rains at 1 − exp(−c0 Δz) per layer,
c0 = `plumeRainRate` 3·10⁻³ m⁻¹. The plume ends where w² reaches zero; one
whose top interface lies above σ 0.7 is deep, any other goes to the
shallow cumulus mass flux unchanged. Deep mass flux per unit base flux:
growing by exp((ε − δ)Δz), δ = max(0, ε − `plumeMassGrowth` 0), up to the
height of neutral buoyancy interpolated linearly in B between layer
midpoints, then falling linearly to zero at the top. The downdraft starts
at the layer of least moist static energy between the top and the
cloud-base layer, `downdraftShare` α = 0.3 of the base flux, saturated by
evaporating rain, entraining at `downdraftEntrainment` 10⁻⁴ m⁻¹ to the
cloud-base layer and detraining below it in proportion to mass, α lowered
where the rain made above a level does not cover what it evaporates to
there. s_l and q_t move in the shallow plume's flux form, both drafts'
fluxes at each interface, the rain made and evaporated as layer sources.
Closure: base flux (CAPE − `plumeCape` 70 J/kg) / (τ F), τ =
`plumeRelaxation` 1 h, CAPE the positive work of the plume's cloudy layers
(`plumeCapeParcel` 'plume'; 'undilute': of the source air), F the change
of the plume's net work over its layers per second and unit base flux from
the scheme's own tendencies with the plume held fixed, times the deck's
opening and the inhibition ramp about `inhibitionThreshold`; at most
`cumulusBoundaryLoss` of the source a step and (M + α|M_d|) gΔt/Δp ≤ 1.
`plumeClosure` 'separate' (the default) runs the shallow plume beside the
deep one on its own closure, 'cape' leaves it out of deep columns,
'maximum' gives the deep plume max(shallow closure, CAPE closure) as first
specified. The rain left after the downdraft falls from where it formed;
below the cloud-base layer it evaporates into cloud-free subsaturated
layers at 1 − exp(−`plumeRainEvaporation` (1 − q/q_s) Δz) of what falls,
10⁻³ m⁻¹, within `rainEvaporation` of the deficit; what reaches the ground
is the cell's convective rain. The cumulus cover below the shallow top is
the larger of the two plumes' M/(ρ w_u), w_u ≥ w0. `plumeMomentum` (off)
moves each edge's normal velocity by the mean of its cells' deep updraft
and downdraft in the same flux form, its kinetic-energy loss returned as
heat. Column enthalpy, water and (with momentum) each edge's momentum
close to rounding.

Tests (`test/convection.test.mjs`): the Jordan column lifts a deep plume
to 211 hPa (CAPE 273 J/kg, base flux 0.049 kg/m²/s, downdraft 0.015 from
518 hPa) that rains 34.7 mm/d after 1.7 mm/d evaporate below cloud base,
heats most at 440 hPa above the cloud-base layer (22.6 K/d; 30.7 K/d in
the cloud-base layer at 949 hPa), cools every subcloud layer (−23.4 K/d
over the lowest 100 m) and keeps enthalpy to 5·10⁻¹⁶ and water to 1·10⁻¹⁶;
the trade-wind column gets the shallow plume bit for bit and no rain; the
Jordan humidity above 850 hPa × 1, 0.8, 0.6, 0.4 tops the plume at 211,
292, 344 hPa and leaves it shallow (CAPE 273, 162, 89, 65 J/kg); the base
flux is (CAPE − CAPE0)/(τF) to 10⁻¹²; on 362 random columns the engines
agree under each closure, source and CAPE parcel to 1.8·10⁻⁴ K and
3·10⁻⁷ kg/kg, the base flux to 5.2·10⁻⁵ of 0.133 kg/m²/s, and with momentum
the winds to 1.6·10⁻⁴ m/s on 245 of 1080 edges (1.1·10⁻⁴ with
`downdraftShare` 0, no downdraft), each edge's column momentum to
2.2·10⁻¹¹ of layer momenta up to 2·10⁴; 'bettsMiller' reproduces the
pinned digests. `test/gpuModel.test.mjs` runs the rain accumulation under
both schemes (24 steps at N=6: plume per-cell rms 8.0·10⁻⁵ convective,
3.3·10⁻⁴ large-scale, last step's rain to 1.2·10⁻⁴ kg/m²).

`scripts/verticalAudit.mjs` prints the plume tops over 15S–15N by 100 hPa,
the deep and shallow shares of the column-steps, and ω500 at 5S–5N by 20°
of longitude. On the N=128 states after 8 CPU steps (plume / bettsMiller):
nine128_day0183 global rain 0.52 / 2.13 mm/d, Pacific ITCZ 0.01 / 4.33,
deep plumes on 0.000 of the tropical column-steps, shallow 0.520;
eight128_day0183 0.41 / 0.59, ITCZ 0.03 / 0.04 (Betts–Miller's activity
starts undecided there). The dilute plume finds no CAPE above 70 J/kg in
an atmosphere the Betts–Miller has held near its 5·10⁻⁵ m⁻¹ parcel, so
deep convection starts only as the troposphere destabilises: global rain
0.70, 1.14, 1.71, 2.02 mm/d on days 184–187 at N=64.

Ten-day N=64 GPU runs from eight64_day0183, the day-193 audit (8 steps;
ITCZ rain, its ω500, the firing columns' heating peak and lowest 100 m,
global rain, zonal-mean peak, SE Pacific rain and firing); runs marked †
on the build before the neutral-height interpolation and the continuous
CAPE consumption; all with c0 2·10⁻³ unless given but the defaults
(3·10⁻³):

| run | albedo 184 / 193 | ASR − OLR 188–193 | ITCZ mm/d; ω500 | peak hPa; lowest 100 m K/d | global | zonal peak | SE Pacific; firing |
|---|---|---|---|---|---|---|---|
| bettsMiller (= `mfv64`) | 0.285 / 0.313 | −1.7 | 4.75; −0.032 | 438; −1.4 | 2.37 | 4.95 (0.5S) | 0.80; 0.018 |
| † as first built | 0.304 / 0.340 | −3.0 | 5.58; −0.038 | 438; −11.2 | 2.26 | 6.20 (7.5N) | 1.66; 0.063 |
| † plumeSource 'lowest' | 0.301 / 0.345 | −1.9 | 5.43; −0.035 | 438; −12.9 | 2.19 | 5.07 (1.5N) | 2.13; 0.050 |
| † plumeClosure 'maximum' | 0.296 / 0.331 | +0.9 | 5.97; −0.038 | 944; −23.6 | 2.17 | 5.20 (2.5S) | 1.39; 0.074 |
| † plumeMassGrowth 2·10⁻⁴ | 0.303 / 0.355 | −5.0 | 7.52; −0.050 | 438; −4.9 | 2.23 | 5.60 (7.5N) | 1.42; 0.085 |
| † plumeMassGrowth −1·10⁻⁴ | 0.303 / 0.337 | −3.8 | 8.60; −0.057 | 944; −19.5 | 2.41 | 7.57 (11.5N) | 1.28; 0.032 |
| † plumeCape 120 | 0.304 / 0.345 | −4.2 | 6.25; −0.039 | 516; −5.8 | 2.30 | 6.17 (10.5N) | 1.97; 0.025 |
| † plumeCape 30 | 0.301 / 0.349 | −1.6 | 5.67; −0.035 | 438; −15.1 | 2.25 | 5.28 (7.5N) | 1.65; 0.040 |
| † c_ε 0.05 | 0.302 / 0.345 | −2.6 | 6.10; −0.041 | 438; −14.4 | 2.14 | 5.12 (4.5S) | 0.93; 0.025 |
| † c_ε 0.2, floor 2·10⁻⁴ | 0.305 / 0.336 | −6.1 | 4.77; −0.019 | 962; −4.6 | 2.41 | 7.49 (6.5N) | 1.00; 0.017 |
| † τ 2 h | 0.304 / 0.343 | −3.7 | 7.12; −0.053 | 438; −6.7 | 2.28 | 5.94 (11.5N) | 1.10; 0.043 |
| † α 0.5 | 0.304 / 0.342 | −3.7 | 6.95; −0.049 | 516; −7.9 | 2.37 | 7.10 (6.5N) | 1.34; 0.039 |
| † no downdraft, no evaporation | 0.302 / 0.345 | −3.1 | 4.90; −0.006 | 787; +0.4 | 2.27 | 6.22 (0.5N) | 1.40; 0.043 |
| † plumeClosure 'cape' | 0.303 / 0.344 | −3.7 | 5.61; −0.034 | 516; −10.4 | 2.37 | 6.06 (27.5S) | 2.04; 0.029 |
| as first built | 0.303 / 0.341 | −3.2 | 5.40; −0.037 | 438; −10.1 | 2.34 | 7.50 (8.5N) | 1.61; 0.046 |
| c0 4·10⁻³, no downdraft entrainment | 0.303 / 0.338 | −3.0 | 6.18; −0.031 | 438; −12.8 | 2.33 | 7.71 (3.5S) | 1.50; 0.075 |
| c0 5·10⁻³, α 0.2 | 0.303 / 0.338 | −2.9 | 6.83; −0.034 | 438; −10.7 | 2.37 | 6.84 (11.5N) | 2.14; 0.057 |
| w0 2 m/s | 0.303 / 0.348 | −3.0 | 6.20; −0.039 | 438; −9.9 | 2.18 | 4.83 (8.5N) | 2.27; 0.072 |
| floor 5·10⁻⁵ | 0.303 / 0.340 | −2.9 | 6.75; −0.043 | 438; −12.7 | 2.33 | 6.78 (3.5S) | 2.12; 0.060 |
| 'undilute' CAPE | 0.301 / 0.343 | −1.3 | 5.91; −0.041 | 921; −59.1 | 2.25 | 4.84 (8.5N) | 0.48; 0.016 |
| 'undilute', τ 2 h, α 0.15 | 0.301 / 0.346 | −2.3 | 6.93; −0.046 | 943; −33.8 | 2.33 | 5.86 (1.5N) | 1.59; 0.018 |
| shallow plume raining above 10⁻³ kg/kg | 0.303 / 0.338 | −3.0 | 5.85; −0.040 | 438; −9.9 | 2.41 | 8.11 (8.5N) | 1.65; 0.055 |
| **the defaults (`fin64`)** | 0.303 / 0.340 | −2.8 | 7.34; −0.057 | 438; −10.1 | 2.31 | 7.99 (7.5N) | 1.63; 0.036 |

The defaults (`fin64`): planetary albedo 0.303, 0.313, 0.317, 0.325,
0.340, 0.347, 0.339, 0.339, 0.337, 0.340 on days 184–193 (`mfv64` 0.285 …
0.313; `base10d64` 0.288 … 0.337); ASR − OLR +1.0, −0.1, −0.2, 0.0, −3.2,
−4.9, −2.7, −1.9, −1.4, −2.6 W/m² (mean over days 188–193 −2.8;
`mfv64` −1.7, `base10d64` −11.5); global rain 0.70, 1.14, 1.71, 2.02,
2.20, 2.29, 2.38, 2.35, 2.22, 2.25 mm/d. The day-193 audit: Pacific ITCZ
7.34 mm/d (convective share 0.81), ω500 −0.057 Pa/s, firing columns'
heating peak 438 hPa (8.9 K/d), −10.1 K/d over the lowest 100 m,
large-scale heating below 1 km at most 1.7 K/d (`mfv64` 1.28); global
2.31 mm/d (convective 0.42; 15S–15N 0.72); zonal-mean peak 7.99 mm/d at
7.5N; SE Pacific 1.63 mm/d (convective 0.26), firing 0.036 of the
column-steps, low cloud 0.108, deck runs 0.005, EIS 1.59 K, resolved
inversion 1955 m (θv jump 3.82 K) (`mfv64`: 0.80 (0.96), 0.018, 0.137,
0.014, 3.69 K, 1614 m); Hadley −180 / 59·10⁹ kg/s; deep plumes on 0.205
of the tropical column-steps (base flux 0.043 kg/m²/s), shallow on 0.449,
tops 200–300 hPa on 0.125 and 800–1000 hPa on 0.405; equatorial ω500
from 180W by 20°: +0.004, −0.013, +0.031, +0.013, +0.030, −0.021, +0.022,
+0.021, −0.016, −0.001, −0.003, −0.124, −0.177, −0.046, −0.017, −0.038,
−0.039, −0.019 Pa/s. The day's own means from the run's log: SE Pacific
0.48 mm/d (convective 0.16), Pacific ITCZ 3.27 mm/d. Day 193 against
`mfv64`, 15S–15N: relative humidity 0.72–0.79 at σ 0.84–0.88 against
0.55–0.63 and 0.72–0.74 at σ 0.31–0.37 against 0.57–0.65; temperature
1.0–1.8 K lower from σ 0.2 to 0.7; 10S–10N cloud water at σ 0.4–0.68
22.3 against 1.6 g/m², above σ 0.4 33.4 against 16.9; the frozen state's
albedo over four sun positions 0.322 against 0.308 (10S–10N 0.285 against
0.273, 10–30N 0.228 against 0.211, 30–90N 0.379 against 0.349), its cloud
above σ 0.68 alone 0.242 against 0.204 and below it 0.224 against 0.237;
convective rain 0.94 against 1.55 mm/d, large-scale 1.31 against 0.92
(40–60N convective 0.19 against 0.63).

Five days at N=128 on the GPU from eight128_day0183 (`fin128e`; the
shallow-only branch from the same state, `mf128e`, in brackets): planetary
albedo 0.267, 0.277, 0.293, 0.298, 0.302 (0.257, 0.275, 0.298, 0.301,
0.300); ASR − OLR +10.6, +8.9, +5.7, +4.8, +5.0 W/m², mean +7.0 (+11.5,
+6.5, +0.3, −0.2, +0.4, mean +3.7); OLR 239.1 → 232.8 W/m² (241.6 →
238.0); global rain 0.67, 1.32, 1.92, 2.33, 2.50 mm/d (1.73 … 2.66). The
day-188 audit (`mf128e` audited under its own scheme): Pacific ITCZ 3.56
mm/d (5.33), ω500 −0.009 Pa/s (−0.018), firing columns' heating peak 515
hPa (438), −9.9 K/d over the lowest 100 m (−1.8); global 2.59 mm/d (2.69),
convective share 0.34 (0.62); zonal-mean peak 7.72 mm/d at 10.5N (5.79
at 9.5N); SE Pacific 0.05 mm/d (0.03), firing 0.003 (0.000), low cloud
0.157 (0.157), deck runs 0.003 (0.004), resolved inversion 1704 m (1691);
Hadley −128 / 58·10⁹ kg/s (−113 / 49); deep plumes on 0.209 of the
tropical column-steps (base flux 0.042), tops 200–300 hPa on 0.141;
equatorial ω500 from 180W by 20°: +0.006, +0.046, +0.021, +0.011, +0.005,
+0.003, +0.049, +0.036, +0.033, −0.025, +0.031, +0.030, −0.020, −0.021,
−0.052, −0.134, −0.089, −0.043 Pa/s (ascent 100E–180, descent
180–80W). At day 188 against `mf128e`, 15S–15N: relative humidity
0.69–0.79 at σ 0.84–0.88 (0.55–0.64), 0.68–0.69 at σ 0.27–0.37
(0.58–0.60), temperature 1.5–2.3 K lower from σ 0.23 to 0.6; convective
rain 0.93 mm/d (1.69), large-scale 1.57 (0.97). From nine128_day0183
(`fin128n`, an M21 state out of balance; `mfv128` in brackets): albedo
0.282, 0.279, 0.287, 0.295, 0.300 (0.286, 0.278, 0.274, 0.275, 0.269);
ASR − OLR +9.5, +11.8, +11.2, +9.5, +8.5, mean +10.1 (+9.0); global rain
1.00, 1.88, 2.24, 2.42, 2.41 mm/d; day 188: Pacific ITCZ 8.36 mm/d
(convective 0.77), ω500 −0.056 Pa/s, heating peak 438 hPa, −12.2 K/d
over the lowest 100 m; global 2.38 mm/d (convective 0.41); zonal-mean
peak 7.46 mm/d at 7.5N; SE Pacific 1.15 mm/d (convective 0.03), firing
0.011, low cloud 0.292, deck runs 0.108, resolved inversion 1554 m;
equatorial ω500 −0.094 at 180W, +0.028 to +0.062 over 160W–120W, −0.071
and −0.115 over 140E–180.

Three days at N=64 on the GPU from nine64_day0091 (`fin64i`): 60–90N ice
9.191 → 8.659·10³ km³, 0.177·10³ km³ a day (package 3 with the Arctic
cover on the CPU 0.173, the shallow branch on the GPU 0.178).

Pace: `js/gpu/profile.module.js` over 128 steps from nine128_day0183
after 64, the plume, b4cc733 and this code under 'bettsMiller' (the
shallow branch) alternated twice with nothing else on the GPU: a step's median 87.0 and 85.2 ms
against 83.8 and 83.8 for b4cc733 (+2.7 %) and 85.7 and 85.7 for the
shallow branch; the adjust kernel's pass 16.1 and 15.7 ms against 14.7
and 14.7 and 16.2 and 16.3, the physics kernel's 7.60 against 7.24 and
7.60. Two to eight spin-ups from another worktree shared the GPU during
the N=128 runs, which took 2.5 to 11.7 min a model day (1.7 min on days
187–188 of `fin128n`); the uncontended ten-day N=64 runs took 6.6 s a
model day (`fin64`) against 6.0 for the shallow branch (`pt`, the same
code under 'bettsMiller').

Review (Sept 30). Without a downdraft the deep plume's downdraft momentum
mixing factors are 1 in every layer; before that, `plumeMomentum` left NaN
winds on 204 edges after one CPU step from fin64_day0193. One CPU step of
`moist.adjust` over every column of fin64_day0193 with `plumeMomentum`:
20,199 plume columns (5,191 deep; the 2 where the filler lost water left
out) keep column enthalpy to 8.9·10⁻¹⁶ and water with the rain to
8.2·10⁻¹⁶; the transport moves 17,943 edges, each edge's column momentum to
2.9·10⁻¹⁰ kg/m/s against layer momenta up to 4.6·10⁴, the heat matching
the kinetic energy lost to 1.2·10⁻¹¹; 27 edges gain kinetic energy, which
no heat pays for. The worker-thread CPU step reproduces the single thread
bit for bit over two steps from fin64_day0193 with `plumeMomentum`. The
closure's F on the 4,063 deep columns of the second CPU step from
fin64_day0193 under `plumeClosure` 'cape': positive in every column, the
buoyant cloud layers (the CAPE's own) giving a median 0.88 of it (5–95 %
0.66–1.00), the negatively buoyant cloud layers a median 0.11 (95 % 0.33),
the layers between the cloud-base layer and the source at most 0.05 (95 %);
the source layers, which the downdraft cools, are outside F; the base flux
falls below (CAPE − CAPE0)/(τF) through the boundary-loss or Courant
limit on 32. No candidate with CAPE above 70 J/kg has F ≤ 0; 10,737 of the
14,802 plumes topping above σ 0.7 have CAPE at most 70 J/kg and go to the
shallow plume. Ten days at N=64 from eight64_day0183 rerun (`adv64`)
reproduce `fin64` byte for byte, its log and its day-193 audit;
`pt_day0193` is byte-identical to `mfv64_day0193`. Equatorial ω500 on day
193 by 20° from 120E to 160W: −0.038, −0.039, −0.019, +0.004, −0.013 Pa/s
(`mfv64` audited under 'bettsMiller': −0.067, −0.065, −0.030, −0.042,
−0.038), from 40E to 80E −0.124, −0.177 (+0.028, −0.022); SH Hadley −180
(−131), NH 59 (45) ·10⁹ kg/s. twin64_day0900 and five64_day2281 (27
layers), six64_day1004, seven64_day0639, eight64_day0183, nine64_day0091
and m21b64_day0183 load and take two CPU steps under the plume with no
non-finite value, 23 to 217 deep plumes and 0.14–0.52 mm/d of rain. Pace
with nothing else on the GPU: two days at N=128 from nine128_day0183
(`adv128`) at 1.0 and 2.0 min, saved after 2.2 min, against 0.9, 1.9 and
2.0 for b4cc733 (`adv128b`); `profileGpu` alternated twice, a step's median
87.4 and 85.1 ms against 84.1 and 84.0 for b4cc733 (+2.6 %) and 85.7 and
85.9 under 'bettsMiller', the adjust pass 15.90 and 15.53 ms against 14.52
and 14.72 and 16.24 and 16.27.
**Item 8, the integration (Sept 30).** The integrated code carries
items 6 and 7 and the boundary-layer entrainment of item 6. Runs from
copies of eight64_day0183 (ten days, N=64), eight128_day0183 (five days,
N=128) and nine64_day0091 (three days, N=64), GPU, `everySteps` 8.

The deep plume on the integrated code: N=64 albedo 0.296, 0.305, 0.310,
0.321, 0.332, 0.342, 0.337, 0.334, 0.335, 0.333 on days 184–193, ASR − OLR
over days 188–193 −0.8, −3.2, −1.7, −0.1, −0.7, −1.0 W/m² (mean −1.3);
N=128 albedo 0.262, 0.273, 0.287, 0.293, 0.301, ASR − OLR +12.2, +10.1,
+7.7, +6.7, +5.5 (mean +8.4), OLR 239.2 → 232.6. The Betts–Miller deep
branch with the shallow plume: N=64 0.274 … 0.296, mean +4.8; N=128 0.247,
0.259, 0.278, 0.283, 0.283, mean +8.7. Three-day N=64 screens, day-186
albedo and ASR − OLR (plume 0.310, +2.2): `criticalHumidity` 0.85 0.315,
+1.3; 0.9 0.324, −0.4; `plumeEntrainment` 0.15 0.314, +1.4; 0.2 0.317,
+0.4; `plumeCape` 30 0.304, +3.9; `cloudLifetime` 1 h 0.271, +4.8.
`upperCloudLifetime` (both engines; the lifetime of cloud water above
`shallowTop`) 1 h: day 187 0.286, +2.2. Ten days at N=64 with
`cloudLifetime` 2 h: albedo 0.284 … 0.302, ASR − OLR over days 188–193
−0.3, −2.4, −2.9, −1.2, +0.6, +2.8 (mean −0.6). At N=128 only the
plume's defaults were run; `convection` defaults to 'bettsMiller'.

Entrainment (`js/physics/boundaryLayer.module.js`, both engines):
w_e = o (1 − s) min(cap, (A B0 + A_s r u*³/h) / max(Δb, b_min)), o the
deck's opening (1 at a gate of 0.5, 0 at 0.6), s the radiation's
stratiform share (0 at EIS 8 K, 1 at 12 K; `radiation.stratiform`, PH
`STRAT`), r = min(1, B0 / `shearOnset`), `shearOnset` 5·10⁻⁵ m²/s³; A 0.2
and A_s 5 are the constants of Driedonks (1982) for the two sources of
Tennekes (1973). Tests: a gate of 0.55 or a share of 0.5 halves w_e to
10⁻¹², 0.6 or a share of 1 gives 0, both together a quarter; w_e is linear
in B0 below the onset (0.080, 0.040, 0.016, 0.002, 0.000 mm/s at B0 5·10⁻⁵,
2.5·10⁻⁵, 10⁻⁵, 10⁻⁶, 10⁻⁸ against 0.080, 0.053, 0.037, 0.027, 0.026
without it); 362 random columns, 78 tapered: engines agree on w_e to
5.8·10⁻⁴ relative, θ to 1.2·10⁻⁴ K, q to 3.3·10⁻⁸; one N=6 step under 6
and 9 K inversions: the share on the ramp in 35 and 50 columns, engines to
7.5·10⁻⁶, w_e to 2.1·10⁻³ mm/s. One CPU step from acc64_day0193, sea
30S–30N with B0 > 0 (13,206 columns), area means: w_e 3.26 mm/s, 3.43
without the onset, 4.16 untapered and 3.36 under the cut at a gate of 0.5
without the share or the onset; 1,960 columns below the onset, 630 with a
share, 499 on the gate's ramp. Three days at N=64: albedo and ASR − OLR on day 186 0.310, +2.2
(plume) and 0.284, +4.6 (Betts–Miller) against 0.310, +2.2 and 0.283, +5.1
before.

The deck's height under cumulus. On mfv64_day0193 (eight CPU steps, SE
Pacific, 4,480 column-steps): Richardson depth 811 m, carried height 818
m, start height 863 m (within 20 m of the floor on 0.620), ceiling 1861 m
on 1.000, jump at h 0.75 K (2 K passed on 0.057), resolved inversion 1919
m (7.14 K), mixed layer cloud-free on 0.535 (cloud base 783 m), gate 0.041,
deck runs 0.014; Peru: 765, 850, 864 m (0.628), ceiling 1417 m, 1.40 K
(0.110), 1271 m (5.95 K), cloud-free 0.390, deck runs 0.226. The binding
rule is the rest: a deck that does not run relaxes to the Richardson depth
and is tested there. `deckRest` 'inversion' (the default; 'depth' keeps the
old rest) starts an unset height and relaxes a resting one toward the
ceiling. Three days at N=64 under the Betts–Miller deep branch, day 186:
'depth' 0.284, +4.6, gate open over 0.061 of the globe; 'inversion' with
`minimumInversion` 2 K 0.326, −7.9, 0.278 (10S–10N 0.206); 3 K 0.325,
−7.9; 4 K 0.306, −2.3, 0.150 (10S–10N 0.069). `minimumInversion` defaults
to 4 K. Under the plume, 'inversion' with 2 K: SE Pacific low cloud 0.614,
deck runs 0.634, Peru 0.573, 0.755 on day 186. Tests: a column whose
Richardson depth is 114 m under a 1200 m inversion forms its deck after
17 h at 1200 m resting at the ceiling and none in 48 h resting at the
depth; the engines agree on the carried height to 1.2·10⁻⁷ (3.2·10⁻⁷ at
rest) and the gate to 4·10⁻⁸.

Acceptance, the defaults (Betts–Miller deep, shallow plume, the deck at
the inversion, 4 K), on the reviewed code (the share kept off land, below):

| | value | asked |
|---|---|---|
| N=64 albedo days 184–193 | 0.273, 0.299, 0.306, 0.304, 0.305, 0.309, 0.313, 0.318, 0.312, 0.312 | 0.30–0.32 from 186 |
| N=64 ASR − OLR days 188–193 (mean) | −0.2, −0.8, −1.8, −2.8, −1.1, −0.8 (−1.3) | ±4 |
| N=64 global rain days 187–193, mm/d | 2.42–2.67 | 2.4–2.8 |
| SE Pacific rain, mm/d; firing; low cloud; deck runs | 1.29; 0.038; 0.264; 0.296 | < 0.5; < 0.02; 0.4–0.7; ≥ 0.3 |
| SE Pacific EIS; deck height where it runs; resolved inversion | 4.39 K; 2045 m; 1908 m (5.01 K) | |
| Peru rain; firing; low cloud; deck runs | 0.00; 0.000; 0.424; 0.429 | |
| Pacific ITCZ rain, mm/d; ω500, Pa/s | 3.71; −0.015 | 6–9 |
| ITCZ firing columns' heating peak; lowest 100 m | 438 hPa; −1.07 K/d | 400–500; −10 to +5 |
| zonal-mean rain peak | 5.64 mm/d at 1.5N | 5–7 at 5–10N |
| ω700 grid-scale share | 0.235 | < 0.3 |
| N=128 albedo days 184–188 | 0.247, 0.265, 0.288, 0.305, 0.308 | 0.29–0.32 on 186–188 |
| N=128 ASR − OLR (mean) | +14.7, +9.8, +3.8, −1.7, −2.4 (+4.8) | ±4 |
| N=128 day 188: SE Pacific rain; low cloud; deck runs; ITCZ rain | 0.12; 0.287; 0.361; 4.31 | |
| N=128 equatorial ω500 100E–180 / 160W–80W, Pa/s | −0.027 to −0.062 / +0.026, +0.045, −0.012, +0.010, +0.031 | ascent / descent |
| 60–90N ice loss, 10³ km³/day | 0.183 (9.191 → 8.643) | ≤ 0.18 |
| fresh start, days 1–30 albedo | 0.276, 0.298, 0.353, 0.394, 0.420, 0.439 … 0.289 (day 17) … 0.281 | 0.29–0.33 by 30 |
| fresh start ASR − OLR, day 30 | +14.1 (+7.0 to +17.1 over days 18–30) | ±10 |

On 94f7e7b, before the review, the same runs gave: N=64 albedo 0.303 …
0.320, ASR − OLR mean −1.6; SE Pacific 1.15 mm/d, firing 0.031, low cloud
0.346, deck runs 0.356 (deck at 2067 m, resolved inversion 1937 m); Peru
low cloud 0.464, deck runs 0.467; ITCZ 3.96 mm/d; zonal peak 5.15 mm/d at
1.5S; N=128 mean +4.8; ice loss 0.183 (9.191 → 8.641); fresh start day 30
0.277 and +13.2. The SE Pacific and ITCZ numbers come from one eight-step
window of the day-193 state.

Ice loss on 94f7e7b with 'depth' and 2 K 0.179, under the plume 0.183.
The fresh start: no NaN, clamped 0 on every day; convection after 30
days: share 0.76 global, 0.99 15S–15N; SE Pacific 1.89 mm/d (convective
1.87) on 0.36 of its column-days; Pacific ITCZ 8.66 mm/d; equator after
30 days: surface current −0.23 m/s (160E–100W) and −0.77 m/s
(140W–100W), undercurrent +0.19 m/s at 93 m, stress −0.040 N/m², mixed
layer 52 m, 1024 class top 185 m and 51 m. Ten days at N=64: convection
after 10 days 0.67 global, 0.99 15S–15N, SE Pacific 0.83 mm/d (0.72) on
0.09, Pacific ITCZ 3.95 mm/d; equator −0.16 and −0.36 m/s, +0.06 m/s at
88 m, −0.019 N/m², 51 m, 144 and 90 m. Pace on 94f7e7b at N=128: days 187 and 188 0.9 and 1.0 min (two
N=64 runs sharing the GPU); `profileGpu` over 128 steps from
eight128_day0183 after 64, alternated twice with the merged parent
(84d7f40) alone on the GPU: a step's median 86.1 and 86.2 ms against 84.8
and 84.5 (+1.8 %), the adjust pass 15.97 and 15.99 ms against 14.58 and
14.59, the physics and boundary-layer passes 7.67 against 7.62.

Review (Sept 30). Rerun from a copy of eight64_day0183, 94f7e7b
reproduces the ten-day N=64 run byte for byte (state and audit). One CPU
step of every column of acc64_day0193: the boundary-layer mix keeps each
column's mass-weighted θ to 7.7·10⁻¹⁶ and water to 8.0·10⁻¹⁶ relative;
`moist.adjust` keeps c_p T + L q to 5.9·10⁻¹⁶ and water with the rain to
7.2·10⁻¹⁶ except in 59 columns, each of which gains no more water than the
negative water it entered with (the filler); the same on the reviewed
code's day 193 (61 columns), and with `upperCloudLifetime` 1 h or the
plume. States from twin64_day0810 (27 layers), seven64_day0365 and
five64_day2190 load and step on the CPU, twin64_day0810 a day on the GPU.

The share's taper acted over land: one CPU step from acc64_day0193, 2,047
of the 7,868 land columns with B0 > 0 tapered, removing 1.22 of their
6.73 mm/s area-mean w_e, and 33 columns over ice. Both engines now keep
`stratiform` at 0 over land and full ice, where no deck forms (land
parity test in `test/landGpu.test.mjs`: under a 9 K inversion the share
would be positive on 36 land cells; engines to 6.1·10⁻⁶, w_e to 3.7·10⁻³
mm/s). On the reviewed code's day 193 (rv3b_day0193), sea 30S–30N with
B0 > 0 (13,139 columns): w_e 3.22 mm/s against 4.12 untapered; the gate's
ramp (0.5 < G < 0.6) tapers 469 entraining columns by 0.08 mm/s of the
area mean, a gate of 0.6 or more stops 1,101 (0.44 mm/s), the share
tapers 409 (34 with a share of 1) by 0.11 mm/s, the onset 1,703 by 0.17
mm/s. Mean w_e by B0 with and without the onset: 0.02 and 2.02 mm/s for B0
up to 10⁻⁶ m²/s³ (29 columns), 0.17 and 2.60 to 5·10⁻⁶, 0.49 and 2.91 to
10⁻⁵, 1.02 and 2.75 to 2.5·10⁻⁵, 1.96 and 2.48 to 5·10⁻⁵; 2,121 columns
have B0 ≤ 0.

The deck against package 2's record (eight CPU steps from each day-193
state; running column-steps): with the inversion rest the SE Pacific deck
is within 20 m of the floor on 0.014 (acc64) and 0.000 (rv3b) of its
running steps, against 0.44–0.54 in package 2, and within 20 m of the
ceiling on 0.186 and 0.175, against 0.03–0.06. It mixes the layer above
the resolved inversion (the interface of largest dθv/dz) on 0.118 and
0.172 of them, every one under a jump below 4 K (mean 2.92 and 2.66 K);
over all open sea on 0.141 and 0.136 (1.98 and 1.90 K). Package 2 and
t_bm_day0193 under the old rule ('depth', 2 K): never in the SE Pacific,
0.049 over the sea (1.06 K). The 4 K default sets the ceiling as well as
the regime test, so the ceiling is the lowest 4 K jump, not package 2's
2 K. `ceilingInversion` (null: `minimumInversion`) separates them; ten
days at N=64 with a 2 K ceiling (rvc2b): albedo 0.269, 0.288, 0.295,
0.296, 0.300, 0.300, 0.303, 0.298, 0.297, 0.298, ASR − OLR over days
188–193 +1.7, +1.9, +1.2, +2.8, +3.8, +3.8 (mean +2.5); day 193: the deck
no longer mixes past the resolved inversion, but the SE Pacific deck runs
on 0.162 (rests within 20 m of the ceiling on 0.447, where the 4 K test
passes on 0.064), low cloud 0.229, rain 1.03 mm/d, Peru deck runs 0.378.
The default keeps the 4 K ceiling. A resting height relaxes toward the
ceiling where the column has one and toward the Richardson depth where
it has none (0.53 of the SE Pacific column-steps have one on rv3b).

**Betts–Miller retired (Oct 1).** The plume of item 7 is the only
convection scheme on both engines. Removed: the deep relaxation, its
parcel, reference profile and trigger, the shallow vent and mixing
line, the `convection`, `shallowScheme` and `cumulusWithDeep` switches,
and the per-cell activity (state field `convectiveActivity`, PH
`CONVACT`). The retired options (`RETIRED_OPTIONS` in
`js/physics/moist.module.js`) throw on either engine; states saved with
the activity load without it. The plume path is unchanged bit for bit:
24 CPU steps at N=4 and 48 GPU steps at N=6 give the parent's digests
under `convection` 'plume', with and without `plumeMomentum`. The
defaults are those of the plume runs above (N=64 albedo 0.333 on day
193, ASR − OLR over days 188–193 −1.3; N=128 albedo 0.301 on day 188,
mean +8.4), which ran before the review's land fix. Digest tests
re-pinned to 3ca002d1, 3d0c610f and da3ea94c. Parity on 362 random
columns under seven option sets: θ to 7.6·10⁻⁵ K, q to 3.6·10⁻⁸, q_c to
1.4·10⁻⁸, base mass flux to 6.5·10⁻⁶ kg/m²/s, cumulus fraction to
6.1·10⁻⁶, no column differing in whether it rains, whether a plume rises
or where it tops.
The radiation-parity fixture is one step after the 10 K inversion: its
layer-24 heating differs between the engines by 7.9·10⁻⁵, 1.4·10⁻⁴,
1.6·10⁻⁴ and 1.5·10⁻⁴ K/day after one to four steps under the plume,
against a tolerance of 10⁻⁴ (9.8·10⁻⁵ after two steps under
Betts–Miller).

Review (Oct 1). From an archive of 7186bdd, 24 CPU steps at N=4 and 48
GPU steps at N=6 under `convection` 'plume' give the same state digest as
the defaults here, with and without `plumeMomentum`. One CPU step of every
column of rv3b_day0193: the boundary-layer mix keeps each column's
mass-weighted θ to 8.4·10⁻¹⁶, q to 8.4·10⁻¹⁶ and q_t to 8.5·10⁻¹⁶
relative, and each edge's column momentum to 1.2·10⁻¹⁵ of Σm|u|;
θ_l = θ − L q_c/(c_p Π), which it does not mix as one field, moves by up to
6.3·10⁻⁷ relative. `moist.adjust` keeps c_p T + L q to 9.0·10⁻¹⁶ and water with
the rain to 8.1·10⁻¹⁶ in every column. nine64_day0091, saved with the
activity: the same to 9.1·10⁻¹⁶ and 8.4·10⁻¹⁶ but for two columns that
gain the negative cloud water the filler removes (5.1·10⁻⁶ and
1.6·10⁻⁶ kg/m²); it loads and takes four GPU steps. The heating-parity
layer is σ 0.93 under a deck of cover 0.41–0.67, cooling 2.7–3.7 K/day.
The full suite (44 files, run concurrently) passes with nothing skipped.

**The cumulus cloud's memory (Oct 3).** The plumes' radiative cloud
(`moist.cumulusCover` and `cumulusWater`, PH `CUCOVER` and `CUWATER`, the
plume layers from `cumulusK0` down) was diagnosed afresh each step from
the mass flux, f = M/(ρ w_u) with the plume's condensate at the layer's
midpoint, and was zero on any step the plume did not run; the page's cloud
column counts it, and a plume firing on alternate steps made a cloud that
blinked (95 % of the overlay's blinks at main, M22). `cumulusMemory` τ
(seconds, default 1800; 0 the instantaneous cloud) carries it: after the
plume stage of every step, in `adjust` on the CPU and around the adjust
kernel's plume on the GPU, each layer's cover f and path P = f × water
relax toward the step's diagnosed f' and P' (0 where no plume ran) as
X ← X' + (X − X') e^(−Δt/τ), the water is P/f, and both are 0 below
CUMULUS_TRACE (10⁻⁶ of cover). The deep plume's cloud, merged into the
same fields, takes the same memory. The filter is linear in f and in P, so
the time means of the cover and of the cover × water the radiation sees
are kept; filtering the water itself would give a 0/x plume a quarter of
its path. The cloud carries no water mass (the plume's condensate rains
or detrains as before), so there is no budget. The scheme is Tiedtke's
(1993), kept by the IFS, in which convective detrainment is the source of
a cloud that then decays on its own timescale, reduced to a first-order
decay of the diagnosed cloud. The value: in a tracked LES shallow-cumulus
ensemble (RICO, 25 m LES over 50 km) active clouds live about 20 min on
average and passive ones about 5 min, 3–7 min over all clouds (Sakradzija,
Seifert and Heus 2015, Nonlin. Processes Geophys. 22, 65–85, Table 2); τ = 30 min lies past the active clouds' mean life by the decay of
what they leave. The persistence of a grid cell's cumulus field beyond a
cloud's life is the plume's forcing persisting, which the closure already
reads each step, so the memory is not set to the field's decorrelation
time. The fields are saved as before, so a state saved before this
change loads with its last instantaneous cloud as the memory's start
(asyncSpinup.test.mjs: one GPU step at N=6 from a day saved with τ = 0
lies within 2.4·10⁻⁷ of X' + (X − X') e^(−Δt/τ) from the saved X in
every plume layer, 710 relaxing toward the plume's and 64 keeping the
saved cloud where the plume made none; the runs split at and inside a
day still end byte for byte on the uninterrupted run's files). With
τ = 0 both engines hash as at f319996 (GPU: state, the whole PH buffer
and the frame after 16 steps from eleven64_day1825 and eleven128_day1825;
CPU: 24 steps at N=6). Over one GPU day at N=64 from eleven64_day1825
the time mean of the cumulus path is 5.513 g/m² at f319996 and 5.507
with the memory, the share of cell-steps with a plume 0.600 and 0.599;
the mean of each column's largest layer cover falls from 0.0071 to
0.0066, the largest of smoothed fields. Over three N=64 days and two
N=128 days with the memory alone the day means move by at most 0.2 W/m²
in SWCRE (−62.7, −61.4, −60.9 → −62.8, −61.5, −61.0 at N=64; −53.3,
−53.6 → −53.5, −53.8 at N=128), 0.1 in LWCRE, ASR and OLR and 0.01 mm/d
in rain, against 0.1 W/m² between replicates; the blinking it removes is
in M22.

Three tests moved with the memory, none through an engine difference.
The pinned twelve-step digest of the moist defaults (physics.test.mjs)
is the parent's with cumulusMemory 0 and has its own with the memory.
The rain split at the diagnostics (gpuModel.test.mjs, 24 steps at N=6
under convectionType 'top') counted 8 cells apart against its bound of
0.02 C (7.2): the CPU against itself under ±10⁻⁴ K of θ noise each step
(16 seeds) parts 2–11 cells by the test's measure with the memory and
1–6 without, every parting a resolved cloud's or the dry adjustment's
decision or beside one, the cumulus traces alike; the test now takes
the shared rule, leaving out the columns within two cells of a parted
decision of test/helpers/decisions.mjs (CPU: 2–15 decisions, 35–190
columns with their neighbours, 0–4 apart outside them; GPU 25, 242 and
1), with the shared bounds C/6 and 0.8 C and the old 0.02 C on the cells
apart outside. The soil-carbon test (soilCarbon.test.mjs, 24 steps at
N=6 with the carbon accelerated 3·10⁵ times) parted by 3.0·10⁻² kg/m²
against 1.8·10⁻²: in cell 18 the shallow plume fires at step 3 on one
run only and its rain wets the soil; the CPU against itself parts the
same cell by 3.0·10⁻² in 8 of 16 seeds with the memory and in none of
8 without. It now leaves out the land cells where the plume fired on one
engine only at some step (CPU: at most 4 of 151, the rest within
3.0·10⁻³; GPU 1, the rest within 5.4·10⁻⁴), bound 0.05 of the land
cells, the carbon's 10⁻³ of its largest change unchanged.

### M22 — A moist boundary layer — done (first tuning; acceptance partly met)

The boundary-layer scheme was the dry Troen–Mahrt K-profile of M14
with the explicit top entrainment of M21: it mixed temperature and
vapour separately, felt no cloud-top cooling, and handed the
stratocumulus regime to the mixed-layer deck model through a gate.
The milestone asked for a moist turbulence closure in conserved
variables that takes cloud-top longwave cooling as a source of
turbulence, entrains at the inversion from its own closure, treats the
stratocumulus-topped layer as one of its regimes, gives the layers it
mixes a cover from its own variance, and takes the equatorial
boundary-layer momentum budget. Acceptance: low cloud of 0.5–0.7 over
the SE Pacific, Peru, Namibia and California boxes with the inversion
at 1–1.5 km, the equatorial stress toward Earth's 0.04–0.06 N/m², the
ten-day and thirty-day tests of M21 balanced, and the pace.

Built (`js/physics/boundaryLayer.module.js`, `pblDiagnose` and the
adjust kernel of `js/gpu/physics.gpu.js`; `turbulence` 'moist', the
default, 'dry' the M14/M21 scheme, which with the parent's moist
defaults, `cloudLifetime` 3 h and `plumeCape` 70, reproduces the parent
80b6bb9 byte for byte over three GPU days at N=64 from
eight64_day0183; under the defaults below it does not):

- Conserved variables. The implicit solve mixes θ_l = θ − L q_c/(c_p Π),
  q_t = q + q_c and the edges' momentum with one set of interface
  coefficients; each layer it touches leaves with θ = θ_l, q = q_t and
  no cloud water, the saturation adjustment that follows returns the
  cloud, and untouched layers keep their values exactly.
- Two profiles. The surface-driven K-profile of M14 over h_s, the
  Richardson depth, or where B0 > 0 the top of a surface parcel (θ_l
  with Holtslag and Boville's excess 8.5 B0 θ_v/(g w_m), q_t, rising with
  its condensate in equilibrium, stopping where its θ_v falls short of
  the layer's by more than 0.5 K) when it has condensed, stops within
  400 m of the base of its first saturated layer and below 3 km, and
  lies above the Richardson depth (Lock et al. 2000's parcel test for a
  stratocumulus-capped layer; a parcel that rises further is cumulus and
  keeps the Richardson depth). The cloud-top profile, where the lowest
  run of cloudy layers (q_c > 10⁻⁶) tops out below 3 km and cools:
  ΔF the longwave cooling summed over the run's layers (the radiation
  keeps each layer's longwave heating, `longwave`, PH `LWH`),
  V³ = (g/θ_v) ΔF/(ρ c_p) z_ml, K = 0.85 κ V z_ml x² (1 − x)^½ with
  x = (z − z_b)/z_ml over z_b < z < h_c, h_c the cloud top's upper
  interface, z_b where a parcel of the cloud top's θ_l less 0.2 K and
  q_t stops sinking (0 when it reaches the lowest layer). The two K add.
- Entrainment across the interface above the mixed layer:
  w_e = min(5 cm/s, (A (w_s³ + V³) + 5 r u*³)/(h max(Δb, 0.015 m/s²))),
  w_s³ = B0 h for a surface-driven top, Δb from the θ_v of the layer
  above or the one above that, whichever is warmer (the inversion's own
  grid layer holds part of the jump), A Nicholls and Turton's
  0.2 [1 + 25 χ* (1 − Δθ_vs/Δθ_v)] at most 1 from the cloudy top layer
  and the jumps in θ_l and q_t as the deck model takes it, 0.2 under a
  clear top (the M21 form); untapered (`entrainment.taper` restores the
  M21 taper; `jumpLayers` 1 the one-layer jump). A decoupled column also
  entrains across its surface-driven top.
- Regimes, per cell (`regime`, PH `REGIME`, saved as `boundaryRegime`
  with the mixing top `mixingTop`): stable, surface-driven, decoupled
  (a cloud-top layer whose z_b lies above h_s: the plume runs into the
  cloud layer) and coupled (z_b at the surface or within h_s: the
  surface profile reaches h_c as well, which is then `depth`).
  `deckRegime` 'boundaryLayer' gates the mixed-layer deck by the coupled
  regime instead of the 4 K jump, and `deckBypass` leaves those columns
  to the resolved cloud; `coupledVeto` (moist) stops the plume in
  coupled columns.
- Cover (`boundaryCover` 'variance', the radiation, both engines): a
  cloudy layer below the mixing top covers ½ [1 + erf(Q₁/√2)],
  Q₁ = a_l (q_t − q_sl(T_l))/σ_s, σ_s = max(0.002 q_sl, 5 l a_l |∂q_t/∂z −
  Π q_s' ∂θ_l/∂z|), l = κz/(1 + κz/λ) with λ 300 m, or 30 m where the
  surface buoyancy flux is not positive, the gradients those to the
  neighbouring mixed layers, erf by Abramowitz and Stegun 7.1.26; on
  the EIS ramp of M21 it blends into the bounded cover as the PDF does.
  The layers above keep the uniform PDF.

Tests (`test/moistBoundaryLayer.test.mjs`, 100 m layers to 2 km): a
stratocumulus column over a 26 °C sea under an 8 K inversion at 1.3 km
with 60 W/m² of cloud-top cooling is coupled with its cloud top at
1243 m, V 1.30 m/s, the cloud-top profile adding K on all 11 interfaces
with its largest at 1030 m; A 0.463 (χ* 0.036, Δθ_v 7.16 K) and w_e
3.84 mm/s equal the Nicholls–Turton rate to 10⁻¹²; one 300 s step keeps
θ_l and q_t to −2.7·10⁻¹⁶ and 1.4·10⁻¹⁶ and the edge momentum to
10⁻¹⁵·K; after 6 h under the cooling θ_l spreads by 0.271 K and q_t by
0.106 g/kg from the surface to the cloud top, the cloud is 106 m thick
(top 1240 m) and its variance cover 1.000. A clear convective column
gives the dry scheme's coefficients to 4·10⁻¹⁶ and its θ after a step to
5.7·10⁻¹⁴ K (w_e 1.57 mm/s, 1.13 with the two-layer jump). A cloud layer
over a 1 K stable layer at 500–700 m is decoupled: surface-driven to
552 m, cloud-top layer 601–1246 m, V 1.05 m/s. The variance cover is
1.35·10⁻³ at −3σ, 0.5 at zero deficit (5·10⁻¹⁰ off), 0.9772 at 2σ; the
humidity PDF gives the zero-deficit layer 0.405. On 362 random bl34
columns (39 stable, 71 surface-driven, 43 decoupled, 209 coupled) the
engines agree on every regime, V to 2.3·10⁻⁷, w_e to 1.8·10⁻⁴, the
coefficients to 1.2·10⁻⁴ of each column's largest, the depths to
1.4·10⁻² m, θ after the step to 1.3·10⁻⁴ K, q to 6.9·10⁻⁸, q_c to
2.4·10⁻⁸, the wind to 9.6·10⁻⁵ m/s, under the defaults, the one-layer
jump, the M21 taper, no surface parcel and `coupledVeto` (no plume in a
coupled column; at most 1 of 362 columns parts where the dry
adjustment merges a near-neutral lowest layer in one engine only).
`test/gpuModel.test.mjs`: the variance cover moves the layer heating by
up to 2.80 K/day from the PDF's and its overcast blend by 4.04, the
engines by 4.6·10⁻⁴ against 24.4 K/day, the longwave heating by
2.7·10⁻³ against 135 W/m²; four N=6 steps of the deck under the moist
layer, jump-gated, regime-gated and bypassed: gates to 4.4·10⁻⁸, deck
water rms 6·10⁻⁴, 2–4 of 362 columns parting where a cloud top or parcel
crosses its threshold in one engine. `plumeConsumption` 'buoyant' (F from
the CAPE's layers only) agrees on 362 random columns to 6.5·10⁻⁵ K.

The deck. Three-day N=64 GPU screens from eight64_day0183 (day 186,
audited; low cloud the radiative cover below 680 hPa, `lowCover`):
gating the mixed-layer deck by the coupled regime (no bypass) gives SE
Pacific 0.147, Peru 0.014, Namibia 0.127; bypassing it 0.196 and 0.112
(cloud water below 680 hPa present) at albedo 0.362; under the 4 K jump gate 0.369, 0.481, 0.713. The resolved
boundary layer of the deck boxes is 0.6–0.9 humid at its top under
1.1–1.9 km inversions with 5–10 % of the columns coupled, so the
resolved stratus does not stand for the deck; the default keeps the
mixed-layer deck behind the 4 K jump (`deckRegime` 'inversion').

The plume on the new layer (ten-day N=64 GPU runs from eight64_day0183;
albedo days 186–193, ASR − OLR mean over days 188–193, 60–90N ice loss
over three days from nine64_day0091 in 10³ km³/day where run):

| run | albedo 186–193 | ASR − OLR | ice |
|---|---|---|---|
| parent (dry, plume, 3 h, CAPE0 70) | 0.330–0.361 | −7.5 | 0.183 |
| σ_s scale 3.2, λ 300, upper cloud 1 h, F buoyant | 0.302–0.343 | −10.7 | |
| scale 5, λ 300, upper cloud 1 h | 0.300–0.338 | −9.2 | 0.217 (with the blend) |
| scale 3.2, λ 300, cloud 1 h, upper 3 h | 0.337–0.369 | −8.4 | |
| scale 3.2, λ 300, cloud 1 h, upper 2 h | 0.320–0.356 | −8.9 | |
| scale 5, cloud 1 h, F buoyant | 0.333–0.367 | −7.8 | |
| scale 3.2, cloud 1 h, upper 2 h, coupled veto | 0.352–0.388 | −19.0 | |
| scale 5, λ 300, cloud 1 h | 0.292–0.330 | −7.3 | |
| scale 5, λ 300, cloud 1 h, CAPE0 120 | 0.304–0.334 | −9.1 | 0.247 (0.277 without the blend) |
| scale 1, λ 150, cloud 3 h, upper 1 h | 0.336–0.362 | −18.5 | 0.177 |
| scale 1, λ 150, cloud 2 h, upper 1 h | 0.331–0.361 | −17.5 | 0.189 |
| scale 5, λ 300 / 30 m stable, cloud 3 h, upper 1 h | 0.321–0.348 | −13.8 | 0.189 |
| scale 5, λ 300 / 30 m stable, cloud 2 h, upper 1 h | 0.310–0.345 | −12.4 | 0.202 |
| **defaults: scale 5, λ 300 / 30 m stable, cloud 1 h, CAPE0 120** | 0.307–0.340 | −10.5 | 0.233 |

Short cloud lifetimes cool less but thin the Arctic stratus; a short
upper lifetime lowers the albedo and raises the OLR as much. The coupled
veto gives the best deck cover (SE Pacific 0.343, Peru 0.442, Namibia
0.513, California 0.329 on day 193) at −19 W/m².

Acceptance on the defaults (copies of the states, GPU, `everySteps` 8):

| | value | asked |
|---|---|---|
| N=64 albedo days 184–193 | 0.287, 0.312, 0.313, 0.319, 0.326, 0.332, 0.340, 0.323, 0.316, 0.307 | 0.30–0.32 from 186 |
| N=64 ASR − OLR days 188–193 (mean) | −11.1, −13.1, −15.6, −10.5, −7.9, −4.9 (−10.5) | ±4 |
| N=64 global rain days 187–193, mm/d | 2.33–2.76 | 2.4–2.8 |
| day 193 low cloud: SE Pacific, Peru, Namibia, California | 0.244, 0.354, 0.525, 0.071 | 0.5–0.7 |
| their resolved inversion, m | 1885, 1863, 1702, 762 | 1000–1500 |
| their low-cloud water in cloud, g/m² | 230, 152, 101, 49 | 50–150 |
| SE Pacific rain, mm/d (10-day log; day-193 window) | 0.53 (2.12) | < 0.3 |
| Pacific ITCZ rain; heating peak | 8.41 mm/d; 438 hPa | 6–9; 400–500 |
| zonal-mean rain peak | 7.99 mm/d at 10.5N | 5–7 at 5–10N |
| N=128 albedo days 184–188 | 0.253, 0.271, 0.286, 0.301, 0.310 | 0.29–0.32 on 186–188 |
| N=128 ASR − OLR (mean) | +7.5, +3.0, −0.7, −5.1, −7.5 (−0.6) | ±4 |
| N=128 day 188 low cloud: SE Pacific, Peru, Namibia, California | 0.247, 0.432, 0.430, 0.397 | 0.5–0.7 |
| N=128 equatorial stress 160E–100W after 5 days, N/m² | −0.039 (parent −0.026) | toward −0.04 |
| 60–90N ice loss, 10³ km³/day | 0.233 (9.191 → 8.491) | ≤ 0.18 |
| fresh start day 30: albedo; ASR − OLR | 0.323; −7.7 (days 25–30 −5.2 to −7.7) | 0.29–0.33; ±10 |

The parent on the same runs: N=64 albedo 0.330–0.361 from day 186, mean
ASR − OLR −7.5, day-193 low cloud 0.280, 0.371, 0.359, 0.197; N=128
albedo 0.298, 0.314, 0.325 on days 186–188, mean +4.4, stress −0.026.
N=128 day 188: ITCZ 6.89 mm/d at 438 hPa, global rain 2.95 mm/d, the
zonal peak 14.6 mm/d at 10.5N (one eight-step window). Fresh start: no
NaN, clamped 0; equator after 30 days −0.58 m/s at 160E–100W, stress
−0.050 N/m².

The equator. The budget of Oct 1 named no boundary-layer term, so none
was changed. On nine128_day0183 after one CPU step the stress over
160E–100W by 20° is −0.0031, +0.0038, −0.0096, −0.0509, −0.0800,
−0.0586 N/m² against −0.0031, +0.0038, −0.0096, −0.0506, −0.0789,
−0.0578 under the dry scheme; w_e there is 6–25 mm/s against 2–11.

Pace (`js/gpu/profile.module.js`, 128 steps from eight128_day0183 after
64, alternated twice with the parent, nothing else on the GPU): a step's
median 86.9 and 86.9 ms against 84.5 and 84.6 (+2.8 %); the physics and
boundary-layer passes 9.36 against 7.65 ms, the adjust pass 15.15
against 14.6. The five N=128 days took 1.5–1.8 min a day beside two
N=64 runs.

What still misses: the N=64 balance (−10.5 against the parent's −7.5)
and its albedo above 0.32 on four of eight days; the deck boxes' low
cloud (0.24–0.53 at N=64, 0.25–0.43 at N=128) with the inversion at
1.7–1.9 km; the Arctic ice loss (0.233); the SE Pacific drizzle.

Review (Oct 1). Rerun from a copy of eight64_day0183, the ten-day N=64
acceptance run reproduces `m22acc64` byte for byte, its log and its
day-193 state. Twenty-four CPU steps at N=4 on bl34 give the parent's
state digest under `turbulence` 'dry' with `cloudLifetime` 3 h and
`plumeCape` 70, and a different one under 'dry' alone. The full suite
(45 files, run concurrently) passes with nothing skipped.

One CPU step of every column: on m22acc64_day0193 the boundary-layer
mix keeps each mixed column's mass-weighted θ_l to 7.4·10⁻¹⁶ and q_t to
8.4·10⁻¹⁶ relative (40,948 columns), leaves every layer outside the
coefficients' reach exactly as it was, and keeps each edge's column
momentum to 1.3·10⁻¹⁵ of Σm|u| (122,880 edges); eight64_day0183 the
same to 7.6·10⁻¹⁶, 8.6·10⁻¹⁶ and 1.0·10⁻¹⁵, nine128_day0183 (163,772
columns) to 9.0·10⁻¹⁶, 9.3·10⁻¹⁶ and 1.1·10⁻¹⁵. After the mix, the
single linearized saturation adjustment of the moist step condenses
34.05 g/m² per mixed column in the mixed layers against 33.82 when it
is iterated to convergence (+0.7 %), the driest cloudy layer at 0.961
humidity.

The cloud-top scale on m22acc64_day0193: in the 11,798 columns with a
cooling cloudy run, ΔF and V recomputed from the radiation's `longwave`
(W/m², positive heating) agree with the diagnosis to 1.9·10⁻¹⁴; ΔF is
3.3, 51.8 and 105.2 W/m² and V 0.32, 0.99 and 1.46 m/s at the 10th,
50th and 90th percentiles; 1,642 runs whose longwave sum warms get no
profile. In coupled columns w_e is 0.77, 6.89 and 32.75 mm/s at the
same percentiles, above 10 mm/s in 0.40 of them and at the 5 cm/s cap
in 0.030.

The variance cover on the same state, over the 46,520 cloudy layers
below the mixing top: none non-finite, all within the 0.01 floor and 1,
0.028 overcast, 369 at the floor and 812 below one half, 0.483
between 0.5 and 0.6; σ_s 0.027, 0.21 and 1.07 g/kg and Q₁ 0.020, 0.225 and 1.40 at the
10th, 50th and 90th percentiles; the floor binds on 0.025 and the 30 m
length on 0.550 of them. nine128_day0183: 283,462 layers, 0.053
overcast, 0.427 between 0.5 and 0.6, σ_s median 0.32 g/kg.

The regime against the deck's gate, the same step, open sea (27,581
columns): stable 0.145, surface-driven 0.587, decoupled 0.165,
coupled 0.103; the gate is open on 0.27, 0.08, 0.30 and 0.12 of
each, so of the open gates 0.082 are coupled. SE Pacific coupled
0.039 (gate open on 0.14 of them), Peru 0.080 (0.53), Namibia 0.253
(0.88), California 0.022 (none). The deck's water is not in q_c, so
where the deck runs the boundary layer finds no cloud top of its own.

The radiation's variance cover reads the boundary layer's surface
buoyancy flux of the step before; states did not carry it and the GPU
model's load zeroed it, so the first step after a load used the 30 m
length everywhere. Spin-up states now carry it (`boundaryBuoyancy`)
and the GPU model mirrors it (`boundaryLayer.buoyancyFlux`). The
in-day restart of `test/asyncSpinup.test.mjs` at N=6 is no closer for
it: the largest daily-forcing difference is 3.1 % with it and 3.0 %
without (netFlux, day 3), 2.4 % under `boundaryCover` 'pdf' and 1.5 %
for the parent, so the test keeps 4 %. `scripts/verticalAudit.mjs` now
starts its window from the state's depth, mixing top, regime and
buoyancy flux, moving the saved heights, which the GPU measures from
the surface, to the CPU's heights above sea level (on
m22acc64_day0193, 2,254 of the 2,927 cells above 1 km would otherwise
have their mixing top below their lowest layer); its first window step
had run with no mixing top. Day 193 of `m22acc64` re-audited: low-cloud
radiative cover 0.256, 0.367, 0.564, 0.072 (SE Pacific, Peru, Namibia,
California), in-cloud water 226, 158, 103, 48 g/m², global rain 2.75
mm/d, the gate replica exact over 6,696 column-steps.

five64_day2281 (27 layers), six64_day1004, seven64_day0639,
nine64_day0091, m21b64_day0183, eight64_day0183 and nine128_day0183
load and take two CPU steps with nothing non-finite; five64_day2281
and nine64_day0091 take 32 GPU steps and save regimes and mixing tops
with nothing non-finite.

**Tuning (Oct 1).** Run ten (a fresh atlas start on 976b0bc, N=128) rested
its decks at the inversion ceiling the cumulus lift, near 2.0–2.1 km, with
cloud layers 800–900 m thick and 400–900 g/m² of water (the radiation takes
at most 150), lost its Arctic ice by day 150 and drove trades of 0.08–0.11
N/m². Three changes, each on both engines with parity tests. Measured with
three-day N=64 GPU runs (`everySteps` 8) from copies of eight64_day0183 and,
for the Arctic, nine64_day0091, the changes taken in turn on top of each
other; day 186's albedo and ASR − OLR are the means of eight samples a day
(eight evenly spaced samples a day; the day line's own value is one step
at a fixed UTC and moves by up to 0.04 in albedo against the day mean), the
deck boxes and the global rain and evaporation from `scripts/verticalAudit.mjs`
on each day-186 state.

- The deck's rest by regime (`deckRest` 'regime', the default; 'inversion'
  and 'depth' remain). A coupled column rests at its boundary-layer top; a
  surface-driven or decoupled column whose ceiling's θ_v jump lies above
  `cumulusCeiling` (2000 m above the surface), or that has no ceiling,
  shuts its gate (G = 0), runs no deck and rests at the boundary-layer top;
  every other column rests at the ceiling. With the deck running, its water
  is not in q_c, so its column reads as surface-driven, and the test falls
  on the height of the inversion. Screens of the stand-down height, end of
  day 186 (albedo, ASR − OLR): 1500 m 0.290, −1.1; 2000 m 0.295, −2.4;
  2500 m 0.307, −6.4; unmodified 0.313, −8.1. At 1500 m the deck runs on
  0.038, 0.212 and 0.137 of the SE Pacific, Peru and Namibia column-steps.
- A lifetime for stratiform cloud (`stratiformLifetime` 3 h, the moist
  physics; null: none). A layer's lifetime moves from `cloudLifetime` (or
  the upper one) to it by a share that is 0 at and below the top of a
  plume that ran in the column this step, and elsewhere the larger of the
  cell's sea-ice cover and, under the moist boundary layer, 1 below the
  mixing top of a coupled column, 0 below that of any other, and the EIS
  share above it.
- The sea's drag coefficient 1.2·10⁻³ at the lowest layer's midpoint
  (`SEA_DRAG`, about 1.4·10⁻³ at 10 m; land keeps 1.5·10⁻³), the heat and
  vapour exchange following it; the GPU ocean's stress takes each cell's
  coefficient as the CPU's does. `SURFACE` (JSON) passes it to the
  spin-up and the audit.

| day 186, N=64 from eight64_day0183 | unmodified | regime rest | + stratiform lifetime | + sea drag |
|---|---|---|---|---|
| albedo, 8 samples (end of day) | 0.294 (0.312) | 0.283 (0.296) | 0.287 (0.298) | 0.279 (0.298) |
| ASR − OLR, W/m² | −2.0 (−8.0) | +1.3 (−2.8) | +0.8 (−2.8) | +3.2 (−2.9) |
| global rain; evaporation, mm/d | 2.25; 2.33 | 2.20; 2.28 | 2.21; 2.28 | 2.07; 2.12 |
| deck runs, SE Pacific / Peru / Namibia | 0.461 / 0.657 / 0.704 | 0.182 / 0.500 / 0.514 | 0.187 / 0.531 / 0.545 | 0.191 / 0.582 / 0.797 |
| deck height where it runs, m | 1694 / 1345 / 1431 | 1458 / 1299 / 1374 | 1460 / 1339 / 1378 | 1385 / 1236 / 1313 |
| its cloud layer, m | 505 / 341 / 450 | 296 / 295 / 387 | 310 / 333 / 396 | 261 / 232 / 303 |
| its water, g/m² (capped; share capped) | 267 / 124 / 225 (82 / 76 / 129; 0.41 / 0.28 / 0.71) | 97 / 107 / 180 (69 / 74 / 108; 0.21 / 0.18 / 0.50) | 107 / 131 / 184 (76 / 77 / 114; 0.26 / 0.24 / 0.55) | 61 / 60 / 113 (51 / 46 / 81; 0.11 / 0.10 / 0.27) |
| low cloud, radiative | 0.414 / 0.523 / 0.694 | 0.254 / 0.448 / 0.637 | 0.281 / 0.447 / 0.666 | 0.233 / 0.416 / 0.711 |
| SE Pacific rain, mm/d | 0.70 | 0.58 | 0.58 | 0.48 |
| 60–90N ice loss from nine64_day0091, 10³ km³/day | 0.233 | 0.225 | 0.175 (9.191 → 8.666) | 0.174 |

The equatorial sea after one CPU step from ten128_day0183, 2S–2N
160E–100W: stress −0.1102 → −0.0881 N/m² (5S–5N −0.1062 → −0.0849), the
lowest layer's wind −7.61 m/s in both, the 10 m wind by the log law −7.08
→ −7.13 m/s. Three N=64 days from run ten's own ten64_day0183, unmodified
and with all three: albedo 0.333 → 0.313, ASR − OLR −12.1 → −5.6 W/m²,
global rain 2.83 → 2.63 and evaporation 2.87 → 2.61 mm/d; the SE Pacific
deck runs on 0.532 → 0.494 of the column-steps (stood down on 0.243) at
1617 → 1567 m, 546 → 485 m thick with 337 → 273 g/m² (capped on 0.631 →
0.545); Namibia 0.908 → 0.905 at 1192 → 1112 m, 92 → 47 g/m².

After the three, one CPU step of tf3m_day0186: `moist.adjust` keeps each
column's c_p T + L q to 1.1·10⁻¹⁵ and its water with the rain to 8.5·10⁻¹⁶
relative but for 27 columns the filler touches. Parity: the regime rest
under moist turbulence at `cumulusCeiling` 3000 and 200 m (height rms
8.7·10⁻⁷, gate 4·10⁻⁸); the stratiform lifetime on 362 random columns with
random regimes, mixing tops, EIS shares and sea ice (θ to 6.5·10⁻⁵ K, q to
3.6·10⁻⁸, q_c to 1.3·10⁻⁸; 338 of 1025 cloudy layers keep more cloud under
3 h than 1 h in both engines); with the sea's drag at zero no stress
reaches the ocean across any of the 670 sea–sea edges at N=6 on either
engine. The in-day restart of `test/asyncSpinup.test.mjs` now differs by
4.6 % (netFlux, day 3; under 4 % with any one change undone), and the test
keeps 5 %.

What still misses: at 2000 m the SE Pacific deck runs on 0.19 of its
column-steps and its low cloud is 0.23–0.28: on 0.55 of the steps the
4 K jump above its boundary layer lies above 2000 m or is missing (the
resolved inversion, the interface of largest dθ_v/dz, at 1.56–1.71 km);
on run ten's state the SE Pacific deck is still 485 m thick at 1.57 km with
273 g/m² of water, a mixed layer moist enough that its cloud base sits near
1.1 km. The sea drag lowers the global rain and evaporation of the
eight64 runs from 2.21 and 2.28 to 2.07 and 2.12 mm/d (Earth 2.6–2.8).

**Sweep (Oct 1).** A scored perturbed-parameter sweep over eleven
parameters on short runs, `scripts/sweep/` (`sweep.mjs` the design and the
screens, `score.mjs` the score, `fit.py` the response surface,
`candidates.mjs` the longer tests; the runs' outputs under
`runs/sweep/`). Every albedo and ASR − OLR is a day mean of eight samples
(eight evenly spaced samples a day).

The score is Σ wᵢ eᵢ², eᵢ = (xᵢ − targetᵢ)/toleranceᵢ, over: from three
N=64 days from eight64_day0183, day 186's ASR − OLR (0 ± 3 W/m², weight
4) and albedo (0.30 ± 0.015, 2); from `scripts/verticalAudit.mjs` on the
day-186 state, the global rain (2.7 ± 0.2 mm/d, 1), the SE Pacific and
Peru radiative low cloud (0.6 ± 0.1, 1 each), their deck water as the
radiation takes it (100 ± 50 g/m², 1 each) and their rain (0.2 ± 0.2
mm/d, 1 each), the Pacific ITCZ's rain (7.5 ± 1.5 mm/d, 1) and heating
peak (450 ± 50 hPa, 0.5), the zonal-mean rain peak's latitude (7.5N ±
2.5, 0.5); the day line's equatorial stress over 160E–100W (−0.05 ±
0.015 N/m², 1); and from three N=64 days from nine64_day0091 the 60–90N
sea-ice volume lost per day (0.15 ± 0.03 10³ km³/day, 2). The balance
weight is 4, not 3: the attribution of run ten's deficit found it about
90 % OLR (the weaker longwave effect of the lost high cloud), which only
the ASR − OLR term sees. The global rain is the audit's at the end of day
186; the day mean of day 186 (1.7–1.9 mm/d) still climbs from the
eight64 start's 0.7 mm/d on day 184 and is kept beside it.

The design: a maximin Latin hypercube of 40 points (seed 20261001, least
distance 0.772 in the unit cube) and the defaults as point 0, over the
variance cover's `varianceScale` (2–10), `mixingLength` (150–600 m) and
`stableMixingLength` (10–60 m), `stratiformLifetime` (1–6 h),
`cloudLifetime` (0.5–2 h), `plumeEntrainment` (0.05–0.2), `plumeCape`
(40–200 J/kg), `minimumInversion` (2–6 K), `criticalHumidity` (0.7–0.9),
`cumulusCeiling` (1500–2500 m) and the sea's drag (1.0–1.5·10⁻³). The
screen scores run from 18.3 (point 17) to 148.4, the defaults 49.0. A
full quadratic in the parameters coded to [−1, 1] (78 coefficients) by
ridge regression, the penalty 0.79 by leave-one-out: R² 0.991,
leave-one-out RMSE 21.8 against a spread of 30.6 (Q² 0.49). Along each
range with the others at the defaults (slope over the half range;
surface at low / default / high):

| parameter | slope | low / default / high |
|---|---|---|
| minimumInversion | +25.7 | 32.3 / 52.5 / 83.8 |
| sea drag | −20.5 | 72.6 / 52.5 / 36.4 |
| cumulusCeiling | −14.2 | 71.1 / 52.5 / 42.7 |
| mixingLength | +14.0 | 41.6 / 52.5 / 64.9 |
| cloudLifetime | −12.5 | 62.6 / 52.5 / 43.2 |
| criticalHumidity | +8.4 | 47.8 / 52.5 / 64.5 |
| varianceScale | +4.6 | 49.8 / 52.5 / 60.5 |
| plumeCape | +2.1 | 59.4 / 52.5 / 63.7 |
| plumeEntrainment | −1.8 | 59.1 / 52.5 / 71.9 |
| stableMixingLength | −1.7 | 59.5 / 52.5 / 63.1 |
| stratiformLifetime | −1.2 | 54.4 / 52.5 / 53.2 |

By term (linear, change of e over the full range): the deck boxes' low
cloud is the minimum jump's (SE Pacific −6.1, Peru −8.1); the albedo
the short lifetime's (+2.6) and the jump's (−1.5); the balance the
jump's (+2.5) and the short lifetime's (−1.8); the Arctic loss the
stratiform lifetime's (−2.2); the global and ITCZ rain the entrainment's
(−1.6, −2.4) and the drag's (+1.3, +1.3); the stress the drag's (−0.5).
The heating peak and the zonal rain peak are fitted with R² 0.22 and 0.24.

The candidates: the four best screens (points 17, 34, 16, 40), the
surface's minimum over the box (qmin: ten of its eleven parameters at an
edge of their range, cloudLifetime inside it at 1.579 h; predicted
−61.7), the minimum of the score composed from each
term's own fit (cmin, composed 11.1) and the defaults (base), each on ten
N=64 days from eight64_day0183 (albedo over days 186–193, ASR − OLR over
188–193, the audit of day 193), five N=128 days from eight128_day0183
(186–188) and thirty N=64 days from a fresh atlas start on bl34 (day 30),
scored by the same terms (the Arctic from the screen, the N=128 runs by
balance, albedo and stress, the fresh starts by balance and albedo).
qminr and baser are qmin and base with the drag 10⁻⁴ of itself larger,
run to measure the scores' noise.

| | qmin | qminr | base | baser | p17 | cmin | p16 | p40 | p34 |
|---|---|---|---|---|---|---|---|---|---|
| score: total (day 193 / N=128 / day 30) | 122 (70 / 11 / 41) | 283 (179 / 11 / 93) | 129 (52 / 23 / 54) | 193 (57 / 24 / 112) | 131 (58 / 4 / 68) | 149 (69 / 14 / 66) | 162 (88 / 1 / 72) | 204 (102 / 1 / 101) | 207 (123 / 2 / 81) |
| screen score | 41.7 | 41.4 | 49.0 | 47.8 | 18.3 | 39.5 | 32.5 | 32.7 | 22.2 |
| ASR − OLR 188–193, W/m² | −3.5 | −3.3 | +1.5 | +1.3 | +0.2 | −4.4 | −3.6 | −3.9 | −0.8 |
| albedo 186–193 | 0.315 | 0.314 | 0.288 | 0.289 | 0.300 | 0.310 | 0.314 | 0.322 | 0.301 |
| global rain, audit day 193 (day means 188–193), mm/d | 2.57 (2.64) | 2.60 (2.64) | 2.42 (2.51) | 2.42 (2.50) | 2.55 (2.58) | 2.74 (2.70) | 2.68 (2.53) | 2.51 (2.48) | 2.55 (2.54) |
| evaporation, mm/d | 2.61 | 2.64 | 2.48 | 2.49 | 2.51 | 2.67 | 2.60 | 2.47 | 2.49 |
| low cloud SE Pacific / Peru / Namibia | 0.32 / 0.49 / 0.41 | 0.38 / 0.50 / 0.38 | 0.25 / 0.39 / 0.39 | 0.22 / 0.40 / 0.33 | 0.23 / 0.44 / 0.38 | 0.42 / 0.60 / 0.29 | 0.23 / 0.36 / 0.59 | 0.22 / 0.30 / 0.40 | 0.30 / 0.50 / 0.41 |
| deck water as radiated, SE Pacific / Peru, g/m² | 143 / 146 | 142 / 143 | 124 / 137 | 138 / 142 | 144 / 145 | 144 / 147 | 146 / 147 | 149 / 147 | 148 / 147 |
| SE Pacific deck runs; resolved inversion, m | 0.234; 1955 | 0.270; 1952 | 0.014; 1889 | 0.023; 1887 | 0.024; 2002 | 0.276; 1884 | 0.062; 1941 | 0.022; 1911 | 0.118; 1951 |
| SE Pacific / Peru rain, mm/d | 1.38 / 0.00 | 1.38 / 0.00 | 1.18 / 0.00 | 1.23 / 0.00 | 1.37 / 0.00 | 1.60 / 0.03 | 1.60 / 0.00 | 1.73 / 0.01 | 1.33 / 0.02 |
| Pacific ITCZ rain, mm/d; heating peak, hPa | 8.26; 438 | 7.90; 438 | 6.81; 438 | 6.52; 439 | 8.41; 438 | 6.52; 438 | 7.62; 439 | 9.93; 438 | 7.10; 438 |
| zonal rain peak, °N | 2.5 | −30.5 | 10.5 | 10.5 | 7.5 | 3.5 | 11.5 | 10.5 | 37.5 |
| equatorial stress day 193, N/m² | −0.026 | −0.025 | −0.024 | −0.021 | −0.024 | −0.026 | −0.025 | −0.027 | −0.023 |
| 60–90N ice loss (screen), 10³ km³/day | 0.219 | 0.219 | 0.174 | 0.174 | 0.171 | 0.159 | 0.197 | 0.152 | 0.180 |
| N=128: albedo / ASR − OLR 186–188; stress | 0.308 / −4.7; −0.034 | 0.307 / −4.5; −0.034 | 0.270 / +5.6; −0.035 | 0.270 / +5.7; −0.034 | 0.287 / +2.0; −0.036 | 0.304 / −5.2; −0.032 | 0.296 / −0.6; −0.037 | 0.302 / −0.5; −0.036 | 0.289 / +0.8; −0.035 |
| fresh day 30: albedo / ASR − OLR (mean of days 25–30) | 0.338 / −8.0 (0.336 / −7.2) | 0.354 / −12.3 (0.350 / −10.3) | 0.339 / −9.6 (0.334 / −9.4) | 0.350 / −14.2 (0.348 / −14.0) | 0.349 / −10.3 (0.349 / −10.7) | 0.343 / −10.6 (0.352 / −14.3) | 0.347 / −10.9 (0.354 / −12.5) | 0.364 / −12.1 (0.360 / −11.5) | 0.348 / −11.7 (0.343 / −10.4) |

No run made a NaN or clamped the ocean. The replicates show the noise:
the ten-day and N=128 terms but one repeat to a few points, while the
zonal rain peak of a single eight-step audit moved from 2.5N to 30.5S in
qmin's pair (+113 points), and the fresh start's day 30 fell by 4–5 W/m²
in both replicates (+40 and +58 points). With the zonal peak left out and
day 30 replaced by the mean of days 25–30, qmin and qminr score 113 and
144, base and baser 123 and 188. The four best screens did not keep
their lead: the three-day screens score a state still recovering from the
eight64 start (day 184's day-mean albedo 0.25–0.27, rain 0.7–0.8 mm/d),
and by day 193 their SE Pacific decks ran on 0.02–0.12 of the
column-steps with 1.3–1.7 mm/d of rain, p16 and p40 had brightened to
0.314 and 0.322, and p34 took a zonal rain peak at 37.5N. The winner is
qmin, the lowest score on the candidates' runs and on both readings of
the replicates, set as the defaults on both engines: `varianceScale` 10,
`mixingLength` 150 m, `stableMixingLength` 10 m, `minimumInversion` 2 K,
`criticalHumidity` 0.7, `cumulusCeiling` 2500 m, `cloudLifetime` 1.579 h,
`stratiformLifetime` 1 h, `plumeEntrainment` 0.05, `plumeCape` 40 J/kg,
`SEA_DRAG` 1.5·10⁻³. Three N=64 days from eight64_day0183 under the new
defaults reproduce qmin's screen byte for byte.

What still misses with the winner: the fresh start's deficit (day 30
−8.0 W/m², albedo 0.338; −12.3 and 0.354 in the replicate); the Arctic
loss, 0.219 against 0.15 (the stratiform lifetime is now shorter than
the short one, and its term is the Arctic's); the SE Pacific deck (low
cloud 0.32, deck on 0.23 of the column-steps under a 1.95 km inversion,
rain 1.38 mm/d against 0.1–0.3); the deck water as the radiation takes it
near the 150 g/m² cap in both boxes; the zonal rain peak at 2.5N; the
equatorial stress at N=64 of −0.026 N/m² against −0.05 (and −0.034 at
N=128 from eight128). Ten of the winner's eleven parameters sit on an
edge of their range, where the fit predicted −61.7 against the 122 the
runs scored, and the 2 K jump is the weak test the deck's own rest was
moved away from: it opens the deck's gate (mlmGate > 0.5) over 0.208 of
the globe and 0.112 of 10S–10N on day 186 and over 0.182 and 0.155 on
day 193, against the old defaults' 0.033 and 0.021, 0.044 and 0.046.
The deck water the score reads is the radiation's, capped at 150 g/m²;
the deck's own water path where it runs on day 193 rose from base's 381
(SE Pacific) and 472 g/m² (Peru) to qmin's 411 and 673, its cloud layer
from 546 and 638 m to 653 and 778 m. Scored on the deck's own water path
with the same target and tolerance, the candidates total base 215, p17
276, p16 283, cmin 288, qmin 291, p34 357 and p40 376.

The winner was not adopted. Its margin over the fixes' defaults (122
against 129 on the specified score) is inside the replicate noise (one
replicate scored 283), ten of its eleven values sit on an edge of the
design, it regresses the Arctic loss to 0.219 from 0.174 by making the
stratiform lifetime shorter than the convective one, and its 2 K
inversion test opens the deck's gate over a fifth of the globe. The
score also could not see the deck's own water path, which the 150 g/m²
cap hides, and every candidate missed the fresh start's deficit of
−8 to −14 W/m² at day 30, so the eleven parameters do not reach that
term. The defaults stay at the fixes' values (variance scale 5, mixing
lengths 300 and 30 m, cloud lifetime 1 h and stratiform 3 h, minimum
inversion 4 K, critical humidity 0.8, sea drag 1.2·10⁻³, cumulus
ceiling 2000 m). The sweep's tools stay: the next score should use the
deck's uncapped water path and thickness, average the fresh start's
days 25–30 and several audit windows, constrain the stratiform
lifetime to at least the convective one, and add the parameters of the
cloud above the boundary layer, since the fresh start's deficit is
longwave.

**Cloud optics (Oct 1).** The radiation gave every cloud the gray optics
τ' = 95 m²/kg × W and κ = 130 m²/kg. `scripts/cloudClasses.mjs` splits the
day-mean cloud effects of a state by class: one CPU step, then the state
held and lit at 24 instants of its day, once with all cloud and once
without each class (the radiation's cloud mask, the deck held at its
diagnosis at the day's start; holding it moves the global ASR − OLR by
0.01 W/m²); resolved cloud classed by the top of its run of cloudy layers.
Before: the gray optics on the day-186 state of the three-day N=64 run
from eight64_day0183 (the run's day means: SWCRE −56.0, LWCRE 17.6). After:
the phase optics below on the day-186 state of the same run under them
(SWCRE −36.0, LWCRE 15.7); on the before state they give −38.2 and 15.6.
Cover, effects W/m², grid-mean path g/m² warmer than 273 K / 273–235 K /
colder, in-cloud path, in-cloud visible τ of the water's physical optics
with its cover shares below 3.6 / 3.6–23 / above 23, and the two-stream's
in-cloud depth τ':

| global | cover | SW | LW | grid path | in cloud | τ (shares) | τ' |
|---|---|---|---|---|---|---|---|
| deck, before | 0.015 | −2.0 | 0.1 | 1.2 / 0.4 / 0.0 | 101 | 12.4 (0.18 / 0.82 / 0.00) | 9.6 |
| deck, after | 0.015 | −1.3 | 0.1 | 1.1 / 0.3 / 0.0 | 100 | 12.2 | 1.7 |
| cumulus, before | 0.010 | −2.0 | 0.2 | 3.2 / 1.9 / 0.0 | 501 | 65.5 (0.00 / 0.21 / 0.79) | 47.6 |
| cumulus, after | 0.010 | −1.9 | 0.2 | 3.3 / 1.9 / 0.0 | 504 | 66.4 | 9.7 |
| low (top below 680 hPa), before | 0.236 | −28.5 | 3.3 | 11.1 / 20.1 / 0.1 | 132 | 15.6 (0.14 / 0.66 / 0.20) | 12.6 |
| low, after | 0.237 | −20.5 | 3.5 | 11.5 / 20.3 / 0.0 | 135 | 15.9 | 2.3 |
| middle (680–440 hPa), before | 0.067 | −5.5 | 1.9 | 3.4 / 11.5 / 0.1 | 224 | 26.2 (0.27 / 0.38 / 0.35) | 21.3 |
| middle, after | 0.067 | −4.3 | 2.0 | 3.5 / 11.7 / 0.1 | 230 | 27.1 | 4.0 |
| high (above 440 hPa), before | 0.144 | −13.4 | 10.7 | 5.0 / 26.9 / 2.3 | 236 | 25.4 (0.44 / 0.26 / 0.30) | 22.5 |
| high, after | 0.135 | −8.0 | 8.6 | 5.2 / 27.3 / 2.0 | 255 | 27.6 (0.43 / 0.25 / 0.32) | 4.3 |
| all, before (sum of classes) | 0.403 | −57.2 (−51.4) | 17.0 (16.2) | 23.9 / 60.7 / 2.4 | | | |
| all, after | 0.399 | −38.1 (−36.1) | 15.1 (14.4) | 24.6 / 61.6 / 2.1 | | | |
| Earth | 0.65–0.68; high 0.2–0.3 | −47 | +26 | | | high: about 0.6 / 0.3 / 0.1 | |

By phase (the ramp below), before: liquid 63.4 and ice 23.5 g/m² globally,
66.2 and 24.8 over the sea (Earth: liquid over the oceans 50–90, O'Dell et
al. 2008; ice 20–70, CloudSat). By region, SWCRE / LWCRE / cover before →
after: 30S–30N −43.7 / 11.0 / 0.223 → −25.2 / 9.2 / 0.221 (high cover 0.104
→ 0.097, low −16.6 → −10.7); 30–60N −55.7 / 15.9 / 0.392 → −39.8 / 15.0 /
0.391; 30–60S −110.0 / 25.5 / 0.629 → −80.5 / 23.8 / 0.627 (low −70.0 →
−53.0 at cover 0.48 and 173 g/m² in cloud); 60–90N −40.1 / 29.2 / 0.741 →
−29.8 / 25.5 / 0.712; 60–90S −35.2 / 30.0 / 0.825 → −22.9 / 26.1 / 0.812;
land −39.9 / 15.6 / 0.383 → −29.0 / 14.7 / 0.381; sea −64.2 / 17.6 / 0.411
→ −41.8 / 15.3 / 0.407. The shortwave excess was low cloud's, half of the
total and −70 of −110 over 30–60S: its 132 g/m² in cloud had τ' 12.6 under
the gray optics against 2.3 from its water's own optics (τ 15.6, g 0.86),
reflecting 0.93 of a beam at μ = 0.5 instead of 0.70. The missing
longwave is high cloud's: 0.144 of cover against 0.2–0.3 (0.104 in the
tropics), and the total cover 0.40 (tropics 0.22) against 0.65–0.68. The
high cloud is thick: 236 g/m² and τ 25 in cloud, 0.30 of its cover above
τ 23 against ISCCP's deep-convective tenth (Rossow and Schiffer 1999), its
gray τ' 22.5 that of τ ≈ 160 liquid. The condensate's totals are Earth's;
the cloud is too little in cover and too thick, and ice colder than
235 K holds 2.4 g/m².

The optics (radiation header; both engines; `cloudOptics`): liquid share
linear in the layer's temperature from 0 at 235.15 K to 1 at 273.15 K
(half at −19 °C; CALIPSO's supercooled half near −20 °C, Hu et al. 2010);
liquid τ = 3W/(2ρ_w r_e), r_e 11.8 µm over sea and ice sheets and 8.5 µm
over land (Han et al. 1994), g from Slingo (1989), κ = 1.66 × 0.090361
m²/g (CAM3); ice r_e = D_e/2 from Ou and Liou (1995) at −60 to −20 °C
(15.55–73.55 µm), τ and g from Ebert and Curry (1992), κ = 1.66 (0.005 +
1/r_e) m²/g; the two-stream takes τ' = (1 − g)τ (Coakley and Chýlek 1975,
upscatter (1 − g)/2). Per kg/m²: liquid τ' 18.0 (sea) and 26.5 m²/kg
(land), κ 150 m²/kg; ice τ' 35.9 and κ 115 at −60 °C, τ' 7.0 and κ 31 at
−20 °C. The deck and the cumulus take their layer's optics.
`cloudScattering` and `cloudAbsorption` now default to null and, set,
restore the gray optics: with 95 and 130 the GPU's day lines from
eight64_day0183 repeat the parent's to every printed digit, and the fresh
start's days 6–10 M21's. Tests (`test/cloudOptics.test.mjs`): the optics
against hand-computed values; an overcast layer over a black surface at
μ 0.5 reflects 0.64304 (100 g/m² of liquid at 285 K) and 0.41806 (20 g/m²
of ice at 213.15 K) with emissivities 0.77687 (10 g/m² of liquid) and
0.89985, each to 10⁻⁹; physics alone at N=6 with 234 warm, 253 mixed and
275 cold cloudy layers and 280 empirical decks, the engines' layer heating
apart by 8.0·10⁻³ against 917 K/day, the surface flux, surface sunlight
and both effects by 6.8·10⁻³, 7.5·10⁻³, 7.0·10⁻³ and 9.1·10⁻⁴ W/m², the
phase optics moving the effects by up to 319 and 12.5; 179 sunlit columns
close to 3.4·10⁻¹⁶ of the beam and their layers to 2.8·10⁻¹⁵.
`test/cloudEffect.test.mjs`'s 48-step run keeps the gray optics: under
the phase optics one cell parts by 1.04 W/m² (rms 1.1·10⁻³).

Review (Oct 1). By hand from the sources, an overcast layer over a black
surface at μ 0.5, reflectance τ'/(τ' + 2μ) and emissivity 1 − exp(−κW):
100 g/m² of liquid at 285 K, sea τ 12.71, g 0.858, τ' 1.801, 0.64304;
land τ 17.65, τ' 2.645, 0.72568; emissivity 1.00000 in both. Ice at
220 K: r_e 21.26 µm, τ 117.8 and τ' 26.09 m²/kg, κ 86.4 m²/kg; at 250 K:
r_e 64.73 µm, 74.7, 11.94 and 79.3. 20 g/m² reflect 0.34289 and 0.19272
and emit 0.82234 and 0.79525; 100 g/m² reflect 0.72292 and 0.54414 and
emit 0.99982 and 0.99964. Both engines' columns give each of these to
the sixth digit, the GPU's from its own `cloudOptics` and `stream`. From
180 to 300 K in steps of 0.1 K the GPU's optics follow the CPU's to
1.6·10⁻⁶ relative, and no step between neighbours moves τ' by more than
0.20 or κ by more than 0.59 m²/kg. One step from eight64_day0183 on both
engines against the parent tree (15754a7) on the same state: with 95 and
130 the CPU's absorbed, reflected and outgoing fluxes, clear-sky fluxes
and surface sunlight are bit-identical cell by cell and the GPU's state
and day sums hash alike; under the phase optics the clear-sky ASR and
OLR (284.82 and 264.13 W/m² for that step) are identical in every cell;
in the 20480 lit columns absorbed and reflected sunlight add up to the
beam to 4.2·10⁻¹⁶ (CPU) and 2.0·10⁻⁷ (GPU, f32), no column has a
negative absorbed, reflected, atmospheric or surface term, and the
night columns have none. `scripts/cloudClasses.mjs` on the day-186 state
of the phase-optics run with every class taken away at once returns the
clear-sky fluxes to 0 W/m² in every cell and the total effects −38.098
and 15.130 W/m², against the single classes' sum −36.1 and 14.4 (30S–30N
−25.6 against −25.2; 30–60S −73.9 against −80.5). The three-day run from
eight64_day0183 and the Arctic run from nine64_day0091, repeated, give
the table's day-186 line and audit and 9.191 → 8.261 (0.310) again.
`scripts/sweep/runs.mjs` reads each parameter's base from
`PHYSICS_DEFAULTS`, where `cloudScattering` and `cloudAbsorption` are
now null: their base in `PARAMETERS2` is NaN.

Not built: an ice fall speed. The standard remedy for thick, short-lived
high cloud is sedimenting ice that sublimates below (Heymsfield and
Donner 1990), a moist-physics change on both engines. Three N=64 days from
eight64_day0183 under the phase optics with `upperCloudLifetime` (the
anvils' layers above 700 hPa) 1 h (the default, as `cloudLifetime`),
3 h and 6 h: SWCRE −36.0, −42.9, −49.4; LWCRE 15.7, 23.5, 31.3; ASR − OLR
+10.2, +10.7, +11.9 W/m²: the upper cloud now adds more longwave than
shortwave effect, where under the gray optics it traded them one for one.

Runs on the phase optics (N=64 GPU, `everySteps` 8; day means; the audit on
each last state):

| | gray | phase |
|---|---|---|
| eight64 day 186: albedo; ASR; OLR; ASR − OLR | 0.312; 234.3; 242.1; −7.8 | 0.253; 254.4; 244.2; +10.2 |
| SWCRE; LWCRE; clear-sky ASR; clear-sky OLR | −56.0; 17.6; 290.3; 259.7 | −36.0; 15.7; 290.4; 259.9 |
| rain day mean (audit), mm/d | 1.70 (2.04) | 1.78 (2.12) |
| SE Pacific / Peru / Namibia low cloud | 0.238 / 0.385 / 0.690 | 0.239 / 0.373 / 0.707 |
| SE Pacific deck water as radiated (its own), g/m²; rain, mm/d | 54.4 (65.0); 0.48 | 50.7 (62.7); 0.49 |
| Pacific ITCZ rain, mm/d; heating peak, hPa | 5.05; 439 | 5.05; 439 |
| ten64 day 186: albedo; ASR − OLR; SWCRE; LWCRE | 0.344; −15.8; −66.6; 19.0 | 0.278; +4.4; −44.3; 17.1 |
| rain day mean (audit); SE Pacific low cloud; its deck water (own) | 2.69 (2.62); 0.527; 104.1 (262.9) | 2.72 (2.67); 0.534; 107.9 (280.6) |
| SE Pacific rain; ITCZ rain; heating peak | 1.02; 3.74; 788 | 0.99; 3.56; 788 |
| 60–90N ice loss from nine64_day0091, 10³ km³/day | 0.166 (9.191 → 8.693) | 0.310 (→ 8.261) |
| its day-94 60–90N SWCRE; LWCRE; cover | −140.3; 26.2; 0.684 | −96.5; 22.8; 0.676 |
| fresh atlas start days 6–10: albedo; ASR; OLR; ASR − OLR | 0.459; 184.2; 220.5; −36.4 | 0.368; 215.2; 225.3; −10.1 |
| SWCRE; LWCRE; rain; sea surface sunlight | −104.3; 36.4; 4.53; 108.8 | −73.8; 32.5; 4.56; 143.1 |

No run made a NaN or clamped the ocean but nine64's known day-94 step. Pace
(`js/gpu/profile.module.js`, 128 steps after 64 from eight64_day0183,
alternated twice with the parent, nothing else on the GPU): step median
19.8 and 19.8 ms against 19.9 and 19.9; the physics and boundary-layer
passes 2.21 against 2.27 ms.

What still misses: the cover (0.40 against 0.65–0.68; high 0.135 against
0.2–0.3) and the high cloud's thickness; SWCRE −36.0 and LWCRE 15.7 against
−47 and +26; ASR − OLR +10.2 W/m² on eight64's day 186, with a clear-sky
OLR of 259.9 against CERES's 265–266; the Arctic loss 0.310 against the
0.15–0.18 set under the gray clouds, whose June SWCRE over 60–90N was
−140; the deck boxes, unchanged. A later sweep should vary the cloud's
amount and spread — `upperCloudLifetime` (or the ice fall speed that would
replace it), `cloudLifetime`, `stratiformLifetime`, `criticalHumidity`,
`boundaryCriticalHumidity` — and the droplet radii at most within the
observed 10–14 µm over sea and 7–10 µm over land; not `cloudScattering`
or `cloudAbsorption` (set, they turn the optics gray), the phase ramp, the
ice radius fit or the longwave coefficients.

**Cloud amount and spread (Oct 1).** `scripts/cloudRegimes.mjs` takes a
state one CPU step on and gives, for seventeen regimes, the total cover
under the radiation's overlap (and under maximum-random and
exponential-random overlap), each class's cover, grid-mean and in-cloud
path and ISCCP optical-depth shares, the relative humidity over water and
ice at 150–350, 350–700 hPa and below, the shares of humid layers that
hold no cloud, the cumulus updraught's area and the liquid and ice paths,
beside Earth's values (Klein and Hartmann 1993 for the deck boxes, the
others from memory as the script marks them). The diagnosis on the
day-186 state of the three-day N=64 run from eight64_day0183 under the
phase optics: total cover 0.40, warm pool 0.15, Pacific ITCZ 0.31, the
trades 0.04–0.21 with low cover 0.02–0.06 and a cumulus cover of
0.002–0.008, the storm tracks 0.50–0.75, land 0.29 and sea 0.36. A layer
held cloud only once its grid mean saturated over liquid water:
46 % of the layers below 700 hPa above 90 % humidity held none, and the
upper troposphere (150–350 hPa) stood at a relative humidity over ice of
0.73 globally, 0.99 over the warm pool and 1.14 over the ITCZ, with 47 and
68 % of their layer area above ice saturation and 88 and 81 % of that
cloud-free: at −40 to −60 °C water saturation lies 1.4–1.6 times ice
saturation, so cirrus could not form below it. The cloud that did form
filled its grid box (in-cloud high-cloud path 255 g/m², 2.0 g/m² of it
colder than 235 K) and was taken away by the 1 h lifetime.

The schemes, both engines (`js/physics/moist.module.js`, the radiation's
cover and overlap, the adjust and physics kernels of
`js/gpu/physics.gpu.js`):

- Saturation over ice (`iceSaturation`): the saturation vapour pressure of
  cloud mixes Bolton's over water and the IFS form over ice,
  e_i = 611.21 exp(22.587 (T − 273.16)/(T + 0.7)) Pa, by the optics' liquid
  share, linear from 235.15 K (all ice) to 273.15 K (all liquid), in the
  condensation, the evaporation of large-scale rain and the variance cover.
  Every phase change takes the latent heat of vaporisation, the fusion heat
  of snow is released at the surface as before.
- A uniform total-water distribution (`condensation` 'uniform'; LeTreut and
  Li 1991, whose fixed width is Sundqvist et al.'s 1989 cover): above the
  moist boundary layer's mixing top a layer holds the condensate
  (Q + b)²/(4b) of a deficit Q = a (q_t − q_s(T_l)) over a half-width
  b = a (1 − RH_c) q_s(T_l), all of Q above b, none below −b, and the
  radiation covers sqrt(q_c/b) of it, so that in-cloud water is sqrt(b q_c)
  and thin cloud is thin. RH_c = 0.75 + 0.225 exp(1 − (p_s/p)²), ECHAM6's
  crs 0.975, crt 0.75 and nex 2 at T63 (CAM3 uses 0.70–0.80 for its high
  clouds, the IFS 0.8 above σ 0.8). These are not observed values: ECHAM6
  sets crs, crt, nex and cvtfall by truncation in mo_echam_cloud_params
  (Stevens et al. 2013, with cvtfall 2.5, crs 0.975, crt 0.75 and nex 2 at
  T63), the values of its tuning (Mauritsen et al. 2012). Below the mixing top the variance cover
  and the adjustment to saturation stay. A clear layer at 600 hPa and
  0.95 of water saturation now holds 0.10 g/kg (`test/cloudIce.test.mjs`).
- Falling ice (`iceFall` 2.5 m/s, `iceFallExponent` 0.16): the ice share of
  each layer's cloud falls at v = 2.5 (ρ q_i/f)^0.16, the in-cloud content
  over the uniform cover, Heymsfield and Donner's (1990) form with ECHAM6's
  coefficient at T63 (Heymsfield and Donner 3.29, ECHAM6 3.0 at other
  resolutions), implicitly from the top down within the step: each layer
  keeps 1/(1 + v Δt/Δz) of its ice with what fell into it, the layer below
  takes the rest as ice in its ice share and as precipitation in its
  liquid share (at any temperature on the phase ramp, not melting),
  and the column is adjusted again so that ice falling into dry air
  sublimates there; only the liquid share converts over the lifetimes.
  30 mg/kg of ice at 193 hPa falls at 0.39 m/s and keeps 0.817 of itself
  over 600 s. Heymsfield and Donner's 3.29 is the fit to observed cirrus;
  2.5 was taken over it because three N=64 days from ten64_day0183 put
  the upper troposphere drier and its high cloud sparser at 3.29, an
  outcome of the missing anvil source rather than an observation of the
  fall speed.
- Exponential-random overlap (`cloudOverlap` 'exponentialRandom'):
  adjacent cloudy layers overlap with α = exp(−Δz/z₀) between maximum and
  random (Hogan and Illingworth 2000), z₀ = 2899 − 27.59 |latitude°| m
  (Shonk et al. 2010, CloudSat and CALIPSO), separated layers randomly. On
  the diagnosis state it gives a total cover of 0.43 for the radiation's
  0.40.
- `iceNucleation` (off): a clear layer colder than 235.15 K forms ice only
  above min(q_sw, (2.583 − T/207.8) q_si) (Kärcher and Lohmann 2002, as the
  IFS Cy43r1 takes it, §7.2.4c). Three days from ten64_day0183 at the fall
  coefficient 3.29: upper-tropospheric RH_i 0.59 over the warm pool with or
  without it, below the threshold, and the warm pool's high cover 0.131
  without, 0.111 with it.

Screens, day 186 of three N=64 days (SWCRE, LWCRE, W/m²), each change on
top of those before it unless marked: from eight64_day0183, the
saturation adjustment over water with the 1 h lifetimes −36.0, 15.7; ice
saturation alone −34.4, 11.3; falling ice alone −34.3, 14.1; ice
saturation, the uniform distribution and falling ice (3.29) −44.8, 21.0,
without the fall −45.0, 21.0; with exponential-random overlap −50.1, 20.9;
with melting of the falling ice −48.3, 20.1; on the overlap's code with
the convective cloud fraction below (−50.5, 21.3), the plumes' rain rate
10⁻³ m⁻¹ −51.8, 22.0 and the fall coefficient 2.5 −51.7, 23.1. From ten64_day0183 on the code with melting: −55.2,
19.2; with nucleation −54.9, 18.4; RH_c's exponent 4 −59.8, 19.8 (earlier
code, −56.7, 19.8 at 2); RH_c aloft 0.70 −58.0, 20.1 (same code); the fall
coefficient 2.5 without nucleation (the defaults) −56.5, 21.0, with it
−56.0, 20.1.

A Xu and Krueger (1991) convective cloud fraction as CAM3 takes it
(k₁ ln(1 + 500 M), k₁ 0.07 shallow and 0.14 deep) was built and removed:
it spreads only a layer's resolved condensate, and three days from
eight64_day0183 moved the total cover by 0.01 and either effect by
0.4 W/m². With `condensation` 'saturation', `iceSaturation` false and
`iceFall` null the CPU engine reproduces the parent's 12-step digests
(3ca002d1, 3d0c610f, da3ea94c); the defaults pin 04f4251c. One moist step
over every column of day 193 of the ten-day run below keeps c_p T + L q to
1.1·10⁻¹⁵ and water with the precipitation to 8.3·10⁻¹⁶ relative (ice fell
in 31,623 of 40,962 columns; eight64_day0183 1.0·10⁻¹⁵ and 7.7·10⁻¹⁶); on
362 random N=6 columns the GPU keeps them to 3.1·10⁻⁸ and 3.4·10⁻⁸ and
agrees with the CPU to 1.4·10⁻⁴ K and 1.4·10⁻⁷ kg/kg, and on the convection's random columns, under the defaults,
the old switches and with nucleation, to 1.2·10⁻⁴ K, 3.6·10⁻⁸ and
1.7·10⁻⁸ kg/kg; the radiation's layer heating under the uniform cover, over
ice and with either overlap agrees to 1.5·10⁻⁴ against 27 K/day.

The regimes, day 186, before → after (the three-day N=64 runs from
eight64_day0183 and, total only, ten64_day0183):

| day 186 | total, eight64 (ten64) | high (its share below τ 3.6) | low | deck | upper-tropospheric RH_i (layer area above ice saturation) | Earth |
|---|---|---|---|---|---|---|
| global | 0.40 → 0.47 (0.44 → 0.49) | 0.135 (0.43) → 0.216 (0.34) | 0.237 → 0.231 | 0.015 → 0.014 | 0.73 (0.27) → 0.54 (0.01) | 0.65–0.68; high 0.2–0.3, 0.6 of it thin |
| warm pool sea, 10S–10N 120–170E | 0.15 → 0.43 (0.16 → 0.35) | 0.134 (0.65) → 0.341 (0.38) | 0.014 → 0.105 | 0.000 → 0.000 | 0.99 (0.47) → 0.68 (0.01) | 0.80–0.90; high 0.55–0.70 |
| Pacific ITCZ | 0.31 → 0.53 (0.30 → 0.52) | 0.278 (0.82) → 0.421 (0.39) | 0.027 → 0.168 | 0.000 → 0.000 | 1.14 (0.68) → 0.76 (0.01) | 0.70–0.85; high 0.45–0.60 |
| N Pacific trades, 15–25N 170–130W | 0.14 → 0.31 (0.21 → 0.33) | 0.089 (0.80) → 0.143 (0.47) | 0.056 → 0.184 | 0.000 → 0.000 | 0.96 (0.43) → 0.72 (0.00) | 0.35–0.55; low 0.2–0.4 |
| S Pacific trades, 10–20S 160–120W | 0.21 → 0.28 (0.15 → 0.28) | 0.119 (0.46) → 0.162 (0.33) | 0.035 → 0.090 | 0.032 → 0.020 | 0.69 (0.17) → 0.54 (0.00) | as above |
| Atlantic trades, 10–20N 50–25W | 0.04 → 0.26 (0.25 → 0.54) | 0.006 (1.00) → 0.142 (0.66) | 0.023 → 0.130 | 0.005 → 0.002 | 0.94 (0.45) → 0.72 (0.00) | as above |
| SE Pacific | 0.29 → 0.35 (0.55 → 0.62) | 0.035 (0.97) → 0.104 (0.97) | 0.087 → 0.104 | 0.133 → 0.125 | 0.42 (0.04) → 0.42 (0.00) | low 0.6–0.8 |
| Peru | 0.41 → 0.56 (0.50 → 0.50) | 0.048 (0.75) → 0.292 (0.94) | 0.053 → 0.072 | 0.368 → 0.371 | 0.76 (0.26) → 0.74 (0.00) | low 0.6–0.8 |
| Namibia | 0.71 → 0.68 (0.34 → 0.34) | 0.000 (n/a) → 0.043 (0.85) | 0.316 → 0.383 | 0.555 → 0.451 | 0.88 (0.51) → 0.72 (0.00) | low 0.6–0.8 |
| California | 0.13 → 0.36 (0.46 → 0.40) | 0.000 (n/a) → 0.056 (0.81) | 0.095 → 0.282 | 0.000 → 0.000 | 0.77 (0.13) → 0.67 (0.00) | low 0.5–0.7 |
| Southern Ocean, 40–60S | 0.75 → 0.81 (0.77 → 0.81) | 0.197 (0.32) → 0.297 (0.24) | 0.570 → 0.502 | 0.009 → 0.019 | 0.66 (0.23) → 0.42 (0.01) | 0.80–0.90; low 0.5–0.7 |
| N Atlantic, 40–60N 50–10W | 0.50 → 0.57 (0.62 → 0.67) | 0.142 (0.50) → 0.162 (0.22) | 0.345 → 0.389 | 0.011 → 0.013 | 0.56 (0.15) → 0.41 (0.00) | 0.75–0.85 |
| N Pacific, 40–60N 150E–140W | 0.57 → 0.65 (0.77 → 0.77) | 0.132 (0.22) → 0.230 (0.12) | 0.436 → 0.400 | 0.008 → 0.003 | 0.48 (0.12) → 0.38 (0.00) | 0.80–0.90 |
| 60–90N | 0.71 → 0.76 (0.68 → 0.71) | 0.253 (0.27) → 0.359 (0.12) | 0.435 → 0.336 | 0.011 → 0.011 | 0.72 (0.24) → 0.49 (0.00) | 0.80–0.90 |
| 60–90S | 0.81 → 0.39 (0.76 → 0.34) | 0.314 (0.44) → 0.221 (0.18) | 0.443 → 0.145 | 0.019 → 0.015 | 1.08 (0.50) → 0.65 (0.00) | 0.65–0.80 |
| land 60S–60N | 0.29 → 0.41 (0.34 → 0.40) | 0.131 (0.37) → 0.223 (0.36) | 0.127 → 0.144 | 0.000 → 0.000 | 0.72 (0.28) → 0.55 (0.01) | 0.50–0.60 |
| sea 60S–60N | 0.36 → 0.47 (0.41 → 0.51) | 0.106 (0.49) → 0.198 (0.39) | 0.233 → 0.260 | 0.020 → 0.020 | 0.70 (0.25) → 0.52 (0.00) | 0.68–0.75 |

The classes (`scripts/cloudClasses.mjs`, the same states, cover under the
radiation's overlap; SW, LW W/m²; grid path g/m²; in-cloud path; share of
the cover below τ 3.6 / 3.6–23 / above 23):

| eight64 day 186 | before | after |
|---|---|---|
| low | 0.237, −20.5, 3.5; 31.8; 135; 0.14 / 0.65 / 0.22 | 0.231, −22.7, 3.6; 29.4; 127; 0.11 / 0.66 / 0.23 |
| middle | 0.067, −4.3, 2.0; 15.3; 230; 0.27 / 0.37 / 0.36 | 0.063, −6.8, 1.8; 15.3; 243; 0.03 / 0.48 / 0.49 |
| high | 0.135, −8.0, 8.6; 34.5; 255; 0.43 / 0.25 / 0.32 | 0.216, −17.6, 15.3; 54.6; 253; 0.34 / 0.32 / 0.35 |
| cumulus, deck | 0.010, −1.9, 0.2; 0.015, −1.3, 0.1 | 0.013, −2.1, 0.1; 0.014, −1.1, 0.1 |
| all: SWCRE, LWCRE; cover; liquid / ice g/m² | −38.1, 15.1; 0.399; 64.8 / 23.6 | −53.4, 21.8; 0.473; 81.1 / 24.5 |
| by region: 30S–30N; 30–60S; 60–90S (SWCRE, LWCRE, cover) | −25.2, 9.2, 0.221; −80.5, 23.8, 0.627; −22.9, 26.1, 0.812 | −47.9, 22.3, 0.368; −92.4, 27.2, 0.703; −19.9, 15.2, 0.394 |
| ten64 day 186: all; liquid / ice | −45.9, 16.4; 0.436; 85.9 / 28.6 | −58.6, 20.2; 0.489; 102.6 / 30.1 |

The runs (N=64 and N=128 GPU, `everySteps` 8; day means; Earth: albedo
about 0.29, SWCRE −47 ± 4, LWCRE +26 ± 3, rain 2.6–2.8 mm/d):

| | before | after |
|---|---|---|
| eight64 day 186: albedo; ASR − OLR; SWCRE; LWCRE; rain day mean (audit) | 0.253; +10.2; −36.0; 15.7; 1.78 (2.12) | 0.293; −0.8; −49.7; 21.9; 1.84 (2.19) |
| its audit: ITCZ rain, mm/d (convective share); global convective share | 5.05 (0.99); 0.45 | 3.97 (0.93); 0.33 |
| SE Pacific / California low cloud, radiative; SE Pacific rain | 0.239 / 0.181; 0.49 | 0.258 / 0.335; 0.43 |
| ten64 day 186: albedo; ASR − OLR; SWCRE; LWCRE; rain (audit) | 0.278; +4.4; −44.3; 17.1; 2.72 (2.67) | 0.314; −7.9; −56.5; 21.0; 2.83 (2.77) |
| its audit: ITCZ rain (convective share); SE Pacific low cloud; zonal rain peak | 3.56 (0.93); 0.534; 11.2 at 7.5N | 4.84 (0.62); 0.574; 12.9 at 9.5N |
| eight128 day 186: albedo; ASR − OLR; SWCRE; LWCRE; rain | 0.243; +12.1; −32.5; 15.1; 2.00 | 0.282; 0.0; −45.8; 19.7; 2.19 |
| its cover: global; warm pool; ITCZ | 0.38; 0.38; 0.16 | 0.44; 0.56; 0.34 |
| 60–90N ice loss from nine64_day0091, 10³ km³/day | 0.310 (9.191 → 8.261) | 0.306 (→ 8.272) |
| its day-94 60–90N SWCRE; LWCRE; cover | −96.5; 22.8; 0.676 | −109.0; 17.8; 0.734 |
| fresh atlas start days 6–10: albedo; ASR − OLR; SWCRE; LWCRE; rain | 0.368; −10.1; −73.8; 32.5; 4.56 | 0.396; −19.3; −82.7; 36.8; 4.71 |

Ten days at N=64 from eight64_day0183, day by day 184–193 (the
upper-tropospheric RH_i of 150–350 hPa, global and over the warm pool;
high cover; resolved liquid / ice path; from each day's state):

| | before | after |
|---|---|---|
| ASR − OLR | +15.7, +11.3, +10.1, +10.5, +10.6, +9.6, +7.6, +7.1, +8.1, +8.2 | +14.6, +3.2, −0.8, −2.3, −2.8, −4.7, −6.5, −7.1, −7.5, −7.2 |
| SWCRE; LWCRE | −28.0 … −38.8; 16.0 … 16.5 | −42.8, −45.2, −49.7, −51.7, −51.2, −51.5, −51.9, −52.0, −51.6, −51.1; 31.8, 22.2, 21.9, 22.5, 21.8, 21.0, 20.0, 19.6, 19.2, 19.1 |
| UT RH_i global; warm pool | 0.70 → 0.77; 0.90 → 0.96 (1.03 on day 188) | 0.55, 0.54, 0.54, 0.53, 0.52, 0.51, 0.50, 0.50, 0.49, 0.49; 0.69 … 0.56 |
| high cover; its grid path, g/m² | 0.156 → 0.139; 33 → 38 | 0.242, 0.218, 0.216, 0.208, 0.197, 0.182, 0.176, 0.170, 0.162, 0.157; 48 → 57 (67 on day 188) |
| liquid / ice path, g/m² | 49.5 / 21.5 → 66.0 / 25.0 | 56.6 / 22.7 → 84.5 / 26.6 (86.8 / 29.1 on day 188) |
| rain, mm/d | 0.73 … 2.57 | 0.87, 1.19, 1.84, 2.34, 2.68, 2.80, 2.74, 2.80, 2.66, 2.63 |

The upper troposphere dries over the ten days and its high cloud thins
with it: the ice that forms above the critical humidity falls out, and the
plumes, which top at 230–270 hPa over the warm pool, do not resupply it.
No run made a NaN or clamped the ocean but nine64's known day-94 step;
five64_day2281 and twin64_day0900 (27 layers), six64_day1004,
seven64_day0639, nine64_day0091 and m21b64_day0183 load and take two CPU
steps with nothing non-finite. Pace (`js/gpu/profile.module.js`, 128 steps
after 64 from eight64_day0183 and eight128_day0183, alternated twice with
the parent alone on the GPU): step median 20.2 and 20.2 ms against 19.7
and 19.7 at N=64 (+2.5 %), 89.1 and 89.2 against 87.4 and 87.4 at N=128
(+1.9 %); the adjust pass 16.5 against 15.2 ms and the physics pass 10.2
against 9.7 at N=128. The full suite (50 files, concurrently) passes,
parallel.test.mjs alone.

What still misses: the total cover (0.47 against 0.65–0.68) and its spread:
the warm pool 0.43 and the ITCZ 0.53 against 0.7–0.9, with their high
cloud too thick (38–39 % of it thin against ISCCP's 60 %) and their upper
troposphere drying through the ten days; the anvil, which a diagnostic
cover cannot hold, wants the detrained condensate and cloud-fraction
sources of Tiedtke (1993) (IFS Cy43r1: source (1 − a) D_up, erosion
a K (q_s − q), K 3·10⁻⁶ s⁻¹), which need an advected cloud fraction; the
trades' low cover (0.09–0.18 against 0.2–0.4) under a cumulus layer at
0.8–0.85 humidity below its RH_c of 0.90–0.93; 60–90S at 0.39 against
0.65–0.80; the storm tracks' deep cloud (in-cloud 300–460 g/m², 12–24 %
thin); the deck boxes, unchanged. Global outcomes, eight64 day 186: SWCRE
−49.7 and LWCRE 21.9 against −47 and +26, albedo 0.293, ASR − OLR −0.8
(ten64: −56.5, 21.0, 0.314, −7.9). A later sweep may vary RH_c aloft
within 0.70–0.80 (CAM3, ECHAM6, IFS) and at the surface within 0.95–0.994
(ECHAM6 across resolutions), the profile's exponent about ECHAM6's 2 (no
published range found), the fall
coefficient within 2.5–3.29 (ECHAM6, Heymsfield and Donner) and its
exponent at 0.16, and the decorrelation length by a factor of 0.5–1.5
about Shonk et al.'s (Hogan and Illingworth 2000 found 1.6 km, Barker 2008
about 2 km); not the phase ramp, the ice saturation or the optics.

Review of the cloud amount and falling ice (Oct 1). The three N=64 days
from eight64_day0183 rerun to the same day lines (day 186: albedo 0.293,
ASR − OLR −0.8, SWCRE −49.7, LWCRE 21.9), the regime and class tables to
the values above, and the Arctic three days to 9.191 → 8.272·10³ km³
(0.306 a day; no ice lies south of 60N on day 91, so 50–90N and 60–90N
agree). The moist step alone, 256 steps of 337.5 s on the day-186 state,
keeps every column's c_p T + L q and water with the precipitation to
3.7·10⁻¹⁵ and 3.1·10⁻¹⁴ relative on the CPU; on the GPU (fusion heat of
snow off, it is released into the lowest layer by design) to 1.2·10⁻⁷ and
1.6·10⁻⁷ after one step and 2.7·10⁻⁵ and 1.8·10⁻⁵ after the day, the
global water changing by 3.3·10⁻⁷ and the enthalpy by 582 J/m² (0.007 W/m²),
f32 round-off accumulated (with the three switches off 8.6·10⁻⁶ and
1.2·10⁻⁵); both engines precipitate 0.175 kg/m² over the day. One layer of
20 mg/kg of ice at 266 hPa (236.4 K, cover 0.453, Δz 1122 m) falls at
0.43 m/s; its mass-weighted fall over 10, 337.5 and 3600 s is 4.25, 127.7
and 652.7 m against v Δt/(1 + v Δt/Δz) of the ice share 0.967 plus the
liquid share's conversion over the 1 h lifetime, 4.25, 127.7 and 652.7 m. Ice of 10⁻⁹ to
3·10⁻³ kg/kg above 250 K stays non-negative and conserved at steps up to
10⁶ s. Ice falling into a layer below sublimates only while it is under
ice saturation; under 'saturation' adjustment one step leaves a mixed-phase
layer at most 1.001 of saturation (at 236 K the slope omits
(e_w − e_i) dα/dT), corrected by the next step; a slope with that term
moved the gpuModel rain-split parity (large-scale rms 5.7·10⁻⁴ to 6.3·10⁻³
on 10⁻⁶ kg/m² of onset drizzle) and was not kept. The uniform cover equals
Sundqvist's 1 − sqrt((1 − RH)/(1 − RH_c)) of the adjusted grid humidity to
0.0004 at 314 hPa and 0.004 at 510 hPa for RH_t 0.80–0.99, is continuous
and bounded, and the radiation's column cover is the decomposition's (0.473
and 0.394 for the globe and 60–90S in both). The in-cloud path the
shortwave takes, the column's path over its visible cover, has a median of
102 g/m², 99 % below 1.48 kg/m², and 67 of 40,962 columns above
10 kg/m² (largest 2.6·10³ kg/m²): plume cumulus of cover below 10⁻⁶ with
no cover floor and a visibility weight of the grid-mean path, over column
covers of 10⁻¹⁰, which carry no flux. Step cost under the exclusive lock,
two alternations: 20.21 and 20.20 ms against 19.88 and 19.74 at N=64
(+1.7 and +2.3 %), 89.65 and 89.46 against 87.50 and 87.57 at N=128
(+2.5 and +2.2 %). five64_day2281, twin64_day0900 (27 layers),
six64_day1004, seven64_day0639, m21b64_day0183 and ten64_day0183 load and
take two CPU steps with nothing non-finite.

**Integration (Oct 1).** The surface by class (M21, branch clear-sky) and
the cloud optics, uniform cover, falling ice and overlap above (branch
cloud-optics) merged on sweep2 (e07e39a), both sides' physics kept. PH
gains only the surface's SNOWALB, CANOPY, SEASONL and SEASONW; the GPU's
upward-absorption escapes take the column's own cloud depth (phase or
gray, the deck's included), as the CPU's do through the same shortwave
call. With grey ice, gray optics and the saturation adjustment the CPU
digests are both parents' (b892e42f, e228ab4c, d8b73e96 and 3ca002d1,
3d0c610f, da3ea94c); the defaults' twelve-step digest is 97902663. Two
cover rules fixed on both engines (`test/cloudIce.test.mjs`,
`test/gpuModel.test.mjs`): a plume layer is seen at least as its plume
(fraction times the visibility of the plume's own path), so a cumulus of
cover 10⁻⁸ and a 262 g/m² plume takes 262 g/m² in cloud where it took
10⁵ kg/m², and one of 0.005 and 52 g/m² takes 52 where it took 227; under
the uniform condensation the stratiform blend takes a (q_t − q_s(T_l)) at
the cover's saturation (the ice ramp), its half-width at most
a (1 − RH_c) q_s(T_l): a layer at 434 hPa and 260 K at 1.02 of the
ice-ramp saturation is covered 0.639 under the full blend where the
liquid saturation at T gave the floor 0.01. That second rule had emptied
60–90S: the day-186 state of the eight64 run below, diagnosed under the
code before it, covers 0.42 there (global 0.47), under it 0.68 (0.50).
The full suite (52 files, concurrently) passes. `scripts/sweep/runs.mjs`
drops `cloudScattering` and `cloudAbsorption` from PARAMETERS2 (their base
was NaN; set, they turn the optics gray); TERMS2 holds no clear-sky albedo
term.

Runs (N=64 and N=128 GPU, bl34, `everySteps` 8, from copies of the
states; the last day's means; cover, audit and classes on the last state;
"surface" is the clear-sky branch under gray cloud and "cloud" the
cloud-optics branch, both from their own paragraphs; Earth: albedo about
0.29, SWCRE −47 ± 4 and LWCRE +26 ± 3 W/m² (CERES EBAF), rain 2.6–2.8
mm/d (GPCP), cover 0.65–0.68 (ISCCP, MODIS, CALIPSO)):

| | albedo | ASR; OLR; ASR − OLR | SWCRE; LWCRE | clear-sky albedo | rain, day (audit); convective share | total cover; 60–90S |
|---|---|---|---|---|---|---|
| eight64 + 3 d (day 186) | 0.296 | 239.7; 240.8; −1.1 | −49.9; 22.8 | 0.149 | 1.83 (2.19); 0.33 | 0.50; 0.68 |
| surface alone | 0.312 | 234.2; 242.2; −8.0 | −55.5; 17.6 | 0.149 | (2.05) | |
| cloud alone | 0.293 | −0.8 | −49.7; 21.9 | | 1.84 (2.19); 0.33 | 0.47; 0.39 |
| ten64 + 3 d (day 186) | 0.318 | 232.4; 240.5; −8.1 | −57.2; 22.1 | 0.150 | 2.83 (2.78); 0.20 | 0.52; 0.64 |
| cloud alone | 0.314 | −7.9 | −56.5; 21.0 | | 2.83 (2.77); 0.20 | 0.49; 0.34 |
| nine64 + 3 d (day 94, June) | 0.304 | 237.1; 241.3; −4.2 | −51.3; 23.3 | 0.153 | 2.42 (2.52); 0.30 | 0.50; 0.66 |
| surface alone | 0.326 | 229.5; 240.6; −11.1 | −58.7 | 0.154 | | |
| nine64_day0365 + 3 d (day 368, March) | 0.312 | 234.3; 229.3; +5.0 | −51.6; 28.6 | 0.160 | 2.45 (2.66); 0.35 | 0.55; 0.70 |
| surface alone | 0.328 | 228.8; 233.3; −4.5 | −57.2 | 0.160 | | |
| eight128 + 3 d (day 186) | 0.285 | 243.5; 243.6; −0.1 | −46.2; 20.8 | 0.149 | 2.19 (2.46); 0.22 | 0.47; 0.64 |
| cloud alone | 0.282 | 0.0 | −45.8; 19.7 | | 2.19 | 0.44 |

The audit (`scripts/verticalAudit.mjs`) on the same states, eight64,
ten64, June, March and eight128: the Pacific ITCZ's firing-column heating
peaks at 975, 974, 975, 517 and 439 hPa (Earth 400–500), its rain 4.02,
4.89, 6.65, 3.21 and 2.74 mm/d (6–9); radiative low cloud over the SE
Pacific 0.26, 0.56, 0.37, 0.20 and 0.28, over Peru 0.39, 0.38, 0.05, 0.03
and 0.34 (0.6–0.7); global evaporation 2.13–2.65 mm/d. Regimes
(`scripts/cloudRegimes.mjs`) on eight64 day 186: warm pool 0.43, ITCZ
0.53, high cover 0.232 (0.36 of it below τ 3.6), upper-tropospheric RH_i
0.54; classes (`scripts/cloudClasses.mjs`) there SWCRE −53.7 and LWCRE
22.7, cumulus cover 0.014 at 335 g/m² in cloud (0.013 in the cloud
branch), liquid / ice 81.0 / 24.8 g/m². The surface classes
(`scripts/clearSkyBudget.mjs`) are M21's to 0.01 on every state: the open
sea matches at every latitude but 0–30° in June (0.102 against
0.08–0.10), partly vegetated land matches (0.199–0.218), dense vegetation
reads 0.156–0.166 (0.12–0.15), cold snow on open land 0.71–0.81
(0.80–0.85; matches in March), wet snow on sea ice 0.81–0.82
(0.65–0.75), melting bare sea ice 0.481 (0.45–0.55), cold snow on sea ice
0.83–0.84 and the ice sheets 0.80 (both match).

Ten days from eight64_day0183, days 184–193: ASR − OLR +14.8, +3.0, −1.1,
−2.6, −3.1, −4.9, −7.1, −7.7, −8.3, −7.5 W/m²; LWCRE 33.2, 23.1, 22.8,
23.3, 22.6, 21.7, 21.0, 20.6, 20.1, 20.0; SWCRE −43.2 to −52.8 (−51.5 on
day 193); the upper troposphere's RH_i 0.55, 0.54, 0.54, 0.53, 0.52,
0.51, 0.50, 0.50, 0.49, 0.49 (warm pool 0.69 → 0.56); high cover 0.264,
0.235, 0.232, 0.223, 0.208, 0.196, 0.193, 0.185, 0.179, 0.175; rain 0.86
→ 2.69 mm/d: the cloud branch's drift (−7.2 W/m² by day 193).

The Arctic test against Earth. The three-day loss from nine64_day0091
(the June solstice) had the target 0.15–0.18·10³ km³/day, set under the
gray cloud. Earth: PIOMAS v2.1 monthly volumes (Schweiger et al. 2011,
JGR 116, C00D06; Polar Science Center), 2011–2025 means June 16.81 and
July 10.07·10³ km³: mid-June to mid-July 0.221·10³ km³/day, 1.6 % of
the standing volume a day near the end of June, about 1.2 % at the
solstice (mid-May to mid-June 0.152, 0.8 %); 1979–1988 0.204·10³
km³/day, 0.8 %. SHEBA, the 20-day block about the solstice (Intrieri et
al. 2002, JGR 107(C10), 8039, read from their Figures 3 and 5–10):
downwelling sunlight about 300 W/m², net sunlight 85 (100 in the next
block), downwelling longwave about 280, net longwave −40 (−15 next),
sensible heat near 0 and latent about 5 upward, a net surface gain of
about 40 W/m² (85 next); tower albedo 0.70 and survey-line albedo
0.5–0.55, the line a mix of melting snow at 0.7, bare ice and new ponds at
0.3 by 15 June and 0.4 on the mean by the end of July (Perovich et al.
2002, JGR 107(C10), 8044); cloud cover 0.77; net surface cloud forcing
−10 W/m² then, −49 in early July. Ocean heat flux into the ice a few
W/m² through June, 16.8 W/m² for July and 33 at its peak (Perovich and
Elder 2002, GRL 29). Over the Arctic Ocean clouds warm the surface on the
annual mean (+10 W/m²) and cool the top of the atmosphere (−12; Kay and
L'Ecuyer 2013, JGR 118, 7219).

The model over three and ten days from nine64_day0091 (a GPU rerun with
the run's options; volumes Σ A·h over the sea cells north of 60N, 9.191
at the start as above; fluxes the means over every step over the cells
iced at the start, leads included, positive into the surface):

| | 70–80N | 80–90N | 60–90N | 60–90N, ten days |
|---|---|---|---|---|
| volume, 10³ km³: start → end | 5.057 → 4.489 | 4.134 → 3.759 | 9.191 → 8.248 | → 5.879 |
| loss a day, 10³ km³; share of the start | 0.189; 3.7 % | 0.125; 3.0 % | 0.314; 3.4 % | 0.331; 3.6 % |
| area, 10⁶ km²; mean thickness, m | 3.81 → 3.54; 1.33 → 1.27 | 2.84 → 2.70; 1.46 → 1.39 | 6.65 → 6.24; 1.38 → 1.32 | → 5.10; → 1.15 |
| snow on the ice, kg/m² | 0.1 | 0.0 | 0.1 | 0.0 |
| downwelling sunlight; surface albedo; absorbed | 190.5; 0.372; 119.6 | 235.7; 0.404; 140.4 | 208.7; 0.387; 128.0 | 201.4; 0.360; 129.0 |
| downwelling longwave; net longwave | 306.3; −7.0 | 298.0; −16.1 | 303.0; −10.7 | 303.4; −9.7 |
| sensible; latent | +17.9; −1.9 | +7.0; −1.6 | +13.5; −1.8 | +21.8; −2.1 |
| net surface flux; ocean flux into the ice | 128.6; 7.2 | 129.7; 1.0 | 129.0; 4.7 | 138.9; 5.5 |
| melt they imply, 10³ km³/day | 0.194 | 0.126 | 0.320 | 0.346 |

No ice lies at 60–70N on day 91 (Earth's lies in Hudson and Baffin Bays
and on the Barents and Bering margins); the iced cells are 0.76 covered,
and the pack is 6.6·10⁶ km² at 1.38 m where PIOMAS holds about 16·10³
km³ at the solstice. The surface and ocean fluxes account for the loss
(0.320 against 0.314). Over the ice itself (the leads' 24 % absorbing
about 0.93 of their sunlight) the absorbed sunlight is about 107 W/m² at
an albedo of 0.49 (melting bare ice 0.48 with the thickness ramp), against
SHEBA's 85 at 0.55–0.70, under downwelling sunlight of 209 against 300
(the June 60–90N SWCRE is −105 W/m², classes on day 94). The ice's net
gain is about 108 W/m² against SHEBA's 40: +22 of the difference
sunlight (no snow left at the solstice, so no melting-snow albedo near
0.7), +30 longwave (net −11 against −40; the downwelling 303 against about
280 makes 23 of it), +14 sensible heat
from air warmer than the melting surface (SHEBA near 0), +3 latent; the
ocean's 4.7 W/m² is SHEBA's. The model's 3.4 % a day is 2–3 times
Earth's 1.2–1.6 %: the surface flux makes about 2 of it (4.7 cm of ice a
day over the pack, against PIOMAS's 0.15–0.22·10³ km³ a day over the June
ice area of the NSIDC Sea Ice Index v4, 8.53·10⁶ km² for 2011–2025
(extent 10.78), 1.8–2.6 cm and about 2.3 at the solstice, and SHEBA's
1.1 cm at 40 W/m², 2.4 at 85), the thin pack about 1.4 (1.38 m against
PIOMAS's June 16.81·10³ km³ over that area, 1.97 m), the ocean none. The
test's reference is the fractional loss: 1.2–1.6 % of the volume a day at
the June solstice (PIOMAS v2.1, 2011–2025 monthly means: the rate
interpolated to 21 June, 0.20·10³ km³/day, over the volume then, 15.8,
and mid-June to mid-July over its mean volume), 0.11–0.15·10³ km³/day for
this state's 9.19, which the second sweep's score now takes (0.129 ±
0.018 in `scripts/sweep/score.mjs`); the melt per area of the pack, about
2.3 cm of ice a day (PIOMAS over the NSIDC area) and 1.1–2.4 cm (SHEBA),
is the second check.

Cost under the exclusive lock (`js/gpu/profile.module.js`, 128 steps
after 16 from eight64_day0183 and eight128_day0183, alternated twice with
the pre-merge tree 15754a7 built by `git archive`): step median 20.49 and
20.43 against 19.59 and 19.56 ms at N=64 (+4.6 and +4.4 %), 90.03 and
90.16 against 86.94 and 86.79 at N=128 (+3.6 and +3.9 %); the physics and
boundary-layer passes 2.79 against 2.25 ms and 11.70 against 9.55 ms. A
day at N=128 takes 60 s of wall time (one day from eight128_day0183 after
6 s of setup). No run made a NaN; nine64's known day-94 clamp is the
only one.

Review of the integration (Oct 1). Each conflicted file's merge equals
the union of both sides' changes but for the GPU escapes above; PH's slots
are named in one sequence, so no two fields share one. The suite (52
files, concurrently) passes. On real states, one step (N=64; eight64_day0183
and nine64_day0091): the CPU with gray optics, maximum-random overlap and
the saturation adjustment is bit-identical to the clear-sky branch in
every flux and field; the clear-sky ASR and OLR of every cell are
bit-identical on each engine under the defaults, the gray optics, either
overlap and the old condensation; `scripts/clearSkyBudget.mjs` gives the
clear-sky branch's table to the digit, with and without the gray cloud.
The shortwave closes per column to 4·10⁻¹⁶ of the beam on the CPU and
2.5·10⁻⁷ on the GPU (also on five64_day2190, nine64_day0365 and
eight128_day0183, which load on both engines); the moist adjustment keeps
column water and c_p T + L q to 8·10⁻¹⁶ on the CPU and 1.6·10⁻⁶ and
8.5·10⁻⁷ on the GPU over 16 steps. Engines under the defaults, one step:
ASR within 0.03, OLR 0.004 and the net surface flux 0.04 W/m² in every
cell (the clear-sky branch's engines differ by up to 165 W/m² in about
400 cells, as do 15754a7's: the saturation adjustment's cover; the new
defaults do not carry it); nine to fifteen columns differ by more than
0.1 K in their lowest ten layers, as on 15754a7. A GPU run repeats
itself bit for bit, but a restart from a saved state is not bit-exact:
three days in one segment and three one-day segments from eight64_day0183
part in every field by day 186 (global day means alike to the log's
digits), on both parents too. The eight64 rerun, day by day as above,
matches the table to every digit; the June rerun, day by day where the
table's ran in one segment, gives ASR 236.9, ASR − OLR −4.4 and SWCRE
−51.4 on day 94 and the same loss, 0.314·10³ km³/day (9.191 → 8.249).

**The spectral gases with the cloud and the surface (Oct 1).** Branch
gas-benchmark (1f61d5c, M21's "The gases on fixed profiles") merged on
sweep2 (5985713), both sides kept. PH gains LWSFCSUM after SNOWALB,
CANOPY, SEASONL and SEASONW in one sequence, LV gains OZS; the spin-up's
daily line carries the sea's net surface longwave and the iced cells'
sunlight. Where the gases meet cloud the gas branch already went through
the arrays the phase optics now fill, so on the CPU nothing else changed:
each g-point combines the layer's cloud emissivity f (1 − exp(−κ W/f)), κ
by phase, as 1 − (1 − ε_gas)(1 − ε_cloud) (in the longwave every layer's
cover enters through its emissivity, under either gas scheme); the CLIRAD
vapour, O₂ and CO₂ take the beam along the whole column's path before the
cloud's two-stream, cloud or none, as the Lacis–Hansen vapour did; the
light the surface sends up loses the near infrared the gases' path down
plus 5/3 of the column absorbs beyond the first, times the all-sky escape
of the exponential-random, deck-blended streams; the deck's light is what
the gases and aerosol leave. On the GPU two conflicted lines needed both
sides: the deck's light and the upward escapes take the visible light the
gases leave (visibleTaken) and the column's own cloud depth (phase or
gray). Fixed on both engines: a dry-air share 1 − q below zero (q above 1
at a 110 Pa layer of `test/cloudOptics.test.mjs`) made the O₂ and CO₂
square roots and the longwave's well-mixed paths NaN; it is held at 0, the
same bits wherever q < 1. Digests (CPU, 12 steps at N=4): under the gray
gases the defaults give the integration's 97902663 and the overcast gray
cloud on grey ice 3ca002d1; under the spectral gases that cloud gives the
gas branch's 41e2f59e, the defaults fbd9be1f.

`test/gasCloud.test.mjs`: under transparent spectral gases a 100 g/m²
liquid layer at 285 K reflects 0.64304 and 10 g/m² emits 0.77687, 20 g/m²
of ice at 213.15 K 0.41806 and 0.89985, the hand values of
`test/cloudOptics.test.mjs` (to 10⁻¹² and 10⁻⁹), as under the gray gases
with an open window; O₂ takes 0.63 % of the beam. In a column at 60 %
relative humidity below 300 hPa the OLR with the cloud is an independent
g-point sum with the hand emissivity in the cloud's layer to 10⁻¹², the
cloud reflects its hand share of the light the gases leave, and the gases
take the same light with and without it. What the gases change is the
overlap: LWCRE of the liquid layer at 800 hPa 5.73 spectral, 5.93 gray,
26.95 for the cloud alone; of the ice at 250 hPa 125.4, 136.4 and 255.6
W/m²; SWCRE −370.8 against −378.3 and −251.5 against −254.7, the light
reaching the cloud 576.6 against 588.2 and 601.6 against 609.2 of 681.0
W/m². With a 100 g/m² layer near 282 K and a 20 g/m² layer near 222 K in
every N=6 column and the defaults (the empirical deck in place of the
mixed layer's), the engines' layer heating differs by 1.4·10⁻³ K/day of
1256, the fluxes at the top and surface and both cloud effects by at most
10⁻³ W/m², and the gases' change of them alike. On eight64_day0183 and
nine64_day0091 (one N=64 step) absorbed plus reflected is the beam to
4.6·10⁻¹³ W/m² on the CPU and 1.6·10⁻⁷ of it on the GPU, the layers'
shortwave heating the atmosphere's absorption to 3.4·10⁻¹³, their
longwave heating σT_s⁴ less the downward longwave less the OLR to
5.1·10⁻¹³ (GPU 5.6·10⁻⁴); CPU against GPU per cell, ASR rms 1.6·10⁻⁷ and
1.3·10⁻⁶ of the field, OLR 6.8·10⁻⁷ and 7.5·10⁻⁷; five64_day2190,
seven64_day0365 and m21a64_day0365 load and close alike. On the
midlatitude-summer profile with 100 g/m² of liquid at 285 K (790 hPa)
under 20 g/m² of ice at 220 K (195 hPa) in every N=4 column, no ozone,
the change of the OLR and of the downward longwave is the sum over
g-points of each layer's phase emissivity acting on the clear-sky fluxes
through the gas transmittances: LWCRE 129.4 at the top and 63.7 at the
surface overcast, 128.2 and 39.0 under the uniform cover (the liquid
layer's 0.521), the CPU to 10⁻¹³ and the GPU to 7·10⁻⁵ W/m²
(`test/gasCloud.test.mjs`).
`scripts/radiationBenchmark.mjs` prints the gas branch's table to the
character. Two deck parity tests (`test/frameGpu.test.mjs`,
`test/gpuModel.test.mjs`) run the gray gases: under the spectral ones a
shallow cumulus fires on the fifth step in one of 362 columns on the CPU
alone, and two steps part one layer's condensate above the deck by
2.9·10⁻⁷ kg/kg, where a 10⁻⁷ perturbation of the state moves the CPU's own
cloud field by 9.3·10⁻³ kg/m² under either scheme. The suite (54 files,
concurrently) passes.

The clear-sky classes (`scripts/clearSkyBudget.mjs`) on the four N=64
day-186, 94 and 368 states below, spectral against gray gases on the same
state (gray: on eight64 the integration's table to the digit): the surface albedos
and their verdicts do not move; the clear-sky reflection at the top falls
by 0.004 over open sea at 0–30°, 0.005–0.007 at 30–50°, 0.009–0.011 at
50–70°, 0.009–0.036 at 70–90°, 0.011–0.012 over partly vegetated land,
0.010–0.012 over snow under forest, 0.019–0.036 among sparse trees,
0.022–0.033 over thin sea ice, 0.028–0.039 over bare ice, 0.035 over cold
snow on open land, 0.055–0.066 over cold snow on sea ice (0.134 on June's
low sun) and 0.039–0.061 over the ice sheets: the near-infrared vapour
absorbs on the way down and up. The open sea stays in CERES's clear-sky
ranges at 0–30°, 30–50° and 50–70° on every state (June's 0–30° moves
from 0.102, high by 0.002, to 0.098). eight64 day 186, global, W/m²:
reflected 50.7 → 48.0, ozone 10.2 → 10.8, vapour 48.8 → 58.5, O₂ and CO₂
0 → 3.5, aerosol 1.4 → 1.3, the surface 229.4 → 218.4.

Runs (GPU, bl34, `everySteps` 8, three days from copies of the states, the
last day's means; rain from the log and, in brackets, the audit's; cover
and 60–90S from `scripts/cloudRegimes.mjs` on the last state; integration
the "Integration (Oct 1)" rows above; Earth as there, the atmosphere's
absorption about 80 W/m² all-sky, Wild et al. 2019):

| | albedo | ASR; OLR; ASR − OLR | atmosphere SW | SWCRE; LWCRE | clear-sky albedo | rain (audit); convective share | cover; 60–90S | sea surface SW; net LW |
|---|---|---|---|---|---|---|---|---|
| eight64 + 3 d | 0.298 | 238.9; 233.2; +5.7 | 82.7 | −53.5; 26.6 | 0.141 | 1.68 (2.12); 0.14 | 0.54; 0.68 | 168.6; −44.3 |
| integration | 0.296 | 239.7; 240.8; −1.1 | 68.6 | −49.9; 22.8 | 0.149 | 1.83 (2.19); 0.33 | 0.50; 0.68 | 184.0; |
| gas branch | | 235.1; 236.9; −1.8 | 80.3 | −58.7; 20.0 | 0.137 | 1.46 | | 166.6; −47.6 |
| ten64 + 3 d | 0.317 | 232.7; 231.9; +0.8 | 83.4 | −59.6; 26.8 | 0.142 | 2.57 (2.54); 0.12 | 0.56; 0.63 | 157.9; −44.2 |
| integration | 0.318 | 232.4; 240.5; −8.1 | 70.1 | −57.2; 22.1 | 0.150 | 2.83 (2.78); 0.20 | 0.52; 0.64 | 172.7; |
| nine64 + 3 d (June) | 0.309 | 235.2; 234.3; +0.9 | 85.4 | −56.0; 26.4 | 0.145 | 2.21 (2.32); 0.17 | 0.55; 0.65 | 148.1; −44.7 |
| integration | 0.304 | 237.1; 241.3; −4.2 | 71.1 | −51.3; 23.3 | 0.153 | 2.42 (2.52); 0.30 | 0.50; 0.66 | |
| nine64_day0365 + 3 d (March) | 0.324 | 230.3; 220.2; +10.1 | 81.8 | −58.5; 33.5 | 0.152 | 2.28 (2.46); 0.16 | 0.60; 0.70 | 161.6; −48.4 |
| integration | 0.312 | 234.3; 229.3; +5.0 | 67.5 | −51.6; 28.6 | 0.160 | 2.45 (2.66); 0.35 | 0.55; 0.70 | 181.8; |
| eight128 + 3 d | 0.282 | 244.6; 236.0; +8.6 | 82.4 | −47.7; 24.8 | 0.141 | 2.06 (2.36); 0.07 | 0.51; 0.62 | 175.2; −48.3 |
| integration | 0.285 | 243.5; 243.6; −0.1 | 68.8 | −46.2; 20.8 | 0.149 | 2.19 (2.46); 0.22 | 0.47; 0.64 | 188.4; |

The audit on the same five states: the Pacific ITCZ's firing-column
heating peaks at 975, 974, 975, 977 and 974 hPa (integration 975, 974,
975, 517, 439; Earth 400–500), its rain 3.23, 6.15, 8.38, 3.43 and 0.91
mm/d (4.02, 4.89, 6.65, 3.21, 2.74; Earth 6–9); radiative low cloud over
the SE Pacific 0.33, 0.61, 0.44, 0.30 and 0.34 (0.26, 0.56, 0.37, 0.20,
0.28), over Peru 0.42, 0.52, 0.31, 0.20 and 0.42 (0.39, 0.38, 0.05,
0.03, 0.34; Earth 0.6–0.7); global evaporation 2.30–2.71 mm/d. Regimes
on eight64 day 186: warm pool 0.52, ITCZ 0.50 (integration 0.43, 0.53),
high cover 0.266 (0.232), upper-tropospheric RH_i 0.59 (0.54).

Ten days from eight64_day0183, days 184–193 (one segment; a GPU rerun of
the same days through `scripts/spinup.mjs` and a driver summing the sea's
fluxes agree to the log's digits): ASR − OLR +22.4, +11.8, +5.6, +4.5,
+3.4, +1.9, +1.0, +1.3, +0.4, +0.9 W/m² (integration +14.8 to −7.5); rain
0.81, 1.03, 1.68, 2.28, 2.54, 2.52, 2.63, 2.61, 2.60, 2.58 mm/d; SWCRE
−40.6 to −57.1, LWCRE 34.5 to 26.5; the global mean surface temperature
17.06, 17.02, 16.93, 16.81, 16.70, 16.63, 16.65, 16.65, 16.64, 16.62 °C
(integration 16.86 to 15.91); the net surface flux over the sea (ice
included) 51.4, 51.2, 40.2, 33.0, 26.9, 23.2, 21.7, 23.5, 23.4, 23.9 W/m²
into it, on day 193 shortwave 163.7, net longwave −43.2, sensible −8.5
and latent −88.1. From day-by-day segments of the same run: the upper
troposphere's RH_i 0.57, 0.58, 0.59, 0.59, 0.59, 0.59, 0.59, 0.58, 0.58,
0.58 (integration 0.55 → 0.49), the warm pool's 0.73, 0.70, 0.71, 0.73,
0.72, 0.70, 0.68, 0.70, 0.73, 0.73 (0.69 → 0.56); high cover 0.281,
0.265, 0.266, 0.264, 0.257, 0.255, 0.251, 0.252, 0.250, 0.251 (0.264 →
0.175); total cover 0.48 → 0.54. The clear-sky OLR (OLR + LWCRE) 263.8,
261.6, 259.8, 259.2, 259.4, 259.7, 259.9, 260.3, 260.5 and 260.7 W/m²
(integration 263.6 on day 186; CERES EBAF about 266, Loeb et al. 2018),
where the gases on the four standard atmospheres give RRTMG's clear-sky
OLR to 0.4 W/m² (M21) and eight64_day0183 itself gives 264.5 at its
first step.

The Arctic three days from nine64_day0091 (the scratchpad diagnostic,
means over every step over the cells iced at the start, leads included,
positive into the surface; integration in brackets): 60–90N 9.191 → 8.272
·10³ km³, 0.306·10³ km³/day, 3.3 % of the standing volume a day (0.314,
3.4 %; PIOMAS 1.2–1.6 %); downwelling sunlight 199.7 (208.7) at a surface
albedo 0.387 (0.387), absorbed 122.4 (128.0; SHEBA 85); downwelling
longwave 305.3 (303.0; SHEBA about 280), net longwave −8.3 (−10.7; −40);
sensible +13.2 (+13.5; near 0), latent −1.7 (−1.8; about −5); net surface
flux 125.6 (129.0; about 40); the ocean's 4.6 (4.7).

Cost under the exclusive lock (`js/gpu/profile.module.js`, 128 steps after
16 from eight64_day0183 and eight128_day0183, alternated twice with
7c6a4a3 built by `git archive`): step median 21.51 and 21.55 against 20.44
and 20.43 ms at N=64 (+5.2 and +5.5 %), 95.35 and 95.44 against 90.31 and
90.35 ms at N=128 (+5.6 and +5.6 %); the physics and boundary-layer passes
3.90 against 2.80 ms and 16.74 against 11.69 ms. A day at N=128 from
eight128_day0183 takes 64 s of wall time (70 s with the 6 s setup). No run made
a NaN; nine64's day-94 clamp is the only one.

Parity (Oct 1), from both engines stepped phase by phase from one state
(nine64_day0091, the slab surface) and every GPU buffer compared across a
save and load. (1) Where the surface parcel tops out at the interface the
stratocumulus's descending parcel stops at, the coupling test compares
two computations of one height: equal on the CPU (coupled, the written
intent), an ulp apart in f32 on the GPU, which decoupled 13 of the 82
such columns after one step (44 of the 351 built for
`test/boundaryLayer.test.mjs`). After one step 13 regimes and 15 columns' lowest ten layers
(by up to 5.2 K) were apart, after four 58 columns (> 0.1 K) and 82 cells'
net surface flux (> 1 W/m², 17 of them poleward of 60°); with the
decision on the parcel top's level index, 0 after one step, 4 and 30 (3
poleward of 60°) after four. The Arctic test is unchanged by it: 9.191 →
8.248·10³ km³ over three days, 3.4 % a day, net surface flux 129.0 and
ocean flux 4.7 W/m², before and after; the southern pack from
eight64_day0183 gains 0.008·10³ km³/day (12.012 → 12.037; net −20.3,
ocean 22.4 W/m²) before and after. The CPU reference melts the Arctic
pack alike over the first day (9.191 → 8.896 against the GPU's 8.895;
net 120.2 and ocean 4.4 W/m² against 120.2 and 4.5): the 3.4 % a day is
the scheme's, not the GPU's. (2) A
reload lost the lowest wind speed (the first step's drag), the
evaporation (the ocean's next freshwater), the shallow cumulus (the next
radiation) and the ocean between its steps (Q and W rebuilt as h·T, the
flux into the ice zero until its next step); every snapshot now carries
them and the ocean's restart arrays (20.3 MB more on an N=64 day's file of
102.4 MB), and three one-day segments from eight64_day0183 end on day 186
byte for byte as one three-day segment (N=6 in
`test/asyncSpinup.test.mjs`); a split inside a day parted through the
rain and runoff the ocean takes as differences of running totals (carried
since, below). (3) Under the saturation adjustment's
cover, gray optics and maximum-random overlap, 1 − e^(−W/W_vis) is 0 in
f32 for trace cloud, so the GPU split the overlap blocks the CPU's expm1
joined (one column covered 0.797 against 0.560, every layer's cover
equal): 11 cells' absorbed sunlight apart by up to 146 W/m² after one
step, 1 cell (20 W/m², the deck's fog cover below) with the small-rate form.
Remaining: the mixed-layer deck's cover falls from 1 to 0.3 as its cloud
base leaves the surface in a fog whose buoyancy flux is negative
throughout, so the engines' rounding of saturation decides it; land
evaporated at the potential rate from the soil under any snow (below);
the dry adjustment's 10⁻⁶ tolerance and plume onsets part a few columns a
step; the ocean flux into the ice is f32 heat content less h·T_f, a
quantum of about 24 W/m² under a 600 m mixed layer (per-cell |Δ| 10 W/m²
over the southern pack, means within 2.2).

Review of the parity work (Oct 2), with both engines' own `step` from
the saved states (deck, rain, boundary-layer and cumulus fields loaded,
the ocean on at `everySteps` 8). Cells whose net surface flux parts by
more than 1 W/m² after 1, 4 and 16 steps: nine64_day0091 0, 32, 117;
eight64_day0183 0, 73, 141; nine64_day0365 0, 111, 216 (largest after
four 93.9, 148 and 224 W/m²); columns more than 0.1 K apart in the lowest
ten layers 0, 3, 28; 0, 9, 30; 0, 17, 66. Most of the cells after four
steps were cold land, 51 of eight64's 73 poleward of 60°: a trace of
snow set the land's wetness to 1 and drew a step's potential evaporation
from the bucket, 10 to 60 times the bucket's own (1.0·10⁻⁴ against
1.8·10⁻⁶ kg/m²/s at 45.7N 85.0E in March). The trace is a plume's rain
evaporating on its way down, 10⁻⁶²–10⁻³⁷ kg/m² on the CPU (244, 66 and 29
cells after the first step of the three states) and 10⁻³⁶–10⁻²³ on the
GPU (17 cells in eight64). Both engines now take the potential rate only
above `TRACE_SNOW`, 10⁻⁶ kg/m² (no GPU state holds land snow between 0 and
10⁻³ kg/m² at a day's end): 0, 24, 107; 0, 19, 68; 0, 26, 86 cells
(largest after four steps 22.4, 7.7 and 18.3 W/m²; poleward of 60° 2, 0
and 0), the columns 0, 3, 25; 0, 7, 21; 0, 7, 34. On the slab surface of
the paragraph above, after four steps: 82 cells (17 poleward) at
7c6a4a3, 36 (3) with the level decision, 28 (2, largest 5.6 W/m²) with
`TRACE_SNOW`. The cells left after four steps are thresholds that single
precision decides: the dry adjustment merging one layer more (40735,
the lowest layers' q 4.6·10⁻⁵ apart a step later), the cloud-top
threshold of 10⁻⁶ kg/kg crossed by a decaying deck (33635, depth 677
against 373 m), a shallow plume stopping on one engine (356). Under the
old cloud options (gray optics 55/130, maximum-random overlap, the
saturation adjustment, no falling ice) one step parts ASR by more than
1 W/m² in 11 cells (146 W/m²) at 7c6a4a3 and in 1 now (14086, 20 W/m²,
the deck's fog cover). Over iced cells one step closes the ice's surface
energy against dt (net + ocean flux) to 10⁻⁷ J/m² on the CPU and to
2·10⁻⁷ of the stored energy on the GPU, and the CPU's absorbed sunlight
equals (1 − α) times the light at the surface to 10⁻¹³ W/m². The CPU
melts the Arctic pack over the first day as above (9.191 → 8.896, 3.21 %
a day; net 120.2, ocean 4.4 W/m²). The Arctic test on the GPU gives the
table's numbers at 7c6a4a3 and with the parity fixes (9.191 → 8.248,
3.42 % a day; downwelling sunlight 208.7, albedo 0.387, absorbed 128.0,
net 129.0, ocean 4.7 W/m²), and with `TRACE_SNOW` 9.191 → 8.247 (0.315·10³
km³ a day, 3.42 %; 208.5, 0.387, 127.9, net 129.1, ocean 4.7); the
southern pack from eight64_day0183 12.012 → 12.037 (net −20.2, ocean
22.3 W/m²). Split runs through
`scripts/spinup.mjs` (N=64, bl34, `everySteps` 8): three one-day
segments from eight64_day0183 end day 186 byte for byte as one three-day
segment. A segment stopped after 128 steps and continued ends day 184
with every state field of the uninterrupted day (the day means alone
cover its last 128 steps) once the in-day snapshot carries the ocean's
running rain and runoff totals and what it has taken of them
(`rainTotal`, `runoffTotal`, `rainSeen`, `runoffSeen`, written back after
the segment's first frame); without them the core and ocean state parted
16 steps after the reload.

**Parity and the gas branch merged (Oct 2).** Branch integrate-c: parity
(300f8e1) on sweep2 (e695112), then gas-benchmark 5970f4f, both sides
kept in every hunk. `scripts/radiationBenchmark.mjs` prints 5970f4f's
output character for character; the gas parent's digest on grey ice is
a73d33f3 on both trees, the defaults' a64fbb13, the cloud parent's
97902663 under the gray gases with the earlier visible split and
Rayleigh bands. Tests: the deck's cloud field (`test/frameGpu.test.mjs`),
the six-step deck and the regime-gated deck (`test/gpuModel.test.mjs`)
pass under the spectral gases and moist defaults; the cloud effects
after 16 steps (`test/cloudEffect.test.mjs`) pass under the defaults,
compared over the columns whose dry adjustment merged the same layers
on both engines (0 of 362 part; under the gray optics and saturation
adjustment column 201 parts at step 14, the CPU leaving its layer 22
8.8·10⁻⁷ of θ under the 10⁻⁶ merge tolerance and the GPU merging it,
1.9 W/m²), as is the iced-cell run (0 of 362 part here; after the
parity merge alone column 356, its layer 24 left 9.4·10⁻⁷ under the
tolerance on the CPU);
the trace-cloud join compares the engines' change by the join,
1.3·10⁻⁵ K/day against 0.78 moved (the joined heating itself 2.3·10⁻⁴
apart against the 2.29·10⁻⁴ it was held to); the treeline warmth,
4.9·10⁻³ K apart, is bounded by the run's largest lowest-air difference,
2.3·10⁻² K (column 172, 10.3N, at step 46, snow-free on both engines).
Still failing: the cloudy columns' layer heating, 1.26·10⁻⁴ against
10⁻⁴ K/day at cell 196 layer 24 (the deck's layer, cover 1.6·10⁻⁶
apart), as on 5970f4f itself (1.28·10⁻⁴; 9.1·10⁻⁵ on 1f61d5c, 1.0·10⁻⁴
on 608c4e5, 1.3·10⁻⁴ from a069a68). The longwave carries 1.32·10⁻⁴ of
it (−36.48518 against −36.48601 W/m²); under gray longwave the worst is
1.3·10⁻⁴ at cell 340, with ozone 'idealized' (the test's state then
differs) 7.1·10⁻⁵; one f32 ulp more in every θ and q moves the CPU's own
heating there by 4.9·10⁻⁵, and the GPU's Exner function there is
5.4·10⁻⁸ below the CPU's (no threshold in it). By hand on the
midlatitude-summer profile (bl34, 76 lit N=4 columns, no aerosol) with
100 g/m² of liquid at 285 K and 20 g/m² of ice at 220 K: the visible
light's ozone loss on its way out (at μ 0.99 8.93 W/m² from the light
the cloud reflects, 0.14 from the surface's) and the near-infrared
Rayleigh's change of the reflected light (0.556 W/m²) are the engines'
to 3·10⁻¹⁶ of the beam on the CPU and 1.7·10⁻⁷ on the GPU, overcast and
under the uniform cover (column cover 1, and 0.749 with 2 g/m² of ice),
and on the CPU with the empirical deck blended in.
Closure: on eight64_day0183 and
nine64_day0091 the CPU closes the shortwave to 4.5·10⁻¹³, the layers'
shortwave to 4.5·10⁻¹³ and their longwave to 9.1·10⁻¹³ W/m², the GPU
to 1.8·10⁻⁴ and 2.3·10⁻⁴ W/m²; five64_day2190, seven64_day0365 and
m21a64_day0365 load and close alike. Parity from nine64_day0091 after
1, 4 and 16 steps: cells whose net surface flux parts by more than
1 W/m² 0, 25, 97 (parity branch 0, 24, 107), columns more than 0.1 K
apart in the lowest ten layers 1, 6, 36 (0, 3, 25; after one step
column 11539, 23.0N 160.9W, decoupled at 334 m on the CPU, coupled to
939 m on the GPU). Three one-day segments from eight64_day0183 end day
186 byte for byte as one three-day segment. Three GPU days (bl34,
`everySteps` 8), the last day's means, sweep2 in brackets: from
eight64_day0183 albedo 0.294 (0.298), ASR 240.6 (238.9), OLR 234.2
(233.2), atmosphere 84.6 (82.7), SWCRE −52.1 (−53.5), LWCRE 26.4 (26.6),
clear-sky reflectance 0.1405 (0.141), rain 1.68 (1.68), cover 0.54
(0.54), 60–90S 0.67 (0.68), sea surface 168.3 (168.6) and −44.5
(−44.3) W/m²; from nine64_day0091 0.305 (0.309), 236.8 (235.2), 235.0
(234.3), 87.3 (85.4), −54.6 (−56.0), 26.2 (26.4), 0.1443 (0.145), 2.21
(2.21), 0.55 (0.55), 0.65 (0.65), 147.8 (148.1) and −44.8 (−44.7). The
Arctic pack 9.191 → 8.272·10³ km³ over days 91–94, 3.3 % a day (3.3);
from eight64_day0183 1.672 → 1.677 and the southern 12.012 → 12.011.
The 0–2.2 hPa layer 259.9, 259.7, 259.2 K over days 184–186 and 260.6,
260.1, 259.5 over days 92–94 in the global mean; the 3.5 hPa layer
267.3 → 260.2 and 267.3 → 260.2, the tropical 85–101 hPa layers 208.8
and 211.4 (day 186), 206.6 and 209.3 (day 94). Cost under the exclusive
lock (128 steps after 16, alternated twice with e695112 built by `git
archive`): N=64 21.55 and 21.72 → 22.28 and 22.51 ms (+3.5 %), the
physics pass 3.89 → 4.67 ms; N=128 97.24 and 95.56 → 98.38 and 98.38 ms
(+2.0 %), the physics pass 16.74–17.22 → 19.74 ms. A day at N=128 from
eight128_day0183 takes 66 s (72 s with setup; sweep2 64 and 70).

**The surface layer by roughness (Oct 2).** `js/physics/exchange.module.js`
(CPU) and `js/gpu/exchange.gpu.js` (WGSL), parity and hand-value tests
in `test/exchange.test.mjs`.

Before: one constant per cell for momentum, heat and vapour alike,
`SEA_DRAG` 1.2·10⁻³ over sea and sea ice and `LAND_DRAG` 1.5·10⁻³ over
land and the ice sheets, at the lowest layer's midpoint (z₁ 16.6 m over
the ice sheets, 18–21 m elsewhere on the day-94 means). It entered the
explicit RK4 drag ρ C max(|v|, 3) |v|/m (`surface.module.js`, the core's
D_DRAG), the ocean's stress, the boundary layer's u* = √C max(|v|, 3) and
its surface buoyancy flux, the radiation column's sensible heat
ρ C U c_p (T_s − T₁) and evaporation, the land's aerodynamic conductance
and the FAO-56 reference evaporation. No stability dependence: the
boundary layer's `stability` changes only the K-profile's velocity scale
(Holtslag and Boville) in unstable columns, and the gustiness is the
3 m/s floor on the wind. The sensible heat used T_s − T₁, not the dry
static energy, about 0.2 K of difference at 20 m (IFS eq. 8.6 carries
g z₁/c_p). At 10 m the sea's constant is a neutral C_D10N of 1.36·10⁻³
(z₀ 1.9·10⁻⁴ m) at every wind. The day means of three N=64 days from
nine64_day0091 (bulk Richardson number of the surface layer, 10th / 50th /
90th percentile of the 8-step samples): forest (trees ≥ 0.5) wind 5.9 m/s,
stress 0.075 N/m², Ri −0.25 / −0.02 / 0.03; bare soil 6.8 m/s, 0.094,
−0.27 / −0.04 / 0.005; ice sheets 9.2 m/s, 0.174, −0.001 / 0.010 / 0.042;
the Arctic pack (≥ 70N, A ≥ 0.8) 5.3 m/s, 0.054, −0.007 / 0.004 / 0.035;
the open sea 20S–20N 5.2 m/s, 0.048, −0.14 / −0.03 / −0.002, 40–60S
10.1 m/s, 0.185, −0.05 / −0.007 / 0.006.

The scheme. Per cell, from the state at the start of the physics:

| surface | z₀m (m) | z₀h (m) | source |
|---|---|---|---|
| forest (the trees' share) | 2.0 | 2.0 | IFS Cy47r3 Part IV Table 8.3 (all tree types), calibrated so that the 10 m wind's error against SYNOP vanishes per type (Sandu et al. 2011) |
| grass (cover less trees) | 0.1 | 0.001 | same, short grass (the table's low vegetation spans 0.034 tundra to 0.5 irrigated crops and bogs) |
| bare soil | 0.013 | 1.3·10⁻⁴ | same, desert |
| snow on grass and bare soil, ice sheets | 1.3·10⁻³ | Andreas (1987) | same, ice caps and glaciers; snow covers the short tiles over min(1, S/30 kg/m²), the IFS's c_sn with D_cr 0.1 m at 300 kg/m³; the trees stand above it |
| sea ice of concentration A | max(10⁻³, 0.93·10⁻³ (1 − A) + 6.05·10⁻³ e^(−17 (A − 0.5)²)) | Andreas (1987) | IFS eq. 3.30 (Andreas et al. 2010, Bidlot et al. 2014) |
| open sea (1 − A) | α u*²/g + 0.11 ν/u*, α = 0.0017 U10N − 0.005 (U10N ≤ 19 m/s, α ≥ 0) | min(1.6·10⁻⁴, 5.8·10⁻⁵ Rr^−0.72) | COARE 3.5 (Edson et al. 2013; the scalar fit as in its coare35vn.m, where COARE 3.0's was min(1.15·10⁻⁴, 5.5·10⁻⁵ Rr^−0.6)), ν of the air at its temperature (COARE's fit), four fixed-point steps from the neutral u* |

Andreas (1987, Table 2 of Andreas 2002): ln(z₀h/z₀m) = b₀ + b₁ ln R* +
b₂ (ln R*)², R* = u* z₀m/ν with the tile's neutral u*, (1.25, 0, 0) for
R* ≤ 0.135, (0.149, −0.55, 0) to 2.5, (0.317, −0.565, −0.183) above (R*
capped at 1000). The tiles blend as the IFS aggregates its tiles'
roughness for its 10 m wind: the shares' neutral C_D and C_H at 10 m are
summed and the cell's z₀m and z₀h backed out (a forest–grass half and
half: z₀m 1.17 m, z₀h 1.03 m). Stability: Monin–Obukhov with the IFS's
surface-layer functions (eqs. 3.16–3.26: Dyer–Hicks integrated by
Paulson when unstable; Holtslag and De Bruin 1988 with a 1, b 2/3, c 5,
d 0.35 when stable; no cap on z/L), ζ = z₁/L from the bulk Richardson
number g z₁ (θv₁ − θv_s)/(θ̄v U²) by five steps (fixed point when
unstable, Newton in ln ζ when stable; worst error against 200 steps
over Ri −10…10 at z 16–40 m 4·10⁻⁴ for the tiles above and their blends
and 1.6·10⁻³ for any z₀m ≤ 2 m with z₀h/z₀m ≥ 10⁻³; below z 10 m with
z₀m 2 m five steps miss by up to 99 %, a height the model's z₁ of
16–21 m does not reach), U the
wind with the 3 m/s floor, θv_s with the skin's saturation humidity over
sea and sea ice and dry over land. z₁ = c_p θv₁ (Π_s − Π₁)/g. C_D and
C_H go to the momentum, the ocean's stress and u*; C_H to the sensible
heat and evaporation, the mixed-layer deck's surface fluxes and the
land's aerodynamic conductance. The reference evapotranspiration takes
FAO-56's own reference grass at z₁ (eq. 4: h 0.12 m, d 2/3 h, z₀m
0.123 h, z₀h 0.1 z₀m and FAO-56's κ 0.41; 2.45·10⁻³ at 20 m, 208/u₂ at
2 m).

Three consequences built with it, on both engines:
- The surface drag is the lower boundary of the boundary layer's
  implicit edge solve (`implicitDrag`, the default under roughness): a
  forest's C_D of 0.02–0.03 on the 40 m lowest layer gives the explicit
  RK4 drag λΔt up to 3–6 at N=64, past RK4's 2.8; at N=6 with Δt 900 s
  the explicit form went NaN within three steps. The solve keeps the
  stress it applies (`surfaceStress`), which the ocean takes (on the GPU
  `PH_STRESS`), and its kinetic energy goes to the lowest layer's
  dissipation heat; the RK4 drag and its heat are off.
- The sensible heat is ρ C_H U (c_p T_s − c_p T₁ − g z₁) (IFS eq. 8.6),
  and the boundary layer's surface buoyancy flux is the radiation's,
  g/θ (H/(ρ c_p Π) + 0.61 θ E/ρ) from the fluxes it applied
  (`surfaceBuoyancy`, on the GPU written to `PH_BUOY` by the physics
  kernel), moist over land too: the boundary layer diagnoses after the
  flux has heated the lowest layer, and with C_H U Δt/Δz of order one
  over forest a buoyancy from the heated layer read as stable.
- The land's season means read the lowest air before this step's
  surface and radiative heating.

Options (`SURFACE`): `exchange` 'roughness' (default) or 'fixed', the
constant coefficients bit for bit; an option set that names
`dragCoefficient` (surface or land) and no `exchange` is fixed, so the
sweep's sea-drag dial still runs as it did; `roughness` {forest, grass,
bare, snow}, `snowCover` 30 kg/m², `blendingHeight` 10 m, `charnock`
[0.0017, −0.005, 19], `smoothFlow` 0.11, `iterations` 5,
`referenceCrop` 0.12 m, `implicitDrag` true.

The sea's neutral 10 m coefficients against the 10 m neutral wind
(×10⁻³; the model's four steps and COARE 3.5's converged relation agree
to 0.4 %): 3 m/s C_D 0.92, C_H 1.10; 5: 0.92, 1.10; 7: 1.05, 1.10; 10:
1.32, 1.11; 15: 1.83, 1.17; 20: 2.35, 1.23; 25: 2.73, 1.28. Large and
Yeager (2004): 1.27, 1.06, 1.06, 1.18, 1.47, 1.81, 2.16. The constant
before: 1.36 at every wind.

C_D / C_H (×10⁻³) at z 20 m by the bulk Richardson number:

| surface | Ri −0.3 | −0.05 | 0 | 0.05 | 0.2 |
|---|---|---|---|---|---|
| forest | 60.8 / 88.3 | 36.4 / 41.4 | 27.8 / 27.8 | 16.0 / 16.0 | 2.29 / 2.17 |
| grass | 8.35 / 4.41 | 6.51 / 3.49 | 5.69 / 3.04 | 4.22 / 2.41 | 1.17 / 0.82 |
| bare soil | 4.22 / 2.61 | 3.40 / 2.11 | 2.97 / 1.83 | 2.12 / 1.38 | 0.65 / 0.48 |
| snow, ice sheet (u* 0.25) | 2.39 / 1.86 | 1.98 / 1.51 | 1.72 / 1.29 | 1.15 / 0.90 | 0.40 / 0.32 |
| sea ice A = 1 | 2.27 / 1.85 | 1.88 / 1.50 | 1.63 / 1.27 | 1.08 / 0.87 | 0.38 / 0.31 |
| sea ice A = 0.5 | 3.39 / 1.79 | 2.80 / 1.50 | 2.48 / 1.33 | 1.85 / 1.06 | 0.63 / 0.43 |

Runs, N=64 GPU, OCEAN `{"everySteps":8}`, from copies of the states,
before (HEAD 3db505a) → after, day means of 8-step samples over the
three days (wind at z₁; stress ρ C_D max(|v|, 3) |v|; H, LE upward
positive; the Arctic pack ≥ 70N with A ≥ 0.8):

| class, nine64 day 91 + 3 | share | wind m/s | stress N/m² | C_D / C_H ×10⁻³ | H W/m² | LE W/m² | skin °C |
|---|---|---|---|---|---|---|---|
| forest (trees ≥ 0.5) | 0.068 | 5.92 → 2.63 | 0.075 → 0.226 | 1.5 → 19.2 / 19.8 | 18.9 → 12.2 | 67.2 → 86.4 | 20.40 → 17.86 |
| grass (≥ 0.5) | 0.008 | 6.72 → 4.69 | 0.071 → 0.122 | 1.5 → 5.80 / 3.26 | 39.0 → 43.3 | 42.4 → 45.7 | 11.40 → 10.27 |
| bare (≥ 0.5) | 0.058 | 6.83 → 4.22 | 0.094 → 0.222 | 1.5 → 9.38 / 8.48 | 42.1 → 54.2 | 26.9 → 36.6 | 27.65 → 25.10 |
| mixed land | 0.126 | 6.46 → 3.30 | 0.086 → 0.227 | 1.5 → 13.8 / 13.9 | 30.7 → 31.0 | 58.1 → 78.2 | 24.12 → 21.28 |
| ice sheets | 0.030 | 9.17 → 9.60 | 0.174 → 0.215 | 1.5 → 1.59 / 1.09 | −19.9 → −19.1 | 2.8 → 2.3 | −36.23 → −36.23 |
| Arctic pack | 0.009 | 5.34 → 5.06 | 0.054 → 0.058 | 1.2 → 1.35 / 1.11 | −5.8 → −5.0 | 2.5 → 2.3 | −0.06 → −0.06 |
| Antarctic sea ice (A ≥ 0.8) | 0.015 | 9.93 → 9.80 | 0.202 → 0.271 | 1.2 → 1.58 / 1.16 | 3.1 → −0.2 | 5.5 → 5.4 | −19.72 → −19.63 |
| open sea 60–40S | 0.107 | 10.11 → 10.46 | 0.185 → 0.236 | 1.2 → 1.24 / 1.08 | 22.3 → 19.2 | 57.8 → 54.8 | 5.67 |
| open sea 20S–20N | 0.262 | 5.16 → 5.52 | 0.048 → 0.048 | 1.2 → 1.01 / 1.19 | 3.7 → 3.2 | 106.2 → 108.2 | 25.25 |
| open sea 40–60N | 0.052 | 5.90 → 6.41 | 0.066 → 0.061 | 1.2 → 0.80 / 0.85 | −14.0 → −10.5 | 15.6 → 14.4 | 10.04 |

| outcome | nine64 day 91 + 3 | eight64 day 183 + 3 | eight64 day 183 + 10 | nine64 day 365 + 3 | eight128 day 183 + 3 |
|---|---|---|---|---|---|
| forest wind, m/s; stress, N/m² | 5.92 → 2.63; 0.075 → 0.226 | 5.91 → 2.56; 0.080 → 0.255 | 6.99 → 2.68; 0.112 → 0.274 | 5.94 → 2.57; 0.086 → 0.291 | 6.03 → 2.50; 0.081 → 0.232 |
| forest LE, W/m²; skin, °C | 67.2 → 86.4; 20.40 → 17.86 | 46.3 → 62.9; 16.07 → 14.34 | 42.8 → 49.6; 14.59 → 12.81 | 59.6 → 78.1; 18.07 → 15.55 | 44.6 → 59.2; 15.29 → 13.57 |
| land skin, all land, last day °C | 17.53 → 14.85 | 15.89 → 13.92 | 13.70 → 11.50 | 7.74 → 6.12 | 15.59 → 13.61 |
| global surface, last day °C | 16.52 → 15.75 | 16.60 → 16.02 | 15.99 → 15.37 | 13.64 → 13.16 | 16.48 → 15.92 |
| global evaporation (samples), mm/d | 2.40 → 2.57 | 1.94 → 2.00 | 2.20 → 2.18 | 2.22 → 2.33 | 2.00 → 2.02 |
| global rain, last day, mm/d | 2.19 → 2.30 | 1.70 → 1.73 | 2.54 → 2.42 | 2.33 → 2.38 | 1.93 → 1.94 |
| global sensible heat, W/m² | 13.03 → 12.96 | 9.23 → 9.44 | 11.67 → 11.21 | 12.80 → 11.59 | 10.13 → 10.62 |
| sea stress magnitude, N/m² | 0.083 → 0.095 | 0.081 → 0.090 | 0.091 → 0.101 | 0.080 → 0.089 | 0.085 → 0.090 |
| equatorial Pacific τx 2S–2N 160E–100W, N/m² (Earth −0.04…−0.06) | −0.038 → −0.038 | −0.022 → −0.019 | −0.031 → −0.027 | −0.006 → −0.006 | −0.020 → −0.018 |
| trades τx 5–20N / 5–20S, N/m² | −0.006 / −0.036 → −0.007 / −0.035 | −0.009 / −0.053 → −0.008 / −0.053 | −0.018 / −0.056 → −0.016 / −0.056 | −0.035 / −0.015 → −0.032 / −0.015 | −0.017 / −0.058 → −0.016 / −0.055 |
| Southern Ocean τx 40–60S mean; peak band, N/m² (Earth about 0.17, SCOW) | 0.132 → 0.170; 0.170 → 0.215 at 47.5S | 0.111 → 0.137; 0.153 → 0.191 at 47.5S | 0.119 → 0.138; 0.136 at 52.5S → 0.154 at 47.5S | 0.093 → 0.114; 0.131 → 0.160 at 52.5S | 0.096 → 0.113; 0.169 → 0.213 at 57.5S |
| Arctic pack H, W/m² (SHEBA near 0 in June) | −5.8 → −5.0 | — | — | 2.0 → 0.3 | — |

The stresses on the sea are the ocean's (`OD_STRESS`, the stress it
receives under the ice's transmission), averaged over the samples.
Earth's Southern Ocean: the 40–60S decadal-mean zonal stress was about
0.17 N/m² in the 1990s, and the westerlies exceed 0.25 N/m² in July
(Risien and Chelton 2008, SCOW, citing Huang et al. 2006). The audit of
the end states (`scripts/verticalAudit.mjs`, each on its own code): day
94 global rain 2.36 → 2.49 and evaporation 2.41 → 2.50 mm/d (Earth
2.6–2.8), the Pacific ITCZ's rain 6.78 → 5.75 mm/d, SE Pacific low cloud
0.38 → 0.42, the zonal rain peak 9.12 → 9.81 mm/d, cloud effects −58.3 →
−59.7 and 20.1 → 20.1 W/m²; day 368 rain 2.46 → 2.63, evaporation 2.34 →
2.43, ITCZ 2.77 → 2.48, SE Pacific low cloud 0.24 → 0.30, zonal peak 5.51
→ 6.08 mm/d, shortwave effect −56.7 → −60.0; day 193 rain 2.48 → 2.44,
evaporation 2.44 → 2.37, ITCZ 6.07 → 9.93, SE Pacific low cloud 0.29 →
0.34, zonal peak 5.73 → 8.93 mm/d, cloud effects −58.8 → −59.4 and 19.1 →
19.9 W/m², the day's ASR − OLR −10.3 → −9.1 W/m²; N=128 day 186 rain
2.15 → 2.34, evaporation 2.21 → 2.25, ITCZ 3.60 → 3.75, SE Pacific low
cloud 0.21 → 0.19, zonal peak 4.97 → 5.67 mm/d, cloud effects −50.4 →
−50.8 and 16.8 → 16.7 W/m². The land's skin falls by 1.7–2.7 K and its
air much less: at day 193's end over the land off the ice sheets the
lowest air 17.50 → 17.15 °C and the skin 19.34 → 17.15 °C (the skin
was 1.8 K above the air and now sits on it), the global lowest air 15.77
→ 15.65 °C.

Step cost under the exclusive lock (each engine's own day-193 and day-186
state, 128 steps after 16, twice): N=64 27.31–27.33 → 27.31–27.61 ms a
step (7.0 → 7.0–7.1 s a model day), N=128 116.8–117.3 → 117.2–117.3 ms
(59.8–60.0 → 60.0–60.1 s a model day).

Parity (N=6, one step from random land, snow, trees and sea ice, 361
cells): C_D and C_H agree to 2–8·10⁻⁶ relative over open sea, forest,
other land, snow and the ice sheets and 6.5·10⁻⁵ over sea ice, the
reference to 8·10⁻⁶, latent heat to 0.04 W/m² and the net surface flux to
0.06 W/m²; eight steps over the continent on bl34: Ts 2.5·10⁻⁴ K, θ rms
4.4·10⁻⁷, column water rms 9.7·10⁻⁷. `exchange` 'fixed' reproduces the
digests unchanged. Tests that feed the boundary layer a hand-made
surface or compare deck regimes, rain sums and cloud effects that one
switching column moves run on 'fixed', the surface they were written for.

The review's checks (Oct 2). On the real state eight64 day 183, both
engines from one copy: the CPU's implicit drag changes each edge
column's momentum by Δt τ to 8·10⁻¹⁴ of the total |Δt τ| (worst edge
1.5·10⁻¹⁰ kg/m/s against a largest Δt τ of 2.2·10³), and the ocean receives the stored
stress unchanged (on the GPU the ocean's stress is the same PH_STRESS
times the ice's transmission). Engine parity (rms over the grid of the
lowest layer's θ and u and of Ts; the 'fixed' exchange's own in
brackets): after one step θ 1.2·10⁻² K (9.6·10⁻⁴) and u 3.5·10⁻³ m/s
(7.6·10⁻⁴), the largest differences in moist-boundary-layer columns of
the coupled regime switching between engines; after 4 steps θ 6.2·10⁻³
(3.8·10⁻³), u 3.3·10⁻³ (3.6·10⁻³), Ts 5.7·10⁻² K (5.7·10⁻²); after 16
θ 3.2·10⁻² (2.7·10⁻²), u 1.9·10⁻² (1.7·10⁻²), Ts 0.18 K (0.18); C_D and C_H after
one step 1.7·10⁻⁵ relative (largest 1.9·10⁻⁴), the stored stress
5.5·10⁻⁵ N/m², latent heat 2.3·10⁻³ W/m². 'fixed' on the CPU reproduces
3db505a bit for bit over nine N=64 steps from that state (every state
array, the land, the ocean and its stress), both by default and with
the sweep's `dragCoefficient` 1.3·10⁻³; on the GPU it differs from
3db505a in the last bits (357 θ values after one step, at most 3 f32
ulps), the shader compiler's rounding of the edited kernels. Workers
reproduce the serial CPU step bit for bit with land (N=8, six steps,
both exchanges). The drag and heat exchange over one step,
λ = ρ C U Δt/m₁ with m₁ the lowest layer's mass, after four steps
(N=64 Δt 337.5 s, N=128 Δt 168.75 s): forest λ_D mean 0.88 (largest
2.3; 4.7 on the first step, in the winds of the fixed drag), λ_H 0.95 (largest 2.2, above 1 in 47 % of the forest cells, above
2 in 0.6 %) at N=64 and λ_H 0.52 (largest 1.24) at N=128; the open sea
λ_H at most 0.32. The drag is implicit; the heat and vapour fluxes stay
explicit sources into the lowest layer. Over 96 N=64 steps the forest's
sensible heat has a mean step-to-step second difference of 0.7 W/m²
against 0.05 on 'fixed' (|H| 37–42 W/m²) and reverses in 0.4–0.5 % of
step pairs by more than 5 W/m² each way (0.0 % on 'fixed'), with no
growth; at N=6 and Δt 900 s, where λ_H reaches 13 over a continent of
forest, 192 steps stay finite with the lowest layer's θ reversing in
0.09 % of step pairs. The three-day runs from nine64 day 91 and eight64
day 183, re-run from copies (3db505a against cf1cdea), reproduce
the tables above to their last digit.

Subgrid orography (reported, not built). Climate models add two
stresses over mountains that the resolved terrain and the roughness
lengths above do not give. (1) Turbulent orographic form drag (Beljaars
et al. 2004, IFS eqs. 3.54–3.57): a drag ∂U/∂t = −C_tofd(z) |U| U spread
over the lowest kilometre or two, C_tofd = α β C_md C_corr 2.109
e^(−(z/1500)^1.5) a₂ z^(−1.2) with α 35, β 1, C_md 0.005, C_corr 0.6 and
a₂ from σ_flt, the standard deviation of the orography filtered to the
3–22 km band; for 10 m/s from the lowest level (20 m) up, σ_flt 25, 50,
100 and 200 m give 0.03, 0.13, 0.52 and 2.1 N/m² (an effective C_D at
20 m of 0.3, 1.2, 4.7 and 19·10⁻³), so over hills of σ_flt near 100 m it
doubles the grass's drag. σ_flt needs 1 km topography (the IFS's from
1 km data; GMTED2010 or SRTM30 would serve); the model's 0.25° file
cannot resolve the band. (2) Gravity-wave drag and low-level blocking
(Lott and Miller 1997, IFS eqs. 4.1–4.19), from the standard deviation
μ, anisotropy γ, orientation and mean slope σ of the subgrid orography
between the 0.25° data and the cell: the wave stress ρ U N (H_eff²/4)
(σ/μ) G (B cos²ψ + C sin²ψ), taken unblocked (H_eff = 2μ) with G 1 and
B π/4, at 10 m/s and N 0.01/s, is over N=64 land: μ < 50 m on 0.37 of the
land (0.003 N/m²), 50–100 m 0.21 (0.020), 100–200 m 0.19 (0.086),
200–400 m 0.16 (0.31), above 400 m 0.07 (1.37), 0.17 N/m² in the land
mean; at N=128 0.12 N/m² (less of the variance is subgrid). Blocking adds
C_d ρ (σ/2μ) U² over the blocked depth where N H/U > 1. These are the
order of the vegetation's turbulent stress itself over the mountains
(0.12–0.29 N/m² after); the μ, γ, θ and σ fields can be computed from the
0.25° file at the model's start.

What still misses: the gustiness is the 3 m/s floor, not a free
convection velocity (COARE's 1.2 w*, the IFS's w* with z_i 1000 m), so
calm stable nights keep a flux; the tiles share one skin, so a forest's
strong coupling and the bare soil's weak one average into one
temperature (the IFS solves each tile's skin); the grass is one IFS type
(0.1 m) where the table spans 0.034 (tundra) to 0.47 (tall grass); there
is no displacement height and no z/L cap; the snow's density is fixed;
the sea's Charnock coefficient takes the neutral u* and no waves; the
equatorial Pacific stress stays at −0.018…−0.038 against −0.04…−0.06 (C_D10N
at 5–7 m/s fell to 0.92–1.05 from 1.36 and the winds there rose by
0.3–0.4 m/s, so the stress stands; its deficit is the trades' own); and the stress over mountains lacks the orographic terms
above.

**The land and roughness branches merged (Oct 2).** Branch integrate-c:
clear-sky 8acba74, then roughness 3f57d8f, both sides kept in every
hunk; the GPU's physics parameters carry the ozone's year fraction in
slot 1 and the land's season and moisture weights and hold in 5–7.
Tests pass but the cloudy columns' layer heating (1.26·10⁻⁴ against
10⁻⁴ K/day). On exchange 'fixed' now: the frame's deck cloud (two
columns coupled on the GPU, decoupled on the CPU, 4.6·10⁻² kg/m²), the
boundary layer under a stratocumulus, the 48-step season means and trees
(1.2·10⁻⁴ in tree cover under roughness, 5.3·10⁻⁵ fixed); the engines'
flux check allows 10⁻⁵ of each cell's own flux (snow columns of 2·10⁴
W/m², 4.6·10⁻⁶ apart). Proofs: 'fixed' with `soilCarbon`, `grassland`
and `treeline` off reproduces b913e99 (`treeline` off) bit for bit over
nine CPU steps from eight64_day0183, ocean step included; the benchmark
prints b913e99's output character for character. On eight64_day0183 and
nine64_day0091 the CPU closes the shortwave to 4.5·10⁻¹³, the layers'
shortwave to 1.1·10⁻¹² and longwave to 5.7·10⁻¹³, the surface flux
against its terms to 9.1·10⁻¹³, the column heating applied against
atmosphere SW + LW + SH to 1.8·10⁻¹² and the vapour against LE to
9.1·10⁻¹³ W/m²; each edge column's momentum changes by Δt τ to
3.2·10⁻¹¹ of 1.3·10³ kg/m/s, and the ocean takes the stored stress on
every edge. The GPU closes to 1.8·10⁻⁴ (SW), 2.2·10⁻⁴ (LW) and 1.7·10⁻⁴
(surface) W/m², one physics kernel's column heating to 8.3·10⁻⁴ and its
vapour to 1.3·10⁻⁴ W/m², and the ocean's stress at its step is the
step before's stored stress times the ice's transmission to 1.2·10⁻⁷
N/m². Parity from nine64_day0091 after 1, 4 and 16 steps, rms of the
lowest layer's T and u and of Ts ('fixed' in brackets): 4.2·10⁻⁵,
4.9·10⁻⁵, 9.3·10⁻⁶ (8.6·10⁻⁴, 1.3·10⁻⁴, 8.7·10⁻⁶); 1.2·10⁻³, 5.5·10⁻⁴,
2.6·10⁻⁵ (5.7·10⁻⁴, 2.7·10⁻⁴, 2.4·10⁻⁵); 1.9·10⁻³, 1.0·10⁻³, 1.2·10⁻³
(1.8·10⁻³, 1.1·10⁻³, 1.6·10⁻³); C_D 1.5·10⁻⁵, 3.8·10⁻⁵, 1.0·10⁻⁴ rms
relative; regime flips 0, 3, 23 (1, 2, 24). Three one-day segments from
eight64_day0183 end day 186 byte for byte as one; a fresh ten-day start
with `LAND_JUMPS=5` split at day 5 ends day 10 byte for byte as unsplit.
five64_day2190, nine64_day0365 and ten64_day0183 load and step. Classes
on day 186: partly vegetated 0.200, dense vegetation 0.134, open sea
0.091 / 0.110 / 0.149, all land snow 0.637; day 94 partly vegetated
0.159.

Three GPU days, b913e99 in brackets: eight64 day 186 albedo 0.296
(0.294), ASR 239.8 (240.6), OLR 233.6 (234.2), SWCRE −53.9 (−52.1),
LWCRE 26.2 (26.4), clear-sky reflectance 0.1374 (0.1405), rain 1.76
(1.68), cover 0.54 (0.54), 60–90S 0.67 (0.67), Ts 16.31 (16.94), land
14.73 (16.87) °C; nine64 day 94 0.307 (0.305), 235.9 (236.8), 234.1
(235.0), −58.5 (−54.6), 26.1 (26.2), 0.1354 (0.1443), 2.39 (2.21), 0.56
(0.55), 0.64 (0.65), 15.98 (16.76), 15.63 (18.31); ten64 day 186 0.312
(0.311), 234.4 (234.7), 232.5 (232.7), −59.1 (−57.8), 26.1 (26.6),
0.1380 (0.1409), 2.54 (2.56), 0.56 (0.56), 0.61 (0.63), 14.61 (15.10),
11.72 (13.45). Day 94 by class, 8-step samples: forest wind 6.17 → 2.66
m/s, stress 0.079 → 0.232 N/m², LE 64.3 → 93.2 W/m², skin − air 3.10 →
0.27 K; grass 6.69 → 4.70, 5.76 → 3.94 K; bare 6.81 → 4.21, 4.50 → 1.55
K; global evaporation 2.54 → 2.77 mm/d. Southern Ocean τx 40–60S 0.130 →
0.167, 0.110 → 0.135, 0.131 → 0.163 N/m²; 2S–2N 160E–100W −0.037,
−0.023, −0.072 → −0.037, −0.023, −0.074. Arctic pack 9.191 → 8.281·10³
km³ over days 91–94 (8.272), the last day's surface SW at ≥ 70N 125.7
(119.7), net LW −7.6 (−6.6) W/m². Pacific ITCZ on ten64 day 186 + 1
CPU day (`scripts/tropicalHeating.mjs`): rain 6.48 (6.02) mm/d,
convective 0.24 (0.26), lowest-layer RH 0.81 (0.81; Jordan 0.88), T −1.3
(−1.4) K against Jordan at 1008 hPa, −2.8 (−2.8) at 848, +1.1 (+1.0) at
516; deep plume fired 0.160 (0.175), fired CAPE 170 (165) J/kg, the
firing columns' peak 973 (974) hPa. The land's reference coefficient
2.44·10⁻³ against 1.5·10⁻³ on bl34 (1.89·10⁻³ on five64's cam26): the
state's day-mean PET 3.31 → 4.00 (eight64), 4.14 → 4.78 (nine64), 2.73
→ 2.93 (five64) mm/d; land with P/PET < 0.2 by the estimated rain 0.138
→ 0.241, 0.218 → 0.290, 0.112 → 0.131, the mean moisture factor 0.433 →
0.369, 0.385 → 0.345, 0.616 → 0.596. Ten days from the atlas: day 10
albedo 0.347, ASR − OLR +7.2, rain 3.87 mm/d; the record 864000 s, its
land means season length 0.624 (land branch 0.546), warmth 7.24 (5.46)
K, rain 1.166 (1.219), demand 2.560 (1.777) mm/d, litter 0.215 (0.177),
decay 0.747 (0.572). Cost under the exclusive lock (alternated twice
with b913e99 built by `git archive`): N=64 22.47 and 22.51 → 22.53 and
22.60 ms, the physics pass 4.67 → 4.75; N=128 98.38 and 98.49 → 99.45
and 99.53 ms, 19.76 → 20.37; a day at N=128 from eight128_day0183 73 s
with setup on both.

Review of the two merges (Oct 2). Every conflict rebuilt with `git
merge-tree`: each line either side added since its merge base is in the
merged tree or in the hunk that joins both, and no line either side
removed is back. The suite passes in concurrent per-file runners but
the cloudy columns' layer heating (1.2592·10⁻⁴ K/day, alone as in the
suite). On the CPU, 'fixed' with `soilCarbon`, `grassland` and
`treeMoisture` off gives b913e99's state and land digests after each of
nine steps from eight64_day0183 with the treeline on as well as off. The
GPU under the same options is not bit for bit: after one step the
sensible heat differs in 20051 of 40962 cells by at most 6.1·10⁻⁵ W/m²
and θ in 475 of 1392708 values by at most 2.4·10⁻⁴ K, and each merge
alone does as much (the land merge: 22691 cells, 6.1·10⁻⁵ W/m²). Parity
from eight64_day0183 after 1, 4 and 16 steps, rms of the lowest layer's
T and u and of Ts (b913e99 in brackets): 4.0·10⁻⁵, 4.9·10⁻⁵, 1.5·10⁻⁵
(3.4·10⁻⁴, 5.4·10⁻⁵, 1.5·10⁻⁵); 1.8·10⁻⁴, 8.8·10⁻⁴, 7.4·10⁻⁴ (1.8·10⁻⁴,
8.9·10⁻⁴, 6.9·10⁻⁴); 1.7·10⁻³, 2.7·10⁻³, 1.3·10⁻³ (2.1·10⁻³, 1.7·10⁻³,
2.2·10⁻³); regime flips 0, 2, 14 (0, 2, 12). One step from
nine64_day0091: each edge column's momentum changes by Δt τ to
3.4·10⁻¹¹ kg/m/s on the CPU; the GPU's stored stress is the CPU's to
4.8·10⁻⁶ rms and 1.5·10⁻⁴ N/m² at most; the GPU's net surface flux is
its shortwave, longwave, sensible and latent terms to 4.5·10⁻⁴ W/m²
(forest) and the engines' net surface fluxes agree to 0.16 W/m² there;
the ocean's stress at the eighth step equals the stored stress of the
seventh on all 79133 open edges. Splits from eight64_day0183 to day 185:
at the day with `LAND_JUMPS=184` byte for byte; inside day 185 (a
step-48 snapshot) with and without that jump, and from a fresh start to
day 2 with `LAND_JUMPS=1` split inside day 2, every field but the
saved last-interval means (the two rains, ASR, OLR, albedo and the two
cloud effects, which cover the segment's part of the day) is byte for
byte; b913e99 splits inside a day with the same seven fields apart. The
three days from eight64_day0183 repeat the day lines above to every
printed digit.

**The apparent heat source, the replicate spread and the longwave
overlap (Oct 2).** The tropical boxes are now measured by the box-mean
apparent heat source of the physics Q1, the same without the shortwave
and longwave Q1R, and Q2 = −(L/c_p) dq_t/dt of the physics, each in K/day
with its peak over 50 hPa bins of the mass-weighted mean and its centroid
Σ p Q dp / Σ Q dp over its positive part (`heatingProfile` in
`js/audit.module.js`); `scripts/tropicalHeating.mjs` prints them with the
column integrals of Q1R and Q2 by process, the total water's budget by
process and layer, the stratiform share of the rain (melted falling ice
and conversion above 700 hPa) and an Amazon land box (10S–2N 70–50W) with
the local solar hour of its convective rain, and `scripts/verticalAudit.mjs`
gives the ITCZ and warm-pool Q1R peak and centroid over its window in
place of the firing columns' convective heating peak. That row was the
layer maximum of the convective trace in K/day without mass weights: on
the diagnosis state (tdb64 rebuilt from b2310ca by `git archive`, its
three-day log identical line for line) it reads 974 hPa, 13.67 K/day, of
which the shallow plume is 12.45; on the gray-gas state 974 hPa again
(6.84). On the same diagnosis state the ITCZ Q1 peaks at 600–650 hPa
(1.81 K/day), Q1R at 704 hPa (2.64, bin 700–750) with 1.29 K/day at 439
hPa and its centroid at 677 hPa; the warm pool's Q1R peaks at 607 hPa
(3.35) with its centroid at 628; the gray-gas state's ITCZ Q1R at 439 hPa
(2.23), centroid 606; convective share 0.13, firing 0.099, mean dilute
and undilute CAPE 43 and 133 J/kg; the replay matches the model in every
column-step, the heat closes to 1.1·10⁻¹³ K and q_t to 4.3·10⁻¹⁹ kg/kg a
step. Over the same eight steps both scripts give the ITCZ 850–900 hPa
(2.14 K/day), centroid 724 hPa, and the warm pool 750–800 hPa (3.34),
centroid 654. Replicates of the three-day N=64 run from eight64_day0183
at this tree's physics (b913e99), the sea drag ×(1 ± 10⁻⁴), budgets over
day 186 → 187, base and range of the three: ITCZ convective share 0.1327
(0.0051), firing 0.0979 (0.0022), Q1R peak 704 hPa, bin 700–750 at 2.643
K/day (0.024), centroid 680.4 hPa (2.3); T − Jordan at 848 / 704 / 516 /
439 hPa −2.40 / −0.09 / +1.71 / +1.71 K (0.002 / 0.009 / 0.005 / 0.005),
RH 0.937 / 0.655 / 0.667 / 0.634 (≤ 0.0006); large-scale rain converted
below 700 hPa 2.248 mm/d (0.030); warm pool share 0.1689 (0.0024), firing
0.2161 (0.0022), centroid 627.6 hPa (0.8); day-186 global rain 1.681 mm/d
(0.0024), SWCRE −52.11 (0.030), LWCRE 26.37 (0.008), ASR − OLR 6.35
(0.024) W/m². The longwave's exponential-random overlap
(`longwaveOverlap`, above, and `physics.gpu.js`): the in-model fluxes
repeat `scripts/longwaveOverlap.mjs` to 8.3·10⁻¹⁶ (OLR) and 6.3·10⁻¹⁶
(surface downward) relative in every column; on the rebuilt diagnosis
state after one CPU step OLR +2.20 W/m² above random, LWCRE 26.22 →
24.02, surface downward −2.82, ITCZ +6.29, warm pool +2.71; on the day-186
state of the base +2.20, 25.77 → 23.57. Three days from eight64_day0183,
against the base: LWCRE 26.37 → 24.07, OLR 234.20 → 236.43, ASR − OLR
6.35 → 4.51 W/m², rain 1.680 mm/d; ITCZ convective share 0.1327 → 0.1210,
firing 0.0979 → 0.0872, Q1R centroid 680.4 → 684.6 hPa, the 300–500 hPa
longwave heating −1.838 → −1.822 K/day (warm pool −1.830 → −1.820, share
0.1689 → 0.1526). Cost under the exclusive lock (128 steps after 16 from
nine64/nine128_day0183, alternated twice with `longwaveOverlap: 'random'`,
which is bit-identical to b913e99 over 32 GPU steps): N=64 22.60 and 22.51
→ 25.90 and 26.34 ms (+16 %), the physics pass 4.68 → 8.00 ms; N=128
99.39 and 99.13 → 111.24 ms (+12 %; the second 125.03, minimum 109.85),
the physics pass 19.93 → 32.07 ms. Layers clear or overcast on both sides
taken as one region, three loops per layer, measured 28.60 and 28.92 ms
and 120.52 and 120.68 ms and is not kept. Tests: the layers' longwave
closes on the surface emission less the back radiation and the OLR to
3·10⁻¹³ W/m² on the CPU and 1.8·10⁻⁴ on the GPU; the cloudy columns'
layer heating parts by 1.26·10⁻⁴ K/day at cell 196 layer 24 as at
b913e99; over the treeline test's 48 GPU steps the lowest air parts by
7.7·10⁻² K (2.3·10⁻² under 'random') and the tree cover by 1.2·10⁻⁴
(2.5·10⁻⁵), above its 10⁻⁴.

**Independent checks of the Q1 metric and the overlap (Oct 2).** The
longwave's chain against the expectation over every sub-column pattern
(each an overcast or clear column, weighted by the chain's pair
transitions), at N=6 on 335 ice-free columns with two adjacent partly
cloudy layers near 500 hPa, 210 of them under a deck in layer 24, half
with resolved water in the deck's layer and the one above it: OLR and back
radiation agree to 2.1·10⁻¹⁵ relative and the layer longwave to 10⁻¹³ W/m²
on the CPU, to 5·10⁻⁶ relative and 1.1·10⁻³ W/m² on the GPU (1.05·10⁻³
under 'random'), for z₀ the default, 2000 m, 10⁻⁶ m (random) and 10¹⁵ m
(maximum); at z₀ 10⁻⁶ m the CPU's OLR and back radiation are bit-identical
to 'random' in 244 of 362 columns and within 4.1·10⁻¹⁶ in the rest. One
step from eight64_day0183 on each engine: the CPU's layers close on the
surface and top fluxes to 5.4·10⁻¹³ W/m², the GPU's to 1.4·10⁻³
(2.3·10⁻⁴ under 'random'); the engines' OLR differ by at most 2.1·10⁻³
W/m² and the layer longwave by 8.2·10⁻³ of 121 W/m² under either overlap.
On the default run's day-186 state (three days from eight64_day0183)
99.5 % of the columns have two adjacent layers of cover strictly between 0
and 1, and every 32-cell group has one, so no collapse to one region can
take back much of the cost; one reused vec2 array for the two regions'
fluxes made the physics pass slower (43.5 ms at N=128 against 32.0). The
three-day runs repeat ceb0 (under 'random') and ceb1 line for line, and
tropicalHeating.mjs repeats their budgets; a separate replay of the ITCZ
from phase snapshots gives the same Q1 (600–650 hPa, 1.761 K/day), Q1R
(700–750 hPa, 2.643, centroid 680.4 hPa) and Q2 (700–750 hPa, 2.491) on
the base, with Q1R's column 144.10 W/m² = L P 138.72 + sensible 3.26 +
the physics' condensate gain 0.87 + dissipation 0.91 + 0.34, and Q2's
column L(P − E) to the printed digits.

**The merged tree (237263d) proved (Oct 2).** Of 59 test files three
failed. `bulkSensible` (`js/audit.module.js`) still took the fixed
`SEA_DRAG`/`LAND_DRAG` and c_p (T_s − T₁): tropicalHeating.mjs's global
check on the audit test's N=12 state read 62.4969 against the model's
92.8688 W/m², the difference landing in the lowest layer between the
sensible and the shortwave terms; it now takes the exchange's C_H and,
under roughness, the dry static energy with the lowest layer's height read
before the physics phase (`lowestHeight`): 92.8688 against 92.8688, and on
the base's day-186 state 15.7736 against 15.7736. The bl34 checkpoint test
(`test/levels.test.mjs`) went NaN on day 2 at N=6: at cell 51 (42N 116E,
skin 237–247 K) C_H max(U, 3) Δt/Δz of the lowest layer reached 9.1 and
then 28 as C_H swung between 2.4·10⁻⁴ and 0.12 with the surface layer's
stability, the lowest air fell to 157.7 K and both engines were
non-finite by step 17–19; the same run is NaN on day 5 at the land
parent and day 6 at the overlap parent (winds 198 m/s), so the test now
runs its spin-ups on the fixed exchange. Three GPU–CPU cases of
`test/gpuModel.test.mjs` parted at single columns: at cell 325 (rain
accumulation, step 23) the GPU stepped from the CPU's exact state gave OLR
241.54 against 257.66 W/m², having found cloud water in layers 21 and
23–26 that the CPU had not (the cover floor 0.01, and the boundary layer's
variance cover 0.33 in layer 21 under 9.5·10⁻⁴ kg/kg of cumulus water); the
cases now leave out the columns whose OLR or absorbed sunlight part by more
than 1 W/m² at a step (every other column within 0.11 W/m²) and assert
their share. The sunlit cloudy columns' heating limit is 1.5·10⁻⁴ K/day
(1.26·10⁻⁴ at cell 196 layer 24 on both parents). The treeline test passed
on this tree (tree cover 1.6·10⁻⁵).

**The base on this tree and its replicate spread (Oct 2).** Three days at
N=64 from eight64_day0183 (`cvb0`) and two replicates with the sea's
Charnock coefficients ×(1 ± 10⁻⁴) (under the surface layer by roughness the
sea has no single drag coefficient), budgets over day 186 → 187 by
tropicalHeating.mjs, base (spread of the three): Pacific ITCZ rain 5.12
(0.11) mm/d, convective share 0.1141 (0.0017), firing 0.0868 (0.0003),
Q1R peak bin 700–750 hPa at 2.776 K/day (0.073), centroid 675.0 hPa
(0.64), large-scale rain converted below 700 hPa 2.388 (0.054) mm/d, T −
Jordan at 848 / 704 / 516 / 439 hPa −2.32 / +0.19 / +1.91 / +1.92 K
(0.003 / 0.007 / 0.011 / 0.008), RH 0.926 / 0.637 / 0.647 / 0.614
(≤ 0.002), RH at 946 / 963 hPa and the lowest layer 0.863 / 0.840 / 0.779
(≤ 0.0004), the shallow plume's export below 950 hPa 1.858 (0.005) mm/d,
mean dilute and undilute CAPE of the cloudy plumes 40.4 (0.13) and 132.7
(0.24) J/kg, undilute plumes stopping at 700–800 hPa 0.582 (0.002), fired
tops above 300 hPa 0.525 (0.020), wettest cell 157 (3) mm/d; warm pool
rain 9.79 (0.08), share 0.145 (0.004), firing 0.206 (0.003), centroid
620.2 (1.0); N Pacific trades rain 1.616 (0.021), share 0.202 (0.007),
firing 0.176 (0.003); Amazon rain 0.279 (0.005), all convective, peaking
at 14 LT; the replayed day's global rain 2.406 (0.0003) mm/d. GPU day-186
means: ASR − OLR 4.36 / 4.28 / 4.26, SWCRE −53.46 / −53.55 / −53.56,
LWCRE 23.87 / 23.88 / 23.87 W/m².

**The deep plume's closure of Bechtold et al. (2014) (Oct 2).**
`capeClosure` 'bechtold' (the default; 'threshold' keeps `plumeCape` 120
J/kg over `plumeRelaxation` 1 h bit for bit on the CPU) on both engines:
M_b = max(0, PCAPE − PCAPE_bl)/(τ F_P), PCAPE = Σ (T_v,u − T_v)/T_v Δp over
the layers the CAPE counts, F_P its change per unit base flux from the
scheme's own tendencies, τ = α_x H/w̄ within 720–10800 s, α_x = 1 + 1.66
dx/125 km, PCAPE_bl = τ_bl/T* Σ dT_v/dt|nc Δp below the plume's base (at
most the lowest 12 layers), T* = 1 K, τ_bl = z_base/max(ū_bl, 2 m/s) over
sea and sea ice and H/w̄ over land; IFS Cy43r1 eqs 6.22–6.29, all IFS
choices. dT_v/dt|nc is each subcloud layer's change of T_v since the end
of the previous adjustment (`subcloudVirtual`, PH SUBTV), saved in
spin-up states. On Jordan's column: PCAPE 199.05 Pa, H 11383 m, w̄ 6.301
m/s, τ 6631 s at N=32 and 4214 s at N=64, the flux scaling as 1/α_x; a sea
column under +2 K/d below cloud base with ū_bl 5 m/s: PCAPE_bl 23.9 Pa;
a land column under +10 K/d: 1583 Pa over τ_bl 1806 s, no deep flux.
Three days from eight64_day0183 against the base: ITCZ convective share
0.114 → 0.618, firing 0.087 → 0.375, large-scale rain below 700 hPa 2.39 →
0.78 mm/d, T − Jordan at 516 / 439 hPa +1.91 / +1.92 → +1.01 / +1.10 K,
wettest ITCZ cell 157 → 106 mm/d, median τ 85 min, Q1R centroid 675 → 703
hPa, fired tops above 300 hPa 0.53 → 0.15; warm pool share 0.145 → 0.446;
trades convective rain 0.33 → 1.34 mm/d, their Q1R peak bin 1000–1050 →
850–900 hPa; Amazon convective rain none at 11–16 LT, starting at 17–18 LT
and peaking at 0 LT; the replayed day's global rain 2.41 → 2.16 mm/d; GPU
day 186 ASR − OLR 4.36 → 9.93, SWCRE −53.46 → −48.46, LWCRE 23.87 → 24.33
W/m². From ten64_day0183: ITCZ share 0.588, warm pool 0.918, zonal-mean
rain peak 7.62 mm/d at 8–10N. Acceptance met: ITCZ share, large-scale rain
below 700 hPa, the 439–516 hPa bias, the wettest cell, the trades' Q1R
peak, τ, the ten64 shares and peak; not met: warm-pool share, ITCZ firing,
global rain, trades' convective rain; the Amazon peak later but past
15–18 LT.

**The deep plume's source: the lowest 50 hPa with the IFS surface-flux
excess (Oct 2).** `plumeSourceDepth` 'surface50' (the default;
'boundaryLayer' the previous source, bit for bit on the CPU): the deep
plume leaves the layers whose midpoints lie within 50 hPa of the surface
with their mean s_l and q_t plus ΔT = min(3 K, 1.5 J_s/(ρ c_p w*)) and
Δq = min(2 g/kg, 1.5 J_q/(ρ L w*)), w* = max((B₀ h)^⅓, u*), from the cell's
surface sensible and latent fluxes (Cy43r1 §6.5, eqs 6.19–6.21; IFS
coefficients). On Jordan with 10 and 130 W/m² at w* 0.585 m/s: ΔT 0.0217 K,
Δq 0.113 g/kg, CAPE 301.5 J/kg plain and 321.2 with the excess (273.0 from
the boundary layer). On element 2's day-186 state after one step the ITCZ
source's h/cp rises 0.73 K (0.50 the cut, 0.23 the excess), CAPE on the 257
columns with CAPE from both sources 46.8 → 48.0 → 53.0 J/kg and the
columns with CAPE 259 → 397 of 687 (warm pool 52.6 → 60.2 J/kg, 506 → 699
of 851). Three days against element 2: mean candidate dilute CAPE 32.1 →
29.3 J/kg (the new weak candidates and the stronger convection's
consumption), undilute plumes stopping at 700–800 hPa 0.407 → 0.386, fired
tops above 300 hPa 0.148 → 0.133, trades deep firing 0.518 → 0.255 and
convective rain 1.34 → 1.51 mm/d, ITCZ share 0.618 → 0.685, firing 0.375 →
0.332, large-scale rain below 700 hPa 0.78 → 0.62 mm/d, T − Jordan at 516 /
439 hPa +0.93 / +1.01 K, warm-pool share 0.446 → 0.523; Amazon peak 0 →
18 LT; global rain 2.16 → 2.11 mm/d; GPU day 186 ASR − OLR 10.00, SWCRE
−47.96, LWCRE 23.96 W/m². From ten64_day0183: ITCZ share 0.658, warm pool
0.976, zonal peak 7.40 mm/d at 8–10N. Acceptance met: the 700–800 hPa
stops, the trades' deep firing, the Amazon peak not earlier; not met: the
candidate CAPE (+20 %), the tops above 300 hPa (+0.05), the trades'
convective rain.

**The shallow cumulus base flux at Grant's 0.03 (Oct 2).**
`cumulusClosure` 0.06 → 0.03: M_b = ρ_LCL c w exp(−CIN/w²) with c the
coefficient Grant (2001, QJRMS 127, 407–421) fitted to LES (M = 0.03 w*)
and the inhibition factor of Bretherton, McCaa and Grenier (2004) kept, a
combination of the two published forms; 0.06 restores the previous
closure bit for bit on the CPU. The trade-wind column lifts 0.02087 against
0.04174 kg/m²/s. Three days against element 3: the ITCZ shallow plume's
export below 950 hPa 1.892 → 1.440 mm/d, RH at 946 / 963 hPa 0.896 / 0.866
→ 0.898 / 0.882, the lowest layer's 0.800 → 0.815, convective share 0.685
→ 0.823, firing 0.332 → 0.476, large-scale rain below 700 hPa 0.62 → 0.29
mm/d, Q1R centroid 701 → 715 hPa, wettest cell 105 → 46 mm/d; warm-pool
share 0.523 → 0.666; N Pacific trades low-cloud cover (the radiation's
lowCover over 8 CPU steps) 0.164 → 0.128; SE Pacific low cloud
(verticalAudit.mjs) 0.449 → 0.394, radiative 0.360 → 0.318, rain 0.28 →
0.26 mm/d; global rain 2.11 → 2.05 mm/d; GPU day 186 ASR − OLR 11.80,
SWCRE −45.85, LWCRE 23.79 W/m². From ten64_day0183: export 2.31 → 1.50
mm/d, RH at 946 / 963 hPa 0.904 / 0.876 → 0.918 / 0.901. Acceptance met:
the lowest layer's RH; not met: the export (≤ 1.1), the 946–963 hPa RH
(+0.025), the trades' low cover (≥ 0.18) and the SE Pacific guard (no more
than 0.03 below element 3).

**Elements 2–4 together (Oct 2).** Cost under the exclusive lock, 128
steps after 16 from nine64 and nine128_day0183, alternated twice with the
previous convection (`capeClosure` 'threshold', `plumeSourceDepth`
'boundaryLayer', `cumulusClosure` 0.06): N=64 27.05 and 27.12 → 27.54 and
27.55 ms (+1.7 %), the adjust pass 4.90 → 5.31 ms; N=128 113.66 and 112.76
→ 114.17 and 113.85 ms (+0.7 %), the adjust pass 17.23 → 18.38 ms. Over
the Amazon the boundary-layer part is negative at night: the fired
columns' mean PCAPE_bl is −1823 Pa against a PCAPE of 1.3 Pa (element 3;
−1884 and 2.8 with element 4), τ_bl = H/w̄ ≈ 30 min turning a cooling
subcloud layer into a flux far beyond what the PCAPE asks; 0.42 of the box's
convective rain falls at 0–6 LT (0.41 with element 4) and none at 11–16 LT. Global rain of the
replayed day 2.41 (base) → 2.05 mm/d, of which convective 0.27 → 0.79.

**PCAPE_bl at least 0 (review, Oct 2).** `pcapeBoundary` 'positive' (the
default; 'signed' is elements 2–4 bit for bit on the CPU) takes
max(0, Σ dT_v/dt|nc Δp) in PCAPE_bl. Bechtold et al. (2014, ECMWF Tech.
Memo. 705, §2b) define PCAPE_bl as the boundary-layer production of PCAPE
that shallow convection takes up, and report that the closure barely
changes the convection at night; the IFS-derived scheme of WRF
(`module_cu_ntiedtke.F`, its non-equilibrium branch) sets
`zcape2 = max(0, zcape2)`. Signed, a cooling subcloud layer adds to the
PCAPE: on Jordan's land column under −10 K/d PCAPE_bl is −1124 Pa against
a PCAPE of 220 Pa and the base flux 0.0963 against 0.0243 kg/m²/s with no
tendency (bounded: 0 Pa and 0.0243). Over one CPU day from element 4's
day-186 state the deep plume fired on 235 106 column-steps with no PCAPE
(none bounded); 0.76 of the fired land column-steps were held at the
boundary-loss or Courant limit (0.02 bounded) and 0.97 of the tropical
land base flux was beyond what the PCAPE alone asks; over sea 0.29 of the
fired column-steps had PCAPE_bl < 0, their flux a median 1.5 times the
PCAPE's. Three days from eight64_day0183 against element 4: Amazon mean
fired PCAPE_bl −1884 → +1.6 Pa, dilute CAPE 6.6 → 25.7 J/kg, firing 0.39 →
0.21, rain 0.32 → 0.47 mm/d, fired tops above 300 hPa 0.01 → 0.16, its
convective rain still none at 11–15 LT with its maximum at 22 LT; ITCZ
convective share 0.823 → 0.826, firing 0.476 → 0.496, Q1R centroid 714.7
→ 718.4 hPa; warm-pool share 0.666 → 0.678; N Pacific trades deep firing
0.704 → 0.757, convective rain 1.60 → 1.59 mm/d; replayed day's global
rain 2.049 → 2.041 mm/d, wettest cell 258 → 242 mm/d; GPU day 186 ASR −
OLR 11.8 → 11.9, SWCRE −45.8 → −44.9, LWCRE 23.8 → 23.1 W/m².

**The source excess over the IFS's own w\* (review, Oct 2).** Element 3
took eq. 6.19's coefficients with the shallow closure's w* = max((B₀
h)^⅓, u*); the IFS forms w* at the lowest model level (Cy43r1 eq. 6.20):
w* = 1.2 (u*³ + 1.5 g z κ/T (J_s/(ρ c_p) + 0.61 T J_q/(ρ L)))^⅓ with u*
0.1 m/s, and its code (WRF's IFS-derived `module_cu_ntiedtke.F`) gives the
parcel an excess only under an upward buoyancy flux, each part at least 0.
`excessVelocity` 'surfaceLayer' (the default; 'convective' is the previous
excess bit for bit on the CPU) does so on both engines. On Jordan's column
with 10 and 130 W/m² and the lowest layer at 21.0 m: w* 0.238 against
0.585 m/s, ΔT 0.0532 against 0.0217 K, Δq 0.278 against 0.113 g/kg, CAPE
344.4 against 321.2 J/kg (301.5 plain). Three days from eight64_day0183
against the bounded PCAPE_bl above: ITCZ mean candidate dilute CAPE 39.4 →
40.4 J/kg, undilute 306 → 342 J/kg, undilute plumes stopping at 700–800
hPa 0.261 → 0.229, fired tops above 300 hPa 0.072 → 0.070, convective
share 0.826 → 0.855, firing 0.496 → 0.504, large-scale rain below 700 hPa
0.291 → 0.228 mm/d, Q1R centroid 718.4 → 717.5 hPa, wettest cell 41 → 32
mm/d; warm-pool share 0.678 → 0.717; N Pacific trades deep firing 0.757 →
0.759, convective rain 1.59 → 1.62 mm/d; Amazon convective rain none at
11–15 LT, maximum 22 LT; replayed day's global rain 2.041 → 2.036 mm/d,
wettest cell 242 → 221 mm/d; GPU day 186 ASR − OLR 11.9 → 12.4, SWCRE
−44.9 → −44.4, LWCRE 23.1 → 23.0 W/m². From ten64_day0183 against element
4: ITCZ share 0.750 → 0.786, warm pool 0.991 → 0.988, zonal-mean rain peak
7.15 mm/d at 8–10N, Amazon mean fired PCAPE_bl −1556 → +23 Pa and dilute
CAPE 44 → 81 J/kg.

**The mountains' drag (Oct 2).** `js/physics/orography.module.js`
(CPU) and `js/gpu/orography.gpu.js` (WGSL) with the fields of
`subgridOrography` in `js/geography.module.js`; the surface layer's
gustiness and land humidity in `exchange.module.js` and
`exchange.gpu.js`; tests in `test/orography.test.mjs` and
`test/exchange.test.mjs`.

The diagnosis, before (3f57d8f), ten N=64 GPU days with 8-step samples
(u the zonal-mean eastward wind, the lowest layer for the surface and ln p
interpolation for the levels, where a level lies below a column's lowest
midpoint the lowest wind; SLP as the frames reduce it; stress the implicit
drag's on the land cells):

| nine64 day 274 + 10 (21–31 Dec) | 30–35N | 45–50N | 50–55N | 55–60N | 60–65N | 70–75N |
|---|---|---|---|---|---|---|
| u lowest layer, all / land, m/s | 1.0 / 1.2 | 4.0 / 3.0 | 3.5 / 2.0 | 2.4 / 1.7 | 1.8 / 1.6 | 0.3 / 1.3 |
| u 850 / 500 / 200 hPa, m/s | 4.1 / 14.2 / 19.7 | 8.9 / 16.9 / 28.6 | 8.1 / 17.4 / 33.0 | 6.1 / 15.6 / 32.3 | 4.5 / 11.9 / 24.7 | 0.9 / 2.6 / 7.2 |
| SLP zonal mean, hPa | 1021.0 | 1012.2 | 1007.5 | 1003.9 | 1001.0 | 996.9 |
| SLP stationary wave rms; wave 1 / 2, hPa | — | 6.2; 1.7 / 6.3 | 7.1; 4.6 / 7.2 | 7.5; 6.2 / 7.2 | 7.4; 6.4 / 6.4 | — |
| land stress magnitude; eastward, N/m² | 0.178; 0.071 | 0.374; 0.238 | 0.407; 0.221 | 0.422; 0.208 | 0.374; 0.144 | 0.212; 0.050 |

The 10° bins at 45–65N have their highest pressure at 45W (1016–1026 hPa)
and the 45–50N band its lowest at 5W (1000 hPa); the Aleutian bin
(60–65N, 175W) holds 984.8 hPa. Over the ranges (turbulent stress
magnitude, N/m²; lowest-layer speed, m/s): the Rockies (30–60N,
125–100W, above 1000 m) 0.301, 4.6; the Andes (55S–10N, above 1000 m)
0.128, 2.6; the Himalaya and Tibet (25–45N, 70–105E, above 2000 m) 0.186,
5.8; Greenland 0.236, 9.9; Antarctica below 2500 m 0.137, 7.5. Nine64 day
91 + 10 (21 June–1 July): the southern winter jet at 200 hPa 44.2 m/s at
40–45S, 40.4 at 35–40S, 29.4 at 30–35S; the lowest-layer westerly 7.5 m/s
at 45–50S; the northern summer 200 hPa maximum 22.0 at 55–60N. References:
ERA5 (Hersbach et al. 2020) northern-winter zonal means, the surface
westerlies over land near 3–5 m/s at 45–55N and the subtropical jet near
40 m/s at 200 hPa and 30N; the NCEP–NCAR reanalysis's DJF means (Kalnay
et al. 1996), the Aleutian low near 1000 hPa, the Icelandic low near
997 hPa, the Siberian high near 1035 hPa at 50N 100E and the zonal-mean
SLP near 1010 hPa at 60–70N. The signature of missing orographic drag is
too strong and too zonal northern midlatitude westerlies with a too deep
polar low (Palmer et al. 1986, McFarlane 1987). The state shows half of
it: the polar low is 9–13 hPa too deep in the zonal mean at 60–75N and
the Aleutian low about 15 hPa too deep, but the westerlies over land at
45–55N are 2.0–3.0 m/s, at or below the reference, and the 200 hPa jet is
not too strong but too far poleward (19.7 m/s at 30–35N, its maximum
33.0 at 50–55N); the stationary pattern is out of phase (no Siberian high,
a high at 45W). The states come from a spin-up with older physics.

The fields. Per cell, from the 0.25° raster clamped at sea level less the
resolved orography (the smoothed surface geopotential over g,
interpolated linearly on the triangle of the three nearest cell centres
to every raster point): μ² = ⟨h²⟩ − ⟨h⟩², γ² = (K − √(L² + M²))/(K +
√(L² + M²)), θ = ½ atan2(M, L) and σ² = K + √(L² + M²), with K, L, M
from the raster's central-difference gradients, area-weighted over the
points nearest the cell (Baines and Palmer 1990; IFS Cy47r3 Part IV
§11.3.4, where the source is 5 km data less the target orography). Over
land, μ below 50 / 50–100 / 100–200 / 200–400 / above 400 m on 0.44 /
0.19 / 0.19 / 0.14 / 0.05 of it at N=64 (land-mean μ 118 m, σ 0.0043,
γ 0.53, 30 raster points a cell) and 0.62 / 0.16 / 0.14 / 0.07 / 0.02 at
N=128 (μ 68 m, σ 0.0037, γ 0.42, 7.4 points a cell, 0.04 of the land
with fewer than 4); 0.19 s to compute at either N, on both engines at
the model's start. A cell the land mask makes sea has no subgrid
orography (at N=64, 2608 of the 29080 sea cells hold a μ from their land
points and the resolved terrain interpolated across the coast: mean 55 m,
883 above 50 m, 861 m on the cell holding Hawaii's Big Island; their drag
would slow air over a sea surface whose ocean never receives it). The raster's values are 0.25° area means, so it holds
the scales between its 28 km spacing and the cell (112 km at N=64, 56 km
at N=128), and loses everything below 28 km. Its land structure function
over 56–222 km rises as r^1.07 north–south and r^0.87 east–west within
30° of the equator, a spectral slope of −2.07 and −1.87 against Beljaars
et al.'s (2004) −1.9; on that slope the raster holds 76 % (N=64) and
53 % (N=128) of the variance between the cell and 5 km, the IFS's lower
bound for these fields, so μ is ×0.87 and ×0.73 of the IFS's, and 12 % and
9 % of the mean-square slope, so σ is ×0.35 and ×0.30. At N=128 a cell's
statistics come from about 7 points and the gradients from the points'
own neighbours, across the cell's edge; nothing is extrapolated.

The scheme, Lott and Miller (1997) as IFS Cy47r3 Part IV Chapter 4
documents it, on the edges after the boundary layer's momentum mixing,
from the state after the physics: the incident wind, density and N
(N² = g Δθ/(θ̄ Δz) between midpoints, eq. 4.26) as layer-mass means over
μ < z < 2μ (4.27); the blocking height Z_b the highest level below 3μ
where ∫ N/U_p dz from it to 3μ reaches H_n,crit, a level where U_p ≤ 0
blocking all below it (4.9, 4.32); below Z_b, ∂u/∂t = −C_d max(2 − 1/r,
0) (σ/2μ) √((Z_b − z)/(z + μ)) (B cos²ψ + C sin²ψ) |U| u/2 (4.14, 4.40),
implicit with |U| at the step's start (4.41); above it the wave stress
τ_0 = ρ_L (H_eff²/9) (σ/μ) G |U_L| √(D1² + D2²) N_L, H_eff = 3μ − Z_b
(4.37; for Z_b = 0 LM97's form with H = 2μ), along (D1, D2) in the frame
of U_L (4.31), constant to Z_b, then cut level by level wherever the
wave Richardson number N²(1 − α)/(S + Nα)² (S the shear of the wind in
the stress's plane, α = N δz/V, ρ N V δz² ∝ τ) falls below Ri_crit, all
of it at a critical level V ≤ 0; breaking below Z_b + Δz (∫ N/U_p dz =
π/2 above Z_b, at least 4μ) spreads linearly in pressure over that depth
(4.33, 4.35), the weight clamped to [0, 1] because p(Z_b) is interpolated
in ln p between midpoints while the interfaces sit at the midpoints' mean
height (up to 80 m apart in the lowest 2 km of bl34); the stress left at
the top goes into the top layer (σ < 0.0022; the IFS spreads what is
left above 9.9 Pa between there and its top, Cy33r1 eq. 4.39).
B = 1 − 0.18γ − 0.04γ², C = 0.48γ + 0.3γ² (Phillips 1984). Constants,
LM97's: C_d 1, G 1, H_n,crit 0.5, Ri_crit 0.25. The IFS documents
H_n,crit 0.5 and Ri_crit 0.25 too, but C_d 2 (eq. 4.17) and H_eff doubled
since Cy32r2 (eq. 4.8), with its own G that the chapter does not state
(G ≈ 1.23 for the elliptical mountain of eq. 4.2); those are not taken.
Each layer's wave drag in a step is at most what stops its own wind along
the stress, the rest passing to the layer above, a limit the IFS's
explicit tendencies do not have. An edge steps as u ← (u + Δt a)/(1 + Δt β); the kinetic energy
removed goes to the dissipation heat, and the momentum the column loses
is a stress on the ground per edge (`stress`, PH_OSTRESS), which the
ocean does not receive. Option `orography` (false off): `blockingDrag`,
`waveDrag`, `criticalHeight`, `criticalRichardson`.

Checks: over a uniform 10 m/s, N 0.01/s flow normal to a ridge (μ 200 m,
σ 0.02, γ 0) the column gives Z_b 100 m, τ_0 1/3 N/m² (ρ 1.2), the
blocking rate at each midpoint and the stress profile τ_0 min(1, (ρ/ρ_L)
(α_c/α_0)², α_c = (√2 − 1)/(2 Ri_crit)) by hand, to 2·10⁻³ of τ_0 (N from
finite differences); an oblique flow over γ 0.5 turns the stress towards
the cross-ridge axis by (D1, D2); after three N=8 steps one application
of the drag conserves each edge's column momentum against its stress to
10⁻¹² and returns the kinetic energy to the heat to 10⁻¹²; the engines agree to 2·10⁻⁴ of the largest
value (N=16, one step, real topography).

Over land at day 186 (eight64 + 3, day means): blocking on 0.33 of the
land, Z_b 447 m where it blocks; the launched wave stress 0.0135 N/m² in
the land mean, all of it taken in the column, 0.0112 below 500 hPa,
0.0013 at 500–100, 0.0009 at 100–10 and 0.0001 above 10 hPa: with H_n
near one most waves break just above the blocked layer (IFS §4.2.2);
the blocked flow's drag 0.014 N/m².

The surface layer. (1) Gustiness: U² = |v|² + u_g², u_g = max(0.2 m/s,
β w*) where B₀ > 0 and 0.2 m/s otherwise, w*³ = B₀ z_i with B₀ = −Ri_b U³
C_H/z the surface buoyancy flux of the step's own coefficients and z_i
the boundary layer's depth of the step before above the lowest midpoint
(at least z), found with the coefficients in four fixed-point passes as COARE
iterates its gust; β 1.2 and u_g 0.2 m/s in stable air over sea and sea
ice (COARE 3.5: Fairall et al. 1996, 2003, Edson et al. 2013; its code,
coare35vn.m and coare36vn, takes β w* alone where B₀ > 0, which drops
the gust from 0.2 m/s to 0 as B₀ crosses zero, so the floor is kept
there to make u_g continuous), β 1 over land (IFS eqs. 3.19–3.20, Beljaars 1994, with the model's z_i for
the IFS's 1000 m), in place of max(|v|, 3 m/s), for every surface flux,
the implicit drag and u*. A gust from the step before's flux went NaN in
seven hourly N=6 steps of a fresh start: a forest's C_H U grows without
bound as U → 0 in Dyer–Hicks free convection (z/z₀ 10: 0.23 m/s at 3 m/s
and 3.96 at 0.2 m/s for a 33 K contrast), and the in-step fixed point
holds it to 0.23–0.93. (2) Over land the bulk Richardson number takes the
humidity the evaporation implies, q₁ + w max(0, q_s(T_s) − q₁), w the
land's wetness from the step before's coefficients (IFS eq. 3.27 takes
the surface's own q), in place of q₁. Options (`SURFACE`):
`convectiveGust` [1.2, 1, 0.2] (false: the floor), `gustIterations` 4,
`landHumidity` 'wetness' or 'air'. With both off and `orography` false
the CPU reproduces 3f57d8f bit for bit (six N=8 steps).

Their effect, eight64 day 183 + 3 (8-step samples; U the speed the
fluxes take, C_H U in mm/s from ρ C_H U/1.2):

| samples | before | after | gust off | humidity 'air' |
|---|---|---|---|---|
| calm tropical sea, \|v\| < 3 (0.28 of 20S–20N): U m/s; C_H ×10⁻³; C_H U; LE W/m² | 3.00; 1.21; 3.54; 37.7 | 1.86; 1.51; 2.54; 27.9 | 3.00; 1.21; 3.55; 37.9 | 1.86; 1.50; 2.54; 27.9 |
| \|v\| < 1.5 (0.10) | 3.00; 1.20; 3.51; 35.4 | 1.02; 1.86; 1.73; 18.3 | 3.00; 1.20; 3.52; 35.6 | 1.02; 1.86; 1.73; 18.3 |
| 20S–20N open sea, LE W/m² | 74.8 | 72.2 | 74.8 | 72.3 |
| forest by day: U; C_H; C_H U; H; LE | 3.42; 25.9; 82.9; 35.1; 108.0 | 2.73; 32.1; 74.7; 35.2; 107.5 | 3.38; 26.9; 85.0; 35.1; 108.0 | 2.74; 29.9; 71.2; 35.4; 106.8 |
| forest by night | 3.56; 18.5; 65.1; −14.4; 35.1 | 2.85; 15.9; 48.8; −9.1; 27.4 | 3.51; 19.2; 66.1; −14.6; 35.4 | 2.92; 14.6; 47.6; −8.7; 26.8 |
| forest by night, \|v\| < 3 | 3.00; 17.4; 49.6; −13.9; 38.6 | 1.90; 13.3; 24.1; −5.2; 26.5 | 3.00; 18.4; 52.3; −14.4; 39.3 | 1.99; 11.2; 22.2; −4.3; 24.6 |

The calm sea's gust is 1.2 w* of its weak buoyancy flux, about 0.4 m/s,
where the floor gave 3 m/s; the transpiring forest's surface humidity
raises its daytime C_H by 7 % against the air's.

Runs, GPU, OCEAN `{"everySteps":8}`, from copies of the states, 3f57d8f →
these commits:

| outcome | nine64 day 274 + 10 | nine64 day 91 + 10 |
|---|---|---|
| orographic stress over land, N/m²: Rockies; Andes; Himalaya–Tibet; Greenland; Antarctica < 2500 m; land mean; land 45–55N | 0.109; 0.068; 0.111; 0.094; 0.056; 0.035; 0.052 | 0.022; 0.151; 0.061; 0.043; 0.123; 0.025; 0.014 |
| turbulent stress over land, N/m²: Rockies; Himalaya–Tibet; land mean; land 45–55N | 0.301 → 0.248; 0.186 → 0.149; 0.219 → 0.198; 0.390 → 0.349 | 0.171 → 0.164; 0.134 → 0.116; 0.194 → 0.182; 0.245 → 0.239 |
| lowest-layer speed, m/s: Rockies; Himalaya–Tibet; Andes | 4.6 → 4.2; 5.8 → 5.0; 2.6 → 2.5 | 2.7 → 2.7; 4.1 → 3.7; 3.2 → 2.9 |
| u lowest layer 45–50N / 50–55N, all; land, m/s | 4.0 / 3.5 → 3.8 / 3.6; 3.0 / 2.0 → 2.7 / 1.9 | 1.6 / 2.8 → 1.8 / 2.8; 1.1 / 1.5 → 1.0 / 1.5 |
| u 850 hPa 45–50N; u 200 hPa 50–55N, m/s | 8.9 → 8.5; 33.0 → 32.9 | 3.4 → 3.6; 19.8 → 20.1 |
| SH: u lowest layer 50–55S; u 200 hPa 40–45S, m/s | 7.1 → 7.0; 27.1 → 27.2 | 6.4 → 5.9; 44.2 → 44.5 |
| SLP zonal mean 60–65N; 70–75N; 80–85N, hPa | 1001.0 → 1001.7; 996.9 → 998.2; 1001.1 → 1002.7 | 1003.4 → 1003.6; 998.2 → 999.4; 1001.0 → 1001.4 |
| SLP stationary rms 45–50N; 60–65N; Aleutian bin, hPa | 6.2 → 5.8; 7.4 → 8.0; 984.8 → 983.4 | 9.9 → 9.9; 4.7 → 4.7; — |
| SH SLP 55–50S rms; wave 2, hPa | 4.5 → 4.5; 1.4 → 1.4 | 6.9 → 6.0; 4.7 → 2.0 |

With the gust and humidity off (the orography alone) the December run
gives the same zonal means to 0.43 m/s and 0.28 hPa in every band (land
45–55N turbulent 0.357, orographic 0.051; SLP 60–65N 1001.7).

| outcome | eight64 day 183 + 3 | eight128 day 183 + 1 |
|---|---|---|
| global evaporation (samples), mm/d | 2.000 → 1.948 (gust off 2.000, humidity 'air' 1.945) | 2.087 → 2.042 |
| global sensible heat, W/m² | 9.44 → 9.92 | 10.90 → 11.35 |
| sea stress magnitude, N/m² | 0.090 → 0.089 | 0.086 → 0.086 |
| equatorial Pacific τx 2S–2N 160E–100W, N/m² | −0.019 → −0.018 | −0.022 → −0.022 |
| Southern Ocean τx 40–60S mean; peak band, N/m² | 0.137 → 0.137; 0.191 → 0.191 at 47.5S | 0.102 → 0.102; 0.211 → 0.211 at 57.5S |
| global surface; land skin, last day, °C | 16.022 → 16.048; 13.91 → 14.01 | 16.251 → 16.261; 14.92 → 14.95 |
| forest night H; LE, W/m² | −14.4 → −9.1; 35.1 → 27.4 | −16.2 → −11.1; 41.3 → 34.5 |
| calm tropical sea LE, W/m² | 37.7 → 27.9 | 44.6 → 33.8 |

The audit of the end states (`scripts/verticalAudit.mjs`, each on its
own code): day 186 global rain 2.03 → 1.97 and evaporation 2.09 → 2.04
mm/d, the Pacific ITCZ's rain 5.27 → 4.71 mm/d, SE Pacific low cloud
0.268 → 0.263, the zonal rain peak 6.35 → 6.61 mm/d, cloud effects
−55.5 → −55.0 and 17.4 → 17.4 W/m²; N=128 day 184 global rain 0.92 →
0.92, evaporation 2.05 → 2.02, SE Pacific low cloud 0.132 → 0.133, the
zonal peak 3.31 → 3.38 mm/d, cloud effects −44.2 → −44.0 and 16.7 →
16.7 W/m².

The turbulent orographic form drag (Beljaars et al. 2004, IFS eqs.
3.54–3.57) needs σ_flt, the standard deviation of the orography in the
3–22 km band (the IFS's 30″ field filtered with Δ 2 km and 20 km, §11.3.3),
from 1 km topography: GMTED2010 (Danielson and Gesch 2011, USGS, public
domain; its 30″ mean product 43200 × 21600 16-bit values, about 1.9 GB),
SRTM30_PLUS (Becker et al. 2009) or GTOPO30 (USGS 1996, public domain,
the same size); the per-cell σ_flt it gives is 4 bytes a cell (0.66 MB
at N=128). Not fetched. On the measured slope the 0.25° band gives σ_flt
= μ (k_flt^n I_H (−n − 1)/(k_c^{n+1} − k_s^{n+1}))^½ (IFS eqs.
11.11–11.14, k_flt 0.00035 m⁻¹, I_H 0.00102 m⁻¹, k_c and k_s the cell's
and the raster's π/Δ): 0.58–0.62 μ at N=64 and 0.99–1.04 μ at N=128,
land means 68–73 m and 67–71 m, and a TOFD stress at 10 m/s from the
lowest level up of 0.59–0.68 N/m² (N=64) and 0.67–0.74 (N=128) in the
land mean, 0.6–0.7 off the ice sheets: three times the vegetation's
turbulent stress, from an extrapolation over a factor 5–10 in scale on
an assumed power law; an order of magnitude for the land as a whole, not
a cell's value. It is not built.

Step cost under the exclusive lock (`profileGpu`, each code on its copy of
eight64 day 183 and eight128 day 183, 128 steps after 16, twice): N=64
20.44–20.48 → 21.34 ms a step (5.24 → 5.46 s a model day), the physics,
boundary-layer and orography kernels 2.64 → 3.42 ms; N=128 89.17–89.36 →
92.06–92.13 ms (45.7 → 47.1 s a model day), those kernels 11.0 → 13.4 ms.
The orography and the gust add 2KC + 9C + E floats to the GPU's physics
buffer (52 MB at N=128).

What still misses: σ from the 0.25° raster is ×0.30–0.35 of the 5 km
data's, so both stresses, each ∝ σ, are low by about that factor; the
TOFD (above); C_d and G are LM97's, not the IFS's calibration; the gust
has no deep-convective (mesoscale downdraft) part (Redelsperger et al.
2000), so the calm tropical sea's evaporation fell by 26 % where
convection would gust it; z_i is the step before's diagnosed depth, so a
column that switches its boundary-layer regime switches its gust; the
land humidity takes the step before's coefficients in its stomatal
factor; the scheme's stress does not reach the ocean (correctly) and its
dissipation heats the layer where the momentum goes; in ten days of
December the zonal-mean wind at 30–70N moved by at most 0.45 m/s at the
lowest layer, 0.65 at 850 hPa and 0.44 at 200 hPa and the SLP by at most
1.9 hPa, so the signature's answer to the drag needs a season, which
needs a spin-up.

The review of the mountains' drag (Oct 2), on 3066df7 and the fixes after it:

- Fields: recomputed for four N=64 cells by a separate script (brute-force
  nearest cells, the resolved terrain on the triangle that contains each
  point, its own differences): the West Siberian plain (61N 75E: 32
  points, μ 10.76 m, γ 0.631, θ 13.0°, σ 5.24·10⁻⁴), the Great Plains
  (39N 99W: μ 20.4 m, σ 7.79·10⁻⁴), the Himalayan front (28N 84E: 17
  points, μ 1335.5 m, γ 0.445, θ 85.6°, σ 0.0350) and the Andes (32.6S
  70W: 20 points, μ 922.7 m, γ 0.238, θ 7.9°, σ 0.0276), all equal to
  `subgridOrography`'s to the digits printed. The land bins, means and
  point counts above reproduce. Sea cells held fields (above); they no
  longer do.
- The column against hand values, a stated mountain (μ 300 m, γ 0.5,
  θ 30°, σ 0.015) in a uniform 10 m/s eastward flow with N 0.01/s on a
  hydrostatic bl34 column: Z_b 399.9999 m (hand 400), τ_0 0.123060 N/m²
  (hand ρ_L × 0.110441 with ρ_L 1.1143), the direction 18.5759° (hand
  atan(D2/D1)), the blocking rate at the six midpoints below Z_b to
  6·10⁻⁷, the stress profile against τ_0 min(1, (ρ/ρ_L)(α_c/α_0)²) to
  5·10⁻⁶ of τ_0, all of it taken by the top. On a synthetic N=16 state
  with those fields everywhere both engines give the same Z_b and τ_0 to
  the digits printed (6 % from the hand values, from the cell
  reconstruction of the edge winds at N=16).
- The stress never grew upward in uniform columns, but over a wind that
  falls with height 4376 of 142104 swept columns had a layer the wave
  drag accelerated, by up to 2.6 % of τ_0 (the low-level breaking's
  weight, above); with the clamp none does.
- Momentum: each edge's column loss against its stress to 1.7·10⁻¹⁶
  relative in ten CPU steps from nine64 day 274 and from eight128 day 183;
  no sea–sea edge takes any. Of the launched stress 0.9995–0.9996 (N=64)
  and 0.9972–0.9977 (N=128) is taken in the column; the rest, in 277–342 and
  774–957 columns, is what the per-step limit keeps from the top layers.
- Stability over those ten steps: the largest blocking rate times the step
  0.20 (N=64, 85.7S 151.5W, μ 422 m, σ 0.038) and 0.21 (N=128, 86.1S
  156.3W), solved implicitly; the largest wave tendency times the step
  2.9 m/s (N=64) and 1.0 m/s (N=128); the largest edge change from the
  drag in a step 2.1 and 2.8 m/s; everything finite, the largest edge wind
  93.6 and 84.0 m/s.
- Engines on nine64 day 274 (N=64, real state, no ocean): after 1, 4 and
  16 steps the launched stress apart by 7.5·10⁻⁵, 2.0·10⁻² and 1.1·10⁻²
  of its largest value, the blocking height by 2.7·10⁻⁴, 2.0·10⁻² and
  6.4·10⁻², the edge stress by 1.2·10⁻³, 9.0·10⁻³ and 8.8·10⁻³, the gusty
  wind by 4·10⁻⁶, 9.5·10⁻³ and 8.1·10⁻²; the state's u apart by
  2.8·10⁻³, 3.2·10⁻³ and 6.1·10⁻³ m/s rms, against 2.8·10⁻³, 4.1·10⁻³ and
  8.9·10⁻³ with the scheme off (the largest single edge 1.8 m/s after one
  step either way).
- With everything off the CPU reproduces 3f57d8f bit for bit (N=8 six
  steps, N=16 eight steps), on 3066df7 and on the fixes.
- The gust: u_g continuous across B₀ = 0 (above). z_i is the depth above
  the lowest midpoint, not the ground (z, 21 m in bl34, short: about 1 %
  of w* at z_i 500 m); measuring it from the ground is right but was not
  kept, because the two 48-step f32/f64 parity tests turn on it by luck:
  with z_i from the ground the treeline test passes under the default
  surface layer (tree cover 8.5·10⁻⁵ apart against 1.8·10⁻⁴, the lowest
  air 1.8·10⁻³ K against 9.6·10⁻³) and the snow-albedo test fails (rms
  3.3·10⁻⁴ against 4.4·10⁻⁵, its tolerance 10⁻⁴); scaling z_i by 0.99,
  1.01 and 1.05 gives that test 1.7·10⁻⁵, 2.1·10⁻⁴ and 2.2·10⁻⁴, and
  with the gust off it gives 1.0·10⁻⁴. The source is the gust's z_i: the
  boundary layer's diagnosed depth jumps where a column changes regime,
  so either engine's regime flip moves the surface fluxes. The treeline
  test keeps 3066df7's premise (gust and land humidity off). The land
  humidity in Ri_b is
  continuous in w and T_s and the same on both engines. The exchange
  flux parity's relative measure stands: the fresh test state's snow
  cells carry |LE| near 2·10⁴ W/m² on 3f57d8f as well (19623 W/m²).
- Sources read for this review: IFS Cy33r1 Part IV Chapter 4 confirms
  H_n,crit 0.5 (§4.1), Ri_crit 0.25 (§4.4.2), H_eff = 2(H − Z_blk) since
  Cy32r2 (eq. 4.8), G ≈ 1.23 for the mountain of eq. 4.2, the /9 of eq.
  4.37 with H = 3μ, eqs. 4.26–4.41 as used, and gives C_d only as close to
  1 by free-streamline theory and nearer 2 with suction behind the
  obstacle; C_d 2 as Cy47r3's value and LM97's C_d 1 and G 1 were not
  checked against those texts. Its gust is |U|² = u² + v² + w*² with
  z_i 1000 m (Cy33r1 eqs. 3.17–3.18) and its Ri_b takes q_surf (eqs.
  3.24–3.25). COARE's gust from coare35vn.m and coare36vn as above.
- Ten days from nine64 day 274 on 3066df7 reproduce the run above to the
  digits printed. A twin with θ perturbed by 10⁻⁷ (relative, random per
  value) gives, against the documented 3f57d8f → 3066df7 changes over the
  5° bands at 30–85N (rms; largest): SLP 0.15 against 0.98 hPa (0.30;
  1.79), the lowest wind over land 0.06 against 0.41 m/s (0.12; 0.90), the
  lowest wind 0.14 against 0.26 (0.31; 0.55), 850 hPa 0.18 against 0.35
  (0.36; 0.68), 500 hPa 0.25 against 0.41 (0.60; 0.81), 200 hPa 0.22
  against 0.25 (0.64; 0.44); the stationary wave's rms at 55–65N 0.25–0.31
  against 0.17–0.59 hPa and the Aleutian bin 0.3 against 1.4 hPa. The
  SLP's rise at 60–85N and the slower land winds are 5–7 times the twin's
  spread; the 200 hPa jet's and the stationary waves' changes are within
  it. The fixed code against 3066df7 over the same ten days differs by
  about the twin's spread (rms SLP 0.15 against 0.15 hPa, the lowest
  wind over land 0.07 against 0.06, 500 hPa 0.34 against 0.25, 200 hPa
  0.20 against 0.22 m/s; the orographic stress over land 0.052 N/m² at
  45–55N in both, 0.034 against 0.035 in the land mean).

**The model top and the mountains merged (Oct 2).** Branch integrate-c:
gas-benchmark 9b7b476, then roughness bd42fa5. The closure (∇⁴,
divergence damping, sponge) sets the dissipation; the implicit drag
with the boundary layer's mixing, the plume's transport, the orographic
drag and the gravity waves add to it, in that order; the mountains'
stress is its own array. Fixes after the merges: the GPU's prior land
wetness takes TRACE_SNOW; the spin-up saves `exchangeHeat` and
`exchangeWind`; four engine-parity tests run with `convectiveGust`
false. The twelve-step digests re-pinned; with `topDragDays` 5,
`spongeDays` 0 and `gravityWaves` false they are 32417c2's. Tests pass
but the cloudy columns' heating (1.26·10⁻⁴ K/day) and the cloud-effect
parity (one column 0.66 against 0.5 W/m² at step 8; passes with
`spongeDays` 0, but with `spongeDays` 0 fails the same way, 0.665
W/m², from initial θ perturbed by 3·10⁻⁷ or 10⁻⁶ relative: a column at
a threshold, not the sponge). Proofs, nine64_day0274, the second CPU
step after loading: each edge
column's momentum + Δt τ 1.9·10⁻¹¹ (implicit drag) and 2.3·10⁻¹³
(mountains) of 1.6·10³ kg/m/s, each cell column's wave force 6.5·10⁻¹⁹
of 4.2·10⁻³ Pa, the energy each gives the heat against its own loss
≤ 2.9·10⁻¹¹ per column; the sponge 3.8·10¹⁵ and 1.1·10¹⁵ J at 1.1
and 3.6 hPa, its angular momentum change 0.033 and 0.020 of a Rayleigh
drag's at its rate; heat applied 9.31·10¹⁷ J against closure 1.43·10¹⁷
+ drag 7.45·10¹⁷ + mountains 4.27·10¹⁶ + waves 5.3·10¹⁴ to 9.6·10⁻¹⁵;
no sea–sea edge touched; the ocean's stress is the implicit drag's
(0 apart). Each GPU kernel: 2.1·10⁻⁵ (drag), 1.6·10⁻⁵ (mountains),
4.2·10⁻⁸ (waves) relative, heat against the sinks 6.5·10⁻⁴. With
`topDragDays` 5, `spongeDays` 0, `gravityWaves` false, `orography`
false, `convectiveGust` false, `landHumidity` 'air' the CPU gives
32417c2's state, land and ocean digests after each of nine steps from
eight64_day0183. Parity from nine64_day0091 after 1, 4, 16 steps: top six
layers' u 7.4·10⁻⁵, 1.7·10⁻⁴, 3.4·10⁻⁴ m/s rms; lowest layer T 1.0·10⁻⁴,
4.0·10⁻⁴, 1.9·10⁻³ K, u 4.9·10⁻⁵, 2.5·10⁻⁴, 1.5·10⁻² m/s (gust off
4.9·10⁻⁵, 5.5·10⁻⁴, 1.6·10⁻³), Ts 9.4·10⁻⁶, 2.6·10⁻⁵, 1.3·10⁻³ K;
regime flips 0, 2, 27; launched stress 1.0·10⁻⁴, 9.6·10⁻⁴, 1.5·10⁻³ of
its largest; wave acceleration 2.7·10⁻⁴, 9.0·10⁻⁴, 1.1·10⁻³ m/s/day rms.
The lowest layer's largest edge apart at step 16, 2.84 m/s, is a sea
edge at 65.8S 43.3E whose cell's boundary layer is in regime 0 on the
CPU and 3 on the GPU (depth 829 against 36 m) under the same gust wind
(27.7 m/s); the gust wind is 9.3·10⁻⁴ m/s rms apart. From the same
state bd42fa5 gives 1.95·10⁻³, 2.05·10⁻³, 8.2·10⁻³ m/s and 32417c2
4.9·10⁻⁵, 5.5·10⁻⁴, 1.0·10⁻³ m/s for the lowest layer's u.
Three one-day segments from eight64_day0183 end day 186 byte for byte as
one. 16 steps from eight128_day0183: finite, top six layers' largest
wind 73 m/s, Courant 0.24 / 0.04. five64_day2190 (cam26) loads and steps.

Three GPU days (32417c2 in brackets): eight64 day 186 albedo 0.293
(0.296), ASR 240.9 (239.8), OLR 233.8 (233.6), SWCRE −52.8 (−53.9),
LWCRE 26.1 (26.2), rain 1.72 (1.76); ten64 day 186 0.310 (0.312), 235.0
(234.4), 232.6 (232.5), −58.5 (−59.1), 26.0 (26.1), 2.53 (2.54); nine64
day 94 0.305 (0.307), 236.8 (235.9), 234.3 (234.1), −57.6 (−58.5), 26.0
(26.1), 2.36 (2.39). eight64, 8-step samples: evaporation 2.219 (2.270;
bd42fa5 1.962) mm/d, sensible 12.01 (11.59) W/m²; calm tropical sea
(|v| < 3 m/s, 0.25 of 20S–20N sea) U 1.92 (3.00) m/s, C_H 1.68 (1.33)
·10⁻³, LE 45.1 (58.0; bd42fa5 28.0) W/m². Pacific ITCZ, ten64 day 186 + 1
CPU day: rain 6.51 (6.48), convective 0.24 (0.24), lowest-layer RH 0.81
(0.81), T −1.3 / −2.8 / +1.0 (−1.3 / −2.8 / +1.1) K against Jordan at
1008 / 848 / 516 hPa, plume fired 0.167 (0.160), CAPE 169 (170) J/kg,
peak 973 (973) hPa. nine64 day 94, 70–90S at 1.1 / 3.5 / 7.4 / 14 / 24
hPa: 222.9 / 237.1 / 229.1 / 222.0 / 213.5 (224.5 / 240.9 / 234.0 /
226.5 / 216.8) K. Ten days from nine64_day0274, day 284: 70–90N at 1.1 /
3.5 / 7.4 / 14 / 24 / 37 / 53 hPa 220.7 / 225.2 / 218.9 / 213.1 / 208.1
/ 205.0 / 203.0 (221.9 / 229.8 / 226.5 / 221.0 / 214.2 / 209.3 / 206.1)
K; zonal-mean u at 65N, 1.1–37 hPa, 49.7 / 39.9 / 40.3 / 39.8 / 37.8 /
35.0 (51.0 / 36.1 / 34.9 / 35.1 / 34.8 / 33.2) m/s; at 1.1 hPa 69.7 m/s
at 35N (40.3), largest wind 86 (76) m/s; 8-step SLP 60–65 / 70–75 /
80–85N 1002.4 / 997.5 / 1002.2 (1001.5 / 996.2 / 1001.0; bd42fa5
1001.9 / 998.0 / 1002.1) hPa; orographic / turbulent stress over land,
N/m², Rockies 0.108 / 0.231 (bd42fa5 0.109 / 0.249), Andes 0.062 /
0.112 (0.068 / 0.113), Himalaya–Tibet 0.093 / 0.108 (0.114 / 0.150),
land 0.033 / 0.189 (0.034 / 0.199; 32417c2 turbulent 0.210), land
45–55N 0.052 / 0.354 (0.052 / 0.351; 0.399). Cost under the exclusive
lock, 128 steps after 16, twice: N=64 22.54 / 22.65 → 24.18 / 24.32 ms
(+7.3 %; physics pass 4.71 → 5.96, mixing pass 4.55 → 4.82); N=128
99.22 / 99.24 → 104.65 / 104.70 ms (+5.5 %; 20.09 → 24.47, 16.43 →
17.24); a day at N=128 1.1 min of steps on both.

**The terrain's fine scales, from data (Oct 2).** `scripts/subgridTerrain.py`
(numpy and Pillow: the 30″ filters need FFTs), `scripts/subgridTerrain.mjs`,
the per-mesh files `data/subgrid_N16.bin`, `_N32`, `_N64`, `_N128`, their
codec and loader in `js/geography.module.js` (`encodeSubgrid`,
`decodeSubgrid`, `meshSubgrid`), `js/physics/formDrag.module.js`, the form
drag in `boundaryLayer.module.js` and `physics.gpu.js`; tests in
`test/subgridTerrain.test.mjs`, `test/formDrag.test.mjs` and
`test/orography.test.mjs`.

The data. GMTED2010's 30-arc-second mean elevation (Danielson, J.J., and
Gesch, D.B., 2011, Global multi-resolution terrain elevation data 2010
(GMTED2010): U.S. Geological Survey Open-File Report 2011–1073, 26 p.), a
USGS product and so in the public domain (a work of the U.S. Government;
cite the report). The 108 tiles of 20° × 30° (`<lat><lon>_20101117_gmted_mea300.tif`,
lower-left corners 90S…70N and 180W…150E), fetched on Oct 2 2026 with
curl from
`https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/topo/downloads/GMTED/Global_tiles_GMTED/300darcsec/mea/<W180…E150>/`
into `/Users/jlhawn/git_repos/jlhawn/geodesic-terrain-cache` (outside every
repository, never committed): 17,339,044 to 17,339,080 bytes each,
1,872,620,024 in all, each equal to the server's Content-Length. The host
publishes no checksum for them; in its place each tile's own GDAL
statistics (minimum, maximum, mean, standard deviation of the valid values)
are reproduced from its data to 10⁻⁶ relative. Each tile is an
uncompressed GeoTIFF of 3600 × 2400 little-endian int16 (Pillow and a raw
read agree), WGS 84 (EPSG 4326), pixel scale 1/120°, RasterPixelIsArea, its
tie point at its named north-west corner less 0.5″ (the grid's half-arc-second
registration offset, 15 m, not carried further), nodata −32768 on the 720
rows north of 84N (31,104,000 points) and the sea at 0 m. Spot heights:
Everest 8625 m at 27.9874N 86.9290E (the summit 8849 m; a 30″ mean), K2
8190, Aconcagua 6678, Denali 5977, Mont Blanc 4447, the Dead Sea −412 m,
Vostok 3500 m (the station 3488), Greenland's Summit 3206 (3216), the plains
at 39N 99W 546 m, the Pacific at 0N 140W 0. `assemble` writes the global
21600 × 43200 grid (1.87 GB) into the cache in 7 s.

The fields, as IFS Cy47r3 Part IV §11.3 defines them, every smoothing
with its operator (eq. 11.4: a radial top hat of width Δ with cosine edges
of half-width δ 1 km), applied by FFT in latitude bands of at most 1° on
the 30″ grid, the sea and land below sea level at 0 m (`filter`, 77 s):
- μ, γ, θ, σ (§11.3.4): the 30″ orography smoothed at Δ 5 km and sampled on a
  2′30″ grid (4.6 km), less the orography the model resolves, its own
  surface geopotential over g (the 0.25° ETOPO1 raster's cell means, two
  smoothing passes on the mesh) interpolated linearly on the triangle of
  cell centres, so the fields hold the scales between 5 km and what the
  mesh carries; gradients by central differences 5 km each way (the
  spacing of the IFS's 5 km grid), bilinear between the 2′30″ points; K, L,
  M, h, h² area-weighted over the 2′30″ points nearest the cell
  (`subgridOrography` on this grid with `spacing` 5000). The IFS subtracts its
  own target-resolution smoothing; this model's resolved terrain and land
  mask are left as they are. Over the land off the ice sheets the
  residual's cell mean is 19.9 m at N=64 (rms 151 m) and 10.1 m at N=128
  (rms 117 m), mostly the mesh's two smoothing passes; μ removes it. At
  0.25° GMTED and the ETOPO1 raster differ by 1.9 m on average (rms 62 m)
  over land outside Antarctica, 6.5 m (rms 76 m) over Antarctica.
- σ_flt (§11.3.3): the 30″ orography smoothed at Δ 2 km and at Δ 20 km; the
  square of the difference averaged over each 2′30″ block (5 × 5 points),
  the blocks area-weighted over the cell, the square root.
Sea cells hold zeros (the land mask, as before). Files: a 16-byte header
and per field one uint16 per cell (μ in steps of 0.1 m, γ 1/65535, θ
π/65535 from −π/2, σ 10⁻⁵, σ_flt 0.05 m): 25,636 (N=16), 102,436 (N=32),
409,636 (N=64) and 1,638,436 bytes (N=128), 5–10 s each. The quantization's
worst error over land, relative with floors of 1 m, 0.01, 0.01 rad, 10⁻³
and 1 m: μ 0.3–5 %, γ ≤ 1.3·10⁻⁴, θ 2.4·10⁻³ rad, σ 0.5 %, σ_flt 1.1–2.5 %
(at the floors). Both engines and the page's worker (which fetches
`data/subgrid_N<N>.bin` for the run's and the device test's N) use them,
with any value on a sea cell of the run's land mask set to zero; a mesh
without a file (or whose fetched file does not decode), a run without
terrain (the files hold only the scales below the resolved terrain) and a
land mask on which more than 0.5 % of the land cells would hold nothing
compute the fields from the topography's raster, say why once in the log,
take Lott and Miller's constants with them and have no form drag.

| land | N=16 | N=32 | N=64 | N=128 |
|---|---|---|---|---|
| cells; 2′30″ points a cell | 731; 17148 | 2950; 4271 | 11882; 1069 | 47357; 267 |
| mean μ, m (0.25° raster) | 292.2 | 202.7 | 141.7 (118) | 100.1 (68) |
| mean σ (0.25° raster) | 0.0119 | 0.0117 | 0.0113 (0.0043) | 0.0109 (0.0037) |
| mean γ (0.25° raster) | 0.795 | 0.750 | 0.704 (0.53) | 0.651 (0.42) |
| μ < 50 / 50–100 / 100–200 / 200–400 / ≥ 400 m | 0.08 / 0.20 / 0.24 / 0.24 / 0.25 | 0.24 / 0.20 / 0.21 / 0.22 / 0.14 | 0.39 / 0.18 / 0.17 / 0.17 / 0.08 | 0.53 / 0.16 / 0.16 / 0.12 / 0.04 |
| σ_flt mean; rms; off the ice sheets, m | 57.6 | 56.0 | 53.8; 86.5; 57.8 | 51.7; 87.3; 55.8 |
| σ_flt < 10 / 10–25 / 25–50 / 50–100 / 100–200 / ≥ 200 m | | | 0.27 / 0.24 / 0.15 / 0.16 / 0.14 / 0.05 | 0.30 / 0.24 / 0.14 / 0.14 / 0.13 / 0.05 |

σ from the 5 km data is 2.62 (N=64) and 2.96 (N=128) times the 0.25°
raster's (the spectral estimate above gave 2.9 and 3.3), μ 1.20 and 1.48
times (1.15 and 1.37). The land-mean σ_flt is 54 and 52 m against 68–73 m
extrapolated from the 0.25° band.

Checks (`test/subgridTerrain.test.mjs`):
- Analytic orography through `filter` on a 5° band at the equator: a ridge
  of 300 m amplitude at 40 km keeps 207.45 m rms after the 5 km smoothing
  (the kernel's response, ∫ h J₀(kr) r dr / ∫ h r dr, gives 207.47) and
  58.39 m in the 3–22 km band (58.30); at 8 km 115.76 (115.95) and 173.18
  (172.58); an isotropic field sin kx sin ky at 40 km 143.96 (143.46) and
  74.48 (74.11).
- An analytic ridge and an isotropic field (amplitude 300 m, 400 waves round
  the equator) through the per-mesh script on N=16: μ within 2 % of A/√2 and
  A/2, σ within 3–4 % of the central-differenced slope, γ < 0.03 for the
  ridge and > 0.95 cos φ for the isotropic field, θ east–west, σ_flt exact.
- A north–south ridge of 20 km wavelength on 2′30″ rows (N=16): σ within
  1 % of A k F/√2, F = 0.594 the response of the bilinear 5 km difference
  (the neighbouring rows' would be 0.683).
- Three N=64 cells recomputed by `scripts/subgridTerrainHand.py` from the
  30″ grid by direct sums (each smoothing laid out at its own row's
  latitude, the cell's members by brute force against its two rings of
  neighbours, the resolved terrain on the containing triangle, the
  gradients 5 km each way bilinear between the 2′30″ points, σ_flt over
  the cell's own 30″ points as the IFS defines it), against the file:

| cell | μ, m | γ | θ, rad | σ | σ_flt, m |
|---|---|---|---|---|---|
| the Great Plains, 39.4N 98.7W (754 points; 18732 at 30″) | 30.376 / 30.4 | 0.5388 / 0.5388 | 1.3842 / 1.3842 | 0.003232 / 0.00323 | 15.01 / 15.00 |
| the Himalayan front, 28.5N 84.4E | 1421.24 / 1421.2 | 0.7991 / 0.7990 | 1.1338 / 1.1338 | 0.08725 / 0.08725 | 461.4 / 462.2 |
| the Andes, 32.7S 70.2W | 940.09 / 940.1 | 0.9272 / 0.9272 | 0.0603 / 0.0604 | 0.05927 / 0.05927 | 359.3 / 358.8 |

- The spectrum (`spectrum`; land segments of 1024 points wholly above 0 m
  within 60° of the equator, Hann window, linear trend removed; 918
  east–west and 6973 north–south segments): slope −1.85 east–west and
  −1.85 north–south over the filter's band (k 0.00014–0.00112 m⁻¹,
  wavelengths 45–5.6 km), against Beljaars et al.'s (2004) n₁ −1.9 that the
  IFS's I_H and k_flt assume (eq. 11.14); −1.68 and −1.69 over 63–10 km;
  −2.63 and −2.68 over k₀–k₁ (10–2.1 km), where the 30″ cell means'
  own averaging steepens the spectrum toward the grid scale. The fit
  recovers −1.91 from a synthetic −1.9 line.

The mountains' drag on these fields. The fields now follow the IFS's
definition (5 km data less the target orography), so the scheme takes the
constants IFS Cy47r3 Part IV Chapter 4 documents for them (read from that
document; Cy43r1's states the same C_d): C_d 2 (eq. 4.17), H_eff = 2(H − Z_blk)
(eq. 4.8, the factor added in Cy32r2 "because diagnostics indicated that
without the factor 2, the gravity wave activity was too weak"), G 1.23 (the
only value the chapter gives, for the elliptical mountain of eq. 4.2 that
the scheme assumes; it calls G "a function of the mountain sharpness" and
states no operational value), H_n,crit 0.5, Ri_crit 0.25
(`OROGRAPHY_DEFAULTS`; `LOTT_MILLER` keeps C_d 1, G 1, factor 1). The
doubled H_eff enters the wave's launch amplitude δz as well as its stress.
Eq. 4.37's remark that its /9 equals eq. 4.19's /4 at Z_blk = 0 holds only
without the factor; taken as eq. 4.8 states it. Against Lott and Miller's
set the launched stress is 4 × 1.23 = 4.9 times larger for one column, the
blocking twice; by hand for μ 200 m, σ 0.02, γ 0, 10 m/s and N 0.01/s,
τ₀ = 4 × 1.23/3 = 1.64 N/m² (ρ 1.2) against 1/3, and the waves, launched at
α = N δz/U = 1.0 against the critical 0.83, break just above the blocked
layer.

The turbulent orographic form drag (Beljaars, Brown and Wood 2004) as IFS
Cy47r3 Part IV §3.4 documents it: ∂U/∂t = −C_tofd(z) |U| U with C_tofd =
α β C_md C_corr 2.109 e^(−(z/1500)^1.5) a₂ z^(−1.2), a₂ = a₁ k₁^(n₁−n₂),
a₁ = σ_flt² (I_H k_flt^n₁)⁻¹ and the published α 35, β 1, C_md 0.005, C_corr
0.6, n₁ −1.9, n₂ −2.8, k₁ 0.003 m⁻¹, k_flt 0.00035 m⁻¹, I_H 0.00102 m⁻¹ (eqs.
3.55–3.57); C_tofd = 3.161·10⁻⁷ σ_flt² e^(−(z/1500)^1.5) z^(−1.2). It enters
the boundary layer's implicit edge solve with the surface drag, on every
layer the solve spans (σ > 0.5, the lowest 15 of bl34's 34, to 5.4 km):
C_tofd(z)|U| on each layer's diagonal, z the midpoint's height above the
model's ground and |U| the cell's wind in that layer at the diagnosis (the
IFS's |U| from the old time level, U implicit; fully implicit here where
the IFS weights 1.5), averaged onto the edge. The momentum it removes is a
stress on the ground per edge (`formStress`, on the GPU `PH_FSTRESS`), not
on the ocean, which takes the surface drag's stress alone; its kinetic
energy goes to the layers' dissipation heat with the solve's own; σ_flt is
zero over sea cells. Option `orography.formDrag` (false off, or the
constants). For σ_flt 100 m, C_tofd is 8.668·10⁻⁵, 1.505·10⁻⁶ and
7.412·10⁻⁸ m⁻¹ at 20 m, 500 m and 2 km, and ρ ∫ C_tofd U² dz from 20 m at
10 m/s is 0.565 N/m²; on the files' σ_flt (∝ σ_flt²) that is 0.42 N/m² in
the land mean at N=64 and 0.43 at N=128 (0.46–0.47 off the ice sheets),
against the 0.6–0.7 estimated from the 0.25° band.

Checks: a column of the solve with σ_flt 100 m on every cell, a uniform
10 m/s and no mixing or surface drag gives each land cell's rate as the
formula by hand times its wind to 10⁻⁹, none over sea cells, each edge
u/(1 + Δt C_tofd |U|) to 10⁻¹², the stress Σ m C_tofd |U| u and the heat
equal to the kinetic energy removed to 10⁻¹²; on a real N=16 state each
edge column loses (τ_surface + τ_form) Δt to 10⁻¹¹ and no sea–sea edge takes
any form stress; the engines lay the same rates, form stress, surface
stress and u after one N=16 step to 5.7·10⁻⁶, 5.7·10⁻⁶, 1.2·10⁻⁵ and
1.8·10⁻⁶ of their largest values. On nine64 day 274 (real state, no
ocean), the engines after 1, 4 and 16 steps, rms of the difference over
rms of the CPU's (largest difference over largest value): launched stress
2.4·10⁻³ (8.8·10⁻³), 1.8·10⁻² (9.8·10⁻²), 1.4·10⁻² (5.9·10⁻²); blocking
height 1.3·10⁻³ (1.8·10⁻²), 1.9·10⁻² (0.22), 2.0·10⁻² (0.23); orographic
stress 1.6·10⁻³, 6.1·10⁻³, 9.0·10⁻³; form stress 1.7·10⁻³, 2.3·10⁻³,
4.1·10⁻³; form rate 4.7·10⁻⁴, 6.2·10⁻³, 4.0·10⁻³; the state's u 8.2·10⁻³,
6.8·10⁻³ and 1.2·10⁻² m/s rms apart against 8.3·10⁻³, 7.8·10⁻³ and
1.3·10⁻² with the scheme off, θ 1.0–1.7·10⁻² K either way; the GPU's σ_flt,
μ and σ equal the file's to float32 in every cell and are zero on every sea
cell. Ten CPU steps from that state: the largest blocking rate times the
step 0.54 (50.0N 87.6E), the form drag's 9.0 (30.0N 94.5E, the lowest
layer), both implicit; the largest wave tendency times the step 6.3 m/s,
in the top layer (53.7N 159.1E), within the per-step limit; of the
launched stress 0.9997 taken in the column; the boundary layer's solve
balances each edge column against its two stresses to 3·10⁻¹¹ of its
momentum change and heats by the kinetic energy it removes to 6·10⁻¹³; the
largest change of an edge in that solve 19.5 m/s; every value finite, the
largest wind 96.1 m/s. From nine64 day 91 (the strongest lowest-layer wind
of the states, 37.2 m/s): 6.6 and 0.55, balance 1.5·10⁻¹⁰; ten GPU steps
from eight128 day 183: 6.3 (29.7S 70.2W) and 0.45, every value finite. A
mesh without a file (N=8, six CPU steps; the raster's fields on N=64, two
CPU steps from nine64 day 274), the scheme off and the terrain off (N=64,
two steps) reproduce bd42fa5 bit for bit.

Runs, N=64 GPU, OCEAN `{"everySteps":8}`, ten days from copies of the
states with 8-step samples, as the diagnosis above (stresses as
magnitudes of each sample's cell vector, the wave stress the launched
τ₀ along its direction, the blocking the orographic stress less it):
bd42fa5 (0.25° fields, Lott and Miller's constants) → the GMTED fields
with Lott and Miller's constants and no form drag → the IFS's constants →
with the form drag (these commits). The twin is the last with θ perturbed
by 10⁻⁷ (relative, random per value).

| nine64 day 274 + 10 | 30–35N | 45–50N | 50–55N | 55–60N | 60–65N | 70–75N | 80–85N |
|---|---|---|---|---|---|---|---|
| SLP zonal mean, hPa | 1021.0 → 1020.9 → 1020.3 → 1020.6 | 1012.9 → 1013.3 → 1014.7 → 1014.8 | 1008.1 → 1008.8 → 1011.2 → 1011.3 | 1004.6 → 1005.9 → 1008.9 → 1009.5 | 1001.7 → 1003.2 → 1006.5 → 1007.5 | 998.0 → 1000.0 → 1004.7 → 1004.0 | 1002.8 → 1005.0 → 1011.0 → 1010.1 |
| u lowest layer, land, m/s | 1.2 → 1.0 → 0.6 → 0.5 | 2.7 → 2.4 → 1.6 → 1.3 | 1.8 → 1.5 → 1.0 → 0.6 | 1.5 → 1.3 → 1.0 → 0.7 | 1.6 → 1.4 → 1.1 → 1.0 | 1.0 → 1.2 → 0.6 → 0.8 | 2.5 → 2.1 → 1.1 → 0.9 |
| u lowest layer, all | 1.0 → 0.9 → 0.8 → 0.6 | 3.8 → 3.6 → 2.7 → 2.8 | 3.6 → 3.2 → 2.4 → 2.3 | 2.3 → 2.2 → 2.2 → 1.8 | 1.9 → 1.6 → 1.6 → 1.5 | −0.1 → −0.0 → −0.8 → −0.1 | −2.3 → −2.8 → −2.8 → −3.0 |
| u 850 hPa | 4.2 → 4.0 → 3.8 → 3.6 | 8.6 → 8.3 → 6.7 → 6.9 | 8.2 → 7.5 → 6.2 → 6.0 | 5.9 → 5.7 → 5.5 → 4.8 | 4.6 → 4.2 → 3.9 → 3.8 | 0.5 → 0.6 → −0.5 → 0.6 | −2.0 → −2.1 → −2.0 → −2.3 |
| u 200 hPa | 19.6 → 19.4 → 19.1 → 19.0 | 28.4 → 28.3 → 27.7 → 27.6 | 33.0 → 32.9 → 32.5 → 32.7 | 32.1 → 31.4 → 31.1 → 31.2 | 24.1 → 23.3 → 22.8 → 22.5 | 7.6 → 7.9 → 6.7 → 7.5 | 1.7 → 2.7 → 2.6 → 1.9 |

The stationary SLP at 45–50N has its maximum at 45W (1025.0 hPa) in
bd42fa5 and with the fields alone (1024.2), at 95E (1025.5 and 1027.1 hPa)
with the IFS's constants: a Siberian high where there was none (NCEP–NCAR
DJF near 1035 hPa at 50N 100E). The 60–65N minimum (the Aleutian bins,
175W) 983.7 → 983.3 → 988.0 → 987.6 hPa. The southern bands (30–55S) move
by at most 1.1 hPa and 0.5 m/s.

| nine64 day 91 + 10 | 45–50N | 50–55N | 60–65N | 70–75N | 80–85N | 45–40S | 50–45S | 55–50S |
|---|---|---|---|---|---|---|---|---|
| SLP, hPa | 1013.1 → 1013.0 → 1012.8 → 1012.8 | 1010.5 → 1010.5 → 1011.0 → 1011.1 | 1003.5 → 1003.9 → 1004.1 → 1004.5 | 999.2 → 1000.1 → 1001.5 → 1002.1 | 1001.9 → 1002.7 → 1005.4 → 1005.6 | 1012.4 → 1012.2 → 1012.4 → 1012.3 | 1004.9 → 1004.6 → 1004.9 → 1004.9 | 997.2 → 996.9 → 997.4 → 997.6 |
| u lowest, land; all, m/s | 1.0 → 1.0 → 0.7 → 0.6; 1.8 → 1.8 → 1.5 → 1.3 | 1.5 → 1.4 → 1.3 → 1.1; 2.9 → 2.8 → 2.7 → 2.6 | 1.5 → 1.4 → 1.3 → 1.1; 2.5 → 2.3 → 2.1 → 2.1 | 0.2 → 0.1 → −0.2 → −0.4; −0.1 → −0.1 → −0.5 → −0.8 | 0.5 → 0.5 → −0.8 → −0.3 | 2.7 → 2.4 → 2.0 → 1.6; 6.5 → 6.6 → 6.5 → 6.4 | 3.9 → 3.5 → 2.9 → 2.0; 7.7 → 7.8 → 7.6 → 7.5 | 4.8 → 4.1 → 3.3 → 2.3; 5.8 → 5.7 → 5.4 → 5.2 |
| u 200 hPa | 12.8 → 12.7 → 12.7 → 12.5 | 19.9 → 20.3 → 20.1 → 19.8 | 19.4 → 18.7 → 18.7 → 18.9 | 7.7 → 7.3 → 7.1 → 7.5 | 1.5 → 1.8 → 0.7 → 0.8 | 44.6 → 44.7 → 44.6 → 44.4 | 38.7 → 38.7 → 38.2 → 38.1 | 28.8 → 28.7 → 28.4 → 28.3 |

Change against the twin's spread over the 5° bands 30–85N (rms; largest):

| | SLP, hPa | u lowest, land | u lowest | u 850 | u 500 | u 200, m/s |
|---|---|---|---|---|---|---|
| Dec: bd42fa5 → these commits | 4.71 (7.31) | 1.02 (1.62) | 0.78 (1.28) | 1.21 (2.24) | 1.12 (1.84) | 0.86 (1.56) |
| Dec: the fields alone | 1.34 (2.20) | 0.28 (0.45) | 0.27 (0.46) | 0.37 (0.68) | 0.51 (1.29) | 0.54 (0.98) |
| Dec: the form drag alone (against the IFS's constants without it) | 0.68 (1.20) | 0.26 (0.41) | 0.26 (0.64) | 0.47 (1.13) | 0.70 (1.45) | 0.38 (0.77) |
| Dec: twin | 0.15 (0.27) | 0.05 (0.09) | 0.12 (0.27) | 0.19 (0.41) | 0.19 (0.29) | 0.16 (0.31) |
| June: bd42fa5 → these commits | 1.95 (3.73) | 0.44 (0.81) | 0.40 (0.78) | 0.57 (0.96) | 0.61 (1.15) | 0.34 (0.72) |
| June: the fields alone | 0.51 (0.91) | 0.08 (0.15) | 0.15 (0.32) | 0.21 (0.47) | 0.38 (0.91) | 0.38 (0.76) |
| June: the form drag alone | 0.39 (0.87) | 0.21 (0.55) | 0.23 (0.64) | 0.28 (0.74) | 0.25 (0.49) | 0.29 (0.57) |
| June: twin | 0.10 (0.24) | 0.04 (0.11) | 0.03 (0.07) | 0.06 (0.12) | 0.10 (0.18) | 0.12 (0.23) |

The December SLP rise at 60–85N (6–7 hPa) is 31 times the twin's rms, the
land wind's fall 19 times, 200 hPa 5 times; the form drag alone moves the
zonal means by 2–7 times the twin's spread. The two December twins (the
fields of 9f1faf7 and of a8105bd) differ in spread by a factor 2.5 (SLP
0.37 and 0.15 hPa), so these ratios carry that much uncertainty.

Stresses over land, N/m² (bd42fa5 → these commits; lowest-layer speed
m/s; veg the turbulent surface stress):

| nine64 day 274 + 10 | speed | veg | form | blocking | waves | orographic | total |
|---|---|---|---|---|---|---|---|
| the Rockies | 4.22 → 2.38 | 0.254 → 0.086 | 0 → 0.214 | 0.060 → 0.265 | 0.061 → 0.358 | 0.111 → 0.579 | 0.363 → 0.874 |
| the Andes | 2.49 → 1.33 | 0.113 → 0.051 | 0 → 0.124 | 0.058 → 0.176 | 0.012 → 0.070 | 0.069 → 0.235 | 0.179 → 0.401 |
| the Himalaya and Tibet | 5.08 → 1.87 | 0.150 → 0.033 | 0 → 0.258 | 0.086 → 0.247 | 0.033 → 0.142 | 0.113 → 0.369 | 0.258 → 0.638 |
| the Alps (43–48N, 5–17E, above 500 m) | 3.51 → 1.18 | 0.337 → 0.031 | 0 → 0.299 | 0.123 → 0.359 | 0.107 → 0.320 | 0.211 → 0.617 | 0.541 → 0.938 |
| Greenland | 8.84 → 6.31 | 0.184 → 0.099 | 0 → 0.107 | 0.053 → 0.145 | 0.050 → 0.165 | 0.086 → 0.267 | 0.267 → 0.466 |
| Antarctica below 2500 m | 6.85 → 5.21 | 0.115 → 0.072 | 0 → 0.027 | 0.035 → 0.076 | 0.032 → 0.109 | 0.056 → 0.146 | 0.169 → 0.243 |
| land mean | 3.83 → 2.94 | 0.199 → 0.125 | 0 → 0.058 | 0.021 → 0.073 | 0.019 → 0.099 | 0.034 → 0.148 | 0.232 → 0.328 |
| land 45–55N | 4.33 → 3.11 | 0.351 → 0.192 | 0 → 0.092 | 0.025 → 0.106 | 0.035 → 0.176 | 0.052 → 0.248 | 0.402 → 0.528 |
| land mean, eastward | | 0.034 → 0.004 | 0 → 0.021 | 0.004 → 0.014 | 0.007 → 0.030 | | 0.045 → 0.069 |
| land 45–55N, eastward | | 0.201 → 0.064 | 0 → 0.061 | 0.011 → 0.036 | 0.027 → 0.107 | | 0.239 → 0.268 |
| **nine64 day 91 + 10** | | | | | | | |
| the Rockies | 2.67 → 1.79 | 0.163 → 0.090 | 0 → 0.100 | 0.014 → 0.068 | 0.010 → 0.073 | 0.022 → 0.125 | 0.184 → 0.313 |
| the Andes | 2.88 → 1.32 | 0.125 → 0.038 | 0 → 0.166 | 0.120 → 0.306 | 0.036 → 0.184 | 0.152 → 0.462 | 0.273 → 0.658 |
| the Himalaya and Tibet | 3.70 → 1.62 | 0.115 → 0.035 | 0 → 0.143 | 0.049 → 0.133 | 0.013 → 0.040 | 0.060 → 0.165 | 0.172 → 0.337 |
| land mean | 3.63 → 3.01 | 0.183 → 0.135 | 0 → 0.044 | 0.017 → 0.054 | 0.011 → 0.061 | 0.024 → 0.099 | 0.207 → 0.276 |
| **eight128 day 183 + 1** | | | | | | | |
| the Andes | 4.16 → 1.79 | 0.207 → 0.065 | 0 → 0.294 | 0.047 → 0.227 | 0.050 → 0.332 | 0.078 → 0.490 | 0.281 → 0.837 |
| the Himalaya and Tibet | 4.54 → 2.19 | 0.142 → 0.049 | 0 → 0.197 | 0.024 → 0.141 | 0.013 → 0.071 | 0.034 → 0.197 | 0.174 → 0.439 |
| land mean | 4.07 → 3.50 | 0.185 → 0.144 | 0 → 0.058 | 0.009 → 0.044 | 0.010 → 0.069 | 0.015 → 0.093 | 0.200 → 0.292 |

In the December land mean the IFS's constants without the form drag give
blocking 0.086 and waves 0.117 N/m² (the fields alone 0.049 and 0.043); the
form drag, slowing the low-level wind, takes 0.058 and leaves blocking
0.073 and waves 0.099. Over the mountains the vegetation's stress falls to
9–54 % of its value as the lowest layer's wind halves or more; the total
over land rises by 41 % (December), 34 % (June) and 46 % (N=128).

Land by cover and the global fluxes (8-step samples, bd42fa5 → these
commits; wind m/s at the lowest layer, stress the turbulent surface
stress ρ C_D U |v| in N/m², H and LE upward W/m², skin °C):

| class | eight64 day 183 + 3: wind; stress; H; LE; skin | eight128 day 183 + 1 |
|---|---|---|
| forest (trees ≥ 0.5) | 2.70 → 2.18; 0.250 → 0.185; 8.0 → 10.8; 58.5 → 56.0; 14.46 → 14.59 | 2.76 → 2.29; 0.243 → 0.187; 9.1 → 12.2; 67.4 → 64.3; 13.99 → 14.05 |
| grass (≥ 0.5) | 5.27 → 2.71; 0.173 → 0.062; 22.0 → 19.0; 24.4 → 21.8; 5.02 → 5.47 | 5.31 → 2.67; 0.171 → 0.059; 25.9 → 21.8; 31.2 → 27.8; 4.48 → 4.84 |
| bare (≥ 0.5) | 3.69 → 3.35; 0.137 → 0.117; 51.0 → 50.6; 23.6 → 23.1; 28.49 → 28.64 | 3.80 → 3.55; 0.142 → 0.126; 61.3 → 60.6; 25.6 → 25.2; 29.62 → 29.68 |
| mixed land | 3.45 → 2.91; 0.233 → 0.177; 23.9 → 25.8; 41.5 → 39.8; 16.78 → 16.93 | 3.47 → 3.03; 0.233 → 0.190; 33.4 → 34.8; 48.2 → 46.5; 17.57 → 17.61 |
| snow-covered land | 8.48 → 3.94; 0.250 → 0.073; −17.7 → −6.4; 6.9 → 2.7; −3.24 → −3.66 | 7.95 → 3.37; 0.234 → 0.061; −15.2 → −4.7; 8.5 → 3.4; −4.53 → −4.62 |
| ice sheets | 10.43 → 8.81; 0.241 → 0.173; −23.2 → −17.6; 2.2 → 1.5; −34.83 → −35.29 | 9.58 → 8.58; 0.203 → 0.162; −17.8 → −15.2; 1.8 → 1.4; −35.30 → −35.43 |
| global evaporation, mm/d; sensible heat, W/m² | 1.949 → 1.930; 9.92 → 10.36 | 2.042 → 2.021; 11.35 → 11.67 |
| sea stress; equatorial Pacific τx; Southern Ocean peak, N/m² | 0.089 → 0.087; −0.018 → −0.018; 0.191 → 0.186 at 47.5S | 0.086 → 0.085; −0.022 → −0.022; 0.211 → 0.211 at 57.5S |
| global surface; land, last day, °C | 16.050 → 16.082; 14.02 → 14.15 | 16.261 → 16.280; 14.95 → 15.02 |

The observed references for the stress over mountains and the
near-surface wind over hilly land are Beljaars et al.'s (2004); the paper
was not read here (the approval covered the USGS data only), so no
observed number is set against these. The circulation references are the
diagnosis's above: the polar low, 9–13 hPa too deep at 60–75N in
December, is now 2.5–6 hPa too deep against the NCEP–NCAR zonal mean near
1010 hPa; the land westerlies at 45–55N, 2.0–3.0 m/s before against ERA5's
3–5, are now 0.6–1.3 m/s.

Step cost under the exclusive lock (`profileGpu`, bd42fa5 and 35c9fae on
copies of eight64 day 183 and eight128 day 183, 128 steps after 16,
twice): N=64 20.81–20.89 → 21.21–21.42 ms a step (5.34 → 5.46 s a model
day), the adjustment and momentum kernels 4.34–4.36 → 4.57–4.67 ms and the
physics, boundary-layer and orography kernels 3.34–3.35 → 3.50–3.57 ms;
N=128 91.02 → 92.61–92.69 ms (46.6 → 47.4 s a model day), those kernels
15.90–15.95 → 16.63–16.71 and 13.36–13.41 → 13.94–14.03 ms. The form drag
adds (K − k_top + 1) C + E floats to the GPU's physics buffer (16 C + E,
12.4 MB at N=128). The four files add 2,176,144 bytes to `data/`; the page
fetches one of them (1.64 MB at N=128, 410 kB at N=64).

The review (Oct 2), on aa98bbf and the fixes after it. The cache holds
the 108 tiles (each the size the fetch log and the server's Content-Length
give, all from edcintl.cr.usgs.gov), the fetch log and listing probes, and
the derived grids (`gmted_mea300.npy`, four `fine_*.f32`), 4.1 GB; no file
in it is executable and none is in the repository. Spot heights read from
the tiles by their own tie points: Everest 8625 m, K2 8190, Aconcagua 6678,
Kilimanjaro 5778 (the next pixel), Denali 6035 (the next pixel), Mont Blanc
4505 (the next pixel), each the 7 × 7 maximum within a pixel of the
summit's coordinates; lake surfaces Titicaca 3815 m (3812), Baikal 449
(456), the Dead Sea −412, the Caspian 0 (−28: GMTED holds it at 0); Dome A
4086 (4093), Greenland's Summit 3206 (3216); 30N 40W and 0N 160W 0 m.
Six cells (N=64 and N=128 at the Alps 46.5N 9.5E, the Ethiopian highlands
9N 38.5E across the 10N tile edge, the Alaska Range 63N 149W) recomputed
from the tiles by a script of the review's own (each smoothing by direct
sums with great-circle distances, σ_flt over each 30″ point): μ to 0.03 %,
σ_flt to 0.8 % (the 2′30″ blocks' version to 0.05 %), and the slope
from neighbouring 2′30″ points to 0.1 % of the files of aa98bbf; against
exact differences 5 km each way that σ was 26–30 % too large at 46.5N
and 63N (3 % at 9N) and the Alps' and the Alaska Range's θ −0.16 and −0.27 rad,
steepest east–west across ridges that run east–west (−1.08 and −1.15 with
5 km each way). a8105bd takes the IFS's 5 km spacing; with exact 5 km
offsets on the 30″ field σ is 2–4 % above the files' bilinear version, γ
within 0.05, θ within 0.13 rad. The band filter's integral reproduces I_H
0.00102 m⁻¹; eq. 11.13 as printed gives k_flt 0.00051 m⁻¹ for that filter
against the documented 0.00035 (taken, with the α 35 tuned with it). The
form drag by hand (σ_flt 150 m, 12 m/s, N=16, no blocking): each rate to
3·10⁻¹⁴ on the CPU, also inside a full step on the state it diagnosed;
ρ∫ C_tofd U² dz from 20 m 1.815 N/m² (0.565 × 2.25 × 1.44 = 1.831); the
layer sum Σ m C_tofd U² with the lowest midpoint at 18 m 3.0 N/m² before
the implicit reduction; each edge u/(1 + Δt C_tofd |U|) exactly. The
page's path (fetch, `arrayBuffer`, `decodeSubgrid`) equals the node loader
to float32 for N=16 and 64. The terrain-off and other-land-mask fallbacks
(b3ba2d8) and the worker's decode fallback (35c9fae) are the review's. The
December ten-day run and its twin on aa98bbf reproduced the reports above
bit for bit before a8105bd; every number in the run tables is from
a8105bd's fields. Full suite 55 files, six at a time under the shared lock:
all pass but parallel.test (its speedup, 1.0× under the load of other
jobs), which passes alone (3.7×).

What still misses: the IFS's G for its operational fields is not stated
in the documentation read, so G 1.23 is the textbook mountain's and the
doubled H_eff rests on eq. 4.8 against eq. 4.37's remark; with them the
December land westerlies at 45–55N fall to 0.6–1.3 m/s, below ERA5's 3–5,
while the polar SLP rises toward the reanalysis: the two references now
pull opposite ways and ten days cannot say which constant is wrong; the
5 km field is sampled at 2′30″ and its differences are bilinear between
those points (σ 2–4 % below exact 5 km offsets on six cells); σ_flt
averages 2′30″ blocks by their centres' cells, not each 30″ point (0.2–0.8 %
on six cells); the form drag's |U| is the cell's at the diagnosis, averaged
onto the edge; the wave stress that the per-step limit carries to the top
layer reaches 6.3 m/s a step there (the top layer's own wind is the
bound); the form drag's 2.109 and z^(−1.2) fit is the IFS's for n₁ −1.9
where these data give −1.85; the summer and winter responses need a
season, which needs a spin-up.

**The terrain's fields merged (Oct 2).** Branch integrate-c: roughness
fe5de7f onto 9ba8c38. The form drag is in the boundary layer's edge
solve with the surface drag, its stress the ground's (`formStress`,
`PH_FSTRESS`), not the ocean's, its energy in the solve's dissipation
once. Tests: the sunlit cloudy layers' heating may part by 1.5·10⁻⁴
K/day (1.26·10⁻⁴ at one cell; the CPU's own heating there moves by
7.7·10⁻⁵ when one input changes by one float32 ulp, 1.3–1.5·10⁻⁴ when
every input of the column changes by ±1 ulp); the cloud-effect parity
leaves out the columns whose dry adjustment, shallow plume top or
cloudy layers differ between the engines after any step, at most 2 %
(4 of 362: 1 merged, 3 plume tops, 3 cloudy-layer patterns; the rest
at most 1.9·10⁻² W/m² apart in the shortwave, 1.2·10⁻² in the
longwave). No digest moved; 63 files pass.
Proofs, nine64_day0274, the second CPU step after loading: each edge
column's momentum + Δt (τ_s + τ_form) 2.3·10⁻¹¹ of 4.8·10⁵ kg/m/s (Δt
τ_form up to 3.3·10³), blocking alone 2.3·10⁻¹³, the orographic waves
alone 4.6·10⁻¹³, each cell column's wave force 6.5·10⁻¹⁹ of 4.2·10⁻³
Pa; the sponge 3.8·10¹⁵ and 1.1·10¹⁵ J, its angular momentum 0.033 and
0.020 of a Rayleigh drag's; heat applied 1.36·10¹⁸ J against closure
1.42·10¹⁷ + surface and form drag 8.52·10¹⁷ (the form drag's 7.5·10¹⁶)
+ mountains 3.64·10¹⁷ + waves 5.2·10¹⁴ to 5·10⁻¹⁴; 37998 edges take
form stress, none sea–sea; no sea cell holds σ_flt, a form or blocking
rate or a wave tendency; the ocean's stress is the surface stress (0
apart). GPU kernels: drag 1.8·10⁻⁵, mountains 1.6·10⁻⁵, waves 4.2·10⁻⁸
relative, heat against the sinks 7.8·10⁻⁴. With `subgrid` false the CPU
gives 9ba8c38's state, land and ocean digests after each of nine steps
from eight64_day0183; N=8 (no file) and N=16 and N=64 with the terrain
off take the raster with its note and give 9ba8c38's six-step digests,
the GPU the raster, Lott and Miller's constants and no form layers.
Parity after 1, 4, 16 steps, lowest layer u rms (9ba8c38; fe5de7f): from
nine64_day0091 4.8·10⁻⁵, 1.4·10⁻³, 2.7·10⁻³ m/s (4.9·10⁻⁵, 2.5·10⁻⁴,
1.5·10⁻²; 1.9·10⁻³, 2.0·10⁻³, 3.7·10⁻³), launched stress 2.4·10⁻⁴,
1.3·10⁻³, 3.8·10⁻² of its largest (1.0·10⁻⁴, 9.6·10⁻⁴, 1.5·10⁻³;
4.7·10⁻⁴, 9.9·10⁻², 2.2·10⁻²), form stress 6.6·10⁻⁶, 8.1·10⁻⁴,
4.8·10⁻⁴ rms relative, regime flips 0, 3, 32 (0, 2, 27; 9, 11, 73);
from nine64_day0274 3.8·10⁻⁴, 1.9·10⁻³, 8.1·10⁻³ (4.9·10⁻⁵, 7.0·10⁻³,
5.5·10⁻³; 7.9·10⁻⁴, 1.9·10⁻³, 5.0·10⁻³). Three one-day segments end day
186 byte for byte as one. Ten CPU steps: the form rate × Δt 9.02 (30.0N
94.5E, lowest layer) from nine64_day0274, 6.58 from nine64_day0091,
6.30 (29.7S 70.2W) at N=128; no column's largest |u| grows and no
column's kinetic energy rises in the solve; balance 3.5·10⁻¹¹ to
1.3·10⁻⁹; all finite, largest wind 95.8, 85.8, 84.0 m/s.
five64_day2190 (cam26) loads and steps on both engines.

GPU runs, 8-step samples, one script on the three trees (it gives
fe5de7f's December numbers above to the digits printed), merge
(9ba8c38; fe5de7f): eight64 day 186 albedo 0.291 (0.293; 0.312), ASR
241.4 (240.9; 234.4), OLR 234.1 (233.8; 241.7), SWCRE −52.4 (−52.8;
−54.3), LWCRE 25.9 (26.1; 17.3), rain 1.70 (1.72; 1.67); nine64 day 94
0.303 (0.305; 0.329), 237.3 (236.8; 228.4), 234.3 (234.3; 239.6),
−57.1 (−57.6; −58.7), 26.1 (26.0; 20.0), 2.36 (2.36; 2.26). Land by
cover, eight64 day 183 + 3, wind m/s; stress N/m²; H; LE W/m²: forest
2.14; 0.197; 14.1; 64.0 (2.65; 0.258; 11.8; 66.4 / 2.18; 0.185; 10.8;
56.0), grass 2.71; 0.062; 22.0; 24.7 (5.19; 0.171; 25.5; 27.3 / 2.71;
0.062; 19.0; 21.8), bare 3.37; 0.123; 61.0; 25.3 (3.68; 0.142; 61.6;
25.7 / 3.35; 0.117; 50.6; 23.1), snow 4.00; 0.074; −6.2; 2.7 (8.43;
0.245; −16.8; 6.7 / 3.94; 0.073; −6.4; 2.7), ice sheets 8.67; 0.163;
−20.9; 2.0 (10.35; 0.233; −26.7; 2.8 / 8.81; 0.173; −17.6; 1.5);
evaporation 2.196 (2.217; 1.930) mm/d, sensible 12.50 (12.11; 10.36)
W/m², calm tropical sea LE 45.5 (45.1; 28.3); nine64 day 91 + 3 forest
wind 2.10 (2.67; 2.11), grass 2.21 (4.35; 2.12), evaporation 2.711
(2.726; 2.505). nine64 day 274 + 10: SLP 60–65 / 70–75 / 80–85N
1006.6 / 1003.3 / 1008.4 (1002.4 / 997.8 / 1002.4; 1007.5 / 1004.0 /
1010.1) hPa, 30–85N rms against 9ba8c38 4.26 hPa against a 10⁻⁷ twin's
0.12; lowest-layer u over land at 45–50 / 50–55N 1.1 / 0.8 (2.7 / 1.9;
1.3 / 0.6) m/s; 200 hPa at 50–55N 32.5 (32.6; 32.7), rms 0.76 against
9ba8c38 (twin 0.08); land-mean stress vegetation / form / blocking /
waves / total 0.112 / 0.052 / 0.068 / 0.085 / 0.296 (0.188 / 0 / 0.020
/ 0.017 / 0.220; 0.125 / 0.058 / 0.073 / 0.099 / 0.328) N/m²; day 284
at 70–90N, 1.1 / 3.5 / 7.4 / 14 / 24 / 37 / 53 hPa, 225.6 / 229.8 /
223.1 / 216.3 / 210.4 / 206.7 / 204.2 K (220.7 / 225.2 / 218.9 / 213.1
/ 208.1 / 205.0 / 203.0; 240.8 / 249.6 / 239.6 / 229.9 / 218.4 / 209.0
/ 204.3; twin within 0.1 K), zonal-mean u at 65–67.5N 38.6 / 34.2 /
36.0 / 35.6 / 34.6 / 32.5 / 29.4 (48.1 / 39.5 / 40.1 / 38.6 / 36.1 /
32.9 / 29.8) m/s; without the form drag 225.2 / 229.8 / 223.1 / 216.5 K
and 39.0 / 34.2 / 35.9 / 36.5 m/s. eight128 day 183 + 1: albedo 0.252
(0.253; 0.281), land-mean wind 3.48 (4.04; 3.50) m/s, stress 0.149 /
0.057 / 0.042 / 0.064 (0.189 / 0 / 0.008 / 0.010; 0.144 / 0.058 / 0.044
/ 0.069). Cost under the exclusive lock, 128 steps after 16, twice:
N=64 25.18 / 25.04 → 25.65 / 26.14 ms (+3.4 %); N=128 105.90 / 106.76
→ 108.13 / 106.60 ms (+0.9 %; physics pass 42.24 / 42.34 → 43.94 /
43.48, mixing pass 17.65 / 17.74 → 18.65 / 18.22); a day at N=128 54.4
→ 55.0 s of steps, 77.1 → 78.1 s to the saved file.

*Review.* The four conflicted files rebuilt with git merge-tree: in each
of the 13 code, test and doc files both sides touched, the lines the
merge adds and removes against 9ba8c38 are fe5de7f's against bd42fa5
line for line, apart from the signatures of layoutFor, createGpuModel
and createModel and the layoutFor call, which carry both sides'
parameters; PH keeps 9ba8c38's 97 slots in order and adds OFLT, TOFD
and FSTRESS; the GPU dispatches pblDiagnose, gravityWaves, orography,
then mixMomentum, orographyApply, gravityWaveDrag, dissipationHeat, the
CPU's order. 63 files, 537 tests pass. With `subgrid` false the GPU too
gives 9ba8c38's state digests after each of nine steps from
eight64_day0183, and the note now says the files are turned off. Own
harnesses, nine64_day0274: the CPU's second step closure 1.421·10¹⁷,
surface and form drag 8.521·10¹⁷, mountains 3.645·10¹⁷, waves
5.244·10¹⁴ J, heat 1.359·10¹⁸ J to 2.4·10⁻¹⁴; the GPU's kernels, run
again after one step, edge momentum against Δt (τ_s + τ_form) to
4.3·10⁻⁶ of the largest, mountains 1.5·10⁻⁷, heat against the
dissipation 3.5·10⁻⁴, 11882 land cells with σ_flt and none at sea. The
CPU's heating at cell 196 layer 24 under ±1 ulp on every input of the
column, another generator: 1.47·10⁻⁴ (20 draws) to 1.64·10⁻⁴ K/day (200
draws), so 1.5·10⁻⁴ lies inside that response, not above it. The
cloud-effect columns 48 and 63 (0.66 W/m²) are first a shallow plume
top apart: at step 9 it is 399 hPa on the CPU and 339 hPa on the GPU
(base flux 4.48·10⁻² against 4.06·10⁻²; both 399 at step 7, 339 at step
8), which cools layer 18 by 0.14 K on one engine and layer 17 by 0.16 K
on the other; their cloudy layers part at step 11. Column 4 (1.59
W/m², merged) is a plume top apart at step 13. The test now compares
the plume top after every step as well; it leaves out the same 4
columns. From nine64_day0274 after one step the launched stress parts
by 6.4·10⁻² of its largest at 49.9N 122.1W, where the CPU blocks to
589.5 m and the GPU not at all; ±10⁻⁶ on the CPU column's winds and
±10⁻⁷ on its θ switch the blocking in 1 of 40 draws (9ba8c38 2.7·10⁻⁴,
fe5de7f 9.2·10⁻⁵; after 4 and 16 steps 2.7·10⁻³ and 5.0·10⁻³, 9ba8c38
2.7·10⁻³ and 8.4·10⁻³, fe5de7f 7.5·10⁻³ and 5.9·10⁻²). The three days
from eight64_day0183 repeat day 186's line and the land-by-cover report
above, one segment byte for byte as three. Ten steps from
eight128_day0183: finite on both engines, largest wind 84.0 then
77.5–82.3 m/s, the form rate × Δt 6.30 at 29.7S 70.2W on the first step
and 1.33, then 0.74 at the tenth on the GPU; the CPU's solve never
raises an edge column's kinetic energy or its largest |u|.

**The gravity waves' sources and bl36's top merged (Oct 2).** Branch
integrate-c: gas-benchmark 4a4eb83 onto 0063c54 (86d02b9). The new wave
source and lid sit beside the mountains' drag and the form drag, the
level set's longwave table beside the cloud optics and the clamped dry
fraction, the remap beside the land's regridding; fresh runs start on
bl36. The twelve-step digests all moved with the fresh stratosphere read
by pressure (the Rayleigh-top ones too), re-pinned; on grey ice the
three deck digests equal 4a4eb83's. Tests: 61 of 63 files pass, the five
that run bl36 among them. dayMeans (reflected sunlight per-cell rms
1.26·10⁻⁴ against 10⁻⁴; one column, 2.7 W/m² of summed reflection, whose
plume base flux parts by 5 % at step 7 with no decision the cloud-effect
test watches apart) and gpuModel's uniform-condensation heating
(2.68·10⁻⁴ K/day at cell 109 layer 24 against 1e-5 × 26.6; the CPU's
own response to ±1 ulp on every input of that column 2.40·10⁻⁴ in 200
draws) fail; both pass with 0063c54's init.module.js on the merged tree,
as did cloudEffect, whose one column at 0.544 W/m² neighbours column 299
(dry adjustment, plume top and cloudy layers apart from step 9): its
neighbours are now left out too (abcbea6; 14 of 362, the rest within
0.048 W/m²). Proofs: with `speedStep` 4 and `sourceDescent` false the
CPU gives 0063c54's state from eight64_day0183 to 7.7·10⁻¹⁷ of the
largest u after one step (the top layer, the spectrum's flux taken out
as a factor) and to 6.5·10⁻¹² (qc), 2.4·10⁻¹² (q), 5.6·10⁻¹³ (u),
3.9·10⁻¹⁴ (θ) after nine; the defaults move u by 3.0·10⁻⁴ in one step.
bl34 states load and step on both engines (eight64_day0183 three days:
day 186 albedo 0.291, ASR 241.3, OLR 234.0, SWCRE −52.4, LWCRE 25.9, rain
1.70, 0063c54's within 0.1 W/m²). On nine64_day0274 remapped to bl36,
the second CPU step: each edge column's momentum + Δt (τ_s + τ_form)
2.6·10⁻¹¹ of 4.8·10⁵ kg/m/s, blocking 2.3·10⁻¹³, orographic waves
4.6·10⁻¹³, each cell column's wave force 1.1·10⁻¹⁸ of 4.3·10⁻³ Pa; the
sponge 5.7·10¹⁴ and 2.6·10¹⁴ J at 0.15 and 0.65 hPa, its angular
momentum 0.033 of a Rayleigh drag's; heat 1.15·10¹⁸ J against the sinks
to 7.9·10⁻¹⁴; GPU kernels drag 1.8·10⁻⁵, mountains 1.6·10⁻⁵, waves
4.7·10⁻⁸, heat 1.3·10⁻³. Parity on bl36 after 1, 4, 16 steps, from
nine64_day0091 (in parentheses nine64_day0274): top eight layers u 8.5·10⁻⁵, 2.1·10⁻⁴,
8.4·10⁻⁴ m/s rms (9.1·10⁻⁵, 1.1·10⁻³, 1.1·10⁻³; largest 0.46 m/s at step
4 in the top layer at 47.6N 142.5E, at step 16 in layer 6 at 49.2N
119.1W), θ 1.8·10⁻⁴, 6.4·10⁻⁴, 2.6·10⁻³ K; lowest layer T 1.0·10⁻⁴,
4.6·10⁻⁴, 2.2·10⁻³ K, u 4.8·10⁻⁵, 4.3·10⁻⁴, 6.4·10⁻³ m/s (4.0·10⁻⁴,
1.9·10⁻⁴, 2.1·10⁻³; the launched stress 6.0·10⁻² of its largest at step
1, the blocking column); wave acceleration 2.6·10⁻⁴, 9.7·10⁻⁴,
1.8·10⁻³ m/s/day rms; regime flips 0, 0, 21. The remap
(eight64_day0183, nine64_day0274): each column's Σ θΠ dσ to 2.4·10⁻¹⁶,
q and qc exactly, u to 3.4·10⁻¹⁶; the θ integral 6.4–6.6·10⁻⁴; 0.15 /
0.65 / 1.6 hPa 242–243 / 257 / 266 K under 272 K at 3.6 hPa, θ rising
upward in every column; the round trip to 5·10⁻¹⁶. Three one-day
segments from the remapped eight64_day0183 end day 186 byte for byte as
one. radiationBenchmark on bl34 is 0063c54's line for line, on bl36
4a4eb83's line for line (M21's bl36 rows above as the script prints
them: TROP −5.05 / −11.08 / −8.83, MLW −8.56 / −12.01 / −7.40 K/day). A paired spin-up on bl36: `NS="64 128" PREFIX=<new>
LEVELS=bl36 OCEAN='{"everySteps":8}' scripts/pairedSpinup.sh` with no
`<PREFIX><N>_day*.bin` in OUT (`STRATOSPHERE=1` for the upper lines);
the land starts neutral and jumps at days 365 and 730 by default. Two
fresh days so: N=64 largest wind 72 / 93 m/s, Courant 0.23 / 0.30, day
2 ASR 204.5, OLR 200.1, albedo 0.399; N=128 77 / 97 m/s, 0.25 / 0.31,
215.8, 202.2, 0.366; both read their terrain files (no fallback note),
so the mountains' drag and the form drag act from the first step.

The top's pre-flight: 60 days at N=64 from a fresh atlas start on bl36,
the defaults (`runs/igpre64.log` in the worktree), day 0 the March
equinox. Top six layers (0.15 / 0.64 / 1.6 / 3.5 / 7.4 / 14 hPa); 5S-5N,
57.5–62.5S and N zonal-mean u in m/s; the largest wind and horizontal
Courant number there; 70-90N and 70-90S layer-mean T in K:

| day | 5S-5N | 60S | 60N | largest, Courant | 70-90N | 70-90S |
|---|---|---|---|---|---|---|
| 10 | +1 / +4 / 0 / 0 / 0 / 0 | 22 / 25 / 22 / 17 / 14 / 12 | 22 / 22 / 19 / 16 / 13 / 12 | 53, 0.17 | 232 / 251 / 248 / 238 / 230 / 223 | 223 / 241 / 241 / 236 / 230 / 223 |
| 20 | +12 / +12 / +12 / −1 / −2 / −3 | 47 / 39 / 29 / 22 / 18 / 16 | 10 / 11 / 8 / 3 / 0 / −1 | 75, 0.24 | 235 / 254 / 249 / 234 / 224 / 220 | 218 / 236 / 234 / 227 / 219 / 215 |
| 30 | 0 / +10 / +34 / 0 / −3 / −3 | 67 / 49 / 34 / 25 / 20 / 17 | 8 / 10 / 7 / 2 / −1 / −3 | 99, 0.32 | 237 / 258 / 255 / 239 / 227 / 221 | 213 / 231 / 230 / 222 / 215 / 211 |
| 40 | −38 / +2 / +49 / +13 / −6 / −8 | 70 / 56 / 42 / 33 / 28 / 25 | −2 / 2 / 0 / −4 / −6 / −4 | 100, 0.32 | 240 / 263 / 260 / 244 / 230 / 224 | 212 / 228 / 224 / 218 / 209 / 205 |
| 50 | −86 / −16 / +47 / +25 / −4 / −11 | 95 / 72 / 56 / 47 / 41 / 37 | −18 / −9 / −7 / −7 / −6 / −5 | 121, 0.39 | 243 / 267 / 266 / 248 / 234 / 226 | 213 / 226 / 223 / 215 / 206 / 201 |
| 60 | −100 / −27 / +42 / +32 / +6 / −13 | 146 / 98 / 70 / 56 / 47 / 40 | −35 / −18 / −12 / −10 / −8 / −6 | 156, 0.50 | 245 / 270 / 270 / 252 / 236 / 228 | 210 / 223 / 218 / 211 / 202 / 197 |

Largest over the run 159 m/s and horizontal Courant 0.51 (day 59),
vertical 0.09 (day 1). Global-mean T
on day 60 of the layers above 100 hPa, 0.15 … 85 hPa: 227.7 / 252.8 /
256.6 / 242.5 / 227.4 / 219.0 / 214.1 / 210.2 / 207.6 / 205.8 / 204.7 K;
trends over days 31-60 −0.018 / +0.009 / −0.013 / −0.013 / −0.020 /
−0.014 / +0.001 / +0.019 / +0.033 / +0.044 / +0.060 K/day. The caps over
days 31-60: 70-90N +0.27 to +0.50 K/day, 70-90S −0.03 / −0.21 / −0.29 /
−0.33 / −0.37 / −0.41 K/day (the season's). The 5S-5N wind's trend in
the top four layers, days 1-20 against 41-60: +1.01 / +0.96 / +0.54 /
−0.03 against −2.80 / −1.46 / −0.38 / +0.90 m/s/day; over days 41-50
−4.84 / −1.61 / −0.43 / +1.11 and 51-60 −0.95 / −0.81 / −0.60 / +0.55.
The 60S wind at 0.15 hPa gains 3.8 m/s/day over days 41-60, the
day's largest wind 2.3 m/s/day over 31-60. Verdict: above 2 hPa the
equatorial winds drift, one way, for 60 days, the 0.15 hPa layer
easterly to −100 m/s, slowing from −4.8 to −0.95 m/s/day but not level,
the 0.64 hPa layer still at −0.8 m/s/day, the 1.6 hPa westerly level at
+42 to +51 since day 35; no reversal, so no oscillation is seen. Nothing
breaks in 60 days; for a multi-year run the open risks are the top
layer's −100 m/s easterly still growing, the winter jet and Courant
number still rising toward June (0.50 on day 60, 0.62 seen at July in
M21), and the winter cap at 210 / 223 / 218 K over 0.15-1.6 hPa, 30-40 K
under AFGL's subarctic winter.

Troposphere and cost. Three N=64 days from eight64_day0183, bl34 /
remapped bl36, days 184-186: ASR 253.3 / 253.3, 249.1 / 249.1, 241.3 /
241.3; OLR 228.9 / 228.6, 235.8 / 235.5, 234.0 / 233.7 W/m²; albedo
0.256 / 0.256, 0.269 / 0.268, 0.291 / 0.291; LWCRE 34.6 / 34.4, 25.7 /
25.6, 25.9 / 25.8; rain 0.81, 1.05, 1.70 mm/d on both. One N=128 day
from eight128_day0183: ASR 254.7 / 254.7, OLR 231.7 / 231.5, albedo
0.252, SWCRE −38.8, LWCRE 32.1 / 32.0, rain 0.77, largest wind 78.0 /
84.8 m/s. Step cost under the exclusive lock, 128 steps after 16,
twice: N=64 bl34 25.69 / 26.70, bl36 27.69 / 29.07 ms (+8 %); N=128 bl34
108.75 / 108.59, bl36 117.05 / 116.28 ms (+7.4 %). A day at N=128 on
bl36: 59.7 s of steps, 77.7 s from the end of setup (20 s) to the saved
file (bl34 55.6 and 74.2 s).

The merge's review (Oct 2). The five conflicted files rebuilt with `git
merge-tree` from 0063c54 and 4a4eb83: the committed merge differs from
git's own only in the conflict hunks, each the union of both sides (the
imports, `longwaveTable` beside the cloud optics, 0063c54's clamped dry
fraction passing 4a4eb83's `gasTable`, both sets of regrid imports,
0063c54's digest tests re-pinned, the three on grey ice equal to
4a4eb83's); `GWS` and `GWF` close the mesh buffers' lists and the kernel
order is 0063c54's. The suite, six files at a time under the shared lock:
61 of 63 pass, dayMeans and gpuModel fail as above. dayMeans' column 135
takes a plume base flux of 4.287·10⁻² on the CPU and 4.504·10⁻² on the
GPU at step 7 (0063c54: 4.508 and 4.508·10⁻²); the CPU alone, with π, θ,
u, the surface temperature and q each scaled by 1 ± 5·10⁻⁷ at random,
gives 4.504·10⁻² in 3 of 7 draws and 4.287·10⁻² in the others, so the
column sits on a threshold that float32 noise crosses. At 4a4eb83 itself
dayMeans passes (reflected rms 3.2·10⁻⁵) and gpuModel's sunlit-heating
test fails at 1.34·10⁻⁴ K/day against that branch's 10⁻⁴. Parity on bl36
with a loader of the review's own, after 1 / 4 / 16 steps: from
nine64_day0091 the top eight layers' u 8.5·10⁻⁵ / 2.1·10⁻⁴ / 9.0·10⁻⁴
m/s rms and θ 1.8·10⁻⁴ / 6.4·10⁻⁴ / 3.0·10⁻³ K, from nine64_day0274
9.1·10⁻⁵ / 1.1·10⁻³ / 1.3·10⁻³ and 2.6·10⁻⁴ / 1.5·10⁻³ / 3.5·10⁻³; the
lowest layer's T 1.0·10⁻⁴ / 1.7·10⁻³ / 1.9·10⁻² K, as 0063c54 gives on
bl34 with the same loader (1.0·10⁻⁴ / 1.9·10⁻³ / 1.9·10⁻²); the wave
accelerations 4.1·10⁻⁵ / 7.0·10⁻⁴ / 1.8·10⁻³ m/s/day rms apart, each
cell column's Σ dσ·a 1.3-1.6·10⁻¹⁶ (CPU) and 5.1-6.2·10⁻⁸ (GPU) of its
deposit, the two lid layers' accelerations equal to 9.4·10⁻¹⁵ m/s/day.
The orographic waves deposit 0.99972 of the launched stress on
nine64_day0274 on bl34 and bl36 alike, 5.9·10⁻³ of it in bl36's two top
layers (bl34's two 7.3·10⁻³). The remap of eight64_day0183 and
nine64_day0274: column enthalpy 2.5·10⁻¹⁶, q and qc exact, u 3.4·10⁻¹⁶,
no θ inversion in any column's top seven layers, the 34 shared layers
bit for bit, the round trip 5.1·10⁻¹⁶. The three bl34 days from
eight64_day0183 give the numbers above to the last digit; on bl36 three
one-day segments equal one three-day segment byte for byte.

Thirty July days on the merged tree from nine64_day0091 remapped
(`STRATOSPHERE=1`), against the second round's review on 4a4eb83 (its
twin in brackets): day 121, the 5S-5N wind at 0.15 / 0.64 / 1.6 hPa −2 /
+15 / −35 m/s (+8 / +17 / −31 [+5 / +40 / −29]), the 70-90S cap 211.5 /
226 / 222 K (216 / 232 / 225 [210 / 226 / 221]), the top six layers'
largest wind 140 m/s and horizontal Courant number 0.45 (152, 0.49
[185, 0.59]), vertical 0.09; the 0.64 hPa wind gains 2.8 m/s/day over
days 112-121. The merge leaves bl36's top as 4a4eb83 had it.

The pre-flight's log read again, least-squares trends over days 1-20,
21-40 and 41-60 (m/s/day): the 5S-5N wind at 0.15 hPa +1.01, −2.56,
−2.80; 0.64 hPa +0.96, −0.54, −1.46; 1.6 hPa +0.54, +1.77, −0.38; 3.5
hPa −0.03, +0.71, +0.90. Over the last ten days the 0.15 hPa wind rises on
one day only (−86 on day 50, −100 on day 60) and the 0.64 hPa wind on
none (−16, −27). The
global means above 30 hPa are flat because the tropics and the winter
cap cool while the summer cap warms: over days 41-60, 20S-20N −0.35 /
−0.27 / −0.17 K/day at 0.15 / 0.64 / 1.6 hPa (229.0 → 220.1, 255.0 →
247.3, 261.8 → 257.0 K from day 20 to day 60; −0.42 / −0.37 / −0.24
over days 51-60), 70-90S −0.22 / −0.35 / −0.37, 70-90N +0.24 / +0.33 /
+0.46. The 60S wind at 0.15 hPa gains 3.8 m/s/day over days 41-60 and
5.2 over 51-60 (95 → 146 m/s); the top six layers' largest wind 2.25
m/s/day over days 31-60 and 4.3 over 51-60; the top layer's horizontal
Courant number 0.39 on day 50, 0.50 on day 60. The July runs from the
spun-up state hold the top layer's 5S-5N wind between −25 and −2 m/s;
the fresh start takes it to −100 by day 60. Verdict: bl36's top is not
shown safe for a multi-year run from a fresh start. Nothing breaks in 60
days, but on day 60 the top layer's easterly, the tropical top's cooling,
the winter jet and the Courant number are all still moving, the last two
faster over the last ten days than over the twenty before, with the
solstice 31 days off.

**The convection and the model top merged (Oct 2).** de82fce merges
0063c54 (the eddy sponge and gravity-wave drag, the mountains' and form
drag on GMTED fields, the gust) into elements 0–4. Every parent digest at
N=4 reproduces with its parent's options (the convection side's under the
Rayleigh top, the model top's with the longwave's random overlap);
re-pinned where both sides moved the inputs: the moist defaults 1dbd465a,
the model top under elements 2–4 01214faf. Three of the suite's files
failed on the merge and were settled: the audit's bulk sensible heat took
max(wind, 3 m/s) where the model now takes the gust wind (82.67 against
94.21 W/m², exact after 388f6b6); the bl34 continent case parted at one
land cell (168, at 947 m) whose absorbed sunlight at step 7 was 646.8
W/m² on the CPU and 546.6 on the GPU (reflected 636.8 against 736.9, OLR
184.8 against 159.1), with the same mixing top (6046 m above sea level),
the cloud water of its cloudy layers 19–21 within 0.9 % and their vapour
within 0.2 %, but their longwave heating apart (layer 21 −58.9 against
−69.2, layer 22 −10.2 against +9.2 K/d): the cover of those layers
decided apart in both bands; cell 189 parted so at step 6 (OLR 185.3
against 161.4). The case now leaves out such columns (2 of 362 at
3f3d915, none at dc14743);
the snow-albedo case's premise failed on the new closure's onsets (16 of
63 snow cells apart; 3 on the previous convection, on which it now runs).
The reference on this tree, three N=64 days from eight64_day0183 with two
Charnock replicates (×(1 ± 1e-4)), values with the replicate spread:
ITCZ convective share 0.871 (0.001), firing 0.448 (0.003), fired tops
above 300 hPa 0.024 (0.0001), Q1R centroid 725.8 hPa (0.4), T − Jordan at
516 / 439 / 848 hPa +0.90 / +0.96 / −1.91 K (0.0001 / 0.0006 / 0.009),
large-scale rain below 700 hPa 0.11 mm/d; warm pool share 0.744 (0.003),
centroid 682.7 hPa (0.25); N Pacific trades deep firing 0.758 (0.004),
convective rain 1.61 mm/d (0.006); SE Pacific convective share (8 audit
steps) 0.72, low cloud 0.459 (0.03); California radiative low cover 0.215
(0.007); trades low cover (cloudRegimes) 0.037 (0.005); the replayed
day's global rain 1.962 mm/d (0.003); GPU day 186 ASR − OLR 13.1, SWCRE
−43.3, LWCRE 22.8 W/m². From ten64_day0183: ITCZ share 0.747, fired tops
above 300 hPa 0.143, centroid 757.9 hPa, ITCZ rain 4.73 mm/d.

**What the deep plume did after elements 2–4 (Oct 2).** A day of CPU
steps from the reference's day-186 state (scratchpad ce/probe.mjs:
before each box column's moist step, twins replay the shallow plume
alone, the deep plume alone, the shallow plume after it, and the same
deep plume unentrained), over fired column-steps:
- Trades (fired 0.758): cloud base 934 hPa, cloud depth 200–300 hPa
  0.27, 300–400 hPa 0.45, 400–500 hPa 0.18; tops 500–600 hPa 0.45,
  600–700 hPa 0.27; PCAPE 24.0 Pa (PCAPE_bl 7.5), τ 92 min, base flux
  0.0097 kg/m²/s. The same columns' shallow plume alone tops at 800–900
  hPa (0.69). Together at 849 / 893 hPa the deep plume heats 3.78 / 4.63
  K/d and the shallow one −0.61 / −2.40, and the moistening (L/c_p
  dq/dt) is −1.04 / −0.61 K/d; the shallow plume alone cools −1.90 /
  −2.25 K/d and moistens +4.67 / +9.04 K/d: the deep plume heats and
  dries the trade-cumulus layer that the shallow plume alone moistens.
  The entraining plume's buoyancy is +0.4 to +0.8 K at 850–893 hPa and
  −0.04 to −0.8 K from 789 to 608 hPa, where it coasts to its top; the
  same plume unentrained keeps +0.6 to +1.2 K to 600 hPa. On the
  old-convection base state (cx0b) the trade plumes' clouds are as deep
  (300–400 hPa: 0.79 of the fired ones) but the 120 J/kg threshold fired
  0.19 of them; the 200 hPa criterion cannot separate them.
- SE Pacific (fired 0.082): depth 200–300 hPa 0.66, tops 600–700 hPa
  0.72, base flux 0.0185, PCAPE 23.7 Pa; California (fired 0.618): depth
  300–500 hPa 0.71, tops 400–600 hPa 0.85, PCAPE 33.0 Pa.
- ITCZ (fired 0.448): depth 200–300 / 300–400 / 400–500 hPa 0.17 / 0.25
  / 0.35, tops 400–500 hPa 0.51, 500–700 hPa 0.42, above 300 hPa 0.024;
  PCAPE 38.3 Pa, τ 88 min, base flux 0.0097. The 2.4 % topping above 300
  hPa have PCAPE 117 Pa, CAPE 137 J/kg; the rest 36 Pa and 37 J/kg. Their
  buoyancy, entraining: +0.26 / +0.46 / +0.75 / +0.95 / +0.95 / +0.65 /
  +0.28 / −0.01 / −0.22 K at 944 / 922 / 848 / 787 / 704 / 606 / 515 / 438
  / 372 hPa for the deep-topping ones; for the rest +0.38 / +0.57 / +0.29
  / +0.04 / −0.17 / −0.54 at 922 / 848 / 788 / 704 / 607 / 516 hPa,
  topping at 400–600 hPa on momentum; unentrained the same plumes keep
  +0.7 to +1.2 K from 850 to 370 hPa. Entrainment is 1.8–3.9 e-4 /m at
  790–920 hPa (Gregory's 0.1 B/w²) and the floor 1e-4 above; the source
  (lowest 50 hPa with the excess) is 1.0–2.5 K of h/c_p below the lowest
  layer, and h*/c_p is least (336 K) at 600–850 hPa.
- Why the tops fell: on the same old-convection base state the old
  closure fires 0.074 of ITCZ column-steps (CAPE 146 J/kg) with 0.44 of
  its tops above 300 hPa; the Bechtold closure fires 0.352 (CAPE 43
  J/kg) with 0.049 above 300 hPa, and the plumes topping above 300 hPa
  number 0.017 of the column-steps against 0.033. Most of the fall is
  selection (many weak plumes fired), the rest the state they leave
  (0.011 of column-steps after three days). The stable 600–850 hPa layer
  limits what entrainment of the cloud-base layers leaves: the plume's
  margin is 1 K.

**The type of convection by the cloud's depth (Oct 2).**
`convectionType` 'cloudDepth' (the default; 'top' the previous test on
the top above σ 0.7 with the shallow plume beside the deep one, bit for
bit on the CPU): the plume is deep if its cloud from base interface to top
interface is deeper than 200 hPa, and a column convects as one type: the
deep plume alone, or the shallow plume alone where the deep plume is not
deep or its closure gives no flux (IFS Cy43r1 §6.4 and §6.4.2; the IFS
text does not say what a deep column with no closure flux does). On 362
random columns 51 are deep by the depth only and none by the top only;
the type removes the shallow plume beside 48. Three days against the
reference (eight64 / ten64): trades deep firing 0.758 → 0.698,
convective rain 1.61 → 1.66 mm/d, low cover (cloudRegimes) 0.037 →
0.133; SE Pacific convective share 0.72 → 0.76 (CPU day 0.67 → 0.72), low
cloud 0.459 → 0.464; California deep firing 0.618 → 0.583, convective
share (CPU day) 0.73 → 0.87, radiative low cover 0.215 → 0.306; ITCZ
convective share 0.871 → 0.951 / 0.747 → 0.888, firing 0.448 → 0.539,
fired tops above 300 hPa 0.024 → 0.009 / 0.143 → 0.203, Q1R centroid
725.8 → 738.8 / 757.9 → 765.9 hPa, T − Jordan 516 / 439 / 848 +0.95 /
+0.99 / −1.48 K; warm pool share 0.744 → 0.850, centroid 682.7 → 692.2
hPa; wettest ITCZ cell 46 → 24 mm/d; replayed global rain 1.962 → 1.944
mm/d; GPU day 186 ASR − OLR 13.1 → 11.9, SWCRE −43.3 → −44.5 W/m². The
trade, deck and California plumes' clouds are 200–500 hPa deep: the type
cannot keep them shallow. It stays on: it is the IFS's rule, it moves its
own observations (trades firing, low cover) the right way, if little.

**The IFS entrainment (Oct 2).** `plumeEntrainmentLaw` 'ifs' (the
default; 'gregory' bit for bit on the CPU): ε = 1.75e-3 /m (1.3 − RH)
(q_s(T)/q_s(T_base))³ where the layer below is buoyant, δ = 0.75e-4 /m
(1.6 − RH), the w² drag 1 + βC_d = 1.949 with mixing rate ε or δ, the
mass flux growing by exp((ε − δ)Δz) where buoyant and falling by
exp(−δΔz) min(1, (1.6 − RH)√(w²_above/w²_below)) elsewhere (IFS Cy43r1
eqs 6.7, 6.8, 6.10, 6.12; Bechtold et al. 2008; all IFS tuning). Not
built: the IFS test parcel's 0.4 ε and its 50 % condensate removal. On
Jordan's column: top 211 hPa, CAPE 209 J/kg; ε 5.7, 5.0, 4.4 e-4 /m in
the first cloud layers. Three days against the type (eight64 / ten64):
trades deep firing 0.698 → 0.711, convective rain 1.66 → 1.62 mm/d, low
cover 0.133 → 0.199; SE Pacific convective share 0.76 → 0.72 (CPU day
0.72 → 0.71), low cloud 0.464 → 0.450; California deep firing 0.583 →
0.594, convective share 0.87 → 0.85, radiative low cover 0.306 → 0.387;
ITCZ convective share 0.951 → 0.871 / 0.888 → 0.891, firing 0.539 →
0.426, fired tops above 300 hPa 0.009 → 0.010 / 0.203 → 0.219, Q1R
centroid 738.8 → 739.6 / 765.9 → 714.8 hPa, ITCZ rain 1.77 → 1.70 / 4.57
→ 6.85 mm/d; T − Jordan 516 / 439 / 848 +0.92 / +1.03 / −1.85 K; warm
pool share 0.850 → 0.797, centroid 692.2 → 680.9 / 797.9 → 771.3 hPa;
wettest ITCZ cell 37 mm/d; replayed global rain 1.972 mm/d; GPU day 186
ASR − OLR 10.2, SWCRE −46.2, LWCRE 22.8 W/m². In the trades the plume
now entrains 5.4–6.1 e-4 /m at 850–924 hPa, loses its buoyancy at 705
hPa (−0.47 K) where it was −0.15, and coasts to 600–700 hPa (0.59 of
tops; clouds 200–300 hPa deep 0.61, 300–400 hPa 0.29): the depth counts
the whole layer in which w² vanishes, about 97 hPa there on bl34. In the ITCZ 0.427 of fired plumes are
200–300 hPa deep with tops at 600–800 hPa; the plumes topping above 300
hPa have PCAPE 147 Pa (ten64: 179 Pa, 0.219 of fired). Kept: it moves the
trades' and California's low cloud and the ten64 ITCZ's heating centroid
(−51 hPa) and rain the right way, and does not move the trades' firing.

**The regimes after the type and the entrainment, and what still misses
(Oct 2).** Acceptance (three N=64 days; the reference's spread above):
N Pacific trades deep firing 0.711 (≤ 0.15) and convective rain 1.62 mm/d
(≤ 1.0): not met; trades low cover 0.199 against the reference's 0.037:
met. SE Pacific convective share 0.72 (CPU day 0.71) and California 0.85
(≤ 0.15): not met; their low cloud 0.450 (within the reference's spread
of 0.459) and radiative 0.387 against 0.215: met. ITCZ convective share
0.871 (0.5–0.8): not met (ten64 0.891); warm pool 0.797 (≥ 0.5): met.
Fired ITCZ tops above 300 hPa 0.010 (≥ 0.4): not met (ten64 0.219). ITCZ
Q1R centroid 739.6 hPa against the base's 675 (eight64): not met; ten64
714.8, and the warm pool's 680.9 / 771.3. T − Jordan at 439–516 hPa
+0.92 / +1.03 K (≤ +1.3): met; 848 hPa −1.85 K against the reference's
−1.91: met. No ITCZ cell above 150 mm/d (37): met. Global rain (replayed
day) 1.972 mm/d. What is in the way of the plume's depth, from the
buoyancy profiles: the entraining plume keeps 0.3–1.0 K of buoyancy below
700 hPa and none above about 650 hPa, against +1.5–1.8 K for the same
plume unentrained (liquid-only thermodynamics, all condensate held as
loading or rained at 3e-3 /m); the plumes that reach above 300 hPa are
those with PCAPE above about 100 Pa. The ice in the plume (element 5:
on Jordan CAPE 177 → 222 J/kg, top 247 → 210 hPa, from the fusion heat of
the condensate frozen above the 0 °C level near 560 hPa) and the plume's
own microphysics (element 6) act exactly in the 300–650 hPa layer where
the fired plumes now stop; they come next. In the trades and the decks the
plume's cloud is deeper than 200 hPa because the sounding has no trade
inversion near 800 hPa (trades h*/c_p 335.0 / 335.9 / 336.1 / 335.1 /
336.9 K at 893 / 850 / 789 / 705 / 608 hPa, no rise): neither element
stops convection there, and the old 120 J/kg threshold had.

**Cost (Oct 2).** Under the exclusive lock, 128 GPU steps from nine64 / nine128_day0183,
alternated twice, merged head (388f6b6) against the type and the
entrainment (16529fa): N=64 step median 28.77, 28.76 → 28.50, 29.86 ms,
the adjust group 5.75, 5.76 → 5.44, 5.64 ms; N=128 123.03, 122.64 →
122.12, 129.64 ms, adjust group 20.72, 20.69 → 19.97, 21.21 ms (the
second pair's physics group, which neither element touches, also rose
38.02 → 39.46 ms). The deep column no longer runs the shallow plume
beside the deep one, which pays for the IFS entrainment's saturation
calls: no measurable cost; 122–123 ms is 62.8 s per model day at N=128.

**Review of the merge and the regime elements (Oct 2).** The merge's
seven conflicted files, rebuilt with git merge-tree, are the only files
where de82fce differs from the automatic merge; every line either side
added is kept or combined there, and 74d7ecb's default digest under the
Rayleigh top, which the merge had stopped asserting, reproduces
(52b46a99). On eight64_day0183 three CPU steps of the full model with
`convectionType` 'top' and `plumeEntrainmentLaw` 'gregory' equal 3f3d915
bit for bit, and with 'gregory' alone 115389e. The IFS formulas at every
cloud layer of Jordan's column and of the same column with its free
troposphere ×0.6 and ×0.3 equal the code exactly (ε 1.1e-8 to 7.3e-4 /m,
δ 5.0e-5 to 1.0e-4 /m); the ×0.3 plume is deep with a 203 hPa cloud
whose w² vanishes inside its 756–826 hPa top layer. On the three-day
state, a day of moist steps keeps each column's water and c_pT + Lq to
1.1e-15 / 1.3e-15 on the CPU and 2.3e-7 / 1.4e-7 on the GPU (the
boundary layer's mixing of θ_l and the dry adjustment are not part of
it). The full models from that state: 104, 153 and 515 of 40962 columns
convect as another type on the two engines at steps 1, 4 and 16 (84 at
step 4 on the previous type and entrainment). The worker-thread engine
parted from the single thread under load (parallel.test, from f6cae8b):
the cells' adjust, which reads the subcloud wind since the Bechtold
closure, shared a phase with the edges' momentum mixing; they are now
two phases. A day of CPU steps from the three-day state, with a twin
of each column every eighth step: deep firing 0.71 in the N Pacific
trades, 0.38 in the S Atlantic trades, 0.19 / 0.24 over the northern /
southern 35–55° oceans, 0.16 over tropical land; with the cloud's top
at the height where w² vanishes inside the top layer instead of that
layer's top interface the fired deep plumes' share would be 0.51, 0.005,
0.10 and 0.13 (ITCZ 0.42 → 0.33, warm pool 0.57 → 0.51), and the IFS's
own test parcel (0.4 ε, half the condensate removed), not built, would
deepen its clouds instead. Nine cells rain more than 150 mm/d, all land,
with convective shares 0.00–0.11.

**The type of convection by the IFS's test parcel (Oct 2).**
`convectionType` 'testParcel' (the default; 'cloudDepth' the previous type
bit for bit on the CPU): the IFS's first-guess deep updraught (Cy43r1 §6.4,
eqs 6.18–6.21 with the w² equation 6.10) types the column before the full
ascent. A test parcel of the deep source's s_l and q_t (the lowest 50 hPa
with the eq. 6.19 excess) leaves the source's top interface at 1 m/s and
mixes toward each layer's air at ε = 0.4 · 1.75e-3 /m (q_s(T)/q_s(T_lowest))³,
below and in its cloud, its w² following the IFS form with that mixing rate,
and keeps half its condensate at each upper interface. Its cloud runs from
the lower interface of its first cloudy layer to the height inside a layer
where its w² vanishes, solved exactly for that layer's buoyancy and mixing
and placed in ln p. Deeper than 200 hPa, the column runs the deep plume
alone; otherwise the shallow plume alone, without the deep ascent. One
departure level only: the IFS repeats the test from higher levels up to
350 hPa above the ground, and the model's deep plume has one source. On
Jordan's column the test cloud passes 200 hPa at 756 hPa and its w² vanishes
at 124 hPa (the plume tops at 211); with the free troposphere ×0.3 and ×0.1
the test parcel is still deep (w² zero at 215 and 266 hPa) while the plume
tops at 756 and 659 hPa: the test parcel's entrainment does not depend on
the environment's humidity. The trade-wind and trade-cumulus columns are
shallow (test clouds 939–833 and 958–808 hPa). Three N=64 days from
eight64_day0183 against the regime elements' reference (cxb; the reference's
replicate spread in brackets): N Pacific trades typed deep on 0.950 of the
column-steps, deep firing 0.711 → 0.749 (0.004), convective rain 1.62 →
1.58 mm/d, low cover (cloudRegimes) 0.199 → 0.306; SE Pacific convective
share (8 audit steps) 0.72 → 0.77, CPU day 0.71 → 0.73, low cloud 0.450 →
0.446 (0.03); California deep firing 0.594 → 0.608, convective share 0.849
→ 0.836, radiative low cover 0.387 → 0.413; ITCZ typed deep 0.703, firing
0.426 → 0.545, convective share 0.871 → 0.881 (0.001), fired tops above 300
hPa 0.010 → 0.010, Q1R centroid 739.6 → 743.0 hPa; warm pool share 0.797 →
0.829; replayed global rain 1.972 → 1.979 mm/d; GPU day 186 SWCRE −46.3,
LWCRE 22.8 W/m². From ten64_day0183: trades firing 0.554 → 0.653, ITCZ share
0.891 → 0.876. Acceptance: trades firing (≤ 0.15) and convective rain
(≤ 1.0 mm/d), the SE Pacific and California shares (≤ 0.15) and the ITCZ
share (0.5–0.8) not met; the low-cloud guards and the warm-pool share met.
Deep firing over a CPU day from the day-186 state (twin map, ce/rv/mapday.mjs):
35–55° oceans 0.195 → 0.177 (north) and 0.237 → 0.211 (south), S Atlantic
trades 0.376 → 0.300, tropical land 0.156 → 0.098.

What the trades' sounding lacks. The test parcel types 0.95 of the trade
column-steps deep because nothing near 800 hPa stops a weakly entraining
parcel: the box's T − Jordan is −2.2 K at 850–924 hPa, −1.4 at 705 and +0.6
at 608 hPa, its humidity falls from 0.95 at 924–946 hPa to 0.66 at 705 and
0.46 at 608 hPa, a decline spread over 250 hPa with no temperature jump,
where the observed trade inversion is a 2–5 K jump in θ over 100–300 m near
850–800 hPa (BOMEX, ATEX). The budget of the 800–900 hPa layer (the 850 and
893 hPa layers, W/m²): longwave −24.4, shortwave +9.0, dynamics −5.7, the
boundary layer's mixing −8.0, condensation +4.3, the deep plume's rain
+23.1, its transport +13.0, its downdraft's evaporation −6.8, the shallow
plume −4.4, recondensation −2.4, rain evaporation −1.6. The deep plume
heats the layer by 2.5–3.2 K/day and the dynamics cool it (−0.5 K/day, net
ascent), where the observed trade-cumulus layer is warmed by subsidence and
cooled and moistened at its top by the shallow cumulus' detrainment. The
inversion is made by the subsidence above the shallow plume's detrainment
and cloud-top cooling; here the deep plume's own heating replaces the
subsidence warming, and the bl34 layers there (789, 850, 893 hPa) are
40–60 hPa (about 500 m) thick, wider than the observed inversion. The
element stays on: it is the IFS's procedure, and the wrong-way move of the
trades' firing (+0.04) comes from the sounding it is given, not from the
test; 'cloudDepth' gives the previous type.

**Ice in the plumes and the melting of their frozen rain (Oct 2).**
`plumePhase` 'mixed' (the default; 'liquid' the previous plumes bit for bit
on the CPU), element 5 of the specification (IFS Cy43r1 §6.6.2 and §6.6.6
on the model's own linear 235.15–273.15 K ramp): both plumes and the test
parcel carry s_li = c_p T + g z − L q_l − (L + L_f) q_i, saturated over
cloudSaturation's mix, q_i = (1 − α(T)) l, T by four Newton steps (within
1.5·10⁻⁶ K of bisection over 180–900 hPa and 205–290 K). A rain's ice share
takes L_f more into s_li, falls as its own stream and melts in the first
layer below it at 273.15 K or warmer (the IFS relaxes toward 0 °C over a
few layers, eqs 6.47–6.48); its evaporation below cloud base takes L + L_f;
what reaches the ground frozen is `convectiveSnow`, to which the surface adds
no fusion heat (and from which it takes L_f where the lowest air is not
freezing); the downdraft evaporates only the liquid and melted rain. The
flux form carries s_li, so frozen condensate that detrains returns to the
environment's L-only convention in the detraining layer, and column
c_p T + L q stays exact (a 265 K column whose plume snows: L_f times the
precipitation under either phase, to 2.3·10⁻¹⁶). On Jordan's column the deep
plume's top rises 211 → 153 hPa and its CAPE 209 → 311 J/kg, its frozen
rain melting at 610 hPa (275.7 K). Three N=64 days from eight64_day0183
against the test parcel's run (ten64_day0183 in brackets): fired ITCZ tops
above 300 hPa 0.010 → 0.025 (0.217 → 0.277), warm pool 0.198 → 0.288
(0.100 → 0.149); the 0 °C layer's melting (607 hPa) −0.02 K/day in the ITCZ,
0.013 K/day per mm/day of deep rain (−0.32 and 0.049; warm pool −0.19 and
0.042), the freezing's heating 0.01–0.14 K/day near 373 hPa: the plume rains
most of its condensate below the freezing level at 3·10⁻³ /m, so little
freezes; ITCZ high cover 0.413 → 0.407 (0.382 → 0.347), warm pool 0.367 →
0.363 (0.301 → 0.291); T − Jordan at 607 / 516 hPa 0.1 / 0.9 → 0.1 / 0.9 K
(0.3 / 0.5 → 0.0 / 0.4); ITCZ convective share 0.881 → 0.903 (0.876 →
0.952), Q1R centroid 743.0 → 744.1 hPa (708.7 → 708.5), warm-pool centroid
688.5 → 671.8 (779.8 → 754.1); trades deep firing 0.749 → 0.759; replayed
global rain 1.979 → 1.989 mm/d; GPU day 186 SWCRE −46.0, LWCRE 22.5 W/m².
Acceptance: the tops (≥ 0.6 or +0.08) not met on eight64 (+0.015) nor
ten64 (+0.060), met for the warm pool (+0.090); the melting (0.1–0.3 K/day
per mm/day) not met, an order of magnitude weaker; high cover not lower:
not met by 0.004–0.035; the 516–600 hPa bias: met.

**The IFS updraught conversion (Oct 2).** `plumeConversion` 'sundqvist'
(the default; 'zhangMcFarlane' the previous 1 − exp(−3·10⁻³ /m Δz) bit for
bit on the CPU), element 6 of the specification (IFS Cy43r1 eqs 6.38–6.40
after Sundqvist 1978; IFS tuning, which the IFS says probably still
overestimates the updraught condensate): where the plume's condensate l at
an upper interface exceeds 0.3 g/kg over sea or 0.5 g/kg over land,
l (1 − exp(−a Δz)) rains, a = c0/(0.75 w)(1 − exp(−(l/l_crit)²)),
c0 = 1.4·10⁻³ /s (1.3 α + 1 − α) on the model's phase ramp, l_crit
0.5 g/kg, w the plume's speed there within 1–10 m/s, and below 268.16 K c0
times and l_crit over 1 + 0.5 √min(268.16 − T, 18). What the plume keeps
detrains where its mass flux falls. The in-updraught fallout (eqs
6.41–6.42) is not built. On Jordan's column the rain of each layer equals
the analytic integral to 1.9·10⁻¹⁶; the detrained condensate (Σ
max(0, M_k+1 − M_k) times the condensate carried in) is 0.065 of the rain
made, centred at 346 hPa (0.002 at 378 hPa before; the specification's
single column on the previous plume expected 0.15–0.35). Three N=64 days
from eight64_day0183 against the mixed-phase run (ten64_day0183 in
brackets): detrained condensate per unit of convective rain, ITCZ 0.327 at
676 hPa, 0.023 of it above 400 hPa (0.187 at 566 hPa), warm pool 0.167 at
602 hPa (0.260 at 630); convective share, ITCZ 0.903 → 0.711 (0.952 →
0.833), warm pool 0.897 → 0.681 (0.988 → 0.955); stratiform share (melted
falling ice and conversion above 700 hPa) ITCZ 0.071 → 0.203 (0.010 →
0.070), warm pool 0.064 → 0.213 (0.002 → 0.017); Q1R peak bin ITCZ 925 →
925 hPa, warm pool 775 → 775 hPa; Q1R centroid ITCZ 744.1 → 730.5 (708.5 →
706.6), warm pool 671.8 → 665.5 (754.1 → 778.0) hPa against spreads of 0.4
and 0.25; warm-pool high cover 0.363 → 0.361 (0.291 → 0.307), its thin
share (τ < 3.6) 0.44 (0.90); ITCZ high cover 0.407 → 0.417; trades low cover
0.282 → 0.318, their convective rain 1.62 → 1.25 mm/d, deep firing 0.759 →
0.735; SE Pacific convective share (audit) 0.74 → 0.60, California (CPU day)
0.826 → 0.696; fired ITCZ tops above 300 hPa 0.025 → 0.034; replayed global
rain 1.989 → 1.967 mm/d; GPU day 186 SWCRE −47.6, LWCRE 23.2 W/m².
Acceptance on eight64: the detrained share (0.15–0.4) met, its centroid
above 400 hPa not met (676 and 602 hPa); the convective share (0.55–0.75)
met; the stratiform share (≥ 0.25) not met (0.20–0.21); the Q1R peak at or
above 600 hPa not met; the centroid's rise beyond the spread met; the warm
pool's high cover not up (within its 0.002 spread), its thin share (≥ 0.38)
met; the trades' low cover (≥ 0.18) met. The detrainment centroid lies low
because the plumes still stop at 500–700 hPa: what the plume carries
detrains where its mass flux falls, and it falls there.

**The ice fall at 3.29 and homogeneous nucleation, tested on the anvil
source (Oct 2).** Three N=64 days from eight64_day0183 on the tree with
elements 5 and 6, the regimes one step on (cloudRegimes.mjs; the reference's
replicate spread of high cover and its thin share is at most 0.002 / 0.01):

| day 186 | defaults (fall 2.5) | `iceFall` 3.29 | `iceNucleation` | both |
|---|---|---|---|---|
| global high cover; thin share | 0.278; 0.49 | 0.249; 0.46 | 0.267; 0.47 | 0.240; 0.45 |
| warm pool high cover; thin share | 0.361; 0.44 | 0.333; 0.43 | 0.358; 0.43 | 0.326; 0.42 |
| ITCZ high cover; thin share | 0.417; 0.30 | 0.377; 0.31 | 0.414; 0.28 | 0.379; 0.31 |
| RH 150–350 hPa over water / ice: global; warm pool; ITCZ | 0.39/0.59; 0.48/0.70; 0.55/0.80 | 0.39/0.59; 0.49/0.71; 0.54/0.80 | 0.40/0.60; 0.49/0.71; 0.55/0.81 | 0.40/0.60; 0.50/0.73; 0.55/0.80 |
| upper-tropospheric layer area above ice saturation (of it clear) | 0.01 (0.00) | 0.00 | 0.02 (0.69) | 0.02 (0.77) |
| GPU day 186 SWCRE; LWCRE | −47.6; 23.2 | −46.4; 20.6 | −47.4; 22.2 | −46.2; 19.6 |

Against the references in the regime table (ISCCP and CALIPSO: global high
0.2–0.3 with about 0.6 of it thin, warm pool high 0.55–0.70 with about 0.5
thin, ITCZ 0.45–0.60), the fall coefficient of 3.29 lowers the high cover
where it is already short (warm pool −0.028, ITCZ −0.040, global −0.029) and
thins nothing; nucleation leaves the high cover within 0.011, lets 0.02 of
the upper-tropospheric layer area stand supersaturated over ice, mostly
clear, and raises the RH over ice by 0.01. The docs give no observed
upper-tropospheric humidity to judge either by. Neither change is
supported by the regime observations; the defaults stay at 2.5 and off. The
anvil source of element 6 is too low to test them fairly: its detrainment
centres at 600–680 hPa, below the cirrus levels where the fall speed and
nucleation act.

**The type, the ice and the conversion together (Oct 2).** Cost under the
exclusive lock, 128 GPU steps from nine64 / nine128_day0183, alternated
twice, 49b2ceb against 994f72f: N=64 step median 29.39, 29.25 → 30.61,
30.64 ms (+4.4 %), the adjust group 5.53 → 6.94 ms; N=128 120.98, 119.87 →
125.86, 123.11 ms (+3.0 %), the adjust group 19.73 → 23.47 ms, 61.7 → 63.7
s per model day. The specification's budget for elements 5 and 6 was 1.3 %;
most of the rest is the mixed-phase saturation, two exponentials in every
Newton step below 273.15 K, and the test parcel's ascent. Deep firing over a
CPU day from the day-186 state of element 6's run against the reference: N
Pacific trades 0.711 → 0.730, S Atlantic trades 0.376 → 0.335, SE Pacific
0.086 → 0.107, California 0.594 → 0.601, 35–55° oceans 0.195 / 0.237 →
0.188 / 0.227, ITCZ 0.426 → 0.559, warm pool 0.571 → 0.664, tropical land
0.156 → 0.093. The full suite (64 files, concurrently) passes; two N=6
coupled parity cases (day means, rain accumulation) run on the previous
type, phase and conversion, with the measured reasons in 218e104.

**Review of the type, the ice and the conversion (Oct 2).** Re-run from
eight64_day0183, the three days, tropicalHeating.mjs and cloudRegimes.mjs
reproduce every number above (replay mismatches 0 of 722432). Three CPU
steps from that day-186 state under convectionType 'cloudDepth', plumePhase
'liquid' and plumeConversion 'zhangMcFarlane' (and under each later pair)
give the digests of 49b2ceb, 8e10cff and e76c952 bit for bit.
- The test parcel against the IFS first guess (Cy43r1 §6.4) integrated by
  hand in 20–40 sub-steps a layer, its base at its first cloudy height:
  Jordan's column deep with w² vanishing at 107 hPa (the code's discrete
  form 106), ×0.3 at 161 (159); a trades sounding built from the model's
  box (T − Jordan −2.2 K below 850 hPa, RH 0.95 → 0.46 over 924–608 hPa)
  without the surface excess shallow, its parcel stopping at 949 hPa below
  its condensation level. On a CPU day from the day-186 state the hand
  parcel and the code type alike on 0.952 (ITCZ), 0.977 (warm pool), 0.995
  (N Pacific trades, 0.964 deep by hand), 0.923 (SE Pacific), 0.942
  (California) and 0.92–0.96 (35–55° oceans) of the column-steps; the GPU
  types, tops and base fluxes match the CPU's on all 362 columns of
  Jordan, Jordan ×0.3, both trades soundings and a cold-based column.
- The mixed-phase plume air at 268, 254 and 238 K from guesses 15 K below
  to 8 K above: CPU within 2.1·10⁻⁶ K of bisection, GPU within 1.2·10⁻⁵ K;
  T(s_li) continuous through 273.15 and 235.15 K to 2·10⁻⁹ K (CPU) and
  3.5·10⁻⁵ K (GPU, single precision).
- The conversion equals eqs 6.38–6.40 by hand on every layer of Jordan's
  column to 5·10⁻¹⁶; a Δz there is 0.3–1.3 in the bl34 layers. Converting
  the condensate at the layer's top over the whole layer keeps a steady
  condensate of aΔz e^(−aΔz)/(1 − e^(−aΔz)) = 0.47–0.86 of the IFS's
  analytic solution with the condensation source (l = l₀e^(−az) + b/a (1 −
  e^(−az)), §6.6.3); with that solution Jordan's detrained condensate per
  unit of rain made is 0.097 (net mass-flux decrease) or 0.189 (with the
  turbulent detrainment δM of the buoyant layers) against 0.065 and 0.126.
  The tropicalHeating.mjs detrainment counts the net decrease only.
- The condensate the deep plume detrains (a CPU day, twin every eighth
  step): the condensation that follows the plume keeps 0.278 (ITCZ) and
  0.293 (warm pool) of it as cloud, 0 where the layer's RH was below 0.7,
  0.12–0.15 at 0.7–0.9, 0.53–0.60 above 0.9; after the step's conversion
  to rain and ice fall the cloud has gained 0.008 and 0.022 of it.
- The fusion heat: where the cloud base lies above the 0 °C layer the
  frozen rain that sublimated below the base, or that the downdraft took,
  was melted again at the 0 °C layer. Over 32 CPU steps from the day-186
  state 8248 column-steps (high-latitude land and sea) were off by up to
  5.5·10⁻³ kg/m² of L_f a step (global mean −1.6·10⁻³ W/m²); fixed in
  8615f9e (16 column-steps left, the same non-convective cells as under
  plumePhase 'liquid'); on 362 cold-based columns 1.3·10⁻⁶ → 5.8·10⁻¹⁶
  relative on the CPU, 9.8·10⁻³ → 4.0·10⁻⁴ kg/m² of L_f on the GPU. Water
  exact to 4·10⁻¹⁵ per column-step with the filler's loss counted.
- Engine parity of the full models from the fixed tree's day-186 state:
  columns of a different convection type 123 / 212 / 577 of 40962 after 1 /
  4 / 16 steps (81 / 118 / 507 under the previous elements), θ rms relative
  2.5 / 5.7 / 13.5·10⁻⁶. The fixed tree's three days give the numbers above
  within the trajectory's noise (ITCZ convective share 0.714, warm pool
  0.682, stratiform share 0.200 / 0.214, detrained 0.323 at 675 hPa and
  0.166 at 602 hPa); the full suite passes (64 files, 623 tests).
- The frozen share of the deep plume's rain, from the column budgets:
  0.05 (ITCZ) and 0.13 (warm pool); the specification's 0.1–0.3 K/day per
  mm/day at the 0 °C layer implies 0.15–0.45.

**Spin-up eleven's tree (Oct 2).** Branch integrate-c at 276d936 (the
surface by class, the land, the surface layer by roughness, the cloud
optics, the spectral gases, the model top, the GMTED terrain drags, bl36
with its lid friction and longwave table) merged with sweep2 at e0d95f0
(the convection-cloud elements 0–6) in 438491c; 306badc settles two parity
checks, 9b2141a gives scripts/longwaveOverlap.mjs the level set's table.
- The merge: createRadiation takes `longwaveOverlap` beside `longwaveTable`,
  the overlap's chain running over the level set's table; cloudEffect and
  dayMeans leave out the union of both sides' parted columns. At N=4 every
  276d936 digest reproduces under its own options, the spectral-gas ones with
  `longwaveOverlap` 'random' (5d1fa897, fc689d9b, 5108787c, 8cb36fdf,
  ce584ac4, e6a0e8e2); e0d95f0's ran from 0063c54's fresh start and are
  re-pinned (the moist defaults c3d3550c, elements 2–6 under the Rayleigh
  top 7545db91, b131be14, d36518ff, ee85badd, then 6baaa821, 9155cbcb,
  6174bf4c, bf86f179, a8f4d734, 1fcfa81c). From e0d95f0's own N=4 fresh
  state the merged tree gives its Rayleigh-top digests bit for bit
  (52b46a99, a64fbb13 with the random overlap, c70751a7, 6761a720, c5362ea5,
  6c479c25); its defaults digest 1dbd465a does not reproduce with the
  gravity waves at their older values (`speedStep` 4, `sourceDescent` false,
  `lidPressure` 0: e7cec906), the flux being scaled after the spectrum's sum
  since b6e8854: after one CPU step from eight64_day0183 2046 of 8478930
  state values differ, by at most 7.1·10⁻¹⁵ (3.9·10⁻¹⁴ relative).
- Suite: 64 files, 567 tests, concurrently; 62 files passed at once.
  dayMeans' per-cell albedo (0.00223 against 10⁻³ over all columns) and
  gpuModel's deckRest 'regime' water (rms 3.0·10⁻³ against 10⁻³; one column
  topping at 84404 Pa on the CPU and 75604 Pa on the GPU, 1.6·10⁻⁴ over the
  other 361; 1.7·10⁻⁴ and 2.0·10⁻⁴ on the parents) now leave out the
  columns whose shallow plume parted (at most 2 % asserted).
- Proofs, CPU, ocean off: three N=64 steps from eight64_day0183 remapped to
  bl36 with `capeClosure` 'threshold', `pcapeBoundary` 'signed',
  `excessVelocity` 'convective', `plumeSourceDepth` 'boundaryLayer',
  `cumulusClosure` 0.06, `convectionType` 'top', `plumeEntrainmentLaw`
  'gregory', `plumePhase` 'liquid', `plumeConversion` 'zhangMcFarlane' and
  `longwaveOverlap` 'random' give 276d936's state digests bit for bit
  (f10554f5, 18df3a04, 621da31e); three steps from eight64_day0183 (bl34)
  with `gravityWaves` false give e0d95f0's (e386cf8d, 0bf5043a, e7e8c953).
- Closure on the bl36 state (scratch closure script, 40962 columns at four
  times of day): CPU |ASR + reflected − incoming| 4.5·10⁻¹³, layers' SW
  6.8·10⁻¹³, LW 5.4·10⁻¹³, surface flux 9.1·10⁻¹³, column heating 1.4·10⁻¹²,
  vapour 1.4·10⁻¹² W/m²; GPU 1.8·10⁻⁴, 1.2·10⁻³ (LW), 1.5·10⁻⁴, one physics
  kernel's heating 8.4·10⁻⁴ and vapour 1.05·10⁻⁴ W/m². The longwave
  overlap's recomputation repeats the OLR to 6.6·10⁻¹⁶ and the surface
  downward longwave to 5.2·10⁻¹⁶. The moist step alone over 16 steps: CPU
  per column c_p T + L q (with L_f of the convective snow) to 8.4·10⁻¹⁶ and
  water to 8.0·10⁻¹⁶; GPU, one adjust kernel from the same state with L_f
  of the CPU's convective snow counted, c_p T + L q to 6.1·10⁻⁸ and water
  to 1.7·10⁻⁷ in every column.
- Engine parity from that state after 1, 4 and 16 steps: convection type
  (none, shallow, deep by a top above 700 hPa) apart on 2, 3 and 21 of 40962
  columns (8351, 7776, 8566 fired), the same type with tops an interface
  apart on 0, 3, 16, regimes apart on 0, 2, 10; lowest layer T rms 4.0·10⁻⁵,
  6.1·10⁻⁴, 1.4·10⁻³ K; top eight layers' u rms 8.6·10⁻⁵, 2.1·10⁻⁴,
  4.8·10⁻⁴ m/s; gravity-wave acceleration rms 7.5·10⁻⁴ to 1.6·10⁻³ m/s/day
  of up to 25.
- Three one-day GPU segments equal one three-day segment byte for byte
  (N=64, bl36, ocean every 8 steps, STRATOSPHERE=1; the snapshot carries
  `subcloudVirtual`, 491544 values). bl34 states load and step:
  eight64_day0183 for three days (day 186 ASR 246.0, OLR 236.8, rain 1.58
  mm/d, SWCRE −47.7, LWCRE 23.2; e0d95f0 246.1, 236.9, 1.58, −47.6, 23.2;
  276d936 on the bl36 remap 241.3, 233.7, 1.70, −52.4, 25.8),
  eight128_day0183 for 64 steps, finite, largest |u| 76.0 m/s.
- Smoke test of eleven's configuration: `NS="64 128" PREFIX=smoke
  LEVELS=bl36 PER_YEAR=36 KEEP=200 OCEAN='{"everySteps":8}' STRATOSPHERE=1
  scripts/pairedSpinup.sh`, fresh from the atlas (WOA, 29078 and 116471 sea
  cells), STOP_smoke placed while N=64's first segment ran, so the driver
  ended after it (day 10); N=128's first segment then ran as the driver
  runs it with DAYS=2. Both log the land's neutral start with its record at
  0 days, read data/subgrid_N64.bin and data/subgrid_N128.bin and the bl36
  longwave table (34 g-points); clamped 0 on every day, no NaN; the day-10
  and day-2 states reload and step 16 steps finite (largest |u| 96.9 and
  82.1 m/s). Day means:

| N, day | ASR | OLR | rain mm/d | albedo | SWCRE | LWCRE |
|---|---|---|---|---|---|---|
| 64, 1 | 218.4 | 191.8 | 2.91 | 0.359 | −73.2 | 53.3 |
| 64, 2 | 214.9 | 185.7 | 4.63 | 0.369 | −76.8 | 58.3 |
| 64, 3 | 224.1 | 195.9 | 4.44 | 0.342 | −67.6 | 51.0 |
| 64, 4 | 228.8 | 199.5 | 4.46 | 0.328 | −62.8 | 49.6 |
| 64, 5 | 230.9 | 202.2 | 4.41 | 0.322 | −60.7 | 48.4 |
| 64, 6 | 230.6 | 205.1 | 4.21 | 0.323 | −61.1 | 46.8 |
| 64, 7 | 230.0 | 207.6 | 3.85 | 0.324 | −61.7 | 45.2 |
| 64, 8 | 230.1 | 210.9 | 3.54 | 0.324 | −61.6 | 42.3 |
| 64, 9 | 230.6 | 213.8 | 3.37 | 0.323 | −61.1 | 39.8 |
| 64, 10 | 230.5 | 215.1 | 3.30 | 0.323 | −61.1 | 38.9 |
| 128, 1 | 227.1 | 194.6 | 3.37 | 0.333 | −64.5 | 51.4 |
| 128, 2 | 226.3 | 192.8 | 4.78 | 0.335 | −65.4 | 52.7 |

  Top six layers (0.148–14 hPa), largest wind and Courant number
  horizontal/vertical over the days: N=64 75 m/s, 0.24/0.09 (day 10: 47,
  0.15/0.03); N=128 85 m/s, 0.27/0.10. scripts/tropicalHeating.mjs on
  smoke64_day0010 (replay mismatches 0 of 722432; heat closes to 1.1·10⁻¹³
  K and q_t to 2.2·10⁻¹⁹ a step): rain, convective share, firing, Q1R
  centroid, stratiform share: Pacific ITCZ 7.06 mm/d, 0.96, 0.79, 654 hPa,
  0.03; warm pool 7.16, 0.95, 0.60, 552, 0.04; SPCZ 7.73, 0.58, 0.42, 608,
  0.30; N Pacific trades 3.82, 0.94, 0.66, 659, 0.04; Amazon 0.42, 0.45,
  0.04, 751, 0.25; global rain 3.16 mm/d (convective 1.43).
- Pace (Apple GPU of this Mac, exclusive lock, 128 steps after 16 from
  eight64/eight128_day0183 on bl36, ocean every 8 steps, alternated twice):
  profiled median N=64 30.09, 30.07 ms a step against 276d936's 26.67, 26.62
  (+13 %), the physics pass 15.3 against 11.9 ms; N=128 127.27, 126.97
  against 115.19, 115.30 ms (+10 %), the physics pass 58.4 against 46.6 ms.
  128 steps back to back, each awaited as the spin-up awaits them: N=64
  35.97, 35.91 against 32.55, 32.57 ms; N=128 151.27, 151.44 against 139.59,
  139.57 ms, 77.5 s for a model day's 512 steps; a one-day segment from
  eight128_day0183 on bl36 (STRATOSPHERE=1) logs the day at 1.3 min and
  saves day 184 82.4 s after its 8 s of setup, 90.8 s of process wall.
- On by default: the correlated longwave on the level set's table with the
  exponential-random overlap; CLIRAD gases with the AFGL ozone and
  near-infrared Rayleigh; the cloud optics; the PDF cloud cover with
  exponential-random overlap; the mixed-layer deck by regime
  (`deckRest` 'regime'); the moist boundary layer with the implicit drag;
  the surface layer by roughness with the convective gust and land
  humidity by wetness; the land's trees by moisture, the treeline,
  grassland, soil darkening, snow ageing and masking, the neutral start and
  its jumps at days 365 and 730; snow and ice albedo by temperature and
  age; the uniform condensation with ice saturation, the falling ice at
  2.5 and the stratiform lifetime; the Bechtold closure with PCAPE_bl at
  least 0, the deep source from the lowest 50 hPa with the IFS excess,
  Grant's shallow closure (0.03), the test-parcel type, the IFS
  entrainment, the mixed-phase plume with melting and the Sundqvist
  conversion; the eddy sponge from 78 Pa and the GISS lid friction on bl36;
  the gravity waves launched as cg_drag does (source descending with
  latitude, 2 m/s phase speeds, the lid's flux spread); the mountains'
  blocking and wave drag and the turbulent form drag on the GMTED fields.
  Off by evidence: the latitude-dependent gravity-wave flux (`northFlux`,
  `southFlux` 0, `equatorialFlux` = `flux`); the ice fall at 3.29 with
  homogeneous nucleation (`iceFall` 2.5, `iceNucleation` false). Also off: the
  deep plume's momentum transport (`plumeMomentum` false).
- Review (Oct 2): the suite again 64 files, 567 tests, 0 failures, in 293 s
  six files at a time; the ten-day N=64 smoke run through pairedSpinup.sh
  repeats every daily line of smoke64.log; the bit-for-bit proofs against
  both parents, the N=4 Rayleigh-top digests from e0d95f0's fresh state,
  the split-run byte equality, the CPU closure and the engine parity above
  reproduce. e0d95f0's defaults digest under the older gravity waves
  differs after 12 steps at N=4 in 287 of 26568 values, by at most
  2.2·10⁻¹⁵. From smoke64_day0010, which carries `subcloudVirtual` and the
  cumulus cloud, the engines part faster than from the remap: convection
  type apart on 12, 51 and 313 of 40962 columns after 1, 4 and 16 steps
  (276d936 on the same state 1, 9, 40), lowest layer T rms 3.7·10⁻⁴,
  8.9·10⁻⁴, 6.2·10⁻³ K; with `subcloudVirtual` zeroed on both, 1, 40, 225.
- Launch (Oct 2, 20:05 PDT; Oct 3 03:05Z): run eleven started on the Verda
  H100 spot instance gcm-eleven (1H100.80S.30V, FIN-02) at 85539ef, main
  fast-forwarded to it, with `NS="64 128" PREFIX=eleven LEVELS=bl36
  PER_YEAR=36 KEEP=1000 OCEAN='{"everySteps":8}' STRATOSPHERE=1 UNTIL=1095`.
  The instance's benchmark (fresh atlas starts on bl36): N=128 13.19 s a
  model day steady, 18.1 with a segment's setup and finish; N=64 2.41 and
  4.50; the three years project to 7.0 h and $13.28 at 1.8911 $/h with
  58 GB of states on the volume. Its suite (30 runners, node 22, driver
  580.178.04, Dawn on Vulkan): every file passes but gpuModel, whose
  three cloud-cover and overlap parity tests fail by 3.6·10⁻⁴ and
  4.1·10⁻⁴ K/day against largest heatings of 23.6 and 26.5 (limits
  1·10⁻⁵ and 1.5·10⁻⁵ of the largest, the second already at the
  settlement rule's ceiling). The differences sit in the same columns of
  the top-but-one layer as the Mac's largest (1.2·10⁻⁴), two to four
  times bigger, with the same rms over all sunlit layers (1.2·10⁻⁵
  against 0.96·10⁻⁵ K/day) and the same per-layer maxima below: the
  transcendental precision of that driver, not a decision or a defect,
  and the limits stay as they are. Before the launch, the merged tree
  from the remapped eight64 day-183 state at N=64 ran 30 days: ASR − OLR
  +22.7 on day 184, +9.5 on 186, +3.9 on 199 and +0.4 over days 204–213
  (ASR 236.6, OLR 236.2, albedo 0.305, SWCRE −55.4, LWCRE 24.0, rain
  2.64 mm/d), Ts falling 1.3 K over the 30 days on that state;
  `scripts/toaBalance.mjs` folds any log's daily lines into the balance
  over its last day, week, month and year and every whole year. The land
  jump was rehearsed on the 150-day bl36 record with a jump forced at day
  151: land albedo 0.260 → 0.234, trees 0.248 → 0.266, topsoil carbon
  1.6 → 5.2 kg/m², and the run continued; that half-year record is a
  northern summer, so the rehearsal put trees in 0–20N and stripped
  10–50S, which a whole year's record does not.

**The cloud overlay's blinking (Oct 3).** On the page's Cloud cover
overlay low cloud blinked on and off from one step to the next.
`scripts/cloudFlicker.mjs` (from eleven128_day1825, 64 steps of 168.75 s,
and eleven64_day1825, 128 of 337.5 s) counts a blink where a cell's
overlay opacity moves by more than 0.4 and back within three steps: 2.51
% of the cells a step at N=128 and 2.14 % at N=64, 97 % of the cloudy
runs one step long and a lag-1 autocorrelation of −0.95 in the step
changes of the blinking cells, a two-step limit cycle. Its mechanism, from
the adjust kernel's stages: a layer near 800–850 hPa at 0.95–0.99
humidity above a shallow mixing top (median 0.75 km) holds the uniform
distribution's condensate (50–140 g/m²); the next step's diagnosis finds
it as a cooling cloud top, the regime goes decoupled or coupled and the
mixing top rises to it (median 2.5 km); the mixing returns its water to
vapour and the layers below the mixing top adjusted to grid saturation,
which at 0.96 is none (99.8 % of the vanishings); with no cloud the top
falls back and the uniform distribution condenses it again. The rain's
conversion emptied a layer in 0.004 % of them and the plume and transport
in under 0.1 %, so neither a limit on the conversion nor a lifetime on
the cloud would stop it: a cloud given a lifetime would decay inside the
mixed layer and pop back above it on a slower cycle. A second, smaller
cycle (1.5 % of the blinks, 96–99 % of those in the Peru and Namibia
boxes) was the deck's: its carried height h rested within 1–2 m of a
layer's midpoint, and the slab of whole layers whose midpoints lie below
h took that layer in on one step and out on the next, the water path
going between 150 and 0 g/m² with the gate above 0.6 throughout.
Each cloud type's blinking on its own overlay is in the per-type table
at the end of M22 (`types` in the script's JSON).

The scheme. `boundaryCondensation` 'cloudLayer' (the default here, until
'uniform' replaced it the same day, below): the run of
cloudy layers that makes a column cloud-topped (the cloud top's layer
and the cloudy layers below it whose cooling the diagnosis sums,
`cloudLayer`, PH `CLOUDK`, diagnosed each step and not saved) holds the
uniform distribution's condensate as the free troposphere does, and the
other mixed layers adjust to saturation as before. The distribution is
ECHAM6's Sundqvist one with RH_c from crs 0.975, crt 0.75 and nex 2
(tuned at T63, not observed); no new parameter. Neither source draws the
line at the cloud-top run: ECHAM6's scheme, and the Unified Model's
large-scale scheme under Lock et al. (2000)'s mixing, condense every
layer of the boundary layer, which is 'uniform' here, so 'cloudLayer' is
this model's restriction of the published scheme to the layers whose
cloud the boundary layer's diagnosis reads. 'uniform' gives every mixed
layer that distribution, 'saturation' is the scheme before. 'uniform'
is the arrangement the sources describe and the default by this model's
rule (below); under it five parity tests failed at their limits, and not
through an engine difference. Stage by stage (`scripts/adjustStages.mjs`
on eleven64_day1825, both engines given the CPU physics phase's state
in single precision and its boundary-layer fields): after the mixing and
after the condensation no layer of the 355,944 below the mixing top or
the 1,118,688 above parts by 10⁻³ K in θ or 10⁻⁷ in q or qc (qc within
2.7·10⁻⁸); after the plumes and after the whole adjust step 31 and 33
of 40,962 columns part in θ and 824 plumes' base flux by 1 % (43, 44
and 849 under 'cloudLayer'); the boundary layer's diagnosis of the
adjusted state parts in no regime and no mixing top by 1 m over 23,903
cloud-topped columns (one of each under 'cloudLayer'). The physics
kernel alone on N=6 states stepped 12 and 18 steps under 'uniform' (295
and 674 cloudy layers below the mixing top) keeps OLR within 7·10⁻⁴
W/m² and the longwave heating within 5·10⁻⁴ K/day, as under
'cloudLayer'. What parts the engines over many steps is the scheme's own
response to their single-precision differences, which reach 10⁻⁴–10⁻³
K of θ in the lowest layers within a few steps. On the rain-accumulation
test's setup (`scripts/perturbedRain.mjs`) the CPU against itself with
±10⁻⁴ K of noise on θ before every step parts the rain of 4 to 23 of
362 cells by 10⁻³ of the largest cell's under 'uniform' over six noise
seeds (20 with the script's default seed) and of 0 to 9 under
'cloudLayer'; ±3·10⁻⁵ K parts none with three seeds of four and 13
with the fourth, all of them among the GPU's 20; one ulp of θ, or
rounding the state to single precision every step, parts none; ±10⁻³ K
once parts 16 and 1. 18 of the GPU's 20 cells part under the CPU's
noise with some seed, and 14 lie within two cells of cells 90 and 223,
where the deep plume fires on one engine only at steps 19 and 20. With
the default seed 11 of the 20 follow a discrete decision (7 a plume's
base flux near its onset, 1 a plume firing, 2 a layer's cloud at the
cloud-top threshold, 1 a merge of the dry adjustment) and 9 none of
their own, each beside a cell whose decision parted: 8 of them part in
the large-scale rain, which every mixed layer's distribution condensate
feeds at each step, by 1.1–7.4·10⁻³ of the largest cell's, one in the
convective rain by 1.1·10⁻³. The
other failures: the cloud effects test finds 4 columns whose dry
adjustment merges apart against a bound of 3.6 and leaves out 37
columns against 32.6; 5 parted columns move the day-mean absorbed
sunlight by 1.6·10⁻⁴ of itself against 10⁻⁴; the treeline's season
length parts by 5.4·10⁻² against 10⁻⁴; and in the stratiform-lifetime
test the engines agree exactly (263 and 263 layers) but the 3 h
lifetime keeps more cloud in 9.6 % of the 2,727 cloudy layers against
its floor of 10 % (Oct 3). Making 'uniform' the default needed those
bounds restated against the CPU's response to perturbations of the
engines' size (below). Every cloudy
layer below the mixing top keeps the variance cover in the radiation,
so the run's condensate is the distribution's and its cover is not: at
N=64 after 8 steps from eleven64_day1825, over the 62,637 run layers
below the mixing top in 22,933 columns, the radiation's cover averages
0.61 (weighted by condensate; 0.1 % at the cover floor) against the
distribution's own 0.71, an in-cloud water of 0.20 against 0.16 g/kg.
`deckSlab` 'fraction' (the default): the deck's slab takes the layer its
height lies in by the share of that layer's height below h, so that the
slab moves with h continuously, except the inversion ceiling's layer,
which stays the free troposphere as before; 'midpoint' is the slab
before. Both conserve θ_l and q_t layer by layer (the condensation is
the existing linearized step). Under 'saturation' and 'midpoint' both
engines reproduce the parent bit for bit (GPU: 16 steps from
eleven64_day1825, state and frame digests; CPU: 24 steps at N=6).
Tested and not kept: the Gaussian condensate of the variance cover
(Sommeria and Deardorff 1977) below the mixing top, 1.37 % at N=128
with a new cycle where the surface-driven top moved between 0.26 and
1.25 km, and 4.6 mm/d of rain and SWCRE +8 W/m² on its first day;
'uniform' with the midpoint slab, 0.034 % at N=128 but SWCRE −3.7 W/m²
over the three days below against −3.0 under 'cloudLayer' ('uniform'
became the default after, below).

After, the same runs: 0.010 % a step at N=128 (934 blink onsets against
246,813) and 0.006 % at N=64 (321 against 108,484), 0.01 % and 0.10 %
at the page's cadence (2.29 and 0.35 before); cells blinking at least
once 0.06 % and 0.24 % (6.8 and 12.3); the lag-1 autocorrelation of the
step change +0.06 and +0.15 over all cells (−0.92 and −0.87), −0.60 and
−0.38 over those still blinking; one-step cloudy runs 56 % and 38 % of
511 and 169 (97 % and 94 % of 128,778 and 56,121); of the blinks left
93 % and 73 % the deck's. Under 'uniform' with the midpoint slab 0.034
% and 0.11 % were left, the deck's in 99 % (Namibia 10.9 % at N=64),
which 'fraction' took to 0.16 %. Three GPU days at N=64 from
eleven64_day1825, days 1826–1828, before → after: albedo 0.318, 0.314, 0.311 → 0.326, 0.323, 0.321; SWCRE
−59.3, −57.7, −56.6 → −61.8, −60.7, −60.0 W/m²; LWCRE 29.2, 29.6, 28.8 →
29.5, 29.9, 29.2; rain 2.77, 2.77, 2.70 → 2.79, 2.78, 2.73 mm/d; total
cover at the end of day 1828 (`scripts/cloudRegimes.mjs`, one CPU step)
0.57 → 0.60, SE Pacific 0.28 → 0.35, Peru 0.33 → 0.43, Namibia 0.37 →
0.50, the trades 0.52–0.55 → 0.60–0.62. A replicate (the cover floor
moved by 10⁻⁷) stays within 0.1 W/m² of the run before, so the −3.0
W/m² mean in SWCRE and +0.3 in LWCRE are the fix's: the cycle had held
these clouds on half the steps, mostly over the tropical and subtropical
seas (−2.17 of the −3.43 W/m² of day 1828 from the seas within 30° of
the equator, −5.8 W/m² there). A step costs the same within the
timing's noise (162.2–164.4 ms at N=128 under either, two runs each).

Review of the same (Oct 3). The measurement reproduces at 930f066 and
under 'saturation' with 'midpoint' to the count (246,813 onsets at
N=128), and with both off the GPU state and the deck, mixing-top,
regime, rain and flux fields hash the same as the parent's after 16
steps at N=64 and N=128. Over those 16 steps the saturation stage keeps
each layer's θ_l to 3.1·10⁻⁵ K and q_t to 1.9·10⁻⁹ (one f32 ulp, as in
the parent) and the global water budget's residual is unchanged
(2.29·10⁻⁶ against 2.28·10⁻⁶ kg/m² a step at N=64). A three-day run
saved after its second day and continued ends byte for byte on the
file of the run that did not stop. Two days at N=128 from
eleven128_day1825: SWCRE −49.3, −49.4 → −52.4, −52.7 W/m², LWCRE 25.5,
25.7 → 25.8, 26.0, rain 2.68, 2.69 → 2.70, 2.73 mm/d, ASR 241.7, 241.6
→ 238.6, 238.2, OLR 235.1, 235.0 → 234.8, 234.6. A step costs
161.4–164.2 ms before and 161.7–163.7 after (two runs each, exclusive
lock). In cloudEffect.test.mjs 30 of 362 columns are left out against
20 at the parent (5 parted on a cloud decision against 3), so it passes
only with the bound on the excluded share at 0.09 (30/362 = 0.083); the
precision bounds on the columns kept are unchanged and met (SWCRE rms
4.0·10⁻⁵, max 0.07 W/m², against 2.2·10⁻⁵ and 0.02 at the parent). The
deck blinks left are a third cycle of the same kind: at 34.8S 127.9E
the carried height sits at 1465–1467 m across a layer's midpoint, the
water path holds at 54–57 g/m² and the cover goes 1, 0.3, 1, 0.3 (57
blinks in 60 steps); the free troposphere the deck entrains
(θ_l and q_t above) is still the first layer whose midpoint lies above
h, so the jump and with it the decoupling ratio switch with the
midpoint. Of the 1,730 deck-driven transitions at N=128 the cover moves
with the blink in all, the water path in 411, the count of layers below
h changed in 1,408.

**'uniform' the default (Oct 3).** `boundaryCondensation` defaults to
'uniform': every mixed layer holds the uniform distribution's
condensate, the arrangement ECHAM6 and the Unified Model's large-scale
scheme under Lock et al. (2000) publish, and the user's choice, for
realistic schemes with parameters tuned to observations; 'cloudLayer' and
'saturation' stay as options with their code untouched (CPU, 24 steps at
N=6: the parent's default and the new explicit 'cloudLayer' hash the
same, as do 'saturation' on both and the parent's explicit 'uniform'
against the new default; the GPU takes its defaults from MOIST_DEFAULTS).
No pinned digest changed. The five tests are restated on one principle:
a parity test compares the engines where the physics is deterministic
and leaves out the columns where a discrete decision parted, with their
neighbours within two cells, asserting that share against a bound from
the CPU's own sensitivity. `test/helpers/decisions.mjs` reads, after each
step on either engine, the regime, the plume's firing, its base flux (by
1 %) and top, each layer's cloud water against the cloud-top threshold
(10⁻⁶) and the dry adjustment's merges; `scripts/perturbedRain.mjs` and
`scripts/perturbedDecisions.mjs` (SETUP cloudEffect, dayMeans, treeline)
count the same on the CPU against itself with ±10⁻⁴ K of θ noise before
every step. The rain accumulation (gpuModel.test.mjs, 24 steps at N=6):
over 18 seeds a decision parts in 36–49 cells, 191–255 within two cells
(0.53–0.70 of 362), and no cell whose rain parts by 10⁻³ of the largest
lies outside them; the GPU 38 and 191 (all 20 parted cells within one
cell of a decision); bounds C/6 and 0.8 C; over the 171 kept cells
convective rms 2.0·10⁻⁵, large-scale 7.6·10⁻⁵ and the last step's rain
9.5·10⁻⁶ kg/m² against the unchanged 10⁻³, 10⁻³ and 3·10⁻⁴, and no kept
cell's rain apart by 10⁻³ of the largest cell's. The cloud
effects (cloudEffect.test.mjs, its own rules joined to the shared one):
over 12 seeds 12–27 columns parted, 0–5 whose bottom block of q merged
apart, 143–248 left out; the GPU 24, 4 and 199; bounds C/10, C/50 and
0.8 C; over the 163 kept SWCRE rms 9.8·10⁻⁶ (max 5.3·10⁻³ W/m²) against
the unchanged 2·10⁻³ and 0.5. The day means (dayMeans.test.mjs): over 16
seeds 37–58 parted, 181–247 left out, the day-mean ASR over all columns
apart by up to 1.7·10⁻⁴ of itself (above 10⁻⁴ in 5 seeds, as between the
engines) and over the kept by at most 6.6·10⁻⁶; the GPU 50 and 229, its
day means over the 133 kept ASR 344.534 / 344.536 W/m² (5·10⁻⁶ of
itself), OLR to 4·10⁻⁷ of itself, albedo to 4·10⁻⁶; bounds C/5 and 0.8 C, the 10⁻⁴ limits now on
the kept columns' means. The treeline (48 steps with ocean and land):
the atmosphere's decisions part in 160–198 columns and cover every land
cell, so the test leaves out the land's own decisions: what parted is
cell 351, whose step fell in season on the GPU only at steps 36 and 40
(each moving the 6-hour mean by 0.042); at ±10⁻⁴ K, whose lowest-air
drift (up to 3.9·10⁻² K) matches the engines', the season or snow parts
in 0–1 of 151 land cells over 12 seeds (at most 16 left out), at ±10⁻³ K
in 1–5 (12–58); bounds 0.04 and 0.4 of the land cells; over the 133 kept
season length 1.5·10⁻⁷, trees within their targets' gap. Snow covers 31
land cells on both engines (24 under 'cloudLayer'). The stratiform
lifetime (convection.test.mjs) counted 263 of 2,727 cloudy layers kept
on both engines against 272.7, the denominator grown by cloud in mixed
layers whose long share is 0; it now counts the layers its rule gives a
positive share: 257 of 876 (29 %; 196 of 723, 27 %, under 'cloudLayer'),
floor a fifth.

The blinking measured again at main, where the page's cloud column has
counted the shallow cumulus' cover × water since 3da313c (merged with the
branch, not in its measurement): 0.162 % a step at N=128 (15,918 onsets)
and 0.091 % at N=64 (4,619), at the page's cadence 0.143 % and 0.117 %,
2.25 % and 2.75 % of cells blinking at least once, lag-1 autocorrelation
of the step change −0.44 and −0.36 over all cells (−0.53 and −0.52 over
the blinking), 86 % and 88 % of cloudy runs one step long; in 95 % of the
blinks the cumulus' change is larger than any grid-scale part's. Under
'cloudLayer' at main the same: 0.161 % and 0.086 % (15,787 and 4,383),
−0.45 and −0.36, the cumulus larger in 94 % and 93 %; the 0.010 % and
0.006 % above were measured without the cumulus in the frame. Three GPU
days at N=64 from eleven64_day1825 (days 1826–1828, the parent → 'cloudLayer'
→ 'uniform'): albedo 0.318, 0.314, 0.311 → 0.326, 0.323, 0.321 → 0.328,
0.325, 0.324; SWCRE −59.3, −57.7, −56.6 → −61.8, −60.7, −60.0 → −62.7,
−61.4, −60.9 W/m² (over the three days −3.8 W/m² from the parent under
'uniform', −3.0 under 'cloudLayer'); LWCRE
29.2, 29.6, 28.8 → 29.5, 29.9, 29.2 → 29.6, 30.0, 29.3; rain 2.77, 2.77,
2.70 → 2.79, 2.78, 2.73 → 2.81, 2.81, 2.75 mm/d; ASR 232.1, 233.6, 234.6
→ 229.5, 230.6, 231.2 → 228.7, 229.9, 230.3; OLR 229.9, 229.7, 230.6 →
229.6, 229.3, 230.1 → 229.5, 229.1, 229.9 W/m².

**Blinking by cloud type (Oct 3).** `scripts/cloudFlicker.mjs` also
counts the blinks of each type's own frame field through its own
overlay (the page's curve stretched to the type's range): the baseline
at main's physics (f319996), from the same states and steps as above,
N=128 / N=64. The `cloud` row reproduces the combined numbers.

| type (range, g/m²) | cells visible | blinks per step | cells ever | lag-1 of the step change: all cells, blinking cells | one-step runs: all, cloudy |
|---|---|---|---|---|---|
| all clouds (100) | 79.8 / 85.4 % | 0.162 / 0.091 % | 2.25 / 2.75 % | −0.44 / −0.36, −0.53 / −0.52 | 64 / 62 %, 86 / 88 % |
| low (200) | 52.1 / 63.3 % | 5 / 8 onsets | 0.002 / 0.017 % | 0.31 / 0.20 | — |
| mid (500) | 31.5 / 33.0 % | none | none | 0.58 / 0.48 | — |
| high (400) | 23.2 / 25.1 % | 0 / 1 onset | 0 / 0.002 % | 0.44 / 0.37 | — |
| cumulus (40) | 43.8 / 51.6 % | 1.18 / 1.79 % | 13.2 / 30.4 % | −0.48 / −0.49, −0.51 / −0.49 | 58 / 53 %, 69 / 65 % |
| deck (150) | 4.8 / 4.3 % | 0.0085 / 0.0082 % | 0.07 / 0.17 % | −0.33 / −0.21, −0.54 / −0.48 | 64 / 55 %, 61 / 48 % |

By region (blinks per cell and step, N=128 / N=64), the cumulus: sea
1.37 / 1.96 %, land 0.71 / 1.37, tropics 1.63 / 2.56, ITCZ box 2.11 /
2.94, storm tracks 0.66 / 1.08 (north) and 0.98 / 0.94 (south), SE
Pacific 0.88 / 1.34, Peru 0.23 / 0.45, Namibia 0 / 0.03, California 0 /
0.57, N Pacific 0.54 / 0.82; the deck: Peru 0 / 0.059, Namibia 0.105 /
0.085, California 0.167 / 0.027, SE Pacific 0.010 / 0.010, none on land
or in the ITCZ. The resolved cloud does not blink at any height on its
own scale; the cumulus does, 30 % of all cells blinking at least once
in the 128 steps at N=64 (13 % in 64 at N=128), and the combined overlay
sees only the part of it that crosses 0.4 on the 40 g/m² curve.

**The cumulus memory and the deck's reference (Oct 3).** The blinks
left at main were the plumes' cloud (95 %, M21's memory) and a deck cycle
at 34.8S 127.9E (above): the free troposphere the deck entrains (θ_l and
q_t above h, as `thetaLAbove` and `qtAbove`) was the first layer whose
midpoint lies above the carried height, so it switched layers as h
crossed a midpoint and the jump, and through the decoupling ratio the
cover, went 1 ↔ 0.3. `deckReference` 'interpolate' moves
that layer's θ_l and q_t toward the next layer up's by the share of the
height from the midpoint below h to the first one above that lies below
h, on both engines and in the host replica of
`scripts/figures/mlmdeck.mjs`: continuous in h, the layer's own value
with h on the midpoint below, the next layer's as h reaches the midpoint
above, as the slab already weights the layer h lies in; 'layer' (the
default) is the reference before, under which both engines hash as at
f319996 (GPU as for M21's memory; CPU: the pinned deck digests in
physics.test.mjs). 'layer' stays the default because the interpolation
as built is a climate change (below), while the cycle it removes is
0.008 % of cells a step; the reference that belongs to the physics is
the free troposphere's air at h itself, extrapolated down from the
layers above, which the deck's tuning against observed cover, water
path and thickness (minimumInversion, the entrainment efficiency and a
drizzle sink) is to set together with the gate's threshold. The
subsidence's bracket was already the two interfaces about h and is
unchanged (only the density it divides by is the layer's that holds h).
By layer the reference sat from one layer's spacing above h down to h
itself as h rose toward a midpoint, half a spacing on average; the
interpolated reference sits one spacing above h wherever h lies (exactly
so on even spacing), so the jump the gate tests and the entrainment
reads is larger. Over the ice-free sinking sea columns the host replica
tests at the fourth step from day 1825, the virtual jump grows by a
median 1.2 K at N=64 (mean 2.2 K, 10 % of columns by more than 5 K) and
1.6 K at N=128 (mean 2.6 K), and the columns passing the 4 K test go
from 3,608 to 5,776 of 12,655 (N=128: 15,821 to 25,898 of 45,260); the
shift is not uniform, so no single minimumInversion restores the old
pass count (6.5 K still passes 4,276 at N=64). The
host replica agrees with the GPU on the night side as before (eleven64:
964 decks on both, none on one only, LWP ≥ 1 g/m² within 2.4·10⁻⁴
relative against 2.8·10⁻⁴ at f319996; eleven128: 3,916, 2.8·10⁻⁴ against
2.9·10⁻⁴; no gate decision apart).

Measured with `scripts/cloudFlicker.mjs` (BOX 10–20N 160–140W, the same
steps as before) from eleven128_day1825 (64 steps) and eleven64_day1825
(128), f319996 → this branch: blink onsets 15,918 → 18 at N=128 (0.162 %
→ 0.0002 % of cells a step) and 4,619 → 44 at N=64 (0.091 % → 0.0009 %);
at the page's cadence 0.143 % → 0.002 % and 0.117 % → 0.040 %; cells
blinking at least once 2.25 % → 0.01 % and 2.75 % → 0.11 %; the lag-1
autocorrelation of the step change over all cells −0.44 → +0.31 and
−0.36 → +0.26 (over the blinking cells −0.53 → +0.06 and −0.52 → +0.06);
one-step cloudy runs 86 % of 10,167 → none of 3 and 88 % of 2,886 → 2 of
2. Before, the cumulus' change was larger than any grid-scale part's in
95 % and 96 % of the blinks; after, in none, and the few left are
resolved low cloud (56 % and 80 %) and the deck (33 % and 5 %). Apart:
the memory alone leaves 760 and 194 onsets (0.008 % and 0.004 %), 98 %
and 82 % of them the deck's; the reference alone 15,253 and 4,506, the
cumulus larger in 99.9 % and 98.7 %, the deck-driven transitions 1,485
and 298 → 13 and 4. Over the trades box the eight-step strip loses the
pale speckle the cumulus put on and off. Three N=64 days from
eleven64_day1825 (days 1826–1828), f319996 → both, with a replicate
(cumulusMemory 1801 s): SWCRE −62.7, −61.4, −60.9 → −63.9, −64.5, −64.1
(replicate −63.9, −64.4, −64.1) W/m²; LWCRE 29.6, 30.0, 29.3 → 29.7,
30.3, 29.6; albedo 0.328, 0.325, 0.324 → 0.332, 0.334, 0.333; ASR 228.7,
229.9, 230.3 → 227.4, 226.8, 227.1; OLR 229.5, 229.1, 229.9 → 229.4,
228.9, 229.6 (ASR − OLR −0.8, +0.8, +0.4 → −2.0, −2.1, −2.5); rain 2.81, 2.81,
2.75 → 2.81, 2.81, 2.76 mm/d. Two N=128
days (1826–1827): SWCRE −53.3, −53.6 → −54.4, −56.5 (replicate −54.5,
−56.5); LWCRE 25.9, 26.1 → 26.0, 26.4; albedo 0.302, 0.303 → 0.305,
0.311; ASR 237.7, 237.3 → 236.5, 234.5; OLR 234.7, 234.5 → 234.6, 234.2;
rain 2.72, 2.75 → 2.72, 2.76. The memory moves the cloud effects by
0.1–0.2 W/m² (M21); the −2.5 W/m² of SWCRE at N=64 (−2.0 at N=128) is
the reference's, run alone: −63.8, −64.2, −64.0 and −54.3, −56.4. Its
larger jump passes the gate's 4 K test (minimumInversion, set against
the layer reference) more often: the gate stands open over 0.100 of the
globe at the end of day 1828 against 0.060 (0.073 of 30S–30N against
0.040; N=128 day 1827: 0.097 against 0.061), the carried height there
1449 against 1498 m. At the end of day 1828 (`scripts/cloudRegimes.mjs`,
one CPU step, f319996's code for its state) the deck's cover is 0.059
against 0.028 globally, total cover SE Pacific 10–30S 110–80W 0.41 → 0.52
(SWCRE −35.4 → −54.3 W/m²), Peru 0.45 → 0.63 (−35.2 → −64.8), Namibia
0.52 → 0.68 (−66.9 → −115.0), California 0.84 → 0.87, toward the decks'
observed 0.6–0.8 low cover; but the Atlantic trades 10–20N 50–25W go
0.64 → 0.77 against Earth's 0.35–0.55, the deck there 0.16 → 0.39, and
the N Pacific trades stay 0.83. minimumInversion has not been restated
against the interpolated reference.

### M23 — The equatorial ocean — in progress

What the atmosphere's changes will not fix on their own. The M21
drag divisor released an undercurrent that overshoots: in the paired
run nine it surfaced as a 0.5 m/s eastward jet on the equator in June,
collapsed to a 0.35 m/s westward current by September and surfaced
again at the March equinox when the trades fell to zero stress, so the
cold tongue swung between 3.9 K and 0.3 K instead of settling. The
diagnosis of M21 item 4 named the sinks: the ∇⁴ closure pulls the
thermocline classes toward the westward mixed layer at their token
edges (−0.9 to −3·10⁻⁷ m/s² within two rings of a token edge), the
50 m floors on the mixed layer hold the eastern equatorial mixed layer
at 50 m where Earth has 20–40, and the thermocline classes at 0.25
kg/m³ spacing are patchy on the equator. The ocean also hits its 5 m/s
speed limit at 1S 99–100E off Sumatra at N=64 in every recent run.
The shearMixing and closureFill options of M21 exist to try here.
Acceptance: a 0.2 m/s westward surface current east of 140W, an
undercurrent of 0.5–1 m/s near 100 m within two degrees of the
equator, the eastern 1024 class top at 40–60 m, 2 K between the warm
pool and the cold tongue held through a year, and no speed-limit
clamps.

**Diagnosis (Sept 30).** On the nine128 states (`scripts/equatorialOcean.mjs`,
`scripts/closureCoupling.mjs`, `scripts/waveSpeeds.mjs`), 140–120W / 120–100W,
2S–2N:

- Mixed-layer budget, day 91 (10⁻⁷ m/s²): stress −5.42/−6.29, −g∇η
  +8.06/+5.83 (8.2/5.9 cm per 1000 km), baroclinic −0.48/−0.40, drag
  below +0.01/−1.31, closure −0.14/+0.01; u₀ +0.39/+0.69 m/s. Day 365:
  stress −3.60/−4.08, −g∇η +5.26/+6.94; u₀ +0.67/+0.19 m/s. Stress over
  pressure force in the 50–55 m mixed layer 0.67–1.08.
- The atlas start's tilt: 1024 top 172 m west, 59 m east, sea level 45 cm
  higher in the west; by day 91 the east 45–67 m deeper (1024 top 117 m at
  120–100W) and 31 cm. c₁ = 2.5–2.7 m/s, Kelvin crossing 64–70 days,
  Rossby 190–210.
- The classes 1022.25–1024.0 hold water in 12–43 % of the equatorial
  cells, 70–96 % of their edges tokens (0–25 % at the atlas start). The
  closure ties them to the mixed layer at 0.5–3.1×10⁻⁵ s⁻¹ (drag r/h
  4×10⁻⁶). Column closure against stress (10⁻⁵ m²/s²): day 91
  +9.76/+3.61 against −2.95/−3.41, day 183 −7.84/−0.36 against
  −7.13/−5.44, day 274 +4.81/+2.68 against −3.03/−2.96, day 365
  +10.99/+0.51 against −2.09/−2.12. Day 91, 140–120W: 1024.00 at
  +33 cm/s with pressure force −22.05 and closure +21.68×10⁻⁷ m/s²,
  1024.50 at +44 cm/s with −16.31 and +24.26.
- 140–100W surface current +0.54, −0.36, +0.27, +0.42 m/s on days 91,
  183, 274, 365; a 175 m slab spun up by the stress in about 47 days.
- Day 365 stress: −0.021 to −0.034 N/m² over 180–100W, +0.08 and +0.13
  over 160E–180 and 140E–160E.
- Cell 6777 (1S 100E, N=64) started with 1020.5 from 50 to 1194 m under a
  50 m atlas column; 165 such columns within 25° at N=64.

**Changed (Sept 30).** The atlas fill from the nearest column that reaches
the bottom (**Start from the World Ocean Atlas**, M18): 165 light columns
to 13. The interior fill of the closure with its transpose, the default
(**Closure at token edges**, M18): fastest-mode growth per N=128 step
1.12–1.30 to 0.99; coupling to the mixed layer 0.0–47.6 to −0.9–2.5×10⁻⁶
s⁻¹ (19.9 in 1021.5 on day 274). `interiorShearMixing` (**Drag**, M18),
off.

60 coupled N=64 days from `nine64_day0365` (day 365: u(5 m) −0.03/−0.02
m/s, strongest eastward 0.02 m/s at 69 m, 1024 top 128/91 m, W−E 0.5 K)
and `eight64_day0365` (W−E 2.0 K), `{"everySteps":8}`, the unmodified
code and the change; 2S–2N, u(5 m) 160E–100W / 140W–100W, the strongest
eastward class over 180–100W with its depth and the band where it holds
half its peak, τ 160E–100W, h₀ 140–100W, 1024 top 150E–180 / 120–90W,
warm pool / cold tongue / W−E from the spin-up log, clamps as edge-days
to that day:

| Run | Day | u(5 m) (m/s) | strongest eastward | τ (N/m²) | h₀ | 1024 top | SST (°C), W−E | clamps |
|---|---|---|---|---|---|---|---|---|
| nine, before | 15 | −0.07 / −0.03 | 0.04 at 161 m (1024.75), 2S–8N | −0.022 | 51 m | 122 / 88 m | 25.8 / 25.9, 0.3 K | 0 |
| nine, after | 15 | −0.06 / −0.02 | 0.04 at 162 m (1024.75), 3S–8N | −0.019 | 51 | 124 / 88 | 25.8 / 25.9, 0.3 | 0 |
| nine, before | 30 | −0.03 / −0.01 | 0.04 at 68 m (1022.25), 0–2N | −0.024 | 51 | 137 / 84 | 25.6 / 26.0, 0.0 | 0 |
| nine, after | 30 | −0.00 / −0.02 | 0.04 at 159 m (1024.75), 3S–7N | −0.013 | 51 | 132 / 85 | 25.7 / 26.0, 0.2 | 0 |
| nine, before | 45 | +0.05 / −0.01 | 0.23 at 74 m (1022), 3S–1N | −0.007 | 50 | 143 / 83 | 25.7 / 26.1, 0.0 | 0 |
| nine, after | 45 | +0.00 / +0.08 | 0.07 at 132 m (1024.25), 1S–6N | −0.010 | 50 | 131 / 85 | 25.6 / 26.0, −0.0 | 0 |
| nine, before | 60 | +0.16 / +0.00 | 0.34 at 76 m (1022), 2S–2N | −0.005 | 51 | 133 / 83 | 25.7 / 25.9, 0.2 | 0 |
| nine, after | 60 | +0.01 / +0.10 | 0.09 at 161 m (1024.75), 3S–3N | −0.013 | 51 | 132 / 88 | 25.4 / 26.1, −0.2 | 0 |
| eight, before | 15 | +0.00 / +0.11 | 0.08 at 169 m (1024.75), 2S–7N | −0.013 | 50 | 148 / 89 | 26.7 / 25.7, 1.6 | 2 |
| eight, after | 15 | −0.00 / +0.09 | 0.08 at 169 m (1024.75), 3S–7N | −0.015 | 50 | 148 / 88 | 26.7 / 25.7, 1.5 | 2 |
| eight, before | 30 | −0.09 / −0.03 | 0.02 at 207 m (1025.25), 6S–8N | −0.019 | 51 | 148 / 89 | 26.6 / 25.9, 1.0 | 2 |
| eight, after | 30 | −0.02 / +0.04 | 0.04 at 152 m (1024.5), 4S–8N | −0.010 | 50 | 149 / 89 | 26.7 / 26.0, 1.2 | 2 |
| eight, before | 45 | −0.08 / −0.09 | 0.01 at 186 m (1025), 2–8N | −0.017 | 51 | 150 / 83 | 26.6 / 26.0, 0.9 | 8 |
| eight, after | 45 | −0.03 / −0.02 | 0.03 at 170 m (1024.75), 0–8N | −0.010 | 51 | 149 / 90 | 26.7 / 26.1, 1.1 | 2 |
| eight, before | 60 | −0.01 / +0.06 | 0.05 at 95 m (1023), 3S–4N | −0.013 | 51 | 150 / 87 | 26.5 / 26.1, 0.8 | 10 |
| eight, after | 60 | −0.03 / +0.04 | 0.02 at 169 m (1024.75), 1–8N | −0.013 | 51 | 149 / 88 | 26.6 / 26.1, 1.0 | 2 |

The clamps: 13.6S 144.3E on days 379–380 in both eight runs, 1S 130.7E
on days 409–411 before; none at 1S 99–100E in any run. The interior fill
without its transpose, from the same states, day 60: u(5 m) +0.15/+0.20
and +0.02/+0.09 m/s, strongest eastward 0.19 at 87 m and 0.06 at 88 m,
clamps 0 and 8. With `interiorShearMixing` as well: +0.02/+0.10 (nine,
0.21 at 181 m) and −0.09/−0.06 (eight, 0.09 at 186 m), 1024 top east 87
and 77 m, W−E 0.4 and 0.8 K, and 1.5–3.0 m/s in thermocline classes at
0–3N 35–42W, 1S 145E, 7.7N 84.6E and 20N 125E; with backgroundViscosity
10⁻³ m²/s +0.12/+0.10 and +0.11/+0.29, and 1.6–2.3 m/s at 3S–3N
32–41W, 140E and 118E.

5 coupled N=128 days from `nine128_day0091` (day 91: +0.18/+0.54 m/s,
0.47 m/s at 108 m in 1023, W−E 3.9 K), the unmodified code / the change:

| Day | u(5 m) 160E–100W / 140W–100W (m/s) | strongest eastward | 1024 top 150E–180 / 120–90W | W−E | clamps |
|---|---|---|---|---|---|
| 92 | +0.17 / +0.55, +0.16 / +0.53 | 0.46 at 107 m (1023), 0.36 at 107 m (1023) | 157 / 104, 157 / 104 m | 3.9, 3.9 K | 0, 0 |
| 93 | +0.17 / +0.56, +0.15 / +0.51 | 0.45 at 106, 0.31 at 106 | 157 / 104, 157 / 103 | 3.8, 3.8 | 0, 0 |
| 94 | +0.16 / +0.56, +0.13 / +0.49 | 0.43 at 105, 0.28 at 105 | 158 / 103, 158 / 102 | 3.8, 3.8 | 0, 0 |
| 95 | +0.14 / +0.55, +0.11 / +0.45 | 0.42 at 107 (1022.75), 0.25 at 106 (1022.75) | 158 / 103, 158 / 101 | 3.7, 3.7 | 0, 0 |
| 96 | +0.13 / +0.52, +0.09 / +0.40 | 0.41 at 107, 0.23 at 105 | 156 / 105, 156 / 102 | 3.7, 3.7 | 0, 0 |

Day 96, the band holding half the peak: 1S–2N (0.53 m/s at 0°) before,
1S–3N (0.29 at 0°) after; τ −0.024 to −0.026 N/m², h₀ 53–54 m in both.
The interior fill without its transpose clamped 585, 497, 546, 573 and
584 edges on days 92–96.

Against the acceptance, day 60 at N=64 under τ −0.005 to −0.013 N/m²
and day 96 at N=128 under −0.025: u(5 m) 140–100W +0.10, +0.04 and
+0.40 m/s (target −0.2); strongest eastward 0.09, 0.02 and 0.23 m/s at
161, 169 and 105 m (target 0.5–1 near 100 m); eastern 1024 top 88, 88
and 102 m (40–60); W−E 0.5 → −0.2, 2.0 → 1.0 and 3.9 → 3.7 K; clamps 0,
2 (13.6S 144.3E) and 0.

**Review (Sept 30).** The closure with its transpose is a drag on the
undercurrent classes of the states it starts from. On `nine128_day0091`
at 140–120W and 120–100W, 2S–2N, the classes 1022.5–1024.0 at +23 to
+77 cm/s take −11.3 to −50.9×10⁻⁷ m/s² from it, against +0.3 to
+21.7 through the tokens, the pressure force −22.1 to +4.1 and the
interfacial drag −10.9 to +1.2 (`scripts/equatorialOcean.mjs`). As a
damping rate of each class's own flow over 2S–2N 180–100W
(`scripts/closureDrag.mjs`, 1022.25–1024.5, 10⁻⁶ s⁻¹):

| State | tokens | default | first ring alone |
|---|---|---|---|
| nine128 days 91, 183, 274, 365 | −4.4 to +0.2 | 6.6 to 95 | −15.8 to +10.6 |
| eqB128_day0096 (5 days, unmodified code) | −2.0 to +0.4 | 7.3 to 60 | −2.9 to +3.4 |
| eqE128_day0096 (5 days, the change) | −32.7 to −0.9 | 0.2 to 13.3 | −37.8 to +0.5 |
| nine64_day0365 | −0.8 to +0.4 | 1.0 to 7.8 | −0.1 to +1.9 |
| eqEnine_day0425 (60 days, the change) | −5.5 to +0.2 | 0.3 to 4.8 | −2.8 to +0.3 |

In the N=128 pair the strongest eastward class fell from 0.47 to 0.36
m/s on the first day with the change and to 0.46 without it. The 60-day
N=64 pair from `nine64_day0365`, repeated, matched its saved states bit
for bit at days 380, 395, 410 and 425. Over 10S–10N the
grid-scale share of each class's thickness-flux divergence (cell minus
the mean of its thick neighbours over the variance, white noise 1.17;
`scripts/divergenceNoise.mjs`)
after the 5 N=128 days was 0.80, 0.89, 0.41 and 0.28 in the mixed layer,
1022–1024.75, 1025–1026.5 and the deeper classes away from tokens, 0.98
and 1.20 beside the fitted tokens and 1.45 and 1.24 beside the tokens
over the sea floor, against 1.02, 0.84, 0.56, 0.26, 1.07, 1.11, 1.59
and 1.11 with the unmodified code and 1.26, 1.05, 1.05, 0.92, 1.87,
1.56, 1.53 and 1.50 after one day of the fill without its transpose; at
day 60 of the N=64 pair 0.48, 0.77, 0.38, 0.39, 0.75, 0.85, 1.37 and
1.22 with the change, 0.64, 0.68, 0.37, 0.33, 0.94, 0.95, 1.88 and 1.26
without it.

The default stays at `closureFill` 0 (the token edges carry the layer above's velocity as before): the two-ring fill with its transpose damps the thermocline classes' own flow by 0.2–95·10⁻⁶ /s on the audited states and took the N=128 undercurrent from 0.47 to 0.36 m/s in a day, while the first ring alone with its transpose showed no systematic drag but was not run for stability. The fill and its transpose remain options; the first-ring form is the next step here. The atlas deep fill is on.

### M24 — Performance — in progress (seven levers done)

The goal is a model day in a minute at N=128 on the M1 Max with every
scheme of M21–M23 in place. Three levers, each an option whose off value
steps bit for bit as 086dccf did (checked on both engines, and for two
model days at N=128 against 086dccf itself), with its realism cost
measured against a replicate (every θ scaled by 1 ± 5·10⁻⁷) and its time
saving measured on the M1 Max at N=128 (`scripts/profileGpu.mjs`, GPU
timestamps per pass, and `scripts/paceGpu.mjs`, whole model days):

| | 086dccf, page (ocean every 4) | 086dccf, spin-ups (ocean every 8) | now (radiation every 4, ocean every 8) |
|---|---|---|---|
| dynamics: three RK4 stages with their advance | 42.8 ms | 42.9 | 43.3 |
| dynamics: the fourth stage | 13.1 | 13.1 | 13.2 |
| physics (radiation, surface, deck; PBL diagnosis, gravity waves, orography) | 41.1 | 41.1 | 23.7 |
| adjust (moist, mixing, drag, dissipation) | 24.8 | 24.8 | 25.5 |
| ocean, amortised over its calls | ≈ 48 | 25.3 | 26.2 |
| ∇⁴ closures, divergence damping, sponge, combine | 10.7 | 10.8 | 11.1 |
| GPU time per step | 180.5 ms | 158.0 | 142.8 |
| a model day, steps awaited one at a time | 94.2 s | 81.6 | 74.2 |
| a model day, queued 8 to a submission | 94.6 s | 81.6 | 73.7 |

(N=128 from runs/eleven128_day1825.bin, 24 steps profiled before and
32 after; the page's worker loop in Node, `scripts/pageLoop.mjs`, runs
at 19.9 simulated hours a minute, 72.3 s a model day, with the frame's
fields; N=64 goes from 11.8 to 9.5 s a model day.)

Two model days at N=128 from eleven128_day1826 with these defaults,
against 086dccf at the spin-ups' cadence (`scripts/spinup.mjs`, whose
log lines the radiation every step and the ocean every 8 reproduce
exactly): day 1 Ts 14.79 against 14.78 °C, ASR 237.1 alike, OLR
234.4 against 234.5 W/m², SWCRE −53.8 and LWCRE 26.2 alike; day 2 Ts
14.82 and ASR 235.5 alike, OLR 233.9 against 234.0, LWCRE 26.7 against
26.6, precipitation 2.75 alike; no clamped ocean edges, currents at most
1.19 m/s in both. The ocean every 4 steps with these defaults gave the
same day lines.

**The radiation held between full calls** (`radiationEvery`, both
engines; the drivers' default from `js/cadence.module.js`, RADIATION_MINUTES
11.25: every 4 steps at N=128, every 2 at N=64). As in every GCM (CAM and
the IFS call their radiation hourly) the full longwave and shortwave,
clear-sky pass included, run every k steps. Between calls the shortwave
heating of the layers and the surface and top shortwave fluxes scale
with the cosine of the zenith angle over the call's, and a change in the
surface's emission since the call, σT⁴ less the call's, is absorbed by the
layers in the shares the call found and the rest escapes (the approximate
updates of Hogan and Bozzo 2015, J. Adv. Model. Earth Syst. 7, 1401: the
upward fluxes' derivative with respect to the surface's emission, one
upward pass per g-point at the call); the longwave heating, the downward
longwave and the column's own emission stay the call's. The call takes
the insolation-weighted mean cosine of the steps it covers, Σμ²/Σμ of the
positive ones: a cell whose sun rises between calls has sunlight to scale
(the mean over the sunlit part, Hogan and Hirahara 2016, Geophys. Res.
Lett. 43, 482, does that too but left the surface's sunlight 0.3 and 1.0
W/m² further below every step's at 45 and 90 minutes than the weighted
mean, whose first-order error in μ averages out over the interval). The
deck, the cloud its sunlight reads and the surface fluxes of heat and
vapour step every step under the current sun, and the day-mean sums take
what each step received: over 2–3 days at N=64 (every 1, 2 and 4 steps)
and 2 at N=128 (every 1, 4 and 8), the accumulated ASR less OLR equals the heating the
physics pass applied to the layers and the surface to 0.0004 W/m² in the
global mean and 0.002 in any cell (f32 sums). Measured with
`scripts/radiationInterval.mjs` against every step, from
eleven64_day1826 over 3 days (two replicates) and eleven128_day1826 over
2 (W/m², K):

| | ASR | OLR | LWCRE | surface SW | Ts | land diurnal Ts amplitude | rms of per-cell day-mean OLR, day 1 |
|---|---|---|---|---|---|---|---|
| N=64 replicates | −0.03, −0.02 | +0.01, 0.00 | −0.01, 0.00 | −0.04, −0.02 | 0.000 | 0.00, −0.03 % | 0.77, 0.82 |
| N=64 every 2 (11 min) | −0.04 | −0.07 | +0.06 | −0.02 | 0.000 | −0.01 % | 1.15 |
| N=64 every 4 (22.5 min) | −0.10 | −0.15 | +0.14 | −0.05 | +0.002 | −0.02 % | 1.42 |
| N=64 every 8 (45 min) | −0.28 | −0.25 | +0.24 | −0.13 | +0.004 | −0.07 % | 2.00 |
| N=64 every 16 (90 min) | −0.43 | −0.31 | +0.30 | −0.02 | +0.010 | −0.21 % | 3.32 |
| N=128 replicate | −0.01 | 0.00 | 0.00 | −0.01 | 0.000 | +0.02 % | 0.81 |
| N=128 every 4 (11.25 min) | −0.03 | −0.07 | +0.07 | −0.02 | +0.001 | +0.01 % | 1.13 |
| N=128 every 8 (22.5 min) | −0.08 | −0.12 | +0.12 | −0.02 | +0.003 | +0.01 % | 1.42 |

The cost grows with the interval in minutes rather than in steps: the
OLR falls and the longwave cloud effect rises as the cloud the longwave
sees ages (by half the interval on average), and the shortwave cloud
effect strengthens likewise. 11.25 minutes is the longest interval whose
means stay within 0.1 W/m² of every step; the peak of the diurnal cycle
of Ts over land (15.6 h local at N=64, 15.4 at N=128) moves by under
0.03 h at any interval, its amplitude by under 0.25 %. The physics pass
falls from 41.1 to 23.7 ms at N=128 (19.9 every 8), the radiation's
24 ms falling by (k−1)/k.

**The ocean's step** (the ocean's `everySteps`; the drivers' default from
`js/cadence.module.js`: OCEAN_MINUTES 45 up to N=64 and shorter in
proportion to the cell spacing at finer N, so every 8 steps at both N=64
and N=128, the eleven spin-up's cadence; coarser grids keep the ocean's
own every 4 steps). Coupled models exchange with
the ocean every 30–60 minutes (CESM every 30), but the layered ocean at
N=128 does not take a 45-minute step: every 16 steps from
eleven128_day1826 its currents reach the 5 m/s cap on the first day
(305279 clamped edges) and the run is NaN on the fourth; every 12 (33.75
minutes) held for the two days it was run. At N=64 every 16 (90 minutes)
stayed within a replicate's noise over ten days. Against every 8, over
the run (scripts/spinup.mjs, radiation every step; K, W/m², m):

| | Ts | ASR | OLR | sea surface SW | warm pool / cold tongue SST | mixed layer, global / equatorial | SST rms | clamped, largest current |
|---|---|---|---|---|---|---|---|---|
| N=64 replicate, 10 days | 0.00 | 0.0 | +0.1 | +0.06 | +0.001 / −0.005 | −0.24 / +0.08 | 0.028 | 0, 1.60 m/s |
| N=64 every 16 (90 min) | 0.00 | +0.02 | 0.0 | +0.07 | −0.001 / −0.011 | −0.30 / −0.10 | 0.030 | 0, 1.69 |
| N=64 every 4 (22.5 min) | 0.00 | +0.05 | +0.01 | +0.18 | +0.003 / −0.009 | +0.07 / +0.01 | 0.031 | 0, 1.51 |
| N=128 replicate, 5 days | 0.00 | +0.04 | +0.02 | +0.08 | +0.002 / −0.009 | +0.01 / +0.01 | 0.010 | 0, 1.20 |
| N=128 every 4 (11.25 min) | 0.00 | +0.02 | 0.00 | +0.02 | +0.002 / −0.005 | +0.01 / −0.01 | 0.011 | 0, 1.21 |

The page, which took the ocean engine's own default of every 4 steps,
now steps it as the spin-ups did, halving its ocean (48 to 25 ms a step
at N=128).

**The page's steps queued as one batch** (`model.stepBatch` in the
worker's loop, with `wait` false so that two submissions stay in flight;
the frame's capture queued ahead as before; a pause still takes effect at
the end of the frame's steps). On Metal the steps run back to back either
way, the GPU time per step equalling the wall time, so the rate is the
same: N=128 19.9 simulated hours a minute per step and 19.9/19.8 batched,
N=64 155.9/158.4 and 158.2/158.2 (two 120 s runs of
`scripts/pageLoop.mjs` each). It sends one submission a frame in place of
about ten a step, which on Vulkan saved 16 % at N=128 (127 against
151 ms a step, the Verda benchmark); the spin-up's BATCH default stays 1.
The page itself in Chrome on the M1 Max at N=128 from eleven128_day1825
(the browser pane hidden, so the globe drew few frames) ran at 15.9
simulated hours a minute on 086dccf and 20.2 with all three levers (90.6
and 71.3 s a model day); pausing, resuming and changing the overlay
while running behaved as before. With a pane hidden the page reports no
late frames, so whether one frame's steps in one submission delay the
globe's drawing more than ten submissions a step did is not measured;
the pacer's pause now comes once a frame rather than once a step.

What is left between 74 s and a minute a model day at N=128: the
dynamics' four RK4 stages (56 ms, 39 %), the adjust pass (25.5), the
ocean (26) and the physics pass (24, of which the deck's own cloud and
sunlight, the surface fluxes and the boundary-layer diagnosis are now
most). Candidates are the ocean momentum kernel at L=45, the adjust
kernel's second saturation adjustment after the plume, a longer stable
ocean step at N=128 (every 12 held two days), the ∇⁴ closure passes and
the deck's ring passes, and on the page the frame and overlay costs.

**The atmosphere's dynamics kernels one element a thread, and D kept to
what is read** (branch perf2-atmos-core, four commits, each bit-identical:
after 16 GPU steps from eleven64_day1826 and 8 from eleven128_day1826 the
core's S and PH and the ocean's S and OD are word for word those of
0bcd938, and so is every remaining field of D compared by name but the
RK stages' scratch, which after the fourth holds the closures' values at
the end of a step). The profile's split mode (`SPLIT=1` in
`scripts/profileGpu.mjs`) found the stencil kernels latency-bound at one
thread per layer and element, each layer's loads waiting on its
connectivity's. In turn:

1. flux, divergence, pvVertex, pvEdge and dissipationHeat take one edge,
   cell or vertex a thread and loop over the 36 layers with the element's
   connectivity and geometry held in registers; divCurl tests whether it
   differences U or the first Laplacian once, outside its layer loops.
2. The column diagnosis stores the lowest interfaces' Exner function EXL
   only where it is read: in the fourth RK stage (a kernel of its own,
   columnLast) for the physics pass that follows the stages, and in
   adjust. It stores no interface values of θ, q and qc: cellTendency
   forms them from the column's own layers in each layer's iteration.
   Carrying a layer's lower interface into the next layer's upper gave
   the same values but let the Metal compiler contract the vertical flux
   differently (one tendency from eleven64_day1826 differed at rounding
   level in 90 % of the θ tendencies); formed twice, they are bit for bit
   the stored ones. DEX, which nothing read, and LAPB, which no kernel
   referenced, go.
3. The column and cellTendency each sum the layers' divergence from FLUX
   in the divergence kernel's order (cellTendency in a loop of its own:
   inside its other edge sums the compiler reordered it), and the
   divergence dispatch and DIV go.
4. The closures' scratch fields (LAPA, DIVS, CURLS, LAP1) share memory
   with the RK stages' (FLUX, QV, QE, PHI); neither set is read outside
   its own passes.

The device was never quiet while this was measured: the page in the
Claude app's pane took 10–50 % of the GPU, and in eight attempts over an
hour and a quarter no two streaming probes in a row reached 350 GB/s.
The times are therefore each
kernel's least dispatch time over three alternated rounds of
`WARM=8 STEPS=32 SPLIT=1 scripts/profileGpu.mjs` from eleven128_day1826
(µs per dispatch, which the timestamps resolve to 65.5 µs):

| | 0bcd938 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| flux (×4) | 852 | 393 | 393 | 328 | 393 |
| divergence (×4) | 1376 | 393 | 393 | — | — |
| column, stages 1–3 (×3) | 2097 | 1638 | 786 | 983 | 1049 |
| column, stage 4 | 2163 | 1770 | 983 | 1376 | 1376 |
| pvVertex (×4) | 1770 | 459 | 459 | 459 | 459 |
| pvEdge (×4) | 786 | 459 | 459 | 459 | 393 |
| cellTendency (×4) | 3277 | 3146 | 2949 | 2621 | 2556 |
| divCurl (×2, ∇⁴) | 786 | 459 | 459 | 459 | 459 |
| dissipationHeat | 2097 | 524 | 524 | 524 | 524 |
| pblDiagnose | 4719 | 4194 | 2818 | 2818 | 2752 |
| adjust | 24183 | 23855 | 23069 | 23069 | 23134 |
| every kernel's least times its dispatches, ms a step | 155.3 | 135.0 | 128.5 | 126.0 | 126.5 |
| D at N=128, MB | 806.1 | 806.1 | 641.0 | 617.4 | 405.0 |

The least times overstate 0bcd938's 142.6 ms a step alone by 9 %; scaled
to it, the four commits save about 26 ms a step (−19 in the first, −6 in
the second, −2 in the third, nothing measurable in the fourth), the four
RK stages' kernels falling from 58.5 to 37.6 ms of least time. Whole
model days on the shared device (`scripts/paceGpu.mjs`, two rounds
alternated): N=128 164.2 and 162.6 s a model day on 0bcd938 against
125.7 and 125.7 (−23 %), N=64 25.4 and 21.3 against 22.0 and 18.6, with
identical day means. At N=64, where the kernels launch a quarter of the
threads, no kernel became slower by more than one timestamp tick
(cellTendency, and the unchanged momentum, by one; the least times' sum
36.6 to 31.4 ms a step), so every size takes the new forms. D at N=64
falls from 201.5 to 101.3 MB. Bit-identity also shows in
`scripts/figures/mlmdeck.mjs`, whose 7.5 MB of JSON from eleven64_day1826
is byte for byte 0bcd938's, and the suite passes (607 tests) at the
fourth commit. Not measured alone: the projection of about 116 ms a step
and 60 s a model day at N=128.

**The layered ocean's kernels one element a thread**
(`js/gpu/layeredOcean.gpu.js`). The ocean's tendency ran one thread per
(class, element), so at its 45 classes (the mixed layer and 44 density
classes) every class loaded the element's connectivity again and each
gather waited on an index load first; a class cost 2.4 times what an
atmospheric layer does. In two steps:

- oDivCurl, oLapVelocity, oVertexVort, oEdgePV and oFlux run one cell,
  vertex or edge a thread over the classes with the connectivity loaded
  once, and the edge thickness and its sill cap are one kernel. Bit for
  bit: the device buffers after 16 steps from eleven64_day1826 and 8
  from eleven128_day1826 are identical to 0bcd938's (5d979f7).
- oMomentum, oCellTendency and oKineticPhi likewise, oMomentum with its
  neighbours' PV weights and masks in registers and the thick classes of
  its edge as a 45-bit mask. The expressions are the same, but Metal's
  compiler fuses and orders them differently once the operands sit in
  registers: after one ocean call the velocities differ by at most
  2.2·10⁻⁷ of the largest, h, h·T and h·S by 1.2·10⁻⁷ or less (49814b5).

Times a dispatch at N=128 (`scripts/profileGpu.mjs` with SPLIT=1, 96
steps from eleven128_day1826, before and after alternated twice; the
Claude app's GPU process took about half the GPU throughout, so each is
the least of a kernel's 48–96 dispatches, in µs): oMomentum 23986 →
4325, oCellTendency 7078 → 2621, oKineticPhi 4325 → 2621, oDivCurl 5046
and 5112 → 1049, oLapVelocity 2032 → 786 (twice), oVertexVort 2097 →
590, oEdgePV 1770 → 1049, oFlux 918 → 721, the edge thickness 852 + 328
→ 721. Amortised over the ocean's every 8 steps the ocean's kernels
fall from 35.0 to 14.2 ms a step by these least times, the first
step's kernels 7.1 of it and the second's 12.9. Scaled by 0.68–0.81,
each ocean pass's time on the quiet device (55b18b8 at 142.6 ms a
step) over its kernels' least times on the shared one, that is about
15 ms, a projected 127 ms a step and 65 s a model day.
At N=64 every changed kernel is as fast or faster (oMomentum 5964 →
1311 µs, oCellTendency 1704 → 590, oDivCurl 1180 → 262). Whole model
days on the shared GPU (`scripts/paceGpu.mjs`, 0bcd938 and 49814b5
alternated four times) say the same with the sharing's noise: at N=128
150.5–184.8 s a model day before (median 161.5) and 129.8–152.4 after
(130.0), at N=64 21.1–25.5 (22.4) and 18.5–22.8 (20.7).

Against 0bcd938 and its replicates, run means less 0bcd938's
(`scripts/spinup.mjs` with STRATOSPHERE=1 at the drivers' cadence; K,
W/m², mm/d; the SST rms is the end state's over the sea cells, the max
wind the run's largest daily maximum, m/s):

| | Ts | ASR | OLR | LWCRE | SWCRE | precipitation | warm pool / cold tongue SST | SST rms | clamped, largest current, max wind |
|---|---|---|---|---|---|---|---|---|---|
| N=64 replicates, 3 days | 0.000, 0.000 | 0.00, +0.03 | −0.03, 0.00 | 0.00, 0.00 | −0.03, +0.03 | 0.000, +0.003 | 0.0 / 0.0 | 0.0047, 0.0048 | 0, 1.02, 100.0 and 99.9 |
| N=64 one element a thread | 0.000 | +0.03 | +0.03 | 0.00 | +0.03 | +0.003 | 0.0 / 0.0 | 0.0046 | 0, 1.02, 100.1 |
| N=128 replicates, 2 days | 0.000, 0.000 | 0.00, 0.00 | 0.00, 0.00 | 0.00, 0.00 | 0.00, 0.00 | −0.005, −0.005 | 0.0 / 0.0 | 0.0035, 0.0034 | 0, 1.19, 89.0 and 88.9 |
| N=128 one element a thread | 0.000 | 0.00 | 0.00 | 0.00 | 0.00 | 0.000 | 0.0 / 0.0 | 0.0033 | 0, 1.19, 88.9 |

(0bcd938: warm pool 27.6 and cold tongue 24.2 °C at N=64, 28.3 and
23.9 at N=128, no clamped edges, max wind 99.9 and 88.9 m/s.) The
mixed layer at 140–100W (51 m), the equatorial surface current and the
undercurrent are alike in all four, the daily energy balance agrees to
0.1 W/m² as the replicates' does, and the ocean's tests and the suite
pass with no bound changed.

**The drags laid every 11.25 minutes** (`dragEvery`, both engines; the
drivers' default from `js/cadence.module.js`, DRAG_MINUTES 11.25, the
radiation's interval: every 4 steps at N=128, every 2 at N=64, every
step at N=32 and coarser). The non-orographic gravity-wave drag and the
subgrid orography's blocking and gravity waves follow the large-scale
wind, which changes over hours, yet their columns took 3.6 and 2.1 ms a
step at N=128, a quarter of the physics pass. They are laid at a
column's first step (after a model is built or its physics uploaded) and
at the steps whose number is a multiple of k, and applied at every step
from what was laid last: the gravity-wave accelerations as pushes, the
blocking implicitly against the current wind, u ← (u + Δt a)/(1 + Δt β),
so that it still stops a wind that has slowed or turned since. The
orographic waves' limiter, at most what stops a layer's wind along the
stress, takes the interval kΔt. On the GPU the gravityWaves and
orography kernels are dispatched on one step in k; the CPU keeps the
step each column last laid them in shared memory, as the held radiation
does, so that its workers decide alike (test/gpuModel.test.mjs holds the
engines together with the drags laid every 4 steps, and the limiter over
16). At the drivers' cadences the drags are laid at the full
radiation's calls. With these the column kernels load a column's six
edges once and form each layer's cell wind from them (pblDiagnose,
orography and gravityWaves; bit for bit, by the saved device buffers
after 16 steps at N=64 and 8 at N=128).

Against the every-step drags (0bcd938, which dragEvery 1 reproduces bit
for bit on both engines), with the full-precision day means of each day
(`model.diagnostics()` at the day's end, the numbers
`scripts/spinup.mjs` logs, unrounded) and two replicates of it, every θ
scaled by 1 ± 5·10⁻⁷ (W/m², K, mm/d):

| | ASR | OLR | SWCRE | LWCRE | Ts | precipitation |
|---|---|---|---|---|---|---|
| N=64 replicates, 3 days | −0.027, +0.025 | −0.008, +0.001 | −0.027, +0.026 | +0.006, −0.002 | +0.0005, −0.0003 | +0.0003, +0.0011 |
| N=64 every 2 (11.25 min) | +0.003 | +0.009 | +0.003 | −0.010 | +0.0004 | +0.0009 |
| N=128 replicates, 2 days | −0.003, +0.006 | −0.003, +0.004 | −0.003, +0.006 | +0.001, −0.006 | +0.0002, +0.0001 | −0.0007, −0.0001 |
| N=128 every 4 (11.25 min) | −0.013 | +0.004 | −0.013 | −0.007 | 0.0000 | −0.0004 |

Every difference is within 1.5 times the largest spread among the run
and its replicates (the N=128 ASR and SWCRE, 0.013 against 0.009, the
most). The spin-ups' logs of the same runs against the replicates of
0bcd938: no clamped ocean edges in any; the largest current 1.02 m/s at
N=64 and 1.19 at N=128 in all; the day's largest wind, as a mean over
the days, 91.93 against 92.03, 92.00 and 91.97 m/s at N=64 and 87.80
against 87.80, 88.05 and 87.65 at N=128; the top six layers' largest
wind (100 m/s at N=64, 85 at N=128) and horizontal Courant number (0.32,
0.27) those of the replicates. Over 30 days at N=64 from
eleven64_day1826, against 0bcd938: no NaN and no clamped edge in either,
the top six layers' largest wind 100 m/s in both (each layer's at most
the base's) and their Courant number at most 0.32 in both, the largest
wind 103.8 against 106.3 m/s and the largest current 1.78 against 1.79.

The saving, at N=128 with every dispatch timed in its own pass
(`scripts/profileGpu.mjs` SPLIT=1, base and change alternated): the
desktop app's GPU process held the device through every measurement
(probe medians 158–390 GB/s against 390 alone, about 290 ms of GPU time
a step against 142.6), so each kernel's tight cluster of dispatch times
is taken. gravityWaves (3.85 ms a dispatch) and orography (2.18) run on
one step in four, 4.5 ms a step less, 4.3 against the alone profile's
rows; the edges loaded once take pblDiagnose from 5.44 to 3.80 ms and
the drags' kernels to 2.8 and 1.85 ms a dispatch, 2.0 ms a step less.
Together about 6.5 ms a step here and 6.2 against the alone profile (the
three kernels ran 4–6 % slower here than alone); 0bcd938 against the
change in one session, three rounds each, gives pblDiagnose 5.5 to 3.8
ms, gravityWaves 3.87 to 0.70 ms a step and orography 2.16 to 0.46, 6.6
ms a step less. Whole model days under that load (`scripts/paceGpu.mjs`,
151–166 s a day at N=128 for both) do not resolve it.

**The adjust kernel's shallow cumulus fed from the deep plume** (GPU
only, bit for bit). The deep plume (`diagnosedPlume`) falls back to the
shallow cumulus in five places after building the column's temperature,
pressure, layer mass, height and plume environment from a state it has
not yet changed; the shallow cumulus (`cumulusFrom`) now takes those six
arrays instead of building them again, and builds its own only for the
separate shallow plume after the deep plume (with the deck's gate closed
it returns before building any). The deep plume's per-layer rain
arrays are cleared only on its own path. Every device buffer matches
0bcd938 word for word after 16
steps from eleven64_day1826 and 8 from eleven128_day1826. The adjust
kernel at N=128 falls from 23.8 to 21.7 ms a step (−2.1: 2.0 the
hand-over, 0.1 the clearing), at N=64 from 6.6 to 6.0. These are its
times on a shared device: another GPU process took nearly all of the GPU
through the night (the streaming probe at 170–390 GB/s against 390), so
the variants were interleaved step by step in one process (64
split-profiled dispatches each in each of two runs at N=128, 192 at
N=64) and each figure is the tight lower cluster of a variant's dispatch
times, its 10th percentile. Alone the kernel runs about a fifth faster
(19.0 ms), so the saving there is likely nearer 1.7 ms. A model day at
N=128 under the same load: 294.0 ms a step before and 292.7 after (two
pairs of `scripts/paceGpu.mjs`); at N=64 the pairs scattered by ±10 %.

**Columns sorted by convective class** (tried, not kept). A counting
sort of the column order within blocks, so that each 32-lane SIMD group
of the adjust kernel runs columns of one class, costs under 0.07 ms but
gains less than 2 ms. With last step's deep convection as the costliest
class (then the deck's gate open with a positive buoyancy flux, open,
closed), adjust took 1.7–2.9 ms longer than in the natural order with
blocks of 128, 256 and 1024 columns. At N=64 26 % of columns convect
deeply and 92 % of those did so the step before, but the other 8 % and
the columns the test parcel types deep without a deep closure (the full
plume runs for 39 %) land in most groups: in blocks of 256, 61 % of the
sorted groups hold a deep column, against 85 % unsorted and 32 % were
the class known. With the test parcel's typing as the class (99 %
persistent) and two classes, adjust took 1.7–1.8 ms less (blocks of
1024), and 1.2–1.4 ms less on top of the hand-over (blocks of 1024 and
4096); a random order within blocks of 1024 took 15 ms more.

**The four levers together** (branch perf-night: 55b18b8 with
perf2-atmos-core, perf2-ocean-kernels, perf2-slow-drags and
perf2-adjust-order merged in that order, the largest saving first). The
code merged without a conflict; pblDiagnose takes both the atmosphere's
column diagnosis that stores only what is read and the drags' six edges
loaded once. The suite passes (71 files, 609 tests); the only test
changes are the drags' two new tests and their dragEvery key in two
exact assertions, no bound changed. On the quiet device (the streaming
probe at 356–372 GB/s throughout), 55b18b8 and perf-night alternated,
from eleven128_day1826 and eleven64_day1826 (`scripts/profileGpu.mjs`
WARM=8 STEPS=32, three profiles of 55b18b8 and four of perf-night;
`scripts/paceGpu.mjs`, four model days of each at N=128 awaited and one
queued, four two-day runs at N=64; each column's runs within 0.4 ms
and 0.2 s of each other):

| | 55b18b8 | perf-night | less |
|---|---|---|---|
| dynamics: three RK4 stages with their advance | 43.1 ms | 27.8 | 15.3 |
| dynamics: the fourth stage | 13.2 | 8.5 | 4.7 |
| physics (radiation, surface, deck; PBL diagnosis, gravity waves, orography) | 23.4 | 15.7 | 7.7 |
| adjust (moist, mixing, drag, dissipation) | 24.8 | 19.7 | 5.1 |
| ocean, amortised over its calls | 25.2 | 11.9 | 13.3 |
| ∇⁴ closures, divergence damping, sponge, combine | 10.7 | 9.3 | 1.4 |
| GPU time per step | 140.5 ms | 92.9 | 47.6 (34 %) |
| step median, awaited one at a time | 110.1 ms | 74.9 | 35.2 |
| a model day at N=128, steps awaited one at a time | 72.6 s | 48.4 | 24.2 (33 %) |
| a model day at N=128, queued 8 to a submission | 72.6 s | 48.3 | 24.3 |
| a model day at N=64, steps awaited one at a time | 9.2 s | 7.0 | 2.2 (24 %) |
| GPU buffers at N=128 | 3668 MiB | 3286 | 382 |
| GPU buffers at N=64 | 917 MiB | 822 | 96 |

The four entries' own estimates sum to about 48–49 ms a step (the
atmosphere's kernels about 26, the ocean's 15, the drags 6.2, the
hand-over 1.7–2.0); together they save 47.6, the sum within the
estimates' error. By kernel (two SPLIT=1 profiles each of 0bcd938 and
perf-night, ms a step): the column −10.1 (+1.2 for columnLast),
oMomentum −8.2, cellTendency −4.4, adjust −4.1 (its own column
diagnosis stores fewer fields too, besides the hand-over), pblDiagnose
−4.0, divergence −3.5, oDivCurl −2.3, pvVertex −2.3, gravityWaves −2.2,
divCurl −1.4, orography −1.3, oCellTendency −1.2, dissipationHeat −1.0.
The GPU buffers are those alive after two steps, counted by wrapping
the device's createBuffer; all of the 382 MiB at N=128 is D (768.8 to
386.3 MiB).

Not every lever is bit for bit (the ocean's kernels round differently,
and the drags are laid every 4 steps at N=128 and 2 at N=64), so the day
lines are compared with 0bcd938's, whose model is 55b18b8's, and its two
replicates (`scripts/spinup.mjs` with STRATOSPHERE=1). Over two days at
N=128 from eleven128_day1826 every figure of the day lines is 55b18b8's
but the day's largest wind (89.0 and 86.1 m/s against 88.9 and 86.7, the
replicates 89.0 and 87.1, 88.9 and 86.4), the sea surface's sunlight
(163.2 and 160.3 W/m² against 163.1 and 160.4, the replicates 163.2 and
160.3, 163.2 and 160.4), that of the iced cells north of 60° on the
first day (12.7 against 12.6, both replicates 12.7) and the sea
surface's net longwave on the second (−47.9 against −48.0): Ts 14.79 and
14.82 °C, ASR 237.1 and 235.5, OLR 234.4 and 233.9, SWCRE −53.8 and
−55.4, LWCRE 26.2 and 26.7 W/m², precipitation 2.76 and 2.75 mm/d, no
clamped edge, currents at most 1.19 m/s, the balance lines alike. In the
end lines the regions' rain differs by at most 0.1 mm/d, as the
replicates' does, the Pacific ITCZ's by 0.01 (4.83 against 4.82), the
Kalahari's surface reads 37 against 36 °C (36 in both replicates, the
log's whole degrees), and the equatorial Pacific's lines are alike. Over
three days at N=64 the run means less 55b18b8's are ASR +0.03, OLR 0.00,
SWCRE +0.03, LWCRE −0.03 (the third day's 28.8 against 28.9 at the log's
rounding) and precipitation +0.003, the replicates' sizes; no clamped
edge, the largest current 1.02 m/s alike, and the top six layers'
largest wind and Courant number the replicates'.

The page from perf-night in the hidden Browser pane, with another page
holding the GPU as in the first two timing rounds (the probe at 167–390
GB/s, 55b18b8 at 289–297 ms of GPU time a step): at N=64 14–15 s a
simulated day; the overlays (wind, temperature and rain at N=64, OLR and
low cloud at N=128) switched with their legends, a pause held the clock
and a resume continued it, and the console stayed empty. At N=128 it
read 81.8–82.8 and 82.7–87.1 s a simulated day in two runs of two to
three minutes, 55b18b8's page between them 118.9–120.4 (−31 %); the page
on the quiet device was not measured.

What is left at 92.9 ms a step (SPLIT=1, ms a step): the momentum kernel
15.1 (one edge a thread over the layers on both), adjust 15.6, the
physics kernel 13.2 (the radiation 8.7 of it, amortised), cellTendency
9.4, advance and combine 5.8. Of the 3286 MiB at N=128 the RK4 registers
S, T and K1–K4 hold 1834 (the atmosphere's 6 × 136.9, the ocean's 6 ×
168.8), the ocean's OD 642, D 386 and PH 283.

**The load writes only data, in 1 MiB pieces** (memory; bit for bit).
Once the default state had loaded at N=64, the page's tab held about
1.9 GB in its renderer and 2.3 GB in Chrome's GPU process against the
822 MiB of the model's GPU buffers, and Safari on the iPhone 17 Pro,
which ends a tab at around 1–1.5 GB, crashed on it. Most of the
difference was shared memory mapped into both processes alike, 1258 MB.
Chrome gives a queue.writeBuffer larger than a few MiB, and a buffer
mapped at creation, shared memory of its own in both processes and
keeps it (vmmap shows one region per such write, sized to it), while
writes of 1–2 MiB go through its 16 MiB transfer ring, which waits for
the GPU process when it is full: in a test page 128 MiB written in 1
MiB pieces without a wait left 16 MB of shared memory, one 63 MiB write
64 MB, and writes of 16, 32 and 64 MiB awaited in turn 113 MB. The
worker's load wrote 1355 MiB at N=64 (5421 at N=128), most of it zeros
and copies: the core's and the ocean's K1–K4, D and OD as zeros, T as a
copy of S, the whole of PH three times and the ocean twice, first from
the analytic start that model.load laid and the saved ocean replaced.
Now:

- `core.upload` writes S, and one encoder copies it into T and clears
  K1–K4 and D; the ocean's `uploadArrays` does the same for its S, T,
  K1–K4 and the whole of OD before its sections' writes; `clearFrame`
  clears FR on the device, and `uploadPhysics` clears PH and writes only
  the 1 MiB blocks of it that hold a nonzero word. Each clear is
  submitted before the queue writes it must not wipe, and a load throws
  inside a batch, whose encoder would be submitted after them.
- `model.load({ ocean: false })` leaves the ocean to the ocean's load or
  initialize, which the worker's start, restore and device probe always
  follow it with; `model.load()` alone still initializes it.
- `storageBuffer` (the mesh's MI and MF, LV, the parameter buffers)
  writes through the queue instead of mapping at creation, and it, S
  and the ocean's S go in 1 MiB pieces (`writeInPieces` in
  `js/gpu/device.module.js`), as PH's blocks do.
- The worker waits for the queue after building the model, after
  model.load, after the ocean's load and after placing the land.

The worker's load now writes 147 MiB at N=64 and 537 at N=128. In
headless Chrome for Testing from the page's default parts, with the
phys_footprint 15 s after ready, its peak and the shared VM_ALLOCATE in
each process (MB; ready in seconds from navigation):

| | renderer | peak | GPU process | peak | shared | ready |
|---|---|---|---|---|---|---|
| N=64, 3610d2a (4 runs) | 1811–2000 | 2835–2998 | 2262–2422 | 2414–2546 | 1258–1259 | 3.2–5.5 |
| N=64, the device's clears | 746 | 1307 | 1384 | 1401 | 129 | 3.4 |
| N=64, and no throwaway ocean | 744 | 1111 | 1161 | 1178 | 129 | 2.8 |
| N=64, and the 1 MiB pieces (3 runs) | 637–644 | 991–1005 | 1044–1059 | 1067–1084 | 23–24 | 2.8–3.0 |
| N=128, 3610d2a | 7214 | 9892 | 9136 | 9822 | 5518 | 12.6 |
| N=128, all but the pieces | 2250 | 3496 | 4062 | 4097 | 458 | 8.4 |
| N=128, all | 1783 | 3109 | 3653 | 3683 | 39–40 | 8.2 |

With the clears alone the shared memory was the regions of the writes
of MI (13.6 MiB), MF (16.3), S (34.2) and the ocean's S (42.2) beside
the ring; at N=128 those four were 425 MiB. Of what is left at N=64,
the GPU process holds the model's 865 MB of Metal buffers, and the
renderer 484 MB of the worker's JavaScript arrays (the second model
`sourceFor` builds, the GPU model's CPU core and mirrors) and 117 MB of
other heap. The device after the load is word for word 3610d2a's: the
sha256 of every one of the 23 GPU buffers after the worker's start
(Node, model.worker.js, both default states), after 10 steps at N=64
and 5 at N=128 through model.step and through model.stepBatch with the
worker's snapshot of them, and after a snapshot's restore into the live
model; `scripts/spinup.mjs` from eleven64_day1826 and eleven128_day1826,
BATCH 1 and 8, saves states that cmp equal. `test/gpuReload.test.mjs`
loads a stepped model's state into a fresh model at N=16, steps it and
loads the state again, and holds every buffer but the parameters and
the ocean's frame fields to the fresh load word for word. The step is
untouched: `scripts/paceGpu.mjs` alone on the device, 3610d2a and the
change alternated, gives 48.5 and 48.4 s a model day at N=128 against
48.3 and 48.4, and 7.0 s at N=64 in all four runs, with the same day
means. In the headless Chrome the page loads and runs at both defaults
without a console message (6.7–8.1 s a simulated day at N=64, 47.9 at
N=128), and a snapshot saved at N=32 restores into the running model,
which steps on from it. There a save or a 'Download and restore' at
N=64, as at 3610d2a, and the download at N=128 fail in IndexedDB, which
aborts the write without an error (the page's only console messages).
Not measured: Safari, whose WebKit may hold uploads differently.

### M25 — The long spin-up — planned

The asynchronous schedule of M18 (`scripts/asyncSpinup.sh`: a hundred
ocean-only years, then ten coupled, repeated toward a thousand years)
on a rented GPU, once the coupled model holds a year within a few W/m²
of balance with the physics of M21–M23. The World Ocean Atlas start
makes the first century a drift from a measured state rather than from
an analytic one. The Verda spot prices and the measured H100 and A100
paces of Sept 28 put the schedule near $300–400; the Verda tooling
(`scripts/verdaRelaunch.sh`, `scripts/verdaInstances.mjs`) is in place.
The page's default states are then taken from the spun-up ocean.
## 7. Module layout in this repo

```
js/
  grid.module.js            existing: ISEA cells, neighbors, Voronoi vertices
  isea.module.js            existing: projection
  mesh.module.js            NEW  M0: edges, vertices, kites, TRiSK weights → typed arrays
  dynamics/
    operators.module.js     NEW  M0: div, grad, curl, KE, uPerp, ∇², reconstruct
    shallowWater.module.js  NEW  M1: single-layer test core
    sigmaCore.module.js     NEW  M2: K-layer hydrostatic core (steps 1–8 on arrays)
    sponge.module.js        M21: the top sponge on the zonally asymmetric wind
  geography.module.js       M16: land mask, land fraction, elevation and coast from a raster; M17: smoothed surface geopotential; M22: subgrid orography fields and the data/subgrid_N<N>.bin codec
  physics/
    land.module.js          M16: bucket, snow, land albedo and wetness
    radiation.module.js     ported from sim.js RadiationColumn
    longwave.module.js      M21: the longwave gases' g-points (table in longwaveTable.module.js, written by scripts/longwaveFit.mjs)
    shortwaveGases.module.js M21: ozone, water vapour, O2 and CO2 absorption of sunlight after CLIRAD-SW
    surface.module.js       ported: surface and top drag, ocean wind stress, convective adjustment
    exchange.module.js      M22: the surface layer's C_D and C_H by roughness and stability
    gravityWaves.module.js  M21: Alexander & Dunkerton (1999) non-orographic gravity-wave drag
    boundaryLayer.module.js M14: K-profile boundary layer, implicit column mixing of θ, q, qc and u; M22: the implicit surface and form drag
    orography.module.js     M22: Lott and Miller's blocking and gravity-wave drag (js/gpu/orography.gpu.js its port)
    formDrag.module.js      M22: the turbulent orographic form drag's coefficient
    init.module.js          ported: thermal init, balance, seed, geostrophic winds
    regrid.module.js        barycentric interpolation of a state between meshes; ice, snow and soil by source tile
    moist.module.js         M7/M8: saturation adjustment, cloud water, autoconversion, Betts–Miller, filler; M21: its triggered entraining parcel and shallow branch, the shallow cumulus mass flux and the convective plume
    ice.module.js           M9/M11: zero-layer sea ice over the mixed layer, its concentration, zenith albedo
  ocean/
    layered.module.js       M18: 45-layer hybrid isopycnal ocean with a split free surface, the mixed layer coupled through the sea-ice cell update
    seawater.module.js      M18: the Roquet et al. (2015) simplified equation of state, shared with the GPU
    climatology.module.js   M18: the World Ocean Atlas file a fresh ocean starts from, and its columns
    (js/gpu/layeredOcean.gpu.js is its WebGPU port)
  gpu/
    device.module.js         M15: WebGPU device (Dawn in Node, navigator.gpu in the page) and buffer helpers
    core.gpu.js              M15: layouts, dynamics kernels, RK4 and closures, full-step orchestration
    physics.gpu.js           M15: column physics and adjustment kernels
    forcing.gpu.js           M18: the forcing recorder and the ocean and sea ice stepped alone under it
    model.gpu.js             M15: the GPU model behind the CPU model's interface
    profile.module.js        the model dialog's GPU profile: step times and kernel timestamps
  model.module.js           assembles core + physics, RK4 step, diagnostics
  cadence.module.js         M24: the radiation's, the drags' and the ocean's intervals the drivers turn into steps
  audit.module.js           M21: the audit's boxes and the spin-up's convection and equator lines
  forcing.module.js         M18: one recorded day of the ocean's surface forcing, encoded and decoded
  oceanHandOff.module.js    M18: a coupled state with the ocean, sea ice and sea surface of another
  parallel.module.js        M6: the same model stepped on worker threads
  parallel.worker.js        M6: one worker's block of every phase
  threads.module.js         M6: worker_threads / Web Worker primitives behind the engine
  model.worker.js           browser worker: steps the model and serves the page's subscription
  pace.module.js            the pacer: an idle pause between GPU steps only while it cuts the page's late frames
  deviceChoice.module.js    the engine and resolution a first visit's device test picks
  gestures.module.js        two-finger pan, pinch and twist from a pair of touch points
  frames.module.js          the page–worker protocol and the catalogue of fields a frame can carry
  climate.module.js         live model page (climate.html): control panel, overlays, wind layers
  levels.module.js          fields on a pressure surface, comfort measures; season phrase for the model time
  stats.module.js           the optional frame-rate panel (mrdoob's stats.js)
  windParticles.module.js   wind shown by particles advected by the field, their opacity rising with speed
  charts.module.js          synoptic charts from saved states (charts.html)
  unifiedViewer.module.js   existing, gains model-overlay mode; colours per-cell values on the GPU
scripts/
  spinup.mjs, spinup.sh     one spin-up segment on the GPU from the newest snapshot, and a loop of them
  pairedSpinup.sh           resolutions spun up in step, one at a time, compared every EVERY days
  oceanSpinup.mjs           the ocean and sea ice spun up alone under a recorded year of forcing
  asyncSpinup.sh            coupled and ocean-only phases in turn at one resolution, supervised
  runControl.mjs            the spin-ups' stop on SIGTERM and their SYNC_CMD queue
  verdaRelaunch.sh          a Verda spot instance recreated on its OS volume after each eviction
  verdaInstances.mjs        the verda CLI's JSON as verdaRelaunch.sh reads it
  compareStates.mjs         saved states side by side as a markdown table
  verticalAudit.mjs         M21: the vertical-motion and convection audit's headline numbers from one state, on the CPU
  radiationBenchmark.mjs    M21: the gas radiation of one CPU column on the standard atmospheres against its references
  longwaveFit.mjs           M21: the longwave spectral model's fit and its reduction to the g-points
  standardAtmospheres.mjs   M21: the benchmark's profiles laid onto a model column
  upperAtmosphere.mjs       M21: the layers above 200 hPa of saved states, the top treatment's force, the spin-up's wind line
  remapState.mjs            a saved state carried onto another sigma grid
  packWoa.py                data/woa_annual_1deg.bin from the World Ocean Atlas 2023 NetCDF files
  subgridTerrain.py         M22: GMTED2010's 30″ tiles assembled and filtered to 2′30″ grids in the download cache
  subgridTerrain.mjs        M22: data/subgrid_N<N>.bin from those grids per mesh
  subgridTerrainHand.py, subgridTerrainCells.mjs  M22: a few cells recomputed by direct sums from the 30″ grid
  splitState.mjs            a saved state gzipped into parts for the page
test/
  mesh.test.mjs, operators.test.mjs, trisk.test.mjs, sw_tc2.mjs, sw_tc6.mjs,
  sw_galewsky.mjs, rest_state.mjs, held_suarez.mjs, baseline_compare.mjs
```

Performance envelope at N=32, K=20: per step, roughly 6 flops/cell for
divergence, 1/edge for gradient, 3/vertex for curl, 10/edge for TRiSK,
6/cell for KE — order 50 flops per edge per layer, ~3×10⁷ per step, on
flat typed arrays. Comparable to or faster than the current object-graph
core (which does more work per cell through the adjoint gather lists).

---

## 8. Risks and open questions

1. **TRiSK weight convention.** The single most error-prone formula.
   Mitigation: implement against RTSK eq. 22–24 and the MPAS mesh spec,
   certify with T1–T3 before any dynamics, and treat TC2 as the final
   oracle. Do not reason about signs from memory (Section 4.2 note).
2. **ISEA cell asymmetry (measured).** TRiSK's reconstruction is exact
   for a uniform field only where the dual edge bisects the primal edge.
   On the raw ISEA grid the crossing point is up to 0.19 `l_e` off the
   primal midpoint on ~30% of edges at every N, and the tangential
   velocity — hence the Coriolis force — is 12% wrong on the worst of
   them (T2, Section 3.5). Conservation (T1, T3) is unaffected. Lloyd
   relaxation (centers → spherical-polygon centroids, circumcenters
   recomputed, topology unchanged) is the standard remedy — MPAS meshes
   are centroidal Voronoi tessellations for this reason — and measured
   at N=32 it takes the worst edge from 12% to 1.6% after 4 iterations
   and 0.9% after 8, the mean from 1.5% to 0.2%, with the cell-area ratio
   unchanged (1.245 → 1.244). The residual sits at the pentagons. `Grid`
   now relaxes by default (`relax: 10`); the raw grid remains available
   with `relax: 0` and both are covered by the test suites.
3. **Extra velocity degrees of freedom.** E = 3C edges carry more
   information than a 2-component vector field needs, so hexagonal
   C-grids have a computational mode branch in the divergent part.
   TRiSK keeps it from growing; the ∇⁴ and, if needed, a weak divergence
   damping keep it quiet. Watch for it in TC6 as small-scale divergence
   noise.
4. **Interaction with AB4.** The C-grid has no null mode but its
   gravity-wave phase speeds at 2Δx are *higher* than the smoothed A-grid
   operators' (the A-grid's low gain at grid scale was accidentally
   softening the CFL). Expect the same `dt` at the same spacing but less
   margin; RK3 is the fallback.
5. **Coriolis placement.** `f` at vertices for the vorticity, at edges for
   the geostrophic init. Poles are ordinary pentagon cells — no
   singularity, no arbitrary basis — but `f` at the pole *vertex ring* is
   the maximum; check TC2 error near the poles specifically.
6. **Base-state physics is not fixed by this work.** The Eady analysis
   showed the free troposphere is ~1.65× too statically stable (gray
   radiation equilibrium), capping eddy growth at ~3× slower than Earth
   even with perfect numerics. That is a radiation-scheme item for after
   M5.
7. **Two large changes remain stacked** (new mesh + new core) by decision.
   The milestone structure — M0/M1 certify the mesh and horizontal core
   against exact solutions before any column physics is attached — is the
   mitigation: if M1 passes TC2/TC6/Galewsky, the mesh and operators are
   right, and any later regression is in M2/M3.
8. **The layered ocean's eddies are parameterized, not resolved.** Its
   cells are 112 km at N=64 and 56 km at N=128, while the ocean's eddies
   are 10–30 km across (eddy-permitting begins near 25 km), so the
   eddy transport that flattens isopycnals across the Antarctic
   Circumpolar Current and sets the Southern Ocean overturning — the path
   by which warm deep water reaches the surface under the pack ice —
   comes from the Gent–McWilliams transport of M18 (**Eddy transport**),
   with along-layer (Redi) mixing implicit in the coordinate. Sixty days
   show the interfaces smoothing but cannot show the overturning: judge
   it over the next multi-century spin-up by the Southern Ocean column
   (60–70S, 200–500 m, Earth +1 °C) and the Antarctic pack volume. κ,
   nearly uniform poleward of 20°, is the parameter to revisit, for
   example with the Eady growth-rate dependence of Visbeck et al. (1997).

---

## 9. Decisions recorded

- C-grid built directly on the ISEA `Grid`; no intermediate port of the
  A-grid model to the new mesh.
- The integrated model lives in this repository.
- Second-order centered transport and AB4 first; higher-order transport
  and RK3 are upgrade paths, not prerequisites.
- The A-grid model in `~/Desktop/climate_model` is retired; it is not a
  validation baseline.
