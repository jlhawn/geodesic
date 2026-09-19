# C-Grid Dynamical Core on the ISEA Icosahedral Mesh

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
drag, sensible heat flux), by the Rayleigh sponge, and by the viewer.
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
- `F_e`: Rayleigh drag terms (PBL `σ>0.7`, top `σ<0.05`), bulk surface drag
  on the lowest layer using `|v|_e = ½(|v_i| + |v_j|)` from 3.7, and the
  ∇⁴ closure `−K₄ (∇⁴u)_e`.

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
bounds it, which is why removing the drag was not an option.

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
`climate.html?snapshot=`. The trough metric is the zonal-mean surface
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
  70%. When the implied rain is positive the reference temperature is
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
  follows the model's own humidity: 2 m²/kg times each layer's water
  mass (`vaporCoupling`; 0 restores the prescribed Frierson profile,
  which the single-column initial profile still uses). This is the
  water-vapour feedback the prescribed profile lacked.

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
decays over 3 h; both are exact exponential decays per step, and the
rain falls out at once. Betts–Miller still rains directly. Dry
adjustment mixes `qc` with `q`. The Betts–Miller reference profile uses
Bolton's closed-form LCL and a two-substep moist ascent that stops once
the parcel is 10 K colder than the air (87 → 5 ms per N=16 step).
Column cloud water (TCW) is a diagnostic and a page overlay; the step
costs 126 ms serial at N=16 against 98 ms dry. Cloud–radiation
coupling (albedo and longwave emissivity of cloud) is the next step.

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
  physics/
    radiation.module.js     ported from sim.js RadiationColumn
    surface.module.js       ported: drag, sensible heat, slab ocean, convective adjustment
    init.module.js          ported: thermal init, balance, seed, bands, geostrophic winds
    regrid.module.js        barycentric interpolation of a state between meshes
    moist.module.js         M7/M8: saturation adjustment, cloud water, autoconversion, Betts–Miller, filler
  model.module.js           assembles core + physics, RK4 step, diagnostics
  parallel.module.js        M6: the same model stepped on worker threads
  parallel.worker.js        M6: one worker's block of every phase
  threads.module.js         M6: worker_threads / Web Worker primitives behind the engine
  model.worker.js           browser worker: steps the model, fills shared buffers
  climate.module.js         live model page (climate.html): control panel, overlays, wind layers
  levels.module.js          fields on a pressure surface; season phrase for the model time
  windParticles.module.js   wind traced by particles with fading trails over the globe
  charts.module.js          synoptic charts from saved states (charts.html)
  unifiedViewer.module.js   existing, gains model-overlay mode
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

---

## 9. Decisions recorded

- C-grid built directly on the ISEA `Grid`; no intermediate port of the
  A-grid model to the new mesh.
- The integrated model lives in this repository.
- Second-order centered transport and AB4 first; higher-order transport
  and RK3 are upgrade paths, not prerequisites.
- The A-grid model in `~/Desktop/climate_model` is retired; it is not a
  validation baseline.
