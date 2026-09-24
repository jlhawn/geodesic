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
decays over 3 h; both are exact exponential decays per step, and the
rain falls out at once. Betts–Miller still rains directly. Dry
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
The initial state carries 0.5 m of ice wherever the initial surface is
below freezing (poleward of ~60°). A prescribed ocean heat transport
(`oceanHeatFlux`, Q₀ = 20 W/m²) converges Q₀(3 sin²φ − 1) into the
mixed layer — zero in the global mean, cooling the tropics, warming
the poles by up to 2 Q₀ — and melts ice from below where it is
covered. It stands in for the poleward heat carried by ocean currents,
without which the ice–albedo feedback runs the ice edge to ~45°.

### M10 — Clouds in radiation — done

Each layer's cloud water path gives it a gray emissivity
1 − exp(−130 m²/kg × path) that joins every longwave band: the vapour
and gas bands as 1 − (1 − ε_gas)(1 − ε_cloud), and the window, which is
now an exchange band of its own that is transparent only where there is
no cloud. In the shortwave the column's cloud optical depth
(`cloudScattering` = 60 m²/kg × path) reflects the beam
with the two-stream reflectance τ/(τ + 2μ); what passes is absorbed by
the surface with its per-cell albedo, with the multiple reflections
between surface and cloud summed. The fixed planetary albedo of 0.3 is
gone: it is now produced by clouds and ice, and diagnosed. Radiation's
closure still holds exactly (with the latent heat of evaporation
counted as leaving the surface). The cloud–albedo, cloud–longwave and
ice–albedo feedbacks are all live from here.

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
  (`test/ice.test.mjs`). It runs on the main thread in both engines
  (`phases.ocean`), which keeps the workers' cell updates free of
  neighbour reads and the parallel step bit-identical to the serial
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
atmosphere steps on the main thread, 7–11 % of the 253 ms step.

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
cold climate); the drag and exchange coefficients are 3e-3 over land against
1.5e-3 over water, as per-cell arrays that the surface drag, the
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
own temperature and salinity over seven interior layers of fixed
reference density 1022.0, 1023.0, 1024.0, 1025.0, 1026.0, 1026.6 and
1026.95 kg/m³ (`LAYER_DENSITIES`; 26.3, 23.1, 19.4, 15.2, 10.0, 5.9 and
2.65 °C at 35 psu), on the real bathymetry `D` (the cell-mean ETOPO depth,
at least 50 m, zero on land). Every layer is a TRiSK shallow-water layer
in the vector-invariant form carrying thickness, edge velocity, heat h·T
and salt h·S. Density is the simplified equation of state of Roquet et
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
(ζ̄ + f̄)/max(h_e, 20 m) with the same edge thickness, which keeps the PV
term bounded where a layer thins to nothing. A layer thinner than 5 m at
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
their transport equals the averaged U. Velocities are finally clamped
to 5 m/s and the count of clamps reported (`oceanLimited`, zero in every
run so far).

**Mixed layer.** After each step, per column: the mixed layer swallows
any interior layer lighter than itself (convection), taking only as much
as keeps it within 200 m and at most 100 m a day — unlimited swallowing
with the neutral return below mixed 150 m of 272.4 K deep water into the
polar mixed layer every 22-minute step, a heat supply no winter cooling
could beat, and the sea ice was gone by day 90; entrains the first layer below at the
Kraus–Turner wind-stirring rate w = 2 m u*³/(h₀ Δb) with
m = 0.8 exp(−h₀/100 m), so the wind's stirring fades below the depth it
can reach, and Δb the buoyancy step to that layer, at least 10⁻³ m/s²
(a 0.1 kg/m³ step; the KT rate is unbounded as Δb → 0);
detrains any depth beyond 200 m and, when it is at least as dense as
the water beneath it (within 0.005 kg/m³, convectively neutral),
everything below the 50 m floor (also the minimum thickness) — the swallow and the return each step
still mix the column, but a kilometre-deep mixed layer beside 40 m ones
opened a 5.7 m sea-level hole within an hour off Cape Farewell, and
500–900 m ones that had exhausted the 1027.2 layer and sat, stably, on
the 1027.7 layer ran the Drake Passage at 5 m/s;
detrained water goes to the interior layer whose density is nearest
its own, so water swallowed from a layer returns to that layer (the
first rule, the first layer at least as dense, put 900 m of 1027.24
water into the 1027.7 layer and made a 0.5 m sea-level step with
2.5 m/s currents around it; Bleck's mass-conserving split between the
two bracketing layers had no step at the bottom but ratcheted a tenth
of a metre of every swallow-and-return cycle into the denser class,
and in 90 days 2400 m of the 1027.2 layer had become 1027.7 water in
some Southern Ocean columns and not their neighbours, with 5 m/s
mid-depth flows between); water lighter than the first interior layer
goes there whole (the tropical and summer mixed layers, a known slow
drift);
when the surface buoyancy flux implied by its temperature change since
the last ocean step is stabilising (below −10⁻⁹ m²/s³) and it is deeper
than the Monin–Obukhov depth 2 m u*³/(−B), it detrains the excess over a
day into the first interior layer at least as dense as itself that
already holds water in the column, or the deepest such layer, or not at
all on a shelf that has only mixed layer; and it is kept at least 50 m
thick by entraining from below. The shallowest depth detrainment leaves
is 50 m: with 20 m the tropical mixed layer sat on that floor (every
warming step collapses the Monin–Obukhov depth) and the SST fell 2 K a
month.

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
each interior layer starts at its class salinity (35 psu through the
1025.0 class, then 34.9, 34.85 and 34.8, `LAYER_SALINITIES`) and the
temperature that gives its label density there (the deepest class at
about 0.8 °C), blending poleward of 50° over 20° of latitude toward
0.5 °C at the salinity that keeps the density (the 1026.6 class at
34.33 psu, the deepest at 34.78), as polar oceans hold cold, fresh water
on the density surfaces of the warm subtropical thermocline; so its
density is its label; interior layer bases at 90, 170, 300, 500, 700 and
1100 m in the subtropics, scaled by 0.7 + 0.6 cos²φ toward the poles;
every layer lighter than the local surface water outcropped; the deepest
layer filling to the bottom. Polar surface water at the freezing point
and 34 psu (1026.47) floats on the deepest class (1026.95); brine that
raises it past about 34.6 psu sinks into it, as bottom water forms. With
35 psu in every class and at the poles, a first spin-up year lost all its
sea ice: polar mixed layers denser than the deepest class convected
2.4 °C water up all winter. A
saved ocean with a different number of layers loads as this climatology.

**Interfaces.** `advance(surfaceT, ice, oceanFlux, stress, dt)` as in M13
plus `accumulate(evaporation, rain, dt, runoff)` each atmosphere step;
`fields()` gives the frame its mixed-layer depth, SST, SSS, surface
velocity, thermocline depth (the boundary between the 23.1 and 19.4 °C
classes, close to the 20 °C isotherm) and η;
`serialize()` is {h, u, T, S, eta} flattened layer-major, and
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

**Drag.** Interfacial drag is the linear stress r Δu with r = 2×10⁻⁴
m/s; bottom drag is quadratic, C_D |u| u with C_D = 3×10⁻³, applied to
the deepest layer with water at the edge. A linear bottom drag of
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

**Density-consistent interior.** After the mixed-layer exchanges, an
interior layer more than 0.01 kg/m³ from its label mixes in water from
the nearest layer lying clearly (by more than 0.01) on the other side of
the label, a fraction dt/2 days of the full correction per step. The move
is the same conservative transfer of mass, heat and salt as detrainment;
the correction is gradual and dead-banded because the curvature of the
equation of state makes a mixture denser than the linear estimate, and
correcting in full every step overshot and pulled water back the other
way. Two cases have no donor and are left alone: the shallowest water
under the mixed layer when too dense, and the deepest class when too
light. The second matters: deep water formed like NADW (1026.88) and
AABW (1026.93) is lighter than the 1026.95 label, so the deepest class
can hold water up to about 0.17 kg/m³ light (about 2 K warm), bounded by
the water detrained into it; referencing every class to the surface also
merges NADW- and AABW-like water into that one class. The GPU ocean
carries its free surface from the barotropic solve rather than
re-summing the layers each step: in single precision the sum rounded
about 4×10⁻⁶ m low in the same way every step, a steady loss of volume.

**Open.** The freezing point ignores salinity; the deepest class has no restoring when light; the Kraus–Turner constants
and the 50 m minimum depth are first guesses; the barotropic mode has
no explicit filter beyond the sub-step average.

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
  geography.module.js       M16: land mask, land fraction, elevation and coast from a raster; M17: smoothed surface geopotential
  physics/
    land.module.js          M16: bucket, snow, land albedo and wetness
    radiation.module.js     ported from sim.js RadiationColumn
    surface.module.js       ported: surface and top drag, ocean wind stress, convective adjustment
    boundaryLayer.module.js M14: K-profile boundary layer, implicit column mixing of θ, q, qc and u
    init.module.js          ported: thermal init, balance, seed, geostrophic winds
    regrid.module.js        barycentric interpolation of a state between meshes; ice, snow and soil by source tile
    moist.module.js         M7/M8: saturation adjustment, cloud water, autoconversion, Betts–Miller, filler
    ice.module.js           M9/M11: zero-layer sea ice over the mixed layer, zenith albedo
  ocean/
    layered.module.js       M18: eight-layer hybrid isopycnal ocean with a split free surface, the mixed layer coupled through the sea-ice cell update
    seawater.module.js      M18: the Roquet et al. (2015) simplified equation of state, shared with the GPU
    (js/gpu/layeredOcean.gpu.js is its WebGPU port)
  gpu/
    device.module.js         M15: WebGPU device (Dawn in Node, navigator.gpu in the page) and buffer helpers
    core.gpu.js              M15: layouts, dynamics kernels, RK4 and closures, full-step orchestration
    physics.gpu.js           M15: column physics and adjustment kernels
    model.gpu.js             M15: the GPU model behind the CPU model's interface
  model.module.js           assembles core + physics, RK4 step, diagnostics
  parallel.module.js        M6: the same model stepped on worker threads
  parallel.worker.js        M6: one worker's block of every phase
  threads.module.js         M6: worker_threads / Web Worker primitives behind the engine
  model.worker.js           browser worker: steps the model, fills shared buffers
  climate.module.js         live model page (climate.html): control panel, overlays, wind layers
  levels.module.js          fields on a pressure surface; season phrase for the model time
  windParticles.module.js   wind shown by particles advected by the field, their opacity rising with speed
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
