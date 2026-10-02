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
  `sigmaInterfaces('bl34')` is the same grid above its 2.4 km
  interface (σ 0.7444) with the five layers below it replaced by
  twelve: interfaces at 40, 83, 133, 192, 265, 359, 482,
  649, 882, 1213 and 1698 m in the standard atmosphere, each layer up
  to 1.5 times as thick as the one below, 34 layers in all. A saved
  state records its interfaces (`levels`), and every loader builds the
  model on them; a state without them is on cam26. `remapLevels`
  (`js/physics/regrid.module.js`) carries a state onto another grid
  conservatively in σ, each layer the σ-weighted mean of the layers it
  overlaps, so cam26 and bl34 exchange their layers above 2.4 km
  unchanged; a spin-up seeded from a state on the other grid carries
  its atmosphere across with it. The page starts a run from nothing
  (`climate.html?from=none`) on bl34; a run started from a saved state
  is on that state's grid.
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
instance alive from an always-on machine with the `verda` CLI (1.8)
logged in. Every POLL (120) s it lists the instances and looks at the
one called NAME: a running one is left alone and its OS volume
remembered in STATE_FILE (`$HOME/.verda-relaunch-NAME`); one starting
or stopping is waited for; an offline one is started; and when there is
none, or it was discontinued, it creates a spot instance again with
`verda --agent vm create --kind gpu --instance-type INSTANCE_TYPE
--location LOCATION --is-spot --os <OS volume or image>
--os-volume-size OS_VOLUME_SIZE --os-volume-on-spot-discontinue
keep_detached --ssh-key SSH_KEY --hostname NAME --startup-script
STARTUP_SCRIPT --wait -o json`, on OS_VOLUME when it is set and
otherwise on the remembered OS volume, once `verda volume list` shows it
detached, and on the image OS only while no volume is known. A remembered volume that is no longer
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
   `rayleighDepth` 0, both aerosol depths 0 and `skylight` 0.15 both
   engines are the previous ones bit for bit. Over a black surface the
   clear column reflects τ/(τ + 2μ) of the visible beam less ozone and
   aerosol absorption to 1e-14 at three zenith angles, every column of a
   real-geography state closes to 5e-15, and on a random set of sunlit
   columns the engines differ by 1.1e-5 of the beam
   (`test/clearScattering.test.mjs`). The daily line gives the clear-sky
   reflectance and the sea's surface sunlight; the audit and the second
   sweep's `clearAlbedo` term (target 0.15, tolerance 0.01, weight 3) read
   the clear-sky albedo.
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

### M24 — Performance — planned

The M21 physics costs about 2.5 % a step for the plume and more for
the 45-class ocean (its momentum kernel's thick-layer search grows
faster than the layer count), and the machine measured 57–78 s a model
day at N=128 on the evening of Sept 30 against 54 before. The goal is
a model day in a minute at N=128 with all of the above in place.
Profile with `js/gpu/profile.module.js` back to back; the candidates
are the ocean momentum kernel at L=45, the adjust kernel (16 of an
86 ms step), the second saturation adjustment after the plume, the ∇⁴
closure passes and the deck's ring passes, and on the page the frame
and overlay costs.

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
  geography.module.js       M16: land mask, land fraction, elevation and coast from a raster; M17: smoothed surface geopotential
  physics/
    land.module.js          M16: bucket, snow, land albedo and wetness
    radiation.module.js     ported from sim.js RadiationColumn
    surface.module.js       ported: surface and top drag, ocean wind stress, convective adjustment
    boundaryLayer.module.js M14: K-profile boundary layer, implicit column mixing of θ, q, qc and u
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
  packWoa.py                data/woa_annual_1deg.bin from the World Ocean Atlas 2023 NetCDF files
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
