# Crusst Architecture Reference

> Deep technical reference for the Crusst SDF geometry kernel.
> Covers the computation pipeline, meshing internals, blend mathematics,
> and the DAG representation that ties everything together.

---

## Table of Contents

1. [System Overview](#1-system-overview)
2. [Module Map](#2-module-map)
3. [The SDF Expression DAG](#3-the-sdf-expression-dag)
4. [Three Computation Modes](#4-three-computation-modes)
5. [Primitive SDF Functions](#5-primitive-sdf-functions)
6. [CSG Operations](#6-csg-operations)
7. [Blend Profile Mathematics](#7-blend-profile-mathematics)
8. [Multi-Face Fillet/Chamfer System](#8-multi-face-filletchamfer-system)
9. [Adaptive Octree Construction](#9-adaptive-octree-construction)
10. [Dual Contouring Pipeline](#10-dual-contouring-pipeline)
11. [QEF Vertex Placement](#11-qef-vertex-placement)
12. [Face Generation and Winding](#12-face-generation-and-winding)
13. [STEP Export Tiers](#13-step-export-tiers)
14. [Transport and Parametric Paths](#14-transport-and-parametric-paths)
15. [Voxelization](#15-voxelization)
16. [Feature Introspection System](#16-feature-introspection-system)
17. [Interval Arithmetic](#17-interval-arithmetic)
18. [Core Data Types](#18-core-data-types)
19. [Architectural Dependency Graph](#19-architectural-dependency-graph)

---

## 1. System Overview

Crusst is a constructive solid geometry (CSG) kernel built on **signed distance
fields** (SDFs). Every geometric shape is represented as a function
`f: R^3 -> R` where negative values are inside, zero is the surface, and
positive values are outside.

The kernel is structured in three layers:

```
                  +-----------------------------------+
                  |         Builder API                |  builder.rs
                  |   Shape::sphere().union(...)       |
                  +-----------------+-----------------+
                                    | Arc<SdfNode>
                  +-----------------v-----------------+
                  |      Expression DAG                |  dag.rs
                  |  point / interval / gradient       |
                  +-----------------+-----------------+
                                    |
          +-------------------------+-------------------------+
          v                         v                         v
   +--------------+   +-------------------+    +--------------+
   |   Meshing    |   |   STEP Export     |    | Voxelization |
   |  Octree+DC   |   |  Exact or Tess   |    |  Rayon grid  |
   +--------------+   +-------------------+    +--------------+
```

All geometry flows through the DAG. The DAG supports three computation modes
(point, interval, gradient) that are consumed by different downstream stages.

---

## 2. Module Map

| Module | File | Role |
|---|---|---|
| `dag` | `src/dag.rs` | `SdfNode` / `SdfNode2d` enum -- the core intermediate representation |
| `builder` | `src/builder.rs` | `Shape` wrapper with fluent API over `Arc<SdfNode>` |
| `primitives` | `src/primitives.rs` | Pure SDF computation functions for all primitive shapes |
| `csg` | `src/csg.rs` | Scalar CSG combinators (sharp, smooth, round, chamfer) |
| `blend` | `src/blend.rs` | Blend profile curves and Newton-iteration solvers |
| `feature` | `src/feature.rs` | `FaceInfo`, `EdgeInfo`, `FeatureTarget` for targeted blends |
| `shape` | `src/shape.rs` | Trait-based `Sdf` / `Sdf2d` layer with generic wrappers |
| `types` | `src/types.rs` | `BBox3`, `Interval`, `MeshSettings`, `TriangleMesh` |
| `octree` | `src/octree.rs` | Adaptive octree with interval-arithmetic pruning |
| `dual_contouring` | `src/dual_contouring.rs` | Full DC pipeline: octree -> crossings -> QEF -> triangles |
| `qef` | `src/qef.rs` | Least-squares QEF solver via SVD |
| `mesh` | `src/mesh.rs` | Compatibility shim: `extract_mesh()` delegates to DC |
| `export` | `src/export.rs` | Binary STL writer |
| `obj_export` | `src/obj_export.rs` | Indexed Wavefront OBJ writer |
| `step_export` | `src/step_export.rs` | Tiered AP203 STEP exporter (exact BRep or tessellated) |
| `voxel` | `src/voxel.rs` | `VoxelGrid` struct with Rayon-parallel construction |
| `path` | `src/path.rs` | `Path` trait + `LinePath`, `HelixPath`, `SpiralPath` |
| `frame` | `src/frame.rs` | Frenet frame computation for transport paths |
| `transport` | `src/transport.rs` | order0-order3 swept-shape generation |

---

## 3. The SDF Expression DAG

`SdfNode` is the central data structure -- an `enum` wrapped in `Arc` for
cheap cloning and structural sharing. Every CSG tree, transform chain, and
blend operation composes into a single `SdfNode` tree.

### 3.1 Node Variants

**Primitives** (10 variants):
`Sphere`, `Box3`, `Cylinder`, `CappedCone`, `Torus`, `RoundedBox`,
`Capsule`, `Ellipsoid`, `RoundedCylinder`, `HalfSpace`

**Sharp CSG** (3):
`Union(Arc<SdfNode>, Arc<SdfNode>)`,
`Intersection(...)`,
`Difference(...)`

**Smooth CSG** (3):
`SmoothUnion(..., k)`, `SmoothIntersection(..., k)`,
`SmoothDifference(..., k)` -- polynomial blend with parameter `k`

**Round CSG** (3):
`RoundUnion(..., r)`, `RoundIntersection(..., r)`,
`RoundDifference(..., r)` -- inscribed-circle fillet with radius `r`

**Chamfer CSG** (3):
`ChamferUnion(..., k)`, `ChamferIntersection(..., k)`,
`ChamferDifference(..., k)` -- 45-degree chamfer of size `k`

**Targeted Blends** (2):
`Fillet { inner, profile, targets }`,
`Chamfer { inner, profile, targets }` -- per-edge blend with any `BlendProfile`

**Transforms** (6):
`Translate`, `Rotate`, `Scale` (uniform), `Mirror`, `Shell`, `Round`

**2D-to-3D** (2):
`Revolve(Arc<SdfNode2d>)`, `Extrude(Arc<SdfNode2d>, half_height)`

**Opaque**:
`Custom(Arc<dyn Sdf>)` -- escape hatch for user-defined SDFs

### 3.2 Thread Safety

Both `SdfNode` and `SdfNode2d` are statically asserted `Send + Sync`,
enabling safe parallel octree traversal with Rayon.

### 3.3 2D Sub-DAG

`SdfNode2d` is a parallel enum for 2D cross-sections used by `Revolve` and
`Extrude`:

```
Circle2d { center, radius }
Rect2d   { center, half_extents }
Union2d(Arc<SdfNode2d>, Arc<SdfNode2d>)
Difference2d(Arc<SdfNode2d>, Arc<SdfNode2d>)
Custom2d(Arc<dyn Sdf2d>)
```

---

## 4. Three Computation Modes

Every `SdfNode` supports three functions. Each serves a different stage of
the pipeline.

### 4.1 Point Computation: `compute(point) -> f64`

Returns the exact signed distance at a single point. This is the fundamental
operation -- meshing samples millions of these.

Dispatch is a `match` over all variants:
- **Primitives** delegate to `primitives::sdf_*` functions
- **CSG** applies the scalar combinator from `csg.rs`
- **Transforms** transform the point into local space, run inner, then
  scale/transform the result back
- **Fillet/Chamfer** runs the multi-face Lp blend (Section 8)

### 4.2 Interval Computation: `interval_compute(bbox) -> Interval`

Returns a **conservative bound** on all possible SDF values within a bounding
box. If the returned interval is entirely positive, the entire box is outside
the surface. If entirely negative, entirely inside.

This powers octree pruning (Section 9): cells proved empty are never
subdivided, giving O(surface area) memory instead of O(volume).

Implementation per variant:
- **Sphere**: Per-axis intervals -> L2 norm interval minus radius
- **Box3**: Per-axis `|p-c| - half` intervals, composed for outside/inside
- **HalfSpace**: Linear dot-product interval arithmetic
- **Complex primitives** (Cylinder, Torus, etc.): Corner-sampling heuristic
  with slight interval widening
- **Union**: `min(a_iv, b_iv)` using `Interval::min`
- **SmoothUnion**: Sharp union interval widened by `k/4` on the low end
- **RoundUnion**: Sharp union interval widened by `r` on the low end
- **Fillet/Chamfer**: Inner interval widened by `profile.max_deviation()`
- **Transforms**: Transform the bbox, recurse on inner
- **Revolve, Extrude, Custom**: Return `Interval::entire()` (conservative
  fallback -- no pruning possible)

### 4.3 Gradient Computation: `gradient(point) -> Vector3<f64>`

Returns the unit-length outward surface normal direction. Used by the QEF
solver to place mesh vertices (Section 11).

Per variant:
- **Sphere**: `(point - center).normalize()` -- exact analytical gradient
- **HalfSpace**: Returns `normal` directly
- **Box3, Cylinder, Torus, CappedCone, etc.**: Central differences with
  `eps = 1e-6`
- **Union**: `csg_gradient_min` -- gradient of whichever operand has the
  smaller SDF value, with linear blending in a narrow band (`eps = 0.05`)
  to avoid discontinuous normals at the seam
- **Intersection**: `csg_gradient_max` -- gradient of the larger operand
- **Difference**: Blends gradient of `a` and negated gradient of `b`
- **All smooth/round/chamfer CSG**: Central differences
- **Fillet/Chamfer**: Central differences
- **Translate**: `inner.gradient(point - offset)`
- **Rotate**: `rotation * inner.gradient(rotation.inverse() * point)`
- **Scale**: `inner.gradient(point / factor)` -- normalized
- **Mirror**: Reflects gradient through the mirror plane when the point is
  on the reflected side
- **Shell**: `inner.gradient(point)`, negated when `inner.compute(point) < 0`
- **Round**: `inner.gradient(point)`

The CSG gradient blending uses `eps = 0.05` -- a relatively wide transition
band that produces smooth normals at Union/Intersection seams without visible
faceting in the final mesh.

---

## 5. Primitive SDF Functions

All primitives live in `src/primitives.rs` as pure functions
`fn(point: Vector3<f64>, ...) -> f64`. Negative inside, zero on surface,
positive outside.

| Function | Formula / Notes |
|---|---|
| `sdf_sphere` | `\|p - c\| - r` |
| `sdf_box` | Exact: outside = L2 norm of max(q, 0); inside = max component of q (where q = \|p-c\| - half) |
| `sdf_cylinder` | Arbitrary-axis cylinder. Projects to local frame, exact corner distance |
| `sdf_capped_cone` | Inigo Quilez exact formulation. Handles caps, tip, conical wall |
| `sdf_torus` | `sqrt((sqrt(x^2 + z^2) - R)^2 + y^2) - r` |
| `sdf_rounded_box` | `sdf_box - radius` (Minkowski offset of box) |
| `sdf_capsule` | Sphere-swept segment: project point to [0,1] parameter, `\|p - closest\| - r` |
| `sdf_ellipsoid` | Approximate (Quilez bound-corrected): `k0 * (k0 - 1) / k1` |
| `sdf_rounded_cylinder` | Exact Quilez: corner rounding in 2D (radial, axial) space |

**Design note**: `sdf_ellipsoid` is the only approximate primitive. True
ellipsoid distance requires solving a quartic, so the Quilez approximation
is used with its bound-correction factor.

---

## 6. CSG Operations

All CSG operations are scalar functions on two distance values, defined in
`src/csg.rs`:

### 6.1 Sharp CSG

```
union(d1, d2)        = min(d1, d2)
intersection(d1, d2) = max(d1, d2)
difference(d1, d2)   = max(d1, -d2)
```

### 6.2 Smooth CSG (Polynomial Blend)

```
smooth_union(d1, d2, k):
    h = clamp(0.5 + 0.5*(d2-d1)/k, 0, 1)
    result = d2*(1-h) + d1*h - k*h*(1-h)

smooth_intersection(d1, d2, k):
    h = clamp(0.5 - 0.5*(d2-d1)/k, 0, 1)
    result = d2*(1-h) + d1*h + k*h*(1-h)
```

The `k` parameter controls the blend region width. The blending term
`k*h*(1-h)` is maximal (equal to `k/4`) at the midpoint `h = 0.5` where
`d1 = d2`.

### 6.3 Round CSG (Inscribed Circle)

```
round_union(d1, d2, r):
    u = max(r-d1, 0)
    v = max(r-d2, 0)
    result = max(r, min(d1,d2)) - sqrt(u^2 + v^2)
```

This produces an exact circular arc fillet of radius `r` at the seam.
The formula computes the distance to a circle inscribed in the corner of
the `(d1, d2)` distance field.

### 6.4 Chamfer CSG

```
chamfer_union(d1, d2, k):
    result = min(d1, d2, (d1 + d2 - k) / sqrt(2))
```

The third term is the signed distance to the 45-degree plane
`d1 + d2 = k` in `(d1, d2)` space, creating a flat chamfer.

---

## 7. Blend Profile Mathematics

The `BlendProfile` enum in `src/blend.rs` defines ten distinct edge-rounding
curves. All operate in normalized `(d1, d2)` distance space where `d1` and
`d2` are signed distances to the two faces meeting at an edge.

### 7.1 Profile Overview

| Profile | Lp Exponent | Continuity | Shape | Parameterized By |
|---|---|---|---|---|
| G1 | -- (Newton) | G1 (tangent) | Cubic Bezier curve | `radius` |
| G2 | p = 2 | G2 (curvature) | Circular arc | `radius` |
| G3 | p = 4 | G3 | Squircle (smootherstep) | `radius` |
| Chord | p = 2 | G2 | Circular arc (chord-specified) | `chord_length` |
| Cycloidal | p = 1.5 | G1 | Between chamfer and circle | `radius` |
| Parabolic | p = 3 | G2 | Between circle and squircle | `radius` |
| Hyperbolic | p = 6 | G2 | Between squircle and box | `radius, asymptote` |
| EqualChamfer | -- (exact) | G0 | 45-degree flat | `distance` |
| TwoDistChamfer | -- (exact) | G0 | Asymmetric flat | `distance1, distance2` |
| AngleChamfer | -- (exact) | G0 | Angle-specified flat | `distance, angle_rad` |

### 7.2 Exact Closed-Form Profiles

**G2 (Circular Arc)**:
```
u = max(d1 + r, 0)
v = max(d2 + r, 0)
result = max(max(d1, d2), sqrt(u^2 + v^2) - r)
```
This is the distance to a circle of radius `r` centered at `(-r, -r)` in
`(d1, d2)` space. The `max(d1, d2)` outer bound ensures the blend only
activates near the edge.

**EqualChamfer**:
```
result = max(max(d1, d2), (d1 + d2 + k) / sqrt(2))
```
Distance to the line from `(-k, 0)` to `(0, -k)`.

**TwoDistChamfer** (asymmetric):
```
result = max(max(d1, d2), (d1*k2 + d2*k1 + k1*k2) / sqrt(k1^2 + k2^2))
```
Distance to the line from `(-k1, 0)` to `(0, -k2)`.

**AngleChamfer**: Converts `(distance, angle_rad)` to
`(k1 = distance, k2 = distance * tan(angle))` and delegates to TwoDistChamfer.

### 7.3 Newton-Iteration Profiles

The remaining profiles (G1, G3, Parabolic, Cycloidal, Hyperbolic) use a
parametric curve `P(t) = (d1(t), d2(t))` for `t in [0, 1]`. The distance
from a query point `(d1, d2)` to the nearest point on this curve is found
via Newton's method:

1. **Initial guess**: 16-sample uniform scan over `[0, 1]` via `best_initial_t()`
2. **Newton iteration**: 5 steps minimizing `(P(t) - q) . P'(t) = 0`
3. **Sign determination**: Cross product of the curve tangent and the
   point-to-curve vector determines inside/outside

**G1 (Cubic Bezier)**:
Control points: `P0 = (-r, 0)`, `P1 = (-r + r*kappa, 0)`,
`P2 = (0, -r + r*kappa)`, `P3 = (0, -r)` with `kappa = 0.5523`
(standard circular approximation constant).

**G3 (Quintic Smootherstep)**:
```
f(t) = 6t^5 - 15t^4 + 10t^3
P(t) = (-r + r*f(t), -r + r*f(1-t))
```
Derivatives computed analytically through the chain rule.

**Parabolic**:
```
P(t) = (-r + r*t, -r*t^2)
```

**Cycloidal**:
```
theta = u * pi
d1(u) = -r + r * (theta - sin(theta)) / pi
d2(u) = -r * (1 - cos(theta)) / 2
```
Special handling at the cusps (`u ~ 0`, `u ~ 1`) where the tangent
degenerates -- a nudged fallback tangent is used.

**Hyperbolic (Superellipse)**:
```
p = clamp(asymptote / r, 0.1, 0.99)
d1(t) = -r * (1 - t^p)^(1/p)
d2(t) = -r * t
```
Numerical second derivative via finite difference for the Newton step.

### 7.4 Profile Spectrum

The profiles form a continuous spectrum from sharp to boxy:

```
Chamfer <-- Cycloidal <-- G2 (Circle) <-- Parabolic <-- G3 (Squircle) <-- Hyperbolic --> Box
  L1          L1.5            L2              L3             L4               L6
```

Lower Lp exponents produce sharper, more concave blends. Higher exponents
produce flatter, more convex blends approaching a box corner.

---

## 8. Multi-Face Fillet/Chamfer System

The `Fillet` and `Chamfer` DAG nodes use a **multi-face Lp superellipse**
approach that generalizes blending to any number of faces meeting at a point.
This replaces an earlier per-edge approach that failed at corners where three
or more faces meet.

### 8.1 Process Steps

Given a point `p` and `Fillet { inner, profile, targets }`:

1. **Run inner**: `d_sharp = inner.compute(p)`
2. **Retrieve edge info**: `inner.edge_info()` returns all edges with their
   two face indices
3. **Collect targeted faces**: For each targeted edge, add both `face_a` and
   `face_b` to a deduplicated face set
4. **Compute face distances**: For each face `i` in the set, compute
   `d_i = inner.face_distance(p, i)`
5. **Select Lp exponent and clamping**:
   - G2/Chord: `p = 2`, hard clamp
   - Cycloidal: `p = 1.5`, softplus clamp
   - Parabolic: `p = 3`, softplus clamp
   - G3: `p = 4`, softplus clamp
   - Hyperbolic: `p = 6`, softplus clamp
6. **Compute per-face contributions**:
   - Hard clamp: `u_i = max(d_i + r, 0)`
   - Softplus clamp: `u_i = ln(1 + exp(k * (d_i + r))) / k` with `k = 3/r`
7. **Lp norm**: `result = (sum_i u_i^p)^(1/p) - r`
8. **Final**: `max(d_sharp, lp_result)` -- the blend never makes the shape
   smaller than the sharp version

### 8.2 Why Lp Norms?

The standard G2 circular fillet is:

```
sqrt(max(d1+r, 0)^2 + max(d2+r, 0)^2) - r
```

This is exactly the **L2 norm** of the clamped, shifted face distances minus
the radius. Generalizing from L2 to Lp:

```
(max(d1+r, 0)^p + max(d2+r, 0)^p)^(1/p) - r
```

produces different cross-section shapes:
- `p = 1`: chamfer (L1 is the taxicab norm, giving a flat diagonal)
- `p = 1.5`: cycloidal (between chamfer and circle)
- `p = 2`: circular arc (the standard fillet)
- `p = 3`: parabolic (flatter than circle)
- `p = 4`: squircle (G3 continuous)
- `p = 6`: hyperbolic (nearly rectangular)
- `p -> infinity`: box (no rounding)

### 8.3 Softplus Clamping

For higher-order profiles (G3+), the hard `max(d+r, 0)` clamp produces
only G0 continuity at the blend boundary. The **softplus** function:

```
softplus_k(x) = ln(1 + exp(k*x)) / k
```

with steepness `k = 3/r` smoothly transitions from 0 to `x`, providing
G3+ continuity where the blend meets the flat face. As `k -> infinity`,
softplus converges to `max(x, 0)`.

### 8.4 Corner Handling

At a box corner where three faces meet, the Lp norm naturally extends:

```
(u_x^p + u_y^p + u_z^p)^(1/p) - r
```

For `p = 2` this is a sphere cap. For `p = 4` it is a superellipsoid.
This is the key advantage over the earlier per-edge approach, which had
no principled way to handle three-face intersections.

---

## 9. Adaptive Octree Construction

The octree (`src/octree.rs`) recursively subdivides a bounding box, refining
only near the surface.

### 9.1 Data Structures

```rust
struct OctreeCell {
    bbox:     BBox3,
    depth:    u8,
    corners:  [f64; 8],                    // SDF at 8 corners
    children: Option<Box<[OctreeCell; 8]>>, // None for leaves
}

struct Octree {
    root: OctreeCell,
}
```

### 9.2 Build Algorithm

`Octree::build(node, bbox, settings)` calls `build_cell()` recursively:

```
build_cell(node, bbox, depth, corners, settings):
    has_sign_change = any corner differs in sign from corner[0]

    // Interval arithmetic pruning
    if depth >= min_depth AND NOT has_sign_change:
        interval = node.interval_compute(bbox)
        if interval.definitely_positive() OR interval.definitely_negative():
            return LEAF  // Entirely inside or outside -- prune
        else:
            interval_spans_zero = true  // Surface nearby but between corners

    // Subdivision decision
    if depth < max_depth AND (has_sign_change OR depth < min_depth
                              OR interval_spans_zero):
        subdivide into 8 octants
        run corners of each child
        recurse for each child

    return LEAF
```

### 9.3 The `interval_spans_zero` Heuristic

When interval arithmetic reports that the surface *might* pass through a cell
but no corner shows a sign change, the surface is between corners. This
occurs at acute concave Union edges where the narrow positive gap is smaller
than the cell size. Setting `interval_spans_zero = true` forces subdivision
so finer cells can eventually resolve the sign change.

### 9.4 Complexity

- **Best case**: O(N^2) cells where N is the resolution along one axis
  (surface area scaling)
- **Worst case**: O(N^3) when interval arithmetic cannot prune (e.g., Custom
  nodes returning `Interval::entire()`)
- **min_depth** guarantees a minimum resolution everywhere (prevents
  missing thin features at coarse levels)
- **max_depth** caps the recursion (default 8 gives cells as small as
  bbox_size / 256)

---

## 10. Dual Contouring Pipeline

`src/dual_contouring.rs` implements the full mesh extraction pipeline. Two
parallel code paths exist:

### 10.1 Adaptive Path (`extract_mesh_adaptive`)

For `SdfNode` inputs -- uses interval-arithmetic octree pruning and
analytical gradients.

```
SdfNode -> Octree::build() -> surface_cells() -> compute_cell_vertex()
         -> face generation -> TriangleMesh
```

### 10.2 Trait-Object Path (`extract_mesh_from_sdf`)

For `&dyn Sdf` inputs -- no interval arithmetic available. Uses a
near-surface heuristic instead.

```
&dyn Sdf -> build_cell_from_sdf() -> collect surface cells
          -> compute_cell_vertex_sdf() -> face generation -> TriangleMesh
```

### 10.3 Near-Surface Heuristic

`build_cell_from_sdf` lacks access to `interval_compute`, so it substitutes
a geometric heuristic: when no corner has a sign change but the closest corner
value to zero is less than the cell diagonal, the surface likely passes through
the cell. This forces subdivision.

```
near_surface = (NOT has_sign_change AND depth >= min_depth)
               AND min(|corner_i|) < cell_diagonal
```

This heuristic resolves blank patches at acute Union edges where all 8 corners
are inside the union but the narrow positive gap between the constituent shapes
is sub-cell.

### 10.4 Pipeline Phases

Both paths share the same downstream phases:

**Phase 1**: Build octree / recursive cell tree
**Phase 2**: For each surface leaf cell, compute a QEF vertex (Section 11)
**Phase 3**: Generate faces from shared edges (Section 12)
**Phase 4**: Compute per-vertex normals (analytical or central differences)

---

## 11. QEF Vertex Placement

The Quadratic Error Function (QEF) places one vertex per surface cell at
the point that best satisfies all the surface constraints from edge crossings.

### 11.1 Edge Crossing Detection

Each cell has 12 edges (4 parallel to each axis). The canonical edge table
`CELL_EDGES` maps to corner index pairs:

```
X-edges: (0,1), (2,3), (4,5), (6,7)
Y-edges: (0,2), (1,3), (4,6), (5,7)
Z-edges: (0,4), (1,5), (2,6), (3,7)
```

For each edge where the two corner SDF values have opposite signs, a
**bisection search** (32 iterations) finds the zero-crossing point to
within `edge_tolerance` (default `1e-6`).

### 11.2 QEF Formulation

Given `n` edge crossings with positions `p_i` and surface normals `n_i`,
minimize:

```
E(v) = sum_i (n_i . (v - p_i))^2
```

This is a linear least-squares problem. The normal equation is:

```
ATA * v = ATb
```

where `ATA = sum n_i * n_i^T` (3x3 matrix of outer products) and
`ATb = sum n_i * (n_i . p_i)`.

### 11.3 Regularization and Clamping

**Tikhonov regularization**: `ATA += I * 1e-3`, `ATb += mass_point * 1e-3`.
This biases the solution toward the centroid of the crossing positions,
preventing vertex drift when the normal system is poorly conditioned
(e.g., all normals nearly parallel).

**SVD solver**: `nalgebra::SVD` with threshold `1e-10` for numerical
stability. Falls back to the mass point if the SVD fails.

**Bounding box clamp**: The final vertex is clamped to the cell's bounding
box, preventing vertices from drifting outside their cell.

---

## 12. Face Generation and Winding

### 12.1 Edge Key System

Each edge crossing is registered in an `EdgeKey -> Vec<(CellKey, vertex_index)>`
map. The `EdgeKey` consists of the two endpoint coordinates discretized to
a `2^20` grid, ordered lexicographically. This ensures that the same edge
shared by multiple cells maps to a single key.

### 12.2 CellKey

`CellKey` is the discretized minimum corner of a cell plus its depth:
`(ix, iy, iz, depth)`. The `2^20` discretization scale provides enough
precision for depth-8 octrees while fitting in integers.

### 12.3 Polygon Generation

For each `EdgeKey` with >= 2 unique cells:

1. Collect the QEF vertices from all cells sharing this edge
2. Project each vertex onto the plane perpendicular to the edge at its
   midpoint
3. Sort by `atan2` angle around the edge axis
4. Emit a triangle fan from the sorted ring via `emit_fan()`

### 12.4 Winding Order

The first corner of the `EdgeKey` determines winding:
- If the first corner's SDF value is **positive** (outside): emit triangles
  CCW giving outward-facing normals
- If **negative** (inside): emit CW giving normals still point outward

This convention ensures consistent outward normals across the entire mesh.

---

## 13. STEP Export Tiers

`src/step_export.rs` implements a tiered export strategy for AP203 STEP files.

### 13.1 Classification

```rust
pub fn classify(node: &SdfNode) -> ExportTier {
    // Tier 1 (Exact): Sphere, Box3, Cylinder
    //   -- optionally wrapped in Translate/Scale
    //   -- Rotate only for Sphere (rotational symmetry)
    // Tier 3 (Tessellated): Everything else
}
```

`unwrap_transforms()` peels `Translate` and `Scale` layers from the DAG,
accumulating an offset and scale factor. The inner primitive is then matched
against the exact-export whitelist.

### 13.2 Exact BRep Export

**Sphere**: `SPHERICAL_SURFACE` with `AXIS2_PLACEMENT_3D`. Topology: two
pole vertices, one seam circle edge, two oriented edges forming a closed
loop, one `ADVANCED_FACE`.

**Box**: 8 `CARTESIAN_POINT` / `VERTEX_POINT`, 12 `EDGE_CURVE` with `LINE`
geometry, 6 `ADVANCED_FACE` each with `FACE_OUTER_BOUND`, `PLANE`, and
`AXIS2_PLACEMENT_3D`. Edge orientation (`.T.`/`.F.`) encodes face winding.

**Cylinder**: `CYLINDRICAL_SURFACE` for the lateral face, two `PLANE`
surfaces for caps. Two seam vertices, circle edge curves, a seam line edge,
three `ADVANCED_FACE` entities.

### 13.3 Tessellated Fallback

For all other shapes: extract mesh via `extract_mesh_adaptive` with default
settings, then write each triangle as an `ADVANCED_FACE` with a `PLANE`
surface. This produces valid STEP geometry but loses the parametric surface
information.

### 13.4 STEP Metadata

All exports include:
- `GEOMETRIC_REPRESENTATION_CONTEXT(3)` with uncertainty `1e-7`
- `LENGTH_UNIT` as `SI_UNIT(.MILLI., .METRE.)`
- `PRODUCT`, `PRODUCT_DEFINITION`, `PRODUCT_DEFINITION_FORMATION`
- `SHAPE_DEFINITION_REPRESENTATION`
- Schema: `AUTOMOTIVE_DESIGN` (AP203)

---

## 14. Transport and Parametric Paths

The transport system generates swept shapes by moving a cross-section along
a parametric path.

### 14.1 Path Trait

```rust
pub trait Path: Send + Sync {
    fn point(&self, t: f64) -> Vector3<f64>;    // t in [0, 1]
    fn tangent(&self, t: f64) -> Vector3<f64>;  // default: central differences
}
```

Implementations:
- `LinePath { start, end }` -- linear interpolation
- `HelixPath { radius, pitch, turns }` -- circular helix
- `SpiralPath { radius_fn, height_fn, turns }` -- closure-based variable spiral

### 14.2 Transport Orders

| Order | Name | Capabilities |
|---|---|---|
| 0 | Static | No sweep -- the shape IS the section |
| 1 | Rigid + Scale | Sweeps a circle along a path with per-parameter taper |
| 2 | Frenet Frame | Constant-radius tube (delegates to order1 with scale = 1.0) |
| 3 | Section Morphing | Taper + twist (for circular sections, twist is invisible) |

All transport functions return `TransportShape` which implements `Sdf`.

### 14.3 Capsule-Sweep Algorithm (Order 1/2/3)

The path is discretized into `samples` segments. For each segment, a capsule
SDF is computed:

1. Project the query point onto the segment
2. Linearly interpolate the radius between the segment endpoints
   (accounting for the scale function)
3. Return `|p - closest_point| - interpolated_radius`

Adjacent capsules are joined with polynomial smooth-min:

```
smin(a, b, k=0.3) = min(a, b) - (k - |a-b|)^2 / (4*k)
```

This provides C1 continuity at segment boundaries with a maximum shape
deviation of `k/4 = 0.075`.

### 14.4 Frenet Frames

`FrenetFrame` (`src/frame.rs`) provides a local coordinate system at each
path point:

```rust
struct FrenetFrame {
    tangent:  Vector3<f64>,
    normal:   Vector3<f64>,
    binormal: Vector3<f64>,
}
```

`from_tangent()` constructs a robust frame using an arbitrary-perpendicular
fallback when the tangent is axis-aligned.

---

## 15. Voxelization

`src/voxel.rs` provides a regular grid SDF representation.

```rust
struct VoxelGrid {
    resolution: [usize; 3],    // [nx, ny, nz]
    voxel_size: f64,
    origin:     Vector3<f64>,  // world position of grid minimum corner
    data:       Vec<f32>,      // flat row-major [x][y][z]
}
```

Index mapping: `index = ix * ny * nz + iy * nz + iz`

Construction (`Shape::voxelize`):
1. Compute bounding box, pad by one voxel on each side
2. Compute grid dimensions from `bbox_size / voxel_size`
3. Parallel computation via Rayon `into_par_iter()` over all voxel centers
4. Store `node.compute(center) as f32`

Methods:
- `index_at(world) -> usize`: World coordinate to flat index (clamped)
- `world_at(ix, iy, iz) -> Vector3`: Grid indices to world center

---

## 16. Feature Introspection System

The feature system (`src/feature.rs`) provides topological metadata for
primitives, enabling targeted fillet/chamfer operations.

### 16.1 Data Types

```rust
struct FaceInfo {
    index:  usize,
    label:  String,              // "+X", "-Y", "side", "top", etc.
    normal: Option<Vector3<f64>>, // None for curved faces
}

struct EdgeInfo {
    index:  usize,
    face_a: usize,
    face_b: usize,
    label:  String,              // "+X/+Y", "side/top", etc.
}

struct FeatureTarget {
    component: usize,
    body:      usize,
    kind:      FeatureKind,      // Face or Edge
    indices:   Vec<usize>,       // empty = all
}
```

### 16.2 Per-Primitive Feature Tables

**Box3**: 6 faces (`+X`, `-X`, `+Y`, `-Y`, `+Z`, `-Z`), 12 edges
(all face-pair combinations)

**Cylinder**: 3 faces (`side`, `top`, `bottom`), 2 edges
(`side/top`, `side/bottom`)

**Sphere**: 1 face (`surface`), 0 edges

### 16.3 Feature Target Builder

```rust
ft(component, body).edges(&[0, 1, 2])   // specific edges
ft(component, body).faces(&[0, 2])       // specific faces
ft(component, body).all_edges()          // all edges
ft(component, body).all_faces()          // all faces
```

Transforms and Fillet/Chamfer nodes delegate `face_info()` and `edge_info()`
to their inner node, preserving the feature table through the DAG.

---

## 17. Interval Arithmetic

`Interval` (`src/types.rs`) is a closed interval `[lo, hi]` with conservative
arithmetic for bounding SDF values over regions.

### 17.1 Operations

| Operation | Result |
|---|---|
| `a + b` | `[a.lo + b.lo, a.hi + b.hi]` |
| `a - b` | `[a.lo - b.hi, a.hi - b.lo]` |
| `a * b` | `[min(products), max(products)]` over all 4 endpoint combinations |
| `abs(a)` | `[0, max(\|lo\|, \|hi\|)]` if `a` spans zero, else tighter bound |
| `sqrt(a)` | `[sqrt(max(lo, 0)), sqrt(hi)]` |
| `min(a, b)` | `[min(a.lo, b.lo), min(a.hi, b.hi)]` |
| `max(a, b)` | `[max(a.lo, b.lo), max(a.hi, b.hi)]` |
| `neg(a)` | `[-a.hi, -a.lo]` |
| `union(a, b)` | `[min(a.lo, b.lo), max(a.hi, b.hi)]` |

### 17.2 Predicates

- `definitely_positive()`: `lo > 0` -- entire interval is positive
- `definitely_negative()`: `hi < 0` -- entire interval is negative
- `contains_zero()`: `lo <= 0 && hi >= 0`

### 17.3 Conservative Guarantees

Interval arithmetic is always **conservative**: the true range of `f` over
a region is always contained within the computed interval. This means:
- If `definitely_positive()`, the surface provably does not enter the region
- If `definitely_negative()`, the region is provably entirely inside
- If `contains_zero()`, the surface *might* be present (may be a false positive)

The octree exploits this: false positives cause unnecessary subdivision
(wasting time), but there are never false negatives (missing geometry).

---

## 18. Core Data Types

### 18.1 TriangleMesh

```rust
struct TriangleMesh {
    vertices: Vec<Vector3<f64>>,  // vertex positions
    normals:  Vec<Vector3<f64>>,  // per-vertex unit normals
    indices:  Vec<u32>,           // every 3 = one triangle
}
```

`to_binary()` serializes to a compact format for WebSocket transmission:
```
[nv: u32][vertices: f32 x 3 x nv][normals: f32 x 3 x nv][ni: u32][indices: u32 x ni]
```

### 18.2 BBox3

```rust
struct BBox3 {
    min: Vector3<f64>,
    max: Vector3<f64>,
}
```

Key methods:
- `center()`, `size()`, `contains(point)`
- `corners() -> [Vector3<f64>; 8]` -- canonical corner ordering
- `octants() -> [BBox3; 8]` -- subdivide into 8 child boxes

### 18.3 MeshSettings

```rust
struct MeshSettings {
    max_depth:      u8,    // default 8
    min_depth:      u8,    // default 3
    edge_tolerance: f64,   // default 1e-6
}
```

`max_depth = 8` yields cells as small as `bbox_size / 256` along each axis.
`min_depth = 3` ensures at least 8 subdivisions before interval pruning
activates, preventing thin features from being missed at coarse levels.

---

## 19. Architectural Dependency Graph

```
builder.rs -----------------> dag.rs
  Shape                        SdfNode
  |                             |
  |-- .mesh()                   |-- compute()        <-- primitives.rs
  |    \-> dual_contouring.rs   |                        csg.rs
  |         |-- octree.rs       |-- interval_compute <-- types.rs (Interval)
  |         |    \- interval    |-- gradient()
  |         |-- qef.rs (SVD)    |-- face_info()      <-- feature.rs
  |         \-- TriangleMesh    \-- face_distance()   <-- blend.rs
  |                                                       BlendProfile
  |-- .export_obj() --> obj_export.rs
  |-- .export_stl() --> export.rs
  |-- .export_step() -> step_export.rs
  |                      |-- Exact BRep (Sphere/Box/Cylinder)
  |                      \-- Tessellated fallback
  |
  |-- .voxelize() ----> voxel.rs (Rayon parallel)
  |
  \-- .distance()/.contains() -- direct compute()

transport.rs ---> path.rs (LinePath, HelixPath, SpiralPath)
  order0/1/2/3    frame.rs (FrenetFrame)
  \-- TransportShape: impl Sdf

shape.rs --- trait Sdf / Sdf2d
  Generic wrappers for all primitives, CSG, transforms
  Used by extract_mesh_from_sdf (dyn Sdf path)
```

The DAG is the nexus: every upstream API produces `SdfNode` trees, and every
downstream consumer (meshing, export, voxelization) reads them.
