# Crusst: From SDF Library to Geometry Kernel

**Date:** 2026-02-28
**Author:** Thread & Signal LLC — Engineering Roadmap
**Status:** Honest Assessment + Actionable Plan

---

## 0. What Crusst Is Right Now

6,682 lines of Rust across 19 modules. Here is what each module actually does
versus what a geometry kernel requires it to do.

| Module | Lines | What It Does | What a Kernel Needs |
|--------|-------|-------------|-------------------|
| `primitives.rs` | 152 | Quilez SDF distance functions | Exact analytic surface representations |
| `csg.rs` | 84 | Scalar min/max/smooth blend | Topology-tracked boolean operations |
| `dag.rs` | 1208 | SDF expression tree with eval/interval/gradient | This stays — it's the computational engine |
| `dual_contouring.rs` | 692 | Adaptive DC mesh extraction | Manifold-guaranteed, T-junction-free meshing |
| `qef.rs` | 59 | SVD least-squares vertex placement | Fine as-is, possibly needs biased QEF |
| `octree.rs` | 196 | Adaptive subdivision with interval pruning | Fine as-is |
| `blend.rs` | 682 | 10 blend profiles with Newton solvers | **The differentiator. Keep and extend.** |
| `builder.rs` | 863 | Fluent API wrapping Arc\<SdfNode\> | Needs topology-aware construction |
| `shape.rs` | 656 | Trait-based SDF API | Keep for backwards compat |
| `step_export.rs` | 1249 | Triangle soup STEP + 3 exact primitives | Real NURBS STEP with trimmed surfaces |
| `feature.rs` | 99 | Hardcoded face/edge indices per primitive | Topological feature graph |
| `transport.rs` | 159 | Capsule sweep with smin | Frenet frame sweep with section evolution |
| `types.rs` | 235 | BBox3, Interval, MeshSettings, TriangleMesh | Half-edge mesh, topological types |
| `export.rs` | 54 | STL writer | Fine |
| `obj_export.rs` | 39 | OBJ writer | Fine |
| `mesh.rs` | 46 | TriangleMesh binary serialization | Fine |
| `path.rs` | 107 | Line, Helix, Spiral parametric curves | Needs B-spline curves |
| `frame.rs` | 39 | Frenet frame stub | Needs rotation-minimizing frame |
| `voxel.rs` | 44 | Rayon-parallel voxelization | Fine |

### The Core Confusion

Crusst is an **SDF evaluation and meshing library**. A geometry kernel is a
system that maintains **exact geometric representations** (surfaces, curves)
and **topological relationships** (which faces share which edges, which
edges bound which faces) and can perform **boolean operations** that produce
new topology.

The SDF engine is not the problem — it's the right computational foundation.
The problem is that there is no topology layer, no exact surface
representation, and no way to go from "I computed a distance field" to
"here is a BREP solid with trimmed NURBS faces."

The path forward is **not** to throw away the SDF engine and rewrite
OpenCASCADE. It is to build a topology and surface-fitting layer **on top
of** the SDF engine.

```
What Crusst has:          What a kernel needs:

  SdfNode DAG               SdfNode DAG          <-- keep
      |                         |
  evaluate(point)           evaluate(point)       <-- keep
      |                         |
  Octree + DC               Octree + DC           <-- fix (manifold)
      |                         |
  TriangleMesh              HalfEdgeMesh          <-- new
      |                         |
  triangle soup STEP        Surface Fitting       <-- new
                                |
                            Trimmed NURBS Faces   <-- new
                                |
                            BREP Solid            <-- new
                                |
                            Real STEP AP214       <-- new
```

---

## Phase 1: Fix the Mesh (Weeks 1–3)

The dual contouring output is the foundation everything else builds on.
Right now it has three problems that must be fixed before any topology
layer can trust it.

### 1.1 T-Junction Elimination

**Problem:** When adjacent octree cells differ in depth (one is subdivided,
its neighbor is not), the shared face has vertices on the fine side that
have no matching vertex on the coarse side. This creates T-junctions —
edges that terminate in the middle of a neighboring face. The mesh is not
manifold.

**What manifold means:**

```
manifold_condition = for every edge E in mesh:
    count(faces sharing E) == 2

non_manifold_edge = edge shared by 1 face (boundary) or 3+ faces (junction)
```

**Algorithm — Dual Contouring with Restricted Octree:**

The standard fix is to enforce a **2:1 balance constraint** on the octree:
no cell may be adjacent to a cell that differs by more than one level of
refinement.

```
balance_octree(cell):
    for each neighbor N of cell:
        if N.depth < cell.depth - 1:
            subdivide(N)
            // recurse: subdividing N may violate balance for N's neighbors
            balance_octree(N)
```

After balancing, the DC face generation step already handles the 2:1 case
correctly because each shared edge has at most 4 cells around it (not
arbitrarily many).

**Where to implement:** `octree.rs`, new function `Octree::balance()` called
after `Octree::build()` and before DC extraction.

**Test:** Generate a sphere at depth 8. Count non-manifold edges (edges
shared by != 2 triangles). Must be zero.

### 1.2 Watertight Guarantee

**Problem:** The current `emit_fan` function skips edges with only 2 cells
(line 588: "Degenerate: skip (boundary of the octree domain)"). This
creates holes at the domain boundary.

**Fix:** Expand the meshing bounding box to guarantee the surface is fully
contained. `compute_bbox` in `builder.rs` already does a 5% pad for STEP
export. Apply this universally.

For edges at the domain boundary with exactly 2 cells: emit a single
triangle (degenerate quad). The boundary should only be reached if the
bounding box is too tight, which the padding prevents.

**Test:** For every edge in the output mesh, assert exactly 2 incident
triangles. For every vertex, assert the incident triangles form a closed
fan (Euler characteristic check: V - E + F = 2 for a closed genus-0 mesh).

### 1.3 Consistent Winding

**Problem:** The winding-order logic in `extract_mesh_adaptive` (lines
116–135) searches for matching edge keys by looping over all 12 cell edges.
This is fragile — if the edge key quantization creates collisions, winding
can flip.

**Fix:** Compute winding from the SDF gradient at the edge midpoint. The
gradient points outward (positive SDF direction). The triangle normal
(computed via cross product of two edges) should agree with the gradient.
If it doesn't, flip the triangle.

```
for each triangle (v0, v1, v2):
    centroid = (v0 + v1 + v2) / 3
    triangle_normal = (v1 - v0).cross(v2 - v0).normalize()
    sdf_gradient = node.gradient(centroid)
    if triangle_normal.dot(sdf_gradient) < 0:
        swap(v1, v2)   // flip winding
```

**Where to implement:** Post-processing pass in `dual_contouring.rs` after
face generation, before returning the `TriangleMesh`.

### 1.4 Edge Key Precision

**Problem:** Edge keys are quantized to 20-bit integers (line 521:
`(1u64 << 20) as f64`). For a bounding box of size 10, the resolution is
10 / 2^20 = 9.5e-6. This is fine for most cases but can create collisions
at depth 8+ where cell sizes approach 10 / 2^8 = 0.039 — only ~4000
quantization steps per cell edge.

**Fix:** Increase to 30-bit quantization (resolution 9.3e-9) or switch to
exact rational edge keys using the cell's octree address.

---

## Phase 2: Half-Edge Topology (Weeks 4–7)

Once the mesh is manifold and watertight, build a topological data structure
on top of it.

### 2.1 The Half-Edge Data Structure

This is the standard representation for manifold meshes in computational
geometry. Every edge in the mesh is stored as **two directed half-edges**
pointing in opposite directions.

```rust
/// Unique identifier for topological entities.
type VertexId = u32;
type HalfEdgeId = u32;
type FaceId = u32;

/// A directed half-edge in the mesh.
struct HalfEdge {
    /// Vertex this half-edge points TO.
    target: VertexId,
    /// The oppositely-directed half-edge on the neighboring face.
    twin: HalfEdgeId,
    /// Next half-edge around this face (counter-clockwise).
    next: HalfEdgeId,
    /// Previous half-edge around this face (counter-clockwise).
    prev: HalfEdgeId,
    /// The face this half-edge borders.
    face: FaceId,
}

/// A vertex in the topological mesh.
struct Vertex {
    /// Position in 3D space.
    position: Vector3<f64>,
    /// Normal vector (from SDF gradient).
    normal: Vector3<f64>,
    /// Any outgoing half-edge from this vertex.
    half_edge: HalfEdgeId,
}

/// A face (triangle or polygon) in the topological mesh.
struct Face {
    /// Any half-edge on the boundary of this face.
    half_edge: HalfEdgeId,
    /// Which surface region this face belongs to (see Phase 3).
    surface_id: Option<SurfaceId>,
}

/// The complete topological mesh.
struct HalfEdgeMesh {
    vertices: Vec<Vertex>,
    half_edges: Vec<HalfEdge>,
    faces: Vec<Face>,
}
```

**Why half-edge and not winged-edge:** Half-edge is simpler to implement,
has better cache locality for traversal, and is what CGAL, libigl, and
OpenMesh use. Winged-edge stores less data but requires branch logic to
determine traversal direction.

**Construction from TriangleMesh:**

```
build_half_edge_mesh(tri_mesh):
    // For each triangle, create 3 half-edges
    for each triangle (i0, i1, i2) with face_id F:
        he_a = new HalfEdge { target: i1, face: F, ... }
        he_b = new HalfEdge { target: i2, face: F, ... }
        he_c = new HalfEdge { target: i0, face: F, ... }
        he_a.next = he_b; he_b.next = he_c; he_c.next = he_a
        he_a.prev = he_c; he_b.prev = he_a; he_c.prev = he_b

    // Link twins: for each pair of half-edges sharing the same
    // undirected edge (v_a, v_b), set them as twins.
    // Use a HashMap<(min(v_a,v_b), max(v_a,v_b)), HalfEdgeId>
    // to find matching pairs.
    edge_map: HashMap<(VertexId, VertexId), HalfEdgeId>
    for each half_edge he:
        key = (min(he.source, he.target), max(he.source, he.target))
        if key in edge_map:
            twin_id = edge_map[key]
            he.twin = twin_id
            half_edges[twin_id].twin = he.id
        else:
            edge_map[key] = he.id
```

**Where to implement:** New module `src/topology.rs`.

**Invariant tests:**
- Every half-edge has a twin (no boundary for closed meshes)
- `he.twin.twin == he` for all half-edges
- Following `next` from any half-edge returns to itself (closed face loop)
- `he.twin.target == he.prev.target` (source vertex consistency)

### 2.2 Face Region Classification

After building the half-edge mesh, classify faces into **surface regions**:
groups of connected faces that lie on the same geometric surface.

For SDF primitives, you already know the surface type from the DAG node.
The challenge is that after boolean operations, a single primitive's
surface gets split into regions.

**Algorithm — Region Growing with Normal Coherence:**

```
classify_faces(mesh, sdf_node):
    // Step 1: For each face, determine which primitive generated it
    for each face F:
        centroid = face_centroid(F)
        F.source_primitive = identify_primitive(sdf_node, centroid)
        // identify_primitive walks the DAG and finds which leaf
        // primitive has the smallest absolute SDF value at centroid

    // Step 2: Region-grow connected faces with same source primitive
    // and similar normals (within angular tolerance)
    regions = []
    visited = set()
    for each face F not in visited:
        region = flood_fill(F, same_source_and_normal_tolerance)
        regions.append(region)
        visited.add_all(region)

    return regions
```

The `identify_primitive` function is the key piece. It walks the SDF DAG:

```
identify_primitive(node, point) -> (PrimitiveId, f64):
    match node:
        Sphere/Box3/Cylinder/...:
            return (node.id, node.evaluate(point).abs())
        Union(a, b):
            (id_a, dist_a) = identify_primitive(a, point)
            (id_b, dist_b) = identify_primitive(b, point)
            return if dist_a < dist_b then (id_a, dist_a) else (id_b, dist_b)
        Difference(a, b):
            (id_a, dist_a) = identify_primitive(a, point)
            (id_b, dist_b) = identify_primitive(b, point)
            return if dist_a < dist_b then (id_a, dist_a) else (id_b, dist_b)
        Translate(inner, offset):
            return identify_primitive(inner, point - offset)
        // ... etc for other transforms
```

**Where to implement:** New module `src/surface_classify.rs`.

### 2.3 Edge Classification

Once faces are grouped into surface regions, the edges between regions
become **topological edges** — the curves where two surfaces meet.

```
for each half_edge he in mesh:
    face_a = he.face
    face_b = he.twin.face
    if face_a.surface_id != face_b.surface_id:
        mark he as a BOUNDARY_EDGE between regions
    else:
        mark he as an INTERIOR_EDGE within a region
```

Chains of connected boundary half-edges form **edge loops** — the trimming
curves that bound each surface region. This is exactly the data that STEP
BREP export needs.

---

## Phase 3: Surface Fitting (Weeks 8–14)

This is the hard part. Each surface region from Phase 2 needs to be
represented as an exact geometric surface, not a bag of triangles.

### 3.1 Analytic Surface Fitting

For faces that came from known SDF primitives, the exact surface is already
known:

| SDF Primitive | Exact Surface Type | STEP Entity |
|--------------|-------------------|-------------|
| Sphere | `SphericalSurface(center, radius)` | `SPHERICAL_SURFACE` |
| Box3 (each face) | `Plane(point, normal)` | `PLANE` |
| Cylinder (side) | `CylindricalSurface(axis, radius)` | `CYLINDRICAL_SURFACE` |
| Cylinder (caps) | `Plane(point, normal)` | `PLANE` |
| CappedCone (side) | `ConicalSurface(axis, half_angle)` | `CONICAL_SURFACE` |
| CappedCone (caps) | `Plane(point, normal)` | `PLANE` |
| Torus | `ToroidalSurface(center, major, minor)` | `TOROIDAL_SURFACE` |

These are the **Tier 1** surfaces — closed-form, exact, no approximation.

```rust
enum AnalyticSurface {
    Plane { point: Vector3<f64>, normal: Vector3<f64> },
    Sphere { center: Vector3<f64>, radius: f64 },
    Cylinder { origin: Vector3<f64>, axis: Vector3<f64>, radius: f64 },
    Cone { apex: Vector3<f64>, axis: Vector3<f64>, half_angle: f64 },
    Torus { center: Vector3<f64>, axis: Vector3<f64>, major: f64, minor: f64 },
}
```

**Where to implement:** New module `src/surface.rs`.

For regions that trace back to an analytic primitive, fitting is trivial:
read the parameters directly from the SdfNode.

### 3.2 B-Spline Surface Fitting (Blend Regions)

Faces in blend zones — fillets, chamfers, smooth CSG transitions — do not
correspond to any analytic surface. These need B-spline (NURBS) surface
approximation.

**This is the hardest single problem in the entire roadmap.**

The algorithm:

```
fit_bspline_to_region(region_faces, sdf_node, tolerance):
    // 1. Collect sample points on the SDF zero-surface within the region
    points = []
    normals = []
    for each face F in region_faces:
        for each vertex V of F:
            // Project vertex onto exact zero-surface via Newton iteration
            p = project_to_surface(V.position, sdf_node)
            n = sdf_node.gradient(p).normalize()
            points.append(p)
            normals.append(n)

    // 2. Parameterize: map 3D points to (u, v) parameter space
    // Use conformal parameterization (LSCM or ABF++) of the
    // triangulated region
    (u, v) = parameterize_region(region_faces)

    // 3. Fit B-spline surface to (u, v) -> (x, y, z) mapping
    // Least-squares fit with knot insertion for error control
    surface = fit_nurbs_surface(points, normals, u, v, degree=3)

    // 4. Refine: insert knots where fitting error > tolerance
    while max_error(surface, points) > tolerance:
        insert_knots_at_worst_error(surface)
        refit(surface, points, u, v)

    return surface
```

**Key sub-problems:**

**Parameterization** — Mapping a 3D surface patch to a 2D (u,v) domain.
Least Squares Conformal Maps (LSCM) is the standard approach. It minimizes
angular distortion and is implemented in libigl. The Rust crate
`nalgebra` provides the sparse linear algebra needed. This is a
sparse-matrix eigenvalue problem:

```
minimize ||Jf - Rotation * Jg||^2

where:
    Jf = Jacobian of the 3D -> 2D map
    Jg = identity (conformal = locally rotation + scale)
```

**NURBS fitting** — Given parameterized points, fit a B-spline surface.
This is a linear least-squares problem per control point:

```
surface(u,v) = sum_i sum_j N_i(u) * N_j(v) * P_ij

minimize sum_k ||surface(u_k, v_k) - point_k||^2

// This is a linear system: A * P = B
// where A_ki = N_i(u_k) * N_j(v_k)
//       B_k  = point_k
```

The `nalgebra` SVD or QR solver handles this directly. No new dependencies
needed.

**Knot insertion** — When fitting error exceeds tolerance, insert new
knots in the B-spline knot vector at the parameter location of maximum
error, then refit. This is the Boehm knot insertion algorithm.

**Crate options:**
- Implement from scratch using `nalgebra` (recommended — you need to own
  this math for the kernel to be genuine)
- `nurbs` crate exists but is minimal
- `bspline` crate exists but is 2D only

**Where to implement:** New modules `src/nurbs.rs` (B-spline evaluation,
knot insertion, fitting) and `src/parameterize.rs` (LSCM).

**Estimated scope:** ~2000–3000 lines for a working NURBS surface fitter.
This is the single largest new module.

### 3.3 Surface Projection

For both analytic and B-spline surfaces, you need a function that projects
an arbitrary point onto the surface (closest point query). This is used
for:

- Refining mesh vertices to lie exactly on the surface
- Computing trimming curve intersections
- Validating fitting accuracy

```rust
trait Surface {
    /// Evaluate surface position at parameter (u, v).
    fn evaluate(&self, u: f64, v: f64) -> Vector3<f64>;

    /// Surface normal at parameter (u, v).
    fn normal(&self, u: f64, v: f64) -> Vector3<f64>;

    /// Partial derivatives du, dv at parameter (u, v).
    fn partials(&self, u: f64, v: f64) -> (Vector3<f64>, Vector3<f64>);

    /// Find the parameter (u, v) closest to a 3D point.
    /// Newton iteration on the distance function.
    fn project(&self, point: Vector3<f64>) -> (f64, f64, Vector3<f64>);
}
```

For analytic surfaces, `project` has closed-form solutions. For B-spline
surfaces, it requires Newton iteration on:

```
minimize ||S(u,v) - P||^2

gradient:  [dS/du . (S - P),  dS/dv . (S - P)]
hessian:   standard Newton-Raphson on the 2x2 system
```

---

## Phase 4: Trimming Curves and BREP Assembly (Weeks 15–20)

### 4.1 Trimming Curves

Each surface region from Phase 2 has a boundary — a set of edge loops
where it meets neighboring surfaces. These boundaries need to be
represented as **trimming curves** in the surface's (u, v) parameter space.

```
for each surface region R:
    boundary_loops = find_boundary_edge_loops(R, half_edge_mesh)
    for each loop L:
        // Project 3D boundary vertices onto the surface
        uv_points = []
        for each vertex V in L:
            (u, v, _) = R.surface.project(V.position)
            uv_points.append((u, v))

        // Fit a B-spline curve to the (u, v) points
        trim_curve = fit_bspline_curve(uv_points, degree=3)
        R.trim_curves.append(trim_curve)
```

A BREP face is then:

```rust
struct BrepFace {
    /// The underlying geometric surface.
    surface: Box<dyn Surface>,
    /// Outer boundary loop (CCW in parameter space).
    outer_loop: Vec<TrimCurve>,
    /// Inner loops (holes, CW in parameter space).
    inner_loops: Vec<Vec<TrimCurve>>,
    /// Orientation: does the face normal agree with the surface normal?
    same_sense: bool,
}

struct TrimCurve {
    /// B-spline curve in (u, v) parameter space.
    curve_2d: BSplineCurve2d,
    /// Corresponding B-spline curve in 3D (the edge geometry).
    curve_3d: BSplineCurve3d,
    /// Reference to the topological edge.
    edge_id: EdgeId,
}
```

### 4.2 BREP Solid Assembly

A complete BREP solid is:

```rust
struct BrepSolid {
    /// All faces with their surfaces and trimming curves.
    faces: Vec<BrepFace>,
    /// All topological edges (shared between faces).
    edges: Vec<BrepEdge>,
    /// All topological vertices.
    vertices: Vec<BrepVertex>,
    /// Shells: groups of faces forming closed boundaries.
    shells: Vec<Shell>,
}

struct BrepEdge {
    /// The 3D curve along this edge.
    curve: BSplineCurve3d,
    /// Start vertex.
    start: VertexId,
    /// End vertex.
    end: VertexId,
    /// Faces on each side.
    left_face: FaceId,
    right_face: FaceId,
}

struct Shell {
    /// Indices of faces forming this shell.
    faces: Vec<FaceId>,
    /// Is this an outer shell or an inner void?
    orientation: ShellOrientation,
}
```

**Validation — the Euler-Poincaré formula:**

For a valid BREP solid with genus g (number of handles/holes through the
solid):

```
V - E + F = 2 * (S - G)

where:
    V = number of vertices
    E = number of edges
    F = number of faces
    S = number of shells
    G = total genus across all shells
```

For a simple solid (one shell, no holes): V - E + F = 2.

### 4.3 The SDF-to-BREP Pipeline (Complete)

Putting Phases 1–4 together:

```
sdf_to_brep(node: &SdfNode, tolerance: f64) -> BrepSolid:

    // Phase 1: Extract manifold mesh
    mesh = extract_mesh_adaptive(node, bbox, settings)
    mesh = fix_winding(mesh, node)

    // Phase 2: Build topology and classify
    he_mesh = build_half_edge_mesh(mesh)
    regions = classify_face_regions(he_mesh, node)

    // Phase 3: Fit surfaces
    for each region R:
        if R.source is analytic primitive:
            R.surface = extract_analytic_surface(R.source)
        else:
            R.surface = fit_bspline_surface(R.faces, node, tolerance)

    // Phase 4: Extract trimming curves and assemble
    for each region R:
        R.trim_loops = extract_trim_curves(R, he_mesh)

    solid = assemble_brep(regions)
    validate_euler_poincare(solid)

    return solid
```

---

## Phase 5: Real STEP Export (Weeks 21–25)

### 5.1 STEP AP214 Entity Mapping

The current `step_export.rs` writes AP203. Real CAD interop needs AP214
(automotive design) or at minimum AP203 Edition 2 with proper BREP.

**Entity hierarchy for a BREP face:**

```
ADVANCED_BREP_SHAPE_REPRESENTATION
  └── MANIFOLD_SOLID_BREP
        └── CLOSED_SHELL
              └── ADVANCED_FACE (one per BrepFace)
                    ├── FACE_BOUND / FACE_OUTER_BOUND
                    │     └── EDGE_LOOP
                    │           └── ORIENTED_EDGE (one per trim curve)
                    │                 └── EDGE_CURVE
                    │                       ├── VERTEX_POINT (start)
                    │                       ├── VERTEX_POINT (end)
                    │                       └── B_SPLINE_CURVE_WITH_KNOTS (3D edge curve)
                    └── (surface geometry)
                          ├── PLANE / SPHERICAL_SURFACE / CYLINDRICAL_SURFACE (analytic)
                          └── B_SPLINE_SURFACE_WITH_KNOTS (fitted)
```

**What changes in `step_export.rs`:**

The current file is 1249 lines. Roughly 800 of those are the tessellated
fallback (one `ADVANCED_FACE` per triangle with a `PLANE` surface).
That entire path gets replaced by the BREP export path above.

The exact-primitive path (lines 28–56, `classify()`) gets absorbed into
Phase 3's analytic surface fitting — you no longer need a separate code
path because every face goes through the same BREP pipeline.

Keep the tessellated fallback as a debug/preview mode, but the primary
export path is the BREP pipeline.

### 5.2 STEP Validation

Use an external validator to confirm the output is spec-compliant:

- **FreeCAD** — open the STEP file, check that faces render correctly
- **STEP File Analyzer** (NIST) — free tool that validates AP203/AP214
  compliance
- **OpenCASCADE STEP reader** — if it can read the file without errors,
  it's probably correct

Write integration tests that:
1. Export a sphere, box, cylinder as STEP
2. Read back with a STEP parser (write a minimal one or use `step21` crate)
3. Verify entity counts match expected topology
4. Verify vertex positions match within tolerance

---

## Phase 6: Boolean Operations with Topology (Weeks 26–32)

Right now, CSG in Crusst is purely scalar: `min(d1, d2)` for union.
This works for SDF evaluation but produces no topology. A geometry kernel
needs boolean operations that track which faces survive, which edges are
created at intersections, and which regions are inside/outside.

### 6.1 SDF-Guided Surface-Surface Intersection

The traditional BREP boolean approach (Requicha, Mantyla) intersects every
pair of surfaces, trims them, and rebuilds topology. This is thousands of
lines of intersection code and the primary reason OpenCASCADE is 7 million
lines of C++.

Crusst can take a shortcut because it has the SDF:

```
boolean_union(solid_a, solid_b):
    // 1. Create the combined SDF node (this already works)
    combined = SdfNode::Union(solid_a.sdf, solid_b.sdf)

    // 2. Extract a new manifold mesh from the combined SDF
    mesh = extract_mesh_adaptive(combined, bbox, settings)

    // 3. Build topology and classify faces
    //    Face classification now uses the COMBINED dag,
    //    so it automatically determines which faces of A and B
    //    survive the boolean.
    he_mesh = build_half_edge_mesh(mesh)
    regions = classify_face_regions(he_mesh, combined)

    // 4. For surviving faces from A: reuse A's surface fits
    //    For surviving faces from B: reuse B's surface fits
    //    For NEW faces (intersection curves): fit new surfaces
    for each region R:
        if R traces back to solid_a's primitive:
            R.surface = solid_a.find_surface(R.source_id)
        else if R traces back to solid_b's primitive:
            R.surface = solid_b.find_surface(R.source_id)
        else:
            R.surface = fit_bspline_surface(R, combined, tolerance)

    // 5. Extract new trimming curves at boolean boundaries
    // 6. Assemble new BREP
```

**This is the key insight:** The SDF does the hard computational work
(determining what's inside/outside). The BREP layer just needs to
**classify and label** the result, not compute it from scratch.

This avoids the combinatorial explosion of surface-surface intersection
that makes traditional BREP booleans so complex. The cost is that the
result is approximate (mesh-resolution-dependent), not exact. But for
3D printing, visualization, and most engineering applications, a
tolerance of 1e-4 mm is more than sufficient.

### 6.2 Intersection Curve Extraction

Where two surfaces meet after a boolean, there's a new edge that didn't
exist on either input solid. This curve lies on both surfaces
simultaneously.

```
extract_intersection_curve(region_a, region_b, boundary_edges):
    // The boundary edges between region_a and region_b define
    // a polyline approximation of the intersection curve.
    polyline_3d = [edge.midpoint for edge in boundary_edges]

    // Project onto both surfaces to get (u,v) curves
    uv_a = [region_a.surface.project(p) for p in polyline_3d]
    uv_b = [region_b.surface.project(p) for p in polyline_3d]

    // Fit B-spline curves
    curve_3d = fit_bspline_curve(polyline_3d, degree=3)
    pcurve_a = fit_bspline_curve(uv_a, degree=3)  // trim curve on surface A
    pcurve_b = fit_bspline_curve(uv_b, degree=3)  // trim curve on surface B

    return IntersectionCurve { curve_3d, pcurve_a, pcurve_b }
```

---

## Phase 7: What Makes Crusst Crusst (Ongoing)

Everything above is what **any** geometry kernel needs. Here's what makes
Crusst different from "another OpenCASCADE clone" — and these should be
developed in parallel with the infrastructure phases.

### 7.1 Blend Profiles as First-Class BREP Features

The 10 blend profiles in `blend.rs` are the genuine differentiator. No
other kernel offers cycloidal, parabolic, or hyperbolic blends as first-
class operations. In BREP-land, these become **rolling-ball blends** with
non-circular cross-sections.

When a blend region is classified in Phase 2, instead of generic B-spline
fitting, use the known blend profile to generate an exact parametric
surface:

```
for cycloidal blend of radius r between faces A and B:
    // The blend surface is a sweep of the cycloidal profile
    // along the intersection curve of A and B.
    //
    // Cross-section at parameter t along the edge:
    //   cycloid(theta) for theta in [0, pi/2]
    //   mapped into the normal plane of the edge curve at t
    //
    // This gives an exact parametric surface, not an approximation.

    edge_curve = intersection_curve(A, B)
    profile = cycloidal_cross_section(r)
    surface = sweep_profile_along_curve(profile, edge_curve)
```

This means Crusst's STEP output would contain blend surfaces that are
**more accurate** than what you'd get from fitting a generic B-spline to
the mesh. That's a real competitive advantage.

### 7.2 Transport System with Real Frenet Frames

The current `transport.rs` has Order 2 and Order 3 as stubs that
delegate to Order 1 (capsule sweep). A real sweep needs:

**Rotation-Minimizing Frame (RMF):** The Frenet frame twists
unnaturally at inflection points. The RMF (Wang et al., 2008) maintains
consistent orientation by parallel-transporting the frame along the curve.

```rust
/// Compute a rotation-minimizing frame along a parametric curve.
fn compute_rmf(path: &dyn Path, samples: usize) -> Vec<Frame> {
    let mut frames = Vec::with_capacity(samples);

    // Initial frame at t=0
    let t0 = path.tangent(0.0).normalize();
    let r0 = arbitrary_perpendicular(t0);
    let s0 = t0.cross(&r0);
    frames.push(Frame { tangent: t0, normal: r0, binormal: s0 });

    // Double reflection method (Hanson & Ma, 1995)
    for i in 1..samples {
        let t = i as f64 / (samples - 1) as f64;
        let prev = &frames[i - 1];
        let ti = path.tangent(t).normalize();

        // Reflect previous frame through bisector plane
        let v1 = path.point(t) - path.point((i-1) as f64 / (samples-1) as f64);
        let c1 = v1.dot(&v1);
        let ri_l = prev.normal - v1 * (2.0 * v1.dot(&prev.normal) / c1);
        let ti_l = prev.tangent - v1 * (2.0 * v1.dot(&prev.tangent) / c1);

        // Second reflection to align with actual tangent
        let v2 = ti - ti_l;
        let c2 = v2.dot(&v2);
        let ri = ri_l - v2 * (2.0 * v2.dot(&ri_l) / c2);

        frames.push(Frame {
            tangent: ti,
            normal: ri.normalize(),
            binormal: ti.cross(&ri).normalize(),
        });
    }

    frames
}
```

**Non-circular cross-sections for Order 3:**

Replace the capsule sweep with a proper SDF evaluation in the local frame:

```
order3_evaluate(point, path, frames, section_sdf, scale_fn, twist_fn):
    // Find closest point on path
    t_closest = closest_parameter(path, point)
    frame = interpolate_frame(frames, t_closest)

    // Transform point into local frame
    local = frame.inverse() * (point - path.point(t_closest))

    // Apply inverse scale and twist
    s = scale_fn(t_closest)
    theta = twist_fn(t_closest)
    local = rotate_2d(local.xy, -theta) / s

    // Evaluate section SDF in local 2D
    return section_sdf.evaluate(local.xy) * s
```

### 7.3 The Orders of Topology as Kernel Architecture

Your philosophical framework maps directly to kernel capabilities:

| Order | Topology | Geometry | Crusst Module |
|-------|----------|----------|--------------|
| 0 | Point | Vertex | `types.rs` — exists |
| 1 | Curve (1 DOF) | Edge, Path | `path.rs` + new `curve.rs` |
| 2 | Surface (2 DOF) | Face, Region | new `surface.rs` + `nurbs.rs` |
| 3 | Volume (3 DOF) | Solid, Shell | new `brep.rs` |

Each order introduces one degree of freedom. The kernel operations at each
order:

- **Order 0 → 1:** Sweep a point along a curve = edge creation
- **Order 1 → 2:** Sweep a curve along a curve = surface creation (ruled,
  swept, revolution)
- **Order 2 → 3:** Bound surfaces into a closed shell = solid creation
- **Cross-order:** Boolean operations intersect Order-3 solids, creating
  new Order-1 edges at surface intersections

This isn't just philosophy — it's the actual type hierarchy for the kernel.

---

## Phase 8: Performance (Parallel to All Phases)

### 8.1 Tape Compilation (libfive approach)

The DAG evaluator currently does dynamic dispatch via `match` on every
`evaluate()` call. For meshing (millions of evaluations), compile the DAG
into a linear instruction tape:

```rust
enum TapeOp {
    LoadPoint,                    // push query point xyz
    Sphere { center, radius },    // push distance to sphere
    Box3 { center, half },        // push distance to box
    Union,                        // pop 2, push min
    Intersection,                 // pop 2, push max
    SmoothUnion { k },            // pop 2, push smooth_min
    Translate { offset },         // modify point register
    // ...
}

struct Tape {
    ops: Vec<TapeOp>,
}

impl Tape {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        let mut stack: Vec<f64> = Vec::with_capacity(32);
        let mut pt = point;
        for op in &self.ops {
            match op {
                TapeOp::Sphere { center, radius } => {
                    stack.push((pt - *center).norm() - radius);
                }
                TapeOp::Union => {
                    let b = stack.pop().unwrap();
                    let a = stack.pop().unwrap();
                    stack.push(a.min(b));
                }
                // ...
            }
        }
        stack[0]
    }
}
```

Expected speedup: 2–5x over dynamic dispatch due to better branch
prediction and cache locality. The `libfive` paper reports up to 10x with
region-based tape pruning (removing ops that are provably irrelevant for a
given spatial region).

### 8.2 GPU Evaluation (wgpu)

The existing v2 design doc mentions this. The implementation path:

1. WGSL codegen from the tape (not the DAG — tape is linear, maps to GPU
   trivially)
2. Batch dispatch: upload a grid of query points, get back distances +
   gradients
3. GPU-assisted octree: evaluate all 8 child corners in one dispatch

This is a large effort (~3000+ lines) and should come **after** the BREP
pipeline is working on CPU. The BREP pipeline needs correct meshes, not
fast meshes.

### 8.3 Rayon Parallelism

The octree construction evaluates child cells sequentially. Since children
at the same depth are independent, use `rayon::join` for parallel
evaluation at shallow depths:

```
if depth < parallel_threshold:
    let (children_0_3, children_4_7) = rayon::join(
        || build_children(0..4),
        || build_children(4..8),
    );
```

---

## Module Map After All Phases

```
crusst/src/
├── lib.rs                    # module declarations
│
├── ── SDF Engine (exists, keep) ──
├── primitives.rs             # SDF distance functions
├── csg.rs                    # scalar CSG operations
├── dag.rs                    # expression DAG
├── blend.rs                  # 10 blend profiles ← THE DIFFERENTIATOR
├── shape.rs                  # Sdf/Sdf2d traits
├── builder.rs                # fluent API
├── feature.rs                # feature targeting
│
├── ── Meshing (exists, fix) ──
├── octree.rs                 # adaptive octree + balance constraint  [MODIFY]
├── qef.rs                    # QEF solver
├── dual_contouring.rs        # manifold DC extraction                [MODIFY]
│
├── ── Topology (new) ──
├── topology.rs               # HalfEdgeMesh data structure           [NEW]
├── surface_classify.rs       # face region classification            [NEW]
│
├── ── Geometry (new) ──
├── surface.rs                # Surface trait + analytic surfaces      [NEW]
├── nurbs.rs                  # B-spline curves and surfaces           [NEW]
├── parameterize.rs           # LSCM conformal parameterization        [NEW]
├── curve.rs                  # B-spline curve fitting                 [NEW]
│
├── ── BREP (new) ──
├── brep.rs                   # BrepSolid, BrepFace, Shell types       [NEW]
├── trim.rs                   # trimming curve extraction              [NEW]
├── boolean.rs                # SDF-guided BREP booleans               [NEW]
│
├── ── Export (exists, extend) ──
├── step_export.rs            # BREP STEP AP214 export                [REWRITE]
├── export.rs                 # STL
├── obj_export.rs             # OBJ
│
├── ── Paths & Sweeps (exists, extend) ──
├── path.rs                   # parametric curves + B-spline paths     [EXTEND]
├── frame.rs                  # rotation-minimizing frames             [REWRITE]
├── transport.rs              # sweep with real frame transport        [REWRITE]
│
├── ── Performance (new) ──
├── tape.rs                   # compiled instruction tape              [NEW]
│
├── ── Infrastructure (exists, keep) ──
├── types.rs                  # BBox3, Interval, MeshSettings
├── mesh.rs                   # TriangleMesh
└── voxel.rs                  # voxelization
```

Estimated new code: ~8,000–12,000 lines.
Total kernel after all phases: ~16,000–20,000 lines.

For reference: `truck` (Rust BREP kernel) is ~25,000 lines. OpenCASCADE is
~7,000,000 lines. Crusst at 20,000 lines with SDF-guided BREP would be a
lean, focused kernel — which is the right target for a one-person project.

---

## Dependency Decisions

| Need | Current | Recommendation |
|------|---------|---------------|
| Linear algebra | nalgebra 0.33 | Keep. Sufficient for everything including NURBS fitting. |
| Sparse matrices | none | Add `nalgebra-sparse` or `sprs` for LSCM parameterization. |
| Parallelism | rayon 1.10 | Keep. |
| NURBS evaluation | none | Implement from scratch. ~500 lines. You need to own this. |
| Half-edge mesh | none | Implement from scratch. ~400 lines. Standard algorithm. |
| STEP writing | custom | Extend current. The entity writing code is fine; the topology mapping is what changes. |
| GPU (future) | none | `wgpu` when ready. Don't add until CPU pipeline is solid. |

**Zero new heavy dependencies.** The kernel stays pure Rust + nalgebra.
This is a feature — it's what makes it embeddable and trustworthy.

---

## Milestone Deliverables

| Milestone | Deliverable | How to Verify |
|-----------|------------|--------------|
| M1 (Phase 1) | Manifold mesh with zero non-manifold edges | Euler check: V - E + F = 2 for genus-0 |
| M2 (Phase 2) | Half-edge mesh with face regions labeled | Visual: color faces by region, verify grouping |
| M3 (Phase 3) | B-spline surface fits within 1e-4 tolerance | Max deviation of sample points from fitted surface |
| M4 (Phase 4) | BREP solid assembly passing Euler-Poincaré | V - E + F = 2(S - G) |
| M5 (Phase 5) | STEP file opens in FreeCAD with correct faces | Visual + NIST STEP validator |
| M6 (Phase 6) | Boolean union/difference produces valid BREP | Export boolean result, open in Fusion 360 |
| M7 (Phase 7) | Cycloidal fillet exports as exact swept surface | Compare mesh fillet vs BREP fillet in CAD viewer |
| M8 (Phase 8) | Tape compiler with 3x+ speedup on benchmarks | criterion benchmarks vs current DAG evaluator |

---

## What Not To Do

1. **Don't try to implement exact BREP booleans.** Surface-surface
   intersection with all the degenerate cases (tangencies, coincident
   faces, edge-on-edge) is a multi-decade research problem. Use the SDF
   to resolve containment and classify the results.

2. **Don't add OCCT as a dependency.** The whole point is pure Rust. If
   you link OCCT, you're not a kernel — you're a wrapper.

3. **Don't chase GPU before BREP works on CPU.** Fast wrong answers are
   still wrong answers. GPU acceleration is a multiplier on a working
   pipeline, not a substitute for one.

4. **Don't implement IGES.** It's a dead format. STEP is the standard.
   Anyone who needs IGES can convert from STEP.

5. **Don't try to handle non-manifold geometry.** Non-manifold (sheets,
   wires, mixed-dimension) is a massive scope expansion. Crusst should
   target manifold solids only. That's what 95% of CAD and 100% of
   3D printing needs.

6. **Don't undervalue the SDF engine.** The hybrid approach (SDF for
   computation, BREP for representation) is genuinely novel for an
   open-source Rust kernel. It's not a compromise — it's an
   architecture decision that sidesteps the hardest problems in
   traditional BREP while still delivering CAD-quality output.

---

## Reading List

These are the papers and references behind each phase. Not textbook
recommendations — these are the specific algorithms you'll implement.

**Dual Contouring:**
- Ju et al., "Dual Contouring of Hermite Data" (2002) — the original
- Schaefer & Warren, "Dual Marching Cubes" (2004) — manifold variant

**Half-Edge:**
- Kettner, "Using Generic Programming for Designing a Data Structure
  for Polyhedral Surfaces" (1999) — the CGAL approach

**LSCM Parameterization:**
- Lévy et al., "Least Squares Conformal Maps" (2002) — the standard

**B-Spline Fitting:**
- Piegl & Tiller, "The NURBS Book" (1997) — chapters 9–10 (fitting)
- Park & Lee, "B-spline surface fitting" (2007) — knot placement

**STEP Format:**
- ISO 10303-21 (STEP physical file format)
- ISO 10303-42 (geometric and topological representation)
- The NIST STEP File Analyzer documentation (free, includes entity
  relationship diagrams)

**Rotation-Minimizing Frames:**
- Wang et al., "Computation of RMF" (2008)

**SDF-based BREP (the hybrid approach):**
- Fayolle & Pasko, "An evolutionary approach to the generation of
  SDF representations" (2016)
- libfive documentation (implicit CAD kernel architecture)
