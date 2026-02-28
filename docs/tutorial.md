# Crusst Tutorial — From First Principles to Watertight Meshes

This tutorial walks through the Crusst geometry kernel from foundational concepts to advanced usage. Each section builds on the previous one and includes runnable code examples.

## Table of Contents

1. [What is a Signed Distance Function?](#1-what-is-a-signed-distance-function)
2. [Evaluating Primitives](#2-evaluating-primitives)
3. [Constructive Solid Geometry](#3-constructive-solid-geometry)
4. [The Builder API](#4-the-builder-api)
5. [Spatial Transforms](#5-spatial-transforms)
6. [Mesh Extraction](#6-mesh-extraction)
7. [Export Formats](#7-export-formats)
8. [Fillet & Chamfer Operations](#8-fillet--chamfer-operations)
9. [Round & Chamfer CSG](#9-round--chamfer-csg)
10. [2D Shapes, Revolve, and Extrude](#10-2d-shapes-revolve-and-extrude)
11. [Voxelization](#11-voxelization)
12. [Parametric Paths & Transport](#12-parametric-paths--transport)
13. [Custom SDF Functions](#13-custom-sdf-functions)
14. [The DAG Layer](#14-the-dag-layer)
15. [Performance Considerations](#15-performance-considerations)

---

## 1. What is a Signed Distance Function?

A **signed distance function** (SDF) is a mathematical function `f : R³ → R` that, given any point in 3D space, returns the shortest distance to the nearest surface of a solid. The sign encodes containment:

```
f(p) < 0    →   p is inside the solid
f(p) = 0    →   p is on the surface (the zero-level set)
f(p) > 0    →   p is outside the solid
```

**Example — a sphere of radius `r` centered at the origin:**

```
f(p) = |p| - r
```

At the center (`p = 0`): `f = -r` (negative, inside). On the surface (`|p| = r`): `f = 0`. Far away (`|p| = 2r`): `f = r` (positive, outside).

SDFs are powerful because they compose algebraically:
- **Union** of two solids A and B: `f(p) = min(fA(p), fB(p))`
- **Intersection**: `f(p) = max(fA(p), fB(p))`
- **Difference** (A minus B): `f(p) = max(fA(p), -fB(p))`

This means arbitrarily complex geometry can be built from simple primitives using just min, max, and negation.

---

## 2. Evaluating Primitives

Crusst provides two layers for SDF evaluation:

### Functional layer (`primitives` module)

Pure functions that take parameters and return a distance value:

```rust
use crusst::primitives::{sdf_sphere, sdf_box, sdf_cylinder};
use nalgebra::Vector3;

// Sphere: center at origin, radius 5
let d = sdf_sphere(Vector3::new(3.0, 4.0, 0.0), Vector3::zeros(), 5.0);
// |p| = 5, radius = 5 → d = 0 (on surface)
assert!(d.abs() < 1e-6);

// Box: center at origin, half-extents (3, 2, 1)
let d = sdf_box(Vector3::new(4.0, 0.0, 0.0), Vector3::zeros(), Vector3::new(3.0, 2.0, 1.0));
// 4 - 3 = 1 unit outside in X
assert!((d - 1.0).abs() < 1e-6);

// Cylinder: base at origin, along Y-axis, radius 2, height 10
let d = sdf_cylinder(
    Vector3::new(1.0, 5.0, 0.0),
    Vector3::zeros(),
    Vector3::y(),
    2.0, 10.0,
);
assert!(d < 0.0); // inside
```

### Object layer (`shape` module)

Structs that implement the `Sdf` trait, enabling composition via generics:

```rust
use crusst::shape::{Sphere, Box3, Union, Difference, Sdf};
use nalgebra::Vector3;

let sphere = Sphere::new(Vector3::zeros(), 5.0);
let box_shape = Box3::new(Vector3::zeros(), Vector3::new(3.0, 3.0, 3.0));

// Union: the shape is inside if you're inside EITHER primitive
let combined = Union::new(sphere, box_shape);
assert!(combined.evaluate(Vector3::new(4.5, 0.0, 0.0)) < 0.0); // inside sphere
assert!(combined.evaluate(Vector3::new(0.0, 2.5, 2.5)) < 0.0); // inside box
```

All shape structs implement `Send + Sync`, so they can be safely shared across threads.

**Available primitives:** `Sphere`, `Box3`, `CappedCone`, `Cylinder`, `Torus`, `RoundedBox`, `Capsule`, `Ellipsoid`, `RoundedCylinder`, `HalfSpace`.

---

## 3. Constructive Solid Geometry

CSG operations combine two solids into a new one. Crusst provides three families:

### Sharp CSG (exact boundaries)

```rust
use crusst::csg::{union, intersection, difference};

let d1 = -2.0; // inside solid A
let d2 =  1.0; // outside solid B

let u = union(d1, d2);        // min(-2, 1) = -2 (inside A)
let i = intersection(d1, d2); // max(-2, 1) = 1  (outside B, so outside A∩B)
let d = difference(d1, d2);   // max(-2, -1) = -1 (inside A and outside B)
```

### Smooth CSG (blended boundaries)

Smooth operations blend the transition between two solids using a polynomial kernel. The smoothing radius `k` controls the extent of the blend:

```rust
use crusst::csg::smooth_union;

// Two overlapping spheres
let d1 = -0.5; // slightly inside A
let d2 = -0.3; // slightly inside B
let blended = smooth_union(d1, d2, 1.0); // k=1.0

// The blended value is less than min(d1, d2), creating a smooth fill
assert!(blended < d1.min(d2));
```

At `k = 0`, smooth operations degenerate to their sharp counterparts. As `k` increases, the blend zone widens and the transition becomes more gradual.

### Round and Chamfer CSG

These apply geometric blend profiles at the CSG boundary:

```rust
use crusst::csg::{round_union, chamfer_union};

let d1 = 0.1;
let d2 = 0.1;

let round = round_union(d1, d2, 1.0);   // inscribed-circle fillet
let chamfer = chamfer_union(d1, d2, 1.0); // 45° flat bevel
```

Round CSG uses the inscribed-circle formula: the fillet profile is an exact circular arc tangent to both surfaces.

---

## 4. The Builder API

The `Shape` builder is the recommended API for most use cases. It wraps an `Arc<SdfNode>` and provides a fluent, chainable interface:

```rust
use crusst::builder::Shape;

// A pipe tee: main pipe + branch pipe, smoothly joined
let main_pipe = Shape::cylinder(3.0, 20.0);
let branch = Shape::cylinder(2.0, 12.0)
    .rotate_x(std::f64::consts::FRAC_PI_2)
    .translate(0.0, 5.0, 0.0);

let pipe_tee = main_pipe
    .smooth_union(branch, 1.5) // smooth blend at junction
    .shell(0.5);               // hollow with 0.5-unit wall thickness
```

**Key properties of `Shape`:**
- **Immutable** — every operation returns a new `Shape`, the original is unchanged
- **Cheap to clone** — internally `Arc`-shared, so `clone()` is just a reference count increment
- **Composable** — operations chain naturally because each returns `Shape`

### Querying geometry

```rust
let shape = Shape::sphere(5.0);

// Point queries
let d = shape.distance([3.0, 0.0, 0.0].into());
let inside = shape.contains([0.0, 0.0, 0.0].into());

// Feature queries (for fillet targeting)
let faces = shape.faces();   // Vec<FaceInfo>
let edges = shape.edges();   // Vec<EdgeInfo>

// Bounding box
let bbox = shape.bounding_box(); // BBox3 { min, max }
```

---

## 5. Spatial Transforms

Transforms modify the coordinate system of a shape without changing its SDF definition:

```rust
use crusst::builder::Shape;

let shape = Shape::box3(5.0, 2.0, 3.0);

// Translation: move 10 units along X
let moved = shape.translate(10.0, 0.0, 0.0);

// Rotation: 45° around the Z axis (radians)
let rotated = shape.rotate_z(std::f64::consts::FRAC_PI_4);

// Scale: double the size uniformly
let scaled = shape.scale(2.0);

// Mirror: reflect across the YZ plane (X → -X)
let mirrored = shape.mirror_x();

// Shell: hollow out with wall thickness 0.3
let hollow = shape.shell(0.3);

// Round: offset the surface inward by 0.2 (shrinks the shape)
let rounded = shape.round(0.2);
```

**How transforms work internally:**

Transforms in SDF space operate by transforming the query point into the shape's local coordinate system, then evaluating the original SDF. For example, `Translate` subtracts the offset from the query point before evaluating. `Rotate` applies the inverse rotation. `Scale` divides the query point by the scale factor and multiplies the result.

This "inverse transform" approach means transforms are exact — no approximation is introduced.

---

## 6. Mesh Extraction

Crusst converts SDFs to triangle meshes using **adaptive dual contouring** (DC). DC is preferred over marching cubes because it preserves sharp features and produces higher-quality meshes with fewer triangles.

### Basic usage

```rust
use crusst::builder::Shape;
use crusst::types::MeshSettings;

let shape = Shape::box3(5.0, 5.0, 5.0);
let mesh = shape.mesh(MeshSettings::default());

println!("Vertices: {}", mesh.vertices.len());
println!("Triangles: {}", mesh.indices.len() / 3);
```

### Mesh settings

`MeshSettings` controls the resolution and quality of mesh extraction:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_depth` | 8 | Maximum octree subdivision depth. Each level doubles resolution. Depth 7 → ~128 cells/axis. |
| `min_depth` | 6 | Minimum depth before interval pruning can skip cells. |
| `edge_tolerance` | 1e-6 | Bisection tolerance for finding surface crossings on cell edges. |

```rust
// High-resolution mesh
let hq = MeshSettings {
    max_depth: 9,   // ~512 cells/axis
    min_depth: 7,
    edge_tolerance: 1e-7,
};

// Fast preview
let preview = MeshSettings {
    max_depth: 6,   // ~64 cells/axis
    min_depth: 4,
    edge_tolerance: 1e-5,
};
```

### Using the compatibility wrapper

The `extract_mesh` function accepts any `&dyn Sdf` and a resolution count instead of mesh settings:

```rust
use crusst::mesh::extract_mesh;
use crusst::shape::{Sphere, Sdf};
use nalgebra::Vector3;

let sphere = Sphere::new(Vector3::zeros(), 5.0);
let mesh = extract_mesh(
    &sphere,
    Vector3::new(-6.0, -6.0, -6.0),
    Vector3::new(6.0, 6.0, 6.0),
    128, // resolution → converted to octree depth 7
);
```

### The TriangleMesh struct

```rust
pub struct TriangleMesh {
    pub vertices: Vec<Vector3<f64>>,  // vertex positions
    pub normals: Vec<Vector3<f64>>,   // per-vertex normals (from SDF gradient)
    pub indices: Vec<u32>,            // triangle indices (length = 3 × num_triangles)
}
```

The mesh uses **indexed geometry** with shared vertices — each unique vertex position appears once, and triangles reference vertices by index. This produces watertight topology suitable for 3D printing and CAD interchange.

`TriangleMesh::to_binary()` serializes the mesh to a compact binary format for WebSocket/GPU transmission.

---

## 7. Export Formats

### Wavefront OBJ (indexed)

```rust
shape.export_obj("output.obj").unwrap();
```

Writes `v` (vertex), `vn` (normal), and `f v//vn` (face) lines. The indexed format preserves shared vertex topology, which is important for watertight meshes.

### Binary STL

```rust
shape.export_stl("output.stl").unwrap();
```

Standard 80-byte header + triangle array format. Compatible with all major CAD and 3D printing software.

### STEP AP203 (smart tiered)

```rust
shape.export_step("output.step").unwrap();
```

The STEP exporter uses a smart classification system:
- **Exact BRep** for primitives with known analytical surfaces (sphere → `SPHERICAL_SURFACE`, box → 6 `PLANE`s, cylinder → `CYLINDRICAL_SURFACE`)
- **Tessellated** for everything else (booleans, fillets, smooth ops)

Rigid transforms (translate, rotate, scale) are unwrapped and applied to the underlying exact geometry when possible.

```rust
use crusst::step_export::classify;
use crusst::step_export::ExportTier;

let tier = classify(shape.node());
match tier {
    ExportTier::Exact => println!("Analytical BRep export"),
    ExportTier::Tessellated => println!("Mesh-based fallback"),
}
```

---

## 8. Fillet & Chamfer Operations

The targeted fillet/chamfer system is Crusst's most advanced feature. It applies blend profiles to specific edges of a shape, controlling geometric continuity from G0 (positional) through G3 (curvature-rate).

### Concepts

A **fillet** rounds a sharp edge by replacing the intersection of two faces with a smooth blend curve. The blend profile determines the cross-sectional shape of the rounding:

- **G0 (Chamfer):** Positional continuity only — a flat bevel that meets adjacent faces at a visible angle
- **G1 (Bezier):** Tangent continuity — the blend meets faces tangentially, eliminating visible edges
- **G2 (Circular):** Curvature continuity — the blend has constant curvature, matching the curvature of adjacent flat faces (zero)
- **G3 (Squircle/Parabolic/Hyperbolic):** Curvature-rate continuity — the rate of curvature change is also continuous, producing the smoothest possible transition

### Basic usage

```rust
use crusst::builder::Shape;
use crusst::blend;
use crusst::feature::ft;

let box_shape = Shape::box3(5.0, 5.0, 5.0);

// G2 circular fillet on all edges, radius 1.0
let filleted = box_shape.fillet(blend::g2(1.0), vec![ft(0, 0).all_edges()]);
```

### Blend profiles

```rust
use crusst::blend;

// Chamfers (G0)
let chamfer = blend::equal_chamfer(1.0);          // symmetric 45° bevel
let asym = blend::two_dist_chamfer(0.5, 1.5);     // asymmetric bevel
let angled = blend::angle_chamfer(1.0, 0.6);      // custom angle

// Tangent-continuous (G1)
let g1 = blend::g1(1.0);                          // cubic Bezier arc

// Curvature-continuous (G2)
let g2 = blend::g2(1.0);                          // exact circular arc
let chord = blend::chord(1.4);                    // arc by chord length

// Curvature-rate-continuous (G3)
let g3 = blend::g3(1.0);                          // squircle (L4 norm)
let para = blend::parabolic(1.0);                 // sin-curvature (L3 norm)
let hyper = blend::hyperbolic(1.0, 0.5);          // sinh-curvature (L6 norm)
let cyclo = blend::cycloidal(1.0);                // cycloid approx (L1.5 norm)
```

### Targeting specific edges

The `ft()` function creates a `FeatureTarget` builder. A `Box3` has 6 faces and 12 edges:

```
Faces: 0=+X, 1=-X, 2=+Y, 3=-Y, 4=+Z, 5=-Z
Edges: 0-3 (bottom), 4-7 (top), 8-11 (vertical pillars)
```

```rust
use crusst::feature::ft;

// All edges on body 0 of component 0
let all = ft(0, 0).all_edges();

// Only the top 4 edges
let top = ft(0, 0).edges(&[4, 5, 6, 7]);

// Only the vertical edges
let vertical = ft(0, 0).edges(&[8, 9, 10, 11]);
```

### How the math works

For G2 and above, the fillet evaluation uses a **multi-face Lp norm**. Given face distances `d₁, d₂, ...` to the faces adjacent to the targeted edges:

```
uᵢ = clamp(dᵢ + r, 0)       (hard clamp for G2)
uᵢ = softplus(dᵢ + r)        (soft clamp for G3+)
blend = (Σ uᵢᵖ)^(1/p) − r
```

The Lp exponent `p` determines the profile shape:
- `p = 2` → circle (G2)
- `p = 1.5` → cycloid (between chamfer and circle)
- `p = 3` → parabolic
- `p = 4` → squircle (G3)
- `p = 6` → hyperbolic

The softplus function `ln(1 + exp(k·x)) / k` (where `k = 3/r`) replaces the hard `max(0, x)` for higher-order profiles, causing the face surfaces themselves to curve smoothly into the blend zone — this is what produces G3 continuity.

Per-edge Newton iteration is used for profiles that don't map to Lp norms (G1 Bezier, chamfer variants). The Newton solver finds the closest point on the parametric blend curve and computes a signed distance via cross-product orientation.

---

## 9. Round & Chamfer CSG

While targeted fillets apply to specific edges of a single shape, round/chamfer CSG applies a blend at the boundary between two shapes during a boolean operation:

```rust
use crusst::builder::Shape;

let a = Shape::sphere(5.0).translate(-3.0, 0.0, 0.0);
let b = Shape::sphere(5.0).translate(3.0, 0.0, 0.0);

// Round union: fillet at the junction
let filleted_union = a.clone().round_union(b.clone(), 1.0);

// Round subtraction: fillet at the cut boundary
let filleted_sub = a.clone().round_subtract(b.clone(), 1.0);

// Chamfer union: beveled junction
let beveled = a.chamfer_union(b, 1.0);
```

**When to use which:**
- **Targeted fillet/chamfer** — when you need control over specific edges and continuity class
- **Round/Chamfer CSG** — when you want a blend everywhere two shapes meet during a boolean

---

## 10. 2D Shapes, Revolve, and Extrude

Crusst supports 2D signed distance functions that can be lifted to 3D via revolution or extrusion.

### 2D Primitives

```rust
use crusst::shape::{Circle2d, Rect2d, Sdf2d};
use nalgebra::Vector2;

let circle = Circle2d::new(Vector2::zeros(), 5.0);
let rect = Rect2d::new(Vector2::zeros(), Vector2::new(3.0, 2.0));

assert!(circle.evaluate(Vector2::new(3.0, 0.0)) < 0.0); // inside
```

### Revolve

`Revolve` sweeps a 2D profile around the Y axis to create a 3D solid of revolution:

```rust
use crusst::shape::{Rect2d, Revolve, Sdf};
use nalgebra::{Vector2, Vector3};

// A 2D rectangle revolved around Y → a cylindrical ring
let profile = Rect2d::new(Vector2::new(5.0, 0.0), Vector2::new(1.0, 3.0));
let ring = Revolve::new(profile);

// Query point at (5, 0, 0) is on the revolution centerline of the profile
assert!(ring.evaluate(Vector3::new(5.0, 0.0, 0.0)) < 0.0);
```

### Extrude

`Extrude` extends a 2D profile along the Z axis:

```rust
use crusst::shape::{Circle2d, Extrude, Sdf};
use nalgebra::{Vector2, Vector3};

// A circle extruded ±5 units along Z → a cylinder
let profile = Circle2d::new(Vector2::zeros(), 3.0);
let cylinder = Extrude::new(profile, 5.0);

assert!(cylinder.evaluate(Vector3::new(0.0, 0.0, 4.0)) < 0.0); // inside
assert!(cylinder.evaluate(Vector3::new(0.0, 0.0, 6.0)) > 0.0); // outside
```

---

## 11. Voxelization

Voxelization samples the SDF on a regular 3D grid, producing a `VoxelGrid` of `f32` distance values:

```rust
use crusst::builder::Shape;

let shape = Shape::sphere(5.0);
let grid = shape.voxelize(0.5); // 0.5-unit voxel size

println!("Grid resolution: {:?}", grid.resolution); // e.g., [22, 22, 22]
println!("Origin: {:?}", grid.origin);

// Access the distance value at a specific voxel
let idx = grid.index_at(nalgebra::Vector3::zeros());
let distance = grid.data[idx];
assert!(distance < 0.0); // center of sphere is inside
```

Voxelization is parallelized with rayon. The grid stores `f32` values (not `f64`) to reduce memory usage — for a 128^3 grid, this saves 16 MB vs. `f64`.

**Use cases:**
- Volume rendering
- Physics simulation (collision grids)
- Machine learning training data
- Constructive voxel operations

---

## 12. Parametric Paths & Transport

The transport system generates SDFs by sweeping a circular cross-section along a parametric path. The four "orders" correspond to increasing degrees of freedom:

### Order 0 — Identity

The cross-section itself, unmodified:

```rust
use crusst::shape::Sphere;
use crusst::transport::order0;

let sphere = order0(Sphere::new(nalgebra::Vector3::zeros(), 10.0));
// Result is just the sphere
```

### Order 1 — Rigid sweep with scaling (Cone eigenform)

A circle is swept along a path while its radius changes according to a scale function:

```rust
use crusst::path::LinePath;
use crusst::transport::order1;
use nalgebra::Vector3;

// Cone: radius tapers from 12 to 0 along a 30-unit line
let path = LinePath::new(Vector3::zeros(), Vector3::new(0.0, 0.0, 30.0));
let cone = order1(path, 12.0, |t| 1.0 - t, 128);
```

### Order 2 — Frame transport (Helix eigenform)

A circle is swept along a curved path using Frenet frame transport:

```rust
use crusst::path::HelixPath;
use crusst::transport::order2;

let helix = HelixPath::new(15.0, 8.0, 3.0); // radius, pitch, turns
let tube = order2(helix, 2.5, 256);          // section radius 2.5, 256 samples
```

### Order 3 — Full morphing (Horn eigenform)

Scale and twist vary along the path:

```rust
use crusst::path::SpiralPath;
use crusst::transport::order3;

let spiral = SpiralPath::new(
    |t| 20.0 * (1.0 - t * 0.7), // radius shrinks from 20 to 6
    |t| 40.0 * t,                // height rises to 40
    2.0,                          // 2 turns
);
let horn = order3(spiral, 3.0, |t| 1.0 - t * 0.8, |_| 0.0, 256);
```

The `samples` parameter controls the number of discrete path evaluations. Higher values produce smoother surfaces but are more expensive. 128–256 works well for most cases.

---

## 13. Custom SDF Functions

You can create shapes from arbitrary closures using `FnSdf`:

```rust
use crusst::shape::{FnSdf, Sdf};
use nalgebra::Vector3;

// A gyroid — a triply periodic minimal surface
let gyroid = FnSdf::new(|p: Vector3<f64>| {
    let scale = 0.5;
    let x = p.x * scale;
    let y = p.y * scale;
    let z = p.z * scale;
    (x.sin() * y.cos() + y.sin() * z.cos() + z.sin() * x.cos()) / scale
});

assert!(gyroid.evaluate(Vector3::zeros()) < 1e-6);
```

`FnSdf` implements `Sdf` so the result can be composed, meshed, and exported like any built-in shape.

For 2D custom shapes, use `FnSdf2d`.

---

## 14. The DAG Layer

The `SdfNode` enum is the internal representation of all geometry in Crusst. While the `Shape` builder is the recommended user API, understanding the DAG is useful for advanced applications.

### Structure

Each `SdfNode` variant holds its parameters and child nodes (via `Arc`):

```rust
use crusst::dag::SdfNode;
use nalgebra::Vector3;

// Manually construct a DAG (usually done via Shape builder)
let sphere = SdfNode::Sphere {
    center: Vector3::zeros(),
    radius: 5.0,
};
let translated = SdfNode::Translate(
    std::sync::Arc::new(sphere),
    Vector3::new(10.0, 0.0, 0.0),
);
```

### Three evaluation modes

**Point evaluation** — the standard SDF query:

```rust
let d = node.evaluate(Vector3::new(15.0, 0.0, 0.0));
```

**Interval evaluation** — computes conservative bounds over a bounding box, used for octree pruning:

```rust
use crusst::types::BBox3;

let bbox = BBox3::new(
    Vector3::new(20.0, 20.0, 20.0),
    Vector3::new(30.0, 30.0, 30.0),
);
let interval = node.interval_evaluate(&bbox);
if interval.definitely_positive() {
    // Entire box is outside the shape — no surface here
}
```

**Gradient evaluation** — the SDF gradient `∇f` at a point, used for dual contouring vertex placement and normal computation:

```rust
let normal = node.gradient(Vector3::new(5.0, 0.0, 0.0));
// For a sphere: points radially outward, unit length
```

Gradients are computed analytically where possible (primitives, CSG via gradient blending, transforms via chain rule) and fall back to central differences for complex operations.

### Feature introspection

The DAG exposes geometric feature information for fillet/chamfer targeting:

```rust
let box_node = SdfNode::Box3 {
    center: Vector3::zeros(),
    half: Vector3::new(5.0, 5.0, 5.0),
};

let faces = box_node.face_info().unwrap(); // 6 faces
let edges = box_node.edge_info().unwrap(); // 12 edges

for face in &faces {
    println!("{}: normal={:?}", face.label, face.normal);
}
for edge in &edges {
    println!("{}: faces ({}, {})", edge.label, edge.face_a, edge.face_b);
}
```

---

## 15. Performance Considerations

### Octree depth vs. quality

The octree depth exponentially affects both resolution and computation time:

| `max_depth` | Cells/axis | Approximate triangles | Relative time |
|-------------|------------|----------------------|---------------|
| 5 | 32 | ~10K | 1x |
| 6 | 64 | ~50K | 4x |
| 7 | 128 | ~100–200K | 16x |
| 8 | 256 | ~400–600K | 64x |
| 9 | 512 | ~1–2M | 256x |

The adaptive octree only subdivides near the surface, so the actual cost scales with surface area, not volume.

### SDF evaluation cost

SDF evaluation is the dominant cost in meshing. Complex DAGs (deep boolean trees, many transforms) evaluate more slowly than simple primitives. Tips:

- **Keep DAGs shallow** — fewer operations between the root and each primitive means faster evaluation
- **Avoid redundant transforms** — compose transforms manually when possible
- **Use `smooth_union` judiciously** — smooth operations evaluate both children at every point (no pruning possible)

### Parallelism

- **Voxelization** is parallelized with rayon (automatic multi-core)
- **Mesh extraction** is currently single-threaded but the octree structure is designed for future parallelization
- **SDF evaluation** is inherently parallel — all `Sdf` implementations are `Send + Sync`

### Export format size

| Format | Size (typical sphere, depth 7) | Notes |
|--------|-------------------------------|-------|
| OBJ | ~15 MB | Human-readable, indexed vertices |
| STL | ~5 MB | Binary, per-triangle normals (no sharing) |
| STEP | ~1 KB (exact) / ~5 MB (tessellated) | Exact BRep for simple primitives |

For shapes that qualify for exact STEP export (sphere, box, cylinder + rigid transforms), the file size is orders of magnitude smaller than tessellated formats.
