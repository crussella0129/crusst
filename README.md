# Crusst — SDF Geometry Kernel in Rust

A standalone **Signed Distance Function (SDF)** geometry kernel written in pure Rust. Crusst provides composable shape primitives, constructive solid geometry, targeted fillet/chamfer blending with continuity control (G0–G3), adaptive dual contouring mesh extraction, parametric path sweeping, and multi-format export — everything needed to go from mathematical shape definition to watertight triangle mesh.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Core Concepts](#core-concepts)
- [API Overview](#api-overview)
  - [Builder API](#builder-api)
  - [Primitives](#primitives)
  - [CSG Operations](#csg-operations)
  - [Transforms](#transforms)
  - [Fillet & Chamfer System](#fillet--chamfer-system)
  - [Round & Chamfer CSG](#round--chamfer-csg)
  - [Meshing](#meshing)
  - [Export](#export)
  - [Voxelization](#voxelization)
  - [Parametric Paths & Transport](#parametric-paths--transport)
- [Architecture](#architecture)
- [Documentation](#documentation)
- [Dependencies](#dependencies)
- [License](#license)

## Installation

```toml
[dependencies]
crusst = { git = "https://github.com/crussella0129/crusst.git", branch = "v1.1" }
```

All spatial types use `nalgebra::Vector3<f64>`. Enable the `mint` feature for interoperability with glam, cgmath, or ultraviolet:

```toml
crusst = { git = "https://github.com/crussella0129/crusst.git", branch = "v1.1", features = ["mint"] }
```

## Quick Start

### Evaluate a shape

```rust
use crusst::builder::Shape;

let sphere = Shape::sphere(5.0);
assert!(sphere.distance([0.0, 0.0, 0.0].into()) < 0.0);       // inside
assert!(sphere.distance([5.0, 0.0, 0.0].into()).abs() < 1e-6); // on surface
assert!(sphere.distance([10.0, 0.0, 0.0].into()) > 0.0);       // outside
```

### Boolean operations

```rust
use crusst::builder::Shape;

let base = Shape::box3(5.0, 5.0, 5.0);
let hole = Shape::cylinder(2.0, 12.0).translate(0.0, 0.0, 0.0);
let part = base.subtract(hole);
```

### Fillet edges

```rust
use crusst::builder::Shape;
use crusst::blend;
use crusst::feature::ft;

let box_shape = Shape::box3(5.0, 5.0, 5.0);
let filleted = box_shape.fillet(blend::g2(1.0), vec![ft(0, 0).all_edges()]);
```

### Mesh and export

```rust
use crusst::builder::Shape;
use crusst::types::MeshSettings;

let shape = Shape::sphere(5.0);
let mesh = shape.mesh(MeshSettings::default());

shape.export_obj("sphere.obj").unwrap();
shape.export_stl("sphere.stl").unwrap();
shape.export_step("sphere.step").unwrap();
```

## Core Concepts

### Signed Distance Functions

A **signed distance function** (SDF) maps every point in 3D space to a scalar value representing the shortest distance to the nearest surface:

| Region | SDF value | Convention |
|--------|-----------|------------|
| Inside the solid | Negative | `f(p) < 0` |
| On the surface | Zero | `f(p) = 0` |
| Outside the solid | Positive | `f(p) > 0` |

SDFs compose algebraically: `union(A, B) = min(A, B)`, `intersection(A, B) = max(A, B)`, `difference(A, B) = max(A, -B)`. This makes complex geometry trivially constructible from simple primitives.

### The DAG

Internally, Crusst represents geometry as a **directed acyclic graph** (`SdfNode`) of operations. Each node in the DAG is one of:
- A **primitive** (sphere, box, cylinder, ...)
- A **CSG operation** (union, intersection, difference, smooth variants)
- A **transform** (translate, rotate, scale, mirror, shell)
- A **targeted blend** (fillet or chamfer with a specific profile and feature targets)

The DAG enables three evaluation modes:
1. **Point evaluation** — `f(p) → f64` for SDF queries
2. **Interval evaluation** — `f([x₀,x₁]×[y₀,y₁]×[z₀,z₁]) → [lo, hi]` for conservative octree pruning
3. **Gradient evaluation** — `∇f(p) → Vector3` for dual contouring vertex placement

## API Overview

### Builder API

The `Shape` struct is the primary user-facing API. It wraps an `Arc<SdfNode>` and provides a fluent interface for constructing, querying, and exporting geometry. Shapes are cheaply cloneable (shared ownership via `Arc`).

```rust
use crusst::builder::Shape;

// Build a mounting bracket
let base = Shape::box3(10.0, 2.0, 6.0);
let hole = Shape::cylinder(1.5, 6.0).translate(3.0, 0.0, 0.0);
let bracket = base
    .subtract(hole.clone())
    .subtract(hole.translate(-6.0, 0.0, 0.0))
    .smooth_union(Shape::box3(2.0, 4.0, 6.0).translate(-4.0, 3.0, 0.0), 0.5);

bracket.export_obj("bracket.obj").unwrap();
```

**Primitive constructors:**

| Method | Description |
|--------|-------------|
| `Shape::sphere(radius)` | Sphere at origin |
| `Shape::box3(hx, hy, hz)` | Axis-aligned box (half-extents) |
| `Shape::cylinder(radius, height)` | Y-axis cylinder |
| `Shape::torus(major, minor)` | Torus in the XZ plane |
| `Shape::capsule(a, b, radius)` | Sphere-swept line segment |
| `Shape::capped_cone(a, b, ra, rb)` | Cone frustum between endpoints |
| `Shape::rounded_box(hx, hy, hz, r)` | Box with uniform edge rounding |
| `Shape::ellipsoid(rx, ry, rz)` | Axis-aligned ellipsoid |
| `Shape::rounded_cylinder(r, rr, hh)` | Cylinder with rounded edges |
| `Shape::half_space(normal, d)` | Infinite half-space |

### Primitives

The `primitives` module provides low-level functional SDF distance functions that operate on raw `Vector3<f64>` values:

```rust
use crusst::primitives::{sdf_sphere, sdf_box, sdf_cylinder, sdf_torus};
use nalgebra::Vector3;

let d = sdf_sphere(Vector3::new(3.0, 0.0, 0.0), Vector3::zeros(), 1.0);
assert!((d - 2.0).abs() < 1e-6); // distance 2 from a radius-1 sphere
```

All nine primitives: `sdf_sphere`, `sdf_box`, `sdf_capped_cone`, `sdf_cylinder`, `sdf_torus`, `sdf_rounded_box`, `sdf_capsule`, `sdf_ellipsoid`, `sdf_rounded_cylinder`.

### CSG Operations

**Sharp CSG** — exact min/max operations:

```rust
use crusst::builder::Shape;

let a = Shape::sphere(5.0).translate(-2.0, 0.0, 0.0);
let b = Shape::sphere(5.0).translate(2.0, 0.0, 0.0);

let union = a.clone().union(b.clone());        // A ∪ B
let inter = a.clone().intersect(b.clone());    // A ∩ B
let diff  = a.subtract(b);                     // A \ B
```

**Smooth CSG** — polynomial blending with radius `k`:

```rust
let smooth = a.smooth_union(b, 1.0);  // blended with k=1.0
```

The smoothing radius `k` controls how far from the intersection the blend extends. Larger values create more gradual transitions.

### Transforms

```rust
use crusst::builder::Shape;
use std::f64::consts::PI;

let shape = Shape::box3(5.0, 2.0, 3.0)
    .translate(10.0, 0.0, 0.0)
    .rotate_z(PI / 4.0)        // 45° around Z
    .scale(2.0)                 // uniform scale
    .mirror_x()                 // reflect across YZ plane
    .shell(0.5)                 // hollow with wall thickness 0.5
    .round(0.2);                // offset the surface inward by 0.2
```

| Method | Effect |
|--------|--------|
| `translate(x, y, z)` | Rigid translation |
| `rotate_x/y/z(angle)` | Rotation around axis (radians) |
| `scale(factor)` | Uniform scaling |
| `mirror_x/y/z()` | Reflection across a coordinate plane |
| `shell(thickness)` | Hollow shell (onion-skin) |
| `round(radius)` | Offset rounding |

### Fillet & Chamfer System

Crusst provides a targeted fillet/chamfer system that applies blend profiles to specific edges of a shape. The system uses a **feature ID** mechanism to identify faces and edges on primitives, then applies blend operations via multi-face Lp norm generalization.

**10 blend profiles** spanning G0 through G3 continuity:

| Profile | Constructor | Continuity | Shape |
|---------|-------------|------------|-------|
| Equal Chamfer | `blend::equal_chamfer(d)` | G0 | 45° flat bevel |
| Two-Distance Chamfer | `blend::two_dist_chamfer(d1, d2)` | G0 | Asymmetric flat bevel |
| Angle Chamfer | `blend::angle_chamfer(d, θ)` | G0 | Angled flat bevel |
| G1 Bezier | `blend::g1(r)` | G1 | Cubic Bezier quarter-circle |
| G2 Circular | `blend::g2(r)` | G2 | Exact circular arc (L2 norm) |
| Chord | `blend::chord(c)` | G2 | Arc specified by chord length |
| G3 Squircle | `blend::g3(r)` | G3 | Superellipse (L4 norm) |
| Parabolic | `blend::parabolic(r)` | G3 | Sin-based curvature (L3 norm) |
| Hyperbolic | `blend::hyperbolic(r, a)` | G3 | Sinh-based curvature (L6 norm) |
| Cycloidal | `blend::cycloidal(r)` | G3 | Cycloid approximation (L1.5 norm) |

**Feature targeting** with the `ft()` builder:

```rust
use crusst::builder::Shape;
use crusst::blend;
use crusst::feature::ft;

let box_shape = Shape::box3(5.0, 5.0, 5.0);

// Fillet all edges with a G2 circular arc, radius 1.0
let filleted = box_shape.fillet(blend::g2(1.0), vec![ft(0, 0).all_edges()]);

// Chamfer only the top 4 edges (edges 4-7 on a Box3)
let chamfered = box_shape.chamfer(
    blend::equal_chamfer(0.8),
    vec![ft(0, 0).edges(&[4, 5, 6, 7])],
);
```

**How it works:** The fillet evaluation collects the face indices adjacent to each targeted edge, then computes a multi-face Lp blend:

```
blend = (Σ clamp(dᵢ + r)^p)^(1/p) − r
```

where `dᵢ` are signed distances to each face, `r` is the blend radius, and `p` is the Lp exponent (2 for G2/circular, 4 for G3, etc.). Higher-order profiles (G3, parabolic, hyperbolic) use softplus clamping instead of hard `max(0, ·)` to curve the face surfaces smoothly into the blend zone.

### Round & Chamfer CSG

Global round and chamfer operations apply a blend at all CSG boundaries:

```rust
use crusst::builder::Shape;

let a = Shape::sphere(5.0).translate(-3.0, 0.0, 0.0);
let b = Shape::sphere(5.0).translate(3.0, 0.0, 0.0);

let round_u = a.clone().round_union(b.clone(), 1.0);     // filleted union
let round_i = a.clone().round_intersect(b.clone(), 1.0);  // filleted intersection
let round_d = a.clone().round_subtract(b.clone(), 1.0);   // filleted subtraction

let chamfer_u = a.chamfer_union(b, 1.0);                  // chamfered union
```

These use the inscribed-circle (round) or 45° flat (chamfer) blend formulas from the `csg` module.

### Meshing

Crusst extracts triangle meshes via **adaptive dual contouring** — a feature-preserving isosurface algorithm that places vertices at optimal positions using a Quadratic Error Function (QEF) solver.

```rust
use crusst::builder::Shape;
use crusst::types::MeshSettings;

let shape = Shape::sphere(5.0);

// Default settings: max_depth=8, min_depth=6, edge_tolerance=1e-6
let mesh = shape.mesh(MeshSettings::default());

// Custom settings for higher resolution
let mesh_hq = shape.mesh(MeshSettings {
    max_depth: 9,
    min_depth: 7,
    edge_tolerance: 1e-7,
});

println!("Triangles: {}", mesh.indices.len() / 3);
println!("Vertices: {}", mesh.vertices.len());
```

The pipeline:
1. **Adaptive octree** — subdivides only near the surface, using interval arithmetic to prune cells that are entirely inside or outside
2. **Edge crossing detection** — bisects cell edges to find surface crossings within `edge_tolerance`
3. **QEF vertex placement** — solves a least-squares system from crossing positions and SDF gradients to place each vertex optimally (preserves sharp features)
4. **Face generation** — connects QEF vertices across shared sign-changing edges into triangles

The compatibility function `extract_mesh(sdf, bbox_min, bbox_max, resolution)` accepts any `&dyn Sdf` and converts the resolution count to an octree depth.

### Export

Three export formats:

```rust
shape.export_obj("output.obj").unwrap();   // Wavefront OBJ (indexed, shared vertices)
shape.export_stl("output.stl").unwrap();   // Binary STL
shape.export_step("output.step").unwrap(); // STEP AP203 (smart tiered)
```

**STEP export** uses a smart tiered system:
- **Tier 1 (Exact):** Spheres → `SPHERICAL_SURFACE`, boxes → 6 `PLANE`s, cylinders → `CYLINDRICAL_SURFACE`. Rigid transforms are unwrapped and applied analytically.
- **Tier 3 (Tessellated):** Everything else (booleans, smooth ops, fillets) is meshed and exported as tessellated `ADVANCED_FACE` geometry.

The `classify(node)` function walks the DAG to determine which tier applies.

### Voxelization

Convert any shape to a regular 3D grid of SDF samples:

```rust
let grid = shape.voxelize(0.5); // 0.5-unit voxel size
println!("Resolution: {:?}", grid.resolution); // [nx, ny, nz]
println!("Origin: {:?}", grid.origin);
let value = grid.data[grid.index_at(some_world_point)];
```

Voxelization is parallelized with rayon. The resulting `VoxelGrid` stores `f32` distance values in row-major `[x][y][z]` order.

### Parametric Paths & Transport

The transport system sweeps a circular cross-section along a parametric path, producing SDF shapes of increasing complexity:

| Order | Eigenform | Description |
|-------|-----------|-------------|
| 0 | Sphere | Identity — the shape is the cross-section |
| 1 | Cone | Rigid sweep with scaling along a path |
| 2 | Helix tube | Frame transport along a curved path |
| 3 | Horn | Full morphing: scale + twist along path |

```rust
use crusst::path::{LinePath, HelixPath};
use crusst::transport::{order1, order2};
use nalgebra::Vector3;

// Cone: circle tapers to a point along a line
let path = LinePath::new(Vector3::zeros(), Vector3::new(0.0, 0.0, 30.0));
let cone = order1(path, 12.0, |t| 1.0 - t, 128);

// Helix tube: circle swept along a helical path
let helix = HelixPath::new(15.0, 8.0, 3.0);
let tube = order2(helix, 2.5, 256);
```

**Path types:**
- `LinePath` — straight line between two points
- `HelixPath` — constant-curvature helix (radius, pitch, turns)
- `SpiralPath` — closure-driven spiral with variable radius and height

All transport results implement `Sdf` and can be composed, meshed, or exported like any other shape.

## Architecture

```
crusst
├── primitives          9 functional SDF distance functions
├── csg                 12 scalar CSG operations (sharp, smooth, round, chamfer)
├── shape               Sdf/Sdf2d traits + generic composable shape structs
│   ├── primitives      10 primitive shapes (Sphere, Box3, Cylinder, ...)
│   ├── csg             6 boolean/smooth operations (generic over Sdf)
│   ├── transforms      6 spatial transforms (Translate, Rotate, Scale, Mirror, Shell, Round)
│   └── 2d              2D primitives + Revolve/Extrude lifts
├── dag                 SdfNode expression DAG
│   ├── evaluate()      Point SDF evaluation
│   ├── interval_evaluate()   Interval arithmetic (octree pruning)
│   ├── gradient()      Analytical + central-difference gradients
│   └── face/edge_info  Feature ID introspection for fillet targeting
├── builder             Shape: fluent API wrapping Arc<SdfNode>
├── blend               10 blend profiles (G0–G3) + Newton iteration solvers
├── feature             Feature ID types (FaceInfo, EdgeInfo, FeatureTarget)
├── octree              Adaptive octree with interval pruning
├── qef                 QEF solver (SVD-based, mass-point regularized)
├── dual_contouring     Mesh extraction: octree → edge crossings → QEF → triangles
├── mesh                Compatibility layer (resolution → depth conversion)
├── voxel               Regular SDF voxel grid (f32, rayon-parallel)
├── path                Parametric paths (Line, Helix, Spiral)
├── frame               Frenet frame computation
├── transport           Sweep operations (Orders 0–3 eigenforms)
├── export              Binary STL writer
├── obj_export          Wavefront OBJ writer (indexed)
└── step_export         STEP AP203 writer (exact BRep + tessellated fallback)
```

The design separates **distance evaluation** (primitives, shapes, DAG) from **mesh generation** (octree, dual contouring, QEF) and **file I/O** (export). The SDF layer can be used independently for ray marching, collision detection, or any other distance-field application.

For detailed technical documentation, see:
- **[Tutorial](docs/tutorial.md)** — step-by-step guide from first principles
- **[Architecture Guide](docs/architecture.md)** — deep dive into the meshing pipeline, blend math, and DAG internals

## Dependencies

| Crate | Version | Purpose |
|-------|---------|---------|
| [`nalgebra`](https://nalgebra.org) | 0.33 | Linear algebra (`Vector3`, `Matrix3`, `Rotation3`, SVD) |
| [`rayon`](https://docs.rs/rayon) | 1.10 | Parallel voxelization |
| [`approx`](https://docs.rs/approx) | 0.5 | Floating-point assertions (dev only) |
| [`tiny_http`](https://docs.rs/tiny_http) | 0.12 | HTTP server for the live viewer example (dev only) |

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT License ([LICENSE-MIT](LICENSE-MIT))

at your option.
