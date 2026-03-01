# Crusst — B-Rep Geometry Kernel in Rust

A pure Rust **Boundary Representation (B-Rep)** geometry kernel. Crusst stores shapes as exact mathematical surfaces with explicit topology — vertices, edges, faces, shells, and solids — enabling precise CAD operations, adaptive tessellation, and lossless STEP export.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Core Concepts](#core-concepts)
- [API Overview](#api-overview)
  - [Builder API](#builder-api)
  - [Primitives](#primitives)
  - [Transforms](#transforms)
  - [Profiles & Sketch Operations](#profiles--sketch-operations)
  - [Tessellation](#tessellation)
  - [Export](#export)
  - [Topology Validation](#topology-validation)
- [Architecture](#architecture)
- [Interactive Viewer](#interactive-viewer)
- [Dependencies](#dependencies)
- [License](#license)

## Installation

```toml
[dependencies]
crusst = { git = "https://github.com/crussella0129/crusst.git", branch = "v1.4" }
```

All spatial types use `nalgebra::Vector3<f64>`. Enable the `mint` feature for interoperability with glam, cgmath, or ultraviolet:

```toml
crusst = { git = "https://github.com/crussella0129/crusst.git", branch = "v1.4", features = ["mint"] }
```

## Quick Start

### Create and mesh a shape

```rust
use crusst::builder::Shape;
use crusst::types::TessSettings;

let mesh = Shape::sphere(5.0)
    .translate(10.0, 0.0, 0.0)
    .auto_mesh();

println!("Triangles: {}", mesh.indices.len() / 3);
println!("Vertices: {}", mesh.vertices.len());
```

### Profile-based extrusion

```rust
use crusst::builder::Shape;
use crusst::profile::Profile;

let bracket = Shape::extrude(&Profile::rect(10.0, 4.0), 20.0)
    .translate(0.0, 0.0, -10.0);
```

### Export to multiple formats

```rust
use crusst::builder::Shape;
use crusst::types::TessSettings;
use std::fs::File;

let shape = Shape::cylinder(3.0, 10.0);
let settings = TessSettings::default();

// Tessellation-based formats
shape.write_stl(&settings, &mut File::create("part.stl").unwrap()).unwrap();
shape.write_obj(&settings, &mut File::create("part.obj").unwrap()).unwrap();
shape.write_3mf(&settings, &mut File::create("part.3mf").unwrap()).unwrap();

// Exact B-Rep geometry (no tessellation)
shape.write_step(&mut File::create("part.step").unwrap()).unwrap();
```

## Core Concepts

### Boundary Representation

A B-Rep kernel represents solids as a hierarchy of topological entities, each carrying exact mathematical geometry:

| Entity | Geometry | Role |
|--------|----------|------|
| **Vertex** | `Point3` | Corner point of the solid |
| **Edge** | `Curve3` (line, circle, ellipse, NURBS) | Boundary between two faces |
| **CoEdge** | `Curve2` (parametric trim curve) | Oriented half-edge with face association |
| **Wire** | Closed loop of CoEdges | Boundary loop of a face |
| **Face** | `Surface` (plane, cylinder, sphere, ...) | Geometric surface patch |
| **Shell** | Collection of faces | Closed surface envelope |
| **Solid** | Outer shell + optional inner shells | Watertight volume |

This is the same architecture used by STEP, IGES, and professional CAD kernels — shapes are stored as the actual mathematical equations, not approximations.

### Surfaces

Six analytic surface types plus full NURBS support:

| Surface | Parametric Domain | Description |
|---------|-------------------|-------------|
| `Plane` | `(u, v) → origin + u·e1 + v·e2` | Infinite flat surface |
| `Cylinder` | `(θ, h) → origin + R·cos(θ)·e1 + R·sin(θ)·e2 + h·axis` | Circular cylinder |
| `Cone` | `(θ, h) → apex + h·(axis + tan(α)·(cos(θ)·e1 + sin(θ)·e2))` | Conical surface |
| `Sphere` | `(θ, φ) → center + R·cos(φ)·cos(θ)·e1 + ...` | Sphere |
| `Torus` | `(θ, φ) → center + (R + r·cos(φ))·cos(θ)·e1 + ...` | Ring torus |
| `NurbsSurface` | de Boor evaluation | Freeform surface |

All surfaces implement `evaluate(u, v) → Point3`, `normal(u, v) → Vector3`, `derivative_u/v`, and `min_curvature_radius()`.

### NURBS Math

The `nurbs` module provides the mathematical foundation:
- **Cox-de Boor** basis function evaluation
- **de Boor algorithm** for curve and surface point evaluation
- **Knot operations** — span finding, knot insertion, multiplicity queries
- **`NurbsCurve3`** and **`NurbsSurface`** structs with weighted control points

## API Overview

### Builder API

The `Shape` struct wraps a `(TopoStore, SolidId)` pair and provides a fluent interface for constructing, transforming, and exporting geometry.

```rust
use crusst::builder::Shape;
use crusst::types::TessSettings;

let bracket = Shape::box3(10.0, 2.0, 6.0)
    .translate(0.0, 5.0, 0.0)
    .rotate_z(std::f64::consts::FRAC_PI_4);

let mesh = bracket.auto_mesh();
bracket.write_step(&mut std::fs::File::create("bracket.step").unwrap()).unwrap();
```

### Primitives

| Method | Description | Faces |
|--------|-------------|-------|
| `Shape::box3(hx, hy, hz)` | Axis-aligned box (half-extents) | 6 planar |
| `Shape::sphere(radius)` | Sphere at origin | 1 spherical |
| `Shape::cylinder(radius, height)` | Z-axis cylinder | 3 (barrel + 2 caps) |
| `Shape::cone(r1, r2, height)` | Cone frustum on Z axis | 3 (conical + 2 caps) |
| `Shape::torus(major_r, minor_r)` | Torus on Z axis | 1 toroidal |
| `Shape::wedge(dx, dy, dz)` | Triangular prism | 5 planar |
| `Shape::capsule(radius, height)` | Hemisphere-capped cylinder | 3 (barrel + 2 spherical caps) |

### Transforms

All transforms modify geometry in place (no lazy DAG — vertex positions and surface equations are updated directly):

```rust
let shape = Shape::box3(5.0, 2.0, 3.0)
    .translate(10.0, 0.0, 0.0)
    .rotate_z(std::f64::consts::FRAC_PI_4)
    .scale(2.0)
    .mirror_x();
```

| Method | Effect |
|--------|--------|
| `translate(x, y, z)` | Rigid translation |
| `rotate_x/y/z(angle)` | Rotation around axis (radians) |
| `scale(factor)` | Uniform scaling |
| `mirror_x/y/z()` | Reflection across a coordinate plane |

### Profiles & Sketch Operations

2D profiles define cross-sections for extrude and revolve operations:

```rust
use crusst::builder::Shape;
use crusst::profile::Profile;
use crusst::math::Point2;

// Built-in profiles
let rect = Profile::rect(10.0, 4.0);
let circle = Profile::circle(3.0);

// Custom polygon
let triangle = Profile::polygon(&[
    Point2::new(0.0, 0.0),
    Point2::new(4.0, 0.0),
    Point2::new(2.0, 3.0),
]);

// Shape-generating operations
let prism = Shape::extrude(&rect, 20.0);
let donut = Shape::revolve(&circle, std::f64::consts::TAU);
```

Profiles support line and arc segments. `Profile::evaluate(t)` returns points along the closed wire at normalized parameter `t ∈ [0, 1]`.

### Tessellation

Crusst uses **adaptive parametric tessellation** — each face is sampled in its `(u, v)` parameter space, and triangles are refined where the surface deviates from the linear approximation.

```rust
use crusst::builder::Shape;
use crusst::types::TessSettings;

let shape = Shape::torus(10.0, 3.0);

// Auto-tuned settings (recommended)
let mesh = shape.auto_mesh();

// Manual settings
let mesh = shape.mesh(&TessSettings {
    chord_tolerance: 0.01,  // max deviation from true surface
    max_edge_length: 5.0,   // max triangle edge length
    min_subdivisions: 32,   // initial grid density per face
});
```

**Auto-tuning** computes tessellation parameters from the shape's actual surface curvature, not just its bounding box:

- **Chord tolerance** = 0.3% of the minimum curvature radius across all faces
- **Max edge length** = 10% of the bounding diagonal (loose safety net)
- **Min subdivisions** = 32 for curved shapes, 4 for all-planar shapes

This means a tall thin cylinder (R=2, H=50) gets the same quality as a squat one (R=2, H=2) — the tessellation tracks surface math, not spatial extent.

**Pipeline:**
1. Initial `n × n` grid in `(u, v)` space per face
2. Up to 3 adaptive refinement passes — triangles split where chord deviation exceeds tolerance or edges exceed max length
3. Analytical normals from `surface.normal(u, v)` — no gradient approximation
4. Planar faces use ear-clipping (no over-tessellation of flat geometry)

### Export

Four export formats, two strategies:

| Format | Method | Strategy |
|--------|--------|----------|
| **STEP AP203** | `write_step(writer)` | Exact B-Rep geometry — surfaces, curves, topology map directly to STEP entities |
| **STL** | `write_stl(settings, writer)` | Binary triangle mesh from tessellation |
| **OBJ** | `write_obj(settings, writer)` | Wavefront text format with indexed vertices and normals |
| **3MF** | `write_3mf(settings, writer)` | XML-based mesh with metadata |

STEP export writes exact mathematical geometry — `Plane` becomes `PLANE`, `Sphere` becomes `SPHERICAL_SURFACE`, `Face` becomes `ADVANCED_FACE`, `Solid` becomes `MANIFOLD_SOLID_BREP`. No tessellation required, no precision loss.

### Topology Validation

Every primitive is validated against the genus-aware Euler formula:

```
V - E + F = 2 - 2g
```

where `g` is the genus of the solid (0 for sphere-like, 1 for torus).

```rust
let shape = Shape::torus(10.0, 3.0);
let result = shape.validate();
assert!(result.is_valid());
```

Validation checks:
- Euler characteristic matches genus
- Every wire is closed (CoEdge chain loops back)
- Every edge has exactly 2 CoEdges (manifold guarantee)
- Vertex/edge/face counts are consistent

## Architecture

```
crusst/src/
├── math/              Linear algebra helpers (Point2, Point3, Vector2, Vector3)
├── nurbs/             NURBS math (Cox-de Boor, de Boor, knot operations)
│   ├── basis.rs       Basis function evaluation
│   └── knot.rs        Knot span finding, insertion, multiplicity
├── curve/             Curve3 enum (Line, Circle, Ellipse, NurbsCurve)
│                      Curve2 enum for parametric trim curves (PCurves)
├── surface/           Surface enum (Plane, Cylinder, Cone, Sphere, Torus, NURBS)
│                      evaluate(), normal(), derivatives, min_curvature_radius()
├── topo/              Arena-based topology (no Arc/Rc — typed indices into Vec<T>)
│   ├── store.rs       TopoStore arena with add/get/mut accessors
│   ├── types.rs       Vertex, Edge, CoEdge, Wire, Face, Shell, Solid
│   └── validate.rs    Genus-aware Euler formula validation
├── primitive/         Shape constructors → SolidId in a TopoStore
│   ├── box3.rs        6 planar faces, 12 edges, 8 vertices
│   ├── sphere.rs      1 spherical face with pole trim
│   ├── cylinder.rs    Barrel + 2 circular caps
│   ├── cone.rs        Conical frustum + 2 caps
│   ├── torus.rs       1 toroidal face (genus 1)
│   ├── wedge.rs       Triangular prism (5 faces)
│   └── capsule.rs     Hemisphere-capped cylinder
├── profile/           2D sketch profiles (line/arc segments)
├── tessellate/        Adaptive parametric tessellation
│   ├── mod.rs         Face dispatch (planar → ear-clip, curved → adaptive)
│   └── adaptive.rs    Grid + refinement passes with chord deviation control
├── builder.rs         Shape: fluent API over (TopoStore, SolidId)
├── types.rs           TriangleMesh, BBox3, TessSettings (curvature-based auto-tuning)
├── export/            Multi-format export
│   ├── step.rs        STEP AP203 (exact B-Rep entities)
│   ├── stl.rs         Binary STL
│   ├── obj.rs         Wavefront OBJ with normals
│   └── threemf.rs     3MF XML model
├── ops/               Operations (fillet, chamfer, shell — stubs)
└── boolean/           Boolean trait + stub (deferred — hardest part of any B-Rep kernel)
```

**Key design decisions:**

1. **Arena-based topology** — `TopoStore` holds `Vec<T>` per entity type, accessed via typed indices (`VertexId(usize)`, `EdgeId(usize)`, etc.). This avoids the `Arc`/`Rc` reference cycles inherent in topology's cyclic graph.

2. **Geometry/topology separation** — Geometric carriers (curves, surfaces) are value types stored on topology entities. Topology handles adjacency and orientation. This maps directly to STEP entities.

3. **Direct transforms** — Translations, rotations, and scales modify vertex positions and surface equations in place. No lazy evaluation DAG — what you see is what's stored.

4. **Curvature-driven tessellation** — Quality is controlled by surface math (curvature radius), not bounding box size. A cylinder always looks smooth regardless of how tall it is.

## Interactive Viewer

The included viewer serves a gallery of shapes in browser via a local HTTP server:

```bash
cargo run --example viewer
# Open http://localhost:8080 in your browser
```

Features:
- Real-time 3D preview of all primitives with Three.js
- Auto/manual tessellation toggle with live parameter controls
- Shape info panel (vertices, triangles, bounding diagonal, curvature radius)
- One-click STL, OBJ, STEP, and 3MF export per shape

## Dependencies

| Crate | Version | Purpose |
|-------|---------|---------|
| [`nalgebra`](https://nalgebra.org) | 0.33 | Linear algebra (`Vector3`, `Point3`, `Rotation3`) |
| [`rayon`](https://docs.rs/rayon) | 1.10 | Parallel face tessellation |
| [`approx`](https://docs.rs/approx) | 0.5 | Floating-point assertions (dev only) |
| [`tiny_http`](https://docs.rs/tiny_http) | 0.12 | HTTP server for the viewer example (dev only) |

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT License ([LICENSE-MIT](LICENSE-MIT))

at your option.
