# Crusst — Rust Based Geometry Kernel

A standalone **Signed Distance Function (SDF)** geometry kernel written in pure Rust. Crusst provides composable shape primitives, boolean operations, parametric path sweeping, parallel mesh extraction, and STL export — everything needed to go from mathematical shape definition to triangle mesh.

## Installation

```toml
[dependencies]
crusst = { git = "https://github.com/crussella0129/crusst.git" }
```

Crusst depends on [`nalgebra`](https://nalgebra.org) for linear algebra. All spatial types use `Vector3<f64>`.

## Quick Start

```rust
use crusst::shape::{Sphere, Sdf};
use crusst::mesh::extract_mesh;
use crusst::export::write_stl;
use nalgebra::Vector3;

// Define a sphere with radius 5 centered at the origin
let sphere = Sphere::new(Vector3::zeros(), 5.0);

// Query the signed distance at any point
assert!(sphere.evaluate(Vector3::zeros()) < 0.0);       // inside
assert!(sphere.evaluate(Vector3::new(5.0, 0.0, 0.0)).abs() < 1e-6); // on surface
assert!(sphere.evaluate(Vector3::new(10.0, 0.0, 0.0)) > 0.0);      // outside

// Extract a triangle mesh via marching cubes
let mesh = extract_mesh(
    &sphere,
    Vector3::new(-6.0, -6.0, -6.0), // bounding box min
    Vector3::new(6.0, 6.0, 6.0),    // bounding box max
    64,                               // grid resolution
);

// Export to binary STL
write_stl(&mesh, "sphere.stl").unwrap();
```

## Modules

### `primitives` — Functional SDF Primitives

Low-level distance functions that take a query point and shape parameters, returning a signed distance value. Negative inside, zero on surface, positive outside.

| Function | Description |
|---|---|
| `sdf_sphere(point, center, radius)` | Distance to a sphere |
| `sdf_box(point, center, half_extents)` | Distance to an axis-aligned box |
| `sdf_capped_cone(point, a, b, ra, rb)` | Distance to a cone frustum between endpoints `a` and `b` with radii `ra` and `rb` |
| `sdf_cylinder(point, base, axis, radius, height)` | Distance to a cylinder |

```rust
use crusst::primitives::sdf_sphere;
use nalgebra::Vector3;

let d = sdf_sphere(Vector3::new(3.0, 0.0, 0.0), Vector3::zeros(), 1.0);
assert!((d - 2.0).abs() < 1e-6); // 3 units from center, radius 1 → distance 2
```

### `csg` — Constructive Solid Geometry Operations

Boolean operations on SDF distance values. Sharp variants use min/max; smooth variants blend with a smoothing radius `k`.

| Function | Operation |
|---|---|
| `union(d1, d2)` | A ∪ B — minimum of two distances |
| `intersection(d1, d2)` | A ∩ B — maximum of two distances |
| `difference(d1, d2)` | A \ B — subtracts B from A |
| `smooth_union(d1, d2, k)` | Blended union with radius `k` |
| `smooth_intersection(d1, d2, k)` | Blended intersection with radius `k` |
| `smooth_difference(d1, d2, k)` | Blended difference with radius `k` |

```rust
use crusst::primitives::sdf_sphere;
use crusst::csg::smooth_union;
use nalgebra::Vector3;

let d1 = sdf_sphere(Vector3::zeros(), Vector3::new(-1.0, 0.0, 0.0), 2.0);
let d2 = sdf_sphere(Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0), 2.0);
let blended = smooth_union(d1, d2, 0.5);
```

### `shape` — Composable Shape Objects

The `Sdf` trait provides an object-oriented interface for building shape trees. All shapes implement `Send + Sync` for safe parallel evaluation.

**Trait:**

```rust
pub trait Sdf: Send + Sync {
    fn evaluate(&self, point: Vector3<f64>) -> f64;
}
```

**Shapes:**

| Type | Description |
|---|---|
| `Sphere` | Sphere defined by center and radius |
| `Box3` | Axis-aligned box defined by center and half-extents |
| `CappedCone` | Cone frustum between two endpoints with independent radii |
| `HalfSpace` | Infinite half-space defined by a normal and offset |
| `Union<A, B>` | Boolean union of two shapes |
| `Intersection<A, B>` | Boolean intersection of two shapes |
| `Difference<A, B>` | Boolean difference (A minus B) |

Shapes compose naturally via generics:

```rust
use crusst::shape::{Sphere, Difference, Sdf};
use nalgebra::Vector3;

let outer = Sphere::new(Vector3::zeros(), 5.0);
let inner = Sphere::new(Vector3::zeros(), 3.0);
let shell = Difference::new(outer, inner);

// Point at radius 4 is between the two spheres → inside the shell
assert!(shell.evaluate(Vector3::new(4.0, 0.0, 0.0)) < 0.0);

// Point at the origin is inside the inner sphere → outside the shell
assert!(shell.evaluate(Vector3::zeros()) > 0.0);
```

### `path` — Parametric Paths

Paths are parametric curves over `t ∈ [0, 1]` used by the transport module to sweep shapes through space.

**Trait:**

```rust
pub trait Path: Send + Sync {
    fn point(&self, t: f64) -> Vector3<f64>;
    fn tangent(&self, t: f64) -> Vector3<f64>;
}
```

**Path types:**

| Type | Description | Parameters |
|---|---|---|
| `LinePath` | Straight line between two points | `start`, `end` |
| `HelixPath` | Circular helix rising along the Z axis | `radius`, `pitch`, `turns` |
| `SpiralPath` | Generalized spiral with dynamic radius and height | `radius_fn(t)`, `height_fn(t)`, `turns` |

```rust
use crusst::path::{HelixPath, Path};

let helix = HelixPath::new(
    10.0, // radius
    5.0,  // pitch (vertical rise per turn)
    3.0,  // number of turns
);
let start = helix.point(0.0);  // (10, 0, 0)
let end = helix.point(1.0);    // back to x=10, z = pitch * turns = 15
```

`SpiralPath` accepts closures for full control over the curve shape:

```rust
use crusst::path::SpiralPath;

// A horn-like spiral: radius shrinks from 20 to 6, height rises to 40
let spiral = SpiralPath::new(
    |t| 20.0 * (1.0 - t * 0.7), // radius_fn
    |t| 40.0 * t,                // height_fn
    2.0,                          // turns
);
```

### `frame` — Frenet Frames

`FrenetFrame` computes a local coordinate system (tangent, normal, binormal) along a path. Used internally by the transport module but available for custom sweep operations.

```rust
use crusst::frame::FrenetFrame;
use nalgebra::Vector3;

let frame = FrenetFrame::from_tangent(Vector3::new(0.0, 0.0, 1.0));
let rotation = frame.to_matrix(); // 3x3 rotation: local → world
```

### `transport` — Shape Sweeping (Eigenforms)

Transport operations sweep a circular cross-section along a path, producing increasingly complex geometry. Each "order" adds a degree of freedom:

| Order | Function | Description | Eigenform |
|---|---|---|---|
| 0 | `order0(section)` | Identity — the shape *is* the cross-section | Sphere |
| 1 | `order1(path, radius, scale_fn, samples)` | Rigid sweep with scaling along a path | Cone |
| 2 | `order2(path, radius, samples)` | Tube sweep with frame transport | Helix |
| 3 | `order3(path, radius, scale_fn, twist_fn, samples)` | Full morphing: scale + twist along path | Horn |

All transport functions return a `TransportShape` that implements `Sdf`, so the result can be meshed, composed with CSG, or nested in further operations.

```rust
use crusst::shape::Sphere;
use crusst::path::{LinePath, HelixPath};
use crusst::transport::{order0, order1, order2};
use crusst::mesh::extract_mesh;
use nalgebra::Vector3;

// Order 0: Just a sphere
let sphere = order0(Sphere::new(Vector3::zeros(), 10.0));

// Order 1: Cone — circle swept along a line, tapering to a point
let path = LinePath::new(Vector3::zeros(), Vector3::new(0.0, 0.0, 30.0));
let cone = order1(path, 12.0, |t| 1.0 - t, 128);

// Order 2: Helix tube — circle swept along a helical path
let helix = HelixPath::new(15.0, 8.0, 3.0);
let tube = order2(helix, 2.5, 256);
let mesh = extract_mesh(
    &tube,
    Vector3::new(-20.0, -20.0, -4.0),
    Vector3::new(20.0, 20.0, 28.0),
    100,
);
```

The `samples` parameter controls how many discrete points along the path are evaluated during the sweep. Higher values produce smoother surfaces but are more expensive. A value of 256 works well for most cases.

### `mesh` — Mesh Extraction

`extract_mesh` runs marching cubes over a bounding box at the given grid resolution, producing a `TriangleMesh` with per-vertex normals computed via central differences.

```rust
pub struct TriangleMesh {
    pub vertices: Vec<Vector3<f64>>,
    pub normals: Vec<Vector3<f64>>,
    pub indices: Vec<u32>,
}
```

`TriangleMesh::to_binary()` serializes the mesh to a compact binary format for WebSocket/GPU transmission:

```
[vertex_count: u32 LE]
[vertices: f32 × 3 × vertex_count LE]
[normals: f32 × 3 × vertex_count LE]
[index_count: u32 LE]
[indices: u32 × index_count LE]
```

Resolution determines the number of grid cells per axis. Higher resolution = more triangles and finer detail. Typical values:

| Resolution | Grid Points | Use Case |
|---|---|---|
| 32 | ~33K | Fast preview |
| 64 | ~262K | Development |
| 100 | ~1M | Production |
| 160 | ~4M | High quality |
| 200 | ~8M | Maximum detail |

Mesh extraction is parallelized with [`rayon`](https://docs.rs/rayon) for multi-core performance.

### `export` — STL Export

`write_stl` writes a binary STL file from a `TriangleMesh`:

```rust
use crusst::export::write_stl;

write_stl(&mesh, "output.stl").unwrap();
```

The output follows the standard 80-byte header + triangle array binary STL format, compatible with all major CAD and 3D printing software.

## Architecture

```
crusst
├── primitives     Functional SDF distance functions
├── csg            Boolean operations on distance values
├── shape          Sdf trait + composable shape structs
├── path           Parametric curves (Line, Helix, Spiral)
├── frame          Frenet frame computation
├── transport      Sweep operations (Orders 0–3)
├── mesh           Marching cubes extraction + binary serialization
└── export         Binary STL file writing
```

The design separates **distance evaluation** (primitives, shapes) from **mesh generation** (marching cubes) and **file I/O** (export). This means you can use the SDF layer for ray marching, collision detection, or any other distance-field application without pulling in mesh dependencies.

## Dependencies

| Crate | Version | Purpose |
|---|---|---|
| [`nalgebra`](https://nalgebra.org) | 0.33 | Linear algebra (`Vector3`, `Matrix3`) |
| [`rayon`](https://docs.rs/rayon) | 1.10 | Parallel mesh extraction |
| [`isosurface`](https://docs.rs/isosurface) | 0.1.0-alpha.0 | Marching cubes algorithm |
| [`approx`](https://docs.rs/approx) | 0.5 | Floating-point test assertions (dev only) |

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT License ([LICENSE-MIT](LICENSE-MIT))

at your option.
