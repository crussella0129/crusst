# Crusst — Rust Based Geometry Kernel

A standalone SDF (Signed Distance Function) geometry kernel for constructive solid geometry, mesh extraction, and parametric transport operations.

## Features

- **Primitive SDFs** — sphere, box, cylinder
- **CSG operations** — union, intersection, difference, smooth union
- **Shape system** — trait-based composable SDF shapes (Sphere, Box3, HalfSpace, CappedCone, Union, Intersection, Difference)
- **Path types** — LinePath, HelixPath, SpiralPath with tangent/frame computation
- **Transport orders** — Order 0 (identity), Order 1 (scaled sweep), Order 2 (tube sweep), Order 3 (tapered + twisted sweep)
- **Mesh extraction** — parallel marching cubes via `isosurface` with central-difference normals
- **STL export** — binary STL file writing

## Usage

```toml
[dependencies]
crusst = { git = "https://github.com/crussella0129/crusst.git" }
```

```rust
use crusst::shape::Sphere;
use crusst::mesh::extract_mesh;
use nalgebra::Vector3;

let sphere = Sphere::new(Vector3::zeros(), 5.0);
let mesh = extract_mesh(
    &sphere,
    Vector3::new(-6.0, -6.0, -6.0),
    Vector3::new(6.0, 6.0, 6.0),
    64,
);
println!("{} vertices, {} triangles", mesh.vertices.len(), mesh.indices.len() / 3);
```

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT License ([LICENSE-MIT](LICENSE-MIT))

at your option.
