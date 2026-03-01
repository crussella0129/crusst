//! Tests for B-Rep tessellation — converts exact geometry to triangle meshes.

use crusst::primitive::*;
use crusst::tessellate::tessellate_solid;
use crusst::topo::TopoStore;
use crusst::types::TessSettings;

fn default_settings() -> TessSettings {
    TessSettings {
        chord_tolerance: 0.05,
        max_edge_length: 5.0,
        min_subdivisions: 8,
    }
}

fn assert_manifold_mesh(mesh: &crusst::types::TriangleMesh, label: &str) {
    assert!(!mesh.vertices.is_empty(), "{label}: no vertices");
    assert!(!mesh.indices.is_empty(), "{label}: no indices");
    assert_eq!(mesh.indices.len() % 3, 0, "{label}: indices not multiple of 3");
    assert_eq!(mesh.vertices.len(), mesh.normals.len(), "{label}: vertices/normals mismatch");

    // All indices in bounds
    let nv = mesh.vertices.len() as u32;
    for &idx in &mesh.indices {
        assert!(idx < nv, "{label}: index {idx} out of bounds (nv={nv})");
    }

    // All normals should be unit length
    for (i, n) in mesh.normals.iter().enumerate() {
        let len = n.norm();
        assert!(
            (len - 1.0).abs() < 0.1,
            "{label}: normal {i} has length {len}"
        );
    }
}

#[test]
fn tessellate_box() {
    let mut store = TopoStore::new();
    let solid = make_box(&mut store, 5.0, 3.0, 8.0);
    let mesh = tessellate_solid(&store, solid, &default_settings());
    assert_manifold_mesh(&mesh, "box");

    // A box should have at least 12 triangles (2 per face × 6 faces minimum)
    let n_tris = mesh.indices.len() / 3;
    assert!(n_tris >= 12, "Box should have >= 12 triangles, got {n_tris}");
}

#[test]
fn tessellate_sphere() {
    let mut store = TopoStore::new();
    let solid = make_sphere(&mut store, 10.0);
    let mesh = tessellate_solid(&store, solid, &default_settings());
    assert_manifold_mesh(&mesh, "sphere");

    // All vertices should be approximately on the sphere
    for (i, v) in mesh.vertices.iter().enumerate() {
        let r = v.norm();
        assert!(
            (r - 10.0).abs() < 0.5,
            "Sphere vertex {i} at distance {r}, expected ~10.0"
        );
    }
}

#[test]
fn tessellate_cylinder() {
    let mut store = TopoStore::new();
    let solid = make_cylinder(&mut store, 5.0, 20.0);
    let mesh = tessellate_solid(&store, solid, &default_settings());
    assert_manifold_mesh(&mesh, "cylinder");

    let n_tris = mesh.indices.len() / 3;
    assert!(n_tris >= 20, "Cylinder should have >= 20 triangles, got {n_tris}");
}

#[test]
fn tessellate_cone() {
    let mut store = TopoStore::new();
    let solid = make_cone(&mut store, 8.0, 2.0, 20.0);
    let mesh = tessellate_solid(&store, solid, &default_settings());
    assert_manifold_mesh(&mesh, "cone");
}

#[test]
fn tessellate_torus() {
    let mut store = TopoStore::new();
    let solid = make_torus(&mut store, 10.0, 3.0);
    let mesh = tessellate_solid(&store, solid, &default_settings());
    assert_manifold_mesh(&mesh, "torus");

    // Torus vertices should be at the correct radial distance
    for v in &mesh.vertices {
        let r_xy = (v.x * v.x + v.y * v.y).sqrt();
        let tube_dist = ((r_xy - 10.0).powi(2) + v.z * v.z).sqrt();
        assert!(
            (tube_dist - 3.0).abs() < 0.5,
            "Torus vertex tube distance {tube_dist}, expected ~3.0"
        );
    }
}

#[test]
fn tessellate_wedge() {
    let mut store = TopoStore::new();
    let solid = make_wedge(&mut store, 5.0, 4.0, 10.0);
    let mesh = tessellate_solid(&store, solid, &default_settings());
    assert_manifold_mesh(&mesh, "wedge");
}

#[test]
fn tessellate_capsule() {
    let mut store = TopoStore::new();
    let solid = make_capsule(&mut store, 3.0, 12.0);
    let mesh = tessellate_solid(&store, solid, &default_settings());
    assert_manifold_mesh(&mesh, "capsule");
}

#[test]
fn tessellate_binary_roundtrip() {
    // Verify mesh serialization still works with new tessellation output
    let mut store = TopoStore::new();
    let solid = make_box(&mut store, 1.0, 1.0, 1.0);
    let mesh = tessellate_solid(&store, solid, &default_settings());
    let binary = mesh.to_binary();
    assert!(!binary.is_empty(), "Binary output should not be empty");

    // Decode vertex count
    let nv = u32::from_le_bytes([binary[0], binary[1], binary[2], binary[3]]);
    assert_eq!(nv as usize, mesh.vertices.len());
}
