//! Tests for the fluent builder API.

use crusst::builder::Shape;
use crusst::profile::Profile;
use crusst::types::TessSettings;

fn settings() -> TessSettings {
    TessSettings {
        chord_tolerance: 0.05,
        max_edge_length: 5.0,
        min_subdivisions: 8,
    }
}

#[test]
fn builder_box_mesh() {
    let mesh = Shape::box3(5.0, 3.0, 8.0).mesh(&settings());
    assert!(!mesh.vertices.is_empty());
    assert!(mesh.indices.len() >= 36, "Box needs at least 12 triangles");
}

#[test]
fn builder_sphere_mesh() {
    let mesh = Shape::sphere(10.0).mesh(&settings());
    assert!(!mesh.vertices.is_empty());
}

#[test]
fn builder_cylinder_mesh() {
    let mesh = Shape::cylinder(5.0, 20.0).mesh(&settings());
    assert!(!mesh.vertices.is_empty());
}

#[test]
fn builder_cone_mesh() {
    let mesh = Shape::cone(8.0, 2.0, 20.0).mesh(&settings());
    assert!(!mesh.vertices.is_empty());
}

#[test]
fn builder_torus_mesh() {
    let mesh = Shape::torus(10.0, 3.0).mesh(&settings());
    assert!(!mesh.vertices.is_empty());
}

#[test]
fn builder_wedge_mesh() {
    let mesh = Shape::wedge(5.0, 4.0, 10.0).mesh(&settings());
    assert!(!mesh.vertices.is_empty());
}

#[test]
fn builder_capsule_mesh() {
    let mesh = Shape::capsule(3.0, 12.0).mesh(&settings());
    assert!(!mesh.vertices.is_empty());
}

#[test]
fn builder_translate() {
    let shape = Shape::box3(1.0, 1.0, 1.0).translate(10.0, 0.0, 0.0);
    let mesh = shape.mesh(&settings());

    // All vertices should be shifted +10 in X
    for v in &mesh.vertices {
        assert!(v.x > 8.0, "Translated box vertex x should be > 8: got {}", v.x);
    }
}

#[test]
fn builder_scale() {
    let mesh_small = Shape::box3(1.0, 1.0, 1.0).mesh(&settings());
    let mesh_big = Shape::box3(1.0, 1.0, 1.0).scale(3.0).mesh(&settings());

    // The scaled box should have vertices 3x further from origin
    let max_small: f64 = mesh_small.vertices.iter().map(|v| v.norm()).fold(0.0, f64::max);
    let max_big: f64 = mesh_big.vertices.iter().map(|v| v.norm()).fold(0.0, f64::max);

    assert!(
        (max_big / max_small - 3.0).abs() < 0.1,
        "Scale 3x: max_big/max_small = {}",
        max_big / max_small
    );
}

#[test]
fn builder_rotate_z() {
    use std::f64::consts::FRAC_PI_2;
    let mesh = Shape::box3(5.0, 1.0, 1.0)
        .rotate_z(FRAC_PI_2) // 90° around Z: X becomes Y
        .mesh(&settings());

    // After 90° rotation, the long dimension (±5) should now be along Y
    let max_y: f64 = mesh.vertices.iter().map(|v| v.y.abs()).fold(0.0, f64::max);
    let max_x: f64 = mesh.vertices.iter().map(|v| v.x.abs()).fold(0.0, f64::max);
    assert!(max_y > max_x, "After 90° Z rotation, Y extent ({max_y}) should exceed X ({max_x})");
}

#[test]
fn builder_mirror_x() {
    let shape = Shape::box3(5.0, 3.0, 8.0).translate(10.0, 0.0, 0.0).mirror_x();
    let mesh = shape.mesh(&settings());

    // Box half-extents (5,3,8), translated +10 in X → vertices at x ∈ [5, 15].
    // After mirror_x → x ∈ [-15, -5].
    let x_min: f64 = mesh.vertices.iter().map(|v| v.x).fold(f64::MAX, f64::min);
    let x_max: f64 = mesh.vertices.iter().map(|v| v.x).fold(f64::MIN, f64::max);
    assert!(x_max < -4.5, "Mirrored box x_max should be ≈ -5: got {x_max}");
    assert!(x_min > -15.5, "Mirrored box x_min should be ≈ -15: got {x_min}");

    // Centroid should be near x = -10
    let cx: f64 = mesh.vertices.iter().map(|v| v.x).sum::<f64>() / mesh.vertices.len() as f64;
    assert!((cx + 10.0).abs() < 1.5, "Mirrored box centroid x should be ≈ -10: got {cx}");
}

#[test]
fn builder_extrude_rect() {
    let profile = Profile::rect(4.0, 2.0);
    let shape = Shape::extrude(&profile, 10.0);
    let mesh = shape.mesh(&settings());

    assert!(!mesh.vertices.is_empty(), "Extruded rect should produce vertices");
    assert!(mesh.indices.len() >= 12, "Should have some triangles");

    // Z values should span [0, 10]
    let z_min: f64 = mesh.vertices.iter().map(|v| v.z).fold(f64::MAX, f64::min);
    let z_max: f64 = mesh.vertices.iter().map(|v| v.z).fold(f64::MIN, f64::max);
    assert!(z_min < 0.01, "Extrude should start at z≈0: got {z_min}");
    assert!(z_max > 9.99, "Extrude should end at z≈10: got {z_max}");
}

#[test]
fn builder_revolve_circle() {
    let profile = Profile::rect(2.0, 1.0);
    let shape = Shape::revolve(&profile, std::f64::consts::TAU);
    let mesh = shape.mesh(&settings());
    assert!(!mesh.vertices.is_empty(), "Revolved rect should produce vertices");
}

#[test]
fn builder_validate_all_primitives() {
    let shapes = vec![
        ("box", Shape::box3(5.0, 3.0, 8.0)),
        ("sphere", Shape::sphere(10.0)),
        ("cylinder", Shape::cylinder(5.0, 20.0)),
        ("cone", Shape::cone(8.0, 2.0, 20.0)),
        ("torus", Shape::torus(10.0, 3.0)),
        ("wedge", Shape::wedge(5.0, 4.0, 10.0)),
        ("capsule", Shape::capsule(3.0, 12.0)),
    ];

    for (name, shape) in shapes {
        let result = shape.validate();
        // Torus has Euler characteristic 0 (genus 1), so it'll report as invalid
        // for our simple genus-0 validator. All others should be valid.
        if name != "torus" {
            assert!(result.valid, "{name} should have valid topology: {:?}", result.errors);
        }
    }
}
