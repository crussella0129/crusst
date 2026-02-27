use crusst::qef::solve_qef;
use crusst::types::BBox3;
use nalgebra::Vector3;

#[test]
fn qef_single_plane_intersection() {
    // One crossing at (5, 0, 0) with normal (1, 0, 0)
    // QEF solution should be at (5, 0, 0)
    let positions = vec![Vector3::new(5.0, 0.0, 0.0)];
    let normals = vec![Vector3::new(1.0, 0.0, 0.0)];
    let bounds = BBox3::new(Vector3::new(0.0, -5.0, -5.0), Vector3::new(10.0, 5.0, 5.0));
    let v = solve_qef(&positions, &normals, &bounds);
    assert!((v.x - 5.0).abs() < 1e-4);
}

#[test]
fn qef_two_perpendicular_planes() {
    // Edge: crossing at (5, 0, z) with normal (1,0,0) and (0, 0, z) at normal (0,0,1)
    // QEF should place vertex near the edge intersection
    let positions = vec![
        Vector3::new(5.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 5.0),
    ];
    let normals = vec![
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ];
    let bounds = BBox3::new(Vector3::new(0.0, -5.0, 0.0), Vector3::new(10.0, 5.0, 10.0));
    let v = solve_qef(&positions, &normals, &bounds);
    assert!((v.x - 5.0).abs() < 1e-2);
    assert!((v.z - 5.0).abs() < 1e-2);
}

#[test]
fn qef_vertex_clamped_to_bounds() {
    // Adversarial case: normals that would pull vertex far outside bounds
    let positions = vec![Vector3::new(1.0, 0.0, 0.0)];
    let normals = vec![Vector3::new(0.0, 1.0, 0.0)]; // perpendicular to crossing — degenerate
    let bounds = BBox3::new(Vector3::zeros(), Vector3::new(2.0, 2.0, 2.0));
    let v = solve_qef(&positions, &normals, &bounds);
    assert!(bounds.contains(v));
}
