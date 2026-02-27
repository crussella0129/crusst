use approx::assert_relative_eq;
use crusst::path::{HelixPath, LinePath, Path, SpiralPath};
use nalgebra::Vector3;

#[test]
fn line_path_endpoints() {
    let line = LinePath::new(Vector3::zeros(), Vector3::new(0.0, 0.0, 10.0));
    let start = line.point(0.0);
    let end = line.point(1.0);
    assert_relative_eq!(start, Vector3::zeros(), epsilon = 1e-6);
    assert_relative_eq!(end, Vector3::new(0.0, 0.0, 10.0), epsilon = 1e-6);
}

#[test]
fn line_path_midpoint() {
    let line = LinePath::new(Vector3::zeros(), Vector3::new(10.0, 0.0, 0.0));
    let mid = line.point(0.5);
    assert_relative_eq!(mid, Vector3::new(5.0, 0.0, 0.0), epsilon = 1e-6);
}

#[test]
fn line_path_tangent_is_constant() {
    let line = LinePath::new(Vector3::zeros(), Vector3::new(0.0, 0.0, 10.0));
    let t0 = line.tangent(0.0);
    let t1 = line.tangent(0.5);
    assert_relative_eq!(t0, t1, epsilon = 1e-6);
    assert_relative_eq!(t0.norm(), 1.0, epsilon = 1e-6);
}

#[test]
fn helix_path_start_on_circle() {
    let helix = HelixPath::new(10.0, 5.0, 3.0);
    let start = helix.point(0.0);
    // At t=0, helix starts at (radius, 0, 0)
    assert_relative_eq!(start.x, 10.0, epsilon = 1e-6);
    assert_relative_eq!(start.y, 0.0, epsilon = 1e-6);
    assert_relative_eq!(start.z, 0.0, epsilon = 1e-6);
}

#[test]
fn helix_path_rises_linearly() {
    let helix = HelixPath::new(10.0, 5.0, 2.0);
    let end = helix.point(1.0);
    // Total height = pitch * turns = 5 * 2 = 10
    assert_relative_eq!(end.z, 10.0, epsilon = 1e-6);
}

#[test]
fn spiral_path_radius_decreases() {
    let spiral = SpiralPath::new(
        |t| 20.0 * (1.0 - t), // radius decreases from 20 to 0
        |t| 50.0 * t,         // height increases from 0 to 50
        3.0,
    );
    let start = spiral.point(0.0);
    let end = spiral.point(1.0);
    assert!(start.x.abs() > end.x.abs() + end.y.abs());
}

#[test]
fn helix_tangent_is_unit_length() {
    let helix = HelixPath::new(10.0, 5.0, 3.0);
    for &t in &[0.0, 0.25, 0.5, 0.75, 1.0] {
        let tan = helix.tangent(t);
        assert_relative_eq!(tan.norm(), 1.0, epsilon = 1e-4);
    }
}

#[test]
fn helix_tangent_has_positive_z() {
    let helix = HelixPath::new(10.0, 5.0, 3.0);
    let tan = helix.tangent(0.5);
    assert!(tan.z > 0.0, "helix tangent should have upward z component");
}
