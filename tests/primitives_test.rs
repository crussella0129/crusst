use approx::assert_relative_eq;
use nalgebra::Vector3;
use crusst::primitives::{sdf_sphere, sdf_box, sdf_cylinder};

#[test]
fn sphere_center_is_negative() {
    let d = sdf_sphere(Vector3::zeros(), Vector3::zeros(), 1.0);
    assert_relative_eq!(d, -1.0, epsilon = 1e-6);
}

#[test]
fn sphere_surface_is_zero() {
    let d = sdf_sphere(Vector3::new(1.0, 0.0, 0.0), Vector3::zeros(), 1.0);
    assert_relative_eq!(d, 0.0, epsilon = 1e-6);
}

#[test]
fn sphere_outside_is_positive() {
    let d = sdf_sphere(Vector3::new(3.0, 0.0, 0.0), Vector3::zeros(), 1.0);
    assert_relative_eq!(d, 2.0, epsilon = 1e-6);
}

#[test]
fn box_center_is_negative() {
    let d = sdf_box(Vector3::zeros(), Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
    assert!(d < 0.0);
}

#[test]
fn box_face_center_is_zero() {
    // Point on the center of the +X face of a unit box
    let d = sdf_box(Vector3::new(1.0, 0.0, 0.0), Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
    assert_relative_eq!(d, 0.0, epsilon = 1e-6);
}

#[test]
fn box_outside_is_positive() {
    let d = sdf_box(Vector3::new(3.0, 0.0, 0.0), Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
    assert_relative_eq!(d, 2.0, epsilon = 1e-6);
}

#[test]
fn cylinder_on_axis_is_negative() {
    let d = sdf_cylinder(
        Vector3::new(0.0, 0.0, 5.0),
        Vector3::zeros(),
        Vector3::new(0.0, 0.0, 1.0),
        5.0,
        10.0,
    );
    assert!(d < 0.0);
}

#[test]
fn cylinder_outside_radially() {
    let d = sdf_cylinder(
        Vector3::new(10.0, 0.0, 5.0),
        Vector3::zeros(),
        Vector3::new(0.0, 0.0, 1.0),
        5.0,
        10.0,
    );
    assert_relative_eq!(d, 5.0, epsilon = 1e-6);
}
