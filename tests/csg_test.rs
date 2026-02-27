use approx::assert_relative_eq;
use crusst::csg::{difference, intersection, smooth_union, union};
use crusst::primitives::sdf_sphere;
use nalgebra::Vector3;

#[test]
fn union_takes_minimum() {
    let a = sdf_sphere(Vector3::new(2.0, 0.0, 0.0), Vector3::zeros(), 1.0);
    let b = sdf_sphere(
        Vector3::new(2.0, 0.0, 0.0),
        Vector3::new(5.0, 0.0, 0.0),
        1.0,
    );
    assert_relative_eq!(union(a, b), 1.0, epsilon = 1e-6);
}

#[test]
fn intersection_takes_maximum() {
    let a = sdf_sphere(Vector3::new(0.5, 0.0, 0.0), Vector3::zeros(), 1.0);
    let b = sdf_sphere(
        Vector3::new(0.5, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        1.0,
    );
    assert_relative_eq!(intersection(a, b), -0.5, epsilon = 1e-6);
}

#[test]
fn difference_subtracts() {
    let a = sdf_sphere(Vector3::zeros(), Vector3::zeros(), 2.0);
    let b = sdf_sphere(Vector3::zeros(), Vector3::zeros(), 1.0);
    let d = difference(a, b);
    assert!(d > 0.0, "Point inside both should be outside A-B");
}

#[test]
fn smooth_union_blends() {
    let a = 1.0_f64;
    let b = 1.0_f64;
    let k = 0.5;
    let d = smooth_union(a, b, k);
    assert!(d < a.min(b), "Smooth union should blend below minimum");
}
