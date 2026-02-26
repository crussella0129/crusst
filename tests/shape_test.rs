use approx::assert_relative_eq;
use nalgebra::Vector3;
use crusst::shape::{Sdf, Sphere, Box3, Union, Difference};

#[test]
fn sphere_shape_evaluates() {
    let s = Sphere::new(Vector3::zeros(), 5.0);
    assert_relative_eq!(s.evaluate(Vector3::zeros()), -5.0, epsilon = 1e-6);
    assert_relative_eq!(s.evaluate(Vector3::new(5.0, 0.0, 0.0)), 0.0, epsilon = 1e-6);
}

#[test]
fn union_shape_combines() {
    let a = Sphere::new(Vector3::zeros(), 1.0);
    let b = Sphere::new(Vector3::new(3.0, 0.0, 0.0), 1.0);
    let u = Union::new(a, b);
    assert_relative_eq!(u.evaluate(Vector3::zeros()), -1.0, epsilon = 1e-6);
}

#[test]
fn difference_shape_subtracts() {
    let outer = Sphere::new(Vector3::zeros(), 2.0);
    let inner = Sphere::new(Vector3::zeros(), 1.0);
    let shell = Difference::new(outer, inner);
    assert!(shell.evaluate(Vector3::zeros()) > 0.0);
    assert!(shell.evaluate(Vector3::new(1.5, 0.0, 0.0)) < 0.0);
}
