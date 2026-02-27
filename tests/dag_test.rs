use crusst::dag::SdfNode;
use nalgebra::Vector3;
use std::sync::Arc;

#[test]
fn dag_sphere_at_center_is_negative() {
    let node = SdfNode::Sphere { center: Vector3::zeros(), radius: 10.0 };
    assert!((node.evaluate(Vector3::zeros()) - (-10.0)).abs() < 1e-10);
}

#[test]
fn dag_sphere_on_surface_is_zero() {
    let node = SdfNode::Sphere { center: Vector3::zeros(), radius: 10.0 };
    assert!(node.evaluate(Vector3::new(10.0, 0.0, 0.0)).abs() < 1e-10);
}

#[test]
fn dag_union_is_min() {
    let a = Arc::new(SdfNode::Sphere { center: Vector3::new(-5.0, 0.0, 0.0), radius: 3.0 });
    let b = Arc::new(SdfNode::Sphere { center: Vector3::new(5.0, 0.0, 0.0), radius: 3.0 });
    let u = SdfNode::Union(a, b);
    let d = u.evaluate(Vector3::zeros());
    assert!((d - 2.0).abs() < 1e-10);
}

#[test]
fn dag_translate_moves_shape() {
    let s = Arc::new(SdfNode::Sphere { center: Vector3::zeros(), radius: 5.0 });
    let t = SdfNode::Translate(s, Vector3::new(10.0, 0.0, 0.0));
    assert!((t.evaluate(Vector3::new(10.0, 0.0, 0.0)) - (-5.0)).abs() < 1e-10);
}

#[test]
fn dag_difference_subtracts() {
    let a = Arc::new(SdfNode::Box3 {
        center: Vector3::zeros(), half_extents: Vector3::new(10.0, 10.0, 10.0)
    });
    let b = Arc::new(SdfNode::Sphere { center: Vector3::zeros(), radius: 7.0 });
    let diff = SdfNode::Difference(a, b);
    assert!(diff.evaluate(Vector3::zeros()) > 0.0);
    assert!(diff.evaluate(Vector3::new(9.0, 9.0, 0.0)) < 0.0);
}
