use crusst::dag::SdfNode;
use crusst::feature::{FeatureKind, ft};
use nalgebra::Vector3;
use std::sync::Arc;

fn unit_box() -> SdfNode {
    SdfNode::Box3 {
        center: Vector3::new(0.0, 0.0, 0.0),
        half_extents: Vector3::new(1.0, 1.0, 1.0),
    }
}

fn unit_cylinder() -> SdfNode {
    SdfNode::Cylinder {
        base: Vector3::new(0.0, 0.0, 0.0),
        axis: Vector3::new(0.0, 1.0, 0.0),
        radius: 1.0,
        height: 2.0,
    }
}

#[test]
fn box3_face_count() {
    let b = unit_box();
    let faces = b.face_info().unwrap();
    assert_eq!(faces.len(), 6);
}

#[test]
fn box3_edge_count() {
    let b = unit_box();
    let edges = b.edge_info().unwrap();
    assert_eq!(edges.len(), 12);
}

#[test]
fn box3_face_labels() {
    let b = unit_box();
    let faces = b.face_info().unwrap();
    assert_eq!(faces[0].label, "+X");
    assert_eq!(faces[5].label, "-Z");
}

#[test]
fn box3_closest_face() {
    let b = unit_box();
    // Point at (2, 0, 0): closest to the +X face (index 0)
    let face = b.closest_face(Vector3::new(2.0, 0.0, 0.0)).unwrap();
    assert_eq!(face, 0);
}

#[test]
fn box3_face_distance() {
    let b = unit_box();
    // Face 0 (+X) distance at (2, 0, 0) = 2.0 - 1.0 = 1.0
    let dist = b.face_distance(Vector3::new(2.0, 0.0, 0.0), 0).unwrap();
    assert!((dist - 1.0).abs() < 1e-10);
}

#[test]
fn cylinder_face_count() {
    let c = unit_cylinder();
    let faces = c.face_info().unwrap();
    assert_eq!(faces.len(), 3);
}

#[test]
fn cylinder_edge_count() {
    let c = unit_cylinder();
    let edges = c.edge_info().unwrap();
    assert_eq!(edges.len(), 2);
}

#[test]
fn translate_preserves_faces() {
    let b = unit_box();
    let translated = SdfNode::Translate(Arc::new(b), Vector3::new(5.0, 0.0, 0.0));
    let faces = translated.face_info().unwrap();
    assert_eq!(faces.len(), 6);
}

#[test]
fn ft_builder() {
    let target = ft(1, 1).edges(&[2, 3, 4]);
    assert_eq!(target.component, 1);
    assert_eq!(target.body, 1);
    assert_eq!(target.kind, FeatureKind::Edge);
    assert_eq!(target.indices, vec![2, 3, 4]);
}
