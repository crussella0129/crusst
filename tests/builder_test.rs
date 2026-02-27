use crusst::builder::Shape;
use crusst::types::MeshSettings;
use nalgebra::Vector3;

#[test]
fn builder_sphere_evaluates() {
    let s = Shape::sphere(5.0);
    assert!((s.distance(Vector3::zeros()) - (-5.0)).abs() < 1e-10);
}

#[test]
fn builder_chain_union() {
    let s = Shape::sphere(5.0).union(Shape::sphere(5.0).translate(10.0, 0.0, 0.0));
    assert!(s.contains(Vector3::zeros()));
    assert!(s.contains(Vector3::new(10.0, 0.0, 0.0)));
    assert!(!s.contains(Vector3::new(100.0, 0.0, 0.0)));
}

#[test]
fn builder_subtract() {
    let s = Shape::box3(10.0, 10.0, 10.0).subtract(Shape::sphere(7.0));
    assert!(!s.contains(Vector3::zeros()));
    assert!(s.contains(Vector3::new(9.0, 9.0, 0.0)));
}

#[test]
fn builder_mesh_produces_triangles() {
    let s = Shape::sphere(5.0);
    let mesh = s.mesh(MeshSettings {
        max_depth: 5,
        min_depth: 3,
        edge_tolerance: 1e-6,
    });
    assert!(mesh.indices.len() > 30);
}

#[test]
fn builder_clone_is_cheap() {
    let s = Shape::sphere(5.0).translate(1.0, 2.0, 3.0);
    let s2 = s.clone();
    assert!((s.distance(Vector3::zeros()) - s2.distance(Vector3::zeros())).abs() < 1e-15);
}

#[test]
fn box_has_six_faces() {
    let b = Shape::box3(1.0, 1.0, 1.0);
    let faces = b.faces().expect("box should have faces");
    assert_eq!(faces.len(), 6);
}

#[test]
fn box_has_twelve_edges() {
    let b = Shape::box3(1.0, 1.0, 1.0);
    let edges = b.edges().expect("box should have edges");
    assert_eq!(edges.len(), 12);
}

#[test]
fn translated_box_preserves_face_count() {
    let b = Shape::box3(1.0, 1.0, 1.0).translate(5.0, 0.0, 0.0);
    let faces = b.faces().expect("translated box should have faces");
    assert_eq!(faces.len(), 6);
}
