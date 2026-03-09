//! Tests for primitive B-Rep solid constructors.

use crusst::primitive::*;
use crusst::topo::*;

#[test]
fn box_euler_formula() {
    let mut store = TopoStore::new();
    let solid = make_box(&mut store, 5.0, 3.0, 8.0);

    let v = store.solid_vertices(solid).len();
    let e = store.solid_edges(solid).len();
    let f = store.solid_face_count(solid);

    assert_eq!(v, 8, "Box should have 8 vertices");
    assert_eq!(e, 12, "Box should have 12 edges");
    assert_eq!(f, 6, "Box should have 6 faces");
    assert_eq!(v as i64 - e as i64 + f as i64, 2, "Euler V-E+F should be 2");
}

#[test]
fn box_topology_valid() {
    let mut store = TopoStore::new();
    let solid = make_box(&mut store, 5.0, 3.0, 8.0);
    let result = validate_solid(&store, solid);
    assert!(result.valid, "Box topology should be valid: {:?}", result.errors);
}

#[test]
fn box_vertex_positions() {
    let mut store = TopoStore::new();
    let solid = make_box(&mut store, 5.0, 3.0, 8.0);
    let verts = store.solid_vertices(solid);

    for &vid in &verts {
        let p = store.vertex(vid).point;
        assert!((p.x.abs() - 5.0).abs() < 1e-12 || p.x.abs() < 1e-12,
            "Box vertex x should be ±5.0");
        assert!((p.y.abs() - 3.0).abs() < 1e-12 || p.y.abs() < 1e-12);
        assert!((p.z.abs() - 8.0).abs() < 1e-12 || p.z.abs() < 1e-12);
    }
}

#[test]
fn sphere_euler_formula() {
    let mut store = TopoStore::new();
    let solid = make_sphere(&mut store, 10.0);

    let v = store.solid_vertices(solid).len();
    let e = store.solid_edges(solid).len();
    let f = store.solid_face_count(solid);

    assert_eq!(v, 8, "Cube-sphere should have 8 vertices");
    assert_eq!(e, 12, "Cube-sphere should have 12 edges");
    assert_eq!(f, 6, "Cube-sphere should have 6 faces");
    assert_eq!(v as i64 - e as i64 + f as i64, 2);
}

#[test]
fn sphere_topology_valid() {
    let mut store = TopoStore::new();
    let solid = make_sphere(&mut store, 10.0);
    let result = validate_solid(&store, solid);
    assert!(result.valid, "Sphere topology should be valid: {:?}", result.errors);
}

#[test]
fn sphere_vertices_on_sphere() {
    let mut store = TopoStore::new();
    let r = 10.0;
    let solid = make_sphere(&mut store, r);
    let verts = store.solid_vertices(solid);

    for &vid in &verts {
        let p = store.vertex(vid).point;
        let dist = p.coords.norm();
        assert!(
            (dist - r).abs() < 1e-10,
            "Sphere vertex at distance {dist} from origin, expected {r}"
        );
    }
}

#[test]
fn cylinder_euler_formula() {
    let mut store = TopoStore::new();
    let solid = make_cylinder(&mut store, 5.0, 20.0);

    let v = store.solid_vertices(solid).len();
    let e = store.solid_edges(solid).len();
    let f = store.solid_face_count(solid);

    assert_eq!(v, 4, "Cylinder should have 4 vertices");
    assert_eq!(e, 6, "Cylinder should have 6 edges");
    assert_eq!(f, 4, "Cylinder should have 4 faces");
    assert_eq!(v as i64 - e as i64 + f as i64, 2);
}

#[test]
fn cylinder_topology_valid() {
    let mut store = TopoStore::new();
    let solid = make_cylinder(&mut store, 5.0, 20.0);
    let result = validate_solid(&store, solid);
    assert!(result.valid, "Cylinder topology should be valid: {:?}", result.errors);
}

#[test]
fn cone_euler_formula() {
    let mut store = TopoStore::new();
    let solid = make_cone(&mut store, 8.0, 2.0, 20.0);

    let v = store.solid_vertices(solid).len();
    let e = store.solid_edges(solid).len();
    let f = store.solid_face_count(solid);

    assert_eq!(v, 4);
    assert_eq!(e, 6);
    assert_eq!(f, 4);
    assert_eq!(v as i64 - e as i64 + f as i64, 2);
}

#[test]
fn cone_topology_valid() {
    let mut store = TopoStore::new();
    let solid = make_cone(&mut store, 8.0, 2.0, 20.0);
    let result = validate_solid(&store, solid);
    assert!(result.valid, "Cone topology should be valid: {:?}", result.errors);
}

#[test]
fn torus_euler_formula() {
    let mut store = TopoStore::new();
    let solid = make_torus(&mut store, 10.0, 3.0);

    let v = store.solid_vertices(solid).len();
    let e = store.solid_edges(solid).len();
    let f = store.solid_face_count(solid);

    assert_eq!(v, 4, "Torus should have 4 vertices");
    assert_eq!(e, 8, "Torus should have 8 edges");
    assert_eq!(f, 4, "Torus should have 4 faces");
    // Torus is genus-1: V-E+F = 0 (not 2)
    assert_eq!(v as i64 - e as i64 + f as i64, 0, "Torus Euler should be 0 (genus 1)");
}

#[test]
fn torus_wires_closed() {
    let mut store = TopoStore::new();
    let solid = make_torus(&mut store, 10.0, 3.0);
    let shell = store.solid(solid);
    let faces = &store.shell(shell.outer_shell).faces;

    for &face_id in faces {
        let face = store.face(face_id);
        let coedges = store.wire_coedges(face.outer_wire);
        assert!(coedges.len() >= 3, "Torus face wire should have at least 3 edges");
    }
}

#[test]
fn torus_vertices_on_torus() {
    let mut store = TopoStore::new();
    let major_r = 10.0;
    let minor_r = 3.0;
    let solid = make_torus(&mut store, major_r, minor_r);
    let verts = store.solid_vertices(solid);

    for &vid in &verts {
        let p = store.vertex(vid).point;
        let r_xy = (p.x * p.x + p.y * p.y).sqrt();
        let dist_from_tube_center = ((r_xy - major_r) * (r_xy - major_r) + p.z * p.z).sqrt();
        assert!(
            (dist_from_tube_center - minor_r).abs() < 1e-10,
            "Torus vertex tube distance {dist_from_tube_center}, expected {minor_r}"
        );
    }
}

// --- Wedge tests ---

#[test]
fn wedge_euler_formula() {
    let mut store = TopoStore::new();
    let solid = make_wedge(&mut store, 5.0, 4.0, 10.0);

    let v = store.solid_vertices(solid).len();
    let e = store.solid_edges(solid).len();
    let f = store.solid_face_count(solid);

    assert_eq!(v, 6, "Wedge should have 6 vertices");
    assert_eq!(e, 9, "Wedge should have 9 edges");
    assert_eq!(f, 5, "Wedge should have 5 faces");
    assert_eq!(v as i64 - e as i64 + f as i64, 2, "Euler V-E+F should be 2");
}

#[test]
fn wedge_topology_valid() {
    let mut store = TopoStore::new();
    let solid = make_wedge(&mut store, 5.0, 4.0, 10.0);
    let result = validate_solid(&store, solid);
    assert!(result.valid, "Wedge topology should be valid: {:?}", result.errors);
}

// --- Capsule tests ---

#[test]
fn capsule_euler_formula() {
    let mut store = TopoStore::new();
    let solid = make_capsule(&mut store, 3.0, 12.0);

    let v = store.solid_vertices(solid).len();
    let e = store.solid_edges(solid).len();
    let f = store.solid_face_count(solid);

    assert_eq!(v, 8, "Capsule should have 8 vertices");
    assert_eq!(e, 12, "Capsule should have 12 edges");
    assert_eq!(f, 6, "Capsule should have 6 faces");
    assert_eq!(v as i64 - e as i64 + f as i64, 2);
}

#[test]
fn capsule_topology_valid() {
    let mut store = TopoStore::new();
    let solid = make_capsule(&mut store, 3.0, 12.0);
    let result = validate_solid(&store, solid);
    assert!(result.valid, "Capsule topology should be valid: {:?}", result.errors);
}

#[test]
fn capsule_symmetric() {
    // Capsule should be symmetric about the XY plane (z=0)
    let mut store = TopoStore::new();
    let solid = make_capsule(&mut store, 3.0, 12.0);
    let verts = store.solid_vertices(solid);

    for &vid in &verts {
        let p = store.vertex(vid).point;
        // Each vertex at z should have a mirror at -z
        let has_mirror = verts.iter().any(|&vid2| {
            let p2 = store.vertex(vid2).point;
            (p.x - p2.x).abs() < 1e-10
                && (p.y - p2.y).abs() < 1e-10
                && (p.z + p2.z).abs() < 1e-10
        });
        assert!(has_mirror, "Capsule vertex {:?} has no mirror", p);
    }
}
