//! Tests for all export formats.

use crusst::builder::Shape;
use crusst::types::TessSettings;

fn settings() -> TessSettings {
    TessSettings {
        chord_tolerance: 0.05,
        max_edge_length: 5.0,
        min_subdivisions: 8,
    }
}

// --- STL tests ---

#[test]
fn stl_all_primitives() {
    let shapes: Vec<(&str, Shape)> = vec![
        ("box", Shape::box3(5.0, 3.0, 8.0)),
        ("sphere", Shape::sphere(10.0)),
        ("cylinder", Shape::cylinder(5.0, 20.0)),
        ("cone", Shape::cone(8.0, 2.0, 20.0)),
        ("torus", Shape::torus(10.0, 3.0)),
        ("wedge", Shape::wedge(5.0, 4.0, 10.0)),
        ("capsule", Shape::capsule(3.0, 12.0)),
    ];

    for (name, shape) in &shapes {
        let mut buf = Vec::new();
        shape.write_stl(&settings(), &mut buf).unwrap();
        assert!(buf.len() > 84, "{name} STL should be larger than header");
    }
}

// --- OBJ tests ---

#[test]
fn obj_all_primitives() {
    let shapes: Vec<(&str, Shape)> = vec![
        ("box", Shape::box3(5.0, 3.0, 8.0)),
        ("sphere", Shape::sphere(10.0)),
        ("cylinder", Shape::cylinder(5.0, 20.0)),
    ];

    for (name, shape) in &shapes {
        let mut buf = Vec::new();
        shape.write_obj(&settings(), &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("v "), "{name} OBJ should contain vertex lines");
        assert!(text.contains("f "), "{name} OBJ should contain face lines");
    }
}

// --- STEP tests ---

#[test]
fn step_all_primitives() {
    let shapes: Vec<(&str, Shape)> = vec![
        ("box", Shape::box3(5.0, 3.0, 8.0)),
        ("sphere", Shape::sphere(10.0)),
        ("cylinder", Shape::cylinder(5.0, 20.0)),
        ("cone", Shape::cone(8.0, 2.0, 20.0)),
        ("torus", Shape::torus(10.0, 3.0)),
    ];

    for (name, shape) in &shapes {
        let mut buf = Vec::new();
        shape.write_step(&mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("MANIFOLD_SOLID_BREP"), "{name} STEP should contain MANIFOLD_SOLID_BREP");
        assert!(text.contains("END-ISO-10303-21;"), "{name} STEP should be complete");
    }
}

#[test]
fn step_surface_types() {
    // Box → PLANE only
    let mut buf = Vec::new();
    Shape::box3(5.0, 3.0, 8.0).write_step(&mut buf).unwrap();
    let text = String::from_utf8(buf).unwrap();
    assert!(text.contains("PLANE"));

    // Sphere → SPHERICAL_SURFACE
    let mut buf = Vec::new();
    Shape::sphere(10.0).write_step(&mut buf).unwrap();
    let text = String::from_utf8(buf).unwrap();
    assert!(text.contains("SPHERICAL_SURFACE"));

    // Cylinder → CYLINDRICAL_SURFACE + PLANE (caps)
    let mut buf = Vec::new();
    Shape::cylinder(5.0, 20.0).write_step(&mut buf).unwrap();
    let text = String::from_utf8(buf).unwrap();
    assert!(text.contains("CYLINDRICAL_SURFACE"));
    assert!(text.contains("PLANE"));

    // Torus → TOROIDAL_SURFACE
    let mut buf = Vec::new();
    Shape::torus(10.0, 3.0).write_step(&mut buf).unwrap();
    let text = String::from_utf8(buf).unwrap();
    assert!(text.contains("TOROIDAL_SURFACE"));
}

// --- 3MF tests ---

#[test]
fn threemf_all_primitives() {
    let shapes: Vec<(&str, Shape)> = vec![
        ("box", Shape::box3(5.0, 3.0, 8.0)),
        ("sphere", Shape::sphere(10.0)),
        ("cylinder", Shape::cylinder(5.0, 20.0)),
    ];

    for (name, shape) in &shapes {
        let mut buf = Vec::new();
        shape.write_3mf(&settings(), &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("<model"), "{name} 3MF should contain model tag");
        assert!(text.contains("<vertex"), "{name} 3MF should contain vertices");
        assert!(text.contains("<triangle"), "{name} 3MF should contain triangles");
    }
}
