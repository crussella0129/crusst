use approx::assert_relative_eq;
use nalgebra::Vector3;
use crusst::shape::{Sdf, Sphere};
use crusst::path::{LinePath, HelixPath, SpiralPath};
use crusst::transport::{order0, order1, order2, order3};
use crusst::mesh::extract_mesh;

#[test]
fn order0_is_just_the_section() {
    // Order 0: S = C(s) — shape IS the cross-section
    let section = Sphere::new(Vector3::zeros(), 5.0);
    let shape = order0(section);
    assert_relative_eq!(shape.evaluate(Vector3::zeros()), -5.0, epsilon = 1e-6);
    assert_relative_eq!(shape.evaluate(Vector3::new(5.0, 0.0, 0.0)), 0.0, epsilon = 1e-6);
}

#[test]
#[ignore]
fn order1_cone_eigenform() {
    // Order 1 eigenform: Cone — linear scaling along path
    let path = LinePath::new(Vector3::zeros(), Vector3::new(0.0, 0.0, 20.0));
    let section_radius = 10.0;
    let scale_fn = |t: f64| 1.0 - t;  // Taper to point

    let shape = order1(path, section_radius, scale_fn, 128);

    // At base (z=0), radius should be ~10
    let base_edge = shape.evaluate(Vector3::new(10.0, 0.0, 0.0));
    assert!(base_edge.abs() < 1.0, "Base edge should be near surface, got {base_edge}");

    // At tip (z=20), radius should be ~0 — point should be outside
    let tip = shape.evaluate(Vector3::new(5.0, 0.0, 20.0));
    assert!(tip > 0.0, "Should be outside at tip, got {tip}");

    // Can extract mesh
    let mesh = extract_mesh(
        &shape,
        Vector3::new(-12.0, -12.0, -1.0),
        Vector3::new(12.0, 12.0, 21.0),
        64,
    );
    assert!(mesh.vertices.len() > 50, "Cone mesh should have vertices");
}

#[test]
#[ignore]
fn order2_helix_eigenform() {
    // Order 2 eigenform: Helix — circular section swept along helix path
    let path = HelixPath::new(15.0, 8.0, 3.0);  // NOTE: turns is f64 now
    let section_radius = 2.0;

    let shape = order2(path, section_radius, 128);

    // Point on helix path should be inside (negative SDF)
    let on_path = Vector3::new(15.0, 0.0, 0.0); // t=0 start of helix
    let d = shape.evaluate(on_path);
    assert!(d < 0.0, "Point on path should be inside, got {d}");

    // Point far from helix should be outside
    let far = shape.evaluate(Vector3::new(0.0, 0.0, 0.0));
    assert!(far > 0.0, "Point at origin should be outside helix, got {far}");

    let mesh = extract_mesh(
        &shape,
        Vector3::new(-20.0, -20.0, -2.0),
        Vector3::new(20.0, 20.0, 26.0),
        64,
    );
    assert!(mesh.vertices.len() > 100, "Helix mesh should have many vertices");
}

#[test]
#[ignore]
fn order3_horn_eigenform() {
    // Order 3 eigenform: Horn — section tapers and rotates along a spiral path
    let path = SpiralPath::new(
        |t| 20.0 * (1.0 - t * 0.8),  // Radius decreases
        |t| 50.0 * t,                  // Height increases
        2.0,  // NOTE: turns is f64
    );
    let section_radius = 8.0;
    let scale_fn = |t: f64| 1.0 - 0.9 * t;    // Taper to near-point
    let twist_fn = |t: f64| t * std::f64::consts::PI;  // Half-turn twist

    let shape = order3(path, section_radius, scale_fn, twist_fn, 200);

    // Base of horn should be wide and inside
    let base = shape.evaluate(Vector3::new(20.0, 0.0, 0.0));
    assert!(base < 0.0, "Base of horn should be inside, got {base}");

    // Tip should be narrow — point away from path should be outside
    let near_tip = shape.evaluate(Vector3::new(10.0, 0.0, 45.0));
    assert!(near_tip > 0.0, "Near tip but offset should be outside, got {near_tip}");

    // Should produce a mesh
    let mesh = extract_mesh(
        &shape,
        Vector3::new(-30.0, -30.0, -5.0),
        Vector3::new(30.0, 30.0, 55.0),
        64,
    );
    assert!(mesh.vertices.len() > 100, "Horn mesh should have many vertices");
}
