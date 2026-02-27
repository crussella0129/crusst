//! Showcase test suite — demonstrates engineering-grade geometry kernel capabilities.
//!
//! Each test exercises primitives, boolean operations, transforms, and complex
//! composed shapes, exporting to both STL and STEP for visual inspection and
//! CAD interop analysis.

use std::f64::consts::PI;
use std::path::PathBuf;

use nalgebra::{Rotation3, Unit, Vector3};

use crusst::csg;
use crusst::export::write_stl;
use crusst::mesh::extract_mesh;
use crusst::shape::*;
use crusst::step_export::write_step;

fn output_dir() -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("showcase");
    std::fs::create_dir_all(&dir).ok();
    dir
}

fn export_both(mesh: &crusst::mesh::TriangleMesh, name: &str) {
    let dir = output_dir();
    write_stl(mesh, &dir.join(format!("{}.stl", name))).unwrap();
    write_step(mesh, &dir.join(format!("{}.step", name))).unwrap();
}

// ===========================================================================
// Section 1: All Primitives
// ===========================================================================

#[test]
#[ignore]
fn primitive_sphere() {
    let shape = Sphere::new(Vector3::zeros(), 10.0);
    // Distance accuracy: center should be -radius
    assert!((shape.evaluate(Vector3::zeros()) - (-10.0)).abs() < 1e-10);
    // Surface point should be ~0
    assert!(shape.evaluate(Vector3::new(10.0, 0.0, 0.0)).abs() < 1e-10);
    // Outside should be positive
    assert!(shape.evaluate(Vector3::new(15.0, 0.0, 0.0)) > 0.0);

    let mesh = extract_mesh(&shape, Vector3::from_element(-12.0), Vector3::from_element(12.0), 64);
    assert!(mesh.indices.len() / 3 > 100);
    export_both(&mesh, "01_sphere");
}

#[test]
#[ignore]
fn primitive_box() {
    let shape = Box3::new(Vector3::zeros(), Vector3::new(5.0, 3.0, 8.0));
    // Center is inside
    assert!(shape.evaluate(Vector3::zeros()) < 0.0);
    // Face center should be ~0
    assert!(shape.evaluate(Vector3::new(5.0, 0.0, 0.0)).abs() < 1e-10);
    // Corner distance: sqrt(0 + 0 + 0) = 0 but at corner it's more complex
    assert!(shape.evaluate(Vector3::new(10.0, 0.0, 0.0)) > 0.0);

    let mesh = extract_mesh(
        &shape,
        Vector3::new(-7.0, -5.0, -10.0),
        Vector3::new(7.0, 5.0, 10.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 10);
    export_both(&mesh, "02_box");
}

#[test]
#[ignore]
fn primitive_cylinder() {
    let shape = Cylinder::new(
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
        5.0,
        20.0,
    );
    // On axis, inside
    assert!(shape.evaluate(Vector3::new(0.0, 0.0, 10.0)) < 0.0);
    // Outside radially
    assert!(shape.evaluate(Vector3::new(8.0, 0.0, 10.0)) > 0.0);

    let mesh = extract_mesh(
        &shape,
        Vector3::new(-7.0, -7.0, -2.0),
        Vector3::new(7.0, 7.0, 22.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 50);
    export_both(&mesh, "03_cylinder");
}

#[test]
#[ignore]
fn primitive_capped_cone() {
    let shape = CappedCone::new(
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 20.0),
        8.0,
        2.0,
    );
    // On axis near base, inside
    assert!(shape.evaluate(Vector3::new(0.0, 0.0, 2.0)) < 0.0);
    // Outside
    assert!(shape.evaluate(Vector3::new(15.0, 0.0, 10.0)) > 0.0);

    let mesh = extract_mesh(
        &shape,
        Vector3::new(-10.0, -10.0, -2.0),
        Vector3::new(10.0, 10.0, 22.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 50);
    export_both(&mesh, "04_capped_cone");
}

#[test]
#[ignore]
fn primitive_torus() {
    let shape = Torus::new(Vector3::zeros(), 10.0, 3.0);
    // Inside the tube at (10, 0, 0)
    assert!(shape.evaluate(Vector3::new(10.0, 0.0, 0.0)).abs() < 3.1);
    // Center of the torus (hole) should be positive
    assert!(shape.evaluate(Vector3::zeros()) > 0.0);
    // Far outside
    assert!(shape.evaluate(Vector3::new(20.0, 0.0, 0.0)) > 0.0);
    // On the tube surface at (10+3, 0, 0) should be ~0
    assert!(shape.evaluate(Vector3::new(13.0, 0.0, 0.0)).abs() < 0.1);

    let mesh = extract_mesh(
        &shape,
        Vector3::new(-15.0, -5.0, -15.0),
        Vector3::new(15.0, 5.0, 15.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 100);
    export_both(&mesh, "05_torus");
}

#[test]
#[ignore]
fn primitive_rounded_box() {
    let shape = RoundedBox::new(
        Vector3::zeros(),
        Vector3::new(5.0, 3.0, 8.0),
        1.0, // 1mm rounding
    );
    // Center inside
    assert!(shape.evaluate(Vector3::zeros()) < 0.0);
    // The rounded box is slightly larger than the sharp box
    // Face center at x=6 (5+1) should be ~0
    assert!(shape.evaluate(Vector3::new(6.0, 0.0, 0.0)).abs() < 0.2);

    let mesh = extract_mesh(
        &shape,
        Vector3::new(-8.0, -6.0, -11.0),
        Vector3::new(8.0, 6.0, 11.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 50);
    export_both(&mesh, "06_rounded_box");
}

#[test]
#[ignore]
fn primitive_capsule() {
    let shape = Capsule::new(
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 20.0),
        4.0,
    );
    // On axis, inside
    assert!(shape.evaluate(Vector3::new(0.0, 0.0, 10.0)) < 0.0);
    // At endpoint hemisphere, surface
    assert!(shape.evaluate(Vector3::new(4.0, 0.0, 0.0)).abs() < 0.1);

    let mesh = extract_mesh(
        &shape,
        Vector3::new(-6.0, -6.0, -6.0),
        Vector3::new(6.0, 6.0, 26.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 50);
    export_both(&mesh, "07_capsule");
}

#[test]
#[ignore]
fn primitive_ellipsoid() {
    let shape = Ellipsoid::new(
        Vector3::zeros(),
        Vector3::new(10.0, 5.0, 3.0),
    );
    // Center inside
    assert!(shape.evaluate(Vector3::zeros()) < 0.0);
    // Outside on long axis
    assert!(shape.evaluate(Vector3::new(15.0, 0.0, 0.0)) > 0.0);

    let mesh = extract_mesh(
        &shape,
        Vector3::new(-12.0, -7.0, -5.0),
        Vector3::new(12.0, 7.0, 5.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 50);
    export_both(&mesh, "08_ellipsoid");
}

#[test]
#[ignore]
fn primitive_rounded_cylinder() {
    let shape = RoundedCylinder::new(
        Vector3::zeros(),
        6.0,  // radius
        1.0,  // round radius
        10.0, // half-height
    );
    // On axis, inside
    assert!(shape.evaluate(Vector3::zeros()) < 0.0);
    // Far outside
    assert!(shape.evaluate(Vector3::new(12.0, 0.0, 0.0)) > 0.0);

    let mesh = extract_mesh(
        &shape,
        Vector3::new(-9.0, -13.0, -9.0),
        Vector3::new(9.0, 13.0, 9.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 50);
    export_both(&mesh, "09_rounded_cylinder");
}

// ===========================================================================
// Section 2: Boolean Operations (Sharp)
// ===========================================================================

#[test]
#[ignore]
fn boolean_union() {
    let a = Sphere::new(Vector3::new(-3.0, 0.0, 0.0), 5.0);
    let b = Sphere::new(Vector3::new(3.0, 0.0, 0.0), 5.0);
    let shape = Union::new(a, b);

    // Both centers inside
    assert!(shape.evaluate(Vector3::new(-3.0, 0.0, 0.0)) < 0.0);
    assert!(shape.evaluate(Vector3::new(3.0, 0.0, 0.0)) < 0.0);
    // Midpoint inside (overlap region)
    assert!(shape.evaluate(Vector3::zeros()) < 0.0);
    // Far outside
    assert!(shape.evaluate(Vector3::new(15.0, 0.0, 0.0)) > 0.0);

    let mesh = extract_mesh(
        &shape,
        Vector3::from_element(-10.0),
        Vector3::from_element(10.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 100);
    export_both(&mesh, "10_boolean_union");
}

#[test]
#[ignore]
fn boolean_intersection() {
    let a = Sphere::new(Vector3::new(-2.0, 0.0, 0.0), 5.0);
    let b = Sphere::new(Vector3::new(2.0, 0.0, 0.0), 5.0);
    let shape = Intersection::new(a, b);

    // Midpoint (in overlap) is inside
    assert!(shape.evaluate(Vector3::zeros()) < 0.0);
    // Center of sphere A — outside intersection if not inside B
    // A center at (-2,0,0): dist to B = sqrt(16) - 5 = -1, inside B
    assert!(shape.evaluate(Vector3::new(-2.0, 0.0, 0.0)) < 0.0);
    // Far side of A only (not in B)
    assert!(shape.evaluate(Vector3::new(-6.0, 0.0, 0.0)) > 0.0);

    let mesh = extract_mesh(
        &shape,
        Vector3::from_element(-8.0),
        Vector3::from_element(8.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 50);
    export_both(&mesh, "11_boolean_intersection");
}

#[test]
#[ignore]
fn boolean_difference() {
    // Box with a sphere hole
    let block = Box3::new(Vector3::zeros(), Vector3::new(10.0, 10.0, 10.0));
    let hole = Sphere::new(Vector3::zeros(), 7.0);
    let shape = Difference::new(block, hole);

    // Center (inside sphere, outside result)
    assert!(shape.evaluate(Vector3::zeros()) > 0.0);
    // Corner of box (outside sphere, inside box -> inside result)
    assert!(shape.evaluate(Vector3::new(9.0, 9.0, 9.0)) < 0.0);

    let mesh = extract_mesh(
        &shape,
        Vector3::from_element(-12.0),
        Vector3::from_element(12.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 50);
    export_both(&mesh, "12_boolean_difference");
}

// ===========================================================================
// Section 3: Boolean Operations (Smooth / Blended)
// ===========================================================================

#[test]
#[ignore]
fn smooth_union_blended() {
    let a = Sphere::new(Vector3::new(-4.0, 0.0, 0.0), 5.0);
    let b = Sphere::new(Vector3::new(4.0, 0.0, 0.0), 5.0);
    let shape = SmoothUnion::new(a, b, 2.0);

    // Both centers inside
    assert!(shape.evaluate(Vector3::new(-4.0, 0.0, 0.0)) < 0.0);
    assert!(shape.evaluate(Vector3::new(4.0, 0.0, 0.0)) < 0.0);
    // Midpoint is inside, and the blend makes it MORE negative than sharp union
    let smooth_val = shape.evaluate(Vector3::zeros());
    let sharp_val = csg::union(
        Sphere::new(Vector3::new(-4.0, 0.0, 0.0), 5.0).evaluate(Vector3::zeros()),
        Sphere::new(Vector3::new(4.0, 0.0, 0.0), 5.0).evaluate(Vector3::zeros()),
    );
    assert!(smooth_val <= sharp_val);

    let mesh = extract_mesh(
        &shape,
        Vector3::from_element(-11.0),
        Vector3::from_element(11.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 100);
    export_both(&mesh, "13_smooth_union");
}

#[test]
#[ignore]
fn smooth_intersection_blended() {
    let a = Sphere::new(Vector3::new(-2.0, 0.0, 0.0), 6.0);
    let b = Sphere::new(Vector3::new(2.0, 0.0, 0.0), 6.0);
    let shape = SmoothIntersection::new(a, b, 2.0);

    // Midpoint inside
    assert!(shape.evaluate(Vector3::zeros()) < 0.0);

    let mesh = extract_mesh(
        &shape,
        Vector3::from_element(-8.0),
        Vector3::from_element(8.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 50);
    export_both(&mesh, "14_smooth_intersection");
}

#[test]
#[ignore]
fn smooth_difference_blended() {
    let block = Box3::new(Vector3::zeros(), Vector3::new(8.0, 8.0, 8.0));
    let hole = Sphere::new(Vector3::zeros(), 6.0);
    let shape = SmoothDifference::new(block, hole, 1.5);

    // Corner of box, far from sphere -> inside result
    assert!(shape.evaluate(Vector3::new(7.5, 7.5, 7.5)) < 0.0);

    let mesh = extract_mesh(
        &shape,
        Vector3::from_element(-10.0),
        Vector3::from_element(10.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 50);
    export_both(&mesh, "15_smooth_difference");
}

// ===========================================================================
// Section 4: Transforms
// ===========================================================================

#[test]
#[ignore]
fn transform_translate() {
    let shape = Translate::new(
        Sphere::new(Vector3::zeros(), 5.0),
        Vector3::new(20.0, 0.0, 0.0),
    );
    // Original center should now be outside
    assert!(shape.evaluate(Vector3::zeros()) > 0.0);
    // Translated center should be inside
    assert!((shape.evaluate(Vector3::new(20.0, 0.0, 0.0)) - (-5.0)).abs() < 1e-10);

    let mesh = extract_mesh(
        &shape,
        Vector3::new(13.0, -7.0, -7.0),
        Vector3::new(27.0, 7.0, 7.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 50);
    export_both(&mesh, "16_translate");
}

#[test]
#[ignore]
fn transform_rotate() {
    // A box rotated 45 degrees around Y axis
    let shape = Rotate::new(
        Box3::new(Vector3::zeros(), Vector3::new(10.0, 5.0, 2.0)),
        Rotation3::from_axis_angle(&Unit::new_normalize(Vector3::y_axis().into_inner()), PI / 4.0),
    );
    // Center still inside
    assert!(shape.evaluate(Vector3::zeros()) < 0.0);
    // A point that was on the X face of the original box, now rotated
    // Original face at (10,0,0) -> rotated to (10*cos45, 0, -10*sin45)
    let rotated_face = Vector3::new(10.0 * (PI / 4.0).cos(), 0.0, -10.0 * (PI / 4.0).sin());
    assert!(shape.evaluate(rotated_face).abs() < 0.5);

    let mesh = extract_mesh(
        &shape,
        Vector3::new(-10.0, -7.0, -10.0),
        Vector3::new(10.0, 7.0, 10.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 20);
    export_both(&mesh, "17_rotate");
}

#[test]
#[ignore]
fn transform_scale() {
    let shape = Scale::new(
        Sphere::new(Vector3::zeros(), 5.0),
        3.0, // scale up 3x -> effective radius 15
    );
    // At distance 12, should be inside (15-12 = 3 units inside)
    assert!(shape.evaluate(Vector3::new(12.0, 0.0, 0.0)) < 0.0);
    // At distance 18, should be outside
    assert!(shape.evaluate(Vector3::new(18.0, 0.0, 0.0)) > 0.0);
    // Surface should be at 15
    assert!(shape.evaluate(Vector3::new(15.0, 0.0, 0.0)).abs() < 0.1);

    let mesh = extract_mesh(
        &shape,
        Vector3::from_element(-18.0),
        Vector3::from_element(18.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 100);
    export_both(&mesh, "18_scale");
}

#[test]
#[ignore]
fn transform_mirror() {
    // Capsule from (5,0,0) to (5,0,20), mirrored across YZ plane
    let shape = Mirror::new(
        Capsule::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(5.0, 0.0, 20.0),
            3.0,
        ),
        Vector3::new(1.0, 0.0, 0.0), // mirror across X=0 (YZ plane)
    );
    // Original position inside
    assert!(shape.evaluate(Vector3::new(5.0, 0.0, 10.0)) < 0.0);
    // Mirror position also inside
    assert!(shape.evaluate(Vector3::new(-5.0, 0.0, 10.0)) < 0.0);

    let mesh = extract_mesh(
        &shape,
        Vector3::new(-10.0, -5.0, -5.0),
        Vector3::new(10.0, 5.0, 25.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 50);
    export_both(&mesh, "19_mirror");
}

#[test]
#[ignore]
fn transform_shell() {
    // Hollow sphere: shell of a sphere
    let shape = Shell::new(Sphere::new(Vector3::zeros(), 10.0), 1.0);
    // At surface (r=10): |0| - 1 = -1, inside shell
    assert!(shape.evaluate(Vector3::new(10.0, 0.0, 0.0)) < 0.0);
    // At center (r=0): |(-10)| - 1 = 9, outside shell
    assert!(shape.evaluate(Vector3::zeros()) > 0.0);
    // At r=5: |(-5)| - 1 = 4, outside shell
    assert!(shape.evaluate(Vector3::new(5.0, 0.0, 0.0)) > 0.0);

    let mesh = extract_mesh(
        &shape,
        Vector3::from_element(-13.0),
        Vector3::from_element(13.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 100);
    export_both(&mesh, "20_shell");
}

// ===========================================================================
// Section 5: Complex Composed Shapes (Engineering-Grade Demonstrations)
// ===========================================================================

#[test]
#[ignore]
fn composed_mounting_bracket() {
    // A mounting bracket: base plate with two cylindrical mounting posts
    // and a large hole through the center

    // Base plate
    let plate = Box3::new(Vector3::zeros(), Vector3::new(20.0, 2.0, 15.0));

    // Left mounting post
    let post_left = Capsule::new(
        Vector3::new(-12.0, 2.0, 0.0),
        Vector3::new(-12.0, 15.0, 0.0),
        3.0,
    );

    // Right mounting post
    let post_right = Capsule::new(
        Vector3::new(12.0, 2.0, 0.0),
        Vector3::new(12.0, 15.0, 0.0),
        3.0,
    );

    // Center hole through plate
    let hole = Cylinder::new(
        Vector3::new(0.0, -5.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        5.0,
        14.0,
    );

    // Compose: (plate + posts) - hole
    let with_posts = Union::new(Union::new(plate, post_left), post_right);
    let bracket = Difference::new(with_posts, hole);

    // Inside plate, away from hole
    assert!(bracket.evaluate(Vector3::new(15.0, 0.0, 10.0)) < 0.0);
    // Inside hole
    assert!(bracket.evaluate(Vector3::new(0.0, 0.0, 0.0)) > 0.0);

    let mesh = extract_mesh(
        &bracket,
        Vector3::new(-18.0, -4.0, -18.0),
        Vector3::new(18.0, 18.0, 18.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 100);
    export_both(&mesh, "21_mounting_bracket");
}

#[test]
#[ignore]
fn composed_pipe_tee() {
    // T-shaped pipe fitting: two cylinders intersecting at right angles
    // with the interior hollowed out

    let main_pipe = Capsule::new(
        Vector3::new(-15.0, 0.0, 0.0),
        Vector3::new(15.0, 0.0, 0.0),
        5.0,
    );
    let branch = Capsule::new(
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 15.0, 0.0),
        5.0,
    );

    // Outer shell
    let outer = SmoothUnion::new(main_pipe, branch, 1.5);

    // Inner bore (smaller radius)
    let main_bore = Capsule::new(
        Vector3::new(-16.0, 0.0, 0.0),
        Vector3::new(16.0, 0.0, 0.0),
        3.0,
    );
    let branch_bore = Capsule::new(
        Vector3::new(0.0, -1.0, 0.0),
        Vector3::new(0.0, 16.0, 0.0),
        3.0,
    );
    let inner = Union::new(main_bore, branch_bore);

    // Hollow pipe: outer - inner
    let tee = Difference::new(outer, inner);

    // On the wall (between inner and outer radius)
    assert!(tee.evaluate(Vector3::new(10.0, 4.0, 0.0)) < 0.0);
    // Inside the bore
    assert!(tee.evaluate(Vector3::new(10.0, 0.0, 0.0)) > 0.0);

    let mesh = extract_mesh(
        &tee,
        Vector3::new(-18.0, -3.0, -8.0),
        Vector3::new(18.0, 18.0, 8.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 100);
    export_both(&mesh, "22_pipe_tee");
}

#[test]
#[ignore]
fn composed_gasket_ring() {
    // A gasket: torus with flat top and bottom (intersection with slab)

    let torus = Torus::new(Vector3::zeros(), 12.0, 4.0);
    let slab = Box3::new(Vector3::zeros(), Vector3::new(20.0, 2.0, 20.0));
    let gasket = Intersection::new(torus, slab);

    // Inside the gasket ring
    assert!(gasket.evaluate(Vector3::new(12.0, 0.0, 0.0)) < 0.0);
    // Center hole
    assert!(gasket.evaluate(Vector3::zeros()) > 0.0);
    // Above the slab
    assert!(gasket.evaluate(Vector3::new(12.0, 5.0, 0.0)) > 0.0);

    let mesh = extract_mesh(
        &gasket,
        Vector3::new(-18.0, -4.0, -18.0),
        Vector3::new(18.0, 4.0, 18.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 50);
    export_both(&mesh, "23_gasket_ring");
}

#[test]
#[ignore]
fn composed_rounded_enclosure() {
    // Electronics enclosure: rounded box with interior cavity and cable port

    let outer = RoundedBox::new(
        Vector3::zeros(),
        Vector3::new(15.0, 8.0, 10.0),
        2.0,
    );
    // Interior cavity (slightly smaller)
    let cavity = Box3::new(
        Vector3::new(0.0, 1.0, 0.0), // offset up to leave floor
        Vector3::new(13.0, 7.0, 8.0),
    );
    // Cable port on the side
    let cable_port = Cylinder::new(
        Vector3::new(15.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        3.0,
        10.0,
    );

    let enclosure = Difference::new(Difference::new(outer, cavity), cable_port);

    // In the wall
    assert!(enclosure.evaluate(Vector3::new(14.0, 0.0, 0.0)) < 0.0);
    // In the cavity
    assert!(enclosure.evaluate(Vector3::new(0.0, 3.0, 0.0)) > 0.0);

    let mesh = extract_mesh(
        &enclosure,
        Vector3::new(-19.0, -12.0, -14.0),
        Vector3::new(19.0, 12.0, 14.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 100);
    export_both(&mesh, "24_rounded_enclosure");
}

#[test]
#[ignore]
fn composed_smooth_organic_blob() {
    // Organic shape: multiple spheres smoothly blended
    // Demonstrates the SDF kernel's natural ability to create organic forms

    let s1 = Sphere::new(Vector3::new(0.0, 0.0, 0.0), 6.0);
    let s2 = Sphere::new(Vector3::new(5.0, 4.0, 0.0), 4.0);
    let s3 = Sphere::new(Vector3::new(-4.0, 5.0, 3.0), 3.5);
    let s4 = Sphere::new(Vector3::new(2.0, -3.0, 5.0), 3.0);

    let blob = SmoothUnion::new(
        SmoothUnion::new(s1, s2, 3.0),
        SmoothUnion::new(s3, s4, 3.0),
        3.0,
    );

    // Center should be well inside
    assert!(blob.evaluate(Vector3::zeros()) < -3.0);

    let mesh = extract_mesh(
        &blob,
        Vector3::from_element(-12.0),
        Vector3::from_element(12.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 200);
    export_both(&mesh, "25_organic_blob");
}

#[test]
#[ignore]
fn composed_rotated_multi_body() {
    // Multiple rotated boxes demonstrating transform composition

    let box1 = Box3::new(Vector3::zeros(), Vector3::new(8.0, 2.0, 2.0));
    let box2 = Rotate::new(
        Box3::new(Vector3::zeros(), Vector3::new(8.0, 2.0, 2.0)),
        Rotation3::from_axis_angle(&Unit::new_normalize(Vector3::new(0.0, 0.0, 1.0)), PI / 3.0),
    );
    let box3 = Rotate::new(
        Box3::new(Vector3::zeros(), Vector3::new(8.0, 2.0, 2.0)),
        Rotation3::from_axis_angle(&Unit::new_normalize(Vector3::new(0.0, 0.0, 1.0)), 2.0 * PI / 3.0),
    );

    let star = Union::new(Union::new(box1, box2), box3);

    // Center inside all three
    assert!(star.evaluate(Vector3::zeros()) < 0.0);

    let mesh = extract_mesh(
        &star,
        Vector3::from_element(-11.0),
        Vector3::from_element(11.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 50);
    export_both(&mesh, "26_rotated_star");
}

#[test]
#[ignore]
fn composed_bearing_housing() {
    // Bearing housing: outer cylinder with inner bore and snap ring groove

    // Outer housing
    let housing = Capsule::new(
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 20.0, 0.0),
        12.0,
    );

    // Inner bore
    let bore = Capsule::new(
        Vector3::new(0.0, -2.0, 0.0),
        Vector3::new(0.0, 22.0, 0.0),
        8.0,
    );

    // Snap ring groove (torus cut into inner wall)
    // Rotate torus to lie in XY plane for the groove
    let groove_rotated = Rotate::new(
        Torus::new(Vector3::zeros(), 8.0, 1.5),
        Rotation3::from_axis_angle(&Unit::new_normalize(Vector3::new(1.0, 0.0, 0.0)), PI / 2.0),
    );
    let groove_positioned = Translate::new(groove_rotated, Vector3::new(0.0, 10.0, 0.0));

    let housing_with_bore = Difference::new(housing, bore);
    let final_housing = Difference::new(housing_with_bore, groove_positioned);

    // In the wall between bore and outer
    assert!(final_housing.evaluate(Vector3::new(10.0, 10.0, 0.0)) < 0.0);
    // Inside the bore
    assert!(final_housing.evaluate(Vector3::new(0.0, 10.0, 0.0)) > 0.0);

    let mesh = extract_mesh(
        &final_housing,
        Vector3::new(-15.0, -3.0, -15.0),
        Vector3::new(15.0, 23.0, 15.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 100);
    export_both(&mesh, "27_bearing_housing");
}

// ===========================================================================
// Section 6: Distance Accuracy Validation
// ===========================================================================

#[test]
fn accuracy_sphere_distance() {
    let shape = Sphere::new(Vector3::new(1.0, 2.0, 3.0), 5.0);

    // Test 100 random-ish points and verify distance
    for i in 0..100 {
        let angle1 = (i as f64) * 0.1;
        let angle2 = (i as f64) * 0.07;
        let r = 2.0 + (i as f64) * 0.15;
        let point = Vector3::new(
            1.0 + r * angle1.cos() * angle2.cos(),
            2.0 + r * angle1.sin() * angle2.cos(),
            3.0 + r * angle2.sin(),
        );
        let expected = (point - Vector3::new(1.0, 2.0, 3.0)).norm() - 5.0;
        let actual = shape.evaluate(point);
        assert!(
            (actual - expected).abs() < 1e-10,
            "Sphere distance mismatch at {:?}: expected {}, got {}",
            point, expected, actual
        );
    }
}

#[test]
fn accuracy_torus_on_surface() {
    let shape = Torus::new(Vector3::zeros(), 10.0, 3.0);

    // Points on the surface of the torus (angle around major circle)
    for i in 0..36 {
        let theta = (i as f64) * PI / 18.0;
        // Point on outer surface: major_radius + minor_radius
        let point = Vector3::new(13.0 * theta.cos(), 0.0, 13.0 * theta.sin());
        let d = shape.evaluate(point);
        assert!(d.abs() < 0.1, "Torus surface point should be near zero, got {}", d);
    }
}

#[test]
fn accuracy_boolean_containment() {
    // Verify boolean containment property for many sample points
    let a = Sphere::new(Vector3::new(-2.0, 0.0, 0.0), 5.0);
    let b = Box3::new(Vector3::new(1.0, 0.0, 0.0), Vector3::new(4.0, 4.0, 4.0));

    let union_shape = Union::new(
        Sphere::new(Vector3::new(-2.0, 0.0, 0.0), 5.0),
        Box3::new(Vector3::new(1.0, 0.0, 0.0), Vector3::new(4.0, 4.0, 4.0)),
    );

    // Test grid of points
    for ix in -10..=10 {
        for iy in -8..=8 {
            for iz in -8..=8 {
                let p = Vector3::new(ix as f64, iy as f64, iz as f64);
                let da = a.evaluate(p);
                let db = b.evaluate(p);
                let du = union_shape.evaluate(p);

                // Union is inside iff at least one operand is inside
                let a_inside = da < 0.0;
                let b_inside = db < 0.0;
                let u_inside = du < 0.0;

                if a_inside || b_inside {
                    assert!(u_inside,
                        "Union should be inside at {:?} (a={}, b={}, u={})",
                        p, da, db, du);
                }
                // Note: the converse isn't exactly true due to SDF approximation
                // near the surface, so we only test one direction with a margin
            }
        }
    }
}

// ===========================================================================
// Section 7: Mesh Quality Checks
// ===========================================================================

#[test]
#[ignore]
fn mesh_quality_manifold_check() {
    // Verify basic mesh topology: every edge should appear in exactly 2 triangles
    // for a watertight mesh (best-effort with marching cubes)
    let shape = Sphere::new(Vector3::zeros(), 8.0);
    let mesh = extract_mesh(&shape, Vector3::from_element(-10.0), Vector3::from_element(10.0), 32);

    // Count edge usage
    let mut edge_counts: HashMap<(u32, u32), u32> = HashMap::new();
    for chunk in mesh.indices.chunks(3) {
        let tri = [chunk[0], chunk[1], chunk[2]];
        for i in 0..3 {
            let (a, b) = (tri[i], tri[(i + 1) % 3]);
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_counts.entry(key).or_insert(0) += 1;
        }
    }

    // For marching cubes, most edges should be shared by exactly 2 triangles
    let manifold_edges = edge_counts.values().filter(|&&c| c == 2).count();
    let total_edges = edge_counts.len();
    let manifold_ratio = manifold_edges as f64 / total_edges as f64;

    assert!(
        manifold_ratio > 0.95,
        "Manifold ratio {:.2}% is too low (expected >95%)",
        manifold_ratio * 100.0
    );
}

use std::collections::HashMap;

#[test]
#[ignore]
fn mesh_quality_normal_consistency() {
    // Verify all normals are unit length
    let shape = Torus::new(Vector3::zeros(), 8.0, 3.0);
    let mesh = extract_mesh(
        &shape,
        Vector3::new(-13.0, -5.0, -13.0),
        Vector3::new(13.0, 5.0, 13.0),
        32,
    );

    for (i, n) in mesh.normals.iter().enumerate() {
        let len = n.norm();
        assert!(
            (len - 1.0).abs() < 0.01,
            "Normal {} has non-unit length: {}",
            i, len
        );
    }
}

#[test]
#[ignore]
fn mesh_quality_triangle_area() {
    // Verify no degenerate (zero-area) triangles
    let shape = Box3::new(Vector3::zeros(), Vector3::new(5.0, 5.0, 5.0));
    let mesh = extract_mesh(
        &shape,
        Vector3::new(-7.0, -7.0, -7.0),
        Vector3::new(7.0, 7.0, 7.0),
        32,
    );

    let mut degenerate_count = 0;
    for chunk in mesh.indices.chunks(3) {
        let va = &mesh.vertices[chunk[0] as usize];
        let vb = &mesh.vertices[chunk[1] as usize];
        let vc = &mesh.vertices[chunk[2] as usize];
        let e1 = vb - va;
        let e2 = vc - va;
        let area = e1.cross(&e2).norm() * 0.5;
        if area < 1e-10 {
            degenerate_count += 1;
        }
    }

    let total = mesh.indices.len() / 3;
    let degenerate_pct = degenerate_count as f64 / total as f64 * 100.0;
    assert!(
        degenerate_pct < 5.0,
        "{}% degenerate triangles ({}/{})",
        degenerate_pct, degenerate_count, total
    );
}
