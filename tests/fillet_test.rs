use approx::assert_relative_eq;
use crusst::blend;
use crusst::builder::Shape;
use crusst::feature::ft;
use nalgebra::Vector3;

/// A box with G2 fillet on ALL edges: the corner should be carved away.
#[test]
fn filleted_box_corner_removed() {
    let b = Shape::box3(5.0, 5.0, 5.0).fillet(blend::g2(1.0), vec![ft(0, 0).all_edges()]);
    // At the sharp corner (5, 5, 5), the sharp SDF is 0.
    // With fillet, the corner is carved away, so SDF should be positive (outside).
    let d = b.distance(Vector3::new(5.0, 5.0, 5.0));
    assert!(d > 0.0, "filleted corner should be outside: got {}", d);
}

/// Fillet at the face center: should match the sharp result (on surface).
#[test]
fn fillet_face_center_unchanged() {
    let b = Shape::box3(5.0, 5.0, 5.0).fillet(blend::g2(1.0), vec![ft(0, 0).all_edges()]);
    // At (5, 0, 0) -- on the +X face center, far from any edge -- should be ~0 (on surface).
    let d = b.distance(Vector3::new(5.0, 0.0, 0.0));
    assert_relative_eq!(d, 0.0, epsilon = 0.1);
}

/// Chamfer on all edges removes the corner.
#[test]
fn chamfered_box_corner_removed() {
    let b =
        Shape::box3(5.0, 5.0, 5.0).chamfer(blend::equal_chamfer(1.0), vec![ft(0, 0).all_edges()]);
    let d = b.distance(Vector3::new(5.0, 5.0, 5.0));
    assert!(d > 0.0, "chamfered corner should be outside: got {}", d);
}

/// Selective fillet: only fillet edges on the +X face.
#[test]
fn selective_fillet_top_edges() {
    // Edge indices for +X face: edges 0 (+X/+Y), 1 (+X/-Y), 2 (+X/+Z), 3 (+X/-Z)
    let b = Shape::box3(5.0, 5.0, 5.0).fillet(blend::g2(1.0), vec![ft(0, 0).edges(&[0, 1, 2, 3])]);

    // +X/+Y edge at (5, 5, 0) should be filleted (edge 0)
    let d_filleted = b.distance(Vector3::new(5.0, 5.0, 0.0));
    assert!(
        d_filleted > -0.1,
        "+X edge should be filleted: got {}",
        d_filleted
    );

    // -X/+Y edge at (-5, 5, 0) (edge 4: -X/+Y) should remain sharp
    let d_sharp = b.distance(Vector3::new(-5.0, -5.0, 0.0));
    assert!(
        d_sharp < 0.05,
        "-X/-Y edge should remain sharp: got {}",
        d_sharp
    );
}

/// Interior point should remain inside after fillet.
#[test]
fn fillet_interior_point_still_inside() {
    let b = Shape::box3(5.0, 5.0, 5.0).fillet(blend::g2(1.0), vec![ft(0, 0).all_edges()]);
    // The origin is the box center, should be well inside.
    let d = b.distance(Vector3::zeros());
    assert!(d < -4.0, "origin should be inside filleted box: got {}", d);
}
