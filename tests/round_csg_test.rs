use crusst::builder::Shape;
use nalgebra::Vector3;

#[test]
fn round_union_is_continuous() {
    // Two overlapping spheres with round union
    let a = Shape::sphere(5.0);
    let b = Shape::sphere(5.0).translate(6.0, 0.0, 0.0);
    let shape = a.round_union(b, 1.0);
    // Point at the midpoint between sphere centers should be inside
    let d = shape.distance(Vector3::new(3.0, 0.0, 0.0));
    assert!(d < 0.0, "midpoint should be inside round union, got {d}");
}

#[test]
fn round_intersection_fillet() {
    // Two half-spaces meeting at a right angle
    // HalfSpace solid region: normal.dot(p) + d <= 0
    // For x <= 2: normal=(1,0,0), d=-2 => x - 2 <= 0
    // For y <= 2: normal=(0,1,0), d=-2 => y - 2 <= 0
    let a = Shape::half_space(Vector3::new(1.0, 0.0, 0.0), -2.0); // x <= 2
    let b = Shape::half_space(Vector3::new(0.0, 1.0, 0.0), -2.0); // y <= 2
    let shape = a.round_intersect(b, 1.0);
    // At the origin (well inside both half-spaces), should be inside
    let d = shape.distance(Vector3::new(0.0, 0.0, 0.0));
    assert!(d < 0.0, "origin should be inside intersection, got {d}");
    // Near the corner (1.5, 1.5): should still be inside the rounded intersection
    let d2 = shape.distance(Vector3::new(1.5, 1.5, 0.0));
    assert!(d2 < 0.0, "near corner should be inside, got {d2}");
}

#[test]
fn chamfer_union_cuts_corner() {
    let a = Shape::sphere(5.0);
    let b = Shape::sphere(5.0).translate(6.0, 0.0, 0.0);
    let sharp = a.clone().union(b.clone());
    let chamfered = a.chamfer_union(b, 1.0);
    // Chamfer removes less material than round at the same size
    let d_sharp = sharp.distance(Vector3::new(3.0, 0.0, 0.0));
    let d_chamfer = chamfered.distance(Vector3::new(3.0, 0.0, 0.0));
    assert!(
        d_chamfer <= d_sharp + 0.01,
        "chamfer should be at least as inside as sharp: chamfer={d_chamfer}, sharp={d_sharp}"
    );
}

#[test]
fn round_method_on_builder() {
    let shape = Shape::box3(5.0, 5.0, 5.0).round(1.0);
    // At (6,6,6) for a box with he=(5,5,5), d_box = sqrt(1+1+1) = sqrt(3)
    // With round(1.0): d_rounded = d_box - 1.0 = sqrt(3) - 1.0
    let d_corner = shape.distance(Vector3::new(6.0, 6.0, 6.0));
    let expected = (3.0_f64).sqrt() - 1.0; // sqrt(3) - 1 ≈ 0.73
    assert!(
        (d_corner - expected).abs() < 0.01,
        "rounded corner distance: got {d_corner}, expected {expected}"
    );
}
