use approx::assert_relative_eq;
use crusst::blend::{
    angle_chamfer, blend_difference, blend_intersection, blend_union, chord, cycloidal,
    equal_chamfer, g1, g2, g3, hyperbolic, parabolic, two_dist_chamfer,
};

// ---------------------------------------------------------------------------
// G2 circular arc
// ---------------------------------------------------------------------------

/// G2 blend at the tangent point on face 1: (d1, d2) = (-r, 0).
/// The arc meets face 1 here, so the SDF should be exactly 0.
#[test]
fn g2_tangent_point_d1() {
    let r = 1.0;
    let profile = g2(r);
    let result = blend_intersection(-r, 0.0, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-6);
}

/// G2 blend at the tangent point on face 2: (d1, d2) = (0, -r).
/// By symmetry with the d1 case, the SDF should be exactly 0.
#[test]
fn g2_tangent_point_d2() {
    let r = 1.0;
    let profile = g2(r);
    let result = blend_intersection(0.0, -r, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-6);
}

/// G2 blend at the arc midpoint on the surface (45-degree point).
/// The midpoint of the arc lies at (-r + r/sqrt(2), -r + r/sqrt(2)).
/// The SDF should be 0 because this point is on the blend surface.
#[test]
fn g2_surface_midpoint() {
    let r = 1.0;
    let m = -r + r / std::f64::consts::SQRT_2;
    let profile = g2(r);
    let result = blend_intersection(m, m, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-6);
}

/// G2 blend far outside one face: (d1, d2) = (5r, -5r).
/// When one distance is deep inside and the other far outside, the blend
/// should match the sharp result max(d1, d2).
#[test]
fn g2_far_outside() {
    let r = 1.0;
    let profile = g2(r);
    let result = blend_intersection(5.0 * r, -5.0 * r, &profile);
    let sharp = (5.0_f64 * r).max(-5.0 * r);
    assert_relative_eq!(result, sharp, epsilon = 1e-6);
}

/// G2 blend deep inside both faces: (d1, d2) = (-5r, -5r).
/// Deep inside, the arc SDF clamps to -r (the arc center depth).
/// Result = max(sharp, arc_sdf) = max(-5r, -r) = -r.
#[test]
fn g2_deep_inside() {
    let r = 1.0;
    let profile = g2(r);
    let result = blend_intersection(-5.0 * r, -5.0 * r, &profile);
    assert_relative_eq!(result, -r, epsilon = 1e-6);
}

// ---------------------------------------------------------------------------
// Equal chamfer
// ---------------------------------------------------------------------------

/// EqualChamfer at the tangent point on face 1: (d1, d2) = (-k, 0).
/// The chamfer line meets face 1 here, so SDF = 0.
#[test]
fn chamfer_tangent_d1() {
    let k = 1.0;
    let profile = equal_chamfer(k);
    let result = blend_intersection(-k, 0.0, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-6);
}

/// EqualChamfer at the tangent point on face 2: (d1, d2) = (0, -k).
/// By symmetry, SDF = 0.
#[test]
fn chamfer_tangent_d2() {
    let k = 1.0;
    let profile = equal_chamfer(k);
    let result = blend_intersection(0.0, -k, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-6);
}

/// EqualChamfer at the midpoint of the chamfer line: (-k/2, -k/2).
/// This point lies exactly on the chamfer, so SDF = 0.
#[test]
fn chamfer_midpoint() {
    let k = 1.0;
    let profile = equal_chamfer(k);
    let result = blend_intersection(-k / 2.0, -k / 2.0, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-6);
}

// ---------------------------------------------------------------------------
// G1 cubic Bezier (Newton iteration)
// ---------------------------------------------------------------------------

/// G1 blend at the tangent point on face 1: (d1, d2) = (-r, 0).
/// The Bezier starts here, so the SDF should be ~0.
#[test]
fn g1_tangent_d1() {
    let r = 1.0;
    let profile = g1(r);
    let result = blend_intersection(-r, 0.0, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-4);
}

/// G1 blend at the tangent point on face 2: (d1, d2) = (0, -r).
/// The Bezier ends here, so the SDF should be ~0.
#[test]
fn g1_tangent_d2() {
    let r = 1.0;
    let profile = g1(r);
    let result = blend_intersection(0.0, -r, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-4);
}

// ---------------------------------------------------------------------------
// blend_union symmetry
// ---------------------------------------------------------------------------

/// blend_union(d1, d2) should equal -blend_intersection(-d1, -d2) by
/// definition.  Verify with a concrete query point.
#[test]
fn blend_union_symmetry() {
    let r = 1.0;
    let profile = g2(r);
    let d1 = 0.3;
    let d2 = -0.5;
    let u = blend_union(d1, d2, &profile);
    let i = blend_intersection(-d1, -d2, &profile);
    assert_relative_eq!(u, -i, epsilon = 1e-12);

    // Also verify blend_difference is blend_intersection(d1, -d2).
    let diff = blend_difference(d1, d2, &profile);
    let i2 = blend_intersection(d1, -d2, &profile);
    assert_relative_eq!(diff, i2, epsilon = 1e-12);
}

// ---------------------------------------------------------------------------
// G3 quintic smoothstep (Newton iteration)
// ---------------------------------------------------------------------------

/// G3 blend at the tangent point on face 1: (d1, d2) = (-r, 0).
/// The quintic curve starts here, so SDF should be ~0.
#[test]
fn g3_tangent_d1() {
    let r = 1.0;
    let profile = g3(r);
    let result = blend_intersection(-r, 0.0, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-4);
}

/// G3 blend at the tangent point on face 2: (d1, d2) = (0, -r).
/// The quintic curve ends here, so SDF should be ~0.
#[test]
fn g3_tangent_d2() {
    let r = 1.0;
    let profile = g3(r);
    let result = blend_intersection(0.0, -r, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-4);
}

// ---------------------------------------------------------------------------
// Two-distance chamfer (closed-form)
// ---------------------------------------------------------------------------

/// TwoDistChamfer at the tangent point on face 1: (d1, d2) = (-k1, 0).
/// The chamfer line meets face 1 here, so SDF = 0.
#[test]
fn two_dist_chamfer_tangent_d1() {
    let k1 = 2.0;
    let k2 = 1.0;
    let profile = two_dist_chamfer(k1, k2);
    let result = blend_intersection(-k1, 0.0, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-6);
}

/// TwoDistChamfer at the tangent point on face 2: (d1, d2) = (0, -k2).
/// The chamfer line meets face 2 here, so SDF = 0.
#[test]
fn two_dist_chamfer_tangent_d2() {
    let k1 = 2.0;
    let k2 = 1.0;
    let profile = two_dist_chamfer(k1, k2);
    let result = blend_intersection(0.0, -k2, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-6);
}

// ---------------------------------------------------------------------------
// Angle chamfer (delegates to two-distance)
// ---------------------------------------------------------------------------

/// AngleChamfer at the tangent point on face 1: (d1, d2) = (-d, 0).
/// The chamfer line meets face 1 here, so SDF = 0.
#[test]
fn angle_chamfer_tangent_d1() {
    let d = 1.0;
    let angle = std::f64::consts::FRAC_PI_4; // 45 degrees => k2 = d * tan(45) = d
    let profile = angle_chamfer(d, angle);
    let result = blend_intersection(-d, 0.0, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-6);
}

// ---------------------------------------------------------------------------
// Parabolic blend (Newton iteration)
// ---------------------------------------------------------------------------

/// Parabolic blend at the tangent point on face 1: (d1, d2) = (-r, 0).
/// The parabola starts here, so SDF should be ~0.
#[test]
fn parabolic_tangent_d1() {
    let r = 1.0;
    let profile = parabolic(r);
    let result = blend_intersection(-r, 0.0, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-4);
}

/// Parabolic blend at the tangent point on face 2: (d1, d2) = (0, -r).
/// The parabola ends here, so SDF should be ~0.
#[test]
fn parabolic_tangent_d2() {
    let r = 1.0;
    let profile = parabolic(r);
    let result = blend_intersection(0.0, -r, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-4);
}

// ---------------------------------------------------------------------------
// Cycloidal blend (Newton iteration)
// ---------------------------------------------------------------------------

/// Cycloidal blend at the tangent point on face 1: (d1, d2) = (-r, 0).
/// The cycloid starts here, so SDF should be ~0.
#[test]
fn cycloidal_tangent_d1() {
    let r = 1.0;
    let profile = cycloidal(r);
    let result = blend_intersection(-r, 0.0, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-4);
}

// ---------------------------------------------------------------------------
// Hyperbolic blend (Newton iteration)
// ---------------------------------------------------------------------------

/// Hyperbolic blend at the tangent point on face 1: (d1, d2) = (-r, 0).
/// The curve starts here, so SDF should be ~0.
#[test]
fn hyperbolic_tangent_d1() {
    let r = 1.0;
    let a = 0.5; // asymptote parameter
    let profile = hyperbolic(r, a);
    let result = blend_intersection(-r, 0.0, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-4);
}

// ---------------------------------------------------------------------------
// Chord blend (delegates to G2)
// ---------------------------------------------------------------------------

/// Chord blend at the tangent point on face 1.
/// Chord delegates to G2 with r = chord_length / sqrt(2),
/// so at d1 = -r the SDF should be ~0.
#[test]
fn chord_tangent_d1() {
    let cl = 2.0;
    let r = cl / std::f64::consts::SQRT_2;
    let profile = chord(cl);
    let result = blend_intersection(-r, 0.0, &profile);
    assert_relative_eq!(result, 0.0, epsilon = 1e-6);
}
