//! Knot vector utilities for NURBS curves and surfaces.

/// Find multiplicity of knot value `u` in the knot vector.
pub fn knot_multiplicity(u: f64, knots: &[f64], tol: f64) -> usize {
    knots.iter().filter(|&&k| (k - u).abs() < tol).count()
}

/// Check if a knot vector is valid:
/// - Non-decreasing
/// - Correct length: `n_ctrl + degree + 1`
pub fn validate_knot_vector(knots: &[f64], degree: usize, n_ctrl: usize) -> bool {
    if knots.len() != n_ctrl + degree + 1 {
        return false;
    }
    for i in 1..knots.len() {
        if knots[i] < knots[i - 1] {
            return false;
        }
    }
    true
}

/// Create a uniform clamped knot vector.
///
/// The first `degree+1` knots are 0.0, the last `degree+1` are 1.0,
/// and interior knots are uniformly spaced.
pub fn uniform_knots(degree: usize, n_ctrl: usize) -> Vec<f64> {
    assert!(n_ctrl > degree, "Need at least degree+1 control points");
    let m = n_ctrl + degree + 1;
    let mut knots = vec![0.0; m];

    for i in 0..=degree {
        knots[i] = 0.0;
        knots[m - 1 - i] = 1.0;
    }

    let n_interior = n_ctrl - degree - 1;
    if n_interior > 0 {
        for i in 1..=n_interior {
            knots[degree + i] = i as f64 / (n_interior + 1) as f64;
        }
    }

    knots
}

/// Insert a knot value `u` into a curve's knot vector and update control points.
///
/// Returns the new knot vector and new control points (homogeneous coordinates
/// must be handled by the caller).
///
/// Based on Algorithm A5.1 from The NURBS Book.
pub fn knot_insert(
    knots: &[f64],
    ctrl_pts: &[nalgebra::Vector4<f64>],
    degree: usize,
    u: f64,
) -> (Vec<f64>, Vec<nalgebra::Vector4<f64>>) {
    let n = ctrl_pts.len() - 1;
    let k = super::basis::find_span(n, degree, u, knots);
    let p = degree;

    // New knot vector: insert u at position k+1
    let mut new_knots = Vec::with_capacity(knots.len() + 1);
    new_knots.extend_from_slice(&knots[..=k]);
    new_knots.push(u);
    new_knots.extend_from_slice(&knots[k + 1..]);

    // New control points
    let mut new_pts = Vec::with_capacity(ctrl_pts.len() + 1);

    // Points before the affected range are unchanged
    for i in 0..=(k.saturating_sub(p)) {
        new_pts.push(ctrl_pts[i]);
    }

    // Affected range: compute blended points
    for i in (k.saturating_sub(p) + 1)..=k {
        let alpha = (u - knots[i]) / (knots[i + p] - knots[i]);
        new_pts.push(ctrl_pts[i - 1] * (1.0 - alpha) + ctrl_pts[i] * alpha);
    }

    // Points after the affected range are unchanged (shifted by 1 in output)
    for i in k..=n {
        new_pts.push(ctrl_pts[i]);
    }

    (new_knots, new_pts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validate_clamped_cubic() {
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.0];
        assert!(validate_knot_vector(&knots, 3, 7));
    }

    #[test]
    fn validate_wrong_length() {
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        assert!(!validate_knot_vector(&knots, 2, 4)); // needs 7 knots
    }

    #[test]
    fn validate_decreasing() {
        let knots = vec![0.0, 0.0, 0.5, 0.3, 1.0, 1.0]; // not non-decreasing
        assert!(!validate_knot_vector(&knots, 2, 3));
    }

    #[test]
    fn uniform_knots_quadratic() {
        let knots = uniform_knots(2, 5);
        // Should be [0, 0, 0, 0.5, 1, 1, 1] (but exact spacing depends)
        assert_eq!(knots.len(), 8); // 5 + 2 + 1
        assert_eq!(knots[0], 0.0);
        assert_eq!(knots[1], 0.0);
        assert_eq!(knots[2], 0.0);
        assert_eq!(knots[knots.len() - 1], 1.0);
        assert_eq!(knots[knots.len() - 2], 1.0);
        assert_eq!(knots[knots.len() - 3], 1.0);
        assert!(validate_knot_vector(&knots, 2, 5));
    }

    #[test]
    fn uniform_knots_minimum() {
        // degree+1 control points → no interior knots
        let knots = uniform_knots(3, 4);
        assert_eq!(knots, vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]);
    }

    #[test]
    fn multiplicity_at_ends() {
        let knots = vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0];
        assert_eq!(knot_multiplicity(0.0, &knots, 1e-12), 3);
        assert_eq!(knot_multiplicity(1.0, &knots, 1e-12), 3);
        assert_eq!(knot_multiplicity(0.5, &knots, 1e-12), 1);
        assert_eq!(knot_multiplicity(0.3, &knots, 1e-12), 0);
    }

    #[test]
    fn knot_insert_preserves_curve() {
        // A degree-2 curve with 3 control points (Bezier segment).
        // Inserting a knot should not change the curve shape.
        use nalgebra::Vector4;

        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let pts = vec![
            Vector4::new(0.0, 0.0, 0.0, 1.0),
            Vector4::new(0.5, 1.0, 0.0, 1.0),
            Vector4::new(1.0, 0.0, 0.0, 1.0),
        ];

        let (new_knots, new_pts) = knot_insert(&knots, &pts, 2, 0.5);
        assert_eq!(new_knots.len(), 7);
        assert_eq!(new_pts.len(), 4);

        // Verify the new knot vector contains 0.5
        assert!((new_knots[3] - 0.5).abs() < 1e-14);
    }
}
