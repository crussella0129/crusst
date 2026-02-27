// Blend profiles for SDF-based filleting and chamfering.
//
// Each profile defines a cross-section curve in (d1, d2) space that replaces
// the sharp edge where two half-spaces meet.  The core function
// `blend_intersection` computes the signed distance for an intersection
// (convex edge) case; union and difference are derived by sign flips.

/// A blend profile describes the cross-section curve used to round or
/// chamfer the edge between two SDF half-spaces.
#[derive(Debug, Clone, PartialEq)]
pub enum BlendProfile {
    /// G1 tangent-continuous cubic Bezier blend.
    G1 { radius: f64 },
    /// G2 curvature-continuous circular arc blend (exact, closed-form).
    G2 { radius: f64 },
    /// G3 blend (placeholder radius; blend math falls back to G2 for now).
    G3 { radius: f64 },
    /// Chord-length specified blend.
    Chord { chord_length: f64 },
    /// Cycloidal blend.
    Cycloidal { radius: f64 },
    /// Parabolic blend.
    Parabolic { radius: f64 },
    /// Hyperbolic blend with asymptote offset.
    Hyperbolic { radius: f64, asymptote: f64 },
    /// Equal-distance chamfer (45-degree line from (-k,0) to (0,-k)).
    EqualChamfer { distance: f64 },
    /// Two-distance chamfer (different setbacks on each face).
    TwoDistChamfer { distance1: f64, distance2: f64 },
    /// Angle chamfer (one distance plus angle).
    AngleChamfer { distance: f64, angle_rad: f64 },
}

// ---------------------------------------------------------------------------
// Shorthand constructors
// ---------------------------------------------------------------------------

/// G1 tangent-continuous blend with the given radius.
pub fn g1(radius: f64) -> BlendProfile {
    BlendProfile::G1 { radius }
}

/// G2 curvature-continuous circular arc blend with the given radius.
pub fn g2(radius: f64) -> BlendProfile {
    BlendProfile::G2 { radius }
}

/// G3 blend with the given radius.
pub fn g3(radius: f64) -> BlendProfile {
    BlendProfile::G3 { radius }
}

/// Chord-length specified blend.
pub fn chord(chord_length: f64) -> BlendProfile {
    BlendProfile::Chord { chord_length }
}

/// Cycloidal blend with the given radius.
pub fn cycloidal(radius: f64) -> BlendProfile {
    BlendProfile::Cycloidal { radius }
}

/// Parabolic blend with the given radius.
pub fn parabolic(radius: f64) -> BlendProfile {
    BlendProfile::Parabolic { radius }
}

/// Hyperbolic blend with the given radius and asymptote offset.
pub fn hyperbolic(radius: f64, asymptote: f64) -> BlendProfile {
    BlendProfile::Hyperbolic { radius, asymptote }
}

/// Equal-distance chamfer.
pub fn equal_chamfer(distance: f64) -> BlendProfile {
    BlendProfile::EqualChamfer { distance }
}

/// Two-distance chamfer.
pub fn two_dist_chamfer(distance1: f64, distance2: f64) -> BlendProfile {
    BlendProfile::TwoDistChamfer {
        distance1,
        distance2,
    }
}

/// Angle chamfer.
pub fn angle_chamfer(distance: f64, angle_rad: f64) -> BlendProfile {
    BlendProfile::AngleChamfer {
        distance,
        angle_rad,
    }
}

// ---------------------------------------------------------------------------
// Methods on BlendProfile
// ---------------------------------------------------------------------------

impl BlendProfile {
    /// Effective blend radius for the profile.
    pub fn radius(&self) -> f64 {
        match self {
            BlendProfile::G1 { radius }
            | BlendProfile::G2 { radius }
            | BlendProfile::G3 { radius }
            | BlendProfile::Cycloidal { radius }
            | BlendProfile::Parabolic { radius } => *radius,
            BlendProfile::Chord { chord_length } => chord_length / std::f64::consts::SQRT_2,
            BlendProfile::Hyperbolic { radius, .. } => *radius,
            BlendProfile::EqualChamfer { distance } => *distance,
            BlendProfile::TwoDistChamfer {
                distance1,
                distance2,
            } => distance1.max(*distance2),
            BlendProfile::AngleChamfer { distance, .. } => *distance,
        }
    }

    /// Maximum deviation from the sharp (unblended) result.
    /// Useful for interval arithmetic bounds.
    pub fn max_deviation(&self) -> f64 {
        match self {
            BlendProfile::G1 { radius }
            | BlendProfile::G2 { radius }
            | BlendProfile::G3 { radius } => *radius,
            BlendProfile::Cycloidal { radius } => 0.36 * radius,
            BlendProfile::Parabolic { radius } => radius / 4.0,
            BlendProfile::Hyperbolic { asymptote, .. } => *asymptote,
            BlendProfile::EqualChamfer { distance } => distance * 0.293,
            BlendProfile::TwoDistChamfer {
                distance1,
                distance2,
            } => distance1.max(*distance2) * 0.293,
            BlendProfile::AngleChamfer { distance, .. } => distance * 0.293,
            BlendProfile::Chord { chord_length } => {
                // Effective radius is chord_length / sqrt(2); max_deviation
                // follows the G2 rule (equals the effective radius).
                chord_length / std::f64::consts::SQRT_2
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Core blend functions
// ---------------------------------------------------------------------------

/// Signed-distance blend for an **intersection** (convex-edge) case.
///
/// `d1` and `d2` are SDF values for the two half-spaces.  The blend
/// replaces the sharp edge `max(d1, d2)` with a curved profile of the
/// specified type.
pub fn blend_intersection(d1: f64, d2: f64, profile: &BlendProfile) -> f64 {
    match profile {
        BlendProfile::G2 { radius } => blend_intersection_g2(d1, d2, *radius),
        BlendProfile::EqualChamfer { distance } => blend_intersection_chamfer(d1, d2, *distance),
        BlendProfile::G1 { radius } => blend_intersection_g1(d1, d2, *radius),
        // Profiles without dedicated math fall back to G2 with the effective radius.
        other => blend_intersection_g2(d1, d2, other.radius()),
    }
}

/// Signed-distance blend for a **union** (concave-edge) case.
/// Derived from intersection by negating both inputs and the result.
pub fn blend_union(d1: f64, d2: f64, profile: &BlendProfile) -> f64 {
    -blend_intersection(-d1, -d2, profile)
}

/// Signed-distance blend for a **difference** (A minus B) case.
/// Derived from intersection by negating d2.
pub fn blend_difference(d1: f64, d2: f64, profile: &BlendProfile) -> f64 {
    blend_intersection(d1, -d2, profile)
}

// ---------------------------------------------------------------------------
// G2 circular arc blend (exact, closed-form)
// ---------------------------------------------------------------------------

/// Arc centered at `(-r, -r)` in `(d1, d2)` space with radius `r`.
fn blend_intersection_g2(d1: f64, d2: f64, r: f64) -> f64 {
    let sharp = d1.max(d2);
    let u = (d1 + r).max(0.0);
    let v = (d2 + r).max(0.0);
    let arc_sdf = (u * u + v * v).sqrt() - r;
    sharp.max(arc_sdf)
}

// ---------------------------------------------------------------------------
// Equal chamfer blend (exact, closed-form)
// ---------------------------------------------------------------------------

/// Line from `(-k, 0)` to `(0, -k)`.
fn blend_intersection_chamfer(d1: f64, d2: f64, k: f64) -> f64 {
    let sharp = d1.max(d2);
    let chamfer_sdf = (d1 + d2 + k) / std::f64::consts::SQRT_2;
    sharp.max(chamfer_sdf)
}

// ---------------------------------------------------------------------------
// G1 cubic Bezier blend (Newton iteration)
// ---------------------------------------------------------------------------

/// Cubic Bezier control points for a G1 tangent-continuous quarter-curve
/// from `(-r, 0)` to `(0, -r)`.
///
/// The kappa value 0.5523 gives tangent continuity at the endpoints.
const KAPPA: f64 = 0.5523;

/// Number of Newton iterations for nearest-point search.
const NEWTON_ITERS: usize = 5;

/// Cubic Bezier at parameter `t` given four 2D control points.
#[inline]
fn bezier_point(p0: [f64; 2], p1: [f64; 2], p2: [f64; 2], p3: [f64; 2], t: f64) -> [f64; 2] {
    let s = 1.0 - t;
    let s2 = s * s;
    let t2 = t * t;
    [
        s2 * s * p0[0] + 3.0 * s2 * t * p1[0] + 3.0 * s * t2 * p2[0] + t2 * t * p3[0],
        s2 * s * p0[1] + 3.0 * s2 * t * p1[1] + 3.0 * s * t2 * p2[1] + t2 * t * p3[1],
    ]
}

/// First derivative of cubic Bezier at parameter `t`.
#[inline]
fn bezier_deriv(p0: [f64; 2], p1: [f64; 2], p2: [f64; 2], p3: [f64; 2], t: f64) -> [f64; 2] {
    let s = 1.0 - t;
    [
        3.0 * (s * s * (p1[0] - p0[0]) + 2.0 * s * t * (p2[0] - p1[0]) + t * t * (p3[0] - p2[0])),
        3.0 * (s * s * (p1[1] - p0[1]) + 2.0 * s * t * (p2[1] - p1[1]) + t * t * (p3[1] - p2[1])),
    ]
}

/// Second derivative of cubic Bezier at parameter `t`.
#[inline]
fn bezier_deriv2(p0: [f64; 2], p1: [f64; 2], p2: [f64; 2], p3: [f64; 2], t: f64) -> [f64; 2] {
    let s = 1.0 - t;
    [
        6.0 * (s * (p2[0] - 2.0 * p1[0] + p0[0]) + t * (p3[0] - 2.0 * p2[0] + p1[0])),
        6.0 * (s * (p2[1] - 2.0 * p1[1] + p0[1]) + t * (p3[1] - 2.0 * p2[1] + p1[1])),
    ]
}

#[inline]
fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

#[inline]
fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn blend_intersection_g1(d1: f64, d2: f64, r: f64) -> f64 {
    let sharp = d1.max(d2);

    // Only blend in the neighbourhood of the edge.
    // If both distances are far outside the blend zone, return sharp.
    if d1 > r || d2 > r {
        return sharp;
    }
    // If deep inside both faces, the blend cannot produce a larger value.
    if d1 < -r && d2 < -r {
        return sharp;
    }

    // Control points for the G1 Bezier in (d1, d2) space.
    let p0 = [-r, 0.0];
    let p1 = [-r, -r * KAPPA];
    let p2 = [-r * KAPPA, -r];
    let p3 = [0.0, -r];

    let q = [d1, d2];

    // Initial guess: linear interpolation based on which face is closer.
    let t_init = if (d1 + r).abs() + (d2 + r).abs() < 1e-12 {
        0.5
    } else {
        // Map query position to a reasonable starting t.
        ((d2 + r) / ((d1 + r) + (d2 + r)).max(1e-12)).clamp(0.0, 1.0)
    };

    let mut t = t_init;

    // Newton iteration to find nearest point on the Bezier.
    for _ in 0..NEWTON_ITERS {
        let b = bezier_point(p0, p1, p2, p3, t);
        let bp = bezier_deriv(p0, p1, p2, p3, t);
        let bpp = bezier_deriv2(p0, p1, p2, p3, t);

        let diff = [b[0] - q[0], b[1] - q[1]];
        let num = dot2(diff, bp);
        let den = dot2(bp, bp) + dot2(diff, bpp);

        if den.abs() < 1e-15 {
            break;
        }
        t = (t - num / den).clamp(0.0, 1.0);
    }

    // Signed distance from query to nearest point on the Bezier.
    let nearest = bezier_point(p0, p1, p2, p3, t);
    let tangent = bezier_deriv(p0, p1, p2, p3, t);
    let diff = [q[0] - nearest[0], q[1] - nearest[1]];
    let dist = (diff[0] * diff[0] + diff[1] * diff[1]).sqrt();

    // Sign: positive if outside the blend curve (toward the sharp edge),
    // negative if inside.  The cross product of the tangent with the
    // difference vector gives the orientation.
    let sign = cross2(tangent, diff);
    let signed_dist = if sign >= 0.0 { dist } else { -dist };

    sharp.max(signed_dist)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shorthand_constructors_roundtrip() {
        assert_eq!(g1(2.0).radius(), 2.0);
        assert_eq!(g2(3.0).radius(), 3.0);
        assert_eq!(g3(1.5).radius(), 1.5);
        assert_eq!(equal_chamfer(4.0).radius(), 4.0);
    }
}
