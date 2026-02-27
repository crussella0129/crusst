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
        BlendProfile::G3 { radius } => blend_intersection_g3(d1, d2, *radius),
        BlendProfile::Chord { chord_length } => {
            blend_intersection_g2(d1, d2, chord_length / std::f64::consts::SQRT_2)
        }
        BlendProfile::TwoDistChamfer {
            distance1,
            distance2,
        } => blend_intersection_two_dist_chamfer(d1, d2, *distance1, *distance2),
        BlendProfile::AngleChamfer {
            distance,
            angle_rad,
        } => blend_intersection_angle_chamfer(d1, d2, *distance, *angle_rad),
        BlendProfile::Parabolic { radius } => blend_intersection_parabolic(d1, d2, *radius),
        BlendProfile::Cycloidal { radius } => blend_intersection_cycloidal(d1, d2, *radius),
        BlendProfile::Hyperbolic { radius, asymptote } => {
            blend_intersection_hyperbolic(d1, d2, *radius, *asymptote)
        }
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

    // Control points for the G1 Bezier quarter-circle approximation in (d1, d2)
    // space.  The arc goes from (-r, 0) to (0, -r) and is centered at (-r, -r).
    // In local coords (centered at (-r,-r)), the standard Bezier quarter-circle
    // from (0, r) to (r, 0) has control points (κr, r), (r, κr).
    // Translating back: p1 = (-r + κr, 0), p2 = (0, -r + κr).
    let p0 = [-r, 0.0];
    let p1 = [-r + r * KAPPA, 0.0];
    let p2 = [0.0, -r + r * KAPPA];
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

// ---------------------------------------------------------------------------
// Sampling-based initial guess for Newton iteration
// ---------------------------------------------------------------------------

/// Find the best starting parameter by sampling `N` equally-spaced points
/// along a parametric curve and returning the t that minimizes distance
/// squared to `q`.
fn best_initial_t<F: Fn(f64) -> [f64; 2]>(q: [f64; 2], curve: &F, n: usize) -> f64 {
    let mut best_t = 0.5;
    let mut best_d2 = f64::MAX;
    for i in 0..=n {
        let t = i as f64 / n as f64;
        let p = curve(t);
        let dx = p[0] - q[0];
        let dy = p[1] - q[1];
        let d2 = dx * dx + dy * dy;
        if d2 < best_d2 {
            best_d2 = d2;
            best_t = t;
        }
    }
    best_t
}

// ---------------------------------------------------------------------------
// G3 quintic smoothstep blend (Newton iteration)
// ---------------------------------------------------------------------------

/// Quintic Hermite (smootherstep) parametric curve from `(-r, 0)` to `(0, -r)`.
///
/// Uses `f(t) = 6t^5 - 15t^4 + 10t^3` so that f(0)=0, f(1)=1 with first
/// and second derivatives zero at endpoints, giving G3 (curvature-rate)
/// continuity with both faces.
fn blend_intersection_g3(d1: f64, d2: f64, r: f64) -> f64 {
    let sharp = d1.max(d2);
    if d1 > r || d2 > r {
        return sharp;
    }
    if d1 < -r && d2 < -r {
        return sharp;
    }

    let q = [d1, d2];

    // Quintic smoothstep: f(t) = 6t^5 - 15t^4 + 10t^3
    // f'(t) = 30t^4 - 60t^3 + 30t^2 = 30t^2(t-1)^2
    // f''(t) = 120t^3 - 180t^2 + 60t = 60t(2t^2 - 3t + 1) = 60t(2t-1)(t-1)
    //
    // Curve: P(t) = (-r + r*f(t), -r + r*f(1-t))
    // P(0) = (-r, -r + r*1) = (-r, 0)
    // P(1) = (-r + r*1, -r + r*0) = (0, -r)
    #[inline]
    fn smootherstep(t: f64) -> f64 {
        t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
    }
    #[inline]
    fn smootherstep_d(t: f64) -> f64 {
        30.0 * t * t * (1.0 - t) * (1.0 - t)
    }
    #[inline]
    fn smootherstep_dd(t: f64) -> f64 {
        60.0 * t * (2.0 * t * t - 3.0 * t + 1.0)
    }

    let curve = |t: f64| -> [f64; 2] {
        let s = 1.0 - t;
        [-r + r * smootherstep(t), -r + r * smootherstep(s)]
    };

    let mut t = best_initial_t(q, &curve, 16);

    for _ in 0..NEWTON_ITERS {
        let s = 1.0 - t;
        let px = -r + r * smootherstep(t);
        let py = -r + r * smootherstep(s);
        let dpx = r * smootherstep_d(t);
        let dpy = -r * smootherstep_d(s); // d/dt of r*f(1-t) = -r*f'(1-t)
        let d2px = r * smootherstep_dd(t);
        let d2py = r * smootherstep_dd(s); // d2/dt2 of r*f(1-t) = r*f''(1-t)

        let diff = [px - q[0], py - q[1]];
        let num = diff[0] * dpx + diff[1] * dpy;
        let den = dpx * dpx + dpy * dpy + diff[0] * d2px + diff[1] * d2py;
        if den.abs() < 1e-15 {
            break;
        }
        t = (t - num / den).clamp(0.0, 1.0);
    }

    let s = 1.0 - t;
    let nearest = [-r + r * smootherstep(t), -r + r * smootherstep(s)];
    let tangent = [r * smootherstep_d(t), -r * smootherstep_d(s)];
    let diff = [q[0] - nearest[0], q[1] - nearest[1]];
    let dist = (diff[0] * diff[0] + diff[1] * diff[1]).sqrt();

    let sign = cross2(tangent, diff);
    let signed_dist = if sign >= 0.0 { dist } else { -dist };
    sharp.max(signed_dist)
}

// ---------------------------------------------------------------------------
// Two-distance chamfer (exact, closed-form)
// ---------------------------------------------------------------------------

/// Line from `(-k1, 0)` to `(0, -k2)` in `(d1, d2)` space.
///
/// The line equation is `d1/k1 + d2/k2 = -1`, or `d1*k2 + d2*k1 + k1*k2 = 0`.
/// Signed distance from a point to this line gives the chamfer SDF.
fn blend_intersection_two_dist_chamfer(d1: f64, d2: f64, k1: f64, k2: f64) -> f64 {
    let sharp = d1.max(d2);
    let chamfer_sdf = (d1 * k2 + d2 * k1 + k1 * k2) / (k1 * k1 + k2 * k2).sqrt();
    sharp.max(chamfer_sdf)
}

// ---------------------------------------------------------------------------
// Angle chamfer (delegates to two-distance chamfer)
// ---------------------------------------------------------------------------

/// Converts angle + distance to two distances, then delegates.
///
/// `distance` is the setback along face 1, `angle_rad` is the angle from
/// face 1 toward face 2.  The setback along face 2 is `distance * tan(angle)`.
fn blend_intersection_angle_chamfer(d1: f64, d2: f64, distance: f64, angle_rad: f64) -> f64 {
    let k1 = distance;
    let k2 = distance * angle_rad.tan();
    blend_intersection_two_dist_chamfer(d1, d2, k1, k2)
}

// ---------------------------------------------------------------------------
// Parabolic blend (Newton iteration)
// ---------------------------------------------------------------------------

/// Parabolic arc from `(-r, 0)` to `(0, -r)`.
///
/// Parametric curve: `P(t) = (-r + r*t, -r*t^2)` for `t in [0, 1]`.
/// At t=0: (-r, 0).  At t=1: (0, -r).
fn blend_intersection_parabolic(d1: f64, d2: f64, r: f64) -> f64 {
    let sharp = d1.max(d2);
    if d1 > r || d2 > r {
        return sharp;
    }
    if d1 < -r && d2 < -r {
        return sharp;
    }

    let q = [d1, d2];
    let curve = |t: f64| -> [f64; 2] { [-r + r * t, -r * t * t] };
    let mut t = best_initial_t(q, &curve, 16);
    for _ in 0..NEWTON_ITERS {
        let px = -r + r * t;
        let py = -r * t * t;
        let dpx = r;
        let dpy = -2.0 * r * t;
        let d2px = 0.0;
        let d2py = -2.0 * r;

        let diff = [px - q[0], py - q[1]];
        let num = diff[0] * dpx + diff[1] * dpy;
        let den = dpx * dpx + dpy * dpy + diff[0] * d2px + diff[1] * d2py;
        if den.abs() < 1e-15 {
            break;
        }
        t = (t - num / den).clamp(0.0, 1.0);
    }

    let nearest = [-r + r * t, -r * t * t];
    let tangent = [r, -2.0 * r * t];
    let diff = [q[0] - nearest[0], q[1] - nearest[1]];
    let dist = (diff[0] * diff[0] + diff[1] * diff[1]).sqrt();

    let sign = cross2(tangent, diff);
    let signed_dist = if sign >= 0.0 { dist } else { -dist };
    sharp.max(signed_dist)
}

// ---------------------------------------------------------------------------
// Cycloidal blend (Newton iteration)
// ---------------------------------------------------------------------------

/// Cycloidal arc from `(-r, 0)` to `(0, -r)`.
///
/// Parametric curve over `t in [0, pi]`:
///   `d1(t) = -r + r * (t - sin(t)) / pi`
///   `d2(t) = -r * (1 - cos(t)) / 2`
/// At t=0: (-r, 0).  At t=pi: (0, -r).
fn blend_intersection_cycloidal(d1: f64, d2: f64, r: f64) -> f64 {
    let sharp = d1.max(d2);
    if d1 > r || d2 > r {
        return sharp;
    }
    if d1 < -r && d2 < -r {
        return sharp;
    }

    let q = [d1, d2];
    let pi = std::f64::consts::PI;

    // Reparametrize: u in [0, 1], physical parameter theta = u * pi.
    let curve = |u: f64| -> [f64; 2] {
        let theta = u * pi;
        [
            -r + r * (theta - theta.sin()) / pi,
            -r * (1.0 - theta.cos()) / 2.0,
        ]
    };

    let mut u = best_initial_t(q, &curve, 16);
    for _ in 0..NEWTON_ITERS {
        let theta = u * pi;
        let sin_t = theta.sin();
        let cos_t = theta.cos();

        let px = -r + r * (theta - sin_t) / pi;
        let py = -r * (1.0 - cos_t) / 2.0;

        // d/du derivatives (chain rule: dtheta/du = pi)
        let dpx = r * (1.0 - cos_t); // r * (1 - cos(theta)) * pi / pi
        let dpy = -r * sin_t * pi / 2.0; // -r * sin(theta) * pi / 2

        let d2px = r * sin_t * pi; // r * sin(theta) * pi
        let d2py = -r * cos_t * pi * pi / 2.0;

        let diff = [px - q[0], py - q[1]];
        let num = diff[0] * dpx + diff[1] * dpy;
        let den = dpx * dpx + dpy * dpy + diff[0] * d2px + diff[1] * d2py;
        if den.abs() < 1e-15 {
            break;
        }
        u = (u - num / den).clamp(0.0, 1.0);
    }

    let theta = u * pi;
    let nearest = [
        -r + r * (theta - theta.sin()) / pi,
        -r * (1.0 - theta.cos()) / 2.0,
    ];
    let tangent = [r * (1.0 - theta.cos()), -r * theta.sin() * pi / 2.0];
    let diff = [q[0] - nearest[0], q[1] - nearest[1]];
    let dist = (diff[0] * diff[0] + diff[1] * diff[1]).sqrt();

    // At cusps (u ≈ 0 or u ≈ 1), the tangent is zero.  Use a slightly offset
    // parameter to obtain a non-degenerate tangent for sign determination.
    let tangent_len_sq = tangent[0] * tangent[0] + tangent[1] * tangent[1];
    let sign = if tangent_len_sq > 1e-20 {
        cross2(tangent, diff)
    } else {
        let nudge = if u < 0.5 { 1e-4 } else { 1.0 - 1e-4 };
        let nt = nudge * pi;
        let fallback = [r * (1.0 - nt.cos()), -r * nt.sin() * pi / 2.0];
        cross2(fallback, diff)
    };
    let signed_dist = if sign >= 0.0 { dist } else { -dist };
    sharp.max(signed_dist)
}

// ---------------------------------------------------------------------------
// Hyperbolic blend (Newton iteration)
// ---------------------------------------------------------------------------

/// Hyperbolic-like blend from `(-r, 0)` to `(0, -r)` using a superellipse
/// with exponent `p < 1`.
///
/// The exponent `p = (asymptote / r).clamp(0.1, 0.99)` controls how much the
/// curve "bows in" toward the sharp edge (lower p = more concave).
///
/// Parametric curve:
///   `d1(t) = -r * (1 - t^p)^(1/p)`,  `d2(t) = -r * t`
/// At t=0: d1=-r, d2=0.  At t=1: d1=0, d2=-r.
fn blend_intersection_hyperbolic(d1: f64, d2: f64, r: f64, asymptote: f64) -> f64 {
    let sharp = d1.max(d2);
    if d1 > r || d2 > r {
        return sharp;
    }
    if d1 < -r && d2 < -r {
        return sharp;
    }

    let q = [d1, d2];
    let p = (asymptote / r).clamp(0.1, 0.99);
    let inv_p = 1.0 / p;

    // Curve helper that clamps t away from singularities.
    let curve_pt = |t: f64| -> [f64; 2] {
        let tc = t.clamp(1e-8, 1.0 - 1e-8);
        let tp = tc.powf(p);
        let base = (1.0 - tp).max(1e-15).powf(inv_p);
        [-r * base, -r * tc]
    };

    let mut t = best_initial_t(q, &|t| curve_pt(t), 16);

    for _ in 0..NEWTON_ITERS {
        let tc = t.clamp(1e-8, 1.0 - 1e-8);
        let tp = tc.powf(p);
        let one_minus_tp = (1.0 - tp).max(1e-15);
        let base = one_minus_tp.powf(inv_p);

        let px = -r * base;
        let py = -r * tc;

        // d1'(t) = -r * (1/p) * (1-t^p)^(1/p-1) * (-p * t^(p-1))
        //        =  r * t^(p-1) * (1-t^p)^(1/p-1)
        let dpx = r * tc.powf(p - 1.0) * one_minus_tp.powf(inv_p - 1.0);
        let dpy = -r;

        // Numerical second derivative for d1 (avoids complex analytic form).
        let dt = 1e-6;
        let t_plus = (tc + dt).min(1.0 - 1e-8);
        let t_minus = (tc - dt).max(1e-8);
        let px_plus = -r * (1.0 - t_plus.powf(p)).max(1e-15).powf(inv_p);
        let px_minus = -r * (1.0 - t_minus.powf(p)).max(1e-15).powf(inv_p);
        let d2px = (px_plus - 2.0 * px + px_minus) / (dt * dt);
        let d2py = 0.0;

        let diff = [px - q[0], py - q[1]];
        let num = diff[0] * dpx + diff[1] * dpy;
        let den = dpx * dpx + dpy * dpy + diff[0] * d2px + diff[1] * d2py;
        if den.abs() < 1e-15 {
            break;
        }
        t = (t - num / den).clamp(0.01, 0.99);
    }

    let tc = t.clamp(1e-8, 1.0 - 1e-8);
    let tp = tc.powf(p);
    let one_minus_tp = (1.0 - tp).max(1e-15);
    let base = one_minus_tp.powf(inv_p);
    let nearest = [-r * base, -r * tc];
    let tangent = [r * tc.powf(p - 1.0) * one_minus_tp.powf(inv_p - 1.0), -r];
    let diff = [q[0] - nearest[0], q[1] - nearest[1]];
    let dist = (diff[0] * diff[0] + diff[1] * diff[1]).sqrt();

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
