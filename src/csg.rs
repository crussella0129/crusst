// CSG (Constructive Solid Geometry) operations
//
// These combine SDF distance values to produce boolean combinations of shapes:
// - union:        the merged shape (min)
// - intersection: the overlapping region (max)
// - difference:   A minus B (max of A and negated B)
//
// Smooth variants blend the transition region, producing fillet-like surfaces
// controlled by the blending radius `k`.

/// Boolean union of two SDF values (logical OR).
/// Returns the minimum distance — the point is inside whichever shape is closer.
pub fn union(d1: f64, d2: f64) -> f64 {
    d1.min(d2)
}

/// Boolean intersection of two SDF values (logical AND).
/// Returns the maximum distance — the point must be inside both shapes.
pub fn intersection(d1: f64, d2: f64) -> f64 {
    d1.max(d2)
}

/// Boolean difference: shape A minus shape B.
/// Keeps points inside A that are outside B.
pub fn difference(d1: f64, d2: f64) -> f64 {
    d1.max(-d2)
}

/// Smooth (polynomial) union with blending radius `k`.
/// Produces a fillet-like blend where the two surfaces meet.
/// When `k == 0`, this degenerates to a sharp union.
pub fn smooth_union(d1: f64, d2: f64, k: f64) -> f64 {
    let h = (0.5 + 0.5 * (d2 - d1) / k).clamp(0.0, 1.0);
    d2 * (1.0 - h) + d1 * h - k * h * (1.0 - h)
}

/// Smooth (polynomial) intersection with blending radius `k`.
/// Produces a fillet-like blend in the intersection region.
pub fn smooth_intersection(d1: f64, d2: f64, k: f64) -> f64 {
    let h = (0.5 - 0.5 * (d2 - d1) / k).clamp(0.0, 1.0);
    d2 * (1.0 - h) + d1 * h + k * h * (1.0 - h)
}

/// Smooth difference: shape A minus shape B with blending radius `k`.
/// Equivalent to `smooth_intersection(d1, -d2, k)`.
pub fn smooth_difference(d1: f64, d2: f64, k: f64) -> f64 {
    smooth_intersection(d1, -d2, k)
}

/// Round union (inscribed-circle blend) with radius `r`.
/// Produces an exact circular fillet where two surfaces meet in a union.
pub fn round_union(d1: f64, d2: f64, r: f64) -> f64 {
    let u = (r - d1).max(0.0);
    let v = (r - d2).max(0.0);
    r.max(d1.min(d2)) - (u * u + v * v).sqrt()
}

/// Round intersection with radius `r`.
/// Produces an exact circular fillet in the intersection region.
pub fn round_intersection(d1: f64, d2: f64, r: f64) -> f64 {
    let u = (r + d1).max(0.0);
    let v = (r + d2).max(0.0);
    (-r).min(d1.max(d2)) + (u * u + v * v).sqrt()
}

/// Round difference with radius `r`.
pub fn round_difference(d1: f64, d2: f64, r: f64) -> f64 {
    round_intersection(d1, -d2, r)
}

/// Chamfer union with chamfer size `k`.
pub fn chamfer_union(d1: f64, d2: f64, k: f64) -> f64 {
    d1.min(d2).min((d1 + d2 - k) / std::f64::consts::SQRT_2)
}

/// Chamfer intersection with chamfer size `k`.
pub fn chamfer_intersection(d1: f64, d2: f64, k: f64) -> f64 {
    d1.max(d2).max((d1 + d2 + k) / std::f64::consts::SQRT_2)
}

/// Chamfer difference with chamfer size `k`.
pub fn chamfer_difference(d1: f64, d2: f64, k: f64) -> f64 {
    chamfer_intersection(d1, -d2, k)
}
