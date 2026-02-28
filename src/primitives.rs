use nalgebra::Vector3;

/// Signed distance to a sphere.
/// Negative inside, zero on surface, positive outside.
pub fn sdf_sphere(point: Vector3<f64>, center: Vector3<f64>, radius: f64) -> f64 {
    (point - center).norm() - radius
}

/// Signed distance to an axis-aligned box.
/// `half_extents` is the half-size in each dimension.
pub fn sdf_box(point: Vector3<f64>, center: Vector3<f64>, half_extents: Vector3<f64>) -> f64 {
    let d = (point - center).abs() - half_extents;
    let outside = Vector3::new(d.x.max(0.0), d.y.max(0.0), d.z.max(0.0)).norm();
    let inside = d.x.max(d.y).max(d.z).min(0.0);
    outside + inside
}

/// Signed distance to a capped cylinder.
/// `base` is the center of the bottom cap, `axis` is the unit direction,
/// `radius` is the cross-section radius, `height` is the length along axis.
/// Exact signed distance to a capped cone (or truncated cone).
///
/// Connects two circles: one centered at `a` with radius `ra`, and another
/// at `b` with radius `rb`. For a sharp cone, set `rb = 0`. Handles the
/// base caps, tip, and conical wall correctly with no approximation.
pub fn sdf_capped_cone(
    point: Vector3<f64>,
    a: Vector3<f64>,
    b: Vector3<f64>,
    ra: f64,
    rb: f64,
) -> f64 {
    sdf_capped_cone_with_normal(point, a, b, ra, rb).0
}

/// Exact signed distance to a capped cone, returning both distance and the
/// analytical surface normal (unit gradient of the SDF).
///
/// Uses a closest-point decomposition: finds the nearest point on the cone
/// surface (wall, cap faces, rims, or tip) and derives the gradient from
/// the direction between the query point and that closest point.
///
/// For exterior points: gradient = (P - closest).normalize()
/// For interior points: gradient = (closest - P).normalize()
/// Both cases yield the outward surface normal at the closest point.
pub fn sdf_capped_cone_with_normal(
    point: Vector3<f64>,
    a: Vector3<f64>,
    b: Vector3<f64>,
    ra: f64,
    rb: f64,
) -> (f64, Vector3<f64>) {
    let ab = b - a;
    let ab_len_sq = ab.dot(&ab);
    let ab_len = ab_len_sq.sqrt();

    // Degenerate: a == b → sphere of max radius
    if ab_len < 1e-15 {
        let r = ra.max(rb);
        let d = (point - a).norm() - r;
        let dir = point - a;
        let n = if dir.norm() > 1e-15 {
            dir.normalize()
        } else {
            Vector3::new(0.0, 1.0, 0.0)
        };
        return (d, n);
    }

    let axis = ab / ab_len; // unit axis direction A→B

    // Project point onto cone axis
    let pa = point - a;
    let t_raw = pa.dot(&axis); // signed distance along axis from A
    let t = t_raw / ab_len; // normalized: 0 at A, 1 at B

    // Radial vector (perpendicular to axis)
    let radial_vec = pa - axis * t_raw;
    let r = radial_vec.norm(); // radial distance from axis

    // Cone radius at axial parameter t (clamped to [0,1])
    let r_at_t = ra + t.clamp(0.0, 1.0) * (rb - ra);

    // Inside test: point is between caps AND inside the cone radius
    let inside = t >= 0.0 && t <= 1.0 && r < r_at_t;

    // Radial unit vector (direction from axis toward point)
    let radial_unit = if r > 1e-15 {
        radial_vec / r
    } else {
        // On axis: pick an arbitrary perpendicular direction
        let arb = if axis.x.abs() < 0.9 {
            Vector3::new(1.0, 0.0, 0.0)
        } else {
            Vector3::new(0.0, 1.0, 0.0)
        };
        (arb - axis * axis.dot(&arb)).normalize()
    };

    // --- Find the closest point on the cone surface ---
    // We track (closest_point, geometric_normal) for each candidate feature.
    // geometric_normal is the outward surface normal at the closest point,
    // used as fallback when the query point is very close to the surface.

    let mut best_dist_sq = f64::INFINITY;
    let mut best_closest = a;
    let mut best_geo_normal = Vector3::new(0.0, 1.0, 0.0);

    // Wall outward normal (constant along the entire conical wall)
    let slant_len = ((ra - rb) * (ra - rb) + ab_len_sq).sqrt();
    let wall_n_radial = ab_len / slant_len;
    let wall_n_axial = (ra - rb) / slant_len;
    let wall_geo_normal = radial_unit * wall_n_radial + axis * wall_n_axial;

    // Project onto the conical wall line in (axial, radial) space:
    // Wall line: from (0, ra) to (ab_len, rb) parametrized by s ∈ [0, 1]
    let dr = rb - ra;
    let s_unclamped = (t_raw * ab_len + (r - ra) * dr) / (ab_len_sq + dr * dr);
    let s_wall = s_unclamped.clamp(0.0, 1.0);
    let wall_r = ra + s_wall * dr;
    let wall_axial = s_wall * ab_len;
    let wall_closest = a + axis * wall_axial + radial_unit * wall_r;

    // Wall candidate — only use wall geometric normal when the projection
    // didn't clamp (i.e., the closest point is on the wall interior, not a rim).
    {
        let diff = point - wall_closest;
        let d_sq = diff.dot(&diff);
        if d_sq < best_dist_sq {
            best_dist_sq = d_sq;
            best_closest = wall_closest;
            // If s_wall was clamped, the closest point is actually on a rim/tip,
            // so the geometric wall normal is wrong. We'll use direction-based
            // normal for those cases (handled by the final normal computation below).
            best_geo_normal = wall_geo_normal;
        }
    }

    // Cap A candidates
    if ra > 1e-15 {
        // Cap A face (flat disk at a, radius ra)
        if r <= ra {
            let cap_closest = a + radial_vec; // project onto cap plane
            let diff = point - cap_closest;
            let d_sq = diff.dot(&diff);
            if d_sq < best_dist_sq {
                best_dist_sq = d_sq;
                best_closest = cap_closest;
                best_geo_normal = -axis;
            }
        }
        // Cap A rim (circle at a, radius ra)
        {
            let rim_point = a + radial_unit * ra;
            let diff = point - rim_point;
            let d_sq = diff.dot(&diff);
            if d_sq < best_dist_sq {
                best_dist_sq = d_sq;
                best_closest = rim_point;
                best_geo_normal = -axis; // fallback; direction-based will override
            }
        }
    } else {
        // ra == 0: point tip at a
        let diff = point - a;
        let d_sq = diff.dot(&diff);
        if d_sq < best_dist_sq {
            best_dist_sq = d_sq;
            best_closest = a;
            best_geo_normal = -axis;
        }
    }

    // Cap B candidates
    if rb > 1e-15 {
        // Cap B face (flat disk at b, radius rb)
        if r <= rb {
            let cap_closest = b + radial_vec; // project onto cap plane
            let diff = point - cap_closest;
            let d_sq = diff.dot(&diff);
            if d_sq < best_dist_sq {
                best_dist_sq = d_sq;
                best_closest = cap_closest;
                best_geo_normal = axis;
            }
        }
        // Cap B rim
        {
            let rim_point = b + radial_unit * rb;
            let diff = point - rim_point;
            let d_sq = diff.dot(&diff);
            if d_sq < best_dist_sq {
                best_dist_sq = d_sq;
                best_closest = rim_point;
                best_geo_normal = axis;
            }
        }
    } else {
        // rb == 0: point tip at b
        let diff = point - b;
        let d_sq = diff.dot(&diff);
        if d_sq < best_dist_sq {
            best_dist_sq = d_sq;
            best_closest = b;
            best_geo_normal = axis;
        }
    }

    let dist = best_dist_sq.sqrt();
    let signed_dist = if inside { -dist } else { dist };

    // Compute the gradient (outward surface normal).
    // For exterior: gradient = (P - closest).normalize()
    // For interior: gradient = (closest - P).normalize()
    // When the point is very close to the surface, this direction becomes
    // numerically unstable, so we fall back to the geometric normal.
    let diff = point - best_closest;
    let diff_len = diff.norm();

    let normal = if diff_len > 1e-10 {
        if inside {
            -diff / diff_len // points outward (from interior toward surface)
        } else {
            diff / diff_len // points outward (from surface toward exterior)
        }
    } else {
        // On or very near the surface: use geometric normal
        best_geo_normal.normalize()
    };

    (signed_dist, normal)
}

pub fn sdf_cylinder(
    point: Vector3<f64>,
    base: Vector3<f64>,
    axis: Vector3<f64>,
    radius: f64,
    height: f64,
) -> f64 {
    let p = point - base;
    let along = p.dot(&axis);
    let perp = (p - axis * along).norm();
    let d_radial = perp - radius;
    let d_axial = (along - height / 2.0).abs() - height / 2.0;
    if d_radial > 0.0 && d_axial > 0.0 {
        (d_radial * d_radial + d_axial * d_axial).sqrt()
    } else {
        d_radial.max(d_axial)
    }
}

/// Signed distance to a torus centered at `center`, lying in the XZ plane.
/// `major_radius` is the distance from center to the tube center.
/// `minor_radius` is the tube radius.
/// Exact signed distance (Quilez).
pub fn sdf_torus(
    point: Vector3<f64>,
    center: Vector3<f64>,
    major_radius: f64,
    minor_radius: f64,
) -> f64 {
    let p = point - center;
    let q_x = (p.x * p.x + p.z * p.z).sqrt() - major_radius;
    let q_y = p.y;
    (q_x * q_x + q_y * q_y).sqrt() - minor_radius
}

/// Signed distance to a box with rounded edges.
/// `half_extents` is the sharp box half-size, `radius` is the rounding radius.
/// The total size is `half_extents + radius` in each dimension.
/// Exact signed distance (Quilez).
pub fn sdf_rounded_box(
    point: Vector3<f64>,
    center: Vector3<f64>,
    half_extents: Vector3<f64>,
    radius: f64,
) -> f64 {
    let d = (point - center).abs() - half_extents;
    let outside = Vector3::new(d.x.max(0.0), d.y.max(0.0), d.z.max(0.0)).norm();
    let inside = d.x.max(d.y).max(d.z).min(0.0);
    outside + inside - radius
}

/// Signed distance to a capsule (sphere-swept line segment).
/// `a` and `b` are the endpoints, `radius` is the swept sphere radius.
/// Exact signed distance.
pub fn sdf_capsule(point: Vector3<f64>, a: Vector3<f64>, b: Vector3<f64>, radius: f64) -> f64 {
    let pa = point - a;
    let ba = b - a;
    let h = (pa.dot(&ba) / ba.dot(&ba)).clamp(0.0, 1.0);
    (pa - ba * h).norm() - radius
}

/// Signed distance to an ellipsoid centered at `center` with semi-axis lengths `radii`.
/// This is an approximation — the true ellipsoid SDF has no closed-form solution.
/// Uses the bound-corrected approximation from Quilez.
pub fn sdf_ellipsoid(point: Vector3<f64>, center: Vector3<f64>, radii: Vector3<f64>) -> f64 {
    let p = point - center;
    // Normalize by radii
    let k0 = Vector3::new(p.x / radii.x, p.y / radii.y, p.z / radii.z).norm();
    let k1 = Vector3::new(
        p.x / (radii.x * radii.x),
        p.y / (radii.y * radii.y),
        p.z / (radii.z * radii.z),
    )
    .norm();
    if k1 > 1e-10 {
        k0 * (k0 - 1.0) / k1
    } else {
        -radii.x.min(radii.y).min(radii.z) // degenerate: inside at center
    }
}

/// Signed distance to a rounded cylinder (cylinder with rounded edges at top/bottom).
/// `radius` is the cylinder radius, `round_radius` is the edge rounding,
/// `height` is the half-height. Axis-aligned along Y.
/// Exact signed distance (Quilez).
pub fn sdf_rounded_cylinder(
    point: Vector3<f64>,
    center: Vector3<f64>,
    radius: f64,
    round_radius: f64,
    half_height: f64,
) -> f64 {
    let p = point - center;
    let d_radial = (p.x * p.x + p.z * p.z).sqrt() - radius + round_radius;
    let d_axial = p.y.abs() - half_height + round_radius;
    let outside = Vector3::new(d_radial.max(0.0), d_axial.max(0.0), 0.0).norm();
    let inside = d_radial.max(d_axial).min(0.0);
    outside + inside - round_radius
}
