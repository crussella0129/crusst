use crate::path::Path;
use crate::shape::Sdf;
use nalgebra::Vector3;

/// Polynomial smooth minimum. Blends two SDF values smoothly within
/// radius `k`, making the combined gradient C1 continuous. At most shifts
/// the surface by k/4 at segment boundaries.
#[inline]
fn smin(a: f64, b: f64, k: f64) -> f64 {
    let h = (k - (a - b).abs()).max(0.0) / k;
    a.min(b) - h * h * k * 0.25
}

/// Order 0: S = C(s). Shape IS the cross-section. No transport.
pub fn order0<S: Sdf + 'static>(section: S) -> TransportShape {
    TransportShape {
        eval: Box::new(move |point| section.evaluate(point)),
    }
}

/// Compute the SDF for a sphere-swept line segment (capsule) with varying radius.
///
/// For each segment between two path samples, finds the closest point on the
/// segment via projection, interpolates the radius, and returns the distance.
/// Adjacent capsules share endpoints, guaranteeing a smooth, continuous surface
/// with no gaps or bead artifacts.
#[inline]
fn capsule_segment_sdf(
    point: Vector3<f64>,
    p0: Vector3<f64>,
    p1: Vector3<f64>,
    r0: f64,
    r1: f64,
) -> f64 {
    let seg = p1 - p0;
    let seg_len_sq = seg.dot(&seg);

    if seg_len_sq < 1e-12 {
        // Degenerate segment: treat as a sphere
        let dist = (point - p0).norm();
        return dist - r0;
    }

    let offset = point - p0;
    let alpha = (offset.dot(&seg) / seg_len_sq).clamp(0.0, 1.0);

    let closest = p0 + seg * alpha;
    let radius = r0 + (r1 - r0) * alpha;

    let dist = (point - closest).norm();
    dist - radius
}

/// Order 1: Rigid transport with scaling.
/// Circle of `section_radius` swept along path with `scale_fn(t)` applied uniformly.
/// Produces cones (scale linear), cylinders (scale constant), etc.
///
/// Uses capsule sweep with polynomial smooth-min: each pair of adjacent path
/// samples forms a capsule (sphere-swept line segment) with linearly
/// interpolated radius. The smooth union of all capsules gives a C1-smooth
/// surface with no staircase normal artifacts at segment boundaries.
pub fn order1<P: Path + 'static>(
    path: P,
    section_radius: f64,
    scale_fn: impl Fn(f64) -> f64 + Send + Sync + 'static,
    samples: usize,
) -> TransportShape {
    // Smooth-min radius: small enough to preserve surface accuracy (max shift
    // = k/4 ≈ 0.075 units) but large enough that central-difference normals
    // sample the blend zone, producing C1 gradients.
    let k = 0.3;
    TransportShape {
        eval: Box::new(move |point| {
            let mut min_dist = f64::MAX;

            for i in 0..samples {
                let t0 = i as f64 / samples as f64;
                let t1 = (i + 1) as f64 / samples as f64;

                let p0 = path.point(t0);
                let p1 = path.point(t1);

                let s0 = scale_fn(t0);
                let s1 = scale_fn(t1);

                let r0 = (section_radius * s0).max(0.0);
                let r1 = (section_radius * s1).max(0.0);

                let d = capsule_segment_sdf(point, p0, p1, r0, r1);
                min_dist = smin(min_dist, d, k);
            }

            min_dist
        }),
    }
}

/// Order 2: Frame transport.
/// Circle of `section_radius` swept along path with Frenet frame following curvature.
/// The section is rigid — only the frame rotates.
pub fn order2<P: Path + 'static>(path: P, section_radius: f64, samples: usize) -> TransportShape {
    order1(path, section_radius, |_| 1.0, samples)
}

/// Order 3: Section evolution (morph transport).
/// The cross-section can change scale and twist along the path.
/// `scale_fn(t)` controls taper, `twist_fn(t)` controls rotation about the path tangent.
/// The horn is the eigenform: monotonic taper + smooth twist along a spiral.
///
/// Note: For circular cross-sections, twist has no visible effect (rotating a
/// circle about its center is an identity). Twist will become visible when
/// non-circular cross-sections are supported in a future phase.
pub fn order3<P: Path + 'static>(
    path: P,
    section_radius: f64,
    scale_fn: impl Fn(f64) -> f64 + Send + Sync + 'static,
    _twist_fn: impl Fn(f64) -> f64 + Send + Sync + 'static,
    samples: usize,
) -> TransportShape {
    // For circular cross-sections, twist is invisible (rotation symmetry).
    // Use the same capsule sweep as order1 — twist will be implemented
    // when non-circular cross-sections are added.
    let k = 0.3;
    TransportShape {
        eval: Box::new(move |point| {
            let mut min_dist = f64::MAX;

            for i in 0..samples {
                let t0 = i as f64 / samples as f64;
                let t1 = (i + 1) as f64 / samples as f64;

                let p0 = path.point(t0);
                let p1 = path.point(t1);

                let s0 = scale_fn(t0);
                let s1 = scale_fn(t1);

                let r0 = (section_radius * s0).max(0.0);
                let r1 = (section_radius * s1).max(0.0);

                let d = capsule_segment_sdf(point, p0, p1, r0, r1);
                min_dist = smin(min_dist, d, k);
            }

            min_dist
        }),
    }
}

/// The output of a transport operation: a shape you can evaluate.
pub struct TransportShape {
    eval: Box<dyn Fn(Vector3<f64>) -> f64 + Send + Sync>,
}

impl Sdf for TransportShape {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        (self.eval)(point)
    }
}
