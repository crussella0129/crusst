//! SdfNode expression DAG — an enum-based composable SDF tree.
//!
//! Unlike the trait-based API in `shape.rs` (which uses generics + type erasure),
//! the DAG uses `Arc<SdfNode>` for children. This enables:
//! - Interval arithmetic for octree pruning
//! - Analytical gradients for QEF
//! - DAG introspection for smart STEP export

use crate::feature::{EdgeInfo, FaceInfo};
use crate::types::{BBox3, Interval};
use crate::{csg, primitives};
use nalgebra::{Rotation3, Vector2, Vector3};
use std::sync::Arc;

// ---------------------------------------------------------------------------
// SdfNode2d — 2D profile sub-enum (for Revolve / Extrude)
// ---------------------------------------------------------------------------

/// A 2D signed distance function node, used as the cross-section profile
/// for `SdfNode::Revolve` and `SdfNode::Extrude`.
pub enum SdfNode2d {
    /// 2D circle.
    Circle2d { center: Vector2<f64>, radius: f64 },
    /// 2D axis-aligned rectangle.
    Rect2d {
        center: Vector2<f64>,
        half_extents: Vector2<f64>,
    },
    /// Boolean union of two 2D shapes.
    Union2d(Arc<SdfNode2d>, Arc<SdfNode2d>),
    /// Boolean difference of two 2D shapes (A minus B).
    Difference2d(Arc<SdfNode2d>, Arc<SdfNode2d>),
    /// Opaque 2D SDF — wraps any `Sdf2d` trait object.
    Custom2d(Arc<dyn crate::shape::Sdf2d>),
}

impl SdfNode2d {
    /// Evaluate the 2D signed distance at a point.
    pub fn evaluate(&self, point: Vector2<f64>) -> f64 {
        match self {
            SdfNode2d::Circle2d { center, radius } => (point - center).norm() - radius,
            SdfNode2d::Rect2d {
                center,
                half_extents,
            } => {
                let d = (point - center).abs() - half_extents;
                let outside = Vector2::new(d.x.max(0.0), d.y.max(0.0)).norm();
                let inside = d.x.max(d.y).min(0.0);
                outside + inside
            }
            SdfNode2d::Union2d(a, b) => a.evaluate(point).min(b.evaluate(point)),
            SdfNode2d::Difference2d(a, b) => a.evaluate(point).max(-b.evaluate(point)),
            SdfNode2d::Custom2d(sdf) => sdf.evaluate(point),
        }
    }
}

// ---------------------------------------------------------------------------
// SdfNode — the main 3D DAG
// ---------------------------------------------------------------------------

/// An expression DAG node representing a signed distance function.
///
/// Each variant stores its parameters inline. Children are held behind
/// `Arc<SdfNode>` so the DAG can be cheaply cloned and shared.
pub enum SdfNode {
    // -- Primitives (10) ---------------------------------------------------
    /// Sphere at `center` with given `radius`.
    Sphere { center: Vector3<f64>, radius: f64 },

    /// Axis-aligned box at `center` with `half_extents`.
    Box3 {
        center: Vector3<f64>,
        half_extents: Vector3<f64>,
    },

    /// Capped cylinder with `base` center, unit `axis`, `radius`, and `height`.
    Cylinder {
        base: Vector3<f64>,
        axis: Vector3<f64>,
        radius: f64,
        height: f64,
    },

    /// Capped cone (truncated cone) from `a` (radius `ra`) to `b` (radius `rb`).
    CappedCone {
        a: Vector3<f64>,
        b: Vector3<f64>,
        ra: f64,
        rb: f64,
    },

    /// Torus lying in the XZ plane, centered at `center`.
    Torus {
        center: Vector3<f64>,
        major_radius: f64,
        minor_radius: f64,
    },

    /// Box with rounded edges. Total size is `half_extents + radius`.
    RoundedBox {
        center: Vector3<f64>,
        half_extents: Vector3<f64>,
        radius: f64,
    },

    /// Capsule (sphere-swept segment) from `a` to `b` with `radius`.
    Capsule {
        a: Vector3<f64>,
        b: Vector3<f64>,
        radius: f64,
    },

    /// Ellipsoid with semi-axis lengths `radii`. Approximate SDF.
    Ellipsoid {
        center: Vector3<f64>,
        radii: Vector3<f64>,
    },

    /// Rounded cylinder aligned along Y axis.
    RoundedCylinder {
        center: Vector3<f64>,
        radius: f64,
        round_radius: f64,
        half_height: f64,
    },

    /// Half-space: solid behind the plane `normal . p + d = 0`.
    HalfSpace { normal: Vector3<f64>, d: f64 },

    // -- CSG Operations (6) ------------------------------------------------
    /// Boolean union (logical OR) — min of children.
    Union(Arc<SdfNode>, Arc<SdfNode>),

    /// Boolean intersection (logical AND) — max of children.
    Intersection(Arc<SdfNode>, Arc<SdfNode>),

    /// Boolean difference: A minus B.
    Difference(Arc<SdfNode>, Arc<SdfNode>),

    /// Smooth union with blending radius `k`.
    SmoothUnion(Arc<SdfNode>, Arc<SdfNode>, f64),

    /// Smooth intersection with blending radius `k`.
    SmoothIntersection(Arc<SdfNode>, Arc<SdfNode>, f64),

    /// Smooth difference with blending radius `k`.
    SmoothDifference(Arc<SdfNode>, Arc<SdfNode>, f64),

    // -- Transforms (6) ----------------------------------------------------
    /// Translate by an offset vector.
    Translate(Arc<SdfNode>, Vector3<f64>),

    /// Rotate by a rotation matrix.
    Rotate(Arc<SdfNode>, Rotation3<f64>),

    /// Uniform scale by a factor.
    Scale(Arc<SdfNode>, f64),

    /// Mirror across a plane through the origin with given normal.
    Mirror(Arc<SdfNode>, Vector3<f64>),

    /// Shell (onion): hollow with wall `thickness`.
    Shell(Arc<SdfNode>, f64),

    /// Round (offset): expand/contract surface by `radius`.
    Round(Arc<SdfNode>, f64),

    // -- 2D → 3D (2) ------------------------------------------------------
    /// Revolve a 2D profile around the Y axis.
    Revolve(Arc<SdfNode2d>),

    /// Extrude a 2D profile along the Z axis by `half_height`.
    Extrude(Arc<SdfNode2d>, f64),

    // -- Opaque ------------------------------------------------------------
    /// Wraps any existing `Sdf` trait object for interop.
    Custom(Arc<dyn crate::shape::Sdf>),
}

// Compile-time assertions: SdfNode and SdfNode2d must be Send + Sync
// so the DAG can be shared across threads during parallel octree traversal.
const _: () = {
    #[allow(dead_code)]
    fn assert_send_sync<T: Send + Sync>() {}
    #[allow(dead_code)]
    fn check() {
        assert_send_sync::<SdfNode>();
        assert_send_sync::<SdfNode2d>();
    }
};

impl SdfNode {
    /// Evaluate the signed distance at a 3D point.
    pub fn evaluate(&self, point: Vector3<f64>) -> f64 {
        match self {
            // -- Primitives ------------------------------------------------
            SdfNode::Sphere { center, radius } => primitives::sdf_sphere(point, *center, *radius),
            SdfNode::Box3 {
                center,
                half_extents,
            } => primitives::sdf_box(point, *center, *half_extents),
            SdfNode::Cylinder {
                base,
                axis,
                radius,
                height,
            } => primitives::sdf_cylinder(point, *base, *axis, *radius, *height),
            SdfNode::CappedCone { a, b, ra, rb } => {
                primitives::sdf_capped_cone(point, *a, *b, *ra, *rb)
            }
            SdfNode::Torus {
                center,
                major_radius,
                minor_radius,
            } => primitives::sdf_torus(point, *center, *major_radius, *minor_radius),
            SdfNode::RoundedBox {
                center,
                half_extents,
                radius,
            } => primitives::sdf_rounded_box(point, *center, *half_extents, *radius),
            SdfNode::Capsule { a, b, radius } => primitives::sdf_capsule(point, *a, *b, *radius),
            SdfNode::Ellipsoid { center, radii } => {
                primitives::sdf_ellipsoid(point, *center, *radii)
            }
            SdfNode::RoundedCylinder {
                center,
                radius,
                round_radius,
                half_height,
            } => primitives::sdf_rounded_cylinder(
                point,
                *center,
                *radius,
                *round_radius,
                *half_height,
            ),
            SdfNode::HalfSpace { normal, d } => normal.dot(&point) + d,

            // -- CSG -------------------------------------------------------
            SdfNode::Union(a, b) => csg::union(a.evaluate(point), b.evaluate(point)),
            SdfNode::Intersection(a, b) => csg::intersection(a.evaluate(point), b.evaluate(point)),
            SdfNode::Difference(a, b) => csg::difference(a.evaluate(point), b.evaluate(point)),
            SdfNode::SmoothUnion(a, b, k) => {
                csg::smooth_union(a.evaluate(point), b.evaluate(point), *k)
            }
            SdfNode::SmoothIntersection(a, b, k) => {
                csg::smooth_intersection(a.evaluate(point), b.evaluate(point), *k)
            }
            SdfNode::SmoothDifference(a, b, k) => {
                csg::smooth_difference(a.evaluate(point), b.evaluate(point), *k)
            }

            // -- Transforms ------------------------------------------------
            SdfNode::Translate(inner, offset) => inner.evaluate(point - offset),
            SdfNode::Rotate(inner, rotation) => inner.evaluate(rotation.inverse() * point),
            SdfNode::Scale(inner, factor) => {
                let inv = 1.0 / factor;
                inner.evaluate(point * inv) * factor
            }
            SdfNode::Mirror(inner, normal) => {
                // `normal` is expected to be unit-length; the builder API
                // normalizes at construction time so we skip it on the hot path.
                let d = point.dot(normal);
                let reflected = point - normal * (2.0 * d.min(0.0));
                inner.evaluate(reflected)
            }
            SdfNode::Shell(inner, thickness) => inner.evaluate(point).abs() - thickness,
            SdfNode::Round(inner, radius) => inner.evaluate(point) - radius,

            // -- 2D → 3D --------------------------------------------------
            SdfNode::Revolve(profile) => {
                let r = (point.x * point.x + point.z * point.z).sqrt();
                profile.evaluate(Vector2::new(r, point.y))
            }
            SdfNode::Extrude(profile, half_height) => {
                let d_2d = profile.evaluate(Vector2::new(point.x, point.y));
                let d_z = point.z.abs() - half_height;
                if d_2d > 0.0 && d_z > 0.0 {
                    (d_2d * d_2d + d_z * d_z).sqrt()
                } else {
                    d_2d.max(d_z)
                }
            }

            // -- Opaque ----------------------------------------------------
            SdfNode::Custom(sdf) => sdf.evaluate(point),
        }
    }

    /// Evaluate a conservative interval bound of the SDF over an axis-aligned
    /// bounding box. The returned interval is guaranteed to contain every
    /// value that `evaluate(p)` can take for any `p` inside `bbox`.
    ///
    /// The interval may be wider than the true range (less pruning), but is
    /// **never narrower** (would miss surface crossings).
    pub fn interval_evaluate(&self, bbox: &BBox3) -> Interval {
        match self {
            // -- Primitives ------------------------------------------------
            SdfNode::Sphere { center, radius } => {
                // Per-axis interval of (p - center)
                let x_iv = Interval::new(bbox.min.x - center.x, bbox.max.x - center.x);
                let y_iv = Interval::new(bbox.min.y - center.y, bbox.max.y - center.y);
                let z_iv = Interval::new(bbox.min.z - center.z, bbox.max.z - center.z);
                // |p - center|^2 = x^2 + y^2 + z^2 (interval of each squared, summed)
                let dist_sq = x_iv.mul(x_iv).add(y_iv.mul(y_iv)).add(z_iv.mul(z_iv));
                // |p - center| = sqrt of that, then subtract radius
                dist_sq.sqrt().scalar_sub(*radius)
            }

            SdfNode::Box3 {
                center,
                half_extents,
            } => {
                // sdf_box: d_i = |p_i - center_i| - half_i
                // outside = norm(max(d, 0)), inside = max(d_x, d_y, d_z).min(0)
                // result = outside + inside
                //
                // Conservative: compute interval of each d_i, then bound.
                let x_iv = Interval::new(bbox.min.x - center.x, bbox.max.x - center.x)
                    .abs()
                    .scalar_sub(half_extents.x);
                let y_iv = Interval::new(bbox.min.y - center.y, bbox.max.y - center.y)
                    .abs()
                    .scalar_sub(half_extents.y);
                let z_iv = Interval::new(bbox.min.z - center.z, bbox.max.z - center.z)
                    .abs()
                    .scalar_sub(half_extents.z);

                // Outside distance: norm of max(d, 0) per axis
                let x_clamped = Interval::new(x_iv.lo.max(0.0), x_iv.hi.max(0.0));
                let y_clamped = Interval::new(y_iv.lo.max(0.0), y_iv.hi.max(0.0));
                let z_clamped = Interval::new(z_iv.lo.max(0.0), z_iv.hi.max(0.0));
                let outside_sq = x_clamped
                    .mul(x_clamped)
                    .add(y_clamped.mul(y_clamped))
                    .add(z_clamped.mul(z_clamped));
                let outside_iv = outside_sq.sqrt();

                // Inside distance: min(max(d_x, d_y, d_z), 0)
                let max_d = x_iv.max(y_iv).max(z_iv);
                let inside_iv = Interval::new(max_d.lo.min(0.0), max_d.hi.min(0.0));

                outside_iv.add(inside_iv)
            }

            SdfNode::HalfSpace { normal, d } => {
                // normal.dot(p) + d
                let x_iv = Interval::new(bbox.min.x, bbox.max.x).scalar_mul(normal.x);
                let y_iv = Interval::new(bbox.min.y, bbox.max.y).scalar_mul(normal.y);
                let z_iv = Interval::new(bbox.min.z, bbox.max.z).scalar_mul(normal.z);
                x_iv.add(y_iv).add(z_iv).scalar_add(*d)
            }

            // For complex primitives, use a conservative bounding approach:
            // evaluate all 8 corners and expand slightly for safety.
            SdfNode::Cylinder { .. }
            | SdfNode::CappedCone { .. }
            | SdfNode::Torus { .. }
            | SdfNode::RoundedBox { .. }
            | SdfNode::Capsule { .. }
            | SdfNode::Ellipsoid { .. }
            | SdfNode::RoundedCylinder { .. } => self.interval_from_corners(bbox),

            // -- CSG -------------------------------------------------------
            SdfNode::Union(a, b) => {
                let a_iv = a.interval_evaluate(bbox);
                let b_iv = b.interval_evaluate(bbox);
                a_iv.min(b_iv)
            }
            SdfNode::Intersection(a, b) => {
                let a_iv = a.interval_evaluate(bbox);
                let b_iv = b.interval_evaluate(bbox);
                a_iv.max(b_iv)
            }
            SdfNode::Difference(a, b) => {
                let a_iv = a.interval_evaluate(bbox);
                let b_iv = b.interval_evaluate(bbox);
                a_iv.max(b_iv.neg())
            }

            // Smooth CSG: start from sharp bounds, then widen by the maximum
            // smooth deviation (k/4) to stay conservative.
            SdfNode::SmoothUnion(a, b, k) => {
                let a_iv = a.interval_evaluate(bbox);
                let b_iv = b.interval_evaluate(bbox);
                // Smooth union <= sharp union, so lo can be lower by up to k/4
                let sharp = a_iv.min(b_iv);
                Interval::new(sharp.lo - k * 0.25, sharp.hi)
            }
            SdfNode::SmoothIntersection(a, b, k) => {
                let a_iv = a.interval_evaluate(bbox);
                let b_iv = b.interval_evaluate(bbox);
                // Smooth intersection >= sharp intersection, so hi can be higher by up to k/4
                let sharp = a_iv.max(b_iv);
                Interval::new(sharp.lo, sharp.hi + k * 0.25)
            }
            SdfNode::SmoothDifference(a, b, k) => {
                let a_iv = a.interval_evaluate(bbox);
                let b_iv = b.interval_evaluate(bbox);
                // Smooth difference: similar to smooth intersection of a and -b
                let sharp = a_iv.max(b_iv.neg());
                Interval::new(sharp.lo, sharp.hi + k * 0.25)
            }

            // -- Transforms ------------------------------------------------
            SdfNode::Translate(inner, offset) => {
                // Shift bbox by -offset, evaluate inner
                let shifted = BBox3::new(bbox.min - offset, bbox.max - offset);
                inner.interval_evaluate(&shifted)
            }
            SdfNode::Rotate(inner, rotation) => {
                // Conservative: compute AABB of the inverse-rotated bbox corners
                let inv = rotation.inverse();
                let corners = bbox.corners();
                let mut new_min = Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
                let mut new_max =
                    Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
                for c in &corners {
                    let rc = inv * c;
                    new_min.x = new_min.x.min(rc.x);
                    new_min.y = new_min.y.min(rc.y);
                    new_min.z = new_min.z.min(rc.z);
                    new_max.x = new_max.x.max(rc.x);
                    new_max.y = new_max.y.max(rc.y);
                    new_max.z = new_max.z.max(rc.z);
                }
                inner.interval_evaluate(&BBox3::new(new_min, new_max))
            }
            SdfNode::Scale(inner, factor) => {
                debug_assert!(*factor > 0.0, "Scale factor must be positive");
                // Scale bbox by 1/factor, evaluate inner, multiply result by factor
                let inv = 1.0 / factor;
                let scaled = BBox3::new(bbox.min * inv, bbox.max * inv);
                inner.interval_evaluate(&scaled).scalar_mul(*factor)
            }
            SdfNode::Mirror(inner, normal) => {
                // Mirror is its own inverse. Conservative: evaluate on the
                // AABB of all 8 reflected corners.
                let corners = bbox.corners();
                let mut new_min = Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
                let mut new_max =
                    Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
                for c in &corners {
                    let d = c.dot(normal);
                    let reflected = c - normal * (2.0 * d.min(0.0));
                    new_min.x = new_min.x.min(reflected.x);
                    new_min.y = new_min.y.min(reflected.y);
                    new_min.z = new_min.z.min(reflected.z);
                    new_max.x = new_max.x.max(reflected.x);
                    new_max.y = new_max.y.max(reflected.y);
                    new_max.z = new_max.z.max(reflected.z);
                }
                inner.interval_evaluate(&BBox3::new(new_min, new_max))
            }
            SdfNode::Shell(inner, thickness) => {
                // |inner(p)| - thickness
                inner.interval_evaluate(bbox).abs().scalar_sub(*thickness)
            }
            SdfNode::Round(inner, radius) => {
                // inner(p) - radius
                inner.interval_evaluate(bbox).scalar_sub(*radius)
            }

            // -- 2D -> 3D --------------------------------------------------
            SdfNode::Revolve(_) | SdfNode::Extrude(_, _) => Interval::entire(),

            // -- Opaque ----------------------------------------------------
            SdfNode::Custom(_) => Interval::entire(),
        }
    }

    /// Compute the analytical gradient (surface normal direction) at a point.
    ///
    /// The gradient of an SDF gives the direction of steepest increase.
    /// For a point on the surface, this is the outward normal. The returned
    /// vector is normalized to unit length for use in QEF solvers.
    pub fn gradient(&self, point: Vector3<f64>) -> Vector3<f64> {
        match self {
            // -- Primitives ------------------------------------------------
            SdfNode::Sphere { center, .. } => {
                let d = point - center;
                let len = d.norm();
                if len > 1e-10 {
                    d / len
                } else {
                    Vector3::new(0.0, 1.0, 0.0)
                }
            }

            SdfNode::HalfSpace { normal, .. } => *normal,

            // Complex primitives: central differences fallback
            SdfNode::Box3 { .. }
            | SdfNode::Cylinder { .. }
            | SdfNode::CappedCone { .. }
            | SdfNode::Torus { .. }
            | SdfNode::RoundedBox { .. }
            | SdfNode::Capsule { .. }
            | SdfNode::Ellipsoid { .. }
            | SdfNode::RoundedCylinder { .. } => central_diff_gradient(self, point),

            // -- CSG -------------------------------------------------------
            SdfNode::Union(a, b) => {
                if a.evaluate(point) <= b.evaluate(point) {
                    a.gradient(point)
                } else {
                    b.gradient(point)
                }
            }
            SdfNode::Intersection(a, b) => {
                if a.evaluate(point) >= b.evaluate(point) {
                    a.gradient(point)
                } else {
                    b.gradient(point)
                }
            }
            SdfNode::Difference(a, b) => {
                let da = a.evaluate(point);
                let db_neg = -b.evaluate(point);
                if da > db_neg {
                    a.gradient(point)
                } else {
                    -b.gradient(point)
                }
            }

            // Smooth CSG: central differences (blending makes analytical complex)
            SdfNode::SmoothUnion(_, _, _)
            | SdfNode::SmoothIntersection(_, _, _)
            | SdfNode::SmoothDifference(_, _, _) => central_diff_gradient(self, point),

            // -- Transforms ------------------------------------------------
            SdfNode::Translate(inner, offset) => inner.gradient(point - offset),
            SdfNode::Rotate(inner, rotation) => {
                let local_grad = inner.gradient(rotation.inverse() * point);
                let g = rotation * local_grad;
                let len = g.norm();
                if len > 1e-10 {
                    g / len
                } else {
                    Vector3::new(0.0, 1.0, 0.0)
                }
            }
            SdfNode::Scale(inner, factor) => {
                let g = inner.gradient(point / *factor);
                let len = g.norm();
                if len > 1e-10 {
                    g / len
                } else {
                    Vector3::new(0.0, 1.0, 0.0)
                }
            }
            SdfNode::Mirror(inner, normal) => {
                let d = point.dot(normal);
                let reflected = point - normal * (2.0 * d.min(0.0));
                let g = inner.gradient(reflected);
                if d < 0.0 {
                    // Reflect gradient back through the mirror plane
                    let g_ref = g - normal * (2.0 * g.dot(normal));
                    let len = g_ref.norm();
                    if len > 1e-10 {
                        g_ref / len
                    } else {
                        Vector3::new(0.0, 1.0, 0.0)
                    }
                } else {
                    g
                }
            }
            SdfNode::Shell(inner, _thickness) => {
                if inner.evaluate(point) >= 0.0 {
                    inner.gradient(point)
                } else {
                    -inner.gradient(point)
                }
            }
            SdfNode::Round(inner, _radius) => inner.gradient(point),

            // -- 2D -> 3D --------------------------------------------------
            SdfNode::Revolve(_) | SdfNode::Extrude(_, _) => central_diff_gradient(self, point),

            // -- Opaque ----------------------------------------------------
            SdfNode::Custom(_) => central_diff_gradient(self, point),
        }
    }

    // -----------------------------------------------------------------------
    // Feature enumeration
    // -----------------------------------------------------------------------

    /// Get face information for this primitive.
    /// Returns None for non-primitive nodes.
    pub fn face_info(&self) -> Option<Vec<FaceInfo>> {
        match self {
            SdfNode::Box3 { .. } => Some(vec![
                FaceInfo {
                    index: 0,
                    label: "+X".into(),
                    normal: Some(Vector3::new(1.0, 0.0, 0.0)),
                },
                FaceInfo {
                    index: 1,
                    label: "-X".into(),
                    normal: Some(Vector3::new(-1.0, 0.0, 0.0)),
                },
                FaceInfo {
                    index: 2,
                    label: "+Y".into(),
                    normal: Some(Vector3::new(0.0, 1.0, 0.0)),
                },
                FaceInfo {
                    index: 3,
                    label: "-Y".into(),
                    normal: Some(Vector3::new(0.0, -1.0, 0.0)),
                },
                FaceInfo {
                    index: 4,
                    label: "+Z".into(),
                    normal: Some(Vector3::new(0.0, 0.0, 1.0)),
                },
                FaceInfo {
                    index: 5,
                    label: "-Z".into(),
                    normal: Some(Vector3::new(0.0, 0.0, -1.0)),
                },
            ]),
            SdfNode::Cylinder { .. } => Some(vec![
                FaceInfo {
                    index: 0,
                    label: "side".into(),
                    normal: None,
                },
                FaceInfo {
                    index: 1,
                    label: "top".into(),
                    normal: None, // depends on axis direction
                },
                FaceInfo {
                    index: 2,
                    label: "bottom".into(),
                    normal: None, // depends on axis direction
                },
            ]),
            SdfNode::Sphere { .. } => Some(vec![FaceInfo {
                index: 0,
                label: "surface".into(),
                normal: None,
            }]),
            // Transforms delegate to inner
            SdfNode::Translate(inner, _) | SdfNode::Rotate(inner, _) | SdfNode::Scale(inner, _) => {
                inner.face_info()
            }
            _ => None,
        }
    }

    /// Get edge information for this primitive.
    /// Returns None for non-primitive nodes.
    pub fn edge_info(&self) -> Option<Vec<EdgeInfo>> {
        match self {
            SdfNode::Box3 { .. } => Some(vec![
                EdgeInfo {
                    index: 0,
                    face_a: 0,
                    face_b: 2,
                    label: "+X/+Y".into(),
                },
                EdgeInfo {
                    index: 1,
                    face_a: 0,
                    face_b: 3,
                    label: "+X/-Y".into(),
                },
                EdgeInfo {
                    index: 2,
                    face_a: 0,
                    face_b: 4,
                    label: "+X/+Z".into(),
                },
                EdgeInfo {
                    index: 3,
                    face_a: 0,
                    face_b: 5,
                    label: "+X/-Z".into(),
                },
                EdgeInfo {
                    index: 4,
                    face_a: 1,
                    face_b: 2,
                    label: "-X/+Y".into(),
                },
                EdgeInfo {
                    index: 5,
                    face_a: 1,
                    face_b: 3,
                    label: "-X/-Y".into(),
                },
                EdgeInfo {
                    index: 6,
                    face_a: 1,
                    face_b: 4,
                    label: "-X/+Z".into(),
                },
                EdgeInfo {
                    index: 7,
                    face_a: 1,
                    face_b: 5,
                    label: "-X/-Z".into(),
                },
                EdgeInfo {
                    index: 8,
                    face_a: 2,
                    face_b: 4,
                    label: "+Y/+Z".into(),
                },
                EdgeInfo {
                    index: 9,
                    face_a: 2,
                    face_b: 5,
                    label: "+Y/-Z".into(),
                },
                EdgeInfo {
                    index: 10,
                    face_a: 3,
                    face_b: 4,
                    label: "-Y/+Z".into(),
                },
                EdgeInfo {
                    index: 11,
                    face_a: 3,
                    face_b: 5,
                    label: "-Y/-Z".into(),
                },
            ]),
            SdfNode::Cylinder { .. } => Some(vec![
                EdgeInfo {
                    index: 0,
                    face_a: 0,
                    face_b: 1,
                    label: "side/top".into(),
                },
                EdgeInfo {
                    index: 1,
                    face_a: 0,
                    face_b: 2,
                    label: "side/bottom".into(),
                },
            ]),
            SdfNode::Sphere { .. } => Some(vec![]),
            // Transforms delegate to inner
            SdfNode::Translate(inner, _) | SdfNode::Rotate(inner, _) | SdfNode::Scale(inner, _) => {
                inner.edge_info()
            }
            _ => None,
        }
    }

    /// Find the closest face index at the given point.
    /// Returns None for non-primitive nodes.
    pub fn closest_face(&self, point: Vector3<f64>) -> Option<usize> {
        match self {
            SdfNode::Box3 {
                center,
                half_extents,
            } => {
                let d = point - center;
                // Normalized distance to each face pair
                let nx = d.x.abs() / half_extents.x;
                let ny = d.y.abs() / half_extents.y;
                let nz = d.z.abs() / half_extents.z;
                // The face whose normalized distance is largest is the closest face
                if nx >= ny && nx >= nz {
                    // X faces: 0 for +X, 1 for -X
                    if d.x >= 0.0 { Some(0) } else { Some(1) }
                } else if ny >= nz {
                    // Y faces: 2 for +Y, 3 for -Y
                    if d.y >= 0.0 { Some(2) } else { Some(3) }
                } else {
                    // Z faces: 4 for +Z, 5 for -Z
                    if d.z >= 0.0 { Some(4) } else { Some(5) }
                }
            }
            SdfNode::Cylinder {
                base, axis, height, ..
            } => {
                // Project point onto cylinder axis
                let rel = point - base;
                let t = rel.dot(axis);
                // Distance along axis: 0 = bottom, height = top
                // Radial distance handled by face_distance
                if t < 0.0 {
                    Some(2) // below bottom
                } else if t > *height {
                    Some(1) // above top
                } else {
                    // Compare radial vs axial proximity
                    let axial_to_bottom = t;
                    let axial_to_top = height - t;
                    let axial_min = axial_to_bottom.min(axial_to_top);
                    let radial = (rel - axis * t).norm();
                    if radial > axial_min {
                        Some(0) // side
                    } else if axial_to_top < axial_to_bottom {
                        Some(1) // top
                    } else {
                        Some(2) // bottom
                    }
                }
            }
            SdfNode::Sphere { .. } => Some(0),
            // Transforms: transform point and delegate
            SdfNode::Translate(inner, offset) => inner.closest_face(point - offset),
            SdfNode::Rotate(inner, rotation) => inner.closest_face(rotation.inverse() * point),
            SdfNode::Scale(inner, factor) => inner.closest_face(point / *factor),
            _ => None,
        }
    }

    /// Compute the signed distance from a point to a specific face.
    /// Returns None if the face doesn't exist or this isn't a primitive.
    pub fn face_distance(&self, point: Vector3<f64>, face_index: usize) -> Option<f64> {
        match self {
            SdfNode::Box3 {
                center,
                half_extents,
            } => {
                let d = point - center;
                match face_index {
                    0 => Some(d.x - half_extents.x),  // +X
                    1 => Some(-d.x - half_extents.x), // -X
                    2 => Some(d.y - half_extents.y),  // +Y
                    3 => Some(-d.y - half_extents.y), // -Y
                    4 => Some(d.z - half_extents.z),  // +Z
                    5 => Some(-d.z - half_extents.z), // -Z
                    _ => None,
                }
            }
            SdfNode::Cylinder {
                base,
                axis,
                radius,
                height,
            } => {
                let rel = point - base;
                let t = rel.dot(axis);
                match face_index {
                    0 => {
                        // Side: radial distance from axis
                        let proj = rel - axis * t;
                        Some(proj.norm() - radius)
                    }
                    1 => {
                        // Top cap: signed distance above top
                        Some(t - height)
                    }
                    2 => {
                        // Bottom cap: signed distance below bottom
                        Some(-t)
                    }
                    _ => None,
                }
            }
            SdfNode::Sphere { center, radius } => {
                if face_index == 0 {
                    Some((point - center).norm() - radius)
                } else {
                    None
                }
            }
            // Transforms: transform point and delegate
            SdfNode::Translate(inner, offset) => inner.face_distance(point - offset, face_index),
            SdfNode::Rotate(inner, rotation) => {
                inner.face_distance(rotation.inverse() * point, face_index)
            }
            SdfNode::Scale(inner, factor) => inner
                .face_distance(point / *factor, face_index)
                .map(|d| d * factor),
            _ => None,
        }
    }

    /// Conservative interval from evaluating all 8 corners of the bbox.
    /// The true range over the box is contained within the range of corner
    /// values (for Lipschitz-1 SDFs), expanded by the diagonal / 2 for safety.
    fn interval_from_corners(&self, bbox: &BBox3) -> Interval {
        let corners = bbox.corners();
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for c in &corners {
            let v = self.evaluate(*c);
            lo = lo.min(v);
            hi = hi.max(v);
        }
        // Expand by half-diagonal to be conservative (Lipschitz bound)
        let half_diag = (bbox.max - bbox.min).norm() * 0.5;
        Interval::new(lo - half_diag, hi + half_diag)
    }
}

// ---------------------------------------------------------------------------
// Central differences fallback for gradient computation
// ---------------------------------------------------------------------------

/// Compute the gradient via central finite differences.
///
/// Used for complex primitives (Box3, Cylinder, etc.), smooth CSG, 2D->3D
/// operations, and opaque Custom nodes where analytical gradients are either
/// impractical or error-prone.
fn central_diff_gradient(node: &SdfNode, point: Vector3<f64>) -> Vector3<f64> {
    let eps = 1e-6;
    let dx = node.evaluate(point + Vector3::new(eps, 0.0, 0.0))
        - node.evaluate(point - Vector3::new(eps, 0.0, 0.0));
    let dy = node.evaluate(point + Vector3::new(0.0, eps, 0.0))
        - node.evaluate(point - Vector3::new(0.0, eps, 0.0));
    let dz = node.evaluate(point + Vector3::new(0.0, 0.0, eps))
        - node.evaluate(point - Vector3::new(0.0, 0.0, eps));
    let g = Vector3::new(dx, dy, dz);
    let len = g.norm();
    if len > 1e-10 {
        g / len
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    }
}
