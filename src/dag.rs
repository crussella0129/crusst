//! SdfNode expression DAG — an enum-based composable SDF tree.
//!
//! Unlike the trait-based API in `shape.rs` (which uses generics + type erasure),
//! the DAG uses `Arc<SdfNode>` for children. This enables:
//! - Interval arithmetic for octree pruning
//! - Analytical gradients for QEF
//! - DAG introspection for smart STEP export

use std::sync::Arc;
use nalgebra::{Rotation3, Vector2, Vector3};
use crate::{csg, primitives};

// ---------------------------------------------------------------------------
// SdfNode2d — 2D profile sub-enum (for Revolve / Extrude)
// ---------------------------------------------------------------------------

/// A 2D signed distance function node, used as the cross-section profile
/// for `SdfNode::Revolve` and `SdfNode::Extrude`.
pub enum SdfNode2d {
    /// 2D circle.
    Circle2d { center: Vector2<f64>, radius: f64 },
    /// 2D axis-aligned rectangle.
    Rect2d { center: Vector2<f64>, half_extents: Vector2<f64> },
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
            SdfNode2d::Circle2d { center, radius } => {
                (point - center).norm() - radius
            }
            SdfNode2d::Rect2d { center, half_extents } => {
                let d = (point - center).abs() - half_extents;
                let outside = Vector2::new(d.x.max(0.0), d.y.max(0.0)).norm();
                let inside = d.x.max(d.y).min(0.0);
                outside + inside
            }
            SdfNode2d::Union2d(a, b) => {
                a.evaluate(point).min(b.evaluate(point))
            }
            SdfNode2d::Difference2d(a, b) => {
                a.evaluate(point).max(-b.evaluate(point))
            }
            SdfNode2d::Custom2d(sdf) => {
                sdf.evaluate(point)
            }
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
    Box3 { center: Vector3<f64>, half_extents: Vector3<f64> },

    /// Capped cylinder with `base` center, unit `axis`, `radius`, and `height`.
    Cylinder { base: Vector3<f64>, axis: Vector3<f64>, radius: f64, height: f64 },

    /// Capped cone (truncated cone) from `a` (radius `ra`) to `b` (radius `rb`).
    CappedCone { a: Vector3<f64>, b: Vector3<f64>, ra: f64, rb: f64 },

    /// Torus lying in the XZ plane, centered at `center`.
    Torus { center: Vector3<f64>, major_radius: f64, minor_radius: f64 },

    /// Box with rounded edges. Total size is `half_extents + radius`.
    RoundedBox { center: Vector3<f64>, half_extents: Vector3<f64>, radius: f64 },

    /// Capsule (sphere-swept segment) from `a` to `b` with `radius`.
    Capsule { a: Vector3<f64>, b: Vector3<f64>, radius: f64 },

    /// Ellipsoid with semi-axis lengths `radii`. Approximate SDF.
    Ellipsoid { center: Vector3<f64>, radii: Vector3<f64> },

    /// Rounded cylinder aligned along Y axis.
    RoundedCylinder { center: Vector3<f64>, radius: f64, round_radius: f64, half_height: f64 },

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

impl SdfNode {
    /// Evaluate the signed distance at a 3D point.
    pub fn evaluate(&self, point: Vector3<f64>) -> f64 {
        match self {
            // -- Primitives ------------------------------------------------

            SdfNode::Sphere { center, radius } => {
                primitives::sdf_sphere(point, *center, *radius)
            }
            SdfNode::Box3 { center, half_extents } => {
                primitives::sdf_box(point, *center, *half_extents)
            }
            SdfNode::Cylinder { base, axis, radius, height } => {
                primitives::sdf_cylinder(point, *base, *axis, *radius, *height)
            }
            SdfNode::CappedCone { a, b, ra, rb } => {
                primitives::sdf_capped_cone(point, *a, *b, *ra, *rb)
            }
            SdfNode::Torus { center, major_radius, minor_radius } => {
                primitives::sdf_torus(point, *center, *major_radius, *minor_radius)
            }
            SdfNode::RoundedBox { center, half_extents, radius } => {
                primitives::sdf_rounded_box(point, *center, *half_extents, *radius)
            }
            SdfNode::Capsule { a, b, radius } => {
                primitives::sdf_capsule(point, *a, *b, *radius)
            }
            SdfNode::Ellipsoid { center, radii } => {
                primitives::sdf_ellipsoid(point, *center, *radii)
            }
            SdfNode::RoundedCylinder { center, radius, round_radius, half_height } => {
                primitives::sdf_rounded_cylinder(point, *center, *radius, *round_radius, *half_height)
            }
            SdfNode::HalfSpace { normal, d } => {
                normal.dot(&point) + d
            }

            // -- CSG -------------------------------------------------------

            SdfNode::Union(a, b) => {
                csg::union(a.evaluate(point), b.evaluate(point))
            }
            SdfNode::Intersection(a, b) => {
                csg::intersection(a.evaluate(point), b.evaluate(point))
            }
            SdfNode::Difference(a, b) => {
                csg::difference(a.evaluate(point), b.evaluate(point))
            }
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

            SdfNode::Translate(inner, offset) => {
                inner.evaluate(point - offset)
            }
            SdfNode::Rotate(inner, rotation) => {
                inner.evaluate(rotation.inverse() * point)
            }
            SdfNode::Scale(inner, factor) => {
                let inv = 1.0 / factor;
                inner.evaluate(point * inv) * factor
            }
            SdfNode::Mirror(inner, normal) => {
                let n = normal.normalize();
                let d = point.dot(&n);
                let reflected = point - n * (2.0 * d.min(0.0));
                inner.evaluate(reflected)
            }
            SdfNode::Shell(inner, thickness) => {
                inner.evaluate(point).abs() - thickness
            }
            SdfNode::Round(inner, radius) => {
                inner.evaluate(point) - radius
            }

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

            SdfNode::Custom(sdf) => {
                sdf.evaluate(point)
            }
        }
    }
}
