use nalgebra::{Rotation3, Vector2, Vector3};
use crate::primitives;
use crate::csg;

/// Trait for any object that can evaluate a signed distance.
pub trait Sdf: Send + Sync {
    fn evaluate(&self, point: Vector3<f64>) -> f64;
}

/// Trait for a 2D signed distance function.
/// Used as the cross-section profile for Revolve and Extrude operations.
pub trait Sdf2d: Send + Sync {
    fn evaluate(&self, point: Vector2<f64>) -> f64;
}

// ---------------------------------------------------------------------------
// Primitives
// ---------------------------------------------------------------------------

pub struct Sphere {
    pub center: Vector3<f64>,
    pub radius: f64,
}

impl Sphere {
    pub fn new(center: Vector3<f64>, radius: f64) -> Self {
        Self { center, radius }
    }
}

impl Sdf for Sphere {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        primitives::sdf_sphere(point, self.center, self.radius)
    }
}

pub struct Box3 {
    pub center: Vector3<f64>,
    pub half_extents: Vector3<f64>,
}

impl Box3 {
    pub fn new(center: Vector3<f64>, half_extents: Vector3<f64>) -> Self {
        Self { center, half_extents }
    }
}

impl Sdf for Box3 {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        primitives::sdf_box(point, self.center, self.half_extents)
    }
}

/// Exact capped cone (or truncated cone). Flat base at `a` with radius `ra`,
/// tip/cap at `b` with radius `rb`. Set `rb = 0` for a sharp cone.
pub struct CappedCone {
    pub a: Vector3<f64>,
    pub b: Vector3<f64>,
    pub ra: f64,
    pub rb: f64,
}

impl CappedCone {
    pub fn new(a: Vector3<f64>, b: Vector3<f64>, ra: f64, rb: f64) -> Self {
        Self { a, b, ra, rb }
    }
}

impl Sdf for CappedCone {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        primitives::sdf_capped_cone(point, self.a, self.b, self.ra, self.rb)
    }
}

/// Cylinder shape wrapper.
pub struct Cylinder {
    pub base: Vector3<f64>,
    pub axis: Vector3<f64>,
    pub radius: f64,
    pub height: f64,
}

impl Cylinder {
    pub fn new(base: Vector3<f64>, axis: Vector3<f64>, radius: f64, height: f64) -> Self {
        Self { base, axis: axis.normalize(), radius, height }
    }
}

impl Sdf for Cylinder {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        primitives::sdf_cylinder(point, self.base, self.axis, self.radius, self.height)
    }
}

/// Torus lying in the XZ plane, centered at `center`.
pub struct Torus {
    pub center: Vector3<f64>,
    pub major_radius: f64,
    pub minor_radius: f64,
}

impl Torus {
    pub fn new(center: Vector3<f64>, major_radius: f64, minor_radius: f64) -> Self {
        Self { center, major_radius, minor_radius }
    }
}

impl Sdf for Torus {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        primitives::sdf_torus(point, self.center, self.major_radius, self.minor_radius)
    }
}

/// Box with rounded edges. Total size is `half_extents + radius`.
pub struct RoundedBox {
    pub center: Vector3<f64>,
    pub half_extents: Vector3<f64>,
    pub radius: f64,
}

impl RoundedBox {
    pub fn new(center: Vector3<f64>, half_extents: Vector3<f64>, radius: f64) -> Self {
        Self { center, half_extents, radius }
    }
}

impl Sdf for RoundedBox {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        primitives::sdf_rounded_box(point, self.center, self.half_extents, self.radius)
    }
}

/// Capsule (sphere-swept line segment) from `a` to `b` with given `radius`.
pub struct Capsule {
    pub a: Vector3<f64>,
    pub b: Vector3<f64>,
    pub radius: f64,
}

impl Capsule {
    pub fn new(a: Vector3<f64>, b: Vector3<f64>, radius: f64) -> Self {
        Self { a, b, radius }
    }
}

impl Sdf for Capsule {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        primitives::sdf_capsule(point, self.a, self.b, self.radius)
    }
}

/// Ellipsoid with semi-axis lengths `radii`. Approximate SDF.
pub struct Ellipsoid {
    pub center: Vector3<f64>,
    pub radii: Vector3<f64>,
}

impl Ellipsoid {
    pub fn new(center: Vector3<f64>, radii: Vector3<f64>) -> Self {
        Self { center, radii }
    }
}

impl Sdf for Ellipsoid {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        primitives::sdf_ellipsoid(point, self.center, self.radii)
    }
}

/// Rounded cylinder aligned along Y axis.
pub struct RoundedCylinder {
    pub center: Vector3<f64>,
    pub radius: f64,
    pub round_radius: f64,
    pub half_height: f64,
}

impl RoundedCylinder {
    pub fn new(center: Vector3<f64>, radius: f64, round_radius: f64, half_height: f64) -> Self {
        Self { center, radius, round_radius, half_height }
    }
}

impl Sdf for RoundedCylinder {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        primitives::sdf_rounded_cylinder(
            point, self.center, self.radius, self.round_radius, self.half_height,
        )
    }
}

/// Half-space SDF: the solid occupies the region behind the plane.
/// The plane is defined by `normal . p + d = 0` where `normal` points
/// outward (into empty space). Inside is `normal . p + d < 0`.
pub struct HalfSpace {
    pub normal: Vector3<f64>,
    pub d: f64,
}

impl HalfSpace {
    pub fn new(normal: Vector3<f64>, d: f64) -> Self {
        Self { normal, d }
    }
}

impl Sdf for HalfSpace {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        self.normal.dot(&point) + self.d
    }
}

// ---------------------------------------------------------------------------
// CSG Operations
// ---------------------------------------------------------------------------

pub struct Union<A: Sdf, B: Sdf> {
    pub a: A,
    pub b: B,
}

impl<A: Sdf, B: Sdf> Union<A, B> {
    pub fn new(a: A, b: B) -> Self { Self { a, b } }
}

impl<A: Sdf, B: Sdf> Sdf for Union<A, B> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        csg::union(self.a.evaluate(point), self.b.evaluate(point))
    }
}

pub struct Intersection<A: Sdf, B: Sdf> {
    pub a: A,
    pub b: B,
}

impl<A: Sdf, B: Sdf> Intersection<A, B> {
    pub fn new(a: A, b: B) -> Self { Self { a, b } }
}

impl<A: Sdf, B: Sdf> Sdf for Intersection<A, B> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        csg::intersection(self.a.evaluate(point), self.b.evaluate(point))
    }
}

pub struct Difference<A: Sdf, B: Sdf> {
    pub a: A,
    pub b: B,
}

impl<A: Sdf, B: Sdf> Difference<A, B> {
    pub fn new(a: A, b: B) -> Self { Self { a, b } }
}

impl<A: Sdf, B: Sdf> Sdf for Difference<A, B> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        csg::difference(self.a.evaluate(point), self.b.evaluate(point))
    }
}

/// Smooth union with blending radius `k`.
pub struct SmoothUnion<A: Sdf, B: Sdf> {
    pub a: A,
    pub b: B,
    pub k: f64,
}

impl<A: Sdf, B: Sdf> SmoothUnion<A, B> {
    pub fn new(a: A, b: B, k: f64) -> Self { Self { a, b, k } }
}

impl<A: Sdf, B: Sdf> Sdf for SmoothUnion<A, B> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        csg::smooth_union(self.a.evaluate(point), self.b.evaluate(point), self.k)
    }
}

/// Smooth intersection with blending radius `k`.
pub struct SmoothIntersection<A: Sdf, B: Sdf> {
    pub a: A,
    pub b: B,
    pub k: f64,
}

impl<A: Sdf, B: Sdf> SmoothIntersection<A, B> {
    pub fn new(a: A, b: B, k: f64) -> Self { Self { a, b, k } }
}

impl<A: Sdf, B: Sdf> Sdf for SmoothIntersection<A, B> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        csg::smooth_intersection(self.a.evaluate(point), self.b.evaluate(point), self.k)
    }
}

/// Smooth difference with blending radius `k`.
pub struct SmoothDifference<A: Sdf, B: Sdf> {
    pub a: A,
    pub b: B,
    pub k: f64,
}

impl<A: Sdf, B: Sdf> SmoothDifference<A, B> {
    pub fn new(a: A, b: B, k: f64) -> Self { Self { a, b, k } }
}

impl<A: Sdf, B: Sdf> Sdf for SmoothDifference<A, B> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        csg::smooth_difference(self.a.evaluate(point), self.b.evaluate(point), self.k)
    }
}

// ---------------------------------------------------------------------------
// Transforms
// ---------------------------------------------------------------------------

/// Translate a shape by an offset. Inverse: shift the query point backwards.
pub struct Translate<S: Sdf> {
    pub shape: S,
    pub offset: Vector3<f64>,
}

impl<S: Sdf> Translate<S> {
    pub fn new(shape: S, offset: Vector3<f64>) -> Self { Self { shape, offset } }
}

impl<S: Sdf> Sdf for Translate<S> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        self.shape.evaluate(point - self.offset)
    }
}

/// Rotate a shape. Inverse: rotate the query point by the inverse rotation.
pub struct Rotate<S: Sdf> {
    pub shape: S,
    pub rotation: Rotation3<f64>,
    inv: Rotation3<f64>,
}

impl<S: Sdf> Rotate<S> {
    pub fn new(shape: S, rotation: Rotation3<f64>) -> Self {
        let inv = rotation.inverse();
        Self { shape, rotation, inv }
    }
}

impl<S: Sdf> Sdf for Rotate<S> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        self.shape.evaluate(self.inv * point)
    }
}

/// Uniform scale. Non-uniform scaling is not supported (breaks distance property).
/// The SDF value must be scaled by the same factor.
pub struct Scale<S: Sdf> {
    pub shape: S,
    pub factor: f64,
    inv: f64,
}

impl<S: Sdf> Scale<S> {
    pub fn new(shape: S, factor: f64) -> Self {
        Self { shape, factor, inv: 1.0 / factor }
    }
}

impl<S: Sdf> Sdf for Scale<S> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        self.shape.evaluate(point * self.inv) * self.factor
    }
}

/// Mirror a shape across a plane through the origin with given normal.
pub struct Mirror<S: Sdf> {
    pub shape: S,
    pub normal: Vector3<f64>,
}

impl<S: Sdf> Mirror<S> {
    pub fn new(shape: S, normal: Vector3<f64>) -> Self {
        Self { shape, normal: normal.normalize() }
    }
}

impl<S: Sdf> Sdf for Mirror<S> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        // Reflect point across the mirror plane
        let d = point.dot(&self.normal);
        let reflected = point - self.normal * (2.0 * d.min(0.0));
        self.shape.evaluate(reflected)
    }
}

/// Shell (onion) operation: hollows a shape with wall thickness.
pub struct Shell<S: Sdf> {
    pub shape: S,
    pub thickness: f64,
}

impl<S: Sdf> Shell<S> {
    pub fn new(shape: S, thickness: f64) -> Self { Self { shape, thickness } }
}

impl<S: Sdf> Sdf for Shell<S> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        self.shape.evaluate(point).abs() - self.thickness
    }
}

/// Round (offset) operation: expands or contracts a shape's surface.
pub struct Round<S: Sdf> {
    pub shape: S,
    pub radius: f64,
}

impl<S: Sdf> Round<S> {
    pub fn new(shape: S, radius: f64) -> Self { Self { shape, radius } }
}

impl<S: Sdf> Sdf for Round<S> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        self.shape.evaluate(point) - self.radius
    }
}

// ---------------------------------------------------------------------------
// Closure-based SDF (for dynamic/user-defined shapes)
// ---------------------------------------------------------------------------

/// An SDF defined by a closure. Useful for custom/dynamic shapes.
pub struct FnSdf<F: Fn(Vector3<f64>) -> f64 + Send + Sync> {
    pub func: F,
}

impl<F: Fn(Vector3<f64>) -> f64 + Send + Sync> FnSdf<F> {
    pub fn new(func: F) -> Self { Self { func } }
}

impl<F: Fn(Vector3<f64>) -> f64 + Send + Sync> Sdf for FnSdf<F> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        (self.func)(point)
    }
}

/// A 2D SDF defined by a closure.
pub struct FnSdf2d<F: Fn(Vector2<f64>) -> f64 + Send + Sync> {
    pub func: F,
}

impl<F: Fn(Vector2<f64>) -> f64 + Send + Sync> FnSdf2d<F> {
    pub fn new(func: F) -> Self { Self { func } }
}

impl<F: Fn(Vector2<f64>) -> f64 + Send + Sync> Sdf2d for FnSdf2d<F> {
    fn evaluate(&self, point: Vector2<f64>) -> f64 {
        (self.func)(point)
    }
}

// ---------------------------------------------------------------------------
// 2D Primitives (for use with Revolve and Extrude)
// ---------------------------------------------------------------------------

/// 2D circle SDF.
pub struct Circle2d {
    pub center: Vector2<f64>,
    pub radius: f64,
}

impl Circle2d {
    pub fn new(center: Vector2<f64>, radius: f64) -> Self { Self { center, radius } }
}

impl Sdf2d for Circle2d {
    fn evaluate(&self, point: Vector2<f64>) -> f64 {
        (point - self.center).norm() - self.radius
    }
}

/// 2D rectangle SDF.
pub struct Rect2d {
    pub center: Vector2<f64>,
    pub half_extents: Vector2<f64>,
}

impl Rect2d {
    pub fn new(center: Vector2<f64>, half_extents: Vector2<f64>) -> Self {
        Self { center, half_extents }
    }
}

impl Sdf2d for Rect2d {
    fn evaluate(&self, point: Vector2<f64>) -> f64 {
        let d = (point - self.center).abs() - self.half_extents;
        let outside = Vector2::new(d.x.max(0.0), d.y.max(0.0)).norm();
        let inside = d.x.max(d.y).min(0.0);
        outside + inside
    }
}

/// 2D union.
pub struct Union2d<A: Sdf2d, B: Sdf2d> {
    pub a: A,
    pub b: B,
}

impl<A: Sdf2d, B: Sdf2d> Union2d<A, B> {
    pub fn new(a: A, b: B) -> Self { Self { a, b } }
}

impl<A: Sdf2d, B: Sdf2d> Sdf2d for Union2d<A, B> {
    fn evaluate(&self, point: Vector2<f64>) -> f64 {
        self.a.evaluate(point).min(self.b.evaluate(point))
    }
}

/// 2D difference.
pub struct Difference2d<A: Sdf2d, B: Sdf2d> {
    pub a: A,
    pub b: B,
}

impl<A: Sdf2d, B: Sdf2d> Difference2d<A, B> {
    pub fn new(a: A, b: B) -> Self { Self { a, b } }
}

impl<A: Sdf2d, B: Sdf2d> Sdf2d for Difference2d<A, B> {
    fn evaluate(&self, point: Vector2<f64>) -> f64 {
        self.a.evaluate(point).max(-self.b.evaluate(point))
    }
}

// ---------------------------------------------------------------------------
// Revolve and Extrude (2D → 3D)
// ---------------------------------------------------------------------------

/// Revolve a 2D profile around the Y axis to create a solid of revolution.
///
/// The 2D profile is evaluated in the (r, y) plane where r = sqrt(x² + z²).
/// This maps every 3D point to its corresponding position on the profile,
/// sweeping the profile 360 degrees around Y.
///
/// The profile's X coordinate corresponds to the radial distance from Y axis,
/// and the profile's Y coordinate corresponds to the world Y coordinate.
pub struct Revolve<P: Sdf2d> {
    pub profile: P,
}

impl<P: Sdf2d> Revolve<P> {
    pub fn new(profile: P) -> Self { Self { profile } }
}

impl<P: Sdf2d> Sdf for Revolve<P> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        let r = (point.x * point.x + point.z * point.z).sqrt();
        self.profile.evaluate(Vector2::new(r, point.y))
    }
}

/// Extrude a 2D profile along the Z axis.
///
/// The 2D profile is evaluated in the XY plane. The extrusion extends
/// `half_height` in both the +Z and -Z directions.
pub struct Extrude<P: Sdf2d> {
    pub profile: P,
    pub half_height: f64,
}

impl<P: Sdf2d> Extrude<P> {
    pub fn new(profile: P, half_height: f64) -> Self { Self { profile, half_height } }
}

impl<P: Sdf2d> Sdf for Extrude<P> {
    fn evaluate(&self, point: Vector3<f64>) -> f64 {
        let d_2d = self.profile.evaluate(Vector2::new(point.x, point.y));
        let d_z = point.z.abs() - self.half_height;
        if d_2d > 0.0 && d_z > 0.0 {
            (d_2d * d_2d + d_z * d_z).sqrt()
        } else {
            d_2d.max(d_z)
        }
    }
}
