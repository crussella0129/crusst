use nalgebra::Vector3;
use crate::primitives;
use crate::csg;

/// Trait for any object that can evaluate a signed distance.
pub trait Sdf: Send + Sync {
    fn evaluate(&self, point: Vector3<f64>) -> f64;
}

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

/// Half-space SDF: the solid occupies the region behind the plane.
/// The plane is defined by `normal · p + d = 0` where `normal` points
/// outward (into empty space). Inside is `normal · p + d < 0`.
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
