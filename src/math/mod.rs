//! Linear algebra type aliases and geometric tolerances.

pub type Point3 = nalgebra::Point3<f64>;
pub type Point2 = nalgebra::Point2<f64>;
pub type Vector3 = nalgebra::Vector3<f64>;
pub type Vector2 = nalgebra::Vector2<f64>;
pub type Matrix4 = nalgebra::Matrix4<f64>;

/// Geometric tolerance for point coincidence tests (distance in model units).
pub const TOLERANCE: f64 = 1e-9;

/// Parametric tolerance for curve/surface parameter comparisons.
pub const PARAM_TOL: f64 = 1e-12;

/// Angular tolerance (radians) for tangent/normal comparisons.
pub const ANGLE_TOL: f64 = 1e-6;
