use nalgebra::Vector3;
use std::f64::consts::PI;

/// A parametric path P(t) where t in [0, 1].
pub trait Path: Send + Sync {
    /// Position at parameter t.
    fn point(&self, t: f64) -> Vector3<f64>;

    /// Unit tangent at parameter t.
    fn tangent(&self, t: f64) -> Vector3<f64> {
        let eps = 1e-6;
        let t0 = (t - eps).max(0.0);
        let t1 = (t + eps).min(1.0);
        let dir = self.point(t1) - self.point(t0);
        let len = dir.norm();
        if len > 1e-10 {
            dir / len
        } else {
            Vector3::new(0.0, 0.0, 1.0)
        }
    }
}

/// Straight line from start to end.
pub struct LinePath {
    pub start: Vector3<f64>,
    pub end: Vector3<f64>,
}

impl LinePath {
    pub fn new(start: Vector3<f64>, end: Vector3<f64>) -> Self {
        Self { start, end }
    }
}

impl Path for LinePath {
    fn point(&self, t: f64) -> Vector3<f64> {
        self.start + (self.end - self.start) * t
    }

    fn tangent(&self, _t: f64) -> Vector3<f64> {
        let dir = self.end - self.start;
        let len = dir.norm();
        if len > 1e-10 {
            dir / len
        } else {
            Vector3::new(0.0, 0.0, 1.0)
        }
    }
}

/// Helix with constant radius, pitch, and number of turns.
/// The eigenform of Order 2 (constant curvature and torsion).
pub struct HelixPath {
    pub radius: f64,
    pub pitch: f64,
    pub turns: f64,
}

impl HelixPath {
    pub fn new(radius: f64, pitch: f64, turns: f64) -> Self {
        Self {
            radius,
            pitch,
            turns,
        }
    }
}

impl Path for HelixPath {
    fn point(&self, t: f64) -> Vector3<f64> {
        let angle = t * 2.0 * PI * self.turns;
        let height = t * self.pitch * self.turns;
        Vector3::new(self.radius * angle.cos(), self.radius * angle.sin(), height)
    }
}

/// Spiral with variable radius and height, parameterized by closures.
/// Used for Order 3 eigenform (horn/tusk) paths.
pub struct SpiralPath {
    radius_fn: Box<dyn Fn(f64) -> f64 + Send + Sync>,
    height_fn: Box<dyn Fn(f64) -> f64 + Send + Sync>,
    turns: f64,
}

impl SpiralPath {
    pub fn new(
        radius_fn: impl Fn(f64) -> f64 + Send + Sync + 'static,
        height_fn: impl Fn(f64) -> f64 + Send + Sync + 'static,
        turns: f64,
    ) -> Self {
        Self {
            radius_fn: Box::new(radius_fn),
            height_fn: Box::new(height_fn),
            turns,
        }
    }
}

impl Path for SpiralPath {
    fn point(&self, t: f64) -> Vector3<f64> {
        let angle = t * 2.0 * PI * self.turns;
        let r = (self.radius_fn)(t);
        let h = (self.height_fn)(t);
        Vector3::new(r * angle.cos(), r * angle.sin(), h)
    }
}
