//! 3D and 2D curve types for the B-Rep kernel.
//!
//! `Curve3` represents edges in 3D space.
//! `Curve2` represents parametric-space curves (PCurves) on face surfaces.

use crate::math::{Point2, Point3, Vector2, Vector3};
use crate::nurbs::{NurbsCurve2, NurbsCurve3};

/// A geometric curve in 3D space. Carried by topology `Edge` entities.
#[derive(Clone, Debug)]
pub enum Curve3 {
    Line {
        origin: Point3,
        dir: Vector3,
    },
    Circle {
        center: Point3,
        axis: Vector3,
        radius: f64,
    },
    Ellipse {
        center: Point3,
        major: Vector3,
        minor: Vector3,
    },
    NurbsCurve(NurbsCurve3),
}

impl Curve3 {
    /// Evaluate the curve at parameter `t`.
    pub fn evaluate(&self, t: f64) -> Point3 {
        match self {
            Curve3::Line { origin, dir } => origin + dir * t,
            Curve3::Circle {
                center,
                axis,
                radius,
            } => {
                let (u, v) = circle_frame(axis);
                center + (u * t.cos() + v * t.sin()) * *radius
            }
            Curve3::Ellipse {
                center,
                major,
                minor,
            } => center + major * t.cos() + minor * t.sin(),
            Curve3::NurbsCurve(nurbs) => nurbs.evaluate(t),
        }
    }

    /// First derivative (tangent direction, not necessarily unit length).
    pub fn derivative(&self, t: f64) -> Vector3 {
        match self {
            Curve3::Line { dir, .. } => *dir,
            Curve3::Circle {
                axis, radius, ..
            } => {
                let (u, v) = circle_frame(axis);
                (-u * t.sin() + v * t.cos()) * *radius
            }
            Curve3::Ellipse { major, minor, .. } => -major * t.sin() + minor * t.cos(),
            Curve3::NurbsCurve(nurbs) => {
                let (_, d) = nurbs.derivative(t);
                d
            }
        }
    }

    /// Unit tangent vector at parameter `t`.
    pub fn tangent(&self, t: f64) -> Vector3 {
        let d = self.derivative(t);
        let len = d.norm();
        if len > 1e-15 {
            d / len
        } else {
            Vector3::new(1.0, 0.0, 0.0)
        }
    }

    /// Find the parameter closest to the given point (point inversion).
    ///
    /// For lines, this is a direct projection. For others, uses Newton-Raphson.
    pub fn closest_parameter(&self, point: &Point3) -> f64 {
        match self {
            Curve3::Line { origin, dir } => {
                let d = point - origin;
                d.dot(dir) / dir.dot(dir)
            }
            Curve3::Circle {
                center,
                axis,
                radius: _,
            } => {
                let (u, v) = circle_frame(axis);
                let d = point - center;
                d.dot(&v).atan2(d.dot(&u))
            }
            Curve3::Ellipse {
                center,
                major,
                minor,
            } => {
                let d = point - center;
                d.dot(minor).atan2(d.dot(major))
            }
            Curve3::NurbsCurve(nurbs) => {
                // Newton-Raphson: minimize |C(t) - P|²
                // f(t) = C'(t) · (C(t) - P) = 0
                let (t_min, t_max) = nurbs.domain();
                let mut t = (t_min + t_max) * 0.5;

                // Coarse search first
                let n_samples = 20;
                let mut best_dist = f64::MAX;
                for i in 0..=n_samples {
                    let ti = t_min + (t_max - t_min) * i as f64 / n_samples as f64;
                    let dist = (nurbs.evaluate(ti) - point).norm();
                    if dist < best_dist {
                        best_dist = dist;
                        t = ti;
                    }
                }

                // Newton-Raphson refinement
                for _ in 0..50 {
                    let (c, dc) = nurbs.derivative(t);
                    let diff = c - point;
                    let f = dc.dot(&diff);
                    // f'(t) ≈ C'(t)·C'(t) + C''(t)·(C(t)-P) ≈ |C'|² (ignoring second derivative)
                    let df = dc.dot(&dc);
                    if df.abs() < 1e-30 {
                        break;
                    }
                    let dt = f / df;
                    t -= dt;
                    t = t.clamp(t_min, t_max);
                    if dt.abs() < 1e-12 {
                        break;
                    }
                }
                t
            }
        }
    }
}

/// A geometric curve in 2D parameter space. Used as PCurves for trimming.
#[derive(Clone, Debug)]
pub enum Curve2 {
    Line {
        origin: Point2,
        dir: Vector2,
    },
    Circle {
        center: Point2,
        radius: f64,
    },
    NurbsCurve(NurbsCurve2),
}

impl Curve2 {
    /// Evaluate the curve at parameter `t`.
    pub fn evaluate(&self, t: f64) -> Point2 {
        match self {
            Curve2::Line { origin, dir } => origin + dir * t,
            Curve2::Circle { center, radius } => {
                Point2::new(center.x + radius * t.cos(), center.y + radius * t.sin())
            }
            Curve2::NurbsCurve(nurbs) => nurbs.evaluate(t),
        }
    }

    /// First derivative.
    pub fn derivative(&self, t: f64) -> Vector2 {
        match self {
            Curve2::Line { dir, .. } => *dir,
            Curve2::Circle { radius, .. } => Vector2::new(-radius * t.sin(), radius * t.cos()),
            Curve2::NurbsCurve(_nurbs) => {
                // Finite difference for 2D NURBS (derivative not implemented on NurbsCurve2)
                let h = 1e-7;
                let p_plus = self.evaluate(t + h);
                let p_minus = self.evaluate(t - h);
                (p_plus - p_minus) / (2.0 * h)
            }
        }
    }

    /// Create a straight-line PCurve from `a` to `b`, parameterized on `[0, 1]`.
    pub fn line_segment(a: Point2, b: Point2) -> Self {
        Curve2::Line {
            origin: a,
            dir: b - a,
        }
    }
}

/// Compute an orthonormal frame (u, v) in the plane perpendicular to `axis`.
fn circle_frame(axis: &Vector3) -> (Vector3, Vector3) {
    let a = axis.normalize();
    // Pick a vector not parallel to axis
    let seed = if a.x.abs() < 0.9 {
        Vector3::new(1.0, 0.0, 0.0)
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    };
    let u = a.cross(&seed).normalize();
    let v = a.cross(&u);
    (u, v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn line_evaluate() {
        let line = Curve3::Line {
            origin: Point3::new(1.0, 0.0, 0.0),
            dir: Vector3::new(2.0, 0.0, 0.0),
        };
        let p = line.evaluate(0.5);
        assert!((p - Point3::new(2.0, 0.0, 0.0)).norm() < 1e-14);
    }

    #[test]
    fn line_derivative_is_constant() {
        let line = Curve3::Line {
            origin: Point3::origin(),
            dir: Vector3::new(3.0, 4.0, 0.0),
        };
        let d0 = line.derivative(0.0);
        let d1 = line.derivative(1.0);
        assert!((d0 - d1).norm() < 1e-14);
    }

    #[test]
    fn circle_evaluate_on_circle() {
        let circle = Curve3::Circle {
            center: Point3::origin(),
            axis: Vector3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        };
        for i in 0..=20 {
            let t = 2.0 * PI * i as f64 / 20.0;
            let p = circle.evaluate(t);
            let dist = p.coords.norm();
            assert!(
                (dist - 5.0).abs() < 1e-12,
                "Circle point at t={t} has radius {dist}"
            );
        }
    }

    #[test]
    fn circle_tangent_perpendicular_to_radius() {
        let circle = Curve3::Circle {
            center: Point3::origin(),
            axis: Vector3::new(0.0, 0.0, 1.0),
            radius: 3.0,
        };
        for i in 1..20 {
            let t = 2.0 * PI * i as f64 / 20.0;
            let p = circle.evaluate(t);
            let d = circle.derivative(t);
            let dot = p.coords.dot(&d);
            assert!(dot.abs() < 1e-10, "Tangent should be perpendicular at t={t}");
        }
    }

    #[test]
    fn line_closest_parameter() {
        let line = Curve3::Line {
            origin: Point3::origin(),
            dir: Vector3::new(1.0, 0.0, 0.0),
        };
        let p = Point3::new(3.0, 5.0, 0.0);
        let t = line.closest_parameter(&p);
        assert!((t - 3.0).abs() < 1e-14);
    }

    #[test]
    fn circle_closest_parameter_roundtrip() {
        let circle = Curve3::Circle {
            center: Point3::origin(),
            axis: Vector3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        };
        // Evaluate at a known parameter, then invert — should recover the same parameter.
        let t_orig = 1.2;
        let p = circle.evaluate(t_orig);
        let t_inv = circle.closest_parameter(&p);
        assert!(
            (t_orig - t_inv).abs() < 1e-10,
            "Expected {t_orig}, got {t_inv}"
        );
    }

    #[test]
    fn curve2_line_segment() {
        let c = Curve2::line_segment(Point2::new(0.0, 0.0), Point2::new(1.0, 1.0));
        let p = c.evaluate(0.5);
        assert!((p - Point2::new(0.5, 0.5)).norm() < 1e-14);
    }

    #[test]
    fn ellipse_evaluate() {
        let ellipse = Curve3::Ellipse {
            center: Point3::origin(),
            major: Vector3::new(3.0, 0.0, 0.0),
            minor: Vector3::new(0.0, 2.0, 0.0),
        };
        let p0 = ellipse.evaluate(0.0);
        assert!((p0 - Point3::new(3.0, 0.0, 0.0)).norm() < 1e-14);
        let p_half = ellipse.evaluate(FRAC_PI_2);
        assert!((p_half - Point3::new(0.0, 2.0, 0.0)).norm() < 1e-12);
    }
}
