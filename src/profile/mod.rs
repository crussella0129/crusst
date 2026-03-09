//! 2D sketch profiles for extrude, revolve, and sweep operations.
//!
//! A `Profile` is a closed 2D wire in the XY plane, composed of line
//! and arc segments. Profiles serve as the cross-section for shape-generating
//! operations.

use crate::math::{Point2, Vector2};

/// A 2D sketch profile — a closed wire of line/arc segments.
#[derive(Clone, Debug)]
pub struct Profile {
    /// Ordered list of segments forming the closed profile.
    pub segments: Vec<Segment2>,
}

/// A 2D segment in a profile.
#[derive(Clone, Debug)]
pub enum Segment2 {
    Line { start: Point2, end: Point2 },
    Arc { center: Point2, radius: f64, start_angle: f64, end_angle: f64 },
}

impl Profile {
    /// Create a rectangular profile centered at the origin.
    pub fn rect(width: f64, height: f64) -> Self {
        let hw = width / 2.0;
        let hh = height / 2.0;
        let p0 = Point2::new(-hw, -hh);
        let p1 = Point2::new(hw, -hh);
        let p2 = Point2::new(hw, hh);
        let p3 = Point2::new(-hw, hh);

        Profile {
            segments: vec![
                Segment2::Line { start: p0, end: p1 },
                Segment2::Line { start: p1, end: p2 },
                Segment2::Line { start: p2, end: p3 },
                Segment2::Line { start: p3, end: p0 },
            ],
        }
    }

    /// Create a circular profile centered at the origin.
    pub fn circle(radius: f64) -> Self {
        use std::f64::consts::{PI, TAU};
        Profile {
            segments: vec![
                Segment2::Arc {
                    center: Point2::origin(),
                    radius,
                    start_angle: 0.0,
                    end_angle: PI,
                },
                Segment2::Arc {
                    center: Point2::origin(),
                    radius,
                    start_angle: PI,
                    end_angle: TAU,
                },
            ],
        }
    }

    /// Create a profile from a list of vertex positions (polygon).
    pub fn polygon(points: &[Point2]) -> Self {
        assert!(points.len() >= 3, "Polygon needs at least 3 points");
        let mut segments = Vec::with_capacity(points.len());
        for i in 0..points.len() {
            segments.push(Segment2::Line {
                start: points[i],
                end: points[(i + 1) % points.len()],
            });
        }
        Profile { segments }
    }

    /// Evaluate a point on the profile at a normalized parameter t ∈ [0, 1].
    /// t=0 and t=1 are the same point (closed profile).
    pub fn evaluate(&self, t: f64) -> Point2 {
        let n = self.segments.len();
        let t_scaled = t * n as f64;
        let seg_idx = (t_scaled.floor() as usize).min(n - 1);
        let seg_t = t_scaled - seg_idx as f64;

        match &self.segments[seg_idx] {
            Segment2::Line { start, end } => {
                start + (end - start) * seg_t
            }
            Segment2::Arc { center, radius, start_angle, end_angle } => {
                let angle = start_angle + (end_angle - start_angle) * seg_t;
                Point2::new(
                    center.x + radius * angle.cos(),
                    center.y + radius * angle.sin(),
                )
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rect_profile_closed() {
        let p = Profile::rect(4.0, 2.0);
        let start = p.evaluate(0.0);
        let end = p.evaluate(1.0);
        assert!((start - end).norm() < 1e-12, "Profile should be closed");
    }

    #[test]
    fn circle_profile_on_circle() {
        let p = Profile::circle(5.0);
        for i in 0..20 {
            let t = i as f64 / 20.0;
            let pt = p.evaluate(t);
            let r = pt.coords.norm();
            assert!((r - 5.0).abs() < 1e-10, "Circle profile at t={t}: r={r}");
        }
    }

    #[test]
    fn polygon_profile_triangle() {
        let pts = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.5, 1.0),
        ];
        let p = Profile::polygon(&pts);
        assert_eq!(p.segments.len(), 3);
        let start = p.evaluate(0.0);
        assert!((start - pts[0]).norm() < 1e-12);
    }
}
