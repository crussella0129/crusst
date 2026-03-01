//! Sweep operation.
//!
//! Creates a solid by sweeping a 2D profile along a 3D path curve.
//! The profile is kept perpendicular to the path tangent (Frenet frame).

use crate::topo::*;

/// Sweep a profile along a path curve to create a solid.
///
/// Returns the new solid (or error if not yet implemented).
pub fn sweep_profile(
    _store: &mut TopoStore,
    _profile_points: &[crate::math::Point2],
    _path: &crate::curve::Curve3,
    _t_start: f64,
    _t_end: f64,
) -> Result<SolidId, SweepError> {
    // Sweep requires:
    //   - Frenet frame computation along the path
    //   - Profile transformation at each step
    //   - Face generation between profile stations
    //   - Cap face generation at endpoints
    Err(SweepError::NotImplemented)
}

/// Errors that can occur during sweep operations.
#[derive(Debug)]
pub enum SweepError {
    /// Sweep operation is not yet implemented.
    NotImplemented,
    /// The path curve has zero length.
    DegeneratePath,
    /// The profile is degenerate (zero area).
    DegenerateProfile,
}

impl std::fmt::Display for SweepError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SweepError::NotImplemented => write!(f, "Sweep operation not yet implemented"),
            SweepError::DegeneratePath => write!(f, "Sweep path has zero length"),
            SweepError::DegenerateProfile => write!(f, "Sweep profile is degenerate"),
        }
    }
}

impl std::error::Error for SweepError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Point2;
    use crate::curve::Curve3;
    use crate::math::Vector3;

    #[test]
    fn sweep_returns_not_implemented() {
        let mut store = TopoStore::new();
        let profile = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.5, 1.0),
        ];
        let path = Curve3::Line {
            origin: crate::math::Point3::origin(),
            dir: Vector3::new(0.0, 0.0, 10.0),
        };
        let result = sweep_profile(&mut store, &profile, &path, 0.0, 1.0);
        assert!(matches!(result, Err(SweepError::NotImplemented)));
    }
}
