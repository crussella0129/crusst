//! Loft operation.
//!
//! Creates a solid by blending between two or more 2D profiles.
//! The profiles are connected by ruled or interpolated surfaces.

use crate::topo::*;

/// Loft between two or more profiles to create a solid.
///
/// Each profile is a set of 3D points representing a cross-section.
/// The profiles are connected by surfaces that interpolate between them.
///
/// Returns the new solid (or error if not yet implemented).
pub fn loft_profiles(
    _store: &mut TopoStore,
    _profiles: &[Vec<crate::math::Point3>],
) -> Result<SolidId, LoftError> {
    // Loft requires:
    //   - Profile alignment (matching vertex counts, orientation)
    //   - Surface interpolation between profiles
    //   - Cap face generation at endpoints
    //   - Potentially NURBS surface fitting
    Err(LoftError::NotImplemented)
}

/// Errors that can occur during loft operations.
#[derive(Debug)]
pub enum LoftError {
    /// Loft operation is not yet implemented.
    NotImplemented,
    /// Need at least two profiles.
    TooFewProfiles,
    /// Profiles have incompatible vertex counts.
    IncompatibleProfiles,
}

impl std::fmt::Display for LoftError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LoftError::NotImplemented => write!(f, "Loft operation not yet implemented"),
            LoftError::TooFewProfiles => write!(f, "Need at least two profiles"),
            LoftError::IncompatibleProfiles => write!(f, "Profiles have incompatible vertex counts"),
        }
    }
}

impl std::error::Error for LoftError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Point3;

    #[test]
    fn loft_returns_not_implemented() {
        let mut store = TopoStore::new();
        let profile1 = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
        ];
        let profile2 = vec![
            Point3::new(0.0, 0.0, 5.0),
            Point3::new(2.0, 0.0, 5.0),
            Point3::new(1.0, 2.0, 5.0),
        ];
        let result = loft_profiles(&mut store, &[profile1, profile2]);
        assert!(matches!(result, Err(LoftError::NotImplemented)));
    }
}
