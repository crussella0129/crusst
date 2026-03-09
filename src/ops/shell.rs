//! Shell (thin wall) operation.
//!
//! Creates a hollow solid by offsetting all faces inward by a thickness,
//! optionally removing specified faces to create openings.

use crate::topo::*;

/// Hollow out a solid by offsetting faces inward.
///
/// `thickness` is the wall thickness. `open_faces` are faces to remove
/// (creating openings in the shell).
///
/// Returns the modified solid (or error if not yet implemented).
pub fn shell_solid(
    store: &mut TopoStore,
    solid_id: SolidId,
    _thickness: f64,
    _open_faces: &[FaceId],
) -> Result<SolidId, ShellError> {
    // Shell requires:
    //   - Surface offsetting (each face offset by -thickness)
    //   - Topology reconstruction (inner shell creation)
    //   - Intersection resolution at corners
    let _ = store;
    let _ = solid_id;
    Err(ShellError::NotImplemented)
}

/// Errors that can occur during shell operations.
#[derive(Debug)]
pub enum ShellError {
    /// Shell operation is not yet implemented.
    NotImplemented,
    /// The thickness is too large for the geometry.
    ThicknessTooLarge,
}

impl std::fmt::Display for ShellError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ShellError::NotImplemented => write!(f, "Shell operation not yet implemented"),
            ShellError::ThicknessTooLarge => write!(f, "Shell thickness too large"),
        }
    }
}

impl std::error::Error for ShellError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitive;

    #[test]
    fn shell_returns_not_implemented() {
        let mut store = TopoStore::new();
        let solid = primitive::make_box(&mut store, 5.0, 3.0, 8.0);
        let result = shell_solid(&mut store, solid, 0.5, &[]);
        assert!(matches!(result, Err(ShellError::NotImplemented)));
    }
}
