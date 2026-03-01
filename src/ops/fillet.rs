//! Edge fillet operation.
//!
//! Replaces sharp edges with smooth rolling-ball fillet surfaces.
//! The fillet radius determines the size of the blend surface.

use crate::topo::*;

/// Apply a constant-radius fillet to the specified edges of a solid.
///
/// The rolling-ball algorithm:
/// 1. For each edge, compute the two adjacent face normals
/// 2. Compute the fillet arc center offset from the edge
/// 3. Create a new cylindrical/toroidal blend surface
/// 4. Trim the adjacent faces and insert the fillet face
///
/// Returns the modified solid (or the original if no edges were filleted).
pub fn fillet_edges(
    store: &mut TopoStore,
    solid_id: SolidId,
    _edge_ids: &[EdgeId],
    _radius: f64,
) -> Result<SolidId, FilletError> {
    // Fillet is one of the most complex operations in a B-Rep kernel.
    // It requires:
    //   - Surface-surface intersection for trimming
    //   - Rolling ball offset surface computation
    //   - Topology surgery (splitting faces, inserting new faces)
    //   - Handling of multi-edge chains and vertex blends
    //
    // This is a stub — full implementation requires boolean-level
    // surface intersection capabilities.
    let _ = store;
    let _ = solid_id;
    Err(FilletError::NotImplemented)
}

/// Errors that can occur during fillet operations.
#[derive(Debug)]
pub enum FilletError {
    /// Fillet operation is not yet implemented.
    NotImplemented,
    /// The fillet radius is too large for the given edges.
    RadiusTooLarge,
    /// The specified edge is not part of the solid.
    EdgeNotFound,
}

impl std::fmt::Display for FilletError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FilletError::NotImplemented => write!(f, "Fillet operation not yet implemented"),
            FilletError::RadiusTooLarge => write!(f, "Fillet radius too large for geometry"),
            FilletError::EdgeNotFound => write!(f, "Edge not found in solid"),
        }
    }
}

impl std::error::Error for FilletError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitive;

    #[test]
    fn fillet_returns_not_implemented() {
        let mut store = TopoStore::new();
        let solid = primitive::make_box(&mut store, 5.0, 3.0, 8.0);
        let edges = store.solid_edges(solid);
        let result = fillet_edges(&mut store, solid, &edges[..1], 0.5);
        assert!(matches!(result, Err(FilletError::NotImplemented)));
    }
}
