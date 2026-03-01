//! Edge chamfer operation.
//!
//! Replaces sharp edges with flat planar cuts at 45° (symmetric chamfer)
//! or at a specified distance/angle (asymmetric chamfer).

use crate::topo::*;

/// Apply a symmetric chamfer to the specified edges of a solid.
///
/// A chamfer replaces each edge with a flat planar face whose width
/// is determined by `distance` (measured along both adjacent faces).
///
/// Returns the modified solid (or error if not yet implemented).
pub fn chamfer_edges(
    store: &mut TopoStore,
    solid_id: SolidId,
    _edge_ids: &[EdgeId],
    _distance: f64,
) -> Result<SolidId, ChamferError> {
    // Chamfer is simpler than fillet (planar cuts instead of blend surfaces)
    // but still requires face trimming and topology surgery.
    let _ = store;
    let _ = solid_id;
    Err(ChamferError::NotImplemented)
}

/// Errors that can occur during chamfer operations.
#[derive(Debug)]
pub enum ChamferError {
    /// Chamfer operation is not yet implemented.
    NotImplemented,
    /// The chamfer distance is too large for the given edges.
    DistanceTooLarge,
    /// The specified edge is not part of the solid.
    EdgeNotFound,
}

impl std::fmt::Display for ChamferError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ChamferError::NotImplemented => write!(f, "Chamfer operation not yet implemented"),
            ChamferError::DistanceTooLarge => write!(f, "Chamfer distance too large for geometry"),
            ChamferError::EdgeNotFound => write!(f, "Edge not found in solid"),
        }
    }
}

impl std::error::Error for ChamferError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitive;

    #[test]
    fn chamfer_returns_not_implemented() {
        let mut store = TopoStore::new();
        let solid = primitive::make_box(&mut store, 5.0, 3.0, 8.0);
        let edges = store.solid_edges(solid);
        let result = chamfer_edges(&mut store, solid, &edges[..1], 0.5);
        assert!(matches!(result, Err(ChamferError::NotImplemented)));
    }
}
