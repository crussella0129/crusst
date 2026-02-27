use nalgebra::Vector3;

use crate::shape::Sdf;

// Re-export so that existing `crate::mesh::TriangleMesh` paths keep working.
pub use crate::types::TriangleMesh;

/// Extract a triangle mesh from a signed distance field.
///
/// # Arguments
/// * `sdf` - The signed distance field to extract the surface from.
/// * `bbox_min` - Minimum corner of the bounding box to sample.
/// * `bbox_max` - Maximum corner of the bounding box to sample.
/// * `resolution` - Number of grid cells along each axis.
///
/// # Returns
/// A `TriangleMesh` with vertices on the zero-isosurface, per-vertex normals,
/// and triangle indices.
///
/// # Panics
/// This is a stub — the adaptive dual contouring mesher will replace this.
pub fn extract_mesh(
    _sdf: &dyn Sdf,
    _bbox_min: Vector3<f64>,
    _bbox_max: Vector3<f64>,
    _resolution: usize,
) -> TriangleMesh {
    todo!("extract_mesh: adaptive dual contouring mesher not yet implemented")
}
