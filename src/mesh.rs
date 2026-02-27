use nalgebra::Vector3;

use crate::shape::Sdf;

// Re-export so that existing `crate::mesh::TriangleMesh` paths keep working.
pub use crate::types::TriangleMesh;

// Re-export the adaptive DC entry point for convenience.
pub use crate::dual_contouring::extract_mesh_adaptive;

/// Extract a triangle mesh from a signed distance field.
///
/// This is a compatibility wrapper that delegates to the dual contouring
/// mesher. It converts the resolution parameter to an octree depth and
/// uses central-difference gradients (since we only have `&dyn Sdf`, not
/// an `SdfNode` with analytical gradients).
///
/// # Arguments
/// * `sdf` - The signed distance field to extract the surface from.
/// * `bbox_min` - Minimum corner of the bounding box to sample.
/// * `bbox_max` - Maximum corner of the bounding box to sample.
/// * `resolution` - Number of grid cells along each axis (converted to octree depth).
///
/// # Returns
/// A `TriangleMesh` with vertices on the zero-isosurface, per-vertex normals,
/// and triangle indices.
pub fn extract_mesh(
    sdf: &dyn Sdf,
    bbox_min: Vector3<f64>,
    bbox_max: Vector3<f64>,
    resolution: usize,
) -> TriangleMesh {
    use crate::types::{BBox3, MeshSettings};

    // Convert resolution to depth: 2^depth ~ resolution
    let depth = (resolution as f64).log2().ceil() as u8;

    let bbox = BBox3::new(bbox_min, bbox_max);
    let settings = MeshSettings {
        max_depth: depth,
        min_depth: depth.saturating_sub(2),
        edge_tolerance: 1e-6,
    };

    crate::dual_contouring::extract_mesh_from_sdf(sdf, &bbox, &settings)
}
