use nalgebra::Vector3;

use isosurface::{
    distance::Signed,
    extractor::IndexedInterleavedNormals,
    math::Vec3 as IsoVec3,
    sampler::Sampler,
    source::{CentralDifference, ScalarSource},
    MarchingCubes,
};

use crate::shape::Sdf;

/// A triangle mesh extracted from a signed distance field.
pub struct TriangleMesh {
    /// Vertex positions.
    pub vertices: Vec<Vector3<f64>>,
    /// Per-vertex normals (unit length).
    pub normals: Vec<Vector3<f64>>,
    /// Triangle indices (every 3 consecutive values form one triangle).
    pub indices: Vec<u32>,
}

impl TriangleMesh {
    /// Serialize the mesh to a compact binary format suitable for WebSocket
    /// transmission to a Three.js frontend.
    ///
    /// Layout:
    /// ```text
    /// [nv: u32 LE]
    /// [vertices: f32*3*nv LE]
    /// [normals: f32*3*nv LE]
    /// [ni: u32 LE]
    /// [indices: u32*ni LE]
    /// ```
    pub fn to_binary(&self) -> Vec<u8> {
        let nv = self.vertices.len() as u32;
        let ni = self.indices.len() as u32;

        // Pre-allocate: 4 (nv) + 12*nv (verts) + 12*nv (normals) + 4 (ni) + 4*ni (indices)
        let capacity = 4 + 12 * nv as usize + 12 * nv as usize + 4 + 4 * ni as usize;
        let mut buf = Vec::with_capacity(capacity);

        // Vertex count
        buf.extend_from_slice(&nv.to_le_bytes());

        // Vertices as f32 triples
        for v in &self.vertices {
            buf.extend_from_slice(&(v.x as f32).to_le_bytes());
            buf.extend_from_slice(&(v.y as f32).to_le_bytes());
            buf.extend_from_slice(&(v.z as f32).to_le_bytes());
        }

        // Normals as f32 triples
        for n in &self.normals {
            buf.extend_from_slice(&(n.x as f32).to_le_bytes());
            buf.extend_from_slice(&(n.y as f32).to_le_bytes());
            buf.extend_from_slice(&(n.z as f32).to_le_bytes());
        }

        // Index count
        buf.extend_from_slice(&ni.to_le_bytes());

        // Indices
        for &i in &self.indices {
            buf.extend_from_slice(&i.to_le_bytes());
        }

        buf
    }
}

/// Adapter that maps the isosurface crate's [0,1] coordinate space to world
/// space, then evaluates our `Sdf` trait.
struct SdfAdapter<'a> {
    sdf: &'a dyn Sdf,
    bbox_min: Vector3<f64>,
    bbox_size: Vector3<f64>,
}

impl<'a> ScalarSource for SdfAdapter<'a> {
    fn sample_scalar(&self, p: IsoVec3) -> Signed {
        // Map from [0,1] to world space
        let world = Vector3::new(
            self.bbox_min.x + p.x as f64 * self.bbox_size.x,
            self.bbox_min.y + p.y as f64 * self.bbox_size.y,
            self.bbox_min.z + p.z as f64 * self.bbox_size.z,
        );
        Signed(self.sdf.evaluate(world) as f32)
    }
}

/// Extract a triangle mesh from a signed distance field using marching cubes.
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
pub fn extract_mesh(
    sdf: &dyn Sdf,
    bbox_min: Vector3<f64>,
    bbox_max: Vector3<f64>,
    resolution: usize,
) -> TriangleMesh {
    let bbox_size = bbox_max - bbox_min;

    let adapter = SdfAdapter {
        sdf,
        bbox_min,
        bbox_size,
    };

    let hermite = CentralDifference::new(adapter);
    let sampler = Sampler::new(&hermite);

    let mut raw_vertices: Vec<f32> = Vec::new();
    let mut raw_indices: Vec<u32> = Vec::new();

    let mut mc = MarchingCubes::new(resolution);
    {
        let mut extractor =
            IndexedInterleavedNormals::new(&mut raw_vertices, &mut raw_indices, &hermite);
        mc.extract(&sampler, &mut extractor);
    }

    // raw_vertices is interleaved: [x, y, z, nx, ny, nz, x, y, z, nx, ny, nz, ...]
    // Each vertex has 6 floats.
    let num_vertices = raw_vertices.len() / 6;
    let mut vertices = Vec::with_capacity(num_vertices);
    let mut normals = Vec::with_capacity(num_vertices);

    for i in 0..num_vertices {
        let base = i * 6;
        // Positions are in [0,1] space — map back to world space
        let px = bbox_min.x + raw_vertices[base] as f64 * bbox_size.x;
        let py = bbox_min.y + raw_vertices[base + 1] as f64 * bbox_size.y;
        let pz = bbox_min.z + raw_vertices[base + 2] as f64 * bbox_size.z;
        vertices.push(Vector3::new(px, py, pz));

        // The isosurface crate computes normals via central differences in [0,1]
        // coordinate space. Since our SDF maps [0,1] → world space non-uniformly
        // (when the bounding box isn't a cube), the gradient needs to be divided
        // by the bounding box size per axis to get the correct world-space normal.
        let nx = raw_vertices[base + 3] as f64 / bbox_size.x;
        let ny = raw_vertices[base + 4] as f64 / bbox_size.y;
        let nz = raw_vertices[base + 5] as f64 / bbox_size.z;
        let normal = Vector3::new(nx, ny, nz);
        let len = normal.norm();
        if len > 1e-10 {
            normals.push(normal / len);
        } else {
            normals.push(Vector3::new(0.0, 1.0, 0.0)); // fallback
        }
    }

    TriangleMesh {
        vertices,
        normals,
        indices: raw_indices,
    }
}
