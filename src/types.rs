use nalgebra::Vector3;

/// A triangle mesh with per-vertex normals.
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

        let capacity = 4 + 12 * nv as usize + 12 * nv as usize + 4 + 4 * ni as usize;
        let mut buf = Vec::with_capacity(capacity);

        buf.extend_from_slice(&nv.to_le_bytes());

        for v in &self.vertices {
            buf.extend_from_slice(&(v.x as f32).to_le_bytes());
            buf.extend_from_slice(&(v.y as f32).to_le_bytes());
            buf.extend_from_slice(&(v.z as f32).to_le_bytes());
        }

        for n in &self.normals {
            buf.extend_from_slice(&(n.x as f32).to_le_bytes());
            buf.extend_from_slice(&(n.y as f32).to_le_bytes());
            buf.extend_from_slice(&(n.z as f32).to_le_bytes());
        }

        buf.extend_from_slice(&ni.to_le_bytes());

        for &i in &self.indices {
            buf.extend_from_slice(&i.to_le_bytes());
        }

        buf
    }
}

/// Axis-aligned bounding box.
#[derive(Clone, Copy, Debug)]
pub struct BBox3 {
    pub min: Vector3<f64>,
    pub max: Vector3<f64>,
}

impl BBox3 {
    pub fn new(min: Vector3<f64>, max: Vector3<f64>) -> Self {
        Self { min, max }
    }
    pub fn center(&self) -> Vector3<f64> {
        (self.min + self.max) * 0.5
    }
    pub fn size(&self) -> Vector3<f64> {
        self.max - self.min
    }
    pub fn contains(&self, p: Vector3<f64>) -> bool {
        p.x >= self.min.x
            && p.x <= self.max.x
            && p.y >= self.min.y
            && p.y <= self.max.y
            && p.z >= self.min.z
            && p.z <= self.max.z
    }
    pub fn corners(&self) -> [Vector3<f64>; 8] {
        let mn = self.min;
        let mx = self.max;
        [
            Vector3::new(mn.x, mn.y, mn.z),
            Vector3::new(mx.x, mn.y, mn.z),
            Vector3::new(mn.x, mx.y, mn.z),
            Vector3::new(mx.x, mx.y, mn.z),
            Vector3::new(mn.x, mn.y, mx.z),
            Vector3::new(mx.x, mn.y, mx.z),
            Vector3::new(mn.x, mx.y, mx.z),
            Vector3::new(mx.x, mx.y, mx.z),
        ]
    }
    pub fn octants(&self) -> [BBox3; 8] {
        let c = self.center();
        let mn = self.min;
        let mx = self.max;
        [
            BBox3::new(Vector3::new(mn.x, mn.y, mn.z), Vector3::new(c.x, c.y, c.z)),
            BBox3::new(Vector3::new(c.x, mn.y, mn.z), Vector3::new(mx.x, c.y, c.z)),
            BBox3::new(Vector3::new(mn.x, c.y, mn.z), Vector3::new(c.x, mx.y, c.z)),
            BBox3::new(Vector3::new(c.x, c.y, mn.z), Vector3::new(mx.x, mx.y, c.z)),
            BBox3::new(Vector3::new(mn.x, mn.y, c.z), Vector3::new(c.x, c.y, mx.z)),
            BBox3::new(Vector3::new(c.x, mn.y, c.z), Vector3::new(mx.x, c.y, mx.z)),
            BBox3::new(Vector3::new(mn.x, c.y, c.z), Vector3::new(c.x, mx.y, mx.z)),
            BBox3::new(Vector3::new(c.x, c.y, c.z), Vector3::new(mx.x, mx.y, mx.z)),
        ]
    }
}

/// Settings controlling adaptive tessellation of B-Rep faces.
#[derive(Clone, Copy, Debug)]
pub struct TessSettings {
    /// Maximum chord deviation — the maximum allowed distance between
    /// the tessellated surface and the true mathematical surface.
    pub chord_tolerance: f64,
    /// Maximum edge length in the tessellation.
    pub max_edge_length: f64,
    /// Minimum number of subdivisions per face in each parametric direction.
    pub min_subdivisions: u32,
}

impl Default for TessSettings {
    fn default() -> Self {
        Self {
            chord_tolerance: 0.01,
            max_edge_length: 5.0,
            min_subdivisions: 4,
        }
    }
}

impl TessSettings {
    /// Compute tessellation settings automatically from the bounding diagonal
    /// of a shape. This is a fallback when curvature info is not available.
    pub fn from_bounding_diagonal(diag: f64) -> Self {
        Self::auto_tune(diag, f64::INFINITY)
    }

    /// Compute tessellation settings from bounding diagonal and minimum
    /// curvature radius across all faces. This produces smooth results by
    /// scaling chord tolerance to the smallest curved feature, not the
    /// overall bounding box (which unfairly penalizes elongated shapes).
    ///
    /// - `diag`: bounding box diagonal length
    /// - `min_curvature_r`: smallest curvature radius (INFINITY if all planar)
    pub fn auto_tune(diag: f64, min_curvature_r: f64) -> Self {
        let diag = diag.max(1e-10);

        // Chord tolerance drives visual smoothness for curved surfaces.
        // For a circle of radius R with n segments, chord deviation =
        // R(1 - cos(π/n)). At 0.3% of R, a radius-5 cylinder gets:
        //   tol = 0.015, initial 32 segs dev=0.024 → refines once → 64 segs
        //   dev=0.006 < 0.015 ✓ → smooth 128 segments around full barrel.
        let chord_tol = if min_curvature_r.is_finite() {
            min_curvature_r * 0.003
        } else {
            diag * 0.001
        };

        // Max edge length is a loose safety net to prevent excessively long
        // triangles on flat regions. Chord deviation is the real quality
        // driver — keep this loose so it doesn't cause over-refinement.
        let max_edge = diag * 0.10;

        // Min subdivisions: 32 for curved shapes (good starting grid for
        // adaptive refinement), 4 for all-planar shapes.
        let min_sub = if min_curvature_r.is_finite() { 32 } else { 4 };

        Self {
            chord_tolerance: chord_tol,
            max_edge_length: max_edge,
            min_subdivisions: min_sub,
        }
    }
}
