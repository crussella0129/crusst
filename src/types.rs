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
    /// of a shape. This scales chord tolerance and max edge length relative
    /// to the geometry's size, so small shapes get enough triangles and large
    /// shapes don't get too many.
    ///
    /// `diag` is the bounding box diagonal length (or any characteristic size).
    pub fn from_bounding_diagonal(diag: f64) -> Self {
        let diag = diag.max(1e-10); // avoid division by zero
        Self {
            // 0.5% of diagonal — tight enough for smooth curves,
            // loose enough that flat faces stay coarse
            chord_tolerance: diag * 0.005,
            // 10% of diagonal — prevents overly long triangles
            max_edge_length: diag * 0.10,
            // At least 8 subdivisions per parametric direction
            min_subdivisions: 8,
        }
    }
}
