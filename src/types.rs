use nalgebra::Vector3;

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

/// Axis-aligned bounding box.
#[derive(Clone, Copy, Debug)]
pub struct BBox3 {
    pub min: Vector3<f64>,
    pub max: Vector3<f64>,
}

impl BBox3 {
    pub fn new(min: Vector3<f64>, max: Vector3<f64>) -> Self { Self { min, max } }
    pub fn center(&self) -> Vector3<f64> { (self.min + self.max) * 0.5 }
    pub fn size(&self) -> Vector3<f64> { self.max - self.min }
    pub fn contains(&self, p: Vector3<f64>) -> bool {
        p.x >= self.min.x && p.x <= self.max.x &&
        p.y >= self.min.y && p.y <= self.max.y &&
        p.z >= self.min.z && p.z <= self.max.z
    }
    /// Split into 8 octants.
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

/// Closed interval [lo, hi] for interval arithmetic.
#[derive(Clone, Copy, Debug)]
pub struct Interval {
    pub lo: f64,
    pub hi: f64,
}

impl Interval {
    pub fn new(lo: f64, hi: f64) -> Self { Self { lo, hi } }
    pub fn entire() -> Self { Self { lo: f64::NEG_INFINITY, hi: f64::INFINITY } }
    pub fn definitely_positive(&self) -> bool { self.lo > 0.0 }
    pub fn definitely_negative(&self) -> bool { self.hi < 0.0 }

    // -- Interval-interval arithmetic --

    pub fn add(self, other: Interval) -> Interval {
        Interval::new(self.lo + other.lo, self.hi + other.hi)
    }
    pub fn sub(self, other: Interval) -> Interval {
        Interval::new(self.lo - other.hi, self.hi - other.lo)
    }
    pub fn mul(self, other: Interval) -> Interval {
        let a = self.lo * other.lo;
        let b = self.lo * other.hi;
        let c = self.hi * other.lo;
        let d = self.hi * other.hi;
        Interval::new(a.min(b).min(c).min(d), a.max(b).max(c).max(d))
    }
    pub fn union(self, other: Interval) -> Interval {
        Interval::new(self.lo.min(other.lo), self.hi.max(other.hi))
    }

    // -- Unary operations --

    pub fn abs(self) -> Interval {
        if self.lo >= 0.0 {
            self
        } else if self.hi <= 0.0 {
            Interval::new(-self.hi, -self.lo)
        } else {
            Interval::new(0.0, self.lo.abs().max(self.hi))
        }
    }
    pub fn sqrt(self) -> Interval {
        Interval::new(self.lo.max(0.0).sqrt(), self.hi.max(0.0).sqrt())
    }
    pub fn neg(self) -> Interval {
        Interval::new(-self.hi, -self.lo)
    }

    // -- Interval min/max --

    pub fn max(self, other: Interval) -> Interval {
        Interval::new(self.lo.max(other.lo), self.hi.max(other.hi))
    }
    pub fn min(self, other: Interval) -> Interval {
        Interval::new(self.lo.min(other.lo), self.hi.min(other.hi))
    }

    // -- Scalar operations --

    pub fn scalar_add(self, s: f64) -> Interval {
        Interval::new(self.lo + s, self.hi + s)
    }
    pub fn scalar_sub(self, s: f64) -> Interval {
        Interval::new(self.lo - s, self.hi - s)
    }
    pub fn scalar_mul(self, s: f64) -> Interval {
        if s >= 0.0 {
            Interval::new(self.lo * s, self.hi * s)
        } else {
            Interval::new(self.hi * s, self.lo * s)
        }
    }

    // -- Predicates --

    pub fn contains_zero(&self) -> bool {
        self.lo <= 0.0 && self.hi >= 0.0
    }
}

/// Settings controlling adaptive mesh extraction.
pub struct MeshSettings {
    pub max_depth: u8,
    pub min_depth: u8,
    pub edge_tolerance: f64,
}

impl Default for MeshSettings {
    fn default() -> Self {
        Self { max_depth: 8, min_depth: 3, edge_tolerance: 1e-6 }
    }
}
