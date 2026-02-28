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

impl TriangleMesh {
    /// Split vertices at sharp edges so each face gets correct normals.
    ///
    /// At a sharp crease (e.g., the flat cap meeting the conical wall of a
    /// CappedCone), DC produces one vertex shared by both sides with a single
    /// normal. This causes visible faceting ("chewy" edges) because the
    /// renderer interpolates the wrong normal across faces on the opposite
    /// side of the crease.
    ///
    /// This function recomputes normals per-face, then merges vertex normals
    /// only across edges where adjacent face normals agree (angle < threshold).
    /// Vertices at sharp creases are duplicated so each side gets its own normal.
    ///
    /// `angle_threshold_deg`: edges with adjacent face normal angle above this
    /// are treated as sharp creases (typical: 30-45 degrees).
    pub fn split_sharp_edges(&self, angle_threshold_deg: f64) -> TriangleMesh {
        use std::collections::HashMap;

        let cos_threshold = angle_threshold_deg.to_radians().cos();
        let ntri = self.indices.len() / 3;

        // Phase 1: Compute face normals
        let mut face_normals = Vec::with_capacity(ntri);
        for tri in self.indices.chunks(3) {
            let v0 = self.vertices[tri[0] as usize];
            let v1 = self.vertices[tri[1] as usize];
            let v2 = self.vertices[tri[2] as usize];
            let e1 = v1 - v0;
            let e2 = v2 - v0;
            let n = e1.cross(&e2);
            let len = n.norm();
            if len > 1e-15 {
                face_normals.push(n / len);
            } else {
                face_normals.push(Vector3::new(0.0, 1.0, 0.0));
            }
        }

        // Phase 2: Build edge → face adjacency
        // Key: (min_vertex_idx, max_vertex_idx), Value: list of face indices
        let mut edge_faces: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        for (fi, tri) in self.indices.chunks(3).enumerate() {
            let edges = [
                (tri[0].min(tri[1]), tri[0].max(tri[1])),
                (tri[1].min(tri[2]), tri[1].max(tri[2])),
                (tri[0].min(tri[2]), tri[0].max(tri[2])),
            ];
            for edge in &edges {
                edge_faces.entry(*edge).or_default().push(fi);
            }
        }

        // Phase 3: Build smooth groups via flood fill.
        // Two adjacent faces are in the same smooth group if their normals
        // agree (dot product > cos_threshold).
        // Each face gets assigned a group ID per vertex it uses.
        // face_vertex_group[face_idx][local_vert_idx] = group_id
        // group = (original_vertex_idx, smooth_group_id)

        // For each original vertex, track which faces use it
        let mut vertex_faces: HashMap<u32, Vec<(usize, usize)>> = HashMap::new(); // vert -> [(face_idx, local_idx)]
        for (fi, tri) in self.indices.chunks(3).enumerate() {
            for (li, &vi) in tri.iter().enumerate() {
                vertex_faces.entry(vi).or_default().push((fi, li));
            }
        }

        // For each vertex, group its faces into smooth groups
        let mut new_vertices: Vec<Vector3<f64>> = Vec::new();
        let mut new_normals: Vec<Vector3<f64>> = Vec::new();
        let mut new_indices = vec![0u32; self.indices.len()];

        for (&vi, face_refs) in &vertex_faces {
            // Build adjacency graph among faces sharing this vertex
            // Two faces are smooth-connected if they share an edge through vi
            // AND their face normals agree.
            let nf = face_refs.len();
            let mut visited = vec![false; nf];
            let mut groups: Vec<Vec<usize>> = Vec::new(); // groups of local indices into face_refs

            for start in 0..nf {
                if visited[start] {
                    continue;
                }
                let mut group = vec![start];
                visited[start] = true;
                let mut queue = vec![start];

                while let Some(cur) = queue.pop() {
                    let (cur_fi, _) = face_refs[cur];
                    let cur_fn = face_normals[cur_fi];

                    for next in 0..nf {
                        if visited[next] {
                            continue;
                        }
                        let (next_fi, _) = face_refs[next];
                        let next_fn = face_normals[next_fi];

                        // Check if faces share an edge through vi and normals agree
                        if cur_fn.dot(&next_fn) > cos_threshold && faces_share_edge_through_vertex(
                            &self.indices, cur_fi, next_fi, vi,
                        ) {
                            visited[next] = true;
                            group.push(next);
                            queue.push(next);
                        }
                    }
                }
                groups.push(group);
            }

            // For each smooth group, create one new vertex with averaged normal
            for group in &groups {
                let new_vi = new_vertices.len() as u32;
                new_vertices.push(self.vertices[vi as usize]);

                // Average face normals in this group (area-weighted would be better
                // but uniform is good enough for visual quality)
                let mut avg_n = Vector3::zeros();
                for &local_idx in group {
                    let (fi, _) = face_refs[local_idx];
                    avg_n += face_normals[fi];
                }
                let len = avg_n.norm();
                if len > 1e-15 {
                    new_normals.push(avg_n / len);
                } else {
                    new_normals.push(self.normals[vi as usize]);
                }

                // Update indices for all faces in this group
                for &local_idx in group {
                    let (fi, li) = face_refs[local_idx];
                    new_indices[fi * 3 + li] = new_vi;
                }
            }
        }

        TriangleMesh {
            vertices: new_vertices,
            normals: new_normals,
            indices: new_indices,
        }
    }
}

/// Check if two faces share an edge that passes through vertex `vi`.
fn faces_share_edge_through_vertex(indices: &[u32], fi_a: usize, fi_b: usize, vi: u32) -> bool {
    let ta = &indices[fi_a * 3..fi_a * 3 + 3];
    let tb = &indices[fi_b * 3..fi_b * 3 + 3];
    // Both faces use vi. They share an edge through vi if they share
    // another vertex besides vi.
    for &va in ta {
        if va == vi {
            continue;
        }
        for &vb in tb {
            if va == vb {
                return true;
            }
        }
    }
    false
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
    /// Return the 8 corner vertices of this bounding box.
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
    pub fn new(lo: f64, hi: f64) -> Self {
        debug_assert!(lo <= hi, "Interval lo ({lo}) must be <= hi ({hi})");
        Self { lo, hi }
    }
    pub fn entire() -> Self {
        Self {
            lo: f64::NEG_INFINITY,
            hi: f64::INFINITY,
        }
    }
    pub fn definitely_positive(&self) -> bool {
        self.lo > 0.0
    }
    pub fn definitely_negative(&self) -> bool {
        self.hi < 0.0
    }

    // -- Interval-interval arithmetic --
    // Named methods rather than std::ops traits because interval semantics
    // differ from scalar arithmetic (e.g., sub swaps hi/lo of the rhs).
    #[allow(clippy::should_implement_trait)]
    pub fn add(self, other: Interval) -> Interval {
        Interval::new(self.lo + other.lo, self.hi + other.hi)
    }
    #[allow(clippy::should_implement_trait)]
    pub fn sub(self, other: Interval) -> Interval {
        Interval::new(self.lo - other.hi, self.hi - other.lo)
    }
    #[allow(clippy::should_implement_trait)]
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
    #[allow(clippy::should_implement_trait)]
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
#[derive(Clone, Copy)]
pub struct MeshSettings {
    pub max_depth: u8,
    pub min_depth: u8,
    pub edge_tolerance: f64,
}

impl Default for MeshSettings {
    fn default() -> Self {
        Self {
            max_depth: 8,
            min_depth: 3,
            edge_tolerance: 1e-6,
        }
    }
}
