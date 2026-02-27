//! Feature ID types for addressing faces and edges on SDF primitives.
//!
//! The feature system provides a way to identify and query specific geometric
//! features (faces, edges) on primitives. This is essential for targeted
//! fillet/chamfer operations that apply blends only to selected edges.

/// Kind of geometric feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FeatureKind {
    Face,
    Edge,
}

/// Information about a face on a primitive.
#[derive(Debug, Clone)]
pub struct FaceInfo {
    /// Face index within the primitive.
    pub index: usize,
    /// Human-readable label (e.g., "+X", "-Y", "side").
    pub label: String,
    /// Outward normal (if planar; None for curved).
    pub normal: Option<nalgebra::Vector3<f64>>,
}

/// Information about an edge on a primitive.
#[derive(Debug, Clone)]
pub struct EdgeInfo {
    /// Edge index within the primitive.
    pub index: usize,
    /// Indices of the two faces that share this edge.
    pub face_a: usize,
    pub face_b: usize,
    /// Human-readable label (e.g., "+X/+Y").
    pub label: String,
}

/// Feature Target -- addresses specific faces or edges on a shape.
///
/// Based on the FT(component, body, kind, indices) system.
#[derive(Debug, Clone)]
pub struct FeatureTarget {
    pub component: usize,
    pub body: usize,
    pub kind: FeatureKind,
    pub indices: Vec<usize>,
}

/// Builder for FeatureTarget.
pub struct FtBuilder {
    component: usize,
    body: usize,
}

/// Create a feature target builder for the given component and body.
pub fn ft(component: usize, body: usize) -> FtBuilder {
    FtBuilder { component, body }
}

impl FtBuilder {
    /// Target specific edges by index.
    pub fn edges(self, indices: &[usize]) -> FeatureTarget {
        FeatureTarget {
            component: self.component,
            body: self.body,
            kind: FeatureKind::Edge,
            indices: indices.to_vec(),
        }
    }

    /// Target specific faces by index.
    pub fn faces(self, indices: &[usize]) -> FeatureTarget {
        FeatureTarget {
            component: self.component,
            body: self.body,
            kind: FeatureKind::Face,
            indices: indices.to_vec(),
        }
    }

    /// Target all edges (empty indices = all).
    pub fn all_edges(self) -> FeatureTarget {
        FeatureTarget {
            component: self.component,
            body: self.body,
            kind: FeatureKind::Edge,
            indices: vec![],
        }
    }

    /// Target all faces (empty indices = all).
    pub fn all_faces(self) -> FeatureTarget {
        FeatureTarget {
            component: self.component,
            body: self.body,
            kind: FeatureKind::Face,
            indices: vec![],
        }
    }
}
