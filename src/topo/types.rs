//! Topology entity types and typed index handles.

use crate::curve::{Curve2, Curve3};
use crate::math::Point3;
use crate::surface::Surface;

// --- Typed index handles ---
// These are cheap to copy, store, and compare.

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct VertexId(pub usize);

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct EdgeId(pub usize);

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct CoEdgeId(pub usize);

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct WireId(pub usize);

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct FaceId(pub usize);

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ShellId(pub usize);

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SolidId(pub usize);

// --- Topology entities ---

/// A topological vertex at a specific 3D point.
#[derive(Clone, Debug)]
pub struct Vertex {
    pub point: Point3,
}

/// A topological edge — a bounded curve segment between two vertices.
///
/// The geometric carrier is `curve` parameterized from `t_start` to `t_end`.
/// `start` and `end` are the bounding vertices.
#[derive(Clone, Debug)]
pub struct Edge {
    pub curve: Curve3,
    pub t_start: f64,
    pub t_end: f64,
    pub start: VertexId,
    pub end: VertexId,
}

/// An oriented half-edge (coedge). Each Edge has exactly two CoEdges
/// with opposite orientations, one per adjacent Face.
///
/// CoEdges form a linked list around the wire boundary of a face.
/// The `pcurve` is the edge's representation in the face's (u,v) parameter space.
#[derive(Clone, Debug)]
pub struct CoEdge {
    pub edge: EdgeId,
    /// True if the coedge traverses the edge from start→end; false for end→start.
    pub forward: bool,
    /// The edge curve in the face's parameter space (for trimming).
    pub pcurve: Curve2,
    /// Next coedge in the wire loop (forms a circular linked list).
    pub next: CoEdgeId,
    /// The face this coedge belongs to.
    pub face: FaceId,
}

/// A closed loop of coedges forming a boundary of a face.
///
/// The outer wire has coedges oriented counter-clockwise when viewed
/// from outside the solid. Inner wires (holes) go clockwise.
#[derive(Clone, Debug)]
pub struct Wire {
    /// Any coedge in the loop (entry point for traversal).
    pub first_coedge: CoEdgeId,
}

/// A topological face — a bounded region on a surface.
///
/// The face carries one geometric `Surface` and is bounded by one outer wire
/// and zero or more inner wires (holes).
#[derive(Clone, Debug)]
pub struct Face {
    pub surface: Surface,
    /// The outer boundary wire.
    pub outer_wire: WireId,
    /// Inner boundary wires (holes in the face).
    pub inner_wires: Vec<WireId>,
    /// True if the face normal agrees with the surface normal; false if reversed.
    pub outward: bool,
}

/// A connected set of faces forming a closed (or open) surface.
#[derive(Clone, Debug)]
pub struct Shell {
    pub faces: Vec<FaceId>,
}

/// A solid bounded by one outer shell and zero or more inner shells (cavities).
#[derive(Clone, Debug)]
pub struct Solid {
    pub outer_shell: ShellId,
    pub inner_shells: Vec<ShellId>,
}
