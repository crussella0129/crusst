//! Arena-based topology store.
//!
//! All topology entities live in the `TopoStore`. Entities reference each other
//! via typed indices (e.g., `VertexId`, `EdgeId`). This avoids Rc/Arc reference
//! cycles in the inherently cyclic topology graph.

use super::types::*;

/// Arena-based storage for all topology entities.
#[derive(Clone, Debug, Default)]
pub struct TopoStore {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub coedges: Vec<CoEdge>,
    pub wires: Vec<Wire>,
    pub faces: Vec<Face>,
    pub shells: Vec<Shell>,
    pub solids: Vec<Solid>,
}

impl TopoStore {
    pub fn new() -> Self {
        Self::default()
    }

    // --- Add entities ---

    pub fn add_vertex(&mut self, vertex: Vertex) -> VertexId {
        let id = VertexId(self.vertices.len());
        self.vertices.push(vertex);
        id
    }

    pub fn add_edge(&mut self, edge: Edge) -> EdgeId {
        let id = EdgeId(self.edges.len());
        self.edges.push(edge);
        id
    }

    pub fn add_coedge(&mut self, coedge: CoEdge) -> CoEdgeId {
        let id = CoEdgeId(self.coedges.len());
        self.coedges.push(coedge);
        id
    }

    pub fn add_wire(&mut self, wire: Wire) -> WireId {
        let id = WireId(self.wires.len());
        self.wires.push(wire);
        id
    }

    pub fn add_face(&mut self, face: Face) -> FaceId {
        let id = FaceId(self.faces.len());
        self.faces.push(face);
        id
    }

    pub fn add_shell(&mut self, shell: Shell) -> ShellId {
        let id = ShellId(self.shells.len());
        self.shells.push(shell);
        id
    }

    pub fn add_solid(&mut self, solid: Solid) -> SolidId {
        let id = SolidId(self.solids.len());
        self.solids.push(solid);
        id
    }

    // --- Get entities ---

    pub fn vertex(&self, id: VertexId) -> &Vertex {
        &self.vertices[id.0]
    }

    pub fn edge(&self, id: EdgeId) -> &Edge {
        &self.edges[id.0]
    }

    pub fn coedge(&self, id: CoEdgeId) -> &CoEdge {
        &self.coedges[id.0]
    }

    pub fn wire(&self, id: WireId) -> &Wire {
        &self.wires[id.0]
    }

    pub fn face(&self, id: FaceId) -> &Face {
        &self.faces[id.0]
    }

    pub fn shell(&self, id: ShellId) -> &Shell {
        &self.shells[id.0]
    }

    pub fn solid(&self, id: SolidId) -> &Solid {
        &self.solids[id.0]
    }

    // --- Get mutable entities ---

    pub fn coedge_mut(&mut self, id: CoEdgeId) -> &mut CoEdge {
        &mut self.coedges[id.0]
    }

    pub fn face_mut(&mut self, id: FaceId) -> &mut Face {
        &mut self.faces[id.0]
    }

    // --- Traversal helpers ---

    /// Iterate over all coedges in a wire loop.
    pub fn wire_coedges(&self, wire_id: WireId) -> Vec<CoEdgeId> {
        let first = self.wire(wire_id).first_coedge;
        let mut result = vec![first];
        let mut current = self.coedge(first).next;
        while current != first {
            result.push(current);
            current = self.coedge(current).next;
        }
        result
    }

    /// Count edges in a wire loop.
    pub fn wire_edge_count(&self, wire_id: WireId) -> usize {
        self.wire_coedges(wire_id).len()
    }

    /// Collect all unique edges referenced by a solid.
    pub fn solid_edges(&self, solid_id: SolidId) -> Vec<EdgeId> {
        let mut edges = Vec::new();
        let shell = self.solid(solid_id);
        for &face_id in &self.shell(shell.outer_shell).faces {
            let face = self.face(face_id);
            for &coedge_id in &self.wire_coedges(face.outer_wire) {
                let edge_id = self.coedge(coedge_id).edge;
                if !edges.contains(&edge_id) {
                    edges.push(edge_id);
                }
            }
            for &inner_wire in &face.inner_wires {
                for &coedge_id in &self.wire_coedges(inner_wire) {
                    let edge_id = self.coedge(coedge_id).edge;
                    if !edges.contains(&edge_id) {
                        edges.push(edge_id);
                    }
                }
            }
        }
        edges
    }

    /// Collect all unique vertices referenced by a solid.
    pub fn solid_vertices(&self, solid_id: SolidId) -> Vec<VertexId> {
        let mut verts = Vec::new();
        for &edge_id in &self.solid_edges(solid_id) {
            let edge = self.edge(edge_id);
            if !verts.contains(&edge.start) {
                verts.push(edge.start);
            }
            if !verts.contains(&edge.end) {
                verts.push(edge.end);
            }
        }
        verts
    }

    /// Count faces in a solid's outer shell.
    pub fn solid_face_count(&self, solid_id: SolidId) -> usize {
        let shell = self.solid(solid_id);
        self.shell(shell.outer_shell).faces.len()
    }
}
