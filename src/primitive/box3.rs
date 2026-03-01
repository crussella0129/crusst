//! Box primitive: 6 planar faces, 12 edges, 8 vertices.
//!
//! Creates an axis-aligned box centered at the origin with half-extents (hx, hy, hz).
//! Euler: V(8) - E(12) + F(6) = 2 ✓

use crate::curve::{Curve2, Curve3};
use crate::math::{Point2, Point3, Vector2, Vector3};
use crate::surface::Surface;
use crate::topo::*;

/// Create an axis-aligned box centered at the origin.
///
/// `hx`, `hy`, `hz` are the half-extents along each axis.
pub fn make_box(store: &mut TopoStore, hx: f64, hy: f64, hz: f64) -> SolidId {
    // 8 vertices
    let v = [
        store.add_vertex(Vertex { point: Point3::new(-hx, -hy, -hz) }), // 0: ---
        store.add_vertex(Vertex { point: Point3::new( hx, -hy, -hz) }), // 1: +--
        store.add_vertex(Vertex { point: Point3::new( hx,  hy, -hz) }), // 2: ++-
        store.add_vertex(Vertex { point: Point3::new(-hx,  hy, -hz) }), // 3: -+-
        store.add_vertex(Vertex { point: Point3::new(-hx, -hy,  hz) }), // 4: --+
        store.add_vertex(Vertex { point: Point3::new( hx, -hy,  hz) }), // 5: +-+
        store.add_vertex(Vertex { point: Point3::new( hx,  hy,  hz) }), // 6: +++
        store.add_vertex(Vertex { point: Point3::new(-hx,  hy,  hz) }), // 7: -++
    ];

    // 12 edges (each edge is a line segment between two vertices)
    // Bottom face edges (z = -hz)
    let e0 = add_line_edge(store, v[0], v[1]);  // bottom: 0→1
    let e1 = add_line_edge(store, v[1], v[2]);  // bottom: 1→2
    let e2 = add_line_edge(store, v[2], v[3]);  // bottom: 2→3
    let e3 = add_line_edge(store, v[3], v[0]);  // bottom: 3→0

    // Top face edges (z = +hz)
    let e4 = add_line_edge(store, v[4], v[5]);  // top: 4→5
    let e5 = add_line_edge(store, v[5], v[6]);  // top: 5→6
    let e6 = add_line_edge(store, v[6], v[7]);  // top: 6→7
    let e7 = add_line_edge(store, v[7], v[4]);  // top: 7→4

    // Vertical edges
    let e8  = add_line_edge(store, v[0], v[4]); // 0→4
    let e9  = add_line_edge(store, v[1], v[5]); // 1→5
    let e10 = add_line_edge(store, v[2], v[6]); // 2→6
    let e11 = add_line_edge(store, v[3], v[7]); // 3→7

    // 6 faces with planar surfaces
    // Each face needs: surface, wire (loop of coedges), outward flag

    // Bottom face (z = -hz): loop v0→v3→v2→v1 (CW from +Z = CCW from -Z outward)
    // Each edge reversed so orientations oppose the side faces.
    let bottom = make_planar_face(
        store,
        Surface::Plane {
            origin: Point3::new(0.0, 0.0, -hz),
            normal: Vector3::new(0.0, 0.0, -1.0),
        },
        &[(e3, false), (e2, false), (e1, false), (e0, false)],
        true,
    );

    // Top face (z = +hz): vertices 4,5,6,7 (CCW when viewed from +Z = outward)
    let top = make_planar_face(
        store,
        Surface::Plane {
            origin: Point3::new(0.0, 0.0, hz),
            normal: Vector3::new(0.0, 0.0, 1.0),
        },
        &[(e4, true), (e5, true), (e6, true), (e7, true)],
        true,
    );

    // Front face (y = -hy): vertices 0,1,5,4
    let front = make_planar_face(
        store,
        Surface::Plane {
            origin: Point3::new(0.0, -hy, 0.0),
            normal: Vector3::new(0.0, -1.0, 0.0),
        },
        &[(e0, true), (e9, true), (e4, false), (e8, false)],
        true,
    );

    // Back face (y = +hy): vertices 2,3,7,6
    let back = make_planar_face(
        store,
        Surface::Plane {
            origin: Point3::new(0.0, hy, 0.0),
            normal: Vector3::new(0.0, 1.0, 0.0),
        },
        &[(e2, true), (e11, true), (e6, false), (e10, false)],
        true,
    );

    // Right face (x = +hx): vertices 1,2,6,5
    let right = make_planar_face(
        store,
        Surface::Plane {
            origin: Point3::new(hx, 0.0, 0.0),
            normal: Vector3::new(1.0, 0.0, 0.0),
        },
        &[(e1, true), (e10, true), (e5, false), (e9, false)],
        true,
    );

    // Left face (x = -hx): vertices 3,0,4,7
    let left = make_planar_face(
        store,
        Surface::Plane {
            origin: Point3::new(-hx, 0.0, 0.0),
            normal: Vector3::new(-1.0, 0.0, 0.0),
        },
        &[(e3, true), (e8, true), (e7, false), (e11, false)],
        true,
    );

    let shell = store.add_shell(Shell {
        faces: vec![bottom, top, front, back, right, left],
    });

    store.add_solid(Solid {
        outer_shell: shell,
        inner_shells: vec![],
    })
}

/// Add a line edge between two vertices.
fn add_line_edge(store: &mut TopoStore, start: VertexId, end: VertexId) -> EdgeId {
    let p0 = store.vertex(start).point;
    let p1 = store.vertex(end).point;
    let dir = p1 - p0;

    store.add_edge(Edge {
        curve: Curve3::Line {
            origin: p0,
            dir,
        },
        t_start: 0.0,
        t_end: 1.0,
        start,
        end,
    })
}

/// Create a planar face with a single wire loop.
///
/// `edge_dirs` is a list of (EdgeId, forward) pairs defining the wire loop.
/// The face is created with a dummy pcurve (unit square mapping).
fn make_planar_face(
    store: &mut TopoStore,
    surface: Surface,
    edge_dirs: &[(EdgeId, bool)],
    outward: bool,
) -> FaceId {
    // Pre-allocate the face so coedges can reference it
    let face_id = store.add_face(Face {
        surface,
        outer_wire: WireId(0), // placeholder, will be updated
        inner_wires: vec![],
        outward,
    });

    let n = edge_dirs.len();
    let mut coedge_ids = Vec::with_capacity(n);

    // Create all coedges first (with placeholder `next`)
    for (i, &(edge_id, forward)) in edge_dirs.iter().enumerate() {
        // Simple linear pcurve around the unit square
        let t = i as f64 / n as f64;
        let t_next = (i + 1) as f64 / n as f64;
        let pcurve = Curve2::Line {
            origin: Point2::new(t, 0.0),
            dir: Vector2::new(t_next - t, 0.0),
        };

        let coedge_id = store.add_coedge(CoEdge {
            edge: edge_id,
            forward,
            pcurve,
            next: CoEdgeId(0), // placeholder
            face: face_id,
        });
        coedge_ids.push(coedge_id);
    }

    // Link coedges into a circular list
    for i in 0..n {
        let next = coedge_ids[(i + 1) % n];
        store.coedge_mut(coedge_ids[i]).next = next;
    }

    // Create wire and update face
    let wire = store.add_wire(Wire {
        first_coedge: coedge_ids[0],
    });
    store.face_mut(face_id).outer_wire = wire;

    face_id
}
