//! Wedge primitive: a triangular prism (extruded triangle).
//!
//! 5 faces (2 triangular caps + 3 rectangular sides), 9 edges, 6 vertices.
//! V(6) - E(9) + F(5) = 2 ✓

use crate::curve::{Curve2, Curve3};
use crate::math::{Point2, Point3, Vector2, Vector3};
use crate::surface::Surface;
use crate::topo::*;

/// Create a wedge (triangular prism) along the Z axis.
///
/// The triangular cross-section has vertices at:
/// - `(dx, 0)`, `(-dx, 0)`, `(0, dy)` in the XY plane
///
/// The prism extends from `z=0` to `z=dz`.
pub fn make_wedge(store: &mut TopoStore, dx: f64, dy: f64, dz: f64) -> SolidId {
    // 6 vertices: 3 bottom, 3 top
    let v = [
        store.add_vertex(Vertex { point: Point3::new(dx, 0.0, 0.0) }),   // 0: bottom-right
        store.add_vertex(Vertex { point: Point3::new(-dx, 0.0, 0.0) }),  // 1: bottom-left
        store.add_vertex(Vertex { point: Point3::new(0.0, dy, 0.0) }),   // 2: bottom-apex
        store.add_vertex(Vertex { point: Point3::new(dx, 0.0, dz) }),    // 3: top-right
        store.add_vertex(Vertex { point: Point3::new(-dx, 0.0, dz) }),   // 4: top-left
        store.add_vertex(Vertex { point: Point3::new(0.0, dy, dz) }),    // 5: top-apex
    ];

    // 9 edges
    // Bottom triangle
    let e0 = add_line_edge(store, v[0], v[1]); // 0→1
    let e1 = add_line_edge(store, v[1], v[2]); // 1→2
    let e2 = add_line_edge(store, v[2], v[0]); // 2→0

    // Top triangle
    let e3 = add_line_edge(store, v[3], v[4]); // 3→4
    let e4 = add_line_edge(store, v[4], v[5]); // 4→5
    let e5 = add_line_edge(store, v[5], v[3]); // 5→3

    // Vertical edges
    let e6 = add_line_edge(store, v[0], v[3]); // 0→3
    let e7 = add_line_edge(store, v[1], v[4]); // 1→4
    let e8 = add_line_edge(store, v[2], v[5]); // 2→5

    // Bottom face: loop v0→v2→v1 (reversed edges for opposite orientation from sides)
    let bottom = make_planar_face(
        store,
        Surface::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vector3::new(0.0, 0.0, -1.0),
        },
        &[(e2, false), (e1, false), (e0, false)],
        true,
    );

    // Top face: loop v3→v4→v5
    let top = make_planar_face(
        store,
        Surface::Plane {
            origin: Point3::new(0.0, 0.0, dz),
            normal: Vector3::new(0.0, 0.0, 1.0),
        },
        &[(e3, true), (e4, true), (e5, true)],
        true,
    );

    // Front face (y=0 side): v0→v1→v4→v3
    let front = make_planar_face(
        store,
        Surface::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vector3::new(0.0, -1.0, 0.0),
        },
        &[(e0, true), (e7, true), (e3, false), (e6, false)],
        true,
    );

    // Left face: v1→v2→v5→v4
    // Normal must point outward (away from solid interior)
    let n_left = compute_face_normal(
        store.vertex(v[1]).point,
        store.vertex(v[5]).point,
        store.vertex(v[2]).point,
    );
    let left = make_planar_face(
        store,
        Surface::Plane {
            origin: store.vertex(v[1]).point,
            normal: n_left,
        },
        &[(e1, true), (e8, true), (e4, false), (e7, false)],
        true,
    );

    // Right face: v2→v0→v3→v5
    // Normal must point outward (away from solid interior)
    let n_right = compute_face_normal(
        store.vertex(v[2]).point,
        store.vertex(v[3]).point,
        store.vertex(v[0]).point,
    );
    let right = make_planar_face(
        store,
        Surface::Plane {
            origin: store.vertex(v[2]).point,
            normal: n_right,
        },
        &[(e2, true), (e6, true), (e5, false), (e8, false)],
        true,
    );

    let shell = store.add_shell(Shell {
        faces: vec![bottom, top, front, left, right],
    });

    store.add_solid(Solid {
        outer_shell: shell,
        inner_shells: vec![],
    })
}

fn compute_face_normal(p0: Point3, p1: Point3, p2: Point3) -> Vector3 {
    let e1 = p1 - p0;
    let e2 = p2 - p0;
    let n = e1.cross(&e2);
    let len = n.norm();
    if len > 1e-15 { n / len } else { Vector3::new(0.0, 0.0, 1.0) }
}

fn add_line_edge(store: &mut TopoStore, start: VertexId, end: VertexId) -> EdgeId {
    let p0 = store.vertex(start).point;
    let p1 = store.vertex(end).point;
    store.add_edge(Edge {
        curve: Curve3::Line { origin: p0, dir: p1 - p0 },
        t_start: 0.0,
        t_end: 1.0,
        start,
        end,
    })
}

fn make_planar_face(
    store: &mut TopoStore,
    surface: Surface,
    edge_dirs: &[(EdgeId, bool)],
    outward: bool,
) -> FaceId {
    let face_id = store.add_face(Face {
        surface,
        outer_wire: WireId(0),
        inner_wires: vec![],
        outward,
    });

    let n = edge_dirs.len();
    let mut coedge_ids = Vec::with_capacity(n);

    for (i, &(edge_id, forward)) in edge_dirs.iter().enumerate() {
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
            next: CoEdgeId(0),
            face: face_id,
        });
        coedge_ids.push(coedge_id);
    }

    for i in 0..n {
        store.coedge_mut(coedge_ids[i]).next = coedge_ids[(i + 1) % n];
    }

    let wire = store.add_wire(Wire { first_coedge: coedge_ids[0] });
    store.face_mut(face_id).outer_wire = wire;
    face_id
}
