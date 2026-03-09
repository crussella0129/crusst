//! Cylinder primitive: 3 faces (barrel + 2 caps), 2 circles + 2 seam edges, 2 vertices.
//!
//! Actually, for clean Euler formula:
//! 2 circle edges (top/bottom rims) + 2 seam lines = 4 edges
//! 2 seam vertices (where seam meets the circles) = 2 vertices per rim = 4 total
//! Wait, a seam line connects a top vertex to a bottom vertex.
//!
//! Simplest valid: 2 seam vertices on each rim = 4 vertices total.
//! 2 half-circle edges per rim = 4 circle arcs + 2 seam lines = 6 edges.
//! 3 faces: barrel, top cap, bottom cap.
//! But barrel face needs to split into 2 for manifold (each seam line is shared).
//!
//! Better: use the same approach as box — split barrel into 2 half-cylinder faces.
//! 4 vertices (2 top, 2 bottom), 6 edges (2 top arcs, 2 bottom arcs, 2 seams), 4 faces (2 barrel, 2 caps)
//! V(4)-E(6)+F(4) = 2 ✓
//!
//! But that's more faces than needed. Even simpler:
//! 2 vertices (1 top seam, 1 bottom seam), 3 edges (top circle, bottom circle, seam line)
//! Each circle edge appears in 2 coedges (cap + barrel)
//! Seam edge appears in 2 coedges (barrel left + barrel right)
//! But barrel is only 1 face → the seam would have 2 coedges on the same face → not standard manifold.
//!
//! For proper manifold topology:
//! Split with 2 seam lines → 4 faces is cleanest.
//! 4 vertices, 8 edges (2 top half-circles, 2 bottom half-circles, 2 seams, 2 more seams... no)
//!
//! Let me just use: 2 seam vertices per rim (4 total), 2 top arcs + 2 bottom arcs + 2 seam lines = 6 edges
//! 2 barrel faces (front/back half), 1 top cap, 1 bottom cap = 4 faces
//! V(4) - E(6) + F(4) = 2 ✓

use crate::curve::{Curve2, Curve3};
use crate::math::{Point2, Point3, Vector2, Vector3};
use crate::surface::Surface;
use crate::topo::*;
use std::f64::consts::PI;

/// Create a cylinder centered on the Z axis.
///
/// `radius` is the cylinder radius, `height` is the full height.
/// The cylinder extends from z=0 to z=height.
pub fn make_cylinder(store: &mut TopoStore, radius: f64, height: f64) -> SolidId {
    // 4 vertices: 2 on bottom circle, 2 on top circle (at seam positions)
    let v_b0 = store.add_vertex(Vertex { point: Point3::new(radius, 0.0, 0.0) });
    let v_b1 = store.add_vertex(Vertex { point: Point3::new(-radius, 0.0, 0.0) });
    let v_t0 = store.add_vertex(Vertex { point: Point3::new(radius, 0.0, height) });
    let v_t1 = store.add_vertex(Vertex { point: Point3::new(-radius, 0.0, height) });

    let axis = Vector3::new(0.0, 0.0, 1.0);

    // Bottom circle: 2 half-circle arcs
    let e_b_front = store.add_edge(Edge {
        curve: Curve3::Circle { center: Point3::new(0.0, 0.0, 0.0), axis, radius },
        t_start: 0.0, t_end: PI,
        start: v_b0, end: v_b1,
    });
    let e_b_back = store.add_edge(Edge {
        curve: Curve3::Circle { center: Point3::new(0.0, 0.0, 0.0), axis, radius },
        t_start: PI, t_end: 2.0 * PI,
        start: v_b1, end: v_b0,
    });

    // Top circle: 2 half-circle arcs
    let e_t_front = store.add_edge(Edge {
        curve: Curve3::Circle { center: Point3::new(0.0, 0.0, height), axis, radius },
        t_start: 0.0, t_end: PI,
        start: v_t0, end: v_t1,
    });
    let e_t_back = store.add_edge(Edge {
        curve: Curve3::Circle { center: Point3::new(0.0, 0.0, height), axis, radius },
        t_start: PI, t_end: 2.0 * PI,
        start: v_t1, end: v_t0,
    });

    // 2 seam lines (vertical)
    let e_seam0 = store.add_edge(Edge {
        curve: Curve3::Line {
            origin: Point3::new(radius, 0.0, 0.0),
            dir: Vector3::new(0.0, 0.0, height),
        },
        t_start: 0.0, t_end: 1.0,
        start: v_b0, end: v_t0,
    });
    let e_seam1 = store.add_edge(Edge {
        curve: Curve3::Line {
            origin: Point3::new(-radius, 0.0, 0.0),
            dir: Vector3::new(0.0, 0.0, height),
        },
        t_start: 0.0, t_end: 1.0,
        start: v_b1, end: v_t1,
    });

    // Bottom cap face: circle edges reversed (inward normal)
    let bottom = make_face_with_wire(
        store,
        Surface::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vector3::new(0.0, 0.0, -1.0),
        },
        &[(e_b_back, false), (e_b_front, false)],
        true,
    );

    // Top cap face: circle edges forward (outward normal)
    let top = make_face_with_wire(
        store,
        Surface::Plane {
            origin: Point3::new(0.0, 0.0, height),
            normal: Vector3::new(0.0, 0.0, 1.0),
        },
        &[(e_t_front, true), (e_t_back, true)],
        true,
    );

    // Front barrel face: bottom_front → seam1 → top_front_rev → seam0_rev
    let barrel_front = make_face_with_wire(
        store,
        Surface::Cylinder {
            origin: Point3::origin(),
            axis: Vector3::new(0.0, 0.0, 1.0),
            radius,
        },
        &[(e_b_front, true), (e_seam1, true), (e_t_front, false), (e_seam0, false)],
        true,
    );

    // Back barrel face: bottom_back → seam0 → top_back_rev → seam1_rev
    let barrel_back = make_face_with_wire(
        store,
        Surface::Cylinder {
            origin: Point3::origin(),
            axis: Vector3::new(0.0, 0.0, 1.0),
            radius,
        },
        &[(e_b_back, true), (e_seam0, true), (e_t_back, false), (e_seam1, false)],
        true,
    );

    let shell = store.add_shell(Shell {
        faces: vec![bottom, top, barrel_front, barrel_back],
    });

    store.add_solid(Solid {
        outer_shell: shell,
        inner_shells: vec![],
    })
}

fn make_face_with_wire(
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

    let wire = store.add_wire(Wire {
        first_coedge: coedge_ids[0],
    });
    store.face_mut(face_id).outer_wire = wire;

    face_id
}
