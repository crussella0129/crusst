//! Cone (frustum) primitive.
//!
//! A capped cone with bottom radius r1 and top radius r2.
//! Same topology as cylinder: 4 vertices, 6 edges, 4 faces.
//! V(4) - E(6) + F(4) = 2 ✓

use crate::curve::{Curve2, Curve3};
use crate::math::{Point2, Point3, Vector2, Vector3};
use crate::surface::Surface;
use crate::topo::*;
use std::f64::consts::PI;

/// Create a cone (frustum) centered on the Z axis.
///
/// `r1` is the bottom radius, `r2` is the top radius, `height` is the full height.
/// If `r2 == 0`, this is a pointed cone; otherwise it's a frustum.
pub fn make_cone(store: &mut TopoStore, r1: f64, r2: f64, height: f64) -> SolidId {
    let half_angle = ((r1 - r2) / height).atan();

    // 4 vertices
    let v_b0 = store.add_vertex(Vertex { point: Point3::new(r1, 0.0, 0.0) });
    let v_b1 = store.add_vertex(Vertex { point: Point3::new(-r1, 0.0, 0.0) });
    let v_t0 = store.add_vertex(Vertex { point: Point3::new(r2, 0.0, height) });
    let v_t1 = store.add_vertex(Vertex { point: Point3::new(-r2, 0.0, height) });

    let axis = Vector3::new(0.0, 0.0, 1.0);

    // Bottom circle: 2 half-circle arcs
    let e_b_front = store.add_edge(Edge {
        curve: Curve3::Circle { center: Point3::new(0.0, 0.0, 0.0), axis, radius: r1 },
        t_start: 0.0, t_end: PI,
        start: v_b0, end: v_b1,
    });
    let e_b_back = store.add_edge(Edge {
        curve: Curve3::Circle { center: Point3::new(0.0, 0.0, 0.0), axis, radius: r1 },
        t_start: PI, t_end: 2.0 * PI,
        start: v_b1, end: v_b0,
    });

    // Top circle: 2 half-circle arcs
    let e_t_front = store.add_edge(Edge {
        curve: Curve3::Circle { center: Point3::new(0.0, 0.0, height), axis, radius: r2 },
        t_start: 0.0, t_end: PI,
        start: v_t0, end: v_t1,
    });
    let e_t_back = store.add_edge(Edge {
        curve: Curve3::Circle { center: Point3::new(0.0, 0.0, height), axis, radius: r2 },
        t_start: PI, t_end: 2.0 * PI,
        start: v_t1, end: v_t0,
    });

    // 2 seam lines (slanted)
    let e_seam0 = store.add_edge(Edge {
        curve: Curve3::Line {
            origin: Point3::new(r1, 0.0, 0.0),
            dir: Point3::new(r2, 0.0, height) - Point3::new(r1, 0.0, 0.0),
        },
        t_start: 0.0, t_end: 1.0,
        start: v_b0, end: v_t0,
    });
    let e_seam1 = store.add_edge(Edge {
        curve: Curve3::Line {
            origin: Point3::new(-r1, 0.0, 0.0),
            dir: Point3::new(-r2, 0.0, height) - Point3::new(-r1, 0.0, 0.0),
        },
        t_start: 0.0, t_end: 1.0,
        start: v_b1, end: v_t1,
    });

    // When r1 ≈ r2, this is a cylinder, not a cone.
    // Use Cylinder surface to avoid degenerate half_angle=0.
    let barrel_surface = if (r1 - r2).abs() < 1e-10 {
        Surface::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::new(0.0, 0.0, 1.0),
            radius: r1,
        }
    } else {
        // Compute cone apex: the point where all slant lines converge.
        let (apex_point, cone_axis) = if r1 > r2 {
            let apex_z = r1 * height / (r1 - r2);
            (Point3::new(0.0, 0.0, apex_z), Vector3::new(0.0, 0.0, -1.0))
        } else {
            let apex_z = -(r1 * height / (r2 - r1));
            (Point3::new(0.0, 0.0, apex_z), Vector3::new(0.0, 0.0, 1.0))
        };
        Surface::Cone {
            apex: apex_point,
            axis: cone_axis,
            half_angle: half_angle.abs(),
        }
    };

    // Faces
    let bottom = make_cone_face(
        store,
        Surface::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vector3::new(0.0, 0.0, -1.0),
        },
        &[(e_b_back, false), (e_b_front, false)],
        true,
    );

    let top = make_cone_face(
        store,
        Surface::Plane {
            origin: Point3::new(0.0, 0.0, height),
            normal: Vector3::new(0.0, 0.0, 1.0),
        },
        &[(e_t_front, true), (e_t_back, true)],
        true,
    );

    let barrel_front = make_cone_face(
        store,
        barrel_surface.clone(),
        &[(e_b_front, true), (e_seam1, true), (e_t_front, false), (e_seam0, false)],
        true,
    );

    let barrel_back = make_cone_face(
        store,
        barrel_surface,
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

fn make_cone_face(
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
