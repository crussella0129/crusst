//! Torus primitive.
//!
//! A torus needs seam edges to make a valid manifold topology.
//! With 2 seam circles (one in each parametric direction):
//! 1 vertex (intersection of seams), 2 edges (u-seam circle, v-seam circle), 1 face
//! V(1) - E(2) + F(1) = 0 ≠ 2 — genus-1 surface has Euler characteristic 0. ✓ for torus!
//!
//! But our validator expects V-E+F=2 (genus 0). For a torus, the correct
//! Euler characteristic is 0. We should handle this properly.
//!
//! Alternative: split the torus into 4 faces (quadrants) to avoid the genus issue
//! in the simple validator. This gives:
//! 4 vertices, 8 edges, 4 faces → V-E+F = 0 (still genus 1)
//!
//! The Euler characteristic for genus-1 is indeed 0, so our validator needs
//! to account for this. For now, let's build the simplest correct topology
//! and update the validator to accept genus-1.

use crate::curve::{Curve2, Curve3};
use crate::math::{Point2, Point3, Vector2, Vector3};
use crate::surface::Surface;
use crate::topo::*;
use std::f64::consts::PI;

/// Create a torus centered at the origin with the given major and minor radii.
///
/// The torus is centered on the Z axis.
/// `major_r` is the distance from the center to the tube center.
/// `minor_r` is the tube radius.
///
/// Uses a 4-face decomposition for cleaner tessellation.
pub fn make_torus(store: &mut TopoStore, major_r: f64, minor_r: f64) -> SolidId {
    let axis = Vector3::new(0.0, 0.0, 1.0);

    // 4 vertices at the intersections of 2 seam circles
    // Seam positions: u = 0 and u = π, v = 0 and v = π
    let v00 = store.add_vertex(Vertex { point: torus_point(major_r, minor_r, 0.0, 0.0) });
    let v10 = store.add_vertex(Vertex { point: torus_point(major_r, minor_r, PI, 0.0) });
    let v01 = store.add_vertex(Vertex { point: torus_point(major_r, minor_r, 0.0, PI) });
    let v11 = store.add_vertex(Vertex { point: torus_point(major_r, minor_r, PI, PI) });

    // 8 edges: 4 in u-direction (half major circles), 4 in v-direction (half minor circles)
    // u-edges at v=0: v00→v10, v10→v00 (wrapping)
    let eu0_a = add_torus_u_edge(store, v00, v10, major_r, minor_r, 0.0, 0.0, PI);
    let eu0_b = add_torus_u_edge(store, v10, v00, major_r, minor_r, 0.0, PI, 2.0 * PI);

    // u-edges at v=π: v01→v11, v11→v01
    let eu1_a = add_torus_u_edge(store, v01, v11, major_r, minor_r, PI, 0.0, PI);
    let eu1_b = add_torus_u_edge(store, v11, v01, major_r, minor_r, PI, PI, 2.0 * PI);

    // v-edges at u=0: v00→v01, v01→v00
    let ev0_a = add_torus_v_edge(store, v00, v01, major_r, minor_r, 0.0, 0.0, PI);
    let ev0_b = add_torus_v_edge(store, v01, v00, major_r, minor_r, 0.0, PI, 2.0 * PI);

    // v-edges at u=π: v10→v11, v11→v10
    let ev1_a = add_torus_v_edge(store, v10, v11, major_r, minor_r, PI, 0.0, PI);
    let ev1_b = add_torus_v_edge(store, v11, v10, major_r, minor_r, PI, PI, 2.0 * PI);

    let torus_surface = Surface::Torus {
        center: Point3::origin(),
        axis,
        major_r,
        minor_r,
    };

    // 4 faces (quadrants of the torus surface)
    // Face 0: u ∈ [0, π], v ∈ [0, π]
    // Wire: eu0_a → ev1_a → eu1_a_rev → ev0_a_rev
    let f0 = make_torus_face(
        store,
        torus_surface.clone(),
        &[(eu0_a, true), (ev1_a, true), (eu1_a, false), (ev0_a, false)],
    );

    // Face 1: u ∈ [π, 2π], v ∈ [0, π]
    // Wire: eu0_b → ev0_a → eu1_b_rev → ev1_a_rev
    let f1 = make_torus_face(
        store,
        torus_surface.clone(),
        &[(eu0_b, true), (ev0_a, true), (eu1_b, false), (ev1_a, false)],
    );

    // Face 2: u ∈ [0, π], v ∈ [π, 2π]
    // Wire: eu1_a → ev1_b → eu0_a_rev... no, need to think about orientation
    // Actually: eu1_a → ev1_b → eu0_a_rev → ev0_b_rev... hmm
    // Let me reconsider. For v ∈ [π, 2π]:
    // Wire: eu1_a (forward, from v01→v11) at bottom of this patch,
    //        ev1_b (forward, from v11→v10) on right side,
    //        eu0_b (reversed) at top,
    //        ev0_b (reversed) on left side
    // Actually no, eu0_b goes from v10→v00 at v=0.
    // Let me trace this more carefully.
    //
    // The patch v ∈ [π, 2π] is "below" v=π and wraps back to v=0.
    // Its boundary:
    //   Bottom (v=π): eu1_a from (0,π)→(π,π) ... hmm, eu1_a goes v01→v11
    //   Right (u=π): ev1_b from v11 (u=π,v=π) → v10 (u=π,v=2π≡0)
    //   Top (v=2π≡0): eu0_a reversed... eu0_a goes v00→v10, reversed: v10→v00
    //                  Actually eu0_b goes v10→v00 forward
    //   Wait, eu0_b goes from v10→v00 at v=0 (which is also v=2π).
    //   Left (u=0): ev0_b from v01(u=0,v=π)→v00(u=0,v=2π≡0), reversed: v00→v01
    //
    // So the wire (CCW): eu1_a(fwd) → ev1_b(fwd) → eu0_b(rev) → ev0_b(rev)
    // Hmm, eu0_b(rev) goes v00→v10, but we need v10→v00. eu0_b forward goes v10→v00.
    // Let me just fix the winding:

    let f2 = make_torus_face(
        store,
        torus_surface.clone(),
        &[(eu1_a, true), (ev1_b, true), (eu0_a, false), (ev0_b, false)],
    );

    // Face 3: u ∈ [π, 2π], v ∈ [π, 2π]
    let f3 = make_torus_face(
        store,
        torus_surface,
        &[(eu1_b, true), (ev0_b, true), (eu0_b, false), (ev1_b, false)],
    );

    let shell = store.add_shell(Shell {
        faces: vec![f0, f1, f2, f3],
    });

    store.add_solid(Solid {
        outer_shell: shell,
        inner_shells: vec![],
    })
}

fn torus_point(major_r: f64, minor_r: f64, u: f64, v: f64) -> Point3 {
    let r = major_r + minor_r * v.cos();
    Point3::new(r * u.cos(), r * u.sin(), minor_r * v.sin())
}

fn add_torus_u_edge(
    store: &mut TopoStore,
    start: VertexId,
    end: VertexId,
    major_r: f64,
    minor_r: f64,
    v: f64,
    u_start: f64,
    u_end: f64,
) -> EdgeId {
    let r = major_r + minor_r * v.cos();
    let z = minor_r * v.sin();
    store.add_edge(Edge {
        curve: Curve3::Circle {
            center: Point3::new(0.0, 0.0, z),
            axis: Vector3::new(0.0, 0.0, 1.0),
            radius: r,
        },
        t_start: u_start,
        t_end: u_end,
        start,
        end,
    })
}

fn add_torus_v_edge(
    store: &mut TopoStore,
    start: VertexId,
    end: VertexId,
    major_r: f64,
    minor_r: f64,
    u: f64,
    v_start: f64,
    v_end: f64,
) -> EdgeId {
    // The v-edge is a circle in the plane containing the axis and the point at angle u
    let tube_center = Point3::new(major_r * u.cos(), major_r * u.sin(), 0.0);
    // The tube circle axis is the tangent to the major circle at u,
    // which is (-sin(u), cos(u), 0) × (0,0,1) ... actually the tube circle lies
    // in the plane spanned by the radial direction and Z.
    // Radial direction at u: (cos(u), sin(u), 0)
    // So the tube circle axis is perpendicular to both radial and Z:
    // axis = radial × Z... no, the tube circle is in the plane of radial and Z.
    // Its axis is (-sin(u), cos(u), 0) (tangent to major circle).
    let tube_axis = Vector3::new(-u.sin(), u.cos(), 0.0);

    store.add_edge(Edge {
        curve: Curve3::Circle {
            center: tube_center,
            axis: tube_axis,
            radius: minor_r,
        },
        t_start: v_start,
        t_end: v_end,
        start,
        end,
    })
}

fn make_torus_face(
    store: &mut TopoStore,
    surface: Surface,
    edge_dirs: &[(EdgeId, bool)],
) -> FaceId {
    let face_id = store.add_face(Face {
        surface,
        outer_wire: WireId(0),
        inner_wires: vec![],
        outward: true,
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
