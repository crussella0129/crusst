//! Capsule primitive: a cylinder with hemispherical end caps.
//!
//! Topology is like a sphere (cube-sphere approach) but stretched along Z.
//! 8 vertices, 12 edges, 6 faces → V(8) - E(12) + F(6) = 2 ✓
//!
//! The 4 middle faces are cylindrical-ish (stretched sphere), the top and
//! bottom faces are hemispherical caps.

use crate::curve::{Curve2, Curve3};
use crate::math::{Point2, Point3, Vector2, Vector3};
use crate::surface::Surface;
use crate::topo::*;

/// Create a capsule centered on the Z axis.
///
/// `radius` is the hemisphere/cylinder radius.
/// `height` is the total height (must be >= 2*radius; the cylindrical section
/// has height `height - 2*radius`).
pub fn make_capsule(store: &mut TopoStore, radius: f64, height: f64) -> SolidId {
    let half_cyl = (height - 2.0 * radius).max(0.0) / 2.0;

    // Use cube-sphere topology: 8 vertices at the "corners" of the capsule.
    // Bottom hemisphere corners at z = -half_cyl, top at z = +half_cyl
    let s = radius / 2.0_f64.sqrt(); // 45° on the hemisphere

    let corners = [
        // Bottom hemisphere (z offset = -half_cyl)
        capsule_corner(-s, -s, -half_cyl, radius),
        capsule_corner(s, -s, -half_cyl, radius),
        capsule_corner(s, s, -half_cyl, radius),
        capsule_corner(-s, s, -half_cyl, radius),
        // Top hemisphere (z offset = +half_cyl)
        capsule_corner(-s, -s, half_cyl, radius),
        capsule_corner(s, -s, half_cyl, radius),
        capsule_corner(s, s, half_cyl, radius),
        capsule_corner(-s, s, half_cyl, radius),
    ];

    let v: Vec<VertexId> = corners
        .iter()
        .map(|p| store.add_vertex(Vertex { point: *p }))
        .collect();

    // 12 edges connecting vertices (same topology as box)
    let edge_pairs = [
        (0, 1), (1, 2), (2, 3), (3, 0), // bottom ring
        (4, 5), (5, 6), (6, 7), (7, 4), // top ring
        (0, 4), (1, 5), (2, 6), (3, 7), // vertical
    ];

    let edges: Vec<EdgeId> = edge_pairs
        .iter()
        .map(|&(a, b)| add_line_edge(store, v[a], v[b], corners[a], corners[b]))
        .collect();

    // Use sphere surfaces for the caps, cylinder for the barrel faces
    let bottom_sphere = Surface::Sphere {
        center: Point3::new(0.0, 0.0, -half_cyl),
        radius,
    };
    let top_sphere = Surface::Sphere {
        center: Point3::new(0.0, 0.0, half_cyl),
        radius,
    };
    let barrel = Surface::Cylinder {
        origin: Point3::new(0.0, 0.0, -half_cyl),
        axis: Vector3::new(0.0, 0.0, 1.0),
        radius,
    };

    // Faces (same winding as box/sphere)
    let face_defs: [(Surface, bool, [(usize, bool); 4]); 6] = [
        (bottom_sphere, true, [(3, false), (2, false), (1, false), (0, false)]),
        (top_sphere, true, [(4, true), (5, true), (6, true), (7, true)]),
        (barrel.clone(), true, [(0, true), (9, true), (4, false), (8, false)]),
        (barrel.clone(), true, [(2, true), (11, true), (6, false), (10, false)]),
        (barrel.clone(), true, [(1, true), (10, true), (5, false), (9, false)]),
        (barrel, true, [(3, true), (8, true), (7, false), (11, false)]),
    ];

    let mut face_ids = Vec::new();
    for (surface, outward, ref edge_refs) in face_defs {
        let face_edges: Vec<(EdgeId, bool)> = edge_refs
            .iter()
            .map(|&(idx, fwd)| (edges[idx], fwd))
            .collect();
        let face_id = make_face(store, surface, &face_edges, outward);
        face_ids.push(face_id);
    }

    let shell = store.add_shell(Shell { faces: face_ids });
    store.add_solid(Solid {
        outer_shell: shell,
        inner_shells: vec![],
    })
}

/// Project a cube corner onto the hemisphere of a capsule.
fn capsule_corner(x: f64, y: f64, z_center: f64, radius: f64) -> Point3 {
    // Project (x, y, ±s) onto hemisphere centered at (0, 0, z_center)
    let z_local = if z_center >= 0.0 {
        (radius * radius - x * x - y * y).max(0.0).sqrt()
    } else {
        -(radius * radius - x * x - y * y).max(0.0).sqrt()
    };
    // Actually, we want the corner projected onto the sphere centered at z_center
    let dir = Vector3::new(x, y, z_local).normalize();
    Point3::new(dir.x * radius, dir.y * radius, z_center + dir.z * radius)
}

fn add_line_edge(store: &mut TopoStore, start: VertexId, end: VertexId, p0: Point3, p1: Point3) -> EdgeId {
    store.add_edge(Edge {
        curve: Curve3::Line { origin: p0, dir: p1 - p0 },
        t_start: 0.0,
        t_end: 1.0,
        start,
        end,
    })
}

fn make_face(
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
