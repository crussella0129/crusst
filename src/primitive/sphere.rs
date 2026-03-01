//! Sphere primitive.
//!
//! A sphere has 1 spherical face, 2 vertices (poles), 2 edges (meridian seams),
//! and 1 wire (one seam edge traversed forward + backward).
//!
//! Topology: V(2) - E(2) + F(1) ... wait, that's 2-2+1 = 1, not 2.
//! Actually, the standard B-Rep sphere uses:
//! - 2 pole vertices
//! - 3 edges: 2 meridian half-seams + 1 equator (or similar decomposition)
//!
//! Simplest valid decomposition:
//! - 1 face, 1 seam edge (self-loop), 1 vertex, 1 wire
//! - But V(1) - E(1) + F(1) = 1, not 2
//!
//! For Euler V-E+F=2, we need the right decomposition.
//! Standard: 2 hemisphere faces, 1 equator edge, 2 pole vertices, 2 meridian edges
//! V(2) - E(3) + F(2) = 1 — still not 2!
//!
//! The issue is that some edges are seam edges (they appear twice in the same face).
//! For Euler counting on manifolds with seam edges, we need to be careful.
//!
//! Simplest valid closed manifold: split into 4 faces (quadrants) meeting at poles.
//! Actually, the simplest that works for tessellation:
//! - 2 faces (upper/lower hemisphere), 1 equator edge, 2 pole vertices, 2 seam edges
//! - Each hemisphere is bounded by: equator + seam, with seam appearing as 2 coedges
//!
//! For simplicity and correct Euler (accounting for degenerate edges at poles):
//! Use a 6-face decomposition like a cube mapped to sphere — "cubed sphere"
//! V(8) - E(12) + F(6) = 2 ✓

use crate::curve::{Curve2, Curve3};
use crate::math::{Point2, Point3, Vector2};
use crate::surface::Surface;
use crate::topo::*;

/// Create a sphere centered at the origin with the given radius.
///
/// Uses a 2-hemisphere decomposition with pole vertices and seam edge.
/// For tessellation, each hemisphere face is parameterized as (u,v) on the sphere.
pub fn make_sphere(store: &mut TopoStore, radius: f64) -> SolidId {
    // Use a latitude-band decomposition:
    // 2 pole vertices, 2 equatorial vertices, 3 latitude circles (edges),
    // 2 meridian edges as seams
    //
    // Actually, simplest correct approach: single face with periodic seam.
    // But for valid Euler formula with our topology, use a 2-face decomposition:
    //
    // Top cap face: north pole to equator
    // Bottom cap face: south pole to equator
    // 1 equator edge (full circle), 2 seam edges (each a half-meridian)
    // 2 poles + 0 equator vertices (seam edges connect poles to themselves... no)
    //
    // Let's use the simplest that works:
    // 2 pole vertices, 1 equator circle edge (from equator vertex to itself)
    // Actually, let's just use 4 triangular faces like an octahedron topology:

    // Octahedron topology: 6 vertices, 12 edges, 8 faces → V-E+F = 6-12+8 = 2 ✓
    // But that's more complexity than needed.

    // Pragmatic approach: 2 hemisphere faces, with the equator as a shared edge
    // and a seam meridian edge. Each face is a "half-sphere" bounded by the equator
    // and the seam.
    //
    // 2 vertices (north pole, south pole — though equator points exist on the boundary)
    // We need equator vertices where the seam meets the equator.
    //
    // Simplest: 2 vertices (equator-seam points), 3 edges, 2 faces
    // v0 = seam point 1 on equator, v1 = seam point 2 on equator (opposite side)
    // e0 = equator arc from v0 to v1 (front half)
    // e1 = equator arc from v1 to v0 (back half)
    // e2 = seam meridian from v0 over north pole to v1
    // e3 = seam meridian from v1 under south pole to v0
    //
    // Face 0 (front): wire = e0 → e2_rev → e1_rev → e3
    // Actually this gets complicated. Let me use the cube-sphere approach.
    //
    // Cube-sphere: project a cube onto a sphere. 8 vertices, 12 edges, 6 faces.
    // V(8) - E(12) + F(6) = 2 ✓
    // Each face is a spherical quadrilateral.

    let s = radius / 3.0_f64.sqrt(); // cube vertex at distance radius from center

    // 8 cube vertices projected onto sphere
    let corners = [
        Point3::new(-s, -s, -s).coords.normalize() * radius,
        Point3::new( s, -s, -s).coords.normalize() * radius,
        Point3::new( s,  s, -s).coords.normalize() * radius,
        Point3::new(-s,  s, -s).coords.normalize() * radius,
        Point3::new(-s, -s,  s).coords.normalize() * radius,
        Point3::new( s, -s,  s).coords.normalize() * radius,
        Point3::new( s,  s,  s).coords.normalize() * radius,
        Point3::new(-s,  s,  s).coords.normalize() * radius,
    ];

    let v: Vec<VertexId> = corners
        .iter()
        .map(|p| store.add_vertex(Vertex { point: Point3::from(*p) }))
        .collect();

    // 12 edges — great circle arcs between adjacent cube-sphere vertices
    let edge_pairs = [
        (0, 1), (1, 2), (2, 3), (3, 0), // bottom ring
        (4, 5), (5, 6), (6, 7), (7, 4), // top ring
        (0, 4), (1, 5), (2, 6), (3, 7), // vertical
    ];

    let edges: Vec<EdgeId> = edge_pairs
        .iter()
        .map(|&(a, b)| {
            let pa = Point3::from(corners[a]);
            let pb = Point3::from(corners[b]);
            add_arc_edge(store, v[a], v[b], pa, pb, radius)
        })
        .collect();

    // 6 faces — each is a spherical surface, bounded by 4 arcs
    let sphere_surface = Surface::Sphere {
        center: Point3::origin(),
        radius,
    };

    // Face indices using same convention as box: bottom, top, front, back, right, left
    let face_defs: [(bool, [(usize, bool); 4]); 6] = [
        // Bottom (edges 0,1,2,3 all reversed for opposite orientation from sides)
        (false, [(3, false), (2, false), (1, false), (0, false)]),
        // Top (edges 4,5,6,7): 4→5→6→7
        (true, [(4, true), (5, true), (6, true), (7, true)]),
        // Front (0, 9, 4r, 8r)
        (true, [(0, true), (9, true), (4, false), (8, false)]),
        // Back (2, 11, 6r, 10r)
        (true, [(2, true), (11, true), (6, false), (10, false)]),
        // Right (1, 10, 5r, 9r)
        (true, [(1, true), (10, true), (5, false), (9, false)]),
        // Left (3, 8, 7r, 11r)
        (true, [(3, true), (8, true), (7, false), (11, false)]),
    ];

    let mut face_ids = Vec::new();
    for &(outward, ref edge_refs) in &face_defs {
        let face_edges: Vec<(EdgeId, bool)> = edge_refs
            .iter()
            .map(|&(idx, fwd)| (edges[idx], fwd))
            .collect();

        let face_id = make_spherical_face(store, sphere_surface.clone(), &face_edges, outward);
        face_ids.push(face_id);
    }

    let shell = store.add_shell(Shell { faces: face_ids });
    store.add_solid(Solid {
        outer_shell: shell,
        inner_shells: vec![],
    })
}

/// Add a great circle arc edge between two points on a sphere.
fn add_arc_edge(
    store: &mut TopoStore,
    start: VertexId,
    end: VertexId,
    p0: Point3,
    p1: Point3,
    radius: f64,
) -> EdgeId {
    // For a great circle arc, use a Circle curve with appropriate parameterization
    let axis = p0.coords.cross(&p1.coords);
    let axis_len = axis.norm();

    if axis_len < 1e-15 {
        // Degenerate: use a line
        let dir = p1 - p0;
        return store.add_edge(Edge {
            curve: Curve3::Line {
                origin: p0,
                dir,
            },
            t_start: 0.0,
            t_end: 1.0,
            start,
            end,
        });
    }

    let axis_n = axis / axis_len;
    let angle = (p0.coords.dot(&p1.coords) / (radius * radius)).clamp(-1.0, 1.0).acos();

    store.add_edge(Edge {
        curve: Curve3::Circle {
            center: Point3::origin(),
            axis: axis_n,
            radius,
        },
        t_start: 0.0,
        t_end: angle,
        start,
        end,
    })
}

/// Create a face on a sphere surface.
fn make_spherical_face(
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
