use crusst::dag::SdfNode;
use crusst::dual_contouring::extract_mesh_adaptive;
use crusst::types::{BBox3, MeshSettings, TriangleMesh};
use nalgebra::Vector3;
use std::collections::HashMap;
use std::sync::Arc;

// ---------------------------------------------------------------------------
// Manifold assertion helpers
// ---------------------------------------------------------------------------

/// Quantize a vertex position to a canonical integer key for position-based
/// edge matching. Uses 1e-8 resolution which is well above floating-point
/// noise but fine enough for mesh-scale geometry.
fn quant_pos(v: &Vector3<f64>) -> (i64, i64, i64) {
    const SCALE: f64 = 1e8;
    (
        (v.x * SCALE).round() as i64,
        (v.y * SCALE).round() as i64,
        (v.z * SCALE).round() as i64,
    )
}

/// Position-based edge key: two quantized vertex positions, ordered.
type PosEdge = ((i64, i64, i64), (i64, i64, i64));

fn pos_edge(a: &Vector3<f64>, b: &Vector3<f64>) -> PosEdge {
    let qa = quant_pos(a);
    let qb = quant_pos(b);
    if qa <= qb { (qa, qb) } else { (qb, qa) }
}

/// Assert the mesh is strictly manifold and watertight:
/// 1. Every undirected edge (by position) is shared by exactly 2 triangles.
/// 2. Euler characteristic V - E + F = 2 (genus-0 closed surface).
///
/// Uses position-based edge matching so that `split_sharp_edges` (which
/// duplicates vertices at creases) does not break the check — two triangles
/// that share a geometric edge but have different vertex indices are still
/// recognized as sharing that edge.
fn assert_manifold(mesh: &TriangleMesh, label: &str) {
    assert!(
        !mesh.indices.is_empty(),
        "{}: mesh has no triangles",
        label
    );
    assert_eq!(
        mesh.indices.len() % 3,
        0,
        "{}: index count not divisible by 3",
        label
    );

    // Count edges by geometric position, not by vertex index.
    let mut edge_count: HashMap<PosEdge, usize> = HashMap::new();
    for tri in mesh.indices.chunks(3) {
        for i in 0..3 {
            let a = &mesh.vertices[tri[i] as usize];
            let b = &mesh.vertices[tri[(i + 1) % 3] as usize];
            let key = pos_edge(a, b);
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }

    let non_manifold: Vec<_> = edge_count
        .iter()
        .filter(|&(_, &count)| count != 2)
        .collect();

    assert!(
        non_manifold.is_empty(),
        "{}: {} non-manifold edges out of {} total. Samples: {:?}",
        label,
        non_manifold.len(),
        edge_count.len(),
        non_manifold.iter().take(5).collect::<Vec<_>>(),
    );

    // For Euler characteristic after split_sharp_edges, count unique vertex
    // positions rather than vertex array length (duplicated verts are the
    // same topological vertex).
    let unique_verts: std::collections::HashSet<(i64, i64, i64)> =
        mesh.vertices.iter().map(|v| quant_pos(v)).collect();
    let v = unique_verts.len() as i64;
    let e = edge_count.len() as i64;
    let f = (mesh.indices.len() / 3) as i64;
    let euler = v - e + f;
    assert_eq!(
        euler, 2,
        "{}: Euler characteristic V - E + F = {} (expected 2). V={}, E={}, F={}",
        label, euler, v, e, f
    );
}

/// Assert all triangle normals point outward (agree with SDF gradient).
fn assert_outward_winding(mesh: &TriangleMesh, node: &SdfNode, label: &str) {
    let mut flipped = 0usize;
    let total = mesh.indices.len() / 3;
    for tri in mesh.indices.chunks(3) {
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];
        let face_normal = (v1 - v0).cross(&(v2 - v0));
        if face_normal.norm_squared() < 1e-30 {
            continue;
        }
        let centroid = (v0 + v1 + v2) / 3.0;
        let gradient = node.gradient(centroid);
        if face_normal.dot(&gradient) < 0.0 {
            flipped += 1;
        }
    }
    assert_eq!(
        flipped, 0,
        "{}: {}/{} triangles have wrong winding",
        label, flipped, total
    );
}

fn mesh_node(node: &SdfNode, bbox: &BBox3, max_depth: u8) -> TriangleMesh {
    let settings = MeshSettings {
        max_depth,
        min_depth: (max_depth / 2).max(2),
        edge_tolerance: 1e-6,
    };
    extract_mesh_adaptive(node, bbox, &settings)
}

// ---------------------------------------------------------------------------
// Preserved existing tests
// ---------------------------------------------------------------------------

#[test]
fn dc_sphere_produces_manifold_mesh() {
    let node = SdfNode::Sphere {
        center: Vector3::zeros(),
        radius: 5.0,
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let settings = MeshSettings {
        max_depth: 5,
        min_depth: 3,
        edge_tolerance: 1e-6,
    };
    let mesh = extract_mesh_adaptive(&node, &bbox, &settings);
    assert!(mesh.vertices.len() > 10);
    assert!(mesh.indices.len() > 10);
    assert_eq!(mesh.indices.len() % 3, 0);
}

#[test]
fn dc_box_has_sharp_edges() {
    let node = SdfNode::Box3 {
        center: Vector3::zeros(),
        half_extents: Vector3::new(5.0, 5.0, 5.0),
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let settings = MeshSettings {
        max_depth: 6,
        min_depth: 3,
        edge_tolerance: 1e-6,
    };
    let mesh = extract_mesh_adaptive(&node, &bbox, &settings);

    let corner = Vector3::new(5.0, 5.0, 5.0);
    let min_dist = mesh
        .vertices
        .iter()
        .map(|v| (v - corner).norm())
        .fold(f64::INFINITY, f64::min);
    let cell_size = 14.0 / 64.0;
    assert!(
        min_dist < cell_size * 2.0,
        "Nearest vertex to corner is {:.4}, expected < {:.4}",
        min_dist,
        cell_size * 2.0
    );
}

#[test]
fn dc_higher_depth_produces_more_triangles() {
    let node = SdfNode::Sphere {
        center: Vector3::zeros(),
        radius: 5.0,
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let s4 = MeshSettings {
        max_depth: 4,
        min_depth: 2,
        edge_tolerance: 1e-6,
    };
    let s6 = MeshSettings {
        max_depth: 6,
        min_depth: 2,
        edge_tolerance: 1e-6,
    };
    let m4 = extract_mesh_adaptive(&node, &bbox, &s4);
    let m6 = extract_mesh_adaptive(&node, &bbox, &s6);
    assert!(m6.indices.len() > m4.indices.len());
}

// ---------------------------------------------------------------------------
// Strict manifold tests — Phase 1 additions
// ---------------------------------------------------------------------------

#[test]
fn manifold_sphere() {
    let node = SdfNode::Sphere {
        center: Vector3::zeros(),
        radius: 5.0,
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let mesh = mesh_node(&node, &bbox, 6);
    assert_manifold(&mesh, "sphere");
    assert_outward_winding(&mesh, &node, "sphere");
}

#[test]
fn manifold_box() {
    let node = SdfNode::Box3 {
        center: Vector3::zeros(),
        half_extents: Vector3::new(5.0, 5.0, 5.0),
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let mesh = mesh_node(&node, &bbox, 6);
    assert_manifold(&mesh, "box");
    assert_outward_winding(&mesh, &node, "box");
}

#[test]
fn manifold_cylinder() {
    let node = SdfNode::Cylinder {
        base: Vector3::new(0.0, -5.0, 0.0),
        axis: Vector3::new(0.0, 1.0, 0.0),
        radius: 3.0,
        height: 10.0,
    };
    let bbox = BBox3::new(Vector3::new(-5.0, -7.0, -5.0), Vector3::new(5.0, 7.0, 5.0));
    let mesh = mesh_node(&node, &bbox, 6);
    assert_manifold(&mesh, "cylinder");
    assert_outward_winding(&mesh, &node, "cylinder");
}

#[test]
fn manifold_union_two_spheres() {
    let a = Arc::new(SdfNode::Sphere {
        center: Vector3::new(-2.0, 0.0, 0.0),
        radius: 5.0,
    });
    let b = Arc::new(SdfNode::Sphere {
        center: Vector3::new(2.0, 0.0, 0.0),
        radius: 5.0,
    });
    let node = SdfNode::Union(a, b);
    let bbox = BBox3::new(Vector3::from_element(-9.0), Vector3::from_element(9.0));
    let mesh = mesh_node(&node, &bbox, 6);
    assert_manifold(&mesh, "union_spheres");
    assert_outward_winding(&mesh, &node, "union_spheres");
}

#[test]
fn manifold_smooth_union_two_spheres() {
    let a = Arc::new(SdfNode::Sphere {
        center: Vector3::new(-2.0, 0.0, 0.0),
        radius: 5.0,
    });
    let b = Arc::new(SdfNode::Sphere {
        center: Vector3::new(2.0, 0.0, 0.0),
        radius: 5.0,
    });
    let node = SdfNode::SmoothUnion(a, b, 1.5);
    let bbox = BBox3::new(Vector3::from_element(-10.0), Vector3::from_element(10.0));
    let mesh = mesh_node(&node, &bbox, 6);
    assert_manifold(&mesh, "smooth_union_spheres");
    assert_outward_winding(&mesh, &node, "smooth_union_spheres");
}

#[test]
fn manifold_thin_box() {
    let node = SdfNode::Box3 {
        center: Vector3::zeros(),
        half_extents: Vector3::new(1.0, 1.0, 0.01),
    };
    let bbox = BBox3::new(Vector3::new(-2.0, -2.0, -2.0), Vector3::new(2.0, 2.0, 2.0));
    let mesh = mesh_node(&node, &bbox, 8);
    assert_manifold(&mesh, "thin_box");
    assert_outward_winding(&mesh, &node, "thin_box");
}
