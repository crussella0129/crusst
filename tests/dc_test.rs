use crusst::dag::SdfNode;
use crusst::dual_contouring::extract_mesh_adaptive;
use crusst::types::{BBox3, MeshSettings};
use nalgebra::Vector3;

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
fn dc_sphere_is_watertight() {
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

    // Check: every edge appears exactly twice (manifold)
    use std::collections::HashMap;
    let mut edge_count: HashMap<(u32, u32), usize> = HashMap::new();
    for tri in mesh.indices.chunks(3) {
        for i in 0..3 {
            let a = tri[i];
            let b = tri[(i + 1) % 3];
            let key = (a.min(b), a.max(b));
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    let manifold = edge_count.values().filter(|&&c| c == 2).count();
    let total = edge_count.len();
    let ratio = manifold as f64 / total as f64;
    assert!(
        ratio > 0.95,
        "Manifold ratio: {:.2}% ({}/{})",
        ratio * 100.0,
        manifold,
        total
    );
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
    let cell_size = 14.0 / 64.0; // bbox_size / 2^max_depth
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
