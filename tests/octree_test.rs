use crusst::dag::SdfNode;
use crusst::octree::Octree;
use crusst::types::{BBox3, MeshSettings};
use nalgebra::Vector3;

#[test]
fn octree_sphere_has_leaf_cells() {
    let node = SdfNode::Sphere {
        center: Vector3::zeros(),
        radius: 5.0,
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let settings = MeshSettings {
        max_depth: 4,
        min_depth: 2,
        edge_tolerance: 1e-6,
    };
    let tree = Octree::build(&node, &bbox, &settings);
    assert!(tree.leaf_count() > 0);
}

#[test]
fn octree_empty_region_has_no_surface_cells() {
    let node = SdfNode::Sphere {
        center: Vector3::new(100.0, 0.0, 0.0),
        radius: 1.0,
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let settings = MeshSettings {
        max_depth: 6,
        min_depth: 2,
        edge_tolerance: 1e-6,
    };
    let tree = Octree::build(&node, &bbox, &settings);
    assert_eq!(tree.surface_cell_count(), 0);
}

#[test]
fn octree_deeper_produces_more_cells() {
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
    let t4 = Octree::build(&node, &bbox, &s4);
    let t6 = Octree::build(&node, &bbox, &s6);
    assert!(t6.surface_cell_count() > t4.surface_cell_count());
}

#[test]
fn octree_find_leaf_containing_returns_correct_cell() {
    use crusst::dag::SdfNode;
    use crusst::octree::Octree;
    use crusst::types::{BBox3, MeshSettings};
    use nalgebra::Vector3;

    let node = SdfNode::Sphere {
        center: Vector3::zeros(),
        radius: 5.0,
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let settings = MeshSettings {
        max_depth: 4,
        min_depth: 2,
        edge_tolerance: 1e-6,
    };
    let tree = Octree::build(&node, &bbox, &settings);

    // A point near the surface should land in a deep cell
    let probe = Vector3::new(4.9, 0.0, 0.0);
    let leaf = tree.find_leaf_containing(probe);
    assert!(leaf.is_some());
    let cell = leaf.unwrap();
    assert!(cell.bbox.contains(probe));
    assert!(cell.depth >= 3);

    // A point far outside should land in a shallow cell (pruned early)
    let far = Vector3::new(6.5, 6.5, 6.5);
    let leaf_far = tree.find_leaf_containing(far);
    assert!(leaf_far.is_some());
    assert!(leaf_far.unwrap().bbox.contains(far));
}

#[test]
fn octree_balance_enforces_two_to_one() {
    use crusst::dag::SdfNode;
    use crusst::octree::Octree;
    use crusst::types::{BBox3, MeshSettings};
    use nalgebra::Vector3;

    let node = SdfNode::Sphere {
        center: Vector3::zeros(),
        radius: 5.0,
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let settings = MeshSettings {
        max_depth: 6,
        min_depth: 2,
        edge_tolerance: 1e-6,
    };
    let mut tree = Octree::build(&node, &bbox, &settings);
    tree.balance(&node);

    // Verify: for every leaf cell, all 6 face-neighbors differ by at most 1 depth level
    let leaves = tree.collect_all_leaves();
    for cell in &leaves {
        let half = cell.bbox.size() * 0.5;
        let center = cell.bbox.center();
        // Step just past the face to land in the immediate neighbor
        let eps = half.x * 1e-3;
        let dx = half.x + eps;
        let dy = half.y + eps;
        let dz = half.z + eps;
        let probes = [
            center + Vector3::new(dx, 0.0, 0.0),
            center - Vector3::new(dx, 0.0, 0.0),
            center + Vector3::new(0.0, dy, 0.0),
            center - Vector3::new(0.0, dy, 0.0),
            center + Vector3::new(0.0, 0.0, dz),
            center - Vector3::new(0.0, 0.0, dz),
        ];
        for probe in &probes {
            if let Some(neighbor) = tree.find_leaf_containing(*probe) {
                let depth_diff = (cell.depth as i16 - neighbor.depth as i16).abs();
                assert!(
                    depth_diff <= 1,
                    "Balance violation: cell at depth {} (center {:?}) has neighbor at depth {} (center {:?})",
                    cell.depth, center, neighbor.depth, neighbor.bbox.center(),
                );
            }
        }
    }
}
