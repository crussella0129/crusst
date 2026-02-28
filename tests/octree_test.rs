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
