use crusst::builder::Shape;
use nalgebra::Vector3;

#[test]
fn voxelize_sphere_center_is_negative() {
    let s = Shape::sphere(5.0);
    let grid = s.voxelize(1.0); // 1mm voxels
    // Center voxel should be negative (inside)
    let center_idx = grid.index_at(Vector3::zeros());
    assert!(grid.data[center_idx] < 0.0);
}

#[test]
fn voxelize_sphere_far_is_positive() {
    let s = Shape::sphere(5.0);
    let grid = s.voxelize(1.0);
    // Far corner should be positive (outside)
    let far = grid.index_at(Vector3::new(grid.origin.x, grid.origin.y, grid.origin.z));
    assert!(grid.data[far] > 0.0);
}

#[test]
fn voxelize_resolution_matches() {
    let s = Shape::sphere(5.0);
    let grid = s.voxelize(0.5);
    // Bounding box ~14 units wide, 0.5mm voxels -> ~28 cells per axis
    assert!(grid.resolution[0] >= 20);
    assert!(grid.resolution[0] <= 40);
}
