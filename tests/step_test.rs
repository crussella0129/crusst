use crusst::builder::Shape;

#[test]
fn step_export_sphere_contains_spherical_surface() {
    let s = Shape::sphere(10.0);
    let path = "target/tmp/test_sphere.step";
    std::fs::create_dir_all("target/tmp").unwrap();
    s.export_step(path).unwrap();
    let content = std::fs::read_to_string(path).unwrap();
    assert!(content.contains("SPHERICAL_SURFACE"));
    assert!(content.contains("CLOSED_SHELL"));
}

#[test]
fn step_export_box_contains_planes() {
    let s = Shape::box3(5.0, 3.0, 8.0);
    let path = "target/tmp/test_box.step";
    std::fs::create_dir_all("target/tmp").unwrap();
    s.export_step(path).unwrap();
    let content = std::fs::read_to_string(path).unwrap();
    assert!(content.contains("PLANE"));
    assert!(content.contains("CLOSED_SHELL"));
}

#[test]
fn step_export_smooth_union_uses_tessellated() {
    let s = Shape::sphere(5.0).smooth_union(Shape::sphere(5.0).translate(8.0, 0.0, 0.0), 2.0);
    let path = "target/tmp/test_smooth.step";
    std::fs::create_dir_all("target/tmp").unwrap();
    s.export_step(path).unwrap();
    let content = std::fs::read_to_string(path).unwrap();
    // Smooth union cannot be exact BRep -> should fallback to tessellated
    assert!(content.contains("ADVANCED_FACE"));
}
