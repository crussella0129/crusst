use nalgebra::Vector3;
use crusst::shape::Sphere;
use crusst::mesh::extract_mesh;
use crusst::export::write_stl;
use std::path::PathBuf;

#[test]
#[ignore]
fn write_stl_produces_valid_file() {
    let sphere = Sphere::new(Vector3::zeros(), 5.0);
    let mesh = extract_mesh(
        &sphere,
        Vector3::new(-6.0, -6.0, -6.0),
        Vector3::new(6.0, 6.0, 6.0),
        32,
    );

    let path = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("test_sphere.stl");
    write_stl(&mesh, &path).unwrap();

    // STL binary format: 80 byte header + 4 byte triangle count + 50 bytes per triangle
    let metadata = std::fs::metadata(&path).unwrap();
    let expected_tris = mesh.indices.len() / 3;
    let expected_size = 80 + 4 + expected_tris * 50;
    assert_eq!(metadata.len() as usize, expected_size);

    std::fs::remove_file(&path).ok();
}
