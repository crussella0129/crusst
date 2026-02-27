use nalgebra::Vector3;
use crusst::shape::Sphere;
use crusst::path::{LinePath, HelixPath, SpiralPath};
use crusst::transport::{order0, order1, order2, order3};
use crusst::mesh::extract_mesh;
use crusst::export::write_stl;
use std::path::PathBuf;

fn output_dir() -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("eigenforms");
    std::fs::create_dir_all(&dir).ok();
    dir
}

#[test]
#[ignore]
fn eigenform_order0_sphere() {
    let shape = order0(Sphere::new(Vector3::zeros(), 10.0));
    let mesh = extract_mesh(&shape, Vector3::from_element(-12.0), Vector3::from_element(12.0), 64);
    assert!(mesh.indices.len() / 3 > 100);
    write_stl(&mesh, &output_dir().join("order0_sphere.stl")).unwrap();
}

#[test]
#[ignore]
fn eigenform_order1_cone() {
    let path = LinePath::new(Vector3::zeros(), Vector3::new(0.0, 0.0, 30.0));
    let shape = order1(path, 12.0, |t| 1.0 - t, 128);
    let mesh = extract_mesh(
        &shape,
        Vector3::new(-14.0, -14.0, -2.0),
        Vector3::new(14.0, 14.0, 32.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 100, "Cone should have >100 triangles");
    write_stl(&mesh, &output_dir().join("order1_cone.stl")).unwrap();
}

#[test]
#[ignore]
fn eigenform_order2_helix() {
    let path = HelixPath::new(15.0, 8.0, 3.0);
    let shape = order2(path, 2.5, 256);
    let mesh = extract_mesh(
        &shape,
        Vector3::new(-20.0, -20.0, -4.0),
        Vector3::new(20.0, 20.0, 28.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 200, "Helix should have >200 triangles");
    write_stl(&mesh, &output_dir().join("order2_helix.stl")).unwrap();
}

#[test]
#[ignore]
fn eigenform_order3_horn() {
    let path = SpiralPath::new(
        |t| 20.0 * (1.0 - t * 0.7),
        |t| 40.0 * t,
        2.0,
    );
    let shape = order3(path, 7.0, |t| 1.0 - 0.85 * t, |t| t * std::f64::consts::PI, 256);
    let mesh = extract_mesh(
        &shape,
        Vector3::new(-29.0, -29.0, -9.0),
        Vector3::new(29.0, 29.0, 49.0),
        64,
    );
    assert!(mesh.indices.len() / 3 > 200, "Horn should have >200 triangles");
    write_stl(&mesh, &output_dir().join("order3_horn.stl")).unwrap();
}
