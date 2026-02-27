use crusst::mesh::extract_mesh;
use crusst::shape::Sphere;
use nalgebra::Vector3;

#[test]
fn sphere_mesh_has_triangles() {
    let sphere = Sphere::new(Vector3::zeros(), 5.0);
    let bbox_min = Vector3::new(-6.0, -6.0, -6.0);
    let bbox_max = Vector3::new(6.0, 6.0, 6.0);
    let mesh = extract_mesh(&sphere, bbox_min, bbox_max, 32);
    assert!(
        mesh.vertices.len() > 10,
        "Mesh should have many vertices, got {}",
        mesh.vertices.len()
    );
    assert!(
        mesh.indices.len() > 10,
        "Mesh should have many triangle indices"
    );
    assert_eq!(
        mesh.indices.len() % 3,
        0,
        "Index count must be multiple of 3"
    );
}

#[test]
fn sphere_mesh_vertices_near_surface() {
    let sphere = Sphere::new(Vector3::zeros(), 5.0);
    let bbox_min = Vector3::new(-6.0, -6.0, -6.0);
    let bbox_max = Vector3::new(6.0, 6.0, 6.0);
    let mesh = extract_mesh(&sphere, bbox_min, bbox_max, 64);
    for v in &mesh.vertices {
        let dist = v.norm();
        assert!(
            (dist - 5.0).abs() < 0.5,
            "Vertex at distance {dist} should be near radius 5.0"
        );
    }
}

#[test]
fn mesh_normals_computed() {
    let sphere = Sphere::new(Vector3::zeros(), 5.0);
    let bbox_min = Vector3::new(-6.0, -6.0, -6.0);
    let bbox_max = Vector3::new(6.0, 6.0, 6.0);
    let mesh = extract_mesh(&sphere, bbox_min, bbox_max, 32);
    assert_eq!(mesh.normals.len(), mesh.vertices.len(), "Normal per vertex");
    for n in &mesh.normals {
        let len = n.norm();
        assert!(
            (len - 1.0).abs() < 0.01,
            "Normals should be unit length, got {len}"
        );
    }
}

#[test]
fn mesh_binary_serialization_roundtrip() {
    let sphere = Sphere::new(Vector3::zeros(), 3.0);
    let mesh = extract_mesh(
        &sphere,
        Vector3::from_element(-4.0),
        Vector3::from_element(4.0),
        16,
    );
    let binary = mesh.to_binary();
    // Should start with vertex count as u32 LE
    let nv = u32::from_le_bytes([binary[0], binary[1], binary[2], binary[3]]);
    assert_eq!(nv as usize, mesh.vertices.len());
}
