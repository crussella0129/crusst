use crusst::mesh::extract_mesh;
use crusst::shape::{Box3, Sphere};
use nalgebra::Vector3;
use std::collections::HashMap;

fn diagnose(name: &str, mesh: &crusst::mesh::TriangleMesh) {
    let tri_count = mesh.indices.len() / 3;
    let vert_count = mesh.vertices.len();

    let mut edge_counts: HashMap<(u32, u32), u32> = HashMap::new();
    for chunk in mesh.indices.chunks(3) {
        let tri = [chunk[0], chunk[1], chunk[2]];
        for i in 0..3 {
            let (a, b) = (tri[i], tri[(i + 1) % 3]);
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_counts.entry(key).or_insert(0) += 1;
        }
    }

    let boundary = edge_counts.values().filter(|&&c| c == 1).count();
    let manifold = edge_counts.values().filter(|&&c| c == 2).count();
    let nonmanifold = edge_counts.values().filter(|&&c| c > 2).count();
    let total = edge_counts.len();

    let mut pos_map: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for (i, v) in mesh.vertices.iter().enumerate() {
        let key = [(v.x * 1e6) as i64, (v.y * 1e6) as i64, (v.z * 1e6) as i64];
        pos_map.entry(key).or_default().push(i);
    }
    let dup_positions = pos_map.values().filter(|v| v.len() > 1).count();

    println!("=== {} ===", name);
    println!("  Triangles: {}, Vertices: {}", tri_count, vert_count);
    println!("  Edges: {} total", total);
    println!(
        "  Boundary (1 tri):    {} ({:.1}%)",
        boundary,
        boundary as f64 / total as f64 * 100.0
    );
    println!(
        "  Manifold (2 tri):    {} ({:.1}%)",
        manifold,
        manifold as f64 / total as f64 * 100.0
    );
    println!(
        "  Non-manifold (>2):   {} ({:.1}%)",
        nonmanifold,
        nonmanifold as f64 / total as f64 * 100.0
    );
    println!("  Duplicate positions: {}", dup_positions);
    println!(
        "  WATERTIGHT: {}",
        if boundary == 0 && nonmanifold == 0 {
            "YES"
        } else {
            "NO"
        }
    );
    println!();
}

fn main() {
    let sphere = Sphere::new(Vector3::zeros(), 10.0);
    let mesh = extract_mesh(
        &sphere,
        Vector3::from_element(-12.0),
        Vector3::from_element(12.0),
        32,
    );
    diagnose("Sphere (res=32)", &mesh);

    let box3 = Box3::new(Vector3::zeros(), Vector3::new(5.0, 5.0, 5.0));
    let mesh = extract_mesh(
        &box3,
        Vector3::from_element(-7.0),
        Vector3::from_element(7.0),
        32,
    );
    diagnose("Box (res=32)", &mesh);
}
