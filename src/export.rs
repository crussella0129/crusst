use crate::mesh::TriangleMesh;
use std::io::Write;
use std::path::Path;

/// Write a binary STL file from a triangle mesh.
pub fn write_stl(mesh: &TriangleMesh, path: &Path) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;

    // 80-byte header
    let header = [0u8; 80];
    file.write_all(&header)?;

    // Triangle count
    let tri_count = (mesh.indices.len() / 3) as u32;
    file.write_all(&tri_count.to_le_bytes())?;

    // Triangles
    for chunk in mesh.indices.chunks(3) {
        let (a, b, c) = (chunk[0] as usize, chunk[1] as usize, chunk[2] as usize);
        let va = &mesh.vertices[a];
        let vb = &mesh.vertices[b];
        let vc = &mesh.vertices[c];

        // Face normal
        let e1 = vb - va;
        let e2 = vc - va;
        let normal = e1.cross(&e2).normalize();

        // Normal (3 floats)
        file.write_all(&(normal.x as f32).to_le_bytes())?;
        file.write_all(&(normal.y as f32).to_le_bytes())?;
        file.write_all(&(normal.z as f32).to_le_bytes())?;

        // Vertex A
        file.write_all(&(va.x as f32).to_le_bytes())?;
        file.write_all(&(va.y as f32).to_le_bytes())?;
        file.write_all(&(va.z as f32).to_le_bytes())?;

        // Vertex B
        file.write_all(&(vb.x as f32).to_le_bytes())?;
        file.write_all(&(vb.y as f32).to_le_bytes())?;
        file.write_all(&(vb.z as f32).to_le_bytes())?;

        // Vertex C
        file.write_all(&(vc.x as f32).to_le_bytes())?;
        file.write_all(&(vc.y as f32).to_le_bytes())?;
        file.write_all(&(vc.z as f32).to_le_bytes())?;

        // Attribute byte count (unused)
        file.write_all(&0u16.to_le_bytes())?;
    }

    Ok(())
}
