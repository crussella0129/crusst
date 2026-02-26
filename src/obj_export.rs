use crate::mesh::TriangleMesh;
use std::io::Write;
use std::path::Path;

/// Write a Wavefront OBJ file from a triangle mesh.
///
/// Unlike STL, OBJ is an indexed format — vertices are written once and
/// triangles reference them by index. This preserves the shared-vertex
/// topology of the mesh, producing sealed (watertight) solids that viewers
/// and downstream tools can correctly interpret as closed surfaces.
pub fn write_obj(mesh: &TriangleMesh, path: &Path) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;

    writeln!(file, "# Crusst Geometry Kernel — OBJ Export")?;
    writeln!(file, "# Vertices: {}", mesh.vertices.len())?;
    writeln!(file, "# Triangles: {}", mesh.indices.len() / 3)?;
    writeln!(file)?;

    // Vertices (shared — written once, referenced by index)
    for v in &mesh.vertices {
        writeln!(file, "v {:.8} {:.8} {:.8}", v.x, v.y, v.z)?;
    }
    writeln!(file)?;

    // Per-vertex normals
    for n in &mesh.normals {
        writeln!(file, "vn {:.8} {:.8} {:.8}", n.x, n.y, n.z)?;
    }
    writeln!(file)?;

    // Faces (1-indexed in OBJ format), referencing shared vertices
    // Format: f v1//vn1 v2//vn2 v3//vn3
    for chunk in mesh.indices.chunks(3) {
        let (a, b, c) = (chunk[0] + 1, chunk[1] + 1, chunk[2] + 1);
        writeln!(file, "f {}//{} {}//{} {}//{}", a, a, b, b, c, c)?;
    }

    Ok(())
}
