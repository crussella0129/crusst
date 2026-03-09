//! Wavefront OBJ export.
//!
//! Writes vertex positions (`v`), vertex normals (`vn`), and faces (`f`)
//! referencing both position and normal indices.

use crate::types::TriangleMesh;
use std::io::{self, Write};

/// Write a `TriangleMesh` as Wavefront OBJ text to the given writer.
pub fn write_obj<W: Write>(mesh: &TriangleMesh, writer: &mut W) -> io::Result<()> {
    writeln!(writer, "# Crusst B-Rep kernel OBJ export")?;
    writeln!(writer, "# Vertices: {}, Triangles: {}", mesh.vertices.len(), mesh.indices.len() / 3)?;

    // Vertex positions
    for v in &mesh.vertices {
        writeln!(writer, "v {:.6} {:.6} {:.6}", v.x, v.y, v.z)?;
    }

    // Vertex normals
    for n in &mesh.normals {
        writeln!(writer, "vn {:.6} {:.6} {:.6}", n.x, n.y, n.z)?;
    }

    // Faces (OBJ indices are 1-based)
    let n_tris = mesh.indices.len() / 3;
    for t in 0..n_tris {
        let i0 = mesh.indices[t * 3] + 1;
        let i1 = mesh.indices[t * 3 + 1] + 1;
        let i2 = mesh.indices[t * 3 + 2] + 1;
        writeln!(writer, "f {i0}//{i0} {i1}//{i1} {i2}//{i2}")?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::Shape;
    use crate::types::TessSettings;

    #[test]
    fn obj_contains_vertices_and_faces() {
        let mesh = Shape::box3(5.0, 3.0, 8.0).mesh(&TessSettings::default());
        let mut buf = Vec::new();
        write_obj(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        let v_count = text.lines().filter(|l| l.starts_with("v ")).count();
        let vn_count = text.lines().filter(|l| l.starts_with("vn ")).count();
        let f_count = text.lines().filter(|l| l.starts_with("f ")).count();

        assert_eq!(v_count, mesh.vertices.len());
        assert_eq!(vn_count, mesh.normals.len());
        assert_eq!(f_count, mesh.indices.len() / 3);
    }

    #[test]
    fn obj_indices_one_based() {
        let mesh = Shape::box3(1.0, 1.0, 1.0).mesh(&TessSettings::default());
        let mut buf = Vec::new();
        write_obj(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        // No face index should reference 0
        for line in text.lines().filter(|l| l.starts_with("f ")) {
            for part in line.split_whitespace().skip(1) {
                let idx: u32 = part.split("//").next().unwrap().parse().unwrap();
                assert!(idx >= 1, "OBJ indices must be 1-based, got {idx}");
            }
        }
    }
}
