//! Binary STL export.
//!
//! Writes a standard 80-byte header + triangle facets in little-endian binary format.
//! Each triangle stores a facet normal followed by three vertex positions.

use crate::types::TriangleMesh;
use std::io::{self, Write};

/// Write a `TriangleMesh` as a binary STL to the given writer.
pub fn write_stl<W: Write>(mesh: &TriangleMesh, writer: &mut W) -> io::Result<()> {
    // 80-byte header
    let header = b"Binary STL from Crusst B-Rep kernel\0";
    let mut header_buf = [0u8; 80];
    let len = header.len().min(80);
    header_buf[..len].copy_from_slice(&header[..len]);
    writer.write_all(&header_buf)?;

    // Number of triangles
    let n_tris = mesh.indices.len() / 3;
    writer.write_all(&(n_tris as u32).to_le_bytes())?;

    // Each triangle: normal (3×f32) + 3 vertices (3×3×f32) + attribute (u16)
    for tri in 0..n_tris {
        let i0 = mesh.indices[tri * 3] as usize;
        let i1 = mesh.indices[tri * 3 + 1] as usize;
        let i2 = mesh.indices[tri * 3 + 2] as usize;

        let v0 = &mesh.vertices[i0];
        let v1 = &mesh.vertices[i1];
        let v2 = &mesh.vertices[i2];

        // Compute facet normal from vertex positions
        let e1 = v1 - v0;
        let e2 = v2 - v0;
        let n = e1.cross(&e2);
        let len = n.norm();
        let n = if len > 1e-15 { n / len } else { nalgebra::Vector3::new(0.0, 0.0, 1.0) };

        // Write normal
        write_f32(writer, n.x as f32)?;
        write_f32(writer, n.y as f32)?;
        write_f32(writer, n.z as f32)?;

        // Write 3 vertices
        for v in [v0, v1, v2] {
            write_f32(writer, v.x as f32)?;
            write_f32(writer, v.y as f32)?;
            write_f32(writer, v.z as f32)?;
        }

        // Attribute byte count (unused)
        writer.write_all(&0u16.to_le_bytes())?;
    }

    Ok(())
}

fn write_f32<W: Write>(writer: &mut W, val: f32) -> io::Result<()> {
    writer.write_all(&val.to_le_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::Shape;
    use crate::types::TessSettings;

    #[test]
    fn stl_box_valid_size() {
        let mesh = Shape::box3(5.0, 3.0, 8.0).mesh(&TessSettings::default());
        let mut buf = Vec::new();
        write_stl(&mesh, &mut buf).unwrap();

        let n_tris = mesh.indices.len() / 3;
        // Expected: 80 header + 4 count + n_tris * 50
        let expected = 80 + 4 + n_tris * 50;
        assert_eq!(buf.len(), expected, "STL file size mismatch");
    }

    #[test]
    fn stl_header_correct() {
        let mesh = Shape::box3(1.0, 1.0, 1.0).mesh(&TessSettings::default());
        let mut buf = Vec::new();
        write_stl(&mesh, &mut buf).unwrap();

        assert!(buf.len() >= 84);
        assert!(buf[..35].starts_with(b"Binary STL from Crusst"));
    }
}
