//! 3MF export.
//!
//! Writes a 3MF file as uncompressed XML. In a production implementation,
//! this would be wrapped in a ZIP container, but for now we output the
//! core model XML which can be manually packaged.

use crate::types::TriangleMesh;
use std::io::{self, Write};

/// Write a `TriangleMesh` as 3MF model XML to the given writer.
///
/// Note: A complete 3MF file is a ZIP archive containing this XML at
/// `3D/3dmodel.model`. This function writes just the XML payload.
pub fn write_3mf<W: Write>(mesh: &TriangleMesh, writer: &mut W) -> io::Result<()> {
    writeln!(writer, r#"<?xml version="1.0" encoding="UTF-8"?>"#)?;
    writeln!(writer, r#"<model unit="millimeter" xml:lang="en-US""#)?;
    writeln!(writer, r#"  xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">"#)?;
    writeln!(writer, r#"  <resources>"#)?;
    writeln!(writer, r#"    <object id="1" type="model">"#)?;
    writeln!(writer, r#"      <mesh>"#)?;

    // Vertices
    writeln!(writer, r#"        <vertices>"#)?;
    for v in &mesh.vertices {
        writeln!(writer, r#"          <vertex x="{:.6}" y="{:.6}" z="{:.6}" />"#, v.x, v.y, v.z)?;
    }
    writeln!(writer, r#"        </vertices>"#)?;

    // Triangles
    writeln!(writer, r#"        <triangles>"#)?;
    let n_tris = mesh.indices.len() / 3;
    for t in 0..n_tris {
        let v1 = mesh.indices[t * 3];
        let v2 = mesh.indices[t * 3 + 1];
        let v3 = mesh.indices[t * 3 + 2];
        writeln!(writer, r#"          <triangle v1="{v1}" v2="{v2}" v3="{v3}" />"#)?;
    }
    writeln!(writer, r#"        </triangles>"#)?;

    writeln!(writer, r#"      </mesh>"#)?;
    writeln!(writer, r#"    </object>"#)?;
    writeln!(writer, r#"  </resources>"#)?;
    writeln!(writer, r#"  <build>"#)?;
    writeln!(writer, r#"    <item objectid="1" />"#)?;
    writeln!(writer, r#"  </build>"#)?;
    writeln!(writer, r#"</model>"#)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::Shape;
    use crate::types::TessSettings;

    #[test]
    fn threemf_valid_xml() {
        let mesh = Shape::box3(5.0, 3.0, 8.0).mesh(&TessSettings::default());
        let mut buf = Vec::new();
        write_3mf(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        assert!(text.contains("<?xml"));
        assert!(text.contains("<model"));
        assert!(text.contains("<vertices>"));
        assert!(text.contains("<triangles>"));
        assert!(text.contains("</model>"));
    }

    #[test]
    fn threemf_correct_counts() {
        let mesh = Shape::box3(1.0, 1.0, 1.0).mesh(&TessSettings::default());
        let mut buf = Vec::new();
        write_3mf(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        let vert_count = text.lines().filter(|l| l.contains("<vertex ")).count();
        let tri_count = text.lines().filter(|l| l.contains("<triangle ")).count();

        assert_eq!(vert_count, mesh.vertices.len());
        assert_eq!(tri_count, mesh.indices.len() / 3);
    }
}
