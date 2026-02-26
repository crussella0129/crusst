use crate::mesh::TriangleMesh;
use std::io::Write;
use std::path::Path;
use std::collections::HashMap;

/// Write a STEP AP203 file from a triangle mesh.
///
/// This produces a tessellated CLOSED_SHELL representation — every triangle
/// becomes an ADVANCED_FACE with a planar surface. The result can be opened
/// in any STEP-compatible CAD tool (FreeCAD, Fusion 360, SolidWorks) to
/// inspect the quality of the SDF kernel's mesh approximation.
pub fn write_step(mesh: &TriangleMesh, path: &Path) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;

    // STEP header
    writeln!(file, "ISO-10303-21;")?;
    writeln!(file, "HEADER;")?;
    writeln!(file, "FILE_DESCRIPTION(('Crusst SDF Geometry Kernel Output'),'2;1');")?;
    writeln!(file, "FILE_NAME('{}','2026-02-26',('Crusst'),('Crusst Geometry Kernel'),'Crusst 0.1.0','Crusst','');",
        path.file_name().unwrap_or_default().to_string_lossy())?;
    writeln!(file, "FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));")?;
    writeln!(file, "ENDSEC;")?;
    writeln!(file, "DATA;")?;

    let mut id: u64 = 1;
    let mut next_id = || { let current = id; id += 1; current };

    // Deduplicate vertices and build entity IDs
    let mut point_ids: Vec<u64> = Vec::with_capacity(mesh.vertices.len());
    let mut vertex_ids: Vec<u64> = Vec::with_capacity(mesh.vertices.len());

    // Write CARTESIAN_POINTs and VERTEX_POINTs
    for v in &mesh.vertices {
        let pid = next_id();
        writeln!(file, "#{}=CARTESIAN_POINT('',({}E0,{}E0,{}E0));", pid, v.x, v.y, v.z)?;
        point_ids.push(pid);

        let vid = next_id();
        writeln!(file, "#{}=VERTEX_POINT('',#{});", vid, pid)?;
        vertex_ids.push(vid);
    }

    // Direction for Z axis and X axis (used for planes)
    let dir_z_id = next_id();
    writeln!(file, "#{}=DIRECTION('',(0.E0,0.E0,1.E0));", dir_z_id)?;
    let dir_x_id = next_id();
    writeln!(file, "#{}=DIRECTION('',(1.E0,0.E0,0.E0));", dir_x_id)?;

    // Build edges (deduplicated by vertex pair)
    // Key: (min_vertex_idx, max_vertex_idx) -> edge entity ID
    let mut edge_map: HashMap<(usize, usize), u64> = HashMap::new();

    let tri_count = mesh.indices.len() / 3;
    let mut face_ids: Vec<u64> = Vec::with_capacity(tri_count);

    for chunk in mesh.indices.chunks(3) {
        let (ai, bi, ci) = (chunk[0] as usize, chunk[1] as usize, chunk[2] as usize);

        // Create or reuse edges
        let edges = [(ai, bi), (bi, ci), (ci, ai)];
        let mut edge_entity_ids = Vec::with_capacity(3);

        for &(v0, v1) in &edges {
            let key = if v0 < v1 { (v0, v1) } else { (v1, v0) };
            let eid = *edge_map.entry(key).or_insert_with(|| {
                // LINE for edge geometry
                let line_dir_id = next_id();
                let va = &mesh.vertices[key.0];
                let vb = &mesh.vertices[key.1];
                let dx = vb.x - va.x;
                let dy = vb.y - va.y;
                let dz = vb.z - va.z;
                let len = (dx*dx + dy*dy + dz*dz).sqrt();
                if len > 1e-12 {
                    writeln!(file, "#{}=DIRECTION('',({:.15E},{:.15E},{:.15E}));",
                        line_dir_id, dx/len, dy/len, dz/len).unwrap();
                } else {
                    writeln!(file, "#{}=DIRECTION('',(1.E0,0.E0,0.E0));", line_dir_id).unwrap();
                }

                let vec_id = next_id();
                writeln!(file, "#{}=VECTOR('',#{},1.E0);", vec_id, line_dir_id).unwrap();

                let line_id = next_id();
                writeln!(file, "#{}=LINE('',#{},#{});",
                    line_id, point_ids[key.0], vec_id).unwrap();

                let edge_id = next_id();
                writeln!(file, "#{}=EDGE_CURVE('',#{},#{},#{},.T.);",
                    edge_id, vertex_ids[key.0], vertex_ids[key.1], line_id).unwrap();
                edge_id
            });
            edge_entity_ids.push((eid, v0 == key.0)); // true if same orientation
        }

        // ORIENTED_EDGEs for this face
        let mut oe_ids = Vec::with_capacity(3);
        for (eid, same_dir) in &edge_entity_ids {
            let oe_id = next_id();
            let orient = if *same_dir { ".T." } else { ".F." };
            writeln!(file, "#{}=ORIENTED_EDGE('',*,*,#{},{});", oe_id, eid, orient)?;
            oe_ids.push(oe_id);
        }

        // EDGE_LOOP
        let el_id = next_id();
        writeln!(file, "#{}=EDGE_LOOP('',(#{},#{},#{}));",
            el_id, oe_ids[0], oe_ids[1], oe_ids[2])?;

        // FACE_BOUND
        let fb_id = next_id();
        writeln!(file, "#{}=FACE_BOUND('',#{},.T.);", fb_id, el_id)?;

        // Compute face normal for the PLANE
        let va = &mesh.vertices[ai];
        let vb = &mesh.vertices[bi];
        let vc = &mesh.vertices[ci];
        let e1x = vb.x - va.x; let e1y = vb.y - va.y; let e1z = vb.z - va.z;
        let e2x = vc.x - va.x; let e2y = vc.y - va.y; let e2z = vc.z - va.z;
        let nx = e1y * e2z - e1z * e2y;
        let ny = e1z * e2x - e1x * e2z;
        let nz = e1x * e2y - e1y * e2x;
        let nlen = (nx*nx + ny*ny + nz*nz).sqrt();
        let (fnx, fny, fnz) = if nlen > 1e-12 {
            (nx/nlen, ny/nlen, nz/nlen)
        } else {
            (0.0, 0.0, 1.0)
        };

        // Direction for face normal
        let fn_dir_id = next_id();
        writeln!(file, "#{}=DIRECTION('',({:.15E},{:.15E},{:.15E}));",
            fn_dir_id, fnx, fny, fnz)?;

        // Ref direction (pick one perpendicular to normal)
        let (rx, ry, rz) = if fnx.abs() < 0.9 {
            let cx = 0.0 - fnz * 0.0 + fny * 0.0;
            let cy = fnz * 1.0 - fnx * 0.0;
            let cz = fnx * 0.0 - fny * 1.0;
            let clen = (cx*cx + cy*cy + cz*cz).sqrt();
            if clen > 1e-12 { (cy/clen, cz/clen, cx/clen) } else { (0.0, 1.0, 0.0) }
        } else {
            let cx = fny * 1.0 - fnz * 0.0;
            let cy = fnz * 0.0 - fnx * 1.0;
            let cz = fnx * 0.0 - fny * 0.0;
            let clen = (cx*cx + cy*cy + cz*cz).sqrt();
            if clen > 1e-12 { (cx/clen, cy/clen, cz/clen) } else { (1.0, 0.0, 0.0) }
        };

        let ref_dir_id = next_id();
        writeln!(file, "#{}=DIRECTION('',({:.15E},{:.15E},{:.15E}));",
            ref_dir_id, rx, ry, rz)?;

        // AXIS2_PLACEMENT_3D for the plane
        let axis_id = next_id();
        writeln!(file, "#{}=AXIS2_PLACEMENT_3D('',#{},#{},#{});",
            axis_id, point_ids[ai], fn_dir_id, ref_dir_id)?;

        // PLANE
        let plane_id = next_id();
        writeln!(file, "#{}=PLANE('',#{});", plane_id, axis_id)?;

        // ADVANCED_FACE
        let face_id = next_id();
        writeln!(file, "#{}=ADVANCED_FACE('',(#{}),#{},.T.);",
            face_id, fb_id, plane_id)?;

        face_ids.push(face_id);
    }

    // CLOSED_SHELL
    let shell_id = next_id();
    write!(file, "#{}=CLOSED_SHELL('Crusst_Shell',(", shell_id)?;
    for (i, fid) in face_ids.iter().enumerate() {
        if i > 0 { write!(file, ",")?; }
        write!(file, "#{}", fid)?;
    }
    writeln!(file, "));")?;

    // MANIFOLD_SOLID_BREP
    let brep_id = next_id();
    writeln!(file, "#{}=MANIFOLD_SOLID_BREP('Crusst_Body',#{});", brep_id, shell_id)?;

    // Axis placement for shape representation
    let origin_id = next_id();
    writeln!(file, "#{}=CARTESIAN_POINT('',(0.E0,0.E0,0.E0));", origin_id)?;
    let axis_place_id = next_id();
    writeln!(file, "#{}=AXIS2_PLACEMENT_3D('',#{},#{},#{});",
        axis_place_id, origin_id, dir_z_id, dir_x_id)?;

    // Geometric context
    let geo_context_id = next_id();
    writeln!(file,
        "#{}=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{}))GLOBAL_UNIT_ASSIGNED_CONTEXT((#{})));",
        geo_context_id, geo_context_id + 1, geo_context_id + 2)?;

    let uncertainty_id = next_id();
    writeln!(file, "#{}=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-7),#{},'distance accuracy');",
        uncertainty_id, uncertainty_id + 1)?;

    let length_unit_id = next_id();
    writeln!(file, "#{}=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));", length_unit_id)?;

    // Shape representation
    let shape_rep_id = next_id();
    writeln!(file, "#{}=SHAPE_REPRESENTATION('Crusst_Shape',(#{},#{}),#{});",
        shape_rep_id, axis_place_id, brep_id, geo_context_id)?;

    // Product definition
    let product_id = next_id();
    writeln!(file, "#{}=PRODUCT('Crusst_Part','Crusst SDF Part','',(#{}));",
        product_id, product_id + 1)?;

    let product_context_id = next_id();
    writeln!(file, "#{}=PRODUCT_CONTEXT('',#{},'mechanical');",
        product_context_id, product_context_id + 1)?;

    let app_context_id = next_id();
    writeln!(file, "#{}=APPLICATION_CONTEXT('automotive design');", app_context_id)?;

    let pdf_id = next_id();
    writeln!(file, "#{}=PRODUCT_DEFINITION_FORMATION('','',#{});", pdf_id, product_id)?;

    let pdc_id = next_id();
    writeln!(file, "#{}=PRODUCT_DEFINITION_CONTEXT('part definition',#{},'design');",
        pdc_id, app_context_id)?;

    let pd_id = next_id();
    writeln!(file, "#{}=PRODUCT_DEFINITION('design','',#{},#{});",
        pd_id, pdf_id, pdc_id)?;

    let pds_id = next_id();
    writeln!(file, "#{}=PRODUCT_DEFINITION_SHAPE('','',#{});", pds_id, pd_id)?;

    let sdr_id = next_id();
    writeln!(file, "#{}=SHAPE_DEFINITION_REPRESENTATION(#{},#{});",
        sdr_id, pds_id, shape_rep_id)?;

    writeln!(file, "ENDSEC;")?;
    writeln!(file, "END-ISO-10303-21;")?;

    Ok(())
}
