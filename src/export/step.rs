//! STEP AP203 export.
//!
//! Maps B-Rep topology and geometry directly to STEP entities without tessellation.
//! This produces exact geometry interchange — the core value proposition of B-Rep.
//!
//! Entity mapping:
//! - `Vertex` → `VERTEX_POINT` + `CARTESIAN_POINT`
//! - `Edge` → `EDGE_CURVE` + curve entity
//! - `Face` → `ADVANCED_FACE` + `FACE_OUTER_BOUND` + surface entity
//! - `Solid` → `MANIFOLD_SOLID_BREP` + `CLOSED_SHELL`
//!
//! Surface mapping:
//! - `Plane` → `PLANE` + `AXIS2_PLACEMENT_3D`
//! - `Cylinder` → `CYLINDRICAL_SURFACE`
//! - `Cone` → `CONICAL_SURFACE`
//! - `Sphere` → `SPHERICAL_SURFACE`
//! - `Torus` → `TOROIDAL_SURFACE`
//! - `NurbsSurface` → `B_SPLINE_SURFACE_WITH_KNOTS`

use crate::surface::Surface;
use crate::topo::*;
use std::io::{self, Write};

/// Write a solid as STEP AP203 to the given writer.
pub fn write_step<W: Write>(store: &TopoStore, solid_id: SolidId, writer: &mut W) -> io::Result<()> {
    let mut eid = EntityCounter::new();

    // Header
    writeln!(writer, "ISO-10303-21;")?;
    writeln!(writer, "HEADER;")?;
    writeln!(writer, "FILE_DESCRIPTION(('Crusst B-Rep export'),'2;1');")?;
    writeln!(writer, "FILE_NAME('shape.stp','2026-02-28',(''),(''),'',")?;
    writeln!(writer, "  'Crusst B-Rep kernel','');")?;
    writeln!(writer, "FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));")?;
    writeln!(writer, "ENDSEC;")?;
    writeln!(writer, "DATA;")?;

    // Collect all entities
    let solid = store.solid(solid_id);
    let shell = store.shell(solid.outer_shell);

    // Build vertex entities
    let vert_ids: Vec<VertexId> = store.solid_vertices(solid_id);
    let mut vert_entity_map: Vec<(VertexId, u64, u64)> = Vec::new(); // (vid, point_eid, vertex_eid)

    for &vid in &vert_ids {
        let pt = &store.vertex(vid).point;
        let pt_eid = eid.next();
        writeln!(writer, "#{pt_eid}=CARTESIAN_POINT('',({:.6},{:.6},{:.6}));", pt.x, pt.y, pt.z)?;
        let vx_eid = eid.next();
        writeln!(writer, "#{vx_eid}=VERTEX_POINT('',#{pt_eid});")?;
        vert_entity_map.push((vid, pt_eid, vx_eid));
    }

    // Build edge entities
    let edge_ids: Vec<EdgeId> = store.solid_edges(solid_id);
    let mut edge_entity_map: Vec<(EdgeId, u64)> = Vec::new(); // (eid, edge_curve_eid)

    for &edge_id in &edge_ids {
        let edge = store.edge(edge_id);

        // Write the curve geometry
        let curve_eid = write_curve_entity(writer, &edge.curve, &mut eid)?;

        // Find vertex entity IDs
        let v_start_eid = vert_entity_map.iter().find(|(vid, _, _)| *vid == edge.start).unwrap().2;
        let v_end_eid = vert_entity_map.iter().find(|(vid, _, _)| *vid == edge.end).unwrap().2;

        let ec_eid = eid.next();
        writeln!(writer, "#{ec_eid}=EDGE_CURVE('',#{v_start_eid},#{v_end_eid},#{curve_eid},.T.);")?;
        edge_entity_map.push((edge_id, ec_eid));
    }

    // Build face entities
    let mut face_entity_ids: Vec<u64> = Vec::new();

    for &face_id in &shell.faces {
        let face = store.face(face_id);

        // Write surface entity
        let surf_eid = write_surface_entity(writer, &face.surface, &mut eid)?;

        // Build edge loop from wire
        let coedge_ids = store.wire_coedges(face.outer_wire);
        let mut oriented_edge_eids = Vec::new();

        for &ce_id in &coedge_ids {
            let ce = store.coedge(ce_id);
            let ec_eid = edge_entity_map.iter().find(|(eid, _)| *eid == ce.edge).unwrap().1;
            let oe_eid = eid.next();
            let orient = if ce.forward { ".T." } else { ".F." };
            writeln!(writer, "#{oe_eid}=ORIENTED_EDGE('',*,*,#{ec_eid},{orient});")?;
            oriented_edge_eids.push(oe_eid);
        }

        let loop_eid = eid.next();
        let oe_refs: Vec<String> = oriented_edge_eids.iter().map(|e| format!("#{e}")).collect();
        writeln!(writer, "#{loop_eid}=EDGE_LOOP('',({oe_list}));", oe_list = oe_refs.join(","))?;

        let bound_eid = eid.next();
        let orient = if face.outward { ".T." } else { ".F." };
        writeln!(writer, "#{bound_eid}=FACE_OUTER_BOUND('',#{loop_eid},{orient});")?;

        let face_eid = eid.next();
        writeln!(writer, "#{face_eid}=ADVANCED_FACE('', (#{bound_eid}),#{surf_eid},.T.);")?;
        face_entity_ids.push(face_eid);
    }

    // Closed shell
    let shell_eid = eid.next();
    let face_refs: Vec<String> = face_entity_ids.iter().map(|e| format!("#{e}")).collect();
    writeln!(writer, "#{shell_eid}=CLOSED_SHELL('',({face_list}));", face_list = face_refs.join(","))?;

    // Manifold solid brep
    let brep_eid = eid.next();
    writeln!(writer, "#{brep_eid}=MANIFOLD_SOLID_BREP('Shape',#{shell_eid});")?;

    // ─── AP203 Product Context (required for STEP readers to find geometry) ───

    // Units context
    let len_unit = eid.next();
    writeln!(writer, "#{len_unit}=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));")?;

    let angle_unit = eid.next();
    writeln!(writer, "#{angle_unit}=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));")?;

    let solid_angle_unit = eid.next();
    writeln!(writer, "#{solid_angle_unit}=(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT());")?;

    // Uncertainty measure for length
    let uncertainty = eid.next();
    writeln!(writer, "#{uncertainty}=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#{len_unit},'distance_accuracy_value','confusion accuracy');")?;

    // Representation context
    let rep_context = eid.next();
    writeln!(writer, "#{rep_context}=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{uncertainty}))GLOBAL_UNIT_ASSIGNED_CONTEXT((#{len_unit},#{angle_unit},#{solid_angle_unit}))REPRESENTATION_CONTEXT('Context3D','3D Context with 1e-7 uncertainty'));")?;

    // ADVANCED_BREP_SHAPE_REPRESENTATION links geometry to context
    let shape_rep = eid.next();
    writeln!(writer, "#{shape_rep}=ADVANCED_BREP_SHAPE_REPRESENTATION('Shape',(#{brep_eid}),#{rep_context});")?;

    // Application context
    let app_ctx = eid.next();
    writeln!(writer, "#{app_ctx}=APPLICATION_CONTEXT('core data for automotive mechanical design processes');")?;

    let app_proto = eid.next();
    writeln!(writer, "#{app_proto}=APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#{app_ctx});")?;

    // Product context
    let prod_ctx = eid.next();
    writeln!(writer, "#{prod_ctx}=PRODUCT_CONTEXT('',#{app_ctx},'mechanical');")?;

    // Product
    let product = eid.next();
    writeln!(writer, "#{product}=PRODUCT('Shape','Shape','',(#{prod_ctx}));")?;

    // Product definition formation
    let pdf = eid.next();
    writeln!(writer, "#{pdf}=PRODUCT_DEFINITION_FORMATION('','',#{product});")?;

    // Product definition context
    let pdc = eid.next();
    writeln!(writer, "#{pdc}=PRODUCT_DEFINITION_CONTEXT('part definition',#{app_ctx},'design');")?;

    // Product definition
    let prod_def = eid.next();
    writeln!(writer, "#{prod_def}=PRODUCT_DEFINITION('design','',#{pdf},#{pdc});")?;

    // Product definition shape
    let pds = eid.next();
    writeln!(writer, "#{pds}=PRODUCT_DEFINITION_SHAPE('','',#{prod_def});")?;

    // Shape definition representation — connects product to geometry
    let sdr = eid.next();
    writeln!(writer, "#{sdr}=SHAPE_DEFINITION_REPRESENTATION(#{pds},#{shape_rep});")?;

    writeln!(writer, "ENDSEC;")?;
    writeln!(writer, "END-ISO-10303-21;")?;

    Ok(())
}

struct EntityCounter(u64);

impl EntityCounter {
    fn new() -> Self { EntityCounter(0) }
    fn next(&mut self) -> u64 {
        self.0 += 1;
        self.0
    }
}

fn write_curve_entity<W: Write>(
    writer: &mut W,
    curve: &crate::curve::Curve3,
    eid: &mut EntityCounter,
) -> io::Result<u64> {
    use crate::curve::Curve3;
    match curve {
        Curve3::Line { origin, dir } => {
            let pt_eid = eid.next();
            writeln!(writer, "#{pt_eid}=CARTESIAN_POINT('',({:.6},{:.6},{:.6}));",
                     origin.x, origin.y, origin.z)?;
            let dir_eid = eid.next();
            let d = dir.normalize();
            writeln!(writer, "#{dir_eid}=DIRECTION('',({:.6},{:.6},{:.6}));", d.x, d.y, d.z)?;
            let vec_eid = eid.next();
            writeln!(writer, "#{vec_eid}=VECTOR('',#{dir_eid},{:.6});", dir.norm())?;
            let line_eid = eid.next();
            writeln!(writer, "#{line_eid}=LINE('',#{pt_eid},#{vec_eid});")?;
            Ok(line_eid)
        }
        Curve3::Circle { center, axis, radius } => {
            let placement_eid = write_axis2_placement(writer, center, axis, eid)?;
            let circ_eid = eid.next();
            writeln!(writer, "#{circ_eid}=CIRCLE('',#{placement_eid},{radius:.6});")?;
            Ok(circ_eid)
        }
        Curve3::Ellipse { center, major, minor } => {
            let axis = major.cross(minor);
            let len = axis.norm();
            let axis_n = if len > 1e-15 { axis / len } else { nalgebra::Vector3::new(0.0, 0.0, 1.0) };
            let placement_eid = write_axis2_placement(writer, center, &axis_n, eid)?;
            let ell_eid = eid.next();
            writeln!(writer, "#{ell_eid}=ELLIPSE('',#{placement_eid},{:.6},{:.6});",
                     major.norm(), minor.norm())?;
            Ok(ell_eid)
        }
        Curve3::NurbsCurve(nurbs) => {
            // B_SPLINE_CURVE_WITH_KNOTS
            let n_pts = nurbs.control_points.len();
            let mut pt_eids = Vec::with_capacity(n_pts);
            for pt in &nurbs.control_points {
                let pe = eid.next();
                writeln!(writer, "#{pe}=CARTESIAN_POINT('',({:.6},{:.6},{:.6}));", pt.x, pt.y, pt.z)?;
                pt_eids.push(pe);
            }

            let pt_refs: Vec<String> = pt_eids.iter().map(|e| format!("#{e}")).collect();
            let knot_mults = compute_knot_multiplicities(&nurbs.knots);
            let mult_vals: Vec<String> = knot_mults.iter().map(|(_, m)| m.to_string()).collect();
            let unique_knots: Vec<String> = knot_mults.iter().map(|(k, _)| format!("{k:.6}")).collect();

            let bsc_eid = eid.next();
            writeln!(writer, "#{bsc_eid}=B_SPLINE_CURVE_WITH_KNOTS('',{degree},({pts}),.UNSPECIFIED.,.F.,.F.,({mults}),({knots}),.UNSPECIFIED.);",
                     degree = nurbs.degree,
                     pts = pt_refs.join(","),
                     mults = mult_vals.join(","),
                     knots = unique_knots.join(","))?;
            Ok(bsc_eid)
        }
    }
}

fn write_surface_entity<W: Write>(
    writer: &mut W,
    surface: &Surface,
    eid: &mut EntityCounter,
) -> io::Result<u64> {
    match surface {
        Surface::Plane { origin, normal } => {
            let placement_eid = write_axis2_placement(writer, origin, normal, eid)?;
            let plane_eid = eid.next();
            writeln!(writer, "#{plane_eid}=PLANE('',#{placement_eid});")?;
            Ok(plane_eid)
        }
        Surface::Cylinder { origin, axis, radius } => {
            let placement_eid = write_axis2_placement(writer, origin, axis, eid)?;
            let cyl_eid = eid.next();
            writeln!(writer, "#{cyl_eid}=CYLINDRICAL_SURFACE('',#{placement_eid},{radius:.6});")?;
            Ok(cyl_eid)
        }
        Surface::Cone { apex, axis, half_angle } => {
            let placement_eid = write_axis2_placement(writer, apex, axis, eid)?;
            let cone_eid = eid.next();
            writeln!(writer, "#{cone_eid}=CONICAL_SURFACE('',#{placement_eid},0.0,{:.6});",
                     half_angle.to_degrees())?;
            Ok(cone_eid)
        }
        Surface::Sphere { center, radius } => {
            let axis = nalgebra::Vector3::new(0.0, 0.0, 1.0);
            let placement_eid = write_axis2_placement(writer, center, &axis, eid)?;
            let sph_eid = eid.next();
            writeln!(writer, "#{sph_eid}=SPHERICAL_SURFACE('',#{placement_eid},{radius:.6});")?;
            Ok(sph_eid)
        }
        Surface::Torus { center, axis, major_r, minor_r } => {
            let placement_eid = write_axis2_placement(writer, center, axis, eid)?;
            let tor_eid = eid.next();
            writeln!(writer, "#{tor_eid}=TOROIDAL_SURFACE('',#{placement_eid},{major_r:.6},{minor_r:.6});")?;
            Ok(tor_eid)
        }
        Surface::NurbsSurface(nurbs) => {
            // B_SPLINE_SURFACE_WITH_KNOTS
            let n_u = nurbs.control_points.len();
            let _n_v = if n_u > 0 { nurbs.control_points[0].len() } else { 0 };

            let mut pt_eids = Vec::new();
            for row in &nurbs.control_points {
                let mut row_eids = Vec::new();
                for pt in row {
                    let pe = eid.next();
                    writeln!(writer, "#{pe}=CARTESIAN_POINT('',({:.6},{:.6},{:.6}));", pt.x, pt.y, pt.z)?;
                    row_eids.push(pe);
                }
                pt_eids.push(row_eids);
            }

            let rows: Vec<String> = pt_eids.iter().map(|row| {
                let refs: Vec<String> = row.iter().map(|e| format!("#{e}")).collect();
                format!("({})", refs.join(","))
            }).collect();

            let knot_mults_u = compute_knot_multiplicities(&nurbs.knots_u);
            let knot_mults_v = compute_knot_multiplicities(&nurbs.knots_v);

            let bss_eid = eid.next();
            writeln!(writer, "#{bss_eid}=B_SPLINE_SURFACE_WITH_KNOTS('',{du},{dv},({pts}),.UNSPECIFIED.,.F.,.F.,.F.,({mu}),({mv}),({ku}),({kv}),.UNSPECIFIED.);",
                     du = nurbs.degree_u,
                     dv = nurbs.degree_v,
                     pts = rows.join(","),
                     mu = knot_mults_u.iter().map(|(_, m)| m.to_string()).collect::<Vec<_>>().join(","),
                     mv = knot_mults_v.iter().map(|(_, m)| m.to_string()).collect::<Vec<_>>().join(","),
                     ku = knot_mults_u.iter().map(|(k, _)| format!("{k:.6}")).collect::<Vec<_>>().join(","),
                     kv = knot_mults_v.iter().map(|(k, _)| format!("{k:.6}")).collect::<Vec<_>>().join(","),
            )?;
            Ok(bss_eid)
        }
    }
}

fn write_axis2_placement<W: Write>(
    writer: &mut W,
    origin: &nalgebra::Point3<f64>,
    axis: &nalgebra::Vector3<f64>,
    eid: &mut EntityCounter,
) -> io::Result<u64> {
    let pt_eid = eid.next();
    writeln!(writer, "#{pt_eid}=CARTESIAN_POINT('',({:.6},{:.6},{:.6}));", origin.x, origin.y, origin.z)?;

    let dir_eid = eid.next();
    let a = axis.normalize();
    writeln!(writer, "#{dir_eid}=DIRECTION('',({:.6},{:.6},{:.6}));", a.x, a.y, a.z)?;

    // Compute ref_direction perpendicular to axis
    let seed = if a.x.abs() < 0.9 { nalgebra::Vector3::new(1.0, 0.0, 0.0) }
               else { nalgebra::Vector3::new(0.0, 1.0, 0.0) };
    let ref_dir = a.cross(&seed).normalize();
    let ref_eid = eid.next();
    writeln!(writer, "#{ref_eid}=DIRECTION('',({:.6},{:.6},{:.6}));", ref_dir.x, ref_dir.y, ref_dir.z)?;

    let placement_eid = eid.next();
    writeln!(writer, "#{placement_eid}=AXIS2_PLACEMENT_3D('',#{pt_eid},#{dir_eid},#{ref_eid});")?;

    Ok(placement_eid)
}

/// Compute knot multiplicities from a flat knot vector.
fn compute_knot_multiplicities(knots: &[f64]) -> Vec<(f64, usize)> {
    let mut result = Vec::new();
    if knots.is_empty() {
        return result;
    }
    let mut current = knots[0];
    let mut count = 1;
    for &k in &knots[1..] {
        if (k - current).abs() < 1e-10 {
            count += 1;
        } else {
            result.push((current, count));
            current = k;
            count = 1;
        }
    }
    result.push((current, count));
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::Shape;

    #[test]
    fn step_box_contains_entities() {
        let shape = Shape::box3(5.0, 3.0, 8.0);
        let mut buf = Vec::new();
        write_step(&shape.store, shape.solid, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        assert!(text.contains("ISO-10303-21;"));
        assert!(text.contains("MANIFOLD_SOLID_BREP"));
        assert!(text.contains("CLOSED_SHELL"));
        assert!(text.contains("ADVANCED_FACE"));
        assert!(text.contains("PLANE"));
        assert!(text.contains("CARTESIAN_POINT"));
        assert!(text.contains("END-ISO-10303-21;"));
    }

    #[test]
    fn step_sphere_contains_spherical_surface() {
        let shape = Shape::sphere(10.0);
        let mut buf = Vec::new();
        write_step(&shape.store, shape.solid, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        assert!(text.contains("SPHERICAL_SURFACE"));
    }

    #[test]
    fn step_cylinder_contains_cylindrical_surface() {
        let shape = Shape::cylinder(5.0, 20.0);
        let mut buf = Vec::new();
        write_step(&shape.store, shape.solid, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        assert!(text.contains("CYLINDRICAL_SURFACE"));
        assert!(text.contains("PLANE")); // caps
    }

    #[test]
    fn step_entity_count_for_box() {
        let shape = Shape::box3(5.0, 3.0, 8.0);
        let mut buf = Vec::new();
        write_step(&shape.store, shape.solid, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        // Count ADVANCED_FACE entries — should be 6 for a box
        let face_count = text.lines().filter(|l| l.contains("ADVANCED_FACE")).count();
        assert_eq!(face_count, 6, "Box should have 6 ADVANCED_FACE entities");
    }
}
