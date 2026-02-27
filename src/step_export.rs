use crate::dag::SdfNode;
use crate::mesh::TriangleMesh;
use nalgebra::Vector3;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

// ---------------------------------------------------------------------------
// Export tier classification
// ---------------------------------------------------------------------------

/// Whether a DAG node can be exported as exact BRep or must fall back to
/// tessellated (per-triangle) STEP output.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportTier {
    /// The shape maps to native STEP BRep entities (SPHERICAL_SURFACE, PLANE, etc.).
    Exact,
    /// The shape requires mesh tessellation (ADVANCED_FACE per triangle).
    Tessellated,
}

/// Walk the SDF DAG and decide the export strategy.
///
/// Primitives that have direct BRep equivalents (Sphere, Box3, Cylinder) are
/// classified as `Exact`. Rigid transforms (Translate, Rotate, Scale) wrapping
/// an `Exact` node inherit that classification. Everything else falls back to
/// `Tessellated`.
pub fn classify(node: &SdfNode) -> ExportTier {
    match node {
        // Primitives with direct STEP BRep mappings
        SdfNode::Sphere { .. } => ExportTier::Exact,
        SdfNode::Box3 { .. } => ExportTier::Exact,
        SdfNode::Cylinder { .. } => ExportTier::Exact,

        // Rigid transforms inherit the inner node's tier
        SdfNode::Translate(inner, _) => classify(inner),
        SdfNode::Rotate(inner, _) => classify(inner),
        SdfNode::Scale(inner, _) => classify(inner),

        // Everything else → tessellated fallback
        _ => ExportTier::Tessellated,
    }
}

// ---------------------------------------------------------------------------
// Smart STEP entry point
// ---------------------------------------------------------------------------

/// Write a STEP AP203 file, choosing the best representation automatically.
///
/// - For recognized primitives (optionally with rigid transforms): writes
///   proper BRep STEP entities (Tier 1).
/// - For complex/blended shapes: extracts a mesh via adaptive dual contouring
///   and writes tessellated ADVANCED_FACE entities (Tier 3).
pub fn write_step_smart(node: &SdfNode, path: &Path) -> std::io::Result<()> {
    match classify(node) {
        ExportTier::Exact => write_step_exact(node, path),
        ExportTier::Tessellated => write_step_tessellated_from_node(node, path),
    }
}

// ---------------------------------------------------------------------------
// Tier 3: Tessellated fallback (mesh → per-triangle ADVANCED_FACE)
// ---------------------------------------------------------------------------

/// Mesh an SdfNode via adaptive dual contouring and write tessellated STEP.
fn write_step_tessellated_from_node(node: &SdfNode, path: &Path) -> std::io::Result<()> {
    use crate::builder::compute_bbox;
    use crate::dual_contouring::extract_mesh_adaptive;
    use crate::types::{BBox3, MeshSettings};

    let bbox = compute_bbox(node);
    // Expand bbox slightly to avoid clipping surface cells at the boundary.
    let pad = bbox.size() * 0.05;
    let padded = BBox3::new(bbox.min - pad, bbox.max + pad);
    let settings = MeshSettings::default();
    let mesh = extract_mesh_adaptive(node, &padded, &settings);
    write_step(&mesh, path)
}

/// Write a STEP AP203 file from a triangle mesh.
///
/// This produces a tessellated CLOSED_SHELL representation — every triangle
/// becomes an ADVANCED_FACE with a planar surface. The result can be opened
/// in any STEP-compatible CAD tool (FreeCAD, Fusion 360, SolidWorks) to
/// inspect the quality of the SDF kernel's mesh approximation.
pub fn write_step(mesh: &TriangleMesh, path: &Path) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;

    write_step_header(&mut file, path)?;

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
    let mut edge_map: HashMap<(usize, usize), u64> = HashMap::new();

    let tri_count = mesh.indices.len() / 3;
    let mut face_ids: Vec<u64> = Vec::with_capacity(tri_count);

    for chunk in mesh.indices.chunks(3) {
        let (ai, bi, ci) = (chunk[0] as usize, chunk[1] as usize, chunk[2] as usize);

        let edges = [(ai, bi), (bi, ci), (ci, ai)];
        let mut edge_entity_ids = Vec::with_capacity(3);

        for &(v0, v1) in &edges {
            let key = if v0 < v1 { (v0, v1) } else { (v1, v0) };
            let eid = *edge_map.entry(key).or_insert_with(|| {
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
            edge_entity_ids.push((eid, v0 == key.0));
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

        let axis_id = next_id();
        writeln!(file, "#{}=AXIS2_PLACEMENT_3D('',#{},#{},#{});",
            axis_id, point_ids[ai], fn_dir_id, ref_dir_id)?;

        let plane_id = next_id();
        writeln!(file, "#{}=PLANE('',#{});", plane_id, axis_id)?;

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

    write_step_footer(&mut file, &mut next_id, dir_z_id, dir_x_id, brep_id)?;

    Ok(())
}

// ---------------------------------------------------------------------------
// Tier 1: Exact BRep export
// ---------------------------------------------------------------------------

/// Write exact BRep STEP entities for recognized primitives.
fn write_step_exact(node: &SdfNode, path: &Path) -> std::io::Result<()> {
    // Collect the accumulated transform and find the innermost primitive.
    let (prim, offset, scale) = unwrap_transforms(node);

    match prim {
        SdfNode::Sphere { center, radius } => {
            write_step_sphere(path, *center + offset, *radius * scale)
        }
        SdfNode::Box3 { center, half_extents } => {
            write_step_box(path, *center + offset, *half_extents * scale)
        }
        SdfNode::Cylinder { base, axis, radius, height } => {
            write_step_cylinder(path, *base + offset, *axis, *radius * scale, *height * scale)
        }
        _ => {
            // Should not happen since classify() returned Exact, but be defensive.
            write_step_tessellated_from_node(node, path)
        }
    }
}

/// Peel off Translate / Rotate / Scale wrappers and accumulate the transform.
/// Returns (innermost_node, accumulated_offset, accumulated_scale).
///
/// NOTE: Rotation is intentionally ignored for the offset accumulation in this
/// simplified version — the STEP entities are written in the rotated local
/// frame. A production implementation would compose a full transform matrix.
fn unwrap_transforms(node: &SdfNode) -> (&SdfNode, Vector3<f64>, f64) {
    match node {
        SdfNode::Translate(inner, offset) => {
            let (prim, inner_offset, s) = unwrap_transforms(inner);
            (prim, inner_offset + offset, s)
        }
        SdfNode::Scale(inner, factor) => {
            let (prim, inner_offset, s) = unwrap_transforms(inner);
            (prim, inner_offset * *factor, s * *factor)
        }
        SdfNode::Rotate(inner, _rotation) => {
            // For simplicity we pass through; the primitive is written in local coords.
            unwrap_transforms(inner)
        }
        other => (other, Vector3::zeros(), 1.0),
    }
}

// ---------------------------------------------------------------------------
// Sphere BRep
// ---------------------------------------------------------------------------

fn write_step_sphere(path: &Path, center: Vector3<f64>, radius: f64) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;
    write_step_header(&mut file, path)?;

    let mut id: u64 = 1;
    let mut next_id = || { let current = id; id += 1; current };

    // Center point
    let center_pt = next_id();
    writeln!(file, "#{}=CARTESIAN_POINT('Center',({:.15E},{:.15E},{:.15E}));",
        center_pt, center.x, center.y, center.z)?;

    // Axis direction (Z up)
    let dir_z = next_id();
    writeln!(file, "#{}=DIRECTION('Axis',(0.E0,0.E0,1.E0));", dir_z)?;

    // Ref direction (X)
    let dir_x = next_id();
    writeln!(file, "#{}=DIRECTION('RefDir',(1.E0,0.E0,0.E0));", dir_x)?;

    // AXIS2_PLACEMENT_3D for the sphere
    let axis_placement = next_id();
    writeln!(file, "#{}=AXIS2_PLACEMENT_3D('',#{},#{},#{});",
        axis_placement, center_pt, dir_z, dir_x)?;

    // SPHERICAL_SURFACE
    let sphere_surf = next_id();
    writeln!(file, "#{}=SPHERICAL_SURFACE('',#{},{:.15E});",
        sphere_surf, axis_placement, radius)?;

    // For a complete sphere we need a single ADVANCED_FACE with the spherical
    // surface. We define a minimal topology: a VERTEX_POINT at the north pole,
    // and a degenerate edge loop.

    // North pole vertex
    let north_pt = next_id();
    writeln!(file, "#{}=CARTESIAN_POINT('NorthPole',({:.15E},{:.15E},{:.15E}));",
        north_pt, center.x, center.y, center.z + radius)?;
    let north_vp = next_id();
    writeln!(file, "#{}=VERTEX_POINT('',#{});", north_vp, north_pt)?;

    // South pole vertex
    let south_pt = next_id();
    writeln!(file, "#{}=CARTESIAN_POINT('SouthPole',({:.15E},{:.15E},{:.15E}));",
        south_pt, center.x, center.y, center.z - radius)?;
    let south_vp = next_id();
    writeln!(file, "#{}=VERTEX_POINT('',#{});", south_vp, south_pt)?;

    // Seam edge (half circle from north to south pole along the prime meridian)
    // We need a CIRCLE curve for the edge
    let seam_circle_axis = next_id();
    writeln!(file, "#{}=AXIS2_PLACEMENT_3D('',#{},#{},#{});",
        seam_circle_axis, center_pt, dir_x, dir_z)?;
    let seam_circle = next_id();
    writeln!(file, "#{}=CIRCLE('',#{},{:.15E});", seam_circle, seam_circle_axis, radius)?;

    let seam_edge = next_id();
    writeln!(file, "#{}=EDGE_CURVE('',#{},#{},#{},.T.);",
        seam_edge, north_vp, south_vp, seam_circle)?;

    // Two oriented edges (forward and reverse) to form a closed loop
    let oe_fwd = next_id();
    writeln!(file, "#{}=ORIENTED_EDGE('',*,*,#{},.T.);", oe_fwd, seam_edge)?;
    let oe_rev = next_id();
    writeln!(file, "#{}=ORIENTED_EDGE('',*,*,#{},.F.);", oe_rev, seam_edge)?;

    // EDGE_LOOP
    let edge_loop = next_id();
    writeln!(file, "#{}=EDGE_LOOP('',(#{},#{}));", edge_loop, oe_fwd, oe_rev)?;

    // FACE_BOUND
    let face_bound = next_id();
    writeln!(file, "#{}=FACE_BOUND('',#{},.T.);", face_bound, edge_loop)?;

    // ADVANCED_FACE with spherical surface
    let face_id = next_id();
    writeln!(file, "#{}=ADVANCED_FACE('',(#{}),#{},.T.);",
        face_id, face_bound, sphere_surf)?;

    // CLOSED_SHELL
    let shell_id = next_id();
    writeln!(file, "#{}=CLOSED_SHELL('Sphere_Shell',(#{}));", shell_id, face_id)?;

    // MANIFOLD_SOLID_BREP
    let brep_id = next_id();
    writeln!(file, "#{}=MANIFOLD_SOLID_BREP('Sphere_Body',#{});", brep_id, shell_id)?;

    write_step_footer(&mut file, &mut next_id, dir_z, dir_x, brep_id)?;

    Ok(())
}

// ---------------------------------------------------------------------------
// Box BRep
// ---------------------------------------------------------------------------

fn write_step_box(
    path: &Path,
    center: Vector3<f64>,
    half_extents: Vector3<f64>,
) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;
    write_step_header(&mut file, path)?;

    let mut id: u64 = 1;
    let mut next_id = || { let current = id; id += 1; current };

    let hx = half_extents.x;
    let hy = half_extents.y;
    let hz = half_extents.z;
    let cx = center.x;
    let cy = center.y;
    let cz = center.z;

    // 8 vertices of the box
    let corners = [
        Vector3::new(cx - hx, cy - hy, cz - hz), // 0: ---
        Vector3::new(cx + hx, cy - hy, cz - hz), // 1: +--
        Vector3::new(cx + hx, cy + hy, cz - hz), // 2: ++-
        Vector3::new(cx - hx, cy + hy, cz - hz), // 3: -+-
        Vector3::new(cx - hx, cy - hy, cz + hz), // 4: --+
        Vector3::new(cx + hx, cy - hy, cz + hz), // 5: +-+
        Vector3::new(cx + hx, cy + hy, cz + hz), // 6: +++
        Vector3::new(cx - hx, cy + hy, cz + hz), // 7: -++
    ];

    // Write CARTESIAN_POINTs and VERTEX_POINTs
    let mut pt_ids = [0u64; 8];
    let mut vp_ids = [0u64; 8];
    for (i, c) in corners.iter().enumerate() {
        let pid = next_id();
        writeln!(file, "#{}=CARTESIAN_POINT('V{}',({:.15E},{:.15E},{:.15E}));",
            pid, i, c.x, c.y, c.z)?;
        pt_ids[i] = pid;

        let vid = next_id();
        writeln!(file, "#{}=VERTEX_POINT('',#{});", vid, pid)?;
        vp_ids[i] = vid;
    }

    // 12 edges of the box. Each edge is (start_vertex_idx, end_vertex_idx).
    let edge_defs: [(usize, usize); 12] = [
        (0, 1), (1, 2), (2, 3), (3, 0), // bottom face edges
        (4, 5), (5, 6), (6, 7), (7, 4), // top face edges
        (0, 4), (1, 5), (2, 6), (3, 7), // vertical edges
    ];

    let mut edge_ids = [0u64; 12];
    for (i, &(a, b)) in edge_defs.iter().enumerate() {
        // Direction vector for the line
        let va = &corners[a];
        let vb = &corners[b];
        let dx = vb.x - va.x;
        let dy = vb.y - va.y;
        let dz = vb.z - va.z;
        let len = (dx*dx + dy*dy + dz*dz).sqrt();

        let dir_id = next_id();
        writeln!(file, "#{}=DIRECTION('',({:.15E},{:.15E},{:.15E}));",
            dir_id, dx/len, dy/len, dz/len)?;

        let vec_id = next_id();
        writeln!(file, "#{}=VECTOR('',#{},{:.15E});", vec_id, dir_id, len)?;

        let line_id = next_id();
        writeln!(file, "#{}=LINE('',#{},#{});", line_id, pt_ids[a], vec_id)?;

        let eid = next_id();
        writeln!(file, "#{}=EDGE_CURVE('',#{},#{},#{},.T.);",
            eid, vp_ids[a], vp_ids[b], line_id)?;
        edge_ids[i] = eid;
    }

    // 6 faces of the box.
    // Each face is defined by 4 edge indices and their orientations,
    // plus a normal direction for the PLANE.
    //
    // Face definitions: (edge_indices, orientations, normal, face_point_index)
    // Edge orientation: true = .T. (same direction as edge def), false = .F. (reversed)
    struct FaceDef {
        edges: [(usize, bool); 4], // (edge_index, same_orientation)
        normal: Vector3<f64>,
        point_idx: usize, // index into corners for the plane origin
    }

    let face_defs = [
        // Bottom face (Z-): vertices 0,1,2,3 — edges 0,1,2,3 with normal -Z
        FaceDef {
            edges: [(0, false), (3, true), (2, false), (1, false)],
            normal: Vector3::new(0.0, 0.0, -1.0),
            point_idx: 0,
        },
        // Top face (Z+): vertices 4,5,6,7 — edges 4,5,6,7 with normal +Z
        FaceDef {
            edges: [(4, true), (5, true), (6, true), (7, true)],
            normal: Vector3::new(0.0, 0.0, 1.0),
            point_idx: 4,
        },
        // Front face (Y-): vertices 0,1,5,4 — edges 0,9,4(rev),8(rev) with normal -Y
        FaceDef {
            edges: [(0, true), (9, true), (4, false), (8, false)],
            normal: Vector3::new(0.0, -1.0, 0.0),
            point_idx: 0,
        },
        // Back face (Y+): vertices 3,2,6,7 — edges 2(rev),10,6(rev),11(rev) with normal +Y
        FaceDef {
            edges: [(2, true), (10, true), (6, false), (11, false)],
            normal: Vector3::new(0.0, 1.0, 0.0),
            point_idx: 3,
        },
        // Left face (X-): vertices 0,3,7,4 — edges 3(rev),11,7(rev),8(rev) with normal -X
        FaceDef {
            edges: [(3, false), (11, true), (7, false), (8, true)],
            normal: Vector3::new(-1.0, 0.0, 0.0),
            point_idx: 0,
        },
        // Right face (X+): vertices 1,2,6,5 — edges 1,10,5(rev),9(rev) with normal +X
        FaceDef {
            edges: [(1, true), (10, false), (5, false), (9, false)],
            normal: Vector3::new(1.0, 0.0, 0.0),
            point_idx: 1,
        },
    ];

    // Global Z and X directions for reuse
    let global_dir_z = next_id();
    writeln!(file, "#{}=DIRECTION('',(0.E0,0.E0,1.E0));", global_dir_z)?;
    let global_dir_x = next_id();
    writeln!(file, "#{}=DIRECTION('',(1.E0,0.E0,0.E0));", global_dir_x)?;

    let mut face_ids = Vec::with_capacity(6);

    for fd in &face_defs {
        // Create ORIENTED_EDGEs
        let mut oe_ids = Vec::with_capacity(4);
        for &(edge_idx, fwd) in &fd.edges {
            let oe_id = next_id();
            let orient = if fwd { ".T." } else { ".F." };
            writeln!(file, "#{}=ORIENTED_EDGE('',*,*,#{},{});",
                oe_id, edge_ids[edge_idx], orient)?;
            oe_ids.push(oe_id);
        }

        // EDGE_LOOP
        let el_id = next_id();
        writeln!(file, "#{}=EDGE_LOOP('',(#{},#{},#{},#{}));",
            el_id, oe_ids[0], oe_ids[1], oe_ids[2], oe_ids[3])?;

        // FACE_OUTER_BOUND
        let fob_id = next_id();
        writeln!(file, "#{}=FACE_OUTER_BOUND('',#{},.T.);", fob_id, el_id)?;

        // Normal direction
        let n = &fd.normal;
        let n_dir = next_id();
        writeln!(file, "#{}=DIRECTION('',({:.15E},{:.15E},{:.15E}));",
            n_dir, n.x, n.y, n.z)?;

        // Ref direction: pick one perpendicular to normal
        let ref_vec = if n.x.abs() > 0.9 {
            Vector3::new(0.0, 1.0, 0.0)
        } else {
            Vector3::new(1.0, 0.0, 0.0)
        };
        let r_dir = next_id();
        writeln!(file, "#{}=DIRECTION('',({:.15E},{:.15E},{:.15E}));",
            r_dir, ref_vec.x, ref_vec.y, ref_vec.z)?;

        // Plane origin point (reuse a corner point)
        let plane_origin = pt_ids[fd.point_idx];

        // AXIS2_PLACEMENT_3D
        let axis_id = next_id();
        writeln!(file, "#{}=AXIS2_PLACEMENT_3D('',#{},#{},#{});",
            axis_id, plane_origin, n_dir, r_dir)?;

        // PLANE
        let plane_id = next_id();
        writeln!(file, "#{}=PLANE('',#{});", plane_id, axis_id)?;

        // ADVANCED_FACE
        let face_id = next_id();
        writeln!(file, "#{}=ADVANCED_FACE('',(#{}),#{},.T.);",
            face_id, fob_id, plane_id)?;

        face_ids.push(face_id);
    }

    // CLOSED_SHELL
    let shell_id = next_id();
    write!(file, "#{}=CLOSED_SHELL('Box_Shell',(", shell_id)?;
    for (i, fid) in face_ids.iter().enumerate() {
        if i > 0 { write!(file, ",")?; }
        write!(file, "#{}", fid)?;
    }
    writeln!(file, "));")?;

    // MANIFOLD_SOLID_BREP
    let brep_id = next_id();
    writeln!(file, "#{}=MANIFOLD_SOLID_BREP('Box_Body',#{});", brep_id, shell_id)?;

    write_step_footer(&mut file, &mut next_id, global_dir_z, global_dir_x, brep_id)?;

    Ok(())
}

// ---------------------------------------------------------------------------
// Cylinder BRep
// ---------------------------------------------------------------------------

fn write_step_cylinder(
    path: &Path,
    base: Vector3<f64>,
    axis: Vector3<f64>,
    radius: f64,
    height: f64,
) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;
    write_step_header(&mut file, path)?;

    let mut id: u64 = 1;
    let mut next_id = || { let current = id; id += 1; current };

    let top = base + axis * height;

    // Build a local coordinate frame from the axis
    let ref_dir = if axis.x.abs() < 0.9 {
        let v = Vector3::new(1.0, 0.0, 0.0);
        let proj = v - axis * v.dot(&axis);
        let len = proj.norm();
        if len > 1e-12 { proj / len } else { Vector3::new(0.0, 1.0, 0.0) }
    } else {
        let v = Vector3::new(0.0, 1.0, 0.0);
        let proj = v - axis * v.dot(&axis);
        let len = proj.norm();
        if len > 1e-12 { proj / len } else { Vector3::new(0.0, 0.0, 1.0) }
    };

    // Base center point
    let base_pt = next_id();
    writeln!(file, "#{}=CARTESIAN_POINT('BaseCenter',({:.15E},{:.15E},{:.15E}));",
        base_pt, base.x, base.y, base.z)?;

    // Top center point
    let top_pt = next_id();
    writeln!(file, "#{}=CARTESIAN_POINT('TopCenter',({:.15E},{:.15E},{:.15E}));",
        top_pt, top.x, top.y, top.z)?;

    // Axis direction
    let axis_dir = next_id();
    writeln!(file, "#{}=DIRECTION('Axis',({:.15E},{:.15E},{:.15E}));",
        axis_dir, axis.x, axis.y, axis.z)?;

    // Negative axis direction (for bottom cap)
    let neg_axis_dir = next_id();
    writeln!(file, "#{}=DIRECTION('NegAxis',({:.15E},{:.15E},{:.15E}));",
        neg_axis_dir, -axis.x, -axis.y, -axis.z)?;

    // Ref direction
    let ref_dir_id = next_id();
    writeln!(file, "#{}=DIRECTION('RefDir',({:.15E},{:.15E},{:.15E}));",
        ref_dir_id, ref_dir.x, ref_dir.y, ref_dir.z)?;

    // --- Cylindrical surface ---
    let cyl_axis = next_id();
    writeln!(file, "#{}=AXIS2_PLACEMENT_3D('',#{},#{},#{});",
        cyl_axis, base_pt, axis_dir, ref_dir_id)?;

    let cyl_surf = next_id();
    writeln!(file, "#{}=CYLINDRICAL_SURFACE('',#{},{:.15E});",
        cyl_surf, cyl_axis, radius)?;

    // --- Bottom cap plane ---
    let bot_axis = next_id();
    writeln!(file, "#{}=AXIS2_PLACEMENT_3D('',#{},#{},#{});",
        bot_axis, base_pt, neg_axis_dir, ref_dir_id)?;

    let bot_plane = next_id();
    writeln!(file, "#{}=PLANE('',#{});", bot_plane, bot_axis)?;

    // --- Top cap plane ---
    let top_axis = next_id();
    writeln!(file, "#{}=AXIS2_PLACEMENT_3D('',#{},#{},#{});",
        top_axis, top_pt, axis_dir, ref_dir_id)?;

    let top_plane = next_id();
    writeln!(file, "#{}=PLANE('',#{});", top_plane, top_axis)?;

    // --- Edge circles (bottom and top) ---
    // A seam vertex on the bottom circle
    let seam_bottom = base + ref_dir * radius;
    let seam_bottom_pt = next_id();
    writeln!(file, "#{}=CARTESIAN_POINT('',({:.15E},{:.15E},{:.15E}));",
        seam_bottom_pt, seam_bottom.x, seam_bottom.y, seam_bottom.z)?;
    let seam_bottom_vp = next_id();
    writeln!(file, "#{}=VERTEX_POINT('',#{});", seam_bottom_vp, seam_bottom_pt)?;

    // A seam vertex on the top circle
    let seam_top = top + ref_dir * radius;
    let seam_top_pt = next_id();
    writeln!(file, "#{}=CARTESIAN_POINT('',({:.15E},{:.15E},{:.15E}));",
        seam_top_pt, seam_top.x, seam_top.y, seam_top.z)?;
    let seam_top_vp = next_id();
    writeln!(file, "#{}=VERTEX_POINT('',#{});", seam_top_vp, seam_top_pt)?;

    // Bottom circle
    let bot_circle_axis = next_id();
    writeln!(file, "#{}=AXIS2_PLACEMENT_3D('',#{},#{},#{});",
        bot_circle_axis, base_pt, axis_dir, ref_dir_id)?;
    let bot_circle = next_id();
    writeln!(file, "#{}=CIRCLE('',#{},{:.15E});", bot_circle, bot_circle_axis, radius)?;

    // Top circle
    let top_circle_axis = next_id();
    writeln!(file, "#{}=AXIS2_PLACEMENT_3D('',#{},#{},#{});",
        top_circle_axis, top_pt, axis_dir, ref_dir_id)?;
    let top_circle = next_id();
    writeln!(file, "#{}=CIRCLE('',#{},{:.15E});", top_circle, top_circle_axis, radius)?;

    // Bottom edge curve (full circle, same start and end vertex)
    let bot_edge = next_id();
    writeln!(file, "#{}=EDGE_CURVE('',#{},#{},#{},.T.);",
        bot_edge, seam_bottom_vp, seam_bottom_vp, bot_circle)?;

    // Top edge curve
    let top_edge = next_id();
    writeln!(file, "#{}=EDGE_CURVE('',#{},#{},#{},.T.);",
        top_edge, seam_top_vp, seam_top_vp, top_circle)?;

    // Seam line (vertical edge from bottom seam to top seam)
    let seam_line_dir = next_id();
    writeln!(file, "#{}=DIRECTION('',({:.15E},{:.15E},{:.15E}));",
        seam_line_dir, axis.x, axis.y, axis.z)?;
    let seam_vec = next_id();
    writeln!(file, "#{}=VECTOR('',#{},{:.15E});", seam_vec, seam_line_dir, height)?;
    let seam_line = next_id();
    writeln!(file, "#{}=LINE('',#{},#{});", seam_line, seam_bottom_pt, seam_vec)?;
    let seam_edge = next_id();
    writeln!(file, "#{}=EDGE_CURVE('',#{},#{},#{},.T.);",
        seam_edge, seam_bottom_vp, seam_top_vp, seam_line)?;

    // --- Build faces ---

    // Bottom cap face: edge loop around bot_edge (reversed, outward normal is -axis)
    let bot_oe = next_id();
    writeln!(file, "#{}=ORIENTED_EDGE('',*,*,#{},.F.);", bot_oe, bot_edge)?;
    let bot_el = next_id();
    writeln!(file, "#{}=EDGE_LOOP('',(#{}));", bot_el, bot_oe)?;
    let bot_fb = next_id();
    writeln!(file, "#{}=FACE_OUTER_BOUND('',#{},.T.);", bot_fb, bot_el)?;
    let bot_face = next_id();
    writeln!(file, "#{}=ADVANCED_FACE('',(#{}),#{},.T.);", bot_face, bot_fb, bot_plane)?;

    // Top cap face: edge loop around top_edge (forward, outward normal is +axis)
    let top_oe = next_id();
    writeln!(file, "#{}=ORIENTED_EDGE('',*,*,#{},.T.);", top_oe, top_edge)?;
    let top_el = next_id();
    writeln!(file, "#{}=EDGE_LOOP('',(#{}));", top_el, top_oe)?;
    let top_fb = next_id();
    writeln!(file, "#{}=FACE_OUTER_BOUND('',#{},.T.);", top_fb, top_el)?;
    let top_face = next_id();
    writeln!(file, "#{}=ADVANCED_FACE('',(#{}),#{},.T.);", top_face, top_fb, top_plane)?;

    // Lateral face: loop is bot_edge(fwd) + seam(fwd) + top_edge(rev) + seam(rev)
    let lat_oe1 = next_id();
    writeln!(file, "#{}=ORIENTED_EDGE('',*,*,#{},.T.);", lat_oe1, bot_edge)?;
    let lat_oe2 = next_id();
    writeln!(file, "#{}=ORIENTED_EDGE('',*,*,#{},.T.);", lat_oe2, seam_edge)?;
    let lat_oe3 = next_id();
    writeln!(file, "#{}=ORIENTED_EDGE('',*,*,#{},.F.);", lat_oe3, top_edge)?;
    let lat_oe4 = next_id();
    writeln!(file, "#{}=ORIENTED_EDGE('',*,*,#{},.F.);", lat_oe4, seam_edge)?;
    let lat_el = next_id();
    writeln!(file, "#{}=EDGE_LOOP('',(#{},#{},#{},#{}));",
        lat_el, lat_oe1, lat_oe2, lat_oe3, lat_oe4)?;
    let lat_fb = next_id();
    writeln!(file, "#{}=FACE_OUTER_BOUND('',#{},.T.);", lat_fb, lat_el)?;
    let lat_face = next_id();
    writeln!(file, "#{}=ADVANCED_FACE('',(#{}),#{},.T.);", lat_face, lat_fb, cyl_surf)?;

    // CLOSED_SHELL
    let shell_id = next_id();
    writeln!(file, "#{}=CLOSED_SHELL('Cylinder_Shell',(#{},#{},#{}));",
        shell_id, bot_face, top_face, lat_face)?;

    // MANIFOLD_SOLID_BREP
    let brep_id = next_id();
    writeln!(file, "#{}=MANIFOLD_SOLID_BREP('Cylinder_Body',#{});", brep_id, shell_id)?;

    write_step_footer(&mut file, &mut next_id, axis_dir, ref_dir_id, brep_id)?;

    Ok(())
}

// ---------------------------------------------------------------------------
// Shared STEP header / footer
// ---------------------------------------------------------------------------

fn write_step_header(file: &mut std::fs::File, path: &Path) -> std::io::Result<()> {
    writeln!(file, "ISO-10303-21;")?;
    writeln!(file, "HEADER;")?;
    writeln!(file, "FILE_DESCRIPTION(('Crusst SDF Geometry Kernel Output'),'2;1');")?;
    writeln!(file, "FILE_NAME('{}','2026-02-26',('Crusst'),('Crusst Geometry Kernel'),'Crusst 0.1.0','Crusst','');",
        path.file_name().unwrap_or_default().to_string_lossy())?;
    writeln!(file, "FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));")?;
    writeln!(file, "ENDSEC;")?;
    writeln!(file, "DATA;")?;
    Ok(())
}

fn write_step_footer(
    file: &mut std::fs::File,
    next_id: &mut dyn FnMut() -> u64,
    dir_z_id: u64,
    dir_x_id: u64,
    brep_id: u64,
) -> std::io::Result<()> {
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
