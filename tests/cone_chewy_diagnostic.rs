/// Diagnostic test for "chewy" (jagged/faceted) CappedCone mesh at base rim and tip.
///
/// Analyzes geometry, normals, and triangle quality for the Monad Order 1 eigenform cone.
use crusst::mesh::extract_mesh;
use crusst::shape::{CappedCone, Sdf};
use nalgebra::Vector3;
use std::collections::HashMap;

#[test]
fn diagnose_chewy_cone() {
    // -----------------------------------------------------------------------
    // 1. Create the cone and mesh (Monad Order 1 eigenform parameters)
    // -----------------------------------------------------------------------
    let cone = CappedCone::new(
        Vector3::zeros(),
        Vector3::new(0.0, 0.0, 30.0),
        12.0,
        0.0,
    );

    let bbox_min = Vector3::new(-14.0, -14.0, -1.0);
    let bbox_max = Vector3::new(14.0, 14.0, 32.0);

    // resolution 128 => depth = log2(128) = 7
    let mesh = extract_mesh(&cone, bbox_min, bbox_max, 128);

    let n_verts = mesh.vertices.len();
    let n_tris = mesh.indices.len() / 3;

    println!("========================================================");
    println!("  CAPPED CONE MESH DIAGNOSTIC (Monad Order 1 eigenform)");
    println!("========================================================");
    println!();
    println!("--- 1. OVERALL MESH STATS ---");
    println!("  Vertices:  {}", n_verts);
    println!("  Triangles: {}", n_tris);
    println!("  Normals:   {}", mesh.normals.len());
    println!();

    // -----------------------------------------------------------------------
    // 2a. Vertices near the base rim (z ~ 0, r ~ 12)
    // -----------------------------------------------------------------------
    println!("--- 2a. BASE RIM VERTICES (z ~ 0, r ~ 12) ---");
    let mut rim_verts: Vec<(usize, f64, f64, f64, f64)> = Vec::new(); // (idx, x, y, z, r)
    for (i, v) in mesh.vertices.iter().enumerate() {
        let r = (v.x * v.x + v.y * v.y).sqrt();
        if v.z.abs() < 1.0 && (r - 12.0).abs() < 2.0 {
            rim_verts.push((i, v.x, v.y, v.z, r));
        }
    }
    println!("  Found {} vertices near base rim", rim_verts.len());

    if !rim_verts.is_empty() {
        let radii: Vec<f64> = rim_verts.iter().map(|v| v.4).collect();
        let z_vals: Vec<f64> = rim_verts.iter().map(|v| v.3).collect();
        let avg_r = radii.iter().sum::<f64>() / radii.len() as f64;
        let min_r = radii.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_r = radii.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_z = z_vals.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_z = z_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let var = radii.iter().map(|r| (r - avg_r).powi(2)).sum::<f64>() / radii.len() as f64;
        let std_dev = var.sqrt();

        println!("  Radius: avg={avg_r:.4}, min={min_r:.4}, max={max_r:.4}, std_dev={std_dev:.4}");
        println!("  Z range: [{min_z:.4}, {max_z:.4}]");
        println!();

        // Print a sample of rim vertices sorted by angle
        let mut by_angle: Vec<(f64, usize, f64, f64, f64, f64)> = rim_verts
            .iter()
            .map(|(i, x, y, z, r)| (y.atan2(*x), *i, *x, *y, *z, *r))
            .collect();
        by_angle.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        println!("  First 20 rim verts (sorted by angle):");
        println!("  {:>6} {:>8} {:>8} {:>8} {:>8} {:>8}", "idx", "x", "y", "z", "r", "angle");
        for (angle, idx, x, y, z, r) in by_angle.iter().take(20) {
            println!(
                "  {:>6} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.3}",
                idx, x, y, z, r, angle.to_degrees()
            );
        }
    }
    println!();

    // -----------------------------------------------------------------------
    // 2b. Vertices near the tip (z ~ 30)
    // -----------------------------------------------------------------------
    println!("--- 2b. TIP VERTICES (z ~ 30) ---");
    let mut tip_verts: Vec<(usize, f64, f64, f64, f64)> = Vec::new();
    for (i, v) in mesh.vertices.iter().enumerate() {
        if (v.z - 30.0).abs() < 2.0 {
            let r = (v.x * v.x + v.y * v.y).sqrt();
            tip_verts.push((i, v.x, v.y, v.z, r));
        }
    }
    println!("  Found {} vertices near tip", tip_verts.len());

    if !tip_verts.is_empty() {
        let z_vals: Vec<f64> = tip_verts.iter().map(|v| v.3).collect();
        let r_vals: Vec<f64> = tip_verts.iter().map(|v| v.4).collect();
        let min_z = z_vals.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_z = z_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let max_r = r_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        println!("  Z range: [{min_z:.4}, {max_z:.4}]");
        println!("  Max radial offset from axis: {max_r:.4}");

        // Check for outliers
        let outliers: Vec<_> = tip_verts.iter().filter(|v| v.4 > 2.0 || (v.3 - 30.0).abs() > 1.0).collect();
        if outliers.is_empty() {
            println!("  No outlier tip vertices.");
        } else {
            println!("  WARNING: {} outlier tip vertices!", outliers.len());
            for (idx, x, y, z, r) in outliers.iter().take(10) {
                println!("    idx={idx} pos=({x:.3}, {y:.3}, {z:.3}) r={r:.3}");
            }
        }

        println!("  All tip verts:");
        println!("  {:>6} {:>8} {:>8} {:>8} {:>8}", "idx", "x", "y", "z", "r");
        for (idx, x, y, z, r) in tip_verts.iter().take(30) {
            println!("  {:>6} {:>8.3} {:>8.3} {:>8.3} {:>8.3}", idx, x, y, z, r);
        }
    }
    println!();

    // -----------------------------------------------------------------------
    // 3. Degenerate triangles (area < 1e-8)
    // -----------------------------------------------------------------------
    println!("--- 3. DEGENERATE TRIANGLES (area < 1e-8) ---");
    let mut degenerate_count = 0;
    let mut total_area = 0.0;
    let mut areas: Vec<f64> = Vec::with_capacity(n_tris);

    for tri in mesh.indices.chunks(3) {
        if tri.len() < 3 { continue; }
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];
        let e1 = v1 - v0;
        let e2 = v2 - v0;
        let cross = e1.cross(&e2);
        let area = cross.norm() * 0.5;
        areas.push(area);
        total_area += area;
        if area < 1e-8 {
            degenerate_count += 1;
        }
    }

    let near_zero_count = areas.iter().filter(|&&a| a < 1e-6).count();
    let small_count = areas.iter().filter(|&&a| a < 1e-4).count();

    areas.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let min_area = areas.first().copied().unwrap_or(0.0);
    let max_area = areas.last().copied().unwrap_or(0.0);
    let median_area = if !areas.is_empty() { areas[areas.len() / 2] } else { 0.0 };

    println!("  Degenerate (area < 1e-8): {degenerate_count}");
    println!("  Near-zero  (area < 1e-6): {near_zero_count}");
    println!("  Small      (area < 1e-4): {small_count}");
    println!("  Total area:  {total_area:.4}");
    println!("  Min area:    {min_area:.2e}");
    println!("  Median area: {median_area:.4}");
    println!("  Max area:    {max_area:.4}");

    // Show the 10 smallest triangles
    if n_tris > 0 {
        println!();
        println!("  10 smallest triangles:");
        let mut indexed_areas: Vec<(usize, f64)> = areas.iter().enumerate().map(|(i, &a)| (i, a)).collect();
        indexed_areas.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        // We lost the original triangle index mapping because areas was sorted. Redo.
        let mut tri_areas: Vec<(usize, f64)> = Vec::with_capacity(n_tris);
        for (ti, tri) in mesh.indices.chunks(3).enumerate() {
            if tri.len() < 3 { continue; }
            let v0 = mesh.vertices[tri[0] as usize];
            let v1 = mesh.vertices[tri[1] as usize];
            let v2 = mesh.vertices[tri[2] as usize];
            let e1 = v1 - v0;
            let e2 = v2 - v0;
            let area = e1.cross(&e2).norm() * 0.5;
            tri_areas.push((ti, area));
        }
        tri_areas.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        for (ti, area) in tri_areas.iter().take(10) {
            let tri = &mesh.indices[ti * 3..ti * 3 + 3];
            let v0 = mesh.vertices[tri[0] as usize];
            let v1 = mesh.vertices[tri[1] as usize];
            let v2 = mesh.vertices[tri[2] as usize];
            println!(
                "    tri[{ti}]: area={area:.2e}, verts=[{},{},{}] at ({:.2},{:.2},{:.2}),({:.2},{:.2},{:.2}),({:.2},{:.2},{:.2})",
                tri[0], tri[1], tri[2],
                v0.x, v0.y, v0.z,
                v1.x, v1.y, v1.z,
                v2.x, v2.y, v2.z
            );
        }
    }
    println!();

    // -----------------------------------------------------------------------
    // 4. Thin/long triangles (aspect ratio > 100)
    // -----------------------------------------------------------------------
    println!("--- 4. THIN/LONG TRIANGLES (aspect ratio) ---");
    let mut thin_count = 0;
    let mut very_thin_count = 0;
    let mut max_aspect = 0.0f64;
    let mut aspect_ratios: Vec<(usize, f64)> = Vec::with_capacity(n_tris);

    for (ti, tri) in mesh.indices.chunks(3).enumerate() {
        if tri.len() < 3 { continue; }
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];

        let a_len = (v1 - v0).norm();
        let b_len = (v2 - v1).norm();
        let c_len = (v0 - v2).norm();
        let longest = a_len.max(b_len).max(c_len);
        let area = (v1 - v0).cross(&(v2 - v0)).norm() * 0.5;

        // Aspect ratio: longest_edge^2 / (4 * sqrt(3) * area) for equilateral = 1
        // Simpler: longest_edge / height where height = 2*area/longest
        let height = if longest > 1e-15 { 2.0 * area / longest } else { 0.0 };
        let aspect = if height > 1e-15 { longest / height } else { f64::INFINITY };

        aspect_ratios.push((ti, aspect));
        if aspect > 100.0 { very_thin_count += 1; }
        if aspect > 10.0 { thin_count += 1; }
        if aspect > max_aspect && aspect.is_finite() { max_aspect = aspect; }
    }

    println!("  Thin triangles (aspect > 10):  {thin_count}");
    println!("  Very thin     (aspect > 100): {very_thin_count}");
    println!("  Max finite aspect ratio:       {max_aspect:.2}");

    aspect_ratios.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    println!("  10 worst aspect ratio triangles:");
    for (ti, ar) in aspect_ratios.iter().take(10) {
        let tri = &mesh.indices[ti * 3..ti * 3 + 3];
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];
        let avg_z = (v0.z + v1.z + v2.z) / 3.0;
        let avg_r = ((v0.x.powi(2) + v0.y.powi(2)).sqrt()
            + (v1.x.powi(2) + v1.y.powi(2)).sqrt()
            + (v2.x.powi(2) + v2.y.powi(2)).sqrt()) / 3.0;
        println!(
            "    tri[{ti}]: aspect={ar:.1}, avg_z={avg_z:.2}, avg_r={avg_r:.2}, verts=[{},{},{}]",
            tri[0], tri[1], tri[2]
        );
    }
    println!();

    // -----------------------------------------------------------------------
    // 5. Vertex normals near the base rim
    // -----------------------------------------------------------------------
    println!("--- 5. VERTEX NORMALS AT BASE RIM ---");
    println!("  Checking normals at rim vertices (z~0, r~12)...");
    println!("  Expected: wall normal has radial + upward components; cap normal is (0,0,-1)");
    println!();

    // Compute expected wall normal direction
    // Wall slant: from (0, ra=12) to (30, rb=0) in (axial, radial) space
    // Wall normal points outward: perpendicular to the wall line
    // Wall line direction = (30, -12) => outward normal = (12, 30) normalized = (12/sqrt(1044), 30/sqrt(1044))
    let slant_len = (12.0_f64.powi(2) + 30.0_f64.powi(2)).sqrt();
    let wall_n_radial = 30.0 / slant_len; // outward radial component
    let wall_n_axial = 12.0 / slant_len;  // upward axial component
    println!("  Expected wall normal (radial, axial): ({wall_n_radial:.4}, {wall_n_axial:.4})");
    println!("  Expected cap-A normal: (0, 0, -1) = downward");
    println!();

    // For each rim vertex, decompose its normal into radial and axial components
    let mut wall_count = 0;
    let mut cap_count = 0;
    let mut ambiguous_count = 0;

    if !rim_verts.is_empty() {
        println!("  {:>6} {:>8} {:>8} {:>8} | {:>8} {:>8} {:>8} | {:>8} {:>8} {:>6}",
            "idx", "x", "y", "z", "nx", "ny", "nz", "n_rad", "n_axial", "type");
        for (idx, x, y, z, r) in &rim_verts {
            let n = &mesh.normals[*idx];
            // Decompose into axial (z) and radial components
            let n_axial = n.z; // component along cone axis (Z)
            let radial_dir = if *r > 1e-10 {
                Vector3::new(*x, *y, 0.0) / *r
            } else {
                Vector3::new(1.0, 0.0, 0.0)
            };
            let n_radial = n.dot(&radial_dir); // outward radial component

            // Classify: wall normal should have n_radial ~ 0.93, n_axial ~ 0.37
            // Cap normal should have n_axial ~ -1.0, n_radial ~ 0.0
            let is_wall = n_radial > 0.5 && n_axial > 0.0;
            let is_cap = n_axial < -0.5;
            let label = if is_wall && !is_cap {
                wall_count += 1;
                "WALL"
            } else if is_cap && !is_wall {
                cap_count += 1;
                "CAP"
            } else {
                ambiguous_count += 1;
                "MIXED!"
            };

            println!(
                "  {:>6} {:>8.3} {:>8.3} {:>8.3} | {:>8.4} {:>8.4} {:>8.4} | {:>8.4} {:>8.4} {:>6}",
                idx, x, y, z, n.x, n.y, n.z, n_radial, n_axial, label
            );
        }
        println!();
        println!("  Normal classification at rim: WALL={wall_count}, CAP={cap_count}, MIXED={ambiguous_count}");
        if ambiguous_count > 0 {
            println!("  *** MIXED normals indicate per-vertex averaging across the sharp edge! ***");
            println!("  *** This is a primary cause of the 'chewy' appearance. ***");
        }
    }
    println!();

    // -----------------------------------------------------------------------
    // 6. Triangles straddling the base rim sharp edge
    // -----------------------------------------------------------------------
    println!("--- 6. TRIANGLES SPANNING THE BASE RIM SHARP EDGE ---");
    println!("  Looking for triangles where some vertices are on the wall side");
    println!("  and others are on the cap side of the rim edge...");
    println!();

    // For each vertex near the rim, classify it
    // wall-side: z > 0 and r ~ 12 (on the cone surface above the base)
    // cap-side: z <= 0 or (z ~ 0 and r < 12) (on the base cap or below)
    let mut vertex_side: HashMap<usize, &str> = HashMap::new();
    for (idx, _x, _y, z, r) in &rim_verts {
        // A vertex is "wall-side" if it's above z=0 with radius close to the cone surface
        // Expected cone radius at z: r_cone(z) = 12 * (1 - z/30) for z in [0, 30]
        let r_cone_at_z = if *z >= 0.0 && *z <= 30.0 { 12.0 * (1.0 - z / 30.0) } else { 0.0 };
        let on_wall = *z > 0.2 || (*z > -0.2 && (*r - r_cone_at_z).abs() < 1.0 && *r > 10.0);
        let on_cap = *z < 0.2 && *r < 11.5;

        if on_wall && !on_cap {
            vertex_side.insert(*idx, "WALL");
        } else if on_cap && !on_wall {
            vertex_side.insert(*idx, "CAP");
        } else {
            vertex_side.insert(*idx, "EDGE"); // right on the rim edge
        }
    }

    let _rim_vert_set: std::collections::HashSet<usize> = rim_verts.iter().map(|v| v.0).collect();
    let mut straddling_count = 0;
    let mut straddling_examples: Vec<String> = Vec::new();

    for (ti, tri) in mesh.indices.chunks(3).enumerate() {
        if tri.len() < 3 { continue; }
        let i0 = tri[0] as usize;
        let i1 = tri[1] as usize;
        let i2 = tri[2] as usize;

        // Check if this triangle has vertices on different sides
        let sides: Vec<Option<&&str>> = vec![
            vertex_side.get(&i0),
            vertex_side.get(&i1),
            vertex_side.get(&i2),
        ];

        let has_wall = sides.iter().any(|s| s.copied() == Some("WALL"));
        let has_cap = sides.iter().any(|s| s.copied() == Some("CAP"));

        if has_wall && has_cap {
            straddling_count += 1;
            if straddling_examples.len() < 10 {
                let v0 = mesh.vertices[i0];
                let v1 = mesh.vertices[i1];
                let v2 = mesh.vertices[i2];
                straddling_examples.push(format!(
                    "    tri[{ti}]: [{i0}({:?}),{i1}({:?}),{i2}({:?})] at z=[{:.2},{:.2},{:.2}] r=[{:.2},{:.2},{:.2}]",
                    sides[0].unwrap_or(&"?"), sides[1].unwrap_or(&"?"), sides[2].unwrap_or(&"?"),
                    v0.z, v1.z, v2.z,
                    (v0.x.powi(2)+v0.y.powi(2)).sqrt(),
                    (v1.x.powi(2)+v1.y.powi(2)).sqrt(),
                    (v2.x.powi(2)+v2.y.powi(2)).sqrt()
                ));
            }
        }
    }

    println!("  Straddling triangles (wall+cap vertices): {straddling_count}");
    for ex in &straddling_examples {
        println!("{ex}");
    }
    println!();

    // -----------------------------------------------------------------------
    // 7. Normal consistency check: angle between face normal and vertex normals
    // -----------------------------------------------------------------------
    println!("--- 7. FACE vs VERTEX NORMAL CONSISTENCY (rim region) ---");
    println!("  Large deviations indicate normals pointing wrong direction...");
    println!();

    let mut big_deviation_count = 0;
    let mut deviation_examples: Vec<String> = Vec::new();

    for (ti, tri) in mesh.indices.chunks(3).enumerate() {
        if tri.len() < 3 { continue; }
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];

        let centroid = (v0 + v1 + v2) / 3.0;
        let cz = centroid.z;
        let cr = (centroid.x.powi(2) + centroid.y.powi(2)).sqrt();

        // Only check triangles near the base rim
        if cz.abs() > 2.0 || (cr - 12.0).abs() > 3.0 { continue; }

        let face_n = (v1 - v0).cross(&(v2 - v0));
        let face_n_len = face_n.norm();
        if face_n_len < 1e-12 { continue; }
        let face_n = face_n / face_n_len;

        // Check each vertex normal against the face normal
        for &vi in tri {
            let vn = &mesh.normals[vi as usize];
            let dot = face_n.dot(vn);
            if dot < 0.0 {
                big_deviation_count += 1;
                if deviation_examples.len() < 5 {
                    deviation_examples.push(format!(
                        "    tri[{ti}] vert[{vi}]: dot(face,vert)={dot:.4}, face_n=({:.3},{:.3},{:.3}), vert_n=({:.3},{:.3},{:.3})",
                        face_n.x, face_n.y, face_n.z, vn.x, vn.y, vn.z
                    ));
                }
            }
        }
    }
    println!("  Rim vertex normals opposing face normal: {big_deviation_count}");
    for ex in &deviation_examples {
        println!("{ex}");
    }
    println!();

    // -----------------------------------------------------------------------
    // 8. Resolution cell size analysis
    // -----------------------------------------------------------------------
    println!("--- 8. RESOLUTION ANALYSIS ---");
    let bbox_size = bbox_max - bbox_min;
    // resolution 128 => depth = 7, so 2^7 = 128 cells per axis
    let depth = (128.0_f64).log2().ceil() as u32;
    let cells_per_axis = 2u32.pow(depth);
    let cell_size_x = bbox_size.x / cells_per_axis as f64;
    let cell_size_y = bbox_size.y / cells_per_axis as f64;
    let cell_size_z = bbox_size.z / cells_per_axis as f64;

    println!("  BBox: ({:.1},{:.1},{:.1}) to ({:.1},{:.1},{:.1})",
        bbox_min.x, bbox_min.y, bbox_min.z, bbox_max.x, bbox_max.y, bbox_max.z);
    println!("  BBox size: ({:.1}, {:.1}, {:.1})", bbox_size.x, bbox_size.y, bbox_size.z);
    println!("  Resolution: 128 => octree depth={depth}, cells_per_axis={cells_per_axis}");
    println!("  Cell size: ({cell_size_x:.4}, {cell_size_y:.4}, {cell_size_z:.4})");
    println!();

    // Base rim circumference = 2 * pi * 12 ~ 75.4
    // With cell size ~0.22, we get ~75.4/0.22 ~ 343 cells around the rim
    let circumference = 2.0 * std::f64::consts::PI * 12.0;
    let min_cell = cell_size_x.min(cell_size_y).min(cell_size_z);
    let cells_around_rim = circumference / min_cell;
    println!("  Base rim circumference: {circumference:.1}");
    println!("  Estimated cells around rim: {cells_around_rim:.0}");
    println!("  Angular resolution: {:.2} degrees per cell", 360.0 / cells_around_rim);
    println!();

    // -----------------------------------------------------------------------
    // 9. Gradient / normal smoothness at rim via direct SDF evaluation
    // -----------------------------------------------------------------------
    println!("--- 9. ANALYTICAL GRADIENT BEHAVIOR AT BASE RIM ---");
    println!("  Testing gradient continuity along the rim circle...");
    println!();

    // Sample points around the base rim at various offsets
    let test_offsets = [
        ("On rim (z=0, r=12)", 12.0, 0.0),
        ("Just above (z=0.1, r=12)", 12.0, 0.1),
        ("Just below (z=-0.1, r=12)", 12.0, -0.1),
        ("Just inside (z=0, r=11.9)", 11.9, 0.0),
        ("Just outside (z=0, r=12.1)", 12.1, 0.0),
        ("Wall side (z=0.5, r=11.8)", 11.8, 0.5),
        ("Cap side (z=0, r=11)", 11.0, 0.0),
    ];

    for (label, r, z) in &test_offsets {
        let p = Vector3::new(*r, 0.0, *z);
        let grad = cone.gradient(p).unwrap();
        let sdf_val = cone.evaluate(p);
        println!("  {label}:");
        println!("    SDF={sdf_val:.6}, gradient=({:.4}, {:.4}, {:.4}), |g|={:.6}",
            grad.x, grad.y, grad.z, grad.norm());
    }
    println!();

    // Check for gradient discontinuity by sampling a fine sweep across the rim
    println!("  Gradient sweep across rim edge (r=12, z from -0.5 to +0.5):");
    println!("  {:>6} {:>8} {:>8} {:>8} {:>8} {:>8}", "z", "sdf", "gx", "gy", "gz", "|g|");
    for i in -10..=10 {
        let z = i as f64 * 0.05;
        let p = Vector3::new(12.0, 0.0, z);
        let grad = cone.gradient(p).unwrap();
        let sdf_val = cone.evaluate(p);
        println!("  {:>6.2} {:>8.5} {:>8.4} {:>8.4} {:>8.4} {:>8.5}",
            z, sdf_val, grad.x, grad.y, grad.z, grad.norm());
    }
    println!();

    // -----------------------------------------------------------------------
    // 10. Summary / Diagnosis
    // -----------------------------------------------------------------------
    println!("========================================================");
    println!("  DIAGNOSIS SUMMARY");
    println!("========================================================");
    println!();
    println!("  (a) Geometry issues:");
    println!("      Degenerate tris (area<1e-8): {degenerate_count}");
    println!("      Very thin tris (aspect>100): {very_thin_count}");
    println!("      Straddling tris (wall+cap):  {straddling_count}");
    println!();
    println!("  (b) Normal issues:");
    println!("      Mixed/ambiguous rim normals:  {ambiguous_count}");
    println!("      Normals opposing face normal: {big_deviation_count}");
    println!();
    println!("  (c) Resolution:");
    println!("      Cell size:     {min_cell:.4}");
    println!("      Cells at rim:  {cells_around_rim:.0}");
    println!("      Total vertices: {n_verts}");
    println!("      Total tris:     {n_tris}");
    println!();

    if ambiguous_count > 0 {
        println!("  >>> PRIMARY SUSPECT: Per-vertex normals at the sharp rim edge are");
        println!("  >>> averaging the wall normal and cap normal, creating a 'smoothed'");
        println!("  >>> appearance where there should be a sharp crease. DC produces");
        println!("  >>> one vertex per cell, so rim vertices serve BOTH the wall and");
        println!("  >>> cap faces — but get only ONE normal. This is the 'chewy' look.");
    }
    if very_thin_count > 0 {
        println!("  >>> CONTRIBUTING FACTOR: {very_thin_count} very thin triangles may cause");
        println!("  >>> shading artifacts (faceting).");
    }
    if straddling_count > 0 {
        println!("  >>> CONTRIBUTING FACTOR: {straddling_count} triangles span the wall/cap");
        println!("  >>> boundary, which creates visible creasing artifacts.");
    }
    println!();
    println!("  End of diagnostic.");
}
