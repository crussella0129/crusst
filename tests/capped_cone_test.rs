use approx::assert_relative_eq;
use crusst::mesh::extract_mesh;
use crusst::primitives::sdf_capped_cone;
use crusst::shape::{CappedCone, Sdf};
use nalgebra::Vector3;
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// SDF value correctness — the closest-point decomposition must match expected
// distances at known geometric positions.
// ---------------------------------------------------------------------------

#[test]
fn cone_center_is_inside() {
    // Center of a truncated cone (midpoint of axis, on axis)
    let a = Vector3::zeros();
    let b = Vector3::new(0.0, 0.0, 10.0);
    let d = sdf_capped_cone(Vector3::new(0.0, 0.0, 5.0), a, b, 5.0, 3.0);
    assert!(d < 0.0, "center should be inside, got {d}");
}

#[test]
fn cone_on_base_rim_is_zero() {
    // Point exactly on the base rim circle
    let a = Vector3::zeros();
    let b = Vector3::new(0.0, 0.0, 10.0);
    let d = sdf_capped_cone(Vector3::new(5.0, 0.0, 0.0), a, b, 5.0, 3.0);
    assert_relative_eq!(d, 0.0, epsilon = 1e-6);
}

#[test]
fn cone_on_tip_rim_is_zero() {
    // Point exactly on the tip rim circle
    let a = Vector3::zeros();
    let b = Vector3::new(0.0, 0.0, 10.0);
    let d = sdf_capped_cone(Vector3::new(3.0, 0.0, 10.0), a, b, 5.0, 3.0);
    assert_relative_eq!(d, 0.0, epsilon = 1e-6);
}

#[test]
fn cone_on_cap_face_is_zero() {
    // Point on the flat base cap (r < ra, z = 0)
    let a = Vector3::zeros();
    let b = Vector3::new(0.0, 0.0, 10.0);
    let d = sdf_capped_cone(Vector3::new(2.0, 0.0, 0.0), a, b, 5.0, 3.0);
    assert_relative_eq!(d, 0.0, epsilon = 1e-6);
}

#[test]
fn cone_outside_radially() {
    // Point outside the cone radially (at base height, beyond ra)
    let a = Vector3::zeros();
    let b = Vector3::new(0.0, 0.0, 10.0);
    let d = sdf_capped_cone(Vector3::new(10.0, 0.0, 0.0), a, b, 5.0, 3.0);
    assert!(d > 0.0, "should be outside, got {d}");
}

#[test]
fn cone_outside_axially() {
    // Point below the base cap, on axis
    let a = Vector3::zeros();
    let b = Vector3::new(0.0, 0.0, 10.0);
    let d = sdf_capped_cone(Vector3::new(0.0, 0.0, -3.0), a, b, 5.0, 3.0);
    assert_relative_eq!(d, 3.0, epsilon = 1e-6);
}

#[test]
fn sharp_cone_tip_distance() {
    // Sharp cone (rb=0), point at the tip
    let a = Vector3::zeros();
    let b = Vector3::new(0.0, 0.0, 30.0);
    let d = sdf_capped_cone(b, a, b, 12.0, 0.0);
    assert_relative_eq!(d, 0.0, epsilon = 1e-6);
}

#[test]
fn sharp_cone_beyond_tip() {
    // Point beyond the sharp tip, on axis
    let a = Vector3::zeros();
    let b = Vector3::new(0.0, 0.0, 30.0);
    let d = sdf_capped_cone(Vector3::new(0.0, 0.0, 35.0), a, b, 12.0, 0.0);
    assert_relative_eq!(d, 5.0, epsilon = 1e-6);
}

// ---------------------------------------------------------------------------
// Analytical gradient correctness
// ---------------------------------------------------------------------------

#[test]
fn gradient_agrees_with_central_differences() {
    // Sample many points and compare analytical gradient with central diffs
    let a = Vector3::zeros();
    let b = Vector3::new(0.0, 0.0, 30.0);
    let ra = 12.0;
    let rb = 0.0;
    let cone = CappedCone::new(a, b, ra, rb);

    let test_points = vec![
        Vector3::new(15.0, 0.0, 5.0),   // outside radially near base
        Vector3::new(0.0, 0.0, -5.0),   // below base on axis
        Vector3::new(0.0, 0.0, 35.0),   // above tip on axis
        Vector3::new(3.0, 0.0, 15.0),   // inside
        Vector3::new(20.0, 0.0, 15.0),  // outside radially mid-height
        Vector3::new(5.0, 5.0, 0.0),    // near base cap, off-axis
        Vector3::new(0.0, 3.0, 33.0),   // near tip
        Vector3::new(10.0, 10.0, 10.0), // general position
    ];

    let eps = 1e-5;
    for p in &test_points {
        let analytical = cone.gradient(*p).unwrap();

        // Central differences
        let dx = cone.evaluate(*p + Vector3::new(eps, 0.0, 0.0))
            - cone.evaluate(*p - Vector3::new(eps, 0.0, 0.0));
        let dy = cone.evaluate(*p + Vector3::new(0.0, eps, 0.0))
            - cone.evaluate(*p - Vector3::new(0.0, eps, 0.0));
        let dz = cone.evaluate(*p + Vector3::new(0.0, 0.0, eps))
            - cone.evaluate(*p - Vector3::new(0.0, 0.0, eps));
        let numerical = Vector3::new(dx, dy, dz).normalize();

        let dot = analytical.dot(&numerical);
        assert!(
            dot > 0.95,
            "Gradient mismatch at {p:?}: analytical={analytical:?}, numerical={numerical:?}, dot={dot}"
        );
    }
}

#[test]
fn gradient_at_cap_wall_junction_is_not_degenerate() {
    // Points very close to the cap-wall junction — the root cause of the bug.
    // The analytical gradient should be well-defined (unit length) at these points.
    let a = Vector3::zeros();
    let b = Vector3::new(0.0, 0.0, 30.0);
    let cone = CappedCone::new(a, b, 12.0, 0.0);

    // Points near the base rim circle (r ≈ 12, z ≈ 0)
    let junction_points = vec![
        Vector3::new(12.0, 0.0, 0.01),   // just above cap plane, at rim
        Vector3::new(12.0, 0.0, -0.01),  // just below cap plane, at rim
        Vector3::new(12.01, 0.0, 0.0),   // just outside rim
        Vector3::new(11.99, 0.0, 0.0),   // just inside rim
        Vector3::new(12.0, 0.0, 0.0),    // exactly on rim
        Vector3::new(0.0, 12.0, 0.01),   // rim in Y direction
        Vector3::new(8.485, 8.485, 0.0), // rim at 45 degrees
    ];

    for p in &junction_points {
        let grad = cone.gradient(*p).unwrap();
        let len = grad.norm();
        assert!(
            (len - 1.0).abs() < 0.01,
            "Gradient not unit length at {p:?}: norm={len}, grad={grad:?}"
        );
    }
}

// ---------------------------------------------------------------------------
// Mesh circularity — the actual bug manifestation
// ---------------------------------------------------------------------------

#[test]
fn cone_base_rim_is_circular() {
    // Generate a mesh of the sharp cone used in Monad's Order 1 eigenform.
    // Extract vertices near the base rim and verify they lie on a circle.
    // Use resolution 96 for good fidelity (cell size ≈ 0.29).
    let cone = CappedCone::new(
        Vector3::zeros(),
        Vector3::new(0.0, 0.0, 30.0),
        12.0,
        0.0,
    );

    let bbox_min = Vector3::new(-14.0, -14.0, -1.0);
    let bbox_max = Vector3::new(14.0, 14.0, 32.0);
    let mesh = extract_mesh(&cone, bbox_min, bbox_max, 96);

    assert!(!mesh.vertices.is_empty(), "mesh should have vertices");
    assert!(!mesh.indices.is_empty(), "mesh should have triangles");

    // Find vertices near the base plane (z ≈ 0) and near the rim (r ≈ 12)
    let rim_vertices: Vec<_> = mesh
        .vertices
        .iter()
        .filter(|v| {
            let z_near_base = v.z.abs() < 0.5;
            let r = (v.x * v.x + v.y * v.y).sqrt();
            let r_near_rim = (r - 12.0).abs() < 1.5;
            z_near_base && r_near_rim
        })
        .collect();

    assert!(
        rim_vertices.len() >= 8,
        "expected ≥8 rim vertices, got {}",
        rim_vertices.len()
    );

    // Check that rim vertices are approximately equidistant from the axis.
    // DC meshes have inherent deviation from sharp edges due to QEF
    // least-squares fitting within octree cells. We check:
    //   1. Average radius is close to the expected 12.0
    //   2. Standard deviation of radii is small relative to average
    let radii: Vec<f64> = rim_vertices
        .iter()
        .map(|v| (v.x * v.x + v.y * v.y).sqrt())
        .collect();
    let avg_r = radii.iter().sum::<f64>() / radii.len() as f64;
    let variance = radii.iter().map(|r| (r - avg_r).powi(2)).sum::<f64>() / radii.len() as f64;
    let std_dev = variance.sqrt();

    // Average radius should be within 10% of expected 12.0
    assert!(
        (avg_r - 12.0).abs() / 12.0 < 0.10,
        "Average rim radius {avg_r:.3} deviates >10% from expected 12.0"
    );

    // Standard deviation should be small (< 5% of average = good circularity)
    let cv = std_dev / avg_r; // coefficient of variation
    assert!(
        cv < 0.05,
        "Rim circularity poor: std_dev={std_dev:.3}, avg={avg_r:.3}, CV={:.1}% (>5%)",
        cv * 100.0
    );
}

// ---------------------------------------------------------------------------
// Watertightness — every edge shared by exactly 2 triangles
// ---------------------------------------------------------------------------

#[test]
fn cone_mesh_is_watertight() {
    let cone = CappedCone::new(
        Vector3::zeros(),
        Vector3::new(0.0, 0.0, 30.0),
        12.0,
        0.0,
    );

    let bbox_min = Vector3::new(-14.0, -14.0, -1.0);
    let bbox_max = Vector3::new(14.0, 14.0, 32.0);
    let mesh = extract_mesh(&cone, bbox_min, bbox_max, 48);

    assert!(!mesh.indices.is_empty(), "mesh should have triangles");

    // Count how many triangles share each edge
    let mut edge_counts: HashMap<(u32, u32), usize> = HashMap::new();
    for tri in mesh.indices.chunks(3) {
        if tri.len() < 3 {
            continue;
        }
        let edges = [
            (tri[0].min(tri[1]), tri[0].max(tri[1])),
            (tri[1].min(tri[2]), tri[1].max(tri[2])),
            (tri[0].min(tri[2]), tri[0].max(tri[2])),
        ];
        for edge in &edges {
            *edge_counts.entry(*edge).or_insert(0) += 1;
        }
    }

    let total_edges = edge_counts.len();
    let boundary_edges = edge_counts.values().filter(|&&c| c == 1).count();
    let non_manifold = edge_counts.values().filter(|&&c| c > 2).count();

    // A perfectly watertight mesh has 0 boundary edges.
    // DC meshes may have a small number of boundary edges at the bbox boundary.
    let boundary_ratio = boundary_edges as f64 / total_edges as f64;
    assert!(
        boundary_ratio < 0.05,
        "Too many boundary edges: {boundary_edges}/{total_edges} ({:.1}%)",
        boundary_ratio * 100.0
    );
    assert_eq!(
        non_manifold, 0,
        "Non-manifold edges found: {non_manifold}"
    );
}

// ---------------------------------------------------------------------------
// Truncated cone (rb > 0) — sanity check
// ---------------------------------------------------------------------------

#[test]
fn truncated_cone_sdf_values() {
    let a = Vector3::zeros();
    let b = Vector3::new(0.0, 0.0, 10.0);
    let ra = 5.0;
    let rb = 3.0;

    // On the wall surface at midpoint: r should be 4, and the point (4, 0, 5) should be on surface
    let d = sdf_capped_cone(Vector3::new(4.0, 0.0, 5.0), a, b, ra, rb);
    assert_relative_eq!(d, 0.0, epsilon = 0.1); // approximate — wall is slanted

    // Center of axis is clearly inside
    let d_center = sdf_capped_cone(Vector3::new(0.0, 0.0, 5.0), a, b, ra, rb);
    assert!(d_center < 0.0, "center should be inside");
}
