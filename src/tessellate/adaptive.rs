//! Adaptive parametric tessellation.
//!
//! Subdivides a surface's (u,v) domain into triangles, refining where the
//! chord deviation exceeds the tolerance or edge lengths exceed the maximum.

use crate::surface::Surface;
use nalgebra::Vector3;

/// Adaptively tessellate a surface patch.
///
/// Returns (vertices, normals, triangle_indices).
pub fn adaptive_tessellate(
    surface: &Surface,
    u_range: (f64, f64),
    v_range: (f64, f64),
    base_n: usize,
    chord_tol: f64,
    max_edge_len: f64,
    outward: bool,
) -> (Vec<Vector3<f64>>, Vec<Vector3<f64>>, Vec<u32>) {
    // Phase 1: Build initial grid
    let nu = base_n.max(2);
    let nv = base_n.max(2);

    let mut params: Vec<(f64, f64)> = Vec::with_capacity((nu + 1) * (nv + 1));
    let mut verts: Vec<Vector3<f64>> = Vec::with_capacity((nu + 1) * (nv + 1));
    let mut normals: Vec<Vector3<f64>> = Vec::with_capacity((nu + 1) * (nv + 1));

    for iv in 0..=nv {
        let v = v_range.0 + (v_range.1 - v_range.0) * iv as f64 / nv as f64;
        for iu in 0..=nu {
            let u = u_range.0 + (u_range.1 - u_range.0) * iu as f64 / nu as f64;
            let pt = surface.evaluate(u, v);
            let mut n = surface.normal(u, v);
            if !outward {
                n = -n;
            }
            params.push((u, v));
            verts.push(Vector3::new(pt.x, pt.y, pt.z));
            normals.push(n);
        }
    }

    // Phase 2: Build initial triangulation from grid
    let mut indices: Vec<u32> = Vec::new();
    for iv in 0..nv {
        for iu in 0..nu {
            let i00 = (iv * (nu + 1) + iu) as u32;
            let i10 = i00 + 1;
            let i01 = i00 + (nu + 1) as u32;
            let i11 = i01 + 1;

            if outward {
                indices.extend_from_slice(&[i00, i10, i11]);
                indices.extend_from_slice(&[i00, i11, i01]);
            } else {
                indices.extend_from_slice(&[i00, i11, i10]);
                indices.extend_from_slice(&[i00, i01, i11]);
            }
        }
    }

    // Phase 3: Adaptive refinement — subdivide triangles with excessive chord deviation
    // We do a fixed number of refinement passes.
    let max_passes = 3;
    for _ in 0..max_passes {
        let mut new_indices: Vec<u32> = Vec::new();
        let mut any_refined = false;

        let mut tri_idx = 0;
        while tri_idx < indices.len() {
            let i0 = indices[tri_idx] as usize;
            let i1 = indices[tri_idx + 1] as usize;
            let i2 = indices[tri_idx + 2] as usize;
            tri_idx += 3;

            let needs_refine = should_refine(
                &verts, &params, i0, i1, i2, surface, chord_tol, max_edge_len,
            );

            if needs_refine {
                any_refined = true;
                // Subdivide by inserting midpoints on all 3 edges
                let m01 = get_or_add_midpoint(
                    &mut verts,
                    &mut normals,
                    &mut params,
                    i0, i1,
                    surface,
                    outward,
                );
                let m12 = get_or_add_midpoint(
                    &mut verts,
                    &mut normals,
                    &mut params,
                    i1, i2,
                    surface,
                    outward,
                );
                let m20 = get_or_add_midpoint(
                    &mut verts,
                    &mut normals,
                    &mut params,
                    i2, i0,
                    surface,
                    outward,
                );

                // 4 sub-triangles
                new_indices.extend_from_slice(&[i0 as u32, m01, m20]);
                new_indices.extend_from_slice(&[m01, i1 as u32, m12]);
                new_indices.extend_from_slice(&[m20, m12, i2 as u32]);
                new_indices.extend_from_slice(&[m01, m12, m20]);
            } else {
                new_indices.extend_from_slice(&[i0 as u32, i1 as u32, i2 as u32]);
            }
        }

        indices = new_indices;
        if !any_refined {
            break;
        }
    }

    (verts, normals, indices)
}

/// Check if a triangle needs refinement based on chord deviation and edge length.
fn should_refine(
    verts: &[Vector3<f64>],
    params: &[(f64, f64)],
    i0: usize,
    i1: usize,
    i2: usize,
    surface: &Surface,
    chord_tol: f64,
    max_edge_len: f64,
) -> bool {
    let edges = [(i0, i1), (i1, i2), (i2, i0)];

    for &(ia, ib) in &edges {
        // Check edge length
        let edge_len = (verts[ia] - verts[ib]).norm();
        if edge_len > max_edge_len {
            return true;
        }

        // Check chord deviation: compare midpoint on surface vs linear midpoint
        let (ua, va) = params[ia];
        let (ub, vb) = params[ib];
        let u_mid = (ua + ub) * 0.5;
        let v_mid = (va + vb) * 0.5;

        let surface_mid = surface.evaluate(u_mid, v_mid);
        let linear_mid = (verts[ia] + verts[ib]) * 0.5;

        let deviation = (Vector3::new(surface_mid.x, surface_mid.y, surface_mid.z) - linear_mid).norm();
        if deviation > chord_tol {
            return true;
        }
    }

    false
}

/// Get or add the midpoint vertex between two existing vertices.
/// Uses a simple approach: always adds a new vertex (no dedup for simplicity).
fn get_or_add_midpoint(
    verts: &mut Vec<Vector3<f64>>,
    normals: &mut Vec<Vector3<f64>>,
    params: &mut Vec<(f64, f64)>,
    i0: usize,
    i1: usize,
    surface: &Surface,
    outward: bool,
) -> u32 {
    let (u0, v0) = params[i0];
    let (u1, v1) = params[i1];
    let u_mid = (u0 + u1) * 0.5;
    let v_mid = (v0 + v1) * 0.5;

    let pt = surface.evaluate(u_mid, v_mid);
    let mut n = surface.normal(u_mid, v_mid);
    if !outward {
        n = -n;
    }

    let idx = verts.len() as u32;
    verts.push(Vector3::new(pt.x, pt.y, pt.z));
    normals.push(n);
    params.push((u_mid, v_mid));
    idx
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Point3;

    #[test]
    fn plane_tessellation_flat() {
        let plane = Surface::Plane {
            origin: Point3::origin(),
            normal: crate::math::Vector3::new(0.0, 0.0, 1.0),
        };

        let (verts, normals, indices) = adaptive_tessellate(
            &plane,
            (0.0, 1.0),
            (0.0, 1.0),
            4,
            0.01,
            10.0,
            true,
        );

        assert!(!verts.is_empty());
        assert_eq!(verts.len(), normals.len());
        assert!(indices.len() % 3 == 0, "Indices must be multiple of 3");

        // All normals should point in +Z for a flat plane
        for n in &normals {
            assert!(n.z > 0.99, "Plane normal should be +Z: got {:?}", n);
        }

        // All Z coordinates should be 0 for a Z=0 plane
        for v in &verts {
            assert!(v.z.abs() < 1e-12, "Plane vertex z should be 0");
        }
    }

    #[test]
    fn sphere_tessellation_on_sphere() {
        let sphere = Surface::Sphere {
            center: Point3::origin(),
            radius: 5.0,
        };

        let (verts, _, indices) = adaptive_tessellate(
            &sphere,
            (0.0, std::f64::consts::TAU),
            (-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2),
            8,
            0.1,
            5.0,
            true,
        );

        assert!(!verts.is_empty());
        assert!(indices.len() > 100, "Sphere should have many triangles");

        // All vertices should be on the sphere (distance = radius from origin)
        for v in &verts {
            let r = v.norm();
            assert!(
                (r - 5.0).abs() < 1e-10,
                "Sphere vertex at distance {r}, expected 5.0"
            );
        }
    }

    #[test]
    fn cylinder_tessellation_on_cylinder() {
        let cyl = Surface::Cylinder {
            origin: Point3::origin(),
            axis: crate::math::Vector3::new(0.0, 0.0, 1.0),
            radius: 3.0,
        };

        let (verts, normals, _) = adaptive_tessellate(
            &cyl,
            (0.0, std::f64::consts::TAU),
            (0.0, 10.0),
            8,
            0.01,
            5.0,
            true,
        );

        // All vertices should be on the cylinder (xy distance = radius)
        for v in &verts {
            let r = (v.x * v.x + v.y * v.y).sqrt();
            assert!(
                (r - 3.0).abs() < 1e-10,
                "Cylinder vertex at xy-distance {r}, expected 3.0"
            );
        }

        // Normals should be radially outward (z component ≈ 0)
        for n in &normals {
            assert!(n.z.abs() < 0.01, "Cylinder normal z should be ~0: got {}", n.z);
        }
    }
}
