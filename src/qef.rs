use nalgebra::{Matrix3, Vector3};
use crate::types::BBox3;

/// Solve the Quadratic Error Function: find the point v that minimizes
/// Σ (nᵢ · (v - pᵢ))²
/// where pᵢ are edge crossing positions and nᵢ are surface normals.
///
/// Uses nalgebra's SVD to solve the 3×3 least-squares system.
/// Result is clamped to the cell bounds to prevent vertex drift.
pub fn solve_qef(
    positions: &[Vector3<f64>],
    normals: &[Vector3<f64>],
    bounds: &BBox3,
) -> Vector3<f64> {
    assert_eq!(positions.len(), normals.len());

    let n = positions.len();
    if n == 0 {
        return bounds.center();
    }

    // Compute mass point (centroid of crossing positions) for regularization.
    let mass_point: Vector3<f64> = positions.iter().sum::<Vector3<f64>>() / n as f64;

    // Build the normal equations: ATA * v = ATb
    // Each constraint row is nᵢᵀ, and each RHS element is nᵢ · pᵢ.
    // ATA = Σ nᵢ * nᵢᵀ   (3×3 outer product sum)
    // ATb = Σ nᵢ * (nᵢ · pᵢ)
    let mut ata = Matrix3::<f64>::zeros();
    let mut atb = Vector3::<f64>::zeros();

    for i in 0..n {
        let ni = &normals[i];
        let pi = &positions[i];
        let d = ni.dot(pi);

        // Accumulate outer product nᵢ * nᵢᵀ
        ata += ni * ni.transpose();
        // Accumulate nᵢ * (nᵢ · pᵢ)
        atb += ni * d;
    }

    // Mass point regularization to prevent degeneracy when normals are
    // parallel or nearly so.
    let regularization_weight = 1e-3;
    ata += Matrix3::identity() * regularization_weight;
    atb += mass_point * regularization_weight;

    // Solve via SVD (handles rank-deficient matrices gracefully).
    let svd = ata.svd(true, true);
    let v = svd.solve(&atb, 1e-10).unwrap_or(mass_point);

    // Clamp to cell bounds to prevent vertex drift.
    Vector3::new(
        v.x.clamp(bounds.min.x, bounds.max.x),
        v.y.clamp(bounds.min.y, bounds.max.y),
        v.z.clamp(bounds.min.z, bounds.max.z),
    )
}
