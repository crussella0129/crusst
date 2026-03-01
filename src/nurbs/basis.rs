//! Cox-de Boor B-spline basis function evaluation.
//!
//! Implements Algorithms A2.1–A2.3 from "The NURBS Book" (Piegl & Tiller).
//! These are the foundational routines for all NURBS curve and surface evaluation.

/// Find the knot span index such that `knots[span] <= u < knots[span+1]`.
///
/// `n` is the index of the last control point (= num_control_points - 1).
/// `p` is the degree. Uses binary search for O(log n) performance.
pub fn find_span(n: usize, p: usize, u: f64, knots: &[f64]) -> usize {
    // Clamp to domain boundaries
    if u >= knots[n + 1] {
        return n;
    }
    if u <= knots[p] {
        return p;
    }

    let mut lo = p;
    let mut hi = n + 1;
    let mut mid = (lo + hi) / 2;

    while u < knots[mid] || u >= knots[mid + 1] {
        if u < knots[mid] {
            hi = mid;
        } else {
            lo = mid;
        }
        mid = (lo + hi) / 2;
    }

    mid
}

/// Evaluate all non-zero basis functions at parameter `u`.
///
/// Returns a vector of `p+1` values: `N[span-p], N[span-p+1], ..., N[span]`.
/// Algorithm A2.2 from The NURBS Book.
pub fn basis_funs(span: usize, u: f64, p: usize, knots: &[f64]) -> Vec<f64> {
    let mut n = vec![0.0; p + 1];
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];

    n[0] = 1.0;

    for j in 1..=p {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;

        for r in 0..j {
            let temp = n[r] / (right[r + 1] + left[j - r]);
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }

    n
}

/// Evaluate basis functions and their derivatives up to order `n_ders`.
///
/// Returns `ders[k][j]` where `k` is the derivative order (0..=n_ders)
/// and `j` indexes the `p+1` non-zero basis functions.
/// Algorithm A2.3 from The NURBS Book.
pub fn ders_basis_funs(
    span: usize,
    u: f64,
    p: usize,
    n_ders: usize,
    knots: &[f64],
) -> Vec<Vec<f64>> {
    // ndu stores both the basis functions (upper triangle)
    // and the knot differences (lower triangle).
    let mut ndu = vec![vec![0.0; p + 1]; p + 1];
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];

    ndu[0][0] = 1.0;

    for j in 1..=p {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;

        for r in 0..j {
            // Lower triangle: knot differences
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = ndu[r][j - 1] / ndu[j][r];
            // Upper triangle: basis functions
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    // Load the basis functions (0th derivative)
    let n_out = n_ders.min(p);
    let mut ders = vec![vec![0.0; p + 1]; n_out + 1];
    for j in 0..=p {
        ders[0][j] = ndu[j][p];
    }

    // Compute derivatives using the recurrence relation.
    // a[s1][..] and a[s2][..] alternate as working rows.
    let mut a = vec![vec![0.0; p + 1]; 2];

    for r in 0..=p {
        let mut s1: usize = 0;
        let mut s2: usize = 1;
        a[0][0] = 1.0;

        for k in 1..=n_out {
            let mut d = 0.0;
            let rk = r as isize - k as isize;
            let pk = p as isize - k as isize;

            if rk >= 0 {
                a[s2][0] = a[s1][0] / ndu[pk as usize + 1][rk as usize];
                d = a[s2][0] * ndu[rk as usize][pk as usize];
            }

            let j1 = if rk >= -1 { 1usize } else { (-rk) as usize };
            let j2 = if (r as isize - 1) <= pk {
                k - 1
            } else {
                p - r
            };

            for j in j1..=j2 {
                a[s2][j] =
                    (a[s1][j] - a[s1][j - 1]) / ndu[pk as usize + 1][(rk + j as isize) as usize];
                d += a[s2][j] * ndu[(rk + j as isize) as usize][pk as usize];
            }

            if r <= pk as usize {
                a[s2][k] = -a[s1][k - 1] / ndu[pk as usize + 1][r];
                d += a[s2][k] * ndu[r][pk as usize];
            }

            ders[k][r] = d;
            std::mem::swap(&mut s1, &mut s2);
        }
    }

    // Multiply by correct factorial prefactors: p! / (p-k)!
    let mut factor = p as f64;
    for k in 1..=n_out {
        for j in 0..=p {
            ders[k][j] *= factor;
        }
        factor *= (p - k) as f64;
    }

    ders
}

#[cfg(test)]
mod tests {
    use super::*;

    // Clamped cubic knot vector: [0,0,0,0, 1,2,3, 4,4,4,4]
    // 7 control points, degree 3, 11 knots
    fn cubic_knots() -> Vec<f64> {
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.0]
    }

    #[test]
    fn find_span_interior() {
        let knots = cubic_knots();
        // n=6 (7 control points), p=3
        assert_eq!(find_span(6, 3, 0.0, &knots), 3);
        assert_eq!(find_span(6, 3, 0.5, &knots), 3);
        assert_eq!(find_span(6, 3, 1.0, &knots), 4);
        assert_eq!(find_span(6, 3, 2.5, &knots), 5);
        assert_eq!(find_span(6, 3, 4.0, &knots), 6); // at end
    }

    #[test]
    fn basis_partition_of_unity() {
        // B-spline basis functions must sum to 1 at any parameter value.
        let knots = cubic_knots();
        let n = 6;
        let p = 3;
        for &u in &[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0] {
            let span = find_span(n, p, u, &knots);
            let basis = basis_funs(span, u, p, &knots);
            let sum: f64 = basis.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-14,
                "Partition of unity failed at u={u}: sum={sum}"
            );
        }
    }

    #[test]
    fn basis_non_negative() {
        let knots = cubic_knots();
        let n = 6;
        let p = 3;
        for &u in &[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0] {
            let span = find_span(n, p, u, &knots);
            let basis = basis_funs(span, u, p, &knots);
            for (i, &val) in basis.iter().enumerate() {
                assert!(
                    val >= -1e-15,
                    "Negative basis function at u={u}, i={i}: {val}"
                );
            }
        }
    }

    #[test]
    fn ders_zeroth_matches_basis() {
        let knots = cubic_knots();
        let n = 6;
        let p = 3;
        let u = 1.5;
        let span = find_span(n, p, u, &knots);
        let basis = basis_funs(span, u, p, &knots);
        let ders = ders_basis_funs(span, u, p, 2, &knots);

        for j in 0..=p {
            assert!(
                (ders[0][j] - basis[j]).abs() < 1e-14,
                "Zeroth derivative should match basis function at j={j}"
            );
        }
    }

    #[test]
    fn ders_first_derivative_sum_zero() {
        // The sum of first derivatives of all non-zero basis functions is 0,
        // because the sum of basis functions is constant (1).
        let knots = cubic_knots();
        let n = 6;
        let p = 3;
        for &u in &[0.5, 1.5, 2.5, 3.5] {
            let span = find_span(n, p, u, &knots);
            let ders = ders_basis_funs(span, u, p, 1, &knots);
            let sum: f64 = ders[1].iter().sum();
            assert!(
                sum.abs() < 1e-12,
                "Sum of first derivatives should be 0 at u={u}: got {sum}"
            );
        }
    }

    // Simple linear basis (degree 1) — easy to verify by hand
    #[test]
    fn linear_basis_interpolation() {
        // Knots: [0, 0, 1, 2, 3, 3], degree=1, 4 control points
        let knots = vec![0.0, 0.0, 1.0, 2.0, 3.0, 3.0];
        let n = 3;
        let p = 1;

        // At u=0 (clamped start): N_0(0)=1, N_1(0)=0
        let span = find_span(n, p, 0.0, &knots);
        let b = basis_funs(span, 0.0, p, &knots);
        assert!((b[0] - 1.0).abs() < 1e-14);
        assert!((b[1] - 0.0).abs() < 1e-14);

        // At u=0.5: N[0]=0.5, N[1]=0.5 (linear interpolation)
        let span = find_span(n, p, 0.5, &knots);
        let b = basis_funs(span, 0.5, p, &knots);
        assert!((b[0] - 0.5).abs() < 1e-14);
        assert!((b[1] - 0.5).abs() < 1e-14);
    }
}
