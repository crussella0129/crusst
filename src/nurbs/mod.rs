//! NURBS (Non-Uniform Rational B-Spline) curves and surfaces.
//!
//! This module provides the mathematical foundation for the B-Rep kernel:
//! basis function evaluation, knot vector operations, and rational
//! curve/surface evaluation with derivatives.

pub mod basis;
pub mod knot;

use crate::math::{Point2, Point3, Vector3};

/// A Non-Uniform Rational B-Spline curve in 3D.
#[derive(Clone, Debug)]
pub struct NurbsCurve3 {
    pub degree: usize,
    pub knots: Vec<f64>,
    pub control_points: Vec<Point3>,
    pub weights: Vec<f64>,
}

impl NurbsCurve3 {
    pub fn new(
        degree: usize,
        knots: Vec<f64>,
        control_points: Vec<Point3>,
        weights: Vec<f64>,
    ) -> Self {
        assert_eq!(control_points.len(), weights.len());
        assert!(
            knot::validate_knot_vector(&knots, degree, control_points.len()),
            "Invalid knot vector: {} knots for degree {} with {} control points",
            knots.len(),
            degree,
            control_points.len()
        );
        Self {
            degree,
            knots,
            control_points,
            weights,
        }
    }

    /// Evaluate the curve at parameter `u`.
    pub fn evaluate(&self, u: f64) -> Point3 {
        let n = self.control_points.len() - 1;
        let span = basis::find_span(n, self.degree, u, &self.knots);
        let b = basis::basis_funs(span, u, self.degree, &self.knots);

        let mut point = Vector3::zeros();
        let mut w = 0.0;

        for i in 0..=self.degree {
            let idx = span - self.degree + i;
            let wi = self.weights[idx];
            let bi = b[i] * wi;
            point += bi * self.control_points[idx].coords;
            w += bi;
        }

        Point3::from(point / w)
    }

    /// Evaluate the curve and its first derivative at parameter `u`.
    ///
    /// Uses the quotient rule on the rational form:
    /// C(u) = A(u) / w(u) where A(u) = Σ N_i·w_i·P_i and w(u) = Σ N_i·w_i
    /// C'(u) = (A'(u) - w'(u)·C(u)) / w(u)
    pub fn derivative(&self, u: f64) -> (Point3, Vector3) {
        let n = self.control_points.len() - 1;
        let span = basis::find_span(n, self.degree, u, &self.knots);
        let ders = basis::ders_basis_funs(span, u, self.degree, 1, &self.knots);

        let mut sw = 0.0;
        let mut dsw = 0.0;
        let mut sc = Vector3::zeros();
        let mut dsc = Vector3::zeros();

        for i in 0..=self.degree {
            let idx = span - self.degree + i;
            let wi = self.weights[idx];
            let pi = self.control_points[idx].coords;

            sc += ders[0][i] * wi * pi;
            sw += ders[0][i] * wi;
            dsc += ders[1][i] * wi * pi;
            dsw += ders[1][i] * wi;
        }

        let c = sc / sw;
        let point = Point3::from(c);
        let deriv = (dsc - dsw * c) / sw;

        (point, deriv)
    }

    /// Unit tangent vector at parameter `u`.
    pub fn tangent(&self, u: f64) -> Vector3 {
        let (_, d) = self.derivative(u);
        let len = d.norm();
        if len > 1e-15 {
            d / len
        } else {
            Vector3::new(1.0, 0.0, 0.0)
        }
    }

    /// The parameter domain `[u_min, u_max]`.
    pub fn domain(&self) -> (f64, f64) {
        (
            self.knots[self.degree],
            self.knots[self.knots.len() - self.degree - 1],
        )
    }

    /// Number of control points.
    pub fn num_control_points(&self) -> usize {
        self.control_points.len()
    }
}

/// A Non-Uniform Rational B-Spline curve in 2D (for parametric-space trim curves).
#[derive(Clone, Debug)]
pub struct NurbsCurve2 {
    pub degree: usize,
    pub knots: Vec<f64>,
    pub control_points: Vec<Point2>,
    pub weights: Vec<f64>,
}

impl NurbsCurve2 {
    pub fn new(
        degree: usize,
        knots: Vec<f64>,
        control_points: Vec<Point2>,
        weights: Vec<f64>,
    ) -> Self {
        assert_eq!(control_points.len(), weights.len());
        assert!(
            knot::validate_knot_vector(&knots, degree, control_points.len()),
            "Invalid knot vector for 2D NURBS curve"
        );
        Self {
            degree,
            knots,
            control_points,
            weights,
        }
    }

    /// Evaluate the curve at parameter `u`.
    pub fn evaluate(&self, u: f64) -> Point2 {
        let n = self.control_points.len() - 1;
        let span = basis::find_span(n, self.degree, u, &self.knots);
        let b = basis::basis_funs(span, u, self.degree, &self.knots);

        let mut point = nalgebra::Vector2::zeros();
        let mut w = 0.0;

        for i in 0..=self.degree {
            let idx = span - self.degree + i;
            let wi = self.weights[idx];
            let bi = b[i] * wi;
            point += bi * self.control_points[idx].coords;
            w += bi;
        }

        Point2::from(point / w)
    }

    /// The parameter domain `[u_min, u_max]`.
    pub fn domain(&self) -> (f64, f64) {
        (
            self.knots[self.degree],
            self.knots[self.knots.len() - self.degree - 1],
        )
    }
}

/// A Non-Uniform Rational B-Spline surface.
#[derive(Clone, Debug)]
pub struct NurbsSurface {
    pub degree_u: usize,
    pub degree_v: usize,
    pub knots_u: Vec<f64>,
    pub knots_v: Vec<f64>,
    /// Control point grid, indexed `[u_index][v_index]`.
    pub control_points: Vec<Vec<Point3>>,
    /// Weight grid, same layout as control points.
    pub weights: Vec<Vec<f64>>,
}

impl NurbsSurface {
    pub fn new(
        degree_u: usize,
        degree_v: usize,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        control_points: Vec<Vec<Point3>>,
        weights: Vec<Vec<f64>>,
    ) -> Self {
        let nu = control_points.len();
        let nv = control_points[0].len();
        assert!(
            knot::validate_knot_vector(&knots_u, degree_u, nu),
            "Invalid U knot vector"
        );
        assert!(
            knot::validate_knot_vector(&knots_v, degree_v, nv),
            "Invalid V knot vector"
        );
        assert_eq!(weights.len(), nu);
        for (i, row) in weights.iter().enumerate() {
            assert_eq!(
                row.len(),
                nv,
                "Weight row {i} has {} elements, expected {nv}",
                row.len()
            );
        }
        Self {
            degree_u,
            degree_v,
            knots_u,
            knots_v,
            control_points,
            weights,
        }
    }

    /// Evaluate the surface at parameters `(u, v)`.
    pub fn evaluate(&self, u: f64, v: f64) -> Point3 {
        let nu = self.control_points.len() - 1;
        let nv = self.control_points[0].len() - 1;

        let span_u = basis::find_span(nu, self.degree_u, u, &self.knots_u);
        let span_v = basis::find_span(nv, self.degree_v, v, &self.knots_v);
        let bu = basis::basis_funs(span_u, u, self.degree_u, &self.knots_u);
        let bv = basis::basis_funs(span_v, v, self.degree_v, &self.knots_v);

        let mut point = Vector3::zeros();
        let mut w = 0.0;

        for i in 0..=self.degree_u {
            let ui = span_u - self.degree_u + i;
            for j in 0..=self.degree_v {
                let vj = span_v - self.degree_v + j;
                let wij = self.weights[ui][vj];
                let bij = bu[i] * bv[j] * wij;
                point += bij * self.control_points[ui][vj].coords;
                w += bij;
            }
        }

        Point3::from(point / w)
    }

    /// Evaluate surface and partial derivatives at `(u, v)`.
    ///
    /// Returns `(point, dS/du, dS/dv)`. Uses the quotient rule on the
    /// rational bivariate form.
    pub fn derivatives(&self, u: f64, v: f64) -> (Point3, Vector3, Vector3) {
        let nu = self.control_points.len() - 1;
        let nv = self.control_points[0].len() - 1;

        let span_u = basis::find_span(nu, self.degree_u, u, &self.knots_u);
        let span_v = basis::find_span(nv, self.degree_v, v, &self.knots_v);
        let ders_u = basis::ders_basis_funs(span_u, u, self.degree_u, 1, &self.knots_u);
        let ders_v = basis::ders_basis_funs(span_v, v, self.degree_v, 1, &self.knots_v);

        let mut sw = 0.0;
        let mut sc = Vector3::zeros();
        let mut dsw_u = 0.0;
        let mut dsc_u = Vector3::zeros();
        let mut dsw_v = 0.0;
        let mut dsc_v = Vector3::zeros();

        for i in 0..=self.degree_u {
            let ui = span_u - self.degree_u + i;
            for j in 0..=self.degree_v {
                let vj = span_v - self.degree_v + j;
                let wij = self.weights[ui][vj];
                let pij = self.control_points[ui][vj].coords;

                let b00 = ders_u[0][i] * ders_v[0][j] * wij;
                let b10 = ders_u[1][i] * ders_v[0][j] * wij;
                let b01 = ders_u[0][i] * ders_v[1][j] * wij;

                sc += b00 * pij;
                sw += b00;
                dsc_u += b10 * pij;
                dsw_u += b10;
                dsc_v += b01 * pij;
                dsw_v += b01;
            }
        }

        let c = sc / sw;
        let point = Point3::from(c);
        let du = (dsc_u - dsw_u * c) / sw;
        let dv = (dsc_v - dsw_v * c) / sw;

        (point, du, dv)
    }

    /// Compute the surface normal at `(u, v)`.
    pub fn normal(&self, u: f64, v: f64) -> Vector3 {
        let (_, du, dv) = self.derivatives(u, v);
        let n = du.cross(&dv);
        let len = n.norm();
        if len > 1e-15 {
            n / len
        } else {
            Vector3::new(0.0, 0.0, 1.0)
        }
    }

    /// Parameter domain in U.
    pub fn domain_u(&self) -> (f64, f64) {
        (
            self.knots_u[self.degree_u],
            self.knots_u[self.knots_u.len() - self.degree_u - 1],
        )
    }

    /// Parameter domain in V.
    pub fn domain_v(&self) -> (f64, f64) {
        (
            self.knots_v[self.degree_v],
            self.knots_v[self.knots_v.len() - self.degree_v - 1],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Point3;

    /// Helper: create a degree-2 NURBS curve that is a straight line from (0,0,0) to (2,0,0).
    fn make_line_curve() -> NurbsCurve3 {
        NurbsCurve3::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
            ],
            vec![1.0, 1.0, 1.0],
        )
    }

    /// Helper: create a NURBS quarter-circle arc in XY plane.
    /// Uses the standard rational Bezier representation with weight 1/sqrt(2) on the middle point.
    fn make_quarter_circle(radius: f64) -> NurbsCurve3 {
        let w = std::f64::consts::FRAC_1_SQRT_2; // 1/√2
        NurbsCurve3::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![
                Point3::new(radius, 0.0, 0.0),
                Point3::new(radius, radius, 0.0),
                Point3::new(0.0, radius, 0.0),
            ],
            vec![1.0, w, 1.0],
        )
    }

    #[test]
    fn line_curve_endpoints() {
        let curve = make_line_curve();
        let p0 = curve.evaluate(0.0);
        let p1 = curve.evaluate(1.0);

        assert!((p0.x - 0.0).abs() < 1e-14);
        assert!((p0.y).abs() < 1e-14);
        assert!((p1.x - 2.0).abs() < 1e-14);
    }

    #[test]
    fn line_curve_midpoint() {
        let curve = make_line_curve();
        let pm = curve.evaluate(0.5);
        assert!((pm.x - 1.0).abs() < 1e-14);
    }

    #[test]
    fn quarter_circle_endpoints() {
        let r = 5.0;
        let curve = make_quarter_circle(r);
        let p0 = curve.evaluate(0.0);
        let p1 = curve.evaluate(1.0);

        assert!((p0.x - r).abs() < 1e-12, "Start should be (r, 0, 0)");
        assert!(p0.y.abs() < 1e-12);
        assert!(p1.x.abs() < 1e-12, "End should be (0, r, 0)");
        assert!((p1.y - r).abs() < 1e-12);
    }

    #[test]
    fn quarter_circle_on_circle() {
        // Every point on the NURBS quarter circle should have distance = radius from origin.
        let r = 5.0;
        let curve = make_quarter_circle(r);

        for i in 0..=20 {
            let u = i as f64 / 20.0;
            let p = curve.evaluate(u);
            let dist = (p.x * p.x + p.y * p.y).sqrt();
            assert!(
                (dist - r).abs() < 1e-12,
                "Point at u={u} has distance {dist} from origin, expected {r}"
            );
        }
    }

    #[test]
    fn quarter_circle_tangent_perpendicular() {
        // Tangent at any point on a circle is perpendicular to the radius vector.
        let r = 5.0;
        let curve = make_quarter_circle(r);

        for i in 1..20 {
            let u = i as f64 / 20.0;
            let (p, d) = curve.derivative(u);
            let radius_dir = p.coords;
            let dot = radius_dir.dot(&d);
            assert!(
                dot.abs() < 1e-10,
                "Tangent should be perpendicular to radius at u={u}: dot={dot}"
            );
        }
    }

    #[test]
    fn line_derivative_constant() {
        let curve = make_line_curve();
        // For a straight line from (0,0,0) to (2,0,0), the derivative should be constant.
        let (_, d0) = curve.derivative(0.25);
        let (_, d1) = curve.derivative(0.75);
        assert!(
            (d0 - d1).norm() < 1e-10,
            "Derivative of a line should be constant"
        );
    }

    #[test]
    fn derivative_finite_difference() {
        // Check that the analytical derivative matches a finite-difference approximation.
        let curve = make_quarter_circle(3.0);
        let u = 0.4;
        let h = 1e-7;
        let (_, analytic) = curve.derivative(u);
        let p_plus = curve.evaluate(u + h);
        let p_minus = curve.evaluate(u - h);
        let fd = (p_plus - p_minus) / (2.0 * h);

        assert!(
            (analytic - fd).norm() < 1e-5,
            "Analytical derivative should match finite difference: analytic={analytic:?}, fd={fd:?}"
        );
    }

    #[test]
    fn curve_domain() {
        let curve = make_line_curve();
        assert_eq!(curve.domain(), (0.0, 1.0));
    }

    // --- Surface tests ---

    /// Helper: create a bilinear NURBS surface patch (a flat quad).
    fn make_bilinear_surface() -> NurbsSurface {
        NurbsSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
                vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
            ],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        )
    }

    #[test]
    fn bilinear_surface_corners() {
        let surf = make_bilinear_surface();
        let p00 = surf.evaluate(0.0, 0.0);
        let p10 = surf.evaluate(1.0, 0.0);
        let p01 = surf.evaluate(0.0, 1.0);
        let p11 = surf.evaluate(1.0, 1.0);

        assert!((p00 - Point3::new(0.0, 0.0, 0.0)).norm() < 1e-14);
        assert!((p10 - Point3::new(1.0, 0.0, 0.0)).norm() < 1e-14);
        assert!((p01 - Point3::new(0.0, 1.0, 0.0)).norm() < 1e-14);
        assert!((p11 - Point3::new(1.0, 1.0, 0.0)).norm() < 1e-14);
    }

    #[test]
    fn bilinear_surface_center() {
        let surf = make_bilinear_surface();
        let pc = surf.evaluate(0.5, 0.5);
        assert!((pc - Point3::new(0.5, 0.5, 0.0)).norm() < 1e-14);
    }

    #[test]
    fn bilinear_surface_normal_is_z() {
        let surf = make_bilinear_surface();
        let n = surf.normal(0.5, 0.5);
        // A flat quad in XY should have normal pointing in +Z or -Z
        assert!(
            n.z.abs() > 0.99,
            "Normal of flat XY surface should be ±Z: got {n:?}"
        );
    }

    #[test]
    fn surface_derivatives_finite_difference() {
        let surf = make_bilinear_surface();
        let u = 0.3;
        let v = 0.7;
        let h = 1e-7;

        let (_, du, dv) = surf.derivatives(u, v);

        let fd_u = (surf.evaluate(u + h, v) - surf.evaluate(u - h, v)) / (2.0 * h);
        let fd_v = (surf.evaluate(u, v + h) - surf.evaluate(u, v - h)) / (2.0 * h);

        assert!(
            (du - fd_u).norm() < 1e-5,
            "dS/du analytical vs FD: {du:?} vs {fd_u:?}"
        );
        assert!(
            (dv - fd_v).norm() < 1e-5,
            "dS/dv analytical vs FD: {dv:?} vs {fd_v:?}"
        );
    }
}
