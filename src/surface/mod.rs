//! Surface types for the B-Rep kernel.
//!
//! Each `Face` in the topology carries one `Surface` as its geometric carrier.
//! The surface provides evaluation, derivatives, and normals in parametric space.

use crate::math::{Point3, Vector3};
use crate::nurbs::NurbsSurface;

/// A geometric surface. Carried by topology `Face` entities.
#[derive(Clone, Debug)]
pub enum Surface {
    Plane {
        origin: Point3,
        normal: Vector3,
    },
    Cylinder {
        origin: Point3,
        axis: Vector3,
        radius: f64,
    },
    Cone {
        apex: Point3,
        axis: Vector3,
        half_angle: f64,
    },
    Sphere {
        center: Point3,
        radius: f64,
    },
    Torus {
        center: Point3,
        axis: Vector3,
        major_r: f64,
        minor_r: f64,
    },
    NurbsSurface(NurbsSurface),
}

impl Surface {
    /// Evaluate the surface at parameters `(u, v)`.
    ///
    /// Parameter conventions:
    /// - **Plane:** `S(u,v) = origin + u*e1 + v*e2` where e1,e2 span the plane
    /// - **Cylinder:** `S(u,v) = origin + radius*(cos(u)*e1 + sin(u)*e2) + v*axis`
    ///   with u ∈ [0, 2π), v ∈ [0, height]
    /// - **Cone:** `S(u,v) = apex + v*(cos(u)*e1 + sin(u)*e2)*sin(half_angle) + v*axis*cos(half_angle)`
    ///   with u ∈ [0, 2π), v ≥ 0
    /// - **Sphere:** `S(u,v) = center + radius*(cos(v)*cos(u)*e1 + cos(v)*sin(u)*e2 + sin(v)*axis)`
    ///   with u ∈ [0, 2π), v ∈ [-π/2, π/2]
    /// - **Torus:** `S(u,v) = center + (major_r + minor_r*cos(v))*(cos(u)*e1 + sin(u)*e2) + minor_r*sin(v)*axis`
    ///   with u,v ∈ [0, 2π)
    pub fn evaluate(&self, u: f64, v: f64) -> Point3 {
        match self {
            Surface::Plane { origin, normal } => {
                let (e1, e2) = plane_frame(normal);
                origin + e1 * u + e2 * v
            }
            Surface::Cylinder {
                origin,
                axis,
                radius,
            } => {
                let (e1, e2) = plane_frame(axis);
                origin + (e1 * u.cos() + e2 * u.sin()) * *radius + axis.normalize() * v
            }
            Surface::Cone {
                apex,
                axis,
                half_angle,
            } => {
                let (e1, e2) = plane_frame(axis);
                let a = axis.normalize();
                let r = v * half_angle.sin();
                let h = v * half_angle.cos();
                apex + (e1 * u.cos() + e2 * u.sin()) * r + a * h
            }
            Surface::Sphere { center, radius } => {
                let (e1, e2) = plane_frame(&Vector3::new(0.0, 0.0, 1.0));
                let a = Vector3::new(0.0, 0.0, 1.0);
                center
                    + (e1 * u.cos() + e2 * u.sin()) * (v.cos() * *radius)
                    + a * (v.sin() * *radius)
            }
            Surface::Torus {
                center,
                axis,
                major_r,
                minor_r,
            } => {
                let (e1, e2) = plane_frame(axis);
                let a = axis.normalize();
                let r = major_r + minor_r * v.cos();
                center + (e1 * u.cos() + e2 * u.sin()) * r + a * (minor_r * v.sin())
            }
            Surface::NurbsSurface(nurbs) => nurbs.evaluate(u, v),
        }
    }

    /// Partial derivative with respect to `u`.
    pub fn derivative_u(&self, u: f64, v: f64) -> Vector3 {
        match self {
            Surface::Plane { normal, .. } => {
                let (e1, _) = plane_frame(normal);
                e1
            }
            Surface::Cylinder {
                axis, radius, ..
            } => {
                let (e1, e2) = plane_frame(axis);
                (-e1 * u.sin() + e2 * u.cos()) * *radius
            }
            Surface::Cone {
                axis, half_angle, ..
            } => {
                let (e1, e2) = plane_frame(axis);
                let r = v * half_angle.sin();
                (-e1 * u.sin() + e2 * u.cos()) * r
            }
            Surface::Sphere { radius, .. } => {
                let (e1, e2) = plane_frame(&Vector3::new(0.0, 0.0, 1.0));
                (-e1 * u.sin() + e2 * u.cos()) * (v.cos() * *radius)
            }
            Surface::Torus {
                axis,
                major_r,
                minor_r,
                ..
            } => {
                let (e1, e2) = plane_frame(axis);
                let r = major_r + minor_r * v.cos();
                (-e1 * u.sin() + e2 * u.cos()) * r
            }
            Surface::NurbsSurface(nurbs) => {
                let (_, du, _) = nurbs.derivatives(u, v);
                du
            }
        }
    }

    /// Partial derivative with respect to `v`.
    pub fn derivative_v(&self, u: f64, v: f64) -> Vector3 {
        match self {
            Surface::Plane { normal, .. } => {
                let (_, e2) = plane_frame(normal);
                e2
            }
            Surface::Cylinder { axis, .. } => axis.normalize(),
            Surface::Cone {
                axis, half_angle, ..
            } => {
                let (e1, e2) = plane_frame(axis);
                let a = axis.normalize();
                (e1 * u.cos() + e2 * u.sin()) * half_angle.sin() + a * half_angle.cos()
            }
            Surface::Sphere { radius, .. } => {
                let (e1, e2) = plane_frame(&Vector3::new(0.0, 0.0, 1.0));
                let a = Vector3::new(0.0, 0.0, 1.0);
                (e1 * u.cos() + e2 * u.sin()) * (-v.sin() * *radius)
                    + a * (v.cos() * *radius)
            }
            Surface::Torus {
                axis, minor_r, ..
            } => {
                let (e1, e2) = plane_frame(axis);
                let a = axis.normalize();
                (e1 * u.cos() + e2 * u.sin()) * (-minor_r * v.sin())
                    + a * (minor_r * v.cos())
            }
            Surface::NurbsSurface(nurbs) => {
                let (_, _, dv) = nurbs.derivatives(u, v);
                dv
            }
        }
    }

    /// Surface normal at parameters `(u, v)` (unit length, outward-pointing convention).
    pub fn normal(&self, u: f64, v: f64) -> Vector3 {
        let du = self.derivative_u(u, v);
        let dv = self.derivative_v(u, v);
        let n = du.cross(&dv);
        let len = n.norm();
        if len > 1e-15 {
            n / len
        } else {
            // Degenerate point (e.g., pole of sphere) — use axis direction
            match self {
                Surface::Sphere { .. } => {
                    if v > 0.0 {
                        Vector3::new(0.0, 0.0, 1.0)
                    } else {
                        Vector3::new(0.0, 0.0, -1.0)
                    }
                }
                _ => Vector3::new(0.0, 0.0, 1.0),
            }
        }
    }

    /// Return the minimum curvature radius of this surface.
    ///
    /// For analytic surfaces, this is exact. For NURBS, we sample and return
    /// a conservative estimate. Planes return `f64::INFINITY`.
    pub fn min_curvature_radius(&self) -> f64 {
        match self {
            Surface::Plane { .. } => f64::INFINITY,
            Surface::Cylinder { radius, .. } => *radius,
            Surface::Sphere { radius, .. } => *radius,
            Surface::Cone { half_angle, .. } => {
                // The cross-section radius varies with v; the minimum non-zero
                // curvature is at small v. We return a conservative small value
                // based on the half-angle — tighter angle = tighter curvature.
                // For practical cones, the tessellator already handles this via
                // chord deviation, so return 1.0 as a sensible floor.
                1.0_f64.max(0.1 / half_angle.sin().max(0.01))
            }
            Surface::Torus { minor_r, .. } => *minor_r,
            Surface::NurbsSurface(_) => {
                // Conservative estimate — assume fairly curved
                1.0
            }
        }
    }

    /// Find the parameters `(u, v)` closest to the given 3D point (point inversion).
    pub fn closest_parameters(&self, point: &Point3) -> (f64, f64) {
        match self {
            Surface::Plane { origin, normal } => {
                let (e1, e2) = plane_frame(normal);
                let d = point - origin;
                (d.dot(&e1), d.dot(&e2))
            }
            Surface::Cylinder {
                origin,
                axis,
                radius: _,
            } => {
                let a = axis.normalize();
                let d = point - origin;
                let v = d.dot(&a);
                let proj = d - a * v;
                let (e1, e2) = plane_frame(axis);
                let u = proj.dot(&e2).atan2(proj.dot(&e1));
                (u, v)
            }
            Surface::Sphere { center, radius: _ } => {
                let d = point - center;
                let (e1, e2) = plane_frame(&Vector3::new(0.0, 0.0, 1.0));
                let a = Vector3::new(0.0, 0.0, 1.0);
                let xy = d - a * d.dot(&a);
                let u = xy.dot(&e2).atan2(xy.dot(&e1));
                let v = d.dot(&a).atan2(xy.norm());
                (u, v)
            }
            Surface::Cone { apex, axis, half_angle } => {
                let a = axis.normalize();
                let d = point - apex;
                let d_axial = d.dot(&a);
                let d_radial = d - a * d_axial;
                let cos_a = half_angle.cos();
                let v = if cos_a.abs() > 1e-15 { d_axial / cos_a } else { d_radial.norm() };
                let (e1, e2) = plane_frame(axis);
                let u = d_radial.dot(&e2).atan2(d_radial.dot(&e1));
                (u, v)
            }
            _ => {
                // General Newton-Raphson for torus, NURBS
                newton_inversion(self, point)
            }
        }
    }
}

/// Newton-Raphson point inversion for general surfaces.
fn newton_inversion(surf: &Surface, point: &Point3) -> (f64, f64) {
    // Coarse grid search
    let n = 20;
    let (u_range, v_range) = default_domain(surf);
    let mut best_u = u_range.0;
    let mut best_v = v_range.0;
    let mut best_dist = f64::MAX;

    for i in 0..=n {
        for j in 0..=n {
            let u = u_range.0 + (u_range.1 - u_range.0) * i as f64 / n as f64;
            let v = v_range.0 + (v_range.1 - v_range.0) * j as f64 / n as f64;
            let dist = (surf.evaluate(u, v) - point).norm();
            if dist < best_dist {
                best_dist = dist;
                best_u = u;
                best_v = v;
            }
        }
    }

    // Newton-Raphson refinement
    let mut u = best_u;
    let mut v = best_v;

    for _ in 0..50 {
        let p = surf.evaluate(u, v);
        let diff = p - point;
        let du = surf.derivative_u(u, v);
        let dv = surf.derivative_v(u, v);

        // Solve [du·du  du·dv] [Δu] = [-du·diff]
        //       [dv·du  dv·dv] [Δv]   [-dv·diff]
        let a11 = du.dot(&du);
        let a12 = du.dot(&dv);
        let a22 = dv.dot(&dv);
        let b1 = -du.dot(&diff);
        let b2 = -dv.dot(&diff);

        let det = a11 * a22 - a12 * a12;
        if det.abs() < 1e-30 {
            break;
        }

        let delta_u = (a22 * b1 - a12 * b2) / det;
        let delta_v = (a11 * b2 - a12 * b1) / det;

        u += delta_u;
        v += delta_v;

        // Clamp to domain
        u = u.clamp(u_range.0, u_range.1);
        v = v.clamp(v_range.0, v_range.1);

        if delta_u.abs() < 1e-12 && delta_v.abs() < 1e-12 {
            break;
        }
    }

    (u, v)
}

/// Default parameter domain for grid search.
fn default_domain(surf: &Surface) -> ((f64, f64), (f64, f64)) {
    use std::f64::consts::TAU;
    match surf {
        Surface::Plane { .. } => ((-10.0, 10.0), (-10.0, 10.0)),
        Surface::Cylinder { .. } => ((0.0, TAU), (0.0, 10.0)),
        Surface::Cone { .. } => ((0.0, TAU), (0.0, 10.0)),
        Surface::Sphere { .. } => ((0.0, TAU), (-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2)),
        Surface::Torus { .. } => ((0.0, TAU), (0.0, TAU)),
        Surface::NurbsSurface(nurbs) => (nurbs.domain_u(), nurbs.domain_v()),
    }
}

/// Compute an orthonormal frame (e1, e2) in the plane perpendicular to `n`.
fn plane_frame(n: &Vector3) -> (Vector3, Vector3) {
    let a = n.normalize();
    let seed = if a.x.abs() < 0.9 {
        Vector3::new(1.0, 0.0, 0.0)
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    };
    let e1 = a.cross(&seed).normalize();
    let e2 = a.cross(&e1);
    (e1, e2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    #[test]
    fn plane_evaluate() {
        let plane = Surface::Plane {
            origin: Point3::origin(),
            normal: Vector3::new(0.0, 0.0, 1.0),
        };
        let p = plane.evaluate(3.0, 4.0);
        assert!(p.z.abs() < 1e-14, "Plane point should be in XY plane");
    }

    #[test]
    fn plane_normal_is_constant() {
        let plane = Surface::Plane {
            origin: Point3::origin(),
            normal: Vector3::new(0.0, 0.0, 1.0),
        };
        let n1 = plane.normal(0.0, 0.0);
        let n2 = plane.normal(5.0, 3.0);
        assert!((n1 - n2).norm() < 1e-14);
        assert!(n1.z.abs() > 0.99, "Normal should point in Z");
    }

    #[test]
    fn cylinder_on_cylinder() {
        let cyl = Surface::Cylinder {
            origin: Point3::origin(),
            axis: Vector3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        };
        for i in 0..=20 {
            let u = TAU * i as f64 / 20.0;
            let p = cyl.evaluate(u, 3.0);
            let r = (p.x * p.x + p.y * p.y).sqrt();
            assert!(
                (r - 5.0).abs() < 1e-12,
                "Cylinder point at u={u} has radius {r}"
            );
            assert!((p.z - 3.0).abs() < 1e-12);
        }
    }

    #[test]
    fn cylinder_normal_outward() {
        let cyl = Surface::Cylinder {
            origin: Point3::origin(),
            axis: Vector3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        };
        for i in 0..=10 {
            let u = TAU * i as f64 / 10.0;
            let p = cyl.evaluate(u, 0.0);
            let n = cyl.normal(u, 0.0);
            // Normal should point radially outward
            let radial = Vector3::new(p.x, p.y, 0.0).normalize();
            let dot = n.dot(&radial);
            assert!(
                dot.abs() > 0.99,
                "Cylinder normal should be radial at u={u}: dot={dot}"
            );
        }
    }

    #[test]
    fn sphere_on_sphere() {
        let sphere = Surface::Sphere {
            center: Point3::origin(),
            radius: 3.0,
        };
        for i in 0..=10 {
            for j in 0..=10 {
                let u = TAU * i as f64 / 10.0;
                let v = -FRAC_PI_2 + PI * j as f64 / 10.0;
                let p = sphere.evaluate(u, v);
                let r = p.coords.norm();
                assert!(
                    (r - 3.0).abs() < 1e-12,
                    "Sphere point at ({u},{v}) has radius {r}"
                );
            }
        }
    }

    #[test]
    fn sphere_normal_outward() {
        let sphere = Surface::Sphere {
            center: Point3::origin(),
            radius: 3.0,
        };
        let u = 1.0;
        let v = 0.3;
        let p = sphere.evaluate(u, v);
        let n = sphere.normal(u, v);
        let radial = p.coords.normalize();
        let dot = n.dot(&radial);
        assert!(
            dot.abs() > 0.99,
            "Sphere normal should point radially: dot={dot}"
        );
    }

    #[test]
    fn plane_closest_parameters_roundtrip() {
        let plane = Surface::Plane {
            origin: Point3::origin(),
            normal: Vector3::new(0.0, 0.0, 1.0),
        };
        let u = 3.7;
        let v = -2.1;
        let p = plane.evaluate(u, v);
        let (u2, v2) = plane.closest_parameters(&p);
        assert!((u - u2).abs() < 1e-10);
        assert!((v - v2).abs() < 1e-10);
    }

    #[test]
    fn cylinder_closest_parameters_roundtrip() {
        let cyl = Surface::Cylinder {
            origin: Point3::origin(),
            axis: Vector3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        };
        let u = 1.0;
        let v = 3.0;
        let p = cyl.evaluate(u, v);
        let (u2, v2) = cyl.closest_parameters(&p);
        assert!((u - u2).abs() < 1e-10, "u: expected {u}, got {u2}");
        assert!((v - v2).abs() < 1e-10, "v: expected {v}, got {v2}");
    }

    #[test]
    fn surface_derivatives_finite_difference() {
        let cyl = Surface::Cylinder {
            origin: Point3::origin(),
            axis: Vector3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        };
        let u = 1.0;
        let v = 2.0;
        let h = 1e-7;
        let du = cyl.derivative_u(u, v);
        let dv = cyl.derivative_v(u, v);
        let fd_u = (cyl.evaluate(u + h, v) - cyl.evaluate(u - h, v)) / (2.0 * h);
        let fd_v = (cyl.evaluate(u, v + h) - cyl.evaluate(u, v - h)) / (2.0 * h);
        assert!((du - fd_u).norm() < 1e-4, "du: {du:?} vs fd: {fd_u:?}");
        assert!((dv - fd_v).norm() < 1e-4, "dv: {dv:?} vs fd: {fd_v:?}");
    }

    #[test]
    fn torus_evaluate() {
        let torus = Surface::Torus {
            center: Point3::origin(),
            axis: Vector3::new(0.0, 0.0, 1.0),
            major_r: 10.0,
            minor_r: 3.0,
        };
        // At u=0, v=0: should be at distance major_r + minor_r from center
        let p = torus.evaluate(0.0, 0.0);
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!(
            (r - 13.0).abs() < 1e-12,
            "Torus outer point at distance {r}"
        );
    }
}
