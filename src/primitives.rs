use nalgebra::Vector3;

/// Signed distance to a sphere.
/// Negative inside, zero on surface, positive outside.
pub fn sdf_sphere(point: Vector3<f64>, center: Vector3<f64>, radius: f64) -> f64 {
    (point - center).norm() - radius
}

/// Signed distance to an axis-aligned box.
/// `half_extents` is the half-size in each dimension.
pub fn sdf_box(point: Vector3<f64>, center: Vector3<f64>, half_extents: Vector3<f64>) -> f64 {
    let d = (point - center).abs() - half_extents;
    let outside = Vector3::new(d.x.max(0.0), d.y.max(0.0), d.z.max(0.0)).norm();
    let inside = d.x.max(d.y).max(d.z).min(0.0);
    outside + inside
}

/// Signed distance to a capped cylinder.
/// `base` is the center of the bottom cap, `axis` is the unit direction,
/// `radius` is the cross-section radius, `height` is the length along axis.
/// Exact signed distance to a capped cone (or truncated cone).
///
/// Connects two circles: one centered at `a` with radius `ra`, and another
/// at `b` with radius `rb`. For a sharp cone, set `rb = 0`. Handles the
/// base caps, tip, and conical wall correctly with no approximation.
///
/// Based on Inigo Quilez's exact capped cone SDF.
pub fn sdf_capped_cone(
    point: Vector3<f64>,
    a: Vector3<f64>,
    b: Vector3<f64>,
    ra: f64,
    rb: f64,
) -> f64 {
    let rba = rb - ra;
    let ba = b - a;
    let baba = ba.dot(&ba);
    let pa = point - a;
    let papa = pa.dot(&pa);
    let paba = pa.dot(&ba) / baba;
    let x = (papa - paba * paba * baba).max(0.0).sqrt();
    let cax = (x - if paba < 0.5 { ra } else { rb }).max(0.0);
    let cay = (paba - 0.5).abs() - 0.5;
    let k = rba * rba + baba;
    let f = ((rba * (x - ra) + paba * baba) / k).clamp(0.0, 1.0);
    let cbx = x - ra - f * rba;
    let cby = paba - f;
    let s = if cbx < 0.0 && cay < 0.0 { -1.0 } else { 1.0 };
    s * (cax * cax + cay * cay * baba)
        .min(cbx * cbx + cby * cby * baba)
        .sqrt()
}

pub fn sdf_cylinder(
    point: Vector3<f64>,
    base: Vector3<f64>,
    axis: Vector3<f64>,
    radius: f64,
    height: f64,
) -> f64 {
    let p = point - base;
    let along = p.dot(&axis);
    let perp = (p - axis * along).norm();
    let d_radial = perp - radius;
    let d_axial = (along - height / 2.0).abs() - height / 2.0;
    if d_radial > 0.0 && d_axial > 0.0 {
        (d_radial * d_radial + d_axial * d_axial).sqrt()
    } else {
        d_radial.max(d_axial)
    }
}
