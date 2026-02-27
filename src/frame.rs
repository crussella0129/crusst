use nalgebra::{Matrix3, Vector3};

/// A local coordinate frame at a point on a path.
/// The tangent (T), normal (N), and binormal (B) vectors.
pub struct FrenetFrame {
    pub tangent: Vector3<f64>,
    pub normal: Vector3<f64>,
    pub binormal: Vector3<f64>,
}

impl FrenetFrame {
    /// Compute a robust frame for a tangent direction.
    /// Uses an arbitrary perpendicular for the normal if curvature is zero.
    pub fn from_tangent(tangent: Vector3<f64>) -> Self {
        let t = tangent.normalize();

        // Find a vector not parallel to tangent
        let arbitrary = if t.x.abs() < 0.9 {
            Vector3::new(1.0, 0.0, 0.0)
        } else {
            Vector3::new(0.0, 1.0, 0.0)
        };

        let normal = t.cross(&arbitrary).normalize();
        let binormal = t.cross(&normal).normalize();

        Self {
            tangent: t,
            normal,
            binormal,
        }
    }

    /// Build a rotation matrix that maps local frame to world frame.
    /// Local X -> normal, local Y -> binormal, local Z -> tangent.
    pub fn to_matrix(&self) -> Matrix3<f64> {
        Matrix3::from_columns(&[self.normal, self.binormal, self.tangent])
    }
}
