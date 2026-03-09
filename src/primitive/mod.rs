//! Primitive solid constructors.
//!
//! Each function builds a complete B-Rep solid (topology + geometry) in a `TopoStore`.

mod box3;
mod capsule;
mod cone;
mod cylinder;
mod sphere;
mod torus;
mod wedge;

pub use box3::make_box;
pub use capsule::make_capsule;
pub use cone::make_cone;
pub use cylinder::make_cylinder;
pub use sphere::make_sphere;
pub use torus::make_torus;
pub use wedge::make_wedge;
