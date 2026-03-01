//! Shape operations: fillet, chamfer, shell, sweep, loft.
//!
//! These operations modify B-Rep topology by adding, removing, or reshaping
//! faces and edges. They are the bread and butter of parametric CAD modeling.

pub mod fillet;
pub mod chamfer;
pub mod shell;
pub mod sweep;
pub mod loft;

pub use fillet::fillet_edges;
pub use chamfer::chamfer_edges;
pub use shell::shell_solid;
pub use sweep::sweep_profile;
pub use loft::loft_profiles;
