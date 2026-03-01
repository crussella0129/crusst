//! Arena-based topology for the B-Rep kernel.
//!
//! Topology entities (Vertex, Edge, CoEdge, Wire, Face, Shell, Solid)
//! are stored in a central `TopoStore` and referenced via typed index handles.

pub mod store;
pub mod types;
pub mod validate;

pub use store::TopoStore;
pub use types::*;
pub use validate::{validate_solid, ValidationResult};
