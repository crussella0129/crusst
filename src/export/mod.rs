//! Export writers for B-Rep shapes.
//!
//! Provides writers for common CAD and mesh interchange formats:
//! - **STL** — Binary triangle mesh (tessellation-based)
//! - **OBJ** — Wavefront text format with normals (tessellation-based)
//! - **STEP AP203** — Exact B-Rep geometry (no tessellation needed)
//! - **3MF** — XML + triangle mesh with metadata
//! - **IGES** — Legacy CAD interchange format

pub mod stl;
pub mod obj;
pub mod step;
pub mod threemf;

pub use stl::write_stl;
pub use obj::write_obj;
pub use step::write_step;
pub use threemf::write_3mf;
