//! Boolean operations on B-Rep solids.
//!
//! Boolean operations (union, intersect, subtract) are the most complex
//! algorithms in a B-Rep kernel. They require:
//! - Surface-surface intersection curve computation
//! - Face classification (inside/outside/on)
//! - Topology surgery (splitting faces at intersection curves)
//! - Result assembly from classified face fragments
//!
//! This module defines the `BooleanEngine` trait and provides a stub
//! implementation. Full implementation is deferred as it requires
//! robust surface intersection — the hardest problem in computational geometry.

use crate::topo::*;

/// Trait for boolean operations on B-Rep solids.
///
/// Implementations must handle:
/// - Surface-surface intersection computation
/// - Face splitting at intersection curves
/// - Inside/outside classification
/// - Result topology assembly
pub trait BooleanEngine {
    /// Compute the union of two solids (A ∪ B).
    fn union(
        &self,
        store: &mut TopoStore,
        a: SolidId,
        b: SolidId,
    ) -> Result<SolidId, BooleanError>;

    /// Compute the intersection of two solids (A ∩ B).
    fn intersect(
        &self,
        store: &mut TopoStore,
        a: SolidId,
        b: SolidId,
    ) -> Result<SolidId, BooleanError>;

    /// Subtract solid B from solid A (A \ B).
    fn subtract(
        &self,
        store: &mut TopoStore,
        a: SolidId,
        b: SolidId,
    ) -> Result<SolidId, BooleanError>;
}

/// Errors that can occur during boolean operations.
#[derive(Debug)]
pub enum BooleanError {
    /// Boolean operations are not yet implemented.
    NotImplemented,
    /// The solids do not overlap.
    NoOverlap,
    /// Surface-surface intersection failed.
    IntersectionFailed,
    /// The result would not be a valid manifold solid.
    NonManifoldResult,
}

impl std::fmt::Display for BooleanError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BooleanError::NotImplemented => write!(f, "Boolean operations not yet implemented"),
            BooleanError::NoOverlap => write!(f, "Solids do not overlap"),
            BooleanError::IntersectionFailed => write!(f, "Surface intersection computation failed"),
            BooleanError::NonManifoldResult => write!(f, "Result would not be manifold"),
        }
    }
}

impl std::error::Error for BooleanError {}

/// Stub boolean engine that returns `NotImplemented` for all operations.
///
/// This serves as a placeholder until surface-surface intersection
/// algorithms are implemented.
pub struct StubBooleanEngine;

impl BooleanEngine for StubBooleanEngine {
    fn union(
        &self,
        _store: &mut TopoStore,
        _a: SolidId,
        _b: SolidId,
    ) -> Result<SolidId, BooleanError> {
        Err(BooleanError::NotImplemented)
    }

    fn intersect(
        &self,
        _store: &mut TopoStore,
        _a: SolidId,
        _b: SolidId,
    ) -> Result<SolidId, BooleanError> {
        Err(BooleanError::NotImplemented)
    }

    fn subtract(
        &self,
        _store: &mut TopoStore,
        _a: SolidId,
        _b: SolidId,
    ) -> Result<SolidId, BooleanError> {
        Err(BooleanError::NotImplemented)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitive;

    #[test]
    fn stub_union_returns_not_implemented() {
        let engine = StubBooleanEngine;
        let mut store = TopoStore::new();
        let a = primitive::make_box(&mut store, 5.0, 3.0, 8.0);
        let b = primitive::make_sphere(&mut store, 4.0);
        let result = engine.union(&mut store, a, b);
        assert!(matches!(result, Err(BooleanError::NotImplemented)));
    }

    #[test]
    fn stub_intersect_returns_not_implemented() {
        let engine = StubBooleanEngine;
        let mut store = TopoStore::new();
        let a = primitive::make_box(&mut store, 5.0, 5.0, 5.0);
        let b = primitive::make_cylinder(&mut store, 3.0, 10.0);
        let result = engine.intersect(&mut store, a, b);
        assert!(matches!(result, Err(BooleanError::NotImplemented)));
    }

    #[test]
    fn stub_subtract_returns_not_implemented() {
        let engine = StubBooleanEngine;
        let mut store = TopoStore::new();
        let a = primitive::make_box(&mut store, 10.0, 10.0, 10.0);
        let b = primitive::make_sphere(&mut store, 5.0);
        let result = engine.subtract(&mut store, a, b);
        assert!(matches!(result, Err(BooleanError::NotImplemented)));
    }
}
