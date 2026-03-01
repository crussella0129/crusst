//! Topology validation for B-Rep solids.
//!
//! Checks invariants that must hold for a valid manifold solid:
//! - Every edge has exactly 2 coedges (manifold condition)
//! - Every wire is closed (circular linked list)
//! - Euler formula: V - E + F = 2 (for genus-0 solids)

use super::store::TopoStore;
use super::types::*;
use std::collections::HashMap;

/// Result of topology validation.
#[derive(Debug)]
pub struct ValidationResult {
    pub valid: bool,
    pub errors: Vec<String>,
}

/// Validate that a solid satisfies B-Rep invariants.
pub fn validate_solid(store: &TopoStore, solid_id: SolidId) -> ValidationResult {
    let mut errors = Vec::new();

    let solid = store.solid(solid_id);
    let shell = store.shell(solid.outer_shell);

    // 1. Collect all coedges and check each edge has exactly 2
    let mut edge_coedge_count: HashMap<EdgeId, usize> = HashMap::new();
    for &face_id in &shell.faces {
        let face = store.face(face_id);
        for &coedge_id in &store.wire_coedges(face.outer_wire) {
            let coedge = store.coedge(coedge_id);
            *edge_coedge_count.entry(coedge.edge).or_insert(0) += 1;
        }
        for &inner_wire in &face.inner_wires {
            for &coedge_id in &store.wire_coedges(inner_wire) {
                let coedge = store.coedge(coedge_id);
                *edge_coedge_count.entry(coedge.edge).or_insert(0) += 1;
            }
        }
    }

    for (&edge_id, &count) in &edge_coedge_count {
        if count != 2 {
            errors.push(format!(
                "Edge {:?} has {count} coedges (expected 2 for manifold)",
                edge_id
            ));
        }
    }

    // 2. Check that paired coedges have opposite orientations
    let mut edge_orientations: HashMap<EdgeId, Vec<bool>> = HashMap::new();
    for &face_id in &shell.faces {
        let face = store.face(face_id);
        for &coedge_id in &store.wire_coedges(face.outer_wire) {
            let coedge = store.coedge(coedge_id);
            edge_orientations
                .entry(coedge.edge)
                .or_default()
                .push(coedge.forward);
        }
        for &inner_wire in &face.inner_wires {
            for &coedge_id in &store.wire_coedges(inner_wire) {
                let coedge = store.coedge(coedge_id);
                edge_orientations
                    .entry(coedge.edge)
                    .or_default()
                    .push(coedge.forward);
            }
        }
    }

    for (&edge_id, orients) in &edge_orientations {
        if orients.len() == 2 && orients[0] == orients[1] {
            errors.push(format!(
                "Edge {:?}: both coedges have same orientation (should be opposite for manifold)",
                edge_id
            ));
        }
    }

    // 3. Check all wires are closed
    for &face_id in &shell.faces {
        let face = store.face(face_id);
        if !is_wire_closed(store, face.outer_wire) {
            errors.push(format!(
                "Face {:?}: outer wire {:?} is not closed",
                face_id, face.outer_wire
            ));
        }
        for &inner_wire in &face.inner_wires {
            if !is_wire_closed(store, inner_wire) {
                errors.push(format!(
                    "Face {:?}: inner wire {:?} is not closed",
                    face_id, inner_wire
                ));
            }
        }
    }

    // 4. Euler formula: V - E + F = 2 (genus 0, no holes)
    let v = store.solid_vertices(solid_id).len() as i64;
    let e = store.solid_edges(solid_id).len() as i64;
    let f = store.solid_face_count(solid_id) as i64;
    let euler = v - e + f;

    if euler != 2 {
        errors.push(format!(
            "Euler formula V-E+F = {v}-{e}+{f} = {euler} (expected 2 for genus-0)"
        ));
    }

    ValidationResult {
        valid: errors.is_empty(),
        errors,
    }
}

/// Check that a wire forms a closed loop (the linked list returns to the start).
fn is_wire_closed(store: &TopoStore, wire_id: WireId) -> bool {
    let first = store.wire(wire_id).first_coedge;
    let mut current = store.coedge(first).next;
    let mut count = 1;
    let max_iter = 10000;

    while current != first {
        current = store.coedge(current).next;
        count += 1;
        if count > max_iter {
            return false; // Infinite loop — not closed
        }
    }
    true
}
