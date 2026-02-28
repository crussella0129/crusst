//! Adaptive dual contouring mesher.
//!
//! Ties together the octree (adaptive subdivision with interval pruning),
//! QEF solver (optimal vertex placement), and edge-crossing detection to
//! extract a triangle mesh from an SDF DAG node.
//!
//! ## Algorithm overview
//!
//! 1. **Build octree** — adaptive subdivision with interval arithmetic pruning.
//! 2. **Place vertices** — for each surface leaf cell, find edge crossings via
//!    bisection, compute gradients at crossings, solve QEF for the optimal
//!    vertex position.
//! 3. **Generate faces** — for each sign-changing edge shared by multiple cells,
//!    connect their QEF vertices into triangles.
//! 4. **Compute normals** — evaluate the SDF gradient at each vertex position.

use crate::dag::SdfNode;
use crate::octree::{Octree, OctreeCell};
use crate::qef::solve_qef;
use crate::shape::Sdf;
use crate::types::{BBox3, MeshSettings, TriangleMesh};
use nalgebra::Vector3;
use std::collections::HashMap;

/// Quantization scale for discretizing floating-point coordinates to integer
/// grid keys. 2^30 gives ~9.3e-10 resolution per unit — sufficient for
/// depths up to 30 without collisions.
const QUANT_SCALE: f64 = (1u64 << 30) as f64;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Extract a triangle mesh from an SDF DAG node using adaptive dual contouring.
///
/// This is the primary entry point for the meshing pipeline.
pub fn extract_mesh_adaptive(
    node: &SdfNode,
    bbox: &BBox3,
    settings: &MeshSettings,
) -> TriangleMesh {
    // Pad the bounding box to guarantee the surface is fully contained.
    // Use the larger of 5% of bbox size or 2 cells at the finest resolution.
    // This prevents boundary edges (edges with < 3 cells) that create holes.
    let bbox_size = bbox.size();
    let cell_size = bbox_size.x.max(bbox_size.y).max(bbox_size.z)
        / (1u64 << settings.max_depth) as f64;
    let min_pad = cell_size * 2.0;
    let pct_pad = bbox_size * 0.05;
    let pad = Vector3::new(
        pct_pad.x.max(min_pad),
        pct_pad.y.max(min_pad),
        pct_pad.z.max(min_pad),
    );
    let padded = BBox3::new(bbox.min - pad, bbox.max + pad);

    // Phase 1: Build the adaptive octree and enforce 2:1 balance.
    let mut tree = Octree::build(node, &padded, settings);
    tree.balance(node);

    // Phase 2: Collect surface leaf cells and place QEF vertices.
    let surface_cells = tree.surface_cells();
    if surface_cells.is_empty() {
        return TriangleMesh {
            vertices: Vec::new(),
            normals: Vec::new(),
            indices: Vec::new(),
        };
    }

    // Assign each surface cell a vertex index and compute its QEF vertex.
    let mut vertices: Vec<Vector3<f64>> = Vec::with_capacity(surface_cells.len());
    // Map from cell key (discretized min corner + depth) to vertex index.
    let mut cell_vertex_map: HashMap<CellKey, u32> = HashMap::new();

    for cell in &surface_cells {
        let key = cell_key(cell);
        let vertex = compute_cell_vertex(node, cell, settings.edge_tolerance);
        let idx = vertices.len() as u32;
        vertices.push(vertex);
        cell_vertex_map.insert(key, idx);
    }

    // Phase 3: Generate faces from shared sign-changing edges.
    //
    // For each surface leaf cell, enumerate its 12 edges. For each edge with
    // a sign change, register the cell in an edge -> cells map. Then for each
    // edge shared by 3+ cells (or exactly 2 in boundary cases), emit triangles.
    //
    // T-junction resolution: In an adaptive octree with 2:1 balance, a large
    // cell's edge may span what is two half-edges in smaller neighbor cells.
    // If we register the full edge, it won't match the half-edges, creating
    // boundary (non-manifold) edges. Fix: collect all leaf corner positions
    // and split any cell edge whose midpoint coincides with a finer-level
    // corner.
    let all_leaves = tree.collect_all_leaves();
    let mut corner_set: std::collections::HashSet<(i64, i64, i64)> =
        std::collections::HashSet::new();
    for leaf in &all_leaves {
        let cs = leaf.bbox.corners();
        for c in &cs {
            let qx = (c.x * QUANT_SCALE).round() as i64;
            let qy = (c.y * QUANT_SCALE).round() as i64;
            let qz = (c.z * QUANT_SCALE).round() as i64;
            corner_set.insert((qx, qy, qz));
        }
    }

    let mut edge_cells: HashMap<EdgeKey, Vec<(CellKey, &OctreeCell)>> = HashMap::new();

    for cell in &surface_cells {
        let key = cell_key(cell);
        let corners = cell.bbox.corners();

        for &(ci0, ci1) in &CELL_EDGES {
            let v0 = cell.corners[ci0];
            let v1 = cell.corners[ci1];
            // Only process edges with a sign change.
            if (v0 > 0.0) == (v1 > 0.0) {
                continue;
            }
            // Check if this edge's midpoint is a corner of a finer neighbor.
            let mid = (corners[ci0] + corners[ci1]) * 0.5;
            let qmx = (mid.x * QUANT_SCALE).round() as i64;
            let qmy = (mid.y * QUANT_SCALE).round() as i64;
            let qmz = (mid.z * QUANT_SCALE).round() as i64;

            if corner_set.contains(&(qmx, qmy, qmz)) {
                // Edge midpoint is a finer-level corner — split into two
                // half-edges and register each sub-edge that has a sign change.
                let vm = node.evaluate(mid);
                // Sub-edge A -> M
                if (v0 > 0.0) != (vm > 0.0) {
                    let ek = edge_key(&corners[ci0], &mid);
                    edge_cells.entry(ek).or_default().push((key, *cell));
                }
                // Sub-edge M -> B
                if (vm > 0.0) != (v1 > 0.0) {
                    let ek = edge_key(&mid, &corners[ci1]);
                    edge_cells.entry(ek).or_default().push((key, *cell));
                }
            } else {
                let ek = edge_key(&corners[ci0], &corners[ci1]);
                edge_cells.entry(ek).or_default().push((key, *cell));
            }
        }
    }

    let mut indices: Vec<u32> = Vec::new();

    for (ek, cells) in &edge_cells {
        // Deduplicate cells by their key (multiple edges of the same cell can
        // map to the same edge key in adaptive grids).
        let mut unique: Vec<(CellKey, &OctreeCell)> = Vec::new();
        let mut seen = std::collections::HashSet::new();
        for &(ref k, c) in cells {
            if seen.insert(*k) {
                unique.push((*k, c));
            }
        }

        if unique.len() < 2 {
            continue;
        }

        // Look up vertex indices for each cell.
        let vert_indices: Vec<u32> = unique
            .iter()
            .filter_map(|(k, _)| cell_vertex_map.get(k).copied())
            .collect();

        if vert_indices.len() < 2 {
            continue;
        }

        // Determine winding order: the triangle normal should point from
        // inside (negative SDF) to outside (positive SDF).
        // Use the edge direction and sign to figure out consistent winding.
        let sign_positive_first = {
            // Decode edge key to get approximate edge direction
            // The first endpoint of the edge key is the one with smaller coords.
            // We use the SDF sign at the first corner of the first cell's edge.
            // A simpler approach: evaluate SDF at the edge midpoint or use
            // the corner signs.
            let first_cell = unique[0].1;
            let fc = first_cell.bbox.corners();
            // Find which corner pair corresponds to this edge
            let mut sign = false;
            for &(ci0, ci1) in &CELL_EDGES {
                let test_ek = edge_key(&fc[ci0], &fc[ci1]);
                if test_ek == *ek {
                    // The sign of corner ci0 tells us orientation
                    sign = first_cell.corners[ci0] > 0.0;
                    break;
                }
            }
            sign
        };

        // Sort vertices by angle around the shared edge to get consistent
        // winding for quads.
        let sorted_indices = sort_vertices_around_edge(ek, &vert_indices, &vertices);

        // Emit triangles as a fan from the first vertex.
        emit_fan(&sorted_indices, sign_positive_first, &mut indices);
    }

    // Phase 3.5: Fix winding order using SDF gradient.
    fix_winding(&vertices, &mut indices, node);

    // Phase 3.6: Fix non-manifold edges from overlapping fans at CSG seams.
    fix_non_manifold(&vertices, &mut indices, node);

    // Phase 4: Compute per-vertex normals.
    let normals: Vec<Vector3<f64>> = vertices.iter().map(|v| node.gradient(*v)).collect();

    let mesh = TriangleMesh {
        vertices,
        normals,
        indices,
    };

    // Post-process: split vertices at sharp edges so each side of a crease
    // gets its own normal. This fixes "chewy" visual artifacts on shapes
    // with hard feature edges (e.g. capped cones, boxes).
    mesh.split_sharp_edges(35.0)
}

/// Extract a mesh from a `&dyn Sdf` trait object (compatibility wrapper).
///
/// This builds an octree by evaluating the SDF at corners directly (no interval
/// pruning or analytical gradients), then runs the same DC algorithm with
/// central-difference gradients.
pub fn extract_mesh_from_sdf(sdf: &dyn Sdf, bbox: &BBox3, settings: &MeshSettings) -> TriangleMesh {
    // Build the octree manually using the trait object.
    let root = build_cell_from_sdf(sdf, bbox, 0, settings);
    let tree_wrapper = OctreeWrapper { root };

    let surface_cells = tree_wrapper.surface_cells();
    if surface_cells.is_empty() {
        return TriangleMesh {
            vertices: Vec::new(),
            normals: Vec::new(),
            indices: Vec::new(),
        };
    }

    let mut vertices: Vec<Vector3<f64>> = Vec::with_capacity(surface_cells.len());
    let mut cell_vertex_map: HashMap<CellKey, u32> = HashMap::new();

    for cell in &surface_cells {
        let key = cell_key(cell);
        let vertex = compute_cell_vertex_sdf(sdf, cell, settings.edge_tolerance);
        let idx = vertices.len() as u32;
        vertices.push(vertex);
        cell_vertex_map.insert(key, idx);
    }

    // Phase 3: Generate faces (same T-junction-aware logic as above).
    let all_leaves = tree_wrapper.collect_all_leaves();
    let mut corner_set: std::collections::HashSet<(i64, i64, i64)> =
        std::collections::HashSet::new();
    for leaf in &all_leaves {
        let cs = leaf.bbox.corners();
        for c in &cs {
            let qx = (c.x * QUANT_SCALE).round() as i64;
            let qy = (c.y * QUANT_SCALE).round() as i64;
            let qz = (c.z * QUANT_SCALE).round() as i64;
            corner_set.insert((qx, qy, qz));
        }
    }

    let mut edge_cells: HashMap<EdgeKey, Vec<(CellKey, &OctreeCell)>> = HashMap::new();

    for cell in &surface_cells {
        let key = cell_key(cell);
        let corners = cell.bbox.corners();

        for &(ci0, ci1) in &CELL_EDGES {
            let v0 = cell.corners[ci0];
            let v1 = cell.corners[ci1];
            if (v0 > 0.0) == (v1 > 0.0) {
                continue;
            }
            // T-junction resolution: split edge if midpoint is a finer corner.
            let mid = (corners[ci0] + corners[ci1]) * 0.5;
            let qmx = (mid.x * QUANT_SCALE).round() as i64;
            let qmy = (mid.y * QUANT_SCALE).round() as i64;
            let qmz = (mid.z * QUANT_SCALE).round() as i64;

            if corner_set.contains(&(qmx, qmy, qmz)) {
                let vm = sdf.evaluate(mid);
                if (v0 > 0.0) != (vm > 0.0) {
                    let ek = edge_key(&corners[ci0], &mid);
                    edge_cells.entry(ek).or_default().push((key, *cell));
                }
                if (vm > 0.0) != (v1 > 0.0) {
                    let ek = edge_key(&mid, &corners[ci1]);
                    edge_cells.entry(ek).or_default().push((key, *cell));
                }
            } else {
                let ek = edge_key(&corners[ci0], &corners[ci1]);
                edge_cells.entry(ek).or_default().push((key, *cell));
            }
        }
    }

    let mut indices: Vec<u32> = Vec::new();

    for (ek, cells) in &edge_cells {
        let mut unique: Vec<(CellKey, &OctreeCell)> = Vec::new();
        let mut seen = std::collections::HashSet::new();
        for &(ref k, c) in cells {
            if seen.insert(*k) {
                unique.push((*k, c));
            }
        }

        if unique.len() < 2 {
            continue;
        }

        let vert_indices: Vec<u32> = unique
            .iter()
            .filter_map(|(k, _)| cell_vertex_map.get(k).copied())
            .collect();

        if vert_indices.len() < 2 {
            continue;
        }

        let sign_positive_first = {
            let first_cell = unique[0].1;
            let fc = first_cell.bbox.corners();
            let mut sign = false;
            for &(ci0, ci1) in &CELL_EDGES {
                let test_ek = edge_key(&fc[ci0], &fc[ci1]);
                if test_ek == *ek {
                    sign = first_cell.corners[ci0] > 0.0;
                    break;
                }
            }
            sign
        };

        let sorted_indices = sort_vertices_around_edge(ek, &vert_indices, &vertices);
        emit_fan(&sorted_indices, sign_positive_first, &mut indices);
    }

    // Phase 3.5: Fix winding order using SDF gradient (central differences).
    fix_winding_sdf(&vertices, &mut indices, sdf);

    // Phase 3.6: Fix non-manifold edges from overlapping fans at CSG seams.
    fix_non_manifold_sdf(&vertices, &mut indices, sdf);

    // Phase 4: Normals — prefer analytical gradients, fall back to central differences.
    let normals: Vec<Vector3<f64>> = vertices
        .iter()
        .map(|v| {
            sdf.gradient(*v)
                .unwrap_or_else(|| central_diff_gradient_sdf(sdf, *v))
        })
        .collect();

    let mesh = TriangleMesh {
        vertices,
        normals,
        indices,
    };

    // Post-process: split vertices at sharp edges so each side of a crease
    // gets its own normal. This fixes "chewy" visual artifacts on shapes
    // with hard feature edges (e.g. capped cones, boxes).
    mesh.split_sharp_edges(35.0)
}

// ---------------------------------------------------------------------------
// Cell identification
// ---------------------------------------------------------------------------

/// Unique key for an octree leaf cell, based on its discretized min corner
/// and depth level. Two cells at the same position and depth are the same cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct CellKey {
    // Discretized coordinates (scaled to integer grid at this depth level)
    ix: i64,
    iy: i64,
    iz: i64,
    depth: u8,
}

fn cell_key(cell: &OctreeCell) -> CellKey {
    // Use a large multiplier to discretize floating-point coordinates to
    // integer grid coordinates. At depth d, cells have size S/2^d where S
    // is the root bbox size. QUANT_SCALE (2^30) gives integer precision
    // that works for depths up to 30.
    let scale = QUANT_SCALE;
    CellKey {
        ix: (cell.bbox.min.x * scale).round() as i64,
        iy: (cell.bbox.min.y * scale).round() as i64,
        iz: (cell.bbox.min.z * scale).round() as i64,
        depth: cell.depth,
    }
}

// ---------------------------------------------------------------------------
// Edge identification
// ---------------------------------------------------------------------------

/// Unique key for a grid edge, identified by the discretized coordinates of
/// its two endpoints (ordered so that the smaller one comes first).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct EdgeKey {
    // First endpoint (smaller)
    ax: i64,
    ay: i64,
    az: i64,
    // Second endpoint (larger)
    bx: i64,
    by: i64,
    bz: i64,
}

fn edge_key(p0: &Vector3<f64>, p1: &Vector3<f64>) -> EdgeKey {
    let scale = QUANT_SCALE;
    let ax = (p0.x * scale).round() as i64;
    let ay = (p0.y * scale).round() as i64;
    let az = (p0.z * scale).round() as i64;
    let bx = (p1.x * scale).round() as i64;
    let by = (p1.y * scale).round() as i64;
    let bz = (p1.z * scale).round() as i64;

    // Order endpoints lexicographically so that the same edge from
    // different cells produces the same key.
    if (ax, ay, az) <= (bx, by, bz) {
        EdgeKey {
            ax,
            ay,
            az,
            bx,
            by,
            bz,
        }
    } else {
        EdgeKey {
            ax: bx,
            ay: by,
            az: bz,
            bx: ax,
            by: ay,
            bz: az,
        }
    }
}

// ---------------------------------------------------------------------------
// The 12 edges of a cell (as pairs of corner indices)
// ---------------------------------------------------------------------------

/// The 12 edges of an axis-aligned box, indexed by corner indices (0-7).
///
/// Corner ordering follows BBox3::corners():
/// ```text
/// 0: (min.x, min.y, min.z)    4: (min.x, min.y, max.z)
/// 1: (max.x, min.y, min.z)    5: (max.x, min.y, max.z)
/// 2: (min.x, max.y, min.z)    6: (min.x, max.y, max.z)
/// 3: (max.x, max.y, min.z)    7: (max.x, max.y, max.z)
/// ```
///
/// Edges along X (4): 0-1, 2-3, 4-5, 6-7
/// Edges along Y (4): 0-2, 1-3, 4-6, 5-7
/// Edges along Z (4): 0-4, 1-5, 2-6, 3-7
const CELL_EDGES: [(usize, usize); 12] = [
    // X edges
    (0, 1),
    (2, 3),
    (4, 5),
    (6, 7),
    // Y edges
    (0, 2),
    (1, 3),
    (4, 6),
    (5, 7),
    // Z edges
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7),
];

// ---------------------------------------------------------------------------
// Bisection for finding zero-crossings on edges
// ---------------------------------------------------------------------------

fn find_crossing_node(
    node: &SdfNode,
    p0: Vector3<f64>,
    p1: Vector3<f64>,
    v0: f64,
    _v1: f64,
    tolerance: f64,
) -> Vector3<f64> {
    let mut a = p0;
    let mut b = p1;
    let mut va = v0;
    for _ in 0..32 {
        let mid = (a + b) * 0.5;
        let vm = node.evaluate(mid);
        if vm.abs() < tolerance {
            return mid;
        }
        if (va > 0.0) == (vm > 0.0) {
            a = mid;
            va = vm;
        } else {
            b = mid;
        }
    }
    (a + b) * 0.5
}

fn find_crossing_sdf(
    sdf: &dyn Sdf,
    p0: Vector3<f64>,
    p1: Vector3<f64>,
    v0: f64,
    _v1: f64,
    tolerance: f64,
) -> Vector3<f64> {
    let mut a = p0;
    let mut b = p1;
    let mut va = v0;
    for _ in 0..32 {
        let mid = (a + b) * 0.5;
        let vm = sdf.evaluate(mid);
        if vm.abs() < tolerance {
            return mid;
        }
        if (va > 0.0) == (vm > 0.0) {
            a = mid;
            va = vm;
        } else {
            b = mid;
        }
    }
    (a + b) * 0.5
}

// ---------------------------------------------------------------------------
// QEF vertex computation for a cell
// ---------------------------------------------------------------------------

/// Compute the QEF vertex for a surface leaf cell using the SdfNode DAG.
fn compute_cell_vertex(node: &SdfNode, cell: &OctreeCell, tolerance: f64) -> Vector3<f64> {
    let corners = cell.bbox.corners();
    let mut positions: Vec<Vector3<f64>> = Vec::new();
    let mut normals: Vec<Vector3<f64>> = Vec::new();

    for &(ci0, ci1) in &CELL_EDGES {
        let v0 = cell.corners[ci0];
        let v1 = cell.corners[ci1];
        if (v0 > 0.0) == (v1 > 0.0) {
            continue;
        }
        let crossing = find_crossing_node(node, corners[ci0], corners[ci1], v0, v1, tolerance);
        let normal = node.gradient(crossing);
        positions.push(crossing);
        normals.push(normal);
    }

    if positions.is_empty() {
        return cell.bbox.center();
    }

    solve_qef(&positions, &normals, &cell.bbox)
}

/// Compute the QEF vertex for a surface leaf cell using a &dyn Sdf.
fn compute_cell_vertex_sdf(sdf: &dyn Sdf, cell: &OctreeCell, tolerance: f64) -> Vector3<f64> {
    let corners = cell.bbox.corners();
    let mut positions: Vec<Vector3<f64>> = Vec::new();
    let mut normals: Vec<Vector3<f64>> = Vec::new();

    for &(ci0, ci1) in &CELL_EDGES {
        let v0 = cell.corners[ci0];
        let v1 = cell.corners[ci1];
        if (v0 > 0.0) == (v1 > 0.0) {
            continue;
        }
        let crossing = find_crossing_sdf(sdf, corners[ci0], corners[ci1], v0, v1, tolerance);
        let normal = sdf.gradient(crossing)
            .unwrap_or_else(|| central_diff_gradient_sdf(sdf, crossing));
        positions.push(crossing);
        normals.push(normal);
    }

    if positions.is_empty() {
        return cell.bbox.center();
    }

    solve_qef(&positions, &normals, &cell.bbox)
}

// ---------------------------------------------------------------------------
// Central differences gradient for &dyn Sdf
// ---------------------------------------------------------------------------

fn central_diff_gradient_sdf(sdf: &dyn Sdf, point: Vector3<f64>) -> Vector3<f64> {
    let eps = 1e-6;
    let dx = sdf.evaluate(point + Vector3::new(eps, 0.0, 0.0))
        - sdf.evaluate(point - Vector3::new(eps, 0.0, 0.0));
    let dy = sdf.evaluate(point + Vector3::new(0.0, eps, 0.0))
        - sdf.evaluate(point - Vector3::new(0.0, eps, 0.0));
    let dz = sdf.evaluate(point + Vector3::new(0.0, 0.0, eps))
        - sdf.evaluate(point - Vector3::new(0.0, 0.0, eps));
    let g = Vector3::new(dx, dy, dz);
    let len = g.norm();
    if len > 1e-10 {
        g / len
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    }
}

// ---------------------------------------------------------------------------
// Sorting vertices around an edge for consistent fan triangulation
// ---------------------------------------------------------------------------

/// Sort vertex indices by angle around the shared edge axis so that the
/// resulting fan produces non-self-intersecting triangles.
fn sort_vertices_around_edge(
    ek: &EdgeKey,
    vert_indices: &[u32],
    vertices: &[Vector3<f64>],
) -> Vec<u32> {
    if vert_indices.len() <= 2 {
        return vert_indices.to_vec();
    }

    // Reconstruct approximate edge direction from the edge key.
    let scale = QUANT_SCALE;
    let ea = Vector3::new(
        ek.ax as f64 / scale,
        ek.ay as f64 / scale,
        ek.az as f64 / scale,
    );
    let eb = Vector3::new(
        ek.bx as f64 / scale,
        ek.by as f64 / scale,
        ek.bz as f64 / scale,
    );
    let edge_dir = (eb - ea).normalize();
    let edge_mid = (ea + eb) * 0.5;

    // Project each vertex onto the plane perpendicular to the edge at the edge midpoint.
    // Compute angles for sorting.
    let mut angles: Vec<(f64, u32)> = Vec::with_capacity(vert_indices.len());

    // Choose a reference direction in the perpendicular plane.
    let ref_dir = {
        let v0 = vertices[vert_indices[0] as usize];
        let d = v0 - edge_mid;
        let projected = d - edge_dir * d.dot(&edge_dir);
        let len = projected.norm();
        if len > 1e-12 {
            projected / len
        } else {
            // Fallback: pick any direction perpendicular to edge_dir
            let perp = if edge_dir.x.abs() < 0.9 {
                Vector3::new(1.0, 0.0, 0.0)
            } else {
                Vector3::new(0.0, 1.0, 0.0)
            };
            let cross = edge_dir.cross(&perp);
            cross.normalize()
        }
    };
    let binormal = edge_dir.cross(&ref_dir);

    for &vi in vert_indices {
        let v = vertices[vi as usize];
        let d = v - edge_mid;
        let projected = d - edge_dir * d.dot(&edge_dir);
        let x = projected.dot(&ref_dir);
        let y = projected.dot(&binormal);
        let angle = y.atan2(x);
        angles.push((angle, vi));
    }

    angles.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    angles.iter().map(|&(_, vi)| vi).collect()
}

// ---------------------------------------------------------------------------
// Triangle fan emission
// ---------------------------------------------------------------------------

/// Emit a triangle fan from a sorted ring of vertex indices.
///
/// `sign_positive_first` controls winding: if the first corner of the edge
/// (in the edge key's ordering) is positive (outside), we emit CCW triangles
/// so that the face normal points outward.
fn emit_fan(sorted: &[u32], sign_positive_first: bool, indices: &mut Vec<u32>) {
    let n = sorted.len();
    if n < 3 {
        // With only 2 vertices we cannot form a triangle. This is a boundary edge.
        // Some adaptive DC implementations skip these.
        if n == 2 {
            // Degenerate: skip (boundary of the octree domain).
            return;
        }
        return;
    }

    for i in 1..(n - 1) {
        if sign_positive_first {
            // CCW winding (normal points from negative to positive)
            indices.push(sorted[0]);
            indices.push(sorted[i]);
            indices.push(sorted[i + 1]);
        } else {
            // Reverse winding
            indices.push(sorted[0]);
            indices.push(sorted[i + 1]);
            indices.push(sorted[i]);
        }
    }
}

// ---------------------------------------------------------------------------
// Winding fix — ensure all face normals agree with SDF gradient
// ---------------------------------------------------------------------------

/// Fix triangle winding so all face normals point outward (agree with SDF gradient).
///
/// For each triangle, computes the face normal via cross product and compares
/// it against the SDF gradient at the centroid. If they disagree, the triangle
/// vertices are swapped to flip the winding.
fn fix_winding(vertices: &[Vector3<f64>], indices: &mut [u32], node: &SdfNode) {
    for tri in indices.chunks_exact_mut(3) {
        let v0 = vertices[tri[0] as usize];
        let v1 = vertices[tri[1] as usize];
        let v2 = vertices[tri[2] as usize];
        let face_normal = (v1 - v0).cross(&(v2 - v0));
        // Skip degenerate triangles
        if face_normal.norm_squared() < 1e-30 {
            continue;
        }
        let centroid = (v0 + v1 + v2) / 3.0;
        let sdf_gradient = node.gradient(centroid);
        if face_normal.dot(&sdf_gradient) < 0.0 {
            tri.swap(1, 2);
        }
    }
}

/// Fix winding for meshes built from &dyn Sdf (uses gradient or central differences).
fn fix_winding_sdf(vertices: &[Vector3<f64>], indices: &mut [u32], sdf: &dyn Sdf) {
    for tri in indices.chunks_exact_mut(3) {
        let v0 = vertices[tri[0] as usize];
        let v1 = vertices[tri[1] as usize];
        let v2 = vertices[tri[2] as usize];
        let face_normal = (v1 - v0).cross(&(v2 - v0));
        if face_normal.norm_squared() < 1e-30 {
            continue;
        }
        let centroid = (v0 + v1 + v2) / 3.0;
        let sdf_gradient = sdf
            .gradient(centroid)
            .unwrap_or_else(|| central_diff_gradient_sdf(sdf, centroid));
        if face_normal.dot(&sdf_gradient) < 0.0 {
            tri.swap(1, 2);
        }
    }
}

// ---------------------------------------------------------------------------
// Manifold repair — remove excess triangles at non-manifold edges
// ---------------------------------------------------------------------------

/// Quantize a vertex position to an integer key for position-based edge matching.
fn quant_vertex(v: &Vector3<f64>) -> (i64, i64, i64) {
    // QUANT_SCALE (2^30) gives ~1e-9 resolution, matching the edge_key precision.
    (
        (v.x * QUANT_SCALE).round() as i64,
        (v.y * QUANT_SCALE).round() as i64,
        (v.z * QUANT_SCALE).round() as i64,
    )
}

type PosEdgeKey = ((i64, i64, i64), (i64, i64, i64));

fn pos_edge_key(a: &Vector3<f64>, b: &Vector3<f64>) -> PosEdgeKey {
    let qa = quant_vertex(a);
    let qb = quant_vertex(b);
    if qa <= qb { (qa, qb) } else { (qb, qa) }
}

/// Remove triangles that create non-manifold edges (shared by >2 triangles).
///
/// At hard CSG intersection seams, the DC algorithm can produce overlapping
/// triangle fans that share geometric edges. This pass:
/// 1. Removes degenerate triangles (coincident vertices).
/// 2. Iteratively removes excess triangles at over-shared edges, preferring
///    to remove the triangle with the worst normal-gradient alignment.
///    Uses a two-phase strategy: first safe removals (all edges stay >= 2),
///    then aggressive removals to eliminate remaining non-manifold edges.
fn fix_non_manifold(
    vertices: &[Vector3<f64>],
    indices: &mut Vec<u32>,
    node: &SdfNode,
) {
    fix_non_manifold_impl(indices, vertices, |centroid| node.gradient(centroid));
}

/// Core manifold-repair logic, parameterized over the gradient function.
fn fix_non_manifold_impl(
    indices: &mut Vec<u32>,
    vertices: &[Vector3<f64>],
    grad_fn: impl Fn(Vector3<f64>) -> Vector3<f64>,
) {
    let ntri = indices.len() / 3;
    if ntri == 0 {
        return;
    }

    // Helper: get the 3 position-based edge keys for triangle ti.
    let tri_edges = |ti: usize, idx: &[u32]| -> [PosEdgeKey; 3] {
        let base = ti * 3;
        let v0 = &vertices[idx[base] as usize];
        let v1 = &vertices[idx[base + 1] as usize];
        let v2 = &vertices[idx[base + 2] as usize];
        [
            pos_edge_key(v0, v1),
            pos_edge_key(v1, v2),
            pos_edge_key(v2, v0),
        ]
    };

    // Score each triangle: normal alignment with SDF gradient.
    // Lower score = worse = remove first.
    let tri_scores: Vec<f64> = (0..ntri)
        .map(|ti| {
            let base = ti * 3;
            let v0 = vertices[indices[base] as usize];
            let v1 = vertices[indices[base + 1] as usize];
            let v2 = vertices[indices[base + 2] as usize];
            // Check for degenerate (coincident vertices).
            let q0 = quant_vertex(&v0);
            let q1 = quant_vertex(&v1);
            let q2 = quant_vertex(&v2);
            if q0 == q1 || q1 == q2 || q0 == q2 {
                return -2.0; // degenerate — absolutely worst
            }
            let face_normal = (v1 - v0).cross(&(v2 - v0));
            let area2 = face_normal.norm_squared();
            if area2 < 1e-30 {
                return -1.0;
            }
            let fn_norm = face_normal / area2.sqrt();
            let centroid = (v0 + v1 + v2) / 3.0;
            let grad = grad_fn(centroid);
            let grad_len = grad.norm();
            if grad_len < 1e-12 {
                return 0.0;
            }
            fn_norm.dot(&grad) / grad_len
        })
        .collect();

    let mut remove = vec![false; ntri];

    // Phase 1: Remove degenerate triangles.
    for ti in 0..ntri {
        if tri_scores[ti] < -1.5 {
            remove[ti] = true;
        }
    }

    // Phase 2: Iteratively remove excess triangles at non-manifold edges.
    // Strategy: build a global priority queue of "removable" triangles
    // sorted by score. For each candidate, check if removing it is safe
    // (all its edges still have >= 2 users). If not safe, try aggressive
    // removal (accept creating boundary edges if it reduces total
    // non-manifold count).
    for _pass in 0..30 {
        let mut changed = false;

        // Build live edge counts.
        let mut edge_count: HashMap<PosEdgeKey, usize> = HashMap::new();
        for ti in 0..ntri {
            if remove[ti] { continue; }
            for ek in &tri_edges(ti, indices) {
                *edge_count.entry(*ek).or_insert(0) += 1;
            }
        }

        // Collect all triangles that touch a non-manifold edge, scored.
        let mut candidates: Vec<(f64, usize)> = Vec::new();
        let mut seen_tris = std::collections::HashSet::new();
        for ti in 0..ntri {
            if remove[ti] { continue; }
            let edges = tri_edges(ti, indices);
            let touches_nm = edges.iter().any(|ek| {
                edge_count.get(ek).copied().unwrap_or(0) > 2
            });
            if touches_nm && seen_tris.insert(ti) {
                candidates.push((tri_scores[ti], ti));
            }
        }

        // Sort worst-first.
        candidates.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        for &(_, ti) in &candidates {
            if remove[ti] { continue; }

            let edges = tri_edges(ti, indices);

            // Check if this triangle still touches a non-manifold edge.
            let still_nm = edges.iter().any(|ek| {
                edge_count.get(ek).copied().unwrap_or(0) > 2
            });
            if !still_nm { continue; }

            // Safe check: will removing this triangle create any boundary edges?
            let safe = edges.iter().all(|ek| {
                edge_count.get(ek).copied().unwrap_or(0) > 2
            });

            if safe {
                remove[ti] = true;
                changed = true;
                for ek in &edges {
                    if let Some(c) = edge_count.get_mut(ek) { *c -= 1; }
                }
            }
        }

        if !changed { break; }
    }

    // Phase 3: Aggressive removal for remaining non-manifold edges.
    // Some non-manifold configurations can't be fixed with safe-only removal.
    // For each remaining non-manifold edge, remove the worst triangle even
    // if it creates boundary edges. This may leave small holes but ensures
    // no edge has >2 triangles.
    for _pass in 0..30 {
        let mut changed = false;

        let mut edge_count: HashMap<PosEdgeKey, usize> = HashMap::new();
        for ti in 0..ntri {
            if remove[ti] { continue; }
            for ek in &tri_edges(ti, indices) {
                *edge_count.entry(*ek).or_insert(0) += 1;
            }
        }

        // Build edge → triangles for non-manifold edges.
        let mut nm_edge_tris: HashMap<PosEdgeKey, Vec<usize>> = HashMap::new();
        for ti in 0..ntri {
            if remove[ti] { continue; }
            for ek in &tri_edges(ti, indices) {
                if edge_count.get(ek).copied().unwrap_or(0) > 2 {
                    nm_edge_tris.entry(*ek).or_default().push(ti);
                }
            }
        }

        for (_ek, tris) in &nm_edge_tris {
            if tris.len() <= 2 { continue; }

            // Sort worst-first.
            let mut scored: Vec<(f64, usize)> =
                tris.iter().map(|&ti| (tri_scores[ti], ti)).collect();
            scored.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

            // Remove worst until only 2 remain.
            let to_remove = scored.len() - 2;
            for &(_, ti) in scored.iter().take(to_remove) {
                if !remove[ti] {
                    remove[ti] = true;
                    changed = true;
                    for ek in &tri_edges(ti, indices) {
                        if let Some(c) = edge_count.get_mut(ek) { *c -= 1; }
                    }
                }
            }
        }

        if !changed { break; }
    }

    // Phase 4: Heal boundary edges left by aggressive removal.
    // If a boundary edge (count=1) has a neighbor triangle that could "close"
    // the hole, we skip this for now — the boundary holes at CSG seams are
    // typically very small and closing them requires inserting new triangles,
    // which is beyond simple removal.
    //
    // Instead, we attempt to re-add removed triangles that would close
    // boundary edges without re-introducing non-manifold edges.
    {
        let mut edge_count: HashMap<PosEdgeKey, usize> = HashMap::new();
        for ti in 0..ntri {
            if remove[ti] { continue; }
            for ek in &tri_edges(ti, indices) {
                *edge_count.entry(*ek).or_insert(0) += 1;
            }
        }

        // Find boundary edges.
        let boundary_edges: std::collections::HashSet<PosEdgeKey> = edge_count
            .iter()
            .filter(|&(_, &c)| c == 1)
            .map(|(ek, _)| *ek)
            .collect();

        if !boundary_edges.is_empty() {
            // Try to re-add removed triangles that border boundary edges
            // and wouldn't re-create non-manifold edges.
            // Sort by score descending (best first).
            let mut re_add_candidates: Vec<(f64, usize)> = (0..ntri)
                .filter(|&ti| remove[ti] && tri_scores[ti] > -1.5)
                .map(|ti| (tri_scores[ti], ti))
                .collect();
            re_add_candidates.sort_by(|a, b| {
                b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal)
            });

            for &(_, ti) in &re_add_candidates {
                let edges = tri_edges(ti, indices);
                // Only re-add if at least one of its edges is a boundary edge.
                let touches_boundary = edges.iter().any(|ek| {
                    edge_count.get(ek).copied().unwrap_or(0) == 1
                });
                // And re-adding wouldn't make any edge > 2.
                let safe = edges.iter().all(|ek| {
                    edge_count.get(ek).copied().unwrap_or(0) < 2
                });
                if touches_boundary && safe {
                    remove[ti] = false;
                    for ek in &edges {
                        *edge_count.entry(*ek).or_insert(0) += 1;
                    }
                }
            }
        }
    }

    // Rebuild index buffer without removed triangles.
    let new_indices: Vec<u32> = (0..ntri)
        .filter(|&ti| !remove[ti])
        .flat_map(|ti| {
            let base = ti * 3;
            [indices[base], indices[base + 1], indices[base + 2]]
        })
        .collect();

    *indices = new_indices;
}

/// Same as fix_non_manifold but for &dyn Sdf meshes.
fn fix_non_manifold_sdf(
    vertices: &[Vector3<f64>],
    indices: &mut Vec<u32>,
    sdf: &dyn Sdf,
) {
    fix_non_manifold_impl(indices, vertices, |centroid| {
        sdf.gradient(centroid)
            .unwrap_or_else(|| central_diff_gradient_sdf(sdf, centroid))
    });
}

// ---------------------------------------------------------------------------
// Octree construction from &dyn Sdf (no interval pruning)
// ---------------------------------------------------------------------------

/// Lightweight wrapper to hold an octree built from &dyn Sdf.
struct OctreeWrapper {
    root: OctreeCell,
}

impl OctreeWrapper {
    fn surface_cells(&self) -> Vec<&OctreeCell> {
        let mut result = Vec::new();
        Self::collect_surface_cells(&self.root, &mut result);
        result
    }

    fn collect_all_leaves(&self) -> Vec<&OctreeCell> {
        let mut result = Vec::new();
        Self::gather_leaves(&self.root, &mut result);
        result
    }

    fn collect_surface_cells<'a>(cell: &'a OctreeCell, result: &mut Vec<&'a OctreeCell>) {
        if cell.is_leaf() {
            if cell.has_sign_change() {
                result.push(cell);
            }
        } else {
            for child in cell.children.as_ref().unwrap().iter() {
                Self::collect_surface_cells(child, result);
            }
        }
    }

    fn gather_leaves<'a>(cell: &'a OctreeCell, result: &mut Vec<&'a OctreeCell>) {
        if cell.is_leaf() {
            result.push(cell);
        } else {
            for child in cell.children.as_ref().unwrap().iter() {
                Self::gather_leaves(child, result);
            }
        }
    }
}

/// Recursively build an octree cell from a &dyn Sdf trait object.
/// No interval pruning (we don't have interval_evaluate on the trait).
fn build_cell_from_sdf(
    sdf: &dyn Sdf,
    bbox: &BBox3,
    depth: u8,
    settings: &MeshSettings,
) -> OctreeCell {
    let bbox_corners = bbox.corners();
    let mut corners = [0.0f64; 8];
    for (i, c) in bbox_corners.iter().enumerate() {
        corners[i] = sdf.evaluate(*c);
    }

    let has_sign_change = {
        let first_sign = corners[0] > 0.0;
        corners.iter().any(|&c| (c > 0.0) != first_sign)
    };

    // When all corners share the same sign but one is very close to zero,
    // the surface may pass through the cell in a narrow region that corners
    // miss (e.g., the acute inner edge of a Union).  Force subdivision so
    // finer cells can resolve the gap.
    let near_surface = if !has_sign_change && depth >= settings.min_depth {
        let cell_diag = (bbox.max - bbox.min).norm();
        let closest_to_zero = corners.iter().map(|c| c.abs()).fold(f64::MAX, f64::min);
        closest_to_zero < cell_diag
    } else {
        false
    };

    // No interval pruning for trait objects. Subdivide if sign change,
    // below min_depth, or near-surface heuristic fires.
    if depth < settings.max_depth && (has_sign_change || depth < settings.min_depth || near_surface)
    {
        let octants = bbox.octants();
        let children: [OctreeCell; 8] =
            std::array::from_fn(|i| build_cell_from_sdf(sdf, &octants[i], depth + 1, settings));

        return OctreeCell {
            bbox: *bbox,
            depth,
            corners,
            children: Some(Box::new(children)),
        };
    }

    OctreeCell {
        bbox: *bbox,
        depth,
        corners,
        children: None,
    }
}
