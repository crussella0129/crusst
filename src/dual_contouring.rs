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
    // Phase 1: Build the adaptive octree.
    let tree = Octree::build(node, bbox, settings);

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
            let ek = edge_key(&corners[ci0], &corners[ci1]);
            edge_cells.entry(ek).or_default().push((key, *cell));
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

    // Phase 3: Generate faces (same as above).
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
            let ek = edge_key(&corners[ci0], &corners[ci1]);
            edge_cells.entry(ek).or_default().push((key, *cell));
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
    // is the root bbox size. We multiply by 2^20 to get integer precision
    // that works for depths up to 20.
    let scale = (1u64 << 20) as f64;
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
    let scale = (1u64 << 20) as f64;
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
    let scale = (1u64 << 20) as f64;
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
