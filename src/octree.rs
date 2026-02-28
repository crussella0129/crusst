//! Adaptive octree construction with interval arithmetic pruning.
//!
//! The octree recursively subdivides an axis-aligned bounding box, refining
//! only near the surface of the SDF. Cells that are entirely inside or
//! entirely outside (proved via `interval_evaluate`) are pruned early,
//! giving O(surface area) memory instead of O(volume).
//!
//! The resulting tree is consumed by dual contouring (Task 7) which calls
//! `surface_cells()` to find cells that need mesh vertices.

use crate::dag::SdfNode;
use crate::types::{BBox3, MeshSettings};
use nalgebra::Vector3;

/// A single cell in the adaptive octree.
pub struct OctreeCell {
    /// Axis-aligned bounding box for this cell.
    pub bbox: BBox3,
    /// Depth of this cell in the tree (root = 0).
    pub depth: u8,
    /// SDF values at the 8 corners of the bounding box.
    pub corners: [f64; 8],
    /// Children (8 octants) if this cell was subdivided, or `None` for leaves.
    pub children: Option<Box<[OctreeCell; 8]>>,
}

/// An adaptive octree built from a signed distance field.
pub struct Octree {
    /// The root cell spanning the entire bounding box.
    pub root: OctreeCell,
}

impl OctreeCell {
    /// Does this cell have sign changes among its 8 corner SDF values?
    ///
    /// A sign change indicates the surface passes through this cell.
    pub fn has_sign_change(&self) -> bool {
        let first_sign = self.corners[0] > 0.0;
        self.corners.iter().any(|&c| (c > 0.0) != first_sign)
    }

    /// Is this a leaf cell (no children)?
    pub fn is_leaf(&self) -> bool {
        self.children.is_none()
    }
}

impl Octree {
    /// Build an adaptive octree for the given SDF node within the bounding box.
    ///
    /// The tree is refined adaptively: cells near the surface are subdivided
    /// up to `settings.max_depth`, while cells proved entirely inside or
    /// outside by interval arithmetic are pruned at `settings.min_depth`.
    pub fn build(node: &SdfNode, bbox: &BBox3, settings: &MeshSettings) -> Self {
        // Evaluate SDF at the 8 corners of the root cell
        let bbox_corners = bbox.corners();
        let mut corners = [0.0f64; 8];
        for (i, c) in bbox_corners.iter().enumerate() {
            corners[i] = node.evaluate(*c);
        }

        // Recursively build the tree starting from depth 0
        let root = Self::build_cell(node, bbox, 0, corners, settings);
        Octree { root }
    }

    /// Recursively build a single octree cell.
    ///
    /// Decision logic:
    /// 1. If interval arithmetic proves the cell is entirely inside or outside
    ///    AND we have reached min_depth -> leaf (no surface).
    /// 2. If there are sign changes AND depth < max_depth -> subdivide into 8.
    /// 3. If depth >= max_depth -> leaf (surface cell if sign changes, else not).
    /// 4. If no sign changes AND depth >= min_depth -> leaf (no surface).
    fn build_cell(
        node: &SdfNode,
        bbox: &BBox3,
        depth: u8,
        corners: [f64; 8],
        settings: &MeshSettings,
    ) -> OctreeCell {
        let has_sign_change = {
            let first_sign = corners[0] > 0.0;
            corners.iter().any(|&c| (c > 0.0) != first_sign)
        };

        // Pruning via interval arithmetic: if the interval over this cell
        // is entirely positive (outside) or entirely negative (inside),
        // and we've reached min_depth, stop subdividing.
        //
        // When the interval spans zero but no corner shows a sign change,
        // the surface passes through the cell between corners (e.g., at an
        // acute concave Union edge).  Force subdivision so finer cells can
        // catch the sign change.
        let mut interval_spans_zero = false;
        if depth >= settings.min_depth && !has_sign_change {
            let iv = node.interval_evaluate(bbox);
            if iv.definitely_positive() || iv.definitely_negative() {
                return OctreeCell {
                    bbox: *bbox,
                    depth,
                    corners,
                    children: None,
                };
            }
            interval_spans_zero = true;
        }

        // Subdivide if there are sign changes (or we haven't reached min_depth,
        // or interval arithmetic detects a possible surface that corners missed)
        // and we haven't hit the maximum depth yet.
        if depth < settings.max_depth
            && (has_sign_change || depth < settings.min_depth || interval_spans_zero)
        {
            let octants = bbox.octants();
            let children: [OctreeCell; 8] = std::array::from_fn(|i| {
                let child_bbox = &octants[i];
                let child_corners_pts = child_bbox.corners();
                let mut child_corners = [0.0f64; 8];
                for (j, c) in child_corners_pts.iter().enumerate() {
                    child_corners[j] = node.evaluate(*c);
                }
                Self::build_cell(node, child_bbox, depth + 1, child_corners, settings)
            });

            return OctreeCell {
                bbox: *bbox,
                depth,
                corners,
                children: Some(Box::new(children)),
            };
        }

        // Leaf cell: either at max_depth or no sign change beyond min_depth
        OctreeCell {
            bbox: *bbox,
            depth,
            corners,
            children: None,
        }
    }

    /// Count all leaf cells in the tree.
    pub fn leaf_count(&self) -> usize {
        Self::count_leaves(&self.root)
    }

    /// Count leaf cells that have sign changes (surface-crossing cells).
    pub fn surface_cell_count(&self) -> usize {
        Self::count_surface_cells(&self.root)
    }

    /// Collect references to all leaf cells that have sign changes.
    pub fn surface_cells(&self) -> Vec<&OctreeCell> {
        let mut result = Vec::new();
        Self::collect_surface_cells(&self.root, &mut result);
        result
    }

    fn count_leaves(cell: &OctreeCell) -> usize {
        if cell.is_leaf() {
            1
        } else {
            cell.children
                .as_ref()
                .unwrap()
                .iter()
                .map(Self::count_leaves)
                .sum()
        }
    }

    fn count_surface_cells(cell: &OctreeCell) -> usize {
        if cell.is_leaf() {
            if cell.has_sign_change() { 1 } else { 0 }
        } else {
            cell.children
                .as_ref()
                .unwrap()
                .iter()
                .map(Self::count_surface_cells)
                .sum()
        }
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

    /// Find the leaf cell containing the given point.
    /// Returns `None` if the point is outside the root bounding box.
    pub fn find_leaf_containing(&self, point: Vector3<f64>) -> Option<&OctreeCell> {
        Self::find_leaf_in(&self.root, point)
    }

    fn find_leaf_in<'a>(cell: &'a OctreeCell, point: Vector3<f64>) -> Option<&'a OctreeCell> {
        if !cell.bbox.contains(point) {
            return None;
        }
        if cell.is_leaf() {
            return Some(cell);
        }
        // Descend into the child whose octant contains the point
        for child in cell.children.as_ref().unwrap().iter() {
            if child.bbox.contains(point) {
                return Self::find_leaf_in(child, point);
            }
        }
        // Fallback: point is exactly on an octant boundary
        Some(cell.children.as_ref().unwrap().first().unwrap())
    }

    /// Collect all leaf cells in the tree.
    pub fn collect_all_leaves(&self) -> Vec<&OctreeCell> {
        let mut result = Vec::new();
        Self::gather_leaves(&self.root, &mut result);
        result
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

    /// Enforce 2:1 balance constraint on the octree.
    ///
    /// After this call, no leaf cell has a face-neighbor that differs by more
    /// than one depth level. This eliminates T-junctions in the dual contouring
    /// mesh. Uses a work-queue approach: immutable scan finds violations, then
    /// mutable pass subdivides. Repeats until no violations remain.
    pub fn balance(&mut self, node: &SdfNode) {
        loop {
            let violations = self.find_balance_violations();
            if violations.is_empty() {
                break;
            }
            for target_bbox in violations {
                Self::subdivide_at(&mut self.root, node, &target_bbox);
            }
        }
    }

    /// Scan all leaves and return the bboxes of coarser neighbors that violate 2:1.
    fn find_balance_violations(&self) -> Vec<BBox3> {
        let leaves = self.collect_all_leaves();
        let mut targets = Vec::new();
        let mut seen = std::collections::HashSet::new();

        for cell in &leaves {
            let half = cell.bbox.size() * 0.5;
            let center = cell.bbox.center();
            // Step just past the face (half-width + small nudge) to land in the
            // immediate face-neighbor cell rather than overshooting into a more
            // distant cell.
            let eps = half.x * 1e-3;
            let dx = half.x + eps;
            let dy = half.y + eps;
            let dz = half.z + eps;
            let probes = [
                center + Vector3::new(dx, 0.0, 0.0),
                center - Vector3::new(dx, 0.0, 0.0),
                center + Vector3::new(0.0, dy, 0.0),
                center - Vector3::new(0.0, dy, 0.0),
                center + Vector3::new(0.0, 0.0, dz),
                center - Vector3::new(0.0, 0.0, dz),
            ];
            for probe in &probes {
                if let Some(neighbor) = self.find_leaf_containing(*probe) {
                    if neighbor.depth + 1 < cell.depth {
                        // neighbor is too coarse — needs subdivision
                        let key = (
                            (neighbor.bbox.min.x * 1e9).round() as i64,
                            (neighbor.bbox.min.y * 1e9).round() as i64,
                            (neighbor.bbox.min.z * 1e9).round() as i64,
                        );
                        if seen.insert(key) {
                            targets.push(neighbor.bbox);
                        }
                    }
                }
            }
        }
        targets
    }

    /// Subdivide the leaf cell at the given bbox for balance purposes.
    /// This is NOT the build subdivision — no interval pruning, no recursion.
    /// Just evaluates corners and creates 8 leaf children.
    fn subdivide_at(cell: &mut OctreeCell, node: &SdfNode, target: &BBox3) {
        if cell.is_leaf() {
            // Is this the cell we're looking for?
            let eps = cell.bbox.size().x * 1e-6;
            if (cell.bbox.min - target.min).norm() < eps
                && (cell.bbox.max - target.max).norm() < eps
            {
                let octants = cell.bbox.octants();
                let children: [OctreeCell; 8] = std::array::from_fn(|i| {
                    let child_bbox = octants[i];
                    let corner_pts = child_bbox.corners();
                    let mut corners = [0.0f64; 8];
                    for (j, c) in corner_pts.iter().enumerate() {
                        corners[j] = node.evaluate(*c);
                    }
                    OctreeCell {
                        bbox: child_bbox,
                        depth: cell.depth + 1,
                        corners,
                        children: None,
                    }
                });
                cell.children = Some(Box::new(children));
            }
            return;
        }
        // Recurse into children to find the target
        if let Some(ref mut children) = cell.children {
            for child in children.iter_mut() {
                if child.bbox.contains(target.center()) {
                    Self::subdivide_at(child, node, target);
                    return;
                }
            }
        }
    }
}
