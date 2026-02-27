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
        }

        // Subdivide if there are sign changes (or we haven't reached min_depth)
        // and we haven't hit the maximum depth yet.
        if depth < settings.max_depth && (has_sign_change || depth < settings.min_depth) {
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
}
