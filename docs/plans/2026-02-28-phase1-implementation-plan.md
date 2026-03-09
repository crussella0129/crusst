# Phase 1: Manifold Mesh Pipeline — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make the dual contouring pipeline produce strictly manifold, watertight meshes with consistent outward-facing winding. Zero non-manifold edges. V - E + F = 2 for genus-0 shapes.

**Architecture:** Four sequential changes to the existing meshing pipeline: (1) increase edge key quantization precision, (2) add 2:1 octree balancing to eliminate T-junctions, (3) move bbox padding into the mesher for universal watertight guarantee, (4) add gradient-based winding correction as a post-process.

**Tech Stack:** Rust, nalgebra 0.33, existing crusst crate modules (`octree.rs`, `dual_contouring.rs`, `builder.rs`, `types.rs`). Tests use `cargo test`.

---

## Task 1: Edge Key Precision — Extract Constant and Increase to 30-bit

**Files:**
- Modify: `src/dual_contouring.rs:291` (cell_key), `src/dual_contouring.rs:319` (edge_key), `src/dual_contouring.rs:535` (sort_vertices_around_edge decode — 2 locations)

**Step 1: Add the module-level constant**

At the top of `src/dual_contouring.rs`, after the `use` block (after line 23), add:

```rust
/// Quantization scale for discretizing floating-point coordinates to integer
/// grid keys. 2^30 gives ~9.3e-10 resolution per unit — sufficient for
/// depths up to 30 without collisions.
const QUANT_SCALE: f64 = (1u64 << 30) as f64;
```

**Step 2: Replace all 4 usage sites**

In `cell_key` (line 291), replace:
```rust
    let scale = (1u64 << 20) as f64;
```
with:
```rust
    let scale = QUANT_SCALE;
```

In `edge_key` (line 319), replace:
```rust
    let scale = (1u64 << 20) as f64;
```
with:
```rust
    let scale = QUANT_SCALE;
```

In `sort_vertices_around_edge` (line 535), replace:
```rust
    let scale = (1u64 << 20) as f64;
```
with:
```rust
    let scale = QUANT_SCALE;
```

**Step 3: Run existing tests to verify no regression**

Run: `cargo test --lib --test dc_test --test octree_test -- --nocapture 2>&1 | tail -20`
Expected: All existing tests pass. The quantization change is transparent to consumers.

**Step 4: Commit**

```bash
git add src/dual_contouring.rs
git commit -m "refactor: extract QUANT_SCALE constant, increase to 30-bit precision

Edge keys and cell keys now use 2^30 quantization (was 2^20).
Resolution improves from 9.5e-7 to 9.3e-10 per unit, eliminating
potential collisions at octree depth 8+."
```

---

## Task 2: Octree 2:1 Balancing — Neighbor Finding

**Files:**
- Modify: `src/octree.rs` (add helper methods to `OctreeCell` and `Octree`)

This task adds the two building blocks for balancing: finding the leaf cell that contains a point, and collecting all leaf cells.

**Step 1: Write the failing test for find_leaf_containing**

Add to `tests/octree_test.rs`:

```rust
#[test]
fn octree_find_leaf_containing_returns_correct_cell() {
    use crusst::dag::SdfNode;
    use crusst::octree::Octree;
    use crusst::types::{BBox3, MeshSettings};
    use nalgebra::Vector3;

    let node = SdfNode::Sphere {
        center: Vector3::zeros(),
        radius: 5.0,
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let settings = MeshSettings {
        max_depth: 4,
        min_depth: 2,
        edge_tolerance: 1e-6,
    };
    let tree = Octree::build(&node, &bbox, &settings);

    // A point near the surface should land in a deep cell
    let probe = Vector3::new(4.9, 0.0, 0.0);
    let leaf = tree.find_leaf_containing(probe);
    assert!(leaf.is_some());
    let cell = leaf.unwrap();
    assert!(cell.bbox.contains(probe));
    // Near-surface cells should be at or near max_depth
    assert!(cell.depth >= 3);

    // A point far outside should land in a shallow cell (pruned early)
    let far = Vector3::new(6.5, 6.5, 6.5);
    let leaf_far = tree.find_leaf_containing(far);
    assert!(leaf_far.is_some());
    assert!(leaf_far.unwrap().bbox.contains(far));
}
```

**Step 2: Run test to verify it fails**

Run: `cargo test octree_find_leaf_containing -- --nocapture 2>&1 | tail -5`
Expected: FAIL — `find_leaf_containing` does not exist.

**Step 3: Implement find_leaf_containing and collect_all_leaves**

Add to `src/octree.rs`, inside the `impl Octree` block (after `collect_surface_cells`, around line 195):

```rust
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

    /// Collect all leaf cells in the tree with their bounding boxes.
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
```

Add the missing import at the top of `src/octree.rs` (after line 2):
```rust
use nalgebra::Vector3;
```

**Step 4: Run test to verify it passes**

Run: `cargo test octree_find_leaf_containing -- --nocapture 2>&1 | tail -5`
Expected: PASS

**Step 5: Commit**

```bash
git add src/octree.rs tests/octree_test.rs
git commit -m "feat: add find_leaf_containing and collect_all_leaves to Octree

Building blocks for 2:1 balance enforcement. find_leaf_containing
walks from root to the leaf cell containing a probe point (O(depth)).
collect_all_leaves gathers every leaf for the violation scan."
```

---

## Task 3: Octree 2:1 Balancing — Balance Pass

**Files:**
- Modify: `src/octree.rs` (add `balance` method and `subdivide_for_balance` helper)

**Step 1: Write the failing test for balance**

Add to `tests/octree_test.rs`:

```rust
#[test]
fn octree_balance_enforces_two_to_one() {
    use crusst::dag::SdfNode;
    use crusst::octree::Octree;
    use crusst::types::{BBox3, MeshSettings};
    use nalgebra::Vector3;

    let node = SdfNode::Sphere {
        center: Vector3::zeros(),
        radius: 5.0,
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let settings = MeshSettings {
        max_depth: 6,
        min_depth: 2,
        edge_tolerance: 1e-6,
    };
    let mut tree = Octree::build(&node, &bbox, &settings);
    tree.balance(&node);

    // Verify: for every leaf cell, all 6 face-neighbors differ by at most 1 depth level
    let leaves = tree.collect_all_leaves();
    for cell in &leaves {
        let size = cell.bbox.size();
        let center = cell.bbox.center();
        // Check 6 face neighbors (±x, ±y, ±z)
        let probes = [
            center + Vector3::new(size.x, 0.0, 0.0),
            center - Vector3::new(size.x, 0.0, 0.0),
            center + Vector3::new(0.0, size.y, 0.0),
            center - Vector3::new(0.0, size.y, 0.0),
            center + Vector3::new(0.0, 0.0, size.z),
            center - Vector3::new(0.0, 0.0, size.z),
        ];
        for probe in &probes {
            if let Some(neighbor) = tree.find_leaf_containing(*probe) {
                let depth_diff = (cell.depth as i16 - neighbor.depth as i16).abs();
                assert!(
                    depth_diff <= 1,
                    "Balance violation: cell at depth {} (center {:?}) has neighbor at depth {} (center {:?})",
                    cell.depth,
                    center,
                    neighbor.depth,
                    neighbor.bbox.center(),
                );
            }
            // If probe is outside root bbox, no neighbor — that's fine
        }
    }
}
```

**Step 2: Run test to verify it fails**

Run: `cargo test octree_balance_enforces_two_to_one -- --nocapture 2>&1 | tail -5`
Expected: FAIL — `balance` does not exist.

**Step 3: Implement the balance method**

Add to `src/octree.rs`, inside the `impl Octree` block:

```rust
    /// Enforce 2:1 balance constraint on the octree.
    ///
    /// After this call, no leaf cell has a face-neighbor that differs by more
    /// than one depth level. This eliminates T-junctions in the dual contouring
    /// mesh. Uses a work-queue approach: immutable scan finds violations, then
    /// mutable pass subdivides. Repeats until no violations remain.
    pub fn balance(&mut self, node: &SdfNode) {
        loop {
            // Immutable pass: find all leaf cells that need a coarser neighbor subdivided.
            let violations = self.find_balance_violations();
            if violations.is_empty() {
                break;
            }
            // Mutable pass: subdivide each violation target.
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
            let size = cell.bbox.size();
            let center = cell.bbox.center();
            // Check 6 face neighbors
            let probes = [
                center + Vector3::new(size.x, 0.0, 0.0),
                center - Vector3::new(size.x, 0.0, 0.0),
                center + Vector3::new(0.0, size.y, 0.0),
                center - Vector3::new(0.0, size.y, 0.0),
                center + Vector3::new(0.0, 0.0, size.z),
                center - Vector3::new(0.0, 0.0, size.z),
            ];
            for probe in &probes {
                if let Some(neighbor) = self.find_leaf_containing(*probe) {
                    if neighbor.depth + 1 < cell.depth {
                        // neighbor is too coarse — needs subdivision
                        // Use quantized bbox min as dedup key
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
                // Subdivide: create 8 leaf children with evaluated corners
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
            // Not the target cell — leave it alone
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
```

Add the `std::collections::HashSet` import is already available via the existing `use std::collections::HashMap` pattern — Rust resolves it from the qualified path in the function body.

**Step 4: Run test to verify it passes**

Run: `cargo test octree_balance_enforces_two_to_one -- --nocapture 2>&1 | tail -10`
Expected: PASS

**Step 5: Run ALL existing tests to verify no regression**

Run: `cargo test 2>&1 | tail -5`
Expected: All tests pass.

**Step 6: Commit**

```bash
git add src/octree.rs tests/octree_test.rs
git commit -m "feat: add 2:1 octree balance pass to eliminate T-junctions

Octree::balance() enforces that no leaf cell has a face-neighbor
differing by more than one depth level. Uses a work-queue approach:
immutable violation scan, then mutable forced subdivision, looping
until convergence. Forced subdivision evaluates corners only — no
interval pruning or recursive build logic."
```

---

## Task 4: Integrate Balancing into DC Pipeline

**Files:**
- Modify: `src/dual_contouring.rs:37-38` (call balance after build in `extract_mesh_adaptive`)

**Step 1: Add balance call to extract_mesh_adaptive**

In `src/dual_contouring.rs`, replace lines 37-38:

```rust
    // Phase 1: Build the adaptive octree.
    let tree = Octree::build(node, bbox, settings);
```

with:

```rust
    // Phase 1: Build the adaptive octree and enforce 2:1 balance.
    let mut tree = Octree::build(node, bbox, settings);
    tree.balance(node);
```

**Step 2: Run all tests**

Run: `cargo test 2>&1 | tail -5`
Expected: All tests pass. The balanced octree should produce identical or better meshes.

**Step 3: Commit**

```bash
git add src/dual_contouring.rs
git commit -m "feat: integrate octree balancing into DC pipeline

extract_mesh_adaptive now calls tree.balance(node) after build,
ensuring the octree satisfies 2:1 face-neighbor constraint before
face generation. This eliminates T-junctions in the output mesh."
```

---

## Task 5: Watertight Bbox Padding

**Files:**
- Modify: `src/dual_contouring.rs:36-38` (add padding in `extract_mesh_adaptive`)
- Modify: `src/builder.rs:770-776` (remove redundant padding from `Shape::mesh()`)

**Step 1: Add smart padding to extract_mesh_adaptive**

In `src/dual_contouring.rs`, after the function signature (line 36) and before the octree build, add padding logic. Replace the block at lines 36-39:

```rust
) -> TriangleMesh {
    // Phase 1: Build the adaptive octree and enforce 2:1 balance.
    let mut tree = Octree::build(node, bbox, settings);
    tree.balance(node);
```

with:

```rust
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
```

**Step 2: Remove redundant padding from Shape::mesh()**

In `src/builder.rs`, replace lines 770-776:

```rust
    pub fn mesh(&self, settings: MeshSettings) -> TriangleMesh {
        let bbox = self.bounding_box();
        // Expand bbox slightly to avoid clipping surface cells at the boundary.
        let pad = bbox.size() * 0.05;
        let padded = BBox3::new(bbox.min - pad, bbox.max + pad);
        extract_mesh_adaptive(&self.node, &padded, &settings)
    }
```

with:

```rust
    pub fn mesh(&self, settings: MeshSettings) -> TriangleMesh {
        let bbox = self.bounding_box();
        extract_mesh_adaptive(&self.node, &bbox, &settings)
    }
```

**Step 3: Run all tests**

Run: `cargo test 2>&1 | tail -5`
Expected: All tests pass.

**Step 4: Commit**

```bash
git add src/dual_contouring.rs src/builder.rs
git commit -m "feat: move bbox padding into extract_mesh_adaptive for watertight guarantee

All callers now get smart padding: max(5% of bbox, 2 cells at finest
resolution). Removes redundant 5% pad from Shape::mesh() to avoid
double-padding. Ensures the surface never touches the domain boundary."
```

---

## Task 6: Gradient-Based Winding Fix

**Files:**
- Modify: `src/dual_contouring.rs` (add `fix_winding` function, call it before `split_sharp_edges`)

**Step 1: Add the fix_winding function**

Add after the `emit_fan` function (after line 622 in the current file — adjust for earlier edits) in `src/dual_contouring.rs`:

```rust
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
```

**Step 2: Call fix_winding in extract_mesh_adaptive**

In `extract_mesh_adaptive`, between the face generation loop and the normals computation, add the winding fix call. Find the line that says:

```rust
    // Phase 4: Compute per-vertex normals.
    let normals: Vec<Vector3<f64>> = vertices.iter().map(|v| node.gradient(*v)).collect();
```

Insert before it:

```rust
    // Phase 3.5: Fix winding order using SDF gradient.
    fix_winding(&vertices, &mut indices, node);

```

**Step 3: Also add fix_winding to extract_mesh_from_sdf**

In `extract_mesh_from_sdf`, find the similar normals computation block:

```rust
    // Phase 4: Normals — prefer analytical gradients, fall back to central differences.
    let normals: Vec<Vector3<f64>> = vertices
```

Insert before it:

```rust
    // Phase 3.5: Fix winding order using SDF gradient (central differences).
    fix_winding_sdf(&vertices, &mut indices, sdf);

```

And add the companion function for `&dyn Sdf`:

```rust
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
```

**Step 4: Run all tests**

Run: `cargo test 2>&1 | tail -5`
Expected: All tests pass.

**Step 5: Commit**

```bash
git add src/dual_contouring.rs
git commit -m "feat: add gradient-based winding fix for consistent outward normals

fix_winding() compares each triangle's face normal against the SDF
gradient at the centroid. If they disagree, the triangle is flipped.
Applied as a post-process after face generation, before split_sharp_edges.
Much more robust than the previous edge-key-search approach."
```

---

## Task 7: Strict Manifold Tests

**Files:**
- Modify: `tests/dc_test.rs` (add `assert_manifold` helper and 6 strict tests)

**Step 1: Write the assert_manifold helper and all 6 test functions**

Replace the entire contents of `tests/dc_test.rs` with:

```rust
use crusst::dag::SdfNode;
use crusst::dual_contouring::extract_mesh_adaptive;
use crusst::types::{BBox3, MeshSettings, TriangleMesh};
use nalgebra::Vector3;
use std::collections::HashMap;
use std::sync::Arc;

// ---------------------------------------------------------------------------
// Manifold assertion helper
// ---------------------------------------------------------------------------

/// Assert the mesh is strictly manifold and watertight:
/// 1. Every undirected edge is shared by exactly 2 triangles.
/// 2. Euler characteristic V - E + F = 2 (genus-0 closed surface).
fn assert_manifold(mesh: &TriangleMesh, label: &str) {
    assert!(
        !mesh.indices.is_empty(),
        "{}: mesh has no triangles",
        label
    );
    assert_eq!(
        mesh.indices.len() % 3,
        0,
        "{}: index count not divisible by 3",
        label
    );

    // Count how many triangles share each undirected edge
    let mut edge_count: HashMap<(u32, u32), usize> = HashMap::new();
    for tri in mesh.indices.chunks(3) {
        for i in 0..3 {
            let a = tri[i];
            let b = tri[(i + 1) % 3];
            let key = (a.min(b), a.max(b));
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }

    let non_manifold: Vec<_> = edge_count
        .iter()
        .filter(|(_, &count)| count != 2)
        .collect();

    assert!(
        non_manifold.is_empty(),
        "{}: {} non-manifold edges out of {} total (expected all edges shared by exactly 2 triangles). \
         Samples: {:?}",
        label,
        non_manifold.len(),
        edge_count.len(),
        non_manifold.iter().take(5).collect::<Vec<_>>(),
    );

    // Euler characteristic: V - E + F = 2 for genus-0
    let v = mesh.vertices.len() as i64;
    let e = edge_count.len() as i64;
    let f = (mesh.indices.len() / 3) as i64;
    let euler = v - e + f;
    assert_eq!(
        euler, 2,
        "{}: Euler characteristic V - E + F = {} (expected 2). V={}, E={}, F={}",
        label, euler, v, e, f
    );
}

/// Assert all triangle normals point outward (agree with SDF gradient).
fn assert_outward_winding(mesh: &TriangleMesh, node: &SdfNode, label: &str) {
    let mut flipped = 0usize;
    let total = mesh.indices.len() / 3;
    for tri in mesh.indices.chunks(3) {
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];
        let face_normal = (v1 - v0).cross(&(v2 - v0));
        if face_normal.norm_squared() < 1e-30 {
            continue;
        }
        let centroid = (v0 + v1 + v2) / 3.0;
        let gradient = node.gradient(centroid);
        if face_normal.dot(&gradient) < 0.0 {
            flipped += 1;
        }
    }
    assert_eq!(
        flipped, 0,
        "{}: {}/{} triangles have wrong winding (face normal disagrees with SDF gradient)",
        label, flipped, total
    );
}

// ---------------------------------------------------------------------------
// Helper to build a mesh from an SdfNode
// ---------------------------------------------------------------------------

fn mesh_node(node: &SdfNode, bbox: &BBox3, max_depth: u8) -> TriangleMesh {
    let settings = MeshSettings {
        max_depth,
        min_depth: (max_depth / 2).max(2),
        edge_tolerance: 1e-6,
    };
    extract_mesh_adaptive(node, bbox, &settings)
}

// ---------------------------------------------------------------------------
// Existing tests (preserved)
// ---------------------------------------------------------------------------

#[test]
fn dc_sphere_produces_manifold_mesh() {
    let node = SdfNode::Sphere {
        center: Vector3::zeros(),
        radius: 5.0,
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let settings = MeshSettings {
        max_depth: 5,
        min_depth: 3,
        edge_tolerance: 1e-6,
    };
    let mesh = extract_mesh_adaptive(&node, &bbox, &settings);
    assert!(mesh.vertices.len() > 10);
    assert!(mesh.indices.len() > 10);
    assert_eq!(mesh.indices.len() % 3, 0);
}

#[test]
fn dc_box_has_sharp_edges() {
    let node = SdfNode::Box3 {
        center: Vector3::zeros(),
        half_extents: Vector3::new(5.0, 5.0, 5.0),
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let settings = MeshSettings {
        max_depth: 6,
        min_depth: 3,
        edge_tolerance: 1e-6,
    };
    let mesh = extract_mesh_adaptive(&node, &bbox, &settings);

    let corner = Vector3::new(5.0, 5.0, 5.0);
    let min_dist = mesh
        .vertices
        .iter()
        .map(|v| (v - corner).norm())
        .fold(f64::INFINITY, f64::min);
    let cell_size = 14.0 / 64.0;
    assert!(
        min_dist < cell_size * 2.0,
        "Nearest vertex to corner is {:.4}, expected < {:.4}",
        min_dist,
        cell_size * 2.0
    );
}

#[test]
fn dc_higher_depth_produces_more_triangles() {
    let node = SdfNode::Sphere {
        center: Vector3::zeros(),
        radius: 5.0,
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let s4 = MeshSettings {
        max_depth: 4,
        min_depth: 2,
        edge_tolerance: 1e-6,
    };
    let s6 = MeshSettings {
        max_depth: 6,
        min_depth: 2,
        edge_tolerance: 1e-6,
    };
    let m4 = extract_mesh_adaptive(&node, &bbox, &s4);
    let m6 = extract_mesh_adaptive(&node, &bbox, &s6);
    assert!(m6.indices.len() > m4.indices.len());
}

// ---------------------------------------------------------------------------
// Strict manifold tests — Phase 1 additions
// ---------------------------------------------------------------------------

#[test]
fn manifold_sphere() {
    let node = SdfNode::Sphere {
        center: Vector3::zeros(),
        radius: 5.0,
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let mesh = mesh_node(&node, &bbox, 6);
    assert_manifold(&mesh, "sphere");
    assert_outward_winding(&mesh, &node, "sphere");
}

#[test]
fn manifold_box() {
    let node = SdfNode::Box3 {
        center: Vector3::zeros(),
        half_extents: Vector3::new(5.0, 5.0, 5.0),
    };
    let bbox = BBox3::new(Vector3::from_element(-7.0), Vector3::from_element(7.0));
    let mesh = mesh_node(&node, &bbox, 6);
    assert_manifold(&mesh, "box");
    assert_outward_winding(&mesh, &node, "box");
}

#[test]
fn manifold_cylinder() {
    let node = SdfNode::Cylinder {
        base: Vector3::new(0.0, -5.0, 0.0),
        axis: Vector3::new(0.0, 1.0, 0.0),
        radius: 3.0,
        height: 10.0,
    };
    let bbox = BBox3::new(Vector3::new(-5.0, -7.0, -5.0), Vector3::new(5.0, 7.0, 5.0));
    let mesh = mesh_node(&node, &bbox, 6);
    assert_manifold(&mesh, "cylinder");
    assert_outward_winding(&mesh, &node, "cylinder");
}

#[test]
fn manifold_union_two_spheres() {
    let a = Arc::new(SdfNode::Sphere {
        center: Vector3::new(-2.0, 0.0, 0.0),
        radius: 5.0,
    });
    let b = Arc::new(SdfNode::Sphere {
        center: Vector3::new(2.0, 0.0, 0.0),
        radius: 5.0,
    });
    let node = SdfNode::Union(a, b);
    let bbox = BBox3::new(Vector3::from_element(-9.0), Vector3::from_element(9.0));
    let mesh = mesh_node(&node, &bbox, 6);
    assert_manifold(&mesh, "union_spheres");
    assert_outward_winding(&mesh, &node, "union_spheres");
}

#[test]
fn manifold_smooth_union_two_spheres() {
    let a = Arc::new(SdfNode::Sphere {
        center: Vector3::new(-2.0, 0.0, 0.0),
        radius: 5.0,
    });
    let b = Arc::new(SdfNode::Sphere {
        center: Vector3::new(2.0, 0.0, 0.0),
        radius: 5.0,
    });
    let node = SdfNode::SmoothUnion(a, b, 1.5);
    let bbox = BBox3::new(Vector3::from_element(-10.0), Vector3::from_element(10.0));
    let mesh = mesh_node(&node, &bbox, 6);
    assert_manifold(&mesh, "smooth_union_spheres");
    assert_outward_winding(&mesh, &node, "smooth_union_spheres");
}

#[test]
fn manifold_thin_box() {
    let node = SdfNode::Box3 {
        center: Vector3::zeros(),
        half_extents: Vector3::new(1.0, 1.0, 0.01),
    };
    let bbox = BBox3::new(Vector3::new(-2.0, -2.0, -2.0), Vector3::new(2.0, 2.0, 2.0));
    let mesh = mesh_node(&node, &bbox, 8);
    assert_manifold(&mesh, "thin_box");
    assert_outward_winding(&mesh, &node, "thin_box");
}
```

**Step 2: Run just the new strict tests**

Run: `cargo test manifold_ -- --nocapture 2>&1 | tail -20`
Expected: All 6 `manifold_*` tests PASS.

**Step 3: Run the full test suite**

Run: `cargo test 2>&1 | tail -5`
Expected: All tests pass (new + existing).

**Step 4: Commit**

```bash
git add tests/dc_test.rs
git commit -m "test: add strict manifold + winding tests for 6 shape types

assert_manifold: every edge shared by exactly 2 triangles, Euler V-E+F=2.
assert_outward_winding: every triangle normal agrees with SDF gradient.
Test shapes: sphere, box, cylinder, union, smooth union, thin box."
```

---

## Task 8: Remove Old Relaxed Watertight Test

**Files:**
- Modify: `tests/dc_test.rs`

The old `dc_sphere_is_watertight` test uses a relaxed 95% threshold. It is now superseded by the strict `manifold_sphere` test.

**Step 1: Remove dc_sphere_is_watertight**

Delete the `dc_sphere_is_watertight` function from `tests/dc_test.rs` (it was the second test in the original file, now reproduced in the rewrite above — verify it is NOT present in the Task 7 rewrite above, which it is not).

This was already handled in Task 7 — the full rewrite dropped `dc_sphere_is_watertight` and replaced it with `manifold_sphere`.

**Step 2: Run full suite**

Run: `cargo test 2>&1 | tail -5`
Expected: All tests pass.

**Step 3: Commit (only if dc_sphere_is_watertight was not already removed in Task 7)**

If the Task 7 rewrite already removed it, this task is a no-op. Otherwise:

```bash
git add tests/dc_test.rs
git commit -m "test: remove relaxed 95% watertight test, superseded by strict manifold"
```

---

## Summary

| Task | What | Files | Estimated effort |
|------|------|-------|-----------------|
| 1 | Edge key precision | `dual_contouring.rs` | 3 min |
| 2 | Neighbor finding helpers | `octree.rs` | 10 min |
| 3 | 2:1 balance pass | `octree.rs` | 20 min |
| 4 | Integrate balancing into DC | `dual_contouring.rs` | 3 min |
| 5 | Watertight bbox padding | `dual_contouring.rs`, `builder.rs` | 5 min |
| 6 | Gradient winding fix | `dual_contouring.rs` | 10 min |
| 7 | Strict manifold tests | `tests/dc_test.rs` | 10 min |
| 8 | Remove old relaxed test | `tests/dc_test.rs` | 1 min |

Total: 8 tasks, 8 commits, ~60 min of implementation.
