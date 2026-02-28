# Phase 1: Manifold Mesh Pipeline

**Date:** 2026-02-28
**Status:** Approved
**Goal:** Produce strictly manifold, watertight meshes with consistent winding from the dual contouring pipeline. Zero non-manifold edges. V - E + F = 2 for genus-0 shapes.

---

## Implementation Order

```
1.4 Edge key precision       (const extraction, 4 locations)
  → 1.1 Octree 2:1 balancing (structural foundation)
    → 1.2 Watertight bbox     (padding in extract_mesh_adaptive)
      → 1.3 Winding fix       (gradient-based post-process)
        → Tests               (strict manifold + Euler validation)
```

---

## 1.4 Edge Key Precision

**Files:** `src/dual_contouring.rs`

- Extract `(1u64 << 20)` into `const QUANT_SCALE: f64 = (1u64 << 30) as f64;`
- Replace all 4 usage sites atomically: `cell_key`, `edge_key`, `sort_vertices_around_edge` (2 locations in the decode)
- Resolution: 9.3e-10 per unit (was 9.5e-7)

---

## 1.1 Octree 2:1 Balancing

**Files:** `src/octree.rs`

### Borrow Checker Strategy: Work-Queue

The recursive ownership (`children: Option<Box<[OctreeCell; 8]>>`) prevents simultaneous immutable traversal + mutation. Use a two-pass loop:

```
loop:
    violations = find_balance_violations(&tree)  // immutable borrow
    if violations.is_empty(): break
    for bbox in violations:
        subdivide_at(&mut tree, &node, bbox)     // mutable borrow
```

Converges in at most `max_depth` iterations (each pass creates violations at most one level deeper).

### Neighbor Finding: Probe-Based

No addresses needed. Each cell stores its bbox:

```
face_neighbor(root, cell, axis, direction):
    probe = cell.center()
    probe[axis] += direction * cell.size()[axis]
    return find_leaf_containing(root, probe)
```

O(max_depth) per lookup. Only need 6 face-neighbors per cell (not edge/corner — face 2:1 implies sufficient edge constraint for T-junction-free DC).

### Forced Subdivision (Not Build Subdivision)

When a pruned cell is subdivided for balance, use a dedicated function that:
- Evaluates corners only (no interval pruning)
- Creates 8 leaf children regardless of sign
- Does NOT recurse deeper

DC correctly ignores these cells (no sign change = no vertex). The outer balance loop handles cascading.

### Public API

```rust
impl Octree {
    pub fn balance(&mut self, node: &SdfNode) { ... }
}
```

Called between `Octree::build()` and DC extraction. The `extract_mesh_adaptive` function calls it automatically.

---

## 1.2 Watertight Guarantee

**Files:** `src/dual_contouring.rs`

### Bbox Padding in extract_mesh_adaptive

Move padding into `extract_mesh_adaptive` itself so all callers get it:

```rust
let cell_size = bbox.size().x / (1u64 << settings.max_depth) as f64;
let min_pad = cell_size * 2.0;
let pct_pad = bbox.size() * 0.05;
let pad = Vector3::new(
    pct_pad.x.max(min_pad),
    pct_pad.y.max(min_pad),
    pct_pad.z.max(min_pad),
);
let padded = BBox3::new(bbox.min - pad, bbox.max + pad);
```

This guarantees at least 2 cells of padding at finest resolution, regardless of bbox size.

### Remove padding from Shape::mesh()

Since `extract_mesh_adaptive` now handles padding, remove the 5% pad in `builder.rs:772-774` to avoid double-padding.

### emit_fan 2-cell edge handling

After padding, 2-cell edges at domain boundary should never occur for closed SDFs. Keep the `n < 3` skip as-is (it's correct for the impossible case). No degenerate triangle emission needed.

---

## 1.3 Consistent Winding

**Files:** `src/dual_contouring.rs`

### Gradient-Based Post-Processing

Replace the fragile edge-key-search winding (`dual_contouring.rs:116-135`) with a post-processing pass after all faces are generated, before `split_sharp_edges`:

```rust
fn fix_winding(vertices: &[Vector3<f64>], indices: &mut [u32], node: &SdfNode) {
    for tri in indices.chunks_exact_mut(3) {
        let v0 = vertices[tri[0] as usize];
        let v1 = vertices[tri[1] as usize];
        let v2 = vertices[tri[2] as usize];
        let centroid = (v0 + v1 + v2) / 3.0;
        let face_normal = (v1 - v0).cross(&(v2 - v0));
        let sdf_gradient = node.gradient(centroid);
        if face_normal.dot(&sdf_gradient) < 0.0 {
            tri.swap(1, 2);
        }
    }
}
```

The existing `sign_positive_first` logic in the face generation loop can be simplified or removed — the post-process handles all winding. However, keeping a best-effort attempt in emit_fan reduces the number of flips needed.

---

## Testing Strategy

**File:** `tests/dc_test.rs`

### Helper: strict manifold assertion

```rust
fn assert_manifold(mesh: &TriangleMesh) {
    // Every undirected edge shared by exactly 2 triangles
    // V - E + F = 2 for genus-0
}
```

### Test shapes (6 total)

1. **Sphere** (depth 6) — smooth, genus-0
2. **Box** (depth 6) — sharp edges
3. **Cylinder** (depth 6) — mixed flat/curved
4. **Union of two overlapping spheres** (depth 6) — CSG, tests sign changes at intersection
5. **SmoothUnion of two overlapping spheres** (depth 6) — blend zone stress test, rapid gradient rotation
6. **Thin box** (half_extents [1.0, 1.0, 0.01], depth 8) — features thinner than coarse cell size

All must pass `assert_manifold` with zero non-manifold edges.
