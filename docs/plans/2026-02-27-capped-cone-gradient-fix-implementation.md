# CappedCone Gradient Fix — Implementation Plan

> **Goal:** Fix jagged CappedCone mesh edges by replacing the Quilez SDF with a
> closest-point formulation + analytical gradient, and extending the Sdf trait.

**Architecture:** Extend `Sdf` trait with optional `gradient()`. Rewrite
CappedCone SDF as closest-point decomposition. Update DC mesher to prefer
analytical gradients.

**Tech Stack:** Rust, nalgebra

---

### Task 1: Extend Sdf Trait with gradient()

**Files:** `src/shape.rs`

Add `fn gradient(&self, point: Vector3<f64>) -> Option<Vector3<f64>> { None }`
to the `Sdf` trait. Default implementation returns `None`.

### Task 2: Rewrite CappedCone SDF

**Files:** `src/primitives.rs`

Replace `sdf_capped_cone()` with a closest-point decomposition that returns
both distance and gradient. The function signature becomes:

```rust
pub fn sdf_capped_cone_with_normal(
    point: Vector3<f64>, a: Vector3<f64>, b: Vector3<f64>, ra: f64, rb: f64,
) -> (f64, Vector3<f64>)
```

Keep the old function as a thin wrapper (returns just the distance) for API
compatibility.

### Task 3: Implement CappedCone::gradient()

**Files:** `src/shape.rs`

Override `gradient()` on `CappedCone` to call the new closest-point function
and return the analytical gradient.

### Task 4: Update DC Mesher

**Files:** `src/dual_contouring.rs`

In `compute_cell_vertex_sdf()` and final normal computation, check
`sdf.gradient(point)` first. If `Some(n)`, use it. Otherwise fall back to
`central_diff_gradient_sdf()`.

### Task 5: Tests

**Files:** `tests/capped_cone_test.rs`

- Gradient accuracy: compare analytical vs numerical at random points
- Circularity: mesh cone, check base rim vertices lie on circle
- Watertightness: every edge shared by exactly 2 triangles
- Regression: existing tests still pass

### Task 6: Commit and verify

Commit, run full test suite, verify all pass.
