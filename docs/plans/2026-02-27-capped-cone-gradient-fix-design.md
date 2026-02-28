# CappedCone Analytical Gradient Fix — Design

## Problem

The Quilez CappedCone SDF formula produces gradient discontinuities at the
cap-wall junction. The dual contouring mesher computes gradients via central
finite differences (eps = 1e-6). When a stencil straddles the discontinuity,
the estimated gradient is a wrong average of two constraint surfaces. This
feeds incorrect normals to the QEF solver, which misplaces vertices — creating
jagged "chewy" edges instead of perfectly circular cap rims.

## Requirements

- Perfectly circular base cap edge on the meshed cone
- Watertight manifold mesh output
- Pure Rust — no new dependencies
- All existing 61 tests must pass

## Approach: Closest-Point SDF + Analytical Gradient Trait

### 1. Sdf Trait Extension

Add an optional `gradient()` method:

```rust
pub trait Sdf: Send + Sync {
    fn evaluate(&self, point: Vector3<f64>) -> f64;
    fn gradient(&self, point: Vector3<f64>) -> Option<Vector3<f64>> { None }
}
```

Default returns `None` — fully backward-compatible.

### 2. CappedCone Closest-Point Decomposition

Replace the Quilez formula with a closest-point formulation. For axis A→B,
radii ra/rb, project point P onto the axis to get (t, r):

| Voronoi Region | Closest Feature | Gradient Direction |
|----------------|----------------|--------------------|
| Conical wall   | Perpendicular projection onto slant | Outward normal to cone surface |
| Cap A face     | Projection onto cap A disk plane | -axis |
| Cap B face     | Projection onto cap B disk plane | +axis |
| Cap A rim      | Nearest point on base circle | Toward P from circle point |
| Cap B rim      | Nearest point on tip circle | Toward P from circle point |

The gradient is `(P - closest) / |P - closest|` for exterior, negated for interior.

### 3. DC Mesher Update

In `dual_contouring.rs`, check `sdf.gradient(point)` before calling
`central_diff_gradient_sdf()`. If `Some(n)` is returned, use it directly.

### 4. Files

| File | Change |
|------|--------|
| `src/shape.rs` | Add `gradient()` to Sdf trait |
| `src/primitives.rs` | New `sdf_capped_cone_closest_point()` + gradient |
| `src/dual_contouring.rs` | Prefer analytical gradient over central diffs |
| `tests/capped_cone_test.rs` | Circularity, watertightness, gradient accuracy |

### 5. Testing

- **Circularity:** Mesh a cone, extract base rim vertices, verify they lie on a circle within tolerance
- **Watertightness:** Every mesh edge shared by exactly 2 triangles
- **Gradient accuracy:** Analytical vs central-diff agreement at 10K random points
- **Regression:** All 61 existing tests pass
