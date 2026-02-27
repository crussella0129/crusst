# Fillet & Chamfer System Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add engineering-grade fillet and chamfer operations with exact profile curves, selective edge targeting via Feature IDs (FT), and full G1-G3 + exotic profile support.

**Architecture:** Two-layer system: (1) global CSG blend variants (RoundUnion, ChamferUnion, etc.) applied everywhere two surfaces meet, and (2) targeted Fillet/Chamfer nodes that use face identification and Feature Targeting (FT) to apply profiles to specific edges. All profiles compute exact signed distance to the profile curve in the (d1, d2) face-distance plane.

**Tech Stack:** Rust, nalgebra, rayon (parallel computation). No new dependencies.

**Design doc:** `docs/plans/2026-02-27-fillet-chamfer-design.md`

---

### Task 1: Blend profile math — core profiles (G2, EqualChamfer, G1)

Create `src/blend.rs` with the `BlendProfile` enum and exact blend functions for the three most important profiles. Create `tests/blend_test.rs` with ~11 tests covering tangent points, corner behavior, surface midpoints, and deep-inside fallback. Add `pub mod blend;` to `src/lib.rs`.

Key implementation details:

**BlendProfile enum:** 10 variants (G1, G2, G3, Chord, Cycloidal, Parabolic, Hyperbolic, EqualChamfer, TwoDistChamfer, AngleChamfer). Each has a `radius()` and `max_deviation()` method.

**Core blend function:** `blend_intersection(d1, d2, r, profile) -> f64` computes the SDF for an intersection (convex edge) case. Union and difference are derived by sign flips.

**G2 Circular (exact, closed-form):** Arc centered at `(-r, -r)` in `(d1, d2)` space with radius `r`. Formula:
```
u = max(d1 + r, 0); v = max(d2 + r, 0);
arc_sdf = sqrt(u*u + v*v) - r;
result = max(sharp, arc_sdf)
```

**EqualChamfer (exact, closed-form):** Line from `(-k, 0)` to `(0, -k)`. Formula:
```
chamfer_sdf = (d1 + d2 + k) / sqrt(2);
result = max(sharp, chamfer_sdf)
```

**G1 Hermite (Newton iteration):** Cubic Bezier from `(-r, 0)` to `(0, -r)` with control points giving tangent continuity. Nearest-point via Newton iteration (~5 iterations).

Also provide shorthand constructors: `g1(r)`, `g2(r)`, `equal_chamfer(k)`, etc.

---

### Task 2: Blend profile math — remaining 7 profiles

Add G3, Chord, Cycloidal, Parabolic, Hyperbolic, TwoDistChamfer, and AngleChamfer to `blend_intersection`. Add ~8 more tests.

**G3 Quintic:** Quintic Bezier with curvature-rate-matching control points. Newton iteration for nearest point.

**Chord:** Delegated to G2 with `r = chord_length / sqrt(2)`.

**TwoDistChamfer:** Line from `(-k1, 0)` to `(0, -k2)`. Closed-form: `(d1/k1 + d2/k2 + 1) * k1*k2 / sqrt(k1^2 + k2^2)`.

**AngleChamfer:** Converts angle to two distances, then delegates to TwoDistChamfer.

**Parabolic:** Parabolic arc from `(-r, 0)` to `(0, -r)`. Nearest point via Cardano's cubic formula (closed-form).

**Cycloidal:** Parametric cycloid `d1=r(t-sin t), d2=r(1-cos t)`. Newton iteration.

**Hyperbolic:** Hyperbolic arc `d1*d2 = a^2` scaled to fit. Newton iteration.

---

### Task 3: Feature ID types and face/edge enumeration

Create `src/feature.rs` with `FeatureTarget`, `FeatureKind`, `FaceInfo`, `EdgeInfo`, and the `ft()` builder helper. Add `face_info()`, `edge_info()`, `closest_face()`, and `face_distance()` methods to `SdfNode` in `src/dag.rs`. Create `tests/feature_test.rs` with ~9 tests.

**Box3 faces:** 6 faces (+X=0, -X=1, +Y=2, -Y=3, +Z=4, -Z=5).

**Box3 edges:** 12 edges, each connecting two adjacent faces.

**Box3 face_distance (exact half-space decomposition):**
```
face 0 (+X): q.x - half_extents.x
face 1 (-X): -q.x - half_extents.x
...
```

**Cylinder faces:** 3 (side=0, top=1, bottom=2). **Edges:** 2 (side/top=0, side/bottom=1).

**closest_face:** Uses normalized coordinate comparison (Box3) or distance comparison (Cylinder).

**Transforms pass through:** Translate, Rotate, Scale transform the point but preserve face numbering.

**ft() helper:** `ft(component, body).edges(&[2, 3, 4])` produces a `FeatureTarget`.

---

### Task 4: Global CSG blend functions (Round + Chamfer)

Add 6 functions to `src/csg.rs`: `round_union`, `round_intersection`, `round_difference`, `chamfer_union`, `chamfer_intersection`, `chamfer_difference`.

Add 6 matching `SdfNode` variants to `src/dag.rs` with `eval()`, `interval_eval()`, and `gradient()` arms.

Add 7 builder methods to `src/builder.rs`: `.round_union()`, `.round_intersect()`, `.round_subtract()`, `.chamfer_union()`, `.chamfer_intersect()`, `.chamfer_subtract()`, plus `.round(r)` (already in DAG, just needs builder exposure).

Add `compute_bbox` arms (same pattern as SmoothUnion — expand by radius).

Create `tests/round_csg_test.rs` with ~4 tests.

**Key formulas (csg.rs):**

Round union (Quilez inscribed-circle):
```
u = max(r - d1, 0); v = max(r - d2, 0);
result = max(r, min(d1, d2)) - sqrt(u*u + v*v)
```

Round intersection:
```
u = max(r + d1, 0); v = max(r + d2, 0);
result = min(-r, max(d1, d2)) + sqrt(u*u + v*v)
```

Chamfer union:
```
result = min(min(d1, d2), (d1 + d2 - k) / sqrt(2))
```

Chamfer intersection:
```
result = max(max(d1, d2), (d1 + d2 + k) / sqrt(2))
```

**Interval arithmetic:** Same pattern as SmoothUnion — sharp bounds expanded by max_deviation. Round expands by `r`, chamfer by `k * 0.3`.

**Gradient:** Central finite differences (consistent with existing smooth CSG).

---

### Task 5: Builder query methods (.faces(), .edges())

Add `.faces()` and `.edges()` methods to `Shape` in `builder.rs` that delegate to `self.node.face_info()` and `self.node.edge_info()`. Add 3 tests verifying box face/edge counts and transform pass-through.

---

### Task 6: Targeted Fillet/Chamfer SdfNode variants

Add `Fillet` and `Chamfer` variants to `SdfNode` with fields: `inner: Arc<SdfNode>`, `profile: BlendProfile`, `targets: Vec<FeatureTarget>`.

Add `eval()` logic:
1. Compute sharp distance via `inner.eval(point)`
2. Resolve target edges: convert FeatureTargets to (face_a, face_b) pairs using `edge_info()`
3. For each target edge: compute `d1 = face_distance(point, face_a)`, `d2 = face_distance(point, face_b)`
4. If in blend zone (`d1 > -r AND d2 > -r`): return `blend_intersection(d1, d2, r, profile)`
5. Otherwise: return sharp

Add `interval_eval()`: sharp bounds expanded by `profile.max_deviation()`.

Add `gradient()`: central finite differences.

Add `compute_bbox()`: inner bbox expanded by profile radius.

Add `.fillet(profile, target)` and `.chamfer(profile, target)` to Shape builder.

Create `tests/fillet_test.rs` with ~5 dimensional accuracy tests: corner removal, tangent points, arc midpoint, selective targeting, face-based targeting.

---

### Task 7: Fillet/chamfer mesh quality tests

Create `tests/fillet_mesh_test.rs` with ~4 integration tests verifying that filleted/chamfered shapes produce valid meshes through the dual contouring pipeline. Tests: filleted box meshes without panic, round union meshes, chamfered box meshes, rounded box has reasonable vertex count.

---

### Task 8: Showcase models for viewer

Add 6 new showcase models to `examples/viewer.rs`:
- `28_filleted_box`: Box3 with G2 fillet on all edges
- `29_chamfered_box`: Box3 with equal chamfer on all edges
- `30_round_union`: Two spheres with round union
- `31_round_subtract`: Box minus cylinder with round fillet
- `32_selective_fillet`: Box with fillets only on top face edges
- `33_rounded_box_simple`: Box3.round(0.5)

---

### Task 9: Final validation

Run `cargo fmt`, `cargo clippy -- -D warnings`, and `cargo test --release`. Fix any issues. Verify final test count (~150+, 0 failures, 0 clippy warnings, clean fmt).
