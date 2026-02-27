# Fillet & Chamfer System Design

**Goal:** Engineering-grade fillet and chamfer operations with exact profile curves, selective edge targeting via a Feature ID (FT) addressing system, and full G1-G3 + exotic profile support.

**Approach:** Hybrid A+C — multi-channel SDF face tracking for general identification, plus half-space decomposition for primitives with planar faces (box, cylinder) for mathematically exact face distances.

**Key Constraint:** All profiles must be dimensionally exact (computed to f64 machine precision), not graphics approximations. The SDF blend at each point is the exact signed distance to the profile curve in the (d1, d2) plane.

---

## 1. Unified Blend Profile System

Every fillet and chamfer is a **profile curve** in the (d1, d2) plane, where d1 and d2 are signed distances to the two meeting faces. The SDF blend value is the exact signed distance to that profile curve.

### BlendProfile Enum

```rust
pub enum BlendProfile {
    // --- Fillets (smooth profiles) ---
    G1 { radius: f64 },          // Tangent-continuous (cubic Hermite)
    G2 { radius: f64 },          // Curvature-continuous (circular arc)
    G3 { radius: f64 },          // Curvature-rate-continuous (quintic Hermite)
    Chord { chord_length: f64 },  // Chord-length fillet (G2 with r = chord/sqrt(2))
    Cycloidal { radius: f64 },    // Cycloid profile (stress-optimized)
    Parabolic { radius: f64 },    // Parabolic arc (aerospace)
    Hyperbolic { radius: f64, asymptote: f64 }, // Hyperbolic arc

    // --- Chamfers (flat/angled profiles) ---
    EqualChamfer { distance: f64 },
    TwoDistChamfer { distance1: f64, distance2: f64 },
    AngleChamfer { distance: f64, angle_rad: f64 },
}
```

### Profile Math

Each profile implements `exact_blend(d1, d2) -> f64` for the intersection case. Union and difference derived by sign flips.

| Profile | Curve in (d1,d2) plane | Nearest-point method |
|---------|----------------------|---------------------|
| G2 Circular | Quarter-circle: `(r-d1)^2 + (r-d2)^2 = r^2` | Radial projection (closed-form, O(1)) |
| G1 Hermite | Cubic tangent to both faces at endpoints | Quintic root-finding (Newton, ~3-5 iters) |
| G3 Quintic | Quintic with curvature-rate boundary conditions | 9th degree polynomial (Newton, ~4-6 iters) |
| Chord | Circle with `r = chord/sqrt(2)` | Same as G2 (closed-form, O(1)) |
| Cycloidal | `d1=r(t-sin t), d2=r(1-cos t)` | Newton on parametric curve (~3-5 iters) |
| Parabolic | `d2 = r - d1^2/(4r)` | Cubic equation (Cardano, closed-form, O(1)) |
| Hyperbolic | `d1*d2 = a^2` scaled to radius r | Quartic or Newton (~3 iters) |
| EqualChamfer | Line: `d1 + d2 = k` | Perpendicular projection (closed-form, O(1)) |
| TwoDistChamfer | Line: `d1/k1 + d2/k2 = 1` | Perpendicular projection (closed-form, O(1)) |
| AngleChamfer | Line at angle theta | Perpendicular projection (closed-form, O(1)) |

### Key Invariant

All profiles satisfy: `blend(r, 0) == 0` and `blend(0, r) == 0` (curve touches both faces at the blend radius boundary).

### Interval Arithmetic Bounds

For octree pruning, each profile has a known maximum deviation from the sharp result. The interval is expanded conservatively: `sharp +/- max_deviation`. This does not affect profile accuracy.

| Profile | Max deviation from sharp |
|---------|-------------------------|
| G2 Circular | r |
| G1 Hermite | r |
| G3 Quintic | r |
| Cycloidal | ~0.36r |
| Parabolic | r/4 |
| Hyperbolic | a (asymptote) |
| All chamfers | k * 0.293 |

---

## 2. Feature ID System (FT)

Topological naming for faces and edges. Each primitive has a fixed, well-defined set of faces and edges with stable indices.

### Core Types

```rust
pub struct FeatureTarget {
    pub component: u32,
    pub body: u32,
    pub kind: FeatureKind,
    pub indices: Vec<u32>,
}

pub enum FeatureKind {
    Edge,  // 0
    Face,  // 1
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct FaceId {
    pub body: u32,
    pub face_index: u32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct EdgeId {
    pub body: u32,
    pub edge_index: u32,
}
```

### Face/Edge Enumeration Per Primitive

| Primitive | Faces | Edges |
|-----------|-------|-------|
| Box3 | 6: +X(0), -X(1), +Y(2), -Y(3), +Z(4), -Z(5) | 12: pairs of adjacent faces |
| Cylinder | 3: side(0), top(1), bottom(2) | 2: side-top(0), side-bottom(1) |
| CappedCone | 3: side(0), top(1), bottom(2) | 2: side-top(0), side-bottom(1) |
| Sphere | 1: surface(0) | 0 |
| Torus | 1: surface(0) | 0 |
| Capsule | 1: surface(0) | 0 |
| Ellipsoid | 1: surface(0) | 0 |
| RoundedBox | 6: same as Box3 | 12: same as Box3 |
| RoundedCylinder | 3: same as Cylinder | 2: same as Cylinder |
| HalfSpace | 1: plane(0) | 0 |

### Face Identification

Primitives implement `closest_face(point) -> u32` analytically:

- **Box3**: Largest normalized coordinate `|q_i / half_extent_i|` determines the face.
- **Cylinder**: Compare radial distance-to-side vs distance-to-caps.

### Transforms

Transforms preserve face/edge numbering:
- Translate, Rotate, Scale: face IDs pass through unchanged.
- Mirror: face IDs pass through, mirrored face pairs swap.

### Boolean Edges

CSG boolean operations create new edges where operand surfaces intersect. For V1, these are identified by body pair: "the edge where body A's surface meets body B's surface."

---

## 3. Multi-Channel SDF & SdfNode Integration

### New SdfNode Variants

```rust
// Targeted fillet/chamfer (requires feature tracking):
Fillet { inner: Arc<SdfNode>, profile: BlendProfile, targets: Vec<FeatureTarget> },
Chamfer { inner: Arc<SdfNode>, profile: BlendProfile, targets: Vec<FeatureTarget> },

// Global CSG blends (no feature tracking):
RoundUnion(Arc<SdfNode>, Arc<SdfNode>, f64),
RoundIntersection(Arc<SdfNode>, Arc<SdfNode>, f64),
RoundDifference(Arc<SdfNode>, Arc<SdfNode>, f64),
ChamferUnion(Arc<SdfNode>, Arc<SdfNode>, f64),
ChamferIntersection(Arc<SdfNode>, Arc<SdfNode>, f64),
ChamferDifference(Arc<SdfNode>, Arc<SdfNode>, f64),
```

### Face-Tracking Evaluation

```rust
impl SdfNode {
    pub fn evaluate(&self, point: Vector3<f64>) -> f64 { ... }           // unchanged
    fn evaluate_with_face(&self, point: Vector3<f64>) -> (f64, FaceId) { ... } // new, internal
}
```

- `evaluate()` (hot path) remains unchanged. Zero overhead for non-fillet shapes.
- `evaluate_with_face()` only activated by Fillet/Chamfer nodes on their inner subtree.
- Face IDs propagate through transforms (unchanged), CSG (dominating operand), Shell/Round (pass-through).

### Half-Space Decomposition (Approach C)

For Box3 and Cylinder caps, `face_distance(point, face_index) -> f64` returns the exact signed distance to a specific face. This feeds directly into the blend function for engineering-grade accuracy.

### Fillet Evaluation Logic

1. Compute sharp distance via `inner.evaluate(point)`.
2. For each targeted edge: compute face-pair distances (exact via decomposition or multi-channel).
3. If point is within blend radius of a targeted edge: return `profile.exact_blend(d1, d2)`.
4. Otherwise: return sharp distance.

### Touch Points Per New Variant

| Location | What |
|----------|------|
| `dag.rs::evaluate()` | Blend logic / CSG function call |
| `dag.rs::interval_evaluate()` | Conservative bounds (sharp +/- max_deviation) |
| `dag.rs::gradient()` | Central finite differences |
| `builder.rs::compute_bbox()` | Expand by blend radius |
| `step_export.rs::classify()` | All -> ExportTier::Tessellated (existing fallback handles it) |

---

## 4. Builder API

### Global Blends (No FT)

```rust
shape.round(radius)                   // global offset rounding (SdfNode::Round, already in DAG)
shape.round_union(other, radius)      // G2 circular fillet union
shape.round_subtract(other, radius)
shape.round_intersect(other, radius)
shape.chamfer_union(other, distance)  // 45 degree bevel union
shape.chamfer_subtract(other, distance)
shape.chamfer_intersect(other, distance)
```

### Targeted Fillets/Chamfers (With FT)

```rust
shape.fillet(profile, target)    // any profile, targeted edges/faces
shape.chamfer(profile, target)   // any profile, targeted edges/faces
```

### Feature Targeting Helper

```rust
ft(component, body).edges(&[2, 3, 4])   // specific edges
ft(component, body).edge(0)              // single edge
ft(component, body).face(4)              // all edges of a face
ft(component, body).all_edges()          // every edge
```

### Profile Shorthand Functions

```rust
g1(radius), g2(radius), g3(radius)
chord(length), cycloidal(radius), parabolic(radius), hyperbolic(radius, asymptote)
equal_chamfer(distance), two_dist_chamfer(d1, d2), angle_chamfer(distance, angle_rad)
```

### Example

```rust
// Fil(G3, 2.5) -> FT(3, 1, 0, <2,3,4>)
shape.fillet(g3(2.5), ft(3, 1).edges(&[2, 3, 4]))
```

### Query Methods

```rust
shape.faces() -> Vec<FaceInfo>   // enumerate faces with normals and labels
shape.edges() -> Vec<EdgeInfo>   // enumerate edges with face pairs and labels
```

---

## 5. Testing Strategy

### Test Levels

1. **Profile math** (~20 tests): Each blend function verified against known analytical values. Key invariant: `blend(r, 0) == 0` and `blend(0, r) == 0` for all profiles.

2. **Face/edge identification** (~12 tests): Primitives correctly report closest face, enumerate edges, return exact face distances.

3. **Integrated fillet geometry** (~15 tests): Full pipeline dimensional accuracy. Verify filleted surface passes through geometrically correct locations (tangent points, midpoints).

4. **Mesh quality** (~5 tests): Filleted shapes produce valid watertight meshes through dual contouring.

5. **Showcase models** (visual): 8 new viewer models demonstrating each profile type.

**Total: ~52 new tests, bringing the suite to ~158.**

---

## 6. Implementation Phases

| Phase | What | Dependencies |
|-------|------|-------------|
| 1 | `blend.rs` - BlendProfile enum + all 10 exact blend functions | None |
| 2 | `feature.rs` - FeatureTarget, FaceId, EdgeId, ft(), face/edge enumeration | None |
| 3 | Global CSG blends - 6 csg.rs functions + 6 SdfNode variants + builder | Phase 1 |
| 4 | `.round()` exposure + `.faces()`/`.edges()` query methods | Phase 2 |
| 5 | `evaluate_with_face()` - multi-channel SDF, face tracking through DAG | Phase 2 |
| 6 | Targeted Fillet/Chamfer SdfNode variants with half-space decomposition | Phase 1 + 5 |
| 7 | Builder API `.fillet()`/`.chamfer()` with FT targeting | Phase 6 |
| 8 | Showcase models + viewer integration | Phase 7 |

Phases 1 and 2 are independent and can run in parallel. Each phase has its own tests and commit.
