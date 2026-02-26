# Crusst V2: Mission-Critical Geometry Kernel Design

**Date:** 2026-02-26
**Status:** Approved
**Target:** Full-stack SDF kernel for parametric CAD, simulation, and additive manufacturing

## 1. Architecture: Hybrid Trait + Expression Graph

The `Sdf` trait remains the user-facing API. Internally, all built-in shapes and operations
construct an `SdfNode` expression DAG (Directed Acyclic Graph) enabling:

- **CPU computation** via trait `evaluate()` (backwards compatible)
- **GPU compilation** via `compile_wgsl()` to wgpu compute shaders
- **Interval arithmetic** via `interval_evaluate(bbox)` for spatial acceleration
- **Automatic differentiation** via `gradient(point)` for exact normals
- **Tape optimization** via linear bytecode compilation with region-based pruning

### SdfNode Enum

Primitives: Sphere, Box3, RoundedBox, Cylinder, CappedCone, Torus, Capsule, Ellipsoid, HalfSpace
CSG: Union, Intersection, Difference, SmoothUnion, SmoothIntersection, SmoothDifference
Transforms: Translate, Rotate, UniformScale, Mirror
Domain: Repeat, Onion (shell), Round (offset), Elongate
Advanced: Extrude2D, Revolution, Twist, Bend

## 2. Primitive Expansion

### Phase 1 (Foundation)
- Torus: exact SDF from Quilez
- Rounded Box: exact SDF with edge radius
- Capsule: sphere-swept line segment
- Ellipsoid: approximate via axis-scaled sphere
- Rounded Cylinder: cylinder with rounded edges

### Phase 2 (Engineering)
- 2D profile extrusion along axis
- Revolution surface from 2D profile
- Bezier curve SDF (cubic, distance solver)
- Infinite cylinder/plane

### Phase 3 (Advanced)
- TPMS surfaces (Gyroid, Schwarz-P, Diamond)
- Onion (shell) operation
- Annular/sector shapes

## 3. Transform System

Composable transforms via inverse-transform of query points.
Builder methods: translate, rotate, scale (uniform only), mirror, union, intersect, subtract.
Non-uniform scaling is forbidden (breaks distance property).

## 4. Adaptive Dual Contouring Mesher

Replace `isosurface` crate with custom pipeline:

### Stage 1: Adaptive Octree
Root cell covers bounding box. Interval arithmetic prunes empty cells.
Recursive subdivision to target resolution near surface only.

### Stage 2: Dual Contouring + QEF
Edge crossing detection via sign changes. Gradient computation at crossings.
QEF (Quadratic Error Function) vertex placement. Sharp edge and corner recovery.

### Stage 3: Mesh Quality Guarantees
Manifold dual contouring (multiple vertices per cell).
Watertight guarantee (every edge has exactly 2 faces).
Self-intersection detection and resolution.

## 5. GPU Pipeline (wgpu, all Rust)

Phase 1: WGSL code generation from SdfNode DAG
Phase 2: Batch compute dispatch (grid points to distances + gradients)
Phase 3: Sphere-traced real-time preview renderer
Phase 4: GPU-accelerated meshing (parallel MC/DC)

## 6. Performance Acceleration (CPU)

1. Interval Arithmetic: skip cells where interval excludes zero
2. Lipschitz Pruning: reduce CSG nodes per-region (up to 629x, Barbier 2025)
3. BVH: O(log n) culling for large unions
4. Tape Optimization: compile DAG to linear tape, prune per region (libfive approach)

## 7. Export Formats

1. Binary STL (exists)
2. OBJ (ASCII mesh)
3. STEP AP203 (tessellated shell for CAD interop)
4. glTF 2.0 (web visualization)
5. 3MF (modern AM format)
6. PLY (research/scanning)

## 8. Testing Framework

Layer 1: Distance accuracy (random-sample verification within tolerance)
Layer 2: Boolean correctness (containment, co-planar stress, tangent geometries)
Layer 3: Mesh quality (manifold, watertight, self-intersection, triangle quality)
Layer 4: Regression (golden-file comparison for standard models)
Layer 5: Benchmarks (criterion: throughput, meshing time, GPU vs CPU)
Layer 6: Showcase (complex engineering models, STL + STEP export)

## Decision Log

| Decision | Rationale |
|----------|-----------|
| Hybrid trait + DAG | Preserves ergonomic API while enabling GPU/acceleration |
| Own mesher over isosurface | Alpha dependency, no sharp features, no adaptive refinement |
| GPU from start (wgpu) | All-Rust, dramatic performance, modern kernel design |
| STEP via tessellation | Shows approximation quality; full NURBS export deferred |
| Uniform scale only | Non-uniform breaks distance field; document clearly |
