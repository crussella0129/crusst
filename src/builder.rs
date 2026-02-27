//! Declarative Shape builder API.
//!
//! `Shape` is the user-facing entry point for composing SDF geometry.
//! It wraps an `Arc<SdfNode>` so it is cheaply cloneable and immutable.
//!
//! # Example
//!
//! ```rust,no_run
//! use crusst::builder::Shape;
//!
//! let bracket = Shape::box3(20.0, 3.0, 15.0)
//!     .union(Shape::cylinder(5.0, 20.0).translate(0.0, 10.0, 0.0))
//!     .subtract(Shape::cylinder(3.0, 25.0))
//!     .translate(0.0, 5.0, 0.0);
//! bracket.export_obj("bracket.obj").unwrap();
//! ```

use std::sync::Arc;
use nalgebra::{Rotation3, Vector3};
use rayon::prelude::*;
use crate::dag::SdfNode;
use crate::types::{TriangleMesh, MeshSettings, BBox3};
use crate::voxel::VoxelGrid;
use crate::dual_contouring::extract_mesh_adaptive;

/// A composable, immutable shape backed by an SDF expression DAG.
///
/// `Shape` is cheaply cloneable (`Arc` under the hood) and can be combined
/// with CSG operations and transforms to build complex geometry.
#[derive(Clone)]
pub struct Shape {
    node: Arc<SdfNode>,
}

// ---------------------------------------------------------------------------
// Primitives — all centered at origin unless otherwise noted
// ---------------------------------------------------------------------------

impl Shape {
    /// Sphere centered at the origin.
    pub fn sphere(radius: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::Sphere {
                center: Vector3::zeros(),
                radius,
            }),
        }
    }

    /// Axis-aligned box centered at the origin.
    /// `hx`, `hy`, `hz` are the half-extents along each axis.
    pub fn box3(hx: f64, hy: f64, hz: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::Box3 {
                center: Vector3::zeros(),
                half_extents: Vector3::new(hx, hy, hz),
            }),
        }
    }

    /// Capped cylinder aligned along the Y axis, base at the origin.
    pub fn cylinder(radius: f64, height: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::Cylinder {
                base: Vector3::zeros(),
                axis: Vector3::new(0.0, 1.0, 0.0),
                radius,
                height,
            }),
        }
    }

    /// Torus lying in the XZ plane, centered at the origin.
    pub fn torus(major: f64, minor: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::Torus {
                center: Vector3::zeros(),
                major_radius: major,
                minor_radius: minor,
            }),
        }
    }

    /// Capsule (sphere-swept segment) from `a` to `b` with given `radius`.
    pub fn capsule(a: Vector3<f64>, b: Vector3<f64>, radius: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::Capsule { a, b, radius }),
        }
    }

    /// Capped cone (truncated cone) from `a` (radius `ra`) to `b` (radius `rb`).
    pub fn capped_cone(a: Vector3<f64>, b: Vector3<f64>, ra: f64, rb: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::CappedCone { a, b, ra, rb }),
        }
    }

    /// Box with rounded edges. Total size is `half_extents + radius`.
    pub fn rounded_box(hx: f64, hy: f64, hz: f64, radius: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::RoundedBox {
                center: Vector3::zeros(),
                half_extents: Vector3::new(hx, hy, hz),
                radius,
            }),
        }
    }

    /// Ellipsoid with semi-axis lengths `rx`, `ry`, `rz`.
    pub fn ellipsoid(rx: f64, ry: f64, rz: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::Ellipsoid {
                center: Vector3::zeros(),
                radii: Vector3::new(rx, ry, rz),
            }),
        }
    }

    /// Rounded cylinder aligned along the Y axis, centered at origin.
    pub fn rounded_cylinder(radius: f64, round_radius: f64, half_height: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::RoundedCylinder {
                center: Vector3::zeros(),
                radius,
                round_radius,
                half_height,
            }),
        }
    }

    /// Half-space: the solid region where `normal . p + d <= 0`.
    pub fn half_space(normal: Vector3<f64>, d: f64) -> Self {
        let n = normal.normalize();
        Self {
            node: Arc::new(SdfNode::HalfSpace { normal: n, d }),
        }
    }
}

// ---------------------------------------------------------------------------
// CSG operations — consume self, return new Shape
// ---------------------------------------------------------------------------

impl Shape {
    /// Boolean union: the volume of either shape.
    pub fn union(self, other: Shape) -> Self {
        Self {
            node: Arc::new(SdfNode::Union(self.node, other.node)),
        }
    }

    /// Boolean subtraction: self minus other.
    pub fn subtract(self, other: Shape) -> Self {
        Self {
            node: Arc::new(SdfNode::Difference(self.node, other.node)),
        }
    }

    /// Boolean intersection: the volume shared by both shapes.
    pub fn intersect(self, other: Shape) -> Self {
        Self {
            node: Arc::new(SdfNode::Intersection(self.node, other.node)),
        }
    }

    /// Smooth union (fillet blend) with blending radius `k`.
    pub fn smooth_union(self, other: Shape, k: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::SmoothUnion(self.node, other.node, k)),
        }
    }

    /// Smooth intersection with blending radius `k`.
    pub fn smooth_intersect(self, other: Shape, k: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::SmoothIntersection(self.node, other.node, k)),
        }
    }

    /// Smooth subtraction with blending radius `k`.
    pub fn smooth_subtract(self, other: Shape, k: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::SmoothDifference(self.node, other.node, k)),
        }
    }
}

// ---------------------------------------------------------------------------
// Transforms
// ---------------------------------------------------------------------------

impl Shape {
    /// Translate by `(x, y, z)`.
    pub fn translate(self, x: f64, y: f64, z: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::Translate(self.node, Vector3::new(x, y, z))),
        }
    }

    /// Rotate around the X axis by `angle` radians.
    pub fn rotate_x(self, angle: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::Rotate(
                self.node,
                Rotation3::from_axis_angle(&Vector3::x_axis(), angle),
            )),
        }
    }

    /// Rotate around the Y axis by `angle` radians.
    pub fn rotate_y(self, angle: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::Rotate(
                self.node,
                Rotation3::from_axis_angle(&Vector3::y_axis(), angle),
            )),
        }
    }

    /// Rotate around the Z axis by `angle` radians.
    pub fn rotate_z(self, angle: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::Rotate(
                self.node,
                Rotation3::from_axis_angle(&Vector3::z_axis(), angle),
            )),
        }
    }

    /// Uniform scale by `factor`.
    pub fn scale(self, factor: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::Scale(self.node, factor)),
        }
    }

    /// Mirror across the YZ plane (reflect X).
    pub fn mirror_x(self) -> Self {
        Self {
            node: Arc::new(SdfNode::Mirror(self.node, Vector3::new(1.0, 0.0, 0.0))),
        }
    }

    /// Mirror across the XZ plane (reflect Y).
    pub fn mirror_y(self) -> Self {
        Self {
            node: Arc::new(SdfNode::Mirror(self.node, Vector3::new(0.0, 1.0, 0.0))),
        }
    }

    /// Mirror across the XY plane (reflect Z).
    pub fn mirror_z(self) -> Self {
        Self {
            node: Arc::new(SdfNode::Mirror(self.node, Vector3::new(0.0, 0.0, 1.0))),
        }
    }

    /// Shell (hollow) with wall `thickness`.
    pub fn shell(self, thickness: f64) -> Self {
        Self {
            node: Arc::new(SdfNode::Shell(self.node, thickness)),
        }
    }
}

// ---------------------------------------------------------------------------
// Query
// ---------------------------------------------------------------------------

impl Shape {
    /// Evaluate the signed distance at `point`.
    /// Negative inside, zero on surface, positive outside.
    pub fn distance(&self, point: Vector3<f64>) -> f64 {
        self.node.evaluate(point)
    }

    /// Returns `true` if `point` is inside the shape (SDF < 0).
    pub fn contains(&self, point: Vector3<f64>) -> bool {
        self.node.evaluate(point) < 0.0
    }
}

// ---------------------------------------------------------------------------
// Bounding box computation — walk the DAG
// ---------------------------------------------------------------------------

impl Shape {
    /// Compute a conservative axis-aligned bounding box by walking the DAG.
    pub fn bounding_box(&self) -> BBox3 {
        compute_bbox(&self.node)
    }
}

/// Recursively compute a conservative AABB for an SdfNode.
fn compute_bbox(node: &SdfNode) -> BBox3 {
    match node {
        // -- Primitives --------------------------------------------------------

        SdfNode::Sphere { center, radius } => {
            let r = Vector3::new(*radius, *radius, *radius);
            BBox3::new(center - r, center + r)
        }

        SdfNode::Box3 { center, half_extents } => {
            BBox3::new(center - half_extents, center + half_extents)
        }

        SdfNode::Cylinder { base, axis, radius, height } => {
            // The cylinder goes from base to base + axis * height.
            // Compute a conservative AABB from the two end-cap discs.
            let top = base + axis * *height;

            // For each axis, the disc extends by radius * sin(angle between
            // cylinder axis and that coordinate axis). Conservative: use
            // radius for each perpendicular extent.
            let perp_extent = Vector3::new(
                radius * (1.0 - axis.x * axis.x).max(0.0).sqrt(),
                radius * (1.0 - axis.y * axis.y).max(0.0).sqrt(),
                radius * (1.0 - axis.z * axis.z).max(0.0).sqrt(),
            );

            let min1 = base - perp_extent;
            let max1 = base + perp_extent;
            let min2 = top - perp_extent;
            let max2 = top + perp_extent;

            BBox3::new(
                Vector3::new(min1.x.min(min2.x), min1.y.min(min2.y), min1.z.min(min2.z)),
                Vector3::new(max1.x.max(max2.x), max1.y.max(max2.y), max1.z.max(max2.z)),
            )
        }

        SdfNode::CappedCone { a, b, ra, rb } => {
            // Conservative: AABB that encloses both end-cap spheres.
            let r_max = ra.max(*rb);
            let ext = Vector3::new(r_max, r_max, r_max);
            let min_a = a - ext;
            let max_a = a + ext;
            let min_b = b - ext;
            let max_b = b + ext;
            BBox3::new(
                Vector3::new(min_a.x.min(min_b.x), min_a.y.min(min_b.y), min_a.z.min(min_b.z)),
                Vector3::new(max_a.x.max(max_b.x), max_a.y.max(max_b.y), max_a.z.max(max_b.z)),
            )
        }

        SdfNode::Torus { center, major_radius, minor_radius } => {
            // Torus lies in XZ plane. Extent in X and Z is major + minor,
            // extent in Y is just minor.
            let r_xz = major_radius + minor_radius;
            BBox3::new(
                Vector3::new(center.x - r_xz, center.y - minor_radius, center.z - r_xz),
                Vector3::new(center.x + r_xz, center.y + minor_radius, center.z + r_xz),
            )
        }

        SdfNode::RoundedBox { center, half_extents, radius } => {
            let total = half_extents + Vector3::new(*radius, *radius, *radius);
            BBox3::new(center - total, center + total)
        }

        SdfNode::Capsule { a, b, radius } => {
            let ext = Vector3::new(*radius, *radius, *radius);
            let min_a = a - ext;
            let max_a = a + ext;
            let min_b = b - ext;
            let max_b = b + ext;
            BBox3::new(
                Vector3::new(min_a.x.min(min_b.x), min_a.y.min(min_b.y), min_a.z.min(min_b.z)),
                Vector3::new(max_a.x.max(max_b.x), max_a.y.max(max_b.y), max_a.z.max(max_b.z)),
            )
        }

        SdfNode::Ellipsoid { center, radii } => {
            BBox3::new(center - radii, center + radii)
        }

        SdfNode::RoundedCylinder { center, radius, round_radius, half_height } => {
            let r_xz = radius + round_radius;
            let h = half_height + round_radius;
            BBox3::new(
                Vector3::new(center.x - r_xz, center.y - h, center.z - r_xz),
                Vector3::new(center.x + r_xz, center.y + h, center.z + r_xz),
            )
        }

        SdfNode::HalfSpace { .. } => {
            // Half-space is unbounded; use a large fallback.
            let big = 1000.0;
            BBox3::new(
                Vector3::new(-big, -big, -big),
                Vector3::new(big, big, big),
            )
        }

        // -- CSG ---------------------------------------------------------------

        SdfNode::Union(a, b) => {
            let ba = compute_bbox(a);
            let bb = compute_bbox(b);
            BBox3::new(
                Vector3::new(ba.min.x.min(bb.min.x), ba.min.y.min(bb.min.y), ba.min.z.min(bb.min.z)),
                Vector3::new(ba.max.x.max(bb.max.x), ba.max.y.max(bb.max.y), ba.max.z.max(bb.max.z)),
            )
        }

        SdfNode::Intersection(a, b) => {
            // Intersection is bounded by the overlap of both boxes.
            let ba = compute_bbox(a);
            let bb = compute_bbox(b);
            BBox3::new(
                Vector3::new(ba.min.x.max(bb.min.x), ba.min.y.max(bb.min.y), ba.min.z.max(bb.min.z)),
                Vector3::new(ba.max.x.min(bb.max.x), ba.max.y.min(bb.max.y), ba.max.z.min(bb.max.z)),
            )
        }

        SdfNode::Difference(a, _b) => {
            // Difference is bounded by the first operand.
            compute_bbox(a)
        }

        SdfNode::SmoothUnion(a, b, k) => {
            let ba = compute_bbox(a);
            let bb = compute_bbox(b);
            // Expand by k to account for the blending region.
            let ext = Vector3::new(*k, *k, *k);
            BBox3::new(
                Vector3::new(ba.min.x.min(bb.min.x), ba.min.y.min(bb.min.y), ba.min.z.min(bb.min.z)) - ext,
                Vector3::new(ba.max.x.max(bb.max.x), ba.max.y.max(bb.max.y), ba.max.z.max(bb.max.z)) + ext,
            )
        }

        SdfNode::SmoothIntersection(a, b, k) => {
            let ba = compute_bbox(a);
            let bb = compute_bbox(b);
            let ext = Vector3::new(*k, *k, *k);
            BBox3::new(
                Vector3::new(ba.min.x.max(bb.min.x), ba.min.y.max(bb.min.y), ba.min.z.max(bb.min.z)) - ext,
                Vector3::new(ba.max.x.min(bb.max.x), ba.max.y.min(bb.max.y), ba.max.z.min(bb.max.z)) + ext,
            )
        }

        SdfNode::SmoothDifference(a, _b, k) => {
            let ba = compute_bbox(a);
            let ext = Vector3::new(*k, *k, *k);
            BBox3::new(ba.min - ext, ba.max + ext)
        }

        // -- Transforms --------------------------------------------------------

        SdfNode::Translate(inner, offset) => {
            let b = compute_bbox(inner);
            BBox3::new(b.min + offset, b.max + offset)
        }

        SdfNode::Rotate(inner, rotation) => {
            let b = compute_bbox(inner);
            // Rotate all 8 corners and compute the new AABB.
            let corners = b.corners();
            let mut new_min = Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
            let mut new_max = Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
            for c in &corners {
                let rc = rotation * c;
                new_min.x = new_min.x.min(rc.x);
                new_min.y = new_min.y.min(rc.y);
                new_min.z = new_min.z.min(rc.z);
                new_max.x = new_max.x.max(rc.x);
                new_max.y = new_max.y.max(rc.y);
                new_max.z = new_max.z.max(rc.z);
            }
            BBox3::new(new_min, new_max)
        }

        SdfNode::Scale(inner, factor) => {
            let b = compute_bbox(inner);
            BBox3::new(b.min * *factor, b.max * *factor)
        }

        SdfNode::Mirror(inner, normal) => {
            let b = compute_bbox(inner);
            // The mirror includes the original and the reflected geometry.
            let corners = b.corners();
            let mut new_min = b.min;
            let mut new_max = b.max;
            for c in &corners {
                let d = c.dot(normal);
                // Only reflect points on the negative side
                if d < 0.0 {
                    let reflected = c - normal * (2.0 * d);
                    new_min.x = new_min.x.min(reflected.x);
                    new_min.y = new_min.y.min(reflected.y);
                    new_min.z = new_min.z.min(reflected.z);
                    new_max.x = new_max.x.max(reflected.x);
                    new_max.y = new_max.y.max(reflected.y);
                    new_max.z = new_max.z.max(reflected.z);
                }
            }
            BBox3::new(new_min, new_max)
        }

        SdfNode::Shell(inner, thickness) => {
            let b = compute_bbox(inner);
            let ext = Vector3::new(*thickness, *thickness, *thickness);
            BBox3::new(b.min - ext, b.max + ext)
        }

        SdfNode::Round(inner, radius) => {
            let b = compute_bbox(inner);
            let ext = Vector3::new(*radius, *radius, *radius);
            BBox3::new(b.min - ext, b.max + ext)
        }

        // -- 2D -> 3D ----------------------------------------------------------

        SdfNode::Revolve(_) | SdfNode::Extrude(_, _) => {
            // Conservative fallback for 2D->3D operations.
            let big = 100.0;
            BBox3::new(
                Vector3::new(-big, -big, -big),
                Vector3::new(big, big, big),
            )
        }

        // -- Opaque ------------------------------------------------------------

        SdfNode::Custom(_) => {
            let big = 100.0;
            BBox3::new(
                Vector3::new(-big, -big, -big),
                Vector3::new(big, big, big),
            )
        }
    }
}

// ---------------------------------------------------------------------------
// Meshing
// ---------------------------------------------------------------------------

impl Shape {
    /// Extract a triangle mesh using adaptive dual contouring.
    pub fn mesh(&self, settings: MeshSettings) -> TriangleMesh {
        let bbox = self.bounding_box();
        // Expand bbox slightly to avoid clipping surface cells at the boundary.
        let pad = bbox.size() * 0.05;
        let padded = BBox3::new(bbox.min - pad, bbox.max + pad);
        extract_mesh_adaptive(&self.node, &padded, &settings)
    }
}

// ---------------------------------------------------------------------------
// Export
// ---------------------------------------------------------------------------

impl Shape {
    /// Mesh with default settings and write to a Wavefront OBJ file.
    pub fn export_obj(&self, path: &str) -> std::io::Result<()> {
        let mesh = self.mesh(MeshSettings::default());
        crate::obj_export::write_obj(&mesh, std::path::Path::new(path))
    }

    /// Mesh with default settings and write to a binary STL file.
    pub fn export_stl(&self, path: &str) -> std::io::Result<()> {
        let mesh = self.mesh(MeshSettings::default());
        crate::export::write_stl(&mesh, std::path::Path::new(path))
    }
}

// ---------------------------------------------------------------------------
// Voxelization
// ---------------------------------------------------------------------------

impl Shape {
    /// Sample the SDF onto a regular 3D grid with the given `voxel_size`.
    ///
    /// The grid covers the shape's bounding box plus one voxel of padding
    /// on each side. Evaluation is parallelized across voxels with rayon.
    pub fn voxelize(&self, voxel_size: f64) -> VoxelGrid {
        let bbox = self.bounding_box();
        let pad = Vector3::new(voxel_size, voxel_size, voxel_size);
        let origin = bbox.min - pad;
        let extent = bbox.max + pad - origin;

        let nx = (extent.x / voxel_size).ceil() as usize;
        let ny = (extent.y / voxel_size).ceil() as usize;
        let nz = (extent.z / voxel_size).ceil() as usize;
        let total = nx * ny * nz;

        let node = &self.node;
        let data: Vec<f32> = (0..total)
            .into_par_iter()
            .map(|idx| {
                let ix = idx / (ny * nz);
                let iy = (idx / nz) % ny;
                let iz = idx % nz;
                let world = origin + Vector3::new(
                    (ix as f64 + 0.5) * voxel_size,
                    (iy as f64 + 0.5) * voxel_size,
                    (iz as f64 + 0.5) * voxel_size,
                );
                node.evaluate(world) as f32
            })
            .collect();

        VoxelGrid {
            resolution: [nx, ny, nz],
            voxel_size,
            origin,
            data,
        }
    }
}

// ---------------------------------------------------------------------------
// Access
// ---------------------------------------------------------------------------

impl Shape {
    /// Expose the inner DAG node (for STEP export and other introspection).
    pub fn node(&self) -> &SdfNode {
        &self.node
    }
}
