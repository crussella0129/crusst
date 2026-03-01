//! Fluent builder API for constructing B-Rep shapes.
//!
//! `Shape` wraps a `(TopoStore, SolidId)` pair and provides chainable methods
//! for creating primitives, applying transforms, and generating meshes.
//!
//! # Examples
//! ```
//! use crusst::builder::Shape;
//! use crusst::types::TessSettings;
//!
//! let mesh = Shape::box3(5.0, 3.0, 8.0)
//!     .translate(10.0, 0.0, 0.0)
//!     .mesh(&TessSettings::default());
//! ```

use crate::math::{Point3, Vector3};
use crate::primitive;
use crate::profile::Profile;
use crate::tessellate;
use crate::topo::*;
use crate::types::{TessSettings, TriangleMesh};

/// A B-Rep shape with fluent API.
pub struct Shape {
    pub store: TopoStore,
    pub solid: SolidId,
}

impl Shape {
    // --- Primitive constructors ---

    /// Axis-aligned box centered at origin with half-extents (hx, hy, hz).
    pub fn box3(hx: f64, hy: f64, hz: f64) -> Self {
        let mut store = TopoStore::new();
        let solid = primitive::make_box(&mut store, hx, hy, hz);
        Shape { store, solid }
    }

    /// Sphere centered at origin.
    pub fn sphere(radius: f64) -> Self {
        let mut store = TopoStore::new();
        let solid = primitive::make_sphere(&mut store, radius);
        Shape { store, solid }
    }

    /// Cylinder on the Z axis from z=0 to z=height.
    pub fn cylinder(radius: f64, height: f64) -> Self {
        let mut store = TopoStore::new();
        let solid = primitive::make_cylinder(&mut store, radius, height);
        Shape { store, solid }
    }

    /// Cone (frustum) on the Z axis. `r1` bottom radius, `r2` top radius.
    pub fn cone(r1: f64, r2: f64, height: f64) -> Self {
        let mut store = TopoStore::new();
        let solid = primitive::make_cone(&mut store, r1, r2, height);
        Shape { store, solid }
    }

    /// Torus centered at origin on the Z axis.
    pub fn torus(major_r: f64, minor_r: f64) -> Self {
        let mut store = TopoStore::new();
        let solid = primitive::make_torus(&mut store, major_r, minor_r);
        Shape { store, solid }
    }

    /// Wedge (triangular prism) on the Z axis.
    pub fn wedge(dx: f64, dy: f64, dz: f64) -> Self {
        let mut store = TopoStore::new();
        let solid = primitive::make_wedge(&mut store, dx, dy, dz);
        Shape { store, solid }
    }

    /// Capsule (hemisphere-capped cylinder) on the Z axis.
    pub fn capsule(radius: f64, height: f64) -> Self {
        let mut store = TopoStore::new();
        let solid = primitive::make_capsule(&mut store, radius, height);
        Shape { store, solid }
    }

    // --- Transform operations ---

    /// Translate the shape by (dx, dy, dz).
    pub fn translate(mut self, dx: f64, dy: f64, dz: f64) -> Self {
        let offset = Vector3::new(dx, dy, dz);
        for v in &mut self.store.vertices {
            v.point += offset;
        }
        // Update surface origins
        for face in &mut self.store.faces {
            translate_surface(&mut face.surface, &offset);
        }
        // Update edge curve origins
        for edge in &mut self.store.edges {
            translate_curve(&mut edge.curve, &offset);
        }
        self
    }

    /// Rotate around the X axis by `angle` radians.
    pub fn rotate_x(mut self, angle: f64) -> Self {
        let rot = nalgebra::Rotation3::from_axis_angle(&nalgebra::Vector3::x_axis(), angle);
        apply_rotation(&mut self.store, &rot);
        self
    }

    /// Rotate around the Y axis by `angle` radians.
    pub fn rotate_y(mut self, angle: f64) -> Self {
        let rot = nalgebra::Rotation3::from_axis_angle(&nalgebra::Vector3::y_axis(), angle);
        apply_rotation(&mut self.store, &rot);
        self
    }

    /// Rotate around the Z axis by `angle` radians.
    pub fn rotate_z(mut self, angle: f64) -> Self {
        let rot = nalgebra::Rotation3::from_axis_angle(&nalgebra::Vector3::z_axis(), angle);
        apply_rotation(&mut self.store, &rot);
        self
    }

    /// Uniform scale.
    pub fn scale(mut self, factor: f64) -> Self {
        for v in &mut self.store.vertices {
            v.point.coords *= factor;
        }
        for face in &mut self.store.faces {
            scale_surface(&mut face.surface, factor);
        }
        for edge in &mut self.store.edges {
            scale_curve(&mut edge.curve, factor);
        }
        self
    }

    /// Mirror across the YZ plane (negate X).
    pub fn mirror_x(self) -> Self {
        self.mirror(Vector3::new(-1.0, 1.0, 1.0))
    }

    /// Mirror across the XZ plane (negate Y).
    pub fn mirror_y(self) -> Self {
        self.mirror(Vector3::new(1.0, -1.0, 1.0))
    }

    /// Mirror across the XY plane (negate Z).
    pub fn mirror_z(self) -> Self {
        self.mirror(Vector3::new(1.0, 1.0, -1.0))
    }

    fn mirror(mut self, scale: Vector3) -> Self {
        for v in &mut self.store.vertices {
            v.point.x *= scale.x;
            v.point.y *= scale.y;
            v.point.z *= scale.z;
        }
        // Flip face orientations since mirroring reverses winding
        for face in &mut self.store.faces {
            face.outward = !face.outward;
            mirror_surface(&mut face.surface, &scale);
        }
        for edge in &mut self.store.edges {
            mirror_curve(&mut edge.curve, &scale);
        }
        self
    }

    // --- Profile-based operations ---

    /// Extrude a 2D profile along the Z axis by `depth`.
    ///
    /// Creates a solid by extruding the profile from z=0 to z=depth,
    /// capping both ends with planar faces.
    pub fn extrude(profile: &Profile, depth: f64) -> Self {
        // Build a prismatic solid from the profile
        // For now, approximate the profile as a polygon and build topology
        let n_samples = profile.segments.len().max(4) * 4;
        let mut bottom_pts: Vec<Point3> = Vec::with_capacity(n_samples);
        let mut top_pts: Vec<Point3> = Vec::with_capacity(n_samples);

        for i in 0..n_samples {
            let t = i as f64 / n_samples as f64;
            let p2 = profile.evaluate(t);
            bottom_pts.push(Point3::new(p2.x, p2.y, 0.0));
            top_pts.push(Point3::new(p2.x, p2.y, depth));
        }

        build_extrusion(bottom_pts, top_pts, depth)
    }

    /// Revolve a 2D profile around the Z axis by `angle` radians.
    pub fn revolve(profile: &Profile, angle: f64) -> Self {
        let n_profile = profile.segments.len().max(4) * 4;
        let n_sweep = ((angle / std::f64::consts::FRAC_PI_4).ceil() as usize).max(4) * 2;

        let mut store = TopoStore::new();
        let mut all_verts: Vec<Vec<VertexId>> = Vec::new();

        // Create vertex grid: profile × sweep
        for j in 0..=n_sweep {
            let theta = angle * j as f64 / n_sweep as f64;
            let cos_t = theta.cos();
            let sin_t = theta.sin();
            let mut ring = Vec::with_capacity(n_profile);

            for i in 0..n_profile {
                let t = i as f64 / n_profile as f64;
                let p2 = profile.evaluate(t);
                // Revolve: (x, y, 0) in profile → (x·cos(θ), x·sin(θ), y) in 3D
                let pt = Point3::new(p2.x * cos_t, p2.x * sin_t, p2.y);
                ring.push(store.add_vertex(Vertex { point: pt }));
            }
            all_verts.push(ring);
        }

        // Build quad faces between adjacent rings
        let mut face_ids = Vec::new();
        for j in 0..n_sweep {
            for i in 0..n_profile {
                let i_next = (i + 1) % n_profile;

                let v00 = all_verts[j][i];
                let v10 = all_verts[j][i_next];
                let v01 = all_verts[j + 1][i];
                let v11 = all_verts[j + 1][i_next];

                let face_id = build_quad_face(&mut store, v00, v10, v11, v01);
                face_ids.push(face_id);
            }
        }

        // Cap the ends if it's not a full revolution
        let is_full = (angle - std::f64::consts::TAU).abs() < 1e-6;
        if !is_full {
            // Start cap
            let start_face = build_fan_face(&mut store, &all_verts[0], true);
            face_ids.push(start_face);
            // End cap
            let end_face = build_fan_face(&mut store, &all_verts[n_sweep], false);
            face_ids.push(end_face);
        }

        let shell = store.add_shell(Shell { faces: face_ids });
        let solid = store.add_solid(Solid {
            outer_shell: shell,
            inner_shells: vec![],
        });
        Shape { store, solid }
    }

    // --- Output ---

    /// Compute the bounding box diagonal of this shape from its topology vertices.
    pub fn bounding_diagonal(&self) -> f64 {
        let vert_ids = self.store.solid_vertices(self.solid);
        if vert_ids.is_empty() {
            return 1.0;
        }
        let mut min = [f64::MAX; 3];
        let mut max = [f64::MIN; 3];
        for &vid in &vert_ids {
            let p = self.store.vertex(vid).point;
            min[0] = min[0].min(p.x);
            min[1] = min[1].min(p.y);
            min[2] = min[2].min(p.z);
            max[0] = max[0].max(p.x);
            max[1] = max[1].max(p.y);
            max[2] = max[2].max(p.z);
        }
        let dx = max[0] - min[0];
        let dy = max[1] - min[1];
        let dz = max[2] - min[2];
        (dx * dx + dy * dy + dz * dz).sqrt().max(1e-10)
    }

    /// Generate a triangle mesh from this shape with explicit settings.
    pub fn mesh(&self, settings: &TessSettings) -> TriangleMesh {
        tessellate::tessellate_solid(&self.store, self.solid, settings)
    }

    /// Generate a triangle mesh with settings auto-tuned to this shape's size.
    pub fn auto_mesh(&self) -> TriangleMesh {
        let settings = TessSettings::from_bounding_diagonal(self.bounding_diagonal());
        self.mesh(&settings)
    }

    /// Validate the topology of this shape.
    pub fn validate(&self) -> ValidationResult {
        validate_solid(&self.store, self.solid)
    }

    // --- Export ---

    /// Export as binary STL.
    pub fn write_stl<W: std::io::Write>(&self, settings: &TessSettings, writer: &mut W) -> std::io::Result<()> {
        let mesh = self.mesh(settings);
        crate::export::write_stl(&mesh, writer)
    }

    /// Export as Wavefront OBJ.
    pub fn write_obj<W: std::io::Write>(&self, settings: &TessSettings, writer: &mut W) -> std::io::Result<()> {
        let mesh = self.mesh(settings);
        crate::export::write_obj(&mesh, writer)
    }

    /// Export as STEP AP203 (exact geometry, no tessellation).
    pub fn write_step<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        crate::export::write_step(&self.store, self.solid, writer)
    }

    /// Export as 3MF model XML.
    pub fn write_3mf<W: std::io::Write>(&self, settings: &TessSettings, writer: &mut W) -> std::io::Result<()> {
        let mesh = self.mesh(settings);
        crate::export::write_3mf(&mesh, writer)
    }
}

// --- Internal helpers ---

fn translate_surface(surface: &mut crate::surface::Surface, offset: &Vector3) {
    use crate::surface::Surface;
    match surface {
        Surface::Plane { origin, .. } => *origin += offset,
        Surface::Cylinder { origin, .. } => *origin += offset,
        Surface::Cone { apex, .. } => *apex += offset,
        Surface::Sphere { center, .. } => *center += offset,
        Surface::Torus { center, .. } => *center += offset,
        Surface::NurbsSurface(nurbs) => {
            for row in &mut nurbs.control_points {
                for pt in row {
                    *pt += offset;
                }
            }
        }
    }
}

fn translate_curve(curve: &mut crate::curve::Curve3, offset: &Vector3) {
    use crate::curve::Curve3;
    match curve {
        Curve3::Line { origin, .. } => *origin += offset,
        Curve3::Circle { center, .. } => *center += offset,
        Curve3::Ellipse { center, .. } => *center += offset,
        Curve3::NurbsCurve(nurbs) => {
            for pt in &mut nurbs.control_points {
                *pt += offset;
            }
        }
    }
}

fn apply_rotation(store: &mut TopoStore, rot: &nalgebra::Rotation3<f64>) {
    for v in &mut store.vertices {
        v.point = Point3::from(rot * v.point.coords);
    }
    for face in &mut store.faces {
        rotate_surface(&mut face.surface, rot);
    }
    for edge in &mut store.edges {
        rotate_curve(&mut edge.curve, rot);
    }
}

fn rotate_surface(surface: &mut crate::surface::Surface, rot: &nalgebra::Rotation3<f64>) {
    use crate::surface::Surface;
    match surface {
        Surface::Plane { origin, normal } => {
            *origin = Point3::from(rot * origin.coords);
            *normal = rot * *normal;
        }
        Surface::Cylinder { origin, axis, .. } => {
            *origin = Point3::from(rot * origin.coords);
            *axis = rot * *axis;
        }
        Surface::Cone { apex, axis, .. } => {
            *apex = Point3::from(rot * apex.coords);
            *axis = rot * *axis;
        }
        Surface::Sphere { center, .. } => {
            *center = Point3::from(rot * center.coords);
        }
        Surface::Torus { center, axis, .. } => {
            *center = Point3::from(rot * center.coords);
            *axis = rot * *axis;
        }
        Surface::NurbsSurface(nurbs) => {
            for row in &mut nurbs.control_points {
                for pt in row {
                    *pt = Point3::from(rot * pt.coords);
                }
            }
        }
    }
}

fn rotate_curve(curve: &mut crate::curve::Curve3, rot: &nalgebra::Rotation3<f64>) {
    use crate::curve::Curve3;
    match curve {
        Curve3::Line { origin, dir } => {
            *origin = Point3::from(rot * origin.coords);
            *dir = rot * *dir;
        }
        Curve3::Circle { center, axis, .. } => {
            *center = Point3::from(rot * center.coords);
            *axis = rot * *axis;
        }
        Curve3::Ellipse { center, major, minor } => {
            *center = Point3::from(rot * center.coords);
            *major = rot * *major;
            *minor = rot * *minor;
        }
        Curve3::NurbsCurve(nurbs) => {
            for pt in &mut nurbs.control_points {
                *pt = Point3::from(rot * pt.coords);
            }
        }
    }
}

fn scale_surface(surface: &mut crate::surface::Surface, factor: f64) {
    use crate::surface::Surface;
    match surface {
        Surface::Plane { origin, .. } => origin.coords *= factor,
        Surface::Cylinder { origin, radius, .. } => {
            origin.coords *= factor;
            *radius *= factor;
        }
        Surface::Cone { apex, .. } => apex.coords *= factor,
        Surface::Sphere { center, radius } => {
            center.coords *= factor;
            *radius *= factor;
        }
        Surface::Torus { center, major_r, minor_r, .. } => {
            center.coords *= factor;
            *major_r *= factor;
            *minor_r *= factor;
        }
        Surface::NurbsSurface(nurbs) => {
            for row in &mut nurbs.control_points {
                for pt in row {
                    pt.coords *= factor;
                }
            }
        }
    }
}

fn scale_curve(curve: &mut crate::curve::Curve3, factor: f64) {
    use crate::curve::Curve3;
    match curve {
        Curve3::Line { origin, dir } => {
            origin.coords *= factor;
            *dir *= factor;
        }
        Curve3::Circle { center, radius, .. } => {
            center.coords *= factor;
            *radius *= factor;
        }
        Curve3::Ellipse { center, major, minor } => {
            center.coords *= factor;
            *major *= factor;
            *minor *= factor;
        }
        Curve3::NurbsCurve(nurbs) => {
            for pt in &mut nurbs.control_points {
                pt.coords *= factor;
            }
        }
    }
}

fn mirror_surface(surface: &mut crate::surface::Surface, s: &Vector3) {
    use crate::surface::Surface;
    match surface {
        Surface::Plane { origin, normal } => {
            origin.x *= s.x; origin.y *= s.y; origin.z *= s.z;
            normal.x *= s.x; normal.y *= s.y; normal.z *= s.z;
        }
        Surface::Cylinder { origin, axis, .. } => {
            origin.x *= s.x; origin.y *= s.y; origin.z *= s.z;
            axis.x *= s.x; axis.y *= s.y; axis.z *= s.z;
        }
        Surface::Cone { apex, axis, .. } => {
            apex.x *= s.x; apex.y *= s.y; apex.z *= s.z;
            axis.x *= s.x; axis.y *= s.y; axis.z *= s.z;
        }
        Surface::Sphere { center, .. } => {
            center.x *= s.x; center.y *= s.y; center.z *= s.z;
        }
        Surface::Torus { center, axis, .. } => {
            center.x *= s.x; center.y *= s.y; center.z *= s.z;
            axis.x *= s.x; axis.y *= s.y; axis.z *= s.z;
        }
        Surface::NurbsSurface(nurbs) => {
            for row in &mut nurbs.control_points {
                for pt in row {
                    pt.x *= s.x; pt.y *= s.y; pt.z *= s.z;
                }
            }
        }
    }
}

fn mirror_curve(curve: &mut crate::curve::Curve3, s: &Vector3) {
    use crate::curve::Curve3;
    match curve {
        Curve3::Line { origin, dir } => {
            origin.x *= s.x; origin.y *= s.y; origin.z *= s.z;
            dir.x *= s.x; dir.y *= s.y; dir.z *= s.z;
        }
        Curve3::Circle { center, axis, .. } => {
            center.x *= s.x; center.y *= s.y; center.z *= s.z;
            axis.x *= s.x; axis.y *= s.y; axis.z *= s.z;
        }
        Curve3::Ellipse { center, major, minor } => {
            center.x *= s.x; center.y *= s.y; center.z *= s.z;
            major.x *= s.x; major.y *= s.y; major.z *= s.z;
            minor.x *= s.x; minor.y *= s.y; minor.z *= s.z;
        }
        Curve3::NurbsCurve(nurbs) => {
            for pt in &mut nurbs.control_points {
                pt.x *= s.x; pt.y *= s.y; pt.z *= s.z;
            }
        }
    }
}

/// Build an extruded solid from bottom and top polygon rings.
fn build_extrusion(bottom: Vec<Point3>, top: Vec<Point3>, depth: f64) -> Shape {
    let n = bottom.len();
    let mut store = TopoStore::new();

    let bottom_vids: Vec<VertexId> = bottom
        .iter()
        .map(|p| store.add_vertex(Vertex { point: *p }))
        .collect();
    let top_vids: Vec<VertexId> = top
        .iter()
        .map(|p| store.add_vertex(Vertex { point: *p }))
        .collect();

    let mut face_ids = Vec::new();

    // Side faces (quads)
    for i in 0..n {
        let i_next = (i + 1) % n;
        let face_id = build_quad_face(
            &mut store,
            bottom_vids[i],
            bottom_vids[i_next],
            top_vids[i_next],
            top_vids[i],
        );
        face_ids.push(face_id);
    }

    // Bottom cap (fan)
    let bottom_face = build_fan_face(&mut store, &bottom_vids, false);
    face_ids.push(bottom_face);

    // Top cap (fan)
    let top_face = build_fan_face(&mut store, &top_vids, true);
    face_ids.push(top_face);

    let shell = store.add_shell(Shell { faces: face_ids });
    let solid = store.add_solid(Solid {
        outer_shell: shell,
        inner_shells: vec![],
    });
    Shape { store, solid }
}

/// Build a quad face from 4 vertices (v0→v1→v2→v3 in CCW order).
fn build_quad_face(
    store: &mut TopoStore,
    v0: VertexId,
    v1: VertexId,
    v2: VertexId,
    v3: VertexId,
) -> FaceId {
    use crate::curve::{Curve2, Curve3};
    use crate::math::{Point2, Vector2};

    let p0 = store.vertex(v0).point;
    let p1 = store.vertex(v1).point;
    let p2 = store.vertex(v2).point;
    let p3 = store.vertex(v3).point;

    // Compute face normal from vertices
    let e1 = p1 - p0;
    let e2 = p3 - p0;
    let normal = e1.cross(&e2);
    let len = normal.norm();
    let normal = if len > 1e-15 { normal / len } else { Vector3::new(0.0, 0.0, 1.0) };

    let surface = crate::surface::Surface::Plane {
        origin: p0,
        normal,
    };

    // 4 edges
    let edges = [
        store.add_edge(Edge {
            curve: Curve3::Line { origin: p0, dir: p1 - p0 },
            t_start: 0.0, t_end: 1.0, start: v0, end: v1,
        }),
        store.add_edge(Edge {
            curve: Curve3::Line { origin: p1, dir: p2 - p1 },
            t_start: 0.0, t_end: 1.0, start: v1, end: v2,
        }),
        store.add_edge(Edge {
            curve: Curve3::Line { origin: p2, dir: p3 - p2 },
            t_start: 0.0, t_end: 1.0, start: v2, end: v3,
        }),
        store.add_edge(Edge {
            curve: Curve3::Line { origin: p3, dir: p0 - p3 },
            t_start: 0.0, t_end: 1.0, start: v3, end: v0,
        }),
    ];

    let face_id = store.add_face(Face {
        surface,
        outer_wire: WireId(0),
        inner_wires: vec![],
        outward: true,
    });

    let mut coedge_ids = Vec::with_capacity(4);
    for (i, &edge_id) in edges.iter().enumerate() {
        let pcurve = Curve2::Line {
            origin: Point2::new(i as f64 / 4.0, 0.0),
            dir: Vector2::new(0.25, 0.0),
        };
        coedge_ids.push(store.add_coedge(CoEdge {
            edge: edge_id,
            forward: true,
            pcurve,
            next: CoEdgeId(0),
            face: face_id,
        }));
    }

    for i in 0..4 {
        store.coedge_mut(coedge_ids[i]).next = coedge_ids[(i + 1) % 4];
    }

    let wire = store.add_wire(Wire { first_coedge: coedge_ids[0] });
    store.face_mut(face_id).outer_wire = wire;
    face_id
}

/// Build a fan-triangulated planar face from a vertex ring.
fn build_fan_face(store: &mut TopoStore, vids: &[VertexId], outward_up: bool) -> FaceId {
    use crate::curve::{Curve2, Curve3};
    use crate::math::{Point2, Vector2};

    let n = vids.len();
    let p0 = store.vertex(vids[0]).point;
    let p1 = store.vertex(vids[1]).point;
    let p_last = store.vertex(vids[n - 1]).point;

    let e1 = p1 - p0;
    let e2 = p_last - p0;
    let mut normal = e1.cross(&e2);
    let len = normal.norm();
    if len > 1e-15 {
        normal /= len;
    }
    if !outward_up {
        normal = -normal;
    }

    let surface = crate::surface::Surface::Plane {
        origin: p0,
        normal,
    };

    // Create edges around the polygon
    let mut edges = Vec::with_capacity(n);
    for i in 0..n {
        let i_next = (i + 1) % n;
        let pa = store.vertex(vids[i]).point;
        let pb = store.vertex(vids[i_next]).point;
        edges.push(store.add_edge(Edge {
            curve: Curve3::Line { origin: pa, dir: pb - pa },
            t_start: 0.0, t_end: 1.0,
            start: vids[i], end: vids[i_next],
        }));
    }

    let face_id = store.add_face(Face {
        surface,
        outer_wire: WireId(0),
        inner_wires: vec![],
        outward: outward_up,
    });

    let mut coedge_ids = Vec::with_capacity(n);
    for (i, &edge_id) in edges.iter().enumerate() {
        let pcurve = Curve2::Line {
            origin: Point2::new(i as f64 / n as f64, 0.0),
            dir: Vector2::new(1.0 / n as f64, 0.0),
        };
        coedge_ids.push(store.add_coedge(CoEdge {
            edge: edge_id,
            forward: outward_up,
            pcurve,
            next: CoEdgeId(0),
            face: face_id,
        }));
    }

    for i in 0..n {
        store.coedge_mut(coedge_ids[i]).next = coedge_ids[(i + 1) % n];
    }

    let wire = store.add_wire(Wire { first_coedge: coedge_ids[0] });
    store.face_mut(face_id).outer_wire = wire;
    face_id
}
