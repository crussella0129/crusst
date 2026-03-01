//! Boundary-first tessellation of B-Rep faces and solids.
//!
//! Two strategies based on surface type:
//! - **Planar faces**: Discretize wire boundary → ear-clip triangulation.
//!   Correctly handles triangles, circles, arbitrary polygons.
//! - **Curved faces**: Adaptive parametric grid tessellation.
//!   Handles periodic surfaces (cylinder, sphere, torus) robustly.

mod adaptive;

use crate::curve::Curve3;
use crate::math::{Point3, Vector3};
use crate::surface::Surface;
use crate::topo::*;
use crate::types::{TessSettings, TriangleMesh};
use nalgebra::Vector3 as NVec3;

/// Tessellate an entire solid into a single triangle mesh.
pub fn tessellate_solid(store: &TopoStore, solid_id: SolidId, settings: &TessSettings) -> TriangleMesh {
    let solid = store.solid(solid_id);
    let shell = store.shell(solid.outer_shell);

    let mut all_verts: Vec<NVec3<f64>> = Vec::new();
    let mut all_normals: Vec<NVec3<f64>> = Vec::new();
    let mut all_indices: Vec<u32> = Vec::new();

    for &face_id in &shell.faces {
        let face = store.face(face_id);
        let (verts, normals, indices) = if matches!(&face.surface, Surface::Plane { .. }) {
            tessellate_planar_face(store, face_id, &face.surface, face.outward, settings)
        } else {
            tessellate_curved_face(store, face_id, &face.surface, face.outward, settings)
        };

        let base = all_verts.len() as u32;
        all_verts.extend_from_slice(&verts);
        all_normals.extend_from_slice(&normals);
        for &idx in &indices {
            all_indices.push(base + idx);
        }
    }

    TriangleMesh {
        vertices: all_verts,
        normals: all_normals,
        indices: all_indices,
    }
}

/// Tessellate a single face, returning (vertices, normals, triangle_indices).
pub fn tessellate_face(
    store: &TopoStore,
    face_id: FaceId,
    settings: &TessSettings,
) -> (Vec<Point3>, Vec<Vector3>, Vec<[u32; 3]>) {
    let face = store.face(face_id);
    let (verts, normals, flat_indices) = if matches!(&face.surface, Surface::Plane { .. }) {
        tessellate_planar_face(store, face_id, &face.surface, face.outward, settings)
    } else {
        tessellate_curved_face(store, face_id, &face.surface, face.outward, settings)
    };

    let pts: Vec<Point3> = verts.iter().map(|v| Point3::new(v.x, v.y, v.z)).collect();
    let norms: Vec<Vector3> = normals.to_vec();
    let tris: Vec<[u32; 3]> = flat_indices.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

    (pts, norms, tris)
}

// ─── Planar face tessellation (ear clipping) ────────────────────────────────

/// Tessellate a planar face by discretizing its wire boundary and ear-clipping.
fn tessellate_planar_face(
    store: &TopoStore,
    face_id: FaceId,
    surface: &Surface,
    outward: bool,
    settings: &TessSettings,
) -> (Vec<NVec3<f64>>, Vec<NVec3<f64>>, Vec<u32>) {
    // Discretize wire boundary into 3D points
    let boundary_3d = discretize_wire_boundary_3d(store, face_id, settings);
    if boundary_3d.len() < 3 {
        return (vec![], vec![], vec![]);
    }

    // Project to (u,v)
    let boundary_uv: Vec<(f64, f64)> = boundary_3d.iter()
        .map(|pt| surface.closest_parameters(pt))
        .collect();

    // Compute constant normal for plane
    let mut normal = surface.normal(boundary_uv[0].0, boundary_uv[0].1);
    if !outward { normal = -normal; }

    // Map to 3D vertex positions (use evaluated points for consistency)
    let mut verts: Vec<NVec3<f64>> = Vec::with_capacity(boundary_3d.len());
    let mut normals: Vec<NVec3<f64>> = Vec::with_capacity(boundary_3d.len());
    for pt in &boundary_3d {
        verts.push(NVec3::new(pt.x, pt.y, pt.z));
        normals.push(normal);
    }

    // Ear clipping triangulation
    let indices = ear_clip(&boundary_uv, outward);

    (verts, normals, indices)
}

// ─── Curved face tessellation (adaptive parametric grid) ────────────────────

/// Tessellate a curved face using adaptive parametric grid subdivision.
fn tessellate_curved_face(
    store: &TopoStore,
    face_id: FaceId,
    surface: &Surface,
    outward: bool,
    settings: &TessSettings,
) -> (Vec<NVec3<f64>>, Vec<NVec3<f64>>, Vec<u32>) {
    let (u_range, v_range) = compute_face_domain(store, face_id, surface);
    let base_n = settings.min_subdivisions as usize;
    adaptive::adaptive_tessellate(
        surface,
        u_range,
        v_range,
        base_n,
        settings.chord_tolerance,
        settings.max_edge_length,
        outward,
    )
}

/// Compute the (u,v) parameter domain for a curved face.
fn compute_face_domain(store: &TopoStore, face_id: FaceId, surface: &Surface) -> ((f64, f64), (f64, f64)) {
    use std::f64::consts::{FRAC_PI_2, TAU};

    match surface {
        Surface::Cylinder { .. } | Surface::Cone { .. } => {
            let samples = collect_face_uv_samples(store, face_id, surface);
            if samples.is_empty() {
                return ((0.0, TAU), (0.0, 1.0));
            }
            let mut v_min = f64::MAX;
            let mut v_max = f64::MIN;
            for &(_, v) in &samples {
                v_min = v_min.min(v);
                v_max = v_max.max(v);
            }
            let u_angles: Vec<f64> = samples.iter().map(|&(u, _)| u).collect();
            let u_range = compute_periodic_range(&u_angles, TAU);
            (u_range, (v_min, v_max))
        }
        Surface::Sphere { .. } => ((0.0, TAU), (-FRAC_PI_2, FRAC_PI_2)),
        Surface::Torus { .. } => ((0.0, TAU), (0.0, TAU)),
        Surface::NurbsSurface(nurbs) => (nurbs.domain_u(), nurbs.domain_v()),
        Surface::Plane { .. } => unreachable!("Plane faces use ear clipping"),
    }
}

// ─── Wire boundary discretization ──────────────────────────────────────────

/// Discretize a face's outer wire into 3D boundary points.
/// Circular arcs are subdivided adaptively; lines are just endpoints.
fn discretize_wire_boundary_3d(
    store: &TopoStore,
    face_id: FaceId,
    settings: &TessSettings,
) -> Vec<Point3> {
    let face = store.face(face_id);
    let coedge_ids = store.wire_coedges(face.outer_wire);
    let mut pts = Vec::new();

    for &ce_id in &coedge_ids {
        let ce = store.coedge(ce_id);
        let edge = store.edge(ce.edge);

        let (start_vid, end_vid) = if ce.forward {
            (edge.start, edge.end)
        } else {
            (edge.end, edge.start)
        };
        let p_start = store.vertex(start_vid).point;
        let p_end = store.vertex(end_vid).point;

        match &edge.curve {
            Curve3::Circle { center, axis, radius } => {
                // Adaptively subdivide arc
                let n_segs = arc_subdivision_count(*radius, settings.chord_tolerance)
                    .max(settings.min_subdivisions as usize / 2)
                    .max(4);

                let arc_pts = sample_arc_3d(p_start, p_end, *center, *axis, *radius, n_segs);
                // Add all but last point (next edge starts there)
                for i in 0..arc_pts.len() - 1 {
                    pts.push(arc_pts[i]);
                }
            }
            _ => {
                // Line edge: just the start vertex
                pts.push(p_start);
            }
        }
    }

    pts
}

/// Sample points along a circular arc from p_start to p_end.
fn sample_arc_3d(
    p_start: Point3, p_end: Point3,
    center: Point3, axis: Vector3, radius: f64,
    n_segs: usize,
) -> Vec<Point3> {
    let a = axis.normalize();
    // Build a local frame for the circle
    let d_start = p_start - center;
    let d_start_radial = d_start - a * d_start.dot(&a);
    let axial_offset = d_start.dot(&a);

    let r_len = d_start_radial.norm();
    if r_len < 1e-15 {
        return vec![p_start, p_end];
    }
    let e1 = d_start_radial / r_len;
    let e2 = a.cross(&e1);

    // Find angle to end point
    let d_end = p_end - center;
    let d_end_radial = d_end - a * d_end.dot(&a);
    let angle_start = 0.0; // e1 direction
    let mut angle_end = d_end_radial.dot(&e2).atan2(d_end_radial.dot(&e1));

    // Ensure we go the short way around (< π for half-circles, otherwise correct)
    if angle_end <= 1e-10 {
        angle_end += std::f64::consts::TAU;
    }

    let mut pts = Vec::with_capacity(n_segs + 1);
    for i in 0..=n_segs {
        let t = angle_start + (angle_end - angle_start) * i as f64 / n_segs as f64;
        let pt = center + (e1 * t.cos() + e2 * t.sin()) * radius + a * axial_offset;
        pts.push(pt);
    }
    pts
}

/// Compute arc subdivision count from chord tolerance.
fn arc_subdivision_count(radius: f64, tolerance: f64) -> usize {
    if radius < 1e-15 || tolerance <= 0.0 { return 8; }
    let ratio = (tolerance / radius).min(1.0);
    let max_angle = 2.0 * (1.0 - ratio).acos();
    if max_angle < 1e-10 { return 64; }
    let total_angle = std::f64::consts::PI; // half-circle
    ((total_angle / max_angle).ceil() as usize).max(4)
}

// ─── UV sampling for curved face domains ────────────────────────────────────

/// Collect (u, v) parameter samples from a face's wire — both vertex positions
/// AND geometric midpoints of circular arc edges.
fn collect_face_uv_samples(store: &TopoStore, face_id: FaceId, surface: &Surface) -> Vec<(f64, f64)> {
    let face = store.face(face_id);
    let coedge_ids = store.wire_coedges(face.outer_wire);
    let mut samples = Vec::new();

    for &ce_id in &coedge_ids {
        let ce = store.coedge(ce_id);
        let edge = store.edge(ce.edge);

        let vid = if ce.forward { edge.start } else { edge.end };
        let pt = store.vertex(vid).point;
        samples.push(surface.closest_parameters(&pt));

        let p_start = store.vertex(edge.start).point;
        let p_end = store.vertex(edge.end).point;

        match &edge.curve {
            Curve3::Circle { center, axis, radius } => {
                let chord = p_end - p_start;
                let cross = chord.cross(axis);
                let cross_len = cross.norm();
                if cross_len > 1e-10 {
                    let mid_dir = cross / cross_len;
                    let arc_mid = center + mid_dir * *radius;
                    samples.push(surface.closest_parameters(&arc_mid));
                }
            }
            _ => {
                let mid = nalgebra::center(&p_start, &p_end);
                samples.push(surface.closest_parameters(&mid));
            }
        }
    }

    samples
}

/// Compute the angular range using the "largest gap" algorithm.
fn compute_periodic_range(angles: &[f64], period: f64) -> (f64, f64) {
    if angles.is_empty() { return (0.0, period); }

    let mut normalized: Vec<f64> = angles.iter().map(|&a| {
        let mut v = a % period;
        if v < 0.0 { v += period; }
        v
    }).collect();
    normalized.sort_by(|a, b| a.partial_cmp(b).unwrap());
    normalized.dedup_by(|a, b| (*a - *b).abs() < 1e-10);

    if normalized.len() <= 1 { return (0.0, period); }

    let mut max_gap = 0.0_f64;
    let mut gap_end_idx = 0;
    for i in 0..normalized.len() {
        let next = (i + 1) % normalized.len();
        let gap = if next > i {
            normalized[next] - normalized[i]
        } else {
            (period - normalized[i]) + normalized[0]
        };
        if gap > max_gap { max_gap = gap; gap_end_idx = next; }
    }

    let range_start = normalized[gap_end_idx];
    let gap_start_idx = if gap_end_idx == 0 { normalized.len() - 1 } else { gap_end_idx - 1 };
    let range_end = normalized[gap_start_idx];

    if range_end >= range_start {
        (range_start, range_end)
    } else {
        (range_start, range_end + period)
    }
}

// ─── Ear clipping triangulation ─────────────────────────────────────────────

fn ear_clip(polygon: &[(f64, f64)], outward: bool) -> Vec<u32> {
    let n = polygon.len();
    if n < 3 { return vec![]; }

    let area = signed_area(polygon);
    let ccw = area > 0.0;

    // For any polygon (including triangles), winding must account for
    // both the 2D polygon orientation (ccw) AND the outward flag.
    // Keep order when both agree; reverse when they disagree.
    let keep_order = (outward && ccw) || (!outward && !ccw);

    if n == 3 {
        return if keep_order { vec![0, 1, 2] } else { vec![0, 2, 1] };
    }

    let mut indices: Vec<usize> = (0..n).collect();
    let mut result = Vec::with_capacity((n - 2) * 3);

    let mut safety = n * n;
    while indices.len() > 3 && safety > 0 {
        safety -= 1;
        let len = indices.len();
        let mut found_ear = false;

        for i in 0..len {
            let prev = indices[(i + len - 1) % len];
            let curr = indices[i];
            let next = indices[(i + 1) % len];

            if is_ear(polygon, &indices, prev, curr, next, ccw) {
                if keep_order {
                    result.extend_from_slice(&[prev as u32, curr as u32, next as u32]);
                } else {
                    result.extend_from_slice(&[prev as u32, next as u32, curr as u32]);
                }
                indices.remove(i);
                found_ear = true;
                break;
            }
        }

        if !found_ear { break; }
    }

    if indices.len() == 3 {
        let (a, b, c) = (indices[0], indices[1], indices[2]);
        if keep_order {
            result.extend_from_slice(&[a as u32, b as u32, c as u32]);
        } else {
            result.extend_from_slice(&[a as u32, c as u32, b as u32]);
        }
    }

    result
}

fn signed_area(polygon: &[(f64, f64)]) -> f64 {
    let n = polygon.len();
    let mut area = 0.0;
    for i in 0..n {
        let (x1, y1) = polygon[i];
        let (x2, y2) = polygon[(i + 1) % n];
        area += x1 * y2 - x2 * y1;
    }
    area * 0.5
}

fn is_ear(polygon: &[(f64, f64)], indices: &[usize], prev: usize, curr: usize, next: usize, ccw: bool) -> bool {
    let (ax, ay) = polygon[prev];
    let (bx, by) = polygon[curr];
    let (cx, cy) = polygon[next];

    let cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
    if (ccw && cross <= 0.0) || (!ccw && cross >= 0.0) {
        return false;
    }

    for &idx in indices {
        if idx == prev || idx == curr || idx == next { continue; }
        let (px, py) = polygon[idx];
        if point_in_triangle(px, py, ax, ay, bx, by, cx, cy) {
            return false;
        }
    }

    true
}

fn point_in_triangle(px: f64, py: f64, ax: f64, ay: f64, bx: f64, by: f64, cx: f64, cy: f64) -> bool {
    let d1 = (px - bx) * (ay - by) - (ax - bx) * (py - by);
    let d2 = (px - cx) * (by - cy) - (bx - cx) * (py - cy);
    let d3 = (px - ax) * (cy - ay) - (cx - ax) * (py - ay);
    let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
    let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);
    !(has_neg && has_pos)
}
