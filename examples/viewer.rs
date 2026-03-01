//! B-Rep Kernel Test Server
//!
//! Serves a Three.js viewer that showcases every capability of the Crusst v1.4
//! B-Rep kernel: all primitives, transforms, profile operations, and exports.
//!
//! Run: `cargo run --example viewer`
//! Open: http://localhost:8080

use crusst::builder::Shape;
use crusst::profile::Profile;
use crusst::types::TessSettings;
use std::f64::consts::TAU;
use tiny_http::{Header, Response, Server};

/// All demo shapes the server can produce.
fn build_shape(name: &str) -> Option<Shape> {
    Some(match name {
        // ─── Primitives ───
        "box"      => Shape::box3(5.0, 3.0, 8.0),
        "sphere"   => Shape::sphere(10.0),
        "cylinder" => Shape::cylinder(5.0, 20.0),
        "cone"     => Shape::cone(8.0, 2.0, 20.0),
        "torus"    => Shape::torus(10.0, 3.0),
        "wedge"    => Shape::wedge(5.0, 4.0, 10.0),
        "capsule"  => Shape::capsule(3.0, 12.0),

        // ─── Transforms ───
        "translate" => Shape::box3(5.0, 3.0, 8.0).translate(10.0, 5.0, 0.0),
        "rotate"    => Shape::box3(8.0, 2.0, 2.0).rotate_z(std::f64::consts::FRAC_PI_4),
        "scale"     => Shape::box3(3.0, 3.0, 3.0).scale(3.0),
        "mirror"    => Shape::box3(5.0, 3.0, 8.0).translate(10.0, 0.0, 0.0).mirror_x(),

        // ─── Profiles ───
        "extrude_rect" => Shape::extrude(&Profile::rect(6.0, 4.0), 15.0),
        "extrude_tri"  => Shape::extrude(
            &Profile::polygon(&[
                nalgebra::Point2::new(0.0, 0.0),
                nalgebra::Point2::new(8.0, 0.0),
                nalgebra::Point2::new(4.0, 6.0),
            ]),
            12.0,
        ),
        "revolve_rect" => Shape::revolve(&Profile::rect(3.0, 1.5), TAU),
        "revolve_half" => Shape::revolve(&Profile::rect(3.0, 1.5), std::f64::consts::PI),

        // ─── Combined / stress tests ───
        "tall_cylinder" => Shape::cylinder(2.0, 50.0),
        "flat_box"      => Shape::box3(20.0, 20.0, 0.5),
        "thin_torus"    => Shape::torus(15.0, 0.5),
        "fat_torus"     => Shape::torus(6.0, 5.0),
        "tiny_sphere"   => Shape::sphere(0.5),
        "big_sphere"    => Shape::sphere(50.0),
        "needle_cone"   => Shape::cone(10.0, 0.1, 30.0),
        "disk_cone"     => Shape::cone(10.0, 10.0, 1.0),

        _ => return None,
    })
}

/// Shape metadata for the gallery.
struct ShapeMeta {
    name: &'static str,
    label: &'static str,
    category: &'static str,
}

fn shape_gallery() -> Vec<ShapeMeta> {
    vec![
        ShapeMeta { name: "box",        label: "Box",          category: "Primitives" },
        ShapeMeta { name: "sphere",     label: "Sphere",       category: "Primitives" },
        ShapeMeta { name: "cylinder",   label: "Cylinder",     category: "Primitives" },
        ShapeMeta { name: "cone",       label: "Cone",         category: "Primitives" },
        ShapeMeta { name: "torus",      label: "Torus",        category: "Primitives" },
        ShapeMeta { name: "wedge",      label: "Wedge",        category: "Primitives" },
        ShapeMeta { name: "capsule",    label: "Capsule",      category: "Primitives" },

        ShapeMeta { name: "translate",  label: "Translate",    category: "Transforms" },
        ShapeMeta { name: "rotate",     label: "Rotate Z 45°", category: "Transforms" },
        ShapeMeta { name: "scale",      label: "Scale 3×",     category: "Transforms" },
        ShapeMeta { name: "mirror",     label: "Mirror X",     category: "Transforms" },

        ShapeMeta { name: "extrude_rect", label: "Extrude Rect", category: "Profiles" },
        ShapeMeta { name: "extrude_tri",  label: "Extrude Tri",  category: "Profiles" },
        ShapeMeta { name: "revolve_rect", label: "Revolve 360°", category: "Profiles" },
        ShapeMeta { name: "revolve_half", label: "Revolve 180°", category: "Profiles" },

        ShapeMeta { name: "tall_cylinder", label: "Tall Cylinder",  category: "Stress" },
        ShapeMeta { name: "flat_box",      label: "Flat Box",       category: "Stress" },
        ShapeMeta { name: "thin_torus",    label: "Thin Torus",     category: "Stress" },
        ShapeMeta { name: "fat_torus",     label: "Fat Torus",      category: "Stress" },
        ShapeMeta { name: "tiny_sphere",   label: "Tiny Sphere",    category: "Stress" },
        ShapeMeta { name: "big_sphere",    label: "Big Sphere",     category: "Stress" },
        ShapeMeta { name: "needle_cone",   label: "Needle Cone",    category: "Stress" },
        ShapeMeta { name: "disk_cone",     label: "Disk (r1=r2)",   category: "Stress" },
    ]
}

fn parse_settings(query: &str) -> TessSettings {
    let mut settings = TessSettings {
        chord_tolerance: 0.02,
        max_edge_length: 5.0,
        min_subdivisions: 8,
    };
    for pair in query.split('&') {
        let mut kv = pair.splitn(2, '=');
        let key = kv.next().unwrap_or("");
        let val = kv.next().unwrap_or("");
        match key {
            "chord" => if let Ok(v) = val.parse::<f64>() { settings.chord_tolerance = v; },
            "maxedge" => if let Ok(v) = val.parse::<f64>() { settings.max_edge_length = v; },
            "minsub" => if let Ok(v) = val.parse::<u32>() { settings.min_subdivisions = v; },
            _ => {}
        }
    }
    settings
}

fn main() {
    let addr = "0.0.0.0:8080";
    let server = Server::http(addr).expect("Failed to start HTTP server");
    println!("╔════════════════════════════════════════════╗");
    println!("║   Crusst v1.4 B-Rep Kernel Test Server    ║");
    println!("║   http://localhost:8080                    ║");
    println!("║   Press Ctrl+C to stop                    ║");
    println!("╚════════════════════════════════════════════╝");

    for request in server.incoming_requests() {
        let url = request.url().to_string();
        let (path, query) = url.split_once('?').unwrap_or((&url, ""));

        match path {
            "/" => {
                let html = include_str!("../viewer/index.html");
                let header = Header::from_bytes("Content-Type", "text/html; charset=utf-8").unwrap();
                let _ = request.respond(Response::from_string(html).with_header(header));
            }

            "/api/gallery" => {
                // Return the shape gallery as JSON
                let gallery = shape_gallery();
                let json: String = format!("[{}]", gallery.iter().map(|m| {
                    format!(r#"{{"name":"{}","label":"{}","category":"{}"}}"#, m.name, m.label, m.category)
                }).collect::<Vec<_>>().join(","));
                let header = Header::from_bytes("Content-Type", "application/json").unwrap();
                let _ = request.respond(Response::from_string(json).with_header(header));
            }

            p if p.starts_with("/api/mesh/") => {
                let name = &p["/api/mesh/".len()..];
                match build_shape(name) {
                    Some(shape) => {
                        let mesh = if query.is_empty() {
                            shape.auto_mesh()
                        } else {
                            shape.mesh(&parse_settings(query))
                        };
                        let binary = mesh.to_binary();
                        let header = Header::from_bytes("Content-Type", "application/octet-stream").unwrap();
                        let _ = request.respond(Response::from_data(binary).with_header(header));
                    }
                    None => {
                        let _ = request.respond(Response::from_string("Shape not found").with_status_code(404));
                    }
                }
            }

            p if p.starts_with("/api/info/") => {
                let name = &p["/api/info/".len()..];
                match build_shape(name) {
                    Some(shape) => {
                        let validation = shape.validate();
                        let mesh = if query.is_empty() {
                            shape.auto_mesh()
                        } else {
                            shape.mesh(&parse_settings(query))
                        };

                        // Topology counts
                        let n_verts = shape.store.solid_vertices(shape.solid).len();
                        let n_edges = shape.store.solid_edges(shape.solid).len();
                        let n_faces = shape.store.solid_face_count(shape.solid);

                        // Surface types used
                        let shell = shape.store.shell(shape.store.solid(shape.solid).outer_shell);
                        let mut surf_types: Vec<String> = Vec::new();
                        for &face_id in &shell.faces {
                            let face = shape.store.face(face_id);
                            let t = match &face.surface {
                                crusst::surface::Surface::Plane { .. } => "Plane",
                                crusst::surface::Surface::Cylinder { .. } => "Cylinder",
                                crusst::surface::Surface::Cone { .. } => "Cone",
                                crusst::surface::Surface::Sphere { .. } => "Sphere",
                                crusst::surface::Surface::Torus { .. } => "Torus",
                                crusst::surface::Surface::NurbsSurface(_) => "NURBS",
                            };
                            if !surf_types.contains(&t.to_string()) {
                                surf_types.push(t.to_string());
                            }
                        }

                        let genus_str = match validation.genus {
                            Some(g) => format!("{}", g),
                            None => "null".to_string(),
                        };

                        let json = format!(
                            r#"{{"valid":{},"errors":[{}],"vertices":{},"edges":{},"faces":{},"euler":{},"genus":{},"surfaces":[{}],"mesh_vertices":{},"mesh_triangles":{}}}"#,
                            validation.valid,
                            validation.errors.iter().map(|e| format!("\"{}\"", e.replace('"', "\\\""))).collect::<Vec<_>>().join(","),
                            n_verts,
                            n_edges,
                            n_faces,
                            validation.euler,
                            genus_str,
                            surf_types.iter().map(|s| format!("\"{}\"", s)).collect::<Vec<_>>().join(","),
                            mesh.vertices.len(),
                            mesh.indices.len() / 3,
                        );
                        let header = Header::from_bytes("Content-Type", "application/json").unwrap();
                        let _ = request.respond(Response::from_string(json).with_header(header));
                    }
                    None => {
                        let _ = request.respond(Response::from_string("Shape not found").with_status_code(404));
                    }
                }
            }

            p if p.starts_with("/api/export/") => {
                let rest = &p["/api/export/".len()..];
                let (format, name) = rest.split_once('/').unwrap_or(("", rest));
                let settings = parse_settings(query);
                match build_shape(name) {
                    Some(shape) => {
                        let mut buf = Vec::new();
                        let (ct, filename) = match format {
                            "stl" => {
                                shape.write_stl(&settings, &mut buf).unwrap();
                                ("application/octet-stream", format!("{name}.stl"))
                            }
                            "obj" => {
                                shape.write_obj(&settings, &mut buf).unwrap();
                                ("text/plain", format!("{name}.obj"))
                            }
                            "step" => {
                                shape.write_step(&mut buf).unwrap();
                                ("text/plain", format!("{name}.stp"))
                            }
                            "3mf" => {
                                shape.write_3mf(&settings, &mut buf).unwrap();
                                ("application/xml", format!("{name}.3mf"))
                            }
                            _ => {
                                let _ = request.respond(Response::from_string("Unknown format").with_status_code(400));
                                continue;
                            }
                        };
                        let ct_header = Header::from_bytes("Content-Type", ct).unwrap();
                        let disp_header = Header::from_bytes(
                            "Content-Disposition",
                            format!("attachment; filename=\"{}\"", filename)
                        ).unwrap();
                        let _ = request.respond(
                            Response::from_data(buf)
                                .with_header(ct_header)
                                .with_header(disp_header)
                        );
                    }
                    None => {
                        let _ = request.respond(Response::from_string("Shape not found").with_status_code(404));
                    }
                }
            }

            _ => {
                let _ = request.respond(Response::from_string("404").with_status_code(404));
            }
        }
    }
}
