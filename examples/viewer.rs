use nalgebra::Rotation3;
use nalgebra::Unit;
use nalgebra::{Vector2, Vector3};
use std::f64::consts::PI;
use std::path::PathBuf;

use crusst::blend;
use crusst::builder::Shape;
use crusst::feature::ft;
use crusst::mesh::extract_mesh;
use crusst::obj_export::write_obj;
use crusst::shape::*;
use crusst::types::MeshSettings;

use tiny_http::{Header, Response, Server};

/// Generate all showcase OBJ files into the output directory.
fn generate_models(dir: &PathBuf) {
    println!("Generating models...");
    std::fs::create_dir_all(dir).unwrap();

    #[allow(clippy::type_complexity)]
    let models: Vec<(&str, Box<dyn Sdf>, Vector3<f64>, Vector3<f64>)> = vec![
        // Primitives
        (
            "01_sphere",
            Box::new(Sphere::new(Vector3::zeros(), 10.0)),
            Vector3::from_element(-12.0),
            Vector3::from_element(12.0),
        ),
        (
            "02_box",
            Box::new(Box3::new(Vector3::zeros(), Vector3::new(5.0, 3.0, 8.0))),
            Vector3::new(-7.0, -5.0, -10.0),
            Vector3::new(7.0, 5.0, 10.0),
        ),
        (
            "03_cylinder",
            Box::new(Cylinder::new(
                Vector3::zeros(),
                Vector3::new(0.0, 0.0, 1.0),
                5.0,
                20.0,
            )),
            Vector3::new(-7.0, -7.0, -2.0),
            Vector3::new(7.0, 7.0, 22.0),
        ),
        (
            "04_capped_cone",
            Box::new(CappedCone::new(
                Vector3::zeros(),
                Vector3::new(0.0, 0.0, 20.0),
                8.0,
                2.0,
            )),
            Vector3::new(-10.0, -10.0, -2.0),
            Vector3::new(10.0, 10.0, 22.0),
        ),
        (
            "05_torus",
            Box::new(Torus::new(Vector3::zeros(), 10.0, 3.0)),
            Vector3::new(-15.0, -5.0, -15.0),
            Vector3::new(15.0, 5.0, 15.0),
        ),
        (
            "06_rounded_box",
            Box::new(RoundedBox::new(
                Vector3::zeros(),
                Vector3::new(5.0, 3.0, 8.0),
                1.0,
            )),
            Vector3::new(-8.0, -6.0, -11.0),
            Vector3::new(8.0, 6.0, 11.0),
        ),
        (
            "07_capsule",
            Box::new(Capsule::new(
                Vector3::zeros(),
                Vector3::new(0.0, 0.0, 20.0),
                4.0,
            )),
            Vector3::new(-6.0, -6.0, -6.0),
            Vector3::new(6.0, 6.0, 26.0),
        ),
        (
            "08_ellipsoid",
            Box::new(Ellipsoid::new(
                Vector3::zeros(),
                Vector3::new(10.0, 5.0, 3.0),
            )),
            Vector3::new(-12.0, -7.0, -5.0),
            Vector3::new(12.0, 7.0, 5.0),
        ),
        (
            "09_rounded_cylinder",
            Box::new(RoundedCylinder::new(Vector3::zeros(), 6.0, 1.0, 10.0)),
            Vector3::new(-9.0, -13.0, -9.0),
            Vector3::new(9.0, 13.0, 9.0),
        ),
        // Boolean operations
        (
            "10_boolean_union",
            Box::new(Union::new(
                Sphere::new(Vector3::new(-3.0, 0.0, 0.0), 5.0),
                Sphere::new(Vector3::new(3.0, 0.0, 0.0), 5.0),
            )),
            Vector3::from_element(-10.0),
            Vector3::from_element(10.0),
        ),
        (
            "11_boolean_intersection",
            Box::new(Intersection::new(
                Sphere::new(Vector3::new(-2.0, 0.0, 0.0), 5.0),
                Sphere::new(Vector3::new(2.0, 0.0, 0.0), 5.0),
            )),
            Vector3::from_element(-8.0),
            Vector3::from_element(8.0),
        ),
        (
            "12_boolean_difference",
            Box::new(Difference::new(
                Box3::new(Vector3::zeros(), Vector3::new(10.0, 10.0, 10.0)),
                Sphere::new(Vector3::zeros(), 7.0),
            )),
            Vector3::from_element(-12.0),
            Vector3::from_element(12.0),
        ),
        // Smooth booleans
        (
            "13_smooth_union",
            Box::new(SmoothUnion::new(
                Sphere::new(Vector3::new(-4.0, 0.0, 0.0), 5.0),
                Sphere::new(Vector3::new(4.0, 0.0, 0.0), 5.0),
                2.0,
            )),
            Vector3::from_element(-11.0),
            Vector3::from_element(11.0),
        ),
        (
            "14_smooth_intersection",
            Box::new(SmoothIntersection::new(
                Sphere::new(Vector3::new(-2.0, 0.0, 0.0), 6.0),
                Sphere::new(Vector3::new(2.0, 0.0, 0.0), 6.0),
                2.0,
            )),
            Vector3::from_element(-8.0),
            Vector3::from_element(8.0),
        ),
        (
            "15_smooth_difference",
            Box::new(SmoothDifference::new(
                Box3::new(Vector3::zeros(), Vector3::new(8.0, 8.0, 8.0)),
                Sphere::new(Vector3::zeros(), 6.0),
                1.5,
            )),
            Vector3::from_element(-10.0),
            Vector3::from_element(10.0),
        ),
        // Transforms — each shows original (at origin) + transformed copy side by side

        // Translate: original box at origin + translated copy offset to the right
        (
            "16_translate",
            Box::new(Union::new(
                Box3::new(Vector3::zeros(), Vector3::new(4.0, 4.0, 4.0)),
                Translate::new(
                    Box3::new(Vector3::zeros(), Vector3::new(4.0, 4.0, 4.0)),
                    Vector3::new(18.0, 0.0, 0.0),
                ),
            )),
            Vector3::new(-6.0, -6.0, -6.0),
            Vector3::new(24.0, 6.0, 6.0),
        ),
        // Rotate: two boxes unioned at 45° — demonstrates the Rotate transform
        (
            "17_rotate",
            Box::new(Union::new(
                Box3::new(Vector3::zeros(), Vector3::new(12.0, 2.0, 2.0)),
                Rotate::new(
                    Box3::new(Vector3::zeros(), Vector3::new(12.0, 2.0, 2.0)),
                    Rotation3::from_axis_angle(
                        &Unit::new_normalize(Vector3::new(0.0, 0.0, 1.0)),
                        PI / 4.0,
                    ),
                ),
            )),
            Vector3::new(-14.0, -14.0, -4.0),
            Vector3::new(14.0, 14.0, 4.0),
        ),
        // Scale: small sphere at origin + large scaled sphere offset so both visible
        (
            "18_scale",
            Box::new(Union::new(
                Sphere::new(Vector3::new(-10.0, 0.0, 0.0), 4.0),
                Translate::new(
                    Scale::new(Sphere::new(Vector3::zeros(), 4.0), 2.5),
                    Vector3::new(10.0, 0.0, 0.0),
                ),
            )),
            Vector3::new(-16.0, -12.0, -12.0),
            Vector3::new(22.0, 12.0, 12.0),
        ),
        // Mirror: L-bracket (asymmetric!) + its mirror — chirality is clearly visible
        // The L-bracket is a union of two boxes forming an L shape
        (
            "19_mirror",
            {
                // L-bracket: vertical arm + horizontal foot (asymmetric shape)
                let l_bracket = Union::new(
                    Box3::new(Vector3::new(7.0, 5.0, 0.0), Vector3::new(2.0, 5.0, 2.0)), // vertical arm
                    Box3::new(Vector3::new(10.0, 1.0, 0.0), Vector3::new(5.0, 1.0, 2.0)), // horizontal foot
                );
                // Show original + mirror across YZ plane (x=0)
                let mirrored: Box<dyn Sdf> = Box::new(Mirror::new(
                    Union::new(
                        Box3::new(Vector3::new(7.0, 5.0, 0.0), Vector3::new(2.0, 5.0, 2.0)),
                        Box3::new(Vector3::new(10.0, 1.0, 0.0), Vector3::new(5.0, 1.0, 2.0)),
                    ),
                    Vector3::new(1.0, 0.0, 0.0),
                ));
                let combined: Box<dyn Sdf> = Box::new(Union::new(
                    l_bracket,
                    FnSdf::new(move |p| mirrored.evaluate(p)),
                ));
                combined
            },
            Vector3::new(-17.0, -2.0, -4.0),
            Vector3::new(17.0, 12.0, 4.0),
        ),
        (
            "20_shell",
            Box::new(Shell::new(Sphere::new(Vector3::zeros(), 10.0), 1.0)),
            Vector3::from_element(-13.0),
            Vector3::from_element(13.0),
        ),
        // Revolve: a 2D L-profile revolved around Y axis → a vase/cup shape
        (
            "20b_revolve",
            {
                // 2D profile: rectangle forming the wall of a cup
                // Outer wall at x=8..10, y=0..15 (r=8..10, height 15)
                // Bottom at x=0..10, y=-1..0
                let profile = Union2d::new(
                    Rect2d::new(Vector2::new(9.0, 7.5), Vector2::new(1.0, 7.5)), // wall
                    Rect2d::new(Vector2::new(5.0, -0.5), Vector2::new(5.0, 0.5)), // bottom
                );
                let cup: Box<dyn Sdf> = Box::new(Revolve::new(profile));
                cup
            },
            Vector3::new(-13.0, -3.0, -13.0),
            Vector3::new(13.0, 18.0, 13.0),
        ),
        // Extrude: a 2D cross profile extruded along Z → a plus-sign beam
        (
            "20c_extrude",
            {
                let cross = Union2d::new(
                    Rect2d::new(Vector2::new(0.0, 0.0), Vector2::new(8.0, 2.0)), // horizontal bar
                    Rect2d::new(Vector2::new(0.0, 0.0), Vector2::new(2.0, 8.0)), // vertical bar
                );
                let beam: Box<dyn Sdf> = Box::new(Extrude::new(cross, 12.0));
                beam
            },
            Vector3::new(-10.0, -10.0, -14.0),
            Vector3::new(10.0, 10.0, 14.0),
        ),
        // Composed shapes
        (
            "21_mounting_bracket",
            Box::new(Difference::new(
                Union::new(
                    Union::new(
                        Box3::new(Vector3::zeros(), Vector3::new(20.0, 2.0, 15.0)),
                        Capsule::new(
                            Vector3::new(-12.0, 2.0, 0.0),
                            Vector3::new(-12.0, 15.0, 0.0),
                            3.0,
                        ),
                    ),
                    Capsule::new(
                        Vector3::new(12.0, 2.0, 0.0),
                        Vector3::new(12.0, 15.0, 0.0),
                        3.0,
                    ),
                ),
                Cylinder::new(
                    Vector3::new(0.0, -5.0, 0.0),
                    Vector3::new(0.0, 1.0, 0.0),
                    5.0,
                    14.0,
                ),
            )),
            Vector3::new(-18.0, -4.0, -18.0),
            Vector3::new(18.0, 18.0, 18.0),
        ),
        (
            "22_pipe_tee",
            Box::new(Difference::new(
                SmoothUnion::new(
                    Capsule::new(
                        Vector3::new(-15.0, 0.0, 0.0),
                        Vector3::new(15.0, 0.0, 0.0),
                        5.0,
                    ),
                    Capsule::new(
                        Vector3::new(0.0, 0.0, 0.0),
                        Vector3::new(0.0, 15.0, 0.0),
                        5.0,
                    ),
                    1.5,
                ),
                Union::new(
                    Capsule::new(
                        Vector3::new(-16.0, 0.0, 0.0),
                        Vector3::new(16.0, 0.0, 0.0),
                        3.0,
                    ),
                    Capsule::new(
                        Vector3::new(0.0, -1.0, 0.0),
                        Vector3::new(0.0, 16.0, 0.0),
                        3.0,
                    ),
                ),
            )),
            Vector3::new(-18.0, -3.0, -8.0),
            Vector3::new(18.0, 18.0, 8.0),
        ),
        (
            "23_gasket_ring",
            Box::new(Intersection::new(
                Torus::new(Vector3::zeros(), 12.0, 4.0),
                Box3::new(Vector3::zeros(), Vector3::new(20.0, 2.0, 20.0)),
            )),
            Vector3::new(-18.0, -4.0, -18.0),
            Vector3::new(18.0, 4.0, 18.0),
        ),
        (
            "24_rounded_enclosure",
            Box::new(Difference::new(
                Difference::new(
                    RoundedBox::new(Vector3::zeros(), Vector3::new(15.0, 8.0, 10.0), 2.0),
                    Box3::new(Vector3::new(0.0, 1.0, 0.0), Vector3::new(13.0, 7.0, 8.0)),
                ),
                Cylinder::new(
                    Vector3::new(15.0, 0.0, 0.0),
                    Vector3::new(1.0, 0.0, 0.0),
                    3.0,
                    10.0,
                ),
            )),
            Vector3::new(-19.0, -12.0, -14.0),
            Vector3::new(19.0, 12.0, 14.0),
        ),
        (
            "25_organic_blob",
            Box::new(SmoothUnion::new(
                SmoothUnion::new(
                    Sphere::new(Vector3::new(0.0, 0.0, 0.0), 6.0),
                    Sphere::new(Vector3::new(5.0, 4.0, 0.0), 4.0),
                    3.0,
                ),
                SmoothUnion::new(
                    Sphere::new(Vector3::new(-4.0, 5.0, 3.0), 3.5),
                    Sphere::new(Vector3::new(2.0, -3.0, 5.0), 3.0),
                    3.0,
                ),
                3.0,
            )),
            Vector3::from_element(-12.0),
            Vector3::from_element(12.0),
        ),
        (
            "26_rotated_star",
            Box::new(Union::new(
                Union::new(
                    Box3::new(Vector3::zeros(), Vector3::new(8.0, 2.0, 2.0)),
                    Rotate::new(
                        Box3::new(Vector3::zeros(), Vector3::new(8.0, 2.0, 2.0)),
                        Rotation3::from_axis_angle(
                            &Unit::new_normalize(Vector3::new(0.0, 0.0, 1.0)),
                            PI / 3.0,
                        ),
                    ),
                ),
                Rotate::new(
                    Box3::new(Vector3::zeros(), Vector3::new(8.0, 2.0, 2.0)),
                    Rotation3::from_axis_angle(
                        &Unit::new_normalize(Vector3::new(0.0, 0.0, 1.0)),
                        2.0 * PI / 3.0,
                    ),
                ),
            )),
            Vector3::from_element(-11.0),
            Vector3::from_element(11.0),
        ),
        (
            "27_bearing_housing",
            Box::new(Difference::new(
                Difference::new(
                    Capsule::new(Vector3::zeros(), Vector3::new(0.0, 20.0, 0.0), 12.0),
                    Capsule::new(
                        Vector3::new(0.0, -2.0, 0.0),
                        Vector3::new(0.0, 22.0, 0.0),
                        8.0,
                    ),
                ),
                Translate::new(
                    Rotate::new(
                        Torus::new(Vector3::zeros(), 8.0, 1.5),
                        Rotation3::from_axis_angle(
                            &Unit::new_normalize(Vector3::new(1.0, 0.0, 0.0)),
                            PI / 2.0,
                        ),
                    ),
                    Vector3::new(0.0, 10.0, 0.0),
                ),
            )),
            Vector3::new(-15.0, -3.0, -15.0),
            Vector3::new(15.0, 23.0, 15.0),
        ),
    ];

    for (name, shape, bbox_min, bbox_max) in &models {
        print!("  {} ... ", name);
        let mesh = extract_mesh(shape.as_ref(), *bbox_min, *bbox_max, 128);
        let path = dir.join(format!("{}.obj", name));
        write_obj(&mesh, &path).unwrap();
        println!(
            "{} tris, {} verts",
            mesh.indices.len() / 3,
            mesh.vertices.len()
        );
    }

    // Fillet / chamfer showcase models (using the builder API)
    let settings = MeshSettings::default();
    let builder_models: Vec<(&str, Shape)> = vec![
        // Box3 with G2 fillet on all edges
        (
            "30_filleted_box",
            Shape::box3(5.0, 5.0, 5.0).fillet(blend::g2(1.0), vec![ft(0, 0).all_edges()]),
        ),
        // Box3 with equal chamfer on all edges
        (
            "31_chamfered_box",
            Shape::box3(5.0, 5.0, 5.0)
                .chamfer(blend::equal_chamfer(1.0), vec![ft(0, 0).all_edges()]),
        ),
        // Two spheres with round union
        (
            "32_round_union",
            Shape::sphere(5.0).round_union(Shape::sphere(5.0).translate(7.0, 0.0, 0.0), 1.5),
        ),
        // Box minus cylinder with round fillet
        (
            "33_round_subtract",
            Shape::box3(5.0, 5.0, 5.0)
                .round_subtract(Shape::cylinder(3.0, 15.0).translate(0.0, -2.5, 0.0), 1.0),
        ),
        // Box with fillet only on top face edges (+Y edges: 0, 4, 8, 9)
        (
            "34_selective_fillet",
            Shape::box3(5.0, 5.0, 5.0).fillet(blend::g2(1.0), vec![ft(0, 0).edges(&[0, 4, 8, 9])]),
        ),
        // Box3 with simple round (Minkowski offset)
        (
            "35_rounded_box_simple",
            Shape::box3(5.0, 5.0, 5.0).round(0.5),
        ),
        // ─── Rotated union (builder API — analytical gradients for clean edges) ───
        (
            "17b_rotated_union",
            Shape::box3(12.0, 2.0, 2.0).union(Shape::box3(12.0, 2.0, 2.0).rotate_z(PI / 4.0)),
        ),
        // ─── Profile comparison series: each profile on the same 5×5×5 box ───
        // G1 tangent-continuous fillet
        (
            "36_g1_fillet_box",
            Shape::box3(5.0, 5.0, 5.0).fillet(blend::g1(1.0), vec![ft(0, 0).all_edges()]),
        ),
        // G3 curvature-rate continuous fillet
        (
            "37_g3_fillet_box",
            Shape::box3(5.0, 5.0, 5.0).fillet(blend::g3(1.0), vec![ft(0, 0).all_edges()]),
        ),
        // Parabolic fillet
        (
            "38_parabolic_fillet_box",
            Shape::box3(5.0, 5.0, 5.0).fillet(blend::parabolic(1.0), vec![ft(0, 0).all_edges()]),
        ),
        // Hyperbolic fillet (asymptote 0.5 gives a concave profile)
        (
            "39_hyperbolic_fillet_box",
            Shape::box3(5.0, 5.0, 5.0)
                .fillet(blend::hyperbolic(1.0, 0.5), vec![ft(0, 0).all_edges()]),
        ),
        // Cycloidal fillet
        (
            "40_cycloidal_fillet_box",
            Shape::box3(5.0, 5.0, 5.0).fillet(blend::cycloidal(1.0), vec![ft(0, 0).all_edges()]),
        ),
    ];

    let total_models = models.len() + builder_models.len();

    for (name, shape) in &builder_models {
        print!("  {} ... ", name);
        let mesh = shape.mesh(settings);
        let path = dir.join(format!("{}.obj", name));
        write_obj(&mesh, &path).unwrap();
        println!(
            "{} tris, {} verts",
            mesh.indices.len() / 3,
            mesh.vertices.len()
        );
    }

    println!("Done — {} models generated.\n", total_models);
}

fn content_type(path: &str) -> &'static str {
    if path.ends_with(".html") {
        "text/html; charset=utf-8"
    } else if path.ends_with(".js") {
        "application/javascript"
    } else if path.ends_with(".css") {
        "text/css"
    } else if path.ends_with(".obj") {
        "text/plain"
    } else {
        "application/octet-stream"
    }
}

fn main() {
    let out_dir = PathBuf::from("viewer");

    // Generate OBJ files
    generate_models(&out_dir);

    let addr = "0.0.0.0:8080";
    let server = Server::http(addr).expect("Failed to start server");
    println!("Crusst Viewer running at http://localhost:8080");
    println!("Press Ctrl+C to stop.\n");

    for request in server.incoming_requests() {
        let url = request.url().to_string();
        let file_path = if url == "/" {
            out_dir.join("index.html")
        } else {
            out_dir.join(url.trim_start_matches('/'))
        };

        if file_path.exists() && file_path.is_file() {
            let data = std::fs::read(&file_path).unwrap();
            let ct = content_type(&file_path.to_string_lossy());
            let header = Header::from_bytes("Content-Type", ct).unwrap();
            let response = Response::from_data(data).with_header(header);
            let _ = request.respond(response);
        } else {
            let response = Response::from_string("404 Not Found").with_status_code(404);
            let _ = request.respond(response);
        }
    }
}
