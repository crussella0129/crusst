use crusst::blend;
use crusst::builder::Shape;
use crusst::feature::ft;
use crusst::types::MeshSettings;

/// Filleted box meshes without panic and produces triangles.
#[test]
fn filleted_box_meshes() {
    let shape = Shape::box3(5.0, 5.0, 5.0).fillet(blend::g2(1.0), vec![ft(0, 0).all_edges()]);
    let settings = MeshSettings {
        max_depth: 5,
        ..Default::default()
    };
    let mesh = shape.mesh(settings);
    let tri_count = mesh.indices.len() / 3;
    let vert_count = mesh.vertices.len();
    assert!(
        tri_count > 10,
        "should produce reasonable triangle count: got {}",
        tri_count
    );
    assert!(
        vert_count > 10,
        "should produce reasonable vertex count: got {}",
        vert_count
    );
}

/// Round union meshes without panic.
#[test]
fn round_union_meshes() {
    let a = Shape::sphere(5.0);
    let b = Shape::sphere(5.0).translate(6.0, 0.0, 0.0);
    let shape = a.round_union(b, 1.0);
    let settings = MeshSettings {
        max_depth: 5,
        ..Default::default()
    };
    let mesh = shape.mesh(settings);
    let tri_count = mesh.indices.len() / 3;
    assert!(
        tri_count > 10,
        "round union should produce mesh: got {} tris",
        tri_count
    );
}

/// Chamfered box meshes without panic.
#[test]
fn chamfered_box_meshes() {
    let shape =
        Shape::box3(5.0, 5.0, 5.0).chamfer(blend::equal_chamfer(1.0), vec![ft(0, 0).all_edges()]);
    let settings = MeshSettings {
        max_depth: 5,
        ..Default::default()
    };
    let mesh = shape.mesh(settings);
    let tri_count = mesh.indices.len() / 3;
    assert!(
        tri_count > 10,
        "chamfered box should mesh: got {} tris",
        tri_count
    );
}

/// Rounded box (simple .round()) has reasonable vertex count.
#[test]
fn rounded_box_vertex_count() {
    let sharp = Shape::box3(5.0, 5.0, 5.0);
    let rounded = Shape::box3(5.0, 5.0, 5.0).round(1.0);
    let settings_sharp = MeshSettings {
        max_depth: 5,
        ..Default::default()
    };
    let settings_round = MeshSettings {
        max_depth: 5,
        ..Default::default()
    };
    let sharp_mesh = sharp.mesh(settings_sharp);
    let round_mesh = rounded.mesh(settings_round);
    // Rounded version should have comparable or more vertices (more surface detail at edges)
    assert!(
        round_mesh.vertices.len() >= sharp_mesh.vertices.len() / 2,
        "rounded should have comparable vertex count: sharp={}, rounded={}",
        sharp_mesh.vertices.len(),
        round_mesh.vertices.len()
    );
}
