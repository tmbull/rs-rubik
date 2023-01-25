extern crate core;

mod array_cube;
mod cube;

use crate::array_cube::ArrayCube;
use crate::cube::CubeFace::{Front, Left};
use crate::cube::Direction::{Clockwise, Counterclockwise};
use crate::cube::{Color, CUBE_SIZE, NUM_SIDES};
use kiss3d::camera::ArcBall;
use kiss3d::light::Light;
use kiss3d::nalgebra::{Point3, Translation3, UnitQuaternion, Vector3};
use kiss3d::resource::Mesh;
use kiss3d::scene::SceneNode;
use kiss3d::window::Window;
use std::cell::RefCell;
use std::rc::Rc;

/// We draw each face of the cube in the x-y plane and then move and rotate it into position. This
/// array stores the distance to move (translate) each face. The cube is 6x6 units, so each side is
/// 3 units away from the origin in one direction.
const TRANSLATIONS: [Translation3<f32>; NUM_SIDES] = [
    Translation3::new(0.0, 3.0, 0.0),
    Translation3::new(-3.0, 0.0, 0.0),
    Translation3::new(0.0, 0.0, 3.0),
    Translation3::new(3.0, 0.0, 0.0),
    Translation3::new(0.0, 0.0, -3.0),
    Translation3::new(0.0, -3.0, 0.0),
];

/// A simple struct for storing the normal direction and "up" direction of each face.
/// We draw each face in the x-y plane and then move, and rotate it into position. This struct
/// is used to store that rotation for each side of the cube.
struct FaceDirection {
    direction: Vector3<f32>,
    up: Vector3<f32>,
}

/// The direction for each side of the cube. The array is indexed in the same order as documented
/// in [Cube].
const ROTATIONS: [FaceDirection; NUM_SIDES] = [
    FaceDirection {
        direction: Vector3::new(0.0, 1.0, 0.0),
        up: Vector3::new(0.0, 0.0, -1.0),
    },
    FaceDirection {
        direction: Vector3::new(-1.0, 0.0, 0.0),
        up: Vector3::new(0.0, 1.0, 0.0),
    },
    FaceDirection {
        direction: Vector3::new(0.0, 0.0, 1.0),
        up: Vector3::new(0.0, 1.0, 0.0),
    },
    FaceDirection {
        direction: Vector3::new(1.0, 0.0, 0.0),
        up: Vector3::new(0.0, 1.0, 0.0),
    },
    FaceDirection {
        direction: Vector3::new(0.0, 0.0, -1.0),
        up: Vector3::new(0.0, 1.0, 0.0),
    },
    FaceDirection {
        direction: Vector3::new(0.0, -1.0, 0.0),
        up: Vector3::new(0.0, 0.0, 1.0),
    },
];

/// A helper function for looking up RGB color for each sticker [Color]. Since `set_color` takes RGB
/// params instead of a struct, we take in a [SceneNode] and set the color directly.
fn set_color(scene_node: &mut SceneNode, color: &Color) {
    match color {
        Color::White => scene_node.set_color(1.0, 1.0, 1.0),
        Color::Yellow => scene_node.set_color(1.0, 1.0, 0.0),
        Color::Blue => scene_node.set_color(0.0, 0.0, 1.0),
        Color::Green => scene_node.set_color(0.0, 1.0, 0.0),
        Color::Red => scene_node.set_color(1.0, 0.0, 0.0),
        Color::Orange => scene_node.set_color(1.0, 0.65, 0.0),
    }
}

const CUBELET_SIZE: f32 = 1.0;
const STICKER_BORDER: f32 = 0.1;
const STICKER_SIZE: f32 = CUBELET_SIZE - STICKER_BORDER;

/// Add one face to the scene. This is called for each face by [add_cube].
fn add_face(
    window: &mut Window,
    face: &[[Color; CUBE_SIZE]; CUBE_SIZE],
    translation: &Translation3<f32>,
    rotation: &FaceDirection,
) {
    for (row_idx, row) in face.iter().enumerate() {
        for (col_idx, color) in row.iter().enumerate() {
            window.draw_line(
                &Point3::new(-CUBELET_SIZE, -CUBELET_SIZE, 0.0),
                &Point3::new(CUBELET_SIZE, -CUBELET_SIZE, 0.0),
                &Point3::new(0.0, 0.0, 0.0),
            );
            let a = Point3::new(-STICKER_SIZE, -STICKER_SIZE, 0.0);
            let b = Point3::new(STICKER_SIZE, -STICKER_SIZE, 0.0);
            let c = Point3::new(STICKER_SIZE, STICKER_SIZE, 0.0);
            let d = Point3::new(-STICKER_SIZE, STICKER_SIZE, 0.0);
            let vertices = vec![a, b, c, d];
            let indices = vec![Point3::new(0, 1, 2), Point3::new(2, 3, 0)];

            let m = Mesh::new(vertices, indices, None, None, false);

            let mesh = Rc::new(RefCell::new(m));
            let mut c = window.add_mesh(mesh, Vector3::new(1.0, 1.0, 1.0));
            set_color(&mut c, color);
            c.enable_backface_culling(true);
            c.append_translation(&Translation3::new(
                col_idx as f32 * 2.0 - CUBE_SIZE as f32 + 1.0,
                CUBE_SIZE as f32 - 1.0 - (row_idx as f32 * 2.0),
                0.0,
            ));
            let q = UnitQuaternion::face_towards(&rotation.direction, &rotation.up);
            c.append_rotation(&q);
            c.append_translation(translation);
        }
    }
}

fn add_cube(window: &mut Window, cube: &ArrayCube) {
    for (idx, face) in cube.facelets.iter().enumerate() {
        add_face(window, &face, &TRANSLATIONS[idx], &ROTATIONS[idx]);
    }
}

fn main() {
    let mut window = Window::new("Kiss3d: custom_mesh");
    // let mut origin = window.add_sphere(0.1);
    // origin.set_color(0.0, 1.0, 0.0);

    let mut cube = ArrayCube::new([
        [
            [Color::White, Color::White, Color::White],
            [Color::White, Color::White, Color::White],
            [Color::White, Color::White, Color::White],
        ],
        [
            [Color::Green, Color::Green, Color::Green],
            [Color::Blue, Color::Green, Color::Yellow],
            [Color::Yellow, Color::Yellow, Color::Blue],
        ],
        [
            [Color::Red, Color::Red, Color::Red],
            [Color::Red, Color::Red, Color::Green],
            [Color::Yellow, Color::Orange, Color::Green],
        ],
        [
            [Color::Blue, Color::Blue, Color::Blue],
            [Color::Red, Color::Blue, Color::Yellow],
            [Color::Red, Color::Blue, Color::Orange],
        ],
        [
            [Color::Orange, Color::Orange, Color::Orange],
            [Color::Orange, Color::Orange, Color::Red],
            [Color::Yellow, Color::Orange, Color::Green],
        ],
        [
            [Color::Red, Color::Blue, Color::Yellow],
            [Color::Green, Color::Yellow, Color::Yellow],
            [Color::Orange, Color::Green, Color::Blue],
        ],
    ]);
    // let mut cube = Cube::default();
    cube.rotate_face(Front, Clockwise);
    cube.rotate_face(Left, Counterclockwise);
    // cube.randomize(50, 100);
    // println!("{:?}", cube.solve_recursive());
    add_cube(&mut window, &cube);
    window.set_light(Light::StickToCamera);
    let mut cam = ArcBall::new(Point3::new(10.0, 7.0, 10.0), Point3::origin());
    // let rot = UnitQuaternion::from_axis_angle(&Vector3::y_axis(), 0.014);

    while window.render_with_camera(&mut cam) {
        // println!("Cam location: {:?}", cam.eye());
    }
}
