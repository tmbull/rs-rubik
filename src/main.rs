use kiss3d::camera::{ArcBall, Camera};
use kiss3d::light::Light;
use kiss3d::nalgebra::{Point3, Translation2, Translation3, UnitComplex, UnitQuaternion, Vector3};
use kiss3d::resource::Mesh;
use kiss3d::scene::SceneNode;
use kiss3d::window::Window;
use std::borrow::BorrowMut;
use std::cell::RefCell;
use std::rc::Rc;

#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(u8)]
enum Color {
    White = 1,
    Yellow,
    Blue,
    Green,
    Red,
    Orange,
}

/// To keep things simple, we are considering only face turns as moves, and are using the
/// "quarter turn metric".
#[derive(PartialEq, Clone, Copy)]
enum Rotation {
    Clockwise,
    Counterclockwise,
}

/// A 'sticker' for lack of a better term is one square on the face
/// The grid is organized such that 0, 0 = bottom, left.
///
/// We use a 2-d array so that we can support different sizes of cubes later.
#[derive(Clone, Debug, PartialEq)]
struct Face {
    stickers: [[Color; 3]; 3], /*
                                  top_left: Color,
                                  top_center: Color,
                                  top_right: Color,
                                  center_left: Color,
                                  center_center: Color,
                                  center_right: Color,
                                  bottom_left: Color,
                                  bottom_center: Color,
                                  bottom_right: Color
                               */
}

impl Face {
    pub fn new(stickers: [[Color; 3]; 3]) -> Self {
        Self { stickers }
    }

    pub fn rotate(&mut self, rotation: Rotation) {
        let n: usize = self.stickers.len();
        let x: usize = n / 2;
        let y: usize = n - 1;
        for i in 0..x {
            for j in i..y - i {
                match rotation {
                    Rotation::Clockwise => {
                        let k = self.stickers[i][j];
                        self.stickers[i][j] = self.stickers[y - j][i];
                        self.stickers[y - j][i] = self.stickers[y - i][y - j];
                        self.stickers[y - i][y - j] = self.stickers[j][y - i];
                        self.stickers[j][y - i] = k;
                    }
                    Rotation::Counterclockwise => {
                        let k = self.stickers[i][j];
                        self.stickers[i][j] = self.stickers[j][y - i];
                        self.stickers[j][y - i] = self.stickers[y - i][y - j];
                        self.stickers[y - i][y - j] = self.stickers[y - j][i];
                        self.stickers[y - j][i] = k;
                    }
                }
            }
        }
    }
}

#[derive(PartialEq, Clone, Copy)]
enum CubeFace {
    Top = 0,
    Left = 1,
    Front = 2,
    Bottom = 3,
    Right = 4,
    Back = 5,
}

/// The cube has 6 faces. They are organized in the array such that opposite faces are separated
/// by 3. In pictures:
///    T         0
///  L F R Ba  1 2 4 5    where T = top, L = left, F = front, Ba = back, Bo = bottom
///    Bo        3
struct Cube {
    faces: [Face; 6], // top: Face,
                      // left: Face,
                      // front: Face,
                      // bottom: Face,
                      // right: Face
                      // back: Face,
}

impl Cube {
    pub fn rotate_face(mut self, face: CubeFace, rotation: Rotation) -> Self {
        let idx = face as usize;
        self.faces[idx].rotate(rotation);
        self
    }
}
const NUM_SIDES: usize = 6;

/// This is used to "look up" the surrounding faces when rotating a face.
/// There might be a way to calculate this (linear algebra?), but it's hard-coded for now.
const SURROUNDING_FACES: [[usize; 4]; NUM_SIDES] = [
    [1, 2, 4, 5],
    [0, 2, 3, 5],
    [0, 4, 3, 1],
    [5, 4, 2, 1],
    [1, 3, 4, 0],
    [1, 3, 4, 0],
];

/// This is used to "look up" the opposite face on the cube.
/// TODO: Delete if unused.
const OPPOSITE_FACES: [usize; NUM_SIDES] = [3, 4, 5, 0, 1, 2];

/// We draw each face of the cube in the x-y plane and then move and rotate it into position. This
/// array stores the distance to move (translate) each face. The cube is 6x6 units, so each side is
/// 3 units away from the origin in one direction.
const TRANSLATIONS: [Translation3<f32>; NUM_SIDES] = [
    Translation3::new(0.0, 3.0, 0.0),
    Translation3::new(-3.0, 0.0, 0.0),
    Translation3::new(0.0, 0.0, 3.0),
    Translation3::new(0.0, -3.0, 0.0),
    Translation3::new(3.0, 0.0, 0.0),
    Translation3::new(0.0, 0.0, -3.0),
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
        direction: Vector3::new(0.0, -1.0, 0.0),
        up: Vector3::new(0.0, 0.0, 1.0),
    },
    FaceDirection {
        direction: Vector3::new(1.0, 0.0, 0.0),
        up: Vector3::new(0.0, 1.0, 0.0),
    },
    FaceDirection {
        direction: Vector3::new(0.0, 0.0, -1.0),
        up: Vector3::new(0.0, 1.0, 0.0),
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
    face: &Face,
    translation: &Translation3<f32>,
    rotation: &FaceDirection,
) {
    window.set_line_width(STICKER_BORDER);
    for (row_idx, row) in face.stickers.iter().enumerate() {
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
                col_idx as f32 * 2.0 - face.stickers.len() as f32 + 1.0,
                face.stickers[0].len() as f32 - 1.0 - (row_idx as f32 * 2.0),
                0.0,
            ));
            let q = UnitQuaternion::face_towards(&rotation.direction, &rotation.up);
            c.append_rotation(&q);
            c.append_translation(translation);
        }
    }
}

fn add_cube(window: &mut Window, cube: &Cube) {
    for (idx, face) in cube.faces.iter().enumerate() {
        add_face(window, &face, &TRANSLATIONS[idx], &ROTATIONS[idx]);
    }
}

fn main() {
    let mut window = Window::new("Kiss3d: custom_mesh");
    let mut origin = window.add_sphere(0.1);
    origin.set_color(0.0, 1.0, 0.0);

    let mut cube = Cube {
        faces: [
            Face::new([
                [Color::White, Color::White, Color::White],
                [Color::White, Color::White, Color::White],
                [Color::White, Color::White, Color::White],
            ]),
            Face::new([
                [Color::Red, Color::Red, Color::Red],
                [Color::Red, Color::Red, Color::Green],
                [Color::Yellow, Color::Orange, Color::Green],
            ]),
            Face::new([
                [Color::Blue, Color::Blue, Color::Blue],
                [Color::Red, Color::Blue, Color::Yellow],
                [Color::Red, Color::Blue, Color::Orange],
            ]),
            Face::new([
                [Color::Yellow, Color::Yellow, Color::Blue],
                [Color::Blue, Color::Yellow, Color::Green],
                [Color::Red, Color::Green, Color::Orange],
            ]),
            Face::new([
                [Color::Orange, Color::Orange, Color::Orange],
                [Color::Orange, Color::Orange, Color::Red],
                [Color::Yellow, Color::Orange, Color::Green],
            ]),
            Face::new([
                [Color::Green, Color::Green, Color::Green],
                [Color::Blue, Color::Green, Color::Yellow],
                [Color::Yellow, Color::Yellow, Color::Blue],
            ]),
        ],
    };
    cube = cube.rotate_face(CubeFace::Front, Rotation::Clockwise);
    add_cube(&mut window, &cube);
    window.set_light(Light::StickToCamera);
    let mut cam = ArcBall::new(Point3::new(0.0f32, 0.0, 15.0), Point3::origin());
    // let rot = UnitQuaternion::from_axis_angle(&Vector3::y_axis(), 0.014);

    while window.render_with_camera(&mut cam) {
        // println!("Cam location: {:?}", cam.eye());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Rotation::{Clockwise, Counterclockwise};
    use arr_macro::arr;
    use quickcheck::{quickcheck, Arbitrary, Gen};
    // use rstest::rstest;

    impl Arbitrary for Color {
        fn arbitrary(g: &mut Gen) -> Self {
            let val = (u8::arbitrary(g) % 6) + 1;
            unsafe { std::mem::transmute(val as u8) }
        }
    }

    impl Arbitrary for Face {
        fn arbitrary(g: &mut Gen) -> Face {
            Face {
                stickers: {
                    let mut mk_inner = || arr![Color::arbitrary(g); 3];
                    arr![mk_inner(); 3]
                },
            }
        }
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_face_clockwise_90_degrees(face: Face) -> bool {
        let expected = face.clone();
        let mut result = face;
        result.rotate(Clockwise);
        println!("{:?}", &result);

        // Verify that iterating the rotated face vertically from top right == iterating the
        // original face horizontally from top left
        for (row_idx, row) in expected.stickers.iter().enumerate() {
            for (col_idx, color) in row.iter().enumerate() {
                if *color != result.stickers[col_idx][result.stickers.len() - row_idx - 1] {
                    return false;
                }
            }
        }

        return true;
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_face_counterclockwise_90_degrees(face: Face) -> bool {
        let expected = face.clone();
        let mut result = face;
        result.rotate(Counterclockwise);

        // Verify that iterating the rotated face vertically from bottom left == iterating the
        // original face horizontally from top left
        for (row_idx, row) in expected.stickers.iter().enumerate() {
            for (col_idx, color) in row.iter().enumerate() {
                if *color != result.stickers[result.stickers[0].len() - col_idx - 1][row_idx] {
                    return false;
                }
            }
        }

        return true;
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_face_clockwise_180_degrees(face: Face) -> bool {
        let expected = face.clone();
        let mut result = face;
        result.rotate(Clockwise);
        result.rotate(Clockwise);

        // Verify that iterating the rotated face from bottom right == iterating the original face
        // from top left
        for (row_idx, row) in expected.stickers.iter().enumerate() {
            for (col_idx, color) in row.iter().enumerate() {
                if *color
                    != result.stickers[result.stickers.len() - row_idx - 1]
                        [result.stickers[0].len() - col_idx - 1]
                {
                    return false;
                }
            }
        }

        return true;
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_face_180_degrees_clockwise_equals_counterclockwise(face: Face) -> bool {
        let mut clockwise = face.clone();
        clockwise.rotate(Clockwise);
        clockwise.rotate(Clockwise);

        let mut counterclockwise = face;
        counterclockwise.rotate(Counterclockwise);
        counterclockwise.rotate(Counterclockwise);

        clockwise == counterclockwise
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_face_clockwise_four_times_is_identity(face: Face) -> bool {
        let expected = face.clone();
        let mut result = face;
        result.rotate(Clockwise);
        result.rotate(Clockwise);
        result.rotate(Clockwise);
        result.rotate(Clockwise);
        expected == result
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_face_counterclockwise_four_times_is_identity(face: Face) -> bool {
        let expected = face.clone();
        let mut result = face;
        result.rotate(Counterclockwise);
        result.rotate(Counterclockwise);
        result.rotate(Counterclockwise);
        result.rotate(Counterclockwise);
        expected == result
    }
}
