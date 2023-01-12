use kiss3d::camera::ArcBall;
use kiss3d::light::Light;
use kiss3d::nalgebra::{Point3, Translation3, UnitQuaternion, Vector3};
use kiss3d::resource::Mesh;
use kiss3d::scene::SceneNode;
use kiss3d::window::Window;
use std::borrow::BorrowMut;
use std::cell::RefCell;
use std::mem;
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
#[derive(Debug, PartialEq, Clone, Copy)]
enum Rotation {
    Clockwise = 0,
    Counterclockwise,
}

impl Rotation {
    fn get_opposite(self) -> Self {
        match self {
            Rotation::Clockwise => Rotation::Counterclockwise,
            Rotation::Counterclockwise => Rotation::Clockwise,
        }
    }
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

    pub fn get_row(&self, row_idx: usize) -> [Color; 3] {
        self.stickers[row_idx]
    }

    pub fn get_col(&self, col_idx: usize) -> [Color; 3] {
        [
            self.stickers[0][col_idx],
            self.stickers[1][col_idx],
            self.stickers[2][col_idx],
        ]
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
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
#[derive(Debug, PartialEq, Clone)]
struct Cube {
    faces: [Face; 6], // top: Face,
                      // left: Face,
                      // front: Face,
                      // bottom: Face,
                      // right: Face
                      // back: Face,
}

impl Cube {
    pub fn rotate_face(&mut self, face: CubeFace, rotation: Rotation) {
        let idx = face as usize;
        self.faces[idx].rotate(rotation);
        let [top, right, bottom, left] = SURROUNDING_FACES[idx];
        match rotation {
            Rotation::Clockwise => {
                for i in 0..=2 {
                    let tmp = self.swap_stickers(
                        top,
                        SURROUNDING_INDICES[0][i],
                        right,
                        SURROUNDING_INDICES[1][i],
                    );
                    self.swap_stickers(
                        left,
                        SURROUNDING_INDICES[3][i],
                        top,
                        SURROUNDING_INDICES[0][i],
                    );
                    self.swap_stickers(
                        bottom,
                        SURROUNDING_INDICES[2][i],
                        left,
                        SURROUNDING_INDICES[3][i],
                    );
                    self.set_sticker(bottom, SURROUNDING_INDICES[2][i], tmp);
                }
            }
            Rotation::Counterclockwise => {
                for i in 0..=2 {
                    let tmp = self.swap_stickers(
                        top,
                        SURROUNDING_INDICES[0][i],
                        left,
                        SURROUNDING_INDICES[3][i],
                    );
                    self.swap_stickers(
                        right,
                        SURROUNDING_INDICES[1][i],
                        top,
                        SURROUNDING_INDICES[0][i],
                    );
                    self.swap_stickers(
                        bottom,
                        SURROUNDING_INDICES[2][i],
                        right,
                        SURROUNDING_INDICES[1][i],
                    );
                    self.set_sticker(bottom, SURROUNDING_INDICES[2][i], tmp);
                }
            }
        }
    }

    fn swap_stickers(
        &mut self,
        src_face: usize,
        src_sticker: (usize, usize),
        dest_face: usize,
        dest_sticker: (usize, usize),
    ) -> Color {
        let src_color = self.faces[src_face].stickers[src_sticker.0][src_sticker.1];
        let dest_face_ref = self.faces[dest_face].borrow_mut();
        let dest_color = dest_face_ref.stickers[dest_sticker.0][dest_sticker.1];
        dest_face_ref.stickers[dest_sticker.0][dest_sticker.1] = src_color;
        dest_color
    }

    fn set_sticker(&mut self, dest_face: usize, dest_sticker: (usize, usize), color: Color) {
        let dest_face_ref = self.faces[dest_face].borrow_mut();
        dest_face_ref.stickers[dest_sticker.0][dest_sticker.1] = color;
    }
}
const NUM_SIDES: usize = 6;

/// This is used to "look up" the surrounding faces when rotating a face.
/// There might be a way to calculate this (linear algebra?), but it's hard-coded for now.
/// These are indexed in the following order when looking at the face head on:
/// (top, right, bottom, left)
/// The specific order is not important except that they are ordered clockwise. We will iterate
/// through this array in order when rotating clockwise, and in reverse order when rotating
/// counterclockwise.
const SURROUNDING_FACES: [[usize; 4]; NUM_SIDES] = [
    [5, 4, 2, 1],
    [0, 2, 3, 5],
    [0, 4, 3, 1],
    [2, 4, 5, 1],
    [0, 5, 3, 2],
    [0, 1, 3, 4],
];

/// This is used to "look up" the surrounding indexes of a face. This is used in conjunction
/// with [SURROUNDING_FACES] to rotate the stickers that surround a face. It is important that the
/// order if this array is the same as the order of the above array (i.e. top, right, bottom, left).
/// This is because when rotating a face clockwise, the bottom row of the face "above" will become
/// the left-most column of the face "to the right", the left-most column of the face "to the right"
/// will become the top row on the face "below", the top row on the face "below" will become the
/// right-most column on the face "to the left", and the right-most column on the face "to the left"
/// will become the bottom row on the face "above". The reverse is true when rotating
/// counterclockwise.
const SURROUNDING_INDICES: [[(usize, usize); 3]; 4] = [
    [(2, 0), (2, 1), (2, 2)],
    [(0, 0), (1, 0), (2, 0)],
    [(0, 2), (0, 1), (0, 0)],
    [(2, 2), (1, 2), (0, 2)],
];

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
    // let mut origin = window.add_sphere(0.1);
    // origin.set_color(0.0, 1.0, 0.0);

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
    cube.rotate_face(CubeFace::Front, Rotation::Counterclockwise);
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
    use crate::CubeFace::{Back, Bottom, Front, Left, Right, Top};
    use crate::Rotation::{Clockwise, Counterclockwise};
    use arr_macro::arr;
    use quickcheck::{Arbitrary, Gen};
    // use rstest::rstest;

    impl Arbitrary for Rotation {
        fn arbitrary(g: &mut Gen) -> Self {
            let val = u8::arbitrary(g) % 2;
            unsafe { mem::transmute(val as u8) }
        }
    }

    impl Arbitrary for CubeFace {
        fn arbitrary(g: &mut Gen) -> Self {
            let val = u8::arbitrary(g) % 6;
            unsafe { mem::transmute(val as u8) }
        }
    }

    impl Arbitrary for Color {
        fn arbitrary(g: &mut Gen) -> Self {
            let val = (u8::arbitrary(g) % 6) + 1;
            unsafe { mem::transmute(val as u8) }
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

    impl Arbitrary for Cube {
        fn arbitrary(g: &mut Gen) -> Cube {
            Cube {
                faces: { arr![Face::arbitrary(g); 6] },
            }
        }
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_cube_face_four_times_is_identity(cube: Cube, face: CubeFace) -> bool {
        let expected = cube.clone();
        let mut result = cube;
        result.rotate_face(face, Clockwise);
        result.rotate_face(face, Clockwise);
        result.rotate_face(face, Clockwise);
        result.rotate_face(face, Clockwise);

        result == expected
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_clockwise_then_counterclockwise_is_identity(cube: Cube, face: CubeFace) -> bool {
        let expected = cube.clone();
        let mut result = cube;
        result.rotate_face(face, Clockwise);
        result.rotate_face(face, Counterclockwise);

        result == expected
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_cube_face_180_degrees_clockwise_equals_counterclockwise(
        cube: Cube,
        face: CubeFace,
    ) -> bool {
        let mut expected = cube.clone();
        expected.rotate_face(face, Counterclockwise);
        expected.rotate_face(face, Counterclockwise);
        let mut result = cube;
        result.rotate_face(face, Clockwise);
        result.rotate_face(face, Clockwise);

        result == expected
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_front_face_counterclockwise(cube: Cube) -> bool {
        let expected = cube.clone();
        let mut result = cube;
        result.rotate_face(Front, Counterclockwise);
        let mut expected_front = expected.faces[Front as usize].clone();
        expected_front.rotate(Counterclockwise);
        let mut expected_right = expected.faces[Bottom as usize].get_row(0);
        expected_right.reverse();
        let expected_bottom = expected.faces[Left as usize].get_col(2);
        let mut expected_left = expected.faces[Top as usize].get_row(2);
        expected_left.reverse();
        let expected_top = expected.faces[Right as usize].get_col(0);

        result.faces[Front as usize] == expected_front
            && expected_right == result.faces[Right as usize].get_col(0)
            && expected_bottom == result.faces[Bottom as usize].get_row(0)
            && expected_left == result.faces[Left as usize].get_col(2)
            && expected_top == result.faces[Top as usize].get_row(2)
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_front_face_clockwise(cube: Cube) -> bool {
        let expected = cube.clone();
        let mut result = cube;
        result.rotate_face(Front, Clockwise);
        let mut expected_front = expected.faces[Front as usize].clone();
        expected_front.rotate(Clockwise);
        let expected_right = expected.faces[Top as usize].get_row(2);
        let mut expected_bottom = expected.faces[Right as usize].get_col(0);
        expected_bottom.reverse();
        let expected_left = expected.faces[Bottom as usize].get_row(0);
        let mut expected_top = expected.faces[Left as usize].get_col(2);
        expected_top.reverse();

        result.faces[Front as usize] == expected_front
            && expected_right == result.faces[Right as usize].get_col(0)
            && expected_bottom == result.faces[Bottom as usize].get_row(0)
            && expected_left == result.faces[Left as usize].get_col(2)
            && expected_top == result.faces[Top as usize].get_row(2)
    }

    #[quickcheck_macros::quickcheck]
    fn cube_rotations_can_be_undone(cube: Cube, mut rotations: Vec<(CubeFace, Rotation)>) -> bool {
        let expected = cube.clone();
        let mut result = cube;
        for (face, rotation) in rotations.iter() {
            result.rotate_face(*face, *rotation);
        }

        rotations.reverse();
        for (face, rotation) in rotations {
            result.rotate_face(face, rotation.get_opposite());
        }

        result == expected
    }

    #[quickcheck_macros::quickcheck]
    fn swap_cube_stickers(
        cube: Cube,
        src_face: u8,
        (src_x, src_y): (u8, u8),
        dest_face: u8,
        (dest_x, dest_y): (u8, u8),
    ) -> bool {
        let src_face = (src_face % NUM_SIDES as u8) as usize;
        let src_sticker = ((src_x % 3) as usize, (src_y % 3) as usize);
        let dest_sticker = ((dest_x % 3) as usize, (dest_y % 3) as usize);
        let dest_face = (dest_face % NUM_SIDES as u8) as usize;
        let expected = cube.clone();
        let mut result = cube;
        let old_dest = result.swap_stickers(src_face, src_sticker, dest_face, dest_sticker);

        old_dest == expected.faces[dest_face].stickers[dest_sticker.0][dest_sticker.1]
            && result.faces[dest_face].stickers[dest_sticker.0][dest_sticker.1]
                == expected.faces[src_face].stickers[src_sticker.0][src_sticker.1]
    }

    #[quickcheck_macros::quickcheck]
    fn set_cube_sticker(
        cube: Cube,
        dest_face: u8,
        (dest_x, dest_y): (u8, u8),
        color: Color,
    ) -> bool {
        let dest_sticker = ((dest_x % 3) as usize, (dest_y % 3) as usize);
        let dest_face = (dest_face % NUM_SIDES as u8) as usize;
        let mut result = cube;
        result.set_sticker(dest_face, dest_sticker, color);

        result.faces[dest_face].stickers[dest_sticker.0][dest_sticker.1] == color
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_face_clockwise_90_degrees(face: Face) -> bool {
        let expected = face.clone();
        let mut result = face;
        result.rotate(Clockwise);

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
