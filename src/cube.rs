//! This module contains the software representation of the Rubik's cube, and code related to
//! manipulating and solving the cube.
use crate::cube::Color::{Blue, Green, Orange, Red, White, Yellow};
use std::borrow::BorrowMut;
use std::slice::Iter;

/// The number of sides of the puzzle. In the case of a cube, the number of sides is 6. It is
/// unlikely that this solution can be generalized for other shapes, but it would be an interesting
/// exercise some day.
pub const NUM_SIDES: usize = 6;

/// This is used to "look up" the surrounding faces when rotating a face.
/// There might be a way to calculate this (linear algebra?), but it's hard-coded for now.
/// These are indexed in the following order when looking at the face head on:
/// (top, right, bottom, left)
/// The specific order is not important except that they are ordered clockwise. We will iterate
/// through this array in order when rotating clockwise, and in reverse order when rotating
/// counterclockwise.
const SURROUNDING_FACES: [[usize; NUM_SIDES - 2]; NUM_SIDES] = [
    [5, 4, 2, 1],
    [0, 2, 3, 5],
    [0, 4, 3, 1],
    [2, 4, 5, 1],
    [0, 5, 3, 2],
    [0, 1, 3, 4],
];

/// This is used to "look up" the surrounding indexes of a face. This is used in conjunction
/// with [SURROUNDING_FACES] to rotate the facelets that surround a face. It is important that the
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

#[derive(Clone, Copy, Debug, PartialEq)]
#[repr(u8)]
pub enum Color {
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
pub enum Rotation {
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

/// A 'facelet' for lack of a better term is one square on the face
/// The grid is organized such that 0, 0 = bottom, left.
///
/// We use a 2-d array so that we can support different sizes of cubes later.
#[derive(Clone, Debug, PartialEq)]
pub struct Face {
    facelets: [[Color; 3]; 3], /*
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
    pub fn new(facelets: [[Color; 3]; 3]) -> Self {
        Self { facelets }
    }

    pub fn new_solid_color(color: Color) -> Self {
        Self {
            facelets: [[color; 3]; 3],
        }
    }

    pub fn rotate(&mut self, rotation: Rotation) {
        let n: usize = self.facelets.len();
        let x: usize = n / 2;
        let y: usize = n - 1;
        for i in 0..x {
            for j in i..y - i {
                match rotation {
                    Rotation::Clockwise => {
                        let k = self.facelets[i][j];
                        self.facelets[i][j] = self.facelets[y - j][i];
                        self.facelets[y - j][i] = self.facelets[y - i][y - j];
                        self.facelets[y - i][y - j] = self.facelets[j][y - i];
                        self.facelets[j][y - i] = k;
                    }
                    Rotation::Counterclockwise => {
                        let k = self.facelets[i][j];
                        self.facelets[i][j] = self.facelets[j][y - i];
                        self.facelets[j][y - i] = self.facelets[y - i][y - j];
                        self.facelets[y - i][y - j] = self.facelets[y - j][i];
                        self.facelets[y - j][i] = k;
                    }
                }
            }
        }
    }

    pub fn get_row(&self, row_idx: usize) -> [Color; 3] {
        self.facelets[row_idx]
    }

    pub fn get_col(&self, col_idx: usize) -> [Color; 3] {
        [
            self.facelets[0][col_idx],
            self.facelets[1][col_idx],
            self.facelets[2][col_idx],
        ]
    }

    pub fn get_row_iter(&self) -> Iter<'_, [Color; 3]> {
        self.facelets.iter()
    }

    /// The 'size' of the cube is defined as the number of facelets in a row. In other words, a
    /// standard 3x3 Rubik's cube has a size of '3'. This is defined on both cube and face for
    /// convenience.
    pub fn get_size(&self) -> usize {
        self.facelets.len()
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum CubeFace {
    Top = 0,
    Left = 1,
    Front = 2,
    Bottom = 3,
    Right = 4,
    Back = 5,
}

/// The cube has NUM_SIDES faces. They are organized in the array such that opposite faces are separated
/// by 3. In pictures:
///    T         0
///  L F R Ba  1 2 4 5    where T = top, L = left, F = front, Ba = back, Bo = bottom
///    Bo        3
#[derive(Debug, PartialEq, Clone)]
pub struct Cube {
    faces: [Face; NUM_SIDES], // top: Face,
                              // left: Face,
                              // front: Face,
                              // bottom: Face,
                              // right: Face
                              // back: Face,
}

impl Default for Cube {
    fn default() -> Self {
        Cube::new([
            Face::new_solid_color(White),
            Face::new_solid_color(Green),
            Face::new_solid_color(Red),
            Face::new_solid_color(Yellow),
            Face::new_solid_color(Blue),
            Face::new_solid_color(Orange),
        ])
    }
}

impl Cube {
    pub fn new(faces: [Face; NUM_SIDES]) -> Self {
        Self { faces }
    }

    pub fn get_face_iter(&self) -> Iter<'_, Face> {
        self.faces.iter()
    }

    pub fn rotate_face(&mut self, face: CubeFace, rotation: Rotation) {
        let idx = face as usize;
        self.faces[idx].rotate(rotation);
        let [top, right, bottom, left] = SURROUNDING_FACES[idx];
        match rotation {
            Rotation::Clockwise => {
                for i in 0..=2 {
                    let tmp = self.swap_facelets(
                        top,
                        SURROUNDING_INDICES[0][i],
                        right,
                        SURROUNDING_INDICES[1][i],
                    );
                    self.swap_facelets(
                        left,
                        SURROUNDING_INDICES[3][i],
                        top,
                        SURROUNDING_INDICES[0][i],
                    );
                    self.swap_facelets(
                        bottom,
                        SURROUNDING_INDICES[2][i],
                        left,
                        SURROUNDING_INDICES[3][i],
                    );
                    self.set_facelet(bottom, SURROUNDING_INDICES[2][i], tmp);
                }
            }
            Rotation::Counterclockwise => {
                for i in 0..=2 {
                    let tmp = self.swap_facelets(
                        top,
                        SURROUNDING_INDICES[0][i],
                        left,
                        SURROUNDING_INDICES[3][i],
                    );
                    self.swap_facelets(
                        right,
                        SURROUNDING_INDICES[1][i],
                        top,
                        SURROUNDING_INDICES[0][i],
                    );
                    self.swap_facelets(
                        bottom,
                        SURROUNDING_INDICES[2][i],
                        right,
                        SURROUNDING_INDICES[1][i],
                    );
                    self.set_facelet(bottom, SURROUNDING_INDICES[2][i], tmp);
                }
            }
        }
    }

    fn swap_facelets(
        &mut self,
        src_face: usize,
        src_facelet: (usize, usize),
        dest_face: usize,
        dest_facelet: (usize, usize),
    ) -> Color {
        let src_color = self.faces[src_face].facelets[src_facelet.0][src_facelet.1];
        let dest_face_ref = self.faces[dest_face].borrow_mut();
        let dest_color = dest_face_ref.facelets[dest_facelet.0][dest_facelet.1];
        dest_face_ref.facelets[dest_facelet.0][dest_facelet.1] = src_color;
        dest_color
    }

    fn set_facelet(&mut self, dest_face: usize, dest_facelet: (usize, usize), color: Color) {
        let dest_face_ref = self.faces[dest_face].borrow_mut();
        dest_face_ref.facelets[dest_facelet.0][dest_facelet.1] = color;
    }

    /// The 'size' of the cube is defined as the number of facelets in a row. In other words, a
    /// standard 3x3 Rubik's cube has a size of '3'.
    pub fn get_size(&self) -> usize {
        self.faces[0].get_size()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cube;
    use crate::cube::Rotation::{Clockwise, Counterclockwise};
    use crate::CubeFace::{Back, Bottom, Front, Left, Right, Top};
    use arr_macro::arr;
    use quickcheck::{Arbitrary, Gen};
    use std::mem;

    const NUM_SIDES: u8 = cube::NUM_SIDES as u8;

    impl Arbitrary for Rotation {
        fn arbitrary(g: &mut Gen) -> Self {
            let val = u8::arbitrary(g) % 2;
            unsafe { mem::transmute(val as u8) }
        }
    }

    impl Arbitrary for CubeFace {
        fn arbitrary(g: &mut Gen) -> Self {
            let val = u8::arbitrary(g) % NUM_SIDES;
            unsafe { mem::transmute(val as u8) }
        }
    }

    impl Arbitrary for Color {
        fn arbitrary(g: &mut Gen) -> Self {
            let val = (u8::arbitrary(g) % NUM_SIDES) + 1;
            unsafe { mem::transmute(val as u8) }
        }
    }

    impl Arbitrary for Face {
        fn arbitrary(g: &mut Gen) -> Face {
            Face {
                facelets: {
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

    fn face_and_surround_match(
        result: &Cube,
        expected_front: &Face,
        expected_right: &[Color; 3],
        expected_bottom: &[Color; 3],
        expected_left: &[Color; 3],
        expected_top: &[Color; 3],
    ) -> bool {
        &result.faces[Front as usize] == expected_front
            && expected_right == &result.faces[Right as usize].get_col(0)
            && expected_bottom == &result.faces[Bottom as usize].get_row(0)
            && expected_left == &result.faces[Left as usize].get_col(2)
            && expected_top == &result.faces[Top as usize].get_row(2)
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

        face_and_surround_match(
            &result,
            &expected_front,
            &expected_right,
            &expected_bottom,
            &expected_left,
            &expected_top,
        )
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

        face_and_surround_match(
            &result,
            &expected_front,
            &expected_right,
            &expected_bottom,
            &expected_left,
            &expected_top,
        )
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
    fn swap_cube_facelets(
        cube: Cube,
        src_face: u8,
        (src_x, src_y): (u8, u8),
        dest_face: u8,
        (dest_x, dest_y): (u8, u8),
    ) -> bool {
        let src_face = (src_face % NUM_SIDES as u8) as usize;
        let src_facelet = ((src_x % 3) as usize, (src_y % 3) as usize);
        let dest_facelet = ((dest_x % 3) as usize, (dest_y % 3) as usize);
        let dest_face = (dest_face % NUM_SIDES as u8) as usize;
        let expected = cube.clone();
        let mut result = cube;
        let old_dest = result.swap_facelets(src_face, src_facelet, dest_face, dest_facelet);

        old_dest == expected.faces[dest_face].facelets[dest_facelet.0][dest_facelet.1]
            && result.faces[dest_face].facelets[dest_facelet.0][dest_facelet.1]
                == expected.faces[src_face].facelets[src_facelet.0][src_facelet.1]
    }

    #[quickcheck_macros::quickcheck]
    fn set_cube_facelet(
        cube: Cube,
        dest_face: u8,
        (dest_x, dest_y): (u8, u8),
        color: Color,
    ) -> bool {
        let dest_facelet = ((dest_x % 3) as usize, (dest_y % 3) as usize);
        let dest_face = (dest_face % NUM_SIDES as u8) as usize;
        let mut result = cube;
        result.set_facelet(dest_face, dest_facelet, color);

        result.faces[dest_face].facelets[dest_facelet.0][dest_facelet.1] == color
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_face_clockwise_90_degrees(face: Face) -> bool {
        let expected = face.clone();
        let mut result = face;
        result.rotate(Clockwise);

        // Verify that iterating the rotated face vertically from top right == iterating the
        // original face horizontally from top left
        for (row_idx, row) in expected.facelets.iter().enumerate() {
            for (col_idx, color) in row.iter().enumerate() {
                if *color != result.facelets[col_idx][result.facelets.len() - row_idx - 1] {
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
        for (row_idx, row) in expected.facelets.iter().enumerate() {
            for (col_idx, color) in row.iter().enumerate() {
                if *color != result.facelets[result.facelets[0].len() - col_idx - 1][row_idx] {
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
        for (row_idx, row) in expected.facelets.iter().enumerate() {
            for (col_idx, color) in row.iter().enumerate() {
                if *color
                    != result.facelets[result.facelets.len() - row_idx - 1]
                        [result.facelets[0].len() - col_idx - 1]
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
