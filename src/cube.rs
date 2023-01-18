//! This module contains the software representation of the Rubik's cube, and code related to
//! manipulating and solving the cube.
use crate::cube::Color::{Blue, Green, Orange, Red, White, Yellow};
use crate::cube::CubeFace::{Back, Down, Front, Left, Right, Up};
use crate::cube::Rotation::{Clockwise, Counterclockwise};
use enum_iterator::Sequence;
use itertools::Itertools;
use kiss3d::nalgebra::DimAdd;
use num_enum::TryFromPrimitive;
use rand::distributions::{DistIter, DistMap, Standard};
use rand::prelude::*;
use std::borrow::{Borrow, BorrowMut};
use std::cell::{Cell, RefCell};
use std::collections::{HashSet, VecDeque};
use std::convert::TryFrom;
use std::rc::Rc;
use std::slice::{Iter, IterMut};

/// The number of sides of the puzzle. In the case of a cube, the number of sides is 6. It is
/// unlikely that this solution can be generalized for other shapes, but it would be an interesting
/// exercise some day.
pub const NUM_SIDES: usize = 6;

/// In Rubik's cube terminology, 'size' is frequently used to refer to length of each side of the
/// cube in terms of number of facelets. For example, a standard 3x3x3 cube is said to be of size
/// '3'.
pub const CUBE_SIZE: usize = 3;

/// This is used to "look up" the surrounding faces when rotating a face.
/// There might be a way to calculate this (linear algebra?), but it's hard-coded for now.
/// These are indexed in the following order when looking at the face head on:
/// (up, right, down, left)
/// The specific order is not important except that they are ordered clockwise. We will iterate
/// through this array in order when rotating clockwise, and in reverse order when rotating
/// counterclockwise.
const SURROUNDING_FACES: [[usize; NUM_SIDES - 2]; NUM_SIDES] = [
    [5, 4, 2, 1],
    [0, 2, 3, 5], // Left
    [0, 4, 3, 1], // Front
    [2, 4, 5, 1],
    [0, 5, 3, 2],
    [0, 1, 3, 4],
];

fn copy_facelets(mut src: Iter<Color>, dest: IterMut<Color>) {
    for facelet in dest {
        *facelet = *(src.next().unwrap());
    }
}

/// This is used to "look up" the surrounding indexes of a face. This is used in conjunction
/// with [SURROUNDING_FACES] to rotate the facelets that surround a face. It is important that the
/// order if this array is the same as the order of the above array (i.e. up, right, down, left).
/// This is because when rotating a face clockwise, the down row of the face "above" will become
/// the left-most column of the face "to the right", the left-most column of the face "to the right"
/// will become the up row on the face "below", the up row on the face "below" will become the
/// right-most column on the face "to the left", and the right-most column on the face "to the left"
/// will become the down row on the face "above". The reverse is true when rotating
/// counterclockwise.
const SURROUNDING_INDICES: [[(usize, usize); CUBE_SIZE]; 4] = [
    [(2, 0), (2, 1), (2, 2)],
    [(0, 0), (1, 0), (2, 0)],
    [(0, 2), (0, 1), (0, 0)],
    [(2, 2), (1, 2), (0, 2)],
];

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Sequence)]
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
#[derive(Debug, PartialEq, Clone, Copy, TryFromPrimitive)]
#[repr(u8)]
pub enum Rotation {
    Clockwise = 0,
    Counterclockwise,
}

impl Rotation {
    pub fn get_opposite(self) -> Self {
        match self {
            Clockwise => Counterclockwise,
            Counterclockwise => Clockwise,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CubeMove {
    face: CubeFace,
    rotation: Rotation,
}

impl CubeMove {
    pub fn new(face: CubeFace, rotation: Rotation) -> Self {
        Self { face, rotation }
    }
}

/// A 'facelet' for lack of a better term is one square on the face
/// The grid is organized such that 0, 0 = down, left.
///
/// We use a 2-d array so that we can support different sizes of cubes later.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Face {
    facelets: [[Color; CUBE_SIZE]; CUBE_SIZE],
    /*
       up_left: Color,
       up_center: Color,
       up_right: Color,
       center_left: Color,
       center_center: Color,
       center_right: Color,
       down_left: Color,
       down_center: Color,
       down_right: Color
    */
}

impl Face {
    pub fn new(facelets: [[Color; CUBE_SIZE]; CUBE_SIZE]) -> Self {
        Self { facelets }
    }

    pub fn new_solid_color(color: Color) -> Self {
        Self {
            facelets: [[color; CUBE_SIZE]; CUBE_SIZE],
        }
    }

    /// Checks if the face is a single color. Returns that color if so. Otherwise, returns None.
    pub fn is_solid(&self) -> Option<Color> {
        let init = self.facelets[0][0];
        for row in 0..self.facelets.len() {
            for color in self.facelets[row] {
                if color != init {
                    return None;
                }
            }
        }

        return Some(init);
    }

    pub fn rotate(&mut self, rotation: Rotation) {
        let n: usize = self.facelets.len();
        let x: usize = n / 2;
        let y: usize = n - 1;
        for i in 0..x {
            for j in i..y - i {
                match rotation {
                    Clockwise => {
                        let k = self.facelets[i][j];
                        self.facelets[i][j] = self.facelets[y - j][i];
                        self.facelets[y - j][i] = self.facelets[y - i][y - j];
                        self.facelets[y - i][y - j] = self.facelets[j][y - i];
                        self.facelets[j][y - i] = k;
                    }
                    Counterclockwise => {
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

    pub fn get_row(&self, row_idx: usize) -> [Color; CUBE_SIZE] {
        self.facelets[row_idx]
    }

    pub fn get_col(&self, col_idx: usize) -> [Color; CUBE_SIZE] {
        [
            self.facelets[0][col_idx],
            self.facelets[1][col_idx],
            self.facelets[2][col_idx],
        ]
    }

    pub fn get_row_iter(&self) -> Iter<'_, [Color; CUBE_SIZE]> {
        self.facelets.iter()
    }

    /// The 'size' of the cube is defined as the number of facelets in a row. In other words, a
    /// standard 3x3 Rubik's cube has a size of '3'. This is defined on both cube and face for
    /// convenience.
    pub fn get_size(&self) -> usize {
        self.facelets.len()
    }
}

#[derive(Debug, PartialEq, Clone, Copy, TryFromPrimitive)]
#[repr(u8)]
pub enum CubeFace {
    Up = 0,
    Left = 1,
    Front = 2,
    Down = 3,
    Right = 4,
    Back = 5,
}

impl CubeFace {
    #[inline]
    fn us(self) -> usize {
        self as usize
    }
}

/// The cube has NUM_SIDES faces. They are organized in the array such that opposite faces are separated
/// by 3. In pictures:
///    U        0
///  L F R B  1 2 4 5    where U = up, L = left, F = front, B = back, D = down
///    D        3
#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct Cube {
    faces: [Face; NUM_SIDES], // up: Face,
                              // left: Face,
                              // front: Face,
                              // down: Face,
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
    /// Create a new cube with the provided faces
    pub fn new(faces: [Face; NUM_SIDES]) -> Self {
        Self { faces }
    }

    /// A 'solved' cube is one where all sides are a single color, and all colors are represented.
    /// We do not care which color is on which side, only that each side is a different color.
    pub fn is_solved(&self) -> bool {
        let mut colors = HashSet::new();
        for face in self.get_face_iter() {
            if let Some(color) = face.is_solid() {
                colors.insert(color);
            } else {
                return false;
            };
        }

        return colors.len() == NUM_SIDES;
    }

    /// Get an iterator of the faces on the cube.
    pub fn get_face_iter(&self) -> Iter<'_, Face> {
        self.faces.iter()
    }

    /// Rotate a face in a specific direction. This takes a [CubeFace] rather than an index.
    pub fn rotate_face(&mut self, mv: CubeMove) {
        let idx = mv.face as usize;
        self.faces[idx].rotate(mv.rotation);
        let [up, right, down, left] = SURROUNDING_FACES[idx];
        match (mv.face, mv.rotation) {
            (Up, Clockwise) => {
                let mut tmp = self.swap_row_to_row(Back.us(), 0, Right.us(), 0, false);
                tmp.reverse();
                self.swap_row_to_row(Left.us(), 0, Back.us(), 0, true);
                self.swap_row_to_row(Front.us(), 0, Left.us(), 0, false);
                self.set_row(Front.us(), 0, tmp);
            }
            (Up, Counterclockwise) => {
                let mut tmp = self.swap_row_to_row(Back.us(), 0, Left.us(), 0, true);
                self.swap_row_to_row(Right.us(), 0, Back.us(), 0, false);
                self.swap_row_to_row(Front.us(), 0, Right.us(), 0, true);
                self.set_row(Front.us(), 0, tmp);
            }
            (Left, Clockwise) => {}
            (Left, Counterclockwise) => {}
            (Front, Clockwise) => {
                let mut tmp = self.swap_row_to_col(Up.us(), 2, Right.us(), 0, false);
                tmp.reverse();
                self.swap_col_to_row(Left.us(), 2, Up.us(), 2, true);
                self.swap_row_to_col(Down.us(), 0, Left.us(), 2, false);
                self.set_row(Down.us(), 0, tmp);
            }
            (Front, Counterclockwise) => {
                let tmp = self.swap_row_to_col(Up.us(), 2, Left.us(), 2, true);
                self.swap_col_to_row(Right.us(), 0, Up.us(), 2, false);
                self.swap_row_to_col(Down.us(), 0, Right.us(), 0, true);
                self.set_row(Down.us(), 0, tmp);
            }
            (Down, Clockwise) => {}
            (Down, Counterclockwise) => {}
            (Right, Clockwise) => {}
            (Right, Counterclockwise) => {}
            (Back, Clockwise) => {
                let mut tmp = self.swap_row_to_col(Up.us(), 0, Left.us(), 0, false);
                tmp.reverse();
                self.swap_col_to_row(Right.us(), 2, Up.us(), 0, true);
                self.swap_row_to_col(Down.us(), 2, Right.us(), 2, false);
                self.set_row(Down.us(), 2, tmp);
            }
            (Back, Counterclockwise) => {
                let tmp = self.swap_row_to_col(Up.us(), 0, Right.us(), 2, true);
                self.swap_col_to_row(Left.us(), 0, Up.us(), 0, false);
                self.swap_row_to_col(Down.us(), 2, Left.us(), 0, true);
                self.set_row(Down.us(), 2, tmp);
            }
        }
        // match mv.rotation {
        //     Clockwise => {
        //         for i in 0..=2 {
        //             let tmp = self.swap_facelets(
        //                 up,
        //                 SURROUNDING_INDICES[0][i],
        //                 right,
        //                 SURROUNDING_INDICES[1][i],
        //             );
        //             self.swap_facelets(
        //                 left,
        //                 SURROUNDING_INDICES[3][i],
        //                 up,
        //                 SURROUNDING_INDICES[0][i],
        //             );
        //             self.swap_facelets(
        //                 down,
        //                 SURROUNDING_INDICES[2][i],
        //                 left,
        //                 SURROUNDING_INDICES[3][i],
        //             );
        //             self.set_facelet(down, SURROUNDING_INDICES[2][i], tmp);
        //         }
        //     }
        //     Counterclockwise => {
        //         for i in 0..=2 {
        //             let tmp = self.swap_facelets(
        //                 up,
        //                 SURROUNDING_INDICES[0][i],
        //                 left,
        //                 SURROUNDING_INDICES[3][i],
        //             );
        //             self.swap_facelets(
        //                 right,
        //                 SURROUNDING_INDICES[1][i],
        //                 up,
        //                 SURROUNDING_INDICES[0][i],
        //             );
        //             self.swap_facelets(
        //                 down,
        //                 SURROUNDING_INDICES[2][i],
        //                 right,
        //                 SURROUNDING_INDICES[1][i],
        //             );
        //             self.set_facelet(down, SURROUNDING_INDICES[2][i], tmp);
        //         }
        //     }
        // }
    }

    /// A convenience function for swapping a source face to a destination on the cube. The original
    /// value of the destination face is returned.
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

    /// A convenience function for swapping a row on one face to a column on another face.
    /// The original of the destination face is returned.
    fn swap_row_to_col(
        &mut self,
        src_face: usize,
        src_row: usize,
        dest_face: usize,
        dest_col: usize,
        reverse: bool,
    ) -> [Color; CUBE_SIZE] {
        let mut orig_src = self.faces[src_face].get_row(src_row);
        if reverse {
            orig_src.reverse();
        }
        let mut orig_dest = [White; CUBE_SIZE];

        for idx in 0..CUBE_SIZE {
            orig_dest[idx] = self.faces[dest_face].facelets[idx][dest_col];
            self.faces[dest_face].facelets[idx][dest_col] = orig_src[idx];
        }
        orig_dest
    }

    /// A convenience function for swapping a column on one face to a row on another face.
    /// The original row of the destination face is returned.
    fn swap_col_to_row(
        &mut self,
        src_face: usize,
        src_col: usize,
        dest_face: usize,
        dest_row: usize,
        reverse: bool,
    ) -> [Color; CUBE_SIZE] {
        let mut orig_src = self.faces[src_face].get_col(src_col);
        if reverse {
            orig_src.reverse();
        }
        let mut orig_dest = self.faces[dest_face].facelets[dest_row];

        for idx in 0..CUBE_SIZE {
            self.faces[dest_face].facelets[dest_row][idx] = orig_src[idx];
        }
        orig_dest
    }

    /// A convenience function for swapping a row on one face to a column on another face.
    /// The original of the destination face is returned.
    fn swap_row_to_row(
        &mut self,
        src_face: usize,
        src_row: usize,
        dest_face: usize,
        dest_row: usize,
        reverse: bool,
    ) -> [Color; CUBE_SIZE] {
        let mut orig_dest = self.faces[dest_face].facelets[dest_row];
        self.faces[dest_face].facelets[dest_row] = self.faces[src_face].facelets[src_row];
        if reverse {
            self.faces[dest_face].facelets[dest_row].reverse();
        }
        orig_dest
    }

    /// Set a face to a specific color
    fn set_facelet(&mut self, dest_face: usize, dest_facelet: (usize, usize), color: Color) {
        let dest_face_ref = self.faces[dest_face].borrow_mut();
        dest_face_ref.facelets[dest_facelet.0][dest_facelet.1] = color;
    }

    /// Replaces a row on the cube.
    #[inline]
    fn set_row(&mut self, dest_face: usize, dest_row: usize, row: [Color; CUBE_SIZE]) {
        self.faces[dest_face].facelets[dest_row] = row;
    }

    // Replaces a column on the cube
    #[inline]
    fn set_col(&mut self, dest_face: usize, dest_col: usize, col: [Color; CUBE_SIZE]) {
        for idx in 0..CUBE_SIZE {
            self.faces[dest_face].facelets[idx][dest_col] = col[idx];
        }
    }

    /// The 'size' of the cube is defined as the number of facelets in a row. In other words, a
    /// standard 3x3 Rubik's cube has a size of '3'.
    pub fn get_size(&self) -> usize {
        self.faces[0].get_size()
    }

    /// Randomize a cube by performing a series of random rotations on the cube. The method takes a
    /// minimum and maximum for the number of moves. The number of moves will be a random value
    /// between the minimum and maximum. The face and rotation direction for each move will be
    /// random.
    pub fn randomize(&mut self, min_moves: usize, max_moves: usize) {
        let mut rng = thread_rng();
        let num_moves = rng.gen_range(min_moves..=max_moves);
        // let mut moves = Vec::with_capacity(num_moves);
        for _ in 0..num_moves {
            let face = CubeFace::try_from(rng.gen_range(Up as u8..=Back as u8)).unwrap();
            let rotation =
                Rotation::try_from(rng.gen_range(Clockwise as u8..=Counterclockwise as u8))
                    .unwrap();

            self.rotate_face(CubeMove { face, rotation });
        }
    }

    /// The naive solution is a simple BFS. The cube has 12 possible moves (rotate each face in
    /// either direction). At each iteration, we check if the cube is solved. Then we perform all
    /// 12 moves, and recurse into the resulting cubes.
    pub fn solve_recursive(&self) -> Vec<CubeMove> {
        let mut queue: VecDeque<(Cube, Vec<CubeMove>)> = VecDeque::new();
        let moves: Vec<CubeMove> = Vec::new();
        let mut tested: HashSet<Cube> = HashSet::new();
        tested.insert(self.clone());
        queue.push_front((self.clone(), moves));
        while let Some((cube, moves)) = queue.pop_front() {
            if cube.is_solved() {
                return moves;
            }
            for mv in ALL_MOVES {
                let mut cube = cube.clone();
                cube.rotate_face(mv);
                if !tested.contains(&cube) {
                    tested.insert(cube.clone());
                    let mut moves = moves.clone();
                    moves.push(mv);
                    queue.push_back((cube, moves));
                }
            }
        }

        return Vec::new();
    }
}

const ALL_MOVES: [CubeMove; 12] = [
    CubeMove {
        face: Up,
        rotation: Clockwise,
    },
    CubeMove {
        face: Up,
        rotation: Counterclockwise,
    },
    CubeMove {
        face: Left,
        rotation: Clockwise,
    },
    CubeMove {
        face: Left,
        rotation: Counterclockwise,
    },
    CubeMove {
        face: Front,
        rotation: Clockwise,
    },
    CubeMove {
        face: Front,
        rotation: Counterclockwise,
    },
    CubeMove {
        face: Down,
        rotation: Clockwise,
    },
    CubeMove {
        face: Down,
        rotation: Counterclockwise,
    },
    CubeMove {
        face: Right,
        rotation: Clockwise,
    },
    CubeMove {
        face: Right,
        rotation: Counterclockwise,
    },
    CubeMove {
        face: Back,
        rotation: Clockwise,
    },
    CubeMove {
        face: Back,
        rotation: Counterclockwise,
    },
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cube;
    use crate::cube::CubeFace::{Back, Down, Front, Left, Right, Up};
    use crate::cube::Rotation::{Clockwise, Counterclockwise};
    use arr_macro::arr;
    use enum_iterator::all;
    use quickcheck::{Arbitrary, Gen};
    use rstest::rstest;
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

    #[quickcheck_macros::quickcheck]
    fn rotate_face_clockwise_90_degrees(face: Face) -> bool {
        let expected = face.clone();
        let mut result = face;
        result.rotate(Clockwise);

        // Verify that iterating the rotated face vertically from up right == iterating the
        // original face horizontally from up left
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

        // Verify that iterating the rotated face vertically from down left == iterating the
        // original face horizontally from up left
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

        // Verify that iterating the rotated face from down right == iterating the original face
        // from up left
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

    #[test]
    fn all_solid_faces_are_solid() {
        for color in all::<Color>() {
            let face = Face::new_solid_color(color);
            assert_eq!(color, face.is_solid().unwrap())
        }
    }

    #[quickcheck_macros::quickcheck]
    fn only_solid_faces_are_solid(face: Face) -> bool {
        let colors = face.get_row_iter().flatten().collect::<HashSet<&Color>>();
        (colors.len() == 1) == face.is_solid().is_some()
    }

    fn face_and_surround_match(
        result: &Cube,
        expected_front: &Face,
        expected_right: &[Color; CUBE_SIZE],
        expected_down: &[Color; CUBE_SIZE],
        expected_left: &[Color; CUBE_SIZE],
        expected_up: &[Color; CUBE_SIZE],
    ) -> bool {
        &result.faces[Front as usize] == expected_front
            && expected_right == &result.faces[Right as usize].get_col(0)
            && expected_down == &result.faces[Down as usize].get_row(0)
            && expected_left == &result.faces[Left as usize].get_col(2)
            && expected_up == &result.faces[Up as usize].get_row(2)
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_cube_face_four_times_is_identity(cube: Cube, face: CubeFace) -> bool {
        let expected = cube.clone();
        let mut result = cube;
        result.rotate_face(CubeMove {
            face,
            rotation: Clockwise,
        });
        result.rotate_face(CubeMove {
            face,
            rotation: Clockwise,
        });
        result.rotate_face(CubeMove {
            face,
            rotation: Clockwise,
        });
        result.rotate_face(CubeMove {
            face,
            rotation: Clockwise,
        });

        result == expected
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_clockwise_then_counterclockwise_is_identity(cube: Cube, face: CubeFace) -> bool {
        let expected = cube.clone();
        let mut result = cube;
        result.rotate_face(CubeMove {
            face,
            rotation: Clockwise,
        });
        result.rotate_face(CubeMove {
            face,
            rotation: Counterclockwise,
        });

        result == expected
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_cube_face_180_degrees_clockwise_equals_counterclockwise(
        cube: Cube,
        face: CubeFace,
    ) -> bool {
        let mut expected = cube.clone();
        expected.rotate_face(CubeMove {
            face,
            rotation: Counterclockwise,
        });
        expected.rotate_face(CubeMove {
            face,
            rotation: Counterclockwise,
        });
        let mut result = cube;
        result.rotate_face(CubeMove {
            face,
            rotation: Clockwise,
        });
        result.rotate_face(CubeMove {
            face,
            rotation: Clockwise,
        });

        result == expected
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_front_face_counterclockwise(cube: Cube) -> bool {
        let expected = cube.clone();
        let mut result = cube;
        result.rotate_face(CubeMove {
            face: Front,
            rotation: Counterclockwise,
        });
        let mut expected_front = expected.faces[Front as usize].clone();
        expected_front.rotate(Counterclockwise);
        let mut expected_right = expected.faces[Down as usize].get_row(0);
        expected_right.reverse();
        let expected_down = expected.faces[Left as usize].get_col(2);
        let mut expected_left = expected.faces[Up as usize].get_row(2);
        expected_left.reverse();
        let expected_up = expected.faces[Right as usize].get_col(0);

        face_and_surround_match(
            &result,
            &expected_front,
            &expected_right,
            &expected_down,
            &expected_left,
            &expected_up,
        )
    }

    #[rstest]
    #[case(Front, Up, 2, Right, 0, Down, 0, Left, 2)]
    // #[case(Left, Up, Front, Down, Back)]
    // #[case(Back, Up, 0, Left, 0, Down, Right, 2)]
    // #[case(Right, Up, Back, Down, Front)]
    // #[case(Up, Back, Right, Front, Left)]
    // #[case(Down, Front, Right, Back, Left)]
    fn rotate_face_clockwise(
        #[case] front: CubeFace,
        #[case] up: CubeFace,
        #[case] up_row: usize,
        #[case] right: CubeFace,
        #[case] right_col: usize,
        #[case] down: CubeFace,
        #[case] down_row: usize,
        #[case] left: CubeFace,
        #[case] left_col: usize,
    ) {
        let cube = Cube::default();
        let expected = cube.clone();
        let mut result = cube;
        result.rotate_face(CubeMove {
            face: Front,
            rotation: Clockwise,
        });
        let mut expected_front = expected.faces[Front as usize].clone();
        expected_front.rotate(Clockwise);
        let expected_right = expected.faces[up as usize].get_row(up_row);
        let mut expected_down = expected.faces[right as usize].get_col(right_col);
        expected_down.reverse();
        let expected_left = expected.faces[down as usize].get_row(down_row);
        let mut expected_up = expected.faces[left as usize].get_col(left_col);
        expected_up.reverse();

        assert!(face_and_surround_match(
            &result,
            &expected_front,
            &expected_right,
            &expected_down,
            &expected_left,
            &expected_up,
        ))
    }

    #[quickcheck_macros::quickcheck]
    fn cube_rotations_can_be_undone(cube: Cube, mut rotations: Vec<(CubeFace, Rotation)>) -> bool {
        let expected = cube.clone();
        let mut result = cube;
        for (face, rotation) in rotations.iter() {
            result.rotate_face(CubeMove {
                face: *face,
                rotation: *rotation,
            });
        }

        rotations.reverse();
        for (face, rotation) in rotations {
            result.rotate_face(CubeMove {
                face,
                rotation: rotation.get_opposite(),
            });
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
        let dest_facelet = (
            (dest_x % CUBE_SIZE as u8) as usize,
            (dest_y % CUBE_SIZE as u8) as usize,
        );
        let dest_face = (dest_face % NUM_SIDES as u8) as usize;
        let mut result = cube;
        result.set_facelet(dest_face, dest_facelet, color);

        result.faces[dest_face].facelets[dest_facelet.0][dest_facelet.1] == color
    }

    #[test]
    fn default_cube_is_solved() {
        assert!(Cube::default().is_solved())
    }

    #[test]
    fn solid_cube_contains_all_colors() {
        let mut colors: Vec<Color> = all().collect();
        for _ in 0..colors.len() {
            let cube = Cube::new([
                Face::new_solid_color(colors[0]),
                Face::new_solid_color(colors[1]),
                Face::new_solid_color(colors[2]),
                Face::new_solid_color(colors[3]),
                Face::new_solid_color(colors[4]),
                Face::new_solid_color(colors[5]),
            ]);
            assert!(cube.is_solved());
            colors.rotate_left(1);
        }
    }

    #[test]
    fn solid_cube_missing_a_color_is_not_solved() {
        let mut colors: Vec<Color> = all().collect();
        for _ in 0..colors.len() {
            let cube = Cube::new([
                Face::new_solid_color(colors[0]),
                Face::new_solid_color(colors[1]),
                Face::new_solid_color(colors[2]),
                Face::new_solid_color(colors[3]),
                Face::new_solid_color(colors[4]),
                Face::new_solid_color(colors[3]),
            ]);
            assert!(!cube.is_solved());
            colors.rotate_left(1);
        }
    }

    #[quickcheck_macros::quickcheck]
    fn random_cube_is_not_solved(cube: Cube) -> bool {
        // It's highly unlikely that we will get a solved cube randomly, but this will account
        // for it. A solved cube will has 6 sides with unique colors
        let solid_faces = cube
            .get_face_iter()
            .map(|face| face.is_solid())
            .flatten()
            .collect::<HashSet<Color>>();
        cube.is_solved() == (solid_faces.len() == NUM_SIDES as usize)
    }

    #[quickcheck_macros::quickcheck]
    fn swap_row_to_column(
        mut cube: Cube,
        src_face: usize,
        src_row: usize,
        dest_face: usize,
        dest_col: usize,
        reverse: bool,
    ) -> bool {
        let src_face = src_face.checked_rem(NUM_SIDES as usize).unwrap();
        let src_row = src_row.checked_rem(CUBE_SIZE).unwrap();
        let dest_face = dest_face.checked_rem(NUM_SIDES as usize).unwrap();
        let dest_col = dest_col.checked_rem(CUBE_SIZE).unwrap();

        let orig = cube.clone();
        let result = cube.swap_row_to_col(src_face, src_row, dest_face, dest_col, reverse);
        let result_is_correct = result == orig.faces[dest_face].get_col(dest_col);
        let dest_face_is_correct = (0..CUBE_SIZE).all(|idx| {
            if idx == dest_col {
                let mut expected = orig.faces[src_face].get_row(src_row);
                if reverse {
                    expected.reverse();
                }
                cube.faces[dest_face].get_col(idx) == expected
            } else {
                cube.faces[dest_face].get_col(idx) == orig.faces[dest_face].get_col(idx)
            }
        });
        result_is_correct && dest_face_is_correct
    }

    #[quickcheck_macros::quickcheck]
    fn swap_column_to_row(
        mut cube: Cube,
        src_face: usize,
        src_col: usize,
        dest_face: usize,
        dest_row: usize,
        reverse: bool,
    ) -> bool {
        let src_face = src_face.checked_rem(NUM_SIDES as usize).unwrap();
        let src_col = src_col.checked_rem(CUBE_SIZE).unwrap();
        let dest_face = dest_face.checked_rem(NUM_SIDES as usize).unwrap();
        let dest_row = dest_row.checked_rem(CUBE_SIZE).unwrap();

        let orig = cube.clone();
        let mut result = cube.swap_col_to_row(src_face, src_col, dest_face, dest_row, reverse);
        let result_is_correct = result == orig.faces[dest_face].get_row(dest_row);
        let dest_face_is_correct = (0..CUBE_SIZE).all(|idx| {
            if idx == dest_row {
                let mut col = orig.faces[src_face].get_col(src_col);
                if reverse {
                    col.reverse();
                }
                cube.faces[dest_face].get_row(idx) == col
            } else {
                cube.faces[dest_face].get_row(idx) == orig.faces[dest_face].get_row(idx)
            }
        });
        result_is_correct && dest_face_is_correct
    }

    #[quickcheck_macros::quickcheck]
    fn swap_row_to_row(
        mut cube: Cube,
        src_face: usize,
        src_col: usize,
        dest_face: usize,
        dest_row: usize,
        reverse: bool,
    ) -> bool {
        let src_face = src_face.checked_rem(NUM_SIDES as usize).unwrap();
        let src_row = src_col.checked_rem(CUBE_SIZE).unwrap();
        let dest_face = dest_face.checked_rem(NUM_SIDES as usize).unwrap();
        let dest_row = dest_row.checked_rem(CUBE_SIZE).unwrap();

        let orig = cube.clone();
        let mut result = cube.swap_row_to_row(src_face, src_row, dest_face, dest_row, reverse);
        let result_is_correct = result == orig.faces[dest_face].get_row(dest_row);
        let dest_face_is_correct = (0..CUBE_SIZE).all(|idx| {
            if idx == dest_row {
                let mut expected = orig.faces[src_face].get_row(src_row);
                if reverse {
                    expected.reverse();
                }
                cube.faces[dest_face].get_row(idx) == expected
            } else {
                cube.faces[dest_face].get_row(idx) == orig.faces[dest_face].get_row(idx)
            }
        });
        result_is_correct && dest_face_is_correct
    }

    // #[ignore]
    // fn solve_works() {
    //     let mut cube = Cube::default();
    //     let num_moves = 1;
    //     cube.randomize(num_moves, num_moves);
    //     cube.rotate_face(Front,
    //         Clockwise,
    //     });
    //     let soln = cube.solve_recursive();
    //     println!("{:?}", cube);
    //     println!("{:?}", soln);
    //     assert_eq!(num_moves, soln.len());
    // }
}
