use crate::cube::Color::*;
use crate::cube::Rotation::*;
use crate::cube::{Color, Rotation, CUBE_SIZE, NUM_SIDES};
use crate::cube3d::CubeFace::*;
use num_enum::TryFromPrimitive;
use rand::{thread_rng, Rng};
use std::collections::{HashSet, VecDeque};
use std::convert::TryFrom;

#[derive(Debug, PartialEq, Clone, Copy, TryFromPrimitive)]
#[repr(u8)]
pub enum CubeFace {
    Up = 0,
    Left = 1,
    Front = 2,
    Right = 3,
    Back = 4,
    Down = 5,
}

/// The cube is represented as a 3-dimensional array. Each element in the array represents a
/// 'facelet' or 'sticker' on one face.
/// If we consider the following array indices: `cube[x][y][z]`
/// * The 'x' index represents the 'face'. Faces are ordered according to [CubeFace].
/// * The 'y' index represents the 'row' in the face.
/// * The 'z' index represents the 'column' in the face.
#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct Cube {
    pub(crate) facelets: [[[Color; CUBE_SIZE]; CUBE_SIZE]; NUM_SIDES],
}

impl Default for Cube {
    fn default() -> Self {
        Self {
            facelets: [
                [[White; CUBE_SIZE]; CUBE_SIZE],
                [[Green; CUBE_SIZE]; CUBE_SIZE],
                [[Red; CUBE_SIZE]; CUBE_SIZE],
                [[Blue; CUBE_SIZE]; CUBE_SIZE],
                [[Orange; CUBE_SIZE]; CUBE_SIZE],
                [[Yellow; CUBE_SIZE]; CUBE_SIZE],
            ],
        }
    }
}

impl Cube {
    pub fn new(facelets: [[[Color; CUBE_SIZE]; CUBE_SIZE]; NUM_SIDES]) -> Self {
        Self { facelets }
    }

    /// A 'solved' cube is one where all sides are a single color, and all colors are represented.
    /// We do not care which color is on which side, only that each side is a different color.
    pub fn is_solved(&self) -> bool {
        let mut colors = HashSet::new();
        for face in self.facelets.iter() {
            if let Some(color) = is_solid(face) {
                colors.insert(color);
            } else {
                return false;
            };
        }

        return colors.len() == NUM_SIDES;
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
            let face = CubeFace::try_from(rng.gen_range(Up as u8..=Down as u8)).unwrap();
            let rotation =
                Rotation::try_from(rng.gen_range(Clockwise as u8..=Counterclockwise as u8))
                    .unwrap();

            self.rotate_face(face, rotation);
        }
    }

    /// Perform a single quarter-turn move on one face of the cube. Takes a face and a direction.
    /// The direction is with respect to the front of the face.
    pub fn rotate_face(&mut self, face: CubeFace, direction: Rotation) {
        match face {
            Up => self.rotate_row(0, direction),
            Left => self.rotate_column_lr(0, direction),
            Front => self.rotate_column_fb(CUBE_SIZE - 1, direction),
            Right => self.rotate_column_lr(CUBE_SIZE - 1, direction.get_opposite()),
            Back => self.rotate_column_fb(0, direction.get_opposite()),
            Down => self.rotate_row(CUBE_SIZE - 1, direction.get_opposite()),
        }
    }

    /// Rotate a horizontal 'row' in the cube
    /// Note: This method supports rotating any row, but our 'public' API ([rotate_face]) only
    /// supports rotating faces.
    fn rotate_row(&mut self, row: usize, direction: Rotation) {
        match direction {
            Clockwise => {
                let tmp = self.facelets[1][row];
                self.facelets[1][row] = self.facelets[2][row];
                self.facelets[2][row] = self.facelets[3][row];
                self.facelets[3][row] = self.facelets[4][row];
                self.facelets[4][row] = tmp;
            }
            Counterclockwise => {
                let tmp = self.facelets[1][row];
                self.facelets[1][row] = self.facelets[4][row];
                self.facelets[4][row] = self.facelets[3][row];
                self.facelets[3][row] = self.facelets[2][row];
                self.facelets[2][row] = tmp;
            }
        }

        if row == 0 {
            rotate_face_only(&mut self.facelets[0], direction);
        } else if row == CUBE_SIZE - 1 {
            rotate_face_only(&mut self.facelets[5], direction.get_opposite());
        }
    }

    /// Rotate a vertical 'column' in the cube. This method rotates columns aligned from left to
    /// right when looking at the front of the cube. (For a 3x3x3 cube, this would be the left,
    /// center, and right columns.)
    /// Note: This method supports rotating any row, but our 'public' API ([rotate_face]) only
    /// supports rotating faces.
    fn rotate_column_lr(&mut self, col: usize, direction: Rotation) {
        match direction {
            Clockwise => {
                for i in 0..CUBE_SIZE {
                    let tmp = self.facelets[5][i][col];
                    self.facelets[5][i][col] = self.facelets[2][i][col]; // down = front
                    self.facelets[2][i][col] = self.facelets[0][i][col]; // front = up
                    self.facelets[0][i][col] =
                        self.facelets[4][CUBE_SIZE - 1 - i][CUBE_SIZE - 1 - col]; // up = back
                    self.facelets[4][CUBE_SIZE - 1 - i][CUBE_SIZE - 1 - col] = tmp;
                    // back = down
                }
            }
            Counterclockwise => {
                for i in 0..CUBE_SIZE {
                    let tmp = self.facelets[0][i][col];
                    self.facelets[0][i][col] = self.facelets[2][i][col]; // up = front
                    self.facelets[2][i][col] = self.facelets[5][i][col]; // front = down
                    self.facelets[5][i][col] =
                        self.facelets[4][CUBE_SIZE - 1 - i][CUBE_SIZE - 1 - col]; // down = back
                    self.facelets[4][CUBE_SIZE - 1 - i][CUBE_SIZE - 1 - col] = tmp;
                    // back = up
                }
            }
        }

        if col == 0 {
            rotate_face_only(&mut self.facelets[1], direction);
        } else if col == CUBE_SIZE - 1 {
            rotate_face_only(&mut self.facelets[3], direction.get_opposite());
        }
    }

    /// Rotate a vertical 'column' in the cube. This method rotates columns aligned from front to
    /// back when looking at the front of the cube. (For a 3x3x3 cube, this would be the front,
    /// middle, and back columns.)
    /// Note: This method supports rotating any column, but our 'public' API ([rotate_face]) only
    /// supports rotating faces.
    fn rotate_column_fb(&mut self, col: usize, direction: Rotation) {
        match direction {
            Clockwise => {
                for i in 0..CUBE_SIZE {
                    let tmp = self.facelets[3][i][CUBE_SIZE - 1 - col];
                    self.facelets[3][i][CUBE_SIZE - 1 - col] = self.facelets[0][col][i]; // right = up
                    self.facelets[0][col][i] = self.facelets[1][i][col]; // up = left
                    self.facelets[1][i][col] = self.facelets[5][CUBE_SIZE - 1 - col][i]; // left = down (row)
                    self.facelets[5][CUBE_SIZE - 1 - col][i] = tmp;
                    // down = right
                }
                self.facelets[0][col].reverse();
                self.facelets[NUM_SIDES - 1][CUBE_SIZE - 1 - col].reverse();
            }
            Counterclockwise => {
                for i in 0..CUBE_SIZE {
                    let tmp = self.facelets[1][i][col];
                    self.facelets[1][i][col] = self.facelets[0][col][CUBE_SIZE - 1 - i]; // left = up
                    self.facelets[0][col][CUBE_SIZE - 1 - i] =
                        self.facelets[3][CUBE_SIZE - 1 - i][CUBE_SIZE - 1 - col]; // up = right
                    self.facelets[3][CUBE_SIZE - 1 - i][CUBE_SIZE - 1 - col] =
                        self.facelets[5][CUBE_SIZE - 1 - col][i]; // right = down (row)
                    self.facelets[5][CUBE_SIZE - 1 - col][i] = tmp; // down = left
                }
            }
        }

        if col == 0 {
            rotate_face_only(&mut self.facelets[4], direction.get_opposite());
        } else if col == CUBE_SIZE - 1 {
            rotate_face_only(&mut self.facelets[2], direction);
        }
    }

    /// The naive solution is a simple BFS. The cube has 12 possible moves (rotate each face in
    /// either direction). At each iteration, we check if the cube is solved. Then we perform all
    /// 12 moves, and recurse into the resulting cubes.
    pub fn solve(&self) -> Vec<(CubeFace, Rotation)> {
        let mut queue: VecDeque<(Cube, Vec<(CubeFace, Rotation)>)> = VecDeque::new();
        let moves: Vec<(CubeFace, Rotation)> = Vec::new();
        let mut tested: HashSet<Cube> = HashSet::new();
        tested.insert(self.clone());
        queue.push_front((self.clone(), moves));
        while let Some((cube, moves)) = queue.pop_front() {
            if cube.is_solved() {
                return moves;
            }
            for (face, dir) in ALL_MOVES {
                let mut cube = cube.clone();
                cube.rotate_face(face, dir);
                if !tested.contains(&cube) {
                    tested.insert(cube.clone());
                    let mut moves = moves.clone();
                    moves.push((face, dir));
                    queue.push_back((cube, moves));
                }
            }
        }

        return Vec::new();
    }

    fn get_col(&self, face_idx: usize, col_idx: usize) -> [Color; CUBE_SIZE] {
        [
            self.facelets[face_idx][0][col_idx],
            self.facelets[face_idx][1][col_idx],
            self.facelets[face_idx][2][col_idx],
        ]
    }

    fn get_row(&self, face_idx: usize, row_idx: usize) -> [Color; CUBE_SIZE] {
        self.facelets[face_idx][row_idx]
    }
}

const ALL_MOVES: [(CubeFace, Rotation); 12] = [
    (Up, Clockwise),
    (Up, Counterclockwise),
    (Left, Clockwise),
    (Left, Counterclockwise),
    (Front, Clockwise),
    (Front, Counterclockwise),
    (Down, Clockwise),
    (Down, Counterclockwise),
    (Right, Clockwise),
    (Right, Counterclockwise),
    (Back, Clockwise),
    (Back, Counterclockwise),
];

fn rotate_face_only(face: &mut [[Color; CUBE_SIZE]; CUBE_SIZE], direction: Rotation) {
    let n: usize = CUBE_SIZE;
    let x: usize = n / 2;
    let y: usize = n - 1;
    for i in 0..x {
        for j in i..y - i {
            match direction {
                Clockwise => {
                    let k = face[i][j];
                    face[i][j] = face[y - j][i];
                    face[y - j][i] = face[y - i][y - j];
                    face[y - i][y - j] = face[j][y - i];
                    face[j][y - i] = k;
                }
                Counterclockwise => {
                    let k = face[i][j];
                    face[i][j] = face[j][y - i];
                    face[j][y - i] = face[y - i][y - j];
                    face[y - i][y - j] = face[y - j][i];
                    face[y - j][i] = k;
                }
            }
        }
    }
}

fn new_solid_color(color: Color) -> [[Color; CUBE_SIZE]; CUBE_SIZE] {
    [[color; CUBE_SIZE]; CUBE_SIZE]
}

/// Checks if the face is a single color. Returns that color if so. Otherwise, returns None.
fn is_solid(face: &[[Color; CUBE_SIZE]; CUBE_SIZE]) -> Option<Color> {
    let init = face[0][0];
    for row in 0..CUBE_SIZE {
        for color in face[row] {
            if color != init {
                return None;
            }
        }
    }

    return Some(init);
}

#[cfg(test)]
mod tests {
    use crate::cube::Rotation::{Clockwise, Counterclockwise};
    use crate::cube::{Color, Face, Rotation, CUBE_SIZE};
    use crate::cube3d;
    use crate::cube3d::CubeFace::{Back, Down, Front, Left, Right, Up};
    use crate::cube3d::{is_solid, new_solid_color, rotate_face_only, Cube, CubeFace};
    use arr_macro::arr;
    use enum_iterator::all;
    use quickcheck::{Arbitrary, Gen};
    use quickcheck_macros::quickcheck;
    use rstest::rstest;
    use std::collections::HashSet;
    use std::mem;

    const NUM_SIDES: u8 = crate::cube::NUM_SIDES as u8;

    impl Arbitrary for CubeFace {
        fn arbitrary(g: &mut Gen) -> Self {
            let val = u8::arbitrary(g) % NUM_SIDES;
            unsafe { mem::transmute(val as u8) }
        }
    }

    impl Arbitrary for Cube {
        fn arbitrary(g: &mut Gen) -> Cube {
            Cube {
                facelets: { arr![arr![arr![Color::arbitrary(g); 3]; 3]; 6] },
            }
        }
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_face_only_rotates_the_face(cube: Cube, direction: Rotation) -> bool {
        let expected = cube.clone().facelets[0];
        let mut result = cube.facelets[0];
        rotate_face_only(&mut result, direction);

        result != expected || is_solid(&result).is_some()
    }

    fn face_and_surround_match(
        result: &Cube,
        expected_front: &[[Color; CUBE_SIZE]; CUBE_SIZE],
        expected_right: &[Color; CUBE_SIZE],
        expected_down: &[Color; CUBE_SIZE],
        expected_left: &[Color; CUBE_SIZE],
        expected_up: &[Color; CUBE_SIZE],
    ) -> bool {
        &result.facelets[Front as usize] == expected_front
            && expected_right == &result.get_col(Right as usize, 0)
            && expected_down == &result.get_row(Down as usize, 0)
            && expected_left == &result.get_col(Left as usize, 2)
            && expected_up == &result.get_row(Up as usize, 2)
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
        let mut expected_front = expected.facelets[Front as usize].clone();
        rotate_face_only(&mut expected_front, Counterclockwise);
        let mut expected_right = expected.get_row(Down as usize, 0);
        expected_right.reverse();
        let expected_down = expected.get_col(Left as usize, 2);
        let mut expected_left = expected.get_row(Up as usize, 2);
        expected_left.reverse();
        let expected_up = expected.get_col(Right as usize, 0);

        face_and_surround_match(
            &result,
            &expected_front,
            &expected_right,
            &expected_down,
            &expected_left,
            &expected_up,
        )
    }

    fn get_next_face(face: CubeFace) -> usize {
        let mut next = CubeFace::try_from((face as u8 + 1) % NUM_SIDES).unwrap();
        if next == Up || next == Down {
            next = Left;
        }
        next as usize
    }

    fn get_prev_face(face: CubeFace) -> usize {
        let mut next = CubeFace::try_from((face as u8 - 1) % NUM_SIDES).unwrap();
        if next == Up || next == Down {
            next = Back;
        }
        next as usize
    }

    #[rstest]
    #[case(Front, |cube: &Cube| cube.get_row(0, 2), |cube: &Cube| cube.get_row(5, 0), true, false, true, false)]
    #[case(Back, |cube: &Cube| cube.get_row(0, 0), |cube: &Cube| cube.get_row(5, 2), false, true, false, true)]
    #[case(Left, |cube: &Cube| cube.get_col(0, 0), |cube: &Cube| cube.get_col(5, 0), true, false, false, true)]
    #[case(Right, |cube: &Cube| cube.get_col(0, 2), |cube: &Cube| cube.get_col(5, 2), false, true, true, false)]
    fn rotate_face_clockwise(
        #[case] front: CubeFace,
        #[case] up: fn(&Cube) -> [Color; 3],
        #[case] down: fn(&Cube) -> [Color; 3],
        #[case] reverse_up: bool,
        #[case] reverse_right: bool,
        #[case] reverse_down: bool,
        #[case] reverse_left: bool,
    ) {
        let mut cube = Cube::default();
        cube.randomize(20, 20);
        let expected = cube.clone();
        let mut result = cube;
        result.rotate_face(front, Clockwise);
        let mut expected_front = expected.facelets[front as usize].clone();
        rotate_face_only(&mut expected_front, Clockwise);
        assert_eq!(expected_front, result.facelets[front as usize]);
        let mut expected_right = up(&expected);
        if reverse_right {
            expected_right.reverse();
        }
        assert_eq!(expected_right, result.get_col(get_next_face(front), 0));
        let mut expected_down = expected.get_col(get_next_face(front), 0);
        if reverse_down {
            expected_down.reverse();
        }
        assert_eq!(expected_down, down(&result));
        let mut expected_left = down(&expected);
        if reverse_left {
            expected_left.reverse();
        }
        assert_eq!(expected_left, result.get_col(get_prev_face(front), 2));
        let mut expected_up = expected.get_col(get_prev_face(front), 2);
        if reverse_up {
            expected_up.reverse();
        }
        assert_eq!(expected_up, up(&result));
    }

    #[test]
    fn rotate_face_clockwise_up() {
        let mut cube = Cube::default();
        cube.randomize(20, 20);
        let expected = cube.clone();
        let mut result = cube;
        result.rotate_face(Up, Clockwise);
        let mut expected_front = expected.facelets[Up as usize].clone();
        rotate_face_only(&mut expected_front, Clockwise);
        assert_eq!(expected_front, result.facelets[Up as usize]);
        let mut expected_right = expected.get_row(Back as usize, 0);
        assert_eq!(expected_right, result.get_row(Right as usize, 0));
        let mut expected_down = expected.get_row(Right as usize, 0);
        assert_eq!(expected_down, result.get_row(Front as usize, 0));
        let mut expected_left = expected.get_row(Front as usize, 0);
        assert_eq!(expected_left, result.get_row(Left as usize, 0));
        let mut expected_up = expected.get_row(Left as usize, 0);
        assert_eq!(expected_up, result.get_row(Back as usize, 0));
    }

    #[test]
    fn rotate_face_clockwise_down() {
        let mut cube = Cube::default();
        cube.randomize(20, 20);
        let expected = cube.clone();
        let mut result = cube;
        result.rotate_face(Down, Clockwise);
        let mut expected_front = expected.facelets[Down as usize].clone();
        rotate_face_only(&mut expected_front, Clockwise);
        assert_eq!(expected_front, result.facelets[Down as usize]);
        let mut expected_right = expected.get_row(Front as usize, 2);
        assert_eq!(expected_right, result.get_row(Right as usize, 2));
        let mut expected_down = expected.get_row(Right as usize, 2);
        assert_eq!(expected_down, result.get_row(Back as usize, 2));
        let mut expected_left = expected.get_row(Back as usize, 2);
        assert_eq!(expected_left, result.get_row(Left as usize, 2));
        let mut expected_up = expected.get_row(Left as usize, 2);
        assert_eq!(expected_up, result.get_row(Front as usize, 2));
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

    #[test]
    fn default_cube_is_solved() {
        assert!(Cube::default().is_solved())
    }

    #[test]
    fn solid_cube_contains_all_colors() {
        let mut colors: Vec<Color> = all().collect();
        for _ in 0..colors.len() {
            let cube = Cube::new([
                new_solid_color(colors[0]),
                new_solid_color(colors[1]),
                new_solid_color(colors[2]),
                new_solid_color(colors[3]),
                new_solid_color(colors[4]),
                new_solid_color(colors[5]),
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
                new_solid_color(colors[0]),
                new_solid_color(colors[1]),
                new_solid_color(colors[2]),
                new_solid_color(colors[3]),
                new_solid_color(colors[4]),
                new_solid_color(colors[3]),
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
            .facelets
            .iter()
            .map(|face| is_solid(&face))
            .flatten()
            .collect::<HashSet<Color>>();
        cube.is_solved() == (solid_faces.len() == NUM_SIDES as usize)
    }

    #[test]
    fn solve_works() {
        let mut cube = Cube::default();
        let num_moves = 3;
        cube.randomize(num_moves, num_moves);
        let soln = cube.solve();
        println!("{:?}", cube);
        println!("{:?}", soln);
        assert!(soln.len() <= num_moves);
    }
}
