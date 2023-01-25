use crate::cube::CubeFace::*;
use crate::cube::Direction::{Clockwise, Counterclockwise};
use enum_iterator::Sequence;
use num_enum::TryFromPrimitive;
use rand::{thread_rng, Rng};

/// The number of sides of the puzzle. In the case of a cube, the number of sides is 6. It is
/// unlikely that this solution can be generalized for other shapes, but it would be an interesting
/// exercise some day.
pub const NUM_SIDES: usize = 6;

/// In Rubik's cube terminology, 'size' is frequently used to refer to length of each side of the
/// cube in terms of number of facelets. For example, a standard 3x3x3 cube is said to be of size
/// '3'.
pub const CUBE_SIZE: usize = 3;

/// The cube has 6 colors, one for each side. Arbitrarily, we define whi
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

/// We are considering only face turns as moves, and are using the
/// "quarter turn metric".
#[derive(Debug, PartialEq, Clone, Copy, TryFromPrimitive)]
#[repr(u8)]
pub enum Direction {
    Clockwise = 0,
    Counterclockwise,
}

impl Direction {
    pub fn get_opposite(self) -> Self {
        match self {
            Clockwise => Counterclockwise,
            Counterclockwise => Clockwise,
        }
    }
}

/// One of six faces or sides to the cube.
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

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct CubeMove {
    pub face: CubeFace,
    pub direction: Direction,
}

impl CubeMove {
    pub fn new(face: CubeFace, direction: Direction) -> Self {
        Self { face, direction }
    }
}

pub const ALL_MOVES: [CubeMove; 12] = [
    CubeMove {
        face: Up,
        direction: Clockwise,
    },
    CubeMove {
        face: Up,
        direction: Counterclockwise,
    },
    CubeMove {
        face: Left,
        direction: Clockwise,
    },
    CubeMove {
        face: Left,
        direction: Counterclockwise,
    },
    CubeMove {
        face: Front,
        direction: Clockwise,
    },
    CubeMove {
        face: Front,
        direction: Counterclockwise,
    },
    CubeMove {
        face: Down,
        direction: Clockwise,
    },
    CubeMove {
        face: Down,
        direction: Counterclockwise,
    },
    CubeMove {
        face: Right,
        direction: Clockwise,
    },
    CubeMove {
        face: Right,
        direction: Counterclockwise,
    },
    CubeMove {
        face: Back,
        direction: Clockwise,
    },
    CubeMove {
        face: Back,
        direction: Counterclockwise,
    },
];

pub trait Cube {
    /// A 'solved' cube is one where all sides are a single color, and all colors are represented.
    /// We do not care which color is on which side, only that each side is a different color.
    fn is_solved(&self) -> bool;

    fn cube_move(&mut self, cube_move: CubeMove);
    /// Randomize a cube by performing a series of random rotations on the cube. The method takes a
    /// minimum and maximum for the number of moves. The number of moves will be a random value
    /// between the minimum and maximum. The face and rotation direction for each move will be
    /// random.
    fn randomize(&mut self, min_moves: usize, max_moves: usize) {
        let mut rng = thread_rng();
        let num_moves = rng.gen_range(min_moves..=max_moves);
        // let mut moves = Vec::with_capacity(num_moves);
        for _ in 0..num_moves {
            let face = CubeFace::try_from(rng.gen_range(Up as u8..=Down as u8)).unwrap();
            let rotation =
                Direction::try_from(rng.gen_range(Clockwise as u8..=Counterclockwise as u8))
                    .unwrap();

            self.cube_move(CubeMove::new(face, rotation));
        }
    }
    /// The naive solution is a simple BFS. The cube has 12 possible moves (rotate each face in
    /// either direction). At each iteration, we check if the cube is solved. Then we perform all
    /// 12 moves, and recurse into the resulting cubes.
    fn solve(&self) -> Vec<CubeMove>;
}

#[cfg(test)]
mod tests {
    use crate::cube::{Color, CubeFace, Direction};
    use quickcheck::{Arbitrary, Gen};
    use std::mem;

    const NUM_SIDES: u8 = crate::cube::NUM_SIDES as u8;

    impl Arbitrary for Direction {
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
}
