use crate::cube::CubeFace::*;
use crate::cube::Direction::{Clockwise, Counterclockwise};
use enum_iterator::Sequence;
use num_enum::TryFromPrimitive;
use rand::{thread_rng, Rng};
use std::collections::{HashSet, VecDeque};
use std::fmt::{Debug, Formatter};
use std::hash::Hash;
use termion::color;

/// The number of sides of the puzzle. In the case of a cube, the number of sides is 6. It is
/// unlikely that this solution can be generalized for other shapes, but it would be an interesting
/// exercise some day.
pub const NUM_SIDES: usize = 6;

/// In Rubik's cube terminology, 'size' is frequently used to refer to length of each side of the
/// cube in terms of number of facelets. For example, a standard 3x3x3 cube is said to be of size
/// '3'.
pub const CUBE_SIZE: usize = 3;

/// The cube has 6 colors, one for each side. Arbitrarily, we define whi
#[derive(Clone, Copy, PartialEq, Eq, Hash, Sequence)]
#[repr(u8)]
pub enum Color {
    White = 1,
    Yellow,
    Blue,
    Green,
    Red,
    Orange,
}

impl Debug for Color {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Color::White => f.write_fmt(format_args!("{}W", color::Fg(color::White))),
            Color::Yellow => f.write_fmt(format_args!("{}Y", color::Fg(color::Yellow))),
            Color::Blue => f.write_fmt(format_args!("{}B", color::Fg(color::Blue))),
            Color::Green => f.write_fmt(format_args!("{}G", color::Fg(color::Green))),
            Color::Red => f.write_fmt(format_args!("{}R", color::Fg(color::Red))),
            Color::Orange => f.write_fmt(format_args!("{}O", color::Fg(color::Rgb(255, 165, 0)))),
        }
    }
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
#[derive(Debug, PartialEq, Clone, Copy, Eq, Hash, Sequence, TryFromPrimitive)]
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
    U_MOVE,
    U_PRIME_MOVE,
    L_MOVE,
    L_PRIME_MOVE,
    F_MOVE,
    F_PRIME_MOVE,
    D_MOVE,
    D_PRIME_MOVE,
    R_MOVE,
    R_PRIME_MOVE,
    B_MOVE,
    B_PRIME_MOVE,
];

pub const U_MOVE: CubeMove = CubeMove {
    face: Up,
    direction: Clockwise,
};

pub const U_PRIME_MOVE: CubeMove = CubeMove {
    face: Up,
    direction: Counterclockwise,
};

pub const L_MOVE: CubeMove = CubeMove {
    face: Left,
    direction: Clockwise,
};

pub const L_PRIME_MOVE: CubeMove = CubeMove {
    face: Left,
    direction: Counterclockwise,
};

pub const F_MOVE: CubeMove = CubeMove {
    face: Front,
    direction: Clockwise,
};

pub const F_PRIME_MOVE: CubeMove = CubeMove {
    face: Front,
    direction: Counterclockwise,
};

pub const D_MOVE: CubeMove = CubeMove {
    face: Down,
    direction: Clockwise,
};
pub const D_PRIME_MOVE: CubeMove = CubeMove {
    face: Down,
    direction: Counterclockwise,
};

pub const R_MOVE: CubeMove = CubeMove {
    face: Right,
    direction: Clockwise,
};

pub const R_PRIME_MOVE: CubeMove = CubeMove {
    face: Right,
    direction: Counterclockwise,
};

pub const B_MOVE: CubeMove = CubeMove {
    face: Back,
    direction: Clockwise,
};

pub const B_PRIME_MOVE: CubeMove = CubeMove {
    face: Back,
    direction: Counterclockwise,
};

pub trait Cube: Sized + Clone + Eq + Hash {
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
    fn solve(&self) -> Vec<CubeMove> {
        let mut queue: VecDeque<(Self, Vec<CubeMove>)> = VecDeque::new();
        let moves: Vec<CubeMove> = Vec::new();
        let mut tested: HashSet<Self> = HashSet::new();
        tested.insert(self.clone());
        queue.push_front((self.clone(), moves));
        while let Some((cube, moves)) = queue.pop_front() {
            if cube.is_solved() {
                return moves;
            }
            for cubeMove in ALL_MOVES {
                let mut cube = cube.clone();
                cube.cube_move(cubeMove);
                if !tested.contains(&cube) {
                    tested.insert(cube.clone());
                    let mut moves = moves.clone();
                    moves.push(cubeMove);
                    queue.push_back((cube, moves));
                }
            }
        }

        return Vec::new();
    }

    /// The naive solution is a simple BFS. The cube has 12 possible moves (rotate each face in
    /// either direction). At each iteration, we check if the cube is solved. Then we perform all
    /// 12 moves, and recurse into the resulting cubes.
    fn solve_bfs_nocache(&self) -> Vec<CubeMove> {
        let mut queue: VecDeque<(Self, Vec<CubeMove>)> = VecDeque::new();
        let moves: Vec<CubeMove> = Vec::new();

        queue.push_front((self.clone(), moves));
        while let Some((cube, moves)) = queue.pop_front() {
            if cube.is_solved() {
                return moves;
            }
            for cubeMove in ALL_MOVES {
                let mut cube = cube.clone();
                cube.cube_move(cubeMove);
                let mut moves = moves.clone();
                moves.push(cubeMove);
                queue.push_back((cube, moves));
            }
        }

        return Vec::new();
    }

    fn solve_iddfs(&self) -> Option<Vec<CubeMove>> {
        for depth in 0..26 {
            let (found, remaining) = self.dls(Vec::new(), depth);
            if found != None {
                return found;
            } else if !remaining {
                return None;
            }
        }

        None
    }

    fn dls(&self, moves: Vec<CubeMove>, depth: u8) -> (Option<Vec<CubeMove>>, bool) {
        return if depth == 0 {
            if self.is_solved() {
                (Some(moves), true)
            } else {
                (None, true) //   (Not found, but may have children)
            }
        } else {
            let mut any_remaining = false;
            for cubeMove in ALL_MOVES {
                let mut moves = moves.clone();
                moves.push(cubeMove);
                let mut cube = self.clone();
                cube.cube_move(cubeMove);
                let (found, remaining) = cube.dls(moves, depth - 1);
                if found != None {
                    return (found, true);
                }
                any_remaining = remaining; // (At least one node found at depth, let IDDFS deepen)
            }
            (None, any_remaining)
        };
    }
}

#[cfg(test)]
mod tests {
    use crate::cube::{Color, CubeFace, CubeMove, Direction};
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

    impl Arbitrary for CubeMove {
        fn arbitrary(g: &mut Gen) -> Self {
            Self {
                face: CubeFace::arbitrary(g),
                direction: Direction::arbitrary(g),
            }
        }
    }

    impl Arbitrary for Color {
        fn arbitrary(g: &mut Gen) -> Self {
            let val = (u8::arbitrary(g) % NUM_SIDES) + 1;
            unsafe { mem::transmute(val as u8) }
        }
    }
}
