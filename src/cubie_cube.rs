//! A 'cubie'-centric implementation of the cube.
//!
//!
//!

use crate::cube::{Color, Cube, CubeFace, CubeMove, Direction};

const NUM_EDGES: usize = 12;
const NUM_CORNERS: usize = 8;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
enum EdgeOrientation {
    Oriented = 0,
    Flipped = 1,
}

/// An 'edge' of the cube. Each edge has two facelets. They facelet corresponding to index 0
/// represents the primary direction of the edge. For edges along the up/down faces the primary
/// direction is up/down. For edges along the left/right faces, the primary direction is left/right.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Edge {
    facelets: [Color; 2],
    orientation: EdgeOrientation,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
enum CornerOrientation {
    Oriented = 0,
    RotatedCW = 1,
    RotatedCCW = 2,
}

/// A 'corner' of the cube. Each corner has 3 facelets. The facelet corresponding to index 0
/// represents the primary direction of the corner. The primary direction of each corner is always
/// up/down. The other facelets are arranged in clockwise order in the array.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Corner {
    facelets: [Color; 3],
    orientation: CornerOrientation,
}

struct CubieCube {
    edges: [Edge; NUM_EDGES],
    corners: [Corner; NUM_CORNERS],
}

impl Cube for CubieCube {
    fn is_solved(&self) -> bool {
        todo!()
    }

    fn cube_move(&mut self, cube_move: CubeMove) {
        todo!()
    }

    fn randomize(&mut self, min_moves: usize, max_moves: usize) {
        todo!()
    }

    fn solve(&self) -> Vec<CubeMove> {
        todo!()
    }
}
