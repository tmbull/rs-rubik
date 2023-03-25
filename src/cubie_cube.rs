//! A 'cubie'-centric implementation of the cube.
//!
//!
//!

use crate::array_cube::ArrayCube;
use crate::cube::Color::{Blue, Green, Orange, Red, White, Yellow};
use crate::cube::CubeFace::*;
use crate::cube::Direction::{Clockwise, Counterclockwise};
use crate::cube::{Color, Cube, CubeFace, CubeMove, Direction, CUBE_SIZE, NUM_SIDES};
use crate::cubie_cube::CornersIndex::{DLB, DLF, DRB, DRF, ULB, ULF, URB, URF};
use crate::cubie_cube::EdgesIndex::{BL, BR, DB, DF, DL, DR, FL, FR, UB, UF, UL, UR};
use enum_iterator::{all, Sequence};
use num_enum::TryFromPrimitive;

const NUM_EDGES: usize = 12;
const NUM_CORNERS: usize = 8;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, TryFromPrimitive)]
#[repr(u8)]
enum EdgeOrientation {
    Oriented = 0,
    Flipped = 1,
}

/// An 'edge' of the cube. Each edge has two facelets. They facelet corresponding to index 0
/// represents the primary direction of the edge. For edges along the up/down faces the primary
/// direction is up/down. For edges along the left/right faces, the primary direction is left/right.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Edge {
    index: EdgesIndex,
    orientation: EdgeOrientation,
}

impl Default for Edge {
    fn default() -> Self {
        Self {
            index: UB,
            orientation: EdgeOrientation::Oriented,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Sequence, TryFromPrimitive)]
#[repr(u8)]
pub enum EdgesIndex {
    UB = 0,
    UR = 1,
    UF = 2,
    UL = 3,
    FR = 4,
    FL = 5,
    BL = 6,
    BR = 7,
    DF = 8,
    DL = 9,
    DB = 10,
    DR = 11,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, TryFromPrimitive)]
#[repr(u8)]
pub enum CornerOrientation {
    Oriented = 0,
    RotatedCW = 1,
    RotatedCCW = 2,
}

/// A 'corner' of the cube. Each corner has 3 facelets. The facelet corresponding to index 0
/// represents the primary direction of the corner. The primary direction of each corner is always
/// up/down. The other facelets are arranged in clockwise order in the array.
#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
pub struct Corner {
    index: CornersIndex,
    orientation: CornerOrientation,
}

impl Corner {
    pub fn get_index(&self) -> CornersIndex {
        self.index
    }

    pub fn get_orientation(&self) -> CornerOrientation {
        self.orientation
    }
}

impl Default for Corner {
    fn default() -> Self {
        Corner {
            index: ULB,
            orientation: CornerOrientation::Oriented,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Sequence, TryFromPrimitive)]
#[repr(u8)]
pub enum CornersIndex {
    ULB = 0,
    URB = 1,
    URF = 2,
    ULF = 3,
    DLF = 4,
    DLB = 5,
    DRB = 6,
    DRF = 7,
}

#[derive(Clone, Eq, Hash, PartialEq)]
pub struct CubieCube {
    corners: [Corner; NUM_CORNERS],
    edges: [Edge; NUM_EDGES],
    centers: [CubeFace; NUM_SIDES],
}

impl Default for CubieCube {
    fn default() -> Self {
        let mut corners: [Corner; NUM_CORNERS] = [Corner::default(); NUM_CORNERS];
        for (idx, corner_idx) in all::<CornersIndex>().enumerate() {
            corners[idx].index = corner_idx;
        }
        let mut edges: [Edge; NUM_EDGES] = [Edge::default(); NUM_EDGES];
        for (idx, edge_idx) in all::<EdgesIndex>().enumerate() {
            edges[idx].index = edge_idx;
        }
        let mut centers: [CubeFace; NUM_SIDES] = [Up; NUM_SIDES];
        for (idx, face) in all::<CubeFace>().enumerate() {
            centers[idx] = face;
        }

        Self {
            corners,
            edges,
            centers,
        }
    }
}

impl CubieCube {
    pub fn get_edge(&self, edge_index: EdgesIndex) -> Edge {
        self.edges[edge_index as usize]
    }

    pub fn get_corner(&self, corner_index: CornersIndex) -> Corner {
        self.corners[corner_index as usize]
    }

    pub fn to_array_cube(&self) -> ArrayCube {
        let mut colors = [[[White; CUBE_SIZE]; CUBE_SIZE]; NUM_SIDES];
        for face in all::<CubeFace>() {
            for row in 0..CUBE_SIZE {
                for col in 0..CUBE_SIZE {
                    colors[face as usize][row][col] = self.get_facelet_color(face, row, col);
                }
            }
        }
        ArrayCube::new(colors)
    }

    fn get_facelet_color(&self, face_idx: CubeFace, row_idx: usize, col_idx: usize) -> Color {
        match (face_idx, row_idx, col_idx) {
            (Up, 0, 0) => self.get_corner_colors(ULB)[0],
            (Up, 0, 1) => self.get_edge_colors(UB)[0],
            (Up, 0, 2) => self.get_corner_colors(URB)[0],
            (Up, 1, 0) => self.get_edge_colors(UL)[0],
            (Up, 1, 1) => self.get_center_color(Up),
            (Up, 1, 2) => self.get_edge_colors(UR)[0],
            (Up, 2, 0) => self.get_corner_colors(ULF)[0],
            (Up, 2, 1) => self.get_edge_colors(UF)[0],
            (Up, 2, 2) => self.get_corner_colors(URF)[0],
            (Left, 0, 0) => self.get_corner_colors(ULB)[1],
            (Left, 0, 1) => self.get_edge_colors(UL)[1],
            (Left, 0, 2) => self.get_corner_colors(ULF)[2],
            (Left, 1, 0) => self.get_edge_colors(BL)[1],
            (Left, 1, 1) => self.get_center_color(Left),
            (Left, 1, 2) => self.get_edge_colors(FL)[1],
            (Left, 2, 0) => self.get_corner_colors(DLB)[2],
            (Left, 2, 1) => self.get_edge_colors(DL)[1],
            (Left, 2, 2) => self.get_corner_colors(DLF)[1],
            (Front, 0, 0) => self.get_corner_colors(ULF)[1],
            (Front, 0, 1) => self.get_edge_colors(UF)[1],
            (Front, 0, 2) => self.get_corner_colors(URF)[2],
            (Front, 1, 0) => self.get_edge_colors(FL)[0],
            (Front, 1, 1) => self.get_center_color(Front),
            (Front, 1, 2) => self.get_edge_colors(FR)[0],
            (Front, 2, 0) => self.get_corner_colors(DLF)[2],
            (Front, 2, 1) => self.get_edge_colors(DF)[1],
            (Front, 2, 2) => self.get_corner_colors(DRF)[1],
            (Right, 0, 0) => self.get_corner_colors(URF)[1],
            (Right, 0, 1) => self.get_edge_colors(UR)[1],
            (Right, 0, 2) => self.get_corner_colors(URB)[2],
            (Right, 1, 0) => self.get_edge_colors(FR)[1],
            (Right, 1, 1) => self.get_center_color(Right),
            (Right, 1, 2) => self.get_edge_colors(BR)[1],
            (Right, 2, 0) => self.get_corner_colors(DRF)[2],
            (Right, 2, 1) => self.get_edge_colors(DR)[1],
            (Right, 2, 2) => self.get_corner_colors(DRB)[1],
            (Back, 0, 0) => self.get_corner_colors(URB)[1],
            (Back, 0, 1) => self.get_edge_colors(UB)[1],
            (Back, 0, 2) => self.get_corner_colors(ULB)[2],
            (Back, 1, 0) => self.get_edge_colors(BR)[0],
            (Back, 1, 1) => self.get_center_color(Back),
            (Back, 1, 2) => self.get_edge_colors(BL)[0],
            (Back, 2, 0) => self.get_corner_colors(DRB)[2],
            (Back, 2, 1) => self.get_edge_colors(DB)[1],
            (Back, 2, 2) => self.get_corner_colors(DLB)[1],
            (Down, 0, 0) => self.get_corner_colors(DLF)[0],
            (Down, 0, 1) => self.get_edge_colors(DF)[0],
            (Down, 0, 2) => self.get_corner_colors(DRF)[0],
            (Down, 1, 0) => self.get_edge_colors(DL)[0],
            (Down, 1, 1) => self.get_center_color(Down),
            (Down, 1, 2) => self.get_edge_colors(DR)[0],
            (Down, 2, 0) => self.get_corner_colors(DLB)[0],
            (Down, 2, 1) => self.get_edge_colors(DB)[0],
            (Down, 2, 2) => self.get_corner_colors(DRB)[0],
            _ => panic!("Uncovered facelet index on cube!"),
        }
    }

    fn get_face(&self, face_idx: CubeFace) -> [[Color; CUBE_SIZE]; CUBE_SIZE] {
        let mut face = [[White; CUBE_SIZE]; CUBE_SIZE];
        for row in 0..CUBE_SIZE {
            for col in 0..CUBE_SIZE {
                face[row][col] = self.get_facelet_color(face_idx, row, col);
            }
        }

        face
    }

    pub fn get_corner_colors(&self, corner_idx: CornersIndex) -> [Color; 3] {
        let corner = self.corners[corner_idx as usize];
        let mut colors = match corner.index as CornersIndex {
            ULB => [White, Green, Orange],
            URB => [White, Orange, Blue],
            URF => [White, Blue, Red],
            ULF => [White, Red, Green],
            DLF => [Yellow, Green, Red],
            DLB => [Yellow, Orange, Green],
            DRB => [Yellow, Blue, Orange],
            DRF => [Yellow, Red, Blue],
        };

        colors.rotate_left(corner.orientation as usize);

        colors
    }

    fn get_edge_colors(&self, edge_idx: EdgesIndex) -> [Color; 2] {
        let edge = self.edges[edge_idx as usize];
        let mut colors = match edge.index {
            UB => [White, Orange],
            UR => [White, Blue],
            UF => [White, Red],
            UL => [White, Green],
            FR => [Red, Blue],
            FL => [Red, Green],
            BL => [Orange, Green],
            BR => [Orange, Blue],
            DF => [Yellow, Red],
            DL => [Yellow, Green],
            DB => [Yellow, Orange],
            DR => [Yellow, Blue],
        };

        colors.rotate_left(edge.orientation as usize);

        colors
    }
    fn get_center_color(&self, center_idx: CubeFace) -> Color {
        let face = self.centers[center_idx as usize];
        match face {
            Up => White,
            Left => Green,
            Front => Red,
            Right => Blue,
            Back => Orange,
            Down => Yellow,
        }
    }

    fn update_corner_orientation(&mut self, corner_idx: CornersIndex, amount: u8) {
        let idx = corner_idx as usize;
        let updated = (self.corners[idx].orientation as u8 + amount) % 3;
        self.corners[idx].orientation =
            CornerOrientation::try_from(updated).expect("Invalid orientation specified.")
    }

    fn update_edge_orientatiopn(&mut self, edge_idx: EdgesIndex) {
        let idx = edge_idx as usize;
        let updated = if self.edges[idx].orientation == EdgeOrientation::Flipped {
            EdgeOrientation::Oriented
        } else {
            EdgeOrientation::Flipped
        };
        self.edges[idx].orientation = updated;
    }

    fn f_move(&mut self) {
        let tmp = self.corners[ULF as usize];
        self.corners[ULF as usize] = self.corners[DLF as usize];
        self.corners[DLF as usize] = self.corners[DRF as usize];
        self.corners[DRF as usize] = self.corners[URF as usize];
        self.corners[URF as usize] = tmp;

        let tmp = self.edges[UF as usize];
        self.edges[UF as usize] = self.edges[FL as usize];
        self.edges[FL as usize] = self.edges[DF as usize];
        self.edges[DF as usize] = self.edges[FR as usize];
        self.edges[FR as usize] = tmp;

        self.update_corner_orientation(ULF, 1);
        self.update_corner_orientation(URF, 2);
        self.update_corner_orientation(DRF, 1);
        self.update_corner_orientation(DLF, 2);

        self.update_edge_orientatiopn(UF);
        self.update_edge_orientatiopn(FL);
        self.update_edge_orientatiopn(DF);
        self.update_edge_orientatiopn(FR);
    }

    fn f_prime_move(&mut self) {
        let tmp = self.corners[ULF as usize];
        self.corners[ULF as usize] = self.corners[URF as usize];
        self.corners[URF as usize] = self.corners[DRF as usize];
        self.corners[DRF as usize] = self.corners[DLF as usize];
        self.corners[DLF as usize] = tmp;

        let tmp = self.edges[UF as usize];
        self.edges[UF as usize] = self.edges[FR as usize];
        self.edges[FR as usize] = self.edges[DF as usize];
        self.edges[DF as usize] = self.edges[FL as usize];
        self.edges[FL as usize] = tmp;

        self.update_corner_orientation(ULF, 1);
        self.update_corner_orientation(URF, 2);
        self.update_corner_orientation(DRF, 1);
        self.update_corner_orientation(DLF, 2);

        self.update_edge_orientatiopn(UF);
        self.update_edge_orientatiopn(FL);
        self.update_edge_orientatiopn(DF);
        self.update_edge_orientatiopn(FR);
    }

    fn b_move(&mut self) {
        let tmp = self.corners[ULB as usize];
        self.corners[ULB as usize] = self.corners[URB as usize];
        self.corners[URB as usize] = self.corners[DRB as usize];
        self.corners[DRB as usize] = self.corners[DLB as usize];
        self.corners[DLB as usize] = tmp;

        let tmp = self.edges[UB as usize];
        self.edges[UB as usize] = self.edges[BR as usize];
        self.edges[BR as usize] = self.edges[DB as usize];
        self.edges[DB as usize] = self.edges[BL as usize];
        self.edges[BL as usize] = tmp;

        self.update_corner_orientation(ULB, 2);
        self.update_corner_orientation(URB, 1);
        self.update_corner_orientation(DRB, 2);
        self.update_corner_orientation(DLB, 1);

        self.update_edge_orientatiopn(UB);
        self.update_edge_orientatiopn(BL);
        self.update_edge_orientatiopn(DB);
        self.update_edge_orientatiopn(BR);
    }

    fn b_prime_move(&mut self) {
        let tmp = self.corners[ULB as usize];
        self.corners[ULB as usize] = self.corners[DLB as usize];
        self.corners[DLB as usize] = self.corners[DRB as usize];
        self.corners[DRB as usize] = self.corners[URB as usize];
        self.corners[URB as usize] = tmp;

        let tmp = self.edges[UB as usize];
        self.edges[UB as usize] = self.edges[BL as usize];
        self.edges[BL as usize] = self.edges[DB as usize];
        self.edges[DB as usize] = self.edges[BR as usize];
        self.edges[BR as usize] = tmp;

        self.update_corner_orientation(ULB, 2);
        self.update_corner_orientation(URB, 1);
        self.update_corner_orientation(DRB, 2);
        self.update_corner_orientation(DLB, 1);

        self.update_edge_orientatiopn(UB);
        self.update_edge_orientatiopn(BL);
        self.update_edge_orientatiopn(DB);
        self.update_edge_orientatiopn(BR);
    }

    fn l_move(&mut self) {
        let tmp = self.corners[ULB as usize];
        self.corners[ULB as usize] = self.corners[DLB as usize];
        self.corners[DLB as usize] = self.corners[DLF as usize];
        self.corners[DLF as usize] = self.corners[ULF as usize];
        self.corners[ULF as usize] = tmp;

        let tmp = self.edges[UL as usize];
        self.edges[UL as usize] = self.edges[BL as usize];
        self.edges[BL as usize] = self.edges[DL as usize];
        self.edges[DL as usize] = self.edges[FL as usize];
        self.edges[FL as usize] = tmp;

        self.update_corner_orientation(ULB, 1);
        self.update_corner_orientation(ULF, 2);
        self.update_corner_orientation(DLF, 1);
        self.update_corner_orientation(DLB, 2);
    }

    fn l_prime_move(&mut self) {
        let tmp = self.corners[ULB as usize];
        self.corners[ULB as usize] = self.corners[ULF as usize];
        self.corners[ULF as usize] = self.corners[DLF as usize];
        self.corners[DLF as usize] = self.corners[DLB as usize];
        self.corners[DLB as usize] = tmp;

        let tmp = self.edges[UL as usize];
        self.edges[UL as usize] = self.edges[FL as usize];
        self.edges[FL as usize] = self.edges[DL as usize];
        self.edges[DL as usize] = self.edges[BL as usize];
        self.edges[BL as usize] = tmp;

        self.update_corner_orientation(ULB, 1);
        self.update_corner_orientation(ULF, 2);
        self.update_corner_orientation(DLF, 1);
        self.update_corner_orientation(DLB, 2);
    }

    fn r_move(&mut self) {
        let tmp = self.corners[URB as usize];
        self.corners[URB as usize] = self.corners[URF as usize];
        self.corners[URF as usize] = self.corners[DRF as usize];
        self.corners[DRF as usize] = self.corners[DRB as usize];
        self.corners[DRB as usize] = tmp;

        let tmp = self.edges[UR as usize];
        self.edges[UR as usize] = self.edges[FR as usize];
        self.edges[FR as usize] = self.edges[DR as usize];
        self.edges[DR as usize] = self.edges[BR as usize];
        self.edges[BR as usize] = tmp;

        self.update_corner_orientation(URB, 2);
        self.update_corner_orientation(URF, 1);
        self.update_corner_orientation(DRF, 2);
        self.update_corner_orientation(DRB, 1);
    }

    fn r_prime_move(&mut self) {
        let tmp = self.corners[URB as usize];
        self.corners[URB as usize] = self.corners[DRB as usize];
        self.corners[DRB as usize] = self.corners[DRF as usize];
        self.corners[DRF as usize] = self.corners[URF as usize];
        self.corners[URF as usize] = tmp;

        let tmp = self.edges[UR as usize];
        self.edges[UR as usize] = self.edges[BR as usize];
        self.edges[BR as usize] = self.edges[DR as usize];
        self.edges[DR as usize] = self.edges[FR as usize];
        self.edges[FR as usize] = tmp;

        self.update_corner_orientation(URB, 2);
        self.update_corner_orientation(URF, 1);
        self.update_corner_orientation(DRF, 2);
        self.update_corner_orientation(DRB, 1);
    }

    fn u_move(&mut self) {
        let tmp = self.corners[URB as usize];
        self.corners[URB as usize] = self.corners[ULB as usize];
        self.corners[ULB as usize] = self.corners[ULF as usize];
        self.corners[ULF as usize] = self.corners[URF as usize];
        self.corners[URF as usize] = tmp;

        let tmp = self.edges[UR as usize];
        self.edges[UR as usize] = self.edges[UB as usize];
        self.edges[UB as usize] = self.edges[UL as usize];
        self.edges[UL as usize] = self.edges[UF as usize];
        self.edges[UF as usize] = tmp;
    }

    fn u_prime_move(&mut self) {
        let tmp = self.corners[URB as usize];
        self.corners[URB as usize] = self.corners[URF as usize];
        self.corners[URF as usize] = self.corners[ULF as usize];
        self.corners[ULF as usize] = self.corners[ULB as usize];
        self.corners[ULB as usize] = tmp;

        let tmp = self.edges[UR as usize];
        self.edges[UR as usize] = self.edges[UF as usize];
        self.edges[UF as usize] = self.edges[UL as usize];
        self.edges[UL as usize] = self.edges[UB as usize];
        self.edges[UB as usize] = tmp;
    }

    fn d_move(&mut self) {
        let tmp = self.corners[DRB as usize];
        self.corners[DRB as usize] = self.corners[DRF as usize];
        self.corners[DRF as usize] = self.corners[DLF as usize];
        self.corners[DLF as usize] = self.corners[DLB as usize];
        self.corners[DLB as usize] = tmp;

        let tmp = self.edges[DR as usize];
        self.edges[DR as usize] = self.edges[DF as usize];
        self.edges[DF as usize] = self.edges[DL as usize];
        self.edges[DL as usize] = self.edges[DB as usize];
        self.edges[DB as usize] = tmp;
    }

    fn d_prime_move(&mut self) {
        let tmp = self.corners[DRB as usize];
        self.corners[DRB as usize] = self.corners[DLB as usize];
        self.corners[DLB as usize] = self.corners[DLF as usize];
        self.corners[DLF as usize] = self.corners[DRF as usize];
        self.corners[DRF as usize] = tmp;

        let tmp = self.edges[DR as usize];
        self.edges[DR as usize] = self.edges[DB as usize];
        self.edges[DB as usize] = self.edges[DL as usize];
        self.edges[DL as usize] = self.edges[DF as usize];
        self.edges[DF as usize] = tmp;
    }
}

impl Cube for CubieCube {
    fn is_solved(&self) -> bool {
        for (idx, edge) in self.edges.iter().enumerate() {
            if edge.orientation != EdgeOrientation::Oriented
                || edge.index != EdgesIndex::try_from(idx as u8).expect("Unexpected edge index!")
            {
                return false;
            }
        }

        for (idx, corner) in self.corners.iter().enumerate() {
            if corner.orientation != CornerOrientation::Oriented
                || corner.index
                    != CornersIndex::try_from(idx as u8).expect("Unexpected corner index!")
            {
                return false;
            }
        }

        return true;
    }

    fn cube_move(&mut self, cube_move: CubeMove) {
        match (cube_move.face, cube_move.direction) {
            (Front, Clockwise) => self.f_move(),
            (Front, Counterclockwise) => self.f_prime_move(),
            (Back, Clockwise) => self.b_move(),
            (Back, Counterclockwise) => self.b_prime_move(),
            (Left, Clockwise) => self.l_move(),
            (Left, Counterclockwise) => self.l_prime_move(),
            (Right, Clockwise) => self.r_move(),
            (Right, Counterclockwise) => self.r_prime_move(),
            (Up, Clockwise) => self.u_move(),
            (Up, Counterclockwise) => self.u_prime_move(),
            (Down, Clockwise) => self.d_move(),
            (Down, Counterclockwise) => self.d_prime_move(),
        }
    }

    fn solve(&self) -> Vec<CubeMove> {
        todo!()
    }
}

#[cfg(test)]
mod test {
    use crate::array_cube::ArrayCube;
    use crate::cube::Cube;
    use crate::cube::CubeFace::*;
    use crate::cube::Direction::*;
    use crate::cube::{CubeFace, CubeMove, Direction};
    use crate::cubie_cube::CubieCube;
    use quickcheck_macros::quickcheck;
    use rstest::rstest;

    #[test]
    fn default_cube_matches_default_arraycube() {
        let array_cube = ArrayCube::default();
        let cubie_cube = CubieCube::default();

        assert_eq!(array_cube, cubie_cube.to_array_cube());
    }

    #[rstest]
    #[case(Front, Clockwise)]
    #[case(Front, Counterclockwise)]
    #[case(Left, Clockwise)]
    #[case(Left, Counterclockwise)]
    #[case(Right, Clockwise)]
    #[case(Right, Counterclockwise)]
    #[case(Back, Clockwise)]
    #[case(Back, Counterclockwise)]
    #[case(Up, Clockwise)]
    #[case(Up, Counterclockwise)]
    #[case(Down, Clockwise)]
    #[case(Down, Counterclockwise)]
    fn rotate_face(#[case] face: CubeFace, #[case] direction: Direction) {
        let cube_move = CubeMove::new(face, direction);
        let mut array_cube = ArrayCube::default();
        array_cube.cube_move(cube_move);
        let mut cubie_cube = CubieCube::default();
        cubie_cube.cube_move(cube_move);

        assert_eq!(array_cube, cubie_cube.to_array_cube());
    }

    #[test]
    fn double_rotate() {
        let mut array_cube = ArrayCube::default();
        array_cube.cube_move(CubeMove::new(Front, Clockwise));
        array_cube.cube_move(CubeMove::new(Left, Clockwise));

        let mut cubie_cube = CubieCube::default();
        cubie_cube.cube_move(CubeMove::new(Front, Clockwise));
        cubie_cube.cube_move(CubeMove::new(Left, Clockwise));

        assert_eq!(array_cube, cubie_cube.to_array_cube());
    }

    #[quickcheck]
    fn randomized_cubes_match(rotations: Vec<CubeMove>) -> bool {
        let mut array_cube = ArrayCube::default();
        let mut cubie_cube = CubieCube::default();

        for mv in rotations {
            array_cube.cube_move(mv);
            cubie_cube.cube_move(mv)
        }

        let result = cubie_cube.to_array_cube();

        result == array_cube
    }
}
