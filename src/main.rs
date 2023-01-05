use itertools::{Either, Itertools, rev};
use crate::Rotation::Clockwise;

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

/*
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
/// A 'sticker' for lack of a better term is one square on the face
/// The grid is organized such that 0, 0 = bottom, left.
///
/// We use a 2-d array so that we can support different sizes of cubes later.
#[derive(Clone, Debug, PartialEq)]
struct Face {
    stickers: [[Color; 3]; 3]
}

impl Face {
    pub fn rotate(mut self, rotation: Rotation) -> Self {
        let n: usize = self.stickers.len();
        let x: usize = n / 2;
        let y: usize = n - 1;
        for i in 0..x {
            for j in i..y - i {
                match rotation {
                    Clockwise => {
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

        self
    }
}

/// The cube has 6 faces. They are organized in the array such that opposite faces are separated
/// by 3. In pictures:
///    T         0
///  L F R Ba  1 2 4 5    where T = top, L = left, F = front, Ba = back, Bo = bottom
///    Bo        3
/// TODO: Figure out if this works for shapes with more than 6 sides
struct Cube {
    faces: [Face; 6]
    // bottom: Face,
    // top: Face,
    // front: Face,
    // back: Face,
    // left: Face,
    // right: Face
}

/// This is used to "look up" the surrounding faces when rotating a face.
const surrounding_faces: [[usize; 4]; 6] = [
    [1,2,4,5],
    [0,2,3,5],
    [0,4,3,1],
    [5,4,2,1],
    [1,3,4,0],
    [1,3,4,0],
];

fn main() {
    let val = 8;
    let color: Color = unsafe { std::mem::transmute(val as u8) };
    println!("{:?}", color)
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck::{Arbitrary, Gen, quickcheck};
    use arr_macro::arr;
    use crate::Rotation::Counterclockwise;
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

    // #[quickcheck_macros::quickcheck]
    // fn rotate_face_clockwise(face: Face) -> bool {
    //     let expected = face.clone();
    //     let result = face
    //         .rotate(Clockwise)
    //         .rotate(Clockwise)
    //         .rotate(Clockwise)
    //         .rotate(Clockwise);
    //     expected == result
    // }
    #[quickcheck_macros::quickcheck]
    fn rotate_face_clockwise_90_degrees(face: Face) -> bool {
        let expected = face.clone();
        let result = face
            .rotate(Clockwise);

        // Verify that iterating the rotated face vertically from top right == iterating the
        // original face horizontally from top left
        for (row_idx, row) in expected.stickers.iter().enumerate() {
            for (col_idx, color) in row.iter().enumerate() {
                if *color != result.stickers[col_idx][result.stickers.len() - row_idx - 1] {
                    return false
                }
            }
        }

        return true
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_face_counterclockwise_90_degrees(face: Face) -> bool {
        let expected = face.clone();
        let result = face
            .rotate(Counterclockwise);

        // Verify that iterating the rotated face vertically from bottom left == iterating the
        // original face horizontally from top left
        for (row_idx, row) in expected.stickers.iter().enumerate() {
            for (col_idx, color) in row.iter().enumerate() {
                if *color != result.stickers[result.stickers[0].len() - col_idx - 1][row_idx] {
                    return false
                }
            }
        }

        return true
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_face_clockwise_180_degrees(face: Face) -> bool {
        let expected = face.clone();
        let result = face
            .rotate(Clockwise)
            .rotate(Clockwise);

        // Verify that iterating the rotated face from bottom right == iterating the original face
        // from top left
        for (row_idx, row) in expected.stickers.iter().enumerate() {
            for (col_idx, color) in row.iter().enumerate() {
                if *color != result.stickers[result.stickers.len() - row_idx - 1][result.stickers[0].len() - col_idx - 1] {
                    return false
                }
            }
        }

        return true
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_face_180_degrees_clockwise_equals_counterclockwise(face: Face) -> bool {
        let clockwise = face.clone()
            .rotate(Clockwise)
            .rotate(Clockwise);

        let counterclockwise = face
            .rotate(Counterclockwise)
            .rotate(Counterclockwise);

        clockwise == counterclockwise
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_face_clockwise_four_times_is_identity(face: Face) -> bool {
        let expected = face.clone();
        let result = face
            .rotate(Clockwise)
            .rotate(Clockwise)
            .rotate(Clockwise)
            .rotate(Clockwise);
        expected == result
    }

    #[quickcheck_macros::quickcheck]
    fn rotate_face_counterclockwise_four_times_is_identity(face: Face) -> bool {
        let expected = face.clone();
        let result = face
            .rotate(Counterclockwise)
            .rotate(Counterclockwise)
            .rotate(Counterclockwise)
            .rotate(Counterclockwise);
        expected == result
    }

    // #[rstest]
    // #[case::transmit(Command::Transmit, "hello world!!11!".as_bytes())]
    // #[case::ack(Command::Ack, & [])]
    // #[case::nack(Command::Nack, & [])]
    // #[case::start(Command::Start, & [])]
    // #[case::stop(Command::Stop, & [])]
    // #[case::reset(Command::Reset, & [])]
    // fn to_from_bytes_valid_input_returns_ok(#[case] cmd: Command, #[case] data: &[u8]) {
    //     // arrange
    //     let expected = Packet::new(cmd, data);
    //     let bytes = expected.iter();
    //     static BUFF: [u8; 1024] = [0; 1024];
    //
    //     // act
    //     let (header, bytes) = Header::from_iter(bytes).unwrap();
    //     let result = header.to_packet(bytes, &mut BUFF[0..size]);
    //
    //     // assert
    //     assert_eq!(result, Ok(expected));
    // }
}