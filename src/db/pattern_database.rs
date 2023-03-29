use crate::cube::Cube;
use crate::cubie_cube::CubieCube;
use derivative::Derivative;
use std::error::Error;
use std::fmt::{Debug, Formatter};
use std::fs::File;
use std::io::{ErrorKind, Read, Write};
const DEFAULT_VAL: u8 = 0xFF;

#[derive(Derivative)]
#[derivative(Debug)]
struct PatternDatabase<const N: usize> {
    db: [u8; N],
    num_items: usize,
    #[derivative(Debug = "ignore")]
    get_index: fn(&CubieCube) -> usize,
}

impl<const N: usize> PatternDatabase<N> {
    fn new(get_index: fn(&CubieCube) -> usize) -> Self {
        Self {
            db: [DEFAULT_VAL; N],
            num_items: 0,
            get_index,
        }
    }

    fn set_num_moves(&mut self, cube: &CubieCube, num_moves: u8) -> bool {
        let idx = (self.get_index)(cube);
        self.set_num_moves_for_index(idx, num_moves)
    }

    /// Set a value for a given index in the DB
    /// Returns false if the value was previously set. True otherwise.
    fn set_num_moves_for_index(&mut self, index: usize, num_moves: u8) -> bool {
        if self.db[index] != DEFAULT_VAL {
            false
        } else {
            self.db[index] = num_moves;
            self.num_items += 1;
            true
        }
    }

    fn get_num_moves(&self, cube: &CubieCube) -> u8 {
        let idx = (self.get_index)(cube);

        self.get_num_moves_for_index(idx)
    }

    fn get_num_moves_for_index(&self, index: usize) -> u8 {
        self.db[index]
    }

    fn get_size(&self) -> usize {
        self.db.len()
    }

    fn is_full(&self) -> bool {
        self.num_items == self.db.len()
    }

    fn get_num_items(&self) -> usize {
        self.num_items
    }

    fn to_file(&self, file_path: &str) -> Result<(), std::io::Error> {
        let mut file = File::create(file_path)?;
        file.write_all(N.to_le_bytes().as_slice())?;
        file.write_all(self.db.as_slice())?;
        file.write_all(self.num_items.to_le_bytes().as_slice())
    }

    fn from_file(
        file_path: &str,
        get_index: fn(&CubieCube) -> usize,
    ) -> Result<Self, std::io::Error> {
        File::open(file_path).and_then(|mut file| {
            let mut size_bytes = [0; 8];
            file.read(&mut size_bytes)?;
            let size = usize::from_le_bytes(size_bytes);
            if size != N {
                return Err(std::io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Size in file ({} bytes) did not match expected size of database ({} bytes)", size, N))
                )
            }
            let mut data = [0xFF; N];
            file.read(&mut data)?;
            let mut num_items_bytes = [0; 8];
            file.read(&mut num_items_bytes)?;
            let num_items = usize::from_le_bytes(size_bytes);
            Ok(Self {
                db: data,
                num_items,
                get_index,
            })
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::cube::{Cube, U_PRIME_MOVE};
    use crate::cube::{CubeMove, L_MOVE, L_PRIME_MOVE, R_MOVE, R_PRIME_MOVE, U_MOVE};
    use crate::cubie_cube::CornersIndex::ULB;
    use crate::cubie_cube::{CornersIndex, CubieCube};
    use crate::db::pattern_database::PatternDatabase;
    use enum_iterator::{all, cardinality};
    use rand::prelude::SliceRandom;
    use rand::thread_rng;
    use tempfile::{tempfile, NamedTempFile};

    // A Simple DB with 8 elements (since the  corner can be in 8 different positions)
    fn get_ulb_index(cube: &CubieCube) -> usize {
        cube.get_corner(CornersIndex::ULB).get_index() as usize
    }

    fn cube_moves(cube: &CubieCube, moves: Vec<CubeMove>) -> CubieCube {
        let mut result = cube.clone();
        for mv in moves {
            result.cube_move(mv);
        }
        result
    }

    #[test]
    fn smoke_test() {
        let db = PatternDatabase::<8>::new(get_ulb_index);
        assert_eq!(8, db.get_size());
        assert_eq!(0, db.get_num_items());
    }

    #[test]
    fn get_set_by_index_works() {
        let mut db = PatternDatabase::<8>::new(get_ulb_index);
        let mut vals: Vec<u8> = (0u8..db.get_size() as u8).collect();
        let mut rng = thread_rng();
        vals.shuffle(&mut rng);
        for i in 0..db.get_size() {
            assert!(db.set_num_moves_for_index(i, vals[i]));
        }

        assert_eq!(8, db.get_size());
        assert_eq!(8, db.get_num_items());

        for i in 0..db.get_size() {
            assert_eq!(vals[i], db.get_num_moves_for_index(i));
        }
    }

    #[test]
    fn set_index_twice_returns_false() {
        let mut db = PatternDatabase::<8>::new(get_ulb_index);
        assert!(db.set_num_moves_for_index(0, 5));
        assert!(!db.set_num_moves_for_index(0, 7));
        assert_eq!(5, db.get_num_moves_for_index(0))
    }

    #[test]
    fn get_set_by_cube_works() {
        let mut db = PatternDatabase::<8>::new(get_ulb_index);
        let cube = CubieCube::default();
        // There are infinite ways to generate this list. We're just rotating all of the corners
        // through the 'ULB' position to generate all of the indexes.
        let moves = [
            vec![],
            vec![L_MOVE],
            vec![L_MOVE, L_MOVE],
            vec![L_PRIME_MOVE],
            vec![U_PRIME_MOVE],
            vec![U_MOVE, U_MOVE],
            vec![R_MOVE, U_MOVE, U_MOVE],
            vec![R_PRIME_MOVE, U_PRIME_MOVE],
        ];
        let positions = moves.clone().map(|mvs| cube_moves(&cube, mvs));

        for i in 0..db.get_size() {
            assert!(db.set_num_moves(&positions[i], moves[i].len() as u8));
        }

        assert_eq!(8, db.get_size());
        assert_eq!(8, db.get_num_items());

        for i in 0..db.get_size() {
            assert_eq!(moves[i].len() as u8, db.get_num_moves(&positions[i]));
        }
    }

    #[test]
    fn to_from_file_works() {
        let f = NamedTempFile::new().unwrap();

        let mut db = PatternDatabase::<8>::new(get_ulb_index);
        let mut vals: Vec<u8> = (0u8..db.get_size() as u8).collect();
        let mut rng = thread_rng();
        vals.shuffle(&mut rng);
        for i in 0..db.get_size() {
            assert!(db.set_num_moves_for_index(i, vals[i]));
        }

        let write_result = db.to_file(f.path().to_str().unwrap());
        assert!(write_result.is_ok());
        let read_result =
            PatternDatabase::<8>::from_file(f.path().to_str().unwrap(), get_ulb_index);
        assert!(read_result.is_ok());
        let read_db = read_result.unwrap();

        assert_eq!(db.get_size(), read_db.get_size());
        assert_eq!(db.get_num_items(), read_db.get_num_items());
        assert_eq!(db.db, read_db.db);
    }
}
