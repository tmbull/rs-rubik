pub struct BasicBitSet<const N: usize> {
    bits: u64,
    mask: u64,
}

impl<const N: usize> BasicBitSet<N> {
    pub fn new() -> Self {
        Self::from(0)
    }

    pub fn from(val: u64) -> Self {
        if N > 64 || N < 1 {
            panic!("Size of BasicBitSet (N) must be 64 or fewer bits")
        }

        let mask = if N == 64 { u64::MAX } else { (1 << N) - 1 };

        BasicBitSet { bits: val, mask }
    }

    pub fn get_val(&self) -> u64 {
        self.bits & self.mask
    }

    pub fn set_bit(&mut self, bit: usize) {
        if N <= bit {
            panic!("Bit argument must be less than the size of the bit set")
        }
        self.bits |= 1 << bit;
    }

    pub fn unset_bit(&mut self, bit: usize) {
        if N <= bit {
            panic!("Bit argument must be less than the size of the bit set")
        }
        self.bits &= !(1 << bit);
    }
}

#[cfg(test)]
mod tests {
    use crate::db::basic_bitset::BasicBitSet;
    use quickcheck_macros::quickcheck;

    #[quickcheck]
    fn from_random_works(val: u64) -> bool {
        let result = BasicBitSet::<64>::from(val);

        result.get_val() == val && result.mask.count_ones() == 64
    }

    fn my_args() -> Vec<u64> {
        (0..10).collect::<Vec<_>>()
    }

    #[test]
    #[should_panic]
    fn size_less_than_1_panics() {
        let result = BasicBitSet::<0>::new();
    }

    #[test]
    #[should_panic]
    fn size_greater_than_64_panics() {
        let result = BasicBitSet::<65>::new();
    }

    #[test]
    fn mask_bits_equal_size() {
        let result = BasicBitSet::<10>::new();

        assert_eq!(result.mask.count_ones(), 10)
    }

    #[test]
    fn mask_works_as_expected() {
        let mut result = BasicBitSet::<2>::from(u64::MAX);

        assert_eq!(result.get_val(), 3);
    }

    #[test]
    fn mask_works_as_expected_max_u64() {
        let mut result = BasicBitSet::<64>::from(u64::MAX);

        assert_eq!(result.get_val(), u64::MAX);
    }

    #[test]
    #[should_panic]
    fn set_bit_greater_than_size_panics() {
        let mut result = BasicBitSet::<10>::new();
        result.set_bit(11);
    }

    #[test]
    fn set_bit_works() {
        let mut result = BasicBitSet::<10>::from(3);
        result.set_bit(5);

        assert_eq!(result.get_val(), 35)
    }

    #[test]
    fn set_msb_works() {
        let mut result = BasicBitSet::<2>::new();
        result.set_bit(1);

        assert_eq!(result.get_val(), 2)
    }

    #[test]
    fn set_all_bits_works() {
        let mut result = BasicBitSet::<2>::new();
        result.set_bit(0);
        result.set_bit(1);

        assert_eq!(result.get_val(), 3);
    }

    #[test]
    #[should_panic]
    fn unset_bit_greater_than_size_panics() {
        let mut result = BasicBitSet::<10>::new();
        result.unset_bit(11);
    }

    #[test]
    fn unset_bit_works() {
        let mut result = BasicBitSet::<10>::from(35);
        result.unset_bit(5);

        assert_eq!(result.get_val(), 3)
    }

    #[test]
    fn unset_msb_works() {
        let mut result = BasicBitSet::<2>::from(3);
        result.unset_bit(1);

        assert_eq!(result.get_val(), 1)
    }

    #[test]
    fn unset_all_bits_works() {
        let mut result = BasicBitSet::<2>::from(3);
        result.unset_bit(0);
        result.unset_bit(1);

        assert_eq!(result.get_val(), 0);
    }
}
