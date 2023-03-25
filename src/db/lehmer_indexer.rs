//! This module was translated directly from Ben Botto's PermutationIndexer C++ class. See his
//! excellent blog post [here](https://medium.com/@benjamin.botto/sequentially-indexing-permutations-a-linear-algorithm-for-computing-lexicographic-rank-a22220ffd6e3)
//!
use crate::db::basic_bitset::BasicBitSet;
fn factorial(n: u64) -> u64 {
    if n <= 1 {
        1
    } else {
        n * factorial(n - 1)
    }
}

// Calculate nPk: n!/(n-k)!.
fn pick(n: u64, k: u64) -> u64 {
    factorial(n) / factorial(n - k)
}

struct PermutationIndexer<const N: usize, const K: usize = N>
where
    [(); ({ 1 << N }) - 1]: Sized,
{
    ones_count_lookup: [usize; ({ 1 << N }) - 1],

    // Precomputed table of factorials (or "picks" if N != K).  They're in
    // reverse order.
    factorials: [u64; K],
}

impl<const N: usize, const K: usize> PermutationIndexer<N, K>
where
    [(); ({ 1 << N }) - 1]: Sized,
{
    fn new() -> Self {
        let mut indexer = PermutationIndexer {
            ones_count_lookup: [0; (1 << N) - 1],
            factorials: [0; K],
        };
        for i in 0..((1 << N) - 1) {
            indexer.ones_count_lookup[i] = i.count_ones() as usize;
        }

        for i_size in 0..K {
            let i = i_size as u64;
            indexer.factorials[i_size] = pick(N as u64 - 1 - i, K as u64 - 1 - i)
        }

        indexer
    }

    /**
     * Calculate the lexicographic rank (the index) of a permutation in O(n)
     * complexity.
     */
    fn rank(&self, perm: &[u64; K]) -> usize {
        // This will hold the Lehmer code (in a factorial number system).
        let mut lehmer: [u64; K] = [0; K];

        // Set of "seen" digits in the permutation.
        let mut seen = BasicBitSet::<N>::new();

        // The first digit of the Lehmer code is always the first digit of
        // the permutation.
        lehmer[0] = perm[0];

        // Mark the digit as seen (bitset uses right-to-left indexing).
        seen.set_bit(N - 1 - perm[0] as usize);

        for i in 1..K {
            seen.set_bit(N - 1 - perm[i] as usize);

            // The number of "seen" digits to the left of this digit is the
            // count of ones left of this digit.
            let num_ones =
                self.ones_count_lookup[seen.get_val() as usize >> (N - perm[i] as usize)] as u64;

            lehmer[i] = perm[i] - num_ones;
        }

        // Convert the Lehmer code to base-10.
        let mut index = 0;

        for i in 0..K {
            index += lehmer[i] * self.factorials[i];
        }
        return index as usize;
    }
}

#[cfg(test)]
mod tests {
    use crate::db::lehmer_indexer::PermutationIndexer;
    use rstest::rstest;

    #[test]
    fn it_works() {
        let indexer = PermutationIndexer::<4, 4>::new();

        assert_eq!(indexer.rank(&[2, 3, 0, 1]), 16)
    }

    #[rstest]
    #[case(0, [0, 1, 2])]
    #[case(1, [0, 2, 1])]
    #[case(2, [1, 0, 2])]
    #[case(3, [1, 2, 0])]
    #[case(4, [2, 0, 1])]
    #[case(5, [2, 1, 0])]
    fn three_bit_indexer_works(#[case] index: usize, #[case] permutation: [u64; 3]) {
        let indexer = PermutationIndexer::<3, 3>::new();

        assert_eq!(indexer.rank(&permutation), index)
    }

    #[rstest]
    #[case(0, [0, 1])]
    #[case(1, [0, 2])]
    #[case(2, [0, 3])]
    #[case(3, [1, 0])]
    #[case(4, [1, 2])]
    #[case(5, [1, 3])]
    #[case(6, [2, 0])]
    #[case(7, [2, 1])]
    #[case(8, [2, 3])]
    #[case(9, [3, 0])]
    #[case(10, [3, 1])]
    #[case(11, [3, 2])]
    fn four_choose_two_works(#[case] index: usize, #[case] permutation: [u64; 2]) {
        let indexer = PermutationIndexer::<4, 2>::new();

        assert_eq!(indexer.rank(&permutation), index)
    }
}
