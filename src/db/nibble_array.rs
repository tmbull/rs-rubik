use std::ops::Index;

struct NibbleArray<const N: usize>
where
    [(); N / 2]: Sized,
{
    bytes: [u8; N / 2],
}

impl<const N: usize> NibbleArray<N>
where
    [(); N / 2]: Sized,
{
    fn new() -> Self {
        Self { bytes: [0; N / 2] }
    }

    fn get(&self, index: usize) -> u8 {
        let idx = index / 2;
        let val = self.bytes[idx];

        if index % 2 == 0 {
            &val & 0xF
        } else {
            &val >> 4
        }
    }

    // fn set(&self, index: usize, value: u8) {
    //     let idx = index / 2;
    //     if index % 2 == 0 {
    //         self.bytes[idx] =
    //     }
    // }
}
