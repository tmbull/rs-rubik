mod array_cube;
mod cube;

use crate::array_cube::ArrayCube;
use crate::cube::Cube;

fn main() {
    let mut cube = ArrayCube::default();
    cube.randomize(5, 5);

    let soln = cube.solve();
    println!("{:?}", soln)
}
