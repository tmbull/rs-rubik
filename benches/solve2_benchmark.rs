use criterion::{black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use kiss3d::resource::BufferType::Array;
use rs_rubik::array_cube;
use rs_rubik::array_cube::ArrayCube;
use rs_rubik::cube::Cube;
use rs_rubik::cubie_cube::CubieCube;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("solve2_impl_comparison_1-6");
    for num_moves in 1usize..=6 {
        group.bench_with_input(
            BenchmarkId::new("ArrayCube", num_moves),
            &num_moves,
            |b, num_moves| {
                b.iter_batched(
                    || {
                        let mut cube = ArrayCube::default();
                        cube.randomize(*num_moves, *num_moves);
                        cube
                    },
                    |cube| cube.solve_iddfs(),
                    BatchSize::SmallInput,
                )
            },
        );
        group.bench_with_input(
            BenchmarkId::new("CubieCube", num_moves),
            &num_moves,
            |b, num_moves| {
                b.iter_batched(
                    || {
                        let mut cube = CubieCube::default();
                        cube.randomize(*num_moves, *num_moves);
                        cube
                    },
                    |cube| cube.solve_iddfs(),
                    BatchSize::SmallInput,
                )
            },
        );
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
