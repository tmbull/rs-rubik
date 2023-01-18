use criterion::{black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use rs_rubik::cube3d;
use rs_rubik::cube3d::Cube;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("naive_solve_1-4");
    for num_moves in 1usize..=4 {
        group.bench_with_input(
            BenchmarkId::from_parameter(num_moves),
            &num_moves,
            |b, num_moves| {
                b.iter_batched(
                    || {
                        let mut cube = Cube::default();
                        cube.randomize(*num_moves, *num_moves);
                        cube
                    },
                    |cube| cube.solve(),
                    BatchSize::SmallInput,
                )
            },
        );
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
