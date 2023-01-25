use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use rs_rubik::array_cube::ArrayCube;
use rs_rubik::cube::Cube;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("randomize_comparison_10-50");
    for num_moves in (10usize..=50).step_by(10) {
        group.bench_with_input(
            BenchmarkId::new("Arrays", num_moves),
            &num_moves,
            |b, num_moves| {
                b.iter_batched(
                    || {
                        let cube = ArrayCube::default();
                        cube
                    },
                    |mut cube| cube.randomize(*num_moves, *num_moves),
                    BatchSize::SmallInput,
                )
            },
        );
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
