// TODO: Phase 3/4 — BAM query benchmarks
use criterion::{Criterion, criterion_group, criterion_main};

fn bench_placeholder(_c: &mut Criterion) {}

criterion_group!(benches, bench_placeholder);
criterion_main!(benches);
