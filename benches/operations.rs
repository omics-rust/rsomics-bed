use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use criterion::{Criterion, criterion_group, criterion_main};
use tempfile::TempDir;

const RECORDS: usize = 50_000;
const DENSE_COUNTS: [usize; 4] = [500, 1_000, 2_000, 4_000];

struct DenseFixtures {
    count: usize,
    a: PathBuf,
    b: PathBuf,
}

struct Fixtures {
    _directory: TempDir,
    unsorted: PathBuf,
    sorted: PathBuf,
    a: PathBuf,
    b: PathBuf,
    genome: PathBuf,
    dense: Vec<DenseFixtures>,
}

impl Fixtures {
    fn new() -> Self {
        let directory = tempfile::tempdir().expect("create benchmark directory");
        let unsorted = directory.path().join("unsorted.bed");
        let sorted = directory.path().join("sorted.bed");
        let a = directory.path().join("a.bed");
        let b = directory.path().join("b.bed");
        let genome = directory.path().join("genome.tsv");

        let mut unsorted_writer = writer(&unsorted);
        let mut sorted_writer = writer(&sorted);
        let mut a_writer = writer(&a);
        let mut b_writer = writer(&b);
        for index in 0..RECORDS {
            let reverse = RECORDS - index - 1;
            writeln!(
                unsorted_writer,
                "chr1\t{}\t{}\tU{reverse}",
                reverse * 4,
                reverse * 4 + 1
            )
            .unwrap();
            writeln!(
                sorted_writer,
                "chr1\t{}\t{}\tS{index}",
                index * 4,
                index * 4 + 2
            )
            .unwrap();
            writeln!(
                a_writer,
                "chr1\t{}\t{}\tA{index}",
                index * 4 + 2,
                index * 4 + 3
            )
            .unwrap();
            writeln!(b_writer, "chr1\t{}\t{}\tB{index}", index * 4, index * 4 + 1).unwrap();
        }
        std::fs::write(&genome, format!("chr1\t{}\n", RECORDS * 4 + 100)).unwrap();

        let dense = DENSE_COUNTS
            .into_iter()
            .map(|count| {
                let a = directory.path().join(format!("dense-a-{count}.bed"));
                let b = directory.path().join(format!("dense-b-{count}.bed"));
                let mut a_writer = writer(&a);
                let mut b_writer = writer(&b);
                for index in 0..count {
                    writeln!(a_writer, "chr1\t10000\t10001\tA{index}").unwrap();
                    writeln!(b_writer, "chr1\t{index}\t{}\tB{index}", 20_000 + index).unwrap();
                }
                DenseFixtures { count, a, b }
            })
            .collect();

        Self {
            _directory: directory,
            unsorted,
            sorted,
            a,
            b,
            genome,
            dense,
        }
    }
}

fn writer(path: &Path) -> BufWriter<File> {
    BufWriter::new(File::create(path).unwrap())
}

fn run(command: &mut Command) {
    let status = command
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(status.success());
}

fn benchmark(c: &mut Criterion) {
    let version = Command::new("bedtools")
        .arg("--version")
        .output()
        .expect("bedtools 2.31.1 is required for the performance oracle");
    assert!(version.status.success());
    assert!(
        String::from_utf8_lossy(&version.stdout).contains("v2.31.1"),
        "performance oracle must be bedtools 2.31.1"
    );

    let fixtures = Fixtures::new();
    let ours = env!("CARGO_BIN_EXE_rsomics-bed");
    let mut group = c.benchmark_group(format!("bed_operations/{RECORDS}"));
    group.sample_size(10);

    group.bench_function("sort/rsomics", |bench| {
        bench.iter(|| run(Command::new(ours).arg("sort").arg(&fixtures.unsorted)));
    });
    group.bench_function("sort/bedtools", |bench| {
        bench.iter(|| {
            run(Command::new("bedtools")
                .args(["sort", "-i"])
                .arg(&fixtures.unsorted));
        });
    });
    group.bench_function("merge/rsomics", |bench| {
        bench.iter(|| run(Command::new(ours).arg("merge").arg(&fixtures.sorted)));
    });
    group.bench_function("merge/bedtools", |bench| {
        bench.iter(|| {
            run(Command::new("bedtools")
                .args(["merge", "-i"])
                .arg(&fixtures.sorted));
        });
    });
    group.bench_function("intersect/rsomics", |bench| {
        bench.iter(|| {
            run(Command::new(ours)
                .args(["intersect", "-a"])
                .arg(&fixtures.a)
                .arg("-b")
                .arg(&fixtures.b));
        });
    });
    group.bench_function("intersect/bedtools", |bench| {
        bench.iter(|| {
            run(Command::new("bedtools")
                .args(["intersect", "-a"])
                .arg(&fixtures.a)
                .arg("-b")
                .arg(&fixtures.b));
        });
    });
    group.bench_function("subtract/rsomics", |bench| {
        bench.iter(|| {
            run(Command::new(ours)
                .args(["subtract", "-a"])
                .arg(&fixtures.a)
                .arg("-b")
                .arg(&fixtures.b));
        });
    });
    group.bench_function("subtract/bedtools", |bench| {
        bench.iter(|| {
            run(Command::new("bedtools")
                .args(["subtract", "-a"])
                .arg(&fixtures.a)
                .arg("-b")
                .arg(&fixtures.b));
        });
    });
    group.bench_function("complement/rsomics", |bench| {
        bench.iter(|| {
            run(Command::new(ours)
                .arg("complement")
                .arg(&fixtures.sorted)
                .arg("-g")
                .arg(&fixtures.genome));
        });
    });
    group.bench_function("complement/bedtools", |bench| {
        bench.iter(|| {
            run(Command::new("bedtools")
                .args(["complement", "-i"])
                .arg(&fixtures.sorted)
                .arg("-g")
                .arg(&fixtures.genome));
        });
    });
    group.finish();

    let mut dense_group = c.benchmark_group("dense_subtract");
    dense_group.sample_size(10);
    for fixture in &fixtures.dense {
        dense_group.bench_function(format!("{}/rsomics", fixture.count), |bench| {
            bench.iter(|| {
                run(Command::new(ours)
                    .args(["subtract", "-a"])
                    .arg(&fixture.a)
                    .arg("-b")
                    .arg(&fixture.b));
            });
        });
        dense_group.bench_function(format!("{}/bedtools", fixture.count), |bench| {
            bench.iter(|| {
                run(Command::new("bedtools")
                    .args(["subtract", "-a"])
                    .arg(&fixture.a)
                    .arg("-b")
                    .arg(&fixture.b));
            });
        });
    }
    dense_group.finish();
}

criterion_group!(benches, benchmark);
criterion_main!(benches);
