use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use criterion::{Criterion, criterion_group, criterion_main};
use tempfile::TempDir;

const RECORDS: usize = 50_000;
const CHROMOSOMES: usize = 10;
const DENSE_COUNTS: [usize; 4] = [500, 1_000, 2_000, 4_000];
const DENSE_RELATION_RECORDS: usize = 2_000;

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
    cluster: PathBuf,
    dense_relation_a: PathBuf,
    dense_relation_b: PathBuf,
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
        let cluster = directory.path().join("cluster.bed");
        let dense_relation_a = directory.path().join("dense-relation-a.bed");
        let dense_relation_b = directory.path().join("dense-relation-b.bed");
        let genome = directory.path().join("genome.tsv");

        let mut unsorted_writer = writer(&unsorted);
        let mut sorted_writer = writer(&sorted);
        let mut a_writer = writer(&a);
        let mut b_writer = writer(&b);
        let mut cluster_writer = writer(&cluster);
        let mut genome_writer = writer(&genome);
        assert_eq!(RECORDS % CHROMOSOMES, 0);
        let per_chromosome = RECORDS / CHROMOSOMES;
        for chromosome_index in 1..=CHROMOSOMES {
            let chromosome = format!("chr{chromosome_index:02}");
            for index in 0..per_chromosome {
                let reverse = per_chromosome - index - 1;
                writeln!(
                    unsorted_writer,
                    "{chromosome}\t{}\t{}\tU{chromosome_index}-{reverse}",
                    reverse * 100 + 5,
                    reverse * 100 + 23
                )
                .unwrap();

                let group = index / 5;
                let member = index % 5;
                let merge_start = group * 100 + member * 10 + 1;
                writeln!(
                    sorted_writer,
                    "{chromosome}\t{merge_start}\t{}\tM{chromosome_index}-{index}",
                    merge_start + 15
                )
                .unwrap();
                let strand = if index % 2 == 0 { '+' } else { '-' };
                writeln!(
                    cluster_writer,
                    "{chromosome}\t{merge_start}\t{}\tC{chromosome_index}-{index}\t0\t{strand}",
                    merge_start + 15
                )
                .unwrap();

                let start = index * 100 + 5;
                writeln!(
                    a_writer,
                    "{chromosome}\t{start}\t{}\tA{chromosome_index}-{index}",
                    start + 18
                )
                .unwrap();
                let b_start = index * 100 + 10;
                writeln!(
                    b_writer,
                    "{chromosome}\t{b_start}\t{}\tB{chromosome_index}-{index}",
                    b_start + 6
                )
                .unwrap();
                if index % 50 == 0 {
                    writeln!(
                        b_writer,
                        "{chromosome}\t{b_start}\t{}\tB{chromosome_index}-{index}",
                        b_start + 6
                    )
                    .unwrap();
                }
            }
            writeln!(
                genome_writer,
                "{chromosome}\t{}",
                per_chromosome * 100 + 200
            )
            .unwrap();
        }

        let mut dense_relation_a_writer = writer(&dense_relation_a);
        let mut dense_relation_b_writer = writer(&dense_relation_b);
        for index in 0..DENSE_RELATION_RECORDS {
            writeln!(dense_relation_a_writer, "chr1\t10000\t10001\tA{index}").unwrap();
            writeln!(
                dense_relation_b_writer,
                "chr1\t{index}\t{}\tB{index}",
                20_000 + index
            )
            .unwrap();
        }

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
            cluster,
            dense_relation_a,
            dense_relation_b,
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

fn checked_output(command: &mut Command) -> Vec<u8> {
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

fn assert_same(label: &str, ours: &mut Command, upstream: &mut Command) {
    assert_eq!(
        checked_output(ours),
        checked_output(upstream),
        "{label} output differs from bedtools 2.31.1"
    );
}

fn verify_outputs(fixtures: &Fixtures, ours: &str) {
    assert_same(
        "sort",
        Command::new(ours).arg("sort").arg(&fixtures.unsorted),
        Command::new("bedtools")
            .args(["sort", "-i"])
            .arg(&fixtures.unsorted),
    );
    assert_same(
        "merge",
        Command::new(ours).arg("merge").arg(&fixtures.sorted),
        Command::new("bedtools")
            .args(["merge", "-i"])
            .arg(&fixtures.sorted),
    );
    assert_same(
        "intersect",
        Command::new(ours)
            .args(["intersect", "-a"])
            .arg(&fixtures.a)
            .arg("-b")
            .arg(&fixtures.b),
        Command::new("bedtools")
            .args(["intersect", "-a"])
            .arg(&fixtures.a)
            .arg("-b")
            .arg(&fixtures.b),
    );
    assert_same(
        "subtract",
        Command::new(ours)
            .args(["subtract", "-a"])
            .arg(&fixtures.a)
            .arg("-b")
            .arg(&fixtures.b),
        Command::new("bedtools")
            .args(["subtract", "-a"])
            .arg(&fixtures.a)
            .arg("-b")
            .arg(&fixtures.b),
    );
    assert_same(
        "complement",
        Command::new(ours)
            .arg("complement")
            .arg(&fixtures.sorted)
            .arg("-g")
            .arg(&fixtures.genome),
        Command::new("bedtools")
            .args(["complement", "-i"])
            .arg(&fixtures.sorted)
            .arg("-g")
            .arg(&fixtures.genome),
    );
    assert_same(
        "cluster",
        Command::new(ours).arg("cluster").arg(&fixtures.cluster),
        Command::new("bedtools")
            .args(["cluster", "-i"])
            .arg(&fixtures.cluster),
    );
    assert_same(
        "cluster same-strand",
        Command::new(ours)
            .args(["cluster", "--strand", "same"])
            .arg(&fixtures.cluster),
        Command::new("bedtools")
            .args(["cluster", "-s", "-i"])
            .arg(&fixtures.cluster),
    );
    assert_same(
        "window pairs",
        Command::new(ours)
            .args(["window", "--window", "25", "-a"])
            .arg(&fixtures.a)
            .arg("-b")
            .arg(&fixtures.b),
        Command::new("bedtools")
            .args(["window", "-w", "25", "-a"])
            .arg(&fixtures.a)
            .arg("-b")
            .arg(&fixtures.b),
    );
    assert_same(
        "window dense count",
        Command::new(ours)
            .args(["window", "--window", "0", "--report", "count", "-a"])
            .arg(&fixtures.dense_relation_a)
            .arg("-b")
            .arg(&fixtures.dense_relation_b),
        Command::new("bedtools")
            .args(["window", "-w", "0", "-c", "-a"])
            .arg(&fixtures.dense_relation_a)
            .arg("-b")
            .arg(&fixtures.dense_relation_b),
    );
    assert_same(
        "closest",
        Command::new(ours)
            .args(["closest", "-a"])
            .arg(&fixtures.a)
            .arg("-b")
            .arg(&fixtures.b),
        Command::new("bedtools")
            .args(["closest", "-a"])
            .arg(&fixtures.a)
            .arg("-b")
            .arg(&fixtures.b),
    );
    assert_same(
        "closest distance",
        Command::new(ours)
            .args(["closest", "--distance", "unsigned", "-a"])
            .arg(&fixtures.a)
            .arg("-b")
            .arg(&fixtures.b),
        Command::new("bedtools")
            .args(["closest", "-d", "-a"])
            .arg(&fixtures.a)
            .arg("-b")
            .arg(&fixtures.b),
    );
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
    verify_outputs(&fixtures, ours);
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
    group.bench_function("cluster/rsomics", |bench| {
        bench.iter(|| run(Command::new(ours).arg("cluster").arg(&fixtures.cluster)));
    });
    group.bench_function("cluster/bedtools", |bench| {
        bench.iter(|| {
            run(Command::new("bedtools")
                .args(["cluster", "-i"])
                .arg(&fixtures.cluster));
        });
    });
    group.bench_function("cluster_same_strand/rsomics", |bench| {
        bench.iter(|| {
            run(Command::new(ours)
                .args(["cluster", "--strand", "same"])
                .arg(&fixtures.cluster));
        });
    });
    group.bench_function("cluster_same_strand/bedtools", |bench| {
        bench.iter(|| {
            run(Command::new("bedtools")
                .args(["cluster", "-s", "-i"])
                .arg(&fixtures.cluster));
        });
    });
    group.bench_function("window_pairs/rsomics", |bench| {
        bench.iter(|| {
            run(Command::new(ours)
                .args(["window", "--window", "25", "-a"])
                .arg(&fixtures.a)
                .arg("-b")
                .arg(&fixtures.b));
        });
    });
    group.bench_function("window_pairs/bedtools", |bench| {
        bench.iter(|| {
            run(Command::new("bedtools")
                .args(["window", "-w", "25", "-a"])
                .arg(&fixtures.a)
                .arg("-b")
                .arg(&fixtures.b));
        });
    });
    group.bench_function("closest/rsomics", |bench| {
        bench.iter(|| {
            run(Command::new(ours)
                .args(["closest", "-a"])
                .arg(&fixtures.a)
                .arg("-b")
                .arg(&fixtures.b));
        });
    });
    group.bench_function("closest/bedtools", |bench| {
        bench.iter(|| {
            run(Command::new("bedtools")
                .args(["closest", "-a"])
                .arg(&fixtures.a)
                .arg("-b")
                .arg(&fixtures.b));
        });
    });
    group.bench_function("closest_distance/rsomics", |bench| {
        bench.iter(|| {
            run(Command::new(ours)
                .args(["closest", "--distance", "unsigned", "-a"])
                .arg(&fixtures.a)
                .arg("-b")
                .arg(&fixtures.b));
        });
    });
    group.bench_function("closest_distance/bedtools", |bench| {
        bench.iter(|| {
            run(Command::new("bedtools")
                .args(["closest", "-d", "-a"])
                .arg(&fixtures.a)
                .arg("-b")
                .arg(&fixtures.b));
        });
    });
    group.finish();

    let mut relation_group = c.benchmark_group("dense_relation");
    relation_group.sample_size(10);
    relation_group.bench_function("window_count/rsomics", |bench| {
        bench.iter(|| {
            run(Command::new(ours)
                .args(["window", "--window", "0", "--report", "count", "-a"])
                .arg(&fixtures.dense_relation_a)
                .arg("-b")
                .arg(&fixtures.dense_relation_b));
        });
    });
    relation_group.bench_function("window_count/bedtools", |bench| {
        bench.iter(|| {
            run(Command::new("bedtools")
                .args(["window", "-w", "0", "-c", "-a"])
                .arg(&fixtures.dense_relation_a)
                .arg("-b")
                .arg(&fixtures.dense_relation_b));
        });
    });
    relation_group.finish();

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
