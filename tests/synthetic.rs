//! End to end: three genomes of very different size sharing one block,
//! reads mixed 90 / 9 / 1 %, plus reads from the shared block and random reads.
use std::collections::HashMap;
use std::fmt::Write as _;
use std::process::Command;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn seq(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[(self.next() >> 20 & 3) as usize]).collect()
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|b| b"TGCA"[b"ACGT".iter().position(|x| x == b).unwrap()]).collect()
}

fn run(args: &[&str]) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_bycatch")).args(args).output().unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

/// first column -> second column
fn parse(report: &str) -> HashMap<String, f64> {
    report
        .lines()
        .skip(1)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].to_string(), f[1].parse().unwrap())
        })
        .collect()
}

#[test]
fn mixed_sample() {
    let dir = std::env::temp_dir().join(format!("bycatch_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = |name: &str| dir.join(name).to_str().unwrap().to_string();
    let mut rng = Rng(0x2545F4914F6CDD1D);

    let block = rng.seq(500);
    let mut list = String::new();
    let mut genomes = Vec::new();
    for (name, size) in [("big", 1_000_000), ("mid", 100_000), ("small", 10_000)] {
        let mut g = rng.seq(size);
        g[size / 2..size / 2 + 500].copy_from_slice(&block);
        std::fs::write(p(name), format!(">{name}\n{}\n", String::from_utf8_lossy(&g))).unwrap();
        writeln!(list, "{name}\t{}", p(name)).unwrap();
        genomes.push(g);
    }
    std::fs::write(p("genomes.tsv"), list).unwrap();

    let mut fastq = String::new();
    let mut add = |read: Vec<u8>, flip: bool| {
        let read = if flip { revcomp(&read) } else { read };
        writeln!(fastq, "@r\n{}\n+\n{}", String::from_utf8_lossy(&read), "I".repeat(100)).unwrap();
    };
    for (g, n) in genomes.iter().zip([18_000, 1_800, 200]) {
        for i in 0..n {
            let start = rng.next() as usize % (g.len() - 100);
            add(g[start..start + 100].to_vec(), i % 2 == 0);
        }
    }
    for i in 0..500 {
        let start = rng.next() as usize % 400;
        add(block[start..start + 100].to_vec(), i % 2 == 0);
        add(rng.seq(100), false);
    }
    std::fs::write(p("reads.fq"), fastq).unwrap();

    run(&["build", &p("genomes.tsv"), "--out", &p("idx")]);
    let report = run(&["scan", "-i", &p("idx"), "--host", "big", "--every", "1", &p("reads.fq")]);
    let r = parse(&report);
    println!("{report}");

    // a few genome reads fall inside the shared block and are rightly not assigned
    assert!((17_950.0..=18_000.0).contains(&r["big"]));
    assert!((1_760.0..=1_800.0).contains(&r["mid"]));
    assert!((170.0..=200.0).contains(&r["small"]));
    assert!((500.0..=560.0).contains(&r["shared"]));
    assert_eq!(r["unmatched"], 500.0);
    assert!(!report.contains("low_diversity"), "random reads over whole genomes must not be flagged");
    assert!(r["conflict"] <= 5.0);
    assert_eq!(r["# sampled_fragments"], 21_000.0);
    let contam = 100.0 * (r["mid"] + r["small"]) / 21_000.0;
    assert!((r["# contamination_percent"] - contam).abs() < 1e-3);
    // the CI holds the true 2,000 of 21,000 fragments; EM moves at most the weak + conflict reads
    let ci: Vec<f64> = report.lines().find(|l| l.starts_with("# contamination_ci95")).unwrap()
        .split('\t').skip(1).map(|x| x.parse().unwrap()).collect();
    let truth = 100.0 * 2_000.0 / 21_000.0;
    assert!(ci[0] < truth && truth < ci[1] && ci[1] - ci[0] < 1.0, "{ci:?}");
    let em = r["# em_contamination_percent"];
    assert!(em >= r["# contamination_percent"] - 1e-9 && em - contam <= 100.0 * (r["weak"] + r["conflict"]) / 21_000.0 + 1e-3);

    // any number of threads gives the same report and the same per-read file, in input order
    let threaded = |t: &str, out: &str| {
        let rep = run(&["scan", "-i", &p("idx"), "--every", "3", "-t", t, "--reads", &p(out), &p("reads.fq"), &p("reads.fq")]);
        (rep, std::fs::read_to_string(p(out)).unwrap())
    };
    assert_eq!(threaded("1", "t1.tsv"), threaded("8", "t8.tsv"));

    let args = ["scan", "-i", &p("idx"), "--every", "10", "--reads", &p("reads.tsv"), &p("reads.fq")];
    let r = parse(&run(&args));
    assert_eq!(r["# sampled_fragments"], 2_100.0);
    assert_eq!(r["# total_fragments"], 21_000.0);
    // per-read file: one line per sampled read, categories add up to the report
    let dump = std::fs::read_to_string(p("reads.tsv")).unwrap();
    let mut per_cat: HashMap<String, f64> = HashMap::new();
    for line in dump.lines().skip(1) {
        *per_cat.entry(line.split('\t').nth(1).unwrap().to_string()).or_default() += 1.0;
    }
    assert_eq!(per_cat.values().sum::<f64>(), 2_100.0);
    for (cat, n) in &per_cat {
        assert_eq!(r[cat], *n, "{cat}");
    }

    // --max-reads: stops after the 100th sampled read, record 1 + 99 * 10
    let report = run(&["scan", "-i", &p("idx"), "--every", "10", "--max-reads", "100", &p("reads.fq")]);
    let r = parse(&report);
    assert_eq!(r["# sampled_fragments"], 100.0);
    assert_eq!(r["# total_fragments"], 991.0);
    assert_eq!(r["# stopped_at_max_reads"], 100.0);
    // a limit above the number of sampled reads changes nothing
    let full = run(&["scan", "-i", &p("idx"), "--every", "10", &p("reads.fq")]);
    let capped = run(&["scan", "-i", &p("idx"), "--every", "10", "--max-reads", "5000", &p("reads.fq")]);
    assert_eq!(full, capped);

    // paired: the same file as R1 and R2 gives the same minimizers, so the same report
    let single = run(&["scan", "-i", &p("idx"), "--every", "7", &p("reads.fq")]);
    let paired = run(&["scan", "-i", &p("idx"), "--every", "7", &p("reads.fq"), &p("reads.fq")]);
    assert_eq!(single, paired);
    // R2 shorter or longer than R1 is an error, not a silent truncation
    let all = std::fs::read_to_string(p("reads.fq")).unwrap();
    let short: String = all.lines().take(4 * 20_000).map(|l| format!("{l}\n")).collect();
    std::fs::write(p("short.fq"), short).unwrap();
    for (a, b) in [("reads.fq", "short.fq"), ("short.fq", "reads.fq")] {
        let out = Command::new(env!("CARGO_BIN_EXE_bycatch"))
            .args(["scan", "-i", &p("idx"), "--every", "7", &p(a), &p(b)])
            .output()
            .unwrap();
        assert!(!out.status.success(), "{a} + {b} must fail");
    }
    // stopping early on a pair: R2 is not read to the end, so no count check, no hang
    let r = parse(&run(&["scan", "-i", &p("idx"), "--every", "7", "--max-reads", "50", &p("reads.fq"), &p("short.fq")]));
    assert_eq!(r["# sampled_fragments"], 50.0);

    // adapter read-through, poly-G and random fragments each get their own category
    let mut fq = String::new();
    let mut put = |read: Vec<u8>| writeln!(fq, "@x{} extra\n{}\n+\n{}", fq.len(), String::from_utf8_lossy(&read), "I".repeat(read.len())).unwrap();
    for _ in 0..300 {
        let mut r = rng.seq(66);
        r.extend_from_slice(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC");
        put(r);
    }
    (0..200).for_each(|_| put(vec![b'G'; 100]));
    (0..100).for_each(|_| put(rng.seq(100)));
    std::fs::write(p("extra.fq"), fq).unwrap();
    let args = ["scan", "-i", &p("idx"), "--every", "1", "--unmatched-out", &p("un.fa"), "--unmatched-n", "50", &p("extra.fq")];
    let r = parse(&run(&args));
    assert_eq!((r["adapter"], r["low_complexity"], r["unmatched"]), (300.0, 200.0, 100.0));
    let fasta = std::fs::read_to_string(p("un.fa")).unwrap();
    assert_eq!(fasta.lines().filter(|l| l.starts_with('>')).count(), 50);
    assert!(fasta.lines().all(|l| !l.contains(' ')), "ids cut at the first space");

    // --precision: checked after each batch of 4,096, from 10,000 fragments on
    let r = parse(&run(&["scan", "-i", &p("idx"), "--host", "big", "--every", "1", "--precision", "2", &p("reads.fq")]));
    assert_eq!(r["# stopped_at_precision"], 12_288.0);
    assert_eq!(r["# total_fragments"], 12_288.0);
    let out = Command::new(env!("CARGO_BIN_EXE_bycatch"))
        .args(["scan", "-i", &p("idx"), "--precision", "2", &p("reads.fq")])
        .output()
        .unwrap();
    assert!(!out.status.success(), "--precision without --host must fail");

    // add: big + mid, then small, gives the same index as all three at once
    let all = std::fs::read_to_string(p("genomes.tsv")).unwrap();
    let lines: Vec<&str> = all.lines().collect();
    std::fs::write(p("first.tsv"), format!("{}\n{}\n", lines[0], lines[1])).unwrap();
    std::fs::write(p("last.tsv"), format!("{}\n", lines[2])).unwrap();
    run(&["build", &p("first.tsv"), "--out", &p("idx12")]);
    run(&["add", "-i", &p("idx12"), &p("last.tsv"), "--out", &p("idx123")]);
    for ext in ["kmers", "table", "meta"] {
        let (a, b) = (std::fs::read(p(&format!("idx.{ext}"))).unwrap(), std::fs::read(p(&format!("idx123.{ext}"))).unwrap());
        assert!(a == b, "add differs from build in .{ext}");
    }

    // run: per-sample reports are the scan reports; matrix, flags and MultiQC files
    let list = format!("s1\tL1\t{}\ns2\tL1\t{}\t{}\n# comment\ns3\tL2\t{}\n", p("reads.fq"), p("reads.fq"), p("reads.fq"), p("extra.fq"));
    std::fs::write(p("samples.tsv"), list).unwrap();
    run(&["run", "-i", &p("idx"), "--host", "big", "--every", "7", "-t", "2", "-j", "2", &p("samples.tsv"), "--out", &p("run")]);
    let single = run(&["scan", "-i", &p("idx"), "--host", "big", "--every", "7", "--mqc", &p("one"), "--name", "s1", &p("reads.fq")]);
    assert_eq!(std::fs::read_to_string(p("run.s1.L1.tsv")).unwrap(), single);
    let matrix = std::fs::read_to_string(p("run.matrix.tsv")).unwrap();
    assert_eq!(matrix.lines().count(), 4);
    assert!(matrix.lines().nth(3).unwrap().starts_with("s3\tL2\t"));
    assert!(std::fs::read_to_string(p("run.flags.tsv")).unwrap().starts_with("sample\tlane\torganism"));
    let mqc = std::fs::read_to_string(p("run_bycatch_mqc.tsv")).unwrap();
    assert!(mqc.contains("# plot_type: 'bargraph'") && mqc.lines().filter(|l| !l.starts_with('#')).count() == 4);
    let one = std::fs::read_to_string(p("one_bycatch_stats_mqc.tsv")).unwrap();
    assert!(one.lines().last().unwrap().starts_with("s1\t"));

    std::fs::remove_dir_all(&dir).unwrap();
}
