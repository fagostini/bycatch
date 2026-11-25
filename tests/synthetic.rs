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

/// Writes `name` genomes as FASTA in `dir`; returns the list lines `name<TAB>path[<TAB>group]`.
fn write_genomes(dir: &std::path::Path, genomes: &[(&str, &[u8], Option<&str>)]) -> String {
    let mut list = String::new();
    for (name, g, group) in genomes {
        let path = dir.join(name).to_str().unwrap().to_string();
        std::fs::write(&path, format!(">{name}\n{}\n", String::from_utf8_lossy(g))).unwrap();
        writeln!(list, "{name}\t{path}{}", group.map_or(String::new(), |gr| format!("\t{gr}"))).unwrap();
    }
    list
}

fn reads_fq(reads: &[Vec<u8>]) -> String {
    let mut fq = String::new();
    for (i, r) in reads.iter().enumerate() {
        writeln!(fq, "@r{i}\n{}\n+\n{}", String::from_utf8_lossy(r), "I".repeat(r.len())).unwrap();
    }
    fq
}

/// Two organisms of one group share a block: reads from it go to the group, not to
/// `shared`, and, with one of them as host, the group is not contamination. `add`
/// gives the same index as a single build.
#[test]
fn groups() {
    let dir = std::env::temp_dir().join(format!("bycatch_groups_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = |name: &str| dir.join(name).to_str().unwrap().to_string();
    let mut rng = Rng(0x9E3779B97F4A7C15);
    let block = rng.seq(3000);
    let (mut a, mut b, c) = (rng.seq(100_000), rng.seq(100_000), rng.seq(100_000));
    a[1000..4000].copy_from_slice(&block);
    b[5000..8000].copy_from_slice(&block);
    let mut reads = Vec::new();
    for (g, n) in [(&a, 800), (&c, 100)] {
        for _ in 0..n {
            let s = rng.next() as usize % (g.len() - 100);
            reads.push(g[s..s + 100].to_vec());
        }
    }
    for _ in 0..100 {
        let s = rng.next() as usize % (block.len() - 100);
        reads.push(block[s..s + 100].to_vec());
    }
    std::fs::write(p("reads.fq"), reads_fq(&reads)).unwrap();

    let scan = |idx: &str| run(&["scan", "-i", &p(idx), "--host", "a", "--every", "1", &p("reads.fq")]);
    std::fs::write(p("plain.tsv"), write_genomes(&dir, &[("a", &a, None), ("b", &b, None), ("c", &c, None)])).unwrap();
    run(&["build", &p("plain.tsv"), "--out", &p("plain")]);
    let plain = parse(&scan("plain"));
    assert!(plain["shared"] >= 95.0, "block reads are shared without groups: {plain:?}");

    let grouped = write_genomes(&dir, &[("a", &a, Some("ab")), ("b", &b, Some("ab")), ("c", &c, None)]);
    std::fs::write(p("grouped.tsv"), &grouped).unwrap();
    run(&["build", &p("grouped.tsv"), "--out", &p("grouped")]);
    let report = scan("grouped");
    let r = parse(&report);
    assert!(r["ab"] >= 95.0 && r["shared"] <= 5.0, "block reads go to the group: {r:?}");
    // a's reads that fall wholly inside the block (~3 % of a) also go to the group
    assert!(r["a"] + r["ab"] >= 895.0 && r["a"] >= 760.0 && r["c"] >= 99.0, "{r:?}");
    assert!(report.lines().any(|l| l.starts_with("ab\t") && l.contains("group,host_group")), "{report}");
    assert!(report.contains("# absent_below_percent\t0.38"), "Wilson upper end for 0 of 1000: {report}");
    let info = run(&["info", "-i", &p("grouped")]);
    assert!(info.contains("# organisms\t3\n# groups\t1\n") && info.contains("\nab\t") && info.contains("\ta,b\n"), "{info}");
    // contamination: c only, not the host's group
    let contam: f64 = report.lines().find(|l| l.starts_with("# contamination_percent")).unwrap().split('\t').nth(1).unwrap().parse().unwrap();
    assert!((contam - 100.0 * r["c"] / 1000.0).abs() < 1e-3, "{contam} {r:?}");

    // add: a first, then b (same group) and c
    let mut lines = grouped.lines();
    std::fs::write(p("first.tsv"), format!("{}\n", lines.next().unwrap())).unwrap();
    std::fs::write(p("rest.tsv"), lines.map(|l| format!("{l}\n")).collect::<String>()).unwrap();
    run(&["build", &p("first.tsv"), "--out", &p("first")]);
    run(&["add", "-i", &p("first"), &p("rest.tsv"), "--out", &p("added")]);
    assert_eq!(scan("added"), report, "add gives the same index as one build");
    std::fs::remove_dir_all(&dir).unwrap();
}

/// Over 255 organisms the index switches to wider ids; `add` across the limit gives the
/// same result as one build.
#[test]
fn more_than_255_organisms() {
    let dir = std::env::temp_dir().join(format!("bycatch_wide_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = |name: &str| dir.join(name).to_str().unwrap().to_string();
    let mut rng = Rng(0xD1B54A32D192ED03);
    let genomes: Vec<(String, Vec<u8>)> = (0..300).map(|i| (format!("g{i}"), rng.seq(3000))).collect();
    let entries: Vec<(&str, &[u8], Option<&str>)> = genomes.iter().map(|(n, g)| (n.as_str(), g.as_slice(), None)).collect();
    let list = write_genomes(&dir, &entries);
    std::fs::write(p("all.tsv"), &list).unwrap();
    let mut reads = Vec::new();
    for g in [0, 150, 299] {
        for k in 0..20 {
            reads.push(genomes[g].1[k * 100..k * 100 + 100].to_vec());
        }
    }
    std::fs::write(p("reads.fq"), reads_fq(&reads)).unwrap();
    run(&["build", &p("all.tsv"), "--out", &p("all")]);
    assert!(std::fs::read_to_string(p("all.meta")).unwrap().starts_with("contam3\t"));
    let report = run(&["scan", "-i", &p("all"), "--every", "1", &p("reads.fq")]);
    let r = parse(&report);
    assert_eq!((r["g0"], r["g150"], r["g299"]), (20.0, 20.0, 20.0), "{report}");

    let lines: Vec<&str> = list.lines().collect();
    std::fs::write(p("first.tsv"), lines[..250].join("\n") + "\n").unwrap();
    std::fs::write(p("rest.tsv"), lines[250..].join("\n") + "\n").unwrap();
    run(&["build", &p("first.tsv"), "--out", &p("first")]);
    assert!(std::fs::read_to_string(p("first.meta")).unwrap().starts_with("contam2\t"));
    run(&["add", "-i", &p("first"), &p("rest.tsv"), "--out", &p("added")]);
    assert_eq!(run(&["scan", "-i", &p("added"), "--every", "1", &p("reads.fq")]), report);
    std::fs::remove_dir_all(&dir).unwrap();
}

/// `--depleted-out` writes exactly the pairs not assigned to the host, in input order,
/// with ids, sequences and qualities intact.
#[test]
fn depleted_out() {
    use std::io::Read;
    let dir = std::env::temp_dir().join(format!("bycatch_depleted_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = |name: &str| dir.join(name).to_str().unwrap().to_string();
    let mut rng = Rng(0xA0761D6478BD642F);
    let (host, other) = (rng.seq(200_000), rng.seq(50_000));
    std::fs::write(p("list.tsv"), write_genomes(&dir, &[("host", &host, None), ("other", &other, None)])).unwrap();
    run(&["build", &p("list.tsv"), "--out", &p("idx")]);
    let (mut fq1, mut fq2, mut expect1, mut expect2) = (String::new(), String::new(), String::new(), String::new());
    for i in 0..600 {
        // every 3rd pair from `other`, every 7th random; the rest host
        let (r1, r2) = match (i % 3, i % 7) {
            (_, 0) => (rng.seq(100), rng.seq(100)),
            (0, _) => {
                let s = rng.next() as usize % (other.len() - 300);
                (other[s..s + 100].to_vec(), revcomp(&other[s + 200..s + 300]))
            }
            _ => {
                let s = rng.next() as usize % (host.len() - 300);
                (host[s..s + 100].to_vec(), revcomp(&host[s + 200..s + 300]))
            }
        };
        let q: String = (0..100).map(|k| (b'!' + ((i + k) % 40) as u8) as char).collect();
        let (a, b) = (format!("@p{i} 1:N\n{}\n+\n{q}\n", String::from_utf8_lossy(&r1)), format!("@p{i} 2:N\n{}\n+\n{q}\n", String::from_utf8_lossy(&r2)));
        fq1 += &a;
        fq2 += &b;
        if i % 3 == 0 || i % 7 == 0 {
            expect1 += &a;
            expect2 += &b;
        }
    }
    std::fs::write(p("r1.fq"), fq1).unwrap();
    std::fs::write(p("r2.fq"), fq2).unwrap();
    let report = run(&["scan", "-i", &p("idx"), "--host", "host", "--every", "1", "--depleted-out", &p("dep"), &p("r1.fq"), &p("r2.fq")]);
    let gunzip = |f: &str| {
        let mut s = String::new();
        flate2::read::MultiGzDecoder::new(std::fs::File::open(p(f)).unwrap()).read_to_string(&mut s).unwrap();
        s
    };
    assert_eq!(gunzip("dep_R1.fastq.gz"), expect1, "{report}");
    assert_eq!(gunzip("dep_R2.fastq.gz"), expect2);
    // refused without --every 1
    let out = Command::new(env!("CARGO_BIN_EXE_bycatch"))
        .args(["scan", "-i", &p("idx"), "--host", "host", "--depleted-out", &p("x"), &p("r1.fq"), &p("r2.fq")])
        .output()
        .unwrap();
    assert!(!out.status.success());
    std::fs::remove_dir_all(&dir).unwrap();
}

/// `build --window 7` keeps about twice the minimizers: 50 bp reads get more of them,
/// fewer stay unassigned, and the index records its window.
#[test]
fn smaller_window() {
    let dir = std::env::temp_dir().join(format!("bycatch_window_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = |name: &str| dir.join(name).to_str().unwrap().to_string();
    let mut rng = Rng(0x94D049BB133111EB);
    let (a, b) = (rng.seq(200_000), rng.seq(200_000));
    std::fs::write(p("list.tsv"), write_genomes(&dir, &[("a", &a, None), ("b", &b, None)])).unwrap();
    let reads: Vec<Vec<u8>> = (0..2000)
        .map(|_| {
            let s = rng.next() as usize % (a.len() - 50);
            a[s..s + 50].to_vec()
        })
        .collect();
    std::fs::write(p("reads.fq"), reads_fq(&reads)).unwrap();
    run(&["build", &p("list.tsv"), "--out", &p("w15")]);
    run(&["build", &p("list.tsv"), "--out", &p("w7"), "--window", "7"]);
    assert!(run(&["info", "-i", &p("w7")]).contains("# w\t7\n"));
    let scan = |idx: &str| parse(&run(&["scan", "-i", &p(idx), "--every", "1", "--min-frac", "0", &p("reads.fq")]));
    let (r15, r7) = (scan("w15"), scan("w7"));
    assert!(r7["a"] > r15["a"] + 100.0, "w=7 assigns more 50 bp reads: {r15:?} {r7:?}");
    assert_eq!(r7["b"], 0.0);
    // out of range
    let out = Command::new(env!("CARGO_BIN_EXE_bycatch")).args(["build", &p("list.tsv"), "--out", &p("x"), "--window", "3"]).output().unwrap();
    assert!(!out.status.success());
    std::fs::remove_dir_all(&dir).unwrap();
}

/// BGZF (as bcl-convert writes), plain gzip and plain FASTQ give the same report.
#[test]
fn bgzf_input() {
    use std::io::Write;
    let dir = std::env::temp_dir().join(format!("bycatch_bgzfin_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = |name: &str| dir.join(name).to_str().unwrap().to_string();
    let mut rng = Rng(0xBF58476D1CE4E5B9);
    let (a, b) = (rng.seq(100_000), rng.seq(100_000));
    std::fs::write(p("list.tsv"), write_genomes(&dir, &[("a", &a, None), ("b", &b, None)])).unwrap();
    run(&["build", &p("list.tsv"), "--out", &p("idx")]);
    let reads: Vec<Vec<u8>> = (0..20_000)
        .map(|i| {
            let g = if i % 5 == 0 { &b } else { &a };
            let s = rng.next() as usize % (g.len() - 150);
            g[s..s + 150].to_vec()
        })
        .collect();
    let fq = reads_fq(&reads);
    std::fs::write(p("r.fq"), &fq).unwrap();
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(fq.as_bytes()).unwrap();
    std::fs::write(p("r.fq.gz"), gz.finish().unwrap()).unwrap();
    // BGZF: gzip members of at most 64 KB, each with a BC subfield giving its size
    let mut bgzf = Vec::new();
    for part in fq.as_bytes().chunks(50_000).chain(std::iter::once(&[][..])) {
        let mut e = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::default());
        e.write_all(part).unwrap();
        let cdata = e.finish().unwrap();
        bgzf.extend_from_slice(&[0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0]);
        bgzf.extend_from_slice(&((18 + cdata.len() + 8 - 1) as u16).to_le_bytes());
        bgzf.extend_from_slice(&cdata);
        let mut c = flate2::Crc::new();
        c.update(part);
        bgzf.extend_from_slice(&c.sum().to_le_bytes());
        bgzf.extend_from_slice(&(part.len() as u32).to_le_bytes());
    }
    std::fs::write(p("r.bgzf.gz"), bgzf).unwrap();
    let scan = |f: &str| run(&["scan", "-i", &p("idx"), "--host", "a", "--every", "1", "-t", "8", &p(f)]);
    let plain = scan("r.fq");
    assert_eq!(scan("r.fq.gz"), plain);
    assert_eq!(scan("r.bgzf.gz"), plain);
    assert!(parse(&plain)["b"] >= 3900.0, "{plain}");
    std::fs::remove_dir_all(&dir).unwrap();
}
