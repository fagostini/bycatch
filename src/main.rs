mod bgzf;
mod index;
mod kmer;
mod run;
mod scan;

use clap::{Args, Parser};

type Res<T> = Result<T, Box<dyn std::error::Error + Send + Sync>>;

/// Quantify contamination in FASTQ files with k-mers unique to each reference genome.
#[derive(Parser)]
enum Cmd {
    /// Build the reference index from a TSV of `name<TAB>fasta_path`
    Build {
        genomes: String,
        /// Output prefix: writes <OUT>.kmers, <OUT>.table and <OUT>.meta
        #[arg(short, long)]
        out: String,
        /// Genomes read at once. Each needs about 1-2 GB of RAM per Gb of genome
        #[arg(short, long, default_value_t = 4)]
        jobs: usize,
        /// Minimizer window in k-mers (5-31): smaller keeps more k-mers, for short reads,
        /// at the cost of a larger index (density about 2 / (w + 1))
        #[arg(long, default_value_t = kmer::W)]
        window: usize,
    },
    /// Add genomes (TSV of `name<TAB>fasta_path`) to an index, without reading its genomes again
    Add {
        /// Existing index prefix
        #[arg(short, long)]
        index: String,
        genomes: String,
        /// Output prefix of the new index (must differ from --index)
        #[arg(short, long)]
        out: String,
        /// Genomes read at once
        #[arg(short, long, default_value_t = 4)]
        jobs: usize,
    },
    /// Print an index's layout, organisms (with unique markers) and groups
    Info {
        /// Index prefix
        #[arg(short, long)]
        index: String,
    },
    /// Rewrite an index in the current, smaller layout; results stay the same
    Repack {
        /// Existing index prefix
        #[arg(short, long)]
        index: String,
        /// Output prefix of the new index (must differ from --index)
        #[arg(short, long)]
        out: String,
    },
    /// Scan a FASTQ file (plain or gzipped, `-` for stdin), or a read pair
    Scan {
        #[command(flatten)]
        common: Common,
        /// Also write one line per classified read (id, category, hits, sequence) to this TSV
        #[arg(long)]
        reads: Option<String>,
        /// Write a random sample of unmatched fragments to this FASTA (for BLAST, sylph...)
        #[arg(long)]
        unmatched_out: Option<String>,
        /// Write the reads not assigned to the host (or its group) to PREFIX_R1/_R2.fastq.gz
        /// (PREFIX.fastq.gz for single-end); needs --host and --every 1
        #[arg(long, value_name = "PREFIX")]
        depleted_out: Option<String>,
        /// Fragments in --unmatched-out
        #[arg(long, default_value_t = 1000)]
        unmatched_n: usize,
        /// Write MultiQC custom-content files <MQC>_bycatch_mqc.tsv and <MQC>_bycatch_stats_mqc.tsv
        #[arg(long)]
        mqc: Option<String>,
        /// Sample name in the MultiQC files (default: from the R1 file name)
        #[arg(long)]
        name: Option<String>,
        r1: String,
        r2: Option<String>,
    },
    /// Scan all samples of a run from a TSV of `sample<TAB>lane<TAB>R1[<TAB>R2]`, and
    /// compare them: matrix, background / elevated / index-hopping flags, MultiQC files
    Run {
        #[command(flatten)]
        common: Common,
        samples: String,
        /// Output prefix: <OUT>.<sample>.<lane>.tsv, <OUT>.matrix.tsv, <OUT>.flags.tsv, MultiQC files
        #[arg(short, long)]
        out: String,
        /// Samples scanned at once, each with --threads classifier threads
        #[arg(short, long, default_value_t = 4)]
        jobs: usize,
    },
}

#[derive(Args)]
struct Common {
    /// Index prefix given to `build --out`
    #[arg(short, long)]
    index: String,
    /// Organism the sample is expected to be; all others count as contamination
    #[arg(long)]
    host: Option<String>,
    /// Classify 1 read in N
    #[arg(long, default_value_t = 50, value_parser = clap::value_parser!(u64).range(1..))]
    every: u64,
    /// Stop after classifying this many reads (or pairs), even if the file goes on
    #[arg(long, value_parser = clap::value_parser!(u64).range(1..))]
    max_reads: Option<u64>,
    /// Stop once the contamination 95 % CI is within ± this many percentage points (needs --host)
    #[arg(long)]
    precision: Option<f64>,
    /// Share of a fragment's minimizers the top organism must hit to be assigned (0 = off)
    #[arg(long, default_value_t = 0.15, value_parser = unit_fraction)]
    min_frac: f64,
    /// Threads classifying reads (plus one reading R1 and one reading R2)
    #[arg(short = 't', long, default_value_t = 4, value_parser = clap::value_parser!(u64).range(1..))]
    threads: u64,
}

impl Common {
    fn opts(&self) -> scan::Opts {
        scan::Opts {
            every: self.every,
            max_reads: self.max_reads,
            threads: self.threads as usize,
            precision: self.precision,
            reads_out: None,
            unmatched_out: None,
            unmatched_n: 0,
            min_frac: self.min_frac,
            depleted_out: None,
        }
    }
}

fn unit_fraction(s: &str) -> Result<f64, String> {
    match s.parse::<f64>() {
        Ok(f) if (0.0..=1.0).contains(&f) => Ok(f),
        _ => Err(format!("{s} is not a number from 0 to 1")),
    }
}

/// Sample name from a FASTQ path: `x/S1_L001_R1_001.fastq.gz` -> `S1_L001`.
fn sample_name(path: &str) -> String {
    let mut n = path.rsplit('/').next().unwrap_or(path).to_string();
    for ext in [".gz", ".fastq", ".fq", ".fasta", ".fa", "_001", "_R1", ".R1", "_1"] {
        n = n.strip_suffix(ext).map(str::to_string).unwrap_or(n);
    }
    n
}

fn main() {
    let result = match Cmd::parse() {
        Cmd::Build { genomes, out, jobs, window } => index::build(&genomes, &out, jobs, window),
        Cmd::Add { index, genomes, out, jobs } => index::add(&index, &genomes, &out, jobs),
        Cmd::Repack { index, out } => index::repack(&index, &out),
        Cmd::Info { index } => index::info(&index),
        Cmd::Scan { common, reads, unmatched_out, unmatched_n, depleted_out, mqc, name, r1, r2 } => {
            let opts = scan::Opts { reads_out: reads, unmatched_out, unmatched_n, depleted_out, ..common.opts() };
            let name = name.unwrap_or_else(|| sample_name(&r1));
            index::Index::load(&common.index).and_then(|idx| {
                let mqc = mqc.as_deref().map(|m| (m, name.as_str()));
                scan::scan(&idx, common.host.as_deref(), &opts, &r1, r2.as_deref(), mqc)
            })
        }
        Cmd::Run { common, samples, out, jobs } => index::Index::load(&common.index)
            .and_then(|idx| run::run(&idx, &samples, &out, common.host.as_deref(), &common.opts(), jobs)),
    };
    if let Err(e) = result {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}
