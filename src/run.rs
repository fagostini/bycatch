//! Many samples of one sequencing run: a report per sample, a samples × organisms
//! matrix, flags that compare each sample with the rest of the run, and MultiQC files.
use crate::index::Index;
use crate::scan::{host_id, scan_sample, Opts, Report, SPECIAL};
use crate::Res;
use std::fs;
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

pub struct Sample {
    pub name: String,
    pub lane: String,
    r1: String,
    r2: Option<String>,
}

/// `sample<TAB>lane<TAB>R1[<TAB>R2]`; blank lines and `#` lines are skipped.
fn read_list(path: &str) -> Res<Vec<Sample>> {
    let text = fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let mut samples: Vec<Sample> = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').map(str::trim).collect();
        if !(3..=4).contains(&f.len()) || f[..3].iter().any(|x| x.is_empty()) {
            return Err(format!("{path}: expected sample<TAB>lane<TAB>R1[<TAB>R2], got '{line}'").into());
        }
        if samples.iter().any(|s| s.name == f[0] && s.lane == f[1]) {
            return Err(format!("{path}: sample '{}' lane '{}' listed twice", f[0], f[1]).into());
        }
        let r2 = f.get(3).filter(|x| !x.is_empty()).map(|x| x.to_string());
        samples.push(Sample { name: f[0].into(), lane: f[1].into(), r1: f[2].into(), r2 });
    }
    if samples.is_empty() {
        return Err(format!("{path}: no samples").into());
    }
    Ok(samples)
}

/// Reads per sample below which an organism counts as not detected.
const MIN_READS: u64 = 2;
/// Background: detected in at least this share of the samples (and at least this many samples)...
const BACKGROUND_PREVALENCE: f64 = 0.8;
const BACKGROUND_MIN_SAMPLES: usize = 5;
/// ...and no more than this many times the run median in this sample.
const BACKGROUND_FOLD: f64 = 3.0;
/// Elevated: lower CI bound above this many times the run median, with at least MIN_ELEVATED reads.
const ELEVATED_FOLD: f64 = 10.0;
const MIN_ELEVATED: u64 = 5;
/// Index hopping moves up to about 2 % of a sample's reads to other samples on its lane
/// (patterned flow cells, combinatorial indexes); look for sources of at least 1 %.
const HOP_RATE: f64 = 0.02;
const HOP_SOURCE_PCT: f64 = 1.0;

pub struct Flag {
    pub sample: usize,
    pub organism: usize,
    pub median: f64,
    pub prevalence: f64,
    pub flags: Vec<String>,
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

/// One row per (sample, organism) detected, host excluded. Flags:
/// - `background`: in most samples of the run at a similar level. Reagent or kit
///   contamination, or sequence missing from the host assembly, not this sample's own.
/// - `elevated`: well above the run median. Most likely contamination of this sample.
/// - `index_hopping:<sample>`: at most ~2 % of a sample on the same lane that has a lot
///   of this organism. Possibly reads moved from there by index hopping.
pub fn flags(samples: &[Sample], reports: &[&Report], n_org: usize, host: Option<usize>) -> Vec<Flag> {
    let mut rows = Vec::new();
    let n = reports.len();
    for g in (0..n_org).filter(|&g| Some(g) != host) {
        let pct: Vec<f64> = reports.iter().map(|r| r.pct(r.reads[g] as f64)).collect();
        let med = median(pct.clone());
        let detected = reports.iter().filter(|r| r.reads[g] >= MIN_READS).count();
        let prevalence = detected as f64 / n as f64;
        for (s, r) in reports.iter().enumerate() {
            if r.reads[g] < MIN_READS {
                continue;
            }
            let mut f = Vec::new();
            if n >= BACKGROUND_MIN_SAMPLES && prevalence >= BACKGROUND_PREVALENCE && pct[s] <= BACKGROUND_FOLD * med {
                f.push("background".to_string());
            }
            if r.reads[g] >= MIN_ELEVATED && r.ci(r.reads[g]).0 > ELEVATED_FOLD * med {
                f.push("elevated".to_string());
            }
            let source = (0..n)
                .filter(|&o| o != s && samples[o].lane == samples[s].lane)
                .filter(|&o| pct[o] >= HOP_SOURCE_PCT && pct[s] <= HOP_RATE * pct[o])
                .max_by(|&a, &b| pct[a].total_cmp(&pct[b]));
            if let Some(o) = source {
                f.push(format!("index_hopping:{}", samples[o].name));
            }
            rows.push(Flag { sample: s, organism: g, median: med, prevalence, flags: f });
        }
    }
    rows.sort_by_key(|f| (f.sample, f.organism));
    rows
}

fn create(path: &str) -> Res<BufWriter<fs::File>> {
    Ok(BufWriter::new(fs::File::create(path).map_err(|e| format!("{path}: {e}"))?))
}

/// Distinct colours (Tableau 20 without its dark grey) for the most abundant columns of a
/// bar graph, in Tableau's order: each hue dark then light. The rest are pooled into one
/// "other" column in Tableau's dark grey, so no colour stands for two things (MultiQC's
/// default palette repeats after about 10).
const PALETTE: [&str; 19] = [
    "#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F", "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
    "#E15759", "#FF9D9A", "#BAB0AC", "#D37295", "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6",
];
const OTHER: (&str, &str) = ("other", "#79706E");

/// A bar-graph column: (key, display name, fragments per row).
type Column = (String, String, Vec<u64>);

/// Columns with at least one fragment, most abundant first (mean % of sampled fragments
/// over the rows, so deep samples do not outweigh the others; ties keep input order).
/// The top PALETTE.len() keep their own colour, the rest are summed as "<other> (N)".
fn top_columns(rows: &[(String, &Report)], cols: Vec<Column>, other: &str) -> Vec<(Column, &'static str)> {
    let mean = |c: &Column| c.2.iter().zip(rows).map(|(&x, (_, r))| r.pct(x as f64)).sum::<f64>();
    let mut cols: Vec<Column> = cols.into_iter().filter(|c| c.2.iter().any(|&x| x > 0)).collect();
    cols.sort_by(|a, b| mean(b).total_cmp(&mean(a)));
    let rest = cols.split_off(cols.len().min(PALETTE.len()));
    let mut out: Vec<(Column, &str)> = cols.into_iter().zip(PALETTE).collect();
    if !rest.is_empty() {
        let sum = (0..rows.len()).map(|i| rest.iter().map(|c| c.2[i]).sum()).collect();
        out.push(((OTHER.0.into(), format!("{other} ({})", rest.len()), sum), OTHER.1));
    }
    out
}

/// `# categories:` block. Every category gets a name: MultiQC 1.35 crashes on an empty one.
fn write_categories(w: &mut impl Write, cats: &[(Column, &str)]) -> Res<()> {
    writeln!(w, "# categories:")?;
    for ((key, name, _), color) in cats {
        writeln!(w, "#     {key}:")?;
        writeln!(w, "#         name: '{name}'")?;
        writeln!(w, "#         color: '{color}'")?;
    }
    Ok(())
}

/// Data rows of a bar graph, one per sample, values from `cell(report, fragments)`.
fn write_bars(w: &mut impl Write, rows: &[(String, &Report)], cats: &[(Column, &str)], cell: impl Fn(&Report, u64) -> String) -> Res<()> {
    let head: Vec<&str> = cats.iter().map(|c| c.0 .0.as_str()).collect();
    writeln!(w, "Sample\t{}", head.join("\t"))?;
    for (i, (name, r)) in rows.iter().enumerate() {
        let v: Vec<String> = cats.iter().map(|c| cell(r, c.0 .2[i])).collect();
        writeln!(w, "{name}\t{}", v.join("\t"))?;
    }
    Ok(())
}

/// MultiQC custom content, one file per report section:
/// - `<prefix>_bycatch_mqc.tsv`: all fragments, organisms (host included) and special
///   categories ranked together, shown as percentages by default;
/// - `<prefix>_bycatch_contaminants_mqc.tsv`: non-host organisms only, as % of sampled
///   fragments, so traces are visible next to a host at ~90 %;
/// - `<prefix>_bycatch_stats_mqc.tsv`: General Statistics columns.
///
/// In both bar graphs the columns are sorted by abundance and only the top 19 get their
/// own colour; the rest are summed as one "other" column.
pub fn write_mqc(prefix: &str, idx: &Index, rows: &[(String, &Report)]) -> Res<()> {
    let host = rows.first().and_then(|r| r.1.host);
    let organism = |g: usize, name: String| -> Column { (idx.names[g].clone(), name, rows.iter().map(|r| r.1.reads[g]).collect()) };

    let mut all: Vec<Column> = (0..idx.names.len())
        .map(|g| organism(g, if Some(g) == host { format!("{} (host)", idx.names[g]) } else { idx.names[g].clone() }))
        .collect();
    all.extend(SPECIAL.iter().enumerate().map(|(c, &n)| (n.to_string(), n.to_string(), rows.iter().map(|r| r.1.special[c]).collect())));
    let cats = top_columns(rows, all, "other");
    let mut w = create(&format!("{prefix}_bycatch_mqc.tsv"))?;
    writeln!(w, "# id: 'bycatch'")?;
    writeln!(w, "# section_name: 'bycatch: all fragments'")?;
    writeln!(
        w,
        "# description: 'Sampled fragments per organism and category, most abundant first; the least abundant pooled as other.'"
    )?;
    writeln!(w, "# plot_type: 'bargraph'")?;
    writeln!(w, "# pconfig:")?;
    writeln!(w, "#     id: 'bycatch_bargraph'")?;
    writeln!(w, "#     title: 'bycatch: sampled fragments'")?;
    writeln!(w, "#     ylab: 'Fragments'")?;
    writeln!(w, "#     cpswitch_c_active: False")?;
    write_categories(&mut w, &cats)?;
    write_bars(&mut w, rows, &cats, |_, x| x.to_string())?;
    w.flush()?;

    // MultiQC orders custom sections by id: bycatch, bycatch_contaminants, bycatch_flags
    let contaminants: Vec<Column> =
        (0..idx.names.len()).filter(|&g| Some(g) != host).map(|g| organism(g, idx.names[g].clone())).collect();
    let cats = top_columns(rows, contaminants, "other organisms");
    let mut w = create(&format!("{prefix}_bycatch_contaminants_mqc.tsv"))?;
    writeln!(w, "# id: 'bycatch_contaminants'")?;
    writeln!(w, "# section_name: 'bycatch: contaminants'")?;
    writeln!(
        w,
        "# description: 'Non-host organisms, in % of sampled fragments (strict rule), most abundant first; the least abundant pooled as other organisms.'"
    )?;
    writeln!(w, "# plot_type: 'bargraph'")?;
    writeln!(w, "# pconfig:")?;
    writeln!(w, "#     id: 'bycatch_organisms_bargraph'")?;
    writeln!(w, "#     title: 'bycatch: non-host organisms'")?;
    writeln!(w, "#     ylab: '% of sampled fragments'")?;
    writeln!(w, "#     cpswitch: False")?;
    writeln!(w, "#     tt_decimals: 4")?;
    write_categories(&mut w, &cats)?;
    write_bars(&mut w, rows, &cats, |r, x| format!("{:.4}", r.pct(x as f64)))?;
    w.flush()?;

    let mut w = create(&format!("{prefix}_bycatch_stats_mqc.tsv"))?;
    writeln!(w, "# id: 'bycatch_stats'")?;
    writeln!(w, "# plot_type: 'generalstats'")?;
    writeln!(w, "# headers:")?;
    for (col, title, desc, scale) in [
        ("contam_percent", "Contam", "Non-host fragments, % of sampled (bycatch)", "OrRd"),
        ("contam_ci95_high", "Contam CI high", "Upper end of the 95 % CI of Contam", "OrRd"),
        ("unmatched_percent", "Unmatched", "Fragments matching no reference, % of sampled", "Greys"),
        ("adapter_percent", "Adapter", "Adapter dimers and read-through, % of sampled", "Greys"),
    ] {
        writeln!(w, "#     {col}:")?;
        writeln!(w, "#         title: '{title}'")?;
        writeln!(w, "#         description: '{desc}'")?;
        writeln!(w, "#         suffix: '%'")?;
        writeln!(w, "#         format: '{{:,.3f}}'")?;
        writeln!(w, "#         min: 0")?;
        writeln!(w, "#         scale: '{scale}'")?;
    }
    writeln!(w, "Sample\tcontam_percent\tcontam_ci95_high\tunmatched_percent\tadapter_percent")?;
    for (name, r) in rows {
        let (c, hi) = match r.contamination() {
            Some(c) => (format!("{:.4}", r.pct(c as f64)), format!("{:.4}", r.ci(c).1)),
            None => (String::new(), String::new()),
        };
        let um = r.pct(r.special[5] as f64);
        writeln!(w, "{name}\t{c}\t{hi}\t{um:.4}\t{:.4}", r.pct(r.special[3] as f64))?;
    }
    w.flush()?;
    Ok(())
}

/// MultiQC table of a run's `elevated` and `index_hopping` rows: `<prefix>_bycatch_flags_mqc.tsv`.
/// `background` rows are many (one per sample) and expected, so the description only
/// names those organisms; MultiQC also turns tables of over 500 rows into violin plots.
fn write_mqc_flags(prefix: &str, idx: &Index, samples: &[Sample], reports: &[&Report], flags: &[Flag]) -> Res<()> {
    let mut background: Vec<(usize, usize)> = Vec::new(); // (organism, samples flagged)
    for f in flags.iter().filter(|f| f.flags.iter().any(|x| x == "background")) {
        match background.iter_mut().find(|b| b.0 == f.organism) {
            Some(b) => b.1 += 1,
            None => background.push((f.organism, 1)),
        }
    }
    let background: Vec<String> =
        background.iter().map(|(g, n)| format!("{} ({n} of {})", idx.names[*g], samples.len())).collect();
    let background = if background.is_empty() { "none".to_string() } else { background.join(", ") };
    let mut w = create(&format!("{prefix}_bycatch_flags_mqc.tsv"))?;
    writeln!(w, "# id: 'bycatch_flags'")?;
    writeln!(w, "# section_name: 'bycatch: flags'")?;
    writeln!(
        w,
        "# description: 'Organisms well above the run median (elevated), or possibly from a sample on the same lane (index_hopping). Background in this run, at a similar level in most samples: {background}. All rows are in the run''s flags.tsv.'"
    )?;
    writeln!(w, "# plot_type: 'table'")?;
    writeln!(w, "# pconfig:")?;
    writeln!(w, "#     id: 'bycatch_flags_table'")?;
    writeln!(w, "#     title: 'bycatch: flags'")?;
    writeln!(w, "#     only_defined_headers: False")?;
    writeln!(w, "# headers:")?;
    for (col, title, format, suffix) in [
        ("reads", "Reads", "{:,.0f}", ""),
        ("percent", "Percent", "{:,.4f}", "%"),
        ("run_median_percent", "Run median", "{:,.4f}", "%"),
        ("prevalence", "Prevalence", "{:,.2f}", ""),
    ] {
        writeln!(w, "#     {col}:")?;
        writeln!(w, "#         title: '{title}'")?;
        writeln!(w, "#         format: '{format}'")?;
        writeln!(w, "#         suffix: '{suffix}'")?;
    }
    writeln!(w, "Sample / organism\treads\tpercent\trun_median_percent\tprevalence\tflags")?;
    for f in flags.iter().filter(|f| f.flags.iter().any(|x| x != "background")) {
        let (s, r) = (&samples[f.sample], reports[f.sample]);
        let x = r.reads[f.organism];
        writeln!(
            w,
            "{}_{} {}\t{x}\t{:.4}\t{:.4}\t{:.3}\t{}",
            s.name,
            s.lane,
            idx.names[f.organism],
            r.pct(x as f64),
            f.median,
            f.prevalence,
            f.flags.join(",")
        )?;
    }
    w.flush()?;
    Ok(())
}

/// `bycatch run`: scans `jobs` samples at once, each with `opts.threads` classifier threads.
/// Writes `<out>.<sample>.<lane>.tsv`, `<out>.matrix.tsv`, `<out>.flags.tsv` and MultiQC files.
pub fn run(idx: &Index, list: &str, out: &str, host: Option<&str>, opts: &Opts, jobs: usize) -> Res<()> {
    let samples = read_list(list)?;
    let host = host_id(idx, host)?;
    let next = AtomicUsize::new(0);
    let results: Mutex<Vec<Option<Res<Report>>>> = Mutex::new((0..samples.len()).map(|_| None).collect());
    std::thread::scope(|scope| {
        for _ in 0..jobs.clamp(1, samples.len()) {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                let Some(s) = samples.get(i) else { return };
                let r = scan_sample(idx, opts, host, &s.r1, s.r2.as_deref()).and_then(|r| {
                    let path = format!("{out}.{}.{}.tsv", s.name, s.lane);
                    let mut w = create(&path)?;
                    r.write_tsv(idx, &mut w)?;
                    w.flush()?;
                    Ok(r)
                });
                match &r {
                    Ok(r) => eprintln!("{}\t{}\t{} fragments sampled", s.name, s.lane, r.sampled),
                    Err(e) => eprintln!("{}\t{}\terror: {e}", s.name, s.lane),
                }
                results.lock().unwrap()[i] = Some(r);
            });
        }
    });
    let results = results.into_inner().unwrap();
    let mut failed = Vec::new();
    let (mut ok_samples, mut reports) = (Vec::new(), Vec::new());
    for (s, r) in samples.into_iter().zip(results) {
        match r {
            Some(Ok(r)) => {
                ok_samples.push(s);
                reports.push(r);
            }
            _ => failed.push(format!("{} {}", s.name, s.lane)),
        }
    }
    let refs: Vec<&Report> = reports.iter().collect();

    let mut w = create(&format!("{out}.matrix.tsv"))?;
    let names: Vec<&str> = idx.names.iter().map(String::as_str).chain(SPECIAL).collect();
    writeln!(
        w,
        "sample\tlane\tsampled_fragments\tduplicate_design_effect\t{}\tcontamination_percent\tcontamination_ci95_low\tcontamination_ci95_high\tem_contamination_percent",
        names.join("\t")
    )?;
    for (s, r) in ok_samples.iter().zip(&refs) {
        let pct: Vec<String> = r.reads.iter().chain(&r.special).map(|&x| format!("{:.4}", r.pct(x as f64))).collect();
        let contam = match (r.contamination(), r.em_contamination()) {
            (Some(c), Some(e)) => {
                let (lo, hi) = r.ci(c);
                format!("{:.4}\t{lo:.4}\t{hi:.4}\t{:.4}", r.pct(c as f64), r.pct(e))
            }
            _ => "\t\t\t".to_string(),
        };
        writeln!(w, "{}\t{}\t{}\t{:.3}\t{}\t{contam}", s.name, s.lane, r.sampled, r.design, pct.join("\t"))?;
    }
    w.flush()?;

    let mut w = create(&format!("{out}.flags.tsv"))?;
    writeln!(w, "sample\tlane\torganism\treads\tpercent\tci95_low\tci95_high\trun_median_percent\tprevalence\tflags")?;
    let flag_rows = flags(&ok_samples, &refs, idx.names.len(), host);
    for f in &flag_rows {
        let (s, r) = (&ok_samples[f.sample], refs[f.sample]);
        let x = r.reads[f.organism];
        let (lo, hi) = r.ci(x);
        writeln!(
            w,
            "{}\t{}\t{}\t{x}\t{:.4}\t{lo:.4}\t{hi:.4}\t{:.4}\t{:.3}\t{}",
            s.name,
            s.lane,
            idx.names[f.organism],
            r.pct(x as f64),
            f.median,
            f.prevalence,
            f.flags.join(",")
        )?;
    }
    w.flush()?;

    let rows: Vec<(String, &Report)> =
        ok_samples.iter().zip(&refs).map(|(s, &r)| (format!("{}_{}", s.name, s.lane), r)).collect();
    write_mqc(out, idx, &rows)?;
    write_mqc_flags(out, idx, &ok_samples, &refs, &flag_rows)?;
    if failed.is_empty() {
        Ok(())
    } else {
        Err(format!("{} sample(s) failed (results of the others are written): {}", failed.len(), failed.join(", ")).into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn report(reads: Vec<u64>, sampled: u64) -> Report {
        let n = reads.len();
        Report {
            host: Some(0),
            est: reads.iter().map(|&x| x as f64).collect(),
            reads,
            hits: vec![0; n],
            distinct: vec![0; n],
            special: [0; 6],
            sampled,
            total: sampled,
            stop: None,
            design: 1.0,
            eps: 0.01,
        }
    }

    fn sample(name: &str, lane: &str) -> Sample {
        Sample { name: name.into(), lane: lane.into(), r1: String::new(), r2: None }
    }

    #[test]
    fn columns_ranked_by_mean_percent() {
        let col = |k: &str, v: Vec<u64>| -> Column { (k.into(), k.into(), v) };
        // b is most abundant on average (a's deep sample does not outweigh), z never seen
        let (r1, r2) = (report(vec![0], 1_000), report(vec![0], 10_000));
        let rows = vec![("1".to_string(), &r1), ("2".to_string(), &r2)];
        let cats = top_columns(&rows, vec![col("a", vec![5, 100]), col("b", vec![50, 200]), col("z", vec![0, 0])], "other");
        let keys: Vec<&str> = cats.iter().map(|c| c.0 .0.as_str()).collect();
        assert_eq!(keys, ["b", "a"]);
        // more columns than colours: the least abundant 3 are summed into "other (3)"
        let many: Vec<Column> = (0..22u64).map(|i| col(&format!("c{i}"), vec![100 - i, 0])).collect();
        let cats = top_columns(&rows, many, "other");
        assert_eq!(cats.len(), 20);
        assert_eq!(cats[0].0 .0, "c0");
        assert_eq!(cats[19].0 .1, "other (3)");
        assert_eq!(cats[19].0 .2, vec![81 + 80 + 79, 0]);
        let colors: std::collections::HashSet<&str> = cats.iter().map(|c| c.1).collect();
        assert_eq!(colors.len(), 20, "no colour repeats");
    }

    #[test]
    fn run_flags() {
        // organism 1 at ~0.05 % in all 6 samples (background), 2 % in sample 5;
        // organism 2 at 5 % in sample 0, and 0.03 % in sample 1 on the same lane (hopping)
        let n = 100_000;
        let mut reads = [[99_000, 50, 0], [99_000, 45, 30], [99_000, 55, 0], [99_000, 50, 0], [99_000, 48, 0], [97_000, 2_000, 0]];
        reads[0][2] = 5_000;
        let samples: Vec<Sample> =
            (0..6).map(|i| sample(&format!("s{i}"), if i < 3 { "L1" } else { "L2" })).collect();
        let reports: Vec<Report> = reads.iter().map(|r| report(r.to_vec(), n)).collect();
        let refs: Vec<&Report> = reports.iter().collect();
        let rows = flags(&samples, &refs, 3, Some(0));
        let get = |s: usize, g: usize| rows.iter().find(|f| f.sample == s && f.organism == g).map(|f| f.flags.join(","));
        assert_eq!(get(0, 1).unwrap(), "background");
        assert_eq!(get(5, 1).unwrap(), "elevated");
        assert_eq!(get(1, 2).unwrap(), "elevated,index_hopping:s0");
        assert!(get(0, 0).is_none(), "host has no rows");
        assert!(get(2, 2).is_none(), "not detected: no row");
    }
}
