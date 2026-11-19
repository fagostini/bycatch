//! 2-bit canonical k-mers and window minimizers. Build and scan must use the
//! same function, or the index and the reads would disagree on which k-mers exist.
use std::collections::VecDeque;

/// 27 bases = 54 bits, which leaves 8 bits of a u64 for the organism id.
pub const K: usize = 27;
/// Window of 15 k-mers: about 1 k-mer in 8 kept, at least one per 41 bp.
pub const W: usize = 15;
const MASK: u64 = (1 << (2 * K)) - 1;

fn code(b: u8) -> Option<u64> {
    match b | 0x20 {
        b'a' => Some(0),
        b'c' => Some(1),
        b'g' => Some(2),
        b't' => Some(3),
        _ => None,
    }
}

/// splitmix64 finalizer. Bijective, so equal hash means equal k-mer.
pub fn mix(mut x: u64) -> u64 {
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

/// Homopolymers and simple repeats: fewer than 8 distinct 3-mers out of 25.
pub fn low_complexity(kmer: u64) -> bool {
    let mut seen = 0u64;
    for i in 0..K - 2 {
        seen |= 1 << ((kmer >> (2 * i)) & 63);
    }
    seen.count_ones() < 8
}

/// Appends the canonical minimizer k-mers of `seq` to `out` (duplicates possible).
/// Any non-ACGT byte breaks the sequence.
pub fn minimizers(seq: &[u8], out: &mut Vec<u64>) {
    let mut emit = |kmer: u64| {
        if !low_complexity(kmer) {
            out.push(kmer)
        }
    };
    // (hash, k-mer, k-mer position in the current ACGT run), hashes increasing
    let mut q: VecDeque<(u64, u64, usize)> = VecDeque::new();
    let (mut fwd, mut rev, mut len) = (0u64, 0u64, 0usize);
    let mut last = usize::MAX;
    for c in seq.iter().map(|&b| code(b)).chain(std::iter::once(None)) {
        let Some(c) = c else {
            // run shorter than one window: keep its single minimum
            if len >= K && len - K + 1 < W {
                emit(q[0].1);
            }
            q.clear();
            len = 0;
            last = usize::MAX;
            continue;
        };
        fwd = ((fwd << 2) | c) & MASK;
        rev = (rev >> 2) | ((3 - c) << (2 * (K - 1)));
        len += 1;
        if len < K {
            continue;
        }
        let n = len - K;
        let kmer = fwd.min(rev);
        let h = mix(kmer);
        while q.back().is_some_and(|b| b.0 >= h) {
            q.pop_back();
        }
        q.push_back((h, kmer, n));
        if q[0].2 + W <= n {
            q.pop_front();
        }
        if n + 1 >= W && q[0].2 != last {
            last = q[0].2;
            emit(q[0].1);
        }
    }
}

/// Calls `f` with every canonical k-mer of `seq` (k <= 32), low-complexity ones included.
pub fn kmers(seq: &[u8], k: usize, mut f: impl FnMut(u64)) {
    kmers_at(seq, k, |_, kmer| f(kmer))
}

/// As `kmers`, also passing the k-mer's start position in `seq`.
pub fn kmers_at(seq: &[u8], k: usize, mut f: impl FnMut(usize, u64)) {
    let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let (mut fwd, mut rev, mut len) = (0u64, 0u64, 0usize);
    for (i, &b) in seq.iter().enumerate() {
        let Some(c) = code(b) else {
            len = 0;
            continue;
        };
        fwd = ((fwd << 2) | c) & mask;
        rev = (rev >> 2) | ((3 - c) << (2 * (k - 1)));
        len += 1;
        if len >= k {
            f(i + 1 - k, fwd.min(rev));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set(seq: &[u8]) -> Vec<u64> {
        let mut v = Vec::new();
        minimizers(seq, &mut v);
        v.sort_unstable();
        v.dedup();
        v
    }

    fn random_seq(n: usize) -> Vec<u8> {
        let mut x = 42u64;
        (0..n)
            .map(|_| {
                x = mix(x.wrapping_add(0x9e3779b97f4a7c15));
                b"ACGT"[(x & 3) as usize]
            })
            .collect()
    }

    #[test]
    fn strand_independent() {
        for n in [30, 500] {
            let seq = random_seq(n);
            let rc: Vec<u8> = seq
                .iter()
                .rev()
                .map(|b| match b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'A',
                })
                .collect();
            assert!(!set(&seq).is_empty());
            assert_eq!(set(&seq), set(&rc));
            assert_eq!(set(&seq), set(&seq.to_ascii_lowercase()));
        }
    }

    #[test]
    fn n_breaks_kmers() {
        let mut seq = random_seq(26);
        seq.push(b'N');
        seq.extend(random_seq(26));
        assert!(set(&seq).is_empty());
    }

    #[test]
    fn all_kmers_canonical() {
        let seq = random_seq(100);
        let (mut a, mut b) = (Vec::new(), Vec::new());
        kmers(&seq, K, |k| a.push(k));
        let rc: Vec<u8> = seq.iter().rev().map(|&c| b"TGCA"[b"ACGT".iter().position(|&x| x == c).unwrap()]).collect();
        kmers(&rc, K, |k| b.push(k));
        b.reverse();
        assert_eq!(a.len(), 100 - K + 1);
        assert_eq!(a, b);
        // every minimizer is one of the k-mers
        assert!(set(&seq).iter().all(|m| a.contains(m)));
    }

    #[test]
    fn low_complexity_dropped() {
        assert!(set(&[b'A'; 200]).is_empty());
        assert!(set(&b"AC".repeat(100)).is_empty());
    }

    #[test]
    fn every_window_covered() {
        // any 41 bp stretch must contain a minimizer of the full sequence
        let seq = random_seq(2000);
        let all = set(&seq);
        for start in 0..seq.len() - (K + W - 1) {
            let win = set(&seq[start..start + K + W - 1]);
            assert!(win.iter().all(|k| all.binary_search(k).is_ok()));
        }
    }
}
