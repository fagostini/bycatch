//! Parallel reader for BGZF, the blocked gzip that bcl-convert and samtools write: a
//! series of gzip members of at most 64 KB, each naming its own size. Blocks are
//! independent, so batches of them are inflated on several threads and handed over in
//! order. Plain gzip (one long member) cannot be split this way and keeps the usual reader.
use crate::Res;
use flate2::{Crc, Decompress, FlushDecompress, Status};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufReader, Read};
use std::sync::mpsc::{sync_channel, Receiver};
use std::thread::JoinHandle;

/// Blocks inflated per batch: about 16 MB of FASTQ.
const BATCH: usize = 256;

/// The size field of a BGZF header (`BC` subfield), if `head` starts one.
fn block_size(head: &[u8]) -> Option<usize> {
    if head.len() < 12 || head[..4] != [0x1f, 0x8b, 8, 4] {
        return None;
    }
    let xlen = u16::from_le_bytes([head[10], head[11]]) as usize;
    let extra = head.get(12..12 + xlen)?;
    let mut i = 0;
    while i + 4 <= extra.len() {
        let slen = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
        if extra[i] == b'B' && extra[i + 1] == b'C' && slen == 2 {
            let bsize = u16::from_le_bytes([*extra.get(i + 4)?, *extra.get(i + 5)?]) as usize;
            return Some(bsize + 1);
        }
        i += 4 + slen;
    }
    None
}

/// Whether the file starts with a BGZF block.
pub fn is_bgzf(path: &str) -> bool {
    let mut head = [0u8; 64];
    let n = File::open(path).and_then(|mut f| f.read(&mut head)).unwrap_or(0);
    block_size(&head[..n]).is_some()
}

/// One whole compressed block: header, deflate data, CRC32, input size.
fn inflate(block: &[u8]) -> Res<Vec<u8>> {
    let xlen = u16::from_le_bytes([block[10], block[11]]) as usize;
    let n = block.len();
    if n < 12 + xlen + 8 {
        return Err("BGZF block too short".into());
    }
    let crc = u32::from_le_bytes(block[n - 8..n - 4].try_into().unwrap());
    let isize = u32::from_le_bytes(block[n - 4..].try_into().unwrap()) as usize;
    let mut out = Vec::with_capacity(isize);
    let mut d = Decompress::new(false);
    let status = d.decompress_vec(&block[12 + xlen..n - 8], &mut out, FlushDecompress::Finish)?;
    let mut c = Crc::new();
    c.update(&out);
    if !matches!(status, Status::StreamEnd) || out.len() != isize || c.sum() != crc {
        return Err("corrupt BGZF block (size or CRC32 mismatch)".into());
    }
    Ok(out)
}

/// Reads raw blocks on a thread, inflates each batch on `threads` threads, and yields
/// the bytes in file order.
pub struct BgzfReader {
    rx: Receiver<Result<Vec<u8>, String>>,
    buf: Vec<u8>,
    pos: usize,
    thread: Option<JoinHandle<()>>,
}

impl BgzfReader {
    pub fn open(path: &str, threads: usize) -> Res<BgzfReader> {
        let mut r = BufReader::with_capacity(1 << 20, File::open(path).map_err(|e| format!("{path}: {e}"))?);
        let pool = rayon::ThreadPoolBuilder::new().num_threads(threads.max(1)).build()?;
        let (tx, rx) = sync_channel::<Result<Vec<u8>, String>>(2);
        let path = path.to_string();
        let thread = std::thread::spawn(move || {
            let mut next_batch = || -> Res<Option<Vec<u8>>> {
                let mut blocks = Vec::with_capacity(BATCH);
                while blocks.len() < BATCH {
                    let mut head = [0u8; 18];
                    match r.read_exact(&mut head[..12]) {
                        Ok(()) => {}
                        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                        Err(e) => return Err(e.into()),
                    }
                    let xlen = u16::from_le_bytes([head[10], head[11]]) as usize;
                    let mut block = head[..12].to_vec();
                    block.resize(12 + xlen, 0);
                    r.read_exact(&mut block[12..])?;
                    let size = block_size(&block).ok_or("not a BGZF block: a plain gzip member in the middle?")?;
                    if size < block.len() {
                        return Err("bad BGZF block size".into());
                    }
                    let start = block.len();
                    block.resize(size, 0);
                    r.read_exact(&mut block[start..])?;
                    blocks.push(block);
                }
                if blocks.is_empty() {
                    return Ok(None);
                }
                let parts: Vec<Res<Vec<u8>>> = pool.install(|| blocks.par_iter().map(|b| inflate(b)).collect());
                let mut out = Vec::with_capacity(parts.iter().map(|p| p.as_ref().map_or(0, Vec::len)).sum());
                for p in parts {
                    out.extend_from_slice(&p?);
                }
                Ok(Some(out))
            };
            loop {
                let msg = match next_batch() {
                    Ok(Some(b)) => Ok(b),
                    Ok(None) => return,
                    Err(e) => Err(format!("{path}: {e}")),
                };
                let stop = msg.is_err();
                if tx.send(msg).is_err() || stop {
                    return; // reader dropped, or error sent
                }
            }
        });
        Ok(BgzfReader { rx, buf: Vec::new(), pos: 0, thread: Some(thread) })
    }
}

impl Read for BgzfReader {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        while self.pos == self.buf.len() {
            match self.rx.recv() {
                Ok(Ok(b)) => (self.buf, self.pos) = (b, 0),
                Ok(Err(e)) => return Err(std::io::Error::other(e)),
                Err(_) => {
                    if let Some(t) = self.thread.take() {
                        let _ = t.join();
                    }
                    return Ok(0); // end of file
                }
            }
        }
        let n = out.len().min(self.buf.len() - self.pos);
        out[..n].copy_from_slice(&self.buf[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use std::io::Write;

    /// A BGZF file of `data`, in blocks of `chunk` bytes, with the standard EOF block.
    pub fn bgzf(data: &[u8], chunk: usize) -> Vec<u8> {
        let mut out = Vec::new();
        for part in data.chunks(chunk).chain(std::iter::once(&[][..])) {
            let mut e = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::default());
            e.write_all(part).unwrap();
            let cdata = e.finish().unwrap();
            let bsize = (18 + cdata.len() + 8 - 1) as u16;
            out.extend_from_slice(&[0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0]);
            out.extend_from_slice(&bsize.to_le_bytes());
            out.extend_from_slice(&cdata);
            let mut c = Crc::new();
            c.update(part);
            out.extend_from_slice(&c.sum().to_le_bytes());
            out.extend_from_slice(&(part.len() as u32).to_le_bytes());
        }
        out
    }

    #[test]
    fn round_trip_and_corruption() {
        let dir = std::env::temp_dir().join(format!("bycatch_bgzf_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("x.gz").to_str().unwrap().to_string();
        // more blocks than one batch, so batches and their order are exercised
        let data: Vec<u8> = (0..3_000_000u32).map(|i| b"ACGT\n"[(i.wrapping_mul(2654435761) >> 7) as usize % 5]).collect();
        let file = bgzf(&data, 10_000);
        std::fs::write(&path, &file).unwrap();
        assert!(is_bgzf(&path));
        let mut back = Vec::new();
        BgzfReader::open(&path, 3).unwrap().read_to_end(&mut back).unwrap();
        assert!(back == data, "round trip");
        // flip a byte inside the second block's deflate data
        let mut bad = file.clone();
        let second = 18 + u16::from_le_bytes([file[16], file[17]]) as usize + 1 - 18;
        bad[second + 30] ^= 0x55;
        std::fs::write(&path, &bad).unwrap();
        assert!(BgzfReader::open(&path, 2).unwrap().read_to_end(&mut Vec::new()).is_err());
        // plain gzip is not BGZF
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        e.write_all(&data[..1000]).unwrap();
        std::fs::write(&path, e.finish().unwrap()).unwrap();
        assert!(!is_bgzf(&path));
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
