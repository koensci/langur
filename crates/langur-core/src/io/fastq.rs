//! Streaming FASTQ reader backed by optional `igzip -dc` (auto-detected on Unix)
//! with seamless fallback to `flate2` (`zlib-ng`) + `needletail`, delivering
//! fixed-size flat-buffer read batches for deterministic `rayon` waves.

use crate::error::{Result, LangurError};
use flate2::read::MultiGzDecoder;
use needletail::{parse_fastx_reader, FastxReader};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;
#[cfg(unix)]
use std::path::PathBuf;
#[cfg(unix)]
use std::process::{Child, ChildStdout, Command, Stdio};
#[cfg(unix)]
use std::sync::OnceLock;

/// Compact owned FASTQ record inside a batch.
#[derive(Clone, Debug)]
pub struct FastqRec {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

/// Flat-buffer chunk of FASTQ reads eliminating per-read `Vec<u8>` heap allocations
/// (2 heap allocations per chunk instead of `N` allocations; `ids`/`quals` stay empty unless `--sam`).
#[derive(Clone, Debug)]
pub struct FastqChunk {
    pub seqs: Vec<u8>,
    pub seq_offsets: Vec<u32>,
    pub ids: Vec<u8>,
    pub id_offsets: Vec<u32>,
    pub quals: Vec<u8>,
    pub qual_offsets: Vec<u32>,
}

impl FastqChunk {
    pub fn with_capacity(max_records: usize, keep_id_qual: bool) -> Self {
        let mut seq_offsets = Vec::with_capacity(max_records + 1);
        seq_offsets.push(0);
        let (ids, id_offsets, quals, qual_offsets) = if keep_id_qual {
            let mut io = Vec::with_capacity(max_records + 1);
            io.push(0);
            let mut qo = Vec::with_capacity(max_records + 1);
            qo.push(0);
            (
                Vec::with_capacity(max_records * 32),
                io,
                Vec::with_capacity(max_records * 104),
                qo,
            )
        } else {
            (Vec::new(), Vec::new(), Vec::new(), Vec::new())
        };
        Self {
            seqs: Vec::with_capacity(max_records * 104),
            seq_offsets,
            ids,
            id_offsets,
            quals,
            qual_offsets,
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.seq_offsets.len().saturating_sub(1)
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[inline]
    pub fn seq(&self, idx: usize) -> &[u8] {
        let s = self.seq_offsets[idx] as usize;
        let e = self.seq_offsets[idx + 1] as usize;
        &self.seqs[s..e]
    }

    #[inline]
    pub fn id(&self, idx: usize) -> &[u8] {
        if self.id_offsets.len() <= idx + 1 {
            return &[];
        }
        let s = self.id_offsets[idx] as usize;
        let e = self.id_offsets[idx + 1] as usize;
        &self.ids[s..e]
    }

    #[inline]
    pub fn qual(&self, idx: usize) -> &[u8] {
        if self.qual_offsets.len() <= idx + 1 {
            return &[];
        }
        let s = self.qual_offsets[idx] as usize;
        let e = self.qual_offsets[idx + 1] as usize;
        &self.quals[s..e]
    }
}

#[cfg(unix)]
fn detect_igzip_path() -> Option<&'static PathBuf> {
    static IGZIP_PATH: OnceLock<Option<PathBuf>> = OnceLock::new();
    IGZIP_PATH
        .get_or_init(|| {
            for cand in ["igzip", "/opt/homebrew/bin/igzip", "/usr/local/bin/igzip"] {
                if let Ok(status) = Command::new(cand)
                    .arg("--version")
                    .stdin(Stdio::null())
                    .stdout(Stdio::null())
                    .stderr(Stdio::null())
                    .status()
                {
                    if status.success() {
                        return Some(PathBuf::from(cand));
                    }
                }
            }
            None
        })
        .as_ref()
}

#[cfg(unix)]
struct IgzipReader {
    stdout: BufReader<ChildStdout>,
    child: Option<Child>,
}

#[cfg(unix)]
impl Read for IgzipReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.stdout.read(buf)?;
        if n == 0 {
            if let Some(mut child) = self.child.take() {
                let status = child.wait()?;
                if !status.success() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("igzip decompression failed with status {status}"),
                    ));
                }
            }
        }
        Ok(n)
    }
}

#[cfg(unix)]
impl Drop for IgzipReader {
    fn drop(&mut self) {
        if let Some(mut child) = self.child.take() {
            let _ = child.kill();
            let _ = child.wait();
        }
    }
}

fn open_raw_stream(path: &Path, force_gz: bool) -> Result<BufReader<Box<dyn Read + Send>>> {
    let mut f = File::open(path).map_err(|source| LangurError::IoPath {
        path: path.to_path_buf(),
        source,
    })?;
    let mut magic = [0u8; 2];
    let n = f.read(&mut magic)?;
    drop(f);

    let is_gz = force_gz || (n == 2 && magic == [0x1F, 0x8B]);

    if is_gz {
        #[cfg(unix)]
        if let Some(igzip_bin) = detect_igzip_path() {
            if let Ok(mut child) = Command::new(igzip_bin)
                .arg("-dc")
                .arg(path)
                .stdin(Stdio::null())
                .stdout(Stdio::piped())
                .stderr(Stdio::null())
                .spawn()
            {
                if let Some(stdout) = child.stdout.take() {
                    let reader = IgzipReader {
                        stdout: BufReader::with_capacity(512 * 1024, stdout),
                        child: Some(child),
                    };
                    return Ok(BufReader::with_capacity(512 * 1024, Box::new(reader)));
                }
                let _ = child.kill();
                let _ = child.wait();
            }
        }

        // Fallback: flate2 with `zlib-ng` SIMD decompression
        let f = File::open(path).map_err(|source| LangurError::IoPath {
            path: path.to_path_buf(),
            source,
        })?;
        let dec = MultiGzDecoder::new(BufReader::with_capacity(512 * 1024, f));
        Ok(BufReader::with_capacity(512 * 1024, Box::new(dec)))
    } else {
        let f = File::open(path).map_err(|source| LangurError::IoPath {
            path: path.to_path_buf(),
            source,
        })?;
        Ok(BufReader::with_capacity(512 * 1024, Box::new(f)))
    }
}

fn open_fastx_reader(path: &Path, force_gz: bool) -> Result<Box<dyn FastxReader>> {
    let stream = open_raw_stream(path, force_gz)?;
    parse_fastx_reader(stream)
        .map_err(|e| LangurError::Fasta(format!("FASTQ open error {}: {e}", path.display())))
}

/// Zero-copy 4-line FASTQ block scanner directly over `BufReader` (skips header/plus/qual
/// copies completely when `!keep_id_qual`, validating `@`, `+`, and `seq.len() == qual.len()`).
struct FastFastqReader {
    buf: BufReader<Box<dyn Read + Send>>,
}

impl FastFastqReader {
    /// Consume one line, returning `(trimmed_len, first_byte)` without copying bytes anywhere.
    #[inline]
    fn skip_line(&mut self) -> Result<Option<(usize, u8)>> {
        let mut total = 0usize;
        let mut first = 0u8;
        let mut saw_any = false;
        let mut last_was_cr = false;
        loop {
            let avail = self.buf.fill_buf()?;
            if avail.is_empty() {
                if !saw_any {
                    return Ok(None);
                }
                if last_was_cr {
                    total = total.saturating_sub(1);
                }
                return Ok(Some((total, first)));
            }
            if !saw_any {
                first = avail[0];
                saw_any = true;
            }
            if let Some(pos) = avail.iter().position(|&b| b == b'\n') {
                let mut line_len = pos;
                if line_len > 0 && avail[line_len - 1] == b'\r' {
                    line_len -= 1;
                } else if line_len == 0 && last_was_cr {
                    total = total.saturating_sub(1);
                }
                total += line_len;
                self.buf.consume(pos + 1);
                return Ok(Some((total, first)));
            }
            let len = avail.len();
            last_was_cr = avail[len - 1] == b'\r';
            total += len;
            self.buf.consume(len);
        }
    }

    /// Append one line (stripped of `\r\n`) directly into `dst`, returning `Some(len_appended)` or `None` on EOF.
    #[inline]
    fn append_line(&mut self, dst: &mut Vec<u8>) -> Result<Option<usize>> {
        let start_len = dst.len();
        let mut saw_any = false;
        loop {
            let avail = self.buf.fill_buf()?;
            if avail.is_empty() {
                if !saw_any {
                    return Ok(None);
                }
                if dst.len() > start_len && dst[dst.len() - 1] == b'\r' {
                    dst.pop();
                }
                return Ok(Some(dst.len() - start_len));
            }
            saw_any = true;
            if let Some(pos) = avail.iter().position(|&b| b == b'\n') {
                let slice = &avail[..pos];
                dst.extend_from_slice(slice);
                if dst.len() > start_len && dst[dst.len() - 1] == b'\r' {
                    dst.pop();
                }
                self.buf.consume(pos + 1);
                return Ok(Some(dst.len() - start_len));
            }
            dst.extend_from_slice(avail);
            let len = avail.len();
            self.buf.consume(len);
        }
    }

    fn next_chunk(&mut self, max_records: usize, keep_id_qual: bool) -> Result<FastqChunk> {
        let mut chunk = FastqChunk::with_capacity(max_records, keep_id_qual);
        while chunk.len() < max_records {
            // Line 1: @header
            if keep_id_qual {
                let id_start = chunk.ids.len();
                let Some(id_len) = self.append_line(&mut chunk.ids)? else {
                    break;
                };
                if id_len == 0 || chunk.ids[id_start] != b'@' {
                    return Err(LangurError::Fasta(
                        "FASTQ parse error: expected '@' header line".into(),
                    ));
                }
                // Remove leading '@' and keep up to first whitespace (matching needletail `rec.id()`)
                chunk.ids.remove(id_start);
                if let Some(ws) = chunk.ids[id_start..]
                    .iter()
                    .position(|&b| b == b' ' || b == b'\t')
                {
                    chunk.ids.truncate(id_start + ws);
                }
                chunk.id_offsets.push(chunk.ids.len() as u32);
            } else {
                let Some((hdr_len, first)) = self.skip_line()? else {
                    break;
                };
                if hdr_len == 0 || first != b'@' {
                    return Err(LangurError::Fasta(
                        "FASTQ parse error: expected '@' header line".into(),
                    ));
                }
            }

            // Line 2: sequence
            let Some(seq_len) = self.append_line(&mut chunk.seqs)? else {
                return Err(LangurError::Fasta(
                    "FASTQ parse error: truncated record before sequence line".into(),
                ));
            };
            chunk.seq_offsets.push(chunk.seqs.len() as u32);

            // Line 3: '+' separator
            let Some((plus_len, plus_first)) = self.skip_line()? else {
                return Err(LangurError::Fasta(
                    "FASTQ parse error: truncated record before '+' line".into(),
                ));
            };
            if plus_len == 0 || plus_first != b'+' {
                return Err(LangurError::Fasta(
                    "FASTQ parse error: expected '+' separator line".into(),
                ));
            }

            // Line 4: quality
            let qual_len = if keep_id_qual {
                let Some(ql) = self.append_line(&mut chunk.quals)? else {
                    return Err(LangurError::Fasta(
                        "FASTQ parse error: truncated record before quality line".into(),
                    ));
                };
                chunk.qual_offsets.push(chunk.quals.len() as u32);
                ql
            } else {
                let Some((ql, _)) = self.skip_line()? else {
                    return Err(LangurError::Fasta(
                        "FASTQ parse error: truncated record before quality line".into(),
                    ));
                };
                ql
            };
            if qual_len != seq_len {
                return Err(LangurError::Fasta(format!(
                    "FASTQ parse error: sequence length ({seq_len}) != quality length ({qual_len})"
                )));
            }
        }
        Ok(chunk)
    }
}

enum BatchReaderInner {
    FastFastq(FastFastqReader),
    Needletail(Box<dyn FastxReader>),
}

/// Single-end streaming batch FASTQ/FASTA reader.
pub struct FastqBatchReader {
    inner: BatchReaderInner,
    keep_id_qual: bool,
}

impl FastqBatchReader {
    pub fn open(path: &Path, force_gz: bool, keep_id_qual: bool) -> Result<Self> {
        let mut stream = open_raw_stream(path, force_gz)?;
        let first_byte = {
            let peek = stream.fill_buf()?;
            peek.first().copied()
        };
        let inner = if first_byte == Some(b'@') {
            BatchReaderInner::FastFastq(FastFastqReader { buf: stream })
        } else {
            let reader = parse_fastx_reader(stream).map_err(|e| {
                LangurError::Fasta(format!("FASTQ open error {}: {e}", path.display()))
            })?;
            BatchReaderInner::Needletail(reader)
        };
        Ok(Self {
            inner,
            keep_id_qual,
        })
    }

    pub fn next_chunk(&mut self, max_records: usize) -> Result<FastqChunk> {
        match &mut self.inner {
            BatchReaderInner::FastFastq(r) => r.next_chunk(max_records, self.keep_id_qual),
            BatchReaderInner::Needletail(reader) => {
                let mut chunk = FastqChunk::with_capacity(max_records, self.keep_id_qual);
                while chunk.len() < max_records {
                    let Some(rec_res) = reader.next() else {
                        break;
                    };
                    let rec = rec_res
                        .map_err(|e| LangurError::Fasta(format!("FASTQ parse error: {e}")))?;
                    chunk.seqs.extend_from_slice(&rec.seq());
                    chunk.seq_offsets.push(chunk.seqs.len() as u32);
                    if self.keep_id_qual {
                        chunk.ids.extend_from_slice(rec.id());
                        chunk.id_offsets.push(chunk.ids.len() as u32);
                        if let Some(q) = rec.qual() {
                            chunk.quals.extend_from_slice(q);
                        }
                        chunk.qual_offsets.push(chunk.quals.len() as u32);
                    }
                }
                Ok(chunk)
            }
        }
    }

    pub fn next_batch(&mut self, max_records: usize) -> Result<Vec<FastqRec>> {
        let chunk = self.next_chunk(max_records)?;
        let mut batch = Vec::with_capacity(chunk.len());
        for i in 0..chunk.len() {
            batch.push(FastqRec {
                id: chunk.id(i).to_vec(),
                seq: chunk.seq(i).to_vec(),
                qual: chunk.qual(i).to_vec(),
            });
        }
        Ok(batch)
    }
}

/// Synchronized paired-end streaming batch FASTQ reader (R1 + R2).
pub struct PairedFastqBatchReader {
    r1: Box<dyn FastxReader>,
    r2: Box<dyn FastxReader>,
    keep_id_qual: bool,
}

impl PairedFastqBatchReader {
    pub fn open(r1_path: &Path, r2_path: &Path, force_gz: bool, keep_id_qual: bool) -> Result<Self> {
        Ok(Self {
            r1: open_fastx_reader(r1_path, force_gz)?,
            r2: open_fastx_reader(r2_path, force_gz)?,
            keep_id_qual,
        })
    }

    pub fn next_batch(&mut self, max_pairs: usize) -> Result<Vec<(FastqRec, FastqRec)>> {
        let mut batch = Vec::with_capacity(max_pairs);
        while batch.len() < max_pairs {
            match (self.r1.next(), self.r2.next()) {
                (Some(res1), Some(res2)) => {
                    let rec1 = res1
                        .map_err(|e| LangurError::Fasta(format!("R1 FASTQ parse error: {e}")))?;
                    let rec2 = res2
                        .map_err(|e| LangurError::Fasta(format!("R2 FASTQ parse error: {e}")))?;
                    let f1 = FastqRec {
                        id: if self.keep_id_qual {
                            rec1.id().to_vec()
                        } else {
                            Vec::new()
                        },
                        seq: rec1.seq().into_owned(),
                        qual: if self.keep_id_qual {
                            rec1.qual().map(|q| q.to_vec()).unwrap_or_default()
                        } else {
                            Vec::new()
                        },
                    };
                    let f2 = FastqRec {
                        id: if self.keep_id_qual {
                            rec2.id().to_vec()
                        } else {
                            Vec::new()
                        },
                        seq: rec2.seq().into_owned(),
                        qual: if self.keep_id_qual {
                            rec2.qual().map(|q| q.to_vec()).unwrap_or_default()
                        } else {
                            Vec::new()
                        },
                    };
                    batch.push((f1, f2));
                }
                (None, None) => break,
                _ => {
                    return Err(LangurError::Fasta(
                        "paired FASTQ files have different numbers of reads".into(),
                    ));
                }
            }
        }
        Ok(batch)
    }
}
