//! BGZF (block-gzip) offset scanner.
//!
//! Produces the same `(virtual_offset, cds_count)` pairs as the Python
//! `utils.get_genbank_offsets`/`get_fasta_offsets` BgzfReader scanners, but
//! decompresses blocks in Rust. Virtual offsets use the standard
//! `block_start << 16 | within_block` encoding (matching Bio.bgzf and
//! `parser::open_at`), so the existing seek/partition path consumes them
//! unchanged.
//!
//! Each BGZF block is an independent gzip member whose total size is carried in
//! the gzip extra field's `BC` subfield (htslib/bgzip convention). We read the
//! header, slice out the member, and decompress it with libdeflater (already in
//! the dependency tree via the `bgzf` crate). Any structural surprise raises a
//! Python error, which the caller treats as "fall back to the Python scanner".

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::fs::File;
use std::io::{BufReader, Read};

use crate::locus_unit_is_aa;

/// Cap on bytes retained from the start of each line. Enough to classify a line
/// ("LOCUS "/"     CDS ") and reach a LOCUS line's molecule-unit token (`parts[3]`)
/// for any standard-width LOCUS line, without buffering long ORIGIN/qualifier lines.
const PREFIX_CAP: usize = 256;

/// Read up to `buf.len()` bytes, returning how many were actually read (0 at clean
/// EOF). Unlike `read_exact`, a short read is reported rather than erroring, so the
/// caller can distinguish a clean end-of-file from a truncated block.
fn read_full<R: Read>(r: &mut R, buf: &mut [u8]) -> std::io::Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match r.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        }
    }
    Ok(filled)
}

/// Find the BGZF `BC` subfield (SI1='B'=66, SI2='C'=67) in a gzip extra field and
/// return its BSIZE value (total block size minus one).
fn bc_bsize(extra: &[u8]) -> Option<u16> {
    let mut i = 0;
    while i + 4 <= extra.len() {
        let si1 = extra[i];
        let si2 = extra[i + 1];
        let slen = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
        if si1 == 66 && si2 == 67 && slen >= 2 && i + 6 <= extra.len() {
            return Some(u16::from_le_bytes([extra[i + 4], extra[i + 5]]));
        }
        i += 4 + slen;
    }
    None
}

/// Walk the BGZF blocks of `path`, invoking `f(block_start, &decompressed)` for
/// each non-empty block, where `block_start` is the block's compressed-file offset.
fn for_each_block<F>(path: &str, mut f: F) -> PyResult<()>
where
    F: FnMut(u64, &[u8]),
{
    let mut reader = BufReader::new(File::open(path)?);
    let mut dec = libdeflater::Decompressor::new();
    let mut block_start: u64 = 0;
    let mut header = [0u8; 12];
    let mut block: Vec<u8> = Vec::new();
    let mut out: Vec<u8> = Vec::new();

    loop {
        // Fixed gzip header through XLEN. A clean EOF here ends the file.
        match read_full(&mut reader, &mut header)? {
            0 => break,
            12 => {}
            _ => return Err(PyValueError::new_err("truncated BGZF header")),
        }
        if header[0] != 0x1f || header[1] != 0x8b || header[2] != 8 || (header[3] & 0x04) == 0 {
            return Err(PyValueError::new_err(
                "not a BGZF block (gzip magic or FEXTRA flag missing)",
            ));
        }
        let xlen = u16::from_le_bytes([header[10], header[11]]) as usize;
        let mut extra = vec![0u8; xlen];
        if read_full(&mut reader, &mut extra)? != xlen {
            return Err(PyValueError::new_err("truncated BGZF extra field"));
        }
        let bsize = bc_bsize(&extra)
            .ok_or_else(|| PyValueError::new_err("BGZF block missing BC subfield"))?;
        let block_size = bsize as usize + 1;
        let consumed = 12 + xlen;
        if block_size < consumed + 8 {
            return Err(PyValueError::new_err("BGZF block size too small"));
        }

        // Reassemble the full gzip member (header + extra + body) for decompression.
        block.clear();
        block.extend_from_slice(&header);
        block.extend_from_slice(&extra);
        let body = block_size - consumed;
        let base = block.len();
        block.resize(base + body, 0);
        if read_full(&mut reader, &mut block[base..])? != body {
            return Err(PyValueError::new_err("truncated BGZF block body"));
        }

        // ISIZE (decompressed size) is the last 4 bytes of the member.
        let isize_dec = u32::from_le_bytes([
            block[block_size - 4],
            block[block_size - 3],
            block[block_size - 2],
            block[block_size - 1],
        ]) as usize;

        if isize_dec == 0 {
            // BGZF EOF marker (or an empty member): nothing to scan.
            block_start += block_size as u64;
            continue;
        }

        out.clear();
        out.resize(isize_dec, 0);
        dec.gzip_decompress(&block, &mut out)
            .map_err(|e| PyValueError::new_err(format!("BGZF block decompress failed: {e:?}")))?;
        f(block_start, &out);
        block_start += block_size as u64;
    }
    Ok(())
}

/// GenBank line scanner with state carried across block boundaries. Mirrors the
/// state machine in `genbank_offsets` (lib.rs) but keyed on virtual offsets.
#[derive(Default)]
struct GbScan {
    out: Vec<(u64, u64)>,
    last_offset: u64,
    cdss: u64,
    is_aa: bool,
    locus_seen: bool,
    line_pending: bool,
    line_voffset: u64,
    prefix: Vec<u8>,
}

impl GbScan {
    fn feed(&mut self, block_start: u64, data: &[u8]) {
        for (i, &b) in data.iter().enumerate() {
            if !self.line_pending {
                self.line_voffset = (block_start << 16) | (i as u64);
                self.prefix.clear();
                self.line_pending = true;
            }
            if b == b'\n' {
                self.classify();
                self.line_pending = false;
            } else if self.prefix.len() < PREFIX_CAP {
                self.prefix.push(b);
            }
        }
    }

    fn classify(&mut self) {
        if self.prefix.starts_with(b"LOCUS ") {
            if self.locus_seen {
                self.out.push((self.last_offset, self.cdss));
            }
            self.locus_seen = true;
            self.is_aa = locus_unit_is_aa(&self.prefix);
            self.cdss = if self.is_aa { 1 } else { 0 };
            self.last_offset = self.line_voffset;
        } else if self.prefix.starts_with(b"     CDS ") && !self.is_aa {
            self.cdss += 1;
        }
    }

    fn finish(mut self) -> Vec<(u64, u64)> {
        if self.line_pending {
            self.classify();
        }
        if self.locus_seen {
            self.out.push((self.last_offset, self.cdss));
        }
        self.out
    }
}

/// FASTA line scanner: record the virtual offset of every line beginning with '>'.
#[derive(Default)]
struct FaScan {
    out: Vec<(u64, u64)>,
    line_pending: bool,
}

impl FaScan {
    fn feed(&mut self, block_start: u64, data: &[u8]) {
        for (i, &b) in data.iter().enumerate() {
            if !self.line_pending {
                self.line_pending = true;
                if b == b'>' {
                    self.out.push(((block_start << 16) | (i as u64), 1));
                }
            }
            if b == b'\n' {
                self.line_pending = false;
            }
        }
    }
}

/// Scan a BGZF-compressed GenBank file, returning `(virtual_offset, cds_count)`.
#[pyfunction]
pub fn genbank_offsets_bgzf(path: &str) -> PyResult<Vec<(u64, u64)>> {
    let mut scan = GbScan::default();
    for_each_block(path, |bs, data| scan.feed(bs, data))?;
    Ok(scan.finish())
}

/// Scan a BGZF-compressed FASTA file, returning `(virtual_offset, 1)` per record.
#[pyfunction]
pub fn fasta_offsets_bgzf(path: &str) -> PyResult<Vec<(u64, u64)>> {
    let mut scan = FaScan::default();
    for_each_block(path, |bs, data| scan.feed(bs, data))?;
    Ok(scan.out)
}
