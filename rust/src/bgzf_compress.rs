//! Parallel BGZF (block-gzip) compression.
//!
//! `bgzf_compress_file` reads an input file and writes a standard BGZF file
//! (gunzip-readable, and seekable/partitionable by `bgzf_scan` and `bgzf::Reader`).
//! BGZF blocks are independent gzip members, so blocks are compressed in parallel
//! across threads and concatenated **in input order** -- that ordering is the key
//! correctness invariant (a reordered block stream decompresses to scrambled bytes).
//!
//! Memory is bounded (~`threads * 64 KiB`), independent of file size: the input is
//! processed one super-chunk at a time, never slurped whole. The GIL is released
//! for the whole read/compress/write loop so this composes with Python-level
//! multiprocessing.

use bgzf::{CompressionLevel, Compressor, BGZF_BLOCK_SIZE};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::thread;

/// How many BGZF blocks each worker thread compresses per super-chunk. Larger
/// values amortize thread-spawn overhead; the super-chunk held in memory is
/// `threads * BLOCKS_PER_THREAD * BGZF_BLOCK_SIZE` bytes.
const BLOCKS_PER_THREAD: usize = 16;

/// Read up to `buf.len()` bytes, returning how many were read (0 at clean EOF).
fn read_upto<R: Read>(r: &mut R, buf: &mut [u8]) -> std::io::Result<usize> {
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

/// Compress one super-chunk: split `data` into BGZF_BLOCK_SIZE pieces, divide those
/// pieces into `nseg` contiguous segments compressed on separate threads, and
/// return the per-segment compressed bytes IN ORDER. Each thread owns its own
/// `Compressor` (not `Sync`) and compresses its blocks sequentially, so within and
/// across segments the block order matches the input.
fn compress_superchunk(data: &[u8], level: CompressionLevel, nseg: usize) -> Result<Vec<Vec<u8>>, String> {
    let n_blocks = data.len().div_ceil(BGZF_BLOCK_SIZE);
    let nseg = nseg.max(1).min(n_blocks.max(1));
    let blocks_per_seg = n_blocks.div_ceil(nseg);

    let results: Vec<Result<Vec<u8>, String>> = thread::scope(|s| {
        let mut handles = Vec::with_capacity(nseg);
        for seg in 0..nseg {
            let start_block = seg * blocks_per_seg;
            if start_block >= n_blocks {
                break;
            }
            let end_block = ((seg + 1) * blocks_per_seg).min(n_blocks);
            let byte_start = start_block * BGZF_BLOCK_SIZE;
            let byte_end = (end_block * BGZF_BLOCK_SIZE).min(data.len());
            let slice = &data[byte_start..byte_end];
            handles.push(s.spawn(move || -> Result<Vec<u8>, String> {
                let mut compressor = Compressor::new(level);
                let mut out = Vec::new();
                let mut block_buf = Vec::new();
                let mut pos = 0;
                while pos < slice.len() {
                    let end = (pos + BGZF_BLOCK_SIZE).min(slice.len());
                    compressor
                        .compress(&slice[pos..end], &mut block_buf)
                        .map_err(|e| format!("BGZF block compression failed: {e:?}"))?;
                    out.extend_from_slice(&block_buf);
                    pos = end;
                }
                Ok(out)
            }));
        }
        handles
            .into_iter()
            .map(|h| h.join().unwrap_or_else(|_| Err("BGZF compression thread panicked".to_string())))
            .collect()
    });

    let mut ordered = Vec::with_capacity(results.len());
    for r in results {
        ordered.push(r?);
    }
    Ok(ordered)
}

/// Compress `in_path` to `out_path` as standard BGZF using `threads` worker threads
/// (0 = autodetect). `level` is the deflate compression level (libdeflate range).
#[pyfunction]
#[pyo3(signature = (in_path, out_path, level=6, threads=0))]
pub fn bgzf_compress_file(
    py: Python<'_>,
    in_path: &str,
    out_path: &str,
    level: u8,
    threads: usize,
) -> PyResult<()> {
    let clevel = CompressionLevel::try_from(level)
        .map_err(|_| PyValueError::new_err(format!("invalid BGZF compression level {level}")))?;
    let nthreads = if threads == 0 {
        thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
    } else {
        threads
    };
    let superchunk_size = nthreads * BLOCKS_PER_THREAD * BGZF_BLOCK_SIZE;

    let result: Result<(), String> = py.allow_threads(|| {
        let mut reader = BufReader::new(File::open(in_path).map_err(|e| format!("opening {in_path}: {e}"))?);
        let mut writer = BufWriter::new(File::create(out_path).map_err(|e| format!("creating {out_path}: {e}"))?);
        let mut buf = vec![0u8; superchunk_size];
        loop {
            let got = read_upto(&mut reader, &mut buf).map_err(|e| format!("reading {in_path}: {e}"))?;
            if got == 0 {
                break;
            }
            for block in compress_superchunk(&buf[..got], clevel, nthreads)? {
                writer.write_all(&block).map_err(|e| format!("writing {out_path}: {e}"))?;
            }
        }
        // The BGZF EOF marker (empty block) is appended exactly once, after all data.
        let mut eof = Vec::new();
        Compressor::append_eof(&mut eof);
        writer.write_all(&eof).map_err(|e| format!("writing {out_path}: {e}"))?;
        writer.flush().map_err(|e| format!("flushing {out_path}: {e}"))?;
        Ok(())
    });

    result.map_err(PyValueError::new_err)
}
