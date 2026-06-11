//! Native acceleration for Domainator's GenBank/FASTA hot path.
//!
//! This crate is an OPTIONAL extension (`domainator._gbfast`). The Python code
//! always works without it (pure-Python / gb-io fallbacks), so anything here is a
//! speed optimization, never a correctness dependency.
//!
//! Increment 2 (this file): the offset scanner that replaces the pure-Python
//! `get_genbank_offsets` / `get_fasta_offsets` byte iteration. It returns the
//! same `(byte_offset, cds_count)` pairs, so the partition scheme is unchanged.
//! Only uncompressed files are handled here; compressed (BGZF) inputs stay on the
//! Python virtual-offset scanner.

use pyo3::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};

mod parser;

/// Return the 4th whitespace-delimited token of a LOCUS line, if present.
/// Mirrors the Python scanner's `parts[3]` (the molecule-unit field: "aa"/"bp").
fn locus_unit_is_aa(line: &[u8]) -> bool {
    let mut tokens = line.split(|&b| b == b' ' || b == b'\t' || b == b'\r' || b == b'\n').filter(|t| !t.is_empty());
    // tokens: [LOCUS, name, length, unit, ...]
    matches!(tokens.nth(3), Some(b"aa"))
}

/// Scan a GenBank file, returning (byte_offset, cds_count) for each record.
///
/// Equivalent to utils.get_genbank_offsets: cds_count is 1 for amino-acid (aa)
/// records and the number of `     CDS ` feature lines for nucleotide (bp) records.
#[pyfunction]
fn genbank_offsets(path: &str) -> PyResult<Vec<(u64, u64)>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut out: Vec<(u64, u64)> = Vec::new();

    let mut offset: u64 = 0;
    let mut last_offset: u64 = 0;
    let mut cdss: u64 = 0;
    let mut is_aa = false;
    let mut locus_seen = false;
    let mut line: Vec<u8> = Vec::new();

    // Advance to the first LOCUS line.
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        if line.starts_with(b"LOCUS ") {
            locus_seen = true;
            is_aa = locus_unit_is_aa(&line);
            cdss = if is_aa { 1 } else { 0 };
            last_offset = offset;
            offset += n as u64;
            break;
        }
        offset += n as u64;
    }

    // Remaining records.
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        if line.starts_with(b"LOCUS ") {
            out.push((last_offset, cdss));
            is_aa = locus_unit_is_aa(&line);
            cdss = if is_aa { 1 } else { 0 };
            last_offset = offset;
        } else if line.starts_with(b"     CDS ") {
            if !is_aa {
                cdss += 1;
            }
        }
        offset += n as u64;
    }

    if locus_seen {
        out.push((last_offset, cdss));
    }
    Ok(out)
}

/// Scan a FASTA file, returning (byte_offset, 1) for each record.
/// Equivalent to utils.get_fasta_offsets.
#[pyfunction]
fn fasta_offsets(path: &str) -> PyResult<Vec<(u64, u64)>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut out: Vec<(u64, u64)> = Vec::new();

    let mut offset: u64 = 0;
    let mut line: Vec<u8> = Vec::new();
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        if line.first() == Some(&b'>') {
            out.push((offset, 1));
        }
        offset += n as u64;
    }
    Ok(out)
}

#[pymodule]
fn _gbfast(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(genbank_offsets, m)?)?;
    m.add_function(wrap_pyfunction!(fasta_offsets, m)?)?;
    m.add_function(wrap_pyfunction!(parser::parse_lean, m)?)?;
    m.add_function(wrap_pyfunction!(parser::parse_lean_search, m)?)?;
    m.add_function(wrap_pyfunction!(parser::parse_fasta_search, m)?)?;
    m.add_class::<parser::LeanSearchContig>()?;
    m.add_class::<parser::LeanFastaContig>()?;
    Ok(())
}
