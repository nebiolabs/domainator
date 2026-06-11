//! Native GenBank record parser (increments 3a/3b).
//!
//! Parses a partition of a GenBank file with the `gb-io` crate. Two modes:
//!
//! * `parse_lean` (3a): builds the full per-feature `LeanFeature` objects in Rust
//!   and returns header tuples. Used by `iter_lean_genbank_native`.
//!
//! * `parse_lean_search` (3b): returns `LeanSearchContig` pyclasses that HOLD the
//!   parsed record in Rust and expose only what `domainate`'s scan needs
//!   (molecule_type, seq, cds peptides). No per-feature Python objects are built
//!   for the bulk of records; the full record is materialized (`.materialize()`)
//!   only when a contig produces a hit. This removes the per-feature
//!   object-construction floor from the ~99% of records that don't hit.
//!
//! Records whose LOCUS line gb-io could not parse (molecule_type None), parse
//! errors, unmodelled location kinds, or (for search mode) a CDS lacking a
//! `/translation` all stop the run early so the Python side falls back to the
//! Biopython parser for the remaining records.

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};

use gb_io::reader::SeqReader;
use gb_io::seq::{Date, Location, Seq, Topology};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};

/// Open `path` and position a reader at `seek_to`, decompressing BGZF on the fly.
///
/// For uncompressed input, `seek_to` is a plain byte offset. For BGZF, `seek_to`
/// is a 64-bit virtual offset (`block_start << 16 | within_block`, the same
/// encoding Bio.bgzf produces): we seek the underlying file to the block boundary
/// (each BGZF block is an independent gzip member), decompress from there, and
/// discard the within-block bytes. The returned reader yields the decompressed
/// stream from the record boundary onward.
fn open_at(path: &str, seek_to: u64, compressed: bool) -> io::Result<Box<dyn Read>> {
    let mut file = File::open(path)?;
    if !compressed {
        if seek_to > 0 {
            file.seek(SeekFrom::Start(seek_to))?;
        }
        return Ok(Box::new(BufReader::new(file)));
    }
    let block_start = seek_to >> 16;
    let within = (seek_to & 0xFFFF) as usize;
    file.seek(SeekFrom::Start(block_start))?;
    let mut reader = bgzf::Reader::new(file);
    let mut remaining = within;
    let mut buf = [0u8; 8192];
    while remaining > 0 {
        let want = remaining.min(buf.len());
        let n = reader.read(&mut buf[..want])?;
        if n == 0 {
            break;
        }
        remaining -= n;
    }
    Ok(Box::new(reader))
}

const MONTHS: [&str; 12] = [
    "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
];

fn fmt_date(d: &Date) -> String {
    let m = d.month();
    let mon = if (1..=12).contains(&m) { MONTHS[(m - 1) as usize] } else { "JAN" };
    format!("{:02}-{}-{}", d.day(), mon, d.year())
}

/// Replicate lean_record._join_multiline_qualifier + _normalize_qualifier_value.
fn normalize_qualifier(key: &str, value: &str) -> String {
    let joined = if value.contains('\n') {
        if key == "translation" {
            value.replace('\n', "")
        } else {
            value.split('\n').map(|p| p.trim()).collect::<Vec<_>>().join(" ")
        }
    } else {
        value.to_string()
    };
    let stripped = joined.trim_end();
    let bytes = stripped.as_bytes();
    if bytes.len() >= 2 && bytes[0] == b'"' && bytes[bytes.len() - 1] == b'"' {
        stripped[1..stripped.len() - 1].to_string()
    } else {
        joined
    }
}

type Part = (i64, i64, i32, bool, bool); // (start, end, strand, before, after)

/// Convert a gb-io Location to lean (operator, between, parts). None for kinds we
/// don't model (OneOf/Bond/External/Gap).
fn lean_location(loc: &Location, strand: i32) -> Option<(Option<&'static str>, bool, Vec<Part>)> {
    match loc {
        Location::Range((start, before), (end, after)) => {
            Some((None, false, vec![(*start, *end, strand, before.0, after.0)]))
        }
        Location::Between(_s, e) => Some((None, true, vec![(*e, *e, strand, false, false)])),
        Location::Complement(inner) => lean_location(inner, -strand),
        Location::Join(subs) | Location::Order(subs) => {
            let op = if matches!(loc, Location::Join(_)) { "join" } else { "order" };
            let mut groups: Vec<Vec<Part>> = Vec::with_capacity(subs.len());
            for sub in subs {
                let (_o, _b, parts) = lean_location(sub, strand)?;
                groups.push(parts);
            }
            if strand == -1 {
                groups.reverse();
            }
            Some((Some(op), false, groups.into_iter().flatten().collect()))
        }
        _ => None,
    }
}

/// Translate one codon with the NCBI standard table (table 1). Any codon with an
/// ambiguous/unknown base -> 'X'; '*' is a stop. Used only to produce SEARCH
/// peptides for pyhmmer (the materialized record's /translation comes from the
/// file or, when absent, from Biopython at hit time). Codon-by-codon so the
/// peptide length matches Biopython's, which keeps hit coordinates consistent.
fn translate_codon(c0: u8, c1: u8, c2: u8) -> u8 {
    fn idx(b: u8) -> Option<usize> {
        match b.to_ascii_uppercase() {
            b'T' | b'U' => Some(0),
            b'C' => Some(1),
            b'A' => Some(2),
            b'G' => Some(3),
            _ => None,
        }
    }
    // TCAG-ordered standard genetic code.
    const TABLE: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    match (idx(c0), idx(c1), idx(c2)) {
        (Some(a), Some(b), Some(c)) => TABLE[a * 16 + b * 4 + c],
        _ => b'X',
    }
}

/// IUPAC reverse-complement, preserving case, matching Biopython's complement.
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T', b'T' => b'A', b'U' => b'A', b'G' => b'C', b'C' => b'G',
            b'a' => b't', b't' => b'a', b'u' => b'a', b'g' => b'c', b'c' => b'g',
            b'R' => b'Y', b'Y' => b'R', b'S' => b'S', b'W' => b'W', b'K' => b'M', b'M' => b'K',
            b'r' => b'y', b'y' => b'r', b's' => b's', b'w' => b'w', b'k' => b'm', b'm' => b'k',
            b'B' => b'V', b'V' => b'B', b'D' => b'H', b'H' => b'D',
            b'b' => b'v', b'v' => b'b', b'd' => b'h', b'h' => b'd',
            b'N' => b'N', b'n' => b'n',
            other => other,
        })
        .collect()
}

fn translate_dna(seq: &[u8]) -> String {
    let mut out = Vec::with_capacity(seq.len() / 3);
    let mut i = 0;
    while i + 3 <= seq.len() {
        out.push(translate_codon(seq[i], seq[i + 1], seq[i + 2]));
        i += 3;
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// Extract the spliced nucleotide sequence for a feature location from the contig,
/// honoring per-part strand (reverse-complementing minus parts) and join order.
/// Returns None for unmodelled locations.
fn extract_location_seq(loc: &Location, contig: &[u8]) -> Option<Vec<u8>> {
    let (_op, _between, parts) = lean_location(loc, 1)?;
    let mut out = Vec::new();
    for (start, end, strand, _b, _a) in parts {
        let s = start.max(0) as usize;
        let e = (end.max(0) as usize).min(contig.len());
        if s >= e {
            continue;
        }
        let mut piece = contig[s..e].to_vec();
        if strand == -1 {
            piece = revcomp(&piece);
        }
        out.extend_from_slice(&piece);
    }
    Some(out)
}

fn is_searchable_cds(feature: &gb_io::seq::Feature) -> bool {
    feature.kind.as_ref() == "CDS"
        && !feature.qualifiers.iter().any(|(k, _)| k == "pseudo" || k == "pseudogene")
}

fn translation_value(feature: &gb_io::seq::Feature) -> Option<&str> {
    feature
        .qualifiers
        .iter()
        .find(|(k, _)| k == "translation")
        .and_then(|(_, v)| v.as_deref())
        .filter(|s| !s.is_empty())
}

fn build_feature<'py>(
    py: Python<'py>,
    lean_feature_cls: &Bound<'py, PyAny>,
    feature: &gb_io::seq::Feature,
) -> PyResult<Option<Bound<'py, PyAny>>> {
    let (operator, between, parts) = match lean_location(&feature.location, 1) {
        Some(v) => v,
        None => return Ok(None),
    };
    let part_objs: Vec<Bound<'py, PyTuple>> = parts
        .iter()
        .map(|(s, e, st, b, a)| {
            PyTuple::new_bound(py, [
                (*s).into_py(py), (*e).into_py(py), (*st).into_py(py), (*b).into_py(py), (*a).into_py(py),
            ])
        })
        .collect();
    let parts_tuple = PyTuple::new_bound(py, &part_objs);

    let quals = PyDict::new_bound(py);
    for (key, value_opt) in &feature.qualifiers {
        let value = match value_opt {
            Some(v) => normalize_qualifier(key, v),
            None => String::new(),
        };
        match quals.get_item(key.as_ref())? {
            Some(existing) => existing.downcast::<PyList>()?.append(value)?,
            None => quals.set_item(key.as_ref(), PyList::new_bound(py, [value]))?,
        }
    }
    let op_obj = match operator {
        Some(s) => s.into_py(py),
        None => py.None(),
    };
    Ok(Some(lean_feature_cls.call1((feature.kind.as_ref(), parts_tuple, op_obj, between, quals))?))
}

/// Build the 14-field header tuple + features list for a record, dropping any
/// features whose kind is in `dropped`. Returns None if a feature location is
/// unmodelled (caller falls back).
fn record_fields<'py>(
    py: Python<'py>,
    seq: &Seq,
    dropped: &HashSet<String>,
) -> PyResult<Option<Bound<'py, PyTuple>>> {
    let lean_mod = py.import_bound("domainator.lean_record")?;
    let lean_feature_cls = lean_mod.getattr("LeanFeature")?;

    let features = PyList::empty_bound(py);
    for feature in &seq.features {
        if !dropped.is_empty() && dropped.contains(feature.kind.as_ref()) {
            continue;
        }
        match build_feature(py, &lean_feature_cls, feature)? {
            Some(f) => features.append(f)?,
            None => return Ok(None),
        }
    }
    Ok(Some(header_tuple(py, seq, features)))
}

fn header_tuple<'py>(py: Python<'py>, seq: &Seq, features: Bound<'py, PyList>) -> Bound<'py, PyTuple> {
    let circular = matches!(seq.topology, Topology::Circular);
    let date_str = seq.date.as_ref().map(fmt_date);
    let (source_name, source_organism) = match &seq.source {
        Some(s) => (Some(s.source.clone()), s.organism.clone()),
        None => (None, None),
    };
    let seq_str = String::from_utf8_lossy(&seq.seq.to_ascii_uppercase()).into_owned();
    PyTuple::new_bound(py, [
        seq.name.clone().into_py(py),
        seq.accession.clone().into_py(py),
        seq.version.clone().into_py(py),
        seq.definition.clone().into_py(py),
        seq.molecule_type.clone().into_py(py),
        circular.into_py(py),
        seq.division.clone().into_py(py),
        date_str.into_py(py),
        seq.keywords.clone().into_py(py),
        source_name.into_py(py),
        source_organism.into_py(py),
        seq.dblink.clone().into_py(py),
        seq_str.into_py(py),
        features.into_py(py),
    ])
}

#[pyfunction]
#[pyo3(signature = (path, seek_to, max_recs, compressed=false))]
pub fn parse_lean<'py>(
    py: Python<'py>,
    path: &str,
    seek_to: u64,
    max_recs: i64,
    compressed: bool,
) -> PyResult<(Bound<'py, PyList>, bool, usize)> {
    let reader = SeqReader::new(open_at(path, seek_to, compressed)?);
    let records = PyList::empty_bound(py);
    let mut n_emitted = 0usize;
    let mut stopped_early = false;
    let empty: HashSet<String> = HashSet::new();
    for rec_result in reader {
        if max_recs >= 0 && (n_emitted as i64) >= max_recs {
            break;
        }
        let seq = match rec_result {
            Ok(s) => s,
            Err(_) => { stopped_early = true; break; }
        };
        if seq.molecule_type.is_none() {
            stopped_early = true;
            break;
        }
        match record_fields(py, &seq, &empty)? {
            Some(t) => records.append(t)?,
            None => { stopped_early = true; break; }
        }
        n_emitted += 1;
    }
    Ok((records, stopped_early, n_emitted))
}

/// A parsed GenBank record held in Rust, exposing only what domainate's bulk scan
/// needs. The full record is built lazily via `materialize` (hit-only).
#[pyclass]
pub struct LeanSearchContig {
    seq: Seq,
    default_molecule_type: Option<String>,
    seq_upper: String,
}

#[pymethods]
impl LeanSearchContig {
    /// Post-swap display id (the GenBank LOCUS name): id == record.id used in the
    /// rest of domainator's in-memory convention.
    #[getter]
    fn id(&self) -> String {
        self.seq.name.clone()
            .or_else(|| self.seq.accession.clone())
            .or_else(|| self.seq.version.clone())
            .unwrap_or_else(|| "<unknown id>".to_string())
    }

    #[getter]
    fn name(&self) -> String {
        self.seq.version.clone()
            .or_else(|| self.seq.accession.clone())
            .or_else(|| self.seq.name.clone())
            .unwrap_or_default()
    }

    #[getter]
    fn molecule_type(&self) -> Option<String> {
        self.seq.molecule_type.clone().or_else(|| self.default_molecule_type.clone())
    }

    #[getter]
    fn topology(&self) -> &'static str {
        if matches!(self.seq.topology, Topology::Circular) { "circular" } else { "linear" }
    }

    #[getter]
    fn seq(&self) -> &str {
        &self.seq_upper
    }

    /// Minimal annotations dict so the object duck-types where parse_seqfiles and
    /// get_contig_search_sequence read record.annotations.
    #[getter]
    fn annotations<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new_bound(py);
        if let Some(mt) = self.molecule_type() {
            d.set_item("molecule_type", mt)?;
        }
        d.set_item("topology", self.topology())?;
        Ok(d)
    }

    fn __len__(&self) -> usize {
        self.seq_upper.len()
    }

    /// CDS peptides for the search as (surviving_feature_index, peptide). Indices
    /// are over the post-drop feature list so they match `materialize(dropped)`.
    /// Peptides are alnum-filtered to match get_prot_list.
    fn cds_peptides(&self, dropped: HashSet<String>) -> Vec<(usize, String)> {
        let mut out = Vec::new();
        let mut idx = 0usize;
        for feature in &self.seq.features {
            if !dropped.is_empty() && dropped.contains(feature.kind.as_ref()) {
                continue;
            }
            if is_searchable_cds(feature) {
                let raw = match translation_value(feature) {
                    Some(tr) => tr.to_string(),
                    None => match extract_location_seq(&feature.location, &self.seq.seq) {
                        Some(nt) => translate_dna(&nt),
                        None => String::new(),
                    },
                };
                // alnum filter matches get_prot_list (strips '*' stop chars etc.)
                let pep: String = raw.chars().filter(|c| c.is_ascii_alphanumeric()).collect();
                out.push((idx, pep));
            }
            idx += 1;
        }
        out
    }

    /// Build the full record header tuple + LeanFeature list (hit path only),
    /// dropping cleared feature kinds. Returns the same 14-tuple as parse_lean.
    fn materialize<'py>(&self, py: Python<'py>, dropped: HashSet<String>) -> PyResult<Bound<'py, PyTuple>> {
        match record_fields(py, &self.seq, &dropped)? {
            Some(t) => Ok(t),
            // Shouldn't happen (we validated at parse time), but be safe.
            None => Err(PyValueError::new_err("unmodelled feature location")),
        }
    }
}

/// Parse a partition into LeanSearchContig objects. A record is only accepted if
/// gb-io parsed its LOCUS (molecule_type Some), all its locations are modelled,
/// and every searchable CDS has a /translation (otherwise we stop early and the
/// Python side falls back to Biopython, which can translate).
#[pyfunction]
#[pyo3(signature = (path, seek_to, max_recs, default_molecule_type=None, compressed=false))]
pub fn parse_lean_search<'py>(
    py: Python<'py>,
    path: &str,
    seek_to: u64,
    max_recs: i64,
    default_molecule_type: Option<String>,
    compressed: bool,
) -> PyResult<(Bound<'py, PyList>, bool, usize)> {
    let reader = SeqReader::new(open_at(path, seek_to, compressed)?);
    let records = PyList::empty_bound(py);
    let mut n_emitted = 0usize;
    let mut stopped_early = false;

    for rec_result in reader {
        if max_recs >= 0 && (n_emitted as i64) >= max_recs {
            break;
        }
        let seq = match rec_result {
            Ok(s) => s,
            Err(_) => { stopped_early = true; break; }
        };
        if seq.molecule_type.is_none() {
            stopped_early = true;
            break;
        }
        // Validate: all locations are modelled (CDS translation, when absent, is
        // computed in Rust by cds_peptides, so it does not force a fallback).
        let mut ok = true;
        for feature in &seq.features {
            if lean_location(&feature.location, 1).is_none() {
                ok = false;
                break;
            }
        }
        if !ok {
            stopped_early = true;
            break;
        }
        let seq_upper = String::from_utf8_lossy(&seq.seq.to_ascii_uppercase()).into_owned();
        let contig = LeanSearchContig {
            seq,
            default_molecule_type: default_molecule_type.clone(),
            seq_upper,
        };
        records.append(Py::new(py, contig)?)?;
        n_emitted += 1;
    }
    Ok((records, stopped_early, n_emitted))
}

// ---------------------------------------------------------------------------
// FASTA fast path
// ---------------------------------------------------------------------------

/// A parsed FASTA record held in Rust. A FASTA record behaves like a contig with
/// no features: for protein input the whole sequence is the search peptide, for
/// nucleotide input it is searched whole by nhmmer. Materializes to a minimal
/// record (annotations = {molecule_type}) matching Biopython's FASTA SeqRecord.
#[pyclass]
pub struct LeanFastaContig {
    record_id: String,
    description: String,
    seq: String,
    molecule_type: Option<String>,
}

#[pymethods]
impl LeanFastaContig {
    #[getter]
    fn id(&self) -> &str {
        &self.record_id
    }

    #[getter]
    fn name(&self) -> &str {
        // Biopython sets FASTA name == id (the first whitespace token).
        &self.record_id
    }

    #[getter]
    fn description(&self) -> &str {
        &self.description
    }

    #[getter]
    fn seq(&self) -> &str {
        &self.seq
    }

    #[getter]
    fn molecule_type(&self) -> Option<String> {
        self.molecule_type.clone()
    }

    /// Minimal annotations dict (only molecule_type), matching a Biopython FASTA
    /// SeqRecord after parse_seqfiles sets molecule_type.
    #[getter]
    fn annotations<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new_bound(py);
        if let Some(mt) = &self.molecule_type {
            d.set_item("molecule_type", mt)?;
        }
        Ok(d)
    }

    fn __len__(&self) -> usize {
        self.seq.len()
    }

    /// FASTA records have no CDS features; the protein sequence (if any) is taken
    /// from `.seq` by get_prot_list, and nucleotide records are searched whole.
    fn cds_peptides(&self, _dropped: HashSet<String>) -> Vec<(usize, String)> {
        Vec::new()
    }

    /// Materialize to the fields a minimal LeanContig is built from on the Python
    /// side: (id, description, seq, molecule_type). No features, no header
    /// normalization (matches Biopython's FASTA SeqRecord).
    fn materialize<'py>(&self, py: Python<'py>, _dropped: HashSet<String>) -> Bound<'py, PyTuple> {
        PyTuple::new_bound(py, [
            self.record_id.clone().into_py(py),
            self.description.clone().into_py(py),
            self.seq.clone().into_py(py),
            self.molecule_type.clone().into_py(py),
        ])
    }
}

/// Parse a FASTA partition into LeanFastaContig objects. `seek_to`/`compressed`
/// follow the same convention as the GenBank parser. Sequence case is preserved
/// (Biopython does not upper-case FASTA), and the id is the first whitespace
/// token while the description is the full title line.
#[pyfunction]
#[pyo3(signature = (path, seek_to, max_recs, default_molecule_type=None, compressed=false))]
pub fn parse_fasta_search<'py>(
    py: Python<'py>,
    path: &str,
    seek_to: u64,
    max_recs: i64,
    default_molecule_type: Option<String>,
    compressed: bool,
) -> PyResult<(Bound<'py, PyList>, bool, usize)> {
    let mut reader = BufReader::new(open_at(path, seek_to, compressed)?);
    let records = PyList::empty_bound(py);
    let mut n_emitted = 0usize;

    let mut title: Option<String> = None;
    let mut seq = String::new();
    let mut line = Vec::new();

    // Helper to push the in-progress record.
    macro_rules! flush {
        () => {{
            if let Some(t) = title.take() {
                if max_recs >= 0 && (n_emitted as i64) >= max_recs {
                    return Ok((records, false, n_emitted));
                }
                let record_id = t.split_whitespace().next().unwrap_or("").to_string();
                let contig = LeanFastaContig {
                    record_id,
                    description: t,
                    seq: std::mem::take(&mut seq),
                    molecule_type: default_molecule_type.clone(),
                };
                records.append(Py::new(py, contig)?)?;
                n_emitted += 1;
            } else {
                seq.clear();
            }
        }};
    }

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        if line.first() == Some(&b'>') {
            flush!();
            // title = everything after '>', trailing newline/CR stripped.
            let raw = &line[1..];
            let end = raw.iter().rposition(|&b| b != b'\n' && b != b'\r').map_or(0, |i| i + 1);
            title = Some(String::from_utf8_lossy(&raw[..end]).into_owned());
        } else if title.is_some() {
            for &b in &line[..] {
                if b != b'\n' && b != b'\r' {
                    seq.push(b as char);
                }
            }
        }
    }
    flush!();

    Ok((records, false, n_emitted))
}
