"""Lean record model for the domain_search hot path.

The domain_search bottleneck is constructing Biopython ``SeqRecord`` /
``SeqFeature`` / ``FeatureLocation`` / position objects for *every* record of a
large database (measured: ~66% of the per-record construction cost is the
feature/location/position object graph, only ~34% is qualifiers). See
fast_parser_plan.md.

``LeanContig`` / ``LeanFeature`` are slotted, allocation-light stand-ins built
during parsing. They carry exactly what ``domainate`` needs to scan a record and
run pyhmmer, plus enough to faithfully rebuild a full Biopython ``SeqRecord`` for
the rare records that produce a hit (the "hit-conversion boundary"):

    bulk path (all records)  : parse -> LeanContig, scan CDSs, search
    hit path (rare)          : LeanContig.to_seqrecord() -> existing machinery

A lean location is stored inline on ``LeanFeature`` as primitive data (no
``FeatureLocation``/position objects):

    parts:    tuple of (start:int, end:int, strand:int, before:bool, after:bool),
              already in final Biopython part order (complement(join) parts
              pre-reversed); strand is per-part so mixed-strand (trans-spliced)
              joins are represented faithfully
    operator: "join" / "order" / None
    between:  True for a GenBank ``x^y`` between-site

Records are produced in the *pre-(name/id)-swap* convention (``id`` = VERSION,
``name`` = LOCUS) exactly like the other parse_seqfiles backends, so the existing
``swap_name_id`` step applies uniformly.
"""

from domainator.Bio.Seq import Seq
from domainator.Bio.SeqRecord import SeqRecord
from domainator.Bio.SeqFeature import (
    SeqFeature,
    FeatureLocation,
    CompoundLocation,
    ExactPosition,
    BeforePosition,
    AfterPosition,
)
try:
    from domainator import _gbfast  # optional native (Rust) acceleration
except ImportError:
    _gbfast = None

# The Rust-held search record types (None when the native extension is absent).
LeanSearchContig = getattr(_gbfast, "LeanSearchContig", None)  # GenBank
LeanFastaContig = getattr(_gbfast, "LeanFastaContig", None)    # FASTA


class LeanParseError(Exception):
    """Raised when the native (Rust) parser stops early on a record it cannot
    handle (an unparseable LOCUS line, an unmodelled location, a rejected
    sequence), carrying how many records it already emitted so callers can fall
    back to the Biopython parser for the remainder."""


# --- GenBank header-field normalization (mirrors Bio.GenBank's parsing rules) ---
# These operate on the raw header strings the native parser hands back, so the
# per-record normalization is single-sourced here.

def _strip_trailing_period(value):
    """Remove a single trailing '.' (GenBank treats it as syntax, not data)."""
    if value and value.endswith("."):
        return value[:-1]
    return value


def _split_keywords(keyword_string):
    """Mirror Bio.GenBank's _split_keywords for the KEYWORDS line."""
    if keyword_string is None:
        return [""]
    if keyword_string == "" or keyword_string == ".":
        keywords = ""
    elif keyword_string[-1] == ".":
        keywords = keyword_string[:-1]
    else:
        keywords = keyword_string
    return [x.strip() for x in keywords.split(";")]


def _split_source(source_name):
    """Mirror Bio.GenBank consumer.source: strip one trailing '.', '' stays ''."""
    if not source_name:
        return ""
    if source_name[-1] == ".":
        return source_name[:-1]
    return source_name


def _split_organism_block(organism_block):
    """Split an ORGANISM block into (organism, taxonomy_list).

    The native parser concatenates the ORGANISM line and the indented lineage
    continuation lines with newlines: the first line is the organism name, the
    remainder is the taxonomy lineage (';'-separated, terminated by '.').
    """
    if not organism_block:
        return "", []
    lines = organism_block.split("\n")
    organism = lines[0].strip()
    taxonomy_string = " ".join(line.strip() for line in lines[1:]).strip()
    if not taxonomy_string or taxonomy_string == ".":
        return organism, []
    if taxonomy_string[-1] == ".":
        taxonomy_string = taxonomy_string[:-1]
    taxonomy = [x.strip() for x in taxonomy_string.split(";") if x.strip()]
    return organism, taxonomy


def _sequence_version(version_string):
    """Extract the integer sequence version from a 'ACCESSION.N' VERSION token."""
    if not version_string or "." not in version_string:
        return None
    tail = version_string.rsplit(".", 1)[1]
    try:
        return int(tail)
    except ValueError:
        return None


class LeanFeature:
    """An allocation-light feature: type, inline location primitives, qualifiers."""

    __slots__ = ("type", "parts", "operator", "between", "qualifiers")

    def __init__(self, type, parts, operator, between, qualifiers):
        self.type = type
        self.parts = parts          # tuple of (start, end, strand, before, after)
        self.operator = operator    # "join" | "order" | None
        self.between = between       # bool
        self.qualifiers = qualifiers # dict[str, list[str]]

    def __len__(self):
        # Location length (sum of parts); used by get_taxid's source-sorting.
        return sum(end - start for start, end, _s, _b, _a in self.parts)

    @property
    def strand(self):
        # Feature-level strand: the first part's strand (parts are uniform-strand
        # except for rare trans-spliced joins).
        return self.parts[0][2] if self.parts else 0

    @property
    def start(self):
        return min(p[0] for p in self.parts)

    @property
    def end(self):
        return max(p[1] for p in self.parts)


class LeanContig:
    """An allocation-light contig/record: the bulk-path stand-in for SeqRecord.

    Exposes ``.id`` / ``.name`` / ``.description`` / ``.seq`` / ``.annotations`` /
    ``.features`` / ``.dbxrefs`` so it duck-types where domainate's scan and
    get_taxid expect a SeqRecord. ``swap_name_id`` works on it (it just swaps the
    two string attributes). ``.seq`` is a plain ``str`` (uppercased), not a ``Seq``.
    """

    __slots__ = ("id", "name", "description", "seq", "annotations", "features", "dbxrefs", "_length")

    def __init__(self, id, name, description, seq, annotations, features, dbxrefs=None):
        self.id = id
        self.name = name
        self.description = description
        self.seq = seq                # str, uppercased
        self.annotations = annotations
        self.features = features
        self.dbxrefs = dbxrefs if dbxrefs is not None else []
        self._length = len(seq)

    def __len__(self):
        return self._length


def build_lean_contig(name, accession, version, definition, molecule_type, circular,
                      division, date_str, keywords, source_name, source_organism,
                      dblink, seq, features, default_molecule_type=None):
    """Assemble a LeanContig from already-extracted header fields + a features list.

    Shared by the native (Rust) producer paths so the per-record header
    normalization is single-sourced (no Python/Rust divergence).
    ``features`` is a prebuilt list of LeanFeature; ``seq`` is the uppercased str;
    ``date_str`` is the already-formatted "%d-%b-%Y" string or None.
    """
    rec_name = name or accession or version or ""
    record_id = version or accession or name or "<unknown id>"
    if definition is not None and "\n" in definition:
        definition = " ".join(part.strip() for part in definition.split("\n"))
    description = _strip_trailing_period(definition) or ""

    annotations = {}
    mol = molecule_type or default_molecule_type
    if mol is not None:
        annotations["molecule_type"] = mol
    annotations["topology"] = "circular" if circular else "linear"
    if division:
        annotations["data_file_division"] = division
    if date_str is not None:
        annotations["date"] = date_str
    if accession:
        annotations["accessions"] = [accession]
    sequence_version = _sequence_version(version)
    if sequence_version is not None:
        annotations["sequence_version"] = sequence_version
    annotations["keywords"] = _split_keywords(keywords)
    if source_name is not None or source_organism is not None:
        annotations["source"] = _split_source(source_name)
        organism, taxonomy = _split_organism_block(source_organism)
        annotations["organism"] = organism
        annotations["taxonomy"] = taxonomy

    dbxrefs = []
    if dblink:
        for entry in dblink.split("\n"):
            entry = entry.strip()
            if entry:
                dbxrefs.append(entry.replace(": ", ":"))

    return LeanContig(record_id, rec_name, description, seq, annotations, features, dbxrefs)


# --- Biopython SeqRecord <-> lean (fallback path + hit-conversion boundary) ---

def seqrecord_to_lean(record):
    """Convert an already-parsed Biopython SeqRecord to a LeanContig.

    Used on the fallback path (records the native parser cannot parse are parsed by Biopython
    then converted). Assumes ``record`` is in the same convention parse_seqfiles
    yields it.
    """
    features = []
    for f in record.features:
        loc = f.location
        parts = tuple(
            (
                int(p.start),
                int(p.end),
                p.strand if p.strand is not None else 0,
                type(p.start).__name__ == "BeforePosition",
                type(p.end).__name__ == "AfterPosition",
            )
            for p in loc.parts
        )
        operator = getattr(loc, "operator", None) if len(loc.parts) > 1 else None
        qualifiers = {k: list(v) for k, v in f.qualifiers.items()}
        features.append(LeanFeature(f.type, parts, operator, False, qualifiers))
    return LeanContig(
        record.id,
        record.name,
        record.description,
        str(record.seq),
        dict(record.annotations),
        features,
        list(record.dbxrefs),
    )


def _position(value, before=False, after=False):
    if before:
        return BeforePosition(value)
    if after:
        return AfterPosition(value)
    return ExactPosition(value)


def _lean_feature_to_seqfeature(feature):
    if feature.between:
        start, end, strand, _b, _a = feature.parts[0]
        pos = ExactPosition(end)
        location = FeatureLocation(pos, pos, strand=strand)
    elif len(feature.parts) == 1:
        start, end, strand, before, after = feature.parts[0]
        location = FeatureLocation(_position(start, before=before), _position(end, after=after), strand=strand)
    else:
        sub = [
            FeatureLocation(_position(s, before=b), _position(e, after=a), strand=st)
            for (s, e, st, b, a) in feature.parts
        ]
        location = CompoundLocation(sub, operator=feature.operator or "join")
    return SeqFeature(location=location, type=feature.type, qualifiers={k: list(v) for k, v in feature.qualifiers.items()})


def lean_translate(feature, contig_seq):
    """Translate a lean CDS feature against the contig sequence.

    Reuses Biopython's ``SeqFeature.translate`` (via a one-off reconstructed
    feature) so the result is byte-identical to what ``clean_rec`` produces on the
    Biopython path. The Rust core (increment 3) will replace this with an in-Rust
    translation; this keeps the pure-Python increment correct.
    """
    seqfeature = _lean_feature_to_seqfeature(feature)
    return seqfeature.translate(Seq(contig_seq), cds=False)


def lean_to_seqrecord(lean):
    """Rebuild a full Biopython SeqRecord from a LeanContig (hit-conversion boundary)."""
    record = SeqRecord(
        Seq(lean.seq),
        id=lean.id,
        name=lean.name,
        description=lean.description,
        dbxrefs=list(lean.dbxrefs),
        annotations=dict(lean.annotations),
        features=[_lean_feature_to_seqfeature(f) for f in lean.features],
    )
    return record


def iter_lean_genbank_native(path, seek_to=0, max_recs=-1, default_molecule_type=None, compressed=False):
    """Yield LeanContigs for a GenBank partition using the native (Rust) parser.

    The Rust side parses up to ``max_recs`` records into prebuilt LeanFeature lists
    + header fields; this wraps each into a LeanContig (header normalization stays
    in Python). When the Rust parser stops early (an unparseable/unreliable record),
    this raises LeanParseError carrying the count already yielded, so the caller
    falls back to Biopython for the remaining records. ``compressed`` indicates a
    BGZF input, in which case ``seek_to`` is a 64-bit virtual offset.
    """
    records, stopped_early, n_emitted = _gbfast.parse_lean(str(path), int(seek_to), int(max_recs), bool(compressed))
    for rec in records:
        yield _build_lean_from_fields(rec, default_molecule_type)
    if stopped_early:
        raise LeanParseError(
            f"native parser stopped after {n_emitted} record(s); falling back to Biopython"
        )


def _build_lean_from_fields(fields, default_molecule_type):
    (name, accession, version, definition, molecule_type, circular, division,
     date_str, keywords, source_name, source_organism, dblink, seq, features) = fields
    return build_lean_contig(
        name, accession, version, definition, molecule_type, circular, division,
        date_str, keywords, source_name, source_organism, dblink, seq, features,
        default_molecule_type=default_molecule_type,
    )


def iter_lean_search_native(path, seek_to=0, max_recs=-1, default_molecule_type=None, compressed=False):
    """Yield LeanSearchContig objects for a GenBank partition (the bulk fast path).

    Each yielded object holds the parsed record in Rust and is materialized to a
    full record only on a hit. Raises LeanParseError (with the count already
    yielded) when the native parser stops early, so the caller falls back to
    Biopython for the remaining records. ``compressed`` indicates a BGZF input,
    in which case ``seek_to`` is a 64-bit virtual offset.
    """
    records, stopped_early, n_emitted = _gbfast.parse_lean_search(
        str(path), int(seek_to), int(max_recs), default_molecule_type, bool(compressed)
    )
    yield from records
    if stopped_early:
        raise LeanParseError(
            f"native parser stopped after {n_emitted} record(s); falling back to Biopython"
        )


def iter_lean_fasta_native(path, seek_to=0, max_recs=-1, default_molecule_type=None, compressed=False):
    """Yield LeanFastaContig objects for a FASTA partition (the bulk fast path).

    Like the GenBank search path but for FASTA: records hold (id, description,
    sequence) and are materialized to minimal records only on a hit. ``compressed``
    indicates a BGZF input (``seek_to`` is then a virtual offset).
    """
    records, _stopped, _n = _gbfast.parse_fasta_search(
        str(path), int(seek_to), int(max_recs), default_molecule_type, bool(compressed)
    )
    yield from records


def search_contig_to_lean(search_contig, dropped_types, default_molecule_type=None):
    """Materialize a LeanSearchContig (a GenBank hit) into a full post-swap LeanContig."""
    fields = search_contig.materialize(set(dropped_types) if dropped_types else set())
    lean = _build_lean_from_fields(fields, default_molecule_type)
    # The search object's id/name are already post-swap; build_lean_contig produced
    # pre-swap, so swap to match the in-memory convention used everywhere else.
    lean.id, lean.name = lean.name, lean.id
    return lean


def search_fasta_to_lean(fasta_contig, dropped_types, default_molecule_type=None):
    """Materialize a LeanFastaContig (a hit) into a minimal LeanContig.

    Matches a Biopython-parsed FASTA SeqRecord: id == name == first title token,
    description == full title, no features, annotations == {molecule_type}.
    """
    record_id, description, seq, molecule_type = fasta_contig.materialize(set())
    mol = molecule_type or default_molecule_type
    annotations = {"molecule_type": mol} if mol is not None else {}
    return LeanContig(record_id, record_id, description, seq, annotations, [])


# Unified handling of the native "search" record types (GenBank + FASTA). domainate
# and parse_seqfiles treat these the same: bulk-scan cheaply, materialize on hit.
LEAN_SEARCH_TYPES = tuple(t for t in (LeanSearchContig, LeanFastaContig) if t is not None)


def materialize_lean_search(search_contig, dropped_types, default_molecule_type=None):
    """Materialize any native search record (GenBank or FASTA) into a LeanContig."""
    if LeanFastaContig is not None and isinstance(search_contig, LeanFastaContig):
        return search_fasta_to_lean(search_contig, dropped_types, default_molecule_type)
    return search_contig_to_lean(search_contig, dropped_types, default_molecule_type)
