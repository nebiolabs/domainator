"""Shared plumbing for structure-based search and annotation.

This module has no CLI. It is imported by structure_to_genbank.py,
structure_domainate.py and structure_search.py.

Three things live here:

  * Input resolution. A structure input may be a prebuilt aligner database, a
    directory of structures, a glob, or individual files. resolve_structure_input()
    normalizes all of those.

  * The backend seam. StructureAligner is an abstract structural aligner; FoldseekAligner
    is the only implementation today. Backends declare which alignment types and which
    optional metrics they support, so an unsupported request produces a clear error
    instead of a backend crash. Adding reseek or a LoLalign-capable foldseek means adding
    a subclass and an ALIGNERS entry.

  * Translation into domainator's annotation vocabulary. Hits become
    domainate.SearchResult objects keyed by input record name, and input entries become
    protein SeqRecords, so the existing annotation writers can be reused unchanged.

Alignment direction is fixed: references are the QUERY and the input database is the
TARGET. That makes E-values comparable to the rest of domainator -- they are always
computed against the database being searched, not against however many references the
user happened to supply. It also means the coordinates of a hit *on the input* come from
the target columns, which is the opposite of domainator's older sequence-driven foldseek
path (see structure_hits_to_search_results).
"""

import contextlib
import glob as glob_module
import os
import shutil
import subprocess
import sys
import tempfile
import warnings
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, FrozenSet, Iterable, Iterator, List, NamedTuple, Optional, Sequence, Tuple

from domainator.Bio.Seq import Seq
from domainator.Bio.SeqRecord import SeqRecord
from domainator import utils
from domainator.domainate import SearchResult
from domainator.output_guardrails import (
    make_temporary_output_path, OutputSizeLimitExceeded, write_records_with_limit,
)
from domainator.utils import open_if_is_name_for_write

# Structure file extensions foldseek and reseek can read. Deliberately NOT added to
# utils.EXTENSION_TO_TYPE: that map feeds parse_seqfiles, which cannot read structures.
STRUCTURE_EXTENSIONS = {"pdb", "cif", "mmcif", "ent", "bcif", "pdb1", "ent"}

# Optional per-hit metrics a backend may be asked to compute. These map onto the
# same-named optional fields of domainate.SearchResult.
METRIC_CHOICES = ("tmscore", "lddt", "rmsd", "prob")

# Metrics that require C-alpha coordinates in the databases being compared.
COORDINATE_METRICS = frozenset({"tmscore", "lddt", "rmsd"})

# foldseek --alignment-type values that require C-alpha coordinates.
COORDINATE_ALIGNMENT_TYPES = frozenset({1})


class StructureInputError(ValueError):
    """Raised when a -i/-r value cannot be interpreted as structures or a database."""


def is_structure_path(path) -> bool:
    """Whether path looks like a structure file, ignoring one compression suffix."""
    suffixes = [s[1:].lower() for s in Path(path).suffixes]
    if suffixes and suffixes[-1] in utils.COMPRESSION_EXTENSIONS:
        suffixes = suffixes[:-1]
    return bool(suffixes) and suffixes[-1] in STRUCTURE_EXTENSIONS


def is_aligner_db(value) -> bool:
    """Whether value is the prefix of an existing mmseqs/foldseek-style database.

    Both files are required: a lone .dbtype would also match some unrelated paths.
    """
    return Path(f"{value}.dbtype").is_file() and Path(f"{value}.index").is_file()


def resolve_structure_input(value) -> Tuple[str, List[str]]:
    """Interpret one -i/-r value.

    Returns:
        ("db", [prefix]) if value is a prebuilt aligner database prefix, or
        ("structures", [paths...]) for a directory, glob, or individual file.

    Directories are searched recursively. Globs are expanded here rather than relying on
    the shell, because values coming from a --config file are never shell-expanded.
    """
    value = str(value)

    if is_aligner_db(value):
        return "db", [value]

    path = Path(value)
    if path.is_dir():
        found = sorted(str(p) for p in path.rglob("*") if p.is_file() and is_structure_path(p))
        if not found:
            raise StructureInputError(
                f"No structure files found under directory '{value}'. "
                f"Recognized extensions: {', '.join(sorted(STRUCTURE_EXTENSIONS))} "
                "(optionally .gz)."
            )
        return "structures", found

    if path.is_file():
        if not is_structure_path(path):
            raise StructureInputError(
                f"'{value}' is not a recognized structure file. Recognized extensions: "
                f"{', '.join(sorted(STRUCTURE_EXTENSIONS))} (optionally .gz)."
            )
        return "structures", [value]

    if any(char in value for char in "*?["):
        found = sorted(p for p in glob_module.glob(value, recursive=True)
                       if os.path.isfile(p) and is_structure_path(p))
        if not found:
            raise StructureInputError(f"Glob '{value}' matched no structure files.")
        return "structures", found

    raise StructureInputError(
        f"Could not interpret '{value}'. It is not an aligner database prefix "
        f"('{value}.dbtype' and '{value}.index' do not both exist), not a directory, "
        "not an existing structure file, and not a glob that matched anything."
    )


def resolve_structure_inputs(values: Sequence[str]) -> Tuple[str, List[str]]:
    """Resolve several -i/-r values into one database prefix or one list of structures.

    Mixing a prebuilt database with loose structure files in a single argument is
    rejected: the two need different handling and silently ignoring one would be worse
    than an error. Multiple databases are likewise rejected, since a search runs against
    exactly one target database.
    """
    kinds = []
    structures: List[str] = []
    dbs: List[str] = []
    for value in values:
        kind, resolved = resolve_structure_input(value)
        kinds.append(kind)
        if kind == "db":
            dbs.append(resolved[0])
        else:
            structures.extend(resolved)

    if dbs and structures:
        raise StructureInputError(
            "Cannot mix a prebuilt aligner database with structure files in the same "
            f"argument. Got database(s) {dbs} and {len(structures)} structure file(s)."
        )
    if len(dbs) > 1:
        raise StructureInputError(
            f"Only one aligner database can be given per argument, got {dbs}."
        )
    if dbs:
        return "db", dbs
    # Preserve order but drop duplicates, which overlapping globs make easy to produce.
    deduplicated = list(dict.fromkeys(structures))
    return "structures", deduplicated


def database_label(prefix) -> str:
    """Label for the /database qualifier.

    Path.name rather than Path.stem: a database prefix has no extension to strip, so
    .stem would collapse 'refs.v1' and 'refs.v2' onto the same label.
    """
    name = Path(prefix).name
    # Strip a compression suffix then one structure suffix, so 'ref.pdb.gz' -> 'ref'.
    for _ in range(2):
        if is_structure_path(name) or Path(name).suffix[1:].lower() in utils.COMPRESSION_EXTENSIONS:
            name = Path(name).stem
    return name



def input_label(values: Sequence[str], kind: str, resolved: Sequence[str]) -> str:
    """A human-meaningful name for a resolved -i/-r argument.

    Used for the /database qualifier, which reaches enum_report, summary_report and
    extract_domains --databases, so it has to describe what the user supplied rather than
    the temporary prefix a database happened to be built at.
    """
    if kind == "db":
        return database_label(resolved[0])
    if len(values) == 1 and Path(values[0]).is_dir():
        return Path(values[0]).name or "references"
    if len(resolved) == 1:
        return database_label(resolved[0])
    parents = {str(Path(path).parent) for path in resolved}
    if len(parents) == 1:
        return Path(next(iter(parents))).name or "references"
    return "references"


def strip_structure_extension(name: str) -> str:
    """'1abc.pdb_A' -> '1abc_A'; leaves names without a structure extension alone.

    foldseek names database entries '<file basename>_<chain>', so the file extension ends
    up in the middle of the entry name. Domain names reach plots and reports, so trimming
    it is usually what you want.
    """
    for extension in sorted(STRUCTURE_EXTENSIONS, key=len, reverse=True):
        needle = f".{extension}_"
        index = name.lower().rfind(needle)
        if index != -1:
            return name[:index] + "_" + name[index + len(needle):]
        suffix = f".{extension}"
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return name


class StructureHit(NamedTuple):
    """One alignment, in backend-neutral form.

    query is a reference, target is an input record: see the module docstring on
    alignment direction.
    """
    query: str
    target: str
    qstart: int
    qend: int
    qlen: int
    tstart: int
    tend: int
    tlen: int
    evalue: float
    bits: float
    fident: float
    alnlen: int
    tseq: str = ""
    qheader: str = ""
    theader: str = ""
    tmscore: Optional[float] = None
    lddt: Optional[float] = None
    rmsd: Optional[float] = None
    prob: Optional[float] = None


class StructureDB(NamedTuple):
    """A structural database on disk.

    Attributes:
        prefix: path prefix shared by the database's files.
        owned: True when this process built it in a temporary location, so callers know
            whether it is safe to delete or must be left alone.
        source_paths: the ordered input paths it was built from, empty for a prebuilt
            database. Used to resolve provenance out of the .lookup file.
    """
    prefix: str
    owned: bool = False
    source_paths: Tuple[str, ...] = ()

    def has_coordinates(self) -> bool:
        """Whether a C-alpha sub-database is present.

        Databases built by foldseek createdb from real structures have one. Databases
        built from amino acid + 3Di FASTA (domainator's older sequence-driven path) do
        not, so they cannot do TM-alignment or any coordinate-based metric.
        """
        return Path(f"{self.prefix}_ca").is_file()


def resolve_binary(name: str, override: Optional[str] = None, install_hint: str = "") -> str:
    """Resolve an external executable, raising a descriptive error when it is missing."""
    if override is not None:
        if not (os.path.isfile(override) and os.access(override, os.X_OK)) and shutil.which(override) is None:
            raise RuntimeError(f"'{override}' is not an executable file.")
        return override
    resolved = shutil.which(name)
    if resolved is None:
        hint = f" {install_hint}" if install_hint else ""
        raise RuntimeError(
            f"Could not find '{name}' on PATH.{hint} "
            f"Alternatively pass an explicit path with --{name}_path."
        )
    return resolved


def run_command(argv: Sequence[str], description: str, expect: Optional[str] = None, env=None) -> None:
    """Run an external command, raising RuntimeError with stderr on failure.

    A nonzero exit code alone is not treated as failure when `expect` names an artifact
    that was produced: foldseek's multi-step workflows exit nonzero during their cleanup
    phase even when the run succeeded.

    subprocess.run rather than Popen followed by stderr.read(), which deadlocks once the
    child writes more than the pipe buffer.
    """
    argv = [str(a) for a in argv]
    completed = subprocess.run(argv, capture_output=True, env=env)
    if completed.returncode == 0:
        return
    if expect is not None and Path(expect).exists():
        return
    stderr = completed.stderr.decode("utf-8", errors="replace")
    tail = "\n".join(stderr.strip().split("\n")[-25:])
    raise RuntimeError(
        f"{description} failed with exit code {completed.returncode}.\n"
        f"Command: {' '.join(argv)}\n{tail}"
    )


class StructureAligner(ABC):
    """A structural aligner backend.

    Subclasses declare their capabilities so that callers can reject an unsupported
    request with a clear message rather than letting the backend fail obscurely.
    """

    name: str = "abstract"
    supports_alignment_types: FrozenSet[int] = frozenset()
    supports_metrics: FrozenSet[str] = frozenset()
    path_option: str = "--aligner_path"

    def __init__(self, bin_path: Optional[str] = None, cpu: int = 0, device: Optional[str] = None,
                 extra_args: Optional[Sequence[str]] = None, tmp_dir: Optional[str] = None):
        self.cpu = cpu
        self.device = device
        self.extra_args = list(extra_args) if extra_args else []
        self.tmp_dir = tmp_dir
        self.bin_path = self._resolve_bin(bin_path)

    @abstractmethod
    def _resolve_bin(self, bin_path: Optional[str]) -> str:
        """Resolve the backend executable. Called from __init__, before any real work."""

    @abstractmethod
    def build_db(self, inputs: Sequence[str], out_prefix: str) -> StructureDB:
        """Build a database from structure files."""

    @abstractmethod
    def search(self, query_db: StructureDB, target_db: StructureDB, *, evalue: float,
               max_seqs: int, want_tseq: bool, alignment_type: int,
               metrics: Sequence[str], work_dir: str,
               sort_by_target: bool = False) -> Iterator[StructureHit]:
        """Align every entry of query_db against target_db.

        When sort_by_target is set, hits are yielded grouped by target, so a caller can
        stream one input record at a time instead of holding every hit in memory.
        """

    @abstractmethod
    def iter_sequences(self, db: StructureDB) -> Iterator[Tuple[str, str, Optional[str]]]:
        """Yield (name, amino_acid_sequence, source_path_or_None) for every entry."""

    def iter_3di(self, db: StructureDB) -> Iterator[Tuple[str, str]]:
        """Yield (name, 3Di_sequence) for every entry, if the backend has 3Di states."""
        raise NotImplementedError(f"{self.name} does not expose a 3Di alphabet.")

    def check_capabilities(self, alignment_type: int, metrics: Sequence[str],
                           databases: Sequence[StructureDB]) -> None:
        """Fail early on an unsupported or impossible request."""
        if alignment_type not in self.supports_alignment_types:
            raise RuntimeError(
                f"Backend '{self.name}' does not support alignment type {alignment_type}. "
                f"Supported: {sorted(self.supports_alignment_types)}."
            )
        unsupported = sorted(set(metrics) - set(self.supports_metrics))
        if unsupported:
            raise RuntimeError(
                f"Backend '{self.name}' cannot compute {unsupported}. "
                f"Supported metrics: {sorted(self.supports_metrics)}."
            )

        needs_coordinates = sorted(set(metrics) & COORDINATE_METRICS)
        if alignment_type in COORDINATE_ALIGNMENT_TYPES:
            needs_coordinates.append(f"alignment type {alignment_type}")
        if not needs_coordinates:
            return
        for db in databases:
            if not db.has_coordinates():
                verb = "requires" if len(needs_coordinates) == 1 else "require"
                raise RuntimeError(
                    f"{', '.join(needs_coordinates)} {verb} C-alpha coordinates, but the "
                    f"database '{db.prefix}' has no '{Path(db.prefix).name}_ca' file. "
                    "Databases built from amino acid and 3Di sequences (rather than from "
                    "structure files) carry no coordinates. Rebuild it from the original "
                    "structures, or drop the coordinate-dependent options."
                )


class FoldseekAligner(StructureAligner):
    """foldseek backend."""

    name = "foldseek"
    supports_alignment_types = frozenset({0, 1, 2})
    supports_metrics = frozenset(METRIC_CHOICES)
    path_option = "--foldseek_path"

    # convertalis column name for each metric we expose.
    METRIC_COLUMNS = {
        "tmscore": "alntmscore",
        "lddt": "lddt",
        "rmsd": "rmsd",
        "prob": "prob",
    }

    BASE_COLUMNS = ("query", "target", "qheader", "fident", "alnlen",
                    "qstart", "qend", "qlen", "tstart", "tend", "tlen", "evalue", "bits")

    def _resolve_bin(self, bin_path):
        return resolve_binary(
            "foldseek", bin_path,
            install_hint="Install it with 'conda install -c conda-forge -c bioconda foldseek'.",
        )

    def _thread_args(self) -> List[str]:
        # cpu == 0 means "let foldseek decide", which is its own default.
        return ["--threads", str(self.cpu)] if self.cpu and self.cpu > 0 else []

    def _gpu_args(self) -> Tuple[List[str], Optional[dict]]:
        """foldseek selects a GPU through CUDA_VISIBLE_DEVICES, not a command line flag."""
        if self.device is None or self.device == "cpu":
            return [], None
        env = None
        if ":" in self.device:
            env = os.environ.copy()
            env["CUDA_VISIBLE_DEVICES"] = self.device.split(":", 1)[1]
        return ["--gpu", "1"], env

    def build_db(self, inputs, out_prefix) -> StructureDB:
        """foldseek createdb, driven by a TSV list of input paths.

        The TSV form is used rather than passing paths as arguments because it has no
        ARG_MAX limit, it sidesteps foldseek's "only one directory can be given"
        restriction, and it fixes the file numbering that .lookup refers to.
        """
        inputs = list(inputs)
        if not inputs:
            raise StructureInputError("No structure files to build a database from.")
        list_path = f"{out_prefix}.inputs.tsv"
        Path(list_path).parent.mkdir(parents=True, exist_ok=True)
        with open(list_path, "w") as handle:
            for path in inputs:
                handle.write(f"{os.path.abspath(path)}\n")

        argv = [self.bin_path, "createdb", list_path, out_prefix,
                "--chain-name-mode", "1", "--write-lookup", "1", "-v", "1"]
        argv += self._thread_args()
        run_command(argv, "foldseek createdb", expect=f"{out_prefix}.dbtype")
        if not Path(f"{out_prefix}.dbtype").is_file():
            raise RuntimeError(
                f"foldseek createdb produced no database at '{out_prefix}'. "
                "Check that the input files are readable structures."
            )
        return StructureDB(prefix=out_prefix, owned=True, source_paths=tuple(inputs))

    @staticmethod
    def _source_key(path) -> str:
        """Reproduce the name foldseek records in <prefix>.source for an input file.

        foldseek strips a compression suffix, then up to two more suffixes: 'x.pdb' -> 'x',
        'x.pdb.gz' -> 'x', 'some.long.name.pdb' -> 'some.long', and
        'some.long.name.pdb.gz' -> 'some.long'. (Note this is one suffix more aggressive
        than the entry naming in .lookup, which keeps 'some.long.name'.)
        """
        name = Path(path).name
        if Path(name).suffix[1:].lower() in utils.COMPRESSION_EXTENSIONS:
            name = Path(name).stem
        for _ in range(2):
            name = Path(name).stem
        return name

    def _lookup_sources(self, db: StructureDB) -> Dict[str, str]:
        """Map entry name -> originating input path.

        <prefix>.lookup gives key -> (entry name, file index) and <prefix>.source gives
        file index -> a truncated source name. The file index follows foldseek's own
        ordering of the inputs, NOT the order they were listed in, so .source has to be
        consulted rather than indexing into source_paths directly.

        .source truncates dotted names ('some.long.name.pdb' is recorded as 'some.long'),
        so the truncated name is matched back to a full input path. Provenance is
        cosmetic, so an ambiguous match is dropped rather than guessed at.
        """
        lookup_path = Path(f"{db.prefix}.lookup")
        source_path = Path(f"{db.prefix}.source")
        if not lookup_path.is_file() or not db.source_paths:
            return {}

        by_key: Dict[str, List[str]] = {}
        for path in db.source_paths:
            by_key.setdefault(self._source_key(path), []).append(path)

        index_to_path: Dict[int, str] = {}
        if source_path.is_file():
            with open(source_path) as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 2:
                        continue
                    try:
                        file_index = int(fields[0])
                    except ValueError:
                        continue
                    candidates = by_key.get(fields[1], [])
                    if len(candidates) == 1:
                        index_to_path[file_index] = candidates[0]
        elif len(db.source_paths) == 1:
            # Single input: the index can only refer to it.
            index_to_path[0] = db.source_paths[0]

        sources: Dict[str, str] = {}
        with open(lookup_path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                try:
                    file_index = int(fields[2])
                except ValueError:
                    continue
                if file_index in index_to_path:
                    sources[fields[1]] = index_to_path[file_index]
        return sources

    @staticmethod
    def _iter_flat_db(prefix: str) -> Iterator[Tuple[int, str]]:
        """Stream (key, sequence) from an mmseqs/foldseek flat sequence database.

        The data file holds '\n\0'-terminated records and the .index file holds
        'key<TAB>offset<TAB>length' rows. Read directly rather than through
        'foldseek convert2fasta' because the _ss (3Di) sub-database has no matching header
        sub-database, so convert2fasta refuses to read it.

        Entries are seeked to one at a time rather than reading the data file in one go,
        so memory does not scale with the size of the database.
        """
        index_path = Path(f"{prefix}.index")
        if not index_path.is_file() or not Path(prefix).is_file():
            raise RuntimeError(f"'{prefix}' is not a readable foldseek database.")
        with open(index_path) as index_handle, open(prefix, "rb") as data_handle:
            for line in index_handle:
                fields = line.split()
                if len(fields) < 3:
                    continue
                key, offset, length = int(fields[0]), int(fields[1]), int(fields[2])
                data_handle.seek(offset)
                record = data_handle.read(length)
                yield key, record.rstrip(b"\x00").rstrip(b"\n").decode()

    @staticmethod
    def _read_names(prefix: str) -> Dict[int, str]:
        """Read {key: entry name} from <prefix>.lookup."""
        names = {}
        lookup_path = Path(f"{prefix}.lookup")
        if not lookup_path.is_file():
            raise RuntimeError(
                f"Database '{prefix}' has no .lookup file, so entry names cannot be "
                "recovered. Rebuild it with 'foldseek createdb --write-lookup 1'."
            )
        with open(lookup_path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2:
                    names[int(fields[0])] = fields[1]
        return names

    def iter_sequences(self, db: StructureDB):
        sources = self._lookup_sources(db)
        names = self._read_names(db.prefix)
        for key, sequence in self._iter_flat_db(db.prefix):
            name = names.get(key, str(key))
            yield name, sequence, sources.get(name)

    def iter_3di(self, db: StructureDB):
        if not Path(f"{db.prefix}_ss").is_file():
            raise RuntimeError(
                f"Database '{db.prefix}' has no 3Di sub-database ('{Path(db.prefix).name}_ss')."
            )
        names = self._read_names(db.prefix)
        for key, sequence in self._iter_flat_db(f"{db.prefix}_ss"):
            yield names.get(key, str(key)), sequence

    def search(self, query_db, target_db, *, evalue, max_seqs, want_tseq,
               alignment_type, metrics, work_dir, sort_by_target=False):
        columns = list(self.BASE_COLUMNS)
        if want_tseq:
            columns.append("tseq")
        metric_columns = [self.METRIC_COLUMNS[m] for m in metrics]
        columns += metric_columns

        aln_path = os.path.join(work_dir, "aln")
        foldseek_tmp = os.path.join(work_dir, "foldseek_tmp")
        results_path = os.path.join(work_dir, "results.tsv")
        os.makedirs(foldseek_tmp, exist_ok=True)

        gpu_args, env = self._gpu_args()
        search_argv = [self.bin_path, "search", query_db.prefix, target_db.prefix,
                       aln_path, foldseek_tmp,
                       "-e", str(evalue),
                       "--max-seqs", str(max_seqs),
                       "--alignment-type", str(alignment_type)]
        if metric_columns:
            # lddt/rmsd/alntmscore are computed from the alignment backtrace, which
            # foldseek search only stores when -a is given.
            search_argv.append("-a")
        search_argv += self._thread_args() + gpu_args + self.extra_args
        run_command(search_argv, "foldseek search", expect=f"{aln_path}.dbtype", env=env)

        convert_argv = [self.bin_path, "convertalis", query_db.prefix, target_db.prefix,
                        aln_path, results_path,
                        "--format-mode", "0",
                        "--format-output", ",".join(columns)]
        convert_argv += self._thread_args()
        run_command(convert_argv, "foldseek convertalis", expect=results_path)

        if sort_by_target:
            # convertalis groups its output by query (i.e. by reference), but records are
            # keyed by target. Sorting on the target column with coreutils sort keeps the
            # grouping pass O(one target's hits) in memory instead of O(all hits), which
            # matters when the target database has hundreds of millions of entries.
            target_column = columns.index("target") + 1
            sorted_path = os.path.join(work_dir, "results.by_target.tsv")
            run_command(
                ["sort", "-t", "\t", f"-k{target_column},{target_column}",
                 "-T", work_dir, "-o", sorted_path, results_path],
                "sorting foldseek results by target", expect=sorted_path,
            )
            results_path = sorted_path

        yield from self._parse_results(results_path, columns)

    @staticmethod
    def _parse_results(results_path: str, columns: Sequence[str]) -> Iterator[StructureHit]:
        inverse_metrics = {v: k for k, v in FoldseekAligner.METRIC_COLUMNS.items()}
        with open(results_path) as handle:
            for line in handle:
                line = line.rstrip("\n")
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) != len(columns):
                    raise RuntimeError(
                        f"Unexpected foldseek output: got {len(fields)} columns, "
                        f"expected {len(columns)} ({','.join(columns)}) in line: {line[:200]}"
                    )
                row = dict(zip(columns, fields))
                values = {
                    "query": row["query"],
                    "target": row["target"],
                    "qstart": int(row["qstart"]),
                    "qend": int(row["qend"]),
                    "qlen": int(row["qlen"]),
                    "tstart": int(row["tstart"]),
                    "tend": int(row["tend"]),
                    "tlen": int(row["tlen"]),
                    "evalue": float(row["evalue"]),
                    "bits": float(row["bits"]),
                    "fident": float(row["fident"]),
                    "alnlen": int(row["alnlen"]),
                    "tseq": row.get("tseq", ""),
                    "qheader": row.get("qheader", ""),
                }
                for column, metric in inverse_metrics.items():
                    if column in row:
                        try:
                            values[metric] = float(row[column])
                        except ValueError:
                            values[metric] = None
                yield StructureHit(**values)


ALIGNERS = {"foldseek": FoldseekAligner}


def structure_hits_to_search_results(
    hits: Iterable[StructureHit],
    db_name: str,
    evalue: float,
    min_evalue: float = 0.0,
    program: str = "foldseek",
    strip_extension: bool = True,
) -> Dict[str, List[SearchResult]]:
    """Group hits by input record, as domainate.SearchResult objects.

    References are the query and inputs are the target, so the annotation coordinates
    (start/end) come from the *target* columns and the reference coordinates
    (rstart/rend/rlen) come from the *query* columns. That is inverted relative to
    domainate.foldseek_hits_to_search_results, where the input was the query.
    """
    out: Dict[str, List[SearchResult]] = {}
    for hit in hits:
        if not (hit.evalue < evalue and hit.evalue >= min_evalue):
            continue
        name = strip_structure_extension(hit.query) if strip_extension else hit.query
        description = hit.qheader.split(" ", 1)[1] if " " in hit.qheader else ""
        # min/max because a backend may report a descending range; FeatureLocation
        # rejects start > end.
        start = min(hit.tstart, hit.tend) - 1
        end = max(hit.tstart, hit.tend)
        out.setdefault(hit.target, []).append(
            SearchResult(
                name,
                description,
                "",
                hit.evalue,
                hit.bits,
                start,
                end,
                db_name,
                hit.fident * 100.0,
                min(hit.qstart, hit.qend),
                max(hit.qstart, hit.qend),
                hit.qlen,
                program,
                1,
                hit.tmscore,
                hit.lddt,
                hit.rmsd,
                hit.prob,
            )
        )
    return out


def build_protein_record(name: str, sequence: str, source_path: Optional[str] = None,
                         threedi: Optional[str] = None) -> SeqRecord:
    """Build a protein SeqRecord for one structure chain.

    id is set as well as name because utils.write_genbank calls swap_name_id before
    writing, so the GenBank LOCUS ends up coming from id.
    """
    description = f"structure={source_path}" if source_path else ""
    record = SeqRecord(seq=Seq(sequence), id=name, name=name, description=description)
    record.annotations["molecule_type"] = "protein"
    if threedi is not None:
        if len(threedi) != len(sequence):
            warnings.warn(
                f"3Di length ({len(threedi)}) does not match sequence length "
                f"({len(sequence)}) for {name}; not storing 3Di."
            )
        else:
            from domainator.Bio.SeqFeature import SeqFeature, FeatureLocation
            record.features.append(
                SeqFeature(
                    location=FeatureLocation(0, len(sequence), strand=1),
                    type="misc_feature",
                    qualifiers={"threedi": [threedi]},
                )
            )
    return record


def add_backend_arguments(parser, include_search_arguments: bool = True) -> None:
    """Add the aligner-backend arguments shared by every structure tool.

    Defined here rather than in each tool so the three CLIs cannot drift apart, and so a
    new backend only has to be added to ALIGNERS to become selectable everywhere.
    """
    parser.add_argument('--algorithm', default="foldseek", type=str, choices=sorted(ALIGNERS),
                        help="structural alignment backend to use.")
    for name in sorted(ALIGNERS):
        parser.add_argument(f'--{name}_path', default=None, type=str,
                            help=f"path to the {name} executable. If not supplied, {name} is looked up on PATH.")
    parser.add_argument('--aligner_arg', nargs='+', default=None, type=str,
                        help="extra arguments appended verbatim to the backend's search command line, for options domainator does not wrap. Example: --aligner_arg --prefilter-mode 2")
    parser.add_argument('--cpu', type=int, default=0,
                        help="number of threads to give the backend. 0 lets the backend use all available cores.")
    parser.add_argument('--tmp_dir', default=None, type=str,
                        help="parent directory for temporary databases and alignment files. Structural searches can need a lot of scratch space, so point this at a large filesystem for big inputs.")
    if include_search_arguments:
        parser.add_argument('--alignment_type', type=int, default=2, choices=[0, 1, 2],
                            help="how to compute the alignment: 0 = 3Di only, 1 = TM-align, 2 = 3Di+AA. Type 1 requires C-alpha coordinates in both databases.")
        parser.add_argument('--metrics', nargs='+', default=None, type=str, choices=list(METRIC_CHOICES),
                            help="additional per-hit structural metrics to compute and store as feature qualifiers. tmscore, lddt and rmsd require C-alpha coordinates in both databases.")
        parser.add_argument('--max_seqs', type=int, default=1000,
                            help="maximum number of hits the backend keeps per reference. Raise this when searching a large database, otherwise coverage is silently capped.")
        parser.add_argument('--device', default="cpu", type=str,
                            help="where to run the search. 'cpu' for a CPU-only search, or a CUDA device string ('cuda', 'cuda:0', 'cuda:1', ...) for a GPU-accelerated search. GPU search requires a recent foldseek. [default: cpu]")


def build_aligner(params) -> StructureAligner:
    """Instantiate the backend selected by the parsed arguments."""
    aligner_class = ALIGNERS[params.algorithm]
    return aligner_class(
        bin_path=getattr(params, f"{params.algorithm}_path", None),
        cpu=getattr(params, "cpu", 0),
        device=getattr(params, "device", None),
        extra_args=getattr(params, "aligner_arg", None),
        tmp_dir=getattr(params, "tmp_dir", None),
    )


class PreparedDatabases(NamedTuple):
    """The databases and scratch space for one search."""
    input_db: StructureDB
    reference_db: StructureDB
    work_dir: str
    input_label: str
    reference_label: str


@contextlib.contextmanager
def prepared_databases(aligner: StructureAligner, input_values: Sequence[str],
                       reference_values: Sequence[str], tmp_dir: Optional[str] = None,
                       keep_db: Optional[str] = None):
    """Resolve -i and -r into databases, building any that are not prebuilt.

    Yields a PreparedDatabases. Databases built here live in a temporary directory that
    is removed on exit, unless keep_db names a prefix for the input database, in which
    case that one is left in place for reuse.
    """
    input_kind, input_resolved = resolve_structure_inputs(input_values)
    reference_kind, reference_resolved = resolve_structure_inputs(reference_values)

    if input_kind == "db" and keep_db is not None:
        raise RuntimeError(
            "--keep_db writes the database built from input structures, but the input is "
            f"already a database ('{input_resolved[0]}')."
        )

    with tempfile.TemporaryDirectory(dir=tmp_dir) as work_dir:
        if input_kind == "db":
            input_db = StructureDB(prefix=input_resolved[0])
        else:
            prefix = keep_db if keep_db is not None else os.path.join(work_dir, "inputdb")
            input_db = aligner.build_db(input_resolved, prefix)

        if reference_kind == "db":
            reference_db = StructureDB(prefix=reference_resolved[0])
        else:
            reference_db = aligner.build_db(reference_resolved, os.path.join(work_dir, "refdb"))

        yield PreparedDatabases(
            input_db=input_db,
            reference_db=reference_db,
            work_dir=work_dir,
            input_label=input_label(input_values, input_kind, input_resolved),
            reference_label=input_label(reference_values, reference_kind, reference_resolved),
        )


def warn_on_saturation(hit_counts: Dict[str, int], max_seqs: int) -> None:
    """Warn when a reference returned exactly --max_seqs hits, i.e. was likely truncated.

    Silent truncation reads as "this reference has few homologs" when it may have many, so
    it must never pass unreported.
    """
    saturated = sorted(name for name, count in hit_counts.items() if count >= max_seqs)
    if saturated:
        warnings.warn(
            f"--max_seqs ({max_seqs}) was reached by {len(saturated)} reference(s): "
            f"{', '.join(saturated[:5])}{' ...' if len(saturated) > 5 else ''}. "
            "Results for those references are truncated; raise --max_seqs for full coverage.",
            RuntimeWarning,
        )


def group_hits_by_target(hits: Iterable[StructureHit], db_name: str, evalue: float,
                         min_evalue: float = 0.0, program: str = "foldseek",
                         strip_extension: bool = True
                         ) -> Iterator[Tuple[str, List[SearchResult], str]]:
    """Stream (target_name, [SearchResult], target_sequence) groups from sorted hits.

    Requires hits sorted by target (see StructureAligner.search(sort_by_target=True)); a
    run of the same target that reappears later would be emitted as a second group.

    The target sequence travels with its group rather than being accumulated in a lookup
    table, so a caller keeping only the best N groups keeps only N sequences.
    """
    current_target = None
    current: List[SearchResult] = []
    current_sequence = ""
    for hit in hits:
        results = structure_hits_to_search_results([hit], db_name, evalue, min_evalue,
                                                   program, strip_extension)
        if not results:
            continue
        if hit.target != current_target:
            if current_target is not None and current:
                yield current_target, current, current_sequence
            current_target = hit.target
            current = []
            current_sequence = ""
        current.extend(results[hit.target])
        if hit.tseq and not current_sequence:
            current_sequence = hit.tseq
    if current_target is not None and current:
        yield current_target, current, current_sequence


def copy_hits_tsv(work_dir: str, destination: str) -> None:
    """Copy the backend's raw hit table out of the scratch directory.

    The table is written by the backend inside work_dir, which is deleted when the search
    finishes, so it has to be copied out while the search context is still open.
    """
    for candidate in ("results.by_target.tsv", "results.tsv"):
        source = os.path.join(work_dir, candidate)
        if os.path.isfile(source):
            shutil.copyfile(source, destination)
            return
    warnings.warn(f"No raw hit table was found in '{work_dir}'; --hits_tsv not written.")


def write_records(records, output, max_output_bytes, output_description, mitigation_options):
    """Write genbank records to a path or stdout, enforcing the output size limit.

    Writes to a temporary file next to the target and renames on success, so a failed or
    size-limited run does not leave a partial output file behind. stdout is written
    directly, since there is nothing to rename.

    An output of None or "-" means stdout, matching the rest of the suite.
    """
    from domainator.utils import write_genbank

    if output is None or output == "-":
        write_genbank(records, sys.stdout, default_molecule_type="protein")
        return

    temp_output_path = make_temporary_output_path(output)
    handle, _ = open_if_is_name_for_write(temp_output_path)
    try:
        written_bytes = 0
        for record in records:
            written_bytes = write_records_with_limit(
                (record,), handle, written_bytes, max_output_bytes,
                output_description, mitigation_options,
            )
    except OutputSizeLimitExceeded as exc:
        handle.close()
        os.remove(temp_output_path)
        raise SystemExit(str(exc)) from None
    except Exception:
        handle.close()
        if os.path.exists(temp_output_path):
            os.remove(temp_output_path)
        raise
    handle.close()
    os.replace(temp_output_path, output)
