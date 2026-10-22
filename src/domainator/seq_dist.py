"""Calculates similarity scores between protein sequences or HMM profiles

Typically the input and reference will be the same sequence file, creating a pairwise similarity matrix. 
An hmm file can also be used as a reference to make a table of profile scores for each input peptide.
An hmm file can also be used as both input and reference to compare profile similarity.

"""
import logging
import warnings

from jsonargparse import ActionConfigFile
import heapq
import re
import shlex
import shutil
import subprocess
import tempfile
import threading
from domainator.cli.common import add_input, reject_stdin, tool_parser
from Bio import SeqIO
from collections import OrderedDict,namedtuple, defaultdict
import numpy as np
import scipy.sparse
import pyhmmer
import tqdm
from typing import Iterator, List, Optional
from domainator.seqio import get_file_type, parse_seqfiles
from domainator.hmmio import open_hmm_file
from domainator._mp import make_pool, resolve_cpu
from domainator.ssn_edges import StreamingMstKnnAccumulator
from domainator.data_matrix import DataMatrix, symmetrize_matrix
from domainator.hmmer_search import compare_hmmer
from domainator.matrix_output import add_matrix_output_arguments, add_symmetrize_argument, resolve_matrix_outputs, write_matrix_outputs
from domainator.errors import ExternalToolError
from domainator.output_guardrails import add_max_output_gb_argument, format_bytes, max_output_gb_to_bytes
from domainator.transform_matrix import MODES, LOG10_2
from domainator.sparsify import (MAX_HITS_PER_QUERY_SPARSIFY_NOTE, add_sparsify_arguments, sparsify_option_label,
                                 sparsify_requested)

logger = logging.getLogger(__name__)

# Modes whose post-transform values are computable per edge (no global matrix maxes),
# and therefore compatible with streaming --mst/--knn pruning.
SPARSIFY_STREAMABLE_MODES = {"score", "bool", "efi_score"}

CmpResult = namedtuple("CmpResult", ["score", "input", "reference"])


def _progress(iterable, enabled, **kwargs):
    if not enabled:
        return iterable
    return tqdm.tqdm(iterable, leave=True, dynamic_ncols=True, **kwargs)


def _store_top_k_result(result_heaps, result, k):
    heap = result_heaps[result.input]
    if len(heap) < k:
        heapq.heappush(heap, result)
    else:
        heapq.heappushpop(heap, result)


def _iter_top_k_results(result_heaps):
    for input_id in result_heaps:
        for result in sorted(result_heaps[input_id], key=lambda item: (-item.score, item.reference)):
            yield result


def _build_sparse_max_matrix(row_indices, col_indices, values, shape):
    if len(values) == 0:
        return scipy.sparse.csr_array(shape, dtype=np.float64)

    row_indices = np.asarray(row_indices, dtype=np.int64)
    col_indices = np.asarray(col_indices, dtype=np.int64)
    values = np.asarray(values, dtype=np.float64)

    order = np.lexsort((col_indices, row_indices))
    row_indices = row_indices[order]
    col_indices = col_indices[order]
    values = values[order]

    starts = np.empty(values.shape[0], dtype=bool)
    starts[0] = True
    starts[1:] = (row_indices[1:] != row_indices[:-1]) | (col_indices[1:] != col_indices[:-1])
    reduce_starts = np.flatnonzero(starts)

    if reduce_starts.size != values.size:
        values = np.maximum.reduceat(values, reduce_starts)
        row_indices = row_indices[reduce_starts]
        col_indices = col_indices[reduce_starts]

    return scipy.sparse.csr_array((values, (row_indices, col_indices)), shape=shape)


class _SparseMaxResultBuilder:
    def __init__(self, query_name_to_idx, db_name_to_idx):
        self.query_name_to_idx = query_name_to_idx
        self.db_name_to_idx = db_name_to_idx
        self.shape = (len(query_name_to_idx), len(db_name_to_idx))
        self.row_indices = []
        self.col_indices = []
        self.values = []

    def _lookup_query_idx(self, query_id):
        try:
            return self.query_name_to_idx[query_id]
        except KeyError as exc:
            raise RuntimeError(f"Unknown query id returned by similarity search: {query_id}") from exc

    def _lookup_subject_idx(self, subject_id):
        try:
            return self.db_name_to_idx[subject_id]
        except KeyError as exc:
            raise RuntimeError(f"Unknown subject id returned by similarity search: {subject_id}") from exc

    def add_result(self, query_id, subject_id, value):
        self.row_indices.append(self._lookup_query_idx(query_id))
        self.col_indices.append(self._lookup_subject_idx(subject_id))
        self.values.append(value)

    def build(self):
        return _build_sparse_max_matrix(self.row_indices, self.col_indices, self.values, self.shape)


class _GroupedQuerySparseMaxResultBuilder(_SparseMaxResultBuilder):
    def __init__(self, query_name_to_idx, db_name_to_idx, stream_name, expected_query_order=None):
        super().__init__(query_name_to_idx, db_name_to_idx)
        self.stream_name = stream_name
        self.expected_query_positions = None
        if expected_query_order is not None:
            self.expected_query_positions = {query_id: idx for idx, query_id in enumerate(expected_query_order)}
        self.last_query_position = None
        self.current_query_id = None
        self.current_query_scores = dict()
        self.closed_queries = set()

    def _start_query(self, query_id):
        self._lookup_query_idx(query_id)
        if query_id in self.closed_queries:
            raise RuntimeError(
                f"{self.stream_name} yielded query '{query_id}' again after moving past it. "
                "This builder assumes results are grouped by query."
            )

        if self.expected_query_positions is not None:
            try:
                query_position = self.expected_query_positions[query_id]
            except KeyError as exc:
                raise RuntimeError(
                    f"{self.stream_name} yielded query '{query_id}', which is not present in the expected query order."
                ) from exc
            if self.last_query_position is not None and query_position < self.last_query_position:
                raise RuntimeError(
                    f"{self.stream_name} yielded query '{query_id}' out of order. "
                    "This builder assumes queries arrive in nondecreasing expected-query order."
                )
            self.last_query_position = query_position

        self.current_query_id = query_id
        self.current_query_scores = dict()

    def _flush_current_query(self):
        if self.current_query_id is None:
            return

        row_idx = self._lookup_query_idx(self.current_query_id)
        for col_idx, value in sorted(self.current_query_scores.items()):
            self.row_indices.append(row_idx)
            self.col_indices.append(col_idx)
            self.values.append(value)

        self.closed_queries.add(self.current_query_id)
        self.current_query_id = None
        self.current_query_scores = dict()

    def add_result(self, query_id, subject_id, value):
        if self.current_query_id is None:
            self._start_query(query_id)
        elif query_id != self.current_query_id:
            self._flush_current_query()
            self._start_query(query_id)

        subject_idx = self._lookup_subject_idx(subject_id)
        previous_value = self.current_query_scores.get(subject_idx)
        if previous_value is None or value > previous_value:
            self.current_query_scores[subject_idx] = value

    def build(self):
        self._flush_current_query()
        return super().build()


def _drain_stream(stream, sink):
    for line in stream:
        sink.append(line)


def _advance_query_progress(progress_bar, query_positions, current_query_idx, query_id):
    query_idx = query_positions.get(query_id)
    if query_idx is None:
        return current_query_idx

    if current_query_idx is None:
        progress_bar.update(query_idx + 1)
        return query_idx

    if query_idx > current_query_idx:
        progress_bar.update(query_idx - current_query_idx)
        return query_idx

    return current_query_idx


def _finish_query_progress(progress_bar, current_query_idx, total_queries):
    if current_query_idx is None:
        progress_bar.update(total_queries)
    else:
        remaining = total_queries - current_query_idx - 1
        if remaining > 0:
            progress_bar.update(remaining)
    progress_bar.close()

# A large square all-vs-all search, at a sensitivity where diamond collapses its index
# to a single chunk, is the configuration that exhausts scratch space. Warn about it,
# but do not retune it: peak scratch is not predictable from the input size.
PREFLIGHT_SEQUENCE_THRESHOLD = 50_000
# Sensitivities where diamond auto-selects --index-chunks 1 (measured, diamond v2.2.5.185).
SINGLE_INDEX_CHUNK_MODES = frozenset({"vs", "us"})
# --diamond_args flags that mean the caller has already bounded diamond's scratch themselves.
_SCRATCH_TUNING_FLAGS = frozenset({"-b", "--block-size", "-c", "--index-chunks"})


def build_diamond_preflight_warning(mode, n_sequences, total_letters, square,
                                    tmpdir=None, diamond_params=None):
    """Return an advisory for a large all-vs-all diamond search, or None.

    This only advises. It deliberately does not adjust diamond's options: peak scratch
    demand depends on the redundancy of the input and on how much of the filesystem
    other users take during the run, neither of which is knowable up front.
    """
    if not square or n_sequences < PREFLIGHT_SEQUENCE_THRESHOLD:
        return None
    if mode not in SINGLE_INDEX_CHUNK_MODES:
        return None
    if diamond_params and _SCRATCH_TUNING_FLAGS.intersection(diamond_params):
        return None  # caller has already bounded it

    lines = [
        "WARNING: large all-vs-all diamond search; this can exhaust scratch space.",
        f"  {n_sequences:,} sequences ({total_letters:,} aa) compared against themselves,",
        f"  at --algorithm diamond_{mode}, where diamond uses a single index chunk.",
    ]
    if tmpdir is not None:
        free = _describe_free_space(tmpdir)
        lines.append(f"  Scratch directory: {tmpdir}" + (f" ({free})" if free else ""))
        lines.append("  Set TMPDIR to place it on another filesystem. diamond unlinks this")
        lines.append("  scratch while holding it open, so `du` and `ls` will not show it;")
        lines.append("  only the filesystem's free space reflects it.")
    lines.append("  Peak scratch can reach the terabyte range for a redundant input, and is")
    lines.append("  not tuned automatically because it cannot be predicted from input size.")
    lines.append("  To bound it, pass block/index-chunk options through --diamond_args, e.g.:")
    lines.append("      --diamond_args=\"-b 0.07 -c 4\"")
    lines.append("  A smaller -b lowers peak scratch further; -c splits the reference index")
    lines.append("  into that many passes. Neither changes the resulting matrix.")
    lines.append("  Lowering --algorithm (diamond_s, diamond_mids) reduces the work most, but")
    lines.append("  does change which hits are found.")
    return "\n".join(lines)


# Substrings diamond writes when a run dies for a reason the caller can act on.
# Matched case-insensitively against the captured stderr.
_DIAMOND_DISK_PATTERNS = ("no space left on device", "error writing to file",
                          "write error in", "disk full", "enospc")
_DIAMOND_MEMORY_PATTERNS = ("std::bad_alloc", "failed to allocate", "out of memory",
                            "cannot allocate memory", "bad_alloc")
_DIAMOND_OPTION_PATTERNS = ("invalid option", "unknown option", "unrecognized option")
_DIAMOND_INPUT_PATTERNS = ("error opening file", "no such file or directory",
                           "invalid input file format", "empty database",
                           "error detecting input file format", "seems to be empty",
                           "invalid character", "unrecognized sequence type")

_DISK_MITIGATIONS = (
    "Lower the sensitivity: --algorithm diamond_s (16 seed shapes) or diamond_mids (8), "
    "instead of diamond_us (64). This is usually the single largest reduction.",
    "Shrink diamond's query blocks so less scratch is live at once, e.g. "
    "--diamond_args=\"-b 0.01 -c 4\". Peak scratch scales roughly with -b.",
    "Put the scratch on a larger filesystem: set TMPDIR, or "
    "--diamond_args=\"--tmpdir /path/with/space\".",
    "Consider fewer input sequences (for example dereplicating first); the seed-hit "
    "volume grows quadratically with the size of a redundant input.",
)
_MEMORY_MITIGATIONS = (
    "Shrink diamond's query blocks: --diamond_args=\"-b 0.01 -c 4\".",
    "Lower the sensitivity: --algorithm diamond_s or diamond_mids.",
    "Use fewer threads (--cpu N); diamond allocates per-thread buffers.",
)


def _collapse_stderr(text, max_lines=15):
    """Collapse runs of identical lines (diamond repeats write errors once per thread)."""
    collapsed = []
    for line in (candidate.rstrip() for candidate in text.splitlines()):
        if not line:
            continue
        if collapsed and collapsed[-1][0] == line:
            collapsed[-1][1] += 1
        else:
            collapsed.append([line, 1])
    rendered = [f"{line} (x{count})" if count > 1 else line for line, count in collapsed]
    if len(rendered) > max_lines:
        hidden = len(rendered) - max_lines
        rendered = [f"... {hidden} earlier line(s) omitted ..."] + rendered[-max_lines:]
    return "\n".join("    " + line for line in rendered)


def _describe_free_space(path):
    """Return a human-readable free-space note for the filesystem holding path."""
    try:
        usage = shutil.disk_usage(path)
    except OSError:
        return None
    return f"{format_bytes(usage.free)} free of {format_bytes(usage.total)}"


def build_diamond_failure_message(stage, returncode, stderr_text, input_fasta,
                                  tmpdir=None, command=None, diamond_params=None):
    """Explain a failed diamond run and suggest what the caller can change.

    Recognizes the failure modes that are actionable from the seq_dist command line
    (scratch disk exhaustion, memory exhaustion, a rejected --diamond_args flag, unreadable
    input) and falls back to the raw stderr when the cause is not one of those.
    """
    haystack = stderr_text.lower()
    lines = [f"diamond {stage} failed (exit code {returncode}) on: {input_fasta}"]

    def _add_mitigations(mitigations):
        lines.append("")
        lines.append("  Suggested mitigations, roughly most effective first:")
        for number, text in enumerate(mitigations, start=1):
            lines.append(f"    {number}. {text}")

    if any(pattern in haystack for pattern in _DIAMOND_DISK_PATTERNS):
        lines.append("")
        lines.append("  Cause: diamond ran out of disk space writing its seed-search scratch.")
        if tmpdir is not None:
            free = _describe_free_space(tmpdir)
            lines.append(f"  Scratch directory: {tmpdir}" + (f" ({free})" if free else ""))
            lines.append("  Note: diamond unlinks these files while holding them open, so they")
            lines.append("        never appear in `du` or `ls` -- only free space reflects them.")
        lines.append("  This happens before diamond emits any results, so --mst/--knn")
        lines.append("  streaming cannot reduce it: all seed hits for a query block are")
        lines.append("  materialized before that block produces any output.")
        _add_mitigations(_DISK_MITIGATIONS)
    elif returncode in (-9, 137) or any(pattern in haystack for pattern in _DIAMOND_MEMORY_PATTERNS):
        lines.append("")
        if returncode in (-9, 137):
            lines.append("  Cause: diamond was killed with SIGKILL, which usually means the")
            lines.append("         out-of-memory killer stopped it. Check `dmesg` to confirm.")
        else:
            lines.append("  Cause: diamond could not allocate memory.")
        _add_mitigations(_MEMORY_MITIGATIONS)
    elif any(pattern in haystack for pattern in _DIAMOND_OPTION_PATTERNS):
        lines.append("")
        lines.append("  Cause: diamond rejected one of its command line options.")
        if diamond_params:
            rendered = " ".join(flag if value is None else f"{flag} {value}" for flag, value in diamond_params.items())
            lines.append(f"  Options passed via --diamond_args: {rendered}")
            lines.append("  Check these against `diamond blastp --help`; --diamond_args is forwarded as-is.")
        else:
            lines.append("  No --diamond_args were supplied, so this may be a diamond version mismatch.")
    elif any(pattern in haystack for pattern in _DIAMOND_INPUT_PATTERNS):
        lines.append("")
        lines.append("  Cause: diamond could not read an input or database file.")
        lines.append("  Check that the file is non-empty, readable, and valid protein FASTA.")

    if command:
        lines.append("")
        lines.append("  Command: " + " ".join(command))
    lines.append("")
    lines.append("  diamond stderr:")
    lines.append(_collapse_stderr(stderr_text) or "    (no stderr captured)")
    return "\n".join(lines)


# diamond options that seq_dist sets itself. Overriding these through --diamond_args would
# either break parsing of the result stream or silently contradict a seq_dist option,
# so they are rejected with a pointer to the option that does control them.
RESERVED_DIAMOND_PARAMS = {
    "-q": "-i/--input",
    "--query": "-i/--input",
    "-d": "-r/--reference",
    "--db": "-r/--reference",
    "-f": "the fixed tabular output format required to parse the result stream",
    "--outfmt": "the fixed tabular output format required to parse the result stream",
    "-p": "--cpu",
    "--threads": "--cpu",
    "-k": "-k",
    "--max-target-seqs": "-k",
    "--max-hsps": "seq_dist (fixed at 1)",
    "--min-score": "--threshold",
    "--verbose": "seq_dist (required for progress reporting)",
    "--fast": "--algorithm",
    "--mid-sensitive": "--algorithm",
    "--sensitive": "--algorithm",
    "--more-sensitive": "--algorithm",
    "--very-sensitive": "--algorithm",
    "--ultra-sensitive": "--algorithm",
}

# diamond options seq_dist supplies a default for, but which --diamond_args may replace.
def _default_overridable_diamond_params(tmpdir):
    return OrderedDict((("--algo", "0"), ("--tmpdir", tmpdir)))


def parse_diamond_args(args_str):
    """Parse a --diamond_args string into an ordered dict of {flag: value}.

    The string is split like a shell command line, for example: "-c 4 -b 0.4 --target-indexed".
    A flag followed by another flag (or by nothing) takes no argument and maps to None.
    """
    if not args_str:
        return OrderedDict()
    try:
        tokens = shlex.split(args_str)
    except ValueError as exc:
        raise ValueError(f"Could not split --diamond_args into options: {exc}") from exc

    def is_flag(token):
        return token.startswith("-") and not _NUMBER_RE.fullmatch(token)

    parsed = OrderedDict()
    i = 0
    while i < len(tokens):
        flag = tokens[i]
        if not is_flag(flag):
            raise ValueError(f"--diamond_args must be diamond options starting with '-', but got: {flag!r}")
        if flag in RESERVED_DIAMOND_PARAMS:
            raise ValueError(
                f"--diamond_args may not set '{flag}', because seq_dist controls it via "
                f"{RESERVED_DIAMOND_PARAMS[flag]}."
            )
        if i + 1 < len(tokens) and not is_flag(tokens[i + 1]):
            parsed[flag] = tokens[i + 1]
            i += 2
        else:
            parsed[flag] = None
            i += 1
    return parsed


_NUMBER_RE = re.compile(r"-?\d+(\.\d*)?([eE][-+]?\d+)?")


def add_params_to_args_list(arg_list, params):
    """Append {flag: value} pairs to a command line. A null value emits a bare flag."""
    for flag, value in params.items():
        arg_list.append(flag)
        if value is not None:
            arg_list.append(str(value))


def diamond(input_fasta, reference_fasta, max_target_seqs, threads, tmpdir, mode, max_hsps=1, min_score: float = 0.0, query_order=None, progress: bool = False,
            diamond_params=None, diamond_path=None) -> Iterator[CmpResult]:
    """
        runs a search using diamond, returns a list of tuples of (inputid, referenceid, score)

        diamond_path: the diamond executable, if it is not "diamond" on the PATH.
    """
    diamond_bin = diamond_path or "diamond"
    mode_to_arg = { 
                    "f": ["--fast"], 
                    "d":[],
                    "mids":["--mid-sensitive"], 
                    "s": ["--sensitive"], 
                    "mors":["--more-sensitive"], 
                    "vs": ["--very-sensitive"], 
                    "us": ["--ultra-sensitive"], 
                    }

    dbpath =  tmpdir+"/"+"db"
    mkdb_command = [diamond_bin, "makedb", "--in", reference_fasta, "-d", dbpath]
    try:
        mkdb_result = subprocess.run(mkdb_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                encoding='utf-8')
    except FileNotFoundError:
        raise ExternalToolError(f"Could not run '{diamond_bin}'. It is not on PATH.\n"
                                "  Install it (conda install -c bioconda diamond), activate the environment\n"
                                "  that provides it, or give its location with --diamond_path, then re-run.") from None

    if mkdb_result.returncode != 0:
        raise ExternalToolError(build_diamond_failure_message(
            "makedb", mkdb_result.returncode,
            (mkdb_result.stderr or "") + (mkdb_result.stdout or ""),
            reference_fasta, tmpdir=tmpdir, command=mkdb_command))

    dmnd_command = [ "stdbuf", "-oL", "-eL",
        diamond_bin, "blastp", "-q", input_fasta, "-d", dbpath, "-p", str(threads),
        "--outfmt", "6", "qseqid", "sseqid", "bitscore",
        "--max-hsps", str(max_hsps), "--verbose"
    ]
    if max_target_seqs is not None and max_target_seqs > 0:
        dmnd_command.extend(["--max-target-seqs", str(max_target_seqs)])
    else:
        dmnd_command.extend(["--max-target-seqs", str(0)])
    if min_score > 0:
        dmnd_command.extend(["--min-score", str(min_score)])
    dmnd_command.extend(mode_to_arg[mode])

    # Defaults seq_dist supplies but --diamond_args may replace. --tmpdir defaults to the
    # managed temporary directory so that diamond's (large, unlinked-while-open)
    # search scratch lands next to the database instead of in the working directory.
    overridable_params = _default_overridable_diamond_params(tmpdir)
    if diamond_params:
        overridable_params.update(diamond_params)
    add_params_to_args_list(dmnd_command, overridable_params)

    # Echo the exact command so a run can be reproduced or retuned straight from a log.
    logger.info("Running diamond:\n  " + " ".join(dmnd_command))

    try:
        dmnd_process = subprocess.Popen(
            dmnd_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding='utf-8',
            bufsize=1,
        )
    except FileNotFoundError as exc:
        missing = "stdbuf" if "stdbuf" in str(exc) else diamond_bin
        raise ExternalToolError(f"Could not run '{missing}'. It is not on PATH.\n"
                                f"  Command: {' '.join(dmnd_command)}") from None

    stderr_lines = list()
    stderr_thread = None
    if dmnd_process.stderr is not None:
        stderr_thread = threading.Thread(target=_drain_stream, args=(dmnd_process.stderr, stderr_lines), daemon=True)
        stderr_thread.start()

    query_progress = None
    current_query_idx = None
    if progress and query_order is not None:
        query_progress = tqdm.tqdm(total=len(query_order), desc="diamond queries", leave=True, dynamic_ncols=True)
        query_positions = {query_id: idx for idx, query_id in enumerate(query_order)}
    else:
        query_positions = None

    try:
        if dmnd_process.stdout is None:
            raise RuntimeError("diamond did not provide a stdout stream")

        for line in dmnd_process.stdout:
            line = line.strip()
            parts = line.split("\t")
            if len(parts) == 3:
                if query_progress is not None:
                    current_query_idx = _advance_query_progress(query_progress, query_positions, current_query_idx, parts[0])
                yield CmpResult(float(parts[2]), parts[0], parts[1])
    finally:
        if dmnd_process.stdout is not None:
            dmnd_process.stdout.close()
        returncode = dmnd_process.wait()
        if stderr_thread is not None:
            stderr_thread.join()
        if query_progress is not None:
            _finish_query_progress(query_progress, current_query_idx, len(query_order))

    if returncode != 0:
        raise ExternalToolError(build_diamond_failure_message(
            "blastp", returncode, "".join(stderr_lines), input_fasta,
            tmpdir=tmpdir, command=dmnd_command, diamond_params=diamond_params))

def _seq_name_to_idx_dict(names_and_lengths, input_path):
    """name_to_idx, idx_to_name and idx_to_len for sequences, keeping the first of duplicate names."""
    name_to_idx = OrderedDict()
    idx_to_name = list()
    idx_to_len = list()
    for name, length in names_and_lengths:
        if name in name_to_idx:
            warnings.warn(f"Warning: duplicate sequence name found in input {input_path}: {name}")
        else:
            name_to_idx[name] = len(idx_to_name)
            idx_to_len.append(length)
            idx_to_name.append(name)
    return name_to_idx, idx_to_name, idx_to_len


def make_seq_name_to_idx_dict(input_path):
    return _seq_name_to_idx_dict(((rec.id, len(rec.seq)) for rec in parse_seqfiles([input_path])), input_path)

def make_hmm_name_to_idx_dict_from_path(hmm_path):
    name_to_idx = OrderedDict()
    idx_to_name = list()
    idx_to_len = list()

    idx = 0
    for rec in open_hmm_file(hmm_path):
        name = rec.name
        if name in name_to_idx:  
            warnings.warn(f"Warning: duplicate hmm name found: {name}")
        else:
            name_to_idx[name] = idx
            idx_to_name.append(name)
            idx_to_len.append(rec.M)
            idx += 1
    
    return name_to_idx, idx_to_name, idx_to_len

def make_hmm_name_to_idx_dict(hmm_list):
    name_to_idx = OrderedDict()
    idx_to_name = list()
    idx_to_len = list()

    idx = 0
    for rec in hmm_list:
        name = rec.name
        if name in name_to_idx:  
            warnings.warn(f"Warning: duplicate hmm name found: {name}")
        else:
            name_to_idx[name] = idx
            idx_to_name.append(name)
            idx_to_len.append(rec.M)
            idx += 1
    
    return name_to_idx, idx_to_name, idx_to_len

def run_hmmsearch(hmmer_seqs, hmmer_models, k, threads, search_type, min_score: float = 0.0, progress: bool = False, query_count=None):
    if search_type == "hmmsearch":
        hits = pyhmmer.hmmer.hmmsearch(hmmer_models, hmmer_seqs, cpus=threads)
    elif search_type == "phmmer":
        hits = pyhmmer.hmmer.phmmer(hmmer_models, hmmer_seqs, cpus=threads)
    else:
        raise ValueError(f"search_type {search_type} not recognized")

    query_hits = _progress(hits, progress, total=query_count, desc=f"{search_type} queries")
    top_k_by_input = defaultdict(list)
    for top_hits in query_hits:
        for hit in top_hits:
            name = hit.name
            domain = hit.best_domain.alignment.hmm_name
            score = round(hit.score,2)
            if score <= min_score:
                continue
            result = CmpResult(score, name, domain)
            if k is None:
                yield result
            else:
                _store_top_k_result(top_k_by_input, result, k)

    del hits
    if k is not None:
        yield from _iter_top_k_results(top_k_by_input)


class _run_hmmer_compare_worker():
    def __init__(self, hmmer_targets, k=None, min_score: float = 0.0):
        self.k = k
        self.hmmer_targets = hmmer_targets
        self.min_score = min_score
        # Parsed on first use inside the worker process, not here: pyhmmer objects
        # are not picklable, and this instance is sent to every worker. Caching it
        # trades memory (the whole target set stays resident per worker) for not
        # re-parsing the target file once per query profile.
        self._targets = None

    def _target_profiles(self):
        if self._targets is None:
            self._targets = list(open_hmm_file(self.hmmer_targets))
        return self._targets

    def __call__(self, input_profile):
        out_heap = []
        for target_profile in self._target_profiles():
            
            score, _traceback, _max_index, _match_scores = compare_hmmer(input_profile, target_profile)
            score = round(score,2)
            if score <= self.min_score:
                continue
            result = CmpResult(score, input_profile.name, target_profile.name)
            if (self.k is None) or (len(out_heap) < self.k):
                heapq.heappush(out_heap, result)
            else:
                heapq.heappushpop(out_heap, result)
        return out_heap

def run_hmmer_compare(hmmer_queries, hmmer_targets, k, threads, min_score: float = 0.0, progress: bool = False, query_count=None):
    
    worker = _run_hmmer_compare_worker(hmmer_targets, k, min_score=min_score)
    with make_pool(processes=threads) as pool:
        for hits in _progress(pool.imap_unordered(worker, open_hmm_file(hmmer_queries)), progress, total=query_count, desc="hmmer_compare queries"):
            for hit in sorted(hits, key=lambda item: (-item.score, item.reference)):
                yield hit

def convert_to_fasta(input_path, output_path):
    """Write input_path as FASTA, and return its make_seq_name_to_idx_dict (from the same pass)."""
    def records():
        with open(output_path, "w") as outh:
            for rec in parse_seqfiles([input_path]):
                SeqIO.write(rec, outh, "fasta")
                yield rec.id, len(rec.seq)
    return _seq_name_to_idx_dict(records(), input_path)


def _fasta_and_index(input_path, input_format, output_path):
    """A FASTA path for input_path (itself if it is FASTA, else a conversion written to output_path),
    and its make_seq_name_to_idx_dict."""
    if input_format == "fasta":
        return input_path, make_seq_name_to_idx_dict(input_path)
    return output_path, convert_to_fasta(input_path, output_path)


def _make_result_builder(algorithm, query_name_to_idx, db_name_to_idx, search_type=None, query_order=None):
    if algorithm.startswith("diamond"):
        return _GroupedQuerySparseMaxResultBuilder(
            query_name_to_idx,
            db_name_to_idx,
            stream_name="diamond result stream",
            expected_query_order=query_order,
        )

    if algorithm == "hmmer_compare":
        return _GroupedQuerySparseMaxResultBuilder(
            query_name_to_idx,
            db_name_to_idx,
            stream_name="hmmer_compare result stream",
        )

    if algorithm == "hmmer" and search_type == "phmmer":
        return _GroupedQuerySparseMaxResultBuilder(
            query_name_to_idx,
            db_name_to_idx,
            stream_name="phmmer result stream",
            expected_query_order=query_order,
        )

    return _SparseMaxResultBuilder(query_name_to_idx, db_name_to_idx)


def seq_dist(input_path: str, *, input_format: Optional[str] = None, reference: Optional[str] = None,
             reference_format: Optional[str] = None, max_hits_per_query: Optional[int] = None,
             algorithm: str = "diamond_us", mode: str = "score", cpu: int = 0, threshold: float = 0,
             progress: bool = False, mst: bool = False, knn: Optional[int] = None, diamond_args: Optional[List[str]] = None,
             symmetrize: str = "none", diamond_path: Optional[str] = None) -> DataMatrix:
    """Compare every sequence or profile in ``input_path`` with every one in ``reference``.

    Returns a DataMatrix of ``mode`` values with a row per input and a column per reference;
    pass the same path as ``input_path`` and ``reference`` for a self-comparison. The options
    mean the same as the ``domainator seq_dist`` options with the same names; ``diamond_args``
    is the parsed form of ``--diamond_args`` (see :func:`parse_diamond_args`). A format of None
    is read from the file name; a ``reference`` of None compares ``input_path`` with itself; ``cpu``
    0 means one thread per physical core.
    """
    if input_format is None:
        input_format = get_file_type(input_path)
    if reference is None:
        reference, reference_format = input_path, input_format
    elif reference_format is None:
        reference_format = get_file_type(reference)
    cpu = resolve_cpu(cpu)
    sparsify = sparsify_requested(mst, knn)
    if sparsify:
        option = sparsify_option_label(mst, knn)
        if mode not in SPARSIFY_STREAMABLE_MODES:
            raise ValueError(
                f"{option} is only supported with --mode {sorted(SPARSIFY_STREAMABLE_MODES)} "
                f"(modes that normalize against global matrix maxima cannot be streamed). "
                f"For '{mode}', run 'seq_dist --mode score --hdf_out out.hdf5' then "
                f"'transform_matrix -i out.hdf5 --mode {mode} {option} --hdf_out pruned.hdf5'."
            )
        if reference != input_path:
            raise ValueError(
                f"{option} requires a square, symmetric comparison (reference must be the same file as input)."
            )
        if symmetrize != "none":
            raise ValueError(f"{option} output is always max-symmetric, so --symmetrize {symmetrize} cannot be applied on top of it.")
    if symmetrize != "none" and reference != input_path:
        raise ValueError("--symmetrize requires a square, symmetric comparison (reference must be the same file as input).")

    if input_format == "hmm" and reference_format == "hmm":
        if algorithm != "hmmer_compare":
            algorithm = "hmmer_compare"
            warnings.warn(f"Input and reference are both hmm files, but algorithm is not hmmer_compare. Changing algorithm to hmmer_compare.")
    search_type = None
    with tempfile.TemporaryDirectory() as tmpdir:
        # each of these needs to set query_names, db_names, and results_table
        
        if algorithm.startswith("diamond"):
            parts = algorithm.split("_")

            # convert input and db to fasta if not fasta
            dmnd_input, (query_name_to_idx, query_idx_to_name, query_idx_to_len) = _fasta_and_index(input_path, input_format, tmpdir + "/input.fasta")
            db_name_to_idx, db_idx_to_name, db_idx_to_len = query_name_to_idx, query_idx_to_name, query_idx_to_len
            if input_path == reference:
                dmnd_ref = dmnd_input
            else:
                dmnd_ref, (db_name_to_idx, db_idx_to_name, db_idx_to_len) = _fasta_and_index(reference, reference_format, tmpdir + "/ref.fasta")
            preflight = build_diamond_preflight_warning(
                parts[1], len(query_name_to_idx), int(sum(query_idx_to_len)),
                square=(input_path == reference),
                tmpdir=tmpdir, diamond_params=diamond_args,
            )
            if preflight is not None:
                logger.warning(preflight)

            # call diamond
            
            results_iter = diamond(dmnd_input, dmnd_ref, max_hits_per_query, cpu, tmpdir, parts[1], min_score=threshold, query_order=query_idx_to_name, progress=progress, diamond_params=diamond_args, diamond_path=diamond_path)
            score_progress = False
        
        elif algorithm == "hmmer":
            hmmer_seqfile = input_path
            if input_format != "fasta":
                hmmer_seqfile = tmpdir + "/input.fasta"
                convert_to_fasta(input_path, hmmer_seqfile)

            # The hits are named by easel, so the matrix labels come from easel's names too.
            with pyhmmer.easel.SequenceFile(hmmer_seqfile, digital=True) as seq_file:
                hmmer_seqs = list(seq_file)
            query_name_to_idx, query_idx_to_name, query_idx_to_len = _seq_name_to_idx_dict(((seq.name, len(seq)) for seq in hmmer_seqs), input_path)


            if reference_format == "hmm":
                db_data = list(open_hmm_file(reference))
                db_name_to_idx, db_idx_to_name, db_idx_to_len = make_hmm_name_to_idx_dict(db_data)
                search_type = "hmmsearch"
            else: # a sequence type
                search_type = "phmmer"
                hmmer_dbfile = reference
                if reference_format != "fasta":
                    hmmer_dbfile = tmpdir + "/ref.fasta"
                    convert_to_fasta(reference, hmmer_dbfile)
                with pyhmmer.easel.SequenceFile(hmmer_dbfile, digital=True) as seq_file:
                    db_data = list(seq_file)
                db_name_to_idx, db_idx_to_name, db_idx_to_len = _seq_name_to_idx_dict(((seq.name, len(seq)) for seq in db_data), reference)
                       
            results_iter = run_hmmsearch(hmmer_seqs, db_data, max_hits_per_query, cpu, search_type, min_score=threshold, progress=progress, query_count=len(db_data))
            score_progress = False
        
        
        elif algorithm == "hmmer_compare":
            query_name_to_idx, query_idx_to_name, query_idx_to_len = make_hmm_name_to_idx_dict_from_path(input_path)
            db_name_to_idx, db_idx_to_name, db_idx_to_len = query_name_to_idx, query_idx_to_name, query_idx_to_len
            if input_path != reference:
                db_name_to_idx, db_idx_to_name, db_idx_to_len = make_hmm_name_to_idx_dict_from_path(reference)
            results_iter = run_hmmer_compare(input_path, reference, max_hits_per_query, cpu, min_score=threshold, progress=progress, query_count=len(query_name_to_idx))
            score_progress = False

        else: 
            raise RuntimeError(f"Algorithm not recognized: {algorithm}")

        

        if sparsify:
            # Stream edges straight into the MST and/or kNN accumulator with
            # bounded memory; the full all-vs-all matrix is never materialized. mode value
            # is computed per edge (no global maxima needed for score/efi_score). bool
            # selects on the raw score, since all-1 values would leave nothing to rank, and
            # becomes a 0/1 mask afterward. The raw-score --threshold gate runs first (matching the
            # help text), then the accumulator prunes on values > 0.
            log_lengths = np.log10(np.asarray(query_idx_to_len, dtype=np.float64)) if mode == "efi_score" else None
            accumulator = StreamingMstKnnAccumulator(len(query_name_to_idx), mst=mst, knn=knn, lower_bound=0)
            for (score, query_id, subject_id) in _progress(results_iter, score_progress, desc="Scoring pairs"):
                if score <= threshold:
                    continue
                query_idx = query_name_to_idx[query_id]
                subject_idx = db_name_to_idx[subject_id]
                if mode == "efi_score":
                    value = score * LOG10_2 - (log_lengths[query_idx] + log_lengths[subject_idx])
                    if value <= 0:
                        continue
                else:  # score, bool
                    value = score
                accumulator.add_edge(query_idx, subject_idx, value)
            matrix = accumulator.to_csr()
            if mode == "bool":
                matrix.data[:] = 1
        else:
            result_builder = _make_result_builder(
                algorithm,
                query_name_to_idx,
                db_name_to_idx,
                search_type=search_type,
                query_order=query_idx_to_name,
            )

            for (score, query_id, subject_id) in _progress(results_iter, score_progress, desc="Scoring pairs"):
                if score <= threshold:
                    continue
                matrix_score = 1 if mode == "bool" else score
                result_builder.add_result(query_id, subject_id, matrix_score)

            matrix = symmetrize_matrix(result_builder.build(), symmetrize)

    query_idx_to_len = np.array(query_idx_to_len)
    db_idx_to_len = np.array(db_idx_to_len)
    if progress:
        logger.info("Transforming matrix...")
    if not sparsify and MODES[mode] is not None:
        matrix = MODES[mode](matrix, query_idx_to_len, db_idx_to_len)


    return DataMatrix(matrix, query_idx_to_name, db_idx_to_name, query_idx_to_len, db_idx_to_len, data_type=mode, copy_data=False)

def build_parser():
    parser = tool_parser(__doc__)
    add_input(parser, "Input query file. Queries in the input file will be compared to the reference sequences/profiles. This can be a fasta, genbank, or hmm file.", multiple=False, stdin=False)

    parser.add_argument('--input_format', type=str, required=False, default=None, choices=["fasta", "genbank", "hmm"],
                        help="File format of the input, if it cannot be told from the file extension.")

    parser.add_argument('-r', "--reference", type=str, required=False,
                        help="Reference file. If not provided, the input file will be used as the reference, and the output will be a pairwise similarity matrix between the input sequences/profiles. This can be a fasta, genbank, or hmm file.")

    parser.add_argument('--reference_format', type=str, required=False, default=None, choices=["fasta", "genbank", "hmm"],
                        help="File format of the reference, if it cannot be told from the file extension. Ignored for some kinds of algorithms, for example, hmmer expects a .hmm reference file.")

    add_matrix_output_arguments(parser, noun="similarity (or, for *_dist modes, distance) matrix")

    parser.add_argument('--threshold', default=0, type=float, required=False,
                        help="Round any scores lower than this down to zero (this applies to raw scores not processed/normalized).")


    parser.add_argument('--max_hits_per_query', type=int, required=False, default=None,
                        help="Keep at most this many hits for each input sequence, the highest-scoring ones. This caps the search itself (diamond --max-target-seqs, or the hmmer hit list). " + MAX_HITS_PER_QUERY_SPARSIFY_NOTE + " Default: keep all hits.")

    parser.add_argument('--algorithm', type=str, required=False, default="diamond_us", choices=["diamond_f", "diamond_d", "diamond_mids", "diamond_s", "diamond_mors", "diamond_vs", "diamond_us", "hmmer", "hmmer_compare"],
                        help="Which distance metric to use, diamond_*; f: fast, d: default, mids: mid-sensitive, s: sensitive, mors: more-sensitive, vs: very-sensitive, us: ultra-sensitive. default: diamond_us (ultra sensitive)")

    parser.add_argument("--diamond_args", default=None,
                        help="Extra options for diamond blastp, as one quoted string split like a shell command line. Use the = form: "
                             '--diamond_args="-c 4 -b 0.4 --target-indexed". '
                             "Options that seq_dist sets itself (query, database, --outfmt, -p, --max-target-seqs, --max-hsps, --min-score, --verbose, and the sensitivity flags) "
                             "are rejected; use the corresponding seq_dist option instead. "
                             "seq_dist defaults --algo to 0 and --tmpdir to its own managed temporary directory; both can be replaced here. "
                             "Useful for large all-vs-all runs, where diamond's seed-search scratch can exceed available disk: "
                             "a smaller --block-size (-b) and/or more --index-chunks (-c) bound its peak scratch usage.")

    parser.add_argument("--diamond_path", type=str, default=None,
                        help="path to the diamond executable, if it is not on the PATH.")

    add_sparsify_arguments(parser, requirements="Computed while streaming, so the full matrix is never held in memory. Requires the reference to be the same file as the input (a square, symmetric comparison) and --mode in {score, bool, efi_score}; bool selects on the scores and then keeps the selected cells as 1. Best paired with --hdf_out.")

    add_symmetrize_argument(parser, extra_help="Applied to the raw scores, after --threshold and --max_hits_per_query. Requires the reference to be the same file as the input. Cannot be combined with --mst/--knn, whose output is always max-symmetric.")

    parser.add_argument("--mode", type=str, required=False, default="score", choices=list(MODES),
                        help="what kind of values should be in the matrix. score: raw score, bool: 1 if a hit otherwise 0, score_dist: 1 - (score / min(row_max, col_max)), norm_score: score/min(row_max, col_max), efi_score: -log10[2^(-score) * (input_seq_length * reference_seq_length)], efi_score_dist: 1 - (efi_score / min(row_max, col_max)).")

    parser.add_argument('--cpu', type=int, default=0, required=False,
                        help="how many cpu threads to use. 0: one per physical core.")

    parser.add_argument('--progress', action='store_true',
                        help="Show a progress bar for long-running steps.")

    add_max_output_gb_argument(parser)

    parser.add_argument('--config', action=ActionConfigFile)
    return parser


def main(argv):
    #TODO: support for nucleotide sequences? (see domainator compare_contigs), can be implemented easily with nhmmer, without adding new dependencies.
    parser = build_parser()
    params = parser.parse_args(argv)

    ### make output directory

    reject_stdin(parser, [params.input, params.reference] if params.reference else [params.input])
    hdf, dense_text = resolve_matrix_outputs(params)

    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)

    diamond_params = parse_diamond_args(params.diamond_args)
    if diamond_params and not params.algorithm.startswith("diamond"):
        raise ValueError(f"--diamond_args only applies to the diamond_* algorithms, but --algorithm is '{params.algorithm}'.")

    matrix = seq_dist(params.input, input_format=params.input_format, reference=params.reference,
                      reference_format=params.reference_format, max_hits_per_query=params.max_hits_per_query,
                      algorithm=params.algorithm, mode=params.mode, cpu=params.cpu,
                      threshold=params.threshold, progress=params.progress, mst=params.mst, knn=params.knn,
                      diamond_args=diamond_params, symmetrize=params.symmetrize, diamond_path=params.diamond_path)
    write_matrix_outputs(matrix, hdf=hdf, dense_text=dense_text, max_output_bytes=max_output_bytes,
                         mitigation_options=["--max_hits_per_query", "--threshold", "--mst", "--knn"])
