"""Calculates similarity scores between protein sequences or HMM profiles

Typically the input and reference will be the same sequence file, creating a pairwise similarity matrix. 
An hmm file can also be used as a reference to make a table of profile scores for each input peptide.
An hmm file can also be used as both input and reference to compare profile similarity.

"""
import warnings
warnings.filterwarnings("ignore", module='numpy')

from jsonargparse import ArgumentParser, ActionConfigFile
import heapq
import json
import shutil
import sys
import subprocess
import tempfile
import threading
from domainator.Bio import SeqIO
from collections import OrderedDict,namedtuple, defaultdict
import numpy as np
import scipy.stats
import scipy.sparse
import pyhmmer
import tqdm
from typing import Iterator, List
from domainator.utils import get_file_type, parse_seqfiles, make_pool, pyhmmer_decode
from domainator import __version__, RawAndDefaultsFormatter
from domainator.data_matrix import DataMatrix, StreamingMstKnnAccumulator
from domainator.hmmer_search import compare_hmmer
from domainator.output_guardrails import add_max_output_gb_argument, enforce_matrix_output_limit, format_bytes, max_output_gb_to_bytes, OutputSizeLimitExceeded
from domainator.transform_matrix import MODES, _mst_knn_arg, _knn_arg, LOG10_2

# Modes whose post-transform values are computable per edge (no global matrix maxes),
# and therefore compatible with streaming --mst_knn pruning.
MST_KNN_STREAMABLE_MODES = {"score", "bool", "efi_score"}
import psutil

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
# --params keys that mean the caller has already bounded diamond's scratch themselves.
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
    lines.append("  To bound it, pass block/index-chunk options through --params, e.g.:")
    lines.append("      --params '\"-b\":0.07, \"-c\":4'")
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
    "--params '\"-b\":0.01, \"-c\":4'. Peak scratch scales roughly with -b.",
    "Put the scratch on a larger filesystem: set TMPDIR, or "
    "--params '\"--tmpdir\":\"/path/with/space\"'.",
    "Consider fewer input sequences (for example dereplicating first); the seed-hit "
    "volume grows quadratically with the size of a redundant input.",
)
_MEMORY_MITIGATIONS = (
    "Shrink diamond's query blocks: --params '\"-b\":0.01, \"-c\":4'.",
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
    (scratch disk exhaustion, memory exhaustion, a rejected --params flag, unreadable
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
        lines.append("  This happens before diamond emits any results, so --mst_knn/--knn")
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
            rendered = ", ".join(f'"{flag}":{value!r}' for flag, value in diamond_params.items())
            lines.append(f"  Options passed via --params: {rendered}")
            lines.append("  Check these against `diamond blastp --help`; --params is forwarded as-is.")
        else:
            lines.append("  No --params were supplied, so this may be a diamond version mismatch.")
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


# diamond options that seq_dist sets itself. Overriding these through --params would
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
    "--min-score": "--lb",
    "--verbose": "seq_dist (required for progress reporting)",
    "--fast": "--algorithm",
    "--mid-sensitive": "--algorithm",
    "--sensitive": "--algorithm",
    "--more-sensitive": "--algorithm",
    "--very-sensitive": "--algorithm",
    "--ultra-sensitive": "--algorithm",
}

# diamond options seq_dist supplies a default for, but which --params may replace.
def _default_overridable_diamond_params(tmpdir):
    return OrderedDict((("--algo", "0"), ("--tmpdir", tmpdir)))


def parse_diamond_params(params_str):
    """Parse a --params string into an ordered dict of {flag: value}.

    The string is a json dict body without the enclosing braces, matching the
    --params syntax of deduplicate_genbank.py, for example: '"-c":4, "-b":0.4'.
    A null value means the flag takes no argument, for example: '"--target-indexed":null'.
    """
    if not params_str:
        return OrderedDict()
    try:
        parsed = json.loads("{" + params_str + "}", object_pairs_hook=OrderedDict)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Could not parse --params as a json dict: {exc}\n"
            "Expected a json dict body without the outer braces, such as: "
            "--params '\"-c\":4, \"-b\":0.4'"
        ) from exc
    if not isinstance(parsed, dict):
        raise ValueError("--params must be a json dict of diamond options, such as: --params '\"-c\":4'")
    for flag in parsed:
        if not flag.startswith("-"):
            raise ValueError(
                f"--params keys must be diamond option flags starting with '-', but got: {flag!r}"
            )
        if flag in RESERVED_DIAMOND_PARAMS:
            raise ValueError(
                f"--params may not set '{flag}', because seq_dist controls it via "
                f"{RESERVED_DIAMOND_PARAMS[flag]}."
            )
    return parsed


def add_params_to_args_list(arg_list, params):
    """Append {flag: value} pairs to a command line. A null value emits a bare flag."""
    for flag, value in params.items():
        arg_list.append(flag)
        if value is not None:
            arg_list.append(str(value))


def diamond(input_fasta, reference_fasta, max_target_seqs, threads, tmpdir, mode, max_hsps=1, min_score: float = 0.0, query_order=None, progress: bool = False,
            diamond_params=None) -> Iterator[CmpResult]:
    """
        runs a search using diamond, returns a list of tuples of (inputid, referenceid, score)
    """
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
    mkdb_command = ["diamond", "makedb", "--in", reference_fasta, "-d", dbpath]
    try:
        mkdb_result = subprocess.run(mkdb_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                encoding='utf-8')
    except FileNotFoundError:
        print("Could not run 'diamond'. It is not on PATH.\n"
              "  Install it (conda install -c bioconda diamond) or activate the environment\n"
              "  that provides it, then re-run.", file=sys.stderr)
        exit(1)

    if mkdb_result.returncode != 0:
        print(build_diamond_failure_message(
            "makedb", mkdb_result.returncode,
            (mkdb_result.stderr or "") + (mkdb_result.stdout or ""),
            reference_fasta, tmpdir=tmpdir, command=mkdb_command), file=sys.stderr)
        exit(1)

    dmnd_command = [ "stdbuf", "-oL", "-eL",
        "diamond", "blastp", "-q", input_fasta, "-d", dbpath, "-p", str(threads),
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

    # Defaults seq_dist supplies but --params may replace. --tmpdir defaults to the
    # managed temporary directory so that diamond's (large, unlinked-while-open)
    # search scratch lands next to the database instead of in the working directory.
    overridable_params = _default_overridable_diamond_params(tmpdir)
    if diamond_params:
        overridable_params.update(diamond_params)
    add_params_to_args_list(dmnd_command, overridable_params)

    # Echo the exact command so a run can be reproduced or retuned straight from a log.
    print("Running diamond:\n  " + " ".join(dmnd_command), file=sys.stderr, flush=True)

    try:
        dmnd_process = subprocess.Popen(
            dmnd_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding='utf-8',
            bufsize=1,
        )
    except FileNotFoundError as exc:
        missing = "stdbuf" if "stdbuf" in str(exc) else "diamond"
        print(f"Could not run '{missing}'. It is not on PATH.\n"
              f"  Command: {' '.join(dmnd_command)}", file=sys.stderr)
        exit(1)

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
        print(build_diamond_failure_message(
            "blastp", returncode, "".join(stderr_lines), input_fasta,
            tmpdir=tmpdir, command=dmnd_command, diamond_params=diamond_params), file=sys.stderr)
        exit(1)

def make_seq_name_to_idx_dict(input_path):
    name_to_idx = OrderedDict()
    idx_to_name = list()
    idx_to_len = list()
    idx = 0
    for rec in parse_seqfiles([input_path]):
        id = rec.id
        if id in name_to_idx:  
            warnings.warn(f"Warning: duplicate sequence name found in input {input_path}: {rec.id}")
        else:
            name_to_idx[id] = idx
            idx_to_len.append(len(rec.seq))
            idx_to_name.append(id)
            idx += 1
    return name_to_idx, idx_to_name, idx_to_len

def make_hmm_name_to_idx_dict_from_path(hmm_path):
    name_to_idx = OrderedDict()
    idx_to_name = list()
    idx_to_len = list()

    idx = 0
    for rec in pyhmmer.plan7.HMMFile(hmm_path):
        name = pyhmmer_decode(rec.name)
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
        name = pyhmmer_decode(rec.name)
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
            name = pyhmmer_decode(hit.name)
            domain = pyhmmer_decode(hit.best_domain.alignment.hmm_name)
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
    
    def __call__(self, input_profile):
        out_heap = []
        for target_profile in pyhmmer.plan7.HMMFile(self.hmmer_targets):
            
            score, _traceback, _max_index, _match_scores = compare_hmmer(input_profile, target_profile)
            score = round(score,2)
            if score <= self.min_score:
                continue
            result = CmpResult(score, pyhmmer_decode(input_profile.name), pyhmmer_decode(target_profile.name))
            if (self.k is None) or (len(out_heap) < self.k):
                heapq.heappush(out_heap, result)
            else:
                heapq.heappushpop(out_heap, result)
        return out_heap

def run_hmmer_compare(hmmer_queries, hmmer_targets, k, threads, min_score: float = 0.0, progress: bool = False, query_count=None):
    
    worker = _run_hmmer_compare_worker(hmmer_targets, k, min_score=min_score)
    with make_pool(processes=threads) as pool:
        for hits in _progress(pool.imap_unordered(worker, pyhmmer.plan7.HMMFile(hmmer_queries)), progress, total=query_count, desc="hmmer_compare queries"):
            for hit in sorted(hits, key=lambda item: (-item.score, item.reference)):
                yield hit

def convert_to_fasta(input_path, output_path):
    with open(output_path,"w") as outh:
        for rec in parse_seqfiles([input_path]):
            SeqIO.write(rec,outh,"fasta")


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


def seq_dist(input_path, input_type, reference_path, reference_type, k, algorithm, mode, threads, dense, dense_text, sparse, lb, max_output_bytes=None, progress: bool = False, mst_knn=None, knn=None, diamond_params=None):
    """

    """
    if mst_knn is not None and knn is not None:
        raise ValueError("--mst_knn already includes the kNN edges; pass either --mst_knn or --knn, not both.")
    # One streaming accumulator serves both flags; only the spanning tree differs.
    sparsify_k = mst_knn if mst_knn is not None else knn
    sparsify_include_mst = mst_knn is not None
    if sparsify_k is not None:
        option = "--mst_knn" if sparsify_include_mst else "--knn"
        if mode not in MST_KNN_STREAMABLE_MODES:
            raise ValueError(
                f"{option} is only supported with --mode {sorted(MST_KNN_STREAMABLE_MODES)} "
                f"(modes that normalize against global matrix maxima cannot be streamed). "
                f"For '{mode}', run 'seq_dist --mode score --sparse out.hdf5' then "
                f"'transform_matrix -i out.hdf5 --mode {mode} {option} {sparsify_k} --sparse pruned.hdf5'."
            )
        if reference_path != input_path:
            raise ValueError(
                f"{option} requires a square, symmetric comparison (reference must be the same file as input)."
            )
    if sparse is not None and mode == "score_dist":
        warnings.warn(f"Sparse output requested for a score_dist matrix. Distance matrices are usually not sparse. This may use more memory and have a larger output file size than expected.")

    if input_type == "hmm" and reference_type == "hmm":
        if algorithm != "hmmer_compare":
            algorithm = "hmmer_compare"
            warnings.warn(f"Input and reference are both hmm files, but algorithm is not hmmer_compare. Changing algorithm to hmmer_compare.")
    search_type = None
    with tempfile.TemporaryDirectory() as tmpdir:
        # each of these needs to set query_names, db_names, and results_table
        
        if algorithm.startswith("diamond"):
            parts = algorithm.split("_")

            # convert input and db to fasta if not fasta
            dmnd_input = input_path
            if input_type != "fasta":
                dmnd_input = tmpdir + "/input.fasta"
                convert_to_fasta(input_path, dmnd_input) #TODO: better error message when supplied with hmmer file
            query_name_to_idx, query_idx_to_name, query_idx_to_len = make_seq_name_to_idx_dict(dmnd_input) #name to index
            db_name_to_idx, db_idx_to_name, db_idx_to_len = query_name_to_idx, query_idx_to_name, query_idx_to_len
            if input_path == reference_path:
                dmnd_ref = dmnd_input
            else:
                dmnd_ref = reference_path
                if reference_type != "fasta":
                    dmnd_ref = tmpdir + "/ref.fasta"
                    convert_to_fasta(reference_path, dmnd_ref)
                db_name_to_idx, db_idx_to_name, db_idx_to_len = make_seq_name_to_idx_dict(dmnd_ref) #name to index
            preflight = build_diamond_preflight_warning(
                parts[1], len(query_name_to_idx), int(sum(query_idx_to_len)),
                square=(input_path == reference_path),
                tmpdir=tmpdir, diamond_params=diamond_params,
            )
            if preflight is not None:
                print(preflight, file=sys.stderr, flush=True)

            # call diamond
            
            results_iter = diamond(dmnd_input, dmnd_ref, k, threads, tmpdir, parts[1], min_score=lb, query_order=query_idx_to_name, progress=progress, diamond_params=diamond_params)
            score_progress = False
        
        elif algorithm == "hmmer":
            hmmer_seqfile = input_path
            if input_type != "fasta": #might be possible to just go from 
                hmmer_seqfile = tmpdir + "/input.fasta"
                convert_to_fasta(input_path, hmmer_seqfile)

            with pyhmmer.easel.SequenceFile(hmmer_seqfile, digital=True) as seq_file:
                hmmer_seqs = list(seq_file)
            query_name_to_idx, query_idx_to_name, query_idx_to_len = make_seq_name_to_idx_dict(hmmer_seqfile) #name to index


            if reference_type == "hmm":
                db_data = list(pyhmmer.plan7.HMMFile(reference_path))
                db_name_to_idx, db_idx_to_name, db_idx_to_len = make_hmm_name_to_idx_dict(db_data)
                search_type = "hmmsearch"
            else: # a sequence type
                search_type = "phmmer"
                hmmer_dbfile = reference_path
                if reference_type != "fasta":  
                    hmmer_dbfile = tmpdir + "/ref.fasta"
                    convert_to_fasta(reference_path, hmmer_dbfile)
                with pyhmmer.easel.SequenceFile(hmmer_dbfile, digital=True) as seq_file:
                    db_data = list(seq_file)
                db_name_to_idx, db_idx_to_name, db_idx_to_len = make_seq_name_to_idx_dict(hmmer_dbfile)
                
            #TODO: do easel and biopython split names and descriptions in the same way?
                       
            results_iter = run_hmmsearch(hmmer_seqs, db_data, k, threads, search_type, min_score=lb, progress=progress, query_count=len(db_data))
            score_progress = False
        
        
        elif algorithm == "hmmer_compare":
            query_name_to_idx, query_idx_to_name, query_idx_to_len = make_hmm_name_to_idx_dict_from_path(input_path)
            db_name_to_idx, db_idx_to_name, db_idx_to_len = query_name_to_idx, query_idx_to_name, query_idx_to_len
            if input_path != reference_path:
                db_name_to_idx, db_idx_to_name, db_idx_to_len = make_hmm_name_to_idx_dict_from_path(reference_path)
            results_iter = run_hmmer_compare(input_path, reference_path, k, threads, min_score=lb, progress=progress, query_count=len(query_name_to_idx))
            score_progress = False

        else: 
            raise RuntimeError(f"Algorithm not recognized: {algorithm}")

        

        if sparsify_k is not None:
            # Stream edges straight into the kNN (optionally plus MST) accumulator with
            # bounded memory; the full all-vs-all matrix is never materialized. mode value
            # is computed per edge (no global maxima needed for score/bool/efi_score). The
            # raw-score --lb gate runs first (matching the help text), then the
            # accumulator prunes on values > 0.
            log_lengths = np.log10(np.asarray(query_idx_to_len, dtype=np.float64)) if mode == "efi_score" else None
            accumulator = StreamingMstKnnAccumulator(len(query_name_to_idx), sparsify_k, lower_bound=0,
                                                     include_mst=sparsify_include_mst)
            for (score, query_id, subject_id) in _progress(results_iter, score_progress, desc="Scoring pairs"):
                if score <= lb:
                    continue
                query_idx = query_name_to_idx[query_id]
                subject_idx = db_name_to_idx[subject_id]
                if mode == "bool":
                    value = 1.0
                elif mode == "efi_score":
                    value = score * LOG10_2 - log_lengths[query_idx] - log_lengths[subject_idx]
                    if value <= 0:
                        continue
                else:  # score
                    value = score
                accumulator.add_edge(query_idx, subject_idx, value)
            matrix = accumulator.to_csr()
        else:
            result_builder = _make_result_builder(
                algorithm,
                query_name_to_idx,
                db_name_to_idx,
                search_type=search_type,
                query_order=query_idx_to_name,
            )

            for (score, query_id, subject_id) in _progress(results_iter, score_progress, desc="Scoring pairs"):
                if score <= lb:
                    continue
                matrix_score = 1 if mode == "bool" else score
                result_builder.add_result(query_id, subject_id, matrix_score)

            matrix = result_builder.build()

    query_idx_to_len = np.array(query_idx_to_len)
    db_idx_to_len = np.array(db_idx_to_len)
    if progress:
        print("Transforming matrix...", flush=True, file=sys.stderr)
    if sparsify_k is None and MODES[mode] is not None:
        matrix = MODES[mode](matrix, query_idx_to_len, db_idx_to_len)


    try:
        if dense is not None:
            enforce_matrix_output_limit(
                output_type="dense",
                matrix=matrix,
                row_names=query_idx_to_name,
                col_names=db_idx_to_name,
                row_lengths=query_idx_to_len,
                col_lengths=db_idx_to_len,
                data_type=mode,
                max_output_bytes=max_output_bytes,
                output_path=dense,
                mitigation_options=["--sparse", "-k", "--lb"],
            )
            DataMatrix.write_dense(matrix, dense, query_idx_to_name, db_idx_to_name, query_idx_to_len, db_idx_to_len, mode)
        if dense_text is not None:
            enforce_matrix_output_limit(
                output_type="dense_text",
                matrix=matrix,
                row_names=query_idx_to_name,
                col_names=db_idx_to_name,
                row_lengths=query_idx_to_len,
                col_lengths=db_idx_to_len,
                data_type=mode,
                max_output_bytes=max_output_bytes,
                output_path=dense_text,
                mitigation_options=["--sparse", "-k", "--lb"],
            )
            DataMatrix.write_dense_text(matrix, dense_text, query_idx_to_name, db_idx_to_name, query_idx_to_len, db_idx_to_len, mode)
        if sparse:
            enforce_matrix_output_limit(
                output_type="sparse",
                matrix=matrix,
                row_names=query_idx_to_name,
                col_names=db_idx_to_name,
                row_lengths=query_idx_to_len,
                col_lengths=db_idx_to_len,
                data_type=mode,
                max_output_bytes=max_output_bytes,
                output_path=sparse,
                mitigation_options=["-k", "--lb"],
            )
            DataMatrix.write_sparse(matrix, sparse, query_idx_to_name, db_idx_to_name, query_idx_to_len, db_idx_to_len, mode)
    except OutputSizeLimitExceeded as exc:
        raise SystemExit(str(exc)) from None

def main(argv):
    #TODO: support for nucleotide sequences? (see compare_contigs.py), can be implemented easily with nhmmer, without adding new dependencies.
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)
    parser.add_argument('-i', '--input', type=str, required=True,
                        help="Input query file. Queries in the input file will be compared to the reference sequences/profiles. This can be a fasta, genbank, or hmm file.")

    parser.add_argument('--input_type', type=str, required=False, default=None, choices={"fasta", "genbank", "hmm"},
                        help="File type for the sequence file.")

    parser.add_argument('-r', "--reference", type=str, required=False,
                        help="Reference file. If not provided, the input file will be used as the reference, and the output will be a pairwise similarity matrix between the input sequences/profiles. This can be a fasta, genbank, or hmm file.")

    parser.add_argument('--reference_type', type=str, required=False, default=None, choices={"fasta", "genbank", "hmm"},
                        help="File type for the sequence file. Ignored for some kinds of algorithms, for example, hmmer expects a .hmm reference file.")

    parser.add_argument("--dense", type=str, default=None, help="Write a dense distance matrix hdf5 file to this path.")
    
    parser.add_argument("--dense_text", type=str, default=None, help="Write a dense distance matrix tsv file to this path.")
    
    parser.add_argument("--sparse", type=str, default=None, help="Write a sparse distance matrix hdf5 file to this path.")

    parser.add_argument('--lb', default=0, type=float, required=False,
                        help="Round any scores lower than this down to zero (this applies to raw scores not processed/normalized).")


    parser.add_argument('-k', type=int, required=False, default=None,
                        help="Include at most this many non-zero entries in the matrix for each input sequence. Default: Include all hits.")

    parser.add_argument('--algorithm', type=str, required=False, default="diamond_us", choices={"diamond_f", "diamond_d", "diamond_mids", "diamond_s", "diamond_mors", "diamond_vs", "diamond_us", "hmmer", "hmmer_compare"},
                        help="Which distance metric to use, diamond_*; f: fast, d: default, mids: mid-sensitive, s: sensitive, mors: more-sensitive, vs: very-sensitive, us: ultra-sensitive. default: diamond_us (ultra sensitive)")

    parser.add_argument("--params", default=None,
                        help='String of extra parameters to pass to diamond blastp, as a json dict body in single quotes. '
                             'Example: \'"-c":4, "-b":0.4\'. Use a null value for flags that take no argument, such as \'"--target-indexed":null\'. '
                             "Options that seq_dist sets itself (query, database, --outfmt, -p, --max-target-seqs, --max-hsps, --min-score, --verbose, and the sensitivity flags) "
                             "are rejected; use the corresponding seq_dist option instead. "
                             "seq_dist defaults --algo to 0 and --tmpdir to its own managed temporary directory; both can be replaced here. "
                             "Useful for large all-vs-all runs, where diamond's seed-search scratch can exceed available disk: "
                             "a smaller --block-size (-b) and/or more --index-chunks (-c) bound its peak scratch usage.")

    sparsify_group = parser.add_mutually_exclusive_group(required=False)
    sparsify_group.add_argument('--mst_knn', type=_mst_knn_arg, required=False, default=None,
                        help="Prune the output to the maximum spanning tree plus OR-symmetric k-nearest-neighbor edges (integer >= 0, where 0 keeps only the maximum spanning tree), computed as a streaming operation to keep memory and output size small. Requires the reference to be the same file as the input (a square symmetric comparison) and --mode in {score, bool, efi_score}. Best paired with --sparse.")
    sparsify_group.add_argument('--knn', type=_knn_arg, required=False, default=None,
                        help="Prune the output to OR-symmetric k-nearest-neighbor edges only (integer >= 1), computed as a streaming operation. Unlike --mst_knn this does not preserve the connected components of the full comparison. Same requirements as --mst_knn.")

    parser.add_argument("--mode", type=str, required=False, default="score", choices=set(MODES.keys()),
                        help="what kind of values should be in the matrix. score: raw score, bool: 1 if a hit otherwise 0, score_dist: 1 - (score / min(row_max, col_max)), norm_score: score/min(row_max, col_max), efi_score: -log10[2^(-score) * (input_seq_length * reference_seq_length)], efi_score_dist: 1 - (efi_score / min(row_max, col_max)). Default: score")

    parser.add_argument('--cpu', type=int, default=0, required=False,
                        help="how many cpu threads to use. Default: use all available cores.")

    parser.add_argument('--progress', action='store_true',
                        help="Show a progress bar for long-running steps.")

    add_max_output_gb_argument(parser)

    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)

    ### make output directory

    input_type = params.input_type
    if params.input_type is None:
        input_type = get_file_type(params.input)

    if params.reference is None:
        reference_path = params.input
        reference_type = input_type
    else:
        reference_path = params.reference
        reference_type = params.reference_type
        if params.reference_type is None:
            reference_type = get_file_type(reference_path)

    if params.dense is not None and get_file_type(params.dense) != "hdf5":
        raise ValueError("Please use an hdf5 related extension for the --dense output, such as .h5, .hdf5, or .hdf.")
    if params.sparse is not None and get_file_type(params.sparse) != "hdf5":
        raise ValueError("Please use an hdf5 related extension for the --sparse output, such as .h5, .hdf5, or .hdf.")
    dense = params.dense
    dense_text = params.dense_text
    sparse = params.sparse

    if params.cpu <= 0:
        cpus = psutil.cpu_count(logical=False)
    else:
        cpus = params.cpu

    if ((dense is None) and (dense_text is None) and (sparse is None)):
        raise ValueError("No output specified! Please specify at least one of: dense, dense_text, sparse")

    if sparse is not None and params.mode == "score_dist":
        raise ValueError("Sparse distance matrices not implemented.")

    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)

    diamond_params = parse_diamond_params(params.params)
    if diamond_params and not params.algorithm.startswith("diamond"):
        raise ValueError(f"--params only applies to the diamond_* algorithms, but --algorithm is '{params.algorithm}'.")

    seq_dist(params.input, input_type, reference_path, reference_type, params.k, params.algorithm, params.mode, cpus, dense, dense_text, sparse, params.lb, max_output_bytes=max_output_bytes, progress=params.progress, mst_knn=params.mst_knn, knn=params.knn, diamond_params=diamond_params)

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])

