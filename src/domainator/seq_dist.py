"""Calculates similarity scores between protein sequences or HMM profiles

Typically the input and reference will be the same sequence file, creating a pairwise similarity matrix. 
An hmm file can also be used as a reference to make a table of profile scores for each input peptide.
An hmm file can also be used as both input and reference to compare profile similarity.

"""
import warnings
warnings.filterwarnings("ignore", module='numpy')

from jsonargparse import ArgumentParser, ActionConfigFile
import heapq
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
from domainator.data_matrix import DataMatrix
from domainator.hmmer_search import compare_hmmer
from domainator.output_guardrails import add_max_output_gb_argument, enforce_matrix_output_limit, max_output_gb_to_bytes, OutputSizeLimitExceeded
from domainator.transform_matrix import MODES
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

def diamond(input_fasta, reference_fasta, max_target_seqs, threads, tmpdir, mode, max_hsps=1, min_score: float = 0.0, query_order=None, progress: bool = False) -> Iterator[CmpResult]:
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
    mkdb_result = subprocess.run(["diamond", "makedb", "--in", reference_fasta, "-d", dbpath], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            encoding='utf-8')

    if mkdb_result.returncode != 0:
        print(
            f'Error in diamond makedb: {reference_fasta}\n{mkdb_result.stdout}\n{mkdb_result.stderr}', file=sys.stderr)
        exit(1)

    dmnd_command = [ "stdbuf", "-oL", "-eL",
        "diamond", "blastp", "-q", input_fasta, "-d", dbpath, "-p", str(threads),
        "--outfmt", "6", "qseqid", "sseqid", "bitscore", "--algo", "0",
        "--max-hsps", str(max_hsps), "--verbose"
    ]
    if max_target_seqs is not None and max_target_seqs > 0:
        dmnd_command.extend(["--max-target-seqs", str(max_target_seqs)])
    else:
        dmnd_command.extend(["--max-target-seqs", str(0)])
    if min_score > 0:
        dmnd_command.extend(["--min-score", str(min_score)])
    dmnd_command.extend(mode_to_arg[mode])

    dmnd_process = subprocess.Popen(
        dmnd_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='utf-8',
        bufsize=1,
    )

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
        print(
            f'Error in diamond search: {input_fasta}\n{"".join(stderr_lines)}', file=sys.stderr)
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


def seq_dist(input_path, input_type, reference_path, reference_type, k, algorithm, mode, threads, dense, dense_text, sparse, lb, max_output_bytes=None, progress: bool = False):
    """

    """
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
            # call diamond
            
            results_iter = diamond(dmnd_input, dmnd_ref, k, threads, tmpdir, parts[1], min_score=lb, query_order=query_idx_to_name, progress=progress)
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
    if MODES[mode] is not None:
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

    seq_dist(params.input, input_type, reference_path, reference_type, params.k, params.algorithm, params.mode, cpus, dense, dense_text, sparse, params.lb, max_output_bytes=max_output_bytes, progress=params.progress)

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])

