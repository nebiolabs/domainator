"""Calculates similarity scores between protein sequences. 

Typically the input and reference will be the same sequence file, creating a pairwise similarity matrix. 
An hmm file can also be used as a reference to make a table of profile scores for each input peptide.

"""

import argparse
import heapq
import sys
import subprocess
import tempfile
import warnings
from domainator.Bio import SeqIO
from collections import OrderedDict,namedtuple, defaultdict
import numpy as np
import scipy.stats
import scipy.sparse
import pyhmmer
from typing import List
from domainator.utils import get_file_type, parse_seqfiles
from domainator import __version__, RawAndDefaultsFormatter
from domainator.data_matrix import DataMatrix
from domainator.hmmer_search import compare_hmmer
from domainator.transform_matrix import MODES
from multiprocessing import Pool
import psutil

#TODO: maybe merge with compare_contigs.py? Domainatorify more than it is.

#TODO: ESM embeddings! (we could do direct PCA on the embeddings, or convert the embeddings to pairwise distances then do PCA or UMAP)
#TODO: extract scores from domainator annotations
#TODO: BigScape/CORASON-style scoring for domainator/hmm scores, also support anti-smash gb files (See compare_contigs)
#TODO: rank is kind of weird, and also not compatible with sparse data I think it should be 0 as the worst and num_rows as the best.
#TODO: upgrade to scipy 1.8
#TODO: it's more like seq_similarity than seq_dist
#TODO: add support for comparing hmm profiles, and/or MSAs
#TODO: maybe negative log evalue output?

#TODO: with sparse output, I think we are destroying singletons, can we account for that somehow, maybe by having a 1-column line?

CmpResult = namedtuple("CmpResult", ["score", "input", "reference"])

def diamond(input_fasta, reference_fasta, max_target_seqs, threads, tmpdir, mode, max_hsps=1) -> List[CmpResult]:
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
    out_path=tmpdir+"/"+"dmnd_out.txt"
    mkdb_result = subprocess.run(["diamond", "makedb", "--in", reference_fasta, "-d", dbpath], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            encoding='utf-8')

    if mkdb_result.returncode != 0:
        print(
            f'Error in diamond makedb: {reference_fasta}\n{mkdb_result.stdout}\n{mkdb_result.stderr}', file=sys.stderr)
        exit(1)

    dmnd_out = subprocess.run(["diamond", "blastp", "-q", input_fasta, "-d", dbpath, "-p", str(threads), "--outfmt", "6", "qseqid", "sseqid", "bitscore", 
                            "--algo", "0", "--max-hsps", str(max_hsps), "-k", str(max_target_seqs), "-o", out_path] + mode_to_arg[mode], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')

    if dmnd_out.returncode != 0:
        print(
            f'Error in diamond search: {input_fasta}\n{dmnd_out.stdout}\n{dmnd_out.stderr}', file=sys.stderr)
        exit(1)
    
    #out = list()
    with open(out_path) as dmnd_results:
        for line in dmnd_results:
            line = line.strip()
            parts = line.split("\t")
            if len(parts) == 3:
                yield CmpResult(float(parts[2]), parts[0], parts[1])

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
        name = rec.name.decode()
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
        name = rec.name.decode()
        if name in name_to_idx:  
            warnings.warn(f"Warning: duplicate hmm name found: {name}")
        else:
            name_to_idx[name] = idx
            idx_to_name.append(name)
            idx_to_len.append(rec.M)
            idx += 1
    
    return name_to_idx, idx_to_name, idx_to_len

def run_hmmsearch(hmmer_seqs, hmmer_models, k, threads, search_type):
    if search_type == "hmmsearch":
        hits = pyhmmer.hmmer.hmmsearch(hmmer_models, hmmer_seqs, cpus=threads)
    elif search_type == "phmmer":
        hits = pyhmmer.hmmer.phmmer(hmmer_models, hmmer_seqs, cpus=threads)
    else:
        raise ValueError(f"search_type {search_type} not recognized")

    hits_list = list()
    for top_hits in hits:
        for hit in top_hits:
            name = hit.name.decode()
            domain = hit.best_domain.alignment.hmm_name.decode()
            score = hit.score
            result = CmpResult(round(score,2), name, domain)
            if k is None:
                yield result
            else:
                hits_list.append(result)

    #TODO: make this faster/more memory efficient by keeping heaps of size k for each query
    del hits
    if k is not None: #then we need to filter the hits somehow
        seen = defaultdict(lambda x: 0)
        
        hits_list.sort(key=lambda x: x.score)
        for result in hits_list:
            if seen[result.input] < k:
                yield result
                seen[result.input] += 1


class _run_hmmer_compare_worker():
    def __init__(self, hmmer_targets, k=None):
        self.k = k
        self.hmmer_targets = hmmer_targets
    
    def __call__(self, input_profile):
        out_heap = []
        for target_profile in pyhmmer.plan7.HMMFile(self.hmmer_targets):
            
            score, _traceback, _max_index, _match_scores = compare_hmmer(input_profile, target_profile)
            result = CmpResult(round(score,2), input_profile.name.decode(), target_profile.name.decode())
            if (self.k is None) or (len(out_heap) < self.k):
                heapq.heappush(out_heap, result)
            else:
                heapq.heappushpop(out_heap, result)
        return out_heap

def run_hmmer_compare(hmmer_queries, hmmer_targets, k, threads):
    
    worker = _run_hmmer_compare_worker(hmmer_targets,k)
    with Pool(processes=threads) as pool:
        for hits in pool.imap_unordered(worker, pyhmmer.plan7.HMMFile(hmmer_queries)):
            for hit in hits:
                yield hit

def convert_to_fasta(input_path, output_path):
    with open(output_path,"w") as outh:
        for rec in parse_seqfiles([input_path]):
            SeqIO.write(rec,outh,"fasta")


def seq_dist(input_path, input_type, reference_path, reference_type, k, algorithm, mode, threads, dense, dense_text, sparse, lb):
    """

    """
    if sparse is not None and mode == "score_dist":
        warnings.warn(f"Sparse output requested for a score_dist matrix. Distance matrices are usually not sparse. This may use more memory and have a larger output file size than expected.")

    if input_type == "hmm" and reference_type == "hmm":
        if algorithm != "hmmer_compare":
            algorithm = "hmmer_compare"
            warnings.warn(f"Input and reference are both hmm files, but algorithm is not hmmer_compare. Changing algorithm to hmmer_compare.")
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
            
            results_iter = diamond(dmnd_input, dmnd_ref, k, threads, tmpdir, parts[1])
        
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
                       
            results_iter = run_hmmsearch(hmmer_seqs, db_data, k, threads, search_type)
        
        
        elif algorithm == "hmmer_compare":
            query_name_to_idx, query_idx_to_name, query_idx_to_len = make_hmm_name_to_idx_dict_from_path(input_path)
            db_name_to_idx, db_idx_to_name, db_idx_to_len = query_name_to_idx, query_idx_to_name, query_idx_to_len
            if input_path != reference_path:
                db_name_to_idx, db_idx_to_name, db_idx_to_len = make_hmm_name_to_idx_dict_from_path(reference_path)
            results_iter = run_hmmer_compare(input_path, reference_path, k, threads)

        else: 
            raise RuntimeError(f"Algorithm not recognized: {algorithm}")

        

        matrix = scipy.sparse.dok_array((len(query_name_to_idx), len(db_name_to_idx)),dtype=np.float64)

        for (score, query_id, subject_id) in results_iter:
            if mode == "bool": # in bool mode any entry gets converted to a 1
                score = 1
            if score > lb:
                if score > matrix[query_name_to_idx[query_id], db_name_to_idx[subject_id]]: # if there are multiple pairwise hits, keep only the largest.
                    matrix[query_name_to_idx[query_id], db_name_to_idx[subject_id]] = score
        
    query_idx_to_len = np.array(query_idx_to_len)
    db_idx_to_len = np.array(db_idx_to_len)
    if MODES[mode] is not None:
        matrix = MODES[mode](matrix, query_idx_to_len, db_idx_to_len)
    

    if dense is not None:
        DataMatrix.write_dense(matrix, dense, query_idx_to_name, db_idx_to_name, query_idx_to_len, db_idx_to_len, mode)
    if dense_text is not None:
        DataMatrix.write_dense_text(matrix, dense_text, query_idx_to_name, db_idx_to_name, query_idx_to_len, db_idx_to_len, mode)
    if sparse:
        DataMatrix.write_sparse(matrix, sparse, query_idx_to_name, db_idx_to_name, query_idx_to_len, db_idx_to_len, mode)

def main(argv):
    #TODO: support for nucleotide sequences? (see compare_contigs.py)
    parser = argparse.ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)
    parser.add_argument('-i', '--input', type=str, required=True,
                        help="Input query file.")

    parser.add_argument('--input_type', type=str, required=False, default=None, choices={"fasta", "genbank", "hmm"},
                        help="File type for the sequence file")

    parser.add_argument('-r', "--reference", type=str, required=True,
                        help="Reference file.")

    parser.add_argument('--reference_type', type=str, required=False, default=None, choices={"fasta", "genbank", "hmm"},
                        help="File type for the sequence file. Ignored for some kinds of algorithms, for example, hmmer expects a .hmm reference file.")

    # TODO: is it weird to not have a -o parameter?
    # parser.add_argument('-o', '--output', type=str, required=True,
    #                     help="file to write the distance matrix to as a sparse hdf5.")
    parser.add_argument("--dense", type=str, default=None, help="Write a dense distance matrix hdf5 file to this path.")
    
    parser.add_argument("--dense_text", type=str, default=None, help="Write a dense distance matrix tsv file to this path.")
    
    parser.add_argument("--sparse", type=str, default=None, help="Write a sparse distance matrix hdf5 file to this path.")

    # parser.add_argument("--hist", type=str, default=None, help="Write a histogram figure of the data to this path.") #TODO: implement this, maybe as a plotly html?, how should 0s (1s in dist mode) be handled, maybe write a count but don't include in hist, maybe there is other information that would be helpful in a report, like number of singletons?
    
    parser.add_argument('--lb', default=0, type=float, required=False,
                        help="Round any scores lower than this down to zero (this applies to raw scores not processed/normalized).")


    parser.add_argument('-k', type=int, required=False, default=None,
                        help="Include at most this many non-zero entries in the matrix for each input sequence. Default: Include all hits.")

    parser.add_argument('--algorithm', type=str, required=False, default="diamond_us", choices={"diamond_f", "diamond_d", "diamond_mids", "diamond_s", "diamond_mors", "diamond_vs", "diamond_us", "hmmer", "hmmer_compare"},
                        help="Which distance metric to use, diamond_*; f: fast, d: default, mids: mid-sensitive, s: sensitive, mors: more-sensitive, vs: very-sensitive, us: ultra-sensitive. default: diamond_us (ultra sensitive)")

    #TODO: random subset of references

    #TODO: what should the default be for rank, max_rank+1?
    #TODO: add rank back rank: hit scores converted to ranks,
    #TODO: maybe better to use raw score as default and then cluster on cosine distance? Would that be similar to using a normalized score matrix?
    #TODO: would efi_score_dist be useful? 
    parser.add_argument("--mode", type=str, required=False, default="score", choices=set(MODES.keys()), 
                        help="what kind of values should be in the matrix. score: raw score, bool: 1 if a hit otherwise 0, score_dist: 1 - (score / min(row_max, col_max)), norm_score: score/min(row_max, col_max), efi_score: -log10[2^(-score) * (input_seq_length * reference_seq_length)], efi_score_dist: 1 - (efi_score / min(row_max, col_max)). Default: score")

    parser.add_argument('--cpu', type=int, default=0, required=False,
                        help="how many cpu threads to use. Default: use all available cores.")

    params = parser.parse_args(argv)

    ### make output directory

    input_type = params.input_type
    if params.input_type is None:
        input_type = get_file_type(params.input)

    reference_type = params.reference_type
    if params.reference_type is None:
        reference_type = get_file_type(params.reference)

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

    seq_dist(params.input, input_type, params.reference, reference_type, params.k, params.algorithm, params.mode, cpus, dense, dense_text, sparse, params.lb)

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])

