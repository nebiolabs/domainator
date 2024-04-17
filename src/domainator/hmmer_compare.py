"""Aligns and calculates alignment scores between hmmer3 profiles 

Kind of like hhsearch except much slower and for hmmer3 profiles instead of hhsuite profiles.

Adapted from pseudocode in:
Steinegger, Martin, Markus Meier, Milot Mirdita, Harald Vöhringer, Stephan J. Haunsberger, and Johannes Söding. “HH-Suite3 for Fast Remote Homology Detection and Deep Protein Annotation.” BMC Bioinformatics 20, no. 1 (September 14, 2019): 473. https://doi.org/10.1186/s12859-019-3019-7.

and

Söding, Johannes. “Protein Homology Detection by HMM–HMM Comparison.” Bioinformatics 21, no. 7 (April 1, 2005): 951–60. https://doi.org/10.1093/bioinformatics/bti125.

"""
import argparse
import sys
import pyhmmer
from typing import Iterable, TextIO
import heapq
from multiprocessing import Pool
from domainator import __version__, RawAndDefaultsFormatter
from domainator.hmmer_search import read_hmms, compare_hmmer, traceback, HmmerHit

class _hmmer_compare_worker():
    def __init__(self, hmmer_targets, alignment=False, k=None, score_cutoff=float("-inf")):
        self.k = k
        self.hmmer_targets = hmmer_targets
        self.alignment = alignment
        self.score_cutoff = score_cutoff
    
    def __call__(self, input_profile):
        out_heap = []
        for target_dataset in self.hmmer_targets.values():
            for target_profile in target_dataset.values():
                score, backtrace, max_index, match_scores = compare_hmmer(input_profile, target_profile)
                if score >= self.score_cutoff:
                    if self.alignment:
                        alignment = traceback(input_profile,target_profile,backtrace,max_index, match_scores)
                    else:
                        alignment = None
                    
                    result = HmmerHit(score, input_profile.name.decode(), target_profile.name.decode(), alignment)
                    if (self.k is None) or (len(out_heap) < self.k):
                        heapq.heappush(out_heap, result)
                    else:
                        heapq.heappushpop(out_heap, result)
        out_heap.sort(reverse=True)
        return out_heap



def hmmer_compare(query_files:Iterable[str], reference_files:Iterable[str], out_handle:TextIO, score_cutoff:float, alignments:bool, k:int, cpu:int):
    references = read_hmms(hmm_files=reference_files) # list of lists of pyhmmer hmm objects

    worker = _hmmer_compare_worker(references, alignments, k, score_cutoff)

    sep="\t"
    print(sep.join(("query","reference","score")), file=out_handle) #TODO: how to write the alignment?
    
    for file in query_files:
        # file_name = os.path.basename(Path(file).stem)
        with Pool(processes=cpu) as pool:
            for hits in pool.imap(worker, pyhmmer.plan7.HMMFile(file), chunksize=1): # I tested some chunk sizes and it didn't seem to make a difference
                for hit in hits:
                    print(sep.join( (hit.query_name,hit.reference_name,f"{round(hit.score,2):.2f}") ), file=out_handle)
                    if alignments:
                        print(hit.alignment, file=out_handle)
                        print("\n\n", file=out_handle)


def main(argv):
    parser = argparse.ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)
    parser.add_argument('-i', '--input', type=str, required=True, nargs='+',
                        help="Input query files. One or more hmm text files with one or more hmmer3 profiles.")
    parser.add_argument('-r', "--reference", type=str, required=True, nargs='+',
                        help="Reference files. One or more hmm text files with one or more hmmer3 profiles.")

    parser.add_argument('--score_cutoff', type=float, default = 0,
                        help="Report alignments with scores greater than or equal to this.") #TODO: what is a reasonable cutoff?

    parser.add_argument('-k', type=int, required=False, default=None,
                        help="Include at most this many of the top hits for each query. Default: Include all hits.")

    parser.add_argument('-o', '--output', type=str, default=None,
                        help="File to write the scores and alignments to.")

    parser.add_argument('--alignments', action='store_true', default=False,
                        help="when activated, will write the alignments to the output.")
    
    parser.add_argument('--cpu', type=int, default=8, required=False,
                        help="how many cpu threads to use. Default: 8")


    params = parser.parse_args(argv)

    if params.output is None:
        out = sys.stdout
    else:
        out = open(params.output, "w")

    hmmer_compare(params.input, params.reference, out, params.score_cutoff, params.alignments, params.k, params.cpu)

    if params.output is not None:
        out.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
