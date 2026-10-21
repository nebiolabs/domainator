"""Aligns and calculates alignment scores between hmmer3 profiles 

Kind of like hhsearch except much slower and for hmmer3 profiles instead of hhsuite profiles.

Adapted from pseudocode in:
Steinegger, Martin, Markus Meier, Milot Mirdita, Harald Vöhringer, Stephan J. Haunsberger, and Johannes Söding. “HH-Suite3 for Fast Remote Homology Detection and Deep Protein Annotation.” BMC Bioinformatics 20, no. 1 (September 14, 2019): 473. https://doi.org/10.1186/s12859-019-3019-7.

and

Söding, Johannes. "Protein Homology Detection by HMM–HMM Comparison." Bioinformatics 21, no. 7 (April 1, 2005): 951–60. https://doi.org/10.1093/bioinformatics/bti125.

Amino acid, DNA, and RNA profiles are all supported, but the -i and -r sides must use
the same alphabet (DNA and RNA count as different alphabets).
"""
from jsonargparse import ActionConfigFile
from contextlib import ExitStack
from typing import Iterable, Optional, TextIO
import heapq
from domainator.cli.common import add_input, add_output, check_outputs, resolve_cpu, tool_parser
from domainator.outputs import output_file
from domainator.hmmer_search import read_hmms, compare_hmmer, traceback, HmmerHit, check_alphabets_match
from domainator._mp import make_pool
from domainator.hmmio import common_hmm_alphabet, iter_hmms_with_alphabet
from domainator.output_guardrails import add_max_output_gb_argument, enforce_output_limit, max_output_gb_to_bytes


def _format_hmmer_compare_hit(hit: HmmerHit) -> str:
    return "\t".join((hit.query_name, hit.reference_name, f"{round(hit.score, 2):.2f}")) + "\n"

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
                    
                    result = HmmerHit(score, input_profile.name, target_profile.name, alignment)
                    if (self.k is None) or (len(out_heap) < self.k):
                        heapq.heappush(out_heap, result)
                    else:
                        heapq.heappushpop(out_heap, result)
        out_heap.sort(reverse=True)
        return out_heap



def hmmer_compare(query_files:Iterable[str], *, references:Iterable[str], tsv_out_handle:Optional[TextIO], min_score:float, max_hits_per_query:int, cpu:int, alignments_out_handle:Optional[TextIO]=None, max_output_bytes=None, output_description: str = "hmmer_compare output"):
    """For each hit, write a TSV row (query, reference, score) to tsv_out_handle, and the same
    row followed by the alignment and a blank line to alignments_out_handle. Either handle can be None."""
    reference_hmms = read_hmms(hmm_files=references) # list of lists of pyhmmer hmm objects
    reference_alphabet = common_hmm_alphabet(references, role="reference")

    worker = _hmmer_compare_worker(reference_hmms, alignments_out_handle is not None, max_hits_per_query, min_score)

    header = "" if tsv_out_handle is None else "\t".join(("query", "reference", "score")) + "\n"
    header_bytes = len(header.encode("utf-8"))
    enforce_output_limit(
        projected_bytes=header_bytes,
        max_output_bytes=max_output_bytes,
        output_description=output_description,
        mitigation_options=["--max_hits_per_query", "--min_score", "--alignments_out"],
    )
    if tsv_out_handle is not None:
        tsv_out_handle.write(header)
    written_bytes = header_bytes
    
    for file in query_files:
        # Each query file is read once, so it can be standard input.
        file_alphabet, query_profiles = iter_hmms_with_alphabet(file, role="input")
        if file_alphabet is None:  # no profiles in this file
            continue
        check_alphabets_match(file_alphabet, reference_alphabet, str(file), "domainator hmmer_compare")
        with make_pool(processes=cpu) as pool:
            for hits in pool.imap(worker, query_profiles, chunksize=1): # I tested some chunk sizes and it didn't seem to make a difference
                for hit in hits:
                    row = _format_hmmer_compare_hit(hit)
                    alignment = None if alignments_out_handle is None else row + hit.alignment + "\n\n"
                    hit_bytes = (0 if tsv_out_handle is None else len(row.encode("utf-8"))) + (0 if alignment is None else len(alignment.encode("utf-8")))
                    enforce_output_limit(
                        projected_bytes=written_bytes + hit_bytes,
                        max_output_bytes=max_output_bytes,
                        output_description=output_description,
                        mitigation_options=["--max_hits_per_query", "--min_score", "--alignments_out"],
                    )
                    if tsv_out_handle is not None:
                        tsv_out_handle.write(row)
                    if alignment is not None:
                        alignments_out_handle.write(alignment)
                    written_bytes += hit_bytes


def build_parser():
    parser = tool_parser(__doc__)
    add_input(parser, "Input query files. One or more hmm text files with one or more hmmer3 profiles.")
    parser.add_argument('-r', "--reference", type=str, required=True, nargs='+',
                        help="Reference files. One or more hmm text files with one or more hmmer3 profiles.")

    parser.add_argument('--min_score', type=float, default = 0,
                        help="Report alignments with scores greater than or equal to this. Scores are on "
                             "the score scale of the input alphabet: DNA and RNA scores run at roughly "
                             "half the amino acid scale, so a threshold tuned on proteins is about twice "
                             "as strict on nucleotide profiles.")

    parser.add_argument('--max_hits_per_query', type=int, required=False, default=None,
                        help="Include at most this many of the top hits for each query. Default: Include all hits.")

    add_output(parser, "tsv", "Write the scores to this tab-separated file.", primary=True)
    add_output(parser, "alignments", "Also write each hit's alignment to this text file.")
    
    parser.add_argument('--cpu', type=int, default=0, required=False,
                        help="how many cpu threads to use. 0: one per physical core.")

    add_max_output_gb_argument(parser)

    parser.add_argument('--config', action=ActionConfigFile)
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params, ["tsv_out", "alignments_out"])

    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)
    with ExitStack() as stack:
        tsv_out = None if params.tsv_out is None else stack.enter_context(output_file(params.tsv_out))
        alignments_out = None if params.alignments_out is None else stack.enter_context(output_file(params.alignments_out))
        hmmer_compare(params.input, references=params.reference, tsv_out_handle=tsv_out, alignments_out_handle=alignments_out,
                      min_score=params.min_score, max_hits_per_query=params.max_hits_per_query,
                      cpu=resolve_cpu(params.cpu), max_output_bytes=max_output_bytes)
