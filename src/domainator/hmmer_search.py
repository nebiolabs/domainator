"""Filter an hmm file based on similarity to reference hmm profiles

If you want a distance matrix, use seq_dist.py with .hmm files for both inputs and outputs.

If you want something more like hhsearch or blast, use hmmer_compare.py, which uses 
the same algorithm, but a different interface (can output the same input profile multiple times).


Adapted from pseudocode in:
Steinegger, Martin, Markus Meier, Milot Mirdita, Harald Vöhringer, Stephan J. Haunsberger, and Johannes Söding. “HH-Suite3 for Fast Remote Homology Detection and Deep Protein Annotation.” BMC Bioinformatics 20, no. 1 (September 14, 2019): 473. https://doi.org/10.1186/s12859-019-3019-7.

and

Söding, Johannes. “Protein Homology Detection by HMM–HMM Comparison.” Bioinformatics 21, no. 7 (April 1, 2005): 951–60. https://doi.org/10.1093/bioinformatics/bti125.

Amino acid, DNA, and RNA profiles are all supported, but the -i and -r sides must use
the same alphabet (DNA and RNA count as different alphabets). Profile-profile scores are
on the scale of whichever alphabet is in use, so the default --score_cutoff differs
between amino acid and nucleotide profiles.
"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module='numpy') # suppress "UserWarning: The value of the smallest subnormal for <class 'numpy.float64'> type is zero."
from jsonargparse import ArgumentParser, ActionConfigFile
from dataclasses import dataclass
import sys
import io
import os
import numpy as np
import pyhmmer
from typing import List, Iterable, TextIO, Tuple, Union, Dict, Optional, BinaryIO
import heapq
from domainator import __version__, RawAndDefaultsFormatter
from domainator.utils import make_pool, pyhmmer_decode, read_hmms, alphabet_name, alphabet_score_scale, common_hmm_alphabet, iter_hmms_with_alphabet
from domainator.output_guardrails import add_max_output_gb_argument, enforce_output_limit, max_output_gb_to_bytes, OutputSizeLimitExceeded, make_temporary_output_path
from numba import jit
import numba as nb
import psutil


def _input_display_name(input_file) -> str:
    if isinstance(input_file, (str, os.PathLike)):
        return os.fspath(input_file)
    return getattr(input_file, "name", repr(input_file))


def _estimate_hmm_profile_size_bytes(profile: pyhmmer.plan7.HMM) -> int:
    buffer = io.BytesIO()
    profile.write(buffer)
    return len(buffer.getvalue())


DEFAULT_SCORE_CUTOFF_AMINO = 10.0
"""Default --score_cutoff for amino acid profiles."""


def default_score_cutoff(alphabet) -> float:
    """Alphabet-appropriate default for --score_cutoff.

    Profile-profile scores are sums of natural-log-odds column scores, and the
    per-column ceiling depends on the alphabet: about log(4)=1.39 for DNA and RNA
    against about log(20)=3.0 for amino acids. A threshold tuned on proteins is
    therefore roughly twice as strict on nucleotides, where it silently returns an
    empty result that looks the same as "no homologs found". Scaling the default
    keeps it at the same relative position on each alphabet's score scale.

    Only the default is scaled: a --score_cutoff the user supplied explicitly is
    always taken literally.
    """
    if alphabet is None:
        return DEFAULT_SCORE_CUTOFF_AMINO
    return DEFAULT_SCORE_CUTOFF_AMINO * alphabet_score_scale(alphabet)


def check_alphabets_match(input_alphabet, reference_alphabet, input_name, tool_name) -> None:
    """Raise if the -i and -r sides of a profile-profile comparison disagree.

    Called in the parent process, before any worker pool is created, so that a
    mismatch reads as one clear line instead of a multiprocessing traceback.
    """
    if input_alphabet is None or reference_alphabet is None:
        return
    if input_alphabet == reference_alphabet:
        return
    raise ValueError(
        f"Input (-i) profiles in '{input_name}' are {alphabet_name(input_alphabet)}, but the "
        f"reference (-r) profiles are {alphabet_name(reference_alphabet)}. {tool_name} aligns "
        f"profiles column by column, so both sides must use the same alphabet. Rebuild one side "
        f"with hmmer_build.py --alphabet, or use domainate.py / domain_search.py to search "
        f"nucleotide profiles against sequences."
    )


#from Algorithm2 pseudocode in:
# Steinegger, Martin, Markus Meier, Milot Mirdita, Harald Vöhringer, Stephan J. Haunsberger, and Johannes Söding. “HH-Suite3 for Fast Remote Homology Detection and Deep Protein Annotation.” BMC Bioinformatics 20, no. 1 (September 14, 2019): 473. https://doi.org/10.1186/s12859-019-3019-7.

_SMM:int = 0
_SMI:int = 1
_SIM:int = 2
_SDG:int = 3
_SGD:int = 4
_TRACEBACK_DTYPE=np.uint64
_STOP_FLAG = 100000000000000000 #np.iinfo(_TRACEBACK_DTYPE).max
#             18446744073709551615

@jit(nb.types.Tuple((nb.float32,nb.uint64))(nb.float32, nb.float32, nb.uint64, nb.uint64),cache=True, nopython=True)
def max2(sMM:float, sXY:float, layer1:int, layer2:int) -> Tuple[float, int]:
    if sMM > sXY:
        score = sMM
        bt = layer1
    else:
        score = sXY
        bt = layer2
    return score, bt

@jit(nb.types.Tuple((nb.float32, nb.uint64))(nb.float32, nb.float32, nb.float32, nb.float32, nb.float32, nb.int64),cache=True, nopython=True)
def max6(sMM: float, sMI: float, sIM: float, sDG: float, sGD: float, global_mode:int = 0) -> Tuple[float, int]:
    """

    Args:
        sMM (float): 
        sMI (float): 
        sIM (float): 
        sDG (float): 
        sGD (float):         

        sSTOP (float, optional): Minimum possible alignment score. 0 for local alignment, -inf for global alignment. Defaults to 0. 

    Returns:
        Tuple[float, int]: score, backtrace_layer

    """
    score = -1000000000000000000.0
    bt = None
    sSTOP = 0 #if global_mode == 0 else float("-inf")
    
    if sSTOP > sMM and global_mode == 0:
        score = sSTOP
        bt = _STOP_FLAG
    else:
        score = sMM
        bt = _SMM
    
    if sIM > score:
        score = sIM
        bt = _SIM
    if sMI > score:
        score = sMI
        bt = _SMI
    if sDG > score:
        score = sDG
        bt = _SDG
    if sGD > score:
        score = sGD
        bt = _SGD
    
    return  score, bt
    
@jit(nb.float32(nb.float32[:], nb.float32[:], nb.float32[:]),cache=True, nopython=True)
def Saa(q:np.array, r:np.array, background:np.array) -> float: #TODO: test
    """calculate a similarity score for two residue distributions

    Alphabet agnostic: the numba signature leaves the dimension open, so the
    arrays are length alphabet.K, which is 20 for amino acids and 4 for DNA and RNA.

    Args:
        q (np.array): array size K, query residue probabilities
        r (np.array): array size K, reference residue probabilities
        background (np.array): array size K, background residue probabilities

    Returns:
        float: natural-log-odds score for the column pair. Bounded above by
            log(1/min(background)), which is about log(4)=1.39 for a nucleotide
            alphabet against about 3 for amino acids, so scores from different
            alphabets are not comparable. See utils.alphabet_score_scale.
    """

    return np.log(np.sum(np.divide(np.multiply(q, r), background)))

@jit(nb.void(nb.float32[:,:], nb.float32[:,:], nb.float32[:,:], nb.float32[:], nb.float32[:,:], nb.float32[:,:,:], nb.uint64[:,:,:], nb.float32[:,:], nb.int64),cache=True, nopython=True)
def compare_hmmer_inner(query_transitions:np.array, query_match_emissions:np.array, reference_transitions:np.array, background:np.array, reference_match_emissions:np.array, scores:np.array, backtrace:np.array, match_scores:np.array, global_mode=0) -> None:
    TMM=0
    TMI=1
    TMD=2
    TIM=3
    TII=4
    TDM=5
    TDD=6
    # 0: Mn -> Mn+1 
    # 1: Mn -> In+1 
    # 2: Mn -> Dn+1 
    # 3: In -> Mn+1 
    # 4: In -> In+1 
    # 5: Dn -> Mn+1 
    # 6: Dn -> Dn+1 


    backtrace[:,0,:] = _STOP_FLAG # stop backtrace at the left column
    backtrace[0,:,:] = _STOP_FLAG # stop backtrace at the top row
    # SMM = 0
    # SMI = 1
    # SIM = 2
    # SDG = 3
    # SGD = 4

    for qi in range(1,query_match_emissions.shape[0]):
        for ri in range(1,reference_match_emissions.shape[0]):

            # SMM
            scores[qi,ri,_SMM], tb_layer = \
                max6( # find the best previous state
                    scores[qi-1,ri-1,_SMM] - query_transitions[qi-1, TMM] - reference_transitions[ri-1, TMM], # subtracting negative logs is the same as adding logs, which is the same as multiplying unlogs
                    scores[qi-1,ri-1,_SMI] - query_transitions[qi-1, TMM] - reference_transitions[ri-1, TIM],
                    scores[qi-1,ri-1,_SIM] - query_transitions[qi-1, TIM] - reference_transitions[ri-1, TMM], 
                    scores[qi-1,ri-1,_SDG] - query_transitions[qi-1, TDM] - reference_transitions[ri-1, TMM], 
                    scores[qi-1,ri-1,_SGD] - query_transitions[qi-1, TMM] - reference_transitions[ri-1, TDM],
                    0 # global alignment
                 ) 
            backtrace[qi,ri,_SMM] = tb_layer
            #Saa_score = Saa(query_match_emissions[qi,:], reference_match_emissions[ri,:], background)
            #print(Saa_score)
            match_score = Saa(query_match_emissions[qi,:], reference_match_emissions[ri,:], background)
            match_scores[qi,ri] = match_score
            scores[qi,ri,_SMM] += match_score
            # could add a secondary structure term here if available
            
            # SGD
            scores[qi,ri,_SGD], tb_layer = \
                max2(
                    scores[qi,ri-1,_SMM] - reference_transitions[ri-1, TMD],
                    scores[qi,ri-1,_SGD] - reference_transitions[ri-1, TDD],
                    _SMM,
                    _SGD
                )
            backtrace[qi,ri,_SGD] = tb_layer

            # SIM
            scores[qi,ri,_SIM], tb_layer = \
                max2(
                    scores[qi,ri-1,_SMM] - query_transitions[qi, TMI] - reference_transitions[ri-1, TMM],
                    scores[qi,ri-1,_SIM] - query_transitions[qi, TII] - reference_transitions[ri-1, TMM],
                    _SMM,
                    _SIM
                )
            backtrace[qi,ri,_SIM] = tb_layer


            # SDG
            scores[qi,ri,_SDG], tb_layer = \
                max2(
                    scores[qi-1,ri,_SMM] - query_transitions[qi-1, TMD],
                    scores[qi-1,ri,_SDG] - query_transitions[qi-1, TDD],
                    _SMM,
                    _SDG
                )
            backtrace[qi,ri,_SDG] = tb_layer

            # SMI
            scores[qi,ri,_SMI], tb_layer = \
                max2(
                    scores[qi-1,ri,_SMM] - query_transitions[qi-1, TMM] - reference_transitions[ri, TMI],
                    scores[qi-1,ri,_SMI] - query_transitions[qi-1, TMM] - reference_transitions[ri, TII],
                    _SMM,
                    _SMI
                )
            backtrace[qi,ri,_SMI] = tb_layer

def compare_hmmer(qhmm:pyhmmer.plan7.HMM, rhmm:pyhmmer.plan7.HMM) -> Tuple[float, np.array, Tuple[int,int,int]]:
    """Run the viterbi alignment algorithm to compare two plan7 HMM profiles

    Args:
        qhmm (pyhmmer.plan7.HMM): query HMM object
        rhmm (pyhmmer.plan7.HMM): reference (target) HMM object

    Raises:
        ValueError: if hmms are not of the same alphabet

    Returns:
        Tuple[float, np.array, Tuple[int,int,int]]: score, backtrace, max_index (query_idx, ref_idx, layer)
    """
    if qhmm.alphabet != rhmm.alphabet:
        raise ValueError(
            f"Error, cannot compare hmms with different alphabets: query "
            f"'{pyhmmer_decode(qhmm.name)}' is {alphabet_name(qhmm.alphabet)}, but reference "
            f"'{pyhmmer_decode(rhmm.name)}' is {alphabet_name(rhmm.alphabet)}. Note that DNA "
            f"and RNA are distinct alphabets."
        )
    background = np.asarray(pyhmmer.plan7.Background(qhmm.alphabet).residue_frequencies) #TODO: could just store this as a constant

    # raw probabilities
    query_match_emissions=np.asarray(qhmm.match_emissions) # M + 1, K. Row 0 is entry probabilities, so emissions are unused, but it still needs to sum to 1.0, so it is [1.0, 0.0, 0.0, ...], but don't use it for any calcs
    reference_match_emissions=np.asarray(rhmm.match_emissions) # M + 1, K. Row 0 is entry probabilities, so emissions are unused, but it still needs to sum to 1.0, so it is [1.0, 0.0, 0.0, ...], but don't use it for any calcs
    
    # query_insert_emissions=np.asarray(qhmm.insert_emissions) # row 0 is insert emissions for left-gaps (before the first match emission) 
    # reference_insert_emissions=np.asarray(rhmm.insert_emissions) # row 0 is insert emissions for left-gaps (before the first match emission)

    with np.errstate(divide='ignore'):
        query_transitions=-1 * np.log(np.asarray(qhmm.transition_probabilities))  # M + 1, 7 probabilities are negative natural logs
        reference_transitions=-1 * np.log(np.asarray(rhmm.transition_probabilities)) # M + 1, 7 probabilities are negative natural logs
    
    # print(np.asarray(qhmm.transition_probabilities))

    scores = np.zeros((query_match_emissions.shape[0], reference_match_emissions.shape[0],5),dtype=np.float32) # dim3 = SMM, SMI, SIM, SDG, SGD, for SXY, X is query, Y is reference
    backtrace = np.zeros((query_match_emissions.shape[0], reference_match_emissions.shape[0],5),dtype=_TRACEBACK_DTYPE) #score_row, score_column, state_type (dim3 of scores array) TODO: should we account for alternative alignments, like by ties in the max fuctions?. TODO: what order of dimensions is best for caching?
    match_scores = np.zeros((query_match_emissions.shape[0], reference_match_emissions.shape[0]),dtype=np.float32) # scores for match states only, for backtrace. Can probably be np.empty.
    # values are flags to the previous state type. There is only one direction from each previous state, so directionality is implied in the state type


    compare_hmmer_inner(query_transitions, query_match_emissions, reference_transitions, background, reference_match_emissions, scores, backtrace, match_scores, 0)

    max_index = np.unravel_index(np.argmax(scores),scores.shape)  
    score = scores[max_index]
    

    return score, backtrace, max_index, match_scores

def traceback(qhmm,rhmm,backtrace,trace_start,match_scores):
    POSITION_PADDING = 10

    qcons = qhmm.consensus
    rcons = rhmm.consensus
    qend = trace_start[0] # int 
    rend = trace_start[1] # int
    qstart = qend
    rstart = rend

    ptr = trace_start # 3-tuple q_pos, r_pos, tb_layer

    prev = backtrace[ptr]
    qline=list()
    midline=list()
    rline=list()
    # SMM = 0
    # SMI = 1
    # SIM = 2
    # SDG = 3
    # SGD = 4

    # hhsearch style alignment output symbols
    # https://github.com/soedinglab/hh-suite/wiki#hmm-hmm-pairwise-alignments
    # = : column score below -1.5
    # - : column score between -1.5 and -0.5
    # . : column score between -0.5 and +0.5
    # + : column score between +0.5 and +1.5
    # | : column score above   +1.5
    # The published cutoffs above are on the amino acid score scale. Saa scores are
    # natural-log odds against the alphabet background, and a nucleotide column tops
    # out near log(4)=1.39, so on DNA and RNA profiles '|' would be unreachable and
    # the whole midline would collapse into '.' and '+'. Rescaling keeps each symbol
    # at the same relative position on the alphabet's own scale. The scale is exactly
    # 1.0 for amino acids, so protein output is unchanged.
    score_scale = alphabet_score_scale(qhmm.alphabet) # compare_hmmer already checked that rhmm matches
    very_bad_cutoff = -1.5 * score_scale
    bad_cutoff = -0.5 * score_scale
    neutral_cutoff = 0.5 * score_scale
    good_cutoff = 1.5 * score_scale


    while True:
        ptr_layer = ptr[2]
        qpos = ptr[0]-1 #-1 because the backtrace array has an extra row at the front
        rpos = ptr[1]-1 #-1 because the backtrace array has an extra column at the front
        if ptr_layer == _SMM:
            qchar = qcons[qpos]
            rchar = rcons[rpos]
            qline.append(qchar)
            rline.append(rchar)
            next_ptr = (ptr[0]-1, ptr[1]-1, prev)
            qstart = ptr[0]
            rstart = ptr[1]

            score = match_scores[qpos+1,rpos+1]
            if score < very_bad_cutoff:
                midline.append("=")
            elif score < bad_cutoff:
                midline.append("-")
            elif score < neutral_cutoff:
                midline.append(".")
            elif score < good_cutoff:
                midline.append("+")
            else:
                midline.append("|")
        elif ptr_layer == _SMI:
            qline.append(qcons[qpos])
            rline.append("-")
            midline.append(" ")
            qstart = ptr[0]
            next_ptr = (ptr[0]-1, ptr[1], prev)
        elif ptr_layer == _SIM:
            qline.append("-")
            rline.append(rcons[rpos])
            midline.append(" ")
            rstart = ptr[1]
            next_ptr = (ptr[0], ptr[1]-1, prev)
        elif ptr_layer == _SDG:
            qline.append(qcons[qpos])
            rline.append("-")
            midline.append(" ")
            qstart = ptr[0]
            next_ptr = (ptr[0]-1, ptr[1], prev)
        elif ptr_layer == _SGD:
            qline.append("-")
            rline.append(rcons[rpos])
            midline.append(" ")
            rstart = ptr[1]
            next_ptr = (ptr[0], ptr[1]-1, prev)
        
        if prev == _STOP_FLAG:
            break
        ptr = next_ptr
        prev = backtrace[ptr]
    return str(qstart).ljust(POSITION_PADDING) + ''.join(reversed(qline)) + str(qend).rjust(POSITION_PADDING) + "\n" + \
           ' '*POSITION_PADDING + ''.join(reversed(midline)) + ' ' * POSITION_PADDING + "\n" + \
           str(rstart).ljust(POSITION_PADDING) + ''.join(reversed(rline)) + str(rend).rjust(POSITION_PADDING)

@dataclass
class HmmerHit():
    score: float
    query_name: str
    reference_name: str
    alignment: str
    profile: Optional[pyhmmer.plan7.HMM] = None
    
    def __eq__(self, other):
        return self.score == other.score

    def __ne__(self, other):
        return self.score != other.score

    def __lt__(self, other):
        return self.score < other.score

    def __le__(self, other):
        return self.score <= other.score

    def __gt__(self, other):
        return self.score > other.score

    def __ge__(self, other):
        return self.score >= other.score


class _hmmer_search_worker():
    def __init__(self, hmmer_targets, score_cutoff=float("-inf")):
        self.hmmer_targets = hmmer_targets
        self.score_cutoff = score_cutoff
    
    def __call__(self, input_profile):
        best_result = None
        for target_dataset in self.hmmer_targets.values():
            for target_profile in target_dataset.values():
                score, _, _, _ = compare_hmmer(input_profile, target_profile)
                if score >= self.score_cutoff:
                    alignment = None
                    new_result = HmmerHit(score, pyhmmer_decode(input_profile.name), pyhmmer_decode(target_profile.name), alignment, input_profile)
                    if best_result is None or new_result > best_result:
                        best_result = new_result
        return best_result


def hmmer_search(input_files:Iterable[str], reference_files:Iterable[str], hmmer_handle:BinaryIO, score_cutoff:Optional[float], max_hits:int, cpu:int, max_output_bytes: Optional[int] = None, output_description: str = "hmmer_search HMM output"):
    references = read_hmms(hmm_files=reference_files) # list of lists of pyhmmer hmm objects
    reference_alphabet = common_hmm_alphabet(reference_files, role="reference")

    if score_cutoff is None:
        score_cutoff = default_score_cutoff(reference_alphabet)
        alphabet_label = alphabet_name(reference_alphabet) if reference_alphabet is not None else "amino"
        print(f"Using --score_cutoff {score_cutoff:g} (the default for {alphabet_label} profiles).", file=sys.stderr)

    worker = _hmmer_search_worker(references, score_cutoff)

    
    def run_comparison():
        for input_file in input_files:
            # Read the first profile up front so that an alphabet mismatch is reported
            # here, in the parent process, rather than as a traceback out of a pool worker.
            input_alphabet, input_profiles = iter_hmms_with_alphabet(input_file, role="input")
            if input_alphabet is None:
                continue
            check_alphabets_match(input_alphabet, reference_alphabet, _input_display_name(input_file), "hmmer_search.py")
            with make_pool(processes=cpu) as pool:
                for hit in pool.imap(worker, input_profiles, chunksize=1): # I tested some chunk sizes and it didn't seem to make a difference
                    if hit is not None:
                        yield hit

    written_bytes = 0
    mitigation_options = ["--max_hits", "--score_cutoff"]

    if max_hits is None:
        for hit in run_comparison():
            profile_bytes = _estimate_hmm_profile_size_bytes(hit.profile) if max_output_bytes is not None else 0
            enforce_output_limit(
                projected_bytes=written_bytes + profile_bytes,
                max_output_bytes=max_output_bytes,
                output_description=output_description,
                mitigation_options=mitigation_options,
            )
            hit.profile.write(hmmer_handle)
            written_bytes += profile_bytes
    else:
        out_heap = []
        for hit in run_comparison():
            if len(out_heap) < max_hits:
                heapq.heappush(out_heap, hit)
            else:
                heapq.heappushpop(out_heap, hit)
        out_heap.sort(reverse=True)
        for hit in out_heap:
            profile_bytes = _estimate_hmm_profile_size_bytes(hit.profile) if max_output_bytes is not None else 0
            enforce_output_limit(
                projected_bytes=written_bytes + profile_bytes,
                max_output_bytes=max_output_bytes,
                output_description=output_description,
                mitigation_options=mitigation_options,
            )
            hit.profile.write(hmmer_handle)
            written_bytes += profile_bytes


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)
    parser.add_argument('-i', '--input', type=str, required=False, nargs='+',
                        help="Input query files. One or more hmm text files with one or more hmmer3 profiles. Default: stdin.\nIf -i and -r are different, -i will typically be the bigger of the two, as it is what is being filtered.")
    
    parser.add_argument('-r', "--reference", type=str, required=True, nargs='+',
                        help="Reference files. One or more hmm text files with one or more hmmer3 profiles.\nIf -i and -r are different, -r will typically be the smaller of the two, as it is used as a filtering criterion.")
    
    parser.add_argument('-o', '--output', type=str, required=False, default=None,
                        help=".hmm file to write hit profiles to. Default: stdout")

    parser.add_argument('--score_cutoff', type=float, default = None,
                        help=f"Report alignments with scores greater than or equal to this. "
                             f"Default: {DEFAULT_SCORE_CUTOFF_AMINO:g} for amino acid profiles, and the "
                             f"equivalent point on the nucleotide score scale for DNA and RNA profiles "
                             f"(about {default_score_cutoff(pyhmmer.easel.Alphabet.dna()):.2f}). A value supplied here is "
                             f"used as-is, without rescaling.")

    parser.add_argument('--max_hits', type=int, default=None,
                        help="The maximum number of hmms returned by the search. Prioritized by bitscore of best scoring profile. Default: return all hits passing the score cutoff.")
    
    parser.add_argument('--cpu', type=int, default=0, required=False,
                        help="how many cpu threads to use. Default: use all available cores")

    add_max_output_gb_argument(parser)

    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)

    if params.cpu <= 0:
        cpus = psutil.cpu_count(logical=False)
    else:
        cpus = params.cpu

    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)
    output_description = f"hmmer_search HMM output '{params.output if params.output is not None else 'stdout'}'"
    temp_output_path = None
    if params.output is None:
        out = sys.stdout.buffer
    else:
        temp_output_path = make_temporary_output_path(params.output)
        out = open(temp_output_path, "wb")
    
    if params.input is None:
        input_files = [sys.stdin.buffer] # pyhmmer needs a binary stream
    else:
        input_files = params.input

    try:
        hmmer_search(input_files, params.reference, out, params.score_cutoff, params.max_hits, cpus, max_output_bytes=max_output_bytes, output_description=output_description)
    except OutputSizeLimitExceeded as exc:
        if params.output is not None:
            out.close()
            if temp_output_path is not None and os.path.exists(temp_output_path):
                os.remove(temp_output_path)
        raise SystemExit(str(exc)) from None
    except Exception:
        if params.output is not None:
            out.close()
            if temp_output_path is not None and os.path.exists(temp_output_path):
                os.remove(temp_output_path)
        raise

    if params.output is not None:
        out.close()
        os.replace(temp_output_path, params.output)

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
