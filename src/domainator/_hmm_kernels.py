"""numba kernels for profile-profile alignment in domainator hmmer_search (see that module).

They are in their own module so that importing hmmer_search doesn't load numba.
"""
from typing import Tuple

import numba as nb
from numba import jit
import numpy as np

from domainator.hmmer_search import _SDG, _SGD, _SIM, _SMI, _SMM, _STOP_FLAG


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
