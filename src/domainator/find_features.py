"""Annotate sequence files with feature detection algorithms

For each CDS or protein in the input file, run feature detection algorithms to 
identify transmembrane regions, signal peptides, and other sequence features.

Supported algorithms:
    - TMbed: Transmembrane protein prediction using ProtT5 language model embeddings
    - MOTIF: PROSITE-style pattern matching for exact motif search
"""

import sys
import re
import math
import warnings
warnings.filterwarnings("ignore", module='numpy')

from typing import List, Dict, Tuple, Optional, NamedTuple, Iterator
from jsonargparse import ArgumentParser, ActionConfigFile

from domainator.Bio.SeqRecord import SeqRecord
from domainator.Bio.SeqFeature import SeqFeature, FeatureLocation
from domainator import utils, DOMAIN_FEATURE_NAME
from domainator.utils import get_cds_unique_name, parse_seqfiles, write_genbank, BooleanEvaluator
from domainator import __version__, RawAndDefaultsFormatter
from domainator.domainate import clean_rec, prodigal_CDS_annotate

# TMbed prediction class labels
# From TMbed format 3 (undirected):
# B: Transmembrane beta strand
# H: Transmembrane alpha helix
# S: Signal peptide
# i: Non-Transmembrane, inside
# o: Non-Transmembrane, outside
TMBED_LABELS = {
    0: ('B', 'transmembrane_beta_strand', 'transmembrane beta strand'),
    1: ('B', 'transmembrane_beta_strand', 'transmembrane beta strand'),  # b -> B (undirected)
    2: ('H', 'transmembrane_alpha_helix', 'transmembrane alpha helix'),
    3: ('H', 'transmembrane_alpha_helix', 'transmembrane alpha helix'),  # h -> H (undirected)
    4: ('S', 'signal_peptide', 'signal peptide'),
    5: ('i', 'non_transmembrane_inside', 'non-transmembrane, inside'),
    6: ('o', 'non_transmembrane_outside', 'non-transmembrane, outside'),
}


class FeatureResult(NamedTuple):
    """Result from a feature detection algorithm"""
    name: str
    """feature name"""
    desc: str
    """feature description"""
    score: float
    """confidence score (0-100 for probability-based)"""
    evalue: float
    """pseudo-evalue: e^(-score/100 * scale_factor)"""
    start: int
    """start position on sequence (0-based)"""
    end: int
    """end position on sequence (0-based, exclusive)"""
    program: str
    """name of program that made the prediction"""
    identity: float
    """same as score for these algorithms"""
    probabilities: Optional[Dict[str, List[float]]] = None
    """per-residue probabilities for each class"""


def score_to_pseudo_evalue(score: float) -> float:
    """
    Convert a probability-based score (0-100) to a pseudo-evalue.
    
    Uses e^(-score) to produce evalue-like numbers. For example:
    - 90% confidence → e^(-90) ≈ 1.2e-39
    - 50% confidence → e^(-50) ≈ 1.9e-22
    
    Note: These are not true evalues (which depend on database size),
    but provide a familiar scale for interpreting confidence.
    
    Args:
        score: Confidence score (0-100, where 100 is best)
    
    Returns:
        Pseudo-evalue where lower is better confidence
    """
    return math.exp(-score)


def sanitize_motif_name(pattern: str) -> str:
    """
    Convert a PROSITE pattern to a name safe for boolean evaluators.
    
    This is a convenience wrapper around BooleanEvaluator.sanitize_identifier.
    
    Args:
        pattern: The PROSITE pattern string
    
    Returns:
        A sanitized name suitable for boolean expressions
    """
    return BooleanEvaluator.sanitize_identifier(pattern)


def prosite_to_regex(pattern: str) -> str:
    """
    Convert a PROSITE pattern to a Python regex pattern.
    
    PROSITE syntax:
    - x: any amino acid
    - [ABC]: any of A, B, C
    - {ABC}: any except A, B, C  
    - -: separator (optional)
    - (n): repeat n times
    - (n,m): repeat n to m times
    - <: N-terminal anchor
    - >: C-terminal anchor
    - .: end of pattern (optional)
    
    Args:
        pattern: PROSITE pattern string
    
    Returns:
        Python regex pattern string
    """
    # Remove trailing period and whitespace
    pattern = pattern.strip().rstrip('.')
    
    result = []
    i = 0
    
    while i < len(pattern):
        c = pattern[i]
        
        if c == '<':
            # N-terminal anchor
            result.append('^')
            i += 1
        elif c == '>':
            # C-terminal anchor
            result.append('$')
            i += 1
        elif c == '-':
            # Separator, ignore
            i += 1
        elif c == 'x' or c == 'X':
            # Any amino acid
            result.append('.')
            i += 1
        elif c == '[':
            # Character class (allowed residues)
            end = pattern.index(']', i)
            result.append(pattern[i:end+1].upper())
            i = end + 1
        elif c == '{':
            # Negated character class (disallowed residues)
            end = pattern.index('}', i)
            chars = pattern[i+1:end].upper()
            result.append(f'[^{chars}]')
            i = end + 1
        elif c == '(':
            # Repetition
            end = pattern.index(')', i)
            rep = pattern[i+1:end]
            if ',' in rep:
                n, m = rep.split(',')
                result.append(f'{{{n.strip()},{m.strip()}}}')
            else:
                result.append(f'{{{rep.strip()}}}')
            i = end + 1
        elif c.isalpha():
            # Single amino acid
            result.append(c.upper())
            i += 1
        else:
            # Skip unknown characters
            i += 1
    
    return ''.join(result)


def search_motif(
    sequence: str,
    pattern: str,
    pattern_name: Optional[str] = None,
) -> List[FeatureResult]:
    """
    Search for a PROSITE motif in a protein sequence.
    
    Args:
        sequence: Amino acid sequence to search
        pattern: PROSITE pattern string
        pattern_name: Optional name for the pattern (if None, uses sanitized pattern)
    
    Returns:
        List of FeatureResult objects for each match
    """
    results = []
    
    # Convert pattern to regex
    try:
        regex_pattern = prosite_to_regex(pattern)
        regex = re.compile(regex_pattern, re.IGNORECASE)
    except (ValueError, re.error) as e:
        warnings.warn(f"Invalid PROSITE pattern '{pattern}': {e}")
        return results
    
    # Determine name
    if pattern_name is None:
        name = sanitize_motif_name(pattern)
    else:
        name = pattern_name
    
    # Search for all non-overlapping matches
    sequence =sequence.upper().rstrip("*")  # Remove trailing stop codon if present

    for match in regex.finditer(sequence):
        results.append(FeatureResult(
            name=name,
            desc=pattern,  # Use raw pattern as description
            score=100.0,   # Exact match = 100%
            evalue=0.0,    # Exact match = 0
            start=match.start(),
            end=match.end(),
            program="MOTIF",
            identity=100.0,  # Exact match = 100%
            probabilities=None,
        ))
    
    return results


def search_all_motifs(
    sequence: str,
    motifs: List[str],
) -> List[FeatureResult]:
    """
    Search for multiple PROSITE motifs in a protein sequence.
    
    Args:
        sequence: Amino acid sequence to search
        motifs: List of PROSITE pattern strings
    
    Returns:
        List of FeatureResult objects for all matches
    """
    results = []
    for pattern in motifs:
        results.extend(search_motif(sequence, pattern))
    return results


def get_tmbed_predictor(gpu_device: Optional[str] = "cuda:0", cpu_threads: int = 0, model_dir: Optional[str] = None):
    """
    Initialize TMbed models and encoder.
    
    Args:
        gpu_device: Device to use for GPU acceleration (e.g., "cuda:0", "0", or None for CPU)
        cpu_threads: Number of CPU threads (0 = use all available)
        model_dir: Directory to download/load the ProtT5 model from. If None, uses
                   tmbed's default (HF_HOME env var, or tmbed/models/t5/ in package)
    
    Returns:
        Tuple of (encoder, models, decoder, device, use_gpu)
    """
    try:
        import torch
        from tmbed.embed import T5Encoder
        from tmbed.tmbed import load_models, load_encoder, init
        from tmbed.viterbi import Decoder
        from tmbed.utils import seed_all
    except ImportError:
        raise ImportError(
            "TMbed is required for transmembrane prediction. "
            "Install it with: pip install git+https://github.com/BernhoferM/TMbed.git"
        )
    
    # Determine if we should use GPU
    use_gpu = False
    if gpu_device is not None and gpu_device.lower() != "none":
        if torch.cuda.is_available():
            use_gpu = True
            if gpu_device.isdigit():
                gpu_device = f"cuda:{gpu_device}"
        else:
            warnings.warn(f"GPU device {gpu_device} requested but CUDA is not available. Using CPU.")
    
    # Initialize TMbed
    num_threads = cpu_threads if cpu_threads > 0 else (torch.get_num_threads() if hasattr(torch, 'get_num_threads') else 1)
    init(use_gpu, num_threads)
    
    # Set device
    device = torch.device(gpu_device if use_gpu else 'cpu')
    
    # Load models
    encoder = load_encoder(use_gpu, model_dir)  # None uses default model path
    models = load_models()
    decoder = Decoder()
    
    return encoder, models, decoder, device, use_gpu


def run_tmbed_batch(
    sequences: List[Tuple[str, str]],  # List of (id, sequence) tuples
    encoder,
    models,
    decoder,
    device,
    cpu_fallback: bool = True,
    batch_size: int = 4000,
) -> Dict[str, Tuple]:
    """
    Run TMbed prediction on a batch of sequences.
    
    Args:
        sequences: List of (id, sequence) tuples
        encoder: TMbed T5 encoder
        models: TMbed CNN models
        decoder: TMbed Viterbi decoder
        device: PyTorch device
        cpu_fallback: Fall back to CPU if GPU fails
        batch_size: Batch size for processing
    
    Returns:
        Dict mapping sequence id to (prediction_tensor, probabilities_tensor)
    """
    import torch
    from tmbed.utils import Protein, make_batches, make_mask
    from tmbed.tmbed import encode_sequences, predict_sequences
    
    # Convert to Protein objects
    proteins = [Protein(header=f">{seq_id}", sequence=seq) for seq_id, seq in sequences]
    
    # Sort by length for efficient batching
    proteins_sorted = sorted(proteins, key=lambda p: p.length)
    
    # Create batches
    batches = make_batches(proteins_sorted, batch_size)
    
    predictions = {}
    
    for a, b in batches:
        batch = proteins_sorted[a:b]
        
        lengths = [protein.length for protein in batch]
        seqs = [protein.sequence for protein in batch]
        
        try:
            embeddings = encode_sequences(encoder, seqs, cpu_fallback)
        except RuntimeError as e:
            if not cpu_fallback:
                raise
            warnings.warn(f"TMbed encoding failed: {e}")
            continue
        
        embeddings = embeddings.to(device=device, dtype=torch.float32)
        
        mask = make_mask(embeddings, lengths)
        
        probabilities = predict_sequences(models, embeddings, mask)
        
        mask = mask.cpu()
        probabilities = probabilities.cpu()
        
        prediction = decoder(probabilities, mask).byte()
        
        # Permute probabilities to (batch, length, classes)
        probabilities = probabilities.permute(0, 2, 1)
        
        for idx, protein in enumerate(batch):
            length = protein.length
            seq_hash = protein.seq_hash
            predictions[seq_hash] = (
                prediction[idx, :length],
                probabilities[idx, :length]
            )
    
    return predictions


def parse_tmbed_predictions(
    seq_id: str,
    sequence: str,
    prediction_tensor,
    probabilities_tensor,
) -> List[FeatureResult]:
    """
    Parse TMbed predictions into FeatureResult objects.
    
    Merges consecutive residues with the same prediction into contiguous features.
    Only outputs transmembrane and signal peptide regions (not i/o regions).
    
    Args:
        seq_id: Sequence identifier
        sequence: Amino acid sequence
        prediction_tensor: TMbed prediction tensor (per-residue class labels)
        probabilities_tensor: TMbed probabilities tensor (per-residue class probabilities)
    
    Returns:
        List of FeatureResult objects for transmembrane/signal peptide regions
    """
    results = []
    
    predictions = prediction_tensor.tolist()
    probs = probabilities_tensor.tolist()
    
    # Merge consecutive residues with same label into regions
    if len(predictions) == 0:
        return results
    
    current_label = predictions[0]
    region_start = 0
    region_probs = [probs[0]]
    
    for i in range(1, len(predictions)):
        if predictions[i] == current_label:
            region_probs.append(probs[i])
        else:
            # Output the completed region (only for TM/signal peptide)
            if current_label in (0, 1, 2, 3, 4):  # B, b, H, h, S
                result = _create_feature_result(
                    current_label, region_start, i, region_probs
                )
                results.append(result)
            
            # Start new region
            current_label = predictions[i]
            region_start = i
            region_probs = [probs[i]]
    
    # Don't forget the last region
    if current_label in (0, 1, 2, 3, 4):
        result = _create_feature_result(
            current_label, region_start, len(predictions), region_probs
        )
        results.append(result)
    
    return results


def _create_feature_result(
    label: int,
    start: int,
    end: int,
    region_probs: List[List[float]],
) -> FeatureResult:
    """Create a FeatureResult from a region."""
    label_char, name, desc = TMBED_LABELS[label]
    
    # Calculate average probability for the predicted class
    # TMbed probabilities are ordered: [B, H, S, i, o] (5 classes)
    # Map label to probability index
    if label in (0, 1):  # B, b -> beta strand
        prob_idx = 0
    elif label in (2, 3):  # H, h -> alpha helix
        prob_idx = 1
    elif label == 4:  # S -> signal peptide
        prob_idx = 2
    else:
        prob_idx = 0  # fallback
    
    # Average probability for the region
    avg_prob = sum(p[prob_idx] for p in region_probs) / len(region_probs)
    score = avg_prob * 100  # Convert to percentage
    
    # Store all probabilities for the region
    prob_dict = {
        'P(B)': [p[0] for p in region_probs],
        'P(H)': [p[1] for p in region_probs],
        'P(S)': [p[2] for p in region_probs],
        'P(i)': [p[3] for p in region_probs],
        'P(o)': [p[4] for p in region_probs],
    }
    
    return FeatureResult(
        name=name,
        desc=desc,
        score=score,
        evalue=score_to_pseudo_evalue(score),
        start=start,
        end=end,
        program="TMbed",
        identity=score,
        probabilities=prob_dict,
    )


def add_feature_annotations(
    contig: SeqRecord,
    features: List[FeatureResult],
    is_protein: bool,
    cds_feature=None,
    cds_index: int = 0,
):
    """
    Add feature annotations to a contig.
    
    Args:
        contig: SeqRecord to annotate
        features: List of FeatureResult objects
        is_protein: Whether the contig is a protein sequence
        cds_feature: The CDS feature (for nucleotide sequences)
        cds_index: Index of the CDS feature (for nucleotide sequences)
    """
    for feat in features:
        if is_protein:
            location = FeatureLocation(feat.start, feat.end, strand=1)
            cds_id = f'0_1_{len(contig)}'
        else:
            # Map protein coordinates to nucleotide coordinates
            annot_length = (feat.end - feat.start) * 3
            try:
                location = cds_feature.location.overlay(feat.start * 3, annot_length)
            except Exception as e:
                warnings.warn(
                    f"Could not overlay location for {feat.name} on {contig.id}. "
                    f"Skipping annotation: {e}"
                )
                continue
            cds_id = cds_feature.qualifiers.get('cds_id', ['.'])[0]
        
        # Build qualifiers
        qualifiers = {
            'program': [feat.program],
            'database': [feat.program],  # For these algorithms, database = program
            'description': [feat.desc],
            'accession': ['.'],
            'evalue': [f"{feat.evalue:.1e}"],
            'score': [f"{feat.score:.1f}"],
            'name': [feat.name],
            'identity': [f"{feat.identity:.1f}"],
            'cds_id': [cds_id],
            'rstart': ['.'],
            'rend': ['.'],
            'rlen': ['.'],
        }
        
        f = SeqFeature(
            location=location,
            type=DOMAIN_FEATURE_NAME,
            qualifiers=qualifiers
        )
        contig.features.append(f)


def find_features(
    seq_iterator: Iterator[SeqRecord],
    algorithms: List[str] = ["all"],
    gpu_device: Optional[str] = "cuda:0",
    cpu: int = 0,
    batch_size: int = 4000,
    gene_call: Optional[str] = None,
    model_dir: Optional[str] = None,
    tmbed_model_dir: Optional[str] = None,
    motifs: Optional[List[str]] = None,
) -> Iterator[SeqRecord]:
    """
    Main function for feature detection.
    
    Yields annotated SeqRecords.
    
    Args:
        seq_iterator: Iterator yielding SeqRecord objects
        algorithms: List of algorithms to run (["all"] or specific algorithms)
        gpu_device: GPU device to use (None or "None" for CPU only)
        cpu: Number of CPU threads (0 = use all available)
        batch_size: Batch size for processing
        gene_call: Gene calling mode ("all", "unannotated", or None)
        model_dir: Base directory for all algorithm models (can be overridden per-algorithm)
        tmbed_model_dir: Directory for TMbed's ProtT5 model (overrides model_dir if set)
        motifs: List of PROSITE pattern strings to search for
    
    Yields:
        Annotated SeqRecord objects
    """
    # Determine which algorithms to run
    run_tmbed = "all" in algorithms or "TMbed" in algorithms or "tmbed" in algorithms
    
    if gene_call not in (None, 'all', 'unannotated'):
        raise ValueError("gene_call must be one of None, 'all', or 'unannotated'")
    
    clear_CDS_annotations = gene_call == 'all'
    
    # Initialize TMbed if needed
    tmbed_components = None
    get_md5 = None
    if run_tmbed:
        try:
            # Use algorithm-specific dir if set, otherwise fall back to generic model_dir
            effective_tmbed_model_dir = tmbed_model_dir if tmbed_model_dir else model_dir
            tmbed_components = get_tmbed_predictor(gpu_device, cpu, effective_tmbed_model_dir)
            from tmbed.utils import get_md5
        except ImportError as e:
            warnings.warn(str(e))
            run_tmbed = False
    
    # Process sequences in batches
    # We need to collect sequences to process in batches for efficiency
    batch_contigs = []
    batch_sequences = []  # List of (contig_idx, cds_idx, seq_hash, sequence) tuples
    
    # Simple MD5 fallback if tmbed not available
    if get_md5 is None:
        from hashlib import md5
        def get_md5(s):
            return md5(s.encode()).hexdigest()
    
    for rec in seq_iterator:
        CDS_count = clean_rec(rec, clear_CDS_annotations=clear_CDS_annotations)
        
        # Gene calling if needed
        if CDS_count == 0 and rec.annotations.get('molecule_type') != "protein" and gene_call is not None:
            prodigal_CDS_annotate(rec)
        
        contig_idx = len(batch_contigs)
        batch_contigs.append(rec)
        
        # Collect sequences for batch processing
        if rec.annotations.get('molecule_type') == "protein":
            seq = str(rec.seq)
            seq_hash = get_md5(seq)
            batch_sequences.append((contig_idx, -1, seq_hash, seq))  # -1 indicates protein contig
        else:
            for cds_idx, feature in enumerate(rec.features):
                if feature.type == 'CDS' and 'translation' in feature.qualifiers:
                    seq = feature.qualifiers['translation'][0]
                    seq_hash = get_md5(seq)
                    batch_sequences.append((contig_idx, cds_idx, seq_hash, seq))
        
        # Process batch if large enough
        if len(batch_sequences) >= batch_size:
            yield from _process_batch(
                batch_contigs, batch_sequences, run_tmbed, tmbed_components, batch_size, motifs
            )
            batch_contigs = []
            batch_sequences = []
    
    # Process remaining sequences
    if batch_contigs:  # Always process remaining contigs, even if no sequences
        yield from _process_batch(
            batch_contigs, batch_sequences, run_tmbed, tmbed_components, batch_size, motifs
        )


def _process_batch(
    batch_contigs: List[SeqRecord],
    batch_sequences: List[Tuple[int, int, str, str]],
    run_tmbed: bool,
    tmbed_components,
    batch_size: int,
    motifs: Optional[List[str]] = None,
) -> Iterator[SeqRecord]:
    """Process a batch of sequences and yield annotated contigs."""
    
    # Run motif search if patterns provided
    if motifs:
        # Group results by contig
        contig_motif_features: Dict[int, Dict[int, List[FeatureResult]]] = {}
        
        for contig_idx, cds_idx, seq_hash, seq in batch_sequences:
            features = search_all_motifs(seq, motifs)
            if features:
                if contig_idx not in contig_motif_features:
                    contig_motif_features[contig_idx] = {}
                contig_motif_features[contig_idx][cds_idx] = features
        
        # Add motif annotations to contigs
        for contig_idx, contig in enumerate(batch_contigs):
            if contig_idx in contig_motif_features:
                is_protein = contig.annotations.get('molecule_type') == "protein"
                
                for cds_idx, features in contig_motif_features[contig_idx].items():
                    if cds_idx == -1:  # Protein contig
                        add_feature_annotations(contig, features, is_protein=True)
                    else:  # Nucleotide contig with CDS
                        cds_feature = contig.features[cds_idx]
                        add_feature_annotations(
                            contig, features, is_protein=False,
                            cds_feature=cds_feature, cds_index=cds_idx
                        )
    
    if run_tmbed and tmbed_components is not None:
        encoder, models, decoder, device, use_gpu = tmbed_components
        
        # Prepare sequences for TMbed
        seq_data = [(seq_hash, seq) for _, _, seq_hash, seq in batch_sequences]
        
        # Run TMbed
        predictions = run_tmbed_batch(
            seq_data, encoder, models, decoder, device,
            cpu_fallback=True, batch_size=batch_size
        )
        
        # Group results by contig
        contig_features: Dict[int, Dict[int, List[FeatureResult]]] = {}
        
        for contig_idx, cds_idx, seq_hash, seq in batch_sequences:
            if seq_hash in predictions:
                pred_tensor, prob_tensor = predictions[seq_hash]
                features = parse_tmbed_predictions(seq_hash, seq, pred_tensor, prob_tensor)
                
                if contig_idx not in contig_features:
                    contig_features[contig_idx] = {}
                contig_features[contig_idx][cds_idx] = features
        
        # Add annotations to contigs
        for contig_idx, contig in enumerate(batch_contigs):
            if contig_idx in contig_features:
                is_protein = contig.annotations.get('molecule_type') == "protein"
                
                for cds_idx, features in contig_features[contig_idx].items():
                    if cds_idx == -1:  # Protein contig
                        add_feature_annotations(contig, features, is_protein=True)
                    else:  # Nucleotide contig with CDS
                        cds_feature = contig.features[cds_idx]
                        add_feature_annotations(
                            contig, features, is_protein=False,
                            cds_feature=cds_feature, cds_index=cds_idx
                        )
    
    # Yield all contigs
    for contig in batch_contigs:
        yield contig


def main(argv):
    parser = ArgumentParser(
        f"\nversion: {__version__}\n\n" + __doc__,
        formatter_class=RawAndDefaultsFormatter
    )
    
    parser.add_argument(
        '-i', '--input', default=None, required=True,
        nargs='+', type=str,
        help="The genbank or fasta files to annotate. Genbank files can be "
             "nucleotide (with CDS annotations) or peptide. Fasta files must be peptide."
    )
    
    parser.add_argument(
        "--fasta_type", type=str, default="protein",
        choices={"protein", "nucleotide"},
        help="Whether the sequences in fasta files are protein or nucleotide sequences."
    )
    
    parser.add_argument(
        '-o', '--output', default=None, type=str, required=False,
        help="Output genbank filename. If not supplied, writes to stdout."
    )
    
    parser.add_argument(
        '--cpu', type=int, default=0,
        help="The number of CPU threads to use [default: use all available cores]"
    )
    
    parser.add_argument(
        '--gene_call', type=str, default=None,
        choices={"all", "unannotated"}, required=False,
        help="When activated, new CDS annotations will be added with Prodigal in "
             "Metagenomic mode. If 'all', then any existing CDS annotations will be "
             "deleted and all contigs will be re-annotated. If 'unannotated', then "
             "only contigs without CDS annotations will be annotated. [default: None]. "
             "If you supply a nucleotide fasta file, be sure to also specify "
             "--fasta_type nucleotide."
    )
    
    parser.add_argument(
        "--gpu_device", default="cuda:0", type=str, required=False,
        help="For models that can use GPU acceleration, use this device if available. "
             "By default: if not set to 'None' try using the specified device for "
             "algorithms that allow GPU acceleration, if GPU not available, fall back "
             "to CPU mode if available. Choices: ['None', '0', 'cuda:0', ...]"
    )
    
    parser.add_argument(
        "--algorithms", default=["all"], nargs="+", type=str, required=False,
        help="Which algorithms to run. By default, run all algorithms. "
             "Choices: ['all', 'TMbed', 'none']. Use 'none' to run only motif "
             "search without ML algorithms. Note: MOTIF search is controlled "
             "separately via --motif."
    )
    
    parser.add_argument(
        '--motif', type=str, default=None, nargs='+',
        help="PROSITE motifs to search for. More information about PROSITE "
             "motif specification can be found: "
             "https://emboss.sourceforge.net/apps/release/6.5/emboss/apps/fuzzpro.html "
             "Example: '[DE](2)HS{P}X(2)PX(2,4)C' matches two Asp/Glu, then "
             "His, Ser, any except Pro, two of any, Pro, two to four of any, Cys. "
             "Motif names for boolean expressions are auto-generated by sanitizing "
             "the pattern (replacing parentheses and dashes with underscores)."
    )
    
    parser.add_argument(
        "--batch_size", type=int, default=4000, required=False,
        help="Batch size for sequence processing. Larger values may be faster "
             "but require more memory."
    )
    
    parser.add_argument(
        "--model_dir", type=str, default=None, required=False,
        help="Base directory for downloading/loading algorithm model weights. "
             "Individual algorithms will create subdirectories as needed. "
             "Can be overridden per-algorithm with --[algorithm]_model_dir."
    )
    
    parser.add_argument(
        "--tmbed_model_dir", type=str, default=None, required=False,
        help="Directory to download/load TMbed's ProtT5 encoder model (~2.25 GB). "
             "Overrides --model_dir for TMbed. "
             "If neither is specified, uses HF_HOME env var if set, otherwise "
             "downloads to tmbed/models/t5/ inside the installed package."
    )
    
    parser.add_argument('--config', action=ActionConfigFile)
    
    params = parser.parse_args(argv)
    
    # Validate algorithms
    # 'none' allows running only motif search without any ML algorithms
    valid_algorithms = {'all', 'TMbed', 'tmbed', 'none'}
    for alg in params.algorithms:
        if alg not in valid_algorithms:
            raise ValueError(f"Unknown algorithm: {alg}. Valid choices: {valid_algorithms}")
    
    # Handle gpu_device="None" string
    gpu_device = params.gpu_device
    if gpu_device and gpu_device.lower() == "none":
        gpu_device = None
    
    # Setup output
    if params.output is None:
        out = sys.stdout
    else:
        out = open(params.output, "w")
    
    # Run
    write_genbank(
        find_features(
            parse_seqfiles(
                params.input,
                None,
                filetype_override=None,
                default_molecule_type=params.fasta_type
            ),
            algorithms=params.algorithms,
            gpu_device=gpu_device,
            cpu=params.cpu,
            batch_size=params.batch_size,
            gene_call=params.gene_call,
            model_dir=params.model_dir,
            tmbed_model_dir=params.tmbed_model_dir,
            motifs=params.motif,
        ),
        out
    )
    
    if params.output is not None:
        out.close()


def _entrypoint():
    main(sys.argv[1:])


if __name__ == '__main__':
    main(sys.argv[1:])
