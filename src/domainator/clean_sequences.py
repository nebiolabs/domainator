"""Clean and filter sequences in genbank or fasta files

This module provides functionality to:
- Clean sequence names for compatibility with newick format
- Deduplicate sequence names to ensure uniqueness
- Strip non-canonical residues from sequence ends
- Filter out sequences containing non-canonical residues

"""
from dataclasses import dataclass
from jsonargparse import ActionConfigFile
import sys

from domainator.cli.common import add_alphabet, add_input, add_seq_outputs, check_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from domainator.utils import parse_seqfiles, slice_record
from domainator.Bio.Data import IUPACData


@dataclass
class CleaningStats:
    """Track statistics about sequence cleaning operations."""
    input_count: int = 0
    output_count: int = 0
    filtered_count: int = 0
    empty_after_strip_count: int = 0
    names_changed_count: int = 0
    sequences_changed_count: int = 0

# Canonical amino acids (20 standard amino acids)
CANONICAL_PROTEIN = set(IUPACData.protein_letters.upper())
# Canonical nucleotides (GATC for DNA, GAUC for RNA)
CANONICAL_NUCLEOTIDE = set(IUPACData.unambiguous_dna_letters.upper() + IUPACData.unambiguous_rna_letters.upper())


def clean_name_for_newick(name):
    """
    Converts problematic characters for newick format to underscores.
    
    Characters converted: space, semicolon, colon, comma, parentheses, quotes
    
    Args:
        name: The sequence name to clean
        
    Returns:
        Cleaned name string
    """
    bad_chars = " ;:,()'\""
    chars = ["_" if x in bad_chars else x for x in name]
    return "".join(chars)


def deduplicate_name(name, seen_names):
    """
    Ensure a name is unique by appending a suffix if necessary.
    
    Args:
        name: The sequence name
        seen_names: A dict tracking seen names and their counts
        
    Returns:
        A unique version of the name
    """
    candidate = name
    while True:
        if candidate not in seen_names:
            seen_names[candidate] = 1
            return candidate

        count = seen_names[candidate]
        seen_names[candidate] = count + 1
        candidate = f"{candidate}_{count}"


def strip_non_canonical(sequence, canonical_chars):
    """
    Strip non-canonical characters from both ends of a sequence.
    Also strips * characters (stop codons).
    
    Args:
        sequence: The sequence string
        canonical_chars: Set of canonical characters to keep
        
    Returns:
        Stripped sequence string
    """
    start, end = strip_non_canonical_bounds(sequence, canonical_chars)
    return str(sequence)[start:end]


def strip_non_canonical_bounds(sequence, canonical_chars):
    """Return the slice bounds after trimming non-canonical characters from both ends."""
    seq_upper = str(sequence).upper()

    start = 0
    while start < len(seq_upper) and (seq_upper[start] not in canonical_chars or seq_upper[start] == '*'):
        start += 1

    end = len(seq_upper)
    while end > start and (seq_upper[end - 1] not in canonical_chars or seq_upper[end - 1] == '*'):
        end -= 1

    return start, end


def has_non_canonical(sequence, canonical_chars):
    """
    Check if a sequence contains any non-canonical characters.
    Also considers * as non-canonical.
    
    Args:
        sequence: The sequence string
        canonical_chars: Set of canonical characters
        
    Returns:
        True if the sequence contains non-canonical characters, False otherwise
    """
    seq_upper = str(sequence).upper()
    for char in seq_upper:
        if char not in canonical_chars or char == '*':
            return True
    return False


def clean_sequences(contigs, fasta_type="protein", clean_name=False, deduplicate_names=False, 
                    strip=False, filter_by_aa=False, stats=None):
    """
    Clean and filter sequences based on specified criteria.
    
    Args:
        contigs: An iterable/iterator of SeqRecords
        fasta_type: "protein" or "nucleotide" - determines which canonical character set to use
        clean_name: If True, make names compatible with newick format
        deduplicate_names: If True, ensure all sequence names are unique
        strip: If True, strip non-canonical residues from sequence ends
        filter_by_aa: If True, filter out sequences containing non-canonical residues
        stats: Optional CleaningStats object to track statistics
        
    Yields:
        Cleaned/filtered SeqRecords
    """
    if fasta_type == "protein":
        canonical_chars = CANONICAL_PROTEIN
    else:
        canonical_chars = CANONICAL_NUCLEOTIDE
    
    seen_names = dict()
    
    for contig in contigs:
        if stats is not None:
            stats.input_count += 1
        
        original_id = contig.id
        original_seq = str(contig.seq)
        
        # Strip non-canonical characters from ends if requested
        if strip:
            start, end = strip_non_canonical_bounds(contig.seq, canonical_chars)
            if start != 0 or end != len(contig.seq):
                contig = slice_record(contig, start, end)
        
        # Filter out sequences with non-canonical characters if requested
        if filter_by_aa:
            if has_non_canonical(contig.seq, canonical_chars):
                if stats is not None:
                    stats.filtered_count += 1
                continue
        
        # Skip empty sequences
        if len(contig.seq) == 0:
            if stats is not None:
                stats.empty_after_strip_count += 1
            continue
        
        # Clean name for newick compatibility if requested
        if clean_name:
            contig.id = clean_name_for_newick(contig.id)
            contig.name = clean_name_for_newick(contig.name)
        
        # Deduplicate names if requested
        if deduplicate_names:
            new_id = deduplicate_name(contig.id, seen_names)
            if new_id != contig.id:
                contig.id = new_id
                contig.name = new_id
        
        # Track statistics
        if stats is not None:
            if contig.id != original_id:
                stats.names_changed_count += 1
            if str(contig.seq) != original_seq:
                stats.sequences_changed_count += 1
            stats.output_count += 1
        
        yield contig


def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, 'Names of input genbank files.')

    add_alphabet(parser)

    add_seq_outputs(parser)

    parser.add_argument('--clean_name', action='store_true', default=False, 
                        help="if true then make the name compatible with newick format.")

    parser.add_argument('--deduplicate_names', action='store_true', default=False, 
                        help="if true then rename sequences such that output sequences have unique names.")

    parser.add_argument('--strip', action='store_true', default=False, 
                        help="if true then strip non-canonical residues or nucleotides from the ends of sequences (depending on input type). Including * characters.")

    parser.add_argument('--filter_by_aa', action='store_true', default=False, 
                        help="if true then throw out sequences containing non-canonical residues or nucleotides (depending on input type). Including * characters. (applied after strip)")

    parser.add_argument('--log_out', default=None, required=False,
                        help="Write the cleaning statistics to this file. If not supplied, writes them to stderr.")

    parser.add_argument('--config', action=ActionConfigFile)
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)

    check_outputs(parser, params, ["genbank_out", "fasta_out"])
    log_handle = sys.stderr if params.log_out is None else open(params.log_out, "w")

    # Initialize statistics tracking
    stats = CleaningStats()

    # Run
    cleaned_contigs_iterator = clean_sequences(
        parse_seqfiles(params.input, default_molecule_type=params.alphabet),
        fasta_type=params.alphabet,
        clean_name=params.clean_name,
        deduplicate_names=params.deduplicate_names,
        strip=params.strip,
        filter_by_aa=params.filter_by_aa,
        stats=stats
    )

    write_seq_outputs(cleaned_contigs_iterator, genbank=params.genbank_out, fasta=params.fasta_out)

    # Write log
    print(f"Input sequences:            {stats.input_count}", file=log_handle)
    print(f"Output sequences:           {stats.output_count}", file=log_handle)
    print(f"Filtered (non-canonical):   {stats.filtered_count}", file=log_handle)
    print(f"Filtered (empty after strip): {stats.empty_after_strip_count}", file=log_handle)
    print(f"Names changed:              {stats.names_changed_count}", file=log_handle)
    print(f"Sequences changed:          {stats.sequences_changed_count}", file=log_handle)

    if params.log_out is not None:
        log_handle.close()
