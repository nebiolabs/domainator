"""Partial (fragment) status of contigs, CDSs and proteins.
"""
from typing import Optional, Tuple, Union
import warnings

from domainator import lean_record
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, BeforePosition, AfterPosition
from Bio.SeqRecord import SeqRecord
from domainator.model import get_cds_unique_name


# UniProt FASTA headers flag partial proteins with one of these tokens.
FRAGMENT_TOKENS = (" (Fragment)", " (Fragments)")
# UniProt text (.dat) entries flag them with a "Flags:" DE line, e.g. "Flags: Fragment;" or "Flags: Precursor; Fragments;",
# which the swiss parser appends to the description (and so to the DEFINITION of a converted GenBank file).
FRAGMENT_FLAGS_PREFIX = "Flags:"
# MGnify protein FASTA headers say whether a protein is full length with an FL= token: FL=0 is a fragment, FL=1 full length.
MGNIFY_FRAGMENT_TOKEN = "FL=0"
# UniProt marks a missing terminus with a NON_TER feature on the first or last residue.
NON_TERMINAL_FEATURE_NAME = "NON_TER"
PARTIAL_CHOICES = ("include", "exclude", "only")

def location_partial_ends(location) -> Tuple[bool, bool]:
    """
        Returns (five_prime_partial, three_prime_partial) for a location, in transcript orientation.
        A '<' or '>' marker on the outer end of the first part (5') or last part (3') makes that end partial.
    """
    first = location.parts[0]
    last = location.parts[-1]
    if first.strand == -1:
        five_prime = isinstance(first.end, AfterPosition)
    else:
        five_prime = isinstance(first.start, BeforePosition)
    if last.strand == -1:
        three_prime = isinstance(last.start, BeforePosition)
    else:
        three_prime = isinstance(last.end, AfterPosition)
    return five_prime, three_prime

def feature_partial_ends(feature:Union[SeqFeature, "lean_record.LeanFeature"]) -> Tuple[bool, bool]:
    """
        location_partial_ends for a SeqFeature or a LeanFeature.
    """
    if isinstance(feature, lean_record.LeanFeature):
        return feature.partial_ends
    return location_partial_ends(feature.location)

def normalize_partial_codes(record) -> None:
    """
        Converts prodigal-style /partial="10" codes on a record's features to '<' / '>' location markers, and removes the codes.
        See lean_record.apply_partial_code. Called by parse_seqfiles on every GenBank record.
    """
    if isinstance(record, lean_record.LeanContig):
        for feature in record.features:
            lean_record.normalize_lean_partial_code(feature, len(record))
        return
    for feature in record.features:
        ends = lean_record.partial_code_ends(feature.qualifiers)
        if ends is None:
            continue
        del feature.qualifiers[lean_record.PARTIAL_QUALIFIER]
        if any(location_partial_ends(feature.location)) or ends == (False, False):
            continue
        left_partial, right_partial = ends
        parts = list(feature.location.parts)
        leftmost = min(range(len(parts)), key=lambda i: int(parts[i].start))
        rightmost = max(range(len(parts)), key=lambda i: int(parts[i].end))
        left_extension = right_extension = 0
        if left_partial:
            part = parts[leftmost]
            start, _end, left_extension, _r = lean_record.apply_partial_code(int(part.start), int(part.end), True, False, len(record))
            parts[leftmost] = FeatureLocation(BeforePosition(start), part.end, strand=part.strand)
        if right_partial:
            part = parts[rightmost]
            _start, end, _l, right_extension = lean_record.apply_partial_code(int(part.start), int(part.end), False, True, len(record))
            parts[rightmost] = FeatureLocation(part.start, AfterPosition(end), strand=part.strand)
        # the 5' end is the start of the first part on the plus strand, the end of the first part on the minus strand
        if parts[0].strand == -1:
            five_prime_extension = right_extension if rightmost == 0 else 0
        else:
            five_prime_extension = left_extension if leftmost == 0 else 0
        feature.location = parts[0] if len(parts) == 1 else CompoundLocation(parts, operator=feature.location.operator)
        if feature.type == "CDS":
            codon_start = lean_record.new_codon_start(feature.qualifiers, five_prime_extension)
            if codon_start is not None:
                feature.qualifiers["codon_start"] = [codon_start]

def partial_ends_to_status(five_prime:bool, three_prime:bool) -> Optional[str]:
    """
        Converts partial end flags to a fragment status: None (complete), "N", "C", or "NC".
    """
    status = ("N" if five_prime else "") + ("C" if three_prime else "")
    return status if status else None

def get_fragment_status(record_or_feature:Union[SeqRecord, SeqFeature]) -> Optional[str]:
    """
        Returns the fragment status of a CDS feature or a protein record:
            None: complete (or nucleotide record, whose fragment status is per-CDS)
            "N": missing the N-terminus (5' partial)
            "C": missing the C-terminus (3' partial)
            "NC": missing both ends
            "?": a fragment with unknown missing ends (UniProt " (Fragment)" or MGnify "FL=0" header, or a CDS with a /partial qualifier without a value)

        For a protein record, a CDS or Protein feature spanning the record (e.g. from domain_search --translate, or GenPept)
        decides. Otherwise UniProt NON_TER features on the first or last residue (UniProt text entries) decide,
        and otherwise a UniProt fragment token in the description (see description_is_fragment) gives "?".
    """
    if isinstance(record_or_feature, (SeqFeature, lean_record.LeanFeature)):
        status = partial_ends_to_status(*feature_partial_ends(record_or_feature))
        if status is None and lean_record.has_bare_partial(record_or_feature.qualifiers):
            return "?" # legacy /partial qualifier with no value: partial, missing end unknown
        return status
    record = record_or_feature
    if record.annotations.get("molecule_type") != "protein":
        return None
    for feature in record.features:
        if feature.type in ("CDS", "Protein"):
            if isinstance(feature, lean_record.LeanFeature):
                start, end = feature.start, feature.end
            else:
                start, end = int(feature.location.start), int(feature.location.end)
            if start == 0 and end >= len(record):
                return partial_ends_to_status(*feature_partial_ends(feature))
    missing_n_terminus = False
    missing_c_terminus = False
    for feature in record.features:
        if feature.type == NON_TERMINAL_FEATURE_NAME:
            if isinstance(feature, lean_record.LeanFeature):
                start, end = feature.start, feature.end
            else:
                start, end = int(feature.location.start), int(feature.location.end)
            missing_n_terminus = missing_n_terminus or start == 0
            missing_c_terminus = missing_c_terminus or end >= len(record)
    if missing_n_terminus or missing_c_terminus:
        return partial_ends_to_status(missing_n_terminus, missing_c_terminus)
    if description_is_fragment(record.description):
        return "?"
    return None

def description_is_fragment(description:str) -> bool:
    """
        True if a description carries a fragment token: UniProt " (Fragment)" or " (Fragments)" (FASTA headers),
        a UniProt "Flags:" list containing Fragment or Fragments (text entries), or an MGnify FL=0 token (FASTA headers).
    """
    if any(token in description for token in FRAGMENT_TOKENS):
        return True
    if MGNIFY_FRAGMENT_TOKEN in description.split():
        return True
    flags_start = description.find(FRAGMENT_FLAGS_PREFIX)
    return flags_start != -1 and "Fragment" in description[flags_start:]

def contig_has_fragment(record) -> bool:
    """
        Contig-level fragment test used by contig selection and reporting:
        a protein record is a fragment if get_fragment_status says so; a nucleotide contig if it contains any partial CDS.
    """
    if record.annotations.get("molecule_type") == "protein":
        return get_fragment_status(record) is not None
    return any(feature.type == "CDS" and get_fragment_status(feature) is not None for feature in record.features)

def cds_fragment_status(record, feature:Optional[SeqFeature]) -> Optional[str]:
    """
        The fragment status of a CDS on a record: the record's status if it is a protein record, otherwise the CDS feature's.
    """
    if record.annotations.get("molecule_type") == "protein" or feature is None:
        return get_fragment_status(record)
    return get_fragment_status(feature)

def filter_by_partial(records, partial:str):
    """
        Yields the records that pass a --partial filter at the contig level (see contig_has_fragment).
    """
    for record in records:
        if partial == "include" or fragment_status_allowed("?" if contig_has_fragment(record) else None, partial):
            yield record

def fragment_status_allowed(status:Optional[str], partial:str) -> bool:
    """
        Whether a protein/CDS with the given fragment status passes a --partial filter ("include", "exclude", or "only").
    """
    if partial == "include":
        return True
    elif partial == "exclude":
        return status is None
    elif partial == "only":
        return status is not None
    raise ValueError(f"partial must be one of {', '.join(PARTIAL_CHOICES)}, not {partial}")

def codon_offset(feature:SeqFeature) -> int:
    """
        Returns the 0-based reading frame offset of a CDS from its codon_start qualifier (0 if absent or invalid).
    """
    try:
        offset = int(feature.qualifiers.get("codon_start", ["1"])[0]) - 1
    except ValueError:
        return 0
    return offset if 0 <= offset <= 2 else 0

def on_contig_translation(translation:str, location_length:int, offset:int, five_prime_partial:bool, three_prime_partial:bool) -> Optional[str]:
    """
        Restricts a partial CDS's translation to the residues encoded on the contig.

        Some files give partial CDSs the translation of the full protein, including residues off the edge of the contig.
        If the translation is longer than the on-contig coding sequence allows, the excess is removed from the partial end.
        A complete 3' end ends in a stop codon, which GenBank translations omit, but Domainator's own translations include as '*'.

        Returns:
            the on-contig translation, or None if the excess can't be placed (both ends partial).
    """
    expected = (location_length - offset) // 3
    if not three_prime_partial and not translation.endswith("*"):
        expected -= 1 # the stop codon is not in the translation
    excess = len(translation) - expected
    if excess <= 0 or not (five_prime_partial or three_prime_partial):
        return translation
    if five_prime_partial and not three_prime_partial:
        return translation[excess:]
    if three_prime_partial and not five_prime_partial:
        return translation[:expected]
    return None

def normalize_cds_translation(feature:SeqFeature, contig_seq) -> None:
    """
        Sets feature's translation qualifier to the on-contig translation: translated from contig_seq if absent or empty,
        and trimmed if it includes residues off the edge of the contig (see on_contig_translation).
    """
    translation = feature.qualifiers.get("translation", [""])[0]
    if len(translation) != 0:
        five_prime, three_prime = location_partial_ends(feature.location)
        translation = on_contig_translation(translation, len(feature.location), codon_offset(feature), five_prime, three_prime)
        if translation is None:
            warnings.warn(f"Translation of partial CDS {get_cds_unique_name(feature)} is longer than its location, and both ends are partial. Re-translating from the contig.")
            translation = ""
    if len(translation) == 0:
        translation = str(feature.translate(contig_seq, cds=False)) # cds = False, because we don't want to throw exceptions on weird annotations.
    feature.qualifiers["translation"] = [translation]
