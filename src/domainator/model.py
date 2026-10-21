"""CDS naming and the per-CDS view of Domainator annotations.
"""
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Union
import hashlib
import warnings

from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
from domainator.location_ops import copy_location, stranded_end_human_readable, stranded_start_human_readable


def feature_evalue(feature) -> float:
    """The e-value of a Domainator or Domain_Search annotation."""
    return float(feature.qualifiers["evalue"][0])


def feature_score(feature) -> float:
    """The bit score of a Domainator or Domain_Search annotation."""
    return float(feature.qualifiers["score"][0])


def evalue_passes(feature, evalue: Optional[float]) -> bool:
    """True if the annotation's e-value is at most ``evalue`` (None: no cutoff). Every ``-e/--evalue`` filter uses this."""
    return evalue is None or feature_evalue(feature) <= evalue


# Record annotations that domainate, domain_search and select_by_cds set on the records they yield.
DIST_TO_START_ANNOTATION = "dist_to_start"
"""Distance from the start of a region to the start of its focus feature."""
DIST_TO_END_ANNOTATION = "dist_to_end"
"""Distance from the end of the focus feature to the end of the region."""
FOCUS_STRAND_ANNOTATION = "focus_strand"
"""Strand of the focus feature in region coordinates (1 or -1)."""
HIT_SCORES_ANNOTATION = "hit_scores"
"""{CDS index: best domain score} of a searched contig."""
HIT_NAMES_ANNOTATION = "hit_names"
"""{CDS index: best domain name} of a searched contig."""
HIT_BEST_SCORE_ANNOTATION = "hit_best_score"
"""The best domain score of a searched contig or region."""
SEARCH_ANNOTATIONS = frozenset({DIST_TO_START_ANNOTATION, DIST_TO_END_ANNOTATION, FOCUS_STRAND_ANNOTATION,
                                HIT_SCORES_ANNOTATION, HIT_NAMES_ANNOTATION, HIT_BEST_SCORE_ANNOTATION})


def _get_annotation(record, key, what):
    try:
        return record.annotations[key]
    except KeyError:
        raise KeyError(f"{what} not set for {record.id}") from None


def set_dist_to_start(record, dist: int) -> None:
    record.annotations[DIST_TO_START_ANNOTATION] = dist


def get_dist_to_start(record) -> int:
    return _get_annotation(record, DIST_TO_START_ANNOTATION, "dist to start")


def set_dist_to_end(record, dist: int) -> None:
    record.annotations[DIST_TO_END_ANNOTATION] = dist


def get_dist_to_end(record) -> int:
    return _get_annotation(record, DIST_TO_END_ANNOTATION, "dist to end")


def set_focus_strand(record, strand: int) -> None:
    record.annotations[FOCUS_STRAND_ANNOTATION] = strand


def get_focus_strand(record) -> int:
    return _get_annotation(record, FOCUS_STRAND_ANNOTATION, "focus strand")


def set_hit_names(record, hit_names: Dict[int, str]) -> None:
    record.annotations[HIT_NAMES_ANNOTATION] = hit_names


def get_hit_names(record) -> Dict[int, str]:
    return _get_annotation(record, HIT_NAMES_ANNOTATION, "hit names")


def set_hit_scores(record, hit_scores: Dict[int, float]) -> None:
    record.annotations[HIT_SCORES_ANNOTATION] = hit_scores


def get_hit_scores(record) -> Dict[int, float]:
    return _get_annotation(record, HIT_SCORES_ANNOTATION, "hit scores")


def set_hit_best_score(record, score: float) -> None:
    record.annotations[HIT_BEST_SCORE_ANNOTATION] = score


def get_hit_best_score(record) -> float:
    return _get_annotation(record, HIT_BEST_SCORE_ANNOTATION, "hit best score")


ANNOTATION_TYPES = ("domainator", "domain_search", "all")
"""Which annotations a function reads: ``domainator`` is the domain annotations written by
domainate (``Domainator`` features), ``domain_search`` the best hits written by domain_search
(``Domain_Search`` features), and ``all`` both."""


def check_annotation_type(annotation_type: str, allowed=ANNOTATION_TYPES[:2]) -> str:
    """Return ``annotation_type``, or raise ValueError if it isn't one of ``allowed``."""
    if annotation_type not in allowed:
        raise ValueError(f"annotation_type must be one of {', '.join(allowed)}, not {annotation_type!r}.")
    return annotation_type


def get_cds_unique_name(feature):
    """
        If the feature already has a cds_id, then keep it, otherwise generate one based on the position on the contig.
    """
    if "cds_id" in feature.qualifiers:
        return feature.qualifiers["cds_id"][0]
    else:
        # need the strand information to account for circular contigs.
        name_parts = ["_".join( (str(stranded_start_human_readable(p)), str(p.strand), str(stranded_end_human_readable(p))) ) for p in feature.location.parts]
        return " ".join(name_parts) # space so it can be split into multiple lines when writing genbank files

def get_cds_name(feature, precedence=None): #(contig_id, feature):
    if precedence is None:
        precedence = ["gene_id", "locus_tag", "gene", "protein_id"]

    for p in precedence:
        if p in feature.qualifiers:
            return feature.qualifiers[p][0]
    return get_cds_unique_name(feature)
    

@dataclass
class DomainatorCDS():
    name: str
    num: str
    index: int #feature index from the full contig
    feature: SeqFeature
    domain_features: List[SeqFeature] = field(default_factory=lambda: [])
    domain_names: Set[str] = field(default_factory=lambda: set())
    domain_search_feature: SeqFeature = None
    is_nucleic_acid:bool = False

    @classmethod
    def list_from_contig(cls, contig, domain_evalue=float("inf"), domain_score=float("-inf"), skip_pseudo=False, name_precedence=None, include_nucleic_acid_annotations=False):
        """
            returns a list of DomainatorCDS objects from a contig SeqRecord
        """
        cdss: List[DomainatorCDS] = list()
        name_to_idx: Dict[str,int] = dict()
        idx:int = 0
        nuc_idx:int = 0 # nucleotide annotation index, used for naming nucleotide annotations that don't have cds_id qualifiers based on their position on the contig.

        def location_key(feature: SeqFeature) -> str:
            name_parts = ["_".join((str(stranded_start_human_readable(p)), str(p.strand), str(stranded_end_human_readable(p)))) for p in feature.location.parts]
            return " ".join(name_parts)

        for i, feature in enumerate(contig.features): #get the CDS features
            if feature.type == 'CDS':
                if skip_pseudo and ("pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers):
                    continue
                cds_num = get_cds_unique_name(feature)
                name = get_cds_name(feature, name_precedence)
                if cds_num not in name_to_idx:
                    #cds_rec = get_seq_record(start, end, name, strand, contig.seq)
                    cdss.append(DomainatorCDS(name, cds_num, i, feature))
                    name_to_idx[cds_num] = idx
                    idx += 1
                else:
                    warnings.warn(f"multiple cdss found with contig id: {contig.id}, CDS num: {cds_num}, so only the first one is kept.")
        for feature in contig.features: #get the domainator features
            if feature.type == DOMAIN_FEATURE_NAME or feature.type == DOMAIN_SEARCH_BEST_HIT_NAME:
                if feature_evalue(feature) <= domain_evalue and feature_score(feature) > domain_score:
                    if include_nucleic_acid_annotations and feature.qualifiers.get("cds_id", [None])[0] == ".":
                        cds_num = f"nuc_{nuc_idx}"
                        nuc_idx += 1
                        if cds_num not in name_to_idx: # should always be true.
                            name = feature.qualifiers.get("name", [cds_num])[0]
                            pseudo_feature = SeqFeature(location=feature.location, type="CDS", qualifiers={"cds_id": [cds_num]})
                            cdss.append(DomainatorCDS(name, cds_num, contig.features.index(feature), pseudo_feature))
                            cdss[-1].is_nucleic_acid = True
                            name_to_idx[cds_num] = idx
                            idx += 1
                    else:
                        cds_num = feature.qualifiers['cds_id'][0]
                        if cds_num not in name_to_idx:
                            if cds_num != ".": # We don't need to warn about missing cds_id qualifiers for features that are associated with nucleic acid annotations.
                                warnings.warn(f"Domainator annotation found with no associated CDS: {cds_num}")
                            continue
                    
                    if feature.type == DOMAIN_FEATURE_NAME:
                        cdss[name_to_idx[cds_num]].domain_features.append(feature)
                        cdss[name_to_idx[cds_num]].domain_names.add(feature.qualifiers["name"][0])
                    elif feature.type == DOMAIN_SEARCH_BEST_HIT_NAME:
                        cdss[name_to_idx[cds_num]].domain_search_feature = feature
        return cdss
    
    def get_domain_names(self, databases:Optional[Set[str]]=None, evalue:Optional[float]=None):

        out = set()
        for f in self.domain_features:
            keep = True
            if databases is not None and f.qualifiers["database"][0] not in databases:
                keep = False
            if not evalue_passes(f, evalue):
                keep = False

            if keep == True:    
                out.add(f.qualifiers["name"][0])
        return out


def get_domain_cds_annotation_strings(seq_record):
    """
        returns a dict of dicts of strings:
            {cds_id:{domain_db: notes_string}}
    """
    out = dict()
    for feature in seq_record.features:
        try:
            if feature.type == DOMAIN_FEATURE_NAME:
                cds_id = feature.qualifiers['cds_id'][0]
                db = feature.qualifiers['database'][0]
                name = feature.qualifiers['name'][0]
                eval = feature.qualifiers['evalue'][0]
                score = feature.qualifiers['score'][0]
                description = feature.qualifiers['description'][0]
                if cds_id not in out:
                    out[cds_id] = defaultdict(list)
                out[cds_id][db].append(f"{name} ({description}, {eval}, {score})")
        except KeyError as e:
            raise KeyError(f"Domainator feature on record {seq_record.id} is missing the {e} qualifier") from e
    return out

def region_id(parent:SeqRecord, location:Union[FeatureLocation,CompoundLocation], flip:bool=False) -> str:
    """The id of a region of parent, extracted with slice_record_from_location(parent, location), then reverse complemented if flip.

        All tools that extract regions of contigs name them this way:
            {parent.id}_{first}:{last}, followed by 'rc' if the region is reverse complemented relative to the parent.
            first and last are the 1-based coordinates on the parent of the first and last bases of the region, as written.
            So first > last means either that the region is reverse complemented, or that it wraps around the origin of a circular contig.
        If the region is the whole contig, starting at the original origin, then the coordinates are omitted:
            {parent.id} for the forward strand, and {parent.id}_rc for the reverse strand.
        Protein regions are never reverse complemented.

        Args:
            parent: the record the region was extracted from
            location: the location passed to slice_record_from_location
            flip: True if the extracted region was then reverse complemented (e.g. to keep the original direction of a reverse strand location)

        Returns:
            str: the id of the region
    """
    first = stranded_start_human_readable(location)
    last = stranded_end_human_readable(location)
    reverse = parent.annotations.get("molecule_type") != "protein" and location.parts[0].strand == -1
    if flip:
        first, last = last, first
        reverse = not reverse
    contig_length = len(parent)
    if len(location) == contig_length and (first, last) == ((contig_length, 1) if reverse else (1, contig_length)):
        return f"{parent.id}_rc" if reverse else parent.id
    return f"{parent.id}_{first}:{last}" + ("rc" if reverse else "")

def get_sources(rec:SeqRecord):
    """
        returns a list of all source features in a SeqRecord
    """
    sources = []
    for f in rec.features:
        if f.type == "source":
            sources.append(f)
    return sources

def get_non_domainator_features(rec:SeqRecord):
    features = []
    for f in rec.features:
        if f.type != DOMAIN_FEATURE_NAME and f.type != DOMAIN_SEARCH_BEST_HIT_NAME: #and f.type != "CDS":
            features.append(f)
    return features


def copy_feature(feature:SeqFeature, location:Union[FeatureLocation,CompoundLocation]=None) -> SeqFeature:
    """
        Creates a copy of a SeqFeature object, optionally replacing the location with a provided location,

    Args:
        feature (SeqFeature): the SeqFeature to copy
        location (FeatureLocation, optional): Used as the location. If None, a copy of the existing FeatureLocation will be used. Defaults to None.

    Returns:
        SeqFeature: a new SeqFeature object.
    """

    #   https://www.insdc.org/submitting-standards/feature-table/
    if location is None:
        location = copy_location(feature.location)
    return SeqFeature(
                location=location,
                type=feature.type,
                id=feature.id,
                qualifiers=feature.qualifiers.copy(),
            )

def sequence_hashes(seq, both_strands: bool = False) -> List[bytes]:
    """
        SHA-256 digests identifying a sequence for exact-duplicate detection.

        The sequence is uppercased and terminal '*' stripped before hashing, so case and stop codons don't distinguish sequences.

        seq: a Seq or str.
        both_strands: if True, also return the digests of the complement and reverse complement (for nucleotide sequences),
            so a sequence matches its other-strand copies.

        returns: a list of digests (1, or 3 if both_strands). Two sequences are duplicates if their lists share any digest.
    """
    seq_str = str(seq).upper().strip("*")
    out = [hashlib.sha256(seq_str.encode("utf-8")).digest()]
    if both_strands:
        seq_obj = Seq(seq_str)
        out.append(hashlib.sha256(str(seq_obj.complement()).encode("utf-8")).digest())
        out.append(hashlib.sha256(str(seq_obj.reverse_complement()).encode("utf-8")).digest())
    return out
