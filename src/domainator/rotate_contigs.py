"""Rotate circular contigs so that their origin is set relative to a domain annotation.

For each contig, one domain is chosen from new annotations found by searching the references (-r) against the contig,
and/or from pre-existing Domainator (--domains) or domain_search (--annotation_type domain_search) annotations.
If there are several candidates, the one with the highest score is chosen (ties broken by lowest E-value, then by position).
The contig is then rotated so that the origin is at the start, end, or center of the chosen domain (--anchor),
and reverse complemented if necessary to put the chosen domain on the forward strand (unless --keep_direction).

Positions are relative to the chosen domain, in its own orientation:
    start:  the first base of the domain becomes the first base of the contig.
    end:    the last base of the domain becomes the last base of the contig.
    center: the origin is placed in the middle of the domain.

Inputs must be nucleotide sequences. Only circular contigs are rotated. FASTA records and linear GenBank records are
skipped with a warning, unless --circularize is set, in which case they are treated as circular.
Contigs without a qualifying domain are skipped with a warning. Output contigs are always marked as circular.

Protein references are searched against CDS translations, so contigs without CDS annotations (such as all FASTA inputs)
need --gene_call to be found by protein references. Nucleotide references (nucleotide HMMs, nucleotide FASTA, and
infernal CM files) are searched against the whole contig.
"""

import warnings
from itertools import islice
from typing import Iterable, Iterator, List, Optional, Set, Tuple, Union

from jsonargparse import ActionConfigFile
from jsonargparse.typing import PositiveFloat

from domainator.cli.common import STDIN_PATHS, add_annotation_type, add_input, add_seq_outputs, check_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator import utils
from domainator.Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from domainator.Bio.SeqRecord import SeqRecord, DOMAINATOR_DIST_TO_START_ANNOTATION, DOMAINATOR_DIST_TO_END_ANNOTATION, DOMAINATOR_FOCUS_STRAND_ANNOTATION, DOMAINATOR_HIT_SCORES_ANNOTATION, DOMAINATOR_HIT_NAMES_ANNOTATION, DOMAINATOR_HIT_BEST_SCORE_ANNOTATION
from domainator.domainate import domainate, read_references
from domainator.lean_record import LeanContig, lean_to_seqrecord, LEAN_SEARCH_TYPES, materialize_lean_search
from domainator.utils import parse_seqfiles, list_and_file_to_dict_keys, slice_record_from_location, region_id, WHERE_CHOICES

_NUCLEOTIDE_LETTERS = set("ACGTURYSWKMBDHVN-.")

# Temporary annotations set by domainate. They describe the record before rotation, so they are not carried over.
_SEARCH_ANNOTATIONS = {DOMAINATOR_DIST_TO_START_ANNOTATION, DOMAINATOR_DIST_TO_END_ANNOTATION, DOMAINATOR_FOCUS_STRAND_ANNOTATION, DOMAINATOR_HIT_SCORES_ANNOTATION, DOMAINATOR_HIT_NAMES_ANNOTATION, DOMAINATOR_HIT_BEST_SCORE_ANNOTATION}


def _as_seqrecord(rec) -> SeqRecord:
    """Returns rec as a full SeqRecord, converting lean records from the native parser."""
    if LEAN_SEARCH_TYPES and isinstance(rec, LEAN_SEARCH_TYPES):
        return lean_to_seqrecord(materialize_lean_search(rec, None))
    if isinstance(rec, LeanContig):
        return lean_to_seqrecord(rec)
    return rec


def _location_strand(location) -> int:
    """The strand of a location, taken from its first part. Unstranded locations are treated as forward."""
    return -1 if location.parts[0].strand == -1 else 1


def envelope_interval(location, contig_length:int) -> Tuple[int, int, int]:
    """
        The smallest gapless region of a circular contig containing location.

        Returns:
            (start, length, strand): start is the 0-based forward-strand coordinate where the region begins,
                the region may continue across the origin.
    """
    strand = _location_strand(location)
    parts = location.circular_envelope(contig_length, strand).parts
    length = sum(len(part) for part in parts)
    if len(parts) == 1:
        start = int(parts[0].start)
    else: # spans the origin
        start = next((int(part.start) for part in parts if int(part.end) == contig_length), None)
        if start is None:
            raise ValueError(f"Could not find the start of location {location} on a circular contig of length {contig_length}.")
    return start, length, strand


def origin_cut(location, contig_length:int, where:str) -> Tuple[int, int]:
    """
        Where to put the new origin of a circular contig, relative to location.

        Args:
            location: the location of the domain to rotate on
            contig_length: the length of the contig
            where: "start", "end", or "center", relative to the domain in its own orientation

        Returns:
            (cut, strand): cut is a 0-based forward-strand coordinate. The new origin is between bases cut-1 and cut.
                strand is the strand of the domain.
    """
    if where not in WHERE_CHOICES:
        raise ValueError(f"where must be one of {', '.join(WHERE_CHOICES)}, not {where}")
    start, length, strand = envelope_interval(location, contig_length)
    if where == "start":
        offset = 0 if strand == 1 else length
    elif where == "end":
        offset = length if strand == 1 else 0
    else: # center
        offset = length // 2 if strand == 1 else length - length // 2
    return (start + offset) % contig_length, strand


def rotation_location(cut:int, contig_length:int, reverse:bool) -> Union[FeatureLocation, CompoundLocation]:
    """
        A location covering a whole circular contig, starting at a new origin.

        Args:
            cut: the new origin, between bases cut-1 and cut (0-based forward-strand coordinates)
            contig_length: the length of the contig
            reverse: if True, the location is on the reverse strand, starting at base cut-1 and reading toward the old origin.
    """
    if reverse:
        parts = [FeatureLocation(0, cut, -1), FeatureLocation(cut, contig_length, -1)]
    else:
        parts = [FeatureLocation(cut, contig_length, 1), FeatureLocation(0, cut, 1)]
    parts = [part for part in parts if len(part) > 0]
    if len(parts) == 1:
        return parts[0]
    return CompoundLocation(parts, "join")


def rotate_record(rec:SeqRecord, cut:int, reverse:bool, features:Optional[Iterable[SeqFeature]]=None, keep_name:bool=False) -> SeqRecord:
    """
        Rotate a circular contig to a new origin, keeping its annotations and features (features crossing the new origin become joins).

        Args:
            rec: a circular contig
            cut: the new origin, between bases cut-1 and cut (0-based forward-strand coordinates)
            reverse: if True, reverse complement the contig
            features: the features to keep. If None, keep all features.
            keep_name: if False, rename the contig to note how it was rotated.
    """
    location = rotation_location(cut, len(rec), reverse)
    rotated = slice_record_from_location(rec, location, features)
    rotated.annotations = {key: value for key, value in rec.annotations.items() if key not in _SEARCH_ANNOTATIONS}
    rotated.annotations["topology"] = "circular"
    if not keep_name:
        rotated.id = region_id(rec, location)
    return rotated


def _is_protein(rec:SeqRecord) -> bool:
    """True if rec is a protein record, or has letters that aren't IUPAC nucleotide codes (fasta records are read as DNA)."""
    if rec.annotations.get("molecule_type") == "protein":
        return True
    return len(set(str(rec.seq).upper()) - _NUCLEOTIDE_LETTERS) != 0


def _prepare_records(records:Iterable, circularize:bool) -> Iterator[SeqRecord]:
    """Yields the input records that can be rotated, as full circular SeqRecords. Warns about skipped records."""
    for rec in records:
        rec = _as_seqrecord(rec)
        if _is_protein(rec):
            warnings.warn(f"Contig {rec.id} is a protein sequence. Only nucleotide sequences can be rotated. Skipping contig.")
            continue
        if rec.annotations.get("topology") != "circular":
            if not circularize:
                warnings.warn(f"Contig {rec.id} is not circular. Use --circularize to rotate linear contigs and fasta records. Skipping contig.")
                continue
            rec.annotations["topology"] = "circular"
        yield rec


def _search(records:Iterator[SeqRecord], reference_groups, z, evalue, cpu, gene_call, max_mode, batch_size) -> Iterator[Tuple[SeqRecord, Set[int]]]:
    """
        Yields (record, new_hit_ids), where new_hit_ids are the ids of the Domainator features added by searching reference_groups.
    """
    if reference_groups is None:
        for rec in records:
            yield rec, set()
        return
    # domainate buffers all contigs when there are only nucleotide references, so it gets fixed-size chunks.
    while True:
        chunk = list(islice(records, batch_size))
        if len(chunk) == 0:
            break
        # domainate yields the same record objects, keeps their existing features, and appends new hits.
        existing = {id(rec): {id(f) for f in rec.features} for rec in chunk}
        for rec in domainate(chunk, None, z, evalue=evalue, cpu=cpu, pre_parsed_references=reference_groups, gene_call=gene_call, max_mode=max_mode):
            new_hit_ids = {id(f) for f in rec.features if f.type == DOMAIN_FEATURE_NAME and id(f) not in existing[id(rec)]}
            yield rec, new_hit_ids


def _qualifier_float(feature:SeqFeature, name:str, default:float) -> float:
    try:
        return float(feature.qualifiers[name][0])
    except (KeyError, IndexError, ValueError):
        return default


def _candidate_features(rec:SeqRecord, new_hit_ids:Set[int], evalue:float, domains:Optional[Set[str]], databases:Optional[Set[str]], search_hits:bool) -> List[SeqFeature]:
    """The features of rec that the contig can be rotated on."""
    candidates = list()
    for f in rec.features:
        if id(f) in new_hit_ids:
            selected = True
        else: # pre-existing annotation
            if f.type == DOMAIN_FEATURE_NAME:
                selected = domains is not None and f.qualifiers.get("name", [None])[0] in domains
            elif f.type == DOMAIN_SEARCH_BEST_HIT_NAME:
                selected = search_hits
            else:
                selected = False
            if databases is not None and f.qualifiers.get("database", [None])[0] not in databases:
                selected = False
        if selected and _qualifier_float(f, "evalue", float("inf")) <= evalue:
            candidates.append(f)
    return candidates


def _pick_anchor(candidates:List[SeqFeature]) -> SeqFeature:
    """The candidate with the highest score, then lowest evalue, then lowest start."""
    return min(candidates, key=lambda f: (-_qualifier_float(f, "score", float("-inf")), _qualifier_float(f, "evalue", float("inf")), int(f.location.start)))


def rotate_contigs(records:Iterable, references:Optional[List[str]]=None, pre_parsed_references=None, z:int=1000, evalue:float=0.001,
                   domains:Optional[Iterable[str]]=None, databases:Optional[Iterable[str]]=None, search_hits:bool=False, where:str="start",
                   keep_direction:bool=False, circularize:bool=False, keep_name:bool=False, keep_annotation:bool=False,
                   gene_call:Optional[str]=None, cpu:int=0, max_mode:bool=False, batch_size:int=1000) -> Iterator[SeqRecord]:
    """
        Rotate circular contigs so that their origin is set relative to a domain annotation.

        Returns a generator of rotated SeqRecords. Contigs that can't be rotated are skipped with a warning.
        Arguments are validated when this function is called, before any records are read.

        Args:
            records: an iterable of nucleotide SeqRecords
            references: hmm, cm, or fasta files to search against the records, to find new domains to rotate on.
            pre_parsed_references: the output of domainate.read_references, used instead of references.
            z: passed as the -Z parameter to hmmsearch.
            evalue: the maximum evalue of new and pre-existing domains to rotate on.
            domains: names of pre-existing Domainator annotations to rotate on.
            databases: only consider pre-existing annotations from these databases. If None, consider all databases.
            search_hits: if True, rotate on pre-existing domain_search annotations.
            where: "start", "end", or "center". Where to put the origin, relative to the domain in its own orientation.
            keep_direction: if True, don't reverse complement contigs with the domain on the reverse strand.
            circularize: if True, rotate linear contigs and fasta records as though they were circular. If False, skip them.
            keep_name: if True, keep the original contig names. If False, rename contigs to note how they were rotated.
            keep_annotation: if True, keep the new annotation of the domain that the contig was rotated on.
                Other new annotations are never kept.
            gene_call: None, "all", or "unannotated". Passed to domainate to call genes with prodigal before the search.
            cpu: number of cpus to use for the search. 0 means use all available cpus.
            max_mode: run the search in maximum sensitivity mode.
            batch_size: number of contigs to search at a time.
    """
    if where not in WHERE_CHOICES:
        raise ValueError(f"where must be one of {', '.join(WHERE_CHOICES)}, not {where}")
    if gene_call not in (None, "all", "unannotated"):
        raise ValueError("gene_call must be one of None, 'all', or 'unannotated'")
    if batch_size < 1:
        raise ValueError("batch_size must be at least 1.")
    domains = set(domains) if domains else None
    databases = set(databases) if databases is not None else None

    reference_groups = pre_parsed_references
    if references is not None:
        reference_groups = read_references(references, None)
    if reference_groups is None and domains is None and not search_hits:
        raise ValueError("Specify references, domains, or search_hits to choose domains to rotate on.")
    if reference_groups is None and gene_call is not None:
        raise ValueError("gene_call only applies when searching references.")
    if gene_call == "all" and (domains is not None or search_hits):
        raise ValueError("gene_call 'all' deletes pre-existing annotations, so it can't be combined with domains or search_hits.")

    return _rotate_contigs(records, reference_groups, z, evalue, domains, databases, search_hits, where, keep_direction, circularize, keep_name, keep_annotation, gene_call, cpu, max_mode, batch_size)


def _rotate_contigs(records, reference_groups, z, evalue, domains, databases, search_hits, where, keep_direction, circularize, keep_name, keep_annotation, gene_call, cpu, max_mode, batch_size):
    """The generator behind rotate_contigs. It is separate so that rotate_contigs can validate its arguments before any records are read."""
    for rec, new_hit_ids in _search(_prepare_records(records, circularize), reference_groups, z, evalue, cpu, gene_call, max_mode, batch_size):
        candidates = _candidate_features(rec, new_hit_ids, evalue, domains, databases, search_hits)
        if len(candidates) == 0:
            warnings.warn(f"No qualifying domain found on contig {rec.id}. Skipping contig.")
            continue
        anchor = _pick_anchor(candidates)
        cut, strand = origin_cut(anchor.location, len(rec), where)
        reverse = strand == -1 and not keep_direction
        features = [f for f in rec.features if id(f) not in new_hit_ids or (keep_annotation and f is anchor)]
        yield rotate_record(rec, cut, reverse, features, keep_name)


def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, 'Names of input genbank or fasta files. Must be nucleotide sequences.')
    add_seq_outputs(parser)

    parser.add_argument("-r", "--reference", required=False, type=str, default=None, nargs="+",
                        help="reference files to search against the input, to find new domains to rotate on. Supported types are protein HMMs, protein FASTA, nucleotide FASTA, nucleotide HMMs, and infernal CM files.")
    parser.add_argument("-Z", default=1000, type=int,
                        help="Passed as the -Z parameter to hmmsearch: Assert that the total number of targets in your searches is <x>, for the purposes of per-sequence E-value calculations, rather than the actual number of targets seen. "
                        "Set to 0 to use the number of actual target sequences in the input files (not available when reading from stdin).")
    parser.add_argument("-e", "--evalue", type=PositiveFloat, default=0.001,
                        help="threshold E value for domains to rotate on. Both new and pre-existing annotations must pass this threshold. Must be <=10.")

    parser.add_argument("--domains", type=str, default=None, nargs="+",
                        help="Rotate on pre-existing domainator annotations with these names.")
    parser.add_argument("--domains_file", type=str, default=None,
                        help="text file containing names of pre-existing domainator annotations to rotate on.")
    parser.add_argument("--databases", type=str, default=None, nargs="+",
                        help="Consider only pre-existing annotations from these databases. default: all databases.")
    add_annotation_type(parser, help="Which pre-existing annotations to rotate on.")

    parser.add_argument("--anchor", type=str, default="start", choices=WHERE_CHOICES,
                        help="where to put the origin, relative to the chosen domain in its own orientation (not relative to the original orientation of the contig). "
                        "start: the first base of the domain becomes the first base of the contig. end: the last base of the domain becomes the last base of the contig. center: the origin is placed in the middle of the domain.")
    parser.add_argument("--keep_direction", action="store_true", default=False,
                        help="by default, contigs are reverse complemented if necessary to put the chosen domain on the forward strand. Setting this option will keep the original orientation of the contig.")
    parser.add_argument("--circularize", action="store_true", default=False,
                        help="treat fasta records and linear genbank records as circular, before searching and rotating. By default, they are skipped.")
    parser.add_argument("--keep_name", action="store_true", default=False,
                        help="by default, rotated contigs will be renamed to note the original coordinates of their first and last bases, and whether they were reverse complemented. Set this option to keep contigs with their original names.")
    parser.add_argument("--keep_annotation", action="store_true", default=False,
                        help="by default, new annotations from the reference search are not written to the output. Setting this option will keep the new annotation of the domain that the contig was rotated on. Has no effect if the contig was rotated on a pre-existing annotation.")

    parser.add_argument("--gene_call", type=str, default=None, choices={"all", "unannotated"}, required=False,
                        help="When activated, new CDS annotations will be added with Prodigal in Metagenomic mode before the search. If 'all', then any existing CDS and domainator annotations will be deleted and all contigs will be re-annotated. If 'unannotated', then only contigs without CDS annotations will be annotated. "
                        "Needed to search protein references against contigs without CDS annotations, such as fasta records.")
    parser.add_argument("--cpu", type=int, default=0,
                        help="the number of cores of the cpu which are used at a time to run the search (0: one thread per physical core)")
    parser.add_argument("--max_mode", action="store_true", default=False,
                        help="Run hmmsearch, phmmer, or nhmmer in maximum sensitivity mode, which is much slower, but more sensitive.")
    parser.add_argument("--batch_size", type=int, default=1000,
                        help="How many contigs to search at one time in a batch. Increasing this number might improve speeds at the cost of memory.")

    parser.add_argument("--config", action=ActionConfigFile)
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)

    if params.evalue > 10:
        parser.error("-e/--evalue must be <= 10.")

    domains = list_and_file_to_dict_keys(params.domains, params.domains_file, as_set=True)

    input_files = params.input

    z = params.Z
    if z == 0:
        if any(path in STDIN_PATHS for path in params.input):
            parser.error("-Z 0 counts the sequences in the input files, so it can't be used when reading from stdin.")
        z = sum(utils.count_peptides(file_path) for file_path in params.input)

    rotated = rotate_contigs(
        parse_seqfiles(input_files, None, default_molecule_type="DNA"),
        references=params.reference,
        z=z,
        evalue=params.evalue,
        domains=domains,
        databases=params.databases,
        search_hits=params.annotation_type == "domain_search",
        where=params.anchor,
        keep_direction=params.keep_direction,
        circularize=params.circularize,
        keep_name=params.keep_name,
        keep_annotation=params.keep_annotation,
        gene_call=params.gene_call,
        cpu=params.cpu,
        max_mode=params.max_mode,
        batch_size=params.batch_size,
    )

    check_outputs(parser, params, ["genbank_out", "fasta_out"])

    write_seq_outputs(rotated, genbank=params.genbank_out, fasta=params.fasta_out)
