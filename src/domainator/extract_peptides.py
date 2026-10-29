"""Converts annotated CDSs into annotated protein sequences

Takes a domainator-annotated genbank file 
Creates peptide genbank or fasta files with the domainator annotations from the CDSs found in the original
genbank file

By default, extracts peptides from all CDSs in the genbank file.

If any criteria are supplied, then only peptides from CDSs that meet the criteria will be extracted.

Criteria are combined with the following logic:
    (all OR target_CDSs OR target_domains OR unannotated OR domain_search hits) AND (target_contigs OR strand) AND partial

"""

from jsonargparse import ActionConfigFile
from jsonargparse.typing import PositiveFloat
from typing import Optional, Set, Iterable, Iterator
from domainator.cli.common import add_annotation_type, add_has_search_hit, add_input, add_seq_outputs, check_outputs, declare_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, CompoundLocation, BeforePosition, AfterPosition
from Bio.SeqRecord import SeqRecord
from domainator.cli.common import list_and_file_to_dict_keys, partial_help
from domainator.seqio import parse_seqfiles
from domainator.model import DomainatorCDS, get_sources, copy_feature, check_annotation_type, get_hit_names, get_hit_scores, set_hit_best_score
from domainator.locations import slice_record_from_location
from domainator.partial import codon_offset, fragment_status_allowed, get_fragment_status, PARTIAL_CHOICES
import warnings
from itertools import chain
from domainator.cli.common import add_taxonomy_filters, taxonomy_from_args
from domainator.taxonomy import NCBITaxonomy
from domainator.taxonomy.records import compile_taxonomy_filter

# TODO: account for different databases in domain selection
# TODO: how to handle stop codons? Right now they are translated if they are covered by in the CDS annotation.

# def get_seq_record(start, end, id, sequence):
#     peptide = Seq.translate(Seq.reverse_complement(sequence[start:end]))
#     record = SeqRecord(seq=peptide, id=id, name=id)
#     record.annotations["molecule_type"] = "protein"
#     return record

def _dna_to_peptide_position(position, offset):
    """Converts a nucleotide position (relative to the start of the CDS) to a peptide position, keeping '<' and '>' markers."""
    peptide_position = max(0, (int(position) - offset) // 3)
    if isinstance(position, BeforePosition):
        return BeforePosition(peptide_position)
    if isinstance(position, AfterPosition):
        return AfterPosition(peptide_position)
    return peptide_position

def dna_to_peptide_location(location, offset=0):
    """
        Converts a nucleotide location, relative to the start of a CDS, to peptide coordinates.
        offset is the reading frame offset of the CDS (codon_start - 1).
    """
    if isinstance(location, CompoundLocation):
        return CompoundLocation([dna_to_peptide_location(part, offset)
                                    for part in location.parts], operator=location.operator)
    start = _dna_to_peptide_position(location.start, offset)
    end = _dna_to_peptide_position(location.end, offset)
    return FeatureLocation(start, end, strand=location.strand)

def extract_peptides(records: Iterable[SeqRecord], *, evalue: Optional[float] = None, domains: Optional[Set] = None,
                     cds_ids: Optional[Set[str]] = None,
                     keep_cds_feature: bool = True, cds_id_type: str = "name", annotation_type: str = "domainator",
                     strand: Optional[str] = None, unannotated: bool = False, _from_domain_search: bool = False,
                     extract_all: bool = False, invert: bool = False, _domain_search_negatives: Optional[Set[str]] = None,
                     databases: Optional[Set[str]] = None, keep_name: bool = False, name_field: Optional[str] = None,
                     partial: str = "include", has_search_hit: bool = False,
                     ncbi_taxonomy: Optional[NCBITaxonomy] = None,
                     include_taxids: Optional[Iterable[int]] = None, exclude_taxids: Optional[Iterable[int]] = None,
                     taxonomy_expr: Optional[str] = None, taxonomy_unknown: str = "exclude") -> Iterator[SeqRecord]:
    """
        Extracts all peptide sequences from each SeqRecord in records.
        Yields the extracted peptides as SeqRecord objects.
        Args:
            records: an iterable of SeqRecords
            evalue: the evalue cutoff for domain annotations to be put onto the new genbank file (None: no cutoff)
            domains: only extract peptides containing these domains
            cds_ids: only consider peptides from cdss with these names
            keep_cds_feature: if True, keep the CDS feature on the peptide record. If False, remove it.
            cds_id_type: if "name", then cds_ids is a list of CDS names. If "num", then cds_ids is a list of CDS numbers.
            annotation_type: which annotations supply domain names: "domainator", "domain_search", or "all".
            has_search_hit: also select CDSs with search hits passing evalue, ignoring domain
                names and databases. Participates in invert.
            strand: only extract peptides of CDSs on the specified strand. If None, then extract peptides from CDSs on both strands. forward strand is "f", reverse strand is "r".
            unannotated: if True, extract peptides from CDSs with no domain annotations.
            _from_domain_search: if True, then the input records are from a domain search, and we can skip all the annotation propagation.
            extract_all: if True, then extract peptides from all CDSs, regardless of the other criteria.
            invert: if True, then invert the CDS selection criteria. i.e. return peptides from CDSs that don't match the CDS selection criteria. (Only applies to CDS selection, not contig selection)
            _domain_search_negatives: if not None, then this is a set of domain names that are considered negative hits. If a CDS has the best search hit as a negative, then it will be skipped.
            databases: if not None, then only consider domains from these databases when filtering by domain. default: all databases.
            keep_name: if True, then keep the original name of the contig. If False, then the name will be the name of the CDS.
            partial: select CDSs by fragment status: "include" (default, all), "exclude" (skip fragments), or "only" (only fragments).
                A CDS is a fragment if its location has '<' or '>' (see utils.get_fragment_status). Not affected by invert.
            ncbi_taxonomy: database required when taxonomy criteria are supplied.
            include_taxids: include these taxids and descendants (OR).
            exclude_taxids: exclude these taxids and descendants, overriding inclusion.
            taxonomy_expr: boolean lineage expression, exclusive with include/exclude_taxids.
            taxonomy_unknown: unresolved taxonomy policy: exclude (default), include,
                evaluate against an empty lineage, or error. Taxonomy is resolved on the
                original CDS, before translation, and is not affected by invert.
    """
    check_annotation_type(annotation_type)
    taxonomy_filter = compile_taxonomy_filter(ncbi_taxonomy, include_taxids, exclude_taxids,
                                              taxonomy_expr, taxonomy_unknown)
    
    # TODO: if fasta_out is set, we can skip all the annotation propagation, because it won't be used. This should speed things up for that case.
    
    databases = set(databases) if databases is not None else None
    name_precedence = ["protein_id", "gene_id", "locus_tag", "gene"]
    if name_field is not None:
        name_precedence = [name_field] + ["protein_id", "gene_id", "locus_tag", "gene"]

    for contig in records:
        cdss = DomainatorCDS.list_from_contig(contig, evalue, skip_pseudo=True, name_precedence=name_precedence)
        source_features = get_sources(contig)
        if _from_domain_search:
            hit_scores = get_hit_scores(contig)
            hit_names = get_hit_names(contig)
        for cds in cdss:
            keep = False
            target_cds_id = cds.name if cds_id_type == "name" else cds.num
            if _from_domain_search:
                if cds.index in hit_scores:
                    if _domain_search_negatives is None or hit_names[cds.index] not in _domain_search_negatives:
                        # if _domain_search_negatives is not None, then we are doing a negative selection, so skip any cds that has a negative domain
                        keep = True
                        best_score = hit_scores[cds.index]
            else:
                if extract_all:
                    keep = True
                if (cds_ids is not None and target_cds_id in cds_ids):
                    keep = True
                domain_names = cds.get_domain_names(databases, annotation_type=annotation_type)
                if domains is not None and len(domains.intersection(domain_names)) > 0:
                    keep = True
                if has_search_hit and cds.domain_search_feature is not None:
                    keep = True
                if unannotated and len(domain_names) == 0:
                    keep = True
            
            if strand is not None:
                if strand == "f" and cds.feature.location.parts[0].strand == -1: # if strand is forward, but the cds is on the reverse strand, then don't keep it
                    keep = False
                elif strand == "r" and cds.feature.location.parts[0].strand != -1:
                    keep = False
            
            if invert:
                keep = not keep

            if keep and partial != "include" and not fragment_status_allowed(get_fragment_status(cds.feature), partial):
                keep = False
            if keep and taxonomy_filter is not None:
                keep = taxonomy_filter.matches(contig, cds.feature)

            if keep:
                if keep_cds_feature:
                    cds_feature_tmp = copy_feature(cds.feature)
                    cds_feature_tmp.qualifiers["source_contig"] = [contig.id]

                    if "translation" in cds_feature_tmp.qualifiers: # remove the translation qualifier, because it's not necessary in the output
                        del cds_feature_tmp.qualifiers["translation"]
                    cds_feature_tmp.qualifiers.pop("codon_start", None) # the reading frame is meaningless in peptide coordinates
                    
                    features = chain(source_features, (x for x in (cds.domain_search_feature,) if x is not None), cds.domain_features, (cds_feature_tmp,))
                else:
                    features = chain(source_features, (x for x in (cds.domain_search_feature,) if x is not None), cds.domain_features)
                
                contig_slice = slice_record_from_location(contig, cds.feature.location, features)
                if ('translation' not in cds.feature.qualifiers) or (len(cds.feature.qualifiers["translation"][0]) == 0):
                    peptide = Seq.translate(contig_slice.seq[codon_offset(cds.feature):], cds=False) # cds=False because we don't want to throw errors on weird stuff.
                else:
                    peptide = Seq(cds.feature.qualifiers["translation"][0])

                new_name = cds.name
                if keep_name:
                    new_name = contig.id
                rec = SeqRecord(seq=peptide, id=new_name, name=new_name)

                rec.annotations["molecule_type"] = "protein"
                if "source" in contig.annotations:
                    rec.annotations["source"] = contig.annotations["source"]
                if "organism" in contig.annotations:
                    rec.annotations["organism"] = contig.annotations["organism"]
                if "taxonomy" in contig.annotations:
                    rec.annotations["taxonomy"] = contig.annotations["taxonomy"]
                rec.dbxrefs=contig.dbxrefs.copy()
                rec.description = contig.description

                rec.features = list()
                for feature in contig_slice.features:
                    feature.location = dna_to_peptide_location(feature.location, codon_offset(cds.feature))
                    rec.features.append(feature)

                if _from_domain_search:
                    set_hit_best_score(rec, best_score)

                yield rec

def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, 'Genbank filenames.')

    add_seq_outputs(parser)
    add_taxonomy_filters(parser)

    parser.add_argument('--contigs', default=[], nargs='+',
                        help="only look for peptides on contigs with ids in this list. Additive with --contigs_file. default: all contigs.")
    
    parser.add_argument('--contigs_file', default=None, 
                        help="only look for peptides on contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs.")

    parser.add_argument('--domains', type=str, default=None, nargs='+',
                        help="Only extract peptides containing these domains")  # TODO: add some way to specify the database of a domain
    parser.add_argument('--domains_file', type=str, default=None,
                        help="text file containing names of domains, where peptides containing these domains will be extracted.")
    
    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Consider only domains from these databases when filtering by domain. default: all databases.")
    
    parser.add_argument('--unannotated', action='store_true', default=False, help="extract peptides from CDSs with no domain annotations.")

    add_annotation_type(parser, help="Which annotations to select peptides by.")
    add_has_search_hit(parser)

    parser.add_argument('--strand', type=str, default=None, choices=["f", "r"], help="Only extract peptides from CDSs on the specified strand.")
    parser.add_argument('--partial', type=str, default="include", choices=PARTIAL_CHOICES, help=partial_help("Which CDSs to extract peptides from") + " Not affected by --invert.")

    parser.add_argument('--invert', action='store_true',
                        help="Invert the CDS selection criteria. i.e. return peptides from CDSs that don't match the CDS selection criteria. (Only applies to CDS selection, not contig selection)")
    
    parser.add_argument('--cds', type=str, default=None, nargs='+',
                        help="only look for peptides on cdss with these names") 
    parser.add_argument('--cds_file', type=str, default=None,
                        help="only look for peptides on cdss with names listed in this text file.")

    parser.add_argument('--keep_name', action='store_true', default=False,
                       help="by default, extracted regions will be renamed to note their original coordinates and whether "
                       "they were reverse complemented. Set this option to keep contigs with their original names. "
                       "WARNING: if you aren't careful, this can easily result in duplicate names.")

    parser.add_argument('--name_field', type=str, default=None,
                        help='The CDS qualifier to preferentially use as the name of the peptide.'
                        ' If not supplied, then qualifiers will be checked in order of: ["protein_id", "gene_id", "locus_tag", "gene"]')

    parser.add_argument('-e', '--evalue', type=PositiveFloat, default=None,
                        help="the evalue cutoff (max) for domain annotations to be put onto the new genbank file. Default: no cutoff.")
    

    parser.add_argument('--config', action=ActionConfigFile)
    declare_outputs(parser, ["genbank_out", "fasta_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
   
    ### Figure out what input and output files ####
    
    genbanks = params.input
    

    check_outputs(parser, params)

    target_contigs = list_and_file_to_dict_keys(params.contigs, params.contigs_file, as_set=True)
    target_domains = list_and_file_to_dict_keys(params.domains, params.domains_file, as_set=True)
    target_cdss = list_and_file_to_dict_keys(params.cds, params.cds_file, as_set=True)

    extract_all = False
    if target_domains is None and target_cdss is None and not params.unannotated and not params.has_search_hit:
        extract_all = True
        if params.invert:
            warnings.warn("Invert supplied with no selection criteria. Invert will be ignored.")
            params.invert = False

    # Run
    extracted_peptides_iterator = extract_peptides(
                    parse_seqfiles(genbanks, target_contigs),
                    evalue=params.evalue,
                    domains=target_domains,
                    cds_ids=target_cdss,
                    keep_cds_feature=True,
                    annotation_type = params.annotation_type,
                    has_search_hit=params.has_search_hit,
                    strand = params.strand,
                    unannotated = params.unannotated,
                    extract_all=extract_all,
                    invert=params.invert,
                    databases=params.databases,
                    keep_name = params.keep_name,
                    name_field = params.name_field,
                    partial = params.partial,
                    ncbi_taxonomy=taxonomy_from_args(params),
                    include_taxids=params.include_taxids,
                    exclude_taxids=params.exclude_taxids,
                    taxonomy_expr=params.taxonomy_expr,
                    taxonomy_unknown=params.taxonomy_unknown,
                    )

    write_seq_outputs(extracted_peptides_iterator, genbank=params.genbank_out, fasta=params.fasta_out)
