"""Converts annotated CDSs into annotated protein sequences

Takes a domainator-annotated genbank file 
Creates peptide genbank or fasta files with the domainator annotations from the CDSs found in the original
genbank file

By default, extracts peptides from all CDSs in the genbank file.

If any criteria are supplied, then only peptides from CDSs that meet the criteria will be extracted.

Criteria are combined with the following logic:
    (all OR target_CDSs OR target_domains OR unannotated OR search_hits) AND (target_contigs OR strand)

"""

import sys
import argparse
from typing import List, Optional, Dict, Set
from domainator.Bio.Seq import Seq
from domainator.Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation
from domainator.Bio.SeqRecord import SeqRecord
from domainator.Bio.Seq import Seq
from domainator.utils import list_and_file_to_dict_keys, parse_seqfiles, write_genbank, DomainatorCDS, get_sources, copy_feature, slice_record_from_location
import warnings
from domainator.Bio import SeqIO
from domainator import __version__, RawAndDefaultsFormatter
from itertools import chain

# TODO: account for different databases in domain selection
# TODO: how to handle stop codons? Right now they are translated if they are covered by in the CDS annotation.

# def get_seq_record(start, end, id, sequence):
#     peptide = Seq.translate(Seq.reverse_complement(sequence[start:end]))
#     record = SeqRecord(seq=peptide, id=id, name=id)
#     record.annotations["molecule_type"] = "protein"
#     return record

def dna_to_peptide_location(location):

    if isinstance(location, CompoundLocation):
        return CompoundLocation([dna_to_peptide_location(part)
                                    for part in location.parts], operator=location.operator)
    start = int(location.start) // 3
    end = int(location.end) // 3
    return FeatureLocation(start, end, strand=location.strand)

def extract_peptides(records, evalue, target_domains:Optional[Set], target_cdss, keep_cds_feature=True, cds_id_type="name", search_hits=False, strand=None, unannotated=False, _from_domain_search=False, extract_all=False, invert=False, _domain_search_negatives:Optional[Set[str]]=None, databases=None, keep_name=False):
    """
        Extracts all peptide sequences from each SeqRecord in records.
        Yields the extracted peptides as SeqRecord objects.
        Args:
            records: an iterable of SeqRecords
            evalue: the evalue cutoff for domain annotations to be put onto the new genbank file
            target_domains: only extract peptides containing these domains
            target_cdss: only consider peptides from cdss with these names
            keep_cds_feature: if True, keep the CDS feature on the peptide record. If False, remove it.
            cds_id_type: if "name", then target_cdss is a list of CDS names. If "num", then target_cdss is a list of CDS numbers.
            search_hits: if True, extract peptides from the search hits (from domain_search.py).
            strand: only extract peptides of CDSs on the specified strand. If None, then extract peptides from CDSs on both strands. forward strand is "f", reverse strand is "r".
            unannotated: if True, extract peptides from CDSs with no domain annotations.
            _from_domain_search: if True, then the input records are from a domain search, and we can skip all the annotation propagation.
            extract_all: if True, then extract peptides from all CDSs, regardless of the other criteria.
            invert: if True, then invert the CDS selection criteria. i.e. return peptides from CDSs that don't match the CDS selection criteria. (Only applies to CDS selection, not contig selection)
            _domain_search_negatives: if not None, then this is a set of domain names that are considered negative hits. If a CDS has the best search hit as a negative, then it will be skipped.
            databases: if not None, then only consider domains from these databases when filtering by domain. default: all databases.
            keep_name: if True, then keep the original name of the contig. If False, then the name will be the name of the CDS.
    """
    
    # TODO: if fasta_out is set, we can skip all the annotation propagation, because it won't be used. This should speed things up for that case.
    
    databases = set(databases) if databases is not None else None
    for contig in records:
        cdss = DomainatorCDS.list_from_contig(contig, evalue, skip_pseudo=True)
        source_features = get_sources(contig)
        if _from_domain_search:
            hit_scores = contig.get_hit_scores()
            hit_names = contig.get_hit_names()
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
                if (target_cdss is not None and target_cds_id in target_cdss):
                    keep = True
                if target_domains is not None and len(target_domains.intersection(cds.get_domain_names(databases))) > 0:
                    keep = True
                if search_hits and cds.domain_search_feature is not None:
                    keep = True
                if unannotated and len(cds.get_domain_names(databases)) == 0:
                    keep = True
            
            if strand is not None:
                if strand == "f" and cds.feature.location.parts[0].strand == -1: # if strand is forward, but the cds is on the reverse strand, then don't keep it
                    keep = False
                elif strand == "r" and cds.feature.location.parts[0].strand != -1:
                    keep = False
            
            if invert:
                keep = not keep

            if keep:
                if keep_cds_feature:
                    cds_feature_tmp = copy_feature(cds.feature)
                    cds_feature_tmp.qualifiers["source_contig"] = [contig.id]

                    if "translation" in cds_feature_tmp.qualifiers: # remove the translation qualifier, because it's not necessary in the output
                        del cds_feature_tmp.qualifiers["translation"]
                    
                    features = chain(source_features, (x for x in (cds.domain_search_feature,) if x is not None), cds.domain_features, (cds_feature_tmp,))
                else:
                    features = chain(source_features, (x for x in (cds.domain_search_feature,) if x is not None), cds.domain_features)
                
                contig_slice = slice_record_from_location(contig, cds.feature.location, features, truncate_features=True)
                if ('translation' not in cds.feature.qualifiers) or (len(cds.feature.qualifiers["translation"][0]) == 0):
                    peptide = Seq.translate(contig_slice.seq, cds=False) # cds=False because we don't want to throw errors on weird stuff.
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
                    feature.location = dna_to_peptide_location(feature.location)
                    rec.features.append(feature)

                if _from_domain_search:
                    rec.set_hit_best_score(best_score)

                yield rec

def main(argv):
    parser = argparse.ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', nargs='+', required=False, default=None,
                       help="Genbank filenames. If not supplied, reads from stdin.")

    parser.add_argument('-o', '--output', default=None, required=False,
                        help="the name of the output genbank file. If not supplied writes to stdout.")

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

    parser.add_argument('--search_hits', action='store_true', default=False, help="extract peptides from the search hits (from domain_search.py).")

    parser.add_argument('--strand', type=str, default=None, choices=["f", "r"], help="Only extract peptides from CDSs on the specified strand.")

    parser.add_argument('--invert', action='store_true',
                        help="Invert the CDS selection criteria. i.e. return peptides from CDSs that don't match the CDS selection criteria. (Only applies to CDS selection, not contig selection)")
    
    parser.add_argument('--cdss', type=str, default=None, nargs='+',
                        help="only look for peptides on cdss with these names") 
    parser.add_argument('--cdss_file', type=str, default=None,
                        help="only look for peptides on cdss with names listed in this text file.")

    parser.add_argument('--keep_name', action='store_true', default=False,
                       help="by default, extracted regions will be renamed to note their original coordinates and whether "
                       "they were reverse complemented. Set this option to keep contigs with their original names. "
                       "WARNING: if you aren't careful, this can easily result in duplicate names.")

    parser.add_argument('-e', '--evalue', default=100000000,
                        help="the evalue cutoff for domain annotations to be put onto the new genbank file")
    
    parser.add_argument('--fasta_out', action='store_true', default=False,
                        help="makes output a fasta file when activated")

    params = parser.parse_args(argv)
   
    ### Figure out what input and output files ####
    
    filetype_override = None
    if params.input is None:
        genbanks = [sys.stdin]
        filetype_override = "genbank"
    else:
        genbanks = params.input

    if params.output is None:
        out = sys.stdout
    else:
        out = open(params.output, "w")




    target_contigs = list_and_file_to_dict_keys(params.contigs, params.contigs_file, as_set=True)
    target_domains = list_and_file_to_dict_keys(params.domains, params.domains_file, as_set=True)
    target_cdss = list_and_file_to_dict_keys(params.cdss, params.cdss_file, as_set=True)

    extract_all = False
    if target_domains is None and target_cdss is None and not params.unannotated and not params.search_hits:
        extract_all = True
        if params.invert:
            warnings.warn("Invert supplied with no selection criteria. Invert will be ignored.")
            params.invert = False


    # Run
    extracted_peptides_iterator = extract_peptides(
                    parse_seqfiles(genbanks, target_contigs, filetype_override=filetype_override),
                    params.evalue,
                    target_domains,
                    target_cdss,
                    keep_cds_feature=True,
                    search_hits = params.search_hits,
                    strand = params.strand,
                    unannotated = params.unannotated,
                    extract_all=extract_all,
                    invert=params.invert,
                    databases=params.databases,
                    keep_name = params.keep_name,
                    )

    if params.fasta_out:
        SeqIO.write(extracted_peptides_iterator, out, "fasta")
    else:
        write_genbank(extracted_peptides_iterator, out, default_molecule_type="protein")


    if params.output is not None:
        out.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
