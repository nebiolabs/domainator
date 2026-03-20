""" Annotate sequence files with protein or nucleic-acid references.

Protein references are searched against CDSs or protein records with hmmsearch or phmmer.
Nucleic-acid references are searched against whole nucleotide contigs with nhmmer or infernal cmsearch.
Hits passing the threshold are added to the output GenBank file as new Domainator features.

Supported references include protein HMMs, protein FASTA, nucleotide FASTA, nucleotide HMMs, and infernal CM files.
Protein inputs cannot be searched with nucleotide references.

NOTES: genbank files that contain translation annotations in their CDS features will run much faster than those without translation annotations.
    domainate.py stores the entire reference set in memory, so it is not suitable for cases where the reference
    set is larger than your system memory (this may change in future versions).
"""
import warnings
warnings.filterwarnings("ignore", module='numpy')
import sys
from jsonargparse import ArgumentParser, ActionConfigFile
from typing import NamedTuple, List, Dict, Set, Tuple, Optional, Union, Iterable, Iterator
from domainator.Bio.Seq import Seq
from domainator.Bio.SeqRecord import SeqRecord
from domainator.Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from domainator import utils, DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.utils import get_cds_unique_name, parse_seqfiles, write_genbank, read_hmms, read_infernal_cms, get_file_type, read_pyhmmer_fastas, read_pyhmmer_peptide_fastas, filter_by_taxonomy, pyhmmer_decode
import pyhmmer
import pyinfernal
from domainator import __version__, RawAndDefaultsFormatter
import pyrodigal
from pathlib import Path
from domainator.Taxonomy import NCBITaxonomy
from domainator import foldseek as foldseek_lib

#TODO: add support for distinguishing complete domains from cutoff domains. For example see the pfam or hmmer domain graphics. I think this can be done by comparing the size of the domain consensus to the domain border coords (from domain.Alignment) on the hits.
#TODO: also look at: https://github.com/Larofeticus/hpc_hmmsearch, https://www.higithub.com/apcamargo/repo/hpc_pfam_search, https://docs.nersc.gov/performance/case-studies/hmmer3/
#TODO: add progress bar
#TODO: For phmmer hits, show the alignment as part of the annotation
#TODO: add database name as DEFINITION or ACCESSION
#TODO: I think there may be some issue in pyhmmer where it makes a pool of all available cores, and then just uses one of them.



MAX_PROTEIN_SIZE = 100_000
NUCLEIC_ACID_REFERENCE_GROUPS = {"nhmmer", "infernal"}


def _is_nucleic_acid_alphabet(alphabet) -> bool:
    return alphabet.is_dna() or alphabet.is_rna()


def _fasta_is_nucleic_acid(file_name: str) -> bool:
    with pyhmmer.easel.SequenceFile(file_name, digital=True) as seq_file:
        for seq in seq_file:
            return _is_nucleic_acid_alphabet(seq.alphabet)
    return False


def _hmm_is_nucleic_acid(file_name: str) -> bool:
    with pyhmmer.plan7.HMMFile(file_name) as hmm_file:
        for hmm in hmm_file:
            return _is_nucleic_acid_alphabet(hmm.alphabet)
    return False


def has_nucleic_acid_references(reference_groups) -> bool:
    return any(group in reference_groups for group in NUCLEIC_ACID_REFERENCE_GROUPS)


def get_reference_length(reference) -> int:
    try:
        return reference.M
    except AttributeError:
        return len(reference.sequence)


def get_max_reference_length(reference_groups, algorithm: str) -> int:
    if algorithm not in reference_groups:
        return 0
    max_length = 0
    for db_dict in reference_groups[algorithm].values():
        for reference in db_dict.values():
            max_length = max(max_length, get_reference_length(reference))
    return max_length


def get_contig_search_sequence(record: SeqRecord, extension_length: int = 0) -> str:
    sequence = ''.join(filter(str.isalpha, str(record.seq)))
    if record.annotations.get("topology") == "circular" and extension_length > 0:
        sequence += sequence[:extension_length]
    return sequence


def validate_input_molecule_types(molecule_types, reference_groups):
    if "protein" in molecule_types and has_nucleic_acid_references(reference_groups):
        raise ValueError("Cannot search protein inputs with a nucleotide reference.")


def get_input_molecule_types(input_files, default_molecule_type="protein"):
    molecule_types = set()
    for input_file in input_files:
        for rec in parse_seqfiles((input_file,), None, default_molecule_type=default_molecule_type, max_recs=1):
            molecule_types.add(rec.annotations["molecule_type"])
    return molecule_types

#{"name": <hmm_name> as str, "desc": <hmm_description> as str, "evalue": <domain_evalue> as float,
#        "score": <domain_score> as float, "start": <current_left> as int, "end": <current_right> as int}
class SearchResult(NamedTuple):
    name: str
    """hmm profile name (hmm_name)"""
    desc: str
    """hmm profile description (description)"""
    acc: str
    """hmm profile accession (hmm_accession)"""
    evalue: float
    """hit evalue (i_evalue)"""
    score: float
    """hit bitscore (score)"""
    start: int
    """start position on sequence (env_from)"""
    end: int
    """end position on sequence (env_to)"""
    database: str
    """name of hmm file containing profile"""
    identity: float
    """percent identity of hit (based on alignment length)"""
    rstart: int
    """start position on reference hmm (hmm_from)"""
    rend: int
    """end position on reference hmm (hmm_to)"""
    rlen: Union[str, int] = "" #TODO: should default be "" or None?
    """length of reference hmm (hmm_length)"""
    program:str = "."
    """name of program that made the hit (hmmsearch, phmmer, foldseek)"""
    strand: int = 1
    """strand of the hit on the target sequence"""
    target_kind: str = "protein"
    """whether the hit targets a protein/CDS or a nucleic acid contig span"""




def hmmer_hits_to_search_results(hits, references, evalue, db_name, min_evalue, program, out):

    for top_hits in hits: # each top_hits covers a single hmm
        for hit in top_hits:
            hit_name = pyhmmer_decode(hit.name).split(",")
            contig_index = int(hit_name[0])
            cds_index = int(hit_name[1])
            for domain in hit.domains:
                if domain.i_evalue < evalue and domain.i_evalue >= min_evalue:
                    domain_name = pyhmmer_decode(domain.alignment.hmm_name)
                    domain_accession = pyhmmer_decode(domain.alignment.hmm_accession)
                    match_positions = sum(1 for c in domain.alignment.identity_sequence if c.isalpha())
                    
                    # TODO: maybe add an alignment length and envelope length to the SearchResult? Currently, identity is based on alignment size, but start and end are from the envelope, which can be confusing.
                    # print(domain.alignment.hmm_sequence)
                    # print(domain.alignment.identity_sequence)
                    # print(domain.alignment.target_sequence)
                    # print(domain.alignment.hmm_to - (domain.alignment.hmm_from - 1))
                    # print(domain.env_to - (domain.env_from - 1))

                    identity = 100 * match_positions / (domain.alignment.hmm_to - (domain.alignment.hmm_from - 1))

                    if references[domain_name].description is None: #TODO: add test case for hmmer and protein queries that lack descriptions
                        domain_description = ""
                    else:
                        domain_description = pyhmmer_decode(references[domain_name].description) #domain.alignment.hmm_name.decode()
                    # if domain_description == "":
                    #     domain_description = ""
                    domain_evalue = domain.i_evalue
                    score = domain.score
                    start = domain.env_from - 1 #TODO: double check this for off-by-one
                    end = domain.env_to #TODO: double check this for off-by-one
                    ref_start = domain.alignment.hmm_from #TODO: double check this for off-by-one
                    ref_end = domain.alignment.hmm_to #TODO: double check this for off-by-one
                    try:
                        ref_len = references[domain_name].M
                    except AttributeError:
                        ref_len = len(references[domain_name].sequence)

                    if contig_index not in out:
                        out[contig_index] = dict()
                    if cds_index not in out[contig_index]:
                        out[contig_index][cds_index] = list()
                    out[contig_index][cds_index].append(SearchResult(domain_name, domain_description, domain_accession, domain_evalue, score, start, end, db_name, identity, ref_start, ref_end, ref_len, program))

def foldseek_hits_to_search_results(hits, references, evalue, db_name, min_evalue, program, out):
    #["query","target","qheader","theader","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits"]
    for hit in hits:
        hit_name = hit.query.split(",")
        contig_index = int(hit_name[0])
        cds_index = int(hit_name[1])
        domain_name = hit.target
        domain_accession = ""
        domain_description = hit.theader
        domain_evalue = float(hit.evalue)
        score = float(hit.bits)
        start = int(hit.qstart) - 1 #TODO: double check this for off-by-one
        end = int(hit.qend) #TODO: double check this for off-by-one
        ref_start = int(hit.tstart) #TODO: double check this for off-by-one
        ref_end = int(hit.tend) #TODO: double check this for off-by-one
        identity = float(hit.pident)
        ref_len = int(hit.tlen)

        if contig_index not in out: #TODO: deduplicate, code with hmmer_hits_to_search_results
            out[contig_index] = dict()
        if cds_index not in out[contig_index]:
            out[contig_index][cds_index] = list()
        
        out[contig_index][cds_index].append(SearchResult(domain_name, domain_description, domain_accession, domain_evalue, score, start, end, db_name, identity, ref_start, ref_end, ref_len, program))


def nhmmer_hits_to_search_results(hits, references, evalue, db_name, min_evalue, out):

    for top_hits in hits:
        query_name = pyhmmer_decode(top_hits.query.name)
        reference = references[query_name]
        if reference.description is None:
            query_description = ""
        else:
            query_description = pyhmmer_decode(reference.description)
        for hit in top_hits:
            hit_name = pyhmmer_decode(hit.name).split(",")
            contig_index = int(hit_name[0])
            contig_length = int(hit_name[2]) if len(hit_name) > 2 else None
            for domain in hit.domains:
                if domain.i_evalue < evalue and domain.i_evalue >= min_evalue:
                    identity_sequence = getattr(domain.alignment, "identity_sequence", "")
                    match_positions = sum(1 for c in identity_sequence if c.isalpha())
                    alignment_length = domain.alignment.hmm_to - (domain.alignment.hmm_from - 1)
                    identity = 0.0
                    if alignment_length > 0:
                        identity = 100 * match_positions / alignment_length
                    if contig_index not in out:
                        out[contig_index] = dict()
                    if -1 not in out[contig_index]:
                        out[contig_index][-1] = list()
                    try:
                        reference_length = reference.M
                    except AttributeError:
                        reference_length = len(reference.sequence)
                    start = domain.env_from - 1
                    end = domain.env_to
                    if contig_length is not None and start >= contig_length:
                        continue
                    out[contig_index][-1].append(
                        SearchResult(
                            query_name,
                            query_description,
                            "",
                            domain.i_evalue,
                            domain.score,
                            start,
                            end,
                            db_name,
                            identity,
                            domain.alignment.hmm_from,
                            domain.alignment.hmm_to,
                            reference_length,
                            "nhmmer",
                            1,
                            "nucleic_acid",
                        )
                    )


def infernal_hits_to_search_results(hits, references, evalue, db_name, min_evalue, out):

    for top_hits in hits:
        query_name = pyhmmer_decode(top_hits.query.name)
        query_description = pyhmmer_decode(top_hits.query.description) if top_hits.query.description is not None else ""
        query_accession = pyhmmer_decode(top_hits.query.accession) if top_hits.query.accession is not None else ""
        reference = references[query_name]
        for hit in top_hits:
            if hit.evalue < evalue and hit.evalue >= min_evalue:
                hit_name = pyhmmer_decode(hit.name).split(",")
                contig_index = int(hit_name[0])
                contig_length = int(hit_name[2]) if len(hit_name) > 2 else None
                if contig_index not in out:
                    out[contig_index] = dict()
                if -1 not in out[contig_index]:
                    out[contig_index][-1] = list()
                target_from = hit.alignment.target_from
                target_to = hit.alignment.target_to
                strand = -1 if hit.strand == "-" else 1
                start = min(target_from, target_to) - 1
                end = max(target_from, target_to)
                if contig_length is not None and start >= contig_length:
                    continue
                out[contig_index][-1].append(
                    SearchResult(
                        query_name,
                        query_description,
                        query_accession,
                        hit.evalue,
                        hit.score,
                        start,
                        end,
                        db_name,
                        0.0,
                        hit.alignment.cm_from,
                        hit.alignment.cm_to,
                        reference.M,
                        "infernal",
                        strand,
                        "nucleic_acid",
                    )
                )

    
    


def read_references(reference_files: Optional[List[str]], foldseek: Optional[List[str]]) -> Dict[str, Dict[str, Union[ Dict[str, pyhmmer.plan7.HMM], Dict[str, pyhmmer.easel.DigitalSequence], str]]]:
    """
    Args: 
        reference_files: list of str
            a list of hmm or fasta files

    Returns: dict of dicts of dicts of pyhmmer.Plan7.HMM objects (hmmsearch),
                pyhmmer.easel.DigitalSequence objects (phmmer), 
                paths to foldseek databases.
            algorithm: db_name: seq_or_profile_name: HMM or DigitalSequence
            or
            "foldseek": db_name: db_path

    """

    out = {}
    hmm_files = []
    nucleotide_hmm_files = []
    cm_files = []
    fasta_files = []
    nucleotide_fasta_files = []
    if reference_files is not None:
        for file in reference_files:
            file_type = get_file_type(file)
            if file_type not in {"fasta", "hmm", "cm"}:
                raise RuntimeError("Reference files must be fasta, hmmer, or infernal cm.")
            if file_type == "hmm":
                if _hmm_is_nucleic_acid(file):
                    nucleotide_hmm_files.append(file)
                else:
                    hmm_files.append(file)
            elif file_type == "cm":
                cm_files.append(file)
            else:
                if _fasta_is_nucleic_acid(file):
                    nucleotide_fasta_files.append(file)
                else:
                    fasta_files.append(file)
    
    if len(hmm_files) > 0:
        out["hmmsearch"] = read_hmms(hmm_files)
    if len(fasta_files) > 0:
        out["phmmer"] = read_pyhmmer_peptide_fastas(fasta_files)
    if len(nucleotide_hmm_files) > 0:
        out["nhmmer"] = read_hmms(nucleotide_hmm_files)
    if len(nucleotide_fasta_files) > 0:
        out["nhmmer"] = {**out.get("nhmmer", {}), **read_pyhmmer_fastas(nucleotide_fasta_files)}
    if len(cm_files) > 0:
        out["infernal"] = read_infernal_cms(cm_files)

    if foldseek is not None:
        out["foldseek"] = dict()
        for file in foldseek:
            out["foldseek"][Path(file).stem] = file
    
    return out

def run_search(proteins, nucleic_acids, infernal_nucleic_acids, foldseek, references, evalue:int, cpu:int, z=None, min_evalue:float=0.0, max_mode:bool=False):
    """
    Args:
        proteins: list
            a list of pyhmmer.easel.DigitalSequence objects

        foldseek: list
            a list of foldseek foldseek fasta strings

        references: dict of dicts of dicts of pyhmmer.Plan7.HMM objects (hmmsearch) or pyhmmer.easel.DigitalSequence objects (phmmer).
             algorithm: db_name: seq_or_profile_name: HMM or DigitalSequence
            

        evalue: float
            The threshold E value for the hmmsearch hit to be reported

        cpu: int
            The number of CPU cores to be used to run hmmsearch


        z: if set, then pass to hmmsearch as the -Z parameter

        algorithm: hmmsearch, phmmer, or foldseek depending on what the reference type is.

    Returns: dict of dicts of lists of SearchResults
            contig_index: cds_index: [SearchResult]
    """

    out = dict()
    if len(proteins) == 0 and len(nucleic_acids) == 0 and len(infernal_nucleic_acids) == 0:
        return out


    varargs={"sequences": proteins, "cpus": cpu, "Z": z, "E": evalue, "domE": evalue, "incdomE": evalue}
    if max_mode:
        varargs["bias_filter"] = False
        varargs["F1"] = 1.0
        varargs["F2"] = 1.0
        varargs["F3"] = 1.0

    for (algorithm, db_dict) in references.items():

        for db_name in db_dict:
            #hits = pyhmmer.hmmer.hmmsearch(hmms[db_name].values(), proteins, cpu, Z=z, domE=evalue, incdomE=evalue) #yields top_hits for each query
            if algorithm == "hmmsearch":
                hits = pyhmmer.hmmer.hmmsearch(db_dict[db_name].values(), **varargs)
                hmmer_hits_to_search_results(hits, db_dict[db_name], evalue, db_name, min_evalue, "hmmsearch", out)
            elif algorithm == "phmmer":
                hits = pyhmmer.hmmer.phmmer(db_dict[db_name].values(), **varargs)
                hmmer_hits_to_search_results(hits, db_dict[db_name], evalue, db_name, min_evalue, "phmmer", out)
            elif algorithm == "foldseek":
                hits = foldseek_lib.search(db_dict[db_name], proteins = proteins, foldseek = foldseek, cpu = cpu, E = evalue)
                foldseek_hits_to_search_results(hits, db_dict[db_name], evalue, db_name, min_evalue, "foldseek", out)
            elif algorithm == "nhmmer":
                hits = pyhmmer.hmmer.nhmmer(db_dict[db_name].values(), nucleic_acids, cpus=cpu, E=evalue, incE=evalue)
                nhmmer_hits_to_search_results(hits, db_dict[db_name], evalue, db_name, min_evalue, out)
            elif algorithm == "infernal":
                infernal_sequences = pyhmmer.easel.DigitalSequenceBlock(pyhmmer.easel.Alphabet.rna(), infernal_nucleic_acids)
                hits = pyinfernal.cmsearch(db_dict[db_name].values(), infernal_sequences, cpus=cpu, parallel="queries")
                infernal_hits_to_search_results(hits, db_dict[db_name], evalue, db_name, min_evalue, out)
    return out


def clean_rec(rec, clear_best_hit=False, clear_domainator_annotations=False, clear_CDS_annotations=False):
    """
        adds translation and normalized cds_id to SeqRecord objects. Removes best_hit annotations if they exist
        
        rec: Domainator SeqRecord
        
        clear_best_hit: if True, then clear any domainator.DOMAIN_SEARCH_BEST_HIT_NAME annotations

        clear_domainator_annotations: if True, then clear any domainator annotations

        clear_CDS_annotations: if True, then clear any CDS annotations
    """
    new_list = list()
    CDS_counter = 0

    for feature in rec.features:

        if clear_best_hit and feature.type == DOMAIN_SEARCH_BEST_HIT_NAME: # TODO: it would be nice to move this into the genbank parser, like by adding an option to skip features with certain names
            continue
        elif feature.type == DOMAIN_FEATURE_NAME and clear_domainator_annotations:
            continue
        elif feature.type == 'CDS' and clear_CDS_annotations:
            continue

        # add translations and normalized cds_id
        if rec.annotations['molecule_type'] != "protein":
            if feature.type == 'CDS' and "pseudo" not in feature.qualifiers and "pseudogene" not in feature.qualifiers:
                feature.qualifiers['cds_id'] = [get_cds_unique_name(feature)]
                CDS_counter += 1
                if ('translation' not in feature.qualifiers) or (len(feature.qualifiers["translation"][0]) == 0):
                    feature.qualifiers['translation'] = [feature.translate(rec.seq,cds=False)] # cds = False, because we don't want to throw exceptions on weird annotations.

        
        

        new_list.append(feature)

    
    rec.features = new_list
    return CDS_counter

def get_prot_list(contig, unique_id):
    """ 

    Gets a dictionary of every protein to be annotated from the contig
    if a peptide sequence is passed in, then return the peptide sequence
    if a DNA sequence is passed in then extract the CDS translations.

    Args:
        contig: a SeqRecord
        unique_id: a unique id for the contig

    Returns:
        iterator of (name, protein_sequence)
        sequence names are contig_index, cds_index

    """
    #prot_list = list()
    i = unique_id
    if contig.annotations['molecule_type'] == "protein":
        j = 0 #0 because there is only one
        prot = ''.join(filter(str.isalnum, str(contig.seq)))
        if len(prot) > MAX_PROTEIN_SIZE:
            warnings.warn(f"Skipping protein longer than {MAX_PROTEIN_SIZE} aa: {contig.id}")
        else:
            yield (f"{i},{j}", prot)
    else:
        for j, feature in enumerate(contig.features):
            # skip over not cds features
            if feature.type == 'CDS' and "pseudo" not in feature.qualifiers and "pseudogene" not in feature.qualifiers:
                prot = ''.join(filter(str.isalnum, feature.qualifiers['translation'][0]))
                if len(prot) > MAX_PROTEIN_SIZE:
                    warnings.warn(f"Skipping protein longer than {MAX_PROTEIN_SIZE} aa: {contig.id}")
                else:
                    yield (f"{i},{j}", prot)
    
    #return prot_list

def get_pyhmmer_digital_sequence(name, prot):
    return pyhmmer.easel.TextSequence(name=bytes(name, encoding='utf8'),sequence=prot).digitize(pyhmmer.easel.Alphabet.amino())


def get_pyhmmer_digital_nucleotide_sequence(name, sequence, alphabet=None):
    if alphabet is None:
        alphabet = pyhmmer.easel.Alphabet.dna()
    return pyhmmer.easel.TextSequence(name=bytes(name, encoding='utf8'), sequence=sequence).digitize(alphabet)


def circular_hit_segments(hit: SearchResult, contig_length: int):
    if hit.end <= contig_length:
        return [(hit.start, hit.end)]
    return [(hit.start, contig_length), (0, hit.end - contig_length)]


def circular_hit_contains(container: SearchResult, contained: SearchResult, contig_length: int) -> bool:
    container_segments = circular_hit_segments(container, contig_length)
    for contained_segment in circular_hit_segments(contained, contig_length):
        if not any(segment[0] <= contained_segment[0] and segment[1] >= contained_segment[1] for segment in container_segments):
            return False
    return True


def deduplicate_contig_hits(contig: SeqRecord, hits_list: List[SearchResult]) -> List[SearchResult]:
    if contig.annotations.get("topology") != "circular":
        return hits_list
    contig_length = len(contig)
    kept = list()
    for hit in sorted(hits_list, key=lambda x: x.score, reverse=True):
        if any(
            existing.name == hit.name and existing.database == hit.database and existing.program == hit.program and circular_hit_contains(existing, hit, contig_length)
            for existing in kept
        ):
            continue
        kept.append(hit)
    return kept


def build_contig_hit_location(contig: SeqRecord, hit: SearchResult):
    if contig.annotations.get("topology") == "circular" and hit.end > len(contig):
        return CompoundLocation([
            FeatureLocation(hit.start, len(contig), strand=hit.strand),
            FeatureLocation(0, hit.end - len(contig), strand=hit.strand),
        ])
    return FeatureLocation(hit.start, min(hit.end, len(contig)), strand=hit.strand)

def filter_by_overlap(hits, allowed_percent_overlap, presorted=False, by_db=False):
    """
    Removes domains that overlap another domain with a higher bitscore at a percentage
    higher than that of the allowed percent overlap

    Args:
        hits: list of SearchResult
            All of the hits and the information on them for a domain

        allowed_percent_overlap: float
            The highest percentage overlap which will be allowed between domains

    Returns:
        a list of tuples of the hits of a coding sequence with the non-overlapping domains 
        with the highest bitscores
    """
    #TODO: filter different dbs independently.
    no_overlap_hits = []
    if not presorted:
        hits.sort(key=lambda x: x.score, reverse=True)  # sort by bitscore
    
    no_overlap_hits = dict()

    
    for i in range(0, len(hits)):
        if by_db:
            i_db = hits[i].database
        else:
            i_db = None
        
        # Add the first hit as this is the one with the highest bitscore
        if i_db not in no_overlap_hits:
            no_overlap_hits[i_db] = [hits[i]]
            continue
           # Number of non-overlapping hits the current hit does not overlap with
        
        non_overlap_count = 0
        for j in range(0, len(no_overlap_hits[i_db])):
            if not utils.regions_overlap((no_overlap_hits[i_db][j].start, no_overlap_hits[i_db][j].end), (hits[i].start, hits[i].end), allowed_percent_overlap):
                non_overlap_count += 1

        # If it does not overlap all of the kept hits then add it to the list of kept hists
        if non_overlap_count == len(no_overlap_hits[i_db]):
            no_overlap_hits[i_db].append(hits[i])
    # flatten no_overlap_hits
    no_overlap_hits = [hit for sublist in no_overlap_hits.values() for hit in sublist]
    return no_overlap_hits


def add_protein_annotations(contig:SeqRecord, hits_list:List[SearchResult] , max_hits:int, max_overlap:float, no_annotations:bool, best_annotation:bool, overlap_by_db:bool=False):
    """add annoations to a protein contig

    Args:
        contig (SeqRecord): SeqRecord to add annotations to
        hits_list (List[SearchResult]): List of SearchResult objects
        max_hits (int): the maximum number of annotations to add to the contig
        max_overlap (float): the maximum fractional of overlap between domains to be included in the annotated contig. If >= 1, then no overlap filtering will be done.
        no_annotations (bool): if true, then no annotations will be added to the contig
        best_annotation (bool): if True then add a special Domainator annotation to each peptide or CDS with the best hit from this search.
        overlap_by_db (bool): if True, then filter by overlap by database. If False, then filter by overlap by all databases.
    """
    hits_list.sort(key=lambda x: x.score, reverse=True)
    contig.set_hit_best_score(hits_list[0].score)
    contig.set_hit_names({0: hits_list[0].name}) # 0 is a dummy variable, since there is only one CDS in a protein contig.
    if best_annotation:
        hit = hits_list[0]
        location = FeatureLocation(hit.start, hit.end, strand=1)
        #namedtuple("SearchResult", ["name", "desc", "evalue", "score", "start", "end", "database"])
        desc = hit.desc
        if desc.strip() == "":
            desc = "."
        acc = hit.acc
        if acc.strip() == "":
            acc = "."
        f = SeqFeature(location=location, type=DOMAIN_SEARCH_BEST_HIT_NAME,
                        qualifiers={'program': [hit.program], "database":[hit.database], 'description': [desc], 'accession': [acc], 'evalue': [f"{hit.evalue:.1e}"],
                                    'score': [f"{hit.score:.1f}"], 'name':[hit.name], 'identity': [f"{hit.identity:.1f}"], 'cds_id': [f'0_1_{len(contig)}'],
                                    'rstart': [f"{hit.rstart}"], 'rend': [f"{hit.rend}"], 'rlen': [f"{hit.rlen}"]})
        contig.features.append(f)
    if not no_annotations:
        for hit in filter_by_overlap(hits_list[:max_hits], max_overlap,  presorted=True, by_db=overlap_by_db):
            # Create new domain feature with information on the HMMER hit
            location = FeatureLocation(hit.start, hit.end, strand=1)
            #namedtuple("SearchResult", ["name", "desc", "evalue", "score", "start", "end", "database"])
            desc = hit.desc
            if desc.strip() == "":
                desc = "."
            acc = hit.acc
            if acc.strip() == "":
                acc = "."
            f = SeqFeature(location=location, type=DOMAIN_FEATURE_NAME,
                            qualifiers={'program': [hit.program], "database":[hit.database], 'description': [desc], 'accession': [acc], 'evalue': [f"{hit.evalue:.1e}"], 'score': [f"{hit.score:.1f}"],
                                        'name':[hit.name], 'identity': [f"{hit.identity:.1f}"], 'cds_id': [f'0_1_{len(contig)}'], 'rstart': [f"{hit.rstart}"],
                                        'rend': [f"{hit.rend}"], 'rlen': [f"{hit.rlen}"]})
            contig.features.append(f)

#contig.features, cds_index, seq_group_index, contig_index, db_index, database_names, hits, max_hits, max_overlap
#contig, contigs_list[contig_id], max_hits, max_overlap
def add_nucleic_acid_annotations(contig:SeqRecord, hits:Dict[int,List[SearchResult]] , max_hits:int, max_overlap:float, no_annotations:bool, best_annotation:bool, overlap_by_db:bool=False):
    """add annoations to a nucleic acid contig

    Args:
        contig (SeqRecord): SeqRecord to add annotations to
        hits_list (Dict[List[SearchResult]]): Keys are cds_index
        max_hits (int): the maximum number of annotations to add to the contig
        max_overlap (float): the maximum fractional of overlap between domains to be included in the annotated contig. If >= 1, then no overlap filtering will be done.
        no_annotations (bool): if true, then no annotations will be added to the contig
        best_annotation (bool): if True then add a special Domainator annotation to each peptide or CDS with the best hit from this search.
        overlap_by_db (bool): if True, then filter by overlap by database. If False, then filter by overlap by all databases.
    """
    hit_scores = dict()
    best_hits = dict()
    for feature_id in hits: #iterate through the CDS features
        feature = contig.features[feature_id]
        hits[feature_id].sort(key=lambda x: x.score, reverse=True)
        hit_scores[feature_id] = hits[feature_id][0].score #TODO: this might be inconsistent if max_hits is 0, but it would be weird for max_hits to be 0.
        best_hits[feature_id] = hits[feature_id][0].name
        if best_annotation:
            hit = hits[feature_id][0]
            annot_length = (hit.end - hit.start) * 3
            try:
                location = feature.location.overlay(hit.start*3, annot_length)
            except:
                warnings.warn(f"Could not overlay location for {hit.name} on {contig.name}, {feature.qualifiers['cds_id'][0]}. Skipping annotation.")
                continue
            desc = hit.desc
            if desc.strip() == "":
                desc = "."
            acc = hit.acc
            if acc.strip() == "":
                acc = "."
            f = SeqFeature(location=location, type=DOMAIN_SEARCH_BEST_HIT_NAME, 
                            qualifiers={'program': [hit.program], "database":[hit.database], 'description': [desc], 'accession': [acc], 'evalue': [f"{hit.evalue:.1e}"], 'score': [f"{hit.score:.1f}"], 
                                        'name':[hit.name], 'identity': [f"{hit.identity:.1f}"], 'cds_id': [feature.qualifiers['cds_id'][0]], 'rstart': [f"{hit.rstart}"],
                                        'rend': [f"{hit.rend}"], 'rlen': [f"{hit.rlen}"]})

            contig.features.append(f)
        
        if not no_annotations:
            for hit in filter_by_overlap(hits[feature_id][:max_hits], max_overlap, presorted=True, by_db=overlap_by_db):
                # Create new domain feature with information on the HMMER hit
                annot_length = (hit.end - hit.start) * 3
                try:
                    location = feature.location.overlay(hit.start*3, annot_length)
                except:
                    warnings.warn(f"Could not overlay location for {hit.name} on {contig.name}, {feature.qualifiers['cds_id'][0]}. Skipping annotation.")
                    continue

                desc = hit.desc
                if desc.strip() == "":
                    desc = "."

                acc = hit.acc
                if acc.strip() == "":
                    acc = "."
                f = SeqFeature(location=location, type=DOMAIN_FEATURE_NAME, 
                                qualifiers={'program': [hit.program], "database":[hit.database], 'description': [desc], 'accession': [acc], 'evalue': [f"{hit.evalue:.1e}"], 'score': [f"{hit.score:.1f}"], 
                                            'name':[hit.name], 'identity': [f"{hit.identity:.1f}"], 'cds_id': [feature.qualifiers['cds_id'][0]], 'rstart': [f"{hit.rstart}"], 
                                            'rend': [f"{hit.rend}"], 'rlen': [f"{hit.rlen}"]})

                contig.features.append(f)
    
    # makes a temporary annotation in the SeqRecord containing a {cds_index: best_domain_score}
    contig.set_hit_scores(hit_scores)
    contig.set_hit_names(best_hits)


def add_contig_nucleic_acid_annotations(contig: SeqRecord, hits_list: List[SearchResult], max_hits: int, max_overlap: float, no_annotations: bool, best_annotation: bool, overlap_by_db: bool = False):
    hits_list = deduplicate_contig_hits(contig, hits_list)
    hits_list.sort(key=lambda x: x.score, reverse=True)
    contig.set_hit_best_score(hits_list[0].score)
    hit_scores = dict()
    hit_names = dict()

    if best_annotation:
        hit = hits_list[0]
        desc = hit.desc if hit.desc.strip() else "."
        acc = hit.acc if hit.acc.strip() else "."
        feature_index = len(contig.features)
        contig.features.append(
            SeqFeature(
                location=build_contig_hit_location(contig, hit),
                type=DOMAIN_SEARCH_BEST_HIT_NAME,
                qualifiers={
                    'program': [hit.program],
                    'database': [hit.database],
                    'description': [desc],
                    'accession': [acc],
                    'evalue': [f"{hit.evalue:.1e}"],
                    'score': [f"{hit.score:.1f}"],
                    'name': [hit.name],
                    'identity': [f"{hit.identity:.1f}"],
                    'cds_id': ['contig'],
                    'rstart': [f"{hit.rstart}"],
                    'rend': [f"{hit.rend}"],
                    'rlen': [f"{hit.rlen}"],
                    'target_kind': [hit.target_kind],
                },
            )
        )
        hit_scores[feature_index] = hit.score
        hit_names[feature_index] = hit.name

    contig.set_hit_scores(hit_scores)
    contig.set_hit_names(hit_names)

    if no_annotations:
        return

    for hit in filter_by_overlap(hits_list[:max_hits], max_overlap, presorted=True, by_db=overlap_by_db):
        desc = hit.desc if hit.desc.strip() else "."
        acc = hit.acc if hit.acc.strip() else "."
        contig.features.append(
            SeqFeature(
                location=build_contig_hit_location(contig, hit),
                type=DOMAIN_FEATURE_NAME,
                qualifiers={
                    'program': [hit.program],
                    'database': [hit.database],
                    'description': [desc],
                    'accession': [acc],
                    'evalue': [f"{hit.evalue:.1e}"],
                    'score': [f"{hit.score:.1f}"],
                    'name': [hit.name],
                    'identity': [f"{hit.identity:.1f}"],
                    'cds_id': ['contig'],
                    'rstart': [f"{hit.rstart}"],
                    'rend': [f"{hit.rend}"],
                    'rlen': [f"{hit.rlen}"],
                    'target_kind': [hit.target_kind],
                },
            )
        )
        

def domainator_inner(contigs_list, proteins_list, nucleic_acid_list, infernal_nucleic_acid_list, foldseek_list, reference_groups, evalue, cpu, z, hits_only, no_annotations, max_hits, max_overlap, best_annotation, min_evalue=0.0, max_mode=False, overlap_by_db=False):
    """

    Args:
        contigs_list (_type_): 
        proteins_list (_type_): 
        foldseek_list (_type_):
        reference_groups (_type_): 
        evalue (_type_): 
        cpu (_type_): 
        z (_type_): 
        hits_only (_type_): 
        no_annotations (_type_): 
        max_hits (_type_): 
        max_overlap (_type_): 
        best_annotation (_type_): 
        min_evalue (float, optional): . Defaults to 0.0.
        max_mode (bool, optional): . Defaults to False.
        overlap_by_db (bool, optional): Defaults to False.

    Returns:
        : 
    """
    hits =  run_search(proteins_list, nucleic_acid_list, infernal_nucleic_acid_list, foldseek_list, reference_groups, evalue, cpu, z, min_evalue=min_evalue, max_mode=max_mode) # returns a dict of dicts of lists of SearchResults, where keys are contig_index, cds_index
    
    for contig_id in hits:
        contig = contigs_list[contig_id]
        if contig.annotations['molecule_type'] == "protein":
            add_protein_annotations(contig, hits[contig_id][0], max_hits, max_overlap, no_annotations, best_annotation, overlap_by_db=overlap_by_db)
        elif -1 in hits[contig_id]:
            add_contig_nucleic_acid_annotations(contig, hits[contig_id][-1], max_hits, max_overlap, no_annotations, best_annotation, overlap_by_db=overlap_by_db)
        else:
            add_nucleic_acid_annotations(contig, hits[contig_id], max_hits, max_overlap, no_annotations, best_annotation, overlap_by_db=overlap_by_db)

    if hits_only:
        return_contigs = list(hits.keys())
    else:
        return_contigs = list(range(len(contigs_list)))

    return return_contigs

def prodigal_CDS_annotate(rec:SeqRecord):
    """Annotate a SeqRecord with CDS features using prodigal

    Args:
        rec (SeqRecord): SeqRecord to annotate

    """
    orf_finder = pyrodigal.GeneFinder(meta=True) #TODO: maybe make this static?

    i = 1
    for pred in orf_finder.find_genes(bytes(rec.seq)):
        
        if pred.partial_begin or pred.partial_end:
            continue
        translation = pred.translate()
        strand = pred.strand
        start = pred.begin - 1
        end = pred.end
        feature = SeqFeature(location=FeatureLocation(start,end,strand), type="CDS", qualifiers={"translation":[translation], "gene_id":[f"{rec.id}_{i}"]})
        feature.qualifiers['cds_id'] = [get_cds_unique_name(feature)]
        rec.features.append(feature)
        i += 1

def domainate(seq_iterator, references, z, evalue=10, max_hits=sys.maxsize, max_overlap=1, cpu=0,  batch_size=10000, hits_only=False, no_annotations=False, pre_parsed_references=None, best_annotation=False, gene_call=None, min_evalue=0.0, ncbi_taxonomy=None, include_taxids=None, exclude_taxids=None, max_mode=False, foldseek=None, esm2_3Di_weights=None, esm2_3Di_device=None, overlap_by_db=False):
    """
    The main function of the hmmer domain annotation algorithm

    yields annotated SeqRecords

    Args:
        sequence_groups: an iterator yielding SeqRecord objects
            
        references: list of str, path-like
            A list of hmm or protein fasta file paths, the input HMM databases to be searched against the sequences using hmmsearch or phmmer.
            If you pass None to references, you can supply your own pre-parsed references via pre_parsed_references.

        z: Z parameter to pass to hmmsearch

        evalue: float
            The threshold E value for hmmer hit

        max_hits: int
            The maximum number of hits per CDS to be reported

        max_overlap: float
            The maximum fractional of overlap to be allowed between domains

        cpu: int
            The number of CPU cores that should be used to run hmmsearch
        
        batch_size: how many protein sequences at a time to pass to pymmer hmmsearch

        hits_only: if True, then only output contigs with at least one domain annotation

        no_annotations: if True, then don't add any new annotations to the sequences, just return the hits with existing annotations.

        pre_parsed_references: dict, the output of read_references

        best_annotation: if True then add a special Domain_search annotation to each peptide or CDS with the best hit from this search.

        gene_call: If 'all', then any existing CDS annotations will be deleted and all contigs will be re-annotated. If 'unannotated', then only contigs without CDS annotations will be annotated. [default: keep existing annotations]

        min_evalue: only report hits with evalue >= min_evalue. [default: 0.0]

        ncbi_taxonomy: an NCBITaxonomy object (needed for filtering by taxonomy, otherwise not needed)
        
        include_taxids: only keep contigs with a taxonomy id in this list
        
        exclude_taxids: only keep contigs with a taxonomy id not in this list

        max_mode: if True, then run hmmsearch/phmmer in maximum sensitivity mode, which is much slower, but more sensitive.

        foldseek: paths to foldseek database files. [default: None]

        esm2_3Di_weights: path to esm2_3Di checkpoint file. [default: None]

        esm2_3Di_device: device to run esm2_3Di on. [default: None]

    """

    if references is None and foldseek is None:
        reference_groups = pre_parsed_references
    else:
        reference_groups = read_references(references, foldseek=foldseek)

    if "foldseek" in reference_groups:
        foldseek_builder = foldseek_lib.foldseekBuilder(esm2_3Di_device, esm2_3Di_weights)

    if gene_call not in (None, 'all', 'unannotated'):
        raise ValueError("gene_call must be one of None, 'all', or 'unannotated'")
    
    clear_domainator_annotations = False
    clear_CDS_annotations = False
    if gene_call == 'all':
        clear_domainator_annotations = True
        clear_CDS_annotations = True

    if include_taxids is not None:
        include_taxids = set(include_taxids)
    if exclude_taxids is not None:
        exclude_taxids = set(exclude_taxids)
    
    contigs_list = list()
    proteins_list = list()
    nucleic_acid_list = list()
    infernal_nucleic_acid_list = list()
    foldseek_list = list()
    contig_index = 0
    if include_taxids or exclude_taxids:
        seq_iterator = filter_by_taxonomy(seq_iterator, include_taxids, exclude_taxids, ncbi_taxonomy)

    nhmmer_extension = get_max_reference_length(reference_groups, "nhmmer")
    infernal_extension = get_max_reference_length(reference_groups, "infernal")
    
    for rec in seq_iterator:
        CDS_count = clean_rec(rec, best_annotation, clear_domainator_annotations, clear_CDS_annotations)
        if CDS_count == 0 and rec.annotations['molecule_type'] != "protein" and gene_call is not None:
            prodigal_CDS_annotate(rec)
        contigs_list.append(rec)
        if "foldseek" in reference_groups:
            foldseek_list.extend([foldseek_builder(name, prot) for name, prot in get_prot_list(rec, contig_index)])
        if "hmmsearch" in reference_groups or "phmmer" in reference_groups or "foldseek" in reference_groups:
            proteins_list.extend([get_pyhmmer_digital_sequence(name, prot) for (name, prot) in get_prot_list(rec, contig_index)])
        if "nhmmer" in reference_groups and rec.annotations['molecule_type'] != "protein":
            nucleic_acid_list.append(get_pyhmmer_digital_nucleotide_sequence(f"{contig_index},contig,{len(rec)}", get_contig_search_sequence(rec, nhmmer_extension)))
        if "infernal" in reference_groups and rec.annotations['molecule_type'] != "protein":
            infernal_nucleic_acid_list.append(get_pyhmmer_digital_nucleotide_sequence(f"{contig_index},contig,{len(rec)}", get_contig_search_sequence(rec, infernal_extension), pyhmmer.easel.Alphabet.rna()))
        contig_index += 1
        if len(proteins_list) >= batch_size:
            return_contigs = domainator_inner(contigs_list, proteins_list, nucleic_acid_list, infernal_nucleic_acid_list, foldseek_list, reference_groups, evalue, cpu, z, hits_only, no_annotations, max_hits, max_overlap, best_annotation, min_evalue=min_evalue, max_mode=max_mode, overlap_by_db=overlap_by_db)
            for contig_id in range(len(contigs_list)):
                if contig_id in return_contigs:
                    yield contigs_list[contig_id]
            contigs_list = list()
            proteins_list = list()
            nucleic_acid_list = list()
            infernal_nucleic_acid_list = list()
            foldseek_list = list()
            contig_index = 0
    
    return_contigs = domainator_inner(contigs_list, proteins_list, nucleic_acid_list, infernal_nucleic_acid_list, foldseek_list, reference_groups, evalue, cpu, z, hits_only, no_annotations, max_hits, max_overlap, best_annotation, min_evalue=min_evalue, max_mode=max_mode, overlap_by_db=overlap_by_db)
    #for contig_id in range(len(contigs_list)): #TODO: why did I think this more complicated loop was a good idea?
    #    if contig_id in return_contigs:
    for contig_id in return_contigs:
        yield contigs_list[contig_id]

def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)


    parser.add_argument('-i', '--input', default=None, required=True,
                       nargs='+', type=str,
                       help="the genbank or fasta files to annotate. Genbank files can be nucleotide (with CDS annotations) or peptide. Fasta files can be protein or nucleotide; use --fasta_type to specify which.")
    parser.add_argument("--fasta_type", type=str, default="protein", choices={"protein", "nucleotide"}, 
                        help="Whether the sequences in fasta files are protein or nucleotide sequences.")
    
    # TODO: might be better to allow user to have different command line arguments for different kinds of databases? Particularly for precompiled databases like hhsuite and foldseek.
    parser.add_argument('-r', '--references', required=False, type=str,
                        default=None, nargs='+',
                        help="reference files to search against the input. Supported types are protein HMMs, protein FASTA, nucleotide FASTA, nucleotide HMMs, and infernal CM files.")
    
    parser.add_argument("--foldseek", nargs="+", default=None, type=str, required=False,
                        help="Foldseek database files.")
    parser.add_argument("--esm2_3Di_weights", default=None, type=str, required=False,
                        help="checkpoint file for esm2_3Di model.") # TODO: this should be automated.
    parser.add_argument("--esm2_3Di_device", default="cuda:0", type=str, required=False,
                        help="device to use for esm2_3Di model.")
    
    parser.add_argument('-o', '--output', default=None, type=str, required=False,
                        help="output genbank filename. If not supplied, writes to stdout.")


    parser.add_argument('-Z', default=1000, type=int, 
                        help="Passed as the -Z parameter to hmmsearch: Assert that the total number of targets in your searches is <x>, for the purposes of per-sequence E-value calculations, rather than the actual number of targets seen. "
                        "Default: 1000\n"
                        "Set to 0 to use the number actual target sequences (currently this does does not account for taxonomic filtering)." # TODO: account for taxonomy filtering
                        "Supplying -Z is recommended for most use cases of domainate, as it will speed up the search and make comparisons between searches more meaningful.")

    parser.add_argument('-e', '--evalue', type=float, default=0.001, 
                        help="threshold E value for the domain hit. Must be <=10. [default 0.001]")
    
    parser.add_argument('--min_evalue', type=float, default=0.0,
                        help="hits with E value LOWER (as in BETTER) than this will be filtered out. NOT FREQUENTLY USED. Use only if you want to eliminate close matches. Must be >=0. [default 0]")

    parser.add_argument('--max_domains', type=int, default=0,
                        help="the maximum number of domains to be reported per CDS. If not specified, then no max_domains filter is applied.")
    
    parser.add_argument('--max_mode', action="store_true", default=False,
                        help="Run hmmsearch, phmmer, or nhmmer in maximum sensitivity mode, which is much slower, but more sensitive.")

    parser.add_argument("--include_taxids", nargs='+', default=None, type=int, help="Space separated list of taxids to include. Contigs with taxonomy not in this list will be skipped.")
    parser.add_argument("--exclude_taxids", nargs='+', default=None, type=int, help="Space separated list of taxids to exclude. Contigs with taxonomy in this list will be skipped.")
    parser.add_argument("--ncbi_taxonomy_path", type=str,  default="/tmp/ncbi_taxonomy", help="Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist.")
    parser.add_argument("--taxonomy_update", action="store_true", help="If taxonomy database exists, check it against the version on the ncbi server and update if there is a newer version.")

    parser.add_argument('--cpu', type=int, default=0,
                        help="the number of cores of the cpu which are used at a time to run the search [default: use all available cores]")
    parser.add_argument('--max_overlap', type=float, default=1,
                        help="the maximum fractional of overlap between domains to be included in the annotated genbank. If >= 1, then no overlap filtering will be done.")
    parser.add_argument("--overlap_by_db", action='store_true', default=False,
                        help="If activated, then overlap filtering will be done by database, rather than all together.")
    parser.add_argument('--hits_only', action='store_true', default=False,
                        help="when activated, the ouptut will only have contigs with at least one domain annotation. In many cases, domain_search.py will be faster and more appropriate.")
    parser.add_argument('--no_annotations', action='store_true', default=False,
                        help="when activated, new annotations will not be added to the file. This is useful if you are just trying to extract a set of contigs for annotation in a later step.")

    parser.add_argument("--batch_size", type=int, default=10000, required=False,
                        help="How many target sequences to search at one time in a batch. Increasing this number might improve speeds at the cost of memory.")
    #TODO: maybe add an hmm_batch_size parameter in case there are really huge hmm files that don't all fit in memory?
 
    parser.add_argument("--offset", type=int, default=None, required=False,
                       help="File offset to start reading records at. (Note that this usually only makes sense when there is just a single input file). Default: start at the top of the file")
    parser.add_argument("--recs_to_read", type=int, default=None, required=False,
                       help="Stop after reading this many records. Default: read all the records")

    parser.add_argument('--gene_call', type=str, default=None, choices = {"all", "unannotated"}, required=False,
                        help="When activated, new CDS annotations will be added with Prodigal in Metagenomic mode. If 'all', then any existing CDS annotations will be deleted and all contigs will be re-annotated. If 'unannotated', then only contigs without CDS annotations will be annotated. [default: None]. "
                        "If you supply a nucleotide fasta file, be sure to also specify --fasta_type nucleotide. "
                        "Note that if you do gene calling, it is STRONGLY recommended that you also supply -Z, because database size is pre-calcuated at the beginning of the execution, whereas gene-calling is done on the fly. Not supplying Z may become an error in the future.")

    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)

    if params.evalue >10:
        raise ValueError("evalue must be <= 10.")
    if params.min_evalue < 0:
        raise ValueError("min_evalue must be >= 0.")

    if params.references is None and params.foldseek is None:
        raise ValueError("You must specify at least one reference database.")

    reference_groups = read_references(params.references, foldseek=params.foldseek)
    molecule_types = get_input_molecule_types(params.input, default_molecule_type=params.fasta_type)
    validate_input_molecule_types(molecule_types, reference_groups)

    ## figure out if we are using every contig or just certain, named contigs ##
    #contigs_needed = list_and_file_to_dict_keys(params.contigs, params.contigs_file)

    Z = params.Z
    if Z == 0:
        Z = 0
        for file_path in params.input: 
            Z += utils.count_peptides(file_path)

    if params.output is None:
        out = sys.stdout
    else:
        out = open(params.output, "w")

    if params.max_domains == 0:
        max_domains = sys.maxsize
    else:
        max_domains = params.max_domains

    if params.no_annotations and not params.hits_only:
        raise RuntimeError(f"You've specified no_annotations, without specifying hits_only. This seems like a useless thing to do. If you want to select a subset of contigs, please use select_by_contig.py.")
    recs_to_read = float("inf")
    if params.recs_to_read is not None:
        recs_to_read = params.recs_to_read

    ncbi_taxonomy = None
    if params.include_taxids or params.exclude_taxids:
        # create the path to the NCBI taxonomy database
        Path(params.ncbi_taxonomy_path).mkdir(parents=True, exist_ok=True)
        # load the NCBI taxonomy database
        ncbi_taxonomy = NCBITaxonomy(params.ncbi_taxonomy_path, params.taxonomy_update)

    # Run
    write_genbank(
        domainate(
            parse_seqfiles(params.input, None, filetype_override=None, seek_to=params.offset, max_recs=recs_to_read, default_molecule_type=params.fasta_type),
            references=None,
            z=Z,
            evalue=params.evalue,
            min_evalue=params.min_evalue,
            max_hits=max_domains,
            max_overlap=params.max_overlap,
            cpu=params.cpu,
            batch_size=params.batch_size,
            hits_only=params.hits_only,
            no_annotations=params.no_annotations,
            pre_parsed_references=reference_groups,
            gene_call=params.gene_call,
            ncbi_taxonomy=ncbi_taxonomy,
            include_taxids=params.include_taxids,
            exclude_taxids=params.exclude_taxids,
            max_mode=params.max_mode,
            foldseek=None,
            esm2_3Di_weights=params.esm2_3Di_weights,
            esm2_3Di_device=params.esm2_3Di_device,
            overlap_by_db=params.overlap_by_db,
        ),
        out)

    if params.output is not None:
        out.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
