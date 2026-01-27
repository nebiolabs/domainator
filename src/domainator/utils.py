"""Miscellaneous common functions for Domainator

"""
import warnings
warnings.filterwarnings("ignore", module='numpy')
from domainator.Bio import SeqIO, BiopythonParserWarning, BiopythonWarning
from domainator.Bio.SeqRecord import SeqRecord
from domainator.Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, ExactPosition
from domainator.Taxonomy import NCBITaxonomy
from pathlib import Path
from datetime import date
import pyhmmer
import os
import multiprocessing
from array import array

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Set, Union, Tuple

from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from typing import List, Optional
from collections import defaultdict, OrderedDict
import functools
import sys
import seaborn as sns
import pandas as pd

from io import IOBase

EXTENSION_TO_TYPE = {
    "gb":"genbank",
    "gbk":"genbank",
    "genbank":"genbank",
    "gbff":"genbank",
    "fasta":"fasta",
    "fa":"fasta",
    "fna":"fasta",
    "pep":"fasta",
    "faa":"fasta",
    # "gff":"gff",
    # "gtf":"gff",
    # "gff3":"gff",
    "hmmer":"hmm",
    "hmm":"hmm",
    "hdf5":"hdf5",
    "hdf":"hdf5",
    "h5":"hdf5",
}


DEFAULT_BUFFER_SIZE = 5000

_MP_CONTEXT = None


def pyhmmer_decode(value):
    """Decode a pyhmmer attribute that may be bytes or str.
    
    PyHMMER 0.12.0 changed many attributes from bytes to str.
    This helper provides backwards compatibility with both versions.
    
    Args:
        value: A bytes or str value from a pyhmmer object attribute
        
    Returns:
        str: The decoded string value
    """
    if value is None:
        return None
    if isinstance(value, bytes):
        return value.decode()
    return value


def get_multiprocessing_context():
    """Return a cached multiprocessing context that avoids unsafe fork start."""

    global _MP_CONTEXT
    if _MP_CONTEXT is not None:
        return _MP_CONTEXT

    preferred_method = os.environ.get("DOMAINATOR_MP_START_METHOD", "spawn")
    try:
        ctx = multiprocessing.get_context(preferred_method)
    except ValueError:
        warnings.warn(
            f"Requested multiprocessing start method '{preferred_method}' is unavailable; "
            "falling back to Python's default start method.",
            RuntimeWarning,
            stacklevel=2,
        )
        ctx = multiprocessing.get_context()

    _MP_CONTEXT = ctx
    return ctx


def make_pool(*args, **kwargs):
    """Create a multiprocessing Pool using the safe start method."""

    return get_multiprocessing_context().Pool(*args, **kwargs)


def make_manager():
    """Create a multiprocessing Manager using the safe start method."""

    return get_multiprocessing_context().Manager()

warnings.filterwarnings("ignore", category=BiopythonParserWarning)

def parse_simple_list(filename):
	"""
		splits a text file on newlines, returns a set with the stripped lines as entries
        returns a dict where keys are lines from the file and values are None
	"""
	out = dict()
	(infile, input_type) = open_if_is_name(filename)
	
	for l in infile:
		out[l.strip()]  = None
	
	if (input_type == "name"):
		infile.close()
	return out

def read_hmms(hmm_files:Iterable[Union[str,os.PathLike,IOBase]]) -> Dict[str, Dict[str,pyhmmer.plan7.HMM]]:
    """
        hmm_files: a list of paths to .hmm files

        returns:
            a dict of dicts of pyhmmer HMM objects
                db_name: hmm_name: HMM

    """
    out = dict()
    for file in hmm_files:
        if isinstance(file, IOBase):
            name = file.name
        else:
            name = file
        name = os.path.basename(Path(name).stem)

        hmmer_models = OrderedDict() 
        for model in pyhmmer.plan7.HMMFile(file):
            model_name = pyhmmer_decode(model.name)
            if model_name in hmmer_models:
                warnings.warn(f"multiple hmms with the same name ({model_name}) in file: {file}, only one will be used.")
            hmmer_models[model_name] = model
        if name in out:
            raise RuntimeError(f"Multiple hmm files with the same name, please combine the hmms into a single file, or rename one of the files. This is important to avoid searching the same domains twice, and for naming the source databases.")
        out[name] = hmmer_models
    return out


def read_pyhmmer_peptide_fastas(peptide_files):
    out = dict()
    for file_name in peptide_files:
        name = os.path.basename(Path(file_name).stem)
        seqs_dict = dict()
        with pyhmmer.easel.SequenceFile(file_name, digital=True) as seq_file:
            for seq in seq_file:
                seq_name = pyhmmer_decode(seq.name)
                if seq_name in seqs_dict:
                    warnings.warn(f"multiple reference sequences with the same name ({seq_name}) in file: {file_name}, only one will be used.")
                seqs_dict[seq_name] = seq
        if name in out:
            raise RuntimeError(f"Multiple reference sequence files with the same name, please combine the reference sequences into a single file, or rename one of the files. This is important to avoid searching the same domains twice, and for naming the source databases.")
        out[name] = seqs_dict
    return out

def get_file_type(filename):
    out = None
    file_extension = Path(filename).suffix[1:]
    if file_extension in EXTENSION_TO_TYPE:
        out = EXTENSION_TO_TYPE[file_extension]
    return out

def open_if_is_name(filename_or_handle, mode="r"):
    """
        if a file handle is passed, return the file handle
        if a Path object or path string is passed, open and return a file handle to the file.

        returns:
            file_handle, input_type ("name" | "handle")
    """
    out = filename_or_handle
    input_type = "handle"
    try:
        out = open(filename_or_handle, mode)
        input_type = "name"
    except TypeError:
        pass
    except Exception as e:
        raise(e)

    return (out, input_type)


def get_cds_unique_name(feature):
    """
        If the feature already has a cds_id, then keep it, otherwise generate one based on the position on the contig.
    """
    if "cds_id" in feature.qualifiers:
        return feature.qualifiers["cds_id"][0]
    else:
        # need the strand information to account for circular contigs.
        name_parts = ["_".join( (str(p.stranded_start_human_readable), str(p.strand), str(p.stranded_end_human_readable)) ) for p in feature.location.parts]
        return " ".join(name_parts) # space so it can be split into multiple lines when writing genbank files

def get_cds_name(feature, precedence=None): #(contig_id, feature):
    if precedence is None:
        precedence = ["gene_id", "locus_tag", "gene", "protein_id"]

    for p in precedence:
        if p in feature.qualifiers:
            return feature.qualifiers[p][0]
    return get_cds_unique_name(feature)


class TaxonomyData:
    """
    Dataclass for storing pre_computed taxonomy data
    """
    # lineage: List[str]
    # rank: str
    def __init__(self, ncbi_taxonomy:NCBITaxonomy, record:SeqRecord=None, taxid:int=None):
        if record is None and taxid is None:
            raise RuntimeError("Must provide either a SeqRecord or a taxid")
        if record is not None and taxid is not None:
            raise RuntimeError("Must provide either a SeqRecord or a taxid, not both")
        if record is not None:
            self.taxid = get_taxid(record)
        else:
            self.taxid = taxid

        self.taxid = ncbi_taxonomy.normalized_taxid(self.taxid)
        self.lineage = ncbi_taxonomy.lineage(self.taxid)
        self.lineage.reverse() # reverse so that self has the highest index
        self.ranks = [ ncbi_taxonomy.rank(x).lower() for x in self.lineage ]
        self.names = [ ncbi_taxonomy.name(x) for x in self.lineage ]
        self.rank_to_name = dict(zip(self.ranks, self.names))
        self.rank_to_taxid = dict(zip(self.ranks, self.lineage))

    def __str__(self):
        return f"TaxonomyData(taxid={self.taxid}, lineage={self.lineage}, ranks={self.ranks}, names={self.names})"
    def __repr__(self):
        return str(self)
    

@dataclass
class DomainatorCDS():
    name: str
    num: str
    index: int #feature index from the full contig
    feature: SeqFeature
    domain_features: List[SeqFeature] = field(default_factory=lambda: [])
    domain_names: Set[str] = field(default_factory=lambda: set())
    domain_search_feature: SeqFeature = None

    @classmethod
    def list_from_contig(cls, contig, domain_evalue=float("inf"), domain_score=float("-inf"), skip_pseudo=False, name_precedence=None):
        """
            returns a list of DomainatorCDS objects from a contig SeqRecord
        """
        cdss: List[DomainatorCDS] = list()
        name_to_idx: Dict[str,int] = dict()
        idx:int = 0
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
                if float(feature.qualifiers['evalue'][0]) < domain_evalue and float(feature.qualifiers['score'][0]) > domain_score:
                    cds_num = feature.qualifiers['cds_id'][0]
                    if cds_num not in name_to_idx:
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
            if evalue is not None and float(f.qualifiers["evalue"][0]) < evalue:
                keep = False

            if keep == True:    
                out.add(f.qualifiers["name"][0])
        return out
    
        

class BooleanEvaluator():
    def __init__(self, expression):
        # TODO: add support for escaped characters
        tokens = list()
        self.tok_to_domain = dict()
        self.expression_parts = list()
        # tokenize the expression
        elem = ''

        for y in str(expression):
            if(y == '&'):
                if(elem != ''):
                    tokens.append(elem)
                    elem = ''
                tokens.append(y)
            elif(y == '|'):
                if(elem != ''):
                    tokens.append(elem)
                    elem = ''
                tokens.append(y)
            elif(y == '~'):
                tokens.append(y)
            elif(y == '('):
                tokens.append(y)
            elif(y == ')'):
                tokens.append(elem)
                elem = ''
                tokens.append(y)
            elif(y != ' '):
                elem += y
        if(elem != ''):
            tokens.append(elem)


        # translate the expression
        i = 0
        for x in tokens:
            if(x == '(' or x == ')'):
                self.expression_parts.append(x)
            elif(x == '|'): 
                self.expression_parts.append("or")
            elif(x == '&'):
                self.expression_parts.append('and')
            elif(x == '~'):
                self.expression_parts.append('not')
            else:
                self.expression_parts.append(f"tok{i}")
                self.tok_to_domain[f"tok{i}"] = x
                i += 1
        #TODO: check sytax? Maybe just plug in zero for every tok and call eval?

    def check_expression(self, domains):

        domains = set(domains)
        domain_present = dict()
        for tok, domain in self.tok_to_domain.items():
            if domain in domains:
                domain_present[tok] = "1"
            else:
                domain_present[tok] = "0"

        interpolated_expression = list()
        for part in self.expression_parts:
            if part[0:3] == "tok":
                if part in domain_present:
                    interpolated_expression.append(domain_present[part])
                else:
                    interpolated_expression.append(domain_present[part])
            else:
                interpolated_expression.append(part)
        expression = " ".join(interpolated_expression)
        
        return eval(expression)
    
    @classmethod
    def sanitize_identifier(cls, identifier: str) -> str:
        """
        Convert an identifier string to a name safe for boolean expressions.
        
        The BooleanEvaluator uses these special characters:
            & - AND operator
            | - OR operator  
            ~ - NOT operator
            ( - grouping open
            ) - grouping close
            space - token separator
        
        This function replaces these characters to create valid identifiers.
        
        Args:
            identifier: The string to sanitize
        
        Returns:
            A sanitized string suitable for boolean expressions
        """
        # Characters that have special meaning in BooleanEvaluator
        # We use double/triple underscores to avoid collision with single underscores
        # that might already be in the identifier
        replacements = [
            ('(', '__'),    # grouping open
            (')', '___'),   # grouping close
            ('&', '_AND_'), # AND operator
            ('|', '_OR_'),  # OR operator
            ('~', '_NOT_'), # NOT operator
            (' ', '_'),     # space separator
            ('-', '_'),     # common separator in patterns (e.g., PROSITE)
            ('<', 'Nterm_'),  # N-terminal anchor (PROSITE)
            ('>', '_Cterm'),  # C-terminal anchor (PROSITE)
            ('.', ''),      # trailing period (PROSITE)
        ]
        
        result = identifier
        for old, new in replacements:
            result = result.replace(old, new)
        
        # Remove trailing underscores
        result = result.rstrip('_')
        
        return result


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
            print(seq_record)
            raise e
    return out


def write_genbank(seq_records, file_name, default_molecule_type="DNA", add_domain_cds_annotations=True, preserve_original=False, replace_space_with=""):
    """

        seq_records: an iterable of SeqRecords

        file_name can be a file handle

        This function will modify the input seqrecords as a side effect. To stop it from doing that, use "preserve_original=True"
    """
    warnings.filterwarnings("ignore", category=BiopythonWarning, module=".*InsdcIO.*")
    warnings.filterwarnings("ignore", category=BiopythonWarning, module=".*InsdcIO.*")


    outfile, file_type = open_if_is_name(file_name, "w")

    for record in seq_records:
        if preserve_original:
            record = copy_SeqRecord(record)

        if 'molecule_type' not in record.annotations:
            record.annotations['molecule_type'] = default_molecule_type
        if "date" not in record.annotations:
            record.annotations["date"] = date.today().strftime(
                "%d-%b-%Y").upper()

        if add_domain_cds_annotations:
            cds_domain_annotations = get_domain_cds_annotation_strings(record)
            
            for feature in record.features:
                if feature.type == "CDS":
                    seen_domain_dbs = set()
                    if "cds_id" in feature.qualifiers and feature.qualifiers["cds_id"][0] in cds_domain_annotations:
                        
                        for domain_db, annotation_string_list in cds_domain_annotations[feature.qualifiers["cds_id"][0]].items():
                            feature.qualifiers[f"domainator_{domain_db}"] = [
                                " & ".join(annotation_string_list)]
                            seen_domain_dbs.add(f"domainator_{domain_db}")
                    
                    for qual in list(feature.qualifiers.keys()): #delete qualifiers for databases where all of the domains have been filtered out.
                        if qual.startswith("domainator_") and qual not in seen_domain_dbs:
                            del feature.qualifiers[qual]

        swap_name_id(record)
        record.name = record.name.replace(" ",replace_space_with)
        record.features.sort(key=lambda x: x.location.start)
        SeqIO.write(record, outfile, "genbank")
        swap_name_id(record)

    if file_type == "name":
        outfile.close()

def swap_name_id(record):
    """
        necesarry because biopython handles genbank names and ids differently than other programs, such as Geneious
    """
    id = record.id
    record.id = record.name
    record.name = id

def copy_SeqRecord(seqrecord):
    new_record = SeqRecord(seq=seqrecord.seq, id=seqrecord.id, name=seqrecord.name,
                           description=seqrecord.description,
                           dbxrefs=seqrecord.dbxrefs.copy(),
                           features=seqrecord.features.copy(),
                           annotations=seqrecord.annotations.copy(),
                           letter_annotations=seqrecord.letter_annotations,
                           )
    return new_record


def parse_seqfiles(seqfiles, contigs=None, filetype_override=None, seek_to=None, max_recs=float("inf"), default_molecule_type='protein'):
    """
        args:
            seqfiles: a genbank or fasta path
            contigs: a dict where keys are the names of desired contigs, and values are None. Or None.
            filetype_override: if supplied then don't try to guess filetype from the extension, just use the supplied filetype
            seek_to: seek to this position in the file before reading any records
            max_recs: stop after returning this many sequences
        output:
            yields genbank records one at a time
    """
    

    for file in seqfiles:
        file_type = filetype_override
        (infile, input_type) = open_if_is_name(file)
        
        if input_type == 'handle' and file_type is None:
            raise ValueError(f"filetype_override must be specified when passing a handle.")

        if file_type is None:
            file_type = get_file_type(file)
            if file_type not in {"genbank", "fasta"}:
                raise ValueError(f"file extension not recognized: {file}, please specify format.")
        if seek_to is not None:
            infile.seek(seek_to)

        recs_returned = 0
        for rec in SeqIO.parse(infile, file_type):
            if recs_returned >= max_recs:
                break
            if file_type == "genbank":
                swap_name_id(rec)
            if contigs is None or rec.id in contigs:
                if 'molecule_type' not in rec.annotations:  # TODO: is this ok (will a nucleotide sequence always have molecule_type defined and a peptide sequence always not have?)
                    rec.annotations['molecule_type'] = default_molecule_type
                yield rec
                recs_returned += 1
                
        if input_type == "name":
            infile.close()

def split_string_list(names, delimeter="::"):
    """
        input:
            a list of strings, some of which might have internal delimeters.
        output:
            a list of lists of strings where input strings are split by the delimeter
    """
    return [[single_gff.strip() for single_gff in list_of_gffs.strip().split(delimeter)] for list_of_gffs in names]


def list_and_file_to_dict_keys(input_list=None, input_file=None, as_set=False):
    """
        args:
            input_list: a list of items to add, in order, as keys to the output dict.

            input_file: The path to a text file, either simple list of strings, a fasta file, a genbank file, or an hmm file. 
                        If a sequence file is passed, then return the sequence ids,
                        If an hmm file is passed then return the hmm profile names
                        If any other kind of file is passed, then just read the lines and pass the lines individually

        output:
            None if input_list and input_file are both None.
            else, a dict where keys are values read from input_list and input_file, and values are None.
    """

    out = None  # dict() we're treating this as an ordered set where values = None
    if input_list is not None and len(input_list) > 0:
        out = dict()
        for c in input_list:
            out[c] = None
    if input_file is not None:
        if out is None:
            out = dict()
        input_file_type = get_file_type(input_file)
        if input_file_type not in {"hmm", "fasta", "genbank"}: #if it's not an hmm file, a fasta file, or a genbank file, we assume it is a text file with one value per line
            with open(input_file, "r") as input_file_handle:
                for line in input_file_handle:
                    line = line.strip()
                    if len(line) != 0:
                        out[line] = None
        elif input_file_type == "hmm":
            hmms = read_hmms([input_file]) # dict of dicts of hmm objects, we're looking for the keys to the second level dicts, because those are the hmm names
            for hmmdb in hmms:
                for hmm_name in hmms[hmmdb]:
                    out[hmm_name] = None
        elif input_file_type in {"fasta","genbank"}:
            for rec in parse_seqfiles([input_file]):
                out[rec.id] = None

    if as_set and out is not None:
        return set(out.keys())

    return out


def regions_overlap(region1, region2, min_overlap_fraction=0.0):
    """
        regions are tuples of start and stop coordinates
        returns true if a fraction of region2 >= min_overlap_fraction overlaps with region1
        coordinates within regions must be sorted low to high
    """
    if min_overlap_fraction >= 1:
        return False

    # TODO: does biopython coordinate system require us to adjust by 1 somewhere?
    region_1_size = region1[1] - region1[0]
    region_2_size = region2[1] - region2[0]
    if region_1_size == 0:
        region_1_size = 0.1
    if region_2_size == 0:
        region_2_size = 0.1

    r1_s, r1_e = region1
    r2_s, r2_e = region2
    
    overlap_size = max(0, min(r1_e, r2_e) - max(r1_s, r2_s))
    if overlap_size == 0:
        return False
    else:
        return (overlap_size / region_2_size) >= min_overlap_fraction


def i_get_fasta_offsets(input_path): #TODO: replace with cython
    """
        input: a path to a fasta file

        yields: (file_offset, num_proteins)
                note that a fasta file can only have one protein per record, so the num_proteins will always be 1.
    
    """
    with open(input_path,"rb") as infile:
        offset=0
        line = infile.readline()
        while line:
            if line[0] == 62: # 62 = b'>'
                yield offset, 1
            offset = infile.tell()
            line = infile.readline()

def i_get_genbank_offsets(input_path): #TODO: replace with cython
    """
        input: a path to a genbank file

        yields: (file_offset, num_proteins)
    """

    last_offset = 0
    offset = 0
    cdss = 0
    rec_type = ""
    locus_seen = False
    with open(input_path,"rb") as infile:
        line = infile.readline()
        while line:
            if line.startswith(b"LOCUS "):
                parts=line.split()
                locus_seen = True
                if parts[3] == b"aa": #aa
                    rec_type = b"aa"
                    cdss = 1
                else: #bp
                    rec_type = b"bp"
                    cdss = 0

                last_offset=offset
                line = infile.readline()
                break
            offset = infile.tell()
            line = infile.readline()
        while line:
            if line.startswith(b"LOCUS "): # find subsequent "LOCUS" tags
                yield last_offset, cdss

                parts=line.split()
                if parts[3] == b"aa": #aa
                    rec_type = b"aa"
                    cdss = 1
                else: #bp
                    rec_type = b"bp"
                    cdss = 0

                last_offset=offset
            elif line.startswith(b"     CDS "):
                if rec_type == b"bp":
                    cdss += 1
            offset=infile.tell()
            line = infile.readline()
    if locus_seen:
        yield last_offset, cdss


def get_fasta_offsets(input_path): #TODO: replace with cython
    """
        input: a path to a fasta file

        output: two arrays of unsigned long longs, (file_offset, num_proteins)
                note that a fasta file can only have one protein per record, so the num_proteins will always be 1.
    
    """
    offsets = array("Q")
    num_proteins = array("Q")
    with open(input_path,"rb") as infile:
        offset=0
        line = infile.readline()
        while line:
            if line[0] == 62: # 62 = b'>'
                offsets.append(offset)
                num_proteins.append(1)
            offset = infile.tell()
            line = infile.readline()
    return offsets, num_proteins


def get_genbank_offsets(input_path): #TODO: replace with cython
    """
        input: a path to a genbank file

        output: two arrays of unsigned long longs, (file_offset, num_proteins)
    """

    last_offset = 0
    offset = 0
    cdss = 0
    offsets = array("Q")
    num_proteins = array("Q")
    rec_type = ""
    locus_seen = False
    with open(input_path,"rb") as infile:
        line = infile.readline()
        while line:
            if line.startswith(b"LOCUS "):
                parts=line.split()
                locus_seen = True
                if parts[3] == b"aa": #aa
                    rec_type = b"aa"
                    cdss = 1
                else: #bp
                    rec_type = b"bp"
                    cdss = 0

                last_offset=offset
                line = infile.readline()
                break
            offset = infile.tell()
            line = infile.readline()
        while line:
            if line.startswith(b"LOCUS "): # find subsequent "LOCUS" tags
                offsets.append(last_offset)
                num_proteins.append(cdss)
                parts=line.split()
                if parts[3] == b"aa": #aa
                    rec_type = b"aa"
                    cdss = 1
                else: #bp
                    rec_type = b"bp"
                    cdss = 0

                last_offset=offset
            elif line.startswith(b"     CDS "):
                if rec_type == b"bp":
                    cdss += 1
            offset=infile.tell()
            line = infile.readline()
    if locus_seen:
        offsets.append(last_offset)
        num_proteins.append(cdss)
    
    return offsets, num_proteins

def get_offsets(input_path):
    """
        input: a path to a fasta or genbank file

        output: two arrays of unsigned long longs, (offset, num_proteins)
    """

    filetype = get_file_type(input_path)
    if filetype == "fasta":
        offsets, num_proteins = get_fasta_offsets(input_path) # First array is file offsets of the record, second array is protein count (always 1 for fasta)

    elif filetype == "genbank":
        offsets, num_proteins = get_genbank_offsets(input_path) #
    else:
        raise ValueError(f"Filetype not recognized for input file: {input_path}")
    return offsets, num_proteins


def i_get_offsets(input_path):
    """
        input: a path to a fasta or genbank file

        returns an iterator of (offset, cdss)
    """

    filetype = get_file_type(input_path)
    if filetype == "fasta":
        return i_get_fasta_offsets(input_path) # First array is file offsets of the record, second array is protein count (always 1 for fasta)
    elif filetype == "genbank":
        return i_get_genbank_offsets(input_path) #
    else:
        raise ValueError(f"Filetype not recognized for input file: {input_path}")



def count_peptides(file_path): 
    _, num_proteins = get_offsets(file_path)
    return sum(num_proteins)
    
def count_peptides_in_record(rec):
    peptides = 0
    if rec.annotations['molecule_type'] == "protein":
        peptides += 1
    else:
        for feature in rec.features:
            if feature.type == 'CDS':
                peptides += 1
    return peptides

def slice_record(r:SeqRecord, start:int, stop:int, features:Iterable[SeqFeature]=None, truncate_features=True) -> SeqRecord:
    """ Adapted from Biopython SeqRecord::SeqRecord::__getitem__
        
        Unlike __getitem__, this function keeps 'source' features (resizing them to the size of the slice),
            as well as other metadata such as dbxrefs, contig 'source', 'organism', and 'taxonomy'.
        
        It also allows you to specify a subset of features to slice, instead of all features.

        Args:
            r (SeqRecord): The record to slice
            start (int): The start coordinate of the slice
            stop (int): The stop coordinate of the slice
            features (Iterable[SeqFeature]): The features to slice. If None, then all features will be sliced.
            truncate_features (bool): If True, then features that are partially within the slice will be truncated to fit within the slice.
                If False, then features that are partially within the slice will be omitted from the output record.
                source features will always be truncated, regardless of this setting.

        #TODO: maybe make this a method of SeqRecord, like make a slice method and have an option to keep source.
    """

    #TODO: avoid code duplication with slice_record_from_location, by making this function a wrapper for that one, as it is a special case.

    return slice_record_from_location(r, FeatureLocation(start, stop), features, truncate_features=truncate_features)

def pad_location(r:SeqRecord, location:Union[FeatureLocation,CompoundLocation], pad_upstream:int = 0, pad_downstream:int = 0, return_extension_sizes:bool = False) -> Union[Union[FeatureLocation,CompoundLocation],Tuple[Union[FeatureLocation,CompoundLocation],Tuple[int,int]]]:
    """ Pad a location by a specified amount of bases on either side.
        If the location is a CompoundLocation, then the first and last subparts will be padded.
        Padding will take strand information into account.
            if the first subpart is on the reverse strand, then upstream padding will be added to its end coordinate, 
            and if the last subpart is on the reverse strand, then downstream padding will be added to its start coordinate.
        
        Padding will also take into account the fact that the location might extend beyond the end of the sequence.
            If the padded location extends beyond the end of the sequence, then the padding will be reduced to the maximum possible padding.
        
        Finally, padding will take into account circular sequences.
            If padding extends beyond the origin, then padding will be extended beyond the origin, with new parts added to the location as necessary.
            If upstream padding overlaps with downstream padding, then the overlapping region will be divided equally between the two parts.

        Args:
            r (SeqRecord): The record to which the location belongs
            location (Union[FeatureLocation,CompoundLocation]): The location to pad
            pad_upstream (int): The number of bases/amino acids to pad upstream
            pad_downstream (int): The number of bases/amino acids to pad downstream
        
        Returns:
            Union[FeatureLocation,CompoundLocation]: The padded location
    """
    if pad_upstream < 0 or pad_downstream < 0:
        raise ValueError("Padding cannot be negative")
    if pad_upstream == 0 and pad_downstream == 0:
        if return_extension_sizes:
            return location, (0, 0)
        else:
            return location
    
    rec_len = len(r)
    if r.annotations.get("topology", None) != "circular":
        if len(location.parts) == 1:
            location = location.parts[0]
            if location.strand == -1:
                min_coord = max(0, int(location.start)-pad_downstream)
                downstream_extension = int(location.start)-min_coord
                max_coord = min(rec_len, int(location.end)+pad_upstream)
                upstream_extension = max_coord-int(location.end)

            else: # location.strand == 1 or location.strand == 0
                min_coord = max(0, int(location.start)-pad_upstream)
                upstream_extension = int(location.start)-min_coord
                max_coord = min(rec_len, int(location.end)+pad_downstream)
                downstream_extension = max_coord-int(location.end)
            out_location = FeatureLocation(
                    min_coord,
                    max_coord,
                    strand=location.strand
                )
            if return_extension_sizes:
                return out_location, (upstream_extension, downstream_extension)
            else:
                return out_location
        else: # CompoundLocation with multiple parts
            # in the weird case where the first part and last part are on different strands or where a middle part is beyond the bounds of the outer parts, this code could result in covering the same regions multiple times.
            first_part, upstream_extensions = pad_location(r, location.parts[0], pad_upstream, 0, return_extension_sizes=True)
            last_part, downstream_extensions = pad_location(r, location.parts[-1], 0, pad_downstream, return_extension_sizes=True)
            new_parts = [first_part]
            for part in location.parts[1:-1]:
                new_parts.append(part)
            new_parts.append(last_part)
            out_location = CompoundLocation(new_parts, operator=location.operator)
            if return_extension_sizes:
                return out_location, (upstream_extensions[0], downstream_extensions[1])
            else:
                return out_location
    else: # contig is circular
        if len(location.parts) == 1: # simple location (one segment)
            location = location.parts[0]
            l_start = int(location.start)
            l_end = int(location.end)
            #start_to_origin = int(location.start)
            end_to_origin = rec_len-int(location.end)
            
            # find the bounds of the region to be extracted
            if location.strand == -1:
                min_coord = location.start-pad_downstream
                max_coord = location.end+pad_upstream
            else: # location.strand == 1 or location.strand == 0
                min_coord = location.start-pad_upstream
                max_coord = location.end+pad_downstream
            
            if min_coord >= 0 and max_coord <= rec_len:
                # the bounds to not extend beyond the origin, on either side, so  we don't care if they overlap
                out_location = FeatureLocation(
                    min_coord,
                    max_coord,
                    strand=location.strand
                )
                if return_extension_sizes:
                    return out_location, (pad_upstream, pad_downstream)
                else:
                    return out_location
            else:
                # If regions extending beyond the origin run into the FeatureLocation, truncate them there.
                # We prefer to not create a copy of any part of the contig.
                
                # adjust min_coord and max_coord so that they don't wrap around the FeatureLocation bounds.
                if min_coord < 0 and ((-1 * min_coord) > end_to_origin):
                    min_coord = -1 * end_to_origin #TODO: test for off-by-one
                if max_coord > rec_len and (max_coord - rec_len > l_start):
                    max_coord = rec_len + l_start
                
                # If the adjusted regions overlap with each other, then divide the overlap equally between min and max.
                # shift the perspective so that min_coord is at the origin
                min_coord_shifted = 0 # min_coord_shifted = min_coord - min_coord
                max_coord_shifted = max_coord - min_coord
                
                # calculuate overlap size:
                overlap_size = max_coord_shifted - rec_len
                if overlap_size > 0:
                    # divide overlap equally between min and max
                    min_coord_shifted += overlap_size // 2
                    max_coord_shifted -= overlap_size // 2
                    # if overlap_size is odd, then we need to adjust one of the coordinates by 1
                    if overlap_size % 2 == 1:
                         max_coord_shifted -= 1 
                # shift the perspective back to the original coordinate system
                max_coord = max_coord_shifted + min_coord
                min_coord = min_coord_shifted + min_coord

                if location.strand == -1:
                    upstream_extension = max_coord-l_end
                    downstream_extension = l_start-min_coord 
                else: # location.strand == 1 or location.strand == 0
                    upstream_extension = l_start-min_coord
                    downstream_extension = max_coord-l_end 


                parts = list() # we're touching or wrapping around the origin, so we may have multiple parts to extract
                if min_coord < 0: # min_coord is across the origin from the start of the FeatureLocation, so we need to break it up into two parts
                    parts.append(FeatureLocation(rec_len+min_coord, rec_len, strand=location.strand)) #from the start of the extraction window to the origin
                    if l_end > 0: # not touching the origin, so we need to add a part for the region from the origin to the end of the FeatureLocation
                        parts.append(FeatureLocation(0, l_end, strand=location.strand))
                elif min_coord != l_start: # min_coord is not across the origin from the start of the FeatureLocation
                    parts.append(FeatureLocation(min_coord, l_end, strand=location.strand))
                else: # min_coord == l_start
                    parts.append(FeatureLocation(l_start, l_end, strand=location.strand))

                if max_coord > rec_len: # max_coord is across the origin from the end of the FeatureLocation, so we need to break it up into two parts
                    parts[-1]._end = ExactPosition(rec_len)
                    parts.append(FeatureLocation(0, max_coord-rec_len, strand=location.strand))
                else:
                    parts[-1]._end = ExactPosition(max_coord)

                if location.strand == -1:
                    parts.reverse()
                
                if len(parts) == 1:
                    out_location = parts[0]
                else:
                    out_location = CompoundLocation(parts, operator="join")

                if return_extension_sizes:
                    return out_location, (upstream_extension, downstream_extension)
                else:
                    return out_location
        else: # CompoundLocation with more than one part
            # we treat the CompoundLocation as a single FeatureLocation, so we need to calculate the bounds of the CompoundLocation
            # we consider the Strand of the CompoundLocation to be the Strand of the first part
            # TODO: handle the case where the first part is on the opposite strand of the last part.
            # TODO: handle the case where a middle part is beyond the bounds of the outer parts
            
            part_1_strand = location.parts[0].strand
            envelope = location.circular_envelope(len(r), part_1_strand)
            if len(envelope) == len(r): # envelope covers the entire contig, so just return the starting location
                return location
            elif len(envelope.parts) == 1: # envelope is a single FeatureLocation (i.e. it does not cross the origin, so the padding might cross the origin)
                _, extensions = pad_location(r, envelope.parts[0], pad_upstream, pad_downstream, return_extension_sizes=True)
                
                out_parts = list()
                start_parts = pad_location(r, location.parts[0], extensions[0], 0)
                end_parts = pad_location(r, FeatureLocation(location.parts[-1].start, location.parts[-1].end, part_1_strand), 0, extensions[1])
                out_parts.extend(start_parts.parts)
                out_parts.extend(location.parts[1:-1])
                out_parts.extend(end_parts.parts)
                out_location = CompoundLocation(out_parts, operator="join")
                if return_extension_sizes:
                    return out_location, extensions
                else:
                    return out_location
            else: # envelope crosses the origin, implying that the padding will not cross the origin
                # adjust pad_upstream and pad_downstream so that they won't run into the envelope.
                l_start = envelope.stranded_start
                l_end = envelope.stranded_end 
                remaining_span = abs(l_start - l_end) # size of the record minus the size of the envelope
                pad_upstream = min(pad_upstream, remaining_span) # adjust pad_upstream so that it won't run into the envelope
                pad_downstream = min(pad_downstream, remaining_span) # adjust pad_downstream so that it won't run into the envelope
                padding_size = pad_upstream + pad_downstream
                if padding_size > remaining_span:
                    padding_size_diff = padding_size - remaining_span # 
                    # divide the padding_size_diff equally between pad_upstream and pad_downstream
                    pad_upstream -= padding_size_diff // 2
                    pad_downstream -= padding_size_diff // 2
                    # if padding_size_diff is odd, then we need to adjust one of the coordinates by 1
                    if padding_size_diff % 2 == 1:
                        pad_downstream -= 1
                out_parts = list()
                first_part, upstream_extensions = pad_location(r, location.parts[0], pad_upstream, 0, return_extension_sizes=True)
                last_part, downstream_extensions = pad_location(r, FeatureLocation(location.parts[-1].start, location.parts[-1].end, part_1_strand), 0, pad_downstream, return_extension_sizes=True)
                out_parts.extend(first_part.parts) # first part should only have one part
                out_parts.extend(location.parts[1:-1])
                out_parts.extend(last_part.parts) # last part should only have one part
                out_location = CompoundLocation(out_parts, operator="join")
                if return_extension_sizes:
                    return out_location, (upstream_extensions[0], downstream_extensions[1])
                else:
                    return out_location

                
def slice_record_from_location(r:SeqRecord, location:Union[FeatureLocation,CompoundLocation], features:Iterable[SeqFeature]=None, pad_upstream:int = 0, pad_downstream:int = 0, truncate_features=False) -> SeqRecord:
    """Extract a sub-record from supplied parent record using the location of the SeqFeature object.
        Output record will have the same id, name, description, dbxrefs, and annotations as the parent record.
        Output record will have the same features as the parent record, except that features that are not completely contained in the location will be removed.
        Output record will have the same letter annotations as the parent record, except that the letter annotations will be sliced to match the sliced sequence.

        sub-parts of the location will be extracted from the parent record, and then concatenated together to form the output record.
        If the location is a CompoundLocation, then the sub-parts will be extracted in the order they are listed in the CompoundLocation.
        parts on the reverse strand will be extracted in reverse order, and then reversed again to be in the correct orientation.


        

        Adapted from Biopython SeqRecord::SeqRecord::__getitem__
        
        Unlike __getitem__, this function keeps 'source' features (resizing them to the size of the slice),
            as well as other metadata such as dbxrefs, contig 'source', 'organism', and 'taxonomy'.
        
        It also allows you to specify a subset of features to slice, instead of all features.

        r (SeqRecord): parent record
        location (FeatureLocation): location to extract
        features (Iterable[SeqFeature]): subset of features to slice. If None, all features will be sliced.
        pad_upstream (int): number of bases/residues to pad upstream of the location
        pad_downstream (int): number of bases/residues to pad downstream of the location
        truncate_features (bool): if True, then features that are not completely contained in the location will 
            be removed from the output record. If False, then features that are not completely contained in the location will be 
            included in the output record, but will be truncated to the size of the location. "source" features will always be included, regardless of this setting.

        Returns:
            SeqRecord: sliced record

        TODO: should padding be handled here, or should we expect the FeatureLocation to already be padded?
    """

    # TODO: to deal with padding 
    #       linear contigs:   add padding to the first or last part of the location, depending on the strand, and slice the record from there.
    #       circular contigs: add padding to the first and last part of the location, depending on the strand, truncate the padding if there is overlap, 
    #                           and then slice the record from there. 
    #                           No rotation should be necessary!, for example if you are slicing a kb-range across the origin, you just need to define a compound 
    #                           location with one part on either side of the origin, and then add padding to both parts, and then slice the record from there. No rotation needed.)

    if r.seq is None:
        raise ValueError("If the sequence is None, we cannot slice it.")

    if features is None:
        features = r.features

    if pad_upstream < 0 or pad_downstream < 0:
        raise ValueError("Padding must be a positive integer.")
    if pad_upstream != 0 or pad_downstream != 0:
        location = pad_location(r, location, pad_upstream, pad_downstream)
    # get sequence parts
    seq_parts = list()
    seq_parts_len = 0
    new_features = list()
    coord_shifts = list()
    strands = list()
    # location_parts = merge_location_parts(location.parts) #TODO: something like this might fix edge cases where an annotation is completely enclosed by the location, but for some reason the location has a size-zero split with the overlap region.
    for part in location.parts:
        strand = part.strand
        strands.append(strand)
        seq_parts.append(part.extract(r.seq))

        # If the part is on the reverse strand, then to calculate shift for overlaping parts,
        # we need to shift left by the length of the part, flip, then add to the length of the previous parts.
        # 
        # if the part is on the forward strand, then to calculate shift for overlaping parts,
        # we need to shift left by the length of the part, then add to the length of the previous parts.
        # 
        # in either case, we need to add the length of the previous parts to the shift.

        coord_shifts.append(seq_parts_len)         
        
        seq_parts_len += len(seq_parts[-1])
    
    # include features if they completely overlap with the location, or they are a source feature
    for f in features:
        f_parts = f.location.parts
        overlap = [set() for _ in range(len(f_parts))] # list of sets of indices of location parts that overlap with the feature part
        total_overlap = 0
        for feature_part_index, f_part in enumerate(f_parts):
            part_found = False
            f_start = int(f_part.start)
            f_end = int(f_part.end)
            for loc_part_index, loc_part in enumerate(location.parts):
                loc_start = int(loc_part.start)
                loc_end = int(loc_part.end)
                if ( # part is completely contained in loc_part
                    loc_start <= f_start
                    and loc_end >= f_end
                ):
                    overlap[feature_part_index].add(loc_part_index)
                    part_found = True
                elif f.type == "source" or truncate_features: # part is partially contained in loc_part
                    if regions_overlap((f_start, f_end), (loc_start, loc_end)): #TODO: should only be overlap for source features, all other features should be completely contained.
                        overlap[feature_part_index].add(loc_part_index)
                        part_found = True
            if part_found:
                total_overlap += 1
        if total_overlap == len(f_parts) or f.type == "source" or truncate_features:
            new_parts = list()
            for feature_part_index, f_part in enumerate(f_parts):
                f_start = int(f_part.start)
                f_end = int(f_part.end)
                for loc_part_index in overlap[feature_part_index]:
                    loc_part = location.parts[loc_part_index]
                    loc_start = int(loc_part.start)
                    loc_end = int(loc_part.end)
                    loc_part_strand = strands[loc_part_index]
                    
                    new_start = max(f_start - loc_start, 0) # if the feature does not completely overlap with the extracted region, only include the part that does
                    new_end = min(f_end-loc_start, loc_end-loc_start)
                    new_loc = FeatureLocation(new_start, new_end, strand=f_part.strand)
                    if loc_part_strand == -1:
                        new_loc = new_loc._flip(loc_end-loc_start)
                    new_loc = new_loc._shift(coord_shifts[loc_part_index])
                    new_parts.append(new_loc)
            # Merge parts that are adjacent or overlap
            if getattr(f.location, "operator", None) == "join" or getattr(f.location, "operator", None) is None: # only merge parts if the original location was a join or not a CompoundLocation (for example, an origin-spanning source) TODO: should we merge parts for any other operators?
                new_parts = CompoundLocation.merge_parts(new_parts) # TODO: There are some edge cases where this will merge things that maybe users won't want to merge, but it's better than not merging at all. Can fix later if anyone complains.

            if len(new_parts) == 0: # this can happen for source features that don't intersect with the location
                pass
            else:
                if len(new_parts) == 1:
                    new_location = new_parts[0]
                elif len(new_parts) > 1:
                    new_operator = getattr(f.location, "operator", "join")
                    new_location = CompoundLocation(new_parts, new_operator) 
                new_feature = f.__class__(
                    new_location,
                    type=f.type,
                    id=f.id,
                    qualifiers=f.qualifiers,
                )
                new_features.append(new_feature)
    # TODO: handle other partially cut annotations



    answer = r.__class__(
        functools.reduce(lambda x, y: x + y, seq_parts),
        id=r.id,
        name=r.name,
        description=r.description,
        dbxrefs=r.dbxrefs.copy(),
        features=new_features,
    )

    if "molecule_type" in r.annotations:
        # This will still apply, and we need it for GenBank/EMBL etc output
        answer.annotations["molecule_type"] = r.annotations["molecule_type"]

    if "source" in r.annotations:
        answer.annotations["source"] = r.annotations["source"]

    if "organism" in r.annotations:
        answer.annotations["organism"] = r.annotations["organism"]

    if "taxonomy" in r.annotations:
        answer.annotations["taxonomy"] = r.annotations["taxonomy"]

    # # if "references" in r.annotations: #TODO: maybe keep this also?
    # #     answer.annotations["references"] = r.annotations["references"]

    # TODO: handle letter_annotations (genbank files don't have these, so probably not relevant for Domainator)

    return answer

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

def get_taxid(record:SeqRecord) -> Optional[int]:
    """
        returns 32644, "unidentified" if no taxid is found
    """
    taxid = 32644
    if "ncbi_taxid" in record.annotations: # from SwissProt-style text file
        return int(record.annotations["ncbi_taxid"][0])
    elif " OX=" in record.description: # from SwissProt-style fasta file
        return int(record.description.split(" OX=")[-1].split(" ")[0])
    elif " TaxID=" in record.description: # from Uniref-style fasta file
        return int(record.description.split(" TaxID=")[-1].split(" ")[0])
    else: # from GenBank-style genbank file, take taxid from the longest source feature
        sources = get_sources(record)
        sources.sort(key=len, reverse=True)
        if len(sources) > 0:
            if "db_xref" in sources[0].qualifiers:
                for xref in sources[0].qualifiers["db_xref"]:
                    if xref.startswith("taxon:"):
                        try:
                            return int(xref.split(":")[-1])
                        except ValueError:
                            pass
    return taxid

def filter_by_taxonomy(records, include_taxids, exclude_taxids, ncbi_taxonomy):
    for record in records:
        taxid = get_taxid(record)
        if taxid is None: # no taxid, skip
            continue

        lineage = set(ncbi_taxonomy.lineage(taxid))
        if include_taxids and not lineage.intersection(include_taxids): # not in include_taxids, skip
            continue
        if exclude_taxids and lineage.intersection(exclude_taxids): # in exclude_taxids, skip
            continue
        yield record


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
        location = feature.location.copy()
    return SeqFeature(
                location=location,
                type=feature.type,
                id=feature.id,
                qualifiers=feature.qualifiers.copy(),
            )

def circular_dist(start, end, size):
    if end < start:
        end = end + size
    return end - start


def get_palette(values):
    """
        values is a list of unique values

        returns a dict mapping values to hex colors
    """
    out = dict()
    #colors = sns.color_palette('husl', n_colors=len(values)).as_hex()
    colors = sns.color_palette('tab20', n_colors=len(values)).as_hex() + sns.color_palette('tab20b', n_colors=len(values)).as_hex() + sns.color_palette('tab20c', n_colors=len(values)).as_hex()
    if len(values) < len(colors):
        colors = colors[:len(values)]
    # random.Random(12).shuffle(colors)

    for i in range(len(values)):
        color_i = i % len(colors)
        if pd.isna(values[i]):
            out[None] = colors[color_i]
        else:
            out[values[i]] = colors[color_i]
    return out