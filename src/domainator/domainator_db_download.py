"""
    Download sequence databases from the internet in formats appropriate for Domainator
    
    Including the following databases:
    - Uniprot (includes Trembl and SwissProt) (downloads fasta by default, use _gb suffix to download genbank file with addtional annotations)
    - Trembl (downloads fasta by default, use _gb suffix to download genbank file with addtional annotations)
    - SwissProt (downloads fasta by default, use _gb suffix to download genbank file with addtional annotations)
    - ncbi_complete_genome_proks: NCBI genomes where assembly_level is "Complete Genome". (excluding eukaryotes)
    - ncbi_complete_genome_nonredundant_proks: non-redundant NCBI genomes where assembly_level is "Complete Genome". (excluding eukaryotes)
    - ncbi_representative_proks : representative and reference genomes from the NCBI database (excluding eukaryotes)
    - ncbi_representative_all : representative and reference genomes from the NCBI database (all taxa, including eukaryotes)
    - ncbi_nonredundant_proks : non-redundant genomes from the NCBI database (excluding eukaryotes)
    - ncbi_nonredundant_all : non-redundant genomes from the NCBI database (all taxa, including eukaryotes)
    - ncbi_all : all genomes from the NCBI database including duplicates of the same species (all taxa, including eukaryotes). Can produce very large files, up to 10 TB or more, so use with caution.
 

    Key taxa:
        81490 : prokaryotic environmental samples
        408169 : metagenomes
        2 : bacteria
        2157 : archaea
        10239 : viruses
        
        2759 : eukaryota

    WARNING: Downloading large databases can take a long time and use a lot of disk space!
"""

from jsonargparse import ArgumentParser, ActionConfigFile
import gzip
import requests
from domainator.Bio import SeqIO
from domainator.Taxonomy import NCBITaxonomy
from pathlib import Path
import sys
from domainator.Bio.SeqFeature import FeatureLocation, SeqFeature, UnknownPosition
from domainator.utils import write_genbank, parse_seqfiles, filter_by_taxonomy, list_and_file_to_dict_keys
from domainator import __version__, RawAndDefaultsFormatter
from typing import Dict, List, Set, Tuple, Union, Optional, Iterator
from domainator.domainate import clean_rec, prodigal_CDS_annotate
import tqdm
import datetime
import tempfile
from multiprocessing import Pool, Manager
import psutil
import functools
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry, MaxRetryError

# TODO: write tests for genbank downloads.
# TODO: NCBI TSA, mgnify, other databases, maybe with automatic taxonomy assignment if not annotated, that would be cool!

##### NCBI #####

def load_ncbi_assembly_summary(file_path: Optional[str] = None) -> Iterator[Dict[str, str]]:
    """
    Downloads and parses a GenBank assembly summary file (e.g., the file found at https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt).

    Args:
        file_path (str, optional): Path to a local copy of a GenBank assembly summary file. If None, the file will be downloaded from the URL.

    Returns:
        Iterator[Dict[str, str]]: An iterator of dictionaries where keys are column labels and values are row entries.
    """

    def parse_file(file):
        # Discard the first line (it's a comment)
        next(file)

        # Parse the header
        header = next(file).strip().split('\t')
        header[0] = header[0].strip('# ')

        # Parse the rest of the file
        for line in file:
            line = line.strip()
            if len(line) > 0:
                parts = line.split("\t")
                out = dict(zip(header, parts))
                out["versionless_accession"] = out["assembly_accession"].split(".")[0]
                yield out

    print(f"Loading assembly_summary_genbank.txt", file=sys.stderr)
    # Open the file and parse it
    if file_path is not None:
        with open(file_path) as f:
            yield from parse_file(f)
    else:
        url = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
        with requests.get(url, stream=True) as response:
            response.raise_for_status()
            f = (line.decode('utf-8') for line in response.iter_lines()) # returns a generator
            yield from parse_file(f)

def download_ncbi(output_handle, include_taxa=None, exclude_taxa=None, taxonomy_database=None, filter_by_uniqueness=True, filter_by_complete_genome=False, no_na=False, gene_call=None, num_recs=None, ncbi_summary=None, before:Optional[datetime.date]=None, after:Optional[datetime.date]=None, cpus:int=1, success_rec_log=None, exclude_accessions:Optional[Set[str]]=None):
    """
    Downloads nucleotide data from NCBI GenBank and writes it to a genbank file.

    Args:
        output_handle (file-like object): The file handle to which the downloaded GenBank data will be written.
        include_taxa (list, optional): A list of NCBI taxids of clades to download. If supplied, only assemblies with at least one taxid from this list in their lineage will be downloaded.
        exclude_taxa (list, optional): A list of NCBI taxids of clades to exclude from downloading. Assemblies with any taxid from this list in their lineage will not be downloaded.
        taxonomy_database (object, optional): A domainator.Taxonomy.NCBITaxonomy object. If None, both include_taxids and exclude_taxids must be None.
        filter_by_uniqueness (bool, optional): If True, only one genome per taxid will be downloaded, except for metagenomes and prokaryotic environmental samples, where all genomes will be downloaded. Defaults to True.
        filter_by_complete_genome: If True, only genomes with assembly_level == "Complete Genome" will be downloaded. Defaults to False.
        no_na (bool, optional): If True, only genomes without 'na' in the 'refseq_category' column of the assembly_summary_genbank.txt file will be downloaded. Defaults to False.
        gene_call (str, optional): If 'all', then any existing CDS annotations will be deleted and all contigs will be re-annotated. If 'unannotated', then only contigs without CDS annotations will be annotated. If None, existing annotations will be kept. Defaults to None.
        num_recs (int, optional): The maximum number of records to download. If None, all records will be downloaded. Defaults to None.
        ncbi_summary (file-like object, optional): The file handle of the NCBI assembly summary file. If None, the assembly summary file will be downloaded. Defaults to None.
        before (datetime.date, optional): If specified, only assemblies with a 'last_updated' date before this date will be downloaded. Defaults to None.
        after (datetime.date, optional): If specified, only assemblies with a 'last_updated' date after this date will be downloaded. Defaults to None.
        cpus (int, optional): The number of CPUs to use for downloading. Defaults to 1.
        success_rec_log (file-like object, optional): The file handle to which the accessions of successfully downloaded records will be written. Defaults to None.
        exclude_accessions (set, optional): A set of accessions to exclude from downloading. Defaults to None.


    Raises:
        ValueError: If the taxonomy database is None and either include_taxa or exclude_taxa are specified or filter_by_uniqueness is True.
        ValueError: If the gene_call parameter is not one of None, 'all', or 'unannotated'.
    """

    if (taxonomy_database is None and (include_taxa is not None or exclude_taxa is not None or filter_by_uniqueness == True)):
        raise ValueError(f"if include_taxa, exclude_taxa, or filter_by_uniqueness are specified, then taxonomy_database must also be supplied")

    if gene_call not in (None, 'all', 'unannotated'):
        raise ValueError("gene_call must be one of None, 'all', or 'unannotated'")

    if include_taxa is not None:
        include_taxa = set(include_taxa)

    if exclude_taxa is not None:
        exclude_taxa = set(exclude_taxa)
    if exclude_accessions is not None:
        exclude_accessions = set(exclude_accessions)

    # Download the index file
    genbank_accessions = list(load_ncbi_assembly_summary(ncbi_summary))

    genbank_accessions = filter_by_ftp_path(genbank_accessions)
    genbank_accessions = filter_by_accession_exclusion(genbank_accessions, exclude_accessions)
    genbank_accessions = filter_by_date(genbank_accessions, before, after)
    genbank_accessions = filter_genbank_accessions_by_taxonomy(genbank_accessions, include_taxa, exclude_taxa, taxonomy_database)
    genbank_accessions = filter_by_unique_taxid(genbank_accessions, taxonomy_database) if filter_by_uniqueness else genbank_accessions
    genbank_accessions = filter_by_assembly_level(genbank_accessions, categories={"Complete Genome"}) if filter_by_complete_genome else genbank_accessions
    genbank_accessions = filter_by_refseq_category(genbank_accessions) if no_na else genbank_accessions

    process_genbank_accessions(genbank_accessions, output_handle, gene_call, num_recs, cpus, success_rec_log)

def filter_by_accession_exclusion(genbank_accessions, exclude_accessions):
    """
    Filters a list of GenBank accessions based on a set of excluded accessions.

    Args:
        genbank_accessions (list): A list of dictionaries containing GenBank accession information.
        exclude_accessions (set): A set of accessions to exclude from downloading.

    Returns:
        list: A filtered list of GenBank accessions.
    """

    if exclude_accessions is not None:
        genbank_accessions = [r for r in genbank_accessions if r['versionless_accession'] not in exclude_accessions]
        print(f"{len(genbank_accessions)} accessions after filtering by accession exclusion", file=sys.stderr)

    return genbank_accessions

def filter_by_date(genbank_accessions, before:Optional[datetime.date], after:Optional[datetime.date]):
    """
    Filters a list of GenBank accessions based on their 'last_updated' date.

    Args:
        genbank_accessions (list): A list of dictionaries containing GenBank accession information.
        before (datetime.date, optional): If specified, only assemblies with a 'last_updated' date before this date will be downloaded.
        after (datetime.date, optional): If specified, only assemblies with a 'last_updated' date after this date will be downloaded.

    Returns:
        list: A filtered list of GenBank accessions.
    """

    if before is not None:
        genbank_accessions = [r for r in genbank_accessions if datetime.datetime.strptime(r['seq_rel_date'], "%Y/%m/%d").date() < before]
        print(f"{len(genbank_accessions)} accessions before {before}", file=sys.stderr)

    if after is not None:
        genbank_accessions = [r for r in genbank_accessions if datetime.datetime.strptime(r['seq_rel_date'], "%Y/%m/%d").date() > after]
        print(f"{len(genbank_accessions)} accessions after {after}", file=sys.stderr)

    return genbank_accessions

def filter_by_ftp_path(genbank_accessions):
    """
    Filters a list of GenBank accessions based on the presence of 'na' in the FTP path.

    Args:
        genbank_accessions (list): A list of dictionaries containing GenBank accession information.

    Returns:
        list: A filtered list of GenBank accessions without 'na' in their FTP path.
    """

    filtered_gb_accessions = [r for r in genbank_accessions if r['ftp_path'].strip() != "na"]
    print(f"{len(filtered_gb_accessions)} accessions in assembly_summary_genbank.txt", file=sys.stderr)
    return filtered_gb_accessions

def filter_genbank_accessions_by_taxonomy(genbank_accessions, include_taxa, exclude_taxa, taxonomy_database):
    """
    Filters a list of GenBank accessions based on their taxonomy.

    Args:
        genbank_accessions (list): A list of dictionaries containing GenBank accession information.
        include_taxa (list): A list of NCBI taxids of clades to include in the filtered list.
        exclude_taxa (list): A list of NCBI taxids of clades to exclude from the filtered list.
        taxonomy_database (object): A domainator.Taxonomy.NCBITaxonomy object.

    Returns:
        list: A filtered list of GenBank accessions based on the provided taxonomy criteria.
    """
    if include_taxa is None and exclude_taxa is None:
        return genbank_accessions

    filtered_gb_accessions = []
    for r in genbank_accessions:
        taxid = int(r['taxid'])
        lineage = set(taxonomy_database.lineage(taxid))
        should_download = True

        if include_taxa is not None:
            if len(lineage.intersection(include_taxa)) == 0:
                should_download = False
        if exclude_taxa is not None:
            if len(lineage.intersection(exclude_taxa)) > 0:
                should_download = False

        if should_download:
            filtered_gb_accessions.append(r)

    print(f"{len(filtered_gb_accessions)} accessions in assembly_summary_genbank.txt after filtering by taxonomy", file=sys.stderr)
    return filtered_gb_accessions

def filter_by_unique_taxid(genbank_accessions, taxonomy_database):
    """
    Filters a list of GenBank accessions to include only one genome per taxid, except for metagenomes and prokaryotic environmental samples.

    Args:
        genbank_accessions (list): A list of dictionaries containing GenBank accession information.
        taxonomy_database (object): A domainator.Taxonomy.NCBITaxonomy object.

    Returns:
        list: A filtered list of GenBank accessions containing unique taxids, with an exception for metagenomes and prokaryotic environmental samples.
    """

    METAGENOME_TAXID = 81490
    PROKARYOTIC_ENV_SAMPLE_TAXID = 408169
    
    filtered_gb_accessions = []
    genbank_accessions.sort(reverse=True, key=lambda x: x['refseq_category'])
    seen = set()

    for r in genbank_accessions:
        taxid = int(r['taxid'])
        lineage = set(taxonomy_database.lineage(taxid))

        if (taxid not in seen) or (METAGENOME_TAXID in lineage) or (PROKARYOTIC_ENV_SAMPLE_TAXID in lineage):
            seen.add(taxid)
            filtered_gb_accessions.append(r)

    print(f"{len(filtered_gb_accessions)} accessions in assembly_summary_genbank.txt after filtering by unique taxid", file=sys.stderr)
    return filtered_gb_accessions

def filter_by_refseq_category(genbank_accessions):
    """
    Filters a list of GenBank accessions based on the 'refseq_category' field.

    Args:
        genbank_accessions (list): A list of dictionaries containing GenBank accession information.

    Returns:
        list: A filtered list of GenBank accessions with 'refseq_category' field not equal to 'na'.
    """

    filtered_gb_accessions = [r for r in genbank_accessions if r['refseq_category'] != "na"]
    print(f"{len(filtered_gb_accessions)} accessions in assembly_summary_genbank.txt after filtering by refseq_category != na", file=sys.stderr)
    return filtered_gb_accessions

def filter_by_assembly_level(genbank_accessions, categories={"Complete Genome"}):
    """
    Filters a list of GenBank accessions based on the 'assembly_level' field.

    Args:
        genbank_accessions (list): A list of dictionaries containing GenBank accession information.
        categories (set): A set of assembly levels to include in the filtered list.

    Returns:
        list: A filtered list of GenBank accessions.
    """

    filtered_gb_accessions = [r for r in genbank_accessions if r['assembly_level'] in categories]
    print(f"{len(filtered_gb_accessions)} accessions in assembly_summary_genbank.txt after filtering by assembly_level in {categories}", file=sys.stderr)
    return filtered_gb_accessions

def worker_process(r, output_file_path, gene_call, num_recs, records_written, lock):

    if gene_call not in (None, 'all', 'unannotated'):
        raise ValueError("gene_call must be one of None, 'all', or 'unannotated'")

    clear_domainator_annotations = True
    clear_best_annotation = True
    clear_CDS_annotations = False
    if gene_call == 'all':
        clear_best_annotation = True
        clear_domainator_annotations = True
        clear_CDS_annotations = True

    assembly_long_name = r['ftp_path'].split('/')[-1]
    remote_nucleotide_file = r['ftp_path'] + '/' + assembly_long_name + "_genomic.gbff.gz"
    
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=0.1, status_forcelist=[ 500, 502, 503, 504 ])
    adapter = HTTPAdapter(max_retries=retries)
    session.mount('https://', adapter)
    
    try:
        with session.get(remote_nucleotide_file, stream=True) as response:
            response.raise_for_status()
            
            if gene_call is None:
                with gzip.open(response.raw, "rt") as gzipped_file: #TODO: if more than one thread, maybe want to download and then unzip instead of streaming, for better parallelization?
                    with lock:
                        if num_recs is not None and records_written.value >= num_recs:
                            return False
                        with open(output_file_path, "a") as output_handle:
                            for line in gzipped_file:
                                output_handle.write(line)
                        if num_recs is not None:
                            records_written.value += 1
            else:
                with gzip.open(response.raw, "rt") as gzipped_file:
                    with tempfile.TemporaryDirectory() as output_dir:
                        temp_file = output_dir + "/temp.gb"
                        with open(temp_file, "w") as temp_handle:
                            for line in gzipped_file:
                                temp_handle.write(line)
                        records = list(parse_seqfiles([temp_file], filetype_override="genbank"))
                        for record in records:
                            CDS_count = clean_rec(record, clear_best_annotation, clear_domainator_annotations, clear_CDS_annotations)
                            if CDS_count == 0 and gene_call is not None:
                                prodigal_CDS_annotate(record)
                        with lock:
                            with open(output_file_path, "a") as output_handle:  
                                if num_recs is not None and records_written.value >= num_recs:
                                    return False
                                
                                for record in records:                         
                                    write_genbank((record,), output_handle)
                                
                                if num_recs is not None:
                                    records_written.value += 1
        return r["versionless_accession"]

    except requests.RequestException as e:
        print(f"error downloading {remote_nucleotide_file}: {e}", file=sys.stderr)
        return False
    except EOFError as e:
        print(f"error parsing {remote_nucleotide_file}: {e}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"error processing {remote_nucleotide_file}: {e}", file=sys.stderr)
        return False

def process_genbank_accessions(genbank_accessions, output_file_path, gene_call, num_recs, cpus:int=1, success_rec_log=None):
    """
    Processes a list of GenBank accessions, downloading the corresponding nucleotide data and writing it to the provided output handle.

    Args:
        genbank_accessions (list): A list of dictionaries containing GenBank accession information.
        output_handle (file-like object): The file handle to which the downloaded GenBank data will be written.
        gene_call (str): If 'all', then any existing CDS annotations will be deleted and all contigs will be re-annotated. If 'unannotated', then only contigs without CDS annotations will be annotated. If None, existing annotations will be kept.
        num_recs (int): The maximum number of records to download. If None, all records will be downloaded.
        skipped_record_log (file-like object, optional): The file handle to which information about skipped records will be written. Defaults to None.
        cpus (int, optional): The number of CPUs to use for parallelization. Defaults to 1.
    """

    if success_rec_log is not None:
        success_log_handle = open(success_rec_log, "w")
        

    with Manager() as manager:
        lock = manager.Lock()
        records_written = manager.Value('i', 0)
        with Pool(processes=cpus) as pool:
            for success in tqdm.tqdm(pool.imap_unordered(functools.partial(worker_process, output_file_path=output_file_path, gene_call=gene_call, num_recs=num_recs, records_written=records_written, lock=lock), genbank_accessions), total=len(genbank_accessions), desc="GenBank genomes", leave=True, dynamic_ncols=True):
                with lock:
                    if success is not False and success_rec_log is not None:
                        success_log_handle.write(success + "\n")
                    if num_recs is not None and records_written.value >= num_recs:
                        break
    
    if success_rec_log is not None:
        success_log_handle.close()



##### Uniprot #####

def download_and_convert_uniprot_dat(url, output_file_path, include_taxids, exclude_taxids, ncbi_taxonomy, num_recs=None):
    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        with open(output_file_path, "w") as output_handle:
            with gzip.open(response.raw, "rt") as gzipped_file:
                recs_written = 0
                records = tqdm.tqdm(SeqIO.parse(gzipped_file, "swiss"), desc=url.split("/")[-1], leave=True, dynamic_ncols=True)
                if include_taxids or exclude_taxids:
                    records = filter_by_taxonomy(records, include_taxids, exclude_taxids, ncbi_taxonomy)
                
                for record in records:
                    features = [SeqFeature(FeatureLocation(0, len(record)), type="source", qualifiers={"db_xref": f"taxon:{record.annotations['ncbi_taxid'][0]}"})]
                    for feature in record.features:
                        skip = False
                        for part in feature.location.parts:
                            if isinstance(part.start, UnknownPosition) or isinstance(part.end, UnknownPosition):
                                skip = True
                                break
                        if not skip:
                            features.append(feature)
                    record.features = features
                    write_genbank((record,), output_handle)
                    recs_written += 1
                    if num_recs is not None and recs_written >= num_recs:
                        break

def download_uniprot_fasta(url, output_file_path, include_taxids, exclude_taxids, ncbi_taxonomy, num_recs=None):
    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        with open(output_file_path, "w") as output_handle:
            with gzip.open(response.raw, "rt") as gzipped_file:
                recs_written = 0
                records = tqdm.tqdm(parse_seqfiles([gzipped_file], filetype_override="fasta"), desc=url.split("/")[-1], leave=True, dynamic_ncols=True)
                if include_taxids or exclude_taxids:
                    records = filter_by_taxonomy(records, include_taxids, exclude_taxids, ncbi_taxonomy)
                for record in records:
                    SeqIO.write(record, output_handle, "fasta")
                    recs_written += 1
                    if num_recs is not None and recs_written >= num_recs:
                        break


def domainator_db_download(ncbi_taxonomy, output_file_path, include_taxids, exclude_taxids, db, gene_call, num_recs, ncbi_summary, before:datetime, after:datetime, cpus:int=1, success_rec_log=None, exclude_accessions:Optional[Set[str]]=None):

    if include_taxids is not None:
        include_taxids = set(include_taxids)
    if exclude_taxids is not None:
        exclude_taxids = set(exclude_taxids)


    if db.lower() in {"swissprot_gb", "uniprot_gb"}:
        download_and_convert_uniprot_dat(
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz",
            output_file_path,
            include_taxids,
            exclude_taxids,
            ncbi_taxonomy,
            num_recs=num_recs,
        )
    if db.lower() in {"trembl_gb", "uniprot_gb"}:
        download_and_convert_uniprot_dat(
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz",
            output_file_path,
            include_taxids,
            exclude_taxids,
            ncbi_taxonomy,
            num_recs=num_recs,
        )
    if db.lower() in {"swissprot", "uniprot"}:
        download_uniprot_fasta(
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
            output_file_path,
            include_taxids,
            exclude_taxids,
            ncbi_taxonomy,
            num_recs=num_recs,
        )
    if db.lower() in {"trembl", "uniprot"}:
        download_uniprot_fasta(
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz",
            output_file_path,
            include_taxids,
            exclude_taxids,
            ncbi_taxonomy,
            num_recs=num_recs,
        )

    if db.lower() == "ncbi_complete_genome_proks":
        if exclude_taxids is None:
            exclude_taxids = set()
        download_ncbi(
            output_file_path,
            include_taxids,
            exclude_taxids.union({2759}), # exclude Eukaryota
            ncbi_taxonomy,
            filter_by_uniqueness=False,
            no_na=False,
            filter_by_complete_genome=True,
            gene_call=gene_call,
            num_recs=num_recs,
            ncbi_summary=ncbi_summary,
            before=before,
            after=after,
            cpus=cpus,
            success_rec_log=success_rec_log,
            exclude_accessions=exclude_accessions,
        )
    if db.lower() == "ncbi_complete_genome_nonredundant_proks":
            if exclude_taxids is None:
                exclude_taxids = set()
            download_ncbi(
                output_file_path,
                include_taxids,
                exclude_taxids.union({2759}), # exclude Eukaryota
                ncbi_taxonomy,
                filter_by_uniqueness=True,
                no_na=False,
                filter_by_complete_genome=True,
                gene_call=gene_call,
                num_recs=num_recs,
                ncbi_summary=ncbi_summary,
                before=before,
                after=after,
                cpus=cpus,
                success_rec_log=success_rec_log,
                exclude_accessions=exclude_accessions,
            )

    if db.lower() == "ncbi_representative_proks":
        if exclude_taxids is None:
            exclude_taxids = set()
        download_ncbi(
            output_file_path,
            include_taxids,
            exclude_taxids.union({2759}), # exclude Eukaryota
            ncbi_taxonomy,
            filter_by_uniqueness=False,
            no_na=True,
            gene_call=gene_call,
            num_recs=num_recs,
            ncbi_summary=ncbi_summary,
            before=before,
            after=after,
            cpus=cpus,
            success_rec_log=success_rec_log,
            exclude_accessions=exclude_accessions,
        )
    if db.lower() == "ncbi_representative_all":
        download_ncbi(
            output_file_path,
            include_taxids,
            exclude_taxids,
            ncbi_taxonomy,
            filter_by_uniqueness=False,
            no_na=False,
            gene_call=gene_call,
            num_recs=num_recs,
            ncbi_summary=ncbi_summary,
            before=before,
            after=after,
            cpus=cpus,
            success_rec_log=success_rec_log,
            exclude_accessions=exclude_accessions,
        )
    if db.lower() == "ncbi_nonredundant_proks":
        if exclude_taxids is None:
            exclude_taxids = set()
        download_ncbi(
            output_file_path,
            include_taxids,
            exclude_taxids.union({2759}),
            ncbi_taxonomy,
            filter_by_uniqueness=True,
            no_na=False,
            gene_call=gene_call,
            num_recs=num_recs,
            ncbi_summary=ncbi_summary,
            before=before,
            after=after,
            cpus=cpus,
            success_rec_log=success_rec_log,
            exclude_accessions=exclude_accessions,
        )
    if db.lower() == "ncbi_nonredundant_all":
        download_ncbi(
            output_file_path,
            include_taxids,
            exclude_taxids,
            ncbi_taxonomy,
            filter_by_uniqueness=True,
            no_na=False,
            gene_call=gene_call,
            num_recs=num_recs,
            ncbi_summary=ncbi_summary,
            before=before,
            after=after,
            cpus=cpus,
            success_rec_log=success_rec_log,
            exclude_accessions=exclude_accessions,
        )
    if db.lower() == "ncbi_all":
        download_ncbi(
            output_file_path,
            include_taxids,
            exclude_taxids,
            ncbi_taxonomy,
            filter_by_uniqueness=False,
            no_na=False,
            gene_call=gene_call,
            num_recs=num_recs,
            ncbi_summary=ncbi_summary,
            before=before,
            after=after,
            cpus=cpus,
            success_rec_log=success_rec_log,
            exclude_accessions=exclude_accessions,
        )

def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)
    parser.add_argument("--db", type=str.lower, choices={"ncbi_complete_genome_proks", "ncbi_complete_genome_nonredundant_proks", "ncbi_representative_proks", "ncbi_representative_all", "ncbi_nonredundant_proks", "ncbi_nonredundant_all", "ncbi_all", "swissprot", "trembl", "uniprot", "swissprot_gb", "trembl_gb", "uniprot_gb"}, help="Database to download")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output filename")
    parser.add_argument("--include_taxids", nargs='+', default=None, type=int, help="Space separated list of taxids to include")
    parser.add_argument("--exclude_taxids", nargs='+', default=None, type=int, help="Space separated list of taxids to exclude")
    parser.add_argument("--ncbi_taxonomy_path", type=str,  default="/tmp/ncbi_taxonomy", help="Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist.")
    parser.add_argument("--skip_taxonomy_update", action="store_true", help="By default, if taxonomy database exists, check it against the version on the ncbi server and update if there is a newer version. Setting this flag will skip this check.")
    parser.add_argument("--success_rec_log", type=str, default=None, help="Path to log file for successfully written records. If not specified, successfully written records will not be logged.")

    parser.add_argument('--exclude_accessions', type=str, default=None, nargs='+',
                        help="Exclude genbank records with these accessions. Should be versionless (i.e. without the '.' and number after the '.')")
    parser.add_argument('--exclude_accessions_file', type=str, default=None,
                        help="text file containing accessions genbank records to exclude. Should be versionless (i.e. without the '.' and number after the '.')")


    parser.add_argument("--ncbi_summary", type=str, default=None, help="Path to NCBI summary file. Normally it will be downloaded from the Internet on every execution, but you may want to pre-download it, or make custom subsets, the url is: https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt.")
    parser.add_argument('--gene_call', type=str, default=None, choices = {"all", "unannotated"}, required=False,
                        help="When activated, for nucleotide databases, such as Genbank, new CDS annotations will be added with Prodigal in Metagenomic mode. If 'all', then any existing CDS annotations will be deleted and all contigs will be re-annotated. If 'unannotated', then only contigs without CDS annotations will be annotated. [default: None].")
    parser.add_argument("--num_recs", type=int, default=None, help="Stop after downloading this many records. For UniProt, this is the number of proteins. For genbank, the number of genomes, including multi-contig genomes. This is mainly for testing.")
    parser.add_argument("--append", action="store_true", help="Append to existing file") 
    parser.add_argument("--before", type=str, default=None, help="For ncbi databases, only download records before this date. Format is YYYY-MM-DD.")
    parser.add_argument("--after", type=str, default=None, help="For ncbi databases, only download records after this date. Format is YYYY-MM-DD.")
    parser.add_argument('--cpu', type=int, default=0,
                        help="the number of cores of the cpu which are used at a time to download ncbi records [default: use all available cores]")

    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)

    if params.append:
        if not Path(params.output).is_file():
            raise Exception(f"File {params.output} is not a file, so you cannot append to it.")
    else: # overwrite the output file with a blank file
        with open(params.output, "w") as output_handle:
            pass
    
        

    exclude_accessions = list_and_file_to_dict_keys(params.exclude_accessions, params.exclude_accessions_file, as_set=True)

    ncbi_taxonomy = None
    if params.include_taxids or params.exclude_taxids or params.db.lower().endswith("_proks"):
        # create the path to the NCBI taxonomy database
        Path(params.ncbi_taxonomy_path).mkdir(parents=True, exist_ok=True)
        # load the NCBI taxonomy database
        ncbi_taxonomy = NCBITaxonomy(params.ncbi_taxonomy_path, not params.skip_taxonomy_update)


    before = None
    if params.before:
        before = datetime.datetime.strptime(params.before, "%Y-%m-%d").date()

    after = None
    if params.after:
        after = datetime.datetime.strptime(params.after, "%Y-%m-%d").date()
    

    # parameter validation
    if params.cpu <= 0:
        cpus = psutil.cpu_count(logical=False)
    else:
        cpus = params.cpu
    

    domainator_db_download(ncbi_taxonomy, params.output, params.include_taxids, params.exclude_taxids, params.db, params.gene_call, params.num_recs, params.ncbi_summary, before, after, cpus, success_rec_log=params.success_rec_log, exclude_accessions=exclude_accessions)

def _entrypoint():
    main(sys.argv[1:])

if __name__ == "__main__":
    main(sys.argv[1:])
