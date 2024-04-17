"""Partitions the sequence IDs from a sequence file

Given an input sequence file, write new text files consisting of the sequence IDs partitioned into groups.
Also prints the total number of records in the file (without a newline at the end).
The "number of records" is the number of CDSs in nucleotide databases, and the number of contigs in protein databases.

"""
#TODO: maybe deprecate this whole file?

#TODO: change the purpose of the file to split based on some kind of annotation, like SSN_cluster, then it could be piped into a script for generating group-specific HMMS.

import argparse
from domainator.utils import parse_seqfiles, count_peptides_in_record
from domainator import __version__, RawAndDefaultsFormatter
import sys

def partition_seqids(input_path, output_prefix, partitions, ids_per_partition, filetype=None):
    cds_count = 0
    rec_names = list()
    
    if partitions is None and ids_per_partition is None:
        raise ValueError("Error: Must specify either partitions or ids_per_partition")
    
    if partitions is not None and ids_per_partition is not None:
        raise ValueError("Error: Must specify only one of partitions or ids_per_partition")

    for rec in parse_seqfiles(input_path, filetype_override=filetype):
        cds_count += count_peptides_in_record(rec)
        rec_names.append(rec.id)

    num_recs = len(rec_names)

    if num_recs == 0:
        raise ValueError("Error: No sequence records found in input file (maybe you specified the wrong file type or there is a formatting error in the input file).")

    if partitions is None: #partitions is None, so we need to set it based on num_recs and ids_per_partition
        partitions = int(num_recs / ids_per_partition)
        if num_recs % ids_per_partition != 0: #  if rec ids don't divide evenly into ids_per_partition we need an extra partition
            partitions += 1
    if ids_per_partition is None: #ids_per_partition is None, so we need to set it based on num_recs and partitions
        ids_per_partition = int(num_recs / partitions)
        if num_recs % partitions != 0: # if rec ids don't divide evenly into partitions we need an extra partition
            ids_per_partition += 1
    max_digits = len(f"{partitions}")
    out_index = 0
    out_file = None

    for position, id in enumerate(rec_names):
        
        if position % ids_per_partition == 0:
            out_index += 1
            if out_file is not None:
                out_file.close()
            out_file = open(output_prefix + str(out_index).zfill(max_digits)+ ".txt", "w")
        print(id, file=out_file)
    
    out_file.close()

    print(cds_count,end='')

def main(argv):
    parser = argparse.ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', default=None,
                       nargs='+', type=str,
                       help="the genbank or fasta files to split the contig names of. Genbank files can be nucleotide (with CDS annotations) or peptide. Fasta files must be peptide.")
    parser.add_argument('--output_prefix', required=True, type=str,
                        help="Output files will be named [output_prefix][0-9]+.txt")

    overwrite_group = parser.add_mutually_exclusive_group(required=True)
    overwrite_group.add_argument('--partitions', type=int, default=None,
                        help="The number of partitions to divide the ids into, roughly evenly. Actual number of partitions will usually be a little smaller than this number due to bin rounding.")
    overwrite_group.add_argument('--ids_per_partition', type=int, default=None,
                        help="The number of ids to write to each partition.")
    params = parser.parse_args(argv)


    partition_seqids(params.input, params.output_prefix, params.partitions, params.ids_per_partition)

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])