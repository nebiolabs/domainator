"""Find record offsets in a sequence file

Given an input sequence file, writes the total count of CDSs/protein sequences,
followed by a list of offset\trecords pairs to divide the sequence file.

Seeking to each offset and reading the specified number of records will result in reading the entire file once.
"""
from jsonargparse import ArgumentParser, ActionConfigFile
import sys
from domainator import __version__, RawAndDefaultsFormatter
from domainator import utils


def i_partition_seqfiles(input_paths, cdss_per_partition): #TODO: test this!
    """
        input:
            input_paths: a list of paths to genbank (peptide or nucleotide) or fasta (peptide) files.
            cdss_per_partition: how many cdss to try to include in each partition (mutually exclusive with partitions)

        yields:
            (input_path, offset, recs_to_read)
    
    """

    for input_path in input_paths:
        running_sum = 0
        next_offset = None
        recs_in_buffer = 0
        for offset, cds_count in utils.i_get_offsets(input_path):
            if next_offset is None:
                next_offset = offset
            running_sum += cds_count
            recs_in_buffer += 1
            if running_sum >= cdss_per_partition:
                yield (input_path, next_offset, recs_in_buffer)
                running_sum = 0
                recs_in_buffer = 0
                next_offset = None
        if recs_in_buffer > 0:
            yield (input_path, next_offset, recs_in_buffer)


def partition_seqfile(input_path, partitions=None, cdss_per_partition=None):
    """
        input:
            input_path: path to a genbank (peptide or nucleotide) or fasta (peptide) file.
            partitions: desired number of partitions to divide the file into 
            cdss_per_partition: how many cdss to try to include in each partition (mutually exclusive with partitions)

        output:
            total_cds_count, [(input_path, offset, recs_to_read)]
    """

    cds_count = 0
    out_list = list()

    if partitions is None and cdss_per_partition is None:
        raise ValueError("Error: Must specify either partitions or cdss_per_partition")

    if partitions is not None and cdss_per_partition is not None:
        raise ValueError("Error: Must specify only one of partitions or cdss_per_partition")
    
    offsets, num_proteins = utils.get_offsets(input_path)

    num_recs = len(offsets)

    cds_count = sum(num_proteins)
    
    if num_recs == 0:
        raise ValueError("Error: No sequence records found in input file (maybe the file has the wrong extension or there is a formatting error in the input file).")


    if cdss_per_partition is None:
        cdss_per_partition = int(cds_count/partitions)

    running_sum = 0
    next_offset = offsets[0]
    recs_in_buffer = 0
    for i in range(num_recs):
        running_sum += num_proteins[i]
        recs_in_buffer += 1
        if running_sum >= cdss_per_partition:
            out_list.append((input_path, next_offset, recs_in_buffer))
            running_sum = 0
            recs_in_buffer = 0
            if i + 1 < num_recs:
                next_offset = offsets[i+1] # next_offset = the offset of the next record
    if recs_in_buffer > 0:
        out_list.append((input_path, next_offset, recs_in_buffer))

    return cds_count, out_list


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', default=None, type=str, required=True,
                       help="the genbank or fasta file to split the contig names of. Genbank files can be nucleotide (with CDS annotations) or peptide. Fasta files must be peptide.")
    parser.add_argument('-o',"--output", default=None, type=str,
                        help="Name of output file. Default: stdout")

    overwrite_group = parser.add_mutually_exclusive_group(required=True)
    overwrite_group.add_argument('--partitions', type=int, default=None,
                        help="The number of partitions to divide the ids into, roughly evenly.")
    overwrite_group.add_argument('--cdss_per_partition', type=int, default=None,
                        help="The approximate number of ids to write to each partition. Partitioning algorithm is greedy, it adds records until the CDS count is met or exceeded, then goes to the next start pointer.")

    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)

    if params.output is None:
        out = sys.stdout
    else:
        out = open(params.output, "w")

    
    cds_count, splits = partition_seqfile(params.input, params.partitions, params.cdss_per_partition)
    print(cds_count, file=out)
    for _, offset, recs in splits:
        print(f"{offset}\t{recs}", file=out)
    
    if params.output is not None:
        out.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
