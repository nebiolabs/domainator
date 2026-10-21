"""Record offsets in sequence files, and partitioning files for parallel workers.
"""
from array import array
import os

from domainator import db_index
from domainator.seqio import detect_compression, get_file_type, open_seqfile
try:
    from domainator import _gbfast  # optional native (Rust) acceleration
except ImportError:
    _gbfast = None


def i_get_fasta_offsets(input_path): #TODO: replace with cython
    """
        input: a path to a fasta file

        yields: (file_offset, num_proteins)
                note that a fasta file can only have one protein per record, so the num_proteins will always be 1.

    """
    infile, _ = open_seqfile(input_path, binary=True)
    try:
        offset=0
        line = infile.readline()
        while line:
            if line[0] == 62: # 62 = b'>'
                yield offset, 1
            offset = infile.tell()
            line = infile.readline()
    finally:
        infile.close()

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
    with open_seqfile(input_path, binary=True)[0] as infile:
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
    with open_seqfile(input_path, binary=True)[0] as infile:
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
    with open_seqfile(input_path, binary=True)[0] as infile:
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

def _gbfast_usable(input_path):
    """True when a native offset scanner can be used: the extension is installed
    and the file is a real path that is either uncompressed or BGZF. Plain gzip
    (non-blocked) still falls back to the Python scanner."""
    if _gbfast is None:
        return False
    if not isinstance(input_path, (str, os.PathLike)):
        return False
    compression = detect_compression(input_path)
    if compression is None:
        return True
    if compression == "bgzf":
        # BGZF support was added alongside genbank_offsets_bgzf; guard so an older
        # built extension (uncompressed-only) cleanly falls back.
        return hasattr(_gbfast, "genbank_offsets_bgzf")
    return False


def _gbfast_offset_pairs(input_path, filetype):
    """Return native (offset, cds_count) pairs. Offsets are plain byte offsets for
    uncompressed files and BGZF virtual offsets (block_start<<16 | within) for
    BGZF, matching the Python scanners / Bio.bgzf and the parser's seek path."""
    bgzf_input = detect_compression(input_path) == "bgzf"
    if filetype == "fasta":
        scan = _gbfast.fasta_offsets_bgzf if bgzf_input else _gbfast.fasta_offsets
    elif filetype == "genbank":
        scan = _gbfast.genbank_offsets_bgzf if bgzf_input else _gbfast.genbank_offsets
    else:
        raise ValueError(f"Filetype not recognized for input file: {input_path}")
    return scan(str(input_path))


def _index_offsets(input_path, filetype):
    """Return (offsets, num_proteins) from a valid, fresh ``.didx`` sidecar, or None
    to signal the caller to scan. Plain gzip never uses an index (no usable
    random-access offset space)."""
    if not isinstance(input_path, (str, os.PathLike)):
        return None
    compression = detect_compression(input_path)
    if compression == "gzip":
        return None
    return db_index.read_index(input_path, filetype=filetype, compression=compression)


def index_total_cds(input_path):
    """Return the total CDS/protein count from a valid, fresh ``.didx`` header
    (O(1), no scan), or None if there is no usable index for this file."""
    if not isinstance(input_path, (str, os.PathLike)):
        return None
    compression = detect_compression(input_path)
    if compression == "gzip":
        return None
    return db_index.read_total_cds(input_path, filetype=get_file_type(input_path), compression=compression)


def get_offsets(input_path):
    """
        input: a path to a fasta or genbank file

        output: two arrays of unsigned long longs, (offset, num_proteins)
    """

    filetype = get_file_type(input_path)
    indexed = _index_offsets(input_path, filetype)
    if indexed is not None:
        return indexed
    if _gbfast_usable(input_path):
        pairs = _gbfast_offset_pairs(input_path, filetype)
        offsets = array("Q", (p[0] for p in pairs))
        num_proteins = array("Q", (p[1] for p in pairs))
        return offsets, num_proteins

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
    if isinstance(input_path, (str, os.PathLike)):
        compression = detect_compression(input_path)
        if compression != "gzip":
            indexed = db_index.i_read_index(input_path, filetype=filetype, compression=compression)
            if indexed is not None:
                return indexed
    if _gbfast_usable(input_path):
        return iter(_gbfast_offset_pairs(input_path, filetype))
    if filetype == "fasta":
        return i_get_fasta_offsets(input_path) # First array is file offsets of the record, second array is protein count (always 1 for fasta)
    elif filetype == "genbank":
        return i_get_genbank_offsets(input_path) #
    else:
        raise ValueError(f"Filetype not recognized for input file: {input_path}")


def partition_seqfile(input_path, cdss_per_partition):
    """Divide a sequence file into runs of whole records holding about ``cdss_per_partition`` CDSs.

    Args:
        input_path: path to a genbank (peptide or nucleotide) or fasta (peptide) file.
        cdss_per_partition: a partition is closed as soon as it holds at least this many CDSs.

    Returns:
        (total_cds_count, [(input_path, offset, records_to_read), ...]). Seeking to each offset and
        reading that many records reads the whole file exactly once.
    """
    offsets, num_proteins = get_offsets(input_path)
    num_recs = len(offsets)
    if num_recs == 0:
        raise ValueError("Error: No sequence records found in input file (maybe the file has the wrong extension or there is a formatting error in the input file).")

    partitions = list()
    running_sum = 0
    next_offset = offsets[0]
    recs_in_buffer = 0
    for i in range(num_recs):
        running_sum += num_proteins[i]
        recs_in_buffer += 1
        if running_sum >= cdss_per_partition:
            partitions.append((input_path, next_offset, recs_in_buffer))
            running_sum = 0
            recs_in_buffer = 0
            if i + 1 < num_recs:
                next_offset = offsets[i+1]
    if recs_in_buffer > 0:
        partitions.append((input_path, next_offset, recs_in_buffer))

    return sum(num_proteins), partitions


def i_partition_seqfiles(input_paths, cdss_per_partition):
    """Lazily yield the partitions of each file, as :func:`partition_seqfile` would list them.

    Yields:
        (input_path, offset, records_to_read)
    """
    for input_path in input_paths:
        running_sum = 0
        next_offset = None
        recs_in_buffer = 0
        for offset, cds_count in i_get_offsets(input_path):
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
