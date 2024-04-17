from domainator.partition_seqids import partition_seqids
import tempfile
from glob import glob
from io import StringIO
import pytest

@pytest.mark.parametrize("partitions, ids_per_partition, rec_count",
[(1, None, 10),
(3, None, 10),
(6, None, 10),
(9, None, 10),
(10, None, 10),
(None, 5, 10),
(None, 9, 10),
(None, 2, 10),
#  (1, None, ":memory:", 10),
#  (1, None, ":memory:", 10),
])
def test_partition_seqids_fasta(partitions, ids_per_partition, rec_count, shared_datadir, capsys):
    names = [f"seq{x}" for x in range(rec_count)]
    f = StringIO("\n".join([f">{x}\nMAGICCATS" for x in names]))
    # (input_path, output_prefix, partitions, ids_per_partition, index_path=None, file_format="fasta")
    if partitions is not None:
        intended_partitions = partitions
    else:
        intended_partitions = int(rec_count / ids_per_partition)
        if rec_count % ids_per_partition != 0:
            intended_partitions += 1
    with tempfile.TemporaryDirectory() as output_dir:
        output_prefix = output_dir + "/out"
        partition_seqids([f], output_prefix, partitions, ids_per_partition, filetype="fasta")
        captured = capsys.readouterr()
        assert captured.out == str(rec_count)
        assert len(glob(output_prefix+"*.txt")) <= intended_partitions
        assert len(glob(output_prefix+"*.txt")) >= int(intended_partitions/2)
        ids_in_files = set()
        bins_with_too_few = 0
        for f in glob(output_prefix+"*.txt"):
            ids_in_file = 0
            with open(f) as inf:
                for line in inf:
                    ids_in_file += 1
                    line = line.rstrip()
                    ids_in_files.add(line)
            if ids_per_partition is not None:
                assert ids_in_file <= ids_per_partition
                if ids_in_file < ids_per_partition:
                    bins_with_too_few += 1

        assert bins_with_too_few <= 1
        assert len(ids_in_files) == len(names)
        assert len(ids_in_files.intersection(set(names))) == len(ids_in_files)

def test_partition_seqids_genbank(shared_datadir, capsys):
    with tempfile.TemporaryDirectory() as output_dir:
        
        output_prefix = output_dir + "/out"
        partition_seqids([str(shared_datadir / "pDONR201_multi_genemark_domainator.gb")], output_prefix, 2, None, filetype="genbank")
        captured = capsys.readouterr()
        assert captured.out == "24"
        assert len(glob(output_prefix+"*.txt")) == 2


def test_partition_seqids_genbank_peptide(shared_datadir, capsys):
    with tempfile.TemporaryDirectory() as output_dir:
        
        output_prefix = output_dir + "/out"
        partition_seqids([str(shared_datadir / "FeSOD_20.gb")], output_prefix, 2, None, filetype="genbank")
        captured = capsys.readouterr()
        assert captured.out == "20"
        assert len(glob(output_prefix+"*.txt")) == 2
