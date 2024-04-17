from math import prod
from domainator.partition_seqfile import main
import tempfile
from glob import glob
from io import StringIO
import pytest

@pytest.mark.parametrize("file,option,value,num_proteins,offsets,recs_to_read",
[("pDONR201.gb","--partitions",2,3,[0],[1]),
("pDONR201.gb","--cdss_per_partition",1,3,[0],[1]),
("pDONR201_multigenemark_partition.gb","--partitions",2,24,[0,16382],[2,2]),
("pDONR201_multigenemark_partition.gb","--partitions",1,24,[0],[4]),
("pDONR201_multigenemark_partition.gb","--partitions",10,24,[0,8191,16382,24573],[1,1,1,1]),
("pDONR201_multigenemark_partition.gb","--partitions",4,24,[0,8191,16382,24573],[1,1,1,1]),
("pDONR201_multigenemark_partition.gb","--cdss_per_partition",6,24,[0,8191,16382,24573],[1,1,1,1]),
("pDONR201_multigenemark_partition.gb","--cdss_per_partition",7,24,[0,16382],[2,2]),
])
def test_partition_seqfile(file, num_proteins, option, value, offsets, recs_to_read, shared_datadir):
    with tempfile.TemporaryDirectory() as output_dir:
        outfile = output_dir + "/" + "outfile.txt"
        args = ["-i", str(shared_datadir / file), "-o", outfile, option, str(value)]
        main(args)
        with open(outfile) as f:
            proteins = int(f.readline().strip())
            assert proteins == num_proteins
            produced_offsets = list()
            produced_recs_to_read = list()
            for line in f:
                parts = line.strip().split()
                produced_offsets.append(int(parts[0]))
                produced_recs_to_read.append(int(parts[1]))
            assert produced_offsets == offsets
            assert produced_recs_to_read == recs_to_read
