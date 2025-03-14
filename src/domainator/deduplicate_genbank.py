"""Remove redundant sequences from a genbank file

Runs a clustering algorithm, such as cdhit or usearch on sequences from a genbank (or fasta) file to reduce redundancy.
The "hash" algorithm is a faster way to deduplicate exact duplicates, by hashing the sequences and eliminating hash duplicates.
"""

from dataclasses import dataclass, field
import sys
import argparse
from typing import List
from domainator.Bio import SeqIO, Seq
from domainator.utils import parse_seqfiles, write_genbank
import tempfile
import subprocess
import json
from abc import ABC, abstractmethod
from pathlib import Path
import psutil
import hashlib
from domainator import __version__, RawAndDefaultsFormatter

def add_params_to_args_list(arg_list, params):
    for k,v in params.items():
        arg_list.append(k)
        if v is not None:
            arg_list.append(str(v))


class ClusterProgram(ABC):
    search_args: List[str]
    input_fasta_path: str
    cluster_output_fasta_path: str
    cluster_assignment_path: str


    def __init__(self, params, bin_path, log_handle=sys.stderr):

        if bin_path is not None:
            self.bin_path = bin_path
        self.log_handle = log_handle
        self.search_args = [self.bin_path] + self.search_args
        add_params_to_args_list(self.search_args, params)

    def run(self):
        """

        """

        out = subprocess.run(self.search_args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        if out.returncode != 0:
            print(
                f'Error executing clustering program: {self.search_args}', file=self.log_handle)
            exit(1)
    
    @abstractmethod
    def get_cluster_members(self):
        pass

    def get_reduced_names(self):
        #get names of reduced set of sequences
        reduced_names = set()
        for rec in SeqIO.parse(self.cluster_output_fasta_path, "fasta"):
            reduced_names.add(rec.id)
        return reduced_names

class USearch(ClusterProgram):
    def __init__(self, input_fasta_path, id, params, bin_path=None, cpus=1, log_handle=sys.stderr):
        self.bin_path = "usearch"
        self.input_fasta_path = input_fasta_path
        self.cluster_output_fasta_path = Path(input_fasta_path).parents[0] / "output.fasta"
        self.cluster_assignment_path = Path(input_fasta_path).parents[0] / "output.uc"
        self.search_args = ['-sort', 'length', '-threads', str(int(cpus)), '-cluster_fast', input_fasta_path, "-id", str(id), "-centroids", self.cluster_output_fasta_path, "-uc", self.cluster_assignment_path]
        super().__init__(params, bin_path, log_handle=log_handle)
    
    def get_cluster_members(self):
        out = dict()
        with open(self.cluster_assignment_path, "r") as cluster_file:
            for line in cluster_file:
                line=line.strip()
                parts = line.split("\t")
                if parts[0] == "H":
                    if parts[9] not in out:
                        out[parts[9]] = list()
                    out[parts[9]].append(parts[8])
                elif parts[0] == "S":
                    if parts[8] not in out:
                        out[parts[8]] = list()
                    out[parts[8]].append(parts[8])
        return out


class CdHit(ClusterProgram, ABC):
    def __init__(self, input_fasta_path, id, params, bin_path=None, cpus=1, word_size=10, log_handle=sys.stderr):
        self.input_fasta_path = input_fasta_path
        self.cluster_output_fasta_path = Path(input_fasta_path).parents[0] / "output.fasta"
        self.cluster_assignment_path = Path(input_fasta_path).parents[0] / "output.fasta.clstr"
        self.word_size = word_size
        self.search_args = ['-n', str(self.word_size), '-M', "0", '-T', str(int(cpus)), "-d", "0", "-c", str(id), "-i", input_fasta_path, "-o", self.cluster_output_fasta_path]
        super().__init__(params, bin_path, log_handle=log_handle)
    def get_cluster_members(self):
        out = dict()
        with open(self.cluster_assignment_path, "r") as cluster_file:
            center = None
            members = []
            for line in cluster_file:
                if line[0] == ">":
                    if center != None:
                        out[center] = members
                    
                    center = None
                    members = []
                else:
                    parts = line.split()
                    if parts[3] == "*": # cluster center
                        center = parts[2][1:-3]
                    members.append(parts[2][1:-3])

            out[center] = members
        return out


class CdHitProtein(CdHit):
    def __init__(self, input_fasta_path, id, params, bin_path=None, cpus=1, log_handle=sys.stderr):
        self.bin_path = "cd-hit"
        word_size = 5
        if id > 0.7:
            word_size = 5
        elif  id > 0.6:
            word_size = 4
        elif id > 0.5:
            word_size = 3
        elif id >= 0.4:
            word_size = 2
        else:
            raise ValueError("For cd-hit on protein sequences, id must be >= 0.4")

        super().__init__(input_fasta_path, id, params, bin_path=bin_path, cpus=cpus, word_size=word_size, log_handle=log_handle)

class CdHitEST(CdHit):
    def __init__(self, input_fasta_path, id, params, bin_path=None, cpus=1, log_handle=sys.stderr):
        self.bin_path = "cd-hit-est"
        word_size = 10
        if id > 0.95:
            word_size = 10
        elif  id > 0.9:
            word_size = 9
        elif id > 0.88:
            word_size = 7
        elif id > 0.85:
            word_size = 6
        elif id >= 0.8:
            word_size = 5
        else:
            raise ValueError("For cd-hit on nucleotide sequences, id must be >= 0.8")

        super().__init__(input_fasta_path, id, params, bin_path=bin_path, cpus=cpus, word_size=word_size, log_handle=log_handle)
   
# class HashDedup(ClusterProgram):
#     pass
# TODO

def run_clustering(algorithm, input_fasta_path, id, params, bin_path=None, add_count=None, cpus=1, log_handle=sys.stderr):
    """
        algorithm: which program to use to cluster

        input_fasta_path: path to a fasta file
        
        id: the identity threshold, as a float >= 0, <= 1.
        
        params: a dict of additional command line arguments to pass to the clustering algorithm

        bin_path: alternative path to the binary

        add_count: "prefix", "suffix", or None. If not None, then a dict of cluster_center_name: contigs_in_cluster will be returned.

        cpus: number of cpus to use

        log_handle: where to write log messages

    """

    search_classes = {"cd-hit":CdHitProtein, "usearch": USearch, "cd-hit-est": CdHitEST}
    
    search = search_classes[algorithm](input_fasta_path, id, params, bin_path, cpus, log_handle=log_handle)

    search.run()

    reduced_names = search.get_reduced_names()

    cluster_members = None
    if add_count is not None:
        cluster_members = search.get_cluster_members()


    return reduced_names, cluster_members
    
@dataclass
class HashHit():
    rep_name: str
    cluster_members: list[str] = field(default_factory=list)

class DeduplicateGenbank():
    def __init__(self):
        self._cluster_members = None

    @property
    def cluster_members(self):
        if self._cluster_members is None:
            raise ValueError("cluster_members not set. Please run deduplicate_genbank first.")
        else:
            return self._cluster_members

    def deduplicate_genbank(self, genbanks, algorithm, params, id, bin_path, add_count, cpus, fasta_type="protein", log_handle=sys.stderr):
        """
        input: 
                genbanks: list of genbank files to cluster sequences from
                algorithm: which program to use for clustering
                params: a dict of additional command line arguments to pass to the clustering algorithm
        """

        #dict of algorithm names to functions that will run the algorithm, the functions must take three arguments, input_fasta_path,cluster_output_fasta_path,params

        if algorithm == "hash": #TODO: multiprocessing
            if id != 1:
                raise ValueError(f"hash algorithm only valid for id = 1, not {id}.")
            seen = dict()
            for i, rec in enumerate(parse_seqfiles(genbanks, default_molecule_type=fasta_type)):
                seq_hash = hashlib.sha256(str(rec.seq).upper().strip('*').encode("utf-8")).hexdigest()
                if seq_hash not in seen:
                    seen[seq_hash] = HashHit(str(i))
                seen[seq_hash].cluster_members.append(str(i))
            reduced_names = {x.rep_name for x in seen.values()}
            cluster_members = {x.rep_name: x.cluster_members for x in seen.values()}
                

        else: #not hash algorithm
            with tempfile.TemporaryDirectory() as output_dir:
                # output_dir = "test_out"
                input_fasta_path = output_dir + "/input.fasta"

                seqtype = None
                #convert input to fasta
                with open(input_fasta_path, "w") as input_handle:
                    for i, rec in enumerate(parse_seqfiles(genbanks, default_molecule_type=fasta_type)):
                        if seqtype is None:
                            seqtype = rec.annotations['molecule_type']
                        rec.id = str(i)
                        rec.description = ""
                        rec.name = ""
                        rec.seq = Seq.Seq(str(rec.seq).upper().replace("*", "f"))
                        SeqIO.write([rec],input_handle,"fasta") 
                
                #run clustering algorithm
                if algorithm == "cd-hit" and seqtype != "protein":
                    algorithm = "cd-hit-est"
                reduced_names, cluster_members = run_clustering(algorithm, input_fasta_path, id, params, bin_path, add_count, cpus, log_handle=log_handle)

        #write reduced set of sequences to a genbank file
        input_ct = 0
        idx_to_original_name = list()
        idx_to_center_name = dict()

        for i, rec in enumerate(parse_seqfiles(genbanks, default_molecule_type=fasta_type)):
            i_str = str(i)
            input_ct += 1
            idx_to_original_name.append(rec.id)
            if i_str in reduced_names:
                if add_count == "prefix":
                    rec.id = f"{len(cluster_members[i_str])}-" + rec.id
                elif add_count == "suffix":
                    rec.id = rec.id + f"-{len(cluster_members[i_str])}"
                idx_to_center_name[i] = rec.id
                yield rec
        
        self._cluster_members = list()
        for center_i_str, members_i_strs in cluster_members.items():
            center_i = int(center_i_str)
            member_names = [idx_to_original_name[int(i)] for i in members_i_strs]
            self._cluster_members.append((idx_to_center_name[int(center_i)], member_names))


        print(f"Input  size: {input_ct}\nOutput size: {len(reduced_names)}", file=log_handle)
   
def main(argv):
    parser = argparse.ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', nargs='+', default=None, required=True, #reads the input multiple times, so it needs to be a file.
                        help="Genbank or fasta filenames.")

    parser.add_argument('-o', '--output', default=None, required=False,
                        help="The name of the output file. If not supplied, writes to stdout.")
    
    parser.add_argument('--log', default=None, required=False,
                        help="The name of the log file. If not supplied, writes to stderr.")

    parser.add_argument("--fasta_type", type=str, default="protein", choices={"protein", "nucleotide"}, 
                        help="Whether the sequences in fasta files are protein or nucleotide sequences.")

    parser.add_argument("--algorithm", type=str, default="cd-hit", choices={"cd-hit", "usearch", "hash"},
                        help="Which clustering algorithm to use.")

    parser.add_argument("--bin_path", type=str, default=None,
                        help="If the executable for the algorithm you're using isn't in the system path, you can provide the full path to the executable here.")

    parser.add_argument("--id", type=float, required=True, default=None,
                        help="Identity threshold (between 0 and 1), passed to cd-hit as the -c parameter or to usearach as the -id parameter.")

    labels = parser.add_mutually_exclusive_group(required=False)
    labels.add_argument("--prefix_count", action="store_true",
                        help="In the output genbank file, will prepend X- to the names of contigs, where X is the number of contigs from the original file collapsed into this centroid.")
    labels.add_argument("--suffix_count", action="store_true",
                        help="In the output genbank file, will append -X to the names of contigs, where X is the number of contigs from the original file collapsed into this centroid.")
    
    parser.add_argument('--cpu', type=int, default=0,
                    help="The number of threads to use for deduplication [default: use all available cores]")

    #TODO: make the syntax less weird, like removing requirement for quotes, or add common options as distinct command line parameters
    parser.add_argument("--params", default=None,
                        help="string of parameters to pass to the clustering algorithm, in the form of a json dict in single quotes. example: '\"-s\":0.9' ")

    parser.add_argument('--fasta_out', action='store_true', default=False,
                        help="makes output a fasta file when activated")
    
    parser.add_argument('--cluster_table', default=None, required=False,
                        help="If supplied, then write a tab separated table with columns: representative, contigs.")
    
    parser.add_argument('--cluster_sep', default=' ; ', required=False,
                        help="Separator to separate individual contig names in the second column of the cluster table. Default: ' ; '")

    params = parser.parse_args(argv)
   
    if params.cpu <= 0:
        cpus = psutil.cpu_count(logical=False)
    else:
        cpus = params.cpu

    ### Figure out what input and output files ####
    if params.output is None:
        out = sys.stdout
    else:
        out = open(params.output, "w")
    
    if params.log is None:
        log = sys.stderr
    else:
        log = open(params.log, "w")


    genbanks = params.input


    algorithm_parameters = {}
    if params.params:
        new_params = "{" + params.params + "}" 
        
        algorithm_parameters = json.loads(new_params)

    add_count = None
    if params.prefix_count:
        add_count = "prefix"
    elif params.suffix_count:
        add_count = "suffix"

    if params.id < 0 or params.id > 1:
        raise ValueError(f"id must be between 0 and 1, not {params.id}")

    dedup = DeduplicateGenbank()
    # Run
    deduplicate_iterator = dedup.deduplicate_genbank(
            genbanks, 
            params.algorithm, 
            algorithm_parameters, 
            params.id, 
            params.bin_path,
            add_count,
            cpus,
            params.fasta_type,
            log_handle=log)

    if params.fasta_out:
        SeqIO.write(deduplicate_iterator, out, "fasta")
    else:
        write_genbank(deduplicate_iterator, out, default_molecule_type=params.fasta_type)
    
    if params.cluster_table is not None:
        with open(params.cluster_table, "w") as cluster_table_out:
            print("representative\tcontigs", file=cluster_table_out)
            for representative, contigs in dedup.cluster_members:
                print(f"{representative}\t{params.cluster_sep.join(contigs)}", file=cluster_table_out)

    if params.output is not None:
        out.close()
    
    if params.log is not None:
        log.close()

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
