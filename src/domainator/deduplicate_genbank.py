"""Remove redundant sequences from a genbank file

Runs a clustering algorithm, such as cdhit, usearch, or diamond on sequences from a genbank (or fasta) file to reduce redundancy.
The "hash" algorithm is a faster way to deduplicate exact duplicates, by hashing the sequences and eliminating hash duplicates.
"""

from dataclasses import dataclass, field
import sys
from jsonargparse import ActionConfigFile
from typing import List
from domainator.cli.common import add_alphabet, add_input, add_seq_outputs, check_outputs, reject_stdin, resolve_cpu, tool_parser
from domainator.outputs import write_seq_outputs
from domainator.Bio import SeqIO, Seq
from domainator.utils import append_extra_args, parse_seqfiles, sequence_hashes
import shlex
import tempfile
import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
class ClusterProgram(ABC):
    search_args: List[str]
    input_fasta_path: str
    cluster_output_fasta_path: str
    cluster_assignment_path: str


    def __init__(self, extra_args, bin_path, log_handle=sys.stderr):

        if bin_path is not None:
            self.bin_path = bin_path
        self.log_handle = log_handle
        self.search_args = append_extra_args([self.bin_path] + self.search_args, extra_args, f"--{self.args_option}_args")

    def run(self):
        """

        """

        try:
            out = subprocess.run(self.search_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")
        except FileNotFoundError:
            raise RuntimeError(f"Clustering program not found: {self.search_args[0]}. Install it or give its path with --{self.path_option}_path.") from None

        if out.returncode != 0:
            raise RuntimeError(f"Clustering program failed (exit code {out.returncode}): {' '.join(map(str, self.search_args))}\n"
                               f"{out.stdout or ''}{out.stderr or ''}")
    
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
    args_option = "usearch"
    path_option = "usearch"
    def __init__(self, input_fasta_path, id, extra_args, bin_path=None, cpus=1, log_handle=sys.stderr):
        self.bin_path = "usearch"
        self.input_fasta_path = input_fasta_path
        self.cluster_output_fasta_path = Path(input_fasta_path).parents[0] / "output.fasta"
        self.cluster_assignment_path = Path(input_fasta_path).parents[0] / "output.uc"
        self.search_args = ['-sort', 'length', '-threads', str(int(cpus)), '-cluster_fast', input_fasta_path, "-id", str(id), "-centroids", self.cluster_output_fasta_path, "-uc", self.cluster_assignment_path]
        super().__init__(extra_args, bin_path, log_handle=log_handle)
    
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
    args_option = "cdhit"
    path_option = "cdhit"
    def __init__(self, input_fasta_path, id, extra_args, bin_path=None, cpus=1, word_size=10, log_handle=sys.stderr):
        self.input_fasta_path = input_fasta_path
        self.cluster_output_fasta_path = Path(input_fasta_path).parents[0] / "output.fasta"
        self.cluster_assignment_path = Path(input_fasta_path).parents[0] / "output.fasta.clstr"
        self.word_size = word_size
        self.search_args = ['-n', str(self.word_size), '-M', "0", '-T', str(int(cpus)), "-d", "0", "-c", str(id), "-i", input_fasta_path, "-o", self.cluster_output_fasta_path]
        super().__init__(extra_args, bin_path, log_handle=log_handle)
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
    def __init__(self, input_fasta_path, id, extra_args, bin_path=None, cpus=1, log_handle=sys.stderr):
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

        super().__init__(input_fasta_path, id, extra_args, bin_path=bin_path, cpus=cpus, word_size=word_size, log_handle=log_handle)

class CdHitEST(CdHit):
    def __init__(self, input_fasta_path, id, extra_args, bin_path=None, cpus=1, log_handle=sys.stderr):
        self.bin_path = "cd-hit-est"
        self.path_option = "cdhit_est"
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

        super().__init__(input_fasta_path, id, extra_args, bin_path=bin_path, cpus=cpus, word_size=word_size, log_handle=log_handle)


class Diamond(ClusterProgram):
    args_option = "diamond"
    path_option = "diamond"
    def __init__(self, input_fasta_path, id, extra_args, bin_path=None, cpus=1, log_handle=sys.stderr):
        self.bin_path = "diamond"
        self.input_fasta_path = input_fasta_path
        self.cluster_output_fasta_path = None
        self.cluster_assignment_path = Path(input_fasta_path).parents[0] / "output.tsv"
        self._cluster_members = None
        self.search_args = [
            "cluster",
            "-d",
            input_fasta_path,
            "-o",
            self.cluster_assignment_path,
            "--approx-id",
            str(id * 100),
            "--threads",
            str(int(cpus)),
        ]
        super().__init__(extra_args, bin_path, log_handle=log_handle)

    def get_cluster_members(self):
        # diamond writes one "representative<tab>member" line per sequence, including a line
        # where the representative is its own member, so no post-processing is needed.
        if self._cluster_members is None:
            out = dict()
            with open(self.cluster_assignment_path, "r") as cluster_file:
                for line in cluster_file:
                    parts = line.strip().split("\t")
                    if len(parts) < 2:
                        continue
                    representative, member = parts[0], parts[1]
                    if member.lower() == "member":  # header line, only written with --header
                        continue
                    out.setdefault(representative, list()).append(member)
            self._cluster_members = out
        return self._cluster_members

    def get_reduced_names(self):
        return set(self.get_cluster_members().keys())
   
# class HashDedup(ClusterProgram):
#     pass
# TODO

def run_clustering(algorithm, input_fasta_path, id, extra_args, bin_path=None, add_count=None, cpus=1, log_handle=sys.stderr):
    """
        algorithm: which program to use to cluster

        input_fasta_path: path to a fasta file
        
        id: the identity threshold, as a float >= 0, <= 1.
        
        extra_args: a list of additional command line arguments to pass to the clustering algorithm

        bin_path: alternative path to the binary

        add_count: "prefix", "suffix", or None. If not None, then a dict of cluster_center_name: contigs_in_cluster will be returned.

        cpus: number of cpus to use

        log_handle: where to write log messages

    """

    search_classes = {"cd-hit":CdHitProtein, "usearch": USearch, "cd-hit-est": CdHitEST, "diamond": Diamond}
    
    search = search_classes[algorithm](input_fasta_path, id, extra_args, bin_path, cpus, log_handle=log_handle)

    search.run()

    reduced_names = search.get_reduced_names()

    cluster_members = None
    # if add_count is not None:
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

    def deduplicate_genbank(self, genbanks, algorithm, extra_args, id, bin_paths, add_count, cpus, fasta_type="protein", both_strands=False, log_handle=sys.stderr):
        """
        input: 
                genbanks: list of genbank files to cluster sequences from
                algorithm: which program to use for clustering
                extra_args: a list of additional command line arguments to pass to the clustering algorithm
                bin_paths: None, or a dict of program name ("cd-hit", "cd-hit-est", "usearch", "diamond") to the path of its executable,
                           for programs that are not on the PATH
        """

        #dict of algorithm names to functions that will run the algorithm, the functions must take three arguments, input_fasta_path,cluster_output_fasta_path,params

        if both_strands and fasta_type != "nucleotide":
            raise ValueError("--both_strands is only supported for nucleotide input.")

        if algorithm == "hash": #TODO: multiprocessing
            if id != 1:
                raise ValueError(f"hash algorithm only valid for id = 1, not {id}.")
            seen = dict()
            for i, rec in enumerate(parse_seqfiles(genbanks, default_molecule_type=fasta_type)):
                seq_hashes = sequence_hashes(rec.seq, both_strands)

                hit = None
                for seq_hash in seq_hashes:
                    if seq_hash in seen:
                        hit = seen[seq_hash]
                        break

                if hit is None:
                    hit = HashHit(str(i))
                    for seq_hash in seq_hashes:
                        seen[seq_hash] = hit

                hit.cluster_members.append(str(i))
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
                if algorithm == "diamond" and seqtype != "protein":
                    raise ValueError("diamond deduplication is only supported for protein input.")
                extra_args = list(extra_args or ())
                if both_strands:
                    if algorithm == "usearch" and "-strand" not in extra_args:
                        extra_args += ["-strand", "both"]
                    elif algorithm == "cd-hit-est" and "-r" not in extra_args:
                        extra_args += ["-r", "1"]
                bin_path = (bin_paths or {}).get(algorithm)
                reduced_names, cluster_members = run_clustering(algorithm, input_fasta_path, id, extra_args, bin_path, add_count, cpus, log_handle=log_handle)

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
   
def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, "GenBank or FASTA files.", stdin=False)

    add_seq_outputs(parser)
    
    parser.add_argument('--log_out', default=None, required=False,
                        help="Write the log to this file. If not supplied, writes it to stderr.")

    add_alphabet(parser)

    parser.add_argument("--algorithm", type=str, default="cd-hit", choices={"cd-hit", "usearch", "hash", "diamond"},
                        help="Which clustering algorithm to use.")

    for option, program in (("cdhit", "cd-hit"), ("cdhit_est", "cd-hit-est"), ("usearch", "usearch"), ("diamond", "diamond")):
        parser.add_argument(f"--{option}_path", type=str, default=None,
                            help=f"path to the {program} executable, if it is not on the PATH.")

    parser.add_argument("--identity", type=float, required=True, default=None,
                        help="Identity threshold (between 0 and 1), passed to cd-hit as the -c parameter, to usearch as the -id parameter, or to diamond as --approx-id after converting to percent.")

    labels = parser.add_mutually_exclusive_group(required=False)
    labels.add_argument("--prefix_count", action="store_true",
                        help="In the output genbank file, will prepend X- to the names of contigs, where X is the number of contigs from the original file collapsed into this centroid.")
    labels.add_argument("--suffix_count", action="store_true",
                        help="In the output genbank file, will append -X to the names of contigs, where X is the number of contigs from the original file collapsed into this centroid.")
    
    parser.add_argument('--cpu', type=int, default=0,
                    help="The number of threads to use for deduplication. 0: one per physical core.")

    for option, program in (("cdhit", "cd-hit or cd-hit-est"), ("usearch", "usearch"), ("diamond", "diamond cluster")):
        parser.add_argument(f"--{option}_args", default=None, type=str,
                            help=f"extra options for {program}, as one quoted string split like a shell command line. "
                                 f"Use the = form: --{option}_args=\"-s 0.9\". Options that domainator sets itself are refused.")


    parser.add_argument('--both_strands', action='store_true', default=False,
                        help="When set for nucleotide inputs, consider both strands for clustering/hash deduplication.")
    
    parser.add_argument('--cluster_tsv_out', default=None, required=False,
                        help="If supplied, then write a tab separated table with columns: representative, contigs.")
    
    parser.add_argument('--cluster_sep', default=' ; ', required=False,
                        help="Separator to separate individual contig names in the second column of the cluster table.")

    parser.add_argument('--config', action=ActionConfigFile)
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
   
    cpus = resolve_cpu(params.cpu)
    reject_stdin(parser, params.input)

    ### Figure out what input and output files ####
    check_outputs(parser, params, ["genbank_out", "fasta_out"])
    
    log = sys.stderr if params.log_out is None else open(params.log_out, "w")


    genbanks = params.input


    args_options = {"cd-hit": "cdhit", "usearch": "usearch", "diamond": "diamond"}
    for algorithm, option in args_options.items():
        if algorithm != params.algorithm and getattr(params, f"{option}_args"):
            parser.error(f"--{option}_args was given, but --algorithm is {params.algorithm}.")
    selected_args = getattr(params, f"{args_options[params.algorithm]}_args") if params.algorithm in args_options else None
    extra_args = shlex.split(selected_args) if selected_args else []

    add_count = None
    if params.prefix_count:
        add_count = "prefix"
    elif params.suffix_count:
        add_count = "suffix"

    if params.identity < 0 or params.identity > 1:
        parser.error(f"--identity must be between 0 and 1, not {params.identity}")

    dedup = DeduplicateGenbank()
    # Run
    deduplicate_iterator = dedup.deduplicate_genbank(
            genbanks, 
            params.algorithm, 
            extra_args,
            params.identity,
            {"cd-hit": params.cdhit_path, "cd-hit-est": params.cdhit_est_path, "usearch": params.usearch_path, "diamond": params.diamond_path},
            add_count,
            cpus,
            params.alphabet,
            params.both_strands,
            log_handle=log)

    write_seq_outputs(deduplicate_iterator, genbank=params.genbank_out, fasta=params.fasta_out)
    
    if params.cluster_tsv_out is not None:
        with open(params.cluster_tsv_out, "w") as cluster_table_out:
            print("representative\tcontigs", file=cluster_table_out)
            for representative, contigs in dedup.cluster_members:
                print(f"{representative}\t{params.cluster_sep.join(contigs)}", file=cluster_table_out)
    
    if params.log_out is not None:
        log.close()
