"""Calculates similarity metrics for gene neighborhoods

Takes a genbank file annotated with domainator and calculates similarities between contigs and writes reports.

ji, ai, and dss metrics are similar to and inspired by the metrics used in the BigScape method:
Navarro-Muñoz, Jorge C., Nelly Selem-Mojica, Michael W. Mullowney, Satria A. Kautsar, James H. Tryon, Elizabeth I. Parkinson, Emmanuel L. C. De Los Santos, et al. “A Computational Framework to Explore Large-Scale Biosynthetic Diversity from Large-Scale Genomic Data.” Nature Chemical Biology 16, no. 1 (January 2020): 60–68. https://doi.org/10.1038/s41589-019-0400-9.
"""
import warnings
warnings.filterwarnings("ignore", module='numpy')
import sys
import argparse
from domainator.utils import list_and_file_to_dict_keys, parse_seqfiles, write_genbank, copy_SeqRecord, get_file_type
from abc import ABC, abstractmethod
from collections import OrderedDict
import scipy.spatial
import scipy.cluster
import math
from domainator import __version__, DOMAIN_FEATURE_NAME, RawAndDefaultsFormatter
from domainator.data_matrix import DataMatrix
from domainator.filter_domains import filter_domains
from typing import Optional, Set
import numpy as np
import scipy.sparse

#TODO: add threading

class ContigMetric(ABC):
    weight: float = 0.0


    @abstractmethod
    def compute(self, recs, k):
        """
            computes the metric for each of the contigs, returns an iterable of tuples:
            (query, target, score)
            
            k: return scores for this many of the top hits for each query, otherwise return all scores. (zero-scores are always omitted)
        """
        pass

class JaccardIndex(ContigMetric):
    def __init__(self, weight: float = 0.0, databases=None):
        self.weight = weight
        self.databases = databases
    
    def get_contig_data(self, recs, databases=None):
        """
            recs: a list of SeqRecords
            databases: a set of database names to keep, if None, keep all databases

            store the domain content for the contig
        """
        out = OrderedDict()
        for contig in recs:
            if contig.id not in out:
                out[contig.id] = set()
            
            for feature in contig.features:
                if feature.type == DOMAIN_FEATURE_NAME:
                    if databases is not None and feature.qualifiers["database"][0] not in databases:
                        continue
                    domain_name = feature.qualifiers["name"][0]
                    out[contig.id].add(domain_name)

        return out


    def compute(self, recs, k):
        contig_data = self.get_contig_data(recs, self.databases)
        for query_id, query_domains in contig_data.items():
            tmp_out = list()
            for target_id, target_domains in contig_data.items():
                if target_id == query_id:
                    tmp_out.append((query_id, target_id, 1.0))
                    continue
                p_or_q = len(query_domains.union(target_domains))
                p_and_q = len(query_domains.intersection(target_domains))
                if p_or_q == 0:
                    ji = 1.0 # if both contigs have no domains, then they are 100% similar
                else:
                    ji = p_and_q / p_or_q
                if ji > 0:
                    tmp_out.append((query_id,target_id,ji))
            if k is not None:
                tmp_out.sort(key=lambda x: -x[2])
                for i in range(k):
                    yield(tmp_out[i])
            else:
                for i in range(len(tmp_out)):
                    yield tmp_out[i]


class AdjacencyIndex(ContigMetric):
    def __init__(self, weight: float = 0.0, databases=None):
        self.weight = weight
        self.databases = databases
        #self.contig_data = OrderedDict()
    
    def get_contig_data(self, recs, databases=None):
        """
            store the domain content for the contig
        """
        out = OrderedDict()
        for contig in recs:
            if contig.id not in out:
                out[contig.id] = set()
            
            contig.features.sort(key=lambda x: x.location.start) #sort features by start location
            last_feature = None
            for feature in contig.features:
                if feature.type == DOMAIN_FEATURE_NAME:
                    if databases is not None and feature.qualifiers["database"][0] not in databases:
                        continue
                    domain_name = feature.qualifiers["name"][0]
                    if last_feature is not None:
                        #out[contig.id].add(tuple(sorted([domain_name,last_feature])))  
                        out[contig.id].add(tuple([domain_name, last_feature])) # If we don't sort alphabetically, we are assuming that the contigs are oriented in the same way, which might be a bad assumption?
                    last_feature = domain_name
        return out
    
    def compute(self, recs, k):
        contig_data = self.get_contig_data(recs, self.databases)
        for query_id, query_domains in contig_data.items():
            tmp_out = list()
            for target_id, target_domains in contig_data.items():
                if target_id == query_id:
                    tmp_out.append((query_id,target_id,1.0))
                    continue
                p_or_q = len(query_domains.union(target_domains))
                p_and_q = len(query_domains.intersection(target_domains))
                if p_or_q == 0:
                    ai = 1.0
                else:
                    ai = p_and_q / p_or_q
                if ai > 0:
                    tmp_out.append((query_id,target_id,ai))
            if k is not None:
                tmp_out.sort(key=lambda x: -x[2])
                for i in range(k):
                    yield(tmp_out[i])
            else:
                for i in range(len(tmp_out)):
                    yield tmp_out[i]
        
class DomainSequenceSimilarity(ContigMetric):
    pass


class DistanceReport(ABC):
    outpath: str 

    @abstractmethod
    def write(self, recs, scores, options):
        """
            recs: a list of SeqRecords
            scores: a dok_sparse matrix of scores representing how similar each pair of contigs are
            options: a dict of additional options
        """
        pass

class SparseTableOutput(DistanceReport):
    def __init__(self, outpath):
        self.outpath = outpath
    
    def write(self, recs, scores, options):
        row_names = [x.id for x in recs]
        DataMatrix.write_sparse(scores, self.outpath, row_names, row_names, data_type="score")


class ClinkerOutput(DistanceReport):

    def __init__(self, outpath):
        self.outpath = outpath
    
    def write(self, recs, scores, options):
        warnings.warn("Clinker output not implemented")

class SvgTreeOutput(DistanceReport):

    def __init__(self, outpath):
        self.outpath = outpath
    
    def write(self, recs, scores, options):
        warnings.warn("SVG Tree output not implemented")

class GenbankOutput(DistanceReport):

    def __init__(self, outpath):
        self.outpath = outpath
    
    def hierarchical_cluster(self, scores):
        """
            runs hierarchical clustering on score matrix, returns a list of record indices optimally ordered by similarity
        """
        
        #TODO: for really large datasets it might be good to find a way to cluster on sparse matrices.
        dense = scores.toarray()
        dense = 1 - (dense/np.max(dense)) #convert scores to distances, #TODO: could also cluster on the euclidean distance of score matrix?
        np.fill_diagonal(dense,0)
        if np.max(dense) == 0: #if all scores are zero, just return the original order
            return list(range(dense.shape[0]))
        dense = scipy.spatial.distance.squareform(dense)
        linkage_matrix = scipy.cluster.hierarchy.linkage(dense, method='ward',  optimal_ordering=True)
        out = list()
        for r in range(len(linkage_matrix)):
            if linkage_matrix[r,0] < scores.shape[0]:
                out.append(int(linkage_matrix[r,0]))
            if linkage_matrix[r,1] < scores.shape[0]:
                out.append(int(linkage_matrix[r,1]))
        return out

    def write(self, recs, scores, options):
        if len(recs) > 0:
            order = self.hierarchical_cluster(scores)
            digits = math.ceil(math.log10(len(order)+1))
            with open(self.outpath, "w") as outfile:
                print_i = 0
                for i in order:
                    rec = recs[i]
                    if "name_by_order" in options:
                        copy_SeqRecord(rec)
                        rec.id = str(print_i).zfill(digits) + "_" + rec.id
                        print_i += 1
                    write_genbank([rec], outfile)
        else:
            with open(self.outpath, "w") as outfile:
                pass


class DenseTextTableOutput(DistanceReport):

    def __init__(self, outpath):
        self.outpath = outpath
    
    def write(self, recs, scores, options):
        row_names = [x.id for x in recs]
        DataMatrix.write_dense_text(scores, self.outpath, row_names, row_names)

class DenseTableOutput(DistanceReport):

    def __init__(self, outpath):
        self.outpath = outpath
    
    def write(self, recs, scores, options):
        row_names = [x.id for x in recs]

        DataMatrix.write_dense(scores.toarray(), self.outpath, row_names, row_names, data_type="score")

def compare_contigs(genbanks, metrics, reports, k, contigs, name_by_order, databases:Optional[Set[str]]=None):
    
    report_options = dict()
    if name_by_order:
        report_options["name_by_order"] = True


    name_to_idx = dict()
    recs = list()
    #collect data needed for the metrics
    for rec in parse_seqfiles(genbanks, contigs):
        if rec.id not in name_to_idx:
            name_to_idx[rec.id] = len(recs)
            if databases is not None:
                rec = tuple(filter_domains((rec,), evalue=float("inf"), max_overlap=1, databases_keep=databases))[0]
            recs.append(rec)
        else:
            warnings.warn(f"Duplicate contig names !!! :{rec.id}")
            
    results_matrix = scipy.sparse.dok_array((len(recs), len(recs)),dtype=np.float64) #lil_matrix would probably also be fine, maybe a little faster in some cases
    
    # compute metrics
    for metric in metrics:
        for query, target, value in metric.compute(recs, k):
            results_matrix[name_to_idx[query],name_to_idx[target]] +=  metric.weight * value

    # write reports
    for report in reports:
        report.write(recs, results_matrix, report_options)


def main(argv):
    parser = argparse.ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', nargs='+', default=None, required=True, help="Genbank filenames.")

    parser.add_argument('-o', '--output', default=None, required=False,
                        help="Name of genbank output, contigs will be sorted hierarchically within the output genbank. If not supplied then no genbank output will be generated. To write to stdout, use '-'. default: None")
    
    parser.add_argument("--name_by_order", action="store_true", 
                        help="Will rename contigs in the genbank file such that they have a prefix that can be sorted to restore the sorted order.")

    parser.add_argument('--contigs', default=[], nargs='+',
                        help="only annotate contigs with ids in this list. Additive with --contigs_file. default: all contigs.")

    parser.add_argument('--contigs_file', default=None, 
                        help="only annotate contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs.")

    parser.add_argument('-k', default=None, type=int,
                        help="for each metric, return scores for this many of the top hits, and consider all other hits to be zero. default: count all comparisons.")

    # Metrics
    parser.add_argument("--ji", default=0.5, required=False, type=float,
                        help="weighting for jaccard index metric, a comparison of all distinct types of domains in the contigs. When two contigs have no domains, their ji is 1.0.")
    parser.add_argument("--ai", default=0.5, required=False, type=float,
                        help="weighting for adjacency index metric, a comparison of shared domain pairs in the contigs. When two contigs have 0 or 1 domains, their ai is 1.0.")
#   parser.add_argument("--dss", default=0.0, required=False,
#                        help="weighting for domain sequence similarity, a comparison of shared domain pairs in the contigs.")
    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Domain databases to use for domain content comparison. default: all databases.")

    # Output
    # parser.add_argument("--clinker", type=str, default=None, help="Write an interactive html in clinker format to this path.")

    # parser.add_argument("--svg_tree", type=str, default=None, help="Write a tree graphic in svg format to this path.")
    
    parser.add_argument("--dense", type=str, default=None, help="Write a dense distance matrix hdf5 file to this path.")
    
    parser.add_argument("--dense_text", type=str, default=None, help="Write a dense distance matrix tsv file to this path.")
    
    parser.add_argument("--sparse", type=str, default=None, help="Write a sparse distance matrix hdf5 file to this path.")


    params = parser.parse_args(argv)
   
    ### Figure out what input and output files ####

    # if pipe_in:
    #     genbanks = [sys.stdin]
    # elif params.genbank is not None:
    genbanks = params.input


    reports = list()
    # if params.clinker is not None:
    #     reports.append(ClinkerOutput(params.clinker))
    # if params.svg_tree is not None:
    #     reports.append(SvgTreeOutput(params.svg_tree))
    if params.dense is not None:
        if get_file_type(params.dense) != "hdf5":
            raise ValueError("Please use an hdf5 related extension for the --dense output, such as .h5, .hdf5, or .hdf.")
        reports.append(DenseTableOutput(params.dense))
    if params.dense_text is not None:
        reports.append(DenseTextTableOutput(params.dense_text))
    if params.sparse is not None:
        if get_file_type(params.sparse) != "hdf5":
            raise ValueError("Please use an hdf5 related extension for the --sparse output, such as .h5, .hdf5, or .hdf.")
        reports.append(SparseTableOutput(params.sparse))
    if params.output is not None:
        if params.output == "-":
            reports.append(GenbankOutput('/dev/stdout')) # TODO: not cross platform, will not work on windows.
        else:
            reports.append(GenbankOutput(params.output))

    databases=None
    if params.databases is not None:
        databases = set(params.databases)

    metrics = list()
    if params.ji != 0.0:
        metrics.append(JaccardIndex(params.ji, databases=databases))
    if params.ai != 0.0:
        metrics.append(AdjacencyIndex(params.ai, databases=databases))
    # if params.dss != 0.0:
    #     reports["dss"] = params.dss
    if len(metrics) == 0:
        raise ValueError("Please specify a non-zero weight for at least one metric.")


        
    ## figure out if we are using every contig or just certain, named contigs ##
    contigs_needed = list_and_file_to_dict_keys(params.contigs, params.contigs_file)
    
    # Run
    compare_contigs(genbanks, metrics, reports, params.k, contigs=contigs_needed, name_by_order=params.name_by_order)

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
