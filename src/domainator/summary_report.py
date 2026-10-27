""" From a domainator-annotated genbank file, calculates contig stats and domain frequency stats

"""


### number of contigs
### length distribution
### cds_density
### raw domain counts
### domain co-occurence table
### (optional) domain distance density from "focus" CDSs (or focus domains)
### TODO: a ridge plot (joy plot) https://www.python-graph-gallery.com/ridgeline-graph-seaborn
###       https://www.python-graph-gallery.com/ridgeline-graph-plotly
### (optional) hierarchical tree of domain co-occurence
### (optional) score cutoffs

### TODO: add index
### TODO: maybe make the histograms selective, like skip if there is only one of the other domain.

### TODO: maybe make the output just a single html file instead of a directory?
from jsonargparse import ActionConfigFile
import contextlib
from pathlib import Path
import pandas as pd
import math
from domainator.cli.common import add_input, add_output, check_outputs, declare_outputs, tool_parser
from domainator.outputs import output_file
from domainator.cli.common import list_and_file_to_dict_keys, partial_help
from domainator.seqio import parse_seqfiles
from domainator.taxonomy.records import TaxonomyData, get_taxid
from domainator.partial import get_fragment_status, fragment_status_allowed, PARTIAL_CHOICES
from domainator.model import get_cds_unique_name, feature_score
from domainator.data_matrix import DataMatrix
from domainator import DOMAIN_FEATURE_NAME
import numpy as np
from typing import Dict, Iterable, Optional, Set, TextIO
from domainator.text_plot import text_histogram
from domainator.taxonomy import NCBITaxonomy, default_ncbi_taxonomy_path
import json
from domainator.filter_domains import filter_domains
from domainator import report_html
from Bio.SeqRecord import SeqRecord


class SummaryTextWriter():
    def __init__(self, out_handle):
        # self.columns = columns
        # self.column_types = column_types
        self.out_handle = out_handle
    
    def write_header(self, num_contigs, num_cdss, lengths, num_partial_cdss=0, num_fragment_proteins=0):
        if len(lengths) != 0:
            hist_str = text_histogram(lengths, title="Contig Lengths")
        else:
            hist_str = ""
        print(f"""Contig Stats

contigs: {num_contigs}
CDSs: {num_cdss}
partial CDSs: {num_partial_cdss}
protein fragments: {num_fragment_proteins}
CDSs per 10 kb: {num_cdss / ((sum(lengths)+0.01)/10000) : .2f}

{hist_str}""", file=self.out_handle)

    def write_domain_table(self, domain_table_df, domain_cooccurence_df):
        print(f"""Domain Stats
        
        Domain Frequency

{domain_table_df.to_string(justify="left", float_format='{:.0f}'.format, index=False)}
""", file=self.out_handle)
        if domain_cooccurence_df is not None and len(domain_cooccurence_df) > 0:
            print(f"""Domain Co-occurence (what percent of the time when the row-domain occurs does the column domain occur on the same contig)
{domain_cooccurence_df.to_string(justify="left", float_format='{:.1f}'.format)}
""", file=self.out_handle)

    def write_taxonomy(self, taxa_data):
        #taxa_data = dict() # dict of dicts, taxid: {count: int, data: TaxonomyData}
        table_data = {"taxid":[], "name":[], "rank":[], "count":[]}
        for taxid in taxa_data:
            table_data["taxid"].append(taxid)
            table_data["name"].append(taxa_data[taxid]["data"].names[-1])
            table_data["rank"].append(taxa_data[taxid]["data"].ranks[-1])
            table_data["count"].append(taxa_data[taxid]["count"])
            
        taxa_df = pd.DataFrame(table_data)

        print("Taxonomy", file=self.out_handle)
        print(taxa_df.to_string(justify="left", float_format='{:.0f}'.format), file=self.out_handle)

    def write_footer(self):
        pass

# The version of the --json_out document; see docs/file_formats.md.
JSON_SCHEMA_VERSION = 1


class SummaryJSONWriter():
    """Accumulates the summary into a single compact JSON object (agent-friendly).

    Mirrors the text/HTML writer interface; the full document is emitted on
    write_footer so the object is self-contained and easily parsed.
    """
    def __init__(self, out_handle):
        self.out_handle = out_handle
        self.payload = {"schema_version": JSON_SCHEMA_VERSION}

    def write_header(self, num_contigs, num_cdss, lengths, num_partial_cdss=0, num_fragment_proteins=0):
        length_stats = None
        if len(lengths) != 0:
            arr = np.asarray(lengths, dtype=float)
            length_stats = {
                "min": int(np.min(arr)),
                "max": int(np.max(arr)),
                "mean": float(np.mean(arr)),
                "median": float(np.percentile(arr, 50)),
                "total": int(np.sum(arr)),
            }
        self.payload["contig_stats"] = {
            "contigs": num_contigs,
            "cdss": num_cdss,
            "partial_cdss": num_partial_cdss,
            "fragment_proteins": num_fragment_proteins,
            "cds_per_10kb": num_cdss / ((sum(lengths) + 0.01) / 10000),
            "length": length_stats,
        }

    def write_domain_table(self, domain_table_df, domain_cooccurence_df):
        self.payload["domain_frequency"] = domain_table_df.rename(
            columns={"avg score": "avg_score"}
        ).to_dict(orient="records")
        if domain_cooccurence_df is not None and len(domain_cooccurence_df) > 0:
            self.payload["domain_cooccurrence"] = domain_cooccurence_df.to_dict(orient="index")
        else:
            self.payload["domain_cooccurrence"] = None

    def write_taxonomy(self, taxa_data):
        self.payload["taxonomy"] = [
            {
                "taxid": taxid,
                "name": taxa_data[taxid]["data"].names[-1],
                "rank": taxa_data[taxid]["data"].ranks[-1],
                "count": taxa_data[taxid]["count"],
            }
            for taxid in taxa_data
        ]

    def write_footer(self):
        json.dump(self.payload, self.out_handle, separators=(",", ":"), default=str)
        self.out_handle.write("\n")

class SummaryHTMLWriter(SummaryJSONWriter):
    """The interactive page (web/reports/summary_report.*): the --json_out document, and the
    contig length histogram and taxonomy tree the page draws."""
    def __init__(self, out_handle):
        super().__init__(out_handle)
        self.length_histogram = None
        self.taxonomy_tree = None

    def write_header(self, num_contigs, num_cdss, lengths, num_partial_cdss=0, num_fragment_proteins=0):
        super().write_header(num_contigs, num_cdss, lengths, num_partial_cdss=num_partial_cdss,
                             num_fragment_proteins=num_fragment_proteins)
        self.length_histogram = length_histogram(lengths)

    def write_taxonomy(self, taxa_data):
        super().write_taxonomy(taxa_data)
        self.taxonomy_tree = taxonomy_tree(taxa_data)

    def write_footer(self):
        # Read back, so the page holds exactly what --json_out writes.
        summary = json.loads(json.dumps(self.payload, default=str))
        data = {"summary": summary, "length_histogram": self.length_histogram, "taxonomy_tree": self.taxonomy_tree}
        self.out_handle.write(report_html.render_report("summary_report", "Summary Report", "summary_report", data))


# A histogram aims for about sqrt(n) bins, at least Sturges' ceil(log2(n)) + 1 and at most this.
MAX_HISTOGRAM_BINS = 50


def _nice_bin_width(raw: int) -> int:
    """The smallest 1, 2 or 5 times a power of ten that is at least ``raw`` (>= 1)."""
    power = 1
    while True:
        for factor in (1, 2, 5):
            if factor * power >= raw:
                return factor * power
        power *= 10


def length_histogram(lengths) -> Optional[Dict]:
    """Integer lengths binned for the page: ``{"bin_edges": [...], "counts": [...]}``, bin i
    counting the lengths in [bin_edges[i], bin_edges[i + 1]), or None for no lengths.

    Integer arithmetic only, so other implementations can give the same bins: n lengths from lo
    to hi aim for ``target = min(50, max((n - 1).bit_length() + 1, isqrt(n)))`` bins; the bin
    width is the smallest 1, 2 or 5 times a power of ten that is at least
    ``ceil((hi - lo + 1) / target)``, and the first edge is lo rounded down to a multiple of it.
    """
    if len(lengths) == 0:
        return None
    lo, hi, n = min(lengths), max(lengths), len(lengths)
    target = min(MAX_HISTOGRAM_BINS, max((n - 1).bit_length() + 1, math.isqrt(n)))
    width = _nice_bin_width(-(-(hi - lo + 1) // target))
    start = lo // width * width
    bins = (hi - start) // width + 1
    counts = [0] * bins
    for length in lengths:
        counts[(length - start) // width] += 1
    return {"bin_edges": [start + i * width for i in range(bins + 1)], "counts": counts}


# Ranks kept in the taxonomy tree even when they add nothing (one child, no contigs of their own).
TREE_KEPT_RANKS = {"superkingdom", "domain", "realm", "kingdom"}


def taxonomy_tree(taxa_data) -> Dict:
    """The taxa as a tree for the page: nodes ``{"taxid", "name", "rank", "count", "children"}``,
    ``count`` being the contigs assigned to that taxon itself; the root is taxid 1, "All".

    A taxon hangs from its nearest ancestor that has a rank (not "no rank"). A node with one
    child and no contigs of its own is replaced by that child, unless its rank is in
    TREE_KEPT_RANKS. Children are ordered by their total contigs, most first, then by name and
    taxid.

    Args:
        taxa_data: taxid: {"count": int, "data": TaxonomyData}.
    """
    root = {"taxid": 1, "name": "All", "rank": "no rank", "count": 0, "children": {}}
    for count_and_data in taxa_data.values():
        data = count_and_data["data"]
        last = len(data.lineage) - 1
        node = root
        for i, taxid in enumerate(data.lineage):
            if taxid == 1 or (i < last and data.ranks[i] == "no rank"):
                continue
            children = node["children"]
            if taxid not in children:
                children[taxid] = {"taxid": taxid, "name": data.names[i], "rank": data.ranks[i], "count": 0, "children": {}}
            node = children[taxid]
        node["count"] += count_and_data["count"]

    def finish(node):
        """(the node with its children as a sorted list, its total count)"""
        children = []
        for child in node["children"].values():
            while len(child["children"]) == 1 and child["count"] == 0 and child["rank"] not in TREE_KEPT_RANKS:
                child = next(iter(child["children"].values()))
            children.append(finish(child))
        children.sort(key=lambda pair: (-pair[1], pair[0]["name"], pair[0]["taxid"]))
        total = node["count"] + sum(total for _, total in children)
        return {**node, "children": [child for child, _ in children]}, total

    return finish(root)[0]


#TODO: only print correlation table for focus domains. Otherwise it ends up being a huge, sparse table, and the output files are enormous.
def write_cooccurence_matrix(domains:Dict[str, Dict], dense_text_file: str, row_norm: bool):
    """
        calculates and writes a co-occurence frequency matrix.

    Args:
        domains: domain_name: {count: int, description: str, locations} 
        dense_text_file (str): path to write the new file to.
        row_norm: if True, then divide each row by the diagonal before saving.
    """
    domain_names = list(domains.keys())
    domains_ct = len(domain_names)
    corr_mat = np.zeros((domains_ct,domains_ct))
    for domain_1_i, domain1 in enumerate(domain_names):
        for domain_2_i, domain2 in enumerate(domain_names):
            if domain2 in domains[domain1]['correlations']:
                corr_mat[domain_1_i, domain_2_i] = domains[domain1]['correlations'][domain2]

    if row_norm:
        corr_mat = corr_mat / corr_mat.diagonal()[:,np.newaxis]


    DataMatrix.write_dense_text(corr_mat, dense_text_file, domain_names, domain_names)

def domain_dist(domain1_coords, domain2_coords):
    """
        input:
            tuples of start, end, strand

        
        if domain1 and domain2 overlap, return 0
        otherwise return the distance between the closest edges
    """

    assert domain1_coords[0] <= domain1_coords[1]
    assert domain2_coords[0] <= domain2_coords[1]
    
    out = 0
    # domain2 has an edge inside domain1
    if (domain2_coords[0] >= domain1_coords[0] and domain2_coords[0] <= domain1_coords[1]) or (domain2_coords[1] >= domain1_coords[0] and domain2_coords[1] <= domain1_coords[1]):
        pass
    # domain1 has an edge inside domain2
    elif (domain1_coords[0] >= domain2_coords[0] and domain1_coords[0] <= domain2_coords[1]) or (domain1_coords[1] >= domain2_coords[0] and domain1_coords[1] <= domain2_coords[1]):
        pass
    # domain2 is before domain1 on the contig
    elif (domain2_coords[1] < domain1_coords[0]):
        out = -1 * (domain1_coords[0] - domain2_coords[1])
    # domain2 is after domain1 on the contig
    elif (domain2_coords[0] > domain1_coords[1]):
        out = domain2_coords[0] - domain1_coords[1]
    
    if domain1_coords[2] == -1 and out != 0: #TODO: test this
        out = -1 * out
    
    return out

def summary_report(records: Iterable[SeqRecord], *, text_out_handle: Optional[TextIO] = None,
                   html_out_handle: Optional[TextIO] = None,
                   domains: Optional[Set[str]] = None, co_occurrence_fraction_out: Optional[str] = None,
                   co_occurrence_count_out: Optional[str] = None, taxonomy: bool = False,
                   ncbi_taxonomy: Optional[NCBITaxonomy] = None, domains_tsv_out: Optional[str] = None,
                   databases: Optional[Set[str]] = None, json_out_handle: Optional[TextIO] = None,
                   partial: str = "include") -> None:
    """Write a summary of ``records``: contig and CDS counts, contig lengths and a table of the
    domains found, to any of the text, HTML and JSON handles (each can be None).

    Args:
        records: the annotated contigs or proteins (SeqRecords).
        text_out_handle: an open text handle for the text report.
        html_out_handle: an open text handle for the HTML report.
        json_out_handle: an open text handle for the JSON report.
        domains: focus domains: the report gives how often each other domain occurs with them,
            and how far apart they are.
        co_occurrence_fraction_out: a path for the domain-by-domain co-occurrence matrix, as
            fractions of contigs.
        co_occurrence_count_out: a path for the co-occurrence matrix as counts of contigs.
        taxonomy: also summarize the records' taxonomy, using ``ncbi_taxonomy`` (an NCBITaxonomy).
        domains_tsv_out: a path for the domain table as TSV.
        databases: count only domains from these databases.
        partial: which CDSs/proteins to count, by fragment status: "include" (all, the default),
            "exclude" (skip fragments), or "only" (only fragments). Protein records that are
            filtered out are skipped entirely. On nucleotide contigs, filtered CDSs and the
            domains on them are not counted, but the contig is. See partial.get_fragment_status.
    """
    longform_outputs = list()
    if text_out_handle is not None:
        longform_outputs.append(SummaryTextWriter(text_out_handle))
    if html_out_handle is not None:
        longform_outputs.append(SummaryHTMLWriter(html_out_handle))
    if json_out_handle is not None:
        longform_outputs.append(SummaryJSONWriter(json_out_handle))
    

    num_contigs = 0
    num_cdss = 0
    num_partial_cdss = 0 # CDSs with '<' or '>' in their location
    num_fragment_proteins = 0 # protein records that are fragments (see utils.get_fragment_status)
    lengths = list() # contig lengths
    domain_stats = dict() # dict of dicts, domain_name: {count: int, description: str, correlations: dict, average_score: float, database: str} 
    if domains is None:
        domains = list()

    domain_dists = dict() # dict of dicts of lists. focus_domain: other_domain: list_of_dists
    for focus_domain in domains:
        domain_dists[focus_domain] = dict() # domain:list_of_dists

    taxa_data = dict() # dict of dicts, taxid: {count: int, data: TaxonomyData}

    for contig in records:
        if databases is not None:
            contig = tuple(filter_domains((contig,), keep_databases=databases))[0]

        excluded_cds_ids = set() # cds_ids of CDSs filtered out by partial
        if partial != "include":
            if contig.annotations.get("molecule_type") == "protein":
                if not fragment_status_allowed(get_fragment_status(contig), partial):
                    continue
            else:
                excluded_cds_ids = {get_cds_unique_name(f) for f in contig.features if f.type == "CDS" and not fragment_status_allowed(get_fragment_status(f), partial)}
        
        lengths.append(len(contig))
        num_contigs += 1
        if contig.annotations.get("molecule_type") == "protein" and get_fragment_status(contig) is not None:
            num_fragment_proteins += 1
        contig_domains = dict() #dict of lists of locations.
        
        if taxonomy:
            taxid = ncbi_taxonomy.normalized_taxid(get_taxid(contig))
            
            if taxid in taxa_data:
                taxa_data[taxid]['count'] += 1
            else:
                taxa_data[taxid] = {'count': 1, 'data': TaxonomyData(ncbi_taxonomy, taxid=taxid)}

        for feature in contig.features:
            if feature.type == 'CDS' and (not excluded_cds_ids or get_cds_unique_name(feature) not in excluded_cds_ids):
                num_cdss += 1
                if get_fragment_status(feature) is not None:
                    num_partial_cdss += 1
            if feature.type == DOMAIN_FEATURE_NAME and (not excluded_cds_ids or feature.qualifiers.get("cds_id", ["."])[0] not in excluded_cds_ids):
                db = feature.qualifiers['database'][0]
                name = feature.qualifiers['name'][0]
                description = feature.qualifiers['description'][0]
                score = feature_score(feature)
                location = (feature.location.start, feature.location.end, feature.location.strand)
                id = name # TODO: someday distinguish different databases
                
                if id not in contig_domains:
                    contig_domains[id] = list()
                contig_domains[id].append(location)
                
                if id not in domain_stats:
                    domain_stats[id] = {'count': 0, 'description': set(), 'correlations':dict(), 'avg score':0, 'database': set()}
                
                domain_stats[id]['database'].add(db) # TODO: maybe there is a better way to handle the same domain name occuring in different databases.
                domain_stats[id]['description'].add(description)

                domain_stats[id]['count'] += 1                
                
                # compute online avg score
                inv_count = 1/domain_stats[id]['count']
                domain_stats[id]['avg score'] = (inv_count * score) + ((1 - inv_count)*domain_stats[id]['avg score'])
        
        # calculate correlations
        for domain1 in contig_domains:
            for domain2 in contig_domains:
                if domain2 not in domain_stats[domain1]['correlations']:
                    domain_stats[domain1]['correlations'][domain2] = 0
                domain_stats[domain1]['correlations'][domain2] += 1
        
        # calculate domain distances
        for focus_domain in domain_dists:
            if focus_domain in contig_domains:
                for focus_location in contig_domains[focus_domain]:
                    for other_domain in contig_domains:
                        for other_location in contig_domains[other_domain]:
                            if (focus_domain, focus_location) != (other_domain, other_location):
                                if other_domain not in domain_dists[focus_domain]:
                                   domain_dists[focus_domain][other_domain] = list()
                                domain_dists[focus_domain][other_domain].append(domain_dist(focus_location, other_location))
    
    # convert database set to string
    for domain in domain_stats:
        domain_stats[domain]['database'] = '; '.join(domain_stats[domain]['database'])
        domain_stats[domain]['description'] = '; '.join(domain_stats[domain]['description'])
    # Only calculate correlations for focus domains
    correlations_dict = {}
    for domain1 in domains:
        if domain1 in domain_stats:
            # A new dict: the raw counts are still needed for the co-occurrence matrix outputs.
            original_count = domain_stats[domain1]['correlations'][domain1]
            correlations_dict[domain1] = {domain2: 100 * (count / original_count)
                                          for domain2, count in domain_stats[domain1]['correlations'].items()}

    domain_table_df = pd.DataFrame.from_dict(domain_stats, orient='index', columns=['database', 'count', "description", 'avg score']).sort_values("count", ascending=False)
    domain_table_df.reset_index(inplace=True, names='domain')
    domain_cooccurence_df = None
    if len(correlations_dict) > 0:
        domain_cooccurence_df = pd.DataFrame.from_dict(correlations_dict, orient='index').fillna("")

    for output in longform_outputs:
        output.write_header(num_contigs, num_cdss, lengths, num_partial_cdss=num_partial_cdss, num_fragment_proteins=num_fragment_proteins)
        output.write_domain_table(domain_table_df, domain_cooccurence_df)
        if taxonomy:
            output.write_taxonomy(taxa_data)
        output.write_footer()

    if co_occurrence_fraction_out is not None:
        write_cooccurence_matrix(domain_stats, co_occurrence_fraction_out, row_norm=True)
    
    if co_occurrence_count_out is not None:
        write_cooccurence_matrix(domain_stats, co_occurrence_count_out, row_norm=False)

    if domains_tsv_out is not None:
        domain_table_df.to_csv(domains_tsv_out, sep="\t", index=False)

def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, "name of input genbank files.")

    add_output(parser, "text", "Write the report as text to this file.", primary=True)
    add_output(parser, "html", "Write the report as html to this file.")
    add_output(parser, "json", "Write the full summary as a single compact JSON object to this file.")
    add_output(parser, "domains_tsv", "Write a tab-separated table of domain frequencies to this file. Headers are: domain, database, count, description, avg score.", stdout=False)
    add_output(parser, "co_occurrence_fraction", "Write a dense text matrix of domain co-occurrence frequencies to this file. Every row is a domain, the values are the fraction of contigs containing the row-domain that also include the column-domain (so the matrix is asymmetric).", stdout=False)
    add_output(parser, "co_occurrence_count", "Write a dense text matrix of domain co-occurrence frequencies to this file. Every row is a domain, the values are the number of contigs containing both the row-domain and the column-domain (so the matrix is symmetric).", stdout=False)

    parser.add_argument('--taxonomy', default=False, required=False, action="store_true",
                        help="Write a summary of the taxonomic origins of the contigs.")

    parser.add_argument("--taxonomy_update", action="store_true", help="If taxonomy database exists, check it against the version on the ncbi server and update if there is a newer version.")
    parser.add_argument("--ncbi_taxonomy_path", type=str,  default=default_ncbi_taxonomy_path(), help="Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist.")


    parser.add_argument('--contigs', type=str, default=None, nargs='+',
                        help="Only consider contigs with these names.")
    parser.add_argument('--contigs_file', type=str, default=None,
                        help="text file containing names of contigs to select.")
    
    parser.add_argument('--domains', type=str, default=None, nargs='+',
                        help="calculate distance graphs for domains with these names")  # TODO: add some way to specify the database of a domain
    parser.add_argument('--domains_file', type=str, default=None,
                        help="text file containing names of domains to calculate distance graphs for.")

    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Ignore domain annotations not from these databases. default: consider all databases.")
    parser.add_argument('--partial', type=str, default="include", choices=PARTIAL_CHOICES,
                        help=partial_help("Which CDSs/proteins to count") + " Protein records that are filtered out are skipped; on nucleotide contigs, filtered CDSs and their domains are not counted, but the contig is.")

    parser.add_argument('--config', action=ActionConfigFile)

    # TODO: domain_search scores histogram
    declare_outputs(parser, ["text_out", "html_out", "json_out", "domains_tsv_out",
                                   "co_occurrence_fraction_out", "co_occurrence_count_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)

    check_outputs(parser, params)

    target_contigs = list_and_file_to_dict_keys(
        params.contigs, params.contigs_file, as_set=True)
    selected_domains = list_and_file_to_dict_keys(
        params.domains, params.domains_file, as_set=True)
    
    ncbi_taxonomy = None
    if params.taxonomy:
        # create the path to the NCBI taxonomy database
        Path(params.ncbi_taxonomy_path).mkdir(parents=True, exist_ok=True)
        # load the NCBI taxonomy database
        ncbi_taxonomy = NCBITaxonomy(params.ncbi_taxonomy_path, params.taxonomy_update)

    databases=None
    if params.databases is not None:
        databases = set(params.databases)

    ### Run
    with contextlib.ExitStack() as outputs:
        def open_output(path):
            return None if path is None else outputs.enter_context(output_file(path))
        summary_report(
            parse_seqfiles(params.input, target_contigs),
            text_out_handle=open_output(params.text_out),
            html_out_handle=open_output(params.html_out),
            domains=selected_domains,
            co_occurrence_fraction_out=params.co_occurrence_fraction_out,
            co_occurrence_count_out=params.co_occurrence_count_out,
            taxonomy=params.taxonomy,
            ncbi_taxonomy=ncbi_taxonomy,
            domains_tsv_out=params.domains_tsv_out,
            databases=databases,
            json_out_handle=open_output(params.json_out),
            partial=params.partial,
        )
