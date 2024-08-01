""" From a domainator-annotated genbank file, calculates contig stats and domain frequency stats.

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
import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal for <class 'numpy.float64'> type is zero.")
import argparse
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from io import BytesIO
import base64
import html
from domainator.utils import list_and_file_to_dict_keys, parse_seqfiles, TaxonomyData, get_taxid
from domainator.data_matrix import DataMatrix
from domainator import __version__, DOMAIN_FEATURE_NAME, RawAndDefaultsFormatter
import numpy as np
from typing import Dict, List
from bashplotlib.histogram import plot_hist
import io
from contextlib import redirect_stdout
from domainator.Taxonomy import NCBITaxonomy
import json
from collections import deque
from domainator.filter_domains import filter_domains


class SummaryTextWriter():
    def __init__(self, out_handle): #TODO: consolidate with HTMLWriter with a super class?
        # self.columns = columns
        # self.column_types = column_types
        self.out_handle = out_handle
    
    def write_header(self, num_contigs, num_cdss, lengths):
        if len(lengths) != 0:
            hist_str = io.StringIO()
            with redirect_stdout(hist_str):
                plot_hist(lengths, bincount=50, title="Contig Lengths", height=20.0, xlab=True, regular=True, showSummary=False)
            hist_str = hist_str.getvalue() # get the string
            hist_str = hist_str.replace('[39m', '') # remove color codes
        else:
            hist_str = ""
        print(f"""Contig Stats

contigs: {num_contigs}
CDSs: {num_cdss}
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

class SummaryHTMLWriter():
    def __init__(self, out_handle): #TODO: consolidate with TSVWriter with a super class?
        self.out_handle = out_handle
    
    def write_header(self, num_contigs, num_cdss, lengths):
        print(f"""<!doctype html><html lang="en"><head><meta charset="utf-8" name="viewport" content="width=device-width"/><title>Summary Report</title></head><body>
<div>
<h1>Contig Stats</h1>
<table>
<tbody>
    <tr> <th>contigs</th> <td>{num_contigs}</td> </tr> 
    <tr> <th>CDSs</th> <td>{num_cdss}</td> </tr> 
    <tr> <th>CDSs per 10 kb</th> <td>{num_cdss / ((sum(lengths)+0.01)/10000) : .2f}</td> </tr>
</tbody>
</table>
<h2>contig lengths</h2>
    <img alt="Contig lengths histogram" src="data:image/png;base64, {histplot_base64(lengths)}" />
</div>""", file=self.out_handle)
    
    def write_domain_table(self, domain_table_df, domain_cooccurence_df):
        print(f"""<div>
<h1>Domain Stats</h1>
<div>
    <h2>Domain Frequency</h2>
    <table border="1">
    <thead>
    {"<th>" + "</th><th>".join(domain_table_df.columns) + "</th>"}
    </thead>
    <tbody>
"""
        , file=self.out_handle)
        for index, row in domain_table_df.iterrows():
            print(f"""<tr>
        <td>{html.escape(row['domain'])}</td>
        <td>{html.escape(row['database'])}</td>
        <td align="right">{html.escape(str(row['count']))}</td>
        <td>{html.escape(row['description'])}</td>
        <td align="right">{html.escape(str(round(row['avg score']) ))}</td>
        </tr>
"""
            , file=self.out_handle)
        print(f"""</tbody>
            </table>
        </div>""", file=self.out_handle)
        if domain_cooccurence_df is not None and len(domain_cooccurence_df) > 0:
            print(f"""<div>
        <h2>Domain Co-occurence (what percent of the time when the row-domain occurs does the column domain occur on the same contig)</h2>
            {domain_cooccurence_df.to_html(float_format='{:.1f}'.format)}
        </div>
            """
            , file=self.out_handle)
        print("</div>", file=self.out_handle)

    def write_taxonomy(self, taxa_data):

        print("""<div id="taxonomy_chart">
<h2>Taxonomy</h2>
</div>
<script type=module>
// icicle plot code modified from code here: https://observablehq.com/@d3/zoomable-icicle
// Copyright 2018–2023 Observable, Inc.
// Permission to use, copy, modify, and/or distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.

import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";
function chart(data) {
// Specify the chart’s dimensions.
const width = 928;
const height = 1000;

// Create the color scale.
const color = d3.scaleOrdinal(d3.quantize(d3.interpolateRainbow, data.children.length + 1));

// Compute the layout.
const hierarchy = d3.hierarchy(data)
    .sum(d => d.value)
    .sort((a, b) => b.height - a.height || b.value - a.value);
const root = d3.partition()
    .size([height, (hierarchy.height + 1) * width / 3])
    (hierarchy);

// Create the SVG container.
const svg = d3.create("svg")
    .attr("viewBox", [0, 0, width, height])
    .attr("width", width)
    .attr("height", height)
    .attr("style", "max-width: 100%; height: auto; font: 10px sans-serif;");

// Append cells.
const cell = svg
    .selectAll("g")
    .data(root.descendants())
    .join("g")
    .attr("transform", d => `translate(${d.y0},${d.x0})`);

const rect = cell.append("rect")
    .attr("width", d => d.y1 - d.y0 - 1)
    .attr("height", d => rectHeight(d))
    .attr("fill-opacity", 0.6)
    .attr("fill", d => {
        if (!d.depth) return "#ccc";
        while (d.depth > 1) d = d.parent;
        return color(d.data.name);
    })
    .style("cursor", "pointer")
    .on("click", clicked);

const text = cell.append("text")
    .style("user-select", "none")
    .attr("pointer-events", "none")
    .attr("x", 4)
    .attr("y", 13)
    .attr("fill-opacity", d => +labelVisible(d));

text.append("tspan")
    .text(d => d.data.name);

const format = d3.format(",d");
const tspan = text.append("tspan")
    .attr("fill-opacity", d => labelVisible(d) * 0.7)
    .text(d => ` ${format(d.value)}`);

cell.append("title")
    .text(d => `${d.ancestors().map(d => d.data.name).reverse().join("/")}\n${format(d.value)}`);

// On click, change the focus and transitions it into view.
let focus = root;
function clicked(event, p) {
    focus = focus === p ? p = p.parent : p;

    root.each(d => d.target = {
    x0: (d.x0 - p.x0) / (p.x1 - p.x0) * height,
    x1: (d.x1 - p.x0) / (p.x1 - p.x0) * height,
    y0: d.y0 - p.y0,
    y1: d.y1 - p.y0
    });

    const t = cell.transition().duration(750)
        .attr("transform", d => `translate(${d.target.y0},${d.target.x0})`);

    rect.transition(t).attr("height", d => rectHeight(d.target));
    text.transition(t).attr("fill-opacity", d => +labelVisible(d.target));
    tspan.transition(t).attr("fill-opacity", d => labelVisible(d.target) * 0.7);
}

function rectHeight(d) {
    return d.x1 - d.x0 - Math.min(1, (d.x1 - d.x0) / 2);
}

function labelVisible(d) {
    return d.y1 <= width && d.y0 >= 0 && d.x1 - d.x0 > 16;
}

document.getElementById("taxonomy_chart").append(svg.node());
}""" , file=self.out_handle)
        
        
        print(f"const taxa_data = {json.dumps(tax_data_to_hierarchy(taxa_data))};", file=self.out_handle)
        #print(f"const taxa_data = {taxa_data};", file=self.out_handle)

        print(f"chart(taxa_data);\n</script>", file=self.out_handle)
                        
    def write_footer(self):
        print("</body></html>", file=self.out_handle)



#TODO: only print correlation table for focus domains. Otherwise it ends up being a huge, sparse table, and the output files are enormous.
def plt_to_base64_png():
    bytes = BytesIO()
    plt.savefig(bytes, format='png')
    bytes.seek(0)
    png_string = base64.b64encode(bytes.read()).decode("utf-8")
    return png_string

def histplot_base64(items_list, zero_center=False):
    binwidth=None
    if zero_center:
        bound = max(abs(min(items_list)),  abs(max(items_list)))
        if bound < 20:
            bound = 20
        binwidth = int(bound/10)
    if len(items_list) == 1:
        binwidth = None
    
    old_backend = plt.get_backend()
    plt.switch_backend('agg') # this is necessary to run on a headless server        
    sns.histplot(items_list, binwidth=binwidth)
    
    if zero_center:
        plt.xlim(-1*bound, bound)
    out = plt_to_base64_png()
    plt.switch_backend(old_backend)
    plt.close()
    return out

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

def tax_data_to_hierarchy(tax_data):
    """
        creates a data structure that can be used to create a d3 sunburst chart
        input: 
            tax_data: dict of dicts, taxid: {count: int, data: TaxonomyData}
        output:
            list of dicts, each dict has keys 'name', 'children', and 'size'
            where 'name' is the name of the taxon, 'children' is a list of dicts with the same structure, and 'size' is the number of contigs assigned to that taxon.
    """

    hierarchy = dict() # taxid: {name: str, children: set, size: int}

    hierarchy[1] = {'name': 'All', 'children': set(), 'size': 0}

    for taxid, count_and_data in tax_data.items():
        count = count_and_data['count']
        data = count_and_data['data']
        ranks_count = len(data.ranks)
        if taxid not in hierarchy:
            hierarchy[taxid] = {'name': data.names[-1], 'children': set(), 'size': count}
        else:
            hierarchy[taxid]['size'] += count
        prev_taxid = data.lineage[-1]
        for i in range(ranks_count-2, -1, -1): # iterate from the second to last rank to the first rank
            if data.ranks[i] != 'no rank':
                rank_taxid = data.lineage[i]
                if rank_taxid not in hierarchy:
                    hierarchy[rank_taxid] = {'name': data.names[i], 'rank': data.ranks[i], 'children': set(), 'size': 0}
                hierarchy[rank_taxid]['children'].add(prev_taxid)
                prev_taxid = rank_taxid
        hierarchy[1]['children'].add(prev_taxid)
    

    formatted_hierarchy = {"name": "All", "children": list(), "value": hierarchy[1]["size"]} # define the root node
    
    queue = list() #deque()
    for child in hierarchy[1]['children']:
        queue.append((child, formatted_hierarchy))
    #queue.append((1, formatted_hierarchy))

    #depth first traversal of the hierarchy
    while len(queue) > 0:

        taxid, parent_node = queue.pop()
        while len(hierarchy[taxid]['children']) == 1 and hierarchy[taxid]['size'] == 0 and hierarchy[taxid]['rank'] not in {"superkingdom", "kingdom"}: # skip nodes where (len(children) == 1 and size == 0), except for superkingdom and kingdom
            taxid = next(iter(hierarchy[taxid]['children']))
        
        node = {'name': hierarchy[taxid]['name'], 'value': hierarchy[taxid]['size']}
        if len(hierarchy[taxid]['children']) > 0:
            node['children'] = list()
            for child in hierarchy[taxid]['children']:
                queue.append((child, node))
        parent_node['children'].append(node)
    
    return formatted_hierarchy







    

def summary_report(records, out_text_handle, out_html_handle, focus_domains=None, co_occurence_fraction_outfile=None, co_occurence_count_outfile=None, report_taxonomy=False, ncbi_taxonomy=None, domains_table=None, databases=None):
    """

        Args:


        Returns:

    """
    longform_outputs = list()
    if out_text_handle is not None:
        longform_outputs.append(SummaryTextWriter(out_text_handle))
    if out_html_handle is not None:
        longform_outputs.append(SummaryHTMLWriter(out_html_handle))
    

    num_contigs = 0
    num_cdss = 0
    lengths = list() # contig lengths
    domains = dict() # dict of dicts, domain_name: {count: int, description: str, correlations: dict, average_score: float, database: str} 
    if focus_domains is None:
        focus_domains = list()

    domain_dists = dict() # dict of dicts of lists. focus_domain: other_domain: list_of_dists
    for focus_domain in focus_domains:
        domain_dists[focus_domain] = dict() # domain:list_of_dists

    taxa_data = dict() # dict of dicts, taxid: {count: int, data: TaxonomyData}

    for contig in records:
        if databases is not None:
            contig = tuple(filter_domains((contig,), evalue=float("inf"), max_overlap=1, databases_keep=databases))[0]
        
        lengths.append(len(contig))
        num_contigs += 1
        contig_domains = dict() #dict of lists of locations.
        
        if report_taxonomy:
            taxid = ncbi_taxonomy.normalized_taxid(get_taxid(contig))
            
            if taxid in taxa_data:
                taxa_data[taxid]['count'] += 1
            else:
                taxa_data[taxid] = {'count': 1, 'data': TaxonomyData(ncbi_taxonomy, taxid=taxid)}

        for feature in contig.features:
            if feature.type == 'CDS':
                num_cdss += 1
            if feature.type == DOMAIN_FEATURE_NAME:
                db = feature.qualifiers['database'][0]
                name = feature.qualifiers['name'][0]
                description = feature.qualifiers['description'][0]
                score = float(feature.qualifiers['score'][0])
                location = (feature.location.start, feature.location.end, feature.strand)
                id = name # TODO: someday distinguish different databases
                
                if id not in contig_domains:
                    contig_domains[id] = list()
                contig_domains[id].append(location)
                
                if id not in domains:
                    domains[id] = {'count': 0, 'description': set(), 'correlations':dict(), 'avg score':0, 'database': set()}
                
                domains[id]['database'].add(db) # TODO: maybe there is a better way to handle the same domain name occuring in different databases.
                domains[id]['description'].add(description)

                domains[id]['count'] += 1                
                
                # compute online avg score
                inv_count = 1/domains[id]['count']
                domains[id]['avg score'] = (inv_count * score) + ((1 - inv_count)*domains[id]['avg score'])
                
                # # domains[id]['locations'].add(location)
                # if description != domains[id]['description']:
                #     warnings.warn(f"Domain {id} has two different descriptions")
        
        # calculate correlations
        for domain1 in contig_domains:
            for domain2 in contig_domains:
                if domain2 not in domains[domain1]['correlations']:
                    domains[domain1]['correlations'][domain2] = 0
                domains[domain1]['correlations'][domain2] += 1
        
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
    for domain in domains:
        domains[domain]['database'] = '; '.join(domains[domain]['database'])
        domains[domain]['description'] = '; '.join(domains[domain]['description'])
    # Only calculate correlations for focus domains
    correlations_dict = {}
    for domain1 in focus_domains:
        if domain1 in domains:
            original_count =  domains[domain1]['correlations'][domain1]
            for domain2 in domains[domain1]['correlations']:
                domains[domain1]['correlations'][domain2] = 100 * (domains[domain1]['correlations'][domain2] / original_count)
            correlations_dict[domain1] = domains[domain1]["correlations"]

    domain_table_df = pd.DataFrame.from_dict(domains, orient='index', columns=['database', 'count', "description", 'avg score']).sort_values("count", ascending=False)
    domain_table_df.reset_index(inplace=True, names='domain')
    domain_cooccurence_df = None
    if len(correlations_dict) > 0:
        domain_cooccurence_df = pd.DataFrame.from_dict(correlations_dict, orient='index').fillna("")

    for output in longform_outputs:
        output.write_header(num_contigs, num_cdss, lengths)
        output.write_domain_table(domain_table_df, domain_cooccurence_df)
        if report_taxonomy:
            output.write_taxonomy(taxa_data)
        output.write_footer()

    if co_occurence_fraction_outfile is not None:
        write_cooccurence_matrix(domains, co_occurence_fraction_outfile, row_norm=True)
    
    if co_occurence_count_outfile is not None:
        write_cooccurence_matrix(domains, co_occurence_count_outfile, row_norm=False)

    if domains_table is not None:
        domain_table_df.to_csv(domains_table, sep="\t", index=False)

def main(argv):
    parser = argparse.ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", default=None, nargs='+', type=str, help="name of input genbank files.")
    
    parser.add_argument('-o', '--output', default=None, required=False,
                        help="Text file to write output to.")
    
    parser.add_argument('--html', default=None, required=False,
                        help="html file to write output to.")

    parser.add_argument('--domains_table', default=None, required=False,
                        help="Write a tab-separated table of domain frequencies to this file. Headers are: domain, database, count, description, avg score.")

    parser.add_argument('--co_occurence_fraction', default=None, required=False,
                        help="Write a dense text matrix of domain co-occurence frequencies to this file. Every row is a domain, the values are the fraction of contigs containing the row-domain that also include the column-domain (so the matrix is assymetric).")
    
    parser.add_argument('--co_occurence_count', default=None, required=False,
                        help="Write a dense text matrix of domain co-occurence frequencies to this file. Every row is a domain, the values are the number of contigs containing both the row-domain and the column-domain (so the matrix is symmetric).")

    parser.add_argument('--taxonomy', default=False, required=False, action="store_true",
                        help="Write a summary of the taxonomic origins of the contigs.")

    parser.add_argument("--taxonomy_update", action="store_true", help="If taxonomy database exists, check it against the version on the ncbi server and update if there is a newer version.")
    parser.add_argument("--ncbi_taxonomy_path", type=str,  default="/tmp/ncbi_taxonomy", help="Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist.")


    parser.add_argument('--contigs', type=str, default=None, nargs='+',
                        help="Only consider contigs with these names.")
    parser.add_argument('--contigs_file', type=str, default=None,
                        help="text file containing names of contigs to select.")
    
    # parser.add_argument('--domain_dist', action="store_true", default=False,
    #                     help="If set, then domain distance plots will be generated.")

    parser.add_argument('--domains', type=str, default=None, nargs='+',
                        help="calculate distance graphs for domains with these names")  # TODO: add some way to specify the database of a domain
    parser.add_argument('--domains_file', type=str, default=None,
                        help="text file containing names of domains to calculate distance graphs for.")

    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Ignore domain annotations not from these databases. default: consider all databases.")


    # TODO: domain scores histogram #for command line, maybe use bashplotlib
    # TODO: domain_search scores histogram
    params = parser.parse_args(argv)

    if params.output is None and params.html is None and params.co_occurence_fraction is None and params.co_occurence_count is None and params.domains_table is None:
        parser.error("No output file specified. Use --output, --html, --co_occurence_fraction, and/or --co_occurence_count and/or --domains_table to specify at least one output file.")


    out = None
    if params.output is not None:
        out = open(params.output, "w")
    
    out_html_handle = None
    if params.html is not None:
        out_html_handle = open(params.html, "w")
 

    # decide whether input is from stdin
    filetype_override = None
    if params.input is None:
        genbanks = [sys.stdin]
        filetype_override = "genbank"
    else:
        genbanks = params.input

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
    summary_report(
        parse_seqfiles(genbanks, target_contigs, filetype_override=filetype_override), 
        out,
        out_html_handle,
        selected_domains,
        params.co_occurence_fraction,
        params.co_occurence_count,
        report_taxonomy=params.taxonomy,
        ncbi_taxonomy=ncbi_taxonomy,
        domains_table=params.domains_table,
        databases=databases,
    )

    if params.output is not None and params.output != "-":
        out.close()
    
    if params.html is not None:
        out_html_handle.close()

def _entrypoint():
    main(sys.argv[1:])    

if __name__ == '__main__':
    main(sys.argv[1:])
  