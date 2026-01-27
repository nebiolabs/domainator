"""Writes record-wise information from Genbank files into tab-separated or html tables

Takes a genbank file which has been annotated using Domainate, write tabulated descriptions of each contig, CDS, or domain.

"""

import sys
import argparse
from jsonargparse import ArgumentParser, ActionConfigFile
from domainator.utils import get_sources, DomainatorCDS, parse_seqfiles, list_and_file_to_dict_keys, slice_record_from_location, TaxonomyData
from domainator.select_by_cds import get_cds_neighborhood
from domainator import __version__, DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME, RawAndDefaultsFormatter
from pathlib import Path
from typing import List, Tuple, Union, Dict, Any
import json
from domainator.Bio.SeqFeature import FeatureLocation, CompoundLocation
from domainator.Bio.SeqRecord import SeqRecord
from domainator.Taxonomy import NCBITaxonomy
from domainator.filter_domains import filter_domains

SELECT_CATEGORIES = {"contig", "cds", "domain"}
COLS_ARG_NAME = "output_cols"

class DynamicArg(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)
    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest, None)
        if items is None:
            items = []
        items.append( (self.const, values) )
        setattr(namespace, self.dest, items)

def get_strand(loc: Union[FeatureLocation, CompoundLocation]) -> str:
    if loc.parts[0].strand == -1:
        return "-"
    else:
        return "+"

def get_domain_content(rec:SeqRecord, sep:str="; ") -> Tuple[str, str]:
    """
        given a SeqRecord, return the domain content as a semicolon-separated string, in alphabetical order.
    """
    domains = set()

    for feature in rec.features:
        if feature.type == DOMAIN_FEATURE_NAME:
            domains.add((feature.qualifiers["name"][0], feature.qualifiers["description"][0]))

    domains = list(domains)
    descriptions = tuple()
    domains.sort(key=lambda x: x[0])
    if len(domains) > 0:
        domains, descriptions = zip(*domains)
    return sep.join(domains), sep.join(descriptions)

def get_domain_architecture(rec:SeqRecord, sep:str="; ") -> str:
    domainator_features = list()
    for feature in rec.features:
        if feature.type == DOMAIN_FEATURE_NAME:
            domainator_features.append(feature)
    domainator_features.sort(key=lambda x:x.location.start)
    return sep.join([feature.qualifiers["name"][0] for feature in domainator_features])

def get_domain_architecture_detailed(rec:SeqRecord, domain_sep:str="/", CDS_sep:str="; ", no_domains_tag="no_domains") -> str: #TODO: add test for this.
    CDS_features = list()
    
    if rec.annotations["molecule_type"] == "protein":
        domainator_features = list()
        for feature in rec.features:
            if feature.type == DOMAIN_FEATURE_NAME:
                domainator_features.append(feature)
        if len(domainator_features) == 0:
            CDS_features.append(no_domains_tag)
        else:
            domainator_features.sort(key=lambda x:x.location.start)
            CDS_features.append(domain_sep.join([feature.qualifiers["name"][0] for feature in domainator_features]))
        
    else: # nucleic acid
        cdss = DomainatorCDS.list_from_contig(rec, skip_pseudo=True)
        cdss.sort(key=lambda x: x.feature.location.stranded_start)
        for cds in cdss:
            domainator_features = list()
            for feature in cds.domain_features:
                domainator_features.append(feature)
            
            if len(domainator_features) == 0:
                CDS_features.append(no_domains_tag)
            else:
                domainator_features.sort(key=lambda x:x.location.stranded_start)
                if (cds.feature.location.strand == -1):
                    domainator_features.reverse()
                CDS_features.append(domain_sep.join([feature.qualifiers["name"][0] for feature in domainator_features]))

    return CDS_sep.join(CDS_features)


def get_domain_scores_sum(rec:SeqRecord)->float:
    score = 0.0
    for feature in rec.features:
        if feature.type == DOMAIN_FEATURE_NAME:
            score += float(feature.qualifiers["score"][0])
    return score

def get_domain_identities_sum(rec:SeqRecord)->float:
    score = 0.0
    for feature in rec.features:
        if feature.type == DOMAIN_FEATURE_NAME:
            score += float(feature.qualifiers.get("identity",(0,))[0])
    return score

def count_feature(rec:SeqRecord, feature_type:str)->int:
    count = 0
    for feature in rec.features:
        if feature.type == feature_type:
            if feature.type == 'CDS' and ("pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers):
                continue # skip pseudogenes for CDSs
            count += 1
    return count

def count_qualifier(rec:SeqRecord, feature_type:str, qualifier:str)->int:
    count = 0
    for feature in rec.features:
        if feature.type == feature_type:
            if qualifier in feature.qualifiers:
                count += 1
    return count

def domain_search(rec:SeqRecord) -> Tuple[str, str, int, int, int, int, int, str, float]:
    """

    Args:
        rec (SeqRecord): 

    Returns:
        Tuple(search_best_query, best_search_query_description, search_match_start, search_match_end, search_match_score)
    """
    search_best_query = ""
    best_search_query_description = ""
    query_start = 0
    query_end = 0
    query_len = 0
    search_match_start = 0
    search_match_end = 0
    search_match_strand = ""
    search_match_score = 0
    search_match_identity = 0
    
    for feature in rec.features:
        if feature.type == DOMAIN_SEARCH_BEST_HIT_NAME:
            search_best_query = str(feature.qualifiers["name"][0])
            best_search_query_description = str(feature.qualifiers["description"][0])
            query_start = str(feature.qualifiers["rstart"][0])
            query_end = str(feature.qualifiers["rend"][0])
            query_len = str(feature.qualifiers["rlen"][0])
            search_match_start = feature.location.stranded_start_human_readable
            search_match_end = feature.location.stranded_end_human_readable
            search_match_strand = get_strand(feature.location)
            search_match_score = float(feature.qualifiers["score"][0])
            search_match_identity = float(feature.qualifiers.get("identity",(0,))[0])
            break
    
    return (search_best_query, best_search_query_description, query_start, query_end, query_len, search_match_start, search_match_end, search_match_strand, search_match_score, search_match_identity)

class EnumTSVWriter():
    def __init__(self, columns, column_types, out_handle): #TODO: consolidate with HTMLWriter with a super class?
        self.columns = columns
        self.column_types = column_types
        self.out_handle = out_handle
    
    def write_header(self):
        print("\t".join(self.columns), file=self.out_handle)

    def write_row(self, values):
        out_vals = list()
        try:
            for i, v in enumerate(values):
                data_type = self.column_types[i]
                if data_type == "str":
                    if v is None:
                        v = ""
                    v = v.replace("\t", " ")
                    out_vals.append(v)
                if data_type == "int":
                    if v is None:
                        v = ""
                    else:
                        v = str(int(v))
                    out_vals.append(v)
                if data_type == "float":
                    out_vals.append(f"%0.1f" % float(v))
        except TypeError:
            print(values, file=sys.stderr)
            print(self.columns, file=sys.stderr)
            print(self.columns[i], file=sys.stderr)
            print(self.column_types, file=sys.stderr)
            raise
        print("\t".join(out_vals), file=self.out_handle)

    def write_footer(self):
        pass

class EnumHTMLWriter():
    sorters = {"str": "string", "int": "number", "float": "number"}
    def __init__(self, columns, column_types, out_handle, max_height): #TODO: consolidate with TSVWriter with a super class?
        self.columns = columns
        self.column_types = column_types
        self.out_handle = out_handle
        self.max_height = max_height
        self.id_counter = 0
    
    def write_header(self):
        print("""<!doctype html><html lang='en'>
<head>
<meta charset="UTF-8" />
<title>Enum Report</title>
<link href="https://unpkg.com/tabulator-tables@5.4.3/dist/css/tabulator.min.css" rel="stylesheet">
<script type="text/javascript" src="https://unpkg.com/tabulator-tables@5.4.3/dist/js/tabulator.min.js"></script>
</head>

<body>
<div>
    <button id="download-csv">Download CSV</button>
    <button id="download-html">Download HTML</button>
    <button id="download-json">Download JSON</button>
</div>
<div id="enum_report_table"></div>
<script>
var tabledata = [""",
        file=self.out_handle)
        self.id_counter = 0
        

    def write_row(self, values):
        out_vals = list()
        for i, v in enumerate(values):
            data_type = self.column_types[i]
            if data_type == "str":
                out_vals.append(json.dumps(v))
            if data_type == "int":
                out_vals.append(str(int(v)))
            if data_type == "float":
                out_vals.append(f"%0.1f" % float(v))
        self.out_handle.write("{id:" + str(self.id_counter) + ",")
        for i,c in enumerate(out_vals):
            self.out_handle.write(f"""{json.dumps(self.columns[i])}:{c}, """)
        self.out_handle.write("},\n")
        self.id_counter += 1

    def write_footer(self):
        print("""];
    var table = new Tabulator("#enum_report_table", {
    data:tabledata,
    pagination:"local",
    paginationSize:25,
    paginationSizeSelector:[25, 50, 100, 1000],
    movableColumns:true,
    paginationCounter:"rows",""", file=self.out_handle)
        if self.max_height is not None:
            print(f"    maxHeight:{self.max_height},", file=self.out_handle)
        print("    columns:[", file=self.out_handle)
        for i, column in enumerate(self.columns):
            print(json.dumps({"title":column, "field":column, "width":150, "sorter": self.sorters[self.column_types[i]], "tooltip":True}) + ",", file=self.out_handle)
        print("    ],", file=self.out_handle)
        print("""});

//trigger download of data.csv file
document.getElementById("download-csv").addEventListener("click", function(){
    table.download("csv", "data.csv");
});

//trigger download of data.json file
document.getElementById("download-json").addEventListener("click", function(){
    table.download("json", "data.json");
});

//trigger download of data.html file
document.getElementById("download-html").addEventListener("click", function(){
    table.download("html", "data.html", {style:true});
});

</script>
</body></html>""", file=self.out_handle)

def taxid_factory(ranks):
    
    def get_taxids(rec, loc, tax):
        taxids = list()
        for r in ranks:
            r = r.lower()
            if r == "self":
                taxids.append(tax.lineage[-1])
            elif r in tax.rank_to_taxid:
                taxids.append(tax.rank_to_taxid[r])
            else:
                taxids.append(None)
        return taxids

    return {
        "columns": ["taxid_" + r for r in ranks],
        "column_types": ["int"] * len(ranks),
        "function": get_taxids,
    }

def taxname_factory(ranks):
    
    def get_taxnames(rec, loc, tax):
        taxnames = list()
        for r in ranks:
            r = r.lower()
            if r == "self":
                taxnames.append(tax.names[-1])
            elif r in tax.rank_to_name:
                taxnames.append(tax.rank_to_name[r])
            else:
                taxnames.append(None)
        return taxnames

    return {
        "columns": ["taxname_" + r for r in ranks],
        "column_types": ["str"] * len(ranks),
        "function": get_taxnames,
    }

def append_factory(spec):
    
    if len(spec) != 3:
        raise ValueError("append spec must be length 3")
    if spec[1] not in ["int", "float", "str"]:
        raise ValueError("append spec must be one of int, float, str")

    def append(rec, loc, tax):
        return spec[2]

    return {
        "columns": [spec[0]],
        "column_types": [spec[1]],
        "function": append,
    }

def qualifier_factory(spec, sep="; "):
    
    if len(spec) != 2:
        raise ValueError("qualifier spec must be length 2")

    def qualifier(rec, loc, tax):
        rec.features.sort(key=lambda x: x.location.stranded_start)
        hits = list()
        for feature in rec.features:
            if feature.type == spec[0]:
                if feature.type == 'CDS' and ("pseudo" in feature.qualifiers and "pseudogene" in feature.qualifiers):
                    continue # skip pseudogenes
                if spec[1] in feature.qualifiers:
                    hits.append(feature.qualifiers[spec[1]][0])
        return sep.join(hits)

    return {
        "columns": ["qualifier_" + spec[0] + "_" + spec[1]],
        "column_types": ["str"],
        "function": qualifier,
    }

def feature_count_factory(spec): #TODO: deduplicate with count_feature
    
    if len(spec) != 1:
        raise ValueError("feature_count spec must be length 1")

    def qualifier(rec, loc, tax):
        count = 0
        for feature in rec.features:
            if feature.type == spec[0]:
                if feature.type == 'CDS' and ("pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers):
                    continue # skip pseudogenes
                count += 1
        return count

    return {
        "columns": [spec[0] + "_count"],
        "column_types": ["int"],
        "function": qualifier,
    }



def get_analysis_names(analyses):
    analysis_names = list()
    for analysis in analyses:
        if isinstance(analysis, str):
            analysis_names.append(analysis)
        else:
            analysis_names.append(analysis[0])
    return set(analysis_names)

def process_record(rec:SeqRecord, by:str, analyses_to_run:List[Dict[str,Any]], ncbi_taxonomy):
    """
    Args:
        rec (SeqRecord): 
        by (str): one of "contig", "cds", "domain"
        analyses_to_run (list): list of dicts with keys "columns", "column_types", "function"
        ncbi_taxonomy (NCBITaxonomy): 

    Yields:
        list: one line of output
    """
    contig_name = rec.id
    # Preserve the source filename if it exists
    source_filename = rec.annotations.get("_source_filename", "")
    sub_recs = list()
    locs = list()
    if by == "contig":
        sub_recs.append(rec)
        locs.append(FeatureLocation(0, len(rec), 1))
    else: #by cds or domain
        cds_names = list()
        domain_names = list()
        sources = get_sources(rec)
        if rec.annotations['molecule_type'] != "protein": # nucleic acid
            cdss = DomainatorCDS.list_from_contig(rec)
            cdss.sort(key=lambda x: x.feature.location.stranded_start)

            for i in range(len(cdss)):
                if by == "cds":
                    cds_rec = slice_record_from_location(rec, cdss[i].feature.location, sources + [cdss[i].feature] + cdss[i].domain_features, truncate_features=True) #  get_cds_neighborhood(rec, cdss, i) # extract a single CDS.
                    if source_filename:
                        cds_rec.annotations["_source_filename"] = source_filename
                    sub_recs.append(cds_rec)
                    locs.append(cdss[i].feature.location)
                    cds_names.append(cdss[i].name)
                else: # by == "domain"
                    cds_summary = cdss[i]
                    for domain in cds_summary.domain_features:
                        domain_name = domain.qualifiers["name"][0]
                        domain_rec = slice_record_from_location(rec, domain.location, sources + [domain], truncate_features=True)
                        if source_filename:
                            domain_rec.annotations["_source_filename"] = source_filename
                        cds_names.append(cds_summary.name)
                        domain_names.append(domain_name)
                        sub_recs.append(domain_rec)
                        locs.append(domain.location)
        else: # protein
            if by == "cds":
                sub_recs.append(rec)
                cds_names.append(" ")
                locs.append(FeatureLocation(0, len(rec), 1))
            else: # by == "domain"
                sources = get_sources(rec)
                for feature in rec.features:
                    if feature.type == DOMAIN_FEATURE_NAME:
                        domain_rec = slice_record_from_location(rec, feature.location, sources + [feature], truncate_features=True)
                        if source_filename:
                            domain_rec.annotations["_source_filename"] = source_filename
                        cds_names.append(" ")        
                        domain_names.append(feature.qualifiers["name"][0])
                        sub_recs.append(domain_rec)
                        locs.append(feature.location)
            
    
    for i, sub_rec in enumerate(sub_recs):
        out_line = [contig_name]
        if by == "cds" or by == "domain":
            out_line.append(cds_names[i])
            if by == "domain":
                out_line.append(domain_names[i]) # I'm not sure what to call the domain, maybe just what domain it is?
        tax_data = None
        if ncbi_taxonomy is not None:
            tax_data = TaxonomyData(ncbi_taxonomy, record=sub_rec)

        for analysis in analyses_to_run:
            analysis_result = analysis["function"](sub_rec, locs[i], tax_data)
            
            if isinstance(analysis_result, list) or isinstance(analysis_result, tuple): # doesn't handle generators
                out_line.extend(analysis_result)
            else:
                out_line.append(analysis_result)
        yield out_line

def parse_seqfiles_with_filenames(seqfiles, contigs=None, filetype_override=None):
    """
    Wrapper around parse_seqfiles that adds the source filename to each record's annotations.
    
    Args:
        seqfiles: a list of genbank or fasta paths
        contigs: a dict where keys are the names of desired contigs, and values are None. Or None.
        filetype_override: if supplied then don't try to guess filetype from the extension
        
    Yields:
        SeqRecord objects with '_source_filename' added to annotations
    """
    from domainator.utils import open_if_is_name
    
    for file in seqfiles:
        # Determine the filename to store
        filename = None
        try:
            # If it's a file path or Path object, convert to string
            filename = str(file)
        except:
            # If it's already a file handle, try to get the name
            if hasattr(file, 'name'):
                filename = file.name
            else:
                filename = "<stdin>"
        
        # Parse records from this file
        for rec in parse_seqfiles([file], contigs, filetype_override):
            rec.annotations["_source_filename"] = filename
            yield rec

def enum_report(records, by, analyses, tsv_out_handle, html_out_handle, column_names, html_max_height, ncbi_taxonomy, databases=None):
    """
      input: 
        records: an iterator of SeqRecords
        by: One line in output for every by contig, cds or domain. default: by contig
        analyses: what data to report
        tsv_out_handle: a file handle to write the analysis table to.
        html_out_handle: a file handle to write html formatted output to.
        column_names: If supplied, then this list will be used instead of the default column names.
        html_max_height: Max height in pixels to set the html output to.
    """
    
    # {command_line_variable: {"columns":[names_to_appear_in_output], "function": function_mapping_SeqRec_to_table_value} }
    # the function should return a string if "columns" has one member, or a list of strings if "columns" has multiple members.
    STATIC_ANALYSES = {"domains": {"columns": ["domains"], "column_types": ["str"], "function": lambda rec,loc,tax: get_domain_content(rec)[0]}, 
                "domain_descriptions": {"columns": ["domain_descriptions"], "column_types": ["str"], "function": lambda rec,loc,tax: get_domain_content(rec)[1]}, 
                "sequence": {"columns": ["sequence"], "column_types": ["str"], "function": lambda rec,loc,tax: str(rec.seq)},  #TODO: maybe for domains and CDSs convert to forward strand.
                "start": {"columns": ["start"], "column_types": ["int"], "function": lambda rec,loc,tax: loc.stranded_start_human_readable}, #TODO: can be off by one depending on strand
                "end": {"columns": ["end"], "column_types": ["int"], "function": lambda rec,loc,tax: loc.stranded_end_human_readable}, #TODO: can be off by one depending on strand
                "strand": {"columns": ["strand"], "column_types": ["str"], "function": lambda rec,loc,tax: get_strand(loc)},
                "length": {"columns": ["length"], "column_types": ["int"], "function": lambda rec,loc,tax: str(len(rec))}, 
                "topology": {"columns": ["topology"], "column_types": ["str"], "function": lambda rec,loc,tax: rec.annotations.get("topology", "linear")}, 
                "architecture": {"columns": ["architecture"], "column_types": ["str"], "function": lambda rec,loc,tax: get_domain_architecture(rec)},
                "architecture_detailed": {"columns": ["architecture_detailed"], "column_types": ["str"], "function": lambda rec,loc,tax: get_domain_architecture_detailed(rec)},
                "score": {"columns": ["score"], "column_types": ["float"], "function":  lambda rec,loc,tax: get_domain_scores_sum(rec)}, 
                "identity": {"columns": ["identity"], "column_types": ["float"], "function":  lambda rec,loc,tax: get_domain_identities_sum(rec)}, 
                "definition": {"columns": ["definition"], "column_types": ["str"], "function": lambda rec,loc,tax: rec.description}, 
                "filename": {"columns": ["filename"], "column_types": ["str"], "function": lambda rec,loc,tax: rec.annotations.get("_source_filename", "")}, 
                "cds_count": {"columns": ["cds_count"], "column_types": ["int"], "function": lambda rec,loc,tax: count_feature(rec, "CDS")}, 
                "domain_count": {"columns": ["domain_count"], "column_types": ["int"], "function": lambda rec,loc,tax: count_feature(rec, DOMAIN_FEATURE_NAME)}, 
                "domain_search": {"columns": ["best_query", "query_description", "query_start", "query_end", "query_length", "match_start", "match_end", "match_strand", "match_score", "match_identity"], 
                                "column_types": ["str", "str", "int", "int", "int", "int", "int", "str", "float", "float"], "function": lambda rec,loc,tax: domain_search(rec)},
                "taxid_lineage": { "columns": ["taxid_lineage"], "column_types": ["str"], "function": lambda rec,loc,tax: "; ".join([str(x) for x in tax.lineage][1:]) }, # [1:] to remove root
                "taxname_lineage": { "columns": ["taxid_lineage"], "column_types": ["str"], "function": lambda rec,loc,tax: "; ".join([x for x in tax.names][1:]) }, # [1:] to remove root
                "rank_lineage": { "columns": ["rank_lineage"], "column_types": ["str"], "function": lambda rec,loc,tax: "; ".join([x for x in tax.ranks][1:]) }, # [1:] to remove root
                }
    DYNAMIC_ANALYSES = {"taxid": taxid_factory, "taxname": taxname_factory, "qualifier":qualifier_factory, "feature_count": feature_count_factory, "append": append_factory} # values are functions that return dicts of {"columns":[names_to_appear_in_output],  "column_types": [types_of_columns], "function": function taking rec, loc, tax as arguments and returning a scalar or list of scalars}
    
    headers = ["contig"]
    column_types = ["str"]
    analyses_to_run = list()
    if by in {"cds","domain"}:
        headers.append("cds")
        column_types.append("str")
    if by == "domain":
        headers.append("domain")
        column_types.append("str")


    for analysis in analyses:
        if isinstance(analysis, str): # if it's a string, then it's a static analysis
            headers.extend(STATIC_ANALYSES[analysis]["columns"])
            column_types.extend(STATIC_ANALYSES[analysis]["column_types"])
            analyses_to_run.append(STATIC_ANALYSES[analysis])
        else:
            generated_analysis = DYNAMIC_ANALYSES[analysis[0]](analysis[1])
            headers.extend(generated_analysis["columns"])
            column_types.extend(generated_analysis["column_types"])
            analyses_to_run.append(generated_analysis)

    
    if column_names is not None:
        if len(headers) != len(column_names):
            raise ValueError('The number of override column names must be the same as the number of columns.\nColumns:\n' + "; ".join(headers) + "\nOverride columns:\n" + "; ".join(column_names))
        headers = column_names

    writers = list()
    if tsv_out_handle is not None:
        writers.append(EnumTSVWriter(headers, column_types, tsv_out_handle))
    if html_out_handle is not None:
        writers.append(EnumHTMLWriter(headers, column_types, html_out_handle, html_max_height))

    for writer in writers:
        writer.write_header()
    

    for rec in records:
        if databases is not None:
            rec = tuple(filter_domains((rec,), evalue=float("inf"), max_overlap=1, databases_keep=databases))[0]
        for out_line in process_record(rec, by, analyses_to_run, ncbi_taxonomy):
            for writer in writers:
                writer.write_row(out_line)
    
    for writer in writers:
        writer.write_footer()

def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', nargs='+', required=False, default=None,
                        help="Genbank filenames. If not supplied then reads from to stdin.")

    parser.add_argument('-o', '--output', default=None, required=False, 
                        help="The name of the output tsv file.")
    parser.add_argument('--html', default=None, required=False, 
                        help="Write a table in html format to this file.")
    parser.add_argument('--html_max_height', default=None, required=False, type=int,
                        help="Max height in pixels to set the html output to.")

    parser.add_argument('--contigs', default=[], nargs='+',
                        help="only examine contigs with ids in this list. Additive with --contigs_file. default: all contigs.")

    parser.add_argument('--contigs_file', default=None, 
                        help="only examine contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs.")

    parser.add_argument('--by', type=str.lower, default="contig", choices=SELECT_CATEGORIES,
                        help="One line in output for every by contig, cds or domain. For protein genbanks, contig and cds are treated the same. default: by contig")

    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Ignore domain annotations not from these databases. default: consider all databases.")

    parser.add_argument('--domains', action='append_const', dest=COLS_ARG_NAME, const="domains", #It's alphabetical so that different contigs can be compared
                        help="report domain content as alphabetical and semicolon-separated, only listing duplicates once.")

    parser.add_argument('--domain_descriptions', action='append_const', dest=COLS_ARG_NAME, const="domain_descriptions",
                        help="report domain content as alphabetical and semicolon-separated, only listing duplicates once.")
    
    parser.add_argument('--score', action='append_const', dest=COLS_ARG_NAME, const="score",
                        help="report report the sum of domain scores.")

    parser.add_argument('--identity', action='append_const', dest=COLS_ARG_NAME, const="identity",
                        help="report report the sum of domain identities.")

    parser.add_argument('--architecture', action='append_const', dest=COLS_ARG_NAME, const="architecture",
                       help="report domains in order of appearance of start position in the sequence and semicolon-separated.")
    
    parser.add_argument('--architecture_detailed', action='append_const', dest=COLS_ARG_NAME, const="architecture_detailed",
                       help="report domain architecture in order of appearance of start position in the sequence. Within each CDS, domains are separated by a '/'. CDSs are separated by a ';'. CDSs with no domains are represented by 'no_domains'.")

    parser.add_argument('--start', action='append_const', dest=COLS_ARG_NAME, const="start",
                        help="report the start coordinate of the feature.")
    
    parser.add_argument('--end', action='append_const', dest=COLS_ARG_NAME, const="end",
                        help="report the end coordinate of the feature.")
    
    parser.add_argument('--strand', action='append_const', dest=COLS_ARG_NAME, const="strand",
                        help="report the strand of the feature.")

    parser.add_argument('--length', action='append_const', dest=COLS_ARG_NAME, const="length",
                        help="report the (contig, cds, or domain) length as a column in the output.")
    
    parser.add_argument('--topology', action='append_const', dest=COLS_ARG_NAME, const="topology",
                        help="report the topology in the output (typically only makes sense when using --by contig).")
    
    parser.add_argument('--sequence', action='append_const', dest=COLS_ARG_NAME, const="sequence",
                        help="report the (contig, cds, or domain) sequence as a column in the output.")
    
    parser.add_argument('--cds_count', action='append_const', dest=COLS_ARG_NAME, const="cds_count",
                        help="report the number of CDSs as a column in the output")
    
    parser.add_argument('--domain_count', action='append_const', dest=COLS_ARG_NAME, const="domain_count",
                        help="report the number of domains as a column in the output")
    
    parser.add_argument('--definition', action='append_const', dest=COLS_ARG_NAME, const="definition",
                        help="report sequence definition line from the genbank file.")
    
    parser.add_argument('--filename', action='append_const', dest=COLS_ARG_NAME, const="filename",
                        help="report the filename of the input file from which the record originated.")
    
    parser.add_argument('--domain_search', action='append_const', dest=COLS_ARG_NAME, const="domain_search",
                        help="This is a composite metric that will return multiple columns: best_query, query_description, query_start, query_end, query_length, match_start, match_end, match_strand, match_score, match_identity.")
    
    parser.add_argument('--taxid_lineage', action='append_const', dest=COLS_ARG_NAME, const="taxid_lineage", help="a semicolon separated list of taxids from the root to the taxid of the sequence.") 
    parser.add_argument('--taxname_lineage', action='append_const', dest=COLS_ARG_NAME, const="taxname_lineage", help="a semicolon separated list of taxnames from the root to the taxid of the sequence.") 
    parser.add_argument('--rank_lineage', action='append_const', dest=COLS_ARG_NAME, const="rank_lineage", help="a semicolon separated list of ranks from the root to the taxid of the sequence.")
    
    parser.add_argument('--taxid', nargs="+", action=DynamicArg, const="taxid", dest=COLS_ARG_NAME, help="The taxid at the given rank (e.g. superkingdom, genus, species), use 'self' for the lowest level taxid available for the record, the column will be named 'taxid_[rank]") # 
    parser.add_argument('--taxname', nargs="+", action=DynamicArg, const="taxname", dest=COLS_ARG_NAME, help="The taxname at the given rank (e.g. superkingdom, genus, species), use 'self' for the lowest level taxname available for the record, the column will be named 'taxname_[rank]") #

    parser.add_argument("--taxonomy_update", action="store_true", help="If taxonomy database exists, check it against the version on the ncbi server and update if there is a newer version.")
    parser.add_argument("--ncbi_taxonomy_path", type=str,  default="/tmp/ncbi_taxonomy", help="Path to NCBI taxonomy database directory. Will be created and downloaded if it does not exist.")

    parser.add_argument('--qualifier', nargs=2, required=False, action=DynamicArg, dest=COLS_ARG_NAME, const="qualifier",
                        help="Supply two strings, an feature type and a qualifier name. Report the values of all instances of this feature and qualifier combination. If multiple are found, they are separated by a ';'.")
    parser.add_argument('--feature_count', nargs=1, required=False, action=DynamicArg, dest=COLS_ARG_NAME, const="feature_count",
                        help="Reports the number of features of this type.")
    # parser.add_argument('--qualifier_count', type=str, nargs=2, required=False, action=DynamicArg, dest=COLS_ARG_NAME, #TODO: implement this using DynamicArg
    #                     help="Supply two strings, an feature type and a qualifier name. Report the number of instances of this feature-qualifier combination.")
    
    parser.add_argument('--append', nargs=3, required=False, action=DynamicArg, dest=COLS_ARG_NAME, const="append",
                        help="Supply three strings, a column will be added with the first string as the column name, the second string as the column type (str, int, float) and the third string as the value for all rows.")

    parser.add_argument('--column_names', nargs='+', default=None, required=False, type=str,
                        help="If supplied, then this list will be used instead of the default column names.")
    
    parser.add_argument('--config', action=ActionConfigFile)



    params = parser.parse_args(argv)
   
    ### Figure out what input and output files ####
    filetype_override = None
    if params.input is None:
        genbanks = [sys.stdin]
        filetype_override = "genbank"
    else:
        genbanks = params.input
    
    out = None
    if params.output is not None:
        out = open(params.output, "w")

    html_out_handle = None
    if params.html is not None:
        html_out_handle = open(params.html, "w")

    if out is None and html_out_handle is None:
        raise ValueError("You must supply either --output or --html or both.")

    #print(getattr(params, COLS_ARG_NAME))
    # figure out other parameters

    target_contigs = list_and_file_to_dict_keys(
        params.contigs, params.contigs_file, as_set=True)

    
    if hasattr(params, COLS_ARG_NAME): 
        analyses = getattr(params, COLS_ARG_NAME)
    if analyses is None:
        analyses = list()
    analysis_names = get_analysis_names(analyses)
    ncbi_taxonomy = None
    if "taxid_lineage" in analysis_names or "taxname_lineage" in analysis_names or "taxname" in analysis_names or "taxid" in analysis_names:
        # create the path to the NCBI taxonomy database
        Path(params.ncbi_taxonomy_path).mkdir(parents=True, exist_ok=True)
        # load the NCBI taxonomy database
        ncbi_taxonomy = NCBITaxonomy(params.ncbi_taxonomy_path, params.taxonomy_update)

    databases=None
    if params.databases is not None:
        databases = set(params.databases)

    # Determine whether to track filenames
    use_filename_tracking = "filename" in analysis_names
    
    # Choose the appropriate record parser
    if use_filename_tracking:
        records_iter = parse_seqfiles_with_filenames(genbanks, target_contigs, filetype_override=filetype_override)
    else:
        records_iter = parse_seqfiles(genbanks, target_contigs, filetype_override=filetype_override)

    # Run
    enum_report(records_iter,
                params.by,
                analyses,
                out,
                html_out_handle,
                params.column_names,
                params.html_max_height,
                ncbi_taxonomy=ncbi_taxonomy,
                databases=databases
                )

    if params.output is not None:
        out.close()
    
    if params.html is not None:
        html_out_handle.close()


def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
