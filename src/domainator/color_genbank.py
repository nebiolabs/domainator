"""add colors to domains or CDS annotations in genbank files

Takes a genbank file which has been annotated using the domainator program, adds color annotations to the sequence features.

Domains are colored based on domain name.

CDSs are colored based on the associated domains, prioritized by domain score.

If a color_table is provided, then the colors in the table will be used, otherwise a random color will be assigned to each domain.

"""

import sys
import argparse
from domainator.utils import list_and_file_to_dict_keys, parse_seqfiles, write_genbank, open_if_is_name
import warnings
import re
import seaborn as sns
from collections import OrderedDict
import random
from domainator import __version__, DOMAIN_FEATURE_NAME, RawAndDefaultsFormatter, DOMAIN_SEARCH_BEST_HIT_NAME

#TODO: color gradient based on bit score

def normalize_color_hex(color_hex):
    '''
        takes in a string of the format:
            (#?)([a-fA-F0-9]{6})
        raises an exception if string is not of expected format
        otherwise, returns string in format #[A-F0-9]{6}
    '''
    color_hex = color_hex.strip().upper()
    matches = re.match(r'^(#?)([A-F0-9]{6})$',color_hex)
    if not matches:
        raise RuntimeError(f"Color code not recognized: {color_hex}")
    else:
        return "#" + matches.group(2)
    
def read_color_table(color_table):
    """
        input: path or handle to a tsv with two columns and no header
                columns are: domain id, hex color. For example: CCDB   cc0000
        output:
            a dict of domain: hex color
                for example: {"CCDB": "#CC0000"}
    """
    
    file_handle, input_type = open_if_is_name(color_table)
    out = dict()

    for line in file_handle:
        line = line.strip()
        if len(line) >0:
            parts = line.split("\t")
            if len(parts) != 2:
                raise RuntimeError(f"badly formatted line in color table: {line}")
            if parts[0] not in out:
                out[parts[0]] = normalize_color_hex(parts[1])
            else:
                raise RuntimeError(f"domain specified multiple times in color table: {parts[0]}")

    if input_type == "name":
        file_handle.close()
    
    return out

def color_seqrecord_by_domain(record, color_map, color_targets, disable_warnings=False, search_hit_color=None, clear=False):
    """
        input: 
            record: SeqRecord containing CDS and domain annotations from Domainator program
            color_map: dict of {domain_name: hex_color}, 
                    for example: {'RNA_pol': '#cc0000', 'DNA_pol_A': '#0000cc', 'DnaB_C':'#FFA500', 'Terminase':'#00cc00'}
            color_targets: a set or list containing "domains", "cdss", or both. Specifying which type of feature to color.
            disable_warnings: if True, then don't warn when a CDS gets multiple colors. NOT IMPLEMENTED
            search_hit_color: if not None, then color search_hit annotations with this color.
            clear: if True, then remove all color annotations from the genbank file before adding any new ones.

        output:
            modifies the SeqRecord to add "Color" qualifiers
    """

    new_cds_colors = dict()
    new_cds_bitscores = dict()
    new_search_hit_colors = dict()
    ninf = float("-inf")
    for feature in record.features:
        if clear:
            if "Color" in feature.qualifiers:
                del feature.qualifiers['Color']
        if feature.type == DOMAIN_FEATURE_NAME:
            if feature.qualifiers['name'][0] in color_map:
                if "domains" in color_targets:
                    feature.qualifiers['Color'] = [ color_map[feature.qualifiers['name'][0]] ]

                score = float(feature.qualifiers['score'][0])
                if "cdss" in color_targets:
                    cds_id = feature.qualifiers['cds_id'][0]
                    
                    if score > new_cds_bitscores.get(cds_id, ninf):
                        new_cds_bitscores[cds_id] = score
                        new_cds_colors[cds_id] = color_map[feature.qualifiers['name'][0]]
                    # elif not disable_warnings:
                    #     warnings.warn(f"cds {cds_id} in {record.id} has more than one color annotation")
        if feature.type == DOMAIN_SEARCH_BEST_HIT_NAME:
            if search_hit_color is not None:
                if "domains" in color_targets:
                    feature.qualifiers['Color'] = [search_hit_color]
                if "cdss" in color_targets:
                    new_search_hit_colors[feature.qualifiers['cds_id'][0]] = search_hit_color
    
    for cds, color in new_search_hit_colors.items():
        new_cds_colors[cds] = color

    if "cdss" in color_targets:
        for feature in record.features:
            if feature.type == "CDS" and 'cds_id' in feature.qualifiers:
                if feature.qualifiers['cds_id'][0] in new_cds_colors:
                    feature.qualifiers['Color'] = [new_cds_colors[feature.qualifiers['cds_id'][0]]]

def auto_color_map(genbanks, contigs):
    domains = OrderedDict()

    for rec in parse_seqfiles(genbanks, contigs):
        for feature in rec.features:
            if feature.type == DOMAIN_FEATURE_NAME:
                domains[feature.qualifiers['name'][0]] = None

    colors = sns.color_palette('husl', n_colors=len(domains)).as_hex()
    random.Random(12).shuffle(colors)
    for index, domain in enumerate(domains.keys()):
        domains[domain] = colors[index].upper()
    
    return domains

class AutoColormap:
    def __init__(self, num_colors=30):
        self.next_color = 0
        self.num_colors = num_colors
        self.color_list = [x.upper() for x in sns.color_palette('husl', n_colors=self.num_colors).as_hex()]

        self.color_map = dict()
        random.Random(12).shuffle(self.color_list)
    def __getitem__(self, key: str) -> str:
        if key not in self.color_map:
            self.color_map[key] = self.color_list[self.next_color]
            self.next_color = (self.next_color + 1)  % self.num_colors
        return self.color_map[key]
    def __contains__(self, key: str) -> bool:
        x = self[key]
        return True
    def items(self):
        return self.color_map.items()

        

def color_genbank(records, color_table, color_targets, search_hit_color=None, clear=False, color_table_out=None):
    """
      input: 
        records: an iterator of SeqRecords
        contigs: collection if not None, then only . if None then don't filter out any contigs (most common).
        color_table: a dict mapping domain names to color hex.
        color_targets: a set or list containing "domains", "cdss", or both. Specifying which type of feature to color.
        search_hit_color: if not None, then color search_hit annotations with this color.
        clear: if True, then remove all color annotations from the genbank file before adding any new ones.
    """
    disable_multiple_colors_warning = False
    if color_table is not None:
        color_map = color_table
    else:
        disable_multiple_colors_warning = True
        color_map = AutoColormap()

    for rec in records:
        color_seqrecord_by_domain(rec, color_map, color_targets, disable_multiple_colors_warning, search_hit_color=search_hit_color, clear=clear)
        yield rec
    
    if color_table_out is not None:
        with open(color_table_out, "w") as out:
            for domain, color in color_map.items():
                out.write(f"{domain}\t{color}\n")

def main(argv):
    parser = argparse.ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', nargs='+', default=None, required=False,
                        help="Genbank filenames. If not supplied, reads from stdin.")

    parser.add_argument('-o', '--output', default=None, required=False,
                        help="the name of the output genbank file. If not supplied, writes to stdout.")

    parser.add_argument('--contigs', default=[], nargs='+',
                        help="only color contigs with ids in this list (other contigs will be excluded from output). Additive with --contigs_file. default: all contigs.")
    parser.add_argument('--contigs_file', default=None, 
                        help="only color contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs.")

    parser.add_argument("--color_table", required=False, default=None, type=str, help="tab separated file with two columns and no header, columns are: domain id, hex color. For example: CCDB   cc0000")
    parser.add_argument("--color_table_out", required=False, default=None, type=str, help="tab separated file with two columns and no header, columns are: annotation, hex color. Written after the color table is updated with new colors, for example if using --color_by, but not supplying an external color table.")

    parser.add_argument("--search_hit_color", default=None, type=str, help="color to use for search hit annotations.")

    parser.add_argument("--clear", action="store_true", default=False,
                        help="if set, then remove all color annotations from the genbank file before adding any new ones.")

    #TODO: maybe just color the top 16 or so most common domains, and give everything else the same color? Or have like 5 colors split among all of the less common domains.
    color_selection = parser.add_mutually_exclusive_group() # TODO: change this from mutually exclusive to just a list of options on one parameter, 
    color_selection.add_argument("--color_domains", action="store_true", default=False, 
                                help="if set, then color domain annotations instead of CDS annotations (by default CDS annotations will be colored)")
    color_selection.add_argument("--color_both", action="store_true", default=False, 
                                help="if set, then color domain and CDS annotations annotations (by default only CDS annotations will be colored)")
 

    params = parser.parse_args(argv)
   
    ### Figure out what input and output files ####

    filetype_override = None
    if params.input is None:
        genbanks = [sys.stdin]
        filetype_override = "genbank"
    else:
        genbanks = params.input


    if params.output is None:
        out = sys.stdout
    else:
        out = open(params.output, "w")

    # figure out other parameters

    target_contigs = list_and_file_to_dict_keys(
        params.contigs, params.contigs_file, as_set=True)

    color_targets = {}
    if params.color_domains:
        color_targets = {"domains"}
    elif params.color_both:
        color_targets = {"domains", "cdss"}
    else:
        color_targets = {"cdss"}

    color_table = None
    if params.color_table is not None:
        color_table = read_color_table(params.color_table)

    search_hit_color = params.search_hit_color 
    if search_hit_color is not None:
        search_hit_color = normalize_color_hex(search_hit_color)

    # Run
    write_genbank(color_genbank(
            parse_seqfiles(genbanks, target_contigs, filetype_override=filetype_override),
            color_table, 
            color_targets,
            search_hit_color,
            params.clear,
            color_table_out=params.color_table_out
            ), 
        out)
    

    if params.output is not None:
        out.close()

def _entrypoint():
    main(sys.argv[1:])
    
if __name__ == '__main__':
    main(sys.argv[1:])
