"""Add colors to domains or CDS annotations in genbank files

Takes a genbank file which has been annotated using the domainator program, adds color annotations to the sequence features.

Domains are colored based on domain name.

CDSs are colored based on the associated domains, prioritized by domain score.

If a color_table is provided, then the colors in the table will be used, otherwise a random color will be assigned to each domain.

"""

from typing import AbstractSet, Dict, Iterable, Iterator, Optional, Set
from jsonargparse import ActionConfigFile
from domainator.cli.common import add_input, add_seq_outputs, check_outputs, declare_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from domainator.colors import AutoColormap, normalize_color_hex, read_color_table
from domainator.cli.common import list_and_file_to_dict_keys
from domainator.seqio import parse_seqfiles
from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.model import feature_score
from Bio.SeqRecord import SeqRecord

#TODO: color gradient based on bit score

COLOR_MODES = {
    "domains": {"domains"},
    "cdss": {"cdss"},
    "domains_and_cdss": {"domains", "cdss"},
}


def color_seqrecord_by_domain(record, color_map, color_targets, disable_warnings=False, search_hit_color=None, clear=False, databases=None, color_feature="name"):
    """
        input: 
            record: SeqRecord containing CDS and domain annotations from Domainator program
            color_map: dict of {domain_name_or_database: hex_color}, 
                    for example: {'RNA_pol': '#cc0000', 'DNA_pol_A': '#0000cc', 'DnaB_C':'#FFA500', 'Terminase':'#00cc00'}
            color_targets: a set or list containing "domains", "cdss", or both. Specifying which type of feature to color.
            disable_warnings: if True, then don't warn when a CDS gets multiple colors. NOT IMPLEMENTED
            search_hit_color: if not None, then color search_hit annotations with this color.
            clear: if True, then remove all color annotations from the genbank file before adding any new ones.
            databases: if not None, then only use domains from these databases for determining color.
            color_feature: qualifier name to use for assigning colors. Either "name" or "database".

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
            if databases is not None and feature.qualifiers['database'][0] not in databases:
                continue
            color_key = feature.qualifiers[color_feature][0]
            if color_key in color_map:
                if "domains" in color_targets:
                    feature.qualifiers['Color'] = [ color_map[color_key] ]

                score = feature_score(feature)
                if "cdss" in color_targets:
                    cds_id = feature.qualifiers['cds_id'][0]
                    
                    if score > new_cds_bitscores.get(cds_id, ninf):
                        new_cds_bitscores[cds_id] = score
                        new_cds_colors[cds_id] = color_map[color_key]
                    # elif not disable_warnings:
                    #     warnings.warn(f"cds {cds_id} in {record.id} has more than one color annotation")
        if feature.type == DOMAIN_SEARCH_BEST_HIT_NAME:
            # TODO: should "databases" be applied here?
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

def color_genbank(records: Iterable[SeqRecord], *, color_table: Optional[Dict[str, str]] = None,
                  targets: AbstractSet[str] = frozenset({"cdss"}),
                  search_hit_color: Optional[str] = None, clear: bool = False, color_table_out: Optional[str] = None,
                  databases: Optional[Set[str]] = None, color_feature: str = "name") -> Iterator[SeqRecord]:
    """Yield ``records`` with color qualifiers on their domains and/or CDSs.

    Args:
        records: the contigs or proteins (SeqRecords).
        color_table: maps domain names (or database names, with ``color_feature="database"``)
            to hex colors. None: assign colors automatically.
        targets: what to color: a set holding "domains", "cdss" (the default), or both.
        search_hit_color: color the Domain_Search best-hit features this color.
        clear: remove all existing colors first.
        color_table_out: write the color table, with any colors assigned, to this path.
        databases: color only by domains from these databases. None: all.
        color_feature: the qualifier that picks the color: "name" or "database".
    """
    disable_multiple_colors_warning = False
    if color_table is not None:
        color_map = color_table
    else:
        disable_multiple_colors_warning = True
        color_map = AutoColormap()

    for rec in records:
        color_seqrecord_by_domain(rec, color_map, targets, disable_multiple_colors_warning, search_hit_color=search_hit_color, clear=clear, databases=databases, color_feature=color_feature)
        yield rec
    
    if color_table_out is not None:
        with open(color_table_out, "w") as out:
            for annotation, color in color_map.items():
                out.write(f"{annotation}\t{color}\n")

def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, 'Genbank filenames.')

    add_seq_outputs(parser)

    parser.add_argument('--contigs', default=[], nargs='+',
                        help="only color contigs with ids in this list (other contigs will be excluded from output). Additive with --contigs_file. default: all contigs.")
    parser.add_argument('--contigs_file', default=None, 
                        help="only color contigs with ids listed in this file (one per line). Additive with --contigs. default: all contigs.")

    parser.add_argument("--color_table", required=False, default=None, type=str, help="tab separated file with two columns and no header, columns are: domain id, hex color. For example: CCDB   cc0000")
    parser.add_argument("--color_table_out", required=False, default=None, type=str, help="tab separated file with two columns and no header, columns are: annotation, hex color. Written after the color table is updated with new colors, for example if using --color_by, but not supplying an external color table.")

    parser.add_argument("--search_hit_color", default=None, type=str, help="color to use for search hit annotations.")

    parser.add_argument("--clear", action="store_true", default=False,
                        help="if set, then remove all color annotations from the genbank file before adding any new ones.")

    parser.add_argument("--databases", default=None, required=False, type=str, nargs="+",
                        help="Domain databases to use for coloring. default: all databases.")

    parser.add_argument("--target", default="cdss", choices=tuple(COLOR_MODES.keys()),
                        help="what to color. domains: color only domain annotations. cdss: color only CDS annotations. domains_and_cdss: color both domains and CDS annotations.")
    parser.add_argument("--domain_database", action="store_true", default=False,
                        help="if set, then color domains based on their database instead of their name. If used with --color_table, the first column must contain database ids.")

 
    parser.add_argument('--config', action=ActionConfigFile)
    declare_outputs(parser, ["genbank_out", "fasta_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
   
    ### Figure out what input and output files ####

    genbanks = params.input

    check_outputs(parser, params)

    # figure out other parameters

    target_contigs = list_and_file_to_dict_keys(
        params.contigs, params.contigs_file, as_set=True)

    color_targets = COLOR_MODES[params.target]
    color_feature = "database" if params.domain_database else "name"

    color_table = None
    if params.color_table is not None:
        color_table = read_color_table(params.color_table)

    search_hit_color = params.search_hit_color 
    if search_hit_color is not None:
        search_hit_color = normalize_color_hex(search_hit_color)

    databases=None
    if params.databases is not None:
        databases = set(params.databases)

    # Run
    write_seq_outputs(color_genbank(
            parse_seqfiles(genbanks, target_contigs),
            color_table=color_table, 
            targets=color_targets,
            search_hit_color=search_hit_color,
            clear=params.clear,
            color_table_out=params.color_table_out,
            databases=databases,
            color_feature=color_feature
            ), genbank=params.genbank_out, fasta=params.fasta_out)
