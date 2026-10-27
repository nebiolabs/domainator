"""Create graphical representations of the contigs and their annotations

It works best for relatively small datasets, of a few megabases or a few hundred contigs.

"""
# The page is web/reports/plot_contigs.html and .js (see report_html.py). This module builds its
# data, and the page lays the plot out: it stacks each contig's features and spaces the contigs.

import csv
from collections import defaultdict
from os import PathLike
from typing import Dict, List, Optional, Tuple, Union

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from jsonargparse import ActionConfigFile

from domainator import DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME, report_html
from domainator.cli.common import add_input, add_output, check_outputs, declare_outputs, tool_parser
from domainator.colors import OTHER_COLOR, get_palette, normalize_color_hex, read_color_table
from domainator.errors import InputError
from domainator.model import get_cds_name
from domainator.outputs import output_file
from domainator.seqio import open_if_is_name, parse_seqfiles

#: The feature types drawn by default, in stacking order.
DEFAULT_TYPES = ["CDS", "exon", "repeat_region", "ncRNA", "rRNA", "tRNA", DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME]

#: The color of a contig's line when it is not colored by metadata, and of a value a color
#: table does not list.
DEFAULT_CONTIG_COLOR = "#000000"

#: The record annotations a contig's details list, in this order, when they have a value.
#: References and structured comments are left out.
CONTIG_ANNOTATIONS = ("molecule_type", "topology", "data_file_division", "date", "accessions", "sequence_version",
                      "keywords", "source", "organism", "taxonomy", "comment")


def _strand(strand: Optional[int]) -> str:
    return "-" if strand == -1 else "+"


def get_feature_name(feature: SeqFeature) -> str:
    if feature.type == "CDS" or feature.type == "gene":
        return get_cds_name(feature)
    if feature.type == "source" and "organism" in feature.qualifiers:
        return feature.qualifiers["organism"][0]
    for key in ("name", "product", "gene", "locus_tag"):
        if key in feature.qualifiers:
            return feature.qualifiers[key][0]
    return feature.type


def feature_to_dict(feature: SeqFeature, default_color: Optional[Union[str, Dict[str, str]]] = None) -> dict:
    """A feature as the page draws it: {name, type, color, parts, qualifiers}.

    parts are [start, end, strand], 0-based and end-exclusive, strand "+" or "-"; qualifiers
    are [key, [values]] pairs.

    Args:
        feature: The SeqFeature to convert
        default_color: The color of a feature with no color or Color qualifier. Can be a hex
            string, eg "CC00CC", or a dict of feature type to hex string, eg
            {"CDS": "CCCCCC", "exon": "00FF00"}, types missing from it being "#CCCCCC".
    """
    if default_color is None:
        default_color = defaultdict(lambda: "#CCCCCC", {"CDS": "#CCCCCC", "exon": "#CCCCCC", "repeat_region": "#0000FF", DOMAIN_FEATURE_NAME: "#00FFFF", DOMAIN_SEARCH_BEST_HIT_NAME: "#06FF0D", "source": "#2626FF"})
    elif isinstance(default_color, defaultdict):
        pass
    elif isinstance(default_color, dict):
        default_color = {k: normalize_color_hex(v) for k, v in default_color.items()}
        default_color = defaultdict(lambda: "#CCCCCC", default_color)
    else: # assume it's a hex string
        color = normalize_color_hex(default_color)
        default_color = defaultdict(lambda: color)

    if "color" in feature.qualifiers:
        color = normalize_color_hex(feature.qualifiers["color"][0])
    elif "Color" in feature.qualifiers:
        color = normalize_color_hex(feature.qualifiers["Color"][0])
    else:
        color = default_color[feature.type]

    return {
        "name": get_feature_name(feature),
        "type": feature.type,
        "color": color,
        "parts": [[int(part.start), int(part.end), _strand(part.strand)] for part in feature.location.parts],
        "qualifiers": [[key, [str(value) for value in values]] for key, values in feature.qualifiers.items()],
    }


def _annotation_text(value) -> str:
    if isinstance(value, (list, tuple)):
        return "; ".join(str(item) for item in value)
    return str(value)


def contig_details(record: SeqRecord, metadata: Optional[List[Tuple[str, str]]] = None) -> List[List[str]]:
    """What the page lists for a contig, as [label, text] pairs: its name, description and
    length, the CONTIG_ANNOTATIONS that have a value, its cross-references, and ``metadata``."""
    details = [["Name", record.id], ["Description", record.description], ["Length", str(len(record))]]
    for key in CONTIG_ANNOTATIONS:
        text = _annotation_text(record.annotations[key]) if key in record.annotations else ""
        if text != "":
            details.append([key, text])
    for xref in record.dbxrefs:
        db, colon, identifier = xref.partition(":")
        details.append([db.strip(), identifier.strip()] if colon else ["dbxref", xref])
    details.extend([column, value] for column, value in metadata or [])
    return details


def read_metadata(paths: List[Union[str, PathLike]]) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    """Tab-separated tables of contig metadata, as (columns, {contig id: {column: value}}).

    Each table's first column is the contig id and its first row the column names. Values are
    kept as the text in the file, and an empty cell is no value. Within a table the first row
    for an id counts; when tables share a column, a contig gets the value of the first table
    that has one. Columns are listed in the order they are first seen.
    """
    columns = []
    values = {}
    for path in paths:
        with open(path, newline="") as handle:
            reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
            names = next(reader, [])[1:]
            columns.extend(name for name in names if name not in columns)
            seen = set()
            for row in reader:
                if not row or row[0] in seen:
                    continue
                seen.add(row[0])
                contig = values.setdefault(row[0], {})
                for name, value in zip(names, row[1:]):
                    if value != "" and name not in contig:
                        contig[name] = value
    return columns, values


def plot_contigs(records: List[SeqRecord], *,
                 html_out: str,
                 metadata: Optional[List[Union[str, PathLike]]] = None,
                 color_by: Optional[str] = None,
                 color_table: Optional[Dict[str, str]] = None,
                 palette_offset: int = 0,
                 height: int = 2400, width: int = 1200,
                 types: Optional[List[str]] = None,
                 feature_height: int = 50, arrow_width: int = 20,
                 contig_thickness: int = 20,
                 contig_label_font: str = "sans-serif", contig_label_px: int = 32,
                 feature_label_font: str = "sans-serif", feature_font_size: int = 26,
                 details_height: int = 100, details_width: Optional[int] = None
                 ) -> None:
    """Write an interactive HTML plot of ``records`` to ``html_out`` (a path or a text handle).

    The options mean the same as the ``domainator plot_contigs`` options with the same names;
    ``types`` is the feature types to draw, and ``palette_offset`` is ignored when
    ``color_table`` is given. A contig with no ``color_by`` value is drawn in OTHER_COLOR, and
    one whose value ``color_table`` does not list in DEFAULT_CONTIG_COLOR.
    """
    if types is None:
        types = DEFAULT_TYPES
    columns, values = read_metadata(metadata) if metadata is not None else ([], {})
    if color_by is not None and color_by not in columns:
        raise InputError(f"--color_by {color_by} is not a column of the --metadata tables.")

    contig_metadata = [[(column, values[record.id][column]) for column in columns if column in values.get(record.id, {})]
                       for record in records]
    colors = [DEFAULT_CONTIG_COLOR] * len(records)
    if color_by is not None:
        by = [dict(row).get(color_by) for row in contig_metadata]
        present = list(dict.fromkeys(value for value in by if value is not None))
        if color_table is None:
            palette = get_palette(present, offset=palette_offset)
        else:
            palette = {value: color_table.get(value, DEFAULT_CONTIG_COLOR) for value in present}
        colors = [OTHER_COLOR if value is None else palette[value] for value in by]

    data = {
        "types": list(types),
        "color_by": color_by,
        "layout": {
            "width": int(width), "height": int(height), "feature_height": int(feature_height),
            "arrow_width": int(arrow_width), "contig_thickness": int(contig_thickness),
            "contig_label_font": contig_label_font, "contig_label_px": int(contig_label_px),
            "feature_label_font": feature_label_font, "feature_font_size": int(feature_font_size),
            "details_width": int(width if details_width is None else details_width),
            "details_height": int(details_height),
        },
        "contigs": [{
            "name": record.id,
            "length": len(record),
            "color": color,
            "details": contig_details(record, row),
            "features": [feature_to_dict(feature) for feature in record.features if feature.type in types],
        } for record, row, color in zip(records, contig_metadata, colors)],
    }
    page = report_html.render_report("plot_contigs", "Contig Plot", "plot_contigs", data, encoding="gzip+base64")
    outfile_handle, file_type = open_if_is_name(html_out, "w")
    outfile_handle.write(page)
    if file_type == "name":
        outfile_handle.close()


def show_contigs(recs:Union[str,List[SeqRecord]], window_height=600, window_width=1200, **kwargs):
    """
        Show the contigs in a Jupyter notebook.
        NOTE: this function creates some temporary files in the current working director that are not deleted. 
              they will be of the form "domainator_contig_plot_*.html"
        NOTE 2: seems to work in-line in Jupyter lab, but not in vscode notebooks or in Colab.
                For those platforms, you can open the created html file in a web browser,
                or use the plot_contigs function to create the html file with a specified name.
    """
    from IPython.display import IFrame, display
    # capture the output of plot_contigs as a string
    import tempfile
    import os
    
    if isinstance(recs, str):
        recs = list(parse_seqfiles([recs]))
    
    with tempfile.NamedTemporaryFile(mode="w", prefix="domainator_contig_plot_", suffix=".html", dir=os.getcwd(), delete=False) as s:
        plot_contigs(recs, html_out=s.name, height=window_height, width=window_width, **kwargs)
        print(os.path.basename(s.name))
        display(IFrame(os.path.basename(s.name), width=window_width, height=window_height))


def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, "names of input genbank files.")
    add_output(parser, "html", "Write the plot to this html file.", primary=True, required=True)
    parser.add_argument("--height", default=2400, type=int, required=False,
                        help="height of the contig plot, in pixels")
    parser.add_argument("--width", default=1200, type=int, required=False,
                        help="width of the contig plot, in pixels")
    parser.add_argument("--types", default=None, type=str, required=False, nargs="+",
                        help="Feature types to include in the plot. Default is 'CDS exon repeat_region ncRNA rRNA tRNA Domainator Domain_Search'.")
    parser.add_argument("--feature_height", default=50, type=int, required=False,
                        help="height of the feature arrows, in pixels.")
    parser.add_argument("--feature_font_size", default=26, type=int, required=False,
                        help="font size of the feature labels, in pixels.")

    parser.add_argument('--metadata', type=str, nargs="+", required=False, default=None, 
                        help="tab separated files of additional contig metadata.")

    parser.add_argument('--color_by', type=str, required=False, default=None,
                        help="Color the contig lines in the plot based on this column of the metadata table.")
    
    parser.add_argument("--color_table", required=False, default=None, type=str,
                        help="tab separated file with two columns and no header: a --color_by value and a hex color, for example: CCDB   cc0000. "
                             "Without it, colors are assigned automatically.")
    parser.add_argument("--palette_offset", required=False, default=0, type=int,
                        help="Skip this many colors at the start of the palette before assigning colors to --color_by values. "
                             "For example, if a plot colored by one column used the first 5 colors, --palette_offset 5 colors another "
                             "column without reusing them. Ignored with --color_table.")
    # TODO: change defaults for proteins?
    

    parser.add_argument('--config', action=ActionConfigFile)
    declare_outputs(parser, ["html_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)

    ### validate input


    check_outputs(parser, params)

    if params.color_table is not None and params.color_by is None:
        parser.error("--color_table requires --color_by.")
    if params.color_by is not None and params.metadata is None:
        parser.error("--color_by requires --metadata.")
    color_table = read_color_table(params.color_table) if params.color_table is not None else None

    ### Run
    recs = list(parse_seqfiles(params.input))
    
    with output_file(params.html_out) as html_handle:
        plot_contigs(recs, html_out=html_handle, metadata=params.metadata, color_by=params.color_by, color_table=color_table,
                     palette_offset=params.palette_offset,
                     feature_height=params.feature_height, feature_font_size=params.feature_font_size,
                     height=params.height, width=params.width, types=params.types)
