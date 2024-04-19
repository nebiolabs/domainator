"""Create graphical representations of the contigs and their annotations.

It works best for relatively small datasets, of a few megabases or a few hundred contigs.

WARNING: THIS IS A WORK IN PROGRESS. THE API MAY CHANGE.

"""

# TODO: more command line options
# TODO: better support for eukaryotic genomes

from jsonargparse import ArgumentParser, ActionConfigFile
import sys
from typing import Tuple, List, Optional, Set, NamedTuple, Union, Dict
from os import PathLike
from domainator import __version__, RawAndDefaultsFormatter, DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME
from domainator.utils import parse_seqfiles, get_cds_name, open_if_is_name
from domainator.Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from domainator.Bio.SeqRecord import SeqRecord
from collections import defaultdict
import html
import json
from domainator.color_genbank import normalize_color_hex
from pprint import pprint
import pandas as pd
from domainator.utils import get_palette

def normalize_strand(strand:int) -> str:
    if strand is None:
        return "+"
    elif int(strand) == -1:
        return "-"
    else:
        return "+"

def location_to_list_of_dicts(location:Union[SimpleLocation, CompoundLocation]) -> dict:
    out = []
    for part in location.parts:
        out.append({
            "start": part.start,
            "end": part.end,
            "strand": normalize_strand(location.strand)
        })
    return out

def get_feature_mouseover(feature:SeqFeature) -> str:
    return "<br>".join([f'<b style="color:#AA2222;">{html.escape(feature.type)}</b>', 
                        f'<b>start</b>: {html.escape(str(feature.location.start))}',
                        f'<b>end</b>: {html.escape(str(feature.location.end))}',
                        f'<b>strand</b>: {html.escape(normalize_strand(feature.location.strand))}',]
                       + [f"<b>{html.escape(k)}</b>: {html.escape(';'.join(v))}" for k, v in feature.qualifiers.items()]
                       )

def get_contig_mouseover(rec:SeqRecord, metadata=None, color_by=None, palette=None) -> str:
    annots = [("Name", rec.id), 
              ("Description", rec.description),
              ("Length", len(rec)),
             ]
    annots += [(k, v) for k, v in rec.annotations.items()]
    for xref in rec.dbxrefs:
        if ":" in xref:
            db, id = xref.split(":", 1)
            annots.append((db, id))
        else:
            annots.append(("dbxref", xref))
    html_escaped = [f"<b>{html.escape(str(k))}</b>: {html.escape(str(v))}" for (k, v) in annots]


    if metadata is not None:
        if rec.id in metadata.index:
            row = metadata.loc[rec.id]
            for k, v in row.items():
                if color_by is not None and k == color_by:
                    html_escaped.append(f'<b>{html.escape(str(k))}</b>: <span style="color:{palette[v]};">{html.escape(str(v))} </span>')
                else:
                    html_escaped.append(f"<b>{html.escape(k)}</b>: {html.escape(str(v))}")

        #html_escaped.append(f"<b>{html.escape(k)}</b>: {html.escape(str(v))}")
    return "<br>".join(
        [f'<b style="color:#AA2222;">Contig</b>'] + html_escaped
        )

def get_feature_name(feature:SeqFeature) -> str:
    name = ""
    if feature.type == "CDS" or feature.type == "gene":
        name += get_cds_name(feature)
    elif feature.type == "source" and "organism" in feature.qualifiers:
        name += feature.qualifiers["organism"][0]
    elif "name" in feature.qualifiers:
        name += feature.qualifiers["name"][0]
    elif "product" in feature.qualifiers:
        name += feature.qualifiers["product"][0]
    elif "gene" in feature.qualifiers:
        name += feature.qualifiers["gene"][0]
    elif "locus_tag" in feature.qualifiers:
        name += feature.qualifiers["locus_tag"][0]
    else:
        name += feature.type
    return name


# same thing as above but as a dict initialized with the Feature fields
def feature_to_dict(feature: SeqFeature, default_color:Optional[Union[str, Dict[str, str]]]=None) -> dict:
    """
    Convert a SeqFeature to a dictionary that can be serialized to JSON.

    Args:
        feature: The SeqFeature to convert
        default_color: The default color to use for the feature. Can be a hex string, eg "CC00CC". 
                        Or a dict of feature type to hex string, eg {"CDS": "CCCCCC", "exon": "00FF00", "repeat_region": "0000FF",}.
    """
    if default_color is None:
        default_color = defaultdict(lambda: "#CCCCCC", {"CDS": "#CCCCCC", "exon": "#CCCCCC", "repeat_region": "#0000FF", DOMAIN_FEATURE_NAME: "#00FFFF", DOMAIN_SEARCH_BEST_HIT_NAME: "#06FF0D", "source": "#2626FF"})
    elif isinstance(default_color, defaultdict):
        default_color = default_color
    elif isinstance(default_color, dict):
        default_color = {k: normalize_color_hex(v) for k, v in default_color.items()}
        default_color = defaultdict(lambda: "#CCCCCC", default_color)
    else: # assume it's a hex string
        default_color = defaultdict(normalize_color_hex(default_color))

    if "color" in feature.qualifiers:
        color = normalize_color_hex(feature.qualifiers["color"][0])
    elif "Color" in feature.qualifiers:
        color = normalize_color_hex(feature.qualifiers["Color"][0])
    else:
        color = default_color[feature.type]

    return {
        "name": get_feature_name(feature),
        "start": feature.location.start,
        "end": feature.location.end,
        "parts": location_to_list_of_dicts(feature.location),
        "location": feature.location,
        "type": feature.type,
        "color": color,
        "mouseover": get_feature_mouseover(feature),
        "z_order": 0
    }

def add_z_order(features:List[Dict[str, Union[int, str]]], 
                precedence:List[str] = None,
                squash:Set[str] = None,
                ):
    """
        features is a list of features
        precedence is a list of feature types in order of precedence, e.g. ["source", "CDS", "repeat_region", "ncRNA", "rRNA", "tRNA"] will put "CDS" features closer to the bottom of the plot than "repeat_region" features.
        squash is a set of feature types that should be squashed together, e.g. {"CDS"} will put all CDS features on the same level as all other CDS features, regardless of their start and end positions.

        Adds a z_order attribute to each feature in the list, which is used to determine the vertical position of the feature in the plot.
    """

    if precedence is None:
        precedence = ["source", "CDS", "exon", "repeat_region", "ncRNA", "rRNA", "tRNA", DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME]
    if squash is None:
        squash = {"CDS", "exon"}


    type_order = defaultdict(lambda: len(precedence) + 1)
    type_order.update({t: i for i, t in enumerate(precedence)})

    features.sort(key=lambda x: (type_order[x["type"]], x["start"])) # sort by start position and then by type, so that features of the same type are grouped together

    # Initialize z_order calculation
    z_orders = []
    for current_feature in features:
        occupied_levels = set()
        level_override = None
        for other_feature, other_level in z_orders:
            if current_feature["type"] in squash:
                if other_feature["type"] in squash:
                    level_override = other_level
                    break
            if current_feature["location"].overlaps(other_feature["location"]):
                occupied_levels.add(other_level)

        if level_override is not None:
            current_level = level_override
        else:
            current_level = 0
            while current_level in occupied_levels:
                current_level += 1
        current_feature["z_order"] = current_level
        z_orders.append((current_feature, current_level))

def create_d3_html(output_file, contigs:List[SeqRecord],
                   metadata_files:List[Union[str, PathLike]]=None,
                   color_by:str=None,
                   window_height=2400, window_width=1200,
                   feature_types:Optional[List[str]] = None,
                   arrow_height=50, arrow_width=20,
                   contig_thickness=20,
                   contig_label_font="sans-serif", contig_label_px=32,
                   feature_label_font="sans-serif", feature_label_px=26,
                   details_height=100, details_width=None
                   ):
    
    # Convert input to correct types
    arrow_height = int(arrow_height)
    arrow_width = int(arrow_width)
    contig_label_font = html.escape(contig_label_font)
    feature_label_font = html.escape(feature_label_font)
    contig_label_px = int(contig_label_px)
    feature_label_px = int(feature_label_px)
    window_height = int(window_height)
    window_width = int(window_width)
    details_height = int(details_height)
    if details_width is None:
        details_width = window_width
    else:
        details_width = int(details_width)

    if feature_types is None:
        feature_types = ["CDS", "exon", "repeat_region", "ncRNA", "rRNA", "tRNA", DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME]

    offset_padding = contig_label_px + (contig_thickness / 2) + 5 # Adjust this value as needed

    metadata = pd.DataFrame()
    metadata.index = [rec.id for rec in contigs]
    if metadata_files is not None:
        for file in metadata_files:
            df = pd.read_csv(file, sep="\t", index_col=0)
            metadata = metadata.merge(df,how="left", left_index=True, right_index=True)
        # drop duplicate indexes
        metadata = metadata[~metadata.index.duplicated(keep='first')]

    palette = None
    if color_by is not None:
        palette = get_palette(metadata[color_by].unique())

    # Preparing data for D3
    contigs_data = []
    current_offset = 0  # Start with the initial offset for the first contig
    
    for contig in contigs:
        features = [feature_to_dict(feature) for feature in contig.features if feature.type in feature_types]
        add_z_order(features, feature_types)
        max_z_order = max([feature["z_order"] for feature in features], default=0)
        for feature in features:
            del feature["location"]  # Remove the location attribute, which is not JSON serializable

        current_offset += (max_z_order + 1) * arrow_height + offset_padding
        contigs_data.append({
            'name': contig.id,
            'length': len(contig),
            'features': features,
            'offset': current_offset,  # Add an offset for each contig based on its content
            'mouseover': get_contig_mouseover(contig, metadata=metadata, color_by=color_by, palette=palette),
            'color': "#000000",
        })
        if color_by is not None:
            contigs_data[-1]["color"] = palette[metadata.loc[contig.id, color_by]]

        # Update current_offset for the next contig, accounting for the maximum z_order of the current contig's features
        # current_offset += (max_z_order + 1) * arrow_height + offset_padding  # Adjust padding as needed for additional spacing
    
    widest_contig_length = max(contig['length'] for contig in contigs_data)
    scale_factor = min(1, window_width / widest_contig_length)
    if scale_factor < 0.1:
        scale_factor = 0.1


    arrow_width = int(arrow_width)
    arrow_height = int(arrow_height)
    # Write to file
    outfile_handle, file_type = open_if_is_name(output_file, "w")
    outfile_handle.write(f"""
    <html>
    <head>
        <title>Domainator Contigs Plot</title>
        <script src="https://d3js.org/d3.v7.min.js"></script>
        <style>
            .tooltip {{
                position: absolute;
                text-align: left;
                padding: 8px;
                font: 12px sans-serif;
                background: lightblue;
                border: 0px;
                border-radius: 8px;
                pointer-events: none;
                opacity: 0;
            }}
            .contig-label {{
                font: {contig_label_px}px {contig_label_font};
            }}
            .feature-label {{
                font: {feature_label_px}px {feature_label_font};
            }}
        </style>
    </head>
    <body>
        <div id="details" style="padding: 5px; border: 1px solid #ccc; margin-top: 5;  margin-bottom: 15px; width: {details_width}; height: {details_height}; overflow: scroll; font: 10px sans-serif;"></div>
        <div> Zoom with mousewheel, left-click to move annotations to the upper box, left-click and drag to pan </div>
        <div id="tooltip" class="tooltip"></div>
        <div>
        <svg width="{window_width}" height="{window_height}" style="padding: 5px; border: 1px solid #ccc; margin-top: 2px; overflow: scroll;"></svg>
        </div>
        <script>
            var svg = d3.select("svg");
            var detailsDiv = d3.select("#details");
            var tooltip = d3.select("#tooltip");
    
            var zoomContainer = svg.append("g");

            var zoom = d3.zoom()
                .on("zoom", function(event) {{
                    zoomContainer.attr("transform", event.transform);
                    
                }});
            var initialTransform = d3.zoomIdentity.translate(0, 0).scale({scale_factor});
            svg.call(zoom)
                .call(zoom.transform, initialTransform); // Apply the initial zoom transform correctly

            
                    
            var contigs = {json.dumps(contigs_data)};
            var contigGroup = zoomContainer.selectAll(".contig")
                .data(contigs)
                .enter()
                .append("g")
                .attr("class", "contig")
                .attr("transform", function(d, i) {{ return "translate(0," + d.offset + ")"; }});
            
    
            contigGroup.append("line")
                .attr("x1", 0)
                .attr("x2", function(d) {{ return d.length; }})
                .attr("stroke", function(d) {{ return d.color; }})
                .attr("stroke-width", {contig_thickness})
                .on("click", function(event, d) {{
                        detailsDiv.html(d.mouseover);
                }})            
                .on("mouseover", function(event, d) {{
                    tooltip.html(d.mouseover)
                        .style("opacity", 1)
                        .style("left", (event.pageX + 5) + "px")
                        .style("top", (event.pageY + 5) + "px");

                }})
                .on("mousemove", function(event, d) {{
                    tooltip.style("left", (event.pageX + 5) + "px")
                    .style("top", (event.pageY + 5) + "px");
                }})
                .on("mouseout", function() {{
                    tooltip.style("opacity", 0);
                }});

    
            var featureGroup = contigGroup.selectAll(".feature")
                .data(function(d) {{ return d.features; }})
                .enter()
                .append("g")
                .attr("class", "feature")
                .attr("transform", function(d) {{ return "translate(" + d.start + "," + (-{arrow_height} * d.z_order) + ")"; }});
    
            featureGroup.each(function(featureData) {{
                d3.select(this).selectAll(".part")
                    .data(featureData.parts)
                    .enter().append("path")
                    .attr("d", function(d) {{
                        var length = d.end - d.start;
                        var arrow_width = {arrow_width};  // Width of the arrow part
                        var arrow_height = {arrow_height};  // Height of the arrow part
                        if (d.strand === "-") {{
                            return `M${{arrow_width}},0 L${{length}},0 L${{length}},${{-arrow_height}} L${{arrow_width}},${{-arrow_height}} L0,${{-arrow_height/2}} Z`;
                        }} else {{
                            return `M0,0 L${{length - arrow_width}},0 L${{length}},${{-arrow_height/2}} L${{length - arrow_width}},${{-arrow_height}} L0,${{-arrow_height}} Z`;
                        }}
                    }})
                    .attr("fill", d => featureData.color)
                    .attr("stroke", "black")
                    .on("click", function(event, d) {{
                        detailsDiv.html(featureData.mouseover);
                    }})            
                    .on("mouseover", function(event, d) {{
                        tooltip.html(featureData.mouseover)
                            .style("opacity", 1)
                            .style("left", (event.pageX + 5) + "px")
                            .style("top", (event.pageY + 5) + "px");

                    }})
                    .on("mousemove", function(event, d) {{
                        tooltip.style("left", (event.pageX + 5) + "px")
                        .style("top", (event.pageY + 5) + "px");
                    }})
                    .on("mouseout", function() {{
                        tooltip.style("opacity", 0);
                    }});
            }});
    
            featureGroup.append("text")
                .attr("class", "feature-label")
                .attr("x", {arrow_width})
                .attr("y", -7)
                .attr("text-anchor", "start")
                .text(d => d.name)
                .on("click", function(event, d) {{
                    detailsDiv.html(d.mouseover);
                }})            
                .on("mouseover", function(event, d) {{
                    tooltip.html(d.mouseover)
                        .style("opacity", 1)
                        .style("left", (event.pageX + 5) + "px")
                        .style("top", (event.pageY + 5) + "px");
                }})
                .on("mousemove", function(event, d) {{
                    tooltip.style("left", (event.pageX + 5) + "px")
                    .style("top", (event.pageY + 5) + "px");
                }})
                .on("mouseout", function() {{
                    tooltip.style("opacity", 0);
                }});

            contigGroup.append("text")
                .attr("class", "contig-label")
                .attr("y", {contig_label_px} + {contig_thickness / 2})
                .attr("x", 0)
                .attr("text-anchor", "top")
                .text(function(d) {{ return d.name; }})
                .on("click", function(event, d) {{
                    detailsDiv.html(d.mouseover);
                }})
                .on("mouseover", function(event, d) {{
                    tooltip.html(d.mouseover)
                        .style("opacity", 1)
                        .style("left", (event.pageX + 5) + "px")
                        .style("top", (event.pageY + 5) + "px");
                }})
                .on("mousemove", function(event) {{
                    tooltip.style("left", (event.pageX + 5) + "px")
                        .style("top", (event.pageY + 5) + "px");
                }})
                .on("mouseout", function() {{
                    tooltip.style("opacity", 0);
                }});
        </script>
    </body>
    </html>
    """)
    if file_type == "name":
        outfile_handle.close()

def show_contigs(recs:Union[str,List[SeqRecord]], window_height=600, window_width=1200, **kwargs):
    """
        Show the contigs in a Jupyter notebook.
        NOTE: this function creates some temporary files in the current working director that are not deleted. 
              they will be of the form "domainator_contig_plot_*.html"
        NOTE 2: seems to work in-line in Jupyter lab, but not in vscode notebooks or in Colab.
                For those platforms, you can open the created html file in a web browser,
                or use the create_d3_html function to create the html file with a specified name.
    """
    from IPython.display import HTML, IFrame, display
    # capture the output of create_d3_html as a string
    from io import StringIO
    import tempfile
    import os
    
    if isinstance(recs, str):
        recs = list(parse_seqfiles([recs]))
    
    with tempfile.NamedTemporaryFile(mode="w", prefix="domainator_contig_plot_", suffix=".html", dir=os.getcwd(), delete=False) as s:
        create_d3_html(s.name, recs, window_height=window_height, window_width=window_width, **kwargs)
        print(os.path.basename(s.name))
        display(IFrame(os.path.basename(s.name), width=window_width, height=window_height))


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", default=None, nargs="+", type=str, required=False,
                        help="names of input genbank files. If not supplied, reads from stdin.")
    parser.add_argument("--html", default=None, type=str, required=True,
                        help="name of output html file")
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
    
    # TODO: allow user to supply a color table, like build_ssn.py
    #TODO: option to embed d3.js (v7) in the html file for offline use
    # TODO: change defaults for proteins?
    

    parser.add_argument('--config', action=ActionConfigFile)
    
    params = parser.parse_args(argv)

    ### validate input


    filetype_override = None
    if params.input is None:
        input_path = (sys.stdin,)
        filetype_override = "genbank"
    else:
        input_path = params.input

    feature_types = params.types
    if params.types is None:
        feature_types = ["CDS", "exon", "repeat_region", "ncRNA", "rRNA", "tRNA", DOMAIN_FEATURE_NAME, DOMAIN_SEARCH_BEST_HIT_NAME]
        


    ### Run
    recs = list(parse_seqfiles(input_path, filetype_override=filetype_override))
    
    create_d3_html(params.html, recs, metadata_files=params.metadata, color_by=params.color_by,
                   arrow_height=params.feature_height,
                   feature_label_px=params.feature_font_size,
                   window_height=params.height, window_width=params.width, feature_types=feature_types)

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])