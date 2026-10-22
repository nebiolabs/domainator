""" Generates an SVG legend based on a color table

The color table is a tsv file that maps annotations to hex color codes.
The generated SVG file contains rectangles filled with the corresponding colors and text labels for each annotation.
"""


from jsonargparse import ActionConfigFile
import sys
from typing import Dict
from domainator.cli.common import add_input, add_output, check_outputs, declare_outputs, STDIN_PATHS, tool_parser
from domainator.outputs import output_file
from domainator.colors import read_color_table
from collections import OrderedDict
import html


def color_table_to_legend(table: Dict[str, str], *, svg_out: str, title: str) -> None:
    """Write an SVG legend of ``table`` (annotation -> hex color) to ``svg_out``, a path (``-`` for stdout)."""
    # Constants for layout and styling
    title_font_size = 24
    item_font_size = 20
    stroke_width = 2    # Width of the stroke for boxes
    padding = 10        # Padding around the content inside the legend box
    title_space = 40    # Space allocated for the title at the top of the legend

    box_height = item_font_size * 2  # Height of each color box, increased to add whitespace
    box_width = box_height  # Width of color boxes
    text_offset_y = box_height/2   # Adjust to vertically center text in the color box

    longest_key = max(len(k) for k in table.keys())

    text_offset_x = box_width + padding*2  # X offset for text to align it nicely with boxes
    total_height = title_space + len(table) * box_height + 2 * padding
    total_width = max(box_width + padding*3 + (longest_key * item_font_size)/2, padding*3 + (len(title) * title_font_size) /2) # Assumed extra space for text
    

    with output_file(svg_out) as f:
        f.write(f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg width="{total_width}" height="{total_height}" xmlns="http://www.w3.org/2000/svg">
    <g font-family="Arial, sans-serif">
        <rect x="0" y="0" width="{total_width}" height="{total_height}" fill="white" stroke="black" stroke-width="{stroke_width}"/>
        <text x="{padding}" y="{title_font_size}" font-size="{title_font_size}" font-weight="bold">{html.escape(title)}</text>
""")
        y = title_space
        for k, v in table.items():
            # Draw the color box with a black border
            f.write(f"""        <rect x="{padding}" y="{y}" width="{box_width}" height="{box_height - 10}" fill="{v}" stroke="black" stroke-width="{stroke_width}"/>
        <text x="{text_offset_x}" y="{y + text_offset_y}" font-size="{item_font_size}">{html.escape(k)}</text>
""")
            y += box_height
        f.write("""    </g>
</svg>
""")

def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, """names of color table files.
                        Files are tab separated with two columns and no header, columns are: annotation, hex color. For example: CCDB   cc0000""")
    add_output(parser, "svg", "Write the legend to this svg file.", primary=True, required=True)
    parser.add_argument("--title", default="Legend", type=str, required=False,
                        help="Title of the legend.")
    #TODO: font, font size, indicator shape, etc.

    parser.add_argument('--config', action=ActionConfigFile)
    declare_outputs(parser, ["svg_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)

    ### validate input


    check_outputs(parser, params)
    table = OrderedDict()
    for color_table_file in params.input:
        if color_table_file in STDIN_PATHS:
            color_table_file = sys.stdin
        for k, v in read_color_table(color_table_file).items():
            table[k] = v

    ### Run
    
    color_table_to_legend(table, svg_out=params.svg_out, title=params.title)
