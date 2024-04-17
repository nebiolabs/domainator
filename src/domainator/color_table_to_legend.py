""" Generates an SVG legend based on a color table.
The color table is a tsv file that maps annotations to hex color codes.
The generated SVG file contains rectangles filled with the corresponding colors and text labels for each annotation.
"""


from jsonargparse import ArgumentParser, ActionConfigFile
import sys
from typing import Tuple, List, Optional, Set, NamedTuple, Union, Dict
from os import PathLike
from domainator import __version__, RawAndDefaultsFormatter
from domainator.color_genbank import read_color_table
from collections import OrderedDict
import html


def color_table_to_legend(table: Dict[str, str], svg: str, title: str):
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
    

    with open(svg, "w") as f:
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

def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument("-i", "--input", default=None, nargs="+", type=str, required=False,
                        help="""names of color table files. If not supplied, reads from stdin.\n
                        Files are tab separated with two columns and no header, columns are: annotation, hex color. For example: CCDB   cc0000""")
    parser.add_argument("--svg", default=None, type=str, required=True,
                        help="name of output svg file")
    parser.add_argument("--title", default="Legend", type=str, required=False,
                        help="Title of the legend. Default: 'Legend'")
    #TODO: font, font size, indicator shape, etc.

    parser.add_argument('--config', action=ActionConfigFile)
    
    params = parser.parse_args(argv)

    ### validate input


    if params.input is None:
        input_path = (sys.stdin,)
    else:
        input_path = params.input
    
    table = OrderedDict()
    for color_table_file in input_path:
        for k, v in read_color_table(color_table_file).items():
            table[k] = v

    ### Run
    
    color_table_to_legend(table, params.svg, params.title)

def _entrypoint():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])


