"""Convert system calls made by other tools into Domain_System features

Reads GenBank records that another tool has annotated, and adds a Domain_System feature (see
find_systems) for each system it called, so that Domainator's selection, reporting, coloring and
plotting tools work on them.

Sources:
    antismash: antiSMASH GenBank output. --level protocluster (the default) makes one system per
        protocluster, located at its core, named by its product, in the antismash database; its
        members are the core CDSs whose profiles triggered the product. --level candidate_cluster
        and --level region make one system per candidate cluster (with antiSMASH's /kind) or
        region, located at the whole feature. Systems imported before are replaced.

    domainator import_systems -i antismash_output.gbk --source antismash -o systems.gb
    domainator select_by_system -i systems.gb --system_categories PKS --model_neighborhood -o pks_regions.gb
"""
from typing import Iterable, Iterator

from Bio.SeqRecord import SeqRecord
from jsonargparse import ActionConfigFile

from domainator.cli.common import add_input, add_seq_outputs, check_outputs, declare_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from domainator.seqio import parse_seqfiles
from domainator.systems.imports.antismash import LEVELS, import_antismash

SOURCES = ("antismash",)


def import_systems(records: Iterable[SeqRecord], *, source: str = "antismash", level: str = "protocluster",
                   clear_source: bool = False) -> Iterator[SeqRecord]:
    """Yield records with a Domain_System feature for each system that another tool called.

    Args:
        records: the records, as written by the source tool.
        source: the tool that made the calls: "antismash".
        level: antiSMASH only: which antiSMASH feature becomes a system: "protocluster" (at its
            core), "candidate_cluster" or "region".
        clear_source: remove the source tool's cluster features (antiSMASH's region,
            cand_cluster, protocluster and proto_core) after converting them.
    """
    if source not in SOURCES:
        raise ValueError(f"source must be one of {', '.join(SOURCES)}, not {source!r}")
    if level not in LEVELS:
        raise ValueError(f"level must be one of {', '.join(LEVELS)}, not {level!r}")
    for record in records:
        yield import_antismash(record, level=level, clear_source=clear_source)


def build_parser():
    parser = tool_parser(__doc__)
    add_input(parser, "GenBank files written by the source tool.")
    add_seq_outputs(parser)
    parser.add_argument("--source", type=str, default="antismash", choices=SOURCES,
                        help="The tool that made the system calls.")
    parser.add_argument("--level", type=str, default="protocluster", choices=LEVELS,
                        help="antiSMASH only: which antiSMASH feature becomes a system. protocluster: one per "
                             "protocluster, located at its core. candidate_cluster, region: one per candidate cluster "
                             "or region, located at the whole feature.")
    parser.add_argument("--clear_source", action="store_true", default=False,
                        help="Remove the source tool's cluster features (antiSMASH's region, cand_cluster, protocluster "
                             "and proto_core) after converting them.")
    parser.add_argument("--config", action=ActionConfigFile)
    declare_outputs(parser, ["genbank_out", "fasta_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params)
    write_seq_outputs(import_systems(parse_seqfiles(params.input), source=params.source, level=params.level,
                                     clear_source=params.clear_source),
                      genbank=params.genbank_out, fasta=params.fasta_out)
