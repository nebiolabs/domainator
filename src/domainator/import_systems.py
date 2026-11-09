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
    defensefinder: DefenseFinder's genes tables (*_defense_finder_genes.tsv), given with --calls,
        for the records DefenseFinder was run on. Systems are named by their subtype, in the
        defensefinder database (as system_db_download --source defensefinder names its models).
        Members are matched to CDSs by protein name (protein_id, locus_tag, gene_id, gene or
        cds_id), or for DefenseFinder's own gene calls (<contig>_<n>) as the contig's n-th CDS.
    padloc: PADLOC's tables (*_padloc.csv), given with --calls. Systems are named by the model's
        yaml.name, in the padloc database. Members are matched to CDSs by coordinates, then by
        name. PADLOC's table has no scores and doesn't say which genes are core, so systems score
        0 and every member is core. CRISPR arrays and ncRNAs are not imported.
    Systems from the same database already on the records are replaced.

    domainator import_systems -i antismash_output.gbk --source antismash -o systems.gb
    domainator import_systems -i genome.gb --source defensefinder --calls genome_defense_finder_genes.tsv -o systems.gb
    domainator select_by_system -i systems.gb --system_categories PKS --model_neighborhood -o pks_regions.gb
"""
from typing import Iterable, Iterator, Optional, Sequence

from Bio.SeqRecord import SeqRecord
from jsonargparse import ActionConfigFile

from domainator.cli.common import add_input, add_seq_outputs, check_outputs, declare_outputs, tool_parser
from domainator.outputs import write_seq_outputs
from domainator.seqio import parse_seqfiles
from domainator.errors import InputError
from domainator.systems.imports.antismash import LEVELS, import_antismash
from domainator.systems.imports.calltables import (DEFENSEFINDER_DATABASE, PADLOC_DATABASE, import_calls,
                                                   read_defensefinder, read_padloc)

SOURCES = ("antismash", "defensefinder", "padloc")
CALL_TABLE_SOURCES = {"defensefinder": (read_defensefinder, DEFENSEFINDER_DATABASE),
                      "padloc": (read_padloc, PADLOC_DATABASE)}


def import_systems(records: Iterable[SeqRecord], *, source: str = "antismash", level: str = "protocluster",
                   clear_source: bool = False, calls: Optional[Sequence[str]] = None) -> Iterator[SeqRecord]:
    """Yield records with a Domain_System feature for each system that another tool called.

    Args:
        records: the records, as written by the source tool.
        source: the tool that made the calls: "antismash", "defensefinder" or "padloc".
        level: antiSMASH only: which antiSMASH feature becomes a system: "protocluster" (at its
            core), "candidate_cluster" or "region".
        clear_source: remove the source tool's cluster features (antiSMASH's region,
            cand_cluster, protocluster and proto_core) after converting them.
        calls: defensefinder and padloc: the tool's output tables (*_defense_finder_genes.tsv,
            *_padloc.csv), read before the first record is yielded.
    """
    if source not in SOURCES:
        raise ValueError(f"source must be one of {', '.join(SOURCES)}, not {source!r}")
    if level not in LEVELS:
        raise ValueError(f"level must be one of {', '.join(LEVELS)}, not {level!r}")
    if source in CALL_TABLE_SOURCES:
        if not calls:
            raise InputError(f"--source {source} needs --calls: the tables of systems {source} called.")
        read, database = CALL_TABLE_SOURCES[source]
        by_contig = dict()
        for call in read(calls):
            by_contig.setdefault(call.contig, []).append(call)
        for record in records:
            yield import_calls(record, by_contig.get(record.id, []) or by_contig.get(record.name, []), database)
        return
    if calls:
        raise InputError("--calls is for --source defensefinder and padloc; antiSMASH's calls are in its GenBank output.")
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
    parser.add_argument("--calls", type=str, nargs="+", default=None,
                        help="defensefinder and padloc: the tables of the systems they called: DefenseFinder's "
                             "*_defense_finder_genes.tsv, or PADLOC's *_padloc.csv.")
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
    if params.source in CALL_TABLE_SOURCES and not params.calls:
        parser.error(f"--source {params.source} needs --calls: the tables of the systems it called.")
    if params.source not in CALL_TABLE_SOURCES and params.calls:
        parser.error("--calls is for --source defensefinder and padloc.")
    write_seq_outputs(import_systems(parse_seqfiles(params.input), source=params.source, level=params.level,
                                     clear_source=params.clear_source, calls=params.calls),
                      genbank=params.genbank_out, fasta=params.fasta_out)
