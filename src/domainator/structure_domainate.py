"""Annotate protein structures with hits to reference structures.

The structure-input counterpart of domainator domainate. Reads structure files (pdb/cif,
optionally gzipped, given as files, globs, or a directory searched recursively) or a
prebuilt foldseek database, aligns reference structures against them, and writes a protein
genbank file where every input chain is one record carrying Domainator features for the
references that hit it.

Every input record is annotated and written, whether or not it had hits.
Use structure_search to select only matching records.
Use domainator structure_search instead when searching a large database with a few references.

References are aligned as the query and the input database as the target, so --evalue is
always computed against the database being annotated, exactly as it is for the rest of
domainator. The annotated coordinates are positions on the input record; the matching
positions on the reference are recorded in the /rstart, /rend and /rlen qualifiers.

Because the input database is built from real structures it carries C-alpha coordinates,
so --alignment_type 1 (TM-align) and the coordinate-based --metrics are available here.

Record coordinates are indices into the *observed* residue sequence extracted from the
structure, not original PDB residue numbers.
"""

from typing import Iterator, List, Optional
import sys

from jsonargparse import ActionConfigFile
from jsonargparse.typing import NonNegativeFloat, PositiveFloat

from domainator.cli.common import add_input, add_seq_outputs, check_outputs, declare_outputs, tool_parser
from domainator import domainate
from domainator import structure_lib
from domainator.output_guardrails import add_max_output_gb_argument, max_output_gb_to_bytes
from Bio.SeqRecord import SeqRecord


MITIGATION_OPTIONS = ["--max_hits_per_cds", "-e/--evalue"]


def structure_domainate(input_values: List[str], *, references: List[str], aligner: structure_lib.StructureAligner,
                        evalue: float = 0.001, min_evalue: float = 0.0,
                        max_hits_per_cds: Optional[int] = None, max_overlap: float = 1.0, overlap_by_db: bool = False,
                        alignment_type: Optional[int] = None,
                        metrics: Optional[List[str]] = None, max_seqs: Optional[int] = None, tmp_dir: Optional[str] = None,
                        keep_db: Optional[str] = None, database_name: Optional[str] = None,
                        hits_tsv_out: Optional[str] = None) -> Iterator[SeqRecord]:
    """Yield annotated protein SeqRecords, one per input structure chain.

    Args:
        input_values: -i values: structure files, globs, directories, or a database prefix.
        references: -r values, same forms as input_values.
        aligner: a structure_lib.StructureAligner.
        evalue: keep hits with E-value strictly below this.
        min_evalue: keep hits with E-value at or above this.
        max_hits_per_cds: maximum annotations to add per record (None: all).
        max_overlap: maximum fractional overlap between kept annotations. >= 1 disables
            overlap filtering.
        overlap_by_db: filter overlaps within each database rather than across all.
        alignment_type: the foldseek alignment type (None: the backend's default).
        metrics: extra structural metrics to compute and store as qualifiers.
        max_seqs: maximum hits the backend keeps per reference.
        tmp_dir: parent directory for temporary databases and scratch files.
        keep_db: persist the database built from the input structures at this prefix.
        database_name: value for the /database qualifier. Defaults to the reference name.
        hits_tsv_out: copy the backend's raw hit table here, for cross-checking.

    Yields:
        SeqRecord objects with molecule_type "protein".
    """
    metrics = list(metrics or [])

    with structure_lib.prepared_databases(aligner, input_values, references,
                                          tmp_dir=tmp_dir, keep_db=keep_db) as prepared:
        input_db = prepared.input_db
        reference_db = prepared.reference_db
        work_dir = prepared.work_dir
        aligner.check_capabilities(alignment_type, metrics, [input_db, reference_db])

        db_name = database_name or prepared.reference_label

        hits = aligner.search(reference_db, input_db, evalue=evalue, max_seqs=max_seqs,
                              want_tseq=False, alignment_type=alignment_type,
                              metrics=metrics, work_dir=work_dir)
        # Annotate mode has to walk every input record anyway, so hits are held in a dict
        # keyed by input record. The dict is bounded by (references x --max_seqs).
        hit_counts = {}
        materialized = []
        for hit in hits:
            hit_counts[hit.query] = hit_counts.get(hit.query, 0) + 1
            materialized.append(hit)
        structure_lib.warn_on_saturation(hit_counts, aligner.effective_max_seqs(max_seqs))
        if hits_tsv_out is not None:
            structure_lib.copy_hits_tsv(work_dir, hits_tsv_out, aligner)

        by_target = structure_lib.structure_hits_to_search_results(
            materialized, db_name, evalue, min_evalue, aligner.name)

        for name, sequence, source_path in aligner.iter_sequences(input_db):
            record_hits = by_target.get(name)
            record = structure_lib.build_protein_record(name, sequence, source_path)
            if record_hits:
                domainate.add_protein_annotations(
                    record, record_hits, sys.maxsize if max_hits_per_cds is None else max_hits_per_cds, max_overlap,
                    False, False, overlap_by_db)
            yield record


def build_parser():
    parser = tool_parser(__doc__)

    add_input(parser, "structures to annotate: pdb/cif files (optionally gzipped), a directory to search recursively, a glob, or the prefix of a prebuilt foldseek database.", stdin=False)
    parser.add_argument('-r', '--reference', default=None, required=True, nargs='+', type=str,
                        help="reference structures to annotate the input with. Same accepted forms as --input.")
    add_seq_outputs(parser)
    parser.add_argument('-e', '--evalue', type=PositiveFloat, default=0.001,
                        help="threshold E value for a hit. E values are computed against the input database.")
    parser.add_argument('--min_evalue', type=NonNegativeFloat, default=0.0,
                        help="minimum E value for a hit. Useful for excluding near-identical matches, such as a structure's hit to itself.")
    parser.add_argument('--max_hits_per_cds', type=int, default=None,
                        help="maximum number of annotations to add per record. If not supplied, there is no limit.")
    parser.add_argument('--max_overlap', type=float, default=1,
                        help="maximum fraction of overlap allowed between annotations on one record. 1 disables overlap filtering.")
    parser.add_argument('--overlap_by_db', action='store_true', default=False,
                        help="filter overlapping annotations within each database independently, instead of across all databases.")
    parser.add_argument('--database_name', default=None, type=str,
                        help="value for the /database qualifier on the annotations. Defaults to the reference database or directory name.")
    parser.add_argument('--db_out', default=None, type=str,
                        help="also write the aligner database built from the input structures, using this path as its prefix, so it can be reused as --input later.")
    parser.add_argument('--hits_tsv_out', default=None, type=str,
                        help="write the aligner's own hit table here, for debugging or for cross-checking against a hand-run search. With reseek this is the table before the per-reference --max_seqs cap is applied, so it may hold more rows than the output.")
    structure_lib.add_backend_arguments(parser)
    add_max_output_gb_argument(parser)
    parser.add_argument('--config', action=ActionConfigFile)
    declare_outputs(parser, ["genbank_out", "fasta_out"])
    return parser


def main(argv):
    parser = build_parser()
    params = parser.parse_args(argv)
    check_outputs(parser, params)

    if params.max_hits_per_cds is not None and params.max_hits_per_cds < 1:
        parser.error("--max_hits_per_cds must be at least 1. Leave it out to add every annotation.")

    aligner = structure_lib.build_aligner(params)
    # Validate the backend-capability arguments before any database is built, so an
    # unsupported request fails immediately rather than after the expensive step.
    aligner.effective_max_seqs(params.max_seqs)
    aligner.check_capabilities(params.alignment_type, params.metrics or [], [])

    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)

    records = structure_domainate(
        params.input, references=params.reference, aligner=aligner,
        evalue=params.evalue, min_evalue=params.min_evalue, max_hits_per_cds=params.max_hits_per_cds,
        max_overlap=params.max_overlap, overlap_by_db=params.overlap_by_db,
        alignment_type=params.alignment_type, metrics=params.metrics,
        max_seqs=params.max_seqs, tmp_dir=params.tmp_dir, keep_db=params.db_out,
        database_name=params.database_name, hits_tsv_out=params.hits_tsv_out,
    )

    structure_lib.write_records(records, params.genbank_out, params.fasta_out, max_output_bytes, MITIGATION_OPTIONS)
