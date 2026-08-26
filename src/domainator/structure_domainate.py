"""Annotate protein structures with hits to reference structures.

The structure-input counterpart of domainate.py. Reads structure files (pdb/cif,
optionally gzipped, given as files, globs, or a directory searched recursively) or a
prebuilt foldseek database, aligns reference structures against them, and writes a protein
genbank file where every input chain is one record carrying Domainator features for the
references that hit it.

Every input record is written, whether or not it had hits, unless --hits_only is given.
Use structure_search.py instead when searching a large database with a few references.

References are aligned as the query and the input database as the target, so --evalue is
always computed against the database being annotated, exactly as it is for the rest of
domainator. The annotated coordinates are positions on the input record; the matching
positions on the reference are recorded in the /rstart, /rend and /rlen qualifiers.

Because the input database is built from real structures it carries C-alpha coordinates,
so --alignment_type 1 (TM-align) and the coordinate-based --metrics are available here.

Record coordinates are indices into the *observed* residue sequence extracted from the
structure, not original PDB residue numbers.
"""

import sys

from jsonargparse import ArgumentParser, ActionConfigFile

from domainator import __version__, RawAndDefaultsFormatter
from domainator import domainate
from domainator import structure_lib
from domainator.output_guardrails import add_max_output_gb_argument, max_output_gb_to_bytes


MITIGATION_OPTIONS = ["--hits_only", "--max_domains", "-e/--evalue", "--no_annotations"]


def structure_domainate(input_values, reference_values, aligner, evalue=0.001, min_evalue=0.0,
                        max_domains=sys.maxsize, max_overlap=1.0, overlap_by_db=False,
                        hits_only=False, no_annotations=False, alignment_type=2, metrics=None,
                        max_seqs=None, tmp_dir=None, keep_db=None, database_name=None,
                        hits_tsv=None):
    """Yield annotated protein SeqRecords, one per input structure chain.

    Args:
        input_values: -i values: structure files, globs, directories, or a database prefix.
        reference_values: -r values, same forms as input_values.
        aligner: a structure_lib.StructureAligner.
        evalue: keep hits with E-value strictly below this.
        min_evalue: keep hits with E-value at or above this.
        max_domains: maximum annotations to add per record.
        max_overlap: maximum fractional overlap between kept annotations. >= 1 disables
            overlap filtering.
        overlap_by_db: filter overlaps within each database rather than across all.
        hits_only: write only records that had at least one hit.
        no_annotations: do not write Domainator features (only meaningful with hits_only).
        alignment_type: backend alignment type.
        metrics: extra structural metrics to compute and store as qualifiers.
        max_seqs: maximum hits the backend keeps per reference.
        tmp_dir: parent directory for temporary databases and scratch files.
        keep_db: persist the database built from the input structures at this prefix.
        database_name: value for the /database qualifier. Defaults to the reference name.
        hits_tsv: copy the backend's raw hit table here, for cross-checking.

    Yields:
        SeqRecord objects with molecule_type "protein".
    """
    metrics = list(metrics or [])

    with structure_lib.prepared_databases(aligner, input_values, reference_values,
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
        if hits_tsv is not None:
            structure_lib.copy_hits_tsv(work_dir, hits_tsv, aligner)

        by_target = structure_lib.structure_hits_to_search_results(
            materialized, db_name, evalue, min_evalue, aligner.name)

        for name, sequence, source_path in aligner.iter_sequences(input_db):
            record_hits = by_target.get(name)
            if record_hits is None and hits_only:
                continue
            record = structure_lib.build_protein_record(name, sequence, source_path)
            if record_hits:
                domainate.add_protein_annotations(
                    record, record_hits, max_domains, max_overlap,
                    no_annotations, False, overlap_by_db)
            yield record


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', default=None, required=True, nargs='+', type=str,
                        help="structures to annotate: pdb/cif files (optionally gzipped), a directory to search recursively, a glob, or the prefix of a prebuilt foldseek database.")
    parser.add_argument('-r', '--references', default=None, required=True, nargs='+', type=str,
                        help="reference structures to annotate the input with. Same accepted forms as --input.")
    parser.add_argument('-o', '--output', default=None, type=str,
                        help="output genbank filename. If not supplied, or if given as '-', writes to stdout.")
    parser.add_argument('-e', '--evalue', type=float, default=0.001,
                        help="threshold E value for a hit. E values are computed against the input database.")
    parser.add_argument('--min_evalue', type=float, default=0.0,
                        help="minimum E value for a hit. Useful for excluding near-identical matches, such as a structure's hit to itself.")
    parser.add_argument('--max_domains', type=int, default=0,
                        help="maximum number of annotations to add per record. 0 means no limit.")
    parser.add_argument('--max_overlap', type=float, default=1,
                        help="maximum fraction of overlap allowed between annotations on one record. 1 disables overlap filtering.")
    parser.add_argument('--overlap_by_db', action='store_true', default=False,
                        help="filter overlapping annotations within each database independently, instead of across all databases.")
    parser.add_argument('--hits_only', action='store_true', default=False,
                        help="write only records that had at least one hit.")
    parser.add_argument('--no_annotations', action='store_true', default=False,
                        help="do not write Domainator features. Only meaningful together with --hits_only, where it selects hit records without annotating them.")
    parser.add_argument('--database_name', default=None, type=str,
                        help="value for the /database qualifier on the annotations. Defaults to the reference database or directory name.")
    parser.add_argument('--keep_db', default=None, type=str,
                        help="also write the aligner database built from the input structures, using this path as its prefix, so it can be reused as --input later.")
    parser.add_argument('--hits_tsv', default=None, type=str,
                        help="write the aligner's own hit table here, for debugging or for cross-checking against a hand-run search. With reseek this is the table before the per-reference --max_seqs cap is applied, so it may hold more rows than the output.")
    structure_lib.add_backend_arguments(parser)
    add_max_output_gb_argument(parser)
    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)

    if params.evalue <= 0:
        raise RuntimeError("--evalue must be greater than 0.")
    if params.min_evalue < 0:
        raise RuntimeError("--min_evalue must be greater than or equal to 0.")
    if params.no_annotations and not params.hits_only:
        raise RuntimeError("--no_annotations only makes sense together with --hits_only, "
                           "otherwise the output would carry no information at all.")
    max_domains = params.max_domains if params.max_domains > 0 else sys.maxsize

    aligner = structure_lib.build_aligner(params)
    # Validate the backend-capability arguments before any database is built, so an
    # unsupported request fails immediately rather than after the expensive step.
    aligner.effective_max_seqs(params.max_seqs)
    aligner.check_capabilities(params.alignment_type, params.metrics or [], [])

    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)
    output_description = f"structure_domainate genbank output ({params.output or 'stdout'})"

    records = structure_domainate(
        params.input, params.references, aligner,
        evalue=params.evalue, min_evalue=params.min_evalue, max_domains=max_domains,
        max_overlap=params.max_overlap, overlap_by_db=params.overlap_by_db,
        hits_only=params.hits_only, no_annotations=params.no_annotations,
        alignment_type=params.alignment_type, metrics=params.metrics,
        max_seqs=params.max_seqs, tmp_dir=params.tmp_dir, keep_db=params.keep_db,
        database_name=params.database_name, hits_tsv=params.hits_tsv,
    )

    structure_lib.write_records(records, params.output, max_output_bytes,
                                output_description, MITIGATION_OPTIONS)


def _entrypoint():
    main(sys.argv[1:])


if __name__ == '__main__':
    main(sys.argv[1:])
