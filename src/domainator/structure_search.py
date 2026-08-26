"""Search a structure database with a small number of reference structures.

The structure-input counterpart of domain_search.py. For searching large structure
databases (for example a foldseek database built from AlphaFold DB) with a handful of
reference structures.

Only records that had a hit are written. Each one carries a Domain_Search feature marking
the region homologous to the best-matching reference; pass --add_annotations to also write
a Domainator feature for every hit.

References are aligned as the query and the searched database as the target, which is both
what foldseek is optimized for and what makes --evalue mean the same thing it means
elsewhere in domainator: it is computed against the database being searched. The hit
coordinates are positions on the database record; the matching positions on the reference
are recorded in the /rstart, /rend and /rlen qualifiers.

Hit records are built from the sequences the backend reports, so the searched database is
never read end to end and its size does not bound memory. Note that --max_seqs bounds how
many hits the backend keeps per reference; the default is low for a large database, and a
warning is issued whenever a reference saturates it.

Use structure_domainate.py instead to annotate every input structure rather than select
hits, and structure_to_genbank.py to convert structures without searching at all.
"""

import heapq
import sys

from jsonargparse import ArgumentParser, ActionConfigFile

from domainator import __version__, RawAndDefaultsFormatter
from domainator import domainate
from domainator import structure_lib
from domainator.output_guardrails import add_max_output_gb_argument, max_output_gb_to_bytes


MITIGATION_OPTIONS = ["--max_hits", "-e/--evalue", "--no_annotations"]


def structure_search(input_values, reference_values, aligner, evalue=0.001, min_evalue=0.0,
                     max_hits=None, max_overlap=1.0, overlap_by_db=False,
                     add_annotations=False, alignment_type=2, metrics=None, max_seqs=None,
                     tmp_dir=None, keep_db=None, database_name=None, hits_tsv=None):
    """Yield protein SeqRecords for database records that matched a reference.

    Args:
        input_values: -i values: the database to search. Structure files, globs,
            directories, or a prebuilt database prefix.
        reference_values: -r values: the reference structures to search with.
        aligner: a structure_lib.StructureAligner.
        evalue: keep hits with E-value strictly below this.
        min_evalue: keep hits with E-value at or above this. Useful for excluding a
            reference's hit to its own copy in the database.
        max_hits: return at most this many records, the best-scoring ones. None means all.
        max_overlap: maximum fractional overlap between annotations on one record.
        overlap_by_db: filter overlaps within each database rather than across all.
        add_annotations: also write a Domainator feature per hit, not just the
            Domain_Search best-hit feature.
        alignment_type: backend alignment type.
        metrics: extra structural metrics to compute and store as qualifiers.
        max_seqs: maximum hits the backend keeps per reference.
        tmp_dir: parent directory for temporary databases and scratch files.
        keep_db: persist the database built from input structures at this prefix.
        database_name: value for the /database qualifier. Defaults to the reference name.
        hits_tsv: copy the backend's raw hit table here, for cross-checking.

    Yields:
        SeqRecord objects with molecule_type "protein". Sorted best-first when max_hits is
        set, otherwise streamed in whatever order the backend grouped them.
    """
    metrics = list(metrics or [])

    with structure_lib.prepared_databases(aligner, input_values, reference_values,
                                          tmp_dir=tmp_dir, keep_db=keep_db) as prepared:
        aligner.check_capabilities(alignment_type, metrics,
                                   [prepared.input_db, prepared.reference_db])
        db_name = database_name or prepared.reference_label

        hits = aligner.search(prepared.reference_db, prepared.input_db, evalue=evalue,
                              max_seqs=max_seqs, want_tseq=True,
                              alignment_type=alignment_type, metrics=metrics,
                              work_dir=prepared.work_dir, sort_by_target=True)

        hit_counts = {}

        def counted(hit_iterable):
            for hit in hit_iterable:
                hit_counts[hit.query] = hit_counts.get(hit.query, 0) + 1
                yield hit

        groups = structure_lib.group_hits_by_target(
            counted(hits), db_name, evalue, min_evalue, aligner.name)

        if aligner.provides_target_sequence:
            # Search output carries each hit target's full sequence, so records are built
            # without ever reading the searched database. group_hits_by_target hands each
            # group its own sequence, so --max_hits bounds how many are held at once.
            emitted = _emit_groups(groups, max_hits, max_overlap, overlap_by_db,
                                   add_annotations)
        else:
            # The backend reports no target sequence (reseek), so the hits are collected
            # first and the sequences are then resolved in one pass over the database.
            # Memory is bounded by the selected hits, not by the database.
            emitted = _emit_groups_with_db_sequences(
                groups, aligner, prepared.input_db, max_hits, max_overlap, overlap_by_db,
                add_annotations)

        for record in emitted:
            yield record

        structure_lib.warn_on_saturation(hit_counts, aligner.effective_max_seqs(max_seqs))
        if hits_tsv is not None:
            structure_lib.copy_hits_tsv(prepared.work_dir, hits_tsv, aligner)


def _build_record(target, record_hits, sequence, max_overlap, overlap_by_db, add_annotations):
    if not sequence:
        raise RuntimeError(
            f"The backend reported a hit to '{target}' without its sequence, so the "
            "output record cannot be built. This usually means the target database is "
            "missing entries it claims to contain."
        )
    record = structure_lib.build_protein_record(target, sequence)
    domainate.add_protein_annotations(
        record, record_hits, sys.maxsize, max_overlap,
        not add_annotations, True, overlap_by_db)
    return record


def _emit_groups(groups, max_hits, max_overlap, overlap_by_db, add_annotations):
    """Build records from grouped hits, keeping only the best max_hits of them.

    Without max_hits records stream out as their groups complete. With it, a bounded heap
    keyed on best hit score holds the current top max_hits and the survivors are written
    best-first once the search finishes -- the same policy as domain_search.py. Each heap
    entry carries its own sequence, so the heap bounds memory as well as output.
    """
    if max_hits is None:
        for target, record_hits, sequence in groups:
            yield _build_record(target, record_hits, sequence, max_overlap,
                                overlap_by_db, add_annotations)
        return

    heap = []
    counter = 0
    for target, record_hits, sequence in groups:
        best_score = max(hit.score for hit in record_hits)
        counter += 1
        # counter breaks score ties deterministically and stops the comparison before it
        # reaches the payload, which heapq would otherwise try to order.
        entry = (best_score, -counter, target, record_hits, sequence)
        if len(heap) < max_hits:
            heapq.heappush(heap, entry)
        elif entry > heap[0]:
            heapq.heapreplace(heap, entry)

    for best_score, _, target, record_hits, sequence in sorted(heap, reverse=True):
        yield _build_record(target, record_hits, sequence, max_overlap,
                            overlap_by_db, add_annotations)


def _emit_groups_with_db_sequences(groups, aligner, input_db, max_hits, max_overlap,
                                   overlap_by_db, add_annotations):
    """Emit records for backends whose search output has no target sequences.

    The hit groups are collected and reduced to the selected set first, then the database
    is streamed once to attach sequences. Records come out in database order rather than
    score order, which is why the score-ordered case is handled separately.
    """
    selected = {}
    heap = []
    counter = 0
    for target, record_hits, _ in groups:
        best_score = max(hit.score for hit in record_hits)
        counter += 1
        entry = (best_score, -counter, target)
        if max_hits is None:
            selected[target] = record_hits
            continue
        if len(heap) < max_hits:
            heapq.heappush(heap, entry)
            selected[target] = record_hits
        elif entry > heap[0]:
            _, _, evicted = heapq.heapreplace(heap, entry)
            selected.pop(evicted, None)
            selected[target] = record_hits

    if not selected:
        return

    if max_hits is None:
        order = None
    else:
        # emit best-first, matching the behavior of the sequence-carrying path
        order = {target: rank for rank, (_, _, target)
                 in enumerate(sorted(heap, reverse=True))}

    found = {}
    for name, sequence, _ in aligner.iter_sequences(input_db):
        if name in selected:
            found[name] = sequence
            if len(found) == len(selected):
                break

    missing = sorted(set(selected) - set(found))
    if missing:
        raise RuntimeError(
            f"{len(missing)} hit target(s) are absent from the searched database "
            f"(for example {missing[:3]}), so their records cannot be built."
        )

    targets = sorted(selected, key=lambda t: order[t]) if order is not None else list(found)
    for target in targets:
        yield _build_record(target, selected[target], found[target], max_overlap,
                            overlap_by_db, add_annotations)


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', default=None, required=True, nargs='+', type=str,
                        help="the structure database to search: the prefix of a prebuilt foldseek database, or pdb/cif files (optionally gzipped), a directory to search recursively, or a glob.")
    parser.add_argument('-r', '--references', default=None, required=True, nargs='+', type=str,
                        help="reference structures to search with. Same accepted forms as --input.")
    parser.add_argument('-o', '--output', default=None, type=str,
                        help="output genbank filename. If not supplied, or if given as '-', writes to stdout.")
    parser.add_argument('-e', '--evalue', type=float, default=0.001,
                        help="threshold E value for a hit. E values are computed against the searched database.")
    parser.add_argument('--min_evalue', type=float, default=0.0,
                        help="minimum E value for a hit. Useful for excluding a reference's match to its own copy in the database.")
    parser.add_argument('--max_hits', type=int, default=None,
                        help="return at most this many records, the best-scoring ones. Applies to the whole search, not per reference. When set, output is sorted best-first and written after the search completes.")
    parser.add_argument('--max_overlap', type=float, default=1,
                        help="maximum fraction of overlap allowed between annotations on one record. 1 disables overlap filtering.")
    parser.add_argument('--overlap_by_db', action='store_true', default=False,
                        help="filter overlapping annotations within each database independently, instead of across all databases.")
    parser.add_argument('--add_annotations', action='store_true', default=False,
                        help="also write a Domainator feature for every hit, in addition to the Domain_Search feature for the best hit.")
    parser.add_argument('--database_name', default=None, type=str,
                        help="value for the /database qualifier on the annotations. Defaults to the reference database or directory name.")
    parser.add_argument('--keep_db', default=None, type=str,
                        help="when the input is structure files rather than a prebuilt database, also write the database built from them, using this path as its prefix.")
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
    if params.max_hits is not None and params.max_hits < 1:
        raise RuntimeError("--max_hits must be at least 1.")

    aligner = structure_lib.build_aligner(params)
    # Validate the backend-capability arguments before any database is built, so an
    # unsupported request fails immediately rather than after the expensive step.
    aligner.effective_max_seqs(params.max_seqs)
    aligner.check_capabilities(params.alignment_type, params.metrics or [], [])

    max_output_bytes = max_output_gb_to_bytes(params.max_output_gb)
    output_description = f"structure_search genbank output ({params.output or 'stdout'})"

    records = structure_search(
        params.input, params.references, aligner,
        evalue=params.evalue, min_evalue=params.min_evalue, max_hits=params.max_hits,
        max_overlap=params.max_overlap, overlap_by_db=params.overlap_by_db,
        add_annotations=params.add_annotations, alignment_type=params.alignment_type,
        metrics=params.metrics, max_seqs=params.max_seqs, tmp_dir=params.tmp_dir,
        keep_db=params.keep_db, database_name=params.database_name,
        hits_tsv=params.hits_tsv,
    )

    structure_lib.write_records(records, params.output, max_output_bytes,
                                output_description, MITIGATION_OPTIONS)


def _entrypoint():
    main(sys.argv[1:])


if __name__ == '__main__':
    main(sys.argv[1:])
