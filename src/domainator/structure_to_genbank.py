"""Convert protein structures into a Domainator protein genbank file.

Reads structure files (pdb/cif, optionally gzipped, given as files, globs, or a directory
searched recursively) or a prebuilt foldseek database, and writes one protein genbank
record per structure chain. This makes structures ordinary domainator inputs, so the
sequence-based tools (domainate.py, domain_search.py, seq_dist.py, build_ssn.py, the
reports) work on them with no further conversion.

Chain splitting, structure parsing and naming are all done by the structural aligner, so
records are named exactly as they are in a search against the same structures.

Use --keep_db to also persist the aligner database that was built along the way. That
database carries C-alpha coordinates, so it can be reused as the -i or -r argument of
structure_domainate.py / structure_search.py without paying to rebuild it, and it is how
you turn a directory of reference structures into a reusable reference database.

Note that record coordinates are indices into the *observed* residue sequence extracted
from the structure, not original PDB residue numbers: a chain with a disordered loop has a
sequence shorter than its residue-number span.
"""

import os
import sys
import tempfile

from jsonargparse import ArgumentParser, ActionConfigFile

from domainator import __version__, RawAndDefaultsFormatter
from domainator import structure_lib
from domainator.utils import write_genbank


def structure_to_genbank(inputs, aligner, keep_db=None, store_3di=False, tmp_dir=None):
    """Yield one protein SeqRecord per structure chain.

    Args:
        inputs: list of -i values (structure files, globs, directories, or one database
            prefix).
        aligner: a structure_lib.StructureAligner.
        keep_db: if not None, build the database at this prefix and leave it in place.
        store_3di: store each chain's 3Di string in a /threedi qualifier.
        tmp_dir: parent directory for the temporary database, when keep_db is None.

    Yields:
        SeqRecord objects with molecule_type "protein".
    """
    kind, resolved = structure_lib.resolve_structure_inputs(inputs)

    if kind == "db":
        if keep_db is not None:
            raise RuntimeError(
                "--keep_db builds a database from structure files; the input is already "
                f"a database ('{resolved[0]}')."
            )
        yield from _records(aligner, structure_lib.StructureDB(prefix=resolved[0]), store_3di)
        return

    if keep_db is not None:
        db = aligner.build_db(resolved, keep_db)
        yield from _records(aligner, db, store_3di)
        return

    with tempfile.TemporaryDirectory(dir=tmp_dir) as scratch:
        db = aligner.build_db(resolved, os.path.join(scratch, "structuredb"))
        yield from _records(aligner, db, store_3di)


def _records(aligner, db, store_3di):
    threedi = dict(aligner.iter_3di(db)) if store_3di else {}
    for name, sequence, source_path in aligner.iter_sequences(db):
        yield structure_lib.build_protein_record(
            name, sequence, source_path, threedi.get(name) if store_3di else None
        )


def main(argv):
    parser = ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)

    parser.add_argument('-i', '--input', default=None, required=True, nargs='+', type=str,
                        help="structure files (pdb/cif, optionally gzipped), a directory to search recursively, a glob, or the prefix of a prebuilt foldseek database.")
    parser.add_argument('-o', '--output', default=None, type=str,
                        help="output genbank filename. If not supplied, or if given as '-', writes to stdout. Optional when --keep_db is given.")
    parser.add_argument('--keep_db', default=None, type=str,
                        help="also write the aligner database built from the input structures, using this path as its prefix. The database keeps C-alpha coordinates, so it can be reused as an input or reference for structure_domainate.py and structure_search.py.")
    parser.add_argument('--store_3di', action='store_true', default=False,
                        help="store each chain's 3Di structural alphabet string in a /threedi qualifier on a full-length misc_feature. Note that this qualifier is a single opaque string, so it does not survive tools that slice records (extract_domains.py, select_by_coord.py).")
    structure_lib.add_backend_arguments(parser, include_search_arguments=False)
    parser.add_argument('--config', action=ActionConfigFile)

    params = parser.parse_args(argv)

    aligner = structure_lib.build_aligner(params)

    if params.store_3di and not aligner.supports_3di:
        raise RuntimeError(
            f"--store_3di is not available with the '{aligner.name}' backend, which does "
            "not expose a per-residue structural alphabet. Use --algorithm foldseek, or "
            "drop --store_3di."
        )

    # None or "-" means stdout, matching the rest of the suite.
    out = sys.stdout if params.output in (None, "-") else params.output

    write_genbank(
        structure_to_genbank(
            params.input,
            aligner,
            keep_db=params.keep_db,
            store_3di=params.store_3di,
            tmp_dir=params.tmp_dir,
        ),
        out,
        default_molecule_type="protein",
    )


def _entrypoint():
    main(sys.argv[1:])


if __name__ == '__main__':
    main(sys.argv[1:])
