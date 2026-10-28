"""The command functions are the library API (docs/stability.md): each command is
``from domainator.<command> import <command>``, takes its main input positionally, and takes
everything else as keywords named like the command's options. These tests keep that true.
"""
import importlib
import inspect
import typing

import pytest

from domainator.cli.registry import TOOLS

# Keyword parameters that are not one of the command's options. Everything else is named like
# an option (its dest), the plural of a repeatable one (references for -r/--reference), or
# <dest>_handle (an open handle for that output).
NOT_OPTIONS = {
    "domainate": {"z", "ncbi_taxonomy", "allowed_taxids", "pre_parsed_references", "best_annotation", "max_hits_per_contig"},
    "domain_search": {"z", "ncbi_taxonomy", "normalize_direction"},
    "structure_domainate": {"aligner", "keep_db"},
    "structure_search": {"aligner", "keep_db"},
    "clean_sequences": {"stats"},
    "color_genbank": {"color_feature"},
    "deduplicate_genbank": {"extra_args", "bin_paths", "add_count", "clusters"},
    "extract_peptides": {"cds_ids", "keep_cds_feature", "cds_id_type", "extract_all", "_from_domain_search", "_domain_search_negatives"},
    "rotate_contigs": {"z", "pre_parsed_references"},
    "select_by_cds": {"cds_ids", "normalize_direction", "_from_domain_search", "_domain_search_negatives"},
    "select_by_contig": {"ncbi_taxonomy"},
    "select_by_coord": {"coordinate_selections"},
    "structure_to_genbank": {"aligner"},
    "hmmer_select": {"hmm_iterator", "search_name", "search_description", "search_accession", "match_all"},
    "hmmer_search": {"max_output_bytes", "output_description"},
    "hmmer_compare": {"max_output_bytes", "output_description"},
    "structure_dist": {"aligner"},
    "compare_contigs": {"metrics"},
    "transform_matrix": {"row_lengths", "col_lengths"},
    "build_ssn": {"subset_labels", "max_output_bytes"},
    "build_tree": {"subset_labels"},
    "build_projection": {"algorithm_parameters"},
    "enum_report": {"analyses", "ncbi_taxonomy"},
    "hmmer_report": {"analyses"},
    "summary_report": {"ncbi_taxonomy"},
    "matrix_report": {"stage_timings", "progress_callback"},
    "plot_contigs": {"arrow_width", "contig_thickness", "contig_label_font", "contig_label_px", "feature_label_font",
                     "details_height", "details_width"},
    "db_download": {"output", "user_agent"},
}


def _command(name):
    module = importlib.import_module(f"domainator.{TOOLS[name].module}")
    return module, getattr(module, name)


@pytest.mark.parametrize("name", list(TOOLS))
def test_only_the_main_input_is_positional(name):
    _module, function = _command(name)
    parameters = list(inspect.signature(function).parameters.values())
    assert parameters[0].kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
    assert all(p.kind is inspect.Parameter.KEYWORD_ONLY for p in parameters[1:]), name


@pytest.mark.parametrize("name", list(TOOLS))
def test_keywords_are_named_like_options(name):
    module, function = _command(name)
    dests = {action.dest for action in module.build_parser()._actions}

    def named_like_an_option(keyword):
        return (keyword in dests or keyword.removesuffix("s") in dests
                or keyword.removesuffix("_handle") in dests)

    keywords = [p.name for p in inspect.signature(function).parameters.values() if p.kind is inspect.Parameter.KEYWORD_ONLY]
    others = {k for k in keywords if not named_like_an_option(k)}
    assert others == NOT_OPTIONS.get(name, set())


@pytest.mark.parametrize("name", list(TOOLS))
def test_documented_and_annotated(name):
    """The API reference (docs/api/commands.md) shows each function's summary line and types."""
    _module, function = _command(name)
    summary = (inspect.getdoc(function) or "").split("\n", 1)[0]
    assert summary and summary[0].isupper() and not summary.endswith(":"), (name, summary)
    signature = inspect.signature(function)
    assert signature.return_annotation is not inspect.Signature.empty
    assert all(p.annotation is not inspect.Parameter.empty for p in signature.parameters.values())
    import pyhmmer  # hmmer_select's annotations name it without importing it at run time

    hints = typing.get_type_hints(function, localns={"pyhmmer": pyhmmer})
    defaults_none = [p.name for p in signature.parameters.values() if p.default is None]
    assert [p for p in defaults_none if type(None) not in typing.get_args(hints[p])] == []  # Optional[...]


def test_package_root_exports():
    """The package root exports the shared pieces lazily (importing domainator stays cheap)."""
    import domainator

    for name in domainator.__all__:
        assert getattr(domainator, name) is not None
        assert name in dir(domainator)
    from domainator import DataMatrix, DomainatorError, parse_seqfiles, write_genbank  # noqa: F401
    from domainator.errors import InputError
    assert issubclass(InputError, domainator.DomainatorError)
    with pytest.raises(AttributeError, match="no attribute 'utils'"):
        domainator.__getattr__("utils")


# Each command function gives what its command writes. (command, options, keyword arguments)
_EDITOR_CASES = [
    ("select_by_contig", "simple_genpept.gb", ["--domains", "CAT"], {"domains": {"CAT"}}),
    ("select_by_cds", "pDONR201_multi_genemark_domainator.gb", ["--domains", "CAT"], {"domains": {"CAT"}}),
    ("filter_domains", "pDONR201_multi_genemark_domainator.gb", ["-e", "1e-30"], {"evalue": 1e-30}),
    ("extract_domains", "pDONR201_multi_genemark_domainator.gb", ["--domains", "CAT", "APH"], {"domains": {"CAT", "APH"}}),
    ("extract_unannotated", "simple_genpept.gb", [], {}),
    ("trim_contigs", "pDONR201_multi_genemark_domainator.gb", ["--cds_range", "2", "--domains", "CAT"],
     {"cds_range": (2, 2), "domains": {"CAT"}}),
    ("clean_sequences", "FeSOD_20.fasta", ["--clean_name"], {"clean_name": True}),
    ("sort_contigs", "pDONR201_multi_genemark_domainator.gb", ["--domains", "CAT"], {"domains": {"CAT"}}),
]


@pytest.mark.parametrize("name,fixture,options,keywords", _EDITOR_CASES, ids=[c[0] for c in _EDITOR_CASES])
def test_editor_function_matches_its_command(shared_datadir, tmp_path, name, fixture, options, keywords):
    from domainator.seqio import parse_seqfiles, write_genbank

    module, function = _command(name)
    path = str(shared_datadir / fixture)
    command_out, function_out = tmp_path / "command.gb", tmp_path / "function.gb"
    module.main(["-i", path, "-o", str(command_out), *options])
    write_genbank(function(parse_seqfiles([path]), **keywords), str(function_out))
    assert function_out.read_text() == command_out.read_text()
    assert command_out.stat().st_size > 0


def test_matrix_functions_match_their_commands(shared_datadir, tmp_path):
    import numpy as np

    from domainator.data_matrix import DataMatrix
    from domainator.kmer_dist import kmer_dist, main as kmer_dist_main
    from domainator.transform_matrix import main as transform_main, transform_matrix

    path = str(shared_datadir / "kmer_dist_proteins.faa")
    kmer_dist_main(["-i", path, "--hdf_out", str(tmp_path / "k.hdf5")])
    from_command = DataMatrix.from_file(tmp_path / "k.hdf5")
    from_function = kmer_dist(path)
    assert from_function.rows == from_command.rows and from_function.columns == from_command.columns
    np.testing.assert_allclose(from_function.toarray(), from_command.toarray())

    transform_main(["-i", str(tmp_path / "k.hdf5"), "--hdf_out", str(tmp_path / "t.hdf5"), "--mode", "score_dist"])
    np.testing.assert_allclose(transform_matrix(from_command.data, mode="score_dist"),
                               DataMatrix.from_file(tmp_path / "t.hdf5").toarray())
    np.testing.assert_array_equal(transform_matrix(from_command.data).toarray(), from_command.toarray())


def test_report_function_matches_its_command(shared_datadir, tmp_path):
    import io

    from domainator.enum_report import enum_report, main as enum_report_main
    from domainator.seqio import parse_seqfiles

    path = str(shared_datadir / "pDONR201_multi_genemark_domainator.gb")
    enum_report_main(["-i", path, "-o", str(tmp_path / "r.tsv"), "--columns", "domains", "length"])
    handle = io.StringIO()
    enum_report(parse_seqfiles([path]), analyses=["domains", "length"], tsv_out_handle=handle)
    assert handle.getvalue() == (tmp_path / "r.tsv").read_text()
