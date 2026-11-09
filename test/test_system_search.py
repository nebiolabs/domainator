"""system_search: find systems in sequence databases and select the regions around them.

test/data/systems/padloc_db is a miniature padloc-db whose GajA and GajB profiles are copies of
CcdB and CcdA, so each pDONR201 contig's ccdB-ccdA pair is a gabija system.
"""
import pytest

from domainator import DOMAIN_FEATURE_NAME, SYSTEM_FEATURE_NAME, SYSTEM_SEARCH_BEST_HIT_NAME
from domainator import domainate as domainate_module
from domainator import system_search as system_search_module
from domainator.errors import InputError
from domainator.model import get_hit_best_score
from domainator.select_by_system import select_by_system
from domainator.seqio import parse_seqfiles
from domainator.system_db_download import system_db_download
from domainator.system_search import system_search
from domainator.systems.records import MEMBERSHIP_QUALIFIER


@pytest.fixture
def padloc_models(shared_datadir, tmp_path):
    output_dir = tmp_path / "padloc"
    system_db_download(str(shared_datadir / "systems" / "padloc_db"), source="padloc", output_dir=str(output_dir))
    return str(output_dir)


@pytest.fixture
def contigs(shared_datadir):
    return str(shared_datadir / "pDONR201_multi_genemark.gb")


def records(path):
    return list(parse_seqfiles([str(path)]))


def features(record, feature_type):
    return [f for f in record.features if f.type == feature_type]


def search(contigs, **kwargs):
    return list(system_search(parse_seqfiles([contigs]), z=1000, cpu=2, **kwargs))


def test_model_directory(padloc_models, contigs):
    regions = search(contigs, system_models=[padloc_models], cds_range=(1, 1))
    assert [region.id for region in regions] == [f"pDONR201_{n}_2265:1rc" for n in range(1, 5)]
    for region in regions:
        focus, = features(region, SYSTEM_SEARCH_BEST_HIT_NAME)
        assert (focus.qualifiers["name"], focus.qualifiers["database"]) == (["gabija"], ["padloc"])
        # Only the members' hits; the system itself is the System_Search feature.
        assert sorted(hit.qualifiers["name"][0] for hit in features(region, DOMAIN_FEATURE_NAME)) == ["GajA_00001", "GajB_00001"]
        assert not features(region, SYSTEM_FEATURE_NAME)
        assert get_hit_best_score(region) == float(focus.qualifiers["score"][0])


def test_same_systems_as_domainate_then_select_by_system(padloc_models, contigs, shared_datadir):
    """system_search is domainate --system_models followed by select_by_system."""
    annotated = domainate_module.domainate(parse_seqfiles([contigs]), references=None, z=1000, cpu=1,
                                           system_models=[padloc_models])
    expected = list(select_by_system(annotated, cds_range=(2, 1)))
    found = search(contigs, system_models=[padloc_models], cds_range=(2, 1), add_annotations=True)
    assert len(found) == len(expected) == 4
    for region, other in zip(found, expected):
        assert str(region.seq) == str(other.seq)
        focus, = features(region, SYSTEM_SEARCH_BEST_HIT_NAME)
        other_focus, = features(other, SYSTEM_SEARCH_BEST_HIT_NAME)
        assert focus.qualifiers == other_focus.qualifiers and focus.location == other_focus.location
        # add_annotations keeps the region's hits and systems.
        assert len(features(region, DOMAIN_FEATURE_NAME)) == len(features(other, DOMAIN_FEATURE_NAME))
        assert len(features(region, SYSTEM_FEATURE_NAME)) == 1


def test_command_line_model(contigs, shared_datadir):
    hmms = str(shared_datadir / "pdonr_hmms.hmm")
    regions = search(contigs, references=[hmms], domain_expr="CcdB & CcdA", window_cds=2)
    assert len(regions) == 4
    focus, = features(regions[0], SYSTEM_SEARCH_BEST_HIT_NAME)
    assert (focus.qualifiers["name"], focus.qualifiers["database"]) == (["domain_expr"], ["command_line"])
    # The other references' hits (CAT, APH, ...) are dropped without add_annotations.
    assert sorted(hit.qualifiers["name"][0] for hit in features(regions[0], DOMAIN_FEATURE_NAME)) == ["CcdA", "CcdB"]
    annotated = search(contigs, references=[hmms], domain_expr="CcdB & CcdA", window_cds=2, add_annotations=True,
                       whole_contig=True)
    assert len({hit.qualifiers["name"][0] for hit in features(annotated[0], DOMAIN_FEATURE_NAME)}) > 2


def test_system_expr(padloc_models, contigs):
    assert len(search(contigs, system_models=[padloc_models], system_expr="gabija & ~cas_test")) == 4
    assert search(contigs, system_models=[padloc_models], system_expr="gabija & cas_test") == []
    with pytest.raises(InputError, match="names no system positively"):
        search(contigs, system_models=[padloc_models], system_expr="~cas_test")
    with pytest.raises(InputError, match="aren't in the system models"):
        search(contigs, system_models=[padloc_models], system_expr="gabija & ~no_such_system")


def test_max_hits_translate_and_min_context(padloc_models, contigs):
    assert len(search(contigs, system_models=[padloc_models], max_hits=2)) == 2
    proteins = search(contigs, system_models=[padloc_models], translate=True, max_hits=1)
    assert len(proteins) == 2  # one system's two members
    assert all(p.annotations["molecule_type"] == "protein" for p in proteins)
    memberships = sorted(f.qualifiers[MEMBERSHIP_QUALIFIER][0] for p in proteins for f in p.features if f.type == "CDS")
    assert memberships == ["core|GajA|padloc|gabija", "core|GajB|padloc|gabija"]
    # Each contig is 4470 bp, so no system has 5 kb on both sides.
    assert search(contigs, system_models=[padloc_models], min_context_kb=(5, 5)) == []
    assert len(search(contigs, system_models=[padloc_models], min_context_kb=(0.5, 0.5))) == 4


def test_main(padloc_models, contigs, tmp_path):
    out = tmp_path / "out.gb"
    system_search_module.main(["-i", contigs, "--system_models", padloc_models, "--systems", "gabija",
                               "--cds_range", "1", "--pad", "-o", str(out), "--cpu", "2"])
    regions = records(out)
    assert len(regions) == 4 and len({len(region) for region in regions}) == 1


def test_main_errors(padloc_models, contigs, shared_datadir, tmp_path, capsys):
    out = str(tmp_path / "out.gb")
    with pytest.raises(SystemExit):
        system_search_module.main(["-i", contigs, "-o", out])
    assert "--system_models, --domain_expr, or both" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        system_search_module.main(["-i", str(shared_datadir / "pDONR201_multi_genemark.gb"), "--system_models",
                                   padloc_models, "--translate", "--cds_range", "2", "-o", out])
    assert "--translate outputs proteins" in capsys.readouterr().err
