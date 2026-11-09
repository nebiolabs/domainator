"""import_systems: antiSMASH cluster calls as Domain_System features.

test/data/systems/antismash_region.gbk is antiSMASH 8.0.4's example output for region 1 of
Streptomyces coelicolor A3(2) (NC_003888.3), trimmed to the features import_systems reads. It has
two protoclusters, hglE-KS (core: SCO0126-SCO0127) and T1PKS (core: SCO0127), in one
chemical-hybrid candidate cluster.
"""
import pytest

from domainator import SYSTEM_FEATURE_NAME
from domainator import import_systems as import_systems_module
from domainator.import_systems import import_systems
from domainator.select_by_system import select_by_system
from domainator.seqio import parse_seqfiles


@pytest.fixture
def antismash_record(shared_datadir):
    return next(iter(parse_seqfiles([str(shared_datadir / "systems" / "antismash_region.gbk")])))


def systems(record):
    return [f for f in record.features if f.type == SYSTEM_FEATURE_NAME]


def test_protoclusters(antismash_record):
    record, = import_systems([antismash_record])
    hgl, t1pks = systems(record)
    assert (hgl.qualifiers["name"], hgl.qualifiers["database"], hgl.qualifiers["category"]) == (["hglE-KS"], ["antismash"], ["PKS"])
    assert str(hgl.location) == "[20000:33018](+)"  # the core
    assert hgl.qualifiers["component"] == ["core|hglE|20001_1_26249|antismash|hglE", "core|hglD|26296_1_33018|antismash|hglD"]
    assert hgl.qualifiers["score"] == ["1334.0"]  # 720.9 + 613.1, each CDS's triggering profile
    assert hgl.qualifiers["evalue"] == ["7.4e-185"]
    assert hgl.qualifiers["neighborhood_kb"] == ["20"] and hgl.qualifiers["window"] == ["anchor:kb=20"]
    assert hgl.qualifiers["description"] == ["(hglE or hglD)"]
    # One component per profile that triggered the product: a CDS can fill several.
    assert t1pks.qualifiers["component"] == ["core|PKS_KS|26296_1_33018|antismash|PKS_KS",
                                             "core|PKS_AT|26296_1_33018|antismash|PKS_AT"]
    assert t1pks.qualifiers["member_count"] == ["1"] and t1pks.qualifiers["score"] == ["91.4"]
    # The members' cds_ids are stored on their CDSs.
    cds_ids = {f.qualifiers["cds_id"][0] for f in record.features if f.type == "CDS"}
    assert {"20001_1_26249", "26296_1_33018"} <= cds_ids


@pytest.mark.parametrize("level, kind", [("candidate_cluster", ["chemical_hybrid"]), ("region", None)])
def test_candidate_clusters_and_regions(antismash_record, level, kind):
    record, = import_systems([antismash_record], level=level)
    system, = systems(record)
    assert system.qualifiers["name"] == ["hglE-KS+T1PKS"]
    assert str(system.location) == "[0:53018](+)"
    assert system.qualifiers.get("kind") == kind
    assert [value.split("|")[1] for value in system.qualifiers["component"]] == ["hglE", "hglD", "PKS_KS", "PKS_AT"]


def test_contig_edge(antismash_record):
    for feature in antismash_record.features:
        if feature.type == "protocluster":
            feature.qualifiers["contig_edge"] = ["True"]
    record, = import_systems([antismash_record], level="candidate_cluster")
    assert systems(record)[0].qualifiers["system_partial"] == ["<>"]
    record, = import_systems([record])
    assert all("system_partial" not in system.qualifiers for system in systems(record))  # the cores don't touch an end


def test_reimporting_replaces_and_clear_source(antismash_record):
    record, = import_systems([antismash_record])
    record, = import_systems([record], clear_source=True)
    assert len(systems(record)) == 2
    assert not {"region", "cand_cluster", "protocluster", "proto_core"} & {f.type for f in record.features}


def test_imported_systems_select_with_their_neighborhood(antismash_record):
    record, = import_systems([antismash_record])
    regions = list(select_by_system([record], systems={"T1PKS"}, model_neighborhood=True))
    assert len(regions) == 1 and len(regions[0]) == 53018 - 6295  # the core padded by 20 kb, cut at the contig end


def test_main(shared_datadir, tmp_path):
    out = tmp_path / "out.gb"
    import_systems_module.main(["-i", str(shared_datadir / "systems" / "antismash_region.gbk"), "-o", str(out),
                                "--level", "region"])
    record, = parse_seqfiles([str(out)])
    assert [f.qualifiers["name"][0] for f in systems(record)] == ["hglE-KS+T1PKS"]


# DefenseFinder and PADLOC call tables (synthetic, for pDONR201_multi_genemark.gb; its CDSs are
# named pDONR201_1 ... pDONR201_6 by gene_id, with ccdB second and ccdA third).

@pytest.fixture
def pdonr(shared_datadir):
    return str(shared_datadir / "pDONR201_multi_genemark.gb")


def test_defensefinder_calls(pdonr, shared_datadir):
    calls = [str(shared_datadir / "systems" / "defensefinder_genes.tsv")]
    with pytest.warns(UserWarning, match="pDONR201_2: 1 system members match no CDS"):
        records = list(import_systems(parse_seqfiles([pdonr]), source="defensefinder", calls=calls))
    system, = systems(records[0])
    assert (system.qualifiers["name"], system.qualifiers["database"], system.qualifiers["program"]) == (
        ["CcdAB"], ["defensefinder"], ["macsyfinder"])
    assert (system.qualifiers["accession"], system.qualifiers["category"]) == (["DefenseFinder/CcdAB/CcdAB"], ["DefenseFinder"])
    assert str(system.location) == "[958:1391](-)"
    # ccdB by its gene_id, ccdA as DefenseFinder's gene 3 of the replicon.
    assert system.qualifiers["component"] == ["core|CcdB|1264_-1_959|defensefinder|CcdAB__CcdB",
                                              "accessory|CcdA|1391_-1_1266|defensefinder|CcdAB__CcdA_2"]
    assert (system.qualifiers["score"], system.qualifiers["evalue"], system.qualifiers["completeness"]) == (
        ["2.0"], ["3.4e-15"], ["1.000"])
    assert [f.qualifiers["member_count"] for f in systems(records[1])] == [["1"]]
    assert systems(records[2]) == systems(records[3]) == []


def test_padloc_calls(pdonr, shared_datadir):
    calls = [str(shared_datadir / "systems" / "padloc.csv")]
    with pytest.warns(UserWarning, match="1 system members match no CDS"):  # the CRISPR array
        record = next(iter(import_systems(parse_seqfiles([pdonr]), source="padloc", calls=calls)))
    system, = systems(record)
    assert (system.qualifiers["name"], system.qualifiers["database"], system.qualifiers["category"]) == (
        ["gabija"], ["padloc"], ["defense"])
    # GajA by its coordinates, GajB by its name.
    assert system.qualifiers["component"] == ["core|GajA|1264_-1_959|padloc|GajA_00001",
                                              "core|GajB|1391_-1_1266|padloc|GajB_00001"]
    assert (system.qualifiers["score"], system.qualifiers["evalue"], system.qualifiers["completeness"]) == (
        ["0.0"], ["2.5e-12"], ["."])
    # Importing again replaces them.
    with pytest.warns(UserWarning):
        record, = import_systems([record], source="padloc", calls=calls)
    assert len(systems(record)) == 1


def test_call_table_errors(pdonr, shared_datadir, tmp_path, capsys):
    from domainator.errors import InputError
    with pytest.raises(InputError, match="is not a DefenseFinder genes table"):
        list(import_systems(parse_seqfiles([pdonr]), source="defensefinder",
                            calls=[str(shared_datadir / "systems" / "padloc.csv")]))
    with pytest.raises(SystemExit):
        import_systems_module.main(["-i", pdonr, "--source", "padloc", "-o", str(tmp_path / "out.gb")])
    assert "needs --calls" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        import_systems_module.main(["-i", pdonr, "--calls", str(shared_datadir / "systems" / "padloc.csv"),
                                    "-o", str(tmp_path / "out.gb")])
    assert "--calls is for" in capsys.readouterr().err


def test_main_padloc(pdonr, shared_datadir, tmp_path):
    out = tmp_path / "out.gb"
    with pytest.warns(UserWarning):
        import_systems_module.main(["-i", pdonr, "--source", "padloc", "--calls", str(shared_datadir / "systems" / "padloc.csv"),
                                    "-o", str(out)])
    regions = list(select_by_system(parse_seqfiles([str(out)]), systems={"padloc:gabija"}, cds_range=(1, 1)))
    assert len(regions) == 1
