"""Convert MacSyFinder model packages, such as DefenseFinder's defense-finder-models, to native
system models.

Written from MacSyFinder's documented package format (models format version 2.0), and checked
against defense-finder-models 3.1.0. The mapping is in docs/design/systems/definition_formats.md
("MacSyFinder / DefenseFinder packages"). A package is a directory with:

    metadata.yml                       vers, license, ...
    definitions/[<family>/...]<model>.xml
    profiles/<gene>.hmm[.gz]           one HMM per gene, named after the gene
    model_conf.xml                     optional: the package's filtering options

A model's genes become components, named after the gene without its ``<family>__`` prefix;
exchangeable genes become extra references of their gene's component. Each profile is referred
to by its NAME line, which can differ from its file name. MacSyFinder's hit filters become
reference thresholds: i_evalue_sel, coverage_profile, and with cut_ga (the default) the
profile's domain gathering threshold as min_score.
"""
import gzip
import os
from pathlib import Path
import re
from typing import Dict, List, NamedTuple, Optional
import xml.etree.ElementTree as ElementTree

import yaml

from domainator import __version__
from domainator.errors import FormatVersionError, InputError
from domainator.systems.formats import Conversion, Note
from domainator.systems.model import Component, ModelCollection, ReferenceRule, Source, SystemModel, Window

DEFENSEFINDER_DATABASE = "defensefinder"
DEFENSEFINDER_REPOSITORY = "mdmparis/defense-finder-models"
MODELS_VERSION = "2.0"
DEFAULT_I_EVALUE_SEL = 0.001
DEFAULT_COVERAGE_PROFILE = 0.5
DEFENSEFINDER_COVERAGE_PROFILE = 0.4
"""DefenseFinder runs MacSyFinder with --coverage-profile 0.4."""
ROLES = {"mandatory": "core", "accessory": "accessory", "neutral": "neutral", "forbidden": "forbidden"}
_PREFIX_SEPARATOR = "__"


class Profile(NamedTuple):
    path: Path
    name: str
    """The NAME line, which hits are labeled with."""
    gathering: Optional[float]
    """The domain gathering threshold (the second GA value), if the profile has one."""


class Filtering(NamedTuple):
    i_evalue_sel: float
    coverage_profile: float
    cut_ga: bool


def is_package(path: os.PathLike) -> bool:
    path = Path(path)
    return (path / "metadata.yml").is_file() and (path / "definitions").is_dir()


def read_profile(path: Path) -> Profile:
    """The NAME and domain GA of an HMMER3 profile file (only its header is read)."""
    name = None
    gathering = None
    with (gzip.open(path, "rt") if path.name.endswith(".gz") else open(path)) as handle:
        for line in handle:
            if line.startswith("NAME "):
                name = line.split()[1]
            elif line.startswith("GA "):
                values = line[3:].replace(";", " ").split()
                gathering = float(values[1] if len(values) > 1 else values[0])
            elif line.startswith("HMM "):
                break
    if name is None:
        raise InputError(f"{path} is not an HMMER3 profile: it has no NAME line.")
    return Profile(path, name, gathering)


def _profile_paths(path: Path) -> Dict[str, Path]:
    """{gene name: profile file} for the package's profiles/ directory."""
    out = dict()
    profiles = path / "profiles"
    if profiles.is_dir():
        for file in sorted(profiles.iterdir()):
            for suffix in (".hmm.gz", ".hmm"):
                if file.name.endswith(suffix):
                    out.setdefault(file.name[:-len(suffix)], file)
    return out


def _filtering(path: Path, coverage_profile: Optional[float], default_coverage: float) -> Filtering:
    """The hit filters: model_conf.xml's, else MacSyFinder's defaults; coverage_profile, if
    given, overrides both."""
    i_evalue_sel, coverage, cut_ga = DEFAULT_I_EVALUE_SEL, default_coverage, True
    conf = path / "model_conf.xml"
    if conf.exists():
        filtering = _parse_xml(conf).find("filtering")
        if filtering is not None:
            value = filtering.findtext("i_evalue_sel")
            i_evalue_sel = float(value) if value else i_evalue_sel
            value = filtering.findtext("coverage_profile")
            coverage = float(value) if value else coverage
            value = filtering.findtext("cut_ga")
            cut_ga = cut_ga if not value else value.strip().lower() in ("1", "true", "yes")
    if coverage_profile is not None:
        coverage = coverage_profile
    return Filtering(i_evalue_sel, coverage, cut_ga)


def _parse_xml(path: Path) -> ElementTree.Element:
    try:
        return ElementTree.parse(path).getroot()
    except ElementTree.ParseError as e:
        raise InputError(f"{path} is not valid XML: {e}") from e


def _licenses(path: Path) -> Optional[str]:
    """Every license: value in metadata.yml. defense-finder-models lists two, and a YAML reader
    keeps only the last."""
    with open(path / "metadata.yml", encoding="utf-8") as handle:
        values = [line.split(":", 1)[1].strip() for line in handle if line.startswith("license:")]
    return "; ".join(values) or None


def _flag(element: ElementTree.Element, key: str) -> bool:
    return element.get(key, "0").strip().lower() in ("1", "true", "yes")


def _component_names(genes: List[str]) -> Dict[str, str]:
    """Component names: the gene names without their ``<family>__`` prefix, unless that would
    make two the same."""
    short = {gene: gene.split(_PREFIX_SEPARATOR, 1)[-1] for gene in genes}
    if len(set(short.values())) < len(short):
        return {gene: gene for gene in genes}
    return short


def convert(path: os.PathLike, *, database: str, release: Optional[str] = None, url: Optional[str] = None,
            coverage_profile: Optional[float] = None, default_coverage: float = DEFAULT_COVERAGE_PROFILE,
            source_name: Optional[str] = None) -> Conversion:
    """Convert the MacSyFinder package at path into a collection named database.

    coverage_profile overrides the package's (model_conf.xml) and MacSyFinder's default
    (default_coverage) minimum profile coverage. release defaults to metadata.yml's vers.
    """
    path = Path(path)
    if not is_package(path):
        raise InputError(f"{path} is not a MacSyFinder model package: it needs metadata.yml and definitions/.")
    with open(path / "metadata.yml", encoding="utf-8") as handle:
        try:
            metadata = yaml.safe_load(handle) or {}
        except yaml.YAMLError as e:
            raise InputError(f"{path / 'metadata.yml'} is not valid YAML: {e}") from e
    filtering = _filtering(path, coverage_profile, default_coverage)
    profile_paths = _profile_paths(path)
    profiles: Dict[str, Profile] = dict()

    def profile(gene: str) -> Optional[Profile]:
        if gene not in profiles and gene in profile_paths:
            profiles[gene] = read_profile(profile_paths[gene])
        return profiles.get(gene)

    definitions = path / "definitions"
    files = sorted(file for file in definitions.rglob("*.xml"))
    stems: Dict[str, int] = dict()
    for file in files:
        stems[file.stem] = stems.get(file.stem, 0) + 1

    notes: List[Note] = list()
    models: List[SystemModel] = list()
    for file in files:
        relative = file.relative_to(definitions).with_suffix("")
        accession = relative.as_posix()
        name = file.stem
        if stems[name] > 1:
            name = accession
            notes.append(Note(name, "name", f"qualified: {stems[file.stem]} models are named {file.stem}"))
        root = _parse_xml(file)
        if root.tag != "model" or root.get("vers") != MODELS_VERSION:
            version = root.get("vers", "1.0" if root.tag == "system" else "unknown")
            raise FormatVersionError(f"{file} is a MacSyFinder model of format version {version}; "
                                     f"Domainator converts version {MODELS_VERSION}.")
        if root.get("inter_gene_max_space") is None:
            raise InputError(f"{file}: the model has no inter_gene_max_space.")

        gene_elements = root.findall("gene")
        component_names = _component_names([gene.get("name", "") for gene in gene_elements])
        components = list()
        mandatory = 0
        seen_genes = set()  # a gene listed twice is one gene, as in MacSyFinder
        for gene in gene_elements:
            gene_name = gene.get("name")
            presence = gene.get("presence")
            if gene_name is None or presence not in ROLES:
                raise InputError(f"{file}: every gene needs a name and a presence (one of {', '.join(ROLES)}).")
            if gene_name in seen_genes:
                notes.append(Note(name, gene_name, "merged: the gene is listed more than once"))
                continue
            mandatory += presence == "mandatory"
            references = list()
            for reference_gene, exchangeable in [(gene_name, False)] + [
                    (g.get("name"), True) for exchangeables in gene.findall("exchangeables")
                    for g in exchangeables.findall("gene")]:
                if reference_gene in seen_genes:
                    notes.append(Note(name, reference_gene, "merged: the gene is listed more than once"))
                    continue
                seen_genes.add(reference_gene)
                found = profile(reference_gene)
                if found is None:
                    notes.append(Note(name, reference_gene, "dropped: no profile in profiles/ has this gene's name"))
                    continue
                references.append(ReferenceRule(
                    found.name, database=database, evalue=filtering.i_evalue_sel,
                    min_score=found.gathering if filtering.cut_ga else None,
                    min_reference_coverage=filtering.coverage_profile, exchangeable=exchangeable))
            if not references:
                continue
            component_name = component_names[gene_name]
            if "|" in component_name:
                raise InputError(f"{file}: gene name {gene_name!r} contains '|'.")
            gap = gene.get("inter_gene_max_space")
            loner = _flag(gene, "loner")
            if loner:
                notes.append(Note(name, component_name, "approximated: loner is evaluated as an ordinary gene "
                                                        "(it must lie within the cluster)"))
            components.append(Component(component_name, ROLES[presence], tuple(references), loner=loner,
                                        multi_system=_flag(gene, "multi_system"), multi_model=_flag(gene, "multi_model"),
                                        max_gap_cds=int(gap) if gap is not None else None))

        counted = sum(1 for component in components if component.role in ("core", "accessory"))
        if counted == 0:
            notes.append(Note(name, "model", "dropped: none of its mandatory or accessory genes has a profile"))
            continue
        min_core = int(root.get("min_mandatory_genes_required", mandatory))
        min_total = int(root.get("min_genes_required", mandatory))
        core = sum(1 for component in components if component.role == "core")
        if min_core > core:
            notes.append(Note(name, "min_mandatory_genes_required", f"lowered: from {min_core} to {core}, the mandatory genes left"))
            min_core = core
        if min_total > counted:
            notes.append(Note(name, "min_genes_required", f"lowered: from {min_total} to {counted}, the mandatory and accessory genes left"))
            min_total = counted
        multi_loci = _flag(root, "multi_loci")
        if multi_loci:
            notes.append(Note(name, "multi_loci", "approximated: systems are found as single loci"))
        max_components = root.get("max_nb_genes")
        models.append(SystemModel(
            name=name, database=database, program="macsyfinder", components=tuple(components),
            window=Window("chain", max_cds=int(root.get("inter_gene_max_space"))),
            min_core=min_core, min_total=min_total,
            category=relative.parts[0] if len(relative.parts) > 1 else ".",
            accession=accession, max_components=int(max_components) if max_components is not None else None,
            multi_loci=multi_loci))

    if not models:
        raise InputError(f"{path} has no models that can be converted.")
    used = sorted({rule.name for model in models for component in model.components for rule in component.references})
    by_name = {found.name: found for found in profiles.values()}
    reference_file = f"{database}.hmm"
    collection = ModelCollection(
        database=database, models=tuple(models), semantics="macsyfinder",
        source=Source(name=source_name or path.name, url=url,
                      release=release or (str(metadata["vers"]) if metadata.get("vers") is not None else None),
                      converter=f"domainator {__version__}", license=_licenses(path)),
        reference_files=(reference_file,))
    license_file = path / "LICENSE" if (path / "LICENSE").exists() else None
    return Conversion(collection, {reference_file: [by_name[name].path for name in used]}, license_file, notes)


def database_name(path: os.PathLike) -> str:
    """A collection name for a package: its directory's name, without a version suffix."""
    name = Path(path).resolve().name
    return re.sub(r"[-_.]?v?\d+(\.\d+)*$", "", name) or name


def defensefinder_asset(release: str) -> str:
    """The name of a defense-finder-models release's tarball asset."""
    return f"defense-finder-models-{release}.tar.gz"
