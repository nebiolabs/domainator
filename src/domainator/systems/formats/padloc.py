"""Convert padloc-db (https://github.com/padlocbio/padloc-db) to native system models.

The mapping is in docs/design/systems/definition_formats.md ("PADLOC"). It follows padloc's own
reading of the database (bin/padloc.R): hmm_meta.txt is read by column position, protein names
and group ("secondary") names are split on '|', group names in a model's gene lists are
expanded to the proteins in the group, and overlaps between the gene lists are removed in
padloc's order. Models are named by their yaml.name, the name padloc reports.
"""
import csv
from dataclasses import replace
import os
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Tuple

import yaml

from domainator import __version__
from domainator.errors import InputError
from domainator.systems.formats import Conversion, Note
from domainator.systems.model import Component, ModelCollection, ReferenceRule, Source, SystemModel, Window

DATABASE = "padloc"
REPOSITORY = "padlocbio/padloc-db"
DEFAULT_EVALUE = 1e-5
DEFAULT_COVERAGE = 0.3
CRISPR_ARRAY = "CRISPR_array"
NCRNA = "ncRNA"
_ROLES = (("core_genes", "core"), ("secondary_genes", "accessory"), ("neutral_genes", "neutral"),
          ("prohibited_genes", "forbidden"))


class HmmRow(NamedTuple):
    accession: str
    name: str
    proteins: Tuple[str, ...]
    groups: Tuple[str, ...]
    evalue: float
    reference_coverage: float
    target_coverage: float


def _na(value: Optional[str]) -> Optional[str]:
    value = None if value is None else value.strip()
    return None if value in (None, "", "NA") else value


def _names(value: Optional[str]) -> Tuple[str, ...]:
    """A '|'-separated list of names, without spaces, as padloc splits them."""
    value = _na(value)
    if value is None:
        return ()
    return tuple(name for name in value.replace(" ", "").split("|") if name)


def _number(value: Optional[str], default: float) -> float:
    value = _na(value)
    return default if value is None else float(value)


def read_hmm_meta(path: os.PathLike) -> List[HmmRow]:
    """The rows of hmm_meta.txt, read by column position (the header's spelling of the columns
    has varied), with padloc's defaults for missing thresholds."""
    rows = list()
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t", quotechar='"')
        next(reader, None)
        for line_number, fields in enumerate(reader, start=2):
            if not fields or not any(field.strip() for field in fields):
                continue
            if len(fields) < 9:
                raise InputError(f"{path}, line {line_number}: expected at least 9 tab-separated columns, found {len(fields)}.")
            accession, name = _na(fields[0]), _na(fields[1])
            proteins = _names(fields[3])
            if accession is None or name is None or not proteins:
                raise InputError(f"{path}, line {line_number}: hmm.accession, hmm.name and protein.name are required.")
            rows.append(HmmRow(accession, name, proteins, _names(fields[4]), _number(fields[6], DEFAULT_EVALUE),
                               _number(fields[7], DEFAULT_COVERAGE), _number(fields[8], DEFAULT_COVERAGE)))
    return rows


def _read_table(path: os.PathLike) -> List[Dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t", quotechar='"'))


def _gene_list(value) -> List[str]:
    """A gene list from a model file: NA, a missing list and NA entries are empty."""
    if value is None:
        return []
    if not isinstance(value, list):
        value = [value]
    out = list()
    for gene in value:
        gene = _na(None if gene is None else str(gene))
        if gene is not None and gene not in out:
            out.append(gene)
    return out


def _without(genes: List[str], *others: List[str]) -> List[str]:
    excluded = {gene for other in others for gene in other}
    return [gene for gene in genes if gene not in excluded]


def _boolean(value) -> bool:
    return value is True or str(value).upper() == "TRUE"


def convert(path: os.PathLike, *, release: Optional[str] = None, url: Optional[str] = None) -> Conversion:
    """Convert a padloc-db directory (with sys/, hmm/, hmm_meta.txt and sys_meta.txt)."""
    path = Path(path)
    for required in ("sys", "hmm", "hmm_meta.txt", "sys_meta.txt"):
        if not (path / required).exists():
            raise InputError(f"{path} is not a padloc-db directory: it has no {required}.")
    hmm_rows = read_hmm_meta(path / "hmm_meta.txt")
    by_protein: Dict[str, List[HmmRow]] = dict()
    by_group: Dict[str, List[str]] = dict()
    for row in hmm_rows:
        for protein in row.proteins:
            by_protein.setdefault(protein, []).append(row)
        for group in row.groups:
            proteins = by_group.setdefault(group, [])
            proteins.extend(protein for protein in row.proteins if protein not in proteins)
    cm_names = [_na(row.get("cm.name")) for row in _read_table(path / "cm_meta.txt")] if (path / "cm_meta.txt").exists() else []
    cm_names = [name for name in cm_names if name is not None]

    notes: List[Note] = list()
    models: List[SystemModel] = list()
    used_rows: Dict[str, HmmRow] = dict()
    uses_ncrna = False
    for meta in _read_table(path / "sys_meta.txt"):
        name = _na(meta.get("yaml.name"))
        if name is None or not _boolean(meta.get("search")):
            continue
        model_path = path / "sys" / f"{name}.yaml"
        if not model_path.exists():
            notes.append(Note(name, "model", "dropped: sys_meta.txt lists it, but sys/ has no model file"))
            continue
        with open(model_path, encoding="utf-8") as handle:
            params = yaml.safe_load(handle) or {}

        def expand(genes: List[str]) -> List[str]:
            out = list(genes)
            for gene in genes:
                out.extend(protein for protein in by_group.get(gene, []) if protein not in out)
            return out

        genes = {key: expand(_gene_list(params.get(key))) for key, _ in _ROLES}
        genes["secondary_genes"] = _without(genes["secondary_genes"], genes["core_genes"], genes["prohibited_genes"])
        genes["prohibited_genes"] = _without(genes["prohibited_genes"], genes["core_genes"])
        genes["neutral_genes"] = _without(genes["neutral_genes"], genes["core_genes"], genes["secondary_genes"],
                                          genes["prohibited_genes"])

        components = list()
        for key, role in _ROLES:
            for gene in genes[key]:
                if gene == CRISPR_ARRAY:
                    notes.append(Note(name, gene, f"dropped: Domainator has no CRISPR array detector ({role} component)"))
                    continue
                if gene == NCRNA:
                    if not cm_names:
                        notes.append(Note(name, gene, "dropped: the database has no covariance models (cm_meta.txt)"))
                        continue
                    uses_ncrna = True
                    references = tuple(ReferenceRule(cm, database=DATABASE) for cm in cm_names)
                elif gene in by_protein:
                    references = tuple(ReferenceRule(row.name, database=DATABASE, evalue=row.evalue,
                                                     min_reference_coverage=row.reference_coverage,
                                                     min_target_coverage=row.target_coverage)
                                       for row in by_protein[gene])
                    used_rows.update((row.accession, row) for row in by_protein[gene])
                elif gene in by_group:
                    continue  # a group name, expanded to its proteins above
                else:
                    notes.append(Note(name, gene, f"dropped: no profile in hmm_meta.txt has this protein name ({role} component)"))
                    continue
                components.append(Component(gene, role, references))

        core = sum(1 for component in components if component.role == "core")
        total = core + sum(1 for component in components if component.role == "accessory")
        if total == 0:
            notes.append(Note(name, "model", "dropped: none of its core or secondary genes can be detected"))
            continue
        min_core, min_total = int(params.get("minimum_core") or 0), int(params.get("minimum_total") or 0)
        if min_core > core:
            notes.append(Note(name, "minimum_core", f"lowered: from {min_core} to {core}, the core components left"))
            min_core = core
        if min_total > total:
            notes.append(Note(name, "minimum_total", f"lowered: from {min_total} to {total}, the core and secondary components left"))
            min_total = total
        literature = tuple(reference.strip() for reference in (_na(meta.get("references")) or "").split(";") if reference.strip())
        models.append(SystemModel(
            name=name, database=DATABASE, program="padloc", components=tuple(components),
            window=Window("chain", max_cds=int(params.get("maximum_separation") or 0)),
            min_core=min_core, min_total=min_total, category="defense",
            description=_na(meta.get("system")) or ".", accession=name,
            same_strand=_boolean(params.get("force_strand", False)), literature=literature))

    models = _with_superiors(models, notes)
    reference_files = [f"{DATABASE}.hmm"] + ([f"{DATABASE}.cm"] if uses_ncrna else [])
    collection = ModelCollection(
        database=DATABASE, models=tuple(models), semantics="padloc",
        source=Source(name="padloc-db", url=url, release=release, converter=f"domainator {__version__}", license="MIT"),
        reference_files=tuple(reference_files))
    profile_files = {f"{DATABASE}.hmm": _hmm_files(path, used_rows)}
    if uses_ncrna:
        profile_files[f"{DATABASE}.cm"] = _cm_files(path, cm_names)
    license_file = path / "LICENSE" if (path / "LICENSE").exists() else None
    return Conversion(collection, profile_files, license_file, notes)


def _with_superiors(models: List[SystemModel], notes: List[Note]) -> List[SystemModel]:
    """padloc drops an ``*_other`` (catch-all) system where it overlaps a canonical one. Here,
    every canonical model that shares a component with an ``_other`` model is its superior."""
    canonical = [model for model in models if "_other" not in model.name]
    out = list()
    for model in models:
        if "_other" in model.name:
            names = {component.name for component in model.components if component.role != "forbidden"}
            superiors = tuple(other.key for other in canonical
                              if names & {c.name for c in other.components if c.role != "forbidden"})
            if superiors:
                model = replace(model, superiors=superiors)
                notes.append(Note(model.name, "superiors", f"added: the {len(superiors)} canonical models it shares a component with"))
        out.append(model)
    return out


def _hmm_files(path: Path, rows: Dict[str, HmmRow]) -> List[Path]:
    """The profile file of each used hmm_meta row (hmm/<accession>.hmm), or the combined
    hmm/padlocdb.hmm that padloc builds, when the single files are missing."""
    files = [path / "hmm" / f"{accession}.hmm" for accession in sorted(rows)]
    missing = [file for file in files if not file.exists()]
    if not missing:
        return files
    combined = path / "hmm" / "padlocdb.hmm"
    if combined.exists():
        return [combined]
    raise InputError(f"{path}/hmm has no profile file for {len(missing)} of the profiles the models use, "
                     f"for example {missing[0].name}, and no combined padlocdb.hmm.")


def _cm_files(path: Path, cm_names: List[str]) -> List[Path]:
    files = [path / "cm" / f"{name}.cm" for name in cm_names]
    missing = [file.name for file in files if not file.exists()]
    if missing:
        raise InputError(f"{path}/cm has no covariance model file for {', '.join(missing[:5])}.")
    return files
