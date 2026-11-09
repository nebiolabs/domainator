"""Searching for the profiles that system models name (domainate --system_models).

See "Thresholds are applied at detection time" in docs/design/systems/definition_formats.md.
"""
import os
from typing import Dict, List, Sequence, Tuple

from domainator.seqio import db_name_from_path
from domainator.systems.model import ModelCollection
from domainator.systems.native import reference_paths


def model_reference_files(collections: Sequence[ModelCollection]) -> List[str]:
    """The reference files of the collections and their models, without repeats."""
    out = list()
    seen = set()
    for collection in collections:
        for path in reference_paths(collection):
            real = os.path.realpath(path)
            if real not in seen:
                seen.add(real)
                out.append(path)
    return out


def model_profile_evalues(collections: Sequence[ModelCollection], reference_files: Sequence[str],
                          evalue: float) -> Dict[Tuple[str, str], float]:
    """The e-value to search each model profile with, by (database, name): the larger of evalue
    and the loosest e-value threshold of the model rules that name it.

    The database of a profile is the name domainate gives its file (seqio.db_name_from_path). A
    rule with no database names the profile in every reference file.
    """
    databases = [db_name_from_path(path) for path in reference_files]
    out: Dict[Tuple[str, str], float] = dict()
    for collection in collections:
        for model in collection.models:
            for component in model.components:
                for rule in component.references:
                    threshold = evalue if rule.evalue is None else max(evalue, rule.evalue)
                    for database in databases:
                        if rule.database is None or rule.database == database:
                            key = (database, rule.name)
                            out[key] = max(out.get(key, threshold), threshold)
    return out
