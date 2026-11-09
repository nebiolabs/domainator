"""Reading and writing the native system model format (YAML with a ``domainator_systems:`` key).

See docs/design/systems/definition_formats.md and docs/file_formats.md. Commands read only
this format; foreign formats (PADLOC, MacSyFinder, antiSMASH) are refused with a message
saying they must be converted.
"""
from dataclasses import replace
import os
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import yaml

from domainator.errors import FormatVersionError, InputError
from domainator.seqio import db_name_from_path
from domainator.systems.conditions import Condition
from domainator.systems.model import (COMMAND_LINE_DATABASE, ROLES, THRESHOLD_KEYS, WINDOW_MODES, Component, Example,
                                      ModelCollection, ReferenceRule, Source, SystemModel, Window)

FORMAT_KEY = "domainator_systems"
FORMAT_VERSION = (1, 0)
"""The newest format version this module reads and the one it writes."""
DIRECTORY_FILE = "systems.yaml"
"""The file a model directory (written by system_db_download) holds."""
SEMANTICS = ("domainator", "padloc", "macsyfinder", "antismash")

_COLLECTION_KEYS = {FORMAT_KEY, "database", "semantics", "source", "reference_files", "defaults", "competing",
                    "competing_overlap_aa", "systems"}
_SOURCE_KEYS = {"name", "url", "release", "date", "converter", "license"}
_DEFAULTS_KEYS = set(THRESHOLD_KEYS) | {"window", "same_strand"}
_MODEL_KEYS = {"name", "category", "description", "accession", "window", "min_core", "min_total", "condition",
               "components", "superiors", "neighborhood_kb", "neighborhood_cds", "same_strand", "max_components",
               "multi_loci", "orientation_from", "examples", "literature", "reference_files"}
_COMPONENT_KEYS = {"name", "role", "domains", "loner", "multi_system", "multi_model", "max_gap_cds"}
_RULE_KEYS = {"name", "database"} | set(THRESHOLD_KEYS) | {"exchangeable", "attach_to_cds"}
_WINDOW_KEYS = {"mode", "max_cds", "max_kb", "max_span_cds", "max_span_kb"}
_EXAMPLE_KEYS = {"accession", "start", "end", "note"}
_FORMAT_1_1_KEYS = {"define": "define aliases", "extenders": "extenders"}
_INTEGER_THRESHOLDS = {"min_target_length", "max_target_length"}

_OPERATOR_CHARACTERS = set("&|~()<")


class _Mapping(dict):
    """A YAML mapping that remembers the line it started on."""
    line: Optional[int] = None


# libyaml's loader, when PyYAML was built with it: converted collections are megabytes of YAML.
class _Loader(getattr(yaml, "CSafeLoader", yaml.SafeLoader)):
    pass


def _construct_mapping(loader, node):
    mapping = _Mapping(loader.construct_mapping(node, deep=True))
    mapping.line = node.start_mark.line + 1
    return mapping


_Loader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping)


class _Reader:
    """Reads one collection, raising InputError with the file and line of the problem."""

    def __init__(self, path: Optional[str]):
        self.path = path

    def error(self, message: str, node: Any = None) -> InputError:
        where = self.path or "system models"
        line = getattr(node, "line", None)
        if line is not None:
            where = f"{where}, line {line}"
        return InputError(f"{where}: {message}")

    def check_keys(self, mapping, allowed, what):
        if not isinstance(mapping, dict):
            raise self.error(f"{what} must be a mapping, not {mapping!r}.")
        for key in mapping:
            if key in _FORMAT_1_1_KEYS and key not in allowed:
                raise self.error(f"{_FORMAT_1_1_KEYS[key]} ({key!r}) need system model format 1.1; "
                                 f"this version of Domainator reads format 1.", mapping)
            if key not in allowed:
                raise self.error(f"unknown key {key!r} in {what}. Known keys: {', '.join(sorted(allowed))}.", mapping)

    def text(self, value, what, node=None) -> str:
        if not isinstance(value, str) or not value.strip():
            raise self.error(f"{what} must be text, not {value!r} (quote it if YAML reads it as something else).", node)
        return value

    def name(self, value, what, node=None, *, forbidden=_OPERATOR_CHARACTERS | {":"}) -> str:
        value = self.text(value, what, node)
        bad = sorted({c for c in value if c in forbidden or c.isspace()})
        if bad:
            raise self.error(f"{what} {value!r} can't contain {' '.join(repr(c) for c in bad)}.", node)
        return value

    def number(self, value, what, node=None, *, integer=False) -> Union[int, float]:
        # YAML 1.1 reads 1e-5 (no dot) as text, so numbers may arrive as strings.
        try:
            if isinstance(value, bool):
                raise ValueError
            number = float(value)
        except (TypeError, ValueError):
            raise self.error(f"{what} must be a number, not {value!r}.", node) from None
        if number < 0:
            raise self.error(f"{what} must not be negative, not {value!r}.", node)
        if integer:
            if number != int(number):
                raise self.error(f"{what} must be a whole number, not {value!r}.", node)
            return int(number)
        return number

    def boolean(self, value, what, node=None) -> bool:
        if not isinstance(value, bool):
            raise self.error(f"{what} must be true or false, not {value!r}.", node)
        return value

    def text_list(self, value, what, node=None) -> Tuple[str, ...]:
        if not isinstance(value, list):
            raise self.error(f"{what} must be a list.", node)
        return tuple(self.text(item, f"each of {what}", node) for item in value)


def parse_version(value) -> Tuple[int, int]:
    """The format version as (major, minor), from ``1``, ``1.0`` or ``"1.1"``. Raises ValueError."""
    if isinstance(value, bool):
        raise ValueError(value)
    if isinstance(value, int):
        return (value, 0)
    text = str(value).strip()
    major, _, minor = text.partition(".")
    return (int(major), int(minor) if minor else 0)


def _version_text(version: Tuple[int, int]) -> str:
    return f"{version[0]}.{version[1]}"


def _foreign_format(path: Path) -> Optional[str]:
    """The name of the foreign model format ``path`` looks like, if any."""
    if path.is_dir():
        if (path / "sys").is_dir() and (path / "hmm_meta.txt").exists():
            return "padloc-db directory"
        if (path / "metadata.yml").exists() and (path / "definitions").is_dir():
            return "MacSyFinder model package (such as DefenseFinder's)"
        if (path / "cluster_rules").is_dir() or (path / "hmmdetails.txt").exists() or (path / "data" / "hmmdetails.txt").exists():
            return "antiSMASH hmm_detection directory"
        return None
    suffix = path.suffix.lower()
    if suffix == ".xml":
        return "MacSyFinder model definition"
    if suffix == ".txt":
        try:
            with open(path, encoding="utf-8", errors="replace") as handle:
                head = handle.read(65536)
        except OSError:
            return None
        if any(line.lstrip().startswith("RULE ") for line in head.splitlines()):
            return "antiSMASH cluster rule file"
    return None


def _foreign_error(path, label) -> InputError:
    message = (f"{path} looks like a {label}. Domainator reads only its own system model format "
               f"(a YAML file with a '{FORMAT_KEY}:' key, or a directory holding {DIRECTORY_FILE}), "
               f"so a {label} must be converted to it first")
    if label.lower().startswith("padloc"):
        directory = path if label.startswith("padloc-db") else "<its padloc-db directory>"
        return InputError(f"{message}: domainator system_db_download --source padloc --input {directory} --output_dir <directory>")
    return InputError(f"{message}. This version of Domainator can't convert it yet.")


def _load_yaml(path: Path):
    try:
        with open(path, encoding="utf-8") as handle:
            return yaml.load(handle, Loader=_Loader)  # nosec: _Loader is a SafeLoader
    except OSError as e:
        raise InputError(f"Can't read system models from {path}: {e.strerror}.") from e
    except yaml.YAMLError as e:
        raise InputError(f"{path} is not valid YAML: {e}") from e


def read_collection(path: Union[str, os.PathLike]) -> ModelCollection:
    """Read one native model file, or a directory holding systems.yaml."""
    path = Path(path)
    if not path.exists():
        raise InputError(f"System model file or directory not found: {path}")
    foreign = _foreign_format(path)
    if path.is_dir():
        if (path / DIRECTORY_FILE).exists():
            file = path / DIRECTORY_FILE
            default_database = path.resolve().name
        elif foreign is not None:
            raise _foreign_error(path, foreign)
        else:
            raise InputError(f"{path} is a directory without a {DIRECTORY_FILE}; give a system model file, or a directory holding {DIRECTORY_FILE}.")
    else:
        if foreign is not None:
            raise _foreign_error(path, foreign)
        file = path
        default_database = db_name_from_path(path)
    data = _load_yaml(file)
    if isinstance(data, dict) and FORMAT_KEY not in data:
        if {"core_genes", "maximum_separation", "minimum_core"} & set(data):
            raise _foreign_error(file, "PADLOC system model")
        if {"maintainer", "vers", "cite"} & set(data):
            raise _foreign_error(file, "MacSyFinder package metadata file")
    if not isinstance(data, dict) or FORMAT_KEY not in data:
        raise InputError(f"{file} is not a system model file: it has no top-level '{FORMAT_KEY}:' key giving the format version.")
    return collection_from_data(data, path=str(file), default_database=default_database)


def collection_from_data(data: Dict[str, Any], *, path: Optional[str] = None,
                         default_database: str = "systems") -> ModelCollection:
    """A ModelCollection from the parsed contents of a native model file."""
    reader = _Reader(path)
    reader.check_keys(data, _COLLECTION_KEYS, "the file")
    try:
        version = parse_version(data[FORMAT_KEY])
    except (TypeError, ValueError):
        raise reader.error(f"{FORMAT_KEY} must be a version number like 1 or \"1.1\", not {data[FORMAT_KEY]!r}.", data) from None
    if version[0] != FORMAT_VERSION[0] or version > FORMAT_VERSION:
        raise FormatVersionError(f"{path or 'system models'}: system model format {_version_text(version)} can't be read by this "
                                 f"version of Domainator, which reads format {_version_text(FORMAT_VERSION)}.")
    database = reader.name(data.get("database", default_database), "database", data)
    semantics = data.get("semantics", "domainator")
    if semantics not in SEMANTICS:
        raise reader.error(f"semantics must be one of {', '.join(SEMANTICS)}, not {semantics!r}.", data)
    source = None
    if data.get("source") is not None:
        reader.check_keys(data["source"], _SOURCE_KEYS, "source")
        source = Source(**{key: str(value) for key, value in data["source"].items() if value is not None})
    reference_files = reader.text_list(data.get("reference_files", []), "reference_files", data)

    defaults = data.get("defaults") or _Mapping()
    reader.check_keys(defaults, _DEFAULTS_KEYS, "defaults")
    default_thresholds = {key: _threshold(reader, key, defaults[key], defaults) for key in THRESHOLD_KEYS if defaults.get(key) is not None}
    default_window = _window(reader, defaults["window"]) if defaults.get("window") is not None else None
    default_same_strand = reader.boolean(defaults.get("same_strand", False), "same_strand", defaults)

    competing = list()
    for group in data.get("competing") or []:
        competing.append(frozenset(reader.text_list(group, "each competing group", data)))
    competing_overlap_aa = reader.number(data.get("competing_overlap_aa", 20), "competing_overlap_aa", data, integer=True)

    systems = data.get("systems")
    if isinstance(systems, dict):
        entries = list()
        for name, condition in systems.items():
            if not isinstance(condition, str):
                raise reader.error(f"in the short form ('systems:' as a mapping) each value is a condition; {name!r} has {condition!r}.", systems)
            entry = _Mapping(name=name, condition=condition)
            entry.line = getattr(systems, "line", None)
            entries.append(entry)
    elif isinstance(systems, list):
        entries = systems
    else:
        raise reader.error("'systems:' must be a list of models, or a mapping from model name to condition.", data)
    if not entries:
        raise reader.error("'systems:' has no models.", data)

    parsed = list()
    raw_superiors = list()
    for entry in entries:
        model, superiors = _model(reader, entry, database=database, program=semantics, default_thresholds=default_thresholds,
                                  default_window=default_window, default_same_strand=default_same_strand)
        parsed.append(model)
        raw_superiors.append((entry, superiors))
    names = [model.name for model in parsed]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise reader.error(f"duplicate model names: {', '.join(duplicates)}.")
    models = list()
    for model, (entry, superiors) in zip(parsed, raw_superiors):
        keys = list()
        for superior in superiors:
            superior_database, _, superior_name = superior.rpartition(":")
            if not superior_database:
                if superior_name not in names:
                    raise reader.error(f"model {model.name!r} lists superior {superior_name!r}, which is not a model in this file "
                                       "(write database:name for a model in another collection).", entry)
                superior_database = database
            keys.append((superior_database, superior_name))
        models.append(replace(model, superiors=tuple(keys)))
    return ModelCollection(database=database, models=tuple(models), semantics=semantics, source=source,
                           competing=tuple(competing), competing_overlap_aa=competing_overlap_aa,
                           reference_files=reference_files, path=path)


def _threshold(reader, key, value, node):
    return reader.number(value, key, node, integer=key in _INTEGER_THRESHOLDS)


def _window(reader, value) -> Window:
    reader.check_keys(value, _WINDOW_KEYS, "window")
    mode = value.get("mode", WINDOW_MODES[0])
    if mode not in WINDOW_MODES:
        raise reader.error(f"window mode must be one of {', '.join(WINDOW_MODES)}, not {mode!r}.", value)
    limits = {key: reader.number(value[key], key, value, integer=key.endswith("cds"))
              for key in ("max_cds", "max_kb", "max_span_cds", "max_span_kb") if value.get(key) is not None}
    if "max_cds" not in limits and "max_kb" not in limits:
        raise reader.error("a window needs max_cds, max_kb or both.", value)
    if mode == "span" and ("max_span_cds" in limits or "max_span_kb" in limits):
        raise reader.error("max_span_cds and max_span_kb are for chain and anchor windows; in span mode the span is the window.", value)
    return Window(mode=mode, **limits)


def implicit_components(condition: Condition, thresholds: Dict[str, Any]) -> Tuple[Component, ...]:
    """The components of a model with no components list: one per name in the condition, with that
    one reference (``database:name`` limits it to one database). Names used only under ``~`` are
    forbidden; the rest are core."""
    components = list()
    seen = set()
    for token in condition._tokens:
        if token not in condition.names or token in seen:
            continue
        seen.add(token)
        database, _, name = token.rpartition(":")
        rule = ReferenceRule(name=name, database=database or None, **thresholds)
        role = "core" if token in condition.positive_names else "forbidden"
        components.append(Component(name=token, role=role, references=(rule,)))
    return tuple(components)


def _default_min_core(condition: Optional[Condition], components: Sequence[Component]) -> int:
    if condition is not None:
        return 0
    return sum(1 for component in components if component.role == "core")


def _rule(reader, value, thresholds, node) -> ReferenceRule:
    if isinstance(value, str):
        return ReferenceRule(name=reader.text(value, "a domain name", node), **thresholds)
    reader.check_keys(value, _RULE_KEYS, "a domain")
    if "name" not in value:
        raise reader.error("each domain needs a name.", value)
    fields = dict(thresholds)
    for key in THRESHOLD_KEYS:
        if value.get(key) is not None:
            fields[key] = _threshold(reader, key, value[key], value)
    return ReferenceRule(name=reader.text(value["name"], "a domain name", value),
                         database=None if value.get("database") is None else reader.text(value["database"], "database", value),
                         exchangeable=reader.boolean(value.get("exchangeable", False), "exchangeable", value),
                         attach_to_cds=reader.boolean(value.get("attach_to_cds", False), "attach_to_cds", value),
                         **fields)


def _component(reader, value, thresholds) -> Component:
    reader.check_keys(value, _COMPONENT_KEYS, "a component")
    if "name" not in value:
        raise reader.error("each component needs a name.", value)
    name = reader.name(value["name"], "component name", value, forbidden=_OPERATOR_CHARACTERS)
    role = value.get("role", "core")
    if role not in ROLES:
        raise reader.error(f"component {name!r}: role must be one of {', '.join(ROLES)}, not {role!r}.", value)
    domains = value.get("domains")
    if not isinstance(domains, list) or not domains:
        raise reader.error(f"component {name!r} needs a list of domains.", value)
    return Component(name=name, role=role, references=tuple(_rule(reader, d, thresholds, value) for d in domains),
                     loner=reader.boolean(value.get("loner", False), "loner", value),
                     multi_system=reader.boolean(value.get("multi_system", False), "multi_system", value),
                     multi_model=reader.boolean(value.get("multi_model", False), "multi_model", value),
                     max_gap_cds=None if value.get("max_gap_cds") is None else reader.number(value["max_gap_cds"], "max_gap_cds", value, integer=True))


def _model(reader, entry, *, database, program, default_thresholds, default_window, default_same_strand) -> Tuple[SystemModel, List[str]]:
    reader.check_keys(entry, _MODEL_KEYS, "a model")
    if "name" not in entry:
        raise reader.error("each model needs a name.", entry)
    name = reader.name(entry["name"], "model name", entry)

    def error(message):
        return reader.error(f"model {name!r}: {message}", entry)

    window = _window(reader, entry["window"]) if entry.get("window") is not None else default_window
    if window is None:
        raise error("no window (give one, or a window in defaults).")
    condition = None
    if entry.get("condition") is not None:
        try:
            condition = Condition(reader.text(entry["condition"], "condition", entry))
        except InputError as e:
            raise error(str(e)) from None
    if entry.get("components") is not None:
        if not isinstance(entry["components"], list) or not entry["components"]:
            raise error("components must be a non-empty list.")
        components = tuple(_component(reader, c, default_thresholds) for c in entry["components"])
        component_names = [c.name for c in components]
        duplicates = sorted({c for c in component_names if component_names.count(c) > 1})
        if duplicates:
            raise error(f"duplicate components: {', '.join(duplicates)}.")
        if condition is not None:
            unknown = sorted(condition.names - set(component_names))
            if unknown:
                raise error(f"the condition names {', '.join(unknown)}, which are not components.")
    else:
        if condition is None:
            raise error("a model needs components, a condition, or both.")
        bad = sorted(n for n in condition.names if "|" in n)
        if bad:
            raise error(f"names can't contain '|': {', '.join(bad)}.")
        components = implicit_components(condition, default_thresholds)
    min_core = entry.get("min_core")
    min_core = _default_min_core(condition, components) if min_core is None else reader.number(min_core, "min_core", entry, integer=True)
    min_total = entry.get("min_total")
    min_total = min_core if min_total is None else reader.number(min_total, "min_total", entry, integer=True)
    if min_core == 0 and min_total == 0 and condition is None:
        raise error("no requirement: give min_core, min_total or a condition.")
    orientation_from = entry.get("orientation_from")
    if orientation_from is not None and orientation_from not in {c.name for c in components}:
        raise error(f"orientation_from names {orientation_from!r}, which is not a component.")
    superiors = list(reader.text_list(entry.get("superiors", []), "superiors", entry))
    examples = list()
    for example in entry.get("examples") or []:
        reader.check_keys(example, _EXAMPLE_KEYS, "an example")
        if "accession" not in example:
            raise reader.error("each example needs an accession.", example)
        examples.append(Example(accession=str(example["accession"]),
                                start=None if example.get("start") is None else reader.number(example["start"], "start", example, integer=True),
                                end=None if example.get("end") is None else reader.number(example["end"], "end", example, integer=True),
                                note=None if example.get("note") is None else str(example["note"])))

    def optional_number(key, integer=False):
        return None if entry.get(key) is None else reader.number(entry[key], key, entry, integer=integer)

    model = SystemModel(
        name=name, database=database, program=program, components=components, window=window,
        min_core=min_core, min_total=min_total, condition=None if condition is None else condition.text,
        category=reader.text(entry.get("category", "."), "category", entry),
        neighborhood_kb=optional_number("neighborhood_kb"),
        neighborhood_cds=optional_number("neighborhood_cds", integer=True),
        description=reader.text(entry.get("description", "."), "description", entry),
        accession=reader.text(str(entry.get("accession", ".")), "accession", entry),
        same_strand=reader.boolean(entry.get("same_strand", default_same_strand), "same_strand", entry),
        max_components=optional_number("max_components", integer=True),
        multi_loci=reader.boolean(entry.get("multi_loci", False), "multi_loci", entry),
        orientation_from=orientation_from,
        examples=tuple(examples),
        literature=reader.text_list(entry.get("literature", []), "literature", entry),
        reference_files=reader.text_list(entry.get("reference_files", []), "reference_files", entry),
    )
    return model, superiors


def read_system_models(sources: Iterable[Union[str, os.PathLike, ModelCollection]]) -> List[ModelCollection]:
    """Read native model files and directories (and pass ModelCollections through).

    Model (database, name) pairs must be unique across all of them, and every superior must exist.
    """
    collections = [source if isinstance(source, ModelCollection) else read_collection(source) for source in sources]
    keys = dict()
    for collection in collections:
        for model in collection.models:
            if model.key in keys:
                raise InputError(f"Model {model.database}:{model.name} is defined twice ({keys[model.key] or 'system models'} "
                                 f"and {collection.path or 'system models'}).")
            keys[model.key] = collection.path
    for collection in collections:
        for model in collection.models:
            for superior in model.superiors:
                if superior not in keys:
                    raise InputError(f"Model {model.database}:{model.name} lists superior {superior[0]}:{superior[1]}, "
                                     "which is not in any of the system models given.")
    return collections


def select_models(collections: Sequence[ModelCollection], systems: Optional[Iterable[str]] = None,
                  system_categories: Optional[Iterable[str]] = None) -> List[ModelCollection]:
    """The collections, keeping only the models named in ``systems`` (``name`` or
    ``database:name``) and in ``system_categories``; None selects every model. Collections left
    with no models are dropped. A name in systems that names no model is an InputError."""
    systems = set(systems) if systems is not None else None
    system_categories = set(system_categories) if system_categories is not None else None
    unmatched = set(systems) if systems is not None else set()
    out = list()
    for collection in collections:
        models = list()
        for model in collection.models:
            names = {model.name, f"{model.database}:{model.name}"}
            if systems is not None:
                if not names & systems:
                    continue
                unmatched -= names
            if system_categories is not None and model.category not in system_categories:
                continue
            models.append(model)
        if models:
            out.append(replace(collection, models=tuple(models)))
    if unmatched:
        raise InputError(f"--systems names models that aren't in the system models given: {', '.join(sorted(unmatched))}.")
    return out


def reference_paths(collection: ModelCollection) -> List[str]:
    """The reference files of a collection and its models, without repeats. Relative paths are
    relative to the directory of the collection's file."""
    base = Path(collection.path).parent if collection.path else Path(".")
    files = list(collection.reference_files)
    for model in collection.models:
        files.extend(file for file in model.reference_files if file not in files)
    return [str(base / file) for file in files]


def window_model(domain_expr: str, *, window_cds: Optional[int] = None, window_kb: Optional[float] = None,
                 window_mode: str = "span", system_name: str = "domain_expr", same_strand: bool = False) -> ModelCollection:
    """The model --domain_expr and the --window options define: a short-form model in the
    ``command_line`` database, with the expression as its description."""
    if window_cds is None and window_kb is None:
        raise InputError("--domain_expr needs a window: give --window_cds, --window_kb or both.")
    window = _Mapping(mode=window_mode)
    if window_cds is not None:
        window["max_cds"] = window_cds
    if window_kb is not None:
        window["max_kb"] = window_kb
    entry = _Mapping(name=system_name, condition=domain_expr, window=window, description=domain_expr, same_strand=same_strand)
    reader = _Reader("the command line")
    model, _ = _model(reader, entry, database=COMMAND_LINE_DATABASE, program="domainator", default_thresholds={},
                      default_window=None, default_same_strand=False)
    return ModelCollection(database=COMMAND_LINE_DATABASE, models=(model,))


# ---- writing


def _window_data(window: Window) -> Dict[str, Any]:
    data = {"mode": window.mode}
    for key in ("max_cds", "max_kb", "max_span_cds", "max_span_kb"):
        if getattr(window, key) is not None:
            data[key] = getattr(window, key)
    return data


def _common_defaults(models: Sequence[SystemModel]) -> Dict[str, Any]:
    """The thresholds, window and same_strand shared by every rule and model, to write once in ``defaults``."""
    defaults = dict()
    rules = [rule for model in models for component in model.components for rule in component.references]
    for key in THRESHOLD_KEYS:
        values = {getattr(rule, key) for rule in rules}
        if len(values) == 1 and None not in values:
            defaults[key] = values.pop()
    windows = {model.window for model in models}
    if len(windows) == 1:
        defaults["window"] = windows.pop()
    if models and all(model.same_strand for model in models):
        defaults["same_strand"] = True
    return defaults


def _rule_data(rule: ReferenceRule, defaults: Dict[str, Any]):
    data = {"name": rule.name}
    if rule.database is not None:
        data["database"] = rule.database
    for key in THRESHOLD_KEYS:
        value = getattr(rule, key)
        if value is not None and value != defaults.get(key):
            data[key] = value
    for key in ("exchangeable", "attach_to_cds"):
        if getattr(rule, key):
            data[key] = True
    return rule.name if len(data) == 1 else data


def _component_data(component: Component, defaults: Dict[str, Any]) -> Dict[str, Any]:
    data = {"name": component.name, "role": component.role,
            "domains": [_rule_data(rule, defaults) for rule in component.references]}
    for key in ("loner", "multi_system", "multi_model"):
        if getattr(component, key):
            data[key] = True
    if component.max_gap_cds is not None:
        data["max_gap_cds"] = component.max_gap_cds
    return data


def _model_data(model: SystemModel, defaults: Dict[str, Any]) -> Dict[str, Any]:
    data = {"name": model.name}
    for key in ("category", "description", "accession"):
        if getattr(model, key) != ".":
            data[key] = getattr(model, key)
    if model.window != defaults.get("window"):
        data["window"] = _window_data(model.window)
    if model.same_strand != defaults.get("same_strand", False):
        data["same_strand"] = model.same_strand
    condition = None
    if model.condition is not None:
        data["condition"] = model.condition
        condition = Condition(model.condition)
    thresholds = {key: defaults[key] for key in THRESHOLD_KEYS if key in defaults}
    if condition is None or model.components != implicit_components(condition, thresholds):
        data["components"] = [_component_data(component, defaults) for component in model.components]
    if model.min_core != _default_min_core(condition, model.components):
        data["min_core"] = model.min_core
    if model.min_total != model.min_core:
        data["min_total"] = model.min_total
    if model.superiors:
        data["superiors"] = [name if database == model.database else f"{database}:{name}" for database, name in model.superiors]
    for key in ("neighborhood_kb", "neighborhood_cds", "max_components", "orientation_from"):
        if getattr(model, key) is not None:
            data[key] = getattr(model, key)
    if model.multi_loci:
        data["multi_loci"] = True
    if model.examples:
        data["examples"] = [{key: value for key, value in vars(example).items() if value is not None} for example in model.examples]
    for key in ("literature", "reference_files"):
        if getattr(model, key):
            data[key] = list(getattr(model, key))
    return data


def collection_to_data(collection: ModelCollection) -> Dict[str, Any]:
    """The native-format contents for ``collection``. Models that fit the short form are written in it."""
    data: Dict[str, Any] = {FORMAT_KEY: FORMAT_VERSION[0], "database": collection.database}
    if collection.semantics != "domainator":
        data["semantics"] = collection.semantics
    if collection.source is not None:
        data["source"] = {key: value for key, value in vars(collection.source).items() if value is not None}
    if collection.reference_files:
        data["reference_files"] = list(collection.reference_files)
    defaults = _common_defaults(collection.models)
    if defaults:
        data["defaults"] = {key: (_window_data(value) if key == "window" else value) for key, value in defaults.items()}
    if collection.competing:
        data["competing"] = [sorted(group) for group in collection.competing]
    if collection.competing_overlap_aa != 20:
        data["competing_overlap_aa"] = collection.competing_overlap_aa
    entries = [_model_data(model, defaults) for model in collection.models]
    if all(set(entry) == {"name", "condition"} for entry in entries):
        data["systems"] = {entry["name"]: entry["condition"] for entry in entries}
    else:
        data["systems"] = entries
    return data


def write_system_models(collection: ModelCollection, path: Union[str, os.PathLike]) -> None:
    """Write ``collection`` as a native model file."""
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(collection_to_data(collection), handle, sort_keys=False, default_flow_style=None, allow_unicode=True)
