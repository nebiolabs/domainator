"""Convert antiSMASH's cluster rules (detection/hmm_detection) to native system models.

The mapping is in docs/design/systems/antismash.md ("Converting antiSMASH rules"). The rule
language is read as documented in antiSMASH's rule parser (common/hmm_rule_parser/rule_parser.py):

    DEFINE ANY_KS AS (PKS_KS or ene_KS)
    RULE transAT-PKS
        CATEGORY PKS
        DESCRIPTION Trans-AT polyketide
        EXAMPLE NCBI CP006259.1 382565-499230 kirromycin
        RELATED PP-binding
        SUPERIORS ...
        CUTOFF 20
        NEIGHBOURHOOD 20
        CONDITIONS cds(PKS_AT and not ANY_KS) and cds(ATd and ANY_KS)
        EXTENDERS cds(ANY_KS and not PKS_AT)

``#`` starts a comment, and a section runs until the next keyword, across lines. Conditions keep
their functions (``cds``, ``minimum``, ``repeatable``, ``minscore``) and are spelled with ``&``,
``|`` and ``~``; DEFINE aliases become the collection's ``define``. Each rule is an ``anchor``
model with antiSMASH's CUTOFF (whose ``<`` becomes ``<=`` one bp sooner) and its NEIGHBOURHOOD as
``neighborhood_kb``. Profiles are hmmdetails.txt's, with its cutoffs as ``min_score``, and
filterhmmdetails.txt's groups are the collection's competing groups. Dynamic profiles are code,
so their components have no references unless a substitutes file gives some.

antiSMASH's rules and profiles are AGPL-3.0 licensed: they are read from the user's own copy,
and the converted directory carries that license.
"""
import datetime
from pathlib import Path
import re
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

from domainator import __version__
from domainator.errors import InputError
from domainator.expr import tokenize
from domainator.systems.conditions import FUNCTIONS, OPERATORS, Condition
from domainator.systems.formats import Conversion, Note
from domainator.systems.model import Component, Example, ModelCollection, ReferenceRule, Source, SystemModel, Window
from domainator.systems.native import implicit_components, read_substitutes

DATABASE = "antismash"
STRICTNESS = ("strict", "relaxed", "loose")
"""antiSMASH's rule files, loaded cumulatively: relaxed adds relaxed.txt to strict.txt, and loose adds loose.txt."""
DEFAULT_STRICTNESS = "relaxed"
_MARKERS = {"RULE", "CATEGORY", "DESCRIPTION", "EXAMPLE", "RELATED", "SUPERIORS", "CUTOFF", "NEIGHBOURHOOD",
            "CONDITIONS", "EXTENDERS", "DEFINE"}
_REPEATABLE = {"EXAMPLE"}
_OPERATORS = {"and": "&", "or": "|", "not": "~"}
_LICENSE_NAMES = ("LICENSE.txt", "LICENSE", "COPYING")
LICENSE_NOTICE = (
    "This directory was converted by Domainator's system_db_download from antiSMASH's cluster rules and\n"
    "profiles, which are licensed under the GNU Affero General Public License v3.0 (AGPL-3.0).\n"
    "It is a derived work of those files: keep it for your own use, or share it under the AGPL.\n")


class Profile(NamedTuple):
    identifier: str
    name: str
    """The profile's NAME line: what Domainator hits are named by."""
    cutoff: float
    path: Path


class _Rule(NamedTuple):
    sections: Dict[str, List[str]]
    examples: List[str]
    file: Path
    line: int


def find_detection_dir(path) -> Path:
    """The hmm_detection directory: ``path`` itself, or inside an antiSMASH checkout or installation."""
    path = Path(path)
    for candidate in (path, path / "detection" / "hmm_detection", path / "antismash" / "detection" / "hmm_detection"):
        if (candidate / "cluster_rules").is_dir() and (candidate / "data" / "hmmdetails.txt").exists():
            return candidate
    raise InputError(f"{path} is not an antiSMASH hmm_detection directory (with cluster_rules/ and data/hmmdetails.txt), "
                     "nor an antiSMASH checkout or installation holding one.")


def _strip_comments(text: str) -> List[Tuple[int, str]]:
    return [(number, line.split("#", 1)[0]) for number, line in enumerate(text.splitlines(), 1)]


def _translate(words: Sequence[str]) -> str:
    """A condition from antiSMASH's spelling (and, or, not) to Domainator's (&, |, ~)."""
    text = ""
    previous = None
    for token in tokenize(" ".join(words), OPERATORS):
        token = _OPERATORS.get(token, token)
        joined = (previous is None or previous in ("(", "[", "~") or token in (")", "]", ",")
                  or (token == "(" and previous in FUNCTIONS))
        text += token if joined else " " + token
        previous = token
    return text


def _antismash_version(detection: Path) -> Optional[str]:
    """antiSMASH's version, from the main.py of the package holding ``detection`` (read, not imported)."""
    main = detection.parent.parent / "main.py"
    if main.is_file():
        match = re.search(r"^__version__\s*=\s*[\"']([^\"']+)[\"']", main.read_text(encoding="utf-8"), re.M)
        if match:
            return match.group(1)
    return None


def _sections(file: Path) -> List[Tuple[str, List[str], int]]:
    """The file as (keyword, its words, line) sections: each keyword's words run to the next keyword."""
    sections: List[Tuple[str, List[str], int]] = list()
    for number, line in _strip_comments(file.read_text(encoding="utf-8")):
        for word in line.split():
            if word in _MARKERS:
                sections.append((word, [], number))
            elif not sections:
                raise InputError(f"{file}, line {number}: {word!r} before the first RULE or DEFINE.")
            else:
                sections[-1][1].append(word)
    return sections


def _read_rules(file: Path) -> Tuple[List[_Rule], Dict[str, str]]:
    """The rules of one rule file, and its DEFINE aliases (translated)."""
    rules: List[_Rule] = list()
    define: Dict[str, str] = dict()
    rule: Optional[_Rule] = None
    for keyword, words, number in _sections(file):
        if keyword == "DEFINE":
            if len(words) < 3 or words[1] != "AS":
                raise InputError(f"{file}, line {number}: expected DEFINE <label> AS <conditions>.")
            define[words[0]] = _translate(words[2:])
            rule = None
        elif keyword == "RULE":
            rule = _Rule({"RULE": words}, [], file, number)
            rules.append(rule)
        elif rule is None:
            raise InputError(f"{file}, line {number}: {keyword} outside a RULE.")
        elif keyword in _REPEATABLE:
            rule.examples.append(words)
        elif keyword in rule.sections:
            raise InputError(f"{file}, line {number}: rule {' '.join(rule.sections['RULE'])} has two {keyword} sections.")
        else:
            rule.sections[keyword] = words
    return rules, define


def read_profiles(detection: Path, notes: List[Note]) -> Dict[str, Profile]:
    """{identifier: Profile} from data/hmmdetails.txt (identifier, description, cutoff, file)."""
    profiles = dict()
    details = detection / "data" / "hmmdetails.txt"
    for number, line in _strip_comments(details.read_text(encoding="utf-8")):
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 4:
            raise InputError(f"{details}, line {number}: expected identifier, description, cutoff and file, tab-separated.")
        identifier, _, cutoff, file = fields
        path = detection / "data" / file
        if not path.exists():
            notes.append(Note(identifier, "profile", f"dropped: hmmdetails.txt names {file}, which is missing"))
            continue
        name = None
        with open(path, encoding="utf-8") as handle:
            for profile_line in handle:
                if profile_line.startswith("NAME "):
                    name = profile_line.split()[1]
                    break
        if name is None:
            raise InputError(f"{path} is not an HMMER3 profile: it has no NAME line.")
        profiles[identifier] = Profile(identifier, name, float(cutoff), path)
    return profiles


def _example(words: Sequence[str]) -> Optional[Example]:
    """``NCBI CP006259.1 382565-499230 kirromycin`` as an Example."""
    if len(words) < 2:
        return None
    start = end = None
    note_words = list(words[3:])
    if len(words) >= 3:
        first, _, last = words[2].partition("-")
        if first.isdigit() and last.isdigit():
            start, end = int(first), int(last)
        else:
            note_words = list(words[2:])
    return Example(accession=words[1], start=start, end=end, note=" ".join(note_words) or None)


def _ids(words: Sequence[str]) -> Tuple[str, ...]:
    return tuple(part.strip() for part in " ".join(words).split(",") if part.strip())


def convert(path, *, strictness: str = DEFAULT_STRICTNESS, substitutes_file: Optional[str] = None,
            cutoff_multiplier: float = 1.0, neighborhood_multiplier: float = 1.0) -> Conversion:
    """Convert the cluster rules of an antiSMASH hmm_detection directory (or a checkout or
    installation holding one).

    strictness: which rule files to load, cumulatively, as antiSMASH's --hmmdetection-strictness.
    substitutes_file: a file of components (see read_substitutes) whose references replace those
        of the rules' components of the same name, such as find_features motifs for dynamic profiles.
    cutoff_multiplier, neighborhood_multiplier: scale every CUTOFF and NEIGHBOURHOOD, as
        antiSMASH's fungal multipliers do.
    """
    if strictness not in STRICTNESS:
        raise ValueError(f"strictness must be one of {', '.join(STRICTNESS)}, not {strictness}")
    detection = find_detection_dir(path)
    notes: List[Note] = list()
    profiles = read_profiles(detection, notes)
    substitutes = read_substitutes(substitutes_file) if substitutes_file is not None else {}
    dynamic = {file.stem.lower() for file in (detection / "dynamic_profiles").glob("*.py") if not file.stem.startswith("_")} \
        if (detection / "dynamic_profiles").is_dir() else set()

    rules: List[_Rule] = list()
    define: Dict[str, str] = dict()
    for level in STRICTNESS[:STRICTNESS.index(strictness) + 1]:
        file = detection / "cluster_rules" / f"{level}.txt"
        if not file.exists():
            raise InputError(f"{file} is missing: antiSMASH's {level} rules.")
        file_rules, file_define = _read_rules(file)
        rules.extend(file_rules)
        define.update(file_define)
    shadowed = sorted(set(define) & set(profiles))
    if shadowed:
        raise InputError(f"DEFINE labels can't also be profile identifiers: {', '.join(shadowed)}.")

    names = [" ".join(rule.sections.get("RULE", [])) for rule in rules]
    models: List[SystemModel] = list()
    used: Dict[str, Profile] = dict()
    substituted = set()
    for rule, name in zip(rules, names):
        where = f"{rule.file}, line {rule.line}"

        def section(key, required=True):
            words = rule.sections.get(key)
            if words is None and required:
                raise InputError(f"{where}: rule {name or '?'} has no {key}.")
            return words

        if len(rule.sections["RULE"]) != 1:
            raise InputError(f"{where}: RULE takes one identifier.")
        try:
            cutoff = float(" ".join(section("CUTOFF")))
            neighbourhood = float(" ".join(section("NEIGHBOURHOOD")))
        except ValueError:
            raise InputError(f"{where}: rule {name}'s CUTOFF and NEIGHBOURHOOD must be numbers (kb).") from None
        max_kb = (1000 * cutoff * cutoff_multiplier - 1) / 1000
        if max_kb < 0:
            notes.append(Note(name, "CUTOFF", f"approximated: CUTOFF {cutoff:g} as a 0 kb anchor window (CDSs that touch "
                                              "or overlap the anchor)"))
            max_kb = 0.0
        condition_text = _translate(section("CONDITIONS"))
        extenders_words = section("EXTENDERS", required=False)
        extenders_text = _translate(extenders_words) if extenders_words else None
        try:
            condition = Condition(condition_text, define)
            extenders = Condition(extenders_text, define) if extenders_text else None
        except InputError as e:
            raise InputError(f"{where}: rule {name}: {e}") from None
        components = list()
        for component in implicit_components(condition, {}, extenders):
            identifier = component.name
            if identifier in substitutes:
                references = substitutes[identifier]
                substituted.add(identifier)
            elif identifier in profiles:
                profile = profiles[identifier]
                used.setdefault(identifier, profile)
                references = (ReferenceRule(profile.name, database=DATABASE, min_score=profile.cutoff),)
            elif identifier.lower() in dynamic:
                references = ()
                notes.append(Note(name, identifier, "no references: a dynamic profile (code), so it is never found; "
                                                    "give a substitute with --substitutes_file"))
            else:
                references = ()
                notes.append(Note(name, identifier, "no references: no profile in hmmdetails.txt has this identifier"))
            components.append(Component(identifier, component.role, references))
        superiors = list()
        for superior in _ids(section("SUPERIORS", required=False) or []):
            if superior in names:
                superiors.append((DATABASE, superior))
            else:
                notes.append(Note(name, "superiors", f"dropped: {superior} is not among the {strictness} rules"))
        examples = tuple(example for example in (_example(words) for words in rule.examples) if example is not None)
        aliases = condition.aliases | (extenders.aliases if extenders is not None else frozenset())
        models.append(SystemModel(
            name=name, database=DATABASE, program="antismash", components=tuple(components),
            window=Window("anchor", max_kb=max_kb),
            condition=condition_text, extenders=extenders_text,
            define=tuple(sorted((alias, define[alias]) for alias in aliases)),
            superiors=tuple(superiors), category=" ".join(section("CATEGORY")),
            neighborhood_kb=neighbourhood * neighborhood_multiplier,
            description=" ".join(section("DESCRIPTION", required=False) or ["."]),
            accession=f"{DATABASE}:{name}", examples=examples,
            related=_ids(section("RELATED", required=False) or [])))
    for identifier in sorted(set(substitutes) - substituted):
        notes.append(Note(".", identifier, "unused substitute: no rule has a component with this name"))
    for identifier in sorted(substituted):
        notes.append(Note(".", identifier, "substituted: references from the substitutes file"))
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise InputError(f"Rules defined twice: {', '.join(duplicates)}.")

    competing = list()
    filter_file = detection / "filterhmmdetails.txt"
    if filter_file.exists():
        for _, line in _strip_comments(filter_file.read_text(encoding="utf-8")):
            group = frozenset(profiles[identifier].name for identifier in _ids([line]) if identifier in profiles)
            if len(group) > 1:
                competing.append(group)
    if not models:
        raise InputError(f"No rules in {detection / 'cluster_rules'}.")
    collection = ModelCollection(
        database=DATABASE, models=tuple(models), semantics="antismash",
        source=Source(name="antiSMASH cluster rules", url=None, release=_antismash_version(detection),
                      date=datetime.date.today().isoformat(), converter=f"domainator {__version__}",
                      license="AGPL-3.0"),
        competing=tuple(competing), competing_overlap_aa=20, reference_files=(f"{DATABASE}.hmm",))
    license_file = next((parent / name for parent in [detection, *detection.parents][:6] for name in _LICENSE_NAMES
                         if (parent / name).is_file()), None)
    notes.append(Note(".", "license", "AGPL-3.0: the converted directory is derived from antiSMASH's rules and profiles"))
    return Conversion(collection, {f"{DATABASE}.hmm": [profile.path for profile in used.values()]},
                      license_file, notes, license_notice=LICENSE_NOTICE)
