"""Checking that system models still find their own examples (format_systems --check_examples_out).

A model's ``examples`` are known occurrences: an accession, and optionally a range. Each one is
downloaded from NCBI (nuccore, as GenBank with its features), annotated with the model's profiles
by domainate --system_models (calling genes with pyrodigal if the record has no CDSs), and
the model passes when it finds at least one system there. antiSMASH's EXAMPLE lines serve the
same purpose.
"""
from dataclasses import replace
import os
import tempfile
import time
from typing import Callable, Iterator, List, NamedTuple, Optional, Sequence

from Bio.SeqRecord import SeqRecord

from domainator import SYSTEM_FEATURE_NAME, __version__
from domainator.errors import DomainatorError
from domainator.systems.model import Example, ModelCollection, SystemModel

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_REQUEST_INTERVAL = 0.4  # NCBI allows 3 requests a second without an API key
REPORT_COLUMNS = ("database", "system", "accession", "start", "end", "note", "found", "systems", "best_score", "error")


class ExampleCheck(NamedTuple):
    model: SystemModel
    example: Example
    found: bool
    systems: int
    """How many systems of the model were found in the example."""
    best_score: Optional[float]
    error: Optional[str]
    """Why the example couldn't be checked (a download failure, ...), or None."""

    def row(self) -> List[str]:
        example = self.example
        return [self.model.database, self.model.name, example.accession,
                "" if example.start is None else str(example.start), "" if example.end is None else str(example.end),
                example.note or "", "yes" if self.found else "no", str(self.systems),
                "" if self.best_score is None else f"{self.best_score:.1f}", self.error or ""]


class NcbiFetcher:
    """Downloads examples from NCBI nuccore with efetch, one request at a time."""

    def __init__(self, api_key: Optional[str] = None):
        import requests
        self.session = requests.Session()
        self.session.headers["User-Agent"] = f"domainator/{__version__}"
        self.api_key = api_key
        self._last = 0.0

    def __call__(self, example: Example) -> List[SeqRecord]:
        from domainator.seqio import parse_seqfiles
        params = {"db": "nuccore", "id": example.accession, "rettype": "gbwithparts", "retmode": "text",
                  "tool": "domainator"}
        if example.start is not None and example.end is not None:
            # examples give 0- or 1-based starts; a base more on the left does no harm
            params["seq_start"] = str(max(example.start, 1))
            params["seq_stop"] = str(example.end)
        if self.api_key:
            params["api_key"] = self.api_key
        wait = self._last + _REQUEST_INTERVAL - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        response = self.session.get(EFETCH_URL, params=params, timeout=120)
        self._last = time.monotonic()
        response.raise_for_status()
        if not response.text.lstrip().startswith("LOCUS"):
            raise DomainatorError(f"NCBI returned no GenBank record for {example.accession}.")
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "example.gb")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(response.text)
            return list(parse_seqfiles([path]))


def check_examples(collections: Sequence[ModelCollection], *,
                   fetch: Optional[Callable[[Example], List[SeqRecord]]] = None, cpu: int = 0) -> Iterator[ExampleCheck]:
    """Check every example of every model in ``collections``, in order.

    fetch: downloads an example's records (default: from NCBI). An example that can't be fetched
    or annotated is reported with its error, not raised.
    """
    from domainator.domainate import domainate

    fetch = fetch or NcbiFetcher()
    for collection in collections:
        for model in collection.models:
            if not model.examples:
                continue
            single = replace(collection, models=(model,))
            for example in model.examples:
                try:
                    records = fetch(example)
                    annotated = list(domainate(records, references=None, z=1000, system_models=[single],
                                               gene_call="unannotated", cpu=cpu))
                except Exception as e:  # report it, and go on to the next example
                    yield ExampleCheck(model, example, False, 0, None, f"{type(e).__name__}: {e}".replace("\n", " "))
                    continue
                scores = [float(feature.qualifiers["score"][0]) for record in annotated for feature in record.features
                          if feature.type == SYSTEM_FEATURE_NAME and feature.qualifiers["name"][0] == model.name
                          and feature.qualifiers["database"][0] == model.database]
                yield ExampleCheck(model, example, bool(scores), len(scores), max(scores, default=None), None)
