"""GenBank reading and writing on upstream Biopython, with the few format differences Domainator needs.

Domainator reads and writes GenBank through these subclasses of Biopython's scanner, consumer
and writer instead of ``Bio.SeqIO``:

- **Long-name GenPept LOCUS lines.** Biopython misreads a protein LOCUS line whose name is too
  long for its column (it loses the molecule type and date, and then can't write the record).
- **Qualifier values are not wrapped** (except ``/translation``). Biopython wraps long values at
  80 columns, so a long value with no spaces comes back with spaces in it.
- **No NCBI quote-escape check** on qualifier values, which is most of Biopython's parse time.

The subclasses use a few private Biopython names (``_FeatureConsumer``, ``_feed_first_line``,
``_write_feature_qualifier``), so the supported Biopython versions are pinned in pyproject.toml,
and test_genbank.py round-trips every GenBank fixture.

Biopython's parser warns about every quirk of real-world files; those warnings are ignored while
a record is read or written here, without changing the warning filters anywhere else.
"""
import warnings
from typing import Iterable, Iterator

from Bio import BiopythonParserWarning, BiopythonWarning
from Bio.GenBank import _FeatureConsumer
from Bio.GenBank.Scanner import GenBankScanner
from Bio.GenBank.utils import FeatureValueCleaner
from Bio.SeqIO.InsdcIO import GenBankWriter
from Bio.SeqRecord import SeqRecord


class _Consumer(_FeatureConsumer):
    def feature_qualifier(self, key, value):
        if value is None:
            if key not in self._cur_feature.qualifiers:
                self._cur_feature.qualifiers[key] = [""]
            return
        if len(value) > 1 and value[0] == '"' and value[-1] == '"':
            value = value[1:-1]
        value = value.replace('""', '"')
        if key == "translation":
            value = "".join(value.split())
        qualifiers = self._cur_feature.qualifiers
        if key in qualifiers:
            qualifiers[key].append(value)
        else:
            qualifiers[key] = [value]


class _Scanner(GenBankScanner):
    def parse(self, handle, do_features=True):
        consumer = _Consumer(use_fuzziness=1, feature_cleaner=FeatureValueCleaner())
        return consumer.data if self.feed(handle, consumer, do_features) else None

    def _feed_first_line(self, consumer, line):
        tokens = line.split()
        if len(tokens) == 6 and tokens[3] == "aa" and line[29:33] != " aa ":
            # A GenPept LOCUS line whose long name pushed the fields out of their columns.
            consumer.locus(tokens[1])
            consumer.size(tokens[2])
            consumer.residue_type("PROTEIN")
            consumer.data_file_division(tokens[4])
            consumer.date(tokens[5])
            return
        return super()._feed_first_line(consumer, line)


class _Writer(GenBankWriter):
    def _write_feature_qualifier(self, key, value=None, quote=None):
        if key == "translation":
            return super()._write_feature_qualifier(key, value, quote)
        if value is None:
            self.handle.write(f"{self.QUALIFIER_INDENT_STR}/{key}\n")
            return
        if isinstance(value, str):
            value = value.replace('"', '""')
        if quote is None:
            quote = not (isinstance(value, int) or key in self.FTQUAL_NO_QUOTE)
        if quote:
            self.handle.write(f'{self.QUALIFIER_INDENT_STR}/{key}="{value}"\n')
        else:
            self.handle.write(f"{self.QUALIFIER_INDENT_STR}/{key}={value}\n")


def parse_genbank(handle) -> Iterator[SeqRecord]:
    """Yield the records of an open GenBank text handle."""
    records = _Scanner(debug=0).parse_records(handle)
    while True:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            record = next(records, None)
        if record is None:
            return
        yield record


def write_genbank_records(records: Iterable[SeqRecord], handle) -> None:
    """Write records to an open GenBank text handle (no other changes; see seqio.write_genbank)."""
    writer = _Writer(handle)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        for record in records:
            writer.write_record(record)
