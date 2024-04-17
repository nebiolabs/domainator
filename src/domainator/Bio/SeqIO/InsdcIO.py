# Copyright 2007-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for the "genbank" and "embl" file formats.

You are expected to use this module via the Bio.SeqIO functions.
Note that internally this module calls Bio.GenBank to do the actual
parsing of GenBank, EMBL and IMGT files.

See Also:
International Nucleotide Sequence Database Collaboration
http://www.insdc.org/

GenBank
http://www.ncbi.nlm.nih.gov/Genbank/

EMBL Nucleotide Sequence Database
http://www.ebi.ac.uk/embl/

DDBJ (DNA Data Bank of Japan)
http://www.ddbj.nig.ac.jp/

IMGT (use a variant of EMBL format with longer feature indents)
http://imgt.cines.fr/download/LIGM-DB/userman_doc.html
http://imgt.cines.fr/download/LIGM-DB/ftable_doc.html
http://www.ebi.ac.uk/imgt/hla/docs/manual.html

"""
import warnings

from datetime import datetime

from domainator.Bio import BiopythonWarning
from domainator.Bio import SeqFeature
from domainator.Bio import SeqIO
from domainator.Bio.GenBank.Scanner import GenBankScanner
from domainator.Bio.Seq import UndefinedSequenceError
from domainator.Bio.Seq import UnknownSeq

from .Interfaces import _get_seq_string
from .Interfaces import SequenceIterator
from .Interfaces import SequenceWriter


# NOTE
# ====
# The "brains" for parsing GenBank, EMBL and IMGT files (and any
# other flat file variants from the INSDC in future) is in
# Bio.GenBank.Scanner (plus the _FeatureConsumer in Bio.GenBank)
# However, all the writing code is in this file.


class GenBankIterator(SequenceIterator):
    """Parser for GenBank files."""

    def __init__(self, source):
        """Break up a Genbank file into SeqRecord objects.

        Argument source is a file-like object opened in text mode or a path to a file.
        Every section from the LOCUS line to the terminating // becomes
        a single SeqRecord with associated annotation and features.

        Note that for genomes or chromosomes, there is typically only
        one record.

        This gets called internally by Bio.SeqIO for the GenBank file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("GenBank/cor6_6.gb", "gb"):
        ...     print(record.id)
        ...
        X55053.1
        X62281.1
        M81224.1
        AJ237582.1
        L31939.1
        AF297471.1

        Equivalently,

        >>> with open("GenBank/cor6_6.gb") as handle:
        ...     for record in GenBankIterator(handle):
        ...         print(record.id)
        ...
        X55053.1
        X62281.1
        M81224.1
        AJ237582.1
        L31939.1
        AF297471.1

        """
        super().__init__(source, mode="t", fmt="GenBank")

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        records = GenBankScanner(debug=0).parse_records(handle)
        return records


class GenBankCdsFeatureIterator(SequenceIterator):
    """Parser for GenBank files, creating a SeqRecord for each CDS feature."""

    def __init__(self, source):
        """Break up a Genbank file into SeqRecord objects for each CDS feature.

        Argument source is a file-like object opened in text mode or a path to a file.

        Every section from the LOCUS line to the terminating // can contain
        many CDS features.  These are returned as with the stated amino acid
        translation sequence (if given).
        """
        super().__init__(source, mode="t", fmt="GenBank")

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        return GenBankScanner(debug=0).parse_cds_features(handle)

def _insdc_feature_position_string(pos, offset=0):
    """Build a GenBank/EMBL position string (PRIVATE).

    Use offset=1 to add one to convert a start position from python counting.
    """
    if isinstance(pos, SeqFeature.ExactPosition):
        return "%i" % (pos + offset)
    elif isinstance(pos, SeqFeature.WithinPosition):
        # TODO - avoid private variables
        return "(%i.%i)" % (
            pos._left + offset,
            pos._right + offset,
        )
    elif isinstance(pos, SeqFeature.BetweenPosition):
        # TODO - avoid private variables
        return "(%i^%i)" % (
            pos._left + offset,
            pos._right + offset,
        )
    elif isinstance(pos, SeqFeature.BeforePosition):
        return "<%i" % (pos + offset)
    elif isinstance(pos, SeqFeature.AfterPosition):
        return ">%i" % (pos + offset)
    elif isinstance(pos, SeqFeature.OneOfPosition):
        return "one-of(%s)" % ",".join(
            _insdc_feature_position_string(p, offset) for p in pos.position_choices
        )
    elif isinstance(pos, SeqFeature.Position):
        raise NotImplementedError("Please report this as a bug in Biopython.")
    else:
        raise ValueError("Expected a SeqFeature position object.")

def _insdc_location_string_ignoring_strand_and_subfeatures(location, rec_length):
    if location.ref:
        ref = f"{location.ref}:"
    else:
        ref = ""
    assert not location.ref_db
    if (
        isinstance(location.start, SeqFeature.ExactPosition)
        and isinstance(location.end, SeqFeature.ExactPosition)
        and location.start == location.end
    ):
        # Special case, for 12:12 return 12^13
        # (a zero length slice, meaning the point between two letters)
        if location.end == rec_length:
            # Very special case, for a between position at the end of a
            # sequence (used on some circular genomes, Bug 3098) we have
            # N:N so return N^1
            return "%s%i^1" % (ref, rec_length)
        else:
            return "%s%i^%i" % (ref, location.end, location.end + 1)
    if (
        isinstance(location.start, SeqFeature.ExactPosition)
        and isinstance(location.end, SeqFeature.ExactPosition)
        and location.start + 1 == location.end
    ):
        # Special case, for 11:12 return 12 rather than 12..12
        # (a length one slice, meaning a single letter)
        return "%s%i" % (ref, location.end)
    elif isinstance(location.start, SeqFeature.UnknownPosition) or isinstance(
        location.end, SeqFeature.UnknownPosition
    ):
        # Special case for features from SwissProt/UniProt files
        if isinstance(location.start, SeqFeature.UnknownPosition) and isinstance(
            location.end, SeqFeature.UnknownPosition
        ):
            # warnings.warn("Feature with unknown location", BiopythonWarning)
            # return "?"
            raise ValueError("Feature with unknown location")
        elif isinstance(location.start, SeqFeature.UnknownPosition):
            # Treat the unknown start position as a BeforePosition
            return "%s<%i..%s" % (
                ref,
                location.end,
                _insdc_feature_position_string(location.end),
            )
        else:
            # Treat the unknown end position as an AfterPosition
            return "%s%s..>%i" % (
                ref,
                _insdc_feature_position_string(location.start, +1),
                location.start + 1,
            )
    else:
        # Typical case, e.g. 12..15 gets mapped to 11:15
        return (
            ref
            + _insdc_feature_position_string(location.start, +1)
            + ".."
            + _insdc_feature_position_string(location.end)
        )


def _insdc_location_string(location, rec_length):
    """Build a GenBank/EMBL location from a (Compound) SimpleLocation (PRIVATE).

    There is a choice of how to show joins on the reverse complement strand,
    GenBank used "complement(join(1,10),(20,100))" while EMBL used to use
    "join(complement(20,100),complement(1,10))" instead (but appears to have
    now adopted the GenBank convention). Notice that the order of the entries
    is reversed! This function therefore uses the first form. In this situation
    we expect the CompoundLocation and its parts to all be marked as
    strand == -1, and to be in the order 19:100 then 0:10.
    """
    try:
        parts = location.parts
        # CompoundLocation
        if location.strand == -1:
            # Special case, put complement outside the join/order/... and reverse order
            return "complement(%s(%s))" % (
                location.operator,
                ",".join(
                    _insdc_location_string_ignoring_strand_and_subfeatures(
                        p, rec_length
                    )
                    for p in parts[::-1]
                ),
            )
        else:
            return "%s(%s)" % (
                location.operator,
                ",".join(_insdc_location_string(p, rec_length) for p in parts),
            )
    except AttributeError:
        # SimpleLocation
        loc = _insdc_location_string_ignoring_strand_and_subfeatures(
            location, rec_length
        )
        if location.strand == -1:
            return f"complement({loc})"
        else:
            return loc


class _InsdcWriter(SequenceWriter):
    """Base class for GenBank and EMBL writers (PRIVATE)."""

    MAX_WIDTH = 80
    QUALIFIER_INDENT = 21
    QUALIFIER_INDENT_STR = " " * QUALIFIER_INDENT
    QUALIFIER_INDENT_TMP = "     %s                "  # 21 if %s is empty
    FTQUAL_NO_QUOTE = (
        "anticodon",
        "citation",
        "codon_start",
        "compare",
        "direction",
        "estimated_length",
        "mod_base",
        "number",
        "rpt_type",
        "rpt_unit_range",
        "tag_peptide",
        "transl_except",
        "transl_table",
    )

    def _write_feature_qualifier(self, key, value=None, quote=None):
        if key == "translation":
            MAX_WIDTH = 80
        else:
            MAX_WIDTH=10000000000000
        if value is None:
            # Value-less entry like /pseudo
            self.handle.write(f"{self.QUALIFIER_INDENT_STR}/{key}\n")
            return

        if type(value) == str:
            value = value.replace(
                '"', '""'
            )  # NCBI says escape " as "" in qualifier values

        # Quick hack with no line wrapping, may be useful for testing:
        # self.handle.write('%s/%s="%s"\n' % (self.QUALIFIER_INDENT_STR, key, value))
        if quote is None:
            # Try to mimic unwritten rules about when quotes can be left out:
            if isinstance(value, int) or key in self.FTQUAL_NO_QUOTE:
                quote = False
            else:
                quote = True
        if quote:
            line = f'{self.QUALIFIER_INDENT_STR}/{key}="{value}"'
        else:
            line = f"{self.QUALIFIER_INDENT_STR}/{key}={value}"
        if len(line) <= MAX_WIDTH:
            self.handle.write(line + "\n")
            return
        while line.lstrip():
            if len(line) <= MAX_WIDTH:
                self.handle.write(line + "\n")
                return
            # Insert line break...
            for index in range(
                min(len(line) - 1, MAX_WIDTH), self.QUALIFIER_INDENT + 1, -1
            ):
                if line[index] == " ":
                    break
            if line[index] != " ":
                # No nice place to break...
                index = MAX_WIDTH
            assert index <= MAX_WIDTH
            self.handle.write(line[:index] + "\n")
            line = self.QUALIFIER_INDENT_STR + line[index:].lstrip()

    def _wrap_location(self, location):
        """Split a feature location into lines (break at commas) (PRIVATE)."""
        # TODO - Rewrite this not to recurse!
        length = self.MAX_WIDTH - self.QUALIFIER_INDENT
        if len(location) <= length:
            return location
        index = location[:length].rfind(",")
        if index == -1:
            # No good place to split (!)
            warnings.warn(f"Couldn't split location:\n{location}", BiopythonWarning)
            return location
        return (
            location[: index + 1]
            + "\n"
            + self.QUALIFIER_INDENT_STR
            + self._wrap_location(location[index + 1 :])
        )

    def _write_feature(self, feature, record_length):
        """Write a single SeqFeature object to features table (PRIVATE)."""
        assert feature.type, feature
        location = _insdc_location_string(feature.location, record_length)
        f_type = feature.type.replace(" ", "_")
        line = (
            (self.QUALIFIER_INDENT_TMP % f_type)[: self.QUALIFIER_INDENT]
            + self._wrap_location(location)
            + "\n"
        )
        self.handle.write(line)
        # Now the qualifiers...
        # Note as of Biopython 1.69, this is an ordered-dict, don't sort it:
        for key, values in feature.qualifiers.items():
            if isinstance(values, (list, tuple)):
                for value in values:
                    self._write_feature_qualifier(key, value)
            else:
                # String, int, etc - or None for a /pseudo tpy entry
                self._write_feature_qualifier(key, values)

    @staticmethod
    def _get_annotation_str(record, key, default=".", just_first=False):
        """Get an annotation dictionary entry (as a string) (PRIVATE).

        Some entries are lists, in which case if just_first=True the first entry
        is returned.  If just_first=False (default) this verifies there is only
        one entry before returning it.
        """
        try:
            answer = record.annotations[key]
        except KeyError:
            return default
        if isinstance(answer, list):
            if not just_first:
                assert len(answer) == 1
            return str(answer[0])
        else:
            return str(answer)

    @staticmethod
    def _split_multi_line(text, max_len):
        """Return a list of strings (PRIVATE).

        Any single words which are too long get returned as a whole line
        (e.g. URLs) without an exception or warning.
        """
        # TODO - Do the line splitting while preserving white space?
        text = text.strip()
        if len(text) <= max_len:
            return [text]

        words = text.split()
        text = ""
        while words and len(text) + 1 + len(words[0]) <= max_len:
            text += " " + words.pop(0)
            text = text.strip()
        # assert len(text) <= max_len
        answer = [text]
        while words:
            text = words.pop(0)
            while words and len(text) + 1 + len(words[0]) <= max_len:
                text += " " + words.pop(0)
                text = text.strip()
            # assert len(text) <= max_len
            answer.append(text)
        assert not words
        return answer

    def _split_contig(self, record, max_len):
        """Return a list of strings, splits on commas (PRIVATE)."""
        # TODO - Merge this with _write_multi_line method?
        # It would need the addition of the comma splitting logic...
        # are there any other cases where that would be sensible?
        contig = record.annotations.get("contig", "")
        if isinstance(contig, (list, tuple)):
            contig = "".join(contig)
        contig = self.clean(contig)
        answer = []
        while contig:
            if len(contig) > max_len:
                # Split lines at the commas
                pos = contig[: max_len - 1].rfind(",")
                if pos == -1:
                    raise ValueError("Could not break up CONTIG")
                text, contig = contig[: pos + 1], contig[pos + 1 :]
            else:
                text, contig = contig, ""
            answer.append(text)
        return answer


class GenBankWriter(_InsdcWriter):
    """GenBank writer."""

    HEADER_WIDTH = 12
    QUALIFIER_INDENT = 21
    STRUCTURED_COMMENT_START = "-START##"
    STRUCTURED_COMMENT_END = "-END##"
    STRUCTURED_COMMENT_DELIM = " :: "
    LETTERS_PER_LINE = 60
    SEQUENCE_INDENT = 9

    def _write_single_line(self, tag, text):
        """Write single line in each GenBank record (PRIVATE).

        Used in the 'header' of each GenBank record.
        """
        assert len(tag) < self.HEADER_WIDTH
        if len(text) > self.MAX_WIDTH - self.HEADER_WIDTH:
            if tag:
                warnings.warn(
                    f"Annotation {text!r} too long for {tag!r} line", BiopythonWarning
                )
            else:
                # Can't give such a precise warning
                warnings.warn(f"Annotation {text!r} too long", BiopythonWarning)
        self.handle.write(
            "%s%s\n" % (tag.ljust(self.HEADER_WIDTH), text.replace("\n", " "))
        )

    def _write_multi_line(self, tag, text):
        """Write multiple lines in each GenBank record (PRIVATE).

        Used in the 'header' of each GenBank record.
        """
        # TODO - Do the line splitting while preserving white space?
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        lines = self._split_multi_line(text, max_len)
        self._write_single_line(tag, lines[0])
        for line in lines[1:]:
            self._write_single_line("", line)

    def _write_multi_entries(self, tag, text_list):
        # used for DBLINK and any similar later line types.
        # If the list of strings is empty, nothing is written.
        for i, text in enumerate(text_list):
            if i == 0:
                self._write_single_line(tag, text)
            else:
                self._write_single_line("", text)

    @staticmethod
    def _get_date(record):
        default = "01-JAN-1980"
        try:
            date = record.annotations["date"]
        except KeyError:
            return default
        # Cope with a list of one string:
        if isinstance(date, list) and len(date) == 1:
            date = date[0]
        if isinstance(date, datetime):
            date = date.strftime("%d-%b-%Y").upper()

        months = [
            "JAN",
            "FEB",
            "MAR",
            "APR",
            "MAY",
            "JUN",
            "JUL",
            "AUG",
            "SEP",
            "OCT",
            "NOV",
            "DEC",
        ]
        if not isinstance(date, str) or len(date) != 11:
            return default
        try:
            datetime(int(date[-4:]), months.index(date[3:6]) + 1, int(date[0:2]))
        except ValueError:
            date = default
        return date

    @staticmethod
    def _get_data_division(record):
        try:
            division = record.annotations["data_file_division"]
        except KeyError:
            division = "UNK"
        if division in [
            "PRI",
            "ROD",
            "MAM",
            "VRT",
            "INV",
            "PLN",
            "BCT",
            "VRL",
            "PHG",
            "SYN",
            "UNA",
            "EST",
            "PAT",
            "STS",
            "GSS",
            "HTG",
            "HTC",
            "ENV",
            "CON",
            "TSA",
        ]:
            # Good, already GenBank style
            #    PRI - primate sequences
            #    ROD - rodent sequences
            #    MAM - other mammalian sequences
            #    VRT - other vertebrate sequences
            #    INV - invertebrate sequences
            #    PLN - plant, fungal, and algal sequences
            #    BCT - bacterial sequences [plus archaea]
            #    VRL - viral sequences
            #    PHG - bacteriophage sequences
            #    SYN - synthetic sequences
            #    UNA - unannotated sequences
            #    EST - EST sequences (expressed sequence tags)
            #    PAT - patent sequences
            #    STS - STS sequences (sequence tagged sites)
            #    GSS - GSS sequences (genome survey sequences)
            #    HTG - HTGS sequences (high throughput genomic sequences)
            #    HTC - HTC sequences (high throughput cDNA sequences)
            #    ENV - Environmental sampling sequences
            #    CON - Constructed sequences
            #    TSA - Transcriptome Shotgun Assembly
            #
            # (plus UNK for unknown)
            pass
        else:
            # See if this is in EMBL style:
            #    Division                 Code
            #    -----------------        ----
            #    Bacteriophage            PHG - common
            #    Environmental Sample     ENV - common
            #    Fungal                   FUN - map to PLN (plants + fungal)
            #    Human                    HUM - map to PRI (primates)
            #    Invertebrate             INV - common
            #    Other Mammal             MAM - common
            #    Other Vertebrate         VRT - common
            #    Mus musculus             MUS - map to ROD (rodent)
            #    Plant                    PLN - common
            #    Prokaryote               PRO - map to BCT (poor name)
            #    Other Rodent             ROD - common
            #    Synthetic                SYN - common
            #    Transgenic               TGN - ??? map to SYN ???
            #    Unclassified             UNC - map to UNK
            #    Viral                    VRL - common
            #
            # (plus XXX for submitting which we can map to UNK)
            embl_to_gbk = {
                "FUN": "PLN",
                "HUM": "PRI",
                "MUS": "ROD",
                "PRO": "BCT",
                "UNC": "UNK",
                "XXX": "UNK",
            }
            try:
                division = embl_to_gbk[division]
            except KeyError:
                division = "UNK"
        assert len(division) == 3
        return division

    def _get_topology(self, record):
        """Set the topology to 'circular', 'linear' if defined (PRIVATE)."""
        max_topology_len = len("circular")

        topology = self._get_annotation_str(record, "topology", default="")
        if topology and len(topology) <= max_topology_len:
            return topology.ljust(max_topology_len)
        else:
            return " " * max_topology_len

    def _write_the_first_line(self, record):
        """Write the LOCUS line (PRIVATE)."""
        locus = record.name
        if not locus or locus == "<unknown name>":
            locus = record.id
        if not locus or locus == "<unknown id>":
            locus = self._get_annotation_str(record, "accession", just_first=True)
        if len(locus) > 16:
            if len(locus) + 1 + len(str(len(record))) > 28:
                # Locus name and record length to long to squeeze in.
                # Per updated GenBank standard (Dec 15, 2018) 229.0
                # the Locus identifier can be any length, and a space
                # is added after the identifier to keep the identifier
                # and length fields separated
                # warnings.warn(
                #     "Increasing length of locus line to allow "
                #     "long name. This will result in fields that "
                #     "are not in usual positions.",
                #     BiopythonWarning,
                # )
                pass

        if len(locus.split()) > 1:
            raise ValueError(f"Invalid whitespace in {locus!r} for LOCUS line")
        if len(record) > 99999999999:
            # As of the GenBank release notes 229.0, the locus line can be
            # any length. However, long locus lines may not be compatible
            # with all software.
            # warnings.warn(
            #     "The sequence length is very long. The LOCUS "
            #     "line will be increased in length to compensate. "
            #     "This may cause unexpected behavior.",
            #     BiopythonWarning,
            # )
            pass

        # Get the molecule type
        mol_type = self._get_annotation_str(record, "molecule_type", None)
        if mol_type is None:
            raise ValueError("missing molecule_type in annotations")
        if mol_type and len(mol_type) > 7:
            # Deal with common cases from EMBL to GenBank
            mol_type = mol_type.replace("unassigned ", "").replace("genomic ", "")
            if len(mol_type) > 7:
                warnings.warn(f"Molecule type {mol_type!r} too long", BiopythonWarning)
                mol_type = "DNA"
        if mol_type in ["protein", "PROTEIN"]:
            mol_type = ""

        if mol_type == "":
            units = "aa"
        else:
            units = "bp"

        topology = self._get_topology(record)

        division = self._get_data_division(record)

        # Accommodate longer header, with long accessions and lengths
        if len(locus) > 16 and len(str(len(record))) > (11 - (len(locus) - 16)):
            name_length = locus + " " + str(len(record))

        # This is the older, standard 80 position header
        else:
            name_length = str(len(record)).rjust(28)
            name_length = locus + name_length[len(locus) :]
            assert len(name_length) == 28, name_length
            assert " " in name_length, name_length

        assert len(units) == 2
        assert len(division) == 3
        line = "LOCUS       %s %s    %s %s %s %s\n" % (
            name_length,
            units,
            mol_type.ljust(7),
            topology,
            division,
            self._get_date(record),
        )
        # Extra long header
        if len(line) > 80:
            splitline = line.split()
            if splitline[3] not in ["bp", "aa"]:
                raise ValueError(
                    "LOCUS line does not contain size units at "
                    "expected position:\n" + line
                )

            if not (
                splitline[3].strip() == "aa"
                or "DNA" in splitline[4].strip().upper()
                or "RNA" in splitline[4].strip().upper()
            ):
                raise ValueError(
                    "LOCUS line does not contain valid "
                    "sequence type (DNA, RNA, ...):\n" + line
                )

            self.handle.write(line)

        # 80 position header
        else:
            assert len(line) == 79 + 1, repr(line)  # plus one for new line

            # We're bending the rules to allow an identifier over 16 characters
            # if we can steal spaces from the length field:
            # assert line[12:28].rstrip() == locus, \
            #     'LOCUS line does not contain the locus at the expected position:\n' + line
            # assert line[28:29] == " "
            # assert line[29:40].lstrip() == str(len(record)), \
            #     'LOCUS line does not contain the length at the expected position:\n' + line
            assert line[12:40].split() == [locus, str(len(record))], line

            # Tests copied from Bio.GenBank.Scanner
            if line[40:44] not in [" bp ", " aa "]:
                raise ValueError(
                    "LOCUS line does not contain size units at "
                    "expected position:\n" + line
                )
            if line[44:47] not in ["   ", "ss-", "ds-", "ms-"]:
                raise ValueError(
                    "LOCUS line does not have valid strand "
                    "type (Single stranded, ...):\n" + line
                )
            if not (
                line[47:54].strip() == ""
                or "DNA" in line[47:54].strip().upper()
                or "RNA" in line[47:54].strip().upper()
            ):
                raise ValueError(
                    "LOCUS line does not contain valid "
                    "sequence type (DNA, RNA, ...):\n" + line
                )
            if line[54:55] != " ":
                raise ValueError(
                    "LOCUS line does not contain space at position 55:\n" + line
                )
            if line[55:63].strip() not in ["", "linear", "circular"]:
                raise ValueError(
                    "LOCUS line does not contain valid "
                    "entry (linear, circular, ...):\n" + line
                )
            if line[63:64] != " ":
                raise ValueError(
                    "LOCUS line does not contain space at position 64:\n" + line
                )
            if line[67:68] != " ":
                raise ValueError(
                    "LOCUS line does not contain space at position 68:\n" + line
                )
            if line[70:71] != "-":
                raise ValueError(
                    "LOCUS line does not contain - at position 71 in date:\n" + line
                )
            if line[74:75] != "-":
                raise ValueError(
                    "LOCUS line does not contain - at position 75 in date:\n" + line
                )

            self.handle.write(line)

    def _write_references(self, record):
        number = 0
        for ref in record.annotations["references"]:
            if not isinstance(ref, SeqFeature.Reference):
                continue
            number += 1
            data = str(number)
            # TODO - support more complex record reference locations?
            if ref.location and len(ref.location) == 1:
                molecule_type = record.annotations.get("molecule_type")
                if molecule_type and "protein" in molecule_type:
                    units = "residues"
                else:
                    units = "bases"
                data += "  (%s %i to %i)" % (
                    units,
                    ref.location[0].start + 1,
                    ref.location[0].end,
                )
            self._write_single_line("REFERENCE", data)
            if ref.authors:
                # We store the AUTHORS data as a single string
                self._write_multi_line("  AUTHORS", ref.authors)
            if ref.consrtm:
                # We store the consortium as a single string
                self._write_multi_line("  CONSRTM", ref.consrtm)
            if ref.title:
                # We store the title as a single string
                self._write_multi_line("  TITLE", ref.title)
            if ref.journal:
                # We store this as a single string - holds the journal name,
                # volume, year, and page numbers of the citation
                self._write_multi_line("  JOURNAL", ref.journal)
            if ref.medline_id:
                # This line type is obsolete and was removed from the GenBank
                # flatfile format in April 2005. Should we write it?
                # Note this has a two space indent:
                self._write_multi_line("  MEDLINE", ref.medline_id)
            if ref.pubmed_id:
                # Note this has a THREE space indent:
                self._write_multi_line("   PUBMED", ref.pubmed_id)
            if ref.comment:
                self._write_multi_line("  REMARK", ref.comment)

    def _write_comment(self, record):
        # This is a bit complicated due to the range of possible
        # ways people might have done their annotation...
        # Currently the parser uses a single string with newlines.
        # A list of lines is also reasonable.
        # A single (long) string is perhaps the most natural of all.
        # This means we may need to deal with line wrapping.
        lines = []
        if "structured_comment" in record.annotations:
            comment = record.annotations["structured_comment"]
            # Find max length of keys for equal padded printing
            padding = 0
            for key, data in comment.items():
                for subkey, subdata in data.items():
                    padding = len(subkey) if len(subkey) > padding else padding
            # Construct output
            for key, data in comment.items():
                lines.append(f"##{key}{self.STRUCTURED_COMMENT_START}")
                for subkey, subdata in data.items():
                    spaces = " " * (padding - len(subkey))
                    lines.append(
                        f"{subkey}{spaces}{self.STRUCTURED_COMMENT_DELIM}{subdata}"
                    )
                lines.append(f"##{key}{self.STRUCTURED_COMMENT_END}")
        if "comment" in record.annotations:
            comment = record.annotations["comment"]
            if isinstance(comment, str):
                lines += comment.split("\n")
            elif isinstance(comment, (list, tuple)):
                lines += list(comment)
            else:
                raise ValueError("Could not understand comment annotation")
        self._write_multi_line("COMMENT", lines[0])
        for line in lines[1:]:
            self._write_multi_line("", line)

    def _write_contig(self, record):
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        lines = self._split_contig(record, max_len)
        self._write_single_line("CONTIG", lines[0])
        for text in lines[1:]:
            self._write_single_line("", text)

    def _write_sequence(self, record):
        # Loosely based on code from Howard Salis
        # TODO - Force lower case?

        try:
            data = _get_seq_string(record)
        except UndefinedSequenceError:
            # We have already recorded the length, and there is no need
            # to record a long sequence of NNNNNNN...NNN or whatever.
            if "contig" in record.annotations:
                self._write_contig(record)
            else:
                self.handle.write("ORIGIN\n")
            return

        # Catches sequence being None:
        data = data.lower()
        seq_len = len(data)
        self.handle.write("ORIGIN\n")
        for line_number in range(0, seq_len, self.LETTERS_PER_LINE):
            self.handle.write(str(line_number + 1).rjust(self.SEQUENCE_INDENT))
            for words in range(
                line_number, min(line_number + self.LETTERS_PER_LINE, seq_len), 10
            ):
                self.handle.write(f" {data[words:words + 10]}")
            self.handle.write("\n")

    def write_record(self, record):
        """Write a single record to the output file."""
        handle = self.handle
        self._write_the_first_line(record)

        default = record.id
        if default.count(".") == 1 and default[default.index(".") + 1 :].isdigit():
            # Good, looks like accession.version and not something
            # else like identifier.start-end
            default = record.id.split(".", 1)[0]
        accession = self._get_annotation_str(
            record, "accession", default, just_first=True
        )
        acc_with_version = accession
        if record.id.startswith(accession + "."):
            try:
                acc_with_version = "%s.%i" % (
                    accession,
                    int(record.id.split(".", 1)[1]),
                )
            except ValueError:
                pass
        gi = self._get_annotation_str(record, "gi", just_first=True)

        descr = record.description
        if descr == "<unknown description>":
            descr = ""  # Trailing dot will be added later

        # The DEFINITION field must end with a period
        # see ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt [3.4.5]
        # and discussion https://github.com/biopython/biopython/pull/616
        # So let's add a period
        descr += "."
        self._write_multi_line("DEFINITION", descr)

        self._write_single_line("ACCESSION", accession)
        if gi != ".":
            self._write_single_line("VERSION", f"{acc_with_version}  GI:{gi}")
        else:
            self._write_single_line("VERSION", acc_with_version)

        # The NCBI initially expected two types of link,
        # e.g. "Project:28471" and "Trace Assembly Archive:123456"
        #
        # This changed and at some point the formatting switched to
        # include a space after the colon, e.g.
        #
        # LOCUS       NC_000011               1606 bp    DNA     linear   CON 06-JUN-2016
        # DEFINITION  Homo sapiens chromosome 11, GRCh38.p7 Primary Assembly.
        # ACCESSION   NC_000011 REGION: complement(5225466..5227071) GPC_000001303
        # VERSION     NC_000011.10  GI:568815587
        # DBLINK      BioProject: PRJNA168
        #             Assembly: GCF_000001405.33
        # ...
        #
        # Or,
        #
        # LOCUS       JU120277                1044 bp    mRNA    linear   TSA 27-NOV-2012
        # DEFINITION  TSA: Tupaia chinensis tbc000002.Tuchadli mRNA sequence.
        # ACCESSION   JU120277
        # VERSION     JU120277.1  GI:379775257
        # DBLINK      BioProject: PRJNA87013
        #             Sequence Read Archive: SRR433859
        # ...
        dbxrefs_with_space = []
        for x in record.dbxrefs:
            if ": " not in x:
                x = x.replace(":", ": ")
            dbxrefs_with_space.append(x)
        self._write_multi_entries("DBLINK", dbxrefs_with_space)
        del dbxrefs_with_space

        try:
            # List of strings
            # Keywords should be given separated with semi colons,
            keywords = "; ".join(record.annotations["keywords"])
            # with a trailing period:
            if not keywords.endswith("."):
                keywords += "."
        except KeyError:
            # If no keywords, there should be just a period:
            keywords = "."
        self._write_multi_line("KEYWORDS", keywords)

        if "segment" in record.annotations:
            # Deal with SEGMENT line found only in segmented records,
            # e.g. AH000819
            segment = record.annotations["segment"]
            if isinstance(segment, list):
                assert len(segment) == 1, segment
                segment = segment[0]
            self._write_single_line("SEGMENT", segment)

        self._write_multi_line("SOURCE", self._get_annotation_str(record, "source"))
        # The ORGANISM line MUST be a single line, as any continuation is the taxonomy
        org = self._get_annotation_str(record, "organism")
        if len(org) > self.MAX_WIDTH - self.HEADER_WIDTH:
            org = org[: self.MAX_WIDTH - self.HEADER_WIDTH - 4] + "..."
        self._write_single_line("  ORGANISM", org)
        try:
            # List of strings
            # Taxonomy should be given separated with semi colons,
            taxonomy = "; ".join(record.annotations["taxonomy"])
            # with a trailing period:
            if not taxonomy.endswith("."):
                taxonomy += "."
        except KeyError:
            taxonomy = "."
        self._write_multi_line("", taxonomy)

        if "db_source" in record.annotations:
            # Hack around the issue of BioSQL loading a list for the db_source
            db_source = record.annotations["db_source"]
            if isinstance(db_source, list):
                db_source = db_source[0]
            self._write_single_line("DBSOURCE", db_source)

        if "references" in record.annotations:
            self._write_references(record)

        if (
            "comment" in record.annotations
            or "structured_comment" in record.annotations
        ):
            self._write_comment(record)

        handle.write("FEATURES             Location/Qualifiers\n")
        rec_length = len(record)
        for feature in record.features:
            self._write_feature(feature, rec_length)
        self._write_sequence(record)
        handle.write("//\n")


def _genbank_convert_fasta(in_file, out_file):
    """Fast GenBank to FASTA (PRIVATE)."""
    # We don't need to parse the features...
    records = GenBankScanner().parse_records(in_file, do_features=False)
    return SeqIO.write(records, out_file, "fasta")

# if __name__ == "__main__":
#     from Bio._utils import run_doctest

#     run_doctest(verbose=0)
