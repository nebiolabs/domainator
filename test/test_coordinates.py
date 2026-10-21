"""Coordinates of the hits domainate writes: on proteins, on both strands of a contig, on CDSs
split into parts or spanning the origin, and with a codon_start offset.

Each case embeds an HMM's consensus sequence at a known place, so the domain must cover exactly
that place: on proteins, its location and rstart/rend; on contigs, extracting the domain's
location from the contig and translating it must give the consensus back.
"""
import pyhmmer
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from domainator import DOMAIN_FEATURE_NAME
from domainator.domainate import domainate
from domainator.locations import pad_location

_CODONS = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT",
           "K": "AAA", "L": "CTT", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT",
           "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}
_LEFT = "DEKRSTNQDEKRSTNQ"  # protein flanks that don't align to the HMM
_RIGHT = "PPPPPPPPPP"


@pytest.fixture
def reference(shared_datadir):
    path = str(shared_datadir / "pdonr_hmms.hmm")
    with pyhmmer.plan7.HMMFile(path) as handle:
        hmm = next(h for h in handle if h.name == "CcdA")
    return path, hmm.consensus.upper()


def _back_translate(protein):
    return "".join(_CODONS[aa] for aa in protein) + "TAA"


def _domains(record):
    return [f for f in record.features if f.type == DOMAIN_FEATURE_NAME and f.qualifiers["name"][0] == "CcdA"]


def _search(record, path):
    out = list(domainate([record], references=[path], z=1000, evalue=1e-5, cpu=1))
    domains = _domains(out[0])
    assert len(domains) == 1
    return out[0], domains[0]


def test_protein_hit_coordinates(reference):
    path, consensus = reference
    record = SeqRecord(Seq(_LEFT + consensus + _RIGHT), id="p", name="p", annotations={"molecule_type": "protein"})
    _, domain = _search(record, path)
    # The location is 0-based and half-open; rstart and rend are 1-based and inclusive on the reference.
    assert (int(domain.location.start), int(domain.location.end)) == (len(_LEFT), len(_LEFT) + len(consensus))
    assert domain.qualifiers["rstart"] == ["1"]
    assert domain.qualifiers["rend"] == [str(len(consensus))]
    assert domain.qualifiers["rlen"] == [str(len(consensus))]


def _contig(seq, cds_location, topology="linear", codon_start=None):
    qualifiers = {} if codon_start is None else {"codon_start": [str(codon_start)]}
    return SeqRecord(Seq(seq), id="c", name="c",
                     annotations={"molecule_type": "DNA", "topology": topology},
                     features=[SeqFeature(cds_location, type="CDS", qualifiers=qualifiers)])


def _assert_domain_is_consensus(record, path, consensus):
    contig, domain = _search(record, path)
    assert str(domain.location.extract(contig.seq).translate()) == consensus
    return domain


def test_forward_strand_cds(reference):
    path, consensus = reference
    cds = _back_translate(_LEFT + consensus + _RIGHT)
    seq = "C" * 15 + cds + "C" * 20
    domain = _assert_domain_is_consensus(_contig(seq, SimpleLocation(15, 15 + len(cds), 1)), path, consensus)
    assert domain.location.strand == 1
    assert int(domain.location.start) == 15 + 3 * len(_LEFT)


def test_reverse_strand_cds(reference):
    path, consensus = reference
    cds = _back_translate(_LEFT + consensus + _RIGHT)
    seq = str(Seq("C" * 15 + cds + "C" * 20).reverse_complement())
    domain = _assert_domain_is_consensus(_contig(seq, SimpleLocation(20, 20 + len(cds), -1)), path, consensus)
    assert domain.location.strand == -1
    assert int(domain.location.end) == 20 + len(cds) - 3 * len(_LEFT)


@pytest.mark.parametrize("split", [0, 1, 2], ids=["between_codons", "after_first_base", "after_second_base"])
@pytest.mark.parametrize("strand", [1, -1])
def test_cds_split_into_parts(reference, split, strand):
    # A CDS of two parts with a 30-base intron in the middle of the domain.
    path, consensus = reference
    cds = _back_translate(_LEFT + consensus + _RIGHT)
    cut = 3 * (len(_LEFT) + 20) + split
    seq = "C" * 15 + cds[:cut] + "G" * 30 + cds[cut:] + "C" * 20
    first = SimpleLocation(15, 15 + cut, 1)
    second = SimpleLocation(15 + cut + 30, 15 + len(cds) + 30, 1)
    location = CompoundLocation([first, second])
    if strand == -1:
        seq = str(Seq(seq).reverse_complement())
        length = len(seq)
        location = CompoundLocation([SimpleLocation(length - int(p.end), length - int(p.start), -1) for p in (first, second)])
    domain = _assert_domain_is_consensus(_contig(seq, location), path, consensus)
    assert len(domain.location.parts) == 2


@pytest.mark.parametrize("strand", [1, -1])
def test_cds_spanning_the_origin(reference, strand):
    path, consensus = reference
    cds = _back_translate(_LEFT + consensus + _RIGHT)
    head = 3 * (len(_LEFT) + 30) + 1  # bases of the CDS before the origin, cutting a codon
    seq = cds[head:] + "C" * 40 + cds[:head]
    length = len(seq)
    location = CompoundLocation([SimpleLocation(length - head, length, 1), SimpleLocation(0, len(cds) - head, 1)])
    if strand == -1:
        seq = str(Seq(seq).reverse_complement())
        location = CompoundLocation([SimpleLocation(length - int(p.end), length - int(p.start), -1) for p in location.parts])
    domain = _assert_domain_is_consensus(_contig(seq, location, topology="circular"), path, consensus)
    assert len(domain.location.parts) == 2


@pytest.mark.parametrize("codon_start", [2, 3])
def test_codon_start(reference, codon_start):
    # The CDS begins codon_start - 1 bases before its first full codon, as a partial gene does.
    path, consensus = reference
    cds = "A" * (codon_start - 1) + _back_translate(_LEFT + consensus + _RIGHT)
    seq = "C" * 15 + cds + "C" * 20
    location = SimpleLocation(15, 15 + len(cds), 1)
    _assert_domain_is_consensus(_contig(seq, location, codon_start=codon_start), path, consensus)


@pytest.mark.parametrize("strand", [1, -1])
def test_padding_around_a_circular_contig_stops_at_the_feature(strand):
    # Padding that wraps all the way around a circular contig covers each base exactly once.
    record = SeqRecord(Seq("A" * 100), id="c", annotations={"molecule_type": "DNA", "topology": "circular"})
    location = SimpleLocation(40, 60, strand)
    padded = pad_location(record, location, 500, 500)
    covered = [i for part in padded.parts for i in range(int(part.start), int(part.end))]
    assert sorted(covered) == list(range(100))
