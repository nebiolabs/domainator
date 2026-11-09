"""Synthetic contigs with hand-written Domainator features, for the systems tests.

A layout is a list of slots, one per CDS, in contig order:

- ``""``: a CDS with no hits;
- ``"A"``, ``"A+B"``: a CDS with a hit to A (and B: a fusion);
- ``"A:80"``: a hit with bit score 80 (default 50);
- a leading ``-``: the CDS is on the minus strand;
- ``"nc:R"``: no CDS, but a contig-level nucleotide hit to R (``cds_id="."``), 200 bp long.

CDSs are 900 bp with 100 bp between slots, so CDS i (0-based) starts at 1000 * i.
"""
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

GENE_LENGTH = 900
GAP = 100
NC_LENGTH = 200


def hit_feature(location, name, cds_id, *, score=50.0, evalue=1e-10, database="db", program="hmmsearch"):
    return SeqFeature(location, type="Domainator", qualifiers={
        "program": [program], "database": [database], "description": ["."], "accession": ["."],
        "evalue": [f"{evalue:.1e}"], "score": [f"{score:.1f}"], "name": [name], "identity": ["0.0"],
        "cds_id": [cds_id], "rstart": ["1"], "rend": ["100"], "rlen": ["100"]})


def _hits(spec):
    out = list()
    for part in spec.split("+"):
        name, _, score = part.partition(":")
        out.append((name, float(score) if score else 50.0))
    return out


def contig(layout, *, topology="linear", record_id="contig", evalue=1e-10, database="db", tail=GAP):
    """A record built from ``layout`` (see the module docstring)."""
    features = list()
    position = 0
    for i, slot in enumerate(layout):
        strand = -1 if slot.startswith("-") else 1
        body = slot.lstrip("-")
        if body.startswith("nc:"):
            location = FeatureLocation(position, position + NC_LENGTH, strand)
            for name, score in _hits(body[3:]):
                features.append(hit_feature(location, name, ".", score=score, evalue=evalue, database=database, program="nhmmer"))
            position += NC_LENGTH + GAP
            continue
        location = FeatureLocation(position, position + GENE_LENGTH, strand)
        cds_id = f"cds{i}"
        features.append(SeqFeature(location, type="CDS", qualifiers={"cds_id": [cds_id], "locus_tag": [f"L{i}"]}))
        if body:
            for name, score in _hits(body):
                features.append(hit_feature(location, name, cds_id, score=score, evalue=evalue, database=database))
        position += GENE_LENGTH + GAP
    record = SeqRecord(Seq("A" * (position - GAP + tail)), id=record_id, name=record_id, description="synthetic contig")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = topology
    record.features = features
    return record


def layout(positions, length=12, strand=""):
    """A layout of ``length`` empty CDSs, with the domains at the given 1-based positions."""
    slots = [strand] * length
    for name, position in positions.items():
        slots[position - 1] = strand + name
    return slots


def member_ids(hit):
    """The cds_ids of a SystemHit's members, in order, without repeats."""
    out = list()
    for member in hit.members:
        if member.cds_id not in out:
            out.append(member.cds_id)
    return out
