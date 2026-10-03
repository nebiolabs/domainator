from domainator.Bio import SeqIO

def compare_files(f1,f2, skip_lines=0):
    with open(f1,"r") as newfile, open(f2, "r") as oldfile:
        for x in range(skip_lines):
            newfile.readline()
            oldfile.readline()
        assert newfile.read() == oldfile.read()

def compare_iterables(i1, i2):
    assert all(a == b for a,b in zip(i1, i2))

def compare_seqfiles(gb1, gb2, format="genbank", skip_attrs={}, skip_qualifiers={}):
    recs1 = list(SeqIO.parse(gb1, format))
    recs2 = list(SeqIO.parse(gb2, format))
    assert len(recs1) == len(recs2)
    for i in range(len(recs1)):
        compare_seqrecords(recs1[i], recs2[i], skip_attrs=skip_attrs, skip_qualifiers=skip_qualifiers)

def compare_seqrecords(rec1, rec2, skip_attrs={}, skip_qualifiers={}):
    attrs = {"seq", "id", "description", "name"}
    skip_attrs = set(skip_attrs)
    skip_qualifiers = set(skip_qualifiers)
    attrs = attrs.difference(skip_attrs)

    for attr in attrs:
        try:
            assert getattr(rec1, attr) == getattr(rec2, attr)
        except AssertionError as e:
            e.args += (attr, rec1, rec2)
            raise
        

    assert rec1.letter_annotations == rec2.letter_annotations
    for k in rec1.letter_annotations:
        assert rec1.letter_annotations[k] == rec2.letter_annotations[k]
    for k in rec1.annotations:
        if k != "date":
            assert rec1.annotations[k] == rec2.annotations[k]
    assert len(rec1.features) == len(rec2.features)


    for i in range(len(rec1.features)):
        feature1 = rec1.features[i]
        feature2 = rec2.features[i]

        for qualifier in feature1.qualifiers:
            if qualifier in skip_qualifiers:
                continue
            try:
                assert feature1.qualifiers[qualifier] == feature2.qualifiers[qualifier], f"qualifiers not equal in: {rec1}, {rec2}"
            except:
                #print(f"{rec1}, {rec2}")
                print(f"{feature1}, {feature2}")
                print(f"{feature1.qualifiers[qualifier]}, {feature2.qualifiers[qualifier]}")
                raise


def gzip_file(src, dst):
    """Write src to dst as plain gzip. Returns str(dst)."""
    import gzip as _gzip
    with open(src, "rb") as fh, _gzip.open(str(dst), "wb") as w:
        w.write(fh.read())
    return str(dst)


def bgzip_file(src, dst):
    """Write src to dst as BGZF (block gzip). Returns str(dst)."""
    from domainator.Bio import bgzf as _bgzf
    with open(src, "rb") as fh, _bgzf.BgzfWriter(str(dst)) as w:
        w.write(fh.read())
    return str(dst)


# --- partial CDS fixture: JABFVH010000506_extraction.gb has a complete CDS (HOV79_30120) and a 5'-partial CDS (HOV79_30125) ---

PARTIAL_FIXTURE = "JABFVH010000506_extraction.gb"
PARTIAL_CDS_LOCUS_TAG = "HOV79_30125"
COMPLETE_CDS_LOCUS_TAG = "HOV79_30120"

def write_partial_fixture_hmms(shared_datadir, path):
    """Writes single-sequence HMMs built from the translations of both CDSs in the partial fixture, named by locus_tag."""
    import pyhmmer
    record = next(SeqIO.parse(str(shared_datadir / PARTIAL_FIXTURE), "genbank"))
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    with open(path, "wb") as handle:
        for feature in record.features:
            if feature.type == "CDS":
                name = feature.qualifiers["locus_tag"][0].encode()
                sequence = pyhmmer.easel.TextSequence(name=name, sequence=feature.qualifiers["translation"][0]).digitize(alphabet)
                hmm, _, _ = builder.build(sequence, background)
                hmm.write(handle)

def domainate_partial_fixture(shared_datadir, output_dir):
    """Annotates the partial fixture with domains named after the CDSs they hit. Returns the path of the annotated genbank file."""
    from domainator import domainate
    hmms = output_dir + "/partial_fixture.hmm"
    write_partial_fixture_hmms(shared_datadir, hmms)
    out = output_dir + "/partial_fixture_domainated.gb"
    domainate.main(["-i", str(shared_datadir / PARTIAL_FIXTURE), "-r", hmms, "-o", out])
    return out

def write_uniprot_fragment_fasta(shared_datadir, path):
    """Writes swissprot_CuSOD_subset.fasta with ' (Fragment)' added to the first header. Returns the id of the fragment."""
    records = list(SeqIO.parse(str(shared_datadir / "swissprot_CuSOD_subset.fasta"), "fasta"))
    with open(path, "w") as handle:
        for i, record in enumerate(records):
            description = record.description.replace(" OS=", " (Fragment) OS=", 1) if i == 0 else record.description
            handle.write(f">{description}\n{str(record.seq)}\n")
    return records[0].id
