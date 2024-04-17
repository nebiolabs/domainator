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
