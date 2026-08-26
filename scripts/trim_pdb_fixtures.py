"""Trim a PDB to first-model backbone+CB ATOM records, for small test fixtures."""
import sys

KEEP_ATOMS = {"N", "CA", "C", "O", "CB"}

def trim(in_path, out_path, atoms=KEEP_ATOMS, chains=None):
    kept = []
    in_model = True
    for line in open(in_path):
        rec = line[:6]
        if rec == "MODEL ":
            in_model = int(line[10:14]) == 1
            continue
        if rec == "ENDMDL":
            in_model = False
            continue
        if not in_model or rec != "ATOM  ":
            continue
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue
        name = line[12:16].strip()
        chain = line[21]
        if name not in atoms:
            continue
        if chains and chain not in chains:
            continue
        kept.append(line)
    with open(out_path, "w") as fh:
        fh.writelines(kept)
        fh.write("END\n")
    return len(kept)

if __name__ == "__main__":
    atoms = KEEP_ATOMS if len(sys.argv) < 4 else set(sys.argv[3].split(","))
    print(sys.argv[2], trim(sys.argv[1], sys.argv[2], atoms), "atoms")
