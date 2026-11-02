#!/usr/bin/env bash
# Build the sdist and wheel, install each into a clean venv, and smoke-test the installed package
# (see "release checks" in docs/developing_domainator.md).
#
#     conda activate domainator             # python, and the Rust toolchain for the native extension
#     scripts/release_check.sh [work_dir]   # work_dir must be new or empty; default: a new temporary directory
#
# The wheel is built from the sdist, so a file the build needs but MANIFEST.in leaves out fails it.
# Every command runs outside the repository, so it tests the installed package, not src/.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=${1:-$(mktemp -d)}
DATA=$REPO/test/data

mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)
if [ -n "$(ls -A "$WORK")" ]; then
    echo "error: $WORK is not empty" >&2
    exit 1
fi
cd "$WORK"
echo "working in $WORK"

echo "=== build the sdist, then the wheel from the sdist ==="
python -m venv buildenv
buildenv/bin/pip install -q build twine
buildenv/bin/python -m build --outdir dist "$REPO"
ls -l dist
buildenv/bin/twine check --strict dist/*

smoke() {  # smoke <venv> <package file>
    local venv=$WORK/$1 package=$2
    echo "=== install $(basename "$package") into a clean venv ==="
    python -m venv "$venv"
    "$venv/bin/pip" install -q "$package"
    local out=$venv-out && mkdir -p "$out" && cd "$out"
    "$venv/bin/python" -c "import domainator, domainator._gbfast; print('imported from', domainator.__file__)"
    "$venv/bin/domainator" --version
    "$venv/bin/domainator" doctor
    "$venv/bin/domainator" completion bash | bash -n
    "$venv/bin/domainator" describe > manifest.json
    "$venv/bin/domainator" viewer --html_out viewer.html   # needs the viewer page's package data
    # the smoke workflow: domainate -> select_by_cds -> enum_report
    "$venv/bin/domainator" domainate -i "$DATA/pDONR201_multi_genemark.gb" -r "$DATA/pdonr_hmms.hmm" -o annotated.gb
    "$venv/bin/domainator" select_by_cds -i annotated.gb --domains CcdB --cds_up 1 --cds_down 1 -o selected.gb
    "$venv/bin/domainator" enum_report -i selected.gb --by contig --columns length domains -o report.tsv --html_out report.html
    cat report.tsv
    for f in annotated.gb selected.gb report.tsv report.html viewer.html manifest.json; do
        [ -s "$f" ] || { echo "FAIL: $f is empty"; exit 1; }
    done
    grep -q "CcdB" report.tsv || { echo "FAIL: CcdB not in report.tsv"; exit 1; }
    cd "$WORK"
}

smoke wheelenv "$(ls dist/*.whl)"
smoke sdistenv "$(ls dist/*.tar.gz)"
echo "=== all release checks passed ==="
