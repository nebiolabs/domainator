[index](README.md)
# Developing Domainator

## testing

run all tests
```
pytest test
```

Two markers are registered in `pyproject.toml`:

- `requires_binary("diamond", ...)`: the test is skipped unless every named external program is on `PATH` (implemented in `test/conftest.py`). Use it for any test that runs an external program, directly or through a tool's default algorithm.
- `network`: the test downloads from the internet. Skip these with `pytest -m "not network"`.

Domainator's own `FutureWarning`s are turned into errors in the test suite, so tests can't rely on deprecated usage.

CI (`.github/workflows/tests.yml`) runs `pytest -m "not network"` on pull requests, and when started by hand from the Actions tab (Tests > Run workflow). The `pip` jobs have no external programs, so binary-dependent tests skip; the `conda` job installs them from bioconda and runs everything. `.github/workflows/docs.yml` builds the documentation site with `mkdocs build --strict` (a broken link, a missing `#heading` or an unparsable docstring fails it), and on pushes to `main` publishes it to GitHub Pages at https://nebiolabs.github.io/domainator/.

examine test coverage
```
coverage run -m pytest test
coverage report -m
coverage html
```
then open htmlcov/index.html in a browser to see coverage

## parallel workers

Every pool of workers comes from `domainator._mp.make_pool`. On free-threaded Python with the GIL off, its workers are threads, which share memory; otherwise they are processes started with spawn. A task gets large read-only data (parsed references, matrices) through `_mp.Shared`, never through module globals or attributes of the task's function: a process pool sends the function and its arguments with every task, but a `Shared` reaches each worker once, and threads read the parent's copy. A task must not change anything it shares with other tasks.

`DOMAINATOR_MP_START_METHOD` overrides the choice: `threads`, or a multiprocessing start method (`spawn`, `forkserver`, `fork`). Running the tests with `DOMAINATOR_MP_START_METHOD=threads` checks the thread path on an ordinary Python. On free-threaded Python, pyhmmer, pyinfernal and pyrodigal still turn the GIL back on when imported (they haven't declared that they can run without it), so pools are processes there too unless you start Python with `PYTHON_GIL=0`.

## the SSN viewer

The viewer page is built from files in `src/domainator/ssn_viewer/`: `viewer.html`, `viewer.css` and the scripts in `js/`. `js/viewer.js` is the page's one script; each `// @include <file>` line pulls in a file from `js/`, which then shares its scope. Values from Python (versions, palettes, an embedded bundle) reach the page through its `dsnv-config` JSON block (`DSNV_CONFIG` in the script), set in `ssn_viewer_html.py`. `test/test_ssn_viewer_js.py` checks the assembled script with `node --check` and tests the layout worker in Node; `test/test_ssn_viewer_browser.py` drives the page in Chromium (see "Browser/SSN-viewer tests" in CLAUDE.md).

## the HTML reports

The `--html_out` pages are built by `report_html.py` from files in `src/domainator/web/`: `reports/<name>.html` and `reports/<name>.js` for each page, and `common/` for the stylesheet and the scripts the pages share (`payload.js`, `download.js`, `table.js`, and `charts.js`, `line_chart.js`, `hierarchy.js` and `pack.js` for SVG charts drawn without a charting library). As in the viewer, `// @include <path>` (relative to `web/`) pulls a file into the page's script. There is no third-party JavaScript, so the pages work offline. A page's data goes in its `report-data` JSON block (an envelope with `payload_version`, which changes whenever the data a page reads changes shape); a table page streams its rows into a `report-rows` NDJSON block. Anything a page draws that another implementation must reproduce is computed in Python and put in the data with an integer-only rule, not left to the page (summary_report's `length_histogram` bins and `taxonomy_tree` order); layout that only the page uses, like where matrix_report's slider puts each stop or how plot_contigs stacks features, stays in the page. A page that needs a Python constant keeps a copy, checked by a test (matrix_report's `SPLIT_AXIS_LABELS` against `merge_impact_axis_labels`). `report_html.report_shell(name)` renders a page with tokens in place of its title, data and rows, for programs that write these pages without Python (domainatoredux vendors them). `test/test_report_html.py` checks the data and runs `node --check` on each page's script; `test/test_report_browser.py` drives the pages in Chromium with the network blocked.

## benchmarks

`scripts/bench/bench.py` times parsing, writing, `domainate`, `domain_search` and `extract_domains` on a test genome. Save timings before a change and compare after:
```
python scripts/bench/bench.py --json_out before.json
python scripts/bench/bench.py --baseline before.json
```
Each step reports its fastest of `--repeat` runs; `--only domainate` runs a subset.


## Projects that domainator depends on

  - hmmer3
  - usearch
  - python3
  - biopython
  - pytest
  - pytest-datadir
  - pandas
  - seaborn
  - cd-hit
  - scipy
  - pyhmmer
  - umap-learn
  - diamond
  - foldseek (optional, for the structure_* tools; installed with the conda environment)
  - reseek (optional, alternative structure backend; must be built from source, see [Structure search and annotation](structure_search.md))
  - coverage
