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

CI (`.github/workflows/tests.yml`) runs `pytest -m "not network"` on every push and pull request. The `pip` jobs have no external programs, so binary-dependent tests skip; the `conda` job installs them from bioconda and runs everything.

examine test coverage
```
coverage run -m pytest test
coverage report -m
coverage html
```
then open htmlcov/index.html in a browser to see coverage

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
