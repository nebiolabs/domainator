[index](../README.md)
# Python API reference

What 1.x keeps stable in Python is listed in [stability.md](../stability.md): the names the
package root exports, each command's function, and the exceptions.

- [The package root](package.md): `parse_seqfiles`, `write_genbank`, `DataMatrix`, the
  exceptions and the feature-type constants.
- [Command functions](commands.md): `from domainator.<command> import <command>` for every
  command.

The pages render with [mkdocs-material](https://squidfunk.github.io/mkdocs-material/) and
[mkdocstrings](https://mkdocstrings.github.io/), which turn the `:::` lines into the
functions' signatures and docstrings:

```bash
pip install -e ".[docs]"
mkdocs serve      # or: mkdocs build --strict, which writes site/
```

The site is published at https://nebiolabs.github.io/domainator/ from `main` (see
`.github/workflows/docs.yml`).
