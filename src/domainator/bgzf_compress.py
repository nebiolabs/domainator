"""BGZF (block-gzip) compression for domain_search databases.

Produces standard BGZF output (gunzip-readable, and seekable/partitionable by the
native ``genbank_offsets_bgzf`` scanner and ``Bio.bgzf.BgzfReader``). Uses the
native parallel Rust compressor (``_gbfast.bgzf_compress_file``) when available,
falling back to the vendored single-threaded ``BgzfWriter`` otherwise.
"""
import shutil

try:
    from domainator import _gbfast  # optional native (Rust) acceleration
except ImportError:
    _gbfast = None

from domainator.Bio import bgzf


def native_available():
    """True if the native parallel BGZF compressor is built."""
    return _gbfast is not None and hasattr(_gbfast, "bgzf_compress_file")


def compress_to_bgzf(in_path, out_path, level=6, threads=0):
    """Compress ``in_path`` to ``out_path`` as standard BGZF.

    ``threads`` (0 = autodetect) only applies to the native path. The Python
    fallback streams the input through ``BgzfWriter`` so memory stays bounded.
    """
    in_path = str(in_path)
    out_path = str(out_path)
    if native_available():
        _gbfast.bgzf_compress_file(in_path, out_path, int(level), int(threads))
        return out_path
    with open(in_path, "rb") as fh, bgzf.BgzfWriter(out_path, compresslevel=level) as w:
        shutil.copyfileobj(fh, w)
    return out_path
