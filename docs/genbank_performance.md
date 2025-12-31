# GenBank Parsing Performance Optimization Guide

This document outlines performance bottlenecks in GenBank parsing and potential optimizations for Domainator.

## Current Bottlenecks

### 1. Line-by-Line I/O and String Operations
The scanner processes files line-by-line with heavy string operations (slicing, strip, split, replace) on every line. For large files this creates significant overhead.

### 2. Regex Compilation
Several regex patterns in `utils.py` and `SeqFeature.py` are compiled at module load, but some patterns are matched repeatedly without caching results. Location parsing is particularly expensive.

### 3. Object Creation Overhead
Every feature creates multiple objects:
- `SeqFeature` instance
- `SimpleLocation` or `CompoundLocation` 
- 2+ position objects (`ExactPosition`, `BeforePosition`, etc.)
- Qualifier dict and lists

For a GenBank file with 10,000 CDS features, this creates 40,000+ objects.

### 4. Consumer Pattern Indirection
The `_FeatureConsumer` callback architecture adds method call overhead for every parsed element.

### 5. No `__slots__` on Core Classes
`SeqFeature`, `SeqRecord`, and position classes use `__dict__` instead of `__slots__`, increasing memory usage and attribute access time.

## Optimization Strategies

### High Impact (Recommended)

#### 1. Add `__slots__` to SeqFeature and Position Classes
```python
# In SeqFeature.py
class SeqFeature:
    __slots__ = ('location', 'type', 'id', 'qualifiers')
    # ...

class ExactPosition(int):
    __slots__ = ()
    # ...
```
**Expected improvement:** 20-30% memory reduction, 5-10% faster attribute access.

#### 2. Lazy Location Parsing
Store the raw location string and only parse to objects when accessed:
```python
class SeqFeature:
    __slots__ = ('_location', '_location_str', 'type', 'id', 'qualifiers')
    
    @property
    def location(self):
        if self._location is None:
            self._location = parse_location(self._location_str)
        return self._location
```
**Expected improvement:** Significant if many features are filtered out before location access.

#### 3. Buffered Reading
Read larger chunks and parse in memory instead of line-by-line:
```python
def parse_records_buffered(handle, buffer_size=1024*1024):
    buffer = handle.read(buffer_size)
    # Split on record boundaries "//\n" and process
```
**Expected improvement:** 15-30% faster I/O for large files.

#### 4. Use `re.Scanner` for Tokenization
Replace line-based parsing with a regex scanner for the feature table:
```python
import re
feature_scanner = re.Scanner([
    (r'^\s{5}(\w+)\s+', lambda s, t: ('FEATURE_KEY', t)),
    (r'^\s{21}/(\w+)=', lambda s, t: ('QUALIFIER_KEY', t)),
    # ...
])
```
**Expected improvement:** 10-20% faster feature parsing.

### Medium Impact

#### 5. Pre-compiled Location Patterns with Caching
Use `functools.lru_cache` for location parsing:
```python
from functools import lru_cache

@lru_cache(maxsize=10000)
def parse_location_cached(location_str: str, seq_length: int, is_circular: bool):
    return _parse_location(location_str, seq_length, is_circular)
```
**Expected improvement:** Significant when same locations repeat (e.g., gene + CDS pairs).

#### 6. Intern Common Strings
Many qualifier keys and feature types repeat. Use `sys.intern()`:
```python
import sys
feature_key = sys.intern(feature_key.strip())
qualifier_name = sys.intern(qualifier_name)
```
**Expected improvement:** Reduced memory, faster string comparisons.

#### 7. Skip Unnecessary Features
Add early filtering in the scanner:
```python
def parse_features(self, skip=False, feature_filter=None):
    # Skip parsing features that don't match filter
    if feature_filter and feature_key not in feature_filter:
        # Skip qualifier parsing entirely
        continue
```
**Expected improvement:** Linear speedup based on filter selectivity.

### Lower Impact / More Complex

#### 8. Cython/PyPy Compilation
The scanner code is pure Python and would benefit from compilation:
- Cython: Compile `Scanner.py` and `utils.py`
- PyPy: Run Domainator under PyPy for JIT optimization

**Expected improvement:** 2-5x faster with Cython, 3-10x with PyPy.

#### 9. Parallel Record Parsing
Parse independent records in parallel:
```python
from concurrent.futures import ProcessPoolExecutor

def parse_parallel(handle, num_workers=4):
    # Read and split into record chunks
    record_texts = split_on_record_boundary(handle.read())
    with ProcessPoolExecutor(num_workers) as pool:
        yield from pool.map(parse_single_record, record_texts)
```
**Expected improvement:** Near-linear scaling for multi-record files.

#### 10. Alternative Parser (msgpack/protobuf)
For Domainator's internal pipeline, convert GenBank to a binary format:
```python
# First pass: GenBank -> msgpack
# Subsequent passes: msgpack only (10x+ faster)
```
**Expected improvement:** 10-50x faster for repeated reads.

## Single-Pass Parsing for Parallel Workflows

For parallel processing workflows like `domain_search`, Domainator now supports a single-pass parsing mode that eliminates double file reading **when Z is pre-specified**.

### The Problem with Offset-Based Partitioning

The original parallel approach:
1. **First pass:** Scan file with binary search to find byte offsets of record boundaries
2. **Second pass:** Workers seek to offsets and parse their assigned sections

This approach:
- Reads each file twice
- Requires seekable files (incompatible with compressed input)
- Has overhead from binary scanning

### Single-Pass Solution

The new `--single_pass` mode (when `-Z` is specified):
1. Read file once in main process
2. Split raw text on record boundaries (LOCUS/CDS for GenBank, `>` for FASTA)
3. Batch raw text chunks and send to workers
4. Workers parse raw text in parallel

Benefits:
- **No double reading:** File is read exactly once
- **Cheap IPC:** Raw text strings are fast to serialize vs. SeqRecord objects
- **Parallel parsing:** Workers still parse independently

### Limitation: Z=0 Mode

When `-Z 0` is used (auto-count mode), single-pass provides **no benefit**:
- Workers need Z to calculate accurate e-values
- Z isn't known until the entire file has been read
- Buffering all raw text while counting would consume huge memory

In this case, `--single_pass` is ignored and the original offset-based approach is used.

### Usage

```bash
# Enable single-pass mode (requires -Z to be specified)
domain_search -i input.gb -r queries.hmm -o out.gb --single_pass -Z 1000

# With Z=0, single_pass is ignored (falls back to offset-based)
domain_search -i input.gb -r queries.hmm -o out.gb --single_pass -Z 0
```

### API

```python
from domainator.domain_search import domain_search_raw
from domainator.utils import i_partition_raw_records, parse_raw_text

# Low-level: partition files into raw text batches
for raw_text, filetype in i_partition_raw_records(['file.gb'], cdss_per_partition=10000):
    # raw_text: str containing one or more GenBank records
    # filetype: 'genbank' or 'fasta'
    records = list(parse_raw_text(raw_text, filetype))

# High-level: run domain search with single-pass mode (requires Z to be specified)
for hit in domain_search_raw(
    input_paths=['file.gb'],
    references=['query.hmm'],
    cdss_per_partition=10000,
    # ... other args
):
    print(hit.id)
```

## Benchmarking

Before optimizing, establish baselines:

```python
import cProfile
import pstats
from domainator.Bio import SeqIO

# Profile parsing
cProfile.run(
    'list(SeqIO.parse("large_file.gb", "genbank"))',
    'genbank_parse.prof'
)

# Analyze hotspots
stats = pstats.Stats('genbank_parse.prof')
stats.sort_stats('cumulative').print_stats(20)
```

Key metrics to track:
- Records/second
- Memory per record
- Time in location parsing vs I/O vs object creation

## Quick Wins for Domainator

1. **Add `__slots__` to SeqFeature** - Low risk, easy to implement
2. **Cache location parsing** - Add `@lru_cache` to `_loc()` function
3. **Feature filtering in scanner** - Skip non-CDS features when only CDS needed
4. **Intern feature types and qualifier keys** - One-line changes

## Implementation Priority

| Priority | Change | Effort | Impact |
|----------|--------|--------|--------|
| 1 | `__slots__` on SeqFeature/positions | Low | Medium |
| 2 | Feature type filtering | Low | High (use-case dependent) |
| 3 | Location parsing cache | Low | Medium |
| 4 | String interning | Low | Low-Medium |
| 5 | Single-pass partitioning | Medium | High |
| 6 | Buffered reading | Medium | Medium |
| 7 | Lazy location parsing | Medium | High |
| 8 | Cython compilation | High | High |
