# CONTACT

<!-- - Document staged installation and analysis usage. -->
## Implemented stage

`contact_networks` is an importable modern Python module for detailed-pathway
parsing, filtering, weighted residue networks, JSON/CSV exports, and optional
headless PNGs. It also reads prepared PDB/mmCIF models and optionally expands
explicit biological assemblies, preserving original chains and distinct copies.

**Structure-to-pathway generation is not implemented yet.** The archived C++
source conflicts with the paper's chemistry exclusions and has a branch-local
relief-score defect. A scientific-policy decision is needed before that port.
The existing Python 2 script remains untouched as a historical reference.

### Install and test

```sh
uv sync --extra dev --extra plot
uv run pytest
uv build
```

For compute-only installations omit `--extra plot`. Importing the module does
not execute analysis or import Matplotlib. Prepared input must already contain
hydrogens unless `Settings(allow_hydrogen_free=True)` is explicitly selected;
hydrogen presence does not certify correct placement. All coordinates must be
finite. Models are selected by deposited number, not list position.

### Python

```python
from contact_networks import analyze_pathways, read_pathways, write_results
result = analyze_pathways(read_pathways("pathways.txt"))
write_results(result, "results/run-001", plots=False)
```

`read_structure(path, Settings(model=1, assembly="1"))` reads a selected prepared
model/assembly; it does not discover pathways. Omitting `assembly` uses supplied
coordinates. Unknown assemblies or unmatched explicit chain/subchain selectors
are errors, never silent partial assemblies.

### CLI and external schedulers

```sh
uv run contact-networks extract pathways.txt --output results/run-001 --plots
```

Submit one invocation per input to your scheduler/job array, using distinct
output directories. Existing destinations are refused. Scientific extraction
options are `--relief-threshold`, `--min-network-size`, `--min-pathway-length`,
`--model`, and repeatable `--residue` JSON objects, for example:

```sh
uv run contact-networks extract pathways.txt --output results/selected \
  --residue '{"model":1,"chain":"A","number":99}'
```

Exit codes: **0** complete, **1** input/export error, **2** command-line usage
error, **3** truncated. When importing known incomplete upstream output, use
`--input-status truncated --stopping-reason "upstream enumeration limit"`.
This records existing truncation; it does not impose a new enumeration limit.

Each result directory contains `raw_pathways.csv`, `pathways.csv`,
`residues.csv`, `edges.csv`, `networks.csv`, and a final `results.json` manifest.
Optional images are named `network_0001.png`, etc.; colors do not determine
filenames. JSON carries full residue identities and aligned pathway records;
CSV empty measurements mean missing, not zero. The manifest is published
atomically after successful exports. A failed run may leave diagnostic files
but never publishes a manifest for that failure.

`generation_settings_known` is false for imported detailed pathways: the
reported default generation parameters in `Settings` do not establish which
parameters produced an arbitrary input file. No upstream parameters are
inferred from filenames.

<!-- - Explain historical artifacts and scientific-policy blockers. -->
## DHFR reference example

The original CONTACT supplement includes a prepared WT DHFR model and pathway
artifacts. Fetch its checksum-pinned public archive explicitly:

```sh
uv run examples/fetch_dhfr.py --output reference/dhfr
CONTACT_DHFR_REFERENCE=reference/dhfr uv run pytest tests/test_dhfr.py
uv run contact-networks extract reference/dhfr/pathways.txt --output results/dhfr --plots
```

The fetcher extracts selected source/reference files only, not bundled
executables. It records source DOI, archive hash, and per-file hashes. The
reference artifacts are downloaded on demand rather than redistributed here.
Keep their original notices. The original code's MIT notice is retained in
`UPSTREAM_LICENSE.txt`; paper DOI: `10.1038/nmeth.2592`.

From the supplied detailed output at relief threshold 0.90:

| Quantity | Reference |
| --- | ---: |
| Detailed raw pathways | 2,415 |
| Ordered maximal residue sequences, before relief filtering | 157 |
| Relief-passing pathways | 340 |
| Retained networks with at least three residues | 9 |
| Retained residues / weighted graph edges | 56 / 52 |
| Networks with at least four residues | 4 |
| Residues in those four networks | 41 |

These reproduce the **saved artifact**, not a rerun of the generator. Historical
relief scores may reflect its DFS scoring defect. The paper's stated minimum
is four residues; the requested new default is three. Tests explicitly distinguish
the two. Do not treat 157 as the relief-filtered pathway count.

### Documented differences

- Preserve model, chain, insertion code, and assembly copy; residue numbers alone
  no longer merge unrelated residues.
- Flush the last pathway at EOF, reject malformed detailed headers/records, and
  ignore empty pathway separators rather than creating phantom paths.
- Only the first relief value is a missing predecessor marker. Later zero
  values remain real zero values and fail positive relief thresholds.
- Remove undersized connected components from returned graphs, not just plots.
- Count consecutive undirected edge occurrences, including repeated transitions.
- Eliminate only ordered contiguous subpaths; reversal or equal residue sets
  do not imply equivalent pathways.
- Keep raw residue occurrences separate from relief/length/residue filtering.
- Exclude recognized waters, including heavy water; keep other cofactors.

<!-- - Record reproducibility checks and measured optimizations. -->
## Optimizations and audit trail

Graph aggregation uses direct edge lookup rather than scanning every existing
edge per pathway transition. Network pathway counts use component membership
instead of rescanning every pathway for every component. Ordered-subpath
comparisons use an index of potential starting positions. Export occurrence
counts are aggregated once; plots are opt-in.

```sh
uv run benchmarks/benchmark_networks.py --paths 2000 --residues 100
```

The benchmark asserts exact weighted-graph equality before reporting runtime
and Python allocation peaks. Allocation peaks are not total process RSS;
this does not benchmark native chemistry or pathway enumeration.

One local run (Python 3.13.11, NetworkX 3.7; 2,000 four-residue pathways and
100 residues) measured 1.419 seconds for the legacy scan and 0.0455 seconds
for direct aggregation: approximately 31x faster. Python allocation peaks
were 71,222 and 104,561 bytes respectively; this is a speed improvement,
not a demonstrated memory reduction. Re-run on your deployment workload.

Implementation is staged as comment-outline commits followed by executable
commits. Tests enforce bullet comments of at most eight words and code sections
of at most ten nonempty physical lines across authored Python files. Fresh
GPT-6 Luna adversarial reviewers inspect each commit.

The first core/structure implementation reviews found mutable nested records,
permissive header recognition, retained heavy water, and incomplete assembly
selection. These were fixed in `74abba3`, with regression tests and two fresh
reviewers reporting no further findings in that fix. The outline-stage review
also requested an explicit finite-coordinate acceptance criterion, documented
above and covered by tests. No review result is a claim of exhaustive correctness.

<!-- - Verify installed wheel imports and upstream attribution. -->
<!-- - Clarify source truncation versus selected model completeness. -->
