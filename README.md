# CONTACT

<!-- - Document staged installation and analysis usage. -->
## Usage

`contact_networks` is an importable modern Python module for prepared PDB/mmCIF
structures → corrected pathways → weighted residue networks → JSON/CSV and
optional headless PNGs. Optional biological assembly expansion preserves
original chains and distinct operator copies, including split source-chain parts.
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
from contact_networks import analyze, write_results
result = analyze("prepared_structure.cif")
write_results(result, "results/run-001", plots=False)
```

`read_structure(path, Settings(model=1, assembly="1"))` reads a selected prepared
model/assembly. `generate_pathways(residues, settings)` exposes the low-level
corrected generator; `analyze_pathways(read_pathways("pathways.txt"))` imports
existing detailed scores without regenerating them. Omitting `assembly` uses supplied
coordinates. Unknown assemblies or unmatched explicit chain/subchain selectors
are errors, never silent partial assemblies.

### CLI and external schedulers

```sh
uv run contact-networks analyze prepared_structure.cif --output results/run-001 --plots
uv run contact-networks extract pathways.txt --output results/imported
```

Submit one invocation per input to your scheduler/job array, using distinct
output directories. Existing destinations are refused. Scientific extraction
options are `--relief-threshold`, `--min-network-size`, `--min-pathway-length`,
`--model`, and repeatable `--residue` JSON objects, for example:

```sh
uv run contact-networks extract pathways.txt --output results/selected \
  --residue '{"model":1,"chain":"A","number":99}'
```

Structure analysis additionally supports `--clash-percentile`,
`--max-pathway-length`, `--max-candidates`, `--side-chain-only`, `--cross-chain`,
`--assembly ID`, and `--allow-hydrogen-free`. Searches are chain-local by default;
cross-chain discovery and biological assembly expansion are separately opt-in.
Minimum pathway length requires strictly more residues than the supplied value;
maximum pathway length counts transitions.

Exit codes: **0** complete, **1** input/export error, **2** command-line usage
error, **3** truncated. When importing known incomplete upstream output, use
`--input-status truncated --stopping-reason "upstream enumeration limit"`.
This records existing truncation; it does not impose a new enumeration limit.
The manifest's `status_scope` is `source_input`. An empty selected-model subset
does not establish that model's enumeration is complete when the source is
truncated; exit 3 and source truncation are retained even with zero matching
records. No selected-model completeness is inferred from absence of records.
For structure analysis, status scope is `generation` and the state budget is
enforced by the generator itself. JSON records the expanded `candidate_states`,
completed groups, group-specific cutoffs, known generation parameters, and
`science_revision: "contact-corrected-v1"`. A complete result is exhaustive within
its configured chain/model/assembly and transition-length scope.

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
uv run contact-networks analyze reference/dhfr/dhfr.pdb --output results/dhfr-corrected --plots
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

```sh
uv run benchmarks/benchmark_generation.py --copies 100
```

An independently specified workload of 100 three-residue cascades verifies
400 completed pathways and 600 expanded states before reporting timings.
One local run completed in 0.155 seconds with 379,934 Python-tracked peak bytes
and a 54.4 MB whole-process peak RSS. This is a local CPU measurement, not
cluster throughput or a native historical-generator comparison.

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
Wheel tests verify the complete upstream notice and `License-File` metadata
are present in built distributions, not only in this source checkout.
<!-- - Clarify source truncation versus selected model completeness. -->

<!-- - Document corrected chemistry and state-transition criteria. -->
## Scientific corrections

New pathway generation uses corrected behavior only; there is no second engine
that deliberately reproduces archived bugs. These changes restore the stated
CONTACT method and make its numerical decisions explicit. They do **not**
constitute experimental validation of the inferred allosteric networks.

| Area | Archived implementation | Corrected generation |
| --- | --- | --- |
| Backbone hydrogen bonds | The MMDB exclusion was commented out. | Exclude protein main-chain/main-chain contacts when either donor direction satisfies the MMDB geometric criterion, evaluated for the selected conformers. Side-chain-involving contacts remain eligible. |
| Cysteine pairs | The apparent CYS SG–SG exclusion only printed a message. | Exclude comparisons between CYS residues at the residue-pair level, following the paper's explicit “pairs of cysteine residues” wording. This is broader than excluding only SG–SG atom pairs. |
| Backbone classification | Names alone could misclassify ligand atoms; neighboring sequence numbers could suppress inter-chain contacts. | Apply backbone rules to amino-acid residues only. Sequence-neighbor exclusions require the same model, original chain, and assembly copy. Include terminal backbone names and glycine alpha hydrogens. |
| Hydrogen radii | The lysine tests repeated HZ1 and omitted HZ2. Isotope/name variants could receive inconsistent radii. | Treat HZ1/HZ2/HZ3 equally. Map deuterium names to equivalent hydrogen names/radii. Alpha hydrogens retain 1.17 Å; specified polar/aromatic hydrogens use 1.00 Å. |
| Mobile-state population | Static residues and single named states could enter the contact population. | Evaluate residue pairs with at least two actual nonblank alternate states in the selected atom set. Shared blank atoms participate in every conformer; backbone-only alternatives do not create side-chain mobility. |
| Starting and continuing a pathway | A seed's direction was not validated, and propagation could follow an already more severe clash. | Require the before-to-after switch to introduce or worsen each propagated clash: the previous normalized distance must be at least the new one, including equality as in the paper's inequalities. |
| Branch-local relief | A shared relief vector grew without backtracking. | Carry immutable per-branch relief tuples and reconstruct each measurement from that pathway's conformer choices. A sibling branch cannot change another pathway's score. |
| Terminal validation | Terminal alternatives were checked only against the immediate predecessor; reaching the length limit could produce a false terminal path. | Every switch must clear clashes with all earlier selected path states. Emit only when no admissible new clash requires propagation. Discard unresolved branches at the configured transition limit. |
| Missing versus excluded measurements | Uninitialized or sentinel zero values could be interpreted as measurements. | Preserve real zero and distinguish chemically excluded pairs from measured distances. Generated JSON `relief_exclusions` and CSV `relief_excluded` explicitly explain a missing numeric relief after exclusion; no synthetic 0, 1, or infinity is substituted. |
| Percentile endpoints | The 100% rank indexed past the end of the sample. | Use rank `min(floor(p*n), n-1)` without interpolation, retaining ties. Zero percent selects the minimum, 100% the maximum; an empty clash population has cutoff zero and no edges. |
| Resource and length boundaries | A terminal-path counter exceeded its nominal limit; maximum length meant edges rather than residues. | `max_pathway_length` explicitly counts transitions: 100 allows at most 101 residues. `max_candidates` caps expanded DFS states across all groups, stopping before an over-budget expansion, with partial results and a distinct status. |

The main-chain hydrogen-bond criterion is **N–O distance ≤ 3.5 Å and angle
C–O–N ≥ 90°**, with nondegenerate vectors. It is independently implemented
from the criterion in MMDB2 2.0.22's `Residue::isMainchainHBond`; we inspect
both directions and do not arbitrarily choose the first alternate atom.
This classifier is the MMDB convention, not a new physical force field.

The supported element radii are pinned to numerical MMDB2 model values with
the original CONTACT H/O/C overrides and explicit atom corrections. Unknown
elements in participating conformers fail rather than receiving a guessed
radius. Supported elements are H, D, C, N, O, P, S, F, Cl, Br, I, Se, Li, Na,
K, Mg, Ca, Mn, Fe, Co, Ni, Cu, and Zn. Backbone carbonyl C is 1.65 Å; generic
C is 1.75 Å, N 1.55 Å, O 1.40 Å, H/D 1.17 Å, and S/P 1.80 Å.

Exact spatial cells discover potential clash pairs; on-demand cached
measurements still evaluate distant pairs used for relief. Read-only coordinate
arrays and vectorized atom-pair calculations avoid repeated scalar work. Tests
compare spatial discovery with an independently calculated exhaustive oracle,
check hydrogen-bond distance/angle boundaries in both directions, and recompute
emitted branch scores from small explicit coordinate examples.

Scientific references:
- [CONTACT paper, online methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC3760795/),
  DOI `10.1038/nmeth.2592`.
- [MMDB2 2.0.22 source archive](https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/dependencies/mmdb2-2.0.22.tar.gz):
  `mmdb2/mmdb_atom.cpp` and `mmdb2/mmdb_tables.cpp`; archive SHA-256
  `a9646933ce9f34633e1ae4901f2383af0e5318d6490851328f5b6aa06195ab51`.
  This is a verified numerical/geometry reference, not a claim that the
  historical MMDB 1.23.2.2 build was independently reproduced.

<!-- - Separate historical references from corrected generation results. -->
<!-- - Document full API, cluster CLI, and corrected measurements. -->
With default parameters, the supplied prepared DHFR structure completed in one
local run with **19,406 expanded states, 4,800 raw pathways, 260 maximal raw
residue sequences, and 2,057 accepted pathways**. The retained graph had
112 residues, 133 weighted edges, and four networks of sizes 102, 4, 3, and 3.
Chain A's stress cutoff was approximately 0.886573. Low-level generation took
about 0.60 seconds with a 66.6 MB whole-process peak RSS on the local machine.

These corrected measurements deliberately differ from the saved artifact's
2,415/157/340 counts and nine networks. The larger inferred network is an
algorithmic result of the documented corrections, not confirmation of improved
biological accuracy. The optional corrected DHFR test checks completion,
provenance, and structural invariants; historical counts remain separate.

The historical DHFR test checks imported saved scores and the 157 maximal
residue sequences. It remains useful for parser/export compatibility, but it
does **not** validate corrected generation. Corrected cutoff values, pathways,
weights, and network memberships may differ. Do not silently replace historical
expected values with corrected results or label corrected differences as
successful reproduction of the original paper.
