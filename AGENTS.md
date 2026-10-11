<!-- - Describe project scope and module boundaries. -->
# CONTACT agent guide

## Scope and entry points
- `contact_networks` requires Python >=3.11; use `uv` and the `uv_build` backend.
- Analyze **prepared** PDB/mmCIF inputs. qFit, refinement, and hydrogen placement
  remain external; hydrogen-free inputs require explicit opt-in.
- Public API: `analyze(path)` runs generation and network analysis;
  `generate_pathways(residues)` accepts residue records;
  `analyze_pathways(read_pathways(path))` imports saved scores without regeneration.
  CLI `contact-networks analyze` and `extract` preserve that distinction.
- In `src/contact_networks/`: `models.py` owns records/settings, `structures.py`
  input/assembly identity, `contacts.py` chemistry, `engine.py` DFS generation,
  `networks.py` parsing/filtering/graphs, `analysis.py` orchestration,
  `output.py` exports, and `cli.py` argument/status handling.
- Keep `scripts/network_extractor_v1.py` as historical Python 2 reference.
  Do not execute bundled historical binaries.

<!-- - Provide verified setup and validation commands. -->
## Checks from the repository root

```sh
uv sync --extra dev --extra plot
uv run pytest
uv build
```

`dev` supplies pytest; `plot` is optional for compute-only use. Tests mirror module
names. Run affected tests while editing, then the suite; build for packaging changes.
Optional DHFR integration (fetch once; the fetcher refuses an existing destination):

```sh
uv run examples/fetch_dhfr.py --output reference/dhfr
CONTACT_DHFR_REFERENCE=reference/dhfr uv run pytest
```

For performance changes, verify equality/count assertions before comparing timings:

```sh
uv run benchmarks/benchmark_networks.py --paths 2000 --residues 100
uv run benchmarks/benchmark_generation.py --copies 100
```

These are local CPU workloads, not cluster claims; generation RSS uses Unix `resource`.

<!-- - State scientific invariants and reference distinctions. -->
## Scientific changes require explicit evidence
- Read `README.md` → **Scientific corrections** before changing chemistry or DFS;
  use `tests/test_contacts.py` and `tests/test_engine.py` as regression entry points.
- Preserve `(model, chain, number, insertion, copy)` throughout; never key by bare
  residue number. Split original-chain parts under one assembly operator share
  copy identity; distinct operator copies remain distinct (`tests/test_structures.py`).
- Contact ratios are distance / sum of radii, using selected conformers and shared
  blank-altloc atoms. Preserve CYS-pair and selected-conformer backbone H-bond masks.
- Scores belong to each DFS branch's conformer choices. Every switch clears earlier
  selected states; length exhaustion must not manufacture a terminal pathway.
- Subpath elimination means **ordered contiguous** sequences, not sets or reversals.
  Real zero is measured data, never missing/excluded; preserve explicit exclusion flags.
- `max_pathway_length` counts transitions, not residues. `max_candidates` caps
  globally expanded DFS states across groups, not emitted pathways.
- Keep truncation and stopping reasons explicit (CLI exit 3). Imported provenance
  has unknown generation settings and source-input status, even for empty model subsets.
- In `tests/test_dhfr.py`, historical artifact import and corrected generation are
  separate tests. Keep archived 157 maximal / 340 relief-passing expectations;
  never replace them with new generation counts to bless changed science.
  Corrected results are not experimental validation of allosteric networks.

<!-- - Preserve comment, output, and commit conventions. -->
## Editing and delivery contracts
- Authored Python in `src/`, `tests/`, `examples/`, and `benchmarks/` starts with
  `# - `. Bullet comments have <=8 words; sections between bullets have <=10
  nonempty non-comment physical lines. `tests/test_comment_contract.py` is the
  executable rule; use meaningful section comments, not additional filler.
- Use fresh output directories; refuse existing destinations. Publish `results.json`
  atomically **last**, only after successful exports (`tests/test_output.py`).
  Keep Matplotlib optional, lazy-imported, and headless.
- Do not commit generated results, reference downloads, environments, caches, build
  outputs, or `.DS_Store`; preserve upstream notices in source and distributions.
- Follow explicitly requested comments-first/audited commit workflows; the migration's
  reviewer model choices are not permanent requirements. Use Conventional Commits
  (`type(scope): imperative summary`, <=72 characters); never push unless asked.

<!-- - Cite current agent-context research without overloading guidance. -->
## Context design note — sources checked 2026-10-10

Keep this guide focused on commands and costly repository-specific mistakes; link
to detailed science instead of copying the README. This is a maintenance choice,
not a proven performance gain for CONTACT. Evidence is mixed:

- [OpenAI instruction discovery](https://developers.openai.com/codex/guides/agents-md):
  root-to-directory scope and overrides; official living documentation.
- [OpenAI customization](https://developers.openai.com/codex/concepts/customization):
  recommends small, durable guidance and separate reusable workflows.
- [Gloaguen et al., v3, 2026-09-29](https://arxiv.org/abs/2602.11988v3):
  no significant success gain; generated files had small, nonsignificant declines
  and increased costs. Length alone was not a demonstrated cause of poorer results.
- [Lulla et al., v2, 2026-03-30](https://arxiv.org/abs/2601.20404v2):
  124 PRs across 10 repositories showed lower median runtime/output tokens with
  instructions. Different tasks and metrics limit generalization of either study.
