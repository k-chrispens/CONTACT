# - Import corrected contact maps and pathway records.
# - Record bounded generation and completed group provenance.
# - Group selected residues by chain and assembly copy.
# - Build alternate-state lookup without mixing residue identities.
# - Seed clashes introduced by an alternate-state change.
# - Carry branch-local scores through iterative depth-first search.
# - Require every switch to relieve previous path clashes.
# - Propagate only newly introduced or worsened clashes.
# - Emit only terminated chemically valid pathways.
# - Bound transitions without mislabeling unfinished terminal branches.
# - Enforce candidate-state budgets across all groups.
# - Return partial results with explicit stopping reasons.
