# - Import structure reading and immutable residue records.
# - Read PDB and mmCIF with identical semantics.
# - Select explicit model without mixing coordinate sets.
# - Expand only explicitly selected biological assemblies.
# - Preserve original chains and generated copy identity.
# - Normalize alternate locations and insertion codes.
# - Reject duplicate residues and invalid coordinates.
# - Keep shared atoms and relevant cofactors.
# - Reject hydrogen-free inputs unless explicitly allowed.
# - Return deterministic immutable prepared residues.
