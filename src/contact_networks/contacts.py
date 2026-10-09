# - Import exact geometry and conformer record utilities.
# - Pin CONTACT radii with verified MMDB provenance.
# - Apply carbonyl and chemically equivalent hydrogen overrides.
# - Reject unknown elements instead of guessing radii.
# - Identify protein backbone atoms without misclassifying cofactors.
# - Build mobile conformers including shared blank atoms.
# - Detect main-chain hydrogen bonds in both directions.
# - Exclude cysteine pairs according to published methods.
# - Suppress bonded backbone pairs within original chains.
# - Calculate symmetric minimum eligible atom-pair ratios.
# - Represent chemically excluded pairs without fabricated scores.
# - Cache relief measurements independently of clash adjacency.
# - Find exact clash candidates using spatial cells.
# - Select percentile ranks with defined endpoint behavior.
# - Preserve threshold ties and deterministic adjacency.
