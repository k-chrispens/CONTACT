# - Import graph, counting, and pathway utilities.
# - Finalize detailed records at boundaries and EOF.
# - Reject malformed and unaligned pathway fields.
# - Preserve chain, insertion, and conformer identity.
# - Distinguish initial relief markers from zero relief.
# - Stream records without empty phantom pathways.
# - Apply relief, length, and residue filters.
# - Count residue occurrences before scientific filtering.
# - Aggregate repeated undirected transitions directly.
# - Remove connected components below minimum size.
# - Count pathways touching each retained network.
# - Remove duplicate ordered residue sequences.
# - Remove contiguous subpaths without reversing sequences.
