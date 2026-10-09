# - Expose analysis types and network operations.
from .models import Atom, Pathway, Residue, ResidueId, Settings
from .networks import (
    build_graph, filter_pathways, maximal_pathways, network_pathway_counts,
    read_pathways, residue_counts,
)

# - Declare the package version without side effects.
__version__ = "0.1.0"

# - Expose staged analysis and serialization helpers.
# - Expose prepared structure reading without generation.
