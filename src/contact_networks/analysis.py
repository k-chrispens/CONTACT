# - Import shared scientific and graph operations.
from dataclasses import dataclass, field
import networkx as nx
from .models import Pathway, ResidueId, Settings
from .networks import build_graph, filter_pathways, maximal_pathways

# - Record complete and explicitly truncated analysis results.
@dataclass(frozen=True)
class Analysis:
    source: str
    settings: Settings
    raw_pathways: tuple[Pathway, ...]
    pathways: tuple[Pathway, ...]
    graph: nx.Graph
    maximal_paths: tuple[tuple[ResidueId, ...], ...]
    # - Record completion status and input processing provenance.
    status: str = "complete"
    stopping_reason: str | None = None
    thresholds: dict[str, float] = field(default_factory=dict)
    processed_groups: tuple[str, ...] = ()
    input_kind: str = "detailed_pathways"
    # - Prevent incomplete results masquerading as complete analyses.
    def __post_init__(self):
        if self.status not in ("complete", "truncated"):
            raise ValueError("status must be complete or truncated")
        if (self.status == "truncated") != bool(self.stopping_reason):
            raise ValueError("only truncated results require a stopping reason")
        if self.input_kind != "detailed_pathways":
            raise ValueError("this stage analyzes detailed pathways only")

# - Analyze prepared structures through pathway generation.
# - Carry generation thresholds, budgets, and explicit science provenance.
# - Analyze supplied detailed pathways without regeneration.
def analyze_pathways(paths, settings=None, *, source="", status="complete", stopping_reason=None,
                     thresholds=None, processed_groups=()):
    settings = Settings() if settings is None else settings
    raw = tuple(paths)
    if any(not isinstance(path, Pathway) for path in raw):
        raise TypeError("analysis requires Pathway records")

# - Preserve raw counts before scientific filtering.
    accepted = filter_pathways(raw, settings)

# - Build retained networks and maximal residue sequences.
    return Analysis(str(source), settings, raw, accepted, build_graph(accepted, settings.min_network_size),
                    maximal_pathways(raw), status, stopping_reason, dict(thresholds or {}), tuple(processed_groups))
