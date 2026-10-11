# - Import shared scientific and graph operations.
from dataclasses import dataclass, field, replace
import networkx as nx
from .models import Pathway, ResidueId, Settings
from .networks import build_graph, filter_pathways, maximal_pathways
from .engine import generate_pathways
from .structures import read_structure

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
    candidate_states: int | None = None
    # - Prevent incomplete results masquerading as complete analyses.
    def __post_init__(self):
        if self.status not in ("complete", "truncated"):
            raise ValueError("status must be complete or truncated")
        if (self.status == "truncated") != bool(self.stopping_reason):
            raise ValueError("only truncated results require a stopping reason")
        # - Validate input-specific generation provenance.
        if self.input_kind not in ("detailed_pathways", "structure"):
            raise ValueError("input_kind must be detailed_pathways or structure")
        if self.input_kind == "structure":
            if type(self.candidate_states) is not int or self.candidate_states < 0:
                raise ValueError("structure analyses require nonnegative candidate_states")
        elif self.candidate_states is not None:
            raise ValueError("imported pathway analyses have unknown candidate_states")

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

# - Analyze a structure while preserving generation provenance.
def analyze(filename, settings=None):
    settings = Settings() if settings is None else settings
    generation = generate_pathways(read_structure(filename, settings), settings)
    result = analyze_pathways(generation.pathways, settings, source=str(filename),
                              status=generation.status, stopping_reason=generation.stopping_reason,
                              thresholds=generation.thresholds, processed_groups=generation.processed_groups)
    return replace(result, input_kind="structure", candidate_states=generation.candidate_states)
