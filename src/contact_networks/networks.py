# - Import graph, counting, and pathway utilities.
import csv
import re
from collections import Counter, defaultdict
from itertools import pairwise
from math import isfinite
import networkx as nx
from .models import Pathway, ResidueId

# - Validate supported detailed CONTACT headers before parsing.
DETAILED_HEADER = (
    "PATH(ModelNo, ChainID, SeqNum(InsertionCode), AltLoc (clash) --> AltLoc (relief), "
    "clash with prvs, relief, clash atom (prvs rsd), clash atom (current rsd)"
)

# - Finalize detailed records at boundaries and EOF.
def read_pathways(filename):
    with open(filename, encoding="utf-8", newline="") as stream:
        rows, started = [], False
        for line_number, line in enumerate(stream, 1):
            text = line.strip()
            # - Stream records without empty phantom pathways.
            # - Reject malformed pathway headers instead of guessing.
            if text.startswith("PATH") and " ".join(text.split()) != DETAILED_HEADER:
                raise ValueError(f"line {line_number}: unsupported detailed CONTACT header")
            if text.startswith("PATH(") or not text:
                if rows:
                    yield _pathway(rows)
                    rows = []
                started = text.startswith("PATH(")
                continue
            # - Reject records outside a detailed pathway.
            if not started:
                raise ValueError(f"line {line_number}: expected PATH header")
            try:
                rows.append(_record(next(csv.reader([line]))))
            except (ValueError, IndexError) as error:
                raise ValueError(f"line {line_number}: {error}") from error
        if rows:
            yield _pathway(rows)

# - Preserve chain, insertion, and conformer identity.
def _record(fields):
    if len(fields) != 8:
        raise ValueError("expected eight detailed CONTACT fields")
    model, chain, number, transition, clash, relief, _, _ = map(str.strip, fields)
    match = re.fullmatch(r"(-?\d+)\(([^()]*)\)", number)
    if match is None:
        raise ValueError("expected residue number(insertion-code)")
    identity = ResidueId(int(model), chain, int(match[1]), match[2])
    conformers = tuple(part.strip() for part in transition.split("-->"))
    # - Reject malformed and unaligned pathway fields.
    if len(conformers) != 2:
        raise ValueError("expected clashing --> relieving conformers")
    return identity, _number(relief), conformers, _number(clash)

# - Reject nonfinite or negative numeric measures.
def _number(value):
    result = float(value)
    if not isfinite(result) or result < 0:
        raise ValueError("measures must be finite and nonnegative")
    return result

# - Distinguish initial relief markers from zero relief.
def _pathway(rows):
    residues, reliefs, conformers, clashes = map(tuple, zip(*rows))
    if reliefs[0] != 0:
        raise ValueError("initial relief must be the zero start marker")
    return Pathway(residues, (None, *reliefs[1:]), conformers, (None, *clashes[1:]))

# - Apply relief, length, and residue filters.
def filter_pathways(paths, settings):
    selected = frozenset(settings.residue_filter)
    return tuple(
        path for path in paths
        if len(path.residues) > settings.min_pathway_length
        and all(flag or value >= settings.relief_threshold for value, flag in
                zip(path.reliefs[1:], (path.relief_exclusions or (False,) * len(path.reliefs))[1:]))
        and (not selected or any(item in selected for item in path.residues))
    )

# - Count residue occurrences before scientific filtering.
def residue_counts(paths):
    return Counter(item for path in paths for item in path.residues)

# - Aggregate repeated undirected transitions directly.
def build_graph(paths, min_network_size=3):
    if isinstance(min_network_size, bool) or not isinstance(min_network_size, int) or min_network_size < 1:
        raise ValueError("min_network_size must be a positive integer")
    graph = nx.Graph()
    for path in paths:
        graph.add_nodes_from(path.residues)
        for left, right in pairwise(path.residues):
            previous = graph.get_edge_data(left, right, {}).get("weight", 0)
            graph.add_edge(left, right, weight=previous + 1)
    # - Remove connected components below minimum size.
    small = [node for part in nx.connected_components(graph) if len(part) < min_network_size for node in part]
    graph.remove_nodes_from(small)
    return graph

# - Count pathways touching each retained network.
def network_pathway_counts(graph, paths):
    components = sorted(nx.connected_components(graph), key=lambda part: tuple(sorted(part)))
    membership = {node: index for index, part in enumerate(components) for node in part}
    counts = [0] * len(components)
    for path in paths:
        for index in {membership[node] for node in path.residues if node in membership}:
            counts[index] += 1
    return tuple(counts)

# - Remove duplicate ordered residue sequences.
def maximal_pathways(paths):
    unique = tuple(dict.fromkeys(path.residues for path in paths))
    positions = defaultdict(list)
    for sequence in unique:
        for index, node in enumerate(sequence):
            positions[node].append((sequence, index))
    # - Remove contiguous subpaths without reversing sequences.
    return tuple(
        sequence for sequence in unique
        if not any(len(other) > len(sequence) and other[index:index + len(sequence)] == sequence
                   for other, index in positions[sequence[0]])
    )
