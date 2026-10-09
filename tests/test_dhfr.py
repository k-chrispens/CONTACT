# - Locate explicitly downloaded historical DHFR reference artifacts.
import os
from pathlib import Path
import networkx as nx
import pytest
from contact_networks import Settings, build_graph, filter_pathways, maximal_pathways, read_pathways

# - Check raw, maximal, and relief-filtered pathway counts.
def test_historical_dhfr_artifacts():
    directory = os.environ.get("CONTACT_DHFR_REFERENCE")
    if directory is None:
        pytest.skip("set CONTACT_DHFR_REFERENCE to downloaded reference directory")
    raw = tuple(read_pathways(Path(directory) / "pathways.txt"))
    accepted = filter_pathways(raw, Settings())
    assert len(raw) == 2415 and len(maximal_pathways(raw)) == 157
    assert len(accepted) == 340
    graph = build_graph(accepted, 3)

# - Check retained network membership at both minima.
    expected = [
        [7, 8, 14, 18, 44, 63, 64, 65, 98, 99, 102, 105, 113, 115, 125, 162],
        [4, 78, 79, 101, 103, 104, 106, 107, 108, 126, 127, 128, 158],
        [109, 112, 133, 135, 137, 154, 155, 156], [138, 139, 141, 153],
        [5, 94, 161], [25, 26, 147], [34, 36, 38], [40, 52, 54], [58, 60, 84],
    ]
    observed = [sorted(node.number for node in part) for part in nx.connected_components(graph)]
    assert sorted(observed) == sorted(expected)
    assert len(graph) == 56 and graph.number_of_edges() == 52

# - Keep integration tests offline and opt-in.
    larger = build_graph(accepted, 4)
    assert len(larger) == 41
    assert sorted(map(len, nx.connected_components(larger))) == [4, 8, 13, 16]
