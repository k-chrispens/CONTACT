# - Import reusable types and network operations.
import pytest
from contact_networks import Pathway, ResidueId, Settings
from contact_networks.networks import build_graph, filter_pathways, maximal_pathways
from contact_networks.networks import network_pathway_counts, read_pathways, residue_counts

# - Construct independent residue and pathway examples.
def rid(n, chain="A", insertion=""):
    return ResidueId(1, chain, n, insertion)

def path(*numbers, relief=1.0):
    return Pathway(tuple(rid(n) for n in numbers), (None,) + (relief,) * (len(numbers) - 1))

# - Build actual legacy detailed-pathway input records.
HEADER = ("PATH(ModelNo, ChainID, SeqNum(InsertionCode), AltLoc (clash) --> AltLoc (relief), "
          "clash with prvs, relief, clash atom (prvs rsd), clash atom (current rsd)\n")

def record(number, relief=0, chain="A", insertion=""):
    return f"1, {chain}, {number}({insertion}), A --> B, 0.8, {relief}, CA, CB\n"

def read_text(tmp_path, text):
    file = tmp_path / "pathways.txt"
    file.write_text(text)
    return tuple(read_pathways(file))

# - Check threshold equality and genuine zero relief.
def test_whole_path_relief_threshold():
    equal = path(1, 2, relief=0.9)
    bad = Pathway((rid(1), rid(2), rid(3)), (None, 0.99, 0.1))
    assert filter_pathways([equal, bad, path(4, 5, relief=0)], Settings()) == (equal,)
    assert filter_pathways([path(1, 2, relief=0)], Settings(relief_threshold=0)) == (path(1, 2, relief=0),)

# - Check EOF and consecutive empty pathway boundaries.
def test_actual_records_and_eof(tmp_path):
    rows = read_text(tmp_path, "\n" + HEADER + HEADER + record(1) + record(2, 0.9).rstrip())
    assert len(rows) == 1
    assert rows[0].residues == (rid(1), rid(2))
    assert rows[0].reliefs == (None, 0.9)
    assert rows[0].conformers == (("A", "B"), ("A", "B"))
    assert read_text(tmp_path, HEADER + "\n\n" + HEADER) == ()

# - Reject malformed records and nonfinite numeric fields.
@pytest.mark.parametrize("text", [record(1), HEADER + "garbage\n", HEADER + record(1, float("nan")),
                                      HEADER + record(1, 0) + record(2, -1)])
def test_malformed_records(tmp_path, text):
    with pytest.raises(ValueError):
        read_text(tmp_path, text)

# - Keep chain and insertion-code identities separate.
def test_chain_and_insertion_parsing(tmp_path):
    paths = read_text(tmp_path, HEADER + record(-12, insertion="B") + record(2, 0.95, "AB"))
    assert paths[0].residues == (rid(-12, insertion="B"), rid(2, "AB"))
    graph = build_graph([Pathway((rid(2, "A"), rid(2, "AB")), (None, 1))], 1)
    assert graph.number_of_nodes() == 2

# - Check repeated edge weights and network boundaries.
def test_edge_counts_and_component_sizes():
    paths = [path(1, 2, 3), path(3, 2, 1), path(8, 9)]
    graph = build_graph(paths, 3)
    assert graph.number_of_nodes() == 3
    assert graph[rid(1)][rid(2)]["weight"] == 2
    assert build_graph(paths, 4).number_of_nodes() == 0
    assert build_graph([path(1)], 1).number_of_nodes() == 1
    assert residue_counts(paths)[rid(2)] == 2
    assert network_pathway_counts(graph, iter(paths)) == (2,)

# - Apply strict length and selected-residue filters.
def test_length_and_residue_filter():
    short, long = path(1, 2), path(1, 2, 3)
    assert filter_pathways([short, long], Settings(min_pathway_length=2)) == (long,)
    assert filter_pathways([short, long], Settings(residue_filter=(rid(3),))) == (long,)

# - Preserve ordered contiguous-subpath elimination.
def test_ordered_subpaths():
    long, reversed_path = path(1, 2, 3), path(3, 2, 1)
    assert maximal_pathways([path(1, 2), long, long, reversed_path, path(1, 3)]) == (
        long.residues, reversed_path.residues, path(1, 3).residues)

# - Reject invalid scientific settings and records.
@pytest.mark.parametrize("kwargs", [{"max_pathway_length": True}, {"max_candidates": 0},
    {"clash_percentile": float("nan")}, {"clash_percentile": 1.1}, {"assembly": ""}])
def test_invalid_settings(kwargs):
    with pytest.raises(ValueError):
        Settings(**kwargs)

# - Reject unaligned mutable and nonfinite pathways.
def test_invalid_pathways():
    with pytest.raises(ValueError):
        Pathway((rid(1), rid(2)), (None, float("inf")))
    with pytest.raises(TypeError):
        Pathway([rid(1)], (None,))
    with pytest.raises(ValueError):
        Pathway((rid(1), rid(2)), (None,))

# - Reject mutable nested atom and pathway values.
def test_deeply_immutable_records():
    from contact_networks import Atom, Residue
    with pytest.raises(TypeError):
        Atom("CA", "C", [0, 0, 0])
    atom = Atom("CA", "C", (0, 0, 0))
    with pytest.raises(TypeError):
        Residue(rid(1), "ALA", [atom])
    with pytest.raises(ValueError):
        Pathway((rid(1),), (None,), (["A", "B"],))

# - Reject malformed or residue-only pathway headers.
@pytest.mark.parametrize("header", ["PATH(not-the-legacy-header)", "PATH(ModelNo, ChainID, SeqNum(InsertionCode), AltLoc):"])
def test_invalid_headers(tmp_path, header):
    with pytest.raises(ValueError, match="header"):
        read_text(tmp_path, header + "\n" + record(1))
