# - Import real records, serialization, and concurrency utilities.
import csv
import json
from concurrent.futures import ThreadPoolExecutor

import networkx as nx
import pytest

from contact_networks.analysis import Analysis
from contact_networks.models import Pathway, ResidueId, Settings
from contact_networks.networks import maximal_pathways
from contact_networks.output import to_document, write_results


# - Construct distinct chain and assembly-copy identities.
def residue(chain, number, copy=""):
    return ResidueId(1, chain, number, copy=copy)


# - Construct immutable aligned pathway measurement records.
def pathway(*residues, reliefs=None, conformers=(), clashes=()):
    reliefs = reliefs if reliefs is not None else (None,) + (0.0,) * (len(residues) - 1)
    return Pathway(tuple(residues), reliefs, conformers, clashes)


# - Construct explicitly supplied analysis result records.
def analysis(raw, accepted, graph):
    return Analysis("input.pdb", Settings(), tuple(raw), tuple(accepted), graph,
                    maximal_pathways(raw), thresholds={"relief": 0.0}, processed_groups=("1",))


# - Exercise chain collisions, zero relief, and copies.
def fixture_result():
    a, b, c, d = residue("A", 7), residue("B", 7), residue("A", 7, "copy-2"), residue("B", 8)
    first = pathway(a, b, reliefs=(None, 0.0), conformers=(("", ""), ("clash", "relieve")),
                    clashes=(None, 2.5))
    second = pathway(c, d, reliefs=(None, 1.25), clashes=(None, 3.0))
    graph = nx.Graph()
    graph.add_edge(a, b, weight=4)
    graph.add_edge(c, d, weight=2)
    return analysis((first, second), (first,), graph)


# - Verify JSON counts and full identity mappings.
def test_document_and_csv_preserve_counts_identity_and_values(tmp_path):
    result = fixture_result()
    document = to_document(result)
    output = write_results(result, tmp_path / "result")
    saved = json.loads((output / "results.json").read_text())
    assert saved == document
    assert document["counts"] == {"raw_pathways": 2, "accepted_pathways": 1,
                                  "maximal_raw_pathways": 2, "raw_residue_occurrences": 4}
    # - Verify exact weighted edge and residue identities.
    assert [(n["chain"], n["copy"]) for n in document["nodes"]] == [
        ("A", ""), ("A", "copy-2"), ("B", ""), ("B", "")]
    assert [(e["source"], e["target"], e["weight"]) for e in document["edges"]] == [(1, 3, 4), (2, 4, 2)]
    with (output / "raw_pathways.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    # - Preserve missing values and genuine numeric zero.
    assert rows[0]["relief"] == "" and rows[1]["relief"] == "0.0"
    assert rows[1]["clash"] == "2.5"
    assert json.loads(rows[1]["conformers"]) == ["clash", "relieve"]
    assert json.loads(rows[0]["residue_id"]) == {
        "model": 1, "chain": "A", "number": 7, "insertion": "", "copy": ""}


# - Verify repeatability and existing destination refusal.
def test_repeatability_and_existing_destination_refusal(tmp_path):
    result = fixture_result()
    first = write_results(result, tmp_path / "one")
    second = write_results(result, tmp_path / "two")
    assert {p.name: p.read_bytes() for p in first.iterdir()} == {
        p.name: p.read_bytes() for p in second.iterdir()}
    with pytest.raises(FileExistsError):
        write_results(result, first)


# - Verify one writer wins concurrent destination reservation.
def test_concurrent_destination_collision_has_one_winner(tmp_path):
    result = fixture_result()
    destination = tmp_path / "shared"

    def attempt(_):
        try:
            write_results(result, destination)
            return "written"
        except FileExistsError:
            return "refused"

    # - Check the winning writer finalizes complete outputs.
    with ThreadPoolExecutor(max_workers=2) as pool:
        outcomes = list(pool.map(attempt, range(2)))
    assert sorted(outcomes) == ["refused", "written"]
    assert (destination / "results.json").exists()


# - Verify more than nine plots retain unique filenames.
def test_plots_number_more_than_nine_networks(tmp_path):
    pytest.importorskip("matplotlib")
    graph = nx.Graph()
    for index in range(10):
        graph.add_node(residue("A", index + 1))
    result = analysis((), (), graph)
    output = write_results(result, tmp_path / "plots", plots=True)
    assert sorted(p.name for p in output.glob("network_*.png")) == [
        f"network_{index:04d}.png" for index in range(1, 11)]
    assert (output / "results.json").exists()


# - Inject reachable export failure before finalization.
def test_failed_export_never_writes_completion_manifest(tmp_path, monkeypatch):
    import contact_networks.output as output_module

    def fail(*args, **kwargs):
        raise OSError("simulated export failure")

    # - Verify failed exports never expose completion manifests.
    monkeypatch.setattr(output_module, "_write_csv", fail)
    destination = tmp_path / "failed"
    with pytest.raises(OSError, match="simulated"):
        write_results(fixture_result(), destination)
    assert destination.exists()
    assert not (destination / "results.json").exists()
