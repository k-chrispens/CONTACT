# - Run CLI against representative detailed pathways.
import json
import subprocess
import sys
from contact_networks import Settings, analyze_pathways, read_pathways, to_document
from contact_networks.cli import main
from test_networks import HEADER, record

# - Verify matching API and command-line results.
def test_api_cli_equivalence(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text(HEADER + record(1) + record(2, 0.9) + record(3, 0.99))
    output = tmp_path / "result"
    assert main(["extract", str(source), "--output", str(output)]) == 0
    expected = analyze_pathways(read_pathways(source), Settings(), source=str(source))
    assert json.loads((output / "results.json").read_text()) == to_document(expected)

# - Verify explicit output destination and exit statuses.
def test_partial_and_failure_statuses(tmp_path):
    source, output = tmp_path / "input.txt", tmp_path / "partial"
    source.write_text(HEADER + record(1) + record(2, 0.99))
    args = ["extract", str(source), "--output", str(output), "--input-status", "truncated",
            "--stopping-reason", "upstream enumeration limit"]
    assert main(args) == 3
    assert json.loads((output / "results.json").read_text())["status"] == "truncated"
    assert main(args) == 1

# - Verify malformed inputs produce ordinary failure statuses.
def test_invalid_input_does_not_finalize(tmp_path):
    source, output = tmp_path / "bad.txt", tmp_path / "result"
    source.write_text("not detailed CONTACT input")
    assert main(["extract", str(source), "--output", str(output)]) == 1
    assert not output.exists()

# - Verify imports perform no analysis or rendering.
def test_import_and_module_entry_point():
    code = "import sys, contact_networks; assert 'matplotlib' not in sys.modules"
    assert subprocess.run([sys.executable, "-c", code], capture_output=True).returncode == 0
    command = subprocess.run([sys.executable, "-m", "contact_networks", "--help"], capture_output=True)
    assert command.returncode == 0 and b"extract" in command.stdout

# - Preserve source truncation when model selection is empty.
def test_partial_source_without_selected_model(tmp_path):
    source, output = tmp_path / "other-model.txt", tmp_path / "result"
    source.write_text((HEADER + record(1) + record(2, 0.99)).replace("1, A,", "2, A,"))
    args = ["extract", str(source), "--output", str(output), "--model", "1",
            "--input-status", "truncated", "--stopping-reason", "upstream limit"]
    assert main(args) == 3
    document = json.loads((output / "results.json").read_text())
    assert document["counts"]["raw_pathways"] == 0
    assert document["status"] == "truncated" and document["status_scope"] == "source_input"

# - Run prepared models through corrected structure analysis.
# - Preserve generated partial status and ordinary input failures.
from test_analysis import write_cascade


# - Export generated structure analyses with truthful status.
def test_structure_cli_and_partial_export(tmp_path):
    source = write_cascade(tmp_path / "cascade.pdb")
    output = tmp_path / "generated"
    args = ["analyze", str(source), "--allow-hydrogen-free", "--output", str(output)]
    assert main(args) == 0
    document = json.loads((output / "results.json").read_text())
    assert document["input_kind"] == "structure"
    assert document["generation_settings_known"] is True
    # - Confirm partial generation exports and signals truncation.
    partial = tmp_path / "partial-generated"
    partial_args = ["analyze", str(source), "--allow-hydrogen-free",
                    "--max-candidates", "1", "--output", str(partial)]
    assert main(partial_args) == 3
    data = json.loads((partial / "results.json").read_text())
    assert data["candidate_states"] == 1
    assert data["status_scope"] == "generation"
    assert data["stopping_reason"]


# - Reject malformed and missing generated inputs without export.
def test_structure_cli_errors_do_not_export(tmp_path):
    bad = tmp_path / "bad.pdb"
    bad.write_text("invalid structure\n")
    for source in (bad, tmp_path / "missing.pdb"):
        output = tmp_path / (source.stem + "-out")
        assert main(["analyze", str(source), "--output", str(output)]) == 1
        assert not output.exists()
