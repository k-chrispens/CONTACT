# - Write independently specified prepared transition cascades.
# - Compare structure analysis with explicit expected network counts.
# - Validate known science parameters and group provenance.
# - Preserve partial results and exact budget stopping behavior.
# - Compare equivalent prepared PDB and mmCIF analyses.
from contact_networks import Settings, analyze, to_document
import gemmi


# - Write a prepared three-residue alternate-state cascade.
def write_cascade(path):
    rows = []
    serial = 1
    # - Give each residue two alternate conformers.
    for residue, number, states in (
        ("LIG", 1, (("A", 0), ("B", -10))),
        ("CB", 2, (("A", 2), ("B", 10))),
        ("CB", 3, (("A", 12), ("B", 30))),
    ):
        # - Write carbon-only atoms in fixed-column format.
        for altloc, x in states:
            rows.append(
                f"HETATM{serial:5d} {'CB':>4}{altloc}{residue:>3} A{number:4d}    "
                f"{x:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}           C\n"
            )
            serial += 1
    path.write_text("".join(rows) + "END\n")
    return path


# - Check structure provenance and generated network measures.
def test_structure_pipeline_and_budget(tmp_path):
    source = write_cascade(tmp_path / "cascade.pdb")
    settings = Settings(allow_hydrogen_free=True)
    result = analyze(source, settings)
    document = to_document(result)
    cif = tmp_path / "cascade.cif"
    gemmi.read_structure(str(source)).make_mmcif_document().write_file(str(cif))
    # - Ensure coordinate formats yield matching counts.
    assert to_document(analyze(cif, settings))["counts"] == document["counts"]
    assert document["counts"]["raw_pathways"] == 4
    assert document["counts"]["accepted_pathways"] == 4
    assert document["counts"]["maximal_raw_pathways"] == 2
    assert len(document["nodes"]) == 3
    assert [edge["weight"] for edge in document["edges"]] == [3, 3]
    assert document["candidate_states"] == 6
    assert document["generation_settings_known"] is True
    assert document["status_scope"] == "generation"
    assert document["science_revision"] == "contact-corrected-v1"
    # - Verify the candidate-state budget preserves partial truth.
    partial = to_document(analyze(source, Settings(allow_hydrogen_free=True, max_candidates=1)))
    assert partial["status"] == "truncated"
    assert partial["candidate_states"] == 1
    assert partial["counts"]["raw_pathways"] == 0
