# - Build small independently checkable transition geometries.
from math import dist
import pytest
from contact_networks import Atom, Residue, ResidueId, Settings
from contact_networks.engine import generate_pathways

def mobile(number, positions, chain="A"):
    atoms = tuple(Atom("CB", "C", (x, 0, 0), altloc) for altloc, x in zip("ABC", positions))
    return Residue(ResidueId(1, chain, number), "LIG", atoms)

# - Create a directed three-residue clash-relief cascade.
def cascade(chain="A", offset=0, branching=False):
    middle = (2 + offset, 10 + offset, 4 + offset) if branching else (2 + offset, 10 + offset)
    return (mobile(1, (offset, offset - 10), chain), mobile(2, middle, chain),
            mobile(3, (12 + offset, 30 + offset), chain))

def generate(residues, **settings):
    return generate_pathways(residues, Settings(allow_hydrogen_free=True, **settings))

# - Verify source switches introduce rather than relieve clashes.
def test_directed_cascade():
    result = generate(cascade())
    assert result.status == "complete" and result.candidate_states == 6
    assert len(result.pathways) == 4
    assert {(tuple(node.number for node in path.residues)) for path in result.pathways} == {
        (1, 2, 3), (2, 1), (2, 3), (3, 2, 1)}

# - Verify propagation respects before-and-after overlap inequalities.
def test_emitted_transitions_obey_overlap_direction():
    residues = cascade(branching=True)
    coordinates = {(residue.id, atom.altloc): atom.position for residue in residues for atom in residue.atoms}
    for path in generate(residues).pathways:
        for step in range(1, len(path.residues) - 1):
            before, after = path.conformers[step]
            next_state = path.conformers[step + 1][0]
            old = dist(coordinates[path.residues[step], before], coordinates[path.residues[step + 1], next_state])
            new = dist(coordinates[path.residues[step], after], coordinates[path.residues[step + 1], next_state])
            assert old >= new

# - Recompute branch-local scores from emitted conformer choices.
def test_branch_local_relief_scores():
    residues = cascade(branching=True)
    coordinates = {(residue.id, atom.altloc): atom.position for residue in residues for atom in residue.atoms}
    for path in generate(residues).pathways:
        for step in range(1, len(path.residues)):
            active = path.conformers[step - 1][0 if step == 1 else 1]
            previous = coordinates[path.residues[step - 1], active]
            before, after = (coordinates[path.residues[step], label] for label in path.conformers[step])
            assert path.clashes[step] == pytest.approx(dist(previous, before) / 3.5)
            assert path.reliefs[step] == pytest.approx(dist(previous, after) / 3.5)

# - Reject terminal alternatives conflicting with earlier sources.
def test_all_prior_states_are_relieved():
    residues = cascade(branching=True)
    coordinates = {(residue.id, atom.altloc): atom.position for residue in residues for atom in residue.atoms}
    for path in generate(residues).pathways:
        for step in range(1, len(path.residues)):
            target = coordinates[path.residues[step], path.conformers[step][1]]
            for earlier in range(step):
                label = path.conformers[earlier][0 if earlier == 0 else 1]
                assert dist(coordinates[path.residues[earlier], label], target) / 3.5 > 2 / 3.5

# - Check ordered residue identities and absence of cycles.
def test_unique_residues_and_selected_model():
    for path in generate(cascade(branching=True)).pathways:
        assert len(set(path.residues)) == len(path.residues)
    with pytest.raises(ValueError, match="selected model"):
        generate_pathways(cascade(), Settings(model=2, allow_hydrogen_free=True))
    with pytest.raises(ValueError, match="unique"):
        generate((*cascade(), cascade()[0]))

# - Check transition limits and exact candidate-budget boundaries.
def test_limits_and_exact_completion_boundary():
    limited = generate(cascade(), max_pathway_length=1)
    assert len(limited.pathways) == 2
    assert all(len(path.residues) == 2 for path in limited.pathways)
    assert generate(cascade(), max_candidates=6).status == "complete"
    partial = generate(cascade(), max_candidates=1)
    assert partial.status == "truncated" and partial.candidate_states == 1
    assert partial.pathways == () and partial.processed_groups == ()
    assert "max_candidates=1" in partial.stopping_reason

# - Verify deterministic partial results across multiple groups.
def test_group_budget_and_repeatability():
    residues = cascade("A") + cascade("B", offset=100)
    first = generate(residues, max_candidates=6)
    second = generate(tuple(reversed(residues)), max_candidates=6)
    assert first == second
    assert first.status == "truncated" and len(first.processed_groups) == 1
    assert len(first.pathways) == 4 and first.candidate_states == 6
    assert '"B"' in first.stopping_reason

# - Check chain-local and cross-chain generation independently.
def test_cross_chain_scope():
    residues = (cascade("A")[0], *cascade("B")[1:])
    local = generate(residues)
    assert all(len({node.chain for node in path.residues}) == 1 for path in local.pathways)
    cross = generate(residues, cross_chain=True)
    assert any(len({node.chain for node in path.residues}) == 2 for path in cross.pathways)

# - Preserve excluded terminal contacts without invented numeric relief.
def test_chemically_excluded_terminal_score():
    donor = Residue(ResidueId(1, "A", 1), "ALA", (
        Atom("N", "N", (0, 0, 0), "A"), Atom("N", "N", (-10, 0, 0), "B")))
    acceptor = Residue(ResidueId(1, "A", 4), "ALA", (
        Atom("O", "O", (1, 0, 0), "A"), Atom("C", "C", (0, 0, 0), "A"),
        Atom("O", "O", (2, 0, 0), "B"), Atom("C", "C", (3, 0, 0), "B")))
    paths = generate((donor, acceptor)).pathways
    assert any(path.reliefs[-1] is None and path.relief_exclusions[-1] for path in paths)

# - Require explicit opt-in for missing nuclear hydrogens.
def test_hydrogen_policy_and_no_clashes():
    with pytest.raises(ValueError, match="hydrogen"):
        generate_pathways(cascade())
    assert generate((cascade()[0],)).pathways == ()
    assert generate(()).status == "complete"
