# - Import contact construction and fixture records.
from math import dist
from random import Random
import pytest
from contact_networks.contacts import ContactMap, Conformer, build_contacts, _hbond, _radius
from contact_networks.models import Atom, Residue, ResidueId, Settings

# - Build small immutable atom and residue fixtures.
def atom(name, element, x, alt="", y=0, z=0):
    return Atom(name, element, (x, y, z), alt)

# - Create residues with stable identities.
def residue(number, name, atoms, chain="A", copy=""):
    return Residue(ResidueId(1, chain, number, copy=copy), name, tuple(atoms))

# - Make two-state residues around a coordinate.
def mobile(number, x, name="ALA", chain="A", atoms=()):
    base = list(atoms)
    base.extend((atom("CB", "C", x, "A"), atom("CB", "C", x, "B")))
    return residue(number, name, base, chain)

# - Confirm mobile states and shared atoms.
def test_mobile_states_include_shared_atoms():
    shared = atom("C", "C", 0)
    result = build_contacts((mobile(1, 0, atoms=(shared,)),), Settings())
    assert [node.altloc for node in result.nodes] == ["A", "B"]
    assert all(shared in node.atoms for node in result.nodes)
    assert result.neighbors == ((), ())

# - Drop states made mobile only by backbone atoms.
def test_side_chain_filter_precedes_altloc_collection():
    atoms = (atom("CA", "C", 0, "A"), atom("CA", "C", 0, "B"))
    result = build_contacts((residue(1, "ALA", atoms),), Settings(side_chain_only=True))
    assert result.nodes == () and result.cutoff == 0

# - Keep ligand backbone-like names outside protein classification.
def test_ligand_names_are_not_backbone():
    lig = mobile(1, 0, "LIG", atoms=(atom("C", "C", 0), atom("O", "O", 0)))
    result = build_contacts((lig,), Settings(side_chain_only=True))
    assert len(result.nodes) == 2

# - Measure distant interactions on demand symmetrically.
def test_distant_ratio_is_measured_but_not_adjacent():
    left, right = mobile(1, 0), mobile(8, 30)
    result = build_contacts((left, right), Settings())
    assert result.ratio(0, 2) == result.ratio(2, 0)
    assert result.ratio(0, 2) > 1
    assert result.neighbors == ((), (), (), ())

# - Distinguish excluded cysteine interactions from zero.
def test_all_cysteine_residue_pairs_are_excluded():
    left = mobile(1, 0, "CYS", atoms=(atom("SG", "S", 0, "A"), atom("SG", "S", 0, "B")))
    right = mobile(9, 0, "CYS", atoms=(atom("SG", "S", 0, "A"), atom("SG", "S", 0, "B")))
    result = build_contacts((left, right), Settings())
    assert result.ratio(0, 2) is None

# - Preserve real measured zero for unrelated residues.
def test_zero_distance_is_not_an_exclusion():
    result = build_contacts((mobile(1, 0), mobile(9, 0)), Settings())
    assert result.ratio(0, 2) == 0

# - Check exact radius source and hydrogen corrections.
def test_radius_values_and_named_hydrogen_parity():
    lys = residue(1, "LYS", (atom("HZ2", "H", 0),))
    assert _radius(lys, lys.atoms[0]) == 1
    assert _radius(lys, atom("HZ3", "D", 0)) == 1
    assert _radius(lys, atom("HZZ", "H", 0)) == 1.17
    assert _radius(lys, atom("C", "C", 0)) == 1.65

# - Preserve isotope parity and alpha-hydrogen radii.
def test_isotope_and_backbone_hydrogens():
    ala = residue(1, "ALA", (atom("H", "H", 0),))
    assert _radius(ala, atom("D", "D", 0)) == 1
    assert _radius(ala, atom("DA2", "D", 0)) == 1.17
    assert _radius(ala, atom("HA2", "H", 0)) == 1.17
    for name in ("DHZ", "DZ1", "DZ2", "DZ3"):
        lys = residue(1, "LYS", (atom(name, "D", 0),))
        assert _radius(lys, lys.atoms[0]) == (1 if name != "DHZ" else 1.17)

# - Reject unsupported elemental radii explicitly.
def test_unknown_element_fails():
    bad = mobile(1, 0, atoms=(atom("XX", "Xx", 0),))
    with pytest.raises(ValueError, match="unsupported element"):
        build_contacts((bad,), Settings())

# - Exclude adjacent amino-acid backbone contacts only.
def test_adjacent_backbone_mask_respects_chain():
    first = residue(1, "ALA", (atom("N", "N", 0, "A"), atom("N", "N", 0, "B")))
    next_res = residue(2, "ALA", (atom("N", "N", 0, "A"), atom("N", "N", 0, "B")))
    other_chain = residue(1, "ALA", (atom("N", "N", 0, "A"), atom("N", "N", 0, "B")), "B")
    result = build_contacts((first, next_res, other_chain), Settings())
    assert result.ratio(0, 2) is None
    assert result.ratio(0, 4) == 0

# - Keep contact thresholds inclusive and percentile bounded.
def test_percentile_endpoints_and_tied_contacts():
    residues = (mobile(1, 0), mobile(9, 2), mobile(17, 2))
    low = build_contacts(residues, Settings(clash_percentile=0))
    high = build_contacts(residues, Settings(clash_percentile=1))
    assert low.cutoff <= high.cutoff
    assert all(i in low.neighbors[j] for j, row in enumerate(low.neighbors) for i in row)
    assert all(i in high.neighbors[j] for j, row in enumerate(high.neighbors) for i in row)

# - Recognize both hydrogen-bond orientations and boundaries.
def test_hbond_geometry_and_reverse_orientation():
    donor = mobile(1, 0, atoms=(atom("N", "N", 0),))
    acceptor = mobile(8, 0, atoms=(atom("O", "O", 3.5), atom("C", "C", 4.5)))
    assert _hbond(build_contacts((donor,), Settings()).nodes[0], build_contacts((acceptor,), Settings()).nodes[0])
    assert not _hbond(build_contacts((acceptor,), Settings()).nodes[0], build_contacts((donor,), Settings()).nodes[0])

# - Build independent backbone-only geometry classification fixtures.
def bond_pair(distance, carbon):
    left = residue(1, "ALA", (atom("N", "N", 0),))
    right = residue(8, "ALA", (atom("O", "O", distance), Atom("C", "C", carbon)))
    return Conformer(left, "A", left.atoms), Conformer(right, "A", right.atoms)

# - Test precise distance, angle, and degeneracy boundaries.
@pytest.mark.parametrize("distance,carbon,expected", [
    (3.5, (4.5, 0, 0), True), (3.50001, (4.5, 0, 0), False),
    (3.5, (3.5, 1, 0), True), (3.5, (2.5, 0, 0), False),
    (0, (1, 0, 0), False), (3.5, (3.5, 0, 0), False)])
def test_hbond_exact_boundaries(distance, carbon, expected):
    donor, acceptor = bond_pair(distance, carbon)
    assert _hbond(donor, acceptor) is expected

# - Check chemical masking in both canonical index orders.
def test_hbond_mask_is_symmetric():
    donor, acceptor = bond_pair(3.5, (4.5, 0, 0))
    for nodes in ((donor, acceptor), (acceptor, donor)):
        contact = ContactMap(nodes, 0, ((), ()))
        assert contact.ratio(0, 1) is None
        assert contact.ratio(1, 0) is None

# - Generate deterministic cells spanning positive and negative coordinates.
def random_residues():
    random, residues = Random(413), []
    for number in range(1, 13):
        center = tuple(random.uniform(-10, 10) for _ in range(3))
        if number <= 2:
            center = (2 * (number - 1), 0, 0)
        # - Include shared atoms and distinct selected conformers.
        atoms = (Atom("C0", "C", center),
                 Atom("C1", "C", (center[0] - 0.3, center[1], center[2]), "A"),
                 Atom("C1", "C", (center[0] + 0.3, center[1], center[2]), "B"))
        residues.append(residue(number, "LIG", atoms))
    return tuple(residues)

# - Compare spatial discovery against independent exhaustive geometry.
@pytest.mark.parametrize("percentile", [0, 0.3, 1])
def test_exact_spatial_candidate_oracle(percentile):
    inputs = random_residues()
    contact = build_contacts(inputs, Settings(clash_percentile=percentile))
    nodes = [(item.id, alt, tuple(a for a in item.atoms if not a.altloc or a.altloc == alt))
             for item in inputs for alt in ("A", "B")]
    assert [(node.residue.id, node.altloc) for node in contact.nodes] == [(rid, alt) for rid, alt, _ in nodes]
    ratios = {}
    # - Evaluate every distinct-residue pair without spatial pruning.
    for left, (rid, _, atoms) in enumerate(nodes):
        for right in range(left + 1, len(nodes)):
            if nodes[right][0] != rid:
                ratios[left, right] = min(dist(a.position, b.position) / 3.5 for a in atoms for b in nodes[right][2])
    sample = sorted(value for value in ratios.values() if value < 1.0001)
    cutoff = sample[min(int(percentile * len(sample)), len(sample) - 1)]
    assert contact.cutoff == pytest.approx(cutoff)
    # - Compare inclusive threshold adjacency and distant relief values.
    expected = {key for key, value in ratios.items() if value <= cutoff}
    observed = {(left, right) for left, row in enumerate(contact.neighbors) for right in row if left < right}
    assert observed == expected
    for (left, right), value in ratios.items():
        assert contact.ratio(left, right) == pytest.approx(value)
