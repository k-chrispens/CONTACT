# - Import structure readers and synthetic fixture utilities.
import gemmi
import pytest
from contact_networks import Settings
from contact_networks.structures import read_structure

# - Build small prepared multi-conformer structures.
def fixture(model_numbers=(1,)):
    structure = gemmi.Structure()
    for number in model_numbers:
        model = gemmi.Model(number)
        for name in ("A", "B"):
            chain = gemmi.Chain(name)
            residue = gemmi.Residue()
            residue.name, residue.seqid, residue.subchain = "ALA", gemmi.SeqId(7, "I"), name + "x"
            # - Add shared hydrogen and alternate carbon coordinates.
            for atom_name, element, altloc, x in (("H", "H", "\0", 0), ("CB", "C", "A", 1), ("CB", "C", "B", 2)):
                atom = gemmi.Atom()
                atom.name, atom.element, atom.altloc = atom_name, gemmi.Element(element), altloc
                atom.pos = gemmi.Position(x, 0, 0)
                residue.add_atom(atom)
            chain.add_residue(residue)
            model.add_chain(chain)
        structure.add_model(model)
    return structure

# - Serialize equivalent PDB and mmCIF fixtures.
def write(structure, file):
    if file.suffix == ".pdb":
        structure.write_pdb(str(file))
    else:
        structure.make_mmcif_document().write_file(str(file))
    return file

# - Verify equivalent PDB and mmCIF interpretations.
def test_formats_and_identifiers(tmp_path):
    structure = fixture()
    pdb = read_structure(write(structure, tmp_path / "x.pdb"), Settings())
    cif = read_structure(write(structure, tmp_path / "x.cif"), Settings())
    assert pdb == cif
    assert [(residue.id.chain, residue.id.number, residue.id.insertion) for residue in pdb] == [
        ("A", 7, "I"), ("B", 7, "I")]
    assert [(atom.name, atom.altloc) for atom in pdb[0].atoms] == [("H", ""), ("CB", "A"), ("CB", "B")]

# - Select deposited models rather than list offsets.
def test_models(tmp_path):
    path = write(fixture((2, 5)), tmp_path / "models.cif")
    assert {residue.id.model for residue in read_structure(path, Settings(model=5))} == {5}
    with pytest.raises(ValueError, match="model numbered 1"):
        read_structure(path, Settings())

# - Reject hydrogen-free input without explicit permission.
def test_hydrogen_policy(tmp_path):
    structure = fixture()
    for chain in structure[0]:
        del chain[0][0]
    path = write(structure, tmp_path / "no-h.cif")
    with pytest.raises(ValueError, match="hydrogen"):
        read_structure(path, Settings())
    assert len(read_structure(path, Settings(allow_hydrogen_free=True))) == 2

# - Reject nonfinite coordinates in either input format.
@pytest.mark.parametrize("value", [float("nan"), float("inf")])
def test_coordinates(tmp_path, value):
    structure = fixture()
    structure[0][0][0][0].pos.x = value
    path = write(structure, tmp_path / "bad.cif")
    with pytest.raises(ValueError, match="coordinates"):
        read_structure(path, Settings())

# - Preserve ligand identity while excluding water residues.
def test_ligands_and_water(tmp_path):
    structure = fixture()
    for number, name in ((100, "NAP"), (101, "HOH")):
        residue = structure[0][0][0].clone()
        residue.name, residue.seqid = name, gemmi.SeqId(number, " ")
        structure[0][0].add_residue(residue)
    residues = read_structure(write(structure, tmp_path / "ligands.cif"), Settings())
    assert {residue.name for residue in residues} == {"ALA", "NAP"}

# - Build explicit assembly operators and metadata.
def assembly_fixture(subchains=False):
    structure = fixture()
    assembly, generator = gemmi.Assembly("1"), gemmi.Assembly.Gen()
    if subchains:
        generator.subchains = ["Ax"]
    else:
        generator.chains = ["A"]
    # - Generate two distinct translated source-chain copies.
    for index, translation in enumerate((0, 10)):
        operator = gemmi.Assembly.Operator()
        operator.name = str(index + 1)
        operator.transform.vec.x = translation
        generator.operators.append(operator)
    assembly.generators.append(generator)
    structure.assemblies.append(assembly)
    return structure

# - Verify assembly copies preserve chains and transformations.
@pytest.mark.parametrize("extension,subchains", [("pdb", False), ("cif", True)])
def test_assembly(tmp_path, extension, subchains):
    path = write(assembly_fixture(subchains), tmp_path / f"assembly.{extension}")
    residues = read_structure(path, Settings(assembly="1"))
    assert len(residues) == 2
    assert {residue.id.chain for residue in residues} == {"A"}
    assert len({residue.id.copy for residue in residues}) == 2
    assert sorted(residue.atoms[0].position[0] for residue in residues) == [0, 10]

# - Reject missing assemblies without coordinate fallback.
def test_missing_assembly(tmp_path):
    path = write(fixture(), tmp_path / "plain.cif")
    with pytest.raises(ValueError, match="assembly"):
        read_structure(path, Settings(assembly="missing"))

# - Exclude heavy water without discarding cofactors.
def test_heavy_water(tmp_path):
    structure = fixture()
    residue = structure[0][0][0].clone()
    residue.name, residue.seqid = "DOD", gemmi.SeqId(101, " ")
    structure[0][0].add_residue(residue)
    selected = read_structure(write(structure, tmp_path / "heavy.cif"), Settings())
    assert all(residue.name != "DOD" for residue in selected)

# - Reject partially unmatched assembly coordinate selectors.
@pytest.mark.parametrize("subchains", [False, True])
def test_partial_assembly_selection(tmp_path, subchains):
    structure = assembly_fixture(subchains)
    generator = structure.assemblies[0].generators[0]
    if subchains:
        generator.subchains = ["Ax", "missing"]
    else:
        generator.chains = ["A", "Z"]
    # - Verify missing selectors cannot produce partial assemblies.
    path = write(structure, tmp_path / ("incomplete.cif" if subchains else "incomplete.pdb"))
    with pytest.raises(ValueError, match="selectors"):
        read_structure(path, Settings(assembly="1"))

# - Preserve split source-chain parts within shared operator copies.
def test_split_source_chain_parts_share_assembly_copy(tmp_path):
    structure = assembly_fixture()
    part = structure[0][0].clone()
    part[0].seqid, part[0].name = gemmi.SeqId(8, " "), "NAP"
    structure[0].add_chain(part)
    residues = read_structure(write(structure, tmp_path / "parts.pdb"), Settings(assembly="1"))
    assert len(residues) == 4
    assert len({residue.id.copy for residue in residues}) == 2
    assert len({residue.id for residue in residues}) == 4

# - Report malformed structure input through a consistent boundary.
def test_malformed_structure_errors(tmp_path):
    path = tmp_path / "bad.cif"
    path.write_text("not a valid coordinate data block")
    with pytest.raises(ValueError, match="prepared structure"):
        read_structure(path, Settings())

# - Reject empty coordinate documents through the input boundary.
@pytest.mark.parametrize("suffix", ["pdb", "cif"])
@pytest.mark.parametrize("contents", ["", " \n\t"])
def test_empty_coordinate_documents(tmp_path, suffix, contents):
    path = tmp_path / f"empty.{suffix}"
    path.write_text(contents)
    with pytest.raises(ValueError):
        read_structure(path, Settings())
