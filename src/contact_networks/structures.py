# - Import structure reading and immutable residue records.
from math import isfinite
import gemmi
from .models import Atom, Residue, ResidueId

# - Normalize alternate locations and insertion codes.
def _empty(value):
    return "" if value in ("", " ", "\0") else value

# - Read PDB and mmCIF with identical semantics.
def read_structure(filename, settings):
    # - Normalize malformed vendor input errors at the boundary.
    # - Normalize empty vendor documents without masking application bugs.
    try:
        structure = gemmi.read_structure(str(filename))
    except (RuntimeError, ValueError, IndexError) as error:
        raise ValueError(f"unable to read prepared structure: {error}") from error
    models = [model for model in structure if model.num == settings.model]
    if len(models) != 1:
        raise ValueError(f"expected one model numbered {settings.model}")
    model = models[0]
    copies = _copies(structure, model, settings.assembly)
    residues, seen = [], set()
    # - Preserve original chains and generated copy identity.
    for chain, copy, transform, selected in copies:
        for group in selected:
            # - Exclude all recognized water including heavy water.
            if group.is_water():
                continue
            identity = ResidueId(model.num, chain.name, group.seqid.num, _empty(group.seqid.icode), copy)
            if identity in seen:
                raise ValueError(f"duplicate residue identity: {identity}")
            seen.add(identity)
            atoms = _atoms(group, transform)
            residues.append(Residue(identity, group.name, atoms))
    # - Reject hydrogen-free inputs unless explicitly allowed.
    if not residues:
        raise ValueError("selected structure contains no non-water residues")
    has_hydrogen = any(atom.element in ("H", "D") for residue in residues for atom in residue.atoms)
    if not has_hydrogen and not settings.allow_hydrogen_free:
        raise ValueError("structure contains no hydrogen atoms; opt in explicitly")
    return tuple(sorted(residues, key=lambda residue: residue.id))

# - Select explicit biological assemblies without fallback.
def _copies(structure, model, assembly_id):
    if assembly_id is None:
        return [(chain, "", None, tuple(chain)) for chain in model]
    assemblies = [item for item in structure.assemblies if item.name == assembly_id]
    if len(assemblies) != 1:
        raise ValueError(f"expected one assembly named {assembly_id!r}")
    copies = []
    # - Reject unmatched explicit assembly chain and subchain selectors.
    chains = {chain.name for chain in model}
    subchains = {residue.subchain for chain in model for residue in chain}
    for generator in assemblies[0].generators:
        missing = (set(generator.chains) - chains) | (set(generator.subchains) - subchains)
        if missing:
            raise ValueError(f"assembly selectors absent from model: {sorted(missing)}")
    # - Expand chain and mmCIF subchain selections.
    for generator_index, generator in enumerate(assemblies[0].generators):
        for chain_index, chain in enumerate(model):
            selected = tuple(residue for residue in chain
                             if chain.name in generator.chains or residue.subchain in generator.subchains)
            if not selected:
                continue
            for operator_index, operator in enumerate(generator.operators):
                # - Keep source-chain parts in the same assembled copy.
                copy = f"{assembly_id}:{generator_index}:{operator_index}:{operator.name}"
                copies.append((chain, copy, operator.transform, selected))
    # - Reject empty assembly selections without fallback.
    if not copies:
        raise ValueError("assembly metadata selects no coordinate copies")
    return copies

# - Reject duplicate atoms and invalid coordinates.
def _atoms(group, transform):
    atoms, seen = [], set()
    for atom in group:
        position = transform.apply(atom.pos) if transform is not None else atom.pos
        xyz = (position.x, position.y, position.z)
        if not all(isfinite(value) for value in xyz):
            raise ValueError("non-finite atomic coordinates")
        name, altloc = atom.name.strip(), _empty(atom.altloc)
        key = (name, altloc)
        # - Keep shared atoms and relevant cofactors.
        if key in seen:
            raise ValueError(f"duplicate atom identity {key}")
        seen.add(key)
        atoms.append(Atom(name, atom.element.name, xyz, altloc))
    if not atoms:
        raise ValueError("residue contains no atoms")
    return tuple(atoms)
