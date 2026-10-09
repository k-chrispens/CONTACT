# - Import geometry, cache, and record utilities.
from dataclasses import dataclass, field
from math import floor, isfinite
from typing import Optional
import gemmi
import numpy as np
from .models import Atom, Residue, Settings

# - Pin supported elemental van der Waals radii.
_RADII = dict(H=1.17, D=1.17, C=1.75, O=1.40, N=1.55, P=1.80, S=1.80, F=1.47, Cl=1.75, Br=1.85, I=1.98, Se=1.90, Li=1.82, Na=2.27, K=2.75, Mg=1.73, Ca=2.65, Mn=1.73, Fe=1.90, Co=1.75, Ni=1.63, Cu=1.40, Zn=1.39)

# - Identify standard protein backbone atom names.
_BACKBONE = frozenset("N CA C O OXT OT1 OT2 H H1 H2 H3 HN HN1 HN2 HN3 HA HA2 HA3".split())

# - Store immutable conformer and contact records.
@dataclass(frozen=True)
class Conformer:
    residue: Residue
    altloc: str
    atoms: tuple[Atom, ...]
    # - Cache read-only geometry for repeated contact evaluation.
    positions: np.ndarray = field(init=False, repr=False, compare=False)
    radii: np.ndarray = field(init=False, repr=False, compare=False)
    backbone: np.ndarray = field(init=False, repr=False, compare=False)
    def __post_init__(self):
        if len({atom.name for atom in self.atoms}) != len(self.atoms):
            raise ValueError("shared and alternate atoms duplicate a selected atom name")
        # - Freeze independently constructed coordinate and classification arrays.
        arrays = {"positions": np.asarray([atom.position for atom in self.atoms], dtype=float),
                  "radii": np.asarray([_radius(self.residue, atom) for atom in self.atoms]),
                  "backbone": np.asarray([_backbone(self.residue, atom) for atom in self.atoms])}
        for name, values in arrays.items():
            values.flags.writeable = False
            object.__setattr__(self, name, values)

# - Cache on-demand contact ratios symmetrically.
@dataclass
class ContactMap:
    nodes: tuple[Conformer, ...]
    cutoff: float
    neighbors: tuple[tuple[int, ...], ...]
    _ratios: dict = field(default_factory=dict, repr=False, compare=False)

# - Return measured ratio or chemical exclusion.
    def ratio(self, i: int, j: int) -> Optional[float]:
        key = (min(i, j), max(i, j))
        if self.nodes[key[0]].residue.id == self.nodes[key[1]].residue.id:
            return None
# - Honor cached chemical exclusions.
        if key not in self._ratios:
            self._ratios[key] = _ratio(self.nodes[key[0]], self.nodes[key[1]])
# - Return absent chemical measurements explicitly.
        return self._ratios.get(key)

# - Classify backbone only in amino-acid residues.
def _backbone(residue, atom):
    amino = gemmi.find_tabulated_residue(residue.name).is_amino_acid()
    name = atom.name.upper()
    if atom.element.upper() == "D" and name.startswith("D"):
        name = "H" + name[1:]
    return amino and name in _BACKBONE

# - Resolve radii with specific atom corrections.
def _radius(residue, atom):
    element, name, res = atom.element.upper().capitalize(), atom.name.upper(), residue.name.upper()
    if element in ("H", "D") and name.startswith("D"):
        name = "H" + name[1:]
# - Reject elements outside the pinned radii table.
    radius = _RADII.get(element)
    if radius is None:
        raise ValueError(f"unsupported element: {atom.element}")
    if element == "C" and name == "C" and gemmi.find_tabulated_residue(residue.name).is_amino_acid():
        return 1.65
# - Adjust hydrogen radii by atom identity.
    if element in ("H", "D"):
        if (_backbone(residue, atom) and not name.startswith("HA")) or name in _NAMED_H.get(res, ()):
            return 1.00
    return radius

# - Map chemically shortened side-chain hydrogens.
_NAMED_H = {
    "PHE": "HD1 HE1 HZ HE2 HD2", "TRP": "HD1 HE1 HZ2 HH2 HZ3 HE3",
    "TYR": "HD1 HE1 HE2 HD2", "ASN": "HD21 HD22", "GLN": "HE21 HE22",
    "LYS": "HZ1 HZ2 HZ3", "ARG": "HE HH11 HH12 HH21 HH22",
    "SER": "HG", "CYS": "HG", "THR": "HG1",
}
_NAMED_H = {key: frozenset(value.split()) for key, value in _NAMED_H.items()}

# - Measure minimum eligible atom-pair normalized distance.
def _ratio(left, right):
    if left.residue.name.upper() == right.residue.name.upper() == "CYS":
        return None
    mask = np.ones((len(left.atoms), len(right.atoms)), dtype=bool)
# - Ignore chemically masked atom pairs.
    if left.backbone.any() and right.backbone.any() and _backbone_excluded(left, right):
        mask &= ~(left.backbone[:, None] & right.backbone[None, :])
    if not mask.any():
        return None
    # - Vectorize exact eligible atom-pair distance evaluation.
    differences = left.positions[:, None, :] - right.positions[None, :, :]
    distances = np.linalg.norm(differences, axis=2)
    ratios = distances / (left.radii[:, None] + right.radii[None, :])
    value = float(np.min(ratios[mask]))
    if not isfinite(value):
        raise ValueError("normalized contact geometry is not representable")
    return value

# - Exclude published cysteine and backbone interactions.
def _backbone_excluded(left, right):
    li, ri = left.residue.id, right.residue.id
    adjacent = (li.model, li.chain, li.copy) == (ri.model, ri.chain, ri.copy) and abs(li.number-ri.number) <= 1
    return adjacent or _hbond(left, right) or _hbond(right, left)

# - Infer the MMDB directional mainchain hydrogen bond.
def _hbond(donor, acceptor):
    for n in donor.atoms:
        if n.name.upper() != "N":
            continue
        for o in acceptor.atoms:
            if o.name.upper() != "O":
                continue
            no = tuple(x-y for x, y in zip(n.position, o.position))
            if sum(x*x for x in no) > 12.25 or not any(no):
                continue
# - Check acceptor carbon angle without inverse trig.
            for c in acceptor.atoms:
                if c.name.upper() != "C":
                    continue
                co = tuple(x-y for x, y in zip(c.position, o.position))
                if any(co) and sum(x*y for x, y in zip(co, no)) <= 0:
                    return True
    return False

# - Build selected conformers and exact spatial candidates.
def build_contacts(residues: tuple[Residue, ...], settings: Settings) -> ContactMap:
    nodes = []
    for residue in sorted(residues, key=lambda r: r.id):
        atoms = tuple(a for a in residue.atoms if not (settings.side_chain_only and _backbone(residue, a)))
        alts = sorted({a.altloc for a in atoms if a.altloc.strip()})
        if len(alts) < 2:
            continue
        for alt in alts:
            selected = tuple(a for a in atoms if not a.altloc.strip() or a.altloc == alt)
            nodes.append(Conformer(residue, alt, selected))
# - Handle structures without mobile conformers.
    if not nodes:
        return ContactMap((), 0.0, ())
# - Build spatial cell occupancy by selected atoms.
    width = 2 * max(float(np.max(node.radii)) for node in nodes) * 1.0001
    cells = {}
    for i, node in enumerate(nodes):
        for atom in node.atoms:
            cell = tuple(floor(x / width) for x in atom.position)
            cells.setdefault(cell, set()).add(i)
# - Visit neighboring cells once per occupancy cell.
    candidates = set()
    for cell, indices in cells.items():
        nearby = set(indices)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    nearby.update(cells.get((cell[0]+dx, cell[1]+dy, cell[2]+dz), ()))
        candidates.update((min(i,j), max(i,j)) for i in indices for j in nearby if i != j and nodes[i].residue.id != nodes[j].residue.id)
# - Measure candidates and select percentile clashes.
    ratios = {(i, j): _ratio(nodes[i], nodes[j]) for i, j in sorted(candidates)}
# - Select requested rank from sorted clash sample.
    sample = sorted(value for value in ratios.values() if value is not None and value < 1.0001)
    cutoff = sample[min(floor(settings.clash_percentile * len(sample)), len(sample)-1)] if sample else 0.0
# - Keep all ties and sort adjacency deterministically.
    neighbors = [set() for _ in nodes]
    for (i,j), value in ratios.items():
        if value is not None and value <= cutoff:
            neighbors[i].add(j)
            neighbors[j].add(i)
    return ContactMap(tuple(nodes), cutoff, tuple(tuple(sorted(n)) for n in neighbors), ratios)
