# - Import immutable record and validation utilities.
from dataclasses import dataclass
from math import isfinite
from numbers import Real

# - Identify residues without chain or copy collisions.
@dataclass(frozen=True, order=True)
class ResidueId:
    model: int
    chain: str
    number: int
    insertion: str = ""
    copy: str = ""

# - Reject ambiguous residue identity types.
    def __post_init__(self):
        if type(self.model) is not int or self.model < 1 or type(self.number) is not int:
            raise ValueError("model must be positive; residue number integer")
        if any(not isinstance(value, str) for value in (self.chain, self.insertion, self.copy)):
            raise TypeError("chain, insertion and copy must be strings")

# - Represent prepared atoms and conformer residues.
@dataclass(frozen=True)
class Atom:
    name: str
    element: str
    position: tuple[float, float, float]
    altloc: str = ""
    # - Reject mutable positions and nonfinite atom coordinates.
    def __post_init__(self):
        if not isinstance(self.position, tuple) or len(self.position) != 3:
            raise TypeError("atom position must be an immutable coordinate triple")
        if not all(isinstance(x, Real) and not isinstance(x, bool) and isfinite(x) for x in self.position):
            raise ValueError("atom coordinates must be finite real numbers")
        if not all(isinstance(x, str) for x in (self.name, self.element, self.altloc)) or not self.name or not self.element:
            raise ValueError("atom names, elements and altloc must be strings")

# - Represent immutable prepared residue atom groups.
@dataclass(frozen=True)
class Residue:
    id: ResidueId
    name: str
    atoms: tuple[Atom, ...]
    # - Require immutable correctly typed residue atom collections.
    def __post_init__(self):
        if not isinstance(self.id, ResidueId) or not isinstance(self.atoms, tuple):
            raise TypeError("residue requires immutable identity and atom tuple")
        if not self.atoms or not all(isinstance(atom, Atom) for atom in self.atoms):
            raise ValueError("residue requires nonempty Atom records")
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("residue name must be a nonempty string")

# - Represent ordered pathways and relief markers.
@dataclass(frozen=True)
class Pathway:
    residues: tuple[ResidueId, ...]
    reliefs: tuple[float | None, ...]
    conformers: tuple[tuple[str, str], ...] = ()
    clashes: tuple[float | None, ...] = ()

# - Validate aligned finite pathway records.
# - Reject absent or mismatched residue sequences.
    def __post_init__(self):
        if not isinstance(self.residues, tuple) or not all(isinstance(item, ResidueId) for item in self.residues):
            raise TypeError("residues must be an immutable ResidueId tuple")
        if not all(isinstance(value, tuple) for value in (self.reliefs, self.conformers, self.clashes)):
            raise TypeError("pathway fields must be immutable tuples")
        n = len(self.residues)
        if not n or len(self.reliefs) != n or self.reliefs[0] is not None:
            raise ValueError("pathway residues and reliefs must align")
        # - Validate optional aligned conformer and clash fields.
        if self.conformers and len(self.conformers) != n:
            raise ValueError("conformers must align")
        if self.clashes and len(self.clashes) != n:
            raise ValueError("clashes must align")
        if any(not _measure(x) for x in self.reliefs[1:]):
            raise ValueError("reliefs must be finite and nonnegative")
# - Validate optional clash measurements.
        # - Require immutable two-string conformer transition pairs.
        if any(x is not None and not _measure(x) for x in self.clashes):
            raise ValueError("clashes must be finite and nonnegative")
        if any(not isinstance(pair, tuple) or len(pair) != 2 or not all(isinstance(x, str) for x in pair) for pair in self.conformers):
            raise ValueError("conformer pairs must contain two strings")
        if len({item.model for item in self.residues}) != 1:
            raise ValueError("a pathway cannot mix structural models")

# - Validate finite nonnegative scientific measurements.
def _measure(x):
    return isinstance(x, Real) and not isinstance(x, bool) and isfinite(x) and x >= 0

# - Define scientific settings and resource budgets.
@dataclass(frozen=True)
class Settings:
    clash_percentile: float = 0.30
    max_pathway_length: int = 100
    relief_threshold: float = 0.90
    min_network_size: int = 3
    min_pathway_length: int = 1
    residue_filter: tuple[ResidueId, ...] = ()
    side_chain_only: bool = False
# - Preserve model and assembly selection.
    model: int = 1
    cross_chain: bool = False
    assembly: str | None = None
    allow_hydrogen_free: bool = False
    max_candidates: int = 1_000_000

# - Reject invalid thresholds and integer limits.
    def __post_init__(self):
# - Validate bounded scientific thresholds.
        if not _measure(self.clash_percentile) or self.clash_percentile > 1:
            raise ValueError("clash_percentile must be in [0, 1]")
        if not _measure(self.relief_threshold):
            raise ValueError("relief_threshold must be nonnegative")
# - Require positive compute budgets.
        for value in (self.max_pathway_length, self.min_network_size, self.model, self.max_candidates):
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError("integer settings must be positive")
        if isinstance(self.min_pathway_length, bool) or not isinstance(self.min_pathway_length, int) or self.min_pathway_length < 0:
            raise ValueError("min_pathway_length must be nonnegative")

# - Preserve immutable filters and input selection.
        if not isinstance(self.residue_filter, tuple) or any(not isinstance(x, ResidueId) for x in self.residue_filter):
            raise TypeError("residue_filter must contain ResidueId values")
        if any(not isinstance(x, bool) for x in (self.side_chain_only, self.cross_chain, self.allow_hydrogen_free)):
            raise TypeError("boolean settings must be booleans")
        if self.assembly is not None and not isinstance(self.assembly, str):
            raise TypeError("assembly must be a string or None")
        if self.assembly == "":
            raise ValueError("assembly must name an explicit assembly")
