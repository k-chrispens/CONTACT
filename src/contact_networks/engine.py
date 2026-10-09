# - Import corrected contact maps and pathway records.
import json
from collections import defaultdict
from dataclasses import dataclass
from .contacts import build_contacts
from .models import Pathway, Residue, Settings

# - Record bounded generation and completed group provenance.
@dataclass(frozen=True)
class Generation:
    pathways: tuple[Pathway, ...]
    status: str
    stopping_reason: str | None
    thresholds: dict[str, float]
    processed_groups: tuple[str, ...]
    candidate_states: int

# - Group selected residues by chain and assembly copy.
def _groups(residues, settings):
    groups = defaultdict(list)
    for residue in sorted(residues, key=lambda item: item.id):
        key = ("cross_chain", residue.id.model) if settings.cross_chain else (
            "chain", residue.id.model, residue.id.chain, residue.id.copy)
        groups[json.dumps(key, separators=(",", ":"))].append(residue)
    return tuple((label, tuple(group)) for label, group in sorted(groups.items()))

# - Build alternate-state lookup without mixing residue identities.
def _alternatives(contacts):
    lookup = defaultdict(list)
    for index, node in enumerate(contacts.nodes):
        lookup[node.residue.id].append(index)
    return tuple(tuple(lookup[node.residue.id]) for node in contacts.nodes)

# - Seed clashes introduced by an alternate-state change.
def _introduced(contacts, before, after, target):
    old, new = contacts.ratio(before, target), contacts.ratio(after, target)
    return new is not None and (old is None or old >= new)

def _seeds(contacts, alternatives):
    for head, neighbors in enumerate(contacts.neighbors):
        for tail in neighbors:
            if any(old != head and _introduced(contacts, old, head, tail) for old in alternatives[head]):
                yield ((head, tail),), ()

# - Carry branch-local scores through iterative depth-first search.
def _states(contacts, alternatives, maximum, remaining):
    count = 0
    for seed in _seeds(contacts, alternatives):
        stack = [seed]
        while stack:
            # - Stop before expanding any over-budget state.
            if count == remaining:
                yield None
                return
            count += 1
            state = stack.pop()
            terminals, children = _expand(contacts, alternatives, state, maximum)
            yield terminals
            stack.extend(reversed(children))

# - Require every switch to relieve previous path clashes.
def _relieves(contacts, edges, alternate):
    for head, _ in edges:
        score = contacts.ratio(head, alternate)
        if score is not None and score <= contacts.cutoff:
            return False
    return True

# - Propagate only newly introduced or worsened clashes.
def _expand(contacts, alternatives, state, maximum):
    edges, reliefs = state
    head, tail = edges[-1]
    visited = {contacts.nodes[node].residue.id for edge in edges for node in edge}
    terminals, children = [], []
    for alternate in alternatives[tail]:
        if alternate == tail or not _relieves(contacts, edges, alternate):
            continue
        score = contacts.ratio(head, alternate)
        # - Follow only admissible forward conformer transitions.
        targets = tuple(target for target in contacts.neighbors[alternate]
                        if contacts.nodes[target].residue.id not in visited
                        and _introduced(contacts, tail, alternate, target))
        if not targets:
            terminals.append(_pathway(contacts, edges, reliefs, alternate, score))
        elif len(edges) < maximum:
            children.extend((edges + ((alternate, target),), reliefs + (score,)) for target in targets)
    return terminals, children

# - Emit only terminated chemically valid pathways.
def _pathway(contacts, edges, reliefs, alternate, terminal_score):
    nodes = contacts.nodes
    identities = tuple(nodes[head].residue.id for head, _ in edges) + (nodes[edges[-1][1]].residue.id,)
    transitions = ((nodes[edges[0][0]].altloc, "-"),) + tuple(
        (nodes[edges[index - 1][1]].altloc, nodes[edges[index][0]].altloc) for index in range(1, len(edges)))
    transitions += ((nodes[edges[-1][1]].altloc, nodes[alternate].altloc),)
    scores = (None,) + reliefs + (terminal_score,)
    clashes = (None,) + tuple(contacts.ratio(head, tail) for head, tail in edges)
    excluded = tuple(index > 0 and score is None for index, score in enumerate(scores))
    return Pathway(identities, scores, transitions, clashes, excluded)

# - Bound transitions without mislabeling unfinished terminal branches.
def generate_pathways(residues, settings=None):
    settings = Settings() if settings is None else settings
    residues = tuple(residues)
    if any(not isinstance(residue, Residue) or residue.id.model != settings.model for residue in residues):
        raise ValueError("prepared residues must belong to the selected model")
    if len({residue.id for residue in residues}) != len(residues):
        raise ValueError("prepared residue identities must be unique")
    # - Require explicit permission for hydrogen-free prepared records.
    if residues and not settings.allow_hydrogen_free and not any(
        atom.element in ("H", "D") for residue in residues for atom in residue.atoms):
        raise ValueError("prepared residues contain no hydrogen atoms")
    paths, thresholds, processed, count = [], {}, [], 0

# - Enforce candidate-state budgets across all groups.
    for label, group in _groups(residues, settings):
        contacts = build_contacts(group, settings)
        alternatives = _alternatives(contacts)
        thresholds[label] = contacts.cutoff
        # - Share one state budget across independently scored groups.
        for terminals in _states(contacts, alternatives, settings.max_pathway_length, settings.max_candidates - count):
            if terminals is None:
                reason = f"max_candidates={settings.max_candidates} reached in group {label}"
                return Generation(tuple(paths), "truncated", reason, thresholds, tuple(processed), count)
            count += 1
            paths.extend(terminals)
        processed.append(label)

# - Return partial results with explicit stopping reasons.
    return Generation(tuple(paths), "complete", None, thresholds, tuple(processed), count)
