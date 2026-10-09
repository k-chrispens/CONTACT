# - Import deterministic serialization and filesystem utilities.
import csv
import json
from collections import Counter
from dataclasses import asdict
from pathlib import Path

import networkx as nx

from .networks import maximal_pathways, network_pathway_counts


# - Sort full residue identities without delimiter collisions.
def _residue_key(residue):
    return (residue.model, residue.chain, residue.number, residue.insertion, residue.copy)


# - Serialize model, chain, insertion, and copy identity.
def _residue_data(residue):
    return {name: getattr(residue, name) for name in ("model", "chain", "number", "insertion", "copy")}


# - Preserve CSV identities as unambiguous JSON objects.
def _residue_id(residue):
    return json.dumps(_residue_data(residue), sort_keys=True, separators=(",", ":"))


# - Export aligned pathway measurements without losing zero.
def _pathway_data(pathway):
    return {"residues": [_residue_data(r) for r in pathway.residues],
            "reliefs": list(pathway.reliefs),
            "conformers": [list(pair) for pair in pathway.conformers],
            "clashes": list(pathway.clashes) if pathway.clashes is not None else None}


# - Assign stable raw and accepted pathway identifiers.
def _pathway_rows(pathways, name):
    return [{"path_id": f"{name}_{i:06d}", **_pathway_data(p)} for i, p in enumerate(pathways, 1)]


# - Assign stable residue and connected-network identifiers.
def to_document(result):
    occurrences = Counter(node for path in result.raw_pathways for node in path.residues)
    nodes = sorted(set(occurrences) | set(result.graph), key=_residue_key)
    node_ids = {node: i for i, node in enumerate(nodes, 1)}
    components = sorted((sorted(c, key=_residue_key) for c in nx.connected_components(result.graph)),
                       key=lambda c: tuple(_residue_key(n) for n in c))
    membership = {node: f"network_{i:04d}" for i, part in enumerate(components, 1) for node in part}
    raw_maxima = maximal_pathways(result.raw_pathways)
    settings = asdict(result.settings)
    settings["residue_filter"] = list(settings["residue_filter"])
    # - Serialize parameters, provenance, and completion status.
    # - Scope imported truncation to the source input.
    return {
        "source": result.source, "settings": settings,
        "input_kind": result.input_kind, "generation_settings_known": False,
        "version": _version(), "status": result.status,
        "stopping_reason": result.stopping_reason, "thresholds": result.thresholds,
        "processed_groups": list(result.processed_groups),
        "counts": {"raw_pathways": len(result.raw_pathways), "accepted_pathways": len(result.pathways),
                   "maximal_raw_pathways": len(raw_maxima),
                   "raw_residue_occurrences": sum(len(p.residues) for p in result.raw_pathways)},
        # - Export raw and accepted aligned pathway records.
        "raw_pathways": _pathway_rows(result.raw_pathways, "raw"),
        "pathways": _pathway_rows(result.pathways, "accepted"),
        "nodes": [{"node_id": node_ids[n], **_residue_data(n),
                   "raw_occurrences": occurrences[n], "network_id": membership.get(n)}
                  for n in nodes],
        "edges": _edge_rows(result.graph, node_ids),
        "networks": _network_rows(result, components),
    }


# - Report the installed package version.
def _version():
    from importlib.metadata import PackageNotFoundError, version
    try:
        return version("contact-networks")
    except PackageNotFoundError:
        return "unknown"


# - Count each accepted pathway once per network.
def _network_rows(result, components):
    counts = network_pathway_counts(result.graph, result.pathways)
    return [{"network_id": f"network_{i:04d}", "size": len(c),
             "accepted_pathway_count": counts[i - 1]}
            for i, c in enumerate(components, 1)]


# - Export deterministically oriented weighted undirected edges.
def _edge_rows(graph, node_ids):
    edges = []
    for u, v, data in graph.edges(data=True):
        source, target = sorted((node_ids[u], node_ids[v]))
        edges.append({"source": source, "target": target, "weight": data["weight"]})
    return sorted(edges, key=lambda edge: (edge["source"], edge["target"]))


# - Write each CSV with an explicit column schema.
def _write_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


# - Export residue occurrences and weighted network edges.
def _aligned_rows(pathways, label):
    rows = []
    for index, p in enumerate(pathways, 1):
        for step, residue in enumerate(p.residues, 1):
            rows.append({"path_id": f"{label}_{index:06d}", "step": step,
                         "residue_id": _residue_id(residue),
                         "conformers": json.dumps(p.conformers[step - 1])
                         if p.conformers and p.conformers[step - 1] is not None else "",
                         "clash": _blank(p.clashes, step), "relief": _blank(p.reliefs, step)})
    return rows


# - Distinguish missing numeric values from genuine zero.
def _blank(values, step):
    value = values[step - 1] if values else None
    return "" if value is None else value


# - Render optional headless uniquely numbered network PNGs.
def _plot(graph, components, directory, Figure, FigureCanvasAgg):
    for index, component in enumerate(components, 1):
        subgraph = nx.Graph()
        subgraph.add_nodes_from(component)
        edges = [(u, v, data) for u, v, data in graph.subgraph(component).edges(data=True)]
        subgraph.add_edges_from(sorted(edges, key=lambda edge: (_residue_key(edge[0]), _residue_key(edge[1]))))
        positions = nx.spring_layout(subgraph, seed=0)
        # - Label full identities instead of ambiguous residue numbers.
        figure = Figure()
        FigureCanvasAgg(figure)
        labels = {node: f"{node.chain}:{node.number}{node.insertion}:{node.copy}" for node in component}
        nx.draw_networkx(subgraph, positions, ax=figure.subplots(), labels=labels, font_size=7)
        figure.savefig(directory / f"network_{index:04d}.png")


# - Preflight optional plotting before reserving output destinations.
def write_results(result, destination, plots=False):
    destination = Path(destination)
    if plots:
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
    else:
        Figure = FigureCanvasAgg = None
    # - Refuse existing destinations and shared output names.
    document = to_document(result)
    payload = json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n"
    destination.mkdir(parents=True, exist_ok=False)
    raw, accepted = result.raw_pathways, result.pathways
    _write_csv(destination / "raw_pathways.csv", ("path_id", "step", "residue_id", "conformers", "clash", "relief"),
               _aligned_rows(raw, "raw"))
    _write_csv(destination / "pathways.csv", ("path_id", "step", "residue_id", "conformers", "clash", "relief"),
               _aligned_rows(accepted, "accepted"))
    _write_csv(destination / "residues.csv", ("node_id", "model", "chain", "number", "insertion", "copy", "raw_occurrences", "network_id"),
               document["nodes"])
    # - Export network memberships and unique numbered PNGs.
    _write_csv(destination / "edges.csv", ("source", "target", "weight"), document["edges"])
    _write_csv(destination / "networks.csv", ("network_id", "size", "accepted_pathway_count"), document["networks"])
    if plots:
        components = sorted((sorted(c, key=_residue_key) for c in nx.connected_components(result.graph)),
                            key=lambda c: tuple(_residue_key(n) for n in c))
        _plot(result.graph, components, destination, Figure, FigureCanvasAgg)
    # - Write completion manifest only after successful exports.
    manifest = destination / ".results.json.tmp"
    manifest.write_text(payload, encoding="utf-8")
    manifest.replace(destination / "results.json")
    return destination
