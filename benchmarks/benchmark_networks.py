# - Construct deterministic repeated weighted pathway workloads.
import argparse
import json
import platform
import time
import tracemalloc
import networkx as nx
from contact_networks import Pathway, ResidueId, build_graph

# - Compare legacy scans against direct graph aggregation.
def legacy_graph(paths):
    graph = nx.Graph()
    for path in paths:
        for left, right in zip(path.residues, path.residues[1:]):
            found = False
            for edge in graph.edges:
                if left in edge and right in edge:
                    graph[edge[0]][edge[1]]["weight"] += 1
                    found = True
            # - Match legacy insertion without changing edge counts.
            if not found:
                graph.add_edge(left, right, weight=1)
    return graph

# - Assert exact weighted graph equality before timing.
def measure(function, paths):
    tracemalloc.start()
    started = time.perf_counter()
    graph = function(paths)
    seconds = time.perf_counter() - started
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return graph, {"seconds": seconds, "peak_python_bytes": peak}

# - Report workload, runtime, memory, and software versions.
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paths", type=int, default=2000)
    parser.add_argument("--residues", type=int, default=100)
    args = parser.parse_args()
    if args.paths < 1 or args.residues < 4:
        parser.error("require positive paths and at least four residues")
    residues = tuple(ResidueId(1, "A", number) for number in range(args.residues))
    # - Generate reproducible distinct four-residue pathway transitions.
    paths = tuple(Pathway(tuple(residues[(index + step) % args.residues] for step in range(4)),
                          (None, 1.0, 1.0, 1.0)) for index in range(args.paths))
    slow, baseline = measure(legacy_graph, paths)
    fast, optimized = measure(build_graph, paths)
    assert nx.utils.graphs_equal(slow, fast), "scientific graph mismatch"
    print(json.dumps({"paths": args.paths, "residues": args.residues, "python": platform.python_version(),
                      "networkx": nx.__version__, "legacy": baseline, "optimized": optimized,
                      "speedup": baseline["seconds"] / optimized["seconds"]}, indent=2))

# - Run benchmarks only when explicitly invoked.
if __name__ == "__main__":
    main()
