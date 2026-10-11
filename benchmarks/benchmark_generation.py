# - Build exact prepared multi-conformer cascade workloads.
import argparse
import json
import platform
import resource
import sys
import time
import tracemalloc
from contact_networks import Atom, Residue, ResidueId, Settings
from contact_networks.engine import generate_pathways

# - Construct distinct chains with exact independent transition cascades.
def workload(copies):
    residues = []
    for index in range(copies):
        for number, positions in enumerate(((0, -10), (2, 10), (12, 30)), 1):
            atoms = tuple(Atom("CB", "C", (x + 100 * index, 0, 0), alt) for x, alt in zip(positions, "AB"))
            residues.append(Residue(ResidueId(1, f"chain{index:06d}", number), "LIG", atoms))
    return tuple(residues)

# - Measure corrected generation, peak memory, and state counts.
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--copies", type=int, default=100)
    args = parser.parse_args()
    if args.copies < 1:
        parser.error("copies must be positive")
    residues = workload(args.copies)
    settings = Settings(allow_hydrogen_free=True, max_candidates=6 * args.copies)
    # - Profile the exact prepared generation workload only.
    tracemalloc.start()
    started = time.perf_counter()
    result = generate_pathways(residues, settings)
    elapsed = time.perf_counter() - started
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

# - Validate expected completed pathways before reporting timings.
    assert result.status == "complete"
    assert len(result.pathways) == 4 * args.copies
    assert result.candidate_states == 6 * args.copies
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    print(json.dumps({"copies": args.copies, "residues": len(residues), "paths": len(result.pathways),
                      "expanded_states": result.candidate_states, "seconds": elapsed, "peak_python_bytes": peak,
                      "peak_process_rss_bytes": rss, "python": platform.python_version()}, indent=2))

# - Execute measurements only when explicitly invoked.
if __name__ == "__main__":
    main()
