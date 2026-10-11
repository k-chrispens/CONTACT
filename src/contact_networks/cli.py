# - Import command parsing and shared analysis operations.
import argparse
import json
import sys
from .analysis import analyze, analyze_pathways
from .models import ResidueId, Settings
from .networks import read_pathways
from .output import write_results

# - Expose scientific settings and input selection.
# - Offer corrected structure analysis alongside pathway import.
def _parser():
    parser = argparse.ArgumentParser(description="CONTACT analysis")
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("analyze").add_argument("input")
    commands.add_parser("extract").add_argument("input")
    # - Share common filters and output controls.
    for name in ("analyze", "extract"):
        command = commands.choices[name]
        command.add_argument("--model", type=int, default=1)
        command.add_argument("--relief-threshold", type=float, default=0.9)
        command.add_argument("--min-network-size", type=int, default=3)
        command.add_argument("--min-pathway-length", type=int, default=1)
        command.add_argument("--residue", action="append", default=[])
        command.add_argument("--output", required=True)
        command.add_argument("--plots", action="store_true")
    # - Keep generation-only options off pathway imports.
    generated = commands.choices["analyze"]
    generated.add_argument("--clash-percentile", type=float, default=0.30)
    generated.add_argument("--max-pathway-length", type=int, default=100)
    generated.add_argument("--max-candidates", type=int, default=1_000_000)
    generated.add_argument("--side-chain-only", action="store_true")
    generated.add_argument("--cross-chain", action="store_true")
    generated.add_argument("--assembly")
    generated.add_argument("--allow-hydrogen-free", action="store_true")
    # - Preserve upstream completion metadata for imports.
    imported = commands.choices["extract"]
    imported.add_argument("--input-status", choices=("complete", "truncated"), default="complete")
    imported.add_argument("--stopping-reason")
    return parser

# - Dispatch structure or detailed-pathway analysis.
# - Forward scientific scope and shared generation budgets.
def main(argv=None):
    args = _parser().parse_args(argv)
    try:
        filters = tuple(ResidueId(**json.loads(value)) for value in args.residue)
        # - Keep each input mode's scientific settings explicit.
        if args.command == "analyze":
            settings = Settings(model=args.model, relief_threshold=args.relief_threshold,
                                min_network_size=args.min_network_size, min_pathway_length=args.min_pathway_length,
                                residue_filter=filters, clash_percentile=args.clash_percentile,
                                max_pathway_length=args.max_pathway_length, max_candidates=args.max_candidates,
                                side_chain_only=args.side_chain_only, cross_chain=args.cross_chain,
                                assembly=args.assembly, allow_hydrogen_free=args.allow_hydrogen_free)
            result = analyze(args.input, settings)
        else:
            # - Preserve imported model filtering and status.
            settings = Settings(model=args.model, relief_threshold=args.relief_threshold,
                                min_network_size=args.min_network_size, min_pathway_length=args.min_pathway_length,
                                residue_filter=filters)
            paths = (path for path in read_pathways(args.input) if path.residues[0].model == args.model)
            result = analyze_pathways(paths, settings, source=str(args.input),
                                      status=args.input_status, stopping_reason=args.stopping_reason)
        # - Export through the shared deterministic writer.
        # - Finalize isolated exports before reporting completion.
        write_results(result, args.output, plots=args.plots)

# - Report ordinary errors without false successful status.
    except (OSError, ValueError, TypeError, ImportError) as error:
        print(f"contact-networks: {error}", file=sys.stderr)
        return 1

# - Preserve partial outputs and return truncation status.
    return 3 if result.status == "truncated" else 0
