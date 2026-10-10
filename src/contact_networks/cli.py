# - Import command parsing and shared analysis operations.
import argparse
import json
import sys
from .analysis import analyze_pathways
from .models import ResidueId, Settings
from .networks import read_pathways
from .output import write_results

# - Expose scientific settings and input selection.
# - Offer corrected structure analysis alongside pathway import.
def _parser():
    parser = argparse.ArgumentParser(description="CONTACT analysis; pathway generation is a separate stage")
    command = parser.add_subparsers(dest="command", required=True).add_parser("extract")
    command.add_argument("input", help="original detailed CONTACT pathway file")
    command.add_argument("--model", type=int, default=1)
    command.add_argument("--relief-threshold", type=float, default=0.9)
    command.add_argument("--min-network-size", type=int, default=3)
    command.add_argument("--min-pathway-length", type=int, default=1)
    command.add_argument("--residue", action="append", default=[], help="JSON ResidueId object; repeatable")

# - Require explicit per-run output destinations.
    command.add_argument("--output", required=True)
    command.add_argument("--plots", action="store_true")
    command.add_argument("--input-status", choices=("complete", "truncated"), default="complete")
    command.add_argument("--stopping-reason")
    return parser

# - Dispatch structure or detailed-pathway analysis.
# - Forward scientific scope and shared generation budgets.
def main(argv=None):
    args = _parser().parse_args(argv)
    try:
        filters = tuple(ResidueId(**json.loads(value)) for value in args.residue)
        settings = Settings(model=args.model, relief_threshold=args.relief_threshold,
                            min_network_size=args.min_network_size, min_pathway_length=args.min_pathway_length,
                            residue_filter=filters)
        paths = (path for path in read_pathways(args.input) if path.residues[0].model == args.model)
        result = analyze_pathways(paths, settings, source=str(args.input),
                                  status=args.input_status, stopping_reason=args.stopping_reason)
        # - Finalize isolated exports before reporting completion.
        write_results(result, args.output, plots=args.plots)

# - Report ordinary errors without false successful status.
    except (OSError, ValueError, TypeError, ImportError) as error:
        print(f"contact-networks: {error}", file=sys.stderr)
        return 1

# - Preserve partial outputs and return truncation status.
    return 3 if result.status == "truncated" else 0
