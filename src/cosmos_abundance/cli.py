import argparse
import math
import subprocess

from . import __version__


def positive_int(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return value


def nonnegative_float(value):
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise argparse.ArgumentTypeError("must be finite and non-negative")
    return value


def percentage(value):
    value = nonnegative_float(value)
    if value > 100:
        raise argparse.ArgumentTypeError("must be between 0 and 100")
    return value


def proportion(value):
    value = nonnegative_float(value)
    if value > 1:
        raise argparse.ArgumentTypeError("must be between 0 and 1")
    return value


def occupancy(value):
    value = nonnegative_float(value)
    if not 0 < value <= 1:
        raise argparse.ArgumentTypeError("must be greater than 0 and at most 1")
    return value


def add_filter_args(parser):
    parser.add_argument("--min-coverage", type=percentage, default=90, help="Minimum KMA template coverage (%%; default: 90)")
    parser.add_argument("--min-identity", type=percentage, default=97.5, help="Minimum KMA identity across covered length (%%; default: 97.5)")
    parser.add_argument("--min-kma-count", type=nonnegative_float, default=5, help="Minimum KMA read count (default: 5)")
    parser.add_argument("--min-kallisto-count", type=nonnegative_float, default=5, help="Minimum kallisto bootstrap mean count (default: 5)")
    parser.add_argument("--allow-single-tool", action="store_true", help="Disable the both-tools detection requirement; other filters still apply")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Prepare reference sets and estimate read-based relative abundance.")
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare", help="Align, trim and cluster a reference FASTA; draw the Mash-distance tree")
    prepare.add_argument("--fasta", required=True)
    prepare.add_argument("--species", help="Tab-separated id/species file; optional tree_label column with a header")
    prepare.add_argument("--out", required=True)
    prepare.add_argument("--mash-distance-threshold", type=nonnegative_float, default=0.01, help="Complete-linkage clustering threshold (default: 0.01)")
    prepare.add_argument("--min-reference-length", type=positive_int, default=500, help="Minimum ungapped length before and after trimming (default: 500 bp)")
    prepare.add_argument("--max-n-fraction", type=proportion, default=0.05, help="Maximum fraction of N bases before alignment (0 to 1; default: 0.05)")
    prepare.add_argument("--kmer-size", type=positive_int, default=31, help="K-mer size for Mash-like distances (default: 31)")
    prepare.add_argument("--allowed-gap-size", type=int, default=3, help="Maximum intervening gappy columns in the shared region (default: 3)")
    prepare.add_argument("--min-occupancy", type=occupancy, default=0.5, help="Minimum occupied fraction at alignment ends (default: 0.5)")
    prepare.add_argument("--threads", type=positive_int, default=1)
    run = commands.add_parser("run", help="Run KMA, kallisto and abundance estimation")
    references = run.add_mutually_exclusive_group(required=True)
    references.add_argument("--reference", help="Prepared reference directory")
    references.add_argument("--fasta", help="Use an existing FASTA without alignment, trimming or clustering")
    run.add_argument("--species", help="ID/species mapping for --fasta")
    run.add_argument("--samples", required=True, help="TSV with sample, read1 and read2 columns")
    run.add_argument("--out", required=True)
    run.add_argument("--cores", type=positive_int, default=1, help="Total workflow CPU budget")
    run.add_argument("--threads", type=positive_int, default=1, help="Maximum threads per assignment job")
    run.add_argument("--bootstraps", type=positive_int, default=100)
    run.add_argument("--seed", type=positive_int, default=42)
    run.add_argument("--dry-run", action="store_true", help="Write run configuration and show jobs without executing them")
    add_filter_args(run)
    estimate = commands.add_parser("estimate", help="Recompute abundance from a completed run")
    estimate.add_argument("--run-dir", required=True)
    add_filter_args(estimate)
    for subparser in [prepare, run, estimate]:
        subparser.add_argument("--keep-intermediate", action="store_true", help="Retain intermediate files after successful completion")
    args = parser.parse_args(argv)
    if args.command == "prepare":
        if args.allowed_gap_size < 0:
            parser.error("--allowed-gap-size must be non-negative")
        if args.kmer_size > args.min_reference_length:
            parser.error("--kmer-size cannot exceed --min-reference-length")
    return args


def main(argv=None):
    args = parse_args(argv)
    try:
        if args.command == "prepare":
            from .reference.prepare import prepare_reference
            prepare_reference(args)
        elif args.command == "run":
            from .run import run_pipeline
            run_pipeline(args)
        else:
            from .abundance import estimate_abundance
            estimate_abundance(args.run_dir, args.min_coverage, args.min_identity,
                               args.min_kma_count, args.min_kallisto_count, not args.allow_single_tool, args.keep_intermediate)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"cosmos-abundance: {error}") from error
