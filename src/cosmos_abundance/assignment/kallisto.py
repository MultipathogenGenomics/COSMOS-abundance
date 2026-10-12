# Helper commands for Kallisto read-assignment outputs.

import argparse
import json

import polars as pl

from .common import count_reads_in_fastq_gz, get_sample_name_from_subdir


EMPTY_ABUNDANCE_SCHEMA = {
    "target_id": pl.Utf8,
    "length": pl.Int64,
    "eff_length": pl.Float64,
    "est_counts": pl.Float64,
    "tpm": pl.Float64,
}


def make_empty_abundance_file(output_path: str) -> None:
    pl.DataFrame({col: pl.Series([], dtype=dtype) for col, dtype in EMPTY_ABUNDANCE_SCHEMA.items()}).write_csv(output_path, separator="\t")


def make_empty_run_info_file(output_path: str, num_bootstraps: int, num_reads: int) -> None:
    with open(output_path, "w") as handle:
        json.dump({"n_processed": num_reads, "n_pseudoaligned": 0, "n_bootstraps": num_bootstraps}, handle, indent=4)


def make_empty_outputs(abundance_output: str, run_info_output: str, bootstraps_output_dir: str, num_bootstraps: int, read1: str) -> None:
    num_reads = count_reads_in_fastq_gz(read1)
    make_empty_abundance_file(abundance_output)
    make_empty_run_info_file(run_info_output, num_bootstraps, num_reads)
    for bootstrap in range(num_bootstraps):
        make_empty_abundance_file(f"{bootstraps_output_dir}/bs_abundance_{bootstrap}.tsv")


def combine_bootstraps(input_files: list[str], output_file: str) -> None:
    results = []
    for bootstrap, file_path in enumerate(input_files):
        df = pl.read_csv(file_path, separator="\t", schema_overrides={"target_id": pl.String, "sample": pl.String})
        df = df.with_columns(pl.col("est_counts").cast(pl.Float64), pl.lit(str(bootstrap)).alias("bootstrap"))
        results.append(df)
    pl.concat(results, how="diagonal_relaxed").write_csv(output_file, separator="\t")


def parse_kallisto_abundance_file(file_path: str) -> pl.DataFrame:
    df = pl.read_csv(file_path, separator="\t", has_header=True)
    return df.with_columns(pl.lit(get_sample_name_from_subdir(file_path)).alias("sample"))


def compute_average_across_bootstraps(df: pl.DataFrame) -> pl.DataFrame:
    rows = []
    for target in df["target_id"].unique():
        target_df = df.filter(pl.col("target_id") == target)
        row = {
            "sample": target_df.select(pl.col("sample").first()).item(),
            "target_id": target,
            "length": target_df.select(pl.col("length").first()).item(),
        }
        for col in ["eff_length", "est_counts", "tpm"]:
            target_df = target_df.with_columns(pl.col(col).cast(pl.Float64))
            row[f"{col}_mean"] = target_df.select(pl.col(col).mean()).item()
            row[f"{col}_sd"] = target_df.select(pl.col(col).std()).item()
        rows.append(row)
    if rows:
        return pl.DataFrame(rows)
    return pl.DataFrame(
        {
            "sample": pl.Series([], dtype=pl.Utf8),
            "target_id": pl.Series([], dtype=pl.Utf8),
            "length": pl.Series([], dtype=pl.Int64),
            "eff_length_mean": pl.Series([], dtype=pl.Float64),
            "eff_length_sd": pl.Series([], dtype=pl.Float64),
            "est_counts_mean": pl.Series([], dtype=pl.Float64),
            "est_counts_sd": pl.Series([], dtype=pl.Float64),
            "tpm_mean": pl.Series([], dtype=pl.Float64),
            "tpm_sd": pl.Series([], dtype=pl.Float64),
        }
    )


def collate_kallisto_results(input_files: list[str], output_file: str) -> None:
    results = []
    for file_path in input_files:
        results.append(compute_average_across_bootstraps(parse_kallisto_abundance_file(file_path)))
    pl.concat(results, how="diagonal_relaxed").write_csv(output_file, separator="\t")


def parse_run_info(file_path: str) -> tuple[int, int]:
    with open(file_path) as handle:
        data = json.load(handle)
    return data.get("n_processed", 0), data.get("n_pseudoaligned", 0)


def collate_proportion_mapped(input_files: list[str], output_file: str) -> None:
    records = []
    for file_path in input_files:
        n_processed, n_aligned = parse_run_info(file_path)
        proportion_mapped = n_aligned / n_processed if n_processed > 0 else 0
        records.append(
            {
                "sample": get_sample_name_from_subdir(file_path),
                "n_processed": n_processed,
                "n_aligned": n_aligned,
                "proportion_mapped": proportion_mapped,
                "proportion_unmapped": 1 - proportion_mapped if n_processed > 0 else 0,
            }
        )
    pl.DataFrame(records).write_csv(output_file, separator="\t")


def parse_args():
    parser = argparse.ArgumentParser(description="Kallisto read-assignment helper commands.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    empty_parser = subparsers.add_parser("make-empty")
    empty_parser.add_argument("--abundance_output", required=True)
    empty_parser.add_argument("--run_info_output", required=True)
    empty_parser.add_argument("--bootstraps_output_dir", required=True)
    empty_parser.add_argument("--num_bootstraps", type=int, required=True)
    empty_parser.add_argument("--read1", required=True)

    combine_parser = subparsers.add_parser("combine-bootstraps")
    combine_parser.add_argument("--input_files", nargs="+", required=True)
    combine_parser.add_argument("--output_file", required=True)

    collate_parser = subparsers.add_parser("collate-results")
    collate_parser.add_argument("--input_files", nargs="+", required=True)
    collate_parser.add_argument("--output_file", required=True)

    mapped_parser = subparsers.add_parser("proportion-mapped")
    mapped_parser.add_argument("--input_files", nargs="+", required=True)
    mapped_parser.add_argument("--output_file", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.command == "make-empty":
        make_empty_outputs(args.abundance_output, args.run_info_output, args.bootstraps_output_dir, args.num_bootstraps, args.read1)
    elif args.command == "combine-bootstraps":
        combine_bootstraps(args.input_files, args.output_file)
    elif args.command == "collate-results":
        collate_kallisto_results(args.input_files, args.output_file)
    elif args.command == "proportion-mapped":
        collate_proportion_mapped(args.input_files, args.output_file)


if __name__ == "__main__":
    main()
