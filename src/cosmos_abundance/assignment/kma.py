# Helper commands for KMA read-assignment outputs.

import argparse

import polars as pl

from .common import get_sample_name_from_subdir


def parse_kma_res_file(file_path: str) -> pl.DataFrame:
    df = pl.read_csv(file_path, separator="\t", has_header=True, schema_overrides={"#Template": pl.String})
    return df.with_columns(pl.lit(get_sample_name_from_subdir(file_path)).alias("sample"))


def parse_kma_mapstat_file(file_path: str) -> pl.DataFrame:
    df = pl.read_csv(file_path, separator="\t", has_header=True, comment_prefix="##", schema_overrides={"# refSequence": pl.String})
    return df.with_columns(pl.lit(get_sample_name_from_subdir(file_path)).alias("sample"))


def collate_kma_results(input_files: list[str], output_file: str) -> None:
    results = []
    for file_path in input_files:
        df = parse_kma_res_file(file_path)
        mapstat_df = parse_kma_mapstat_file(file_path.replace(".res", ".mapstat"))
        df = df.join(mapstat_df, left_on=["sample", "#Template"], right_on=["sample", "# refSequence"], how="left")
        for col in df.columns:
            if df[col].dtype == pl.String:
                df = df.with_columns(pl.col(col).str.strip_chars().alias(col))
        results.append(df)
    pl.concat(results, how="diagonal_relaxed").write_csv(output_file, separator="\t")


def extract_processed_from_mapstat(file_path: str) -> int | None:
    with open(file_path) as handle:
        for line in handle:
            if line.startswith("## fragmentCount"):
                return int(line.strip().split("\t")[1].strip())
    return None


def compute_n_aligned_from_mapstat(file_path: str) -> float:
    frag_file = file_path.replace(".mapstat", ".frag.gz")
    try:
        df = pl.read_csv(frag_file, separator="\t", has_header=False)
    except pl.exceptions.NoDataError:
        return 0
    df.columns = ["query_seq", "n_mapping_templates", "alignment_score", "start_pos", "end_pos", "template_id", "query_id"]
    df = df.select(pl.col("query_id").str.split(" ").list.get(0).alias("query_id"))
    n_unique_entries = df.select(pl.col("query_id").unique().count()).item()
    if all(df["query_id"].str.ends_with("/1") | df["query_id"].str.ends_with("/2")):
        return n_unique_entries / 2
    return n_unique_entries


def collate_proportion_mapped(input_files: list[str], output_file: str) -> None:
    records = []
    for file_path in input_files:
        n_processed = extract_processed_from_mapstat(file_path)
        n_aligned = compute_n_aligned_from_mapstat(file_path)
        proportion_mapped = n_aligned / n_processed if n_processed and n_processed > 0 else 0
        records.append(
            {
                "sample": get_sample_name_from_subdir(file_path),
                "n_processed": n_processed,
                "n_aligned": n_aligned,
                "proportion_mapped": proportion_mapped,
                "proportion_unmapped": 1 - proportion_mapped if n_processed and n_processed > 0 else 0,
            }
        )
    pl.DataFrame(records).write_csv(output_file, separator="\t")


def parse_args():
    parser = argparse.ArgumentParser(description="KMA read-assignment helper commands.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    collate_parser = subparsers.add_parser("collate-results")
    collate_parser.add_argument("--input_files", nargs="+", required=True)
    collate_parser.add_argument("--output_file", required=True)

    mapped_parser = subparsers.add_parser("proportion-mapped")
    mapped_parser.add_argument("--input_files", nargs="+", required=True)
    mapped_parser.add_argument("--output_file", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.command == "collate-results":
        collate_kma_results(args.input_files, args.output_file)
    elif args.command == "proportion-mapped":
        collate_proportion_mapped(args.input_files, args.output_file)


if __name__ == "__main__":
    main()
