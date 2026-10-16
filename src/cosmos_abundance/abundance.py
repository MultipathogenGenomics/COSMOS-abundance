import json
from pathlib import Path

import polars as pl

from .io import write_json, clean_run_outputs


def _rename_with_prefix(df: pl.DataFrame, prefix: str, keep_cols: list[str]) -> pl.DataFrame:
    return df.rename({col: f"{prefix}:{col}" for col in df.columns if col not in keep_cols})


def _rename_mapped_to_assigned(df: pl.DataFrame) -> pl.DataFrame:
    rename_dict = {
        col: col.replace("mapped", "assigned").replace("aligned", "assigned")
        for col in df.columns
        if col != "sample"
    }
    return df.rename(rename_dict)


def join_kallisto_kma_outputs(read_assignment_dir: str) -> tuple[pl.DataFrame, pl.DataFrame]:
    root = Path(read_assignment_dir) / "reference" / "samples"
    paths = {
        "kallisto_results": root / "kallisto/results/data/kallisto_results.all_samples.tsv",
        "kma_results": root / "kma/results/data/kma_results.all_samples.tsv",
        "kallisto_assigned": root / "kallisto/results/data/kallisto_proportion_mapped.all_samples.tsv",
        "kma_assigned": root / "kma/results/data/kma_proportion_mapped.all_samples.tsv",
    }
    kallisto_df = pl.read_csv(paths["kallisto_results"], separator="\t", schema_overrides={"sample": pl.String, "target_id": pl.String, "est_counts_mean": pl.Float64})
    kma_df = pl.read_csv(paths["kma_results"], separator="\t", schema_overrides={"sample": pl.String, "#Template": pl.String, "readCount": pl.Float64, "Template_Identity": pl.Float64, "Template_Coverage": pl.Float64, "Query_Identity": pl.Float64}).rename({"#Template": "target_id"})

    kallisto_metric = pl.col("est_counts_mean").cast(pl.Float64, strict=False)
    kma_metric = pl.col("readCount").cast(pl.Float64, strict=False)
    kallisto_df = kallisto_df.with_columns(
        (kallisto_metric.is_not_null() & kallisto_metric.is_finite() & (kallisto_metric > 0)).alias("identified")
    )
    kma_df = kma_df.with_columns(
        (kma_metric.is_not_null() & kma_metric.is_finite() & (kma_metric > 0)).alias("identified")
    )

    kallisto_df = _rename_with_prefix(kallisto_df, "kallisto", keep_cols=["sample", "target_id"])
    kma_df = _rename_with_prefix(kma_df, "kma", keep_cols=["sample", "target_id"])

    joined_df = kallisto_df.join(kma_df, on=["sample", "target_id"], how="full", coalesce=True)
    joined_df = joined_df.with_columns(
        (
            pl.col("kallisto:identified").fill_null(False)
            & pl.col("kma:identified").fill_null(False)
        ).alias("identified_by_both_tools"),
    )
    joined_df = joined_df.select(
        ["sample", "target_id"]
        + [col for col in joined_df.columns if col not in ["sample", "target_id"]]
    )

    kallisto_assigned_df = pl.read_csv(paths["kallisto_assigned"], separator="\t", schema_overrides={"sample": pl.String})
    kma_assigned_df = pl.read_csv(paths["kma_assigned"], separator="\t", schema_overrides={"sample": pl.String})
    kallisto_assigned_df = _rename_mapped_to_assigned(
        _rename_with_prefix(kallisto_assigned_df, "kallisto", keep_cols=["sample"])
    )
    kma_assigned_df = _rename_mapped_to_assigned(
        _rename_with_prefix(kma_assigned_df, "kma", keep_cols=["sample"])
    )
    assigned_df = kallisto_assigned_df.join(kma_assigned_df, on="sample", how="full", coalesce=True)
    assigned_df = assigned_df.select(
        ["sample"]
        + [col for col in assigned_df.columns if col not in ["sample"]]
    )
    return joined_df, assigned_df


def add_kma_identity_across_covered_length(df: pl.DataFrame) -> pl.DataFrame:
    rename_dict = {
        "kma:Template_Identity": "kma:template_identity",
        "kma:Template_Coverage": "kma:template_coverage",
        "kma:Query_Identity": "kma:query_identity",
    }
    rename_dict = {old: new for old, new in rename_dict.items() if old in df.columns}
    df = df.rename(rename_dict)
    return df.with_columns(
        (
            pl.col("kma:template_identity")
            / pl.min_horizontal(pl.lit(100.0), pl.col("kma:template_coverage")).clip(lower_bound=1e-9)
            * 100
        ).alias("kma:identity_across_covered_length")
    )


def make_predictor_table(
    df: pl.DataFrame,
    predictor: str,
    assignment_names: dict[str, str],
    predictor_threshold: float = 0,
    observed_threshold: int = 0,
    filter_threshold_dict: dict | None = None,
) -> pl.DataFrame:
    df = df.with_columns(
        pl.col("target_id").replace_strict(assignment_names).alias("species")
    )
    df = df.with_columns(pl.col(predictor).fill_null(0).fill_nan(0).cast(pl.Float64))

    aggregate_cols = [predictor]
    if "kallisto:est_counts_mean" in df.columns and predictor != "kallisto:est_counts_mean":
        aggregate_cols.append("kallisto:est_counts_mean")
    for col in aggregate_cols:
        df = df.with_columns(pl.col(col).sum().over(["sample", "species"]).alias(col))
    for col in [
        "kma:query_identity",
        "kma:template_identity",
        "kma:template_coverage",
        "kma:identity_across_covered_length",
        "identified_by_both_tools",
    ]:
        if col in df.columns:
            df = df.with_columns(pl.col(col).max().over(["sample", "species"]).alias(col))
    df = df.unique(subset=["sample", "species"])
    df = df.drop("target_id").rename({"species": "target_id"})
    if filter_threshold_dict is not None:
        for col, threshold in filter_threshold_dict.items():
            if isinstance(threshold, bool):
                df = df.filter(pl.col(col) == threshold)
            else:
                df = df.filter(pl.col(col) >= threshold)
    if predictor_threshold > 0:
        df = df.filter(pl.col(predictor) >= predictor_threshold)

    if df.is_empty():
        return pl.DataFrame(schema={"sample": pl.String})
    pa_pivot = df.pivot(values=predictor, index="sample", on="target_id").fill_null(0)
    non_zero_cols = [col for col in pa_pivot.columns if col == "sample" or pa_pivot[col].sum() > 0]
    pa_pivot = pa_pivot.select(non_zero_cols)
    sufficient_obs_cols = [
        col
        for col in pa_pivot.columns
        if col == "sample" or pa_pivot.select((pl.col(col) > 0).sum()).item() >= observed_threshold
    ]
    return pa_pivot.select(["sample"] + sorted(col for col in sufficient_obs_cols if col != "sample"))


def convert_counts_to_abundance(count_df: pl.DataFrame) -> pl.DataFrame:
    count_cols = [col for col in count_df.columns if col != "sample"]
    if not count_cols:
        return count_df
    abundance_df = count_df.with_columns(pl.sum_horizontal([pl.col(col) for col in count_cols]).alias("total_reads"))
    abundance_df = abundance_df.with_columns(
        [
            (pl.col(col) / pl.col("total_reads").cast(pl.Float64).clip(lower_bound=1e-9)).alias(col)
            for col in count_cols
        ]
    )
    return abundance_df.drop("total_reads")


def estimate_abundance(run_dir: str, min_coverage: float, min_identity: float,
                       min_kma_count: float, min_kallisto_count: float,
                       require_both: bool = True, keep_intermediate: bool = False) -> None:
    root = Path(run_dir)
    if not (root / "assignments.tsv").is_file() or not (root / "read-assignment").is_dir():
        raise ValueError("Re-estimation requires a run created with --keep-intermediate")
    names_df = pl.read_csv(root / "assignments.tsv", separator="\t", infer_schema=False)
    assignment_names = dict(zip(names_df["target_id"], names_df["assignment_name"]))
    joined_df, assigned_df = join_kallisto_kma_outputs(root / "read-assignment")
    unknown = set(joined_df["target_id"].to_list()) - set(assignment_names)
    if unknown:
        raise ValueError(f"Assignment results contain unknown reference IDs: {sorted(unknown)}")
    joined_df = add_kma_identity_across_covered_length(joined_df)
    filters = {
        "kma:template_coverage": min_coverage,
        "kma:identity_across_covered_length": min_identity,
        "kma:readCount": min_kma_count,
        "kallisto:est_counts_mean": min_kallisto_count,
    }
    if require_both:
        filters["identified_by_both_tools"] = True
    counts = make_predictor_table(joined_df, "kma:readCount", assignment_names,
                                  filter_threshold_dict=filters)
    samples = pl.read_csv(root / "samples.tsv", separator="\t", infer_schema=False).select("sample")
    counts = samples.join(counts, on="sample", how="left").fill_null(0)
    count_cols = [col for col in counts.columns if col != "sample"]
    # Retain the original paired-end convention: assigned fragments * 2 = reads.
    total = pl.sum_horizontal(count_cols) if count_cols else pl.lit(0.0)
    all_counts = counts.with_columns(total.alias("total_filter_pass_reads"))
    all_counts = all_counts.join(assigned_df.select("sample", "kma:n_assigned"), on="sample", how="left")
    if all_counts["kma:n_assigned"].null_count():
        raise ValueError("Missing KMA assignment totals for one or more samples")
    all_counts = all_counts.with_columns(
        (pl.col("kma:n_assigned") * 2 - pl.col("total_filter_pass_reads")).alias("non_confidently_assigned"))
    if all_counts.filter(pl.col("non_confidently_assigned") < -1e-6).height:
        raise ValueError("Filtered read counts exceed assigned fragments * 2; check KMA count units")
    all_counts = all_counts.with_columns(pl.col("non_confidently_assigned").clip(lower_bound=0))
    all_counts = all_counts.drop("total_filter_pass_reads", "kma:n_assigned")
    summary = samples.join(assigned_df, on="sample", how="left")
    summary = summary.with_columns(
        (pl.col("kma:n_assigned").is_null() | pl.col("kallisto:n_assigned").is_null()).alias("missing_read_assignment"))
    if summary["missing_read_assignment"].any():
        raise ValueError("Missing read-assignment results for one or more samples")
    summary = summary.join(counts.select("sample", total.alias("filter_pass_reads")), on="sample", how="left")
    summary = summary.with_columns((pl.col("filter_pass_reads") == 0).alias("no_filter_pass_reads"))
    out = root / "relative_abundance"
    out.mkdir(exist_ok=True)
    joined_df.with_columns(pl.col("target_id").replace_strict(assignment_names).alias("assignment_name")).write_csv(
        out / "kallisto_kma_outputs_joined.tsv", separator="\t")
    counts.write_csv(out / "filter_pass_read_counts.tsv", separator="\t")
    convert_counts_to_abundance(counts).write_csv(out / "abundance_within_filter_pass_reads.tsv", separator="\t")
    convert_counts_to_abundance(all_counts).write_csv(out / "abundance_within_all_assigned_reads.tsv", separator="\t")
    summary.write_csv(out / "sample_read_assignment_summary.tsv", separator="\t")
    provenance = {}
    for name in ["inputs", "config"]:
        path = root / f"{name}.json"
        if path.exists():
            provenance[name] = json.loads(path.read_text())
    write_json(out / "parameters.json", {"predictor": "kma:readCount", "filters": filters,
                                        "aggregation": "sum counts and maximise quality metrics before filtering",
                                        "assigned_fragment_to_read_multiplier": 2, **provenance})
    if not keep_intermediate:
        clean_run_outputs(root, keep_intermediate=False)
    print(f"Abundance tables written to {out}")
