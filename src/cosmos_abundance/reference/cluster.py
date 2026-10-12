# Cluster references from a distance matrix and select representative sequences.

import numpy as np
import polars as pl
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform


def cluster_sequences_by_distance(ids: list[str], dist_matrix: np.ndarray, distance_threshold: float) -> pl.DataFrame:
    if len(ids) == 1:
        return pl.DataFrame({"seq_id": ids, "cluster_id": [1]})
    condensed_distance = squareform(dist_matrix, checks=False)
    linkage_matrix = sch.linkage(condensed_distance, method="complete")
    cluster_labels = sch.fcluster(linkage_matrix, t=distance_threshold, criterion="distance")
    return pl.DataFrame({"seq_id": ids, "cluster_id": cluster_labels.tolist()})


def identify_best_representatives(cluster_df: pl.DataFrame, dist_matrix: np.ndarray) -> pl.DataFrame:
    ids = cluster_df["seq_id"].to_list()
    best_reps = []
    for cluster_id in sorted(cluster_df["cluster_id"].unique().to_list()):
        members = cluster_df.filter(pl.col("cluster_id") == cluster_id)["seq_id"].to_list()
        if len(members) == 1:
            best_reps.append(members[0])
            continue
        member_indices = [ids.index(member) for member in members]
        sub_dist_matrix = dist_matrix[np.ix_(member_indices, member_indices)]
        avg_distances = sub_dist_matrix.mean(axis=1)
        avg_distance_by_member = dict(zip(members, avg_distances))
        best_reps.append(min(members, key=lambda member: (avg_distance_by_member[member], member)))
    return cluster_df.with_columns(
        pl.col("seq_id").is_in(best_reps).alias("is_best_representative"),
    )


def select_representative_records(records: list[dict], cluster_df: pl.DataFrame) -> list[dict]:
    representative_ids = set(cluster_df.filter(pl.col("is_best_representative"))["seq_id"].to_list())
    return [record for record in records if record["id"] in representative_ids]
