# Pure-python Mash-like distances from k-mer Jaccard similarity.

import math

import numpy as np


VALID_BASES = set("ACGT")


def sequence_to_kmers(sequence: str, kmer_size: int) -> set[str]:
    sequence = sequence.upper().replace("-", "")
    kmers = set()
    for i in range(0, len(sequence) - kmer_size + 1):
        kmer = sequence[i : i + kmer_size]
        if set(kmer).issubset(VALID_BASES):
            kmers.add(kmer)
    return kmers


def jaccard_to_mash_distance(jaccard: float, kmer_size: int) -> float:
    if jaccard <= 0:
        return 1.0
    return max(0.0, -(1 / kmer_size) * math.log((2 * jaccard) / (1 + jaccard)))


def compute_mash_distance_matrix(records: list[dict], kmer_size: int) -> tuple[list[str], np.ndarray]:
    ids = [record["id"] for record in records]
    kmer_sets = [sequence_to_kmers(record["sequence"], kmer_size) for record in records]
    n = len(records)
    dist_matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            union_size = len(kmer_sets[i] | kmer_sets[j])
            if union_size == 0:
                distance = 1.0
            else:
                jaccard = len(kmer_sets[i] & kmer_sets[j]) / union_size
                distance = jaccard_to_mash_distance(jaccard, kmer_size)
            dist_matrix[i, j] = distance
            dist_matrix[j, i] = distance
    return ids, dist_matrix


def write_distance_matrix(ids: list[str], dist_matrix: np.ndarray, output_file: str) -> None:
    with open(output_file, "w") as out_file:
        out_file.write("id\t" + "\t".join(ids) + "\n")
        for seq_id, row in zip(ids, dist_matrix):
            out_file.write(seq_id + "\t" + "\t".join(f"{value:.8g}" for value in row) + "\n")
