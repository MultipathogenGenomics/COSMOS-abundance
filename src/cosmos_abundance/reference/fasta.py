# Alignment and shared-region trimming, adapted from COSMOS-analysis.

import subprocess

from Bio import SeqIO
import numpy as np

from ..io import write_fasta


def run_mafft(input_fasta: str, output_fasta: str, threads: int) -> None:
    with open(output_fasta, "w") as out_fasta, open(str(output_fasta) + ".log", "w") as log:
        subprocess.run(["mafft", "--auto", "--thread", str(threads), str(input_fasta)],
                       check=True, stdout=out_fasta, stderr=log)


def find_best_gap_free_region(gap_free_cols: np.ndarray, allow_gap: int) -> tuple[int, int, int]:
    indices = np.flatnonzero(gap_free_cols)
    if len(indices) == 0:
        return 0, -1, 0
    start = end = best_start = best_end = int(indices[0])
    for index in indices[1:]:
        if index - end - 1 <= allow_gap:
            end = int(index)
        else:
            if end - start > best_end - best_start:
                best_start, best_end = start, end
            start = end = int(index)
    if end - start > best_end - best_start:
        best_start, best_end = start, end
    return best_start, best_end, best_end - best_start + 1


def trim_alignment(msa_fasta: str, output_dir, min_length: int, allowed_gap_size: int,
                   min_occupancy: float) -> tuple[list[dict], list[dict], dict]:
    with open(msa_fasta) as handle:
        records = [{"id": record.id, "sequence": str(record.seq).upper()}
                   for record in SeqIO.parse(handle, "fasta")]
    matrix = np.array([list(record["sequence"]) for record in records])
    occupied = np.flatnonzero(np.mean(matrix != "-", axis=0) >= min_occupancy)
    if len(occupied) == 0:
        raise ValueError("No alignment columns meet the minimum occupancy")
    left, right = int(occupied[0]), int(occupied[-1]) + 1
    matrix = matrix[:, left:right]
    kept = list(range(len(records)))
    removed = []
    while kept:
        start, end, length = find_best_gap_free_region(np.all(matrix[kept] != "-", axis=0), allowed_gap_size)
        if length >= min_length:
            break
        worst = max(kept, key=lambda index: (int(np.sum(matrix[index] == "-")), records[index]["id"]))
        removed.append({"id": records[worst]["id"], "reason": "insufficient_shared_region"})
        kept.remove(worst)
    if not kept:
        raise ValueError("No shared region meets --min-reference-length")
    trimmed, ungapped = [], []
    for index in kept:
        sequence = "".join(matrix[index, start:end + 1])
        if len(sequence.replace("-", "")) < min_length:
            removed.append({"id": records[index]["id"], "reason": "short_after_trimming"})
            continue
        trimmed.append({"id": records[index]["id"], "sequence": sequence})
        ungapped.append({"id": records[index]["id"], "sequence": sequence.replace("-", "")})
    if not ungapped:
        raise ValueError("No references meet the minimum length after trimming")
    write_fasta(trimmed, output_dir / "alignment.trimmed.fasta")
    write_fasta(ungapped, output_dir / "references.trimmed.fasta")
    return ungapped, removed, {"alignment_start_1based": left + start + 1,
                              "alignment_end_1based": left + end + 1}
