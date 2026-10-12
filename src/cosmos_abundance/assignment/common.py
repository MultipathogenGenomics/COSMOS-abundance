# Shared helpers for read-assignment pipeline scripts.

import os


def get_sample_name_from_subdir(file_path: str) -> str:
    return os.path.basename(os.path.dirname(file_path))


def count_reads_in_fastq_gz(fastq_gz: str) -> int:
    import gzip

    opener = gzip.open if str(fastq_gz).endswith(".gz") else open
    with opener(fastq_gz, "rt") as handle:
        return sum(1 for _ in handle) // 4
